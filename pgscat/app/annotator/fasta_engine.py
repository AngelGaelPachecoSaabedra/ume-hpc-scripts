"""
FASTA-Aware Variant Consequence Engine
=======================================
Computes biologically correct coding consequences and splice-site classification
using pyfaidx for random-access FASTA queries.

Capabilities:
  - Determines REF vs ALT by matching alleles against the reference FASTA
  - Extracts the reference codon containing a CDS variant
  - Applies the ALT substitution in coding-strand orientation
  - Translates both codons and classifies: synonymous / missense / stop_gained /
    stop_lost / start_lost (SNPs) or frameshift / inframe_insertion/deletion (indels)
  - Classifies intronic variants near splice sites as
    splice_donor_variant or splice_acceptor_variant by reading the canonical
    GT/AG dinucleotides directly from the FASTA

Requires:
  - pyfaidx >= 0.5
  - A GFF3Index with .transcript_cds_map populated (gff3_parser v1.2+)

FASTA chromosome convention:
  - hg38.fa uses 'chr' prefix:  chr1, chr2, ..., chrX, chrY, chrM
  - Our data uses bare names:   1, 2, ..., X, Y, MT
  - This module handles the conversion transparently.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from gff3_parser import GFF3Index

from sequence_utils import (
    reverse_complement,
    complement_nt,
    translate_codon,
    aa_3letter,
    snp_codon_consequence,
    indel_consequence,
)

logger = logging.getLogger(__name__)

# Chromosome name remapping for FASTA lookup
_BARE_TO_FASTA: dict[str, str] = {"MT": "chrM"}


def _to_fasta_chrom(bare_chrom: str) -> str:
    """Convert bare chromosome name to FASTA contig name (adds 'chr' prefix)."""
    if bare_chrom.startswith("chr"):
        return bare_chrom
    return _BARE_TO_FASTA.get(bare_chrom, f"chr{bare_chrom}")


def _strip_chr(name: str) -> str:
    """Remove 'chr' prefix; remap M→MT."""
    bare = name[3:] if name.startswith("chr") else name
    return {"M": "MT"}.get(bare, bare)


class FASTAEngine:
    """
    FASTA-aware variant consequence engine.

    Usage:
        engine = FASTAEngine("/ref/GRCh38.fa")
        result = engine.get_codon_consequence(
            chrom="1", pos=149934520,
            effect_allele="C", other_allele="T",
            transcript_id="ENST00000853977",
            gff3_index=idx,
        )
        splice_type = engine.get_splice_type("1", 12345678, gff3_index)
    """

    def __init__(self, fasta_path: Path) -> None:
        from pyfaidx import Fasta
        self.fasta_path = Path(fasta_path)
        logger.info("Opening FASTA: %s", self.fasta_path)
        self.fa = Fasta(str(self.fasta_path), as_raw=False, sequence_always_upper=True)
        logger.info("FASTA loaded.  Chromosomes available: %d", len(self.fa.keys()))

    # ── Public API ────────────────────────────────────────────────────────────

    def get_ref_allele(self, chrom: str, pos: int, length: int = 1) -> Optional[str]:
        """
        Fetch the reference sequence at 1-based position pos (length bp).
        Returns uppercase sequence string, or None on error.
        """
        fchrom = _to_fasta_chrom(_strip_chr(str(chrom)))
        try:
            seq = str(self.fa[fchrom][pos - 1: pos - 1 + length]).upper()
            return seq if seq else None
        except (KeyError, ValueError, Exception):
            return None

    def determine_ref_alt(
        self, chrom: str, pos: int, allele1: str, allele2: str
    ) -> tuple[Optional[str], Optional[str]]:
        """
        Determine which of allele1/allele2 is the REF allele by checking the FASTA.

        Returns:
            (ref_allele, alt_allele) where ref_allele matches the FASTA,
            or (None, None) if neither matches.
        """
        fasta_ref = self.get_ref_allele(chrom, pos, len(allele1))
        if fasta_ref and fasta_ref.upper() == allele1.upper():
            return allele1.upper(), allele2.upper()

        fasta_ref2 = self.get_ref_allele(chrom, pos, len(allele2))
        if fasta_ref2 and fasta_ref2.upper() == allele2.upper():
            return allele2.upper(), allele1.upper()

        return None, None

    def get_codon_consequence(
        self,
        chrom: str,
        pos: int,
        effect_allele: str,
        other_allele: str,
        transcript_id: str,
        gff3_index: "GFF3Index",
    ) -> dict:
        """
        Compute the coding consequence of a CDS variant.

        Args:
            chrom:          Chromosome (bare, no 'chr')
            pos:            1-based genomic position
            effect_allele:  EFFECT_ALLELE from betamap
            other_allele:   OTHER_ALLELE from betamap
            transcript_id:  ENST ID of the annotating transcript (no version suffix)
            gff3_index:     GFF3Index with transcript_cds_map populated

        Returns:
            dict with keys:
                ref_allele, alt_allele,
                codon_ref, codon_alt,
                aa_ref (1-letter), aa_alt (1-letter),
                aa_ref_3, aa_alt_3 (3-letter),
                consequence,
                codon_index, pos_in_codon,
                error (only present on failure)
        """
        # 1. Determine REF / ALT from FASTA
        ref_allele, alt_allele = self.determine_ref_alt(
            chrom, pos, effect_allele, other_allele
        )
        base = {
            "ref_allele":  ref_allele or effect_allele,
            "alt_allele":  alt_allele or other_allele,
            "codon_ref":   "",
            "codon_alt":   "",
            "aa_ref":      "",
            "aa_alt":      "",
            "aa_ref_3":    "",
            "aa_alt_3":    "",
            "consequence": "coding_sequence_variant",
            "codon_index":  None,
            "pos_in_codon": None,
        }

        if ref_allele is None:
            base["error"] = "cannot_determine_ref_alt"
            return base

        # 2. Indel — frameshift or inframe
        if len(ref_allele) != 1 or len(alt_allele) != 1:
            base["consequence"] = indel_consequence(ref_allele, alt_allele)
            return base

        # 3. SNP — look up CDS offset
        cds_info = gff3_index.get_cds_offset(transcript_id, chrom, pos)
        if cds_info is None:
            base["error"] = "cds_offset_not_found"
            return base

        codon_index  = cds_info["codon_index"]
        pos_in_codon = cds_info["pos_in_codon"]
        strand       = cds_info["strand"]

        base["codon_index"]  = codon_index
        base["pos_in_codon"] = pos_in_codon

        # 4. Extract reference codon from FASTA
        ref_codon = self._get_codon_seq(transcript_id, codon_index, gff3_index)
        if ref_codon is None:
            base["error"] = "codon_extraction_failed"
            return base

        # 5. Build ALT codon
        # ALT allele is always on the + (forward) strand in betamap.
        # The coding codon is on the coding strand:
        #   + strand: same as forward → use ALT directly
        #   - strand: coding strand is RC of forward → complement ALT
        coding_alt = complement_nt(alt_allele) if strand == "-" else alt_allele
        codon_list = list(ref_codon.upper())
        codon_list[pos_in_codon] = coding_alt.upper()
        alt_codon = "".join(codon_list)

        # 6. Translate and classify
        consequence = snp_codon_consequence(ref_codon, alt_codon, codon_index)
        ref_aa = translate_codon(ref_codon)
        alt_aa = translate_codon(alt_codon)

        base.update({
            "codon_ref":  ref_codon,
            "codon_alt":  alt_codon,
            "aa_ref":     ref_aa,
            "aa_alt":     alt_aa,
            "aa_ref_3":   aa_3letter(ref_aa),
            "aa_alt_3":   aa_3letter(alt_aa),
            "consequence": consequence,
        })
        return base

    def get_splice_type(
        self,
        chrom: str,
        pos: int,
        gff3_index: "GFF3Index",
        max_dist: int = 10,
    ) -> Optional[str]:
        """
        For a variant near a splice site, classify as 'donor' or 'acceptor'.

        Uses the canonical GT (donor) and AG (acceptor) dinucleotides at
        exon-intron boundaries, fetched from the FASTA.

        Args:
            chrom:    Bare chromosome name
            pos:      1-based genomic position
            gff3_index: GFF3Index with exon interval tree
            max_dist: Maximum distance from exon boundary to consider (default 10)

        Returns:
            'donor' | 'acceptor' | None
        """
        chrom_bare = _strip_chr(str(chrom))
        if chrom_bare not in gff3_index.trees:
            return None
        if "exon" not in gff3_index.trees[chrom_bare]:
            return None

        exon_tree = gff3_index.trees[chrom_bare]["exon"]
        qp = pos - 1  # 0-based

        # Find nearest exon boundary within max_dist
        best_dist     = None
        best_boundary = None   # "start" or "end"
        best_iv       = None

        for iv in exon_tree[max(0, qp - max_dist): qp + max_dist + 1]:
            d_begin = abs(qp - iv.begin)
            d_end   = abs(qp - (iv.end - 1))
            if d_begin <= d_end:
                d, btype = d_begin, "start"
            else:
                d, btype = d_end, "end"
            if best_dist is None or d < best_dist:
                best_dist, best_boundary, best_iv = d, btype, iv

        if best_iv is None:
            return None

        strand = best_iv.data.strand
        fchrom = _to_fasta_chrom(chrom_bare)

        try:
            return self._classify_splice_boundary(
                fchrom, strand, best_boundary, best_iv
            )
        except Exception as exc:
            logger.debug("splice_type lookup failed at %s:%d: %s", chrom, pos, exc)
            return None

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _get_codon_seq(
        self,
        transcript_id: str,
        codon_index: int,
        gff3_index: "GFF3Index",
    ) -> Optional[str]:
        """
        Extract the 3-nt codon at codon_index from the FASTA.
        Handles multi-exon codons and strand orientation.
        Returns uppercase codon string or None on failure.
        """
        intervals = gff3_index.transcript_cds_map.get(transcript_id, [])
        if not intervals:
            return None

        strand = intervals[0][3]
        if strand == "+":
            sorted_ivs = sorted(intervals, key=lambda x: x[1])   # asc by begin
        else:
            sorted_ivs = sorted(intervals, key=lambda x: -x[2])  # desc by end

        target_s = codon_index * 3
        target_e = target_s + 3

        cumulative = 0
        nts: list[str] = []

        for (iv_chrom, begin, end, _) in sorted_ivs:
            exon_len   = end - begin
            exon_cds_s = cumulative
            exon_cds_e = cumulative + exon_len

            ov_s = max(target_s, exon_cds_s)
            ov_e = min(target_e, exon_cds_e)

            if ov_s < ov_e:
                rel_s = ov_s - exon_cds_s
                rel_e = ov_e - exon_cds_s
                fchrom = _to_fasta_chrom(iv_chrom)

                if strand == "+":
                    g_s = begin + rel_s
                    g_e = begin + rel_e
                    seq = str(self.fa[fchrom][g_s:g_e]).upper()
                else:
                    # - strand: rel_s=0 maps to the rightmost nt of this exon
                    # Genomic range: [end - rel_e : end - rel_s]
                    g_s = end - rel_e
                    g_e = end - rel_s
                    seq = str(self.fa[fchrom][g_s:g_e]).upper()
                    seq = reverse_complement(seq)

                nts.extend(seq)

            cumulative += exon_len
            if cumulative >= target_e:
                break

        return "".join(nts[:3]) if len(nts) >= 3 else None

    def _classify_splice_boundary(
        self,
        fchrom: str,
        strand: str,
        boundary: str,
        iv,
    ) -> Optional[str]:
        """
        Read the canonical dinucleotide at a splice boundary and classify.

        boundary='start' → near iv.begin (genomic start of exon)
        boundary='end'   → near iv.end-1 (genomic end of exon)

        + strand:
          end   boundary → downstream intron → GT = donor
          start boundary → upstream   intron → AG = acceptor

        - strand (transcriptomic direction reversed):
          start boundary → downstream intron (lower coords) → GT = donor
          end   boundary → upstream   intron (higher coords) → AG = acceptor
        """
        if strand == "+":
            if boundary == "end":
                # Intron starts at iv.end (0-based)
                di = str(self.fa[fchrom][iv.end: iv.end + 2]).upper()
                return "donor" if di == "GT" else None
            else:
                # Intron ends just before iv.begin
                if iv.begin < 2:
                    return None
                di = str(self.fa[fchrom][iv.begin - 2: iv.begin]).upper()
                return "acceptor" if di == "AG" else None
        else:  # strand == "-"
            if boundary == "start":
                # Intron at lower coords than iv.begin
                if iv.begin < 2:
                    return None
                di = str(self.fa[fchrom][iv.begin - 2: iv.begin]).upper()
                return "donor" if reverse_complement(di) == "GT" else None
            else:
                # Intron at higher coords than iv.end
                di = str(self.fa[fchrom][iv.end: iv.end + 2]).upper()
                return "acceptor" if reverse_complement(di) == "AG" else None
