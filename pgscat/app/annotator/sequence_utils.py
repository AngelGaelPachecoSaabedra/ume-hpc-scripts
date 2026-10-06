"""
Sequence Utilities
==================
Pure functions for nucleotide manipulation, codon translation, and
coding-consequence prediction.  No external dependencies.

Used by fasta_engine.py for FASTA-aware variant annotation.
"""
from __future__ import annotations

from typing import Optional

# ── Codon table (NCBI Standard genetic code, 1-letter AA) ────────────────────
CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# 1-letter → 3-letter amino acid names
_AA_3LETTER: dict[str, str] = {
    "A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe",
    "G": "Gly", "H": "His", "I": "Ile", "K": "Lys", "L": "Leu",
    "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln", "R": "Arg",
    "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr",
    "*": "Ter", "?": "???",
}

# Complement table (handles N and lowercase)
_COMPLEMENT = str.maketrans("ACGTacgtNnRrYySsWwKkMmBbDdHhVv",
                             "TGCAtgcaNnYyRrSsWwMmKkVvHhDdBb")


# ── Nucleotide helpers ────────────────────────────────────────────────────────

def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


def complement_nt(nt: str) -> str:
    """Return the complement of a single nucleotide."""
    return nt.translate(_COMPLEMENT)


# ── Codon helpers ─────────────────────────────────────────────────────────────

def translate_codon(codon: str) -> str:
    """
    Translate a 3-nt codon to 1-letter amino acid code.
    Returns '*' for stop codons, '?' for unknown/degenerate codons.
    """
    return CODON_TABLE.get(codon.upper(), "?")


def aa_3letter(aa_1: str) -> str:
    """Convert 1-letter AA code to 3-letter name (e.g. 'M' → 'Met')."""
    return _AA_3LETTER.get(aa_1, aa_1)


# ── Consequence logic ─────────────────────────────────────────────────────────

def snp_codon_consequence(
    ref_codon: str,
    alt_codon: str,
    codon_index: int,
) -> str:
    """
    Determine SO coding consequence from REF and ALT codons.

    Args:
        ref_codon:    3-nt reference codon in coding-strand orientation.
        alt_codon:    3-nt alternative codon in coding-strand orientation.
        codon_index:  0-based codon index in CDS (0 = first Met codon).

    Returns:
        SO consequence term:
          synonymous_variant | missense_variant | stop_gained |
          stop_lost | start_lost
    """
    ref_codon = ref_codon.upper()
    alt_codon = alt_codon.upper()

    ref_aa = translate_codon(ref_codon)
    alt_aa = translate_codon(alt_codon)

    # Start codon loss: codon 0 was Met, now is not
    if codon_index == 0 and ref_aa == "M" and alt_aa != "M":
        return "start_lost"

    # Stop → non-stop
    if ref_aa == "*" and alt_aa != "*":
        return "stop_lost"

    # Non-stop → stop
    if ref_aa != "*" and alt_aa == "*":
        return "stop_gained"

    # Synonymous
    if ref_aa == alt_aa:
        return "synonymous_variant"

    # Missense
    return "missense_variant"


def indel_consequence(ref_allele: str, alt_allele: str) -> str:
    """
    Determine consequence for an insertion or deletion in CDS.

    Args:
        ref_allele: REF allele string
        alt_allele: ALT allele string

    Returns:
        frameshift_variant | inframe_insertion | inframe_deletion |
        coding_sequence_variant (for complex MNPs)
    """
    len_diff = len(alt_allele) - len(ref_allele)
    if len_diff == 0:
        return "coding_sequence_variant"   # complex substitution / MNP
    if len_diff % 3 != 0:
        return "frameshift_variant"
    if len_diff > 0:
        return "inframe_insertion"
    return "inframe_deletion"


# ── Loss-of-function helper ───────────────────────────────────────────────────

LOF_CONSEQUENCES = frozenset({
    "stop_gained",
    "stop_lost",
    "start_lost",
    "frameshift_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
})


def is_lof_consequence(consequence: str) -> bool:
    """Return True if the consequence is considered loss-of-function."""
    return consequence in LOF_CONSEQUENCES
