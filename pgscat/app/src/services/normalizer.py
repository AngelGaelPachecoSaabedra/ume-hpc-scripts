"""
Data normalizer
───────────────
Converts polymorphic pgscat API / DB-cache / local JSON data into
clean, typed structures before template rendering.

Rules:
  - Every field templates touch must be a plain Python str / int / bool / list[str] / dict[str,str].
  - No isinstance() or complex expressions in templates.
  - None-safe everywhere: missing / null / empty → safe default.
"""
from __future__ import annotations

import json
import logging
from typing import Optional

logger = logging.getLogger(__name__)


# ── Scalar helpers ────────────────────────────────────────────────────────────

def normalize_str(val) -> str:
    """Any value → string, strip whitespace, empty string if falsy."""
    if val is None:
        return ""
    return str(val).strip()


def normalize_int(val) -> Optional[int]:
    """Any value → int or None."""
    if val is None:
        return None
    try:
        return int(val)
    except (ValueError, TypeError):
        return None


# ── EFO trait IDs ─────────────────────────────────────────────────────────────

def normalize_efo(raw) -> str:
    """
    Normalise the polymorphic trait_efo field from pgscat / local JSON.

    Accepted input shapes:
      None / ""                                    → ""
      "EFO_0001360"                                → "EFO_0001360"
      ["EFO_0001360", ...]                         → "EFO_0001360, ..."
      [{"id": "EFO_...", "label": "..."}, ...]     → "EFO_... (label), ..."
      [{"id": "EFO_..."}, ...]                     → "EFO_..."

    Returns a plain display string (may be empty).
    """
    if not raw:
        return ""
    if isinstance(raw, str):
        return raw.strip()
    if isinstance(raw, list):
        parts = []
        for item in raw:
            if isinstance(item, dict):
                efo_id = normalize_str(item.get("id") or item.get("efo_id") or "")
                label = normalize_str(item.get("label") or item.get("term") or "")
                if efo_id and label:
                    parts.append(f"{efo_id} ({label})")
                elif efo_id:
                    parts.append(efo_id)
                elif label:
                    parts.append(label)
            elif item:
                parts.append(normalize_str(item))
        return ", ".join(parts)
    return normalize_str(raw)


# ── Publication ───────────────────────────────────────────────────────────────

def normalize_publication(pub) -> dict:
    """
    Normalise publication metadata.

    Accepted input:
      None / {}      → all-empty dict
      dict           → normalised str values
      str (JSON)     → parsed then normalised (DB cache edge case)
    """
    if not pub:
        return {"pmid": "", "journal": "", "title": "", "authors": "", "date": ""}
    if isinstance(pub, str):
        try:
            pub = json.loads(pub)
        except Exception:
            return {"pmid": "", "journal": "", "title": pub, "authors": "", "date": ""}
    if not isinstance(pub, dict):
        return {"pmid": "", "journal": "", "title": "", "authors": "", "date": ""}
    return {
        "pmid":    normalize_str(pub.get("pmid") or pub.get("PMID")),
        "journal": normalize_str(pub.get("journal")),
        "title":   normalize_str(pub.get("title")),
        "authors": normalize_str(pub.get("authors")),
        "date":    normalize_str(pub.get("date") or pub.get("date_publication")),
    }


# ── Chromosomes ───────────────────────────────────────────────────────────────

def normalize_chromosomes(chrs) -> list[str]:
    """Any input → list of strings (never list-of-non-string)."""
    if not chrs:
        return []
    if isinstance(chrs, str):
        # JSON-encoded list from PostgreSQL JSONB column
        try:
            chrs = json.loads(chrs)
        except Exception:
            return [chrs] if chrs else []
    if isinstance(chrs, list):
        return [normalize_str(c) for c in chrs if c is not None]
    return []


# ── Composite normalizers ─────────────────────────────────────────────────────

def normalize_local_info(info: dict) -> dict:
    """
    Normalise the dict returned by LocalCatalog.get_pgs_info() before
    passing to dashboard.tpl / pgs_remote.tpl.

    Adds convenience fields:
      efo_display   – ready-to-render EFO string
      chrom_display – comma-joined chromosome list or '—'
    """
    if not info:
        return info

    chrs = normalize_chromosomes(info.get("chromosomes", []))
    info["chromosomes"] = chrs
    info["chrom_display"] = ", ".join(chrs) if chrs else "—"
    info["efo_display"] = normalize_efo(info.get("trait_efo"))
    info["n_variants"] = normalize_int(info.get("n_variants"))

    # total_metadata sub-dict
    tmeta = info.get("total_metadata")
    if isinstance(tmeta, str):
        try:
            tmeta = json.loads(tmeta)
        except Exception:
            tmeta = {}
    info["total_metadata"] = tmeta or {}

    return info


def normalize_remote_info(info: dict) -> dict:
    """
    Normalise the dict returned by PGSCatClient.get_score_with_cache()
    before passing to pgs_remote.tpl / search.tpl.
    """
    if not info or "error" in info:
        return info

    info["efo_display"] = normalize_efo(info.get("trait_efo"))
    info["variants_number"] = normalize_int(info.get("variants_number"))
    info["publication"] = normalize_publication(info.get("publication"))
    return info


def normalize_search_result(result: dict) -> dict:
    """
    Normalise a single search result dict for search.tpl.
    Safe to call even if 'error' key is present.
    """
    if not result:
        return result
    result["variants_number"] = normalize_int(result.get("variants_number"))
    if "publication" in result:
        result["publication"] = normalize_publication(result.get("publication"))
    return result
