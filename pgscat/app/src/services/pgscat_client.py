"""
PGSCatClient service
─────────────────────
Queries the PGS Catalog API for remote score metadata.

Integration modes:
  python (default) → pgscat.pgsclient.PGSCatalogClient (library)
  cli              → subprocess pgscat --json <cmd>

Caching:
  Remote metadata is cached in PostgreSQL (remote_pgs_cache) with a
  configurable TTL (default 24h). This avoids hitting the external API
  repeatedly for the same PGS ID.

Limitations documented:
  - pgscat does NOT support arbitrary free-text search.
  - "text" mode tries PGS ID pattern first, then EFO trait ID.
  - EFO IDs (e.g. EFO_0001360) work; plain English words do not.
  - PMIDs must be numeric strings.
"""
import json
import logging
import subprocess
from typing import Optional, Union

from config import Config

logger = logging.getLogger(__name__)

try:
    from pgscat.pgsclient import PGSCatalogClient as _LibClient  # type: ignore
    _LIB_AVAILABLE = True
    logger.info("pgscat Python library loaded")
except ImportError:
    _LibClient = None  # type: ignore
    _LIB_AVAILABLE = False
    logger.warning("pgscat library not importable; will use CLI fallback")


class PGSCatClient:
    def __init__(self, cfg: Config, db=None) -> None:
        self.cfg = cfg
        self._db = db      # DBPool | None; used for caching + audit
        self._lib: Optional[_LibClient] = None

        if cfg.PGSCAT_MODE == "python" and _LIB_AVAILABLE:
            try:
                self._lib = _LibClient()
                logger.info("pgscat: using Python library")
            except Exception as exc:
                logger.warning("PGSCatalogClient() init failed: %s — using CLI", exc)

    # ── Availability ──────────────────────────────────────────────────────────

    @property
    def available(self) -> bool:
        if not self.cfg.REMOTE_ENABLED:
            return False
        return self._lib is not None or self._cli_probe()

    def _cli_probe(self) -> bool:
        try:
            subprocess.run([self.cfg.PGSCAT_BIN, "--help"],
                           capture_output=True, timeout=5, check=False)
            return True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False

    # ── Public API (with cache) ───────────────────────────────────────────────

    def get_score_with_cache(self, pgs_id: str) -> tuple[dict, bool]:
        """
        Fetch metadata for a single PGS ID.
        Returns (normalised_dict, from_cache).
        """
        if not self.cfg.REMOTE_ENABLED:
            return {"error": "Remote access is disabled", "pgs_id": pgs_id}, False

        # Check DB cache first
        if self._db:
            cached = self._db.get_remote_cache(pgs_id)
            if cached:
                logger.debug("Remote cache hit: %s", pgs_id)
                return self._from_db_cache(cached, pgs_id), True

        # Live query
        result = self.get_score(pgs_id)
        if "error" not in result and self._db:
            self._db.upsert_remote_cache(pgs_id, result, self.cfg.REMOTE_CACHE_TTL)
        return result, False

    def search_with_cache(self, query: str, search_type: str = "text") -> tuple[list[dict], bool]:
        """
        Search with caching (only ID lookups are cacheable per-PGS).
        Returns (results, from_cache).
        """
        # For exact PGS ID searches, use single-item cache
        if search_type == "id":
            result, from_cache = self.get_score_with_cache(query.upper())
            return ([result] if "error" not in result else [result]), from_cache

        results = self.search(query, search_type)
        # Cache individual results that have valid pgs_ids
        if self._db:
            for r in results:
                pid = r.get("pgs_id")
                if pid and "error" not in r:
                    self._db.upsert_remote_cache(pid, r, self.cfg.REMOTE_CACHE_TTL)
        return results, False

    # ── Core query methods ────────────────────────────────────────────────────

    def get_score(self, pgs_id: str) -> dict:
        """Fetch metadata without cache. Returns normalised dict or error dict."""
        if not self.cfg.REMOTE_ENABLED:
            return {"error": "Remote access is disabled", "pgs_id": pgs_id}
        try:
            if self._lib:
                data = self._lib.get_score(pgs_id)
                if not data:
                    return {"error": f"{pgs_id} not found in PGS Catalog", "pgs_id": pgs_id}
                return self._normalize(data)
            return self._cli_get_score(pgs_id)
        except Exception as exc:
            logger.exception("get_score(%s) failed", pgs_id)
            return {"error": str(exc), "pgs_id": pgs_id}

    def search(self, query: str, search_type: str = "text") -> list[dict]:
        """
        Search the PGS Catalog.

        search_type options:
          "id"    → exact PGS ID lookup (reliable)
          "trait" → EFO ontology term, e.g. "EFO_0001360" (reliable)
          "pmid"  → PubMed ID, e.g. "25855707" (reliable)
          "text"  → tries PGS ID pattern first, then EFO; NOT free-text search
        """
        if not self.cfg.REMOTE_ENABLED:
            return [{"error": "Remote access is disabled"}]
        if not query:
            return []
        try:
            if self._lib:
                return self._lib_search(query, search_type)
            return self._cli_search(query, search_type)
        except Exception as exc:
            logger.exception("search(%r, %r) failed", query, search_type)
            return [{"error": str(exc)}]

    # ── Library-mode search ───────────────────────────────────────────────────

    def _lib_search(self, query: str, search_type: str) -> list[dict]:
        client = self._lib
        if search_type == "id":
            r = client.get_score(query.upper())
            return [self._normalize(r)] if r else []
        if search_type == "trait":
            return [self._normalize(r) for r in (client.search_scores(trait_id=query) or [])]
        if search_type == "pmid":
            return [self._normalize(r) for r in (client.search_scores(pmid=query) or [])]

        # "text" mode
        import re
        if re.match(r"^PGS\d{4,6}$", query.upper()):
            r = client.get_score(query.upper())
            if r:
                return [self._normalize(r)]
        results = client.search_scores(trait_id=query) or []
        normalized = [self._normalize(r) for r in results]
        if not normalized:
            normalized = [{
                "error": (
                    f"No results for '{query}'. "
                    "pgscat does not support free-text search. "
                    "Try an EFO ID (e.g. EFO_0001360), a PGS ID (e.g. PGS000004), "
                    "or a PubMed ID."
                )
            }]
        return normalized

    # ── CLI-mode search ───────────────────────────────────────────────────────

    def _cli_get_score(self, pgs_id: str) -> dict:
        raw = self._run_cli(["score", pgs_id])
        if isinstance(raw, dict) and "error" in raw:
            return {**raw, "pgs_id": pgs_id}
        return self._normalize(raw)

    def _cli_search(self, query: str, search_type: str) -> list[dict]:
        if search_type == "id":
            r = self._cli_get_score(query.upper())
            return [r] if "error" not in r else []
        if search_type == "pmid":
            raw = self._run_cli(["search-scores", "--pmid", query])
        else:
            raw = self._run_cli(["search-scores", "--trait", query])
        if isinstance(raw, list):
            return [self._normalize(r) for r in raw]
        if isinstance(raw, dict):
            return [raw] if "error" in raw else [self._normalize(raw)]
        return []

    def _run_cli(self, args: list[str]) -> Union[dict, list]:
        cmd = [self.cfg.PGSCAT_BIN, "--json"] + args
        logger.debug("CLI: %s", " ".join(cmd))
        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True,
                timeout=self.cfg.PGSCAT_TIMEOUT, check=False,
            )
            if proc.returncode != 0:
                err = proc.stderr.strip() or f"pgscat exit {proc.returncode}"
                return {"error": err}
            return json.loads(proc.stdout)
        except subprocess.TimeoutExpired:
            return {"error": f"pgscat CLI timed out after {self.cfg.PGSCAT_TIMEOUT}s"}
        except json.JSONDecodeError as exc:
            return {"error": f"pgscat returned invalid JSON: {exc}"}
        except FileNotFoundError:
            return {"error": f"pgscat binary not found: {self.cfg.PGSCAT_BIN!r}"}

    # ── Normalisation ─────────────────────────────────────────────────────────

    def _normalize(self, raw: dict) -> dict:
        if not raw:
            return {}
        pub = raw.get("publication") or {}
        return {
            "pgs_id":           raw.get("id") or raw.get("pgs_id", ""),
            "name":             raw.get("name", ""),
            "trait_reported":   raw.get("trait_reported", ""),
            "trait_efo":        raw.get("trait_efo", []),
            "variants_number":  raw.get("variants_number"),
            "license":          raw.get("license", ""),
            "ftp_scoring_file": raw.get("ftp_scoring_file", ""),
            "is_harmonized":    (raw.get("harmonization_details") or {}).get("is_harmonized", False),
            "publication": {
                "pmid":    pub.get("PMID") or pub.get("pmid", ""),
                "journal": pub.get("journal", ""),
                "title":   pub.get("title", ""),
                "authors": pub.get("authors", ""),
                "date":    pub.get("date_publication", ""),
            },
            "_raw": raw,
        }

    def _from_db_cache(self, row: dict, pgs_id: str) -> dict:
        """Reconstruct a normalised dict from a remote_pgs_cache row."""
        pub = row.get("publication") or {}
        if isinstance(pub, str):
            import json as _json
            try:
                pub = _json.loads(pub)
            except Exception:
                pub = {}
        return {
            "pgs_id":           pgs_id,
            "name":             row.get("name", ""),
            "trait_reported":   row.get("trait_reported", ""),
            "trait_efo":        row.get("trait_efo") or [],
            "variants_number":  row.get("variants_number"),
            "license":          "",
            "ftp_scoring_file": row.get("ftp_scoring_file", ""),
            "is_harmonized":    row.get("is_harmonized", False),
            "publication":      pub if isinstance(pub, dict) else {},
            "_from_db_cache":   True,
        }
