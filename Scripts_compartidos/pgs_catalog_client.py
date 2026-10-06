#!/usr/bin/env python3
"""
pgs_catalog_client.py
=====================
Cliente Python para la API REST del PGS Catalog (https://www.pgscatalog.org/rest/).

Permite consultar los cinco tipos de entidades del catálogo:
  - Polygenic Scores  (PGS)
  - Publications       (PGP)
  - Traits / EFO       (EFO)
  - Performance Metrics (PPM)
  - Sample Sets         (PSS)

Uso interactivo (menú) o como módulo importable.

Requisitos: requests  (pip install requests)
Autor:      Generado para el proyecto MCPS-HPC
Fecha:      2026-03-24
"""

import json
import csv
import sys
import time
import argparse
from pathlib import Path
from typing import Any

try:
    import requests
except ImportError:
    sys.exit("Error: instala 'requests' con:  pip install requests")


# ── Configuración ────────────────────────────────────────────────────────────
BASE_URL = "https://www.pgscatalog.org/rest"
HEADERS  = {"Accept": "application/json"}
RATE_LIMIT_PAUSE = 0.25          # segundos entre peticiones paginadas
DEFAULT_TIMEOUT  = 30            # segundos


# ── Clase principal ──────────────────────────────────────────────────────────

class PGSCatalogClient:
    """Cliente para la REST API del PGS Catalog."""

    def __init__(self, base_url: str = BASE_URL, timeout: int = DEFAULT_TIMEOUT):
        self.base_url = base_url.rstrip("/")
        self.timeout  = timeout
        self.session  = requests.Session()
        self.session.headers.update(HEADERS)

    # ── Métodos internos ─────────────────────────────────────────────────────

    def _get(self, endpoint: str, params: dict | None = None) -> dict:
        """GET con manejo de errores y rate-limit."""
        url = f"{self.base_url}/{endpoint.lstrip('/')}"
        try:
            r = self.session.get(url, params=params, timeout=self.timeout)
            r.raise_for_status()
            return r.json()
        except requests.exceptions.HTTPError as e:
            print(f"  ⚠ HTTP {r.status_code}: {r.text[:200]}")
            raise
        except requests.exceptions.RequestException as e:
            print(f"  ⚠ Error de conexión: {e}")
            raise

    def _get_all_pages(self, endpoint: str, params: dict | None = None,
                       limit: int | None = None) -> list[dict]:
        """Recorre todas las páginas de un endpoint paginado.

        El PGS Catalog usa el esquema:
          { "count": N, "next": url|null, "previous": url|null, "results": [...] }

        Parameters
        ----------
        limit : int, optional
            Número máximo de registros a devolver (None = todos).
        """
        all_results: list[dict] = []
        data = self._get(endpoint, params)

        # Si la respuesta no es paginada (es un solo objeto), devolver tal cual
        if "results" not in data:
            return [data]

        all_results.extend(data["results"])
        total = data.get("count", "?")
        print(f"  → Página 1  ({len(all_results)}/{total} registros)")

        while data.get("next") and (limit is None or len(all_results) < limit):
            time.sleep(RATE_LIMIT_PAUSE)
            try:
                r = self.session.get(data["next"], timeout=self.timeout)
                r.raise_for_status()
                data = r.json()
            except requests.exceptions.RequestException as e:
                print(f"  ⚠ Error al paginar: {e}")
                break
            all_results.extend(data["results"])
            page = len(all_results) // max(len(data["results"]), 1) + 1
            print(f"  → Página {page}  ({len(all_results)}/{total} registros)")

        if limit:
            all_results = all_results[:limit]
        return all_results

    # ── Scores ───────────────────────────────────────────────────────────────

    def get_score(self, pgs_id: str) -> dict:
        """Obtener metadatos de un PGS específico (e.g. PGS000001)."""
        return self._get(f"score/{pgs_id}")

    def get_all_scores(self, limit: int | None = None) -> list[dict]:
        """Listar todos los scores del catálogo (paginado)."""
        return self._get_all_pages("score/all", limit=limit)

    def search_scores(self, trait_id: str | None = None,
                      pmid: str | None = None,
                      limit: int | None = None) -> list[dict]:
        """Buscar scores por trait EFO ID y/o PubMed ID."""
        params = {}
        if trait_id:
            params["trait_id"] = trait_id
        if pmid:
            params["pmid"] = pmid
        if not params:
            raise ValueError("Debes indicar al menos trait_id o pmid")
        return self._get_all_pages("score/search", params=params, limit=limit)

    # ── Publicaciones ────────────────────────────────────────────────────────

    def get_publication(self, pgp_id: str) -> dict:
        """Obtener metadatos de una publicación (e.g. PGP000001)."""
        return self._get(f"publication/{pgp_id}")

    def get_all_publications(self, limit: int | None = None) -> list[dict]:
        """Listar todas las publicaciones."""
        return self._get_all_pages("publication/all", limit=limit)

    def search_publications(self, pmid: str | None = None,
                            limit: int | None = None) -> list[dict]:
        """Buscar publicaciones por PubMed ID."""
        params = {}
        if pmid:
            params["pmid"] = pmid
        return self._get_all_pages("publication/search", params=params, limit=limit)

    # ── Traits ───────────────────────────────────────────────────────────────

    def get_trait(self, efo_id: str) -> dict:
        """Obtener información de un trait EFO (e.g. EFO_0000305)."""
        return self._get(f"trait/{efo_id}")

    def get_all_traits(self, limit: int | None = None) -> list[dict]:
        """Listar todos los traits del catálogo."""
        return self._get_all_pages("trait/all", limit=limit)

    # ── Performance Metrics ──────────────────────────────────────────────────

    def search_performance(self, pgs_id: str | None = None,
                           pmid: str | None = None,
                           limit: int | None = None) -> list[dict]:
        """Buscar métricas de rendimiento de un PGS o publicación."""
        params = {}
        if pgs_id:
            params["pgs_id"] = pgs_id
        if pmid:
            params["pmid"] = pmid
        return self._get_all_pages("performance/search", params=params, limit=limit)

    # ── Sample Sets ──────────────────────────────────────────────────────────

    def get_sample_set(self, pss_id: str) -> dict:
        """Obtener un sample set específico (e.g. PSS000001)."""
        return self._get(f"sample_set/{pss_id}")

    def get_all_sample_sets(self, limit: int | None = None) -> list[dict]:
        """Listar todos los sample sets."""
        return self._get_all_pages("sample_set/all", limit=limit)

    # ── Información general ──────────────────────────────────────────────────

    def get_ancestry_categories(self) -> dict:
        """Obtener las categorías de ancestría definidas por el catálogo."""
        return self._get("ancestry_categories")

    def get_release_info(self) -> dict:
        """Obtener información de la versión/release actual del catálogo."""
        return self._get("info")


# ── Utilidades de exportación ────────────────────────────────────────────────

def flatten_record(record: dict, parent_key: str = "", sep: str = ".") -> dict:
    """Aplana un dict anidado para poder escribirlo como fila CSV."""
    items: list[tuple[str, Any]] = []
    for k, v in record.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_record(v, new_key, sep).items())
        elif isinstance(v, list):
            # Listas de dicts → concatenar IDs o representar como JSON compacto
            if v and isinstance(v[0], dict):
                items.append((new_key, json.dumps(v, ensure_ascii=False)))
            else:
                items.append((new_key, "; ".join(str(x) for x in v)))
        else:
            items.append((new_key, v))
    return dict(items)


def export_json(data: list[dict] | dict, filepath: str) -> None:
    """Exportar resultados a JSON."""
    path = Path(filepath)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"  ✓ Exportado a {path}  ({path.stat().st_size / 1024:.1f} KB)")


def export_csv(data: list[dict], filepath: str) -> None:
    """Exportar resultados (lista de dicts) a CSV."""
    if not data:
        print("  ⚠ Sin datos para exportar.")
        return
    flat = [flatten_record(r) for r in data]
    # Obtener todas las columnas posibles
    all_keys: list[str] = []
    seen: set[str] = set()
    for row in flat:
        for k in row:
            if k not in seen:
                all_keys.append(k)
                seen.add(k)

    path = Path(filepath)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=all_keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(flat)
    print(f"  ✓ Exportado a {path}  ({path.stat().st_size / 1024:.1f} KB)")


# ── Funciones de presentación ────────────────────────────────────────────────

def print_score_summary(score: dict) -> None:
    """Imprime un resumen legible de un PGS."""
    traits = score.get("trait_efo", [])
    trait_names = ", ".join(t.get("label", "?") for t in traits) if traits else "N/A"
    pub = score.get("publication", {}) or {}
    print(f"  ┌─ {score.get('id', '?')}  {score.get('name', '')}")
    print(f"  │  Trait reportado : {score.get('trait_reported', 'N/A')}")
    print(f"  │  Trait(s) EFO    : {trait_names}")
    print(f"  │  Variantes       : {score.get('variants_number', '?')}")
    print(f"  │  Genome build    : {score.get('variants_genomebuild', 'N/A')}")
    print(f"  │  Weight type     : {score.get('weight_type', 'NR')}")
    print(f"  │  Publicación     : {pub.get('id', '?')} — {pub.get('title', '')[:70]}")
    print(f"  │  PMID            : {pub.get('PMID', 'N/A')}")
    ftp = score.get("ftp_scoring_file", "")
    print(f"  └─ Scoring file   : {ftp}")


def print_publication_summary(pub: dict) -> None:
    """Imprime un resumen de una publicación PGP."""
    print(f"  ┌─ {pub.get('id', '?')}")
    print(f"  │  Título    : {pub.get('title', 'N/A')[:80]}")
    print(f"  │  Autores   : {pub.get('firstauthor', '?')} et al.")
    print(f"  │  Journal   : {pub.get('journal', '?')} ({pub.get('date_publication', '?')})")
    print(f"  │  PMID      : {pub.get('PMID', 'N/A')}")
    print(f"  │  DOI       : {pub.get('doi', 'N/A')}")
    n_scores = len(pub.get("associated_pgs_ids", {}).get("development", []))
    n_eval   = len(pub.get("associated_pgs_ids", {}).get("evaluation", []))
    print(f"  └─ PGS dev/eval: {n_scores} / {n_eval}")


def print_trait_summary(trait: dict) -> None:
    """Imprime un resumen de un trait."""
    print(f"  ┌─ {trait.get('id', '?')}")
    print(f"  │  Label       : {trait.get('label', 'N/A')}")
    print(f"  │  Descripción : {(trait.get('description', '') or '')[:90]}")
    print(f"  │  Categoría   : {trait.get('category', 'N/A')}")
    n_scores = len(trait.get("associated_pgs_ids", []))
    print(f"  └─ # PGS asociados: {n_scores}")


# ── Menú interactivo ─────────────────────────────────────────────────────────

def interactive_menu():
    """Menú interactivo por consola."""
    client = PGSCatalogClient()

    MENU = """
╔══════════════════════════════════════════════════════════════╗
║              PGS CATALOG — Cliente REST API                 ║
╠══════════════════════════════════════════════════════════════╣
║  SCORES                                                     ║
║   1. Consultar score por PGS ID  (e.g. PGS000001)          ║
║   2. Buscar scores por trait     (e.g. EFO_0001360)         ║
║   3. Buscar scores por PMID      (e.g. 25855707)           ║
║   4. Listar todos los scores     (paginado)                 ║
║                                                              ║
║  PUBLICACIONES                                               ║
║   5. Consultar publicación por PGP ID                       ║
║   6. Buscar publicación por PMID                            ║
║   7. Listar todas las publicaciones                         ║
║                                                              ║
║  TRAITS                                                      ║
║   8. Consultar trait por EFO ID                             ║
║   9. Listar todos los traits                                ║
║                                                              ║
║  PERFORMANCE & SAMPLES                                       ║
║  10. Métricas de rendimiento por PGS ID                     ║
║  11. Consultar sample set por PSS ID                        ║
║                                                              ║
║  GENERAL                                                     ║
║  12. Info del catálogo (release actual)                      ║
║  13. Categorías de ancestría                                ║
║                                                              ║
║  EXPORTAR                                                    ║
║  14. Exportar última consulta a JSON                        ║
║  15. Exportar última consulta a CSV                         ║
║                                                              ║
║   0. Salir                                                   ║
╚══════════════════════════════════════════════════════════════╝
"""
    last_results: list[dict] | dict | None = None

    while True:
        print(MENU)
        choice = input("Opción ▸ ").strip()

        try:
            if choice == "0":
                print("¡Hasta luego!")
                break

            elif choice == "1":
                pgs_id = input("  PGS ID (e.g. PGS000001): ").strip().upper()
                result = client.get_score(pgs_id)
                print_score_summary(result)
                last_results = result

            elif choice == "2":
                efo_id = input("  Trait EFO ID (e.g. EFO_0001360): ").strip()
                lim = input("  Límite de resultados (Enter=todos): ").strip()
                limit = int(lim) if lim else None
                results = client.search_scores(trait_id=efo_id, limit=limit)
                print(f"\n  Encontrados: {len(results)} scores")
                for s in results[:20]:
                    print_score_summary(s)
                if len(results) > 20:
                    print(f"  ... y {len(results)-20} más.")
                last_results = results

            elif choice == "3":
                pmid = input("  PubMed ID: ").strip()
                results = client.search_scores(pmid=pmid)
                print(f"\n  Encontrados: {len(results)} scores")
                for s in results:
                    print_score_summary(s)
                last_results = results

            elif choice == "4":
                lim = input("  Límite de resultados (Enter=50): ").strip()
                limit = int(lim) if lim else 50
                results = client.get_all_scores(limit=limit)
                print(f"\n  Descargados: {len(results)} scores")
                for s in results[:10]:
                    print_score_summary(s)
                if len(results) > 10:
                    print(f"  ... mostrando 10 de {len(results)}.")
                last_results = results

            elif choice == "5":
                pgp_id = input("  PGP ID (e.g. PGP000001): ").strip().upper()
                result = client.get_publication(pgp_id)
                print_publication_summary(result)
                last_results = result

            elif choice == "6":
                pmid = input("  PubMed ID: ").strip()
                results = client.search_publications(pmid=pmid)
                print(f"\n  Encontradas: {len(results)} publicaciones")
                for p in results:
                    print_publication_summary(p)
                last_results = results

            elif choice == "7":
                lim = input("  Límite de resultados (Enter=50): ").strip()
                limit = int(lim) if lim else 50
                results = client.get_all_publications(limit=limit)
                print(f"\n  Descargadas: {len(results)} publicaciones")
                for p in results[:10]:
                    print_publication_summary(p)
                last_results = results

            elif choice == "8":
                efo_id = input("  EFO ID (e.g. EFO_0000305): ").strip()
                result = client.get_trait(efo_id)
                print_trait_summary(result)
                last_results = result

            elif choice == "9":
                lim = input("  Límite de resultados (Enter=50): ").strip()
                limit = int(lim) if lim else 50
                results = client.get_all_traits(limit=limit)
                print(f"\n  Descargados: {len(results)} traits")
                for t in results[:15]:
                    print_trait_summary(t)
                last_results = results

            elif choice == "10":
                pgs_id = input("  PGS ID (e.g. PGS000001): ").strip().upper()
                results = client.search_performance(pgs_id=pgs_id)
                print(f"\n  Encontradas: {len(results)} métricas")
                for pm in results[:10]:
                    ppm_id = pm.get("id", "?")
                    effect = pm.get("effect_sizes", [])
                    class_acc = pm.get("class_acc", [])
                    other = pm.get("othermetrics", [])
                    print(f"  ─ {ppm_id}  effect_sizes={effect}  "
                          f"class_acc={class_acc}  other={other}")
                last_results = results

            elif choice == "11":
                pss_id = input("  PSS ID (e.g. PSS000001): ").strip().upper()
                result = client.get_sample_set(pss_id)
                print(f"  Sample set: {result.get('id', '?')}")
                samples = result.get("samples", [])
                for s in samples:
                    anc = s.get("ancestry_broad", "?")
                    n   = s.get("sample_number", "?")
                    coh = s.get("cohorts_additional", "")
                    print(f"    Ancestría: {anc}  N={n}  {coh}")
                last_results = result

            elif choice == "12":
                info = client.get_release_info()
                print(f"  Release  : {info.get('pgs_catalog_release', '?')}")
                print(f"  Fecha    : {info.get('date', '?')}")
                counts = info.get("counts", {})
                for k, v in counts.items():
                    print(f"  {k:20s}: {v}")
                last_results = info

            elif choice == "13":
                cats = client.get_ancestry_categories()
                print("  Categorías de ancestría:")
                if isinstance(cats, dict):
                    for key, val in cats.items():
                        print(f"    {key}: {val}")
                else:
                    print(f"    {cats}")
                last_results = cats

            elif choice == "14":
                if last_results is None:
                    print("  ⚠ No hay datos. Realiza una consulta primero.")
                    continue
                fname = input("  Archivo de salida (Enter=output.json): ").strip()
                fname = fname or "output.json"
                export_json(last_results, fname)

            elif choice == "15":
                if last_results is None:
                    print("  ⚠ No hay datos. Realiza una consulta primero.")
                    continue
                data_list = last_results if isinstance(last_results, list) else [last_results]
                fname = input("  Archivo de salida (Enter=output.csv): ").strip()
                fname = fname or "output.csv"
                export_csv(data_list, fname)

            else:
                print("  Opción no válida.")

        except Exception as e:
            print(f"  ✗ Error: {e}")
        print()


# ── CLI con argparse ─────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Cliente para la REST API del PGS Catalog.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos:
  python pgs_catalog_client.py score PGS000001
  python pgs_catalog_client.py search-scores --trait EFO_0001360 --limit 5
  python pgs_catalog_client.py search-scores --pmid 25855707
  python pgs_catalog_client.py publication PGP000001
  python pgs_catalog_client.py trait EFO_0000305
  python pgs_catalog_client.py performance --pgs PGS000001
  python pgs_catalog_client.py info
  python pgs_catalog_client.py --interactive
""",
    )
    p.add_argument("--interactive", "-i", action="store_true",
                    help="Abrir menú interactivo")
    p.add_argument("--json", type=str, metavar="FILE",
                    help="Exportar resultado a archivo JSON")
    p.add_argument("--csv", type=str, metavar="FILE",
                    help="Exportar resultado a archivo CSV")

    sub = p.add_subparsers(dest="command")

    # score
    sc = sub.add_parser("score", help="Obtener un score por PGS ID")
    sc.add_argument("pgs_id", help="e.g. PGS000001")

    # search-scores
    ss = sub.add_parser("search-scores", help="Buscar scores por trait/PMID")
    ss.add_argument("--trait", type=str, help="EFO ID (e.g. EFO_0001360)")
    ss.add_argument("--pmid", type=str, help="PubMed ID")
    ss.add_argument("--limit", type=int, default=None)

    # publication
    pb = sub.add_parser("publication", help="Obtener una publicación por PGP ID")
    pb.add_argument("pgp_id", help="e.g. PGP000001")

    # search-publications
    sp = sub.add_parser("search-publications", help="Buscar publicaciones por PMID")
    sp.add_argument("--pmid", type=str, required=True)

    # trait
    tr = sub.add_parser("trait", help="Obtener un trait por EFO ID")
    tr.add_argument("efo_id", help="e.g. EFO_0000305")

    # performance
    pf = sub.add_parser("performance", help="Métricas de rendimiento")
    pf.add_argument("--pgs", type=str, help="PGS ID")
    pf.add_argument("--pmid", type=str, help="PubMed ID")

    # sample-set
    sset = sub.add_parser("sample-set", help="Obtener un sample set por PSS ID")
    sset.add_argument("pss_id", help="e.g. PSS000001")

    # info
    sub.add_parser("info", help="Info del release actual del catálogo")

    # ancestry
    sub.add_parser("ancestry", help="Categorías de ancestría")

    return p


def cli_main():
    parser = build_parser()
    args = parser.parse_args()

    if args.interactive:
        interactive_menu()
        return

    if not args.command:
        parser.print_help()
        return

    client = PGSCatalogClient()
    result: Any = None

    if args.command == "score":
        result = client.get_score(args.pgs_id)
        print_score_summary(result)

    elif args.command == "search-scores":
        result = client.search_scores(trait_id=args.trait, pmid=args.pmid,
                                      limit=args.limit)
        print(f"Encontrados: {len(result)} scores\n")
        for s in result:
            print_score_summary(s)

    elif args.command == "publication":
        result = client.get_publication(args.pgp_id)
        print_publication_summary(result)

    elif args.command == "search-publications":
        result = client.search_publications(pmid=args.pmid)
        print(f"Encontradas: {len(result)} publicaciones\n")
        for p in result:
            print_publication_summary(p)

    elif args.command == "trait":
        result = client.get_trait(args.efo_id)
        print_trait_summary(result)

    elif args.command == "performance":
        result = client.search_performance(pgs_id=args.pgs, pmid=args.pmid)
        print(f"Encontradas: {len(result)} métricas\n")
        print(json.dumps(result, indent=2, ensure_ascii=False))

    elif args.command == "sample-set":
        result = client.get_sample_set(args.pss_id)
        print(json.dumps(result, indent=2, ensure_ascii=False))

    elif args.command == "info":
        result = client.get_release_info()
        print(json.dumps(result, indent=2, ensure_ascii=False))

    elif args.command == "ancestry":
        result = client.get_ancestry_categories()
        print(json.dumps(result, indent=2, ensure_ascii=False))

    # Exportar si se pidió
    if result is not None:
        if args.json:
            export_json(result if isinstance(result, (list, dict)) else [result],
                        args.json)
        if args.csv:
            data_list = result if isinstance(result, list) else [result]
            export_csv(data_list, args.csv)


# ── Entry point ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    cli_main()
