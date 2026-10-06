-- ══════════════════════════════════════════════════════════════════════════
-- PGS Dashboard  –  PostgreSQL schema
-- Target: PostgreSQL 16
-- Apply: psql -U pgs_user -d pgs_dashboard -f sql/schema.sql
-- Idempotent: safe to run multiple times (uses IF NOT EXISTS / OR REPLACE).
-- ══════════════════════════════════════════════════════════════════════════

BEGIN;

-- ── Extensions ────────────────────────────────────────────────────────────
CREATE EXTENSION IF NOT EXISTS "pg_stat_statements";

-- ── Local catalog index ────────────────────────────────────────────────────
-- Indexed mirror of the filesystem scan. Populated by SyncService.
-- Source of truth for the catalog listing page (faster than filesystem scan).
-- Does NOT store scores or sample-level data.
CREATE TABLE IF NOT EXISTS local_pgs_index (
    pgs_id          VARCHAR(12)  PRIMARY KEY,
    trait_name      TEXT,
    trait_efo       JSONB        DEFAULT '[]',
    n_variants      INTEGER,
    has_parquet     BOOLEAN      NOT NULL DEFAULT FALSE,
    has_tsv         BOOLEAN      NOT NULL DEFAULT FALSE,
    has_metadata    BOOLEAN      NOT NULL DEFAULT FALSE,
    chromosomes     JSONB        DEFAULT '[]',  -- list of computed chromosomes ["1","2",...]
    n_chromosomes   SMALLINT     NOT NULL DEFAULT 0,
    first_seen_at   TIMESTAMPTZ  NOT NULL DEFAULT NOW(),
    last_synced_at  TIMESTAMPTZ  NOT NULL DEFAULT NOW(),
    sync_notes      TEXT
);
CREATE INDEX IF NOT EXISTS idx_local_pgs_trait ON local_pgs_index (trait_name);
CREATE INDEX IF NOT EXISTS idx_local_pgs_parquet ON local_pgs_index (has_parquet);
CREATE INDEX IF NOT EXISTS idx_local_pgs_synced ON local_pgs_index (last_synced_at DESC);

-- ── Remote PGS Catalog cache ──────────────────────────────────────────────
-- Caches normalised responses from pgscat/PGS Catalog API.
-- TTL enforced by the application (not DB triggers).
CREATE TABLE IF NOT EXISTS remote_pgs_cache (
    pgs_id              VARCHAR(12)  PRIMARY KEY,
    name                TEXT,
    trait_reported      TEXT,
    trait_efo           JSONB        DEFAULT '[]',
    variants_number     INTEGER,
    is_harmonized       BOOLEAN      DEFAULT FALSE,
    ftp_scoring_file    TEXT,
    publication         JSONB        DEFAULT '{}',
    raw_response        JSONB,       -- full API response for future use
    fetched_at          TIMESTAMPTZ  NOT NULL DEFAULT NOW(),
    ttl_expires_at      TIMESTAMPTZ  NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_remote_cache_expires ON remote_pgs_cache (ttl_expires_at);
CREATE INDEX IF NOT EXISTS idx_remote_cache_trait ON remote_pgs_cache (trait_reported);

-- ── File inventory ────────────────────────────────────────────────────────
-- Per-file inventory for each local PGS directory.
-- Populated by SyncService alongside local_pgs_index.
CREATE TABLE IF NOT EXISTS file_inventory (
    id              SERIAL       PRIMARY KEY,
    pgs_id          VARCHAR(12)  NOT NULL REFERENCES local_pgs_index(pgs_id) ON DELETE CASCADE,
    filename        TEXT         NOT NULL,
    file_type       VARCHAR(20),   -- parquet | tsv_total | tsv_chrom | metadata | betamap | json | other
    chrom           VARCHAR(3),    -- NULL for non-chrom files
    size_bytes      BIGINT,
    last_seen_at    TIMESTAMPTZ  NOT NULL DEFAULT NOW(),
    UNIQUE (pgs_id, filename)
);
CREATE INDEX IF NOT EXISTS idx_inventory_pgs ON file_inventory (pgs_id);
CREATE INDEX IF NOT EXISTS idx_inventory_type ON file_inventory (file_type);

-- ── Computed stats cache ───────────────────────────────────────────────────
-- Persists DuckDB-computed statistics so gunicorn workers survive restarts
-- without re-scanning parquet. Invalidated by source file mtime change.
-- Does NOT store sample-level scores.
CREATE TABLE IF NOT EXISTS stats_cache (
    pgs_id          VARCHAR(12)  PRIMARY KEY REFERENCES local_pgs_index(pgs_id) ON DELETE CASCADE,
    score_column    TEXT         NOT NULL,
    data_format     VARCHAR(10),               -- parquet | tsv
    stats_json      JSONB        NOT NULL,     -- n, mean, stddev, min, max, median, pXX, histogram
    source_mtime    BIGINT,                    -- file mtime (seconds) for invalidation
    computed_at     TIMESTAMPTZ  NOT NULL DEFAULT NOW()
);

-- ── Search audit log ──────────────────────────────────────────────────────
-- Records every remote pgscat search for observability.
-- Useful for understanding usage patterns and debugging pgscat failures.
CREATE TABLE IF NOT EXISTS search_audit (
    id              SERIAL       PRIMARY KEY,
    query           TEXT         NOT NULL,
    search_type     VARCHAR(20),               -- id | trait | pmid | text
    n_results       SMALLINT,
    from_cache      BOOLEAN      NOT NULL DEFAULT FALSE,
    error           TEXT,
    duration_ms     INTEGER,
    searched_at     TIMESTAMPTZ  NOT NULL DEFAULT NOW(),
    client_ip       INET
);
CREATE INDEX IF NOT EXISTS idx_audit_searched ON search_audit (searched_at DESC);
CREATE INDEX IF NOT EXISTS idx_audit_query ON search_audit (query);

-- ── Pipeline plans ────────────────────────────────────────────────────────
-- Stores generated pipeline plans. In Phase 2, these can be picked up
-- by an operator or automated submission agent.
-- The web app only writes; Slurm submission happens externally.
CREATE TABLE IF NOT EXISTS pipeline_plans (
    id              SERIAL       PRIMARY KEY,
    pgs_id          VARCHAR(12)  NOT NULL,
    plan_json       JSONB        NOT NULL,
    manifest_path   TEXT,                      -- path to /work manifest file if written
    status          VARCHAR(20)  NOT NULL DEFAULT 'generated',
                                               -- generated | manifest_written | submitted | done | failed
    created_at      TIMESTAMPTZ  NOT NULL DEFAULT NOW(),
    updated_at      TIMESTAMPTZ  NOT NULL DEFAULT NOW(),
    notes           TEXT
);
CREATE INDEX IF NOT EXISTS idx_plans_pgs ON pipeline_plans (pgs_id);
CREATE INDEX IF NOT EXISTS idx_plans_status ON pipeline_plans (status);

-- ── Catalog sync runs ─────────────────────────────────────────────────────
-- History of filesystem→DB sync operations.
CREATE TABLE IF NOT EXISTS sync_runs (
    id              SERIAL       PRIMARY KEY,
    started_at      TIMESTAMPTZ  NOT NULL DEFAULT NOW(),
    finished_at     TIMESTAMPTZ,
    pgs_scanned     INTEGER,
    pgs_added       INTEGER,
    pgs_updated     INTEGER,
    status          VARCHAR(20)  NOT NULL DEFAULT 'running',  -- running | done | failed
    error           TEXT
);

-- ── Variant annotation jobs ───────────────────────────────────────────────
-- Tracks external annotation pipeline runs (executed via Apptainer on HPC).
-- The web app WRITES job records when users request annotation via UI;
-- the actual pipeline runs externally and writes files to ANNOTATIONS_DIR.
-- Status is updated by polling the filesystem (VariantAnnotationService).
-- No sample-level data is stored here — only job metadata.
CREATE TABLE IF NOT EXISTS annotation_jobs (
    id              SERIAL       PRIMARY KEY,
    pgs_id          VARCHAR(12)  NOT NULL,
    requested_at    TIMESTAMPTZ  NOT NULL DEFAULT NOW(),
    requested_by    INET,                              -- client IP (best-effort)
    started_at      TIMESTAMPTZ,
    finished_at     TIMESTAMPTZ,
    status          VARCHAR(20)  NOT NULL DEFAULT 'requested',
                                                       -- requested | running | done | failed | partial
    gff3_path       TEXT,                              -- GFF3 reference used
    output_dir      TEXT,                              -- filesystem output path
    n_variants      INTEGER,                           -- from summary JSON
    n_coding        INTEGER,
    n_intergenic    INTEGER,
    elapsed_seconds NUMERIC(10,1),
    error           TEXT,
    summary_json    JSONB        DEFAULT '{}'          -- copy of annotation_summary.json stats block
);
CREATE INDEX IF NOT EXISTS idx_annot_jobs_pgs    ON annotation_jobs (pgs_id);
CREATE INDEX IF NOT EXISTS idx_annot_jobs_status ON annotation_jobs (status);
CREATE INDEX IF NOT EXISTS idx_annot_jobs_req    ON annotation_jobs (requested_at DESC);

-- ── Schema version ────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS schema_version (
    version         INTEGER      PRIMARY KEY,
    applied_at      TIMESTAMPTZ  NOT NULL DEFAULT NOW(),
    description     TEXT
);
INSERT INTO schema_version (version, description)
VALUES (1, 'Initial schema: local_pgs_index, remote_pgs_cache, file_inventory, stats_cache, search_audit, pipeline_plans, sync_runs')
ON CONFLICT (version) DO NOTHING;

INSERT INTO schema_version (version, description)
VALUES (2, 'Add annotation_jobs table for Variant Annotator module (Mode D)')
ON CONFLICT (version) DO NOTHING;

COMMIT;
