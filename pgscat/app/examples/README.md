# Examples

Minimal dummy datasets for testing and development. **No real cohort data.**

## data/PGS000001/

A minimal PGS000001 directory with synthetic data that matches the expected file layout.

| File | Description |
|------|-------------|
| `PGS000001.json` | Dummy PGS Catalog metadata |
| `PGS000001_PRS_total.tsv` | 10 synthetic samples × 5 chromosomes |
| `PGS000001_PRS_total_metadata.json` | Aggregation report stub |

### Usage

Point the app at the examples directory to run without a real cohort:

```bash
export APP_DATA_DIR=/path/to/app/examples/data
cd src && python3 app.py
```

Then visit `http://localhost:8080/` — PGS000001 should appear in the catalog.

### Convert TSV → Parquet (optional, recommended for dashboards)

```bash
duckdb -c "
  COPY (SELECT * FROM read_csv_auto(
    'examples/data/PGS000001/PGS000001_PRS_total.tsv',
    delim=chr(9), header=true))
  TO 'examples/data/PGS000001/PGS000001_PRS_total.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD)"
```
