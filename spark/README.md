# spark — Utilidades Apache Spark sobre Slurm

Scripts de cómputo distribuido con **Apache Spark** (instalación en CephFS,
`SPARK_HOME=/mnt/cephfs/biocontainers/spark/spark`) lanzado sobre Slurm.

| Script | Qué hace |
|--------|----------|
| `match_cadd.py` | Cruza el catálogo de variantes del paciente contra la base masiva **CADD** (`whole_genome_SNVs.tsv.gz`) para anotar *scores* PHRED por variante. |
| `match_prs.py` | Cruza el catálogo de variantes del Zarr contra los pesos PGS limpios (coincidencia `CHR,POS,EA,OA`) para el PRS de diabetes. |
| `spark_slurm_pi.sh` | Ejemplo mínimo: levanta un clúster Spark multinodo sobre Slurm y corre el cálculo de π (prueba de humo de la integración Spark↔Slurm). |
