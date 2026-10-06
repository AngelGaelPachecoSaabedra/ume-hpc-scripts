# Ordanamiento_Bases — Construcción del fenotipo y pipeline PRS

Dos bloques de trabajo:

1. **Fases 1–5 (raíz):** construcción reproducible de la **cohorte
   caso-control oncológica** de MCPS a partir de los diccionarios de datos
   (PDF), las bases clínicas (BASELINE, MORTALITY) y los *updates* genómicos
   (Exoma/WGS QC). Pipeline en `pandas`.
2. **`scripts/`:** pipeline PRS *end-to-end* (cómputo, comparación contra PGS
   Catalog y gráficas en R) usado como back-end de cálculo de scores.

> Las bases de entrada (`.csv`/`.xlsx`/`.pdf` con datos de sujetos) **no se
> incluyen** en el repo. Los scripts referencian sus rutas en CephFS.

## Fases 1–5 (construcción del fenotipo)

| Fase | Script | Qué hace |
|------|--------|----------|
| 1 | `fase1_mort_dict.py` | Extrae texto del diccionario de datos de *mortalidad* (PDF → texto). |
| 1 | `fase1_mortalidad_p2.py` | Extrae las páginas con los códigos ICD-10 de cáncer del *listing* de mortalidad. |
| 1 | `fase1_pdf_mortalidad.py` | Volcado completo del PDF de mortalidad a texto. |
| 1 | `fase1_pdf_basal.py` | Busca términos oncológicos (cáncer, tumor, carcinoma…) en el diccionario *baseline*. |
| 1 | `fase1_pdf_basal_completo.py` | Volcado completo del PDF *baseline* a texto. |
| 2 | `fase2_inspect_excel.py` | Inspecciona las hojas QC Pass de los Excel de Exoma/WGS (columnas, primeras filas). |
| 2 | `fase2_inspect_flags.py` | Inspecciona la distribución de las *MCPS Flags* en Exoma/WGS. |
| 2 | `fase2_denominador_genomico.py` | Construye el **denominador genómico** uniendo QC Pass de Exoma + WGS. |
| 2 | `fase2_denominador_final.py` | Versión final del denominador (con flags de linaje/consistencia). |
| — | `preproceso_genomic_flag.py` | Deriva `GENOMIC_QC_FLAG` / riesgo de *linkage* a partir de `EXOME_MCPS_FLAG`. |
| 3 | `fase3_head_inspection.py` | Inspección de cabeceras y dtypes de BASELINE y MORTALITY. |
| 3 | `fase3_basal_validation.py` | Valida la presencia de las variables oncológicas obligatorias en BASELINE. |
| 3 | `fase3_id_consistency.py` | Chequea la consistencia de `PATID` entre BASELINE, MORTALITY y el denominador. |
| 4 | `fase4_construccion_fenotipo.py` | Construye el **fenotipo oncológico** (caso/control) uniendo denominador + BASELINE + MORTALITY. |
| 4 | `fase4b_cancer_site.py` | Estratifica los casos por **sitio de cáncer** (pulmón, próstata, mama, etc.). |
| 5 | `fase5_exportacion_final.py` | Exporta la cohorte final (`MCPS_Cohorte_Casos_Controles_Oncologicos_F145K.csv`) + metadata. |

El `scripts/` interno tiene su propio README.
