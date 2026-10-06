# espirometrias — Matching de 4 fuentes de espirometría

Emparejamiento (*record linkage*) de sujetos entre **4 fuentes** de datos de
espirometría (`op.csv`, `Rcv_trials.csv`, `Session.csv`, `Trials.csv`) con una
estrategia híbrida: match por ID + fecha + talla/peso cuando hay identificador,
y por fecha + sexo + edad + características cuando no lo hay. Se ofrecen varias
implementaciones según la escala de cómputo. Entorno gestionado con `uv`
(`pyproject.toml`).

## Implementaciones (de menor a mayor escala)

| Script | Qué hace |
|--------|----------|
| `matching_simple_v2.py` | Versión secuencial de referencia: busca la mejor coincidencia de cada sujeto por múltiples criterios. |
| `matching_paralelo.py` | Versión paralela con `multiprocessing` (reparte el trabajo entre N CPUs de un nodo). |
| `matching_mpi.py` | Versión **MPI** multinodo (`mpirun -np N`): cada nodo procesa un subconjunto de sujetos. |
| `matching_4fuentes_mpi.py` | Variante MPI que integra las 4 fuentes en un solo cruce. |
| `matching_distribuido.py` | Estrategia híbrida distribuida (prioridad 1: ID+fecha+talla+peso; prioridad 2: fecha+sexo+edad+características). |
| `main.py` | *Entry point* mínimo del paquete. |

## Jobs Slurm

| Script | Lanza |
|--------|-------|
| `job_matching.sh` | `matching_paralelo.py` (un nodo, N CPUs). |
| `job_matching_mpi.sh` | `matching_mpi.py` (4 nodos, MPI). |
| `job_4fuentes.sh` | `matching_4fuentes_mpi.py` (4 nodos). |
| `job_distribuido.sh` | `matching_distribuido.py`. |

> Los `.csv` de entrada contienen datos de sujetos y **no se incluyen** en el repo.
