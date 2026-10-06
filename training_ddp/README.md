# training_ddp — Entrenamiento distribuido PyTorch DDP

Ejemplo de entrenamiento **PyTorch DistributedDataParallel (DDP)** multi-GPU
en un nodo, con entorno gestionado por `uv`.

| Script | Qué hace |
|--------|----------|
| `train_ddp_4gpu_real.py` | Entrenamiento DDP sobre 4 GPUs (`backend='nccl'`): inicializa el grupo de procesos, reparte el modelo/datos y entrena. |
| `job_ddp_uv.sh` | Job que crea/activa el venv `uv` (`.venv`) y lanza el entrenamiento con `torchrun`. |

> El mismo `train_ddp_4gpu_real.py` se usó como plantilla distribuida en el
> matching de `espirometrias/`; aquí vive su copia canónica.
