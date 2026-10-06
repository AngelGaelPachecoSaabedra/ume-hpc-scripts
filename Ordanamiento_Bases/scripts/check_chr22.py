import zarr
import sys

zarr_path = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zar_files/chr22.zarr"
print(f"Abriendo {zarr_path}...")
store = zarr.open_group(zarr_path, mode='r')
ds = store["call_DS"]

print(f"Iniciando escaneo de {ds.shape[0]} variantes...")
chunk_size = ds.chunks[0]

for start_idx in range(0, ds.shape[0], chunk_size):
    end_idx = min(start_idx + chunk_size, ds.shape[0])
    try:
        # Extraer el bloque fuerza a Blosc a descomprimirlo
        _ = ds[start_idx:end_idx]
    except Exception as e:
        print(f"\n[¡ERROR DETECTADO!] Chunk corrupto en las variantes: {start_idx} a {end_idx}")
        print(f"Detalle: {e}")
        sys.exit(1)

print("\nEscaneo completado. Todos los chunks están íntegros.")
