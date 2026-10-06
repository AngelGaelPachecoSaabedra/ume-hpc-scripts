import sys
import subprocess
from dask.distributed import Client

def procesar_cromosoma(chrom, pgs_id, script_dir):
    # Usar sys.executable garantiza que estemos en el micromamba vcfarr-hpc
    cmd = [
        sys.executable, f"{script_dir}/compute_prs.py",
        "--pgs-id", pgs_id,
        "--chrom", chrom,
        "--window-size", "10000"
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    
    if res.returncode == 0:
        return f" chr{chrom} OK"
    else:
        # compute_prs.py escupe sus errores en stdout, capturamos ambos:
        error_msg = f" chr{chrom} ERROR (Code {res.returncode}):\n"
        if res.stdout:
            error_msg += f"--- STDOUT ---\n{res.stdout}\n"
        if res.stderr:
            error_msg += f"--- STDERR ---\n{res.stderr}\n"
        return error_msg

if __name__ == "__main__":
    client = Client(sys.argv[1])
    chroms = [str(i) for i in range(1, 23)] + ["X"]
    
    print(f"Repartiendo 23 cromosomas a {len(client.scheduler_info()['workers'])} workers...")
    futuros = client.map(procesar_cromosoma, chroms, pgs_id=sys.argv[2], script_dir=sys.argv[3])
    
    for r in client.gather(futuros):
        print(r)
    client.close()
