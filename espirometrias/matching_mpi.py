#!/usr/bin/env python3
"""
MATCHING DISTRIBUIDO CON MPI PARA MÚLTIPLES NODOS
==================================================

Distribuye el trabajo entre múltiples nodos usando MPI.
Cada nodo procesa un subconjunto de sujetos en paralelo.

USO:
    mpirun -np 4 python matching_mpi.py
    
    o con Slurm:
    sbatch job_matching_mpi.sh
"""

from mpi4py import MPI
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
import time
import warnings
import os

warnings.filterwarnings('ignore')

# ============================================================
# CONFIGURACIÓN
# ============================================================

ARCHIVO_SUJETOS = "op.csv"
ARCHIVO_RCV_TRIALS = "Rcv_trials.csv"
ARCHIVO_SESSION = "Session.csv"
ARCHIVO_SALIDA = "resultados_matching.csv"
ARCHIVO_SIN_MATCH = "sin_coincidencia.csv"

# Tolerancias
TOLERANCIA_PESO = 2.0
TOLERANCIA_TALLA = 2.0
TOLERANCIA_EDAD = 2
TOLERANCIA_FECHA = 7

# Puntajes
SCORE_ID_EXACTO = 100
SCORE_PESO = 25
SCORE_TALLA = 25
SCORE_EDAD = 20
SCORE_FECHA_TRIAL = 35
SCORE_SEXO = 10

UMBRAL_MINIMO = 50

# ============================================================
# FUNCIONES DE MATCHING
# ============================================================

def extraer_digitos(valor):
    if pd.isna(valor):
        return None
    resultado = ''.join(c for c in str(valor) if c.isdigit())
    return resultado if resultado else None

def parsear_fecha(fecha_str):
    if pd.isna(fecha_str) or str(fecha_str).strip() == '':
        return None
    fecha_str = str(fecha_str).strip()
    if ' ' in fecha_str:
        fecha_str = fecha_str.split()[0]
    formatos = ['%d/%m/%Y', '%Y-%m-%d', '%d-%m-%Y', '%m/%d/%Y', '%Y/%m/%d']
    for fmt in formatos:
        try:
            return pd.to_datetime(fecha_str, format=fmt)
        except:
            continue
    try:
        return pd.to_datetime(fecha_str, dayfirst=True)
    except:
        return None

def normalizar_sexo(valor):
    if pd.isna(valor):
        return None
    valor = str(valor).strip().upper()
    if valor in ['M', 'MASCULINO', 'MALE', '1', 'H', 'HOMBRE']:
        return 'M'
    elif valor in ['F', 'FEMENINO', 'FEMALE', '2', 'MUJER']:
        return 'F'
    return valor

def calcular_score(sujeto, trial):
    score = 0
    detalles = []
    
    id_sujeto = extraer_digitos(sujeto.get('p_08070803'))
    id_trial = extraer_digitos(trial.get('id'))
    if id_sujeto and id_trial and id_sujeto == id_trial:
        score += SCORE_ID_EXACTO
        detalles.append("ID")
    
    try:
        peso_suj = float(sujeto.get('p_080612') or 0)
        peso_trial = float(trial.get('weight') or 0)
        if peso_suj > 0 and peso_trial > 0:
            diff = abs(peso_suj - peso_trial)
            if diff <= TOLERANCIA_PESO:
                score += SCORE_PESO * (1 - diff/TOLERANCIA_PESO)
                detalles.append("Peso")
    except:
        pass
    
    try:
        talla_suj = float(sujeto.get('p_080610') or 0)
        talla_trial = float(trial.get('height') or 0)
        if talla_suj > 0 and talla_trial > 0:
            diff = abs(talla_suj - talla_trial)
            if diff <= TOLERANCIA_TALLA:
                score += SCORE_TALLA * (1 - diff/TOLERANCIA_TALLA)
                detalles.append("Talla")
    except:
        pass
    
    try:
        edad_suj = float(sujeto.get('edad_hoy') or sujeto.get('edad') or 0)
        edad_trial = float(trial.get('age') or 0)
        if edad_suj > 0 and edad_trial > 0:
            diff = abs(edad_suj - edad_trial)
            if diff <= TOLERANCIA_EDAD:
                score += SCORE_EDAD * (1 - diff/TOLERANCIA_EDAD)
                detalles.append("Edad")
    except:
        pass
    
    try:
        fecha_suj = parsear_fecha(sujeto.get('fecha'))
        fecha_trial_date = parsear_fecha(trial.get('trialdate'))
        if fecha_suj and fecha_trial_date:
            diff_dias = abs((fecha_trial_date - fecha_suj).days)
            if diff_dias <= TOLERANCIA_FECHA:
                score += SCORE_FECHA_TRIAL * (1 - diff_dias/TOLERANCIA_FECHA)
                detalles.append("Fecha")
    except:
        pass
    
    sexo_suj = normalizar_sexo(sujeto.get('sexo'))
    sexo_trial = normalizar_sexo(trial.get('sex'))
    if sexo_suj and sexo_trial and sexo_suj == sexo_trial:
        score += SCORE_SEXO
        detalles.append("Sexo")
    
    return score, detalles

def determinar_confianza(score, detalles):
    tiene_id = "ID" in detalles
    otros = len([d for d in detalles if d != "ID"])
    if tiene_id and otros >= 3:
        return "MUY ALTA"
    elif tiene_id:
        return "ALTA"
    elif score >= 80:
        return "MEDIA-ALTA"
    elif score >= 60:
        return "MEDIA"
    return "BAJA"

def procesar_sujeto(args):
    """Procesa un sujeto contra todos los trials."""
    sujeto, trials_list, session_dict = args
    nombre = f"{sujeto.get('nombre', '')} {sujeto.get('paterno', '')} {sujeto.get('materno', '')}".strip()
    
    mejor_score = 0
    mejor_trial = None
    mejor_detalles = []
    
    for trial in trials_list:
        score, detalles = calcular_score(sujeto, trial)
        if score > mejor_score:
            mejor_score = score
            mejor_trial = trial
            mejor_detalles = detalles
    
    if mejor_score >= UMBRAL_MINIMO and mejor_trial:
        confianza = determinar_confianza(mejor_score, mejor_detalles)
        session_data = session_dict.get(mejor_trial.get('id_session'), {})
        
        return ('match', {
            "registro": sujeto.get('registro', ''),
            "nombre_completo": nombre,
            "nombre": sujeto.get('nombre', ''),
            "paterno": sujeto.get('paterno', ''),
            "materno": sujeto.get('materno', ''),
            "id_sujeto": extraer_digitos(sujeto.get('p_08070803')),
            "fecha_encuesta": sujeto.get('fecha', ''),
            "peso_sujeto": sujeto.get('p_080612', ''),
            "talla_sujeto": sujeto.get('p_080610', ''),
            "edad_sujeto": sujeto.get('edad_hoy', ''),
            "sexo_sujeto": sujeto.get('sexo', ''),
            "id_trial": mejor_trial.get('id', ''),
            "id_session": mejor_trial.get('id_session', ''),
            "birthdate_trial": mejor_trial.get('birthdate', ''),
            "trialdate": mejor_trial.get('trialdate', ''),
            "peso_trial": mejor_trial.get('weight', ''),
            "talla_trial": mejor_trial.get('height', ''),
            "edad_trial": mejor_trial.get('age', ''),
            "sexo_trial": mejor_trial.get('sex', ''),
            "fvcbest": session_data.get('fvcbest', ''),
            "fev1best": session_data.get('fev1best', ''),
            "pefbest": session_data.get('pefbest', ''),
            "match_score": round(mejor_score, 1),
            "confianza": confianza,
            "criterios_match": "|".join(mejor_detalles),
            "num_criterios": len(mejor_detalles),
            "municipio": sujeto.get('municipio', ''),
        })
    else:
        return ('sin_match', {
            "registro": sujeto.get('registro', ''),
            "nombre_completo": nombre,
            "id_sujeto": extraer_digitos(sujeto.get('p_08070803')),
            "mejor_score": round(mejor_score, 1),
        })

def procesar_chunk_local(sujetos_chunk, trials_list, session_dict, n_threads):
    """Procesa un chunk de sujetos usando multiprocessing local."""
    args_list = [(s, trials_list, session_dict) for s in sujetos_chunk]
    
    with Pool(processes=n_threads) as pool:
        resultados = pool.map(procesar_sujeto, args_list)
    
    matches = [r[1] for r in resultados if r[0] == 'match']
    sin_match = [r[1] for r in resultados if r[0] == 'sin_match']
    
    return matches, sin_match

# ============================================================
# PROGRAMA PRINCIPAL MPI
# ============================================================

def main():
    # Inicializar MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    # Obtener CPUs por nodo
    cpus_por_nodo = int(os.environ.get('SLURM_CPUS_PER_TASK', cpu_count()))
    
    if rank == 0:
        print("="*60)
        print("MATCHING DISTRIBUIDO MPI")
        print("="*60)
        print(f"Nodos MPI: {size}")
        print(f"CPUs por nodo: {cpus_por_nodo}")
        print(f"Total hilos: {size * cpus_por_nodo}")
        print("="*60)
    
    # Rank 0 carga los datos y los distribuye
    if rank == 0:
        print(f"\n[Rank 0] Cargando datos...")
        inicio = time.time()
        
        # Cargar archivos
        df_sujetos = pd.read_csv(ARCHIVO_SUJETOS, dtype=str, low_memory=False)
        print(f"  ✓ Sujetos: {len(df_sujetos):,}")
        
        df_rcv = pd.read_csv(ARCHIVO_RCV_TRIALS, low_memory=False)
        df_rcv = df_rcv.drop_duplicates(subset=['id', 'birthdate'], keep='first')
        print(f"  ✓ Trials: {len(df_rcv):,}")
        
        try:
            df_session = pd.read_csv(ARCHIVO_SESSION, low_memory=False)
            session_dict = df_session.set_index('id_session').to_dict('index')
        except:
            session_dict = {}
        
        # Filtrar pendientes
        col_loc = None
        for col in ['Localizado', 'localizado', 'verificado']:
            if col in df_sujetos.columns:
                col_loc = col
                break
        
        if col_loc:
            df_pendientes = df_sujetos[
                df_sujetos[col_loc].fillna('').str.strip().str.lower() != 'si'
            ]
        else:
            df_pendientes = df_sujetos
        
        sujetos_list = df_pendientes.to_dict('records')
        trials_list = df_rcv.to_dict('records')
        
        print(f"  → {len(sujetos_list):,} sujetos a procesar")
        
        # Dividir sujetos entre nodos
        chunks = np.array_split(sujetos_list, size)
        chunks = [list(c) for c in chunks]
        
        print(f"\n[Rank 0] Distribuyendo trabajo...")
        for i, chunk in enumerate(chunks):
            print(f"  Nodo {i}: {len(chunk):,} sujetos")
    else:
        chunks = None
        trials_list = None
        session_dict = None
    
    # Distribuir datos a todos los nodos
    if rank == 0:
        mi_chunk = chunks[0]
        for i in range(1, size):
            comm.send(chunks[i], dest=i, tag=1)
            comm.send(trials_list, dest=i, tag=2)
            comm.send(session_dict, dest=i, tag=3)
    else:
        mi_chunk = comm.recv(source=0, tag=1)
        trials_list = comm.recv(source=0, tag=2)
        session_dict = comm.recv(source=0, tag=3)
    
    # Sincronizar antes de procesar
    comm.Barrier()
    
    if rank == 0:
        print(f"\n[Todos] Procesando...")
        inicio_proceso = time.time()
    
    # Cada nodo procesa su chunk con multiprocessing local
    mis_matches, mis_sin_match = procesar_chunk_local(
        mi_chunk, trials_list, session_dict, cpus_por_nodo
    )
    
    print(f"  [Rank {rank}] Terminado: {len(mis_matches)} matches, {len(mis_sin_match)} sin match")
    
    # Recolectar resultados en rank 0
    comm.Barrier()
    
    todos_matches = comm.gather(mis_matches, root=0)
    todos_sin_match = comm.gather(mis_sin_match, root=0)
    
    # Rank 0 guarda los resultados
    if rank == 0:
        tiempo_proceso = time.time() - inicio_proceso
        
        # Aplanar listas
        matches_final = [item for sublist in todos_matches for item in sublist]
        sin_match_final = [item for sublist in todos_sin_match for item in sublist]
        
        print(f"\n[Rank 0] Guardando resultados...")
        
        if matches_final:
            df_res = pd.DataFrame(matches_final)
            df_res = df_res.sort_values('match_score', ascending=False)
            df_res.to_csv(ARCHIVO_SALIDA, index=False, encoding='utf-8-sig')
            print(f"  ✓ {ARCHIVO_SALIDA}: {len(df_res):,} coincidencias")
        
        if sin_match_final:
            df_sin = pd.DataFrame(sin_match_final)
            df_sin.to_csv(ARCHIVO_SIN_MATCH, index=False, encoding='utf-8-sig')
            print(f"  ✓ {ARCHIVO_SIN_MATCH}: {len(df_sin):,} sin match")
        
        tiempo_total = time.time() - inicio
        total_sujetos = len(sujetos_list)
        
        print("\n" + "="*60)
        print("RESUMEN")
        print("="*60)
        print(f"Sujetos procesados: {total_sujetos:,}")
        print(f"Coincidencias: {len(matches_final):,} ({len(matches_final)/total_sujetos*100:.1f}%)")
        print(f"Sin match: {len(sin_match_final):,}")
        print(f"Tiempo proceso: {tiempo_proceso:.1f}s")
        print(f"Tiempo total: {tiempo_total:.1f}s")
        print(f"Velocidad: {total_sujetos/tiempo_proceso:.0f} sujetos/seg")
        print("="*60)

if __name__ == "__main__":
    main()
