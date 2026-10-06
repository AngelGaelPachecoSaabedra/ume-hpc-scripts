#!/usr/bin/env python3
"""
MATCHING 4 FUENTES - VERSIÓN MPI DISTRIBUIDA
============================================
"""

from mpi4py import MPI
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
import time
import warnings
import os

warnings.filterwarnings('ignore')

# Archivos
ARCHIVO_OP = "op.csv"
ARCHIVO_RCV = "Rcv_trials.csv"
ARCHIVO_SESSION = "Session.csv"
ARCHIVO_TRIALS = "Trials.csv"
ARCHIVO_SALIDA = "matching_4fuentes_completo.csv"

# Tolerancias y puntajes
TOL_PESO, TOL_TALLA, TOL_EDAD, TOL_FECHA = 2.0, 2.0, 2, 7
SCORE_ID, SCORE_FECHA, SCORE_TALLA, SCORE_PESO, SCORE_EDAD, SCORE_SEXO = 100, 40, 30, 30, 25, 15
UMBRAL = 60

# Columnas
COLS_OP = ['secuencia', 'registro', 'cadena', 'municipio', 'ageb', 'area', 'manzana', 'vivienda',
           'nombre', 'paterno', 'materno', 'sexo', 'edad', 'edad_hoy', 'fecha', 'hora_ini',
           'p_0101', 't_0101', 'p_0102', 't_0102', 'p_0103', 't_0103',
           'p_080608', 'p_080609', 'p_080610', 'p_080612', 'p_08070802', 'p_08070803',
           'p_0811', 't_0811', 'p_0815', 't_0815', 'p72', 'vive', 'tableta', 'encuesto', 'verificado']

COLS_SESSION = ['id_session', 'id_visit', 'sessiondat', 'sessiontim', 'reproducib', 'besttriali',
                'bestfvcfev', 'fvcbest', 'fev1best', 'pefbest', 'vcbest', 'mvvbest',
                'bestfvc', 'bestfev1', 'bestpef', 'bestrestri', 'bestobstru',
                'predictedf', 'predictedv', 'predictedm', 'fef2575_pc', 'fev1_xvc',
                'id_doctor', 'id_operato', 'sessionnot']

COLS_TRIALS = ['id_trial', 'id_session', 'id_trialty', 'id_device', 'trialdate', 'trialtime',
               'trialtempe', 'qualitycod', 'recordqual', 'bestparamv', 'bestfev1', 'bestvc',
               'bestmvv', 'prepost', 'activate', 'fvc_corr', 'fivc_corr']

COLS_RCV = ['id_rcv', 'id_trial', 'id_session', 'id', 'subjectcod', 'subjectsur', 'subjectnam',
            'birthdate', 'sex', 'age', 'height', 'weight', 'id_ethnicg', 'id_smokety',
            'smoker', 'smokeqty', 'trialdate', 'trialtime', 'qualitycod', 'recordqual',
            'symptoms', 'dyspnea_in', 'dyspnea_en', 'fatigue_in', 'fatigue_en']

# Funciones auxiliares
def extraer_digitos(v):
    if pd.isna(v): return None
    r = ''.join(c for c in str(v) if c.isdigit())
    return r if r else None

def parsear_fecha(f):
    if pd.isna(f) or str(f).strip() == '': return None
    f = str(f).strip().split()[0]
    for fmt in ['%d/%m/%Y', '%Y-%m-%d', '%d-%m-%Y', '%m/%d/%Y']:
        try: return pd.to_datetime(f, format=fmt)
        except: continue
    try: return pd.to_datetime(f, dayfirst=True)
    except: return None

def normalizar_sexo(v):
    if pd.isna(v): return None
    v = str(v).strip().upper()
    if v in ['M', 'MASCULINO', 'MALE', '1', 'H', 'HOMBRE']: return 'M'
    if v in ['F', 'FEMENINO', 'FEMALE', '2', 'MUJER']: return 'F'
    return v

def to_float(v):
    try: return float(v) if not pd.isna(v) else None
    except: return None

def calcular_score(sujeto, rcv):
    score = 0
    det = {'id': False, 'fecha': False, 'talla': False, 'peso': False, 'edad': False, 'sexo': False}
    crit, diffs = [], {}
    
    id_s, id_r = extraer_digitos(sujeto.get('p_08070803')), extraer_digitos(rcv.get('id'))
    if id_s and id_r and id_s == id_r:
        score += SCORE_ID; det['id'] = True; crit.append("ID")
    
    f_s, f_r = parsear_fecha(sujeto.get('fecha')), parsear_fecha(rcv.get('trialdate'))
    if f_s and f_r:
        d = abs((f_r - f_s).days); diffs['fecha'] = d
        if d <= TOL_FECHA:
            score += SCORE_FECHA * (1 - d/TOL_FECHA); det['fecha'] = True; crit.append(f"Fecha±{d}d")
    
    t_s, t_r = to_float(sujeto.get('p_080610')), to_float(rcv.get('height'))
    if t_s and t_r and t_s > 0 and t_r > 0:
        d = abs(t_s - t_r); diffs['talla'] = round(d, 1)
        if d <= TOL_TALLA:
            score += SCORE_TALLA * (1 - d/TOL_TALLA); det['talla'] = True; crit.append(f"Talla±{d:.1f}")
    
    p_s, p_r = to_float(sujeto.get('p_080612')), to_float(rcv.get('weight'))
    if p_s and p_r and p_s > 0 and p_r > 0:
        d = abs(p_s - p_r); diffs['peso'] = round(d, 1)
        if d <= TOL_PESO:
            score += SCORE_PESO * (1 - d/TOL_PESO); det['peso'] = True; crit.append(f"Peso±{d:.1f}")
    
    e_s, e_r = to_float(sujeto.get('edad_hoy') or sujeto.get('edad')), to_float(rcv.get('age'))
    if e_s and e_r and e_s > 0 and e_r > 0:
        d = abs(e_s - e_r); diffs['edad'] = int(d)
        if d <= TOL_EDAD:
            score += SCORE_EDAD * (1 - d/TOL_EDAD); det['edad'] = True; crit.append(f"Edad±{int(d)}")
    
    s_s = normalizar_sexo(sujeto.get('sexo') or sujeto.get('p_080608'))
    s_r = normalizar_sexo(rcv.get('sex'))
    if s_s and s_r and s_s == s_r:
        score += SCORE_SEXO; det['sexo'] = True; crit.append("Sexo")
    
    return round(score, 1), det, crit, diffs

def get_confianza(score, det):
    n = sum(det.values())
    if det['id']:
        if n >= 5: return "CERTEZA"
        if n >= 3: return "MUY_ALTA"
        return "ALTA"
    if score >= 100: return "MEDIA_ALTA"
    if score >= 80: return "MEDIA"
    if score >= 60: return "BAJA"
    return "MUY_BAJA"

def procesar_sujeto(args):
    sujeto, rcv_list, session_idx, trials_idx = args
    mejor = (0, None, None, None, None)
    
    for rcv in rcv_list:
        score, det, crit, diffs = calcular_score(sujeto, rcv)
        if score > mejor[0]:
            mejor = (score, rcv, det, crit, diffs)
    
    return sujeto, mejor, session_idx, trials_idx

def procesar_chunk(chunk, rcv_list, session_idx, trials_idx, n_threads):
    args_list = [(s, rcv_list, session_idx, trials_idx) for s in chunk]
    with Pool(n_threads) as pool:
        return pool.map(procesar_sujeto, args_list)

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', cpu_count()))
    
    if rank == 0:
        print("="*60)
        print("MATCHING 4 FUENTES - MPI DISTRIBUIDO")
        print("="*60)
        print(f"Nodos MPI: {size}, CPUs/nodo: {cpus}, Total: {size*cpus}")
        print("="*60)
        
        # Cargar datos
        print("\n[1/4] Cargando datos...")
        df_op = pd.read_csv(ARCHIVO_OP, dtype=str, low_memory=False)
        print(f"  ✓ op.csv: {len(df_op):,}")
        
        df_rcv = pd.read_csv(ARCHIVO_RCV, low_memory=False)
        df_rcv = df_rcv.drop_duplicates(subset=['id', 'birthdate'], keep='first')
        print(f"  ✓ Rcv_trials.csv: {len(df_rcv):,}")
        
        session_idx, trials_idx = {}, {}
        try:
            df_ses = pd.read_csv(ARCHIVO_SESSION, low_memory=False)
            for _, r in df_ses.iterrows(): session_idx[r.get('id_session')] = r.to_dict()
            print(f"  ✓ Session.csv: {len(df_ses):,}")
        except: print("  ⚠ Session.csv no encontrado")
        
        try:
            df_tri = pd.read_csv(ARCHIVO_TRIALS, low_memory=False)
            for _, r in df_tri.iterrows(): trials_idx[(r.get('id_session'), r.get('id_trial'))] = r.to_dict()
            print(f"  ✓ Trials.csv: {len(df_tri):,}")
        except: print("  ⚠ Trials.csv no encontrado")
        
        sujetos = df_op.to_dict('records')
        rcv_list = df_rcv.to_dict('records')
        
        chunks = [list(c) for c in np.array_split(sujetos, size)]
        print(f"\n[2/4] Distribuyendo {len(sujetos):,} sujetos en {size} nodos...")
    else:
        chunks, rcv_list, session_idx, trials_idx = None, None, None, None
    
    # Distribuir
    if rank == 0:
        mi_chunk = chunks[0]
        for i in range(1, size):
            comm.send((chunks[i], rcv_list, session_idx, trials_idx), dest=i)
    else:
        mi_chunk, rcv_list, session_idx, trials_idx = comm.recv(source=0)
    
    comm.Barrier()
    
    if rank == 0:
        print(f"\n[3/4] Procesando...")
        inicio = time.time()
    
    # Procesar
    mis_resultados = procesar_chunk(mi_chunk, rcv_list, session_idx, trials_idx, cpus)
    print(f"  [Rank {rank}] {len(mis_resultados)} procesados")
    
    comm.Barrier()
    todos = comm.gather(mis_resultados, root=0)
    
    if rank == 0:
        tiempo = time.time() - inicio
        resultados = [r for sub in todos for r in sub]
        
        print(f"\n[4/4] Generando salida...")
        
        filas = []
        stats = {'certeza': 0, 'muy_alta': 0, 'alta': 0, 'media_alta': 0, 'media': 0, 'baja': 0, 'sin_match': 0}
        
        for sujeto, (score, rcv, det, crit, diffs), ses_idx, tri_idx in resultados:
            fila = {}
            for col in COLS_OP:
                if col in sujeto: fila[f'op_{col}'] = sujeto.get(col, '')
            
            if score >= UMBRAL and rcv:
                conf = get_confianza(score, det)
                if conf == 'CERTEZA': stats['certeza'] += 1
                elif conf == 'MUY_ALTA': stats['muy_alta'] += 1
                elif conf == 'ALTA': stats['alta'] += 1
                elif conf == 'MEDIA_ALTA': stats['media_alta'] += 1
                elif conf == 'MEDIA': stats['media'] += 1
                else: stats['baja'] += 1
                
                fila['MATCH'] = 'SI'
                fila['match_score'] = score
                fila['match_confianza'] = conf
                fila['match_criterios'] = '|'.join(crit)
                fila['match_n_criterios'] = len(crit)
                for k in ['id', 'fecha', 'talla', 'peso', 'edad', 'sexo']:
                    fila[f'match_{k}'] = 'SI' if det[k] else 'NO'
                for k, v in diffs.items():
                    fila[f'diff_{k}'] = v
                
                for col in COLS_RCV:
                    fila[f'rcv_{col}'] = rcv.get(col, '')
                
                id_ses = rcv.get('id_session')
                if id_ses in ses_idx:
                    for col in COLS_SESSION:
                        fila[f'ses_{col}'] = ses_idx[id_ses].get(col, '')
                
                key = (id_ses, rcv.get('id_trial'))
                if key in tri_idx:
                    for col in COLS_TRIALS:
                        if col != 'id_session':
                            fila[f'tri_{col}'] = tri_idx[key].get(col, '')
            else:
                stats['sin_match'] += 1
                fila['MATCH'] = 'NO'
                fila['match_score'] = score if score > 0 else 0
                fila['match_confianza'] = 'SIN_MATCH'
            
            filas.append(fila)
        
        df_out = pd.DataFrame(filas)
        df_out = df_out.sort_values(['MATCH', 'match_score'], ascending=[False, False])
        df_out.to_csv(ARCHIVO_SALIDA, index=False, encoding='utf-8-sig')
        
        n_match = stats['certeza'] + stats['muy_alta'] + stats['alta'] + stats['media_alta'] + stats['media'] + stats['baja']
        
        print(f"\n  ✓ {ARCHIVO_SALIDA}: {len(df_out):,} filas, {len(df_out.columns)} columnas")
        print("\n" + "="*60)
        print("RESUMEN")
        print("="*60)
        print(f"Total: {len(resultados):,} | Match: {n_match:,} | Sin: {stats['sin_match']:,}")
        print(f"CERTEZA: {stats['certeza']} | MUY_ALTA: {stats['muy_alta']} | ALTA: {stats['alta']}")
        print(f"MEDIA_ALTA: {stats['media_alta']} | MEDIA: {stats['media']} | BAJA: {stats['baja']}")
        print(f"Tiempo: {tiempo:.1f}s | Velocidad: {len(resultados)/tiempo:.0f}/seg")
        print("="*60)

if __name__ == "__main__":
    main()
