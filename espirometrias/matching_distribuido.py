#!/usr/bin/env python3
"""
MATCHING DE ESPIROMETRÍAS - ESTRATEGIA HÍBRIDA
===============================================

PRIORIDAD 1: Si tiene p_08070803, usar ID + fecha + talla + peso
PRIORIDAD 2: Si no tiene p_08070803, usar fecha + sexo + edad + talla/peso

Meta: ~7,000+ matches (3,300 por ID + 4,700 por características)
"""

import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
import argparse
from pathlib import Path
from datetime import datetime
import os
import time
import math

BASE_DIR = "/mnt/cephfs/orgs/home/angel.pacheco/espirometrias"
TEMP_DIR = os.path.join(BASE_DIR, "temp_chunks")

# Variables globales
rcv_by_date = None
rcv_by_id = None
session_dict = None


def init_worker(rcv_date_data, rcv_id_data, ses_data):
    global rcv_by_date, rcv_by_id, session_dict
    rcv_by_date = rcv_date_data
    rcv_by_id = rcv_id_data
    session_dict = ses_data


def parse_fecha_solo_dia(fecha_str):
    if pd.isna(fecha_str):
        return None
    fecha_str = str(fecha_str).strip()
    if ' ' in fecha_str:
        fecha_str = fecha_str.split()[0]
    return fecha_str


def fecha_to_datetime(fecha_str):
    if not fecha_str:
        return None
    try:
        for fmt in ['%d/%m/%Y', '%Y-%m-%d', '%d-%m-%Y']:
            try:
                return datetime.strptime(fecha_str, fmt)
            except:
                continue
        return None
    except:
        return None


def extraer_año(fecha_str):
    try:
        return int(str(fecha_str).split('/')[2].split()[0])
    except:
        return None


def safe_int(val):
    if pd.isna(val) or val == '':
        return None
    try:
        return int(float(val))
    except:
        return None


def convertir_sexo_op(sexo_op):
    """Convierte sexo op (1=H, 2=M) a rcv (1=H, 0=M)"""
    if sexo_op == 1 or sexo_op == '1':
        return 1
    elif sexo_op == 2 or sexo_op == '2':
        return 0
    return None


def buscar_por_id(op_row):
    """
    PRIORIDAD 1: Buscar usando p_08070803 (ID ya asignado manualmente)
    """
    global rcv_by_id
    
    op_id = op_row.get('p_08070803')
    if not op_id or pd.isna(op_id) or str(op_id).strip() == '':
        return None, 0, [], {}, 0
    
    op_id = str(op_id).strip()
    
    if op_id not in rcv_by_id:
        return None, 0, [], {}, 0
    
    candidatos = rcv_by_id[op_id]
    
    op_fecha = op_row.get('fecha_dia')
    op_talla = op_row.get('talla_int')
    op_peso = op_row.get('peso_int')
    op_edad = op_row.get('edad_int')
    
    mejor_match = None
    mejor_score = 0
    mejor_criterios = ['ID_DIRECTO']
    mejor_diffs = {}
    
    for cand in candidatos:
        score = 100  # Base alta por tener ID
        criterios = ['ID_DIRECTO']
        diffs = {}
        
        # Fecha exacta
        rcv_fecha = cand.get('trialdate')
        if op_fecha and rcv_fecha:
            if op_fecha == rcv_fecha:
                score += 50
                criterios.append('FECHA_EXACTA')
                diffs['fecha'] = 0
            else:
                # Intentar fecha cercana
                op_dt = op_row.get('fecha_dt')
                rcv_dt = cand.get('trialdate_dt')
                if op_dt and rcv_dt:
                    diff_dias = abs((rcv_dt - op_dt).days)
                    if diff_dias <= 5:
                        score += 30 - (diff_dias * 5)
                        criterios.append(f'FECHA_±{diff_dias}d')
                        diffs['fecha'] = diff_dias
                    else:
                        continue  # Fecha muy diferente, saltar
        
        # Talla
        rcv_talla = cand.get('height_int')
        if op_talla is not None and rcv_talla is not None:
            diff = abs(rcv_talla - op_talla)
            diffs['talla'] = diff
            if diff == 0:
                score += 30
                criterios.append('TALLA_EXACTA')
            elif diff <= 1:
                score += 20
                criterios.append('TALLA_±1')
            elif diff <= 2:
                score += 10
                criterios.append('TALLA_±2')
            else:
                score -= 20  # Penalizar diferencia grande
        
        # Peso
        rcv_peso = cand.get('weight_int')
        if op_peso is not None and rcv_peso is not None:
            diff = abs(rcv_peso - op_peso)
            diffs['peso'] = diff
            if diff == 0:
                score += 30
                criterios.append('PESO_EXACTO')
            elif diff <= 1:
                score += 20
                criterios.append('PESO_±1')
            elif diff <= 2:
                score += 10
                criterios.append('PESO_±2')
            else:
                score -= 20
        
        # Edad
        rcv_edad = cand.get('age_int')
        if op_edad is not None and rcv_edad is not None:
            diff = abs(rcv_edad - op_edad)
            diffs['edad'] = diff
            if diff == 0:
                score += 20
                criterios.append('EDAD_EXACTA')
            elif diff <= 1:
                score += 10
                criterios.append('EDAD_±1')
            elif diff <= 2:
                score += 5
                criterios.append('EDAD_±2')
        
        if score > mejor_score:
            mejor_score = score
            mejor_match = cand
            mejor_criterios = criterios
            mejor_diffs = diffs
    
    return mejor_match, mejor_score, mejor_criterios, mejor_diffs, len(candidatos)


def buscar_por_caracteristicas(op_row):
    """
    PRIORIDAD 2: Buscar por fecha + sexo + edad + talla/peso
    """
    global rcv_by_date
    
    op_fecha_dt = op_row.get('fecha_dt')
    
    if not op_fecha_dt:
        return None, 0, [], {}, 0
    
    # Buscar candidatos en fechas cercanas (0 a +5 días)
    candidatos = []
    
    for delta in range(0, 6):
        fecha_buscar = op_fecha_dt + pd.Timedelta(days=delta)
        fecha_str = fecha_buscar.strftime('%d/%m/%Y')
        
        if fecha_str in rcv_by_date:
            candidatos.extend(rcv_by_date[fecha_str])
    
    if not candidatos:
        return None, 0, [], {}, 0
    
    # Pre-filtrar por sexo
    op_sexo = convertir_sexo_op(op_row.get('sexo'))
    if op_sexo is not None:
        candidatos = [c for c in candidatos if safe_int(c.get('sex')) == op_sexo]
    
    if not candidatos:
        return None, 0, [], {}, 0
    
    op_talla = op_row.get('talla_int')
    op_peso = op_row.get('peso_int')
    op_edad = op_row.get('edad_int')
    op_año = op_row.get('año_nac')
    
    mejor_match = None
    mejor_score = 0
    mejor_criterios = []
    mejor_diffs = {}
    
    for cand in candidatos:
        score = 10  # Base por sexo
        criterios = ['SEXO']
        diffs = {}
        
        # Fecha
        rcv_fecha_dt = cand.get('trialdate_dt')
        if op_fecha_dt and rcv_fecha_dt:
            diff_dias = (rcv_fecha_dt - op_fecha_dt).days
            diffs['fecha'] = diff_dias
            
            if diff_dias == 0:
                score += 30
                criterios.append('FECHA_EXACTA')
            elif diff_dias <= 2:
                score += 20
                criterios.append(f'FECHA_+{diff_dias}d')
            else:
                score += 10
                criterios.append(f'FECHA_+{diff_dias}d')
        
        # Edad
        rcv_edad = cand.get('age_int')
        if op_edad is not None and rcv_edad is not None and rcv_edad > 0:
            diff = abs(rcv_edad - op_edad)
            diffs['edad'] = diff
            
            if diff == 0:
                score += 30
                criterios.append('EDAD_EXACTA')
            elif diff == 1:
                score += 20
                criterios.append('EDAD_±1')
            elif diff == 2:
                score += 10
                criterios.append('EDAD_±2')
            elif diff > 3:
                continue  # Saltar si edad muy diferente
        
        # Talla
        rcv_talla = cand.get('height_int')
        if op_talla is not None and rcv_talla is not None:
            diff = abs(rcv_talla - op_talla)
            diffs['talla'] = diff
            
            if diff == 0:
                score += 25
                criterios.append('TALLA_EXACTA')
            elif diff == 1:
                score += 15
                criterios.append('TALLA_±1')
            elif diff == 2:
                score += 8
                criterios.append('TALLA_±2')
            elif diff > 3:
                score -= 10
        
        # Peso
        rcv_peso = cand.get('weight_int')
        if op_peso is not None and rcv_peso is not None:
            diff = abs(rcv_peso - op_peso)
            diffs['peso'] = diff
            
            if diff == 0:
                score += 25
                criterios.append('PESO_EXACTO')
            elif diff == 1:
                score += 15
                criterios.append('PESO_±1')
            elif diff == 2:
                score += 8
                criterios.append('PESO_±2')
            elif diff > 3:
                score -= 10
        
        # Año nacimiento
        rcv_año = cand.get('año_nac')
        if op_año is not None and rcv_año is not None:
            diff = abs(rcv_año - op_año)
            diffs['año_nac'] = diff
            
            if diff == 0:
                score += 20
                criterios.append('AÑO_NAC_EXACTO')
            elif diff == 1:
                score += 10
                criterios.append('AÑO_NAC_±1')
            elif diff > 2:
                score -= 15
        
        if score > mejor_score:
            mejor_score = score
            mejor_match = cand
            mejor_criterios = criterios
            mejor_diffs = diffs
    
    return mejor_match, mejor_score, mejor_criterios, mejor_diffs, len(candidatos)


def buscar_match(op_row):
    """Busca el mejor match usando estrategia híbrida."""
    
    # PRIORIDAD 1: Buscar por ID directo
    match, score, criterios, diffs, n_cand = buscar_por_id(op_row)
    
    if match is not None and score >= 150:  # Umbral alto para ID directo
        return match, score, criterios, diffs, n_cand, 'ID_DIRECTO'
    
    # PRIORIDAD 2: Buscar por características
    match2, score2, criterios2, diffs2, n_cand2 = buscar_por_caracteristicas(op_row)
    
    # Elegir el mejor
    if match is not None and match2 is not None:
        if score >= score2:
            return match, score, criterios, diffs, n_cand, 'ID_DIRECTO'
        else:
            return match2, score2, criterios2, diffs2, n_cand2, 'CARACTERISTICAS'
    elif match is not None:
        return match, score, criterios, diffs, n_cand, 'ID_DIRECTO'
    elif match2 is not None:
        return match2, score2, criterios2, diffs2, n_cand2, 'CARACTERISTICAS'
    else:
        return None, 0, [], {}, 0, 'NINGUNO'


def get_confianza(score, criterios, metodo):
    """Determina nivel de confianza."""
    
    if metodo == 'ID_DIRECTO':
        if score >= 200:
            return 'EXACTA_ID'
        elif score >= 150:
            return 'MUY_ALTA_ID'
        else:
            return 'ALTA_ID'
    else:
        tiene_talla = any('TALLA' in c for c in criterios)
        tiene_peso = any('PESO' in c for c in criterios)
        
        if tiene_talla and tiene_peso:
            if score >= 100:
                return 'EXACTA'
            elif score >= 80:
                return 'MUY_ALTA'
            elif score >= 60:
                return 'ALTA'
            else:
                return 'MEDIA'
        else:
            if score >= 70:
                return 'ALTA_SIN_MEDIDAS'
            elif score >= 50:
                return 'MEDIA_SIN_MEDIDAS'
            else:
                return 'BAJA_SIN_MEDIDAS'


def process_chunk(chunk_data):
    """Procesa un chunk de sujetos."""
    chunk_id, op_chunk, op_cols = chunk_data
    global session_dict
    
    results = []
    
    for idx, op_row in op_chunk.iterrows():
        row_dict = op_row.to_dict()
        match, score, criterios, diffs, n_cand, metodo = buscar_match(row_dict)
        
        result = {}
        
        # Columnas de op
        for col in op_cols:
            result[f'op_{col}'] = op_row.get(col, '')
        
        umbral = 150 if metodo == 'ID_DIRECTO' else 40
        
        if match is not None and score >= umbral:
            confianza = get_confianza(score, criterios, metodo)
            
            result['MATCH'] = 'SI'
            result['match_score'] = score
            result['match_confianza'] = confianza
            result['match_metodo'] = metodo
            result['match_criterios'] = '|'.join(criterios)
            result['match_n_candidatos'] = n_cand
            result['diff_fecha'] = diffs.get('fecha', '')
            result['diff_edad'] = diffs.get('edad', '')
            result['diff_talla'] = diffs.get('talla', '')
            result['diff_peso'] = diffs.get('peso', '')
            result['diff_año_nac'] = diffs.get('año_nac', '')
            
            # Datos de RCV
            for col, val in match.items():
                if not col.endswith('_dt') and not col.endswith('_int') and col != 'key':
                    result[f'rcv_{col}'] = val
            
            # Datos de Session
            id_session = match.get('id_session_from_trials')
            if id_session and not pd.isna(id_session):
                id_session = int(float(id_session))
                if id_session in session_dict:
                    for col, val in session_dict[id_session].items():
                        result[f'ses_{col}'] = val
        else:
            result['MATCH'] = 'NO'
            result['match_score'] = score if score else 0
            result['match_confianza'] = 'SIN_MATCH'
            result['match_metodo'] = metodo
            result['match_criterios'] = ''
            result['match_n_candidatos'] = n_cand
        
        results.append(result)
    
    return chunk_id, results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--task-id', type=int, default=0)
    parser.add_argument('--total-tasks', type=int, default=1)
    parser.add_argument('--cpus', type=int, default=None)
    parser.add_argument('--combine-only', action='store_true')
    args = parser.parse_args()
    
    if args.combine_only:
        combine_results(args.total_tasks)
        return
    
    print(f"\n{'='*60}")
    print(f"MATCHING HÍBRIDO - TAREA {args.task_id + 1}/{args.total_tasks}")
    print(f"{'='*60}")
    print("PRIORIDAD 1: ID directo (p_08070803)")
    print("PRIORIDAD 2: Características (fecha+sexo+edad+talla+peso)")
    print(f"{'='*60}\n")
    
    Path(TEMP_DIR).mkdir(parents=True, exist_ok=True)
    
    # 1. Cargar datos
    print("=== CARGANDO DATOS ===")
    
    op = pd.read_csv(os.path.join(BASE_DIR, "op.csv"), dtype=str, low_memory=False)
    op_cols = list(op.columns)
    print(f"op.csv: {len(op):,}")
    
    rcv = pd.read_csv(os.path.join(BASE_DIR, "Rcv_trials.csv"), dtype=str, low_memory=False)
    print(f"Rcv_trials.csv: {len(rcv):,}")
    
    tri = pd.read_csv(os.path.join(BASE_DIR, "Trials.csv"), dtype=str, low_memory=False)
    print(f"Trials.csv: {len(tri):,}")
    
    ses = pd.read_csv(os.path.join(BASE_DIR, "Session.csv"), dtype=str, low_memory=False)
    print(f"Session.csv: {len(ses):,}")
    
    # 2. Relacionar Rcv con Trials
    print("\n=== RELACIONANDO TABLAS ===")
    
    rcv['key'] = rcv['trialdate'].astype(str) + '|' + rcv['trialtime'].astype(str)
    tri['key'] = tri['trialdate'].astype(str) + '|' + tri['trialtime'].astype(str)
    
    tri_subset = tri[['key', 'id_session']].drop_duplicates(subset=['key'], keep='first')
    tri_subset = tri_subset.rename(columns={'id_session': 'id_session_from_trials'})
    
    rcv = rcv.merge(tri_subset, on='key', how='left')
    print(f"  Rcv con id_session: {rcv['id_session_from_trials'].notna().sum():,}")
    
    # 3. Preparar op.csv
    print("\n=== PREPARANDO DATOS ===")
    
    op['fecha_dia'] = op['fecha'].apply(parse_fecha_solo_dia)
    op['fecha_dt'] = op['fecha_dia'].apply(lambda x: fecha_to_datetime(x) if x else None)
    op['fecha_dt'] = pd.to_datetime(op['fecha_dt'])
    
    op['talla_int'] = op['p_080610'].apply(lambda x: math.floor(float(x)) if pd.notna(x) and x != '' else None)
    op['peso_int'] = op['p_080612'].apply(lambda x: math.floor(float(x)) if pd.notna(x) and x != '' else None)
    op['edad_int'] = op['edad_hoy'].apply(lambda x: int(float(x)) if pd.notna(x) and x != '' else None)
    
    op['año_encuesta'] = op['fecha'].apply(extraer_año)
    op['año_nac'] = op.apply(lambda r: r['año_encuesta'] - r['edad_int'] if r['año_encuesta'] and r['edad_int'] else None, axis=1)
    
    n_con_id = (op['p_08070803'].notna() & (op['p_08070803'] != '')).sum()
    n_con_medidas = ((op['talla_int'].notna()) & (op['peso_int'].notna())).sum()
    print(f"  Sujetos con p_08070803: {n_con_id:,}")
    print(f"  Sujetos con talla+peso: {n_con_medidas:,}")
    
    # 4. Preparar Rcv - convertir campos numéricos
    print("\n=== INDEXANDO ESPIROMETRÍAS ===")
    
    rcv['trialdate_dt'] = rcv['trialdate'].apply(lambda x: fecha_to_datetime(str(x)) if pd.notna(x) else None)
    rcv['trialdate_dt'] = pd.to_datetime(rcv['trialdate_dt'])
    rcv['height_int'] = rcv['height'].apply(lambda x: int(float(x)) if pd.notna(x) and x != '' else None)
    rcv['weight_int'] = rcv['weight'].apply(lambda x: int(float(x)) if pd.notna(x) and x != '' else None)
    rcv['age_int'] = rcv['age'].apply(lambda x: int(float(x)) if pd.notna(x) and x != '' else None)
    rcv['año_nac'] = rcv['birthdate'].apply(extraer_año)
    
    # Índice por fecha (para búsqueda por características)
    rcv_by_date = {}
    for idx, row in rcv.iterrows():
        fecha = row.get('trialdate')
        if pd.isna(fecha):
            continue
        fecha = str(fecha).strip()
        if fecha not in rcv_by_date:
            rcv_by_date[fecha] = []
        rcv_by_date[fecha].append(row.to_dict())
    
    print(f"  Fechas únicas: {len(rcv_by_date):,}")
    
    # Índice por ID (para búsqueda directa)
    rcv_by_id = {}
    for idx, row in rcv.iterrows():
        rcv_id = row.get('id')
        if pd.isna(rcv_id):
            continue
        rcv_id = str(rcv_id).strip()
        if rcv_id not in rcv_by_id:
            rcv_by_id[rcv_id] = []
        rcv_by_id[rcv_id].append(row.to_dict())
    
    print(f"  IDs únicos: {len(rcv_by_id):,}")
    
    # 5. Session dict
    session_dict = {}
    for idx, row in ses.iterrows():
        id_ses = row.get('id_session')
        if pd.notna(id_ses):
            session_dict[int(float(id_ses))] = row.to_dict()
    
    print(f"  Sesiones: {len(session_dict):,}")
    
    # 6. Dividir trabajo
    total_rows = len(op)
    rows_per_task = total_rows // args.total_tasks
    start_idx = args.task_id * rows_per_task
    end_idx = total_rows if args.task_id == args.total_tasks - 1 else start_idx + rows_per_task
    
    op_subset = op.iloc[start_idx:end_idx].copy()
    
    print(f"\n=== DIVISIÓN ===")
    print(f"Esta tarea: {len(op_subset):,} sujetos")
    
    # 7. Procesar
    n_cpus = args.cpus or cpu_count()
    chunk_size = max(100, len(op_subset) // (n_cpus * 2))
    
    chunks = []
    for i in range(0, len(op_subset), chunk_size):
        chunk = op_subset.iloc[i:i+chunk_size]
        chunks.append((len(chunks), chunk, op_cols))
    
    print(f"\n=== PROCESANDO ===")
    print(f"CPUs: {n_cpus}")
    
    inicio = time.time()
    
    with Pool(n_cpus, initializer=init_worker, initargs=(rcv_by_date, rcv_by_id, session_dict)) as pool:
        results = []
        for i, (chunk_id, chunk_results) in enumerate(pool.imap_unordered(process_chunk, chunks)):
            results.append((chunk_id, chunk_results))
            if (i + 1) % 10 == 0:
                pct = (i+1) / len(chunks) * 100
                print(f"  {pct:.0f}%")
    
    print(f"  Tiempo: {time.time()-inicio:.1f}s")
    
    # Combinar
    results.sort(key=lambda x: x[0])
    all_results = []
    for _, chunk_results in results:
        all_results.extend(chunk_results)
    
    df_results = pd.DataFrame(all_results)
    
    # Guardar
    output_file = os.path.join(TEMP_DIR, f"match_task_{args.task_id:04d}.csv")
    df_results.to_csv(output_file, index=False)
    
    # Stats
    n_match = (df_results['MATCH'] == 'SI').sum()
    n_total = len(df_results)
    
    print(f"\n{'='*60}")
    print(f"ESTADÍSTICAS TAREA {args.task_id + 1}")
    print(f"{'='*60}")
    print(f"Total: {n_total:,}")
    print(f"Con match: {n_match:,} ({n_match/n_total*100:.1f}%)")
    
    if n_match > 0:
        print(f"\nPor método:")
        for metodo in ['ID_DIRECTO', 'CARACTERISTICAS']:
            n = (df_results['match_metodo'] == metodo).sum()
            if n > 0:
                print(f"  {metodo}: {n:,}")
        
        print(f"\nPor confianza:")
        for conf in df_results['match_confianza'].value_counts().index:
            if conf != 'SIN_MATCH':
                n = (df_results['match_confianza'] == conf).sum()
                print(f"  {conf}: {n:,}")
    
    # Si es última tarea, combinar
    if args.task_id == args.total_tasks - 1:
        time.sleep(15)
        combine_results(args.total_tasks)


def combine_results(total_tasks):
    """Combina resultados."""
    print(f"\n{'='*60}")
    print("COMBINANDO RESULTADOS")
    print(f"{'='*60}\n")
    
    files = sorted(Path(TEMP_DIR).glob("match_task_*.csv"))
    print(f"Archivos: {len(files)}")
    
    if not files:
        print("ERROR: No hay archivos")
        return
    
    dfs = []
    for f in files:
        dfs.append(pd.read_csv(f, dtype=str, low_memory=False))
    
    df = pd.concat(dfs, ignore_index=True)
    
    df['match_score_num'] = pd.to_numeric(df['match_score'], errors='coerce').fillna(0)
    df = df.sort_values(['MATCH', 'match_score_num'], ascending=[False, False])
    df = df.drop(columns=['match_score_num'])
    
    output = os.path.join(BASE_DIR, "matching_4fuentes_completo.csv")
    df.to_csv(output, index=False, encoding='utf-8-sig')
    
    n_match = (df['MATCH'] == 'SI').sum()
    n_total = len(df)
    
    print(f"\n{'='*60}")
    print("ESTADÍSTICAS FINALES")
    print(f"{'='*60}")
    print(f"Archivo: {output}")
    print(f"Columnas: {len(df.columns)}")
    print(f"Total: {n_total:,}")
    print(f"Con match: {n_match:,} ({n_match/n_total*100:.1f}%)")
    
    if n_match > 0:
        print(f"\nPor método:")
        for metodo in ['ID_DIRECTO', 'CARACTERISTICAS']:
            n = (df['match_metodo'] == metodo).sum()
            if n > 0:
                pct = n/n_match*100
                print(f"  {metodo}: {n:,} ({pct:.1f}%)")
        
        print(f"\nPor confianza:")
        for conf in df['match_confianza'].value_counts().index:
            if conf != 'SIN_MATCH':
                n = (df['match_confianza'] == conf).sum()
                pct = n/n_match*100
                print(f"  {conf}: {n:,} ({pct:.1f}%)")
        
        if 'ses_id_session' in df.columns:
            n_con_session = df[df['MATCH'] == 'SI']['ses_id_session'].notna().sum()
            print(f"\nCon datos de Session: {n_con_session:,}")
    
    # Limpiar
    for f in files:
        f.unlink()
    
    print("\n¡COMPLETADO!")


if __name__ == '__main__':
    main()
