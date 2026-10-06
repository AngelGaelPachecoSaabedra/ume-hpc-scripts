#!/usr/bin/env python3
"""
MATCHING AUTOMÁTICO PARALELIZADO PARA HPC/SLURM
================================================

Versión optimizada para clusters HPC usando multiprocessing.
Divide el trabajo entre múltiples CPUs para procesar más rápido.

USO LOCAL:
    python matching_paralelo.py --cpus 8

USO CON SLURM:
    sbatch job_matching.sh
"""

import pandas as pd
import numpy as np
from datetime import datetime
import sys
import warnings
import argparse
from multiprocessing import Pool, cpu_count
from functools import partial
import time

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
SCORE_FECHA_NAC = 80
SCORE_PESO = 25
SCORE_TALLA = 25
SCORE_EDAD = 20
SCORE_FECHA_TRIAL = 35
SCORE_SEXO = 10

UMBRAL_MINIMO = 50

# ============================================================
# FUNCIONES AUXILIARES (deben estar a nivel módulo para pickle)
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
    
    # 1. Comparar ID
    id_sujeto = extraer_digitos(sujeto.get('p_08070803'))
    id_trial = extraer_digitos(trial.get('id'))
    
    if id_sujeto and id_trial and id_sujeto == id_trial:
        score += SCORE_ID_EXACTO
        detalles.append("ID exacto")
    
    # 2. Comparar peso
    try:
        peso_suj = float(sujeto.get('p_080612') or 0)
        peso_trial = float(trial.get('weight') or 0)
        if peso_suj > 0 and peso_trial > 0:
            diff = abs(peso_suj - peso_trial)
            if diff <= TOLERANCIA_PESO:
                score += SCORE_PESO * (1 - diff/TOLERANCIA_PESO)
                detalles.append(f"Peso")
    except:
        pass
    
    # 3. Comparar talla
    try:
        talla_suj = float(sujeto.get('p_080610') or 0)
        talla_trial = float(trial.get('height') or 0)
        if talla_suj > 0 and talla_trial > 0:
            diff = abs(talla_suj - talla_trial)
            if diff <= TOLERANCIA_TALLA:
                score += SCORE_TALLA * (1 - diff/TOLERANCIA_TALLA)
                detalles.append(f"Talla")
    except:
        pass
    
    # 4. Comparar edad
    try:
        edad_suj = float(sujeto.get('edad_hoy') or sujeto.get('edad') or 0)
        edad_trial = float(trial.get('age') or 0)
        if edad_suj > 0 and edad_trial > 0:
            diff = abs(edad_suj - edad_trial)
            if diff <= TOLERANCIA_EDAD:
                score += SCORE_EDAD * (1 - diff/TOLERANCIA_EDAD)
                detalles.append(f"Edad")
    except:
        pass
    
    # 5. Comparar fecha
    try:
        fecha_suj = parsear_fecha(sujeto.get('fecha'))
        fecha_trial_date = parsear_fecha(trial.get('trialdate'))
        
        if fecha_suj and fecha_trial_date:
            diff_dias = abs((fecha_trial_date - fecha_suj).days)
            if diff_dias <= TOLERANCIA_FECHA:
                score += SCORE_FECHA_TRIAL * (1 - diff_dias/TOLERANCIA_FECHA)
                detalles.append(f"Fecha")
    except:
        pass
    
    # 6. Comparar sexo
    sexo_suj = normalizar_sexo(sujeto.get('sexo'))
    sexo_trial = normalizar_sexo(trial.get('sex'))
    
    if sexo_suj and sexo_trial and sexo_suj == sexo_trial:
        score += SCORE_SEXO
        detalles.append("Sexo")
    
    return score, detalles

def determinar_confianza(score, detalles):
    tiene_id = "ID exacto" in detalles
    criterios_adicionales = len([d for d in detalles if d != "ID exacto"])
    
    if tiene_id and criterios_adicionales >= 3:
        return "MUY ALTA"
    elif tiene_id and criterios_adicionales >= 1:
        return "ALTA"
    elif tiene_id:
        return "ALTA"
    elif score >= 80:
        return "MEDIA-ALTA"
    elif score >= 60:
        return "MEDIA"
    else:
        return "BAJA"

def procesar_sujeto(sujeto_dict, trials_list, session_dict):
    """
    Procesa un solo sujeto - función para paralelización.
    """
    nombre = f"{sujeto_dict.get('nombre', '')} {sujeto_dict.get('paterno', '')} {sujeto_dict.get('materno', '')}".strip()
    
    mejor_score = 0
    mejor_trial = None
    mejor_detalles = []
    
    for trial in trials_list:
        score, detalles = calcular_score(sujeto_dict, trial)
        if score > mejor_score:
            mejor_score = score
            mejor_trial = trial
            mejor_detalles = detalles
    
    if mejor_score >= UMBRAL_MINIMO and mejor_trial is not None:
        confianza = determinar_confianza(mejor_score, mejor_detalles)
        
        # Buscar datos de sesión
        session_data = session_dict.get(mejor_trial.get('id_session'), {})
        
        return {
            'tipo': 'match',
            'data': {
                "registro": sujeto_dict.get('registro', ''),
                "nombre_completo": nombre,
                "nombre": sujeto_dict.get('nombre', ''),
                "paterno": sujeto_dict.get('paterno', ''),
                "materno": sujeto_dict.get('materno', ''),
                "id_sujeto": extraer_digitos(sujeto_dict.get('p_08070803')),
                "fecha_encuesta": sujeto_dict.get('fecha', ''),
                "peso_sujeto": sujeto_dict.get('p_080612', ''),
                "talla_sujeto": sujeto_dict.get('p_080610', ''),
                "edad_sujeto": sujeto_dict.get('edad_hoy', ''),
                "sexo_sujeto": sujeto_dict.get('sexo', ''),
                "id_trial": mejor_trial.get('id', ''),
                "id_session": mejor_trial.get('id_session', ''),
                "birthdate_trial": mejor_trial.get('birthdate', ''),
                "trialdate": mejor_trial.get('trialdate', ''),
                "peso_trial": mejor_trial.get('weight', ''),
                "talla_trial": mejor_trial.get('height', ''),
                "edad_trial": mejor_trial.get('age', ''),
                "sexo_trial": mejor_trial.get('sex', ''),
                "fvcbest": session_data.get('fvcbest', mejor_trial.get('fvc_corr', '')),
                "fev1best": session_data.get('fev1best', mejor_trial.get('bestfev1', '')),
                "pefbest": session_data.get('pefbest', ''),
                "bestfvcfev": session_data.get('bestfvcfev', ''),
                "reproducib": session_data.get('reproducib', ''),
                "match_score": round(mejor_score, 1),
                "confianza": confianza,
                "criterios_match": " | ".join(mejor_detalles),
                "num_criterios": len(mejor_detalles),
                "municipio": sujeto_dict.get('municipio', ''),
            }
        }
    else:
        return {
            'tipo': 'sin_match',
            'data': {
                "registro": sujeto_dict.get('registro', ''),
                "nombre_completo": nombre,
                "id_sujeto": extraer_digitos(sujeto_dict.get('p_08070803')),
                "fecha": sujeto_dict.get('fecha', ''),
                "peso": sujeto_dict.get('p_080612', ''),
                "talla": sujeto_dict.get('p_080610', ''),
                "edad": sujeto_dict.get('edad_hoy', ''),
                "mejor_score": round(mejor_score, 1),
            }
        }

def procesar_chunk(args):
    """Procesa un chunk de sujetos."""
    chunk_sujetos, trials_list, session_dict, chunk_id = args
    resultados = []
    sin_match = []
    
    for i, sujeto in enumerate(chunk_sujetos):
        resultado = procesar_sujeto(sujeto, trials_list, session_dict)
        if resultado['tipo'] == 'match':
            resultados.append(resultado['data'])
        else:
            sin_match.append(resultado['data'])
        
        # Progreso cada 1000 registros
        if (i + 1) % 1000 == 0:
            print(f"  [Chunk {chunk_id}] Procesados: {i+1}/{len(chunk_sujetos)}")
    
    return resultados, sin_match

# ============================================================
# PROGRAMA PRINCIPAL
# ============================================================

def main():
    parser = argparse.ArgumentParser(description='Matching paralelo de sujetos')
    parser.add_argument('--cpus', type=int, default=None, 
                        help='Número de CPUs a usar (default: todos disponibles)')
    parser.add_argument('--chunk-size', type=int, default=5000,
                        help='Tamaño de cada chunk de trabajo')
    args = parser.parse_args()
    
    # Determinar CPUs
    n_cpus = args.cpus or cpu_count()
    print("="*60)
    print("MATCHING AUTOMÁTICO PARALELIZADO")
    print("="*60)
    print(f"CPUs disponibles: {cpu_count()}")
    print(f"CPUs a usar: {n_cpus}")
    print(f"Tamaño de chunk: {args.chunk_size}")
    
    # 1. Cargar datos
    print("\n[1/4] Cargando archivos...")
    inicio = time.time()
    
    df_sujetos = pd.read_csv(ARCHIVO_SUJETOS, dtype=str, low_memory=False)
    print(f"  ✓ {ARCHIVO_SUJETOS}: {len(df_sujetos):,} sujetos")
    
    df_rcv = pd.read_csv(ARCHIVO_RCV_TRIALS, low_memory=False)
    df_rcv_clean = df_rcv.drop_duplicates(subset=['id', 'birthdate'], keep='first')
    print(f"  ✓ {ARCHIVO_RCV_TRIALS}: {len(df_rcv_clean):,} registros únicos")
    
    try:
        df_session = pd.read_csv(ARCHIVO_SESSION, low_memory=False)
        session_dict = df_session.set_index('id_session').to_dict('index')
        print(f"  ✓ {ARCHIVO_SESSION}: {len(df_session):,} sesiones")
    except:
        session_dict = {}
        print(f"  ⚠ {ARCHIVO_SESSION} no encontrado")
    
    # 2. Filtrar sujetos pendientes
    print("\n[2/4] Preparando datos...")
    
    col_localizado = None
    for col in ['Localizado', 'localizado', 'verificado']:
        if col in df_sujetos.columns:
            col_localizado = col
            break
    
    if col_localizado:
        df_pendientes = df_sujetos[
            df_sujetos[col_localizado].fillna('').str.strip().str.lower() != 'si'
        ]
    else:
        df_pendientes = df_sujetos
    
    print(f"  → {len(df_pendientes):,} sujetos a procesar")
    
    # Convertir a listas para paralelización
    sujetos_list = df_pendientes.to_dict('records')
    trials_list = df_rcv_clean.to_dict('records')
    
    # 3. Dividir en chunks y procesar en paralelo
    print(f"\n[3/4] Procesando en paralelo con {n_cpus} CPUs...")
    
    # Dividir sujetos en chunks
    chunk_size = args.chunk_size
    chunks = [sujetos_list[i:i+chunk_size] for i in range(0, len(sujetos_list), chunk_size)]
    print(f"  → {len(chunks)} chunks de ~{chunk_size} sujetos cada uno")
    
    # Preparar argumentos para cada chunk
    chunk_args = [(chunk, trials_list, session_dict, i) for i, chunk in enumerate(chunks)]
    
    # Procesar en paralelo
    inicio_proceso = time.time()
    
    with Pool(processes=n_cpus) as pool:
        resultados_chunks = pool.map(procesar_chunk, chunk_args)
    
    tiempo_proceso = time.time() - inicio_proceso
    
    # Combinar resultados
    todos_resultados = []
    todos_sin_match = []
    
    for resultados, sin_match in resultados_chunks:
        todos_resultados.extend(resultados)
        todos_sin_match.extend(sin_match)
    
    # 4. Guardar resultados
    print(f"\n[4/4] Guardando resultados...")
    
    df_resultados = pd.DataFrame(todos_resultados)
    df_sin_match = pd.DataFrame(todos_sin_match)
    
    if len(df_resultados) > 0:
        df_resultados = df_resultados.sort_values('match_score', ascending=False)
        df_resultados.to_csv(ARCHIVO_SALIDA, index=False, encoding='utf-8-sig')
        print(f"  ✓ {ARCHIVO_SALIDA}: {len(df_resultados):,} coincidencias")
    
    if len(df_sin_match) > 0:
        df_sin_match.to_csv(ARCHIVO_SIN_MATCH, index=False, encoding='utf-8-sig')
        print(f"  ✓ {ARCHIVO_SIN_MATCH}: {len(df_sin_match):,} sin coincidencia")
    
    # Resumen
    tiempo_total = time.time() - inicio
    
    print("\n" + "="*60)
    print("RESUMEN")
    print("="*60)
    print(f"Total procesados: {len(sujetos_list):,}")
    print(f"Coincidencias: {len(todos_resultados):,} ({len(todos_resultados)/len(sujetos_list)*100:.1f}%)")
    print(f"Sin match: {len(todos_sin_match):,}")
    print(f"\nTiempo de proceso: {tiempo_proceso:.1f} segundos")
    print(f"Tiempo total: {tiempo_total:.1f} segundos")
    print(f"Velocidad: {len(sujetos_list)/tiempo_proceso:.0f} sujetos/segundo")
    
    if len(df_resultados) > 0:
        print("\nDistribución por confianza:")
        for nivel in ['MUY ALTA', 'ALTA', 'MEDIA-ALTA', 'MEDIA', 'BAJA']:
            count = len(df_resultados[df_resultados['confianza'] == nivel])
            if count > 0:
                print(f"  {nivel}: {count:,}")
    
    print("\n¡COMPLETADO!")
    print("="*60)

if __name__ == "__main__":
    main()
