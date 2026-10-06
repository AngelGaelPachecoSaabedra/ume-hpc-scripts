#!/usr/bin/env python3
"""
MATCHING AUTOMÁTICO DE SUJETOS - VERSIÓN CORREGIDA
===================================================

Busca automáticamente la mejor coincidencia para cada sujeto de op.csv
en el archivo Rcv_trials.csv usando múltiples criterios.

USO:
    python matching_simple_v2.py

ARCHIVOS REQUERIDOS:
    - op.csv: Datos de los sujetos (42,994 registros)
    - Rcv_trials.csv: Datos de trials con info de sujetos (id, birthdate, weight, height, age)
    - Session.csv: Datos adicionales de sesión (se une por id_session)
"""

import pandas as pd
import numpy as np
from datetime import datetime
import sys
import warnings

warnings.filterwarnings('ignore')

print("=" * 60)
print("MATCHING AUTOMÁTICO DE SUJETOS CON SESIONES")
print("Versión 2.0 - Corregida para estructura real de datos")
print("=" * 60)

# ============================================================
# CONFIGURACIÓN
# ============================================================

# Archivos de entrada
ARCHIVO_SUJETOS = "op.csv"
ARCHIVO_RCV_TRIALS = "Rcv_trials.csv"  # Este tiene los datos de sujetos
ARCHIVO_SESSION = "Session.csv"  # Datos adicionales de sesión

# Archivos de salida
ARCHIVO_SALIDA = "resultados_matching.csv"
ARCHIVO_SIN_MATCH = "sin_coincidencia.csv"
ARCHIVO_LOG = "matching_log.txt"

# Tolerancias para considerar una coincidencia
TOLERANCIA_PESO = 2.0  # kg de diferencia permitida
TOLERANCIA_TALLA = 2.0  # cm de diferencia permitida
TOLERANCIA_EDAD = 2  # años de diferencia permitida
TOLERANCIA_FECHA = 7  # días de diferencia permitida

# Puntajes
SCORE_ID_EXACTO = 100  # Match exacto de ID
SCORE_FECHA_NAC = 80  # Match de fecha de nacimiento
SCORE_PESO = 25  # Coincidencia de peso
SCORE_TALLA = 25  # Coincidencia de talla
SCORE_EDAD = 20  # Coincidencia de edad
SCORE_FECHA_TRIAL = 35  # Coincidencia de fecha de intento
SCORE_SEXO = 10  # Coincidencia de sexo

# Umbral mínimo de score para aceptar una coincidencia
UMBRAL_MINIMO = 50


# ============================================================
# FUNCIONES AUXILIARES
# ============================================================

def extraer_digitos(valor):
    """Extrae solo los dígitos de un string."""
    if pd.isna(valor):
        return None
    resultado = ''.join(c for c in str(valor) if c.isdigit())
    return resultado if resultado else None


def parsear_fecha(fecha_str):
    """Intenta parsear una fecha con varios formatos."""
    if pd.isna(fecha_str) or str(fecha_str).strip() == '':
        return None

    fecha_str = str(fecha_str).strip()

    # Si tiene espacio, tomar solo la parte de fecha
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
    """Normaliza valores de sexo a M/F."""
    if pd.isna(valor):
        return None
    valor = str(valor).strip().upper()
    if valor in ['M', 'MASCULINO', 'MALE', '1', 'H', 'HOMBRE']:
        return 'M'
    elif valor in ['F', 'FEMENINO', 'FEMALE', '2', 'MUJER']:
        return 'F'
    return valor


def calcular_score(sujeto, trial):
    """
    Calcula el puntaje de coincidencia entre un sujeto y un trial.
    Retorna (score, detalles_list)
    """
    score = 0
    detalles = []

    # 1. Comparar ID (p_08070803 en op.csv vs id en Rcv_trials)
    id_sujeto = extraer_digitos(sujeto.get('p_08070803'))
    id_trial = extraer_digitos(trial.get('id'))

    if id_sujeto and id_trial and id_sujeto == id_trial:
        score += SCORE_ID_EXACTO
        detalles.append("ID exacto")

    # 2. Comparar fecha de nacimiento
    # En op.csv la fecha parece ser fecha de encuesta, no de nacimiento
    # Buscar si hay fecha de nacimiento calculable
    fecha_trial = parsear_fecha(trial.get('birthdate'))

    # 3. Comparar peso (p_080612 en op.csv vs weight en Rcv_trials)
    try:
        peso_suj = float(sujeto.get('p_080612') or 0)
        peso_trial = float(trial.get('weight') or 0)
        if peso_suj > 0 and peso_trial > 0:
            diff = abs(peso_suj - peso_trial)
            if diff <= TOLERANCIA_PESO:
                score += SCORE_PESO * (1 - diff / TOLERANCIA_PESO)
                detalles.append(f"Peso (Δ{diff:.1f}kg)")
    except:
        pass

    # 4. Comparar talla (p_080610 en op.csv vs height en Rcv_trials)
    try:
        talla_suj = float(sujeto.get('p_080610') or 0)
        talla_trial = float(trial.get('height') or 0)
        if talla_suj > 0 and talla_trial > 0:
            diff = abs(talla_suj - talla_trial)
            if diff <= TOLERANCIA_TALLA:
                score += SCORE_TALLA * (1 - diff / TOLERANCIA_TALLA)
                detalles.append(f"Talla (Δ{diff:.1f}cm)")
    except:
        pass

    # 5. Comparar edad (edad_hoy en op.csv vs age en Rcv_trials)
    try:
        edad_suj = float(sujeto.get('edad_hoy') or sujeto.get('edad') or 0)
        edad_trial = float(trial.get('age') or 0)
        if edad_suj > 0 and edad_trial > 0:
            diff = abs(edad_suj - edad_trial)
            if diff <= TOLERANCIA_EDAD:
                score += SCORE_EDAD * (1 - diff / TOLERANCIA_EDAD)
                detalles.append(f"Edad (Δ{diff:.0f}años)")
    except:
        pass

    # 6. Comparar fecha del trial (fecha en op.csv vs trialdate en Rcv_trials)
    try:
        fecha_suj = parsear_fecha(sujeto.get('fecha'))
        fecha_trial_date = parsear_fecha(trial.get('trialdate'))

        if fecha_suj and fecha_trial_date:
            diff_dias = abs((fecha_trial_date - fecha_suj).days)
            if diff_dias <= TOLERANCIA_FECHA:
                score += SCORE_FECHA_TRIAL * (1 - diff_dias / TOLERANCIA_FECHA)
                detalles.append(f"Fecha (Δ{diff_dias}días)")
    except:
        pass

    # 7. Comparar sexo
    sexo_suj = normalizar_sexo(sujeto.get('sexo'))
    sexo_trial = normalizar_sexo(trial.get('sex'))

    if sexo_suj and sexo_trial and sexo_suj == sexo_trial:
        score += SCORE_SEXO
        detalles.append("Sexo")

    return score, detalles


def determinar_confianza(score, detalles):
    """Determina el nivel de confianza del match."""
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


# ============================================================
# PROGRAMA PRINCIPAL
# ============================================================

def main():
    log_file = open(ARCHIVO_LOG, 'w', encoding='utf-8')

    def log(msg):
        print(msg)
        log_file.write(msg + '\n')

    # 1. Cargar datos
    log("\n[1/5] Cargando archivos...")

    try:
        df_sujetos = pd.read_csv(ARCHIVO_SUJETOS, dtype=str, low_memory=False)
        log(f"  ✓ {ARCHIVO_SUJETOS}: {len(df_sujetos)} sujetos")
    except FileNotFoundError:
        log(f"  ✗ ERROR: No se encontró {ARCHIVO_SUJETOS}")
        sys.exit(1)

    try:
        df_rcv = pd.read_csv(ARCHIVO_RCV_TRIALS, low_memory=False)
        log(f"  ✓ {ARCHIVO_RCV_TRIALS}: {len(df_rcv)} registros")
    except FileNotFoundError:
        log(f"  ✗ ERROR: No se encontró {ARCHIVO_RCV_TRIALS}")
        sys.exit(1)

    try:
        df_session = pd.read_csv(ARCHIVO_SESSION, low_memory=False)
        log(f"  ✓ {ARCHIVO_SESSION}: {len(df_session)} sesiones")
    except FileNotFoundError:
        log(f"  ⚠ ADVERTENCIA: No se encontró {ARCHIVO_SESSION}")
        df_session = pd.DataFrame()

    # 2. Preparar datos de Rcv_trials (eliminar duplicados, quedarnos con el mejor)
    log("\n[2/5] Preparando base de datos de trials...")

    # Eliminar duplicados por id, manteniendo el que tenga más datos
    df_rcv_clean = df_rcv.drop_duplicates(subset=['id', 'birthdate'], keep='first')
    log(f"  → {len(df_rcv_clean)} registros únicos (por id + birthdate)")

    # 3. Filtrar sujetos
    log("\n[3/5] Identificando sujetos a procesar...")

    # Verificar si existe columna 'Localizado' o 'verificado'
    col_localizado = None
    for col in ['Localizado', 'localizado', 'verificado', 'Verificado']:
        if col in df_sujetos.columns:
            col_localizado = col
            break

    if col_localizado:
        df_pendientes = df_sujetos[
            df_sujetos[col_localizado].fillna('').str.strip().str.lower() != 'si'
            ].copy()
        log(f"  → {len(df_pendientes)} sujetos pendientes (columna: {col_localizado})")
    else:
        df_pendientes = df_sujetos.copy()
        log(f"  → Procesando todos los {len(df_pendientes)} sujetos (no se encontró columna de estado)")

    if len(df_pendientes) == 0:
        log("  ¡Todos los sujetos ya están procesados!")
        sys.exit(0)

    # 4. Buscar coincidencias
    log("\n[4/5] Buscando coincidencias...")
    log(f"  Configuración:")
    log(f"    - Tolerancia peso: ±{TOLERANCIA_PESO} kg")
    log(f"    - Tolerancia talla: ±{TOLERANCIA_TALLA} cm")
    log(f"    - Tolerancia edad: ±{TOLERANCIA_EDAD} años")
    log(f"    - Tolerancia fecha: ±{TOLERANCIA_FECHA} días")
    log(f"    - Umbral mínimo: {UMBRAL_MINIMO} puntos")

    resultados = []
    sin_match = []

    total = len(df_pendientes)

    # Convertir df_rcv_clean a lista de diccionarios para mejor rendimiento
    trials_list = df_rcv_clean.to_dict('records')

    for i, (idx, row) in enumerate(df_pendientes.iterrows()):
        sujeto = row.to_dict()
        nombre = f"{sujeto.get('nombre', '')} {sujeto.get('paterno', '')} {sujeto.get('materno', '')}".strip()

        # Mostrar progreso
        if (i + 1) % 500 == 0 or i == 0 or (i + 1) == total:
            log(f"  Procesando: {i + 1:,}/{total:,} ({(i + 1) / total * 100:.1f}%)")

        # Buscar mejor coincidencia
        mejor_score = 0
        mejor_trial = None
        mejor_detalles = []

        for trial in trials_list:
            score, detalles = calcular_score(sujeto, trial)

            if score > mejor_score:
                mejor_score = score
                mejor_trial = trial
                mejor_detalles = detalles

        if mejor_score >= UMBRAL_MINIMO and mejor_trial is not None:
            confianza = determinar_confianza(mejor_score, mejor_detalles)

            # Buscar datos de sesión correspondiente
            session_data = {}
            if len(df_session) > 0 and mejor_trial.get('id_session'):
                session_match = df_session[df_session['id_session'] == mejor_trial['id_session']]
                if len(session_match) > 0:
                    session_data = session_match.iloc[0].to_dict()

            resultado = {
                # === IDENTIFICACIÓN ===
                "registro": sujeto.get('registro', ''),
                "nombre_completo": nombre,
                "nombre": sujeto.get('nombre', ''),
                "paterno": sujeto.get('paterno', ''),
                "materno": sujeto.get('materno', ''),

                # === DATOS DEL SUJETO (op.csv) ===
                "id_sujeto": extraer_digitos(sujeto.get('p_08070803')),
                "fecha_encuesta": sujeto.get('fecha', ''),
                "peso_sujeto": sujeto.get('p_080612', ''),
                "talla_sujeto": sujeto.get('p_080610', ''),
                "edad_sujeto": sujeto.get('edad_hoy', ''),
                "sexo_sujeto": sujeto.get('sexo', ''),

                # === DATOS DEL TRIAL ENCONTRADO (Rcv_trials.csv) ===
                "id_trial": mejor_trial.get('id', ''),
                "id_session": mejor_trial.get('id_session', ''),
                "birthdate_trial": mejor_trial.get('birthdate', ''),
                "trialdate": mejor_trial.get('trialdate', ''),
                "peso_trial": mejor_trial.get('weight', ''),
                "talla_trial": mejor_trial.get('height', ''),
                "edad_trial": mejor_trial.get('age', ''),
                "sexo_trial": mejor_trial.get('sex', ''),

                # === DATOS DE ESPIROMETRÍA ===
                "fvcbest": session_data.get('fvcbest', mejor_trial.get('fvc_corr', '')),
                "fev1best": session_data.get('fev1best', mejor_trial.get('bestfev1', '')),
                "pefbest": session_data.get('pefbest', ''),
                "bestfvcfev": session_data.get('bestfvcfev', ''),
                "reproducib": session_data.get('reproducib', ''),
                "fef2575_pc": session_data.get('fef2575_pc', ''),

                # === METADATOS DEL MATCHING ===
                "match_score": round(mejor_score, 1),
                "confianza": confianza,
                "criterios_match": " | ".join(mejor_detalles),
                "num_criterios": len(mejor_detalles),

                # === DATOS ADICIONALES ===
                "municipio": sujeto.get('municipio', ''),
                "cadena": sujeto.get('cadena', ''),
                "p_0811": sujeto.get('p_0811', ''),
                "subjectcod": mejor_trial.get('subjectcod', ''),
            }

            resultados.append(resultado)
        else:
            sin_match.append({
                "registro": sujeto.get('registro', ''),
                "nombre_completo": nombre,
                "id_sujeto": extraer_digitos(sujeto.get('p_08070803')),
                "fecha": sujeto.get('fecha', ''),
                "peso": sujeto.get('p_080612', ''),
                "talla": sujeto.get('p_080610', ''),
                "edad": sujeto.get('edad_hoy', ''),
                "sexo": sujeto.get('sexo', ''),
                "mejor_score": round(mejor_score, 1) if mejor_score > 0 else 0,
                "municipio": sujeto.get('municipio', ''),
            })

    # 5. Guardar resultados
    log("\n[5/5] Guardando resultados...")

    df_resultados = pd.DataFrame(resultados)
    df_sin_match = pd.DataFrame(sin_match)

    if len(df_resultados) > 0:
        df_resultados = df_resultados.sort_values('match_score', ascending=False)
        df_resultados.to_csv(ARCHIVO_SALIDA, index=False, encoding='utf-8-sig')
        log(f"  ✓ {ARCHIVO_SALIDA}: {len(df_resultados):,} coincidencias")

    if len(df_sin_match) > 0:
        df_sin_match = df_sin_match.sort_values('mejor_score', ascending=False)
        df_sin_match.to_csv(ARCHIVO_SIN_MATCH, index=False, encoding='utf-8-sig')
        log(f"  ✓ {ARCHIVO_SIN_MATCH}: {len(df_sin_match):,} sin coincidencia")

    # Resumen final
    log("\n" + "=" * 60)
    log("RESUMEN FINAL")
    log("=" * 60)
    log(f"Total sujetos` procesados: {total:,}")
    log(f"Coincidencias encontradas: {len(resultados):,}")
    log(f"Sin coincidencia: {len(sin_match):,}")
    log(f"Tasa de match: {len(resultados) / total * 100:.1f}%")

    if len(df_resultados) > 0:
        log("\nDistribución por nivel de confianza:")
        for nivel in ['MUY ALTA', 'ALTA', 'MEDIA-ALTA', 'MEDIA', 'BAJA']:
            count = len(df_resultados[df_resultados['confianza'] == nivel])
            if count > 0:
                pct = count / len(df_resultados) * 100
                log(f"  {nivel:12}: {count:6,} ({pct:5.1f}%)")

        log("\nDistribución por número de criterios coincidentes:")
        for n_crit in sorted(df_resultados['num_criterios'].unique()):
            count = len(df_resultados[df_resultados['num_criterios'] == n_crit])
            log(f"  {n_crit} criterios: {count:,}")

    log("\n" + "=" * 60)
    log("¡PROCESO COMPLETADO!")
    log("=" * 60)

    log_file.close()
    print(f"\nLog guardado en: {ARCHIVO_LOG}")


if __name__ == "__main__":
    main()
