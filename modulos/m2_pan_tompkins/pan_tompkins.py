"""Módulo 2 — Detección de latidos con el método de Pan-Tompkins.

Qué recibe : la señal acondicionada (salida del módulo 1) y fs.
Qué devuelve: las posiciones (números de muestra, enteros) de los picos R, ordenadas, sin repetidos y
            referidas a la señal que entró (si alguna etapa corre la señal en el tiempo, corregirlo).
"""

import numpy as np


def detectar_picos_r(ecg, fs):
    """Función que usan los demás módulos: devuelve sólo los picos R."""
    return pan_tompkins(ecg, fs)["picos_r"]


def pan_tompkins(ecg, fs):
    """Corre todo el algoritmo y devuelve también las señales de cada etapa
    (sirven para la figura del informe que muestra el paso a paso)."""
    filtrada = filtro_pasabanda(ecg, fs)
    derivada = derivar(filtrada, fs)
    al_cuadrado = elevar_al_cuadrado(derivada)
    integrada = integrar_con_ventana_movil(al_cuadrado, fs)
    picos_r = decidir_picos(ecg, integrada, fs)
    return {
        "filtrada": filtrada,
        "derivada": derivada,
        "al_cuadrado": al_cuadrado,
        "integrada": integrada,
        "picos_r": picos_r,
    }


# ---------------------------- Etapas ----------------------------

def filtro_pasabanda(x, fs):


def derivar(x, fs):


def elevar_al_cuadrado(x):


def integrar_con_ventana_movil(x, fs):


def decidir_picos(ecg, integrada, fs):
