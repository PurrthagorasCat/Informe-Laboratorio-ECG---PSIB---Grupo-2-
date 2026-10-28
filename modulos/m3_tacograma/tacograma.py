"""Módulo 3 — Tacograma y frecuencia cardíaca.

Consigna: ítems 6a (cómo varía la FC en cada registro) y 6b (comparación
entre voluntarios).

Qué recibe : las posiciones de los picos R (salida del módulo 2) y fs.
Qué devuelve: t_rr y rr, en SEGUNDOS, con un elemento menos que los picos.
              rr[i]   = tiempo entre el latido i y el latido i+1
              t_rr[i] = instante del latido i+1 (el que cierra el intervalo)
Todo en segundos; pasar a milisegundos sólo para los gráficos.
"""

import numpy as np


def intervalos_rr(picos_r, fs):
    """Arma el tacograma (t_rr, rr) a partir de los picos R."""
    raise NotImplementedError("M3: falta programar intervalos_rr()")


def frecuencia_cardiaca(rr):
    """Frecuencia cardíaca instantánea en latidos por minuto (rr en segundos)."""
    raise NotImplementedError("M3: falta programar frecuencia_cardiaca()")


def resumen(rr):
    """Valores descriptivos del tacograma para la tabla del informe.

    PENDIENTE (M3): definir qué valores reportar para responder 6a y 6b.
    """
    raise NotImplementedError("M3: falta programar resumen()")
