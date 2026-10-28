"""Módulo 4 — Periodograma de Lomb del tacograma.

Consigna: ítems 7a (reposo vs. actividad física) y 7b (comparación entre
voluntarios).
Referencia: Lomb, N. R. (1976). Astrophysics and Space Science, 39, 447-462.

Qué recibe : t_rr y rr (salida del módulo 3, en segundos).
Qué devuelve: f (frecuencias en Hz) y p (potencia), de igual largo, p >= 0.
"""

import numpy as np

# Bandas de frecuencia para el análisis: {"NOMBRE": (f_min, f_max)} en Hz.
# PENDIENTE (M4): completar con bibliografía y citar la fuente en el informe.
BANDAS = {}


def periodograma_lomb(t, x, f=None):
    """Periodograma de Lomb de una serie con muestras no equiespaciadas.

    PENDIENTE (M4):
      - [ ] Elegir las frecuencias a evaluar cuando f no se pasa, y justificarlas.
      - [ ] Decidir qué tratamiento previo recibe el tacograma y documentarlo.
      - [ ] Documentar cómo se normaliza p.
    """
    raise NotImplementedError("M4: falta programar periodograma_lomb()")


def potencia_en_banda(f, p, banda):
    """Área bajo el periodograma entre banda[0] y banda[1] (en Hz)."""
    raise NotImplementedError("M4: falta programar potencia_en_banda()")
