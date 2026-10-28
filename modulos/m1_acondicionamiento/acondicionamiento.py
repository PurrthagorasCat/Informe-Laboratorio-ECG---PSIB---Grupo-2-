"""Módulo 1 — Acondicionamiento de la señal.

Consigna: ítems 4 (¿hace falta acondicionar? ¿qué ruidos asociados al
registro aparecen?) y 5 (acondicionar la señal).

Qué recibe : la señal cruda y fs.
Qué devuelve: la señal acondicionada, con la MISMA cantidad de muestras y SIN
              CORRIMIENTO en el tiempo (un pico que estaba en la muestra n
              sigue en la muestra n). Los módulos 2 y 3 dependen de esto.
"""

import numpy as np


def acondicionar(ecg, fs):
    """Aplica las etapas de acondicionamiento justificadas en el ítem 4.

    PENDIENTE (M1):
      - [ ] Identificar los ruidos presentes en los registros (cuaderno acondicionamiento.ipynb).
      - [ ] Definir y justificar cada etapa.
      - [ ] Programar cada etapa como una función aparte y llamarlas acá en orden.
      - [ ] Comprobar: misma cantidad de muestras y sin corrimiento.
    """
    raise NotImplementedError("M1: falta programar acondicionar()")
