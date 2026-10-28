"""Módulo 0 — Carga de datos.

Consigna: ítems 1 (formato del registro) y 2 (frecuencia de muestreo).

Qué tiene que devolver este módulo: un `Registro` con
    ecg : la señal cruda (array de numpy de una dimensión)
    fs  : la frecuencia de muestreo en Hz
El resto de los módulos sólo usa esos dos datos.
"""

import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

CARPETA_DATOS = Path(__file__).resolve().parents[2] / "datos"   # carpeta datos/ del repo
CARPETA_REGISTROS_CRUDOS = CARPETA_DATOS / "registros_crudos"
LISTA_DE_REGISTROS = CARPETA_DATOS / "lista_de_registros.csv"

VOLUNTARIOS = ("V1", "V2")
CONDICIONES = ("reposo", "post_ejercicio")


@dataclass
class Registro:
    """Un registro de ECG. Se usa así: reg.ecg, reg.fs, reg.voluntario..."""
    ecg: np.ndarray                 # señal cruda
    fs: float                       # frecuencia de muestreo [Hz]
    voluntario: str = ""            # "V1" o "V2"
    condicion: str = ""             # "reposo" o "post_ejercicio"
    archivo: str = ""               # nombre del archivo en datos/registros_crudos
    unidades: str = "sin definir"   # completar: cuentas del conversor, V, mV...

    @property
    def nombre(self):
        return f"{self.voluntario}_{self.condicion}"

    @property
    def duracion(self):
        """Duración del registro en segundos."""
        return len(self.ecg) / self.fs


def cargar_registro(ruta):
    """Lee un archivo tal como lo exporta la interfaz de adquisición.

    PENDIENTE (M0):
      - [ ] Identificar el formato del archivo exportado (ítem 1).
      - [ ] Identificar columnas, separador, encabezado y unidades.
      - [ ] Determinar la frecuencia de muestreo y dejar la justificación
            en el cuaderno carga_de_datos.ipynb (ítem 2).
      - [ ] Devolver un Registro con `ecg` y `fs`.
    """
    raise NotImplementedError("M0: falta programar cargar_registro()")


def leer_lista_de_registros():
    """Lee datos/lista_de_registros.csv como tabla."""
    return pd.read_csv(LISTA_DE_REGISTROS, dtype=str, keep_default_na=False)


def cargar_todos():
    """Carga los registros que figuran en la lista. Saltea las filas sin archivo."""
    registros = []
    for _, fila in leer_lista_de_registros().iterrows():
        if not fila["archivo"]:
            warnings.warn(f"{fila['voluntario']}_{fila['condicion']}: falta el nombre del archivo en la lista")
            continue
        reg = cargar_registro(CARPETA_REGISTROS_CRUDOS / fila["archivo"])
        reg.voluntario = fila["voluntario"]
        reg.condicion = fila["condicion"]
        reg.archivo = fila["archivo"]
        registros.append(reg)
    return registros
