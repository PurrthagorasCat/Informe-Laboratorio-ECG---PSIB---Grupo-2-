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

# Formato del archivo que exporta la interfaz (ítem 1)
COLUMNA_SENAL = "Sample"   # encabezado de la única columna
UNIDADES = "V"             # tensión a la salida del módulo AD8232

# Frecuencia de muestreo (ítem 2). El archivo no la guarda: sale del firmware
# del Arduino (Timer1 en modo CTC, 16 MHz / (64 · (1 + OCR1A)) con OCR1A = 499)
# y es consistente con la duración de los registros: 28 800 muestras = 57,6 s
# (protocolo: 1 min) y 86 600 muestras = 173,2 s (protocolo: 3 min).
# Justificación completa en el informe, sección 2.2.3.
FS_REGISTRO = 500.0        # Hz


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

    Formato (ítem 1, informe sección 2.2.2): CSV de una sola columna con el
    encabezado "Sample" y una muestra por fila. Cada valor es la tensión a la
    salida del AD8232 en voltios (V = 5 n / 1023, redondeada a 10 mV).
    El archivo no guarda la base de tiempo: fs es FS_REGISTRO (ítem 2).
    """
    ruta = Path(ruta)
    tabla = pd.read_csv(ruta)
    columna = COLUMNA_SENAL if COLUMNA_SENAL in tabla.columns else tabla.columns[0]
    ecg = pd.to_numeric(tabla[columna], errors="coerce").to_numpy(dtype=float)

    # Filas vacías o no numéricas (por ejemplo, una última línea cortada si la
    # grabación se interrumpió): se avisa y se descartan.
    malas = ~np.isfinite(ecg)
    if malas.any():
        warnings.warn(f"{ruta.name}: se descartan {malas.sum()} filas no numéricas")
        ecg = ecg[~malas]

    return Registro(ecg=ecg, fs=FS_REGISTRO, archivo=ruta.name, unidades=UNIDADES)


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
