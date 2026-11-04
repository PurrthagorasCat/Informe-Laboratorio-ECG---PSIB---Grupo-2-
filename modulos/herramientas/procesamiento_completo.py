"""Junta los módulos 1 a 4 en orden sobre un registro.

    M1 acondicionar()       señal acondicionada + máscara de muestras válidas
    M2 detectar_picos_r()   picos R sobre la señal acondicionada
    M3 intervalos_rr()      tacograma, sin los intervalos que tocan zonas inválidas
    M4 periodograma_lomb()  periodograma del tacograma

Con `hasta=` se corre sólo hasta una etapa, por ejemplo
procesar(reg, hasta="pan_tompkins").
"""

from dataclasses import dataclass

import numpy as np

from modulos.m0_carga_de_datos import carga_de_datos
from modulos.m0_carga_de_datos.carga_de_datos import Registro
from modulos.m1_acondicionamiento import acondicionamiento
from modulos.m2_pan_tompkins import pan_tompkins
from modulos.m3_tacograma import tacograma
from modulos.m4_periodograma_lomb import periodograma_lomb

ETAPAS = ("acondicionamiento", "pan_tompkins", "tacograma", "lomb")


@dataclass
class Resultado:
    """Todo lo que se obtiene de un registro. Se usa así: res.picos_r, res.rr..."""
    registro: Registro
    ecg_acondicionado: np.ndarray = None
    valida: np.ndarray = None          # True en las muestras aptas (M1, mascara_valida)
    picos_r: np.ndarray = None
    t_rr: np.ndarray = None
    rr: np.ndarray = None
    f: np.ndarray = None
    p: np.ndarray = None


def procesar(reg, hasta="lomb"):
    if hasta not in ETAPAS:
        raise ValueError(f"'hasta' tiene que ser una de estas: {ETAPAS}")
    ultima = ETAPAS.index(hasta)
    res = Resultado(registro=reg)

    res.ecg_acondicionado = acondicionamiento.acondicionar(reg.ecg, reg.fs)
    res.valida = acondicionamiento.mascara_valida(reg.ecg, reg.fs)
    if ultima < 1:
        return res
    res.picos_r = pan_tompkins.detectar_picos_r(res.ecg_acondicionado, reg.fs)
    if ultima < 2:
        return res
    res.t_rr, res.rr = tacograma.intervalos_rr(res.picos_r, reg.fs, valida=res.valida)
    if ultima < 3:
        return res
    res.f, res.p = periodograma_lomb.periodograma_lomb(res.t_rr, res.rr)
    return res


def procesar_todos(hasta="lomb"):
    """Procesa todos los registros de datos/lista_de_registros.csv."""
    return [procesar(reg, hasta=hasta) for reg in carga_de_datos.cargar_todos()]
