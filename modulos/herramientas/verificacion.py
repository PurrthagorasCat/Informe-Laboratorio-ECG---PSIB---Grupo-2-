"""Comprueba que cada módulo devuelve lo que los demás módulos esperan.

Lo usa main.py al arrancar. Resultado por módulo:
    LISTO   cumple con lo acordado
    FALTA   todavía tiene funciones sin programar (o faltan los registros)
    ERROR   está programado pero no cumple; se muestra qué falló
"""

import numpy as np

from modulos.herramientas.ecg_de_prueba import generar_ecg, generar_latidos
from modulos.m0_carga_de_datos import carga_de_datos
from modulos.m1_acondicionamiento import acondicionamiento
from modulos.m2_pan_tompkins import pan_tompkins
from modulos.m3_tacograma import tacograma
from modulos.m4_periodograma_lomb import periodograma_lomb


class FaltanRegistros(Exception):
    pass


def verificar_m0():
    lista = carga_de_datos.leer_lista_de_registros()
    assert {"archivo", "voluntario", "condicion"} <= set(lista.columns), \
        "a lista_de_registros.csv le faltan columnas"
    if not any(lista["archivo"]):
        raise FaltanRegistros("todavía no hay registros en datos/lista_de_registros.csv")
    for reg in carga_de_datos.cargar_todos():
        assert reg.ecg.ndim == 1, f"{reg.nombre}: ecg tiene que ser un array de una dimensión"
        assert np.issubdtype(reg.ecg.dtype, np.floating), f"{reg.nombre}: ecg tiene que ser de tipo float"
        assert np.all(np.isfinite(reg.ecg)), f"{reg.nombre}: ecg tiene valores NaN o infinitos"
        assert reg.fs > 0, f"{reg.nombre}: fs tiene que ser positiva"
        assert reg.duracion > 10, f"{reg.nombre}: dura {reg.duracion:.1f} s, sospechosamente corto"


def verificar_m1():
    for fs in (250.0, 500.0):
        _, ecg, picos = generar_ecg(fs=fs, duracion=20, semilla=0)
        y = acondicionamiento.acondicionar(ecg, fs)
        assert y.shape == ecg.shape, "la salida no tiene la misma cantidad de muestras que la entrada"
        assert np.all(np.isfinite(y)), "la salida tiene valores NaN o infinitos"
        ventana, tolerancia = int(0.05 * fs), int(0.01 * fs)
        for r in picos:
            i0, i1 = max(r - ventana, 0), min(r + ventana, len(y))
            corrimiento = i0 + np.argmax(y[i0:i1]) - r
            assert abs(corrimiento) <= tolerancia, \
                f"los picos R se corren {1000 * corrimiento / fs:.0f} ms (máximo 10 ms)"


def _cantidad_emparejados(a, b, tolerancia):
    """Cuántos elementos de b tienen algún elemento de a a <= tolerancia muestras."""
    return sum(np.any(np.abs(a - x) <= tolerancia) for x in b)


def verificar_m2():
    fs = 500.0
    _, ecg, _ = generar_ecg(fs=fs, duracion=10, semilla=0)
    etapas = pan_tompkins.pan_tompkins(ecg, fs)
    faltan = {"filtrada", "derivada", "al_cuadrado", "integrada", "picos_r"} - set(etapas)
    assert not faltan, f"pan_tompkins() no devuelve: {faltan}"

    for fc in (60, 120):
        _, ecg, reales = generar_ecg(fs=fs, duracion=30, fc_media=fc, ruido_std=0.02, semilla=1)
        det = pan_tompkins.detectar_picos_r(ecg, fs)
        assert np.issubdtype(det.dtype, np.integer), "los picos tienen que ser números enteros (muestras)"
        assert np.all(np.diff(det) > 0), "los picos tienen que estar ordenados y sin repetir"
        assert det.min() >= 0 and det.max() < len(ecg), "hay picos fuera de la señal"
        tol = int(0.05 * fs)
        sensibilidad = _cantidad_emparejados(det, reales, tol) / len(reales)
        vpp = _cantidad_emparejados(reales, det, tol) / max(len(det), 1)
        assert sensibilidad >= 0.95, f"a {fc} lpm detecta sólo el {100 * sensibilidad:.0f} % de los latidos"
        assert vpp >= 0.95, f"a {fc} lpm el {100 * (1 - vpp):.0f} % de las detecciones son falsas"


def verificar_m3():
    fs = 500.0
    picos = np.arange(10) * 400 + 100  # un latido cada 0,8 s
    t_rr, rr = tacograma.intervalos_rr(picos, fs)
    assert len(rr) == len(t_rr) == len(picos) - 1, "t_rr y rr tienen que tener un elemento menos que los picos"
    assert np.allclose(rr, 0.8), "rr tiene que estar en segundos"
    assert np.allclose(t_rr, picos[1:] / fs), "t_rr[i] tiene que ser el instante del latido i+1, en segundos"
    assert np.allclose(tacograma.frecuencia_cardiaca(np.array([1.0, 0.5])), [60, 120]), \
        "la frecuencia cardíaca tiene que estar en latidos por minuto"


def verificar_m4():
    t_r = generar_latidos(duracion=180, fc_media=70, amp_mod=0.05, f_mod=0.2)
    f, p = periodograma_lomb.periodograma_lomb(t_r[1:], np.diff(t_r))
    assert f.shape == p.shape, "f y p tienen que tener el mismo largo"
    assert np.all(p >= 0), "p no puede tener valores negativos"
    rango = (f >= 0.05) & (f <= 0.5)
    assert rango.any(), "las frecuencias evaluadas no cubren 0,05-0,5 Hz"
    f_pico = f[rango][np.argmax(p[rango])]
    assert abs(f_pico - 0.2) < 0.02, f"un tacograma que oscila a 0,2 Hz da el pico en {f_pico:.3f} Hz"

    f_plana = np.linspace(0, 1, 1001)
    area = periodograma_lomb.potencia_en_banda(f_plana, np.ones_like(f_plana), (0.1, 0.3))
    assert abs(area - 0.2) < 0.01, "potencia_en_banda() tiene que devolver el área bajo la curva"


MODULOS = [
    ("M0 Carga de datos", verificar_m0),
    ("M1 Acondicionamiento", verificar_m1),
    ("M2 Pan-Tompkins", verificar_m2),
    ("M3 Tacograma", verificar_m3),
    ("M4 Periodograma de Lomb", verificar_m4),
]


def revisar_modulos():
    """Imprime el estado de cada módulo y devuelve True si están todos LISTOS."""
    todos_listos = True
    for nombre, verificar in MODULOS:
        try:
            verificar()
            estado = "LISTO"
        except (NotImplementedError, FaltanRegistros) as e:
            estado = f"FALTA   {e}"
        except AssertionError as e:
            estado = f"ERROR   {e}"
        except Exception as e:
            estado = f"ERROR   {type(e).__name__}: {e}"
        todos_listos = todos_listos and estado == "LISTO"
        print(f"  {nombre:26s}{estado}")
    return todos_listos
