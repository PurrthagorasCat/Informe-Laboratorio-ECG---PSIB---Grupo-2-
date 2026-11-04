"""Módulo 5 — Resultados finales.

Consigna: figuras y tablas finales de los ítems 6 y 7; material para el ítem 8.

Qué recibe: la lista de resultados de los registros (cada uno con picos_r,
            t_rr, rr, f, p; ver herramientas/procesamiento_completo.py).
Qué hace  : genera todas las figuras de resultados del informe y las tablas
            resumen. Las figuras propias de cada módulo están programadas en
            el .py de ese módulo; acá sólo se las llama y se guardan.

Figuras (informe/figuras/*.pdf):
    pt_etapas            etapas de Pan-Tompkins sobre el registro de reposo (M2)
    deteccion_picos      picos R detectados sobre la señal acondicionada
    tacogramas           tacograma y FC instantánea de todos los registros (M3)
    lomb_periodogramas   periodograma de Lomb de cada registro (M4)
    lomb_ventanas_post   periodograma por minuto de los post ejercicio (M4)

Tablas (informe/figuras/*.csv; los valores se copian a las tablas del informe):
    tabla_deteccion      latidos detectados e intervalos usados
    tabla_tacograma      índices temporales (ítem 6)
    tabla_lomb           índices espectrales de cada registro y por minuto (ítem 7)
"""

import numpy as np
import pandas as pd

from modulos.herramientas import graficos
from modulos.m2_pan_tompkins import pan_tompkins
from modulos.m3_tacograma import tacograma
from modulos.m4_periodograma_lomb import periodograma_lomb as lomb

DURACION_DETALLE = 8.0   # s que se muestran de cada registro en deteccion_picos
TRAMO_ETAPAS = (2.0, 8.0)   # s del registro de reposo que se muestran en pt_etapas


# ------------------------------ Tablas ------------------------------

def tabla_deteccion(resultados):
    """Cuántos latidos se detectaron y cuántos intervalos quedaron para el análisis."""
    filas = {}
    for res in resultados:
        reg = res.registro
        filas[reg.nombre] = {
            "duracion_s": reg.duracion,
            "valida_%": 100 * res.valida.mean(),
            "latidos": len(res.picos_r),
            "RR_posibles": len(res.picos_r) - 1,
            "RR_en_zona_valida": len(res.rr),
            "RR_atipicos": int(tacograma.rr_atipicos(res.rr).sum()),
        }
    return pd.DataFrame.from_dict(filas, orient="index")


def tabla_tacograma(resultados):
    """tacograma.resumen() de cada registro (ítems 6a y 6b)."""
    return pd.DataFrame.from_dict({res.registro.nombre: tacograma.resumen(res.t_rr, res.rr)
                                   for res in resultados}, orient="index")


def tabla_lomb(resultados):
    """Índices espectrales de cada registro completo y, para los post ejercicio, por minuto."""
    filas = []
    for res in resultados:
        reg = res.registro
        validos = ~tacograma.rr_atipicos(res.rr)
        filas.append({"registro": reg.nombre, "tramo": "completo",
                      "n_rr": int(validos.sum()), "T_s": res.t_rr[-1] - res.t_rr[0],
                      "FC_lpm": 60 / res.rr[validos].mean(), **lomb.indices_espectrales(res.f, res.p)})
        if reg.condicion == "post_ejercicio":
            for k, fila in enumerate(lomb.indices_por_ventanas(res.t_rr, res.rr)):
                sel = (res.t_rr >= res.t_rr[0] + fila["inicio_s"]) & (res.t_rr < res.t_rr[0] + fila["inicio_s"] + 60)
                filas.append({"registro": reg.nombre, "tramo": f"minuto {k + 1}",
                              "T_s": res.t_rr[sel][-1] - res.t_rr[sel][0],
                              **{c: v for c, v in fila.items() if c != "inicio_s"}})
    return pd.DataFrame(filas).set_index(["registro", "tramo"])


# ------------------------------ Figuras ------------------------------

def figura_deteccion(resultados, duracion=DURACION_DETALLE):
    """Señal acondicionada con los picos R detectados, un panel por registro.

    En los registros con tramos saturados se muestra el tramo centrado en el
    primer tramo saturado después de los 100 s, para que se vea cómo se excluye
    la zona; si no hay, el tramo del medio del registro. Sombreado: muestras
    inválidas (mascara_valida).
    """
    import matplotlib.pyplot as plt
    from modulos.m1_acondicionamiento import acondicionamiento

    cm = 1 / 2.54
    fig, axs = plt.subplots(len(resultados), 1, figsize=(16 * cm, 4.5 * cm * len(resultados)), squeeze=False)
    for ax, res in zip(axs[:, 0], resultados):
        reg = res.registro
        tramos = [tr for tr in acondicionamiento.tramos_saturados(reg.ecg, reg.fs) if tr[0] / reg.fs > 100]
        centro = (tramos[0][0] + tramos[0][1]) / 2 / reg.fs if tramos else reg.duracion / 2
        t0 = centro - duracion / 2
        t = np.arange(len(reg.ecg)) / reg.fs
        m = (t >= t0) & (t < t0 + duracion)
        picos = res.picos_r[(res.picos_r >= t0 * reg.fs) & (res.picos_r < (t0 + duracion) * reg.fs)]

        est = graficos.estilo_registro(reg.voluntario, reg.condicion)
        ax.plot(t[m], res.ecg_acondicionado[m], lw=0.7, **est)
        ax.plot(t[picos], res.ecg_acondicionado[picos], "o", ms=3.5, mfc="none", mec="k", label="pico R detectado")
        if (~res.valida[m]).any():
            ax.fill_between(t[m], 0, 1, where=~res.valida[m], color="0.85", lw=0,
                            transform=ax.get_xaxis_transform(), label="zona excluida")
        ax.set(ylabel=f"ECG ({reg.unidades})", xlim=(t0, t0 + duracion))
        ax.legend(loc="lower right", bbox_to_anchor=(1, 1), ncol=3, borderaxespad=0.2)
    axs[-1, 0].set_xlabel("Tiempo (s)")
    fig.tight_layout()
    return fig


# ------------------------------ Todo junto ------------------------------

def generar(resultados):
    """Lo llama main.py al final, con los registros ya procesados.

    Guarda las figuras y las tablas en informe/figuras/ y devuelve las tablas
    en un diccionario {nombre: DataFrame} para mostrarlas.
    """
    import matplotlib.pyplot as plt

    graficos.aplicar_estilo()

    figuras = {}
    reposo = next((r for r in resultados if r.registro.condicion == "reposo"), resultados[0])
    etapas = pan_tompkins.pan_tompkins(reposo.ecg_acondicionado, reposo.registro.fs)
    figuras["pt_etapas"] = pan_tompkins.figura_etapas(
        reposo.ecg_acondicionado, reposo.registro.fs, etapas, *TRAMO_ETAPAS,
        titulo=f"{reposo.registro.voluntario} - {reposo.registro.condicion.replace('_', ' ')}")
    figuras["deteccion_picos"] = figura_deteccion(resultados)
    figuras["tacogramas"] = tacograma.figura_tacogramas(resultados)
    figuras["lomb_periodogramas"] = lomb.figura_periodogramas(resultados)
    if any(r.registro.condicion == "post_ejercicio" for r in resultados):
        figuras["lomb_ventanas_post"] = lomb.figura_ventanas(resultados)

    for nombre, fig in figuras.items():
        graficos.guardar_figura(fig, nombre)
        plt.close(fig)

    tablas = {
        "tabla_deteccion": tabla_deteccion(resultados),
        "tabla_tacograma": tabla_tacograma(resultados),
        "tabla_lomb": tabla_lomb(resultados),
    }
    for nombre, tabla in tablas.items():
        tabla.to_csv(graficos.CARPETA_FIGURAS / f"{nombre}.csv", float_format="%.4g")
    return {"figuras": list(figuras), "tablas": tablas}
