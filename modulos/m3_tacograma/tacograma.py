"""Módulo 3 — Tacograma y frecuencia cardíaca.

Consigna: ítems 6a (cómo varía la FC en cada registro) y 6b (comparación
entre voluntarios).

Qué recibe : las posiciones de los picos R (salida del módulo 2), fs y,
             opcionalmente, la máscara de muestras válidas del módulo 1.
Qué devuelve: t_rr y rr, en SEGUNDOS.
              rr[i]   = tiempo entre el latido i y el latido i+1
              t_rr[i] = instante del latido i+1 (el que cierra el intervalo)
              Sin máscara tienen un elemento menos que los picos; con máscara,
              además faltan los intervalos que tocan zonas inválidas.
Todo en segundos; pasar a milisegundos sólo para los gráficos.

Funciones:
  - intervalos_rr()        arma el tacograma (la que usan los demás módulos).
  - frecuencia_cardiaca()  FC instantánea en lpm.
  - rr_atipicos()          marca los RR que no pueden ser fisiológicos
                           (latido perdido o detección falsa). La usa también M4.
  - resumen()              índices del tacograma para la tabla del informe.
  - figura_tacogramas()    figura del informe con los tacogramas y la FC.

Por qué el tacograma puede tener huecos: un intervalo con alguna muestra
inválida (primer segundo del registro, tramo saturado) no se usa. No alcanza
con mirar los dos latidos de los extremos: si el latido de en medio cayó en la
saturación y no se detectó, el intervalo vale el doble aunque sus extremos
estén en zona válida. Los huecos no se rellenan: el periodograma de Lomb (M4)
no necesita muestras equiespaciadas, y los índices de este módulo tienen en
cuenta qué intervalos son consecutivos (ver son_consecutivos()).
"""

import numpy as np

# ------------------------------ Parámetros ------------------------------

VENTANA_ATIPICOS = 11   # latidos. Ancho de la mediana móvil de referencia.
UMBRAL_ATIPICOS = 0.30  # Un RR es atípico si se aparta más de un 30 % de la
                        # mediana de sus vecinos. Un latido perdido duplica el
                        # RR (+100 %) y uno falso lo parte en dos (-50 %). Con
                        # 20 % se perdían crestas reales del reposo de V1 (+28 %).
SEGUNDOS_EXTREMOS = 10.0  # s. La FC "al inicio" y "al final" del registro es
                          # la media de los RR de los primeros / últimos 10 s.
UMBRAL_NN50 = 0.050     # s. Para pNN50: diferencias sucesivas de más de 50 ms
                        # (Task Force, 1996).


# ----------------------------- Tacograma -----------------------------

def intervalos_rr(picos_r, fs, valida=None):
    """Arma el tacograma (t_rr, rr) a partir de los picos R.

    picos_r : números de muestra de los picos R (ordenados).
    fs      : frecuencia de muestreo (Hz).
    valida  : opcional, booleano por muestra (acondicionamiento.mascara_valida).
              Si se pasa, se descartan los intervalos con alguna muestra
              inválida entre sus dos latidos (extremos incluidos).
    """
    picos_r = np.asarray(picos_r, dtype=int)
    t = picos_r / fs                 # instante de cada latido (s)
    rr, t_rr = np.diff(t), t[1:]     # el RR se ubica en el latido que lo cierra
    if valida is None:
        return t_rr, rr

    # Cantidad acumulada de muestras inválidas: invalidas[k] = cuántas hay
    # antes de la muestra k. Así, cuántas hay entre a y b (inclusive) es
    # invalidas[b + 1] - invalidas[a], sin recorrer cada intervalo.
    invalidas = np.r_[0, np.cumsum(~np.asarray(valida, dtype=bool))]
    ok = invalidas[picos_r[1:] + 1] - invalidas[picos_r[:-1]] == 0
    return t_rr[ok], rr[ok]


def frecuencia_cardiaca(rr):
    """Frecuencia cardíaca instantánea en latidos por minuto (rr en segundos)."""
    return 60.0 / np.asarray(rr, dtype=float)


def rr_atipicos(rr, ventana=VENTANA_ATIPICOS, umbral=UMBRAL_ATIPICOS):
    """True en los RR que se apartan más de `umbral` (fracción) de la mediana móvil.

    La mediana de los `ventana` intervalos centrados en cada RR sigue los
    cambios lentos de la FC (la recuperación post ejercicio) pero no se deja
    arrastrar por un valor aislado.
    """
    rr = np.asarray(rr, dtype=float)
    medio = ventana // 2
    referencia = np.array([np.median(rr[max(i - medio, 0): i + medio + 1]) for i in range(len(rr))])
    return np.abs(rr / referencia - 1) > umbral


def son_consecutivos(t_rr, rr, tolerancia=1e-6):
    """True en i si rr[i] y rr[i+1] son intervalos seguidos (sin hueco entre ellos).

    Dos RR son seguidos si el latido que cierra el primero abre el segundo,
    es decir, si t_rr[i+1] - t_rr[i] == rr[i+1]. Tiene un elemento menos que rr.
    """
    t_rr, rr = np.asarray(t_rr, dtype=float), np.asarray(rr, dtype=float)
    return np.abs(np.diff(t_rr) - rr[1:]) < tolerancia


# ----------------------------- Índices -----------------------------

def resumen(t_rr, rr):
    """Valores descriptivos del tacograma para la tabla del informe (ítems 6a y 6b).

    Se calculan sin los RR atípicos. Índices temporales de la Task Force (1996):
      SDNN   desvío estándar de los RR (variabilidad total, incluye la tendencia)
      RMSSD  raíz del promedio de las diferencias sucesivas al cuadrado
             (variabilidad latido a latido, ligada al tono vagal)
      pNN50  % de diferencias sucesivas mayores a 50 ms
    Las diferencias sucesivas se toman sólo entre intervalos consecutivos y
    no atípicos. La pendiente es la de la recta que mejor ajusta la FC
    instantánea en función del tiempo (lpm por minuto).
    """
    t_rr, rr = np.asarray(t_rr, dtype=float), np.asarray(rr, dtype=float)
    atipicos = rr_atipicos(rr)
    t, x = t_rr[~atipicos], rr[~atipicos]
    fc = frecuencia_cardiaca(x)

    # Diferencias sucesivas: ambos RR no atípicos y sin hueco entre ellos
    seguidos = son_consecutivos(t_rr, rr) & ~atipicos[1:] & ~atipicos[:-1]
    dif = np.diff(rr)[seguidos]

    inicio = t <= t[0] + SEGUNDOS_EXTREMOS
    fin = t >= t[-1] - SEGUNDOS_EXTREMOS
    return {
        "N_RR": len(x),
        "atipicos": int(atipicos.sum()),
        "duracion_s": t[-1] - t[0],
        "FC_media_lpm": 60.0 / x.mean(),
        "FC_min_lpm": fc.min(),
        "FC_max_lpm": fc.max(),
        "FC_inicio_lpm": 60.0 / x[inicio].mean(),
        "FC_fin_lpm": 60.0 / x[fin].mean(),
        "pendiente_FC_lpm_min": 60.0 * np.polyfit(t, fc, 1)[0],
        "RR_medio_ms": 1000 * x.mean(),
        "SDNN_ms": 1000 * x.std(ddof=1),
        "RMSSD_ms": 1000 * np.sqrt(np.mean(dif ** 2)),
        "pNN50_%": 100 * np.mean(np.abs(dif) > UMBRAL_NN50),
    }


# ----------------------------- Figura -----------------------------

def figura_tacogramas(resultados, ventana_media=10.0):
    """Tacograma (arriba) y frecuencia cardíaca instantánea (abajo) de todos los registros.

    resultados: objetos con .registro, .t_rr y .rr (procesar_todos()).
    Todos en los mismos ejes, con el tiempo desde el inicio del registro, para
    comparar directamente reposo y post ejercicio. Los RR atípicos se marcan
    con una cruz. En la FC se superpone una media móvil de `ventana_media`
    segundos (trazo grueso) para que se vea la tendencia.
    """
    import matplotlib.pyplot as plt
    from modulos.herramientas import graficos

    cm = 1 / 2.54
    fig, (ax_rr, ax_fc) = plt.subplots(2, 1, sharex=True, figsize=(16 * cm, 10 * cm))
    for res in resultados:
        reg = res.registro
        est = graficos.estilo_registro(reg.voluntario, reg.condicion)
        atipicos = rr_atipicos(res.rr)
        ax_rr.plot(res.t_rr, 1000 * res.rr, marker=".", ms=2, lw=0.7, **est)
        if atipicos.any():
            ax_rr.plot(res.t_rr[atipicos], 1000 * res.rr[atipicos], "kx", ms=4)

        t, fc = res.t_rr[~atipicos], frecuencia_cardiaca(res.rr[~atipicos])
        ax_fc.plot(t, fc, color=est["color"], ls=est["linestyle"], lw=0.5, alpha=0.5)
        media = np.array([fc[np.abs(t - ti) <= ventana_media / 2].mean() for ti in t])
        ax_fc.plot(t, media, color=est["color"], ls=est["linestyle"], lw=1.6, label=est["label"])

    ax_rr.set(ylabel="RR (ms)", title="Tacograma")
    ax_fc.set(ylabel="FC (lpm)", xlabel="Tiempo desde el inicio del registro (s)",
              title=f"Frecuencia cardíaca instantánea (trazo fino) y media móvil de {ventana_media:g} s (grueso)")
    ax_fc.set_xlim(0, None)
    ax_rr.legend(loc="center right")
    fig.tight_layout()
    return fig
