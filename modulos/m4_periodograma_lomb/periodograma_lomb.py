"""Módulo 4 — Periodograma de Lomb del tacograma.

Consigna: ítems 7a (reposo vs. actividad física) y 7b (comparación entre
voluntarios).
Referencia: Lomb, N. R. (1976). Astrophysics and Space Science, 39, 447-462.

Qué recibe : t_rr y rr (salida del módulo 3, en segundos).
Qué devuelve: f (frecuencias en Hz) y p (potencia), de igual largo, p >= 0.

Por qué Lomb: el tacograma está muestreado en los latidos, que no son
equiespaciados, y además le faltan los latidos que M3 descarta (inicio del
registro, tramos saturados). Lomb ajusta senoidales por cuadrados mínimos
directamente sobre esos instantes, sin interpolar ni remuestrear.

Etapas de periodograma_lomb() (justificación en el informe, sección 2.6):
  1. preparar_tacograma     descarta latidos atípicos y quita la tendencia
  2. grilla_de_frecuencias  si no se pasa f: de df a F_MAX, con df = 1/(SOBREMUESTREO*T)
  3. lomb                   periodograma de Lomb, pasado a densidad espectral en s²/Hz

Normalización de p: densidad espectral unilateral en s²/Hz, escalada para que
el área bajo la curva sea la varianza de la serie (Parseval). Así
potencia_en_banda() devuelve potencias en s² (x 1e6 para ms², como en la
bibliografía de VFC).
"""

import numpy as np

# ------------------------------ Parámetros ------------------------------

# Bandas de frecuencia para el análisis: {"NOMBRE": (f_min, f_max)} en Hz.
# Task Force (1996), registros cortos. La banda VLF (< 0,04 Hz) no se incluye:
# la Task Force desaconseja estimarla en registros de menos de 5 min.
BANDAS = {
    "LF": (0.04, 0.15),
    "HF": (0.15, 0.40),
}

F_MAX = 0.5              # Hz. Cubre HF (hasta 0,4 Hz) con margen. Es la "Nyquist
                         # media" de un tacograma a 60 lpm (un latido por segundo).
SOBREMUESTREO = 4        # La resolución natural es 1/T; se evalúa 4 veces más
                         # fino para ubicar bien los picos (no agrega resolución).
VENTANA_ATIPICOS = 11    # latidos. Ancho de la mediana móvil de referencia.
UMBRAL_ATIPICOS = 0.30   # Se descarta un RR que se aparta más de 30 % de la
                         # mediana de sus vecinos. Un latido perdido duplica el
                         # RR (+100 %) y uno falso lo parte en dos. Con 20 % se
                         # perdían crestas reales del reposo de V1 (+28 %).
GRADO_TENDENCIA = 2      # Grado del polinomio que se resta como tendencia. Quita
                         # la recuperación lenta post ejercicio sin tocar LF.


# ---------------------------- Función principal ----------------------------

def periodograma_lomb(t, x, f=None, preparar=True):
    """Periodograma de Lomb de una serie con muestras no equiespaciadas.

    t, x     : instantes (s) y valores (s) del tacograma.
    f        : frecuencias a evaluar (Hz). Si no se pasa, grilla_de_frecuencias(t).
    preparar : si es True, aplica preparar_tacograma() antes (atípicos y tendencia).
    Devuelve f (Hz) y p (s²/Hz).
    """
    t = np.asarray(t, dtype=float)
    x = np.asarray(x, dtype=float)
    if preparar:
        t, x = preparar_tacograma(t, x)
    if f is None:
        f = grilla_de_frecuencias(t)
    f = np.asarray(f, dtype=float)
    return f, lomb(t, x, f)


# --------------------------------- Etapas ---------------------------------

def preparar_tacograma(t, rr):
    """Descarta los RR atípicos y resta la tendencia. Devuelve (t, rr sin tendencia).

    Con Lomb los latidos atípicos se pueden quitar sin más: no hace falta
    rellenar el hueco, porque el método no necesita muestras equiespaciadas.
    """
    ok = ~rr_atipicos(rr)
    t, rr = t[ok], rr[ok]
    return t, quitar_tendencia(t, rr)


def rr_atipicos(rr, ventana=VENTANA_ATIPICOS, umbral=UMBRAL_ATIPICOS):
    """True en los RR que se apartan más de `umbral` (fracción) de la mediana móvil."""
    rr = np.asarray(rr, dtype=float)
    medio = ventana // 2
    referencia = np.array([np.median(rr[max(i - medio, 0): i + medio + 1]) for i in range(len(rr))])
    return np.abs(rr / referencia - 1) > umbral


def quitar_tendencia(t, x, grado=GRADO_TENDENCIA):
    """Resta un polinomio de grado `grado` ajustado por cuadrados mínimos en función de t.

    Con grado 0 sólo se resta la media, que Lomb necesita igual.
    """
    tc = t - t.mean()   # centrar t mejora el condicionamiento del ajuste
    return x - np.polyval(np.polyfit(tc, x, grado), tc)


def grilla_de_frecuencias(t, f_max=F_MAX, sobremuestreo=SOBREMUESTREO):
    """Frecuencias equiespaciadas de df a f_max, con df = 1 / (sobremuestreo * T).

    T es la duración del tacograma. No se incluye f = 0 (la media ya se quitó).
    """
    df = 1.0 / (sobremuestreo * (t[-1] - t[0]))
    return np.arange(df, f_max + df / 2, df)


def lomb(t, x, f):
    """Periodograma de Lomb (1976) de x(t) evaluado en f, como densidad unilateral.

    Para cada frecuencia se ajustan un coseno y un seno por cuadrados mínimos;
    el corrimiento tau hace que los dos ajustes sean independientes:
        tan(2 w tau) = sum sin(2 w t) / sum cos(2 w t)
        P(w) = 1/2 [ (sum x cos w(t-tau))² / sum cos² w(t-tau)
                   + (sum x sin w(t-tau))² / sum sin² w(t-tau) ]
    Con muestras equiespaciadas P coincide con el periodograma clásico |X|²/N.
    Se pasa a densidad unilateral con p = 2 P / fs_media, donde fs_media = N / T
    (latidos por segundo), para que el área bajo p sea la varianza de x.
    """
    x = x - x.mean()
    w = 2 * np.pi * f[:, None]   # una fila por frecuencia, una columna por muestra
    tau = np.arctan2(np.sin(2 * w * t).sum(axis=1), np.cos(2 * w * t).sum(axis=1)) / (2 * w[:, 0])
    fase = w * (t - tau[:, None])
    c, s = np.cos(fase), np.sin(fase)
    P = 0.5 * ((c @ x) ** 2 / (c ** 2).sum(axis=1) + (s @ x) ** 2 / (s ** 2).sum(axis=1))
    fs_media = len(t) / (t[-1] - t[0])
    return 2 * P / fs_media


# ---------------------------- Índices espectrales ----------------------------

def potencia_en_banda(f, p, banda):
    """Área bajo el periodograma entre banda[0] y banda[1] (en Hz)."""
    f = np.asarray(f, dtype=float)
    p = np.asarray(p, dtype=float)
    en_banda = (f >= banda[0]) & (f <= banda[1])
    return np.trapezoid(p[en_banda], f[en_banda])


def frecuencia_pico(f, p, banda):
    """Frecuencia (Hz) del máximo del periodograma dentro de la banda."""
    en_banda = (f >= banda[0]) & (f <= banda[1])
    return f[en_banda][np.argmax(p[en_banda])]


def indices_espectrales(f, p):
    """Valores para la tabla del informe (ítems 7a y 7b).

    LF y HF en ms²; LFnu y HFnu en % de LF + HF (Task Force, 1996); picos en Hz.
    El pico de HF estima la frecuencia respiratoria.
    """
    lf = potencia_en_banda(f, p, BANDAS["LF"]) * 1e6
    hf = potencia_en_banda(f, p, BANDAS["HF"]) * 1e6
    return {
        "LF_ms2": lf,
        "HF_ms2": hf,
        "LF_HF": lf / hf,
        "LFnu_%": 100 * lf / (lf + hf),
        "HFnu_%": 100 * hf / (lf + hf),
        "pico_LF_Hz": frecuencia_pico(f, p, BANDAS["LF"]),
        "pico_HF_Hz": frecuencia_pico(f, p, BANDAS["HF"]),
    }


def indices_por_ventanas(t, rr, duracion=60.0):
    """indices_espectrales() en ventanas consecutivas de `duracion` segundos.

    Sirve para seguir la recuperación post ejercicio, donde el tacograma no es
    estacionario. Cada ventana se trata como un tacograma aparte (atípicos,
    tendencia y grilla propios). Devuelve una lista de diccionarios con
    además el inicio de la ventana (s), la cantidad de RR y la FC media (lpm).
    """
    t = np.asarray(t, dtype=float)
    rr = np.asarray(rr, dtype=float)
    filas = []
    for t0 in np.arange(t[0], t[-1] - duracion / 2, duracion):   # se descarta una última ventana de menos de la mitad
        sel = (t >= t0) & (t < t0 + duracion)
        f, p = periodograma_lomb(t[sel], rr[sel])
        validos = rr[sel][~rr_atipicos(rr[sel])]
        filas.append({"inicio_s": t0 - t[0], "n_rr": len(validos), "FC_lpm": 60 / validos.mean(),
                      **indices_espectrales(f, p)})
    return filas


# --------------------------------- Figura ---------------------------------

# Relleno de cada banda en la figura. VLF se muestra pero no se analiza (ver BANDAS).
COLOR_BANDA = {"VLF": "#c8c8c4", "LF": "#e8836f", "HF": "#8ec5ee"}


def figura_periodogramas(resultados):
    """Un periodograma por registro, en escala lineal, con el área de cada banda pintada.

    resultados: objetos con .registro (voluntario, condicion), .f y .p, como
    los que devuelve herramientas.procesamiento_completo.procesar_todos().
    Cada registro tiene su propio eje y: entre reposo y post ejercicio la
    potencia cambia en dos órdenes de magnitud y con un eje común el post
    quedaría plano.
    """
    import matplotlib.pyplot as plt

    cm = 1 / 2.54
    fig, axs = plt.subplots(len(resultados), 1, figsize=(16 * cm, 5 * cm * len(resultados)),
                            sharex=True, squeeze=False)
    zonas = {"VLF": (0.0, BANDAS["LF"][0]), **BANDAS}
    for ax, res in zip(axs[:, 0], resultados):
        reg = res.registro
        p = res.p * 1e6   # s²/Hz -> ms²/Hz
        for nombre, (f0, f1) in zonas.items():
            sel = (res.f >= f0) & (res.f <= f1)
            etiqueta = f"{nombre}: {f0:.2f}–{f1:.2f} Hz" + (" (no se analiza)" if nombre not in BANDAS else "")
            ax.fill_between(res.f[sel], p[sel], color=COLOR_BANDA[nombre], lw=0, label=etiqueta)
        ax.plot(res.f, p, color="k", lw=0.6)
        for f0 in (BANDAS["LF"][0], BANDAS["HF"][0], BANDAS["HF"][1]):
            ax.axvline(f0, color="0.4", ls="--", lw=0.5)
        ax.set(title=f"{reg.voluntario} - {reg.condicion.replace('_', ' ')}", ylabel="DEP (ms²/Hz)", ylim=(0, None))
        ax.legend(loc="upper right")
    axs[-1, 0].set(xlim=(0, F_MAX), xlabel="Frecuencia (Hz)")
    fig.tight_layout()
    return fig
