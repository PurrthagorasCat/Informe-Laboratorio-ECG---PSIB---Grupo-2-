"""Señales inventadas, con los latidos en posiciones conocidas, para probar los
módulos sin depender de los registros reales.

    generar_latidos()  instantes de los latidos de un tacograma que oscila a una
                       frecuencia conocida (para probar M4).
    generar_ecg()      un ECG sintético con esos latidos, deriva de la línea de
                       base y ruido, más las posiciones exactas de sus picos R
                       (para probar M1, M2 y M3).

Lo usa herramientas/verificacion.py, así que no conviene cambiar los valores
por defecto.
"""

import numpy as np

# Ondas de cada latido como gaussianas: (nombre, retardo respecto de R [s],
# ancho [s], amplitud [mV]). La R está en el retardo 0, así su pico cae
# exactamente en el instante del latido.
ONDAS = (
    ("P", -0.200, 0.025, 0.15),
    ("Q", -0.030, 0.008, -0.15),
    ("R", 0.000, 0.010, 1.00),
    ("S", 0.030, 0.008, -0.25),
    ("T", 0.280, 0.050, 0.30),
)


def generar_latidos(duracion=60.0, fc_media=70.0, amp_mod=0.05, f_mod=0.25, t0=0.0):
    """Instantes (s) de los latidos de un tacograma que oscila a f_mod Hz.

    Cada latido llega un RR después del anterior, con
        RR = (60 / fc_media) * (1 + amp_mod * sin(2 pi f_mod t)),
    donde t es el instante del latido que abre el intervalo.
    """
    rr_medio = 60.0 / fc_media
    t = [t0]
    while t[-1] < t0 + duracion:
        t.append(t[-1] + rr_medio * (1 + amp_mod * np.sin(2 * np.pi * f_mod * t[-1])))
    return np.array(t[:-1])   # el último ya quedó fuera de la duración pedida


def generar_ecg(fs=500.0, duracion=20.0, fc_media=70.0, amp_mod=0.05, f_mod=0.25,
                ruido_std=0.01, deriva=0.2, semilla=None):
    """ECG sintético (mV) con picos R en posiciones conocidas.

    Devuelve (t, ecg, picos_r): tiempo en s, la señal y los números de muestra
    de los picos R. Incluye una deriva de la línea de base de 0,15 Hz (amplitud
    `deriva` en mV), un nivel de continua y ruido blanco gaussiano de desvío
    `ruido_std` (mV).
    """
    rng = np.random.default_rng(semilla)
    n = int(round(duracion * fs))
    t = np.arange(n) / fs

    # Latidos: se deja medio segundo libre en cada borde para que todas las
    # ondas de cada latido queden dentro de la señal.
    latidos = generar_latidos(duracion - 1.0, fc_media, amp_mod, f_mod, t0=0.5)
    picos_r = np.round(latidos * fs).astype(int)

    ecg = np.zeros(n)
    for tr in picos_r / fs:   # se usa el instante redondeado, así R cae justo en la muestra
        for _, retardo, ancho, amp in ONDAS:
            ecg += amp * np.exp(-0.5 * ((t - tr - retardo) / ancho) ** 2)

    ecg += deriva * np.sin(2 * np.pi * 0.15 * t + rng.uniform(0, 2 * np.pi))
    ecg += 0.5   # nivel de continua
    ecg += ruido_std * rng.standard_normal(n)
    return t, ecg, picos_r
