"""Módulo 2 — Detección de latidos con el método de Pan-Tompkins.

Consigna: ítem 6 (detección de latidos).
Referencia: Pan, J., & Tompkins, W. J. (1985). IEEE Trans. Biomed. Eng., 32(3), 230-236.

Qué recibe : la señal acondicionada (salida del módulo 1), fs y, opcionalmente, el umbral relativo (por defecto 0.3).
Qué devuelve: las posiciones (números de muestra, enteros) de los picos R, ordenadas, sin repetidos y referidas a la señal que entró.

Funciones:
  - detectar_picos_r(ecg, fs, umbral_relativo): devuelve sólo los picos R (la que usan los demás módulos).
  - pan_tompkins(ecg, fs, umbral_relativo): devuelve además las señales de cada etapa, el umbral y los QRS
                                            sobre la señal integrada (para las figuras y para encontrar fallas).

Qué hace (etapas):

  1. Pasabanda de 5 a 15 Hz: deja pasar la energía del QRS y atenúa la deriva de la línea de base,
                             las ondas P y T, la interferencia de línea y el ruido muscular.

  2. Derivada de cinco puntos: resalta las pendientes pronunciadas del QRS.

  3. Cuadrado: convierte todos los lóbulos en positivos y acentúa las pendientes grandes.

  4. Integración con ventana móvil de 150 ms: junta en un solo pulso la energía de cada QRS,
                                              y la distingue de las ondas P y T.

  5. Decisión: un máximo local de la señal integrada es un QRS si supera un umbral fijo (una fracción
               del máximo de la integrada) y está a más de 200 ms (período refractario) del anterior.

  6. Ubicación del pico R: en la señal de entrada, se busca el extremo dentro de ±100 ms alrededor de cada QRS.
                           Antes se decide la polaridad del registro, así también funciona si los complejos
                           están invertidos (electrodos al revés).

Diferencias con el algoritmo del artículo:
  - Usamos un único umbral fijo, en lugar de los dos umbrales adaptativos (sobre la señal filtrada y la integrada)
    y la búsqueda hacia atrás de latidos perdidos. Así podemos estudiar cómo cambia la detección al variar el umbral.
  - El pasabanda es un Butterworth (IIR) de orden 2 aplicado ida y vuelta (filtfilt); como el análisis no es en tiempo real,
    usar filtfilt no trae problemas. Por otro lado, la derivada y la ventana móvil están centradas, así no tenemos que corregir el retardo.
"""

import numpy as np
from scipy.signal import butter, filtfilt, find_peaks

# --- Parámetros --------------
BANDA_HZ = (5.0, 15.0)            # pasabanda (en Hz)
ORDEN_PASABANDA = 2               # orden del Butterworth (se aplica dos veces)
VENTANA_INTEGRACION_S = 0.150     # ancho de la ventana móvil (en segundos)
PERIODO_REFRACTARIO_S = 0.200     # no puede haber dos latidos más cerca que esto (en segundos)
SEMIANCHO_BUSQUEDA_S = 0.100      # alrededor de cada QRS se busca el pico R en ±100 ms (en segundos)
UMBRAL_RELATIVO = 0.3             # umbral, como fracción del máximo de la señal integrada (porcentual)


def detectar_picos_r(ecg, fs, umbral_relativo=UMBRAL_RELATIVO):
    """Función que devuelve la ubicación temporal de los picos R."""
    return pan_tompkins(ecg, fs, umbral_relativo)["picos_r"]


def pan_tompkins(ecg, fs, umbral_relativo=UMBRAL_RELATIVO):
    """Corre todo el algoritmo y devuelve las señales de cada etapa (por si hubo una falla, poder identificarla)."""
    ecg = np.asarray(ecg, dtype=float)
    filtrada = filtro_pasabanda(ecg, fs)
    derivada = derivar(filtrada, fs)
    al_cuadrado = elevar_al_cuadrado(derivada)
    integrada = integrar_con_ventana_movil(al_cuadrado, fs)
    qrs, umbral = decidir_picos(integrada, fs, umbral_relativo)
    picos_r = ubicar_picos_r(ecg, qrs, fs)
    return {
        "filtrada": filtrada,
        "derivada": derivada,
        "al_cuadrado": al_cuadrado,
        "integrada": integrada,
        "umbral_integrada": umbral,
        "qrs_integrada": qrs,
        "picos_r": picos_r,
    }


# --- Etapas -----------------------------------------------------------

def filtro_pasabanda(x, fs):
    b, a = butter(ORDEN_PASABANDA, BANDA_HZ, btype="bandpass", fs=fs) # Butterworth pasabanda de 5 a 15 Hz.
    return filtfilt(b, a, x) # Lo aplicamos ida y vuelta (para sacar el desfasaje)


def derivar(x, fs):
    """Derivada de cinco puntos, centrada para no introducir retardo y después tener que compensar"""

    nucleo = np.array([1.0, 2.0, 0.0, -2.0, -1.0]) * fs / 8
    return np.convolve(x, nucleo, mode="same") # Convolucionamos la ventana con la señal, en el plano frecuencial.


def elevar_al_cuadrado(x):
    return np.asarray(x) ** 2


def integrar_con_ventana_movil(x, fs):
    n = max(int(round(VENTANA_INTEGRACION_S * fs)), 1) # Promedio móvil centrado de 150 ms.
    return np.convolve(x, np.ones(n) / n, mode="same") # Convolucionamos con la señal


# --- Etapa 5: decisión ------------------------------------------------------

def decidir_picos(integrada, fs, umbral_relativo=UMBRAL_RELATIVO):
    """Un máximo local de la señal integrada es un QRS si supera el umbral fijo
    (umbral_relativo * máximo de la integrada) y está a más de 200 ms del anterior.
    Devuelve (qrs, umbral), que son las posiciones de los QRS en la señal integrada y el umbral. """

    refractario = max(int(round(PERIODO_REFRACTARIO_S * fs)), 1)
    valor_umbral = umbral_relativo * integrada.max()
    qrs, _ = find_peaks(integrada, height=valor_umbral, distance=refractario)
    return qrs.astype(int), np.full(len(integrada), valor_umbral)


def ubicar_picos_r(ecg, qrs, fs):
    """Busca el pico R en la señal de entrada, en ±100 ms alrededor de cada QRS."""

    ecg = np.asarray(ecg, dtype=float)
    if len(qrs) == 0:
        return np.array([], dtype=int)
    n = len(ecg)
    semiancho = int(round(SEMIANCHO_BUSQUEDA_S * fs))
    tramos = [(max(q - semiancho, 0), min(q + semiancho + 1, n)) for q in qrs]

    # Decide la polaridad del registro. En cada ventana mide cuánto sube la señal por encima de su mediana y cuánto baja
    # por debajo. Si en general baja más que lo que sube, los complejos están invertidos (electrodos al revés) y signo vale −1.
    hacia_arriba, hacia_abajo = [], []
    for i0, i1 in tramos:
        tramo = ecg[i0:i1]
        base = np.median(tramo)
        hacia_arriba.append(tramo.max() - base)
        hacia_abajo.append(base - tramo.min())
    signo = 1.0 if np.median(hacia_arriba) >= np.median(hacia_abajo) else -1.0

    picos = [i0 + int(np.argmax(signo * ecg[i0:i1])) for i0, i1 in tramos] # Ubica los máximos dentro de cada ventana
    return np.unique(np.array(picos, dtype=int)) # ordena los picos y elimina repetidos, por si dos ventanas encontraron el mismo pico.