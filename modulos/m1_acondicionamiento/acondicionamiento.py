"""Módulo 1 — Acondicionamiento de la señal.

Consigna: ítems 4 (¿hace falta acondicionar? ¿qué ruidos asociados al
registro aparecen?) y 5 (acondicionar la señal).

Qué recibe : la señal cruda y fs.
Qué devuelve: la señal acondicionada, con la MISMA cantidad de muestras y SIN
              CORRIMIENTO en el tiempo (un pico que estaba en la muestra n
              sigue en la muestra n). Los módulos 2 y 3 dependen de esto.

Etapas de acondicionar() (justificación en el informe, sección 2.3):
  1. interpolar_saturaciones   tramos largos pegados a los extremos (0 V / 3,3 V)
  2. quitar_offset             nivel de continua (salida referida a 1,65 V)
  3. quitar_linea_de_base      deriva electrodo-piel -> pasaaltos de 0,5 Hz
  4. quitar_red_electrica      notch de 50 Hz (DESACTIVADO: ver APLICAR_NOTCH)
Los filtros se aplican hacia adelante y hacia atrás (fase cero).

Además, mascara_valida() marca qué muestras sirven para el análisis: excluye
el primer segundo (paquete viejo que guarda la interfaz) y los tramos
saturados largos. M2/M3 la pueden usar para descartar latidos ahí.

Organización del archivo:
  - Parámetros: todas las decisiones numéricas, en un solo lugar.
  - acondicionar(): la función que llaman los demás módulos.
  - Etapas: una función por cada paso del acondicionamiento.
  - Calidad del registro: funciones que detectan las zonas que no sirven
    y calculan los números que van al informe.
"""

import numpy as np
from scipy import signal

# ------------------------------ Parámetros ------------------------------
# Todas las decisiones numéricas están acá para no tener números sueltos en
# el código. Si se cambia alguno, actualizar también la tabla del informe
# (sección 2.3) y volver a correr el cuaderno para regenerar las figuras.
#
# OJO: los parámetros en segundos y en Hz dependen de que fs sea la correcta.
# Si fs está mal, un "0,1 s" o un "0,5 Hz" en realidad son otra cosa.

TOLERANCIA_SAT = 0.01   # Una muestra cuenta como saturada si está a menos del
                        # 1 % del rango del registro de su mínimo o su máximo.
DURACION_SAT_MIN = 0.1  # s. Sólo se interpolan los tramos saturados de al menos
                        # 0,1 s. Las puntas de R recortadas duran menos y se
                        # dejan como están, para no mover el pico.
MARGEN_SAT = 0.2        # s. Además del tramo saturado, se excluye este margen
                        # a cada lado (el pasaaltos deja la señal alterada
                        # alrededor del tramo interpolado).
DESCARTE_INICIO = 1.0   # s. El primer paquete que guarda la interfaz (100
                        # muestras) es viejo y no es continuo con el resto.
                        # Se descarta 1 s para tener margen.
FC_PASAALTOS = 0.5      # Hz. Corte del pasaaltos para la deriva de la línea
                        # de base. Queda por debajo de la FC más baja esperable.
ORDEN_PASAALTOS = 2     # Orden del Butterworth. Como se aplica ida y vuelta
                        # (filtfilt), el orden efectivo es el doble: 4.
APLICAR_NOTCH = False   # El espectro no muestra pico en 50 Hz (ver cuaderno),
                        # así que el notch está programado pero apagado.
F_RED, Q_NOTCH = 50.0, 30.0   # Frecuencia de la red en Argentina y factor de
                              # calidad del notch (ancho = 50/30 ≈ 1,7 Hz).


def acondicionar(ecg, fs):
    """Aplica las etapas de acondicionamiento justificadas en el ítem 4.

    Es la única función de este archivo que usan los demás módulos.
    Cada etapa recibe la salida de la anterior.
    """
    # Se convierte a float para poder operar (el CSV podría venir como enteros)
    x = np.asarray(ecg, dtype=float)

    # El orden importa: las saturaciones se arreglan ANTES de filtrar. Si no,
    # el pasaaltos "ve" los escalones de entrada y salida de la meseta y
    # genera oscilaciones alrededor.
    x = interpolar_saturaciones(x, fs)
    x = quitar_offset(x)
    x = quitar_linea_de_base(x, fs)
    if APLICAR_NOTCH:
        x = quitar_red_electrica(x, fs)
    return x


# -------------------------------- Etapas --------------------------------

def quitar_offset(x):
    """Resta la mediana para centrar la señal en cero.

    Se usa la mediana y no la media porque los tramos saturados (en 0 V o en
    3,3 V) arrastran la media hacia un lado; la mediana casi no se mueve.
    El pasaaltos de la etapa siguiente también quitaría la continua, pero
    centrar antes evita un transitorio grande en los bordes del filtro.
    """
    return x - np.median(x)


def interpolar_saturaciones(x, fs, dur_min=DURACION_SAT_MIN):
    """Reemplaza por una recta los tramos saturados de más de dur_min segundos.

    En esos tramos la señal real es desconocida; la meseta produce un escalón
    que el pasaaltos convierte en oscilaciones. Las puntas de R recortadas
    (pocas muestras) no se tocan, para no mover el pico.
    """
    # Copia, para no modificar el array original que nos pasaron
    x = x.copy()

    # tramos_saturados devuelve pares (i0, i1): el tramo va de la muestra i0
    # inclusive a la i1 exclusive, como en los slices de Python (x[i0:i1]).
    for i0, i1 in tramos_saturados(x, fs, dur_min):
        # a = último valor bueno ANTES del tramo.
        # Si el tramo empieza en la muestra 0 no hay valor anterior: se usa el
        # de después (x[i1]), y si el tramo ocupa toda la señal, la mediana.
        a = x[i0 - 1] if i0 > 0 else x[i1] if i1 < len(x) else np.median(x)
        # b = primer valor bueno DESPUÉS del tramo (o a, si el tramo llega al final)
        b = x[i1] if i1 < len(x) else a
        # linspace genera la recta de a a b incluyendo los dos extremos.
        # Como a y b ya existen en la señal, [1:-1] se queda sólo con los
        # puntos del medio, que son exactamente los i1 - i0 que hay que llenar.
        x[i0:i1] = np.linspace(a, b, i1 - i0 + 2)[1:-1]
    return x


def quitar_linea_de_base(x, fs, fc=FC_PASAALTOS, orden=ORDEN_PASAALTOS):
    """Butterworth pasaaltos de fase cero (sosfiltfilt).

    Quita la deriva lenta de la línea de base (potencial electrodo-piel).
    """
    # butter() diseña el filtro. Al pasarle fs=fs, el corte fc se da
    # directamente en Hz (si no, habría que normalizarlo por fs/2).
    # output="sos" lo devuelve como "secciones de segundo orden": un filtro de
    # orden alto partido en etapas de orden 2 en cascada. Es numéricamente más
    # estable que la forma (b, a), sobre todo con cortes tan bajos como 0,5 Hz.
    sos = signal.butter(orden, fc, btype="highpass", fs=fs, output="sos")

    # sosfiltfilt filtra hacia adelante y después hacia atrás. El desfasaje de
    # la ida se cancela con el de la vuelta, así que la fase total es cero y
    # los picos R no se mueven (condición del README). El módulo queda elevado
    # al cuadrado: por eso el orden efectivo es el doble.
    # La salida tiene la misma cantidad de muestras que la entrada.
    return signal.sosfiltfilt(sos, x)


def quitar_red_electrica(x, fs, f0=F_RED, q=Q_NOTCH):
    """Notch de fase cero en f0. No se usa por defecto (APLICAR_NOTCH)."""
    # Si f0 está por encima de Nyquist (fs/2), esa frecuencia no existe en la
    # señal muestreada: no hay nada que filtrar.
    if f0 >= fs / 2:
        return x
    # iirnotch diseña un filtro que elimina una banda angosta alrededor de f0
    # (ancho f0/q). filtfilt, igual que antes, para que la fase sea cero.
    b, a = signal.iirnotch(f0, q, fs=fs)
    return signal.filtfilt(b, a, x)


# --------------------------- Calidad del registro ---------------------------

def muestras_saturadas(x, tol=TOLERANCIA_SAT):
    """Booleano por muestra: True si está pegada al mínimo o al máximo del registro.

    Se usan los extremos del propio registro (0 V y 3,3 V en los post ejercicio)
    para que funcione igual en voltios o en cuentas del conversor. En un registro
    sin saturación los extremos son puntas aisladas y no forman tramos largos.
    """
    x = np.asarray(x, dtype=float)
    # Margen de tolerancia: 1 % del rango total. Con 0 a 3,3 V son ~0,033 V,
    # así que 3,29 V también cuenta como saturado (la interfaz redondea).
    margen = tol * (x.max() - x.min())
    # | es el "o" elemento a elemento entre arrays booleanos
    return (x <= x.min() + margen) | (x >= x.max() - margen)


def tramos_saturados(x, fs, dur_min=DURACION_SAT_MIN):
    """Lista de (inicio, fin) en muestras de los tramos saturados de más de dur_min s."""
    # Ejemplo de cómo funciona, con una señal de 7 muestras:
    #   sat           =    [0, 1, 1, 1, 0, 0, 1]
    #   r_[0, sat, 0] = [0, 0, 1, 1, 1, 0, 0, 1, 0]   (se agrega un 0 en cada punta)
    #   diff          =    [0, 1, 0, 0,-1, 0, 1,-1]
    # Donde diff vale +1 empieza un tramo (índices 1 y 6) y donde vale -1
    # termina (índices 4 y 7). Los tramos son (1, 4) y (6, 7). Los ceros
    # agregados aseguran detectar tramos que tocan el principio o el final.
    sat = muestras_saturadas(x).astype(int)
    d = np.diff(np.r_[0, sat, 0])
    inicios, fines = np.flatnonzero(d == 1), np.flatnonzero(d == -1)

    # La duración mínima se pasa de segundos a muestras: dur_min * fs
    largo_min = int(round(dur_min * fs))
    # Se quedan sólo los tramos suficientemente largos (fin - inicio = largo)
    return [(i0, i1) for i0, i1 in zip(inicios, fines) if i1 - i0 >= largo_min]


def mascara_valida(ecg, fs):
    """True en las muestras aptas para el análisis.

    Excluye el primer DESCARTE_INICIO s y los tramos saturados largos con un
    margen de MARGEN_SAT s a cada lado. Se calcula sobre la señal CRUDA.

    No se borran muestras de la señal (eso cambiaría su longitud y correría
    los tiempos): se devuelve aparte un array de True/False del mismo largo.
    M2/M3 lo usan para descartar los latidos que caen donde vale False.
    """
    x = np.asarray(ecg, dtype=float)
    # Se arranca con todo válido...
    valida = np.ones(len(x), dtype=bool)
    # ...se marca como inválido el primer segundo (pasado a muestras)...
    valida[: int(round(DESCARTE_INICIO * fs))] = False
    # ...y cada tramo saturado largo, ampliado MARGEN_SAT a cada lado.
    m = int(round(MARGEN_SAT * fs))
    for i0, i1 in tramos_saturados(x, fs):
        # max(..., 0) evita índices negativos si el tramo está al principio.
        # Si i1 + m se pasa del final, Python recorta el slice solo.
        valida[max(i0 - m, 0): i1 + m] = False
    return valida


def espectro(x, fs, segundos_ventana=8):
    """DEP por Welch, para mostrar los ruidos en el cuaderno.

    Welch parte la señal en ventanas, calcula el espectro de cada una y los
    promedia: queda menos ruidoso que una sola FFT de toda la señal.
    Ventanas de 8 s dan una resolución de 1/8 = 0,125 Hz, suficiente para
    ver la deriva por debajo de 0,5 Hz.
    """
    x = np.asarray(x, dtype=float)
    # Largo de ventana en muestras (o la señal entera, si es más corta)
    nperseg = min(len(x), int(segundos_ventana * fs))
    # Se resta la media para que la continua no tape el resto del espectro
    return signal.welch(x - np.mean(x), fs=fs, nperseg=nperseg)


def resumen_calidad(ecg, fs):
    """Números para el informe: saturación, tramos largos y fracción válida."""
    x = np.asarray(ecg, dtype=float)
    tramos = tramos_saturados(x, fs)
    return {
        # % de muestras en los extremos (incluye puntas de R recortadas)
        "saturadas_%": 100 * muestras_saturadas(x).mean(),
        # cantidad de tramos saturados largos (los que se interpolan)
        "tramos_largos": len(tramos),
        # cuántos segundos suman esos tramos
        "s_en_tramos_largos": sum(i1 - i0 for i0, i1 in tramos) / fs,
        # % de la señal que queda para el análisis según mascara_valida
        "valida_%": 100 * mascara_valida(x, fs).mean(),
    }
