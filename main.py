"""Corre todo el trabajo práctico de punta a punta.

Uso, desde la carpeta del repo:

    python main.py

1. Revisa el estado de cada módulo (LISTO / FALTA / ERROR) con señales de prueba.
2. Si están todos listos, procesa los registros de datos/lista_de_registros.csv:
   M0 carga -> M1 acondicionamiento -> M2 Pan-Tompkins -> M3 tacograma -> M4 Lomb.
3. Muestra las tablas de resultados (ítems 6 y 7) y guarda las figuras y las
   tablas del informe en informe/figuras/ (M5).
"""

import sys

import matplotlib

matplotlib.use("Agg")   # sin ventanas: las figuras sólo se guardan en PDF

import pandas as pd  # noqa: E402

from modulos.herramientas.procesamiento_completo import procesar_todos  # noqa: E402
from modulos.herramientas.verificacion import revisar_modulos  # noqa: E402
from modulos.m5_resultados_finales import resultados_finales  # noqa: E402

# Columnas que se muestran en pantalla (los CSV guardan todas)
COLUMNAS_TACOGRAMA = ["N_RR", "FC_media_lpm", "FC_inicio_lpm", "FC_fin_lpm", "pendiente_FC_lpm_min",
                      "RR_medio_ms", "SDNN_ms", "RMSSD_ms", "pNN50_%"]
COLUMNAS_LOMB = ["n_rr", "T_s", "FC_lpm", "LF_ms2", "HF_ms2", "LF_HF", "LFnu_%", "pico_HF_Hz"]


def titulo(texto):
    print(f"\n{texto}\n")


def mostrar(tabla, decimales=1):
    with pd.option_context("display.width", 140, "display.max_columns", 20):
        print(tabla.round(decimales).to_string())


def main():
    # Algunas consolas de Windows no tienen "²" o "→": se reemplazan en vez de cortar el programa
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(errors="replace")

    titulo("1. Estado de los módulos")
    if not revisar_modulos():
        print("\nPara correr todo, los módulos M0 a M4 tienen que estar en LISTO.\n")
        return

    titulo("2. Procesamiento de los registros")
    resultados = procesar_todos()
    for res in resultados:
        reg = res.registro
        print(f"  {reg.nombre:20s} fs = {reg.fs:g} Hz  {reg.duracion:6.1f} s  "
              f"{len(res.picos_r):4d} latidos  {len(res.rr):4d} RR en zona válida")

    titulo("3. Resultados finales")
    salida = resultados_finales.generar(resultados)
    tablas = salida["tablas"]

    print("Detección de latidos")
    mostrar(tablas["tabla_deteccion"])
    print("\nTacograma (ítem 6)")
    mostrar(tablas["tabla_tacograma"][COLUMNAS_TACOGRAMA])
    print("\nPeriodograma de Lomb (ítem 7)")
    mostrar(tablas["tabla_lomb"][COLUMNAS_LOMB], decimales=3)

    print("\nFiguras guardadas en informe/figuras/:", ", ".join(f"{n}.pdf" for n in salida["figuras"]))
    print("Tablas guardadas en informe/figuras/:", ", ".join(f"{n}.csv" for n in tablas), "\n")


if __name__ == "__main__":
    main()
