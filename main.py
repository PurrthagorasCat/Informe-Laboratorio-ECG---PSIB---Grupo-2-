"""Corre todo el trabajo práctico de punta a punta.

Uso, desde la carpeta del repo:

    python main.py

1. Revisa el estado de cada módulo (LISTO / FALTA / ERROR).
2. Si están todos listos, procesa los cuatro registros.
3. Genera las figuras y tablas finales en informe/figuras/.

No hace falta modificar este archivo para programar los módulos.
"""

from modulos.herramientas.procesamiento_completo import procesar_todos
from modulos.herramientas.verificacion import revisar_modulos
from modulos.m5_resultados_finales import resultados_finales


def main():
    print("\n1. Estado de los módulos\n")
    if not revisar_modulos():
        print("\nPara correr todo, los módulos M0 a M4 tienen que estar en LISTO.\n")
        return

    print("\n2. Procesamiento de los registros\n")
    resultados = procesar_todos()
    for res in resultados:
        reg = res.registro
        print(f"  {reg.nombre:20s} {reg.duracion:6.1f} s   {len(res.picos_r):4d} latidos")

    print("\n3. Resultados finales\n")
    try:
        resultados_finales.generar(resultados)
        print("  Figuras y tablas guardadas en informe/figuras/\n")
    except NotImplementedError as e:
        print(f"  FALTA   {e}\n")


if __name__ == "__main__":
    main()
