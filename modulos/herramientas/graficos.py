"""Estilo común para que todas las figuras del informe se vean iguales.

Regla para TODAS las figuras:
  - el color indica la condición (reposo / post ejercicio)
  - el tipo de línea indica el voluntario (V1 continua / V2 discontinua)
Así nada depende sólo del color (sirve también si se imprime en blanco y negro).
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

CARPETA_FIGURAS = Path(__file__).resolve().parents[2] / "informe" / "figuras"   # informe/figuras/ del repo

COLOR_CONDICION = {
    "reposo": "#2a78d6",
    "post_ejercicio": "#eb6834",
}
LINEA_VOLUNTARIO = {
    "V1": "-",
    "V2": "--",
}
COLOR_NEUTRO = "#5f5e5a"  # para señales sin condición (p. ej. etapas intermedias)


def aplicar_estilo():
    """Llamarla una vez al principio de cada cuaderno."""
    cm = 1 / 2.54
    mpl.rcParams.update({
        "figure.figsize": (16 * cm, 7 * cm),  # 16 cm = ancho del texto en el informe
        "figure.dpi": 110,
        "savefig.bbox": "tight",
        "font.size": 9,
        "axes.titlesize": 9,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.5,
        "lines.linewidth": 1.0,
    })


def estilo_registro(voluntario, condicion):
    """Color, tipo de línea y etiqueta para plot(). Uso: ax.plot(t, x, **estilo_registro("V1", "reposo"))"""
    return {
        "color": COLOR_CONDICION.get(condicion, COLOR_NEUTRO),
        "linestyle": LINEA_VOLUNTARIO.get(voluntario, "-"),
        "label": f"{voluntario} - {condicion.replace('_', ' ')}",
    }


def guardar_figura(fig, nombre):
    """Guarda la figura en informe/figuras/<nombre>.pdf, lista para el informe."""
    CARPETA_FIGURAS.mkdir(parents=True, exist_ok=True)
    ruta = CARPETA_FIGURAS / f"{nombre}.pdf"
    fig.savefig(ruta)
    return ruta
