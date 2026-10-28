# Laboratorio ECG — Variabilidad de la frecuencia cardíaca

Procesamiento de Señales e Imágenes Biomédicas · ITBA · 2026Q1
**Entrega del informe: miércoles 04/11/2026, 23:59 h.**

## Qué hay en el repositorio

```
lab-ecg-variabilidad/
├── main.py                      corre todo el TP de punta a punta
├── README.md                    este archivo
├── requisitos.txt               librerías de Python que hay que instalar
├── .gitignore                   archivos que Git no sube (no hace falta tocarlo)
│
├── modulos/                     una carpeta por módulo: su código + su cuaderno
│   ├── m0_carga_de_datos/
│   │   ├── carga_de_datos.py          se programa acá
│   │   └── carga_de_datos.ipynb       se prueba, se grafica y se discute acá
│   ├── m1_acondicionamiento/
│   ├── m2_pan_tompkins/
│   ├── m3_tacograma/
│   ├── m4_periodograma_lomb/
│   ├── m5_resultados_finales/         figuras y tablas comparativas finales
│   └── herramientas/                  ya está hecho, no hace falta tocarlo
│       ├── procesamiento_completo.py     corre M1 → M4 en orden
│       ├── verificacion.py               revisa si cada módulo cumple lo acordado
│       ├── ecg_de_prueba.py              ECG inventado para probar sin datos reales
│       └── graficos.py                   estilo común de las figuras
│
├── datos/
│   ├── lista_de_registros.csv       qué archivo es de qué voluntario y condición
│   └── registros_crudos/            los archivos tal cual salen de la interfaz
│
└── informe/                     el informe en LaTeX
    ├── informe.tex                  archivo principal
    ├── secciones/                   una sección por archivo
    ├── figuras/                     las figuras que guardan los cuadernos
    └── referencias.bib              bibliografía
```

**La idea:** cada integrante trabaja **sólo dentro de la carpeta de su módulo**. El `.py` tiene las funciones; el cuaderno las prueba, grafica y discute.

## Para empezar

**Una sola persona** sube esta carpeta a GitHub: crear un repositorio **privado** vacío y, desde esta carpeta, correr:

```bash
git init -b main
git add .
git commit -m "Esqueleto del TP"
git remote add origin <link del repositorio>
git push -u origin main
```

Después suma a los demás como colaboradores del repositorio.

**Cada integrante**, una vez:

```bash
git clone <link del repositorio>
cd lab-ecg-variabilidad
pip install -r requisitos.txt
python main.py
```

## Cómo saber qué falta: `python main.py`

`main.py` primero revisa cada módulo y muestra un renglón por módulo:

```
1. Estado de los módulos

  M0 Carga de datos         FALTA   todavía no hay registros en datos/lista_de_registros.csv
  M1 Acondicionamiento      FALTA   M1: falta programar acondicionar()
  M2 Pan-Tompkins           LISTO
  ...
```

- `LISTO`: el módulo devuelve lo que los demás esperan.
- `FALTA`: hay funciones sin programar (o faltan los registros).
- `ERROR`: está programado pero no cumple; el mensaje dice qué está mal.

Cuando M0 a M4 están en `LISTO`, `main.py` procesa los cuatro registros y genera las figuras y tablas finales (M5) en `informe/figuras/`.

Antes de subir cambios, tu módulo tiene que estar en `LISTO`.

## División del trabajo

| Carpeta | Ítems | Responsable |
|---|---|---|
| `m0_carga_de_datos` | 1, 2 | |
| `m1_acondicionamiento` | 4, 5 | |
| `m2_pan_tompkins` | 6 (detección) | |
| `m3_tacograma` | 6a, 6b | |
| `m4_periodograma_lomb` | 7a, 7b | |
| `m5_resultados_finales` | figuras finales de 6 y 7 | |
| `informe/secciones/01_` y `02_` (teoría) | 1 (otros formatos), 3a, 3b, 3c | |
| `informe/secciones/04_discusion.tex` | 8 | todos |

Con 3 integrantes: **A** → M0 + M1, **B** → M2, **C** → M3 + M4. La teoría (ítems 1 y 3) y M5 los toma quien termine primero. Con 4 integrantes, M3 y M4 van por separado.

### Todos pueden empezar hoy

No hace falta esperar al módulo anterior. `herramientas/ecg_de_prueba.py` genera señales inventadas con los latidos en posiciones conocidas:

| Módulo | Mientras no está lo anterior, prueba con |
|---|---|
| M1 y M2 | `generar_ecg()`: un ECG con sus picos R |
| M3 | los picos R que devuelve `generar_ecg()` |
| M4 | `generar_latidos()`: un tacograma que oscila a una frecuencia conocida |

Lo único que sí hay que tener cuanto antes son **los registros reales y el módulo M0**.

## Lo que cada módulo recibe y devuelve

Esto es lo que permite trabajar por separado. Si alguien necesita cambiarlo, lo avisa al grupo primero.

| Función | Devuelve |
|---|---|
| `cargar_registro(ruta)` | un `Registro` con `ecg` (señal cruda) y `fs` (Hz) |
| `acondicionar(ecg, fs)` | la señal acondicionada, **misma cantidad de muestras y sin corrimiento** |
| `detectar_picos_r(ecg, fs)` | posiciones de los picos R (números de muestra), **referidas a la señal que entró** |
| `intervalos_rr(picos_r, fs)` | `t_rr` y `rr` en segundos; `t_rr[i]` es el latido que cierra el intervalo |
| `periodograma_lomb(t_rr, rr)` | `f` (Hz) y `p` |

Tiempos e intervalos siempre en **segundos**; milisegundos sólo en los ejes de los gráficos.

Desde un cuaderno, `procesar_todos(hasta="pan_tompkins")` corre los cuatro registros sólo hasta la etapa que ya esté programada.

## Cómo trabajar con Git

1. **No subir cambios directo a `main`.** Cada uno trabaja en su propia rama: `m0-carga`, `m1-acondicionamiento`, `m2-pan-tompkins`, `m3-tacograma`, `m4-lomb`, `m5-resultados`, `informe-teoria`.
   ```bash
   git checkout -b m2-pan-tompkins          # crear tu rama (una vez)
   git add modulos/m2_pan_tompkins
   git commit -m "M2: filtro pasabanda"
   git push -u origin m2-pan-tompkins
   ```
2. Cuando el módulo está en `LISTO`, abrir un **Pull Request** en GitHub hacia `main`. Otro integrante lo revisa y lo aprueba.
3. Antes del Pull Request, traer lo último de `main` a tu rama: `git pull origin main`.
4. Cada uno edita **sólo su carpeta de `modulos/` y su sección del informe**. Si necesitás un cambio en el módulo de otro, se lo pedís.
5. Las funciones van en el `.py`; el cuaderno sólo llama, grafica y comenta.
6. Mensajes de commit cortos, empezando por el módulo: `M3: frecuencia cardíaca instantánea`.

### Antes de abrir un Pull Request

- [ ] `python main.py` muestra tu módulo en `LISTO`
- [ ] Las funciones nuevas están en el `.py`, no sueltas en el cuaderno
- [ ] Las figuras del informe se guardan con `graficos.guardar_figura()`
- [ ] No modificaste archivos fuera de tu carpeta (o lo avisaste)

## Datos

1. Copiar cada archivo exportado a `datos/registros_crudos/` **sin modificarlo**.
2. Nombrarlo `V1_reposo`, `V1_post_ejercicio`, `V2_reposo` o `V2_post_ejercicio` (con la extensión original).
3. Completar su fila en `datos/lista_de_registros.csv` (archivo, fecha y observaciones: electrodo que se despegó, movimiento, etc.).

Los voluntarios van como **V1** y **V2**: no poner nombres en el repositorio.

## Figuras

Todas con el mismo estilo: **el color indica la condición** y **el tipo de línea, el voluntario**.

```python
ax.plot(t, x, **graficos.estilo_registro("V1", "reposo"))
graficos.guardar_figura(fig, "tacogramas")   # → informe/figuras/tacogramas.pdf
```

En el informe: `\includegraphics{tacogramas}`.

## Informe

- LaTeX, con una sección por archivo en `informe/secciones/` para no pisarse. En Overleaf: subir la carpeta `informe/` completa y compilar `informe.tex`.
- Lo que falta aparece **en rojo** con `\pendiente{...}`. Para la entrega, en `informe.tex` cambiar `\borradortrue` por `\borradorfalse`.
- Citas con `\parencite{clave}`; la bibliografía sale en formato APA, en orden de aparición. Las fuentes se agregan en `referencias.bib`.
- Cada sección tiene un comentario `% Responde: ítem X`.

| Ítem | Dónde se trabaja | Sección del informe |
|---|---|---|
| 1 Formato del registro / otros formatos | M0 | 2.2 / 1.2 |
| 2 Frecuencia de muestreo | M0 | 2.2 |
| 3a Resolución AD8232 | — | 2.2 |
| 3b Electrodos en miembros | — | 1.1 |
| 3c Otras derivaciones | — | 1.1 |
| 4 ¿Hace falta acondicionar? | M1 | 2.3 |
| 5 Acondicionamiento | M1 | 3.1 |
| 6 Pan-Tompkins y tacograma (6a, 6b) | M2, M3 | 2.4, 2.5, 3.2, 3.3 |
| 7 Periodograma de Lomb (7a, 7b) | M4 | 2.6, 3.4 |
| 8 Conclusiones | M5 | 4, 5 |

## Cronograma sugerido

| Fecha | Objetivo |
|---|---|
| jue 29/10 | Repositorio creado, responsables asignados, registros en `datos/` |
| vie 30/10 | M0 en `main` · borrador de teoría (ítems 1 y 3) |
| dom 01/11 | M1 y M2 en `main` |
| lun 02/11 | M3 y M4 en `main` |
| mar 03/11 | M5 en `main` · informe completo · revisión cruzada · discusión (ítem 8) |
| mié 04/11 | Revisión final, `\borradorfalse` y entrega antes de las 23:59 |
