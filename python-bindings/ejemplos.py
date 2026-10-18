#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ejemplos.py — Recorrido didáctico por la librería alcp (cuerpos finitos).

Diez secciones cortas, de los cuerpos primos F_p a la FFT, pensadas para
leerse en la terminal (o aquí, en el código fuente) de arriba abajo.

Uso::

    python python-bindings/ejemplos.py               recorrido completo
    python python-bindings/ejemplos.py --lista       lista de secciones
    python python-bindings/ejemplos.py --seccion 7   solo la sección 7

Opciones:

    --sin-color   sin colores ANSI (también si existe la variable NO_COLOR o si
                  la salida no es una terminal)
    --color       fuerza los colores (p. ej. para verlos con «less -R»)
    --pausa       espera a que pulses Intro entre secciones
    --sin-pausa   sin pausas (es lo predeterminado)
    --rapido      tamaños reducidos en las medidas de tiempo (para pruebas)

Requiere el módulo nativo compilado: ``python build.py module``.

Cómo está hecho: cada ejemplo se escribe UNA sola vez, como texto, y se
ejecuta con eval/exec.  Así, el código que se ve en pantalla es exactamente
el que se ejecuta.  Solo se evalúan cadenas literales de este archivo, nunca
datos externos.
"""
from __future__ import annotations

import argparse
import contextlib
import errno
import io
import itertools
import os
import re
import shutil
import sys
import textwrap
import time
from fractions import Fraction

# El paquete alcp vive junto a este archivo (python-bindings/alcp).
_AQUI = os.path.dirname(os.path.abspath(__file__))
if _AQUI not in sys.path:
    sys.path.insert(0, _AQUI)

try:
    import alcp
except ImportError as exc:      # normalmente, el módulo nativo no está compilado
    raise SystemExit(f"No se puede importar alcp: {exc}\n"
                     "Compila antes el módulo nativo con:  python build.py module") from None

#: Semilla del generador del núcleo: cada sección empieza con ella, así que la
#: salida es reproducible (salvo los tiempos) y no depende de qué se ejecute antes.
SEMILLA = 2026

#: Límites del ancho de la salida (se adapta a la terminal entre estos valores).
ANCHO_MINIMO, ANCHO_MAXIMO = 64, 100

#: Grados de los polinomios de la sección de multiplicación rápida.
GRADOS_MULTIPLICACION = (256, 512, 1024, 2048)
GRADOS_MULTIPLICACION_RAPIDO = (16, 32, 64)


# ═══════════════════════════════════════════════════════════════════════════
# Salida en la terminal: colores, recuadros, tablas y barras
# ═══════════════════════════════════════════════════════════════════════════
class Salida:
    """Configuración global de la salida (se fija en :func:`main`)."""

    color = False           # ¿se emiten secuencias ANSI?
    ancho = 88              # columnas disponibles
    rapido = False          # ¿tamaños reducidos en las medidas?
    comprobaciones = 0      # comprobaciones realizadas...
    fallos = 0              # ... y cuántas han fallado


class Paleta:
    """Códigos SGR de ANSI (256 colores).

    Tonos medios, legibles con fondo oscuro y claro: azul para la estructura,
    teal para los resultados y naranja para los acentos.  El texto normal y el
    código usan el color por defecto de la terminal.
    """

    MARCO = "38;5;67"           # recuadros, reglas y separadores
    TITULO = "1"                # títulos de sección: negrita
    NUMERO = "1;38;5;173"       # número de sección
    APARTADO = "1;38;5;173"     # subtítulos «▸ ...»
    ENFASIS = "1"               # *palabras destacadas* en las explicaciones
    COMENTARIO_TEXTO = "3"      # comentarios entre ejemplos: cursiva
    NOTA = "38;5;244"           # notas al margen: gris
    PROMPT = "38;5;67"          # «>>> »
    CLAVE = "38;5;133"          # palabras clave de Python
    CADENA = "38;5;107"         # cadenas en el código
    NUMERO_CODIGO = "38;5;173"  # números en el código
    COMENTARIO = "3;38;5;244"   # comentarios en el código
    RESULTADO = "38;5;30"       # resultados
    CABECERA = "38;5;67"        # cabeceras de tabla
    OK = "38;5;71"              # ✔
    ERROR = "38;5;167"          # ✘
    TENUE = "38;5;244"          # ceros de las tablas
    UNO = "1;38;5;71"           # unos de las tablas de multiplicar
    LENTO = "38;5;167"          # barra del algoritmo escolar
    MEDIO = "38;5;179"          # barra de Karatsuba
    VELOZ = "38;5;71"           # barra de la FFT


def pinta(texto: str, sgr: str | None) -> str:
    """``texto`` con el estilo ANSI ``sgr`` (si el color está activado)."""
    if not (Salida.color and sgr and texto):
        return texto
    return f"\x1b[{sgr}m{texto}\x1b[0m"


def decimal(x: float, cifras: int = 1) -> str:
    """Número con coma decimal, como se escribe en español: 3,7."""
    return f"{x:.{cifras}f}".replace(".", ",")


# --- recuadros y textos ------------------------------------------------------
def portada() -> None:
    """Recuadro de bienvenida y leyenda de la salida."""
    ancho = Salida.ancho
    lineas = ["",
              f"alcp {alcp.__version__} · cuerpos finitos en Python",
              f"Un recorrido didáctico en {len(SECCIONES)} pasos",
              ""]
    print(pinta("╭" + "─" * (ancho - 2) + "╮", Paleta.MARCO))
    for i, texto in enumerate(lineas):
        estilo = Paleta.TITULO if i == 1 else Paleta.NOTA
        print(pinta("│", Paleta.MARCO) + pinta(f"   {texto}".ljust(ancho - 2), estilo)
              + pinta("│", Paleta.MARCO))
    print(pinta("╰" + "─" * (ancho - 2) + "╯", Paleta.MARCO))
    print()
    print("  Cómo leer la salida:")
    print("    " + pinta(">>> ", Paleta.PROMPT) + "código    "
          + pinta("se ejecuta tal cual (está en ejemplos.py)", Paleta.NOTA))
    print("    " + pinta("resultado", Paleta.RESULTADO) + "     "
          + pinta("lo que devuelve, como lo escribiría print()", Paleta.NOTA))
    print("    " + pinta("✔", Paleta.OK) + "  " + pinta("✘", Paleta.ERROR) + "          "
          + pinta("comprobación superada  ·  error esperado", Paleta.NOTA))


def recuadro_seccion(n: int, titulo: str) -> None:
    """Título de sección:  ╭────┬──────╮ / │  n │ título │ / ╰────┴──────╯."""
    numero = f" {n:>2} "
    interior = Salida.ancho - len(numero) - 3
    texto = f" {titulo}"[:interior].ljust(interior)
    print()
    print(pinta("╭" + "─" * len(numero) + "┬" + "─" * interior + "╮", Paleta.MARCO))
    print(pinta("│", Paleta.MARCO) + pinta(numero, Paleta.NUMERO) + pinta("│", Paleta.MARCO)
          + pinta(texto, Paleta.TITULO) + pinta("│", Paleta.MARCO))
    print(pinta("╰" + "─" * len(numero) + "┴" + "─" * interior + "╯", Paleta.MARCO))


_NBSP = chr(0xA0)              # espacio no separable (U+00A0)
_OPERADORES_TEXTO = "=+−-·<>≡≠≈×→↦Σ"


def ajusta(texto: str, ancho: int) -> list:
    """Como ``textwrap.wrap``, pero sin partir las fórmulas entre dos líneas.

    Los espacios dentro de paréntesis o llaves y los que rodean a un operador
    («a + 1», «p − 1», «F_p = {…}») se cambian por espacios no separables, en
    los que textwrap no corta, y se restauran después.
    """
    texto = " ".join(texto.split())
    protegido, profundidad = [], 0
    for i, c in enumerate(texto):
        if c in "([{":
            profundidad += 1
        elif c in ")]}":
            profundidad = max(0, profundidad - 1)
        elif c == " " and (profundidad > 0
                           or (i > 0 and texto[i - 1] in _OPERADORES_TEXTO)
                           or (i + 1 < len(texto) and texto[i + 1] in _OPERADORES_TEXTO)):
            c = _NBSP
        protegido.append(c)
    lineas = textwrap.wrap("".join(protegido), ancho)
    return [linea.replace(_NBSP, " ") for linea in lineas]


def _parrafo(texto: str, sangria: str, estilo: str | None) -> None:
    """Imprime ``texto`` ajustado al ancho; ``*así*`` se destaca en negrita."""
    lineas = ajusta(texto, Salida.ancho - len(sangria) - 1)
    bloque = "\n".join(sangria + linea for linea in lineas)
    partes = re.split(r"\*(.+?)\*", bloque, flags=re.DOTALL)
    # Las posiciones impares de «partes» son los trozos destacados.
    print("".join(pinta(p, Paleta.ENFASIS if i % 2 else estilo) for i, p in enumerate(partes)))


def explica(texto: str) -> None:
    """Explicación con la que empieza cada sección (2-4 líneas)."""
    _parrafo(texto, "  ", None)


def comenta(texto: str) -> None:
    """Comentario breve entre dos ejemplos (en cursiva)."""
    _parrafo(texto, "  ", Paleta.COMENTARIO_TEXTO)


def nota(texto: str) -> None:
    """Nota al margen (en gris)."""
    _parrafo(texto, "  ", Paleta.NOTA)


def apartado(texto: str) -> None:
    """Subtítulo dentro de una sección."""
    print()
    print("  " + pinta("▸ " + texto, Paleta.APARTADO))


def comprueba(descripcion: str, condicion: bool) -> bool:
    """Muestra ✔ o ✘ según ``condicion`` y lleva la cuenta de los fallos."""
    Salida.comprobaciones += 1
    if not condicion:
        Salida.fallos += 1
    lineas = ajusta(descripcion, Salida.ancho - 5) or [""]
    marca = pinta("✔", Paleta.OK) if condicion else pinta("✘", Paleta.ERROR)
    estilo = None if condicion else Paleta.ERROR
    print("  " + marca + " " + pinta(lineas[0], estilo))
    for linea in lineas[1:]:
        print("    " + pinta(linea, estilo))
    return condicion


# --- código y resultados -----------------------------------------------------
# Resaltado mínimo de Python: comentarios, cadenas, números y palabras clave.
_TOKENS = re.compile(r"""
      (?P<comentario>\#.*$)
    | (?P<cadena>[rbfu]{0,2}(?:"[^"\\]*(?:\\.[^"\\]*)*"|'[^'\\]*(?:\\.[^'\\]*)*'))
    | (?P<numero>\b\d+\b)
    | (?P<clave>\b(?:and|as|def|else|for|if|import|in|is|lambda|not|or|return|with
                     |True|False|None)\b)
""", re.VERBOSE)

_ESTILO_TOKEN = {"comentario": Paleta.COMENTARIO, "cadena": Paleta.CADENA,
                 "numero": Paleta.NUMERO_CODIGO, "clave": Paleta.CLAVE}


def resalta(linea: str) -> str:
    """Una línea de código Python con colores."""
    if not Salida.color:
        return linea
    trozos, pos = [], 0
    for m in _TOKENS.finditer(linea):
        trozos.append(linea[pos:m.start()])
        trozos.append(pinta(m.group(), _ESTILO_TOKEN[m.lastgroup]))
        pos = m.end()
    trozos.append(linea[pos:])
    return "".join(trozos)


def _separa_comentario(linea: str) -> tuple:
    """``(código, comentario)`` de una línea; el comentario vacío si no lo hay."""
    for m in _TOKENS.finditer(linea):
        if m.lastgroup == "comentario":
            return linea[:m.start()], m.group()
    return linea, ""


def bonito(valor) -> str:
    """Texto de un resultado tal como lo escribiría print().

    A diferencia de ``str(lista)``, que usa el repr de cada elemento
    (``[GF(7)(1), GF(7)(2)]``), aquí los contenedores muestran el str de sus
    elementos: ``[1, 2]``.
    """
    if isinstance(valor, str):
        return valor
    if isinstance(valor, list):
        return "[" + ", ".join(map(bonito, valor)) + "]"
    if isinstance(valor, tuple):
        coma = "," if len(valor) == 1 else ""
        return "(" + ", ".join(map(bonito, valor)) + coma + ")"
    if isinstance(valor, dict):
        return "{" + ", ".join(f"{bonito(k)}: {bonito(v)}" for k, v in valor.items()) + "}"
    return str(valor)


def _corta(linea: str, ancho: int, sangria: int = 4) -> list:
    """Divide una línea larga en trozos de ``ancho`` caracteres como máximo.

    Corta en un espacio y prefiere, por este orden, los que están fuera de
    paréntesis y los que siguen a «,» o «·»: así una factorización larga se
    parte entre dos factores y no dentro de uno.  Las líneas de continuación
    llevan ``sangria`` espacios más.
    """
    trozos = []
    while len(linea) > ancho:
        mejor, clave_mejor, profundidad = None, None, 0
        for i, caracter in enumerate(linea[:ancho + 1]):
            if caracter in "([{":
                profundidad += 1
            elif caracter in ")]}":
                profundidad -= 1
            elif caracter == " " and i > ancho // 2:        # nada de trozos muy cortos
                clave = (profundidad, linea[i - 1] not in ",·|")
                if clave_mejor is None or clave <= clave_mejor:   # el más a la derecha
                    mejor, clave_mejor = i, clave
        if mejor is None:
            break                       # no hay dónde cortar: se deja larga
        trozos.append(linea[:mejor].rstrip())
        linea = " " * sangria + linea[mejor + 1:]
    trozos.append(linea)
    return trozos


def _muestra_resultado(texto: str, estilo: str = Paleta.RESULTADO) -> None:
    """Imprime un resultado (de una o varias líneas) ajustado al ancho."""
    for linea in texto.splitlines():
        for trozo in _corta(linea, Salida.ancho - 3):
            print("  " + pinta(trozo, estilo))


def ejecuta(codigo: str, ns: dict, *, error: bool = False):
    """Muestra ``codigo`` con su «>>> », lo ejecuta en el espacio ``ns`` y
    muestra el resultado.

    Si el código es una expresión, se muestra su valor (con :func:`bonito`);
    si son sentencias, se muestra lo que impriman.  Con ``error=True`` el
    ejemplo DEBE lanzar una excepción, que se muestra como «✘ Tipo: mensaje»;
    sin él, una excepción es un fallo del recorrido y se propaga.
    """
    codigo = textwrap.dedent(codigo).strip("\n")
    for i, linea in enumerate(codigo.splitlines()):
        prompt = "... " if i and linea[:1].isspace() else ">>> "
        instruccion, comentario = _separa_comentario(linea)
        if comentario and len(linea) + 6 > Salida.ancho:
            # Terminal estrecha: el comentario va en su propia línea, encima.
            sangria = len(instruccion) - len(instruccion.lstrip())
            print("      " + " " * sangria + resalta(comentario))
            linea = instruccion.rstrip()
        print("  " + pinta(prompt, Paleta.PROMPT) + resalta(linea))

    capturado = io.StringIO()
    valor = None
    try:
        with contextlib.redirect_stdout(capturado):
            try:
                compilado = compile(codigo, "<ejemplo>", "eval")
            except SyntaxError:         # no es una expresión: sentencias
                exec(compile(codigo, "<ejemplo>", "exec"), ns)
            else:
                valor = eval(compilado, ns)
    except Exception as exc:  # noqa: BLE001  (se muestra o se propaga)
        if not error:
            raise
        _muestra_resultado(capturado.getvalue())
        tipo = type(exc).__name__
        _muestra_resultado(f"✘ {tipo}: {exc}", Paleta.ERROR)
        return exc
    if error:
        raise AssertionError(f"el ejemplo debía lanzar una excepción: {codigo!r}")
    _muestra_resultado(capturado.getvalue())
    if valor is not None:
        _muestra_resultado(bonito(valor))
    return valor


def espacio() -> dict:
    """Espacio de nombres nuevo para los ejemplos de una sección."""
    return {"alcp": alcp, "Fraction": Fraction, "itertools": itertools,
            "mobius": mobius, "divisores": divisores}


# --- tablas, matrices y barras ----------------------------------------------
def tabla(cabecera: list, filas: list, *, color=None, cortes=(1,),
          centrar: bool = False, sangria: int = 4) -> None:
    """Tabla alineada.

    ``cortes``: índices de las columnas precedidas por « │ »; ``color(i, j,
    texto)`` devuelve el estilo de la celda (fila i, columna j) o None.
    """
    anchos = [max(len(t) for t in columna) for columna in zip(cabecera, *filas)]

    def alinea(texto, ancho):
        return texto.center(ancho) if centrar else texto.rjust(ancho)

    def compone(celdas, estilo):
        trozos = []
        for j, (texto, ancho) in enumerate(zip(celdas, anchos)):
            if j:
                trozos.append(pinta(" │ ", Paleta.MARCO) if j in cortes else "  ")
            trozos.append(pinta(alinea(texto, ancho), estilo(j, texto)))
        return " " * sangria + "".join(trozos)

    regla = ""
    for j, ancho in enumerate(anchos):
        if j:
            regla += "─┼─" if j in cortes else "──"
        regla += "─" * ancho
    print(compone(cabecera, lambda j, t: Paleta.CABECERA))
    print(" " * sangria + pinta(regla, Paleta.MARCO))
    for i, fila in enumerate(filas):
        if color is None:
            print(compone(fila, lambda j, t: Paleta.CABECERA if j == 0 else None))
        else:
            print(compone(fila, lambda j, t, i=i: color(i, j, t)))


def tabla_operacion(simbolo: str, elementos: list, resultados: list) -> None:
    """Tabla de una operación del cuerpo: resalta los unos y atenúa los ceros."""
    cabecera = [simbolo] + [str(e) for e in elementos]
    filas = [[str(e)] + [str(v) for v in fila] for e, fila in zip(elementos, resultados)]

    def color(i, j, texto):
        if j == 0:
            return Paleta.CABECERA
        v = resultados[i][j - 1]
        return Paleta.UNO if v == 1 else Paleta.TENUE if v == 0 else None

    tabla(cabecera, filas, color=color, centrar=True)


def matriz(nombre: str, filas: list, sangria: int = 4) -> None:
    """Imprime una matriz entre corchetes ⎡ ⎤ con su nombre a la izquierda."""
    textos = [[str(v) for v in fila] for fila in filas]
    anchos = [max(len(t) for t in columna) for columna in zip(*textos)]
    n = len(textos)
    for i, fila in enumerate(textos):
        izquierda, derecha = (("⎡", "⎤") if i == 0 else ("⎣", "⎦") if i == n - 1
                              else ("⎢", "⎥"))
        etiqueta = f"{nombre} = " if i == (n - 1) // 2 else " " * (len(nombre) + 3)
        cuerpo = "  ".join(t.rjust(w) for t, w in zip(fila, anchos))
        print(" " * sangria + etiqueta + pinta(izquierda, Paleta.MARCO) + " "
              + pinta(cuerpo, Paleta.RESULTADO) + " " + pinta(derecha, Paleta.MARCO))


_OCTAVOS = ("", "▏", "▎", "▍", "▌", "▋", "▊", "▉")


def barra(valor: float, maximo: float, ancho: int) -> str:
    """Barra horizontal de ``ancho`` caracteres como máximo, con precisión de 1/8."""
    octavos = max(1, round(valor / maximo * ancho * 8)) if maximo > 0 else 1
    llenos, resto = divmod(octavos, 8)
    return "█" * llenos + _OCTAVOS[resto]


def mide(funcion, presupuesto: float = 0.25, maximo: int = 5):
    """``(mejor tiempo en segundos, resultado)`` de varias ejecuciones.

    Repite mientras el tiempo acumulado no pase de ``presupuesto`` (al menos
    una vez y como mucho ``maximo``): el mínimo es la medida menos ruidosa.
    """
    mejor, total, veces, resultado = float("inf"), 0.0, 0, None
    while veces < maximo and (veces == 0 or total < presupuesto):
        inicio = time.perf_counter()
        resultado = funcion()
        t = time.perf_counter() - inicio
        mejor, total, veces = min(mejor, t), total + t, veces + 1
    return mejor, resultado


# --- aritmética auxiliar (fórmula de Gauss, matrices) -----------------------
def divisores(n: int) -> list:
    """Divisores positivos de n."""
    return [d for d in range(1, n + 1) if n % d == 0]


def mobius(n: int) -> int:
    """Función de Möbius: 0 si n tiene un factor primo repetido y, si no,
    (−1)^k, con k el número de factores primos de n."""
    k, d = 0, 2
    while d * d <= n:
        if n % d == 0:
            n //= d
            if n % d == 0:
                return 0
            k += 1
        d += 1
    if n > 1:
        k += 1
    return -1 if k % 2 else 1


def toeplitz_completa(t: list, n: int) -> list:
    """Matriz de Toeplitz n×n: primera fila t[0:n] y primera columna t[0], t[n:]."""
    return [[t[j - i] if j >= i else t[n - 1 + i - j] for j in range(n)] for i in range(n)]


def producto_ingenuo(K, filas: list, v: list) -> list:
    """Producto matriz × vector «de libro» (para comparar con los algoritmos rápidos)."""
    return [sum((K(x) * K(y) for x, y in zip(fila, v)), K.zero) for fila in filas]


# ═══════════════════════════════════════════════════════════════════════════
# Secciones
# ═══════════════════════════════════════════════════════════════════════════
#: Lista de (título, función) en orden; el número de sección es la posición + 1.
SECCIONES: list = []


def seccion(titulo: str):
    """Decorador que registra una función como la siguiente sección."""
    def registra(funcion):
        SECCIONES.append((titulo, funcion))
        return funcion
    return registra


# ---------------------------------------------------------------------------
@seccion("Cuerpos primos: aritmética módulo p")
def cuerpos_primos() -> None:
    explica("""Un *cuerpo* es un conjunto en el que se puede sumar, restar,
        multiplicar y dividir (salvo por 0). El más sencillo es *F_p = {0, 1, …,
        p − 1}* con p primo: se opera como con enteros y se toma el resto de
        dividir entre p, como en un reloj de p horas.""")
    ns = espacio()

    apartado("Las operaciones en F₇")
    ejecuta("F = alcp.GF(7)", ns)
    ejecuta("F(3) + F(5)                 # 8 = 7 + 1", ns)
    ejecuta("F(3) * F(5)                 # 15 = 2·7 + 1", ns)
    ejecuta("F(3) - F(5)                 # −2 ≡ 5", ns)
    ejecuta("F(3) == 10                  # 3 y 10 dan el mismo resto", ns)

    apartado("Dividir es multiplicar por el inverso")
    comenta("Como 3·5 = 15 ≡ 1, el 5 es el inverso de 3: dividir entre 3 es "
            "multiplicar por 5.")
    ejecuta("F(3).inverse()", ns)
    ejecuta("F(2) / F(3)                 # 2·5 = 10 ≡ 3", ns)
    ejecuta("F(3) / 0", ns, error=True)

    apartado("Tabla de multiplicar de F₇")
    ejecuta("tabla = [[x * y for y in F] for x in F]", ns)
    tabla_operacion("·", list(ns["F"]), ns["tabla"])
    comenta("Cada fila, salvo la del 0, tiene exactamente un 1: todo elemento no "
            "nulo tiene inverso. Eso ocurre porque 7 es primo; con 6 no hay cuerpo "
            "(2·3 = 0):")
    ejecuta("alcp.GF(6)", ns, error=True)

    apartado("El pequeño teorema de Fermat")
    comenta("Si x ≠ 0, entonces x^(p−1) = 1 en F_p. Por tanto x^(p−2) es el inverso de x.")
    ejecuta("[x ** 6 for x in F if x != 0]", ns)
    comprueba("x⁶ = 1 para los seis elementos no nulos de F₇",
              all(x ** 6 == 1 for x in ns["F"] if x))
    ejecuta("F(3) ** 5, F(3) ** -1       # las dos formas de calcular 3⁻¹", ns)

    apartado("Utilidades sobre primos")
    ejecuta("alcp.is_prime(2**61 - 1)    # un primo de Mersenne", ns)
    ejecuta("alcp.next_prime(1000)", ns)


# ---------------------------------------------------------------------------
@seccion("Extensiones F_q: polinomios módulo un irreducible")
def extensiones() -> None:
    explica("""Si q = pⁿ con n > 1, los restos módulo q no forman un cuerpo. F_q se
        construye con polinomios: se elige uno *irreducible* g de grado n (que no
        se factoriza) y se calcula módulo g. Así, F₉ = F₃[a]/(a² + 1): sus
        elementos son c₀ + c₁·a y se usa la regla a² = −1.""")
    ns = espacio()

    apartado("F₉ y su generador a")
    ejecuta("K = alcp.GF(9)", ns)
    ejecuta("print(K)", ns)
    ejecuta("a = K.gen", ns)
    ejecuta("a ** 2                      # a² = −1 = 2", ns)
    ejecuta("list(K)", ns)
    ejecuta("(a + 1) * (a + 2)           # a² + 3a + 2 = a² + 2 = 1", ns)
    comenta("Así que a + 1 y a + 2 son inversos. Cada elemento es un vector de "
            "coordenadas sobre F₃ y tiene un índice c₀ + c₁·3:")
    ejecuta("(2*a + 1).coeffs, (2*a + 1).index", ns)

    apartado("Tabla de multiplicar de F₄ = F₂[a]/(a² + a + 1)")
    ejecuta("K4 = alcp.GF(4)", ns)
    ejecuta("tabla = [[x * y for y in K4] for x in K4]", ns)
    tabla_operacion("·", list(ns["K4"]), ns["tabla"])
    comenta("En F₄, a² = a + 1 (en característica 2, −1 = 1), así que "
            "a·(a + 1) = a² + a = 1.")

    apartado("El Frobenius x ↦ xᵖ")
    comenta("En característica p se cumple (x + y)ᵖ = xᵖ + yᵖ. Por eso el Frobenius "
            "x ↦ xᵖ respeta sumas y productos: permuta los elementos de F_q y deja "
            "fijos justo los de F_p.")
    ejecuta("[x.frobenius() for x in K]", ns)
    ejecuta("[x for x in K if x.frobenius() == x]", ns)
    K = ns["K"]
    comprueba("(x + y)³ = x³ + y³ para los 81 pares (x, y) de F₉",
              all((x + y) ** 3 == x ** 3 + y ** 3 for x in K for y in K))
    comenta("La traza x + xᵖ + … y la norma x·xᵖ·… no cambian al aplicar el "
            "Frobenius, así que caen en F_p:")
    ejecuta("(a + 1).trace(), (a + 1).norm()", ns)

    apartado("Orden multiplicativo y elemento primitivo")
    comenta("Los q − 1 elementos no nulos forman un grupo cíclico: hay un elemento "
            "primitivo g cuyas potencias los recorren todos.")
    ejecuta("a.multiplicative_order()    # a⁴ = 1: a no es primitivo", ns)
    ejecuta("g = K.primitive_element()", ns)
    ejecuta("g, g.multiplicative_order()", ns)
    ejecuta("[g ** k for k in range(8)]", ns)
    g = ns["g"]
    comprueba("las 8 potencias de g son los 8 elementos no nulos de F₉",
              {g ** k for k in range(8)} == {x for x in K if x})


# ---------------------------------------------------------------------------
@seccion("Conversión de tipos generosa")
def conversion() -> None:
    explica("""Para escribir con comodidad, alcp convierte lo que le des: enteros,
        cadenas como "2a + 1" (con x², multiplicación implícita y ·), listas de
        coeficientes, diccionarios, fracciones… Y si algo no tiene sentido, lo
        dice con un mensaje claro en lugar de dar un resultado erróneo.""")
    ns = espacio()
    ejecuta("K = alcp.GF(9); a = K.gen", ns)

    apartado("En un cuerpo: K(…)")
    ejecuta("K(5)                        # 5 ≡ 2 (mód 3)", ns)
    ejecuta('K("2a + 1")', ns)
    ejecuta('K("(a + 1)²")               # superíndices Unicode', ns)
    ejecuta("K([1, 2])                   # coordenadas c₀, c₁: 1 + 2a", ns)
    ejecuta("K(Fraction(1, 2))           # 1/2 = 2⁻¹ = 2 en F₃", ns)
    ejecuta("alcp.GF(7)(Fraction(3, 4))  # 3·4⁻¹ = 3·2 = 6", ns)

    apartado("En las operaciones")
    ejecuta('2 + a, "a" * a, Fraction(1, 2) * a', ns)
    ejecuta("alcp.GF(3)(2) * a           # F₃ se sumerge en F₉", ns)

    apartado("En un anillo de polinomios: R(…)")
    ejecuta('R = alcp.GF(7)["x"]', ns)
    ejecuta('R("3x² + 2x + 1")', ns)
    ejecuta("R([1, 2, 3])                # orden creciente: 1 + 2x + 3x²", ns)
    ejecuta("R({0: 1, 10: 3})            # disperso: {exponente: coeficiente}", ns)
    ejecuta('R("(x + 1)^3") * "x − 1"    # una cadena como operando', ns)

    apartado("Errores con mensajes claros")
    ejecuta("K(0.5)", ns, error=True)
    ejecuta("alcp.GF(10)", ns, error=True)
    ejecuta('K("2a +* 1")', ns, error=True)
    ejecuta("alcp.GF(5)(1) + alcp.GF(7)(1)", ns, error=True)
    ejecuta('alcp.GF(9, modulus="x^2 + x + 1")', ns, error=True)


# ---------------------------------------------------------------------------
@seccion("Estilos de impresión: Unicode, ASCII y LaTeX")
def impresion() -> None:
    explica("""Un mismo objeto se puede escribir de varias formas: en Unicode para
        leerlo en pantalla, en ASCII para copiarlo y pegarlo (alcp lo vuelve a
        leer) y en LaTeX para documentos y cuadernos Jupyter. Además, repr() da
        una forma inequívoca que reconstruye el objeto.""")
    ns = espacio()
    ejecuta('K = alcp.GF(9); R = K["x"]', ns)
    ejecuta('f = R("(a + 1)x^3 + 2x + a")', ns)

    apartado("Tres estilos")
    ejecuta("print(f)                    # Unicode, el estilo por defecto", ns)
    ejecuta('format(f, "ascii")', ns)
    ejecuta('format(f, "latex")', ns)
    ejecuta('format(K, "latex")', ns)
    ejecuta('R(format(f, "ascii")) == f  # ida y vuelta', ns)

    apartado("repr: inequívoco y evaluable")
    ejecuta("repr(f)", ns)
    ejecuta('eval(repr(f), {"GF": alcp.GF}) == f', ns)

    apartado("f-strings y estilo global")
    ejecuta('f"{f:u}   |   {f:a}"', ns)
    ejecuta('''\
        with alcp.print_style("ascii"):
            print(f.factor())''', ns)
    ejecuta("print(f.factor())           # al salir del with vuelve el Unicode", ns)
    ejecuta('f"{alcp.GF(7)(5):03d}"      # en F_p valen los formatos numéricos', ns)
    nota("alcp.set_print_style(...) cambia el estilo de forma permanente. En Jupyter, "
         "los objetos se muestran solos en LaTeX (_repr_latex_).")


# ---------------------------------------------------------------------------
@seccion("Polinomios: división, MCD y Bézout")
def polinomios() -> None:
    explica("""Los polinomios con coeficientes en un cuerpo se parecen a los
        enteros: se suman, se multiplican y se dividen con resto, f = q·g + r con
        grado(r) < grado(g). Por eso tienen *máximo común divisor*, que se calcula
        con el algoritmo de Euclides, y cumplen la *identidad de Bézout*.""")
    ns = espacio()
    ejecuta('R = alcp.GF(7)["x"]; x = R.gen', ns)

    apartado("Aritmética")
    ejecuta("f = x**4 + 3*x**3 + 2*x + 5", ns)
    ejecuta('g = R("x^2 + 1")', ns)
    ejecuta("f + g", ns)
    ejecuta("f * g", ns)

    apartado("División con resto")
    ejecuta("q, r = divmod(f, g)", ns)
    ejecuta("q, r", ns)
    f, g, q, r = (ns[k] for k in "fgqr")
    comprueba("f = q·g + r  y  grado(r) < grado(g)", f == q * g + r and r.degree < g.degree)
    ejecuta("f / g                       # «/» solo si la división es exacta", ns, error=True)

    apartado("Máximo común divisor e identidad de Bézout")
    ejecuta('f = R("(x + 1)(x + 2)(x^2 + 1)")', ns)
    ejecuta('g = R("(x + 1)(x + 4)")', ns)
    ejecuta("f, g", ns)
    ejecuta("d, s, t = f.xgcd(g)", ns)
    ejecuta("d, s, t", ns)
    ejecuta("s*f + t*g", ns)
    f, g, d, s, t = (ns[k] for k in "fgdst")
    comprueba("s·f + t·g = mcd(f, g) = x + 1", s * f + t * g == d == ns["R"]("x + 1"))

    apartado("Evaluación y composición")
    ejecuta("f(3)", ns)
    ejecuta("[f(c) for c in R.field]     # valores en los 7 elementos", ns)
    ejecuta("[c for c in R.field if f(c) == 0]   # raíces: −2 = 5 y −1 = 6", ns)
    ejecuta('f("x + 1")                  # composición: f(x + 1)', ns)
    ejecuta("f.derivative()", ns)


# ---------------------------------------------------------------------------
@seccion("Potencia modular e irreducibilidad")
def irreducibles() -> None:
    explica("""Un polinomio es *irreducible* si no es producto de otros de grado
        menor: son los «números primos» de los polinomios. Para saberlo no se
        prueban divisores: se calculan potencias x^(pᵏ) módulo f, que salen
        rápido elevando al cuadrado repetidamente (test de Ben-Or).""")
    ns = espacio()
    ejecuta('R = alcp.GF(2)["x"]; x = R.gen', ns)

    apartado("Potencias modulares")
    ejecuta('f = R("x^4 + x + 1")', ns)
    ejecuta("pow(x, 2**4, f)             # x¹⁶ ≡ x, como en todo cuerpo de 16 elementos", ns)
    ejecuta("pow(x, 10**30, f)           # un exponente de 31 cifras, al instante", ns)

    apartado("¿Irreducible?")
    ejecuta("f.is_irreducible()", ns)
    ejecuta('g = R("x^4 + x^2 + 1")', ns)
    ejecuta("g.is_irreducible(), g.factor()", ns)

    apartado("¿Cuántos hay? La fórmula de Gauss")
    comenta("Gauss contó los mónicos irreducibles de grado n sobre F_q: "
            "N_q(n) = (1/n)·Σ μ(d)·q^(n/d), sumando sobre los divisores d de n "
            "(μ es la función de Möbius). Es aproximadamente qⁿ/n. Comprobémoslo "
            "contándolos uno a uno:")
    ejecuta('''\
        def gauss(q, n):
            return sum(mobius(d) * q**(n // d)
                       for d in divisores(n)) // n''', ns)
    ejecuta('''\
        def cuenta(K, n):
            R = K["x"]
            return sum(R([*c, 1]).is_irreducible()
                       for c in itertools.product(K, repeat=n))''', ns)
    ejecuta("cuenta(alcp.GF(3), 2)       # x² + 1, x² + x + 2 y x² + 2x + 2", ns)
    nota("mobius y divisores son funciones auxiliares de este archivo.")

    gauss, cuenta = ns["gauss"], ns["cuenta"]
    cuerpos = (2, 3, 4)
    grados = range(1, 4 if Salida.rapido else 6)
    filas, correctos = [], 0
    for n in grados:
        fila = [str(n)]
        for q in cuerpos:
            teoria, real = gauss(q, n), cuenta(alcp.GF(q), n)
            correctos += teoria == real
            fila += [str(teoria), f"{real} {'✔' if teoria == real else '✘'}"]
        filas.append(fila)
    print()
    cabecera = ["n"] + [t for q in cuerpos for t in (f"Gauss q={q}", "real")]

    def color(i, j, texto):
        if j == 0:
            return Paleta.CABECERA
        if j % 2 == 0:
            return Paleta.OK if texto.endswith("✔") else Paleta.ERROR
        return Paleta.RESULTADO

    tabla(cabecera, filas, color=color, cortes=(1, 3, 5))
    casos = len(grados) * len(cuerpos)
    comprueba(f"la fórmula coincide con el recuento en los {casos} casos", correctos == casos)


# ---------------------------------------------------------------------------
@seccion("Factorización: el algoritmo de Cantor–Zassenhaus")
def factorizacion() -> None:
    explica("""Todo polinomio no nulo es, de forma única, una constante por un
        producto de irreducibles mónicos, igual que un entero es producto de
        primos. alcp lo factoriza con el algoritmo de *Cantor–Zassenhaus*: quita
        los factores repetidos, agrupa los factores por grado y separa cada grupo
        con un truco aleatorio.""")
    ns = espacio()
    ejecuta('R = alcp.GF(7)["x"]', ns)
    ejecuta('f = R("3 (x + 1)^2 (x + 5) (x^2 + 1)")', ns)
    ejecuta("f", ns)

    apartado("La factorización")
    ejecuta("fac = f.factor()", ns)
    ejecuta("print(fac)", ns)
    ejecuta("fac.unit, len(fac), fac[0]  # unidad, nº de factores y el primero", ns)
    comprueba("el producto de los factores reconstruye f", ns["fac"].expand() == ns["f"])
    ejecuta("f.roots()                   # pares (raíz, multiplicidad)", ns)

    apartado("Los tres pasos, uno a uno")
    ejecuta("f.is_squarefree()           # (x + 1)² está repetido", ns)
    ejecuta("u = f.radical()             # 1) sin repeticiones y mónico", ns)
    ejecuta("u", ns)
    ejecuta("u.distinct_degree_factorization()   # 2) producto de los de cada grado", ns)
    ejecuta('R("x^2 + 6x + 5").equal_degree_factorization(1)   # 3) se separan', ns)

    apartado("Un clásico: x^q − x")
    comenta("x^(pⁿ) − x es el producto de todos los mónicos irreducibles de grado "
            "divisor de n. Sobre F₃, x⁹ − x reúne los 3 lineales y los 3 cuadráticos "
            "irreducibles; sobre F₉ se parte en 9 factores lineales (uno por elemento).")
    ejecuta('print(alcp.GF(3)["x"]("x^9 - x").factor())', ns)
    ejecuta('print(alcp.GF(9)["x"]("x^9 - x").factor())', ns)

    apartado("Un polinomio aleatorio sobre F₉")
    ejecuta('h = alcp.GF(9)["x"].random(6)', ns)
    ejecuta("h", ns)
    ejecuta("print(h.factor())", ns)
    comprueba("también aquí el producto de los factores es h",
              ns["h"].factor().expand() == ns["h"])


# ---------------------------------------------------------------------------
@seccion("Multiplicación rápida: escolar, Karatsuba y FFT")
def multiplicacion() -> None:
    explica("""Multiplicar dos polinomios de grado n «como en el colegio» cuesta
        unas n² operaciones. *Karatsuba* cambia 4 productos por 3 al partir por la
        mitad (≈ n^1,58) y la *FFT* evalúa en raíces de la unidad, multiplica
        valores y vuelve a interpolar (≈ n·log n). Lo medimos en F_p con
        p = 998244353 = 119·2²³ + 1, que tiene raíces de la unidad de orden 2²³.""")
    ns = espacio()
    ejecuta('F = alcp.GF(998244353); R = F["x"]', ns)
    ejecuta("f, g = R.random(300), R.random(300)", ns)
    ejecuta('metodos = "school", "karatsuba", "fft"', ns)
    ejecuta("[f.mul(g, algorithm=m) == f * g for m in metodos]", ns)

    apartado("Tiempos (el mejor de varias repeticiones)")
    R = ns["R"]
    algoritmos = (("escolar", "school", Paleta.LENTO),
                  ("Karatsuba", "karatsuba", Paleta.MEDIO),
                  ("FFT", "fft", Paleta.VELOZ))
    # Calentamiento: la primera llamada de cada algoritmo prepara datos internos.
    f0, g0 = R.random(64), R.random(64)
    for _, alg, _ in algoritmos:
        f0.mul(g0, algorithm=alg)

    grados = GRADOS_MULTIPLICACION_RAPIDO if Salida.rapido else GRADOS_MULTIPLICACION
    tiempos, coinciden = {}, True
    for d in grados:
        f, g = R.random(d), R.random(d)
        medidas = {alg: mide(lambda: f.mul(g, algorithm=alg)) for _, alg, _ in algoritmos}
        tiempos[d] = {alg: t for alg, (t, _) in medidas.items()}
        productos = [producto for _, producto in medidas.values()]
        coinciden &= all(p == productos[0] for p in productos)

    ancho_barra = max(12, Salida.ancho - 42)
    nota("Barras relativas al más lento de cada grado; «×» = veces más rápido "
         "que el escolar.")
    for d in grados:
        t = tiempos[d]
        lento = max(t.values())
        print("    " + pinta(f"grado {d}", Paleta.CABECERA))
        for nombre, alg, estilo in algoritmos:
            dibujo = barra(t[alg], lento, ancho_barra).ljust(ancho_barra)
            mejora = "" if alg == "school" else f"×{decimal(t['school'] / t[alg])}"
            print(f"      {nombre:<10}" + pinta(dibujo, estilo)
                  + f" {decimal(t[alg] * 1000):>8} ms " + pinta(f"{mejora:>6}", Paleta.NOTA))
    comprueba("los tres algoritmos dan el mismo producto en todos los grados", coinciden)

    d1, d2 = grados[-2], grados[-1]
    crece = [decimal(tiempos[d2][alg] / tiempos[d1][alg]) for _, alg, _ in algoritmos]
    comenta(f"Al duplicar el grado (de {d1} a {d2}) el escolar tarda ×{crece[0]}, "
            f"Karatsuba ×{crece[1]} y la FFT ×{crece[2]}; la teoría predice ×4, ×3 y "
            "algo más de ×2.")
    nota("Los tiempos dependen del equipo y de su carga: lo importante es cómo crecen.")


# ---------------------------------------------------------------------------
@seccion("DFT y matrices de Toeplitz")
def dft_toeplitz() -> None:
    explica("""La *transformada discreta de Fourier* (DFT) de longitud n evalúa un
        polinomio en las n potencias de una raíz n-ésima primitiva de la unidad w;
        la inversa interpola y recupera los coeficientes. En F₁₇ hay raíces de
        orden 8 porque 8 divide a 17 − 1 = 16.""")
    ns = espacio()
    ejecuta('F = alcp.GF(17); R = F["x"]', ns)
    ejecuta("w = F.root_of_unity(8)", ns)
    ejecuta("w, w.multiplicative_order()", ns)

    apartado("La DFT es una evaluación")
    ejecuta("a = [1, 2, 3, 4, 0, 0, 0, 0]       # f = 1 + 2x + 3x² + 4x³", ns)
    ejecuta("A = alcp.dft(F, a)", ns)
    F, R, w, A = ns["F"], ns["R"], ns["w"], ns["A"]
    f = R(ns["a"])
    print()
    filas = [[str(i), str(w ** i), str(f(w ** i)), str(A[i])] for i in range(8)]
    tabla(["i", "wⁱ", "f(wⁱ)", "A[i]"], filas,
          color=lambda i, j, t: Paleta.CABECERA if j == 0 else Paleta.RESULTADO if j > 1
          else None)
    comprueba("A[i] = f(wⁱ) para i = 0, …, 7", all(f(w ** i) == A[i] for i in range(8)))
    ejecuta("alcp.idft(F, A)                    # la inversa recupera a", ns)

    apartado("Multiplicar = multiplicar valores punto a punto")
    ejecuta("b = [5, 6, 7, 0, 0, 0, 0, 0]       # g = 5 + 6x + 7x²", ns)
    ejecuta("B = alcp.dft(F, b)", ns)
    ejecuta("alcp.idft(F, [x * y for x, y in zip(A, B)])", ns)
    ejecuta("(R(a) * R(b)).coeffs", ns)
    producto = alcp.idft(F, [x * y for x, y in zip(A, ns["B"])])
    comprueba("idft(dft(a)·dft(b)) son los coeficientes de f·g (completados con ceros)",
              producto == (R(ns["a"]) * R(ns["b"])).coeffs + [F.zero] * 2)

    apartado("Matrices de Toeplitz")
    comenta("Una matriz de Toeplitz es constante en cada diagonal, así que la "
            "determinan su primera fila y su primera columna: 2n − 1 números en lugar "
            "de n². Multiplicar por ella es, en el fondo, multiplicar polinomios.")
    ejecuta("t = [1, 2, 3, 4, 5, 6, 7]   # fila 1 2 3 4; columna sin la esquina: 5 6 7",
            ns)
    T = toeplitz_completa(ns["t"], 4)
    print()
    matriz("T", T)
    print()
    ejecuta("v = [1, 0, 2, 1]", ns)
    ejecuta("alcp.toeplitz_vec(F, t, v)  # T·v sin construir la matriz", ns)
    comprueba("coincide con el producto matriz × vector de libro",
              alcp.toeplitz_vec(F, ns["t"], ns["v"]) == producto_ingenuo(F, T, ns["v"]))

    apartado("Inversa de una Toeplitz triangular y división de polinomios")
    ejecuta("c = [2, 1, 0, 3]            # primera columna de L (triangular inferior)", ns)
    ejecuta("inv = alcp.toeplitz_lower_inverse(F, c)", ns)
    ejecuta("alcp.toeplitz_lower_vec(F, c, inv)   # L·(1ª columna de L⁻¹) = e₁", ns)
    comprueba("L·L⁻¹ tiene por primera columna (1, 0, 0, 0)",
              alcp.toeplitz_lower_vec(F, ns["c"], ns["inv"]) == [1, 0, 0, 0])
    ejecuta('f, g = R("x^7 + 3x^5 + 2x + 1"), R("2x^3 + x + 5")', ns)
    ejecuta("alcp.divmod_toeplitz(f, g)  # cociente y resto con sistemas de Toeplitz", ns)
    comprueba("divmod_toeplitz(f, g) = divmod(f, g)",
              alcp.divmod_toeplitz(ns["f"], ns["g"]) == divmod(ns["f"], ns["g"]))


# ---------------------------------------------------------------------------
#: Resumen final: (idea, funciones de alcp con las que se trabaja).
_RESUMEN = (
    ("F_p son los enteros módulo p: todo x ≠ 0 tiene inverso y x^(p−1) = 1.",
     "alcp.GF(7)   x.inverse()   x ** k"),
    ("F_q = F_p[a]/(g): elementos como polinomios en a; Frobenius, traza, norma "
     "y elemento primitivo.",
     "alcp.GF(9)   K.gen   x.frobenius()   K.primitive_element()"),
    ("Conversión generosa (enteros, cadenas, listas, dicts, Fraction) y errores claros.",
     'K("2a + 1")   R([1, 2, 3])   R({0: 1, 10: 3})'),
    ("Tres estilos de impresión y un repr que reconstruye el objeto.",
     'format(f, "latex")   alcp.print_style("ascii")   repr(f)'),
    ("Polinomios: división con resto, MCD, Bézout, evaluación y composición.",
     "divmod(f, g)   f.xgcd(g)   f(3)   f(g)"),
    ("Irreducibilidad con potencias modulares; Gauss cuenta los irreducibles.",
     "pow(x, k, f)   f.is_irreducible()"),
    ("Factorización de Cantor–Zassenhaus y raíces con su multiplicidad.",
     "f.factor()   f.roots()   f.radical()"),
    ("Escolar, Karatsuba y FFT: el mismo producto con costes muy distintos.",
     'f.mul(g, algorithm="fft")'),
    ("DFT = evaluar en raíces de la unidad; Toeplitz = productos de polinomios.",
     "alcp.dft(F, a)   alcp.toeplitz_vec(F, t, v)"),
)


@seccion("Resumen: lo que hemos aprendido")
def resumen() -> None:
    explica("""Un repaso de las ideas del recorrido y de las funciones de alcp con
        las que se ponen en práctica:""")
    print()
    for n, (idea, funciones) in enumerate(_RESUMEN, start=1):
        lineas = ajusta(idea, Salida.ancho - 7)
        print("  " + pinta(f"{n:>2}", Paleta.NUMERO) + "  " + lineas[0])
        for linea in lineas[1:]:
            print("      " + linea)
        print("      " + resalta(funciones))
    print()
    comenta("Para seguir: help(alcp) y los docstrings de cada función, "
            "docs/diseno.md (el diseño de la librería) y python-bindings/tests/ "
            "(cientos de ejemplos más, en forma de pruebas).")


# ═══════════════════════════════════════════════════════════════════════════
# Programa principal
# ═══════════════════════════════════════════════════════════════════════════
def decide_color(sin_color: bool, forzar: bool, entorno, es_terminal: bool) -> bool:
    """¿Usar colores?  --sin-color gana a todo; después, --color; después, la
    variable NO_COLOR (si está definida, sin color); y, por último, solo si la
    salida es una terminal."""
    if sin_color:
        return False
    if forzar:
        return True
    if "NO_COLOR" in entorno:
        return False
    return es_terminal


def _activa_vt_windows() -> bool:
    """En la consola de Windows, activa el procesado de secuencias ANSI.

    Devuelve False solo si la consola no lo admite (Windows anterior a 10).
    Fuera de Windows, o si la salida no es una consola, no hace nada.
    """
    if os.name != "nt":
        return True
    try:
        import ctypes
        from ctypes import wintypes
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.GetStdHandle.restype = wintypes.HANDLE
        kernel32.GetStdHandle.argtypes = (wintypes.DWORD,)
        kernel32.GetConsoleMode.argtypes = (wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD))
        kernel32.SetConsoleMode.argtypes = (wintypes.HANDLE, wintypes.DWORD)
        manejador = kernel32.GetStdHandle(-11 & 0xFFFFFFFF)     # STD_OUTPUT_HANDLE = (DWORD)-11
        modo = wintypes.DWORD()
        if not kernel32.GetConsoleMode(manejador, ctypes.byref(modo)):
            return True         # no es una consola (archivo, tubería, mintty...)
        activar_vt = 0x0004     # ENABLE_VIRTUAL_TERMINAL_PROCESSING
        return bool(kernel32.SetConsoleMode(manejador, modo.value | activar_vt))
    except (OSError, AttributeError, ValueError):
        return True


def _salida_utf8() -> None:
    """Escribe en UTF-8 aunque la consola o la redirección usen otra codificación."""
    codificacion = (getattr(sys.stdout, "encoding", None) or "").lower().replace("-", "")
    if codificacion != "utf8":
        with contextlib.suppress(AttributeError, ValueError, io.UnsupportedOperation):
            sys.stdout.reconfigure(encoding="utf-8")


def _es_terminal() -> bool:
    try:
        return sys.stdout.isatty()
    except (AttributeError, ValueError):
        return False


def configura_salida(args) -> None:
    """Fija color, ancho y modo rápido según las opciones y el entorno."""
    _salida_utf8()
    terminal = _es_terminal()
    Salida.color = decide_color(args.sin_color, args.color, os.environ, terminal)
    if Salida.color and terminal and not _activa_vt_windows():
        Salida.color = False    # consola antigua de Windows sin secuencias ANSI
    columnas = shutil.get_terminal_size((96, 24)).columns
    Salida.ancho = max(ANCHO_MINIMO, min(ANCHO_MAXIMO, columnas - 1))
    Salida.rapido = args.rapido
    Salida.comprobaciones = Salida.fallos = 0


class _FormatoAyuda(argparse.HelpFormatter):
    """Ayuda de argparse con el prefijo «uso:» en lugar de «usage:»."""

    def add_usage(self, usage, actions, groups, prefix=None):
        super().add_usage(usage, actions, groups, "uso: " if prefix is None else prefix)


def analiza_argumentos(argv):
    """Opciones de la línea de órdenes (con la ayuda y los errores en español)."""
    parser = argparse.ArgumentParser(
        prog="ejemplos.py", add_help=False, formatter_class=_FormatoAyuda,
        description="Recorrido didáctico por la librería alcp (cuerpos finitos).")
    opciones = parser.add_argument_group("opciones")
    opciones.add_argument("-h", "--ayuda", action="help", help="muestra esta ayuda y termina")
    opciones.add_argument("--seccion", metavar="N",
                          help="muestra solo la sección N (ver --lista)")
    opciones.add_argument("--lista", action="store_true", help="lista las secciones")
    opciones.add_argument("--sin-color", action="store_true",
                          help="sin colores (también con NO_COLOR o si la salida no es "
                               "una terminal)")
    opciones.add_argument("--color", action="store_true",
                          help="fuerza los colores aunque la salida no sea una terminal")
    opciones.add_argument("--pausa", action="store_true",
                          help="espera a que se pulse Intro entre secciones")
    opciones.add_argument("--sin-pausa", action="store_true",
                          help="sin pausas (es el comportamiento por defecto)")
    opciones.add_argument("--rapido", action="store_true",
                          help="tamaños reducidos en las medidas de tiempo (para pruebas)")
    args = parser.parse_args(argv)
    if args.seccion is not None:
        try:
            args.seccion = int(args.seccion)
        except ValueError:
            parser.error(f"--seccion espera un número entero, no {args.seccion!r}")
        if not 1 <= args.seccion <= len(SECCIONES):
            parser.error(f"la sección debe estar entre 1 y {len(SECCIONES)}")
    if args.sin_pausa:
        args.pausa = False
    return args


def lista() -> None:
    """Lista numerada de las secciones."""
    print()
    print("  " + pinta("Secciones del recorrido", Paleta.TITULO)
          + pinta("  (python ejemplos.py --seccion N)", Paleta.NOTA))
    print()
    for n, (titulo, _) in enumerate(SECCIONES, start=1):
        print("   " + pinta(f"{n:>2}", Paleta.NUMERO) + "  " + titulo)
    print()


def ejecuta_seccion(n: int) -> None:
    titulo, funcion = SECCIONES[n - 1]
    alcp.set_seed(SEMILLA)
    recuadro_seccion(n, titulo)
    funcion()


def _pausa() -> None:
    """Espera a que se pulse Intro (solo si la entrada es una terminal)."""
    if sys.stdin is not None and sys.stdin.isatty():
        with contextlib.suppress(EOFError):
            input("\n  " + pinta("Pulsa Intro para continuar…", Paleta.NOTA))


def pie(segundos: float) -> None:
    """Resumen de las comprobaciones y del tiempo empleado."""
    total, fallos = Salida.comprobaciones, Salida.fallos
    print()
    print(pinta("─" * Salida.ancho, Paleta.MARCO))
    if fallos:
        mensaje = pinta(f"✘ {fallos} de {total} comprobaciones han fallado", Paleta.ERROR)
    elif total == 0:
        mensaje = pinta("✔", Paleta.OK) + " completado (sin comprobaciones automáticas)"
    elif total == 1:
        mensaje = pinta("✔", Paleta.OK) + " 1 comprobación superada"
    else:
        mensaje = pinta("✔", Paleta.OK) + f" {total} comprobaciones superadas"
    print("  " + mensaje + pinta(f"   ·   {decimal(segundos)} s", Paleta.NOTA))
    print()


def main(argv=None) -> int:
    """Punto de entrada.  Devuelve 0 si todas las comprobaciones pasan."""
    args = analiza_argumentos(argv)
    configura_salida(args)
    if args.lista:
        lista()
        return 0
    estilo_anterior = alcp.get_print_style()
    alcp.set_print_style("unicode")
    inicio = time.perf_counter()
    try:
        if args.seccion is not None:
            ejecuta_seccion(args.seccion)
        else:
            portada()
            for n in range(1, len(SECCIONES) + 1):
                ejecuta_seccion(n)
                if args.pausa and n < len(SECCIONES):
                    _pausa()
        pie(time.perf_counter() - inicio)
    finally:
        alcp.set_print_style(estilo_anterior)
    return 1 if Salida.fallos else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print()
        sys.exit(130)
    except OSError as exc:
        # Salida cortada (p. ej. «python ejemplos.py | head»): llega como
        # BrokenPipeError (EPIPE) en Linux y como OSError EINVAL en Windows.
        # Se sale sin traza; cualquier otro error se propaga.
        if exc.errno not in (errno.EPIPE, errno.EINVAL):
            raise
        with contextlib.suppress(OSError):
            os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        sys.exit(1)
