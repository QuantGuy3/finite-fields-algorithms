"""Representación textual de cuerpos, elementos, polinomios y factorizaciones.

Hay tres estilos:

* ``"unicode"`` (por defecto): ``x³ + 2x + 1``, ``(a + 1)x² + a``,
  ``3 · (x + 1)² · (x² + x + 2)``;
* ``"ascii"``: ``x^3 + 2*x + 1``, ``(a + 1)*x^2 + a``,
  ``3 * (x + 1)^2 * (x^2 + x + 2)`` (lo vuelve a leer el analizador);
* ``"latex"``: ``x^{3} + 2 x + 1``, ``\\left(a + 1\\right) x^{2} + a``.

El estilo global se consulta con :func:`get_print_style`, se cambia con
:func:`set_print_style` y se cambia temporalmente con el gestor de contexto
:func:`print_style`.  El estilo es global al proceso (no por hilo).

Convenciones: los términos van por grado descendente; los coeficientes de F_p
se escriben con su representante en [0, p); el coeficiente 1 se omite (salvo en
el término independiente); un coeficiente de una extensión con varios términos
se escribe entre paréntesis.
"""
from __future__ import annotations

import contextlib
from typing import Callable, Iterator, Optional, Sequence, Tuple

ESTILOS = ("unicode", "ascii", "latex")

# Nombres admitidos (también en los especificadores de formato) -> estilo.
_ALIAS = {"unicode": "unicode", "u": "unicode",
          "ascii": "ascii", "a": "ascii",
          "latex": "latex", "l": "latex"}

_estilo = "unicode"

_SUPERINDICES = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")

# Letras griegas que en LaTeX tienen orden propia (\alpha, ...).
_GRIEGAS = frozenset(
    "alpha beta gamma delta epsilon zeta eta theta iota kappa lambda mu nu xi "
    "pi rho sigma tau upsilon phi chi psi omega".split())

#: Texto de un coeficiente no nulo y si tiene varios términos (lleva paréntesis).
Coeficiente = Tuple[str, bool]


def estilo_valido(estilo) -> Optional[str]:
    """Nombre canónico del estilo (admite los alias u, a, l) o None si no lo es."""
    return _ALIAS.get(estilo) if isinstance(estilo, str) else None


def normalizar_estilo(estilo) -> str:
    """Como :func:`estilo_valido`, pero lanza ``ValueError`` si no es válido."""
    e = estilo_valido(estilo)
    if e is None:
        raise ValueError(f"estilo de impresión desconocido: {estilo!r} "
                         f"(válidos: 'unicode', 'ascii', 'latex')")
    return e


def set_print_style(estilo: str) -> None:
    """Fija el estilo global de impresión: ``"unicode"``, ``"ascii"`` o ``"latex"``."""
    global _estilo
    _estilo = normalizar_estilo(estilo)


def get_print_style() -> str:
    """Estilo global de impresión actual."""
    return _estilo


@contextlib.contextmanager
def print_style(estilo: str) -> Iterator[str]:
    """Gestor de contexto que cambia el estilo de impresión temporalmente.

    >>> with print_style("ascii"):
    ...     print(f)          # x^3 + 2*x + 1
    """
    global _estilo
    anterior = _estilo
    _estilo = normalizar_estilo(estilo)
    try:
        yield _estilo
    finally:
        _estilo = anterior


# ---------------------------------------------------------------------------
# Piezas
# ---------------------------------------------------------------------------
def superindice(k: int) -> str:
    """``12 -> '¹²'``."""
    return str(k).translate(_SUPERINDICES)


def exponente(k: int, estilo: str) -> str:
    """Sufijo de potencia (vacío si k == 1)."""
    if k == 1:
        return ""
    if estilo == "unicode":
        return superindice(k)
    if estilo == "ascii":
        return f"^{k}"
    return f"^{{{k}}}"


def variable(nombre: str, estilo: str) -> str:
    """Nombre de variable en el estilo dado (en LaTeX, \\alpha o \\mathit{...})."""
    if estilo == "latex" and len(nombre) > 1:
        return "\\" + nombre if nombre in _GRIEGAS else f"\\mathit{{{nombre}}}"
    return nombre


def parentesis(texto: str, estilo: str) -> str:
    if estilo == "latex":
        return f"\\left({texto}\\right)"
    return f"({texto})"


def _letras_seguidas(texto: str) -> int:
    """Longitud de la racha más larga de letras de ``texto``."""
    mayor = actual = 0
    for c in texto:
        actual = actual + 1 if c.isalpha() else 0
        mayor = max(mayor, actual)
    return mayor


def _separador_producto(estilo: str, nombre_var: str, coeficiente: str) -> str:
    """Lo que va entre un coeficiente y la potencia de la variable."""
    if estilo == "ascii":
        return "*"
    if estilo == "latex":
        return " "
    # Unicode: yuxtaposición (2x, ax, (a + 1)x²), salvo si la variable o el
    # coeficiente usan nombres de varias letras (2·xy, alpha·x).
    return "·" if len(nombre_var) > 1 or _letras_seguidas(coeficiente) > 1 else ""


def polinomio(coefs: Sequence, nombre_var: str, estilo: str,
              formato: Callable[[object], Optional[Coeficiente]]) -> str:
    """Texto de ``sum(coefs[k]·var^k)`` por grado descendente.

    ``formato(c)`` devuelve ``(texto, compuesto)`` para un coeficiente no nulo o
    ``None`` si es cero.  El polinomio sin términos se escribe ``0``.
    """
    var = variable(nombre_var, estilo)
    partes = []
    for k in range(len(coefs) - 1, -1, -1):
        f = formato(coefs[k])
        if f is None:
            continue
        texto, compuesto = f
        if k == 0:
            partes.append(texto)
            continue
        pot = var + exponente(k, estilo)
        if texto == "1":
            partes.append(pot)
        else:
            sep = _separador_producto(estilo, nombre_var, texto)
            partes.append((parentesis(texto, estilo) if compuesto else texto) + sep + pot)
    return " + ".join(partes) if partes else "0"


def formato_entero(c: int) -> Optional[Coeficiente]:
    """Formato de un coeficiente entero (de F_p) para :func:`polinomio`."""
    return (str(c), False) if c else None


def polinomio_entero(coefs: Sequence[int], nombre_var: str, estilo: str) -> str:
    """Polinomio de coeficientes enteros (en [0, p)) en orden creciente."""
    return polinomio(coefs, nombre_var, estilo, formato_entero)


def nombre_cuerpo(q: int, estilo: str) -> str:
    """Nombre corto de F_q: ``GF(9)`` o ``\\mathbb{F}_{9}``."""
    return f"\\mathbb{{F}}_{{{q}}}" if estilo == "latex" else f"GF({q})"


def producto(unidad: Coeficiente, factores: Sequence[Tuple[str, bool, int]],
             estilo: str) -> str:
    """Texto de ``unidad · f1^e1 · f2^e2 ...``.

    La unidad 1 se omite si hay factores; los factores con varios términos van
    entre paréntesis salvo que sean lo único que se escribe (``x + 1``).
    """
    sep = {"unicode": " · ", "ascii": " * ", "latex": " \\cdot "}[estilo]
    texto_u, compuesta = unidad
    con_unidad = texto_u != "1" or not factores
    solo = not con_unidad and len(factores) == 1 and factores[0][2] == 1
    partes = []
    if con_unidad:
        partes.append(parentesis(texto_u, estilo) if compuesta and factores else texto_u)
    for texto, compuesto, e in factores:
        base = parentesis(texto, estilo) if compuesto and not solo else texto
        partes.append(base + exponente(e, estilo))
    return sep.join(partes)


def formatear(obj, spec: str) -> str:
    """Implementación común de ``__format__``.

    ``spec`` vacío: estilo global; ``unicode``/``u``, ``ascii``/``a``,
    ``latex``/``l``: ese estilo; cualquier otro: especificador de cadena
    (anchura, alineación...) aplicado al texto en el estilo global.
    """
    estilo = estilo_valido(spec)
    if estilo is not None:
        return obj._texto(estilo)
    return format(obj._texto(_estilo), spec)


class Imprimible:
    """Mezcla con ``__str__``, ``__format__``, ``_repr_latex_`` y ``_repr_pretty_``.

    Las clases que la usan definen ``_texto(estilo)``.
    """

    __slots__ = ()

    def _texto(self, estilo: str) -> str:  # pragma: no cover - abstracto
        raise NotImplementedError

    def __str__(self) -> str:
        return self._texto(_estilo)

    def __format__(self, spec: str) -> str:
        return formatear(self, spec)

    def _repr_latex_(self) -> str:
        """Representación LaTeX para Jupyter (entre $...$)."""
        return f"${self._texto('latex')}$"

    def _repr_pretty_(self, p, cycle: bool) -> None:
        """Representación para IPython: el texto en el estilo global."""
        p.text(str(self))
