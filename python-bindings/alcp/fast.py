"""Algoritmos rápidos (Práctica 3): DFT, matrices de Toeplitz y división.

Envuelven las funciones del módulo nativo con conversión generosa: los
vectores son listas (o tuplas) de valores convertibles al cuerpo (enteros,
cadenas, elementos...) y los resultados son listas de elementos.

* :func:`dft` / :func:`idft`: transformada de Fourier discreta de longitud
  n = 2^k (Cooley-Tukey).  Si no se da la raíz, se usa
  :func:`root_of_unity` (determinista).
* ``toeplitz_*``: productos e inversas de matrices de Toeplitz.
* :func:`divmod_toeplitz`: división euclídea vía sistemas de Toeplitz.
"""
from __future__ import annotations

from ._nativo import nativo
from ._util import entero
from .field import FiniteField
from .poly import Polynomial


def _cuerpo(K) -> FiniteField:
    if not isinstance(K, FiniteField):
        raise TypeError(f"se esperaba un cuerpo finito (alcp.GF), no {type(K).__name__!r}")
    return K


def _valores(K: FiniteField, a, que: str) -> list:
    if isinstance(a, (str, dict)) or not hasattr(a, "__iter__"):
        raise TypeError(f"{que} debe ser una lista de elementos, no {type(a).__name__!r}")
    return [K._valor(x) for x in a]


def _elementos(K: FiniteField, vs) -> list:
    return [K._crear(v) for v in vs]


def root_of_unity(K: FiniteField, n):
    """Raíz primitiva n-ésima de la unidad de K (ver ``FiniteField.root_of_unity``)."""
    return _cuerpo(K).root_of_unity(n)


def _transformada(funcion, K, a, root):
    K = _cuerpo(K)
    vs = _valores(K, a, "a")
    n = len(vs)
    if n == 0 or n & (n - 1):
        raise ValueError(f"la longitud de la DFT debe ser una potencia de 2, no {n}")
    w = K.root_of_unity(n) if root is None else K(root)
    return _elementos(K, funcion(K._nat, w._v, n.bit_length() - 1, vs))


def dft(K: FiniteField, a, root=None) -> list:
    """DFT de ``a`` (longitud n = 2^k) en K: ``[sum(a_j·w^(ij)) for i in range(n)]``.

    ``root`` debe tener orden exactamente n; por defecto, ``K.root_of_unity(n)``
    (``ValueError`` si K no tiene raíces primitivas n-ésimas, es decir, si n no
    divide a q - 1).
    """
    return _transformada(nativo.dft, K, a, root)


def idft(K: FiniteField, a, root=None) -> list:
    """Transformada inversa: ``idft(K, dft(K, a, w), w) == a`` (misma raíz w)."""
    return _transformada(nativo.idft, K, a, root)


def _n(n, defecto: int) -> int:
    return defecto if n is None else entero(n, "n")


def toeplitz_lower_vec(K: FiniteField, a, b, n=None) -> list:
    """``T·b`` con T triangular inferior de Toeplitz n×n de primera columna ``a``.

    Por defecto n = len(b); ``a`` debe tener al menos n elementos.
    """
    K = _cuerpo(K)
    vb = _valores(K, b, "b")
    vs = nativo.toeplitz_lower_vec(K._nat, _n(n, len(vb)), _valores(K, a, "a"), vb)
    return _elementos(K, vs)


def toeplitz_upper_vec(K: FiniteField, a, b, n=None) -> list:
    """``T·b`` con T triangular superior de Toeplitz n×n de primera fila ``a``."""
    K = _cuerpo(K)
    vb = _valores(K, b, "b")
    vs = nativo.toeplitz_upper_vec(K._nat, _n(n, len(vb)), _valores(K, a, "a"), vb)
    return _elementos(K, vs)


def toeplitz_vec(K: FiniteField, a, b, n=None) -> list:
    """``T·b`` con T de Toeplitz n×n completa.

    ``a`` (longitud 2n - 1) es la primera fila seguida de la primera columna sin
    la esquina: ``T[i][j] = a[j - i]`` si j >= i y ``a[n - 1 + i - j]`` si j < i.
    """
    K = _cuerpo(K)
    vb = _valores(K, b, "b")
    vs = nativo.toeplitz_vec(K._nat, _n(n, len(vb)), _valores(K, a, "a"), vb)
    return _elementos(K, vs)


def _inversa(funcion, K, a, n) -> list:
    K = _cuerpo(K)
    va = _valores(K, a, "a")
    if va and va[0] == 0:
        raise ZeroDivisionError("la matriz de Toeplitz es singular: a[0] = 0")
    return _elementos(K, funcion(K._nat, _n(n, len(va)), va))


def toeplitz_lower_inverse(K: FiniteField, a, n=None) -> list:
    """Primera columna de T⁻¹ (T triangular inferior n×n de primera columna ``a``).

    Por defecto n = len(a).  ``ZeroDivisionError`` si a[0] = 0.
    """
    return _inversa(nativo.toeplitz_lower_inverse, K, a, n)


def toeplitz_upper_inverse(K: FiniteField, a, n=None) -> list:
    """Primera fila de T⁻¹ (T triangular superior n×n de primera fila ``a``).

    Por defecto n = len(a).  ``ZeroDivisionError`` si a[0] = 0.
    """
    return _inversa(nativo.toeplitz_upper_inverse, K, a, n)


def divmod_toeplitz(f: Polynomial, g):
    """División euclídea ``f = q·g + r`` resuelta con sistemas de Toeplitz: ``(q, r)``.

    ``g`` se convierte al anillo de ``f`` (admite enteros, cadenas, listas...).
    """
    if not isinstance(f, Polynomial):
        raise TypeError(f"f debe ser un polinomio, no {type(f).__name__!r}")
    S, a, b = f._otro(g)
    if b.is_zero():
        raise ZeroDivisionError("división por el polinomio cero")
    q, r = nativo.divmod_toeplitz(a, b)
    return S._envolver(q), S._envolver(r)
