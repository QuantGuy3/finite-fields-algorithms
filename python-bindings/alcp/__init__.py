"""alcp: cuerpos finitos, polinomios y algoritmos rápidos (núcleo C++/GMP).

Capa de alto nivel sobre el módulo nativo ``alcp.alcp_native`` con conversión
de tipos generosa y representación cuidada (docs/diseno.md, §6)::

    >>> import alcp
    >>> F = alcp.GF(7)                      # cuerpo primo
    >>> K = alcp.GF(9, modulus="x^2+1")     # extensión F_3[a]/(a² + 1)
    >>> a = K.gen                           # clase de x; se imprime «a»
    >>> R = K["x"]                          # anillo de polinomios K[x]
    >>> f = R("(a+1)x^3 + 2x + a")          # cadenas, listas, dicts, enteros...
    >>> print(f + 2)
    (a + 1)x³ + 2x + a + 2
    >>> print(f.factor())
    (a + 1) · (x + 1) · (x² + 2x + 2a + 2)

Sin ``modulus``, ``GF(p^n)`` usa un módulo determinista (el primer mónico
irreducible en un orden fijo, ver :func:`GF`); ``modulus="random"`` elige uno
aleatorio, reproducible con :func:`set_seed`.

Contenido:

* cuerpos y elementos: :func:`GF`, :class:`FiniteField` (:class:`PrimeField`,
  :class:`ExtensionField`), :class:`FieldElement` (:class:`PrimeFieldElement`,
  :class:`ExtensionFieldElement`);
* polinomios: :class:`PolynomialRing`, :class:`Polynomial`, :class:`Factorization`;
* algoritmos rápidos: :func:`dft`, :func:`idft`, :func:`root_of_unity`,
  ``toeplitz_*``, :func:`divmod_toeplitz`;
* utilidades: :func:`is_prime`, :func:`prime_power_factor`, :func:`next_prime`,
  :func:`set_seed`;
* representación: :func:`set_print_style`, :func:`get_print_style`,
  :func:`print_style` (estilos "unicode", "ascii" y "latex").

El módulo nativo se compila con ``python build.py module``; si falta, importar
el paquete lanza ``ImportError`` con esa indicación.
"""
from __future__ import annotations

from ._nativo import nativo as alcp_native
from ._pretty import get_print_style, print_style, set_print_style
from ._util import is_prime, next_prime, prime_power_factor, set_seed
from .fast import (dft, divmod_toeplitz, idft, root_of_unity, toeplitz_lower_inverse,
                   toeplitz_lower_vec, toeplitz_upper_inverse, toeplitz_upper_vec,
                   toeplitz_vec)
from .field import (GF, ExtensionField, ExtensionFieldElement, FieldElement, FiniteField,
                    PrimeField, PrimeFieldElement)
from .poly import Factorization, Polynomial, PolynomialRing

__version__ = "0.2.0"

__all__ = [
    # cuerpos y elementos
    "GF", "FiniteField", "PrimeField", "ExtensionField",
    "FieldElement", "PrimeFieldElement", "ExtensionFieldElement",
    # polinomios
    "PolynomialRing", "Polynomial", "Factorization",
    # algoritmos rápidos
    "dft", "idft", "root_of_unity", "toeplitz_lower_vec", "toeplitz_upper_vec",
    "toeplitz_vec", "toeplitz_lower_inverse", "toeplitz_upper_inverse", "divmod_toeplitz",
    # utilidades
    "is_prime", "prime_power_factor", "next_prime", "set_seed",
    # representación
    "set_print_style", "get_print_style", "print_style",
    # módulo nativo
    "alcp_native",
]
