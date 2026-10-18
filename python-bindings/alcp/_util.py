"""Utilidades internas: conversión a entero, primos y factorización de enteros."""
from __future__ import annotations

import operator

from ._nativo import nativo

#: Las divisiones de prueba de :func:`factorizar` no pasan de este divisor.
LIMITE_DIVISION = 1 << 20


def entero(x, que: str = "el argumento") -> int:
    """Convierte ``x`` en ``int`` con ``operator.index``.

    Acepta ``int``, ``bool`` y cualquier objeto con ``__index__`` (p. ej.
    ``numpy.int64``).  Si ``x`` no es un entero, lanza ``TypeError`` con un
    mensaje en español; si ``__index__`` lanza otra excepción, se propaga.
    """
    if type(x) is int:
        return x
    try:
        return operator.index(x)
    except TypeError:
        raise TypeError(f"{que} debe ser un entero, no un valor de tipo "
                        f"{type(x).__name__!r}") from None


def is_prime(n) -> bool:
    """¿Es ``n`` primo?  Test probabilístico de GMP (Baillie-PSW + Miller-Rabin)."""
    return nativo.is_prime(entero(n, "n"))


def prime_power_factor(q) -> tuple[int, int]:
    """Devuelve ``(p, n)`` con ``q = p**n`` y ``p`` primo.

    ``ValueError`` si ``q`` no es potencia de un primo (en particular, si q < 2).
    """
    q = entero(q, "q")
    try:
        p, n = nativo.prime_power_factor(q)
    except ValueError:
        raise ValueError(f"{q} no es potencia de un número primo") from None
    return p, n


def next_prime(n) -> int:
    """Menor primo estrictamente mayor que ``n``."""
    c = max(entero(n, "n") + 1, 2)
    while not nativo.is_prime(c):
        c += 1
    return c


def set_seed(seed) -> None:
    """Fija la semilla (entero >= 0) del generador aleatorio global del núcleo.

    Afecta a ``K.random()``, ``R.random(d)``, ``R.irreducible(d)``,
    ``GF(q, modulus="random")`` y a la factorización (Cantor-Zassenhaus), que
    pasan a ser reproducibles.
    """
    seed = entero(seed, "la semilla")
    if seed < 0:
        raise ValueError("la semilla debe ser un entero >= 0")
    nativo.set_seed(seed)


def factorizar(n: int) -> dict[int, int]:
    """Factorización ``{primo: exponente}`` de un entero ``n >= 1``.

    Divide por 2, 3 y los números 6k ± 1 hasta :data:`LIMITE_DIVISION`; cada vez
    que el cofactor cambia se comprueba si ya es primo (test de GMP), así que un
    único factor primo grande no cuesta nada.  Si queda un cofactor compuesto sin
    factores por debajo del límite, lanza ``ValueError``.
    """
    factores: dict[int, int] = {}

    def dividir(d: int) -> bool:
        nonlocal n
        if n % d:
            return False
        while n % d == 0:
            n //= d
            factores[d] = factores.get(d, 0) + 1
        return True

    dividir(2)
    dividir(3)
    d, paso = 5, 2
    primo = nativo.is_prime(n)
    while not primo and d * d <= n:
        if d > LIMITE_DIVISION:
            raise ValueError(
                f"no se puede factorizar el entero: el cofactor {n} es compuesto y "
                f"no tiene factores primos menores que {LIMITE_DIVISION}")
        if dividir(d):
            primo = nativo.is_prime(n)
        d += paso
        paso = 6 - paso
    if n > 1:
        factores[n] = factores.get(n, 0) + 1
    return factores
