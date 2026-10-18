"""Implementaciones de referencia en Python puro para las pruebas del paquete ``alcp``.

Son lentas pero evidentes: aritmética de polinomios con coeficientes enteros
módulo p (listas en orden creciente) y de extensiones F_p[a]/(g) con elementos
como tuplas de coordenadas.  Sirven de oráculo independiente del núcleo C++.
"""
from __future__ import annotations

import itertools


def recortar(c):
    """Quita los ceros finales de una lista de coeficientes."""
    c = list(c)
    while c and c[-1] == 0:
        c.pop()
    return c


def pol_suma(a, b, p):
    n = max(len(a), len(b))
    a, b = list(a) + [0] * (n - len(a)), list(b) + [0] * (n - len(b))
    return recortar([(x + y) % p for x, y in zip(a, b)])


def pol_neg(a, p):
    return recortar([(-x) % p for x in a])


def pol_resta(a, b, p):
    return pol_suma(a, pol_neg(b, p), p)


def pol_mul(a, b, p):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] = (r[i + j] + x * y) % p
    return recortar(r)


def pol_divmod(a, b, p):
    """División euclídea en F_p[x] (b no nulo)."""
    a, b = recortar(a), recortar(b)
    q = [0] * max(len(a) - len(b) + 1, 0)
    inv = pow(b[-1], -1, p)
    r = list(a)
    while len(r) >= len(b) and r:
        c = r[-1] * inv % p
        d = len(r) - len(b)
        q[d] = c
        for i, y in enumerate(b):
            r[i + d] = (r[i + d] - c * y) % p
        r = recortar(r)
    return recortar(q), r


def pol_eval(a, x, p):
    v = 0
    for c in reversed(a):
        v = (v * x + c) % p
    return v


def mobius(n):
    """Función de Möbius."""
    m, d = 1, 2
    while d * d <= n:
        if n % d == 0:
            n //= d
            if n % d == 0:
                return 0
            m = -m
        d += 1
    return -m if n > 1 else m


def irreducibles_gauss(q, n):
    """Número de polinomios mónicos irreducibles de grado n sobre F_q (fórmula de Gauss)."""
    return sum(mobius(d) * q ** (n // d) for d in range(1, n + 1) if n % d == 0) // n


def monicos(p_o_q, n):
    """Todas las listas de coeficientes de los mónicos de grado n (coeficientes en range)."""
    for cola in itertools.product(range(p_o_q), repeat=n):
        yield list(cola) + [1]


class Extension:
    """F_p[a]/(g) con elementos como tuplas (c0, ..., c_{n-1})."""

    def __init__(self, p, g):
        self.p, self.g = p, list(g)
        self.n = len(g) - 1

    def desde_indice(self, i):
        c = []
        for _ in range(self.n):
            i, r = divmod(i, self.p)
            c.append(r)
        return tuple(c)

    def indice(self, c):
        return sum(x * self.p ** k for k, x in enumerate(c))

    def reducir(self, c):
        _, r = pol_divmod(c, self.g, self.p)
        return tuple(r + [0] * (self.n - len(r)))

    def suma(self, x, y):
        return tuple((a + b) % self.p for a, b in zip(x, y))

    def mul(self, x, y):
        return self.reducir(pol_mul(recortar(x), recortar(y), self.p))

    def potencia(self, x, k):
        r = self.reducir([1])
        for _ in range(k):
            r = self.mul(r, x)
        return r
