"""Pruebas de los algoritmos rápidos de alto nivel (alcp.fast).

DFT/IDFT contra la definición ``A_i = sum a_j w^(ij)``, productos e inversas de
Toeplitz contra matrices construidas explícitamente, y la división por
sistemas de Toeplitz contra ``divmod``; en cuerpos primos y en extensiones.
"""
from __future__ import annotations

import random

import pytest

import alcp
from alcp import (dft, divmod_toeplitz, idft, root_of_unity, toeplitz_lower_inverse,
                  toeplitz_lower_vec, toeplitz_upper_inverse, toeplitz_upper_vec, toeplitz_vec)

GF = alcp.GF


@pytest.fixture(autouse=True)
def semilla():
    alcp.set_seed(7)


def matvec(K, T, b):
    return [sum((T[i][j] * K(b[j]) for j in range(len(b))), K.zero) for i in range(len(T))]


def inferior(K, a, n):
    return [[K(a[i - j]) if i >= j else K.zero for j in range(n)] for i in range(n)]


def superior(K, a, n):
    return [[K(a[j - i]) if j >= i else K.zero for j in range(n)] for i in range(n)]


def completa(K, a, n):
    return [[K(a[j - i]) if j >= i else K(a[n - 1 + i - j]) for j in range(n)] for i in range(n)]


def producto(K, A, B):
    n = len(A)
    return [[sum((A[i][k] * B[k][j] for k in range(n)), K.zero) for j in range(n)]
            for i in range(n)]


# ---------------------------------------------------------------------------
# DFT
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("q, n", [(17, 1), (17, 2), (17, 8), (17, 16), (9, 8), (25, 8),
                                  (998244353, 32), (2 ** 4, 1)])
def test_dft_contra_la_definicion(q, n):
    K = GF(q)
    a = [K.random() for _ in range(n)]
    w = root_of_unity(K, n)
    assert w == K.root_of_unity(n) and w.multiplicative_order() == n
    esperado = [sum((a[j] * w ** (i * j) for j in range(n)), K.zero) for i in range(n)]
    A = dft(K, a)
    assert A == esperado and all(x.field is K for x in A)
    assert dft(K, a, w) == esperado and idft(K, A) == a and idft(K, A, w) == a


def test_dft_con_valores_convertibles():
    K = GF(17)
    assert dft(K, [1, "2", 3, K(4)]) == dft(K, [1, 2, 3, 4])
    assert dft(K, (1, 2), root=-1) == [3, 16] and dft(K, [1, 2], root="16") == [3, 16]
    assert dft(K, [1, 1, 1, 1], root=4) == [4, 0, 0, 0]
    assert dft(K, range(4)) == dft(K, [0, 1, 2, 3])


@pytest.mark.parametrize("a", [[], [1, 2, 3], [1] * 6])
def test_dft_longitud_no_potencia_de_2(a):
    with pytest.raises(ValueError, match="la longitud de la DFT debe ser una potencia de 2"):
        dft(GF(17), a)


def test_dft_errores():
    with pytest.raises(ValueError, match="no tiene raíces primitivas de la unidad de orden 4"):
        dft(GF(7), [1, 2, 3, 4])
    with pytest.raises(ValueError, match="no tiene orden exactamente"):
        dft(GF(17), [1, 2, 3, 4], root=2)
    with pytest.raises(TypeError, match="se esperaba un cuerpo finito"):
        dft(17, [1, 2])
    with pytest.raises(TypeError, match="a debe ser una lista de elementos, no 'str'"):
        idft(GF(17), "12")
    with pytest.raises(TypeError, match="no 'int'"):
        dft(GF(17), 5)
    with pytest.raises(TypeError, match="no 'dict'"):
        dft(GF(17), {0: 1})


# ---------------------------------------------------------------------------
# Toeplitz
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("q", [7, 9, 101])
def test_productos_de_toeplitz(q):
    K = GF(q)
    rng = random.Random(q)
    for n in (1, 2, 3, 5, 8):
        a = [K.from_index(rng.randrange(q)) for _ in range(2 * n - 1)]
        b = [K.from_index(rng.randrange(q)) for _ in range(n)]
        assert toeplitz_lower_vec(K, a[:n], b) == matvec(K, inferior(K, a, n), b)
        assert toeplitz_lower_vec(K, a, b, n) == matvec(K, inferior(K, a, n), b)
        assert toeplitz_upper_vec(K, a[:n], b) == matvec(K, superior(K, a, n), b)
        assert toeplitz_upper_vec(K, a, b, n=n) == matvec(K, superior(K, a, n), b)
        assert toeplitz_vec(K, a, b) == matvec(K, completa(K, a, n), b)


@pytest.mark.parametrize("q", [7, 9, 101])
def test_inversas_de_toeplitz(q):
    K = GF(q)
    rng = random.Random(q + 1)
    for n in (1, 2, 3, 4, 7):
        a = [K.from_index(rng.randrange(1, q))] + [K.from_index(rng.randrange(q))
                                                   for _ in range(n)]
        identidad = [[K.one if i == j else K.zero for j in range(n)] for i in range(n)]
        col = toeplitz_lower_inverse(K, a, n)
        assert producto(K, inferior(K, a, n), inferior(K, col, n)) == identidad
        fila = toeplitz_upper_inverse(K, a[:n])
        assert producto(K, superior(K, a, n), superior(K, fila, n)) == identidad


def test_toeplitz_con_valores_convertibles():
    K = GF(7)
    assert toeplitz_lower_vec(K, ["1", 2], (3, "4")) == [3, 10 % 7]
    assert toeplitz_lower_inverse(K, [2]) == [4] and toeplitz_upper_inverse(K, ("3",)) == [5]


def test_toeplitz_errores():
    K = GF(7)
    with pytest.raises(ZeroDivisionError, match="la matriz de Toeplitz es singular: a\\[0\\] = 0"):
        toeplitz_lower_inverse(K, [0, 1])
    with pytest.raises(ZeroDivisionError, match="singular"):
        toeplitz_upper_inverse(K, [7, 1])
    with pytest.raises(ValueError, match="n debe ser >= 1"):
        toeplitz_lower_inverse(K, [])
    with pytest.raises(ValueError, match="debe tener al menos n = 3"):
        toeplitz_lower_vec(K, [1, 2], [1, 2, 3])
    with pytest.raises(ValueError, match="debe tener exactamente"):
        toeplitz_vec(K, [1, 2], [1, 2])
    with pytest.raises(ValueError, match="debe tener exactamente"):
        toeplitz_upper_vec(K, [1, 2, 3], [1, 2], n=3)
    with pytest.raises(TypeError, match="n debe ser un entero"):
        toeplitz_upper_vec(K, [1, 2], [1, 2], n=2.0)
    with pytest.raises(TypeError, match="b debe ser una lista de elementos"):
        toeplitz_vec(K, [1], "1")
    with pytest.raises(TypeError, match="se esperaba un cuerpo finito"):
        toeplitz_upper_inverse("GF(7)", [1])


# ---------------------------------------------------------------------------
# División por Toeplitz
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("q", [7, 9, 998244353])
def test_divmod_toeplitz(q):
    R = GF(q)["x"]
    for df, dg in ((10, 3), (5, 5), (3, 7), (20, 1), (8, 0)):
        f, g = R.random(df), R.random(dg)
        cociente, resto = divmod_toeplitz(f, g)
        assert (cociente, resto) == divmod(f, g)
        assert cociente.ring is R and resto.ring is R


def test_divmod_toeplitz_conversion_y_errores():
    R = GF(3)["x"]
    f = R("x^4 + 2x + 1")
    assert divmod_toeplitz(f, "x + 1") == divmod(f, R("x + 1"))
    q, r = divmod_toeplitz(f, GF(9)["x"]("x + a"))          # se promociona a GF(9)[x]
    assert q.ring is GF(9)["x"] and q * "x + a" + r == f
    with pytest.raises(ZeroDivisionError, match="división por el polinomio cero"):
        divmod_toeplitz(f, 0)
    with pytest.raises(TypeError, match="f debe ser un polinomio, no 'str'"):
        divmod_toeplitz("x", f)
    with pytest.raises(TypeError, match="no se puede convertir un valor de tipo 'NoneType'"):
        divmod_toeplitz(f, None)
