"""Fase 7, propuesta 7: los cuerpos de extensión rápidos FqTab y FqSwar en el binding.

Con ``backend='auto'`` y p < 2^62, ``Fq`` elige ahora FqTab (q <= 2^16) o
FqSwar (p pequeño, 2·n·(p-1)^2 + p < 2^16) cuando caben, y Fq64 en el resto.
Desde Python no debe notarse nada salvo la velocidad: se compara cada cuerpo
con el mismo cuerpo forzado a GMP (``backend='gmp'``, que da por contrato lo
mismo que Fq64) en elementos, índices, coeficientes, polinomios, factorización,
algoritmos rápidos, mensajes de error y consumo del generador aleatorio
(cada prueba parte de una semilla fija, con cuerpos recién creados: «en frío»).
"""
from __future__ import annotations

import itertools
import random

import pytest

from alcp import alcp_native as al

AES = [1, 1, 0, 1, 1, 0, 0, 0, 1]
M3E30 = [2, 0, 2, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0, 2, 2, 1, 0, 1, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1]

# (p, módulo, tipo interno esperado).  F_(97^2) y F_(193^2) admiten FFT de
# 2^6 y 2^7; F_(3^30) y F_(2^20) van por FqSwar; F_(101^5) no cabe en ninguno.
CUERPOS = [
    (2, AES, "tab"),
    (2, [1, 1, 1], "tab"),
    (3, [1, 0, 1], "tab"),
    (7, [3, 4, 1, 2, 1], "tab"),
    (97, [5, 0, 1], "tab"),
    (193, [5, 0, 1], "tab"),
    (3, M3E30, "swar"),
]


def _testigo():
    return al.Fp(2**61 - 1).random()


def par(p, g):
    return al.Fq(p, g), al.Fq(p, g, backend="gmp")


@pytest.mark.parametrize("p, g, _tipo", CUERPOS)
def test_datos_repr_igualdad(p, g, _tipo):
    K, G = par(p, g)
    assert K.backend == "u64" and G.backend == "gmp"
    assert (K.p, K.n, K.q, K.modulus) == (G.p, G.n, G.q, G.modulus)
    assert repr(K) == repr(G).replace(", backend='gmp'", "")
    assert eval(repr(K), {"Fq": al.Fq}) == K
    assert K == al.Fq(p, g) and K != G and hash(K) == hash(G) == hash(al.Fq(p, g))
    assert K.fft_capacity() == G.fft_capacity()


@pytest.mark.parametrize("p, g, _tipo", CUERPOS)
def test_elementos(p, g, _tipo):
    K, G = par(p, g)
    rng = random.Random(p * 1000 + len(g))
    q = K.q
    xs = list(range(q)) if q <= 64 else [0, 1, q - 1] + [rng.randrange(q) for _ in range(40)]
    for x in xs:
        assert K.to_coeffs(x) == G.to_coeffs(x)
        assert K.from_coeffs(G.to_coeffs(x)) == x == K.from_index(x)
        assert K.neg(x) == G.neg(x)
        for k in (0, 1, 3, -1, -2, q + 5, 10**30):
            if x == 0 and k < 0:
                with pytest.raises(ZeroDivisionError) as e1:
                    K.pow(x, k)
                with pytest.raises(ZeroDivisionError) as e2:
                    G.pow(x, k)
                assert str(e1.value) == str(e2.value)
            else:
                assert K.pow(x, k) == G.pow(x, k)
        if x:
            assert K.inv(x) == G.inv(x)
        for y in xs[:20]:
            assert K.add(x, y) == G.add(x, y) and K.sub(x, y) == G.sub(x, y) and K.mul(x, y) == G.mul(x, y)
            if y:
                assert K.div(x, y) == G.div(x, y)
    for n in (0, 1, -1, 5, 10**40):
        assert K.from_int(n) == G.from_int(n)
        assert K.from_index(n) == G.from_index(n)
    assert K.from_coeffs([5, -3, 7, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]) == G.from_coeffs(
        [5, -3, 7, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
    for f in (lambda H: H.inv(0), lambda H: H.div(3, 0)):
        with pytest.raises(ZeroDivisionError) as e1:
            f(K)
        with pytest.raises(ZeroDivisionError) as e2:
            f(G)
        assert str(e1.value) == str(e2.value)


def bateria(K, d1, d2):
    """Cálculos que consumen el generador: todo sale de random/random_poly."""
    out = [K.random() for _ in range(5)]
    a, b = K.random_poly(d1), K.random_poly(d2)
    out += [(a * b).coeffs(), a.mul(b, "school").coeffs(), a.mul(b, "karatsuba").coeffs(),
            a.mul(b, "kronecker").coeffs(), a.mul(b, "fft").coeffs()]
    q, r = divmod(a, b)
    out += [q.coeffs(), r.coeffs(), (a // b).coeffs(), (a % b).coeffs()]
    out += [a.gcd(b).coeffs()]
    out += [[x.coeffs() for x in a.gcd_ext(b)]]
    out += [a.pow_mod(K.q + 3, b).coeffs(), a.compose(b).coeffs(), b.is_irreducible(), b.is_squarefree()]
    out += [a.evaluate(7), a.leading(), a.monic().coeffs(), a.derivative().coeffs(), str(a), repr(a)]
    out += [a.scalar_mul(5).coeffs(), (-a).coeffs(), (a + b - b).coeffs(), (b ** 3).coeffs()]
    out.append(K.fft_capacity() > 0 and K.fft_base_root())
    return out


def factores(K, d):
    f = K.random_poly(d)
    g = K.random_poly(2) * K.random_poly(2) * K.random_poly(1) ** 3 * K.random_poly(3)
    out = []
    for h in (f, g):
        fs, lc = h.factor()
        out.append(([(u.coeffs(), e) for u, e in fs], lc))
        out.append(h.radical().coeffs())
        out.append([(dd, u.coeffs()) for dd, u in h.radical().monic().distinct_degree()])
    for dd, u in g.radical().monic().distinct_degree():
        out.append([x.coeffs() for x in u.equal_degree(dd)])
        out.append(g.multiplicity(u))
    return out


def en_frio(p, g, calc, semilla):
    al.set_seed(semilla)
    ra = calc(al.Fq(p, g))
    ta = _testigo()
    al.set_seed(semilla)
    rb = calc(al.Fq(p, g, backend="gmp"))
    tb = _testigo()
    return ra, ta, rb, tb


@pytest.mark.parametrize("p, g, _tipo", CUERPOS)
@pytest.mark.parametrize("d1, d2", [(9, 4), (40, 33), (90, 70)])
def test_polinomios_en_frio(p, g, _tipo, d1, d2):
    if len(g) > 20 and d1 > 40:
        pytest.skip("F_(3^30) con GMP es lento en grados altos")
    ra, ta, rb, tb = en_frio(p, g, lambda K: bateria(K, d1, d2), 17 + d1)
    assert ra == rb and ta == tb


@pytest.mark.parametrize("p, g, _tipo", CUERPOS)
def test_factor_en_frio(p, g, _tipo):
    d = 5 if len(g) > 20 else 30
    ra, ta, rb, tb = en_frio(p, g, lambda K: factores(K, d), 91)
    assert ra == rb and ta == tb


@pytest.mark.parametrize("q", [2**8, 5**3, 7**4, 2**16, 3**30, 2**20])
def test_modulo_aleatorio(q):
    al.set_seed(5)
    K = al.Fq(q)
    t1 = _testigo()
    al.set_seed(5)
    G = al.Fq(q, backend="gmp")
    assert K.modulus == G.modulus and t1 == _testigo()
    al.set_seed(6)
    a = K.random_poly(6).coeffs()
    al.set_seed(6)
    assert a == G.random_poly(6).coeffs()


def test_rapidos_iguales():
    for p, g in [(97, [5, 0, 1]), (3, M3E30), (2, AES)]:
        K, G = par(p, g)
        cap = K.fft_capacity()
        if cap >= 3:
            al.set_seed(40)  # la raíz base sale del generador: misma semilla
            raiz_k = K.pow(K.fft_base_root(), 2 ** (cap - 3))
            al.set_seed(40)
            raiz_g = G.pow(G.fft_base_root(), 2 ** (cap - 3))
            assert raiz_k == raiz_g
            a = [(7 * i + 3) % K.q for i in range(8)]
            assert al.dft(K, raiz_k, 3, a) == al.dft(G, raiz_g, 3, a)
            assert al.idft(K, raiz_k, 3, a) == al.idft(G, raiz_g, 3, a)
        a, b = [3, 1, 4, 1, 5], [9, 2, 6]
        assert al.toeplitz_lower_vec(K, 3, a, b) == al.toeplitz_lower_vec(G, 3, a, b)
        assert al.toeplitz_upper_vec(K, 3, a, b) == al.toeplitz_upper_vec(G, 3, a, b)
        assert al.toeplitz_vec(K, 3, a, b) == al.toeplitz_vec(G, 3, a, b)
        assert al.toeplitz_lower_inverse(K, 3, [1, 2, 3]) == al.toeplitz_lower_inverse(G, 3, [1, 2, 3])
        f, h = K.poly(list(range(1, 30))), G.poly(list(range(1, 30)))
        d1, d2 = K.poly([1, 2, 3, 1]), G.poly([1, 2, 3, 1])
        assert [x.coeffs() for x in al.divmod_toeplitz(f, d1)] == [x.coeffs() for x in al.divmod_toeplitz(h, d2)]


def test_errores_y_mezclas():
    # Módulo reducible y p no primo: mismos mensajes que con GMP.
    for args in [(2, [1, 0, 1]), (4, [1, 1, 1]), (3, [1]), (3, [])]:
        with pytest.raises(ValueError) as e1:
            al.Fq(*args)
        with pytest.raises(ValueError) as e2:
            al.Fq(*args, backend="gmp")
        assert str(e1.value) == str(e2.value)
    # Polinomios de cuerpos de distinto tipo interno (o distinto backend) no se mezclan.
    K, G = par(2, AES)
    L = al.Fq(3, M3E30)
    for x, y in [(K.poly([1, 2]), G.poly([1, 2])), (K.poly([1, 2]), L.poly([1, 2]))]:
        assert x != y
        with pytest.raises(ValueError):
            x + y
    assert K.poly([1, 2]) == K.poly([1, 2]) and hash(K.poly([1, 2])) == hash(G.poly([1, 2]))


def test_cuerpos_grandes_siguen_en_fq64():
    # F_(101^5): ni FqTab (q > 2^16) ni FqSwar (2·5·100^2 >= 2^16); F_(2^70): n > 64.
    for p, n in [(101, 5), (2, 70)]:
        al.set_seed(3)
        K = al.Fq(p**n)
        al.set_seed(3)
        G = al.Fq(p**n, backend="gmp")
        assert K.modulus == G.modulus
        x, y = K.random(), K.random()
        assert K.mul(x, y) == G.mul(x, y)


def test_enumeracion_completa_pequena():
    for p, g in [(2, [1, 1, 1]), (3, [1, 0, 1]), (5, [2, 1])]:
        K, G = par(p, g)
        for x, y in itertools.product(range(K.q), repeat=2):
            assert (K.add(x, y), K.mul(x, y)) == (G.add(x, y), G.mul(x, y))
