"""Fase 7, propuesta 3: backend de dos palabras (FpN2, ``backend == 'u128'``).

``Fp(p)`` con ``2^62 <= p < 2^128`` usa FpN2 (alcp/fpn.hpp).  Desde Python debe
ser indistinguible del backend GMP (``backend="gmp"``) salvo en la propiedad
``backend``: mismos valores, mensajes, repr, igualdad entre cuerpos y mezclas,
y el mismo consumo del generador aleatorio con la misma semilla.
"""
import copy

import pytest

from alcp import alcp_native as al
from test_native import _bateria

pytestmark = pytest.mark.nativo

M127 = 2**127 - 1
P63 = 2**62 + 135             # el menor primo de 63 bits (primero de FpN2)
P64 = 2**64 - 59              # el mayor primo de 64 bits
P65 = 2**64 + 13
P_NTT128 = 81 * 2**121 + 1    # 128 bits con capacidad FFT 121
P_NTT64 = 2**64 - 2**32 + 1   # 64 bits con capacidad FFT 32
P_MAX = 2**128 - 159          # el mayor primo de 128 bits
PRIMOS = [P63, P64, P65, M127, P_NTT128, P_NTT64, P_MAX]


class Indice:
    def __init__(self, v):
        self.v = v

    def __index__(self):
        return self.v


@pytest.mark.parametrize("p", PRIMOS)
def test_rango_y_backend(p):
    assert al.is_prime(p)
    assert al.Fp(p).backend == "u128" and al.Fp(p, backend="gmp").backend == "gmp"
    assert al.Fp(2**62 - 57).backend == "u64" and al.Fp(2**128 + 51).backend == "gmp"
    assert al.Fq(p, [1, 1]).backend == "gmp"  # Fq no tiene backend de dos palabras


@pytest.mark.parametrize("p", PRIMOS)
def test_misma_bateria_que_gmp(p):
    al.set_seed(23)
    rapido = al.Fp(p)
    al.set_seed(23)
    gmp = al.Fp(p, backend="gmp")
    assert _bateria(rapido, 9) == _bateria(gmp, 9)


@pytest.mark.parametrize("p", PRIMOS)
def test_elementos_y_conversiones(p):
    A, B = al.Fp(p), al.Fp(p, backend="gmp")
    valores = [0, 1, -1, 2, -2, p - 1, p, p + 1, -p, -p - 1, 2**62, 2**63 - 1, -2**63, 2**63, 2**64 - 1, 2**64,
               -2**64, 2**127, 2**128 - 1, 2**128, -(2**200) + 3, 3**300, True, Indice(-7), Indice(2**130)]
    for x in valores:
        assert A.from_int(x) == B.from_int(x)
        for y in (3, p - 2, -5, 2**100 + 1):
            assert A.add(x, y) == B.add(x, y) and A.sub(x, y) == B.sub(x, y)
            assert A.mul(x, y) == B.mul(x, y) and A.div(x, y) == B.div(x, y)
        assert A.neg(x) == B.neg(x) and A.pow(x, 5) == B.pow(x, 5) and A.pow(x, p + 3) == B.pow(x, p + 3)
        if B.from_int(x) != 0:
            assert A.inv(x) == B.inv(x) and A.pow(x, -3) == B.pow(x, -3)
    for K in (A, B):
        with pytest.raises(ZeroDivisionError, match="0 is not invertible"):
            K.inv(p)
        with pytest.raises(ZeroDivisionError):
            K.pow(0, -1)
        with pytest.raises(TypeError):
            K.add(1.5, 1)
    # Listas de coeficientes cortas y largas (con enteros grandes).
    cs = [5, -3, 2**64 + 7, p - 1, -(2**90), 0, 1]
    assert A.poly(cs).coeffs() == B.poly(cs).coeffs() and al.PolyFp(A, cs).coeffs() == B.poly(cs).coeffs()
    assert A.poly([1, 2, 3]).coeffs() == [1, 2, 3]
    assert A.poly(cs).coeff(2) == B.poly(cs).coeff(2) and A.poly(cs).leading() == B.poly(cs).leading()


def test_no_primo_mismo_mensaje():
    n = 2**100 + 1
    for backend in ("auto", "gmp"):
        with pytest.raises(ValueError, match=f"^{n} is not prime; cannot build F_p$"):
            al.Fp(n, backend=backend)


def test_igualdad_hash_repr_y_mezclas():
    a, b = al.Fp(M127), al.Fp(M127, backend="gmp")
    assert a != b and a == al.Fp(M127) and b == al.Fp(M127, backend="gmp") and hash(a) == hash(b)
    assert repr(a) == f"Fp({M127})" and repr(b) == f"Fp({M127}, backend='gmp')"
    assert eval(repr(a), {"Fp": al.Fp}) == a and eval(repr(b), {"Fp": al.Fp}) == b
    assert repr(al.Fp(2**128 + 51)) == f"Fp({2**128 + 51})"
    assert a != al.Fp(P_MAX) and a != al.Fp(7) and copy.copy(a) is a and copy.deepcopy(a) is a
    f, g = a.poly([1, 2]), b.poly([1, 2])
    assert f != g and hash(f) == hash(g) and f == a.poly([1, 2]) and f.field == a
    for x, y in ((f, g), (g, f), (f, al.Fp(P_MAX).poly([1, 2])), (f, al.Fp(7).poly([1, 2]))):
        with pytest.raises(ValueError, match="cuerpos distintos"):
            x + y
        with pytest.raises(ValueError, match="cuerpos distintos"):
            x.gcd(y)
    with pytest.raises(TypeError):
        al.PolyFp(al.Fq(49), [1])
    assert repr(f) == repr(g) and str(f) == str(g)


@pytest.mark.parametrize("p", [M127, P_NTT128, P63])
def test_factor_y_generador(p):
    """factor, distinct_degree, equal_degree y random_poly con la misma semilla."""
    def correr(K):
        al.set_seed(77)
        h = K.random_poly(25)
        r = [h.coeffs(), K.random()]
        fac, lc = h.factor()
        r += [[(u.coeffs(), e) for u, e in fac], lc, h.is_irreducible()]
        rad = h.radical().monic()
        dd = rad.distinct_degree()
        r += [[(d, u.coeffs()) for d, u in dd]]
        r += [[[v.coeffs() for v in u.equal_degree(d)] for d, u in dd]]
        grande = K.random_poly(300) * K.random_poly(200)
        r += [grande.coeffs(), (grande % K.random_poly(150)).coeffs(), K.fft_base_root()]
        r.append(al.Fp(5).random())  # estado del generador al final
        return r
    assert correr(al.Fp(p)) == correr(al.Fp(p, backend="gmp"))


def test_generador_en_frio_con_fft():
    """Un cuerpo recién creado pide la raíz FFT en el mismo producto que GMP."""
    def tras(backend, na, nb):
        al.set_seed(5)
        K = al.Fp(P_NTT128, backend=backend)
        x = K.poly(list(range(1, na + 1))) * K.poly(list(range(2, nb + 2)))
        return x.coeffs(), al.Fp(7).random(), al.Fp(7).random()
    for na, nb in ((20, 20), (40, 40), (5, 200), (20, 200), (100, 100)):
        assert tras("auto", na, nb) == tras("gmp", na, nb)
