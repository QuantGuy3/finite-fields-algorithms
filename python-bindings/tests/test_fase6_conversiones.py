"""Fase 6, propuesta 8: conversión directa ``int`` <-> ``u64`` en el binding.

``bindings.cpp`` convierte ahora sin pasar por ``mpz`` los enteros de Python que
son ``int`` exactos y caben en 64 bits con signo (``K.poly(lista)``,
``PolyFp(K, lista)``, ``coeffs()``, ``coeff``, ``leading``, ``evaluate``,
``scalar_mul`` y las operaciones de elementos), y deja el camino de siempre para
todo lo demás.  Estas pruebas comprueban que nada observable cambia:

* los valores coinciden con la referencia en Python (``x % p``, ``x % q``) y con
  el backend GMP (``backend="gmp"``, que no tiene camino de 64 bits en el
  núcleo), en varios cuerpos y en los bordes de 2^31, 2^32, 2^62, 2^63 y 2^64;
* ``Fq.from_int`` sigue siendo la inmersión ``(a mod p)·1`` y no el índice
  ``a mod q`` (error real encontrado al integrar la propuesta);
* los tipos y mensajes de excepción son los de antes, incluido el orden en que
  se detecta el primer elemento erróneo de una lista;
* bool, subclases de int, ``numpy.int64``, objetos con ``__index__`` y
  secuencias que no son list/tuple siguen el camino de siempre;
* el generador aleatorio se consume igual (cuerpos recién creados, misma
  semilla, mismas extracciones posteriores).
"""
from __future__ import annotations

import random
import re

import pytest

from alcp import alcp_native as al

pytestmark = pytest.mark.nativo

P_NTT = 998244353
M61 = 2**61 - 1
P62 = 2**62 - 57             # mayor primo por debajo de 2^62 (último de Fp64)
M127 = 2**127 - 1            # backend GMP obligatorio

MSG_FLOAT = "se esperaba un entero (int o un objeto con __index__), no un valor de tipo 'float'"
MSG_ARGS = "incompatible (function|constructor) arguments"
RE_FLOAT = "^" + re.escape(MSG_FLOAT) + "$"

# Bordes del camino rápido (long long) y del cambio de mpz_set_si a mpz_import
# (long de 32 bits en Windows).
BORDES = [0, 1, -1, 2, -2, 6, 7, 8, -7, -8, 100, -100,
          2**31 - 1, 2**31, -2**31, -2**31 - 1, 2**32 + 5, -(2**32 + 5),
          M61, M61 + 1, -M61, P62, P62 + 1, -P62, 2**62, -2**62,
          2**63 - 1, -2**63, -2**63 + 1, 2**63, -2**63 - 1, 2**64 - 1, 2**64, -2**64,
          M127, -(2**200) + 3]


class Indice:
    def __init__(self, v):
        self.v = v

    def __index__(self):
        return self.v


class IndiceRoto:
    def __index__(self):
        raise RuntimeError("__index__ roto")


class MiEntero(int):
    pass


def cuerpos_primos():
    """(cuerpo, p) en Fp64 (varios tamaños) y en GMP (forzado y obligatorio)."""
    return [(al.Fp(p), p) for p in (2, 7, 101, P_NTT, M61, P62)] + \
        [(al.Fp(7, backend="gmp"), 7), (al.Fp(P_NTT, backend="gmp"), P_NTT), (al.Fp(M127), M127)]


def recortar(c):
    while c and c[-1] == 0:
        c = c[:-1]
    return c


# ---------------------------------------------------------------------------
# Valores
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("K, p", cuerpos_primos(), ids=lambda x: repr(x))
def test_elementos_contra_referencia(K, p):
    for a in BORDES:
        assert K.from_int(a) == a % p
        assert K.neg(a) == (-a) % p
        for b in (3, -5, 2**63 - 1, -2**63, 2**64 + 1, a):
            assert K.add(a, b) == (a + b) % p
            assert K.sub(a, b) == (a - b) % p
            assert K.mul(a, b) == (a * b) % p
            if b % p:
                assert K.div(a, b) == a * pow(b, -1, p) % p
        if a % p:
            assert K.inv(a) == pow(a, -1, p)
            assert K.pow(a, -3) == pow(a, -3, p)
        assert K.pow(a, 5) == pow(a, 5, p)
    # Las salidas son int exactos (no subclases) en [0, p).
    assert type(K.mul(-1, 1)) is int
    assert type(K.add(2**64, 0)) is int


@pytest.mark.parametrize("K, p", cuerpos_primos(), ids=lambda x: repr(x))
def test_listas_contra_referencia(K, p):
    rnd = random.Random(p % 1000)
    for d in (0, 1, 2, 31, 100, 1000):
        for rango in (p, 2**31, 2**63, 2**80):
            c = [rnd.randrange(-rango, rango) for _ in range(d + 1)]
            esperado = recortar([x % p for x in c])
            f = K.poly(c)
            assert f.coeffs() == esperado
            assert K.poly(tuple(c)).coeffs() == esperado
            assert al.PolyFp(K, c) == f
            assert all(type(x) is int for x in f.coeffs())
            if esperado:
                assert f.leading() == esperado[-1]
                assert f.coeff(0) == esperado[0]
            assert f.coeff(d + 5) == 0
            x = rnd.randrange(-2**70, 2**70)
            assert f.evaluate(x) == sum(ci * pow(x, i, p) for i, ci in enumerate(esperado)) % p
            assert f.scalar_mul(-3).coeffs() == recortar([(-3 * ci) % p for ci in esperado])
    # Bordes en una sola lista (camino rápido y lista con un entero grande).
    cortos = [b for b in BORDES if -2**63 <= b < 2**63]
    assert K.poly(cortos + [1]).coeffs() == [b % p for b in cortos] + [1 % p]
    assert K.poly(BORDES + [1]).coeffs() == [b % p for b in BORDES] + [1 % p]
    assert K.poly([]).coeffs() == []
    assert al.PolyFp(K).coeffs() == []


def test_fp64_igual_que_gmp():
    """Mismo cuerpo con los dos backends: mismos coeficientes, factores y hash."""
    rnd = random.Random(8)
    for p in (7, 101, P_NTT, M61):
        A, B = al.Fp(p), al.Fp(p, backend="gmp")
        for d in (3, 20, 200):
            c = [rnd.randrange(-2**63, 2**63) for _ in range(d + 1)]
            f, g = A.poly(c), B.poly(c)
            assert f.coeffs() == g.coeffs()
            assert hash(f) == hash(g)
            assert f.evaluate(-12345) == g.evaluate(-12345)
        c = [rnd.randrange(-50, 50) for _ in range(13)] + [1]
        fa, fb = A.poly(c).factor(), B.poly(c).factor()
        assert [(u.coeffs(), e) for u, e in fa[0]] == [(u.coeffs(), e) for u, e in fb[0]]
        assert fa[1] == fb[1]


@pytest.mark.parametrize("K", [al.Fq(49), al.Fq(7, [3, 1, 1]), al.Fq(49, backend="gmp"), al.Fq(2**20)],
                         ids=repr)
def test_extensiones(K):
    """En F_q las entradas son índices módulo q y from_int es (a mod p)·1."""
    p, q = K.p, K.q
    for a in BORDES:
        assert K.from_int(a) == a % p          # inmersión, NO a % q
        assert K.from_index(a) == a % q
        assert K.add(a, 1) == K.add(a % q, 1)
        assert K.neg(a) == K.neg(a % q)
        assert K.mul(a, 3) == K.mul(a % q, 3)
    c = [b for b in BORDES if -2**63 <= b < 2**63]
    assert K.poly(c).coeffs() == recortar([b % q for b in c])
    assert K.poly(BORDES).coeffs() == recortar([b % q for b in BORDES])
    assert al.PolyFq(K, c) == K.poly(c)
    f = K.poly([1, 2, 3])
    assert f.evaluate(-1) == f.evaluate(q - 1)
    assert f.scalar_mul(-1) == f.scalar_mul(q - 1)


# ---------------------------------------------------------------------------
# Tipos que no son int exacto: camino de siempre
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("obj, valor", [
    (True, 1), (False, 0), (MiEntero(5), 5), (MiEntero(-2**70), -2**70),
    (Indice(9), 9), (Indice(-2**65), -2**65),
])
def test_no_exactos(obj, valor):
    for K, p in cuerpos_primos():
        assert K.add(obj, 0) == valor % p
        assert K.from_int(obj) == valor % p
        assert K.poly([obj, 1]).coeffs() == [valor % p, 1 % p]
        assert K.poly([1, obj]).coeffs() == recortar([1 % p, valor % p])


def test_numpy():
    np = pytest.importorskip("numpy")
    for K, p in cuerpos_primos():
        for v in (np.int64(-5), np.int64(2**63 - 1), np.int32(-3), np.uint64(2**64 - 1), np.int8(-128)):
            assert K.add(v, 0) == int(v) % p
            assert K.poly([v, 1]).coeffs() == [int(v) % p, 1 % p]
        arr = np.array([-1, 2, -2**63], dtype=np.int64)
        assert K.poly(arr).coeffs() == recortar([int(x) % p for x in arr])
        assert K.poly(list(arr)).coeffs() == recortar([int(x) % p for x in arr])
    F = al.Fp(7)
    with pytest.raises(TypeError, match="float64"):
        F.add(np.float64(1.0), 1)
    with pytest.raises(TypeError, match=MSG_ARGS):
        F.add(np.bool_(True), 1)


def test_secuencias_que_no_son_lista():
    F = al.Fp(101)
    assert F.poly(range(-3, 3)).coeffs() == [x % 101 for x in range(-3, 3)]
    assert F.poly(x for x in (1, 2, 3)).coeffs() == [1, 2, 3]
    # pybind11 3 acepta cualquier iterable (como antes de la fase 6): un set también.
    assert F.poly({5}).coeffs() == [5]
    assert F.poly([Indice(2**80), -1]).coeffs() == [2**80 % 101, 100]


# ---------------------------------------------------------------------------
# Errores: mismos tipos, mensajes y orden de detección
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("K, p", cuerpos_primos(), ids=lambda x: repr(x))
def test_errores_de_elementos(K, p):
    for malo in (1.0, float("nan")):
        with pytest.raises(TypeError, match=RE_FLOAT):
            K.add(malo, 1)
        with pytest.raises(TypeError, match=RE_FLOAT):
            K.mul(1, malo)
        with pytest.raises(TypeError, match=RE_FLOAT):
            K.poly([1, 2]).evaluate(malo)
    for malo in ("3", None, b"1", [1], 1j):
        with pytest.raises(TypeError, match=MSG_ARGS):
            K.add(malo, 1)
        with pytest.raises(TypeError, match=MSG_ARGS):
            K.from_int(malo)
        with pytest.raises(TypeError, match=MSG_ARGS):
            K.poly([1]).scalar_mul(malo)
    # El primer argumento se convierte antes que el segundo.
    with pytest.raises(TypeError, match=RE_FLOAT):
        K.add(2.0, "a")
    with pytest.raises(TypeError, match=MSG_ARGS):
        K.add("a", 2.0)
    with pytest.raises(RuntimeError, match="__index__ roto"):
        K.add(IndiceRoto(), 1)
    with pytest.raises(ZeroDivisionError):
        K.inv(p)
    with pytest.raises(ZeroDivisionError):
        K.div(1, -p)


@pytest.mark.parametrize("K", [al.Fp(7), al.Fp(7, backend="gmp"), al.Fp(M127), al.Fq(49)], ids=repr)
def test_errores_de_listas(K):
    crear = al.PolyFp if K.n == 1 else al.PolyFq
    for fn in (K.poly, lambda c: crear(K, c)):
        for mala in ([1, None], [1, "a"], ["a", 2.0], "123", b"12", 5, None, [[1]], [1, 1j]):
            with pytest.raises(TypeError, match=MSG_ARGS):
                fn(mala)
        # El primer elemento erróneo decide: float antes que str -> mensaje del float.
        for mala in ([1, 2.0], [2.0, "a"], (1, 2**70, 2.5), [2**64, 0.5]):
            with pytest.raises(TypeError, match=RE_FLOAT):
                fn(mala)
        with pytest.raises(RuntimeError, match="__index__ roto"):
            fn([1, IndiceRoto()])
    with pytest.raises(TypeError, match=MSG_ARGS):
        al.PolyFp(al.Fq(49), [1])


# ---------------------------------------------------------------------------
# Generador aleatorio: mismo consumo con y sin conversiones
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("p", [7, P_NTT, M61, M127])
def test_generador_en_frio(p):
    def tras(semilla, con_conversiones):
        K = al.Fp(p)                      # cuerpo recién creado: caché FFT vacía
        al.set_seed(semilla)
        if con_conversiones:
            c = list(range(-600, 600)) + [2**70, -2**63]
            f = K.poly(c)
            assert f.coeffs() == recortar([x % p for x in c])
            K.mul(-1, 2**63 - 1), K.add(5, -5), K.from_int(-3), f.evaluate(-2), f.leading()
            K.poly(tuple(c)).scalar_mul(-1).coeffs()
        return [K.random() for _ in range(8)], K.random_poly(20).coeffs(), K.fft_base_root()

    for semilla in (1, 12345, 2**40 + 7):
        assert tras(semilla, True) == tras(semilla, False)


def test_generador_extension_en_frio():
    def tras(con_conversiones):
        K = al.Fq(7, [3, 1, 1])
        al.set_seed(99)
        if con_conversiones:
            K.poly([-1, 2**63 - 1, 5]).coeffs(), K.from_int(-8), K.mul(-1, 2)
        return [K.random() for _ in range(5)], K.random_poly(10).coeffs()

    assert tras(True) == tras(False)
