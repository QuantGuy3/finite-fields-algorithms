"""Pruebas del módulo nativo ``alcp_native`` (binding pybind11 de python-bindings/bindings.cpp).

Recorren toda la superficie expuesta (docs/diseno.md, §5):

* conversión de enteros Python <-> mpz en ambos sentidos (pequeños, negativos,
  de más de 64 bits, de más de 4300 dígitos, objetos con ``__index__``, bool) y
  rechazo de ``float``;
* traducción de excepciones del núcleo (ZeroDivisionError, ValueError,
  OverflowError) y errores de tipo;
* cuerpos ``Fp`` y ``Fq`` (igualdad, hash, operaciones de elemento, índices);
* polinomios ``PolyFp`` y ``PolyFq`` (aritmética, factorización y demás);
* algoritmos rápidos (DFT, Toeplitz, división) para ``Fp`` y ``Fq``, contra
  versiones ingenuas escritas aquí;
* concordancia con la práctica original ``python/cuerpos_finitos.py`` (solo
  multiplicación, divmod y gcd en F_p con p impar, donde la referencia es
  correcta).

El módulo se importa como ``alcp.alcp_native`` (pyproject.toml añade
python-bindings a sys.path), igual que lo hará el paquete ``alcp``.
"""
from __future__ import annotations

import copy
import fractions
import gc
import itertools
import random
import sys

import pytest

from alcp import alcp_native as al

P_NTT = 998244353            # 2^23 divide a p - 1
M61 = 2**61 - 1              # primo de Mersenne de 61 bits
M127 = 2**127 - 1
M521 = 2**521 - 1


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------
class Indice:
    """Entero «de mentira» con __index__ (como numpy.int64)."""

    def __init__(self, v):
        self.v = v

    def __index__(self):
        return self.v


class IndiceRoto:
    """__index__ que lanza: la excepción debe propagarse tal cual."""

    def __index__(self):
        raise RuntimeError("__index__ roto")


class MiEntero(int):
    """Subclase de int (no es un int exacto)."""


def recortar(c):
    """Quita los ceros finales de una lista de coeficientes."""
    c = list(c)
    while c and c[-1] == 0:
        c.pop()
    return c


def poli_aleatorio(K, rng, grado, rango):
    """Polinomio de grado exacto `grado` con coeficientes en range(rango)."""
    return K.poly([rng.randrange(rango) for _ in range(grado)] + [rng.randrange(1, rango)])


@pytest.fixture(autouse=True)
def semilla():
    """Generador del núcleo con semilla fija: pruebas reproducibles."""
    al.set_seed(20261017)


@pytest.fixture
def rng():
    """Generador de Python con semilla fija para los datos de prueba."""
    return random.Random(7)


@pytest.fixture(scope="module")
def ref(practicas):
    """Práctica original (python/cuerpos_finitos.py), solo como referencia.

    La importa la fixture compartida «practicas» de conftest.py, sin escribir .pyc
    para no dejar rastro en python/.
    """
    return practicas.cf


@pytest.fixture(params=["Fp", "Fq"])
def K(request):
    """Los dos cuerpos de prueba habituales: F_101 y F_9 = F_3[a]/(a^2 + 1)."""
    return al.Fp(101) if request.param == "Fp" else al.Fq(3, [1, 0, 1])


# ---------------------------------------------------------------------------
# Conversión de enteros
# ---------------------------------------------------------------------------
VALORES = [
    0, 1, -1, 7, -7,
    2**31 - 1, 2**31, -2**31, -2**31 - 1,      # frontera de long en Windows (32 bits)
    2**32 + 5, -(2**40) - 3,
    2**62, 2**63 - 1, -(2**63) + 1, -(2**63),  # frontera de 64 bits
    2**63, -(2**63) - 1, 2**64, -(2**64),
    2**100 + 12345, -(2**100) - 12345,
    M521, -M521,
    10**5000 + 3, -(10**5000) - 3,              # más de 4300 dígitos
]


def id_entero(x):
    """Identificador corto (pytest no puede imprimir enteros de más de 4300 dígitos)."""
    return f"{'-' if x < 0 else '+'}{abs(x).bit_length()}bits"


@pytest.mark.parametrize("x", VALORES, ids=id_entero)
def test_entero_ida_y_vuelta(x):
    y = al._int_roundtrip(x)
    assert y == x and type(y) is int


def test_entero_sin_limite_de_digitos():
    # Con el límite de conversión str<->int al mínimo (640 dígitos), un camino
    # que pasara por cadenas decimales fallaría; el binding usa bytes.
    antes = sys.get_int_max_str_digits()
    sys.set_int_max_str_digits(640)
    try:
        x = 7**3000 + 1  # ~2536 dígitos
        assert al._int_roundtrip(x) == x
        assert al._int_roundtrip(-x) == -x
        F = al.Fp(M521)
        assert F.add(x, 0) == x % M521
        assert F.poly([x, -x]).coeffs() == [x % M521, (-x) % M521]
    finally:
        sys.set_int_max_str_digits(antes)


@pytest.mark.parametrize("obj, esperado", [
    (True, 1), (False, 0),
    (MiEntero(5), 5), (MiEntero(-2**70), -2**70),
    (Indice(9), 9), (Indice(-3), -3), (Indice(2**80 + 1), 2**80 + 1),
])
def test_entero_desde_index(obj, esperado):
    assert al._int_roundtrip(obj) == esperado
    F = al.Fp(101)
    assert F.add(obj, 0) == esperado % 101
    assert F.poly([obj, obj]).coeffs() == recortar([esperado % 101] * 2)


def test_entero_numpy():
    np = pytest.importorskip("numpy")
    F = al.Fp(101)
    assert F.mul(np.int64(7), np.int32(3)) == 21
    assert F.poly(np.array([1, 2, 3], dtype=np.int64)).coeffs() == [1, 2, 3]
    with pytest.raises(TypeError):
        F.add(np.float64(1.0), 1)


@pytest.mark.parametrize("malo", [1.0, 2.5, float("nan")])
def test_float_rechazado(malo):
    F = al.Fp(7)
    with pytest.raises(TypeError, match="float"):
        F.add(malo, 1)
    with pytest.raises(TypeError, match="float"):
        F.poly([1, malo])
    with pytest.raises(TypeError, match="float"):
        al._int_roundtrip(malo)


@pytest.mark.parametrize("malo", ["5", fractions.Fraction(1, 2), None, [1], 1j])
def test_no_enteros_rechazados(malo):
    with pytest.raises(TypeError):
        al._int_roundtrip(malo)
    with pytest.raises(TypeError):
        al.Fp(7).add(malo, 1)


def test_index_que_lanza_se_propaga():
    with pytest.raises(RuntimeError, match="__index__ roto"):
        al.Fp(7).add(IndiceRoto(), 1)


def test_lista_de_coeficientes_no_valida():
    F = al.Fp(7)
    with pytest.raises(TypeError):
        F.poly("123")
    with pytest.raises(TypeError):
        F.poly(5)
    assert F.poly((1, 2)).coeffs() == [1, 2]  # cualquier secuencia sirve


# ---------------------------------------------------------------------------
# Enteros: primalidad, potencias de primo y semilla
# ---------------------------------------------------------------------------
def test_is_prime():
    assert al.is_prime(2) and al.is_prime(P_NTT) and al.is_prime(M127)
    assert not al.is_prime(P_NTT * 7)
    assert not al.is_prime(1) and not al.is_prime(0) and not al.is_prime(-7)
    assert al.is_prime(Indice(13)) and not al.is_prime(True)


def test_prime_power_factor():
    assert al.prime_power_factor(3**10) == (3, 10)
    assert al.prime_power_factor(1024) == (2, 10)
    assert al.prime_power_factor(101) == (101, 1)
    assert al.prime_power_factor(M61**3) == (M61, 3)
    for malo in (1, 0, -8, 12, 2**10 * 3):
        with pytest.raises(ValueError):
            al.prime_power_factor(malo)


@pytest.mark.parametrize("s", [0, 5, 2**31 + 1, 2**64 + 3, 10**40])
def test_set_seed_reproducible(s):
    F = al.Fp(M127)
    al.set_seed(s)
    a = [F.random() for _ in range(5)]
    g = al.Fq(3**5).modulus
    al.set_seed(s)
    assert [F.random() for _ in range(5)] == a
    assert al.Fq(3**5).modulus == g
    assert all(0 <= x < M127 for x in a)


def test_set_seed_negativa():
    with pytest.raises(ValueError, match="semilla"):
        al.set_seed(-1)


# ---------------------------------------------------------------------------
# Cuerpo primo Fp
# ---------------------------------------------------------------------------
def test_fp_parametros_y_repr():
    F = al.Fp(7)
    assert (F.p, F.q, F.n) == (7, 7, 1)
    assert repr(F) == "Fp(7)"
    G = al.Fp(M521)
    assert G.p == M521 and repr(G) == f"Fp({M521})"


@pytest.mark.parametrize("malo", [0, 1, 4, 9, -7, 2**64])
def test_fp_no_primo(malo):
    with pytest.raises(ValueError):
        al.Fp(malo)


def test_fp_igualdad_hash_y_copia():
    F, G, H = al.Fp(7), al.Fp(7), al.Fp(11)
    assert F == G and not (F != G) and hash(F) == hash(G)
    assert F != H and F != 7 and F != al.Fq(7)
    assert F.__eq__(7) is NotImplemented
    assert len({F, G, H}) == 2
    assert copy.copy(F) is F and copy.deepcopy(F) is F


def test_fp_operaciones():
    p = 101
    F = al.Fp(p)
    for a, b in itertools.product([0, 1, 5, 100, -3, 250, 10**30], repeat=2):
        assert F.add(a, b) == (a + b) % p
        assert F.sub(a, b) == (a - b) % p
        assert F.mul(a, b) == (a * b) % p
        if b % p:
            assert F.div(a, b) == (a * pow(b, -1, p)) % p
        assert F.neg(a) == (-a) % p
        assert F.from_int(a) == a % p
    for a in range(1, p):
        assert F.mul(a, F.inv(a)) == 1
        assert F.pow(a, p - 1) == 1
        assert F.pow(a, -2) == pow(a, -2, p)
    assert F.pow(0, 0) == 1 and F.pow(5, 0) == 1 and F.pow(3, 2**100) == pow(3, 2**100, p)


def test_fp_division_por_cero():
    F = al.Fp(7)
    for f in (lambda: F.inv(0), lambda: F.inv(14), lambda: F.div(3, 0), lambda: F.pow(0, -1)):
        with pytest.raises(ZeroDivisionError):
            f()


def test_fp_aleatorio():
    F = al.Fp(7)
    vistos = {F.random() for _ in range(300)}
    assert vistos == set(range(7))
    for d in (0, 1, 5):
        f = F.random_poly(d)
        assert isinstance(f, al.PolyFp) and f.degree() == d and f.field is F
    assert F.random_poly(-1).is_zero() and F.random_poly(-5).is_zero()
    for enorme in (2**31 - 1, 2**40):  # d + 1 debe caber en un int del núcleo
        with pytest.raises(OverflowError):
            F.random_poly(enorme)


@pytest.mark.parametrize("p, cap", [(2, 0), (3, 1), (17, 4), (101, 2), (P_NTT, 23)])
def test_fp_fft(p, cap):
    F = al.Fp(p)
    assert F.fft_capacity() == cap
    r = F.fft_base_root()
    assert r == F.fft_base_root()  # fijada en el cuerpo
    assert F.pow(r, 2**cap) == 1
    if cap:
        assert F.pow(r, 2 ** (cap - 1)) != 1


# ---------------------------------------------------------------------------
# Cuerpo de extensión Fq (elementos = índices en base p)
# ---------------------------------------------------------------------------
def test_fq_parametros_y_repr():
    K = al.Fq(3, [1, 0, 1])
    assert (K.p, K.n, K.q, K.modulus) == (3, 2, 9, [1, 0, 1])
    assert repr(K) == "Fq(3, [1, 0, 1])"
    # El módulo se hace mónico: 2x^2 + 2 -> x^2 + 1.
    assert al.Fq(3, [2, 0, 2]).modulus == [1, 0, 1]
    for q, (p, n) in [(9, (3, 2)), (8, (2, 3)), (7, (7, 1)), (5**4, (5, 4))]:
        L = al.Fq(q)
        assert (L.p, L.n, L.q) == (p, n, q)
        assert len(L.modulus) == n + 1 and L.modulus[-1] == 1
        assert al.Fp(p).poly(L.modulus).is_irreducible()
    assert eval(repr(K), {"Fq": al.Fq}) == K


@pytest.mark.parametrize("args", [(6,), (1,), (0,), (-9,), (12,), (3, [1]), (3, []), (3, [2, 0, 1]),
                                  (4, [1, 1, 1]), (3, [1, 0, 3])])
def test_fq_invalido(args):
    with pytest.raises(ValueError):
        al.Fq(*args)


def test_fq_igualdad_y_hash():
    K1, K2 = al.Fq(3, [1, 0, 1]), al.Fq(3, [1, 0, 1])
    K3 = al.Fq(3, [2, 2, 1])  # x^2 + 2x + 2, otro módulo: cuerpo isomorfo, presentación distinta
    assert K1 == K2 and hash(K1) == hash(K2)
    assert K1 != K3 and K1 != al.Fp(3) and K1 != "F9"
    assert K1.__eq__(al.Fp(3)) is NotImplemented
    assert len({K1, K2, K3}) == 2
    assert copy.copy(K1) is K1 and copy.deepcopy(K1) is K1


def test_fq_indices_y_coeficientes():
    K = al.Fq(3, [1, 0, 1])
    for i in range(9):
        c = K.to_coeffs(i)
        assert len(c) == 2 and c == [i % 3, i // 3]
        assert K.from_coeffs(c) == i
        assert K.from_index(i) == i
    assert K.from_index(9) == 0 and K.from_index(-1) == 8 and K.from_index(10**50) == 10**50 % 9
    assert K.to_coeffs(-1) == [2, 2]
    # from_coeffs es generoso: reduce módulo p y módulo g = a^2 + 1.
    assert K.from_coeffs([]) == 0
    assert K.from_coeffs([4]) == 1
    assert K.from_coeffs([-1, 5]) == 2 + 2 * 3
    assert K.from_coeffs([0, 0, 1]) == 2       # a^2 = -1
    assert K.from_coeffs([1, 2, 5, 7]) == 5    # 1 + 2a + 2a^2 + a^3 = 2 + a
    # from_int: inmersión escalar n·1 (no el índice).
    assert [K.from_int(n) for n in (0, 1, 2, 3, 4, -1)] == [0, 1, 2, 0, 1, 2]


def mul_fq_ingenua(K, i, j):
    """Producto en F_p[a]/(g) a partir de los coeficientes (referencia)."""
    p, n, g = K.p, K.n, K.modulus
    a, b = K.to_coeffs(i), K.to_coeffs(j)
    c = [0] * (2 * n - 1)
    for s, x in enumerate(a):
        for t, y in enumerate(b):
            c[s + t] += x * y
    for k in range(len(c) - 1, n - 1, -1):  # a^n = -(g_0 + ... + g_{n-1} a^{n-1})
        for t in range(n):
            c[k - n + t] -= c[k] * g[t]
        c[k] = 0
    return K.from_coeffs([x % p for x in c[:n]])


def test_fq_operaciones():
    K = al.Fq(3, [1, 0, 1])
    for i, j in itertools.product(range(9), repeat=2):
        ci, cj = K.to_coeffs(i), K.to_coeffs(j)
        assert K.add(i, j) == K.from_coeffs([x + y for x, y in zip(ci, cj)])
        assert K.sub(i, j) == K.from_coeffs([x - y for x, y in zip(ci, cj)])
        assert K.mul(i, j) == mul_fq_ingenua(K, i, j)
        if j:
            assert K.mul(K.div(i, j), j) == i
    for i in range(9):
        assert K.add(i, K.neg(i)) == 0
        if i:
            assert K.mul(i, K.inv(i)) == 1
            assert K.pow(i, 8) == 1
            assert K.pow(i, -1) == K.inv(i)
    # índices fuera de [0, q) se reducen módulo q
    assert K.add(9 + 4, -9 + 1) == K.add(4, 1)


def test_fq_division_por_cero():
    K = al.Fq(3, [1, 0, 1])
    for f in (lambda: K.inv(0), lambda: K.inv(9), lambda: K.div(4, 0), lambda: K.pow(0, -3)):
        with pytest.raises(ZeroDivisionError):
            f()


def test_fq_aleatorio_y_fft():
    K = al.Fq(3, [1, 0, 1])
    assert {K.random() for _ in range(400)} == set(range(9))
    f = K.random_poly(4)
    assert isinstance(f, al.PolyFq) and f.degree() == 4
    assert K.fft_capacity() == 3  # 8 = 2^3 divide a q - 1
    r = K.fft_base_root()
    assert K.pow(r, 8) == 1 and K.pow(r, 4) != 1
    assert al.Fq(2, [1, 1, 1]).fft_capacity() == 0  # q - 1 = 3


# ---------------------------------------------------------------------------
# Polinomios: construcción, acceso, igualdad, hash y copias
# ---------------------------------------------------------------------------
def test_poly_construccion(K):
    clase = al.PolyFp if isinstance(K, al.Fp) else al.PolyFq
    f = K.poly([1, 2, 0, 0])
    assert isinstance(f, clase)
    assert f.coeffs() == [1, 2] and f.degree() == 1
    assert clase(K, [1, 2]) == f and clase(K).is_zero() and clase(K, []).degree() == -1
    assert K.poly([]).coeffs() == [] and K.poly([0, 0]).is_zero()
    q = K.q
    assert K.poly([q + 1, -1]).coeffs() == ([1, 100] if isinstance(K, al.Fp) else [1, 8])
    with pytest.raises(TypeError):
        clase(al.Fp(7) if isinstance(K, al.Fq) else al.Fq(9), [1])


def test_poly_accesores(K):
    f = K.poly([3, 0, 5])
    assert f.coeff(0) == 3 and f.coeff(1) == 0 and f.coeff(2) == 5 and f.coeff(10**6) == 0
    assert f.leading() == 5
    assert not f.is_zero() and not f.is_one() and K.poly([1]).is_one() and not K.poly([]).is_one()
    with pytest.raises(ValueError):
        f.coeff(-1)
    with pytest.raises(OverflowError):
        f.coeff(2**70)
    with pytest.raises(ValueError):
        K.poly([]).leading()


def test_poly_field_mismo_objeto(K):
    f = K.poly([1, 1])
    assert f.field is K
    assert (f * f).field is K and f.factor()[0][0][0].field is K


def test_poly_field_sobrevive_al_cuerpo():
    f = al.Fp(7).poly([1, 2, 3])  # el objeto Fp(7) muere tras esta línea
    gc.collect()
    assert f.field == al.Fp(7) and f.field.p == 7
    assert (f * f).coeffs() == [1, 4, 3, 5, 2]
    g = al.Fq(3, [1, 0, 1]).poly([1, 4])
    gc.collect()
    assert g.field == al.Fq(3, [1, 0, 1]) and g.leading() == 4


def test_poly_resultado_reanclado_al_cuerpo_del_primer_operando():
    # gcd(x^2, x) devuelve internamente el segundo operando, que apunta a OTRO
    # objeto cuerpo (igual en parámetros).  El binding lo reancla al del primer
    # operando; tras destruir el segundo cuerpo el resultado sigue siendo válido.
    F1, F2 = al.Fp(7), al.Fp(7)
    a, b = F1.poly([0, 0, 1]), F2.poly([0, 1])
    g = a.gcd(b)
    del F2, b
    gc.collect()
    assert g.field is F1 and g.coeffs() == [0, 1] and (g * g).coeffs() == [0, 0, 1]
    K1, K2 = al.Fq(3, [1, 0, 1]), al.Fq(3, [1, 0, 1])
    h = K1.poly([0, 0, 1]).gcd(K2.poly([0, 1]))
    del K2
    gc.collect()
    assert h.field is K1 and h.coeffs() == [0, 1]


def test_poly_igualdad_y_hash(K):
    f = K.poly([1, 2, 3])
    g = K.poly([1, 2, 3, 0])
    assert f == g and not (f != g) and hash(f) == hash(g)
    assert f != K.poly([1, 2]) and f != 5 and f != "x" and f.__eq__(5) is NotImplemented
    # igualdad entre cuerpos distintos pero con los mismos parámetros
    K2 = al.Fp(101) if isinstance(K, al.Fp) else al.Fq(3, [1, 0, 1])
    h = K2.poly([1, 2, 3])
    assert f == h and hash(f) == hash(h) and len({f, g, h}) == 1
    # cuerpos distintos: distintos (sin excepción)
    otro = al.Fp(103).poly([1, 2, 3]) if isinstance(K, al.Fp) else al.Fq(3, [2, 2, 1]).poly([1, 2, 3])
    assert f != otro
    assert al.Fp(7).poly([1]) != al.Fq(7).poly([1])  # PolyFp frente a PolyFq
    assert {K.poly([]): "cero"}[K.poly([0])] == "cero"


def test_poly_copias(K):
    f = K.poly([1, 2, 3])
    for g in (copy.copy(f), copy.deepcopy(f), copy.deepcopy([f])[0]):
        assert g == f and g is not f and g.field is K


def test_poly_repr_y_to_string():
    F = al.Fp(7)
    f = F.poly([1, 0, 3])
    assert repr(f) == "PolyFp((1) + (3)*x^2)"
    assert f.to_string() == "(1) + (3)*x^2" and f.to_string("y") == "(1) + (3)*y^2"
    assert F.poly([]).to_string() == "0"
    K = al.Fq(3, [1, 0, 1])
    assert repr(K.poly([K.from_coeffs([1, 2]), 1])) == "PolyFq((1+2a) + (1)*x)"


# ---------------------------------------------------------------------------
# Polinomios: aritmética
# ---------------------------------------------------------------------------
def test_poly_aritmetica_basica():
    F = al.Fp(7)
    f, g = F.poly([1, 0, 1]), F.poly([2, 1])
    assert (f + g).coeffs() == [3, 1, 1]
    assert (f - g).coeffs() == [6, 6, 1]
    assert (-g).coeffs() == [5, 6]
    assert (f * g).coeffs() == [2, 1, 2, 1]
    assert (f ** 3).coeffs() == [1, 0, 3, 0, 3, 0, 1] and (f ** 0).is_one()
    q, r = (f * g + F.poly([3])).divmod(g)
    assert q == f and r.coeffs() == [3]
    assert divmod(f * g, g) == (f, F.poly([]))
    assert (f // g) * g + (f % g) == f
    assert f.scalar_mul(3).coeffs() == [3, 0, 3] and f.scalar_mul(7).is_zero()
    assert F.poly([5, 3, 1, 4]).derivative().coeffs() == [3, 2, 5]
    assert F.poly([5]).derivative().is_zero()
    assert F.poly([2, 4]).monic().coeffs() == [4, 1] and F.poly([]).monic().is_zero()


def test_poly_aritmetica_fq():
    K = al.Fq(3, [1, 0, 1])
    a = K.from_coeffs([0, 1])  # la clase de a, con a^2 = -1
    f = K.poly([0, a])          # a·x
    assert (f * f).coeffs() == [0, 0, K.neg(1)]
    assert f.scalar_mul(a).coeffs() == [0, K.neg(1)]
    assert f.evaluate(a) == K.neg(1)
    assert (f + f + f).is_zero()  # característica 3
    assert f.monic().coeffs() == [0, 1]


def test_poly_identidades_algebraicas(K, rng):
    q = K.q
    for _ in range(10):
        f = poli_aleatorio(K, rng, rng.randint(0, 6), q)
        g = poli_aleatorio(K, rng, rng.randint(0, 6), q)
        assert (f - g) + g == f and (-f) + f == K.poly([]) and -(-f) == f
        assert (f * g).derivative() == f.derivative() * g + f * g.derivative()  # regla de Leibniz
        assert f.scalar_mul(1) == f and f.scalar_mul(0).is_zero()
        q_, r_ = divmod(f, g)
        assert q_ * g + r_ == f and (r_.is_zero() or r_.degree() < g.degree())
        assert f.to_string("t").count("t") >= min(f.degree(), 1)


def test_poly_grados_errores(K):
    f = K.poly([1, 1]) * K.poly([2, 1])
    for malo in (lambda: K.poly([]).distinct_degree(), lambda: K.poly([3]).distinct_degree(),
                 lambda: f.equal_degree(0), lambda: f.equal_degree(3),
                 lambda: K.poly([]).equal_degree(1), lambda: K.poly([3]).equal_degree(1)):
        with pytest.raises(ValueError):
            malo()


def test_poly_operandos_de_otro_tipo(K):
    f = K.poly([1, 1])
    for malo in (1, "x", 1.5, None):
        for op in (lambda: f + malo, lambda: f - malo, lambda: f * malo,
                   lambda: f // malo, lambda: f % malo, lambda: divmod(f, malo)):
            with pytest.raises(TypeError):
                op()
    with pytest.raises(TypeError):
        f ** 1.5
    with pytest.raises(TypeError):
        pow(f, 2, f)


def test_poly_cuerpos_distintos():
    f = al.Fp(7).poly([1, 1])
    g = al.Fp(11).poly([1, 1])
    for op in (lambda: f + g, lambda: f - g, lambda: f * g, lambda: f // g, lambda: f % g,
               lambda: divmod(f, g), lambda: f.divmod(g), lambda: f.mul(g, "school"),
               lambda: f.gcd(g), lambda: f.gcd_ext(g), lambda: f.inv_mod(g),
               lambda: f.pow_mod(2, g), lambda: f.multiplicity(g), lambda: f.compose(g),
               lambda: al.divmod_toeplitz(f, g)):
        with pytest.raises(ValueError, match="cuerpos distintos"):
            op()
    K1, K2 = al.Fq(3, [1, 0, 1]), al.Fq(3, [2, 2, 1])
    with pytest.raises(ValueError, match="cuerpos distintos"):
        K1.poly([1, 1]) * K2.poly([1, 1])
    # mismos parámetros, objetos distintos: se opera con normalidad
    assert (f + al.Fp(7).poly([1])).coeffs() == [2, 1]


def test_poly_division_por_cero(K):
    f, cero = K.poly([1, 2, 3]), K.poly([])
    for op in (lambda: f // cero, lambda: f % cero, lambda: divmod(f, cero), lambda: f.divmod(cero),
               lambda: f.pow_mod(3, cero), lambda: al.divmod_toeplitz(f, cero)):
        with pytest.raises(ZeroDivisionError):
            op()


def test_poly_potencia_negativa(K):
    with pytest.raises(ValueError):
        K.poly([1, 1]) ** -1


@pytest.mark.parametrize("algo", ["auto", "school", "schoolbook", "karatsuba", "fft"])
def test_poly_mul_algoritmos(algo, rng):
    F = al.Fp(P_NTT)  # admite FFT de cualquier tamaño de prueba
    a = poli_aleatorio(F, rng, 80, P_NTT)
    b = poli_aleatorio(F, rng, 70, P_NTT)
    assert a.mul(b, algo) == a.mul(b) == a * b == a.mul(b, algo="school")
    G = al.Fp(101)  # FFT imposible a este tamaño: 'fft' recurre a Karatsuba
    c, d = poli_aleatorio(G, rng, 40, 101), poli_aleatorio(G, rng, 35, 101)
    assert c.mul(d, algo) == c.mul(d, "school")
    K = al.Fq(3, [1, 0, 1])
    e, f = poli_aleatorio(K, rng, 3, 9), poli_aleatorio(K, rng, 2, 9)
    assert e.mul(f, algo) == e.mul(f, "school")


def test_poly_mul_algoritmo_desconocido(K):
    f = K.poly([1, 1])
    for malo in ("", "Karatsuba", "ntt", "rapido"):
        with pytest.raises(ValueError, match="desconocido"):
            f.mul(f, malo)


# ---------------------------------------------------------------------------
# Polinomios: evaluación, composición, gcd, potencias modulares
# ---------------------------------------------------------------------------
def test_poly_evaluate_y_compose(K, rng):
    q = K.q
    f = poli_aleatorio(K, rng, 4, q)
    g = poli_aleatorio(K, rng, 2, q)
    h = f.compose(g)
    assert h.degree() == 8
    for x in range(min(q, 30)):
        # Horner a mano con las operaciones de elemento del cuerpo
        acc = 0
        for c in reversed(f.coeffs()):
            acc = K.add(K.mul(acc, x), c)
        assert f.evaluate(x) == acc
        assert h.evaluate(x) == f.evaluate(g.evaluate(x))
    assert K.poly([]).evaluate(3) == 0
    assert f.evaluate(q + 2) == f.evaluate(2)  # x se reduce como cualquier elemento


def test_poly_gcd_bezout_inverso(K, rng):
    q = K.q
    for _ in range(20):
        a = poli_aleatorio(K, rng, rng.randint(1, 6), q)
        b = poli_aleatorio(K, rng, rng.randint(1, 4), q)
        g, x, y = a.gcd_ext(b)
        assert g == x * a + y * b == a.gcd(b)
        assert g.is_zero() or g.leading() == 1
    cero = K.poly([])
    assert cero.gcd(cero).is_zero() and all(t.is_zero() for t in cero.gcd_ext(cero))
    m = K.poly([1, 0, 1, 1])  # x^3 + x^2 + 1, coprimo con x + 1 en F_101 y en F_9
    a = K.poly([1, 1])
    assert a.gcd(m).is_one()
    assert (a * a.inv_mod(m) % m).is_one()
    with pytest.raises(ValueError):
        K.poly([0, 1]).inv_mod(K.poly([0, 0, 1]))  # x no es invertible módulo x^2


def test_poly_pow_mod(K, rng):
    q = K.q
    m = poli_aleatorio(K, rng, 3, q)
    a = poli_aleatorio(K, rng, 5, q)
    for k in (0, 1, 2, 7, 20):
        assert a.pow_mod(k, m) == (a ** k) % m
    while a.gcd(m).degree() != 0:
        a = a + K.poly([1])
    assert (a.pow_mod(-3, m) * a.pow_mod(3, m) % m).is_one()
    # exponente de más de 64 bits: (a^(2^40))^(2^40) = a^(2^80)
    assert a.pow_mod(2**80, m) == a.pow_mod(2**40, m).pow_mod(2**40, m)


# ---------------------------------------------------------------------------
# Polinomios: irreducibilidad y factorización
# ---------------------------------------------------------------------------
def test_poly_irreducible_y_libre_de_cuadrados():
    F = al.Fp(7)
    assert F.poly([1, 0, 1]).is_irreducible()       # -1 no es cuadrado mod 7
    assert not F.poly([6, 0, 1]).is_irreducible()   # (x - 1)(x + 1)
    assert not F.poly([3]).is_irreducible() and not F.poly([]).is_irreducible()
    assert F.poly([6, 0, 1]).is_squarefree() and F.poly([3]).is_squarefree()
    assert not (F.poly([6, 1]) ** 2).is_squarefree() and not F.poly([]).is_squarefree()
    assert not F.poly([1] + [0] * 6 + [1]).is_squarefree()  # x^7 + 1 = (x + 1)^7


def test_poly_multiplicidad_y_radical(K):
    u, v = K.poly([1, 1]), K.poly([1, 0, 1])  # x + 1 y x^2 + 1, coprimos en F_101 y en F_9
    f = u ** 3 * v
    assert f.multiplicity(u) == 3 and f.multiplicity(K.poly([5, 1])) == 0
    r = f.radical()
    assert r.is_squarefree() and r == (u * v).monic()
    assert K.poly([]).radical().is_zero() and K.poly([4]).radical().is_one()
    with pytest.raises(ValueError):
        f.multiplicity(K.poly([3]))
    with pytest.raises(ValueError):
        K.poly([]).multiplicity(u)
    # rama f' = 0 de radical: f = g(x^p)
    p = K.p
    g = K.poly([1] + [0] * (p - 1) + [1])  # x^p + 1 = (x + 1)^p
    assert g.radical() == K.poly([1, 1])


def test_poly_grados_distintos_e_iguales():
    F = al.Fp(7)
    lin = [F.poly([c, 1]) for c in (1, 2, 3)]
    cuad = [F.poly([1, 0, 1]), F.poly([3, 1, 1])]
    assert all(c.is_irreducible() for c in cuad)
    f = lin[0] * lin[1] * lin[2] * cuad[0] * cuad[1]
    dd = f.distinct_degree()
    assert [d for d, _ in dd] == [1, 2]
    assert dd[0][1] == (lin[0] * lin[1] * lin[2]).monic() and dd[1][1] == (cuad[0] * cuad[1]).monic()
    assert isinstance(dd[0][1], al.PolyFp) and dd[0][1].field is F
    lineales = dd[0][1].equal_degree(1)
    assert sorted(t.coeffs() for t in lineales) == sorted(t.coeffs() for t in lin)
    cuadraticos = dd[1][1].equal_degree(2)
    assert sorted(t.coeffs() for t in cuadraticos) == sorted(t.coeffs() for t in cuad)
    for malo in (lambda: F.poly([]).distinct_degree(), lambda: F.poly([3]).distinct_degree(),
                 lambda: f.equal_degree(0), lambda: f.equal_degree(4), lambda: F.poly([3]).equal_degree(1),
                 lambda: F.poly([]).equal_degree(1)):
        with pytest.raises(ValueError):
            malo()
    with pytest.raises(OverflowError):
        f.equal_degree(2**40)


def test_poly_grados_fq():
    K = al.Fq(3, [1, 0, 1])
    f = K.poly([1, 1]) * K.poly([2, 1]) * K.poly([K.from_coeffs([0, 1]), 1])
    dd = f.distinct_degree()
    assert [d for d, _ in dd] == [1] and dd[0][1] == f.monic()
    assert len(dd[0][1].equal_degree(1)) == 3


def reconstruir(K, factores, lider):
    r = K.poly([lider])
    for u, e in factores:
        r = r * u ** e
    return r


def test_factorizacion_fp(rng):
    F = al.Fp(101)
    for _ in range(150):
        f = F.poly([1])
        for _ in range(rng.randint(1, 3)):
            base = poli_aleatorio(F, rng, rng.randint(1, 3), 101).monic()
            f = f * base ** rng.randint(1, 3)
        f = f.scalar_mul(rng.randrange(1, 101))
        facs, lider = f.factor()
        assert lider == f.leading()
        assert all(u.is_irreducible() and u.leading() == 1 and e >= 1 for u, e in facs)
        assert reconstruir(F, facs, lider) == f
        claves = [(u.degree(), list(reversed(u.coeffs()))) for u, _ in facs]
        assert claves == sorted(claves)  # orden canónico


def test_factorizacion_fq(rng):
    K = al.Fq(3, [1, 0, 1])
    for _ in range(80):
        f = K.poly([1])
        for _ in range(rng.randint(1, 2)):
            f = f * poli_aleatorio(K, rng, rng.randint(1, 2), 9).monic() ** rng.randint(1, 2)
        facs, lider = f.factor()
        assert all(u.is_irreducible() for u, _ in facs)
        assert reconstruir(K, facs, lider) == f


def test_factorizacion_casos_limite():
    F2 = al.Fp(2)
    f = F2.poly([0, 0, 0, 0, 0, 0, 1, 1, 1, 1])  # x^6 (x + 1)^3: la práctica original fallaba
    facs, lider = f.factor()
    assert [(u.coeffs(), e) for u, e in facs] == [([0, 1], 6), ([1, 1], 3)] and lider == 1
    assert al.Fp(7).poly([5]).factor() == ([], 5)
    with pytest.raises(ValueError):
        al.Fp(7).poly([]).factor()
    K = al.Fq(9)
    facs, lider = K.poly([3, 0, 1]).factor()
    assert reconstruir(K, facs, lider) == K.poly([3, 0, 1])


# ---------------------------------------------------------------------------
# Algoritmos rápidos (Fp y Fq)
# ---------------------------------------------------------------------------
def dft_ingenua(K, raiz, a):
    n = len(a)
    return [
        _suma(K, [K.mul(a[i], K.pow(raiz, i * j)) for i in range(n)]) for j in range(n)
    ]


def _suma(K, xs):
    acc = 0
    for x in xs:
        acc = K.add(acc, x)
    return acc


def raiz_de_orden(K, k):
    """Raíz primitiva 2^k-ésima de la unidad a partir de la raíz base del cuerpo."""
    return K.pow(K.fft_base_root(), 2 ** (K.fft_capacity() - k))


@pytest.mark.parametrize("K, k", [(al.Fp(17), 4), (al.Fp(P_NTT), 6), (al.Fq(3, [1, 0, 1]), 3),
                                  (al.Fq(5, [2, 0, 1]), 3)])
def test_dft_idft(K, k, rng):
    raiz = raiz_de_orden(K, k)
    a = [rng.randrange(K.q) for _ in range(2**k)]
    A = al.dft(K, raiz, k, a)
    assert A == dft_ingenua(K, raiz, a)
    assert al.idft(K, raiz, k, A) == a
    assert al.dft(K, 1, 0, [5]) == [5] and al.idft(K, 1, 0, [5]) == [5]


def test_dft_compatibilidad_api_original():
    F17 = al.Fp(17)
    a = list(range(16))
    assert al.idft(F17, 3, 4, al.dft(F17, 3, 4, a)) == a  # 3 tiene orden 16 en F_17*


def test_dft_errores():
    F = al.Fp(17)
    with pytest.raises(ValueError):
        al.dft(F, 3, 4, [1] * 15)         # longitud distinta de 2^k
    with pytest.raises(ValueError):
        al.dft(F, 9, 4, [1] * 16)         # 9 tiene orden 8, no 16
    with pytest.raises(ValueError):
        al.idft(F, 16, 4, [1] * 16)       # 16 = -1 tiene orden 2
    with pytest.raises(ValueError):
        al.dft(F, 3, -1, [1])             # k negativo
    with pytest.raises(ValueError):
        al.dft(F, 3, 70, [1])             # k demasiado grande para el núcleo
    with pytest.raises(OverflowError):
        al.dft(F, 3, 2**80, [1])          # k no cabe en un entero de máquina
    with pytest.raises(TypeError):
        al.dft("F17", 3, 4, [1] * 16)     # ninguna sobrecarga admite el primer argumento
    K = al.Fq(3, [1, 0, 1])
    with pytest.raises(ValueError):
        al.dft(K, 1, 2, [1] * 4)          # 1 no tiene orden 4


def toeplitz_ingenua(K, M, b):
    return [_suma(K, [K.mul(M[i][j], b[j]) for j in range(len(b))]) for i in range(len(M))]


def matriz_inferior(a, n):
    return [[a[i - j] if i >= j else 0 for j in range(n)] for i in range(n)]


def matriz_superior(a, n):
    return [[a[j - i] if j >= i else 0 for j in range(n)] for i in range(n)]


def matriz_completa(a, n):
    # a = primera fila (a[0..n-1]) seguida de la primera columna sin la esquina (a[n..2n-2])
    return [[a[j - i] if j >= i else a[n - 1 + i - j] for j in range(n)] for i in range(n)]


def producto(K, A, B):
    n = len(A)
    return [[_suma(K, [K.mul(A[i][t], B[t][j]) for t in range(n)]) for j in range(n)] for i in range(n)]


@pytest.mark.parametrize("n", [1, 2, 5, 8])
def test_toeplitz_productos(K, n, rng):
    q = K.q
    a = [rng.randrange(q) for _ in range(n + 2)]  # a más largo que n: solo se usan n
    b = [rng.randrange(q) for _ in range(n)]
    assert al.toeplitz_lower_vec(K, n, a, b) == toeplitz_ingenua(K, matriz_inferior(a, n), b)
    assert al.toeplitz_upper_vec(K, n, a, b) == toeplitz_ingenua(K, matriz_superior(a, n), b)
    c = [rng.randrange(q) for _ in range(2 * n - 1)]
    assert al.toeplitz_vec(K, n, c, b) == toeplitz_ingenua(K, matriz_completa(c, n), b)


@pytest.mark.parametrize("n", [1, 2, 7, 8])
def test_toeplitz_inversas(K, n, rng):
    q = K.q
    a = [rng.randrange(1, q)] + [rng.randrange(q) for _ in range(n - 1)]
    identidad = [[int(i == j) for j in range(n)] for i in range(n)]
    x = al.toeplitz_lower_inverse(K, n, a)
    assert producto(K, matriz_inferior(a, n), matriz_inferior(x, n)) == identidad
    assert al.toeplitz_lower_vec(K, n, a, x) == [1] + [0] * (n - 1)
    y = al.toeplitz_upper_inverse(K, n, a)
    assert producto(K, matriz_superior(a, n), matriz_superior(y, n)) == identidad


def test_toeplitz_errores(K):
    for f in (lambda: al.toeplitz_lower_vec(K, 0, [], []),
              lambda: al.toeplitz_lower_vec(K, 3, [1, 2], [1, 2, 3]),
              lambda: al.toeplitz_upper_vec(K, 3, [1, 2, 3], [1, 2]),
              lambda: al.toeplitz_vec(K, 2, [1, 2], [1, 2]),
              lambda: al.toeplitz_lower_inverse(K, 0, [1]),
              lambda: al.toeplitz_upper_inverse(K, 3, [1, 2]),
              lambda: al.toeplitz_lower_vec(K, -1, [1], [1])):
        with pytest.raises(ValueError):
            f()
    with pytest.raises(OverflowError):
        al.toeplitz_vec(K, 2**70, [1], [1])
    for f in (lambda: al.toeplitz_lower_inverse(K, 3, [0, 1, 2]),
              lambda: al.toeplitz_upper_inverse(K, 1, [K.q])):
        with pytest.raises(ZeroDivisionError):
            f()


def test_toeplitz_compatibilidad_api_original(rng):
    F = al.Fp(101)
    nn = 8
    col = [rng.randrange(101) for _ in range(nn)]
    vec = [rng.randrange(101) for _ in range(nn)]
    y = al.toeplitz_lower_vec(F, nn, col, vec)
    assert y == [sum(col[i - j] * vec[j] for j in range(i + 1)) % 101 for i in range(nn)]
    col[0] = col[0] or 1
    inv = al.toeplitz_lower_inverse(F, nn, col)
    assert al.toeplitz_lower_vec(F, nn, col, inv) == [1] + [0] * (nn - 1)


def test_divmod_toeplitz(K, rng):
    q = K.q
    for df, dg in [(6, 2), (5, 5), (2, 4), (0, 0), (9, 1)]:
        f = poli_aleatorio(K, rng, df, q)
        g = poli_aleatorio(K, rng, dg, q)
        qt, rt = al.divmod_toeplitz(f, g)
        assert (qt, rt) == f.divmod(g)
        assert qt.field is K and rt.field is K
    cero = K.poly([])
    assert al.divmod_toeplitz(cero, K.poly([1, 1])) == (cero, cero)


# ---------------------------------------------------------------------------
# Concordancia con la práctica original (F_p con p impar)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("p", [7, 101])
def test_concordancia_con_la_practica_original(ref, p, rng):
    fp = ref.cuerpo_fp(p)
    fpx = ref.anillo_fp_x(fp)
    F = al.Fp(p)

    def a_ref(c):
        return fpx.elem_de_tuple(tuple(c))

    def de_ref(poli):
        return recortar(int(c) for c in fpx.conv_a_tuple(poli))

    for _ in range(25):
        fI = [rng.randrange(p) for _ in range(rng.randint(0, 7))] + [rng.randrange(1, p)]
        gI = [rng.randrange(p) for _ in range(rng.randint(0, 4))] + [rng.randrange(1, p)]
        f, g = F.poly(fI), F.poly(gI)
        assert (f * g).coeffs() == de_ref(fpx.mult(a_ref(fI), a_ref(gI)))
        qr, rr = fpx.divmod(a_ref(fI), a_ref(gI))
        q, r = f.divmod(g)
        assert (q.coeffs(), r.coeffs()) == (de_ref(qr), de_ref(rr))
        assert f.gcd(g).coeffs() == de_ref(fpx.gcd(a_ref(fI), a_ref(gI)))
    # Los elementos de la práctica (ElementoFp) tienen __index__: se aceptan tal cual.
    h = fpx.elem_de_tuple((1, 2, 3))
    assert F.poly(fpx.conv_a_tuple(h)).coeffs() == [1, 2, 3]


# ---------------------------------------------------------------------------
# Documentación
# ---------------------------------------------------------------------------
def test_docstrings_en_todo_lo_expuesto():
    assert "Cuerpos finitos" in al.__doc__
    sin_doc = []
    for nombre in dir(al):
        if nombre.startswith("__"):
            continue
        obj = getattr(al, nombre)
        if isinstance(obj, type):
            if not (obj.__doc__ or "").strip():
                sin_doc.append(nombre)
            for meto, valor in vars(obj).items():
                # Métodos especiales propios; se ignoran los atributos internos de
                # Python y de pybind11 (__module__, _pybind11_conduit_v1_, ...).
                if meto.startswith("_") and meto not in {
                    "__init__", "__eq__", "__hash__", "__repr__", "__copy__", "__deepcopy__",
                    "__add__", "__sub__", "__neg__", "__mul__", "__floordiv__", "__mod__",
                    "__divmod__", "__pow__"}:
                    continue
                doc = getattr(valor, "__doc__", None) or ""
                # pybind11 antepone la firma: debe haber texto propio además de ella.
                lineas = [ln for ln in doc.splitlines() if ln.strip()]
                if len(lineas) < 2 and not isinstance(valor, property):
                    sin_doc.append(f"{nombre}.{meto}")
                if isinstance(valor, property) and not (valor.__doc__ or "").strip():
                    sin_doc.append(f"{nombre}.{meto}")
        elif callable(obj):
            lineas = [ln for ln in (obj.__doc__ or "").splitlines() if ln.strip()]
            if len(lineas) < 2:
                sin_doc.append(nombre)
    assert sin_doc == []
