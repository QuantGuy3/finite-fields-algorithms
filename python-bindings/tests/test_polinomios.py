"""Pruebas de anillos de polinomios, polinomios y factorizaciones (alcp.poly).

La aritmética sobre F_p se compara con la referencia en Python puro
(tests/referencia.py); las raíces y la irreducibilidad, con fuerza bruta y con
la fórmula de Gauss para el número de mónicos irreducibles.  Se recorren las
reglas de conversión, las promociones F_p[x] -> F_q[x], los errores y la
igualdad/hash/pickle.
"""
from __future__ import annotations

import copy
import pickle
import random
from fractions import Fraction

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

import alcp
from alcp import Factorization, Polynomial, PolynomialRing
from referencia import irreducibles_gauss, monicos, pol_divmod, pol_eval, pol_mul, pol_resta, pol_suma

GF = alcp.GF


class Indice:
    def __init__(self, v):
        self.v = v

    def __index__(self):
        return self.v


@pytest.fixture(autouse=True)
def semilla():
    alcp.set_seed(424242)
    yield
    alcp.set_print_style("unicode")


def coefs(f):
    """Coeficientes enteros (valores en F_p o índices en F_q) de un polinomio."""
    return [c.index for c in f.coeffs]


# ---------------------------------------------------------------------------
# Anillo
# ---------------------------------------------------------------------------
def test_anillo_cacheado():
    F = GF(7)
    R = F["x"]
    assert PolynomialRing(F) is R and PolynomialRing(F, "x") is R and F.poly_ring() is R
    assert F["y"] is not R and F["y"] == F.poly_ring("y")
    assert R.field is F and R.base_field is F and R.base_ring is F and R.var == "x"
    assert R.gen == R.x and coefs(R.gen) == [0, 1] and R.zero == 0 and R.one == 1
    assert F["t"].t == F["t"].gen
    with pytest.raises(AttributeError, match="'PolynomialRing' no tiene el atributo 'y'"):
        R.y


@pytest.mark.parametrize("cuerpo, var, error, mensaje", [
    (7, "x", TypeError, "el cuerpo base debe ser un cuerpo finito"),
    (GF(7), 1, TypeError, "debe ser una cadena"),
    (GF(7), "x1", ValueError, "solo por letras"),
    (GF(9), "a", ValueError, "la variable 'a' ya es el generador de GF\\(9"),
])
def test_anillo_invalido(cuerpo, var, error, mensaje):
    with pytest.raises(error, match=mensaje):
        PolynomialRing(cuerpo, var)


def test_anillo_igualdad_hash_copia_pickle():
    R, S = GF(7)["x"], GF(9)["x"]
    for A in (R, S, GF(7)["y"]):
        assert copy.copy(A) is A and copy.deepcopy(A) is A
        assert pickle.loads(pickle.dumps(A)) is A
    assert R != S and R != GF(7)["y"] and R.__eq__(7) is NotImplemented and R != "x"
    assert len({R, S, GF(7)["x"], GF(7)["y"]}) == 3 and hash(R) == hash(GF(7)["x"])


def test_monomios():
    R = GF(9)["x"]
    assert R.monomial(3) == "x^3" and R.monomial(0, "a") == "a" and R.monomial(2, 5) == "2x^2"
    assert R.monomial(Indice(1)) == R.x
    with pytest.raises(ValueError, match="el grado de un monomio debe ser >= 0"):
        R.monomial(-1)
    with pytest.raises(ValueError, match="demasiado grande"):
        R.monomial(10**8)
    with pytest.raises(TypeError, match="el grado debe ser un entero"):
        R.monomial(2.0)


def test_aleatorios_e_irreducibles():
    R = GF(9)["x"]
    for d in range(6):
        assert R.random(d).degree == d
    assert R.random(-1) == 0
    alcp.set_seed(3)
    f = R.random(8)
    alcp.set_seed(3)
    assert R.random(8) == f
    for d in (1, 2, 3, 5):
        g = R.irreducible(d)
        assert g.degree == d and g.is_monic() and g.is_irreducible()
    with pytest.raises(ValueError, match="tiene grado >= 1, no 0"):
        R.irreducible(0)
    with pytest.raises(ValueError, match="demasiado grande"):
        R.random(10**8)


# ---------------------------------------------------------------------------
# Conversión
# ---------------------------------------------------------------------------
def test_conversion_constantes():
    F, K = GF(7), GF(9)
    R, S = F["x"], K["x"]
    assert coefs(R()) == [] and coefs(R(10)) == [3] and coefs(R(True)) == [1]
    assert coefs(R(Fraction(1, 2))) == [4] and coefs(R(Indice(-1))) == [6]
    assert coefs(R(F(5))) == [5] and coefs(S(K.gen)) == [3] and coefs(S(GF(3)(2))) == [2]
    assert coefs(R(GF(49)(3))) == [3]
    with pytest.raises(ValueError, match="no pertenece al subcuerpo primo"):
        GF(3)["x"](K.gen)


def test_conversion_cadenas_listas_dicts():
    F, K = GF(7), GF(9)
    R, S = F["x"], K["x"]
    assert coefs(R("x^3 + 2x - 1")) == [6, 2, 0, 1]
    assert coefs(R([1, 2, 0, 0])) == [1, 2] and coefs(R((0, "3", Fraction(1, 2)))) == [0, 3, 4]
    assert coefs(R({5: 1, 0: -1, 2: 0})) == [6, 0, 0, 0, 0, 1] and R({}) == 0
    assert coefs(R(range(3))) == [0, 1, 2] and coefs(R(c for c in [1, 1])) == [1, 1]
    assert S("(a+1)x^2 + a") == S([K.gen, 0, K.gen + 1]) == S({2: "a+1", 0: "a"})
    assert coefs(S([[0, 1], 1])) == [3, 1]          # el coeficiente [0, 1] es a
    with pytest.raises(ValueError, match="exponente negativo en un polinomio: -1"):
        R({-1: 1})
    with pytest.raises(ValueError, match="demasiado grande"):
        R({10**9: 1})
    with pytest.raises(TypeError, match="el exponente debe ser un entero"):
        R({"2": 1})


def test_conversion_entre_anillos():
    F3, K = GF(3), GF(9)
    f = F3["x"]("x^2 + 2")
    assert F3["y"](f) == "y^2 + 2" and F3["y"](f).ring is F3["y"]
    assert K["x"](f).ring is K["x"] and K["x"](f) == "x^2 + 2"
    assert F3["x"](K["x"]("x + 2")) == "x + 2"
    assert F3["x"](f) is f
    with pytest.raises(ValueError, match="no pertenece al subcuerpo primo"):
        F3["x"](K["x"]("a x"))
    with pytest.raises(ValueError, match="cuerpos distintos"):
        GF(5)["x"](f)


@pytest.mark.parametrize("valor, error, mensaje", [
    (1.5, TypeError, "usa enteros o Fraction"),
    (None, TypeError, "no se puede convertir un valor de tipo 'NoneType' en un polinomio"),
    (object(), TypeError, "en un polinomio de GF\\(7\\)\\[x\\]"),
    ({1, 2}, TypeError, "un conjunto no tiene orden"),
    ([1.5], TypeError, "usa enteros o Fraction"),
    ("x + y", ValueError, "variable desconocida 'y'"),
])
def test_conversion_invalida(valor, error, mensaje):
    with pytest.raises(error, match=mensaje):
        GF(7)["x"](valor)


def test_conversion_y_operaciones_con_numpy():
    np = pytest.importorskip("numpy")
    R = GF(7)["x"]
    f = R("x + 1")
    assert R(np.int64(10)) == 3 and R(np.array([1, 2, 3])) == "3x^2 + 2x + 1"
    assert np.int64(2) * f == 2 * f and f + np.int32(1) == "x + 2" and np.int64(3) - f == "2 - x"
    assert f ** np.int64(2) == f * f and f[np.int64(1)] == 1
    with pytest.raises(TypeError, match="usa enteros o Fraction"):
        f * np.float64(2.0)


def test_constructor_polynomial():
    R = GF(7)["x"]
    assert Polynomial(R, "x+1") == "x + 1" and Polynomial(R) == 0
    with pytest.raises(TypeError, match="ring debe ser un PolynomialRing"):
        Polynomial(GF(7), [1])


def test_pertenencia_al_anillo():
    R = GF(7)["x"]
    assert "x^2" in R and 3 in R and [1, 2] in R and GF(7)(1) in R
    assert 1.5 not in R and "y" not in R and None not in R and GF(5)(1) not in R


# ---------------------------------------------------------------------------
# Aritmética contra la referencia
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("p", [2, 7, 13])
def test_aritmetica_contra_referencia(p):
    R = GF(p)["x"]
    rng = random.Random(p)
    for _ in range(60):
        a = [rng.randrange(p) for _ in range(rng.randrange(0, 9))]
        b = [rng.randrange(p) for _ in range(rng.randrange(0, 7))]
        f, g = R(a), R(b)
        assert coefs(f + g) == pol_suma(a, b, p)
        assert coefs(f - g) == pol_resta(a, b, p)
        assert coefs(f * g) == pol_mul(a, b, p)
        assert coefs(-f) == pol_resta([], a, p)
        if g:
            q, r = pol_divmod(a, b, p)
            assert coefs(f // g) == q and coefs(f % g) == r
            assert [coefs(h) for h in divmod(f, g)] == [q, r]
        for x in range(p):
            assert f(x) == pol_eval(a, x, p)


def test_operaciones_con_reflejo():
    K = GF(9)
    a = K.gen
    R = K["x"]
    f = R("x^2 + a")
    assert f + 1 == "x^2 + a + 1" and 1 + f == f + 1 and f - "x" == "x^2 - x + a"
    assert "x" - f == "-x^2 + x - a" and 2 * f == f * 2 and a * f == f * a
    assert f * "x + a" == "x^3 + a x^2 + a x + a^2" and [0, 1] * f == f * R.x
    assert Fraction(1, 2) * f == 2 * f and f // 2 == f * 2 and f % 2 == 0 and f / 2 == 2 * f
    assert 2 // f == 0 and 2 % f == 2 and divmod(f, "x") == (R.x, R(a))
    assert divmod("x^3", f) == (R.x, R("-a x")) and f / R("a") == f * a ** -1
    assert (f * R("x+1")) / "x + 1" == f and R(2) / R(2) == 1 and 2 / R(a) == 2 / a
    assert +f is f and f - f == 0


def test_division_exacta_y_por_cero():
    R = GF(7)["x"]
    f = R("x^2 + 1")
    with pytest.raises(ValueError, match="no es exacta: usa // y %"):
        f / "x + 1"
    with pytest.raises(ValueError, match="no es exacta"):
        2 / R("x")
    for op in (lambda: f // 0, lambda: f % R(0), lambda: divmod(f, 0), lambda: f / 0,
               lambda: 3 // R(0), lambda: divmod(3, R(0))):
        with pytest.raises(ZeroDivisionError, match="división por el polinomio cero"):
            op()


def test_operandos_no_admitidos():
    f = GF(7)["x"]("x")
    for malo in (None, object()):
        with pytest.raises(TypeError):
            f + malo
        with pytest.raises(TypeError):
            malo * f
    assert f.__add__(None) is NotImplemented and f.__rmod__(None) is NotImplemented
    with pytest.raises(TypeError, match="usa enteros o Fraction"):
        f * 1.5


def test_anillos_distintos_no_se_mezclan():
    f = GF(7)["x"]("x")
    for g in (GF(7)["y"]("y"), GF(5)["x"]("x"), GF(49)["y"]("y")):
        with pytest.raises(ValueError, match="son anillos distintos"):
            f + g
    with pytest.raises(ValueError, match="son anillos distintos"):
        GF(9)["x"]("x") * GF(9, var="b")["x"]("x")
    with pytest.raises(ValueError, match="cuerpos distintos"):
        f + GF(5)(1)


def test_promocion_de_fp_a_extension():
    F3, K = GF(3), GF(9)
    a = K.gen
    f, g = F3["x"]("x + 1"), K["x"]("x + a")
    for h in (f + g, g + f, f * g, g - f, f - g, f * a, a * f, f + a, a - f, f // g, g % f,
              f.gcd(g), f.compose(g)):
        assert h.ring is K["x"], h
    assert f + g == "2x + a + 1" and f - g == "1 - a" and a * f == "a x + a"
    assert f == K["x"]("x + 1") and K["x"]("x + 1") == f and hash(f) == hash(K["x"]("x+1"))
    assert g + F3(2) == "x + a + 2" and F3(2) * g == "2x + 2a"


def test_potencias():
    R = GF(7)["x"]
    f = R("x + 2")
    assert f ** 0 == 1 and f ** 3 == f * f * f and f ** Indice(2) == f * f
    assert R(3) ** -1 == 5 and R(3) ** -2 == R(5) ** 2 and R(0) ** 0 == 1
    with pytest.raises(ZeroDivisionError, match="el polinomio cero no tiene inverso"):
        R(0) ** -1
    with pytest.raises(ValueError, match="un polinomio no constante no tiene inverso"):
        f ** -1
    with pytest.raises(TypeError, match="el exponente debe ser un entero"):
        f ** 1.5
    with pytest.raises(TypeError, match="no un elemento de un cuerpo"):
        f ** GF(7)(2)
    with pytest.raises(ValueError, match="demasiado grande"):
        f ** (10**8)


def test_potencia_modular():
    R = GF(7)["x"]
    f, m = R("x + 2"), R("x^3 + x + 1")
    assert pow(f, 10, m) == (f ** 10) % m and f.pow_mod(10, "x^3 + x + 1") == pow(f, 10, m)
    assert (pow(f, -3, m) * f ** 3) % m == 1
    with pytest.raises(ValueError, match="no es invertible módulo"):
        pow(R("x"), -1, R("x^2"))
    with pytest.raises(ZeroDivisionError, match="potencia modular con módulo cero"):
        pow(f, 2, 0)
    with pytest.raises(TypeError, match="no se puede convertir un valor de tipo 'NoneType'"):
        f.pow_mod(2, None)


@pytest.mark.parametrize("q", [7, 9, 998244353])
def test_multiplicacion_por_algoritmos(q):
    R = GF(q)["x"]
    f, g = R.random(70), R.random(65)
    h = f * g
    for algo in ("auto", "school", "schoolbook", "karatsuba", "fft"):
        assert f.mul(g, algo) == h
    assert f.mul("x") == f * R.x
    with pytest.raises(ValueError, match="algoritmo de multiplicación desconocido: 'rapido'"):
        f.mul(g, "rapido")
    with pytest.raises(TypeError, match="algorithm debe ser una cadena"):
        f.mul(g, 1)


# ---------------------------------------------------------------------------
# mcd, inversos, composición y evaluación
# ---------------------------------------------------------------------------
def test_mcd_y_euclides_extendido():
    R = GF(7)["x"]
    f, g = R("(x-1)(x-2)^2"), R("3(x-1)(x-3)")
    assert f.gcd(g) == "x - 1" and f.gcd(0) == f.monic() and R(0).gcd(0) == 0
    for a, b in ((f, g), (R.random(9), R.random(5)), (R("x^2+1"), "x")):
        d, s, t = a.xgcd(b)
        assert d == s * a + t * R(b) and (d.is_monic() or d == 0)
    assert R(0).xgcd(0) == (0, 0, 0)


def test_inverso_modular():
    R = GF(9)["x"]
    m = R.irreducible(4)
    for _ in range(5):
        f = R.random(3)
        assert (f * f.inverse_mod(m)) % m == 1
    with pytest.raises(ValueError, match="no es invertible módulo"):
        R("x^2").inverse_mod("x^3")
    with pytest.raises(ZeroDivisionError, match="inverso módulo el polinomio cero"):
        R("x").inverse_mod(0)


def test_composicion_y_evaluacion():
    F, K = GF(7), GF(9)
    a = K.gen
    R, S = F["x"], K["x"]
    f = R("x^2 + 3x + 1")
    assert f.compose("x + 1") == f(R("x + 1")) == f("x + 1") == "x^2 + 5x + 5"
    assert f(2) == 4 + 6 + 1 and f("2") == f(2) and f(Fraction(1, 2)) == f(4) and f(F(2)) == f(2)
    assert f.evaluate(3) == f(3) and isinstance(f(3), alcp.PrimeFieldElement)
    assert f(F["y"]("y")) == F["y"]("y^2 + 3y + 1")          # cambio de variable
    g = GF(3)["x"]("x^2 + 1")
    assert g(a) == 0 and g(a).field is K                   # F_p[x] evaluado en F_q
    assert g(S("x + a")).ring is S and g(S("x + a")) == "x^2 + 2a x"
    h = S("a x + 1")
    assert h(GF(3)["t"]("t + 1")) == K["t"]("a t + a + 1")   # v sobre F_p se sube a K
    assert h("a") == a * a + 1 and h("x^2") == S("a x^2 + 1")
    with pytest.raises(ValueError, match="cuerpos distintos"):
        h(GF(5)["x"]("x"))
    with pytest.raises(ValueError, match="cuerpos distintos"):
        f(GF(5)(1))


def test_cambio_de_cuerpo():
    f = GF(3)["x"]("x^2 + 1")
    assert f.change_ring(GF(9)).ring is GF(9)["x"]
    with pytest.raises(TypeError, match="se esperaba un cuerpo finito"):
        f.change_ring(3)


# ---------------------------------------------------------------------------
# Datos y secuencia
# ---------------------------------------------------------------------------
def test_datos_del_polinomio():
    K = GF(9)
    a = K.gen
    R = K["x"]
    f = R("a x^3 + 2x + 1")
    assert f.ring is R and f.parent is R and f.field is K and f.base_field is K and f.var == "x"
    assert f.degree == 3 and R(0).degree == -1 and len(f) == 4 and len(R(0)) == 0
    assert f.coeffs == [1, 2, 0, a] and f.coefficients() == f.coeffs
    assert f.coefficients(sparse=True) == [1, 2, a] and f.to_dict() == {0: 1, 1: 2, 3: a}
    assert f.lc == a and f.leading_coefficient == a and R(0).lc == 0
    assert f[0] == 1 and f[3] == a and f[10] == 0 and f[10**30] == 0 and f[-1] == a and f[-4] == 1
    assert f[1:3] == [2, 0] and list(f) == f.coeffs and f[Indice(1)] == 2
    with pytest.raises(IndexError, match="fuera de rango"):
        f[-5]
    with pytest.raises(TypeError, match="el índice debe ser un entero"):
        f[1.0]
    assert bool(f) and not R(0) and R(0).is_zero() and R(1).is_one() and not f.is_one()
    assert f.monic().is_monic() and not f.is_monic() and not R(0).is_monic()
    assert R(3).is_constant() and R(0).is_constant() and not f.is_constant()


def test_derivada_radical_y_libres_de_cuadrados():
    R = GF(3)["x"]
    assert R("x^3 + 2x^2 + x").derivative() == "x + 1"
    assert R("x^6").derivative() == 0
    f = R("(x+1)^3 (x^2+1)^2 x")
    assert f.radical() == R("(x+1)(x^2+1)x") and R(0).radical() == 0 and R(2).radical() == 1
    assert not f.is_squarefree() and f.radical().is_squarefree()
    assert R(2).is_squarefree() and not R(0).is_squarefree()


def test_multiplicidad():
    R = GF(5)["x"]
    f = R("(x-1)^3 (x+1)")
    assert f.multiplicity("x - 1") == 3 and f.multiplicity("x+1") == 1 and f.multiplicity("x") == 0
    assert f.multiplicity(R("(x-1)^2")) == 1
    with pytest.raises(ValueError, match="divisores de grado >= 1"):
        f.multiplicity(2)
    with pytest.raises(ValueError, match="en el polinomio cero es infinita"):
        R(0).multiplicity("x")


# ---------------------------------------------------------------------------
# Factorización
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("q, n", [(2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (3, 2),
                                  (3, 3), (3, 4), (4, 2), (4, 3), (5, 2), (5, 3), (7, 2),
                                  (9, 2)])
def test_numero_de_irreducibles_formula_de_gauss(q, n):
    K = GF(q)
    R = K["x"]
    total = sum(1 for c in monicos(q, n) if R([K.from_index(i) for i in c]).is_irreducible())
    assert total == irreducibles_gauss(q, n)


@pytest.mark.parametrize("q", [2, 3, 4, 7, 9, 25])
def test_factorizacion_aleatoria(q):
    R = GF(q)["x"]
    for d in (1, 2, 5, 9, 12):
        f = R.random(d) * R.random(2) ** 2
        F = f.factor()
        assert F.expand() == f and F.unit == f.lc
        for u, e in F:
            assert u.is_monic() and u.is_irreducible() and e >= 1 and f.multiplicity(u) == e
        claves = [(u.degree, [tuple(c.coeffs) for c in reversed(u.coeffs)]) for u, _ in F]
        claves = [(d_, [[int(x) for x in c] for c in cs]) for d_, cs in claves]
        assert claves == sorted(claves)        # orden canónico


def test_factorizacion_casos():
    R = GF(7)["x"]
    F = R("3(x+1)^2(x^2+1)(x+5)").factor()
    assert [(str(u), e) for u, e in F] == [("x + 1", 2), ("x + 5", 1), ("x² + 1", 1)]
    assert F.unit == 3 and len(F) == 3 and F[0] == (R("x+1"), 2) and F.ring is R
    C = R(5).factor()
    assert len(C) == 0 and C.unit == 5 and C.expand() == 5
    with pytest.raises(ValueError, match="el polinomio cero no se puede factorizar"):
        R(0).factor()


def test_factorizacion_en_grados_distintos_e_iguales():
    R = GF(3)["x"]
    f = R("(x+1)(x+2)(x^2+1)(x^3+2x+1)")
    dd = f.distinct_degree_factorization()
    assert [(d, str(u)) for d, u in dd] == [(1, "x² + 2"), (2, "x² + 1"), (3, "x³ + 2x + 1")]
    assert R("2x^2 + 1").distinct_degree_factorization() == [(1, R("x^2 + 2"))]
    lineales = dd[0][1].equal_degree_factorization(1)
    assert lineales == [R("x + 1"), R("x + 2")]
    g = R("(x^2+1)(x^2+x+2)(x^2+2x+2)")
    assert g.equal_degree_factorization(2) == [R("x^2+1"), R("x^2+x+2"), R("x^2+2x+2")]
    K = GF(9)
    h = K["x"]("(x + a)(x + 1)(x + 2a + 1)")
    # Orden canónico del núcleo: en F_q, coordenadas (c0, c1) en orden lexicográfico.
    assert h.equal_degree_factorization(1) == [K["x"](s) for s in ("x+a", "x+1", "x+2a+1")]
    assert h.equal_degree_factorization(1) == [u for u, _ in h.factor()]
    with pytest.raises(ValueError, match="requiere un polinomio de grado >= 1"):
        R(2).distinct_degree_factorization()
    with pytest.raises(ValueError, match="no es libre de cuadrados"):
        R("(x+1)^2").distinct_degree_factorization()
    with pytest.raises(ValueError, match="no es producto de irreducibles distintos de grado 2"):
        f.equal_degree_factorization(2)
    with pytest.raises(ValueError, match="no es producto de irreducibles distintos de grado 1"):
        R("x^2+1").equal_degree_factorization(1)
    with pytest.raises(ValueError, match="r debe ser >= 1"):
        f.equal_degree_factorization(0)


@pytest.mark.parametrize("q", [5, 7, 8, 9, 16])
def test_raices_contra_fuerza_bruta(q):
    K = GF(q)
    R = K["x"]
    for d in (1, 2, 3, 6):
        for _ in range(4):
            f = R.random(d) * R("x - " + str(K.random().index % K.p)) ** 2
            esperado = []
            for e in K:
                if f(e) == 0:
                    esperado.append((e, f.multiplicity(R.x - e)))
            assert f.roots() == esperado
            assert f.roots(multiplicities=False) == [r for r, _ in esperado]
    assert R(1).roots() == [] and R(K.gen).roots() == []
    with pytest.raises(ValueError, match="el polinomio cero tiene todas las raíces"):
        R(0).roots()


def test_raices_sin_raices_y_una_raiz():
    R = GF(7)["x"]
    assert R("x^2 + 1").roots() == [] and R("3x + 1").roots() == [(GF(7)(2), 1)]
    assert R("(x-3)^4").roots() == [(GF(7)(3), 4)]


# ---------------------------------------------------------------------------
# Igualdad, hash y pickle
# ---------------------------------------------------------------------------
def test_igualdad_y_hash():
    R = GF(7)["x"]
    f = R("x^2 + 1")
    assert f == "x^2 + 1" and f == [1, 0, 1] and f == {2: 1, 0: 1} and f != "x"
    assert R(3) == 3 and 3 == R(3) and R(3) == GF(7)(3) and R(3) == 10 and hash(R(3)) == hash(3)
    assert hash(R(0)) == hash(0) and len({R(3), 3, GF(7)(3)}) == 1
    assert f != None and f != 1.5 and f != GF(7)["y"]("y^2+1") and f != "y"  # noqa: E711
    assert f.__eq__(None) is NotImplemented
    assert hash(f) == hash(R([1, 0, 1])) and hash(f) == hash(f)
    assert len({f, R("x^2+1"), R("x^2+2")}) == 2


def test_copia_y_pickle():
    for f in (GF(7)["x"]("x^3 + 2"), GF(9)["y"]("a y + 1"), GF(9)["x"](0)):
        assert copy.copy(f) is f and copy.deepcopy(f) is f
        g = pickle.loads(pickle.dumps(f))
        assert g == f and g.ring is f.ring


# ---------------------------------------------------------------------------
# Factorization como objeto
# ---------------------------------------------------------------------------
def test_objeto_factorization():
    F = GF(7)
    R = F["x"]
    fa = Factorization([(R("x+1"), 2), ("x + 2", 1)], 3, R)
    assert fa.expand() == R("3(x+1)^2(x+2)") and fa.value() == fa.expand()
    assert fa == R("3(x+1)^2(x+2)").factor() and fa != Factorization([], 3, R)
    assert fa.__eq__(3) is NotImplemented and fa != 3
    assert hash(fa) == hash(R("3(x+1)^2(x+2)").factor())
    assert pickle.loads(pickle.dumps(fa)) == fa
    assert Factorization([(R("x"), 1)], F(1)).ring is R
    vacia = Factorization([], F(3))
    assert vacia.ring is R and vacia.expand() == 3 and list(vacia) == []
    assert Factorization([], 3, F["t"]).ring is F["t"]
    with pytest.raises(TypeError, match="pásalo con ring="):
        Factorization([], 3)
    with pytest.raises(TypeError, match="pásalo con ring="):
        Factorization([("x", 1)], 3)
    with pytest.raises(TypeError, match="ring debe ser un PolynomialRing"):
        Factorization([], 3, F)
    with pytest.raises(TypeError, match="la multiplicidad debe ser un entero"):
        Factorization([(R("x"), 1.5)], 1)


# ---------------------------------------------------------------------------
# Propiedades (hypothesis, pocos ejemplos)
# ---------------------------------------------------------------------------
listas = st.lists(st.integers(min_value=-50, max_value=50), max_size=8)


@settings(max_examples=40, deadline=None)
@given(listas, listas)
def test_propiedad_division_euclidea(a, b):
    R = GF(9)["x"]
    f, g = R([R.field.from_index(x % 9) for x in a]), R([R.field.from_index(x % 9) for x in b])
    if g:
        q, r = divmod(f, g)
        assert f == q * g + r and r.degree < g.degree
    d = f.gcd(g)
    if d:
        assert f % d == 0 and g % d == 0


@settings(max_examples=25, deadline=None)
@given(listas)
def test_propiedad_factorizacion(a):
    f = GF(5)["x"](a)
    if f:
        assert f.factor().expand() == f
