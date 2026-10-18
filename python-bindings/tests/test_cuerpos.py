"""Pruebas de los cuerpos finitos y sus elementos (alcp.field).

Comparan la aritmética con implementaciones de referencia en Python puro
(tests/referencia.py) en cuerpos pequeños y recorren cada regla de conversión,
cada error (tipo y mensaje), la igualdad, el hash, pickle y las funciones de
teoría de cuerpos (traza, norma, orden, raíces cuadradas, polinomio mínimo).
"""
from __future__ import annotations

import copy
import decimal
import itertools
import operator
import pickle
from fractions import Fraction

import pytest

import alcp
from alcp import field as modulo_field
from referencia import Extension, pol_eval

GF = alcp.GF


class Indice:
    """Entero «de mentira» con __index__ (como numpy.int64)."""

    def __init__(self, v):
        self.v = v

    def __index__(self):
        return self.v


class IndiceRoto:
    """__index__ que lanza: la excepción se propaga tal cual."""

    def __index__(self):
        raise RuntimeError("__index__ roto")


@pytest.fixture(autouse=True)
def semilla():
    alcp.set_seed(20261018)
    yield
    alcp.set_print_style("unicode")


# ---------------------------------------------------------------------------
# Construcción y caché
# ---------------------------------------------------------------------------
def test_gf_primo_cacheado_y_con_enteros_generosos():
    F = GF(7)
    assert GF(7) is F
    assert GF(Indice(7)) is F
    assert GF(True + 6) is F
    assert isinstance(F, alcp.PrimeField) and isinstance(F, alcp.FiniteField)
    assert (F.p, F.characteristic, F.q, F.order, F.degree) == (7, 7, 7, 7, 1)
    assert F.var is None and F.is_prime_field
    assert F.prime_subfield is F


@pytest.mark.parametrize("q", [0, 1, 6, 12, -7, 100])
def test_gf_no_potencia_de_primo(q):
    with pytest.raises(ValueError, match="no es potencia de un número primo"):
        GF(q)


@pytest.mark.parametrize("q", [7.0, "7", None, Fraction(7)])
def test_gf_orden_no_entero(q):
    with pytest.raises(TypeError, match="el orden q debe ser un entero"):
        GF(q)


def test_gf_primo_no_admite_modulo():
    with pytest.raises(ValueError, match=r"GF\(7\) es un cuerpo primo: no admite módulo"):
        GF(7, modulus=[1, 1])


def test_gf_extension_cacheada_por_modulo_y_variable():
    K = GF(9)
    assert isinstance(K, alcp.ExtensionField)
    assert (K.p, K.q, K.n if hasattr(K, "n") else K.degree, K.var) == (3, 9, 2, "a")
    F3 = GF(3)
    for m in ("x^2+1", "t^2 + 1", "2x^2 + 2", [1, 0, 1], (2, 0, 2), [4, 3, 1],
              F3["x"]("x^2+1"), F3["y"]("y^2+1")):
        assert GF(9, modulus=m) is K, m
    Kb = GF(9, var="b")
    assert Kb is not K and Kb != K and Kb.var == "b"
    assert GF(9, modulus="x^2+1", var="b") is Kb
    assert K.prime_subfield is F3 and not K.is_prime_field


@pytest.mark.parametrize("q, esperado", [
    (4, [1, 1, 1]), (8, [1, 1, 0, 1]), (9, [1, 0, 1]), (16, [1, 1, 0, 0, 1]),
    (25, [1, 1, 1]), (27, [1, 2, 0, 1]), (49, [1, 0, 1]), (169, [2, 0, 1]),
    (2197, [2, 0, 0, 1]), (1024, [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]),
])
def test_modulo_determinista(q, esperado):
    K = GF(q)
    assert [int(c) for c in K.modulus] == esperado
    assert K.modulus.is_irreducible()


@pytest.mark.parametrize("p, n", [(2, 2), (2, 3), (3, 2), (3, 3), (5, 2), (5, 3),
                                  (7, 3), (13, 2), (11, 3)])
def test_modulo_determinista_es_el_primero_del_orden_documentado(p, n):
    # Para grado <= 3, irreducible <=> sin raíces en F_p (fuerza bruta).
    def irreducible(g):
        return all(pol_eval(g, x, p) for x in range(p))

    def orden():
        for cota in range(2, p + 1):
            for cifras in itertools.product(range(cota), repeat=n):
                if cifras[-1] and max(cifras) == cota - 1:
                    yield list(reversed(cifras)) + [1]

    primero = next(g for g in orden() if irreducible(g))
    assert [int(c) for c in GF(p ** n).modulus] == primero


def test_modulo_aleatorio_reproducible():
    alcp.set_seed(5)
    K1 = GF(3 ** 5, modulus="random")
    alcp.set_seed(5)
    K2 = GF(3 ** 5, modulus="random")
    assert K1 is K2
    assert K1.modulus.degree == 5 and K1.modulus.is_irreducible() and K1.modulus.is_monic()


@pytest.mark.parametrize("modulo, error, mensaje", [
    ("x^2+y", ValueError, "una sola variable"),
    ("2", ValueError, "tiene grado 0, pero q = 3\\^2 exige grado 2"),
    ("x^3+2x+1", ValueError, "tiene grado 3"),
    ("x^2+2", ValueError, r"el módulo x\^2 \+ 2 es reducible sobre GF\(3\)"),
    ([1, 1], ValueError, "tiene grado 1"),
    (5, TypeError, "el módulo debe ser None, 'random', una lista"),
    (1.5, TypeError, "no 'float'"),
])
def test_modulo_invalido(modulo, error, mensaje):
    with pytest.raises(error, match=mensaje):
        GF(9, modulus=modulo)


def test_modulo_polinomio_de_otro_cuerpo():
    with pytest.raises(ValueError, match=r"el módulo debe ser un polinomio sobre GF\(3\)"):
        GF(9, modulus=GF(5)["x"]("x^2+2"))


@pytest.mark.parametrize("var, error", [(1, TypeError), ("a1", ValueError), ("", ValueError),
                                        ("a b", ValueError)])
def test_nombre_de_variable_invalido(var, error):
    with pytest.raises(error, match="el nombre de la variable"):
        GF(9, var=var)


def test_no_se_instancian_directamente():
    with pytest.raises(TypeError, match=r"alcp.GF\(q, ...\)"):
        alcp.PrimeField(7)
    with pytest.raises(TypeError, match=r"K\(x\)"):
        alcp.FieldElement()


# ---------------------------------------------------------------------------
# Datos del cuerpo
# ---------------------------------------------------------------------------
def test_datos_cuerpo_primo():
    F = GF(7)
    assert F.zero == 0 and F.one == 1 and F.gen == 1
    assert str(F.modulus) == "x + 6"
    assert len(F) == 7 and [int(e) for e in F] == list(range(7))
    assert [e.index for e in F.elements()] == list(range(7))
    assert F.from_index(5) == 5 and F.to_index(-1) == 6 and F.to_index("3/2") == 5


def test_datos_extension():
    K = GF(9)
    a = K.gen
    assert a.index == 3 and K.from_index(3) == a and K.to_index("a + 2") == 5
    assert K.zero.index == 0 and K.one.index == 1
    assert str(K.modulus) == "x² + 1" and K.modulus.field is GF(3)
    assert len(K) == 9 and [e.index for e in K] == list(range(9))


@pytest.mark.parametrize("i", [-1, 9, 10**30])
def test_from_index_fuera_de_rango(i):
    with pytest.raises(ValueError, match=r"el índice debe estar en \[0, 9\)"):
        GF(9).from_index(i)


def test_from_index_tipo():
    with pytest.raises(TypeError, match="el índice debe ser un entero"):
        GF(9).from_index(1.0)


def test_len_de_un_cuerpo_enorme():
    K = GF(2 ** 127)
    assert K.order == 2 ** 127
    with pytest.raises(OverflowError, match=r"demasiados para len\(\): usa K.order"):
        len(K)


def test_enumeracion_demasiado_grande():
    with pytest.raises(ValueError, match="demasiados para enumerarlos"):
        GF(1000003).elements()
    with pytest.raises(ValueError, match="demasiados"):
        list(GF(1009 ** 2))


def test_random_en_el_cuerpo_y_reproducible():
    K = GF(3 ** 4)
    alcp.set_seed(1)
    xs = [K.random() for _ in range(20)]
    alcp.set_seed(1)
    assert [K.random() for _ in range(20)] == xs
    assert all(x in K and 0 <= x.index < 81 for x in xs)
    assert len({x.index for x in xs}) > 1


def test_igualdad_hash_copia_y_pickle_de_cuerpos():
    for K in (GF(7), GF(9), GF(9, var="b"), GF(2 ** 8)):
        assert copy.copy(K) is K and copy.deepcopy(K) is K
        assert pickle.loads(pickle.dumps(K)) is K
        assert hash(K) == hash(pickle.loads(pickle.dumps(K)))
    assert GF(7) != GF(11) and GF(7) != 7 and GF(7).__eq__(7) is NotImplemented
    assert GF(9) != GF(9, var="b") and GF(9) != GF(3)
    assert len({GF(7), GF(7), GF(9), GF(9, var="b")}) == 3


# ---------------------------------------------------------------------------
# Aritmética contra la referencia
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("q", [2, 4, 7, 8, 9, 25])
def test_tablas_contra_referencia(q):
    K = GF(q)
    ref = Extension(K.p, [int(c) for c in K.modulus])
    elems = list(K)
    coords = [ref.desde_indice(e.index) for e in elems]
    for x, cx in zip(elems, coords):
        assert tuple(int(c) for c in x.coeffs) == cx
        assert (-x).index == ref.indice(tuple(-c % K.p for c in cx))
        for y, cy in zip(elems, coords):
            assert (x + y).index == ref.indice(ref.suma(cx, cy))
            assert (x - y).index == ref.indice(ref.suma(cx, tuple(-c % K.p for c in cy)))
            assert (x * y).index == ref.indice(ref.mul(cx, cy))
            if y:
                assert (x / y) * y == x
        if x:
            assert x * x.inverse() == 1 and x ** -1 == x.inverse()
            assert x ** (q - 1) == 1


def test_potencias():
    K = GF(9)
    a = K.gen
    ref = Extension(3, [1, 0, 1])
    for k in range(12):
        assert (a ** k).index == ref.indice(ref.potencia((0, 1), k))
    assert a ** -1 == a ** 3 and a ** (10**40) == a ** (10**40 % 8)
    assert a ** Indice(2) == -1 and a ** True == a and +a is a
    assert K.zero ** 0 == 1 and K.zero ** 5 == 0


@pytest.mark.parametrize("k, error, mensaje", [
    (1.0, TypeError, "el exponente debe ser un entero"),
    ("2", TypeError, "el exponente debe ser un entero"),
    (GF(7)(2), TypeError, "no un elemento de un cuerpo"),
])
def test_potencia_exponente_invalido(k, error, mensaje):
    with pytest.raises(error, match=mensaje):
        GF(9).gen ** k


def test_potencia_tres_argumentos():
    with pytest.raises(TypeError, match="pow\\(\\) con tres argumentos"):
        pow(GF(7)(3), 2, 5)


def test_division_por_cero():
    F, K = GF(7), GF(9)
    for cero in (F(0), 0, "0", Fraction(0)):
        with pytest.raises(ZeroDivisionError, match=r"división por cero en GF\(7\)"):
            F(3) / cero
    with pytest.raises(ZeroDivisionError, match="división por cero en GF\\(9"):
        1 / K(0)
    with pytest.raises(ZeroDivisionError, match="no tiene inverso"):
        F(0).inverse()
    with pytest.raises(ZeroDivisionError, match="no tiene inverso"):
        K(0) ** -2


# ---------------------------------------------------------------------------
# Conversión (cada regla de la tabla)
# ---------------------------------------------------------------------------
def test_conversion_enteros():
    F, K = GF(7), GF(9)
    assert F() == 0 and F(10) == 3 and F(-1) == 6 and F(True) == 1 and F(10**50) == 10**50 % 7
    assert F(Indice(9)) == 2 and K(Indice(5)) == 2 and K(-1).index == 2
    assert K(5).index == 2                     # n·1 con n = 5 ≡ 2 (mod 3)
    with pytest.raises(RuntimeError, match="__index__ roto"):
        F(IndiceRoto())


def test_conversion_numpy():
    np = pytest.importorskip("numpy")
    F, K = GF(7), GF(9)
    assert F(np.int64(10)) == 3 and K(np.uint8(4)) == 1
    assert np.int64(2) + F(3) == 5 and F(3) * np.int32(3) == 2
    assert np.int64(2) - K.gen == 2 - K.gen
    assert isinstance(np.int64(2) * K.gen, alcp.ExtensionFieldElement)
    with pytest.raises(TypeError, match="usa enteros o Fraction"):
        F(np.float64(1.0))
    # Un array no escalar tiene __index__ pero es una secuencia de coeficientes.
    assert K(np.array([1, 2])) == 1 + 2 * K.gen and F(np.array(5)) == 5
    with pytest.raises(TypeError, match="solo valen en extensiones"):
        F(np.array([1, 2]))
    with pytest.raises(TypeError):
        F(1) + np.array([1, 2])
    with pytest.raises(TypeError, match="no se puede convertir un valor de tipo 'bool'"):
        F(np.True_)


class IndiceNoEntero:
    """Objeto con __index__ que rechaza la conversión (como numpy.bool_)."""

    def __index__(self):
        raise TypeError("no soy un entero")


class SecuenciaConIndice(IndiceNoEntero):
    """Secuencia con __index__ que falla (como un array de numpy no escalar)."""

    def __iter__(self):
        return iter([1, 2])


def test_conversion_index_que_falla():
    F, K = GF(7), GF(9)
    with pytest.raises(TypeError, match="un valor de tipo 'IndiceNoEntero' en un elemento"):
        F(IndiceNoEntero())
    assert K(SecuenciaConIndice()) == 1 + 2 * K.gen
    assert F(1).__add__(SecuenciaConIndice()) is NotImplemented
    assert GF(7)["x"](SecuenciaConIndice()) == "2x + 1"


def test_conversion_fracciones():
    F, K = GF(7), GF(9)
    assert F(Fraction(1, 2)) == 4 and F(Fraction(-3, 4)) * 4 == -3
    assert K(Fraction(1, 2)) == 2
    with pytest.raises(ZeroDivisionError, match="el denominador de 1/7 es múltiplo de p = 7"):
        F(Fraction(1, 7))


def test_conversion_cadenas():
    F, K = GF(7), GF(9)
    a = K.gen
    assert F("3") == 3 and F("1/2") == 4 and F("2^-1") == 4 and F("-(3+4)*2") == 0
    assert K("a") == a and K("a^2 + 1") == 0 and K("2a+1") == 2 * a + 1 and K("a²") == -1
    assert K("a⁻¹") == a ** -1 and K("(a+1)/(a+2)") == (a + 1) / (a + 2)
    with pytest.raises(ValueError, match="variable desconocida 'x'"):
        F("x")
    with pytest.raises(ZeroDivisionError):
        K("1/(a^2+1)")


def test_conversion_listas_en_extensiones():
    K = GF(9)
    a = K.gen
    assert K([1, 2]) == 1 + 2 * a and K((0, 1)) == a and K([]) == 0
    assert K([1, 0, 1]) == 0                     # se reduce módulo g = x² + 1
    assert K(["2", Fraction(1, 2), 3]) == 2 + 2 * a
    assert K([GF(3)(1), K(2)]) == 1 + 2 * a      # coeficientes del subcuerpo primo
    with pytest.raises(ValueError, match=r"a no pertenece al subcuerpo primo GF\(3\)"):
        K([a])
    with pytest.raises(TypeError, match="las listas de coeficientes solo valen en extensiones"):
        GF(7)([1, 2])


def test_conversion_entre_cuerpos():
    F3, K, Kb = GF(3), GF(9), GF(9, var="b")
    a = K.gen
    assert K(F3(2)) == 2 and K(F3(2)).field is K
    assert F3(K(2)) == 2 and F3(K(2)).field is F3
    assert Kb(K(2)).field is Kb and Kb(K(2)) == 2
    assert K(a) is a
    with pytest.raises(ValueError, match=r"a no pertenece al subcuerpo primo GF\(3\)"):
        F3(a)
    with pytest.raises(ValueError, match="son cuerpos distintos"):
        Kb(a)
    with pytest.raises(ValueError, match=r"GF\(5\) y GF\(7\) son cuerpos distintos"):
        GF(5)(GF(7)(3))
    with pytest.raises(ValueError, match=r"GF\(5\) y GF\(7\) son cuerpos distintos"):
        GF(5)(1) + GF(7)(1)


def test_conversion_polinomios_constantes():
    F, K = GF(7), GF(9)
    R = K["x"]
    assert F(F["x"](3)) == 3 and K(R("a + 1")) == K.gen + 1 and K(R(0)) == 0
    assert K(F3x := GF(3)["x"]("2")) == 2 and F3x.degree == 0
    with pytest.raises(ValueError, match=r"solo un polinomio constante.*no x \+ 1"):
        F(F["x"]("x+1"))


@pytest.mark.parametrize("valor, tipo", [(1.5, "float"), (1j, "complex"),
                                         (decimal.Decimal("1"), "Decimal")])
def test_conversion_numeros_no_enteros(valor, tipo):
    with pytest.raises(TypeError, match=f"un valor de tipo '{tipo}'.*usa enteros o Fraction"):
        GF(7)(valor)


@pytest.mark.parametrize("valor", [None, object(), {1: 2}, {1, 2}])
def test_conversion_tipos_no_admitidos(valor):
    with pytest.raises(TypeError, match="no se puede convertir un valor de tipo"):
        GF(9)(valor)


def test_pertenencia():
    F, K = GF(7), GF(9)
    assert 3 in F and "1/2" in F and Fraction(1, 3) in F and F(2) in F
    assert 1.5 not in F and "x" not in F and Fraction(1, 7) not in F and None not in F
    assert [1, 2] not in F and [1, 2] in K and K.gen not in F and GF(3)(1) in K
    assert GF(5)(1) not in F


# ---------------------------------------------------------------------------
# Operaciones mixtas y por la izquierda
# ---------------------------------------------------------------------------
def test_operaciones_con_reflejo():
    F, K = GF(7), GF(9)
    a = K.gen
    e = F(3)
    assert 2 + e == 5 and 2 - e == 6 and 2 * e == 6 and 2 / e == 3
    assert e + 2 == 5 and e - 2 == 1 and e * 2 == 6 and e / 2 == 5
    assert "3" * e == 2 and Fraction(1, 2) * e == 5 and Indice(2) + e == 5
    assert "a" * a == -1 and a * "a" == -1 and [0, 1] * a == -1 and 1 / a == -a
    assert True + e == 4


def test_operaciones_promocion_de_fp_a_extension():
    F3, K = GF(3), GF(9)
    a = K.gen
    dos = F3(2)
    for r in (dos + a, a + dos, dos - a, a - dos, dos * a, a * dos, dos / a, a / dos):
        assert r.field is K
    assert dos + a == a + 2 and dos - a == 2 - a and dos / a == 2 / a and a / dos == a / 2
    # Las variantes reflejadas también funcionan si se llaman directamente.
    assert dos.__radd__(a) == a + 2 and dos.__rsub__(a) == a - 2
    assert dos.__rtruediv__(a) == a / 2 and dos.__rmul__(a) == 2 * a


def test_operandos_no_admitidos():
    F = GF(7)
    for malo in (None, object(), {1: 2}):
        with pytest.raises(TypeError):
            F(1) + malo
        with pytest.raises(TypeError):
            malo * F(1)
    with pytest.raises(TypeError, match="usa enteros o Fraction"):
        F(1) + 1.5
    with pytest.raises(TypeError, match="usa enteros o Fraction"):
        1.5 * F(1)
    with pytest.raises(TypeError, match="solo valen en extensiones"):
        F(1) + [1]
    assert F(1).__add__(None) is NotImplemented


def test_extensiones_distintas_no_se_mezclan():
    a, b = GF(9).gen, GF(9, var="b").gen
    with pytest.raises(ValueError, match="son cuerpos distintos"):
        a + b
    with pytest.raises(ValueError, match="son cuerpos distintos"):
        b * GF(9)(2)               # ni siquiera con elementos del subcuerpo primo


# ---------------------------------------------------------------------------
# Igualdad, hash, conversiones a int, pickle
# ---------------------------------------------------------------------------
def test_igualdad_con_convertibles():
    F, K = GF(7), GF(9)
    a = K.gen
    assert F(3) == 3 and F(3) == 10 and F(3) == "3" and F(3) == Fraction(3) and 3 == F(3)
    assert F(3) != 4 and F(3) != 1.5 and F(3) != None and F(3) != "x"  # noqa: E711
    assert F(3) != [3] and F(3) != GF(5)(3) and F(3) != Fraction(1, 7) and F(3) != object()
    assert K(2) == GF(3)(2) and GF(3)(2) == K(2) and K(2) == 2
    assert a == "a" and a == [0, 1] and a != GF(3)(1) and a != GF(9, var="b").gen
    assert F(3) == F["x"](3) and F["x"](3) == F(3) and F(3) != F["x"]("x")


def test_hash_coherente():
    F, K = GF(7), GF(9)
    assert hash(F(3)) == hash(3) and len({F(3), 3, F(10)}) == 1
    assert hash(K(2)) == hash(2) == hash(GF(3)(2)) and len({K(2), GF(3)(2), 2}) == 1
    d = {F(3): "tres"}
    assert d[3] == "tres" and d[F(10)] == "tres"
    assert hash(F["x"](3)) == hash(F(3))


def test_bool_int_index():
    F, K = GF(7), GF(9)
    assert not F(0) and F(3) and not K(0) and K.gen
    assert int(F(3)) == 3 and operator.index(F(5)) == 5 and list(range(F(3))) == [0, 1, 2]
    assert ["a", "b", "c"][F(1)] == "b" and float(F(3)) == 3.0
    assert int(K(2)) == 2
    with pytest.raises(TypeError, match=r"a \+ 1 no está en el subcuerpo primo GF\(3\)"):
        int(K.gen + 1)
    with pytest.raises(TypeError):
        operator.index(K(2))


def test_datos_de_elementos():
    K = GF(9)
    e = K("2a + 1")
    assert e.field is K and e.parent is K and e.index == 7
    assert e.coeffs == [1, 2] and all(c.field is GF(3) for c in e.coeffs)
    F = GF(7)
    assert F(3).coeffs == [3] and F(3).coeffs[0].field is F


def test_copia_y_pickle_de_elementos():
    for e in (GF(7)(3), GF(9)("2a+1"), GF(2 ** 8).gen ** 7, GF(9, var="b").gen):
        assert copy.copy(e) is e and copy.deepcopy(e) is e
        f = pickle.loads(pickle.dumps(e))
        assert f == e and f.field is e.field and type(f) is type(e)


# ---------------------------------------------------------------------------
# Teoría de cuerpos
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("q", [7, 8, 9, 16, 27])
def test_frobenius_traza_norma(q):
    K = GF(q)
    p, n = K.p, K.degree
    for x in K:
        conj = [x ** (p ** i) for i in range(n)]
        assert x.frobenius() == x ** p and x.frobenius(n) == x and x.frobenius(-1) == conj[-1]
        t, nm = sum(conj, K.zero), K.one
        for c in conj:
            nm = nm * c
        assert x.trace() == t and x.trace().field is GF(p)
        assert x.norm() == nm and x.norm().field is GF(p)


@pytest.mark.parametrize("q", [2, 7, 16, 25])
def test_orden_multiplicativo(q):
    K = GF(q)
    for x in list(K)[1:]:
        k, y = 1, x
        while y != 1:
            y, k = y * x, k + 1
        assert x.multiplicative_order() == k
    with pytest.raises(ValueError, match="el 0 no tiene orden multiplicativo"):
        K.zero.multiplicative_order()


def test_orden_multiplicativo_con_q_menos_1_no_factorizable(monkeypatch):
    # 443 - 1 = 2·13·17: con un límite de división 10, el cofactor 221 no se factoriza.
    from alcp import _util
    monkeypatch.setattr(_util, "LIMITE_DIVISION", 10)
    with pytest.raises(ValueError, match="no se puede factorizar el entero"):
        GF(443)(3).multiplicative_order()


@pytest.mark.parametrize("q", [2, 3, 7, 8, 9, 13, 16, 17, 25, 41, 97, 81])
def test_raices_cuadradas_contra_fuerza_bruta(q):
    K = GF(q)
    cuadrados = {}
    for y in K:
        cuadrados.setdefault((y * y).index, set()).add(y.index)
    for x in K:
        raices = sorted(cuadrados.get(x.index, ()))
        assert x.is_square() == bool(raices)
        assert [r.index for r in x.sqrt(all=True)] == raices
        if raices:
            assert x.sqrt().index == raices[0] and x.sqrt() ** 2 == x
        else:
            with pytest.raises(ValueError, match="no es un cuadrado en"):
                x.sqrt()


@pytest.mark.parametrize("q", [7, 8, 9, 16, 27])
def test_polinomio_minimo(q):
    K = GF(q)
    for x in K:
        m = x.minimal_polynomial()
        assert m.field is GF(K.p) and m.is_monic() and m.is_irreducible()
        assert K.degree % m.degree == 0 and m(x) == 0
    assert GF(7)(3).minimal_polynomial("t") == GF(7)["t"]("t - 3")


def test_raices_de_la_unidad():
    K = GF(25)
    for n in (1, 2, 3, 4, 6, 8, 12, 24):
        w = K.root_of_unity(n)
        assert w.multiplicative_order() == n
        assert K.root_of_unity(Indice(n)) == w        # determinista
    g = K.primitive_element()
    assert g.multiplicative_order() == 24 and K.primitive_element() is g
    assert GF(2).primitive_element() == 1
    for n in (0, -1, 5, 7):
        with pytest.raises(ValueError, match="no tiene raíces primitivas de la unidad de orden"):
            K.root_of_unity(n)


# ---------------------------------------------------------------------------
# Detalles internos (caché del módulo determinista, reconstrucción)
# ---------------------------------------------------------------------------
def test_reconstruccion_de_elementos_valida_el_indice():
    with pytest.raises(ValueError, match="el índice debe estar en"):
        modulo_field._elemento(GF(9), 9)
    assert modulo_field._elemento(GF(9), 3) == GF(9).gen
