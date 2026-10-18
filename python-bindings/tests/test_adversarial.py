"""Revisión adversarial del paquete ``alcp``: entradas hostiles y casos límite.

Pruebas escritas por un revisor independiente contra la API declarada (no
contra la implementación): conversiones raras, mezclas de cuerpos y anillos,
característica 2, cuerpos grandes, polinomio cero y constantes, exponentes
negativos y enormes, división por cero, pickle/copia, coherencia de == y hash,
textos exactos, determinismo de ``GF(q)`` y comparación exhaustiva con fuerza
bruta (factorización, raíces, órdenes, traza y norma, algoritmos rápidos).

Para que la suite siga siendo rápida, casi todas las pruebas recorren una tabla
de casos dentro de una sola función.
"""
from __future__ import annotations

import copy
import itertools
import math
import pickle
import re
from decimal import Decimal
from fractions import Fraction

import pytest

import alcp
from alcp import GF, Factorization, Polynomial, PolynomialRing
from referencia import pol_divmod, pol_eval

pytestmark = pytest.mark.adversarial

M61 = 2**61 - 1
M127 = 2**127 - 1
M521 = 2**521 - 1


class Indice:
    """Entero «de mentira» con __index__ (como numpy.int64)."""

    def __init__(self, v):
        self.v = v

    def __index__(self):
        return self.v


class IndiceRoto:
    """__index__ que lanza: la excepción debe propagarse (API declarada)."""

    def __index__(self):
        raise RuntimeError("__index__ roto")


def lineas_de_error(exc):
    """(motivo, línea de la expresión, línea del ^) de un error del analizador."""
    motivo, expresion, marca = str(exc.value).split("\n")
    return motivo, expresion, marca


# ---------------------------------------------------------------------------
# Analizador
# ---------------------------------------------------------------------------
def test_cadenas_raras_que_se_aceptan():
    F, K = GF(7), GF(9)
    R, S = F["x"], K["x"]
    casos = [
        (R, "  x ^ 2\t+\n1 ", "x^2 + 1"),
        (R, "x² − 3·x ÷ 2 ⋅ 1 × 1 ∗ 1", "x^2 + 2*x"),
        (R, "2x(x+1)(x-1)", "2*x^3 + 5*x"),
        (R, "((((x))))^(2)", "x^2"),
        (R, "x**2 + x^+2 + x^(-0)", "2*x^2 + 1"),
        (R, "++x - -x + - - -x", "x"),
        (R, "2(3) + x(2) + (x)(x) + 2 x", "x^2 + 4*x + 6"),
        (R, "2⁻¹ + x⁰ + x¹⁰ + x⁺²", "x^10 + x^2 + 5"),
        (R, "0^0 + 2^99999999999999999999", "2"),
        (R, "1/2/2 + 3/4x", "6*x + 2"),
        (R, "-x^2", "6*x^2"),
        (R, "(-x)^2", "x^2"),
        (S, "ax + xa + a x + x a", "a*x"),
        (S, "x/a + x÷(a+1) + a⁻¹x", "(2*a + 2)*x"),
        (S, "(a+1)(x-1)", "(a + 1)*x + 2*a + 2"),
        (S, "3a^2x - a²x²", "x^2"),
        (K, "a^2 + 1", "0"),
        (K, "a/(a+1)", "2*a + 2"),
        (F, "3/4 - -1", "0"),
        (GF(7)["α"], "α² + 2α", "α^2 + 2*α"),
        (GF(7)["X"], "X^2 + 2X", "X^2 + 2*X"),
    ]
    for destino, texto, esperado in casos:
        assert format(destino(texto), "ascii") == esperado, texto


ERRORES_DE_SINTAXIS = [
    # (texto, regex del motivo, columna del ^)
    ("", "expresión vacía", 0),
    ("   ", "expresión vacía", 3),
    ("()", "se esperaba un número", 1),
    ("(x", r"falta '\)'", 2),
    ("x)", r"paréntesis '\)' sin abrir", 1),
    ("2 3", "falta un operador antes de '3'", 2),
    ("x2", "falta un operador antes de '2'", 1),
    ("(x+1)2", "falta un operador antes de '2'", 5),
    ("1.5", "carácter inesperado '.'.*decimales", 1),
    ("1,5", "carácter inesperado ','.*decimales", 1),
    ("y", "variable desconocida 'y' \\(variables válidas: x\\)", 0),
    ("xy", "variable desconocida 'y'", 1),
    ("x^", "el exponente debe ser un número entero", 2),
    ("x^^2", "el exponente debe ser un número entero", 2),
    ("x***2", "el exponente debe ser un número entero", 3),
    ("x^x", "el exponente debe ser un número entero", 2),
    ("x^(2", r"falta '\)'", 4),
    ("x^2^3", "símbolo inesperado", 3),
    ("x⁻", "superíndice no válido: '⁻'", 1),
    ("x²⁻", "superíndice no válido", 1),
    ("x++", "se esperaba un número", 3),
    ("x = 1", "carácter inesperado '='", 2),
    ("２x", "carácter inesperado '２'", 0),
    ("x\t#", "carácter inesperado '#'", 2),
    ("²", "se esperaba un número", 0),
    ("*x", "se esperaba un número", 0),
]


@pytest.mark.parametrize("texto, motivo, columna", ERRORES_DE_SINTAXIS,
                         ids=[repr(c[0]) for c in ERRORES_DE_SINTAXIS])
def test_errores_de_sintaxis_con_su_posicion(texto, motivo, columna):
    with pytest.raises(ValueError) as exc:
        GF(7)["x"](texto)
    linea_motivo, expresion, marca = lineas_de_error(exc)
    assert re.search(motivo, linea_motivo), linea_motivo
    assert expresion == "    " + texto.replace("\t", " ")
    assert marca == "    " + " " * columna + "^"


def test_variables_desconocidas_segun_el_destino():
    with pytest.raises(ValueError, match="'a' \\(variables válidas: ninguna\\)"):
        GF(7)("a")
    with pytest.raises(ValueError, match="'x' \\(variables válidas: a\\)"):
        GF(9)("x")
    with pytest.raises(ValueError, match="'b' \\(variables válidas: a, x\\)"):
        GF(9)["x"]("b")
    # Una variable del cuerpo no vale en un anillo sobre F_p (aunque exista GF(9)).
    with pytest.raises(ValueError, match="variable desconocida 'a'"):
        GF(3)["x"]("x + a")
    assert "a" not in GF(3)["x"] and "a" in GF(9)["x"] and "x" not in GF(9)


def test_conflictos_de_nombre_entre_cuerpo_y_anillo():
    """Nombres que son prefijo unos de otros: el texto unicode y ASCII se relee igual."""
    combinaciones = [("ab", "a"), ("ab", "b"), ("b", "ba"), ("a", "aa"), ("aa", "a"),
                     ("xy", "x"), ("x", "xy"), ("x", "y"), ("A", "a"), ("alpha", "t")]
    for var_cuerpo, var_anillo in combinaciones:
        K = GF(9, var=var_cuerpo)
        R = K[var_anillo]
        a, t = K.gen, R.gen
        for f in (a * t**2 + 2 * a * t + a, (a + 1) * t**3 + a**2 * t + 2,
                  t**11 + a, R(a + 2), R.random(4), R.random(5)):
            for estilo in ("unicode", "ascii"):
                assert R(format(f, estilo)) == f, (var_cuerpo, var_anillo, format(f, estilo))
    # Coincidencia más larga: con a y aa, «aaa» es aa·a (anillo: aa, cuerpo: a).
    R = GF(9)["aa"]
    assert R("aaa") == R.gen * GF(9).gen and R("aaaa") == R.gen**2
    with pytest.raises(ValueError, match="ya es el generador"):
        GF(9, var="t")["t"]


def test_errores_de_calculo_conservan_su_tipo_y_la_posicion():
    K, R = GF(9), GF(7)["x"]
    casos = [
        (lambda: K("1/(a^2+1)"), ZeroDivisionError, "división por cero en GF\\(9", 1),
        (lambda: K("a^2 + (a+1)^-1 * 0^-3"), ZeroDivisionError, "no tiene inverso", 18),
        (lambda: R("x/0"), ZeroDivisionError, "división por el polinomio cero", 1),
        (lambda: R("x + 1/(7)"), ZeroDivisionError, "división por el polinomio cero", 5),
        (lambda: R("(x+1)^-1"), ValueError, "no constante no tiene inverso", 5),
        (lambda: R("2 + x/x"), ValueError, "solo se puede dividir por constantes", 5),
        (lambda: R("x⁻¹"), ValueError, "no constante no tiene inverso", 1),
        (lambda: R("x^10000001"), ValueError, "demasiado grande", 1),
        (lambda: GF(7)("3/7"), ZeroDivisionError, "división por cero en GF\\(7\\)", 1),
    ]
    for funcion, tipo, motivo, columna in casos:
        with pytest.raises(tipo, match=motivo) as exc:
            funcion()
        assert type(exc.value) is tipo
        assert str(exc.value).split("\n")[2] == "    " + " " * columna + "^"


def test_anidamiento_y_expresiones_largas():
    R = GF(7)["x"]
    limite = alcp._parse.PROFUNDIDAD_MAXIMA
    assert R("(" * limite + "x" + ")" * limite) == R.gen
    with pytest.raises(ValueError, match="demasiado anidada"):
        R("(" * (limite + 1) + "x" + ")" * (limite + 1))
    assert R("+".join(["x"] * 3000)) == 3000 * R.gen
    assert R("-" * 20001 + "x") == -R.gen


# ---------------------------------------------------------------------------
# Conversiones y mezclas
# ---------------------------------------------------------------------------
def test_conversiones_raras_de_elementos():
    np = pytest.importorskip("numpy")
    F, K = GF(7), GF(9)
    validos = [
        (F, True, 1), (F, False, 0), (F, Indice(-1), 6), (F, np.int8(-3), 4),
        (F, np.uint64(2**64 - 1), (2**64 - 1) % 7), (F, 10**100, 10**100 % 7),
        (F, -10**100, -10**100 % 7), (F, Fraction(-1, 3), 2), (F, Fraction(14, 3), 0),
        (F, Fraction(10**50 + 1, 3), (10**50 + 1) * 5 % 7), (F, " - 1 ", 6), (F, np.array(3), 3),
        (K, [1, 2, 3], 7), (K, ("1", GF(3)(2)), 7), (K, [Fraction(1, 2)], 2), (K, [], 0),
        (K, np.array([1, 2]), 7), (K, GF(3)(2), 2), (K, GF(27)(2), 2), (K, GF(9, var="b")(1), 1),
        (K, GF(9)["y"](K.gen), 3), (K, GF(3)["x"](2), 2),
    ]
    for destino, valor, indice in validos:
        assert destino(valor).index == indice and valor in destino, (destino, valor)
    invalidos = [
        (F, 1.0, TypeError, "usa enteros o Fraction"), (F, 2j, TypeError, "usa enteros o Fraction"),
        (F, Decimal("1"), TypeError, "usa enteros o Fraction"),
        (F, np.float64(1), TypeError, "usa enteros o Fraction"),
        (F, None, TypeError, "'NoneType'"), (F, b"1", TypeError, "'bytes'"),
        (F, {1}, TypeError, "'set'"), (F, {1: 2}, TypeError, "'dict'"),
        (F, range(2), TypeError, "'range'"), (F, [1], TypeError, "solo valen en extensiones"),
        (F, np.array([1, 2]), TypeError, "solo valen en extensiones"),
        (F, np.bool_(True), TypeError, "no se puede convertir"),
        (K, [[1]], TypeError, "solo valen en extensiones"), (K, [1.5], TypeError, "usa enteros"),
        (K, [K.gen], ValueError, "no pertenece al subcuerpo primo GF\\(3\\)"),
        (F, Fraction(1, 14), ZeroDivisionError, "múltiplo de p = 7"),
        (F, "a", ValueError, "variable desconocida"), (F, GF(5)(1), ValueError, "cuerpos distintos"),
        (GF(3), K.gen, ValueError, "a no pertenece al subcuerpo primo GF\\(3\\)"),
        (K, GF(27).gen, ValueError, "son cuerpos distintos"),
        (K, GF(9, var="b").gen, ValueError, "son cuerpos distintos"),
        (K, GF(9)["x"].gen, ValueError, "solo un polinomio constante"),
    ]
    for destino, valor, tipo, mensaje in invalidos:
        with pytest.raises(tipo, match=mensaje):
            destino(valor)
        assert valor not in destino, (destino, valor)
    # Las excepciones de __index__ se propagan tal cual.
    with pytest.raises(RuntimeError, match="__index__ roto"):
        F(IndiceRoto())
    # El mismo cuerpo devuelve el mismo objeto.
    e = K.gen
    assert K(e) is e and F(F(3)) == 3


def test_mezcla_de_cuerpos_en_las_operaciones():
    F3, K9, K9b, K27, F5 = GF(3), GF(9), GF(9, var="b"), GF(27), GF(5)
    a, b = K9.gen, K9b.gen
    # F_p con una extensión: el resultado vive en la extensión (las 4 operaciones y su reflejo).
    for x, y in ((F3(2), a), (a, F3(2)), (F3(1), K27.gen)):
        L = y.field if x.field.is_prime_field else x.field
        for op in ("__add__", "__sub__", "__mul__", "__truediv__"):
            r = getattr(x, op)(y)
            assert r.field is L and r == getattr(L(x), op)(L(y))
    with pytest.raises(ZeroDivisionError, match="división por cero en GF\\(9"):
        F3(1) / K9(0)
    # Cuerpos distintos: ValueError en las operaciones, False (sin excepción) en ==.
    for x, y in ((a, b), (a, K27.gen), (F3(1), F5(1)), (a, F5(1)), (K27.gen, a)):
        for op in ("__add__", "__radd__", "__sub__", "__rsub__", "__mul__", "__truediv__"):
            with pytest.raises(ValueError, match="cuerpos distintos"):
                getattr(x, op)(y)
        assert (x == y) is False and (x != y) is True
    # El subcuerpo primo de dos extensiones distintas sí es común (y el hash coincide).
    assert K9(2) == K9b(2) == K27(2) == F3(2) == 2
    assert len({K9(2), K9b(2), K27(2), F3(2), 2}) == 1
    # Conversión explícita: del subcuerpo primo, sí; de fuera, no.
    assert K9(K27(2)) == 2 and K9b(K9(2)) == 2


def test_mezcla_de_anillos_en_las_operaciones():
    F3, K9, K9b = GF(3), GF(9), GF(9, var="b")
    f, g = F3["x"]("x^2 + 1"), K9["x"]("x + a")
    for r in (f + g, g + f, f * g, g * f, f - g, g - f, f // g, f % g, f.gcd(g), g.gcd(f)):
        assert r.ring is K9["x"]
    assert divmod(f, g) == divmod(K9["x"](f), g) and f.gcd(g) == g
    assert f + K9.gen == K9["x"]("x^2 + a + 1") == K9.gen + f
    assert F3(2) + g == g + 2 and (F3(2) * g).ring is K9["x"]
    distintos = [
        (F3["y"].gen, g, "anillos distintos"), (K9b["x"].gen, g, "anillos distintos"),
        (GF(5)["x"].gen, f, "anillos distintos"), (f, GF(5)(1), "cuerpos distintos"),
        (K9["y"].gen, K9["x"].gen, "anillos distintos"),
        # F_3[a] no se puede subir a F_9[a]: la variable a es el generador de F_9.
        (F3["a"].gen, K9.gen, "ya es el generador"),
    ]
    for u, v, mensaje in distintos:
        for op in ("__add__", "__sub__", "__mul__", "__floordiv__", "__mod__"):
            with pytest.raises(ValueError, match=mensaje):
                getattr(u, op)(v)
        assert (u == v) is False and (v == u) is False and u != v
    # Reflejo con int, str, elementos, listas y dicts.
    x = F3["x"].gen
    for otro in (2, "x + 1", F3(2), [1, 1], (1, 1), {0: 1, 1: 1}, Fraction(1, 2)):
        R = F3["x"]
        assert otro + x == R(otro) + x and otro * x == R(otro) * x
        assert otro - x == R(otro) - x and x - otro == x - R(otro)
    assert "x^3 + 1" // (x + 1) == R("x^2 + 2x + 1")
    assert divmod("x^3 + 2", x + 1) == (R("x^2 + 2x + 1"), R(1))
    # Con una cadena a la izquierda, % es el formateo de cadenas de Python (no se puede
    # interceptar desde el operando derecho): hay que convertir antes.
    assert R("x^3 + 2") % (x + 1) == 1


# ---------------------------------------------------------------------------
# Característica 2
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("q", [2, 4, 8, 256])
def test_caracteristica_2_aritmetica_de_elementos(q):
    K = GF(q)
    elementos = list(K)
    trazas = [e.trace() for e in elementos]
    for e in elementos:
        assert -e == e and e + e == 0 and e.frobenius() == e * e
        assert e.is_square() and e.sqrt() ** 2 == e and e.sqrt(all=True) == [e.sqrt()]
        assert e.frobenius(K.degree) == e and e.frobenius(-1).frobenius() == e
    assert len({e.frobenius() for e in elementos}) == q              # Frobenius biyectivo
    assert trazas.count(GF(2)(0)) == trazas.count(GF(2)(1)) == q // 2 or q == 2
    for x, y in itertools.islice(itertools.product(elementos, repeat=2), 300):
        assert (x + y).trace() == x.trace() + y.trace() and (x * y).norm() == x.norm() * y.norm()
    for valor in (Fraction(1, 2), "1/2", "a/2" if q > 2 else "1/(1+1)"):
        with pytest.raises(ZeroDivisionError):
            K(valor)
    R = K["x"]
    with pytest.raises(ZeroDivisionError, match="polinomio cero"):
        R.gen / 2
    assert R("x^2 + 1") == R("(x + 1)^2") and (R.gen + 1) ** 2 == R.gen**2 + 1


@pytest.mark.parametrize("q", [2, 4, 8, 16])
def test_caracteristica_2_polinomios(q):
    K = GF(q)
    R = K["x"]
    x = R.gen
    # f' = 0: f = g(x^2) es un cuadrado; radical, libre de cuadrados y factorización.
    f = x**4 + x**2 + 1                                           # (x^2 + x + 1)^2
    assert f.derivative() == 0 and not f.is_squarefree()
    assert f.radical() == (x**2 + x + 1).radical()
    assert f.factor().expand() == f and all(e % 2 == 0 for _, e in f.factor())
    assert (x**8 + x).factor().expand() == x**8 + x
    assert R(K.gen if q > 2 else 1).radical() == 1
    raices = (x**q - x).roots(multiplicities=False)
    assert raices == list(K)                                      # x^q - x = ∏ (x - c)
    assert ((x + 1) ** 4 * x).roots() == [(K(0), 1), (K(1), 4)]


# ---------------------------------------------------------------------------
# Cuerpos grandes
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("p", [M127, M521], ids=["M127", "M521"])
def test_cuerpos_primos_grandes(p):
    F = GF(p)
    menos_uno = F(-1)
    assert menos_uno.index == p - 1 and str(menos_uno) == str(p - 1)
    assert F(p) == 0 and F(2 * p + 3) == 3 and F(-p - 1) == menos_uno
    assert hash(menos_uno) == hash(p - 1) and menos_uno == p - 1 and menos_uno == -1
    for v in (2, 3, 10**30, p - 2):
        e = F(v)
        assert e ** (p - 1) == 1 and e * e.inverse() == 1 and e / e == 1 and e ** p == e
        assert (e * e).sqrt() in (e, -e) and len((e * e).sqrt(all=True)) == 2
    assert F.root_of_unity(2) == menos_uno and F.root_of_unity(1) == 1
    assert eval(repr(menos_uno), {"GF": GF}) == menos_uno
    assert pickle.loads(pickle.dumps(menos_uno)) == menos_uno
    R = F["x"]
    f = R("(x - 2)(x + 3)^2 (x^2 + 1)")
    assert f.factor().expand() == f and f.roots()[0] == (F(2), 1)
    assert alcp.dft(F, [1, 2]) == [F(3), F(-1)]
    if p == M127:
        # (Factorizar p - 1 por división de prueba cuesta: solo una vez.)
        # p - 1 = 2·3^3·7^2·19·43·73·127·337·5419·92737·649657·77158673929
        primos = [2, 3, 7, 19, 43, 73, 127, 337, 5419, 92737, 649657, 77158673929]
        assert 2 * 3**3 * 7**2 * math.prod(primos[3:]) == p - 1
        o = F(3).multiplicative_order()
        assert (p - 1) % o == 0 and F(3) ** o == 1
        assert all(F(3) ** (o // r) != 1 for r in primos if o % r == 0)
    else:
        with pytest.raises(ValueError, match="no se puede factorizar"):
            F(3).multiplicative_order()
    with pytest.raises(OverflowError, match="usa K.order"):
        len(F)
    with pytest.raises(ValueError, match="demasiados para enumerarlos"):
        F.elements()


def test_extension_con_p_grande():
    K = GF(M61**2)
    a = K.gen
    assert format(K.modulus, "ascii") == "x^2 + 1" and a**2 == -1
    assert a.frobenius() == -a and a.trace() == 0 and a.norm() == 1    # p ≡ 3 (mod 4)
    assert a.minimal_polynomial("y") == GF(M61)["y"]("y^2 + 1")
    # (a + 1)^2 = 2a, así que (a + 1)^8 = 16 y 16 = 2^4 tiene orden 61 (2^61 = 1 en F_p).
    assert a.multiplicative_order() == 4 and (a + 1).multiplicative_order() == 8 * 61
    R = K["x"]
    f = R("(x - a)(x - a - 1)(x + 5)")
    assert f.roots(multiplicities=False) == [K(-5), a, a + 1]
    F = R("x^3 - a").factor()
    assert F.expand() == R("x^3 - a") and all(u.is_irreducible() for u, _ in F)
    assert pickle.loads(pickle.dumps(a + 7)) == a + 7 and eval(repr(a + 7), {"GF": GF}) == a + 7


# ---------------------------------------------------------------------------
# Polinomio cero y constantes
# ---------------------------------------------------------------------------
def test_polinomio_cero_en_todas_las_operaciones():
    F = GF(7)
    R = F["x"]
    z, x = R.zero, R.gen
    assert z + z == z * x == -z == z**5 == z.monic() == z.derivative() == z.radical() == 0
    assert z**0 == 1 and z ** (10**30) == 0 and z.compose(x + 1) == 0 and (x + 1).compose(z) == 1
    assert (z.degree, len(z), list(z), z.coeffs, z.coefficients(sparse=True), z.to_dict()) == \
        (-1, 0, [], [], [], {})
    assert z.lc == 0 and z.leading_coefficient == 0 and not z and z.is_zero() and z.is_constant()
    assert not (z.is_one() or z.is_monic() or z.is_irreducible() or z.is_squarefree())
    assert z.gcd(z) == 0 and z.gcd(2 * x + 2) == x + 1 and z.xgcd(z) == (z, z, z)
    assert divmod(z, x) == (0, 0) and z // x == z % x == 0 and z(3) == 0 and z("x+1") == 0
    assert z[0] == z[5] == 0 and z[::-1] == [] and str(z) == "0" and repr(z) == "GF(7)['x']('0')"
    assert hash(z) == hash(0) == hash(F(0)) and z == 0 == F(0) and z == "0" and z == [] and z == {}
    assert pickle.loads(pickle.dumps(z)) == z and copy.deepcopy(z) is z
    errores = [
        (lambda: z ** -1, ZeroDivisionError, "cero no tiene inverso"),
        (lambda: x / z, ZeroDivisionError, "división por el polinomio cero"),
        (lambda: x // z, ZeroDivisionError, "división por el polinomio cero"),
        (lambda: x % z, ZeroDivisionError, "división por el polinomio cero"),
        (lambda: divmod(x, z), ZeroDivisionError, "división por el polinomio cero"),
        (lambda: pow(x, 3, z), ZeroDivisionError, "módulo cero"),
        (lambda: x.inverse_mod(z), ZeroDivisionError, "módulo el polinomio cero"),
        (lambda: z.inverse_mod(x + 1), ValueError, "no es invertible"),
        (lambda: z.factor(), ValueError, "cero no se puede factorizar"),
        (lambda: z.roots(), ValueError, "todas las raíces"),
        (lambda: z.multiplicity(x), ValueError, "infinita"),
        (lambda: x.multiplicity(z), ValueError, "grado >= 1"),
        (lambda: z.distinct_degree_factorization(), ValueError, "grado >= 1"),
        (lambda: z.equal_degree_factorization(1), ValueError, "grado >= 1"),
        (lambda: z[-1], IndexError, "fuera de rango"),
    ]
    for funcion, tipo, mensaje in errores:
        with pytest.raises(tipo, match=mensaje):
            funcion()


def test_constantes_en_todas_las_operaciones():
    F = GF(7)
    R = F["x"]
    c, x = R(3), R.gen
    assert c.degree == 0 and c.is_constant() and c.monic() == 1 and c.derivative() == 0
    assert c.radical() == 1 and c.is_squarefree() and not c.is_irreducible() and not c.is_monic()
    F_c = c.factor()
    assert len(F_c) == 0 and F_c.unit == 3 and F_c.expand() == c and str(F_c) == "3"
    assert c.roots() == [] and c.roots(multiplicities=False) == []
    assert c**-1 == R(5) and c ** -(10**30) == R(F(3) ** -(10**30)) and c ** (10**30) == F(3) ** (10**30)
    assert pow(c, -1, x**2 + 1) == 5 and c.inverse_mod(x**2 + 1) == 5 and c / c == 1
    assert (x**2 + 1) / c == 5 * x**2 + 5 and 6 / c == 2 and c(4) == 3 and c("x+1") == 3
    assert c.compose(x**3) == 3 and (x**2).compose(c) == 2 and c.multiplicity(x) == 0
    assert c.gcd(x**2 + 1) == 1 and c.xgcd(x) == (R(1), R(5), R(0))
    assert hash(c) == hash(3) == hash(F(3)) and c == 3 == F(3) == "3" and c != 4
    with pytest.raises(ValueError, match="grado >= 1"):
        c.distinct_degree_factorization()
    with pytest.raises(ValueError, match="no constante"):
        (x + 1) ** -1


# ---------------------------------------------------------------------------
# Exponentes
# ---------------------------------------------------------------------------
def test_exponentes_negativos_y_enormes():
    np = pytest.importorskip("numpy")
    for q in (7, 9, 2**8, M61):
        K = GF(q)
        e = K.gen + 2 if q != M61 else K(3)
        k = 10**40 + 12345
        assert e**k == e ** (k % (q - 1)) and e**-k == (e**-1) ** k == e ** (-k % (q - 1))
        assert e**0 == 1 and K(0) ** 0 == 1 and K(0) ** k == 0 and e ** True == e
        assert e ** np.int64(3) == e**3 == e ** Indice(3)
        for malo, tipo in ((K(2), TypeError), (2.0, TypeError), ("2", TypeError)):
            with pytest.raises(tipo):
                e**malo
        with pytest.raises(TypeError, match="tres argumentos"):
            pow(e, 2, 5)
        with pytest.raises(ZeroDivisionError, match="no tiene inverso"):
            K(0) ** -k
    R = GF(7)["x"]
    x, m = R.gen, R("x^3 + 3")                            # irreducible: F_7[x]/(m) = F_343
    k = 10**50 + 7
    assert pow(x, k, m) == pow(x, k % 342, m) and pow(x, -k, m) == pow(x.inverse_mod(m), k, m)
    assert pow(x, 0, m) == 1 and pow(x, 0, R(1)) == 0 and pow(x + 1, 342, m) == 1
    with pytest.raises(ValueError, match="no es invertible"):
        pow(x**2, -1, x**3)
    with pytest.raises(ValueError, match="demasiado grande"):
        x ** (10**7 + 1)
    with pytest.raises(TypeError, match="no un elemento"):
        x ** GF(7)(2)


# ---------------------------------------------------------------------------
# División por cero
# ---------------------------------------------------------------------------
def test_division_por_cero_en_todas_partes():
    F, K = GF(7), GF(9)
    R = F["x"]
    x = R.gen
    casos = [
        lambda: F(3) / 0, lambda: F(3) / "0", lambda: F(3) / F(0), lambda: 0 / F(0),
        lambda: 1 / K(0), lambda: K.gen / Fraction(3, 2), lambda: F(0) ** -1, lambda: F(0).inverse(),
        lambda: F(Fraction(1, 7)), lambda: F("1/7"), lambda: K("1/(a^2+1)"),
        lambda: x / 0, lambda: x // 0, lambda: x % 0, lambda: divmod(x, 0), lambda: 3 / R(0),
        lambda: 3 % R(0), lambda: divmod(3, R(0)), lambda: "x" // R(0), lambda: pow(x, 2, 0),
        lambda: x.inverse_mod(0), lambda: alcp.divmod_toeplitz(x, 0), lambda: R("x/0"),
        lambda: R("0^-1"), lambda: R(0) ** -1, lambda: R("1/7"), lambda: x / F(0),
        lambda: alcp.toeplitz_lower_inverse(F, [0, 1]), lambda: alcp.toeplitz_upper_inverse(K, ["0"]),
        lambda: GF(3)(1) / K(0), lambda: GF(3)["x"].gen / K(0),
    ]
    for funcion in casos:
        with pytest.raises(ZeroDivisionError) as exc:
            funcion()
        assert type(exc.value) is ZeroDivisionError
        assert any(palabra in str(exc.value) for palabra in ("cero", "múltiplo", "singular",
                                                             "inverso")), str(exc.value)


# ---------------------------------------------------------------------------
# Pickle y copia
# ---------------------------------------------------------------------------
def test_pickle_y_copia_de_todo():
    alcp.set_seed(99)
    aleatorio = GF(3**5, modulus="random")
    K, Kb, F = GF(9), GF(9, var="beta"), GF(M127)
    R = K["t"]
    objetos = [
        GF(2), K, Kb, aleatorio, F, GF(M61**2), R, Kb["x"], F["x"],
        K.gen, Kb(2), aleatorio.gen ** 7, F(-1), GF(M61**2).gen + 3,
        R.zero, R("(a+1)t^3 + a"), F["x"]("x^2 - 5"),
        R("(a+1)(t+a)^2 t").factor(), R(K.gen + 1).factor(), Kb["x"].one.factor(),
    ]
    for o in objetos:
        for protocolo in range(pickle.HIGHEST_PROTOCOL + 1):
            r = pickle.loads(pickle.dumps(o, protocol=protocolo))
            assert type(r) is type(o) and r == o and hash(r) == hash(o) and repr(r) == repr(o)
            if isinstance(o, (alcp.FiniteField, PolynomialRing)):
                assert r is o                                      # cuerpos y anillos cacheados
        if not isinstance(o, Factorization):
            assert copy.copy(o) is o and copy.deepcopy(o) is o
        else:
            assert copy.copy(o) == o and copy.deepcopy(o) == o
    g = K.gen
    lista = [g, R.gen, {g: R.gen}]
    copia = copy.deepcopy(lista)
    assert copia == lista and copia[0] is g and copia[2][g] is R.gen


# ---------------------------------------------------------------------------
# == y hash
# ---------------------------------------------------------------------------
def _valores_para_igualdad():
    np = pytest.importorskip("numpy")
    F3, K9, K9b, K27, F5 = GF(3), GF(9), GF(9, var="b"), GF(27), GF(5)
    valores = [F3(0), F3(1), F3(2), K9(0), K9(1), K9(2), K9.gen, K9.gen + 1, K9b(2), K9b.gen,
               K27(2), K27.gen, F5(0), F5(2)]
    for R in (F3["x"], K9["x"], K9["y"], K9b["x"], F5["x"]):
        valores += [R(0), R(1), R(2), R([1, 1]), R([2, 0, 1])]
    valores += [K9["x"](K9.gen), K9["x"]("x + a"), F3["x"]("x+1").factor(),
                K9["x"]("x+1").factor(), F3["x"]("2x^2 + 2").factor()]
    # Solo enteros del rango canónico [0, p) de todos los cuerpos de la lista (p >= 3).
    valores += [0, 1, 2, True, False, Fraction(2), np.int64(2), np.uint8(1)]
    return valores


def test_igualdad_implica_el_mismo_hash():
    """a == b ⇒ hash(a) == hash(b), también con enteros del rango canónico [0, p),
    bool, Fraction enteras y enteros de numpy; == es simétrica y nunca lanza."""
    valores = _valores_para_igualdad()
    for a, b in itertools.product(valores, repeat=2):
        if not any(type(v).__module__.startswith("alcp") for v in (a, b)):
            continue                                 # int con numpy, etc.: no es asunto de alcp
        igual = a == b
        assert igual is True or igual is False, (a, b)
        assert bool(igual) == bool(b == a), (a, b)
        assert bool(a != b) == (not igual), (a, b)
        if igual:
            assert hash(a) == hash(b), (a, b)
    # Conjuntos y diccionarios se comportan como con números.
    assert len({GF(7)(3), 3, GF(7)["x"](3), Fraction(3)}) == 1
    assert {GF(9)(2): "x"}[GF(3)(2)] == "x"


def test_igualdad_fuera_del_rango_canonico_documentada():
    """GF(7)(3) == 10 es True (API declarada); el hash es el del representante en
    [0, p), así que hash(10) no coincide: es la excepción documentada (como en SageMath)."""
    e = GF(7)(3)
    for otro in (10, -4, "3", "10", Fraction(13, 2), GF(7)["x"]("10")):
        assert e == otro and otro == e
    assert hash(e) == hash(3) != hash(10)


def test_comparacion_con_tipos_no_convertibles():
    np = pytest.importorskip("numpy")
    K = GF(9)
    R = K["x"]
    alcp_objs = [K.gen, GF(7)(3), R.gen, R(2), R("x+a").factor(), K, R]
    raros = [None, 1.5, 2j, Decimal("1"), object(), {1}, frozenset(), "zzz", "x +", b"x",
             [[1]], {"a": 1}, {-1: 1}, np.float64(1.0), len, float("nan"), [K.gen, 1.5]]
    for o in alcp_objs:
        for r in raros:
            assert (o == r) is False and (r == o) is False and o != r, (o, r)
        with pytest.raises(TypeError):
            o < o  # noqa: B015  (no hay orden)


# ---------------------------------------------------------------------------
# Textos exactos
# ---------------------------------------------------------------------------
def test_textos_exactos_adversariales():
    K4, Ka = GF(4), GF(9, var="alpha")
    R4, Ra = K4["x"], Ka["xx"]
    F3 = GF(3)["x"]
    casos = [
        # (objeto, unicode, ascii, latex, repr)
        (GF(2), "GF(2)", "GF(2)", r"\mathbb{F}_{2}", "GF(2)"),
        (K4, "GF(4) = GF(2)[a]/(a² + a + 1)", "GF(4) = GF(2)[a]/(a^2 + a + 1)",
         r"\mathbb{F}_{4} = \mathbb{F}_{2}[a]/(a^{2} + a + 1)", "GF(4, modulus='x^2 + x + 1')"),
        (Ka, "GF(9) = GF(3)[alpha]/(alpha² + 1)", "GF(9) = GF(3)[alpha]/(alpha^2 + 1)",
         r"\mathbb{F}_{9} = \mathbb{F}_{3}[\alpha]/(\alpha^{2} + 1)",
         "GF(9, modulus='x^2 + 1', var='alpha')"),
        (Ra, "GF(9)[xx]", "GF(9)[xx]", r"\mathbb{F}_{9}[\mathit{xx}]",
         "GF(9, modulus='x^2 + 1', var='alpha')[xx]"),
        (R4("(a+1)x^2 + a x + 1"), "(a + 1)x² + ax + 1", "(a + 1)*x^2 + a*x + 1",
         r"\left(a + 1\right) x^{2} + a x + 1", "GF(4, modulus='x^2 + x + 1')['x']('(a + 1)*x^2 + a*x + 1')"),
        (R4(K4.gen + 1), "a + 1", "a + 1", "a + 1", "GF(4, modulus='x^2 + x + 1')['x']('a + 1')"),
        (Ra("(alpha+1)xx^2 + 2alpha xx + 1"), "(alpha + 1)·xx² + 2·alpha·xx + 1",
         "(alpha + 1)*xx^2 + 2*alpha*xx + 1",
         r"\left(\alpha + 1\right) \mathit{xx}^{2} + 2 \alpha \mathit{xx} + 1",
         "GF(9, modulus='x^2 + 1', var='alpha')['xx']('(alpha + 1)*xx^2 + 2*alpha*xx + 1')"),
        (Ka.gen * 2 + 1, "2·alpha + 1", "2*alpha + 1", r"2 \alpha + 1",
         "GF(9, modulus='x^2 + 1', var='alpha')('2*alpha + 1')"),
        (GF(27)["x"]("a^2 x^2 + (a^2+1) x + 2a^2"), "a²x² + (a² + 1)x + 2a²",
         "a^2*x^2 + (a^2 + 1)*x + 2*a^2", r"a^{2} x^{2} + \left(a^{2} + 1\right) x + 2 a^{2}",
         "GF(27, modulus='x^3 + 2*x + 1')['x']('a^2*x^2 + (a^2 + 1)*x + 2*a^2')"),
        (F3("2(x+1)^3").factor(), "2 · (x + 1)³", "2 * (x + 1)^3", r"2 \cdot \left(x + 1\right)^{3}",
         "Factorization([(GF(3)['x']('x + 1'), 3)], unit=GF(3)(2))"),
        (F3("x + 1").factor(), "x + 1", "x + 1", "x + 1",
         "Factorization([(GF(3)['x']('x + 1'), 1)], unit=GF(3)(1))"),
        (F3("x^2").factor(), "x²", "x^2", "x^{2}", "Factorization([(GF(3)['x']('x'), 2)], unit=GF(3)(1))"),
        (F3("2x").factor(), "2 · x", "2 * x", r"2 \cdot x",
         "Factorization([(GF(3)['x']('x'), 1)], unit=GF(3)(2))"),
        (F3(1).factor(), "1", "1", "1", "Factorization([], unit=GF(3)(1), ring=GF(3)['x'])"),
        (Ra(Ka.gen + 1).factor(), "alpha + 1", "alpha + 1", r"\alpha + 1",
         "Factorization([], unit=GF(9, modulus='x^2 + 1', var='alpha')('alpha + 1'), "
         "ring=GF(9, modulus='x^2 + 1', var='alpha')['xx'])"),
        (GF(9)["x"]("(a+1)(x+1)^2 (x^2+x+a)").factor(), "(a + 1) · (x + 1)² · (x² + x + a)",
         "(a + 1) * (x + 1)^2 * (x^2 + x + a)",
         r"\left(a + 1\right) \cdot \left(x + 1\right)^{2} \cdot \left(x^{2} + x + a\right)", None),
    ]
    for obj, uni, asc, lat, rep in casos:
        assert (str(obj), format(obj, "ascii"), format(obj, "latex")) == (uni, asc, lat), uni
        assert obj._repr_latex_() == f"${lat}$"
        if rep is not None:
            assert repr(obj) == rep
            if not isinstance(obj, (alcp.FiniteField, PolynomialRing)):
                assert eval(rep, {"GF": GF, "Factorization": Factorization}) == obj


def test_especificadores_de_formato():
    e = GF(257)(255)
    numericos = {"d": "255", "+d": "+255", "05d": "00255", "x": "ff", "X": "FF", "b": "11111111",
                 "o": "377", "n": "255", "#x": "0xff"}
    for spec, esperado in numericos.items():
        assert format(e, spec) == esperado == f"{e:{spec}}"
    assert (f"{e:>6}", f"{e:^7}", f"{e:<5}|", f"{e:u}", f"{e:a}", f"{e:l}") == \
        ("   255", "  255  ", "255  |", "255", "255", "255")
    a = GF(9).gen + 1
    assert (f"{a:>7}", f"{a:l}", f"{a:latex}") == ("  a + 1", "a + 1", "a + 1")
    with pytest.raises(ValueError):
        format(a, "d")                              # un elemento de una extensión no es un entero
    f = GF(9)["x"]("(a+1)x^2")
    with alcp.print_style("ascii"):
        assert f"{f}" == "(a + 1)*x^2" and f"{f:>14}" == "   (a + 1)*x^2" and f"{f:u}" == "(a + 1)x²"
    assert f"{f}" == "(a + 1)x²" and f"{f:*^13}" == "**(a + 1)x²**"


def test_estilos_de_impresion_robustos():
    assert alcp.get_print_style() == "unicode"
    for malo in ("UNICODE", "", None, 1, "tex"):
        with pytest.raises(ValueError, match="estilo de impresión desconocido"):
            alcp.set_print_style(malo)
        with pytest.raises(ValueError):
            with alcp.print_style(malo):
                pass                                                   # pragma: no cover
        assert alcp.get_print_style() == "unicode"
    with alcp.print_style("l") as estilo:
        assert estilo == "latex"
        with pytest.raises(KeyError):
            with alcp.print_style("a"):
                assert alcp.get_print_style() == "ascii"
                raise KeyError("x")
        assert alcp.get_print_style() == "latex" and str(GF(7)["x"]("x^2")) == "x^{2}"
    assert alcp.get_print_style() == "unicode"


def test_texto_unicode_y_ascii_se_vuelven_a_leer():
    """parse(str(f)) == f para polinomios, elementos y factorizaciones en varios cuerpos."""
    for q, var, rvar in ((2, "a", "x"), (4, "a", "t"), (9, "alpha", "x"), (27, "a", "xx"),
                         (25, "b", "T"), (M61, "a", "x")):
        K = GF(q, var=var)
        R = K[rvar]
        for grado in range(-1, 7):
            f = R.random(grado)
            for estilo in ("unicode", "ascii"):
                assert R(format(f, estilo)) == f
                if f and grado <= 5 and q < 10**6:
                    assert R(format(f.factor(), estilo)) == f
        for e in [K.random() for _ in range(10)] + [K.zero, K.one, K.gen]:
            assert K(str(e)) == e and K(format(e, "ascii")) == e


# ---------------------------------------------------------------------------
# GF(q) sin módulo: determinismo frente a una enumeración independiente
# ---------------------------------------------------------------------------
def _irreducible_a_mano(g, p):
    """g (lista creciente, mónico) es irreducible sobre F_p: ningún mónico de grado
    1..n/2 lo divide (los de grado 1, por evaluación)."""
    n = len(g) - 1
    if any(pol_eval(g, c, p) == 0 for c in range(p)):
        return False
    for d in range(2, n // 2 + 1):
        for cola in itertools.product(range(p), repeat=d):
            if not pol_divmod(g, list(cola) + [1], p)[1]:
                return False
    return True


def _primer_irreducible_documentado(p, n):
    """Primer mónico irreducible de grado n en el orden documentado de GF(q):
    cotas B = 2, 3, ...; para cada una, cifras < B con alguna igual a B - 1, por
    índice c0 + c1·B + ... creciente."""
    B = 2
    while True:
        for indice in range(B**n):
            cifras = [indice // B**k % B for k in range(n)]
            if max(cifras) == B - 1 and _irreducible_a_mano(cifras + [1], p):
                return cifras + [1]
        B += 1


@pytest.mark.parametrize("q", [4, 8, 16, 32, 64, 128, 256, 9, 27, 81, 243, 729, 25, 125, 625,
                               49, 343, 121, 169, 289, 361, 529, 961, 1009**2])
def test_modulo_determinista_contra_enumeracion_independiente(q):
    p, n = alcp.prime_power_factor(q)
    esperado = _primer_irreducible_documentado(p, n)
    K = GF(q)
    assert [c.index for c in K.modulus.coeffs] == esperado
    assert GF(q) is K is GF(q, modulus=esperado) and pickle.loads(pickle.dumps(K)) is K
    assert K.modulus.is_irreducible() and K.modulus.is_monic()


def test_modulo_determinista_en_cuerpos_grandes():
    # Los primos de Mersenne son ≡ 3 (mod 4): x^2 + 1 (el primer candidato) es irreducible.
    for p in (M61, M127, M521):
        assert format(GF(p**2).modulus, "ascii") == "x^2 + 1"
    assert format(GF(2**127).modulus, "ascii") == "x^127 + x + 1"


# ---------------------------------------------------------------------------
# Factorización y raíces contra fuerza bruta (exhaustivo en cuerpos pequeños)
# ---------------------------------------------------------------------------
def _irreducibles_por_criba(R, grado_maximo):
    """Mónicos irreducibles de grado <= grado_maximo: los que no divide ninguno menor."""
    K = R.field
    irreducibles = []
    for d in range(1, grado_maximo + 1):
        for cola in itertools.product(range(K.order), repeat=d):
            f = R([K.from_index(i) for i in cola] + [1])
            if all(f % u for u in irreducibles if 2 * u.degree <= d):
                irreducibles.append(f)
    return irreducibles


def _clave_canonica(u):
    """Orden canónico: grado y coeficientes desde el de mayor grado, cada uno por coordenadas."""
    return (u.degree, [tuple(c.index for c in coef.coeffs) for coef in reversed(u.coeffs)])


@pytest.mark.parametrize("q, grado_maximo", [(2, 6), (3, 4), (4, 3), (5, 3), (9, 2)])
def test_factorizacion_y_raices_exhaustivas_contra_fuerza_bruta(q, grado_maximo):
    K = GF(q)
    R = K["x"]
    x = R.gen
    irreducibles = _irreducibles_por_criba(R, grado_maximo)
    # Número de mónicos irreducibles de cada grado (fórmula de Gauss, comprobación de la criba).
    for d in range(1, grado_maximo + 1):
        gauss = sum(_mobius(d // k) * q**k for k in range(1, d + 1) if d % k == 0) // d
        assert sum(1 for u in irreducibles if u.degree == d) == gauss
    assert all(u.is_irreducible() for u in irreducibles)
    elementos = list(K)
    for d in range(0, grado_maximo + 1):
        for cola in itertools.product(elementos, repeat=d):
            for principal in (elementos[1:] if d < grado_maximo or q <= 3 else [K.one]):
                f = R(list(cola) + [principal])
                g, factores = f.monic(), []
                for u in irreducibles:
                    e = 0
                    while g.degree >= u.degree and not g % u:
                        g, e = g // u, e + 1
                    if e:
                        factores.append((u, e))
                assert g == 1
                F = f.factor()
                canonica = sorted(factores, key=lambda t: _clave_canonica(t[0]))
                assert F.unit == principal and list(F) == canonica
                assert f.is_irreducible() == (len(factores) == 1 and factores[0][1] == 1)
                assert f.is_squarefree() == all(e == 1 for _, e in factores)
                bruta = []
                for c in elementos:
                    m, h = 0, f
                    while h.degree >= 1 and not h % (x - c):
                        h, m = h // (x - c), m + 1
                    if m:
                        bruta.append((c, m))
                assert f.roots() == bruta and f.roots(multiplicities=False) == [c for c, _ in bruta]


def _mobius(n):
    resultado, d = 1, 2
    while d * d <= n:
        if n % d == 0:
            n //= d
            if n % d == 0:
                return 0
            resultado = -resultado
        d += 1
    return -resultado if n > 1 else resultado


# ---------------------------------------------------------------------------
# Orden multiplicativo, traza y norma (exhaustivo)
# ---------------------------------------------------------------------------
def _phi(n):
    return sum(1 for k in range(1, n + 1) if math.gcd(k, n) == 1)


@pytest.mark.parametrize("q", [2, 3, 4, 5, 8, 9, 16, 25, 27, 32, 49])
def test_orden_multiplicativo_traza_y_norma_exhaustivos(q):
    K = GF(q)
    p, n = K.p, K.degree
    Fp = GF(p)
    ordenes = {}
    for e in list(K)[1:]:
        potencia, orden = e, 1
        while potencia != 1:
            potencia, orden = potencia * e, orden + 1
        assert e.multiplicative_order() == orden
        ordenes[orden] = ordenes.get(orden, 0) + 1
        conjugados = [e ** (p**k) for k in range(n)]
        assert e.trace() == sum(conjugados, K(0)) and e.norm() == math.prod(conjugados, start=K(1))
        assert e.trace().field is Fp and e.norm().field is Fp and e.norm() != 0
    # Para cada d | q - 1 hay φ(d) elementos de orden d.
    assert ordenes == {d: _phi(d) for d in range(1, q) if (q - 1) % d == 0}
    assert K.primitive_element().multiplicative_order() == q - 1
    for d in ordenes:
        assert K.root_of_unity(d).multiplicative_order() == d
    with pytest.raises(ValueError, match="no tiene orden multiplicativo"):
        K(0).multiplicative_order()
    # Traza sobreyectiva y equilibrada; norma sobreyectiva sobre F_p^*.
    trazas = [e.trace().index for e in K]
    normas = [e.norm().index for e in list(K)[1:]]
    assert all(trazas.count(t) == q // p for t in range(p))
    assert all(normas.count(v) == (q - 1) // (p - 1) for v in range(1, p))
    # En F_p (como subcuerpo): traza n·c y norma c^n.
    for c in range(p):
        assert K(c).trace() == n * c and K(c).norm() == c**n


# ---------------------------------------------------------------------------
# Algoritmos rápidos contra versiones ingenuas
# ---------------------------------------------------------------------------
def _matriz_por_vector(T, b):
    return [sum((T[i][j] * b[j] for j in range(len(b))), b[0] * 0) for i in range(len(T))]


@pytest.mark.parametrize("q", [17, 97, 9, 25, 4, M61, M61**2])
def test_algoritmos_rapidos_contra_versiones_ingenuas(q):
    K = GF(q)
    alcp.set_seed(5)
    # DFT de todas las longitudes 2^k que admite el cuerpo (hasta 16).
    n = 1
    while (q - 1) % n == 0 and n <= 16:
        a = [K.random() for _ in range(n)]
        w = K.root_of_unity(n)
        ingenua = [sum((a[j] * w ** (i * j) for j in range(n)), K(0)) for i in range(n)]
        assert alcp.dft(K, a) == alcp.dft(K, a, w) == ingenua
        assert alcp.idft(K, ingenua) == a and alcp.idft(K, alcp.dft(K, a, w**-1), w**-1) == a
        n *= 2
    for n in range(1, 6):
        a = [K.random() for _ in range(2 * n - 1)]
        b = [K.random() for _ in range(n)]
        a[0] = a[0] or K.one
        inferior = [[a[i - j] if i >= j else K(0) for j in range(n)] for i in range(n)]
        superior = [[a[j - i] if j >= i else K(0) for j in range(n)] for i in range(n)]
        completa = [[a[j - i] if j >= i else a[n - 1 + i - j] for j in range(n)] for i in range(n)]
        assert alcp.toeplitz_lower_vec(K, a[:n], b) == _matriz_por_vector(inferior, b)
        assert alcp.toeplitz_upper_vec(K, a[:n], b) == _matriz_por_vector(superior, b)
        assert alcp.toeplitz_vec(K, a, b) == _matriz_por_vector(completa, b)
        for inversa, T in ((alcp.toeplitz_lower_inverse, inferior), (alcp.toeplitz_upper_inverse, superior)):
            v = inversa(K, a[:n])
            # La inversa también es de Toeplitz: T · (primera columna/fila de T^-1) = e_1.
            if inversa is alcp.toeplitz_lower_inverse:
                assert _matriz_por_vector(T, v) == [K(1)] + [K(0)] * (n - 1)
            else:
                assert _matriz_por_vector([list(f) for f in zip(*T)], v) == [K(1)] + [K(0)] * (n - 1)
    R = K["x"]
    for gf, gg in ((7, 3), (3, 7), (5, 0), (0, 2)):
        f, g = R.random(gf), R.random(gg)
        assert alcp.divmod_toeplitz(f, g) == divmod(f, g)


def test_algoritmos_rapidos_entradas_hostiles():
    F, K = GF(17), GF(9)
    errores = [
        (lambda: alcp.dft(F, []), ValueError, "potencia de 2, no 0"),
        (lambda: alcp.dft(F, [1, 2, 3]), ValueError, "potencia de 2, no 3"),
        (lambda: alcp.dft(GF(7), [1, 2, 3, 4]), ValueError, "4 no divide a q - 1 = 6"),
        (lambda: alcp.dft(F, [1, 2, 3, 4], root=16), ValueError, "orden exactamente"),
        (lambda: alcp.dft(F, [1, 2, 3, 4], root=3), ValueError, "orden exactamente"),
        (lambda: alcp.idft(F, [1, 2], root=1), ValueError, "orden exactamente"),
        (lambda: alcp.dft(F, [1, 2], root=GF(5)(4)), ValueError, "cuerpos distintos"),
        (lambda: alcp.dft(7, [1]), TypeError, "cuerpo finito"),
        (lambda: alcp.dft(F, "12"), TypeError, "lista de elementos"),
        (lambda: alcp.dft(F, {0: 1}), TypeError, "lista de elementos"),
        (lambda: alcp.dft(F, 5), TypeError, "lista de elementos"),
        (lambda: alcp.dft(F, [1.5, 2]), TypeError, "usa enteros"),
        (lambda: alcp.toeplitz_lower_vec(F, [1], [1, 1]), ValueError, "al menos n = 2"),
        (lambda: alcp.toeplitz_lower_vec(F, [1, 2], [], 0), ValueError, "n debe ser >= 1"),
        (lambda: alcp.toeplitz_upper_vec(F, [1], [1], -1), ValueError, "n debe ser >= 0"),
        (lambda: alcp.toeplitz_vec(F, [1, 2], [1, 1]), ValueError, "exactamente 3"),
        (lambda: alcp.toeplitz_lower_inverse(K, []), ValueError, "n debe ser >= 1"),
        (lambda: alcp.toeplitz_upper_inverse(K, [1], 3), ValueError, "al menos n = 3"),
        (lambda: alcp.toeplitz_lower_vec(F, [1], [1], 10**30), OverflowError, "demasiado grande"),
        (lambda: alcp.divmod_toeplitz("x^2", GF(7)["x"].gen), TypeError, "debe ser un polinomio"),
        (lambda: alcp.divmod_toeplitz(GF(7)["x"].gen, GF(5)["x"].gen), ValueError, "anillos distintos"),
        (lambda: alcp.root_of_unity(F, 3), ValueError, "3 no divide a q - 1 = 16"),
    ]
    for funcion, tipo, mensaje in errores:
        with pytest.raises(tipo, match=mensaje):
            funcion()
    # Iterables generales, cadenas como elementos y raíces dadas como cadenas.
    assert alcp.dft(F, (i for i in range(4))) == alcp.dft(F, [0, 1, 2, 3])
    # En GF(9) = F_3[a]/(a^2 + 1), a tiene orden 4.
    assert alcp.dft(K, ["a", "1", 0, K.gen], root="a") == alcp.dft(K, [K.gen, 1, 0, K.gen], K.gen)
    assert alcp.dft(GF(4), ["a"]) == [GF(4).gen] and alcp.idft(GF(2), [1]) == [GF(2)(1)]


# ---------------------------------------------------------------------------
# Inmutabilidad y construcción
# ---------------------------------------------------------------------------
def test_inmutabilidad_y_construccion_directa():
    K = GF(9)
    R = K["x"]
    for objeto in (K.gen, GF(7)(3), R.gen, R("x+a").factor()):
        assert not hasattr(objeto, "__dict__")
        with pytest.raises(AttributeError):
            objeto.nuevo = 1
    for propiedad, objeto in (("index", K.gen), ("degree", R.gen), ("unit", R("2x").factor()),
                              ("q", K), ("var", R)):
        with pytest.raises(AttributeError):
            setattr(objeto, propiedad, 1)
    for clase, argumentos in ((alcp.FiniteField, (7,)), (alcp.PrimeField, (7,)),
                              (alcp.ExtensionField, (9,)), (alcp.FieldElement, ()),
                              (alcp.PrimeFieldElement, (3,)), (alcp.ExtensionFieldElement, ())):
        with pytest.raises(TypeError):
            clase(*argumentos)
    assert Polynomial(R) == 0 and Polynomial(R, "x + a") == R("x + a")
    with pytest.raises(TypeError, match="PolynomialRing"):
        Polynomial(K, [1])
    with pytest.raises(TypeError, match="cuerpo finito"):
        PolynomialRing("GF(7)")
