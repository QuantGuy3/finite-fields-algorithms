"""Pruebas del analizador de expresiones (alcp._parse) y de su uso en cuerpos y anillos.

La gramática se prueba primero con enteros y fracciones de Python (operaciones
exactas, independientes del núcleo) y después en cuerpos y anillos: variables,
multiplicación implícita, superíndices Unicode, símbolos alternativos y cada
mensaje de error con su ``^`` señalando la posición.
"""
from __future__ import annotations

import operator
from fractions import Fraction

import pytest

import alcp
from alcp import _parse

GF = alcp.GF


def evaluar(texto, **variables):
    """Evalúa con enteros/fracciones de Python (división exacta)."""
    return _parse.analizar(texto, Fraction, {k: Fraction(v) for k, v in variables.items()})


def mensaje(texto, pos, motivo):
    return f"{motivo}\n    {texto}\n    {' ' * pos}^"


# ---------------------------------------------------------------------------
# Gramática con números racionales
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("texto, valor", [
    ("1 + 2 * 3", 7), ("(1 + 2) * 3", 9), ("2^3", 8), ("2**3", 8),
    ("-2^2", -4), ("(-2)^2", 4), ("2^-1", Fraction(1, 2)), ("2^(-2)", Fraction(1, 4)),
    ("2^(+3)", 8), ("2^+3", 8), ("2²", 4), ("2¹⁰", 1024), ("2⁻¹", Fraction(1, 2)),
    ("2⁺³", 8), ("7 - 2 - 1", 4), ("8 / 2 / 2", 2), ("2 · 3 × 4 ⋅ 5 ∗ 1", 120),
    ("9 ÷ 3", 3), ("5 − 2", 3), ("--3", 3), ("-+-+3", 3), ("+3", 3), ("---3", -3),
    ("2(3)", 6), ("(2)(3)", 6), ("2 (3 + 1)", 8), ("\t1\n+\r2 ", 3),
    ("123456789012345678901234567890 + 1", 123456789012345678901234567891),
    ("x", 5), ("2x", 10), ("2 x", 10), ("x²", 25), ("3x^2y", 3 * 25 * 2), ("x y", 10),
    ("xy", 10), ("yx", 10), ("(x+1)(y-1)", 6), ("x(y)", 10), ("x/2x", Fraction(25, 2)),
    ("2x²y³", 2 * 25 * 8),
])
def test_gramatica(texto, valor):
    assert evaluar(texto, x=5, y=2) == valor


def test_cadena_larga_de_signos_y_anidamiento_permitido():
    assert evaluar("-" * 5001 + "1") == -1
    profundo = "(" * _parse.PROFUNDIDAD_MAXIMA + "1" + ")" * _parse.PROFUNDIDAD_MAXIMA
    assert evaluar(profundo) == 1


def test_reparto_de_letras_en_nombres():
    # Con los nombres «alpha» y «x», «alphax» es alpha·x (coincidencia más larga).
    assert evaluar("alphax", alpha=3, x=5) == 15
    assert evaluar("xalpha", alpha=3, x=5) == 15
    assert evaluar("aax", a=2, x=5) == 20
    assert evaluar("ab", a=2, ab=7) == 7


def test_variables_usadas():
    usadas = set()
    _parse.analizar("x + 2", Fraction, {"x": Fraction(1), "a": Fraction(2)}, usadas=usadas)
    assert usadas == {"x"}


def test_letras():
    assert _parse.letras("x^2 + 1") == {"x"} and _parse.letras("t + a2") == {"t", "a"}
    assert _parse.letras("12") == set()


# ---------------------------------------------------------------------------
# Errores de sintaxis (con la posición)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("texto, pos, motivo", [
    ("", 0, "expresión vacía"),
    ("   ", 3, "expresión vacía"),
    ("2 +", 3, "se esperaba un número, una variable o '('"),
    ("()", 1, "se esperaba un número, una variable o '('"),
    ("(1 + 2", 6, "falta ')'"),
    ("1 + 2)", 5, "paréntesis ')' sin abrir"),
    ("2 3", 2, "falta un operador antes de '3'"),
    ("x 2", 2, "falta un operador antes de '2'"),
    ("2^3^2", 3, "símbolo inesperado"),
    ("x^y", 2, "el exponente debe ser un número entero"),
    ("x^(2", 4, "falta ')'"),
    ("x^²", 2, "el exponente debe ser un número entero"),
    ("1.5", 1, "carácter inesperado '.' (los números decimales no están admitidos: "
               "usa fracciones como 1/2)"),
    ("1,5", 1, "carácter inesperado ',' (los números decimales no están admitidos: "
               "usa fracciones como 1/2)"),
    ("x # 2", 2, "carácter inesperado '#'"),
    ("x⁻", 1, "superíndice no válido: '⁻'"),
    ("x⁺⁻1", 1, "superíndice no válido: '⁺⁻'"),
    ("2 * z", 4, "variable desconocida 'z' (variables válidas: x, y)"),
    ("xzy", 1, "variable desconocida 'zy' (variables válidas: x, y)"),
])
def test_errores_de_sintaxis(texto, pos, motivo):
    with pytest.raises(ValueError) as info:
        evaluar(texto, x=1, y=2)
    assert str(info.value) == mensaje(texto, pos, motivo)


def test_mensaje_alinea_tabuladores():
    with pytest.raises(ValueError) as info:
        evaluar("\t1 +")
    assert str(info.value) == mensaje(" 1 +", 4, "se esperaba un número, una variable o '('")


def test_variable_desconocida_sin_variables():
    with pytest.raises(ValueError, match=r"variable desconocida 'a' \(variables válidas: ninguna\)"):
        evaluar("a")


def test_anidamiento_excesivo():
    n = _parse.PROFUNDIDAD_MAXIMA + 1
    with pytest.raises(ValueError, match="expresión demasiado anidada"):
        evaluar("(" * n + "1" + ")" * n)


def test_errores_de_calculo_conservan_el_tipo_y_la_posicion():
    with pytest.raises(ZeroDivisionError) as info:
        evaluar("1 + 1/0")
    assert str(info.value) == mensaje("1 + 1/0", 5, "Fraction(1, 0)")
    with pytest.raises(ZeroDivisionError) as info:
        GF(7)("2 + 0^-1")
    assert str(info.value).endswith("2 + 0^-1\n         ^")
    with pytest.raises(ValueError, match="solo se puede dividir por constantes") as info:
        GF(7)["x"]("1/x")
    assert str(info.value).endswith("1/x\n     ^")


# ---------------------------------------------------------------------------
# En cuerpos y anillos
# ---------------------------------------------------------------------------
def test_en_cuerpo_primo():
    F = GF(7)
    assert F("3/2 + 1") == 6 and F("-1") == 6 and F("2^10") == 2 and F("(1+2)(3+4)") == 0
    assert F("3⁻¹") == 5 and F("10^100") == pow(10, 100, 7)


def test_en_extension():
    K = GF(9)
    a = K.gen
    assert K("a^2") == -1 and K("2a + 1") == 2 * a + 1 and K("(a+1)(a-1)") == a * a - 1
    assert K("a**3") == a ** 3 and K("a³") == a ** 3 and K("1/a") == a ** -1
    assert K("−a") == -a and K("a·a") == -1 and K("a × 2") == 2 * a
    L = GF(9, var="alpha")
    assert L("alpha^2 + 1") == 0 and L("2alpha") == 2 * L.gen
    with pytest.raises(ValueError, match="variable desconocida 'x'"):
        K("x")


def test_en_anillo_de_polinomios():
    K = GF(9)
    a = K.gen
    R = K["x"]
    x = R.gen
    assert R("(a+1)x^3 + 2x + a") == (a + 1) * x ** 3 + 2 * x + a
    assert R("3a^2x") == 0 and R("2a^2x") == 2 * a ** 2 * x and R("ax") == a * x
    assert R("xa") == a * x and R("(a+1)(x-1)") == (a + 1) * (x - 1) and R("2 x") == 2 * x
    assert R("x²+x¹⁰") == x ** 2 + x ** 10 and R("x/a") == x * a ** -1 and R("x/2") == 2 * x
    assert R("(x+1)^2 / (a+1)") == (x + 1) ** 2 / (a + 1)
    with pytest.raises(ValueError, match="un polinomio no constante no tiene inverso"):
        R("x^-1")
    with pytest.raises(ValueError, match="demasiado grande"):
        R("x^100000000")
    with pytest.raises(ZeroDivisionError, match="división por el polinomio cero"):
        R("x/(a^2+1)")


def test_variables_de_varias_letras_en_anillos():
    L = GF(9, var="alpha")
    R = L["t"]
    assert R("alphat + t") == (L.gen + 1) * R.gen
    S = GF(7)["xy"]
    assert S("2xy^2 + xy") == 2 * S.gen ** 2 + S.gen


def test_operador_division_por_defecto():
    assert _parse.analizar("6/3", Fraction) == 2
    assert _parse.analizar("6/4", Fraction, None, operator.truediv) == Fraction(3, 2)
