"""Pruebas de la representación (alcp._pretty): estilos, formato, LaTeX, IPython y repr.

Cada caso fija el texto exacto en los tres estilos ("unicode", "ascii" y
"latex") y comprueba que ``repr`` es evaluable y que la salida "ascii" la
vuelve a leer el analizador.
"""
from __future__ import annotations

import pytest

import alcp
from alcp import Factorization

GF = alcp.GF


@pytest.fixture(autouse=True)
def estilo_por_defecto():
    alcp.set_print_style("unicode")
    yield
    alcp.set_print_style("unicode")


def casos():
    F, K = GF(7), GF(9)
    a = K.gen
    L = GF(9, var="alpha")
    Rx, Kx = F["x"], K["x"]
    return [
        # (objeto, unicode, ascii, latex, repr)
        (F, "GF(7)", "GF(7)", r"\mathbb{F}_{7}", "GF(7)"),
        (K, "GF(9) = GF(3)[a]/(a² + 1)", "GF(9) = GF(3)[a]/(a^2 + 1)",
         r"\mathbb{F}_{9} = \mathbb{F}_{3}[a]/(a^{2} + 1)", "GF(9, modulus='x^2 + 1')"),
        (L, "GF(9) = GF(3)[alpha]/(alpha² + 1)", "GF(9) = GF(3)[alpha]/(alpha^2 + 1)",
         r"\mathbb{F}_{9} = \mathbb{F}_{3}[\alpha]/(\alpha^{2} + 1)",
         "GF(9, modulus='x^2 + 1', var='alpha')"),
        (GF(9, var="gen"), "GF(9) = GF(3)[gen]/(gen² + 1)", "GF(9) = GF(3)[gen]/(gen^2 + 1)",
         r"\mathbb{F}_{9} = \mathbb{F}_{3}[\mathit{gen}]/(\mathit{gen}^{2} + 1)",
         "GF(9, modulus='x^2 + 1', var='gen')"),
        (F(3), "3", "3", "3", "GF(7)(3)"),
        (a + 1, "a + 1", "a + 1", "a + 1", "GF(9, modulus='x^2 + 1')('a + 1')"),
        (2 * a, "2a", "2*a", "2 a", "GF(9, modulus='x^2 + 1')('2*a')"),
        (K(0), "0", "0", "0", "GF(9, modulus='x^2 + 1')('0')"),
        (GF(27)("a^2 + 2a"), "a² + 2a", "a^2 + 2*a", "a^{2} + 2 a",
         "GF(27, modulus='x^3 + 2*x + 1')('a^2 + 2*a')"),
        (Rx("x^3 + 2x + 1"), "x³ + 2x + 1", "x^3 + 2*x + 1", "x^{3} + 2 x + 1",
         "GF(7)['x']('x^3 + 2*x + 1')"),
        (Rx(0), "0", "0", "0", "GF(7)['x']('0')"),
        (Rx("x^10 - x^2 + 5"), "x¹⁰ + 6x² + 5", "x^10 + 6*x^2 + 5", "x^{10} + 6 x^{2} + 5",
         "GF(7)['x']('x^10 + 6*x^2 + 5')"),
        (Kx("(a+1)x^3 + 2x + a"), "(a + 1)x³ + 2x + a", "(a + 1)*x^3 + 2*x + a",
         r"\left(a + 1\right) x^{3} + 2 x + a",
         "GF(9, modulus='x^2 + 1')['x']('(a + 1)*x^3 + 2*x + a')"),
        (Kx("a x^2 + x"), "ax² + x", "a*x^2 + x", "a x^{2} + x",
         "GF(9, modulus='x^2 + 1')['x']('a*x^2 + x')"),
        (F["xy"]("2xy^2 + xy"), "2·xy² + xy", "2*xy^2 + xy", r"2 \mathit{xy}^{2} + \mathit{xy}",
         "GF(7)['xy']('2*xy^2 + xy')"),
        (L["x"]("alpha x^2 + (alpha+1)x + alpha^3"), "alpha·x² + (alpha + 1)·x + 2·alpha",
         "alpha*x^2 + (alpha + 1)*x + 2*alpha",
         r"\alpha x^{2} + \left(\alpha + 1\right) x + 2 \alpha",
         "GF(9, modulus='x^2 + 1', var='alpha')['x']('alpha*x^2 + (alpha + 1)*x + 2*alpha')"),
        (Rx, "GF(7)[x]", "GF(7)[x]", r"\mathbb{F}_{7}[x]", "GF(7)[x]"),
        (K["t"], "GF(9)[t]", "GF(9)[t]", r"\mathbb{F}_{9}[t]", "GF(9, modulus='x^2 + 1')[t]"),
        (Rx("3(x+1)^2(x^2+1)(x+5)").factor(), "3 · (x + 1)² · (x + 5) · (x² + 1)",
         "3 * (x + 1)^2 * (x + 5) * (x^2 + 1)",
         r"3 \cdot \left(x + 1\right)^{2} \cdot \left(x + 5\right) \cdot \left(x^{2} + 1\right)",
         "Factorization([(GF(7)['x']('x + 1'), 2), (GF(7)['x']('x + 5'), 1), "
         "(GF(7)['x']('x^2 + 1'), 1)], unit=GF(7)(3))"),
        (Rx("x + 1").factor(), "x + 1", "x + 1", "x + 1",
         "Factorization([(GF(7)['x']('x + 1'), 1)], unit=GF(7)(1))"),
        (Rx("(x+1)^2").factor(), "(x + 1)²", "(x + 1)^2", r"\left(x + 1\right)^{2}",
         "Factorization([(GF(7)['x']('x + 1'), 2)], unit=GF(7)(1))"),
        (Rx("3x").factor(), "3 · x", "3 * x", r"3 \cdot x",
         "Factorization([(GF(7)['x']('x'), 1)], unit=GF(7)(3))"),
        (Rx("x^2").factor(), "x²", "x^2", "x^{2}",
         "Factorization([(GF(7)['x']('x'), 2)], unit=GF(7)(1))"),
        (Rx(5).factor(), "5", "5", "5",
         "Factorization([], unit=GF(7)(5), ring=GF(7)['x'])"),
        (Kx("(a+1)(x+1)(x^2+x+a)").factor(), "(a + 1) · (x + 1) · (x² + x + a)",
         "(a + 1) * (x + 1) * (x^2 + x + a)",
         r"\left(a + 1\right) \cdot \left(x + 1\right) \cdot \left(x^{2} + x + a\right)",
         "Factorization([(GF(9, modulus='x^2 + 1')['x']('x + 1'), 1), "
         "(GF(9, modulus='x^2 + 1')['x']('x^2 + x + a'), 1)], "
         "unit=GF(9, modulus='x^2 + 1')('a + 1'))"),
        (Factorization([], a + 1), "a + 1", "a + 1", "a + 1",
         "Factorization([], unit=GF(9, modulus='x^2 + 1')('a + 1'), "
         "ring=GF(9, modulus='x^2 + 1')['x'])"),
    ]


@pytest.mark.parametrize("i", range(len(casos())))
def test_textos_exactos(i):
    obj, uni, asc, lat, rep = casos()[i]
    assert str(obj) == uni and f"{obj}" == uni and format(obj, "u") == uni
    assert format(obj, "unicode") == uni
    assert format(obj, "ascii") == asc and f"{obj:a}" == asc
    assert format(obj, "latex") == lat and f"{obj:l}" == lat
    assert obj._repr_latex_() == f"${lat}$"
    assert repr(obj) == rep


class Impresora:
    """Imita al «pretty printer» de IPython (solo necesita .text)."""

    def __init__(self):
        self.partes = []

    def text(self, s):
        self.partes.append(s)


def test_repr_pretty_de_ipython():
    for obj, uni, asc, *_ in casos():
        p = Impresora()
        obj._repr_pretty_(p, False)
        assert p.partes == [uni]
        with alcp.print_style("ascii"):
            p = Impresora()
            obj._repr_pretty_(p, True)
            assert p.partes == [asc]


def test_repr_evaluable():
    espacio = {"GF": GF, "Factorization": Factorization}
    for obj, *_ in casos():
        if isinstance(obj, alcp.PolynomialRing):
            continue                     # su repr es «GF(7)[x]», como pide el diseño
        assert eval(repr(obj), espacio) == obj


def test_ascii_se_vuelve_a_leer():
    for obj, _, asc, *_ in casos():
        if isinstance(obj, alcp.Polynomial):
            assert obj.ring(asc) == obj
        elif isinstance(obj, alcp.FieldElement):
            assert obj.field(asc) == obj
        elif isinstance(obj, Factorization):
            assert obj.ring(asc) == obj.expand()


def test_estilo_global():
    f = GF(7)["x"]("x^2 + 1")
    assert alcp.get_print_style() == "unicode"
    alcp.set_print_style("ascii")
    assert alcp.get_print_style() == "ascii" and str(f) == "x^2 + 1" and f"{f}" == "x^2 + 1"
    alcp.set_print_style("l")
    assert alcp.get_print_style() == "latex" and str(f) == "x^{2} + 1"
    alcp.set_print_style("u")
    assert str(f) == "x² + 1"
    for malo in ("negrita", "", 3, None):
        with pytest.raises(ValueError, match="estilo de impresión desconocido"):
            alcp.set_print_style(malo)
    assert alcp.get_print_style() == "unicode"


def test_gestor_de_contexto():
    f = GF(7)["x"]("x^2 + 1")
    with alcp.print_style("ascii") as estilo:
        assert estilo == "ascii" and str(f) == "x^2 + 1"
        with alcp.print_style("latex"):
            assert str(f) == "x^{2} + 1"
        assert str(f) == "x^2 + 1"
    assert str(f) == "x² + 1"
    with pytest.raises(RuntimeError):
        with alcp.print_style("a"):
            raise RuntimeError("se restaura igualmente")
    assert alcp.get_print_style() == "unicode"
    with pytest.raises(ValueError, match="estilo de impresión desconocido"):
        with alcp.print_style("raro"):
            pass
    assert alcp.get_print_style() == "unicode"


def test_especificadores_de_formato_genericos():
    f = GF(7)["x"]("x + 1")
    assert f"{f:>8}" == "   x + 1" and f"{f:*^9}" == "**x + 1**"
    with alcp.print_style("ascii"):
        assert f"{GF(7)['x']('x^2'):<6}|" == "x^2   |"
    e = GF(7)(10)
    assert f"{e:03d}" == "003" and f"{e:b}" == "11" and f"{e:x}" == "3" and f"{e:X}" == "3"
    assert f"{e:o}" == "3" and f"{e:n}" == "3" and f"{e:>4}" == "   3" and f"{e:a}" == "3"
    a = GF(9).gen
    assert f"{a + 1:>7}" == "  a + 1"
    with pytest.raises(ValueError):
        format(a, "d")                                   # no es un número
