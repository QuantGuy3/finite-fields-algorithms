"""Propiedades algebraicas comprobadas con hypothesis (pocos ejemplos por prueba).

Cuerpos de característica 2 e impar, primos y extensiones, pequeños y grandes
(p de 61 bits y F_{(2^61-1)^2}).  ``derandomize=True``: los ejemplos son siempre
los mismos (la suite es reproducible y no depende de la base de datos de
hypothesis).  Para una campaña más larga, la variable de entorno
``ALCP_EJEMPLOS`` fija el número de ejemplos por prueba (25 por defecto)::

    ALCP_EJEMPLOS=1000 python build.py pytest -m propiedades
"""
from __future__ import annotations

import os
import pickle

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

import alcp
from alcp import GF

pytestmark = pytest.mark.propiedades

M61 = 2**61 - 1
CUERPOS = [2, 4, 8, 256, 3, 9, 27, 7, 25, M61, M61**2]
PEQUENOS = [2, 4, 8, 3, 9, 7, 25]

POCOS = settings(max_examples=int(os.environ.get("ALCP_EJEMPLOS", "25")), deadline=None,
                 derandomize=True,
                 suppress_health_check=[HealthCheck.too_slow])

cuerpos = st.sampled_from(CUERPOS).map(GF)
cuerpos_pequenos = st.sampled_from(PEQUENOS).map(GF)


def elemento(K):
    return st.integers(min_value=0, max_value=K.order - 1).map(K.from_index)


def polinomio(K, grado_maximo=6):
    return st.lists(elemento(K), max_size=grado_maximo + 1).map(K["x"])


@POCOS
@given(st.data())
def test_axiomas_de_cuerpo_y_frobenius(datos):
    K = datos.draw(cuerpos)
    x, y, z = (datos.draw(elemento(K)) for _ in range(3))
    assert (x + y) + z == x + (y + z) and (x * y) * z == x * (y * z)
    assert x + y == y + x and x * y == y * x and x * (y + z) == x * y + x * z
    assert x - x == 0 and x + (-x) == 0 and x * 1 == x and x + 0 == x
    if y:
        assert (x / y) * y == x and y * y.inverse() == 1 and y ** -3 == (y**3).inverse()
    p = K.p
    assert (x + y) ** p == x**p + y**p and (x * y).frobenius() == x.frobenius() * y.frobenius()
    assert x ** K.order == x and x.frobenius(K.degree) == x
    assert (x + y).trace() == x.trace() + y.trace() and (x * y).norm() == x.norm() * y.norm()
    if x.is_square():
        r = x.sqrt()
        assert r * r == x and r.index <= (-r).index and x.sqrt(all=True)[0] == r
    else:
        assert x.sqrt(all=True) == [] and K.p != 2


@POCOS
@given(st.data())
def test_anillo_de_polinomios(datos):
    K = datos.draw(cuerpos)
    f, g, h = (datos.draw(polinomio(K)) for _ in range(3))
    assert (f + g) * h == f * h + g * h and (f * g) * h == f * (g * h) and f * g == g * f
    assert (f * g).derivative() == f.derivative() * g + f * g.derivative()
    if g:
        q, r = divmod(f, g)
        assert q * g + r == f and r.degree < g.degree and (f * g) / g == f
        assert alcp.divmod_toeplitz(f, g) == (q, r)
    d, s, t = f.xgcd(g)
    assert d == s * f + t * g and d == f.gcd(g)
    if d:
        assert d.is_monic() and f % d == 0 and g % d == 0
    v = datos.draw(elemento(K))
    assert (f * g)(v) == f(v) * g(v) and f(g)(v) == f(g(v)) and f.compose(g) == f(g)
    if f.degree <= 5 and K.order < 10**6:
        F = f.factor() if f else None
        if F is not None:
            assert F.expand() == f and all(u.is_monic() and u.is_irreducible() for u, _ in F)
            assert f.radical() == (f.monic() if f.degree < 1 else
                                   alcp.Factorization([(u, 1) for u, _ in F], 1, f.ring).expand())


@POCOS
@given(st.data())
def test_raices_contra_evaluacion(datos):
    K = datos.draw(cuerpos_pequenos)
    f = datos.draw(polinomio(K, 5))
    if not f:
        return
    raices = f.roots(multiplicities=False)
    assert raices == [c for c in K if f(c) == 0] or f.degree < 1
    x = f.ring.gen
    for c, m in f.roots():
        assert f % (x - c) ** m == 0 and f % (x - c) ** (m + 1) != 0


NOMBRES = [("a", "x"), ("alpha", "t"), ("ab", "a"), ("b", "ba"), ("x", "y"), ("A", "a")]


@POCOS
@given(st.data())
def test_textos_se_vuelven_a_leer(datos):
    q = datos.draw(st.sampled_from([2, 4, 9, 27, 25, M61, M61**2]))
    var, rvar = datos.draw(st.sampled_from(NOMBRES))
    K = GF(q, var=var)
    R = K[rvar]
    f = R(datos.draw(st.lists(elemento(K), max_size=7)))
    for estilo in ("unicode", "ascii"):
        assert R(format(f, estilo)) == f
    e = datos.draw(elemento(K))
    assert K(str(e)) == e and K(format(e, "ascii")) == e
    assert eval(repr(f), {"GF": GF}) == f and eval(repr(e), {"GF": GF}) == e
    assert pickle.loads(pickle.dumps(f)) == f and pickle.loads(pickle.dumps(e)) == e


@POCOS
@given(st.data())
def test_igualdad_y_hash_coherentes_entre_cuerpos_y_anillos(datos):
    """Un polinomio sobre F_p y su copia sobre F_q (y sus constantes, elementos y
    enteros del rango canónico) son iguales y tienen el mismo hash."""
    p, q = datos.draw(st.sampled_from([(2, 4), (2, 8), (3, 9), (3, 27), (5, 25)]))
    Fp, Fq = GF(p), GF(q)
    coefs = datos.draw(st.lists(st.integers(0, p - 1), max_size=6))
    f, g = Fp["x"](coefs), Fq["x"](coefs)
    assert f == g and g == f and hash(f) == hash(g)
    c = datos.draw(st.integers(0, p - 1))
    iguales = [Fp(c), Fq(c), Fp["x"](c), Fq["x"](c), c]
    for a in iguales:
        for b in iguales:
            assert a == b and hash(a) == hash(b)
    assert len(set(iguales)) == 1


# ---------------------------------------------------------------------------
# Fuzzing: entradas arbitrarias nunca producen excepciones inesperadas
# ---------------------------------------------------------------------------
_K_AB = GF(9, var="ab")
DESTINOS_TEXTO = [GF(7), GF(9), GF(7)["x"], GF(9)["x"], _K_AB["a"], _K_AB["b"], GF(2)["t"]]
PIEZAS = ["x", "a", "b", "ab", "t", "y", "2", "0", "7", "10", "+", "-", "−", "*", "**", "^", "/",
          "÷", "·", "(", ")", " ", "²", "⁻", "¹", "⁰", ".", ",", "#", "\t", "99999999999", "x^", "(-"]


@POCOS
@given(st.sampled_from(DESTINOS_TEXTO),
       st.one_of(st.lists(st.sampled_from(PIEZAS), max_size=14).map("".join), st.text(max_size=12)))
def test_analizador_con_cadenas_arbitrarias(destino, texto):
    """El analizador devuelve un valor o lanza ValueError/ZeroDivisionError con la
    expresión y un ^ dentro de ella; nunca otra excepción."""
    try:
        destino(texto)
    except (ValueError, ZeroDivisionError) as exc:
        motivo, expresion, marca = str(exc).split("\n")
        assert expresion.startswith("    ") and marca.strip() == "^"
        assert len(marca) - 5 <= len(texto)


_F7, _K9 = GF(7), GF(9)
_escalares = st.one_of(
    st.none(), st.booleans(), st.integers(), st.floats(), st.complex_numbers(), st.fractions(),
    st.decimals(), st.text(max_size=6), st.binary(max_size=4),
    st.sampled_from([_F7(3), _K9.gen, _F7["x"].gen, _K9["x"]("x+a"), GF(5)(1), GF(9, var="b").gen]))
_valores = st.recursive(_escalares, lambda hijos: st.one_of(
    st.lists(hijos, max_size=4), st.tuples(hijos, hijos),
    st.dictionaries(st.integers(-3, 12), hijos, max_size=3), st.frozensets(st.integers(), max_size=3)),
    max_leaves=8)


@POCOS
@given(st.sampled_from([_F7, _K9, GF(2), _F7["x"], _K9["x"], GF(4)["t"]]), _valores)
def test_conversiones_con_valores_arbitrarios(destino, valor):
    """K(x)/R(x) solo lanzan TypeError, ValueError o ZeroDivisionError; ``x in K``
    coincide con que la conversión funcione; == devuelve siempre un bool."""
    try:
        destino(valor)
        convertible = True
    except (TypeError, ValueError, ZeroDivisionError):
        convertible = False
    assert (valor in destino) is convertible
    e = destino.gen
    igual = e == valor
    assert igual is True or igual is False
    for op in ("__add__", "__mul__", "__sub__", "__truediv__", "__radd__", "__rsub__"):
        try:
            r = getattr(e, op)(valor)
        except (TypeError, ValueError, ZeroDivisionError):
            continue
        assert r is NotImplemented or isinstance(r, (alcp.FieldElement, alcp.Polynomial))


@POCOS
@given(st.data())
def test_algoritmos_rapidos(datos):
    q = datos.draw(st.sampled_from([17, 97, 9, 25, M61**2]))
    K = GF(q)
    k = datos.draw(st.integers(0, 3))
    n = 2**k
    if (q - 1) % n:
        n = 1
    a = [datos.draw(elemento(K)) for _ in range(n)]
    assert alcp.idft(K, alcp.dft(K, a)) == a
    # T·b con T triangular inferior de primera columna a = los n primeros coeficientes de a(x)·b(x).
    m = datos.draw(st.integers(1, 5))
    a = [datos.draw(elemento(K)) for _ in range(m)]
    b = [datos.draw(elemento(K)) for _ in range(m)]
    producto = (K["x"](a) * K["x"](b))
    assert alcp.toeplitz_lower_vec(K, a, b) == [producto[i] for i in range(m)]
    if a[0]:
        inversa = alcp.toeplitz_lower_inverse(K, a)
        assert [(K["x"](a) * K["x"](inversa))[i] for i in range(m)] == [K(1)] + [K(0)] * (m - 1)
