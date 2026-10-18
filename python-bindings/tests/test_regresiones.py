"""Regresiones: un test por cada defecto encontrado en la revisión adversarial.

Cada prueba fallaba (o no terminaba) antes de su corrección:

1. ``sqrt``, ``root_of_unity``, ``primitive_element`` y ``dft`` sin raíz no
   terminaban en una extensión F_{p^n} con p grande (p. ej. p = 2^61 - 1): la
   búsqueda recorría los índices 1, 2, 3... que son los elementos de F_p, y
   ninguno servía (en grado par todos son cuadrados; sus potencias no alcanzan
   órdenes que no dividen a p - 1).  Corrección (alcp/field.py): se saltan los
   elementos de F_p cuando se sabe que ninguno sirve; el resultado sigue siendo
   el de menor índice.
2. El analizador no admitía literales enteros de más de 4300 cifras: el
   ``ValueError`` de CPython (límite de conversión de cadenas a int) salía en
   inglés y sin la posición.  Corrección (alcp/_parse.py): conversión por
   bloques, sin límite, como el resto de la librería.
3. ``bool(K)`` lanzaba ``OverflowError`` en un cuerpo de más de sys.maxsize
   elementos (Python recurría a ``__len__``).  Corrección (alcp/field.py):
   ``FiniteField.__bool__`` devuelve siempre True.
4. Los mensajes de error que citan un entero del usuario fallaban si este
   tenía más de 4300 cifras (``str(int)`` de CPython lanza ValueError): el
   mensaje salía en inglés y ``K(Fraction(1, 7·10^5000))`` en GF(7) lanzaba
   ValueError en lugar de ZeroDivisionError.  Corrección (alcp/_util.py,
   field.py, poly.py): ``texto_entero`` resume los enteros enormes.
5. Comparar un cuerpo, un anillo o una factorización con un escalar o un array
   de numpy lanzaba ValueError o devolvía un array.  Corrección (alcp/field.py,
   poly.py): ``__array_ufunc__ = None`` también en esas clases.
6. En el módulo nativo, ``pow_mod`` con k < 0 y módulo cero lanzaba ValueError
   en lugar del ZeroDivisionError documentado.  Corrección (bindings.cpp).
"""
from __future__ import annotations

from fractions import Fraction

import pytest

import alcp
from alcp import GF

pytestmark = pytest.mark.regresion

M61 = 2**61 - 1
M127 = 2**127 - 1


# ---------------------------------------------------------------------------
# 1. Búsquedas que no terminaban en extensiones con p grande
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("q", [M61**2, M127**2, M61**4], ids=["M61^2", "M127^2", "M61^4"])
def test_raiz_cuadrada_en_extension_con_p_grande(q, potencias_limitadas):
    K = GF(q)
    a = K.gen
    with potencias_limitadas(K, 2000):
        cuadrado = (a + 3) ** 2
        raices = cuadrado.sqrt(all=True)
        assert raices == sorted(raices, key=lambda r: r.index) and len(raices) == 2
        assert {r ** 2 for r in raices} == {cuadrado} and a + 3 in raices
        assert cuadrado.sqrt() == raices[0]
        # Un no cuadrado: el primer elemento de índice >= p que lo es.
        no_cuadrado = next(K.from_index(i) for i in range(K.p, K.p + 50)
                           if not K.from_index(i).is_square())
        assert no_cuadrado.sqrt(all=True) == []
        with pytest.raises(ValueError, match="no es un cuadrado"):
            no_cuadrado.sqrt()


def test_raices_de_la_unidad_en_extension_con_p_grande(potencias_limitadas):
    K = GF(M61**2)                       # q - 1 = (p - 1)(p + 1), p + 1 = 2^61
    with potencias_limitadas(K, 2000):
        for n in (4, 2**62, M61 - 1, M61 + 1, 3 * 2**40):
            w = K.root_of_unity(n)
            assert w ** n == 1 and w.multiplicative_order() == n
        g = K.primitive_element()
    assert g.multiplicative_order() == K.order - 1
    assert g is K.primitive_element()    # se calcula una sola vez


def test_dft_sin_raiz_en_extension_con_p_grande(potencias_limitadas):
    K = GF(M61**2)
    with potencias_limitadas(K, 200):
        w = K.root_of_unity(8)
    datos = [K.gen, 1, 2, "a + 5", 0, 0, 7, K.gen ** 3]
    # Con la búsqueda acotada, dft sin raíz usa la misma raíz w.
    assert alcp.dft(K, datos) == alcp.dft(K, datos, w)
    assert alcp.idft(K, alcp.dft(K, datos)) == [K(v) for v in datos]


@pytest.mark.parametrize("q", [9, 25, 49, 81, 121, 169, 3**5, 2**6])
def test_raices_de_la_unidad_siguen_siendo_las_de_menor_indice(q):
    """La corrección no cambia el resultado: sigue siendo x^((q-1)/n) para el
    x de menor índice con el que tiene orden n (fuerza bruta)."""
    K = GF(q)
    for n in (d for d in range(1, q) if (q - 1) % d == 0):
        e = (q - 1) // n
        esperado = next(K.from_index(i) ** e for i in range(1, q)
                        if (K.from_index(i) ** e).multiplicative_order() == n)
        assert K.root_of_unity(n) == esperado, n
    # Raíces cuadradas: la de menor índice (fuerza bruta).
    for x in list(K)[:40]:
        bruta = [r for r in K if r * r == x]
        assert x.sqrt(all=True) == bruta


# ---------------------------------------------------------------------------
# 2. Literales enteros de más de 4300 cifras
# ---------------------------------------------------------------------------
def test_literales_enteros_de_mas_de_4300_cifras():
    F, R = GF(7), GF(7)["x"]
    N = (10**5000 - 1) // 9                     # 5000 unos
    assert F("1" * 5000) == F(N) and F("-" + "1" * 5000) == F(-N)
    assert R("1" * 5000 + "x + 2") == R([2, N])
    assert GF(9)("1" * 4301 + "a") == GF(9)([0, int("1" * 4300) * 10 + 1])
    # Superíndice de más de 4300 cifras (exponente enorme, reducido en el cuerpo).
    assert F("3" + "¹" * 4400) == F(pow(3, (10**4400 - 1) // 9, 7))
    # Un exponente enorme en un anillo sigue dando el error de grado, con posición.
    with pytest.raises(ValueError, match=r"demasiado grande[\s\S]*\^"):
        R("x^" + "9" * 5000)


# ---------------------------------------------------------------------------
# 3. bool de un cuerpo enorme
# ---------------------------------------------------------------------------
def test_bool_de_un_cuerpo_enorme():
    for K in (GF(2**127), GF(M127), GF(7)):
        assert bool(K) is True and (K or None) is K
    with pytest.raises(OverflowError, match="usa K.order"):
        len(GF(2**127))


# ---------------------------------------------------------------------------
# 4. Mensajes de error con enteros enormes
# ---------------------------------------------------------------------------
def test_mensajes_de_error_con_enteros_de_mas_de_4300_cifras():
    enorme = 10**5000                                # 16610 bits
    F, R = GF(7), GF(7)["x"]
    x = R.gen
    with pytest.raises(ZeroDivisionError, match=r"denominador de 1/<entero de 16613 bits> es "
                                                r"múltiplo de p = 7"):
        F(Fraction(1, 7 * enorme))
    casos = [
        (lambda: x ** enorme, "grado <entero de 16610 bits> demasiado grande"),
        (lambda: R.monomial(enorme), "demasiado grande"),
        (lambda: R.random(enorme), "demasiado grande"),
        (lambda: R.irreducible(enorme), "demasiado grande"),
        (lambda: R.monomial(-enorme), "monomio debe ser >= 0, no -<entero de 16610 bits>"),
        (lambda: R.irreducible(-enorme), "grado >= 1, no -<entero"),
        (lambda: R({-enorme: 1}), "exponente negativo en un polinomio: -<entero"),
        (lambda: (x + 1).equal_degree_factorization(-enorme), "r debe ser >= 1, no -<entero"),
        (lambda: F.from_index(enorme), r"\[0, 7\), no <entero de 16610 bits>"),
        (lambda: F.root_of_unity(enorme), "orden <entero de 16610 bits>: <entero"),
        (lambda: alcp.prime_power_factor(enorme), "<entero de 16610 bits> no es potencia"),
        (lambda: GF(enorme + 1), "no es potencia de un número primo"),
        (lambda: R("x^" + "9" * 5000),
         r"grado <entero de 16610 bits> demasiado grande \(máximo 10000000\)\n    x\^9{5000}\n"
         r"     \^$"),
    ]
    for funcion, mensaje in casos:
        with pytest.raises(ValueError, match=mensaje):
            funcion()
    # Los enteros «normales» se siguen escribiendo completos.
    with pytest.raises(ValueError, match=f"no {10**300}$"):
        F.from_index(10**300)


# ---------------------------------------------------------------------------
# 5. Comparación con escalares y arrays de numpy
# ---------------------------------------------------------------------------
def test_comparacion_de_cuerpos_anillos_y_factorizaciones_con_numpy():
    """Antes, ``factorización == numpy.float64(1)`` lanzaba ValueError (numpy
    intentaba convertir la factorización, una secuencia, en un array) y
    ``GF(9) == numpy.int64(1)`` devolvía un array (comparaba con cada elemento).
    Corrección: ``__array_ufunc__ = None`` en FiniteField, PolynomialRing y
    Factorization, como ya tenían los elementos y los polinomios."""
    np = pytest.importorskip("numpy")
    K = GF(9)
    R = K["x"]
    objetos = [K, GF(2**127), R, R("(x+a)(x+1)").factor(), R(2).factor()]
    for o in objetos:
        for r in (np.float64(1.0), np.int64(1), np.array([1.0, 2.0])):
            assert (o == r) is False and (r == o) is False, (o, r)
            assert (o != r) is True and (r != o) is True, (o, r)


# ---------------------------------------------------------------------------
# 6. Módulo nativo: pow_mod con exponente negativo y módulo cero
# ---------------------------------------------------------------------------
def test_pow_mod_nativo_con_modulo_cero():
    """La documentación de PolyFp/PolyFq.pow_mod promete ZeroDivisionError si el
    módulo es cero; con k < 0 se lanzaba ValueError («not invertible»), porque
    el núcleo intentaba antes invertir módulo 0.  Corrección (bindings.cpp): se
    comprueba el módulo cero antes de llamar al núcleo."""
    al = alcp.alcp_native
    for K in (al.Fp(7), al.Fq(3, [1, 0, 1])):
        f, cero = K.poly([1, 2]), K.poly([])
        for k in (-3, -1, 0, 3):
            with pytest.raises(ZeroDivisionError):
                f.pow_mod(k, cero)
        assert f.pow_mod(-1, K.poly([0, 0, 1])) * f % K.poly([0, 0, 1]) == K.poly([1])
