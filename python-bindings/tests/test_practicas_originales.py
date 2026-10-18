"""Cruce del paquete ``alcp`` con las prácticas ORIGINALES de ``python/``.

Las prácticas (``cuerpos_finitos``, ``factorizacion``, ``algoritmos_rapidos``)
solo se importan, sin escribir ``.pyc`` (fixture ``practicas`` de conftest.py);
no se modifican.  Se compara únicamente donde la referencia es correcta:

* aritmética de elementos, de polinomios, divmod, gcd, Euclides extendido,
  potencia e inverso modulares e irreducibilidad, en F_p y en F_q;
* factorización de polinomios MÓNICOS con p impar (la referencia falla con
  p = 2 y con polinomios no mónicos sobre F_q);
* algoritmos rápidos (Toeplitz, inversas, FFT, división, Karatsuba); en F_q,
  ``fq_ifft`` solo con n | p - 1, que es su precondición documentada (con
  n >= p calcula n^-1 como el elemento de ÍNDICE n, no como n·1).

Los cuerpos de alcp se construyen con el mismo módulo que la referencia.
"""
from __future__ import annotations

import itertools
import random
import sys
from pathlib import Path

import pytest

import alcp
from alcp import GF

pytestmark = pytest.mark.referencia


class Puente:
    """Conversión entre los objetos de alcp y los de la práctica original."""

    def __init__(self, practicas, q):
        cf = practicas.cf
        self.K = GF(q)
        self.R = self.K["x"]
        self.p = self.K.p
        self.fp = cf.cuerpo_fp(self.p)
        if self.K.is_prime_field:
            self.Q, self.anillo = self.fp, cf.anillo_fp_x(self.fp)
        else:
            self.Q = cf.cuerpo_fq(self.fp, [c.index for c in self.K.modulus.coeffs], var="a")
            self.anillo = cf.anillo_fq_x(self.Q)

    # elementos
    def a_ref(self, e):
        e = self.K(e)
        return self.fp.elem(e.index) if self.K.is_prime_field else self.Q.elem_de_int(e.index)

    def de_ref(self, e):
        return self.K(int(e)) if self.K.is_prime_field else self.K.from_index(self.Q.conv_a_int(e))

    # polinomios
    def pol_a_ref(self, f):
        return self.anillo.elem_de_tuple(tuple(self.a_ref(c) for c in f.coeffs))

    def pol_de_ref(self, g):
        return self.R([self.de_ref(c) for c in self.anillo.conv_a_tuple(g)])


CUERPOS_IMPARES = [3, 7, 13, 9, 25, 27, 49]


@pytest.mark.parametrize("q", [3, 7, 13, 9, 25])
def test_aritmetica_de_elementos_exhaustiva(practicas, q):
    P = Puente(practicas, q)
    Q, K = P.Q, P.K
    for x, y in itertools.product(list(K), repeat=2):
        rx, ry = P.a_ref(x), P.a_ref(y)
        assert P.de_ref(Q.suma(rx, ry)) == x + y and P.de_ref(Q.mult(rx, ry)) == x * y
        if y:
            assert P.de_ref(Q.inv_mult(ry)) == y.inverse()
    for x in K:
        assert P.de_ref(Q.inv_adit(P.a_ref(x))) == -x
        assert P.de_ref(Q.pot(P.a_ref(x), 5)) == x**5


@pytest.mark.parametrize("q", CUERPOS_IMPARES + [2, 4, 8])
def test_polinomios_contra_la_referencia(practicas, q):
    P = Puente(practicas, q)
    A, R = P.anillo, P.R
    alcp.set_seed(q)
    rng = random.Random(q)
    m = R.irreducible(3)
    for _ in range(8):
        f, g = R.random(rng.randint(0, 6)), R.random(rng.randint(0, 4))
        rf, rg = P.pol_a_ref(f), P.pol_a_ref(g)
        assert P.pol_de_ref(A.suma(rf, rg)) == f + g and P.pol_de_ref(A.mult(rf, rg)) == f * g
        if g:
            qq, rr = A.divmod(rf, rg)
            assert (P.pol_de_ref(qq), P.pol_de_ref(rr)) == divmod(f, g)
            assert P.pol_de_ref(A.gcd(rf, rg)) == f.gcd(g)
            d, s, t = (P.pol_de_ref(v) for v in A.gcd_ext(rf, rg))
            assert d == s * f + t * g and d.monic() == f.gcd(g)
        assert A.es_irreducible(rg) == g.is_irreducible()
        if f % m:
            rm = P.pol_a_ref(m)
            assert P.pol_de_ref(A.pot_mod(rf, 7, rm)) == pow(f, 7, m)
            assert P.pol_de_ref(A.inv_mod(rf, rm)) == f.inverse_mod(m)


@pytest.mark.parametrize("q", CUERPOS_IMPARES)
def test_factorizacion_monica_contra_la_referencia(practicas, q):
    P = Puente(practicas, q)
    R = P.R
    alcp.set_seed(q + 1)
    rng = random.Random(q + 1)
    casos = [R.random(rng.randint(1, 6)).monic() for _ in range(6)]
    x = R.gen
    casos += [(x + 1) ** 3 * (x**2 + 2), x**P.p - x, (x**2 + x + 2) ** 2 * x]
    for f in casos:
        f = f.monic()
        if f.degree < 1:
            continue
        referencia = {(P.pol_de_ref(u), e) for u, e in P.anillo.factorizar(P.pol_a_ref(f))}
        assert referencia == set(f.factor()), f


@pytest.mark.parametrize("q", [17, 97, 9, 25])
def test_algoritmos_rapidos_contra_la_referencia(practicas, q):
    ar = practicas.ar
    P = Puente(practicas, q)
    K, R = P.K, P.R
    primo = K.is_prime_field
    tupla = lambda v: tuple(P.a_ref(c) for c in v)          # noqa: E731
    lista = lambda t: [P.de_ref(c) for c in t]              # noqa: E731
    fn = {nombre: getattr(ar, ("fp_" if primo else "fq_") + nombre)
          for nombre in ("toep_inf_vec", "toep_sup_vec", "toep_vec", "toep_inf_inv",
                         "toep_sup_inv", "fft", "ifft", "x_divmod", "x_mult_karatsuba")}
    alcp.set_seed(q + 2)
    rng = random.Random(q + 2)
    Q = P.fp if primo else P.Q
    for n in (1, 2, 3, 5):
        a = [K.random() for _ in range(2 * n - 1)]
        b = [K.random() for _ in range(n)]
        a[0] = a[0] or K.one
        assert lista(fn["toep_inf_vec"](Q, n, tupla(a[:n]), tupla(b))) == alcp.toeplitz_lower_vec(K, a[:n], b)
        assert lista(fn["toep_sup_vec"](Q, n, tupla(a[:n]), tupla(b))) == alcp.toeplitz_upper_vec(K, a[:n], b)
        assert lista(fn["toep_vec"](Q, n, tupla(a), tupla(b))) == alcp.toeplitz_vec(K, a, b)
        assert lista(fn["toep_inf_inv"](Q, n, tupla(a[:n]))) == alcp.toeplitz_lower_inverse(K, a[:n])
        assert lista(fn["toep_sup_inv"](Q, n, tupla(a[:n]))) == alcp.toeplitz_upper_inverse(K, a[:n])
    k = 0
    while (q - 1) % 2**k == 0 and k <= 4:
        n = 2**k
        w = K.root_of_unity(n)
        a = [K.random() for _ in range(n)]
        assert lista(fn["fft"](Q, P.a_ref(w), k, tupla(a))) == alcp.dft(K, a, w)
        if (K.p - 1) % n == 0:
            assert lista(fn["ifft"](Q, P.a_ref(w), k, tupla(a))) == alcp.idft(K, a, w)
        k += 1
    for _ in range(5):
        f, g = R.random(rng.randint(0, 9)), R.random(rng.randint(1, 4))
        qq, rr = fn["x_divmod"](P.anillo, P.pol_a_ref(f), P.pol_a_ref(g))
        assert (P.pol_de_ref(qq), P.pol_de_ref(rr)) == alcp.divmod_toeplitz(f, g)
        assert P.pol_de_ref(fn["x_mult_karatsuba"](P.anillo, P.pol_a_ref(f), P.pol_a_ref(g))) == \
            f.mul(g, "karatsuba")


def test_la_referencia_no_se_modifica_ni_deja_rastro(practicas):
    """Las prácticas se importan desde python/ sin escribir .pyc."""
    carpeta = Path(__file__).resolve().parents[2] / "python"
    assert sorted(p.name for p in carpeta.iterdir()) == \
        ["algoritmos_rapidos.py", "cuerpos_finitos.py", "factorizacion.py"]
    assert str(carpeta) not in sys.path
    assert practicas.cf.__file__.startswith(str(carpeta))
