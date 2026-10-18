"""Pruebas de las utilidades del paquete (alcp._util, alcp._nativo) y de su API pública."""
from __future__ import annotations

import doctest

import pytest

import alcp
from alcp import _nativo, _util


class Indice:
    def __init__(self, v):
        self.v = v

    def __index__(self):
        return self.v


def test_api_publica():
    assert alcp.__version__ and isinstance(alcp.__version__, str)
    for nombre in alcp.__all__:
        assert hasattr(alcp, nombre), nombre
    assert alcp.alcp_native is _nativo.nativo
    from alcp import alcp_native
    assert alcp_native.Fp(7).p == 7


def test_ejemplo_de_la_documentacion_del_paquete():
    alcp.set_print_style("unicode")
    resultado = doctest.testmod(alcp, verbose=False)
    assert resultado.attempted >= 6 and resultado.failed == 0


def test_modulo_nativo_ausente():
    with pytest.raises(ImportError, match="python build.py module") as info:
        _nativo.cargar(".no_existe_este_modulo")
    assert isinstance(info.value.__cause__, ImportError)
    assert "alcp.no_existe_este_modulo" in str(info.value)


def test_entero():
    assert _util.entero(5) == 5 and _util.entero(True) == 1 and _util.entero(Indice(-3)) == -3
    with pytest.raises(TypeError, match="el grado debe ser un entero, no un valor de tipo 'float'"):
        _util.entero(1.0, "el grado")


@pytest.mark.parametrize("n, esperado", [(2, True), (3, True), (4, False), (1, False),
                                         (0, False), (-7, False), (2 ** 61 - 1, True),
                                         (2 ** 61 + 1, False), (Indice(13), True)])
def test_is_prime(n, esperado):
    assert alcp.is_prime(n) is esperado


def test_is_prime_tipo():
    with pytest.raises(TypeError, match="n debe ser un entero"):
        alcp.is_prime(7.0)


@pytest.mark.parametrize("q, esperado", [(2, (2, 1)), (9, (3, 2)), (1024, (2, 10)),
                                         (Indice(125), (5, 3)), ((2 ** 61 - 1) ** 3, (2 ** 61 - 1, 3))])
def test_prime_power_factor(q, esperado):
    assert alcp.prime_power_factor(q) == esperado


@pytest.mark.parametrize("q", [0, 1, 6, 12, -9])
def test_prime_power_factor_invalido(q):
    with pytest.raises(ValueError, match=f"{q} no es potencia de un número primo"):
        alcp.prime_power_factor(q)


@pytest.mark.parametrize("n, esperado", [(-5, 2), (0, 2), (1, 2), (2, 3), (13, 17), (24, 29),
                                         (2 ** 61 - 2, 2 ** 61 - 1)])
def test_next_prime(n, esperado):
    assert alcp.next_prime(n) == esperado


def test_set_seed():
    K = alcp.GF(101)
    alcp.set_seed(10**30)
    x = [K.random() for _ in range(5)]
    alcp.set_seed(Indice(10**30))
    assert [K.random() for _ in range(5)] == x
    with pytest.raises(ValueError, match="la semilla debe ser un entero >= 0"):
        alcp.set_seed(-1)
    with pytest.raises(TypeError, match="la semilla debe ser un entero"):
        alcp.set_seed("1")


@pytest.mark.parametrize("n, esperado", [
    (1, {}), (2, {2: 1}), (12, {2: 2, 3: 1}), (97, {97: 1}), (2 * 3 * 5 * 7 * 11, None),
    (2 ** 10 * 3 ** 5 * 1009, None), (998244352, {2: 23, 7: 1, 17: 1}),
    ((2 ** 61 - 1) * 4, {2: 2, 2 ** 61 - 1: 1}),
])
def test_factorizar(n, esperado):
    f = _util.factorizar(n)
    producto = 1
    for p, e in f.items():
        assert alcp.is_prime(p)
        producto *= p ** e
    assert producto == n
    if esperado is not None:
        assert f == esperado


def test_factorizar_cofactor_compuesto_grande(monkeypatch):
    monkeypatch.setattr(_util, "LIMITE_DIVISION", 100)
    with pytest.raises(ValueError, match="el cofactor 10403 es compuesto"):
        _util.factorizar(2 * 101 * 103)
