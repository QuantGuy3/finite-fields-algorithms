"""Prueba de humo del recorrido didáctico ``python-bindings/ejemplos.py``.

Ejecuta el recorrido completo en modo rápido (las medidas de tiempo usan
grados pequeños) con la salida capturada, y comprueba que termina sin
excepciones, que todas sus comprobaciones internas pasan y que, sin color, la
salida no lleva secuencias ANSI.  También prueba --lista, --seccion y la
decisión de usar color.
"""
from __future__ import annotations

import pytest

import alcp
import ejemplos


def test_recorrido_completo_sin_color(capsys):
    assert ejemplos.main(["--sin-color", "--rapido"]) == 0
    salida = capsys.readouterr().out
    assert "\x1b[" not in salida
    for titulo, _ in ejemplos.SECCIONES:
        assert titulo in salida
    assert ejemplos.Salida.comprobaciones > 0 and ejemplos.Salida.fallos == 0
    assert alcp.get_print_style() == "unicode"      # el estilo global queda como estaba


def test_seccion_con_color_forzado(capsys):
    assert ejemplos.main(["--color", "--seccion", "1"]) == 0
    salida = capsys.readouterr().out
    assert "\x1b[" in salida
    assert ejemplos.SECCIONES[0][0] in salida
    assert ejemplos.SECCIONES[1][0] not in salida


def test_lista_y_seccion_invalida(capsys):
    assert ejemplos.main(["--lista", "--sin-color"]) == 0
    salida = capsys.readouterr().out
    for n, (titulo, _) in enumerate(ejemplos.SECCIONES, start=1):
        assert f"{n:>2}  {titulo}" in salida
    for argumento in ("0", str(len(ejemplos.SECCIONES) + 1), "x"):
        with pytest.raises(SystemExit):
            ejemplos.main(["--seccion", argumento])


def test_comprobacion_fallida_da_codigo_1(monkeypatch, capsys):
    def seccion_que_falla():
        ejemplos.comprueba("esto es falso", False)

    monkeypatch.setattr(ejemplos, "SECCIONES", [("Sección de prueba", seccion_que_falla)])
    assert ejemplos.main(["--sin-color", "--seccion", "1"]) == 1
    assert "✘ esto es falso" in capsys.readouterr().out


def test_decide_color():
    decide = ejemplos.decide_color
    assert decide(False, False, {}, True)
    assert not decide(False, False, {}, False)                  # no es una terminal
    assert not decide(False, False, {"NO_COLOR": ""}, True)     # NO_COLOR definida
    assert decide(False, True, {"NO_COLOR": "1"}, False)        # --color la anula
    assert not decide(True, True, {}, True)                     # --sin-color gana a todo
