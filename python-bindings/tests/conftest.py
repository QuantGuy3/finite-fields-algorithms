"""Configuración común de la suite pytest del paquete ``alcp``.

Organización de la suite (python-bindings/tests/):

==============================  ===============================================
``test_native.py``              módulo nativo ``alcp_native`` (binding pybind11)
``test_cuerpos.py``             cuerpos ``GF(q)`` y elementos
``test_polinomios.py``          anillos, polinomios y factorizaciones
``test_rapidos.py``             DFT, Toeplitz y división (alcp.fast)
``test_analizador.py``          analizador de expresiones
``test_representacion.py``      str/repr/format/LaTeX y estilos de impresión
``test_utilidades.py``          utilidades y API pública
``test_adversarial.py``         entradas hostiles y casos límite (revisión adversarial)
``test_propiedades.py``         propiedades algebraicas con hypothesis
``test_practicas_originales.py`` comparación con las prácticas de ``python/``
``test_regresiones.py``         un test por cada defecto encontrado y corregido
==============================  ===============================================

Marcadores (además de ``slow``, declarado en pyproject.toml): ``nativo``,
``adversarial``, ``propiedades``, ``referencia`` y ``regresion``.  Los de
``test_native.py`` y las pruebas con hypothesis se marcan automáticamente; por
ejemplo, ``pytest -m "not propiedades"`` excluye las pruebas de hypothesis.

Fixtures compartidas: ``estado_global`` (automática: semilla fija y estilo de
impresión restaurado tras cada prueba), ``practicas`` (las prácticas originales,
importadas sin escribir ``.pyc``) y ``potencias_limitadas`` (corta una búsqueda
que no termina en lugar de colgar la suite).
"""
from __future__ import annotations

import contextlib
import importlib
import sys
import types
from pathlib import Path

import pytest

import alcp

#: Raíz del repositorio.
RAIZ = Path(__file__).resolve().parents[2]

#: Semilla común del generador del núcleo (pruebas reproducibles).
SEMILLA = 20261018

MARCADORES = {
    "nativo": "pruebas del módulo nativo alcp_native (test_native.py)",
    "adversarial": "entradas hostiles y casos límite que intentan romper la librería",
    "propiedades": "propiedades algebraicas comprobadas con hypothesis",
    "referencia": "comparación con las prácticas originales de python/",
    "regresion": "defectos encontrados y corregidos (cada uno con su prueba)",
}


def pytest_configure(config):
    for nombre, descripcion in MARCADORES.items():
        config.addinivalue_line("markers", f"{nombre}: {descripcion}")


def pytest_collection_modifyitems(config, items):
    """Marca test_native.py como «nativo» y las pruebas de hypothesis como «propiedades»."""
    for item in items:
        if item.path.name == "test_native.py":
            item.add_marker(pytest.mark.nativo)
        if getattr(getattr(item, "obj", None), "is_hypothesis_test", False):
            item.add_marker(pytest.mark.propiedades)


@pytest.fixture(autouse=True)
def estado_global():
    """Semilla fija antes de cada prueba y estilo de impresión «unicode» después.

    Los módulos que necesitan otra semilla la fijan en su propia fixture
    automática, que se ejecuta después de esta.
    """
    alcp.set_seed(SEMILLA)
    yield
    alcp.set_print_style("unicode")


@pytest.fixture(scope="session")
def practicas():
    """Las tres prácticas originales de ``python/``, solo como referencia.

    Se importan sin escribir ``.pyc`` (no se deja rastro en python/) y sin dejar
    la carpeta en ``sys.path``.  ``factorizacion`` añade ``factorizar`` a los
    anillos de ``cuerpos_finitos``.  NO se modifican: los fallos conocidos de la
    referencia (factorización con p = 2 y con polinomios no mónicos sobre F_q,
    ``fq_ifft`` con n >= p) se evitan eligiendo los casos de comparación.
    """
    ruta = str(RAIZ / "python")
    antes = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    sys.path.insert(0, ruta)
    try:
        modulos = {nombre: importlib.import_module(nombre)
                   for nombre in ("cuerpos_finitos", "factorizacion", "algoritmos_rapidos")}
    finally:
        sys.path.remove(ruta)
        sys.dont_write_bytecode = antes
    return types.SimpleNamespace(cf=modulos["cuerpos_finitos"], fz=modulos["factorizacion"],
                                 ar=modulos["algoritmos_rapidos"])


class NativoContado:
    """Envoltorio del cuerpo nativo de un ``FiniteField`` que cuenta las potencias.

    Delega todo en el objeto nativo, pero ``pow`` lanza ``AssertionError`` al
    superar ``maximo`` llamadas: una búsqueda que no termina falla en lugar de
    colgar la suite.
    """

    def __init__(self, nat, maximo: int):
        self._nat = nat
        self._maximo = maximo
        self.llamadas = 0

    def __getattr__(self, nombre):
        return getattr(self._nat, nombre)

    def pow(self, a, k):
        self.llamadas += 1
        if self.llamadas > self._maximo:
            raise AssertionError(f"más de {self._maximo} potencias: la búsqueda no termina "
                                 "en un tiempo razonable")
        return self._nat.pow(a, k)


@pytest.fixture
def potencias_limitadas():
    """Gestor de contexto ``potencias_limitadas(K, maximo)`` (ver :class:`NativoContado`).

    Mientras dura, ``K._nat`` es el envoltorio; al salir (también si falla) se
    restaura el objeto nativo original.  No sirve para funciones que pasan
    ``K._nat`` al módulo nativo (p. ej. ``alcp.dft``): solo para las búsquedas
    en Python (raíces de la unidad, raíces cuadradas...).
    """
    @contextlib.contextmanager
    def limitar(K, maximo: int = 2000):
        original = K._nat
        K._nat = NativoContado(original, maximo)
        try:
            yield K._nat
        finally:
            K._nat = original
    return limitar
