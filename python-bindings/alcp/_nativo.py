"""Carga del módulo nativo ``alcp.alcp_native`` con un mensaje claro si falta.

El módulo nativo (binding pybind11 del núcleo C++) se compila con
``python build.py module`` y se coloca junto a este archivo.  Se importa
SIEMPRE con el nombre ``alcp.alcp_native``: pybind11 no admite inicializar el
mismo binario con dos nombres distintos en un mismo proceso.
"""
from __future__ import annotations

import importlib
from types import ModuleType

#: Orden que compila el módulo nativo (se cita en el mensaje de error).
ORDEN_COMPILAR = "python build.py module"


def cargar(nombre: str = ".alcp_native") -> ModuleType:
    """Importa el módulo nativo ``nombre`` (relativo al paquete ``alcp``).

    Si no está compilado (o no se puede cargar), lanza ``ImportError`` con un
    mensaje en español que indica cómo compilarlo; la excepción original queda
    encadenada en ``__cause__``.
    """
    try:
        return importlib.import_module(nombre, "alcp")
    except ImportError as exc:
        raise ImportError(
            f"no se puede cargar el módulo nativo «alcp{nombre}» ({exc}). "
            f"Compílalo desde la raíz del repositorio con «{ORDEN_COMPILAR}»."
        ) from exc


#: Módulo nativo ya cargado (lo usan el resto de módulos del paquete).
nativo = cargar()
