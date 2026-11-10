// alcp/version.hpp — Versión de la librería (versionado semántico vX.Y.Z).
//
// MAJOR cambia con cambios incompatibles de la API, MINOR con funcionalidad
// nueva compatible y PATCH con correcciones.  Es la única fuente de la versión:
// el módulo pybind11 la expone como alcp_native.__version__ y el paquete Python
// alcp comprueba en sus tests que coincide con su __version__.
#pragma once

#define ALCP_VERSION_MAJOR 1
#define ALCP_VERSION_MINOR 0
#define ALCP_VERSION_PATCH 0

// Cadena "X.Y.Z" construida a partir de los tres números.
#define ALCP_VERSION_CADENA_(x) #x
#define ALCP_VERSION_CADENA(x) ALCP_VERSION_CADENA_(x)
#define ALCP_VERSION_STRING                                                              \
    ALCP_VERSION_CADENA(ALCP_VERSION_MAJOR) "." ALCP_VERSION_CADENA(ALCP_VERSION_MINOR) "." \
        ALCP_VERSION_CADENA(ALCP_VERSION_PATCH)

namespace alcp {

// Versión como valores constantes, para comparar en código C++.
inline constexpr int version_major = ALCP_VERSION_MAJOR;
inline constexpr int version_minor = ALCP_VERSION_MINOR;
inline constexpr int version_patch = ALCP_VERSION_PATCH;
inline constexpr const char* version_string = ALCP_VERSION_STRING;

}  // namespace alcp
