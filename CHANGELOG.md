# Registro de cambios

Formato basado en [Keep a Changelog](https://keepachangelog.com/es-ES/1.1.0/). El proyecto
sigue [versionado semántico](https://semver.org/lang/es/) `vX.Y.Z`:

- `X`: cambios incompatibles de la API.
- `Y`: funcionalidad nueva compatible.
- `Z`: correcciones.

La versión está definida en un único sitio, `include/alcp/version.hpp`. El módulo nativo la
expone como `alcp_native.__version__`, y el paquete `alcp` declara la misma.

## [1.0.0] - 2026-11-10

Primera versión publicada en `main`. Reúne las fases 1 a 8 desarrolladas en `develop`, cuyos
informes están en `docs/optimizacion/`.

### Añadido
- **Núcleo C++23 de solo cabeceras** (`include/alcp/`): cuerpos `Fp` (GMP), `Fp64`,
  `FpN<2>`, `Fq`, `FqTab` y `FqSwar`, y polinomios con producto escolar, Karatsuba,
  Kronecker y FFT/NTT, división de Newton, half-GCD, irreducibilidad y factorización
  (Cantor–Zassenhaus con DDF de Kaltofen–Shoup).
- **Módulo pybind11 `alcp_native` y paquete Python `alcp`**, con cobertura de tests del
  100 %.
- **Herramientas de trabajo:** batería de tests C++, oráculo contra las prácticas Python
  originales, benchmarks (`bench/`), comparativa con python-flint, galois y SymPy, y CI.
- **Versionado semántico:** `version.hpp`, `alcp_native.__version__` y este registro.

### Rendimiento
- **Frente al código original**, en la suite de 128 escenarios: fases 3 a 7. Por ejemplo,
  1,72x en la fase 7 y 1,076x en la fase 8 sobre la anterior.
- **Frente a FLINT:** alcp es más rápido en todos los escenarios de la comparativa salvo
  `factor` de grado 40 en F_101 (0,86x). Ver `docs/optimizacion/fase8.md`.

## [0.1.0] - 2026-06-30

Subida original del repositorio: las prácticas Python de cuerpos finitos, un primer port
C++, sus tests y benchmarks, el binding pybind11 y el oráculo.

[1.0.0]: https://github.com/QuantGuy3/finite-fields-algorithms/compare/v0.1.0...v1.0.0
[0.1.0]: https://github.com/QuantGuy3/finite-fields-algorithms/releases/tag/v0.1.0
