# Flujo de trabajo con git

## Ramas

| Rama | Para qué | Quién escribe en ella |
|---|---|---|
| `main` | Solo versiones publicadas. Muy pocos merges, cada uno con su etiqueta `vX.Y.Z` | Únicamente el merge de una versión desde `develop` |
| `develop` | Integración: lo que está terminado y probado | Únicamente merges `--no-ff` desde `feature/*` |
| `feature/<tema>` | Trabajo en curso, una rama por tema | Commits atómicos |

```text
main     ──●───────────────────────────────●── v1.0.0
           v0.1.0                          ╱
develop    ╲──●──────●──────●──────────────●
               ╲    ╱ ╲    ╱ ╲            ╱
feature/a       ●──●   ╲  ╱   ●──●──●────╱
feature/b               ●●
```

Ciclo de un cambio:

1. `git switch develop && git pull`
2. `git switch -c feature/<tema>`
3. Commits atómicos (ver abajo). `python build.py all` debe acabar en «TODO CORRECTO».
4. `git switch develop && git merge --no-ff feature/<tema>`, con un mensaje
   `merge(feature/<tema>): ...` en el mismo formato.
5. Para publicar una versión: actualizar `include/alcp/version.hpp`, `alcp.__version__` y
   `CHANGELOG.md` en una `feature`, integrarla en `develop`, y después
   `git switch main && git merge --no-ff develop` y `git tag -a vX.Y.Z`.

## Commits

Cada commit hace **un único cambio**: una función, una corrección, un test, un documento. Si
el mensaje necesita «y además», probablemente son dos commits.

Formato del mensaje:

```text
<tipo>(<ámbito>): <resumen de 72 caracteres como mucho>

What:
- qué cambia, concreto.
Why:
- por qué; si no consta, decirlo, no inventarlo.

Files:
- ruta/del/fichero (A|M|D|R)
```

- `<tipo>`: `feat`, `fix`, `perf`, `refactor`, `test`, `docs`, `build`, `ci`, `chore`,
  `bench`, `revert`, `style` o `merge`.
- `<ámbito>`: la parte afectada, por ejemplo `polynomial`, `fft`, `binding`, `paquete`,
  `bench`, `fase8`.
- `Files` lista lo que muestra `git diff --cached --name-status`: añadido (A), modificado
  (M), borrado (D) o renombrado (R).

## Versionado

Se usa [versionado semántico](https://semver.org/lang/es/) `vX.Y.Z`: `X` para cambios
incompatibles de la API, `Y` para funcionalidad compatible y `Z` para correcciones. La
versión vive en `include/alcp/version.hpp`. Los tests comprueban que `alcp.__version__`
coincide con `alcp_native.__version__`. Los cambios de cada versión se anotan en
[CHANGELOG.md](../CHANGELOG.md).

## Historial reescrito

El 10 de noviembre de 2026 se reescribió el historial para que siguiera esta estructura. El
contenido de cada commit no cambió (mismo árbol, autor y fecha), pero sí todos los SHAs. Los
documentos anteriores citan SHAs antiguos; la equivalencia está en
[historial-sha.md](historial-sha.md).
