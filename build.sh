#!/usr/bin/env bash
# Envoltorio de build.py: busca un intérprete de Python y le pasa todos los
# argumentos tal cual.
#
#   ./build.sh tests       ./build.sh oracle      ./build.sh module
#   ./build.sh all         ./build.sh ci          ./build.sh --help
#
# Orden de búsqueda del intérprete:
#   1. la variable de entorno PYTHON, si está definida;
#   2. el entorno virtual del proyecto (.venv);
#   3. python3 / python del PATH (descartando el alias de la Microsoft Store,
#      que no es un intérprete real).
set -euo pipefail

raiz="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ¿Es $1 un Python 3 utilizable?
python_valido() {
    "$1" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 9) else 1)' >/dev/null 2>&1
}

if [[ -n "${PYTHON:-}" ]]; then
    exec "$PYTHON" "$raiz/build.py" "$@"
fi
for candidato in "$raiz/.venv/Scripts/python.exe" "$raiz/.venv/bin/python"; do
    if [[ -x "$candidato" ]]; then
        exec "$candidato" "$raiz/build.py" "$@"
    fi
done
for candidato in python3 python; do
    if command -v "$candidato" >/dev/null 2>&1 && python_valido "$candidato"; then
        exec "$candidato" "$raiz/build.py" "$@"
    fi
done
echo "build.sh: no se encuentra un intérprete de Python 3 (define PYTHON=/ruta/a/python)" >&2
exit 2
