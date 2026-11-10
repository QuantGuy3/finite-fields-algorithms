"""mide_binding.py — Fase 6, propuesta 8: conversiones del binding, base frente a final.

Carga un alcp_native.pyd concreto (el de la base o el de la fase 6) y mide,
en F_998244353 y en un único proceso por módulo:
  * K.poly(lista) de grado 100 y 1000 (conversión de entrada),
  * f.coeffs() de grado 100 y 1000 (conversión de salida),
  * K.mul(a, b) y K.add(a, b) con elementos (enteros pequeños).
Mejor de 7 tandas de >= 50 ms, 16 entradas rotando, gc desactivado.  También
imprime un hash de los resultados para comprobar que ambos módulos coinciden.

Uso (desde la raíz del repositorio, con el Python del .venv):
  .venv/Scripts/python.exe bench/experimentos/fase6/python/mide_binding.py <directorio del .pyd>
"""
from __future__ import annotations

import gc
import hashlib
import importlib.util
import sys
import time
from pathlib import Path

P = 998244353


def cargar(directorio: str):
    """Importa alcp_native desde el directorio dado (sin tocar sys.path global)."""
    pyd = next(Path(directorio).glob("alcp_native*.pyd"))
    spec = importlib.util.spec_from_file_location("alcp_native", pyd)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def mejor_ns(fn, entradas) -> float:
    """Mejor tiempo medio por llamada (ns) en 7 tandas de >= 50 ms."""
    gc.disable()
    try:
        mejor = float("inf")
        n = len(entradas)
        for _ in range(7):
            veces = 0
            t0 = time.perf_counter_ns()
            while True:
                for e in entradas:
                    fn(e)
                veces += n
                dt = time.perf_counter_ns() - t0
                if dt >= 50_000_000:
                    break
            mejor = min(mejor, dt / veces)
        return mejor
    finally:
        gc.enable()


def main() -> None:
    nat = cargar(sys.argv[1])
    K = nat.Fp(P)
    # Entradas deterministas (LCG sencillo): mismas en los dos módulos.
    x = 12345

    def sig() -> int:
        nonlocal x
        x = (x * 6364136223846793005 + 1442695040888963407) % (1 << 64)
        return x % P

    listas = {d: [[sig() for _ in range(d + 1)] for _ in range(16)] for d in (100, 1000)}
    h = hashlib.sha256()
    for d, ls in listas.items():
        polys = [K.poly(c) for c in ls]
        for f in polys:
            h.update(repr(f.coeffs()).encode())
        print(f"K.poly(lista) grado {d:5d}  {mejor_ns(lambda c: K.poly(c), ls):10.0f} ns")
        print(f"coeffs()      grado {d:5d}  {mejor_ns(lambda f: f.coeffs(), polys):10.0f} ns")
    pares = [(sig(), sig()) for _ in range(16)]
    for a, b in pares:
        h.update(repr((K.mul(a, b), K.add(a, b))).encode())
    print(f"K.mul(a, b)                 {mejor_ns(lambda e: K.mul(e[0], e[1]), pares):10.0f} ns")
    print(f"K.add(a, b)                 {mejor_ns(lambda e: K.add(e[0], e[1]), pares):10.0f} ns")
    print("hash", h.hexdigest()[:16])


if __name__ == "__main__":
    main()
