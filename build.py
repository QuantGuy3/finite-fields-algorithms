#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build.py -- driver de compilación y verificación multiplataforma.

Funciona en Windows (g++ de MSYS2/MinGW + CPython oficial de python.org, o el
Python de MSYS2 para todo lo que no sea el módulo) y en Linux (g++ + libgmp-dev).
Requiere Python >= 3.9 y solo usa la biblioteca estándar; pybind11, pytest y
gcovr se buscan en el intérprete con el que se ejecuta este script (normalmente
el del .venv; build.sh lo elige).  Funciona desde cualquier directorio.

Uso:
    python build.py [-v] [--no-color] <subcomando> [opciones]
    (en pytest y bench, lo que sigue al subcomando es para ellos:
     «build.py -v pytest -k foo» es -v de build.py y -k foo de pytest)

Subcomandos:
    tests     Compila cada tests/*.cpp en build/tests/ y lo ejecuta.  crosscheck.cpp
              se compila como «xcheck» (Windows Defender bloquea un ejecutable
              llamado crosscheck.exe) y solo se ejecuta si existe
              build/oracle_vectors.txt (lo genera el subcomando oracle).
    oracle    python-bindings/oracle.py -> build/oracle_vectors.txt y xcheck sobre él.
    module    Compila python-bindings/bindings.cpp ->
              python-bindings/alcp/alcp_native<EXT_SUFFIX> y comprueba que importa.
              Opciones: --debug (-O0 -g), --coverage (-O0 -g --coverage).
    pytest    Ejecuta pytest con la configuración de pyproject.toml; los
              argumentos que siguen al subcomando se pasan tal cual.
    cov-cpp   Cobertura C++ con gcov + gcovr -> build/coverage-cpp/ (ver abajo).
    bench     Delega en bench/run_bench.py (argumentos tal cual); si no existe,
              compila y ejecuta bench/bench.cpp (-O3 -march=native).
    all       tests + oracle + module + pytest (si hay tests).
    ci        Igual que all, sin colores y con grupos de GitHub Actions.
    clean     Borra build/ y los módulos compilados.

Opciones comunes:
    -v / --verbose   muestra cada comando antes de ejecutarlo.
    --no-color       sin colores ANSI (también con la variable NO_COLOR).
    --flags "..."    flags extra del compilador (se añaden tras los de serie).
    -j N             compilaciones en paralelo (por defecto, nº de CPU).
    -B / --force     recompila aunque los binarios estén al día.

Variables de entorno:
    CXX              compilador C++ (por defecto g++).
    GCOV             ejecutable gcov para cov-cpp (por defecto, el que acompaña a CXX).
    ALCP_CXXFLAGS    flags extra del compilador (antes de los de --flags).

Cobertura C++ (cov-cpp):
    Todo lo instrumentado vive en build/cov/: objetos y .gcno/.gcda en
    build/cov/obj/, tests en build/cov/tests/ y el módulo en build/cov/module/.
    El módulo NORMAL de python-bindings/alcp/ no se pisa de forma permanente:
    la suite pytest importa el paquete «alcp», que carga su módulo nativo
    desde su propia carpeta; para ella el módulo instrumentado se coloca
    TEMPORALMENTE en python-bindings/alcp/ y, al terminar (también si pytest
    falla o se interrumpe con Ctrl+C), se restaura el original byte a byte
    desde la copia guardada en build/cov/.
    Si el proceso muriera de forma abrupta en medio de pytest, bastaría con
    «python build.py module» para recompilar el módulo normal.
    El informe (HTML detallado + resumen de texto) filtra include/alcp/ y
    python-bindings/bindings.cpp.

Detalles de Windows:
    * Los ejecutables (tests, bench) se enlazan con -static: así no dependen
      del PATH (en Git Bash, la libstdc++-6.dll de /mingw64/bin, incompatible
      con la de ucrt64, se cargaría antes).  El módulo también enlaza GMP y el
      runtime de C++ estáticamente: solo depende de pythonXY.dll.
    * Si el antivirus retiene un binario recién creado se reintenta; si Smart
      App Control bloquea un binario concreto (WinError 4551), se vuelve a
      enlazar (cambia su hash) y se reintenta.  Un test que se ejecuta y
      falla NUNCA se reintenta: su fallo se propaga siempre.

Compilación incremental: cada binario guarda en build/sellos/ el comando que
lo generó; se recompila si cambia el comando o alguna fuente/cabecera es más
reciente (-B fuerza la recompilación).

Códigos de salida: 0 = todo correcto, 1 = algún paso ha fallado, 2 = error de
uso o de entorno (falta el compilador, pybind11, etc.).
"""
from __future__ import annotations

import argparse
import concurrent.futures
import importlib.util
import os
import shlex
import shutil
import subprocess
import sys
import sysconfig
import time
from dataclasses import dataclass, field, replace
from pathlib import Path

# ---------------------------------------------------------------------------
# Rutas del proyecto (absolutas: el script funciona desde cualquier directorio)
# ---------------------------------------------------------------------------
RAIZ = Path(__file__).resolve().parent
INCLUDE = RAIZ / "include"
TESTS = RAIZ / "tests"
BINDINGS = RAIZ / "python-bindings"
PAQUETE = BINDINGS / "alcp"
PRACTICAS = RAIZ / "python"
BUILD = RAIZ / "build"
VECTORES = BUILD / "oracle_vectors.txt"
COV = BUILD / "cov"
INFORME_COV = BUILD / "coverage-cpp"

ES_WINDOWS = os.name == "nt"
EXE = ".exe" if ES_WINDOWS else ""
# Python de MSYS2/MinGW (sysconfig.get_platform() == "mingw_x86_64_ucrt_gnu"...)
PYTHON_MINGW = sysconfig.get_platform().startswith("mingw")

FLAGS_BASE = ["-std=c++23", "-Wall", "-Wextra"]
LIBS_GMP = ["-lgmpxx", "-lgmp"]
# Ejecutables (tests, bench): en Windows se enlazan estáticamente.  Así no
# dependen del PATH: en Git Bash, /mingw64/bin (el MinGW de Git, con msvcrt)
# va antes que C:/msys64/ucrt64/bin y su libstdc++-6.dll, incompatible, se
# cargaría en lugar de la buena (comprobado: std::ifstream de un fichero
# inexistente NO marcaba failbit con la DLL equivocada).
LIBS_EXE = (["-static"] if ES_WINDOWS else []) + LIBS_GMP

# Reintentos cuando el antivirus bloquea un ejecutable recién creado.
REINTENTOS_AV = 6
ESPERA_AV = 0.5          # segundos; se multiplica por el nº de intento


class ErrorEntorno(Exception):
    """Falta una herramienta o un requisito: código de salida 2."""


# ---------------------------------------------------------------------------
# Salida legible (secciones, OK/FALLO, tiempos)
# ---------------------------------------------------------------------------
def _activar_vt_windows() -> bool:
    """Activa las secuencias ANSI en la consola de Windows (si es posible)."""
    try:
        import ctypes
        k32 = ctypes.windll.kernel32
        h = k32.GetStdHandle(-11)                     # STD_OUTPUT_HANDLE
        modo = ctypes.c_uint32()
        if not k32.GetConsoleMode(h, ctypes.byref(modo)):
            return False
        return bool(k32.SetConsoleMode(h, modo.value | 0x0004))  # VT processing
    except Exception:
        return False


class Consola:
    COLORES = {"OK": "32", "FALLO": "31", "SALTADO": "33", "AL DÍA": "36",
               "AVISOS": "33", "titulo": "1;36", "tenue": "2"}

    def __init__(self) -> None:
        self.verbose = False
        self.color = False
        self.gha = False            # agrupar secciones al estilo GitHub Actions

    def configurar(self, verbose: bool, sin_color: bool, modo_ci: bool) -> None:
        self.verbose = verbose
        self.gha = modo_ci and os.environ.get("GITHUB_ACTIONS") == "true"
        quiere = (not sin_color and not modo_ci and "NO_COLOR" not in os.environ
                  and sys.stdout.isatty())
        self.color = quiere and (not ES_WINDOWS or _activar_vt_windows())

    def _c(self, clave: str, texto: str) -> str:
        if not self.color or clave not in self.COLORES:
            return texto
        return f"\033[{self.COLORES[clave]}m{texto}\033[0m"

    def seccion(self, titulo: str, plegable: bool = True) -> None:
        if self.gha and plegable:
            print(f"::group::{titulo}", flush=True)
        else:
            print("\n" + self._c("titulo", f"== {titulo} " + "=" * max(4, 70 - len(titulo))),
                  flush=True)

    def fin_seccion(self) -> None:
        if self.gha:
            print("::endgroup::", flush=True)

    def linea(self, estado: str, nombre: str, segundos: float | None = None,
              detalle: str = "") -> None:
        t = f"{segundos:7.2f} s" if segundos is not None else " " * 9
        est = self._c(estado, f"{estado:<7}")
        print(f"  {est} {nombre:<22} {t}  {detalle}".rstrip(), flush=True)

    def nota(self, texto: str) -> None:
        print(f"  {self._c('tenue', texto)}", flush=True)

    def bloque(self, texto: str, max_lineas: int | None = 60) -> None:
        """Imprime salida capturada, sangrada (las últimas max_lineas)."""
        lineas = texto.rstrip().splitlines()
        if not lineas:
            return
        if max_lineas is not None and len(lineas) > max_lineas and not self.verbose:
            print(f"    ... ({len(lineas) - max_lineas} líneas omitidas; usa -v para verlas)")
            lineas = lineas[-max_lineas:]
        for ln in lineas:
            print("    " + ln)
        sys.stdout.flush()

    def comando(self, cmd: list[str], cwd: Path | None = None) -> None:
        if self.verbose:
            texto = subprocess.list2cmdline(cmd) if ES_WINDOWS else shlex.join(cmd)
            donde = f"   (en {cwd})" if cwd and cwd != RAIZ else ""
            print(self._c("tenue", f"  $ {texto}{donde}"), flush=True)

    def error(self, texto: str) -> None:
        if self.gha:
            print(f"::error::{texto}", flush=True)
        print(self._c("FALLO", f"  ERROR: {texto}"), flush=True)


UI = Consola()


@dataclass
class Paso:
    """Resultado de un paso de alto nivel (aparece en el resumen final)."""
    nombre: str
    estado: str               # OK | FALLO | SALTADO
    segundos: float = 0.0
    detalle: str = ""


@dataclass
class Resultado:
    rc: int
    salida: str = ""
    segundos: float = 0.0
    al_dia: bool = False      # compilación omitida porque el binario está al día
    avisos: bool = False      # el compilador emitió avisos

    @property
    def ok(self) -> bool:
        return self.rc == 0


@dataclass
class Opciones:
    flags: list[str] = field(default_factory=list)
    jobs: int = os.cpu_count() or 2
    force: bool = False


# ---------------------------------------------------------------------------
# Ejecución de procesos (siempre con listas de argumentos, nunca shell=True)
# ---------------------------------------------------------------------------
def correr(cmd: list[str], *, env: dict | None = None, timeout: float | None = None,
           capturar: bool = True, reintentar_av: bool = False,
           cwd: Path = RAIZ) -> Resultado:
    """Ejecuta `cmd` y devuelve su código de salida, salida combinada y tiempo.

    Con `reintentar_av`, si el sistema impide lanzar el ejecutable (antivirus
    escaneando un binario recién enlazado: PermissionError en Windows o código
    126 en entornos POSIX) se reintenta unas cuantas veces con espera
    creciente.  Cualquier otro código distinto de 0 se devuelve tal cual: los
    fallos de los tests se propagan siempre.
    """
    UI.comando(cmd, cwd)
    intentos = REINTENTOS_AV if reintentar_av else 1
    t0 = time.perf_counter()
    for intento in range(1, intentos + 1):
        try:
            r = subprocess.run(
                cmd, cwd=cwd, env=env, timeout=timeout,
                stdout=subprocess.PIPE if capturar else None,
                stderr=subprocess.STDOUT if capturar else None,
                text=True, encoding="utf-8", errors="replace")
        except FileNotFoundError as ex:
            return Resultado(127, f"no se encuentra el ejecutable {cmd[0]!r}: {ex}",
                             time.perf_counter() - t0)
        except subprocess.TimeoutExpired as ex:
            parcial = ex.output if isinstance(ex.output, str) else ""
            return Resultado(124, (parcial or "") + f"\n[tiempo agotado: {timeout} s]",
                             time.perf_counter() - t0)
        except OSError as ex:
            # PermissionError (WinError 5/32: binario retenido por el antivirus),
            # WinError 4551 (Smart App Control / WDAC aún evaluando un binario
            # nuevo), WinError 225 (Defender)...  Se reintenta y, si persiste,
            # se informa como fallo (nunca se da por bueno).
            if intento < intentos:
                time.sleep(ESPERA_AV * intento)
                continue
            return Resultado(126, f"no se pudo ejecutar {cmd[0]} tras {intentos} intento(s) "
                                  f"(¿bloqueado por el antivirus o por Smart App Control?): {ex}",
                             time.perf_counter() - t0)
        if r.returncode == 126 and intento < intentos:
            time.sleep(ESPERA_AV * intento)
            continue
        return Resultado(r.returncode, r.stdout or "", time.perf_counter() - t0)
    raise AssertionError("inalcanzable")


def env_ejecucion(extra_path: list[Path] | None = None, **variables: str) -> dict:
    """Entorno para lanzar binarios y scripts.

    En Windows se antepone al PATH la carpeta del compilador, para que los
    ejecutables encuentren las DLL de MinGW (libgmp, libstdc++...) aunque el
    script se lance desde PowerShell/cmd sin MSYS2 en el PATH.
    """
    env = dict(os.environ)
    rutas = [str(p) for p in (extra_path or [])]
    if ES_WINDOWS:
        cxx = shutil.which(compilador())
        if cxx:
            rutas.append(str(Path(cxx).parent))
    if rutas:
        env["PATH"] = os.pathsep.join(rutas + [env.get("PATH", "")])
    env.update(variables)
    return env


def borrar(p: Path) -> None:
    """Borra un archivo o un árbol, reintentando unos instantes: en Windows el
    antivirus (o un .pyd recién descargado de memoria) puede retener un archivo
    lo justo para que falle el borrado del directorio que lo contiene."""
    def quitar_solo_lectura(func, ruta, _exc):
        os.chmod(ruta, 0o700)
        func(ruta)
    for intento in range(1, REINTENTOS_AV + 1):
        try:
            if p.is_dir():
                # onexc desde Python 3.12 (onerror queda obsoleto).
                clave = "onexc" if sys.version_info >= (3, 12) else "onerror"
                shutil.rmtree(p, **{clave: quitar_solo_lectura})
            elif p.exists():
                p.unlink()
            return
        except OSError:
            if intento == REINTENTOS_AV:
                raise
            time.sleep(ESPERA_AV * intento)


# ---------------------------------------------------------------------------
# Compilación
# ---------------------------------------------------------------------------
def compilador() -> str:
    return os.environ.get("CXX", "g++")


def comprobar_compilador() -> str:
    cxx = compilador()
    if shutil.which(cxx) is None:
        pista = ("instala MSYS2 y el paquete mingw-w64-ucrt-x86_64-gcc, y añade "
                 "C:/msys64/ucrt64/bin al PATH" if ES_WINDOWS else
                 "instala g++ (p. ej. «sudo apt install g++ libgmp-dev»)")
        raise ErrorEntorno(f"no se encuentra el compilador {cxx!r}: {pista}")
    return cxx


def flags_extra(cli: list[str] | None) -> list[str]:
    """Flags extra: primero ALCP_CXXFLAGS y después --flags (gana el último)."""
    extra = shlex.split(os.environ.get("ALCP_CXXFLAGS", ""))
    for f in cli or []:
        extra += shlex.split(f)
    return extra


def cabeceras() -> list[Path]:
    """Dependencias comunes: todas las cabeceras del núcleo y de los tests."""
    deps = [p for p in INCLUDE.rglob("*") if p.suffix in (".hpp", ".h")]
    return deps + sorted(TESTS.glob("*.hpp"))


def ruta_sello(salida: Path) -> Path:
    """Sello con el comando que generó `salida` (en build/sellos/, para no
    dejar ficheros auxiliares junto al módulo dentro del paquete)."""
    try:
        rel = salida.resolve().relative_to(RAIZ).as_posix()
    except ValueError:
        rel = salida.resolve().as_posix().replace(":", "")
    return BUILD / "sellos" / (rel.replace("/", "__") + ".cmd")


def compilar(cmd: list[str], salida: Path, deps: list[Path], force: bool) -> Resultado:
    """Compila con `cmd` salvo que `salida` esté al día.

    Al día = existe, se generó con el MISMO comando (sello en build/sellos/) y
    es más reciente que todas sus dependencias.  Si el enlazador no puede
    sobrescribir el binario porque el antivirus lo tiene abierto, se reintenta.
    """
    sello = ruta_sello(salida)
    firma = "\n".join(cmd)
    if not force and salida.exists() and sello.exists():
        try:
            mas_nuevo = max(p.stat().st_mtime for p in deps)
            if (sello.read_text(encoding="utf-8") == firma
                    and salida.stat().st_mtime >= mas_nuevo):
                return Resultado(0, al_dia=True)
        except (OSError, ValueError):
            pass
    salida.parent.mkdir(parents=True, exist_ok=True)
    sello.parent.mkdir(parents=True, exist_ok=True)
    sello.unlink(missing_ok=True)
    r = correr(cmd)
    for intento in range(1, 4):
        if r.ok or "Permission denied" not in r.salida:
            break
        time.sleep(ESPERA_AV * intento)        # binario bloqueado por el antivirus
        r = correr(cmd)
    if r.ok:
        sello.write_text(firma, encoding="utf-8")
        r.avisos = bool(r.salida.strip())
    return r


def nombre_test(src: Path) -> str:
    # crosscheck -> xcheck: Defender bloquea por su nombre un «crosscheck.exe».
    return "xcheck" if src.stem == "crosscheck" else src.stem


def fuentes_tests() -> list[Path]:
    return sorted(TESTS.glob("*.cpp"))


def cmd_test(src: Path, exe: Path, extra: list[str]) -> list[str]:
    return [compilador(), *FLAGS_BASE, "-O2", "-I", str(INCLUDE), *extra,
            str(src), "-o", str(exe), *LIBS_EXE]


def informar_compilacion(nombre: str, r: Resultado) -> None:
    if r.al_dia:
        UI.linea("AL DÍA", f"compilar {nombre}")
    elif r.ok and r.avisos:
        UI.linea("AVISOS", f"compilar {nombre}", r.segundos, "compila, pero con avisos:")
        UI.bloque(r.salida)
    elif r.ok:
        UI.linea("OK", f"compilar {nombre}", r.segundos)
    else:
        UI.linea("FALLO", f"compilar {nombre}", r.segundos)
        UI.bloque(r.salida)
        UI.error(f"{nombre} no compila")


def compilar_varios(trabajos: list[tuple[str, callable]], jobs: int) -> dict[str, Resultado]:
    """Ejecuta en paralelo funciones de compilación (subprocess libera el GIL)
    e informa de cada una en el orden de la lista."""
    resultados: dict[str, Resultado] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        futuros = [(nombre, pool.submit(f)) for nombre, f in trabajos]
        for nombre, fut in futuros:
            resultados[nombre] = r = fut.result()
            informar_compilacion(nombre, r)
    return resultados


def ultima_linea(texto: str) -> str:
    lineas = [ln.strip() for ln in texto.strip().splitlines() if ln.strip()]
    return lineas[-1] if lineas else ""


def ejecutar_test(nombre: str, cmd: list[str], timeout: float,
                  reconstruir=None) -> bool:
    """Ejecuta un test e informa.  Si el sistema bloquea el ejecutable de forma
    persistente (rc 126 de correr(): p. ej. Smart App Control, cuyo veredicto
    depende del hash del binario), se vuelve a generar con `reconstruir()` (el
    nuevo enlace cambia la marca de tiempo PE y, con ello, el hash) y se
    reintenta.  Un test que se ejecuta y falla NUNCA se reintenta."""
    r = correr(cmd, env=env_ejecucion(), timeout=timeout, reintentar_av=True)
    for _ in range(2):
        if r.rc != 126 or reconstruir is None:
            break
        UI.nota(f"{nombre}: el sistema bloquea el ejecutable; se regenera y se reintenta")
        if not reconstruir().ok:
            break
        r = correr(cmd, env=env_ejecucion(), timeout=timeout, reintentar_av=True)
    if r.ok:
        UI.linea("OK", f"ejecutar {nombre}", r.segundos, ultima_linea(r.salida))
        if UI.verbose:
            UI.bloque(r.salida, None)
        return True
    UI.linea("FALLO", f"ejecutar {nombre}", r.segundos, f"código de salida {r.rc}")
    UI.bloque(r.salida)
    UI.error(f"el test {nombre} ha fallado (código {r.rc})")
    return False


def hay_vectores() -> bool:
    """¿Existen vectores del oráculo NO vacíos?  (Se comprueba aquí porque un
    xcheck que no puede abrir el archivo no debe darse nunca por bueno.)"""
    try:
        return VECTORES.stat().st_size > 0
    except OSError:
        return False


# ---------------------------------------------------------------------------
# Subcomando: tests
# ---------------------------------------------------------------------------
def paso_tests(op: Opciones, incluir_xcheck: bool = True, timeout: float = 600) -> Paso:
    comprobar_compilador()
    t0 = time.perf_counter()
    UI.seccion("Tests C++ (tests/*.cpp -> build/tests/)")
    destino = BUILD / "tests"
    extra = flags_extra(op.flags)
    deps_comunes = cabeceras()
    trabajos, exes, rehacer = [], {}, {}
    for src in fuentes_tests():
        nombre = nombre_test(src)
        exe = destino / (nombre + EXE)
        exes[nombre] = exe
        cmd = cmd_test(src, exe, extra)
        trabajos.append((nombre, lambda c=cmd, e=exe, s=src:
                         compilar(c, e, [s, *deps_comunes], op.force)))
        rehacer[nombre] = lambda c=cmd, e=exe: compilar(c, e, [], True)
    if not trabajos:
        UI.fin_seccion()
        return Paso("tests", "FALLO", 0, "no hay tests/*.cpp")
    compilados = compilar_varios(trabajos, op.jobs)
    fallos = [n for n, r in compilados.items() if not r.ok]

    for nombre, exe in exes.items():
        if nombre in fallos:
            continue
        if nombre == "xcheck":
            if not incluir_xcheck:
                UI.linea("SALTADO", "ejecutar xcheck", detalle="(lo ejecuta el paso oracle)")
                continue
            if not hay_vectores():
                UI.linea("SALTADO", "ejecutar xcheck",
                         detalle=f"falta {VECTORES.relative_to(RAIZ).as_posix()} "
                                 "(genéralo con «build.py oracle»)")
                continue
            ok = ejecutar_test(nombre, [str(exe), str(VECTORES)], timeout, rehacer[nombre])
        else:
            ok = ejecutar_test(nombre, [str(exe)], timeout, rehacer[nombre])
        if not ok:
            fallos.append(nombre)
    UI.fin_seccion()
    seg = time.perf_counter() - t0
    if fallos:
        return Paso("tests", "FALLO", seg, "fallan: " + ", ".join(fallos))
    return Paso("tests", "OK", seg, f"{len(exes)} tests")


# ---------------------------------------------------------------------------
# Subcomando: oracle
# ---------------------------------------------------------------------------
def paso_oracle(op: Opciones, timeout: float = 600) -> Paso:
    comprobar_compilador()
    t0 = time.perf_counter()
    UI.seccion("Oráculo: prácticas Python originales -> build/oracle_vectors.txt -> xcheck")
    BUILD.mkdir(parents=True, exist_ok=True)
    r = correr([sys.executable, str(BINDINGS / "oracle.py"), str(VECTORES)],
               env=env_ejecucion(PYTHONDONTWRITEBYTECODE="1", PYTHONIOENCODING="utf-8"),
               timeout=timeout)
    if not r.ok or not hay_vectores():
        UI.linea("FALLO", "generar vectores", r.segundos, f"código {r.rc}")
        UI.bloque(r.salida)
        UI.fin_seccion()
        return Paso("oracle", "FALLO", time.perf_counter() - t0, "oracle.py ha fallado")
    UI.linea("OK", "generar vectores", r.segundos)
    UI.bloque(r.salida, None)          # incluye el resumen de oracle.py (stderr)

    src = TESTS / "crosscheck.cpp"
    exe = BUILD / "tests" / ("xcheck" + EXE)
    cmd = cmd_test(src, exe, flags_extra(op.flags))
    rc = compilar(cmd, exe, [src, *cabeceras()], op.force)
    informar_compilacion("xcheck", rc)
    ok = rc.ok and ejecutar_test("xcheck", [str(exe), str(VECTORES)], timeout,
                                 lambda: compilar(cmd, exe, [], True))
    UI.fin_seccion()
    return Paso("oracle", "OK" if ok else "FALLO", time.perf_counter() - t0,
                "" if ok else "xcheck no pasa")


# ---------------------------------------------------------------------------
# Subcomando: module
# ---------------------------------------------------------------------------
def sufijo_modulo() -> str:
    suf = sysconfig.get_config_var("EXT_SUFFIX")
    if not suf:
        suf = ".pyd" if ES_WINDOWS else ".so"
    return suf


def include_pybind11() -> str:
    try:
        import pybind11
    except ImportError:
        raise ErrorEntorno(f"pybind11 no está instalado en {sys.executable} "
                           "(instálalo con: python -m pip install pybind11)") from None
    return pybind11.get_include()


def cmds_modulo(salida: Path, perfil: str, extra: list[str],
                obj_dir: Path | None) -> list[list[str]]:
    """Comandos para compilar el módulo nativo (uno, o compilar+enlazar).

    perfil: "release" (-O3), "debug" (-O0 -g) o "coverage" (-O0 -g --coverage).
    Con `obj_dir` se compila primero a objeto (así .gcno/.gcda quedan allí y no
    junto al módulo) y después se enlaza.
    """
    cxx = compilador()
    opt = {"release": ["-O3"], "debug": ["-O0", "-g"],
           "coverage": ["-O0", "-g", "--coverage"]}[perfil]
    py_inc = sysconfig.get_path("include")
    compilar_f = [*FLAGS_BASE, "-I", str(INCLUDE),
                  # -isystem: los avisos de cabeceras ajenas no cuentan
                  "-isystem", include_pybind11(), "-isystem", py_inc]
    if ES_WINDOWS and not PYTHON_MINGW:
        # g++ (MinGW) + CPython oficial (MSVC): se enlaza contra pythonXY.lib y
        # se enlazan estáticamente GMP, libstdc++, libgcc y winpthread, de modo
        # que el .pyd solo depende de DLL del sistema y de pythonXY.dll.
        xy = f"{sys.version_info.major}{sys.version_info.minor}"
        if sysconfig.get_config_var("Py_GIL_DISABLED"):
            xy += "t"
        libs = Path(sys.base_prefix) / "libs"
        if not (libs / f"python{xy}.lib").exists():
            raise ErrorEntorno(f"no se encuentra {libs / f'python{xy}.lib'}")
        enlazar_f = ["-shared", "-L", str(libs), f"-lpython{xy}", "-static", *LIBS_GMP]
    elif ES_WINDOWS:
        # Python de MSYS2 (MinGW): se enlaza contra libpythonX.Y.dll.a.
        # (Camino no probado en este proyecto; lo habitual es el CPython oficial.)
        enlazar_f = ["-shared", "-L", sysconfig.get_config_var("LIBDIR"),
                     f"-lpython{sysconfig.get_config_var('VERSION')}", *LIBS_GMP]
    else:
        compilar_f.append("-fPIC")
        enlazar_f = ["-shared", *LIBS_GMP]
        if sys.platform == "darwin":                  # no probado
            enlazar_f[1:1] = ["-undefined", "dynamic_lookup"]
    src = str(BINDINGS / "bindings.cpp")
    if obj_dir is None:
        return [[cxx, *compilar_f, *opt, *extra, src, "-o", str(salida), *enlazar_f]]
    obj = obj_dir / "bindings.o"
    cov = ["--coverage"] if perfil == "coverage" else []
    return [[cxx, *compilar_f, *extra, *opt, "-c", src, "-o", str(obj)],
            [cxx, *cov, str(obj), "-o", str(salida), *enlazar_f]]


def construir_modulo(salida: Path, perfil: str, extra: list[str], force: bool,
                     obj_dir: Path | None = None) -> Resultado:
    cmds = cmds_modulo(salida, perfil, extra, obj_dir)
    deps = [BINDINGS / "bindings.cpp", *cabeceras()]
    if len(cmds) == 1:
        return compilar(cmds[0], salida, deps, force)
    obj_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    salidas = []
    for cmd in cmds:                          # compilar a objeto y enlazar
        r = compilar(cmd, Path(cmd[cmd.index("-o") + 1]), deps, True)
        salidas.append(r.salida)
        if not r.ok:
            r.segundos = time.perf_counter() - t0
            return r
    texto = "\n".join(s for s in salidas if s.strip())
    return Resultado(0, texto, time.perf_counter() - t0, avisos=bool(texto.strip()))


def comprobar_import(carpeta: Path) -> Resultado:
    """Importa el módulo recién compilado en un proceso aparte."""
    codigo = ("import sys; sys.path.insert(0, sys.argv[1]); import alcp_native as m; "
              "F = m.Fp(7); print('importado:', m.__file__); "
              "print('Fp(7): (3x^2+2x+1)(x+1) =', F.poly([1, 2, 3]) * F.poly([1, 1]))")
    return correr([sys.executable, "-c", codigo, str(carpeta)], env=env_ejecucion(),
                  timeout=120, reintentar_av=True)


def importar_con_reintento(carpeta: Path, regenerar) -> Resultado:
    """comprobar_import() que, si Smart App Control bloquea ESTE binario
    concreto (su veredicto depende del hash), lo regenera con `regenerar()`
    (un nuevo enlace cambia la marca de tiempo PE y, con ello, el hash) y
    reintenta.  Cualquier otro error se devuelve tal cual."""
    ri = comprobar_import(carpeta)
    for _ in range(2):
        if ri.ok or "Application Control" not in ri.salida:
            break
        UI.nota("el sistema bloquea el módulo recién compilado; se regenera y se reintenta")
        if not regenerar().ok:
            break
        ri = comprobar_import(carpeta)
    return ri


def paso_module(op: Opciones, perfil: str = "release") -> Paso:
    comprobar_compilador()
    include_pybind11()          # requisito: falla aquí, antes de abrir la sección
    t0 = time.perf_counter()
    salida = PAQUETE / ("alcp_native" + sufijo_modulo())
    UI.seccion(f"Módulo nativo ({perfil}) -> {salida.relative_to(RAIZ).as_posix()}")
    PAQUETE.mkdir(parents=True, exist_ok=True)
    obj_dir = BUILD / "module-obj" if perfil == "coverage" else None
    r = construir_modulo(salida, perfil, flags_extra(op.flags), op.force, obj_dir)
    informar_compilacion("alcp_native", r)
    if not r.ok:
        if "Permission denied" in r.salida:
            UI.nota("¿Hay algún proceso Python con el módulo cargado? Ciérralo y repite.")
        UI.fin_seccion()
        return Paso("module", "FALLO", time.perf_counter() - t0, "no compila")
    ri = importar_con_reintento(PAQUETE, lambda: construir_modulo(
        salida, perfil, flags_extra(op.flags), True, obj_dir))
    UI.linea("OK" if ri.ok else "FALLO", "importar alcp_native", ri.segundos)
    UI.bloque(ri.salida)
    UI.fin_seccion()
    return Paso("module", "OK" if ri.ok else "FALLO", time.perf_counter() - t0,
                "" if ri.ok else "no se puede importar")


# ---------------------------------------------------------------------------
# Subcomando: pytest
# ---------------------------------------------------------------------------
def hay_tests_python() -> bool:
    return (BINDINGS / "tests").is_dir()


def paso_pytest(args: list[str], con_cobertura: bool = False) -> Paso:
    if not args and not hay_tests_python():
        return Paso("pytest", "SALTADO", 0, "no existe python-bindings/tests")
    if importlib.util.find_spec("pytest") is None:
        raise ErrorEntorno(f"pytest no está instalado en {sys.executable}")
    UI.seccion("pytest (configuración en pyproject.toml)")
    extra = []
    if con_cobertura:
        if importlib.util.find_spec("pytest_cov") is not None:
            extra = ["--cov", "--cov-report=term-missing"]
        else:
            UI.nota("pytest-cov no está instalado: se ejecuta sin medir cobertura")
    r = correr([sys.executable, "-m", "pytest", *extra, *args], capturar=False,
               env=env_ejecucion(PYTHONIOENCODING="utf-8"))
    UI.linea("OK" if r.ok else "FALLO", "pytest", r.segundos,
             "" if r.ok else f"código {r.rc}")
    UI.fin_seccion()
    return Paso("pytest", "OK" if r.ok else "FALLO", r.segundos)


# ---------------------------------------------------------------------------
# Subcomando: cov-cpp
# ---------------------------------------------------------------------------
def ejecutable_gcov() -> str:
    """gcov a juego con el compilador (g++-14 -> gcov-14), o $GCOV.

    Se devuelve con «/» como separador: gcovr trocea --gcov-executable con
    shlex (admite p. ej. "llvm-cov gcov") y se comería las «\\» de Windows.
    """
    if os.environ.get("GCOV"):
        return os.environ["GCOV"]
    cxx = shutil.which(compilador())
    if cxx:
        p = Path(cxx)
        nombre = p.name.replace("g++", "gcov") if "g++" in p.name else "gcov" + p.suffix
        candidato = p.with_name(nombre)
        if candidato.exists():
            return candidato.as_posix()
    if shutil.which("gcov"):
        return "gcov"
    raise ErrorEntorno("no se encuentra gcov (define la variable GCOV)")


def asegurar_gcovr() -> None:
    if importlib.util.find_spec("gcovr") is not None:
        return
    if sys.prefix == sys.base_prefix:
        raise ErrorEntorno("gcovr no está instalado y no estamos en un entorno virtual; "
                           f"instálalo con: {sys.executable} -m pip install gcovr")
    UI.nota(f"gcovr no está instalado: se instala en el entorno virtual ({sys.prefix})")
    r = correr([sys.executable, "-m", "pip", "install", "--disable-pip-version-check",
                "gcovr"], timeout=600)
    if not r.ok:
        UI.bloque(r.salida)
        raise ErrorEntorno("no se pudo instalar gcovr con pip")


def con_modulo_temporal(instrumentado: Path, accion):
    """Coloca temporalmente el módulo instrumentado en python-bindings/alcp/,
    ejecuta `accion()` y restaura SIEMPRE el módulo normal (o lo retira si no
    había ninguno)."""
    destino = PAQUETE / instrumentado.name
    respaldo = COV / (instrumentado.name + ".normal")
    try:
        if destino.exists():
            os.replace(destino, respaldo)
        shutil.copy2(instrumentado, destino)
    except OSError as ex:
        UI.error(f"no se pudo colocar el módulo instrumentado ({ex}); "
                 "¿hay algún proceso Python con el módulo cargado?")
        if respaldo.exists() and not destino.exists():
            os.replace(respaldo, destino)
        return Paso("pytest", "FALLO", 0, "módulo bloqueado")
    try:
        return accion()
    finally:
        try:
            borrar(destino)                   # con reintentos (antivirus)
            if respaldo.exists():
                os.replace(respaldo, destino)
                UI.nota(f"módulo normal restaurado en {destino.relative_to(RAIZ).as_posix()}")
        except OSError as ex:
            UI.error(f"no se pudo restaurar el módulo normal ({ex}); "
                     "ejecuta «python build.py module»")


def paso_cov_cpp(op: Opciones, timeout: float = 1800) -> Paso:
    comprobar_compilador()
    gcov = ejecutable_gcov()
    asegurar_gcovr()
    t0 = time.perf_counter()
    UI.seccion("Cobertura C++: compilación instrumentada en build/cov/")
    # Siempre desde cero: .gcno y .gcda de compilaciones distintas no casan.
    if COV.exists():
        borrar(COV)
    obj, bin_tests, mod_dir = COV / "obj", COV / "tests", COV / "module"
    for d in (obj, bin_tests, mod_dir):
        d.mkdir(parents=True, exist_ok=True)
    extra = flags_extra(op.flags)
    cov_f = ["-O0", "-g", "--coverage"]
    deps = cabeceras()

    def enlazar_test(src: Path, exe: Path) -> Resultado:
        o = obj / (src.stem + ".o")
        return compilar([compilador(), "--coverage", str(o), "-o", str(exe), *LIBS_EXE],
                        exe, [o], True)

    def construir_test(src: Path, exe: Path) -> Resultado:
        # En dos pasos: así .gcno/.gcda quedan en build/cov/obj/<test>.*
        o = obj / (src.stem + ".o")
        r1 = compilar([compilador(), *FLAGS_BASE, "-I", str(INCLUDE), *extra, *cov_f,
                       "-c", str(src), "-o", str(o)], o, [src, *deps], True)
        if not r1.ok:
            return r1
        r2 = enlazar_test(src, exe)
        r2.segundos += r1.segundos
        r2.salida = (r1.salida + "\n" + r2.salida).strip()
        r2.avisos = bool(r2.salida)
        return r2

    fuentes = {nombre_test(s): s for s in fuentes_tests()}
    exes = {n: bin_tests / (n + EXE) for n in fuentes}
    trabajos = [(n, lambda n=n: construir_test(fuentes[n], exes[n])) for n in fuentes]
    modulo = mod_dir / ("alcp_native" + sufijo_modulo())
    hay_modulo = importlib.util.find_spec("pybind11") is not None
    if hay_modulo:
        trabajos.append(("alcp_native", lambda: construir_modulo(
            modulo, "coverage", extra, True, obj)))
    else:
        UI.linea("SALTADO", "compilar alcp_native", detalle="pybind11 no está instalado")
    compilados = compilar_varios(trabajos, op.jobs)
    fallos = [n for n, r in compilados.items() if not r.ok]
    UI.fin_seccion()

    UI.seccion("Cobertura C++: ejecución de los tests instrumentados")
    if "xcheck" in exes and not hay_vectores():
        UI.nota("generando build/oracle_vectors.txt para xcheck")
        rv = correr([sys.executable, str(BINDINGS / "oracle.py"), str(VECTORES)],
                    env=env_ejecucion(PYTHONDONTWRITEBYTECODE="1"), timeout=600)
        if not rv.ok:
            UI.bloque(rv.salida)
    for nombre, exe in exes.items():
        if nombre in fallos:
            continue
        cmd = [str(exe)]
        if nombre == "xcheck":
            if not hay_vectores():
                UI.linea("SALTADO", f"ejecutar {nombre}", detalle="sin vectores del oráculo")
                continue
            cmd.append(str(VECTORES))
        if not ejecutar_test(nombre, cmd, timeout,
                             lambda n=nombre: enlazar_test(fuentes[n], exes[n])):
            fallos.append(nombre)
    if hay_modulo and "alcp_native" not in fallos:
        # Antes de pytest: el módulo instrumentado debe poder cargarse (si se
        # bloquea, basta con volver a enlazarlo).
        enlace = cmds_modulo(modulo, "coverage", extra, obj)[-1]
        ri = importar_con_reintento(mod_dir, lambda: compilar(enlace, modulo, [], True))
        UI.linea("OK" if ri.ok else "FALLO", "importar alcp_native", ri.segundos)
        if not ri.ok:
            UI.bloque(ri.salida)
            fallos.append("alcp_native")
        elif not hay_tests_python():
            UI.linea("SALTADO", "pytest", detalle="no existe python-bindings/tests")
    UI.fin_seccion()

    # Código Python que ejercita el módulo instrumentado: la suite pytest
    # (incluye python-bindings/tests/test_native.py, que sustituye al antiguo
    # test_module.py).
    if hay_modulo and "alcp_native" not in fallos and hay_tests_python():
        if con_modulo_temporal(modulo, lambda: paso_pytest([])).estado == "FALLO":
            fallos.append("pytest")

    UI.seccion("Cobertura C++: informe gcovr -> build/coverage-cpp/")
    if INFORME_COV.exists():
        borrar(INFORME_COV)
    INFORME_COV.mkdir(parents=True)
    resumen = INFORME_COV / "resumen.txt"
    cmd = [sys.executable, "-m", "gcovr",
           "--root", str(RAIZ),
           "--object-directory", str(obj),
           "--gcov-executable", gcov,
           # Los filtros relativos se interpretan respecto al cwd (RAIZ) y
           # usan siempre «/» como separador, también en Windows.
           "--filter", "include/alcp/",
           "--filter", r"python-bindings/bindings\.cpp",
           # Las cabeceras se instancian en varias unidades de traducción.
           "--merge-mode-functions", "separate",
           "--exclude-throw-branches", "--exclude-unreachable-branches",
           "--html-details", str(INFORME_COV / "index.html"),
           "--html-title", "Cobertura C++ de alcp",
           "--txt", str(resumen),
           "--print-summary",
           str(obj)]
    r = correr(cmd, env=env_ejecucion(), timeout=timeout)
    UI.linea("OK" if r.ok else "FALLO", "gcovr", r.segundos)
    if r.ok and resumen.exists():
        UI.bloque(resumen.read_text(encoding="utf-8", errors="replace"), None)
    UI.bloque(r.salida, None if r.ok else 60)
    if r.ok:
        UI.nota(f"informe HTML: {(INFORME_COV / 'index.html').as_posix()}")
    else:
        fallos.append("gcovr")
    UI.fin_seccion()
    seg = time.perf_counter() - t0
    if fallos:
        return Paso("cov-cpp", "FALLO", seg, "fallan: " + ", ".join(fallos))
    return Paso("cov-cpp", "OK", seg, ultima_linea_cobertura(r.salida))


def ultima_linea_cobertura(salida: str) -> str:
    """Extrae «lines: NN.N% (...)» del resumen de gcovr para el resumen final."""
    for ln in salida.splitlines():
        if ln.strip().lower().startswith("lines:"):
            return "líneas " + ln.split(":", 1)[1].strip()
    return ""


# ---------------------------------------------------------------------------
# Subcomando: bench
# ---------------------------------------------------------------------------
def paso_bench(args: list[str], op: Opciones) -> Paso:
    script = RAIZ / "bench" / "run_bench.py"
    if script.exists():
        UI.seccion("Benchmarks: bench/run_bench.py")
        r = correr([sys.executable, str(script), *args], capturar=False,
                   env=env_ejecucion(PYTHONIOENCODING="utf-8"))
        UI.fin_seccion()
        return Paso("bench", "OK" if r.ok else "FALLO", r.segundos)
    # Alternativa heredada de build.sh: el benchmark C++ original.
    fuente = RAIZ / "bench" / "bench.cpp"
    if not fuente.exists():
        return Paso("bench", "SALTADO", 0, "no existe bench/run_bench.py ni bench/bench.cpp")
    comprobar_compilador()
    t0 = time.perf_counter()
    UI.seccion("Benchmark C++ (bench/bench.cpp, -O3 -march=native)")
    exe = BUILD / "bench" / ("bench" + EXE)
    cmd = [compilador(), "-std=c++23", "-O3", "-march=native", "-Wall", "-Wextra",
           "-I", str(INCLUDE), *flags_extra(op.flags), str(fuente), "-o", str(exe), *LIBS_EXE]
    rc = compilar(cmd, exe, [fuente, *cabeceras()], op.force)
    informar_compilacion("bench", rc)
    if not rc.ok:
        UI.fin_seccion()
        return Paso("bench", "FALLO", time.perf_counter() - t0, "no compila")
    r = correr([str(exe), *args], env=env_ejecucion(), capturar=False, reintentar_av=True)
    UI.fin_seccion()
    return Paso("bench", "OK" if r.ok else "FALLO", time.perf_counter() - t0)


# ---------------------------------------------------------------------------
# Subcomando: clean
# ---------------------------------------------------------------------------
def paso_clean() -> Paso:
    UI.seccion("Limpieza")
    objetivos = [BUILD, *PAQUETE.glob("alcp_native*.pyd"), *PAQUETE.glob("alcp_native*.so"),
                 # datos que generan showcase.py y pytest/coverage
                 BINDINGS / "nemo_in.txt", BINDINGS / "nemo_out.txt",
                 RAIZ / ".pytest_cache", RAIZ / ".coverage",
                 # restos del build.sh antiguo
                 *TESTS.glob("*.exe"), *TESTS.glob("*.o"), *(RAIZ / "bench").glob("*.exe")]
    errores = 0
    for p in objetivos:
        if not p.exists():
            continue
        try:
            borrar(p)
            UI.linea("OK", "borrado", detalle=p.relative_to(RAIZ).as_posix())
        except OSError as ex:
            errores += 1
            UI.linea("FALLO", "borrar", detalle=f"{p.relative_to(RAIZ).as_posix()}: {ex}")
    UI.fin_seccion()
    return Paso("clean", "FALLO" if errores else "OK")


# ---------------------------------------------------------------------------
# Línea de órdenes
# ---------------------------------------------------------------------------
SIN_ANALIZAR = ("pytest", "bench")      # su resto de argumentos se pasa tal cual


def construir_parser() -> argparse.ArgumentParser:
    comunes = argparse.ArgumentParser(add_help=False)
    # SUPPRESS: así la opción vale tanto antes como después del subcomando.
    comunes.add_argument("-v", "--verbose", action="store_true", default=argparse.SUPPRESS,
                         help="muestra los comandos que se ejecutan")
    comunes.add_argument("--no-color", action="store_true", default=argparse.SUPPRESS,
                         help="desactiva los colores")
    compila = argparse.ArgumentParser(add_help=False)
    compila.add_argument("--flags", action="append", default=[], metavar='"..."',
                         help="flags extra del compilador (repetible)")
    compila.add_argument("-j", "--jobs", type=int, default=os.cpu_count() or 2,
                         help="compilaciones en paralelo")
    compila.add_argument("-B", "--force", action="store_true",
                         help="recompila aunque esté al día")

    p = argparse.ArgumentParser(
        prog="build.py", parents=[comunes],
        description="Compilación y verificación de alcp (Windows MinGW y Linux).",
        epilog="Detalles en la cabecera de build.py.")
    sub = p.add_subparsers(dest="orden", metavar="subcomando", required=True)
    sub.add_parser("tests", parents=[comunes, compila], help="compila y ejecuta tests/*.cpp")
    sub.add_parser("oracle", parents=[comunes, compila], help="vectores del oráculo + xcheck")
    m = sub.add_parser("module", parents=[comunes, compila], help="compila el módulo pybind11")
    grupo = m.add_mutually_exclusive_group()
    grupo.add_argument("--debug", action="store_true", help="-O0 -g")
    grupo.add_argument("--coverage", action="store_true", help="-O0 -g --coverage")
    sub.add_parser("pytest", parents=[comunes], help="ejecuta pytest (argumentos tal cual)")
    sub.add_parser("cov-cpp", parents=[comunes, compila], help="cobertura C++ con gcovr")
    sub.add_parser("bench", parents=[comunes], help="benchmarks (argumentos tal cual)")
    sub.add_parser("all", parents=[comunes, compila], help="tests + oracle + module + pytest")
    sub.add_parser("ci", parents=[comunes, compila], help="como all, con salida para CI")
    sub.add_parser("clean", parents=[comunes], help="borra build/ y los módulos")
    return p


def separar_argumentos(argv: list[str]) -> tuple[list[str], list[str]]:
    """Para pytest/bench, todo lo que sigue al subcomando se pasa sin analizar.

    Además, «--flags -O2» se reescribe como «--flags=-O2»: argparse tomaría
    el valor, que empieza por «-», como si fuera otra opción.
    """
    propios, directos = list(argv), []
    for i, a in enumerate(argv):
        if not a.startswith("-"):
            if a in SIN_ANALIZAR:
                propios, directos = list(argv[:i + 1]), list(argv[i + 1:])
            break
    i = 0
    while i < len(propios) - 1:
        if propios[i] == "--flags":
            propios[i:i + 2] = [f"--flags={propios[i + 1]}"]
        i += 1
    return propios, directos


def resumen_final(pasos: list[Paso], t0: float) -> int:
    # El resumen no se pliega en GitHub Actions: debe verse siempre.
    UI.seccion("Resumen", plegable=False)
    for paso in pasos:
        seg = None if paso.estado == "SALTADO" else paso.segundos
        UI.linea(paso.estado, paso.nombre, seg, paso.detalle)
    total = time.perf_counter() - t0
    fallidos = [p.nombre for p in pasos if p.estado == "FALLO"]
    if fallidos:
        UI.error(f"FALLO en: {', '.join(fallidos)} ({total:.1f} s)")
        return 1
    print(f"\n  TODO CORRECTO ({total:.1f} s)", flush=True)
    return 0


def main(argv: list[str]) -> int:
    # Salida en UTF-8 también cuando stdout es una tubería (Git Bash, CI), y
    # sin fallar nunca por un carácter no representable.
    for flujo in (sys.stdout, sys.stderr):
        try:
            flujo.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass
    propios, paso_directo = separar_argumentos(argv)
    a = construir_parser().parse_args(propios)
    UI.configurar(getattr(a, "verbose", False), getattr(a, "no_color", False),
                  a.orden == "ci")
    op = Opciones(flags=getattr(a, "flags", []), jobs=getattr(a, "jobs", os.cpu_count() or 2),
                  force=getattr(a, "force", False))
    t0 = time.perf_counter()
    if UI.verbose:
        UI.nota(f"python: {sys.executable} ({sys.version.split()[0]}); "
                f"CXX: {compilador()}; raíz: {RAIZ}")
    pasos: list[Paso] = []
    try:
        if a.orden == "tests":
            pasos.append(paso_tests(op))
        elif a.orden == "oracle":
            pasos.append(paso_oracle(op))
        elif a.orden == "module":
            perfil = "coverage" if a.coverage else "debug" if a.debug else "release"
            pasos.append(paso_module(op, perfil))
        elif a.orden == "pytest":
            pasos.append(paso_pytest(paso_directo))
        elif a.orden == "cov-cpp":
            pasos.append(paso_cov_cpp(op))
        elif a.orden == "bench":
            pasos.append(paso_bench(paso_directo, op))
        elif a.orden == "clean":
            pasos.append(paso_clean())
        elif a.orden in ("all", "ci"):
            # Cada paso se ejecuta aunque falle el anterior (para ver todos
            # los problemas de una vez); un error de entorno cuenta como FALLO.
            def paso(nombre, funcion):
                try:
                    pasos.append(funcion())
                except ErrorEntorno as ex:
                    UI.error(str(ex))
                    pasos.append(Paso(nombre, "FALLO", 0, str(ex)))
            # xcheck lo ejecuta el paso oracle con vectores recién generados.
            paso("tests", lambda: paso_tests(op, incluir_xcheck=False))
            # (sin -B: xcheck se acaba de compilar en el paso tests)
            paso("oracle", lambda: paso_oracle(replace(op, force=False)))
            paso("module", lambda: paso_module(op))
            if pasos[-1].estado == "OK":
                paso("pytest", lambda: paso_pytest([], con_cobertura=True))
            else:
                pasos.append(Paso("pytest", "SALTADO", 0, "el módulo no compila"))
    except ErrorEntorno as ex:
        UI.error(str(ex))
        return 2
    except KeyboardInterrupt:
        UI.error("interrumpido por el usuario")
        return 130
    return resumen_final(pasos, t0)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
