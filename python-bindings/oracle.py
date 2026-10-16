"""
Oráculo de verificación cruzada (C++ frente a las prácticas Python originales).

Ejecuta las implementaciones Python ORIGINALES (python/cuerpos_finitos.py,
python/factorizacion.py, python/algoritmos_rapidos.py — Prácticas 1/2/3) sobre
entradas pseudoaleatorias con semilla fija y escribe, un registro por línea,
las entradas y las salidas esperadas.  El programa tests/crosscheck.cpp (que se
compila como «xcheck») lee este archivo, recalcula cada resultado con el port
C++ y comprueba que coincide exactamente con la referencia Python (para la
factorización, que es aleatorizada, compara el multiconjunto canónico de
factores).

Las prácticas originales NO se modifican: solo se importan.  Importar
factorizacion / algoritmos_rapidos aplica «monkey-patching» en memoria sobre
las clases de cuerpos_finitos (añade factorizar, mult_fast, divmod_fast...),
que es justo lo que se pretende.

Uso:
    python python-bindings/oracle.py [salida]      (por defecto build/oracle_vectors.txt)

Gramática de los registros (campos separados por '|'):
    OP | field | in... | out...
donde
    field  : "p"            para F_p
             "p;g0,g1,..."  para F_q = F_p[x]/(g), g little-endian y mónico
    poly   : coeficientes separados por comas (little-endian); en F_p son
             enteros, en F_q son índices en base p (== conv_a_int)
    factor : "coefs:mult;coefs:mult;..." en orden canónico

Registros de factorización (FPX_FACTOR / FQX_FACTOR):
    OP | field | poly | factor | py_ok
    py_ok = 1 si la referencia es autoconsistente (el coeficiente principal por
    el producto de los factores reconstruye la entrada y todos los factores son
    irreducibles); py_ok = 0 en caso contrario.  Si la referencia LANZA una
    excepción, el campo factor queda VACÍO y py_ok = 0; el C++ se verifica
    entonces por reconstrucción + irreducibilidad.

La salida es determinista (semilla fija, orden fijo, finales de línea LF).  Al
terminar se imprime por stderr un resumen con el número de registros y los
fallos de la referencia agrupados por motivo.
"""
import os
import random
import sys
from collections import Counter, defaultdict

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)
SALIDA_POR_DEFECTO = os.path.join(RAIZ, "build", "oracle_vectors.txt")

# Las prácticas originales viven en ../python.  No se escriben .pyc junto a
# ellas para no ensuciar ese directorio (es de solo lectura para este proyecto).
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(RAIZ, "python"))
import cuerpos_finitos as cf            # noqa: E402
import factorizacion                    # noqa: E402,F401  parchea anillo_*.factorizar (Cantor-Zassenhaus)
import algoritmos_rapidos as ar         # noqa: E402       parchea mult_fast / divmod_fast, FFT, Toeplitz

random.seed(20260629)
OUT = []

# Estadísticas para el resumen final (stderr).
CUENTA_OPS = Counter()                  # registros emitidos por operación
FALLOS_REF = Counter()                  # (op, motivo) -> veces
CAMPOS_FALLO = defaultdict(set)         # (op, motivo) -> cuerpos en los que ocurrió


def emit(*tokens):
    CUENTA_OPS[str(tokens[0])] += 1
    OUT.append("|".join(str(t) for t in tokens))


def anotar_fallo(op, campo, motivo):
    """Registra un fallo de la referencia Python para el resumen final."""
    FALLOS_REF[(op, motivo)] += 1
    CAMPOS_FALLO[(op, motivo)].add(str(campo))


# ---- F_p helpers -----------------------------------------------------------
def poly_fp(fpx, ints):
    return fpx.elem_de_tuple(tuple(ints))

def coeffs_fp(fpx, poly):
    return [int(c) for c in fpx.conv_a_tuple(poly)]

def s_poly(ints):
    return ",".join(str(i) for i in ints) if ints else "0"

def rand_poly_ints(p, deg):
    if deg < 0:
        return [0]
    c = [random.randrange(p) for _ in range(deg)]
    c.append(random.randrange(1, p))   # nonzero leading
    return c

def canon_factors_fp(fpx, faclist):
    items = []
    for (u, e) in faclist:
        items.append((tuple(coeffs_fp(fpx, u)), e))
    items.sort()
    return ";".join(",".join(str(x) for x in cs) + ":" + str(e) for (cs, e) in items)


# ---- F_q helpers -----------------------------------------------------------
def poly_fq(fqx, fq, idxs):
    return fqx.elem_de_tuple(tuple(fq.elem_de_int(i) for i in idxs))

def coeffs_fq(fqx, fq, poly):
    return [fq.conv_a_int(c) for c in fqx.conv_a_tuple(poly)]

def rand_poly_fq_idxs(q, deg):
    if deg < 0:
        return [0]
    c = [random.randrange(q) for _ in range(deg)]
    c.append(random.randrange(1, q))
    return c

def canon_factors_fq(fqx, fq, faclist):
    items = []
    for (u, e) in faclist:
        items.append((tuple(coeffs_fq(fqx, fq, u)), e))
    items.sort()
    return ";".join(",".join(str(x) for x in cs) + ":" + str(e) for (cs, e) in items)


# ---- factorización de referencia, tolerante a fallos -----------------------
def factor_referencia(op, campo, anillo, comp, cI, constante, coeffs, canon):
    """Factoriza `comp` con la práctica original y devuelve (factor, py_ok).

    - `constante(c)` construye el polinomio constante c (c = coeficiente
      principal codificado como en el registro); `coeffs(poly)` devuelve la
      codificación del polinomio; `canon(faclist)` el campo factor canónico.
    - La reconstrucción multiplica por el coeficiente principal, porque la
      práctica devuelve factores MÓNICOS y la entrada no tiene por qué serlo.
    - Si la referencia lanza una excepción se devuelve ("", 0): el registro
      conserva el formato con el campo de factores vacío.
    """
    try:
        faclist = anillo.factorizar(comp)
        rec = constante(cI[-1])
        for (u, e) in faclist:
            rec = anillo.mult(rec, u ** e)
        reconstruye = coeffs(rec) == cI
        irreducibles = all(u.es_irreducible for (u, _) in faclist)
        factor = canon(faclist)
    except Exception as ex:   # la práctica original tiene fallos conocidos
        anotar_fallo(op, campo, f"excepción {type(ex).__name__}: {ex}")
        return "", 0
    # Un único motivo por registro (puede combinar los dos defectos).
    defectos = []
    if not irreducibles:
        defectos.append("devuelve factores reducibles")
    if not reconstruye:
        defectos.append("no reconstruye la entrada")
    if defectos:
        anotar_fallo(op, campo, " y ".join(defectos))
    return factor, (0 if defectos else 1)


# ===========================================================================
# F_p[x] cases
# ===========================================================================
def casos_fp():
    FP_PRIMES = [2, 3, 5, 7, 13, 101]
    for p in FP_PRIMES:
        fp = cf.cuerpo_fp(p)
        fpx = cf.anillo_fp_x(fp)
        for _ in range(12):
            fI = rand_poly_ints(p, random.randint(0, 7))
            gI = rand_poly_ints(p, random.randint(0, 5))
            f, g = poly_fp(fpx, fI), poly_fp(fpx, gI)

            emit("FPX_MUL", p, s_poly(fI), s_poly(gI), s_poly(coeffs_fp(fpx, fpx.mult(f, g))))
            emit("FPX_KARA", p, s_poly(fI), s_poly(gI), s_poly(coeffs_fp(fpx, fpx.mult_fast(f, g))))

            if cf.PolyK(fpx, gI).grado >= 0 and any(gI):
                q, r = fpx.divmod(f, g)
                emit("FPX_DIVMOD", p, s_poly(fI), s_poly(gI),
                     s_poly(coeffs_fp(fpx, q)), s_poly(coeffs_fp(fpx, r)))
                q2, r2 = fpx.divmod_fast(f, g)
                emit("FPX_DIVMOD_FAST", p, s_poly(fI), s_poly(gI),
                     s_poly(coeffs_fp(fpx, q2)), s_poly(coeffs_fp(fpx, r2)))

            emit("FPX_GCD", p, s_poly(fI), s_poly(gI), s_poly(coeffs_fp(fpx, fpx.gcd(f, g))))
            gg, xx, yy = fpx.gcd_ext(f, g)
            emit("FPX_GCDEXT", p, s_poly(fI), s_poly(gI),
                 s_poly(coeffs_fp(fpx, gg)), s_poly(coeffs_fp(fpx, xx)), s_poly(coeffs_fp(fpx, yy)))

        # irreducibility + factorization
        for _ in range(10):
            fI = rand_poly_ints(p, random.randint(1, 5))
            f = poly_fp(fpx, fI)
            emit("FPX_IRRED", p, s_poly(fI), 1 if f.es_irreducible else 0)
        for _ in range(8):
            # build a known composite from random factors
            comp = poly_fp(fpx, [1])
            for _ in range(random.randint(1, 3)):
                d = random.randint(1, 3)
                e = random.randint(1, 2)
                base = poly_fp(fpx, rand_poly_ints(p, d))
                comp = fpx.mult(comp, base ** e)
            cI = coeffs_fp(fpx, comp)
            if fpx.grado(comp) >= 1:
                factor, py_ok = factor_referencia(
                    "FPX_FACTOR", p, fpx, comp, cI,
                    constante=lambda c: poly_fp(fpx, [c]),
                    coeffs=lambda pol: coeffs_fp(fpx, pol),
                    canon=lambda fl: canon_factors_fp(fpx, fl))
                emit("FPX_FACTOR", p, s_poly(cI), factor, py_ok)

    # pot_mod over a fixed prime
    fp = cf.cuerpo_fp(7); fpx = cf.anillo_fp_x(fp)
    for _ in range(8):
        aI = rand_poly_ints(7, random.randint(1, 4))
        mI = rand_poly_ints(7, random.randint(1, 3))
        k = random.randint(0, 25)
        a, m = poly_fp(fpx, aI), poly_fp(fpx, mI)
        if fpx.grado(m) >= 1:
            emit("FPX_POWMOD", 7, s_poly(aI), s_poly(mI), k, s_poly(coeffs_fp(fpx, fpx.pot_mod(a, k, m))))

    # DFT / IDFT over F_17 (16 = 2^4)
    fp17 = cf.cuerpo_fp(17); k = 4; n = 1 << k; groot = 3   # 3 has order 16 in F_17*
    for _ in range(6):
        aI = [random.randrange(17) for _ in range(n)]
        a_elems = tuple(fp17.elem(v) for v in aI)
        A = ar.fp_fft(fp17, fp17.elem(groot), k, a_elems)
        emit("FP_DFT", 17, groot, k, s_poly(aI), s_poly([int(v) for v in A]))
        B = ar.fp_ifft(fp17, fp17.elem(groot), k, A)
        emit("FP_IDFT", 17, groot, k, s_poly([int(v) for v in A]), s_poly([int(v) for v in B]))

    # Toeplitz lower mat-vec over F_101
    fp101 = cf.cuerpo_fp(101)
    for _ in range(6):
        nn = random.randint(2, 10)
        a = tuple(fp101.elem(random.randrange(101)) for _ in range(nn))
        b = tuple(fp101.elem(random.randrange(101)) for _ in range(nn))
        y = ar.fp_toep_inf_vec(fp101, nn, a, b)
        emit("FP_TOEP_LOWER", 101, nn, s_poly([int(v) for v in a]),
             s_poly([int(v) for v in b]), s_poly([int(v) for v in y]))


# ===========================================================================
# F_q[x] cases (fixed irreducible moduli shared with the C++ side)
# ===========================================================================
def casos_fq():
    FQ_FIELDS = [
        (3, [1, 0, 1]),       # F_9  = F_3[x]/(x^2+1)
        (3, [1, 2, 0, 1]),    # F_27 = F_3[x]/(x^3+2x+1)
        (5, [2, 0, 1]),       # F_25 = F_5[x]/(x^2+2)
    ]
    for (p, gmod) in FQ_FIELDS:
        fp = cf.cuerpo_fp(p)
        fq = cf.cuerpo_fq(fp, tuple(gmod))
        fqx = cf.anillo_fq_x(fq)
        q = fq.fq
        fld = f"{p};{','.join(str(c) for c in gmod)}"
        for _ in range(10):
            fI = rand_poly_fq_idxs(q, random.randint(0, 5))
            gI = rand_poly_fq_idxs(q, random.randint(0, 4))
            f, g = poly_fq(fqx, fq, fI), poly_fq(fqx, fq, gI)
            emit("FQX_MUL", fld, s_poly(fI), s_poly(gI), s_poly(coeffs_fq(fqx, fq, fqx.mult(f, g))))
            if fqx.grado(g) >= 0 and any(gI):
                qq, rr = fqx.divmod(f, g)
                emit("FQX_DIVMOD", fld, s_poly(fI), s_poly(gI),
                     s_poly(coeffs_fq(fqx, fq, qq)), s_poly(coeffs_fq(fqx, fq, rr)))
            emit("FQX_GCD", fld, s_poly(fI), s_poly(gI), s_poly(coeffs_fq(fqx, fq, fqx.gcd(f, g))))
        for _ in range(6):
            fI = rand_poly_fq_idxs(q, random.randint(1, 4))
            f = poly_fq(fqx, fq, fI)
            emit("FQX_IRRED", fld, s_poly(fI), 1 if f.es_irreducible else 0)
        for _ in range(5):
            comp = poly_fq(fqx, fq, [1])
            for _ in range(random.randint(1, 2)):
                d = random.randint(1, 2)
                e = random.randint(1, 2)
                base = poly_fq(fqx, fq, rand_poly_fq_idxs(q, d))
                comp = fqx.mult(comp, base ** e)
            if fqx.grado(comp) >= 1:
                cI = coeffs_fq(fqx, fq, comp)
                factor, py_ok = factor_referencia(
                    "FQX_FACTOR", fld, fqx, comp, cI,
                    constante=lambda c: poly_fq(fqx, fq, [c]),
                    coeffs=lambda pol: coeffs_fq(fqx, fq, pol),
                    canon=lambda fl: canon_factors_fq(fqx, fq, fl))
                emit("FQX_FACTOR", fld, s_poly(cI), factor, py_ok)


def resumen():
    """Resumen legible por stderr (no forma parte del archivo de vectores)."""
    err = sys.stderr
    print(f"oracle: {len(OUT)} registros", file=err)
    for op in sorted(CUENTA_OPS):
        print(f"  {op:<16} {CUENTA_OPS[op]:5d}", file=err)
    n_fact = CUENTA_OPS["FPX_FACTOR"] + CUENTA_OPS["FQX_FACTOR"]
    n_fallo = sum(1 for linea in OUT if "_FACTOR|" in linea and linea.endswith("|0"))
    print(f"oracle: la referencia Python falla en {n_fallo} de {n_fact} factorizaciones "
          "(el C++ se verifica en ellas por reconstrucción + irreducibilidad)", file=err)
    for (op, motivo), veces in sorted(FALLOS_REF.items(), key=lambda kv: (kv[0][0], -kv[1], kv[0][1])):
        campos = ", ".join(sorted(CAMPOS_FALLO[(op, motivo)], key=lambda s: (len(s), s)))
        print(f"  {op:<11} x{veces:<3} {motivo}  [cuerpos: {campos}]", file=err)


def main(argv):
    out_path = argv[1] if len(argv) > 1 else SALIDA_POR_DEFECTO
    casos_fp()
    casos_fq()
    carpeta = os.path.dirname(os.path.abspath(out_path))
    os.makedirs(carpeta, exist_ok=True)
    # newline="\n": finales de línea LF también en Windows (salida reproducible).
    with open(out_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(OUT) + "\n")
    print(f"wrote {len(OUT)} cross-check records to {out_path}", flush=True)
    resumen()
    return 0


if __name__ == "__main__":
    # Salida UTF-8 también cuando stdout/stderr son tuberías en Windows
    # (mismo criterio que showcase.py), para que los acentos salgan bien.
    for flujo in (sys.stdout, sys.stderr):
        try:
            flujo.reconfigure(encoding="utf-8")
        except Exception:
            pass
    sys.exit(main(sys.argv))
