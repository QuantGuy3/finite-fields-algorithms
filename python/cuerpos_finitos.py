"""
PRÁCTICA ALCP 1

Grupo: Ismael Amador, Ángela Ruiz, Miquel Sandonís, Elia Torres y Claudia Vicario
"""



import unittest
import random
from itertools import zip_longest
from typing import List, Tuple, Union, Optional
import builtins
import re

# =========================
# Utilidades básicas
# =========================

def es_primo(p: int) -> bool:
    if p < 2:
        return False
    i = 2
    while i * i <= p:
        if p % i == 0:
            return False
        i += 1
    return True


def factor_potencia(q: int) -> Tuple[int, int]:
    """
    Devuelve (p, n) si q = p^n con p primo.
    """
    if q < 2:
        raise ValueError("q debe ser ≥ 2")
    for p in range(2, int(q**0.5) + 2):
        if es_primo(p):
            n, x = 0, q
            while x % p == 0:
                x //= p
                n += 1
            if x == 1:
                return p, n
    if es_primo(q):
        return q, 1
    raise ValueError(f"{q} no es potencia de un primo")

# =========================================================
# Elementos de Fp + núcleo rígido cuerpo_fp
# =========================================================

class cuerpo_fp:
    """Núcleo rígido para Fp: opera con enteros normalizados [0, p-1]."""
    def __init__(self, p: int):
        if not es_primo(p):
            raise ValueError(f"{p} no es primo. No se puede construir Fp")
        self.p = p

    # --- Primitivas para el cálculo ---
    def _norm(self, a: int) -> int:
        return a % self.p

    def _add(self, a: int, b: int) -> int:
        return (a + b) % self.p

    def _sub(self, a: int, b: int) -> int:
        return (a - b) % self.p

    def _neg(self, a: int) -> int:
        return (-a) % self.p

    def _mul(self, a: int, b: int) -> int:
        return (a * b) % self.p

    def _pow(self, a: int, k: int) -> int:
        a = self._norm(a)
        if k < 0:
            if a == 0:
                raise ZeroDivisionError("0 no es invertible en Fp")
            # a^(p-2) = a^{-1}
            a = self._pow(a, self.p - 2)
            k = -k
        res = 1
        while k:
            if k & 1:
                res = (res * a) % self.p
            a = (a * a) % self.p
            k >>= 1
        return res
    
    def _inv(self, a: int) -> int:
        if a % self.p == 0:
            raise ZeroDivisionError("División por 0 en Fp")
        return self._pow(a, self.p - 2)

    # --- capa objetos opacos --- #
    def elem(self, n: int) -> "ElementoFp":
        return ElementoFp(self, n)
    def cero(self) -> "ElementoFp":
        return ElementoFp(self, 0)
    def uno(self) -> "ElementoFp":
        return ElementoFp(self, 1)
    def elem_de_int(self, n: int) -> "ElementoFp":
        return ElementoFp(self, n)
    def elem_de_str(self, s: str) -> "ElementoFp":
        return ElementoFp(self, int(s))
    def conv_a_int(self, a: "ElementoFp") -> int:
        return int(a)
    def conv_a_str(self, a: "ElementoFp") -> str:
        return str(int(a))
    def suma(self, a, b):
        a, b = ElementoFp(self, a), ElementoFp(self, b)
        return a + b
    def inv_adit(self, a):
        return -ElementoFp(self, a)
    def mult(self, a, b):
        a, b = ElementoFp(self, a), ElementoFp(self, b)
        return a * b
    def pot(self, a, k: int):
        return ElementoFp(self, a) ** k
    def inv_mult(self, a):
        return ElementoFp(self, a) ** -1
    def es_cero(self, a):
        return ElementoFp(self, a) == 0
    def es_uno(self, a):
        return ElementoFp(self, a) == 1
    def es_igual(self, a, b):
        return ElementoFp(self, a) == ElementoFp(self, b)
    def aleatorio(self) -> "ElementoFp":
        return self.elem(random.randrange(self.p))
    def tabla_suma(self):
        return [[(i + j) % self.p for j in range(self.p)] for i in range(self.p)]
    def tabla_mult(self):
        return [[(i * j) % self.p for j in range(self.p)] for i in range(self.p)]
    def tabla_inv_adit(self):
        return [(-i) % self.p for i in range(self.p)]
    def tabla_inv_mult(self):
        t = ["*"]
        for i in range(1, self.p):
            t.append(self._inv(i))
        return t
    def cuadrado_latino(self, a: "ElementoFp"):
        a = a if isinstance(a, ElementoFp) else self.elem_de_int(int(a))
        if self.es_cero(a):
            raise ValueError("cuadrado latino no definido para a=0")
        return [[int(a * self.elem_de_int(i) + self.elem_de_int(j)) for j in range(self.p)] for i in range(self.p)]


class ElementoFp:
    """Elemento opaco de Fp, con operadores sobrecargados
    (todos delegan en primitivas de cuerpo_fp, que trabaja con int).
    """
    __slots__ = ("P", "v")

    def __init__(self, parent: cuerpo_fp, x):
        self.P = parent
        p = parent.p
        if isinstance(x, ElementoFp):
            if x.P.p != p:
                raise TypeError("ElementoFp de cuerpo distinto.")
            self.v = x.v
        elif isinstance(x, int):
            self.v = x % p
        elif isinstance(x, str):
            s = x.strip()
            self.v = int(s) % p
        else:
            try:
                self.v = int(x) % p
            except Exception as e:
                raise TypeError(f"No se puede convertir {type(x)} a ElementoFp") from e

    # coerción
    def _coerce(self, other):
        return ElementoFp(self.P, other)

    # repr/conv
    def __int__(self):
        return self.v
    __index__ = __int__
    def __repr__(self):
        return f"Fp({self.v}; p={self.P.p})"
    def __str__(self):
        return str(self.v)
    def __hash__(self):
        return hash((id(self.P), self.v))

    # comparación
    def __eq__(self, other):
        o = self._coerce(other)
        return self.v == o.v

    # aritmética (delegando en cuerpo_fp)
    def __neg__(self):
        return ElementoFp(self.P, self.P._neg(self.v))

    def __add__(self, other):
        o = self._coerce(other)
        return ElementoFp(self.P, self.P._add(self.v, o.v))
    __radd__ = __add__

    def __sub__(self, other):
        o = self._coerce(other)
        return ElementoFp(self.P, self.P._sub(self.v, o.v))

    def __rsub__(self, other):
        o = self._coerce(other)
        return ElementoFp(self.P, self.P._sub(o.v, self.v))

    def __mul__(self, other):
        o = self._coerce(other)
        return ElementoFp(self.P, self.P._mul(self.v, o.v))
    __rmul__ = __mul__

    def __truediv__(self, other):
        o = self._coerce(other)
        return ElementoFp(self.P, self.P._mul(self.v, self.P._inv(o.v)))

    def __rtruediv__(self, other):
        o = self._coerce(other)
        return ElementoFp(self.P, self.P._mul(o.v, self.P._inv(self.v)))

    def __pow__(self, k: int):
        return ElementoFp(self.P, self.P._pow(self.v, k))

# =======================================================
# Polinomios genéricos PolyK (elemento) + anillos K[var]
# =======================================================

class PolyK:
    """Elemento opaco: polinomio sobre el cuerpo K con variable R.var (coeficientes opacos)."""
    __slots__ = ("R", "coef")

    def __init__(self, ring: "anillo_fp_x | anillo_fq_x", x):
        self.R = ring

        def coerce_list(seq):
            return [self.R._coerce_coef(c) for c in seq]

        # 1) Ya es PolyK
        if isinstance(x, PolyK):
            if x.R is self.R:
                c = list(x.coef)
            else:
                c = coerce_list(x.coef)
            self.coef = self._trim(c)
            return

        # 2) str: si contiene var o '^', parsea; si no, constante
        if isinstance(x, str):
            s = x.replace(" ", "")
            if self.R.var in s or "^" in s:
                poly = self.R.elem_de_str(s)
                self.coef = list(poly.coef)
                return
            else:
                self.coef = [self.R._coerce_coef(s)]
                self.coef = self._trim(self.coef)
                return

        # 3) lista/tupla de coeficientes
        if isinstance(x, (list, tuple)):
            c = coerce_list(list(x))
            if c == []:
                c = [self.R._coerce_coef(0)]
            self.coef = self._trim(c)
            return

        # 4) constante
        try:
            c0 = self.R._coerce_coef(x)
            self.coef = [c0]
            return
        except Exception:
            pass

        raise TypeError(f"No se puede convertir {type(x)} a PolyK sobre {self.R}")

    # util
    def _trim(self, c):
        while len(c) > 1 and c[-1] == 0:
            c.pop()
        return c

    def _coerce_poly(self, other):
        return PolyK(self.R, other)

    # como lista
    def __iter__(self):
        return iter(self.coef)
    def __len__(self):
        return len(self.coef)
    def __getitem__(self, i: int):
        if isinstance(i, slice):
            return type(self)(self.R, self.coef[i])
        return self.coef[i]
    def __setitem__(self, i, val):
        if isinstance(i, slice):
            vals = [self.R._coerce_coef(x) for x in val]
            self.coef[i] = vals
        else:
            self.coef[i] = self.R._coerce_coef(val)
        self.coef = self._trim(self.coef)
    def append(self, x):
        self.coef.append(self.R._coerce_coef(x))
        self.coef = self._trim(self.coef)

    # repr
    def __repr__(self):
        return f"{repr(self.R)}({self.R.conv_a_str(self)})"
    def __str__(self):
        return self.R.conv_a_str(self)
    def __hash__(self):
        return hash((id(self.R), tuple(self.coef)))

    # comparación
    def __eq__(self, other):
        o = self._coerce_poly(other)
        return self.coef == o.coef

    # grado
    @property
    def grado(self) -> int:
        for i in range(len(self.coef) - 1, -1, -1):
            if not (self.coef[i] == 0):
                return i
        return -1

    # -------- aritmética (delegación si hay kernel en el anillo) --------
    def __neg__(self):
        if hasattr(self.R, "inv_adit"):
            return self.R.inv_adit(self)
        return PolyK(self.R, [-a for a in self])

    def __add__(self, other):
        o = self._coerce_poly(other)
        if hasattr(self.R, "suma"):
            return self.R.suma(self, o)
        return PolyK(self.R, [a + b for a, b in zip_longest(self, o, fillvalue=0)])
    __radd__ = __add__

    def __sub__(self, other):
        o = self._coerce_poly(other)
        if hasattr(self.R, "suma") and hasattr(self.R, "inv_adit"):
            return self.R.suma(self, self.R.inv_adit(o))
        return PolyK(self.R, [a - b for a, b in zip_longest(self, o, fillvalue=0)])

    def __rsub__(self, other):
        o = self._coerce_poly(other)
        if hasattr(self.R, "suma") and hasattr(self.R, "inv_adit"):
            return self.R.suma(o, self.R.inv_adit(self))
        return PolyK(self.R, [a - b for a, b in zip_longest(o, self, fillvalue=0)])

    def __mul__(self, other):
        o = self._coerce_poly(other)
        if hasattr(self.R, "mult"):
            return self.R.mult(self, o)
        na, nb = len(self), len(o)
        if na == 0 or nb == 0:
            return PolyK(self.R, [])
        res = [0] * (na + nb - 1)
        for i in range(na):
            ai = self[i]
            if ai == 0:
                continue
            for j in range(nb):
                res[i + j] = res[i + j] + (ai * o[j])
        return PolyK(self.R, res)
    __rmul__ = __mul__

    # División euclídea: //, %, divmod (K cuerpo)
    def divmod(self, other):
        b = self._coerce_poly(other)
        if hasattr(self.R, "divmod"):
            return self.R.divmod(self, b)

        if b == 0:
            raise ZeroDivisionError("División por 0 en K[x]")
        r = PolyK(self.R, self.coef[:])
        db = b.grado
        q = [0] * max(0, self.grado - db + 1)
        inv_lead = b[db] ** -1
        while r.grado >= db:
            dr = r.grado
            pos = dr - db
            factor = r[dr] * inv_lead
            q[pos] = factor
            for j in range(db + 1):
                r.coef[pos + j] = r.coef[pos + j] - factor * b[j]
            r.coef = self._trim(r.coef)
        return PolyK(self.R, q), r

    def __floordiv__(self, other):
        q, _ = PolyK.divmod(self, other)
        return q

    def __mod__(self, other):
        _, r = PolyK.divmod(self, other)
        return r

    def __pow__(self, k: int):
        if k < 0:
            raise TypeError("Potencia negativa no definida en K[var]")
        if k == 0:
            return PolyK(self.R, 1)
        q, r = builtins.divmod(k, 2)
        if r == 0:
            return (self * self) ** q
        else:
            return self * (self * self) ** q

    # -------- utilidades algorítmicas (delegación si existe) --------
    @staticmethod
    def gcd(a: "PolyK", b: "PolyK") -> "PolyK":
        if a.R is not b.R:
            raise TypeError("Polinomios de anillos distintos")
        if hasattr(a.R, "gcd"):
            return a.R.gcd(a, b)
        # fallback genérico
        while b.grado != -1:
            a, b = b, (a % b)
        if a.grado == -1:
            return a
        return a * (a[a.grado] ** -1)

    @staticmethod
    def gcd_ext(a: "PolyK", b: "PolyK"):
        if a.R is not b.R:
            raise TypeError("Polinomios de anillos distintos")
        if hasattr(a.R, "gcd_ext"):
            return a.R.gcd_ext(a, b)
        # fallback genérico
        R = a.R
        x0, x1 = PolyK(R, 1), PolyK(R, 0)
        y0, y1 = PolyK(R, 0), PolyK(R, 1)
        r0, r1 = a, b
        while r1.grado != -1:
            q, r = PolyK.divmod(r0, r1)
            r0, r1 = r1, r
            x0, x1 = x1, x0 - q * x1
            y0, y1 = y1, y0 - q * y1
        inv = r0[r0.grado] ** -1
        return r0 * inv, x0 * inv, y0 * inv

    @staticmethod
    def inv_mod(a: "PolyK", m: "PolyK") -> "PolyK":
        if a.R is not m.R:
            raise TypeError("Polinomios de anillos distintos")
        if hasattr(a.R, "inv_mod"):
            return a.R.inv_mod(a, m)
        g, x, _ = PolyK.gcd_ext(a, m)
        if g.grado != 0 or not (g[0] == 1):
            raise ValueError("No existe inverso módulo m")
        return x % m

    @staticmethod
    def pot_mod(a: "PolyK", k: int, m: "PolyK") -> "PolyK":
        if a.R is not m.R:
            raise TypeError("Polinomios de anillos distintos")
        if hasattr(a.R, "pot_mod"):
            return a.R.pot_mod(a, int(k), m)
        # fallback genérico
        if k < 0:
            ainv = PolyK.inv_mod(a, m)
            return PolyK.pot_mod(ainv, -k, m)
        if k == 0:
            return PolyK(a.R, 1)
        q, r = builtins.divmod(k, 2)
        if r == 0:
            return PolyK.pot_mod((a * a) % m, q, m)
        else:
            return (a * PolyK.pot_mod((a * a) % m, q, m)) % m

    @property
    def es_irreducible(self) -> bool:
        # delegación rápida al anillo
        if hasattr(self.R, "es_irreducible"):
            return self.R.es_irreducible(self)
        # fallback (Rabin genérico)
        f = self
        R = f.R
        K = R.K
        n = f.grado
        q = getattr(K, "fq", None) or getattr(K, "p", None)
        if n <= 0:
            return False
        x = PolyK(R, (0, 1))
        h = x
        for _ in range(1, n // 2 + 1):
            h = PolyK.pot_mod(h, q, f)
            g = PolyK.gcd(h - x, f)
            if g != PolyK(R, 1):
                return False
        h = PolyK.pot_mod(x, q ** n, f)
        return (h - x) % f == 0

    @property
    def factorizar(self):
        # Si el anillo ya tiene kernel rápido, úsalo
        if hasattr(self.R, "factorizar"):
            return self.R.factorizar(self)

        # Fallback genérico SFD → DDF → EDF
        f = self
        L  = f._SDF1()          # [(h, s)]
        L2 = f._DDF1(L)         # [(h, s, d)]
        facs = f._EDF1(L2)      # [(irred, s)]

        # Limpieza: filtra constantes y moniciza
        cleaned = []
        for h, s in facs:
            if h.grado <= 0:
                continue
            lc = h[h.grado]
            h  = h * (lc ** -1)
            cleaned.append((h, s))
        return cleaned


    # ----------------- Fallback SFD/DDF/EDF genéricos -----------------
    def _SDF1(self):
        """Squarefree Decomposition (SFD) estándar."""
        f = self
        R = self.R
        K = R.K
        p = getattr(K, "p", None)
        L = []
        if f == 0:
            return L
        df = PolyK(R, [i * f[i] for i in range(1, len(f))])
        if df == 0:
            if p is None:
                raise ValueError("Característica desconocida para SFD")
            # f = g(x^p), tomar raíz p-ésima: índices múltiplos de p
            g = PolyK(R, [f[i] for i in range(0, len(f), p)])
            for (h, s) in g._SDF1():
                L.append((h, s * p))
            return L
        c = PolyK.gcd(f, df)
        w = f // c
        i = 1
        while w != 1:
            y = PolyK.gcd(w, c)
            z = w // y
            if z != 1:
                L.append((z, i))
            w = y
            c = c // y
            i += 1
        return L

    def _DDF1(self, L):
        """
        L: lista de (f, s) desde SFD.
        Devuelve lista de (h, s, d) donde:
        - h es libre de cuadrados,
        - todos los irreducibles de h tienen grado d,
        - s es la multiplicidad heredada desde SFD.
        """
        result = []
        R = self.R
        K = R.K
        q = getattr(K, "fq", None) or getattr(K, "p", None)
        x = PolyK(R, (0, 1))

        for f, s in L:
            i = 1
            g = f
            while i <= g.grado // 2 and g != PolyK(R, 1):
                h = PolyK.pot_mod(x, pow(q, i), g) - x
                d = PolyK.gcd(g, h)
                if d != PolyK(R, 1):
                    result.append((d, s, i))
                    g = g // d
                else:
                    i += 1
            if g != PolyK(R, 1):
                result.append((g, s, g.grado))
        return result


    def _EDF1(self, L):
        """
        L: lista de (f, s, d) desde DDF.
        Devuelve lista final de (factor_irreducible, s) conservando multiplicidades.
        """
        result = []
        R = self.R
        K = R.K
        q = getattr(K, "fq", None) or getattr(K, "p", None)

        for f, s, d in L:
            g = f
            while g.grado > d:
                found = False
                # limita intentos para evitar bucles infinitos en casos degenerados
                for _ in range(64):
                    coeffs = [K.aleatorio() for _ in range(g.grado)] + [K.cero()]
                    a = (PolyK(R, coeffs) % g)
                    if a == PolyK(R, 0):
                        continue
                    h = PolyK.pot_mod(a, (pow(q, d) - 1) // 2, g)
                    r = PolyK.gcd(g, h - PolyK(R, 1))
                    if r != PolyK(R, 1) and r != g:
                        result.append((r, s))
                        g = g // r
                        found = True
                        break
                if not found:
                    # si no separa, sal del bucle para no estancarte
                    break
            if g.grado > 0:
                result.append((g, s))
        return result


# -----------------------------
# Anillo Fp[x]
# -----------------------------

# -----------------------------
# Anillo Fp[x] (kernel rápido, entradas ya normalizadas)
# -----------------------------

class anillo_fp_x:
    """
    Representación interna (privada) de polinomios: lista de enteros en 0..p-1, sin ceros de cola,
    excepto el cero como [0]. Los métodos _*L asumen entradas YA normalizadas.
    """
    def __init__(self, fp, var='x'):
        self.K = fp if isinstance(fp, cuerpo_fp) else cuerpo_fp(int(fp))
        self.p = self.K.p
        self.var = var

    def __repr__(self):
        return f"F_{self.K.p}[{self.var}]"

    # ---------- Coerciones públicas (aquí sí normalizamos) ----------
    def _coerce_coef(self, c):
        return ElementoFp(self.K, c)  # garantiza 0..p-1

    def _to_int(self, c) -> int:
        # Para coeficientes sueltos, aquí sí normalizamos.
        if isinstance(c, ElementoFp):
            if c.P is not self.K:
                raise TypeError("ElementoFp de cuerpo distinto.")
            return int(c)
        return int(c) % self.p

    def _to_list(self, a) -> list[int]:
        """Convierte a lista de enteros 0..p-1 sin ceros de cola (frontera pública)."""
        if isinstance(a, PolyK) and a.R is self:
            L = [self._to_int(c) for c in a.coef]
        elif isinstance(a, (list, tuple)):
            L = [self._to_int(c) for c in a]
        else:
            L = [self._to_int(a)]
        return self._trim(L)

    def _from_list(self, L: list[int]) -> PolyK:
        """Crea un PolyK desde lista ya normalizada."""
        return PolyK(self, [self.K.elem(ci) for ci in L])

    def _trim(self, L: list[int]) -> list[int]:
        i = len(L) - 1
        while i > 0 and L[i] == 0:
            i -= 1
        return L[:i+1] if L else [0]

    def _gradoL(self, L: list[int]) -> int:
        return -1 if (len(L) == 1 and L[0] == 0) else len(L) - 1

    def _is_unoL(self, L: list[int]) -> bool:
        return len(L) == 1 and L[0] == 1

    # ---------- Núcleo rápido (entradas YA normalizadas) ----------
    def _sumaL(self, A: list[int], B: list[int]) -> list[int]:
        add = self.K._add
        n = max(len(A), len(B))
        R = [0] * n
        for i in range(n):
            ai = A[i] if i < len(A) else 0
            bi = B[i] if i < len(B) else 0
            R[i] = add(ai, bi)
        return self._trim(R)

    def _negL(self, A: list[int]) -> list[int]:
        neg = self.K._neg
        return self._trim([neg(ai) for ai in A])

    def _multL(self, A: list[int], B: list[int]) -> list[int]:
        add, mul = self.K._add, self.K._mul
        da, db = self._gradoL(A), self._gradoL(B)
        if da == -1 or db == -1:
            return [0]
        R = [0] * (len(A) + len(B) - 1)
        for i, ai in enumerate(A):
            if ai == 0:
                continue
            for j, bj in enumerate(B):
                R[i + j] = add(R[i + j], mul(ai, bj))
        return self._trim(R)

    def _divmodL(self, A: list[int], B: list[int]) -> tuple[list[int], list[int]]:
        """Pre: A y B normalizados; B != 0."""
        if self._gradoL(B) == -1:
            raise ZeroDivisionError("División por 0 en Fp[x]")
        db = self._gradoL(B)
        inv_lead = self.K._inv(B[db])
        sub, mul = self.K._sub, self.K._mul
        R = A[:]                               # trabajamos sobre copia
        Q = [0] * max(0, self._gradoL(A) - db + 1)
        while self._gradoL(R) >= db:
            dr = self._gradoL(R)
            pos = dr - db
            factor = mul(R[dr], inv_lead)
            Q[pos] = factor
            # R -= factor * x^pos * B
            for j in range(db + 1):
                R[pos + j] = sub(R[pos + j], mul(factor, B[j]))
            # trim manual barato
            while len(R) > 1 and R[-1] == 0:
                R.pop()
        if not Q: Q = [0]
        if not R: R = [0]
        return self._trim(Q), self._trim(R)

    def _divmodL_fixed(self, A, B, db, inv_lead):
        """Versión fija para potencias modulares (Pre: A,B normalizados)."""
        sub, mul = self.K._sub, self.K._mul
        R = A[:]
        Q = [0] * max(0, self._gradoL(R) - db + 1)
        while self._gradoL(R) >= db:
            dr = self._gradoL(R)
            pos = dr - db
            factor = mul(R[dr], inv_lead)
            Q[pos] = factor
            for j in range(db + 1):
                R[pos + j] = sub(R[pos + j], mul(factor, B[j]))
            while len(R) > 1 and R[-1] == 0:
                R.pop()
        if not Q: Q = [0]
        if not R: R = [0]
        return Q, R

    def _modL_fixed(self, A, B, db, inv_lead):
        _, R = self._divmodL_fixed(A, B, db, inv_lead)
        return R

    def _pot_modL(self, A: list[int], k: int, M: list[int]) -> list[int]:
        """Pre: A,M normalizados; k >= 0."""
        if k < 0:
            raise TypeError("pot_mod con exponente negativo")
        if self._gradoL(M) == -1:
            raise ZeroDivisionError("División por 0 en Fp[x]")
        db = self._gradoL(M)
        inv_lead = self.K._inv(M[db])
        base = self._modL_fixed(A, M, db, inv_lead)
        res  = [1]
        while k:
            if k & 1:
                res = self._modL_fixed(self._multL(res, base), M, db, inv_lead)
            base = self._modL_fixed(self._multL(base, base), M, db, inv_lead)
            k >>= 1
        return res

    def _gcdL(self, A: list[int], B: list[int]) -> list[int]:
        A, B = self._trim(A), self._trim(B)
        while self._gradoL(B) != -1:
            _, R = self._divmodL(A, B)
            A, B = B, R
        d = self._gradoL(A)
        if d == -1:
            return [0]
        inv = self.K._inv(A[d])
        mul = self.K._mul
        return self._trim([mul(ci, inv) for ci in A])

    # ---------- API pública  ------------ #
    def cero(self) -> PolyK:
        return PolyK(self, [])
    
    def uno(self) -> PolyK:
        return PolyK(self, [1])
    
    def elem_de_tuple(self, a: Tuple[Union[int, ElementoFp], ...]) -> PolyK:
        return PolyK(self, list(a))
    
    def elem_de_int(self, a: int) -> PolyK:
        if a == 0:
            return self.cero()
        p = self.p
        coef, x = [], int(a)
        while x:
            coef.append(self.K.elem(x % p))
            x //= p
        return PolyK(self, coef)
    
    def elem_de_str(self, s: str) -> PolyK:
        s = s.replace(" ", "")
        if s in ("", "0"):
            return self.cero()
        if s[0] not in "+-":
            s = "+" + s
        term_pat = re.compile(r"([+-])(\d*)(%s?)(?:\^(\d+))?" % re.escape(self.var))
        pos = 0
        coef: List[ElementoFp] = []
        while pos < len(s):
            m = term_pat.match(s, pos)
            if not m or m.end() == pos:
                raise ValueError(f"Formato inválido cerca de: '{s[pos:]}'")
            sign, num, xflag, exp = m.groups()
            pos = m.end()
            if xflag:
                c_int = 1 if num == "" else int(num)
                d = 1 if exp is None else int(exp)
            else:
                if num == "":
                    continue
                c_int = int(num)
                d = 0
            if sign == "-":
                c_int = -c_int
            c = self._coerce_coef(c_int)
            while len(coef) <= d:
                coef.append(self.K.cero())
            coef[d] = coef[d] + c
        return PolyK(self, coef)

    def conv_a_tuple(self, a: PolyK) -> Tuple[ElementoFp, ...]:
        return tuple(a.coef)

    def conv_a_int(self, a: PolyK) -> int:
        res = 0
        for i, c in enumerate(a.coef):
            res += int(c) * (self.p ** i)
        return res

    def conv_a_str(self, a: PolyK) -> str:
        terms = []
        for i, c in enumerate(a.coef):
            if int(c) == 0:
                continue
            cs = f"({c})"
            if i == 0:
                terms.append(cs)
            elif i == 1:
                terms.append(f"{'' if cs=='(1)' else cs}{self.var}")
            else:
                terms.append(f"{'' if cs=='(1)' else cs}{self.var}^{i}")
        return "+".join(terms) if terms else "0"

    # ----- aritmética ----- #
    def suma(self, a, b) -> PolyK:
        A, B = self._to_list(a), self._to_list(b)
        return self._from_list(self._sumaL(A, B))

    def inv_adit(self, a) -> PolyK:
        A = self._to_list(a)
        return self._from_list(self._negL(A))

    def mult(self, a, b) -> PolyK:
        A, B = self._to_list(a), self._to_list(b)
        return self._from_list(self._multL(A, B))

    def mult_por_escalar(self, a, e) -> PolyK:
        A = self._to_list(a)
        k = self._to_int(e)
        mul = self.K._mul
        return self._from_list([mul(ai, k) for ai in A])

    def divmod(self, a, b):
        A, B = self._to_list(a), self._to_list(b)
        qL, rL = self._divmodL(A, B)
        return self._from_list(qL), self._from_list(rL)

    def div(self, a, b) -> PolyK:
        q, _ = self.divmod(a, b)
        return q

    def mod(self, a, b) -> PolyK:
        _, r = self.divmod(a, b)
        return r

    def grado(self, a):
        return self._gradoL(self._to_list(a))

    def gcd(self, a, b) -> PolyK:
        A, B = self._to_list(a), self._to_list(b)
        return self._from_list(self._gcdL(A, B))

    def gcd_ext(self, a, b):
        A, B = self._to_list(a), self._to_list(b)
        add, mul, inv = self.K._add, self.K._mul, self.K._inv
        X0, X1 = [1], [0]
        Y0, Y1 = [0], [1]
        R0, R1 = A[:], B[:]
        while self._gradoL(R1) != -1:
            qL, rL = self._divmodL(R0, R1)
            R0, R1 = R1, rL
            # X0, X1 := X1, X0 - q*X1
            t = self._sumaL(X0, self._negL(self._multL(qL, X1)))
            X0, X1 = X1, t
            # Y0, Y1 := Y1, Y0 - q*Y1
            t = self._sumaL(Y0, self._negL(self._multL(qL, Y1)))
            Y0, Y1 = Y1, t
        d = self._gradoL(R0)
        if d == -1:
            z = self._from_list([0])
            return z, z, z
        inv_lead = inv(R0[d])
        g = self._trim([mul(ci, inv_lead) for ci in R0])
        x = self._trim([mul(ci, inv_lead) for ci in X0])
        y = self._trim([mul(ci, inv_lead) for ci in Y0])
        return self._from_list(g), self._from_list(x), self._from_list(y)

    def inv_mod(self, a, m):
        g, x, _ = self.gcd_ext(a, m)
        if self.grado(g) != 0 or int(self.conv_a_tuple(g)[0]) != 1:
            raise ValueError("No existe inverso módulo m")
        return self.mod(x, m)

    def pot_mod(self, a, k: int, m):
        A, M = self._to_list(a), self._to_list(m)
        return self._from_list(self._pot_modL(A, int(k), M))

    def es_cero(self, a) -> bool:
        return self._gradoL(self._to_list(a)) == -1

    def es_uno(self, a) -> bool:
        return self._is_unoL(self._to_list(a))

    def es_igual(self, a, b) -> bool:
        return self._trim(self._to_list(a)) == self._trim(self._to_list(b))

    # ----- Test de irreducibilidad (Rabin) -----
    def es_irreducible(self, f):
        F = self._to_list(f)           # normaliza una vez
        n = self._gradoL(F)
        if n <= 0:
            return False
        p = self.p
        xL = [0, 1]                    # x
        H  = self._divmodL(xL, F)[1]   # x mod f
        for _ in range(1, n // 2 + 1):
            H = self._pot_modL(H, p, F)                      # H := H^p (mod f)
            g = self._gcdL(F, self._sumaL(H, self._negL(xL)))# gcd(f, H - x)
            if self._gradoL(g) != 0:                         # != 1
                return False
        Hn = self._pot_modL(xL, pow(p, n), F)                # x^{p^n} (mod f)
        return self._gradoL(self._divmodL(self._sumaL(Hn, self._negL(xL)), F)[1]) == -1


# ===================================================
# Fq = Fp[x]/(g) + núcleo rígido cuerpo_fq
# ===================================================

class cuerpo_fq:
    """
    Núcleo rígido para Fq = Fp[x]/(g), con p primo y g mónico irreducible de grado n.
    * Formato ideal de un elemento: lista de n enteros normalizados (coef. en base x^i).
    """
    def __init__(self, fq, g=1, var: str = 'a'):
        if isinstance(fq,int):
            p, n = factor_potencia(fq)
            self.p = p
            self.n = n
            self.fq = fq
            self.R = anillo_fp_x(p, var='x')
            g_poly = PolyK(self.R, g)
            g_poly = g_poly * (g_poly[g_poly.grado] ** -1)  # mónico
            if not g_poly.es_irreducible or g_poly.grado != n:
                while True:
                    coef = [self.R.K.elem(random.randrange(p)) for _ in range(n)]
                    g_try = PolyK(self.R, coef + [self.R.K.uno()])
                    if g_try.es_irreducible:
                        g_poly = g_try
                        break
        elif isinstance(fq,cuerpo_fp):
            if var=='a':
                a='x'
            else:
                a='a'
            self.R=anillo_fp_x(fq, var=a)
            g_poly = PolyK(self.R, g)
            g_poly = g_poly * (g_poly[g_poly.grado] ** -1)  # mónico
            if not g_poly.es_irreducible:
                raise ValueError(f"g no es irreducible. No se puede construir Fq")
            self.p=fq.p
            self.n=g_poly.grado
            self.fq=self.p**self.n
        self.var = var  # símbolo α
        self.g = g_poly
        # coeficientes de g en formato ideal (enteros) para reducción
        self._g_c = [int(c) for c in g_poly.coef]  # len = n+1, monico en n

    # --- helpers formato ideal ---
    def _vec_from_poly(self, poly: PolyK) -> List[int]:
        v = [0] * self.n
        for i, c in enumerate(poly.coef[:self.n]):
            v[i] = int(c)
        return v

    def _poly_from_vec(self, v: List[int]) -> PolyK:
        coef = [self.R.K.elem(ci) for ci in v]
        return PolyK(self.R, coef)

    # primitivas sobre vectores (formato ideal)
    def _add_vec(self, a: List[int], b: List[int]) -> List[int]:
        p = self.p
        return [(ai + bi) % p for ai, bi in zip(a, b)]

    def _neg_vec(self, a: List[int]) -> List[int]:
        p = self.p
        return [(-ai) % p for ai in a]

    def _sub_vec(self, a: List[int], b: List[int]) -> List[int]:
        p = self.p
        return [(ai - bi) % p for ai, bi in zip(a, b)]

    def _mul_vec(self, a: List[int], b: List[int]) -> List[int]:
        p, n = self.p, self.n
        # convolución
        conv = [0] * (2 * n - 1)
        for i in range(n):
            ai = a[i]
            if ai == 0:
                continue
            for j in range(n):
                conv[i + j] = (conv[i + j] + ai * b[j]) % p
        # reducción módulo g (monico): g(x) = x^n + g0 + g1 x + ... + g_{n-1} x^{n-1}
        g = self._g_c
        for k in range(len(conv) - 1, n - 1, -1):
            c = conv[k]
            if c == 0:
                continue
            shift = k - n
            # conv[shift + t] -= c * g[t] para t=0..n-1
            for t in range(n):
                conv[shift + t] = (conv[shift + t] - c * g[t]) % p
            # conv[k] se anula implícitamente (x^n -> -g_{<n})
        return conv[:n]

    def _pow_vec(self, a: List[int], k: int) -> List[int]:
        if k < 0:
            a = self._pow_vec(a, self.fq - 2)
            k = -k
        # exponenciación binaria en el cuerpo
        n = self.n
        p0 = [0] * n
        res = [0] * n
        res[0] = 1  # 1 en Fq
        base = a[:]
        while k:
            if k & 1:
                res = self._mul_vec(res, base)
            base = self._mul_vec(base, base)
            k >>= 1
        return res

    # --- API opaca (ElementoFq) ---

    def cero(self) -> "ElementoFq":
        return self.elem_de_tuple((0,))
    
    def uno(self) -> "ElementoFq":
        return self.elem_de_tuple((1,))

    def elem_de_tuple(self, a: Tuple[int, ...]) -> "ElementoFq":
        v = [int(ai) % self.p for ai in (list(a) + [0] * self.n)[:self.n]]
        return ElementoFq(self, self._poly_from_vec(v))

    def elem_de_int(self, a: int) -> "ElementoFq":
        t = int(a)
        digits = []
        for _ in range(self.n):
            digits.append(t % self.p)
            t //= self.p
        return ElementoFq(self, self._poly_from_vec(digits))

    def elem_de_str(self, s: str) -> "ElementoFq":
        poly = self.R.elem_de_str(s)
        return ElementoFq(self, poly)

    def conv_a_tuple(self, a: "ElementoFq") -> Tuple[int, ...]:
        return tuple(int(ci) for ci in a)

    def conv_a_int(self, a: "ElementoFq") -> int:
        p = self.p
        acc = 0
        for i, c in enumerate(a.rep.coef):
            acc += int(c) * (p ** i)
        return acc

    def conv_a_str(self, a: "ElementoFq") -> str:
        terms = []
        for i, c in enumerate(a):
            ci = int(c)
            if ci == 0:
                continue
            cs = str(ci)
            if i == 0:
                terms.append(cs)
            elif i == 1:
                terms.append(f"{'' if cs=='1' else cs}{self.var}")
            else:
                terms.append(f"{'' if cs=='1' else cs}{self.var}^{i}")
        return "+".join(terms) if terms else "0"

    def suma(self, a, b) -> "ElementoFq":
        a, b = ElementoFq(self, a), ElementoFq(self, b)
        return a + b
    def inv_adit(self, a) -> "ElementoFq":
        return -ElementoFq(self, a)
    def mult(self, a, b) -> "ElementoFq":
        a, b = ElementoFq(self, a), ElementoFq(self, b)
        return a * b
    def pot(self, a, k: int) -> "ElementoFq":
        return ElementoFq(self, a) ** k
    def inv_mult(self, a) -> "ElementoFq":
        return ElementoFq(self, a) ** -1
    def es_cero(self, a) -> bool:
        return ElementoFq(self, a) == 0
    def es_uno(self, a) -> bool:
        return ElementoFq(self, a) == 1
    def es_igual(self, a, b) -> bool:
        return ElementoFq(self, a) == ElementoFq(self, b)
    
    def aleatorio(self) -> "ElementoFq":
        coef = [self.R.K.elem(random.randrange(self.p)) for _ in range(self.n)]
        return ElementoFq(self, PolyK(self.R, coef))

    def tabla_suma(self):
        return [[self.conv_a_int(self.elem_de_int(i) + self.elem_de_int(j)) for j in range(self.fq)]
                for i in range(self.fq)]
    def tabla_mult(self):
        return [[self.conv_a_int(self.elem_de_int(i) * self.elem_de_int(j)) for j in range(self.fq)]
                for i in range(self.fq)]
    def tabla_inv_adit(self):
        return [self.conv_a_int(-self.elem_de_int(i)) for i in range(self.fq)]
    def tabla_inv_mult(self):
        t = ["*"]
        for i in range(1, self.fq):
            t.append(self.conv_a_int(self.uno() / self.elem_de_int(i)))
        return t
    def cuadrado_latino(self, a: "ElementoFq"):
        a = a if isinstance(a, ElementoFq) else self.elem_de_int(int(a))
        if self.es_cero(a):
            raise ValueError("cuadrado latino no definido para a=0")
        return [[self.conv_a_int(a * self.elem_de_int(i) + self.elem_de_int(j)) for j in range(self.fq)] for i in range(self.fq)]


class ElementoFq:
    """Elemento opaco de Fq = Fp[x]/(g). Operadores sobrecargados.
    La aritmética delega al núcleo rígido cuerpo_fq sobre vectores de enteros.
    """
    __slots__ = ("Q", "rep")  # rep: PolyK sobre R reducido mod g

    def __init__(self, parent: cuerpo_fq, x):
        self.Q = parent
        R = parent.R
        p = parent.p
        # normaliza x a PolyK sobre R y reduce mod g
        if isinstance(x, ElementoFq):
            poly = x.rep if x.Q is parent else PolyK(R, x.rep.coef)
        elif isinstance(x, PolyK):
            poly = x if x.R is R else PolyK(R, x.coef)
        elif isinstance(x, (list, tuple)):
            poly = PolyK(R, list(x))
        elif isinstance(x, str):
            s = x.replace(" ", "")
            if R.var in s or "^" in s:
                poly = R.elem_de_str(s)
            else:
                poly = PolyK(R, [R._coerce_coef(s)])
        elif isinstance(x, ElementoFp):
            if x.P.p != p:
                raise TypeError("ElementoFp con primo distinto al de Fq")
            poly = PolyK(R, [x.v])
        else:
            poly = PolyK(R, [R._coerce_coef(int(x))])
        self.rep = poly % parent.g

    def _coerce(self, other):
        return ElementoFq(self.Q, other)

    def __iter__(self):
        return iter(self.rep.coef)
    def __len__(self):
        return len(self.rep.coef)
    def __getitem__(self, i):
        return self.rep.coef[i]

    def __repr__(self):
        return f"F_{self.Q.fq}({self.Q.conv_a_str(self)})"
    def __str__(self):
        return self.Q.conv_a_str(self)
    def __hash__(self):
        return hash((id(self.Q), tuple(int(c) for c in self.rep.coef)))

    def __eq__(self, other):
        o = self._coerce(other)
        return self.rep == o.rep

    # aritmética via núcleo rígido (vectores)
    def __neg__(self):
        v = self.Q._vec_from_poly(self.rep)
        w = self.Q._neg_vec(v)
        return ElementoFq(self.Q, self.Q._poly_from_vec(w))

    def __add__(self, other):
        o = self._coerce(other)
        v = self.Q._vec_from_poly(self.rep)
        w = self.Q._vec_from_poly(o.rep)
        r = self.Q._add_vec(v, w)
        return ElementoFq(self.Q, self.Q._poly_from_vec(r))
    __radd__ = __add__

    def __sub__(self, other):
        o = self._coerce(other)
        v = self.Q._vec_from_poly(self.rep)
        w = self.Q._vec_from_poly(o.rep)
        r = self.Q._sub_vec(v, w)
        return ElementoFq(self.Q, self.Q._poly_from_vec(r))

    def __rsub__(self, other):
        o = self._coerce(other)
        v = self.Q._vec_from_poly(o.rep)
        w = self.Q._vec_from_poly(self.rep)
        r = self.Q._sub_vec(v, w)
        return ElementoFq(self.Q, self.Q._poly_from_vec(r))

    def __mul__(self, other):
        o = self._coerce(other)
        v = self.Q._vec_from_poly(self.rep)
        w = self.Q._vec_from_poly(o.rep)
        r = self.Q._mul_vec(v, w)
        return ElementoFq(self.Q, self.Q._poly_from_vec(r))
    __rmul__ = __mul__

    def __pow__(self, k: int):
        v = self.Q._vec_from_poly(self.rep)
        r = self.Q._pow_vec(v, k)
        return ElementoFq(self.Q, self.Q._poly_from_vec(r))

    def __truediv__(self, other):
        o = self._coerce(other)
        v = self.Q._vec_from_poly(self.rep)
        w = self.Q._vec_from_poly(o.rep)
        invw = self.Q._pow_vec(w, self.Q.fq - 2)
        r = self.Q._mul_vec(v, invw)
        return ElementoFq(self.Q, self.Q._poly_from_vec(r))

    def __rtruediv__(self, other):
        o = self._coerce(other)
        return o / self

class anillo_fq_x:
    """Anillo Fq[var] con aritmética rígida (listas de vectores) y
    factorización (SFD/DDF/EDF) sin delegar a PolyK.
    """
    def __init__(self, fq, var: str = 'x'):
        """Construye Fq[var].
        - Si `fq` es un `cuerpo_fq`, lo usa.
        - Si `fq` es un entero q = p^n, construye internamente Fq eligiendo
          un polinomio irreducible mónico de grado n en Fp[x].
        """
        if isinstance(fq, cuerpo_fq):
            self.Q = fq
        else:
            self.Q = cuerpo_fq(int(fq), var='a')  # deja que cuerpo_fq elija g irreducible
        self.K = self.Q
        self.var = var

    def __repr__(self):
        return f"F_{self.Q.fq}[{self.var}]"

    # ---------- helpers de coeficiente (ElementoFq <-> vector int) ----------
    def _coerce_coef(self, c):
        return ElementoFq(self.Q, c)

    def _zero_vec(self):
        return [0] * self.Q.n

    def _one_vec(self):
        return [1] + [0]*(self.Q.n-1)

    def _is_zero_vec(self, v):
        return all((x % self.Q.p) == 0 for x in v)

    def _is_uno_vec(self, v):
        p = self.Q.p
        return (v[0] % p) == 1 and all((x % p) == 0 for x in v[1:])

    def _is_unoL(self, L):
        return len(L) == 1 and self._is_uno_vec(L[0])

    def _to_vec(self, c):
        p, n = self.Q.p, self.Q.n
        if isinstance(c, ElementoFq) and c.Q is self.Q:
            v = [int(ci) % p for ci in c.rep.coef]
            return v + [0]*(n - len(v)) if len(v) < n else v[:n]
        if isinstance(c, ElementoFp) and c.P is self.Q.K:
            return [int(c) % p] + [0]*(n-1)
        if isinstance(c, int):
            return [c % p] + [0]*(n-1)
        if isinstance(c, (list, tuple)) and all(isinstance(x, int) for x in c):
            v = [x % p for x in c]
            return v + [0]*(n - len(v)) if len(v) < n else v[:n]
        e = ElementoFq(self.Q, c)
        v = [int(ci) % p for ci in e.rep.coef]
        return v + [0]*(n - len(v)) if len(v) < n else v[:n]

    # ---------- conversión de polinomios ----------
    def _to_list(self, a):
        if isinstance(a, PolyK) and a.R is self:
            return [self._to_vec(c) for c in a.coef]
        if isinstance(a, (list, tuple)):
            return [self._to_vec(c) for c in a]
        return [self._to_vec(a)]

    def _from_list(self, L: list[list[int]]) -> PolyK:
        K0 = self.Q.cero()
        coefs = []
        for v in L:
            if self._is_zero_vec(v):
                coefs.append(K0)
            else:
                coefs.append(ElementoFq(self.Q, self.Q._poly_from_vec(v)))
        return PolyK(self, coefs)

    def _trim(self, L: list[list[int]]) -> list[list[int]]:
        while len(L) > 1 and self._is_zero_vec(L[-1]):
            L.pop()
        if not L:
            L = [self._zero_vec()]
        return L

    def _gradoL(self, L: list[list[int]]) -> int:
        i = len(L) - 1
        while i >= 0 and self._is_zero_vec(L[i]):
            i -= 1
        return i

    # ---------- helpers núcleo sobre listas ----------
    def _sumaL(self, A, B):
        n = max(len(A), len(B))
        R = []
        add = self.Q._add_vec
        zv = self._zero_vec()
        for i in range(n):
            av = A[i] if i < len(A) else zv
            bv = B[i] if i < len(B) else zv
            R.append(add(av, bv))
        return self._trim(R)

    def _negL(self, A):
        return self._trim([self.Q._neg_vec(v) for v in A])

    def _multL(self, A, B):
        if self._gradoL(A) == -1 or self._gradoL(B) == -1:
            return [self._zero_vec()]
        R = [self._zero_vec()[:] for _ in range(len(A)+len(B)-1)]
        add = self.Q._add_vec; mul = self.Q._mul_vec
        for i, ai in enumerate(A):
            if self._is_zero_vec(ai): continue
            for j, bj in enumerate(B):
                R[i+j] = add(R[i+j], mul(ai, bj))
        return self._trim(R)

    def _divmodL(self, A, B):
        B = self._trim([v[:] for v in B]); A = self._trim([v[:] for v in A])
        db = self._gradoL(B)
        if db == -1:
            raise ZeroDivisionError("División por 0 en Fq[x]")
        inv_lead = self.Q._pow_vec(B[db], self.Q.fq - 2)
        R = [v[:] for v in A]
        Q = [self._zero_vec()[:] for _ in range(max(0, self._gradoL(A) - db + 1))]
        sub = self.Q._sub_vec; mul = self.Q._mul_vec
        while self._gradoL(R) >= db:
            dr = self._gradoL(R)
            pos = dr - db
            factor = mul(R[dr], inv_lead)
            Q[pos] = factor
            for j in range(db+1):
                R[pos+j] = sub(R[pos+j], mul(factor, B[j]))
            R = self._trim(R)
        return self._trim(Q), self._trim(R)

    def _modL(self, A, B):
        _, R = self._divmodL(A, B)
        return R

    def _divmodL_fixed(self, A, B, db, inv_lead):
        R = [v[:] for v in self._trim([v[:] for v in A])]
        if db == -1:
            raise ZeroDivisionError("División por 0 en Fq[x]")
        Q = [self._zero_vec()[:] for _ in range(max(0, self._gradoL(R) - db + 1))]
        sub = self.Q._sub_vec; mul = self.Q._mul_vec
        while self._gradoL(R) >= db:
            dr = self._gradoL(R)
            pos = dr - db
            factor = mul(R[dr], inv_lead)
            Q[pos] = factor
            for j in range(db + 1):
                R[pos + j] = sub(R[pos + j], mul(factor, B[j]))
            while len(R) > 0 and self._is_zero_vec(R[-1]):
                R.pop()
        if not Q: Q = [self._zero_vec()]
        if not R: R = [self._zero_vec()]
        return Q, R

    def _modL_fixed(self, A, B, db, inv_lead):
        _, R = self._divmodL_fixed(A, B, db, inv_lead)
        return R

    def _pot_modL(self, A, k, M):
        if k < 0:
            raise TypeError("pot_mod con exponente negativo")
        A = self._trim(A); M = self._trim(M)
        db = self._gradoL(M)
        if db == -1:
            raise ZeroDivisionError("División por 0 en Fq[x]")
        inv_lead = self.Q._pow_vec(M[db], self.Q.fq - 2)
        base = self._modL_fixed(A, M, db, inv_lead)
        res  = [self._one_vec()]   # 1
        while k:
            if k & 1:
                res = self._modL_fixed(self._multL(res, base), M, db, inv_lead)
            base = self._modL_fixed(self._multL(base, base), M, db, inv_lead)
            k >>= 1
        return res

    def _gcdL(self, A: list[list[int]], B: list[list[int]]) -> list[list[int]]:
        A = self._trim(A); B = self._trim(B)
        while self._gradoL(B) != -1:
            _, R = self._divmodL(A, B)
            A, B = B, R
        da = self._gradoL(A)
        if da == -1:
            return [self._zero_vec()]
        inv = self.Q._pow_vec(A[da], self.Q.fq - 2)
        return [self.Q._mul_vec(ci, inv) for ci in A]

    # ---------- API del profe ----------
    def cero(self) -> PolyK:
        return PolyK(self, [])

    def uno(self) -> PolyK:
        return PolyK(self, [self.Q.uno()])

    def elem_de_tuple(self, a):
        return PolyK(self, list(a))

    def elem_de_int(self, a: int) -> PolyK:
        if a == 0:
            return self.cero()
        q = self.Q.fq
        coef, x = [], int(a)
        while x:
            coef.append(self.K.elem_de_int(x % q)) 
            x //= q
        return PolyK(self, coef)

    def elem_de_str(self, s: str):
        s = str(s).replace(" ", "")
        if s in ("", "0"):
            return self.cero()
        if s[0] not in "+-":
            s = "+" + s
        term_pat = re.compile(r"([+-])(.*?)(%s?)(?:\^(\d+))?" % re.escape(self.var))
        pos = 0
        coef = []
        while pos < len(s):
            m = term_pat.match(s, pos)
            if not m or m.end() == pos:
                raise ValueError(f"Formato inválido cerca de: '{s[pos:]}'")
            sign, num, xflag, exp = m.groups()
            pos = m.end()
            if xflag:
                cpart = num
                d = 1 if exp is None else int(exp)
                c = self.Q.uno() if cpart in ("", "1") else ElementoFq(self.Q, cpart)
            else:
                if num == "":
                    continue
                d = 0
                c = ElementoFq(self.Q, num)
            if sign == "-":
                c = -c
            while len(coef) <= d:
                coef.append(self.Q.cero())
            coef[d] = coef[d] + c
        return PolyK(self, coef)

    def conv_a_tuple(self, a):
        return tuple(a.coef)

    def conv_a_int(self, a) -> int:
        res = 0
        for i, c in enumerate(a.coef):
            res += self.Q.conv_a_int(c) * (self.Q.fq ** i)
        return res

    def conv_a_str(self, a) -> str:
        terms = []
        for i, c in enumerate(a.coef):
            if self.Q.es_cero(c):
                continue
            cs = f"({self.Q.conv_a_str(c)})"
            if i == 0:
                terms.append(cs)
            elif i == 1:
                terms.append(f"{'' if cs=='(1)' else cs}{self.var}")
            else:
                terms.append(f"{'' if cs=='(1)' else cs}{self.var}^{i}")
        return "+".join(terms) if terms else "0"

    # ----- aritmética rígida sobre listas de vectores -----
    def suma(self, a, b) -> PolyK:
        A, B = self._to_list(a), self._to_list(b)
        n = max(len(A), len(B))
        R = []
        add = self.Q._add_vec
        zv = self._zero_vec()
        for i in range(n):
            av = A[i] if i < len(A) else zv
            bv = B[i] if i < len(B) else zv
            R.append(add(av, bv))
        return self._from_list(self._trim(R))

    def inv_adit(self, a) -> PolyK:
        A = self._to_list(a)
        R = [self.Q._neg_vec(v) for v in A]
        return self._from_list(self._trim(R))

    def mult(self, a, b) -> PolyK:
        A, B = self._to_list(a), self._to_list(b)
        if self._gradoL(A) == -1 or self._gradoL(B) == -1:
            return self._from_list([self._zero_vec()])
        R = [self._zero_vec()[:] for _ in range(len(A) + len(B) - 1)]
        add = self.Q._add_vec; mul = self.Q._mul_vec
        for i, ai in enumerate(A):
            if self._is_zero_vec(ai):
                continue
            for j, bj in enumerate(B):
                prod = mul(ai, bj)
                R[i + j] = add(R[i + j], prod)
        return self._from_list(self._trim(R))

    def mult_por_escalar(self, a, e) -> PolyK:
        A = self._to_list(a)
        ev = self._to_vec(e)
        R = [self.Q._mul_vec(v, ev) for v in A]
        return self._from_list(self._trim(R))

    def divmod(self, a, b):
        A, B = self._to_list(a), self._to_list(b)
        A, B = self._trim([v[:] for v in A]), self._trim([v[:] for v in B])
        db = self._gradoL(B)
        if db == -1:
            raise ZeroDivisionError("División por 0 en Fq[x]")
        R = [v[:] for v in A]
        Q = [self._zero_vec()[:] for _ in range(max(0, self._gradoL(A) - db + 1))]
        inv_lead = self.Q._pow_vec(B[db], self.Q.fq - 2)
        sub = self.Q._sub_vec; mul = self.Q._mul_vec
        while self._gradoL(R) >= db:
            dr = self._gradoL(R)
            pos = dr - db
            factor = mul(R[dr], inv_lead)
            Q[pos] = factor
            for j in range(db + 1):
                R[pos + j] = sub(R[pos + j], mul(factor, B[j]))
            R = self._trim(R)
        return self._from_list(self._trim(Q)), self._from_list(self._trim(R))

    def div(self, a, b):
        q, _ = self.divmod(a, b)
        return q

    def mod(self, a, b):
        _, r = self.divmod(a, b)
        return r

    def grado(self, a):
        return self._gradoL(self._to_list(a))

    def gcd(self, a, b):
        A, B = self._trim(self._to_list(a)), self._trim(self._to_list(b))
        while self._gradoL(B) != -1:
            _, r = self.divmod(self._from_list(A), self._from_list(B))
            A, B = B, self._trim(self._to_list(r))
        da = self._gradoL(A)
        if da == -1:
            return self._from_list([self._zero_vec()])
        inv = self.Q._pow_vec(A[da], self.Q.fq - 2)
        G = [self.Q._mul_vec(ci, inv) for ci in A]
        return self._from_list(self._trim(G))

    def gcd_ext(self, a, b):
        A, B = self._trim(self._to_list(a)), self._trim(self._to_list(b))
        X0, X1 = [self._one_vec()], [self._zero_vec()]
        Y0, Y1 = [self._zero_vec()], [self._one_vec()]
        R0, R1 = A[:], B[:]
        while self._gradoL(R1) != -1:
            q, r = self.divmod(self._from_list(R0), self._from_list(R1))
            R0, R1 = R1, self._trim(self._to_list(r))
            X0, X1 = X1, self._to_list(
                self.suma(self._from_list(X0),
                          self.inv_adit(self.mult(q, self._from_list(X1))))
            )
            Y0, Y1 = Y1, self._to_list(
                self.suma(self._from_list(Y0),
                          self.inv_adit(self.mult(q, self._from_list(Y1))))
            )
        inv = self.Q._pow_vec(R0[self._gradoL(R0)], self.Q.fq - 2)
        g = [self.Q._mul_vec(ci, inv) for ci in R0]
        x = [self.Q._mul_vec(ci, inv) for ci in X0]
        y = [self.Q._mul_vec(ci, inv) for ci in Y0]
        return self._from_list(self._trim(g)), self._from_list(self._trim(x)), self._from_list(self._trim(y))

    def inv_mod(self, a, b):
        g, x, _ = self.gcd_ext(a, b)
        if self.grado(g) != 0 or not self.Q.es_uno(self.conv_a_tuple(g)[0]):
            raise ValueError("No existe inverso módulo b")
        return self.mod(x, b)

    def pot_mod(self, a, k: int, b):
        if int(k) < 0:
            raise TypeError("pot_mod con exponente negativo")
        A, M = self._to_list(a), self._to_list(b)
        R = self._pot_modL(A, int(k), M)
        return self._from_list(self._trim(R))

    def es_cero(self, a):
        return self.grado(a) == -1

    def es_uno(self, a):
        L = self._to_list(a)
        return len(L) == 1 and self.Q.es_uno(ElementoFq(self.Q, self.Q._poly_from_vec(L[0])))

    def es_igual(self, a, b):
        return self._trim(self._to_list(a)) == self._trim(self._to_list(b))

    # ---------- Irreducibilidad (Rabin) ----------
    def es_irreducible(self, f):
        F = self._to_list(f)
        n = self._gradoL(F)
        if n <= 0:
            return False
        q = self.Q.fq
        xL = [self._zero_vec(), self._one_vec()]
        H  = self._modL(xL, F)
        for _ in range(1, n//2 + 1):
            H = self._pot_modL(H, q, F)
            g = self._gcdL(self._sumaL(H, self._negL(xL)), F)
            if self._gradoL(g) != 0:
                return False
        Hn = self._pot_modL(xL, pow(q, n), F)
        return self._gradoL(self._modL(self._sumaL(Hn, self._negL(xL)), F)) == -1
    

    