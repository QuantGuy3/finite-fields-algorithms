"""Anillos de polinomios K[x], polinomios y factorizaciones.

* :class:`PolynomialRing` (también ``K["x"]`` o ``K.poly_ring("x")``): anillo de
  polinomios en una variable sobre un cuerpo finito; cacheado por (K, variable).
* :class:`Polynomial`: polinomio inmutable (envuelve un PolyFp/PolyFq nativo).
* :class:`Factorization`: resultado de ``f.factor()``.

Conversión en un anillo (``R(x)``, y en las operaciones con un polinomio):

=========================================  ===================================
``int``, ``Fraction``, elemento del cuerpo  polinomio constante
``str``                                     se analiza: ``"x^3 + (a+1)x - 2"``
``list`` / ``tuple``                        coeficientes en orden creciente
                                            (cada uno convertido en el cuerpo)
``dict {exponente: coeficiente}``           polinomio disperso
polinomio del mismo anillo                  él mismo
polinomio sobre el subcuerpo primo          se sumerge (misma variable)
otro iterable (solo en ``R(x)``)            como una lista
``float``                                   ``TypeError``
=========================================  ===================================

Un polinomio sobre F_p operado con un polinomio o un elemento de una extensión
de característica p da un polinomio sobre la extensión.  Anillos con cuerpos o
variables distintas no se mezclan en las operaciones (``ValueError``); la
conversión explícita ``R(f)`` sí cambia de variable.
"""
from __future__ import annotations

import operator

from . import _parse, _pretty
from .field import _PROMOVER, FieldElement, FiniteField, _texto_elemento, validar_nombre
from ._util import entero

#: Grado máximo de un polinomio construido desde Python (protege la memoria).
GRADO_MAXIMO = 10**7

_OPERADORES = {"add": operator.add, "sub": operator.sub, "mul": operator.mul,
               "floordiv": operator.floordiv, "mod": operator.mod,
               "divmod": divmod, "truediv": operator.truediv}


def _poli_desde_indices(R: "PolynomialRing", coefs: list) -> "Polynomial":
    """Reconstruye un polinomio a partir de los valores/índices de sus coeficientes (pickle)."""
    return R._envolver(R._K._nat.poly(coefs))


def _comprobar_grado(d: int) -> None:
    if d > GRADO_MAXIMO:
        raise ValueError(f"grado {d} demasiado grande (máximo {GRADO_MAXIMO})")


# ---------------------------------------------------------------------------
# Anillo
# ---------------------------------------------------------------------------
class PolynomialRing(_pretty.Imprimible):
    """Anillo de polinomios K[var] sobre un cuerpo finito (cacheado).

    ``PolynomialRing(GF(7), "x") is GF(7)["x"]``.
    """

    __slots__ = ("_K", "_var")

    def __new__(cls, field, var: str = "x"):
        if not isinstance(field, FiniteField):
            raise TypeError(f"el cuerpo base debe ser un cuerpo finito (alcp.GF), "
                            f"no {type(field).__name__!r}")
        var = validar_nombre(var)
        R = field._anillos.get(var)
        if R is None:
            if var == field._var:
                raise ValueError(f"la variable {var!r} ya es el generador de {field!r}")
            R = object.__new__(cls)
            R._K = field
            R._var = var
            field._anillos[var] = R
        return R

    # -- datos ------------------------------------------------------------
    @property
    def field(self) -> FiniteField:
        """Cuerpo de coeficientes."""
        return self._K

    base_field = base_ring = field

    @property
    def var(self) -> str:
        """Nombre de la variable."""
        return self._var

    @property
    def gen(self) -> "Polynomial":
        """La variable como polinomio (también ``R.x`` si la variable es x)."""
        return self._envolver(self._K._nat.poly([0, 1]))

    def __getattr__(self, nombre: str):
        # R.x (o R.t...) es el generador si así se llama la variable.
        if not nombre.startswith("_") and nombre == self._var:
            return self.gen
        raise AttributeError(f"{type(self).__name__!r} no tiene el atributo {nombre!r}")

    @property
    def zero(self) -> "Polynomial":
        return self._envolver(self._K._nat.poly([]))

    @property
    def one(self) -> "Polynomial":
        return self._envolver(self._K._nat.poly([1]))

    # -- construcción ----------------------------------------------------
    def _envolver(self, nat) -> "Polynomial":
        f = object.__new__(Polynomial)
        f._R = self
        f._p = nat
        f._h = None
        return f

    def __call__(self, x=0) -> "Polynomial":
        """Convierte ``x`` en un polinomio de este anillo (ver la tabla del módulo)."""
        if x.__class__ is Polynomial and x._R is self:
            return x
        return self._envolver(self._nativo(x))

    def __contains__(self, x) -> bool:
        """``x in R``: ¿se puede convertir ``x`` en un polinomio de R?"""
        try:
            self._nativo(x)
        except (TypeError, ValueError, ArithmeticError):
            return False
        return True

    def monomial(self, deg, c=1) -> "Polynomial":
        """``c·var^deg`` (deg >= 0)."""
        deg = entero(deg, "el grado")
        if deg < 0:
            raise ValueError(f"el grado de un monomio debe ser >= 0, no {deg}")
        _comprobar_grado(deg)
        return self._envolver(self._K._nat.poly([0] * deg + [self._K._valor(c)]))

    def random(self, deg) -> "Polynomial":
        """Polinomio aleatorio de grado exactamente ``deg`` (el cero si deg < 0)."""
        deg = entero(deg, "el grado")
        _comprobar_grado(deg)
        return self._envolver(self._K._nat.random_poly(deg))

    def irreducible(self, deg) -> "Polynomial":
        """Polinomio mónico irreducible aleatorio de grado ``deg`` >= 1 (ver set_seed)."""
        deg = entero(deg, "el grado")
        if deg < 1:
            raise ValueError(f"un polinomio irreducible tiene grado >= 1, no {deg}")
        _comprobar_grado(deg)
        while True:
            f = self._K._nat.random_poly(deg).monic()
            if f.is_irreducible():
                return self._envolver(f)

    # -- conversión --------------------------------------------------------
    def _nativo(self, x, operando: bool = False):
        """Polinomio nativo de ``x`` en este anillo.

        Con ``operando=True`` devuelve ``None`` si el tipo no se admite como
        operando y :data:`_PROMOVER` si ``x`` vive sobre una extensión de este
        cuerpo primo (la operación se hará allí).
        """
        K = self._K
        if isinstance(x, Polynomial):
            return self._desde_polinomio(x, operando)
        if isinstance(x, str):
            return self._analizar(x)._p
        if isinstance(x, dict):
            return self._desde_dict(x)
        if isinstance(x, (list, tuple)):
            return K._nat.poly([K._valor(c) for c in x])
        if isinstance(x, FieldElement):
            v = K._valor(x, operando)
            return v if v is _PROMOVER else K._nat.poly([v])
        v = K._valor(x, True)        # enteros, fracciones...; None si el tipo no vale
        if v is not None:
            return K._nat.poly([v])
        if operando:
            return None
        try:
            coefs = list(x)
        except TypeError:
            raise TypeError(f"no se puede convertir un valor de tipo {type(x).__name__!r} "
                            f"en un polinomio de {self!r}") from None
        return K._nat.poly([K._valor(c) for c in coefs])

    def _desde_polinomio(self, f: "Polynomial", operando: bool):
        S = f._R
        if S is self:
            return f._p
        K, L = self._K, S._K
        if not operando:
            # Conversión explícita: coeficiente a coeficiente (admite otra variable).
            return K._nat.poly([K._valor(c) for c in f])
        if S._var == self._var and L._p == K._p:
            if K._n == 1:
                return _PROMOVER
            if L._n == 1:
                return K._nat.poly(f._p.coeffs())   # F_p se sumerge con el mismo índice
        raise ValueError(f"{S!r} y {self!r} son anillos distintos")

    def _desde_dict(self, d: dict):
        K = self._K
        terminos = {}
        for k, c in d.items():
            k = entero(k, "el exponente")
            if k < 0:
                raise ValueError(f"exponente negativo en un polinomio: {k}")
            _comprobar_grado(k)
            terminos[k] = K._valor(c)
        coefs = [0] * (max(terminos, default=-1) + 1)
        for k, v in terminos.items():
            coefs[k] = v
        return K._nat.poly(coefs)

    def _analizar(self, texto: str, usadas=None) -> "Polynomial":
        variables = {self._var: self.gen}
        if self._K._var:
            variables[self._K._var] = self._envolver(self._K._nat.poly([self._K.gen._v]))
        return _parse.analizar(texto, self._constante, variables, _dividir_constante, usadas)

    def _constante(self, n: int) -> "Polynomial":
        return self._envolver(self._K._nat.poly([n % self._K._p]))

    # -- igualdad, pickle, representación --------------------------------------
    def __eq__(self, otro):
        if self is otro:
            return True
        if not isinstance(otro, PolynomialRing):
            return NotImplemented
        return self._K == otro._K and self._var == otro._var

    def __hash__(self) -> int:
        return hash((self._K, self._var))

    def __reduce__(self):
        return (PolynomialRing, (self._K, self._var))

    def __copy__(self):
        return self

    def __deepcopy__(self, memo):
        return self

    def __repr__(self) -> str:
        return f"{self._K!r}[{self._var}]"

    def _texto(self, estilo: str) -> str:
        return (f"{_pretty.nombre_cuerpo(self._K._q, estilo)}"
                f"[{_pretty.variable(self._var, estilo)}]")


def _dividir_constante(a: "Polynomial", b: "Polynomial") -> "Polynomial":
    """División del analizador: solo por constantes."""
    if b.degree > 0:
        raise ValueError("solo se puede dividir por constantes")
    return a / b


# ---------------------------------------------------------------------------
# Polinomio
# ---------------------------------------------------------------------------
class Polynomial(_pretty.Imprimible):
    """Polinomio inmutable sobre un cuerpo finito.

    Se construye con el anillo: ``R = GF(7)["x"]; f = R("x^2 + 1")`` (o
    ``Polynomial(R, datos)``).  Admite ``+ - * // % divmod **`` (y
    ``pow(f, k, m)``), ``/`` por constantes o divisores exactos, ``==`` con
    cualquier valor convertible, ``f[i]``, ``len(f)`` (grado + 1), iteración por
    los coeficientes y ``f(v)`` (evaluación o composición).
    """

    __slots__ = ("_R", "_p", "_h")
    __array_ufunc__ = None

    def __new__(cls, ring, data=0):
        if not isinstance(ring, PolynomialRing):
            raise TypeError("Polynomial(ring, data): ring debe ser un PolynomialRing")
        return ring(data)

    # -- datos ------------------------------------------------------------
    @property
    def ring(self) -> PolynomialRing:
        return self._R

    parent = ring

    @property
    def field(self) -> FiniteField:
        """Cuerpo de coeficientes."""
        return self._R._K

    base_field = field

    @property
    def var(self) -> str:
        return self._R._var

    @property
    def degree(self) -> int:
        """Grado (-1 para el polinomio cero)."""
        return self._p.degree()

    @property
    def coeffs(self) -> list:
        """Coeficientes (elementos del cuerpo) en orden creciente; [] para el cero."""
        crear = self._R._K._crear
        return [crear(c) for c in self._p.coeffs()]

    def coefficients(self, sparse: bool = False) -> list:
        """Coeficientes en orden creciente; con ``sparse=True``, solo los no nulos."""
        cs = self.coeffs
        return [c for c in cs if c] if sparse else cs

    def to_dict(self) -> dict:
        """``{exponente: coeficiente}`` con los coeficientes no nulos."""
        return {k: c for k, c in enumerate(self.coeffs) if c}

    @property
    def leading_coefficient(self) -> FieldElement:
        """Coeficiente principal (el 0 para el polinomio cero)."""
        p = self._p
        return self._R._K._crear(0 if p.is_zero() else p.leading())

    lc = leading_coefficient

    def is_zero(self) -> bool:
        return self._p.is_zero()

    def is_one(self) -> bool:
        return self._p.is_one()

    def is_monic(self) -> bool:
        """¿Es mónico?  (El polinomio cero no lo es.)"""
        return not self._p.is_zero() and self._p.leading() == 1

    def is_constant(self) -> bool:
        """¿Es constante (grado <= 0)?"""
        return self._p.degree() <= 0

    # -- secuencia ------------------------------------------------------------
    def __len__(self) -> int:
        """Grado + 1 (0 para el polinomio cero)."""
        return self._p.degree() + 1

    def __getitem__(self, i):
        """Coeficiente de x^i (0 si i > grado); índices negativos y cortes como en una lista."""
        if isinstance(i, slice):
            return self.coeffs[i]
        i = entero(i, "el índice")
        if i < 0:
            i += len(self)
            if i < 0:
                raise IndexError("índice de coeficiente fuera de rango")
        K = self._R._K
        return K._crear(self._p.coeff(i) if i <= self._p.degree() else 0)

    def __iter__(self):
        return iter(self.coeffs)

    def __bool__(self) -> bool:
        return not self._p.is_zero()

    # -- operaciones -------------------------------------------------------
    def _par(self, otro):
        """``(anillo, nativo de self, nativo de otro)`` para operar con ``otro``.

        Devuelve None si ``otro`` no es de un tipo admitido como operando.
        """
        R = self._R
        b = R._nativo(otro, True)
        if b is None:
            return None
        if b is _PROMOVER:
            # self está sobre F_p y otro sobre una extensión L: se opera en L[x].
            L = otro._R._K if isinstance(otro, Polynomial) else otro._K
            S = L[R._var]
            return S, S._nativo(self), S._nativo(otro)
        return R, self._p, b

    def _binaria(self, otro, nombre: str, reflejada: bool):
        par = self._par(otro)
        if par is None:
            return NotImplemented
        S, a, b = par
        if reflejada:
            a, b = b, a
        if nombre in ("floordiv", "mod", "divmod", "truediv") and b.is_zero():
            raise ZeroDivisionError("división por el polinomio cero")
        if nombre == "truediv":
            return S._envolver(_division_exacta(S, a, b))
        r = _OPERADORES[nombre](a, b)
        if nombre == "divmod":
            return S._envolver(r[0]), S._envolver(r[1])
        return S._envolver(r)

    def __add__(self, otro):
        if otro.__class__ is Polynomial and otro._R is self._R:
            return self._R._envolver(self._p + otro._p)
        return self._binaria(otro, "add", False)

    def __radd__(self, otro):
        return self._binaria(otro, "add", True)

    def __sub__(self, otro):
        if otro.__class__ is Polynomial and otro._R is self._R:
            return self._R._envolver(self._p - otro._p)
        return self._binaria(otro, "sub", False)

    def __rsub__(self, otro):
        return self._binaria(otro, "sub", True)

    def __mul__(self, otro):
        if otro.__class__ is Polynomial and otro._R is self._R:
            return self._R._envolver(self._p * otro._p)
        return self._binaria(otro, "mul", False)

    def __rmul__(self, otro):
        return self._binaria(otro, "mul", True)

    def __floordiv__(self, otro):
        return self._binaria(otro, "floordiv", False)

    def __rfloordiv__(self, otro):
        return self._binaria(otro, "floordiv", True)

    def __mod__(self, otro):
        return self._binaria(otro, "mod", False)

    def __rmod__(self, otro):
        return self._binaria(otro, "mod", True)

    def __divmod__(self, otro):
        return self._binaria(otro, "divmod", False)

    def __rdivmod__(self, otro):
        return self._binaria(otro, "divmod", True)

    def __truediv__(self, otro):
        """``f / g``: g constante no nula, o divisor exacto de f (si no, ValueError)."""
        return self._binaria(otro, "truediv", False)

    def __rtruediv__(self, otro):
        return self._binaria(otro, "truediv", True)

    def __neg__(self):
        return self._R._envolver(-self._p)

    def __pos__(self):
        return self

    def __pow__(self, k, modulo=None):
        """``f**k`` (k >= 0; k < 0 solo para constantes no nulas) o ``pow(f, k, m)``."""
        if modulo is not None:
            return self.pow_mod(k, modulo)
        if isinstance(k, FieldElement):
            raise TypeError("el exponente debe ser un entero, no un elemento de un cuerpo")
        k = entero(k, "el exponente")
        R, d = self._R, self._p.degree()
        if k < 0:
            if d < 0:
                raise ZeroDivisionError("el polinomio cero no tiene inverso")
            if d > 0:
                raise ValueError("un polinomio no constante no tiene inverso: "
                                 "usa pow(f, k, m) para potencias modulares")
            return R._envolver(R._K._nat.poly([R._K._nat.pow(self._p.coeff(0), k)]))
        if d > 0:
            _comprobar_grado(d * k)
        return R._envolver(self._p ** k)

    def _otro(self, otro):
        """Como :meth:`_par`, pero para los métodos: TypeError si no es convertible."""
        par = self._par(otro)
        if par is None:
            raise TypeError(f"no se puede convertir un valor de tipo {type(otro).__name__!r} "
                            f"en un polinomio de {self._R!r}")
        return par

    def pow_mod(self, k, m) -> "Polynomial":
        """``f**k mod m`` (k < 0: potencia del inverso de f módulo m)."""
        k = entero(k, "el exponente")
        S, a, b = self._otro(m)
        if b.is_zero():
            raise ZeroDivisionError("potencia modular con módulo cero")
        try:
            return S._envolver(a.pow_mod(k, b))
        except ValueError:
            raise ValueError(f"{self:ascii} no es invertible módulo {S._envolver(b):ascii}") from None

    def mul(self, otro, algorithm: str = "auto") -> "Polynomial":
        """Producto con el algoritmo indicado: "auto", "school", "karatsuba" o "fft"."""
        if not isinstance(algorithm, str):
            raise TypeError("algorithm debe ser una cadena: 'auto', 'school', 'karatsuba' o 'fft'")
        S, a, b = self._otro(otro)
        return S._envolver(a.mul(b, algorithm))

    def gcd(self, otro) -> "Polynomial":
        """Máximo común divisor mónico (gcd(0, 0) = 0)."""
        S, a, b = self._otro(otro)
        return S._envolver(a.gcd(b))

    def xgcd(self, otro):
        """Euclides extendido: ``(g, s, t)`` con ``g = s·self + t·otro`` y g mónico."""
        S, a, b = self._otro(otro)
        return tuple(S._envolver(x) for x in a.gcd_ext(b))

    def inverse_mod(self, m) -> "Polynomial":
        """Inverso de f módulo m.  ``ValueError`` si no es invertible."""
        S, a, b = self._otro(m)
        if b.is_zero():
            raise ZeroDivisionError("inverso módulo el polinomio cero")
        try:
            return S._envolver(a.inv_mod(b))
        except ValueError:
            raise ValueError(f"{self:ascii} no es invertible módulo {S._envolver(b):ascii}") from None

    def compose(self, g) -> "Polynomial":
        """Composición ``f(g(x))``."""
        S, a, b = self._otro(g)
        return S._envolver(a.compose(b))

    def __call__(self, v):
        """``f(v)``: evaluación si v es convertible al cuerpo; composición si es un
        polinomio o una cadena que contiene la variable."""
        R, K = self._R, self._R._K
        if isinstance(v, Polynomial):
            S = v._R
            return S._envolver(S._nativo(self).compose(v._p))
        if isinstance(v, str):
            usadas = set()
            g = R._analizar(v, usadas)
            if R._var in usadas:
                return R._envolver(self._p.compose(g._p))
            v = g
        elif isinstance(v, FieldElement) and K._n == 1 and v._K._p == K._p and v._K._n > 1:
            # Polinomio sobre F_p evaluado en un elemento de una extensión.
            return self.change_ring(v._K)(v)
        return K._crear(self._p.evaluate(K._valor(v)))

    evaluate = __call__

    def change_ring(self, field: FiniteField) -> "Polynomial":
        """El mismo polinomio sobre otro cuerpo (coeficientes convertidos uno a uno)."""
        if not isinstance(field, FiniteField):
            raise TypeError(f"se esperaba un cuerpo finito, no {type(field).__name__!r}")
        return field[self._R._var](self)

    # -- propiedades y factorización ----------------------------------------
    def monic(self) -> "Polynomial":
        """Polinomio mónico asociado (el cero queda igual)."""
        return self._R._envolver(self._p.monic())

    def derivative(self) -> "Polynomial":
        """Derivada formal."""
        return self._R._envolver(self._p.derivative())

    def is_irreducible(self) -> bool:
        """Test de irreducibilidad de Ben-Or (False para las constantes)."""
        return self._p.is_irreducible()

    def is_squarefree(self) -> bool:
        """¿Libre de cuadrados?  (El cero no lo es; las constantes no nulas sí.)"""
        return self._p.is_squarefree()

    def radical(self) -> "Polynomial":
        """Parte libre de cuadrados mónica: rad(0) = 0, rad(c) = 1."""
        return self._R._envolver(self._p.radical())

    def multiplicity(self, u) -> int:
        """Mayor e tal que u^e divide a f (``ValueError`` si grado(u) < 1 o f = 0)."""
        S, a, b = self._otro(u)
        return a.multiplicity(b)

    def _exigir_libre_de_cuadrados(self) -> None:
        if self._p.degree() >= 1 and not self._p.is_squarefree():
            raise ValueError(f"{self:ascii} no es libre de cuadrados")

    def distinct_degree_factorization(self) -> list:
        """Factorización en grados distintos de un polinomio libre de cuadrados.

        Lista de ``(d, producto de los factores irreducibles mónicos de grado d)``.
        """
        self._exigir_libre_de_cuadrados()
        return [(d, self._R._envolver(u)) for d, u in self._p.distinct_degree()]

    def equal_degree_factorization(self, r) -> list:
        """Factores irreducibles mónicos de un producto libre de cuadrados de
        irreducibles de grado r (Cantor-Zassenhaus)."""
        r = entero(r, "r")
        self._exigir_libre_de_cuadrados()
        return [self._R._envolver(u) for u in self._p.equal_degree(r)]

    def factor(self) -> "Factorization":
        """Factorización en irreducibles mónicos (orden canónico) y unidad."""
        if self._p.is_zero():
            raise ValueError("el polinomio cero no se puede factorizar")
        factores, unidad = self._p.factor()
        R = self._R
        return Factorization([(R._envolver(u), e) for u, e in factores], R._K._crear(unidad), R)

    def roots(self, multiplicities: bool = True) -> list:
        """Raíces en el cuerpo, ordenadas por índice.

        Lista de ``(raíz, multiplicidad)`` o, con ``multiplicities=False``, de
        raíces.  ``ValueError`` para el polinomio cero.
        """
        if self._p.is_zero():
            raise ValueError("el polinomio cero tiene todas las raíces")
        K = self._R._K
        raices = sorted((K._nat.neg(u._p.coeff(0)), e) for u, e in self.factor() if u.degree == 1)
        if multiplicities:
            return [(K._crear(r), e) for r, e in raices]
        return [K._crear(r) for r, _ in raices]

    # -- comparación, hash, pickle -----------------------------------------
    def __eq__(self, otro):
        if otro.__class__ is Polynomial and otro._R is self._R:
            return self._p == otro._p
        try:
            par = self._par(otro)
        except (TypeError, ValueError, ArithmeticError):
            return False
        if par is None:
            return NotImplemented
        return par[1] == par[2]

    def __hash__(self) -> int:
        if self._h is None:
            cs = self._p.coeffs()
            # Las constantes, como su coeficiente (y este como un entero).
            self._h = hash(cs[0] if cs else 0) if len(cs) <= 1 else hash((self._R._var, tuple(cs)))
        return self._h

    def __reduce__(self):
        return (_poli_desde_indices, (self._R, self._p.coeffs()))

    def __copy__(self):
        return self

    def __deepcopy__(self, memo):
        return self

    # -- representación ------------------------------------------------------
    def _texto(self, estilo: str) -> str:
        K = self._R._K

        def formato(v):
            return _texto_elemento(K, v, estilo) if v else None

        return _pretty.polinomio(self._p.coeffs(), self._R._var, estilo, formato)

    def _compuesto(self) -> bool:
        """¿Tiene más de un término?"""
        return sum(1 for c in self._p.coeffs() if c) > 1

    def __repr__(self) -> str:
        return f"{self._R._K!r}[{self._R._var!r}]({self._texto('ascii')!r})"


def _division_exacta(S: PolynomialRing, a, b):
    """``a / b`` con b no nulo: por el inverso si b es constante; si no, exacta."""
    if b.degree() == 0:
        return a.scalar_mul(S._K._nat.inv(b.coeff(0)))
    q, r = a.divmod(b)
    if not r.is_zero():
        raise ValueError(f"la división de {S._envolver(a):ascii} entre {S._envolver(b):ascii} "
                         "no es exacta: usa // y %")
    return q


# ---------------------------------------------------------------------------
# Factorización
# ---------------------------------------------------------------------------
class Factorization(_pretty.Imprimible):
    """Factorización ``unit · f1^e1 · f2^e2 · ...`` de un polinomio.

    Es una secuencia de pares ``(factor irreducible mónico, multiplicidad)`` en
    orden canónico (por grado y coeficientes), con la unidad en ``.unit``.
    ``expand()`` reconstruye el polinomio.
    """

    __slots__ = ("_factores", "_unidad", "_anillo")

    def __init__(self, factors, unit, ring=None):
        """``factors``: pares (polinomio, multiplicidad); ``unit``: elemento del
        cuerpo; ``ring``: anillo (por defecto, el de los factores o K["x"])."""
        self._factores = tuple((f, entero(e, "la multiplicidad")) for f, e in factors)
        if ring is None:
            ring = self._factores[0][0].ring if self._factores else unit.field["x"]
        self._anillo = ring
        self._unidad = ring.field(unit)

    @property
    def unit(self) -> FieldElement:
        """Unidad (coeficiente principal del polinomio factorizado)."""
        return self._unidad

    def __len__(self) -> int:
        return len(self._factores)

    def __iter__(self):
        return iter(self._factores)

    def __getitem__(self, i):
        return self._factores[i]

    @property
    def ring(self) -> PolynomialRing:
        """Anillo del polinomio factorizado."""
        return self._anillo

    def expand(self) -> Polynomial:
        """Producto de la unidad por los factores con sus multiplicidades."""
        f = self._anillo(self._unidad)
        for u, e in self._factores:
            f = f * u ** e
        return f

    def __eq__(self, otro):
        if not isinstance(otro, Factorization):
            return NotImplemented
        return self._unidad == otro._unidad and self._factores == otro._factores

    def __hash__(self) -> int:
        return hash((self._unidad, self._factores))

    def __repr__(self) -> str:
        R = self._anillo
        # El anillo solo hace falta si no se deduce de los factores ni es K["x"].
        anillo = f", ring={R.field!r}[{R.var!r}]" if not self._factores and R.var != "x" else ""
        return f"Factorization({list(self._factores)!r}, unit={self._unidad!r}{anillo})"

    def _texto(self, estilo: str) -> str:
        u = self._unidad
        unidad = _texto_elemento(u.field, u.index, estilo)
        factores = [(f._texto(estilo), f._compuesto(), e) for f, e in self._factores]
        return _pretty.producto(unidad, factores, estilo)
