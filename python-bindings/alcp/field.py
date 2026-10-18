"""Cuerpos finitos ``GF(q)`` y sus elementos.

* :func:`GF` construye (y cachea) el cuerpo primo F_p o la extensión
  F_q = F_p[a]/(g), q = p^n.
* :class:`PrimeField` y :class:`ExtensionField` son los cuerpos; ambos derivan
  de :class:`FiniteField`.
* :class:`FieldElement` es un elemento (inmutable).  Los de F_p son
  :class:`PrimeFieldElement` (con ``__index__``: se comportan como enteros en
  ``range``, índices, ``hash``...); los de una extensión,
  :class:`ExtensionFieldElement`.

Representación interna: un elemento guarda su cuerpo y un entero, que en F_p es
su valor en [0, p) y en una extensión es su índice en base p
(``c0 + c1·p + ... + c_{n-1}·p^(n-1)`` representa ``c0 + c1·a + ...``), igual que
el módulo nativo.  Así los elementos de F_p sumergidos en F_q tienen el mismo
entero.

Conversión (``K(x)``, de más a menos específica):

=================================  ==========================================
elemento del mismo cuerpo           él mismo
elemento de F_p (mismo p)           se sumerge
elemento de otra extensión          solo si está en el subcuerpo primo
``int``, ``bool``, ``numpy.int*``   n·1 (vía ``operator.index``)
``str``                             se analiza: ``"a^2 + 1"``, ``"2a"``, ``"a²"``
polinomio constante                 su coeficiente
``fractions.Fraction``              num · den⁻¹ (``ZeroDivisionError`` si p | den)
``list`` / ``tuple``                coeficientes sobre F_p en orden creciente
                                    (solo en extensiones)
``float``, ``complex``, ``Decimal``  ``TypeError`` (usa enteros o Fraction)
=================================  ==========================================

En las operaciones (``+ - * / ==``) se aplica la misma conversión al otro
operando; un elemento de F_p operado con uno de una extensión de
característica p da un elemento de la extensión.  Mezclar cuerpos distintos
da ``ValueError`` («GF(5) y GF(7) son cuerpos distintos»).

Igualdad y hash: ``==`` compara tras convertir el otro operando, así que
``GF(7)(3) == 3`` y también ``GF(7)(3) == 10``; si el otro valor no es
convertible, ``==`` devuelve False (nunca lanza ``TypeError``/``ValueError``).
El hash de un elemento es el de su valor/índice entero: ``{GF(7)(3), 3}`` tiene
un solo elemento.  (Como en SageMath, ``hash(GF(7)(3)) != hash(10)`` aunque
sean iguales: la igualdad módulo p no es compatible con el hash de ``int``.)
"""
from __future__ import annotations

import itertools
import math
import numbers
import operator
import sys
from typing import TYPE_CHECKING

from . import _parse, _pretty
from ._nativo import nativo
from ._util import entero, factorizar, prime_power_factor, texto_entero

if TYPE_CHECKING:
    from .poly import Polynomial, PolynomialRing

#: ``K.elements()`` (y ``iter(K)``) solo enumera cuerpos de hasta este cardinal.
ENUMERACION_MAXIMA = 10**6

# Marca interna: «el otro operando vive en una extensión de este cuerpo primo;
# hay que llevar este operando a ella y repetir la operación allí».
_PROMOVER = object()

_OPERADORES = {"add": operator.add, "sub": operator.sub,
               "mul": operator.mul, "div": operator.truediv}

# Cachés: cuerpos por parámetros y módulos deterministas por (p, n).
_cuerpos: dict = {}
_modulos: dict = {}


# ---------------------------------------------------------------------------
# Constructor público
# ---------------------------------------------------------------------------
def validar_nombre(var, que: str = "la variable") -> str:
    """Comprueba que ``var`` es un nombre de variable válido (solo letras)."""
    if not isinstance(var, str):
        raise TypeError(f"el nombre de {que} debe ser una cadena, no {type(var).__name__!r}")
    if not var.isalpha():
        raise ValueError(f"el nombre de {que} debe estar formado solo por letras: {var!r}")
    return var


def GF(q, modulus=None, var: str = "a") -> "FiniteField":
    """Cuerpo finito de ``q`` elementos.

    * ``q`` primo: el cuerpo primo F_q (``modulus`` debe ser None; ``var`` se
      ignora).
    * ``q = p^n`` con n >= 2: la extensión F_p[a]/(g), donde ``a`` (de nombre
      ``var``) es la clase de la variable y ``g`` el módulo:

      - ``modulus=None`` (por defecto): módulo DETERMINISTA, el primer
        polinomio mónico irreducible de grado n en este orden de enumeración:
        se recorren las cotas B = 2, 3, ..., p y, para cada una, los
        polinomios ``x^n + c_{n-1}x^{n-1} + ... + c0`` con todas las cifras
        ``c_i < B`` y alguna igual a B - 1, por índice
        ``c0 + c1·B + ... + c_{n-1}·B^(n-1)`` creciente.  (Así se prueban
        primero los coeficientes pequeños, también cuando p es enorme.)  Con
        p = 2 es el orden por índice creciente.  La irreducibilidad la
        comprueba el núcleo (test de Ben-Or).  Ejemplos: GF(4) usa
        x² + x + 1, GF(8) usa x³ + x + 1, GF(9) usa x² + 1, GF(25) usa
        x² + x + 1 y GF(27) usa x³ + 2x + 1.
      - ``modulus="random"``: un módulo mónico irreducible aleatorio
        (reproducible con ``set_seed``).
      - una lista o tupla de coeficientes en orden creciente, una cadena
        (``"x^2 + 1"``, con cualquier variable de una letra) o un polinomio
        sobre F_p.  Se hace mónico; debe ser irreducible y de grado n.

    Los cuerpos se cachean: ``GF(7) is GF(7)`` y dos llamadas con el mismo
    p, el mismo módulo y la misma variable devuelven el mismo objeto.
    """
    q = entero(q, "el orden q")
    p, n = prime_power_factor(q)
    if n == 1:
        if modulus is not None:
            raise ValueError(f"GF({q}) es un cuerpo primo: no admite módulo")
        K = _cuerpos.get(p)
        if K is None:
            K = _cuerpos[p] = PrimeField._construir(p)
        return K
    var = validar_nombre(var)
    nat = None
    if modulus is None:
        g = _modulo_determinista(p, n)
    elif isinstance(modulus, str) and modulus == "random":
        nat = nativo.Fq(q)
        g = tuple(nat.modulus)
    else:
        g = _modulo_dado(modulus, p, n)
    clave = (p, g, var)
    K = _cuerpos.get(clave)
    if K is None:
        K = _cuerpos[clave] = ExtensionField._construir(p, n, g, var, nat)
    return K


def _modulo_determinista(p: int, n: int) -> tuple:
    """Primer mónico irreducible de grado n en el orden documentado en :func:`GF`."""
    g = _modulos.get((p, n))
    if g is not None:
        return g
    Fp = nativo.Fp(p)
    cota = 2
    while True:
        # product(...) recorre las tuplas en orden lexicográfico; la última
        # cifra es c0, así que es el orden por índice creciente en base `cota`.
        for cifras in itertools.product(range(cota), repeat=n):
            if cifras[-1] == 0 or max(cifras) != cota - 1:
                continue        # x | g (reducible) o ya visto con una cota menor
            g = tuple(reversed(cifras)) + (1,)
            if Fp.poly(list(g)).is_irreducible():
                _modulos[(p, n)] = g
                return g
        cota += 1


def _modulo_dado(modulus, p: int, n: int) -> tuple:
    """Valida y normaliza (mónico) un módulo dado por el usuario."""
    Fp = GF(p)
    if isinstance(modulus, str):
        nombres = _parse.letras(modulus)
        if len(nombres) > 1:
            raise ValueError("el módulo debe estar escrito en una sola variable de una "
                             f"letra: {modulus!r}")
        f = Fp[nombres.pop() if nombres else "x"](modulus)
    elif isinstance(modulus, _poly.Polynomial):
        if modulus.field is not Fp:
            raise ValueError(f"el módulo debe ser un polinomio sobre GF({p}), "
                             f"no sobre {modulus.field!r}")
        f = modulus
    elif isinstance(modulus, (list, tuple)):
        f = Fp["x"](modulus)
    else:
        raise TypeError("el módulo debe ser None, 'random', una lista de coeficientes, "
                        f"una cadena o un polinomio, no {type(modulus).__name__!r}")
    if f.degree != n:
        raise ValueError(f"el módulo {f:ascii} tiene grado {f.degree}, pero "
                         f"q = {p}^{n} exige grado {n}")
    if not f.is_irreducible():
        raise ValueError(f"el módulo {f:ascii} es reducible sobre GF({p})")
    return tuple(f.monic()._p.coeffs())


def _elemento(K: "FiniteField", v: int) -> "FieldElement":
    """Reconstruye un elemento a partir de su valor/índice (pickle)."""
    return K.from_index(v)


# ---------------------------------------------------------------------------
# Cuerpos
# ---------------------------------------------------------------------------
class FiniteField(_pretty.Imprimible):
    """Cuerpo finito (base común).  No se instancia directamente: usa :func:`GF`.

    Atributos: ``p``/``characteristic``, ``q``/``order``, ``degree`` (n),
    ``modulus``, ``gen``, ``zero``, ``one``, ``prime_subfield``.  ``K(x)``
    convierte ``x`` en un elemento; ``x in K`` dice si es convertible;
    ``K["x"]`` es el anillo de polinomios K[x].
    """

    __slots__ = ("_nat", "_p", "_n", "_q", "_var", "_clase", "_anillos",
                 "_fact_q1", "_nocuadrado", "_primitivo")
    # numpy: que no trate el cuerpo (iterable) como un array al compararlo con
    # sus escalares (GF(9) == numpy.int64(1) daba un array de 9 booleanos).
    __array_ufunc__ = None

    def __new__(cls, *args, **kwargs):
        raise TypeError("los cuerpos finitos se construyen con alcp.GF(q, ...)")

    def _iniciar(self, nat, p: int, n: int, var, clase) -> None:
        self._nat = nat
        self._p = p
        self._n = n
        self._q = p ** n
        self._var = var
        self._clase = clase
        self._anillos = {}
        self._fact_q1 = None
        self._nocuadrado = None
        self._primitivo = None

    # -- datos -------------------------------------------------------------
    @property
    def p(self) -> int:
        """Característica p."""
        return self._p

    characteristic = p

    @property
    def q(self) -> int:
        """Cardinal q = p^n."""
        return self._q

    order = q

    @property
    def degree(self) -> int:
        """Grado n de la extensión sobre F_p (1 en un cuerpo primo)."""
        return self._n

    @property
    def var(self):
        """Nombre del generador de la extensión (None en un cuerpo primo)."""
        return self._var

    @property
    def is_prime_field(self) -> bool:
        """¿Es un cuerpo primo?"""
        return self._n == 1

    @property
    def prime_subfield(self) -> "PrimeField":
        """Subcuerpo primo GF(p)."""
        return GF(self._p)

    @property
    def zero(self) -> "FieldElement":
        """El 0 del cuerpo."""
        return self._crear(0)

    @property
    def one(self) -> "FieldElement":
        """El 1 del cuerpo."""
        return self._crear(1)

    # -- elementos -----------------------------------------------------------
    def _crear(self, v: int) -> "FieldElement":
        """Elemento de valor/índice ``v`` (ya normalizado en [0, q))."""
        e = object.__new__(self._clase)
        e._K = self
        e._v = v
        return e

    def __call__(self, x=0) -> "FieldElement":
        """Convierte ``x`` en un elemento de este cuerpo (ver la tabla del módulo)."""
        if x.__class__ is self._clase and x._K is self:
            return x
        return self._crear(self._valor(x))

    def __contains__(self, x) -> bool:
        """``x in K``: ¿se puede convertir ``x`` en un elemento de K?"""
        try:
            self._valor(x)
        except (TypeError, ValueError, ArithmeticError):
            return False
        return True

    def from_index(self, i) -> "FieldElement":
        """Elemento de índice ``i`` en [0, q) (en F_p, el de valor i).

        En una extensión, el índice ``c0 + c1·p + ...`` es el elemento
        ``c0 + c1·a + ...``.  ``ValueError`` si ``i`` está fuera de [0, q).
        """
        i = entero(i, "el índice")
        if not 0 <= i < self._q:
            raise ValueError(f"el índice debe estar en [0, {self._q}), no {texto_entero(i)}")
        return self._crear(i)

    def to_index(self, x) -> int:
        """Índice en [0, q) de ``x`` (convertido antes a este cuerpo)."""
        return self(x)._v

    def random(self) -> "FieldElement":
        """Elemento aleatorio uniforme (generador global del núcleo; ver set_seed)."""
        return self._crear(self._nat.random())

    def elements(self):
        """Iterador sobre todos los elementos, por índice creciente.

        ``ValueError`` si q supera :data:`ENUMERACION_MAXIMA` (10**6).
        """
        if self._q > ENUMERACION_MAXIMA:
            raise ValueError(f"{self!r} tiene {self._q} elementos: demasiados para "
                             f"enumerarlos (máximo {ENUMERACION_MAXIMA})")
        return map(self._crear, range(self._q))

    def __iter__(self):
        return self.elements()

    def __len__(self) -> int:
        """Número de elementos.  ``OverflowError`` si q > sys.maxsize: usa ``K.order``."""
        if self._q > sys.maxsize:
            raise OverflowError(f"{self!r} tiene {self._q} elementos, demasiados para len(): "
                                "usa K.order")
        return self._q

    def __bool__(self) -> bool:
        """Un cuerpo nunca es vacío (sin esto, ``bool(K)`` usaría ``len`` y fallaría
        en los cuerpos de más de sys.maxsize elementos)."""
        return True

    # -- anillos de polinomios ------------------------------------------------
    def poly_ring(self, var: str = "x") -> "PolynomialRing":
        """Anillo de polinomios K[var] (cacheado; también ``K["x"]``)."""
        return _poly.PolynomialRing(self, var)

    def __getitem__(self, var) -> "PolynomialRing":
        return self.poly_ring(var)

    # -- grupo multiplicativo ---------------------------------------------------
    def _factores_q1(self) -> dict:
        """Factorización de q - 1 (se calcula una vez)."""
        if self._fact_q1 is None:
            self._fact_q1 = factorizar(self._q - 1)
        return self._fact_q1

    def root_of_unity(self, n) -> "FieldElement":
        """Raíz primitiva n-ésima de la unidad (de orden exactamente n).

        Determinista: ``w = x^((q-1)/n)`` para el elemento ``x`` de menor índice
        con el que ``w`` tiene orden n.  ``ValueError`` si n no divide a q - 1.
        """
        n = entero(n, "n")
        if n < 1 or (self._q - 1) % n:
            tn = texto_entero(n)
            raise ValueError(f"{self!r} no tiene raíces primitivas de la unidad de orden "
                             f"{tn}: {tn} no divide a q - 1 = {self._q - 1}")
        primos = list(factorizar(n))
        e = (self._q - 1) // n
        pot = self._nat.pow
        # Los índices 1..p-1 son los elementos de F_p^*: sus potencias e-ésimas
        # forman el subgrupo de orden m = (p-1)/mcd(p-1, e), y solo hay alguna de
        # orden n si m == n (en un cuerpo primo, siempre).  Si no, se saltan: en
        # una extensión con p grande, recorrerlos uno a uno no terminaría.  El
        # resultado sigue siendo el de menor índice.
        inicio = 1 if (self._p - 1) // math.gcd(self._p - 1, e) == n else self._p
        # Existe un generador del grupo multiplicativo, así que la búsqueda termina.
        candidatas = (pot(i, e) for i in range(inicio, self._q))
        return self._crear(next(w for w in candidatas
                                if all(pot(w, n // r) != 1 for r in primos)))

    def primitive_element(self) -> "FieldElement":
        """Generador del grupo multiplicativo (determinista; ver :meth:`root_of_unity`)."""
        if self._primitivo is None:
            self._primitivo = self.root_of_unity(self._q - 1)
        return self._primitivo

    def _no_cuadrado(self) -> int:
        """Índice del primer no cuadrado (q impar), para Tonelli-Shanks."""
        if self._nocuadrado is None:
            e = (self._q - 1) // 2
            # En una extensión de grado par todo elemento de F_p (índices < p) es
            # un cuadrado: se empieza por el índice p (el generador a).  En grado
            # impar, un no residuo de F_p tampoco es un cuadrado en F_q.
            inicio = self._p if self._n % 2 == 0 else 2
            self._nocuadrado = next(i for i in range(inicio, self._q)
                                    if self._nat.pow(i, e) != 1)
        return self._nocuadrado

    # -- conversión --------------------------------------------------------------
    def _valor(self, x, operando: bool = False):
        """Valor/índice nativo de ``x`` convertido a este cuerpo.

        Con ``operando=True`` (operadores aritméticos) devuelve ``None`` si el
        tipo no es convertible (el operador devolverá NotImplemented) y
        :data:`_PROMOVER` si ``x`` es de una extensión de este cuerpo primo.
        """
        if type(x) is int:
            return x % self._p
        if isinstance(x, FieldElement):
            L = x._K
            if L is self or L == self:
                return x._v
            if operando and self._n == 1 and L._p == self._p:
                return _PROMOVER
            return self._desde_elemento(x, operando)
        if isinstance(x, str):
            return self._analizar(x)._v
        if isinstance(x, _poly.Polynomial):
            if operando:
                return None
            if x.degree > 0:
                raise ValueError(f"solo un polinomio constante se puede convertir en un "
                                 f"elemento de {self!r}, no {x:ascii}")
            return self._valor(x[0])
        if hasattr(type(x), "__index__"):
            try:
                return operator.index(x) % self._p
            except TypeError:
                # Tienen __index__ sin ser enteros: un array de numpy no escalar
                # (se trata como una secuencia) o numpy.bool_ (se rechaza).
                if not hasattr(x, "__iter__"):
                    raise TypeError(f"no se puede convertir un valor de tipo "
                                    f"{type(x).__name__!r} en un elemento de {self!r}") from None
                return None if operando else self._desde_lista(list(x))
        if isinstance(x, numbers.Rational):
            num, den = x.numerator % self._p, x.denominator % self._p
            if den == 0:
                raise ZeroDivisionError(f"el denominador de {texto_entero(x.numerator)}/"
                                        f"{texto_entero(x.denominator)} es múltiplo de p = {self._p}")
            return num * pow(den, -1, self._p) % self._p
        if isinstance(x, numbers.Number):
            raise TypeError(f"no se puede convertir un valor de tipo {type(x).__name__!r} en "
                            f"un elemento de {self!r}: usa enteros o Fraction")
        if isinstance(x, (list, tuple)):
            return self._desde_lista(x)
        if operando:
            return None
        raise TypeError(f"no se puede convertir un valor de tipo {type(x).__name__!r} "
                        f"en un elemento de {self!r}")

    def _desde_elemento(self, x: "FieldElement", operando: bool) -> int:
        """Elemento de otro cuerpo L -> este cuerpo K (L != K)."""
        L = x._K
        if L._p == self._p:
            if L._n == 1:
                return x._v                 # F_p se sumerge en cualquier extensión
            if not operando:
                if x._v < self._p:
                    return x._v             # está en el subcuerpo primo
                if self._n == 1:
                    raise ValueError(f"{x:ascii} no pertenece al subcuerpo primo {self!r}")
        raise ValueError(f"{self!r} y {L!r} son cuerpos distintos")

    def _analizar(self, texto: str) -> "FieldElement":
        variables = {self._var: self.gen} if self._var else {}
        return _parse.analizar(texto, self._crear_entero, variables)

    def _crear_entero(self, n: int) -> "FieldElement":
        return self._crear(n % self._p)

    # -- igualdad, copia, pickle ---------------------------------------------------
    def _clave(self) -> tuple:
        return (self._p, self._modulo(), self._var)

    def __eq__(self, otro):
        if self is otro:
            return True
        if not isinstance(otro, FiniteField):
            return NotImplemented
        return self._clave() == otro._clave()

    def __hash__(self) -> int:
        return hash(("GF",) + self._clave())

    def __copy__(self):
        return self

    def __deepcopy__(self, memo):
        return self


class PrimeField(FiniteField):
    """Cuerpo primo F_p = Z/pZ.  Se construye con ``GF(p)``."""

    __slots__ = ()

    @classmethod
    def _construir(cls, p: int) -> "PrimeField":
        K = object.__new__(cls)
        K._iniciar(nativo.Fp(p), p, 1, None, PrimeFieldElement)
        return K

    def _modulo(self) -> tuple:
        # F_p = F_p[x]/(x - 1): el generador es 1 (convenio de SageMath).
        return (self._p - 1, 1)

    @property
    def modulus(self) -> "Polynomial":
        """Polinomio x - 1 sobre F_p (F_p = F_p[x]/(x - 1); el generador es 1)."""
        return self["x"]([-1, 1])

    @property
    def gen(self) -> "PrimeFieldElement":
        """Generador como F_p-álgebra: el 1."""
        return self.one

    def _desde_lista(self, x):
        raise TypeError(f"una lista o tupla no representa un elemento de {self!r}: "
                        "las listas de coeficientes solo valen en extensiones")

    def __reduce__(self):
        return (GF, (self._p,))

    def __repr__(self) -> str:
        return f"GF({self._p})"

    def _texto(self, estilo: str) -> str:
        return _pretty.nombre_cuerpo(self._p, estilo)


class ExtensionField(FiniteField):
    """Extensión F_q = F_p[a]/(g), q = p^n, n >= 2.  Se construye con ``GF(q, ...)``."""

    __slots__ = ("_g",)

    @classmethod
    def _construir(cls, p: int, n: int, g: tuple, var: str, nat=None) -> "ExtensionField":
        K = object.__new__(cls)
        K._iniciar(nat if nat is not None else nativo.Fq(p, list(g)), p, n, var,
                   ExtensionFieldElement)
        K._g = g
        return K

    def _modulo(self) -> tuple:
        return self._g

    @property
    def modulus(self) -> "Polynomial":
        """Módulo g (mónico e irreducible) como polinomio sobre F_p en la variable x."""
        return GF(self._p)["x"](list(self._g))

    @property
    def gen(self) -> "ExtensionFieldElement":
        """Generador ``a`` (la clase de x en F_p[x]/(g)); su índice es p."""
        return self._crear(self._p)

    def _desde_lista(self, x) -> int:
        Fp = GF(self._p)
        return self._nat.from_coeffs([Fp._valor(c) for c in x])

    def __reduce__(self):
        return (GF, (self._q, list(self._g), self._var))

    def __repr__(self) -> str:
        g = _pretty.polinomio_entero(self._g, "x", "ascii")
        extra = "" if self._var == "a" else f", var={self._var!r}"
        return f"GF({self._q}, modulus={g!r}{extra})"

    def _texto(self, estilo: str) -> str:
        g = _pretty.polinomio_entero(self._g, self._var, estilo)
        var = _pretty.variable(self._var, estilo)
        return (f"{_pretty.nombre_cuerpo(self._q, estilo)} = "
                f"{_pretty.nombre_cuerpo(self._p, estilo)}[{var}]/({g})")


# ---------------------------------------------------------------------------
# Elementos
# ---------------------------------------------------------------------------
class FieldElement(_pretty.Imprimible):
    """Elemento (inmutable) de un cuerpo finito.  Se crea con ``K(x)``.

    Admite ``+ - * /`` y ``**`` (exponente entero; negativo = potencia del
    inverso) con cualquier valor convertible al cuerpo, también por la
    izquierda (``2 + e``, ``"a" * e``, ``Fraction(1, 2) * e``).
    """

    __slots__ = ("_K", "_v")
    # numpy: que sus operadores cedan el paso a los nuestros (2 + e con numpy.int64).
    __array_ufunc__ = None

    def __new__(cls, *args, **kwargs):
        raise TypeError("los elementos se crean convirtiendo en el cuerpo: K(x)")

    # -- datos ------------------------------------------------------------
    @property
    def field(self) -> FiniteField:
        """Cuerpo al que pertenece."""
        return self._K

    parent = field

    @property
    def index(self) -> int:
        """Índice en [0, q) (en F_p, el valor)."""
        return self._v

    @property
    def coeffs(self) -> list:
        """Coeficientes sobre F_p en orden creciente (longitud n), como elementos de GF(p)."""
        K = self._K
        if K._n == 1:
            return [self]
        Fp = GF(K._p)
        return [Fp._crear(c) for c in K._nat.to_coeffs(self._v)]

    # -- operaciones -------------------------------------------------------
    def _binaria(self, otro, nombre: str, reflejada: bool):
        """Operación con un operando de otro tipo o de otro cuerpo."""
        K = self._K
        v = K._valor(otro, True)
        if v is None:
            return NotImplemented
        if v is _PROMOVER:
            # self es de F_p y otro de una extensión L de F_p: se opera en L.
            a = otro._K._crear(self._v)
            op = _OPERADORES[nombre]
            return op(otro, a) if reflejada else op(a, otro)
        a, b = (v, self._v) if reflejada else (self._v, v)
        if nombre == "div" and b == 0:
            raise ZeroDivisionError(f"división por cero en {K!r}")
        return K._crear(getattr(K._nat, nombre)(a, b))

    def __add__(self, otro):
        K = self._K
        if otro.__class__ is self.__class__ and otro._K is K:
            return K._crear(K._nat.add(self._v, otro._v))
        return self._binaria(otro, "add", False)

    def __radd__(self, otro):
        return self._binaria(otro, "add", True)

    def __sub__(self, otro):
        K = self._K
        if otro.__class__ is self.__class__ and otro._K is K:
            return K._crear(K._nat.sub(self._v, otro._v))
        return self._binaria(otro, "sub", False)

    def __rsub__(self, otro):
        return self._binaria(otro, "sub", True)

    def __mul__(self, otro):
        K = self._K
        if otro.__class__ is self.__class__ and otro._K is K:
            return K._crear(K._nat.mul(self._v, otro._v))
        return self._binaria(otro, "mul", False)

    def __rmul__(self, otro):
        return self._binaria(otro, "mul", True)

    def __truediv__(self, otro):
        K = self._K
        if otro.__class__ is self.__class__ and otro._K is K:
            if not otro._v:
                raise ZeroDivisionError(f"división por cero en {K!r}")
            return K._crear(K._nat.div(self._v, otro._v))
        return self._binaria(otro, "div", False)

    def __rtruediv__(self, otro):
        return self._binaria(otro, "div", True)

    def __neg__(self):
        return self._K._crear(self._K._nat.neg(self._v))

    def __pos__(self):
        return self

    def __pow__(self, k, modulo=None):
        if modulo is not None:
            raise TypeError("pow() con tres argumentos no está definido para elementos "
                            "de un cuerpo finito")
        if isinstance(k, FieldElement):
            raise TypeError("el exponente debe ser un entero, no un elemento de un cuerpo")
        k = entero(k, "el exponente")
        if k < 0 and not self._v:
            raise ZeroDivisionError(f"el 0 de {self._K!r} no tiene inverso")
        return self._K._crear(self._K._nat.pow(self._v, k))

    def inverse(self) -> "FieldElement":
        """Inverso multiplicativo.  ``ZeroDivisionError`` para el 0."""
        return self ** -1

    # -- comparación, hash, conversiones -----------------------------------
    def __eq__(self, otro):
        if otro.__class__ is self.__class__ and otro._K is self._K:
            return self._v == otro._v
        if isinstance(otro, _poly.Polynomial):
            return NotImplemented           # lo decide Polynomial.__eq__
        try:
            v = self._K._valor(otro)
        except (TypeError, ValueError, ArithmeticError):
            return False
        return self._v == v

    def __hash__(self) -> int:
        # Igual al hash del entero en F_p (y en el subcuerpo primo de una
        # extensión): así {GF(7)(3), 3} se comporta como un conjunto de números.
        return hash(self._v)

    def __bool__(self) -> bool:
        return self._v != 0

    def __reduce__(self):
        return (_elemento, (self._K, self._v))

    def __copy__(self):
        return self

    def __deepcopy__(self, memo):
        return self

    # -- teoría de cuerpos ---------------------------------------------------
    def frobenius(self, k=1) -> "FieldElement":
        """Imagen por la potencia k-ésima del Frobenius: ``x^(p^k)`` (k módulo n)."""
        k = entero(k, "k") % self._K._n
        return self ** (self._K._p ** k)

    def _subcuerpo(self, v: int) -> "PrimeFieldElement":
        return GF(self._K._p)._crear(v)

    def trace(self) -> "PrimeFieldElement":
        """Traza absoluta ``x + x^p + ... + x^(p^(n-1))`` (elemento de GF(p))."""
        K, t, c = self._K, self._v, self._v
        for _ in range(K._n - 1):
            c = K._nat.pow(c, K._p)
            t = K._nat.add(t, c)
        return self._subcuerpo(t)

    def norm(self) -> "PrimeFieldElement":
        """Norma absoluta ``x^((q-1)/(p-1))`` (elemento de GF(p))."""
        K = self._K
        return self._subcuerpo(K._nat.pow(self._v, (K._q - 1) // (K._p - 1)))

    def multiplicative_order(self) -> int:
        """Orden en el grupo multiplicativo.

        Factoriza q - 1 por división de prueba (``ValueError`` si q - 1 tiene un
        cofactor compuesto sin factores pequeños) y lanza ``ValueError`` para el 0.
        """
        if not self._v:
            raise ValueError("el 0 no tiene orden multiplicativo")
        K = self._K
        orden = K._q - 1
        for r, e in K._factores_q1().items():
            for _ in range(e):
                if K._nat.pow(self._v, orden // r) != 1:
                    break
                orden //= r
        return orden

    def is_square(self) -> bool:
        """¿Es un cuadrado en su cuerpo?  (En característica 2, siempre.)"""
        K = self._K
        return not self._v or K._p == 2 or K._nat.pow(self._v, (K._q - 1) // 2) == 1

    def sqrt(self, all: bool = False):
        """Raíz cuadrada.

        Devuelve la raíz de menor índice; ``ValueError`` si no es un cuadrado.
        Con ``all=True`` devuelve la lista (ordenada por índice) de todas las
        raíces: ``[]`` si no es un cuadrado.  Tonelli-Shanks si q es impar;
        ``x^(q/2)`` en característica 2.
        """
        K = self._K
        if not self.is_square():
            if all:
                return []
            raise ValueError(f"{self:ascii} no es un cuadrado en {K!r}")
        pot, mul = K._nat.pow, K._nat.mul
        if not self._v or K._p == 2:
            r = pot(self._v, K._q // 2)
        else:
            # Tonelli-Shanks: q - 1 = 2^s · t con t impar.
            s, t = 0, K._q - 1
            while t % 2 == 0:
                s, t = s + 1, t // 2
            c = pot(K._no_cuadrado(), t)
            r = pot(self._v, (t + 1) // 2)
            u = pot(self._v, t)
            m = s
            while u != 1:
                # Menor i con u^(2^i) = 1 (0 < i < m).
                i, w = 0, u
                while w != 1:
                    w, i = mul(w, w), i + 1
                b = pot(c, 1 << (m - i - 1))
                r, c = mul(r, b), mul(b, b)
                u, m = mul(u, c), i
        raices = sorted({r, K._nat.neg(r)})
        if all:
            return [K._crear(x) for x in raices]
        return K._crear(raices[0])

    def minimal_polynomial(self, var: str = "x") -> "Polynomial":
        """Polinomio mínimo sobre F_p: ``∏ (x - c)`` para los conjugados distintos c."""
        K = self._K
        conjugados, c = [], self._v
        while c not in conjugados:
            conjugados.append(c)
            c = K._nat.pow(c, K._p)
        f = K._nat.poly([1])
        for c in conjugados:
            f = f * K._nat.poly([K._nat.neg(c), 1])
        # Los coeficientes están en F_p: sus índices son sus valores.
        Fp = GF(K._p)
        return Fp[var]._envolver(Fp._nat.poly(f.coeffs()))

    # -- representación ------------------------------------------------------
    def _texto(self, estilo: str) -> str:
        return _texto_elemento(self._K, self._v, estilo)[0]

    def __repr__(self) -> str:
        return f"{self._K!r}({self._repr_valor()})"


def _texto_elemento(K: FiniteField, v: int, estilo: str):
    """(texto, compuesto) del elemento de valor/índice ``v`` de K."""
    if K._n == 1:
        return str(v), False
    coefs = K._nat.to_coeffs(v)
    texto = _pretty.polinomio_entero(coefs, K._var, estilo)
    return texto, sum(1 for c in coefs if c) > 1


class PrimeFieldElement(FieldElement):
    """Elemento de F_p: además se comporta como un entero (``int``, ``__index__``)."""

    __slots__ = ()

    def __int__(self) -> int:
        return self._v

    def __index__(self) -> int:
        return self._v

    def _repr_valor(self) -> str:
        return str(self._v)

    def __format__(self, spec: str) -> str:
        # Los especificadores numéricos (d, x, 03d...) se aplican al valor.
        if _pretty.estilo_valido(spec) is None and spec[-1:] in ("b", "d", "o", "x", "X", "n"):
            return format(self._v, spec)
        return _pretty.formatear(self, spec)


class ExtensionFieldElement(FieldElement):
    """Elemento de una extensión F_p[a]/(g)."""

    __slots__ = ()

    def __int__(self) -> int:
        if self._v >= self._K._p:
            raise TypeError(f"{self:ascii} no está en el subcuerpo primo GF({self._K._p}): "
                            "no se puede convertir en int")
        return self._v

    def _repr_valor(self) -> str:
        return repr(self._texto("ascii"))


from . import poly as _poly  # noqa: E402  (importación circular: poly usa este módulo)
