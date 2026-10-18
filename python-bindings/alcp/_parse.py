"""Analizador de expresiones (descendente recursivo) para elementos y polinomios.

Gramática (de menor a mayor precedencia)::

    suma     := termino (('+' | '-') termino)*
    termino  := unario (('*' | '/') unario | implicito)*
    implicito:= potencia            (si empieza por una variable o por '(')
    unario   := ('+' | '-') unario | potencia
    potencia := atomo [('^' | '**') exponente | superindice]
    exponente:= ['+' | '-'] ENTERO | '(' ['+' | '-'] ENTERO ')'
    atomo    := ENTERO | VARIABLE | '(' suma ')'

Admite multiplicación implícita (``2x``, ``3a^2x``, ``(a+1)(x-1)``, ``2 x``),
superíndices Unicode (``x²``, ``x¹⁰``, ``a⁻¹``), ``·``, ``⋅``, ``×`` y ``∗``
como producto, ``÷`` como división, ``−`` (menos Unicode) y espacios.  Un
número no se multiplica implícitamente por otro número (``2 3`` es un error).

Las variables son secuencias de letras; una secuencia se reparte entre los
nombres conocidos tomando en cada paso el más largo que encaja (con
``a`` y ``x``, ``ax`` es ``a·x``).

La expresión se evalúa en la estructura destino con sus propias operaciones:
el llamante da la función que convierte un entero literal, el diccionario de
variables y, opcionalmente, la división.  Los errores de sintaxis son
``ValueError`` con la cadena y un ``^`` bajo la posición del fallo; los errores
de cálculo (división por cero, potencias no válidas) conservan su tipo
(``ZeroDivisionError``, ``ValueError``) y también señalan la posición.
"""
from __future__ import annotations

import operator
from typing import Callable, Dict, List, Mapping, NamedTuple, Optional, Set

# Profundidad máxima de paréntesis (evita RecursionError con entradas hostiles).
PROFUNDIDAD_MAXIMA = 200

_SUPER = {"⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5",
          "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9", "⁻": "-", "⁺": "+"}
_OPERADORES = {"+": "+", "-": "-", "−": "-", "*": "*", "·": "*", "⋅": "*",
               "×": "*", "∗": "*", "/": "/", "÷": "/", "^": "^"}
_DIGITOS = frozenset("0123456789")

# Tipos de token
NUM, VAR, OP, ABRE, CIERRA, SUP, FIN = "número", "variable", "operador", "(", ")", "superíndice", "fin"


class Token(NamedTuple):
    tipo: str
    valor: object
    pos: int


def mensaje(texto: str, pos: int, motivo: str) -> str:
    """``motivo`` seguido de la expresión y un ``^`` bajo la posición ``pos``."""
    # Tabuladores y saltos de línea se muestran como espacios (para alinear el ^).
    linea = "".join(" " if c in "\t\r\n\v\f" else c for c in texto)
    return f"{motivo}\n    {linea}\n    {' ' * pos}^"


def _error(texto: str, pos: int, motivo: str) -> ValueError:
    return ValueError(mensaje(texto, pos, motivo))


def letras(texto: str) -> Set[str]:
    """Conjunto de letras (caracteres alfabéticos) que aparecen en ``texto``."""
    return {c for c in texto if c.isalpha()}


def tokenizar(texto: str, nombres: Mapping[str, object]) -> List[Token]:
    """Divide ``texto`` en tokens; las variables deben estar en ``nombres``."""
    toks: List[Token] = []
    i, n = 0, len(texto)
    # Nombres ordenados del más largo al más corto (coincidencia más larga).
    candidatos = sorted(nombres, key=len, reverse=True)
    while i < n:
        c = texto[i]
        if c.isspace():
            i += 1
        elif c in _DIGITOS:
            j = i
            while j < n and texto[j] in _DIGITOS:
                j += 1
            toks.append(Token(NUM, int(texto[i:j]), i))
            i = j
        elif c in _SUPER:
            j = i
            while j < n and texto[j] in _SUPER:
                j += 1
            cifras = "".join(_SUPER[ch] for ch in texto[i:j])
            try:
                valor = int(cifras)
            except ValueError:
                raise _error(texto, i, f"superíndice no válido: {texto[i:j]!r}") from None
            toks.append(Token(SUP, valor, i))
            i = j
        elif c.isalpha():
            j = i
            while j < n and texto[j].isalpha():
                j += 1
            # Reparto de la secuencia de letras entre los nombres conocidos.
            k = i
            while k < j:
                for nombre in candidatos:
                    if texto.startswith(nombre, k) and k + len(nombre) <= j:
                        toks.append(Token(VAR, nombre, k))
                        k += len(nombre)
                        break
                else:
                    conocidas = ", ".join(sorted(nombres)) or "ninguna"
                    raise _error(texto, k, f"variable desconocida {texto[k:j]!r} "
                                           f"(variables válidas: {conocidas})")
        elif c == "*" and texto.startswith("**", i):
            toks.append(Token(OP, "^", i))
            i += 2
        elif c in _OPERADORES:
            toks.append(Token(OP, _OPERADORES[c], i))
            i += 1
        elif c == "(":
            toks.append(Token(ABRE, c, i))
            i += 1
        elif c == ")":
            toks.append(Token(CIERRA, c, i))
            i += 1
        else:
            pista = " (los números decimales no están admitidos: usa fracciones como 1/2)" \
                if c in ".," else ""
            raise _error(texto, i, f"carácter inesperado {c!r}{pista}")
    toks.append(Token(FIN, None, n))
    return toks


class _Analizador:
    def __init__(self, texto: str, constante: Callable[[int], object],
                 variables: Mapping[str, object],
                 dividir: Callable[[object, object], object],
                 usadas: Optional[Set[str]]):
        self.texto = texto
        self.constante = constante
        self.variables = variables
        self.dividir = dividir
        self.usadas = usadas
        self.toks = tokenizar(texto, variables)
        self.i = 0
        self.profundidad = 0

    # -- utilidades --------------------------------------------------------
    @property
    def actual(self) -> Token:
        return self.toks[self.i]

    def avanzar(self) -> Token:
        t = self.toks[self.i]
        self.i += 1
        return t

    def es_op(self, *ops: str) -> bool:
        t = self.actual
        return t.tipo == OP and t.valor in ops

    def error(self, motivo: str, pos: Optional[int] = None) -> ValueError:
        return _error(self.texto, self.actual.pos if pos is None else pos, motivo)

    def calcular(self, funcion, a, b, pos: int):
        """Aplica una operación que puede fallar; el error señala ``pos``."""
        try:
            return funcion(a, b)
        except (ValueError, ZeroDivisionError) as exc:
            raise type(exc)(mensaje(self.texto, pos, str(exc))) from None

    # -- gramática -----------------------------------------------------------
    def analizar(self):
        if self.actual.tipo == FIN:
            raise self.error("expresión vacía")
        v = self.suma()
        t = self.actual
        if t.tipo != FIN:
            if t.tipo == CIERRA:
                raise self.error("paréntesis ')' sin abrir")
            raise self.error(f"falta un operador antes de {str(t.valor)!r}"
                             if t.tipo == NUM else "símbolo inesperado")
        return v

    def suma(self):
        v = self.termino()
        while self.es_op("+", "-"):
            signo = self.avanzar().valor
            w = self.termino()
            v = v + w if signo == "+" else v - w
        return v

    def termino(self):
        v = self.unario()
        while True:
            t = self.actual
            if t.tipo == OP and t.valor == "*":
                self.avanzar()
                v = v * self.unario()
            elif t.tipo == OP and t.valor == "/":
                self.avanzar()
                v = self.calcular(self.dividir, v, self.unario(), t.pos)
            elif t.tipo in (VAR, ABRE):
                v = v * self.potencia()          # multiplicación implícita
            else:
                return v

    def unario(self):
        # Iterativo (no recursivo): admite cadenas largas de signos.
        negativo = False
        while self.es_op("+", "-"):
            negativo ^= self.avanzar().valor == "-"
        v = self.potencia()
        return -v if negativo else v

    def potencia(self):
        base = self.atomo()
        t = self.actual
        if t.tipo == OP and t.valor == "^":
            self.avanzar()
            return self.calcular(operator.pow, base, self.exponente(), t.pos)
        if t.tipo == SUP:
            self.avanzar()
            return self.calcular(operator.pow, base, t.valor, t.pos)
        return base

    def exponente(self) -> int:
        parentesis = self.actual.tipo == ABRE
        if parentesis:
            self.avanzar()
        signo = 1
        if self.es_op("+", "-"):
            signo = -1 if self.avanzar().valor == "-" else 1
        if self.actual.tipo != NUM:
            raise self.error("el exponente debe ser un número entero")
        k = signo * self.avanzar().valor
        if parentesis:
            self.cerrar()
        return k

    def cerrar(self) -> None:
        if self.actual.tipo != CIERRA:
            raise self.error("falta ')'")
        self.avanzar()

    def atomo(self):
        t = self.actual
        if t.tipo == NUM:
            self.avanzar()
            return self.constante(t.valor)
        if t.tipo == VAR:
            self.avanzar()
            if self.usadas is not None:
                self.usadas.add(t.valor)
            return self.variables[t.valor]
        if t.tipo == ABRE:
            self.profundidad += 1
            if self.profundidad > PROFUNDIDAD_MAXIMA:
                raise self.error("expresión demasiado anidada")
            self.avanzar()
            v = self.suma()
            self.cerrar()
            self.profundidad -= 1
            return v
        raise self.error("se esperaba un número, una variable o '('")


def analizar(texto: str, constante: Callable[[int], object],
             variables: Optional[Mapping[str, object]] = None,
             dividir: Callable[[object, object], object] = operator.truediv,
             usadas: Optional[Set[str]] = None):
    """Evalúa ``texto`` con las operaciones de los valores de la estructura destino.

    ``constante(n)`` convierte un entero literal ``n >= 0``; ``variables``
    asocia cada nombre a su valor; ``dividir(a, b)`` implementa ``/``.  Si se da
    ``usadas``, se le añaden los nombres de variable que aparecen.
    """
    vars_: Dict[str, object] = dict(variables or {})
    return _Analizador(texto, constante, vars_, dividir, usadas).analizar()
