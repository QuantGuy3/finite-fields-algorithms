// alcp/errors.hpp — Jerarquía de excepciones de la librería.
//
// Convención (ver docs/diseno.md, sección 3):
//   * alcp::division_by_zero (deriva de std::domain_error): invertir el 0 de un
//     cuerpo (Fp::inv, Fq::inv, pow con exponente negativo de 0), dividir por el
//     polinomio cero (Poly::divmod, /, %, fast::toeplitz_divmod) o invertir una
//     matriz de Toeplitz triangular singular (a[0] == 0).  El binding pybind11
//     la traduce a ZeroDivisionError.
//   * std::invalid_argument: argumentos inválidos (p no primo, módulo reducible
//     o de grado < 1, longitudes incorrectas, operandos de cuerpos distintos,
//     raíces de la unidad de orden incorrecto, ...).  -> ValueError.
//   * std::domain_error (resto): p. ej. "polinomio no invertible módulo m",
//     potencia negativa de un polinomio o factorizar el polinomio cero.
//
// Como division_by_zero deriva de std::domain_error, el código antiguo que
// capturaba std::domain_error sigue funcionando sin cambios.
#pragma once

#include <stdexcept>

namespace alcp {

// División por cero en un cuerpo finito o en su anillo de polinomios.
struct division_by_zero : std::domain_error {
    using std::domain_error::domain_error;
};

}  // namespace alcp
