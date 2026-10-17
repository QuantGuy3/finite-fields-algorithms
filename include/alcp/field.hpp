// alcp/field.hpp — Concepto C++20 FiniteField: la interfaz de cuerpo finito que
// consumen Poly<Field> (polynomial.hpp), el motor FFT (fft.hpp) y los
// algoritmos rápidos (fast.hpp).
//
// Recoge exactamente la interfaz que implementan Fp y Fq (ver el comentario de
// cabecera de fp.hpp):
//   using Elem;                                   // semirregular (copiable y
//                                                 //   construible por defecto)
//   Elem  zero() / one() / from_int(BigInt) / random();
//   Elem  add/sub/mul/div(Elem, Elem);  Elem neg/inv(Elem);  Elem pow(Elem, BigInt);
//   bool  eq(Elem, Elem) / is_zero(Elem) / is_one(Elem);
//   BigInt cardinality() / characteristic();
//   std::string to_string(Elem);   Elem pth_root(Elem);
//   unsigned long fft_capacity();  const Elem& fft_base_root();
//   FFTCache<Elem> fft_cache_;    // miembro mutable, propiedad del cuerpo
//
// Las comprobaciones `static_assert(FiniteField<Fp>)` y
// `static_assert(FiniteField<Fq>)` están al final de fp.hpp y fq.hpp.
//
// Además se define detail::fields_match, el comparador de cuerpos que usan
// Poly y fast.hpp para rechazar operandos de cuerpos distintos.
#pragma once

#include "bigint.hpp"
#include "fft_tables.hpp"

#include <concepts>
#include <string>

namespace alcp {

template <class F>
concept FiniteField =
    std::semiregular<typename F::Elem> &&
    requires(const F& K, const typename F::Elem& a, const typename F::Elem& b, const BigInt& n) {
        // constantes y coerción
        { K.zero() } -> std::convertible_to<typename F::Elem>;
        { K.one() } -> std::convertible_to<typename F::Elem>;
        { K.from_int(n) } -> std::convertible_to<typename F::Elem>;
        { K.random() } -> std::convertible_to<typename F::Elem>;
        // aritmética
        { K.add(a, b) } -> std::convertible_to<typename F::Elem>;
        { K.sub(a, b) } -> std::convertible_to<typename F::Elem>;
        { K.neg(a) } -> std::convertible_to<typename F::Elem>;
        { K.mul(a, b) } -> std::convertible_to<typename F::Elem>;
        { K.inv(a) } -> std::convertible_to<typename F::Elem>;
        { K.div(a, b) } -> std::convertible_to<typename F::Elem>;
        { K.pow(a, n) } -> std::convertible_to<typename F::Elem>;
        // predicados
        { K.eq(a, b) } -> std::convertible_to<bool>;
        { K.is_zero(a) } -> std::convertible_to<bool>;
        { K.is_one(a) } -> std::convertible_to<bool>;
        // datos del cuerpo
        { K.cardinality() } -> std::convertible_to<BigInt>;
        { K.characteristic() } -> std::convertible_to<BigInt>;
        // conversión y raíz p-ésima (descomposición libre de cuadrados)
        { K.to_string(a) } -> std::convertible_to<std::string>;
        { K.pth_root(a) } -> std::convertible_to<typename F::Elem>;
        // soporte FFT
        { K.fft_capacity() } -> std::convertible_to<unsigned long>;
        { K.fft_base_root() } -> std::convertible_to<const typename F::Elem&>;
        { K.fft_cache_ } -> std::same_as<FFTCache<typename F::Elem>&>;
    };

namespace detail {

// ¿Son a y b el mismo cuerpo?  Camino rápido: el mismo objeto.  Si no, se usa
// la función libre same_field(a, b) del tipo de cuerpo (Fp y Fq la definen; se
// encuentra por ADL en el punto de instanciación), o en su defecto operator==;
// un tipo de cuerpo sin ninguna de las dos solo es igual a sí mismo.
template <class Field>
bool fields_match(const Field& a, const Field& b) {
    if (&a == &b) return true;
    if constexpr (requires { { same_field(a, b) } -> std::convertible_to<bool>; })
        return static_cast<bool>(same_field(a, b));
    else if constexpr (std::equality_comparable<Field>)
        return static_cast<bool>(a == b);
    else
        return false;
}

}  // namespace detail

}  // namespace alcp
