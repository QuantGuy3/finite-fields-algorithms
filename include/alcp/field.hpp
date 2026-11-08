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

// ---------------------------------------------------------------------------
// Ganchos OPCIONALES de acumulación perezosa (reducción diferida).
//
// Un cuerpo puede ofrecer un tipo acumulador Acc capaz de sumar productos a·b
// SIN reducirlos, y reducir una sola vez al final.  Poly<Field> los detecta en
// tiempo de compilación (if constexpr) en sus bucles internos (mul_school,
// divmod): un cuerpo sin ganchos sigue el camino de siempre (mul + add por
// término).  Los ganchos no cambian el algoritmo: la suma de los productos es
// la misma, solo se reduce módulo p (y módulo g en F_q) una vez por
// coeficiente en lugar de una vez por término.
//
//   using Acc;                              // construible por defecto
//   void acc_clear(Acc& s);                 // s = 0 (y lo dimensiona si hace falta)
//   void acc_set(Acc& s, Elem a);           // s = a
//   void acc_fma(Acc& s, Elem a, Elem b);   // s += a·b   (sin reducir)
//   Elem acc_reduce(Acc& s);                // valor reducido de s; puede consumir s
//                                           //   (después solo vale acc_clear/acc_set)
// ---------------------------------------------------------------------------
template <class F>
concept LazyAccumulation =
    std::default_initializable<typename F::Acc> &&
    requires(const F& K, typename F::Acc& s, const typename F::Elem& a) {
        K.acc_clear(s);
        K.acc_set(s, a);
        K.acc_fma(s, a, a);
        { K.acc_reduce(s) } -> std::convertible_to<typename F::Elem>;
    };

namespace detail {

// Fase 6: despacho de acc_fma fuera del bucle.  Un cuerpo puede ofrecer además
//   bool acc_pequeno();                      // ¿vale acc_fma_peq?
//   void acc_fma_peq(Acc& s, Elem a, Elem b);  // s += a·b, versión rápida
// (Fp64 con p < 2^32: producto de 64 bits sin vigilar el acarreo).  `bucle`
// es una lambda genérica que recibe la operación fma; se instancia dos veces
// y la comprobación se hace una sola vez por bucle, no una por término (con la
// comprobación dentro de acc_fma, los bucles de p grande iban un 15 % más
// lentos).  Ambas versiones dan el mismo valor exacto.
template <class F, class Bucle>
inline void con_acc_fma(const F& K, Bucle&& bucle) {
    if constexpr (requires(typename F::Acc& s, const typename F::Elem& a) {
                      { K.acc_pequeno() } -> std::convertible_to<bool>;
                      K.acc_fma_peq(s, a, a);
                  }) {
        if (K.acc_pequeno()) {
            bucle([&K](typename F::Acc& s, const typename F::Elem& a, const typename F::Elem& b) {
                K.acc_fma_peq(s, a, b);
            });
            return;
        }
    }
    bucle([&K](typename F::Acc& s, const typename F::Elem& a, const typename F::Elem& b) { K.acc_fma(s, a, b); });
}

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
