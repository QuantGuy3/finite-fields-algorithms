// alcp/fft_tables.hpp — Precomputed FFT/IFFT data, cached per field.
//
// A field that supports an FFT of size n = 2^k (i.e. 2^k divides q - 1) owns a
// cache of FFTTables keyed by k.  Each entry stores the "ints necessary for the
// (I)FFT" precomputed once and reused across every multiplication of that size:
//   * fwd   : forward twiddle factors  g^0, g^1, ..., g^{n/2 - 1}
//   * inv   : inverse twiddle factors  g^{-0}, ..., g^{-(n/2 - 1)}  (for IFFT)
//   * inv_n : n^{-1} in the field      (the IDFT normalisation factor)
//
// Además (fase 3, solo organización de datos; los valores son los mismos):
//   * fwd_st / inv_st : los twiddles de cada etapa en posiciones CONTIGUAS.  La
//     etapa de longitud len = 2^s (s = 1..k-1) usa w_len^j = tw[j·(n/len)],
//     j = 0..len/2-1; aquí se guardan seguidos a partir del desplazamiento
//     len/2 - 1, de modo que la mariposa los recorre secuencialmente en vez de
//     saltar con paso n/len por toda la tabla.  La última etapa (len = n, paso 1)
//     usa directamente fwd / inv, que ya son contiguos.
//   * swaps : la permutación de inversión de bits precalculada como la lista de
//     pares (i, j) con i < j que hay que intercambiar.
//
// Caché: un vector indexado por k (acceso directo, sin búsqueda en un árbol)
// de punteros compartidos a tablas inmutables.  Compartir es seguro porque una
// tabla no cambia tras construirse, y permite copiar un cuerpo (p. ej. el Fp
// base que guarda Fq) sin duplicar las tablas.
#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace alcp {

// Copia en palabras de máquina de los MISMOS twiddles, para F_p con p < 2^62
// (ver fft.hpp, «NTT en palabras de máquina»).  p == 0 significa «no
// disponible».  Los twiddles de todas las etapas (incluida la última) van
// contiguos: la etapa len empieza en len/2 - 1.  *_sh son los cocientes de
// Shoup floor(w·2^64 / p) que permiten multiplicar por w módulo p sin dividir.
// Fase 7 (convolución AVX2, p < 2^30): copia en 32 bits de los MISMOS
// twiddles de WordNTTTables, con cocientes de Shoup de 32 bits
// floor(w·2^32 / p).  p == 0 significa «no disponible» (CPU sin AVX2,
// p >= 2^30 o n < 16).  f4/f2 (i4/i2) son los twiddles directos (inversos) de
// las etapas half = 4 y half = 2 repetidos para llenar 8 carriles.
struct Word32NTTTables {
    std::uint32_t p = 0, pinv = 0, c = 0, c_sh = 0;  // pinv = -p^{-1} mod 2^32; c = n^{-1}·2^32 mod p
    std::vector<std::uint32_t> fwd, fwd_sh, inv, inv_sh;  // n - 1 entradas, por etapas
    std::uint32_t f4[8] = {}, f4_sh[8] = {}, f2[8] = {}, f2_sh[8] = {};
    std::uint32_t i4[8] = {}, i4_sh[8] = {}, i2[8] = {}, i2_sh[8] = {};
};

struct WordNTTTables {
    std::uint64_t p = 0;
    std::vector<std::uint64_t> fwd, fwd_sh, inv, inv_sh;  // n - 1 entradas
    std::uint64_t inv_n = 0, inv_n_sh = 0;
    // Fase 7 (convolución radix 4): constantes del producto punto a punto de
    // Montgomery fusionado con n^{-1}.  pinv = -p^{-1} mod 2^64 (0 si p es par:
    // entonces se usa el camino anterior); mont_c = n^{-1}·2^64 mod p y su
    // cociente de Shoup.  Se derivan de inv_n: no se calcula ninguna raíz.
    std::uint64_t pinv = 0, mont_c = 0, mont_c_sh = 0;
    Word32NTTTables w32;  // fase 7: camino AVX2 (w32.p == 0 si no se usa)
};

// Versión multi-limb (F_p con p >= 2^62): cada valor ocupa N limbs de 64 bits
// contiguos.  Los twiddles se guardan en forma de Montgomery w·R mod p
// (R = 2^(64N)), de modo que REDC(y · w·R) = y·w mod p ya canónico.  N == 0
// significa «no disponible».
struct LimbNTTTables {
    std::size_t N = 0;                          // limbs por valor
    std::vector<std::uint64_t> p;               // el primo, N limbs
    std::uint64_t pinv = 0;                     // -p^{-1} mod 2^64
    std::vector<std::uint64_t> fwd, inv;        // (n - 1)·N limbs, por etapas
    std::vector<std::uint64_t> inv_n;           // n^{-1}·R mod p, N limbs
};

template <class Elem>
struct FFTTables {
    unsigned long k = 0;       // transform size is n = 2^k
    std::vector<Elem> fwd;     // forward twiddles, length n/2
    std::vector<Elem> inv;     // inverse twiddles, length n/2
    Elem inv_n{};              // n^{-1}, used by the inverse transform
    std::vector<Elem> fwd_st;  // twiddles directos por etapa (contiguos), n/2 - 1
    std::vector<Elem> inv_st;  // twiddles inversos por etapa (contiguos), n/2 - 1
    std::vector<std::pair<std::size_t, std::size_t>> swaps;  // inversión de bits
    WordNTTTables word;        // solo F_p con p < 2^62 (si no, word.p == 0)
    LimbNTTTables limb;        // solo F_p con p >= 2^62 (si no, limb.N == 0)
};

// Caché de tablas de un cuerpo: posición k -> tablas de tamaño 2^k (nullptr si
// todavía no se han construido).
template <class Elem>
using FFTCache = std::vector<std::shared_ptr<const FFTTables<Elem>>>;

}  // namespace alcp
