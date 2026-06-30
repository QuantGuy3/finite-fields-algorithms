// alcp/fft_tables.hpp — Precomputed FFT/IFFT data, cached per field.
//
// A field that supports an FFT of size n = 2^k (i.e. 2^k divides q - 1) owns a
// cache of FFTTables keyed by k.  Each entry stores the "ints necessary for the
// (I)FFT" precomputed once and reused across every multiplication of that size:
//   * fwd   : forward twiddle factors  g^0, g^1, ..., g^{n/2 - 1}
//   * inv   : inverse twiddle factors  g^{-0}, ..., g^{-(n/2 - 1)}  (for IFFT)
//   * inv_n : n^{-1} in the field      (the IDFT normalisation factor)
#pragma once

#include <map>
#include <vector>

namespace alcp {

template <class Elem>
struct FFTTables {
    unsigned long k = 0;       // transform size is n = 2^k
    std::vector<Elem> fwd;     // forward twiddles, length n/2
    std::vector<Elem> inv;     // inverse twiddles, length n/2
    Elem inv_n{};              // n^{-1}, used by the inverse transform
};

template <class Elem>
using FFTCache = std::map<unsigned long, FFTTables<Elem>>;

}  // namespace alcp
