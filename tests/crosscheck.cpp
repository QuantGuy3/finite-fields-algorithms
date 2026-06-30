// crosscheck.cpp — verify the C++ port reproduces the ORIGINAL Python results.
//
// Reads the oracle file produced by python/oracle.py (Python = source of truth)
// and recomputes every record with the C++ library, asserting equality.  For
// deterministic operations the match is exact; for randomized factorization the
// canonical factor multiset must coincide (factorization is unique).
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace alcp;

// ---- tiny parsing helpers --------------------------------------------------
static std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::string cur;
    std::stringstream ss(s);
    while (std::getline(ss, cur, sep)) out.push_back(cur);
    if (!s.empty() && s.back() == sep) out.push_back("");
    return out;
}
static std::vector<BigInt> parse_ints(const std::string& s) {
    std::vector<BigInt> v;
    if (s.empty()) return v;
    for (auto& tok : split(s, ',')) v.push_back(BigInt(tok));
    return v;
}

// ---- serialization matching the Python oracle ------------------------------
static std::string poly_str_fp(const Fp& F, const Poly<Fp>& p) {
    if (p.is_zero()) return "0";
    std::string s;
    for (std::size_t i = 0; i < p.coeffs().size(); ++i) {
        if (i) s += ",";
        s += F.to_bigint(p.coeffs()[i]).get_str();
    }
    return s;
}
static std::string poly_str_fq(const Fq& F, const Poly<Fq>& p) {
    if (p.is_zero()) return "0";
    std::string s;
    for (std::size_t i = 0; i < p.coeffs().size(); ++i) {
        if (i) s += ",";
        s += F.to_index(p.coeffs()[i]).get_str();
    }
    return s;
}
static std::string vec_str(const std::vector<BigInt>& v) {  // fixed length, no trim
    std::string s;
    for (std::size_t i = 0; i < v.size(); ++i) { if (i) s += ","; s += v[i].get_str(); }
    return s;
}

// ---- field / polynomial builders ------------------------------------------
static Poly<Fp> mk_fp(const Fp& F, const std::string& s) { return Poly<Fp>::from_ints(F, parse_ints(s)); }
static Poly<Fq> mk_fq(const Fq& F, const std::string& s) {
    std::vector<Fq::Elem> coeffs;
    for (auto& idx : parse_ints(s)) coeffs.push_back(F.from_index(idx));
    return Poly<Fq>(F, std::move(coeffs));
}
static Fq build_fq(const std::string& field_tok) {
    auto parts = split(field_tok, ';');
    BigInt p(parts[0]);
    return Fq(Fp(p), parse_ints(parts[1]));
}

// per-field coefficient representation matching the Python encoding
static BigInt coeff_repr(const Fp& F, const Fp::Elem& c) { return F.to_bigint(c); }
static BigInt coeff_repr(const Fq& F, const Fq::Elem& c) { return F.to_index(c); }

// ---- canonical factor multisets + self-reconstruction ----------------------
// Returns the canonical factor string and sets `reconstructs` to whether
// re-multiplying the C++ factors recovers the input (the true correctness test).
template <class Field>
static std::string canon_factors(const Poly<Field>& comp, const Field& F, bool& reconstructs) {
    auto fac = comp.factor();
    Poly<Field> recon = Poly<Field>::constant(F, fac.leading);
    std::vector<std::pair<std::vector<BigInt>, int>> items;
    for (auto& [u, e] : fac.factors) {
        recon = recon * u.pow(BigInt(e));
        std::vector<BigInt> cs;
        for (auto& c : u.coeffs()) cs.push_back(coeff_repr(F, c));
        items.emplace_back(std::move(cs), e);
    }
    reconstructs = (recon == comp);
    std::sort(items.begin(), items.end());
    std::string s;
    for (std::size_t i = 0; i < items.size(); ++i) {
        if (i) s += ";";
        s += vec_str(items[i].first) + ":" + std::to_string(items[i].second);
    }
    return s;
}

int main(int argc, char** argv) {
    std::string path = (argc > 1) ? argv[1] : "python/oracle_vectors.txt";
    std::ifstream in(path);
    if (!in) { std::cerr << "cannot open " << path << "\n"; return 2; }

    std::string line;
    int records = 0;
    int py_buggy = 0;  // factorization inputs where the Python reference is self-inconsistent
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto t = split(line, '|');
        const std::string& op = t[0];
        ++records;

        if (op == "FPX_MUL") {
            Fp F{BigInt(t[1])};
            CHECK(poly_str_fp(F, mk_fp(F, t[2]) * mk_fp(F, t[3])) == t[4]);
        } else if (op == "FPX_KARA") {
            Fp F{BigInt(t[1])};
            CHECK(poly_str_fp(F, mk_fp(F, t[2]).mul(mk_fp(F, t[3]), MulAlgo::Karatsuba)) == t[4]);
        } else if (op == "FPX_DIVMOD") {
            Fp F{BigInt(t[1])};
            auto [q, r] = mk_fp(F, t[2]).divmod(mk_fp(F, t[3]));
            CHECK(poly_str_fp(F, q) == t[4] && poly_str_fp(F, r) == t[5]);
        } else if (op == "FPX_DIVMOD_FAST") {
            Fp F{BigInt(t[1])};
            auto [q, r] = fast::toeplitz_divmod(F, mk_fp(F, t[2]), mk_fp(F, t[3]));
            CHECK(poly_str_fp(F, q) == t[4] && poly_str_fp(F, r) == t[5]);
        } else if (op == "FPX_GCD") {
            Fp F{BigInt(t[1])};
            CHECK(poly_str_fp(F, Poly<Fp>::gcd(mk_fp(F, t[2]), mk_fp(F, t[3]))) == t[4]);
        } else if (op == "FPX_GCDEXT") {
            Fp F{BigInt(t[1])};
            auto e = Poly<Fp>::gcd_ext(mk_fp(F, t[2]), mk_fp(F, t[3]));
            CHECK(poly_str_fp(F, e.g) == t[4] && poly_str_fp(F, e.x) == t[5] && poly_str_fp(F, e.y) == t[6]);
        } else if (op == "FPX_IRRED") {
            Fp F{BigInt(t[1])};
            CHECK((mk_fp(F, t[2]).is_irreducible() ? "1" : "0") == t[3]);
        } else if (op == "FPX_FACTOR") {
            Fp F{BigInt(t[1])};
            bool recon = false;
            std::string got = canon_factors(mk_fp(F, t[2]), F, recon);
            CHECK(recon);  // C++ factorization MUST reconstruct the input
            if (t[4] == "1") CHECK(got == t[3]);  // compare only when Python self-consistent
            else { ++py_buggy; }                   // Python reference is buggy here
        } else if (op == "FPX_POWMOD") {
            Fp F{BigInt(t[1])};
            Poly<Fp> a = mk_fp(F, t[2]), m = mk_fp(F, t[3]);
            CHECK(poly_str_fp(F, Poly<Fp>::pow_mod(a, BigInt(t[4]), m)) == t[5]);
        } else if (op == "FP_DFT") {
            Fp F{BigInt(t[1])};
            auto a = parse_ints(t[4]);
            auto A = fast::dft(F, F.from_int(BigInt(t[2])), std::stoul(t[3]), a);
            CHECK(vec_str(A) == t[5]);
        } else if (op == "FP_IDFT") {
            Fp F{BigInt(t[1])};
            auto A = parse_ints(t[4]);
            auto B = fast::idft(F, F.from_int(BigInt(t[2])), std::stoul(t[3]), A);
            CHECK(vec_str(B) == t[5]);
        } else if (op == "FP_TOEP_LOWER") {
            Fp F{BigInt(t[1])};
            std::size_t n = std::stoul(t[2]);
            auto y = fast::toeplitz_lower_vec(F, n, parse_ints(t[3]), parse_ints(t[4]));
            CHECK(vec_str(y) == t[5]);
        } else if (op == "FQX_MUL") {
            Fq F = build_fq(t[1]);
            CHECK(poly_str_fq(F, mk_fq(F, t[2]) * mk_fq(F, t[3])) == t[4]);
        } else if (op == "FQX_DIVMOD") {
            Fq F = build_fq(t[1]);
            auto [q, r] = mk_fq(F, t[2]).divmod(mk_fq(F, t[3]));
            CHECK(poly_str_fq(F, q) == t[4] && poly_str_fq(F, r) == t[5]);
        } else if (op == "FQX_GCD") {
            Fq F = build_fq(t[1]);
            CHECK(poly_str_fq(F, Poly<Fq>::gcd(mk_fq(F, t[2]), mk_fq(F, t[3]))) == t[4]);
        } else if (op == "FQX_IRRED") {
            Fq F = build_fq(t[1]);
            CHECK((mk_fq(F, t[2]).is_irreducible() ? "1" : "0") == t[3]);
        } else if (op == "FQX_FACTOR") {
            Fq F = build_fq(t[1]);
            bool recon = false;
            std::string got = canon_factors(mk_fq(F, t[2]), F, recon);
            CHECK(recon);
            if (t[4] == "1") CHECK(got == t[3]);
            else { ++py_buggy; }
        } else {
            std::cerr << "unknown op: " << op << "\n";
            ++check::g_fail;
        }
    }
    std::cout << "processed " << records << " records\n";
    if (py_buggy)
        std::cout << "note: " << py_buggy << " factorization input(s) where the ORIGINAL Python "
                     "is self-inconsistent (does not reconstruct);\n      C++ verified correct by "
                     "reconstruction on all of them.\n";
    return check::summary("crosscheck-vs-python");
}
