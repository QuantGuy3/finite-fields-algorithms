// crosscheck.cpp — verify the C++ port reproduces the ORIGINAL Python results.
//
// Reads the oracle file produced by python-bindings/oracle.py (Python = source
// of truth) and recomputes every record with the C++ library, asserting
// equality.  For deterministic operations the match is exact; for randomized
// factorization the canonical factor multiset must coincide (factorization is
// unique).
//
// Uso: xcheck [ruta]   (por defecto build/oracle_vectors.txt)
//
// Detalles de robustez:
//   * se eliminan los '\r' finales (vectores con finales de línea CRLF);
//   * un registro con menos campos de los esperados, un py_ok distinto de 0/1
//     o una excepción C++ cuentan como fallo, con el número de línea; un
//     archivo sin ningún registro también es un fallo;
//   * FPX_FACTOR / FQX_FACTOR: "OP|cuerpo|poli|factores|py_ok".  La
//     factorización C++ se verifica SIEMPRE por reconstrucción (coeficiente
//     principal por producto de factores^multiplicidad == entrada) y porque
//     cada factor es mónico, irreducible y distinto de los demás.  Si
//     py_ok == 1 además se compara con el multiconjunto canónico de Python; si
//     py_ok == 0 (la práctica original no es autoconsistente, se cuelga o lanza
//     una excepción) el campo de factores puede venir vacío y no se compara.
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <algorithm>
#include <cstdio>
#include <exception>
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
    if (parts.size() < 2) throw std::invalid_argument("campo de cuerpo F_q mal formado: " + field_tok);
    BigInt p(parts[0]);
    return Fq(Fp(p), parse_ints(parts[1]));
}

// per-field coefficient representation matching the Python encoding
static BigInt coeff_repr(const Fp& F, const Fp::Elem& c) { return F.to_bigint(c); }
static BigInt coeff_repr(const Fq& F, const Fq::Elem& c) { return F.to_index(c); }

// ---- canonical factor multisets + self-verification ------------------------
// Resultado de factorizar con C++: cadena canónica (mismo formato que el
// oráculo) y verificación independiente de Python.
struct FactorCheck {
    std::string canon;         // "coefs:mult;coefs:mult;..." ordenado
    bool reconstructs = false; // leading * prod(u^e) == entrada
    bool all_irreducible = false;
    bool all_monic = false;
    bool distinct = false;
    bool verified() const { return reconstructs && all_irreducible && all_monic && distinct; }
};

template <class Field>
static FactorCheck check_factors(const Poly<Field>& comp, const Field& F) {
    FactorCheck r;
    auto fac = comp.factor();
    Poly<Field> recon = Poly<Field>::constant(F, fac.leading);
    std::vector<std::pair<std::vector<BigInt>, int>> items;
    r.all_irreducible = true;
    r.all_monic = true;
    for (auto& [u, e] : fac.factors) {
        recon = recon * u.pow(BigInt(e));
        if (!u.is_irreducible()) r.all_irreducible = false;
        if (u.is_zero() || !F.is_one(u.lead())) r.all_monic = false;
        std::vector<BigInt> cs;
        for (auto& c : u.coeffs()) cs.push_back(coeff_repr(F, c));
        items.emplace_back(std::move(cs), e);
    }
    r.reconstructs = (recon == comp);
    std::sort(items.begin(), items.end());
    r.distinct = std::adjacent_find(items.begin(), items.end(),
                                    [](const auto& a, const auto& b) { return a.first == b.first; }) ==
                 items.end();
    for (std::size_t i = 0; i < items.size(); ++i) {
        if (i) r.canon += ";";
        r.canon += vec_str(items[i].first) + ":" + std::to_string(items[i].second);
    }
    return r;
}

// Número mínimo de campos (incluido el código de operación) de cada registro.
static std::size_t expected_fields(const std::string& op) {
    if (op == "FPX_MUL" || op == "FPX_KARA" || op == "FQX_MUL" || op == "FQX_GCD") return 5;
    if (op == "FPX_GCD") return 5;
    if (op == "FPX_DIVMOD" || op == "FPX_DIVMOD_FAST" || op == "FQX_DIVMOD") return 6;
    if (op == "FPX_GCDEXT") return 7;
    if (op == "FPX_IRRED" || op == "FQX_IRRED") return 4;
    if (op == "FPX_FACTOR" || op == "FQX_FACTOR") return 5;
    if (op == "FPX_POWMOD") return 6;
    if (op == "FP_DFT" || op == "FP_IDFT" || op == "FP_TOEP_LOWER") return 6;
    return 0;  // desconocida
}

// Contadores de los registros de factorización.
struct FactorStats {
    int total = 0;          // registros de factorización
    int compared = 0;       // py_ok == 1: comparados con Python
    int py_bad = 0;         // py_ok == 0: Python no autoconsistente / sin resultado
    int py_bad_empty = 0;   //   ... de ellos, con el campo de factores vacío
    int py_bad_cpp_ok = 0;  //   ... de ellos, con la factorización C++ verificada
};

template <class Field>
static void handle_factor(const std::vector<std::string>& t, const Poly<Field>& comp, const Field& F,
                          FactorStats& st, int lineno) {
    const std::string& py_ok = t[4];
    if (py_ok != "0" && py_ok != "1") {
        std::printf("  FAIL [línea %d]: py_ok = '%s' (se esperaba 0 o 1)\n", lineno, py_ok.c_str());
        ++check::g_fail;
        return;
    }
    ++st.total;
    FactorCheck r = check_factors(comp, F);
    CHECK_MSG(r.reconstructs, "la factorización C++ reconstruye la entrada");
    CHECK_MSG(r.all_irreducible, "todos los factores C++ son irreducibles");
    CHECK_MSG(r.all_monic, "todos los factores C++ son mónicos");
    CHECK_MSG(r.distinct, "los factores C++ son distintos");
    if (py_ok == "1") {
        ++st.compared;
        CHECK_MSG(r.canon == t[3], "factorización C++ == factorización Python (py_ok = 1)");
        if (r.canon != t[3]) std::printf("    línea %d: C++ '%s' vs Python '%s'\n", lineno, r.canon.c_str(), t[3].c_str());
    } else {
        ++st.py_bad;
        if (t[3].empty()) ++st.py_bad_empty;
        if (r.verified()) ++st.py_bad_cpp_ok;
    }
}

int main(int argc, char** argv) {
    std::string path = (argc > 1) ? argv[1] : "build/oracle_vectors.txt";
    std::ifstream in(path);
    if (!in) { std::cerr << "cannot open " << path << "\n"; return 2; }

    std::string line;
    int records = 0;
    int lineno = 0;
    FactorStats fst;
    while (std::getline(in, line)) {
        ++lineno;
        // Quita los '\r' finales (archivos con CRLF): si no, el último campo de
        // cada registro no coincide nunca (p. ej. py_ok == "1\r").
        while (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        auto t = split(line, '|');
        const std::string& op = t[0];
        ++records;

        const std::size_t need = expected_fields(op);
        if (need == 0) {
            std::printf("  FAIL [línea %d]: unknown op: %s\n", lineno, op.c_str());
            ++check::g_fail;
            continue;
        }
        if (t.size() < need) {
            std::printf("  FAIL [línea %d]: %s con %zu campos (se esperaban %zu)\n", lineno, op.c_str(),
                        t.size(), need);
            ++check::g_fail;
            continue;
        }

        try {
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
                handle_factor(t, mk_fp(F, t[2]), F, fst, lineno);
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
                handle_factor(t, mk_fq(F, t[2]), F, fst, lineno);
            }
        } catch (const std::exception& ex) {
            std::printf("  FAIL [línea %d]: %s lanzó una excepción: %s\n", lineno, op.c_str(), ex.what());
            ++check::g_fail;
        }
    }
    // Todo por printf (como los FAIL de check.hpp): en MinGW mezclar printf y
    // std::cout puede desordenar la salida.
    std::printf("processed %d records\n", records);
    if (records == 0) {
        // Un archivo vacío suele indicar que el oráculo falló: no es un éxito.
        std::printf("  FAIL: %s no contiene ningún registro\n", path.c_str());
        ++check::g_fail;
    }
    std::printf("factorizaciones: %d registros válidos; %d comparados con Python (py_ok = 1); en "
                "todos se verifica la C++ por reconstrucción e irreducibilidad de cada factor\n",
                fst.total, fst.compared);
    if (fst.py_bad) {
        std::printf("nota: en %d factorización(es) la práctica Python original no da un resultado "
                    "autoconsistente (py_ok = 0", fst.py_bad);
        if (fst.py_bad_empty) std::printf("; %d sin resultado: campo de factores vacío", fst.py_bad_empty);
        std::printf(").\n      La factorización C++ se ha verificado (reconstruye la entrada y sus "
                    "factores son mónicos, irreducibles y distintos) en %d de esas %d.\n",
                    fst.py_bad_cpp_ok, fst.py_bad);
    }
    return check::summary("crosscheck-vs-python");
}
