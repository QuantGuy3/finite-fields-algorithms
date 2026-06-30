// bindings.cpp — pybind11 wrapper exposing the C++ algebra to Python.
//
// Module `alcp_native`:
//   is_prime, prime_power_factor
//   Fp(p), Fq(q) / Fq(p, modulus)             -- the finite fields
//   PolyFp, PolyFq                            -- polynomials with full algebra,
//                                                irreducibility and factorization
//   dft/idft, toeplitz_* , divmod_toeplitz    -- the Práctica-3 fast algorithms
//
// Fields are held in shared_ptr so that polynomials (which keep a raw Field*)
// always outlive their field.  F_p coefficients are Python ints; F_q
// coefficients are base-p indices (ints in [0, q)), matching cuerpo_fq.elem_de_int.
#include "alcp/alcp.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <string>
#include <vector>

namespace py = pybind11;
using namespace alcp;

// ---- BigInt <-> Python int -------------------------------------------------
static py::int_ to_py(const BigInt& b) {
    return py::reinterpret_steal<py::int_>(PyLong_FromString(b.get_str().c_str(), nullptr, 10));
}
static BigInt to_big(const py::handle& o) {
    return BigInt(py::str(o).cast<std::string>());
}

// ---- per-field coefficient <-> Python int (overloaded for the template) ----
static Fp::Elem coeff_in(const Fp& K, const BigInt& v) { return K.from_int(v); }
static BigInt coeff_out(const Fp& K, const Fp::Elem& e) { return K.to_bigint(e); }
static Fq::Elem coeff_in(const Fq& K, const BigInt& v) { return K.from_index(v); }
static BigInt coeff_out(const Fq& K, const Fq::Elem& e) { return K.to_index(e); }

// ---- field wrappers (own the field via shared_ptr) -------------------------
struct FpW {
    std::shared_ptr<Fp> K;
    explicit FpW(const py::int_& p) : K(std::make_shared<Fp>(to_big(p))) {}
};
struct FqW {
    std::shared_ptr<Fq> K;
    explicit FqW(const py::int_& q) : K(std::make_shared<Fq>(to_big(q))) {}
    FqW(const py::int_& p, const std::vector<py::int_>& mod) {
        std::vector<BigInt> g;
        for (auto& c : mod) g.push_back(to_big(c));
        K = std::make_shared<Fq>(Fp(to_big(p)), std::move(g));
    }
};

// ---- polynomial wrapper (shares the field's lifetime) ----------------------
template <class Field>
struct PolyW {
    std::shared_ptr<Field> K;
    Poly<Field> p;
    PolyW(std::shared_ptr<Field> k, Poly<Field> poly) : K(std::move(k)), p(std::move(poly)) {}

    static PolyW from_coeffs(std::shared_ptr<Field> k, const std::vector<py::int_>& cs) {
        std::vector<typename Field::Elem> coeffs;
        for (auto& c : cs) coeffs.push_back(coeff_in(*k, to_big(c)));
        return PolyW(std::move(k), Poly<Field>(*k, std::move(coeffs)));
    }
    PolyW same(Poly<Field> q) const { return PolyW(K, std::move(q)); }

    std::vector<py::int_> coeffs() const {
        std::vector<py::int_> out;
        for (auto& c : p.coeffs()) out.push_back(to_py(coeff_out(*K, c)));
        return out;
    }
};

// register PolyFp / PolyFq with an identical method surface
template <class Field>
static py::class_<PolyW<Field>> register_poly(py::module& m, const char* name) {
    using PW = PolyW<Field>;
    py::class_<PW> c(m, name);

    c.def("coeffs", &PW::coeffs, "Coefficient list, little-endian (trimmed).")
        .def("degree", [](const PW& a) { return a.p.degree(); })
        .def("is_zero", [](const PW& a) { return a.p.is_zero(); })
        .def("__eq__", [](const PW& a, const PW& b) { return a.p == b.p; })
        .def("__add__", [](const PW& a, const PW& b) { return a.same(a.p + b.p); })
        .def("__sub__", [](const PW& a, const PW& b) { return a.same(a.p - b.p); })
        .def("__neg__", [](const PW& a) { return a.same(-a.p); })
        .def("__mul__", [](const PW& a, const PW& b) { return a.same(a.p * b.p); })
        .def("__floordiv__", [](const PW& a, const PW& b) { return a.same(a.p / b.p); })
        .def("__mod__", [](const PW& a, const PW& b) { return a.same(a.p % b.p); })
        .def("divmod", [](const PW& a, const PW& b) {
            auto [q, r] = a.p.divmod(b.p);
            return py::make_tuple(a.same(q), a.same(r));
        })
        .def("__pow__", [](const PW& a, const py::int_& k) { return a.same(a.p.pow(to_big(k))); })
        .def("mul", [](const PW& a, const PW& b, const std::string& algo) {
            MulAlgo m = MulAlgo::Auto;
            if (algo == "school") m = MulAlgo::Schoolbook;
            else if (algo == "karatsuba") m = MulAlgo::Karatsuba;
            else if (algo == "fft") m = MulAlgo::FFT;
            return a.same(a.p.mul(b.p, m));
        }, py::arg("other"), py::arg("algo") = "auto",
           "Multiply choosing the algorithm: 'auto' (default), 'school', 'karatsuba', 'fft'.")
        .def("derivative", [](const PW& a) { return a.same(a.p.derivative()); })
        .def("monic", [](const PW& a) { return a.same(a.p.monic()); })
        .def("is_irreducible", [](const PW& a) { return a.p.is_irreducible(); })
        .def("gcd", [](const PW& a, const PW& b) { return a.same(Poly<Field>::gcd(a.p, b.p)); })
        .def("gcd_ext", [](const PW& a, const PW& b) {
            auto e = Poly<Field>::gcd_ext(a.p, b.p);
            return py::make_tuple(a.same(e.g), a.same(e.x), a.same(e.y));
        })
        .def("inv_mod", [](const PW& a, const PW& mod) {
            return a.same(Poly<Field>::inv_mod(a.p, mod.p));
        })
        .def("pow_mod", [](const PW& a, const py::int_& k, const PW& mod) {
            return a.same(Poly<Field>::pow_mod(a.p, to_big(k), mod.p));
        })
        .def("multiplicity", [](const PW& a, const PW& u) { return a.p.multiplicity(u.p); })
        .def("factor", [](const PW& a) {
            auto f = a.p.factor();
            py::list factors;
            for (auto& [u, e] : f.factors) factors.append(py::make_tuple(a.same(u), e));
            return py::make_tuple(factors, to_py(coeff_out(*a.K, f.leading)));
        }, "Returns (list of (irreducible PolyW, multiplicity), leading coefficient).")
        .def("to_string", [](const PW& a, const std::string& v) { return a.p.to_string(v); },
             py::arg("var") = "x")
        .def("__repr__", [name](const PW& a) { return std::string(name) + "(" + a.p.to_string() + ")"; });
    return c;
}

PYBIND11_MODULE(alcp_native, m) {
    m.doc() = "Finite fields, polynomial factorization and fast algorithms (C++/GMP core).";

    m.def("is_prime", [](const py::int_& n) { return is_prime(to_big(n)); });
    m.def("prime_power_factor", [](const py::int_& q) {
        auto [p, n] = prime_power_factor(to_big(q));
        return py::make_tuple(to_py(p), to_py(n));
    });

    // ---- Fp ----
    py::class_<FpW>(m, "Fp")
        .def(py::init<const py::int_&>(), py::arg("p"))
        .def_property_readonly("p", [](const FpW& f) { return to_py(f.K->characteristic()); })
        .def_property_readonly("q", [](const FpW& f) { return to_py(f.K->cardinality()); })
        .def("add", [](const FpW& f, const py::int_& a, const py::int_& b) {
            return to_py(f.K->add(f.K->from_int(to_big(a)), f.K->from_int(to_big(b)))); })
        .def("sub", [](const FpW& f, const py::int_& a, const py::int_& b) {
            return to_py(f.K->sub(f.K->from_int(to_big(a)), f.K->from_int(to_big(b)))); })
        .def("mul", [](const FpW& f, const py::int_& a, const py::int_& b) {
            return to_py(f.K->mul(f.K->from_int(to_big(a)), f.K->from_int(to_big(b)))); })
        .def("neg", [](const FpW& f, const py::int_& a) { return to_py(f.K->neg(f.K->from_int(to_big(a)))); })
        .def("inv", [](const FpW& f, const py::int_& a) { return to_py(f.K->inv(f.K->from_int(to_big(a)))); })
        .def("pow", [](const FpW& f, const py::int_& a, const py::int_& k) {
            return to_py(f.K->pow(f.K->from_int(to_big(a)), to_big(k))); })
        .def("random", [](const FpW& f) { return to_py(f.K->random()); })
        .def("poly", [](const FpW& f, const std::vector<py::int_>& cs) {
            return PolyW<Fp>::from_coeffs(f.K, cs); }, py::arg("coeffs"),
            "Build a PolyFp from a little-endian list of integer coefficients.")
        .def("__repr__", [](const FpW& f) { return "Fp(" + f.K->characteristic().get_str() + ")"; });

    // ---- Fq ----
    py::class_<FqW>(m, "Fq")
        .def(py::init<const py::int_&>(), py::arg("q"),
             "F_q with q = p^n; a random monic irreducible modulus is chosen.")
        .def(py::init<const py::int_&, const std::vector<py::int_>&>(), py::arg("p"), py::arg("modulus"),
             "F_q = F_p[x]/(modulus), modulus little-endian (made monic, must be irreducible).")
        .def_property_readonly("p", [](const FqW& f) { return to_py(f.K->characteristic()); })
        .def_property_readonly("n", [](const FqW& f) { return f.K->extension_degree(); })
        .def_property_readonly("q", [](const FqW& f) { return to_py(f.K->cardinality()); })
        .def_property_readonly("modulus", [](const FqW& f) {
            std::vector<py::int_> v;
            for (auto& c : f.K->modulus()) v.push_back(to_py(c));
            return v; })
        .def("from_index", [](const FqW& f, const py::int_& i) { return to_py(f.K->to_index(f.K->from_index(to_big(i)))); })
        .def("add", [](const FqW& f, const py::int_& a, const py::int_& b) {
            return to_py(f.K->to_index(f.K->add(f.K->from_index(to_big(a)), f.K->from_index(to_big(b))))); })
        .def("mul", [](const FqW& f, const py::int_& a, const py::int_& b) {
            return to_py(f.K->to_index(f.K->mul(f.K->from_index(to_big(a)), f.K->from_index(to_big(b))))); })
        .def("inv", [](const FqW& f, const py::int_& a) {
            return to_py(f.K->to_index(f.K->inv(f.K->from_index(to_big(a))))); })
        .def("pow", [](const FqW& f, const py::int_& a, const py::int_& k) {
            return to_py(f.K->to_index(f.K->pow(f.K->from_index(to_big(a)), to_big(k)))); })
        .def("random", [](const FqW& f) { return to_py(f.K->to_index(f.K->random())); })
        .def("poly", [](const FqW& f, const std::vector<py::int_>& cs) {
            return PolyW<Fq>::from_coeffs(f.K, cs); }, py::arg("coeffs"),
            "Build a PolyFq from a little-endian list of base-p element indices.")
        .def("__repr__", [](const FqW& f) {
            return "Fq(q=" + f.K->cardinality().get_str() + ")"; });

    register_poly<Fp>(m, "PolyFp");
    register_poly<Fq>(m, "PolyFq");

    // ---- Práctica-3 fast algorithms over F_p ----
    m.def("dft", [](const FpW& f, const py::int_& g, unsigned long k, const std::vector<py::int_>& a) {
        std::vector<BigInt> av; for (auto& x : a) av.push_back(f.K->from_int(to_big(x)));
        auto A = fast::dft(*f.K, f.K->from_int(to_big(g)), k, av);
        std::vector<py::int_> out; for (auto& x : A) out.push_back(to_py(x));
        return out;
    }, py::arg("Fp"), py::arg("root"), py::arg("k"), py::arg("a"),
       "DFT_{n,root}(a) with n = 2^k and root a 2^k-th root of unity in F_p.");

    m.def("idft", [](const FpW& f, const py::int_& g, unsigned long k, const std::vector<py::int_>& a) {
        std::vector<BigInt> av; for (auto& x : a) av.push_back(f.K->from_int(to_big(x)));
        auto A = fast::idft(*f.K, f.K->from_int(to_big(g)), k, av);
        std::vector<py::int_> out; for (auto& x : A) out.push_back(to_py(x));
        return out;
    }, py::arg("Fp"), py::arg("root"), py::arg("k"), py::arg("a"));

    m.def("toeplitz_lower_vec", [](const FpW& f, std::size_t n,
                                   const std::vector<py::int_>& a, const std::vector<py::int_>& b) {
        std::vector<BigInt> av, bv;
        for (auto& x : a) av.push_back(f.K->from_int(to_big(x)));
        for (auto& x : b) bv.push_back(f.K->from_int(to_big(x)));
        auto y = fast::toeplitz_lower_vec(*f.K, n, av, bv);
        std::vector<py::int_> out; for (auto& x : y) out.push_back(to_py(x));
        return out;
    });
    m.def("toeplitz_lower_inverse", [](const FpW& f, std::size_t n, const std::vector<py::int_>& a) {
        std::vector<BigInt> av; for (auto& x : a) av.push_back(f.K->from_int(to_big(x)));
        auto y = fast::toeplitz_lower_inverse(*f.K, n, av);
        std::vector<py::int_> out; for (auto& x : y) out.push_back(to_py(x));
        return out;
    });
}
