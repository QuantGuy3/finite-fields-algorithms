// proto_modulo.cpp — Fase 6, tema «python»: módulo pybind11 EXPERIMENTAL mínimo
// (alcp_proto) para medir cuánto cuesta cada capa de una llamada desde Python y
// qué ganan conversiones más baratas.  NO sustituye a python-bindings/bindings.cpp
// ni se usa desde el paquete; solo lo importa mide_python.py.
//
// Contenido (todo sobre Fp64, es decir, F_p con p < 2^62):
//   * Cuerpo(p)                pybind11, guarda std::shared_ptr<const Fp64>.
//       .mul(a, b) / .add(a, b)  elementos con conversión int <-> u64 DIRECTA
//                                (sin mpz): PyLong compacto -> u64.
//       .mul_u64(a, b)           igual, pero con el caster estándar de pybind11
//                                para uint64_t (solo medición: no reduce negativos).
//       .poly(lista)             PolyP con conversión directa (sin vector<BigInt>).
//       .poly_stl(lista)         PolyP con el caster std::vector<uint64_t> de
//                                pybind11 (solo medición: no admite negativos ni
//                                enteros >= 2^64).
//   * PolyP                     pybind11, polinomio opaco: Poly<Fp64> SIN variant
//                                ni comprobación de cuerpo en el constructor.
//       __mul__, __add__, coeffs(), degree(), mul_grado(g) (producto que solo
//       devuelve el grado: aísla el coste de crear el objeto resultado).
//   * PolyH                     HÍBRIDO: clase pybind11 (igual que PolyP) a la que,
//                                tras registrarla, se le sustituye la ranura
//                                nb_multiply por una función C que lee el valor
//                                C++ directamente (sin el despachador de pybind11);
//                                el resultado se sigue creando con py::cast.  Mide
//                                cuánto se gana sin reescribir la clase.  Se crea
//                                con Cuerpo.poly_h(lista).
//   * RawPoly                   tipo CPython «a mano» (PyType_FromSpec), sin
//                                pybind11: nb_multiply / nb_add directos.  Cota
//                                inferior del coste de despacho.  Se crea con
//                                raw_poly(cuerpo, lista).
//   * noop() (pybind11), rnoop() y rmul(a, b) / radd(a, b) (CPython METH_FASTCALL,
//     sobre el cuerpo fijado con fijar_p(p)): suelos del coste de una llamada.
//
// Semántica de la conversión de entrada (idéntica a la de alcp_native para Fp64):
// cualquier entero (o objeto con __index__) se reduce módulo p a [0, p); los float
// se rechazan con TypeError.  Salidas: int en [0, p).
//
// Compilación (desde la raíz del repositorio, con las MISMAS banderas que
// «build.py module» usa para el módulo: -std=c++23 -Wall -Wextra -O3, -isystem de
// pybind11 y de Python, -shared -lpython315 -static -lgmpxx -lgmp):
//   g++ -std=c++23 -Wall -Wextra -I include
//       -isystem .venv/Lib/site-packages/pybind11/include
//       -isystem C:/Users/ismae/AppData/Local/Programs/Python/Python315/Include
//       -O3 bench/experimentos/fase6/python/proto_modulo.cpp
//       -o build/fase6/python/alcp_proto.cp315-win_amd64.pyd
//       -shared -L C:/Users/ismae/AppData/Local/Programs/Python/Python315/libs
//       -lpython315 -static -lgmpxx -lgmp
#include "alcp/alcp.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

namespace py = pybind11;
using alcp::Fp64;
using PolyF = alcp::Poly<Fp64>;
using u64 = std::uint64_t;

namespace {

// ---------------------------------------------------------------------------
// Conversión int de Python -> elemento de F_p (u64 en [0, p)) sin pasar por mpz.
// Devuelve false con una excepción de Python pendiente si falla.
// ---------------------------------------------------------------------------
bool a_elem(PyObject* o, u64 p, u64& out) {
    if (PyLong_CheckExact(o)) {
        // Camino 1: entero «compacto» (|v| < 2^30 en CPython de 64 bits): sin llamadas.
        if (PyUnstable_Long_IsCompact(reinterpret_cast<PyLongObject*>(o))) {
            const Py_ssize_t v = PyUnstable_Long_CompactValue(reinterpret_cast<PyLongObject*>(o));
            if (v >= 0) {
                out = static_cast<u64>(v) % p;
            } else {
                const u64 m = static_cast<u64>(-(v + 1)) + 1;  // |v| sin desbordar
                const u64 r = m % p;
                out = r == 0 ? 0 : p - r;
            }
            return true;
        }
        // Camino 2: cabe en 64 bits con signo.
        int desborde = 0;
        const long long v = PyLong_AsLongLongAndOverflow(o, &desborde);
        if (desborde == 0) {
            if (v == -1 && PyErr_Occurred()) return false;
            if (v >= 0) {
                out = static_cast<u64>(v) % p;
            } else {
                const u64 m = static_cast<u64>(-(v + 1)) + 1;
                const u64 r = m % p;
                out = r == 0 ? 0 : p - r;
            }
            return true;
        }
    } else if (PyFloat_Check(o)) {
        PyErr_Format(PyExc_TypeError,
                     "se esperaba un entero (int o un objeto con __index__), no un valor de tipo '%s'",
                     Py_TYPE(o)->tp_name);
        return false;
    }
    // Camino general (enteros grandes, bool, numpy.int64...): (__index__(o)) % p en
    // Python (resto no negativo, como mpz_fdiv_r en Fp64::from_int).
    PyObject* idx = PyNumber_Index(o);
    if (!idx) return false;
    PyObject* pp = PyLong_FromUnsignedLongLong(p);
    PyObject* r = pp ? PyNumber_Remainder(idx, pp) : nullptr;
    Py_DECREF(idx);
    Py_XDECREF(pp);
    if (!r) return false;
    out = PyLong_AsUnsignedLongLong(r);
    Py_DECREF(r);
    return !PyErr_Occurred();
}

// Lista (o secuencia) de Python -> vector de elementos.
std::vector<u64> lista_a_elems(PyObject* seq, u64 p) {
    PyObject* rapida = PySequence_Fast(seq, "se esperaba una lista de enteros");
    if (!rapida) throw py::error_already_set();
    const Py_ssize_t n = PySequence_Fast_GET_SIZE(rapida);
    PyObject** items = PySequence_Fast_ITEMS(rapida);
    std::vector<u64> v(static_cast<std::size_t>(n));
    for (Py_ssize_t i = 0; i < n; ++i) {
        if (!a_elem(items[i], p, v[static_cast<std::size_t>(i)])) {
            Py_DECREF(rapida);
            throw py::error_already_set();
        }
    }
    Py_DECREF(rapida);
    return v;
}

// Vector de elementos -> list de int (referencia nueva).
PyObject* elems_a_lista(const std::vector<u64>& c) {
    PyObject* l = PyList_New(static_cast<Py_ssize_t>(c.size()));
    if (!l) return nullptr;
    for (std::size_t i = 0; i < c.size(); ++i) {
        PyObject* x = PyLong_FromUnsignedLongLong(c[i]);
        if (!x) {
            Py_DECREF(l);
            return nullptr;
        }
        PyList_SET_ITEM(l, static_cast<Py_ssize_t>(i), x);
    }
    return l;
}

// ---------------------------------------------------------------------------
// Versión pybind11
// ---------------------------------------------------------------------------
// Cuerpo F_p (p < 2^62).  Fp64 valida que p sea primo (ValueError si no).
struct Cuerpo {
    std::shared_ptr<const Fp64> K;
    u64 p;
    explicit Cuerpo(u64 pp) : K(std::make_shared<const Fp64>(alcp::BigInt(std::to_string(pp)))), p(pp) {}
};

// Polinomio opaco: el Poly<Fp64> vive en C++ entre llamadas.
struct PolyP {
    std::shared_ptr<const Fp64> K;  // mantiene vivo el cuerpo al que apunta p
    PolyF p;
};

}  // namespace

// ---------------------------------------------------------------------------
// Versión CPython «a mano» (RawPoly)
// ---------------------------------------------------------------------------
namespace {

struct RawPolyObj {
    PyObject_HEAD
    std::shared_ptr<const Fp64> K;
    PolyF p;
};

PyTypeObject* tipo_raw = nullptr;

PyObject* raw_nuevo(std::shared_ptr<const Fp64> K, PolyF&& q) {
    RawPolyObj* o = PyObject_New(RawPolyObj, tipo_raw);
    if (!o) return nullptr;
    new (&o->K) std::shared_ptr<const Fp64>(std::move(K));
    new (&o->p) PolyF(std::move(q));
    return reinterpret_cast<PyObject*>(o);
}

void raw_dealloc(PyObject* self) {
    auto* o = reinterpret_cast<RawPolyObj*>(self);
    PyTypeObject* t = Py_TYPE(self);
    o->p.~PolyF();
    o->K.~shared_ptr();
    PyObject_Free(self);
    Py_DECREF(t);  // tipo de montículo
}

template <class Op>
PyObject* raw_binaria(PyObject* a, PyObject* b, Op op) {
    if (Py_TYPE(a) != tipo_raw || Py_TYPE(b) != tipo_raw) Py_RETURN_NOTIMPLEMENTED;
    auto* x = reinterpret_cast<RawPolyObj*>(a);
    auto* y = reinterpret_cast<RawPolyObj*>(b);
    try {
        return raw_nuevo(x->K, op(x->p, y->p));
    } catch (const std::exception& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return nullptr;
    }
}

PyObject* raw_mul(PyObject* a, PyObject* b) {
    return raw_binaria(a, b, [](const PolyF& x, const PolyF& y) { return x * y; });
}
PyObject* raw_add(PyObject* a, PyObject* b) {
    return raw_binaria(a, b, [](const PolyF& x, const PolyF& y) { return x + y; });
}
PyObject* raw_coeffs(PyObject* self, PyObject*) { return elems_a_lista(reinterpret_cast<RawPolyObj*>(self)->p.coeffs()); }
PyObject* raw_degree(PyObject* self, PyObject*) { return PyLong_FromLong(reinterpret_cast<RawPolyObj*>(self)->p.degree()); }

PyMethodDef raw_metodos[] = {
    {"coeffs", raw_coeffs, METH_NOARGS, "Coeficientes (lista de int en [0, p))."},
    {"degree", raw_degree, METH_NOARGS, "Grado (-1 para el cero)."},
    {nullptr, nullptr, 0, nullptr},
};

PyType_Slot raw_slots[] = {
    {Py_tp_dealloc, reinterpret_cast<void*>(raw_dealloc)},
    {Py_nb_multiply, reinterpret_cast<void*>(raw_mul)},
    {Py_nb_add, reinterpret_cast<void*>(raw_add)},
    {Py_tp_methods, raw_metodos},
    {Py_tp_doc, const_cast<char*>("Polinomio opaco sobre Fp64 (tipo CPython sin pybind11; experimento).")},
    {0, nullptr},
};

PyType_Spec raw_spec = {"alcp_proto.RawPoly", sizeof(RawPolyObj), 0, Py_TPFLAGS_DEFAULT, raw_slots};

// Cuerpo global de rmul/radd (suelo de coste de una función METH_FASTCALL).
std::shared_ptr<const Fp64> K_global;

PyObject* r_noop(PyObject*, PyObject*) { Py_RETURN_NONE; }

template <class Op>
PyObject* r_elem(PyObject* const* args, Py_ssize_t n, Op op) {
    if (n != 2) {
        PyErr_SetString(PyExc_TypeError, "se esperaban 2 argumentos");
        return nullptr;
    }
    const u64 p = K_global->modulus64();
    u64 a, b;
    if (!a_elem(args[0], p, a) || !a_elem(args[1], p, b)) return nullptr;
    return PyLong_FromUnsignedLongLong(op(*K_global, a, b));
}
PyObject* r_mul(PyObject*, PyObject* const* args, Py_ssize_t n) {
    return r_elem(args, n, [](const Fp64& K, u64 a, u64 b) { return K.mul(a, b); });
}
PyObject* r_add(PyObject*, PyObject* const* args, Py_ssize_t n) {
    return r_elem(args, n, [](const Fp64& K, u64 a, u64 b) { return K.add(a, b); });
}

PyMethodDef r_funciones[] = {
    {"rnoop", r_noop, METH_NOARGS, "No hace nada (CPython METH_NOARGS)."},
    {"rmul", reinterpret_cast<PyCFunction>(reinterpret_cast<void (*)()>(r_mul)), METH_FASTCALL, "a*b en el cuerpo de fijar_p."},
    {"radd", reinterpret_cast<PyCFunction>(reinterpret_cast<void (*)()>(r_add)), METH_FASTCALL, "a+b en el cuerpo de fijar_p."},
    {nullptr, nullptr, 0, nullptr},
};

// --- Híbrido PolyH: clase pybind11 con la ranura nb_multiply sustituida ---
struct PolyH {
    std::shared_ptr<const Fp64> K;
    PolyF p;
};

PyTypeObject* tipo_h = nullptr;

// Puntero al PolyH guardado en una instancia pybind11 (API interna de pybind11:
// instance::get_value_and_holder; válido porque el tipo es exactamente tipo_h).
const PolyH& valor_h(PyObject* o) {
    auto* inst = reinterpret_cast<py::detail::instance*>(o);
    return *static_cast<const PolyH*>(inst->get_value_and_holder().value_ptr());
}

PyObject* h_mul(PyObject* a, PyObject* b) {
    if (Py_TYPE(a) != tipo_h || Py_TYPE(b) != tipo_h) Py_RETURN_NOTIMPLEMENTED;
    try {
        const PolyH& x = valor_h(a);
        const PolyH& y = valor_h(b);
        return py::cast(PolyH{x.K, x.p * y.p}).release().ptr();
    } catch (py::error_already_set& e) {
        e.restore();
        return nullptr;
    } catch (const std::exception& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return nullptr;
    }
}

}  // namespace

PYBIND11_MODULE(alcp_proto, m) {
    m.doc() = "Prototipo de la fase 6 (tema python): conversiones y despacho más baratos. Solo para medir.";

    py::class_<Cuerpo, std::shared_ptr<Cuerpo>> cu(m, "Cuerpo");
    cu.def(py::init<u64>(), py::arg("p"));
    cu.def("mul", [](const Cuerpo& c, py::handle a, py::handle b) {
        u64 x, y;
        if (!a_elem(a.ptr(), c.p, x) || !a_elem(b.ptr(), c.p, y)) throw py::error_already_set();
        return py::reinterpret_steal<py::object>(PyLong_FromUnsignedLongLong(c.K->mul(x, y)));
    });
    cu.def("add", [](const Cuerpo& c, py::handle a, py::handle b) {
        u64 x, y;
        if (!a_elem(a.ptr(), c.p, x) || !a_elem(b.ptr(), c.p, y)) throw py::error_already_set();
        return py::reinterpret_steal<py::object>(PyLong_FromUnsignedLongLong(c.K->add(x, y)));
    });
    cu.def("mul_u64", [](const Cuerpo& c, u64 a, u64 b) { return c.K->mul(a % c.p, b % c.p); });
    cu.def("poly", [](const Cuerpo& c, py::handle lista) {
        return PolyP{c.K, PolyF(*c.K, lista_a_elems(lista.ptr(), c.p))};
    });
    cu.def("poly_stl", [](const Cuerpo& c, std::vector<u64> v) {
        for (auto& x : v) x %= c.p;
        return PolyP{c.K, PolyF(*c.K, std::move(v))};
    });

    py::class_<PolyP> po(m, "PolyP");
    po.def("__mul__", [](const PolyP& a, const PolyP& b) { return PolyP{a.K, a.p * b.p}; }, py::is_operator());
    po.def("__add__", [](const PolyP& a, const PolyP& b) { return PolyP{a.K, a.p + b.p}; }, py::is_operator());
    po.def("mul_grado", [](const PolyP& a, const PolyP& b) { return (a.p * b.p).degree(); });
    po.def("degree", [](const PolyP& a) { return a.p.degree(); });
    po.def("coeffs", [](const PolyP& a) { return py::reinterpret_steal<py::list>(elems_a_lista(a.p.coeffs())); });

    py::class_<PolyH> ph(m, "PolyH");
    ph.def("coeffs", [](const PolyH& a) { return py::reinterpret_steal<py::list>(elems_a_lista(a.p.coeffs())); });
    cu.def("poly_h", [](const Cuerpo& c, py::handle lista) {
        return PolyH{c.K, PolyF(*c.K, lista_a_elems(lista.ptr(), c.p))};
    });
    tipo_h = reinterpret_cast<PyTypeObject*>(ph.ptr());
    tipo_h->tp_as_number->nb_multiply = h_mul;  // tipo de montículo de pybind11
    PyType_Modified(tipo_h);

    m.def("noop", []() {});
    m.def("fijar_p", [](u64 p) { K_global = std::make_shared<const Fp64>(alcp::BigInt(std::to_string(p))); });

    // Tipo y funciones CPython «a mano».
    tipo_raw = reinterpret_cast<PyTypeObject*>(PyType_FromSpec(&raw_spec));
    if (!tipo_raw) throw py::error_already_set();
    m.add_object("RawPoly", py::reinterpret_borrow<py::object>(reinterpret_cast<PyObject*>(tipo_raw)));
    m.def("raw_poly", [](const Cuerpo& c, py::handle lista) {
        PolyF q(*c.K, lista_a_elems(lista.ptr(), c.p));
        PyObject* o = raw_nuevo(c.K, std::move(q));
        if (!o) throw py::error_already_set();
        return py::reinterpret_steal<py::object>(o);
    });
    for (PyMethodDef* d = r_funciones; d->ml_name; ++d) {
        PyObject* f = PyCFunction_New(d, nullptr);
        if (!f) throw py::error_already_set();
        m.add_object(d->ml_name, py::reinterpret_steal<py::object>(f));
    }
}
