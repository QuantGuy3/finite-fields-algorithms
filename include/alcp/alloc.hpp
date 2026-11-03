// alcp/alloc.hpp — Asignador opcional para GMP con caché de bloques pequeños.
//
// Qué es: unas funciones de memoria para GMP (mp_set_memory_functions) que
// guardan en una lista libre por hilo los bloques pequeños que GMP libera y
// los reutilizan en la siguiente reserva del mismo tamaño.  Cada elemento de
// F_p / F_q es un mpz_class, de modo que casi toda operación aritmética crea
// y destruye temporales con 1-4 limbs; con el malloc de UCRT (HeapAlloc) cada
// par reserva/liberación cuesta decenas de ns.
//
// No cambia ningún algoritmo ni ningún resultado: solo de dónde sale la
// memoria de los mpz.
//
// Activación (es OPCIONAL y está desactivada por defecto):
//   * Explícita: llamar a alcp::alloc::instalar() al principio del programa
//     (y, si se quiere, alcp::alloc::desinstalar() al final).
//   * Por macro: compilar con -DALCP_FAST_ALLOC e incluir esta cabecera en
//     alguna unidad de traducción; sin tocar el código fuente basta con
//     «-DALCP_FAST_ALLOC -include alcp/alloc.hpp».  La instalación se hace
//     durante la inicialización dinámica de esa unidad.
//   alcp.hpp NO incluye esta cabecera: quien no la pida no la paga.
//
// Por qué es seguro mezclar bloques de antes y de después de instalar():
//   * Todo bloque que entrega este asignador es un bloque de std::malloc con
//     EXACTAMENTE el tamaño pedido (la caché es por tamaño exacto, sin
//     redondeo).  Así un bloque reservado por la función por defecto de GMP
//     (malloc) antes de instalar() puede acabar en la caché sin riesgo, y un
//     bloque de la caché puede liberarse con free() tras desinstalar().
//   * GMP pasa a la función de liberación el tamaño con el que se reservó el
//     bloque (contrato documentado de mp_set_memory_functions), y ese es el
//     tamaño que indexa la caché.
//   * Solo se guardan tamaños múltiplos de 8 hasta kMaxBytes (los limbs son de
//     8 bytes): lo demás va directo a malloc/free.
//   * La caché es por hilo (thread_local): sin cerrojos.  Un bloque reservado
//     en un hilo y liberado en otro acaba en la caché del segundo, lo cual es
//     correcto porque todos son bloques normales de malloc.
//   * Cada lista tiene un tope (kTopes, en bytes por clase), y al terminar el hilo se
//     devuelven a free() todos sus bloques.  Tras esa limpieza, las
//     liberaciones que lleguen tarde (destructores estáticos) van a free().
#pragma once

#include <gmp.h>

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace alcp::alloc {

namespace detalle {

inline constexpr std::size_t kGranularidad = 8;     // bytes por limb
inline constexpr std::size_t kMaxBytes = 512;       // 64 limbs (4096 bits)
inline constexpr std::size_t kClases = kMaxBytes / kGranularidad;
// Tope de memoria guardada por clase y por hilo (configurable al compilar).
// El tope en bloques es ese número de bytes entre el tamaño del bloque: así
// las clases pequeñas (1-2 limbs, las de casi todos los elementos) guardan
// miles de bloques, lo bastante para absorber la destrucción de un vector de
// coeficientes entero, y el total por hilo queda acotado por
// kClases * ALCP_FAST_ALLOC_BYTES_POR_CLASE (8 MiB con el valor por defecto).
// (Con un tope fijo de 128 bloques, convertir un polinomio de grado 1000
// desde Python era más lento que con malloc: la caché se desbordaba.)
#ifndef ALCP_FAST_ALLOC_BYTES_POR_CLASE
#define ALCP_FAST_ALLOC_BYTES_POR_CLASE (128 * 1024)
#endif
inline constexpr std::size_t kBytesPorClase = ALCP_FAST_ALLOC_BYTES_POR_CLASE;

// Tope de bloques de cada clase (tabla precalculada: evita una división).
struct Topes {
    std::uint32_t v[kClases];
    constexpr Topes() : v{} {
        for (std::size_t i = 0; i < kClases; ++i) {
            const std::size_t t = kBytesPorClase / ((i + 1) * kGranularidad);
            v[i] = static_cast<std::uint32_t>(t > 0xFFFFFFFFu ? 0xFFFFFFFFu : t);
        }
    }
};
inline constexpr Topes kTopes{};

// Nodo de la lista libre: se guarda dentro del propio bloque libre.
struct Nodo {
    Nodo* sig;
};

// Estado por hilo.  Es trivial (sin constructor ni destructor), de modo que
// el acceso no necesita comprobaciones de inicialización.
struct Cache {
    Nodo* cabeza[kClases];
    std::uint32_t n[kClases];
    bool registrada;   // ¿se ha creado ya la Guardia de este hilo?
    bool cerrada;      // la Guardia ya vació la caché (fin del hilo)
};

inline thread_local Cache cache{};

// Devuelve a free() todos los bloques de la caché del hilo actual.
inline void vaciar() noexcept {
    Cache& c = cache;
    for (std::size_t i = 0; i < kClases; ++i) {
        Nodo* p = c.cabeza[i];
        while (p) {
            Nodo* s = p->sig;
            std::free(p);
            p = s;
        }
        c.cabeza[i] = nullptr;
        c.n[i] = 0;
    }
}

// Objeto con destructor que vacía la caché al terminar el hilo.  Solo se
// construye (primer uso en el hilo) cuando la caché va a guardar su primer
// bloque, para no pagar la comprobación en el camino rápido.
struct Guardia {
    Guardia() noexcept {}
    ~Guardia() {
        vaciar();
        cache.cerrada = true;
    }
};

inline thread_local Guardia guardia;

// Índice de clase de un tamaño cacheable, o kClases si no lo es.
inline std::size_t clase(std::size_t bytes) noexcept {
    if (bytes == 0 || bytes > kMaxBytes || bytes % kGranularidad != 0) return kClases;
    return bytes / kGranularidad - 1;
}

[[noreturn]] inline void sin_memoria(std::size_t bytes) noexcept {
    // Mismo comportamiento que las funciones por defecto de GMP.
    std::fprintf(stderr, "GNU MP: Cannot allocate memory (size=%zu)\n", bytes);
    std::abort();
}

inline void* reservar_malloc(std::size_t bytes) noexcept {
    void* p = std::malloc(bytes);
    if (!p) sin_memoria(bytes);
    return p;
}

inline void* reservar(std::size_t bytes) noexcept {
    const std::size_t i = clase(bytes);
    if (i < kClases) {
        Cache& c = cache;
        if (Nodo* p = c.cabeza[i]) {
            c.cabeza[i] = p->sig;
            --c.n[i];
            return p;
        }
    }
    return reservar_malloc(bytes);
}

inline void liberar(void* ptr, std::size_t bytes) noexcept {
    if (!ptr) return;
    const std::size_t i = clase(bytes);
    if (i < kClases) {
        Cache& c = cache;
        if (!c.cerrada && c.n[i] < kTopes.v[i]) {
            if (!c.registrada) {
                c.registrada = true;
                (void)&guardia;   // primer uso: registra la limpieza del hilo
            }
            Nodo* p = static_cast<Nodo*>(ptr);
            p->sig = c.cabeza[i];
            c.cabeza[i] = p;
            ++c.n[i];
            return;
        }
    }
    std::free(ptr);
}

inline void* realocar(void* ptr, std::size_t viejo, std::size_t nuevo) noexcept {
    if (!ptr) return reservar(nuevo);
    if (viejo == nuevo) return ptr;
    // Si alguno de los dos tamaños es cacheable se pasa por la caché (reserva
    // + copia + liberación): así el bloque viejo se reutiliza y el nuevo puede
    // salir de la caché.  Si no, realloc normal (puede crecer en el sitio).
    if (clase(viejo) < kClases || clase(nuevo) < kClases) {
        void* q = reservar(nuevo);
        std::memcpy(q, ptr, viejo < nuevo ? viejo : nuevo);
        liberar(ptr, viejo);
        return q;
    }
    void* q = std::realloc(ptr, nuevo);
    if (!q) sin_memoria(nuevo);
    return q;
}

}  // namespace detalle

// Instala el asignador con caché como funciones de memoria de GMP (para
// todo el proceso).  Se puede llamar en cualquier momento, incluso con mpz
// vivos reservados con las funciones por defecto.
inline void instalar() noexcept {
    mp_set_memory_functions(&detalle::reservar, &detalle::realocar, &detalle::liberar);
}

// Restaura las funciones por defecto de GMP y devuelve a free() la caché del
// hilo actual (las de otros hilos se vacían cuando esos hilos terminan).
inline void desinstalar() noexcept {
    mp_set_memory_functions(nullptr, nullptr, nullptr);
    detalle::vaciar();
}

// ¿Están instaladas las funciones de esta cabecera?
inline bool instalado() noexcept {
    void* (*r)(std::size_t) = nullptr;
    mp_get_memory_functions(&r, nullptr, nullptr);
    return r == &detalle::reservar;
}

#ifdef ALCP_FAST_ALLOC
namespace detalle {
// Instalación automática durante la inicialización dinámica de la unidad de
// traducción que incluya esta cabecera con ALCP_FAST_ALLOC definida.
inline const bool instalado_por_macro = (instalar(), true);
}  // namespace detalle
#endif

}  // namespace alcp::alloc
