// Minimal test harness: CHECK / CHECK_EQ macros with a global pass/fail tally.
#pragma once
#include <cstdio>
#include <string>

namespace check {
inline int g_pass = 0;
inline int g_fail = 0;

inline void report(bool ok, const char* expr, const char* file, int line) {
    if (ok) {
        ++g_pass;
    } else {
        ++g_fail;
        std::printf("  FAIL [%s:%d]: %s\n", file, line, expr);
    }
}

inline int summary(const char* suite) {
    std::printf("[%s] %d passed, %d failed\n", suite, g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
}  // namespace check

#define CHECK(cond) ::check::report((cond), #cond, __FILE__, __LINE__)
#define CHECK_MSG(cond, msg) ::check::report((cond), msg, __FILE__, __LINE__)
