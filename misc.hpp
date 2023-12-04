#pragma once

#ifdef UNIX
#define CLR_SCREEN system("clear");
#endif
#ifdef _WIN32
#define CLR_SCREEN system("cls");
#endif
#define OUT_CONDITION(var) (var % 125 == 0)

inline
constexpr size_t __c_ceil(double n) {
    size_t in = static_cast<size_t>(n);
    if (n == static_cast<double>(in)) {
        return in;
    }

    return in + (n > 0.l ? 1ul : 0ul);
}