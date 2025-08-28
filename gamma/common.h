#pragma once

#ifndef GAMMA_COMMON_H
#define GAMMA_COMMON_H

/** @brief Declaration specifier for (small) functions defined in-header */
#define GAMMA_INLINE static inline

/** @brief Also known as `ARRAYLEN` when other people write it */
#define BUFLEN(buf) (sizeof (buf) / sizeof *(buf))

/** Create all the extern "C"-type macros if not already defined */
#if defined(__cplusplus) && __cplusplus
#   if !defined(EXTERN_C)
#       define EXTERN_C extern "C"
#   endif

#   if !defined(EXTERN_C_BEGIN)
#       define EXTERN_C_BEGIN extern "C" {
#   endif

#   if !defined(EXTERN_C_END)
#       define EXTERN_C_END }
#   endif
#else
#   if !defined(EXTERN_C)
#       define EXTERN_C
#   endif

#   if !defined(EXTERN_C_BEGIN)
#       define EXTERN_C_BEGIN
#   endif

#   if !defined(EXTERN_C_END)
#       define EXTERN_C_END
#   endif
#endif



#if !defined(__cplusplus) && !__cplusplus
GAMMA_INLINE float gamma_sqrf(float x) { return x * x; }
GAMMA_INLINE double gamma_sqr(double x) { return x * x; }
GAMMA_INLINE long double gamma_sqrl(long double x) { return x * x; }

#   define gamma_sqr(x) _Generic(x, \
        float:       gamma_sqrf, \
        double:      gamma_sqr, \
        long double: gamma_sqrl \
    )(x)

#else
template <class T>
GAMMA_INLINE T gamma_sqr(T x) { return x * x; }

#endif


#endif /* GAMMA_COMMON_H */
