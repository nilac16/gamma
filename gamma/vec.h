#pragma once

/** @file Basic linear algebra operations as needed by this implementation of
 *      gamma analysis. Care has been taken with type definitions such that they
 *      might be heavily vectorized by platform compilers given proper
 *      architecture flags, all with the standard C header `stdalign.h`
 */

#ifndef GAMMA_MATH_VEC_H
#define GAMMA_MATH_VEC_H

#include "common.h"
#include <stdalign.h>

EXTERN_C_BEGIN


/** @brief Vector scalar component */
typedef double gamma_scal_t;


/** @brief An array of scalar types, suitably aligned for vector operations,
 *      should the target machine support them
 */
typedef struct gamma_vec {
    alignas (alignof (gamma_scal_t) * 4) gamma_scal_t vec[4];
} gamma_vec_t;


#define GAMMA_VEC_BINOP(name, op) \
GAMMA_INLINE gamma_vec_t gamma_vec_ ##name(const gamma_vec_t *a, const gamma_vec_t *b) \
{ \
    gamma_vec_t res; \
    res.vec[0] = a->vec[0] op b->vec[0]; \
    res.vec[1] = a->vec[1] op b->vec[1]; \
    res.vec[2] = a->vec[2] op b->vec[2]; \
    res.vec[3] = a->vec[3] op b->vec[3]; \
    return res; \
}


#define GAMMA_VEC_BINOPS(name, op) \
GAMMA_INLINE gamma_vec_t gamma_vec_ ##name ##s(const gamma_vec_t *a, gamma_scal_t b) \
{ \
    gamma_vec_t res; \
    res.vec[0] = a->vec[0] op b; \
    res.vec[1] = a->vec[1] op b; \
    res.vec[2] = a->vec[2] op b; \
    res.vec[3] = a->vec[3] op b; \
    return res; \
}


#define GAMMA_VEC_TERNOP(name, op1, op2) \
GAMMA_INLINE gamma_vec_t gamma_vec_ ##name(const gamma_vec_t *a, const gamma_vec_t *b, const gamma_vec_t *c) \
{ \
    gamma_vec_t res; \
    res.vec[0] = c->vec[0] op2 a->vec[0] op1 b->vec[0]; \
    res.vec[1] = c->vec[1] op2 a->vec[1] op1 b->vec[1]; \
    res.vec[2] = c->vec[2] op2 a->vec[2] op1 b->vec[2]; \
    res.vec[3] = c->vec[3] op2 a->vec[3] op1 b->vec[3]; \
    return res; \
}


#define GAMMA_VEC_TERNOPS(name, op1, op2) \
GAMMA_INLINE gamma_vec_t gamma_vec_ ##name ##s(const gamma_vec_t *a, gamma_scal_t b, const gamma_vec_t *c) \
{ \
    gamma_vec_t res; \
    res.vec[0] = c->vec[0] op2 a->vec[0] op1 b; \
    res.vec[1] = c->vec[1] op2 a->vec[1] op1 b; \
    res.vec[2] = c->vec[2] op2 a->vec[2] op1 b; \
    res.vec[3] = c->vec[3] op2 a->vec[3] op1 b; \
    return res; \
}


GAMMA_VEC_BINOP(add, +)
GAMMA_VEC_BINOP(sub, -)
GAMMA_VEC_BINOP(mul, *)
GAMMA_VEC_BINOP(div, /)

GAMMA_VEC_BINOPS(mul, *)
GAMMA_VEC_BINOPS(div, /)

GAMMA_VEC_TERNOP(fmadd, *, +)
GAMMA_VEC_TERNOP(fmsub, *, -)

GAMMA_VEC_TERNOPS(fmadd, *, +)
GAMMA_VEC_TERNOPS(fmsub, *, -)


/** @brief Compute the dot product of two vectors */
GAMMA_INLINE
gamma_scal_t gamma_vec_dp(const gamma_vec_t *a, const gamma_vec_t *b)
{
    return a->vec[0] * b->vec[0]
         + a->vec[1] * b->vec[1]
         + a->vec[2] * b->vec[2]
         + a->vec[3] * b->vec[3];
}


GAMMA_INLINE
gamma_vec_t gamma_vec_cross(const gamma_vec_t *a, const gamma_vec_t *b)
{
    gamma_vec_t res = { 0 };

    res.vec[0] = a->vec[1] * b->vec[2] - a->vec[2] * b->vec[1];
    res.vec[1] = a->vec[2] * b->vec[0] - a->vec[0] * b->vec[2];
    res.vec[2] = a->vec[0] * b->vec[1] - a->vec[1] * b->vec[0];
    return res;
}


EXTERN_C_END

#endif /* GAMMA_MATH_VEC_H */
