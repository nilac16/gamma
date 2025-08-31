#pragma once

#ifndef GAMMA_MATH_IDX_H
#define GAMMA_MATH_IDX_H

#include "common.h"
#include <stdalign.h>
#include <stdint.h>

EXTERN_C_BEGIN


typedef int32_t gamma_iscal_t;


/** @brief A multi-index */
typedef struct gamma_idx {
    alignas (alignof (gamma_iscal_t) * 4) gamma_iscal_t idx[4];
} gamma_idx_t;


#define GAMMA_IDX_BINOP(name, op) \
GAMMA_INLINE gamma_idx_t gamma_idx_ ##name(const gamma_idx_t *a, const gamma_idx_t *b) \
{ \
    gamma_idx_t res; \
    res.idx[0] = a->idx[0] op b->idx[0]; \
    res.idx[1] = a->idx[1] op b->idx[1]; \
    res.idx[2] = a->idx[2] op b->idx[2]; \
    res.idx[3] = a->idx[3] op b->idx[3]; \
    return res; \
}


#define GAMMA_IDX_BOOLOP(name, op) \
GAMMA_INLINE gamma_idx_t gamma_idx_ ##name(const gamma_idx_t *a, const gamma_idx_t *b) \
{ \
    gamma_idx_t res; \
    res.idx[0] = a->idx[0] op b->idx[0] ? (gamma_iscal_t)-1 : 0; \
    res.idx[1] = a->idx[1] op b->idx[1] ? (gamma_iscal_t)-1 : 0; \
    res.idx[2] = a->idx[2] op b->idx[2] ? (gamma_iscal_t)-1 : 0; \
    res.idx[3] = a->idx[3] op b->idx[3] ? (gamma_iscal_t)-1 : 0; \
    return res; \
}


GAMMA_IDX_BINOP(add, +)
GAMMA_IDX_BINOP(sub, -)

GAMMA_IDX_BOOLOP(lt, <)
GAMMA_IDX_BOOLOP(leq, <=)
GAMMA_IDX_BOOLOP(gt, >)
GAMMA_IDX_BOOLOP(geq, >=)
GAMMA_IDX_BOOLOP(eq, ==)
GAMMA_IDX_BOOLOP(neq, !=)

GAMMA_IDX_BINOP(or, |)
GAMMA_IDX_BINOP(and, &)


GAMMA_INLINE gamma_iscal_t gamma_idx_any(const gamma_idx_t *a) { return a->idx[0] | a->idx[1] | a->idx[2] | a->idx[3]; }
GAMMA_INLINE gamma_iscal_t gamma_idx_all(const gamma_idx_t *a) { return a->idx[0] & a->idx[1] & a->idx[2] & a->idx[3]; }


/** @brief Hit-test a set of indices with a rectangular [hyper]volume
 *  @param x
 *      Test index
 *  @param lo
 *      Inclusive lower bound
 *  @param hi
 *      Exclusive upper bound
 *  @returns The hit-test boolean vector, any element nonzero if out-of-bounds.
 *      All nonzero elements are GUARANTEED to have their high (i.e. sign) bit
 *      set
 */
GAMMA_INLINE gamma_idx_t gamma_idx_hittest(const gamma_idx_t *x,
                                           const gamma_idx_t *lo,
                                           const gamma_idx_t *hi)
{
    gamma_idx_t cmplo, cmphi;

    cmplo = gamma_idx_lt(x, lo);
    cmphi = gamma_idx_geq(x, hi);
    return gamma_idx_or(&cmplo, &cmphi);
}


EXTERN_C_END

#endif /* GAMMA_MATH_IDX_H */
