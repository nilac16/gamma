#pragma once

#ifndef GAMMA_INTERP_H
#define GAMMA_INTERP_H

#include <stdalign.h>
#include "common.h"
#include "idx.h"
#include "vec.h"

EXTERN_C_BEGIN


/** @brief Interpolator buffer, fully aligned to petition a clever compiler to
 *      VECTORIZE EVERYTHING that I do with this
 */
struct gamma_interp {
    alignas (alignof (double) * 8) double buf[8];
};


/** @brief Prepare the interpolator for many evaluations
 *  @param interp
 *      Interpolator with corner values already set
 */
GAMMA_INLINE void gamma_interp_prepare(struct gamma_interp *interp)
{
    interp->buf[4] -= interp->buf[0];
    interp->buf[5] -= interp->buf[1];
    interp->buf[6] -= interp->buf[2];
    interp->buf[7] -= interp->buf[3];

    interp->buf[2] -= interp->buf[0];
    interp->buf[3] -= interp->buf[1];
    interp->buf[6] -= interp->buf[4];
    interp->buf[7] -= interp->buf[5];

    interp->buf[1] -= interp->buf[0];
    interp->buf[3] -= interp->buf[2];
    interp->buf[5] -= interp->buf[4];
    interp->buf[7] -= interp->buf[6];
}


/** @brief Evaluate a prepared interpolator
 *  @param interp
 *      Interpolator that has been `gamma_interp_prepare`'d
 *  @param offs
 *      Unit pixel offset into the lattice cell
 *  @returns The interpolated value
 */
GAMMA_INLINE double gamma_interp_eval(const struct gamma_interp *interp,
                                      const gamma_vec_t         *offs)
{
    struct gamma_interp spill = *interp;

    spill.buf[0] += spill.buf[4] * offs->vec[2];
    spill.buf[1] += spill.buf[5] * offs->vec[2];
    spill.buf[2] += spill.buf[6] * offs->vec[2];
    spill.buf[3] += spill.buf[7] * offs->vec[2];

    spill.buf[0] += spill.buf[2] * offs->vec[1];
    spill.buf[1] += spill.buf[3] * offs->vec[1];

    return spill.buf[0] + spill.buf[1] * offs->vec[0];
}


/** @brief Evaluate an interpolator without preparing it first. Use this if you
 *      only need a single value before it will be destroyed
 *  @param interp
 *      Interpolator with corner values set
 *  @param offs
 *      Unit pixel offset
 *  @returns The interpolated value
 */
GAMMA_INLINE double gamma_interp_single(const struct gamma_interp *interp,
                                        const gamma_vec_t         *offs)
{
    struct gamma_interp spill = *interp;
    gamma_vec_t onem;

    onem = gamma_vec_sub(&(const gamma_vec_t){{1, 1, 1, 0}}, offs);

    spill.buf[0] = spill.buf[0] * onem.vec[2] + spill.buf[4] * offs->vec[2];
    spill.buf[1] = spill.buf[1] * onem.vec[2] + spill.buf[5] * offs->vec[2];
    spill.buf[2] = spill.buf[2] * onem.vec[2] + spill.buf[6] * offs->vec[2];
    spill.buf[3] = spill.buf[3] * onem.vec[2] + spill.buf[7] * offs->vec[2];

    spill.buf[0] = spill.buf[0] * onem.vec[1] + spill.buf[2] * offs->vec[1];
    spill.buf[1] = spill.buf[1] * onem.vec[1] + spill.buf[3] * offs->vec[1];

    return spill.buf[0] * onem.vec[0] + spill.buf[1] * offs->vec[0];
}


EXTERN_C_END

#endif /* GAMMA_INTERP_H */
