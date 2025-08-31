#pragma once

#ifndef GAMMA_PSEARCH_H
#define GAMMA_PSEARCH_H

#include "vec.h"


/** @brief Function to be minimized
 *  @param pos
 *      Input coordinates, as a coordinate 4-vector
 *  @param data
 *      Callback data
 *  @return The value at @p pos
 *  @note This function should be pure---the algorithm assumes that any repeated
 *      calls with identical coordinates will yield the same results. Failing to
 *      abide this invalidates the result
 */
typedef double gamma_psrch_func_t(const gamma_vec_t *pos, void *data);


/** @brief Pattern search data */
struct gamma_psfunc {
    gamma_psrch_func_t *func;   /* Function to be minimized */
    void               *data;   /* Function data */
    int                 dims;   /* Dimensions (the amount of basis vectors) */
    const gamma_vec_t  *bases;  /* Basis vectors to be stenciled */
};


/** @brief A pair containing a coordinate vector and a function value */
struct gamma_pspair {
    gamma_vec_t vec;    /* Coordinates */
    double      val;    /* Function value */
};


/** @brief Minimize a function by pattern search
 *  @param func
 *      Function to be minimized
 *  @param[in, out] init
 *      Initial value i.e. the initial coordinates and value of @p func and
 *      where the results are written
 *  @param res
 *      Initial stencil resolution
 *  @param shrinks
 *      The maximum number of times the stencil may shrink before the result is
 *      accepted. This parameter may safely be negative. Runtime is at least
 *      linearly dependent on this parameter
 */
void gamma_pattern_search(const struct gamma_psfunc *func,
                          struct gamma_pspair       *init,
                          gamma_scal_t               res,
                          int                        shrinks);


#endif /* GAMMA_PSEARCH_H */
