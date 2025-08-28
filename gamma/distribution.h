#pragma once

#ifndef GAMMA_DISTRIBUTION_H
#define GAMMA_DISTRIBUTION_H

#include <stddef.h>
#include "common.h"
#include "mat.h"
#include "idx.h"

EXTERN_C_BEGIN


/** @brief All of the information needed for a dose distribution embedded in R3
 */
struct gamma_distribution {
    gamma_mat_t matrix;     /* Pixel-to-physical affine transformation */
    gamma_mat_t inverse;    /* Physical-to-pixelspace inverse transform */
    gamma_idx_t dims;       /* Pixel dimensions */
    size_t      len;        /* Pixel count */
    double      max;        /* Maximum pixel value */
    double     *data;       /* Pixel data */
};


/** @brief Initialize the distribution
 *  @param dist
 *      Distribution. This should not be managing any memory (free the data
 *      pointer before this call if you need to)
 *  @param matr
 *      Affine matrix
 *  @param dims
 *      Pixel dimensions
 *  @param data
 *      Pixel data
 *  @returns true on success, false if @p matr is singular
 *  @note Each pointer may address the relevant member in @p dist without any
 *      pernicious effects
 */
bool gamma_distribution_set(struct gamma_distribution *dist,
                            const gamma_mat_t         *matr,
                            const gamma_idx_t         *dims,
                            double                    *data);


/** @brief Get a value
 *  @param dist
 *      Distribution
 *  @param idx
 *      Index to fetch. If this is out-of-bounds then zero is returned
 *  @returns The value at @p idx without fail
 */
double gamma_distribution_at(const struct gamma_distribution *dist,
                             const gamma_idx_t               *idx);


/** @brief Interpolate a value
 *  @param dist
 *      Dose distribution
 *  @param pos
 *      Real-valued physical coordinates to interpolate data from
 *  @returns The dose value at @p pos or zero if it was out of bounds. This
 *      function does not fail
 */
double gamma_distribution_interp(const struct gamma_distribution *dist,
                                 const gamma_vec_t               *pos);


/** @brief Iterator callback
 *  @param pos
 *      Physical coordinates of this dose value
 *  @param dose
 *      This dose value
 *  @param data
 *      Your callback data
 *  @return true to continue, false to stop iterating
 */
typedef bool gamma_distribution_iterfn_t(const gamma_vec_t *pos,
                                         double             dose,
                                         void              *data);


/** @brief Iterate over a distribution
 *  @param dist
 *      Distribution
 *  @param func
 *      Iterator function
 *  @param data
 *      Iterator function data
 *  @returns true on completion, false if the callback requested early
 *      termination
 */
bool gamma_distribution_foreach(const struct gamma_distribution *dist,
                                gamma_distribution_iterfn_t     *func,
                                void                            *data);


EXTERN_C_END

#endif /* GAMMA_DISTRIBUTION_H */
