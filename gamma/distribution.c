#include <assert.h>
#include <stdio.h>
#include <tgmath.h>
#include "distribution.h"
#include "interp.h"


static size_t gamma_distribution_linearize(const struct gamma_distribution *dis,
                                           const gamma_idx_t               *idx)
{
    return (size_t)idx->idx[0] + dis->dims.idx[0]
        * (idx->idx[1] + dis->dims.idx[1] * idx->idx[2]);
}


bool gamma_distribution_set(struct gamma_distribution *dist,
                            const gamma_mat_t         *matr,
                            const gamma_idx_t         *dims,
                            double                    *data)
{
    size_t i;

    dist->matrix = *matr;
    dist->inverse = *matr;
    if (!gamma_mat_invert(&dist->inverse)) {
        return false;
    }
    dist->dims = *dims;
    dist->dims.idx[3] = INT32_MAX;
    dist->len = (size_t)dims->idx[0] * dims->idx[1] * dims->idx[2];
    dist->max = -HUGE_VAL;
    for (i = 0; i < dist->len; i++) {
        dist->max = fmax(dist->max, data[i]);
    }
    dist->data = data;
    return true;
}


double gamma_distribution_at(const struct gamma_distribution *dist,
                             const gamma_idx_t               *idx)
{
    const gamma_idx_t zero = { 0 };
    double res = 0.0;

    if (gamma_idx_hittest(idx, &zero, &dist->dims)) {
        res = dist->data[gamma_distribution_linearize(dist, idx)];
    }
    return res;
}


/** @brief Split a vector into integral and fractional parts
 *  @param[in, out] vec
 *      The vector on input and on output to receive the fractional parts
 *  @param[out] idx
 *      The index to receive the integral components, rounded towards zero
 */
static void gamma_distribution_modf(gamma_vec_t *vec, gamma_idx_t *idx)
{
    idx->idx[0] = vec->vec[0];
    idx->idx[1] = vec->vec[1];
    idx->idx[2] = vec->vec[2];
    idx->idx[3] = vec->vec[3];
    vec->vec[0] -= idx->idx[0];
    vec->vec[1] -= idx->idx[1];
    vec->vec[2] -= idx->idx[2];
    vec->vec[3] -= idx->idx[3];
}


static void gamma_distribution_corners(const struct gamma_distribution *dist,
                                       struct gamma_interp             *intr,
                                       const gamma_idx_t               *org)
{
    const gamma_idx_t xoffs = {{ 1, 0, 0, 0 }};
    const gamma_idx_t yoffs = {{ 0, 1, 0, 0 }};
    const gamma_idx_t xyoffs = {{ 1, 1, 0, 0 }};
    gamma_idx_t idx, up;

    up = gamma_idx_add(org, &(const gamma_idx_t){{ 0, 0, 1, 0 }});

    idx = *org;
    intr->buf[0] = gamma_distribution_at(dist, &idx);

    idx = gamma_idx_add(org, &xoffs);
    intr->buf[1] = gamma_distribution_at(dist, &idx);

    idx = gamma_idx_add(org, &yoffs);
    intr->buf[2] = gamma_distribution_at(dist, &idx);

    idx = gamma_idx_add(org, &xyoffs);
    intr->buf[3] = gamma_distribution_at(dist, &idx);

    idx = up;
    intr->buf[4] = gamma_distribution_at(dist, &idx);

    idx = gamma_idx_add(&up, &xoffs);
    intr->buf[5] = gamma_distribution_at(dist, &idx);

    idx = gamma_idx_add(&up, &yoffs);
    intr->buf[6] = gamma_distribution_at(dist, &idx);

    idx = gamma_idx_add(&up, &xyoffs);
    intr->buf[7] = gamma_distribution_at(dist, &idx);
}


double gamma_distribution_interp(const struct gamma_distribution *dist,
                                 const gamma_vec_t               *pos)
{
    struct gamma_interp interp;
    gamma_vec_t offs;
    gamma_idx_t lat;

    offs = gamma_matmul_mv(&dist->inverse, pos);
    gamma_distribution_modf(&offs, &lat);
    gamma_distribution_corners(dist, &interp, &lat);
    return gamma_interp_single(&interp, &offs);
}


bool gamma_distribution_foreach(const struct gamma_distribution *dist,
                                gamma_distribution_iterfn_t     *func,
                                void                            *data)
{
    gamma_iscal_t i, j, k;
    gamma_vec_t pos;
    size_t n = 0;

    for (k = 0; k < dist->dims.idx[2]; k++) {
        for (j = 0; j < dist->dims.idx[1]; j++) {
            for (i = 0; i < dist->dims.idx[0]; i++) {
                pos = gamma_matmul_mv(&dist->matrix,
                                      &(const gamma_vec_t){{ i, j, k, 1 }});
                if (!func(&pos, dist->data[n++], data)) {
                    return false;
                }
            }
        }
    }
    return true;
}
