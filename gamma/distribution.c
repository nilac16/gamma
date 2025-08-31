#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <tgmath.h>
#include "distribution.h"
#include "interp.h"


/** @brief Linearize a multi-index
 *  @param dist
 *      Distribution
 *  @param idx
 *      Multi-index
 *  @returns The index into the main data buffer, directly computed
 *  @warning This function does no bounds nor wraparound (overflow) checks
 */
static uint64_t
gamma_distribution_linearize(const struct gamma_distribution *dist,
                             const gamma_idx_t               *idx)
{
    return (uint64_t)idx->idx[0] + (uint64_t)dist->dims.idx[0]
        * (idx->idx[1] + (uint64_t)dist->dims.idx[1] * idx->idx[2]);
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
    gamma_idx_t test;
    double res = 0.0;

    test = gamma_idx_hittest(idx, &zero, &dist->dims);
    if (!gamma_idx_any(&test)) {
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
    static const gamma_idx_t zero = { 0 };
    int64_t gather[8] = {
        gamma_distribution_linearize(dist, org),
        gather[0] + 1,
        gather[0] + dist->dims.idx[0],
        gather[2] + 1,
        gather[0] + dist->dims.idx[0] * dist->dims.idx[1],
        gather[4] + 1,
        gather[4] + dist->dims.idx[0],
        gather[6] + 1
    };
    gamma_idx_t ext, testlo, testhi;

    ext = gamma_idx_add(org, &(const gamma_idx_t){{ 1, 1, 1, 0 }});
    testlo = gamma_idx_hittest(org, &zero, &dist->dims);
    testhi = gamma_idx_hittest(&ext, &zero, &dist->dims);
    gather[0] |= testlo.idx[0] | testlo.idx[1] | testlo.idx[2];
    gather[1] |= testhi.idx[0] | testlo.idx[1] | testlo.idx[2];
    gather[2] |= testlo.idx[0] | testhi.idx[1] | testlo.idx[2];
    gather[3] |= testhi.idx[0] | testhi.idx[1] | testlo.idx[2];
    gather[4] |= testlo.idx[0] | testlo.idx[1] | testhi.idx[2];
    gather[5] |= testhi.idx[0] | testlo.idx[1] | testhi.idx[2];
    gather[6] |= testlo.idx[0] | testhi.idx[1] | testhi.idx[2];
    gather[7] |= testhi.idx[0] | testhi.idx[1] | testhi.idx[2];

    intr->buf[0] = gather[0] < 0 ? 0.0 : dist->data[gather[0]];
    intr->buf[1] = gather[1] < 0 ? 0.0 : dist->data[gather[1]];
    intr->buf[2] = gather[2] < 0 ? 0.0 : dist->data[gather[2]];
    intr->buf[3] = gather[3] < 0 ? 0.0 : dist->data[gather[3]];
    intr->buf[4] = gather[4] < 0 ? 0.0 : dist->data[gather[4]];
    intr->buf[5] = gather[5] < 0 ? 0.0 : dist->data[gather[5]];
    intr->buf[6] = gather[6] < 0 ? 0.0 : dist->data[gather[6]];
    intr->buf[7] = gather[7] < 0 ? 0.0 : dist->data[gather[7]];
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


void gamma_distribution_foreach(const struct gamma_distribution *dist,
                                gamma_distribution_iterfn_t     *func,
                                void                            *data)
{
    gamma_iscal_t i, j, k;
    gamma_vec_t pos;
    size_t n = 0;

#if _OPENMP
#   pragma omp parallel for private(pos, i, j, n)
#endif
    for (k = 0; k < dist->dims.idx[2]; k++) {
        for (j = 0; j < dist->dims.idx[1]; j++) {
            for (i = 0; i < dist->dims.idx[0]; i++) {
                pos = (const gamma_vec_t){{ i, j, k, 1 }};
                pos = gamma_matmul_mv(&dist->matrix, &pos);
                n = i + dist->dims.idx[0] * (j + dist->dims.idx[1] * k);
                func(&pos, dist->data[n], n, data);
            }
        }
    }
}
