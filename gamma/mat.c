#include <tgmath.h>
#include "mat.h"


const gamma_mat_t gamma_mat_identity = {
    (const gamma_vec_t){{ 1, 0, 0, 0 }},
    (const gamma_vec_t){{ 0, 1, 0, 0 }},
    (const gamma_vec_t){{ 0, 0, 1, 0 }},
    (const gamma_vec_t){{ 0, 0, 0, 1 }},
};


/** @brief Swap two vectors
 *  @param a
 *      Vector
 *  @param b
 *      Vector
 */
static void gamma_vecswap(gamma_vec_t *a, gamma_vec_t *b)
{
    gamma_vec_t tmp = *a;

    *a = *b;
    *b = tmp;
}


/** @brief Find a pivot column
 *  @param mat
 *      Matrix
 *  @param idx
 *      The column to be pivot
 *  @returns The index of the pivot column or -1 if none was found
 */
static int gamma_mat_pivot(const gamma_mat_t *mat, const int idx)
{
    gamma_scal_t max, test;
    int piv = idx, i;

    max = fabs(mat->cols[idx].vec[idx]);
    for (i = idx + 1; i < 4; i++) {
        test = fabs(mat->cols[i].vec[idx]);
        if (isgreater(test, max)) {
            piv = i;
            max = test;
        }
    }
    return max ? piv : -1;
}


/** @brief Perform two-way column elimination
 *  @param lhs
 *      Matrix to be inverted
 *  @param rhs
 *      Output
 *  @param idx
 *      Column index for Gaussian elimination
 *  @returns true on success, false if no pivot could be found for @p lhs (i.e.
 *      it is singular)
 */
static bool gamma_mat_fwdelim(gamma_mat_t *lhs, gamma_mat_t *rhs, const int idx)
{
    gamma_scal_t norm, mult;
    int piv, i;

    /* Find a pivot */
    piv = gamma_mat_pivot(lhs, idx);
    if (piv < 0) {
        return false;
    }
    gamma_vecswap(&lhs->cols[idx], &lhs->cols[piv]);
    gamma_vecswap(&rhs->cols[idx], &rhs->cols[piv]);

    /* Normalize the current column */
    norm = lhs->cols[idx].vec[idx];
    lhs->cols[idx] = gamma_vec_divs(&lhs->cols[idx], norm);
    rhs->cols[idx] = gamma_vec_divs(&rhs->cols[idx], norm);

    /* Reduce _all_ remaining columns */
    for (i = 0; i < 4; i++) {
        if (i == idx) {
            continue;
        }
        mult = lhs->cols[i].vec[idx];
        lhs->cols[i] = gamma_vec_fmsubs(&lhs->cols[idx], mult, &lhs->cols[i]);
        rhs->cols[i] = gamma_vec_fmsubs(&rhs->cols[idx], mult, &rhs->cols[i]);
    }
    return true;
}


bool gamma_mat_ldivide(const gamma_mat_t *lhs, gamma_mat_t *rhs)
{
    gamma_mat_t left = *lhs;

    return gamma_mat_fwdelim(&left, rhs, 0)
        && gamma_mat_fwdelim(&left, rhs, 1)
        && gamma_mat_fwdelim(&left, rhs, 2)
        && gamma_mat_fwdelim(&left, rhs, 3);
}


bool gamma_mat_invert(gamma_mat_t *mat)
{
    gamma_mat_t left = *mat;

    *mat = gamma_mat_identity;
    return gamma_mat_ldivide(&left, mat);
}
