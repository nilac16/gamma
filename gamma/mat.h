#pragma once

#ifndef GAMMA_MATH_MAT_H
#define GAMMA_MATH_MAT_H

#include <stdbool.h>
#include "common.h"
#include "vec.h"

EXTERN_C_BEGIN


/** @brief A matrix of 4 4-vectors in column-major orientation */
typedef struct gamma_mat {
    gamma_vec_t cols[4];
} gamma_mat_t;


/** @brief The 4x4 identity matrix */
extern const gamma_mat_t gamma_mat_identity;


/** @brief Matrix left-multiplication with a vector
 *  @param mat
 *      Matrix left operand
 *  @param[in, out] vec
 *      Vector right operand and destination
 */
GAMMA_INLINE gamma_vec_t gamma_matmul_mv(const gamma_mat_t *mat,
                                         const gamma_vec_t *vec)
{
    gamma_vec_t res = { 0 };

    res = gamma_vec_fmadds(&mat->cols[0], vec->vec[0], &res);
    res = gamma_vec_fmadds(&mat->cols[1], vec->vec[1], &res);
    res = gamma_vec_fmadds(&mat->cols[2], vec->vec[2], &res);
    res = gamma_vec_fmadds(&mat->cols[3], vec->vec[3], &res);
    return res;
}


/** @brief Matrix-matrix multiplication
 *  @param lhs
 *      Left operand
 *  @param[in, out] rhs
 *      Right operand and destination
 */
GAMMA_INLINE void gamma_matmul_mm(const gamma_mat_t *lhs, gamma_mat_t *rhs)
{
    rhs->cols[0] = gamma_matmul_mv(lhs, &rhs->cols[0]);
    rhs->cols[1] = gamma_matmul_mv(lhs, &rhs->cols[1]);
    rhs->cols[2] = gamma_matmul_mv(lhs, &rhs->cols[2]);
    rhs->cols[3] = gamma_matmul_mv(lhs, &rhs->cols[3]);
}


/** @brief Divide two matrices
 *  @param lhs
 *      Left operand, to be inverted and left-multiplied
 *  @param[in, out] rhs
 *      Right operand and destination
 *  @returns true on success, false if @p lhs is singular. If @p lhs is found to
 *      be singular, then the contents of @p rhs are undefined
 *  @note This function has similar semantics as its MATLAB equivalent: It is
 *      faster to do this once than to invert @p lhs and then multiply it with
 *      @p rhs once
 */
bool gamma_mat_ldivide(const gamma_mat_t *lhs, gamma_mat_t *rhs);


/** @brief Invert a matrix
 *  @param mat
 *      Matrix
 *  @returns true on success, false if @p mat is singular (i.e. noninvertible).
 *      If @p mat is singular, then its contents are undefined upon return from
 *      this function
 */
bool gamma_mat_invert(gamma_mat_t *mat);


EXTERN_C_END

#endif /* GAMMA_MATH_MAT_H */
