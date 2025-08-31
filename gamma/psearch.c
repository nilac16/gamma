#include <stdbool.h>
#include <tgmath.h>
#include "psearch.h"


/** @brief Invoke the callback
 *  @param func
 *      Function information
 *  @param pos
 *      Coordinates
 *  @returns The result yielded from the optimizer callback
 */
static struct gamma_pspair
gamma_pattern_invoke(const struct gamma_psfunc *func, const gamma_vec_t *pos)
{
    return (struct gamma_pspair){
        .vec = *pos,
        .val = func->func(pos, func->data)
    };
}


/** @brief Test a point against the current stencil candidate
 *  @param func
 *      Optimizer function
 *  @param cand
 *      Candidate point
 *  @param pos
 *      Test coordinates
 *  @returns true if a new minimum was found, false if not
 */
static bool gamma_pattern_test(const struct gamma_psfunc *func,
                               struct gamma_pspair       *cand,
                               const gamma_vec_t         *pos)
{
    struct gamma_pspair test;
    bool res;

    test = gamma_pattern_invoke(func, pos);
    res = test.val < cand->val;
    if (res) {
        *cand = test;
    }
    return res;
}


void gamma_pattern_search(const struct gamma_psfunc *func,
                          struct gamma_pspair       *init,
                          gamma_scal_t               res,
                          int                        shrinks)
{
    struct gamma_pspair cand;
    gamma_vec_t test;
    bool found;
    int i;

    do {
        cand = *init;
        found = false;
        for (i = 0; i < func->dims; i++) {
            test = gamma_vec_fmadds(&func->bases[i], res, &init->vec);
            found = gamma_pattern_test(func, &cand, &test) || found;
            test = gamma_vec_fmsubs(&func->bases[i], res, &init->vec);
            found = gamma_pattern_test(func, &cand, &test) || found;
        }
        if (found) {
            *init = cand;
        } else {
            res /= 2.0;
            shrinks--;
        }
    } while (shrinks >= 0);
}
