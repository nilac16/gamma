#include <stdbool.h>
#include "psearch.h"


void gamma_pattern_search(const struct gamma_psfunc *func,
                          struct gamma_pspair       *init,
                          gamma_scal_t               res,
                          int                        shrinks)
{
    struct gamma_pspair cand, test;
    bool found;
    int i;

    shrinks = shrinks < 0 ? 0 : shrinks;
    do {
        found = false;
        cand = *init;
        for (i = 0; i < func->dims; i++) {
            test.vec = gamma_vec_fmadds(&func->bases[i], res, &init->vec);
            test.val = func->func(&test.vec, func->data);
            if (test.val < cand.val) {
                found = true;
                cand = test;
            }
            test.vec = gamma_vec_fmsubs(&func->bases[i], res, &init->vec);
            test.val = func->func(&test.vec, func->data);
            if (test.val < cand.val) {
                found = true;
                cand = test;
            }
        }
        if (found) {
            *init = cand;
        } else {
            res /= 2.0;
            shrinks--;
        }
    } while (shrinks >= 0);
}
