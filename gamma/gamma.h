#pragma once

#ifndef GAMMA_H
#define GAMMA_H

#include <stdbool.h>
#include "common.h"
#include "distribution.h"
#include "statistics.h"

EXTERN_C_BEGIN


/** @brief Signal value, used for points below threshold when calculating a
 *      distribution
 */
#define GAMMA_SIG (-1.0)


/** @brief % dose difference normalization modes */
typedef enum gamma_normalization {
    GAMMA_NORM_GLOBAL,      /* Use the maximum measured dose value */
    GAMMA_NORM_LOCAL,       /* Use the current measured dose value */
    GAMMA_NORM_ABSOLUTE,    /* Use the supplied dose value */
} gamma_norm_t;


/** @brief Primary gamma parameters */
struct gamma_params {
    double       diff;  /* %difference criterion as a proportion (e.g. 0.03) */
    double       dta;   /* Distance-to-agreement in units of pixel spacing */
    double       thrsh; /* Low-dose threshold as a proportion (e.g. 0.10) */
    gamma_norm_t norm;  /* Normalization mode */
    bool         rel;   /* If true then normalize both distributions */
};


/** @brief Extra options not traditionally considered gamma parameters */
struct gamma_options {
    bool pass_only;     /* Terminate immediately upon finding a pass */
    long shrinks;       /* Pattern search stencil shrink limit */
};


/** @brief Results buffer. Only the pointer must be set by you, and if it is, it
 *      must address a buffer at least the size of the measured dose buffer
 */
struct gamma_results {
    struct gamma_statistics stats;      /* Point statistics */
    long                    pass;       /* Total passing points */
    double                 *dist;       /* The gamma distribution, if nonnull */
};


/** @brief Compute gamma index statistics for two dose distributions
 *  @param params
 *      Gamma parameters
 *  @param options
 *      Extra gamma options
 *  @param ref
 *      Reference/baseline distribution
 *  @param meas
 *      Test distribution
 *  @param[out] res
 *      Results buffer
 */
void gamma_compute(const struct gamma_params       *params,
                   const struct gamma_options      *options,
                   const struct gamma_distribution *ref,
                   const struct gamma_distribution *meas,
                   struct gamma_results            *res);


EXTERN_C_END

#endif /* GAMMA_H */
