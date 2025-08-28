#pragma once

#ifndef GAMMA_STATISTICS_H
#define GAMMA_STATISTICS_H

#include "common.h"
#include <tgmath.h>

EXTERN_C_BEGIN


struct gamma_statistics {
    long   total;   /* Total points */
    double min;     /* Minimum value */
    double max;     /* Maximum value */
    double mean;    /* Arithmetic mean */
    double msqr;    /* Arithmetic mean of squares */
};


#if !defined(__cplusplus) || !__cplusplus
#define gamma_statistics_init() (struct gamma_statistics){ \
    .min = HUGE_VAL, \
    .max = -HUGE_VAL, \
}

#else
GAMMA_INLINE struct gamma_statistics gamma_statistics_init(void)
{
    struct gamma_statistics res{ };

    res.min = HUGE_VAL;
    res.max = -HUGE_VAL;
    return res;
}

#endif


/** @brief Add a new datapoint
 *  @param stat
 *      Statistics
 *  @param x
 *      Data value
 */
GAMMA_INLINE void gamma_statistics_add(struct gamma_statistics *stat, double x)
{
    const double orig_len = stat->total;

    stat->total++;
    stat->min  = fmin(stat->min, x);
    stat->max  = fmax(stat->max, x);
    stat->mean = (orig_len * stat->mean + x) / (double)stat->total;
    stat->msqr = (orig_len * stat->msqr + gamma_sqr(x)) / (double)stat->total;
}


EXTERN_C_END

#endif /* GAMMA_STATISTICS_H */
