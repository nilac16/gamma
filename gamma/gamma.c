#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <tgmath.h>
#include "gamma.h"
#include "psearch.h"

#if _OPENMP
#   include <threads.h>
#   define ADD_MUTEX(name)  mtx_t name
#   define MTX_INIT(mtx)    mtx_init(mtx, mtx_plain)
#   define MTX_DESTROY(mtx) mtx_destroy(mtx)
#   define MTX_LOCK(mtx)    mtx_lock(mtx)
#   define MTX_UNLOCK(mtx)  mtx_unlock(mtx)

#else
#   define ADD_MUTEX(name)
#   define MTX_INIT(mtx)
#   define MTX_DESTROY(mtx)
#   define MTX_LOCK(mtx)
#   define MTX_UNLOCK(mtx)

#endif


struct gamma_objective {
    const struct gamma_distribution *ref;       /* Reference dose */
    double                           ratio;     /* Criteria ratio */
    double                           mdose;     /* (Normalized) measured dose */
    gamma_vec_t                      origin;    /* Measured dose origin */
};


/** @brief Consistently evaluate the objective function given inputs
 *  @param obj
 *      Objective function
 *  @param rdose
 *      Reference dose value
 *  @param rdiff
 *      Displacement vector from the measured dose point (not a coordinate
 *      vector---be certain there is not a one in the last position!)
 *  @returns The value of the objective function for this dose, coordinate tuple
 */
static double gamma_objective_value(const struct gamma_objective *obj,
                                    double                        rdose,
                                    const gamma_vec_t            *rdiff)
{
    return gamma_sqr(obj->ratio * (rdose - obj->mdose))
        + gamma_vec_dp(rdiff, rdiff);
}


/** @brief Evaluate the objective function at a given point
 *  @note Optimizer callback
 *  @param pos
 *      Coordinates
 *  @param data
 *      Objective
 *  @returns The value of the distance-gamma objective
 */
static double gamma_objective_evaluate(const gamma_vec_t *pos, void *data)
{
    const struct gamma_objective *obj = data;
    gamma_vec_t diff;

    diff = gamma_vec_sub(pos, &obj->origin);
    return gamma_objective_value(obj,
                                 gamma_distribution_interp(obj->ref, pos),
                                 &diff);
}


/** @brief The full alphabet of parameters */
struct gamma {
    const struct gamma_params       *parms;     /* Gamma parameters */
    const struct gamma_options      *opts;      /* Gamma options */
    const struct gamma_distribution *ref;       /* Reference dose */
    const struct gamma_distribution *meas;      /* Measured dose */
    struct gamma_results            *res;       /* Results */
    double                           rthrsh;    /* Reference dose threshold */
    double                           mthrsh;    /* Measured dose threshold */
    ADD_MUTEX(mtx);
};


/** @brief Do pointwise gamma
 *  @param gamma
 *      Gamma context
 *  @param pos
 *      Measured dose physical coordinates
 *  @param mdose
 *      Measured dose value
 *  @returns The gamma value at this point
 */
static double gamma_pointwise(const struct gamma *gamma,
                              const gamma_vec_t  *pos,
                              double              mdose)
{
    const gamma_vec_t bases[] = {
        {{ 1, 0, 0, 0 }},
        {{ 0, 1, 0, 0 }},
        {{ 0, 0, 1, 0 }},
    };
    struct gamma_objective obj = {
        .ref    = gamma->ref,
        .ratio  = gamma->parms->dta / gamma->parms->diff,
        .mdose  = mdose,
        .origin = *pos,
    };
    const struct gamma_psfunc func = {
        .func  = gamma_objective_evaluate,
        .data  = &obj,
        .dims  = BUFLEN(bases),
        .bases = bases
    };
    struct gamma_pspair pair;
    double dnorm, rdose;

    /* Check if this point is even above threshold */
    rdose = gamma_distribution_interp(gamma->ref, pos);
    if (rdose < gamma->rthrsh && mdose < gamma->mthrsh) {
        return GAMMA_SIG;
    }

    switch (gamma->parms->norm) {
    case GAMMA_NORM_GLOBAL:
    default:
        dnorm = gamma->meas->max;
        break;
    case GAMMA_NORM_LOCAL:
        dnorm = mdose;
        break;
    case GAMMA_NORM_ABSOLUTE:
        dnorm = 1.0;
        break;
    }
    obj.ratio /= dnorm;

    if (gamma->parms->rel) {
        /* Normalize the measured dose value to the reference dose's range */
        obj.mdose *= gamma->ref->max / gamma->meas->max;
    }

    pair.vec = *pos;
    pair.val = gamma_objective_value(&obj, rdose, &(const gamma_vec_t){ 0 });
    gamma_pattern_search(&func, &pair, gamma->parms->dta, gamma->opts->shrinks);
    return sqrt(pair.val) / gamma->parms->dta;
}


/** @brief Iterator callback
 *  @param pos
 *      Physical coordinates of this position in the measured dose
 *  @param dose
 *      The dose value
 *  @param data
 *      The gamma context pointer
 *  @returns true
 */
static void gamma_iterator(const gamma_vec_t *pos,
                           double             dose,
                           size_t             idx,
                           void              *data)
{
    struct gamma *gamma = data;
    double value;

    value = gamma_pointwise(gamma, pos, dose);
    if (value != GAMMA_SIG) {
        MTX_LOCK(&gamma->mtx);
        gamma->res->pass += value < 1.0;
        gamma_statistics_add(&gamma->res->stats, value);
        MTX_UNLOCK(&gamma->mtx);
    }
    if (gamma->res->dist) {
        gamma->res->dist[idx] = value;
    }
}


void gamma_compute(const struct gamma_params       *params,
                   const struct gamma_options      *options,
                   const struct gamma_distribution *ref,
                   const struct gamma_distribution *meas,
                   struct gamma_results            *res)
{
    struct gamma gamma = {
        .parms  = params,
        .opts   = options,
        .ref    = ref,
        .meas   = meas,
        .res    = res,
        .rthrsh = params->thrsh * ref->max,
        .mthrsh = params->thrsh * meas->max,
    };

    MTX_INIT(&gamma.mtx);

    res->stats = gamma_statistics_init();
    res->pass = 0;
    gamma_distribution_foreach(meas, gamma_iterator, &gamma);

    MTX_DESTROY(&gamma.mtx);
}
