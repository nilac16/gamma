#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <tgmath.h>
#include "gamma.h"
#include "psearch.h"


struct gamma_objective {
    const struct gamma_distribution *ref;       /* Reference dose */
    double                           ratio;     /* Criteria ratio */
    double                           mdose;     /* (Normalized) measured dose */
    double                           rnorm;     /* Reference dose norm */
    gamma_vec_t                      origin;    /* Measured dose origin */
};


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
    struct gamma_objective *obj = data;
    gamma_vec_t diff;
    double dose;

    dose = gamma_distribution_interp(obj->ref, pos);
    diff = gamma_vec_sub(pos, &obj->origin);
    return gamma_sqr(obj->ratio * (dose / obj->rnorm - obj->mdose))
        + gamma_vec_dp(&diff, &diff);
}


/** @brief The full alphabet of parameters */
struct gamma {
    const struct gamma_params       *parms;     /* Gamma parameters */
    const struct gamma_options      *opts;      /* Gamma options */
    const struct gamma_distribution *meas;      /* Measured dose */
    const struct gamma_distribution *ref;       /* Reference dose */
    struct gamma_results            *res;       /* Results */
    double                           mthrsh;    /* Measured dose threshold */
    double                           rthrsh;    /* Reference dose threshold */
    size_t                           idx;       /* Current point index */
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
        .rnorm  = 1.0,
        .origin = *pos,
    };
    struct gamma_psfunc func = {
        .func  = gamma_objective_evaluate,
        .data  = &obj,
        .dims  = BUFLEN(bases),
        .bases = bases
    };
    struct gamma_pspair pair;
    double rdose;

    /* Check if this point is even above threshold */
    rdose = gamma_distribution_interp(gamma->ref, pos);
    if (rdose < gamma->rthrsh && mdose < gamma->mthrsh) {
        return GAMMA_SIG;
    }

    switch (gamma->parms->norm) {
    case GAMMA_NORM_GLOBAL:
    default:
        obj.ratio /= gamma->meas->max;
        break;
    case GAMMA_NORM_LOCAL:
        obj.ratio /= mdose;
        break;
    case GAMMA_NORM_ABSOLUTE:
        obj.ratio /= 1.0;
        break;
    }

    pair.vec = *pos;
    pair.val = gamma_sqr(obj.ratio * (rdose / obj.rnorm - obj.mdose));
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
static bool gamma_iterator(const gamma_vec_t *pos, double dose, void *data)
{
    struct gamma *gamma = data;
    double value;

    value = gamma_pointwise(gamma, pos, dose);
    if (value != GAMMA_SIG) {
        gamma->res->pass += value < 1.0;
        gamma_statistics_add(&gamma->res->stats, value);
    }
    if (gamma->res->dist) {
        gamma->res->dist[gamma->idx] = value;
    }
    gamma->idx++;
    return true;
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
        .idx    = 0,
    };

    res->stats = gamma_statistics_init();
    res->pass = 0;
    gamma_distribution_foreach(meas, gamma_iterator, &gamma);
}
