/*
 * nelder-mead.c
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2018-2024 Bryan Franklin. All rights reserved.
 */

#include <float.h>
#include <stdio.h>
#include "../fnt_util.h"
#include "../fnt_vect.h"

/* ---------------------------- */

typedef struct nm_sample {
    fnt_vect_t parameters;
    double value;
} nm_sample;

void nmSampleCopy(nm_sample *dst, nm_sample *src) {
    fnt_vect_copy(&dst->parameters, &src->parameters);
    dst->value = src->value;
}

typedef struct NMSimplex {
    nm_sample *points;
    int count;
} NMSimplex;

void nmSimplexInit(NMSimplex *simplex, int dimensions) {
    memset(simplex, '\0', sizeof(*simplex));
    simplex->points = calloc(dimensions+1, sizeof(nm_sample));
}

void nmSimplexFree(NMSimplex *simplex) {
    for(int i=0; i<simplex->count; ++i) {
        fnt_vect_free(&simplex->points[i].parameters);
    }
    free(simplex->points); simplex->points=NULL;
    simplex->count = 0;
}

void nmSimplexPrint(NMSimplex *simplex) {
    printf("simplex:\n");
    for(int i=0; i<simplex->count; ++i) {
        printf("\tvalue=%g; ", simplex->points[i].value);
        fnt_vect_println(&simplex->points[i].parameters, "parameters: ", NULL);
    }
}

void nmSimplexAdd(NMSimplex *simplex, nm_sample *sample) {
    fnt_vect_calloc(&simplex->points[simplex->count].parameters, sample->parameters.n);
    nmSampleCopy(&simplex->points[simplex->count], sample);
    simplex->count += 1;
}

void nmSimplexSort(NMSimplex *simplex) {
    if( simplex->count <= 1 )
        return;

    int done = 0;
    nm_sample *points = simplex->points;
    nm_sample tmp;
    fnt_vect_calloc(&tmp.parameters, points[0].parameters.n);
    while( !done ) {
        done = 1;
        for(int i=simplex->count-1; i>0; --i) {
            if( points[i-1].value > points[i].value ) {
                /* swap points i and i-1 */
                nmSampleCopy(&tmp, &points[i-1]);
                nmSampleCopy(&points[i-1], &points[i]);
                nmSampleCopy(&points[i], &tmp);

                done = 0;
            }
        }
    }
    fnt_vect_free(&tmp.parameters);
    #if 0
    nmSimplexPrint(simplex);
    #endif /* 0 */
}

/* ---------------------------- */

enum NMState {
    initial, reflect, expand, contract_out, contract_in, shrink, shrink2
};

/* see: http://www.scholarpedia.org/article/Nelder-Mead_algorithm */
typedef struct NelderMead {
    /* maintain search state */
    int dimensions;
    int iterations;
    NMSimplex simplex;
    fnt_vect_t seed;
    enum NMState state;

    /* sample points */
    nm_sample x_r;
    nm_sample x_e;
    nm_sample x_c;
    fnt_vect_t s_shrink;

    /* hyper-parameters */
    double alpha;
    double beta;
    double gamma;
    double delta;
} NelderMead;

void nm_init(void **nm_ptr, int dimensions) {
    NelderMead *nm = calloc(1, sizeof(NelderMead));
    *nm_ptr = nm;
    memset(nm, '\0', sizeof(*nm));

    nm->dimensions = dimensions;
    nm->iterations = 0;
    nm->state = initial;

    nm->alpha = 1;      /* \alpha > 0 */
    nm->beta = 0.5;     /* 0 < \beta < 1 */
    nm->gamma = 2;      /* \gamma > 1 */
    nm->delta = 0.5;    /* 0 < \delta < 1 */

    fnt_vect_calloc(&nm->seed, dimensions);
    fnt_vect_calloc(&nm->x_r.parameters, dimensions);
    fnt_vect_calloc(&nm->x_e.parameters, dimensions);
    fnt_vect_calloc(&nm->x_c.parameters, dimensions);
    fnt_vect_calloc(&nm->s_shrink, dimensions);

    nmSimplexInit(&nm->simplex, dimensions);
}

void nm_free(void *nm_ptr) {
    NelderMead *nm = nm_ptr;

    fnt_vect_free(&nm->seed);
    fnt_vect_free(&nm->x_r.parameters);
    fnt_vect_free(&nm->x_e.parameters);
    fnt_vect_free(&nm->x_c.parameters);
    fnt_vect_free(&nm->s_shrink);

    nmSimplexFree(&nm->simplex);

    free(nm); nm=NULL;
}

void nm_set_seed(void *nm_ptr, fnt_vect_t *seed) {
    NelderMead *nm = nm_ptr;
    if( nm->state != initial ) { return; }
    fnt_vect_copy(&nm->seed, seed);
}

void nm_best_point(void *nm_ptr, fnt_vect_t *result) {
    NelderMead *nm = nm_ptr;

    int best = 0;
    double min = nm->simplex.points[best].value;
    for(int i=0; i<nm->simplex.count; ++i) {
        if( nm->simplex.points[i].value < min ) {
            min = nm->simplex.points[i].value;
            best = i;
        }
    }

    if( best < nm->simplex.count )
        fnt_vect_copy(result, &nm->simplex.points[best].parameters);
}

void nm_add_result(void *nm_ptr, fnt_vect_t *parameters, double value) {
    NelderMead *nm = nm_ptr;

    nm->iterations += 1;

    nm_sample newSample;
    fnt_vect_calloc(&newSample.parameters, parameters->n);
    fnt_vect_copy(&newSample.parameters, parameters);
    newSample.value = value;

    /* shrink just replaces points, but needs the associated values */
    if( nm->state == shrink2 ) {
        nmSampleCopy(&nm->simplex.points[nm->simplex.count-2], &newSample);
        nm->state = reflect;
        fnt_vect_free(&newSample.parameters);
        return;
    } else if( nm->state == shrink ) {
        nmSampleCopy(&nm->simplex.points[nm->simplex.count-1], &newSample);
        nm->state = shrink2;
        fnt_vect_free(&newSample.parameters);
        return;
    }

    /* check for initialization state */
    if( nm->simplex.count <= nm->dimensions ) {
        nmSimplexAdd(&nm->simplex, &newSample);
        if( nm->simplex.count >= nm->dimensions+1 )
            nm->state = reflect;
        fnt_vect_free(&newSample.parameters);
        return;
    }

    /* sort simple */
    if( nm->state != shrink && nm->state != shrink2 )
        nmSimplexSort(&nm->simplex);

    /* get h, s, l */
    nm_sample h; /* simplex[simplex.count-1] */
    nm_sample s; /* simplex[simplex.count-2] */
    nm_sample l; /* simplex[0] */
    nm_sample r; /* newSample */
    fnt_vect_calloc(&h.parameters, nm->dimensions);
    fnt_vect_calloc(&s.parameters, nm->dimensions);
    fnt_vect_calloc(&l.parameters, nm->dimensions);
    fnt_vect_calloc(&r.parameters, nm->dimensions);
    nmSampleCopy(&h, &nm->simplex.points[nm->simplex.count-1]);
    nmSampleCopy(&s, &nm->simplex.points[nm->simplex.count-2]);
    nmSampleCopy(&l, &nm->simplex.points[0]);
    nmSampleCopy(&r, &newSample);
    fnt_vect_free(&newSample.parameters);

    /* actual parameters can be discarded,
     * as only the output values are needed here */
    fnt_vect_free(&h.parameters);
    fnt_vect_free(&s.parameters);
    fnt_vect_free(&l.parameters);

    #if 0
    printf("f(h) = %g\n", h.value);
    printf("f(s) = %g\n", s.value);
    printf("f(l) = %g\n", l.value);
    printf("f(r) = %g\n", r.value);
    #endif /* 0 */

    /* deal with recently computed point based on state */
    if( nm->state == reflect ) {
        nmSampleCopy(&nm->x_r, &r);

        if( l.value <= nm->x_r.value && nm->x_r.value < s.value ) {
            /* accept x_r and terminate iteration */
            nmSampleCopy(&nm->simplex.points[nm->simplex.count-1], &r);
            fnt_vect_free(&r.parameters);
            return;
        }
    }

    if( nm->state == expand ) {
        nmSampleCopy(&nm->x_e, &r);

        if( nm->x_e.value < nm->x_r.value ) {
            /* accept x_e and terminate iteration */
            nmSampleCopy(&nm->simplex.points[nm->simplex.count-1], &nm->x_e);
        } else {
            /* accept x_r and terminate iteration */
            nmSampleCopy(&nm->simplex.points[nm->simplex.count-1], &nm->x_r);
        }

        nm->state = reflect;
        fnt_vect_free(&r.parameters);
        return;
    }

    if( nm->state == contract_out ) {
        nmSampleCopy(&nm->x_c, &r);

        if( nm->x_c.value < nm->x_r.value ) {
            /* accept x_c and terminate iteration */
            nmSampleCopy(&nm->simplex.points[nm->simplex.count-1], &nm->x_c);
            nm->state = reflect;
            fnt_vect_free(&r.parameters);
            return;
        }
    }

    if( nm->state == contract_in ) {
        nmSampleCopy(&nm->x_c, &r);

        if( nm->x_c.value < h.value ) {
            /* accept x_c and terminate iteration */
            nmSampleCopy(&nm->simplex.points[nm->simplex.count-1], &nm->x_c);
            nm->state = reflect;
            fnt_vect_free(&r.parameters);
            return;
        }
    }
    fnt_vect_free(&r.parameters);

    /* determine next state if new point not accepted */
    if( r.value < l.value ) {
        /* cause x_e to be computed next */
        nm->state = expand;
        return;
    } else if( r.value >= s.value ) {
        /* cause x_c to be computed next */
        if( s.value <= r.value && r.value < h.value )
            nm->state = contract_out;
        else
            nm->state = contract_in;
        return;
    }

    nm->state = shrink;
    return;
}

void nm_next_point(void *nm_ptr, fnt_vect_t *vector) {
    NelderMead *nm = nm_ptr;

    if( nm->state == initial && nm->simplex.count < nm->dimensions+1 ) {
        /* add new initial point */
        if( nm->simplex.count > 0 ) {
            int pos = nm->simplex.count-1;
            fnt_vect_copy(vector, &nm->seed);
            vector->v[pos] += nm->simplex.count;
        } else if( nm->seed.n == vector->n ) {
            fnt_vect_copy(vector, &nm->seed);
        } else {
            fnt_vect_copy(&nm->seed, vector);
        }
        return;
    }

    /* check to see that we have correct number of points */
    if( nm->simplex.count != nm->dimensions+1 ) {
        fnt_vect_copy(vector, &nm->seed);
        return;
    }

    /* sort simplex */
    if( nm->state != shrink && nm->state != shrink2 )
        nmSimplexSort(&nm->simplex);

    /* get h, s, l */
    nm_sample h; /* simplex[simplex.count-1] */
    nm_sample s; /* simplex[simplex.count-2] */
    nm_sample l; /* simplex[0] */
    fnt_vect_calloc(&h.parameters, nm->dimensions);
    fnt_vect_calloc(&s.parameters, nm->dimensions);
    fnt_vect_calloc(&l.parameters, nm->dimensions);
    nmSampleCopy(&h, &nm->simplex.points[nm->simplex.count-1]);
    nmSampleCopy(&s, &nm->simplex.points[nm->simplex.count-2]);
    nmSampleCopy(&l, &nm->simplex.points[0]);

    /* compute centroid */
    fnt_vect_t c, sumVect;
    fnt_vect_calloc(&c, nm->dimensions);
    fnt_vect_calloc(&sumVect, nm->dimensions);
    for(int i=0; i<nm->simplex.count-1; ++i) {
        fnt_vect_add(&sumVect, &nm->simplex.points[i].parameters, &sumVect);
    }
    fnt_vect_scale(&sumVect, 1.0/(nm->simplex.count-1), &c);
    fnt_vect_free(&sumVect);

    /* add appropriate new point(s) */
    fnt_vect_t scaled, tmp;
    fnt_vect_calloc(&scaled, nm->dimensions);
    fnt_vect_calloc(&tmp, nm->dimensions);
    switch( nm->state ) {
        case initial:
            /* initial is handled above */
            printf("This should never happen!\n");
            break;
        case reflect:
            fnt_vect_sub(&c, &h.parameters, &tmp);
            fnt_vect_scale(&tmp, nm->alpha, &scaled);
            fnt_vect_add(&c, &scaled, vector);    /* x_r */
            break;
        case expand:
            fnt_vect_sub(&nm->x_r.parameters, &c, &tmp);
            fnt_vect_scale(&tmp, nm->gamma, &scaled);
            fnt_vect_add(&c, &scaled, vector);    /* x_e */
            break;
        case contract_out:
            fnt_vect_sub(&nm->x_r.parameters, &c, &tmp);
            fnt_vect_scale(&tmp, nm->beta, &scaled);
            fnt_vect_add(&c, &scaled, vector);    /* x_c */
            break;
        case contract_in:
            fnt_vect_sub(&h.parameters, &c, &tmp);
            fnt_vect_scale(&tmp, nm->beta, &scaled);
            fnt_vect_add(&c, &scaled, vector);    /* x_c */
            break;
        case shrink:
            /* store one shrink point for state shrink2 */
            fnt_vect_add(&nm->x_r.parameters, &s.parameters, &tmp);
            fnt_vect_scale(&tmp, 0.5, &nm->s_shrink); /* new s */

            /* return other shrink point */
            fnt_vect_add(&nm->x_r.parameters, &h.parameters, &tmp);
            fnt_vect_scale(&tmp, 0.5, vector);        /* new h */
            break;
        case shrink2:
            /* return second shrink point */
            fnt_vect_copy(vector, &nm->s_shrink);
            fnt_vect_reset(&nm->s_shrink);
            break;
    }
    fnt_vect_free(&tmp);
    fnt_vect_free(&scaled);
    fnt_vect_free(&c);
    fnt_vect_free(&h.parameters);
    fnt_vect_free(&s.parameters);
    fnt_vect_free(&l.parameters);

    return;
}

int nm_simplex_point(void *nm_ptr, int which, fnt_vect_t *point, double *value) {
    NelderMead *nm = nm_ptr;
    if( which >= nm->simplex.count )
        return 0;

    if( point )
        fnt_vect_copy(point, &nm->simplex.points[which].parameters);
    if( value )
        *value = nm->simplex.points[which].value;

    return 1;
}

int nm_done(void *nm_ptr, double threshold, int iterations) {
    NelderMead *nm = nm_ptr;
    if( nm->state == initial )
        return 0;

    /* check maximum iterations */
    if( nm->iterations > iterations ) {
        return 1;
    }

    /* check distance between best and worst parameters */
    double dist;
    if( nm->state != shrink && nm->state != shrink2 )
        nmSimplexSort(&nm->simplex);
    fnt_vect_dist(&nm->simplex.points[0].parameters,
                &nm->simplex.points[nm->simplex.count-1].parameters, &dist);
    #if 0
    printf("iteration: %i; dist: %g\n", nm->iterations, dist);
    #endif /* 0 */
    if( dist < threshold ) {
        return 1;
    }

    return 0;
}
