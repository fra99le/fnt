/*
 * nelder-mead.c
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2018-2024 Bryan Franklin. All rights reserved.
 */

#include <float.h>
#include <stdio.h>
#include "../fnt.h"
#include "../fnt_util.h"
#include "../fnt_vect.h"

extern int fnt_verbose_level;

/* MARK: Sample types */

typedef struct nm_sample {
    fnt_vect_t parameters;
    double value;
} nm_sample_t;


void nm_sample_copy(nm_sample_t *dst, nm_sample_t *src) {
    fnt_vect_copy(&dst->parameters, &src->parameters);
    dst->value = src->value;
}


/* MARK: Simplex types */

typedef struct nm_simplex {
    nm_sample_t *points;
    int count;
} nm_simplex_t;


/* MARK: Simplex functions */

void nm_simplex_init(nm_simplex_t *simplex, int dimensions) {
    memset(simplex, '\0', sizeof(*simplex));
    simplex->points = calloc(dimensions+1, sizeof(nm_sample_t));
}


void nm_simplex_free(nm_simplex_t *simplex) {
    for(int i=0; i<simplex->count; ++i) {
        fnt_vect_free(&simplex->points[i].parameters);
    }
    free(simplex->points); simplex->points=NULL;
    simplex->count = 0;
}


void nm_simplex_print(nm_simplex_t *simplex) {
    printf("simplex:\n");
    for(int i=0; i<simplex->count; ++i) {
        printf("\tvalue=%g; ", simplex->points[i].value);
        fnt_vect_println(&simplex->points[i].parameters, "parameters: ", NULL);
    }
}


void nm_simplex_add(nm_simplex_t *simplex, nm_sample_t *sample) {
    fnt_vect_calloc(&simplex->points[simplex->count].parameters, sample->parameters.n);
    nm_sample_copy(&simplex->points[simplex->count], sample);
    simplex->count += 1;
}


void nm_simplex_sort(nm_simplex_t *simplex) {
    if( simplex->count <= 1 )
        return;

    int done = 0;
    nm_sample_t *points = simplex->points;
    nm_sample_t tmp;
    fnt_vect_calloc(&tmp.parameters, points[0].parameters.n);
    while( !done ) {
        done = 1;
        for(int i=simplex->count-1; i>0; --i) {
            if( points[i-1].value > points[i].value ) {
                /* swap points i and i-1 */
                nm_sample_copy(&tmp, &points[i-1]);
                nm_sample_copy(&points[i-1], &points[i]);
                nm_sample_copy(&points[i], &tmp);

                done = 0;
            }
        }
    }
    fnt_vect_free(&tmp.parameters);

    if( fnt_verbose_level >= 3 ) {
        nm_simplex_print(simplex);
    }
}


/* MARK: Nelder-Mead types */

typedef enum nm_state {
    initial, reflect, expand, contract_out, contract_in, shrink, shrink2
} nm_state_t;


/* see: http://www.scholarpedia.org/article/Nelder-Mead_algorithm */
typedef struct nelder_mead {
    /* maintain search state */
    int dimensions;
    int iterations;
    nm_simplex_t simplex;
    fnt_vect_t seed;
    nm_state_t state;

    /* sample points */
    nm_sample_t x_r;
    nm_sample_t x_e;
    nm_sample_t x_c;
    fnt_vect_t s_shrink;

    /* hyper-parameters */
    double alpha;
    double beta;
    double gamma;
    double delta;

    /* termination criteria */
    double dist_threshold;
    int max_iterations;

} nelder_mead_t;


/* MARK: Nelder-Mead function */

/* \brief Provides the name of the method.
 * \param name Allocated buffer to hold the name.
 * \param size Size of the name buffer in bytes.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_name(char *name, int size) {
    if( strlcpy(name,"nelder-mead",size) >= size ) {
        return FNT_FAILURE;
    }

    return FNT_SUCCESS;
}


int method_init(void **nm_ptr, int dimensions) {
    nelder_mead_t *nm = calloc(1, sizeof(nelder_mead_t));
    *nm_ptr = nm;
    memset(nm, '\0', sizeof(*nm));

    nm->dimensions = dimensions;
    nm->iterations = 0;
    nm->state = initial;

    /* termination critria */
    nm->dist_threshold = 1e-5;
    nm->max_iterations = 30;

    nm->alpha = 1;      /* \alpha > 0 */
    nm->beta = 0.5;     /* 0 < \beta < 1 */
    nm->gamma = 2;      /* \gamma > 1 */
    nm->delta = 0.5;    /* 0 < \delta < 1 */

    fnt_vect_calloc(&nm->seed, dimensions);
    fnt_vect_calloc(&nm->x_r.parameters, dimensions);
    fnt_vect_calloc(&nm->x_e.parameters, dimensions);
    fnt_vect_calloc(&nm->x_c.parameters, dimensions);
    fnt_vect_calloc(&nm->s_shrink, dimensions);

    nm_simplex_init(&nm->simplex, dimensions);

    return FNT_SUCCESS;
}


int method_free(void **nm_ptr) {
    nelder_mead_t *nm = *nm_ptr;

    fnt_vect_free(&nm->seed);
    fnt_vect_free(&nm->x_r.parameters);
    fnt_vect_free(&nm->x_e.parameters);
    fnt_vect_free(&nm->x_c.parameters);
    fnt_vect_free(&nm->s_shrink);

    nm_simplex_free(&nm->simplex);

    free(nm); nm=NULL;

    return FNT_SUCCESS;
}


 
int method_info(void *nm_ptr) {
    /* argument isn't used, but is there for consistency. */

    /* TODO: Write info section for Nelder-Mead. */
    /* don't forget to describe hyper-parameters and their ranges */
    printf("Nelder-Mead info still needs to be filled in.\n");

    return FNT_FAILURE;
}


int method_hparam_set(void *nm_ptr, char *id, void *value_ptr) {
    if( nm_ptr == NULL )    { return FNT_FAILURE; }
    if( id == NULL )        { return FNT_FAILURE; }
    if( value_ptr == NULL ) { return FNT_FAILURE; }
    nelder_mead_t *nm = (nelder_mead_t*)nm_ptr;

    /* TODO: Add bounds checks for hparams */
    if( strncmp("alpha", id, 5) == 0 ) {
        nm->alpha = (double)(*((double*)value_ptr));
        return FNT_SUCCESS;
    }
    if( strncmp("beta", id, 4) == 0 ) {
        nm->beta = (double)(*((double*)value_ptr));
        return FNT_SUCCESS;
    }
    if( strncmp("gamma", id, 5) == 0 ) {
        nm->gamma = (double)(*((double*)value_ptr));
        return FNT_SUCCESS;
    }
    if( strncmp("delta", id, 5) == 0 ) {
        nm->delta = (double)(*((double*)value_ptr));
        return FNT_SUCCESS;
    }

    if( fnt_verbose_level >= 1 ) {
        fprintf(stderr, "No hyper-parameter '%s'.\n", id);
    }

    return FNT_FAILURE;
}


int method_hparam_get(void *nm_ptr, char *id, void *value_ptr) {
    if( nm_ptr == NULL )    { return FNT_FAILURE; }
    if( id == NULL )        { return FNT_FAILURE; }
    if( value_ptr == NULL ) { return FNT_FAILURE; }
    nelder_mead_t *nm = (nelder_mead_t*)nm_ptr;

    if( strncmp("alpha", id, 5) == 0 ) {
        *((double*)value_ptr) = nm->alpha;
        return FNT_SUCCESS;
    }
    if( strncmp("beta", id, 4) == 0 ) {
        *((double*)value_ptr) = nm->beta;
        return FNT_SUCCESS;
    }
    if( strncmp("gamma", id, 5) == 0 ) {
        *((double*)value_ptr) = nm->gamma;
        return FNT_SUCCESS;
    }
    if( strncmp("delta", id, 5) == 0 ) {
        *((double*)value_ptr) = nm->delta;
        return FNT_SUCCESS;
    }

    if( fnt_verbose_level >= 1 ) {
        fprintf(stderr, "No hyper-parameter '%s'.\n", id);
    }

    return FNT_FAILURE;
}


int method_seed(void *nm_ptr, fnt_vect_t *seed) {
    nelder_mead_t *nm = nm_ptr;
    if( nm->state != initial ) { return FNT_FAILURE; }
    fnt_vect_copy(&nm->seed, seed);

    return FNT_SUCCESS;
}


void nm_best_point(void *nm_ptr, fnt_vect_t *result) {
    nelder_mead_t *nm = nm_ptr;

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


int method_value(void *nm_ptr, fnt_vect_t *parameters, double value) {
    if( nm_ptr == NULL )        { return FNT_FAILURE; }
    if( parameters == NULL )    { return FNT_FAILURE; }
    if( parameters->v == NULL ) { return FNT_FAILURE; }

    nelder_mead_t *nm = nm_ptr;

    nm->iterations += 1;

    nm_sample_t new_sample;
    fnt_vect_calloc(&new_sample.parameters, parameters->n);
    fnt_vect_copy(&new_sample.parameters, parameters);
    new_sample.value = value;

    /* shrink just replaces points, but needs the associated values */
    if( nm->state == shrink2 ) {
        nm_sample_copy(&nm->simplex.points[nm->simplex.count-2], &new_sample);
        nm->state = reflect;
        fnt_vect_free(&new_sample.parameters);
        return FNT_SUCCESS;
    } else if( nm->state == shrink ) {
        nm_sample_copy(&nm->simplex.points[nm->simplex.count-1], &new_sample);
        nm->state = shrink2;
        fnt_vect_free(&new_sample.parameters);
        return FNT_SUCCESS;
    }

    /* check for initialization state */
    if( nm->simplex.count <= nm->dimensions ) {
        nm_simplex_add(&nm->simplex, &new_sample);
        if( nm->simplex.count >= nm->dimensions+1 )
            nm->state = reflect;
        fnt_vect_free(&new_sample.parameters);
        return FNT_SUCCESS;
    }

    /* sort simplex */
    if( nm->state != shrink && nm->state != shrink2 )
        nm_simplex_sort(&nm->simplex);

    /* get h, s, l */
    nm_sample_t h; /* simplex[simplex.count-1] */
    nm_sample_t s; /* simplex[simplex.count-2] */
    nm_sample_t l; /* simplex[0] */
    nm_sample_t r; /* new_sample */
    fnt_vect_calloc(&h.parameters, nm->dimensions);
    fnt_vect_calloc(&s.parameters, nm->dimensions);
    fnt_vect_calloc(&l.parameters, nm->dimensions);
    fnt_vect_calloc(&r.parameters, nm->dimensions);
    nm_sample_copy(&h, &nm->simplex.points[nm->simplex.count-1]);
    nm_sample_copy(&s, &nm->simplex.points[nm->simplex.count-2]);
    nm_sample_copy(&l, &nm->simplex.points[0]);
    nm_sample_copy(&r, &new_sample);
    fnt_vect_free(&new_sample.parameters);

    /* actual parameters can be discarded,
     * as only the output values are needed here */
    fnt_vect_free(&h.parameters);
    fnt_vect_free(&s.parameters);
    fnt_vect_free(&l.parameters);

    if( fnt_verbose_level >= 3 ) {
        printf("f(h) = %g\n", h.value);
        printf("f(s) = %g\n", s.value);
        printf("f(l) = %g\n", l.value);
        printf("f(r) = %g\n", r.value);
    }

    /* deal with recently computed point based on state */
    if( nm->state == reflect ) {
        nm_sample_copy(&nm->x_r, &r);

        if( l.value <= nm->x_r.value && nm->x_r.value < s.value ) {
            /* accept x_r and terminate iteration */
            nm_sample_copy(&nm->simplex.points[nm->simplex.count-1], &r);
            fnt_vect_free(&r.parameters);
            return FNT_SUCCESS;
        }
    }

    if( nm->state == expand ) {
        nm_sample_copy(&nm->x_e, &r);

        if( nm->x_e.value < nm->x_r.value ) {
            /* accept x_e and terminate iteration */
            nm_sample_copy(&nm->simplex.points[nm->simplex.count-1], &nm->x_e);
        } else {
            /* accept x_r and terminate iteration */
            nm_sample_copy(&nm->simplex.points[nm->simplex.count-1], &nm->x_r);
        }

        nm->state = reflect;
        fnt_vect_free(&r.parameters);
        return FNT_SUCCESS;
    }

    if( nm->state == contract_out ) {
        nm_sample_copy(&nm->x_c, &r);

        if( nm->x_c.value < nm->x_r.value ) {
            /* accept x_c and terminate iteration */
            nm_sample_copy(&nm->simplex.points[nm->simplex.count-1], &nm->x_c);
            nm->state = reflect;
            fnt_vect_free(&r.parameters);
            return FNT_SUCCESS;
        }
    }

    if( nm->state == contract_in ) {
        nm_sample_copy(&nm->x_c, &r);

        if( nm->x_c.value < h.value ) {
            /* accept x_c and terminate iteration */
            nm_sample_copy(&nm->simplex.points[nm->simplex.count-1], &nm->x_c);
            nm->state = reflect;
            fnt_vect_free(&r.parameters);
            return FNT_SUCCESS;
        }
    }
    fnt_vect_free(&r.parameters);

    /* determine next state if new point not accepted */
    if( r.value < l.value ) {
        /* cause x_e to be computed next */
        nm->state = expand;
        return FNT_SUCCESS;
    } else if( r.value >= s.value ) {
        /* cause x_c to be computed next */
        if( s.value <= r.value && r.value < h.value )
            nm->state = contract_out;
        else
            nm->state = contract_in;
        return FNT_SUCCESS;
    }

    nm->state = shrink;
    return FNT_SUCCESS;
}


int method_next(void *nm_ptr, fnt_vect_t *vector) {
    nelder_mead_t *nm = nm_ptr;

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

        return FNT_SUCCESS;
    }

    /* check to see that we have correct number of points */
    if( nm->simplex.count != nm->dimensions+1 ) {
        fnt_vect_copy(vector, &nm->seed);
        return FNT_SUCCESS;
    }

    /* sort simplex */
    if( nm->state != shrink && nm->state != shrink2 )
        nm_simplex_sort(&nm->simplex);

    /* get h, s, l */
    nm_sample_t h; /* simplex[simplex.count-1] */
    nm_sample_t s; /* simplex[simplex.count-2] */
    nm_sample_t l; /* simplex[0] */
    fnt_vect_calloc(&h.parameters, nm->dimensions);
    fnt_vect_calloc(&s.parameters, nm->dimensions);
    fnt_vect_calloc(&l.parameters, nm->dimensions);
    nm_sample_copy(&h, &nm->simplex.points[nm->simplex.count-1]);
    nm_sample_copy(&s, &nm->simplex.points[nm->simplex.count-2]);
    nm_sample_copy(&l, &nm->simplex.points[0]);

    /* compute centroid */
    fnt_vect_t centroid, sum_vect;
    fnt_vect_calloc(&centroid, nm->dimensions);
    fnt_vect_calloc(&sum_vect, nm->dimensions);
    for(int i=0; i<nm->simplex.count-1; ++i) {
        fnt_vect_add(&sum_vect, &nm->simplex.points[i].parameters, &sum_vect);
    }
    fnt_vect_scale(&sum_vect, 1.0/(nm->simplex.count-1), &centroid);
    fnt_vect_free(&sum_vect);

    /* add appropriate new point(s) */
    fnt_vect_t scaled, tmp;
    fnt_vect_calloc(&scaled, nm->dimensions);
    fnt_vect_calloc(&tmp, nm->dimensions);
    switch( nm->state ) {
        case initial:
            /* initial is handled above */
            if( fnt_verbose_level >= 1 ) {
                fprintf(stderr, "In initial state after bootstapping phase.\n");
                fprintf(stderr, "This should never happen!\n");
            }
            break;
        case reflect:
            fnt_vect_sub(&centroid, &h.parameters, &tmp);
            fnt_vect_scale(&tmp, nm->alpha, &scaled);
            fnt_vect_add(&centroid, &scaled, vector);    /* x_r */
            break;
        case expand:
            fnt_vect_sub(&nm->x_r.parameters, &centroid, &tmp);
            fnt_vect_scale(&tmp, nm->gamma, &scaled);
            fnt_vect_add(&centroid, &scaled, vector);    /* x_e */
            break;
        case contract_out:
            fnt_vect_sub(&nm->x_r.parameters, &centroid, &tmp);
            fnt_vect_scale(&tmp, nm->beta, &scaled);
            fnt_vect_add(&centroid, &scaled, vector);    /* x_c */
            break;
        case contract_in:
            fnt_vect_sub(&h.parameters, &centroid, &tmp);
            fnt_vect_scale(&tmp, nm->beta, &scaled);
            fnt_vect_add(&centroid, &scaled, vector);    /* x_c */
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
    fnt_vect_free(&centroid);
    fnt_vect_free(&h.parameters);
    fnt_vect_free(&s.parameters);
    fnt_vect_free(&l.parameters);

    if( fnt_verbose_level >= 3 ) {
        fnt_vect_println(vector, "next x ", "%.3f");
    }

    return FNT_SUCCESS;
}


int nm_simplex_point(void *nm_ptr, int which, fnt_vect_t *point, double *value) {
    nelder_mead_t *nm = nm_ptr;
    if( which >= nm->simplex.count )
        return 0;

    if( point )
        fnt_vect_copy(point, &nm->simplex.points[which].parameters);
    if( value )
        *value = nm->simplex.points[which].value;

    return 1;
}


int method_done(void *nm_ptr) {
    nelder_mead_t *nm = nm_ptr;
    if( nm->state == initial ) {
        return FNT_CONTINUE;
    }

    /* check maximum iterations */
    if( nm->iterations > nm->max_iterations ) {
        if( fnt_verbose_level >= 3 ) {
            printf("Iteration count (%i) exceeded limit (%i).\n", nm->iterations, nm->max_iterations); 
        }
        return FNT_DONE;
    }

    /* check distance between best and worst parameters */
    double dist;
    if( nm->state != shrink && nm->state != shrink2 )
        nm_simplex_sort(&nm->simplex);
    fnt_vect_dist(&nm->simplex.points[0].parameters,
                &nm->simplex.points[nm->simplex.count-1].parameters, &dist);
    #if 0
    printf("iteration: %i; dist: %g\n", nm->iterations, dist);
    #endif /* 0 */
    if( dist < nm->dist_threshold ) {
        if( fnt_verbose_level >= 3 ) {
            printf("Simplex size limit (%g) reached (%g).\n", nm->dist_threshold, dist); 
        }
        return FNT_DONE;
    }

    return FNT_CONTINUE;
}
