/*
 * de.c
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../fnt.h"

/* MARK: Method type definitions */

typedef enum de_state {
    de_initial, de_running, de_done
} de_state_t;

typedef struct de {

    int dim;    /* number of dimensions in parameter vectors */
    de_state_t state;
    int allocated_NP;

    /* hyper parameters */
    int iterations;
    int NP;
    double F;
    double lambda;
    fnt_vect_t start_point;
    fnt_vect_t lower_bounds;
    fnt_vect_t upper_bounds;
    int has_start_point;
    int has_lower_bounds;
    int has_upper_bounds;

    /* current generation */
    fnt_vect_t *x;
    fnt_vect_t *x_prev;
    double *fx;
    double *fx_prev;
    int best;

    /* trial vector */
    fnt_vect_t v;
    int current;    /* index of vector that v might replace */

    /* results */
    double min_fx;
    fnt_vect_t min_x;
} de_t;


/* MARK: Internal functions */

static int de_print_generation(de_t *ptr) {

    printf("Current generation:\n");
    for(int i=0; i<ptr->NP; ++i) {
        printf("%5d: ", i);
        fnt_vect_print(&ptr->x[i], "f(", "%.4f");
        printf(") -> %g\n", ptr->fx[i]);
    }

    printf("Previous generation:\n");
    for(int i=0; i<ptr->NP; ++i) {
        printf("%5d: ", i);
        fnt_vect_print(&ptr->x_prev[i], "f(", "%.4f");
        printf(") -> %g\n", ptr->fx_prev[i]);
    }

    return FNT_SUCCESS;
}


static int de_allocate_generations(de_t *ptr) {

    /* TODO: Check for calloc failures */
    ptr->x = calloc(ptr->NP, sizeof(fnt_vect_t));
    ptr->x_prev = calloc(ptr->NP, sizeof(fnt_vect_t));
    for(int i=0; i<ptr->NP; ++i) {
        fnt_vect_calloc(&ptr->x[i], ptr->dim);
        fnt_vect_calloc(&ptr->x_prev[i], ptr->dim);
    }
    ptr->fx = calloc(ptr->NP, sizeof(double));
    ptr->fx_prev = calloc(ptr->NP, sizeof(double));

    if( fnt_verbose_level >= FNT_DEBUG ) {
        de_print_generation(ptr);
    }

    return FNT_SUCCESS;
}


static int de_free_generations(de_t *ptr) {

    for(int i=0; i<ptr->NP; ++i) {
        fnt_vect_free(&ptr->x[i]);
        fnt_vect_free(&ptr->x_prev[i]);
    }
    free(ptr->x); ptr->x=NULL;
    free(ptr->x_prev); ptr->x_prev=NULL;
    free(ptr->fx); ptr->fx=NULL;
    free(ptr->fx_prev); ptr->fx_prev=NULL;

    return FNT_SUCCESS;
}


static int de_fill_first_gen(de_t *ptr) {

    int curr = ptr->current;

    /* if an initial value is give, fill based on it. */
    if( ptr->has_start_point ) {
        if( fnt_verbose_level >= FNT_DEBUG ) {
            DEBUG("Filling initial generation using ");
            fnt_vect_print(&ptr->start_point, "start point: ", NULL);
            DEBUG(".\n");
        }

        /* pick a point normally distributed around start point */
        for(int j=0; j<ptr->dim; ++j) {
            /* TODO: rnd should be normally distributed with a hyper-parameter
             * specifying the std. dev. */
            double rnd = FNT_RAND() / (double)FNT_RAND_MAX - 0.5;
            FNT_VECT_ELEM(ptr->v, j) = FNT_VECT_ELEM(ptr->start_point, j)  + rnd;

            /* apply bounds, as available */
            if( ptr->has_lower_bounds
                && FNT_VECT_ELEM(ptr->v, j) < FNT_VECT_ELEM(ptr->lower_bounds, j) ) {
                FNT_VECT_ELEM(ptr->v, j) = FNT_VECT_ELEM(ptr->lower_bounds, j);
            }
            if( ptr->has_upper_bounds
                && FNT_VECT_ELEM(ptr->v, j) > FNT_VECT_ELEM(ptr->upper_bounds, j) ) {
                FNT_VECT_ELEM(ptr->v, j) = FNT_VECT_ELEM(ptr->upper_bounds, j);
            }
        }
    } else {
        DEBUG("Filling initial generation uniformly randomly (curr=%d).\n", curr);

        /* pick uniformly random point, applying bounds as supplied. */
        for(int j=0; j<ptr->dim; ++j) {
            double rnd = FNT_RAND() / (double)FNT_RAND_MAX;
            double lower = -1.0;
            double upper = 1.0;
            if( ptr->has_lower_bounds) {
                lower = FNT_VECT_ELEM(ptr->lower_bounds, j);
                if( !ptr->has_upper_bounds ) {
                    upper = lower + 1.0;
                }
            }
            if( ptr->has_upper_bounds) {
                upper = FNT_VECT_ELEM(ptr->upper_bounds, j);
                if( !ptr->has_lower_bounds ) {
                    lower = lower - 1.0;
                }
            }
            FNT_VECT_ELEM(ptr->v, j) = lower + rnd * (upper-lower);
        }
    }

    return FNT_SUCCESS;
}



/* MARK: Functions called by FNT */

/* \brief Provides the name of the method.
 * \param name Allocated buffer to hold the name.
 * \param size Size of the name buffer in bytes.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_name(char *name, int size) {
    if( snprintf(name, size, "differential evolution") >= size ) {
        return FNT_FAILURE;
    }

    return FNT_SUCCESS;
}


/* \brief Initialize intenal state for method.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_init(void **handle_ptr, int dimensions) {
    de_t *ptr = calloc(1, sizeof(de_t));

    /* record dimensionality */
    ptr->dim = dimensions;
    ptr->state = de_initial;

    /* set up method */
    ptr->iterations = 1000;
    ptr->NP = dimensions * 10;
    ptr->allocated_NP = ptr->NP;
    ptr->F = 0.5;
    ptr->lambda = 0.1;

    /* allocate generations */
    de_allocate_generations(ptr);
    fnt_vect_calloc(&ptr->v, dimensions);
    ptr->current = 0;
    ptr->best = 0;

    /* allocate/initialize results */
    fnt_vect_calloc(&ptr->min_x, dimensions);
    ptr->min_fx = 0.0;

    *handle_ptr = (void*)ptr;

    return FNT_SUCCESS;
}


int method_free(void **handle_ptr) {
    if( handle_ptr == NULL )    { return FNT_FAILURE; }
    if( *handle_ptr == NULL )   { return FNT_FAILURE; }
    de_t *ptr = (de_t*)*handle_ptr;

    /* free generation tracking */
    fnt_vect_free(&ptr->v);
    de_free_generations(ptr);

    /* free vectors, if allocated */
    if( ptr->has_start_point )  { fnt_vect_free(&ptr->start_point);  }
    if( ptr->has_lower_bounds ) { fnt_vect_free(&ptr->lower_bounds); }
    if( ptr->has_upper_bounds ) { fnt_vect_free(&ptr->upper_bounds); }

    /* free results */
    fnt_vect_free(&ptr->min_x);

    free(ptr);  *handle_ptr = ptr = NULL;

    return FNT_SUCCESS;
}


/* \brief Display information about the method to the console.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_info() {
    printf(
"Differential evolution is a minimization method that uses a population of\n"
"randomized guesses that are systematically updated with better guesses until\n"
"a minimum value it found.\n"
"\n"
"Note: crossover is not currently implemented.\n"
"\n"
"Hyper-parameters:\n"
"name\trequired\ttype\t\tDefault\tDescription\n"
"lower\toptional\tfnt_vect_t\tnone\tLower bounds on search region.\n"
"upper\toptional\tfnt_vect_t\tnone\tUpper bounds on search region.\n"
"start\toptional\tfnt_vect_t\tnone\tCenter of initial search region.\n"
"NP\tREQUIRED\tint\t\t10*dims\tNumber of random points.\n"
"F\toptional\tint\t\t0\tScaling factor applied to difference of vectors.\n"
"lambda\toptional\tint\t\t0\tScaling factor applied to best vector difference.\n"
"iterations\toptional\tint\t\t1000\tNumber of iterations to run.\n"
"\n"
"References:\n"
"Storn, R., Price, K. Differential Evolution – A Simple and Efficient\n"
"\tHeuristic for global Optimization over Continuous Spaces.\n"
"\tJournal of Global Optimization 11, 341–359 (1997).\n"
"\thttps://doi.org/10.1023/A:1008202821328\n"
);
    return FNT_SUCCESS;
}

static int validate_hparams(de_t *ptr) {

    if( (ptr->has_lower_bounds && ptr->has_upper_bounds) ) {
        for(int j=0; j<ptr->lower_bounds.n; ++j) {
            double lower = FNT_VECT_ELEM(ptr->lower_bounds, j);
            double upper = FNT_VECT_ELEM(ptr->upper_bounds, j);
            if( upper < lower ) {
                WARN("WARNING: Upper and lower bounds for dimension %i are out of order (lower=%g, upper=%g), swapping them.\n", j, lower, upper);
                FNT_VECT_ELEM(ptr->lower_bounds, j) = upper;
                FNT_VECT_ELEM(ptr->upper_bounds, j) = lower;
            }
        }
    }

    if( ptr->NP < 3 ) {
        ERROR("ERROR: NP must be at least 3, NP was %d, changing it to 3.\n", ptr->NP);
        ptr->NP = 3;
    }

    /* resize generation, if NP changed */
    if( ptr->NP != ptr->allocated_NP ) {
        /* free old generation tracking */
        de_free_generations(ptr);

        /* allocate new generation tracking */
        de_allocate_generations(ptr);
    }

    return FNT_SUCCESS;
}


/* \brief Set any hyper-parameters needed for the method.
 * \param id The name of the hyper-parameter.
 * \param value_ptr A pointer to the value being set.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_hparam_set(void *handle, char *id, void *value_ptr) {
    de_t *ptr = (de_t*)handle;
    if( ptr == NULL )       { return FNT_FAILURE; }
    if( id == NULL )        { return FNT_FAILURE; }
    if( value_ptr == NULL ) { return FNT_FAILURE; }

    FNT_HPARAM_SET("iterations", id, int, value_ptr, ptr->iterations);
    FNT_HPARAM_SET("F", id, double, value_ptr, ptr->F);
    FNT_HPARAM_SET("lambda", id, double, value_ptr, ptr->lambda);
    FNT_HPARAM_SET("NP", id, int, value_ptr, ptr->NP);

    if( strncmp("start", id, 5) == 0 ) {
        if( !ptr->has_start_point ) {
            fnt_vect_calloc(&ptr->start_point, ptr->dim);
        }
        fnt_vect_copy(&ptr->start_point, value_ptr);
        ptr->has_start_point = 1;

        return FNT_SUCCESS;
    }

    if( strncmp("lower", id, 5) == 0 ) {
        if( !ptr->has_lower_bounds ) {
            fnt_vect_calloc(&ptr->lower_bounds, ptr->dim);
        }
        fnt_vect_copy(&ptr->lower_bounds, value_ptr);
        ptr->has_lower_bounds = 1;

        return FNT_SUCCESS;
    }

    if( strncmp("upper", id, 5) == 0 ) {
        if( !ptr->has_upper_bounds ) {
            fnt_vect_calloc(&ptr->upper_bounds, ptr->dim);
        }
        fnt_vect_copy(&ptr->upper_bounds, value_ptr);
        ptr->has_upper_bounds = 1;

        return FNT_SUCCESS;
    }

    ERROR("No hyper-parameter named '%s'.\n", id);

    return FNT_FAILURE;
}


int method_hparam_get(void *handle, char *id, void *value_ptr) {
    de_t *ptr = (de_t*)handle;
    if( ptr == NULL )       { return FNT_FAILURE; }
    if( id == NULL )        { return FNT_FAILURE; }
    if( value_ptr == NULL ) { return FNT_FAILURE; }

    FNT_HPARAM_GET("F", id, double, ptr->F, value_ptr);
    FNT_HPARAM_GET("lambda", id, double, ptr->lambda, value_ptr);
    FNT_HPARAM_GET("NP", id, int, ptr->NP, value_ptr);

    if( strncmp("start", id, 5) == 0 ) {
        if( ptr->has_start_point ) {
            return fnt_vect_copy(value_ptr, &ptr->start_point);
        } else {
            ERROR("Start point requested, but not set.\n");
            return FNT_FAILURE;
        }
    }
    if( strncmp("lower", id, 5) == 0 ) {
        if( ptr->has_lower_bounds ) {
            return fnt_vect_copy(value_ptr, &ptr->lower_bounds);
        } else {
            ERROR("Lower bound requested, but not set.\n");
            return FNT_FAILURE;
        }
    }
    if( strncmp("upper", id, 5) == 0 ) {
        if( ptr->has_upper_bounds ) {
            return fnt_vect_copy(value_ptr, &ptr->upper_bounds);
        } else {
            ERROR("Upper bound requested, but not set.\n");
            return FNT_FAILURE;
        }
    }

    ERROR("No hyper-parameter named '%s'.\n", id);

    return FNT_FAILURE;
}


int method_next(void *handle, fnt_vect_t *vec) {
    de_t *ptr = (de_t*)handle;
    if( ptr == NULL )   { return FNT_FAILURE; }

    int curr = ptr->current;

    /* fill initial generation during initialization phase */
    if( ptr->state == de_initial ) {
        validate_hparams(ptr);

        de_fill_first_gen(ptr);
        return fnt_vect_copy(vec, &ptr->v);
    }

    if( ptr->state != de_running ) {
        ERROR("%s called while in the wrong state.\n", __FUNCTION__);

        return FNT_FAILURE;
    }

    /* pick unique r1, r2, r3 vectors */
    int r1 = FNT_RAND() % ptr->NP;
    int r2 = FNT_RAND() % ptr->NP;
    int r3 = FNT_RAND() % ptr->NP;
    while( r1 == r2 )               { r2 = FNT_RAND() % ptr->NP; }
    while( r1 == r3 || r2 == r3 )   { r3 = FNT_RAND() % ptr->NP; }
    DEBUG("DEBUG: r1, r2, r3 = %d, %d, %d\n", r1, r2, r3);

    /* compute trial vector v */
    fnt_vect_t diff;    /* Note: these could reside in de_t */
    fnt_vect_t scaled;
    fnt_vect_calloc(&diff, ptr->dim);
    fnt_vect_calloc(&scaled, ptr->dim);
    fnt_vect_t *x_prev = ptr->x_prev;
    if( ptr->lambda != 0.0 ) {
        /* scheme DE2 */
        fnt_vect_sub(&x_prev[ptr->best], &x_prev[curr], &diff);
        fnt_vect_scale(&diff, ptr->lambda, &scaled);
        fnt_vect_add(&x_prev[curr], &scaled, &ptr->v);

        fnt_vect_sub(&x_prev[r2], &x_prev[r3], &diff);
        fnt_vect_scale(&diff, ptr->F, &scaled);
        fnt_vect_add(&ptr->v, &scaled, &ptr->v);
    } else if( ptr->F != 0.0 ) {
        /* scheme DE1 */
        fnt_vect_sub(&x_prev[r2], &x_prev[r3], &diff);
        fnt_vect_scale(&diff, ptr->F, &scaled);
        fnt_vect_add(&x_prev[r1], &scaled, &ptr->v);
    }

    /* apply crossover */
    /* TODO: Add crossover */

    /* apply lower and upper bounds */
    if( ptr->has_lower_bounds ) {
        for(int j=0; j<ptr->v.n; ++j) {
            if( FNT_VECT_ELEM(ptr->v, j) < FNT_VECT_ELEM(ptr->lower_bounds, j) ) {
                FNT_VECT_ELEM(ptr->v, j) = FNT_VECT_ELEM(ptr->lower_bounds, j);
            }
        }
    }
    if( ptr->has_upper_bounds ) {
        for(int j=0; j<ptr->v.n; ++j) {
            if( FNT_VECT_ELEM(ptr->v, j) > FNT_VECT_ELEM(ptr->upper_bounds, j) ) {
                FNT_VECT_ELEM(ptr->v, j) = FNT_VECT_ELEM(ptr->upper_bounds, j);
            }
        }
    }

    return fnt_vect_copy(vec, &ptr->v);
}


int method_value(void *handle, fnt_vect_t *vec, double value) {
    de_t *ptr = (de_t*)handle;
    if( ptr == NULL )   { return FNT_FAILURE; }

    /* replace parameter vector with v, if warranted */
    int curr = ptr->current;
    if( value < ptr->fx_prev[curr] || ptr->state == de_initial ) {
        fnt_vect_copy(&ptr->x[curr], vec);
        ptr->fx[curr] = value;
        if( ptr->state == de_initial ) { ptr->state = de_running; }
    } else {
        fnt_vect_copy(&ptr->x[curr], &ptr->x_prev[curr]);
        ptr->fx[curr] = ptr->fx_prev[curr];
    }

    /* fx[curr] and x[curr] are now set correctly */

    /* compare against current best value */
    if( value < ptr->fx[ptr->best] ) {
        if( fnt_verbose_level >= FNT_INFO ) {
            INFO("New best value %g ", value);
            fnt_vect_print(vec, "for input ", NULL);
            INFO(" at position %d.\n", curr);
        }

        /* update best */
        ptr->best = curr;
    }

    /* move to next memeber of the current generation */
    ++ptr->current;

    /* update generation, as needed */
    if( ptr->current >= ptr->NP ) {

        /* leave initial state once first generation is complete */
        if( ptr->state == de_initial ) {
            DEBUG("Finished initial generation of size %d.\n", ptr->NP);
            ptr->state = de_running;
        }

        /* finished current generation, swap */
        void *tmp;
        DEBUG("DEBUG: Swapping generations.\n");
        tmp = ptr->x;   ptr->x = ptr->x_prev;       ptr->x_prev = tmp;
        tmp = ptr->fx;  ptr->fx = ptr->fx_prev;     ptr->fx_prev = tmp;

        ptr->current = 0;

        if( fnt_verbose_level >= FNT_DEBUG ) {
            DEBUG("After swap:\n");
            de_print_generation(ptr);
        }

        --ptr->iterations;
    }

    return FNT_SUCCESS;
}


int method_done(void *handle) {
    de_t *ptr = (de_t*)handle;
    if( ptr == NULL )   { return FNT_FAILURE; }

    if( ptr->state == de_initial ) {
        return FNT_CONTINUE;
    }

    if( ptr->iterations <= 0 ) {

        /* update result fields */
        ptr->min_fx = ptr->fx[ptr->best];
        fnt_vect_copy(&ptr->min_x, &ptr->x[ptr->best]);

        /* mark method as complete */
        ptr->state = de_done;

        return FNT_DONE;
    }

    return FNT_CONTINUE;
}


int method_result(void *handle, char *id, void *value_ptr) {
    de_t *ptr = (de_t*)handle;
    if( ptr == NULL )       { return FNT_FAILURE; }

    FNT_RESULT_GET_VECT("minimum x", id, ptr->min_x, value_ptr);
    FNT_RESULT_GET("minimum f", id, double, ptr->min_fx, value_ptr);

    ERROR("No result named '%s'.\n", id);

    return FNT_FAILURE;
}
