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

/* MARK: Method type defitinions */

typedef enum de_state {
    de_initial, de_running, de_done
} de_state_t;

typedef struct de {
    int dim;    /* number of dimensions in parameter vectors */
    de_state_t state;
    double f_tol;

    /* hyper parameters */
    double F;
    double lambda;

    /* current generation */
    int NP;
    fnt_vect_t *x;
    fnt_vect_t *x_prev;
    double *fx;
    double *fx_prev;
    int best;

    /* trial vector */
    fnt_vect_t v;
    int current;    /* index of vector that v might replace */
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

    /* if bounds are providesd, randomly fill based on those */

    /* otherwise just pick random values */
    if( fnt_verbose_level >= FNT_DEBUG ) {
        printf("Filling initial generation randomly (curr=%d).\n", curr);
    }
    for(int j=0; j<ptr->dim; ++j) {
        ptr->v.v[j] = rand() / (double)RAND_MAX - 0.5;
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
    if( strlcpy(name,"differential evolution",size) >= size ) {
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
    ptr->f_tol = 1e-6;
    ptr->NP = dimensions * 10;
    ptr->F = 0.5;
    ptr->lambda = 0.1;

    /* allocate generations */
    de_allocate_generations(ptr);
    fnt_vect_calloc(&ptr->v, dimensions);
    ptr->current = 0;
    ptr->best = 0;

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

    free(ptr);  *handle_ptr = ptr = NULL;

    return FNT_SUCCESS;
}


int method_info() {
    printf("\n"
           "No information for this method, yet.\n"
           "\n");
    return FNT_SUCCESS;
}


/* \brief Set any hyper-parameters needed for the method.
 * \param id The name of the hyper-parameter.
 * \param value_ptr A pointer to the value being set.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_hparam_set(void *handle, char *id, void *value_ptr) {
    if( id == NULL )        { return FNT_FAILURE; }
    if( value_ptr == NULL ) { return FNT_FAILURE; }
    de_t *ptr = (de_t*)handle;

    int old_NP = ptr->NP;

    FNT_HPARAM_SET("f_tol", id, double, value_ptr, ptr->f_tol);
    FNT_HPARAM_SET("F", id, double, value_ptr, ptr->F);
    FNT_HPARAM_SET("lambda", id, double, value_ptr, ptr->lambda);
    FNT_HPARAM_SET("NP", id, int, value_ptr, ptr->NP);

    if( ptr->NP < 3 ) {
        if( fnt_verbose_level >= FNT_ERROR ) {
            fprintf(stderr, "ERROR: NP must be at least 3.\n");
        }
        ptr->NP = 3;
    }

    /* resize generation, if NP changed */
    if( ptr->NP != old_NP ) {
        /* free old generation tracking */
        de_free_generations(ptr);

        /* allocate new generation tracking */
        de_allocate_generations(ptr);
    }
        
    return FNT_SUCCESS;
}


int method_hparam_get(void *handle, char *id, void *value_ptr) {
    if( id == NULL )        { return FNT_FAILURE; }
    if( value_ptr == NULL ) { return FNT_FAILURE; }
    de_t *ptr = (de_t*)handle;

    FNT_HPARAM_GET("F", id, double, ptr->F, value_ptr);
    FNT_HPARAM_GET("lambda", id, double, ptr->lambda, value_ptr);
    FNT_HPARAM_GET("NP", id, int, ptr->NP, value_ptr);

    return FNT_SUCCESS;
}


int method_next(void *handle, fnt_vect_t *vec) {
    de_t *ptr = (de_t*)handle;
    if( ptr == NULL )   { return FNT_FAILURE; }

    int curr = ptr->current;

    /* fill initial generation on first call */
    if( ptr->state == de_initial ) {
        de_fill_first_gen(ptr);
        return fnt_vect_copy(vec, &ptr->v);
    }

    if( ptr->state != de_running ) {
        if( fnt_verbose_level >= FNT_WARN ) {
            fprintf(stderr, "%s called while in the wrong state.\n", __FUNCTION__);
        }

        return FNT_FAILURE;
    }

    /* pick r1, r2, r3 unique */
    int r1, r2, r3;
    r1 = rand() % ptr->NP;
    r2 = rand() % ptr->NP;
    r3 = rand() % ptr->NP;
    while( r1 == r2 )               { r2 = rand() % ptr->NP; }
    while( r1 == r3 || r2 == r3 )   { r3 = rand() % ptr->NP; }
    if( fnt_verbose_level >= FNT_DEBUG ) {
        printf("r1, r2, r3 = %d, %d, %d\n", r1, r2, r3);
    }

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
    } else {
        fnt_vect_copy(&ptr->x[curr], &ptr->x_prev[curr]);
        ptr->fx[curr] = ptr->fx_prev[curr];
    }

    /* fx[curr] and x[curr] are now set correctly */

    /* compare against current best value */
    if( value < ptr->fx[ptr->best] ) {
        if( fnt_verbose_level >= FNT_INFO ) {
            printf("New best value %g ", value);
            fnt_vect_print(vec, "for input ", NULL);
            printf(" at position %d.\n", curr);
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
            if( fnt_verbose_level >= FNT_DEBUG ) {
                printf("Finished initial generation of size %d.\n", ptr->NP);
            }
            ptr->state = de_running;
        }

        /* finished current generation, swap */
        void *tmp;
        if( fnt_verbose_level >= FNT_DEBUG ) {
            printf("DEBUG: Swapping generations.\n");
        }
        tmp = ptr->x;   ptr->x = ptr->x_prev;       ptr->x_prev = tmp;
        tmp = ptr->fx;  ptr->fx = ptr->fx_prev;     ptr->fx_prev = tmp;

        ptr->current = 0;

        if( fnt_verbose_level >= FNT_DEBUG ) {
            printf("After swap:\n");
            de_print_generation(ptr);
        }
    }

    return FNT_SUCCESS;
}


int method_done(void *handle) {
    de_t *ptr = (de_t*)handle;
    if( ptr == NULL )   { return FNT_FAILURE; }

    if( ptr->state == de_initial ) {
        return FNT_CONTINUE;
    }

    if( ptr->fx[ptr->best] < ptr->f_tol
        && ptr->fx_prev[ptr->best] < ptr->f_tol ) {
        ptr->state = de_done;
        return FNT_DONE;
    }

    return FNT_CONTINUE;
}


int method_result(void *handle, void *extra) {
    return FNT_FAILURE;
}
