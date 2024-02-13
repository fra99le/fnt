/*
 * bisection.c
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../fnt.h"


/* MARK: Method type definitions */

typedef enum bisection_state {
    initial, initial2, running, done
} bisection_state_t;

typedef struct bisection {

    /* method state */
    bisection_state_t state;

    /* hyper-parameters */
    double upper_bound;
    double lower_bound;

    /* termination thresholds */
    double x_tol;
    double f_tol;

    /* current bounds */
    double a;
    double b;

    /* values at each bounds */
    double f_a;
    double f_b;

} bisection_t;


/* MARK: Functions called by FNT */

/* \brief Provides the name of the method.
 * \param name Allocated buffer to hold the name.
 * \param size Size of the name buffer in bytes.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_name(char *name, int size) {
    if( strlcpy(name,"bisection",size) >= size ) {
        return FNT_FAILURE;
    }

    return FNT_SUCCESS;
}


/* \brief Initialize intenal state for method.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_init(void **handle_ptr, int dimensions) {
    if( handle_ptr == NULL )    { return FNT_FAILURE; }
    bisection_t *ptr = calloc(1, sizeof(bisection_t));
    if( ptr == NULL )           { return FNT_FAILURE; }
    *handle_ptr = (void*)ptr;

    /* initialize method here */
    ptr->state = initial;
    ptr->x_tol = 1e-6;
    ptr->f_tol = 1e-6;
    ptr->lower_bound = -1e6;
    ptr->upper_bound = 1e6;

    return FNT_SUCCESS;
}


int method_free(void **handle_ptr) {
    if( handle_ptr == NULL )    { return FNT_FAILURE; }
    if( *handle_ptr == NULL )   { return FNT_FAILURE; }
    bisection_t *ptr = (bisection_t*)*handle_ptr;

    /* free any memory allocated by method */

    free(ptr);  *handle_ptr = ptr = NULL;

    return FNT_SUCCESS;
}


int method_info() {
    printf("\n"
           "This should give usefiul information about the method.\n"
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
    if( handle == NULL )    { return FNT_FAILURE; }
    bisection_t *ptr = (bisection_t*)handle;

    FNT_HPARAM_SET("f_tol", id, double, value_ptr, ptr->f_tol);
    FNT_HPARAM_SET("x_tol", id, double, value_ptr, ptr->x_tol);
    FNT_HPARAM_SET("lower", id, double, value_ptr, ptr->lower_bound);
    FNT_HPARAM_SET("upper", id, double, value_ptr, ptr->upper_bound);

    ptr->a = ptr->lower_bound;
    ptr->b = ptr->upper_bound;

    return FNT_SUCCESS;
}


int method_hparam_get(void *handle, char *id, void *value_ptr) {
    if( id == NULL )        { return FNT_FAILURE; }
    if( value_ptr == NULL ) { return FNT_FAILURE; }
    bisection_t *ptr = (bisection_t*)handle;

    FNT_HPARAM_GET("f_tol", id, double, ptr->f_tol, value_ptr);
    FNT_HPARAM_GET("x_tol", id, double, ptr->x_tol, value_ptr);
    FNT_HPARAM_GET("lower", id, double, ptr->lower_bound, value_ptr);
    FNT_HPARAM_GET("upper", id, double, ptr->upper_bound, value_ptr);

    return FNT_FAILURE;
}


int method_next(void *handle, fnt_vect_t *vec) {
    if( handle == NULL )    { return FNT_FAILURE; }
    if( vec == NULL )       { return FNT_FAILURE; }
    bisection_t *ptr = (bisection_t*)handle;

    if( ptr->state == initial )     { vec->v[0] = ptr->a; return FNT_SUCCESS; }
    if( ptr->state == initial2 )    { vec->v[0] = ptr->b; return FNT_SUCCESS; }

    /* fill vector pointed to by vec with next input to try */
    vec->v[0] = 0.5*ptr->a + 0.5*ptr->b;

    return FNT_SUCCESS;
}


int method_value(void *handle, fnt_vect_t *vec, double value) {
    if( handle == NULL )    { return FNT_FAILURE; }
    bisection_t *ptr = (bisection_t*)handle;

    /* update method using value */
    if( ptr->state == initial2 ) {
        ptr->f_b = value;

        /* ensure that f(a) < f(b) */
        if( ptr->f_b < ptr->f_a ) {
            double tmp;
            tmp = ptr->a;       ptr->b = ptr->a;        ptr->a = tmp;

            tmp = ptr->f_a;     ptr->f_b = ptr->f_a;    ptr->f_a = tmp;
        }

        /* check that endpoints meet bisection precondition, f(a)*f(b) < 0 */
        if( ptr->f_a > 0.0 ) {
            if( fnt_verbose_level >= FNT_ERROR ) {
                fprintf(stderr, "Lower bound is not less than zero (f(%g)=%g)\n", ptr->a, ptr->f_a); 
            }
            return FNT_FAILURE;
        }
        if( ptr->f_b < 0.0 ) {
            if( fnt_verbose_level >= FNT_ERROR ) {
                fprintf(stderr, "Upper bound is not greater than zero (f(%g)=%g)\n", ptr->b, ptr->f_b); 
            }
            return FNT_FAILURE;
        }

        ptr->state = running;
        return FNT_SUCCESS;
    } else if( ptr->state == initial ) {
        ptr->f_a = value;
        ptr->state = initial2;
        return FNT_SUCCESS;
    } 

    if( ptr->state != running ) { return FNT_FAILURE; }

    if( value < 0.0 ) {
        ptr->a = vec->v[0];
        ptr->f_a = value;
    } else if( value > 0.0 ) {
        ptr->b = vec->v[0];
        ptr->f_b = value;
    } else if( value == 0.0 ) {
        ptr->a = vec->v[0];
        ptr->a = vec->v[0];
        ptr->f_a = 0.0;
        ptr->f_b = 0.0;
        ptr->state = done;
    } else {
        if( fnt_verbose_level >= FNT_ERROR ) {
            fprintf(stderr, "Value (%g) is not comparable to zero.\n", value);
        }
        return FNT_FAILURE;
    }

    return FNT_SUCCESS;
}


int method_value_gradient(void *handle, fnt_vect_t *vec, double value, fnt_vect_t gradient) {
    if( handle == NULL )    { return FNT_FAILURE; }

    /* update method using value and derivative/gradient */
    return method_value(handle, vec, value);
}


int method_done(void *handle) {
    if( handle == NULL )    { return FNT_FAILURE; }
    bisection_t *ptr = (bisection_t*)handle;

    /* test for completion
     *  Return FNT_DONE when comlete, or FNT_CONTINUE when not done.
     */
    if( ptr->state == done ) {
        return FNT_DONE;
    }

    if( ptr->b - ptr->a  < ptr->x_tol ) {
        if( fnt_verbose_level >= FNT_INFO ) {
            printf("Upper and lower bound within termination threshold.\n");
        }
        ptr->state = done;
        return FNT_DONE;
    }
    if( ptr->f_b - ptr->f_a  < ptr->f_tol ) {
        if( fnt_verbose_level >= FNT_INFO ) {
            printf("Difference in function's value at upper and lower bound within termination threshold.\n");
        }
        ptr->state = done;
        return FNT_DONE;
    }

    return FNT_CONTINUE;
}


int method_result(void *handle, void *extra) {
    if( handle == NULL )    { return FNT_FAILURE; }
    bisection_t *ptr = (bisection_t*)handle;

    /* Optional method to report any additional results if the method
     * produces such results.
     */

    return FNT_FAILURE;
}
