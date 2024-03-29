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

    /* result */
    double root_x;

} bisection_t;


/* MARK: Functions called by FNT */

/* \brief Provides the name of the method.
 * \param name Allocated buffer to hold the name.
 * \param size Size of the name buffer in bytes.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_name(char *name, int size) {
    if( snprintf(name, size, "bisection") >= size ) {
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


/* \brief Display information about the method to the console.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_info() {
    printf(
"The bisection method is a root finding technique that works by repeatedly "
"dividing a search region in half until it converges on the root."
"\n"
"Hyper-parameters:\n"
"name\trequired\ttype\tDefault\tDescription\n"
"lower\tREQUIRED\tdouble\t-1e6\tLower bound of the region.\n"
"upper\tREQUIRED\tdouble\t1e6\tUpper bound of the region.\n"
"f_tol\toptional\tdouble\t1e-6\tTerminates when |f(x)| < f_tol.\n"
"x_tol\toptional\tdouble\t1e-6\tTerminates when |a-b| < x_tol.\n"
"\n"
"References\n"
"https://en.wikipedia.org/wiki/Bisection_method"
);
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

    ERROR("No hyper-parameter named '%s'.\n", id);

    return FNT_FAILURE;
}


int method_hparam_get(void *handle, char *id, void *value_ptr) {
    if( id == NULL )        { return FNT_FAILURE; }
    if( value_ptr == NULL ) { return FNT_FAILURE; }
    bisection_t *ptr = (bisection_t*)handle;

    FNT_HPARAM_GET("f_tol", id, double, ptr->f_tol, value_ptr);
    FNT_HPARAM_GET("x_tol", id, double, ptr->x_tol, value_ptr);
    FNT_HPARAM_GET("lower", id, double, ptr->lower_bound, value_ptr);
    FNT_HPARAM_GET("upper", id, double, ptr->upper_bound, value_ptr);

    ERROR("No hyper-parameter named '%s'.\n", id);

    return FNT_FAILURE;
}


int method_next(void *handle, fnt_vect_t *vec) {
    if( handle == NULL )    { return FNT_FAILURE; }
    if( vec == NULL )       { return FNT_FAILURE; }
    bisection_t *ptr = (bisection_t*)handle;

    if( ptr->state == initial ) {
        ptr->a = ptr->lower_bound;
        ptr->b = ptr->upper_bound;

        FNT_VECT_ELEM(*vec, 0) = ptr->a;

        return FNT_SUCCESS;
    }
    if( ptr->state == initial2 ) {
        FNT_VECT_ELEM(*vec, 0) = ptr->b;

        return FNT_SUCCESS;
    }

    /* fill vector pointed to by vec with next input to try */
    FNT_VECT_ELEM(*vec, 0) = 0.5*ptr->a + 0.5*ptr->b;

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
            /* swap a & b if order of f(a) & f(b) are awkward. */
            double tmp;
            tmp = ptr->b;       ptr->b = ptr->a;        ptr->a = tmp;
            tmp = ptr->f_b;     ptr->f_b = ptr->f_a;    ptr->f_a = tmp;
        }

        /* check that endpoints meet bisection precondition, f(a)*f(b) < 0 */
        if( ptr->f_a > 0.0 ) {
            ERROR("Lower bound is not less than zero (f(%g)=%g)\n", ptr->a, ptr->f_a); 
            return FNT_FAILURE;
        }
        if( ptr->f_b < 0.0 ) {
            ERROR("Upper bound is not greater than zero (f(%g)=%g)\n", ptr->b, ptr->f_b); 
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
        ptr->a = FNT_VECT_ELEM(*vec, 0);
        ptr->f_a = value;
    } else if( value > 0.0 ) {
        ptr->b = FNT_VECT_ELEM(*vec, 0);
        ptr->f_b = value;
    } else if( value == 0.0 ) {
        ptr->a = FNT_VECT_ELEM(*vec, 0);
        ptr->b = FNT_VECT_ELEM(*vec, 0);
        ptr->f_a = 0.0;
        ptr->f_b = 0.0;
        ptr->root_x = ptr->a;
        ptr->state = done;
    } else {
        ERROR("Value (%g) is not comparable to zero.\n", value);
        return FNT_FAILURE;
    }

    return FNT_SUCCESS;
}


int method_done(void *handle) {
    if( handle == NULL )    { return FNT_FAILURE; }
    bisection_t *ptr = (bisection_t*)handle;

    /* test for completion
     *  Return FNT_DONE when comlete, or FNT_CONTINUE when not done.
     */
    if( ptr->state == initial || ptr->state == initial2 ) {
        return FNT_CONTINUE;
    }
    if( ptr->state == done ) {
        return FNT_DONE;
    }

    if( fabs(ptr->b - ptr->a) < ptr->x_tol ) {
        INFO("Upper and lower bound within termination threshold.\n");
        ptr->root_x = 0.5 * (ptr->b + ptr->a);
        ptr->state = done;
        return FNT_DONE;
    }
    if( fabs(ptr->f_b - ptr->f_a) < ptr->f_tol ) {
        INFO("Difference in function's value at upper and lower bound within termination threshold.\n");
        ptr->root_x = 0.5 * (ptr->b + ptr->a);
        ptr->state = done;
        return FNT_DONE;
    }

    return FNT_CONTINUE;
}


int method_result(void *handle, char *id, void *value_ptr) {
    if( handle == NULL )    { return FNT_FAILURE; }
    bisection_t *ptr = (bisection_t*)handle;

    /* Report any results the method produces. */
    FNT_RESULT_GET("root", id, double, ptr->root_x, value_ptr);

    ERROR("No result named '%s'.\n", id);

    return FNT_FAILURE;
}
