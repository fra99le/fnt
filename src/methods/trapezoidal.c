/*
 * trapezoidal.c
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../fnt.h"


/* MARK: Method type definitions */

typedef enum trapezoidal_states {
    trapezoidal_initial, trapezoidal_running, trapezoidal_done
} trapezoidal_state_t;

typedef struct trapezoidal {

    /* method state */
    trapezoidal_state_t state;
    double first_fx;
    double sum;
    double last_fx;
    int curr_subinterval;

    /* hyper-parameters */
    double x_0;
    double x_1;
    int n;

    /* result */
    double area;

} trapezoidal_t;


/* MARK: Functions called by FNT */

/* \brief Provides the name of the method.
 * \param name Allocated buffer to hold the name.
 * \param size Size of the name buffer in bytes.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_name(char *name, int size) {
    if( snprintf(name, size, "trapezoidal") >= size ) {
        return FNT_FAILURE;
    }

    return FNT_SUCCESS;
}


/* \brief Initialize intenal state for method.
 * \param handle_ptr Pointer to the method handle pointer.
 * \param dimensions Number of dimensions in the objactive function input.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_init(void **handle_ptr, int dimensions) {
    if( handle_ptr == NULL )    { return FNT_FAILURE; }
    trapezoidal_t *ptr = calloc(1, sizeof(trapezoidal_t));
    if( ptr == NULL )           { return FNT_FAILURE; }
    *handle_ptr = (void*)ptr;

    /* initialize method here */
    ptr->state = trapezoidal_initial;
    ptr->curr_subinterval = 0;

    return FNT_SUCCESS;
}


/* \brief Free any resources allocated for the method.
 * \param handle_ptr Pointer to the method handle pointer.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_free(void **handle_ptr) {
    if( handle_ptr == NULL )    { return FNT_FAILURE; }
    if( *handle_ptr == NULL )   { return FNT_FAILURE; }
    trapezoidal_t *ptr = (trapezoidal_t*)*handle_ptr;

    /* free any memory allocated by method */

    free(ptr);  *handle_ptr = ptr = NULL;

    return FNT_SUCCESS;
}


/* \brief Display information about the method to the console.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_info() {
    printf(
"The trapezoidal method is an integration method that samples the interval\n"
"being integrated at regular subintervals and used trapezoids to estimate the\n"
"area under the curve.\n"
"\n"
"Hyper-parameters:\n"
"name\trequired\ttype\tDefault\tDescription\n"
"lower\tREQUIRED\tdouble\t0.0\tLower end of the interval being integrated.\n"
"upper\tREQUIRED\tdouble\t1.0\tUpper end of the interval being integrated.\n"
"n\tREQUIRED\tint\t10\tNumber of subintervals (i.e. trapezoids) to use.\n"
"\n"
"References:\n"
"Fausett, L.V. (2002). Numerical Methods: Algorithms and Applications.\n"
"\tISBN 0-13-031400-5\n"
);
    return FNT_SUCCESS;
}


/* \brief Set any hyper-parameters needed for the method.
 * \param handle Pointer to the method handle.
 * \param id The name of the hyper-parameter.
 * \param value_ptr A pointer to the value being set.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_hparam_set(void *handle, char *id, void *value_ptr) {
    if( id == NULL )        { return FNT_FAILURE; }
    if( value_ptr == NULL ) { return FNT_FAILURE; }
    if( handle == NULL )    { return FNT_FAILURE; }
    trapezoidal_t *ptr = (trapezoidal_t*)handle;

    FNT_HPARAM_SET("lower", id, double, value_ptr, ptr->x_0);
    FNT_HPARAM_SET("upper", id, double, value_ptr, ptr->x_1);
    FNT_HPARAM_SET("subintervals", id, int, value_ptr, ptr->n);
    FNT_HPARAM_SET("n", id, int, value_ptr, ptr->n);

    ERROR("No hyper-parameter named '%s'.\n", id);

    return FNT_FAILURE;
}


/* \brief Get any hyper-parameters values form the method.
 * \param handle Pointer to the method handle.
 * \param id The name of the hyper-parameter.
 * \param value_ptr A pointer to the value being set.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_hparam_get(void *handle, char *id, void *value_ptr) {
    if( id == NULL )        { return FNT_FAILURE; }
    if( value_ptr == NULL ) { return FNT_FAILURE; }
    if( handle == NULL )    { return FNT_FAILURE; }
    trapezoidal_t *ptr = (trapezoidal_t*)handle;

    FNT_HPARAM_GET("lower", id, double, ptr->x_0, value_ptr);
    FNT_HPARAM_GET("upper", id, double, ptr->x_1, value_ptr);
    FNT_HPARAM_GET("subintervals", id, int, ptr->n, value_ptr);
    FNT_HPARAM_GET("n", id, int, ptr->n, value_ptr);

    ERROR("No hyper-parameter named '%s'.\n", id);

    return FNT_FAILURE;
}


int method_next(void *handle, fnt_vect_t *vec) {
    if( handle == NULL )    { return FNT_FAILURE; }
    if( vec == NULL )       { return FNT_FAILURE; }
    if( vec->v == NULL )    { return FNT_FAILURE; }
    trapezoidal_t *ptr = (trapezoidal_t*)handle;

    if( ptr->state == trapezoidal_done ) {
        ERROR("ERROR: Requested next value after method has finished.\n");
        return FNT_FAILURE;
    }

    /* fill vector pointed to by vec with next input to try */
    if( ptr->state == trapezoidal_initial ) {
        FNT_VECT_ELEM(*vec, 0) = ptr->x_0;
        return FNT_SUCCESS;
    }

    /* compute next x to be evaluated */
    FNT_VECT_ELEM(*vec, 0) = ptr->x_0
                                + (double)ptr->curr_subinterval
                                * (ptr->x_1 - ptr->x_0)
                                / (double)ptr->n;

    return FNT_SUCCESS;
}


int method_value(void *handle, fnt_vect_t *vec, double value) {
    if( handle == NULL )    { return FNT_FAILURE; }
    if( vec == NULL )       { return FNT_FAILURE; }
    if( vec->v == NULL )    { return FNT_FAILURE; }
    trapezoidal_t *ptr = (trapezoidal_t*)handle;

    if( ptr->state == trapezoidal_done ) {
        ERROR("Attempting to update method with a value after method completed.\n");
        return FNT_FAILURE;
    }

    /* update method using value */
    if( ptr->state == trapezoidal_initial ) {
        DEBUG("Recording first f(%g)=%g.\n", FNT_VECT_ELEM(*vec, 0), value);

        /* record first point, but no area is computable yet. */
        ptr->first_fx = value;
        ptr->sum = 0.0;
        ptr->curr_subinterval = 1;

        /* update state to running */
        ptr->state = trapezoidal_running;

        return FNT_SUCCESS;
    } else if( ptr->curr_subinterval >= ptr->n ) {
        DEBUG("Recording final f(%g)=%g and computing area.\n", FNT_VECT_ELEM(*vec, 0), value);

        /* compute final result */
        ptr->last_fx = value;
        double h = (ptr->x_1 - ptr->x_0) / (double)ptr->n;
        ptr->area = 0.5 * h * (ptr->first_fx + ptr->last_fx + 2.0 * ptr->sum);

        /* set state to done */
        ptr->state = trapezoidal_done;

        return FNT_SUCCESS;
    }

    DEBUG("Adding f(%g)=%g to sum.\n", FNT_VECT_ELEM(*vec, 0), value);

    /* add value for use in area calculation and update subinteval */
    ptr->sum += value;
    ++ptr->curr_subinterval;

    return FNT_SUCCESS;
}


int method_done(void *handle) {
    if( handle == NULL )    { return FNT_FAILURE; }
    trapezoidal_t *ptr = (trapezoidal_t*)handle;

    if( ptr->state == trapezoidal_done ) {
        return FNT_DONE;
    }

    return FNT_CONTINUE;
}


int method_result(void *handle, char *id, void *value_ptr) {
    if( handle == NULL )    { return FNT_FAILURE; }
    trapezoidal_t *ptr = (trapezoidal_t*)handle;

    if( ptr->state != trapezoidal_done ) {
        ERROR("ERROR: Request for result before method completed.\n");
        return FNT_FAILURE;
    }

    /* report the area under the function */
    FNT_RESULT_GET("area", id, double, ptr->area, value_ptr);

    ERROR("No result named '%s'.\n", id);

    return FNT_FAILURE;
}
