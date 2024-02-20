/*
 * secant.c
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../fnt.h"

#ifndef epsilon
const double epsilon = 1e-6;
#endif /* epsilon */

/* MARK: Method type definitions */

typedef enum secant_state {
    secant_initial, secant_running, secant_done
} secant_state_t;

typedef struct secant {
    secant_state_t state;

    double x_0;
    double x_1;
    double f_tol;

    double x_prev;
    double fx_prev;
    double x_next;
} secant_t;


/* MARK: Functions called by FNT */

/* \brief Provides the name of the method.
 * \param name Allocated buffer to hold the name.
 * \param size Size of the name buffer in bytes.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_name(char *name, int size) {
    if( snprintf(name, size, "secant") >= size ) {
        return FNT_FAILURE;
    }

    return FNT_SUCCESS;
}


/* \brief Initialize intenal state for method.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_init(void **handle_ptr, int dimensions) {
    if( handle_ptr == NULL )    { return FNT_FAILURE; }
    secant_t *ptr = calloc(1, sizeof(secant_t));
    if( ptr == NULL )           { return FNT_FAILURE; }
    *handle_ptr = (void*)ptr;

    if( dimensions > 1 ) {
        ERROR("Newton-Raphson is a single variate method.\n");
        return FNT_FAILURE;
    }

    /* initialize method here */
    ptr->f_tol = 1e-6;
    ptr->x_next = 0.0;

    ptr->state = secant_initial;

    return FNT_SUCCESS;
}


int method_free(void **handle_ptr) {
    if( handle_ptr == NULL )    { return FNT_FAILURE; }
    if( *handle_ptr == NULL )   { return FNT_FAILURE; }
    secant_t *ptr = (secant_t*)*handle_ptr;

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
    secant_t *ptr = (secant_t*)handle;

    FNT_HPARAM_SET("x_0", id, double, value_ptr, ptr->x_0);
    FNT_HPARAM_SET("x_1", id, double, value_ptr, ptr->x_1);
    FNT_HPARAM_SET("f_tol", id, double, value_ptr, ptr->f_tol);

    return FNT_SUCCESS;
}


int method_hparam_get(void *handle, char *id, void *value_ptr) {
    if( id == NULL )        { return FNT_FAILURE; }
    if( value_ptr == NULL ) { return FNT_FAILURE; }
    secant_t *ptr = (secant_t*)handle;

    FNT_HPARAM_GET("x_0", id, double, ptr->x_0, value_ptr);
    FNT_HPARAM_GET("x_1", id, double, ptr->x_1, value_ptr);
    FNT_HPARAM_GET("f_tol", id, double, ptr->f_tol, value_ptr);

    return FNT_SUCCESS;
}


int method_next(void *handle, fnt_vect_t *vec) {
    secant_t *ptr = (secant_t*)handle;
    if( ptr == NULL )       { return FNT_FAILURE; }
    if( vec == NULL )       { return FNT_FAILURE; }
    if( vec->v == NULL )    { return FNT_FAILURE; }

    if( ptr->state == secant_initial ) {
        FNT_VECT_ELEM(*vec, 0) = ptr->x_0;
        return FNT_SUCCESS;
    }

    /* fill vector pointed to by vec with next input to try */
    FNT_VECT_ELEM(*vec, 0) = ptr->x_next;

    return FNT_SUCCESS;
}


int method_value(void *handle, fnt_vect_t *vec, double value) {
    secant_t *ptr = (secant_t*)handle;
    if( ptr == NULL )       { return FNT_FAILURE; }
    if( vec == NULL )       { return FNT_FAILURE; }
    if( vec->v == NULL )    { return FNT_FAILURE; }

    if( ptr->state == secant_initial ) {
        ptr->x_prev = FNT_VECT_ELEM(*vec, 0);
        ptr->fx_prev = value;

        ptr->x_next = ptr->x_1;
        ptr->state = secant_running;

        return FNT_SUCCESS;
    }

    /* update method using value */
    double x = FNT_VECT_ELEM(*vec, 0);
    double fx = value;

    double x_prev = ptr->x_prev;
    double fx_prev = ptr->fx_prev;

    double delta_x = x - x_prev;
    double delta_fx = fx - fx_prev;

    if( fabs(delta_fx) < epsilon ) { return FNT_FAILURE; }

    /* Note: if delta_fx is small, error in x_next can be very large */
    ptr->x_next = x_prev - fx_prev * delta_x / delta_fx;

    ptr->x_prev = x;
    ptr->fx_prev = fx;

    if( ptr->state == secant_initial ) { ptr->state = secant_running; }

    return FNT_SUCCESS;
}


int method_value_gradient(void *handle, fnt_vect_t *vec, double value, fnt_vect_t *gradient) {
    secant_t *ptr = (secant_t*)handle;
    if( ptr == NULL )       { return FNT_FAILURE; }
    if( vec == NULL )       { return FNT_FAILURE; }

    /* update method using value and derivative/gradient */
    return method_value(handle, vec, value);
}


int method_done(void *handle) {
    secant_t *ptr = (secant_t*)handle;
    if( ptr == NULL )       { return FNT_FAILURE; }

    if( ptr->state == secant_initial ) { return FNT_CONTINUE; }

    /* test for completion
     *  Return FNT_DONE when comlete, or FNT_CONTINUE when not done.
     */
    if( ptr->f_tol > fabs(ptr->fx_prev) ) {
        ptr->state = secant_done;
        return FNT_DONE;
    }

    return FNT_CONTINUE;
}


int method_result(void *handle, void *extra) {

    /* Optional method to report any additional results if the method
     * produces such results.
     */

    return FNT_SUCCESS;
}
