/*
 * newton-raphson.c
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

typedef enum nr_state {
    nr_initial, nr_running, nr_done
} nr_state_t;

typedef struct newton_raphson {

    /* method internal state */
    nr_state_t state;
    double last_x;
    double last_fx;
    double next_x;

    /* hyper-parameters */
    double f_tol;

    /* result */
    double root_x;

} newton_raphson_t;


/* MARK: Functions called by FNT */

/* \brief Provides the name of the method.
 * \param name Allocated buffer to hold the name.
 * \param size Size of the name buffer in bytes.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_name(char *name, int size) {
    if( snprintf(name, size, "newton-raphson") >= size ) {
        return FNT_FAILURE;
    }

    return FNT_SUCCESS;
}


/* \brief Initialize intenal state for method.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_init(void **handle_ptr, int dimensions) {
    if( handle_ptr == NULL )    { return FNT_FAILURE; }
    newton_raphson_t *ptr = calloc(1, sizeof(newton_raphson_t));
    if( ptr == NULL )           { return FNT_FAILURE; }
    *handle_ptr = (void*)ptr;

    if( dimensions > 1 ) {
        ERROR("Newton-Raphson is a single variate method.\n");
        return FNT_FAILURE;
    }

    /* initialize method here */
    ptr->f_tol = 1e-6;
    ptr->next_x = 0.0;

    ptr->state = nr_initial;

    return FNT_SUCCESS;
}


int method_free(void **handle_ptr) {
    if( handle_ptr == NULL )    { return FNT_FAILURE; }
    if( *handle_ptr == NULL )   { return FNT_FAILURE; }
    newton_raphson_t *ptr = (newton_raphson_t*)*handle_ptr;

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
    newton_raphson_t *ptr = (newton_raphson_t*)handle;

    FNT_HPARAM_SET("x_0", id, double, value_ptr, ptr->next_x);
    FNT_HPARAM_SET("f_tol", id, double, value_ptr, ptr->f_tol);

    return FNT_SUCCESS;
}


int method_hparam_get(void *handle, char *id, void *value_ptr) {
    if( id == NULL )        { return FNT_FAILURE; }
    if( value_ptr == NULL ) { return FNT_FAILURE; }
    newton_raphson_t *ptr = (newton_raphson_t*)handle;

    FNT_HPARAM_GET("x_0", id, double, ptr->next_x, value_ptr);
    FNT_HPARAM_GET("f_tol", id, double, ptr->f_tol, value_ptr);

    return FNT_SUCCESS;
}


int method_next(void *handle, fnt_vect_t *vec) {
    newton_raphson_t *ptr = (newton_raphson_t*)handle;
    if( ptr == NULL )       { return FNT_FAILURE; }
    if( vec == NULL )       { return FNT_FAILURE; }

    /* fill vector pointed to by vec with next input to try */
    FNT_VECT_ELEM(*vec, 0) = ptr->next_x;

    return FNT_SUCCESS;
}


int method_value(void *handle, fnt_vect_t *vec, double value) {

    /* update method using value */
    ERROR("ERROR: Newton-Raphsom method requires a dervative.\n");

    return FNT_FAILURE;
}


int method_value_gradient(void *handle, fnt_vect_t *vec, double value, fnt_vect_t *gradient) {
    newton_raphson_t *ptr = (newton_raphson_t*)handle;
    if( ptr == NULL )       { return FNT_FAILURE; }
    if( vec == NULL )       { return FNT_FAILURE; }
    if( gradient == NULL )  { return FNT_FAILURE; }

    /* update method using value and derivative/gradient */
    double x = FNT_VECT_ELEM(*vec, 0);
    double fx = value;
    double fx_der = FNT_VECT_ELEM(*gradient, 0);

    if( fabs(fx_der) < epsilon ) { return FNT_FAILURE; }

    ptr->last_x = x;
    ptr->last_fx = fx;
    ptr->next_x = x - fx / fx_der;

    if( ptr->state == nr_initial ) { ptr->state = nr_running; }

    return FNT_SUCCESS;
}


int method_done(void *handle) {
    newton_raphson_t *ptr = (newton_raphson_t*)handle;
    if( ptr == NULL )       { return FNT_FAILURE; }

    if( ptr->state == nr_initial ) { return FNT_CONTINUE; }

    /* test for completion
     *  Return FNT_DONE when comlete, or FNT_CONTINUE when not done.
     */
    if( ptr->f_tol > fabs(ptr->last_fx) ) {

        /* record result */
        ptr->root_x = ptr->last_x;

        return FNT_DONE;
    }

    return FNT_CONTINUE;
}


int method_result(void *handle, char *id, void *value_ptr) {
    if( handle == NULL )    { return FNT_FAILURE; }
    newton_raphson_t *ptr = (newton_raphson_t*)handle;

    FNT_RESULT_GET("root", id, double, ptr->root_x, value_ptr);

    return FNT_SUCCESS;
}
