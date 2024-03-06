/*
 * gradient-estimate.c
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../fnt.h"


/* MARK: Method type definitions */

typedef enum gradient_est_states {
    gradient_est_initial, gradient_est_running, gradient_est_done
} gradient_est_state_t;

typedef struct gradient_est {

    /* hyper-parameters */
    fnt_vect_t x0;
    double step;
    fnt_vect_t steps;
    int has_steps_vec;

    /* method state */
    gradient_est_state_t state;
    double fx0;
    int curr;

    /* results */
    fnt_vect_t gradient;

} gradient_est_t;


/* MARK: Functions called by FNT */

/* \brief Provides the name of the method.
 * \param name Allocated buffer to hold the name.
 * \param size Size of the name buffer in bytes.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_name(char *name, int size) {
    if( snprintf(name, size, "gradient estimate") >= size ) {
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
    gradient_est_t *ptr = calloc(1, sizeof(gradient_est_t));
    if( ptr == NULL )           { return FNT_FAILURE; }
    *handle_ptr = (void*)ptr;

    /* initialize method here */
    ptr->state = gradient_est_initial;

    /* allocate x0, steps and result vectors */
    fnt_vect_calloc(&ptr->x0, dimensions);
    fnt_vect_calloc(&ptr->steps, dimensions);
    fnt_vect_calloc(&ptr->gradient, dimensions);

    /* set default step size */
    ptr->step = 1e-3;
    for(int i=0; i<dimensions; ++i) {
        FNT_VECT_ELEM(ptr->steps, i) = ptr->step;
    }

    return FNT_SUCCESS;
}


/* \brief Free any resources allocated for the method.
 * \param handle_ptr Pointer to the method handle pointer.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_free(void **handle_ptr) {
    if( handle_ptr == NULL )    { return FNT_FAILURE; }
    if( *handle_ptr == NULL )   { return FNT_FAILURE; }
    gradient_est_t *ptr = (gradient_est_t*)*handle_ptr;

    /* free any memory allocated by method */
    fnt_vect_free(&ptr->x0);
    fnt_vect_free(&ptr->steps);
    fnt_vect_free(&ptr->gradient);

    free(ptr);  *handle_ptr = ptr = NULL;

    return FNT_SUCCESS;
}


/* \brief Display information about the method to the console.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_info() {
    printf(
"The gradient estimation method uses small steps in each dimension to\n"
"estimate the gradient of a fucntion at a specified point.\n"
"\n"
"Hyper-parameters:\n"
"name\t\trequired\ttype\t\tDefault\tDescription\n"
"x0\t\tREQUIRED\tfnt_vect_t\tnone\tPoint where the gradient is estimated.\n"
"step\t\toptional\tdouble\t\t1e-3\tStep size to use.\n"
"step_vec\toptional\tfnt_vect_t\tnone\tStep sizes to use in each dimension.\n"
"\n"
"Results:\n"
"name\t\ttype\tDescription\n"
"gradient\tdouble\tEstimated gradient at x0.\n"
"\n"
"References:\n"
"Anton, H. (1992). Calculus with analytic geometry -- 4th ed.\n"
"\tISBN 0-471-50901-9\n"
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
    gradient_est_t *ptr = (gradient_est_t*)handle;

    FNT_HPARAM_SET("step", id, double, value_ptr, ptr->step);
    FNT_HPARAM_SET_VECT("x0", id, value_ptr, &ptr->x0);

    if( strncmp("step_vec", id, 9) == 0 ) {
        ptr->has_steps_vec = 1;
        FNT_HPARAM_SET_VECT("step_vec", id, value_ptr, &ptr->steps);
        /* Note: FNT_HPARAM_SET_VECT will return FNT_SUCCESS */
    }

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
    gradient_est_t *ptr = (gradient_est_t*)handle;
    if( ptr == NULL )       { return FNT_FAILURE; }
    if( id == NULL )        { return FNT_FAILURE; }
    if( value_ptr == NULL ) { return FNT_FAILURE; }

    FNT_HPARAM_GET("step", id, double, ptr->step, value_ptr);
    FNT_HPARAM_GET_VECT("x0", id, &ptr->x0, value_ptr);

    ERROR("No hyper-parameter named '%s'.\n", id);

    return FNT_FAILURE;
}


int method_next(void *handle, fnt_vect_t *vec) {
    gradient_est_t *ptr = (gradient_est_t*)handle;
    if( ptr == NULL )       { return FNT_FAILURE; }
    if( vec == NULL )       { return FNT_FAILURE; }
    if( vec->v == NULL )    { return FNT_FAILURE; }

    /* fill vector pointed to by vec with next input to try */
    if( ptr->state == gradient_est_initial ) {
        return fnt_vect_copy(vec, &ptr->x0);
    }

    fnt_vect_copy(vec, &ptr->x0);
    if( ptr->has_steps_vec ) {
        DEBUG("DEBUG: Updating x0 with element %i of step vector (%g).\n", ptr->curr, FNT_VECT_ELEM(ptr->steps, ptr->curr));
        FNT_VECT_ELEM(*vec, ptr->curr) += FNT_VECT_ELEM(ptr->steps, ptr->curr);
    } else {
        DEBUG("DEBUG: Updating x0 with step (%g).\n", ptr->step);
        FNT_VECT_ELEM(*vec, ptr->curr) += ptr->step;
    }

    return FNT_SUCCESS;
}


int method_value(void *handle, fnt_vect_t *vec, double value) {
    gradient_est_t *ptr = (gradient_est_t*)handle;
    if( ptr == NULL )       { return FNT_FAILURE; }
    if( vec == NULL )       { return FNT_FAILURE; }
    if( vec->v == NULL )    { return FNT_FAILURE; }

    /* update method using value */
    if( ptr->state == gradient_est_initial ) {
        ptr->fx0 = value;
        ptr->curr = 0;
        ptr->state = gradient_est_running;
        return FNT_SUCCESS;
    }

    /* estimate partial derivative with respect to current dimension */
    double step = ptr->step;
    if( ptr->has_steps_vec ) {
        step = FNT_VECT_ELEM(ptr->steps, ptr->curr);
    }
    FNT_VECT_ELEM(ptr->gradient, ptr->curr) = (value-ptr->fx0) / step;


    /* update state */
    ++ptr->curr;
    if( ptr->curr >= ptr->gradient.n ) {
        ptr->state = gradient_est_done;
    }

    return FNT_SUCCESS;
}


int method_done(void *handle) {
    gradient_est_t *ptr = (gradient_est_t*)handle;
    if( ptr == NULL )       { return FNT_FAILURE; }

    /* test for completion
     *  Return FNT_DONE when comlete, or FNT_CONTINUE when not done.
     */
    if( ptr->state == gradient_est_done )
        return FNT_DONE;

    return FNT_CONTINUE;
}


int method_result(void *handle, char *id, void *value_ptr) {
    gradient_est_t *ptr = (gradient_est_t*)handle;
    if( ptr == NULL )       { return FNT_FAILURE; }

    /* Report any results the method produces. */
    FNT_RESULT_GET_VECT("gradient", id, ptr->gradient, value_ptr);

    ERROR("No result named '%s'.\n", id);

    return FNT_FAILURE;
}
