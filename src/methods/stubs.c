/*
 * stub.c
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../fnt.h"


/* MARK: Method type definitions */

typedef enum stub_states {
    stub_initial, stub_running, stub_done
} stub_state_t;

typedef struct stub {
    int placeholder;
} stub_t;


/* MARK: Functions called by FNT */

/* \brief Provides the name of the method.
 * \param name Allocated buffer to hold the name.
 * \param size Size of the name buffer in bytes.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_name(char *name, int size) {
    if( strlcpy(name,"stub",size) >= size ) {
        return FNT_FAILURE;
    }

    return FNT_SUCCESS;
}


/* \brief Initialize intenal state for method.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_init(void **handle_ptr, int dimensions) {
    if( handle_ptr == NULL )    { return FNT_FAILURE; }
    stub_t *ptr = calloc(1, sizeof(stub_t));
    if( ptr == NULL )           { return FNT_FAILURE; }
    *handle_ptr = (void*)ptr;

    /* initialize method here */

    return FNT_SUCCESS;
}


int method_free(void **handle_ptr) {
    if( handle_ptr == NULL )    { return FNT_FAILURE; }
    if( *handle_ptr == NULL )   { return FNT_FAILURE; }
    stub_t *ptr = (stub_t*)*handle_ptr;

    /* free any memory allocated by method */

    free(ptr);  *handle_ptr = ptr = NULL;

    return FNT_SUCCESS;
}


/* \brief Display information about the method to the console.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_info() {
    printf(
"The stub method is not a numerical method, instead it provides a starting\n"
"point for imlementing real numerical methods.\n"
"\n"
"Hyper-parameters:\n"
"name\trequired\ttype\tDefault\tDescription\n"
"placeholder\toptional\tint\t0\tJust an example hyper-parameter.\n"
"\n"
"References:\n"
"None"
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
    stub_t *ptr = (stub_t*)handle;

    FNT_HPARAM_SET("placeholder", id, int, value_ptr, ptr->placeholder);

    return FNT_FAILURE;
}


int method_hparam_get(void *handle, char *id, void *value_ptr) {
    if( id == NULL )        { return FNT_FAILURE; }
    if( value_ptr == NULL ) { return FNT_FAILURE; }
    if( handle == NULL )    { return FNT_FAILURE; }
    stub_t *ptr = (stub_t*)handle;

    FNT_HPARAM_GET("placeholder", id, int, ptr->placeholder, value_ptr);

    return FNT_FAILURE;
}


int method_next(void *handle, fnt_vect_t *vec) {
    if( handle == NULL )    { return FNT_FAILURE; }
    if( vec == NULL )       { return FNT_FAILURE; }
    if( vec->v == NULL )    { return FNT_FAILURE; }
    stub_t *ptr = (stub_t*)handle;

    /* fill vector pointed to by vec with next input to try */

    return FNT_FAILURE;
}


int method_value(void *handle, fnt_vect_t *vec, double value) {
    if( handle == NULL )    { return FNT_FAILURE; }
    if( vec == NULL )       { return FNT_FAILURE; }
    if( vec->v == NULL )    { return FNT_FAILURE; }
    stub_t *ptr = (stub_t*)handle;

    /* update method using value */

    return FNT_FAILURE;
}


int method_value_gradient(void *handle, fnt_vect_t *vec, double value, fnt_vect_t gradient) {
    if( handle == NULL )    { return FNT_FAILURE; }
    if( vec == NULL )       { return FNT_FAILURE; }
    if( vec->v == NULL )    { return FNT_FAILURE; }
    stub_t *ptr = (stub_t*)handle;

    /* update method using value and derivative/gradient */

    return FNT_FAILURE;
}


int method_done(void *handle) {
    if( handle == NULL )    { return FNT_FAILURE; }
    stub_t *ptr = (stub_t*)handle;

    /* test for completion
     *  Return FNT_DONE when comlete, or FNT_CONTINUE when not done.
     */

    return FNT_FAILURE;
}


int method_result(void *handle, void *extra) {
    if( handle == NULL )    { return FNT_FAILURE; }
    stub_t *ptr = (stub_t*)handle;

    /* Optional method to report any additional results if the method
     * produces such results.
     */

    return FNT_FAILURE;
}