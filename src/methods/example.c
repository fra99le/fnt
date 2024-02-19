/*
 * example.c
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../fnt.h"

/* This example computes the L2 norms of each input and returns the sum of
 * those norms.  Since this method does not request additional function
 * evaluations, only the seeded values will be used.  */


/* MARK: Method type definitions */

typedef struct example {
    int count;
    int norm;

    double sum;
    int counter;
} example_t;


/* MARK: Functions called by FNT */

/* \brief Provides the name of the method.
 * \param name Allocated buffer to hold the name.
 * \param size Size of the name buffer in bytes.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_name(char *name, int size) {
    if( strlcpy(name,"example",size) >= size ) {
        return FNT_FAILURE;
    }

    return FNT_SUCCESS;
}


/* \brief Initialize intenal state for method.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_init(void **handle_ptr, int dimensions) {
    example_t *ptr = calloc(1, sizeof(example_t));

    ptr->count = 0;
    ptr-> norm = 0;

    ptr->sum = 0.0;
    ptr->counter = 0;

    *handle_ptr = (void*)ptr;

    return FNT_SUCCESS;
}


int method_free(void **handle_ptr) {
    if( handle_ptr == NULL )    { return FNT_FAILURE; }
    if( *handle_ptr == NULL )   { return FNT_FAILURE; }
    example_t *ptr = (example_t*)*handle_ptr;

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
    example_t *ptr = (example_t*)handle;

    if( strncmp("count", id, 5) == 0 ) {
        ptr->count = (int)(*((int*)value_ptr));
        return FNT_SUCCESS;
    }
    if( strncmp("norm", id, 4) == 0 ) {
        ptr->norm = (int)(*((int*)value_ptr));
        return FNT_SUCCESS;
    }

    return FNT_FAILURE;
}


int method_hparam_get(void *handle, char *id, void *value_ptr) {
    if( id == NULL )        { return FNT_FAILURE; }
    if( value_ptr == NULL ) { return FNT_FAILURE; }
    example_t *ptr = (example_t*)handle;

    if( strncmp("count", id, 5) == 0 ) {
        *((int*)value_ptr) = ptr->count;
        return FNT_SUCCESS;
    }
    if( strncmp("norm", id, 4) == 0 ) {
        (*((int*)value_ptr)) = ptr->norm;
        return FNT_SUCCESS;
    }

    return FNT_FAILURE;
}


int method_next(void *handle, fnt_vect_t *vec) {
    return FNT_FAILURE;
}


int method_value(void *handle, fnt_vect_t *vec, double value) {
    return FNT_FAILURE;
}


int method_done(void *handle) {
    return FNT_FAILURE;
}


int method_result(void *handle, void *extra) {
    return FNT_FAILURE;
}
