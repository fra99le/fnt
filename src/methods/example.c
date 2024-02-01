/*
 * example.c
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */
#include <stdio.h>
#include <string.h>
#include "../fnt.h"

/* This example computes the L2 norms of each input and returns the sum of
 * those norms.  Since this method does not request additional function
 * evaluations, only the seeded values will be used.  */

/* MARK: Method internals */

static double sum = 0.0;

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
int method_init() {
    sum = 0.0;
    return FNT_SUCCESS;
}

int method_info() {
    printf("This example computes the L2 norms of each input and returns the"
           "sum of those norms.  Since this method does not request additional"
           "function evaluations, only the seeded values will be used.\n");
    return FNT_SUCCESS;
}

/* \brief Set any hyper-parameters needed for the method.
 * \param 
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_hparams(double *doubles, void *misc) {
    return FNT_FAILURE;
}

int method_next(double *vec) {
    return FNT_FAILURE;
}

int method_value(double *vec, double value) {
    return FNT_FAILURE;
}

int method_done() {
    return FNT_FAILURE;
}
