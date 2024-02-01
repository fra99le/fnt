/*
 * example_test.c
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../fnt.h"

#ifndef FNT_METHODS_DIR
#define FNT_METHODS_DIR "."
#endif /* FNT_METHODS_DIR */

/** \brief Computes the L_2 norm of the input vector. */
double objective_function(double *x, int len) {
    int sum = 0;
    for(int i=0; i<len; ++i)
        sum += x[i] * x[i];
    return sqrt(sum);
}

#define DIM 3

int main() {
    void *fnt = NULL;

    fnt_verbose(2); /* request informative output */
    fnt_init(&fnt, FNT_METHODS_DIR "/methods");

    /* load example for a DIM parameter problem */
    fnt_set_method(fnt, "example", DIM);

    /* print info, if available. */
    if( fnt_info(fnt) != FNT_SUCCESS ) {
        fprintf(stderr, "No info available.\n");
    }

    double *x = calloc(DIM, sizeof(double));
    while( fnt_done(fnt) != FNT_DONE ) {
        /* get vector to try */
        if( fnt_next(fnt, x) != FNT_SUCCESS ) { break; }

        /* call objective function */
        double fx = objective_function(x, DIM);

        /* update method */
        if( fnt_set_value(fnt, x, fx) != FNT_SUCCESS ) { break; }
    }

    fnt_result(fnt);

    fnt_free(&fnt);

    return 0;
}
