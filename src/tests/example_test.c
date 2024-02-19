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

    fnt_verbose(FNT_INFO); /* request informative output */
    fnt_init(&fnt, FNT_METHODS_DIR "/methods");

    /* load example for a DIM parameter problem */
    fnt_set_method(fnt, "example", DIM);

    /* display info */
    fnt_info(fnt);

    /* set hyper-parameter(s) */
    int count = 5;
    fnt_hparam_set(fnt, "count", (void*)&count);

    /* fetch a hyper-parameter */
    int norm = 0;
    fnt_hparam_get(fnt, "norm", (void*)&norm);
    printf("hyper-parameter 'norm' set to %i.\n", norm);

    /* print info, if available. */
    if( fnt_info(fnt) != FNT_SUCCESS ) {
        fprintf(stderr, "No info available.\n");
    }

    /* allocate input for objective function */
    fnt_vect_t x;
    fnt_vect_calloc(&x, DIM);

    /* loop as long as method is not complete */
    while( fnt_done(fnt) != FNT_DONE ) {

        /* get vector to try */
        if( fnt_next(fnt, &x) != FNT_SUCCESS ) { break; }

        /* call objective function */
        double fx = objective_function(x.v, DIM);

        /* update method */
        if( fnt_set_value(fnt, &x, fx) != FNT_SUCCESS ) { break; }
    }

    /* Get best result. */
    fnt_minimum(fnt, &x, NULL);

    /* Get/report any results beyond best input vector. */
    fnt_result(fnt, NULL);

    /* free the method */
    fnt_free(&fnt);

    return 0;
}
