/*
 * nelder-mead_test.c
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../fnt.h"
#include "../fnt_problems.h"

#ifndef FNT_METHODS_DIR
#define FNT_METHODS_DIR "."
#endif /* FNT_METHODS_DIR */

int main() {

    void *fnt = NULL;

    fnt_verbose(FNT_INFO); /* request informative output */
    fnt_init(&fnt, FNT_METHODS_DIR "/methods");

    /* load nelder-mead to minimize Rosenbrock function */
    if( fnt_set_method(fnt, "nelder-mead", 2) == FNT_FAILURE ) {
        return 1;
    }

    /* display info */
    fnt_info(fnt);

    /* read and report default hyper-parameters */
    double alpha, beta, gamma, delta;
    fnt_hparam_get(fnt, "alpha", &alpha);
    fnt_hparam_get(fnt, "beta", &beta);
    fnt_hparam_get(fnt, "gamma", &gamma);
    fnt_hparam_get(fnt, "delta", &delta);
    printf("\talpha: %g\n\tbeta: %g\n\tgamma: %g\n\tdelta: %g\n", alpha, beta, gamma, delta);

    /* allocate input for objective function */
    fnt_vect_t x;
    fnt_vect_calloc(&x, 2);

    /* loop as long as method is not complete */
    while( fnt_done(fnt) == FNT_CONTINUE ) {

        /* get vector to try */
        if( fnt_next(fnt, &x) != FNT_SUCCESS ) { break; }


        /* call objective function */
        double fx = rosenbrock_2d(FNT_VECT_ELEM(x, 0), FNT_VECT_ELEM(x, 1));

        fnt_vect_print(&x, "f(", "%.3f");
        printf(") -> %g\n", fx);

        /* update method */
        if( fnt_set_value(fnt, &x, fx) != FNT_SUCCESS ) { break; }
    }

    /* Get best result. */
    double min_fx;
    if( fnt_result(fnt, "minimum x", &x) == FNT_SUCCESS
        && fnt_result(fnt, "minimum f", &min_fx) == FNT_SUCCESS ) {
        fnt_vect_print(&x, "Minimum found at f(", NULL);
        printf(") = %g\n", min_fx);
    }

    /* free input vector */
    fnt_vect_free(&x);

    /* free the method */
    fnt_free(&fnt);

    return 0;
}
