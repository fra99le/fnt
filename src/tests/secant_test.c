/*
 * secant_test.c
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

double polynomial(double x) {
    // 3x^3 - 5x^2 - 6x + 5
    return 3*pow(x, 3.0) - 5*pow(x,2.0) - 6*x + 5;
}

int main() {

    void *fnt = NULL;

    fnt_verbose(FNT_DEBUG); /* request informative output */
    fnt_init(&fnt, FNT_METHODS_DIR "/methods");

    /* load nelder-mead to minimize Rosenbrock function */
    if( fnt_set_method(fnt, "secant", 1) == FNT_FAILURE ) {
        return 1;
    }

    /* set threshold for completion */
    double f_tol = 1e-5;
    fnt_hparam_set(fnt, "f_tol", &f_tol);
    double x_0 = 2.0;
    double x_1 = 3.0;
    fnt_hparam_set(fnt, "x_0", &x_0);
    fnt_hparam_set(fnt, "x_1", &x_1);

    /* allocate input for objective function */
    fnt_vect_t x;
    fnt_vect_calloc(&x, 1);

    /* provide an initial guess */
    x.v[0] = 0.0;
    fnt_seed(fnt, &x);

    /* loop as long as method is not complete */
    while( fnt_done(fnt) == FNT_CONTINUE ) {

        /* get vector to try */
        if( fnt_next(fnt, &x) != FNT_SUCCESS ) { break; }

        /* call objective function */
        double fx = polynomial(x.v[0]);

        fnt_vect_print(&x, "f(", "%.3f");
        printf(") -> %g\n", fx);

        /* update method */
        if( fnt_set_value(fnt, &x, fx) != FNT_SUCCESS ) { break; }
    }

    /* Get best result. */
    double x_root;
    if( fnt_root(fnt, &x, &x_root) == FNT_SUCCESS ) {
        fnt_vect_println(&x, "Best result: ", "%.3f");
    }

    /* Get/report any results beyond best input vector. */
    fnt_result(fnt, NULL);

    /* free input vector */
    fnt_vect_free(&x);

    /* free the method */
    fnt_free(&fnt);

    return 0;
}