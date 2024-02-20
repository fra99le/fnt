/*
 * newton-raphson_test.c
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

double derivative(double x) {
    // 9x^2 - 10x - 6
    return 9*pow(x, 2.0) - 10*x - 6;
}

int main() {

    void *fnt = NULL;

    fnt_verbose(FNT_INFO); /* request informative output */
    fnt_init(&fnt, FNT_METHODS_DIR "/methods");

    /* load newton-raphsom to find the root of the polynomial function */
    if( fnt_set_method(fnt, "newton-raphson", 1) == FNT_FAILURE ) {
        return 1;
    }

    /* display info */
    fnt_info(fnt);

    /* set threshold for completion */
    double f_tol = 1e-5;
    fnt_hparam_set(fnt, "f_tol", &f_tol);
    double x_0 = 2.0;
    fnt_hparam_set(fnt, "x_0", &x_0);

    /* allocate input for objective function */
    fnt_vect_t x, dfdx;
    fnt_vect_calloc(&x, 1);
    fnt_vect_calloc(&dfdx, 1);

    /* loop as long as method is not complete */
    while( fnt_done(fnt) == FNT_CONTINUE ) {

        /* get vector to try */
        if( fnt_next(fnt, &x) != FNT_SUCCESS ) { break; }

        /* call objective function */
        double fx = polynomial(FNT_VECT_ELEM(x, 0));
        FNT_VECT_ELEM(dfdx, 0) = derivative(FNT_VECT_ELEM(x, 0));

        fnt_vect_print(&x, "f(", "%.3f");
        printf(") -> %g\tf'(x) -> %g\n", fx, FNT_VECT_ELEM(dfdx, 0));

        /* update method */
        if( fnt_set_value_gradient(fnt, &x, fx, &dfdx) != FNT_SUCCESS ) { break; }
    }

    /* Get best result. */
    if( fnt_root(fnt, &x, NULL) == FNT_SUCCESS )
        fnt_vect_println(&x, "Best result: ", "%.3f");

    /* Get/report any results beyond best input vector. */
    fnt_result(fnt, NULL);

    /* free input vector */
    fnt_vect_free(&x);

    /* free the method */
    fnt_free(&fnt);

    return 0;
}
