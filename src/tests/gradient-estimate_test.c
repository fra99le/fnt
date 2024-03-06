/*
 * gradient-estimate_test.c
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

/* Example 3 from section 16.6 of Anton (page 1032) */
double example3(double x, double y) {
    return 3.0 * (x*x) * y;
}

int main() {

    void *fnt = NULL;

    fnt_verbose(FNT_INFO); /* request informative output */
    fnt_init(&fnt, FNT_METHODS_DIR "/methods");

    /* load gradient-estimate to minimize Rosenbrock function */
    if( fnt_set_method(fnt, "gradient estimate", 2) == FNT_FAILURE ) {
        return 1;
    }

    /* display info */
    fnt_info(fnt);

    /* set hyper-parameters */
    double step = 1e-4;
    fnt_hparam_set(fnt, "step", &step);

    #if 0
    /* use a vector of steps instead */
    fnt_vect_t steps;
    fnt_vect_calloc(&steps, 2);
    FNT_VECT_ELEM(steps, 0) = 1e-6;
    FNT_VECT_ELEM(steps, 1) = 1e-3;
    fnt_hparam_set(fnt, "step_vec", &steps);
    #endif

    fnt_vect_t x0;
    fnt_vect_calloc(&x0, 2);
    FNT_VECT_ELEM(x0, 0) = 1.0;
    FNT_VECT_ELEM(x0, 1) = 2.0;
    fnt_hparam_set(fnt, "x0", &x0);

    /* allocate input for objective function */
    fnt_vect_t x;
    fnt_vect_calloc(&x, 2);

    /* loop as long as method is not complete */
    while( fnt_done(fnt) == FNT_CONTINUE ) {

        /* get vector to try */
        if( fnt_next(fnt, &x) != FNT_SUCCESS ) { break; }

        /* call objective function */
        double fx = example3(FNT_VECT_ELEM(x, 0), FNT_VECT_ELEM(x, 1));

        fnt_vect_print(&x, "f(", "%.4f");
        printf(") -> %g\n", fx);

        /* update method */
        if( fnt_set_value(fnt, &x, fx) != FNT_SUCCESS ) { break; }
    }

    /* Get best result. */
    double min_fx;
    if( fnt_result(fnt, "gradient", &x) == FNT_SUCCESS ) {
        fnt_vect_print(&x0, "Gradient of f at ", NULL);
        fnt_vect_print(&x, " is ", NULL);
        printf(".\n", min_fx);
    }

    /* free input vector */
    fnt_vect_free(&x);
    fnt_vect_free(&x0);
    #if 0
    fnt_vect_free(&steps);
    #endif

    /* free the method */
    fnt_free(&fnt);

    return 0;
}
