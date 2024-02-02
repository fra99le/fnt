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
double rosenbrock(double x, double y) {
    /* Rosenbrock function
     * see: https://en.wikipedia.org/wiki/Rosenbrock_function
     */
    const double a = 1, b =100;   /* minumum should be at (x,y) = (1,1) */
    
    double f = pow(a - x, 2.0) + b * pow(y - x*x , 2.0);

    return f;
}

int main() {
    void *fnt = NULL;

    fnt_verbose(1); /* request informative output */
    fnt_init(&fnt, FNT_METHODS_DIR "/methods");

    /* load nelder-mead to minimize Rosenbrock function */
    if( fnt_set_method(fnt, "nelder-mead", 2) == FNT_FAILURE ) {
        return 1;
    }

    /* set threshold for completion */
    double tolerance = 1e-5;
    fnt_hparam_set(fnt, "tolerence", &tolerance);

    /* read an report default hyper-parameters */
    double alpha, beta, gamma, delta;
    fnt_hparam_get(fnt, "alpha", &alpha);
    fnt_hparam_get(fnt, "beta", &beta);
    fnt_hparam_get(fnt, "gamma", &gamma);
    fnt_hparam_get(fnt, "delta", &delta);
    printf("\talpha: %g\n\tbeta: %g\n\tgamma: %g\n\tdelta: %g\n", alpha, beta, gamma, delta);

    /* allocate input for objective function */
    fnt_vect_t x;
    fnt_vect_calloc(&x, 2);

    /* provide an initial guess */
    x.v[0] = 0.0;
    x.v[1] = 0.0;
    fnt_seed(fnt, &x);

    /* loop as long as method is not complete */
    while( fnt_done(fnt) == FNT_CONTINUE ) {

        /* get vector to try */
        if( fnt_next(fnt, &x) != FNT_SUCCESS ) { break; }


        /* call objective function */
        double fx = rosenbrock(x.v[0], x.v[1]);

        fnt_vect_print(&x, "f(", "%.3f");
        printf(") -> %g\n", fx);

        /* update method */
        if( fnt_set_value(fnt, &x, fx) != FNT_SUCCESS ) { break; }
    }

    /* Get best result. */
    if( fnt_best(fnt, &x) == FNT_SUCCESS )
        fnt_vect_println(&x, "Best result: ", "%.3f");

    /* Get/report any results beyond best input vector. */
    fnt_result(fnt, NULL);

    /* free input vector */
    fnt_vect_free(&x);

    /* free the method */
    fnt_free(&fnt);

    return 0;
}
