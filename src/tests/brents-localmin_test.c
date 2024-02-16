/*
 * de_test.c
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

/* function defined in equation 6.1 of Brent. */
double function_61(double x) {
    double sum = 0.0;
    for(int i=1; i<=20; ++i) {
        sum += pow( (2*i - 5) / (x - i*i), 2.0);
    }
    return sum;
}

int main() {

    void *fnt = NULL;

    fnt_verbose(FNT_INFO); /* request informative output */
    fnt_init(&fnt, FNT_METHODS_DIR "/methods");

    /* load Brent's methos to minimize a function */
    if( fnt_set_method(fnt, "brents-localmin", 1) == FNT_FAILURE ) {
        fprintf(stderr, "Failed to initialize method.\n");
        return 1;
    }

    /* set threshold for completion */
    double tolerance = 1e-30;
    double x_0 = 2;
    double x_1 = 3;
    double eps = 1e-6;
    double t = 1e-6;
    fnt_hparam_set(fnt, "x_0", &x_0);
    fnt_hparam_set(fnt, "x_1", &x_1);
    fnt_hparam_set(fnt, "eps", &eps);
    fnt_hparam_set(fnt, "t", &t);

    /* read and report hyper-parameters */
    fnt_hparam_get(fnt, "x_0", &x_0);
    fnt_hparam_get(fnt, "x_1", &x_1);
    fnt_hparam_get(fnt, "eps", &eps);
    fnt_hparam_get(fnt, "t", &t);
    printf("\ta: %g\n\tb: %g\n\teps: %g\n\tt: %g\n", x_0, x_1, eps, t);

    /* allocate input for objective function */
    fnt_vect_t x;
    fnt_vect_calloc(&x, 1);

    /* loop as long as method is not complete */
    while( fnt_done(fnt) == FNT_CONTINUE ) {

        /* get vector to try */
        if( fnt_next(fnt, &x) != FNT_SUCCESS ) { break; }

        /* call objective function */
        double fx = function_61(x.v[0]);

        if( fnt_verbose_level >= FNT_INFO ) {
            fnt_vect_print(&x, "f(", "%.3f");
            printf(") -> %g\n", fx);
        }

        /* update method */
        if( fnt_set_value(fnt, &x, fx) != FNT_SUCCESS ) { break; }
    }

    /* Get best result. */
    if( fnt_minimum(fnt, &x, NULL) == FNT_SUCCESS )
        fnt_vect_println(&x, "Best result: ", NULL);

    /* Get/report any results beyond best input vector. */
    fnt_result(fnt, NULL);

    /* free input vector */
    fnt_vect_free(&x);

    /* free the method */
    fnt_free(&fnt);

    return 0;
}
