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

    fnt_verbose(FNT_INFO); /* request informative output */
    fnt_verbose(FNT_NONE); /* request no output */
    fnt_init(&fnt, FNT_METHODS_DIR "/methods");

    /* load de to minimize Rosenbrock function */
    if( fnt_set_method(fnt, "differential evolution", 2) == FNT_FAILURE ) {
        fprintf(stderr, "Failed to initialize method.\n");
        return 1;
    }

    /* set threshold for completion */
    double tolerance = 1e-30;
    int NP = 20;
    fnt_hparam_set(fnt, "f_tol", &tolerance);
    fnt_hparam_set(fnt, "NP", &NP);

    #if 0
    /* enable this code to set an initial guess */
    fnt_vect_t start;
    fnt_vect_calloc(&start, 2);
    start.v[0] = 2.0;
    start.v[1] = 2.0;
    fnt_hparam_set(fnt, "start", &start);
    #endif /* 1 */

    #if 0
    /* enable this code to set upper and lower bounds for search */
    fnt_vect_t lower, upper;
    fnt_vect_calloc(&lower, 2);
    fnt_vect_calloc(&upper, 2);
    lower.v[0] = -10;
    lower.v[1] = -10;
    upper.v[0] = 10;
    upper.v[1] = 10;
    fnt_hparam_set(fnt, "lower", &lower);
    fnt_hparam_set(fnt, "upper", &upper);
    #endif /* 1 */

    /* read and report default hyper-parameters */
    double F, lambda;
    fnt_hparam_get(fnt, "f_tol", &tolerance);
    fnt_hparam_get(fnt, "F", &F);
    fnt_hparam_get(fnt, "lambda", &lambda);
    fnt_hparam_get(fnt, "NP", &NP);
    printf("\tf_tol: %g\n\tF: %g\n\tlambda: %g\n\tNP: %d\n", tolerance, F, lambda, NP);

    /* allocate input for objective function */
    fnt_vect_t x;
    fnt_vect_calloc(&x, 2);

    /* provide an initial guess */
    x.v[0] = 1.0;
    x.v[1] = 1.0;
    fnt_seed(fnt, &x);

    /* loop as long as method is not complete */
    while( fnt_done(fnt) == FNT_CONTINUE ) {

        /* get vector to try */
        if( fnt_next(fnt, &x) != FNT_SUCCESS ) { break; }

        /* call objective function */
        double fx = rosenbrock(x.v[0], x.v[1]);

        if( fnt_verbose_level >= FNT_INFO ) {
            fnt_vect_print(&x, "f(", "%.3f");
            printf(") -> %g\n", fx);
        }

        /* update method */
        if( fnt_set_value(fnt, &x, fx) != FNT_SUCCESS ) { break; }
    }

    /* Get best result. */
    if( fnt_best(fnt, &x) == FNT_SUCCESS )
        fnt_vect_println(&x, "Best result: ", NULL);

    /* Get/report any results beyond best input vector. */
    fnt_result(fnt, NULL);

    /* free input vector */
    fnt_vect_free(&x);

    /* free the method */
    fnt_free(&fnt);

    return 0;
}
