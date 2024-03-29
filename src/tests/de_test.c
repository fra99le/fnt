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

int main() {

    void *fnt = NULL;

    fnt_verbose(FNT_INFO); /* request informative output */
    fnt_init(&fnt, FNT_METHODS_DIR "/methods");

    /* load de to minimize Rosenbrock function */
    if( fnt_set_method(fnt, "differential evolution", 2) == FNT_FAILURE ) {
        fprintf(stderr, "Failed to initialize method.\n");
        return 1;
    }

    /* display info */
    fnt_info(fnt);

    /* set threshold for completion */
    int iterations = 1e4;
    int NP = 20;
    fnt_hparam_set(fnt, "iters", &iterations);
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
    fnt_hparam_get(fnt, "iterations", &iterations);
    fnt_hparam_get(fnt, "F", &F);
    fnt_hparam_get(fnt, "lambda", &lambda);
    fnt_hparam_get(fnt, "NP", &NP);
    printf("\titerations: %d\n\tF: %g\n\tlambda: %g\n\tNP: %d\n", iterations, F, lambda, NP);

    /* allocate input for objective function */
    fnt_vect_t x;
    fnt_vect_calloc(&x, 2);

    /* loop as long as method is not complete */
    while( fnt_done(fnt) == FNT_CONTINUE ) {

        /* get vector to try */
        if( fnt_next(fnt, &x) != FNT_SUCCESS ) { break; }

        /* call objective function */
        double fx = ackley(FNT_VECT_ELEM(x, 0), FNT_VECT_ELEM(x, 1));

        if( fnt_verbose_level >= FNT_INFO ) {
            fnt_vect_print(&x, "f(", "%.3f");
            printf(") -> %g\n", fx);
        }

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
