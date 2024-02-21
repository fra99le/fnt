/*
 * bisection_test.c
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
    // 3x^3 - 5x^2 - 6x + 10
    return 3*pow(x, 3.0) - 5*pow(x,2.0) - 6*x + 10;
}

int main() {

    void *fnt = NULL;

    fnt_verbose(FNT_INFO); /* request informative output */
    fnt_init(&fnt, FNT_METHODS_DIR "/methods");

    /* load trapezoidal method to find the area under a polynomial function */
    if( fnt_set_method(fnt, "trapezoidal", 1) == FNT_FAILURE ) {
        return 1;
    }

    /* display info */
    fnt_info(fnt);

    /* place initial bounds for search */
    double x_0 = 2.0;
    double x_1 = 3.0;
    int subintervals = 20;
    fnt_hparam_set(fnt, "lower", &x_0);
    fnt_hparam_set(fnt, "upper", &x_1);
    fnt_hparam_set(fnt, "subintervals", &subintervals);

    /* allocate input for objective function */
    fnt_vect_t x;
    fnt_vect_calloc(&x, 1);

    /* loop as long as method is not complete */
    while( fnt_done(fnt) == FNT_CONTINUE ) {

        /* get vector to try */
        if( fnt_next(fnt, &x) != FNT_SUCCESS ) { break; }

        /* call objective function */
        double fx = polynomial(FNT_VECT_ELEM(x, 0));

        fnt_vect_print(&x, "f(", "%.3f");
        printf(") -> %g\n", fx);

        /* update method */
        if( fnt_set_value(fnt, &x, fx) != FNT_SUCCESS ) { break; }
    }

    /* Get/report any results beyond best input vector. */
    double area = 0.0;
    fnt_result(fnt, &area);
    printf("Area under polynomial is %g\n", area);

    /* free input vector */
    fnt_vect_free(&x);

    /* free the method */
    fnt_free(&fnt);

    return 0;
}
