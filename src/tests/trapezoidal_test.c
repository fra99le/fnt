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

/* Example from page 427 of Fausett. */
double one_over_x(double x) {
    // 1/x
    return 1.0 / x;
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
    double x_0 = 1.0;
    double x_1 = 2.0;
    int subintervals = 8;
    fnt_hparam_set(fnt, "lower", &x_0);
    fnt_hparam_set(fnt, "upper", &x_1);
    fnt_hparam_set(fnt, "n", &subintervals);

    /* allocate input for objective function */
    fnt_vect_t x;
    fnt_vect_calloc(&x, 1);

    /* loop as long as method is not complete */
    while( fnt_done(fnt) == FNT_CONTINUE ) {

        /* get vector to try */
        if( fnt_next(fnt, &x) != FNT_SUCCESS ) { break; }

        /* call objective function */
        double fx = one_over_x(FNT_VECT_ELEM(x, 0));

        fnt_vect_print(&x, "f(", "%.3f");
        printf(") -> %g\n", fx);

        /* update method */
        if( fnt_set_value(fnt, &x, fx) != FNT_SUCCESS ) { break; }
    }

    /* Get/report any results beyond best input vector. */
    double area = 0.0;
    if( fnt_result(fnt, "area", &area) == FNT_SUCCESS ) {
        printf("Area under function is %g\n", area);
    }

    /* free input vector */
    fnt_vect_free(&x);

    /* free the method */
    fnt_free(&fnt);

    return 0;
}
