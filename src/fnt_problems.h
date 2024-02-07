/*
 * fnt_problems.h
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */
#ifndef FNT_PROBLEMS_H
#define FNT_PROBLEMS_H

#include "fnt_vect.h"

/* \brief Computes Rastrigin function.
 * see: https://en.wikipedia.org/wiki/Rastrigin_function
 * Minimum is at (x_0,...,x_n) = (0,...,0).
 */
double rastrigin(fnt_vect_t *x) {
    double A = 10.0;
    double sum = 0.0;
    int n = x->n;
    for(int i=0; i<n; ++i) {
        double x_i = x->v[i];
        sum += x_i * x_i - A * cos(2*M_PI*x_i);
    }

    return A * n + sum;
}


/* \brief Computes Ackley function.
 * see: https://en.wikipedia.org/wiki/Ackley_function
 * Minimum is at (x,y) = (0,0).
 */
double ackley(double x, double y) {
    return (-20.0) * exp(-2.0 * sqrt( 0.5 * (x*x + y*y)) )
            - exp( 0.5 * (cos(2*M_PI*x) + cos(2*M_PI*y)) ) + M_E + 20;
}


/** \brief Computes Rosenbrock function.
 * see: https://en.wikipedia.org/wiki/Rosenbrock_function
 * Minimum is at (x,y) = (1,1).
 */
double rosenbrock(double x, double y) {
    const double a = 1, b =100;   /* minumum should be at (x,y) = (1,1) */
    
    double f = pow(a - x, 2.0) + b * pow(y - x*x , 2.0);

    return f;
}

#endif /* FNT_PROBLEMS_H */
