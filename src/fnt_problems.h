/*
 * fnt_problems.h
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */
#ifndef FNT_PROBLEMS_H
#define FNT_PROBLEMS_H

#include "fnt_vect.h"

/** \brief Computes Rastrigin function.
 * see: https://en.wikipedia.org/wiki/Rastrigin_function
 * Minimum is at (x_0,...,x_n) = (0,...,0).
 * Search range should be x_i \in [-5.12,5.12].
 */
static double rastrigin(fnt_vect_t *x) {
    double A = 10.0;
    double sum = 0.0;
    int n = x->n;

    for(int i=0; i<n; ++i) {
        double x_i = x->v[i];
        sum += x_i * x_i - A * cos(2*M_PI*x_i);
    }

    return A * n + sum;
}


/** \brief Computes Ackley function.
 * see: https://en.wikipedia.org/wiki/Ackley_function
 * Minimum is at (x,y) = (0,0).
 */
static double ackley(double x, double y) {
    return (-20.0) * exp(-2.0 * sqrt( 0.5 * (x*x + y*y)) )
            - exp( 0.5 * (cos(2*M_PI*x) + cos(2*M_PI*y)) ) + M_E + 20;
}


/** \brief Computes Sphere fucntion.
 * see: https://en.wikipedia.org/wiki/Test_functions_for_optimization
 * Minimum is at (x_1,...,x_n) = (0,...,0).
 * Search range is unbounded.
 */
static double sphere(fnt_vect_t *x) {
    double sum = 0.0;
    int n = x->n;

    for(int i=0; i<n; ++i) {
        double x_i = x->v[i];
        sum += x_i * x_i;
    }

    return sum;
}


/** \brief Computes 2D Rosenbrock function.
 * see: https://en.wikipedia.org/wiki/Rosenbrock_function
 * Minimum is at (x,y) = (1,1).
 */
static double rosenbrock_2d(double x, double y) {
    const double a = 1, b =100;

    double f = pow(a - x, 2.0) + b * pow(y - x*x , 2.0);

    return f;
}


/** \brief Computes Rosenbrock function in arbitrary dimensions.
 * see: https://en.wikipedia.org/wiki/Rosenbrock_function
 * Minimum is at (x_1,...,x_n) = (1,...,1).
 */
static double rosenbrock(fnt_vect_t *x) {
    const double a = 1, b =100;
    int n = x->n;
    double sum = 0.0;

    for(int i=0; i<(n-1); ++i) {
        double x_i = x->v[i];
        double x_ip1 = x->v[i+1];
        sum += b * (x_ip1 - x_i*x_i) + pow(a - x_i, 2.0);
    }

    return sum;
}


/** \brief Computes Beale function.
 * see: https://en.wikipedia.org/wiki/Test_functions_for_optimization
 * Minimumm is at (x,y) = (3.0,0.5)
 */
static double beale(double x, double y) {
    return pow(1.5 - x + x*y, 2.0)
            + pow(2.25 - x + x*y*y, 2.0)
            + pow(2.625 - x + x*y*y*y, 2.0);
}


#endif /* FNT_PROBLEMS_H */
