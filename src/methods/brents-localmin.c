/*
 * brent.c
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../fnt.h"


/* MARK: Method type definitions */

typedef enum brent_states {
    brent_initial, brent_starting, brent_running, brent_done
} brent_state_t;

typedef struct brent {
    brent_state_t state;

    /* variables that need to persists between function calls */
    double a;
    double b;
    double c;

    double u;
    double v;
    double w;
    double x;
    double fu;
    double fv;
    double fw;
    double fx;

    double e;
    double d;

    /* hyper-parameters */
    double eps;
    double t;
} brent_t;


/* MARK: Functions called by FNT */

/* \brief Provides the name of the method.
 * \param name Allocated buffer to hold the name.
 * \param size Size of the name buffer in bytes.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_name(char *name, int size) {
    if( strlcpy(name,"brents-localmin",size) >= size ) {
        return FNT_FAILURE;
    }

    return FNT_SUCCESS;
}


/* \brief Initialize intenal state for method.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_init(void **handle_ptr, int dimensions) {
    if( handle_ptr == NULL )    { return FNT_FAILURE; }
    brent_t *ptr = calloc(1, sizeof(brent_t));
    if( ptr == NULL )           { return FNT_FAILURE; }
    *handle_ptr = (void*)ptr;

    /* initialize method here */

    return FNT_SUCCESS;
}


int method_free(void **handle_ptr) {
    if( handle_ptr == NULL )    { return FNT_FAILURE; }
    if( *handle_ptr == NULL )   { return FNT_FAILURE; }
    brent_t *ptr = (brent_t*)*handle_ptr;

    /* free any memory allocated by method */

    free(ptr);  *handle_ptr = ptr = NULL;

    return FNT_SUCCESS;
}


/* \brief Display information about the method to the console.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_info() {
    printf(
"Brent's method is a minimization method, that uses a search strategy\n"
"similar to the Brent-Dekker root finding method.\n"
"\n"
"Hyper-parameters:\n"
"name\trequired\ttype\tDefault\tDescription\n"
"x_0\tREQUIRED\tdouble\tnone\tLower bound of search region.\n"
"x_1\tREQUIRED\tdouble\tnone\tUpper bound of search region.\n"
"eps\toptional\tdouble\t1e-10\tMachine epsilon.\n"
"t\toptional\tdouble\t1e-6\tMachine epsilon.\n"
"\n"
"References:\n"
"R. P. Brent, Algorithms for Minimization without Derivatives,\n\tPrentice-Hall, Englewood Cliffs, New Jersey, 1973, 195 pp.\n\tISBN 0-13-022335-2.\n"
"https://maths-people.anu.edu.au/~brent/pub/pub011.html\n"
);
    return FNT_SUCCESS;
}


/* \brief Set any hyper-parameters needed for the method.
 * \param id The name of the hyper-parameter.
 * \param value_ptr A pointer to the value being set.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_hparam_set(void *handle, char *id, void *value_ptr) {
    if( id == NULL )        { return FNT_FAILURE; }
    if( value_ptr == NULL ) { return FNT_FAILURE; }
    if( handle == NULL )    { return FNT_FAILURE; }
    brent_t *ptr = (brent_t*)handle;

    FNT_HPARAM_SET("x_0", id, double, value_ptr, ptr->a);
    FNT_HPARAM_SET("x_1", id, double, value_ptr, ptr->b);
    FNT_HPARAM_SET("eps", id, double, value_ptr, ptr->eps);
    FNT_HPARAM_SET("t", id, double, value_ptr, ptr->t);

    return FNT_SUCCESS;
}


int method_hparam_get(void *handle, char *id, void *value_ptr) {
    if( id == NULL )        { return FNT_FAILURE; }
    if( value_ptr == NULL ) { return FNT_FAILURE; }
    if( handle == NULL )    { return FNT_FAILURE; }
    brent_t *ptr = (brent_t*)handle;

    FNT_HPARAM_GET("x_0", id, double, ptr->a, value_ptr);
    FNT_HPARAM_GET("x_1", id, double, ptr->b, value_ptr);
    FNT_HPARAM_GET("eps", id, double, ptr->eps, value_ptr);
    FNT_HPARAM_GET("t", id, double, ptr->t, value_ptr);

    return FNT_SUCCESS;
}


int method_next(void *handle, fnt_vect_t *vec) {
    if( handle == NULL )    { return FNT_FAILURE; }
    if( vec == NULL )       { return FNT_FAILURE; }
    if( vec->v == NULL )    { return FNT_FAILURE; }
    brent_t *ptr = (brent_t*)handle;

    /* fill vector pointed to by vec with next input to try */
    if( ptr->state == brent_initial ) {
        double a = ptr->a;
        double b = ptr->b;

        double c = (3.0 - sqrt(5.0)) / 2.0;
        double v, w, x, e;
        v = w = x = a + c * (b - a);   e = 0;
        double d = 0.0;  /* from errata */
        FNT_VECT_ELEM(*vec, 0) = x;  /* f(x) is needed */

        ptr->c = c;
        ptr->v = v;
        ptr->w = w;
        ptr->x = x;
        ptr->e = e;
        ptr->d = d;

        DEBUG("Initializing by requesting f(x) = f(%g).\n", x);

        return FNT_SUCCESS;
    } else {
        FNT_VECT_ELEM(*vec, 0) = ptr->u; /* f(u) is needed */
        DEBUG("Requesting f(u) = f(%g).\n", ptr->u);
        return FNT_SUCCESS;
    }

    return FNT_FAILURE;
}


int method_value(void *handle, fnt_vect_t *vec, double value) {
    if( handle == NULL )    { return FNT_FAILURE; }
    if( vec == NULL )       { return FNT_FAILURE; }
    if( vec->v == NULL )    { return FNT_FAILURE; }
    brent_t *ptr = (brent_t*)handle;

    double a = ptr->a;
    double b = ptr->b;
    double c = ptr->c;
    double u = ptr->v;
    double v = ptr->v;
    double w = ptr->w;
    double x = ptr->x;
    double fu = ptr->fu;
    double fv = ptr->fv;
    double fw = ptr->fw;
    double fx = ptr->fx;
    double e = ptr->e;
    double d = ptr->d;

    if( ptr->state == brent_initial ) {
        v  = w = x = FNT_VECT_ELEM(*vec, 0);
        fv = fw = fx = value;

        DEBUG("Got initial value f(x) = f(%g) = %g.\n", x, fx);

        DEBUG("Setting state to starting.\n");
        ptr->state = brent_starting;
    }

    /* update method using value */

    if( ptr->state == brent_running ) {
        /* skip this on the first iteration,
         * as it is the end of the original loop. */
        u = FNT_VECT_ELEM(*vec, 0);
        fu = value;

        DEBUG("Updating with f(u) = f(%g) = %g.\n", u, fu);

        /* Update a, b, v, w, and x */
        if( fu <= fx ) {
                                        /* FORTRAN: 130 */
            if( u < x ) { b = x; } else { a = x; }
            /* FORTRAN: 140 */
            v = w;  fv = fw;    w = x;  fw = fx;    x = u;  fx = fu;
        } else {
            /* FORTRAN: 150 */          /* FORTRAN: 160 */
            if( u < x ) { a = u; } else { b = u; }
            /* FORTRAN: 170 */
            if( fu <= fw || w == x ) {
                v = w;  fv = fw;    w = u;  fw = fu;
            } else {
                /* FORTRAN: 180 */
                v = u;  fv = fu;
            }
        }
    } else {
        DEBUG("Setting state to running.\n");
        ptr->state = brent_running;
    }

    /* ALGOL: loop */
    /* FORTRAN: 10 */
    double m = 0.5 * (a + b);
    double tol = ptr->eps * fabs(x) + ptr->t;
    double t2 = 2.0 * tol;

    /* Check stopping criterion */
    if( fabs(x - m) > t2 - 0.5 * (b - a) ) {
        double p, q, r;
        p = q = r = 0;
        if( fabs(e) > tol ) {
            DEBUG("Fitting a parabola.\n");
            /* Fit parabola */
            r = (x - w) * (fx -fv);     q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;      q = 2 * (q - r);
                                         /* FORTRAN: 20 */
            if( q > 0 ) { p = -p; } else { q = -q; }
            /* FORTRAN: 30 */
            r = e;      e = d;
        }
        /* FORTRAN: 40 */
        if( fabs(p) < fabs(0.5 * q * r) && p > q * (a - x) && p < q * (b-x) ) {
            DEBUG("Parabolic interpolation.\n");
            /* A "parabolic interpolation" step */
            d = p / q;      u = x + d;
            /* f must not be evaluated too close to a or b */
                                                        /* FORTRAN: 50 */
            if( u - a < t2 || b - u < t2 ) { d = (x < m) ? tol : -tol; }
        } else {
            DEBUG("Golden section step.\n");
            /* A "golden section" step */
            /* FORTRAN: 60 & 70 */          /* FORTRAN: 80 */
            e = ((x < m) ? b : a) - x;      d = c * e;
        }
        /* f must not be evaluated too close to x */
        /* FORTRAN: 90 100 & 110 */
        u = x + ((fabs(d) >= tol ? d : (d > 0 ? tol : -tol)));

        DEBUG("Requesting f(u) = f(%g).\n", u);
        ptr->u = u; /* need f(u) */
        /* FORTRAN: 120 is external to this function. */
    } else {
        /* FORTRAN: 190 */
        DEBUG("Setting state to done.\n");
        ptr->state = brent_done;
        /* local minimum should be in fx */
    }

    ptr->a = a;
    ptr->b = b;
    ptr->c = c;
    ptr->u = u;
    ptr->v = v;
    ptr->w = w;
    ptr->x = x;
    ptr->fu = fu;
    ptr->fv = fv;
    ptr->fw = fw;
    ptr->fx = fx;
    ptr->e = e;
    ptr->d = d;

    return FNT_SUCCESS;
}


int method_done(void *handle) {
    if( handle == NULL )    { return FNT_FAILURE; }
    brent_t *ptr = (brent_t*)handle;

    /* test for completion
     *  Return FNT_DONE when comlete, or FNT_CONTINUE when not done.
     */
    if( ptr->state == brent_done ) {
        return FNT_DONE;
    } else {
        return FNT_CONTINUE;
    }

    return FNT_FAILURE;
}
