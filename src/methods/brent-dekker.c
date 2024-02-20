/*
 * brent_dekker.c
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../fnt.h"


/* MARK: Method type definitions */

typedef enum brent_dekker_state {
    brent_dekker_initial,
    brent_dekker_initial2,
    brent_dekker_starting,
    brent_dekker_running,
    brent_dekker_done
} brent_dekker_state_t;

typedef struct brent_dekker {

    /* execution state */
    brent_dekker_state_t state;

    /* hyper-parameters */
    double macheps; /* machine epsilon */
    double t;       /* a positive tolerance */

    /* method state variables */
    double a;
    double b;
    double c;
    double f_a;
    double f_b;
    double f_c;
    double d;
    double e;

} brent_dekker_t;


/* MARK: Functions called by FNT */

/* \brief Provides the name of the method.
 * \param name Allocated buffer to hold the name.
 * \param size Size of the name buffer in bytes.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_name(char *name, int size) {
    if( snprintf(name,  size,  "brent-dekker") >= size ) {
        return FNT_FAILURE;
    }

    return FNT_SUCCESS;
}


/* \brief Initialize intenal state for method.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_init(void **handle_ptr, int dimensions) {
    if( handle_ptr == NULL )    { return FNT_FAILURE; }
    brent_dekker_t *ptr = calloc(1, sizeof(brent_dekker_t));
    if( ptr == NULL )           { return FNT_FAILURE; }
    *handle_ptr = (void*)ptr;

    if( dimensions > 1 ) {
        ERROR("ERROR: Brent-Dekker is a single variate method, %d dmensions requested.\n", dimensions);
        return FNT_FAILURE;
    }

    /* initialize method here */
    ptr->state = brent_dekker_initial;
    ptr->macheps = 1e-10;
    ptr->t = 1e-6;

    return FNT_SUCCESS;
}


int method_free(void **handle_ptr) {
    if( handle_ptr == NULL )    { return FNT_FAILURE; }
    if( *handle_ptr == NULL )   { return FNT_FAILURE; }
    brent_dekker_t *ptr = (brent_dekker_t*)*handle_ptr;

    /* free any memory allocated by method */

    free(ptr);  *handle_ptr = ptr = NULL;

    return FNT_SUCCESS;
}


/* \brief Display information about the method to the console.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int method_info() {
    printf(
"The Brent-Dekker method is a root finding method, similar to bisection,\n"
"that uses multiple strategies that, in general, reduce the search space\n"
"faster than the bisection method.\n"
"\n"
"Hyper-parameters:\n"
"name\trequired\ttype\tDefault\tDescription\n"
"x_0\tREQUIRED\tdouble\tnone\tLower bound of search region.\n"
"x_1\tREQUIRED\tdouble\tnone\tUpper bound of search region.\n"
"macheps\toptional\tdouble\t1e-10\tMachine epsilon.\n"
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
    if( handle == NULL )    { return FNT_FAILURE; }
    brent_dekker_t *ptr = (brent_dekker_t*)handle;
    if( id == NULL )        { return FNT_FAILURE; }
    if( value_ptr == NULL ) { return FNT_FAILURE; }

    FNT_HPARAM_SET("macheps", id, double, value_ptr, ptr->macheps);
    FNT_HPARAM_SET("t", id, double, value_ptr, ptr->t);
    FNT_HPARAM_SET("x_0", id, double, value_ptr, ptr->a);
    FNT_HPARAM_SET("x_1", id, double, value_ptr, ptr->b);

    return FNT_SUCCESS;
}


int method_hparam_get(void *handle, char *id, void *value_ptr) {
    if( handle == NULL )    { return FNT_FAILURE; }
    brent_dekker_t *ptr = (brent_dekker_t*)handle;
    if( id == NULL )        { return FNT_FAILURE; }
    if( value_ptr == NULL ) { return FNT_FAILURE; }

    FNT_HPARAM_GET("macheps", id, double, ptr->macheps, value_ptr);
    FNT_HPARAM_GET("t", id, double, ptr->t, value_ptr);
    FNT_HPARAM_GET("x_0", id, double, ptr->a, value_ptr);
    FNT_HPARAM_GET("x_1", id, double, ptr->b, value_ptr);

    return FNT_SUCCESS;
}


int method_next(void *handle, fnt_vect_t *vec) {
    if( handle == NULL )    { return FNT_FAILURE; }
    brent_dekker_t *ptr = (brent_dekker_t*)handle;
    if( vec == NULL )       { return FNT_FAILURE; }
    if( vec->v == NULL )    { return FNT_FAILURE; }

    /* fill vector pointed to by vec with next input to try */
    if( ptr->state == brent_dekker_initial ) {
        FNT_VECT_ELEM(*vec, 0) = ptr->a;
        return FNT_SUCCESS;
    }

    /* after initialization, only f(b) is required per iteration. */
    FNT_VECT_ELEM(*vec, 0) = ptr->b;

    return FNT_SUCCESS;
}


int method_value(void *handle, fnt_vect_t *vec, double value) {
    if( handle == NULL )    { return FNT_FAILURE; }
    brent_dekker_t *ptr = (brent_dekker_t*)handle;
    if( vec == NULL )       { return FNT_FAILURE; }
    if( vec->v == NULL )    { return FNT_FAILURE; }

    /* update method using value */
    if( ptr->state == brent_dekker_initial ) {
        ptr->a = FNT_VECT_ELEM(*vec, 0);
        ptr->f_a = value;

        ptr->state = brent_dekker_initial2;

        return FNT_SUCCESS;
    }
    if( ptr->state == brent_dekker_initial2 ) {
        ptr->b = FNT_VECT_ELEM(*vec, 0);
        ptr->f_b = value;

        if( ptr->f_a * ptr->f_b > 0.0 ) {
            ERROR("Objective function must have opposite sign at each end of the search region (f(%g)=%g; f(%g)=%g)\n", ptr->a, ptr->f_a, ptr->b, ptr->f_b);
            ptr->state = brent_dekker_done;
            return FNT_FAILURE;
        } else {
            INFO("f(a) and f(b) have different signs, as required.\n");
        }

        ptr->state = brent_dekker_starting;
    }

    /* check state */
    if( ptr->state != brent_dekker_starting
        && ptr->state != brent_dekker_running ) {
        ERROR("Should be in starting or running state, but is not.\n");
    }

    /* perform Brent-Dekker update to a, b, & c. */

    /* copy common values into local variables */
    double a = ptr->a;
    double b = FNT_VECT_ELEM(*vec, 0);   /* update b */
    double c = ptr->c;
    double f_a = ptr->f_a;
    double f_b = value;     /* update f(b) */
    double f_c = ptr->f_c;
    double d = ptr->d;
    double e = ptr->e;

    /* FORTRAN: 130-ish */
    if( (f_b > 0.0 && f_c > 0.0)
        || (f_b <= 0.0 && f_c <= 0.0)
        || ptr->state == brent_dekker_starting ) {
        /* ALGOL: int */
        /* FORTRAN: 10 */
        c = a;      f_c = f_a;      d = e = b - a;
        if( ptr->state == brent_dekker_starting )
            ptr->state = brent_dekker_running;
    }

    /* ALGOL: ext */
    /* FORTRAN: 20 */
    if( fabs(f_c) < fabs(f_b) ) {
        a = b;      b = c;      c = a;
        f_a = f_b;  f_b = f_c;  f_c = f_a;
    }

    /* FORTRAN: 30 */
    double tol = 2.0 * ptr->macheps * fabs(b) + ptr->t;
    double m = 0.5 * (c - b);

    if( fabs(m) > tol && f_b != 0.0 ) {
        /* see if bisection is forced */
        if( fabs(e) < tol || fabs(f_a) <= fabs(f_b) ) {
            d = e = m;
        } else {
            /* FORTRAN: 40 */
            double p, q;
            double s = f_b / f_a;
            if( a == c ) {
                /* Linear interpolation */
                p = 2.0 * m * s;      q = 1.0 - s;
            } else {
                /* FORTRAN: 50 */
                double r;
                /* Inverse quadratic interpolation */
                q = f_a / f_c;      r = f_b / f_c;
                p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }

            /* FORTRAN: 60 */           /* FORTRAN: 70 */
            if( p > 0 ) { q = -q; } else { p = -p; }

            /* FORTRAN: 80 */
            s = e;      e = d;

            if( (2.0 * p < 3.0 * m * q - fabs(tol * q))
                && (p < fabs(0.5 * s * q)) ) {
                d = p / q;
            } else {
                /* FORTRAN: 90 */
                d = e = m;
            }
        }

        /* FORTRAN: 100 */
        a = b;      f_a = f_b;
                                       /* FORTRAN: 110 & 120 */
        b = b + ((fabs(d) > tol) ? d : ((m > 0) ? tol : -tol));

        /* need updated f(b) */
        /* FORTRAN: 130 is external to this function */
    } else {
        /* FORTRAN: 140 */
        /* method has converged to a solution */
        ptr->state = brent_dekker_done;
        /* b contains the result. */
    }

    /* copy local variables back to presistent struct */
    ptr->a = a;
    ptr->b = b;
    ptr->c = c;
    ptr->f_a = f_a;
    ptr->f_b = f_b;
    ptr->f_c = f_c;
    ptr->d = d;
    ptr->e = e;

    return FNT_SUCCESS;
}


int method_done(void *handle) {
    if( handle == NULL )    { return FNT_FAILURE; }
    brent_dekker_t *ptr = (brent_dekker_t*)handle;

    /* test for completion
     *  Return FNT_DONE when comlete, or FNT_CONTINUE when not done.
     */
    if( ptr->state == brent_dekker_done ) {
        return FNT_DONE;
    } else {
        return FNT_CONTINUE;
    }

    return FNT_FAILURE;
}
