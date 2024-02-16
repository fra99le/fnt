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
    brent_initial, brent_initial2, brent_running, brent_done
} brent_state_t;

typedef struct brent {
    double a;
    double b;
    double c;

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


int method_info() {
    printf("\n"
           "This should give usefiul information about the method.\n"
           "\n");
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

    return FNT_FAILURE;
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

    return FNT_FAILURE;
}


int method_next(void *handle, fnt_vect_t *vec) {
    if( handle == NULL )    { return FNT_FAILURE; }
    if( vec == NULL )       { return FNT_FAILURE; }
    if( vec->v == NULL )    { return FNT_FAILURE; }
    brent_t *ptr = (brent_t*)handle;

    /* fill vector pointed to by vec with next input to try */

    return FNT_FAILURE;
}


int method_value(void *handle, fnt_vect_t *vec, double value) {
    if( handle == NULL )    { return FNT_FAILURE; }
    if( vec == NULL )       { return FNT_FAILURE; }
    if( vec->v == NULL )    { return FNT_FAILURE; }
    brent_t *ptr = (brent_t*)handle;

    if( ptr->state == brent_initial ) {
        ptr->a = vec->v[0];
        ptr->f_a = value;
        ptr->state = brent_initial2;
        return FNT_SUCCESS;
    }
    if( ptr->state == brent_initial2 ) {
        ptr->b = vec->v[0];
        ptr->f_b = value;

        /* Note: some of this could possiblt be moved to the init function,
         * as only one  objective fucntion value is needed. */
        c = (3.0 - sqrt(5.0)) / 2.0;
        v  = w = x = a + c * (b - a);   e = 0;
        fv = fw = fx = f(x); /* THIS WILL NEED TO RETURN AND COME BACK! */

        ptr->state = brent_running;
    }

    /* update method using value */

    double a = ptr->a;
    double b = ptr->b;
    double c = ptr->c;


    /* ALGOL: loop */
    double m = 0.5 * (a + b);
    double tol = ptr->eps * fabs(x) + ptr->t;       t2 = 2.0 * tol;

    /* Check stopping criterion */
    if( fabs(x - m) > t2 - 0.5 * (b - a) ) {
        p = q = r = 0;
        if( fabs(e) > tol ) {
            /* Fit parabola */
            r = (x - w) * (fx -fv);     q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;      q = 2 * (q - r);
            if( q > 0 ) { p = -p; } else { q = -q; }
            r = e;      e = d;
        }
        if( fabs(p) < abs(0.5 * q * r) && p < q * (a - x) && p < q * (b-x) ) {
            /* A "parabolic interpolation" step */
            d = p / q;      u = x + d;
            /* f must not be evaluated too close to a or b */
            if( u - a < t2 || b - u < t2 ) { d = (x < m) ? tol : -tol; }
        } else {
            /* A "golden section" step */
            e = ((x < m) ? b : a) - x;      d = c * e;
        }
        /* f must not be evaluated too close to x */
        u = x + ((fabs(d) >= tol ? d : (d > 0 ? tol : -tol)));
        fu = f(u);  /* THIS WILL NEED TO RETURN AND COME BACK! */
        /* Update a, b, v, w, and x */
        if( fu <= fx ) {
            if( u < x ) { b = x; } else { a = x; }
            v = w;  fv = fw;    w = x;  fw = fx;    x = u;  fx = fu;
        } else {
            if( u < x ) { a = u; } else { b = u; }
            if( fu <= fw || w = x ) {
                v = w;  fv = fw;    w = u;  fw = fu;
            } else {
                v = u;  fv = fu;
            }
        }
    } else {
        ptr->state = brent_done;
        /* local minimum should be in fx */
    }

    return FNT_FAILURE;
}


int method_value_gradient(void *handle, fnt_vect_t *vec, double value, fnt_vect_t gradient) {
    if( handle == NULL )    { return FNT_FAILURE; }
    if( vec == NULL )       { return FNT_FAILURE; }
    if( vec->v == NULL )    { return FNT_FAILURE; }
    brent_t *ptr = (brent_t*)handle;

    /* update method using value and derivative/gradient */

    return FNT_FAILURE;
}


int method_done(void *handle) {
    if( handle == NULL )    { return FNT_FAILURE; }
    brent_t *ptr = (brent_t*)handle;

    /* test for completion
     *  Return FNT_DONE when comlete, or FNT_CONTINUE when not done.
     */

    return FNT_FAILURE;
}
