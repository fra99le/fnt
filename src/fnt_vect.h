/*
 * fnt_vect.h
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fnt_util.h"

#define FNT_VEC_SUCCESS 0
#define FNT_VEC_FAILURE 1


/* MARK: Type and extern declarations */

extern int fnt_verbose_level;

typedef struct fnt_vect {
    double *v;
    size_t n;
} fnt_vect_t;


/* MARK: Vector memory operations */

static int fnt_vect_calloc(fnt_vect_t *vec, int length) {
    if( vec == NULL )   { return FNT_VEC_FAILURE; }

    if( (vec->v = calloc(length, sizeof(double))) == NULL ) {
        if( fnt_verbose_level >= 1 ) {
            perror("calloc");
        }
        return FNT_VEC_FAILURE;
    }
    vec->n = length;

    return FNT_VEC_SUCCESS;
}


static int fnt_vect_free(fnt_vect_t *vec) {
    if( vec == NULL )   { return FNT_VEC_FAILURE; }

    if( vec->v != NULL ) {
        free(vec->v); vec->v = NULL;
        vec->n = 0;
    }

    return FNT_VEC_SUCCESS;
}


static int fnt_vect_reset(fnt_vect_t *vec) {
    if( vec == NULL )       { return FNT_VEC_FAILURE; }
    if( vec->v == NULL )    { return FNT_VEC_FAILURE; }

    memset(vec->v, '\0', vec->n * sizeof(double));

    return FNT_VEC_SUCCESS;
}


static int fnt_vect_copy(fnt_vect_t *dst, fnt_vect_t *src) {
    if( dst == NULL )       { return FNT_VEC_FAILURE; }
    if( src == NULL )       { return FNT_VEC_FAILURE; }
    if( dst->v == NULL )    { return FNT_VEC_FAILURE; }
    if( src->v == NULL )    { return FNT_VEC_FAILURE; }

    if( dst->n != src->n ) {
        if( fnt_verbose_level >= 1 ) {
            fprintf(stderr, "%s: source length (%zu) differant than destination length (%zu).", __FUNCTION__, src->n, dst->n);
        }

        return FNT_VEC_FAILURE;
    }

    return FNT_VEC_SUCCESS;
}


/* MARK: Vector I/O */

static int fnt_vect_print(fnt_vect_t *vec, char *label, char *fmt) {
    if( vec == NULL )   { return FNT_VEC_FAILURE; }
    /* label and fmt can each be  NULL */

    if( label != NULL ) { printf("%s", label); }

    printf("[");
    for(int i=0; i<vec->n; ++i) {
        if( fmt != NULL )       { printf(fmt, vec->v[i]); }
        else                    { printf("%g", vec->v[i]); }

        if( i < vec->n - 1 )    { printf(", "); }
    }
    printf("]");

    return FNT_VEC_SUCCESS;
}


static int fnt_vect_println(fnt_vect_t *vec, char *label, char *fmt) {
    int ret = fnt_vect_print(vec, label, fmt);
    if( ret == FNT_VEC_SUCCESS ) { printf("\n"); }

    return ret;
}


/* MARK: Basic vector operations */

static int fnt_vect_add(fnt_vect_t *sum, fnt_vect_t *a, fnt_vect_t *b) {
    if( sum == NULL )       { return FNT_VEC_FAILURE; }
    if( a == NULL )         { return FNT_VEC_FAILURE; }
    if( b == NULL )         { return FNT_VEC_FAILURE; }
    if( sum->v == NULL )    { return FNT_VEC_FAILURE; }
    if( a->v == NULL )      { return FNT_VEC_FAILURE; }
    if( b->v == NULL )      { return FNT_VEC_FAILURE; }
    if( sum->n != a->n )    { return FNT_VEC_FAILURE; }
    if( b->n != a->n )      { return FNT_VEC_FAILURE; }

    for(int i=0; i<a->n; ++i) {
        sum->v[i] = a->v[i] + b->v[i];
    }

    return FNT_VEC_SUCCESS;
}


static int fnt_vect_sub(fnt_vect_t *diff, fnt_vect_t *a, fnt_vect_t *b) {
    if( diff == NULL )      { return FNT_VEC_FAILURE; }
    if( a == NULL )         { return FNT_VEC_FAILURE; }
    if( b == NULL )         { return FNT_VEC_FAILURE; }
    if( diff->v == NULL )   { return FNT_VEC_FAILURE; }
    if( a->v == NULL )      { return FNT_VEC_FAILURE; }
    if( b->v == NULL )      { return FNT_VEC_FAILURE; }
    if( diff->n != a->n )   { return FNT_VEC_FAILURE; }
    if( b->n != a->n )      { return FNT_VEC_FAILURE; }

    for(int i=0; i<a->n; ++i) {
        diff->v[i] = a->v[i] - b->v[i];
    }

    return FNT_VEC_SUCCESS;
}


static int fnt_vect_scale(fnt_vect_t *result, double scaling, fnt_vect_t *vec) {
    if( vec == NULL )           { return FNT_VEC_FAILURE; }
    if( result == NULL )        { return FNT_VEC_FAILURE; }
    if( vec->v == NULL )        { return FNT_VEC_FAILURE; }
    if( result->v == NULL )     { return FNT_VEC_FAILURE; }
    if( result->n != vec->n )   { return FNT_VEC_FAILURE; }

    for(int i=0; i<vec->n; ++i) {
        result->v[i] = scaling * vec->v[i];
    }

    return FNT_VEC_SUCCESS;
}


/* MARK: Advanced vector operations */

static int fnt_vect_l2norm(fnt_vect_t *vec, double *result) {
    if( vec == NULL )           { return FNT_VEC_FAILURE; }
    if( result == NULL )        { return FNT_VEC_FAILURE; }
    if( vec->v == NULL )        { return FNT_VEC_FAILURE; }

    double sum = 0.0;
    for(int i=0; i<vec->n; ++i) {
        sum += pow(vec->v[i], 2.0);
    }
    *result = sqrt(sum);

    return FNT_VEC_SUCCESS;
}


static int fnt_vect_dist(fnt_vect_t *a, fnt_vect_t *b, double *result) {
    if( result == NULL )    { return FNT_VEC_FAILURE; }
    if( a == NULL )         { return FNT_VEC_FAILURE; }
    if( b == NULL )         { return FNT_VEC_FAILURE; }
    if( a->v == NULL )      { return FNT_VEC_FAILURE; }
    if( b->v == NULL )      { return FNT_VEC_FAILURE; }
    if( b->n != a->n )      { return FNT_VEC_FAILURE; }

    fnt_vect_t diff;
    fnt_vect_calloc(&diff, a->n);
    fnt_vect_sub(&diff, a, b);
    fnt_vect_l2norm(&diff, result);
    fnt_vect_free(&diff);

    return FNT_VEC_SUCCESS;
}
