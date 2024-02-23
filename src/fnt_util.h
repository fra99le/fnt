/*
 * fnt_util.h
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */
#ifndef FNT_UTIL_H
#define FNT_UTIL_H

/* MARK: Constants */

#define FNT_SUCCESS     0
#define FNT_FAILURE     1
#define FNT_CONTINUE    2
#define FNT_DONE        3

#define FNT_NONE    0
#define FNT_ERROR   1
#define FNT_WARN    2
#define FNT_INFO    3
#define FNT_DEBUG   4


/* MARK: Random number generation macros */

#ifndef FNT_RAND
#define FNT_RAND        rand
#endif /* FNT_RAND */
#ifndef FNT_RAND_MAX
#define FNT_RAND_MAX    RAND_MAX
#endif /* FNT_RAND_MAX */

/* MARK: Console output macros */

#define ERROR(...)  if( fnt_verbose_level >= FNT_ERROR ) { fprintf(stderr, __VA_ARGS__); }
#define WARN(...)   if( fnt_verbose_level >= FNT_WARN ) { fprintf(stderr, __VA_ARGS__); }
#define INFO(...)   if( fnt_verbose_level >= FNT_INFO ) { printf(__VA_ARGS__); }
#define DEBUG(...)  if( fnt_verbose_level >= FNT_DEBUG ) { printf(__VA_ARGS__); }

/* MARK: Hyper-parameter accessing macros */

#define FNT_HPARAM_SET(name, id, type, src_ptr, dst) \
    if( strncmp((name), (id), strlen(name)) == 0 ) { \
        (dst) = *(type*)(src_ptr); \
    }

#define FNT_HPARAM_GET(name, id, type, src, dst_ptr) \
    if( strncmp((name), (id), strlen(name)) == 0 ) { \
        *(type*)(dst_ptr) = (src); \
    }


/* MARK: Result accessing macros */

#define FNT_RESULT_GET(name, id, type, src, dst_ptr) \
    if( strncmp((name), (id), strlen(name)) == 0 ) { \
        *(type*)(dst_ptr) = (src); \
    }

#define FNT_RESULT_GET_VECT(name, id, src, dst_ptr) \
    if( strncmp((name), (id), strlen(name)) == 0 ) { \
        fnt_vect_copy(dst_ptr, &src); \
    }


/* MARK: Externed Global Variables */

extern int fnt_verbose_level;

#endif /* FNT_UTIL_H */
