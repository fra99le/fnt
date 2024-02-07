/*
 * fnt_util.h
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */
#ifndef FNT_UTIL_H
#define FNT_UTIL_H

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

/* MARK: Hyper-parameter accessing macros */

#define FNT_HPARAM_SET(name, id, type, src_ptr, dst) \
    if( strncmp((name), (id), strlen(name)) == 0 ) { \
        (dst) = *(type*)(src_ptr); \
    }

#define FNT_HPARAM_GET(name, id, type, src, dst_ptr) \
    if( strncmp((name), (id), strlen(name)) == 0 ) { \
        *(type*)(dst_ptr) = (src); \
    }

extern int fnt_verbose_level;

#endif /* FNT_UTIL_H */
