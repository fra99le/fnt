/*
 * fnt.c
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */
#include <dirent.h>
#include <dlfcn.h>
#include <limits.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "fnt.h"

/* MARK: Internal constants */

#define FNT_MAX_NAME_LENGTH 64

/* MARK: Internal variables */

static int fnt_verbose_level = 1;   /* default to showing errors */

/* MARK: Internal structures */

typedef struct fnt_method_list_entry {
    char name[FNT_MAX_NAME_LENGTH];
    char path[PATH_MAX];
} fnt_method_list_entry_t;

typedef struct fnt_method_list {
    int count;
    int capacity;
    fnt_method_list_entry_t *entries;
} fnt_method_list_t;

typedef struct fnt_method {
    char name[FNT_MAX_NAME_LENGTH];
    int (*init)();
    int (*info)();
    int (*hparam_set)(char *id, void *value_ptr);
    int (*hparam_get)(char *id, void *value_ptr);
    int (*seed)(double *vec);
    int (*next)(double *vec);
    int (*value)(double *vec, double value);
    int (*done)();
    int (*result)(void*);
    int (*free)();
} fnt_method_t;

typedef struct vector_queue_node {
    double *vector;
    int len;
    struct vector_queue_node *next;
} vector_queue_node_t;

typedef struct fnt_context {
    /* list of available methods */
    fnt_method_list_t methods_list;

    /* handle for currently loaded method */
    void *dl_handle;

    /* number of input dimensions */
    int dim;

    /* loaded method, NULL otherwise */
    fnt_method_t method;

    /* list of upcoming inputs that are needed */
    vector_queue_node_t *inputs_head;
} context_t;

/* MARK: Internal functions */

int fnt_method_list_init(fnt_method_list_t *list, int initial_cap) {

    if( (list->entries
            = calloc(initial_cap, sizeof(fnt_method_list_entry_t)))
            == NULL ) {
        if( fnt_verbose_level >= 1 ) {
            perror("calloc");
        }
        return FNT_FAILURE;
    }
    list->capacity = initial_cap;
    list->count = 0;

    return FNT_SUCCESS;
}


int fnt_method_list_free(fnt_method_list_t *list) {
    if( list == NULL )       { return FNT_FAILURE; }

    if( list->entries != NULL ) {
        free(list->entries); list->entries = NULL;
    }
    list->count = list->capacity = 0;

    return FNT_SUCCESS;
}


int fnt_method_list_add(context_t *ctx, fnt_method_list_entry_t *entry) {
    if( ctx == NULL )   { return FNT_FAILURE; }
    if( entry == NULL ) { return FNT_FAILURE; }

    /* allocate list, if needed */
    if( ctx->methods_list.entries == NULL ) {
        fnt_method_list_init(&ctx->methods_list, 10);
    }

    if( ctx->methods_list.count == ctx->methods_list.capacity ) {
        size_t new_size = (ctx->methods_list.capacity*2 + 1)
                                * sizeof(fnt_method_list_entry_t);
        void *ptr = realloc(ctx->methods_list.entries, new_size);
        if( ptr == NULL )   {
            if( fnt_verbose_level >= 1 ) {
                perror("realloc");
            }
            return FNT_FAILURE;
        }
        ctx->methods_list.entries = ptr;
    }

    int pos = ctx->methods_list.count;
    memcpy(&ctx->methods_list.entries[pos], entry, sizeof(fnt_method_list_entry_t));
    ++ctx->methods_list.count;

    return FNT_SUCCESS;
}


int fnt_register_method(context_t *ctx, char *filename) {

    /* open object file */
    void *dl_handle = NULL;
    if( !(dl_handle = dlopen(filename,RTLD_NOW)) ) {
        if( fnt_verbose_level >= 1 ) {
            fprintf(stderr, "ERROR: dlopen: %s\n", dlerror());
        }
        return FNT_FAILURE;
    }

    /* extract name function */
    int (*fnt_method_name)(char *, int) = dlsym(dl_handle, "method_name");
    if( fnt_method_name == NULL ) {
        if( fnt_verbose_level >= 1 ) {
            fprintf(stderr, "ERROR: No fnt_method_name in '%s'.\n", filename);
        }
        dlclose(dl_handle);
        return FNT_FAILURE;
    }
    char name[FNT_MAX_NAME_LENGTH];
    fnt_method_name(name, FNT_MAX_NAME_LENGTH);

    /* set up list entry */
    fnt_method_list_entry_t *entry = calloc(1,sizeof(fnt_method_list_entry_t));
    strlcpy(entry->name, name, sizeof(entry->name));
    strlcpy(entry->path, filename, sizeof(entry->path));

    /* add method entry to list of available methods */
    fnt_method_list_add(ctx, entry);

    if( fnt_verbose_level >= 2 ) {
        printf("\tloaded method '%s' from '%s'.\n", entry->name, filename);
    }

    dlclose(dl_handle);

    return FNT_SUCCESS;
}


int fnt_register_methods(context_t *ctx, char *path) {
    /* check input */
    if( path == NULL ) {
        if( fnt_verbose_level >= 1 ) {
            fprintf(stderr, "ERROR: Methods path is NULL.\n");
        }
        return FNT_FAILURE;
    }

    /* open directory */
    DIR *dirp = NULL;
    if( (dirp = opendir(path)) == NULL ) {
        if( fnt_verbose_level >= 1 ) {
            perror("opendir");
            return FNT_FAILURE;
        }
    }

    /* read contents of the directory */
    if( fnt_verbose_level >= 2 ) {
        printf("Loading methods from '%s'.\n", path);
    }

    struct dirent *dp;
    char filename[PATH_MAX];
    while( (dp = readdir(dirp)) != NULL ) {
        /* check extension */
        int namelen = strlen(dp->d_name);
        if( namelen <= 3
            || dp->d_name[0] == '.'
            || strncasecmp(dp->d_name+namelen-3, ".so", 3) ) {
            if( fnt_verbose_level >= 3 ) {
                printf("DEBUG: Skipping '%s'.\n", dp->d_name);
            }
            continue;
        }

        snprintf(filename,sizeof(filename), "%s/%s", path, dp->d_name);

        if( fnt_register_method(ctx, filename) == FNT_FAILURE ) {
            if( fnt_verbose_level >= 1 ) {
                fprintf(stderr, "ERROR: Failed to load '%s' from '%s'.\n", dp->d_name, path);
            }
        }
    }

    closedir(dirp); dirp = NULL;

    return FNT_SUCCESS;
}


int fnt_method_load(context_t *ctx, char *filename) {
    if( ctx == NULL )       { return FNT_FAILURE; }
    if( filename == NULL )  { return FNT_FAILURE; }

    /* open method file */
    void *dl_handle = NULL;
    if( !(dl_handle = dlopen(filename, RTLD_NOW)) ) {
        if( fnt_verbose_level >= 1 ) {
            fprintf(stderr, "ERROR: dlopen: %s\n", dlerror());
        }
        return FNT_FAILURE;
    }

    printf("fnt_verbose_level: %i\n", fnt_verbose_level);
    if( fnt_verbose_level >= 2 ) {
        printf("Loading method from '%s'.\n", filename);
    }

    /* assign function pointers */
    ctx->dl_handle = dl_handle;
    int (*method_name)(char *, int) = dlsym(dl_handle, "method_name");
    method_name(ctx->method.name, sizeof(ctx->method.name));
    ctx->method.init = dlsym(dl_handle, "method_init");
    ctx->method.info = dlsym(dl_handle, "method_info");
    ctx->method.hparam_get = dlsym(dl_handle, "method_hparam_get");
    ctx->method.hparam_set = dlsym(dl_handle, "method_hparam_set");
    ctx->method.next = dlsym(dl_handle, "method_next");
    ctx->method.value = dlsym(dl_handle, "method_value");
    ctx->method.done = dlsym(dl_handle, "method_done");
    ctx->method.result = dlsym(dl_handle, "method_result");
    ctx->method.free = dlsym(dl_handle, "method_free");

    if( ctx->method.next == NULL
        || ctx->method.value == NULL
        || ctx->method.done == NULL ) {
        fprintf(stderr, "ERROR: '%s' does not have all required methods.\n", filename);
        dlclose(dl_handle); ctx->dl_handle = dl_handle = NULL;
        return FNT_FAILURE;
    }

    /* TODO: Need to set default versions of optional methods. */

    return ctx->method.init();
}

/* MARK: User callable functions */

int fnt_init(void **context, char *path) {
    context_t *ctx = calloc(1,sizeof(context_t));
    *context = ctx;

    /* load list of available methods */
    return fnt_register_methods(ctx, path);
}


int fnt_verbose(int verbosity) {
    fnt_verbose_level = verbosity;
    if( fnt_verbose_level >= 2 ) {
        printf("Verbosity set to %i.\n", fnt_verbose_level);
    }
    return FNT_SUCCESS;
}


int fnt_set_method(void *context, char *name, int dimensions) {
    context_t *ctx = context;
    if( ctx == NULL )   { return FNT_FAILURE; }
    if( name == NULL )  { return FNT_FAILURE; }

    ctx->dim = dimensions;

    /* dynamically load module */
    for(int i=0; i<ctx->methods_list.count; ++i) {
        printf("checking %s\n", ctx->methods_list.entries[i].name);
        if( strncmp(ctx->methods_list.entries[i].name, name,
                sizeof(ctx->methods_list.entries[i].name)) == 0 ) {
            return fnt_method_load(ctx, ctx->methods_list.entries[i].path);
        }
    }

    return FNT_FAILURE;
}


int fnt_free(void **context) {
    if( context == NULL )   { return FNT_FAILURE; }
    context_t *ctx = *(context_t**)context;
    if( ctx == NULL )       { return FNT_FAILURE; }

    int ret = FNT_SUCCESS;
    if( ctx->method.free != NULL )  {
        ret = ctx->method.free();

        if( ret == FNT_SUCCESS && fnt_verbose_level >= 2 ) {
            printf("Freed intenally allocated values for method '%s'.\n", ctx->method.name);
        } else if( ret == FNT_FAILURE && fnt_verbose_level >= 1 ) {
            fprintf(stderr, "ERROR: Call to free for method '%s' failed..\n", ctx->method.name);
        }
    }

    fnt_method_list_free(&ctx->methods_list);
    dlclose(ctx->dl_handle);    ctx->dl_handle = NULL;

    /* TODO: empty eny queued values */

    if( ret == FNT_SUCCESS ) {
        free(*context); *context = ctx = NULL;
    }

    return ret;
}


int fnt_info(void *context) {
    context_t *ctx = (context_t*)context;
    if( ctx == NULL )               { return FNT_FAILURE; }
    if( ctx->method.info == NULL )  { return FNT_FAILURE; }

    if( ctx->method.name[0] == '\0' ) {
        if( fnt_verbose_level >= 1 ) {
            fprintf(stderr, "ERROR: Called %s before setting method.\n", __FUNCTION__);
        }
        return FNT_FAILURE;
    }
    if( ctx->method.info == NULL ) {
        if( fnt_verbose_level >= 1 ) {
            fprintf(stderr, "ERROR: Method '%s' does not provide additional info.\n", ctx->method.name);
        }
        return FNT_FAILURE;
    }

    return ctx->method.info();
}


int fnt_hparam_set(void *context, char *id, void *value_ptr) {
    context_t *ctx = (context_t*)context;
    if( ctx == NULL )                       { return FNT_FAILURE; }
    if( ctx->method.hparam_set == NULL )    { return FNT_FAILURE; }
    if( id == NULL )                        { return FNT_FAILURE; }
    if( value_ptr == NULL )                 { return FNT_FAILURE; }

    int ret = ctx->method.hparam_set(id, value_ptr);

    if( ret == FNT_SUCCESS && fnt_verbose_level >= 2 ) {
        printf("Set hyper-parameter '%s'.\n", id);
    } else if( ret == FNT_FAILURE && fnt_verbose_level >= 1 ) {
        fprintf(stderr, "ERROR: Failed to set hyper-parameter '%s'.\n", id);
    }

    return ret;
}


int fnt_hparam_get(void *context, char *id, void *value_ptr) {
    context_t *ctx = (context_t*)context;
    if( ctx == NULL )                       { return FNT_FAILURE; }
    if( ctx->method.hparam_get == NULL )    { return FNT_FAILURE; }
    if( id == NULL )                        { return FNT_FAILURE; }
    if( value_ptr == NULL )                 { return FNT_FAILURE; }

    int ret = ctx->method.hparam_get(id, value_ptr);

    if( ret == FNT_SUCCESS && fnt_verbose_level >= 2 ) {
        printf("Got hyper-parameter '%s'.\n", id);
    } else if( ret == FNT_FAILURE && fnt_verbose_level >= 1 ) {
        fprintf(stderr, "ERROR: Failed to get hyper-parameter '%s'.\n", id);
    }

    return ret;
}


int fnt_seed(void *context, double *vec) {
    context_t *ctx = (context_t*)context;
    if( ctx == NULL )               { return FNT_FAILURE; }
    if( ctx->method.seed == NULL )  { return FNT_FAILURE; }
    if( vec == NULL )               { return FNT_FAILURE; }

    int ret = ctx->method.seed(vec);

    if( ret == FNT_SUCCESS && fnt_verbose_level >= 2 ) {
        printf("Seeded input vector.\n");
    } else if( ret == FNT_FAILURE && fnt_verbose_level >= 1 ) {
        fprintf(stderr, "ERROR: Failed to seed input vector.\n");
    }

    return ret;
}


int fnt_next(void *context, double *vec) {
    context_t *ctx = (context_t*)context;
    if( ctx == NULL )               { return FNT_FAILURE; }
    if( ctx->method.next == NULL )  { return FNT_FAILURE; }
    if( vec == NULL )               { return FNT_FAILURE; }

    int ret =  ctx->method.next(vec);

    if( ret == FNT_SUCCESS && fnt_verbose_level >= 3 ) {
        printf("DEBUG: Retrieved next input vector.\n");
    } else if( ret == FNT_FAILURE && fnt_verbose_level >= 1 ) {
        fprintf(stderr, "ERROR: Failed to retreive next input vector.\n");
    }

    return ret;
}


int fnt_set_value(void *context, double *vec, double value) {
    context_t *ctx = (context_t*)context;
    if( ctx == NULL )               { return FNT_FAILURE; }
    if( ctx->method.value == NULL ) { return FNT_FAILURE; }
    if( vec == NULL )               { return FNT_FAILURE; }

    int ret = ctx->method.value(vec, value);;

    if( ret == FNT_SUCCESS && fnt_verbose_level >= 3 ) {
        printf("DEBUG: Set objective function value for input vector.\n");
    } else if( ret == FNT_FAILURE && fnt_verbose_level >= 1 ) {
        fprintf(stderr, "ERROR: Failed to set objective value for input vector.\n");
    }

    return ret;
}


int fnt_done(void *context) {
    context_t *ctx = (context_t*)context;
    if( ctx == NULL )               { return FNT_FAILURE; }
    if( ctx->method.done == NULL )  { return FNT_FAILURE; }

    int ret = ctx->method.done();

    if( ret == FNT_DONE && fnt_verbose_level >= 3 ) {
        printf("DEBUG: Method '%s' has finished.\n", ctx->method.name);
    } else if( ret == FNT_FAILURE && fnt_verbose_level >= 1 ) {
        fprintf(stderr, "ERROR: Method completion check failed.\n");
    }

    return ret;
}


int fnt_best(void *context, double *vec) {
    context_t *ctx = (context_t*)context;
    if( ctx == NULL )               { return FNT_FAILURE; }
    if( ctx->method.next == NULL )  { return FNT_FAILURE; }
    if( vec == NULL )               { return FNT_FAILURE; }

    /* TODO: add tracking of the best value. */
    int ret = FNT_FAILURE;

    if( ret == FNT_SUCCESS && fnt_verbose_level >= 3 ) {
        printf("DEBUG: Retrieved best input vector.\n");
    } else if( ret == FNT_FAILURE && fnt_verbose_level >= 1 ) {
        fprintf(stderr, "ERROR: Failed to retreive best input vector.\n");
    }

    return ret;
}


int fnt_result(void *context, void *extra) {
    context_t *ctx = (context_t*)context;
    if( ctx == NULL )                   { return FNT_FAILURE; }
    if( ctx->method.result == NULL )    { return FNT_FAILURE; }
    /* extra can be NULL when method doesn't return a result. */

    int ret = ctx->method.result(extra);

    if( ret == FNT_DONE && fnt_verbose_level >= 3 ) {
        printf("DEBUG: Method '%s' has finished.\n", ctx->method.name);
    } else if( ret == FNT_FAILURE && fnt_verbose_level >= 1 ) {
        fprintf(stderr, "ERROR: Method completion check failed.\n");
    }

    return ret;
}
