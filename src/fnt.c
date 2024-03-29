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
#include <sys/errno.h>
#include "fnt.h"
#include "fnt_util.h"
#include "fnt_vect.h"

/* MARK: Internal constants */

#define FNT_MAX_NAME_LENGTH 64

/* MARK: Internal variables */

int fnt_verbose_level = FNT_WARN;   /* default to showing errors & warnings */

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
    void *handle;
    int (*init)(void **handle, int dimensions);
    int (*free)(void **handle);
    int (*info)();
    int (*hparam_set)(void *handle, char *id, void *value_ptr);
    int (*hparam_get)(void *handle, char *id, void *value_ptr);
    int (*next)(void *handle, fnt_vect_t *vec);
    int (*value)(void *handle, fnt_vect_t *vec, double value);
    int (*value_gradient)(void *handle, fnt_vect_t *vec, double value, fnt_vect_t *gradient);
    int (*done)(void *handle);
    int (*result)(void *handle, char *id, void *value_ptr);
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
        ERROR("calloc: %s\n", strerror(errno));
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
        size_t new_size = (ctx->methods_list.capacity*2 + 1);
        DEBUG("DEBUG: Resizing methods list from %d to %zu\n", ctx->methods_list.capacity, new_size);
        void *ptr = realloc(ctx->methods_list.entries,
                            new_size * sizeof(fnt_method_list_entry_t));
        if( ptr == NULL )   {
            ERROR("realloc: %s\n", strerror(errno));
            return FNT_FAILURE;
        }
        ctx->methods_list.entries = ptr;
    }

    int pos = ctx->methods_list.count;
    memcpy(&ctx->methods_list.entries[pos], entry, sizeof(fnt_method_list_entry_t));
    ++ctx->methods_list.count;

    return FNT_SUCCESS;
}


int fnt_method_list_print(context_t *ctx) {
    if( ctx == NULL )   { return FNT_FAILURE; }

    printf("%24s : %s\n","Method","File");
    printf("------------------------ : ------------------------\n");
    for(int i=0; i<ctx->methods_list.count; ++i) {
        printf("%24s : %s\n", ctx->methods_list.entries[i].name,
                              ctx->methods_list.entries[i].path);
    }
    printf("------------------------ : ------------------------\n");

    return FNT_SUCCESS;
}


int fnt_register_method(context_t *ctx, char *filename) {

    /* open object file */
    void *dl_handle = NULL;
    if( !(dl_handle = dlopen(filename,RTLD_NOW)) ) {
        ERROR("ERROR: dlopen: %s\n", dlerror());
        return FNT_FAILURE;
    }

    /* extract name function */
    int (*fnt_method_name)(char *, int) = dlsym(dl_handle, "method_name");
    if( fnt_method_name == NULL ) {
        ERROR("ERROR: No fnt_method_name in '%s'.\n", filename);
        dlclose(dl_handle);
        return FNT_FAILURE;
    }
    char name[FNT_MAX_NAME_LENGTH];
    fnt_method_name(name, FNT_MAX_NAME_LENGTH);

    /* set up list entry */
    fnt_method_list_entry_t *entry = calloc(1,sizeof(fnt_method_list_entry_t));
    snprintf(entry->name, sizeof(entry->name), "%s", name);
    snprintf(entry->path, sizeof(entry->path), "%s", filename);

    /* add method entry to list of available methods */
    fnt_method_list_add(ctx, entry);
    INFO("\tfound method '%s' in '%s'.\n", entry->name, filename);

    dlclose(dl_handle);

    return FNT_SUCCESS;
}


int fnt_register_methods(context_t *ctx, char *path) {

    /* check input */
    if( path == NULL ) {
        ERROR("ERROR: Methods path is NULL.\n");
        return FNT_FAILURE;
    }

    /* open directory */
    DIR *dirp = NULL;
    if( (dirp = opendir(path)) == NULL ) {
        ERROR("opendir: %s\n", strerror(errno));
        return FNT_FAILURE;
    }

    /* read contents of the directory */
    INFO("Loading methods from '%s'.\n", path);

    struct dirent *dp;
    char filename[PATH_MAX];
    while( (dp = readdir(dirp)) != NULL ) {

        /* check extension */
        int namelen = strlen(dp->d_name);
        if( namelen <= 3
            || dp->d_name[0] == '.'
            || strncasecmp(dp->d_name+namelen-3, ".so", 3) ) {
            DEBUG("DEBUG: Skipping '%s'.\n", dp->d_name);
            continue;
        }

        /* construct filename */
        snprintf(filename,sizeof(filename), "%s/%s", path, dp->d_name);

        /* record method name and filename */
        if( fnt_register_method(ctx, filename) == FNT_FAILURE ) {
            ERROR("ERROR: Failed to load '%s' from '%s'.\n", dp->d_name, path);
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
        ERROR("ERROR: dlopen: %s\n", dlerror());
        return FNT_FAILURE;
    }

    INFO("Loading method from '%s'.\n", filename);

    /* assign function pointers */
    ctx->dl_handle = dl_handle;
    int (*method_name)(char *, int) = dlsym(dl_handle, "method_name");
    method_name(ctx->method.name, sizeof(ctx->method.name));
    ctx->method.init = dlsym(dl_handle, "method_init");
    ctx->method.free = dlsym(dl_handle, "method_free");
    ctx->method.info = dlsym(dl_handle, "method_info");
    ctx->method.hparam_get = dlsym(dl_handle, "method_hparam_get");
    ctx->method.hparam_set = dlsym(dl_handle, "method_hparam_set");
    ctx->method.next = dlsym(dl_handle, "method_next");
    ctx->method.value = dlsym(dl_handle, "method_value");
    ctx->method.value_gradient = dlsym(dl_handle, "method_value_gradient");
    ctx->method.done = dlsym(dl_handle, "method_done");
    ctx->method.result = dlsym(dl_handle, "method_result");

    if( ctx->method.next == NULL
        || ctx->method.value == NULL
        || ctx->method.done == NULL ) {
        ERROR("ERROR: '%s' does not have all required methods.\n", filename);
        if( ctx->method.next == NULL )
            ERROR("\tMISSING method_next(void*, fnt_vect_t*)\n");
        if( ctx->method.value == NULL )
            ERROR("\tMISSING method_value(void*, fnt_vect_t*, double)\n");
        if( ctx->method.done == NULL )
            ERROR("\tMISSING method_done(void*)\n");
        memset(&ctx->method, '\0', sizeof(ctx->method));
        dlclose(dl_handle); ctx->dl_handle = dl_handle = NULL;

        return FNT_FAILURE;
    }

    return FNT_SUCCESS;
}

/* MARK: User callable functions */

int fnt_init(void **context, char *path) {
    context_t *ctx = calloc(1,sizeof(context_t));
    *context = (void*)ctx;

    /* load list of available methods */
    int ret = fnt_register_methods(ctx, path);

    if( fnt_verbose_level >= FNT_DEBUG ) {
        /* print methods/file mapping */
        fnt_method_list_print(ctx);
    }

    return ret;
}


int fnt_verbose(int verbosity) {
    fnt_verbose_level = verbosity;
    INFO("Verbosity set to %i.\n", fnt_verbose_level);
    return FNT_SUCCESS;
}


int fnt_set_method(void *context, char *name, int dimensions) {
    context_t *ctx = (context_t*)context;
    if( ctx == NULL )   { return FNT_FAILURE; }
    if( name == NULL )  { return FNT_FAILURE; }

    ctx->dim = dimensions;
    INFO("Initializing method '%s' for %d dimensions.\n", name, ctx->dim);

    /* dynamically load module */
    for(int i=0; i<ctx->methods_list.count; ++i) {
        DEBUG("DEBUG: checking %s\n", ctx->methods_list.entries[i].name);

        if( strncmp(ctx->methods_list.entries[i].name, name,
                sizeof(ctx->methods_list.entries[i].name)) == 0 ) {
            int ret = fnt_method_load(ctx, ctx->methods_list.entries[i].path);

            if( ret == FNT_SUCCESS ) {
                INFO("Loaded method '%s'.\n", ctx->method.name);
            } else if( ret == FNT_FAILURE ) {
                ERROR("ERROR: Loading method '%s' failed..\n", ctx->method.name);
                continue;   /* keep looking for one that might work */
            }

            ret = ctx->method.init(&ctx->method.handle, dimensions);

            if( ret == FNT_SUCCESS ) {
                INFO("Initialized method '%s' for %i dimensional inputs.\n", ctx->method.name, dimensions);
            } else if( ret == FNT_FAILURE ) {
                ERROR("ERROR: Initialization of method '%s' failed..\n", ctx->method.name);
                continue;   /* keep looking for one that might work */
            }

            return ret;
        }
    }
    ERROR("Failed to find method '%s'.\n", name);

    return FNT_FAILURE;
}


int fnt_free(void **context) {
    if( context == NULL )   { return FNT_FAILURE; }
    context_t *ctx = (context_t*)*context;
    if( ctx == NULL )       { return FNT_FAILURE; }

    int ret = FNT_SUCCESS;
    if( ctx->method.free != NULL )  {
        ret = ctx->method.free(&ctx->method.handle);

        if( ret == FNT_SUCCESS ) {
            DEBUG("DEBUG: Freed internally allocated values for method '%s'.\n", ctx->method.name);
        } else if( ret == FNT_FAILURE ) {
            ERROR("ERROR: Call to free for method '%s' failed..\n", ctx->method.name);
        }
    }

    fnt_method_list_free(&ctx->methods_list);
    dlclose(ctx->dl_handle);    ctx->dl_handle = NULL;

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
        ERROR("ERROR: Called %s before setting method.\n", __FUNCTION__);
        return FNT_FAILURE;
    }
    if( ctx->method.info == NULL ) {
        ERROR("ERROR: Method '%s' does not provide additional info.\n", ctx->method.name);
        return FNT_FAILURE;
    }

    return ctx->method.info();
}


int fnt_hparam_set(void *context, char *id, void *value_ptr) {
    if( context == NULL )                   { return FNT_FAILURE; }
    context_t *ctx = (context_t*)context;
    if( ctx->method.hparam_set == NULL )    { return FNT_FAILURE; }
    if( id == NULL )                        { return FNT_FAILURE; }
    if( value_ptr == NULL )                 { return FNT_FAILURE; }

    int ret = ctx->method.hparam_set(ctx->method.handle, id, value_ptr);

    if( ret == FNT_SUCCESS ) {
        INFO("Set hyper-parameter '%s'.\n", id);
    } else if( ret == FNT_FAILURE ) {
        ERROR("ERROR: Failed to set hyper-parameter '%s'.\n", id);
    }

    return ret;
}


int fnt_hparam_get(void *context, char *id, void *value_ptr) {
    context_t *ctx = (context_t*)context;
    if( ctx == NULL )                       { return FNT_FAILURE; }
    if( ctx->method.hparam_get == NULL )    { return FNT_FAILURE; }
    if( id == NULL )                        { return FNT_FAILURE; }
    if( value_ptr == NULL )                 { return FNT_FAILURE; }

    int ret = ctx->method.hparam_get(ctx->method.handle, id, value_ptr);

    if( ret == FNT_SUCCESS ) {
        INFO("Got hyper-parameter '%s'.\n", id);
    } else if( ret == FNT_FAILURE ) {
        ERROR("ERROR: Failed to get hyper-parameter '%s'.\n", id);
    }

    return ret;
}


int fnt_next(void *context, fnt_vect_t *vec) {
    context_t *ctx = (context_t*)context;
    if( ctx == NULL )               { return FNT_FAILURE; }
    if( ctx->method.next == NULL )  { return FNT_FAILURE; }
    if( vec == NULL )               { return FNT_FAILURE; }

    int ret =  ctx->method.next(ctx->method.handle, vec);

    if( ret == FNT_SUCCESS ) {
        if( fnt_verbose_level >= FNT_DEBUG ) {
            fnt_vect_println(vec, "DEBUG: Retrieved next input vector: ", NULL);
        }
    } else if( ret == FNT_FAILURE ) {
        ERROR("ERROR: Failed to retrieve next input vector.\n");
    }

    return ret;
}


int fnt_set_value(void *context, fnt_vect_t *vec, double value) {
    context_t *ctx = (context_t*)context;
    if( ctx == NULL )               { return FNT_FAILURE; }
    if( ctx->method.value == NULL ) { return FNT_FAILURE; }
    if( vec == NULL )               { return FNT_FAILURE; }
    if( vec->v == NULL )            { return FNT_FAILURE; }

    int ret = ctx->method.value(ctx->method.handle, vec, value);

    if( ret == FNT_SUCCESS ) {
        if( fnt_verbose_level >= FNT_DEBUG ) {
            DEBUG("DEBUG: Set value of objective function");
            fnt_vect_print(vec, " for input ", "%.2f");
            DEBUG(" to %g.\n", value);
        }
    } else if( ret == FNT_FAILURE ) {
        ERROR("ERROR: Failed to set objective value for input vector.\n");
    }

    return ret;
}


int fnt_set_value_gradient(void *context, fnt_vect_t *vec, double value, fnt_vect_t *gradient) {
    context_t *ctx = (context_t*)context;
    if( ctx == NULL )               { return FNT_FAILURE; }
    if( vec == NULL )               { return FNT_FAILURE; }
    if( gradient == NULL )          { return FNT_FAILURE; }

    /* fall back to non-gradient function, if gradient version not supplied. */
    if( ctx->method.value_gradient == NULL ) {
        return fnt_set_value(context, vec, value);
    }

    int ret = ctx->method.value_gradient(ctx->method.handle, vec, value, gradient);

    if( ret == FNT_SUCCESS ) {
        if( fnt_verbose_level >= FNT_DEBUG ) {
            DEBUG("DEBUG: Set value of objective function");
            fnt_vect_print(vec, " for input ", "%.2f");
            DEBUG(" to %g.\n", value);
        }
    } else if( ret == FNT_FAILURE ) {
        ERROR("ERROR: Failed to set objective value for input vector.\n");
    }

    return ret;
}


int fnt_done(void *context) {
    context_t *ctx = (context_t*)context;
    if( ctx == NULL )               { return FNT_FAILURE; }
    if( ctx->method.done == NULL )  { return FNT_FAILURE; }

    int ret = ctx->method.done(ctx->method.handle);

    if( ret == FNT_DONE ) {
        DEBUG("DEBUG: Method '%s' has finished.\n", ctx->method.name);
    } else if( ret == FNT_FAILURE ) {
        ERROR("ERROR: Method completion check failed.\n");
    }

    return ret;
}


int fnt_result(void *context, char *name, void *value_ptr) {
    context_t *ctx = (context_t*)context;
    if( ctx == NULL )                   { return FNT_FAILURE; }
    /* method is optional, so return success if it is not provided. */
    if( ctx->method.result == NULL )    { return FNT_SUCCESS; }
    /* extra can be NULL when method doesn't return a result. */

    if( fnt_done(context) != FNT_DONE ) {
        DEBUG("DEBUG: Method '%s' has not finished yet.\n", ctx->method.name);
        return FNT_FAILURE;
    }

    int ret = ctx->method.result(ctx->method.handle, name, value_ptr);

    if( ret == FNT_FAILURE ) {
        ERROR("ERROR: Method result reporting failed.\n");
    }

    return ret;
}
