/*
 * fnt.h
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */
#ifndef FNT_H
#define FNT_H

#include "fnt_util.h"
#include "fnt_vect.h"

/** \brief Creates an opaque context handle.
 * \param context Pointer to a void* to be assigned to the context.
 * \param dimensions Number of elements in an input vector.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int fnt_init(void **context, char *path);

/** \brief Set the leverl of verbosity.
 * \param verbosity Setls the level of verbosity, higher values are more verbose (default 0).
 *      Level       Description
 *      FNT_NONE    No console output
 *      FNT_ERROR   Errors only on stderr
 *      FNT_WARN    Errors and Warnings on stderr
 *      FNT_INFO    Info on stdout
 *      FNT_DEBUG   Debugging Info on stdout
 * \return Will always return FNT_SUCCESS.
 */
int fnt_verbose(int verbosity);

/** \brief Choose the numeric method to use.
 * \param context FNT context for method.
 * \param name Name of the method being selected.
 * \param dimension Number of dimensions in the input vector.
 */
int fnt_set_method(void *context, char *name, int dimensions);

/** \brief Frees an FNT context.
 * \param context Pointer to the void* to be freed.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int fnt_free(void **context);

/** \brief Display method info, if available.
 * \return FNT_SUCCESS if info available, FNT_FAILURE otherwise.
 */
int fnt_info(void *context);

/** \brief Provide hyper-parameters the method may need.
 * \param id Name of the hyper-parameter.
 * \param value_ptr pointer to the value being set.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int fnt_hparam_set(void *context, char *id, void *value_ptr);

/** \brief Retrieve hyper-parameters from the method.
 * \param id Name of the hyper-parameter.
 * \param value_ptr pointer to the value being set.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int fnt_hparam_get(void *context, char *id, void *value_ptr);

/** \brief Provide initial inputs values.
 * \param context FNT context for method.
 * \param vec Pointer to input vector being seeded.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int fnt_seed(void *context, fnt_vect_t *vec);

/** \brief Get next input vector to try
 * \param context FNT context for method.
 * \param vec Pointer to allocated input vector to be filled in.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int fnt_next(void *context, fnt_vect_t *vec);

/** \brief Provide the value of the objective function for input vector.
 * \param context FNT context for method.
 * \param vec Pointer to the input vector (i.e. v).
 * \param value Value of objective function (i.e., f(v)).
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int fnt_set_value(void *context, fnt_vect_t *vec, double value);

/** \brief Check if method had completed.
 * \param context FNT context to be checked.
 * \return FNT_DONE when complete, zero otherwise.
 */
int fnt_done(void *context);

/** \brief Get best input vector (i.e., the one that produced the lowest objective function value).
 * \param context FNT context for method.
 * \param vec Pointer to allocated input vector to be filled in.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int fnt_best(void *context, fnt_vect_t *vec);

/** \brief Produce final result from method,
 * \param context FNT context for which fnt_done returns FNT_DONE.
 * \param extra Method dependant.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int fnt_result(void *context, void *extra);

#endif /* FNT_H */
