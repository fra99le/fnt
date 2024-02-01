/*
 * fnt.h
 * fnt: Numerical Toolbox
 *
 * Copyright (c) 2024 Bryan Franklin. All rights reserved.
 */

#define FNT_SUCCESS 0
#define FNT_FAILURE 1
#define FNT_DONE    2

/** \brief Creates an opaque context handle.
 * \param context Pointer to a void* to be assigned to the context.
 * \param dimensions Number of elements in an input vector.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int fnt_init(void **context, char *path);

/** \brief Set the leverl of verbosity.
 * \param verbosity Setls the level of verbosity, higher values are more verbose (default 0).
 *      Level   Description
 *      0       No console output
 *      1       Errors only on stderr
 *      2       Info on stdout
 *      3       Debugging Info on stdout
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

/** \brief Provide initial inputs values.
 * \param context FNT context for method.
 * \param vec Pointer to input vector being seeded.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int fnt_seed(void *context, double *vec);

/** \brief Get next input vector to try
 * \param context FNT context for method.
 * \param vec Pointer to allocated input vector.
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int fnt_next(void *context, double *vec);

/** \brief Provide the value of the objective function for input vector.
 * \param context FNT context for method.
 * \param vec Pointer to the input vector (i.e. v).
 * \param value Value of objective function (i.e., f(v)).
 * \return FNT_SUCCESS on success, FNT_FAILURE otherwise.
 */
int fnt_set_value(void *context, double *vec, double value);

/** \brief Check if method had completed.
 * \param context FNT context to be checked.
 * \return FNT_DONE when complete, zero otherwise.
 */
int fnt_done(void *context);

/** \brief Produce final result from method,
 * \param context FNT context for which fnt_done returns FNT_DONE.
 */
int fnt_result(void *context);
