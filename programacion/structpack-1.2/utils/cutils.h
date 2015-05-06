/* cutils.h
 * Function prototypes for C utility functions.
 *
 * Date: 18 february 2012
 * Author: Pablo Martinez Naredo <pmnaredo@gmail.com>
 */

#include <stdio.h>

/* Macro definitions */
#define align_mem(x) x

/* Define boolean values if not already done */
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

int read_input(int n, double **t, FILE *file);

void print_vector(const char *name, double *v, int n);

void print_matrix(const char *name, double *m, int rows, int cols, int lda);

void ctimer_(double *elapsed, double *ucpu, double *scpu, double *gtime);
