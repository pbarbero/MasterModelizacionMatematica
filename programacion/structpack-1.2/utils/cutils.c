/* cutils.c
 * Utility functions for C.
 *
 * Date: 18 february 2012
 * Author: StructPack team
 */

#include <stdio.h>
#include <stdlib.h>


/**
 * Auxiliar function for reading input data from file.
 *
 * \param n Variable that will store input data size.
 * \param t Pointer to variable that will store input data vector.
 * \param file File where contents should be read from.
 * \return \c 1 if successful, \c 0 otherwise.
*/
int read_input(int n, double **t, FILE *file)
{
  /* Check size */
  if (n <= 0)
    return 0;

  /* Read vector */
  int i;
  *t = (double *) malloc(n * sizeof (double));
  for (i=0; i<n; ++i) {
    if (fscanf(file, "%lf", (*t) + i) != 1) {
      free(t);
      return 0;
    }
  }

  return 1;
}

/**
 * Print a given vector using standard error output.
 *
 * \param name Name of the vector (C-style string).
 * \param v Pointer to vector data.
 * \param n Size of the vector.
*/
void print_vector(const char *name, double *v, int n)
{
  int i;
  const int elems = 12;
  
  printf("%s:%s", name, n >= elems ? "\n" : "");
  for (i=0; i<n; ++i) {
    printf("%11.5g ", v[i]);
    if ((i+1) % elems == 0)
      printf("\n");
  }
  if (n >= elems && ((i+1) % elems != 0))
    printf("\n");
  printf("\n");
}

/**
 * Print a given matrix using standard error output.
 *
 * \param name Name of the matrix (C-style string).
 * \param m Pointer to matrix data (column ordering).
 * \param rows Number of rows.
 * \param cols Number of columns.
 * \param lda Leading dimension of the matrix.
*/
void print_matrix(const char *name, double *m, int rows, int cols, int lda)
{
  int i, j;
  printf("%s:\n", name);
  if (rows <= 0 || cols <= 0 || lda == 0)
    return;

  for (j=0; j<rows; ++j) {
    for (i=0; i<cols; ++i)
      printf("%9.5f ", m[i*lda+j]);
    printf("\n");
  }
  printf("\n");
}

