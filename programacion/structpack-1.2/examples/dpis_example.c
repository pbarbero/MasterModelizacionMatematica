/* dpis_example.c
 * A simple C program that solves a tridiagonal symmetric Toeplitz system using
 * the csp_dpis_rojosv function.
 *
 * Author: Pablo Martinez Naredo <pmnaredo@gmail.com>
 * Date: 27 june 2011
 */

/* Solves the system
 *
 *  Tx = b
 *
 * where
 *
 *  T =
 *      (1  2  0  0)
 *      (2  1  2  0)
 *      (0  2  1  2)
 *      (0  0  2  1)
 *
 *  b =
 *      (1  1  1  1)
 */

#include <structpack.h>
#include <csp_dpis.h>

#include <stdio.h>

#define N 4
#define NB 20
#ifndef TRUE
#define TRUE 1
#endif

int main(int argc, char **argv)
{
  const double t0 = 1.;
  const double t1 = 2.;
  double x[N] = {1., 1., 1., 1.};
  int i;

  printf("t0: %11.5g\n", t0);
  printf("t1: %11.5g\n", t1);
  printf("b:");
  for (i = 0; i < N; i++)
    printf(" %11.5g", x[i]);
  printf("\n");

  csp_dpis_rojosv(N, t0, t1, x);

  printf("x:");
  for (i = 0; i < N; i++)
    printf(" %11.5g", x[i]);
  printf("\n");

  return 0;
}
