/* dts_example.c
 * A simple C program that solves a general symmetric Toeplitz system using
 * the csp_dts_sv function.
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
 *      (1  2  3  4)
 *      (2  1  2  3)
 *      (3  2  1  2)
 *      (4  3  2  1)
 *
 *  b =
 *      (1  1  1  1)
 */

#include <structpack.h>
#include <csp_dts.h>

#include <stdio.h>

#define N 4
#define NB 20
#ifndef TRUE
#define TRUE 1
#endif

int main(int argc, char **argv)
{
  const double t[N] = {1., 2., 3., 4.};
  double b[N] = {1., 1., 1., 1.};
  double x[N];
  int i;

  printf("t:");
  for (i = 0; i < N; i++)
    printf(" %11.5g", t[i]);
  printf("\n");

  printf("b:");
  for (i = 0; i < N; i++)
    printf(" %11.5g", b[i]);
  printf("\n");

  csp_dts_sv(N, t, x, b, NB, TRUE);

  printf("x:");
  for (i = 0; i < N; i++)
    printf(" %11.5g", x[i]);
  printf("\n");

  printf("b:");
  for (i = 0; i < N; i++)
    printf(" %11.5g", b[i]);
  printf("\n");

  return 0;
}
