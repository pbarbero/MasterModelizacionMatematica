/* dt_example.c
 * A simple C program that solves a general non symmetric Toeplitz system using
 * the csp_dt_sv function.
 *
 * Author: Daniel Arguelles Martino <daniel.arguelles.martino@gmail.com>
 * Date: 1 August 2012
 */

/* Solves the system
 *
 *  Tx = b
 *
 * where
 *
 *  T =
 *      (1  2  3  4)
 *      (5  1  2  3)
 *      (6  5  1  2)
 *      (7  6  5  1)
 *
 *  b =
 *      (1  1  1  1)
 */
 
 
 
#include <structpack.h>
#include <stdio.h>

#define N 4
#define NB 20
#ifndef TRUE
#define TRUE 1
#endif

int main(int argc, char **argv)
{
  const double u[N] = {1., 2., 3., 4.};
  const double v[N] = {1., 5., 6., 7.};
  double b[N] = {1., 1., 1., 1.};
  double x[N];
  int i;
  int INFO;
  
  
  printf("u:");
  for (i = 0; i < N; i++)
    printf(" %11.5g", u[i]);
  printf("\n");
  
  printf("v:");
  for (i = 0; i < N; i++)
    printf(" %11.5g", v[i]);
  printf("\n");

  printf("b:");
  for (i = 0; i < N; i++)
    printf(" %11.5g", b[i]);
  printf("\n");

  csp_dt_sv(N, u, v, x, b, NB, TRUE, 1, &INFO);

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