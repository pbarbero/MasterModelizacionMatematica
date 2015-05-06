/* dtspg_example.c
 * A simple C program that solves a general symmetric block Toeplitz system using
 * the csp_dtspg_sv function.
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
 *   (  71, 101,  23,  30 )
 *   ( 101, 194,  70,  96 )
 *   (  23,  70,  71, 101 )
 *   (  30,  96, 101, 194 )
 *
 *  b =
 *      (1  1  1  1)
 */
 
 
 
#include <structpack.h>
#include <stdio.h>

#define N 4
#define NB 2
#ifndef TRUE
#define TRUE 1
#endif

int main(int argc, char **argv)
{
  double T[NB*N] = {71.0, 101.0, 101.0,194.0, 23.0, 70.0, 30.0, 96.0};
  double b[N] = {1., 1., 1., 1.};
  double x[N];
  int i;
  
  
  printf("T = \n\
    (  71, 101,  23,  30 ) \n\
    ( 101, 194,  70,  96 ) \n\
    (  23,  70,  71, 101 ) \n\
    (  30,  96, 101, 194 )");
  printf("\n");
  
  printf("b:");
  for (i = 0; i < N; i++)
    printf(" %11.5g", b[i]);
  printf("\n");


  csp_dtspg_sv(N, NB, T, b, x, 1);

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