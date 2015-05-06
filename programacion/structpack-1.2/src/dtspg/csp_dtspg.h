/** @file csp_dtspg.h
 *  @brief Header file for using the csp_dtspg module from C source code
 *
 *  This file contains the suitable function prototypes for calling the
 *  csp_dtspg module routines from C source code.
 *  @author Daniel Arguelles Martino
 */


#ifndef CSP_DTSBG_H
#define CSP_DTSBG_H

/** @brief Solves the linear system
  *   TX=B
  * where T is a symmetric block Toeplitz matrix.
  * @param m    Size of matrix T.
  * @param n    Block size of matrix T.
  * @param T    Pointer to first element of first row block of matrix T.
  * @param B    Pointer to first element of right hand side matrix.
  * @param X    Pointer to first element of solution matrix.
  * @param k    Number of colums of matrix B and matrix X.
  */
void csp_dtspg_sv( const int m, const int n, double *T, double *B, double *X, const int k);

/** @brief Build a random symmetric block Toeplitz matrix.
 *  @param nBlocks  Number of blocks of matrix T.
 *  @param bSize    Block's size of matrix T.
 *  @param seed     Seed for random generation.
 *  @param T        Pointer to first element of matrix T.
 */
void csp_dtspg_rg( const int nBlocks, const int bSize, const int seed, double *T );


/** @brief Calculates the 1-norm of a symmetric block Toeplitz matrix.
 *  @param nBlocks  Number of blocks of matrix T.
 *  @param bSize    Block's size of matrix T.
 *  @param T        Pointer to first element of first row block of matrix T.
 *  @return         The 1-norm of the Toeplitz matrix.
 */
double csp_dtspg_nrm1( const int nBlocks, const int bSize, double *T );


/** @brief Computes the matrix-vector operation defined as
 *        b = beta* b + alpha * T * x,
 *  where T is a symmetric block Toeplitz matrix.
 *  @param nBlocks  Number of blocks of matrix T.
 *  @param bSize    Block's size of matrix T.
 *  @param alpha    Scalar alpha.
 *  @param T        Pointer to first element of first row block of matrix T.
 *  @param x        Pointer to the first element of vector x.
 *  @param beta     Scalar beta.
 *  @param b        Pointer to the first element of vector b. Changed on output.
 */
void csp_dtspg_gemv ( const int nBlocks, const int bSize, const double alfa, const double *T, const double *x, const double beta, double *b);

#endif
