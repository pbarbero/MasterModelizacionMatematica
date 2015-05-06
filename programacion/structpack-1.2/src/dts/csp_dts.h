/** @file csp_dts.h
 *  @brief Header file for using the csp_dts module from C source code
 *
 *  This file contains the suitable function prototypes for calling the
 *  csp_dts module routines from C source code.
 *  @author Pablo Martinez Naredo
 */

#ifndef CSP_TPSYMM_H
#define CSP_TPSYMM_H

/** @brief Solves the linear system
 *        T x = b,
 *  where T is a symmetric Toeplitz matrix.
 *  @param n   Size of matrix.
 *  @param t   Pointer to the first element of the first column (row) of the
 *             Toeplitz matrix.
 *  @param x   Pointer to the first element of solution vector.
 *  @param b   Pointer to the first element of the right hand side vector.
 *             Changed on output.
 *  @param nb  Block size.
 *  @param piv Uses local diagonal pivoting if nonzero.
 *  @param aff Uses local thread/core affinity.
 */
void csp_dts_sv(const int n, const double *t, double *x, double *b,
		       const int nb, const int piv, const int aff );


/** @brief Calculates the 1-norm of a symmetric Toeplitz matrix.
 *  @param n Size of matrix.
 *  @param t Pointer to the first element of the first column (row) of the Toeplitz matrix.
 *  @return  The 1-norm of the Toeplitz matrix.
 */
double csp_dts_nrm1(const int n, const double *t);


/** @brief Computes the matrix-vector operation defined as
 *        b = b + alpha * T * x,
 *  where T is a symmetric Toeplitz matrix.
 *  @param n     Size of matrix.
 *  @param alpha Scalar alpha.
 *  @param t     Pointer to the first element of the first column (row) of the Toeplitz matrix.
 *  @param x     Pointer to the first element of vector x.
 *  @param b     Pointer to the first element of vector b. Changed on output.
 */
void csp_dts_gemv(const int n, const double alpha, const double *t,
		     const double *x, double *b);


#endif /* CSP_TPSYMM_H */
