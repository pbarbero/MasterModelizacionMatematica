/** @file csp_dt.h
 *  @brief Header file for using the csp_dt module from C source code
 *
 *  This file contains the suitable function prototypes for calling the
 *  csp_dt module routines from C source code.
 *  @author Daniel Arguelles Martino
 */

#ifndef CSP_TPGESV_H
#define CSP_TPGESV_H

/** @brief Solves the linear system
 *        T x = b,
 *  where T is a general Toeplitz matrix.
 *  @param n    Size of matrix.
 *  @param u    Pointer to the first element of the first column of the
 *              Toeplitz matrix.
 *  @param v    Pointer to the first element of the row column of the
 *              Toeplitz matrix.               
 *  @param x    Pointer to the first element of solution vector.
 *  @param b    Pointer to the first element of the right hand side vector.
 *              Changed on output.
 *  @param nb   Block size.
 *  @param piv  Uses local diagonal pivoting if nonzero.
 *  @param ref  Number of iterative refinement steps.   
 *  @param INFO = 0: successful exit
 *              = 1: first element of u and first element of v are different
 *              = 2: u and v have different length
 *              = 3: u and b have different length
 */
void csp_dt_sv(const int n, const double *u, const double *v, double *x, double *b,
		       const int nb, const int piv, const int ref, int *INFO);

               
/** @brief Calculates the 1-norm of a general Toeplitz matrix.
 *  @param n Size of matrix.
 *  @param u Pointer to the first element of the first column of the Toeplitz matrix.
 *  @param v Pointer to the first element of the first row of the Toeplitz matrix.
 *  @return  The 1-norm of the Toeplitz matrix.
 */
double csp_dt_nrm1(const int n, const double *u, const double *v );



/** @brief Computes the matrix-vector operation defined as
 *        b = b + alpha * T * x,
 *  where T is a general Toeplitz matrix.
 *  @param n     Size of matrix.
 *  @param alpha Scalar alpha.
 *  @param u     Pointer to the first element of the first column of the Toeplitz matrix.
 *  @param v     Pointer to the first element of the first row of the Toeplitz matrix.
 *  @param x     Pointer to the first element of vector x.
 *  @param b     Pointer to the first element of vector b. Changed on output.
 */
void csp_dt_gemv(const int n, const double alpha, const double *u, const double *v,
		     const double *x, const double beta, double *b);


#endif
