/** @file csp_dpis.h
 *  @brief Header file for using the csp_dpis module from C source code
 *
 *  This file contains the suitable function prototypes for calling the
 *  csp_dpis module routines from C source code.
 *  @author Pablo Martinez Naredo
 */

#ifndef CSP_TPSYTRID_H
#define CSP_TPSYTRID_H

/** @brief Solves the linear system
 *        T x = b,
 *  where T is a tridiagonal symmetric Toeplitz matrix, using the Rojo method.
 *  @param n  Size of matrix.
 *  @param t0 Diagonal of the Toeplitz matrix.
 *  @param t1 Sub/Super-diagonal of the Toeplitz matrix.
 *  @param x  On input points to right hand side vector (b). On output points to the solution vector (x).
 */
void csp_dpis_rojosv(const int n, const double t0, const double t1,
			 double *x);

/** @brief Solves the linear system
 *        T x = b,
 *  where T is a tridiagonal symmetric Toeplitz matrix, using the DST method.
 *  @param n  Size of matrix.
 *  @param t0 Diagonal of the Toeplitz matrix.
 *  @param t1 Sub/Super-diagonal of the Toeplitz matrix.
 *  @param x  On input points to right hand side vector (b). On output points to the solution vector (x).
 */
void csp_dpis_dstsv(const int n, const double t0, const double t1,
			double *x);


/** @brief Solves the linear system
 *        T x = b,
 *  where T is a tridiagonal symmetric Toeplitz matrix, using the LDLt method.
 *  @param n  Size of matrix.
 *  @param t0 Diagonal of the Toeplitz matrix.
 *  @param t1 Sub/Super-diagonal of the Toeplitz matrix.
 *  @param x  On input points to right hand side vector (b). On output points to the solution vector (x).
 *  @return   Zero on success, non zero in case of error.
 */
int csp_dpis_ldltsv(const int n, const double t0, const double t1,
			double *x);


/** @brief Calculates the maximum SVD of a tridiagonal symmetric Toeplitz matrix.
 *  @param n  Size of matrix.
 *  @param t0 Diagonal of the Toeplitz matrix.
 *  @param t1 Sub/Super-diagonal of the Toeplitz matrix.
 *  @return   Maximum SVD of the Toeplitz matrix.
 */
double csp_dpis_maxsvd(const int n, const double t0, const double t1);


/** @brief Calculates the minimum SVD of a tridiagonal symmetric Toeplitz matrix.
 *  @param n  Size of matrix.
 *  @param t0 Diagonal of the Toeplitz matrix.
 *  @param t1 Sub/Super-diagonal of the Toeplitz matrix.
 *  @return   Minimum SVD of the Toeplitz matrix.
 */
double csp_dpis_minsvd(const int n, const double t0, const double t1);


/** @brief Computes the matrix-vector operation defined as
 *        b = b + alpha * T * x,
 *  where T is a tridiagonal symmetric Toeplitz matrix.
 *  @param n     Size of matrix.
 *  @param alpha Scalar alpha.
 *  @param t0    Diagonal of the Toeplitz matrix.
 *  @param t1    Sub/Super-diagonal of the Toeplitz matrix.
 *  @param x     Pointer to the first element of vector x.
 *  @param b     Pointer to the first element of vector b. Changed on output.
 */
void csp_dpis_gemv(const int n, const double alpha, const double t0,
		       const double t1, const double *x, double *b);


#endif /* CSP_TPSYTRID_H */
