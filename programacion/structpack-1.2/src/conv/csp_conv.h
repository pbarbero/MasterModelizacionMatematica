/** @file csp_conv.h
 *  @brief Header file for using the csp_conv module from C source code
 *   
 *  This file contains the suitable function prototypes for calling the
 *  csp_conv module routines from C source code.
 *  @author StructPack Team
 */

#ifndef CSP_CONV_H
#define CSP_CONV_H

    /** @brief Non-Symmetric convolution 
     *       y = beta * y + alpha * T * x, 
     *  where T is   non-symmetric Toeplitz marix. 
     *  @param n     Order of the matrix.
     *  @param alpha alpha Scalar alpha. 
     *  @param t     Pointer to the first element of the first column of the Toeplitz matrix.
     *  @param u     Pointer to the first element of the first row of the Toeplitz matrix.
     *  @param x     Pointer to the first element of vector x.
     *  @param beta  beta Scalar beta. 
     *  @param y     Pointer to the first element of vector y.
     */
    void csp_convn( int n, double alpha, double *t, double *u, double *x, double beta, double *y );

    /** @brief Symmetric convolution 
     *       y = beta * y + alpha * T * x, 
     *  where T is symmetric Toeplitz marix. 
     *  @param n     Order of the matrix.
     *  @param alpha alpha Scalar alpha. 
     *  @param t     Pointer to the first element of the first column (row) of the Toeplitz matrix.
     *  @param x     Pointer to the first element of vector x.
     *  @param beta  beta Scalar beta. 
     *  @param y     Pointer to the first element of vector y.
     */
    void csp_convs( int n, double alpha, double *t, double *x, double beta, double *y );

#endif /* CSP_CONV_H */
