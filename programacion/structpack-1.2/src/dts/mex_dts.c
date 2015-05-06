/* mex_dts.c
 * Source code for the mex_dts.mex file. This MEX file is intended to be
 * invoked the following way:
 *
 *   x = mex_dts(t, b, nb, piv)
 */

#include "mex.h"
#include <csp_dts.h>

/* Entry point */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int     n = mxGetN(prhs[0])*mxGetM(prhs[0]);
  double *t = mxGetPr(prhs[0]);
  double *x;
  double *b = mxGetPr(prhs[1]);
  int    nb = (int) mxGetScalar(prhs[2]);
  int   piv = (int) mxGetScalar(prhs[3]);
  
  plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);

  x = mxGetPr(plhs[0]); 
  
  csp_dts_sv(n, t, x, b, nb, piv);
} 
