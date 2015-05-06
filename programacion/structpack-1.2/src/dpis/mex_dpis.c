/* mex_dpis.c
 * Source code for the mex_dpis.mex file. This MEX file is intended to be
 * invoked the following way:
 *
 *   x = mex_dpis(n, t0, t1, b, method)
 */

#include "mex.h" 
#include <csp_dpis.h>
#include <string.h>


/* Entry point */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int     n = (int) mxGetScalar(prhs[0]);
  double t0 = mxGetScalar(prhs[1]);
  double t1 = mxGetScalar(prhs[2]);
  double *b = mxGetPr(prhs[3]);
  double *x;
  char *method = mxArrayToString(prhs[4]);

  plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);

  x = mxGetPr(plhs[0]); 
  memcpy(x, b, n*sizeof (double));

  if (!strcmp(method,"rojo")) {
    csp_dpis_rojosv(n, t0, t1, x);
  }
  else 
    if (!strcmp(method, "dst")) {
      csp_dpis_dstsv(n, t0, t1, x);
    }
    else
      if (!strcmp(method,"ldlt")) {
	csp_dpis_ldltsv(n, t0, t1, x);
      }
      else
	mexPrintf("error: method %s not recognised.\n", method);
} 
