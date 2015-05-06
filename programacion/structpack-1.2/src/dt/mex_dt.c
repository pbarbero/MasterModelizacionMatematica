/**
  * x = dt(u, v, b, nb, piv )
  */

#include "mex.h"
#include <csp_dt.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int     n = mxGetN(prhs[0])*mxGetM(prhs[0]);
  double *u = mxGetPr(prhs[0]);
  double *v = mxGetPr(prhs[1]);
  double *x;
  double *b = mxGetPr(prhs[2]);
  int    nb = (int) mxGetScalar(prhs[3]);
  int   piv = (int) mxGetScalar(prhs[4]);
  int INFO =0;
  plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);

  x = mxGetPr(plhs[0]); 
  
  csp_dt_sv(n, u, v, x, b, nb, piv, 1, &INFO);
} 
