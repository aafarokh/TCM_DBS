#include <math.h>
#include "mex.h"
#include <omp.h>


static void multMex(double a[], double b[], double c[],int indb[],int n) {
	
	#pragma omp parallel for num_threads(20)
	for(int i=0;i<n;i++)
		c[i] = a[i]*b[indb[i]-1];
	
	
}


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
		//mexErrMsgIdAndTxt("MATLAB:yprime:invalidY", "Second input argument must be a real matrix.");
	    double* a,*b,*c;
		double *n;
		int * indb;
		
		a = mxGetPr(prhs[0]);
		b = mxGetPr(prhs[1]);
		n = mxGetPr(prhs[2]);
		indb = (int*)mxGetPr(prhs[3]);
		
		plhs[0] = mxCreateDoubleMatrix(*n, 1, mxREAL);
		
		c = mxGetPr(plhs[0]);
		
		//c[0] = 10;
		multMex(a, b, c,indb,*n);
		return;
    
}