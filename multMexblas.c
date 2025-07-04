#include <math.h>
#include "mex.h"
#include <omp.h>
#include "blas.h"

static void dotMex(double a[], double b[],unsigned int n, double*c) {
	
    omp_lock_t writelock;

    omp_init_lock(&writelock);
    
    int incx =1;
    int incy =1;
    
    ptrdiff_t m = (ptrdiff_t)n, one1 = 1, one2 = 1;
    
    //*c = ddot(&m,a,&one1,b,&one2);
    double result;
    size_t i=0,chunk;
    //chunk = n/20;
    #pragma omp parallel for num_threads(20) default(shared) private(i) reduction(+:result)  
    for (int i=0;i<n;i++)
    {
       
        result+= a[i]*b[i];

    }
        
	*c=result;
	
}


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
		//mexErrMsgIdAndTxt("MATLAB:yprime:invalidY", "Second input argument must be a real matrix.");
	    double* a,*b,*c;
		unsigned int* n;
		
		a = mxGetPr(prhs[0]);
		b = mxGetPr(prhs[1]);
		n = (unsigned int*)mxGetPr(prhs[2]);
		
		plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
		
		c = mxGetPr(plhs[0]);
		
		//c[0] = 10;
		dotMex(a, b,*n,c);
		return;
    
}