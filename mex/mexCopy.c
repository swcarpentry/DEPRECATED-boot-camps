#include <mex.h>

// 
// In all mex functions, *lhs are the matlab "outputs" and *rhs are the 
// matlab "inputs". Think argc/argv**
//

// you MUST name it this. The file name is the callable name.
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	int i,j,m,n;
	double *dataIn, *dataOut;
	if (nrhs != nlhs)
		mexErrMsgTxt("The number of inputs and outputs had better be the same!");

	// NOTE: if the last line was run, this function has returned.
	
	for (i=0; i < nrhs; i++){ // for each input
		m = mxGetM(prhs[i]); // get size 1
		n = mxGetN(prhs[i]); // get size 2

		plhs[i] = mxCreateDoubleMatrix(m,n,mxREAL); // create an m by n real matrix.
		dataIn = mxGetPtr(plhs[i]); // pointer to input.
		dataOut = mxGetPtr(prhs[i]); // pointer to ouput (for convienence)
		for (j = 0; j < m*n; j++)
			dataOut[j] = 2 * dataIn[j]; // stupid... need better example. 

	}

}


