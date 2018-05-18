/**************************************************************\
*	Columnwise array multiplication between Matrix and Vector  *
*   For myMatAD                                                *
\**************************************************************/

#include "mex.h"
#define	Out	plhs[0]
#define	Vec	prhs[0]
#define	Mat	prhs[1]

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{	
	if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }
	switch(nrhs) {
    default: {
    	mexErrMsgTxt("Incorrect number of arguments.");
        }
    case 0: {
        mexPrintf("\nColumnwise dot-product matrix with vector.\n\n");
        break;
        }
	case 2: {
		int m=mxGetN(Vec);
		int n=mxGetN(Mat);
		Out=mxCreateDoubleMatrix(0,0,mxREAL);
		if (m > 0 && n > 0) {
            int i, j, k=0;
			const int* e=mxGetDimensions(Mat);
			double* outMat;
            double* inMat=mxGetPr(Mat);
            double* inVec=mxGetPr(Vec);
            mxSetDimensions(Out,e,3);
    		mxSetPr(Out,outMat=(double*) mxMalloc(e[0]*e[1]*e[2]*sizeof(double)));
            if (m == 1) { //scalar time matrix
                for (i=0; i<e[0]*e[1]*e[2]; i++) {
                        outMat[i] = inMat[i]*inVec[0];
                }
            } else { // 2D-matrix times 3D
                for (i=0; i<e[2]; i++) {
                    for (j=0; j<e[0]*e[1]; j++) {
                        outMat[k] = inMat[k]*inVec[j];
                        k++;
                    }
                }
            }
		}
        }
	}
}
