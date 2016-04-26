/**************************************************************************
%* 
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>
%*
%* This program is free software: you can redistribute it and/or modify
%* it under the terms of the GNU General Public License as published by
%* the Free Software Foundation, either version 3 of the License, or
%* (at your option) any later version.
%*
%* This program is distributed in the hope that it will be useful,
%* but WITHOUT ANY WARRANTY; without even the implied warranty of
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%* GNU General Public License for more details.
%*
%* You should have received a copy of the GNU General Public License
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.
%* 
%*************************************************************************/

#include "mex.h"
#include "matrix.h"
#include <math.h>       // contains lgamma
#include "lapack.h"

double r8_beta ( double x, double y ) {
    // from: http://people.sc.fsu.edu/~jburkardt/c_src/prob/prob.c
    double value;
    if ( x <= 0.0 || y <= 0.0 ) {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "R8_BETA - Fatal error!\n" );
        fprintf ( stderr, "  Both X and Y must be greater than 0.\n" );
        exit ( 1 );
    }
    value = exp ( 
        lgamma ( x ) 
      + lgamma ( y ) 
      - lgamma ( x + y ) );
    return value;
}

double beta_pdf ( double x, double a, double b ) {
    // computes the beta distribution pdf value at x with parameters a,b
    // from: http://people.sc.fsu.edu/~jburkardt/c_src/prob/prob.c
    double pdf;
    if ( x < 0.0 || 1.0 < x ) {
        pdf = 0.0;
    }
    else {
        pdf = pow ( x, ( a - 1.0 ) ) * pow ( ( 1.0 - x ), ( b - 1.0 ) ) / r8_beta ( a, b );
    }
    return pdf;
}

// The actual computational routine
void empcopulapdf(mwSize M, mwSize D, double *U, mwSize K, double h, double *outMatrix, mwSize outMatLen) {
    mwSignedIndex gridPointsDims[2]; gridPointsDims[0] = 1; gridPointsDims[1] = D;
    mwSignedIndex Kernel_vecDims[2]; Kernel_vecDims[0] = M; gridPointsDims[1] = 1;
    mxArray *gridPoints = mxCreateNumericArray(2, gridPointsDims, mxINT32_CLASS, mxREAL);
    mxArray *Kernel_vec = mxCreateNumericArray(2, Kernel_vecDims, mxDOUBLE_CLASS, mxREAL);
    //for(int ii=0; ii<outMatLen; ii++) {
    //    
    //}
    printf("HELLO :)\n");
    mxDestroyArray(gridPoints);
    mxDestroyArray(Kernel_vec);
}

// The gateway function
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]) {
    double *inMatrix;
    size_t M, D;
    mwSize K;       /* grid spacing */
    double h;       /* h - the kernel bandwidth */
    
    double *outMatrix;
    
    /* check for proper number of arguments */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("Copula:arrayProduct:nrhs","Three inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("Copula:arrayProduct:nlhs","One output required.");
    }
    
    /* make sure the first input argument is type double */
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("Copula:arrayProduct:notDouble","Input matrix must be type double.");
    }
    
    /* make sure the second input argument is scalar */
    if( mxIsComplex(prhs[1]) ||
         mxGetNumberOfElements(prhs[1])!=1 ) {
        mexErrMsgIdAndTxt("Copula:arrayProduct:notScalar","Input K must be a scalar.");
    }
    
    /* make sure the third input argument is scalar */
    if( mxIsComplex(prhs[2]) ||
         mxGetNumberOfElements(prhs[2])!=1 ) {
        mexErrMsgIdAndTxt("Copula:arrayProduct:notScalar","Input h must be a scalar.");
    }
    
    K = (int)mxGetScalar(prhs[1]);
    h = mxGetScalar(prhs[2]);
    
    inMatrix = mxGetPr(prhs[0]);
    M = mxGetM(prhs[0]);
    D = mxGetN(prhs[1]);
    
    /* create the output matrix */
    mwSize outMatLen = pow(K,D);
    plhs[0] = mxCreateDoubleMatrix(1,outMatLen,mxREAL);
    outMatrix = mxGetPr(plhs[0]);
    
    /* call the computational routine */
    empcopulapdf(M, D, inMatrix,K,h,outMatrix,outMatLen);
}