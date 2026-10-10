/*=================================================================
 *
 * ix2dmin.c	.MEX file corresponding to ix2dmin.m
 *	        	Takes 3 1D vectors x,y,z, finds the elements of  y that
 						are relative minimum of the neighboring 8 points
 *
 * The calling syntax is:
 *
 *		y = ix2dmin(degrees, x,y,z);   % degrees is the max number of minima
 *                                   % It is used in the m file but not in this mex file. 
 *
 *  You may also want to look at the corresponding M-code, ix2dmin.m.
 *
 * Author: Ryan King 12-10-05
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"

/* Input Arguments */

#define DEGREES phrs[0]
#define X_IN    prhs[1]
#define Y_IN	 prhs[2]
#define Z_IN    prhs[3]

/* Output Arguments */

#define	OUT     plhs[0]

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif


/* x,y,z vectors and n=length */
static void is2dmin(double *x, double *y, double *z, double *out, int n) {
  int j,next;
  double thisy = y[0];
  next=0;
  if ( thisy < y[1] && thisy < x[0] && thisy < x[1] && thisy < z[0] && thisy < z[1] ) out[next++]=1;

  for (j=1; j < (n-1); j++) {
    thisy=y[j];
    if (thisy < y[j-1] && thisy < y[j+1] && thisy < x[j-1] && thisy < x[j] && thisy < x[j+1] && thisy < z[j-1] && thisy < z[j] && thisy < z[j+1]){
		out[next++] = j+1;
	}
  }

  thisy=y[n-1];
  if (thisy < y[n-2] && thisy < x[n-1] && thisy < x[n-2] && thisy < z[n-2] && thisy < z[n-1]) out[next++]=n;
}


static unsigned int count2dmin(double *x, double *y, double *z, int n) {
	int count=0;
	int j;
	double thisy=y[0];
	count += (thisy < y[1] && thisy < x[0] && thisy < x[1] && thisy < z[0] && thisy < z[1]);
	for (j=1; j < (n-1); j++) {
	    thisy=y[j];
	    count += (thisy < y[j-1] && thisy < y[j+1] && thisy < x[j-1] && thisy < x[j] && thisy < x[j+1] && thisy < z[j-1] && thisy < z[j] && thisy < z[j+1]);
	}
	thisy=y[n-1];
	count+= (thisy < y[n-2] && thisy < x[n-1] && thisy < x[n-2] && thisy < z[n-2] && thisy < z[n-1]);
	return count;
}


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray*prhs[] )
{
    double *outp;
    double *x,*y, *z;
    unsigned int m, xm, ym, zm, newDeg;

    /* Check for proper number of arguments */

    if (! (nrhs == 4)) {
	mexErrMsgTxt("4 input arguments required.");
    } else if (nlhs > 1) {
	mexErrMsgTxt("Too many output arguments.");
    }

    /* Check the dimensions of Y,X,Z */

    ym = MAX(mxGetM(Y_IN), mxGetN(Y_IN));
    zm = MAX(mxGetM(Z_IN), mxGetN(Z_IN));
    xm = MAX(mxGetM(X_IN), mxGetN(X_IN));

    m = MAX(ym, MAX(xm, zm));

    if(!m) {
      mexErrMsgTxt("Arrays must not be empty.");
    }

    if (!mxIsDouble(Y_IN) || mxIsComplex(Y_IN) || ym != m || zm != m || xm != m) {
         mexErrMsgTxt("All 3 arrays must have the same length, and be real valued");
    }


	/* Set up the input arguments for the subroutine */
    x = mxGetPr(X_IN);
    y = mxGetPr(Y_IN);
    z = mxGetPr(Z_IN);

	/* Here we cheat to learn the exact size of the output */
	/* basically, we trade cycles for memory allocation*/
	newDeg=count2dmin(x,y,z,m);

    /* Create a matrix for the return argument */
	OUT = mxCreateDoubleMatrix(newDeg, 1, mxREAL);

    /* Assign pointers to the various parameters */
    outp = mxGetPr(OUT);

    /* Do the actual computations in a subroutine */
    is2dmin(x, y, z, outp, m);

    return;
}
