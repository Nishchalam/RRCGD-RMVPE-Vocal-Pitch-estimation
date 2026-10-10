/* Mex file to call rdft

 Author: James W. Fox

 Round n up to the next largest power of 2.
 If n is a power of 2, n is returned.
 If n==0, 1 is returned.

 Matlab has a routine by the same name
 but it returns the power of 2 exponent; this
 returns the power of 2.
*/
 
static unsigned long nextpow2(unsigned long n)
{ 
  unsigned long N, log2n=0, ncopy=n;

  if(!n) /* Special case, n==0 */
    { return 1; }

  while(ncopy = (ncopy>>1)) /* Remove one bit at a time until it is 0. */
    { ++log2n ; }
  N = 1<<log2n ;    /* N = 2^log2n = 2^floor(log2(n)) */
  if(N != n)        /* n was not a power of 2. */
    { N = (N<<1) ;}
  return N;
}
/*======================================*/
#include "mex.h"
#include <string.h> /* memcpy() */
#include <math.h>   /* floor() */

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#define F   (prhs[0])
#define DIR (prhs[1])
#define LEN (prhs[2])

#define Y   (plhs[0])
/*-------------------------------------------------------*/
void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{ void rdft(int, int, double *);
  mxArray *work;
  unsigned long mF, nF, len, n, nPadded, k;
  double multiplier, *fr, *fi, *yr, *yi, *workP, *workPsave;
  int isComplex, direction, N; /* direction & N are parms to */
                               /* rdft whose prototype wants int */

  /* Check for proper number of arguments */
  if (nlhs > 1)
    { mexErrMsgTxt("Only one output argument is allowed."); }
  else if (nrhs < 2 || nrhs > 3)
    { mexErrMsgTxt("Supply 2 or 3 iput arguments."); }

  isComplex = mxIsComplex(F);
  fr = mxGetPr(F);
  fi = mxGetPi(F);

  /* Find the polynomial dimensions and length */
  mF = mxGetM(F); 
  nF = mxGetN(F);
  n = MAX(mF, nF);
  
  if(!n)
    { Y = mxCreateDoubleMatrix(0, 0, mxREAL);
      return;
    }

  if (!mxIsDouble(F) || (mF > 1) && (nF > 1))
    { mexErrMsgTxt("$1 must be a 1xn or nx1 non-sparse numeric matrix."); }
  else if (!mxIsDouble(DIR) || mxIsComplex(DIR) || (mxGetM(DIR)!= 1) || (mxGetN(DIR)!= 1))
    { mexErrMsgTxt("$2 should equal 1 or -1."); }
  
  direction = (int)floor(.5 + mxGetScalar(DIR)); /* Round to nearest integer. */

  if (direction != -1 && direction != 1)
    { mexErrMsgTxt("$2 must equal 1 or -1."); }
  else if(direction==-1 && (n-1) != nextpow2(n-1))
    { mexErrMsgTxt("Must have length($1)=1+(a power of 2) for inverse transform."); }
  else if (direction == 1 && isComplex)
    { mexErrMsgTxt("$1 must be real for a forward transform."); }

  if(nrhs == 3)
    { if(!mxIsDouble(LEN) || mxIsComplex(LEN) || (mxGetM(LEN)!= 1) || (mxGetN(LEN)!= 1) || mxGetScalar(LEN)<n)
        { mexErrMsgTxt("$3 must be a positive scalar at least as large as length($1)."); }
      else if(direction == -1)
        { mexErrMsgTxt("$3 is not allowed with inverse transforms."); }
      N = nextpow2(MAX(n, mxGetScalar(LEN)));
    }
  else
    { N = nextpow2(n);
    }
  
  if(mF > 1)
    { if(direction == 1)
        Y = mxCreateDoubleMatrix(N/2+1, 1, mxCOMPLEX);
      else
        Y = mxCreateDoubleMatrix(N, 1, mxREAL);
    }
  else
    { if(direction == 1)
        Y = mxCreateDoubleMatrix(1, N/2+1, mxCOMPLEX);
      else
        Y = mxCreateDoubleMatrix(1, N, mxREAL);
    }
  
  yr = mxGetPr(Y);
  yi = mxGetPi(Y);
  if(N == 1)  /* rdft() dies if N < 2. */
    { *yr = *fr;
      return;
    }
  
  work = mxCreateDoubleMatrix(1, N, mxREAL);
  workP = workPsave = mxGetPr(work); /* workPsave never changes. */
  
  if(direction == 1)
   { memcpy(workPsave, fr ,(size_t)n*sizeof(double)); }
  else
   { *workP++ = *fr;
     *workP++ = *(fr + n - 1);
     ++fr;
     if(isComplex)
       { ++fi; } 
     for(k=1; k<=n-2; ++k)
       { *workP++ = *fr++;
         if(isComplex)
           { *workP++ = *fi++; } /* "-*fi" if you want to match Matlab's FFT */
         else
           { *workP++ = 0.0; }
       }
     workP = workPsave;
   }

/*????????????????????????????????
  mexPrintf("\n Before rdft\n");
  for(k=1; k<=N; ++k)
   { mexPrintf(" %lf",*workP++); }
  workP = workPsave;
????????????????????????????????*/

  rdft(N, direction, workPsave);

/*????????????????????????????????
  mexPrintf("\n After rdft\n");
  for(k=1; k<=N; ++k)
   { mexPrintf(" %lf",*workP++); }
  workP = workPsave;
????????????????????????????????*/
  
   if(direction == 1)
     { *yr++ = *workP; /* The first 2 samples of the output of rdft are the */
     *yi++ = 0.0;      /* real parts of 0 Hz. and Nyquist respectively.  */
     workP += 2;       /* Skip past the second sample, the Nyquist info. */
      for(k=1; k <= N/2-1; ++k)
        { *yr++ = *workP++;
	*yi++ = *workP++; /* "-*workP" if you want to match Matlab's FFT */
        }
      *yr = *(workPsave + 1);
    }
  else
    { multiplier = (double)2/N;
      for(k=1; k<=N; ++k)
        { *workP++ *= multiplier; }
      memcpy(yr, workPsave, (size_t)N*sizeof(double));
    }
  mxDestroyArray(work);
}
