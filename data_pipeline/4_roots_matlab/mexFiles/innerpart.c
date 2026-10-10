/* Mex file and C code to compute the inner part of a polynomial evaluation. */
/* Author: James W. Fox */


#include <mex.h>
#include <math.h> /* hypot() */

/*======================================
 [y,err]=innerpart(f,z).
   Given  an array of real polynomial coefficients, *f, and a complex array
   (*zPr,*zPi), find the inner part of the the polynomial evaluations.
   Optionally, estimate the error in each evaluation.
  
   For one specific z, this is the same as:
     polyval(f,z)         if |z|<=1
     polyval(flip(f),1/z) if |z|>1

   Author: James W. Fox
*/

static void innerpart(double *yR, double *yI, double *err,
                      const double *f, const long nCoeff, 
                      double *zPr, double *zPi, unsigned long nPts)
{ double *fp,                  /* Moving pointer to polynomial coefficients. */
         xin, yin,             /* The current evaluation point. */
         xout, yout,           /* (xout,yout) = the answer copied to (*yR++,*yI++). */
         absZSq,               /* |((*zPr),(*zPi))|^2 */
         absZ, u, gamma2, mu,  /* See Higham p. 3, 69, 106 */
         tempR, *fence;

  u = pow(2.0,-53);
  gamma2 = 2*u/(1-2*u);     
  ++nPts;
  while(--nPts)
    { xin = *zPr++; /* See below - in Matlab 5.3 pointers to real parts are never NULL. */
      yin = zPi ? *zPi++ : 0; /* (xin, yin) is the current evaluation point. */
      
      absZ = hypot(xin,yin);
      absZSq = absZ*absZ;
      if(absZSq <= 1)
        { xout = f[0];
          yout = 0;
          fp = (double *) f;              /* Start at the first coefficient. */
          fence = (double *)&f[nCoeff];   /* Stop just before you get here. */
          if(err)    /* Also compute error estimates. */
            { absZ = sqrt(absZSq); /* New code does it above ??????????????? */
              mu = abs(*fp)/2.0;
              while(++fp != fence)
                { tempR = xin*xout - yin*yout;
                  yout  = xin*yout + yin*xout;
                  xout  = tempR + *fp;
                  mu    = absZ*mu + sqrt(xout*xout + yout*yout);
                }
              *err++ = sqrt(2.0)*gamma2*(2*mu - sqrt(xout*xout + yout*yout));
            } 
          else
            { while(++fp != fence)
                { tempR = xin*xout - yin*yout;
                  yout  = xin*yout + yin*xout;
                  xout  = tempR + *fp;
                }
            }
        }
      else
        { xin /= absZSq; /*  Replace (xin,yin) with its complex reciprocal. */
          yin /= -absZSq;
          xout = f[nCoeff-1];
          yout = 0;
          fp =  (double *) &f[nCoeff-1]; /* Start at the last coefficient. */
          fence = (double *) f-1;        /* Stop just before you get here. */
          if(err) /* Also compute error estimates. */
            { absZ = 1.0/sqrt(absZSq);
              mu = abs(*fp)/2.0;
              while(--fp != fence)
                { tempR = xin*xout - yin*yout;
                  yout  = xin*yout + yin*xout;
                  xout  = tempR + *fp;
                  mu    = absZ*mu + sqrt(xout*xout + yout*yout);
                }
              *err++ = sqrt(2.0)*gamma2*(2*mu - sqrt(xout*xout + yout*yout));
            }
          else
            { while(--fp != fence)
                { tempR = xin*xout - yin*yout;
                  yout  = xin*yout + yin*xout;
                  xout  = tempR + *fp;
                }
            }
        }
      *yR++ = xout;
      *yI++ = yout;
    }
}    
/*====================================== */
#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#define Y   (plhs[0])
#define ERR (plhs[1])

#define F (prhs[0])
#define Z (prhs[1])
/*------------------------------------------------------- */
/* By putting in   
      mexPrintf("FPr=%d FPi=%d \n\n",FPr,FPi);
   I discovered something interesting - the pointer to the real part is never NULL.

   Matlab 5.3 stores complex arrays as two separate arrays, the real parts and the
   imaginary parts. If the array is real, the pointer to the imaginary part is NULL.
   However, if the array is purely imaginary, like [i i i], the pointer to the real part
   is not NULL. They create an array of all zeros for the real part. I have no idea why
   this asymmetry or whether this is changed in Matlab 6. If Matlab 6 does something 
   different then this code may not work.
*/
void mexFunction(int nlhs,       mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{ unsigned long mF, nF, mZ, nZ, nCoeff, nPts;
/*  long nCoeff; */
/*  int dataTypes; */
  
  double *YPr, *YPi, *ERRPr, *FPr, *FPi, *ZPr, *ZPi;

  /* Check for proper number of arguments */
  if (nlhs > 2)
    { mexErrMsgTxt("Only one or two output arguments are allowed."); }
  else if (nrhs != 2)
    { mexErrMsgTxt("Supply exactly 2 iput arguments."); }
  
  /* Find the polynomial dimensions and number of coefficients */
  mF = mxGetM(F); 
  nF = mxGetN(F);
  nCoeff = MAX(mF, nF);

  /* Find the dimensions of the array of evaluation points. */
  mZ = mxGetM(Z);
  nZ = mxGetN(Z);
  nPts = MAX(mZ, nZ);

  if (nlhs > 2)
    { mexErrMsgTxt("Only one or two output arguments are allowed."); }
  else if (nrhs != 2)
    { mexErrMsgTxt("Supply exactly 2 iput arguments."); }
  else if (!mxIsDouble(F) || mxIsComplex(F) || (mF != 1) && (nF != 1))
    { mexErrMsgTxt("$1 must be a 1xn or nx1 real matrix."); }
  else if (!mxIsDouble(Z) || (mZ > 1) && (nZ > 1))
    { mexErrMsgTxt("$2 must either be empty or a 1xn or nx1 real or complex matrix."); }
  
  FPr = mxGetPr(F);
  /*  FPi = mxGetPi(F); We checked above that F was not complex */

  ZPr = mxGetPr(Z);
  ZPi = mxGetPi(Z);

  /* Create the array for evaluations.  */
  Y = mxCreateDoubleMatrix(mZ, nZ, mxCOMPLEX);
  YPr = mxGetPr(Y);
  YPi = mxGetPi(Y);

  /* Create the array for error estimates. */
  if (nlhs > 1)
    { ERR = mxCreateDoubleMatrix(mZ, nZ, mxREAL);
      ERRPr = mxGetPr(ERR);
    }
  else
    { ERRPr = NULL;
    }

  innerpart(YPr, YPi, ERRPr, mxGetPr(F), nCoeff, ZPr, ZPi, nPts);
 }
