/* Mex file and C code to deflate of 1 root and its conjugate from a real polynomial.
   Author: James W. Fox

   Deflate a root from a real coefficient polynomial. 
   If the root is complex, also deflate its complex conjugate.

   The original polynomial has degree+1 coefficients. The caller must
   ensure *quoP is a real array of length degree if the root is real,
   or length degree-1 otherwise.

   Here is my code to do polynomial division when both numerator
   and denominator are real. Originally deflate.m used deconvR.dll
   and it worked to degree 70,000. However, this can be optimized. 

   If I am deflating a real root, use Horner's algorithm.
   If I am deflating a complex root unroll the inner for() loop. It will
   be 2 terms long. Also get rid of MIN(LDenom-1 ,LQuo-i) by
   splitting the loop into 2 parts. In one part we always go the full
   amount. In the other part we handle the slide out zone which is only
   one term long. Also, for forward deflation, *denomP==1 so we can get
   rid of the line *qp /= *denomP.

   static void deconvR(double *quoP,   unsigned int LQuo, 
                       double *numerP, unsigned int LNumer,
                       double *denomP, unsigned int LDenom)
   {
     double *qp = quoP ;
     unsigned int   i ,j ,k ;

     memcpy((char *)quoP ,(char *)numerP ,(size_t)(MIN(LQuo,LNumer)*sizeof(double))) ;

     for(i=1 ; i <= LQuo ; i++ , qp++)
       { *qp /= *denomP ;
         if(i < LQuo)
           for(k=i , j=1 ; j <= MIN(LDenom-1 ,LQuo-i) ; j++ , k++)
             quoP[k] -= (*qp) * denomP[j] ;
       }
   }
*/

#include <stddef.h>  /* size_t */
#include <string.h>  /* memcpy -- modern compilers warn/error on implicit declaration without this */
#include <math.h>    /* fabs   -- same */

static void deflate(double *quoP,
                    double *numerP, unsigned long degree,
                    double x, double y)   /* (x,y) is the root */
{
  double *qp, *qFence, mul, c2, c3;
  unsigned long LQuo;

  /* BTW, quoP must not overlap the numerP array or memcpy might fail. */

  if(y)
    LQuo = degree - 1; /* numerP has degree+1 coefficients. */
  else
    LQuo = degree;

  /* The following is required. Consider the line "qFence = quoP + LQuo - 2;"
     when degree=2 and LQuo=1. Then qFence=quoP-1. The while loop is skipped
     but the following 2 lines modify *(quoP+1) which is not part of quo. 
  */
  if(LQuo == 1)
    { *quoP = *numerP;
      return;
    }

  if(y) /* Deflate a complex root and its conjugate; i.e., deflate a quadratic. */
    { /* (z-(x+iy)(*(z-(x-iy)) = z^2 -2xz + (x*x+y*y) = [c1,c2,c3] = [1,c2,c3] */
      c2 = -2.0*x;
      c3 = x*x + y*y;

      if(c3 <= 1) /* Use forward deflation for roots inside the unit circle.*/
        { memcpy((char *)quoP ,(char *)numerP ,(size_t)LQuo*sizeof(double)) ;
          qp = quoP;
          qFence = quoP + LQuo - 3;
          while(qp <= qFence)
            { /* *qp /= c1;         However, c1==1.0 so we can skip this. */
              mul = *qp++;
              *qp++ -= mul*c2;
              *qp-- -= mul*c3;
            }
          mul = *qp++;
          *qp -= mul*c2;
        }
      else /* Use backward deflation. */
        { memcpy((char *)quoP ,(char *)(numerP+2) ,(size_t)LQuo*sizeof(double)) ;
	qp = quoP + LQuo - 1; /* The last element of quo */
	qFence = quoP + 2;    /* The third element. */
          while(qp >= qFence)
            { *qp /= c3;        /* Indexing backwards, c3 is the first coefficient. */
              mul = *qp--;
              *qp-- -= mul*c2;
              *qp++ -= mul;     /* The coefficient is 1.0; no need to multiply by that. */
            }
          *qp /= c3;
          mul = *qp--;
          *qp -= mul*c2;
          *qp /= c3;
        }
    }
  else if(x) /* It is a non-zero real root. */
    { if(fabs(x) <= 1) /* Use forward polynomial evaluation - Horner's method */
        { memcpy((char *)quoP ,(char *)numerP ,(size_t)LQuo*sizeof(double)) ;
          qp = quoP;
          qFence = quoP + LQuo - 2;
          while(qp <= qFence)
            { mul = *qp++;
              *qp += x*mul;              
            }
        }
    else /* Use backward polynomial evaluation on 1/x. */
        { x = 1/x;
          memcpy((char *)quoP ,(char *)(numerP+1) ,(size_t)LQuo*sizeof(double)) ;
          qp = quoP + LQuo - 1; /* The last element of quo */
          qFence = quoP + 1;
          while(qp >= qFence)
            { mul = *qp--;
              *qp += x*mul;
            }
          /* Normalize */
          x = -x;
          qp = quoP;
          qFence = quoP + LQuo - 1; /* The last sample. */
          while(qp <= qFence)
            { *qp++ *= x; }
        }
    }
  else /* Deflate 0 */
    { memcpy((char *)quoP ,(char *)numerP ,(size_t)LQuo*sizeof(double)) ;
    }
}
/*-------------------------------------------------------*/
#include "mex.h"

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#define NUMER prhs[0]
#define ROOT  prhs[1]
#define QUO   plhs[0]

void mexFunction(int nlhs,       mxArray *plhs[], 
         	     int nrhs, const mxArray *prhs[])
{ 
  unsigned long mNumer, nNumer, degree, LQuo;
  double x, y;

  /* Check for proper number of arguments */
  if (nlhs > 1)
    { mexErrMsgTxt("Only one output argument is allowed."); }
  else if(nrhs != 2)
    { mexErrMsgTxt("Supply exactly two input arguments."); }
    
  /* Check the dimensions */
  mNumer = mxGetM(NUMER); 
  nNumer = mxGetN(NUMER);
  degree = MAX(mNumer, nNumer) - 1;
 
  /* Find the root, (x,y) */
  x = (mxGetPr(ROOT) != NULL) ? *mxGetPr(ROOT) : 0;
  y = (mxGetPi(ROOT) != NULL) ? *mxGetPi(ROOT) : 0;

  LQuo = y ? (degree - 1) : degree;

  if(!mxIsDouble(NUMER) || mxIsComplex(NUMER) || mNumer!=1 && nNumer!=1) 
    { mexErrMsgTxt("$1 must be a 1xn or nx1 real matrix."); }
  else if(!mxIsDouble(ROOT) || mxGetM(ROOT)!=1 || mxGetN(ROOT)!=1) 
    { mexErrMsgTxt("$2 must be a scalar."); } 
  else if(degree < 1)
    { mexErrMsgTxt("$1 must have at least 2 values."); }
  if(degree == 1 && y)
    { mexErrMsgTxt("$1 must have at least 3 values to deflate a non-real root."); } 

  /* Create a matrix for the return argument */
  if(nNumer > 1) /* It is a row matrix. */
    { QUO = mxCreateDoubleMatrix(1, LQuo, mxREAL); }
  else
    { QUO = mxCreateDoubleMatrix(LQuo, 1, mxREAL); }

  deflate(mxGetPr(QUO), mxGetPr(NUMER),degree, x,y);
}
/*-------------------------------------------------------*/
