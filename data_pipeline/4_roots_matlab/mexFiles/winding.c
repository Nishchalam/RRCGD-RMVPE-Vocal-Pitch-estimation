/* Mex file and C code to compute the winding number of a polygonal curve about 0.

Author:   James W. Fox

   Find the winding number of a polygon about 0.

   Sometimes when doing the winding number, we wish to automatically close
   the polygon. Sometimes we do not. For example when computing the number of
   roots of a real polynomial inside a circle of radius r, we use the FFT that
   is optimized for real sequences to evaluate the polygon only in the upper 
   half-plane. In that case we do not want to close the polygon.

   I found this algorithm in a post by Ken McElvain decwrl.sci.kenm at the 
   bottom of http://info.acm.org/tog/resources/RTNews/html/rtnv3n4.html#art22
   I modified and slightly optimized  his code which computes the winding
   number about an arbitrary point; this code works only for the winding number
   about 0.
*/

typedef struct DCOMPLEX {double r,i;} dcomplex;

double winding(dcomplex *z, long n, long close)
{ /* close = 1 means close the polygon by repeating the last value at the beginning.
     close = 0 means do not close. This is useful if you want to use a real FFT
             to only evaluate at upper half-plane unit-circle roots.

     The quadrants are numbered:      1  |  0
                                      -------
                                      2  |  3
     For definiteness, suppose a point is in quadrant 0. If the next point is in:
       quadrant 0, do nothing.
       quadrant 1, quadcount += 1
       quadrant 3, quadcount -= 1
       quadrant 2, quadcount += 2 or quadcount -= 2 depending on 
          whether the shortest distance to the new point
          goes in a positive or negative direction. 
     When you are done with all the points, quadcount /= 4.

     If the next point is in quadrant 2, how do you tell if it is +=2 or -=2? 
     Draw a straight line to the next point. If the y-intercept is >0 the arc is 
     traversed in a positive direction. If the y-intercept is <0, the arc is negative.

     If the first point had been in quadrant 2 and the next point in quadrant 0,
     then reverse the sense: y-intercept>0 means the arc is negative. In summary,
     the angle is positive if and only if either:
	  point is in quad 0 & the next point is in quad 2 and y-intercept >0 or
	  point is in quad 2 & the next point is in quad 0 and y-intercept <0.

     Suppose the two points are (x1,y1) and (x2,y2).  y1 = m*x1 + b. Thus
     b = y1 - m*x1 = y1 - ((y2-y1)/(x2-x1))x1 = (y1(x2-x1)-x1(y2-y1))/(x2-x1)
     Let numer= y1(x2-x1)-x1(y2-y1) so b = numer/(x2-x1). Then
       if x2-x1>0 then b<0 if and only if numer < 0.
       if x2-x1<0 then b>0 if and only if numer < 0.
     Thus in either case, the angle is positive if and only if numer < 0.
     Similarly, the angle is negative if and only if numer > 0. 
  */
  short oldquad, newquad ;
  long j, start, wind = 0 ;
  dcomplex oldpt, newpt ;
  double numer ;

  if(close) oldpt = z[n-1] ; /* Repeat the last value at the beginning. */
  else oldpt = z[0] ;

  if(oldpt.r >= 0)
    { if(oldpt.i > 0) 
        oldquad = 0;
      else 
        oldquad = 3;
    }
  else
    { if(oldpt.i >= 0) 
         oldquad = 1;
      else 
         oldquad = 2 ;
    }

  if(close) start = 0 ;
  else start = 1 ;

  for(j=start; j<n; ++j)
    { newpt = z[j] ;

      if(newpt.r >= 0)
        { if(newpt.i > 0) newquad = 0;
          else newquad = 3;
        }
      else
        { if(newpt.i >= 0) newquad = 1;
          else newquad = 2 ;
        }

      if(oldquad != newquad)
        { /* If newquad = oldquad+1 mod 4, then ++wind.
             If oldquad = newquad+1 mod 4, then --wind.
             Note, mod(m,4) == m&3 but m&3 is faster.
          */
          if(newquad == ((oldquad+1) & 3)) 
             ++wind ; 
          else if(oldquad == ((newquad+1) & 3))
             --wind ; 
          else
            { numer = oldpt.i*(newpt.r-oldpt.r) - oldpt.r*(newpt.i-oldpt.i) ; 
		  if(numer < 0)
		    { wind += 2 ; }
		  else                  /* What should I do if numer==0 ??? */
		    { wind -= 2 ; }
            }
        }
	oldquad = newquad ;
	oldpt = newpt ;
    }

  if(close)
    { return (double)wind/4 ; }
  else
    { return (double)wind/2 ; }
}
/*======================================*/
/*======================================*/
#include "mex.h"

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#define Y (plhs[0])

#define Z (prhs[0])
#define F (prhs[1])
/*-------------------------------------------------------*/
void mexFunction(int nlhs,       mxArray *plhs[], 
	           int nrhs, const mxArray *prhs[])
{ unsigned int mZ, nZ;
  long nPts, counter;
  long Flag;
  double *ZPr, *ZPi, *YPr, *wPr;
  mxArray *w;  

  /* Check for proper number of arguments */
  if (nlhs > 1)
    { mexErrMsgTxt("Only one output argument is allowed."); }
  else if(nrhs != 2)
    { mexErrMsgTxt("Supply exactly two input arguments."); }

  /* Find the number of points */
  mZ = mxGetM(Z); 
  nZ = mxGetN(Z);
  nPts = MAX(mZ, nZ);

  ZPr = mxGetPr(Z);
  ZPi = mxGetPi(Z);

  if(ZPi == NULL)
    { mexErrMsgTxt("$1 must have some non-zero imaginary parts."); }
  else if (!mxIsDouble(Z) || ((mZ > 1) && (nZ > 1)))
    { mexErrMsgTxt("$1 must either be empty or a 1xn or nx1 complex matrix."); }

  Flag = (long)mxGetScalar(F);

  if(Flag != 0 && Flag != 1)
    { mexErrMsgTxt("$2 must equal 0 or 1."); }

  /* Copy Z to a true complex array with intertwined real and imaginary parts. */
  w = mxCreateDoubleMatrix(1, 2*nPts, mxREAL);
  wPr = mxGetPr(w);

  counter = nPts + 1;
  while(--counter)
    { *wPr++ = *ZPr++;
    *wPr++ = *ZPi++;   /*   *wPr++ = ZPi ? *ZPi++ : 0;  */
    }
  wPr -= 2*nPts;

  /*  Create the LHS. */
  Y = mxCreateDoubleMatrix(1, 1, mxREAL);
  YPr = mxGetPr(Y);

  *YPr = winding((dcomplex *)wPr, nPts, Flag);
  mxDestroyArray(w);
 }
