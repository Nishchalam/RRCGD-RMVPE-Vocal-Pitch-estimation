/* Copyright (c) James W. Fox 2000. I release this as copyrighted
   freeware. You may make any commercial use of this code without 
   fee or royalty. This is NOT released under the Gnu Public License;
   i.e. you can use this code without publicly releasing the source to 
   your private code. You may make any changes you wish to bitReverse().
   However, if you make changes, document your changes but retain 
   this copyright notice unchanged. End of copyright notice.
*/

/* Of course, you can delete the accompanying driver which tests bitReverse(). */

/*
   PURPOSE
     bitReverse(m,n) computes the bit reversal of 0,1,2,3,...,n-1.
     This is useful for constructing fast Fourier transforms
     and for quasi-monte carlo methods. Cf the van der Corput sequence
   INPUTS
     long *m - the array to be filled in. You must allocate the space.
     long n  - the length of m. It must be a power of 2.
   OUTPUT
     If n is not a power of 2, NULL is returned.
     Otherwise m is returned.
*/

#include <stdio.h>   /* NULL */
unsigned long *bitReverse(unsigned long *m, unsigned long n)
{ /* I used unsigned long because all integers must be non-negative. */
  unsigned long *m1, *m2, *msub1, j, k, log2n, ncopy=n;
  unsigned long decr, incr=2; /* increasing & decreasing powers of 2 */

  /* Check that n is a power of 2 and compute log2n = log2(n). */
  log2n = 0 ;
  while(n = (n>>1)) /* Remove one bit at a time until it is 0. */
    { ++log2n ; }
  n = 1<<log2n ;    /* n = 2^log2n */
  if(n != ncopy)
    { return (unsigned long *)NULL ; } /* n was not a power of 2. */

  m1 = m ;
  msub1 = m + 1 ; /* &m[1]; contains the first non-zero value. */
  /* m and msub1 never change. m1 starts at &m[0] and marches
     straight to the end of m, setting all values along the way. 
     At the beginning of each iteration m2 is reset to &m[1].

     Set m[0:3] = <0, n/2, n/4, n/2 + n/4> 
  */
  *m1++ = 0 ;               if(n==1) return m ;
  *m1++ = n >> 1 ;          if(n==2) return m ;
  *m1++ = decr = (n >> 2) ;
  *m1++ = m[1] + m[2] ;     if(n==4) return m ;

  j = log2n - 1 ;
  while(--j)
    { m2 = msub1 ; /* Prepare for a new run starting at the beginning. */
      *m1++ = (decr >>= 1);   /* *m1 = decr = decr/2 */
      k = (incr <<= 1);       /* k = incr = incr*2   */
	while(--k)  /* Each time is twice as long as the previous time. */
	  { *m1++ = decr | *m2++ ; }
    }
  return m ;
}
/*************************************************/
/* Driver to test bitReverse(m,n)
   It prompts for n, fills in m, and prints the values.
*/

#include <stdlib.h>  /* malloc() */
#include <stdio.h>   /* printf() */

void main()
{ unsigned long *m, *m1, j, k, n;
  
  printf("\n\n\tProgram to test bitReverse(n) for n a power of 2\n") ;
  printf("\t--------------------------------------------\n") ;
  while(1)
    { printf("\nEnter n (or enter 0 to exit). => ") ;
      fscanf(stdin, "%uld", &n) ; /* Aborts if you enter a negative value. */
      if(!n) goto endit ;

      m = malloc(n*sizeof (long));
      m1 = bitReverse(m, n) ;

      /* Print the values generated. */
      if(m1 == NULL)
        { printf("\n\t<This was not a power of 2.>\n") ; 
          free(m) ;
        }
      else
        { for(j = 1; j <= n/10; ++j) /* n/10 is the number of rows. */
	      { for(k=1; k<=10; ++k)   /* Print 10 numbers per row.   */
                { printf("%7d", *m++) ; } /* Hopefully, *m++ < 1E6. */
              printf("\n") ;
            }
          for(j = 1; j <= n%10; ++j)  /* n%10 values in the last partial row. */
            { printf("%7d", *m++) ; }
          printf("\n") ;
          free(m-n) ; /* m-n was the original value of m before all the *m++ */
        }
    }
endit: ;
}