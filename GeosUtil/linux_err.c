/* $Id: linux_err.c,v 1.1 2009/11/20 21:43:04 bmy Exp $ */

#include <math.h>

int is_nan_( double *x ) {

  /*====================================================================
   * C function "is_nan_" is a wrapper for the Linux function isnan(x).  
   * (bmy, 3/22/02)
   * 
   * isnan(x) returns 
   *     non-zero : if x is Not-a-Number (NaN),
   *     zero     : otherwise
   *
   * Note that we must declare x as a pointer (i.e. we write *x
   * instead of x) since FORTRAN will pass by reference.  In other
   * words, FORTRAN will pass the memory location of the variable x
   * to this routine.
   *
   * Also, the underscore in the last character of the routine name
   * "is_nan_" is needed so that we can call this from FORTRAN. 
   *==================================================================== 
   */   

  /* isnan is a library call, return value to calling program */ 
  return isnan( *x );
  
}


int is_inf_( double *x ) {

  /*====================================================================
   * C function "is_inf_" is a wrapper for the Linux function isinf(x).  
   * (bmy, 3/22/02)
   * 
   * isinf(x) returns 
   *     -1 : if x is negative infinity
   *      1 : if x is positive infinity
   *      0 : otherwise
   *
   * Note that we must declare x as a pointer (i.e. we write *x
   * instead of x) since FORTRAN will pass by reference.  In other
   * words, FORTRAN will pass the memory location of the variable x
   * to this routine.
   *
   * Also, the underscore in the last character of the routine name
   * "is_inf_" is needed so that we can call this from FORTRAN. 
   *====================================================================
   */
   
  /* isinf is a library call, return value to calling program */
  return isinf( *x );
  
}


int is_finite_( double *x ) {

  /*====================================================================
   * C function "is_finite_" is a wrapper for the Linux function 
   * finite(x). (bmy, 3/22/02)
   * 
   * finite(x) returns 
   *     non-zero : if x is +/- infinity or NaNis negative infinity
   *     zero     : otherwise
   *
   * Note that we must declare x as a pointer (i.e. we write *x
   * instead of x) since FORTRAN will pass by reference.  In other
   * words, FORTRAN will pass the memory location of the variable x
   * to this routine.
   *
   * Also, the underscore in the last character of the routine name
   * "is_finite_" is needed so that we can call this from FORTRAN. 
   *====================================================================
   */

  /* finite is a library call, return value to calling program */
  return finite( *x );
  
}
