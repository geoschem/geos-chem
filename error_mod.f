! $Id: error_mod.f,v 1.16 2008/01/31 15:41:58 bmy Exp $
      MODULE ERROR_MOD
!
!******************************************************************************
!  Module ERROR_MOD contains error checking routines. (bmy, 3/8/01, 8/14/07)
!
!  Module Routines:
!  ===========================================================================
!  (1 ) NAN_FLOAT        : Checks REAL*4 numbers for IEEE NaN flag
!  (2 ) NAN_DBLE         : Checks REAL*8 numbers for IEEE NaN flag
!  (3 ) FINITE_FLOAT     : Checks REAL*4 numbers for IEEE Infinity flag
!  (4 ) FINITE_DBLE      : Checks REAL*8 numbers for IEEE Infinity flag
!  (5 ) CHECK_REAL_VALUE : Convenience routine for REAL*4 error checking
!  (6 ) CHECK_DBLE_VALUE : Convenience routine for REAL*8 error checking
!  (7 ) ERROR_STOP       : Wrapper for GEOS_CHEM_STOP; also prints error msg
!  (8 ) GEOS_CHEM_STOP   : Deallocates all module arrays and stops the run
!  (9 ) ALLOC_ERR        : Prints error msg for memory allocation errors
!  (10) DEBUG_MSG        : Prints a debug message and flushes stdout buffer
!
!  Module Interfaces:
!  ===========================================================================
!  (1 ) IT_IS_NAN        : Overloads NAN_FLOAT,        NAN_DBLE
!  (2 ) IT_IS_FINITE     : Overloads FINITE_FLOAT,     FINITE_DBLE
!  (3 ) CHECK_VALUE      : Overloads CHECK_REAL_VALUE, CHECK_DBLE_VALUE
!
!  GEOS-CHEM modules referenced by error_mod.f
!  ============================================================================
!  none
!
!  NOTES:
!  (1 ) Added subroutines CHECK_REAL_VALUE and CHECK_DBLE_VALUE, which are
!        overloaded by interface CHECK_VALUE.  This is a convenience
!        so that you don't have to always call IT_IS_NAN directly.
!        (bmy, 6/13/01)
!  (2 ) Updated comments (bmy, 9/4/01)
!  (3 ) Now use correct values for bit masking in FINITE_FLOAT for the
!        ALPHA platform (bmy, 11/15/01)
!  (4 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Also add MODULE INTERFACES section,
!        since we have an interface here. (bmy, 5/28/02)
!  (5 ) Add NaN and infinity error checking for Linux platform (bmy, 3/22/02)
!  (6 ) Added routines ERROR_STOP, GEOS_CHEM_STOP, and ALLOC_ERR to this
!        module.  Also improved CHECK_STT. (bmy, 11/27/02)
!  (7 ) Minor bug fixes in FORMAT statements.   Renamed cpp switch from 
!        DEC_COMPAQ to COMPAQ.  Also added code to trap errors on SUN 
!        platform. (bmy, 3/21/03)
!  (8 ) Added patches for IBM/AIX platform (gcc, bmy, 6/27/03)
!  (9 ) Bug fixes for LINUX platform (bmy, 9/29/03)
!  (10) Now supports INTEL_FC compiler (bmy, 10/24/03)
!  (11) Changed the name of some cpp switches in "define.h" (bmy, 12/2/03)
!  (12) Minor fix for LINUX_IFC and LINUX_EFC (bmy, 1/24/04)
!  (13) Do not flush buffer for LINUX_EFC in ERROR_STOP (bmy, 4/6/04)
!  (14) Move CHECK_STT routine to "tracer_mod.f" (bmy, 7/20/04)
!  (15) Added LINUX_IFORT switch for Intel v8 and v9 compilers (bmy, 10/18/05)
!  (16) Now print IFORT error messages for Intel v8/v9 compiler (bmy, 11/30/05)
!  (17) Cosmetic change in DEBUG_MSG (bmy, 4/10/06)
!  (18) Remove support for LINUX_IFC and LINUX_EFC compilers (bmy, 8/4/06)
!  (19) Now use intrinsic functions for IFORT, remove C routines (bmy, 8/14/07)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "error_mod.f"
      !=================================================================

      ! Make everything PUBLIC ...
      PUBLIC 

      ! ... except these routines
      PRIVATE NAN_FLOAT,        NAN_DBLE
      PRIVATE FINITE_FLOAT,     FINITE_DBLE
      PRIVATE CHECK_REAL_VALUE, CHECK_DBLE_VALUE

      !=================================================================
      ! MODULE INTERFACES -- "bind" two or more routines with different
      ! argument types or # of arguments under one unique name
      !================================================================= 

      ! Interface for NaN-check routines
      INTERFACE IT_IS_NAN
         MODULE PROCEDURE NAN_FLOAT, NAN_DBLE
      END INTERFACE

      ! Interface for finite-check routines
      INTERFACE IT_IS_FINITE
         MODULE PROCEDURE FINITE_FLOAT, FINITE_DBLE
      END INTERFACE

      ! Interface for check-value routines
      INTERFACE CHECK_VALUE
         MODULE PROCEDURE CHECK_REAL_VALUE, CHECK_DBLE_VALUE
      END INTERFACE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      FUNCTION NAN_FLOAT( VALUE ) RESULT( IT_IS_A_NAN )
!
!******************************************************************************
!  Module NAN_FLOAT returns TRUE if a REAL*4 number is equal to the IEEE NaN 
!  (Not-a-Number) flag.  Returns FALSE otherwise. (bmy, 3/8/01, 8/14/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1) VALUE (REAL*4) : Number to be tested for NaN
!
!  NOTES:
!  (1 ) Is overloaded by interface "IT_IS_NAN".
!  (2 ) Now call C routine is_nan(x) for Linux platform (bmy, 6/13/02)
!  (3 ) Eliminate IF statement in Linux section.  Also now trap NaN on
!        the Sun/Sparc platform.  Rename cpp switch from DEC_COMPAQ to
!        COMPAQ. (bmy, 3/23/03)
!  (4 ) Added patches for IBM/AIX platform (gcc, bmy, 6/27/03)
!  (5 ) Use LINUX error-trapping for INTEL_FC (bmy, 10/24/03)
!  (6 ) Renamed SGI to SGI_MIPS, LINUX to LINUX_PGI, INTEL_FC to INTEL_IFC,
!        and added LINUX_EFC. (bmy, 12/2/03)
!  (7 ) Added LINUX_IFORT switch for Intel v8 and v9 compilers (bmy, 10/18/05)
!  (8 ) Remove support for LINUX_IFC & LINUX_EFC compilers (bmy, 8/4/06)
!  (9 ) Now use ISNAN for Linux/IFORT compiler (bmy, 8/14/07)
!******************************************************************************
!
#     include "define.h" ! C-preprocessor switches

#if   defined( IBM_AIX )
      USE IEEE_ARITHMETIC
#endif

      ! Arguments
      REAL*4, INTENT(IN) :: VALUE

      ! Return value
      LOGICAL            :: IT_IS_A_NAN

      !=================================================================
      ! NAN_FLOAT begins here!
      !=================================================================
#if   defined( SGI_MIPS )
      IT_IS_A_NAN = IEEE_IS_NAN( VALUE )   

#elif defined( COMPAQ )
      IT_IS_A_NAN = ISNAN( VALUE )         

!------------------------------------------------------------------------------
! Prior to 8/14/07:
!#elif defined( LINUX_IFORT ) || defined( LINUX_PGI )
!
!      ! Declare IS_NAN as an external function
!      INTEGER, EXTERNAL  :: IS_NAN
!      
!      ! For LINUX or INTEL_FC compilers, use C routine "is_nan" to test if 
!      ! VALUE is NaN.   VALUE must be cast to DBLE since "is_nan" only
!      ! takes doubles.
!      IT_IS_A_NAN = ( IS_NAN( DBLE( VALUE ) ) /= 0 )
!------------------------------------------------------------------------------

#elif defined( LINUX_IFORT )
      IT_IS_A_NAN = ISNAN( VALUE )    

#elif defined( LINUX_PGI )

      ! Declare IS_NAN as an external function
      INTEGER, EXTERNAL  :: IS_NAN
      
      ! For LINUX or INTEL_FC compilers, use C routine "is_nan" to test if 
      ! VALUE is NaN.   VALUE must be cast to DBLE since "is_nan" only
      ! takes doubles.
      IT_IS_A_NAN = ( IS_NAN( DBLE( VALUE ) ) /= 0 )

#elif defined( SPARC )
!-----------------------------------------------------------------------------
! NOTE: If you compile with SunStudio11/12 with the -fast optimization, this 
! will turn on -ftrap=common, which checks for NaN, invalid, division, and 
! inexact IEEE math errors. (bmy, 12/18/07)
!
!      ! Declare IR_ISNAN as an external function
!      INTEGER, EXTERNAL :: IR_ISNAN
!
!      ! Test if VALUE is a NaN
!      IT_IS_A_NAN = ( IR_ISNAN( VALUE ) /= 0 )
!-----------------------------------------------------------------------------
      IT_IS_A_NAN = .FALSE.

#elif defined( IBM_AIX )

      ! For IBM/AIX platform
      IF ( IEEE_SUPPORT_DATATYPE( VALUE ) ) THEN
         IT_IS_A_NAN = IEEE_IS_NAN( VALUE )
      ENDIF

#endif

      ! Return to calling program
      END FUNCTION NAN_FLOAT

!------------------------------------------------------------------------------

      FUNCTION NAN_DBLE( VALUE ) RESULT( IT_IS_A_NAN )
!
!******************************************************************************
!  Module NAN_DBLE returns TRUE if a REAL*8 number is equal to the IEEE NaN 
! (Not-a-Number) flag.  Returns FALSE otherwise. (bmy, 3/8/01, 8/14/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1) VALUE (REAL*8) : Number to be tested for NaN
!
!  NOTES:
!  (1 ) Is overloaded by interface "IT_IS_NAN".
!  (2 ) Now call C routine is_nan(x) for Linux platform (bmy, 6/13/02)
!  (3 ) Eliminate IF statement in Linux section.  Also now trap NaN on
!        the Sun/Sparc platform.  Rename cpp switch from DEC_COMPAQ to
!        COMPAQ. (bmy, 3/23/03)
!  (4 ) Added patches for IBM/AIX (gcc, bmy, 6/27/03)
!  (5 ) Use LINUX error-trapping for INTEL_FC (bmy, 10/24/03)
!  (6 ) Renamed SGI to SGI_MIPS, LINUX to LINUX_PGI, INTEL_FC to INTEL_IFC,
!        and added LINUX_EFC. (bmy, 12/2/03)
!  (7 ) Added LINUX_IFORT switch for Intel v8 and v9 compilers (bmy, 10/18/05)
!  (8 ) Remove support for LINUX_IFC & LINUX_EFC compilers (bmy, 8/4/06)
!  (9 ) Now use ISNAN for Linux/IFORT compiler (bmy, 8/14/07)
!******************************************************************************
!
#     include "define.h" ! C-preprocessor switches

#if   defined( IBM_AIX )
      USE IEEE_ARITHMETIC
#endif

      ! Arguments
      REAL*8, INTENT(IN) :: VALUE

      ! Return value
      LOGICAL            :: IT_IS_A_NAN

      !=================================================================
      ! NAN_DBLE begins here!
      !=================================================================
#if   defined( SGI_MIPS )
      IT_IS_A_NAN = IEEE_IS_NAN( VALUE )   

#elif defined( COMPAQ )
      IT_IS_A_NAN = ISNAN( VALUE )         

!-----------------------------------------------------------------------------
! Prior to 8/14/07:
!#elif defined( LINUX_IFORT ) || defined( LINUX_PGI )
!
!      ! Declare IS_NAN as an external function
!      INTEGER, EXTERNAL  :: IS_NAN
!
!      ! For LINUX or INTEL_FC compilers, use C routine 
!      ! "is_nan" to test if VALUE is NaN.  
!      IT_IS_A_NAN = ( IS_NAN( VALUE ) /= 0 )
!-----------------------------------------------------------------------------

#elif defined( LINUX_IFORT )
      IT_IS_A_NAN = ISNAN( VALUE )      

#elif defined( LINUX_PGI )

      ! Declare IS_NAN as an external function
      INTEGER, EXTERNAL  :: IS_NAN

      ! For LINUX or INTEL_FC compilers, use C routine 
      ! "is_nan" to test if VALUE is NaN.  
      IT_IS_A_NAN = ( IS_NAN( VALUE ) /= 0 )

#elif defined( SPARC )
!-----------------------------------------------------------------------------
! NOTE: If you compile with SunStudio11/12 with the -fast optimization, this 
! will turn on -ftrap=common, which checks for NaN, invalid, division, and 
! inexact IEEE math errors. (bmy, 12/18/07)
!
!      ! Declare ID_ISNAN as an external function
!      INTEGER, EXTERNAL  :: ID_ISNAN
!
!      ! Test if VALUE is NaN
!      IT_IS_A_NAN = ( ID_ISNAN( VALUE ) /= 0 )
!-----------------------------------------------------------------------------
      IT_IS_A_NAN = .FALSE.

#elif defined( IBM_AIX )

       ! For IBM/AIX platform
      IF ( IEEE_SUPPORT_DATATYPE( VALUE ) ) THEN
         IT_IS_A_NAN = IEEE_IS_NAN( VALUE )
      ENDIF

#endif

      ! Return to calling program
      END FUNCTION NAN_DBLE

!-----------------------------------------------------------------------------

      FUNCTION FINITE_FLOAT( VALUE ) RESULT( IT_IS_A_FINITE )
!
!******************************************************************************
!  Module FINITE_FLOAT returns TRUE if a REAL*4 number is equal to the 
!  IEEE Infinity flag.  Returns FALSE otherwise. (bmy, 3/8/01, 8/14/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1) VALUE (REAL*4) : Number to be tested for infinity
!
!  NOTES:
!  (1 ) Is overloaded by interface "IT_IS_FINITE".
!  (2 ) Now use correct values for bit masking (bmy, 11/15/01)
!  (3 ) Eliminate IF statement in Linux section.  Also now trap Infinity on
!        the Sun/Sparc platform.  Rename cpp switch from DEC_COMPAQ to
!        COMPAQ. (bmy, 3/23/03)
!  (4 ) Added patches for IBM/AIX platform (gcc, bmy, 6/27/03)
!  (5 ) Bug fix: now use external C IS_FINITE for PGI/Linux (bmy, 9/29/03)
!  (6 ) Use LINUX error-trapping for INTEL_FC (bmy, 10/24/03)
!  (7 ) Renamed SGI to SGI_MIPS, LINUX to LINUX_PGI, INTEL_FC to INTEL_IFC,
!        and added LINUX_EFC. (bmy, 12/2/03)
!  (8 ) Added LINUX_IFORT switch for Intel v8 and v9 compilers (bmy, 10/18/05)
!  (9 ) Remove support for LINUX_IFC & LINUX_EFC compilers (bmy, 8/4/06)
!  (10) Now use FP_CLASS for IFORT compiler (bmy, 8/14/07)
!******************************************************************************
!
#     include "define.h" ! C-preprocessor switches

#if   defined( IBM_AIX )
      USE IEEE_ARITHMETIC
#endif

      ! Arguments
      REAL*4, INTENT(IN) :: VALUE

      ! Return value
      LOGICAL            :: IT_IS_A_FINITE

      !=================================================================
      ! FINITE_FLOAT begins here!
      !=================================================================  
     
#if   defined( SGI_MIPS )
      IT_IS_A_FINITE = IEEE_FINITE( VALUE )  
 
#elif defined( COMPAQ ) 

      ! Test for finite # using bit masking for DEC/Compaq platform
      ! Now use REAL*4 values (bmy, 11/15/01)
      IF ( VALUE == Z'7F800000' .or. 
     &     VALUE == Z'FF800000' ) THEN
         IT_IS_A_FINITE = .FALSE.
      ELSE
         IT_IS_A_FINITE = .TRUE.
      ENDIF

!-----------------------------------------------------------------------------
! Prior to 8/14/07:
!#elif defined( LINUX_IFORT ) || defined( LINUX_PGI ) 
!
!      ! Declare IS_FINITE as an external function
!      INTEGER, EXTERNAL :: IS_FINITE
!      
!      ! For LINUX or INTEL_FC compilers use C routine "is_finite" to test 
!      ! if VALUE is finite.  VALUE must be cast to DBLE since "is_inf" 
!      ! only takes doubles. 
!      IT_IS_A_FINITE = ( IS_FINITE( DBLE( VALUE ) ) /= 0 )
!-----------------------------------------------------------------------------
      IT_IS_A_NAN = .FALSE.

#elif defined( LINUX_IFORT )

      ! Local variables (parameters copied from "fordef.for")
      INTEGER, PARAMETER :: SNAN=0, QNAN=1, POS_INF=2, NEG_INF=3
      INTEGER            :: FPC

      ! Get the floating point type class for VALUE
      FPC            = FP_CLASS( VALUE )

      ! VALUE is infinite if it is either +Inf or -Inf
      ! Also flag an error if VALUE is a signaling or quiet NaN
      IT_IS_A_FINITE = ( FPC /= POS_INF .and. FPC /= NEG_INF .and. 
     &                   FPC /= SNAN    .and. FPC /= QNAN          )

#elif defined( LINUX_PGI ) 

      ! Declare IS_FINITE as an external function
      INTEGER, EXTERNAL :: IS_FINITE
      
      ! For LINUX or INTEL_FC compilers use C routine "is_finite" to test 
      ! if VALUE is finite.  VALUE must be cast to DBLE since "is_inf" 
      ! only takes doubles. 
      IT_IS_A_FINITE = ( IS_FINITE( DBLE( VALUE ) ) /= 0 )

#elif defined( SPARC )
!-----------------------------------------------------------------------------
! NOTE: If you compile with SunStudio11/12 with the -fast optimization, this 
! will turn on -ftrap=common, which checks for NaN, invalid, division, and 
! inexact IEEE math errors. (bmy, 12/18/07)
!      ! Declare IR_FINITE as an external function
!      INTEGER, EXTERNAL :: IR_FINITE
!
!      ! Test if VALUE is a finite number
!      IT_IS_A_FINITE = ( IR_FINITE( VALUE ) /= 0 )
!-----------------------------------------------------------------------------
      IT_IS_A_FINITE = .TRUE.

#elif defined( IBM_AIX )

      ! For IBM/AIX platform
      IF ( IEEE_SUPPORT_DATATYPE( VALUE ) ) THEN
         IT_IS_A_FINITE = IEEE_IS_FINITE( VALUE )
      ENDIF 
      
#endif
      
      ! Return to calling program
      END FUNCTION FINITE_FLOAT

!-----------------------------------------------------------------------------

      FUNCTION FINITE_DBLE( VALUE ) RESULT( IT_IS_A_FINITE )
!
!*****************************************************************************
!  Module FINITE_DBLE returns TRUE if a REAL*8 number is equal to the 
!  IEEE Infinity flag.  Returns FALSE otherwise. (bmy, 3/8/01, 8/14/07)
!
!  Arguments as Input:
!  ===========================================================================
!  (1) VALUE (REAL*8) : Number to be tested for infinity
!
!  NOTES:
!  (1 ) Is overloaded by interface "IT_IS_FINITE".
!  (2 ) Now call C routine is_nan(x) for Linux platform (bmy, 6/13/02)
!  (3 ) Eliminate IF statement in Linux section.  Also now trap Infinity on
!        the Sun/Sparc platform.  Rename cpp switch from DEC_COMPAQ to
!        COMPAQ. (bmy, 3/23/03)
!  (4 ) Added patches for IBM/AIX platform (gcc, bmy, 6/27/03)
!  (5 ) Bug fix: now use external C IS_FINITE for PGI/Linux (bmy, 9/29/03)
!  (6 ) Use LINUX error-trapping for INTEL_FC (bmy, 10/24/03)
!  (7 ) Renamed SGI to SGI_MIPS, LINUX to LINUX_PGI, INTEL_FC to INTEL_IFC,
!        and added LINUX_EFC. (bmy, 12/2/03)
!  (8 ) Added LINUX_IFORT switch for Intel v8 and v9 compilers (bmy, 10/18/05)
!  (9 ) Remove support for LINUX_IFC & LINUX_EFC compilers (bmy, 8/4/06)
!  (10) Now use FP_CLASS for IFORT compiler (bmy, 8/14/07)
!******************************************************************************
!
#     include "define.h" ! C-preprocessor switches

#if   defined( IBM_AIX )
      USE IEEE_ARITHMETIC
#endif

      ! Arguments
      REAL*8, INTENT(IN) :: VALUE

      ! Return value
      LOGICAL            :: IT_IS_A_FINITE
         
      !=================================================================
      ! FINITE_DBLE begins here!
      !=================================================================
#if   defined( SGI_MIPS )
      IT_IS_A_FINITE = IEEE_FINITE( VALUE )  

#elif defined( COMPAQ ) 

      ! Test for finite # using bit masking for DEC/Compaq F90
      IF ( VALUE == Z'7FF0000000000000'  .or. 
     &     VALUE == Z'FFF0000000000000') THEN
         IT_IS_A_FINITE = .FALSE.
      ELSE
         IT_IS_A_FINITE = .TRUE.
      ENDIF

!-----------------------------------------------------------------------------
! Prior to 8/14/07:
!#elif defined( LINUX_IFORT ) || defined( LINUX_PGI )
!
!      ! Declare IS_FINITE as an external function
!      INTEGER, EXTERNAL :: IS_FINITE
!
!      ! For LINUX or INTEL_FC compilers, use C routine 
!      ! "is_finite" to test if VALUE is infinity
!      IT_IS_A_FINITE = ( IS_FINITE( VALUE ) /= 0 )
!-----------------------------------------------------------------------------

#elif defined( LINUX_IFORT )

      ! Local variables (parameters copied from "fordef.for")
      INTEGER, PARAMETER :: SNAN=0, QNAN=1, POS_INF=2, NEG_INF=3
      INTEGER            :: FPC

      ! Get the floating point type class for VALUE
      FPC            = FP_CLASS( VALUE )

      ! VALUE is infinite if it is either +Inf or -Inf
      ! Also flag an error if VALUE is a signaling or quiet NaN
      IT_IS_A_FINITE = ( FPC /= POS_INF .and. FPC /= NEG_INF .and. 
     &                   FPC /= SNAN    .and. FPC /= QNAN          )

#elif defined( LINUX_PGI )

      ! Declare IS_FINITE as an external function
      INTEGER, EXTERNAL :: IS_FINITE

      ! For LINUX or INTEL_FC compilers, use C routine 
      ! "is_finite" to test if VALUE is infinity
      IT_IS_A_FINITE = ( IS_FINITE( VALUE ) /= 0 )

#elif defined( SPARC )
!-----------------------------------------------------------------------------
! NOTE: If you compile with SunStudio11/12 with the -fast optimization, this 
! will turn on -ftrap=common, which checks for NaN, invalid, division, and 
! inexact IEEE math errors. (bmy, 12/18/07)
! 
!      ! Declare ID_FINITE as an external function
!      INTEGER, EXTERNAL :: ID_FINITE
!
!      ! Test if VALUE is a finite number
!      IT_IS_A_FINITE = ( ID_FINITE( VALUE ) /= 0 )
!-----------------------------------------------------------------------------
      IT_IS_A_FINITE = .TRUE.

#elif defined( IBM_AIX )

      ! For IBM/AIX platform
      IF ( IEEE_SUPPORT_DATATYPE( VALUE ) ) THEN
         IT_IS_A_FINITE = IEEE_IS_FINITE( VALUE )
      ENDIF

#endif
      
      ! Return to calling program
      END FUNCTION FINITE_DBLE
      
!------------------------------------------------------------------------------

      SUBROUTINE CHECK_REAL_VALUE( VALUE, LOCATION, VARNAME, MESSAGE )
!
!******************************************************************************
!  Subroutine CHECK_REAL_VALUE checks to make sure a REAL*4 value is not 
!  NaN or Infinity. This is a wrapper for the interface IT_IS_NAN. 
!  (bmy, 6/13/01)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) VALUE    (REAL*4  ) : Value to be checked
!  (2 ) LOCATION (INTEGER ) : 4-element array w/ box indices (I,J,L,N)
!  (3 ) VARNAME  (CHAR*(*)) : Variable name corresponding to VALUE
!  (4 ) MESSAGE  (CHAR*(*)) : A short descriptive message 
!
!  NOTES:
!  (1 ) Now call GEOS_CHEM_STOP to shutdown safely.  Updated comments,
!        cosmetic changes (bmy, 10/15/02)
!******************************************************************************
!
      ! Arguments
      REAL*4,           INTENT(IN) :: VALUE
      CHARACTER(LEN=*), INTENT(IN) :: VARNAME
      CHARACTER(LEN=*), INTENT(IN) :: MESSAGE
      INTEGER,          INTENT(IN) :: LOCATION(4)

      !=================================================================
      ! CHECK_REAL_VALUE begins here!
      !=================================================================

      ! First check for NaN -- print info & stop run if found
      IF ( IT_IS_NAN( VALUE ) ) THEN 
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, 110   ) TRIM( VARNAME )
         WRITE( 6, 115   ) LOCATION
         WRITE( 6, '(a)' ) TRIM( MESSAGE )
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Next check for infinity -- print info & stop run if found
      IF ( .not. IT_IS_FINITE( VALUE ) ) THEN
         WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
         WRITE( 6, 120       ) TRIM( VARNAME )
         WRITE( 6, 115       ) LOCATION
         WRITE( 6, '(f13.6)' ) VALUE      
         WRITE( 6, '(a)'     ) TRIM ( MESSAGE )
         WRITE( 6, '(a)'     ) REPEAT( '=', 79 )         
         CALL GEOS_CHEM_STOP
      ENDIF

      ! FORMAT statements
 110  FORMAT( 'CHECK_VALUE: ', a, ' is NaN!'        )
 115  FORMAT( 'Grid box (I,J,L,N) : ', 4i4          )
 120  FORMAT( 'CHECK_VALUE: ', a, ' is not finite!' )

      ! Return to calling program
      END SUBROUTINE CHECK_REAL_VALUE

!-----------------------------------------------------------------------------

      SUBROUTINE CHECK_DBLE_VALUE( VALUE, LOCATION, VARNAME, MESSAGE )
!
!******************************************************************************
!  Subroutine CHECK_DBLE_VALUE checks to make sure a REAL*8 value is 
!  not NaN or Infinity.  This is a wrapper for the interface IT_IS_NAN. 
!  (bmy, 6/14/01, 10/15/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) VALUE    (REAL*4  ) : Value to be checked
!  (2 ) LOCATION (INTEGER ) : 4-element array w/ box indices (I,J,L,N)
!  (3 ) VARNAME  (CHAR*(*)) : Variable name corresponding to VALUE
!  (4 ) MESSAGE  (CHAR*(*)) : A short descriptive message 
! 
!  NOTES:
!  (1 ) Now call GEOS_CHEM_STOP to shutdown safely.  Updated comments,
!        cosmetic changes (bmy, 10/15/02)
!******************************************************************************
!
      ! Arguments
      REAL*8,           INTENT(IN) :: VALUE
      CHARACTER(LEN=*), INTENT(IN) :: VARNAME
      CHARACTER(LEN=*), INTENT(IN) :: MESSAGE
      INTEGER,          INTENT(IN) :: LOCATION(4)

      !=================================================================
      ! CHECK_DBLE_VALUE begins here!
      !=================================================================

      ! First check for NaN
      IF ( IT_IS_NAN( VALUE ) )THEN 
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, 110   ) TRIM( VARNAME )
         WRITE( 6, 115   ) LOCATION
         WRITE( 6, '(a)' ) TRIM( MESSAGE )
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Next check for infinity
      IF ( .not. IT_IS_FINITE( VALUE ) ) THEN
         WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
         WRITE( 6, 120       ) TRIM( VARNAME )
         WRITE( 6, 115       ) LOCATION
         WRITE( 6, '(f13.6)' ) VALUE      
         WRITE( 6, '(a)'     ) TRIM ( MESSAGE )
         WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      ! FORMAT statements
 110  FORMAT( 'CHECK_VALUE: ', a, ' is NaN!'        )
 115  FORMAT( 'Grid box (I,J,L,N) : ', 4i4          )
 120  FORMAT( 'CHECK_VALUE: ', a, ' is not finite!' )

      ! Return to calling program
      END SUBROUTINE CHECK_DBLE_VALUE

!------------------------------------------------------------------------------

      SUBROUTINE ERROR_STOP( MESSAGE, LOCATION )
!
!******************************************************************************
!  Subroutine ERROR_STOP is a wrapper for GEOS_CHEM_STOP.  It prints an error
!  message then calls GEOS_CHEM_STOP to free memory and quit. (bmy, 10/15/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) MESSAGE  (CHARACTER) : Error message to be printed
!  (2 ) LOCATION (CHARACTER) : Location where ERROR_STOP is called from
!
!  NOTES:
!******************************************************************************
!

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: MESSAGE, LOCATION

      !=================================================================
      ! ERROR_MSG begins here!
      !=================================================================

!$OMP CRITICAL

      ! Write msg
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, '(a)' ) 'GEOS-CHEM ERROR: ' // TRIM( MESSAGE )
      WRITE( 6, '(a)' ) 'STOP at '          // TRIM( LOCATION )
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

!$OMP END CRITICAL

      ! Deallocate memory and stop the run
      CALL GEOS_CHEM_STOP
      
      ! Return to calling program
      END SUBROUTINE ERROR_STOP

!------------------------------------------------------------------------------

      SUBROUTINE GEOS_CHEM_STOP
!
!******************************************************************************
!  Subroutine GEOS_CHEM_STOP calls CLEANUP to deallocate all module arrays
!  and then stops the run.  (bmy, 10/15/02, 1/24/04) 
!
!  NOTES
!  (1 ) Now EXIT works for LINUX_IFC, LINUX_EFC, so remove #if block.
!        (bmy, 1/24/04)
!******************************************************************************
!
#     include "define.h"

      !=================================================================
      ! GEOS_CHEM_STOP begins here!
      !=================================================================

!$OMP CRITICAL

      ! Deallocate all module arrays
      CALL CLEANUP

      ! Flush all files and stop
      CALL EXIT( 99999 )

!$OMP END CRITICAL

      ! End of program
      END SUBROUTINE GEOS_CHEM_STOP

!------------------------------------------------------------------------------

      SUBROUTINE ALLOC_ERR( ARRAYNAME, AS )
!
!******************************************************************************
!  Subroutine ALLOC_ERR prints an error message if there is not enough
!  memory to allocate a particular allocatable array. (bmy, 6/26/00, 11/30/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ARRAYNAME (CHARACTER) : name of the allocatable array
!
!  NOTES:
!  (1 ) Split off from "ndxx_setup.f" into a separate routine (bmy, 6/26/00)
!  (2 ) Added to "error_mod.f" (bmy, 10/15/02)
!  (3 ) Call special error msg function for IFORT compiler (bmy, 11/30/05)
!******************************************************************************
!
#     include "define.h"

      ! Arguments
      CHARACTER(LEN=*),  INTENT(IN) :: ARRAYNAME
      INTEGER, OPTIONAL, INTENT(IN) :: AS

      ! Local variables
      CHARACTER(LEN=255)            :: ERRMSG

      !=================================================================
      ! ALLOC_ERR begins here!
      !=================================================================

#if   defined( LINUX_IFORT )
     
      !-----------------------
      ! Linux/IFORT compiler 
      !-----------------------
 
      ! More local variables
      CHARACTER(LEN=255) :: IFORT_ERRMSG, MSG

      ! Define error message
      ERRMSG = 'Allocation error in array: ' // TRIM( ARRAYNAME )

      ! If we have passed the allocation status argument ...
      IF ( PRESENT( AS ) ) THEN 

         ! Get IFORT error message
         MSG = IFORT_ERRMSG( AS )

         ! Append IFORT error message 
         ERRMSG = TRIM( ERRMSG ) // ' :: ' // TRIM( MSG ) 

      ENDIF

#else

      !-----------------------
      ! All other compilers
      !-----------------------

      ! Define error message
      ERRMSG = 'Allocation error in array: ' // TRIM( ARRAYNAME )
    
#endif
 
      ! Print error message, deallocate memory, and stop the run
      CALL ERROR_STOP( ERRMSG, 'alloc_err.f' )
      
      ! End of subroutine
      END SUBROUTINE ALLOC_ERR

!------------------------------------------------------------------------------

      SUBROUTINE DEBUG_MSG( MESSAGE )
!
!******************************************************************************
!  Subroutine DEBUG_MSG prints a message to the stdout buffer and flushes.
!  This is useful for determining the exact location where errors occur.
!  (bmy, 1/7/02, 8/4/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) MESSAGE (CHAR*(*)) : Descriptive message string
!
!  NOTES: 
!  (1 ) Now just write the message and flush the buffer (bmy, 7/5/01)
!  (2 ) Renamed from "paftop.f" to "debug_msg.f" (bmy, 1/7/02)
!  (3 ) Bundled into "error_mod.f" (bmy, 11/22/02)
!  (4 ) Now do not FLUSH the buffer for EFC compiler (bmy, 4/6/04)
!  (5 ) Now add a little space for debug output (bmy, 4/10/06)
!  (6 ) Remove support for LINUX_IFC & LINUX_EFC compilers (bmy, 8/4/06)
!******************************************************************************
!
      IMPLICIT NONE

#     include "define.h"

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: MESSAGE

      !=================================================================
      ! DEBUG_MSG begins here!
      !================================================================= 
      WRITE( 6, '(5x,a)' ) MESSAGE

      ! Call FLUSH routine to flush the output buffer
      CALL FLUSH( 6 )

      ! Return to calling program
      END SUBROUTINE DEBUG_MSG

!------------------------------------------------------------------------------

      ! End of module
      END MODULE ERROR_MOD
