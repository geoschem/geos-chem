! $Id: error_mod.f,v 1.2 2003/10/01 20:32:21 bmy Exp $
      MODULE ERROR_MOD
!
!******************************************************************************
!  Module ERROR_MOD contains error checking routines. (bmy, 3/8/01, 9/27/03)
!
!  Module Routines:
!  ===========================================================================
!  (1 ) NAN_FLOAT        : Checks REAL*4 numbers for IEEE NaN flag
!  (2 ) NAN_DBLE         : Checks REAL*8 numbers for IEEE NaN flag
!  (3 ) FINITE_FLOAT     : Checks REAL*4 numbers for IEEE Infinity flag
!  (4 ) FINITE_DBLE      : Checks REAL*8 numbers for IEEE Infinity flag
!  (5 ) CHECK_STT        : Checks STT array for Negatives, NaN's, or Infinities
!  (6 ) CHECK_REAL_VALUE : Convenience routine for REAL*4 error checking
!  (7 ) CHECK_DBLE_VALUE : Convenience routine for REAL*8 error checking
!  (8 ) ERROR_STOP       : Wrapper for GEOS_CHEM_STOP; also prints error msg
!  (9 ) GEOS_CHEM_STOP   : Deallocates all module arrays and stops the run
!  (10) ALLOC_ERR        : Prints error msg for memory allocation errors
!  (11) DEBUG_MSG        : Prints a debug message and flushes stdout buffer
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
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "error_mod.f"
      !=================================================================

      ! PRIVATE module routines
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

!-----------------------------------------------------------------------------

      FUNCTION NAN_FLOAT( VALUE ) RESULT( IT_IS_A_NAN )
!
!*****************************************************************************
!  Module NAN_FLOAT returns TRUE if a REAL*4 number is equal to the IEEE NaN 
!  (Not-a-Number) flag.  Returns FALSE otherwise. (bmy, 3/8/01, 3/23/03)
!
!  Arguments as Input:
!  ===========================================================================
!  (1) VALUE (REAL*4) : Number to be tested for NaN
!
!  NOTES:
!  (1 ) Is overloaded by interface "IT_IS_NAN".
!  (2 ) Now call C routine is_nan(x) for Linux platform (bmy, 6/13/02)
!  (3 ) Eliminate IF statement in Linux section.  Also now trap NaN on
!        the Sun/Sparc platform.  Rename cpp switch from DEC_COMPAQ to
!        COMPAQ. (bmy, 3/23/03)
!  (4 ) Added patches for IBM/AIX platform (gcc, bmy, 6/27/03)
!*****************************************************************************
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
#if   defined( SGI )
      IT_IS_A_NAN = IEEE_IS_NAN( VALUE )   ! for SGI platform

#elif defined( COMPAQ )
      IT_IS_A_NAN = ISNAN( VALUE )         ! for DEC, Compaq platforms

#elif defined( LINUX )

      ! Declare IS_NAN as an external function
      INTEGER, EXTERNAL  :: IS_NAN

      ! For LINUX platform, use C routine "is_nan" to test if VALUE is NaN.  
      ! VALUE must be cast to DBLE since "is_nan" only takes doubles.
      IT_IS_A_NAN = ( IS_NAN( DBLE( VALUE ) ) /= 0 )

#elif defined( SPARC )

      ! Declare IR_ISNAN as an external function
      INTEGER, EXTERNAL :: IR_ISNAN

      ! Test if VALUE is a NaN
      IT_IS_A_NAN = ( IR_ISNAN( VALUE ) /= 0 )

#elif defined( IBM_AIX )

       ! For IBM/AIX platform
       IF ( IEEE_SUPPORT_DATATYPE( VALUE ) ) THEN
          IT_IS_A_NAN = IEEE_IS_NAN( VALUE )
       ENDIF

#endif

      ! Return to calling program
      END FUNCTION NAN_FLOAT

!-----------------------------------------------------------------------------

      FUNCTION NAN_DBLE( VALUE ) RESULT( IT_IS_A_NAN )
!
!*****************************************************************************
!  Module NAN_DBLE returns TRUE if a REAL*8 number is equal to the IEEE NaN 
! (Not-a-Number) flag.  Returns FALSE otherwise. (bmy, 3/8/01, 6/27/03)
!
!  Arguments as Input:
!  ===========================================================================
!  (1) VALUE (REAL*8) : Number to be tested for NaN
!
!  NOTES:
!  (1 ) Is overloaded by interface "IT_IS_NAN".
!  (2 ) Now call C routine is_nan(x) for Linux platform (bmy, 6/13/02)
!  (3 ) Eliminate IF statement in Linux section.  Also now trap NaN on
!        the Sun/Sparc platform.  Rename cpp switch from DEC_COMPAQ to
!        COMPAQ. (bmy, 3/23/03)
!  (4 ) Added patches for IBM/AIX (gcc, bmy, 6/27/03)
!*****************************************************************************
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
#if   defined( SGI )
      IT_IS_A_NAN = IEEE_IS_NAN( VALUE )   ! for SGI platform

#elif defined( COMPAQ )
      IT_IS_A_NAN = ISNAN( VALUE )         ! for DEC, Compaq platforms

#elif defined( LINUX )

      ! Declare IS_NAN as an external function
      INTEGER, EXTERNAL  :: IS_NAN

      ! For LINUX platform, use C routine "is_nan" to test if VALUE is NaN.  
      IT_IS_A_NAN = ( IS_NAN( VALUE ) /= 0 )

#elif defined( SPARC )

      ! Declare ID_ISNAN as an external function
      INTEGER, EXTERNAL  :: ID_ISNAN

      ! Test if VALUE is NaN
      IT_IS_A_NAN = ( ID_ISNAN( VALUE ) /= 0 )

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
!*****************************************************************************
!  Module FINITE_FLOAT returns TRUE if a REAL*4 number is equal to the 
!  IEEE Infinity flag.  Returns FALSE otherwise. (bmy, 3/8/01, 9/29/03)
!
!  Arguments as Input:
!  ===========================================================================
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
!*****************************************************************************
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
#if   defined( SGI )
      IT_IS_A_FINITE = IEEE_FINITE( VALUE )  ! for SGI platform
 
#elif defined( COMPAQ ) 

      ! Test for finite # using bit masking for DEC/Compaq platform
      ! Now use REAL*4 values (bmy, 11/15/01)
      IF ( VALUE == Z'7F800000' .or. 
     &     VALUE == Z'FF800000' ) THEN
         IT_IS_A_FINITE = .FALSE.
      ELSE
         IT_IS_A_FINITE = .TRUE.
      ENDIF

#elif defined( LINUX )

      !---------------------------------------------------------------------
      ! Prior to 9/29/03:
      !! Declare IS_INF as an external function
      !INTEGER, EXTERNAL  :: IS_INF
      !
      !! For LINUX, use C routine "is_inf" to test if VALUE is infinity  
      !! VALUE must be cast to DBLE since "is_inf" only takes doubles.
      !IT_IS_A_FINITE = ( IS_INF( DBLE( VALUE ) ) /= 0 )
      !---------------------------------------------------------------------
      
      ! Declare IS_FINITE as an external function
      INTEGER, EXTERNAL :: IS_FINITE
      
      ! For LINUX, use C routine "is_finite" to test if VALUE is finite  
      ! VALUE must be cast to DBLE since "is_inf" only takes doubles. 
      IT_IS_A_FINITE = ( IS_FINITE( DBLE( VALUE ) ) /= 0 )

#elif defined( SPARC )

      ! Declare IR_FINITE as an external function
      INTEGER, EXTERNAL :: IR_FINITE

      ! Test if VALUE is a finite number
      IT_IS_A_FINITE = ( IR_FINITE( VALUE ) /= 0 )

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
!  IEEE Infinity flag.  Returns FALSE otherwise. (bmy, 3/8/01, 9/29/03)
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
!*****************************************************************************
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
#if   defined( SGI )
      IT_IS_A_FINITE = IEEE_FINITE( VALUE )  ! for SGI platform

#elif defined( COMPAQ ) 

      ! Test for finite # using bit masking for DEC/Compaq F90
      IF ( VALUE == Z'7FF0000000000000'  .or. 
     &     VALUE == Z'FFF0000000000000') THEN
         IT_IS_A_FINITE = .FALSE.
      ELSE
         IT_IS_A_FINITE = .TRUE.
      ENDIF

#elif defined( LINUX )

      !----------------------------------------------------------------------
      ! Prior to 9/29/03:
      !! Declare IS_INF as an external function
      !INTEGER, EXTERNAL  :: IS_INF
      !
      ! For LINUX, use C routine "is_inf" to test if VALUE is infinity  
      !IT_IS_A_FINITE = ( IS_INF( VALUE ) /= 0 )
      !----------------------------------------------------------------------

      ! Declare IS_FINITE as an external function
      INTEGER, EXTERNAL :: IS_FINITE

      ! For LINUX, use C routine "is_finite" to test if VALUE is infinity
      IT_IS_A_FINITE = ( IS_FINITE( VALUE ) /= 0 )

#elif defined( SPARC )
 
      ! Declare ID_FINITE as an external function
      INTEGER, EXTERNAL :: ID_FINITE

      ! Test if VALUE is a finite number
      IT_IS_A_FINITE = ( ID_FINITE( VALUE ) /= 0 )

#elif defined( IBM_AIX )

      ! For IBM/AIX platform
      IF ( IEEE_SUPPORT_DATATYPE( VALUE ) ) THEN
         IT_IS_A_FINITE = IEEE_IS_FINITE( VALUE )
      ENDIF

#endif
      
      ! Return to calling program
      END FUNCTION FINITE_DBLE
      
!-----------------------------------------------------------------------------

      SUBROUTINE CHECK_STT( LOCATION )
!
!******************************************************************************
!  Subroutine CHECK_STT checks the STT tracer array for negative values,
!  NaN values, or Infinity values.  If any of these are found, the code
!  will stop with an error message. (bmy, 3/8/01, 3/23/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1) LOCATION (CHARACTER) : String describing location of error in code
!
!  NOTES:
!  (1 ) CHECK_STT uses the interfaces defined above -- these will do the
!        proper error checking for either SGI or DEC/Compaq platforms.
!        (bmy, 3/8/01)
!  (2 ) Now call GEOS_CHEM_STOP to shutdown safely.  Now use logicals LNAN,
!        LNEG, LINF to flag if we have error conditions, and then stop the
!        run outside of the parallel DO-loop. (bmy, 11/27/02)
!  (3 ) Bug fix in FORMAT statement: replace missing commas (bmy, 3/23/03)
!******************************************************************************
!
#     include "CMN_SIZE"           ! Size parameters
#     include "CMN"                ! STT, NTRACE

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: LOCATION

      ! Local variables
      LOGICAL                      :: LNEG, LNAN, LINF
      INTEGER                      :: I,    J,    L,   N
      
      !=================================================================
      ! CHECK_STT begins here!
      !=================================================================

      ! Initialize
      LNEG = .FALSE.
      LNAN = .FALSE.
      LINF = .FALSE.

      ! Loop over grid boxes
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, NTRACE
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !---------------------------
         ! Check for Negatives
         !---------------------------
         IF ( STT(I,J,L,N) < 0d0 ) THEN 
!$OMP CRITICAL
            LNEG = .TRUE.
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
!$OMP END CRITICAL

         !---------------------------
         ! Check for NaN's
         !---------------------------
         ELSE IF ( IT_IS_NAN( STT(I,J,L,N) ) ) THEN
!$OMP CRITICAL
            LNAN = .TRUE.
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
!$OMP END CRITICAL

         !----------------------------
         ! Check STT's for Infinities
         !----------------------------
         ELSE IF ( .not. IT_IS_FINITE( STT(I,J,L,N) ) ) THEN
!$OMP CRITICAL
            LINF = .TRUE.
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
!$OMP END CRITICAL            

         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Stop the run if any of LNEG, LNAN, LINF is true
      !=================================================================
      IF ( LNEG .or. LNAN .or. LINF ) THEN
         WRITE( 6, 120 ) TRIM( LOCATION )
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! FORMAT statements
      !=================================================================
 100  FORMAT( 'CHECK_STT: STT(',i3,',',i3,',',i3,',',i3,') = ', f13.6 )
 120  FORMAT( 'CHECK_STT: STOP at ', a )

      ! Return to calling program
      END SUBROUTINE CHECK_STT

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
!  and then stops the run.  (bmy, 10/15/02) 
!******************************************************************************
!
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

      SUBROUTINE ALLOC_ERR( ARRAYNAME )
!
!******************************************************************************
!  Subroutine ALLOC_ERR prints an error message if there is not enough
!  memory to allocate a particular allocatable array. (bmy, 6/26/00, 10/15/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ARRAYNAME (CHARACTER) : name of the allocatable array
!
!  NOTES:
!  (1 ) Split off from "ndxx_setup.f" into a separate routine (bmy, 6/26/00)
!  (2 ) Added to "error_mod.f" (bmy, 10/15/02)
!******************************************************************************
!
      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: ARRAYNAME

      ! Local variables
      CHARACTER(LEN=255)           :: ERRMSG

      !=================================================================
      ! ALLOC_ERR begins here!
      !=================================================================

      ! Define error message
      ERRMSG = 'Allocation error in array: ' // TRIM( ARRAYNAME )
      
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
!  (bmy, 1/7/02, 11/22/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) MESSAGE (CHAR*(*)) : Descriptive message string
!
!  NOTES: 
!  (1 ) Now just write the message and flush the buffer (bmy, 7/5/01)
!  (2 ) Renamed from "paftop.f" to "debug_msg.f" (bmy, 1/7/02)
!  (3 ) Bundled into "error_mod.f" (bmy, 11/22/02)
!******************************************************************************
!
      IMPLICIT NONE

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: MESSAGE

      !=================================================================
      ! DEBUG_MSG begins here!
      !================================================================= 
      WRITE( 6, '(a)' ) TRIM( MESSAGE )
      CALL FLUSH( 6 )

      ! Return to calling program
      END SUBROUTINE DEBUG_MSG

!------------------------------------------------------------------------------

      END MODULE ERROR_MOD
