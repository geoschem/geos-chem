!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: error_mod.F90
!
! !DESCRIPTION: Module ERROR\_MOD contains error checking routines.
!\\
!\\
! !INTERFACE:
!
MODULE ERROR_MOD
!
! !USES:
!
  USE ErrCode_Mod
  USE Input_Opt_Mod,      ONLY : OptInput
  USE PRECISION_MOD            ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: ALLOC_ERR
  PUBLIC  :: CHECK_VALUE
  PUBLIC  :: DEBUG_MSG
  PUBLIC  :: ERROR_STOP
  PUBLIC  :: GEOS_CHEM_STOP
  PUBLIC  :: IS_SAFE_DIV
  PUBLIC  :: IS_SAFE_EXP
  PUBLIC  :: IT_IS_NAN
  PUBLIC  :: IT_IS_FINITE
  PUBLIC  :: SAFE_DIV
  PUBLIC  :: SAFE_EXP
  PUBLIC  :: SAFE_LOG
  PUBLIC  :: SAFE_LOG10
  PUBLIC  :: INIT_ERROR
  PUBLIC  :: CLEANUP_ERROR

  ! Interface for NaN-check routines
  INTERFACE IT_IS_NAN
     MODULE PROCEDURE NAN_FLOAT
     MODULE PROCEDURE NAN_DBLE
  END INTERFACE IT_IS_NAN

  ! Interface for finite-check routines
  INTERFACE IT_IS_FINITE
     MODULE PROCEDURE FINITE_FLOAT
     MODULE PROCEDURE FINITE_DBLE
  END INTERFACE IT_IS_FINITE

  ! Interface for check-value routines
  INTERFACE CHECK_VALUE
     MODULE PROCEDURE CHECK_REAL_VALUE
     MODULE PROCEDURE CHECK_DBLE_VALUE
  END INTERFACE CHECK_VALUE

  INTERFACE IS_SAFE_DIV
     MODULE PROCEDURE IS_SAFE_DIV_R4
     MODULE PROCEDURE IS_SAFE_DIV_R8
  END INTERFACE IS_SAFE_DIV
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: CHECK_DBLE_VALUE
  PRIVATE :: FINITE_DBLE
  PRIVATE :: FINITE_FLOAT
  PRIVATE :: NAN_DBLE
  PRIVATE :: NAN_FLOAT
  PRIVATE :: IS_SAFE_DIV_R4
  PRIVATE :: IS_SAFE_DIV_R8
!
! !REVISION HISTORY:
!  08 Mar 2001 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:

  LOGICAL                 :: SHADOW_am_I_Root   ! Shadow for am_I_Root
  TYPE(OptInput), POINTER :: SHADOW_Input_Opt   ! Shadow for Input_Opt

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nan_Float
!
! !DESCRIPTION: Function NAN\_FLOAT returns TRUE if a REAL*4 number is equal
!  to the IEEE NaN (Not-a-Number) flag.  Returns FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
  FUNCTION NAN_FLOAT( VALUE ) RESULT( IT_IS_A_NAN )
!
! !INPUT PARAMETERS:
!
    REAL*4, INTENT(IN) :: VALUE        ! Value to be tested for NaN
!
! !RETURN VALUE:
!
    LOGICAL            :: IT_IS_A_NAN  ! =T if VALUE is NaN; =F otherwise
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    IT_IS_A_NAN = ISNAN( VALUE )

  END FUNCTION NAN_FLOAT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nan_Dble
!
! !DESCRIPTION: Function NAN\_DBLE returns TRUE if a REAL(fp) number is equal
!  to the IEEE NaN (Not-a-Number) flag.  Returns FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
  FUNCTION NAN_DBLE( VALUE ) RESULT( IT_IS_A_NAN )
!
! !INPUT PARAMETERS:
!
    REAL*8, INTENT(IN) :: VALUE        ! Value to be tested for NaN
!
! !RETURN VALUE:
!
    LOGICAL            :: IT_IS_A_NAN  ! =T if VALUE is NaN; =F otherwise
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    IT_IS_A_NAN = ISNAN( VALUE )

  END FUNCTION NAN_DBLE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finite_Float
!
! !DESCRIPTION: Function FINITE\_FLOAT returns FALSE if a REAL*4 number is
!  equal to the IEEE Infinity flag.  Returns TRUE otherwise.
!\\
!\\
! !INTERFACE:
!
  FUNCTION FINITE_FLOAT( VALUE ) RESULT( IT_IS_A_FINITE )
!
! !INPUT PARAMETERS:
!
    REAL*4, INTENT(IN) :: VALUE           ! Value to be tested for infinity
!
! !RETURN VALUE:
!
    LOGICAL            :: IT_IS_A_FINITE  ! =T if VALUE is finite; =F else
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( LINUX_GFORTRAN )

    IT_IS_A_FINITE = ((.not.ISNAN(VALUE)) .and. &
                      (VALUE.lt.HUGE(1.0e0)) .and. &
                      (VALUE.gt.(-1.0e0*HUGE(1.0e0))))

#elif defined( LINUX_IFORT )

    ! Local variables (parameters copied from "fordef.for")
    INTEGER, PARAMETER :: SNAN=0, QNAN=1, POS_INF=2, NEG_INF=3
    INTEGER            :: FPC

    ! Get the floating point type class for VALUE
    FPC            = FP_CLASS( VALUE )

    ! VALUE is infinite if it is either +Inf or -Inf
    ! Also flag an error if VALUE is a signaling or quiet NaN
    IT_IS_A_FINITE = ( FPC /= POS_INF .and. FPC /= NEG_INF .and. &
                       FPC /= SNAN    .and. FPC /= QNAN          )

#endif

  END FUNCTION FINITE_FLOAT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finite_Dble
!
! !DESCRIPTION: Function FINITE\_FLOAT returns FALSE if a REAL(fp) number is
!  equal to the IEEE Infinity flag.  Returns TRUE otherwise.
!\\
!\\
! !INTERFACE:
!
  FUNCTION FINITE_DBLE( VALUE ) RESULT( IT_IS_A_FINITE )
!
! !INPUT PARAMETERS:
!
    REAL*8, INTENT(IN) :: VALUE           ! Value to be tested for infinity
!
! !RETURN VALUE:
!
    LOGICAL            :: IT_IS_A_FINITE  ! =T if VALUE is finite; =F else
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( LINUX_GFORTRAN )

    IT_IS_A_FINITE = ((.not.ISNAN(VALUE)) .and. &
                      (VALUE.lt.HUGE(1.d0)) .and. &
                      (VALUE.gt.(-1.0d0*HUGE(1.d0))))

#elif defined( LINUX_IFORT )

    ! Local variables (parameters copied from "fordef.for")
    INTEGER, PARAMETER :: SNAN=0, QNAN=1, POS_INF=2, NEG_INF=3
    INTEGER            :: FPC

    ! Get the floating point type class for VALUE
    FPC            = FP_CLASS( VALUE )

    ! VALUE is infinite if it is either +Inf or -Inf
    ! Also flag an error if VALUE is a signaling or quiet NaN
    IT_IS_A_FINITE = ( FPC /= POS_INF .and. FPC /= NEG_INF .and. &
                       FPC /= SNAN    .and. FPC /= QNAN          )

#endif

  END FUNCTION FINITE_DBLE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check_Real_Value
!
! !DESCRIPTION: Subroutine CHECK\_REAL\_VALUE checks to make sure a REAL*4
!  value is not NaN or Infinity. This is a wrapper for the interfaces
!  IT\_IS\_NAN and IT\_IS\_FINITE.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHECK_REAL_VALUE( VALUE, LOCATION, VARNAME, MESSAGE )
!
! !INPUT PARAMETERS:
!
    REAL*4,             INTENT(IN) :: VALUE        ! Value to be checked
    CHARACTER(LEN=255), INTENT(IN) :: VARNAME      ! Name of variable
    CHARACTER(LEN=255), INTENT(IN) :: MESSAGE      ! Short descriptive msg
    INTEGER,            INTENT(IN) :: LOCATION(4)  ! (/ I, J, L, N /) indices
!
! !REVISION HISTORY:
!  13 Jun 2001 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    ! First check for NaN -- print info & stop run if found
    IF ( IT_IS_NAN( VALUE ) ) THEN
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       WRITE( 6, 110   ) TRIM( VARNAME )
       WRITE( 6, 115   ) LOCATION
       WRITE( 6, '(a)' ) TRIM( MESSAGE )
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       CALL GEOS_CHEM_STOP()
    ENDIF

    ! Next check for infinity -- print info & stop run if found
    IF ( .not. IT_IS_FINITE( VALUE ) ) THEN
       WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
       WRITE( 6, 120       ) TRIM( VARNAME )
       WRITE( 6, 115       ) LOCATION
       WRITE( 6, '(f13.6)' ) VALUE
       WRITE( 6, '(a)'     ) TRIM ( MESSAGE )
       WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
       CALL GEOS_CHEM_STOP()
    ENDIF

    ! FORMAT statements
110 FORMAT( 'CHECK_VALUE: ', a, ' is NaN!'        )
115 FORMAT( 'Grid box (I,J,L,N) : ', 4i4          )
120 FORMAT( 'CHECK_VALUE: ', a, ' is not finite!' )

  END SUBROUTINE CHECK_REAL_VALUE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check_Dble_Value
!
! !DESCRIPTION: Subroutine CHECK\_DBLE\_VALUE checks to make sure a REAL*4
!  value is not NaN or Infinity. This is a wrapper for the interfaces
!  IT\_IS\_NAN and IT\_IS\_FINITE.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHECK_DBLE_VALUE( VALUE, LOCATION, VARNAME, MESSAGE )
!
! !INPUT PARAMETERS:
!
    REAL*8,             INTENT(IN) :: VALUE        ! Value to be checked
    CHARACTER(LEN=255), INTENT(IN) :: VARNAME      ! Name of variable
    CHARACTER(LEN=255), INTENT(IN) :: MESSAGE      ! Short descriptive msg
    INTEGER,            INTENT(IN) :: LOCATION(4)  ! (/ I, J, L, N /) indices
!
! !REVISION HISTORY:
!  13 Jun 2001 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
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
110 FORMAT( 'CHECK_VALUE: ', a, ' is NaN!'        )
115 FORMAT( 'Grid box (I,J,L,N) : ', 4i4          )
120 FORMAT( 'CHECK_VALUE: ', a, ' is not finite!' )

  END SUBROUTINE CHECK_DBLE_VALUE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Error_Stop
!
! !DESCRIPTION: Subroutine ERROR\_STOP is a wrapper for GEOS\_CHEM\_STOP.  It
!  prints an error message then calls GEOS\_CHEM\_STOP to free memory and quit.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ERROR_STOP( MESSAGE, LOCATION, INSTRUCTIONS )
!
! !USES:
!
    USE CharPak_Mod, ONLY : WordWrapPrint
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: MESSAGE      ! Error msg to print
    CHARACTER(LEN=*), INTENT(IN) :: LOCATION     ! Where ERROR_STOP is called
    CHARACTER(LEN=*), OPTIONAL   :: INSTRUCTIONS ! Further instructions
!
! !REVISION HISTORY:
!  15 Oct 2002 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512) :: TmpMsg

    !$OMP CRITICAL

    ! Write the error message
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    TmpMsg = 'GEOS-CHEM ERROR: ' // TRIM( MESSAGE )
    CALL WordWrapPrint( TmpMsg, 76 )

    ! If the optional INSTRUCTIONS argument is passed, then print it to
    ! the stdout stream.  This is useful for instructing the user to
    ! look for the error in another location (e.g. the HEMCO log file).
    IF ( PRESENT( INSTRUCTIONS ) ) THEN
       WRITE( 6, '(a)' )
       CALL WordWrapPrint( Instructions, 76 )
       WRITE( 6, '(a)' )
    ENDIF

    ! Write the location of the error
    WRITE( 6, '(a)' ) 'STOP at ' // TRIM( LOCATION )
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )

    !$OMP END CRITICAL

    ! Flush text to file before stopping
    CALL FLUSH( 6 )

    ! Deallocate memory and stop the run
    CALL GEOS_CHEM_STOP()

  END SUBROUTINE ERROR_STOP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Geos_Chem_Stop
!
! !DESCRIPTION: Subroutine GEOS\_CHEM\_STOP calls CLEANUP to deallocate all
!  module arrays and then stops the run.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS_CHEM_STOP()
!
! !USES:
!
    USE ErrCode_Mod
    USE Timers_Mod
#if defined( ESMF_ )
    !-----------------------------------------------------------------
    !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
    !
    ! Use GEOS-5 style error reporting when connecting to the GEOS-5
    ! GCM via the ESMF interface (bmy, 3/12/13)
    !-----------------------------------------------------------------
    USE MAPL_Mod
#   include "MAPL_Generic.h"
#endif
!
! !REVISION HISTORY:
!  15 Oct 2002 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    INTEGER :: RC          ! Success / Failure
    LOGICAL :: am_I_Root   ! Is this the root CPU?

#if defined( ESMF_ )
    !-----------------------------------------------------------------
    !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
    !
    ! Use GEOS-5 style error reporting when connecting to the GEOS-5
    ! GCM via the ESMF interface (bmy, 3/12/13)
    !-----------------------------------------------------------------
    __Iam__('GEOSCHEMSTOP')

    ! Only write to stdout if we are on the root CPU
    am_I_Root = .TRUE.

    IF ( SHADOW_Input_Opt%useTimers ) THEN
       CALL Timer_StopAll( RC )
       CALL Timer_PrintAll( SHADOW_Input_Opt, RC )
    ENDIF

#else
    ! Only write to stdout if we are on the root CPU
    am_I_Root = .TRUE.

    IF ( SHADOW_Input_Opt%useTimers ) THEN
       CALL Timer_StopAll( RC )
       CALL Timer_PrintAll( SHADOW_Input_Opt, RC )
    ENDIF

    !-----------------------------------------------------------------
    !         %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
    !
    ! Current practice in the std GEOS-Chem is to call CLEANUP to
    ! deallocate module arrays and then exit (bmy, 3/12/13)
    !-----------------------------------------------------------------
    !$OMP CRITICAL

    ! Deallocate all module arrays
    CALL CLEANUP( SHADOW_am_I_Root, SHADOW_Input_Opt, .TRUE., RC )

    ! Flush all files and stop
    CALL EXIT( 99999 )

    !$OMP END CRITICAL

#endif

  END SUBROUTINE GEOS_CHEM_STOP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Alloc_Err
!
! !DESCRIPTION: Subroutine ALLOC\_ERR prints an error message if there is not
!  enough memory to allocate a particular allocatable array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ALLOC_ERR( ARRAYNAME, AS )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*),  INTENT(IN) :: ARRAYNAME  ! Name of array
    INTEGER, OPTIONAL, INTENT(IN) :: AS         ! Error output from "STAT"
!
! !REVISION HISTORY:
!  26 Jun 2000 - R. Yantosca - Initial version, split off from "ndxx_setup.F90"
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)            :: ERRMSG

    !=================================================================
    ! ALLOC_ERR begins here!
    !=================================================================

#if defined( LINUX_IFORT )

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
    CALL ERROR_STOP( ERRMSG, 'ALLOC_ERR in error_mod.F90' )

  END SUBROUTINE ALLOC_ERR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Debug_Msg
!
! !DESCRIPTION: Subroutine DEBUG\_MSG prints a message to the stdout buffer
!  and flushes.  This is useful for determining the exact location where
!  errors occur.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DEBUG_MSG( MESSAGE )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: MESSAGE   ! Message to print
!
! !REVISION HISTORY:
!  07 Jan 2002 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Print message
    WRITE( 6, '(5x,a)' ) MESSAGE

    ! Call FLUSH routine to flush the output buffer
    CALL FLUSH( 6 )

  END SUBROUTINE DEBUG_MSG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Safe_Div
!
! !DESCRIPTION: Function SAFE\_DIV performs "safe division", that is to
!  prevent overflow, underlow, NaN, or infinity errors.  An alternate value
!  is returned if the division cannot be performed.
!\\
!\\
! !INTERFACE:
!
  FUNCTION SAFE_DIV( N, D, ALT_NAN, ALT_OVER, ALT_UNDER ) &
       RESULT( Q )
!
! !INPUT PARAMETERS:
!
    REAL(fp),           INTENT(IN) :: N         ! Numerator
    REAL(fp),           INTENT(IN) :: D         ! Denominator
    REAL(fp),           INTENT(IN) :: ALT_NAN   ! Alternate value to be
                                                !  returned if the division
                                                !  is either NAN (0/0) or
                                                !  leads to overflow (i.e.,
                                                !  a too large number)
    REAL(fp), OPTIONAL, INTENT(IN) :: ALT_OVER  ! Alternate value to be
                                                !  returned if the division
                                                !  leads to overflow (default
                                                !  is ALT_NAN)
    REAL(fp), OPTIONAL, INTENT(IN) :: ALT_UNDER ! Alternate value to be
                                                !  returned if the division
                                                !  leads to underflow
                                                !  (default is 0, but you
                                                !  could use TINY() if you
                                                !  want a non-zero result).
!
! !RETURN VALUE:
!
    REAL(fp)                       :: Q         ! Output from the division

!
! !REMARKS:
!  For more information, see the discussion on:
!   http://groups.google.com/group/comp.lang.fortran/browse_thread/thread/8b367f44c419fa1d/
!
! !REVISION HISTORY:
!  26 Feb 2008 - P. Le Sager & R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( N==0 .and. D==0 ) THEN

       ! NAN
       Q = ALT_NAN

    ELSE IF ( EXPONENT(N) - EXPONENT(D) >= MAXEXPONENT(N) .OR. D==0 ) THEN

       ! OVERFLOW
       Q = ALT_NAN
       IF ( PRESENT(ALT_OVER) ) Q = ALT_OVER

    ELSE IF ( EXPONENT(N) - EXPONENT(D) <= MINEXPONENT(N) ) THEN

       ! UNDERFLOW
       Q = 0D0
       IF ( PRESENT(ALT_UNDER) ) Q = ALT_UNDER

    ELSE

       ! No problem
       Q = N / D

    ENDIF

  END FUNCTION SAFE_DIV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Is_Safe_Div_r4
!
! !DESCRIPTION: Function IS\_SAFE\_DIV tests for "safe division", that is
!  check if the division will overflow/underflow or hold NaN.  .FALSE. is
!  returned if the division cannot be performed.  The numerator and
!  denominator must be 4-byte floating point.
!\\
!\\
! !INTERFACE:
!
  FUNCTION IS_SAFE_DIV_R4( N, D, R4 ) RESULT( F )
!
! !INPUT PARAMETERS:
!
    REAL(f4), INTENT(IN)           :: N    ! Numerator
    REAL(f4), INTENT(IN)           :: D    ! Denominator
    LOGICAL,  INTENT(IN), OPTIONAL :: R4   ! Logical flag to use the limits
                                           !  of REAL*4 to define underflow
                                           !  or overflow.  Extra defensive.
!
! !OUTPUT PARAMETERS:
!
    LOGICAL                        :: F    ! =F if division isn't allowed
                                           ! =T otherwise
!
! !REMARKS:
!  UnderFlow, OverFlow and NaN are tested for. If you need to
!  differentiate between the three, use the SAFE_DIV (phs, 4/14/09)
!
! !REVISION HISTORY:
!  11 Jun 2008 - P. Le Sager - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: MaxExp, MinExp
    REAL(f4) :: RR

    !==================================================================
    ! IS_SAFE_DIV begins here!
    !==================================================================

    MaxExp = MAXEXPONENT( N )
    MinExp = MINEXPONENT( N )

    IF ( PRESENT( R4 ) ) THEN
       IF ( R4 ) THEN
          MaxExp = MAXEXPONENT( RR )
          MinExp = MINEXPONENT( RR )
       ENDIF
    ENDIF

    IF ( EXPONENT(N) - EXPONENT(D) >= MaxExp .or. D==0 .or. &
         EXPONENT(N) - EXPONENT(D) <= MinExp  ) THEN
       F = .FALSE.
    ELSE
       F = .TRUE.
    ENDIF

  END FUNCTION IS_SAFE_DIV_R4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Is_Safe_Div_r8
!
! !DESCRIPTION: Function IS\_SAFE\_DIV tests for "safe division", that is
!  check if the division will overflow/underflow or hold NaN.  .FALSE. is
!  returned if the division cannot be performed.  The numerator and
!  denominator must be 4-byte floating point.
!\\
!\\
! !INTERFACE:
!
  FUNCTION IS_SAFE_DIV_R8( N, D, R4 ) RESULT( F )
!
! !INPUT PARAMETERS:
!
    REAL(f8), INTENT(IN)           :: N    ! Numerator
    REAL(f8), INTENT(IN)           :: D    ! Denominator
    LOGICAL,  INTENT(IN), OPTIONAL :: R4   ! Logical flag to use the limits
                                           !  of REAL*4 to define underflow
                                           !  or overflow.  Extra defensive.
!
! !OUTPUT PARAMETERS:
!
    LOGICAL                        :: F    ! =F if division isn't allowed
                                           ! =T otherwise
!
! !REMARKS:
!  UnderFlow, OverFlow and NaN are tested for. If you need to
!  differentiate between the three, use the SAFE_DIV (phs, 4/14/09)
!
! !REVISION HISTORY:
!  11 Jun 2008 - P. Le Sager - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: MaxExp, MinExp
    REAL(f4) :: RR

    !==================================================================
    ! IS_SAFE_DIV begins here!
    !==================================================================

    MaxExp = MAXEXPONENT( N )
    MinExp = MINEXPONENT( N )

    IF ( PRESENT( R4 ) ) THEN
       IF ( R4 ) THEN
          MaxExp = MAXEXPONENT( RR )
          MinExp = MINEXPONENT( RR )
       ENDIF
    ENDIF

    IF ( EXPONENT(N) - EXPONENT(D) >= MaxExp .or. D==0 .or. &
         EXPONENT(N) - EXPONENT(D) <= MinExp  ) THEN
       F = .FALSE.
    ELSE
       F = .TRUE.
    ENDIF

  END FUNCTION IS_SAFE_DIV_R8
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Safe_Exp
!
! !DESCRIPTION: Function SAFE\_EXP performs a "safe exponential", that is to
!  prevent overflow, underlow, NaN, or infinity errors when taking the
!  value EXP( x ).  An alternate value is returned if the exponential
!  cannot be performed.
!\\
!\\
! !INTERFACE:
!
  FUNCTION SAFE_EXP( X, ALT ) RESULT( VALUE )
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN) :: X      ! Argument of EXP
    REAL(fp), INTENT(IN) :: ALT    ! Alternate value to be returned
!
! !RETURN VALUE:
!
    REAL(fp)             :: VALUE  ! Output from the exponential
!
! !REVISION HISTORY:
!  04 Jan 2010 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( IS_SAFE_EXP( X ) ) THEN
       VALUE = EXP( X )
    ELSE
       VALUE = ALT
    ENDIF

  END FUNCTION SAFE_EXP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Is_Safe_Exp
!
! !DESCRIPTION: Function IS\_SAFE\_EXP returns TRUE if it is safe to take
!  the value EXP( x ) without encountering a floating point exception.  FALSE
!  is returned if the exponential cannot be performed.
!\\
!\\
! !INTERFACE:
!
  FUNCTION IS_SAFE_EXP( X ) RESULT( F )
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN) :: X    ! Argument to the exponential function
!
! !OUTPUT PARAMETERS:
!
    LOGICAL              :: F    ! =F if exponential isn't allowed
                                 ! =T otherwise
!
! !REMARKS:
!  Empirical testing has revealed that -600 < X < 600 will not result in
!  a floating-point exception on Sun and IFORT compilers.  This is good
!  enough for most purposes.
!
! !REVISION HISTORY:
!  04 Jan 2010 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
#if defined( USE_REAL8 )
    REAL(fp), PARAMETER :: CUTOFF = 600e+0_fp
#else
    REAL(fp), PARAMETER :: CUTOFF = 75e+0_fp
#endif

    ! If -CUTOFF < x < CUTOFF, then it is safe to take EXP( x )
    F = ( ABS( X ) < CUTOFF )

  END FUNCTION IS_SAFE_EXP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Safe_Log
!
! !DESCRIPTION: Function SAFE\_LOG performs a "safe natural logarithm", that
!  is to prevent overflow, underlow, NaN, or infinity errors when taking the
!  value LOG( x ).  An alternate value is returned if the logarithm
!  cannot be performed.
!\\
!\\
! !INTERFACE:
!
  FUNCTION SAFE_LOG( X, ALT ) RESULT( VALUE )
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN) :: X      ! Argument of LOG
    REAL(fp), INTENT(IN) :: ALT    ! Alternate value to be returned
!
! !RETURN VALUE:
!
    REAL(fp)             :: VALUE  ! Output from the natural logarithm
!
! !REVISION HISTORY:
!  04 Jan 2010 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( X > 0e+0_fp ) THEN
       VALUE = LOG( X )          ! Take LOG(x) for positive-definite X
    ELSE
       VALUE = ALT               ! Otherwise return alternate value
    ENDIF

  END FUNCTION SAFE_LOG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Safe_Log10
!
! !DESCRIPTION: Function SAFE\_LOG10 performs a "safe log10", that
!  is to prevent overflow, underlow, NaN, or infinity errors when taking the
!  value LOG10( x ).  An alternate value is returned if the logarithm
!  cannot be performed.
!\\
!\\
! !INTERFACE:
!
  FUNCTION SAFE_LOG10( X, ALT ) RESULT( VALUE )
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN) :: X      ! Argument of LOG10
    REAL(fp), INTENT(IN) :: ALT    ! Alternate value to be returned
!
! !RETURN VALUE:
!
    REAL(fp)             :: VALUE  ! Output from the natural logarithm
!
! !REVISION HISTORY:
!  04 Jan 2010 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( X > 0e+0_fp ) THEN
       VALUE = LOG10( X )        ! Take LOG10(x) for positive-definite X
    ELSE
       VALUE = ALT               ! Otherwise return alternate value
    ENDIF

  END FUNCTION SAFE_LOG10
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Error
!
! !DESCRIPTION: Subroutine INIT\_ERROR stores shadow copies of am\_I\_Root
!  and Input\_Opt.  We need store shadow copies of these variables within
!  error\_mod.F90 to compensate for the removal of logical\_mod.F from
!  GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_ERROR( Input_Opt, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN), TARGET :: Input_Opt ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)        :: RC        ! Success or failure?
!
! !REMARKS:
!  Instead of making a copy of Input_Opt, we use a pointer reference.
!  This should be more efficient memory-wise.
!
! !REVISION HISTORY:
!  04 Jan 2010 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Assume success
    RC               =  GC_SUCCESS

    ! Store a shadow copy of am_I_Root
    SHADOW_am_I_Root =  Input_Opt%amIRoot

    ! Store a shadow copy of Input_Opt (point to it instead of copying)
    SHADOW_Input_Opt => Input_Opt

  END SUBROUTINE INIT_ERROR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Error
!
! !DESCRIPTION: Subroutine CLEANUP\_ERROR finalizes all module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_ERROR()
!
! !REVISION HISTORY:
!  04 Jan 2010 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Free the pointer to Input_Opt
    NULLIFY( SHADOW_Input_Opt )

  END SUBROUTINE CLEANUP_ERROR
!EOC
END MODULE ERROR_MOD
