!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: error_mod
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

      ! Interface for NaN-check routines
      INTERFACE IT_IS_NAN
         MODULE PROCEDURE NAN_FLOAT
         MODULE PROCEDURE NAN_DBLE
      END INTERFACE

      ! Interface for finite-check routines
      INTERFACE IT_IS_FINITE
         MODULE PROCEDURE FINITE_FLOAT
         MODULE PROCEDURE FINITE_DBLE
      END INTERFACE

      ! Interface for check-value routines
      INTERFACE CHECK_VALUE
         MODULE PROCEDURE CHECK_REAL_VALUE
         MODULE PROCEDURE CHECK_DBLE_VALUE
      END INTERFACE
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: CHECK_DBLE_VALUE
      PRIVATE :: CHECK_REAL_VALUE
      PRIVATE :: FINITE_DBLE
      PRIVATE :: FINITE_FLOAT
      PRIVATE :: NAN_DBLE
      PRIVATE :: NAN_FLOAT
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: "LINUX_IFORT" and "LINUX_PGI" ARE CURRENTLY THE ONLY       %%%
!  %%%  SUPPORTED COMPILER OPTIONS.  MOST SYSTEMS NOW USE A UNIX VERSION %%%
!  %%%  BASED ON LINUX, SO OTHER COMPILERS (IBM/AIX, SUN/SPARC, etc.)    %%%
!  %%%  ARE GENERALLY NOT USED ANYMORE.  LEAVE THE OLDER CODE HERE JUST  %%%
!  %%%  IN CASE WE NEED TO REVERT TO IT AGAIN IN THE FUTURE.             %%%
!  %%%     -- Bob Yantosca, 20 Aug 2013                                  %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  08 Mar 2001 - R. Yantosca - Initial version
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
!  (20) Added routine SAFE_DIV (phs, bmy, 2/26/08)
!  (21) Added routine IS_SAFE_DIV (phs, bmy, 6/11/08)
!  (22) Updated routine SAFE_DIV (phs, 4/14/09)
!  (23) Remove support for SGI, COMPAQ compilers (bmy, 7/8/09)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  04 Jan 2010 - R. Yantosca - Added SAFE_EXP and IS_SAFE_EXP functions
!  04 Jan 2010 - R. Yantosca - Added SAVE_LOG and SAFE_LOG10 functions
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: nan_float
!
! !DESCRIPTION: Function NAN\_FLOAT returns TRUE if a REAL*4 number is equal 
!  to the IEEE NaN (Not-a-Number) flag.  Returns FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION NAN_FLOAT( VALUE ) RESULT( IT_IS_A_NAN )
!
! !USES:
!
#if   defined( IBM_AIX ) || defined( IBM_XLF )
      USE IEEE_ARITHMETIC
#endif
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
!  (10) Remove support for SGI, COMPAQ compilers.  Add IBM_XLF switch. 
!        (bmy, 7/8/09)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( LINUX_IFORT )
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

#elif defined( IBM_AIX ) || defined( IBM_XLF )

      ! For IBM/AIX platform
      IF ( IEEE_SUPPORT_DATATYPE( VALUE ) ) THEN
         IT_IS_A_NAN = IEEE_IS_NAN( VALUE )
      ENDIF

#endif

      ! Return to calling program
      END FUNCTION NAN_FLOAT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: nan_dble
!
! !DESCRIPTION: Function NAN\_DBLE returns TRUE if a REAL*8 number is equal 
!  to the IEEE NaN (Not-a-Number) flag.  Returns FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION NAN_DBLE( VALUE ) RESULT( IT_IS_A_NAN )
!
! !USES:
!
#if   defined( IBM_AIX ) || defined( IBM_XLF )
      USE IEEE_ARITHMETIC
#endif
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
!  (10) Remove support for SGI, COMPAQ compilers.  Add IBM_XLF switch. 
!        (bmy, 7/8/09)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
#if   defined( LINUX_IFORT )
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

#elif defined( IBM_AIX ) || defined( IBM_XLF )

       ! For IBM/AIX platform
      IF ( IEEE_SUPPORT_DATATYPE( VALUE ) ) THEN
         IT_IS_A_NAN = IEEE_IS_NAN( VALUE )
      ENDIF

#endif

      ! Return to calling program
      END FUNCTION NAN_DBLE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: finite_float
!
! !DESCRIPTION: Function FINITE\_FLOAT returns FALSE if a REAL*4 number is 
!  equal to the IEEE Infinity flag.  Returns TRUE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION FINITE_FLOAT( VALUE ) RESULT( IT_IS_A_FINITE )
!
! !USES:
!
#if   defined( IBM_AIX ) || defined( IBM_XLF )
      USE IEEE_ARITHMETIC
#endif
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
!  (11) Remove support for SGI, COMPAQ compilers.  Add IBM_XLF switch. 
!        (bmy, 7/8/09)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( LINUX_IFORT )

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

#elif defined( IBM_AIX ) || defined( IBM_XLF )

      ! For IBM/AIX platform
      IF ( IEEE_SUPPORT_DATATYPE( VALUE ) ) THEN
         IT_IS_A_FINITE = IEEE_IS_FINITE( VALUE )
      ENDIF 
      
#endif
      
      ! Return to calling program
      END FUNCTION FINITE_FLOAT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: finite_dble
!
! !DESCRIPTION: Function FINITE\_FLOAT returns FALSE if a REAL*8 number is 
!  equal to the IEEE Infinity flag.  Returns TRUE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION FINITE_DBLE( VALUE ) RESULT( IT_IS_A_FINITE )
!
! !USES:
!
#if   defined( IBM_AIX ) || defined( IBM_XLF )
      USE IEEE_ARITHMETIC
#endif
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
!  (11) Remove support for SGI, COMPAQ compilers.  Add IBM_XLF switch. 
!        (bmy, 7/8/09)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( LINUX_IFORT )

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

#elif defined( IBM_AIX ) || defined( IBM_XLF )

      ! For IBM/AIX platform
      IF ( IEEE_SUPPORT_DATATYPE( VALUE ) ) THEN
         IT_IS_A_FINITE = IEEE_IS_FINITE( VALUE )
      ENDIF

#endif
      
      ! Return to calling program
      END FUNCTION FINITE_DBLE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_real_value
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
!  15 Oct 2002 - R. Yantosca - Now call GEOS_CHEM_STOP to shutdown safely
!  15 Oct 2002 - R. Yantosca - Updated comments, cosmetic changes
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  10 Jun 2013 - R. Yantosca - Avoid array temporaries, use CHAR*255 args
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
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_dble_value
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
!  15 Oct 2002 - R. Yantosca - Now call GEOS_CHEM_STOP to shutdown safely
!  15 Oct 2002 - R. Yantosca - Updated comments, cosmetic changes
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  10 Jun 2013 - R. Yantosca - Avoid array temporaries, use CHAR*255 args
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
 110  FORMAT( 'CHECK_VALUE: ', a, ' is NaN!'        )
 115  FORMAT( 'Grid box (I,J,L,N) : ', 4i4          )
 120  FORMAT( 'CHECK_VALUE: ', a, ' is not finite!' )

      ! Return to calling program
      END SUBROUTINE CHECK_DBLE_VALUE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: error_stop
!
! !DESCRIPTION: Subroutine ERROR\_STOP is a wrapper for GEOS\_CHEM\_STOP.  It 
!  prints an error message then calls GEOS\_CHEM\_STOP to free memory and quit.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ERROR_STOP( MESSAGE, LOCATION )
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN) :: MESSAGE    ! Error msg to print
      CHARACTER(LEN=*), INTENT(IN) :: LOCATION   ! Where ERROR_STOP is called
!
! !REVISION HISTORY:
!  15 Oct 2002 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC

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
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geos_chem_stop
!
! !DESCRIPTION: Subroutine GEOS\_CHEM\_STOP calls CLEANUP to deallocate all 
!  module arrays and then stops the run.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GEOS_CHEM_STOP
!
! !USES:
!
#if defined( ESMF_ )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! Use GEOS-5 style error reporting when connecting to the GEOS-5
      ! GCM via the ESMF interface (bmy, 3/12/13)
      !-----------------------------------------------------------------
      USE MAPL_Mod
#     include "MAPL_Generic.h"
#endif

! !REVISION HISTORY:
!  15 Oct 2002 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Now EXIT works for LINUX_IFC, LINUX_EFC,
!                              so remove #if block.
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  12 Mar 2013 - R. Yantosca - Now use GEOS-5 style traceback when using ESMF
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( ESMF_ )

      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! Use GEOS-5 style error reporting when connecting to the GEOS-5
      ! GCM via the ESMF interface (bmy, 3/12/13)
      !-----------------------------------------------------------------
      INTEGER :: RC      
      __Iam__('GEOSCHEMSTOP')
      STATUS = -1
      VERIFY_(STATUS)

#else
      !-----------------------------------------------------------------
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! Current practice in the std GEOS-Chem is to call CLEANUP to 
      ! deallocate module arrays and then exit (bmy, 3/12/13)
      !-----------------------------------------------------------------
!$OMP CRITICAL

      ! Deallocate all module arrays
      CALL CLEANUP

      ! Flush all files and stop
      CALL EXIT( 99999 )

!$OMP END CRITICAL

#endif

      ! End of program
      END SUBROUTINE GEOS_CHEM_STOP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: alloc_err
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
!  26 Jun 2000 - R. Yantosca - Initial version, split off from "ndxx_setup.f"
!  15 Oct 2002 - R. Yantosca - Added to "error_mod.f"
!  30 Nov 2005 - R. Yantosca - Call IFORT_ERRMSG for Intel Fortran compiler
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
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
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: debug_msg
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
!  (1 ) Now just write the message and flush the buffer (bmy, 7/5/01)
!  (2 ) Renamed from "paftop.f" to "debug_msg.f" (bmy, 1/7/02)
!  (3 ) Bundled into "error_mod.f" (bmy, 11/22/02)
!  (4 ) Now do not FLUSH the buffer for EFC compiler (bmy, 4/6/04)
!  (5 ) Now add a little space for debug output (bmy, 4/10/06)
!  (6 ) Remove support for LINUX_IFC & LINUX_EFC compilers (bmy, 8/4/06)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Print message
      WRITE( 6, '(5x,a)' ) MESSAGE

      ! Call FLUSH routine to flush the output buffer
      CALL FLUSH( 6 )

      ! Return to calling program
      END SUBROUTINE DEBUG_MSG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: safe_div
!
! !DESCRIPTION: Function SAFE\_DIV performs "safe division", that is to 
!  prevent overflow, underlow, NaN, or infinity errors.  An alternate value 
!  is returned if the division cannot be performed.
!\\
!\\
! !INTERFACE:
!
      FUNCTION SAFE_DIV( N,        D, 
     &                   ALT_NAN,  ALT_OVER, 
     &                   ALT_UNDER            ) RESULT( Q )
!
! !INPUT PARAMETERS: 
!
      REAL*8,           INTENT(IN) :: N           ! Numerator
      REAL*8,           INTENT(IN) :: D           ! Denominator
      REAL*8,           INTENT(IN) :: ALT_NAN     ! Alternate value to be 
                                                  !  returned if the division 
                                                  !  is either NAN (0/0) or 
                                                  !  leads to overflow (i.e., 
                                                  !  a too large number)
      REAL*8, OPTIONAL, INTENT(IN) :: ALT_OVER    ! Alternate value to be 
                                                  !  returned if the division
                                                  !  leads to overflow (default
                                                  !  is ALT_NAN)
      REAL*8, OPTIONAL, INTENT(IN) :: ALT_UNDER   ! Alternate value to be 
                                                  !  returned if the division
                                                  !  leads to underflow 
                                                  !  (default is 0, but you 
                                                  !  could use TINY() if you 
                                                  !  want a non-zero result).
!
! !RETURN VALUE:
!
      REAL*8                       :: Q           ! Output from the division

!
! !REMARKS:
!  For more information, see the discussion on:
!   http://groups.google.com/group/comp.lang.fortran/browse_thread/thread/8b367f44c419fa1d/
!
! !REVISION HISTORY:
!  26 Feb 2008 - P. Le Sager & R. Yantosca - Initial version
!  (1) Now can return different alternate values if NAN (that is 0/0),
!      overflow (that is a too large number), or too small (that is greater
!      than 0 but less than smallest possible number). Default value is
!      zero in case of underflow (phs, 4/14/09)
!  (2) Some compiler options flush underflows to zero (-ftz for IFort).
!       To think about it (phs, 4/14/09)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC

      IF ( N==0 .and. D==0 ) THEN

         ! NAN
         Q = ALT_NAN

      ELSE IF ( EXPONENT(N) - EXPONENT(D) >= MAXEXPONENT(N) .OR.
     &          D==0                                        ) THEN

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

      ! Return to calling program
      END FUNCTION SAFE_DIV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: is_safe_div
!
! !DESCRIPTION: Function IS\_SAFE\_DIV tests for "safe division", that is 
!  check if the division will overflow/underflow or hold NaN.  .FALSE. is 
!  returned if the division cannot be performed. (phs, 6/11/08)
!\\
!\\
! !INTERFACE:
!
      FUNCTION IS_SAFE_DIV( N, D, R4 ) RESULT( F )
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN)           :: N    ! Numerator
      REAL*8,  INTENT(IN)           :: D    ! Denominator
      LOGICAL, INTENT(IN), OPTIONAL :: R4   ! Logical flag to use the limits 
                                            !  of REAL*4 to define underflow
                                            !  or overflow.  Extra defensive.
!
! !OUTPUT PARAMETERS:
!
      LOGICAL                       :: F    ! =F if division isn't allowed
                                            ! =T otherwise
!
! !REMARKS:
!  UnderFlow, OverFlow and NaN are tested for. If you need to
!  differentiate between the three, use the SAFE_DIV (phs, 4/14/09)
!
! !REVISION HISTORY:
!  11 Jun 2008 - P. Le Sager - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER MaxExp, MinExp
      REAL*4  RR

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

      IF ( EXPONENT(N) - EXPONENT(D) >= MaxExp .or. D==0 .or.
     &     EXPONENT(N) - EXPONENT(D) <= MinExp  ) THEN
         F = .FALSE.
      ELSE
         F = .TRUE.
      ENDIF

      ! Return to calling program
      END FUNCTION IS_SAFE_DIV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: safe_exp
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
      REAL*8, INTENT(IN) :: X      ! Argument of EXP
      REAL*8, INTENT(IN) :: ALT    ! Alternate value to be returned
!
! !RETURN VALUE:
!
      REAL*8             :: VALUE  ! Output from the exponential
!
! !REVISION HISTORY:
!  04 Jan 2010 - R. Yantosca - Initial version
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
! !IROUTINE: is_safe_exp
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
      REAL*8, INTENT(IN) :: X    ! Argument to the exponential function
!
! !OUTPUT PARAMETERS:
!
      LOGICAL            :: F    ! =F if exponential isn't allowed
                                 ! =T otherwise
!
! !REMARKS:
!  Empirical testing has revealed that -600 < X < 600 will not result in
!  a floating-point exception on Sun and IFORT compilers.  This is good
!  enough for most purposes.
!
! !REVISION HISTORY:
!  04 Jan 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL*8, PARAMETER :: CUTOFF = 600d0

      ! If -CUTOFF < x < CUTOFF, then it is safe to take EXP( x )
      F = ( ABS( X ) < CUTOFF )

      END FUNCTION IS_SAFE_EXP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: safe_log
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
      REAL*8, INTENT(IN) :: X      ! Argument of LOG
      REAL*8, INTENT(IN) :: ALT    ! Alternate value to be returned
!
! !RETURN VALUE:
!
      REAL*8             :: VALUE  ! Output from the natural logarithm
!
! !REVISION HISTORY:
!  04 Jan 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      IF ( X > 0d0 ) THEN
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
! !IROUTINE: safe_log10
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
      REAL*8, INTENT(IN) :: X      ! Argument of LOG10
      REAL*8, INTENT(IN) :: ALT    ! Alternate value to be returned
!
! !RETURN VALUE:
!
      REAL*8             :: VALUE  ! Output from the natural logarithm
!
! !REVISION HISTORY:
!  04 Jan 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      IF ( X > 0d0 ) THEN
         VALUE = LOG10( X )        ! Take LOG10(x) for positive-definite X
      ELSE 
         VALUE = ALT               ! Otherwise return alternate value
      ENDIF

      END FUNCTION SAFE_LOG10
!EOC
      END MODULE ERROR_MOD
