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
! !PUBLIC DATA MEMBERS:
!
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
      PUBLIC  :: PRINT_GLOBAL_SPECIES_KG
#if defined( MODEL_GEOS )
      PUBLIC  :: mpiComm
#endif 

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

      INTERFACE IS_SAFE_DIV
         MODULE PROCEDURE IS_SAFE_DIV_R4
         MODULE PROCEDURE IS_SAFE_DIV_R8
      END INTERFACE
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
!  08 Jul 2014 - R. Yantosca - Added INIT_ERROR and CLEANUP_ERROR subroutines
!  08 Jul 2014 - R. Yantosca - Add shadow variables for am_I_Root, Input_Opt
!                              so that we can pass these to routine CLEANUP
!  02 Dec 2014 - M. Yannetti - Added PRECISION_MOD 
!  19 Dec 2014 - R. Yantosca - Now overload IS_SAFE_DIV w/ REAL*4 and REAL*8 
!                              module procedures to facilitate flex precision
!  13 Aug 2015 - E. Lundgren - Add GIGC_ERROR which sets RC to GIGC_FAILURE
!                              and print msg and location to log
!  22 Jan 2016 - R. Yantosca - Remove SPARC, IBM #ifdefs
!  17 Aug 2016 - M. Sulprizio- Move CHECK_SPC and CHECK_SPC_NESTED here from
!                              obsolete tracer_mod.F
!  16 Sep 2016 - E. Lundgren - Add routine to print global species mass to log
!  26 Jun 2017 - R. Yantosca - Moved GC_ERROR to Headers/errcode_mod.F90
!  15 Jun 2018 - E. Lundgren - Add code for use in GEOS-5 only
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:

#if defined( MODEL_GEOS )
      INTEGER, SAVE           :: mpiComm
#endif
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
! !USES:
!
#if defined( LINUX_PGI )
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
!  22 Jan 2016 - R. Yantosca - Remove SPARC, IBM #ifdefs
!  22 Jan 2016 - R. Yantosca - Use IEEE_IS_NAN for PGI compiler
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( LINUX_IFORT ) || defined( LINUX_GFORTRAN )
      IT_IS_A_NAN = ISNAN( VALUE )    
#elif defined( LINUX_PGI )
      IT_IS_A_NAN = IEEE_IS_NAN( VALUE )
#endif

      ! Return to calling program
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
! !USES:
!
#if defined( LINUX_PGI )
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
!  22 Jan 2016 - R. Yantosca - Removed SPARC, IBM #ifdefs
!  22 Jan 2016 - R. Yantosca - Use IEEE_IS_NAN for PGI compiler
!EOP
!------------------------------------------------------------------------------
!BOC
!
#if defined( LINUX_IFORT ) || defined( LINUX_GFORTRAN )
      IT_IS_A_NAN = ISNAN( VALUE )      
#elif defined( LINUX_PGI )
      IT_IS_A_NAN = IEEE_IS_NAN( VALUE )
#endif

      ! Return to calling program
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
! !USES:
!
#if defined( LINUX_PGI )
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
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  22 Jan 2016 - R. Yantosca - Removed SPARC, IBM #ifdefs
!  22 Jan 2016 - R. Yantosca - Use IEEE_IS_NAN for PGI compiler
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( LINUX_GFORTRAN )

      IT_IS_A_FINITE = ((.not.ISNAN(VALUE)) .and. 
     &                  (VALUE.lt.HUGE(1.0e0)) .and.
     &                  (VALUE.gt.(-1.0e0*HUGE(1.0e0))))

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

      ! PGI compiler
      IT_IS_A_FINITE = IEEE_IS_FINITE( VALUE )

#endif
      
      ! Return to calling program
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
! !USES:
!
#if defined( LINUX_PGI )
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
!  22 Jan 2016 - R. Yantosca - Removed SPARC, IBM #ifdefs
!  22 Jan 2016 - R. Yantosca - Use IEEE_IS_NAN for PGI compiler
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( LINUX_GFORTRAN )

      IT_IS_A_FINITE = ((.not.ISNAN(VALUE)) .and. 
     &                  (VALUE.lt.HUGE(1.d0)) .and.
     &                  (VALUE.gt.(-1.0d0*HUGE(1.d0))))

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

      ! PGI compiler
      IT_IS_A_FINITE = IEEE_IS_FINITE( VALUE )

#endif
      
      ! Return to calling program
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
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  06 Jan 2017 - R. Yantosca - Added optional INSTRUCTIONS argument
!  01 Nov 2017 - R. Yantosca - Flush the stdout buffer before halting
!  13 Nov 2018 - R. Yantosca - Now word-wrap the Instructions string
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
      
      ! Return to calling program
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
      USE Input_Opt_Mod,      ONLY : OptInput
      USE GEOS_TIMERS_MOD
#if defined( ESMF_ )
      !-----------------------------------------------------------------
      !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
      !
      ! Use GEOS-5 style error reporting when connecting to the GEOS-5
      ! GCM via the ESMF interface (bmy, 3/12/13)
      !-----------------------------------------------------------------
      USE MAPL_Mod
#     include "MAPL_Generic.h"
#endif
!
! !REVISION HISTORY:
!  15 Oct 2002 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Now EXIT works for LINUX_IFC, LINUX_EFC,
!                              so remove #if block.
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  12 Mar 2013 - R. Yantosca - Now use GEOS-5 style traceback when using ESMF
!   8 Jul 2014 - R. Yantosca - Now call cleanup.F with shadow variables
!  11 Aug 2015 - M. Yannetti - Now calls all timers to stop.
!  26 Mar 2016 - S.D.Eastham - Re-ordered to allow __Iam__ to contain
!                              variable declarations
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

#if defined( USE_TIMERS )
      CALL GEOS_Timer_StopAll( RC )
      CALL GEOS_Timer_PrintAll( am_I_Root, RC )  
#endif


#else
      ! Only write to stdout if we are on the root CPU
      am_I_Root = .TRUE. 

#if defined( USE_TIMERS )
      CALL GEOS_Timer_StopAll( RC )
      CALL GEOS_Timer_PrintAll( am_I_Root, RC )  
#endif

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

      ! End of program
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
! !IROUTINE: Safe_Div
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
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
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

      IF ( EXPONENT(N) - EXPONENT(D) >= MaxExp .or. D==0 .or.
     &     EXPONENT(N) - EXPONENT(D) <= MinExp  ) THEN
         F = .FALSE.
      ELSE
         F = .TRUE.
      ENDIF

      ! Return to calling program
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
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
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

      IF ( EXPONENT(N) - EXPONENT(D) >= MaxExp .or. D==0 .or.
     &     EXPONENT(N) - EXPONENT(D) <= MinExp  ) THEN
         F = .FALSE.
      ELSE
         F = .TRUE.
      ENDIF

      ! Return to calling program
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
!  06 Feb 2015 - M. Yannetti - Needed to make the CUTOFF smaller for
!                              REAL*4.
!  22 Jan 2016 - R. Yantosca - Use lowercase for #ifdefs; PGI chokes if not.
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
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
!  error\_mod.F to compensate for the removal of logical\_mod.F from 
!  GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_ERROR( am_I_Root, Input_Opt, RC )
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)         :: am_I_Root ! Are we on root CPU?
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Assume success
      RC               =  GC_SUCCESS

      ! Store a shadow copy of am_I_Root
      SHADOW_am_I_Root =  am_I_Root

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
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Free the pointer to Input_Opt
      NULLIFY( SHADOW_Input_Opt )

      END SUBROUTINE CLEANUP_ERROR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Print_Global_Species_Kg
!
! !DESCRIPTION: Subroutine Print\_Global\_Species\_Kg sums up the total 
!   global species mass for species \#1 and prints it to log along with
!   a user-defined location
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Print_Global_Species_Kg( State_Chm, State_Met, LOC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE State_Chm_Mod, ONLY : ChmState
      USE State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN)    :: LOC
      TYPE(ChmState),   INTENT(IN)    :: State_Chm 
      TYPE(MetState),   INTENT(IN)    :: State_Met 
! 
! !REMARKS:
!  This routine is for debugging purposes to trace where species
!  mass is not conserved 
!
! !REVISION HISTORY: 
!  22 Jun 2016 - E. Lundgren - Initial version
!  29 Mar 2017 - R. Yantosca - Now print out sums for up to the 1st 5 species
!                              and use a more efficient algorithm
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS
!
      INTEGER, PARAMETER :: nSums = 5
!
! !LOCAL VARIABLES:
!     
      ! Scalars
      INTEGER            :: I, J, L, N, nSpecies
      REAL(fp)           :: SpcTotal(nSums)

      ! Arrays
      REAL(fp)           :: Species_Kg    (IIPAR,JJPAR,LLPAR      )
      REAL(fp)           :: Sum_Species_Kg(IIPAR,JJPAR,LLPAR,nSums)

      ! Strings
      CHARACTER(LEN=12)  :: SpcName

      !================================================================
      ! Print_Global_Species_Kg begins here!
      !================================================================

      ! Echo info
      WRITE( 6, 100 ) TRIM( LOC )
 100  FORMAT( /, '%%%%% PRINT_GLOBAL_SPECIES_KG at ', a )

      ! Initialize
      Species_Kg     = 0.0_fp
      Sum_Species_Kg = 0.0_fp
      
      ! Check mass for up to the 1st nSums (=5) species
      nSpecies       = MIN( nSums, State_Chm%nSpecies )

      !---------------------------
      ! Compute global sums [kg]
      !---------------------------
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, nSpecies

         ! Compute sum at each grid box
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            Sum_Species_Kg(I,J,L,N) = Sum_Species_Kg(I,J,L,N)
     &                              + ( State_Chm%Species(I,J,L,N)
     &                              *   State_Met%AD     (I,J,L)   ) 
         ENDDO
         ENDDO
         ENDDO

         ! Compute global sum for each species
         SpcTotal(N) = SUM( Sum_Species_Kg(:,:,:,N) )
      ENDDO
!$OMP END PARALLEL DO

      !---------------------------
      ! Display global sums [kg]
      !---------------------------
      DO N = 1, nSpecies
         
         ! Get species name from the species database
         SpcName = TRIM( State_Chm%SpcData(N)%Info%Name )
         
         ! Write formatted output
         WRITE( 6, 110 ) SpcName, SpcTotal(N)
 110     FORMAT( 'Global sum [kg] for ', a12, ' = ', es24.16 )         
         print*, '===> : ', State_Chm%Species(23,34,5,N)
     &                  *   State_Met%AD     (23,34,5)
      ENDDO

      ! Add trailing space
      WRITE( 6, 120 )
 120  FORMAT( / )

      END SUBROUTINE Print_Global_Species_Kg
!EOC
      END MODULE ERROR_MOD
