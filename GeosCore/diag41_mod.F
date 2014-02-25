!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag41_mod
!
! !DESCRIPTION: Module DIAG41\_MOD contains arrays and routines for 
!  archiving the ND41 diagnostic -- Afternoon PBL heights.
!\\
!\\
! !INTERFACE: 
!
      MODULE DIAG41_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      INTEGER, PUBLIC            :: ND41
      INTEGER, PUBLIC, PARAMETER :: PD41 = 2
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC :: ZERO_DIAG41
      PUBLIC :: WRITE_DIAG41
      PUBLIC :: DIAG41
      PUBLIC :: INIT_DIAG41
      PUBLIC :: CLEANUP_DIAG41
!
! !REVISION HISTORY:
!  17 Feb 2005 - R. Yantosca - Initial version
!  (1 ) Updated for GCAP grid (bmy, 6/28/05)
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (3 ) Replace TINY(1d0) with 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Use updated GET_LOCALTIME from time_mod.F
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      INTEGER, ALLOCATABLE :: GOOD_CT(:)     ! Counter of afternoon data pts
      REAL*4,  ALLOCATABLE :: AD41(:,:,:)    ! Array for afternoon PBL ht

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: zero_diag41
!
! !DESCRIPTION: Subroutine ZERO\_DIAG41 zeroes the ND41 diagnostic arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ZERO_DIAG41
!
! !USES:
!
      USE CMN_SIZE_MOD  ! Size parameters
! 
! !REVISION HISTORY: 
!  17 Feb 2005 - R. Yantosca - Initial version
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, N

      !=================================================================
      ! ZERO_DIAG41 begins here!
      !=================================================================

      ! Exit if ND41 is turned off
      IF ( ND41 == 0 ) RETURN

      ! Zero GOOD_CT
      DO I = 1, IIPAR
         GOOD_CT(I) = 0
      ENDDO

      ! Zero AD41
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N )
      DO N = 1, PD41
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         AD41(I,J,N) = 0e0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE ZERO_DIAG41
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_diag41
!
! !DESCRIPTION: Subroutine WRITE\_DIAG41 writes the ND41 diagnostic arrays 
!  to the binary punch file at the proper time.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE_DIAG41
!
! !USES:
!
      USE BPCH2_MOD, ONLY : BPCH2
      USE BPCH2_MOD, ONLY : GET_HALFPOLAR
      USE BPCH2_MOD, ONLY : GET_MODELNAME 
      USE FILE_MOD,  ONLY : IU_BPCH
      USE GRID_MOD,  ONLY : GET_XOFFSET
      USE GRID_MOD,  ONLY : GET_YOFFSET
      USE TIME_MOD,  ONLY : GET_CT_EMIS
      USE TIME_MOD,  ONLY : GET_DIAGb
      USE TIME_MOD,  ONLY : GET_DIAGe

      USE CMN_SIZE_MOD  ! Size parameters
      USE CMN_DIAG_MOD  ! TINDEX
!
! !REMARKS:
!  ND41: Afternoon PBL depth (between 1200 and 1600 Local Time)
!                                                                             .
!   #  Field    : Description             : Units  : Scale factor
!  -------------------------------------------------------------------
!  (1) PBLDEPTH : Afternoon PBL heights   : m      : GOOD_CT
! 
! !REVISION HISTORY:
!  17 Feb 2005 - R. Yantosca - Initial version
!  (1 ) Now call GET_HALFPOLAR from "bpch2_mod.f" to get the HALFPOLAR flag 
!        value for GEOS or GCAP grids. (bmy, 6/28/05)
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (3 ) Replace TINY(1d0) with 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I,         J,           M,      N
      INTEGER            :: CENTER180, HALFPOLAR,   IFIRST
      INTEGER            :: JFIRST,    LFIRST,      LMAX
      REAL*4             :: ARRAY(IIPAR,JJPAR,LLPAR)
      REAL*4             :: LONRES,    LATRES,      EPS
      REAL*8             :: DIAGb,     DIAGe,       SCALE
      CHARACTER(LEN=20)  :: MODELNAME 
      CHARACTER(LEN=40)  :: CATEGORY,  RESERVED,    UNIT

      !=================================================================
      ! WRITE_DIAG41 begins here!
      !=================================================================

      ! Exit if ND41 is turned off
      IF ( ND41 == 0 ) RETURN

      ! Initialize
      CATEGORY  = 'PBLDEPTH'
      CENTER180 = 1
      DIAGb     = GET_DIAGb()
      DIAGe     = GET_DIAGe()
      HALFPOLAR = GET_HALFPOLAR()
      IFIRST    = GET_XOFFSET( GLOBAL=.TRUE. ) + 1
      JFIRST    = GET_YOFFSET( GLOBAL=.TRUE. ) + 1
      LATRES    = DJSIZE
      LFIRST    = 1
      LONRES    = DISIZE
      MODELNAME = GET_MODELNAME()
      RESERVED  = ''
      EPS       = 1d-32
         
      !=================================================================
      ! Write data to the bpch file
      !=================================================================
      
      ! Loop over ND41 diagnostic tracers
      DO M = 1, TMAX(41)
         N = TINDEX(41,M)
         IF ( N > PD41 ) CYCLE

         ! Select proper unit string
         IF ( N == 1 ) UNIT = 'm' 
         IF ( N == 2 ) UNIT = 'level'
                     
         ! Divide by # of afternoon boxes at each longitude
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            SCALE        = DBLE( GOOD_CT(I) ) + EPS
            ARRAY(I,J,1) = AD41(I,J,N)        / SCALE
         ENDDO
         ENDDO

         ! Write to bpch file
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )
      ENDDO

      END SUBROUTINE WRITE_DIAG41
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diag41
!
! !DESCRIPTION: Subroutine DIAG41 produces monthly mean boundary layer 
!  height in meters  between 1200-1600 local time for the U.S. geographical 
!  domain. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DIAG41 
!
! !USES:
!
      USE PBL_MIX_MOD, ONLY : GET_PBL_TOP_L
      USE PBL_MIX_MOD, ONLY : GET_PBL_TOP_m
      USE TIME_MOD,    ONLY : GET_LOCALTIME

      USE CMN_SIZE_MOD    ! Size parameters
! 
! !REVISION HISTORY:
!  18 Nov 1999 - A. Fiore, S. Wu - Initial version
!  (1 ) DIAG41 is written in Fixed-Format F90. 
!  (2 ) XTRA2 must be computed by turning TURBDAY on first.  Also,
!        XTRA2 is a global-size array, so use window offsets IREF, JREF
!        to index it correctly. (bmy, 11/18/99)
!  (3 ) Do a little rewriting so that the DO-loops get executed
!        in the correct order (J first, then I). (bmy, 11/18/99)
!  (4 ) AD41 is now declared allocatable in "diag_mod.f". (bmy, 12/6/99)
!  (5 ) AFTTOT is now declared allocatable in "diag_mod.f". (bmy, 3/17/00)
!  (6 ) Remove NYMD from the argument list -- it wasn't used (bmy, 6/22/00) 
!  (7 ) XTRA2(IREF,JREF,5) is now XTRA2(I,J).  Also updated comments. 
!        (bmy, 9/25/01)
!  (8 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (9 ) Now reference BXHEIGHT from "dao_mod.f".  Also removed obsolete
!        code. (bmy, 9/18/02)
!  (10) Now use function GET_LOCALTIME from "dao_mod.f" (bmy, 2/11/03)
!  (11) Bug fix in DO-loop for calculating local time (bmy, 7/9/03)
!  (12) For GEOS-4, PBL depth is already in meters, so we only have to
!        multiply that by the GOOD array.  Also now references PBL array
!        from "dao_mod.f".  Bug fix: now use barometric law to compute PBL 
!        height in meters for GEOS-1, GEOS-STRAT, GEOS-3.  This eliminates an 
!        overprediction of the PBL height. (swu, bmy, 11/6/03)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_LOCALTIME(I,J,L) from time_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, N, GOOD(IIPAR)
      REAL*8  :: LT, PBLTOP

      !=================================================================
      ! DIAG41 begins here!
      !=================================================================

      !-----------------------------------
      ! Find boxes where it is afternoon
      !-----------------------------------
      DO I = 1, IIPAR

         ! Local time
         LT = GET_LOCALTIME( I, 1, 1 )
   
         ! Find points between 12 and 16 GMT
         IF ( LT >= 12d0 .and. LT <= 16d0 ) THEN
            GOOD(I) = 1
         ELSE
            GOOD(I) = 0
         ENDIF

         ! Increment counter of afternoon boxes
         GOOD_CT(I)  = GOOD_CT(I) + GOOD(I)
      ENDDO

      !-----------------------------------
      ! Archive afternoon PBL heights
      !-----------------------------------
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N, PBLTOP )
      DO N = 1, PD41
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         IF ( N == 1 ) THEN

            ! Afternoon PBL top [m]
            PBLTOP = GET_PBL_TOP_m( I, J ) * GOOD(I)

         ELSE IF ( N == 2 ) THEN

            ! Afternoon PBL top [model layers]
            PBLTOP = GET_PBL_TOP_L( I, J ) * GOOD(I)

         ENDIF
           
         ! Store in AD41 array
         AD41(I,J,N) = AD41(I,J,N) + PBLTOP

      ENDDO    
      ENDDO    
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE DIAG41
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_diag41
!
! !DESCRIPTION: Subroutine CLEANUP\_DIAG41 allocates and zeroes all 
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_DIAG41
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR
   
      USE CMN_SIZE_MOD  ! Size parameters
! 
! !REVISION HISTORY: 
!  17 Feb 2005 - R. Yantosca - Initial version
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS
      
      !=================================================================
      ! INIT_DIAG41 begins here!
      !=================================================================

      ! Exit if ND41 is turned off
      IF ( ND41 == 0 ) RETURN

      ! Counter of afternoon pts
      ALLOCATE( GOOD_CT( IIPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GOOD_CT' )

      ! Diagnostic array
      ALLOCATE( AD41( IIPAR, JJPAR, PD41 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD41' )

      ! Zero arrays
      CALL ZERO_DIAG41

      END SUBROUTINE INIT_DIAG41
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_diag41
!
! !DESCRIPTION: Subroutine CLEANUP\_DIAG41 deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_DIAG41
! 
! !REVISION HISTORY: 
!  17 Feb 2005 - R. Yantosca - Initial version
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_DIAG41 begins here!
      !=================================================================
      IF ( ALLOCATED( AD41    ) ) DEALLOCATE( AD41    ) 
      IF ( ALLOCATED( GOOD_CT ) ) DEALLOCATE( GOOD_CT )

      END SUBROUTINE CLEANUP_DIAG41
!EOC
      END MODULE DIAG41_MOD
