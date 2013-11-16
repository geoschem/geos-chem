!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: geoeng_mod
!     
! !DESCRIPTION: Module GEOENG\_MOD contains routines and variables which
!  are associated with SO2 or SO4 geoengineering emissions
!\\   
!\\   
! !INTERFACE: 
!
      MODULE GEOENG_MOD
!
! !USES:
!
      USE inquireMod, ONLY : findFreeLUN
      USE ERROR_MOD,  ONLY : DEBUG_MSG
      USE LOGICAL_MOD,ONLY : LPRT
      USE CMN_SIZE_MOD

      IMPLICIT NONE
#     include "define.h"
      PRIVATE

!
! !PUBLIC DATA MEMBERS:
!
      REAL*8,             PUBLIC :: GE_TGS ! TgS/year to deliver
      REAL*8,             PUBLIC :: GE_PBASE ! Geoengineering delivery base (hPa)
      REAL*8,             PUBLIC :: GE_PCEIL ! Geoengineering delivery ceiling (hPa)
      REAL*8,             PUBLIC :: GE_MINLAT ! Lower latitude edge (deg)
      REAL*8,             PUBLIC :: GE_MAXLAT ! Upper latitude edge (deg)
      REAL*8,             PUBLIC :: GE_TSURF ! Surface temperature offsett (K)
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: EMIT_GEOENG
      PUBLIC  :: OFFSET_MET_FIELDS
      PUBLIC  :: INIT_GEOENG
      PUBLIC  :: CLEANUP_GEOENG
!
! !REVISION HISTORY: 
!  11 Nov 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      !=================================================================
      ! MODULE PARAMETERS
      !
      ! AVOGADRO        : Avogadro's number (#/mole)
      ! R_UNIV          : Universal gas constant (J/mol K)
      !
      !=================================================================

      REAL*8,  PARAMETER            :: AVOGADRO=6.022d23
      REAL*8, PARAMETER             :: R_UNIV = 8.31448d0

!
! PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! Arrays
      !
      ! GE_INJECT          : Injection mass (kg/box) 
      ! GE_TOFFSET         : Temperature offset (K)
      ! GE_LMIN            : Minimum geoengineering height (level)
      ! GE_LMAX            : Maximum geoengineering height (level)
      ! 
      !=================================================================

      ! Arrays
      REAL*8,DIMENSION(:,:,:),ALLOCATABLE   :: GE_INJECT
      REAL*8,DIMENSION(:,:,:),ALLOCATABLE   :: GE_TOFFSET
      INTEGER,DIMENSION(:,:),ALLOCATABLE    :: GE_LMIN, GE_LMAX

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emit_geoeng
!
! !DESCRIPTION: Subroutine EMIT_GEOENG is used to model sulfate geoengineering
!  using SO2.           
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMIT_GEOENG
!
! !USES:
!
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTSO2
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE DAO_MOD,      ONLY : BXHEIGHT
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE GRID_MOD,     ONLY : GET_YEDGE, GET_AREA_M2
!
! !REMARKS:
! (1) A remark
! 
! !REVISION HISTORY: 
!  28 Jul 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8                   :: KG_PER_STEP
      REAL*8                   :: DTEMIS
      REAL*8                   :: CELL_AREA, CELL_VOL
      REAL*8                   :: CELL_MINLAT, CELL_MAXLAT
      REAL*8                   :: CELL_PCEIL, CELL_PBASE
      REAL*8                   :: CELL_ALTFRAC, CELL_LATFRAC
      REAL*8                   :: GE_GLOBALVOL
      LOGICAL, SAVE            :: FIRST=.TRUE.
      INTEGER, SAVE            :: JMIN, JMAX
      INTEGER                  :: I, J, L
      CHARACTER(LEN=255)       :: OUTMSG

      !=================================================================
      ! EMIT_GEOENG begins here!
      !=================================================================

      ! Timestep length in days
      DTEMIS = GET_TS_EMIS() / (60d0 * 24d0)
      IF (FIRST) THEN
         IF (LPRT) THEN
            CALL DEBUG_MSG( '### GEOENG: Initializing geoeng module')
         ENDIF
         ! Zero out injection mass and limits
         GE_INJECT = 0d0
         GE_LMIN = LLPAR
         GE_LMAX = 1

         ! Calculate how much we need to emit per timestep
         ! Assume that we want to emit an average of X TgS per 365.25 days
         KG_PER_STEP = 1.d9 * GE_TGS * DTEMIS / 365.25d0
         ! Convert from TgS to TgSO2
         KG_PER_STEP = KG_PER_STEP * 64d0 / 32d0
         ! Calculate the total volume into which we are injecting
         DO J=1,JJPAR
         CELL_MINLAT = GET_YEDGE(1,J,1)
         CELL_MAXLAT = GET_YEDGE(1,J+1,1)
         IF (CELL_MAXLAT<GE_MINLAT) THEN
            JMIN = MAX(J,JMIN)
            CYCLE
         ENDIF
         IF (CELL_MINLAT>GE_MAXLAT) THEN
            EXIT
         ENDIF
         JMAX = J
         ! Calculate fraction of cell within limits
         CELL_LATFRAC = MIN(CELL_MAXLAT,GE_MAXLAT) - 
     &                       MAX(CELL_MINLAT,GE_MINLAT)
         CELL_LATFRAC = CELL_LATFRAC/(CELL_MAXLAT-CELL_MINLAT)
         CELL_AREA    = GET_AREA_M2(1,J,1)
         DO I=1,IIPAR
         DO L=1,LLPAR
         CELL_PBASE = GET_PEDGE(I,J,L)
         CELL_PCEIL = GET_PEDGE(I,J+1,L)
         IF (CELL_PCEIL.gt.GE_PBASE) THEN
            GE_LMIN(I,J) = L
            CYCLE
         ENDIF
         IF (CELL_PBASE.lt.GE_PCEIL) THEN
            GE_LMAX(I,J) = L
            EXIT
         ENDIF
         ! As for latitude, calculate limits
         ! Are we using the entire cell?
         IF ((CELL_PBASE.gt.GE_PBASE).or.(CELL_PCEIL.lt.GE_PCEIL)) THEN
            CELL_ALTFRAC = MIN(CELL_PBASE,GE_PBASE)/
     &                     MAX(CELL_PCEIL,GE_PCEIL)
            CELL_ALTFRAC = LOG(CELL_ALTFRAC)/LOG(CELL_PBASE/CELL_PCEIL)
         ELSE
            CELL_ALTFRAC = 1.d0
         ENDIF
         CELL_VOL = CELL_ALTFRAC*CELL_LATFRAC*CELL_AREA*BXHEIGHT(I,J,L)
         GE_GLOBALVOL = GE_GLOBALVOL+CELL_VOL
         GE_INJECT(I,J,L) = CELL_VOL
         ENDDO
         ENDDO
         ENDDO
         ! Divide cell volumes by global volume to get fraction
         ! For each target cell, calculate the volume
         DO J=JMIN,JMAX
         DO I=1,IIPAR
         DO L=GE_LMIN(I,J),GE_LMAX(I,J)
         GE_INJECT(I,J,L) = KG_PER_STEP*GE_INJECT(I,J,L)/GE_GLOBALVOL
         ENDDO
         ENDDO
         ENDDO
         FIRST = .FALSE.
         IF (LPRT) THEN
            WRITE(OUTMSG,'(a,a,1x,E10.4)') 'GE_INJECT grid calculated.',
     &          ' Annual mass: ', SUM(GE_INJECT)*365.25/DTEMIS
            CALL DEBUG_MSG(TRIM(OUTMSG))
         ENDIF
      ENDIF

      ! Simply apply the precalculated matrix
      ! SDE DEBUG
      GE_GLOBALVOL = SUM(SUM(SUM(STT(:,:,:,IDTSO2),3),2),1)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
!$OMP+SCHEDULE( DYNAMIC )
      DO J=JMIN,JMAX
      DO I=1,IIPAR
      DO L=GE_LMIN(I,J),GE_LMAX(I,J)
      STT(I,J,L,IDTSO2) = STT(I,J,L,IDTSO2) + GE_INJECT(I,J,L)
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      IF (LPRT) THEN
          WRITE(OUTMSG,'(a,1x,E10.4,1x,a)') 'Geoengineering emission: ',
     &   (0.5d0*(SUM(SUM(SUM(STT(:,:,:,IDTSO2),3),2),1)-GE_GLOBALVOL)),
     &          ' TgS'
          CALL DEBUG_MSG(TRIM(OUTMSG))
      ENDIF

      END SUBROUTINE EMIT_GEOENG
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: offset_met_fields
!
! !DESCRIPTION: Subroutine OFFSET\_MET\_FIELDS applies offsets to met fields
!  to crudely simulate the effects of climate change when appropriate met
!  forecasts are not available
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OFFSET_MET_FIELDS(FIELDTYPE)
!
! !USES:
!
      USE DAO_MOD,     ONLY : T, TMPU1, TMPU2, TSKIN, TS, MOISTQ
      USE PRESSURE_MOD,ONLY : GET_PCENTER
      USE LOGICAL_MOD, ONLY : LCC_AR5T, LCC_SFCT, LCC_UNIFORMT
      USE TIME_MOD,    ONLY : GET_ELAPSED_MIN
      USE ERROR_MOD,   ONLY : ALLOC_ERR
!
! !INPUT PARAMETERS:
!
      INTEGER,INTENT(IN)          :: FIELDTYPE  ! Type of field (A3, A6, I6, A1, I3)
!
! !REVISION HISTORY: 
!  15 Oct 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER,SAVE                :: LASTMIN=-1
      INTEGER                     :: I, J, L, AS
      REAL*8,PARAMETER            :: PMAX = 50.d0
      REAL*8,PARAMETER            :: PMIN = 300.d0
      LOGICAL,SAVE                :: FIRST=.TRUE.
      REAL*8, SAVE                :: INTERPFAC, TMAX, TMIN, LOGOFFSET
      REAL*8                      :: TLOCAL, PLOCAL

      !=================================================================
      ! OFFSET_MET_FIELDS begins here!
      !=================================================================

      ! Do we need to apply any offsets?
      IF (.not.(LCC_AR5T.or.LCC_SFCT.or.LCC_UNIFORMT)) THEN
          FIRST = .FALSE.
          RETURN
      ENDIF

      ! Do we need to recalculate GE_TOFFSET?
      IF (FIRST.or.(LCC_AR5T.and.(GET_ELAPSED_MIN().gt.LASTMIN))) THEN
         ! Allocate temperature offset array
         IF (FIRST) THEN
            ALLOCATE( GE_TOFFSET( IIPAR, JJPAR, LLPAR ), STAT=AS )
            IF ( AS /= 0 ) CALL ALLOC_ERR( 'GE_TOFFSET' )
         ENDIF
         GE_TOFFSET = 0d0
         IF (LCC_AR5T) THEN
            LASTMIN = GET_ELAPSED_MIN()
            IF (FIRST) THEN
               ! Calculate re-used parameters
               TMAX      = GE_TSURF
               TMIN      = GE_TSURF*(-2.5d0)
               INTERPFAC = (TMAX-TMIN)/LOG(PMIN/PMAX)
               LOGOFFSET = LOG(PMIN)
            ENDIF
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, PLOCAL, TLOCAL )
!$OMP+SCHEDULE( DYNAMIC )
            DO I=1,IIPAR
            DO J=1,JJPAR
            DO L=1,LLPAR
               PLOCAL = GET_PCENTER(I,J,L)
               IF (PLOCAL.gt.PMIN) THEN
                  ! Below 300 hPa
                  TLOCAL = TMAX
               ELSEIF (PLOCAL.lt.PMAX) THEN
                  ! Above 50 hPa
                  TLOCAL = TMIN
               ELSE
                  ! Between the two - altitude-linear scaling
                  TLOCAL = TMAX + (INTERPFAC*(LOG(PLOCAL)-LOGOFFSET))
               ENDIF
               GE_TOFFSET(I,J,L) = TLOCAL
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ELSEIF (LCC_SFCT) THEN
            GE_TOFFSET(:,:,1) = GE_TSURF
         ELSEIF (LCC_UNIFORMT) THEN
            GE_TOFFSET(:,:,:) = GE_TSURF
         ENDIF
      ENDIF

      ! Run through available field types
      SELECT CASE ( FIELDTYPE )
         CASE ( 1 )
            ! A3
#if defined( GEOS_4 ) || defined( GEOS_5 ) || defined( GCAP )
            TS(:,:) = TS(:,:) + GE_TOFFSET(:,:,1)
#if !defined( GCAP )
            TSKIN(:,:) = TSKIN(:,:) + GE_TOFFSET(:,:,1)
#elif defined( MERRA )
            T(:,:,:) = T(:,:,:) + GE_TOFFSET(:,:,:)
#endif
#endif
         CASE ( 2 )
#if defined( GEOS_4 ) || defined( GEOS_5 ) || defined( GCAP )
            ! A6
            T(:,:,:) = T(:,:,:) + GE_TOFFSET(:,:,:)
#endif
         CASE ( 3 )
            ! I6
         CASE ( 4 )
            ! A1
#if defined( MERRA ) || defined( GEOS_57 )
            TS(:,:) = TS(:,:) + GE_TOFFSET(:,:,1)
            TSKIN(:,:) = TSKIN(:,:) + GE_TOFFSET(:,:,1)
#endif
         CASE ( 5 )
            ! I3
#if defined( GEOS_57 )
            TMPU1(:,:,:) = TMPU1(:,:,:) + GE_TOFFSET(:,:,:)
            TMPU2(:,:,:) = TMPU2(:,:,:) + GE_TOFFSET(:,:,:)
#endif
         CASE DEFAULT
            ! Do nothing
      END SELECT

      FIRST = .FALSE.

      END SUBROUTINE OFFSET_MET_FIELDS
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_geoeng
!
! !DESCRIPTION: Subroutine INIT\_GEOENG initializes module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GEOENG( am_I_Root )
!
! !USES:
!
      USE TRACERID_MOD
      USE LOGICAL_MOD
      USE TIME_MOD,     ONLY : GET_YEAR
      USE ERROR_MOD,    ONLY : ALLOC_ERR
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root  ! Is this the root CPU?
!
! !OUTPUT VARIABLES:
!
!
! !REVISION HISTORY: 
!  11 Nov 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: N, AS
      CHARACTER(LEN=255) :: DBGMSG

      !=================================================================
      ! INIT_GEOENG begins here!
      !=================================================================

      WRITE( 6,'(a)') REPEAT( '=', 79 )
      WRITE( 6,'(a)') 'G E O E N G I N E E R I N G'
      WRITE( 6,'(a)') 'Routines written by SEBASTIAN D. EASTHAM'
      WRITE( 6,'(a)') 'Last Modification Date: 11/11/13'
      WRITE( 6,'(a)') REPEAT( '=', 79 )

      ! Allocate geoengineering injection array
      ALLOCATE( GE_INJECT( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GE_INJECT' )
      GE_INJECT = 0d0

      ! Allocate arrays of altitude limits
      ALLOCATE( GE_LMIN( IIPAR,JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GE_LMIN' )
      GE_LMIN = LLPAR

      ALLOCATE( GE_LMAX( IIPAR,JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GE_LMAX' )
      GE_LMAX = 1
      
      END SUBROUTINE INIT_GEOENG
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_geoeng
!
! !DESCRIPTION: Subroutine CLEANUP\_GEOENG deallocates module variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_GEOENG
!
! !REVISION HISTORY: 
!  11 Nov 2013 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      ! CLEANUP_GEOENGAT begins here!
      !=================================================================

      IF ( ALLOCATED( GE_INJECT  ) ) DEALLOCATE( GE_INJECT  )
      IF ( ALLOCATED( GE_TOFFSET ) ) DEALLOCATE( GE_TOFFSET )
      IF ( ALLOCATED( GE_LMIN    ) ) DEALLOCATE( GE_LMIN    )
      IF ( ALLOCATED( GE_LMAX    ) ) DEALLOCATE( GE_LMAX    )
      
      END SUBROUTINE CLEANUP_GEOENG
!
!EOC
      END MODULE GEOENG_MOD
