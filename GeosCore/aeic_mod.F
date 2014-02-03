! $Id: aeic_mod.f,v 1.0 2012/12/12 13:34:07 sde Exp $
      MODULE AEIC_MOD
!
!******************************************************************************
!  Module AEIC_MOD contains variables and routines for aircraft flight
!  emissions into the chemistry and transport grids. (sde, 12/14/12)
!
!  The AEIC fields are stored on 1x1 grids.
!  These fields will be interpolated onto the current horizontal grid.
!
!  Module Variables:
!  ============================================================================
!  (1 ) AEIC_EMIT (REAL*8 )  : Array for emissions on the GEOS-CHEM grid
!  (2 ) OPTIONS_AEIC         : Structure containing options
!
!  Module Routines:
!  ============================================================================
!  (1 ) READ_AEIC_FILE       : Routine to read emissions from disk
!  (2 ) EMIT_AEIC            : Routine to emit pollutants into GEOS-CHEM
!  (3 ) INIT_AEIC            : Routine to allocate/initialize module variables
!  (4 ) CLEANUP_AEIC         : Routine to deallocate module variables
! 
!  GEOS-CHEM modules referenced by biomass_mod.f
!  ============================================================================
!  (1 ) diag_mod.F     : Module containing GEOS-CHEM diagnostic arrays
!  (2 ) error_mod.F    : Module containing NaN and other error check routines
!  (3 ) file_mod.F     : Module containing file unit numbers and error checks
!  (4 ) grid_mod.F     : Module containing horizontal grid information
!  (5 ) pressure_mod.F : Module containing routines to compute P(I,J,L)
!  (6 ) tropopause_mod.F:Module containing routines to find tropopause
! 
!  NOTES:
!  14 Dec 2012 - S. Eastham - Adapted from AIRCRAFT_NOX_MOD. Now emits
!                             aircraft NOx, CO, HC, SO2, SO4, OC and BC
!******************************************************************************
!

      IMPLICIT NONE
#     include "define.h"
#     include "netcdf.inc" ! Needed for netCDF libraries

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "aeic_mod.f"
      !=================================================================

      ! Initially declare all private
      PRIVATE

      ! PUBLIC module routines
      PUBLIC EMIT_AEIC
      PUBLIC CLEANUP_AEIC
      PUBLIC OPTIONS_AEIC

      ! PUBLIC options structure
      TYPE, PUBLIC :: Aircraft_Struct
         LOGICAL   :: LLTO_EMIS
         REAL*8    :: LTO_FB_MULT
         REAL*8    :: LTO_CO_MULT
         REAL*8    :: LTO_NOx_MULT
         REAL*8    :: LTO_HC_MULT
         REAL*8    :: LTO_BC_EI
         REAL*8    :: LTO_OC_EI
         LOGICAL   :: LCRUISE_EMIS
         REAL*8    :: CRUISE_FB_MULT
         REAL*8    :: CRUISE_CO_MULT
         REAL*8    :: CRUISE_NOx_MULT
         REAL*8    :: CRUISE_HC_MULT
         REAL*8    :: CRUISE_BC_EI
         REAL*8    :: CRUISE_OC_EI
         REAL*8    :: FSC
         REAL*8    :: SULFCONV
         INTEGER   :: EXC_CODE
         INTEGER   :: EXC_UR_I
         INTEGER   :: EXC_UR_J
         INTEGER   :: EXC_LL_I
         INTEGER   :: EXC_LL_J
         REAL*8    :: ALT_FLR    ! Altitude cutoff (only emit above)
         REAL*8    :: ALT_OFFSET ! Altitude offset
      END TYPE Aircraft_Struct

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      
      ! Emissions on native 1x1 grid
      REAL*8,  ALLOCATABLE :: AEIC_1x1(:,:,:,:)
      ! Regridded emissions
      REAL*8,  ALLOCATABLE :: AEIC_REGRID(:,:,:,:)
      ! Options structure
      TYPE(Aircraft_Struct) :: OPTIONS_AEIC
 
      ! Number of levels on the 1x1 grid
      ! Standard resolution (~3 GB of data)
      INTEGER, PARAMETER :: NLVL_HD = 5
      REAL*8, PARAMETER  :: DEPTHLVL_HD = 60.96d0
      INTEGER, PARAMETER :: NLVL_SD = 44
      REAL*8, PARAMETER  :: DEPTHLVL_SD = 304.8d0
      INTEGER, PARAMETER :: NLVLAEIC = 49
 
      ! Aircraft emissions on GEOS-CHEM grid
      REAL*8,  ALLOCATABLE :: AEIC_EMIT(:,:,:,:)

      INTEGER, PARAMETER :: FB_INDEX = 1
      INTEGER, PARAMETER :: HC_INDEX = 2
      INTEGER, PARAMETER :: CO_INDEX = 3
      INTEGER, PARAMETER :: NOx_INDEX = 4
      INTEGER, PARAMETER :: BC_INDEX = 5
      INTEGER, PARAMETER :: OC_INDEX = 6
  
      ! THC to TOG conversion factor for aircraft emissions
      REAL, PARAMETER   :: THC2TOG      = 1.16d0
      ! HC speciation is included below
      REAL, PARAMETER   :: ACF_ACET     = 0.003693477d0
      REAL, PARAMETER   :: ACF_ALD2     = 0.042718224d0
      REAL, PARAMETER   :: ACF_ALK4     = 0.213791063d0
      REAL, PARAMETER   :: ACF_C2H6     = 0.005214505d0
      REAL, PARAMETER   :: ACF_C3H8     = 0.000780871d0
      REAL, PARAMETER   :: ACF_CH2O     = 0.123081099d0
      REAL, PARAMETER   :: ACF_PRPE     = 0.178041756d0
      REAL, PARAMETER   :: ACF_MACR     = 0.005362609d0
      REAL, PARAMETER   :: ACF_RCHO     = 0.036769436d0
      ! Note not all aircraft hydrocarbon species modeled in GEOS-Chem

      ! LTO ceiling level (3000 ft = 914.4 m)
      REAL*8            :: LTOCEIL = 9.144d2

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_aeic_file
!
! !DESCRIPTION: Subroutine READ\_AEIC\_FILE reads fuel burn and
!  emissions from AEIC netCDF files (monthly average). (sde, 12/14/12)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_AEIC_FILE( INMONTH )
!
! !USES:
!
      ! NcdfUtil modules for netCDF I/O
      USE m_netcdf_io_open                    ! netCDF open
      USE m_netcdf_io_get_dimlen              ! netCDF dimension queries
      USE m_netcdf_io_read                    ! netCDF data reads
      USE m_netcdf_io_close                   ! netCDF close

      USE DAO_MOD,        ONLY : T
      USE REGRID_A2A_MOD, ONLY : DO_REGRID_A2A
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE TIME_MOD,       ONLY : EXPAND_DATE
      USE ERROR_MOD,      ONLY : ERROR_STOP
      USE PRESSURE_MOD,   ONLY : GET_AP, GET_BP

      USE CMN_SIZE_MOD         ! Size parameters
      USE CMN_GCTM_MOD, ONLY : Rdg0

!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)     :: INMONTH ! Current month number (1-12)
!
! !REVISION HISTORY: 
!  14 Dec 2012 - S. Eastham - Adapted from READ_NONERUP_VOLC
!  01 Aug 2013 - M. Sulprizio- Update to read emissions from netCDF file
!                              created by Christoph Keller. He combined original
!                              monthly files into one file and made file
!                              COARDS compliant.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8, TARGET       :: ACARRAY_IN(360,180,4)
      REAL*8, TARGET       :: ACARRAY_OUT(IIPAR,JJPAR)
      INTEGER              :: I,J,L, AEIC_FID, EMIT_INDEX
      CHARACTER(LEN=255)   :: AEIC_FILE,MONTH_STR
      CHARACTER(LEN=255)   :: LLFILENAME, EMIT_NAME
      CHARACTER(LEN=255)   :: AEIC_GRIDSTR
      REAL*8, POINTER      :: OUTGRID(:,:) => NULL()
      REAL*8, POINTER      :: INGRID(:,:)  => NULL()
      REAL*8               :: LTOFRAC, PBASE, PCEIL
      REAL*8               :: CURRBASE, CURRCEIL
      LOGICAL              :: INCRUISE
      REAL*8,PARAMETER     :: PSEALEV=1013.25d0
      REAL*8,DIMENSION(6)  :: LTOSCALE,CRUISESCALE,EMITSCALE
      REAL*8               :: TEMP ! Local temperature
      INTEGER              :: L_AEIC
      REAL*8               :: FRAC_AEIC,AEIC_CEIL,AEIC_BASE,ACTIVEDEPTH
      REAL*8               :: MINALT,PLUSALT,DEPTHLVLAEIC
      LOGICAL              :: RESCALE_GRID,AEIC_DONE

      !=================================================================
      ! READ_AEIC begins here!
      !=================================================================

      ! Retrieve values from options structure and convert from km to m
      PLUSALT = OPTIONS_AEIC%ALT_OFFSET * 1.0d3
      MINALT = (OPTIONS_AEIC%ALT_FLR * 1.0d3) + PLUSALT

      ! Filename
      AEIC_FILE = TRIM( DATA_DIR_1x1 ) // '/AEIC_201301/' //
     &           'AEIC.nc'

      ! Echo output
      WRITE( 6, 100 ) TRIM( AEIC_FILE )
 100  FORMAT( '     - READ_AEIC: Reading ', a )

      ! Read 3D NetCDF file containing emissions data for this month
      ! Open file
      CALL NcOp_Rd( AEIC_FID, TRIM(AEIC_FILE))

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1 ) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_generic1x1.nc'

      DO EMIT_INDEX = 1, 4
  
         SELECT CASE (EMIT_INDEX)
            CASE (FB_INDEX)
               EMIT_NAME = 'fuelBurn'
            CASE (CO_INDEX)
               EMIT_NAME = 'CO'
            CASE (HC_INDEX)
               EMIT_NAME = 'HC'
            CASE (NOx_INDEX)
               EMIT_NAME = 'NOx'
            CASE DEFAULT
               CALL ERROR_STOP( 'Unrecognized aircraft emission!',
     &               'READ_AEIC (aeic_mod.F90)' )
         END SELECT
  
         CALL NcRd(AEIC_1x1(:,:,:,EMIT_INDEX),AEIC_FID,TRIM(EMIT_NAME),
     &       (/ 1, 1, 1 , INMONTH /),(/360, 180, NLVLAEIC, 1 /))
 
         ! Fuel burn sum
         IF (EMIT_INDEX .eq. FB_INDEX) THEN
            WRITE(6,*) 'Fuelburn start sum: ', SUM(AEIC_1x1), 'kg'
         ENDIF
      ENDDO

      ! Close NetCDF file
      CALL NcCl(AEIC_FID)
 
      ! Split out LTO and cruise emissions
      LTOSCALE(FB_INDEX) = OPTIONS_AEIC%LTO_FB_MULT
      LTOSCALE(HC_INDEX) = OPTIONS_AEIC%LTO_FB_MULT*
     &      OPTIONS_AEIC%LTO_HC_MULT*1.0d-3
      LTOSCALE(NOx_INDEX)= OPTIONS_AEIC%LTO_FB_MULT*
     &      OPTIONS_AEIC%LTO_NOX_MULT*1.0d-3
      LTOSCALE(CO_INDEX) = OPTIONS_AEIC%LTO_FB_MULT*
     &      OPTIONS_AEIC%LTO_CO_MULT*1.0d-3
      LTOSCALE(BC_INDEX) = OPTIONS_AEIC%LTO_FB_MULT*
     &      OPTIONS_AEIC%LTO_BC_EI*1.0d-3
      LTOSCALE(OC_INDEX) = OPTIONS_AEIC%LTO_FB_MULT*
     &      OPTIONS_AEIC%LTO_OC_EI*1.0d-3
      CRUISESCALE(FB_INDEX) = OPTIONS_AEIC%CRUISE_FB_MULT
      CRUISESCALE(HC_INDEX) = OPTIONS_AEIC%CRUISE_FB_MULT*
     &      OPTIONS_AEIC%CRUISE_HC_MULT*1.0d-3
      CRUISESCALE(NOx_INDEX)= OPTIONS_AEIC%CRUISE_FB_MULT*
     &      OPTIONS_AEIC%CRUISE_NOX_MULT*1.0d-3
      CRUISESCALE(CO_INDEX) = OPTIONS_AEIC%CRUISE_FB_MULT*
     &      OPTIONS_AEIC%CRUISE_CO_MULT*1.0d-3
      CRUISESCALE(BC_INDEX) = OPTIONS_AEIC%CRUISE_FB_MULT*
     &      OPTIONS_AEIC%CRUISE_BC_EI*1.0d-3
      CRUISESCALE(OC_INDEX) = OPTIONS_AEIC%CRUISE_FB_MULT*
     &      OPTIONS_AEIC%CRUISE_OC_EI*1.0d-3

      ! First step: rescale grid horizontally at each level
      DO L=1,NLVLAEIC
         ACARRAY_IN = AEIC_1x1(:,:,L,:)
         DO EMIT_INDEX=1,4
            ! Point to input & output arrays
            INGRID  => ACARRAY_IN(:,:,EMIT_INDEX) !Units kg/s
            OUTGRID => ACARRAY_OUT

            ! Regrid
            CALL DO_REGRID_A2A( LLFILENAME, 360,    180,
     &         INGRID,     OUTGRID, PERAREA=1,
     &         netCDF=.TRUE.                   )

            AEIC_REGRID(:,:,L,EMIT_INDEX) = OUTGRID

            ! Free pointers
            NULLIFY( INGRID, OUTGRID )
         ENDDO
      ENDDO

      ! Re-zero output grid
      AEIC_EMIT = 0d0

      DO I=1,IGLOB
         DO J=1,JGLOB
            CURRCEIL = 0d0
            PCEIL = PSEALEV
            INCRUISE = .FALSE.
            EMITSCALE = LTOSCALE
            RESCALE_GRID = .TRUE.
            L_AEIC = 1
            AEIC_CEIL = PLUSALT
            ACARRAY_IN = 0.0d0
            AEIC_DONE = .FALSE.
            DO L=1,LLPAR
               ! Are we in LTO or cruise (or fraction?)
               CURRBASE = CURRCEIL
               PBASE = PCEIL

               ! Calculate ceiling of next box on standard grid
               PCEIL = GET_AP(L+1) + (GET_BP(L+1) * PSEALEV)

               ! Temperature [K]
               TEMP    = T(I,J,L)

               CURRCEIL = CURRBASE + (Rdg0*TEMP*LOG(PBASE/PCEIL))

               IF ((.not. INCRUISE).and.(LTOCEIL .lt. CURRCEIL)) THEN
                  ! Transition box
                  INCRUISE = .TRUE.
                  LTOFRAC = (LTOCEIL-CURRBASE)/(CURRCEIL-CURRBASE)
                  EMITSCALE = LTOSCALE*LTOFRAC
                  EMITSCALE = EMITSCALE+(CRUISESCALE*(1.0d0-LTOFRAC))
               ELSEIF (INCRUISE) THEN
                  EMITSCALE = CRUISESCALE
               ELSE
                  EMITSCALE = LTOSCALE
               ENDIF
 
               ! Scan through layers until we hit a GC boundary
               AEIC_DONE = ((AEIC_CEIL.le.CURRCEIL) .and. 
     &              (L_AEIC.gt.NLVLAEIC))
               ACARRAY_IN = 0d0
               DO WHILE (.not. AEIC_DONE)
                  ! Calculate fraction of AEIC level within GC level
                  IF (L_AEIC .le. NLVL_HD) THEN
                     ! In finely resolved area
                     DEPTHLVLAEIC = DEPTHLVL_HD
                     AEIC_CEIL = PLUSALT + (REAL(L_AEIC)*DEPTHLVLAEIC)
                  ELSE
                     ! Standard resolution region
                     DEPTHLVLAEIC = DEPTHLVL_SD
                     AEIC_CEIL = PLUSALT + (REAL(NLVL_HD)*DEPTHLVL_HD) +
     &                           (REAL(L_AEIC-NLVL_HD)*DEPTHLVL_SD)
                  ENDIF
                  AEIC_BASE = AEIC_CEIL-DEPTHLVLAEIC
                  IF (AEIC_CEIL .gt. MINALT) THEN
                     ! Establish base level for input cell
                     IF (AEIC_BASE.gt. MINALT) THEN                
                        FRAC_AEIC = MAX(CURRBASE,AEIC_BASE)
                        ACTIVEDEPTH = DEPTHLVLAEIC
                     ELSE
                        FRAC_AEIC = MAX(CURRBASE,AEIC_BASE,MINALT)
                        ACTIVEDEPTH=AEIC_CEIL-MAX(MINALT,AEIC_BASE)
                     ENDIF
                  
                     IF (AEIC_CEIL.lt.CURRCEIL) THEN
                        ! Switch input level
                        FRAC_AEIC = (AEIC_CEIL-FRAC_AEIC)/ACTIVEDEPTH
                     ELSE
                        ! Switch output level
                        FRAC_AEIC = (CURRCEIL-FRAC_AEIC)/ACTIVEDEPTH
                     ENDIF
                  ELSE
                     ! Below cutoff
                     FRAC_AEIC = 0.0d0
                  ENDIF

                  ! Add to output grid
                  DO EMIT_INDEX=1,4
                     AEIC_EMIT(I,J,L,EMIT_INDEX) = 
     &                  AEIC_EMIT(I,J,L,EMIT_INDEX) + 
     &                  (AEIC_REGRID(I,J,L_AEIC,EMIT_INDEX)*FRAC_AEIC)
                  ENDDO

                  IF (AEIC_CEIL.lt.CURRCEIL) THEN
                     L_AEIC = L_AEIC + 1
                     IF (L_AEIC.gt.NLVLAEIC) THEN
                        AEIC_DONE = .TRUE.
                     ENDIF
                  ELSE
                     ! Change GC level now
                     AEIC_DONE = .TRUE.
                  ENDIF
               ENDDO

               ! Copy OC/BC from fuelburn, apply EIs in cruise/LTO
               AEIC_EMIT(I,J,L,OC_INDEX) = AEIC_EMIT(I,J,L,FB_INDEX)
               AEIC_EMIT(I,J,L,BC_INDEX) = AEIC_EMIT(I,J,L,FB_INDEX)

               ! Store emissions
               AEIC_EMIT(I,J,L,:) = EMITSCALE*AEIC_EMIT(I,J,L,:)
            ENDDO
         ENDDO
      ENDDO
 
      ! Apply exclusion zones
      IF (OPTIONS_AEIC%EXC_CODE .eq. 0) THEN
         ! Do nothing
      ELSE IF ( OPTIONS_AEIC%EXC_CODE .eq. 1 ) THEN
         ! Only use one zone
         AEIC_EMIT(1:(OPTIONS_AEIC%EXC_LL_I-1),:,:,:) = 0.0d0
         AEIC_EMIT((OPTIONS_AEIC%EXC_UR_I+1):IIPAR,:,:,:) = 0.0d0
         AEIC_EMIT(:,1:(OPTIONS_AEIC%EXC_LL_J-1),:,:) = 0.0d0
         AEIC_EMIT(:,(OPTIONS_AEIC%EXC_UR_J+1):JJPAR,:,:) = 0.0d0
      ELSE IF (OPTIONS_AEIC%EXC_CODE .eq. -1 ) THEN
         ! Exclude one zone
         AEIC_EMIT(OPTIONS_AEIC%EXC_LL_I:OPTIONS_AEIC%EXC_UR_I,
     &         OPTIONS_AEIC%EXC_LL_J:OPTIONS_AEIC%EXC_UR_J,:,:) = 0.0d0
      ELSE
         ! Invalid
         CALL ERROR_STOP( 
     &      'Aircraft emissions exlusion flag invalid!',
     &      'READ_AEIC_FILE ("aeic_mod.F")' )
      END IF
      END SUBROUTINE READ_AEIC_FILE
!EOC
!------------------------------------------------------------------------------

      SUBROUTINE EMIT_AEIC
!
!******************************************************************************
!  Subroutine EMIT_AEIC interpolates AEIC aircraft emissions from the
!  native grid onto the given GEOS-CHEM grid. (sde, 12/14/12)
!
!  NOTES: 
!  (1 ) 
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD,     ONLY : GET_XOFFSET, GET_YOFFSET
      USE TRACER_MOD,   ONLY : STT, XNUMOL
      USE PRESSURE_MOD, ONLY : GET_AP, GET_BP
      USE TRACERID_MOD
      USE CHEMGRID_MOD, ONLY : GET_CHEMGRID_LEVEL
      USE TIME_MOD,     ONLY : GET_MONTH,GET_TS_EMIS
      USE DAO_MOD,      ONLY : T
      USE LOGICAL_MOD,  ONLY : LUCX

      USE CMN_SIZE_MOD  ! Size parameters
      USE CMN_DIAG_MOD  ! Diagnostic switches
      USE CMN_MOD       ! PTOP, SIGE, AVP

      INTEGER             :: I,    J,     IREF, JREF,  L
      INTEGER             :: I0,   J0, MAXLEV
      REAL*8              :: DTSRCE,TMPMULT,TEMP
      INTEGER, SAVE       :: LAST_MONTH = 0
      LOGICAL, SAVE       :: FIRST = .TRUE.
      REAL*8, SAVE        :: FBSUM = 0d0
      
      !=================================================================
      ! EMIT_AEIC begins here!
      !=================================================================

      IF ( FIRST ) THEN
         FIRST = .FALSE.

         ! Echo info
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, 100   )
         WRITE( 6, 110   )
         WRITE( 6, 120   )
         WRITE( 6, 130   )
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

         ! FORMAT strings
 100     FORMAT( 'A I R C R A F T   E M I S S I O N S'   )
 110     FORMAT( 'Routines originally by STEVEN BARRETT'    )
 120     FORMAT( 'Edited by Jamin Koo and S. D. Eastham'   )
 130     FORMAT( 'Last Modification Date: 12/17/12'       )

         CALL INIT_AEIC

         WRITE(6,*) 'Aircraft FSC (ppm) : ',OPTIONS_AEIC%FSC*1.d6
         WRITE(6,*) 'S(VI) -> S(IV) %   : ',OPTIONS_AEIC%SULFCONV*1.d2
         WRITE(6,*) 'LTO included?      : ',OPTIONS_AEIC%LLTO_EMIS
         IF (OPTIONS_AEIC%LLTO_EMIS) THEN
         WRITE(6,*) ' - Fuelburn factor : ',OPTIONS_AEIC%LTO_FB_MULT
         WRITE(6,*) ' - CO factor       : ',OPTIONS_AEIC%LTO_CO_MULT
         WRITE(6,*) ' - NOx factor      : ',OPTIONS_AEIC%LTO_NOx_MULT
         WRITE(6,*) ' - HC factor       : ',OPTIONS_AEIC%LTO_HC_MULT
         WRITE(6,*) ' - BC EI           : ',OPTIONS_AEIC%LTO_BC_EI
         WRITE(6,*) ' - OC EI           : ',OPTIONS_AEIC%LTO_BC_EI
         END IF
         WRITE(6,*) 'Cruise included?   : ',OPTIONS_AEIC%LCRUISE_EMIS
         IF (OPTIONS_AEIC%LCRUISE_EMIS) THEN
         WRITE(6,*) ' - Fuelburn factor : ',OPTIONS_AEIC%CRUISE_FB_MULT
         WRITE(6,*) ' - CO factor       : ',OPTIONS_AEIC%CRUISE_CO_MULT
         WRITE(6,*) ' - NOx factor      : ',OPTIONS_AEIC%CRUISE_NOx_MULT
         WRITE(6,*) ' - HC factor       : ',OPTIONS_AEIC%CRUISE_HC_MULT
         WRITE(6,*) ' - BC EI           : ',OPTIONS_AEIC%CRUISE_BC_EI
         WRITE(6,*) ' - OC EI           : ',OPTIONS_AEIC%CRUISE_BC_EI
         END IF
         WRITE(6,*) 'Offset (km)        : ',OPTIONS_AEIC%ALT_OFFSET
         WRITE(6,*) 'Cutoff (km)        : ',OPTIONS_AEIC%ALT_FLR
#if defined( UCX )
         IF ((LUCX).and.(IDTH2O.ne.0)) THEN
            WRITE(6,*) ' - H2O EI          : ',1260.d0
         ELSE
            WRITE(6,*) 'Aircraft H2O emissions ignored.'
         ENDIF
#endif

         FBSUM = 0d0

      END IF

      ! Read aircraft NOx emissions
      IF (LAST_MONTH/=GET_MONTH()) THEN
         LAST_MONTH = GET_MONTH()
         CALL READ_AEIC_FILE(LAST_MONTH)
         WRITE(6,*) 'Aircraft emissions for month ', LAST_MONTH
         WRITE(6,*) 'Aircraft fuelburn: ',
     &       SUM(AEIC_EMIT(:,:,:,FB_INDEX)), ' kg/s'
         WRITE(6,*) 'Aircraft NOx     : ',
     &       SUM(AEIC_EMIT(:,:,:,NOx_INDEX)), ' kg/s'
         WRITE(6,*) 'Aircraft CO      : ',
     &       SUM(AEIC_EMIT(:,:,:,CO_INDEX)), ' kg/s'
         WRITE(6,*) 'Aircraft BC/PMNV : ',
     &       SUM(AEIC_EMIT(:,:,:,BC_INDEX)), ' kg/s'
         WRITE(6,*) 'Aircraft OC/PMV  : ',
     &       SUM(AEIC_EMIT(:,:,:,OC_INDEX)), ' kg/s'
         WRITE(6,*) 'Aircraft HC      : ',
     &       SUM(AEIC_EMIT(:,:,:,HC_INDEX)), ' kg/s'

         IF (LAST_MONTH .eq. 1) THEN
            FBSUM = 0d0
         ENDIF

      ENDIF

      ! Emissions timestep
      DTSRCE = GET_TS_EMIS() * 60d0

      ! Get nested-grid offsets
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()

!      ! Debug
!      WRITE(6,*) 'PRESO4:  ', SUM(STT(:,:,:,IDTSO4))
!      WRITE(6,*) 'PRENOx:  ', SUM(STT(:,:,:,IDTNOx))

      ! Loop over surface grid boxes
      DO J = 1, JJPAR
         JREF = J + J0
      DO I = 1, IIPAR
         IREF = I + I0
         !==============================================================
         ! Loop over GEOS-CHEM levels 
         !==============================================================
         MAXLEV = GET_CHEMGRID_LEVEL(I,J)
         DO L = 1, LLPAR
             ! First handle fuelburn-multiplied emissinos
             TMPMULT = (OPTIONS_AEIC%FSC*DTSRCE/32.0d0)
             STT(I,J,L,IDTSO4) = STT(I,J,L,IDTSO4) + 
     &          ((AEIC_EMIT(I,J,L,FB_INDEX)*96.0d0*TMPMULT)
     &           *OPTIONS_AEIC%SULFCONV)
             STT(I,J,L,IDTSO2) = STT(I,J,L,IDTSO2) +
     &          ((AEIC_EMIT(I,J,L,FB_INDEX)*64.0d0*TMPMULT)
     &           *(1.d0-OPTIONS_AEIC%SULFCONV))
             STT(I,J,L,IDTNOx) = STT(I,J,L,IDTNOx) +
     &           (AEIC_EMIT(I,J,L,NOx_INDEX)*DTSRCE)
             STT(I,J,L,IDTCO) = STT(I,J,L,IDTCO) +
     &           (AEIC_EMIT(I,J,L,CO_INDEX)*DTSRCE)
             STT(I,J,L,IDTOCPI) = STT(I,J,L,IDTOCPI) +
     &           (AEIC_EMIT(I,J,L,OC_INDEX)*DTSRCE)
             STT(I,J,L,IDTBCPI) = STT(I,J,L,IDTBCPI) +
     &           (AEIC_EMIT(I,J,L,BC_INDEX)*DTSRCE)
#if defined( UCX )
             IF ((LUCX).and.(IDTH2O.ne.0)) THEN
                STT(I,J,L,IDTH2O) = STT(I,J,L,IDTH2O) +
     &              (1.26d0*AEIC_EMIT(I,J,L,FB_INDEX)*DTSRCE)
             ENDIF
#endif
             IF (L.le.MAXLEV) THEN
                ! Only emit HC within chemgrid to avoid buildup
                TMPMULT = THC2TOG*DTSRCE*AEIC_EMIT(I,J,L,HC_INDEX)
                STT(I,J,L,IDTMACR) = STT(I,J,L,IDTMACR) +
     &               (TMPMULT*ACF_MACR)
                STT(I,J,L,IDTRCHO) = STT(I,J,L,IDTRCHO) +
     &               (TMPMULT*ACF_RCHO)
                STT(I,J,L,IDTACET) = STT(I,J,L,IDTACET) +
     &               (TMPMULT*ACF_ACET)
                STT(I,J,L,IDTALD2) = STT(I,J,L,IDTALD2) +
     &               (TMPMULT*ACF_ALD2)
                STT(I,J,L,IDTALK4) = STT(I,J,L,IDTALK4) +
     &               (TMPMULT*ACF_ALK4)
                STT(I,J,L,IDTC2H6) = STT(I,J,L,IDTC2H6) +
     &               (TMPMULT*ACF_C2H6)
                STT(I,J,L,IDTC3H8) = STT(I,J,L,IDTC3H8) +
     &               (TMPMULT*ACF_C3H8)
                STT(I,J,L,IDTCH2O) = STT(I,J,L,IDTCH2O) +
     &               (TMPMULT*ACF_CH2O)
                STT(I,J,L,IDTPRPE) = STT(I,J,L,IDTPRPE) +
     &               (TMPMULT*ACF_PRPE)
             ENDIF
         ENDDO  ! L
      ENDDO     ! I 
      ENDDO     ! J

!      ! Debug
!      WRITE(6,*) 'POSTSO4: ', SUM(STT(:,:,:,IDTSO4))
!      WRITE(6,*) 'POSTNOx: ', SUM(STT(:,:,:,IDTNOx))
!      FBSUM = FBSUM + (DTSRCE*SUM(AEIC_EMIT(:,:,:,FB_INDEX)))
!      WRITE(6,*) 'Fuel burn in current year: ', FBSUM*1.d-9, 'Tg'

      ! Return to calling program
      END SUBROUTINE EMIT_AEIC

!------------------------------------------------------------------------------

      SUBROUTINE INIT_AEIC
!
!******************************************************************************
!  Subroutine INIT_AEIC allocates and initializes module variables.
!  (sde, 12/16/12)
!
!  NOTES:
!  (1 )
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

      USE CMN_SIZE_MOD
      
      ! Local variables
      INTEGER :: AS

      ! AEIC_1x1 is used temporarily for input grids
      ALLOCATE( AEIC_1x1( 360, 180, NLVLAEIC, 4 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEIC_1x1' )
      AEIC_1x1 = 0d0

      ALLOCATE( AEIC_REGRID( IIPAR, JJPAR, NLVLAEIC, 4 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEIC_REGRID' )
      AEIC_REGRID = 0d0
  
      ! AEIC_EMIT holds the output grids
      ALLOCATE( AEIC_EMIT( IIPAR, JJPAR, LLPAR, 6 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AEIC_EMIT' )
      AEIC_EMIT = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_AEIC

!------------------------------------------------------------------------------
!EOC

!BOP
      SUBROUTINE CLEANUP_AEIC
!
!******************************************************************************
!  Subroutine CLEANUP_AEIC deallocates module variables. (sde, 12/16/12)
!
!  NOTES:  
!  (1 ) 
!******************************************************************************
!
      IF ( ALLOCATED( AEIC_1x1    ) ) DEALLOCATE( AEIC_1x1    )
      IF ( ALLOCATED( AEIC_EMIT   ) ) DEALLOCATE( AEIC_EMIT   )

      ! Return to calling program
      END SUBROUTINE CLEANUP_AEIC

!------------------------------------------------------------------------------

      ! End of module
      END MODULE AEIC_MOD
