! $Id: diag48_mod.f,v 1.3 2009/10/15 17:46:24 bmy Exp $
      MODULE DIAG48_MOD
!
!******************************************************************************
!  Module DIAG48_MOD contains variables and routines to save out 3-D 
!  timeseries output to disk (bmy, 7/20/04, 10/13/09)
!
!  Module Variables:
!  ============================================================================
!  (2 ) I0             (INTEGER ) : Lon offset between global & nested grid
!
!  Module Routines:
!  ============================================================================
!  (1 ) DIAG48                 : Main driver routine
!  (2 ) ITS_TIME_FOR_DIAG48    : Returns TRUE if it's time to save to disk
!  (3 ) INIT_DIAG48            : Gets variable values from "input_mod.f"
!
!  GEOS-CHEM modules referenced by "diag48_mod.f" 
!  ============================================================================
!  (1 ) bpch2_mod.f    : Module w/ routines for binary punch file I/O
!  (2 ) dao_mod.f      : Module w/ arrays for DAO met fields
!  (3 ) file_mod.f     : Module w/ file unit numbers & error checks
!  (4 ) grid_mod.f     : Module w/ horizontal grid information   
!  (5 ) pbl_mix_mod.f  : Module w/ routines for PBL height & mixing
!  (6 ) pressure_mod.f : Module w/ routines to compute P(I,J,L)
!  (7 ) time_mod.f     : Module w/ routines for computing time & date 
!  (8 ) tracer_mod.f   : Module w/ GEOS-CHEM tracer array STT etc.  
!  (9 ) tracerid_mod.f : Module w/ pointers to tracers & emissions
!
!  ND48 tracer numbers:
!  ============================================================================
!  1 - nAdvect   : GEOS-CHEM advected species               [v/v      ]
!  501           : OH concentration                         [molec/cm3]
!  502           : NOy concentration                        [v/v      ]
!  503           : Relative Humidity                        [%        ]
!  504           : 3-D Cloud fractions                      [unitless ]
!  505           : Column optical depths                    [unitless ]
!  506           : Cloud top heights                        [hPa      ]
!  507           : Air density                              [molec/cm3]
!  508           : Total seasalt tracer concentration       [unitless ]
!  509           : PBL heights                              [m        ]
!  510           : PBL heights                              [levels   ]
!  511           : Grid box heights                         [m        ]
!  512           : PEDGE-$ (Pressure @ level edges)         [hPa      ]
!                  (lower edge, moist air pressure) 
!  513           : Sea level pressure                       [hPa      ]
!  514           : Zonal wind (a.k.a. U-wind)               [m/s      ]
!  515           : Meridional wind (a.k.a. V-wind)          [m/s      ]
!  516           : Temperature                              [K        ]
!  517           : Sulfate aerosol optical depth            [unitless ]
!  518           : Black carbon aerosol optical depth       [unitless ]
!  519           : Organic carbon aerosol optical depth     [unitless ]
!  520           : Accumulation mode seasalt optical depth  [unitless ]
!  521           : Coarse mode seasalt optical depth        [unitless ]
!  522           : Total dust optical depth                 [unitless ]
!  523-529       : Size resolved dust optical depth         [unitless ]
!  
!  NOTES:
!  (1 ) Now save out cld frac and grid box heights (bmy, 4/20/05)
!  (2 ) Remove TRCOFFSET because it's always zero.  Now call GET_HALFPOLAR
!        to get the value for GEOS or GCAP grids. (bmy, 6/28/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now references XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (5 ) Minor bug fixes in DIAG48 (cdh, bmy, 2/11/08)
!  (6 ) Bug fix: replace "PS-PTOP" with "PEDGE-$" (phs, bmy, 10/7/08)
!  (7 ) Modified to archive O3, NO, NOy as tracers 89, 90, 91  (tmf, 10/22/07)
!  (8 ) Bug fix in DIAG49 for diangostic output of SLP. (bmy, 10/13/09)
!  (9 ) Modify AOD output to wavelength specified in jv_spec_aod.dat 
!       (clh, 05/07/10)
!  (10) Now do not scale AOD output (now recalculated in RDAER and DUST_MOD)
!       Minor bug fix with the accumulation of ODAER across RH bins
!       (skim, 02/03/11)
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  06 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  06 Aug 2012 - R. Yantosca - Now make IU_ND48 a local module variable
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  19 Mar 2014 - M. Sulprizio- Updated to allow for more than 75 tracers
!  10 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE calls with State_Met%PEDGE. 
!                              Remove dependency on pressure_mod.
!  25 Mar 2015 - E. Lundgren - Change tracer units from kg to kg/kg and
!                              remove AD pointer
!******************************************************************************
!
      USE inquireMod,    ONLY : findFreeLUN
      USE PRECISION_MOD       ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "diag48_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE
 
      ! ... and these routines
      PUBLIC :: CLEANUP_DIAG48
      PUBLIC :: DIAG48
      PUBLIC :: INIT_DIAG48 
      PUBLIC :: ITS_TIME_FOR_DIAG48

      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      ! Scalars
      INTEGER              :: HALFPOLAR
      INTEGER, PARAMETER   :: CENTER180=1 
      REAL*4               :: LONRES
      REAL*4               :: LATRES
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: RESERVED = 'ND48 station timeseries'
      CHARACTER(LEN=80)    :: TITLE

      ! LUN for ND48 diagnostic file
      INTEGER              :: IU_ND48

      ! Arrays
      INTEGER, ALLOCATABLE :: ND48_I(:)
      INTEGER, ALLOCATABLE :: ND48_J(:)
      INTEGER, ALLOCATABLE :: ND48_L(:)
      INTEGER, ALLOCATABLE :: ND48_N(:)

      INTEGER              :: id_HNO3, id_N2O5, id_HNO4, id_NO
      INTEGER              :: id_PAN,  id_NPMN, id_PPN,  id_O3
      INTEGER              :: id_R4N2, id_SALA, id_SALC, id_NO2
      INTEGER              :: id_IPMN, id_OH
#endif

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
      SUBROUTINE DIAG48( am_I_Root, Input_Opt, 
     &                   State_Met, State_Chm, RC )
!
!******************************************************************************
!  Subroutine DIAG48 saves station time series diagnostics to disk.
!  (bmy, bey, amf, 6/1/99, 10/13/09)
!
!  NOTES:
!  (1 ) Remove reference to "CMN".  Also now get PBL heights in meters and
!        model layers from GET_PBL_TOP_m and GET_PBL_TOP_L of "pbl_mix_mod.f".
!        (bmy, 2/16/05)
!  (2 ) Now reference CLDF and BXHEIGHT from "dao_mod.f".  Now save 3-D cloud 
!        fraction as tracer #79 and box height as tracer #93.  Now remove
!        reference to PBL from "dao_mod.f" (bmy, 4/20/05)
!  (3 ) Remove references to TRCOFFSET because it's always zero.  Now call 
!        GET_HALFPOLAR from "bpch2_mod.f" to get the HALFPOLAR flag value for 
!        GEOS or GCAP grids.  (bmy, 6/28/05)
!  (4 ) Now do not save SLP data if it is not allocated (bmy, 8/2/05)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Now references XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (7 ) Bug fix: unit for tracer #77 should be "layers".  Also RH should be 
!        tracer #17 under "TIME-SER" category. (cdh, bmy, 2/11/08)
!  (8 ) Bug fix: replace "PS-PTOP" with "PEDGE-$". (phs, bmy, 10/7/08)
!  (9 ) Now save SLP under tracer #18 in "DAO-FLDS" (tai, bmy, 10/13/09)
!  12 Nov 2010 - R. Yantosca - Save out PEDGE-$ (pressure @ level edges) rather
!                              than PSURFACE - PTOP.
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  14 Mar 2013 - M. Payer    - Replace NOx and Ox with NO, NO2, and O3 as part
!                              of removal of NOx-Ox partitioning
!  06 Nov 2014 - R. Yantosca - Now use State_Met%AIRDEN(I,J,L)
!  06 Nov 2014 - R. Yantosca - Now use State_Met%CLDF(I,J,L)
!  06 Nov 2014 - R. Yantosca - Now use State_Met%OPTD(I,J,L)
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  24 Mar 2015 - E. Lundgren - Remove dependency on tracer_mod since
!                              XNUMOLAIR now defined in CMN_GTCM_MOD
!  22 Jun 2016 - M. Yannetti - Replace TCVV with spec db and physical const
!  20 Sep 2016 - R. Yantosca - Define species flags here
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState,           Ind_
      USE State_Met_Mod,      ONLY : MetState
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      USE BPCH2_MOD,          ONLY : BPCH2
      USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR
      USE CMN_FJX_MOD,        ONLY : ODAER, ODMDUST, NRH, NDUST
      USE CMN_FJX_MOD,        ONLY : IWVSELECT, ACOEF_WV, BCOEF_WV
!      USE CMN_O3_MOD               ! SAVEOH, XNUMOLAIR
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GC_GRID_MOD,        ONLY : GET_XOFFSET,        GET_YOFFSET
      USE PBL_MIX_MOD,        ONLY : GET_PBL_TOP_L,      GET_PBL_TOP_m
      USE PhysConstants            ! SCALE_HEIGHT, XNUMOLAIR
      USE TIME_MOD,           ONLY : GET_LOCALTIME,      GET_NYMD
      USE TIME_MOD,           ONLY : GET_NHMS,           GET_TAU
      USE TIME_MOD,           ONLY : EXPAND_DATE,        ITS_A_NEW_DAY
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !LOCAL VARIABLES:
!
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      LOGICAL, SAVE      :: IS_CLDTOPS,  IS_OPTD, IS_SEASALT, IS_SLP
      LOGICAL, SAVE      :: IS_FULLCHEM, IS_NOy
      LOGICAL, SAVE      :: FIRST = .TRUE.
      LOGICAL            :: LINTERP
      INTEGER            :: nAdvect
      INTEGER            :: GMTRC, I, I0, J, J0, L, N, K, R, TMP, W,ISPC
      REAL(fp)           :: LT, Q(LLPAR), SCALEAODnm
      REAL(f8)           :: TAU
      CHARACTER(LEN=40)  :: CATEGORY
      CHARACTER(LEN=40)  :: UNIT
      CHARACTER(LEN=255) :: FILENAME

      ! Aerosol types (rvm, aad, bmy, 7/20/04)
      INTEGER            :: IND(6) = (/ 22, 29, 36, 43, 50, 15 /)

      ! For fields from Input_Opt
      LOGICAL            :: DO_DIAG_WRITE

      ! Pointers
      REAL(fp), POINTER :: Spc(:,:,:,:)
#endif
      !=================================================================
      ! DIAG48 begins here!
      !=================================================================

      ! Assume success
      RC            =  GC_SUCCESS

#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )

      ! Number of advected species
      nAdvect       = State_Chm%nAdvect

      ! Copy values from Input_Opt
      DO_DIAG_WRITE = Input_Opt%DO_DIAG_WRITE

      ! Initialize pointers
      SPC          => State_Chm%Species

      ! Get grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Define species ID's
         ! These have to be done on the first call to DIAG48 so that the
         ! species database will have been already initialized (bmy, 9/20/16)
         id_HNO3 = Ind_('HNO3')
         id_N2O5 = Ind_('N2O5')
         id_HNO4 = Ind_('HNO4')
         id_NO   = Ind_('NO'  )
         id_PAN  = Ind_('PAN' )
         id_NPMN = Ind_('NPMN')
         id_IPMN = Ind_('IPMN')
         id_PPN  = Ind_('PPN' )
         id_O3   = Ind_('O3'  )
         id_R4N2 = Ind_('R4N2')
         id_SALA = Ind_('SALA')
         id_SALC = Ind_('SALC')
         id_NO2  = Ind_('NO2' )
         id_OH   = Ind_('OH'  )

         ! Set logical flags on first timestep
         IS_CLDTOPS  = ASSOCIATED( State_Met%CLDTOPS )
         IS_SLP      = ASSOCIATED( State_Met%SLP     )
         IS_OPTD     = ASSOCIATED( State_Met%OPTD    )
         IS_FULLCHEM = Input_Opt%ITS_A_FULLCHEM_SIM
         IS_SEASALT  = ( id_SALA > 0 .and. id_SALC > 0 )
         IS_NOy      = ( IS_FULLCHEM .and. id_NO   > 0 .and.
     &                   id_NO2  > 0 .and. id_PAN  > 0 .and.
     &                   id_HNO3 > 0 .and. id_NPMN > 0 .and.
     &                   id_PPN  > 0 .and. id_R4N2 > 0 .and.
     &                   id_N2O5 > 0 .and. id_HNO4 > 0 .and.
     &                   id_IPMN > 0 ) 

         ! Replace YYYYMMDD tokens in filename
         FILENAME = TRIM( Input_Opt%ND48_FILE )
         CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )

         IF ( DO_DIAG_WRITE ) THEN
           
            ! Open file for writing
            CALL OPEN_ND48_FILE( FILENAME )

         ENDIF

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Store station timeseries data in the Q array
      !=================================================================

      ! Get TAU value
      TAU = GET_TAU()

      !Determine if optical properties need interpolating
      !The LUT wavelengths in IWVSELECT will match if no interpolation
      !is needed. (output is only for the first requested wavelength)
      !(DAR 10/2013)
      IF(IWVSELECT(1,1).EQ.IWVSELECT(2,1)) THEN
         LINTERP=.FALSE.
      ELSE
         LINTERP=.TRUE.
      ENDIF

      ! Loop over each station
      DO W = 1, Input_Opt%ND48_N_STA

         ! Get lon, lat, alt, tracer indices
         I = ND48_I(W)             
         J = ND48_J(W)  
         K = ND48_L(W) 
         N = ND48_N(W)

         ! Get local time
         LT = GET_LOCALTIME( I, J, K )

         ! Initialize
         Q(:) = 0

         IF ( N <= nAdvect ) THEN

            !------------------------------------
            ! GEOS-CHEM advected species [v/v]
            !------------------------------------
            CATEGORY = 'IJ-AVG-$'
            UNIT     = ''              ! Let GAMAP pick the unit
            GMTRC    = N

            DO L = 1, K
               Q(L) = SPC(I,J,L,N) * ( AIRMW / 
     &                State_Chm%SpcData(N)%Info%emMW_g )
            ENDDO

! NOTE: Restore this when we figure out what the right unit conversion is
! This can be replaced by the netCDF diagnostic (bmy, 11/20/17)
!         ELSE IF ( N == 501 .and. IS_FULLCHEM ) THEN
!
!            !-------------------------------------
!            ! OH CONCENTRATION [molec/cm3]
!            !-------------------------------------
!            CATEGORY = 'TIME-SER'
!            UNIT     = 'molec/cm3'
!            GMTRC    = 1
!
!            DO L = 1, K
!               Q(L) = SAVEOH(I,J,L)
!            ENDDO

         ELSE IF ( N == 502 .and. IS_NOy ) THEN

            !-------------------------------------
            ! NOy CONCENTRATION [v/v]
            !-------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMTRC    = 2

            DO L = 1, K

               ! Initialize
               Q(L) = 0e+0_fp

               ! NO
               Q(L) = Q(L) + ( ( AIRMW 
     &                / State_Chm%SpcData(id_NO)%Info%emMW_g )
     &                         * SPC(I,J,L,id_NO) )

               ! NO2
               Q(L) = Q(L) + ( ( AIRMW 
     &                 / State_Chm%SpcData(id_NO2)%Info%emMW_g )         
     &                         * SPC(I,J,L,id_NO2) )

               ! PAN
               Q(L) = Q(L) + ( ( AIRMW 
     &                 / State_Chm%SpcData(id_PAN)%Info%emMW_g )         
     &                         * SPC(I,J,L,id_PAN) )

               ! HNO3
               Q(L) = Q(L) + ( ( AIRMW 
     &                 / State_Chm%SpcData(id_HNO3)%Info%emMW_g )       
     &                         * SPC(I,J,L,id_HNO3) )
            
               ! NPMN
               Q(L) = Q(L) + ( ( AIRMW 
     &                 / State_Chm%SpcData(id_NPMN)%Info%emMW_g )        
     &                         * SPC(I,J,L,id_NPMN) )

               ! IPMN
               Q(L) = Q(L) + ( ( AIRMW 
     &                 / State_Chm%SpcData(id_IPMN)%Info%emMW_g )        
     &                         * SPC(I,J,L,id_IPMN) )

               ! PPN
               Q(L) = Q(L) + ( ( AIRMW 
     &                 / State_Chm%SpcData(id_PPN)%Info%emMW_g )         
     &                         * SPC(I,J,L,id_PPN) )
 
               ! R4N2
               Q(L) = Q(L) + ( ( AIRMW 
     &                 / State_Chm%SpcData(id_R4N2)%Info%emMW_g )       
     &                         * SPC(I,J,L,id_R4N2) )
            
               ! N2O5
               Q(L) = Q(L) + ( 2d0 * ( AIRMW 
     &                 / State_Chm%SpcData(id_N2O5)%Info%emMW_g ) 
     &                         * SPC(I,J,L,id_N2O5) )
                        
               ! HNO4
               Q(L) = Q(L) + ( ( AIRMW 
     &                 / State_Chm%SpcData(id_HNO4)%Info%emMW_g )       
     &                         * SPC(I,J,L,id_HNO4) )
            ENDDO

         ELSE IF ( N == 503 ) THEN 

            !-------------------------------------
            ! RELATIVE HUMIDITY [%]
            !-------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = '%'
            GMTRC    = 3

            DO L = 1, K
               Q(L) = State_Met%RH(I,J,L)
            ENDDO

         ELSE IF ( N == 504 ) THEN

            !-------------------------------------
            ! 3-D CLOUD FRACTIONS [unitless]
            !-------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'unitless'
            GMTRC    = 4

            DO L = 1, K
               Q(L) = State_Met%CLDF(I,J,L)
            ENDDO

         ELSE IF ( N == 505 .and. IS_OPTD ) THEN
            
            !-------------------------------------
            ! COLUMN OPTICAL DEPTHS [unitless]
            !-------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'unitless'
            GMTRC    = 5

            Q(1) = SUM( State_Met%OPTD(I,J,:) )

         ELSE IF ( N == 506 .and. IS_CLDTOPS ) THEN

            !-------------------------------------
            ! CLOUD TOP HEIGHTS [hPa]
            !-------------------------------------
            CATEGORY = 'TIME_SER'
            UNIT     = 'hPa'
            GMTRC    = 6

            IF ( K == 1 ) THEN
               Q(1) = State_Met%PEDGE(I,J,State_Met%CLDTOPS(I,J))
            ENDIF

         ELSE IF ( N == 507 ) THEN 

            !-------------------------------------
            ! AIR DENSITY [molec/cm3]
            !-------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'molec/cm3'
            GMTRC    = 7

            DO L = 1, K
               Q(L) = State_Met%AIRDEN(I,J,L) * XNUMOLAIR * 1e-6_fp
            ENDDO

         ELSE IF ( N == 508 .and. IS_SEASALT ) THEN

            !-------------------------------------
            ! TOTAL SEASALT TRACER [v/v]
            !-------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMTRC    = 8

            DO L = 1, K
               Q(L) = ( SPC(I,J,L,id_SALA) + SPC(I,J,L,id_SALC) ) *
     &                  (AIRMW / State_Chm%SpcData(id_SALA)%Info%emMW_g)
            ENDDO

         ELSE IF ( N == 509 ) THEN

            !-------------------------------------
            ! PBL HEIGHTS [m] 
            !-------------------------------------
            CATEGORY = 'PBLDEPTH'
            UNIT     = 'm'  
            GMTRC    = 1

            IF ( K == 1 ) THEN
               Q(1) = GET_PBL_TOP_m( I, J )
            ENDIF
               
         ELSE IF ( N == 510 ) THEN

            !-------------------------------------
            ! PBL HEIGHTS [levels] 
            !-------------------------------------   
            CATEGORY = 'PBLDEPTH'
            UNIT     = 'levels'
            GMTRC    = 2

            IF ( K == 1 ) THEN
               Q(1) = GET_PBL_TOP_L( I, J )
            ENDIF

         ELSE IF ( N == 511 ) THEN 

            !-------------------------------------
            ! GRID BOX HEIGHT [m]
            !-------------------------------------
            CATEGORY = 'BXHGHT-$'
            UNIT     = 'm'
            GMTRC    = 1

            DO L = 1, K
               Q(L) = State_Met%BXHEIGHT(I,J,L)
            ENDDO

         ELSE IF ( N == 512 ) THEN

            !-------------------------------------
            ! PEDGE-$ [hPa]
            !-------------------------------------
            CATEGORY = 'PEDGE-$'
            UNIT     = 'hPa'
            GMTRC    = 1            
            
            DO L = 1, K
               Q(L) = State_Met%PEDGE( I, J, L )
            ENDDO

         ELSE IF ( N == 513 .and. IS_SLP ) THEN

            !-------------------------------------
            ! SEA LEVEL PRESSURE [hPa]
            !-------------------------------------
            CATEGORY = 'DAO-FLDS'
            UNIT     = 'hPa'
            GMTRC    = 18

            IF ( K == 1 ) THEN
               Q(1) = State_Met%SLP(I,J)
            ENDIF

         ELSE IF ( N == 514 ) THEN

            !-------------------------------------
            ! ZONAL (U) WIND [m/s]
            !-------------------------------------
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'm/s'
            GMTRC    = 1

            DO L = 1, K
               Q(L) = State_Met%U(I,J,L)
            ENDDO
            
         ELSE IF ( N == 515 ) THEN 

            !-------------------------------------
            ! MERIDIONAL (V) WIND [m/s]
            !-------------------------------------
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'm/s'
            GMTRC    = 2

            DO L = 1, K
               Q(L) = State_Met%V(I,J,L)
            ENDDO

         ELSE IF ( N == 516 ) THEN

            !-------------------------------------
            ! TEMPERATURE [K]
            !-------------------------------------
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'K'
            GMTRC    = 3            

            DO L = 1, K
               Q(L) = State_Met%T(I,J,L)
            ENDDO

         ELSE IF ( N == 517 ) THEN

            !-------------------------------------
            ! SULFATE AOD [unitless]
            ! for wavelengths set in Radiation Menu
            !-------------------------------------
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMTRC    = 6

            ISPC     = 1 !sulfate
            IF ( .not. LINTERP ) THEN
               DO L = 1, K
                  ! Accumulate
                  Q(L) = Q(L) + ODAER(I,J,L,IWVSELECT(1,1),ISPC)
               ENDDO
            ELSE
               DO L = 1, K
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(L) = Q(L) + 
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDDO
            ENDIF

         ELSE IF ( N == 518 ) THEN

            !-------------------------------------
            ! BLACK CARBON AOD [unitless]
            ! for wavelengths set in Radiation Menu
            !-------------------------------------
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMTRC    = 9            

            ISPC     = 2 !BC
            IF ( .not. LINTERP ) THEN
               DO L = 1, K
                  ! Accumulate
                  Q(L) = Q(L) + ODAER(I,J,L,IWVSELECT(1,1),ISPC)
               ENDDO
            ELSE
               DO L = 1, K
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(L) = Q(L) +
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDDO
            ENDIF

         ELSE IF ( N == 519 ) THEN

            !-------------------------------------
            ! ORGANIC CARBON AOD [unitless]
            ! for wavelengths set in Radiation Menu
            !-------------------------------------       
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMTRC    = 12

            ISPC     = 3 !OC
            IF ( .not. LINTERP ) THEN
               DO L = 1, K
                  ! Accumulate
                  Q(L) = Q(L) + ODAER(I,J,L,IWVSELECT(1,1),ISPC)
               ENDDO
            ELSE
               DO L = 1, K
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(L) = Q(L) +
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDDO
            ENDIF

         ELSE IF ( N == 520 ) THEN
            
            !-------------------------------------
            ! ACCUM SEASALT AOD [unitless]
            ! for wavelengths set in Radiation Menu
            !-------------------------------------
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMTRC    = 15

            ISPC     = 4 !SSa
            IF ( .not. LINTERP ) THEN
               DO L = 1, K
                  ! Accumulate
                  Q(L) = Q(L) + ODAER(I,J,L,IWVSELECT(1,1),ISPC)
               ENDDO
            ELSE
               DO L = 1, K
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(L) = Q(L) +
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDDO
            ENDIF

         ELSE IF ( N == 521 ) THEN

            !-------------------------------------
            ! COARSE SEASALT AOD [unitless]
            ! for wavelengths set in Radiation Menu
            !-------------------------------------    
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMTRC    = 18

            ISPC     = 5 !SSc
            IF ( .not. LINTERP ) THEN
               DO L = 1, K
                  ! Accumulate
                  Q(L) = Q(L) + ODAER(I,J,L,IWVSELECT(1,1),ISPC)
               ENDDO
            ELSE
               DO L = 1, K
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(L) = Q(L) +
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDDO
            ENDIF

         ELSE IF ( N == 522 ) THEN

            !-------------------------------------
            ! TOTAL DUST OPT DEPTH [unitless]
            ! for wavelengths set in Radiation Menu
            !-------------------------------------
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMTRC     = 4

            DO R = 1, NDUST

               IF ( .not. LINTERP ) THEN
                  DO L = 1, K
                     ! Accumulate
                     Q(L) = Q(L) + ODMDUST(I,J,L,IWVSELECT(1,1),R)
                  ENDDO
               ELSE
                  DO L = 1, K
                     ! Interpolated using angstrom exponent between
                     ! Closest available wavelengths
                     ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                     !catch any zero values before interpolation
                     IF ((ODMDUST(I,J,L,IWVSELECT(2,1),R).GT.0).AND.
     &                   (ODMDUST(I,J,L,IWVSELECT(1,1),R).GT.0)) THEN
                       Q(L) = Q(L) +
     &                 (ODAER(I,J,L,IWVSELECT(2,1),R)*ACOEF_WV(1)**
     &                 (BCOEF_WV(1)*LOG(ODMDUST(I,J,L,IWVSELECT(1,1),R)/
     &                                ODMDUST(I,J,L,IWVSELECT(2,1),R))))
                     ENDIF
                  ENDDO
               ENDIF

            ENDDO

         ELSE IF ( ( N >= 523 ) .and. ( N <= 529) ) THEN

            !-------------------------------------
            ! DUST BINS 1-7 OPT DEPTH [unitless]
            ! for wavelengths set in Radiation Menu
            !-------------------------------------
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMTRC     = 21+(N-523)

            R = N - 522

            IF ( .not. LINTERP ) THEN
               DO L = 1, K
                  ! Accumulate
                  Q(L) = Q(L) + ODMDUST(I,J,L,IWVSELECT(1,1),R)
               ENDDO
            ELSE
               DO L = 1, K
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODMDUST(I,J,L,IWVSELECT(2,1),R).GT.0).AND.
     &                (ODMDUST(I,J,L,IWVSELECT(1,1),R).GT.0)) THEN
                    Q(L) = Q(L) +
     &              (ODAER(I,J,L,IWVSELECT(2,1),R)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODMDUST(I,J,L,IWVSELECT(1,1),R)/
     &                               ODMDUST(I,J,L,IWVSELECT(2,1),R))))
                  ENDIF
               ENDDO
            ENDIF

         ELSE

            ! Skip other tracers
            CYCLE

         ENDIF
         
         IF ( DO_DIAG_WRITE ) THEN

         !==============================================================
         ! Write each station to a bpch file
         !==============================================================
            CALL BPCH2( IU_ND48,   MODELNAME, LONRES,   LATRES,       
     &           HALFPOLAR, CENTER180, CATEGORY, GMTRC,        
     &           UNIT,      TAU,       TAU,      RESERVED,  
     &           1,         1,         K,        I+I0,         
     &           J+J0,      1,         REAL( Q(1:K),4 ) )
            
         ENDIF

      ENDDO

      !=================================================================
      ! Close the file at the proper time
      !=================================================================
!      IF ( ITS_TIME_TO_CLOSE_FILE() ) THEN
!
!         ! Expand date tokens in the file name
!         FILENAME = TRIM( ND49_OUTPUT_FILE )
!         CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )
!
!         ! Echo info
!         WRITE( 6, 120 ) TRIM( FILENAME )
! 120     FORMAT( '     - DIAG49: Closing file : ', a )
!
!         ! Close file
!         CLOSE( IU_ND49 ) 
!      ENDIF

      IF ( DO_DIAG_WRITE ) THEN

         ! Flush the file once per day
         IF ( ITS_A_NEW_DAY() ) CALL FLUSH( IU_ND48 ) 

      ENDIF

      ! Free pointers
      Spc => NULL()
#endif

      END SUBROUTINE DIAG48

!------------------------------------------------------------------------------

      SUBROUTINE OPEN_ND48_FILE( FILENAME )
!
!******************************************************************************
!  Subroutine OPEN_ND48_FILE opens a binary punch file (version 2.0)
!  for writing GEOS-CHEM ND48 station timeseries data. (bmy, 7/30/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : Name of the file to be opened
!
!  NOTES:
!******************************************************************************
!
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      ! References to F90 modules
      USE FILE_MOD, ONLY : IOERROR
#endif

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME

#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )

      ! Local variables
      INTEGER                      :: IOS
      CHARACTER(LEN=40)            :: FTI
      CHARACTER(LEN=80)            :: TITLE

      !=================================================================
      ! OPEN_ND48_FOR_WRITE begins here!
      !=================================================================

      ! Find a free file LUN
      IU_ND48 = findFreeLUN()

      ! Initialize
      FTI   = 'CTM bin 02 -- GEOS-CHEM station ts' 
      TITLE = 'GEOS-CHEM ND48 station timeseries diagnostic'

!%%% NOTE: The GEOS-5 GCM won't call this routine, so it should be OK
!%%% to use IU_ND48 from the file_mod.F. (bmy, 8/3/12)
!%%%      ! Find a free file LUN
!%%%      IU_ND48 = findFreeLUN()

      ! Open file for output
      OPEN( IU_ND48,    FILE=TRIM( FILENAME ), STATUS='UNKNOWN',
     &      IOSTAT=IOS, FORM='UNFORMATTED',    ACCESS='SEQUENTIAL' )

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_ND48, 'open_nd48_file:1' )

      ! Write file type identifier
      WRITE ( IU_ND48, IOSTAT=IOS ) FTI
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_ND48, 'open_nd48_file:2' )

      ! Write top title
      WRITE ( IU_ND48, IOSTAT=IOS ) TITLE
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_ND48, 'open_nd48_file:3' )
#endif

      ! Return to calling program
      END SUBROUTINE OPEN_ND48_FILE

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_DIAG48( Input_Opt ) RESULT( ITS_TIME )
!
!******************************************************************************
!  Function ITS_TIME_FOR_DIAG48 returns TRUE if it's time for the next
!  timeseries data write. (bmy, 7/20/04)
!
!  NOTES:
!  (1 ) Add a check on the output frequency for validity compared to time 
!        steps used. (ccc, 5/21/09)
!******************************************************************************
!
      ! References to F90 modules
      USE Input_Opt_Mod, ONLY : OptInput
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      USE ERROR_MOD,     ONLY : GEOS_CHEM_STOP
      USE TIME_MOD,      ONLY : GET_ELAPSED_SEC, GET_TS_DIAG
#endif
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options
!
! !RETURN VALUE:
!
      LOGICAL :: ITS_TIME

#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )

      ! Local variables
      INTEGER :: XSEC, TS_DIAG
      LOGICAL, SAVE :: FIRST = .TRUE.

      !=================================================================
      ! ITS_TIME_FOR_DIAG48 begins here!
      !=================================================================
      
      IF ( Input_Opt%DO_ND48 ) THEN

         IF ( FIRST ) THEN
            TS_DIAG = GET_TS_DIAG()

            ! Check if ND48_FREQ is a multiple of TS_DIAG
            IF ( MOD( Input_Opt%ND48_FREQ, TS_DIAG ) /= 0 ) THEN
               WRITE( 6, 100 ) 'ND48', Input_Opt%ND48_FREQ, TS_DIAG
 100           FORMAT( 'The ',a,' output frequency must be a multiple '
     &              'of the largest time step:', i5, i5 )
               CALL GEOS_CHEM_STOP
            ENDIF
            FIRST = .FALSE.
         ENDIF
      
         ! Get elapsed seconds
         XSEC     = GET_ELAPSED_SEC()
         
         ! Is it time to save the next timeseries station?
         ITS_TIME = ( Input_Opt%DO_ND48 .and.
     &                MOD( XSEC, Input_Opt%ND48_FREQ ) == 0 )
      ELSE
         ITS_TIME = Input_Opt%DO_ND48
      ENDIF
#endif

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_DIAG48

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DIAG48( am_I_Root, Input_Opt, RC )
!
!******************************************************************************
!  Subroutine INIT_DIAG48 allocates and zeroes all module arrays (bmy, 7/20/04)
!
!  NOTES:
!  (1 ) Now call GET_HALFPOLAR from "bpch2_mod.f" to get the HALFPOLAR flag 
!        value for GEOS or GCAP grids. (bmy, 6/28/05)
!  22 Jun 2016 - R. Yantosca - Now define species ID's with IND()
!  20 Sep 2016 - R. Yantosca - Move species ID's initialization to DIAG48
!******************************************************************************
!
      ! References to F90 modules
      USE ErrCode_Mod
      USE Input_Opt_Mod,         ONLY : OptInput
      USE State_Chm_Mod,         ONLY : Ind_
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      USE BPCH2_MOD,             ONLY : GET_MODELNAME, GET_HALFPOLAR
      USE CMN_SIZE_MOD                ! Size parameters
      USE ErrCode_Mod
      USE ERROR_MOD,             ONLY : ALLOC_ERR, ERROR_STOP
#endif
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?

#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      ! Local variables
      INTEGER                     :: AS, N
      CHARACTER(LEN=255)          :: LOCATION
#endif

      !=================================================================
      ! INIT_DIAG48 begins here!
      !=================================================================
      
      ! Assume success
      RC = GC_SUCCESS

#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )

      ! Assume success
      RC = GC_SUCCESS

      ! Location string
      LOCATION         = 'INIT_DIAG48 ("diag48_mod.f")'

      ! Return if we are not saving ND48 diagnostics
      IF ( .not. Input_Opt%DO_ND48 ) RETURN

      !------------------------------
      ! Allocate module arrays
      !------------------------------
      ALLOCATE( ND48_I( Input_Opt%ND48_N_STA ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ND48_I' )

      ALLOCATE( ND48_J( Input_Opt%ND48_N_STA ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ND48_J' )

      ALLOCATE( ND48_L( Input_Opt%ND48_N_STA ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ND48_L' )

      ALLOCATE( ND48_N( Input_Opt%ND48_N_STA ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ND48_N' )

      !-------------------------------
      ! Copy values & error check
      !-------------------------------
      DO N = 1, Input_Opt%ND48_N_STA

         ! Error check longitude
         IF ( Input_Opt%ND48_IARR(N) < 1       .or.
     &        Input_Opt%ND48_IARR(N) > IIPAR ) THEN
            CALL ERROR_STOP( 'Bad longitude index!', LOCATION )
         ELSE
            ND48_I(N) = Input_Opt%ND48_IARR(N)
         ENDIF

         ! Error check latitude
         IF ( Input_Opt%ND48_JARR(N) < 1       .or.
     &        Input_Opt%ND48_JARR(N) > JJPAR ) THEN
            CALL ERROR_STOP( 'Bad latitude index!', LOCATION )
         ELSE
            ND48_J(N) = Input_Opt%ND48_JARR(N)
         ENDIF

         ! Error check longitude
         IF ( Input_Opt%ND48_LARR(N) < 1       .or.
     &        Input_Opt%ND48_LARR(N) > LLPAR ) THEN
            CALL ERROR_STOP( 'Bad altitude index!', LOCATION )
         ELSE
            ND48_L(N) = Input_Opt%ND48_LARR(N)
         ENDIF

         ! Tracer array
         ND48_N(N) = Input_Opt%ND48_NARR(N)
      ENDDO

      !-------------------------------
      ! For bpch file output
      !-------------------------------
      LONRES    = DISIZE
      LATRES    = DJSIZE
      MODELNAME = GET_MODELNAME()
      HALFPOLAR = GET_HALFPOLAR()
#endif

      ! Return to calling program
      END SUBROUTINE INIT_DIAG48

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DIAG48()
!
!******************************************************************************
!  Subroutine CLEANUP_DIAG48 deallocates all module arrays (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      !=================================================================
      ! CLEANUP_DIAG48 begins here!
      !=================================================================
      IF ( ALLOCATED( ND48_I ) ) DEALLOCATE( ND48_I )
      IF ( ALLOCATED( ND48_J ) ) DEALLOCATE( ND48_J )
      IF ( ALLOCATED( ND48_L ) ) DEALLOCATE( ND48_L )
      IF ( ALLOCATED( ND48_N ) ) DEALLOCATE( ND48_N )
#endif

      ! Return to calling program
      END SUBROUTINE CLEANUP_DIAG48

!------------------------------------------------------------------------------

      ! Return to calling program
      END MODULE DIAG48_MOD
