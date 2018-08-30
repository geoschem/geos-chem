!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag50_mod.F
!
! !DESCRIPTION: Module DIAG50\_MOD contains variables and routines to 
!  generate 24-hour average timeseries data.
!\\
!\\
! !INTERFACE: 
!
      MODULE DIAG50_MOD
!
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC  :: CLEANUP_DIAG50
      PUBLIC  :: DIAG50
      PUBLIC  :: INIT_DIAG50
!
! !PRIVATE MEMBER FUNCTIONS:
! 
      PRIVATE :: ACCUMULATE_DIAG50
      PRIVATE :: WRITE_DIAG50
      PRIVATE :: GET_I
!
! !REMARKS:
!  ND50 tracer numbers:
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
!  512           : PEDGE-$ (Pressure @ level edges          [hPa      ]
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
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Rewritten for clarity and to save extra quantities (bmy, 7/20/04)
!  (2 ) Added COUNT_CHEM to count the chemistry timesteps per day, since some
!        quantities are only archived after a fullchem call (bmy, 10/25/04)
!  (3 ) Bug fix: Now get I0 and J0 properly for nested grids (bmy, 11/9/04)
!  (4 ) Now only archive AODs once per chemistry timestep (bmy, 1/14/05)
!  (5 ) Now references "pbl_mix_mod.f" (bmy, 2/16/05)
!  (6 ) Now save cloud fractions & grid box heights (bmy, 4/20/05)
!  (7 ) Remove TRCOFFSET since it is always zero.  Also now get HALFPOLAR for
!        both GCAP and GEOS grids. (bmy, 6/24/05)
!  (8 ) Bug fix: do not save SLP unless it is allocated (bmy, 8/2/05)
!  (9 ) Now references XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (10) Modified INIT_DIAG49 to save out transects (cdh, bmy, 11/30/06)
!  (11) Now use 3D timestep counter for full chem in the trop (phs, 1/24/07)
!  (12) Renumber RH diagnostic in WRITE_DIAG50 (bmy, 2/11/08)
!  (13) Bug fix: replace "PS-PTOP" with "PEDGE-$" (bmy, 10/7/08)
!  (14) Modified to archive O3, NO, NOy as tracers 89, 90, 91  (tmf, 9/26/07)
!  (15) Updates & bug fixes in WRITE_DIAG50 (ccc, tai, bmy, 10/13/09)
!  (16) Updates to AOD output.  Also have the option to write to HDF 
!        (amv, bmy, 12/21/09)
!  (17) Modify AOD output to wavelength specified in jv_spec_aod.dat 
!       (clh, 05/07/10)
!  12 Nov 2010 - R. Yantosca - Now save out PEDGE-$ (pressure at level edges)
!                              rather than Psurface - PTOP
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  03 Feb 2011 - S. Kim      - Now do not scale the AOD output
!                              (recalculated in RDAER AND DUST_MOD)
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  06 Aug 2012 - R. Yantosca - Now make IU_ND50 a local module variable
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  19 Mar 2014 - M. Sulprizio- Updated to allow for more than 75 tracers
!  10 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      !=================================================================
      ! MODULE VARIABLES
      !
      ! COUNT            : Counter of timesteps
      ! I0               : Offset between global & nested grid
      ! J0               : Offset between global & nested grid
      ! IOFF             : Longitude offset
      ! JOFF             : Latitude offset
      ! LOFF             : Altitude offset
      ! ND50_NI          : Number of longitudes in DIAG51 region 
      ! ND50_NJ          : Number of latitudes  in DIAG51 region
      ! ND50_NL          : Number of levels     in DIAG51 region
      ! Q                : Accumulator array for various quantities
      ! TAU0             : Starting TAU used to index the bpch file
      ! TAU1             : Ending TAU used to index the bpch file
      ! HALFPOLAR        : Used for bpch file output
      ! CENTER180        : Used for bpch file output
      ! LONRES           : Used for bpch file output
      ! LATRES           : Used for bpch file output
      ! MODELNAME        : Used for bpch file output
      ! RESERVED         : Used for bpch file output
      !=================================================================

      ! Scalars
      INTEGER              :: IU_ND50
      INTEGER              :: COUNT
      INTEGER              :: IOFF,           JOFF   
      INTEGER              :: LOFF,           I0
      INTEGER              :: J0,             ND50_NI
      INTEGER              :: ND50_NJ,        ND50_NL
      INTEGER              :: HALFPOLAR
      INTEGER, PARAMETER   :: CENTER180=1
      REAL*4               :: LONRES,         LATRES
      REAL(f8)             :: TAU0,           TAU1
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE

      ! Arrays
      REAL(fp),  ALLOCATABLE :: Q(:,:,:,:)
      INTEGER, ALLOCATABLE :: COUNT_CHEM3D(:,:,:)
#endif

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DIAG50
!
! !DESCRIPTION: Subroutine DIAG50 generates 24hr average time series.  
!  Output is to binary punch file format or HDF5 file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DIAG50( am_I_Root, Input_Opt, State_Met, State_Chm, RC)
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE State_Met_Mod,      ONLY : MetState
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      USE Time_Mod,           ONLY : ITS_A_NEW_DAY
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
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  06 Nov 2014 - R. Yantosca - Now use State_Met%AIRDEN(I,J,L)
!EOP
!------------------------------------------------------------------------------
!BOC


      !=================================================================
      ! DIAG50 begins here!
      !=================================================================

      ! Assume success
      RC = GC_SUCCESS

#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )

      ! Accumulate data over a 24-hr period in the Q array
      CALL ACCUMULATE_DIAG50
     &    ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

      ! Write data to disk at the end of the day
      IF ( ITS_A_NEW_DAY() .and. 
     &     Input_Opt%DO_DIAG_WRITE    ) THEN
         CALL WRITE_DIAG50( am_I_Root, Input_Opt, State_Chm, RC )
      ENDIF

#endif

      END SUBROUTINE DIAG50
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: accumulate_diag50
!
! !DESCRIPTION: Subroutine ACCUMULATE\_DIAG50 accumulates tracers into the Q 
!  array. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ACCUMULATE_DIAG50( am_I_Root, Input_Opt,
     &                              State_Met, State_Chm, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_
      USE State_Met_Mod,      ONLY : MetState
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR
      USE CMN_FJX_MOD,        ONLY : ODAER, ODMDUST, NRH, NDUST
      USE CMN_FJX_MOD,        ONLY : IWVSELECT, ACOEF_WV, BCOEF_WV
!      USE CMN_O3_MOD               ! SAVEOH
      USE PhysConstants            ! SCALE_HEIGHT, XNUMOLAIR
      USE PBL_MIX_MOD,        ONLY : GET_PBL_TOP_L, GET_PBL_TOP_m
      USE TIME_MOD,           ONLY : GET_ELAPSED_SEC, GET_TS_CHEM
      USE TIME_MOD,           ONLY : TIMESTAMP_STRING, GET_TS_DYN
      USE Aerosol_Mod,        ONLY : PM25
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
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Rewrote to remove hardwiring and for better efficiency.  Added extra
!        diagnostics and updated numbering scheme.  Now scale aerosol & dust
!        optical depths to 400 nm. (rvm, aad, bmy, 7/20/04) 
!  (2 ) Now reference GET_ELAPSED_MIN and GET_TS_CHEM from "time_mod.f".  
!        Also now use extra counter COUNT_CHEM to count the number of
!        chemistry timesteps since NO, NO2, OH, O3 only when a full-chemistry
!        timestep happens. (bmy, 10/25/04)
!  (3 ) Only archive AODs when it is a chem timestep (bmy, 1/14/05)
!  (4 ) Remove reference to "CMN".  Also now get PBL heights in meters and 
!        model layers from GET_PBL_TOP_m and GET_PBL_TOP_L of "pbl_mix_mod.f".
!        (bmy, 2/16/05)
!  (5 ) Now reference CLDF and BXHEIGHT from "dao_mod.f".  Now save 3-D
!        cloud fraction as tracer #79 and box height as tracer #93.  Now 
!        remove references to CLMOSW, CLROSW, and PBL from "dao_mod.f". 
!        (bmy, 4/20/05)
!  (6 ) Remove references to TRCOFFSET because it is always zero (bmy, 6/24/05)
!  (7 ) Now do not save SLP data if it is not allocated (bmy, 8/2/05)
!  (8 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (9 ) Now references XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (10) Now account for time spent in the trop for non-tracers (phs, 1/24/07)
!  (11) IS_CHEM check is not appropriate anymore. Keep COUNT_CHEM3D for 
!       species known in troposphere only (ccc, 8/12/09)
!  (12) Output AOD at 3rd jv_spec.dat row wavelength.  Include all seven dust 
!        bin individual AOD (amv, bmy, 12/21/09)
!  12 Nov 2010 - R. Yantosca - Now save out PEDGE-$ (pressure at level edges)
!                              rather than Psurface - PTOP
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  14 Mar 2013 - M. Payer    - Replace NOx and Ox with NO, NO2, and O3 as part
!                              of removal of NOx-Ox partitioning
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  06 Nov 2014 - R. Yantosca - Now use State_Met%CLDF(I,J,L)
!  06 Nov 2014 - R. Yantosca - Now use State_Met%OPTD(I,J,L)
!  26 Feb 2014 - E. Lundgren - Replace GET_PEDGE with State_Met%PEDGE.
!  24 Mar 2015 - E. Lundgren - Remove dependency on tracer_mod
!  25 Mar 2015 - E. Lundgren - Change tracer units from kg to kg/kg
!  17 May 2016 - M. Sulprizio- Add IS_CHEM logical for incrementing COUNT_CHEM3D
!  16 Jun 2016 - K. Yu       - Now define species IDs with the Ind_ function
!  17 Jun 2016 - R. Yantosca - Only define species IDs on the first call
!  22 Jun 2016 - M. Yannetti - Replace TCVV with spec db and physical const
!  28 Jun 2016 - M. Sulprizio- Replace NPVERT with LLPAR to remove dependence
!                              on comode_loop_mod.F
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  11 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!  06 Feb 2018 - E. Lundgren - Change GET_ELAPSED_MIN to GET_ELAPSED_SEC to
!                              match new timestep unit of seconds
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )

      ! SAVEd scalars
      LOGICAL, SAVE       :: FIRST = .TRUE.
      LOGICAL, SAVE       :: IS_FULLCHEM, IS_SEASALT
      LOGICAL, SAVE       :: IS_CLDTOPS,  IS_NOy,  IS_OPTD, IS_SLP
      INTEGER, SAVE       :: id_HNO3,     id_HNO4, id_N2O5, id_NO
      INTEGER, SAVE       :: id_PAN,      id_NPMN,  id_PPN,  id_O3
      INTEGER, SAVE       :: id_R4N2,     id_SALA, id_SALC, id_NO2
      INTEGER, SAVE       :: id_IPMN,     id_OH

      ! Scalars
      LOGICAL             :: IS_CHEM
      LOGICAL             :: LINTERP
      INTEGER             :: H, I, J, K, L, M, N, nAdvect
      INTEGER             :: PBLINT,  R, X, Y, W
      INTEGER             :: ISPC
      REAL(fp)            :: C1, C2, PBLDEC, TEMPBL, TMP, SCALEAODnm
      CHARACTER(LEN=16)   :: STAMP

      ! Pointers
      REAL(fp), POINTER   :: Spc(:,:,:,:)
#endif
      !=================================================================
      ! ACCUMULATE_DIAG50 begins here!
      !=================================================================

      ! Assume success
      RC        =  GC_SUCCESS

#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )

      ! Number of advected species
      nAdvect = State_Chm%nAdvect

      ! Initialize pointers
      Spc      => State_Chm%Species

      ! First-time setup
      IF ( FIRST ) THEN

         ! Define species ID flags
         id_HNO3 = Ind_('HNO3')
         id_HNO4 = Ind_('HNO4')
         id_N2O5 = Ind_('N2O5')
         id_NO   = Ind_('NO'  )
         id_PAN  = Ind_('PAN' )
         id_NPMN = Ind_('NPMN' )
         id_IPMN = Ind_('IPMN' )
         id_PPN  = Ind_('PPN' )
         id_O3   = Ind_('O3'  )
         id_R4N2 = Ind_('R4N2')
         id_SALA = Ind_('SALA')
         id_SALC = Ind_('SALC')
         id_NO2  = Ind_('NO2' )
         id_OH   = Ind_('OH'  )

         ! Set logical flags
         IS_OPTD     = ASSOCIATED( State_Met%OPTD    )
         IS_CLDTOPS  = ASSOCIATED( State_Met%CLDTOPS )
         IS_SLP      = ASSOCIATED( State_Met%SLP     )
         IS_FULLCHEM = Input_Opt%ITS_A_FULLCHEM_SIM
         IS_SEASALT  = ( id_SALA > 0 .and. id_SALC > 0 )
         IS_NOy      = ( IS_FULLCHEM .and. id_NO   > 0 .and.
     &                   id_NO2  > 0 .and. id_PAN  > 0 .and.
     &                   id_HNO3 > 0 .and. id_NPMN > 0 .and.
     &                   id_PPN  > 0 .and. id_R4N2 > 0 .and.
     &                   id_N2O5 > 0 .and. id_HNO4 > 0 .and.
     &                   id_IPMN > 0 ) 

         ! Reset first-time flag
         FIRST       = .FALSE.
      ENDIF

      ! Is it a chemistry timestep?
      IS_CHEM = ( MOD( GET_ELAPSED_SEC()-GET_TS_DYN(), 
     &                 GET_TS_CHEM() ) == 0 )

      ! Echo time information to the screen
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 100 ) STAMP
 100  FORMAT( '     - DIAG50: Accumulation at ', a )
      
      !=================================================================
      ! Archive tracers into accumulating array Q 
      !=================================================================

      ! Increment counter
      COUNT = COUNT + 1

      !Determine if optical properties need interpolating
      !The LUT wavelengths in IWVSELECT will match if no interpolation
      !is needed. (output is only for the first requested wavelength)
      !(DAR 10/2013)
      IF(IWVSELECT(1,1).EQ.IWVSELECT(2,1)) THEN
         LINTERP=.FALSE.
      ELSE
         LINTERP=.TRUE.
      ENDIF

      ! Also increment 3-D counter for boxes in the chemistry grid
      IF ( IS_FULLCHEM .and. IS_CHEM ) THEN
         
         ! Loop over levels
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( X, Y, K, I, J, L )
!$OMP+SCHEDULE( DYNAMIC )
         DO K = 1, MIN(ND50_NL, LLPAR)
            L = LOFF + K

         ! Loop over latitudes 
         DO Y = 1, ND50_NJ
            J = JOFF + Y

         ! Loop over longitudes
         DO X = 1, ND50_NI
            I = GET_I( X )

            ! Only increment if we were in the chemistry grid
            IF ( State_Met%InChemGrid(I,J,L) ) THEN
               COUNT_CHEM3D(X,Y,K) = COUNT_CHEM3D(X,Y,K) + 1
            ENDIF

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !-----------------------
      ! Accumulate quantities
      !-----------------------

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( W, N, X, Y, K, I, J, L, TMP, H, R, SCALEAODnm )
!$OMP+SCHEDULE( DYNAMIC )
      DO W = 1, Input_Opt%N_ND50

         ! ND50 Tracer number
         N = Input_Opt%ND50_TRACERS(W)

         ! Loop over levels
         DO K = 1, ND50_NL
            L = LOFF + K

         ! Loop over latitudes 
         DO Y = 1, ND50_NJ
            J = JOFF + Y

         ! Loop over longitudes
         DO X = 1, ND50_NI
            I = GET_I( X )

            ! Archive by simulation 
            IF ( N <= nAdvect ) THEN

               !--------------------------------------
               ! GEOS-CHEM advected species [v/v]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( Spc(I,J,L,N) * 
     &                      ( AIRMW 
     &                      / State_Chm%SpcData(N)%Info%emMW_g ) )

! NOTE: We can restore this once we figure out what the proper unit
! conversion is (bmy, 11/20/17)
!            ELSE IF ( N == 501 .and. IS_FULLCHEM ) THEN
!
!               !--------------------------------------
!               ! OH CONCENTRATION [molec/cm3]
!               ! NOTE: Only archive at chem timestep
!               !--------------------------------------
!!               Q(X,Y,K,W) = Q(X,Y,K,W) + SAVEOH(I,J,L)
!               Q(X,Y,K,W) = Q(X,Y,K,W) + Spc(I,J,L,id_OH)

            ELSE IF ( N == 502 .and. IS_NOy ) THEN

               !--------------------------------------
               ! NOy CONCENTRATION [v/v]
               !--------------------------------------
  
               ! Temp variable for accumulation
               TMP = 0e+0_fp
            
               ! NO
               TMP = TMP + ( ( AIRMW 
     &               / State_Chm%SpcData(id_NO)%Info%emMW_g )         
     &                       * Spc(I,J,L,id_NO) )  

               ! NO2
               TMP = TMP + ( ( AIRMW 
     &                / State_Chm%SpcData(id_NO2)%Info%emMW_g )         
     &                       * Spc(I,J,L,id_NO2) )
               ! PAN
               TMP = TMP + ( ( AIRMW 
     &                / State_Chm%SpcData(id_PAN)%Info%emMW_g )         
     &                       * Spc(I,J,L,id_PAN) )

               ! HNO3
               TMP = TMP + ( ( AIRMW 
     &                / State_Chm%SpcData(id_HNO3)%Info%emMW_g )        
     &                       * Spc(I,J,L,id_HNO3) )
            
               ! NPMN
               TMP = TMP + ( ( AIRMW 
     &                / State_Chm%SpcData(id_NPMN)%Info%emMW_g )         
     &                       * Spc(I,J,L,id_NPMN) )
               ! IPMN
               TMP = TMP + ( ( AIRMW 
     &                / State_Chm%SpcData(id_IPMN)%Info%emMW_g )         
     &                       * Spc(I,J,L,id_IPMN) )

               ! PPN
               TMP = TMP + ( ( AIRMW 
     &                / State_Chm%SpcData(id_PPN)%Info%emMW_g )         
     &                       * Spc(I,J,L,id_PPN) )
 
               ! R4N2
               TMP = TMP + ( ( AIRMW 
     &                / State_Chm%SpcData(id_R4N2)%Info%emMW_g )        
     &                       * Spc(I,J,L,id_R4N2) )
            
               ! N2O5
               TMP = TMP + ( 2d0 * ( AIRMW 
     &                / State_Chm%SpcData(id_N2O5)%Info%emMW_g )  
     &                       * Spc(I,J,L,id_N2O5) )
                        
               ! HNO4
               TMP = TMP + ( ( AIRMW 
     &                / State_Chm%SpcData(id_HNO4)%Info%emMW_g )       
     &                       * Spc(I,J,L,id_HNO4) )
 
               ! Accumulate into Q
               Q(X,Y,K,W) = Q(X,Y,K,W) + TMP

            ELSE IF ( N == 503 ) THEN

               !--------------------------------------
               ! RELATIVE HUMIDITY [%]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + State_Met%RH(I,J,L)

            ELSE IF ( N == 504 ) THEN

               !--------------------------------------
               ! 3-D CLOUD FRACTIONS [unitless]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + State_Met%CLDF(I,J,L)

            ELSE IF ( N == 505 .and. IS_OPTD ) THEN

               !--------------------------------------
               ! COLUMN OPTICAL DEPTHS [unitless]
               !--------------------------------------
               Q(X,Y,1,W) = Q(X,Y,1,W) + State_Met%OPTD(I,J,L)

            ELSE IF ( N == 506 .and. IS_CLDTOPS ) THEN

               !--------------------------------------
               ! CLOUD TOP HEIGHTS [mb]
               !--------------------------------------
               IF ( K == 1 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                       State_Met%PEDGE(I,J,State_Met%CLDTOPS(I,J))
               ENDIF
              
            ELSE IF ( N == 507 ) THEN

               !--------------------------------------
               ! AIR DENSITY [molec/cm3] 
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%AIRDEN(I,J,L) *
     &                      XNUMOLAIR  * 1d-6 )

            ELSE IF ( N == 508 .and. IS_SEASALT ) THEN

               !--------------------------------------
               ! TOTAL SEASALT TRACER [v/v]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                      ( Spc(I,J,L,id_SALA) + 
     &                        Spc(I,J,L,id_SALC) ) *
     &                        ( AIRMW 
     &                       / State_Chm%SpcData(id_SALA)%Info%emMW_g )

            ELSE IF ( N == 509 ) THEN

               !--------------------------------------
               ! PBL HEIGHTS [m] 
               !--------------------------------------
               IF ( K == 1 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) + GET_PBL_TOP_m( I, J )
               ENDIF
 
            ELSE IF ( N == 510 ) THEN

               !--------------------------------------
               ! PBL HEIGHTS [levels] 
               !--------------------------------------
               IF ( K == 1 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) + GET_PBL_TOP_L( I, J )
               ENDIF

            ELSE IF ( N == 511 ) THEN

               !--------------------------------------
               ! GRID BOX HEIGHTS [m]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + State_Met%BXHEIGHT(I,J,L)

            ELSE IF ( N == 512 ) THEN

               !--------------------------------------
               ! PEDGE-$ [hPa]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + State_Met%PEDGE( I, J, K )

            ELSE IF ( N == 513 .and. IS_SLP ) THEN

               !--------------------------------------
               ! SEA LEVEL PRESSURE [hPa]
               !--------------------------------------
               IF ( K == 1 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) + State_Met%SLP(I,J)
               ENDIF

            ELSE IF ( N == 514 ) THEN

               !--------------------------------------
               ! ZONAL (U) WIND [m/s]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + State_Met%U(I,J,L)

            ELSE IF ( N == 515 ) THEN

               !--------------------------------------
               ! MERIDIONAL (V) WIND [m/s]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + State_Met%V(I,J,L)

            ELSE IF ( N == 516 ) THEN 

               !--------------------------------------
               ! TEMPERATURE [K]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + State_Met%T(I,J,L)

            ELSE IF ( N == 517 ) THEN

               !--------------------------------------
               ! SULFATE AOD [unitless]
               ! for wavelengths set in Radiation Menu
               !
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               ISPC = 1 !sulfate
               IF ( .not. LINTERP ) THEN
                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                         ODAER(I,J,L,IWVSELECT(1,1),ISPC)
               ELSE
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(X,Y,K,W) = Q(X,Y,K,W) +
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDIF

            ELSE IF ( N == 518 ) THEN

               !--------------------------------------
               ! BLACK CARBON AOD [unitless]
               ! for wavelengths set in Radiation Menu
               !
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               ISPC = 2 !BC
               IF ( .not. LINTERP ) THEN
                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                         ODAER(I,J,L,IWVSELECT(1,1),ISPC)
               ELSE
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(X,Y,K,W) = Q(X,Y,K,W) +
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDIF

            ELSE IF ( N == 519 ) THEN

               !--------------------------------------
               ! ORG CARBON AOD [unitless]
               ! for wavelengths set in Radiation Menu
               !
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               ISPC = 3 !OC
               IF ( .not. LINTERP ) THEN
                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                        ODAER(I,J,L,IWVSELECT(1,1),ISPC)
               ELSE
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(X,Y,K,W) = Q(X,Y,K,W) +
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDIF

            ELSE IF ( N == 520 ) THEN

               !--------------------------------------
               ! ACCUM SEASALT AOD [unitless]
               ! for wavelengths set in Radiation Menu
               !
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               ISPC = 4 !SSa
               IF ( .not. LINTERP ) THEN
                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                         ODAER(I,J,L,IWVSELECT(1,1),ISPC)
               ELSE
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(X,Y,K,W) = Q(X,Y,K,W) +
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDIF

            ELSE IF ( N == 521 ) THEN

               !--------------------------------------
               ! COARSE SEASALT AOD [unitless]
               ! for wavelengths set in Radiation Menu
               !
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               ISPC = 5 !SSc
               IF ( .not. LINTERP ) THEN
                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                         ODAER(I,J,L,IWVSELECT(1,1),ISPC)
               ELSE
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(X,Y,K,W) = Q(X,Y,K,W) +
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDIF

            ELSE IF ( N == 522 ) THEN
               
               !--------------------------------------
               ! TOTAL DUST OPTD [unitless]
               ! for wavelengths set in Radiation Menu
               !
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               DO ISPC = 1, NDUST
                 
                  IF ( .not. LINTERP ) THEN
                     ! Accumulate
                     Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                            ODMDUST(I,J,L,IWVSELECT(1,1),ISPC)
                  ELSE
                   ! Interpolated using angstrom exponent between
                   ! Closest available wavelengths
                   ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                   !catch any zero values before interpolation
                   IF ((ODMDUST(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                 (ODMDUST(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(X,Y,K,W) = Q(X,Y,K,W) +
     &              (ODMDUST(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODMDUST(I,J,L,IWVSELECT(1,1),ISPC)/
     &                             ODMDUST(I,J,L,IWVSELECT(2,1),ISPC))))
                   ENDIF
                  ENDIF

               ENDDO

            ELSE IF ( ( N >= 523 ) .and. ( N <= 529) ) THEN

               !--------------------------------------
               ! DUST OPTD BINS1-7 [unitless]
               ! for wavelengths set in Radiation Menu
               !
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               ISPC = N - 522

               IF ( .not. LINTERP ) THEN
                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                         ODMDUST(I,J,L,IWVSELECT(1,1),ISPC)
               ELSE
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODMDUST(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODMDUST(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(X,Y,K,W) = Q(X,Y,K,W) +
     &              (ODMDUST(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODMDUST(I,J,L,IWVSELECT(1,1),ISPC)/
     &                             ODMDUST(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDIF

            ELSE IF ( N == 530 ) THEN

               !--------------------------------------
               ! PM 2.5, ug/m3 at STP
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + PM25(I,J,L)
               
            ELSE

               ! Skip other tracers
               CYCLE

            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDDO 
!$OMP END PARALLEL DO

      ! Free pointers
      Spc => NULL()

#endif

      END SUBROUTINE ACCUMULATE_DIAG50
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_diag50
!
! !DESCRIPTION: Subroutine WRITE\_DIAG50 computes the 24-hr time-average of 
!  quantities and saves to bpch file format.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE_DIAG50( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      USE BPCH2_MOD,          ONLY : BPCH2
      USE BPCH2_MOD,          ONLY : GET_MODELNAME
      USE BPCH2_MOD,          ONLY : GET_HALFPOLAR
      USE BPCH2_MOD,          ONLY : OPEN_BPCH2_FOR_WRITE
      USE CMN_SIZE_MOD             ! Size Parameters
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GC_GRID_MOD,        ONLY : GET_XOFFSET
      USE GC_GRID_MOD,        ONLY : GET_YOFFSET
      USE inquireMod,         ONLY : findFreeLUN
      USE TIME_MOD,           ONLY : EXPAND_DATE
      USE TIME_MOD,           ONLY : GET_NYMD_DIAG
      USE TIME_MOD,           ONLY : GET_NHMS
      USE TIME_MOD,           ONLY : GET_TAU
      USE TIME_MOD,           ONLY : GET_TS_DYN
      USE TIME_MOD,           ONLY : TIMESTAMP_STRING

#if defined( USE_HDF5 )
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%%% NOTE: THIS IS PROBABLY BROKEN, WE WILL NOT REPLACE IT
      !%%%%% BECAUSE WE ARE GOING TO MOVE TO netCDF OUTPUT (bmy, 4/11/18)
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Only include this if we are linking to HDF5 library (bmy, 12/21/09)
      USE HDF_MOD,            ONLY : OPEN_HDF
      USE HDF_MOD,            ONLY : CLOSE_HDF
      USE HDF_MOD,            ONLY : WRITE_HDF
      USE HDF5,               ONLY : HID_T
      INTEGER(HID_T)              :: IU_ND50_HDF
#endif

#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Rewrote to remove hardwiring and for better efficiency.  Added extra
!        diagnostics and updated numbering scheme. (bmy, 7/20/04)
!  (2 ) Now only archive NO, NO2, OH, O3 on every chemistry timestep (i.e. 
!        only when fullchem is called).  Also remove reference to FIRST. 
!        (bmy, 10/25/04)
!  (3 ) Now divide tracers 82-87 (i.e. various AODs) by GOOD_CT_CHEM since
!        these are only updated once per chemistry timestep (bmy, 1/14/05)
!  (4 ) Now save grid box heights as tracer #93.  Now save 3-D cloud fraction
!        as tracer #79. (bmy, 4/20/05)
!  (5 ) Remove references to TRCOFFSET because it is always zero (bmy, 6/24/05)
!  (6 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (7 ) DIVISOR is now a 3-D array.  Now zero COUNT_CHEM3D.  Now zero Q
!        array with array assignment statement. (phs, 1/24/07)
!  (8 ) RH should be tracer #17 under "TIME-SER" category (bmy, 2/11/08)
!  (9 ) Bug fix: replace "PS-PTOP" with "PEDGE-$" (bmy, 10/7/08)
!  (10) Change timestamp for filename.  Now save SLP under tracer #18 in 
!        "DAO-FLDS".  Also set unit to 'K' for temperature field. 
!        (ccc, tai, bmy, 10/13/09)
!  (11) Now have the option of saving out to HDF5 format.  NOTE: we have to
!        bracket HDF-specific code with an #ifdef statement to avoid problems
!        if the HDF5 libraries are not installed. (amv, bmy, 12/21/09)
!  12 Nov 2010 - R. Yantosca - Now save out PEDGE-$ (pressure at level edges)
!                              rather than Psurface - PTOP
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC args
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      INTEGER            :: DIVISOR(ND50_NI,ND50_NJ,ND50_NL)
      INTEGER            :: I,    J,     L,   W, N, nAdvect
      INTEGER            :: GMNL, GMTRC, IOS, X, Y, K, NHMS
      CHARACTER(LEN=16)  :: STAMP
      CHARACTER(LEN=40)  :: CATEGORY
      CHARACTER(LEN=40)  :: UNIT
      CHARACTER(LEN=255) :: FILENAME
#endif

      !=================================================================
      ! WRITE_DIAG50 begins here!
      !=================================================================

      ! Assume success
      RC        =  GC_SUCCESS

#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )

      ! Nubmer of advected species
      nAdvect   = State_Chm%nAdvect

      ! Replace time & date tokens in the filename
      FILENAME = Input_Opt%ND50_FILE

      ! Change to get the good timestamp: day that was run and not next 
      ! day if saved at midnight
      NHMS = GET_NHMS()
      IF ( NHMS == 0 ) NHMS = 240000

      CALL EXPAND_DATE( FILENAME, GET_NYMD_DIAG(), NHMS )
      
      ! Find a free file LUN
      IU_ND50 = findFreeLun()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME ), IU_ND50
 100  FORMAT( '     - DIAG50: Opening file ', a, ' on unit ', i4 ) 

      ! Open output file
      IF ( Input_Opt%LND50_HDF ) THEN
#if   defined( USE_HDF5 )
         ! Only include this if we are linking to HDF5 library (bmy, 12/21/09)
         CALL OPEN_HDF( IU_ND50_HDF,         FILENAME,
     &                  Input_Opt%ND50_IMAX, Input_Opt%ND50_IMIN,   
     &                  Input_Opt%ND50_JMAX, Input_Opt%ND50_JMIN,
     &                  ND50_NI,             ND50_NJ )
#endif
      ELSE
         CALL OPEN_BPCH2_FOR_WRITE( IU_ND50, FILENAME, TITLE )
      ENDIF

      ! Set ENDING TAU for this bpch write 
      TAU1 = GET_TAU()

      !=================================================================
      ! Compute 24-hr average quantities for bpch file output
      !=================================================================

      ! Echo info
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 110 ) STAMP
 110  FORMAT( '     - DIAG50: Saving to disk at ', a ) 

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( X, Y, K, W, DIVISOR )
      DO W = 1, Input_Opt%N_ND50

         ! ND50 Tracer number
         N = Input_Opt%ND50_TRACERS(W)
         
         ! Pick the proper divisor
         ! OH is defined only in boxes with active chemistry
         ! Others are defined everywhere
         IF ( N == 501 ) THEN
            DIVISOR = COUNT_CHEM3D
         ELSE
            DIVISOR = COUNT
         ENDIF

         ! Loop over grid boxes
         DO K = 1, ND50_NL
         DO Y = 1, ND50_NJ
         DO X = 1, ND50_NI

            ! Avoid division by zero
            IF ( DIVISOR(X,Y,K) > 0 ) THEN
               Q(X,Y,K,W) = Q(X,Y,K,W) / DBLE( DIVISOR(X,Y,K) )
            ELSE
               Q(X,Y,K,W) = 0e+0_fp
            ENDIF

         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      !=================================================================
      ! Write each tracer from "timeseries.dat" to the timeseries file
      !=================================================================
      DO W = 1, Input_Opt%N_ND50

         ! ND50 tracer number
         N = Input_Opt%ND50_TRACERS(W)

         ! Save by simulation
         IF ( N <= nAdvect ) THEN

            !---------------------
            ! GEOS-CHEM species
            !---------------------
            CATEGORY = 'IJ-AVG-$'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND50_NL
            GMTRC    = N

         ELSE IF ( N == 501 ) THEN

            !---------------------
            ! OH 
            !---------------------
            CATEGORY = 'CHEM-L=$'
            UNIT     = 'molec/cm3'
            GMNL     = ND50_NL
            GMTRC    = 1

         ELSE IF ( N == 502 ) THEN

            !---------------------
            ! NOy 
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND50_NL
            GMTRC    = 2

         ELSE IF ( N == 503 ) THEN

            !---------------------
            ! Relative humidity
            !---------------------            
            CATEGORY = 'TIME-SER'
            UNIT     = '%'
            GMNL     = ND50_NL
            GMTRC    = 3

         ELSE IF ( N == 504 ) THEN

            !---------------------
            ! 3-D Cloud fractions
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'unitless'
            GMNL     = ND50_NL
            GMTRC    = 4

         ELSE IF ( N == 505 ) THEN

            !---------------------
            ! Column opt depths 
            !---------------------
            CATEGORY  = 'TIME-SER'
            UNIT      = 'unitless'
            GMNL      = 1
            GMTRC     = 5
            
         ELSE IF ( N == 506 ) THEN
        
            !---------------------
            ! Cloud top heights 
            !---------------------
            CATEGORY  = 'TIME-SER'
            UNIT      = 'hPa'
            GMNL      = 1
            GMTRC     = 6

         ELSE IF ( N == 507 ) THEN

            !---------------------
            ! Air Density 
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'molec/cm3'
            GMNL     = ND50_NL
            GMTRC    = 7

         ELSE IF ( N == 508 ) THEN

            !----------------------
            ! Total Seasalt
            !----------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND50_NL
            GMTRC    = 8

         ELSE IF ( N == 509 ) THEN 

            !---------------------
            ! PBL Height [m] 
            !---------------------
            CATEGORY = 'PBLDEPTH'
            UNIT     = 'm'
            GMNL     = 1
            GMTRC    = 1

         ELSE IF ( N == 510 ) THEN

            !---------------------
            ! PBL Height [levels]
            !---------------------
            CATEGORY = 'PBLDEPTH'
            UNIT     = 'levels'
            GMNL     = 1
            GMTRC    = 2

         ELSE IF ( N == 511 ) THEN

            !---------------------
            ! Grid box heights
            !---------------------            
            CATEGORY = 'BXHGHT-$'
            UNIT     = 'm'
            GMNL     = ND50_NL
            GMTRC    = 1

        ELSE IF ( N == 512 ) THEN

            !---------------------
            ! PEDGE-$ [hPa] 
            !---------------------
            CATEGORY  = 'PEDGE-$'
            UNIT      = 'hPa'
            GMNL      = ND50_NL
            GMTRC     = 1

          ELSE IF ( N == 513 ) THEN

            !---------------------
            ! Sea level prs [hPa]
            !---------------------            
            CATEGORY = 'DAO-FLDS'
            UNIT     = 'hPa'
            GMNL     = 1
            GMTRC    = 18

         ELSE IF ( N == 514 ) THEN

            !---------------------
            ! U-wind [m/s]
            !---------------------            
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'm/s'
            GMNL     = ND50_NL
            GMTRC    = 1

         ELSE IF ( N == 515 ) THEN

            !---------------------
            ! V-wind [m/s]
            !---------------------
            CATEGORY  = 'DAO-3D-$'
            UNIT      = 'm/s'
            GMNL      = ND50_NL
            GMTRC     = 2

          ELSE IF ( N == 516 ) THEN

            !---------------------
            ! Temperature
            !---------------------
            CATEGORY  = 'DAO-3D-$'
            UNIT      = 'K'
            GMNL      = ND50_NL
            GMTRC     = 3

         ELSE IF ( N == 517 ) THEN

            !---------------------
            ! Sulfate AOD
            !---------------------            
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMNL      = ND50_NL
            GMTRC     = 6

         ELSE IF ( N == 518 ) THEN

            !---------------------
            ! Black Carbon AOD
            !---------------------            
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMNL      = ND50_NL
            GMTRC     = 9

         ELSE IF ( N == 519 ) THEN

            !---------------------
            ! Organic Carbon AOD
            !---------------------            
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMNL      = ND50_NL
            GMTRC     = 12
            
         ELSE IF ( N == 520 ) THEN

            !---------------------
            ! SS Accum AOD
            !---------------------            
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMNL      = ND50_NL
            GMTRC     = 15

         ELSE IF ( N == 521 ) THEN

            !---------------------
            ! SS Coarse AOD
            !---------------------            
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMNL      = ND50_NL
            GMTRC     = 18

         ELSE IF ( N == 522 ) THEN

            !---------------------
            ! Total dust OD
            !---------------------   
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMNL      = ND50_NL
            GMTRC     = 4

         ELSE IF ( ( N >= 523 ) .and. ( N <= 529) ) THEN

            !---------------------
            ! Dust OD (bins 1-7)
            !---------------------
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMNL      = ND50_NL
            GMTRC     = 21+(N-523)

         ELSE IF ( N == 530 ) THEN

            !---------------------
            ! PM2.5
            !---------------------   
            CATEGORY  = 'IJ-SOA-$'
            UNIT      = 'ug/m3'
            GMNL      = ND50_NL
            GMTRC     = 18

         ELSE

            ! Otherwise skip
            CYCLE

         ENDIF

         !------------------------
         ! Save to bpch file 
         ! or HDF5 file
         !------------------------
         IF ( Input_Opt%LND50_HDF ) THEN
#if   defined( USE_HDF5 )
            ! Only include this if we are linking to HDF5 library 
            ! (bmy, 12/21/09)
            CALL WRITE_HDF( IU_ND50_HDF,  N,
     &                      CATEGORY,     GMTRC,        UNIT,
     &                      TAU0,         TAU1,         RESERVED,
     &                      ND50_NI,      ND50_NJ,      GMNL,
     &                      Input_Opt%ND50_IMIN+I0,
     &                      Input_Opt%ND50_JMIN+J0,
     &                      Input_Opt%ND50_LMIN,
     &                      REAL( Q(1:ND50_NI, 1:ND50_NJ, 1:GMNL, W)))
#endif
         ELSE
            CALL BPCH2( IU_ND50,      MODELNAME,    LONRES,   
     &                  LATRES,       HALFPOLAR,    CENTER180, 
     &                  CATEGORY,     GMTRC,        UNIT,      
     &                  TAU0,         TAU1,         RESERVED,  
     &                  ND50_NI,      ND50_NJ,      GMNL,     
     &                  Input_Opt%ND50_IMIN+I0,
     &                  Input_Opt%ND50_JMIN+J0,
     &                  Input_Opt%ND50_LMIN, 
     &                  REAL( Q(1:ND50_NI, 1:ND50_NJ, 1:GMNL, W),4))
         ENDIF
      ENDDO

      ! Echo info
      WRITE( 6, 120 ) TRIM( FILENAME )
 120  FORMAT( '     - DIAG50: Closing file ', a )

      ! Close file
      IF ( Input_Opt%LND50_HDF ) THEN
#if   defined( USE_HDF5 )
         ! Only include this if we are linking to HDF5 library (bmy, 12/21/09)
         CALL CLOSE_HDF( IU_ND50_HDF )
#endif
      ELSE 
         CLOSE( IU_ND50 )
      ENDIF

      !=================================================================
      ! Re-initialize quantities for the next diagnostic cycle
      !=================================================================

      ! Echo info
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 130 ) STAMP
 130  FORMAT( '     - DIAG50: Zeroing arrays at ', a )

      ! Set STARTING TAU for the next bpch write
      TAU0       = GET_TAU() + ( GET_TS_DYN() / 3600e+0_fp )

      ! Zero counters
      COUNT        = 0
      COUNT_CHEM3D = 0  
      
      ! Zero accumulating array
      Q            = 0e+0_fp
#endif

      END SUBROUTINE WRITE_DIAG50
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_i
!
! !DESCRIPTION: Function GET\_I returns the absolute longitude index (I), 
!  given the relative longitude index (X).
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_I( X ) RESULT( I )
!
! !USES:
!
      USE CMN_SIZE_MOD         ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: X   ! Relative longitude index
!
! !RETURN VALUE:
!
      INTEGER             :: I   ! Absolute longitude index
!
! !REMARKS:
! 
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! GET_I begins here!
      !=================================================================
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      ! Add the offset to X to get I  
      I = IOFF + X

      ! Handle wrapping around the date line, if necessary
      IF ( I > IIPAR ) I = I - IIPAR
#endif

      END FUNCTION GET_I
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_diag50
!
! !DESCRIPTION: Subroutine INIT\_DIAG50 allocates and zeroes all module arrays.
!  It also gets values for module variables from "input\_mod.f". 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_DIAG50( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      USE BPCH2_MOD,     ONLY : GET_MODELNAME
      USE BPCH2_MOD,     ONLY : GET_HALFPOLAR
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE ERROR_MOD,     ONLY : ALLOC_ERR
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE GC_GRID_MOD,   ONLY : GET_XOFFSET
      USE GC_GRID_MOD,   ONLY : GET_YOFFSET
      USE GC_GRID_MOD,   ONLY : ITS_A_NESTED_GRID
      USE TIME_MOD,      ONLY : GET_TAUb
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
!
! !REMARKS:
! 
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now get I0 and J0 correctly for nested grid simulations (bmy, 11/9/04)
!  (2 ) Now call GET_HALFPOLAR from "bpch2_mod.f" to get the HALFPOLAR flag 
!        value for GEOS or GCAP grids. (bmy, 6/28/05)
!  (3 ) Now allow ND50_IMIN to be equal to ND50_IMAX and ND50_JMIN to be
!        equal to ND50_JMAX.  This will allow us to save out longitude
!        or latitude transects.  Now allocate COUNT_CHEM3D array.
!        (cdh, phs, 1/24/07)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      INTEGER            :: AS
      CHARACTER(LEN=255) :: LOCATION
#endif

      !=================================================================
      ! INIT_DIAG50 begins here!
      !=================================================================

      ! Assume success
      RC = GC_SUCCESS

#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )

      ! Initialize
      LOCATION               = 'INIT_DIAG50 ("diag50_mod.f")'

      ! Exit if ND50 is turned off 
      IF ( .not. Input_Opt%DO_ND50 ) RETURN

      !=================================================================
      ! Error check longitude, latitude, altitude limits
      !=================================================================

      ! Get grid offsets
      IF ( ITS_A_NESTED_GRID() ) THEN
         I0 = GET_XOFFSET()
         J0 = GET_YOFFSET() 
      ELSE
         I0 = GET_XOFFSET( GLOBAL=.TRUE. )
         J0 = GET_YOFFSET( GLOBAL=.TRUE. ) 
      ENDIF

      !-----------
      ! Longitude
      !-----------

      ! Error check ND50_IMIN
      IF ( Input_Opt%ND50_IMIN+I0 < 1       .or.
     &     Input_Opt%ND50_IMIN+I0 > IIPAR ) THEN
         CALL ERROR_STOP( 'Bad ND50_IMIN value!', LOCATION )
      ENDIF

      ! Error check ND50_IMAX
      IF ( Input_Opt%ND50_IMAX+I0 < 1      .or.
     &     Input_Opt%ND50_IMAX+I0 > IIPAR ) THEN
         CALL ERROR_STOP( 'Bad ND50_IMAX value!', LOCATION )
      ENDIF

      ! Compute longitude limits to write to disk 
      ! Also handle wrapping around the date line
      IF ( Input_Opt%ND50_IMAX >= Input_Opt%ND50_IMIN ) THEN
         ND50_NI = ( Input_Opt%ND50_IMAX - Input_Opt%ND50_IMIN ) + 1
      ELSE 
         ND50_NI = ( IIPAR - Input_Opt%ND50_IMIN ) + 1 +
     &               Input_Opt%ND50_IMAX
         WRITE( 6, '(a)' ) 'We are wrapping around the date line!'
      ENDIF

      ! Make sure that ND50_NI <= IIPAR
      IF ( ND50_NI > IIPAR ) THEN
         CALL ERROR_STOP( 'Too many longitudes!', LOCATION )
      ENDIF

      !-----------
      ! Latitude
      !-----------
      
      ! Error check JMIN_AREA
      IF ( Input_Opt%ND50_JMIN+J0 < 1       .or.
     &     Input_Opt%ND50_JMIN+J0 > JJPAR ) THEN
         CALL ERROR_STOP( 'Bad ND50_JMIN value!', LOCATION )
      ENDIF
     
      ! Error check JMAX_AREA
      IF ( Input_Opt%ND50_JMAX+J0 < 1       .or.
     &     Input_Opt%ND50_JMAX+J0 > JJPAR ) THEN
         CALL ERROR_STOP( 'Bad ND50_JMAX value!', LOCATION )
      ENDIF

      ! Compute latitude limits to write to disk (bey, bmy, 3/16/99)
      IF ( Input_Opt%ND50_JMAX >= Input_Opt%ND50_JMIN ) THEN
         ND50_NJ = ( Input_Opt%ND50_JMAX - Input_Opt%ND50_JMIN ) + 1
      ELSE
         CALL ERROR_STOP( 'ND50_JMAX < ND50_JMIN!', LOCATION )
      ENDIF     
  
      !-----------
      ! Altitude
      !-----------

      ! Error check ND50_LMIN, ND50_LMAX
      IF ( Input_Opt%ND50_LMIN < 1       .or.
     &     Input_Opt%ND50_LMAX > LLPAR ) THEN 
         CALL ERROR_STOP( 'Bad ND50 altitude values!', LOCATION )
      ENDIF

      ! # of levels to save in ND50 timeseries
      IF ( Input_Opt%ND50_LMAX >= Input_Opt%ND50_LMIN ) THEN  
         ND50_NL = ( Input_Opt%ND50_LMAX - Input_Opt%ND50_LMIN ) + 1
      ELSE
         CALL ERROR_STOP( 'ND50_LMAX < ND50_LMIN!', LOCATION )
      ENDIF

      !-----------
      ! Offsets
      !-----------
      IOFF      = Input_Opt%ND50_IMIN - 1
      JOFF      = Input_Opt%ND50_JMIN - 1
      LOFF      = Input_Opt%ND50_LMIN - 1

      !------------
      ! Counter
      !------------
      COUNT     = 0

      !------------
      ! For bpch
      !------------
      TAU0      = GET_TAUb()
      TITLE     = 'GEOS-CHEM DIAG50 24hr average time series'
      LONRES    = DISIZE
      LATRES    = DJSIZE
      MODELNAME = GET_MODELNAME()
      HALFPOLAR = GET_HALFPOLAR()

      ! Reset offsets to global values for bpch write
      I0        = GET_XOFFSET( GLOBAL=.TRUE. )
      J0        = GET_YOFFSET( GLOBAL=.TRUE. ) 

      !=================================================================
      ! Allocate arrays
      !=================================================================

      ! Accumulator array
      ALLOCATE( Q( ND50_NI, ND50_NJ, ND50_NL, Input_Opt%N_ND50),STAT=AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Q' )
      Q = 0e+0_fp

      ! 3-D full chemistry timestep counter in troposphere
      ALLOCATE( COUNT_CHEM3D( ND50_NI, ND50_NJ, ND50_NL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COUNT_CHEM3D' )
      COUNT_CHEM3D = 0e+0_fp
#endif

      END SUBROUTINE INIT_DIAG50
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_diag50
!
! !DESCRIPTION: Subroutine CLEANUP\_DIAG50 deallocates all module arrays. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_DIAG50()
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now deallocate COUNT_CHEM3D (phs, 1/24/07)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      !=================================================================
      ! CLEANUP_DIAG50 begins here!
      !=================================================================
      IF ( ALLOCATED( Q            ) ) DEALLOCATE( Q            )
      IF ( ALLOCATED( COUNT_CHEM3D ) ) DEALLOCATE( COUNT_CHEM3D )
#endif

      END SUBROUTINE CLEANUP_DIAG50
!EOC
      END MODULE DIAG50_MOD
