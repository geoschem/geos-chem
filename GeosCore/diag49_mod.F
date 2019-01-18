!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag49_mod.F
!
! !DESCRIPTION: Module DIAG49\_MOD contains variables and routines to save 
!  out 3-D instantaneous timeseries output to disk.
!\\
!\\
! !INTERFACE: 
!
      MODULE DIAG49_MOD
!
! !USES:
!
      USE inquireMod,    ONLY : findFreeLUN ! (ewl, 6/16/15)
      USE PhysConstants
      USE PRECISION_MOD       ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC  :: DIAG49
      PUBLIC  :: ITS_TIME_FOR_DIAG49
      PUBLIC  :: INIT_DIAG49
!
! !PRIVATE MEMBER FUNCTIONS:
! 
      PRIVATE :: ITS_TIME_TO_CLOSE_FILE
      PRIVATE :: GET_I
!
! !REMARKS:
!  NOTE by Melissa Sulprizio, 26 May 2015
!  ----------------------------------------------------------------------------
!  The emission options in the timeseries diagnostics were removed 
!  from v10-01 since HEMCO now handles the emission diagnostics. To 
!  utilize the diagnostics capability from HEMCO and output hourly 
!  isoprene emissions, you can follow these steps:
!
!    1. At the top of your HEMCO_Config.rc file, set DiagnFreq to Hourly 
!       and add a line for DiagnFile:
!
!          DiagnPrefix: HEMCO_Diagnostics
!          DiagnFreq: Hourly
!          DiagnFile: DiagnFile.rc
!
!    2.  Create a new text file in your run directory named DiagnFile.rc 
!        and list the emission fields that you would like to be saved out. 
!        For example:
!
!          # Name        Spec ExtNr Cat Hier Dim OutUnit
!          ISOP_BIOG     ISOP 108    1   1   2   kg/m2/s
!
!        NOTE: The ExtNr, Cat, Hier, and Dim values listed above were 
!        obtained from the MEGAN entries in the HEMCO_Config.rc file.
!
!    3. You can then run GEOS-Chem as usual. HEMCO will write out the 
!       specified diagnostics in a netCDF file named 
!       HEMCO_Diagnostics.YYYYMMDDHHmm.nc. I recommend running a short 
!       1-day simulation to make sure the diagnostic output is what 
!       you expect.
!
!  For more details on the HEMCO diagnostics, please see this post 
!  in the HEMCO User’s Guide:
!
!       http://wiki.geos-chem.org/The_HEMCO_User%27s_Guide#Diagnostics
!
!  ND49 tracer numbers:
!  ============================================================================
!  1 - nAdvect   : GEOS-CHEM advected species               [v/v        ]
!  501           : OH concentration                         [molec/cm3  ]
!  502           : NOy concentration                        [v/v        ]
!  503           : Relative Humidity                        [%          ]
!  504           : 3-D Cloud fractions                      [unitless   ]
!  505           : Column optical depths                    [unitless   ]
!  506           : Cloud top heights                        [hPa        ]
!  507           : Air density                              [molec/cm3  ]
!  508           : Total seasalt tracer concentration       [unitless   ]
!  509           : PBL heights                              [m          ]
!  510           : PBL heights                              [levels     ]
!  511           : Grid box heights                         [m          ]
!  512           : PEDGE-$ (Pressure @ level edges          [hPa        ]
!  513           : Sea level pressure                       [hPa        ]
!  514           : Zonal wind (a.k.a. U-wind)               [m/s        ]
!  515           : Meridional wind (a.k.a. V-wind)          [m/s        ]
!  516           : Temperature                              [K          ]
!  517           : Sulfate aerosol optical depth            [unitless   ]
!  518           : Black carbon aerosol optical depth       [unitless   ]
!  519           : Organic carbon aerosol optical depth     [unitless   ]
!  520           : Accumulation mode seasalt optical depth  [unitless   ]
!  521           : Coarse mode seasalt optical depth        [unitless   ]
!  522           : Total dust optical depth                 [unitless   ]
!  523-529       : Size resolved dust optical depth         [unitless   ]
!  530           : PAR direct                               [hPa        ]
!  531           : PAR diffuse                              [hPa        ]
!  532           : Daily LAI                                [hPa        ]
!  533           : Temperature at 2m                        [K          ]
!  534           : Accumulated precip.                      [mm/day     ]
!  535           : Cosine of Solar zenith angle             [unitless   ]
!  536           : Total organic aerosol                    [ug/m3      ]
!  537           : ISORROPIA aerosol pH                     [unitless   ]
!  538           : ISORROPIA aerosol water                  [ug/m3 air  ]
!  539           : Ozone at 10m                             [ppb        ]
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Bug fix: get I0, J0 properly for nested grids (bmy, 11/9/04)
!  (2 ) Now references "pbl_mix_mod.f" (bmy, 2/16/05)
!  (3 ) Now saves 3-D cld frac & grid box height (bmy, 4/20/05)
!  (4 ) Remove TRCOFFSET since it's always zero  Also now get HALFPOLAR for
!        both GCAP and GEOS grids.  (bmy, 6/28/05)
!  (5 ) Bug fix: do not save SLP if it's not allocated (bmy, 8/2/05)
!  (6 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (7 ) Now references XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (8 ) Modified INIT_DIAG49 to save out transects (cdh, bmy, 11/30/06)
!  (9 ) Bug fix: accumulate into Q(X,Y,K) for dust OD (qli, bmy, 4/30/07)
!  (10) Minor bug fixes in DIAG49 (cdh, bmy, 2/11/08)
!  (11) Bug fix: replace "PS-PTOP" with "PEDGE-$"
!  (12) Modified to archive O3, NO, NOy as tracers 89, 90, 91  (tmf, 9/26/07)
!  (13) Bug fix DIAG49 for diagnostic output of SLP (tai, bmy, 10/13/09)
!  (14) Modify AOD output to wavelength specified in jv_spec_aod.dat 
!       (clh, 05/07/10)
!  (15) Bug fix in ITS_TIME_TO_CLOSE: compare HR1 to 00 not 24. (ccc, 11/11/10)
!  (16) Now do not scale AOD output (recalculated in RDAER AND DUST_MOD)
!       (skim, 02/02/11)
!  12 Nov 2010 - R. Yantosca - Changed tracer 99 to be PEDGE-$ (pressure at
!                              level edges) instead of Psurface-PTOP.
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2013 - M. Sulprizio- Add farnesene emissions for updated SOA (H. Pye)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  19 Mar 2014 - M. Sulprizio- Updated to allow for more than 75 tracers
!  10 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  22 Jun 2016 - M. Yannetti - Replace TCVV with spec db and physical const
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
      ! I0               : Offset between global & nested grid
      ! J0               : Offset between global & nested grid
      ! IOFF             : Longitude offset
      ! JOFF             : Latitude offset
      ! LOFF             : Altitude offset
      ! ND49_NI          : Number of longitudes in DIAG51 region 
      ! ND49_NJ          : Number of latitudes  in DIAG51 region
      ! ND49_NL          : Number of levels     in DIAG51 region
      ! HALFPOLAR        : Used for bpch file output
      ! CENTER180        : Used for bpch file output
      ! LONRES           : Used for bpch file output
      ! LATRES           : Used for bpch file output
      ! MODELNAME        : Used for bpch file output
      ! RESERVED         : Used for bpch file output
      !=================================================================
      INTEGER            :: IOFF,           JOFF,    LOFF
      INTEGER            :: I0,             J0
      INTEGER            :: ND49_NI,        ND49_NJ, ND49_NL
      INTEGER            :: HALFPOLAR
      INTEGER, PARAMETER :: CENTER180=1 
      REAL*4             :: LONRES,         LATRES
      CHARACTER(LEN=20)  :: MODELNAME
      CHARACTER(LEN=40)  :: RESERVED = ''
      CHARACTER(LEN=80)  :: TITLE

      ! Logical file unit for ND49 diagnostic
      INTEGER            :: IU_ND49

      ! Species ID flags
      INTEGER            :: id_HNO3, id_HNO4, id_N2O5, id_NO
      INTEGER            :: id_PAN,  id_NPMN, id_PPN,  id_O3
      INTEGER            :: id_R4N2, id_SALA, id_SALC, id_NO2
      INTEGER            :: id_IPMN, id_OH,   id_Passive
#endif

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diag49 
!
! !DESCRIPTION: Subroutine DIAG49 produces time series (instantaneous fields) 
!  for a geographical domain from the information read in timeseries.dat.  
!  Output will be in binary punch (BPCH) format.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DIAG49( am_I_Root, Input_Opt,  State_Met,
     &                   State_Chm, State_Diag, RC         )
!
! !USES:
!
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_
      USE State_Diag_Mod,     ONLY : DgnState
      USE State_Met_Mod,      ONLY : MetState
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      USE BPCH2_MOD,          ONLY : BPCH2,   OPEN_BPCH2_FOR_WRITE
      USE CMN_FJX_MOD,        ONLY : ODAER, ODMDUST, NRH, NDUST
      USE CMN_FJX_MOD,        ONLY : IWVSELECT, ACOEF_WV, BCOEF_WV
      USE CMN_O3_MOD	           ! SAVEOH  SAVEOA
!------------------------------------------------------------------------------
! This diagnostic requires some fixes. Comment out for now (mps,3/2/18)
!#if defined ( BPCH_DIAG ) || defined( BPCH_TIMESER )
!      USE DIAG_MOD,           ONLY : AD44b, AD_RA, AD_ZL
!#endif
!------------------------------------------------------------------------------
      USE ErrCode_Mod
      USE FILE_MOD,           ONLY : IOERROR
      USE GC_GRID_MOD,        ONLY : GET_XOFFSET,        GET_YOFFSET
      USE inquireMod,         ONLY : findFreeLun
      USE PhysConstants            ! XTRA2, XNUMOLAIR
      USE PBL_MIX_MOD,        ONLY : GET_PBL_TOP_L,      GET_PBL_TOP_m
      USE TIME_MOD,           ONLY : EXPAND_DATE
      USE TIME_MOD,           ONLY : GET_NYMD,           GET_NHMS
      USE TIME_MOD,           ONLY : GET_NYMD_DIAG,      GET_TS_DIAG
      USE TIME_MOD,           ONLY : GET_TAU,            GET_HOUR
      USE TIME_MOD,           ONLY : ITS_A_NEW_DAY
      USE TIME_MOD,           ONLY : TIMESTAMP_STRING
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
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  09 Apr 1999 - I. Bey, R. Martin, R. Yantosca - Initial version
!  (1 ) Now bundled into "diag49_mod.f".  Now reference STT from 
!        "tracer_mod.f".  Now scale aerosol & dust OD's to 400 nm.  
!        (bmy, rvm, aad, 7/9/04)
!  (2 ) Updated tracer # for NO2 (bmy, 10/25/04)
!  (3 ) Remove reference to "CMN".  Also now get PBL heights in meters and 
!        model layers from GET_PBL_TOP_m and GET_PBL_TOP_L of "pbl_mix_mod.f".
!        (bmy, 2/16/05)
!  (4 ) Now reference CLDF and BXHEIGHT from "dao_mod.f".  Now save 3-D cloud 
!        fraction as tracer #79 and box height as tracer #93.  Now remove 
!        reference to PBL from "dao_mod.f"(bmy, 4/20/05)
!  (5 ) Remove references to TRCOFFSET because it is always zero (bmy, 6/24/05)
!  (6 ) Now do not save SLP data if it is not allocated (bmy, 8/2/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Now references XNUMOLAIR from "tracer_mod.f".  Bug fix: now must sum
!        aerosol OD's over all RH bins.  Also zero Q array. (bmy, 11/1/05)
!  (9 ) Bug fix: accumulate into Q(X,Y,K) for dust OD (qli, bmy, 4/30/07)
!  (10) Bug fix: UNIT should be "levels" for tracer 77.  Also RH should be
!        tracer #17 under "TIME-SER" category. (cdh, bmy, 2/11/08)
!  (11) Bug fix: replace "PS-PTOP" with "PEDGE-$" (bmy, phs, 10/7/08)
!  (12) Change the new day condition to open a new file. (ccc, 8/12/09)
!  (13) Change the timestamp for the filename when closing (ccc, 8/12/09)
!  (14) Add outputs for EMISS_BVOC (10 tracers), TS, PARDR, PARDF and ISOLAI
!        (mpb, 11/19/09)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  11 Apr 2012 - R. Yantosca - Replace lai_mod.F with modis_lai_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  14 Mar 2013 - M. Payer    - Replace NOx and Ox with NO, NO2, and O3 as part
!                              of removal of NOx-Ox partitioning
!  13 Aug 2013 - M. Sulprizio- Add farnesene emissions for updated SOA (H. Pye)
!  25 Jul 2014 - R. Yantosca - Remove reference to commsoil_mod.F90
!  06 Nov 2014 - R. Yantosca - Now use State_Met%AIRDEN(I,J,L)
!  06 Nov 2014 - R. Yantosca - Now use State_Met%CLDF(I,J,L)
!  06 Nov 2014 - R. Yantosca - Now use State_Met%OPTD(I,J,L)
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE with State_Met%PEDGE.
!                              Remove dependency on pressure_mod.
!  24 Mar 2015 - E. Lundgren - Remove dependency on tracer_mod
!  25 Mar 2015 - E. Lundgren - Change tracer units from kg to kg/kg and
!                              remove AD pointer
!  15 Apr 2015 - M. Sulprizio- Comment out emission diagnostics for now. These
!                              are now handled by the HEMCO.
!  24 Jun 2015 - E. Lundgren - Add code to print out mass conservation output,
!                              blocked off by #if defined( MASSCONS )
!  17 Nov 2015 - E. Lundgren - Bug fix in masscons kg calculation: use AD
!  16 Jun 2016 - K. Yu       - Now define species ID's with the Ind_ function
!  17 Jun 2016 - R. Yantosca - Only define species ID's on the first call
!  22 Jun 2016 - R. Yantosca - Now make species ID's module variables
!  22 Jun 2016 - R. Yantosca - Rename IDTPASV to id_PASV.  Also move the 
!                              initializations of the id_ flags before the
!                              first mass conservation test.
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  11 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!  03 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!  06 Feb 2018 - E. Lundgren - Change ND49 frequency unit from min to sec
!  18 Jan 2019 - M. Sulprizio- Rename PASV to Passive
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      ! SAVEd scalars
      LOGICAL, SAVE            :: FIRST  = .TRUE.
      LOGICAL, SAVE            :: IS_FULLCHEM
      LOGICAL, SAVE            :: IS_NOy,      IS_CLDTOPS, IS_OPTD
      LOGICAL, SAVE            :: IS_SEASALT,  IS_SLP,     LSOA

      ! Scalars
      LOGICAL                  :: LINTERP
      INTEGER                  :: nAdvect
      INTEGER                  :: IOS,  GMTRC, GMNL, I, J, K, L 
      INTEGER                  :: N,    R,     H,    W, X, Y
      INTEGER                  :: NYMD, NHMS,  TS_DIAG
      INTEGER                  :: ISPC
      REAL(f8)                 :: TAU
      REAL(fp)                 :: TMP, SCALEAODnm
      REAL(fp)                 :: Q( ND49_NI, ND49_NJ, ND49_NL )
      REAL(fp)                 :: TMPO3, CONV
      CHARACTER(LEN=16)        :: STAMP
      CHARACTER(LEN=40)        :: CATEGORY
      CHARACTER(LEN=40)        :: UNIT
      CHARACTER(LEN=255)       :: FILENAME

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%  MASS CONSERVATION TEST:                                  %%%
      !%%%  Define extra quantities for mass conservation printout   %%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER                  :: IU_TRCR
      REAL(fp)                 :: TOTALMASS

      ! Aerosol types (rvm, aad, bmy, 7/20/04)
      INTEGER                  :: IND(6) = (/ 22, 29, 36, 43, 50, 15 /)

      ! For fields from Input_Opt
      LOGICAL                  :: LSOILNOX
      LOGICAL                  :: DO_DIAG_WRITE
 
      ! Pointers
      REAL(fp), POINTER        :: Spc(:,:,:,:)

      ! For mass conservation tests
      LOGICAL                  :: DO_MASSCONS
#endif
      
      !=================================================================
      ! DIAG49 begins here!
      !=================================================================

#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      ! Assume success
      RC            =  GC_SUCCESS

      ! Number of advected species
      nAdvect       = State_Chm%nAdvect

      ! Copy values from Input_Opt
      LSOA          = Input_Opt%LSOA
      LSOILNOX      = Input_Opt%LSOILNOX
      DO_DIAG_WRITE = Input_Opt%DO_DIAG_WRITE

#if defined( MASSCONS )
      DO_MASSCONS   = .TRUE.
#else
      DO_MASSCONS   = .FALSE.
#endif
      
      !----------------------------------------
      ! First-time handling (bmy, 6/22/16)
      !----------------------------------------
      IF ( FIRST ) THEN

         ! Define species ID flags
         id_HNO3 = Ind_('HNO3')
         id_HNO4 = Ind_('HNO4')
         id_N2O5 = Ind_('N2O5')
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
         id_Passive = Ind_('PassiveTracer')

         ! Set logical flags on the first timestep
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

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%  MASS CONSERVATION TEST:                                  %%%
      !%%%  If we are just printing out mass conservation, we can    %%%
      !%%%  prevent writing to the timeseries bpch file.  We are     %%%
      !%%%  only interested in printing out the sum of tracers.      %%%
      !%%%  (ewl, 6/24/15)                                           %%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF ( DO_MASSCONS .or. id_Passive > 0 ) THEN
         DO_DIAG_WRITE = .FALSE.
      ENDIF

      ! Point to
      SPC => State_Chm%Species

      !=================================================================
      ! If it's a new day, open a new BPCH file and write file header
      ! We need to check if it's a new day + 1 ND49 time step (ccc, 8/12/09)
      !=================================================================
!--- Previous to (ccc, 8/12/09)
!      IF ( ITS_A_NEW_DAY() ) THEN
      NHMS    = GET_NHMS()
      TS_DIAG = Input_Opt%ND49_FREQ

      ! To change TS_DIAG to NHMS format
      TS_DIAG = TS_DIAG/3600 * 10000 + (TS_DIAG - (TS_DIAG/3600)*3600)
     &          * 100 + (TS_DIAG - (TS_DIAG/60)*60)

      IF ( NHMS == TS_DIAG ) THEN     ! It's a new day for diagnostics.

         ! Find a free file LUN
         IU_ND49 = findFreeLUN()

         ! Expand date tokens in the file name
         FILENAME = TRIM( Input_Opt%ND49_FILE )
         CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME ), IU_ND49
 100     FORMAT( '     - DIAG49: Opening file ', a, ' on unit ', i4 )
        
         IF ( DO_DIAG_WRITE ) THEN

            ! Open bpch file and write top-of-file header
            CALL OPEN_BPCH2_FOR_WRITE( IU_ND49, FILENAME, TITLE )

         ENDIF

      ENDIF

      !=================================================================
      ! Save tracers to timeseries file
      !=================================================================

      ! Echo info
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 110 ) STAMP
 110  FORMAT( '     - DIAG49: Saving timeseries at ', a )

      ! Time for BPCH file
      TAU  = GET_TAU()

      !Determine if optical properties need interpolating
      !The LUT wavelengths in IWVSELECT will match if no interpolation
      !is needed. (output is only for the first requested wavelength)
      !(DAR 10/2013)
      IF(IWVSELECT(1,1).EQ.IWVSELECT(2,1)) THEN
         LINTERP=.FALSE.
      ELSE
         LINTERP=.TRUE.
      ENDIF

      ! Loop over tracers
      DO W = 1, Input_Opt%N_ND49

         ! ND49 tracer number
         N = Input_Opt%ND49_TRACERS(W)

         ! Zero summing array
         Q = 0e+0_fp

         IF ( N <= nAdvect ) THEN

            !-------------------------------------
            ! GEOS-CHEM advected species [v/v]
            !-------------------------------------
            CATEGORY = 'IJ-AVG-$'
            UNIT     = ''           ! Let GAMAP pick the unit
            GMNL     = ND49_NL
            GMTRC    = N
            
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%  MASS CONSERVATION TEST (ewl, 6/24/15)              %%%
            !%%%  Declare a variable for total mass of tracer.       %%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            IF ( DO_MASSCONS .or. id_Passive > 0 ) THEN
               TOTALMASS = 0.0e+0_fp
            ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )

               !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               !%%%  MASS CONSERVATION TEST (ewl, 6/24/15)           %%%
               !%%%  Archive output in kg (not v/v) for printout     %%%
               !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               IF ( DO_MASSCONS .or. id_Passive > 0 ) THEN
                  Q(X,Y,K) = SPC(I,J,L,N) * State_Met%AD(I,J,L)
               ELSE
                  ! Normal ND49 output, archive in v/v
                  Q(X,Y,K) = SPC(I,J,L,N) * ( AIRMW 
     &                       / State_Chm%SpcData(N)%Info%emMW_g )
               ENDIF
               
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%  MASS CONSERVATION TEST (ewl, 6/16/15)              %%%
            !%%%  Sum the mass and write out to an ASCII file        %%%
            !%%%  Use the filename from Input_Opt%ND49_FILE          %%%
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            IF ( DO_MASSCONS .or. id_Passive > 0 ) THEN

               ! Sum the mass of tracer in kg
               TOTALMASS = SUM( Q )

               ! Find a free file LUN
               IU_TRCR   = findFreeLUN()

               ! Get the timestamp string
               STAMP     = TIMESTAMP_STRING()

               ! Print info that we are opening the file 
               WRITE( 6, 300 ) TRIM( Input_Opt%ND49_FILE ), IU_TRCR
 300           FORMAT( 'DIAG49_MOD: Opening ', a ', on unit', i4 )

               ! Open the output file in append mode
               OPEN( IU_TRCR, FILE   = TRIM( Input_Opt%ND49_FILE ),
     &                        STATUS = 'UNKNOWN', 
     &                        ACCESS = "APPEND", 
     &                        IOSTAT = IOS )

               ! Stop if we can't append to the file
               IF ( IOS > 0 ) THEN
                  CALL IOERROR( IOS, IU_TRCR, 'diag49: masscons #1' )
               ENDIF

               ! Write tracer mass, date, and time
               WRITE( IU_TRCR, *, IOSTAT=IOS ) STAMP, ' : ', TOTALMASS

               ! Stop if we can't write to the file
               IF ( IOS > 0 ) THEN
                  CALL IOERROR( IOS, IU_TRCR, 'diag49: masscons #2' )
               ENDIF

               ! Close file
               CLOSE( IU_TRCR )

            ENDIF

! NOTE: We can restore this once we figure out what the proper unit
! conversion is (bmy, 11/20/17)
!         ELSE IF ( N == 501 .and. IS_FULLCHEM ) THEN
!
!            !--------------------------------------
!            ! OH CONCENTRATION [molec/cm3]
!            !--------------------------------------
!            CATEGORY = 'TIME-SER'
!            UNIT     = 'molec/cm3'
!            GMNL     = ND49_NL
!            GMTRC    = 1
!
!            ! NOTE: OH is a non-advected species, so it will have a 
!            ! molecular weight of -1.  Need to multiply the concentration
!            ! by this to turn it into a positive. (bmy, 11/20/17)
!            TMP      = State_Chm%SpcData(id_OH)%Info%EmMw_g
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L, X, Y, K )
!            DO K = 1, ND49_NL
!               L = LOFF + K
!            DO Y = 1, ND49_NJ
!               J = JOFF + Y
!            DO X = 1, ND49_NI
!               I = GET_I( X )
!               !Q(X,Y,K) = SAVEOH(I,J,L)
!
!               ! We are in units of [kg/kg dry] here.  Follow the same
!               ! unit conversion as in routine ConvertSpc_KgKgDry_to_MND
!               ! in GeosUtil/unitconv_mod.F90 to get [molec/cm3].
!               Q(X,Y,K) = ( TMP * Spc(I,J,L,id_OH)  )
!     &                  * ( State_Met%AIRDEN(I,J,L) )
!     &                  * ( AVO / 17e-3_fp          )
!     &                  / ( 1e+6_fp                 )
!
!            ENDDO 
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!            print*, '@@@@> State_Chm%Spc_Units: ', 
!     &           State_Chm%Spc_units
!            print*, '@@@@> SPECIES(OH)     ; ', 
!     &           State_Chm%Species(23,34,1,id_OH)           
!            print*, '@@@@> Dev OH in DIAG49 (23,34,1): ', 
!     &           Q(23,34,1)
!            print*, '@@@@> Dev AIRDEN in DIAG49 (23,34,1): ', 
!     &           State_Met%AirDen(23,34,1)
!            print*, '@@@@> tmp:                : ', 
!     &           tmp

         ELSE IF ( N == 502 .and. IS_NOy ) THEN

            !--------------------------------------
            ! NOy CONCENTRATION [v/v]
            !--------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND49_NL
            GMTRC    = 2

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K, TMP )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )  

               ! Temp variable for accumulation
               TMP = 0e0_fp
            
               ! NO
               TMP = TMP + ( ( AIRMW 
     &               / State_Chm%SpcData(id_NO)%Info%emMW_g )      
     &               * SPC(I,J,L,id_NO) )

               ! NO2
               TMP = TMP + ( ( AIRMW 
     &               / State_Chm%SpcData(id_NO2)%Info%emMW_g )      
     &               * SPC(I,J,L,id_NO2) )  

               ! PAN
               TMP = TMP + ( ( AIRMW 
     &               / State_Chm%SpcData(id_PAN)%Info%emMW_g )       
     &               * SPC(I,J,L,id_PAN) )

               ! HNO3
               TMP = TMP + ( ( AIRMW
     &               / State_Chm%SpcData(id_HNO3)%Info%emMW_g )     
     &               * SPC(I,J,L,id_HNO3) )
            
               ! NPMN
               TMP = TMP + ( ( AIRMW 
     &               / State_Chm%SpcData(id_NPMN)%Info%emMW_g )      
     &               * SPC(I,J,L,id_NPMN) )

               ! IPMN
               TMP = TMP + ( ( AIRMW 
     &               / State_Chm%SpcData(id_IPMN)%Info%emMW_g )      
     &               * SPC(I,J,L,id_IPMN) )

               ! PPN
               TMP = TMP + ( ( AIRMW 
     &               / State_Chm%SpcData(id_PPN)%Info%emMW_g )      
     &               * SPC(I,J,L,id_PPN) )
 
               ! R4N2
               TMP = TMP + ( ( AIRMW 
     &               / State_Chm%SpcData(id_R4N2)%Info%emMW_g )     
     &               * SPC(I,J,L,id_R4N2) )
            
               ! N2O5
               TMP = TMP + ( ( 2e+0_fp * AIRMW 
     &               / State_Chm%SpcData(id_N2O5)%Info%emMW_g ) 
     &               * SPC(I,J,L,id_N2O5) )
                        
               ! HNO4
               TMP = TMP + ( ( AIRMW 
     &               / State_Chm%SpcData(id_HNO4)%Info%emMW_g )     
     &               * SPC(I,J,L,id_HNO4) )

               ! Save afternoon points
               Q(X,Y,K) = Q(X,Y,K) + TMP

            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 503 ) THEN

            !--------------------------------------
            ! RELATIVE HUMIDITY [%]
            !--------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = '%'
            GMNL     = ND49_NL
            GMTRC    = 3

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = State_Met%RH(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 504 ) THEN 

            !--------------------------------------
            ! 3-D CLOUD FRACTIONS [unitless]
            !--------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'unitless'
            GMNL     = ND49_NL
            GMTRC    = 4

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = State_Met%CLDF(I,J,K)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 505 .and. IS_OPTD ) THEN 

            !--------------------------------------
            ! COLUMN OPTICAL DEPTHS [unitless]
            !--------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'unitless'
            GMNL     = 1
            GMTRC    = 5

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = SUM( State_Met%OPTD(I,J,:) )
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 506 .and. IS_CLDTOPS ) THEN 

            !--------------------------------------
            ! CLOUD TOP HEIGHTS [hPa]
            !--------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'hPa'
            GMNL     = ND49_NL
            GMTRC    = 6

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = State_Met%PEDGE(I,J,State_Met%CLDTOPS(I,J))
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 507 ) THEN

            !--------------------------------------
            ! AIR DENSITY [molec/cm3]
            !--------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'molec/cm3'
            GMNL     = ND49_NL
            GMTRC    = 7

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = State_Met%AIRDEN(I,J,L) * XNUMOLAIR * 1d-6
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 508 .and. IS_SEASALT ) THEN
            
            !--------------------------------------
            ! TOTAL SEASALT TRACER [v/v]
            !--------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''        ! Let GAMAP pick unit
            GMNL     = ND49_NL
            GMTRC    = 8

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
              Q(X,Y,K) = ( SPC(I,J,L,id_SALA) + SPC(I,J,L,id_SALC) ) *
     &                     ( AIRMW
     &                       / State_Chm%SpcData(id_SALA)%Info%emMW_g )     
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 509 ) THEN

            !--------------------------------------
            ! PBL HEIGHTS [m] 
            !--------------------------------------
            CATEGORY = 'PBLDEPTH'
            UNIT     = 'm'  
            GMNL     = 1
            GMTRC    = 1

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y, TMP )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = GET_PBL_TOP_m( I, J )
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 510 ) THEN

            !--------------------------------------
            ! PBL HEIGHTS [levels] 
            !--------------------------------------
            CATEGORY = 'PBLDEPTH'
            UNIT     = 'levels'  
            GMNL     = 1
            GMTRC    = 2

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = GET_PBL_TOP_L( I, J )
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 511 ) THEN

            !--------------------------------------
            ! GRID BOX HEIGHT [m]
            !--------------------------------------
            CATEGORY = 'BXHGHT-$'
            UNIT     = 'm'
            GMNL     = ND49_NL
            GMTRC    = 1

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = State_Met%BXHEIGHT(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 512) THEN 

            !--------------------------------------
            ! PEDGE-$ (Prs @ level edges) [hPa]
            !--------------------------------------
            CATEGORY = 'PEDGE-$'
            UNIT     = 'hPa'
            GMNL     = ND49_NL
            GMTRC    = 1
                  
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = State_Met%PEDGE(I,J,K)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 513 .and. IS_SLP ) THEN

            !--------------------------------------
            ! SEA LEVEL PRESSURE [hPa]
            !--------------------------------------
            CATEGORY = 'DAO-FLDS'
            UNIT     = 'hPa'
            GMNL     = 1
            GMTRC    = 18

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = State_Met%SLP(I,J)
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 514 ) THEN

            !--------------------------------------
            ! ZONAL (U) WIND [m/s]
            !--------------------------------------
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'm/s'
            GMNL     = ND49_NL
            GMTRC    = 1

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = State_Met%U(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 515 ) THEN 

            !--------------------------------------
            ! MERIDIONAL (V) WIND [m/s]
            !--------------------------------------
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'm/s'
            GMNL     = ND49_NL
            GMTRC    = 2

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = State_Met%V(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

        ELSE IF ( N == 516 ) THEN

            !--------------------------------------
            ! TEMPERATURE [K]
            !--------------------------------------
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'K'
            GMNL     = ND49_NL
            GMTRC    = 3

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = State_Met%T(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 517 ) THEN

            !--------------------------------------
            ! SULFATE AOD [unitless]
            ! for wavelengths set in Radiation Menu
            !--------------------------------------
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND49_NL
            GMTRC    = 6

            ISPC     = 1 !sulfate
            IF ( .not. LINTERP ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
               DO K = 1, ND49_NL
                  L = LOFF + K
               DO Y = 1, ND49_NJ
                  J = JOFF + Y
               DO X = 1, ND49_NI
                  I = GET_I( X )
                  Q(X,Y,K) = Q(X,Y,K) +
     $                       ODAER(I,J,L,IWVSELECT(1,1),ISPC)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO
            ELSE
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
               DO K = 1, ND49_NL
                  L = LOFF + K
               DO Y = 1, ND49_NJ
                  J = JOFF + Y
               DO X = 1, ND49_NI
                  I = GET_I( X )
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(X,Y,K) = Q(X,Y,K) + 
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO
            ENDIF

         ELSE IF ( N == 518 ) THEN

            !--------------------------------------
            ! BLACK CARBON AOD [unitless]
            ! for wavelengths set in Radiation Menu
            !--------------------------------------
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND49_NL
            GMTRC    = 9

            ISPC     = 2 !BC
            IF ( .not. LINTERP ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
               DO K = 1, ND49_NL
                  L = LOFF + K
               DO Y = 1, ND49_NJ
                  J = JOFF + Y
               DO X = 1, ND49_NI
                  I = GET_I( X )
                  Q(X,Y,K) = Q(X,Y,K) +
     $                       ODAER(I,J,L,IWVSELECT(1,1),ISPC)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO
            ELSE
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
               DO K = 1, ND49_NL
                  L = LOFF + K
               DO Y = 1, ND49_NJ
                  J = JOFF + Y
               DO X = 1, ND49_NI
                  I = GET_I( X )
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(X,Y,K) = Q(X,Y,K) +
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO
            ENDIF

         ELSE IF ( N == 519 ) THEN

            !--------------------------------------
            ! ORGANIC CARBON AOD [unitless]
            ! for wavelengths set in Radiation Menu
            !--------------------------------------
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND49_NL
            GMTRC    = 12

            ISPC     = 3 !OC
            IF ( .not. LINTERP ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
               DO K = 1, ND49_NL
                  L = LOFF + K
               DO Y = 1, ND49_NJ
                  J = JOFF + Y
               DO X = 1, ND49_NI
                  I = GET_I( X )
                  Q(X,Y,K) = Q(X,Y,K) +
     $                       ODAER(I,J,L,IWVSELECT(1,1),ISPC)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO
            ELSE
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
               DO K = 1, ND49_NL
                  L = LOFF + K
               DO Y = 1, ND49_NJ
                  J = JOFF + Y
               DO X = 1, ND49_NI
                  I = GET_I( X )
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(X,Y,K) = Q(X,Y,K) +
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO
            ENDIF

         ELSE IF ( N == 520 ) THEN
            
            !--------------------------------------
            ! ACCUM SEASALT AOD [unitless]
            ! for wavelengths set in Radiation Menu
            !--------------------------------------
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND49_NL
            GMTRC    = 15

            ISPC     = 4 !SSa
            IF ( .not. LINTERP ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
               DO K = 1, ND49_NL
                  L = LOFF + K
               DO Y = 1, ND49_NJ
                  J = JOFF + Y
               DO X = 1, ND49_NI
                  I = GET_I( X )
                  Q(X,Y,K) = Q(X,Y,K) +
     $                       ODAER(I,J,L,IWVSELECT(1,1),ISPC)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO
            ELSE
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
               DO K = 1, ND49_NL
                  L = LOFF + K
               DO Y = 1, ND49_NJ
                  J = JOFF + Y
               DO X = 1, ND49_NI
                  I = GET_I( X )
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(X,Y,K) = Q(X,Y,K) +
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO
            ENDIF

         ELSE IF ( N == 521 ) THEN

            !--------------------------------------
            ! COARSE SEASALT AOD [unitless]
            ! for wavelengths set in Radiation Menu
            !--------------------------------------
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND49_NL
            GMTRC    = 18

            ISPC     = 5 !SSc
            IF ( .not. LINTERP ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
               DO K = 1, ND49_NL
                  L = LOFF + K
               DO Y = 1, ND49_NJ
                  J = JOFF + Y
               DO X = 1, ND49_NI
                  I = GET_I( X )
                  Q(X,Y,K) = Q(X,Y,K) +
     $                       ODAER(I,J,L,IWVSELECT(1,1),ISPC)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO
            ELSE
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
               DO K = 1, ND49_NL
                  L = LOFF + K
               DO Y = 1, ND49_NJ
                  J = JOFF + Y
               DO X = 1, ND49_NI
                  I = GET_I( X )
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(X,Y,K) = Q(X,Y,K) +
     &              (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/
     &                               ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO
            ENDIF

         ELSE IF ( N == 522 ) THEN

            !--------------------------------------
            ! TOTAL DUST OPT DEPTH [unitless]
            ! for wavelengths set in Radiation Menu
            !--------------------------------------
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMNL      = ND49_NL
            GMTRC     = 4

            DO ISPC = 1, NDUST

               IF ( .not. LINTERP ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
                  DO K = 1, ND49_NL
                     L = LOFF + K
                  DO Y = 1, ND49_NJ
                     J = JOFF + Y
                  DO X = 1, ND49_NI
                     I = GET_I( X )
                     Q(X,Y,K) = Q(X,Y,K) +
     $                          ODMDUST(I,J,L,IWVSELECT(1,1),ISPC)
                  ENDDO
                  ENDDO
                  ENDDO
!$OMP END PARALLEL DO
               ELSE
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
                  DO K = 1, ND49_NL
                     L = LOFF + K
                  DO Y = 1, ND49_NJ
                     J = JOFF + Y
                  DO X = 1, ND49_NI
                     I = GET_I( X )
                   ! Interpolated using angstrom exponent between
                   ! Closest available wavelengths
                   ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                   !catch any zero values before interpolation
                   IF ((ODMDUST(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                 (ODMDUST(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(X,Y,K) = Q(X,Y,K) +
     &              (ODMDUST(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODMDUST(I,J,L,IWVSELECT(1,1),ISPC)/
     &              ODMDUST(I,J,L,IWVSELECT(2,1),ISPC))))
                   ENDIF
                  ENDDO
                  ENDDO
                  ENDDO
!$OMP END PARALLEL DO
               ENDIF
            ENDDO

         ELSE IF ( ( N >= 523 ) .and. ( N <= 529) ) THEN

            !--------------------------------------
            ! DUST OPTD BINS1-7 [unitless]
            ! for wavelengths set in Radiation Menu
            !
            ! NOTE: Only archive at chem timestep
            !--------------------------------------
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMNL      = ND49_NL
            GMTRC     = 21+(N-523)

            ISPC = N - 522
            IF ( .not. LINTERP ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
               DO K = 1, ND49_NL
                  L = LOFF + K
               DO Y = 1, ND49_NJ
                  J = JOFF + Y
               DO X = 1, ND49_NI
                  I = GET_I( X )
                  Q(X,Y,K) = Q(X,Y,K) +
     $                       ODMDUST(I,J,L,IWVSELECT(1,1),ISPC)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO
            ELSE
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
               DO K = 1, ND49_NL
                  L = LOFF + K
               DO Y = 1, ND49_NJ
                  J = JOFF + Y
               DO X = 1, ND49_NI
                  I = GET_I( X )
                  ! Interpolated using angstrom exponent between
                  ! Closest available wavelengths
                  ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                  !catch any zero values before interpolation
                  IF ((ODMDUST(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND.
     &                 (ODMDUST(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                    Q(X,Y,K) = Q(X,Y,K) +
     &              (ODMDUST(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)**
     &              (BCOEF_WV(1)*LOG(ODMDUST(I,J,L,IWVSELECT(1,1),ISPC)/
     &               ODMDUST(I,J,L,IWVSELECT(2,1),ISPC))))
                  ENDIF
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO
            ENDIF

         ELSE IF ( N == 530 ) THEN

            !--------------------------------------
            ! PAR Diffuse [hPa]
            !--------------------------------------
            CATEGORY = 'DAO-FLDS'
            UNIT     = 'W/m2'
            GMNL     = ND49_NL
            GMTRC    = 20

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = State_Met%PARDF(I,J)
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 531 ) THEN

            !--------------------------------------
            ! PAR Direct [hPa]
            !--------------------------------------
            CATEGORY = 'DAO-FLDS'
            UNIT     = 'W/m2'
            GMNL     = ND49_NL
            GMTRC    = 21

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = State_Met%PARDR(I,J)
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 532 ) THEN

            !--------------------------------------
            ! Daily LAI [m2/m2]
            !--------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'm2/m2'
            GMNL     = 1
            GMTRC    = 9
           
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = State_Met%MODISLAI(I,J) 
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 533 ) THEN
            
            !--------------------------------------
            ! T @ 2m (mpb,2008)
            ! [K]
            !--------------------------------------
            CATEGORY = 'DAO-FLDS'
            UNIT     = 'K'     
            GMNL     = 1
            GMTRC    = 5

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = State_Met%TS( I , J )
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

!##############################################################################
!###  HEMCO NOW ARCHIVES THE BIOGENIC EMISSION TIMESERIES DIAGNOSTICS,      ###
!###  AS WELL AS THE SOIL NOx TIMESERIES DIAGNOSTICS.                       ###
!###  SEE THE NOTE ABOVE FOR INSTRUCTIONS ON HOW TO SAVE THE EQUIVALENT     ###
!###  OUTPUTS TO A netCDF FILE. (bmy, mps, 5/26/15)                         ###
!##############################################################################

        ELSE IF ( N == 534 ) THEN

            !--------------------------------------
            ! Accumulated Precip [mm/day]
            !--------------------------------------
            CATEGORY = 'DAO-FLDS'
            UNIT     = 'mm/day'
            GMNL     = 1
            GMTRC    = 3
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) =  State_Met%PRECTOT(I,J)
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 535 ) THEN
            
            !-----------------------------------
            ! Cosine of SZA  (eam, 2014)
            !-----------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'unitless'     
            GMNL     = 1
            GMTRC    = 33

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = State_Met%SUNCOSmid( I, J )
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 536 .and. LSOA ) THEN
            
            !-------------------------------------
            ! Total organic aerosol [ug/m3] (eam, 2014)
            !-------------------------------------
            CATEGORY = 'IJ-SOA-$'
            UNIT     = 'ug/m3'
            GMNL     = ND49_NL
            GMTRC    = 6

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = SAVEOA(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 537 ) THEN
            
            !-------------------------------------
            ! ISORROPIA aerosol pH (unitless) 
            ! (eam, 08/2015) 
            !-------------------------------------
            CATEGORY = 'ISOROP-$'
            UNIT     = 'unitless'
            GMNL     = ND49_NL
            GMTRC    = 1

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = State_Chm%phSav(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 538 ) THEN
            
            !-------------------------------------
            ! ISORROPIA aerosol water (ug/m3 air) 
            ! (eam, 08/2015) 
            !-------------------------------------
            CATEGORY = 'ISOROP-$'
            UNIT     = 'ug/m3'
            GMNL     = ND49_NL
            GMTRC    = 3

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = State_Chm%WaterSav(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

!==============================================================================
! This diagnostic requires some fixes. Comment out for now (mps,3/2/18)
!         ELSE IF ( N == 539 ) THEN
!            
!            !-----------------------------------------
!            ! Adjust ozone from lowest model level to 10m
!            ! (krt, 11/2017) 
!            !-----------------------------------------
!            CATEGORY = 'O3@10m-$'
!            UNIT     = 'ppbv'
!            GMNL     = 1
!            GMTRC    = 2
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, X, Y, TMPO3, CONV )
!            DO Y = 1, ND49_NJ
!               J = JOFF + Y
!            DO X = 1, ND49_NI
!               I = GET_I( X )
!               TMPO3 = SPC(I,J,1,id_O3) * ( AIRMW 
!     &                 / State_Chm%SpcData(id_O3)%Info%emMW_g ) ! v/v
!               
!               CONV = (1-(AD_RA(I,J)/100.)* AD44b(I,J,id_O3) )
!               IF (CONV .lt. 0) CONV = 1.0
!
!               ! Case where Monin-Obhukov similarity does not apply
!               ! and the correction is not valid
!               ! Drydep code says that similarity applies over the
!               ! range -2.5 <= z/L <= 1.5
!               IF ( AD_ZL(I,J) .lt. -2.5   .or.
!     &              AD_ZL(I,J) .gt.  1.5 ) THEN
!                  Q(X,Y,1) = -999e+0_fp
!               ELSE
!                  Q(X,Y,1)= CONV *TMPO3*1e+9_fp
!               ENDIF
!               
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!==============================================================================

         ELSE

            ! Skip other tracers
            CYCLE

         ENDIF

         IF ( DO_DIAG_WRITE ) THEN

         !==============================================================
         ! Save this data block to the ND49 timeseries file
         !==============================================================
            CALL BPCH2( IU_ND49,      MODELNAME,    LONRES,   
     &        LATRES,       HALFPOLAR,    CENTER180, 
     &           CATEGORY,     GMTRC,        UNIT,      
     &           TAU,          TAU,          RESERVED,  
     &           ND49_NI,      ND49_NJ,      GMNL,  
     &           Input_Opt%ND49_IMIN+I0,
     &           Input_Opt%ND49_JMIN+J0,
     &           Input_Opt%ND49_LMIN, 
     &           REAL( Q(1:ND49_NI, 1:ND49_NJ, 1:GMNL),4 ) )

         ENDIF

      ENDDO
            
      !=================================================================
      ! Close the file at the proper time
      !=================================================================

      IF ( ITS_TIME_TO_CLOSE_FILE() .and. DO_DIAG_WRITE ) THEN
         
         ! Expand date tokens in the file name
         FILENAME = TRIM( Input_Opt%ND49_FILE )
!--- Previous to (ccc, 8/12/09)
!         CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )
         CALL EXPAND_DATE( FILENAME, GET_NYMD_DIAG(), GET_NHMS() )

         ! Echo info
         WRITE( 6, 120 ) TRIM( FILENAME )
 120     FORMAT( '     - DIAG49: Closing file : ', a )

         ! Close file
         CLOSE( IU_ND49 ) 
      ENDIF

      ! Free pointers
      Spc => NULL()
#endif

      END SUBROUTINE DIAG49
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_to_close_file
!
! !DESCRIPTION: Function ITS\_TIME\_TO\_CLOSE\_FILE returns TRUE if it's 
!  time to close the ND49 bpch file before the end of the day.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_TO_CLOSE_FILE() RESULT( ITS_TIME )
!
! !USES:
!
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      USE TIME_MOD, ONLY : GET_HOUR
      USE TIME_MOD, ONLY : GET_MINUTE
      USE TIME_MOD, ONLY : GET_SECOND
#endif
!
! !RETURN VALUE:
!
      LOGICAL :: ITS_TIME
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) The time is already updated to the next time step (ccc, 8/12/09)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  13 Mar 2015 - M. Sulprizio- Bug fix: Allow for temporal resolution finer
!                              than 60 min (E. Marais)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      REAL(fp) :: HR1

      !=================================================================
      ! ITS_TIME_TO_CLOSE_FILE begins here!
      !=================================================================

      ! Current hour
      HR1      = GET_HOUR() + ( GET_MINUTE() / 60e+0_fp ) 
     &           + ( GET_SECOND() / 3600e+0_fp )

      ! Allow for temporal resolution finer than 60 min
      HR1      = HR1 * 1e+2_fp

      ITS_TIME = ( INT( HR1 ) == 00 )
#endif

      END FUNCTION ITS_TIME_TO_CLOSE_FILE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_diag49
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_DIAG49 returns TRUE if ND49 is 
!  turned on and it is time to call DIAG49 -- or FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_DIAG49( Input_Opt ) RESULT( ITS_TIME )
!
! !USES:
!
      USE Input_Opt_Mod, ONLY : OptInput
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      USE ERROR_MOD,     ONLY : GEOS_CHEM_STOP
      USE TIME_MOD,      ONLY : GET_ELAPSED_SEC
      USE TIME_MOD,      ONLY : GET_TS_DIAG
#endif
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options
!
! !RETURN VALUE:
!
      LOGICAL :: ITS_TIME
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Add a check on the output frequency for validity compared to time 
!        steps used. (ccc, 5/21/09)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  06 Feb 2018 - E. Lundgren - Change ND49 frequency from min to sec
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      INTEGER       :: XSEC, TS_DIAG
      LOGICAL, SAVE :: FIRST = .TRUE.

      !=================================================================
      ! ITS_TIME_FOR_DIAG49 begins here!
      !=================================================================

      IF ( Input_Opt%DO_ND49 ) THEN
         IF ( FIRST ) THEN
            TS_DIAG = GET_TS_DIAG()
            
            ! Check if ND49_FREQ is a multiple of TS_DIAG
            IF ( MOD( Input_Opt%ND49_FREQ, TS_DIAG ) /= 0 ) THEN
               WRITE( 6, 100 ) 'ND49', Input_Opt%ND49_FREQ, TS_DIAG
 100           FORMAT( 'The ',a,' output frequency must be a multiple '
     &              'of the largest time step:', i5, i5 )
               CALL GEOS_CHEM_STOP
            ENDIF
            FIRST = .FALSE.
         ENDIF
         
         ! Time already elapsed in this run
         XSEC     = GET_ELAPSED_SEC()
         
         ! Is the elapsed time a multiple of ND49_FREQ?
         ITS_TIME = ( Input_Opt%DO_ND49 .and.
     &                MOD( XSEC, Input_Opt%ND49_FREQ ) == 0 )
      ELSE
         ITS_TIME = Input_Opt%DO_ND49
      ENDIF

#endif            
      END FUNCTION ITS_TIME_FOR_DIAG49
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
      INTEGER, INTENT(IN) :: X   ! Relative longitude index (used by Q array)
!
! !RETURN VALUE:
!
      INTEGER             :: I   ! Absolute longitude index
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
! !IROUTINE: init_diag49
!
! !DESCRIPTION: Subroutine INIT\_DIAG49 allocates and zeroes all module 
!  arrays.  It also gets values for module variables from "input\_mod.f". 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_DIAG49( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      USE BPCH2_MOD,     ONLY : GET_MODELNAME
      USE BPCH2_MOD,     ONLY : GET_HALFPOLAR
      USE CMN_SIZE_MOD        ! Size parameters
      USE ErrCode_Mod
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE GC_GRID_MOD,   ONLY : GET_XOFFSET
      USE GC_GRID_MOD,   ONLY : GET_YOFFSET
      USE GC_GRID_MOD,   ONLY : ITS_A_NESTED_GRID
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
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now get I0 and J0 correctly for nested grid simulations (bmy, 11/9/04)
!  (2 ) Now call GET_HALFPOLAR from "bpch2_mod.f" to get the HALFPOLAR flag 
!        value for GEOS or GCAP grids. (bmy, 6/28/05)
!  (3 ) Now allow ND49_IMIN to be equal to ND49_IMAX and ND49_JMIN to be
!        equal to ND49_JMAX.  This will allow us to save out longitude
!        or latitude transects.  (cdh, bmy, 11/30/06)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      CHARACTER(LEN=255) :: LOCATION
#endif

      !=================================================================
      ! INIT_DIAG49 begins here!
      !=================================================================

      ! Assume success
      RC = GC_SUCCESS

#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )

      ! Initialize
      LOCATION               = 'INIT_DIAG49 ("diag49_mod.f")'
     
      ! Return if we are not saving ND49 diagnostics
      IF ( .not. Input_Opt%DO_ND49 ) RETURN

      !=================================================================
      ! Compute lon, lat, alt extents and check for errors
      !=================================================================

      ! Get grid offsets for error checking
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

      ! Error check ND49_IMIN
      IF ( Input_Opt%ND49_IMIN+I0 < 1       .or.
     &     Input_Opt%ND49_IMIN+I0 > IIPAR ) THEN
         CALL ERROR_STOP( 'Bad ND49_IMIN value!', LOCATION )
      ENDIF

      ! Error check ND49_IMAX
      IF ( Input_Opt%ND49_IMAX+I0 < 1       .or.
     &     Input_Opt%ND49_IMAX+I0 > IIPAR ) THEN
         CALL ERROR_STOP( 'Bad ND49_IMAX value!', LOCATION )
      ENDIF

      ! Compute longitude limits to write to disk 
      ! Also handle wrapping around the date line
      IF ( Input_Opt%ND49_IMAX >= Input_Opt%ND49_IMIN ) THEN
         ND49_NI = ( Input_Opt%ND49_IMAX - Input_Opt%ND49_IMIN ) + 1
      ELSE 
         ND49_NI = ( IIPAR - Input_Opt%ND49_IMIN ) + 1 +
     &             Input_Opt%ND49_IMAX
         WRITE( 6, '(a)' ) 'We are wrapping over the date line!'
      ENDIF

      ! Make sure that ND49_NI <= IIPAR
      IF ( ND49_NI > IIPAR ) THEN
         CALL ERROR_STOP( 'Too many longitudes!', LOCATION )
      ENDIF

      !-----------
      ! Latitude
      !-----------
      
      ! Error check JMIN_AREA
      IF ( Input_Opt%ND49_JMIN+J0 < 1       .or.
     &     Input_Opt%ND49_JMIN+J0 > JJPAR ) THEN
         CALL ERROR_STOP( 'Bad ND49_JMIN value!', LOCATION)
      ENDIF
     
      ! Error check JMAX_AREA
      IF ( Input_Opt%ND49_JMAX+J0 < 1       .or.
     &     Input_Opt%ND49_JMAX+J0 > JJPAR ) THEN
         CALL ERROR_STOP( 'Bad ND49_JMAX value!', LOCATION)
      ENDIF

      ! Compute latitude limits to write to disk (bey, bmy, 3/16/99)
      IF ( Input_Opt%ND49_JMAX >= Input_Opt%ND49_JMIN ) THEN      
         ND49_NJ = ( Input_Opt%ND49_JMAX - Input_Opt%ND49_JMIN ) + 1
      ELSE
         CALL ERROR_STOP( 'ND49_JMAX < ND49_JMIN!', LOCATION )
      ENDIF     
  
      !-----------
      ! Altitude
      !-----------

      ! Error check ND49_LMIN, ND49_LMAX
      IF ( Input_Opt%ND49_LMIN < 1       .or.
     &     Input_Opt%ND49_LMAX > LLPAR ) THEN 
         CALL ERROR_STOP( 'Bad ND49 altitude values!', LOCATION )
      ENDIF

      ! # of levels to save in ND49 timeseries
      IF ( Input_Opt%ND49_LMAX >= Input_Opt%ND49_LMIN ) THEN  
         ND49_NL = ( Input_Opt%ND49_LMAX - Input_Opt%ND49_LMIN ) + 1
      ELSE
         CALL ERROR_STOP( 'ND49_LMAX < ND49_LMIN!', LOCATION )
      ENDIF

      !-----------
      ! Offsets
      !-----------
      IOFF      = Input_Opt%ND49_IMIN - 1
      JOFF      = Input_Opt%ND49_JMIN - 1
      LOFF      = Input_Opt%ND49_LMIN - 1

      !-----------
      ! For bpch
      !-----------
      TITLE     = 'GEOS-CHEM DIAG49 instantaneous timeseries'
      LONRES    = DISIZE
      LATRES    = DJSIZE
      MODELNAME = GET_MODELNAME()
      HALFPOLAR = GET_HALFPOLAR()
      
      ! Reset grid offsets to global values for bpch write
      I0        = GET_XOFFSET( GLOBAL=.TRUE. )
      J0        = GET_YOFFSET( GLOBAL=.TRUE. )      
#endif

      END SUBROUTINE INIT_DIAG49
!EOC
      END MODULE DIAG49_MOD
