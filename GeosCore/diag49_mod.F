!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag49_mod
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
      IMPLICIT NONE
#     include "define.h"
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      LOGICAL, PUBLIC :: DO_SAVE_DIAG49
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
!  ND49 tracer numbers:
!  ============================================================================
!  1 - N_TRACERS : GEOS-CHEM transported tracers            [v/v        ]
!  76            : OH concentration                         [molec/cm3  ]
!  77            : NO2 concentration                        [v/v        ]
!  78            : PBL heights                              [m          ]
!  79            : PBL heights                              [levels     ]
!  80            : Air density                              [molec/cm3  ]
!  81            : 3-D Cloud fractions                      [unitless   ]
!  82            : Column optical depths                    [unitless   ]
!  83            : Cloud top heights                        [hPa        ]
!  84            : Sulfate aerosol optical depth            [unitless   ]
!  85            : Black carbon aerosol optical depth       [unitless   ]
!  86            : Organic carbon aerosol optical depth     [unitless   ]
!  87            : Accumulation mode seasalt optical depth  [unitless   ]
!  88            : Coarse mode seasalt optical depth        [unitless   ]
!  89            : Total dust optical depth                 [unitless   ]
!  90            : Total seasalt tracer concentration       [unitless   ]
!  91            : Pure O3 (not Ox) concentration           [v/v        ]
!  92            : NO concentration                         [v/v        ]
!  93            : NOy concentration                        [v/v        ]
!  94            : Grid box heights                         [m          ]
!  95            : Relative Humidity                        [%          ]
!  96            : Sea level pressure                       [hPa        ]
!  97            : Zonal wind (a.k.a. U-wind)               [m/s        ]
!  98            : Meridional wind (a.k.a. V-wind)          [m/s        ]
!  99            : PEDGE-$ (Pressure @ level edges          [hPa        ]
!  100           : Temperature                              [K          ]
!  101           : PAR direct                               [hPa        ]
!  102           : PAR diffuse                              [hPa        ]
!  103           : Daily LAI                                [hPa        ]
!  104           : Temperature at 2m                        [K          ]
!  105           : Isoprene emissions                       [atomC/cm2/s]
!  106           : Total Monoterpene emissions              [atomC/cm2/s]
!  107           : Methyl Butanol emissions                 [atomC/cm2/s]
!  108           : Alpha-Pinene emissions                   [atomC/cm2/s]
!  109           : Beta-Pinene emissions                    [atomC/cm2/s]
!  110           : Limonene emissions                       [atomC/cm2/s]
!  111           : Sabinene emissions                       [atomC/cm2/s]
!  112           : Myrcene emissions                        [atomC/cm2/s]
!  113           : 3-Carene emissions                       [atomC/cm2/s]
!  114           : Ocimene emissions                        [atomC/cm2/s]
!  115-121       : size resolved dust optical depth         [unitless   ]
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES
      !
      ! I0               : Offset between global & nested grid
      ! J0               : Offset between global & nested grid
      ! IOFF             : Longitude offset
      ! JOFF             : Latitude offset
      ! LOFF             : Altitude offset
      ! ND49_IMIN        : Minimum latitude  index for DIAG51 region
      ! ND49_IMAX        : Maximum latitude  index for DIAG51 region
      ! ND49_JMIN        : Minimum longitude index for DIAG51 region
      ! ND49_JMAX        : Maximum longitude index for DIAG51 region
      ! ND49_LMIN        : Minimum altitude  index for DIAG51 region
      ! ND49_LMAX        : Minimum latitude  index for DIAG51 region
      ! ND49_NI          : Number of longitudes in DIAG51 region 
      ! ND49_NJ          : Number of latitudes  in DIAG51 region
      ! ND49_NL          : Number of levels     in DIAG51 region
      ! ND49_N_TRACERS   : Number of tracers for DIAG51
      ! ND49_OUTPUT_FILE : Name of bpch file w  timeseries data
      ! ND49_TRACERS     : Array of DIAG51 tracer numbers
      ! HALFPOLAR        : Used for bpch file output
      ! CENTER180        : Used for bpch file output
      ! LONRES           : Used for bpch file output
      ! LATRES           : Used for bpch file output
      ! MODELNAME        : Used for bpch file output
      ! RESERVED         : Used for bpch file output
      !=================================================================

      INTEGER            :: IOFF,           JOFF,   LOFF
      INTEGER            :: I0,             J0
      INTEGER            :: ND49_N_TRACERS, ND49_TRACERS(120)
      INTEGER            :: ND49_IMIN,      ND49_IMAX
      INTEGER            :: ND49_JMIN,      ND49_JMAX
      INTEGER            :: ND49_LMIN,      ND49_LMAX
      INTEGER            :: ND49_FREQ,      ND49_NI
      INTEGER            :: ND49_NJ,        ND49_NL
      INTEGER            :: HALFPOLAR
      INTEGER, PARAMETER :: CENTER180=1 
      REAL*4             :: LONRES,         LATRES
      CHARACTER(LEN=20)  :: MODELNAME
      CHARACTER(LEN=40)  :: RESERVED = ''
      CHARACTER(LEN=80)  :: TITLE
      CHARACTER(LEN=255) :: ND49_OUTPUT_FILE

      ! Logical file unit for ND49 diagnostic
      INTEGER            :: IU_ND49

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      SUBROUTINE DIAG49
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : BPCH2,   OPEN_BPCH2_FOR_WRITE
      USE DAO_MOD,       ONLY : AD,      AIRDEN, BXHEIGHT, CLDF 
      USE DAO_MOD,       ONLY : CLDTOPS, OPTD,   RH,       SLP     
      USE DAO_MOD,       ONLY : T,       UWND,   VWND
      USE DAO_MOD,       ONLY : TS
      USE DAO_MOD,       ONLY : PARDF, PARDR
      USE MODIS_LAI_MOD, ONLY : ISOLAI => GC_LAI
      USE GRID_MOD,      ONLY : GET_XOFFSET,        GET_YOFFSET
      USE inquireMod,    ONLY : findFreeLun
      USE TIME_MOD,      ONLY : EXPAND_DATE
      USE TIME_MOD,      ONLY : GET_NYMD,           GET_NHMS
      USE TIME_MOD,      ONLY : GET_NYMD_DIAG,      GET_TS_DIAG
      USE TIME_MOD,      ONLY : GET_TAU,            GET_HOUR
      USE TIME_MOD,      ONLY : ITS_A_NEW_DAY,      TIMESTAMP_STRING
      USE PBL_MIX_MOD,   ONLY : GET_PBL_TOP_L,      GET_PBL_TOP_m
      USE TRACER_MOD,    ONLY : ITS_A_FULLCHEM_SIM, N_TRACERS
      USE TRACER_MOD,    ONLY : STT,                TCVV
      USE TRACER_MOD,    ONLY : XNUMOLAIR
      USE PRESSURE_MOD,  ONLY : GET_PEDGE
      USE TRACERID_MOD,  ONLY : IDTHNO3, IDTHNO4, IDTN2O5, IDTNOX  
      USE TRACERID_MOD,  ONLY : IDTPAN,  IDTPMN,  IDTPPN,  IDTOX   
      USE TRACERID_MOD,  ONLY : IDTR4N2, IDTSALA, IDTSALC 
      USE LOGICAL_MOD,   ONLY : DO_DIAG_WRITE

      USE CMN_FJ_MOD,    ONLY : JPMAX,   JPPJ     ! FAST-J stuff
      USE JV_CMN_MOD                              ! ODAER, QAA, QAA_AOD (clh)
      USE CMN_O3_MOD	                          ! Pure O3, SAVENO2
      USE CMN_GCTM_MOD                            ! XTRA2

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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE            :: FIRST  = .TRUE.
      LOGICAL, SAVE            :: IS_FULLCHEM, IS_NOx,     IS_Ox 
      LOGICAL, SAVE            :: IS_NOy,      IS_CLDTOPS, IS_OPTD
      LOGICAL, SAVE            :: IS_SEASALT,  IS_SLP
      LOGICAL                  :: LINTERP
      INTEGER                  :: IOS,  GMTRC, GMNL, I, J, K, L 
      INTEGER                  :: N,    R,     H,    W, X, Y
      INTEGER                  :: NHMS, TS_DIAG
      INTEGER                  :: ISPC
      REAL*8                   :: TAU, TMP,   SCALEAODnm
      REAL*8                   :: Q( ND49_NI, ND49_NJ, ND49_NL )
      CHARACTER(LEN=16)        :: STAMP
      CHARACTER(LEN=40)        :: CATEGORY
      CHARACTER(LEN=40)        :: UNIT
      CHARACTER(LEN=255)       :: FILENAME

      ! Aerosol types (rvm, aad, bmy, 7/20/04)
      INTEGER                  :: IND(6) = (/ 22, 29, 36, 43, 50, 15 /)
      !=================================================================
      ! DIAG49 begins here!
      !=================================================================

      ! Set logical flags on first timestep
      IF ( FIRST ) THEN
         IS_CLDTOPS  = ALLOCATED( CLDTOPS )
         IS_OPTD     = ALLOCATED( OPTD    )
         IS_SLP      = ALLOCATED( SLP     )
         IS_FULLCHEM = ITS_A_FULLCHEM_SIM()
         IS_SEASALT  = ( IDTSALA > 0 .and. IDTSALC > 0 )
         IS_Ox       = ( IS_FULLCHEM .and. IDTOX   > 0 )
         IS_NOx      = ( IS_FULLCHEM .and. IDTNOX  > 0 )
         IS_NOy      = ( IS_FULLCHEM .and. 
     &                   IDTNOX  > 0 .and. IDTPAN  > 0 .and.
     &                   IDTHNO3 > 0 .and. IDTPMN  > 0 .and.
     &                   IDTPPN  > 0 .and. IDTR4N2 > 0 .and.
     &                   IDTN2O5 > 0 .and. IDTHNO4 > 0 ) 
         FIRST       = .FALSE.
      ENDIF

      !=================================================================
      ! If it's a new day, open a new BPCH file and write file header
      ! We need to check if it's a new day + 1 ND49 time step (ccc, 8/12/09)
      !=================================================================
!--- Previous to (ccc, 8/12/09)
!      IF ( ITS_A_NEW_DAY() ) THEN
      NHMS    = GET_NHMS()
      TS_DIAG = ND49_FREQ

      ! To change TS_DIAG to NHMS format
      TS_DIAG = TS_DIAG/60 * 10000 + (TS_DIAG - (TS_DIAG/60)*60) * 100  

      IF ( NHMS == TS_DIAG ) THEN     ! It's a new day for diagnostics.

         ! Find a free file LUN
         IU_ND49 = findFreeLUN()

         ! Expand date tokens in the file name
         FILENAME = TRIM( ND49_OUTPUT_FILE )
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
      DO W = 1, ND49_N_TRACERS

         ! ND49 tracer number
         N = ND49_TRACERS(W)

         ! Zero summing array
         Q = 0d0

         ! Test by tracer number
         IF ( N <= N_TRACERS ) THEN

            !-------------------------------------
            ! GEOS-CHEM tracers [v/v]
            !-------------------------------------
            CATEGORY = 'IJ-AVG-$'
            UNIT     = ''           ! Let GAMAP pick the unit
            GMNL     = ND49_NL
            GMTRC    = N
      
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = STT(I,J,L,N) * TCVV(N) / AD(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 91 .and. IS_Ox ) THEN

            !-------------------------------------
            ! PURE O3 CONCENTRATION [v/v]
            !-------------------------------------
            CATEGORY = 'IJ-AVG-$'
            UNIT     = ''           ! Let GAMAP pick the unit
            GMNL     = ND49_NL
            GMTRC    = N_TRACERS + 1

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = STT(I,J,L,IDTOX) * TCVV(IDTOX)  / 
     &                    AD(I,J,L)        * FRACO3(I,J,L) 
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 92 .and. IS_NOx ) THEN
            
            !-------------------------------------
            ! NO CONCENTRATION [v/v]
            !-------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''           ! Let GAMAP pick the unit
            GMNL     = ND49_NL
            GMTRC    = 9
               
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )   
               Q(X,Y,K) = STT(I,J,L,IDTNOX) * TCVV(IDTNOX) * 
     &                    FRACNO(I,J,L)     / AD(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 93 .and. IS_NOy ) THEN

            !--------------------------------------
            ! NOy CONCENTRATION [v/v]
            !--------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND49_NL
            GMTRC    = 3

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
               TMP = 0d0
            
               ! NOx
               TMP = TMP + ( TCVV(IDTNOX)        *
     &                       STT(I,J,L,IDTNOX)   / AD(I,J,L) )
               ! PAN
               TMP = TMP + ( TCVV(IDTPAN)        * 
     &                       STT(I,J,L,IDTPAN)   / AD(I,J,L) )

               ! HNO3
               TMP = TMP + ( TCVV(IDTHNO3)       *
     &                       STT(I,J,L,IDTHNO3)  / AD(I,J,L) )
            
               ! PMN
               TMP = TMP + ( TCVV(IDTPMN)        *
     &                       STT(I,J,L,IDTPMN)   / AD(I,J,L) )

               ! PPN
               TMP = TMP + ( TCVV(IDTPPN)        * 
     &                       STT(I,J,L,IDTPPN)   / AD(I,J,L) )
 
               ! R4N2
               TMP = TMP + ( TCVV(IDTR4N2)       *
     &                       STT(I,J,L,IDTR4N2)  / AD(I,J,L) )
            
               ! N2O5
               TMP = TMP + ( 2d0 * TCVV(IDTN2O5) *
     &                       STT(I,J,L,IDTN2O5)  / AD(I,J,L) )
                        
               ! HNO4
               TMP = TMP + ( TCVV(IDTHNO4)       *
     &                       STT(I,J,L,IDTHNO4)  / AD(I,J,L) )

               ! Save afternoon points
               Q(X,Y,K) = Q(X,Y,K) + TMP
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 76 .and. IS_FULLCHEM ) THEN

            !-------------------------------------
            ! OH CONCENTRATION [molec/cm3]
            !-------------------------------------              
            CATEGORY = 'TIME-SER'
            UNIT     = 'molec/cm3'
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
               Q(X,Y,K) = SAVEOH(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 77 .and. IS_FULLCHEM ) THEN

            !-------------------------------------
            ! NO2 CONCENTRATION [molec/cm3]
            !-------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''           ! Let GAMAP pick the unit
            GMNL     = ND49_NL
            GMTRC    = 25

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = SAVENO2(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 78 ) THEN

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

         ELSE IF ( N == 79 ) THEN

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

         ELSE IF ( N == 80 ) THEN

            !--------------------------------------
            ! AIR DENSITY [molec/cm3]
            !--------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'molec/cm3'
            GMNL     = ND49_NL
            GMTRC    = 22

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = AIRDEN(L,I,J) * XNUMOLAIR * 1d-6
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 81 ) THEN 

            !--------------------------------------
            ! 3-D CLOUD FRACTIONS [unitless]
            !--------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'unitless'
            GMNL     = ND49_NL
            GMTRC    = 19

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = CLDF(K,I,J)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 82 .and. IS_OPTD ) THEN 

            !--------------------------------------
            ! COLUMN OPTICAL DEPTHS [unitless]
            !--------------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'unitless'
            GMNL     = 1
            GMTRC    = 20

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = SUM( OPTD(:,I,J) )
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 83 .and. IS_CLDTOPS ) THEN 

            !--------------------------------------
            ! CLOUD TOP HEIGHTS [hPa]
            !--------------------------------------
            CATEGORY = 'TIME_SER'
            UNIT     = 'hPa'
            GMNL     = ND49_NL
            GMTRC    = 21

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = GET_PEDGE( I, J, CLDTOPS(I,J) )
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 84 ) THEN

            !--------------------------------------
            ! SULFATE AOD @ jv_spec_aod.dat wavelength [unitless]
            !--------------------------------------
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND49_NL
            GMTRC    = 6

            ISPC=1
            IF (LINTERP.EQ..FALSE.) THEN
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

         ELSE IF ( N == 85 ) THEN

            !--------------------------------------
            ! BLACK CARBON AOD @ jv_spec_aod.dat wavelength [unitless]
            !--------------------------------------
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND49_NL
            GMTRC    = 9

            ISPC=2
            IF (LINTERP.EQ..FALSE.) THEN
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
         ELSE IF ( N == 86 ) THEN

            !--------------------------------------
            ! ORGANIC CARBON AOD @ jv_spec_aod.dat wavelength [unitless]
            !--------------------------------------            
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND49_NL
            GMTRC    = 12

            ISPC=3
            IF (LINTERP.EQ..FALSE.) THEN
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
         ELSE IF ( N == 87 ) THEN
            
            !--------------------------------------
            ! ACCUM SEASALT AOD @ jv_spec_aod.dat wavelength [unitless]
            !--------------------------------------
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND49_NL
            GMTRC    = 15

            ISPC=4
            IF (LINTERP.EQ..FALSE.) THEN
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
         ELSE IF ( N == 88 ) THEN

            !--------------------------------------
            ! COARSE SEASALT AOD @ jv_spec_aod.dat wavelength [unitless]
            !--------------------------------------            
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND49_NL
            GMTRC    = 18

            ISPC=5
            IF (LINTERP.EQ..FALSE.) THEN
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
         ELSE IF ( N == 89 ) THEN

            !-----------------------------------
            ! TOTAL DUST OPT DEPTH [unitless]
            !-----------------------------------
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMNL      = ND49_NL
            GMTRC     = 4

            DO ISPC = 1, NDUST

            IF (LINTERP.EQ..FALSE.) THEN
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
            ENDDO
         ELSE IF (N > 114) THEN

            !--------------------------------------
            ! DUST OPTD BINS1-7 @ jv_spec_aod.dat wavelength[unitless]
            ! NOTE: Only archive at chem timestep
            !--------------------------------------
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMNL      = ND49_NL
            GMTRC     = 21+(N-115)

            ISPC = N - 114
            IF (LINTERP.EQ..FALSE.) THEN
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

         ELSE IF ( N == 90 .and. IS_SEASALT ) THEN
            
            !-----------------------------------
            ! TOTAL SEASALT TRACER [v/v]
            !-----------------------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''        ! Let GAMAP pick unit
            GMNL     = ND49_NL
            GMTRC    = 24

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = ( STT(I,J,L,IDTSALA) + STT(I,J,L,IDTSALC) ) *
     &                      TCVV(IDTSALA)      / AD(I,J,L) 
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 94 ) THEN

            !-----------------------------------
            ! GRID BOX HEIGHT [m]
            !----------------------------------- 
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
               Q(X,Y,K) = BXHEIGHT(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 95 ) THEN

            !-----------------------------------
            ! RELATIVE HUMIDITY [%]
            !----------------------------------- 
            CATEGORY = 'TIME-SER'
            UNIT     = '%'
            GMNL     = ND49_NL
            GMTRC    = 17

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, X, Y, K )
            DO K = 1, ND49_NL
               L = LOFF + K
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,K) = RH(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 96 .and. IS_SLP ) THEN

            !-----------------------------------
            ! SEA LEVEL PRESSURE [hPa]
            !----------------------------------- 
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
               Q(X,Y,1) = SLP(I,J)
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
               
         ELSE IF ( N == 97 ) THEN

            !-----------------------------------
            ! ZONAL (U) WIND [m/s]
            !----------------------------------- 
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
               Q(X,Y,K) = UWND(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 98 ) THEN 

            !-----------------------------------
            ! MERIDIONAL (V) WIND [m/s]
            !----------------------------------- 
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
               Q(X,Y,K) = VWND(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 99 ) THEN 

            !-----------------------------------
            ! PEDGE-$ (Prs @ level edges) [hPa]
            !----------------------------------- 
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
               Q(X,Y,K) = GET_PEDGE(I,J,K)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 100 ) THEN

            !-----------------------------------
            ! TEMPERATURE [K]
            !----------------------------------- 
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
               Q(X,Y,K) = T(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 101 ) THEN

            !-----------------------------------
            ! PAR Diffuse [hPa]
            !----------------------------------- 
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
               Q(X,Y,1) = PARDF(I,J)
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 102 ) THEN

            !-----------------------------------
            ! PAR Direct [hPa]
            !----------------------------------- 
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
               Q(X,Y,1) = PARDR(I,J)
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 103 ) THEN

            !-----------------------------------
            ! Daily LAI [hPa]
            !----------------------------------- 
            CATEGORY = 'TIME-SER'
            UNIT     = 'm2/m2'
            GMNL     = ND49_NL
            GMTRC    = 32
           
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = ISOLAI(I,J) 
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 104 ) THEN
            
            !-----------------------------------
            ! T @ 2m (mpb,2008)
            ! [K]
            !-----------------------------------
            CATEGORY = 'DAO-FLDS'
            UNIT     = 'K'     
            GMNL     = ND49_NL
            GMTRC    = 5

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = TS( I , J )
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 105 ) THEN
            
            !-----------------------------------
            ! ISOPRENE EMISSIONS (mpb,2008)
            ! [atom C/cm2/s]
            !-----------------------------------
            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'    
            GMNL     = ND49_NL
            GMTRC    = 1

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = EMISS_BVOC( I , J , 1 )
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 106 ) THEN
            
            !--------------------------------------
            ! TOTAL MONTERPENE EMISSIONS (mpb,2008)
            ! [atom C/cm2/s]
            !--------------------------------------
            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'       
            GMNL     = ND49_NL
            GMTRC    = 4

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = EMISS_BVOC( I , J , 2 )
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 107 ) THEN
            
            !-----------------------------------
            ! MBO EMISSIONS  (mpb,2008)
            ! [atom C/cm2/s]
            !-----------------------------------
            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'     
            GMNL     = ND49_NL
            GMTRC    = 5

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = EMISS_BVOC( I , J , 3 )
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 108 ) THEN
            
            !-----------------------------------
            ! Alpha-Pinene EMISSIONS  (mpb,2008)
            ! [atom C/cm2/s]
            !-----------------------------------
            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'     
            GMNL     = ND49_NL
            GMTRC    = 7

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = EMISS_BVOC( I , J , 4 )
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 109 ) THEN
            
            !-----------------------------------
            ! Beta-Pinene EMISSIONS  (mpb,2008)
            ! [atom C/cm2/s]
            !-----------------------------------
            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'     
            GMNL     = ND49_NL
            GMTRC    = 8

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = EMISS_BVOC( I , J , 5 )
            ENDDO
            ENDDO
!$OMP END PARALLEL DO


         ELSE IF ( N == 110 ) THEN
            
            !-----------------------------------
            ! Limonene EMISSIONS  (mpb,2008)
            ! [atom C/cm2/s]
            !-----------------------------------
            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'     
            GMNL     = ND49_NL
            GMTRC    = 9

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = EMISS_BVOC( I , J , 6 )
            ENDDO
            ENDDO
!$OMP END PARALLEL DO


         ELSE IF ( N == 111 ) THEN
            
            !-----------------------------------
            ! Sabinene EMISSIONS  (mpb,2008)
            ! [atom C/cm2/s]
            !-----------------------------------
            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'     
            GMNL     = ND49_NL
            GMTRC    = 10

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = EMISS_BVOC( I , J , 7 )
            ENDDO
            ENDDO
!$OMP END PARALLEL DO


         ELSE IF ( N == 112 ) THEN
            
            !-----------------------------------
            ! Myrcene EMISSIONS  (mpb,2008)
            ! [atom C/cm2/s]
            !-----------------------------------
            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'     
            GMNL     = ND49_NL
            GMTRC    = 11

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = EMISS_BVOC( I , J , 8 )
            ENDDO
            ENDDO
!$OMP END PARALLEL DO


         ELSE IF ( N == 113 ) THEN
            
            !-----------------------------------
            ! 3-Carene EMISSIONS  (mpb,2008)
            ! [atom C/cm2/s]
            !-----------------------------------
            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'     
            GMNL     = ND49_NL
            GMTRC    = 12

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = EMISS_BVOC( I , J , 9 )
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE IF ( N == 114 ) THEN
            
            !-----------------------------------
            ! Ocimene EMISSIONS  (mpb,2008)
            ! [atom C/cm2/s]
            !-----------------------------------
            CATEGORY = 'BIOGSRCE'
            UNIT     = 'atomC/cm2/s'     
            GMNL     = ND49_NL
            GMTRC    = 13

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y )
            DO Y = 1, ND49_NJ
               J = JOFF + Y
            DO X = 1, ND49_NI
               I = GET_I( X )
               Q(X,Y,1) = EMISS_BVOC( I , J , 10 )
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE

            ! Skip
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
     &           ND49_IMIN+I0, ND49_JMIN+J0, ND49_LMIN, 
     &           REAL( Q(1:ND49_NI, 1:ND49_NJ, 1:GMNL),4 ) )

         ENDIF

      ENDDO
            
      !=================================================================
      ! Close the file at the proper time
      !=================================================================
      IF ( ITS_TIME_TO_CLOSE_FILE() .and. DO_DIAG_WRITE ) THEN

         ! Expand date tokens in the file name
         FILENAME = TRIM( ND49_OUTPUT_FILE )
!--- Previous to (ccc, 8/12/09)
!         CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )
         CALL EXPAND_DATE( FILENAME, GET_NYMD_DIAG(), GET_NHMS() )

         ! Echo info
         WRITE( 6, 120 ) TRIM( FILENAME )
 120     FORMAT( '     - DIAG49: Closing file : ', a )

         ! Close file
         CLOSE( IU_ND49 ) 
      ENDIF

      END SUBROUTINE DIAG49
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      USE TIME_MOD, ONLY : GET_HOUR
      USE TIME_MOD, ONLY : GET_MINUTE
!
! !RETURN VALUE:
!
      LOGICAL :: ITS_TIME
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) The time is already updated to the next time step (ccc, 8/12/09)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8 :: HR1

      !=================================================================
      ! ITS_TIME_TO_CLOSE_FILE begins here!
      !=================================================================

      ! Current hour
      HR1      = GET_HOUR() + ( GET_MINUTE() / 60d0 )

!--- Previous to (ccc, 8/12/09)
!      ! Hour at the next dynamic timestep
!      HR2      = HR1        + ( ND49_FREQ / 60d0 )

      ! If the next dyn step is the start of a new day, return TRUE
!--- Previous to (ccc, 11/11/10)
!       HR1 varies between 00 and 23:59. So compares to 00 not 24 anymore.
!      ITS_TIME = ( INT( HR1 ) == 24 )

      ITS_TIME = ( INT( HR1 ) == 00 )

      END FUNCTION ITS_TIME_TO_CLOSE_FILE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      FUNCTION ITS_TIME_FOR_DIAG49() RESULT( ITS_TIME )
!
! !USES:
!
      USE TIME_MOD,  ONLY : GET_ELAPSED_MIN
      USE TIME_MOD,  ONLY : GET_TS_DIAG
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER       :: XMIN, TS_DIAG
      LOGICAL, SAVE :: FIRST = .TRUE.

      !=================================================================
      ! ITS_TIME_FOR_DIAG49 begins here!
      !=================================================================

      IF ( DO_SAVE_DIAG49 ) THEN
         IF ( FIRST ) THEN
            TS_DIAG = GET_TS_DIAG()
            
            ! Check if ND49_FREQ is a multiple of TS_DIAG
            IF ( MOD( ND49_FREQ, TS_DIAG ) /= 0 ) THEN
               WRITE( 6, 100 ) 'ND49', ND49_FREQ, TS_DIAG
 100           FORMAT( 'The ',a,' output frequency must be a multiple '
     &              'of the largest time step:', i5, i5 )
               CALL GEOS_CHEM_STOP
            ENDIF
            FIRST = .FALSE.
         ENDIF
         
         ! Time already elapsed in this run
         XMIN     = GET_ELAPSED_MIN()
         
         ! Is the elapsed time a multiple of ND49_FREQ?
         ITS_TIME = ( DO_SAVE_DIAG49 .and. MOD( XMIN, ND49_FREQ ) == 0 )
      ELSE
         ITS_TIME = DO_SAVE_DIAG49
      ENDIF
            
      END FUNCTION ITS_TIME_FOR_DIAG49
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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

      ! Add the offset to X to get I  
      I = IOFF + X

      ! Handle wrapping around the date line, if necessary
      IF ( I > IIPAR ) I = I - IIPAR

      END FUNCTION GET_I
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      SUBROUTINE INIT_DIAG49( DO_ND49, N_ND49, TRACERS, IMIN,    
     &                        IMAX,    JMIN,   JMAX,    LMIN,    
     &                        LMAX,    FREQ,   FILE )
!
! !USES:
!
      USE BPCH2_MOD, ONLY : GET_MODELNAME
      USE BPCH2_MOD, ONLY : GET_HALFPOLAR
      USE GRID_MOD,  ONLY : GET_XOFFSET
      USE GRID_MOD,  ONLY : GET_YOFFSET
      USE GRID_MOD,  ONLY : ITS_A_NESTED_GRID
      USE ERROR_MOD, ONLY : ERROR_STOP

      USE CMN_SIZE_MOD ! Size parameters
!
! !INPUT PARAMETERS: 
!
      ! DO_ND49 : Switch to turn on ND49 timeseries diagnostic
      ! N_ND50  : Number of ND49 read by "input_mod.f"
      ! TRACERS : Array w/ ND49 tracer #'s read by "input_mod.f"
      ! IMIN    : Min longitude index read by "input_mod.f"
      ! IMAX    : Max longitude index read by "input_mod.f" 
      ! JMIN    : Min latitude index read by "input_mod.f" 
      ! JMAX    : Min latitude index read by "input_mod.f" 
      ! LMIN    : Min level index read by "input_mod.f" 
      ! LMAX    : Min level index read by "input_mod.f" 
      ! FREQ    : Frequency for saving to disk [min]
      ! FILE    : ND49 output file name read by "input_mod.f"
      LOGICAL,            INTENT(IN) :: DO_ND49
      INTEGER,            INTENT(IN) :: N_ND49, TRACERS(100)
      INTEGER,            INTENT(IN) :: IMIN,   IMAX 
      INTEGER,            INTENT(IN) :: JMIN,   JMAX      
      INTEGER,            INTENT(IN) :: LMIN,   LMAX 
      INTEGER,            INTENT(IN) :: FREQ
      CHARACTER(LEN=255), INTENT(IN) :: FILE
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255) :: LOCATION
      
      !=================================================================
      ! INIT_DIAG49 begins here!
      !=================================================================

      ! Initialize
      LOCATION               = 'INIT_DIAG49 ("diag49_mod.f")'
      ND49_TRACERS(:)        = 0

      ! Get values from "input_mod.f"
      DO_SAVE_DIAG49         = DO_ND49 
      ND49_N_TRACERS         = N_ND49
      ND49_TRACERS(1:N_ND49) = TRACERS(1:N_ND49)
      ND49_IMIN              = IMIN
      ND49_IMAX              = IMAX
      ND49_JMIN              = JMIN
      ND49_JMAX              = JMAX
      ND49_LMIN              = LMIN
      ND49_LMAX              = LMAX
      ND49_FREQ              = FREQ
      ND49_OUTPUT_FILE       = FILE
     
      ! Return if we are not saving ND49 diagnostics
      IF ( .not. DO_SAVE_DIAG49 ) RETURN

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
      IF ( ND49_IMIN+I0 < 1 .or. ND49_IMIN+I0 > IIPAR ) THEN
         CALL ERROR_STOP( 'Bad ND49_IMIN value!', LOCATION )
      ENDIF

      ! Error check ND49_IMAX
      IF ( ND49_IMAX+I0 < 1 .or. ND49_IMAX+I0 > IIPAR ) THEN
         CALL ERROR_STOP( 'Bad ND49_IMAX value!', LOCATION )
      ENDIF

      ! Compute longitude limits to write to disk 
      ! Also handle wrapping around the date line
      IF ( ND49_IMAX >= ND49_IMIN ) THEN
         ND49_NI = ( ND49_IMAX - ND49_IMIN ) + 1
      ELSE 
         ND49_NI = ( IIPAR - ND49_IMIN ) + 1 + ND49_IMAX
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
      IF ( ND49_JMIN+J0 < 1 .or. ND49_JMIN+J0 > JJPAR ) THEN
         CALL ERROR_STOP( 'Bad ND49_JMIN value!', LOCATION)
      ENDIF
     
      ! Error check JMAX_AREA
      IF ( ND49_JMAX+J0 < 1 .or.ND49_JMAX+J0 > JJPAR ) THEN
         CALL ERROR_STOP( 'Bad ND49_JMAX value!', LOCATION)
      ENDIF

      ! Compute latitude limits to write to disk (bey, bmy, 3/16/99)
      IF ( ND49_JMAX >= ND49_JMIN ) THEN      
         ND49_NJ = ( ND49_JMAX - ND49_JMIN ) + 1
      ELSE
         CALL ERROR_STOP( 'ND49_JMAX < ND49_JMIN!', LOCATION )
      ENDIF     
  
      !-----------
      ! Altitude
      !-----------

      ! Error check ND49_LMIN, ND49_LMAX
      IF ( ND49_LMIN < 1 .or. ND49_LMAX > LLPAR ) THEN 
         CALL ERROR_STOP( 'Bad ND49 altitude values!', LOCATION )
      ENDIF

      ! # of levels to save in ND49 timeseries
      IF ( ND49_LMAX >= ND49_LMIN ) THEN  
         ND49_NL = ( ND49_LMAX - ND49_LMIN ) + 1
      ELSE
         CALL ERROR_STOP( 'ND49_LMAX < ND49_LMIN!', LOCATION )
      ENDIF

      !-----------
      ! Offsets
      !-----------
      IOFF      = ND49_IMIN - 1
      JOFF      = ND49_JMIN - 1
      LOFF      = ND49_LMIN - 1

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

      END SUBROUTINE INIT_DIAG49
!EOC
      END MODULE DIAG49_MOD
