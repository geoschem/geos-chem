! $Id: diag51_mod.f,v 1.19 2005/10/20 14:03:21 bmy Exp $
      MODULE DIAG51_MOD
!
!******************************************************************************
!  Module DIAG51_MOD contains variables and routines to generate save 
!  timeseries data where the local time is between two user-defined limits. 
!  This facilitates comparisons with morning or afternoon-passing satellites
!  such as GOME. (amf, bey, bdf, pip, bmy, 11/30/00, 10/3/05)
!
!  Module Variables:
!  ============================================================================
!  (1 ) DO_SAVE_DIAG51   (LOGICAL ) : Flag to turn on DIAG51 timseries
!  (2 ) GOOD             (INTEGER ) : Array denoting grid boxes w/in LT limits
!  (3 ) GOOD_CT          (INTEGER ) : # of "good" times per grid box
!  (4 ) GOOD_CT_CHEM     (INTEGER ) : # of "good" chemistry timesteps
!  (5 ) ND51_HR_WRITE    (INTEGER ) : Hour at which to save to disk
!  (6 ) I0               (INTEGER ) : Offset between global & nested grid
!  (7 ) J0               (INTEGER ) : Offset between global & nested grid
!  (8 ) IOFF             (INTEGER ) : Longitude offset
!  (9 ) JOFF             (INTEGER ) : Latitude offset
!  (10) LOFF             (INTEGER ) : Altitude offset
!  (11) ND51_HR1         (REAL*8  ) : Starting hour of user-defined LT interval
!  (12) ND51_HR2         (REAL*8  ) : Ending hour of user-defined LT interval
!  (13) ND51_IMIN        (INTEGER ) : Minimum latitude  index for DIAG51 region
!  (14) ND51_IMAX        (INTEGER ) : Maximum latitude  index for DIAG51 region
!  (15) ND51_JMIN        (INTEGER ) : Minimum longitude index for DIAG51 region
!  (16) ND51_JMAX        (INTEGER ) : Maximum longitude index for DIAG51 region
!  (17) ND51_LMIN        (INTEGER ) : Minimum altitude  index for DIAG51 region
!  (18) ND51_LMAX        (INTEGER ) : Minimum latitude  index for DIAG51 region
!  (19) ND51_NI          (INTEGER ) : Number of longitudes in DIAG51 region 
!  (20) ND51_NJ          (INTEGER ) : Number of latitudes  in DIAG51 region
!  (21) ND51_NL          (INTEGER ) : Number of levels     in DIAG51 region
!  (22) ND51_N_TRACERS   (INTEGER ) : Number of tracers for DIAG51
!  (23) ND51_OUTPUT_FILE (CHAR*255) : Name of bpch file w  timeseries data
!  (24) ND51_TRACERS     (INTEGER ) : Array of DIAG51 tracer numbers
!  (25) Q                (REAL*8  ) : Accumulator array for various quantities
!  (26) TAU0             (REAL*8  ) : Starting TAU used to index the bpch file
!  (27) TAU1             (REAL*8  ) : Ending TAU used to index the bpch file
!  (28) HALFPOLAR        (INTEGER ) : Used for bpch file output
!  (29) CENTER180        (INTEGER ) : Used for bpch file output
!  (30) LONRES           (REAL*4  ) : Used for bpch file output
!  (31) LATRES           (REAL*4  ) : Used for bpch file output
!  (32) MODELNAME        (CHAR*20 ) : Used for bpch file output
!  (33) RESERVED         (CHAR*40 ) : Used for bpch file output
!
!  Module Procedures:
!  ============================================================================
!  (1 ) DIAG51                      : Driver subroutine for US grid timeseries 
!  (2 ) GET_LOCAL_TIME              : Computes the local times at each grid box
!  (3 ) WRITE_DIAG51                : Writes timeseries data to a bpch file
!  (4 ) ITS_TIME_FOR_WRITE_DIAG51   : Returns T if it's time to save to disk
!  (5 ) ACCUMULATE_DIAG51           : Accumulates data over for later averaging
!  (6 ) INIT_DIAG51                 : Allocates and zeroes all module arrays 
!  (7 ) CLEANUP_DIAG51              : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by diag51_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f    : Module w/ routines for binary punch file I/O
!  (2 ) dao_mod.f      : Module w/ arrays for DAO met fields
!  (3 ) error_mod.f    : Module w/ NaN and other error check routines
!  (4 ) file_mod.f     : Module w/ file unit numbers and error checks
!  (5 ) grid_mod.f     : Module w/ horizontal grid information
!  (6 ) pbl_mix_mod.f  : Module w/ routines for PBL height & mixing
!  (7 ) pressure_mod.f : Module w/ routines to compute P(I,J,L) 
!  (8 ) time_mod.f     : Module w/ routines to compute date & time
!  (9 ) tracerid_mod.f : Module w/ pointers to tracers & emissions
!
!  ND51 tracer numbers:
!  ============================================================================
!  1 - N_TRACERS : GEOS-CHEM transported tracers            [v/v      ]
!  71            : Pure O3 (not Ox) concentration           [v/v      ]
!  72            : NO concentration                         [v/v      ]
!  73            : NOy concentration                        [v/v      ]
!  74            : OH concentration                         [molec/cm3]
!  75            : NO2 concentration                        [v/v      ]
!  76            : PBL heights                              [m        ]
!  77            : PBL heights                              [levels   ]
!  78            : Air density                              [molec/cm3]
!  79            : 3-D Cloud fractions                      [unitless ]
!  80            : Column optical depths                    [unitless ]
!  81            : Cloud top heights                        [hPa      ]
!  82            : Sulfate aerosol optical depth            [unitless ]
!  83            : Black carbon aerosol optical depth       [unitless ]
!  84            : Organic carbon aerosol optical depth     [unitless ]
!  85            : Accumulation mode seasalt optical depth  [unitless ]
!  86            : Coarse mode seasalt optical depth        [unitless ]
!  87            : Total dust optical depth                 [unitless ]
!  88            : Total seasalt tracer concentration       [unitless ]
!  89 - 92       : RESERVED FOR FUTURE USE
!  93            : Grid box heights                         [m        ]
!  94            : Relative Humidity                        [%        ]
!  95            : Sea level pressure                       [hPa      ]
!  96            : Zonal wind (a.k.a. U-wind)               [m/s      ]
!  97            : Meridional wind (a.k.a. V-wind)          [m/s      ]
!  98            : P(surface) - PTOP                        [hPa      ]
!  99            : Temperature                              [K        ]
!
!  NOTES:
!  (1 ) Rewritten for clarity (bmy, 7/20/04)
!  (2 ) Added extra counters for NO, NO2, OH, O3.  Also all diagnostic counter
!        arrays are 1-D since they only depend on longitude. (bmy, 10/25/04)
!  (3 ) Bug fix: Now get I0 and J0 properly for nested grids (bmy, 11/9/04)
!  (4 ) Now only archive AOD's once per chemistry timestep (bmy, 1/14/05)
!  (5 ) Now references "pbl_mix_mod.f" (bmy, 2/16/05)
!  (6 ) Now save cld frac and grid box heights (bmy, 4/20/05)
!  (7 ) Remove TRCOFFSET since it's always zero  Also now get HALFPOLAR for
!        both GCAP and GEOS grids.  (bmy, 6/28/05)
!  (8 ) Bug fix: do not save SLP if it's not allocated (bmy, 8/2/05)
!  (9 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "diag51_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE 
      
      ! ... except these variables ...
      PUBLIC :: DO_SAVE_DIAG51

      ! ... and these routines
      PUBLIC :: CLEANUP_DIAG51
      PUBLIC :: DIAG51
      PUBLIC :: INIT_DIAG51
      
      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      
      ! Scalars
      LOGICAL              :: DO_SAVE_DIAG51
      INTEGER              :: IOFF,           JOFF,    LOFF
      INTEGER              :: I0,             J0
      INTEGER              :: ND51_N_TRACERS, ND51_TRACERS(100)
      INTEGER              :: ND51_IMIN,      ND51_IMAX
      INTEGER              :: ND51_JMIN,      ND51_JMAX
      INTEGER              :: ND51_LMIN,      ND51_LMAX
      INTEGER              :: ND51_FREQ,      ND51_NI
      INTEGER              :: ND51_NJ,        ND51_NL
      INTEGER              :: HALFPOLAR
      INTEGER, PARAMETER   :: CENTER180=1
      REAL*4               :: LONRES,         LATRES
      REAL*8               :: TAU0,           TAU1
      REAL*8               :: ND51_HR1,       ND51_HR2
      REAL*8               :: ND51_HR_WRITE
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE
      CHARACTER(LEN=255)   :: ND51_OUTPUT_FILE

      ! Arrays
      INTEGER, ALLOCATABLE :: GOOD(:)
      INTEGER, ALLOCATABLE :: GOOD_CT(:)
      INTEGER, ALLOCATABLE :: GOOD_CT_CHEM(:)
      REAL*8,  ALLOCATABLE :: Q(:,:,:,:)

      !=================================================================
      ! Original code from old DIAG51_MOD.  Leave here as a guide to 
      ! figure out when the averaging periods should be and when to
      ! write to disk (bmy, 9/28/04)
      !
      !! For timeseries between 1300 and 1700 LT, uncomment this code:
      !!
      !! Need to write to the bpch file at 12 GMT, since this covers
      !! an entire day over the US grid (amf, bmy, 12/1/00)
      !!
      !INTEGER, PARAMETER   :: NHMS_WRITE = 120000
      !REAL*8,  PARAMETER   :: HR1        = 13d0
      !REAL*8,  PARAMETER   :: HR2        = 17d0
      !CHARACTER(LEN=255)   :: FILENAME   = 'ts1_4pm.bpch'
      !=================================================================
      ! For timeseries between 1000 and 1200 LT, uncomment this code:
      !
      ! Between 10 and 12 has been chosen because the off-polar orbit 
      ! of GOME traverses (westward) through local times between 12 
      ! and 10 over North America, finally crossing the equator at 
      ! 10.30 (local time).
      !
      ! Need to write to the bpch file at 00 GMT, since we will be 
      ! interested in the whole northern hemisphere (pip, 12/1/00)
      !
      !INTEGER, PARAMETER   :: NHMS_WRITE = 000000
      !REAL*8,  PARAMETER   :: HR1        = 10d0
      !REAL*8,  PARAMETER   :: HR2        = 12d0
      !CHARACTER(LEN=255)   :: FILENAME   ='ts10_12pm.bpch'
      !=================================================================

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DIAG51
!
!******************************************************************************
!  Subroutine DIAG51 generates time series (averages from 10am - 12pm LT 
!  or 1pm - 4pm LT) for the US grid area.  Output is to binary punch files.
!  (amf, bey, bdf, pip, bmy, 11/15/99, 9/28/04)
!
!  NOTES:
!  (1 ) Rewritten for clarity (bmy, 7/20/04)
!  (2 ) Added TAU_W as a local variable (bmy, 9/28/04)
!******************************************************************************
!
      ! Local variables
      REAL*8 :: TAU_W

      !=================================================================
      ! DIAG51 begins here!
      !=================================================================
      
      ! Construct array of where local times are between HR1, HR2
      CALL GET_LOCAL_TIME

      ! Accumulate data in the Q array
      CALL ACCUMULATE_DIAG51

      ! Write data to disk at the proper time
      IF ( ITS_TIME_FOR_WRITE_DIAG51( TAU_W ) ) THEN
         CALL WRITE_DIAG51( TAU_W )
      ENDIF

      ! Return to calling program
      END SUBROUTINE DIAG51

!------------------------------------------------------------------------------
 
      SUBROUTINE GET_LOCAL_TIME
!
!******************************************************************************
!  Subroutine GET_LOCAL_TIME computes the local time and returns an array 
!  of points where the local time is between two user-defined limits. 
!  (bmy, 11/29/00, 7/20/04)
!
!  NOTES:
!  (1 ) The 1d-3 in the computation of XLOCTM is to remove roundoff ambiguity 
!        if a the local time should fall exactly on an hour boundary.
!        (bmy, 11/29/00)
!  (2 ) Bug fix: XMID(I) should be XMID(II).  Also updated comments.
!        (bmy, 7/6/01)
!  (3 ) Updated comments (rvm, bmy, 2/27/02)
!  (4 ) Now uses function GET_LOCALTIME of "time_mod.f" (bmy, 3/27/03) 
!  (5 ) Removed reference to CMN (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE TIME_MOD, ONLY : GET_LOCALTIME

#     include "CMN_SIZE"   ! Size parameters

      ! Local variables 
      INTEGER :: I, LT

      !=================================================================
      ! GET_LOCAL_TIME begins here!
      !=================================================================
      DO I = 1, IIPAR

         ! Get local time
         LT = GET_LOCALTIME(I)

         ! GOOD indicates which boxes have local times between HR1 and HR2
         IF ( LT >= ND51_HR1 .and. LT <= ND51_HR2 ) THEN
            GOOD(I) = 1
         ELSE
            GOOD(I) = 0
         ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE GET_LOCAL_TIME

!------------------------------------------------------------------------------

      SUBROUTINE ACCUMULATE_DIAG51
!
!******************************************************************************
!  Subroutine ACCUMULATE_DIAG51 accumulates tracers into the Q array. 
!  (bmy, 8/20/02, 10/3/05)
!
!  NOTES:
!  (1 ) Rewrote to remove hardwiring and for better efficiency.  Added extra
!        diagnostics and updated numbering scheme.  Now scale optical depths
!        to 400 nm (which is usually what QAA(2,*) is.  (bmy, 7/20/04) 
!  (2 ) Now reference GET_ELAPSED_MIN and GET_TS_CHEM from "time_mod.f".  
!        Also now all diagnostic counters are 1-D since they only depend on 
!        longitude. Now only archive NO, NO2, OH, O3 on every chemistry 
!        timestep (i.e. only when fullchem is called). (bmy, 10/25/04)
!  (3 ) Only archive AOD's when it is a chem timestep (bmy, 1/14/05)
!  (4 ) Remove reference to "CMN".  Also now get PBL heights in meters and 
!        model layers from GET_PBL_TOP_m and GET_PBL_TOP_L of "pbl_mix_mod.f".
!        (bmy, 2/16/05)
!  (5 ) Now reference CLDF and BXHEIGHT from "dao_mod.f".  Now save 3-D cloud 
!        fraction as tracer #79 and box height as tracer #93.  Now remove 
!        references to CLMOSW, CLROSW, and PBL from "dao_mod.f". (bmy, 4/20/05)
!  (6 ) Remove TRCOFFSET since it's always zero  Also now get HALFPOLAR for
!        both GCAP and GEOS grids.  (bmy, 6/28/05)
!  (7 ) Now do not save SLP data if it is not allocated (bmy, 8/2/05)
!  (8 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD,      AIRDEN, BXHEIGHT, CLDF 
      USE DAO_MOD,      ONLY : CLDTOPS, OPTD,   RH,       T 
      USE DAO_MOD,      ONLY : UWND,    VWND,   SLP
      USE PBL_MIX_MOD,  ONLY : GET_PBL_TOP_L,   GET_PBL_TOP_m
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TIME_MOD,     ONLY : GET_ELAPSED_MIN, GET_TS_CHEM 
      USE TIME_MOD,     ONLY : TIMESTAMP_STRING
      USE TRACER_MOD,   ONLY : STT, TCVV, ITS_A_FULLCHEM_SIM, N_TRACERS
      USE TRACERID_MOD, ONLY : IDTHNO3, IDTHNO4, IDTN2O5, IDTNOX  
      USE TRACERID_MOD, ONLY : IDTPAN,  IDTPMN,  IDTPPN,  IDTOX   
      USE TRACERID_MOD, ONLY : IDTR4N2, IDTSALA, IDTSALC 

#     include "cmn_fj.h"  ! includes CMN_SIZE
#     include "jv_cmn.h"  ! ODAER
#     include "CMN_O3"    ! FRACO3, FRACNO, SAVEO3, SAVENO2, SAVEHO2, FRACNO2
#     include "CMN_GCTM"  ! SCALE_HEIGHT

      ! Local variables
      LOGICAL, SAVE     :: FIRST = .TRUE.
      LOGICAL, SAVE     :: IS_FULLCHEM, IS_NOx, IS_Ox,   IS_SEASALT
      LOGICAL, SAVE     :: IS_CLDTOPS,  IS_NOy, IS_OPTD, IS_SLP
      LOGICAL           :: IS_CHEM
      INTEGER           :: H, I, J, K, L, M, N
      INTEGER           :: PBLINT,  R, X, Y, W, XMIN
      REAL*8            :: C1, C2, PBLDEC, TEMPBL, TMP, SCALE400nm
      CHARACTER(LEN=16) :: STAMP

      ! Aerosol types (rvm, aad, bmy, 7/20/04)
      INTEGER           :: IND(6) = (/ 22, 29, 36, 43, 50, 15 /)

      !=================================================================
      ! ACCUMULATE_DIAG51 begins here!
      !=================================================================

      ! Set logical flags on first call
      IF ( FIRST ) THEN
         IS_OPTD     = ALLOCATED( OPTD    )
         IS_CLDTOPS  = ALLOCATED( CLDTOPS )
         IS_SLP      = ALLOCATED( SLP     )
         IS_FULLCHEM = ITS_A_FULLCHEM_SIM()
         IS_SEASALT  = ( IDTSALA > 0 .and. IDTSALC > 0 )
         IS_NOx      = ( IS_FULLCHEM .and. IDTNOX  > 0 )
         IS_Ox       = ( IS_FULLCHEM .and. IDTOx   > 0 )
         IS_NOy      = ( IS_FULLCHEM .and. 
     &                   IDTNOX  > 0 .and. IDTPAN  > 0 .and.
     &                   IDTHNO3 > 0 .and. IDTPMN  > 0 .and.
     &                   IDTPPN  > 0 .and. IDTR4N2 > 0 .and.
     &                   IDTN2O5 > 0 .and. IDTHNO4 > 0 ) 
         FIRST       = .FALSE.
      ENDIF

      ! Is it a chemistry timestep?
      IS_CHEM = ( MOD( GET_ELAPSED_MIN(), GET_TS_CHEM() ) == 0 )

      ! Echo info
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 100 ) STAMP
 100  FORMAT( '     - DIAG51: Accumulation at ', a )
      
      !=================================================================
      ! Archive tracers into accumulating array Q 
      !=================================================================

      ! Archive counter array of good points 
      DO X = 1, ND51_NI
         I          = GET_I( X )
         GOOD_CT(X) = GOOD_CT(X) + GOOD(I)
      ENDDO

      ! Archive counter array of good points for chemistry timesteps only
      IF ( IS_CHEM ) THEN
         DO X = 1, ND51_NI
            I               = GET_I( X )
            GOOD_CT_CHEM(X) = GOOD_CT_CHEM(X) + GOOD(I)
         ENDDO
      ENDIF

      ! Accumulate quantities
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( W, N, X, Y, K, I, J, L, TMP, H, R, SCALE400nm ) 
!$OMP+SCHEDULE( DYNAMIC )
      DO W = 1, ND51_N_TRACERS

         ! ND51 Tracer number
         N = ND51_TRACERS(W)

         ! Loop over levels
         DO K = 1, ND51_NL
            L = LOFF + K

         ! Loop over latitudes 
         DO Y = 1, ND51_NJ
            J = JOFF + Y

         ! Loop over longitudes
         DO X = 1, ND51_NI
            I = GET_I( X )

            ! Archive by simulation 
            IF ( N <= N_TRACERS ) THEN

               !--------------------------------------
               ! GEOS-CHEM tracers [v/v]
               !--------------------------------------

               ! Archive afternoon points
               Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                      ( STT(I,J,L,N) * TCVV(N) / 
     &                        AD(I,J,L)    * GOOD(I) )

            ELSE IF ( N == 71 .and. IS_Ox .and. IS_CHEM ) THEN

               !--------------------------------------
               ! Pure O3 [v/v]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------

               ! Accumulate data
               Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &              ( STT(I,J,L,IDTOX) * FRACO3(I,J,L) *
     &                TCVV(IDTOX)      / AD(I,J,L)     * GOOD(I) )

            ELSE IF ( N == 72 .and. IS_NOx .and. IS_CHEM ) THEN

               !--------------------------------------
               ! NO [v/v]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               
               ! Accumulate data
               Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                ( STT(I,J,L,IDTNOX) * FRACNO(I,J,L) *
     &                  TCVV(IDTNOX)      / AD(I,J,L)     * GOOD(I) )

            ELSE IF ( N == 73 .and. IS_NOy ) THEN

               !--------------------------------------
               ! NOy [v/v]
               !--------------------------------------
  
               ! Temp variable for accumulation
               TMP = 0d0
            
               ! NOx
               TMP = TMP + ( TCVV(IDTNOX)        * GOOD(I) *
     &                       STT(I,J,L,IDTNOX)   / AD(I,J,L) )
               ! PAN
               TMP = TMP + ( TCVV(IDTPAN)        * GOOD(I) *
     &                       STT(I,J,L,IDTPAN)   / AD(I,J,L) )

               ! HNO3
               TMP = TMP + ( TCVV(IDTHNO3)       * GOOD(I) *
     &                       STT(I,J,L,IDTHNO3)  / AD(I,J,L) )
            
               ! PMN
               TMP = TMP + ( TCVV(IDTPMN)        * GOOD(I) *
     &                       STT(I,J,L,IDTPMN)   / AD(I,J,L) )

               ! PPN
               TMP = TMP + ( TCVV(IDTPPN)        * GOOD(I) *
     &                       STT(I,J,L,IDTPPN)   / AD(I,J,L) )
 
               ! R4N2
               TMP = TMP + ( TCVV(IDTR4N2)       * GOOD(I) *
     &                       STT(I,J,L,IDTR4N2)  / AD(I,J,L) )
            
               ! N2O5
               TMP = TMP + ( 2d0 * TCVV(IDTN2O5) * GOOD(I) *
     &                       STT(I,J,L,IDTN2O5)  / AD(I,J,L) )
                        
               ! HNO4
               TMP = TMP + ( TCVV(IDTHNO4)       * GOOD(I) *
     &                       STT(I,J,L,IDTHNO4)  / AD(I,J,L) )

               ! Save afternoon points
               Q(X,Y,K,W) = Q(X,Y,K,W) + TMP
    
            ELSE IF ( N == 74 .and. IS_FULLCHEM .and. IS_CHEM ) THEN

               !--------------------------------------
               ! OH [molec/cm3]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------

               ! Accumulate data
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( SAVEOH(I,J,L) * GOOD(I) )
              
            ELSE IF ( N == 75 .and. IS_NOx .and. IS_CHEM ) THEN

               !--------------------------------------
               ! NO2 [v/v]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------     
               Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &              ( STT(I,J,L,IDTNOX)  * FRACNO2(I,J,L) *
     &                TCVV(IDTNOX)       / AD(I,J,L)      * GOOD(I) )
 
            ELSE IF ( N == 76 ) THEN

               !--------------------------------------
               ! PBL HEIGHTS [m] 
               !--------------------------------------
               IF ( K == 1 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                         ( GET_PBL_TOP_m( I, J ) * GOOD(I) )  
               ENDIF

            ELSE IF ( N == 77 ) THEN

               !--------------------------------------
               ! PBL HEIGHTS [layers] 
               !--------------------------------------
               IF ( K == 1 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                         ( GET_PBL_TOP_L( I, J ) * GOOD(I) )
               ENDIF

            ELSE IF ( N == 78 ) THEN

               !--------------------------------------
               ! AIR DENSITY [molec/cm3] 
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &              ( AIRDEN(L,I,J) * XNUMOLAIR * 1d-6 * GOOD(I) )

            ELSE IF ( N == 79 ) THEN

               !--------------------------------------
               ! 3-D CLOUD FRACTION [unitless]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( CLDF(L,I,J) * GOOD(I) )

            ELSE IF ( N == 80 .and. IS_OPTD ) THEN

               !--------------------------------------
               ! COLUMN OPTICAL DEPTH [unitless]
               !--------------------------------------
               Q(X,Y,1,W) = Q(X,Y,1,W) + ( OPTD(L,I,J) * GOOD(I) )

            ELSE IF ( N == 81 .and. IS_CLDTOPS ) THEN

               !--------------------------------------
               ! CLOUD TOP HEIGHTS [mb]
               !--------------------------------------
               IF ( K == 1 ) THEN
                  TMP        = GET_PEDGE( I, J, CLDTOPS(I,J) )
                  Q(X,Y,K,W) = Q(X,Y,K,W) + ( TMP * GOOD(I) )
               ENDIF

            ELSE IF ( N == 82 .and. IS_CHEM ) THEN

               !--------------------------------------
               ! SULFATE AOD @ 400 nm [unitless]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               DO H = 1, NRH
                  
                  ! Scaling factor to 400 nm
                  SCALE400nm = QAA(2,IND(1)+H-1) / QAA(4,IND(1)+H-1) 

                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                         ( ODAER(I,J,L,H) * SCALE400nm * GOOD(I) )
               ENDDO

            ELSE IF ( N == 83 .and. IS_CHEM ) THEN

               !--------------------------------------
               ! BLACK CARBON AOD @ 400 nm [unitless]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               DO R = 1, NRH

                  ! Index for ODAER
                  H          = NRH    + R

                  ! Scaling factor to 400 nm
                  SCALE400nm = QAA(2,IND(2)+R-1) / QAA(4,IND(2)+R-1)

                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                         ( ODAER(I,J,L,H) * SCALE400nm * GOOD(I) )
               ENDDO

            ELSE IF ( N == 84 .and. IS_CHEM ) THEN

               !--------------------------------------
               ! ORG CARBON AOD @ 400 nm [unitless]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               DO R = 1, NRH

                  ! Index for ODAER
                  H          = 2*NRH  + R

                  ! Scaling factor to 400 nm
                  SCALE400nm = QAA(2,IND(3)+R-1) / QAA(4,IND(3)+R-1)

                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                         ( ODAER(I,J,L,H) * SCALE400nm * GOOD(I) )
               ENDDO

            ELSE IF ( N == 85 .and. IS_CHEM ) THEN

               !--------------------------------------
               ! ACCUM SEASALT AOD @ 400 nm [unitless]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               DO R = 1, NRH

                  ! Index for ODAER
                  H          = 3*NRH  + R

                  ! Scaling factor to 400 nm
                  SCALE400nm = QAA(2,IND(4)+R-1) / QAA(4,IND(4)+R-1)

                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                         ( ODAER(I,J,L,H) * SCALE400nm * GOOD(I) ) 
               ENDDO

            ELSE IF ( N == 86 .and. IS_CHEM ) THEN

               !--------------------------------------
               ! COARSE SEASALT AOD 400 nm [unitless]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               DO R = 1, NRH

                  ! Index for ODAER
                  H          = 4*NRH + R

                  ! Scaling factor to 400 nm
                  SCALE400nm = QAA(2,IND(5)+R-1) / QAA(4,IND(5)+R-1)

                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                         ( ODAER(I,J,L,H) * SCALE400nm * GOOD(I) )
               ENDDO

            ELSE IF ( N == 87 .and. IS_CHEM ) THEN               

               !--------------------------------------
               ! TOTAL DUST OPTD @ 400 nm [unitless]
               ! NOTE: Only archive at chem timestep
               !--------------------------------------
               DO R = 1, NDUST 

                  ! Scaling factor to 400 nm
                  SCALE400nm = QAA(2,IND(6)+R-1) / QAA(4,IND(6)+R-1)

                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                       ( ODMDUST(I,J,L,R) * SCALE400nm * GOOD(I) )
               ENDDO

            ELSE IF ( N == 88 .and. IS_SEASALT ) THEN

               !-----------------------------------
               ! TOTAL SEASALT TRACER [v/v]
               !-----------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                      ( STT(I,J,L,IDTSALA) + 
     &                        STT(I,J,L,IDTSALC) ) *
     &                        TCVV(IDTSALA)  / AD(I,J,L) * GOOD(I)

            ELSE IF ( N == 93 ) THEN

               !-----------------------------------
               ! GRID BOX HEIGHTS [m]
               !-----------------------------------               
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( BXHEIGHT(I,J,L) * GOOD(I) )

            ELSE IF ( N == 94 ) THEN

               !-----------------------------------
               ! RELATIVE HUMIDITY [%]
               !-----------------------------------               
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( RH(I,J,L) * GOOD(I) )

            ELSE IF ( N == 95 .and. IS_SLP ) THEN

               !-----------------------------------
               ! SEA LEVEL PRESSURE [hPa]
               !-----------------------------------               
               IF ( K == 1 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) + ( SLP(I,J) * GOOD(I) )
               ENDIF

            ELSE IF ( N == 96 ) THEN

               !-----------------------------------
               ! ZONAL (U) WIND [M/S]
               !-----------------------------------               
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( UWND(I,J,L) * GOOD(I) )

            ELSE IF ( N == 97 ) THEN

               !-----------------------------------
               ! MERIDIONAL (V) WIND [M/S]
               !-----------------------------------               
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( VWND(I,J,L) * GOOD(I) )

            ELSE IF ( N == 98 ) THEN

               !-----------------------------------
               ! SURFACE PRESSURE - PTOP [hPa]
               !-----------------------------------
               IF ( K == 1 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                         ( GET_PEDGE(I,J,K) - PTOP ) * GOOD(I)
               ENDIF

            ELSE IF ( N == 99 ) THEN 

               !-----------------------------------
               ! TEMPERATURE [K]
               !-----------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + ( T(I,J,L) * GOOD(I) )

            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDDO 
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE ACCUMULATE_DIAG51

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_WRITE_DIAG51( TAU_W ) RESULT( ITS_TIME )
!
!******************************************************************************
!  Function ITS_TIME_FOR_WRITE_DIAG51 returns TRUE if it's time to write
!  the ND51 bpch file to disk.  We test the time at the next dynamic
!  timestep so that we can write to disk properly. (bmy, 7/20/04, 9/28/04)
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) TAU_W (REAL*8) : TAU value at time of writing to disk
!
!  NOTES:
!  (1 ) Added TAU_W so to make sure the timestamp is accurate. (bmy, 9/28/04)
!******************************************************************************
!
      ! References to F90 modules
      USE TIME_MOD, ONLY : GET_HOUR, GET_MINUTE, GET_TAU,  
     &                     GET_TAUb, GET_TAUe,   GET_TS_DYN

      ! Arguments
      REAL*8, INTENT(OUT) :: TAU_W

      ! Local variables
      LOGICAL             :: ITS_TIME
      REAL*8              :: TAU, HOUR, DYN

      !=================================================================
      ! ITS_TIME_FOR_WRITE_DIAG51 begins here!
      !=================================================================

      ! Initialize
      ITS_TIME = .FALSE.

      ! Current TAU, Hour, and Dynamic Timestep [hrs]
      TAU      = GET_TAU()
      HOUR     = ( GET_MINUTE() / 60d0 ) + GET_HOUR()
      DYN      = ( GET_TS_DYN() / 60d0 )

      ! If first timestep, return FALSE
      IF ( TAU == GET_TAUb() ) RETURN

      ! If the next dyn timestep is the hour of day
      ! when we have to save to disk, return TRUE
      IF ( MOD( HOUR+DYN, 24d0 ) == ND51_HR_WRITE ) THEN
         ITS_TIME = .TRUE.
         TAU_W    = TAU + DYN
         RETURN
      ENDIF

      ! If the next dyn timestep is the 
      ! end of the run, return TRUE
      IF ( TAU + DYN == GET_TAUe() ) THEN
         ITS_TIME = .TRUE.
         TAU_W    = TAU + DYN
         RETURN
      ENDIF

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_WRITE_DIAG51

!------------------------------------------------------------------------------

      SUBROUTINE WRITE_DIAG51( TAU_W )
!
!******************************************************************************
!  Subroutine WRITE_DIAG51 computes the time-average of quantities between
!  local time limits ND51_HR1 and ND51_HR2 and writes them to a bpch file.
!  Arrays and counters are also zeroed for the next diagnostic interval.
!  (bmy, 12/1/00, 10/3/05)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TAU_W (REAL*8) : TAU value at time of writing to disk 
!
!  NOTES:
!  (1 ) Rewrote to` remove hardwiring and for better efficiency.  Added extra
!        diagnostics and updated numbering scheme. (bmy, 7/20/04) 
!  (2 ) Added TAU_W to the arg list.  Now use TAU_W to set TAU0 and TAU0.
!        Also now all diagnostic counters are 1-D since they only depend on 
!        longitude.  Now only archive NO, NO2, OH, O3 on every chemistry
!        timestep (i.e. only when fullchem is called).  Also remove reference
!        to FIRST. (bmy, 10/25/04)
!  (3 ) Now divide tracers 82-87 (i.e. various AOD's) by GOOD_CT_CHEM since
!        these are only updated once per chemistry timestep (bmy, 1/14/05)
!  (4 ) Now save grid box heights as tracer #93.  Now save 3-D cloud fraction 
!        as tracer #79 (bmy, 4/20/05)
!  (5 ) Remove references to TRCOFFSET because it's always zero (bmy, 6/24/05)
!  (6 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! Reference to F90 modules
      USE BPCH2_MOD,  ONLY : BPCH2,           OPEN_BPCH2_FOR_WRITE
      USE ERROR_MOD,  ONLY : ALLOC_ERR
      USE FILE_MOD,   ONLY : IU_ND51
      USE TIME_MOD,   ONLY : EXPAND_DATE,     GET_NYMD    
      USE TIME_MOD,   ONLY : GET_NHMS,        GET_TAU 
      USE TIME_MOD,   ONLY : TIMESTAMP_STRING
      USE TRACER_MOD, ONLY : N_TRACERS

#     include "CMN_SIZE"  ! Size Parameters

      ! Arguments
      REAL*8, INTENT(IN) :: TAU_W

      ! Local variables
      LOGICAL            :: IS_CHEM
      INTEGER            :: I,   J,  L,  W, N, GMNL, GMTRC
      INTEGER            :: IOS, X, Y, K
      CHARACTER(LEN=16)  :: STAMP
      CHARACTER(LEN=40)  :: CATEGORY
      CHARACTER(LEN=40)  :: UNIT 
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! WRITE_DIAG51 begins here!
      !=================================================================

      ! Replace date tokens in FILENAME
      FILENAME = ND51_OUTPUT_FILE
      CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )
      
      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - DIAG51: Opening file ', a ) 

      ! Open output file
      CALL OPEN_BPCH2_FOR_WRITE( IU_ND51, FILENAME, TITLE )

      ! Set ENDING TAU for this bpch write 
      TAU1 = TAU_W
    
      !=================================================================
      ! Compute time-average of tracers between local time limits
      !=================================================================

      ! Echo info
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 110 ) STAMP
 110  FORMAT( '     - DIAG51: Saving to disk at ', a ) 

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( X, Y, K, W, IS_CHEM )
      DO W = 1, ND51_N_TRACERS
         
         ! Set a flag to denote tracers which are only
         ! accumulated on every chemistry timestep
         IS_CHEM = ( ND51_TRACERS(W) == 71 .or. 
     &               ND51_TRACERS(W) == 72 .or. 
     &               ND51_TRACERS(W) == 74 .or.
     &               ND51_TRACERS(W) == 75 .or.
     &               ND51_TRACERS(W) == 82 .or. 
     &               ND51_TRACERS(W) == 83 .or.
     &               ND51_TRACERS(W) == 84 .or.
     &               ND51_TRACERS(W) == 85 .or.
     &               ND51_TRACERS(W) == 86 .or.
     &               ND51_TRACERS(W) == 87 )

         ! Loop over grid boxes
         DO K = 1, ND51_NL
         DO Y = 1, ND51_NJ
         DO X = 1, ND51_NI

            IF ( IS_CHEM ) THEN 

               ! Avoid division by zero for tracers
               ! which are archived each chem timestep
               IF ( GOOD_CT_CHEM(X) > 0 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) / GOOD_CT_CHEM(X) 
               ELSE
                  Q(X,Y,K,W) = 0d0
               ENDIF

            ELSE

               ! Avoid division by zero for all other tracers
               IF ( GOOD_CT(X) > 0 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) / GOOD_CT(X) 
               ELSE
                  Q(X,Y,K,W) = 0d0
               ENDIF

            ENDIF

         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      !=================================================================
      ! Write each tracer from "timeseries.dat" to the timeseries file
      !=================================================================
      DO W = 1, ND51_N_TRACERS

         ! ND51 tracer number
         N = ND51_TRACERS(W)

         ! Save by simulation
         IF ( N <= N_TRACERS ) THEN

            !---------------------
            ! GEOS-CHEM tracers
            !---------------------
            CATEGORY = 'IJ-AVG-$'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND51_NL
            GMTRC    = N

         ELSE IF ( N == 71 ) THEN

            !---------------------
            ! Pure O3
            !---------------------
            CATEGORY = 'IJ-AVG-$'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND51_NL
            GMTRC    = N_TRACERS + 1

         ELSE IF ( N == 72 ) THEN

            !---------------------
            ! Pure NO [v/v]
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND51_NL
            GMTRC    = 9

         ELSE IF ( N == 73 ) THEN

            !---------------------
            ! NOy 
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND51_NL
            GMTRC    = 3

         ELSE IF ( N == 74 ) THEN

            !---------------------
            ! OH 
            !---------------------
            CATEGORY  = 'CHEM-L=$'
            UNIT      = 'molec/cm3'
            GMNL      = ND51_NL
            GMTRC     = 1

         ELSE IF ( N == 75 ) THEN

            !---------------------
            ! NO2 
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND51_NL
            GMTRC    = 25

         ELSE IF ( N == 76 ) THEN 

            !---------------------
            ! PBL Height [m] 
            !---------------------
            CATEGORY = 'PBLDEPTH'
            UNIT     = 'm'
            GMNL     = 1
            GMTRC    = 1

         ELSE IF ( N == 77 ) THEN

            !---------------------
            ! PBL Height [levels]
            !---------------------
            CATEGORY = 'PBLDEPTH'
            UNIT     = 'levels'
            GMNL     = 1
            GMTRC    = 2

         ELSE IF ( N == 78 ) THEN

            !---------------------
            ! Air Density 
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'molec/cm3'
            GMNL     = ND51_NL
            GMTRC    = 22

         ELSE IF ( N == 79 ) THEN

            !---------------------
            ! 3-D Cloud fractions
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'unitless'
            GMNL     = ND51_NL
            GMTRC    = 19

         ELSE IF ( N == 80 ) THEN

            !---------------------
            ! Column opt depths 
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'unitless'
            GMNL     = 1
            GMTRC    = 20
            
         ELSE IF ( N == 81 ) THEN
        
            !---------------------
            ! Cloud top heights 
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'hPa'
            GMNL     = 1
            GMTRC    = 21

         ELSE IF ( N == 82 ) THEN

            !---------------------
            ! Sulfate AOD
            !---------------------            
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51_NL
            GMTRC    = 6

         ELSE IF ( N == 83 ) THEN

            !---------------------
            ! Black Carbon AOD
            !---------------------            
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51_NL
            GMTRC    = 9

         ELSE IF ( N == 84 ) THEN

            !---------------------
            ! Organic Carbon AOD
            !---------------------            
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51_NL
            GMTRC    = 12
            
         ELSE IF ( N == 85 ) THEN

            !---------------------
            ! SS Accum AOD
            !---------------------            
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51_NL
            GMTRC    = 15

         ELSE IF ( N == 86 ) THEN

            !---------------------
            ! SS Coarse AOD
            !---------------------            
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51_NL
            GMTRC    = 18

         ELSE IF ( N == 87 ) THEN

            !---------------------
            ! Total dust OD
            !---------------------   
            CATEGORY = 'OD-MAP-$'
            UNIT     = 'unitless'
            GMNL     = ND51_NL
            GMTRC    = 4

         ELSE IF ( N == 88 ) THEN

            !---------------------
            ! Total seasalt
            !---------------------            
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND51_NL
            GMTRC    = 24

         ELSE IF ( N == 93 ) THEN

            !---------------------
            ! Grid box heights
            !---------------------            
            CATEGORY = 'BXHGHT-$'
            UNIT     = 'm'
            GMNL     = ND51_NL
            GMTRC    = 1

         ELSE IF ( N == 94 ) THEN

            !---------------------
            ! Relative humidity 
            !---------------------            
            CATEGORY = 'DAO-3D-$'
            UNIT     = '%'
            GMNL     = ND51_NL
            GMTRC    = 11

         ELSE IF ( N == 95 ) THEN

            !---------------------
            ! Sea level prs
            !---------------------            
            CATEGORY = 'DAO-FLDS'
            UNIT     = 'hPa'
            GMNL     = 1
            GMTRC    = 21

         ELSE IF ( N == 96 ) THEN

            !---------------------
            ! U-wind
            !---------------------            
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'm/s'
            GMNL     = ND51_NL
            GMTRC    = 1

         ELSE IF ( N == 97 ) THEN

            !---------------------
            ! V-wind
            !---------------------
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'm/s'
            GMNL     = ND51_NL
            GMTRC    = 2

         ELSE IF ( N == 98 ) THEN

            !---------------------
            ! Psurface - PTOP 
            !---------------------
            CATEGORY = 'PS-PTOP'
            UNIT     = 'hPa'
            GMNL     = 1
            GMTRC    = 1

         ELSE IF ( N == 99 ) THEN

            !---------------------
            ! Temperature
            !---------------------
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'K'
            GMNL     = ND51_NL
            GMTRC    = 3
            
         ELSE

            ! Otherwise skip
            CYCLE

         ENDIF

         !------------------------
         ! Save to bpch file
         !------------------------
         CALL BPCH2( IU_ND51,      MODELNAME,    LONRES,   
     &               LATRES,       HALFPOLAR,    CENTER180, 
     &               CATEGORY,     GMTRC,        UNIT,      
     &               TAU0,         TAU1,         RESERVED,  
     &               ND51_NI,      ND51_NJ,      GMNL,     
     &               ND51_IMIN+I0, ND51_JMIN+J0, ND51_LMIN, 
     &               REAL( Q(1:ND51_NI, 1:ND51_NJ, 1:GMNL, W) ) )
      ENDDO

      ! Echo info
      WRITE( 6, 120 ) TRIM( FILENAME )
 120  FORMAT( '     - DIAG51: Closing file ', a )

      ! Close file
      CLOSE( IU_ND51 )

      !=================================================================
      ! Re-initialize quantities for next diagnostic cycle
      !=================================================================

      ! Echo info
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 130 ) STAMP
 130  FORMAT( '     - DIAG51: Zeroing arrays at ', a )

      ! Set STARTING TAU for the next bpch write
      TAU0 = TAU_W

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( X, Y, K, W )
      DO W = 1, ND51_N_TRACERS
      DO K = 1, ND51_NL
      DO Y = 1, ND51_NJ
      DO X = 1, ND51_NI

         ! Zero accumulating array for tracer
         Q(X,Y,K,W) = 0d0

         ! Zero counters
         IF ( W == 1 .and. K == 1 ) THEN
            GOOD_CT(X)      = 0
            GOOD_CT_CHEM(X) = 0
         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE WRITE_DIAG51

!------------------------------------------------------------------------------

      FUNCTION GET_I( X ) RESULT( I )
!
!******************************************************************************
!  Function GET_I returns the absolute longitude index (I), given the 
!  relative longitude index (X).  (bmy, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) X (INTEGER) : Relative longitude index (used by Q)
!
!  NOTES:
!******************************************************************************
!      
#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: X

      ! Local variables
      INTEGER             :: I
      
      !=================================================================
      ! GET_I begins here!
      !=================================================================

      ! Add the offset to X to get I  
      I = IOFF + X

      ! Handle wrapping around the date line, if necessary
      IF ( I > IIPAR ) I = I - IIPAR

      ! Return to calling program
      END FUNCTION GET_I

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DIAG51( DO_ND51, N_ND51, TRACERS, HR_WRITE, 
     &                        HR1,     HR2,    IMIN,    IMAX,   
     &                        JMIN,    JMAX,   LMIN,    LMAX,  FILE )
!
!******************************************************************************
!  Subroutine INIT_DIAG51 allocates and zeroes all module arrays.  
!  It also gets values for module variables from "input_mod.f". 
!  (bmy, 7/20/04, 6/28/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) DO_ND51  (LOGICAL ) : Switch to turn on ND51 timeseries diagnostic
!  (2 ) N_ND51   (INTEGER ) : Number of ND51 read by "input_mod.f"
!  (3 ) TRACERS  (INTEGER ) : Array w/ ND51 tracer #'s read by "input_mod.f"
!  (4 ) HR_WRITE (REAL*8  ) : GMT hour of day at which to write bpch file
!  (5 ) HR1      (REAL*8  ) : Lower limit of local time averaging bin
!  (6 ) HR2      (REAL*8  ) : Upper limit of local time averaging bin
!  (7 ) IMIN     (INTEGER ) : Min longitude index read by "input_mod.f"
!  (8 ) IMAX     (INTEGER ) : Max longitude index read by "input_mod.f" 
!  (9 ) JMIN     (INTEGER ) : Min latitude index read by "input_mod.f" 
!  (10) JMAX     (INTEGER ) : Min latitude index read by "input_mod.f" 
!  (11) LMIN     (INTEGER ) : Min level index read by "input_mod.f" 
!  (12) LMAX     (INTEGER ) : Min level index read by "input_mod.f" 
!  (13) FILE     (CHAR*255) : ND51 output file name read by "input_mod.f"
!
!  NOTES:
!  (1 ) Diagnostic counter arrays are now only 1-D.  Also add GOOD_CT_CHEM
!        which is the counter array of "good" boxes at each chemistry
!        timesteps.  Now allocate GOOD_CT_CHEM. (bmy, 10/25/04)
!  (2 ) Now get I0 and J0 correctly for nested grid simulations (bmy, 11/9/04)
!  (3 ) Now call GET_HALFPOLAR from "bpch2_mod.f" to get the HALFPOLAR flag 
!        value for GEOS or GCAP grids. (bmy, 6/28/05)
!******************************************************************************
!    
      ! References to F90 modules
      USE BPCH2_MOD,  ONLY : GET_MODELNAME, GET_HALFPOLAR
      USE ERROR_MOD,  ONLY : ALLOC_ERR,   ERROR_STOP
      USE GRID_MOD,   ONLY : GET_XOFFSET, GET_YOFFSET, ITS_A_NESTED_GRID
      USE TIME_MOD,   ONLY : GET_TAUb
      USE TRACER_MOD, ONLY : N_TRACERS
  
#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      LOGICAL,            INTENT(IN) :: DO_ND51
      INTEGER,            INTENT(IN) :: N_ND51, TRACERS(100)
      INTEGER,            INTENT(IN) :: IMIN,   IMAX 
      INTEGER,            INTENT(IN) :: JMIN,   JMAX      
      INTEGER,            INTENT(IN) :: LMIN,   LMAX 
      REAL*8,             INTENT(IN) :: HR1,    HR2
      REAL*8,             INTENT(IN) :: HR_WRITE
      CHARACTER(LEN=255), INTENT(IN) :: FILE

      ! Local variables
      INTEGER                        :: AS
      CHARACTER(LEN=255)             :: LOCATION
      
      !=================================================================
      ! INIT_DIAG51 begins here!
      !=================================================================

      ! Initialize
      LOCATION               = 'INIT_DIAG51 ("diag51_mod.f")'
      ND51_TRACERS(:)        = 0

      ! Get values from "input_mod.f"
      DO_SAVE_DIAG51         = DO_ND51 
      ND51_N_TRACERS         = N_ND51
      ND51_TRACERS(1:N_ND51) = TRACERS(1:N_ND51)
      ND51_HR_WRITE          = HR_WRITE
      ND51_HR1               = HR1
      ND51_HR2               = HR2
      ND51_IMIN              = IMIN
      ND51_IMAX              = IMAX
      ND51_JMIN              = JMIN
      ND51_JMAX              = JMAX
      ND51_LMIN              = LMIN
      ND51_LMAX              = LMAX
      ND51_OUTPUT_FILE       = TRIM( FILE )

      ! Make sure ND51_HR_WRITE is in the range 0-23.999 hrs
      ND51_HR_WRITE = MOD( ND51_HR_WRITE, 24d0 )

      ! Exit if ND51 is turned off 
      IF ( .not. DO_SAVE_DIAG51 ) RETURN

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

      ! Error check ND51_IMIN
      IF ( ND51_IMIN+I0 < 1 .or. ND51_IMIN+I0 > IGLOB ) THEN
         CALL ERROR_STOP( 'Bad ND51_IMIN value!', LOCATION )
      ENDIF

      ! Error check ND51_IMAX
      IF ( ND51_IMAX+I0 < 1 .or. ND51_IMAX+I0 > IGLOB ) THEN
         CALL ERROR_STOP( 'Bad ND51_IMAX value!', LOCATION )
      ENDIF

      ! Compute longitude limits to write to disk
      ! Also handle wrapping around the date line
      IF ( ND51_IMAX > ND51_IMIN ) THEN
         ND51_NI = ( ND51_IMAX - ND51_IMIN ) + 1
      ELSE 
         ND51_NI = ( IIPAR - ND51_IMIN ) + 1 + ND51_IMAX
         WRITE( 6, '(a)' ) 'We are wrapping over the date line!'
      ENDIF

      ! Make sure that ND50_NI <= IIPAR
      IF ( ND51_NI > IIPAR ) THEN
         CALL ERROR_STOP( 'Too many longitudes!', LOCATION )
      ENDIF

      !-----------
      ! Latitude
      !-----------
      
      ! Error check JMIN_AREA
      IF ( ND51_JMIN+J0 < 1 .or. ND51_JMIN+J0 > JGLOB ) THEN
         CALL ERROR_STOP( 'Bad ND51_JMIN value!', LOCATION )
      ENDIF
     
      ! Error check JMAX_AREA
      IF ( ND51_JMAX+J0 < 1 .or.ND51_JMAX+J0 > JGLOB ) THEN
         CALL ERROR_STOP( 'Bad ND51_JMAX value!', LOCATION )
      ENDIF

      ! Compute latitude limits to write to disk (bey, bmy, 3/16/99)
      IF ( ND51_JMAX > ND51_JMIN ) THEN
         ND51_NJ = ( ND51_JMAX - ND51_JMIN ) + 1
      ELSE
         CALL ERROR_STOP( 'ND51_JMAX < ND51_JMIN!', LOCATION )
      ENDIF     
  
      !-----------
      ! Altitude
      !-----------

      ! Error check ND51_LMIN, ND51_LMAX
      IF ( ND51_LMIN < 1 .or. ND51_LMAX > LLPAR ) THEN 
         CALL ERROR_STOP( 'Bad ND51 altitude values!', LOCATION )
      ENDIF

      ! # of levels to save in ND51 timeseries
      IF ( ND51_LMAX >= ND51_LMIN ) THEN  
         ND51_NL = ( ND51_LMAX - ND51_LMIN ) + 1
      ELSE
         CALL ERROR_STOP( 'ND51_LMAX < ND51_LMIN!', LOCATION )
      ENDIF

      !-----------
      ! Offsets
      !-----------
      IOFF      = ND51_IMIN - 1
      JOFF      = ND51_JMIN - 1
      LOFF      = ND51_LMIN - 1

      !-----------
      ! For bpch
      !-----------
      TAU0      = GET_TAUb()
      TITLE     = 'GEOS-CHEM DIAG51 time series'
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

      ! Array denoting where LT is between HR1 and HR2
      ALLOCATE( GOOD( IIPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GOOD' )
      GOOD(:) = 0

      ! Counter of "good" times per day at each grid box
      ALLOCATE( GOOD_CT( ND51_NI ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GOOD_CT' )
      GOOD_CT(:) = 0

      ! Counter of "good" times per day for each chemistry timestep
      ALLOCATE( GOOD_CT_CHEM( ND51_NI ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GOOD_CT_CHEM' )
      GOOD_CT_CHEM(:) = 0

      ! Accumulating array
      ALLOCATE( Q( ND51_NI, ND51_NJ, ND51_NL, ND51_N_TRACERS), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Q' )
      Q(:,:,:,:) = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_DIAG51

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DIAG51
!
!******************************************************************************
!  Subroutine CLEANUP_DIAG51 deallocates all module arrays. 
!  (bmy, 11/29/00, 10/25/04)
!
!  NOTES:
!  (1 ) Now deallocate GOOD_CT_CHEM (bmy, 10/25/04)
!******************************************************************************
! 
      !=================================================================
      ! CLEANUP_DIAG51 begins here!
      !=================================================================
      IF ( ALLOCATED( GOOD         ) ) DEALLOCATE( GOOD         )
      IF ( ALLOCATED( GOOD_CT      ) ) DEALLOCATE( GOOD_CT      )
      IF ( ALLOCATED( GOOD_CT_CHEM ) ) DEALLOCATE( GOOD_CT_CHEM )
      IF ( ALLOCATED( Q            ) ) DEALLOCATE( Q            )

      ! Return to calling program
      END SUBROUTINE CLEANUP_DIAG51

!------------------------------------------------------------------------------

      ! End of module
      END MODULE DIAG51_MOD
