! $Id: diag50_mod.f,v 1.1 2004/09/21 18:04:11 bmy Exp $
      MODULE DIAG50_MOD
!
!******************************************************************************
!  Module DIAG50_MOD contains variables and routines to generate save 
!  timeseries data over the United States where the local time is between 
!  two user-defined limits. (amf, bey, bdf, pip, bmy, 11/30/00, 7/20/04)
!
!  Module Variables:
!  ============================================================================
!  (1 ) COUNT            (INTEGER ) : Counter of timesteps per day
!  (2 ) DO_SAVE_DIAG50   (LOGICAL ) : Flag to turn on DIAG50 timseries
!  (3 ) I0               (INTEGER ) : Lon offset between global & nested grid
!  (4 ) J0               (INTEGER ) : Lat offset between global & nested grid
!  (4 ) IOFF             (INTEGER ) : Offset between relative & absolute lon
!  (5 ) JOFF             (INTEGER ) : Offset between relative & absolute lat
!  (6 ) LOFF             (INTEGER ) : Offset between relative & absolute alt 
!  (7 ) ND50_IMIN        (INTEGER ) : Minimum lat index for DIAG50 region
!  (8 ) ND50_IMAX        (INTEGER ) : Maximum lat index for DIAG50 region
!  (9 ) ND50_JMIN        (INTEGER ) : Minimum lon index for DIAG50 region
!  (10) ND50_JMAX        (INTEGER ) : Maximum lon index for DIAG50 region
!  (11) ND50_LMIN        (INTEGER ) : Minimum alt index for DIAG50 region
!  (12) ND50_LMAX        (INTEGER ) : Minimum alt index for DIAG50 region
!  (13) ND50_NI          (INTEGER ) : Number of longitudes in DIAG50 region 
!  (14) ND50_NJ          (INTEGER ) : Number of latitudes  in DIAG50 region
!  (15) ND50_NL          (INTEGER ) : Number of levels     in DIAG50 region
!  (16) ND50_N_TRACERS   (INTEGER ) : Number of tracers for DIAG50
!  (17) ND50_OUTPUT_FILE (CHAR*255) : Name of output file for timeseries data
!  (18) ND50_TRACERS     (INTEGER ) : Array of DIAG50 tracer numbers
!  (19) Q                (REAL*8  ) : Accumulator array for various quantities
!  (20) TAU0             (REAL*8  ) : Starting TAU used to index the bpch file
!  (21) TAU1             (REAL*8  ) : Ending TAU used to index the bpch file
!  (22) HALFPOLAR        (INTEGER ) : Used for bpch file output
!  (23) CENTER180        (INTEGER ) : Used for bpch file output
!  (24) LONRES           (REAL*4  ) : Used for bpch file output
!  (25) LATRES           (REAL*4  ) : Used for bpch file output
!  (26) MODELNAME        (CHAR*20 ) : Used for bpch file output
!  (27) RESERVED         (CHAR*40 ) : Used for bpch file output
!
!  Module Procedures:
!  ============================================================================
!  (1 ) DIAG50                      : Driver subroutine for 24hr timeseries 
!  (2 ) ACCUMULATE_DIAG50           : Accumulates data for later averaging
!  (3 ) ITS_TIME_FOR_WRITE_DIAG50   : Returns T if it's time to write bpch file
!  (4 ) WRITE_DIAG50                : Writes 24-hr averaged data to a bpch file
!  (5 ) GET_I                       : Converts relative lon index to absolute
!  (5 ) INIT_DIAG50                 : Allocates and zeroes all module arrays 
!  (6 ) CLEANUP_DIAG50              : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by diag50_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f    : Module containing routines for binary punch file I/O
!  (2 ) dao_mod.f      : Module containing arrays for DAO met fields
!  (3 ) error_mod.f    : Module containing NaN and other error check routines
!  (4 ) file_mod.f     : Module containing file unit numbers and error checks
!  (5 ) grid_mod.f     : Module containing horizontal grid information
!  (6 ) pressure_mod.f : Module containing routines to compute P(I,J,L) 
!  (7 ) time_mod.f     : Module containing routines to compute date & time
!  (8 ) tracer_mod.    : Module containing GEOS-CHEM tracer array STT etc.
!  (9 ) tracerid_mod.f : Module containing pointers to tracers & emissions
!
!  ND50 tracer numbers:
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
!  79            : Cloud fractions                          [unitless ]
!  80            : Column optical depths                    [unitless ]
!  81            : Cloud top heights                        [hPa      ]
!  82            : Sulfate aerosol optical depth            [unitless ]
!  83            : Black carbon aerosol optical depth       [unitless ]
!  84            : Organic carbon aerosol optical depth     [unitless ]
!  85            : Accumulation mode seasalt optical depth  [unitless ]
!  86            : Coarse mode seasalt optical depth        [unitless ]
!  87            : Total dust optical depth                 [unitless ]
!  88            : Total seasalt tracer concentration       [unitless ]
!  89 - 93       : RESERVED FOR FUTURE USE
!  94            : Relative humidity                        [%        ]
!  95            : Sea level pressure                       [hPa      ]
!  96            : Zonal wind (a.k.a. U-wind)               [m/s      ]
!  97            : Meridional wind (a.k.a. V-wind)          [m/s      ]
!  98            : P(surface) - PTOP                        [hPa      ]
!  99            : Temperature                              [K        ]
!
!  NOTES:
!  (1 ) Rewritten for clarity and to save extra quantities (bmy, 7/20/04)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "diag50_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE 
      
      ! ... except these variables ...
      PUBLIC :: DO_SAVE_DIAG50

      ! ... and these routines
      PUBLIC :: CLEANUP_DIAG50
      PUBLIC :: DIAG50
      PUBLIC :: INIT_DIAG50
      
      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      
      ! Scalars
      LOGICAL              :: DO_SAVE_DIAG50
      INTEGER              :: COUNT
      INTEGER              :: IOFF,           JOFF   
      INTEGER              :: LOFF,           I0
      INTEGER              :: J0,             ND50_NI
      INTEGER              :: ND50_NJ,        ND50_NL
      INTEGER              :: ND50_N_TRACERS, ND50_TRACERS(100)
      INTEGER              :: ND50_IMIN,      ND50_IMAX
      INTEGER              :: ND50_JMIN,      ND50_JMAX
      INTEGER              :: ND50_LMIN,      ND50_LMAX
      INTEGER, PARAMETER   :: HALFPOLAR=1,    CENTER180=1
      REAL*4               :: LONRES,         LATRES
      REAL*8               :: TAU0,           TAU1
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE
      CHARACTER(LEN=255)   :: ND50_OUTPUT_FILE

      ! Arrays
      REAL*8,  ALLOCATABLE :: Q(:,:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DIAG50
!
!******************************************************************************
!  Subroutine DIAG50 generates 24hr average time series.  Output is to 
!  binary punch file format. (amf, bey, bdf, pip, bmy, 11/15/99, 7/20/04)
!
!  NOTES:
!  (1 ) Rewritten for clarity (bmy, 7/20/04)
!******************************************************************************
!
      !=================================================================
      ! DIAG50 begins here!
      !=================================================================

      ! Accumulate data over a 24-hr period in the Q array
      CALL ACCUMULATE_DIAG50

      ! Write data to disk at the end of the day
      IF ( ITS_TIME_FOR_WRITE_DIAG50() ) CALL WRITE_DIAG50

      ! Return to calling program
      END SUBROUTINE DIAG50

!------------------------------------------------------------------------------

      SUBROUTINE ACCUMULATE_DIAG50
!
!******************************************************************************
!  Subroutine ACCUMULATE_DIAG50 accumulates tracers into the Q array. 
!  (bmy, 8/20/02, 7/20/04)
!
!  NOTES:
!  (1 ) Rewrote to remove hardwiring and for better efficiency.  Added extra
!        diagnostics and updated numbering scheme.  Now scale aerosol & dust
!        optical depths to 400 nm. (rvm, aad, bmy, 7/20/04) 
!******************************************************************************
!
      ! Reference to F90 modules
      USE DAO_MOD,      ONLY : AD,     AIRDEN, BXHEIGHT, CLDTOPS, 
     &                         CLMOSW, CLROSW, OPTD,     RH, T, 
     &                         PBL,    UWND,   VWND,     SLP
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TIME_MOD,     ONLY : TIMESTAMP_STRING
      USE TRACER_MOD,   ONLY : STT, TCVV, ITS_A_FULLCHEM_SIM, N_TRACERS
      USE TRACERID_MOD

#     include "cmn_fj.h"  ! includes CMN_SIZE
#     include "jv_cmn.h"  ! ODAER
#     include "CMN"       ! XTRA2
#     include "CMN_O3"    ! FRACO3, FRACNO, SAVEO3, SAVENO2, SAVEHO2, FRACNO2
#     include "CMN_DIAG"  ! TRCOFFSET
#     include "CMN_GCTM"  ! SCALE_HEIGHT

      ! Local variables
      LOGICAL, SAVE       :: FIRST = .TRUE.
      LOGICAL, SAVE       :: IS_FULLCHEM, IS_NOx, IS_Ox,  IS_SEASALT
      LOGICAL, SAVE       :: IS_CLDTOPS,  IS_NOy, IS_OPTD
      INTEGER             :: H, I, J, K, L, M, N
      INTEGER             :: PBLINT,  R, X, Y, W
      REAL*8              :: C1, C2, PBLDEC, TEMPBL, TMP, SCALE400nm
      CHARACTER(LEN=16)   :: STAMP

      ! Aerosol types (rvm, aad, bmy, 7/20/04)
      INTEGER             :: IND(6) = (/ 22, 29, 36, 43, 50, 15 /)

      !=================================================================
      ! ACCUMULATE_DIAG50 begins here!
      !=================================================================

      ! Set logical flags
      IF ( FIRST ) THEN
         IS_OPTD     = ALLOCATED( OPTD    )
         IS_CLDTOPS  = ALLOCATED( CLDTOPS )
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

      ! Echo time information to the screen
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 100 ) STAMP
 100  FORMAT( '     - DIAG50: Accumulation at ', a )
      
      !=================================================================
      ! Archive tracers into accumulating array Q 
      !=================================================================

      ! Increment counter
      COUNT = COUNT + 1

      ! Accumulate quantities
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( W, N, X, Y, K, I, J, L, TMP, H, R, SCALE400nm ) 
!$OMP+SCHEDULE( DYNAMIC )
      DO W = 1, ND50_N_TRACERS

         ! ND50 Tracer number
         N = ND50_TRACERS(W)

         ! Loop over levels
         DO K = 1, ND50_NL
            L = LOFF + K

         ! Loop over latitudes 
         DO Y = 1, ND50_NJ
            J = JOFF + Y

         ! Loop over altitudes
         DO X = 1, ND50_NI
            I = GET_I( X )

            ! Archive by simulation 
            IF ( N <= N_TRACERS ) THEN

               !--------------------------------------
               ! GEOS-CHEM TRACERS [v/v]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                      ( STT(I,J,L,N) * TCVV(N) / AD(I,J,L) )

            ELSE IF ( N == 71 .and. IS_Ox ) THEN

               !--------------------------------------
               ! PURE O3 CONCENTRATION [v/v]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                      ( STT(I,J,L,IDTOX) * FRACO3(I,J,L) *
     &                        TCVV(IDTOX)      / AD(I,J,L)      )

            ELSE IF ( N == 72 .and. IS_NOx ) THEN

               !--------------------------------------
               ! NO CONCENTRATION [v/v]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                      ( STT(I,J,L,IDTNOX) * FRACNO(I,J,L) *
     &                        TCVV(IDTNOX)      / AD(I,J,L)      )

            ELSE IF ( N == 73 .and. IS_NOy ) THEN

               !--------------------------------------
               ! NOy CONCENTRATION [v/v]
               !--------------------------------------
  
               ! Temp variable for accumulation
               TMP = 0d0
            
               ! NOx
               TMP = TMP + ( TCVV(IDTNOX)       * 
     &                       STT(I,J,L,IDTNOX)  / AD(I,J,L) )

               ! PAN
               TMP = TMP + ( TCVV(IDTPAN)       * 
     &                       STT(I,J,L,IDTPAN)  / AD(I,J,L) )

               ! HNO3
               TMP = TMP + ( TCVV(IDTHNO3)      * 
     &                       STT(I,J,L,IDTHNO3) / AD(I,J,L) )
            
               ! PMN
               TMP = TMP + ( TCVV(IDTPMN)       * 
     &                       STT(I,J,L,IDTPMN)  / AD(I,J,L) )

               ! PPN
               TMP = TMP + ( TCVV(IDTPPN)       * 
     &                       STT(I,J,L,IDTPPN)  / AD(I,J,L) )
 
               ! R4N2
               TMP = TMP + ( TCVV(IDTR4N2)       * 
     &                       STT(I,J,L,IDTR4N2)  / AD(I,J,L) )
            
               ! N2O5
               TMP = TMP + ( 2d0 * TCVV(IDTN2O5) * 
     &                       STT(I,J,L,IDTN2O5)  / AD(I,J,L) )
                        
               ! HNO4
               TMP = TMP + ( TCVV(IDTHNO4)       * 
     &                       STT(I,J,L,IDTHNO4)  / AD(I,J,L) )

               ! Accumulate into Q
               Q(X,Y,K,W) = Q(X,Y,K,W) + TMP
    
            ELSE IF ( N == 74 .and. IS_FULLCHEM ) THEN

               !--------------------------------------
               ! OH CONCENTRATION [molec/cm3]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + SAVEOH(I,J,L)
              
            ELSE IF ( N == 75 .and. IS_NOx ) THEN

               !--------------------------------------
               ! NO2 CONCENTRATION [v/v]
               !--------------------------------------     
               Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                      ( STT(I,J,L,IDTNOX) * FRACNO2(I,J,L) *
     &                        TCVV(IDTNOX)      / AD(I,J,L)       )
 
            ELSE IF ( N == 76 ) THEN

               !--------------------------------------
               ! PBL HEIGHTS [m] 
               !--------------------------------------
#if   defined( GEOS_4 ) 

               ! PBL is in meters, store in PBLM
               IF ( K == 1 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) + PBL(I,J)
               ENDIF

#else
               ! Convert PBL from [hPa] to [m], store in PBLM
               IF ( K == 1 ) THEN
                  
                  ! Surface pressure
                  TMP = GET_PEDGE(I,J,1)

                  ! Convert [hPa] to [m]
                  Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                         ( -SCALE_HEIGHT *
     &                            LOG( ( TMP - PBL(I,J) ) / TMP ) )
               ENDIF

#endif
 
            ELSE IF ( N == 77 ) THEN

               !--------------------------------------
               ! PBL HEIGHTS [layers] 
               !--------------------------------------
               IF ( K == 1 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) + XTRA2(I,J)
               ENDIF

            ELSE IF ( N == 78 ) THEN

               !--------------------------------------
               ! AIR DENSITY [molec/cm3] 
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + 
     &                      ( AIRDEN(L,I,J) * XNUMOLAIR * 1d-6 )

            ELSE IF ( N == 79 ) THEN

               !--------------------------------------
               ! COLUMN CLOUD FRACTION [unitless]
               !--------------------------------------
! Prior to 7/20/04:
! Comment out for now -- need to save for all GEOS models
!               IF ( ALLOCATED( CLMOSW ) .and. ALLOCATED( CLROSW ) ) THEN
!
!                  C1 = CLMOSW(1,I,J)       
!                  C2 = 1.0d0 - CLROSW(1,I,J) 
!
!                  ! Sum over levels
!                  DO Z = 2, ND50_NL
!                     L = LOFF + Z
!                        
!                     IF ( CLMOSW(L,I,J) > CLMOSW(L-1,I,J) ) THEN
!                        C1 = CLMOSW(L,I,J)
!                     ENDIF
!                     
!                     C2 = C2 * ( 1.0d0 - CLROSW(L,I,J) )
!                  ENDDO
!
!                  C1 = 1.0d0 - C1
!                  
!                  ! Save afternoon pionts
!                  Q(X,Y,1,W) = Q(X,Y,1,W) + 
!     &                         ( 1.0d0 - ( C1 * C2 ) ) * GOOD(I)
!               ENDIF

            ELSE IF ( N == 80 .and. IS_OPTD ) THEN

               !--------------------------------------
               ! COLUMN OPTICAL DEPTH [unitless]
               !--------------------------------------
               Q(X,Y,1,W) = Q(X,Y,1,W) + OPTD(L,I,J)

            ELSE IF ( N == 81 .and. IS_CLDTOPS ) THEN

               !--------------------------------------
               ! CLOUD TOP HEIGHTS [mb]
               !--------------------------------------
               IF ( K == 1 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) + GET_PEDGE(I,J,CLDTOPS(I,J))
               ENDIF

            ELSE IF ( N == 82 ) THEN

               !--------------------------------------
               ! SULFATE AOD @ 400 nm [unitless]
               !--------------------------------------
               DO R = 1, NRH
                  
                  ! Scaling factor to 400 nm
                  SCALE400nm = QAA(2,IND(1)+R-1) / QAA(4,IND(1)+R-1) 

                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) + ODAER(I,J,L,R) * SCALE400nm
               ENDDO

            ELSE IF ( N == 83 ) THEN

               !--------------------------------------
               ! BLACK CARBON AOD @ 400 nm [unitless]
               !--------------------------------------
               DO R = 1, NRH

                  ! Index for ODAER
                  H          = NRH + R
                  
                  ! Scaling factor to 400 nm
                  SCALE400nm = QAA(2,IND(2)+R-1) / QAA(4,IND(2)+R-1)

                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) + ODAER(I,J,L,H) * SCALE400nm
               ENDDO

            ELSE IF ( N == 84 ) THEN

               !--------------------------------------
               ! ORG CARBON AOD @ 400 nm [unitless]
               !--------------------------------------
               DO R = 1, NRH

                  ! Index for ODAER
                  H          = 2*NRH + R

                  ! Scaling factor to 400 nm
                  SCALE400nm = QAA(2,IND(3)+R-1) / QAA(4,IND(3)+R-1) 

                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) + ODAER(I,J,L,H) * SCALE400nm
               ENDDO

            ELSE IF ( N == 85 ) THEN

               !--------------------------------------
               ! ACCUM SEASALT AOD @ 400 nm [unitless]
               !--------------------------------------
               DO R = 1, NRH

                  ! Index for ODAER
                  H          = 3*NRH + R

                  ! Scaling factor to 400 nm
                  SCALE400nm = QAA(2,IND(4)+R-1) / QAA(4,IND(4)+R-1)
                 
                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) + ODAER(I,J,L,H) * SCALE400nm
               ENDDO

            ELSE IF ( N == 86 ) THEN

               !--------------------------------------
               ! COARSE SEASALT AOD 400 nm [unitless]
               !--------------------------------------
               DO R = 1, NRH

                  ! Index for ODAER
                  H          = 4*NRH + R

                  ! Scaling factor to 400 nm
                  SCALE400nm = QAA(2,IND(5)+R-1) / QAA(4,IND(5)+R-1)
                  
                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) + ODAER(I,J,L,H) * SCALE400nm
               ENDDO

            ELSE IF ( N == 87 ) THEN
               
               !--------------------------------------
               ! TOTAL DUST OPTD @ 400 nm [unitless]
               !--------------------------------------
               DO R = 1, NDUST

                  ! Scaling factor to 400 nm
                  SCALE400nm = QAA(2,IND(6)+R-1) / QAA(4,IND(6)+R-1)

                  ! Accumulate
                  Q(X,Y,K,W) = Q(X,Y,K,W) + ODMDUST(I,J,L,R)*SCALE400nm
               ENDDO

            ELSE IF ( N == 88 .and. IS_SEASALT ) THEN

               !--------------------------------------
               ! TOTAL SEASALT TRACER [v/v]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) +
     &                      ( STT(I,J,L,IDTSALA) + 
     &                        STT(I,J,L,IDTSALC) ) *
     &                        TCVV(IDTSALA)  / AD(I,J,L) 

            ELSE IF ( N == 94 ) THEN

               !--------------------------------------
               ! RELATIVE HUMIDITY [%]
               !--------------------------------------               
               Q(X,Y,K,W) = Q(X,Y,K,W) + RH(I,J,L)

            ELSE IF ( N == 95 ) THEN

               !--------------------------------------
               ! SEA LEVEL PRESSURE [hPa]
               !--------------------------------------            
               IF ( K == 1 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) + SLP(I,J)
               ENDIF

            ELSE IF ( N == 96 ) THEN

               !--------------------------------------
               ! ZONAL (U) WIND [m/s]
               !--------------------------------------               
               Q(X,Y,K,W) = Q(X,Y,K,W) + UWND(I,J,L)

            ELSE IF ( N == 97 ) THEN

               !--------------------------------------
               ! MERIDIONAL (V) WIND [m/s]
               !--------------------------------------            
               Q(X,Y,K,W) = Q(X,Y,K,W) + VWND(I,J,L)

            ELSE IF ( N == 98 ) THEN

               !--------------------------------------
               ! SURFACE PRESSURE - PTOP [hPa]
               !--------------------------------------
               IF ( K == 1 ) THEN
                  Q(X,Y,K,W) = Q(X,Y,K,W) + ( GET_PEDGE(I,J,K) - PTOP ) 
               ENDIF

            ELSE IF ( N == 99 ) THEN 

               !--------------------------------------
               ! TEMPERATURE [K]
               !--------------------------------------
               Q(X,Y,K,W) = Q(X,Y,K,W) + T(I,J,L)

            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDDO 
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE ACCUMULATE_DIAG50

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_WRITE_DIAG50() RESULT( ITS_TIME )
!
!******************************************************************************
!  Function ITS_TIME_FOR_WRITE_DIAG51 returns TRUE if it's time to write
!  the ND51 bpch file to disk.  We test the time at the next dynamic timestep,
!  so that we can close the file before the end of the run properly.
!  (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE TIME_MOD, ONLY : GET_HOUR, GET_MINUTE, GET_TS_DYN

      ! Local variables
      LOGICAL :: ITS_TIME
      REAL*8  :: HR1, HR2

      !=================================================================
      ! ITS_TIME_FOR_WRITE_DIAG50 begins here!
      !=================================================================

      ! Current hour
      HR1      = GET_HOUR() + ( GET_MINUTE() / 60d0 )

      ! Hour at the next dynamic timestep
      HR2      = HR1        + ( GET_TS_DYN() / 60d0 )

      ! If the next dyn step is the start of a new day, return TRUE
      ITS_TIME = ( INT( HR2 ) == 24 )

      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_WRITE_DIAG50

!------------------------------------------------------------------------------
 
      SUBROUTINE WRITE_DIAG50
!
!******************************************************************************
!  Subroutine WRITE_DIAG50 computes the 24-hr time-average of quantities
!  and saves to bpch file format. (bmy, 12/1/00, 7/20/04)  
!
!  NOTES:
!  (1 ) Rewrote to remove hardwiring and for better efficiency.  Added extra
!        diagnostics and updated numbering scheme. (bmy, 7/20/04) 
!******************************************************************************
!
      ! Reference to F90 modules
      USE BPCH2_MOD
      USE ERROR_MOD,  ONLY : ALLOC_ERR
      USE FILE_MOD,   ONLY : IU_ND50
      USE GRID_MOD,   ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,   ONLY : EXPAND_DATE, GET_NYMD,   GET_NHMS,
     &                       GET_TAU,     GET_TS_DYN, TIMESTAMP_STRING
      USE TRACER_MOD, ONLY : N_TRACERS

#     include "CMN_SIZE"  ! Size Parameters
#     include "CMN_DIAG"  ! TRCOFFSET

      ! Local variables
      LOGICAL            :: FIRST = .TRUE.
      INTEGER            :: I,    J,     L,   W, N  
      INTEGER            :: GMNL, GMTRC, IOS, X, Y, K
      CHARACTER(LEN=16)  :: STAMP
      CHARACTER(LEN=40)  :: CATEGORY
      CHARACTER(LEN=40)  :: UNIT
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! WRITE_DIAG50 begins here!
      !=================================================================

      ! Replace time & date tokens in the filename
      FILENAME = ND50_OUTPUT_FILE
      CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - DIAG50: Opening file ', a ) 

      ! Open output file
      CALL OPEN_BPCH2_FOR_WRITE( IU_ND50, FILENAME, TITLE )

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
!$OMP+PRIVATE( X, Y, K, W )
      DO W = 1, ND50_N_TRACERS
      DO K = 1, ND50_NL
      DO Y = 1, ND50_NJ
      DO X = 1, ND50_NI

         ! Avoid division by zero
         IF ( COUNT > 0 ) THEN
            Q(X,Y,K,W) = Q(X,Y,K,W) / COUNT 
         ELSE
            Q(X,Y,K,W) = 0d0
         ENDIF

      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      !=================================================================
      ! Write each tracer from "timeseries.dat" to the timeseries file
      !=================================================================
      DO W = 1, ND50_N_TRACERS

         ! ND50 tracer number
         N = ND50_TRACERS(W)

         ! Save by simulation
         IF ( N <= N_TRACERS ) THEN

            !---------------------
            ! GEOS-CHEM tracers
            !---------------------
            CATEGORY = 'IJ-AVG-$'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND50_NL
            GMTRC    = N + TRCOFFSET

         ELSE IF ( N == 71 ) THEN

            !---------------------
            ! Pure O3
            !---------------------
            CATEGORY = 'IJ-AVG-$'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND50_NL
            GMTRC    = N_TRACERS + 1 + TRCOFFSET

         ELSE IF ( N == 72 ) THEN

            !---------------------
            ! Pure NO 
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND50_NL
            GMTRC    = 9

         ELSE IF ( N == 73 ) THEN

            !---------------------
            ! NOy 
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND50_NL
            GMTRC    = 3

         ELSE IF ( N == 74 ) THEN

            !---------------------
            ! OH 
            !---------------------
            CATEGORY = 'CHEM-L=$'
            UNIT     = 'molec/cm3'
            GMNL     = ND50_NL
            GMTRC    = 1

         ELSE IF ( N == 75 ) THEN

            !---------------------
            ! NO2 
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''
            GMNL     = ND50_NL
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
            ! PBL Height [layers]
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
            GMNL     = ND50_NL
            GMTRC    = 22

         ELSE IF ( N == 79 ) THEN

            !---------------------
            ! Cloud fractions
            !---------------------
            CATEGORY = 'TIME-SER'
            UNIT     = 'unitless'
            GMNL     = 1
            GMTRC    = 19

         ELSE IF ( N == 80 ) THEN

            !---------------------
            ! Optical depths 
            !---------------------
            CATEGORY  = 'TIME-SER'
            UNIT      = 'unitless'
            GMNL      = 1
            GMTRC     = 20
            
         ELSE IF ( N == 81 ) THEN
        
            !---------------------
            ! Cloud top heights 
            !---------------------
            CATEGORY  = 'TIME-SER'
            GMNL      = 1
            GMTRC     = 21

         ELSE IF ( N == 82 ) THEN

            !---------------------
            ! Sulfate AOD
            !---------------------            
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMNL      = ND50_ NL
            GMTRC     = 6

         ELSE IF ( N == 83 ) THEN

            !---------------------
            ! Black Carbon AOD
            !---------------------            
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMNL      = ND50_NL
            GMTRC     = 9

         ELSE IF ( N == 84 ) THEN

            !---------------------
            ! Organic Carbon AOD
            !---------------------            
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMNL      = ND50_ NL
            GMTRC     = 12
            
         ELSE IF ( N == 85 ) THEN

            !---------------------
            ! SS Accum AOD
            !---------------------            
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMNL      = ND50_NL
            GMTRC     = 15

         ELSE IF ( N == 86 ) THEN

            !---------------------
            ! SS Coarse AOD
            !---------------------            
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMNL      = ND50_NL
            GMTRC     = 18

         ELSE IF ( N == 87 ) THEN

            !---------------------
            ! Total dust OD
            !---------------------   
            CATEGORY  = 'OD-MAP-$'
            UNIT      = 'unitless'
            GMNL      = ND50_NL
            GMTRC     = 4

         ELSE IF ( N == 88 ) THEN

            !----------------------
            ! Total Seasalt
            !----------------------
            CATEGORY = 'TIME-SER'
            UNIT     = ''              ! Let GAMAP pick unit
            GMNL     = ND50_NL
            GMTRC    = 24

         ELSE IF ( N == 94 ) THEN

            !---------------------
            ! Relative humidity
            !---------------------            
            CATEGORY = 'DAO-3D-$'
            UNIT     = '%'
            GMNL     = ND50_NL
            GMTRC    = 11

         ELSE IF ( N == 95 ) THEN

            !---------------------
            ! Sea level prs [hPa]
            !---------------------            
            CATEGORY = 'DAO-FLDS'
            UNIT     = 'hPa'
            GMNL     = 1
            GMTRC    = 21

         ELSE IF ( N == 96 ) THEN

            !---------------------
            ! U-wind [m/s]
            !---------------------            
            CATEGORY = 'DAO-3D-$'
            UNIT     = 'm/s'
            GMNL     = ND50_NL
            GMTRC    = 1

         ELSE IF ( N == 97 ) THEN

            !---------------------
            ! V-wind [m/s]
            !---------------------
            CATEGORY  = 'DAO-3D-$'
            UNIT      = 'm/s'
            GMNL      = ND50_NL
            GMTRC     = 2

         ELSE IF ( N == 98 ) THEN

            !---------------------
            ! Psurface [hPa] 
            !---------------------
            CATEGORY  = 'PS-PTOP'
            UNIT      = 'hPa'
            GMNL      = 1
            GMTRC     = 1

         ELSE IF ( N == 99 ) THEN

            !---------------------
            ! Temperature
            !---------------------
            CATEGORY  = 'DAO-3D-$'
            GMNL      = ND50_NL
            GMTRC     = 3
            
         ELSE

            ! Otherwise skip
            CYCLE

         ENDIF

         !------------------------
         ! Save to bpch file
         !------------------------
         CALL BPCH2( IU_ND50,      MODELNAME,    LONRES,   
     &               LATRES,       HALFPOLAR,    CENTER180, 
     &               CATEGORY,     GMTRC,        UNIT,      
     &               TAU0,         TAU1,         RESERVED,  
     &               ND50_NI,      ND50_NJ,      GMNL,     
     &               ND50_IMIN+I0, ND50_JMIN+J0, ND50_LMIN, 
     &               REAL( Q(1:ND50_NI, 1:ND50_NJ, 1:GMNL, W) ) )
      ENDDO

      ! Echo info
      WRITE( 6, 120 ) TRIM( FILENAME )
 120  FORMAT( '     - DIAG50: Closing file ', a )

      ! Close file
      CLOSE( IU_ND50 )

      !=================================================================
      ! Re-initialize quantities for the next diagnostic cycle
      !=================================================================

      ! Echo info
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 130 ) STAMP
 130  FORMAT( '     - DIAG50: Zeroing arrays at ', a )

      ! Set STARTING TAU for the next bpch write
      TAU0  = GET_TAU() + ( GET_TS_DYN() / 60d0 )

      ! Zero counter
      COUNT = 0
      
      ! Zero accumulating array
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( X, Y, K, W )
      DO W = 1, ND50_N_TRACERS
      DO K = 1, ND50_NL
      DO Y = 1, ND50_NJ
      DO X = 1, ND50_NI
         Q(X,Y,K,W) = 0d0
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      ! Return to calling program
      END SUBROUTINE WRITE_DIAG50

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

      SUBROUTINE INIT_DIAG50( DO_ND50, N_ND50, TRACERS, IMIN, IMAX,    
     &                        JMIN,    JMAX,   LMIN,    LMAX, FILE )
!
!******************************************************************************
!  Subroutine INIT_DIAG50 allocates and zeroes all module arrays.  
!  It also gets values for module variables from "input_mod.f". 
!  (bmy, 7/20/04)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) DO_ND50 (LOGICAL ) : Switch to turn on ND50 timeseries diagnostic
!  (2 ) N_ND50  (INTEGER ) : Number of ND50 read by "input_mod.f"
!  (3 ) TRACERS (INTEGER ) : Array w/ ND50 tracer #'s read by "input_mod.f"
!  (4 ) IMIN    (INTEGER ) : Min longitude index read by "input_mod.f"
!  (5 ) IMAX    (INTEGER ) : Max longitude index read by "input_mod.f" 
!  (6 ) JMIN    (INTEGER ) : Min latitude index read by "input_mod.f" 
!  (7 ) JMAX    (INTEGER ) : Min latitude index read by "input_mod.f" 
!  (8 ) LMIN    (INTEGER ) : Min level index read by "input_mod.f" 
!  (9 ) LMAX    (INTEGER ) : Min level index read by "input_mod.f" 
!  (11) FILE    (CHAR*255) : ND50 output file name read by "input_mod.f"
!
!  NOTES:
!******************************************************************************
!    
      ! References to F90 modules
      USE BPCH2_MOD,  ONLY : GET_MODELNAME
      USE ERROR_MOD,  ONLY : ALLOC_ERR,   ERROR_STOP
      USE GRID_MOD,   ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,   ONLY : GET_TAUb
      USE TRACER_MOD, ONLY : N_TRACERS
  
#     include "CMN_SIZE"

      ! Arguments
      LOGICAL,            INTENT(IN) :: DO_ND50
      INTEGER,            INTENT(IN) :: N_ND50, TRACERS(100)
      INTEGER,            INTENT(IN) :: IMIN,   IMAX 
      INTEGER,            INTENT(IN) :: JMIN,   JMAX      
      INTEGER,            INTENT(IN) :: LMIN,   LMAX 
      CHARACTER(LEN=255), INTENT(IN) :: FILE

      ! Local variables
      INTEGER                        :: AS
      CHARACTER(LEN=255)             :: LOCATION

      !=================================================================
      ! INIT_DIAG50 begins here!
      !=================================================================

      ! Initialize
      LOCATION               = 'INIT_DIAG50 ("diag50_mod.f")'
      ND50_TRACERS(:)        = 0

      ! Get values from "input_mod.f"
      DO_SAVE_DIAG50         = DO_ND50 
      ND50_N_TRACERS         = N_ND50
      ND50_TRACERS(1:N_ND50) = TRACERS(1:N_ND50)
      ND50_IMIN              = IMIN
      ND50_IMAX              = IMAX
      ND50_JMIN              = JMIN
      ND50_JMAX              = JMAX
      ND50_LMIN              = LMIN
      ND50_LMAX              = LMAX
      ND50_OUTPUT_FILE       = TRIM( FILE )

      ! Exit if ND50 is turned off 
      IF ( .not. DO_SAVE_DIAG50 ) RETURN

      !=================================================================
      ! Error check longitude, latitude, altitude limits
      !=================================================================

      ! Get grid offsets between global & nested grid
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !-----------
      ! Longitude
      !-----------

      ! Error check ND50_IMIN
      IF ( ND50_IMIN+I0 < 1 .or. ND50_IMIN+I0 > IGLOB ) THEN
         CALL ERROR_STOP( 'Bad ND50_IMIN value!', LOCATION )
      ENDIF

      ! Error check ND50_IMAX
      IF ( ND50_IMAX+I0 < 1 .or. ND50_IMAX+I0 > IGLOB ) THEN
         CALL ERROR_STOP( 'Bad ND50_IMAX value!', LOCATION )
      ENDIF

      ! Compute longitude limits to write to disk 
      ! Also handle wrapping around the date line
      IF ( ND50_IMAX > ND50_IMIN ) THEN
         ND50_NI = ( ND50_IMAX - ND50_IMIN ) + 1
      ELSE 
         ND50_NI = ( IIPAR - ND50_IMIN ) + 1 + ND50_IMAX
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
      IF ( ND50_JMIN+J0 < 1 .or. ND50_JMIN+J0 > JGLOB ) THEN
         CALL ERROR_STOP( 'Bad ND50_JMIN value!', LOCATION )
      ENDIF
     
      ! Error check JMAX_AREA
      IF ( ND50_JMAX+J0 < 1 .or.ND50_JMAX+J0 > JGLOB ) THEN
         CALL ERROR_STOP( 'Bad ND50_JMAX value!', LOCATION )
      ENDIF

      ! Compute latitude limits to write to disk (bey, bmy, 3/16/99)
      IF ( ND50_JMAX > ND50_JMIN ) THEN
         ND50_NJ = ( ND50_JMAX - ND50_JMIN ) + 1
      ELSE
         CALL ERROR_STOP( 'ND50_JMAX < ND50_JMIN!', LOCATION )
      ENDIF     
  
      !-----------
      ! Altitude
      !-----------

      ! Error check ND50_LMIN, ND50_LMAX
      IF ( ND50_LMIN < 1 .or. ND50_LMAX > LLPAR ) THEN 
         CALL ERROR_STOP( 'Bad ND50 altitude values!', LOCATION )
      ENDIF

      ! # of levels to save in ND50 timeseries
      IF ( ND50_LMAX >= ND50_LMIN ) THEN  
         ND50_NL = ( ND50_LMAX - ND50_LMIN ) + 1
      ELSE
         CALL ERROR_STOP( 'ND50_LMAX < ND50_LMIN!', LOCATION )
      ENDIF

      !-----------
      ! Offsets
      !-----------
      IOFF      = ND50_IMIN - 1
      JOFF      = ND50_JMIN - 1
      LOFF      = ND50_LMIN - 1

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

      !=================================================================
      ! Allocate arrays
      !=================================================================

      ! Accumulator array
      ALLOCATE( Q( ND50_NI, ND50_NJ, ND50_NL, ND50_N_TRACERS), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Q' )
      Q(:,:,:,:) = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_DIAG50

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DIAG50
!
!******************************************************************************
!  Subroutine CLEANUP_DIAG50 deallocates all module arrays. 
!  (bmy, 11/29/00, 7/20/04)
!
!  NOTES:
!******************************************************************************
! 
      !=================================================================
      ! CLEANUP_DIAG50 begins here!
      !=================================================================
      IF ( ALLOCATED( Q ) ) DEALLOCATE( Q )

      ! Return to calling program
      END SUBROUTINE CLEANUP_DIAG50

!------------------------------------------------------------------------------

      ! End of module
      END MODULE DIAG50_MOD
