! $Id: diag49.f,v 1.6 2004/05/03 14:46:15 bmy Exp $
      SUBROUTINE DIAG49
! 
!******************************************************************************
!  Subroutine DIAG49 produces time series (instantaneous fields) for a 
!  geographical domain from the information read in timeseries.dat. Output 
!  will be in binary punch (BPCH) format. (bey, bmy, rvm, 4/9/99, 4/26/04)
!
!  NOTES:
!  (1 ) Now use F90 syntax for declarations (bmy, 3/26/99)
!  (2 ) Now add TRCOFFSET to the tracer number for special chemistry
!        simulations (bmy, 3/26/99)
!  (3 ) Rename "NAMEDIAG" to "CATEGORY", for compatibility with binary
!       punch file format (bmy, 5/27/99)
!  (4 ) Added surface pressure as tracer #98 and temperature as tracer #99.  
!       Also skip tracers #25 (O3) and #26 (NO) and increase the length
!        of FILENAME to 255 characters. (bmy, 1/5/00)  
!  (5 ) Now archive NO2 as tracer #29 and NO as tracer #26.  Also fixed
!        a few unit conversion problems and now we compute the ending
!        indices XI1, XJ1, XL1 properly.  Also added a few cosmetic changes 
!        and updated comments. (rvm, bmy, 5/10/00)
!  (6 ) Now use function NYMD_STRING from "time_mod.f" to generate a
!        Y2K compliant string for all data sets.  Also include code for
!        GEOS-2 and GEOS-3 data sets. (bmy, 6/22/00)
!  (7 ) Reference F90 module "bpch2_mod.f" which contains routines BPCH2_HDR, 
!        BPCH2, and GET_MODELNAME for writing data to binary punch files. 
!        (bmy, 6/22/00)
!  (8 ) Only add TRCOFFSET to actual CTM tracers, and not to other quantities
!        like surface pressure or NO concentrations.   Also skip over tracers
!        #30, #42, #43, and #44.  (bmy, 7/17/00) 
!  (9 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
!  (10) Bug fix: put a CYCLE statement for CASE(27) (amf, bmy, 7/26/01)
!  (11) Replace PW(:,:) with P(:,:) (bmy, 10/3/01)
!  (12) Removed obsolete code from 9/01 and 10/01 (bmy, 10/23/01)
!  (13) Now pass XI0+I0 as IFIRST and XJ0+J0 as JFIRST in call to BPCH2.
!        This will take care of window regions smaller than the globe.
!        (yxw, bmy, 12/18/01)
!  (14) Zero STT_TEMPO for safety's sake.  Also make STT_TEMPO of size
!        (IIPAR,JJPAR,LLPAR).  Also skip OH and NO2 if we are not running
!        the SMVGEAR chemistry simulation. (bmy, 1/7/02)
!  (15) Eliminate obsolete code from 1/02 (bmy, 2/27/02)
!  (16) Now reference IU_ND49 and from "file_mod.f".  Now use IU_ND49 instead 
!        of IUT as the file unit #.  Also call routine OPEN_BPCH2_FOR_WRITE to
!        start writing to the output file. (bmy, 6/26/02)
!  (17) Replaced references to P(I,J) with call to GET_PEDGE(I,J,1) from 
!        "pressure_mod.f".  Also eliminated obsolete, commented-out code
!         from 6/02. (dsa, bdf, bmy, 8/20/02)
!  (18) Now reference AD & T from "dao_mod.f".  Now reference IDTNOX and IDTOX
!        from "tracerid_mod.f". (bmy, 11/6/02)
!  (19) Bug fix: we need to loop over I,J when calling GET_PEDGE, since this
!        is a function and not an array.  (yxw, bmy, 1/30/03)
!  (20) Change tracer #'s for sulfate tracers beyond 24.  Also updated
!        comments (rjp, bmy, 3/23/03)
!  (21) Remove FIRSTDIAG49 from the arg list -- this is now a local variable.
!        Remove NYMD from the arg list, now use function DATE_STRING from
!        the new "time_mod.f".  Now use functions GET_TAU, GET_DAY, and
!        TIMESTAMP_STRING from "time_mod.f".  Now use functions GET_XOFFSET, 
!        GET_YOFFSET from "grid_mod.f".  Now adjust diagnostic tracer numbers 
!        for the extra fullchem sulfate tracers.  Now also add UWND and
!        VWND as tracers 71, 72. (bmy, 3/31/03)
!  (22) LINUX has a problem putting a function call w/in a WRITE statement.
!        Now save output from TIMESTAMP_STRING to STAMP and print that.
!        (bmy, 9/29/03)
!  (23) Readjust tracer numbers for other fields to accommodate 41 regular 
!        CTM tracers -- QUICK FIX.  (bmy, 4/16/04)
!  (24) Added parallel DO-loops.  Replaced calls to O3COMP by references to 
!        the arrays FRACO3 and FRACNO.  Added aerosol optical depths to
!        the CASE statement.  Also save out total dust and seasalt mass.
!        (bmy, 4/26/04)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DAO_MOD,      ONLY : AD, T, UWND, VWND
      USE FILE_MOD,     ONLY : IU_ND49
      USE GRID_MOD,     ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,     ONLY : DATE_STRING, ITS_A_NEW_DAY, 
     &                         GET_TAU,     TIMESTAMP_STRING
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TRACERID_MOD

      IMPLICIT NONE

#     include "cmn_fj.h"        ! FAST-J stuff, includes CMN_SIZE
#     include "jv_cmn.h"        ! ODAER
#     include "CMN_O3"		! Pure O3, SAVENO2
#     include "CMN"             ! STT, T
#     include "CMN_TIMES"       ! STT_TEMPO
#     include "CMN_DIAG"        ! TRCOFFSET (bmy, 3/26/99)

      ! Local variables
      CHARACTER(LEN=255), SAVE :: FILENAME
      CHARACTER(LEN=16)        :: STAMP
      LOGICAL, SAVE            :: FIRSTDIAG49 = .TRUE.
      INTEGER, SAVE            :: DAY_LAST
      INTEGER, SAVE            :: COUNT_IN_DAY
      INTEGER, SAVE            :: NUMBER_OF_ARCHIVAL
      INTEGER                  :: TT,  I_MID, IOS,    XTRAC,  XI0 
      INTEGER                  :: XJ0, XL0,   GMTRC,  NL2,    XI1
      INTEGER                  :: XJ1, XL1,   I,      J,      L
      INTEGER                  :: N,   R,     IFIRST, JFIRST, RH
      REAL*4                   :: XARRAY(IIPAR, JJPAR, LLPAR)
      REAL*8                   :: XTAU1, XTAU2
      REAL*8                   :: STT_TEMPO(IIPAR, JJPAR, LLPAR)  

      ! For binary punch file, version 2.0
      REAL*4                   :: LONRES, LATRES
      INTEGER, PARAMETER       :: HALFPOLAR = 1
      INTEGER, PARAMETER       :: CENTER180 = 1 
      CHARACTER(LEN=20)        :: MODELNAME
      CHARACTER(LEN=40)        :: CATEGORY
      CHARACTER(LEN=40)        :: UNIT
      CHARACTER(LEN=40)        :: RESERVED = ''
      CHARACTER(LEN=80)        :: TITLE
      
      !=================================================================
      ! DIAG49 begins here!
      !=================================================================

      ! Initialize
      TITLE     = 'GEOS-CHEM time series for geographical domain'
      LONRES    = DISIZE
      LATRES    = DJSIZE
      MODELNAME = GET_MODELNAME()

      !=================================================================
      ! If it's a new day, open a new BPCH file and write file header
      !=================================================================
      IF ( ITS_A_NEW_DAY() ) THEN
         COUNT_IN_DAY = 0
         NUMBER_OF_ARCHIVAL = ( 24 * 60 ) / FREQ_AREA

         ! Filename
         FILENAME = 'ts' // DATE_STRING() // '.bpch'

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - DIAG49: Opening file ', a )
        
         ! Open bpch file and write top-of-file header
         CALL OPEN_BPCH2_FOR_WRITE( IU_ND49, FILENAME, TITLE )
      ENDIF

      !=================================================================
      ! Define time and grid information
      !=================================================================

      ! Increment counter
      COUNT_IN_DAY = COUNT_IN_DAY + 1
 
      ! Echo info
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 110 ) STAMP
 110  FORMAT( '     - DIAG49: Saving timeseries at ', a )

      ! Times
      XTAU1 = GET_TAU()
      XTAU2 = XTAU1

      ! Starting indices
      XI0   = IMIN_AREA
      XJ0   = JMIN_AREA
      XL0   = LMIN_AREA

      ! Ending indices
      ! We must subtract 1 from XI0, XJ0, XL0 (bmy, 5/11/00)
      XI1   = ( XI0 - 1 ) + NI
      XJ1   = ( XJ0 - 1 ) + NJ
      XL1   = ( XL0 - 1 ) + NL

      ! Get nested-grid offsets
      IFIRST = XI0 + GET_XOFFSET( GLOBAL=.TRUE. )
      JFIRST = XJ0 + GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Loop over tracers
      !=================================================================
      DO TT = 1, NTRAC_AREA

         ! Tracer number
         XTRAC = TRAC_AREA(TT)

         ! Zero STT_TEMPO for safety's sake (bmy, 1/7/02)
         STT_TEMPO(:,:,:) = 0d0

         ! CASE statement for tracers
         SELECT CASE ( XTRAC ) 

            !===========================================================
            ! Save pure O3 as tracer #71 (bmy, 4/16/04)
            !===========================================================
            CASE ( 71 )
                  
               ! Skip if not a full-chemistry simulation 
               ! or a single-tracer Ox simulation 
               IF ( NSRCX /= 3 .and. NSRCX /= 6 ) CYCLE

               CATEGORY = 'IJ-AVG-$'
               UNIT     = ''
               NL2      = NL
               GMTRC    = NNPAR + 1 + TRCOFFSET

               IF ( IDTOX > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
                  DO L = 1, LLPAR
                  DO J = 1, JJPAR
                  DO I = 1, IIPAR
            
                     ! Archive pure O3 -- convert from [kg] to [v/v]
                     STT_TEMPO(I,J,L) = STT(I,J,L,IDTOX) * 
     &                    TCVV(IDTOX) / AD(I,J,L) * FRACO3(I,J,L) 

                  ENDDO
                  ENDDO
                  ENDDO
!$OMP END PARALLEL DO
               ENDIF

            !===========================================================
            ! Save NO as tracer #72 (bmy, 4/16/04)
            !===========================================================
            CASE ( 72 )

               ! NO is only defined for full chemistry simulations
               IF ( NSRCX /= 3 ) CYCLE

               CATEGORY = 'TIME-SER'
               UNIT     = ''
               NL2      = NL
               GMTRC    = 9
               
               IF ( IDTNOX > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
                  DO L = 1, LLPAR
                  DO J = 1, JJPAR
                  DO I = 1, IIPAR
 
                        ! Archive NO -- convert from [kg] to [v/v]
                        STT_TEMPO(I,J,L) = STT(I,J,L,IDTNOX) *
     &                       TCVV(IDTNOX) * FRACNO(I,J,L) / AD(I,J,L)
                  ENDDO
                  ENDDO
                  ENDDO
!$OMP END PARALLEL DO
               ENDIF

            !===========================================================
            ! Skip over NOy (tracer #73)  (bmy, 4/16/04)
            !===========================================================
            CASE ( 73 )
               CYCLE

            !===========================================================
            ! Store OH as tracer #74 (bmy, 4/16/04)
            !===========================================================
            CASE ( 74 )
                  
               ! Skip if not a full chemistry run (bmy, 1/7/02)
               IF ( NSRCX /= 3 ) CYCLE
                  
               CATEGORY = 'TIME-SER'
               UNIT     = ''
               NL2      = NL
               GMTRC    = 2

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  STT_TEMPO(I,J,L) = SAVEOH(I,J,L)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

            !===========================================================
            ! Store NO2 as tracer #75 (bmy, 4/16/04)
            !===========================================================
            CASE ( 75 )

               ! Skip if not a full chemistry run (bmy, 1/7/02)
               IF ( NSRCX /= 3 ) CYCLE

               CATEGORY = 'TIME-SER'
               UNIT     = ''
               NL2      = NL
               GMTRC    = 19

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L)
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  STT_TEMPO(I,J,L) = SAVENO2(I,J,L)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

            !===========================================================
            ! Leave #76 - #81 blank for now (bmy, 4/16/04)
            !===========================================================
            CASE ( 76:81 )
               CYCLE

            !===========================================================
            ! #82: Sulfate aerosol optical depth [unitless]
            !===========================================================
            CASE ( 82 )
               CATEGORY = 'OD-MAP-$'
               UNIT     = 'unitless'
               GMTRC    = 6
               NL2      = NL
               N        = 1

               DO R = 1, NRH
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR

                  ! Index for type of aerosol and RH value
                  RH = ( (N-1) * NRH ) + R

                  ! Save into STT_TEMPO
                  STT_TEMPO(I,J,L) = STT_TEMPO(I,J,L) + ODAER(I,J,L,RH)
               ENDDO
               ENDDO
               ENDDO
               ENDDO

            !===========================================================
            ! #83: Black Carbon aerosol optical depth [unitless]
            !===========================================================
            CASE ( 83 )
               CATEGORY = 'OD-MAP-$'
               UNIT     = 'unitless'
               GMTRC    = 9
               NL2      = NL
               N        = 2

               DO R = 1, NRH
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR

                  ! Index for type of aerosol and RH value
                  RH = ( (N-1) * NRH ) + R

                  ! Save into STT_TEMPO
                  STT_TEMPO(I,J,L) = STT_TEMPO(I,J,L) + ODAER(I,J,L,RH)
               ENDDO
               ENDDO
               ENDDO
               ENDDO

            !===========================================================
            ! #84: Organic Carbon aerosol optical depth (unitless)
            !===========================================================
            CASE ( 84 )
               CATEGORY = 'OD-MAP-$'
               UNIT     = 'unitless'
               GMTRC    = 12
               NL2      = NL
               N        = 3

               DO R = 1, NRH
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR

                  ! Index for type of aerosol and RH value
                  RH = ( (N-1) * NRH ) + R

                  ! Save into STT_TEMPO
                  STT_TEMPO(I,J,L) = STT_TEMPO(I,J,L) + ODAER(I,J,L,RH)
               ENDDO
               ENDDO
               ENDDO
               ENDDO

            !===========================================================
            ! #85: Accum Sea Salt aerosol optical depth (unitless)
            !===========================================================
            CASE ( 85 )
               CATEGORY = 'OD-MAP-$'
               UNIT     = 'unitless'
               GMTRC    = 15
               NL2      = NL
               N        = 4

               DO R = 1, NRH
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR

                  ! Index for type of aerosol and RH value
                  RH = ( (N-1) * NRH ) + R

                  ! Save into STT_TEMPO
                  STT_TEMPO(I,J,L) = STT_TEMPO(I,J,L) + ODAER(I,J,L,RH)
               ENDDO
               ENDDO
               ENDDO
               ENDDO

            !===========================================================
            ! #86: Coarse Sea Salt aerosol optical depth (unitless)
            !===========================================================
            CASE ( 86 )
               CATEGORY = 'OD-MAP-$'
               UNIT     = 'unitless'
               GMTRC    = 18
               NL2      = NL
               N        = 5

               DO R = 1, NRH
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR

                  ! Index for type of aerosol and RH value
                  RH = ( (N-1) * NRH ) + R

                  ! Save into STT_TEMPO
                  STT_TEMPO(I,J,L) = STT_TEMPO(I,J,L) + ODAER(I,J,L,RH)
               ENDDO
               ENDDO
               ENDDO
               ENDDO

            !===========================================================
            ! #87: Total dust tracer (all size bins) [v/v]
            !===========================================================
            CASE ( 87 )

               ! Skip if DUST isn't defined
               IF ( IDTDST1 + IDTDST2 + IDTDST3 + IDTDST4 > 0 ) THEN

                  CATEGORY = 'TIME-SER'
                  UNIT     = ''          ! Let GAMAP pick unit
                  GMTRC    = 23
                  NL2      = NL

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
                  DO L = 1, LLPAR
                  DO J = 1, JJPAR
                  DO I = 1, IIPAR
                     STT_TEMPO(I,J,L) = ( STT(I,J,L,IDTDST1) + 
     &                                    STT(I,J,L,IDTDST2) +
     &                                    STT(I,J,L,IDTDST3) + 
     &                                    STT(I,J,L,IDTDST4) ) *
     &                                  TCVV(IDTDST1) / AD(I,J,L) 

                  ENDDO
                  ENDDO
                  ENDDO
!$OMP END PARALLEL DO
               ENDIF

            !===========================================================
            ! #88: Total seasalt tracer (accum + coarse) [v/v]
            !===========================================================
            CASE ( 88 )

               ! Skip if SEASALT isn't defined
               IF ( IDTSALA + IDTSALC > 0 ) THEN

                  CATEGORY = 'TIME-SER'
                  UNIT     = ''            ! Let GAMAP pick unit
                  GMTRC    = 24
                  NL2      = NL

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
                  DO L = 1, LLPAR
                  DO J = 1, JJPAR
                  DO I = 1, IIPAR
                     STT_TEMPO(I,J,L) = ( STT(I,J,L,IDTSALA) + 
     &                                    STT(I,J,L,IDTSALC) ) *
     &                                  TCVV(IDTSALA) / AD(I,J,L) 
                  ENDDO
                  ENDDO
                  ENDDO
!$OMP END PARALLEL DO
               ENDIF
               
            !===========================================================
            ! Leave #89-95 blank for now
            !===========================================================
            CASE ( 89:95 )
               CYCLE

            !========================================================
            ! Store UWND as tracer #96 (bmy, 4/16/04)
            !========================================================
            CASE ( 96 )
               CATEGORY = 'DAO-3D-$'
               UNIT     = 'm/s'
               GMTRC    = 1
               NL2      = NL

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L)
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  STT_TEMPO(I,J,L) = UWND(I,J,L)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

            !===========================================================
            ! Store VWND as tracer #97 (bmy, 4/16/04)
            !===========================================================
            CASE ( 97 )
               CATEGORY = 'DAO-3D-$'
               UNIT     = 'm/s'
               GMTRC    = 2
               NL2      = NL

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  STT_TEMPO(I,J,L) = VWND(I,J,L)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO                   

            !===========================================================
            ! Store Psurface - PTOP as tracer #98 (bmy, 4/16/04)
            !===========================================================
            CASE ( 98 )
               CATEGORY = 'PS-PTOP'
               UNIT     = 'hPa'
               NL2      = 1
               XL1      = 1 
               GMTRC    = 1
                  
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  STT_TEMPO(I,J,1) = GET_PEDGE(I,J,1) - PTOP
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

            !===========================================================
            ! Store temperature as tracer #99 (bmy, 4/16/04)
            !===========================================================
            CASE ( 99 )
               CATEGORY = 'DAO-3D-$'
               UNIT     = 'K'
               NL2      = NL
               GMTRC    = 3

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L ) 
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  STT_TEMPO(I,J,L) = T(I,J,L)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

            !===========================================================
            ! Store CTM Tracers [v/v] as tracers 1-41
            !===========================================================
            CASE DEFAULT
               CATEGORY = 'IJ-AVG-$'
               UNIT     = ''
               GMTRC    = XTRAC + TRCOFFSET
               NL2      = NL
                  
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  STT_TEMPO(I,J,L) = STT(I,J,L,XTRAC) *
     &                               TCVV(XTRAC) / AD(I,J,L) 
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

         END SELECT

         !==============================================================
         ! IF (word_wrap) then we are going through the date line
         !==============================================================
         IF ( WORD_WRAP ) THEN      
            I_MID = (IIPAR-XI0) + 1

            XARRAY(1:I_MID,1:NJ,1:NL2) =
     &           STT_TEMPO( XI0:IIPAR, XJ0:XJ1, XL0:XL1 )

            XARRAY(I_MID+1:NI,1:NJ,1:NL2) =
     &           STT_TEMPO( 1:IMAX_AREA, XJ0:XJ1, XL0:XL1 )

         ELSE

            ! Don't wrap around the date line
            XARRAY(1:NI,1:NJ,1:NL2) =
     &           STT_TEMPO( XI0:XI1, XJ0:XJ1, XL0:XL1 )

         ENDIF

         !==============================================================
         ! Write the block data + header for each block
         ! we need to pass twice tau as xtau1, xtau2 for GAMAP
         ! Also write model information (name, resolution, etc) 
         ! to each data block
         !==============================================================
         CALL BPCH2( IU_ND49,   MODELNAME,       LONRES,   
     &               LATRES,    HALFPOLAR,       CENTER180, 
     &               CATEGORY,  GMTRC,           UNIT,      
     &               XTAU1,     XTAU2,           RESERVED,  
     &               NI,        NJ,              NL2,  
     &               IFIRST,    JFIRST,          XL0, 
     &               XARRAY(1:NI, 1:NJ, 1:NL2) )
      ENDDO
            
      !=================================================================
      ! close the file if it is the last time in the day
      ! i.e. if (24*60) / archival frequency = icount_in_day
      !=================================================================
      IF ( NUMBER_OF_ARCHIVAL == COUNT_IN_DAY ) THEN
         CLOSE( IU_ND49 )
 
         WRITE( 6, 120 ) TRIM( FILENAME )
 120     FORMAT( '     - DIAG49: Closing file : ', a )
      ENDIF

      ! Return to calling program
      END SUBROUTINE DIAG49
