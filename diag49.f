! $Id: diag49.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
      SUBROUTINE DIAG49
! 
!******************************************************************************
!  Subroutine DIAG49 produces time series (instantaneous fields) for a 
!  geographical domain from the information read in timeseries.dat. Output 
!  will be in binary punch (BPCH) format. (bey, bmy, rvm, 4/9/99, 3/27/03)
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
!  (20) Remove FIRSTDIAG49 from the arg list -- this is now a local variable.
!        Remove NYMD from the arg list, now use function DATE_STRING from
!        the new "time_mod.f".  Now use functions GET_TAU, GET_DAY, and
!        TIMESTAMP_STRING from "time_mod.f".  Now use functions GET_XOFFSET, 
!        GET_YOFFSET from "grid_mod.f".  Now adjust diagnostic tracer numbers 
!        for the extra fullchem sulfate tracers.  Now also add UWND and
!        VWND as tracers 71, 72. (bmy, 3/31/03)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DAO_MOD,      ONLY : AD, T, UWND, VWND
      USE FILE_MOD,     ONLY : IU_ND49
      USE GRID_MOD,     ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,     ONLY : DATE_STRING, GET_DAY, 
     &                         GET_TAU,     TIMESTAMP_STRING
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TRACERID_MOD, ONLY : IDTNOX, IDTOX

      IMPLICIT NONE

#     include "CMN_SIZE"        ! Size parameters
#     include "CMN_O3"		! Pure O3, SAVENO2
#     include "CMN"             ! STT, T
#     include "CMN_TIMES"       ! STT_TEMPO
#     include "CMN_DIAG"        ! TRCOFFSET (bmy, 3/26/99)

      ! Local variables
      CHARACTER(LEN=255), SAVE :: FILENAME
      LOGICAL, SAVE            :: FIRSTDIAG49 = .TRUE.
      INTEGER, SAVE            :: DAY_LAST
      INTEGER, SAVE            :: COUNT_IN_DAY
      INTEGER, SAVE            :: NUMBER_OF_ARCHIVAL
      INTEGER                  :: TT, I_MID, IOS
      INTEGER                  :: XTRAC, XI0, XJ0, XL0, GMTRC
      INTEGER                  :: NL2,   XI1, XJ1, XL1 
      INTEGER                  :: I,     J,   L, IFIRST, JFIRST
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
      !
      ! Information for the header of the file
      !=================================================================
      TITLE     = 'GEOS-CTM time series for geographical domain'
      LONRES    = DISIZE
      LATRES    = DJSIZE

      ! Get the proper model name for the binary punch file (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      !=================================================================
      ! Test if it is the first time in the subroutine.
      ! Now also set FIRSTDIAG49 = .FALSE. here instead of in MAIN 
      !=================================================================
      IF ( FIRSTDIAG49 ) THEN
         DAY_LAST    = -1
         FIRSTDIAG49 = .FALSE.
      ENDIF

      !=================================================================
      ! Test if it is a new day.
      !
      ! If it is a new day, open a new file tsYYMMDD.bpch and call 
      ! BPCH2_HDR to write the file header (binary punch file v. 2.0).
      !
      !  The header is made of
      !  - fti (character*40)
      !  - title (character*80)
      !=================================================================
      IF ( GET_DAY() /= DAY_LAST ) THEN
         COUNT_IN_DAY = 0
         NUMBER_OF_ARCHIVAL = (24*60)/FREQ_AREA

         ! Use DATE_STRING to insert the proper date (bmy, 3/27/03)
         FILENAME = 'ts' // DATE_STRING() // '.bpch'

         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - DIAG49: Opening file ', a )
        
         ! Open *bpch file and write top-of-file header
         CALL OPEN_BPCH2_FOR_WRITE( IU_ND49, FILENAME, TITLE )

         DAY_LAST = GET_DAY()
      ENDIF

      !=================================================================
      ! Extract tracer concentration from the global fields)
      ! 
      ! For now, we have only the ND45 diag - "IJ-AVG-$"
      !=================================================================
      IF ( DIAG_AREA == 45 ) THEN

         ! Increment counter
         COUNT_IN_DAY = COUNT_IN_DAY + 1
 
         WRITE( 6, 110 ) TIMESTAMP_STRING()
 110     FORMAT( '     - DIAG49: Saving timeseries at ', a )


         ! loop over tracers
         DO TT = 1, NTRAC_AREA

            ! Tracer number
            XTRAC = TRAC_AREA(TT)

            ! Times
            XTAU1 = GET_TAU()
            XTAU2 = GET_TAU()

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

            ! Zero STT_TEMPO for safety's sake (bmy, 1/7/02)
            STT_TEMPO(:,:,:) = 0d0

            ! CASE statement for tracers
            SELECT CASE ( XTRAC ) 

               !========================================================
               ! Save pure O3 as tracer #32 
               !========================================================
               CASE ( 32 )
                  
                  ! Skip if not a full-chemistry simulation 
                  ! or a single-tracer Ox simulation 
                  IF ( NSRCX /= 3 .and. NSRCX /= 6 ) CYCLE

                  CATEGORY = 'IJ-AVG-$'
                  UNIT     = ''
                  GMTRC    = XTRAC + TRCOFFSET
                  NL2      = NL

                  IF ( IDTOX > 0 ) THEN
                     DO L = 1, LLPAR
                     DO J = 1, JJPAR
                     DO I = 1, IIPAR
            
                        ! Get the fraction of Ox that is pure O3 -- FO3(5)
                        CALL O3COMP(I,J,L)
               
                        ! Archive O3 -- convert from [kg] to [v/v]
                        STT_TEMPO(I,J,L) = STT(I,J,L,IDTOX) * 
     &                       TCVV(IDTOX) / AD(I,J,L) * FO3(5)  
                     ENDDO
                     ENDDO
                     ENDDO
                  ENDIF

               !========================================================
               ! Save NO as tracer #33 (bmy, 3/23/03)
               ! GAMAP thinks NO is tracer #9 for the ND48 diagnostic
               !========================================================
               CASE ( 33 )

                  ! NO is only defined for full chemistry simulations
                  IF ( NSRCX /= 3 ) CYCLE

                  CATEGORY = 'TIME-SER'
                  UNIT     = ''
                  GMTRC    = 9
                  NL2      = NL
               
                  IF ( IDTNOX > 0 ) THEN
                     DO L = 1, LLPAR
                     DO J = 1, JJPAR
                     DO I = 1, IIPAR
            
                        ! Get the fraction of NOx that is pure NO -- FNO(5)
                        CALL O3COMP(I,J,L)
               
                        ! Archive NO -- convert from [kg] to [v/v]
                        STT_TEMPO(I,J,L) = STT(I,J,L,IDTNOX) *
     &                       TCVV(IDTNOX) * FNO(5) / AD(I,J,L)
                     ENDDO
                     ENDDO
                     ENDDO
                  ENDIF

               !========================================================
               ! # 34 -- leave blank for now
               !========================================================
               CASE ( 34 )
                  CYCLE

               !========================================================
               ! Store OH as tracer #35 
               ! GAMAP thinks NO2 is tracer #2 for ND48 diagnostic
               !========================================================
               CASE ( 35 )
                  
                  ! Skip if not a full chemistry run (bmy, 1/7/02)
                  IF ( NSRCX /= 3 ) CYCLE
                  
                  CATEGORY         = 'TIME-SER'
                  UNIT             = ''
                  STT_TEMPO(:,:,:) = SAVEOH(:,:,:)
                  GMTRC            = 2
                  NL2              = NL

               !========================================================
               ! Store NO2 as tracer #36
               ! GAMAP thinks NO2 is tracer #19 for ND48 diagnostic
               !========================================================
               CASE ( 36 )

                  ! Skip if not a full chemistry run (bmy, 1/7/02)
                  IF ( NSRCX /= 3 ) CYCLE

                  CATEGORY         = 'TIME-SER'
                  UNIT             = ''
                  STT_TEMPO(:,:,:) = SAVENO2(:,:,:)
                  GMTRC            = 19
                  NL2              = NL

               !========================================================
               ! Leave #37 - #44 blank for now
               !========================================================
               CASE ( 37:44 )
                  CYCLE

               !========================================================
               ! Store UWND as tracer #71 (GAMAP tracer is #1)
               !========================================================
               CASE ( 71 )
                  CATEGORY         = 'DAO-3D-$'
                  UNIT             = 'm/s'
                  STT_TEMPO(:,:,:) = UWND(:,:,:)
                  GMTRC            = 1
                  NL2              = NL

               !========================================================
               ! Store VWND as tracer #72 (GAMAP tracer is #2)
               !========================================================
               CASE ( 72 )
                  CATEGORY         = 'DAO-3D-$'
                  UNIT             = 'm/s'
                  STT_TEMPO(:,:,:) = VWND(:,:,:)
                  GMTRC            = 2
                  NL2              = NL

               !========================================================
               ! Store Psurface - PTOP as tracer #98
               ! GAMAP thinks pressure is tracer 1 for ND31 diagnostic
               !========================================================
               CASE ( 98 )
                  CATEGORY         = 'PS-PTOP'
                  UNIT             = 'mb'
                  
                  ! Need to loop over w/ I and J (yxw, bmy, 1/30/03)
                  DO J = 1, JJPAR
                  DO I = 1, IIPAR
                     STT_TEMPO(I,J,1) = GET_PEDGE(I,J,1) - PTOP
                  ENDDO
                  ENDDO

                  GMTRC            = 1
                  NL2              = 1
                  XL1              = 1 

               !========================================================
               ! Store temperature as tracer #99
               ! GAMAP thinks temp is tracer #3 for ND68 diagnostic
               !========================================================
               CASE ( 99 )
                  CATEGORY         = 'DAO-3D-$'
                  UNIT             = 'K'
                  STT_TEMPO(:,:,:) = T(:,:,:)
                  GMTRC            = 3
                  NL2              = NL

               !========================================================
               ! Store CTM Tracers [v/v] as tracers 1-24
               !========================================================
               CASE DEFAULT
                  CATEGORY         = 'IJ-AVG-$'
                  UNIT             = ''
                  GMTRC            = XTRAC + TRCOFFSET
                  NL2              = NL
                  
                  ! Convert STT from [kg] to [v/v]
                  STT_TEMPO(:,:,:) = STT(:,:,:,XTRAC) * 
     &                               TCVV(XTRAC) / AD(:,:,:) 

            END SELECT

            !===========================================================
            ! IF (word_wrap) then we are going through the date line
            !===========================================================
            IF ( WORD_WRAP ) THEN      
               I_MID = (IIPAR-XI0) + 1

               XARRAY(1:I_MID,1:NJ,1:NL2) =
     &              STT_TEMPO( XI0:IIPAR, XJ0:XJ1, XL0:XL1 )

               XARRAY(I_MID+1:NI,1:NJ,1:NL2) =
     &              STT_TEMPO( 1:IMAX_AREA, XJ0:XJ1, XL0:XL1 )

            ELSE

               ! Don't wrap around the date line
               XARRAY(1:NI,1:NJ,1:NL2) =
     &              STT_TEMPO( XI0:XI1, XJ0:XJ1, XL0:XL1 )

            ENDIF

            !===========================================================
            ! Write the block data + header for each block
            ! we need to pass twice tau as xtau1, xtau2 for GAMAP
            ! Also write model information (name, resolution, etc) 
            ! to each data block
            !===========================================================
            CALL BPCH2( IU_ND49,   MODELNAME,       LONRES,   
     &                  LATRES,    HALFPOLAR,       CENTER180, 
     &                  CATEGORY,  GMTRC,           UNIT,      
     &                  XTAU1,     XTAU2,           RESERVED,  
     &                  NI,        NJ,              NL2,  
     &                  IFIRST,    JFIRST,          XL0, 
     &                  XARRAY(1:NI, 1:NJ, 1:NL2) )
         ENDDO
      ENDIF
            
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
