! $Id: diag52.f,v 1.1 2004/05/03 15:26:53 bmy Exp $
      SUBROUTINE DIAG52
! 
!******************************************************************************
!  Subroutine DIAG52 produces column time series for a geographical domain
!  from the information read in timeseries.dat. Output will be in binary 
!  punch (BPCH) format.  This subroutine is hardwired for ICARTT. 
!  (stu, cas, bmy, 4/22/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE COMODE_MOD,   ONLY : JLOP,        CSPEC
      USE FILE_MOD,     ONLY : IU_ND52
      USE GRID_MOD,     ONLY : GET_XOFFSET, GET_YOFFSET, GET_AREA_CM2
      USE TIME_MOD,     ONLY : DATE_STRING, ITS_A_NEW_DAY, 
     &                         GET_TAU,     TIMESTAMP_STRING
      USE TRACERID_MOD

      IMPLICIT NONE

#     include "CMN_SIZE"        ! Size parameters
#     include "CMN_O3"		! Pure O3, SAVENO2
#     include "CMN"             ! STT
#     include "CMN_TIMES"       ! NI, NJ, etc
#     include "CMN_DIAG"        ! TRCOFFSET (bmy, 3/26/99)

      ! Local variables
      CHARACTER(LEN=255), SAVE :: FILENAME
      CHARACTER(LEN=16)        :: STAMP
      LOGICAL, SAVE            :: FIRST = .TRUE.
      INTEGER, SAVE            :: COUNT
      INTEGER, SAVE            :: NUMBER_OF_ARCHIVAL
      INTEGER                  :: I,     J,     L,        I_MID,  IOS
      INTEGER                  :: XI0,   XJ0,   XL0,      GMTRC,  XI1
      INTEGER                  :: XJ1,   XL1,   IFIRST,   JFIRST
      REAL*8                   :: XTAU1, XTAU2, AREA_CM2, VOL
      REAL*8                   :: O3(IIPAR,JJPAR)
      REAL*8                   :: NO2(IIPAR,JJPAR)
      REAL*8                   :: CH2O(IIPAR,JJPAR)
      REAL*8                   :: CO(IIPAR,JJPAR)
      REAL*4                   :: ARRAY(IIPAR,JJPAR,1)

      ! For binary punch file, version 2.0
      REAL*4                   :: LONRES, LATRES
      INTEGER, PARAMETER       :: HALFPOLAR = 1
      INTEGER, PARAMETER       :: CENTER180 = 1 
      CHARACTER(LEN=20)        :: MODELNAME
      CHARACTER(LEN=40)        :: CATEGORY
      CHARACTER(LEN=40)        :: UNIT
      CHARACTER(LEN=40)        :: RESERVED = ''
      CHARACTER(LEN=80)        :: TITLE
      
      ! External functions
      REAL*8,  EXTERNAL        :: BOXVL

      !=================================================================
      ! DIAG52 begins here!
      !=================================================================

      ! Initialize
      TITLE     = 'GEOS-CHEM column timeseries'
      LONRES    = DISIZE
      LATRES    = DJSIZE
      MODELNAME = GET_MODELNAME()

      !=================================================================
      ! Test if it is a new day.  If it is a new day, open a new file 
      ! ts_column.YYMMDD.bpch and call BPCH2_HDR to write the file hdr
      !=================================================================
      IF ( ITS_A_NEW_DAY() ) THEN
         COUNT              = 0
         NUMBER_OF_ARCHIVAL = ( 24 * 60 ) / FREQ_AREA

         ! Use DATE_STRING to insert the proper date (bmy, 3/27/03)
         FILENAME = 'ts_col.' // DATE_STRING() // '.bpch'

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - DIAG52: Opening file ', a )
        
         ! Open *bpch file and write top-of-file header
         CALL OPEN_BPCH2_FOR_WRITE( IU_ND52, FILENAME, TITLE )
      ENDIF

      !=================================================================
      ! Define variables for this timestep
      !=================================================================

      ! Echo info
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 110 ) STAMP
 110  FORMAT( '     - DIAG52: Saving timeseries at ', a )

      ! Increment counter
      COUNT = COUNT + 1
 
      ! Times
      XTAU1 = GET_TAU()
      XTAU2 = XTAU1

      ! Starting indices
      XI0   = IMIN_AREA
      XJ0   = JMIN_AREA
      XL0   = LMIN_AREA

      ! Ending indices
      XI1   = ( XI0 - 1 ) + NI
      XJ1   = ( XJ0 - 1 ) + NJ
      XL1   = ( XL0 - 1 ) + NL

      ! Get nested-grid offsets
      IFIRST = XI0 + GET_XOFFSET( GLOBAL=.TRUE. )
      JFIRST = XJ0 + GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Compute columns (O3, NO2, CH2O are trop, CO is trop+strat)
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, AREA_CM2, VOL )
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over latitudes
      DO J = 1, JJPAR
            
         ! Grid box surface area
         AREA_CM2 = GET_AREA_CM2( J )

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Loop over levels
            DO L = 1, LLPAR

               ! Box volume [cm3]
               VOL = BOXVL(I,J,L)
            
               ! Tropopsheric columns
               IF ( L < LPAUSE(I,J) ) THEN

                  ! Convert O3 from [molec/cm3] to [molec] and sum
                  O3(I,J)   = O3(I,J)   + CSPEC(JLOP(I,J,L),IDO3)  * VOL

                  ! Convert NO2 from [molec/cm3] to [molec] and sum
                  NO2(I,J)  = NO2(I,J)  + CSPEC(JLOP(I,J,L),IDNO2) * VOL

                  ! Convert CH2O from [kg] to [molec] and sum 
                  CH2O(I,J) = CH2O(I,J) + ( STT(I,J,L,IDTCH2O) *
     &                                         XNUMOL(IDTCH2O) )
               ENDIF

               ! Convert CO from [kg] to [molec] and sum
               CO(I,J) = CO(I,J) + ( STT(I,J,L,IDTCO) * XNUMOL(IDTCO) )
            ENDDO

            ! Divide by grid box area to get [molec/cm2]
            O3(I,J)   = O3(I,J)   / AREA_CM2
            NO2(I,J)  = NO2(I,J)  / AREA_CM2
            CH2O(I,J) = CH2O(I,J) / AREA_CM2
            CO(I,J)   = CO(I,J)   / AREA_CM2

         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Write Tropospheric Column of CO
      !=================================================================      
      CATEGORY = 'INST-COL'
      UNIT     = 'molec/cm2'
      GMTRC    = 4

      ! Save out the NH
      ARRAY( 1:NI, 1:NJ, 1 ) = CO( XI0:XI1, XJ0:XJ1 )

      ! Write to bpch file
      CALL BPCH2( IU_ND52,   MODELNAME, LONRES,    LATRES,    
     &            HALFPOLAR, CENTER180, CATEGORY,  GMTRC,           
     &            UNIT,      XTAU1,     XTAU2,     RESERVED,  
     &            NI,        NJ,        1,         IFIRST,    
     &            JFIRST,    1,         ARRAY(1:NI, 1:NJ, 1) )

            
      !=================================================================
      ! Write tropospheric column of CH2O [molec/cm2]
      !=================================================================      
      CATEGORY = 'INST-COL'
      UNIT     = 'molec/cm2'
      GMTRC    = 20

      ! Save out the NH
      ARRAY( 1:NI, 1:NJ, 1 ) = CH2O( XI0:XI1, XJ0:XJ1 )

      ! Write to bpch file
      CALL BPCH2( IU_ND52,   MODELNAME, LONRES,    LATRES,    
     &            HALFPOLAR, CENTER180, CATEGORY,  GMTRC,           
     &            UNIT,      XTAU1,     XTAU2,     RESERVED,  
     &            NI,        NJ,        1,         IFIRST,    
     &            JFIRST,    1,         ARRAY(1:NI, 1:NJ, 1) )

      !=================================================================
      ! Write Tropospheric Column of O3
      !=================================================================      
      CATEGORY = 'INST-COL'
      UNIT     = 'molec/cm2'
      GMTRC    = 42

      ! Save out the NH
      ARRAY( 1:NI, 1:NJ, 1 ) = O3( XI0:XI1, XJ0:XJ1 )

      ! Write to bpch file
      CALL BPCH2( IU_ND52,   MODELNAME,  LONRES,    LATRES,    
     &            HALFPOLAR, CENTER180,  CATEGORY,  GMTRC,           
     &            UNIT,      XTAU1,      XTAU2,     RESERVED,  
     &            NI,        NJ,         1,         IFIRST,    
     &            JFIRST,    1,          ARRAY(1:NI, 1:NJ, 1) )

      !=================================================================
      ! Write Tropospheric Column of NO2
      !=================================================================      
      CATEGORY = 'INST-COL'
      UNIT     = 'molec/cm2'
      GMTRC    = 45

      ! Save out the NH
      ARRAY( 1:NI, 1:NJ, 1 ) = NO2( XI0:XI1, XJ0:XJ1 )

      ! Write to bpch file
      CALL BPCH2( IU_ND52,   MODELNAME,  LONRES,    LATRES,    
     &            HALFPOLAR, CENTER180,  CATEGORY,  GMTRC,           
     &            UNIT,      XTAU1,      XTAU2,     RESERVED,  
     &            NI,        NJ,         1,         IFIRST,    
     &            JFIRST,    1,          ARRAY(1:NI, 1:NJ, 1) )

      !=================================================================
      ! Close the file if it's the last time in the day
      !=================================================================
      IF ( NUMBER_OF_ARCHIVAL == COUNT ) THEN
         CLOSE( IU_ND52 )

         ! Echo info
         WRITE( 6, 120 ) TRIM( FILENAME )
 120     FORMAT( '     - DIAG52: Closing file : ', a )
      ENDIF

      ! Return to calling program
      END SUBROUTINE DIAG52
