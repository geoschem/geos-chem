! $Id: diag50.f,v 1.3 2003/12/05 21:13:59 bmy Exp $
      SUBROUTINE DIAG50
!
!******************************************************************************
!  Subroutine DIAG50 produces time series (24-h averages) for a geographical 
!  domain from the information read in timeseries.dat Output is in binary 
!  punch (BPCH) file format. (bey, bmy, 6/10/99, 9/29/03)
!
!  NOTES:
!  (1 ) Now use F90 syntax for declarations (bmy, 4/9/99)
!  (2 ) Now add TRCOFFSET to the tracer number for special chemistry
!        simulations (bmy, 3/26/99)
!  (3 ) Now use HR1_NO and HR2_NO from CMN_SETUP as local time limits
!        for NO archival (bmy, 4/9/99) 
!  (4 ) Rename "NAMEDIAG" to "CATEGORY", for compatibility with binary
!        punch file format (bmy, 5/27/99)
!  (5 ) Cosmetic changes.  Also add support for GEOS-2. (bmy, 10/13/99) 
!  (6 ) Now save pure O3 as tracer #25 and NO as tracer #26. 
!        Also declare STT_TEMPO2 as a dynamically allocatable array 
!        in "diag_mod.f". (qli, bmy, 1/5/00)
!  (7 ) Remove NYMD from the argument list, since it is not used anywhere.
!        Added a few cosmetic changes and updated some comments.  Also added
!        code for the GEOS-2 and GEOS-3 data sets.  (bmy, 6/22/00)
!  (8 ) Reference F90 module "bpch2_mod" which contains routines BPCH2_HDR, 
!        BPCH2, and GET_MODELNAME for writing data to binary punch files. 
!        (bmy, 6/22/00)
!  (9 ) For now, skip over unused tracers #42, #43, #44 (bmy, 7/17/00)
!  (10) Remove obsolete code from 5/11/00. (bmy, 8/31/00)
!  (11) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
!  (12) Only save NTRACE tracers instead of NNPAR tracers.  For full
!        chemistry, save O3 as tracer NTRACE+1 and NO as tracer NTRACE+2. 
!        Also added parallel DO-loops. (bmy, 10/13/00) 
!  (13) Eliminated obsolete commented-out code.  Updated comments and
!        made some cosmetic changes. (bmy, 4/20/01)
!  (14) Made ZTAU1 and ZTAU2 into local variables -- removed these from
!        CMN_TIMES to avoid conflicts w/ DIAG50 and DIAG51 (bmy, 7/17/01)
!  (15) Now pass XI0+I0 as IFIRST and XJ0+J0 as JFIRST in call to BPCH2.
!        This will take care of window regions smaller than the globe.
!        (yxw, bmy, 12/18/01)
!  (16) Eliminate obsolete code from 12/01 (bmy, 2/27/02)
!  (17) Now reference IU_ND50 and from "file_mod.f".  Now use IU_ND50 instead 
!        of IUT as the file unit #.  Also call routine OPEN_BPCH2_FOR_WRITE to
!        start writing to the output file.  (bmy, 6/26/02)
!  (18) Now reference AD from "dao_mod.f".  Now reference IDTOX from F90
!        module "tracerid_mod.f". (bmy, 11/6/02)
!  (19) Now remove FIRSTDIAG50 from the arg list -- this is now a local
!        variable.  Now use functions GET_DAY, GET_TAU, GET_TAUE, 
!        GET_LOCALTIME, and TIMESTAMP_STRING from "time_mod.f".  Now use 
!        GET_XOFFSET and GET_YOFFSET from "grid_mod.f". (bmy, 3/14/03)
!  (20) LINUX has a problem putting a function call w/in a WRITE statement.  
!        Now save output from TIMESTAMP_STRING to STAMP and print that.
!        (bmy, 9/29/03)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DAO_MOD,      ONLY : AD
      USE FILE_MOD,     ONLY : IU_ND50
      USE GRID_MOD,     ONLY : GET_XOFFSET, GET_YOFFSET
      USE DIAG_MOD,     ONLY : STT_TEMPO2
      USE TIME_MOD,     ONLY : GET_DAY, GET_TAU, GET_TAUE, 
     &                         GET_LOCALTIME, TIMESTAMP_STRING
      USE TRACERID_MOD, ONLY : IDTOX

      IMPLICIT NONE

#     include "CMN_SIZE"      ! Size parameters
#     include "CMN"           ! STT, NSRCX
#     include "CMN_TIMES"     ! Timeseries parameters
#     include "CMN_O3"        ! FO3, XNUMOL
#     include "CMN_SETUP"     ! HR1_NO, HR2_NO
#     include "CMN_DIAG"      ! TRCOFFSET

      ! Local variables
      CHARACTER(LEN=16), SAVE :: FILENAME
      CHARACTER(LEN=16)       :: STAMP
      LOGICAL, SAVE           :: FIRSTDIAG50 = .TRUE.
      INTEGER, SAVE           :: COUNT_NO(IIPAR,JJPAR)
      INTEGER, SAVE           :: DAY_LAST
      INTEGER, SAVE           :: COUNT_IN_DAY
      INTEGER, SAVE           :: NUMBER_OF_ARCHIVAL
      INTEGER                 :: TT, I_MID, I, J, L, IOS
      INTEGER                 :: XTRAC, XI0, XJ0, XL0, N, XTRAC2
      INTEGER                 :: XI1, XJ1, XL1, IFIRST, JFIRST

      REAL*4                  :: ZARRAY(NI,NJ,NL)

      REAL*8                  :: XLOCTM
      REAL*8,  SAVE           :: ZTAU1, ZTAU2

      ! For binary punch file, version 2.0
      REAL*4                  :: LONRES, LATRES
      INTEGER, PARAMETER      :: HALFPOLAR = 1
      INTEGER, PARAMETER      :: CENTER180 = 1

      CHARACTER(LEN=20)       :: MODELNAME
      CHARACTER(LEN=40)       :: CATEGORY
      CHARACTER(LEN=40)       :: UNIT     = '' 
      CHARACTER(LEN=40)       :: RESERVED = ''
      CHARACTER(LEN=80)       :: TITLE
!
!*****************************************************************************
!  DIAG50 begins here!
!
!  Information for the header of the file
!*****************************************************************************
!
      TITLE     = 'GEOS-CTM time series for geographical domain'
      LONRES    = DISIZE
      LATRES    = DJSIZE

      ! Get the proper model name for the binary punch file (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()
!
!*****************************************************************************
!  If it is the first time in the subroutine, then open the output file
!  for writing.  Also set some variables and allocate the STT_TEMPO2 array. 
!*****************************************************************************
!
      IF ( FIRSTDIAG50 ) THEN

          ! Open file
          FILENAME = 'ts24h.bpch'

          WRITE( 6, 100 ) TRIM( FILENAME )
 100      FORMAT( '     - DIAG50: Opening file ', a )

          ! Open *bpch file for output
          CALL OPEN_BPCH2_FOR_WRITE( IU_ND50, FILENAME, TITLE )

          DAY_LAST            = GET_DAY()
          FIRSTDIAG50         = .FALSE.
          ZTAU1               = GET_TAU()
        
          ! Initialize the STT_TEMPO2 array
          STT_TEMPO2(:,:,:,:) = 0D0

          ! Reset FIRSTDIAG50
          FIRSTDIAG50 = .FALSE.

       ENDIF
!
!*****************************************************************************
!  Test if it is a new day.  If it is a new day :
!  (1) Compute the average of tracers, O3, and NO
!  (2) Write the average values in the newfile.
!  (3) Flush file and reset variables to zero
!  (4) Accumulate in the array for the first time step of the new day
!*****************************************************************************
!
      IF ( GET_DAY() /= DAY_LAST .or. GET_TAU() == GET_TAUE() ) THEN

         !==============================================================
         ! (1) Compute the mean of tracers, O3, NO
         !==============================================================

         ! Now only process NTRACE tracers (bmy, 10/13/00)
         ! Add parallel DO loops (bmy, 10/13/00)
!$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, J, L, N )
         DO N = 1, NTRACE
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Compute the mean of each tracer
            STT_TEMPO2(I,J,L,N) = STT_TEMPO2(I,J,L,N) / COUNT_IN_DAY

            ! Compute the mean of pure O3, and store as tracer NTRACE+1.
            ! Only do this for a full chemistry run (bmy, 10/13/00)
            IF ( NSRCX == 3 .and. N == IDTOX ) THEN
               STT_TEMPO2(I,J,L,NTRACE+1) = 
     &              STT_TEMPO2(I,J,L,NTRACE+1) / COUNT_IN_DAY      
            ENDIF
         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! Compute average of NO -- where COUNT_NO is nonzero.  Also store 
         ! NO as tracer NTRACE+2 -- only for full chemistry  (bmy, 10/13/00)
         IF ( NSRCX == 3 ) THEN

!$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, J, L, N )
            DO L = 1, LLPAR
            DO I = 1, IIPAR 
            DO J = 1, JJPAR
               IF ( COUNT_NO(I,J) .NE. 0 ) THEN
                  STT_TEMPO2(I,J,L,NTRACE+2) = 
     X                 STT_TEMPO2(I,J,L,NTRACE+2) / COUNT_NO(I,J)
               ELSE
                  STT_TEMPO2(I,J,L,NTRACE+2) = 0d0
               ENDIF
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ENDIF
         
         ZTAU2 = GET_TAU()
          
         !==============================================================
         ! (2) ND45 diagnostic -- write to binary punch file
         !     Wrap around the International Date Line if necessary
         !==============================================================
         IF ( DIAG_AREA == 45 ) THEN
            CATEGORY = 'IJ-AVG-$'

            ! Loop over the number of tracers in "timeseries.dat"
            DO TT = 1, NTRAC_AREA

               ! XTRAC = Current tracer number
               XTRAC = TRAC_AREA(TT)
               !WRITE(6,*) 'TRAC_AREA = ', TRAC_AREA(TT)

               ! XTRAC > 90 are for DIAG49 only...skip over (bmy, 1/5/00)
               IF ( XTRAC > 90 ) CYCLE

               ! For now, skip over cloud diagnostics (bmy, 7/17/00)
               IF ( XTRAC == 42 ) CYCLE
               IF ( XTRAC == 43 ) CYCLE
               IF ( XTRAC == 44 ) CYCLE

               ! Get minimum boundaries for I,J,L dimensions 
               XI0 = IMIN_AREA
               XJ0 = JMIN_AREA
               XL0 = LMIN_AREA

               ! Get maximum boundaries for I,J,L dimensions
               XI1 = ( XI0 - 1 ) + NI
               XJ1 = ( XJ0 - 1 ) + NJ
               XL1 = ( XL0 - 1 ) + NL

               ! Get nested-grid offsets
               IFIRST = XI0 + GET_XOFFSET( GLOBAL=.TRUE. )
               JFIRST = XJ0 + GET_YOFFSET( GLOBAL=.TRUE. )

               IF ( WORD_WRAP ) THEN

                  ! Here we are spanning the date line...
                  I_MID = (IIPAR-XI0) + 1
 
                  ! Now avoid subscript errors (bmy, 5/11/00)
                  ZARRAY( 1:I_MID, 1:NJ, 1:NL ) =
     &               STT_TEMPO2( XI0:IIPAR, XJ0:XJ1, XL0:XL1, XTRAC )
                  
                  ZARRAY( I_MID+1:NI, 1:NJ, 1:NL ) =
     &               STT_TEMPO2( 1:IMAX_AREA, XJ0:XJ1, XL0:XL1, XTRAC )

               ELSE

                  ! Here we are NOT spanning the date line...
                  ! Now avoid subscript errors (bmy, 5/11/00)
                  ZARRAY( 1:NI, 1:NJ, 1:NL ) = 
     &               STT_TEMPO2( XI0:XI1, XJ0:XJ1, XL0:XL1, XTRAC )

               ENDIF

               ! Write data block to binary punch file v. 2.0
               CALL BPCH2( IU_ND50,    MODELNAME,       LONRES,   
     &                     LATRES,    HALFPOLAR,       CENTER180, 
     &                     CATEGORY,  XTRAC+TRCOFFSET, UNIT,      
     &                     ZTAU1,     ZTAU2,           RESERVED,  
     &                     NI,        NJ,              NL,   
     &                     IFIRST,    JFIRST,          XL0, 
     &                     ZARRAY )

            ENDDO

         ENDIF 

         !==============================================================
         ! (3) Flush the file and zero variables after writing
         !==============================================================
         IF ( GET_TAU() == GET_TAUE() ) CLOSE ( IU_ND50 )

         DAY_LAST            = GET_DAY()
         COUNT_IN_DAY        = 0
         COUNT_NO(:,:)       = 0
         STT_TEMPO2(:,:,:,:) = 0d0

         !==============================================================
         ! (4) Accumulate for the first time step of this new day...
         !==============================================================
         IF ( DIAG_AREA == 45 ) THEN

            ! Increment counter and reset ZTAU1, if necessary
            COUNT_IN_DAY = COUNT_IN_DAY + 1
            IF ( COUNT_IN_DAY .EQ. 1 ) ZTAU1 = GET_TAU()

            ! Echo output
            STAMP = TIMESTAMP_STRING()
            WRITE( 6, 110 ) STAMP
 110        FORMAT( '     - DIAG50: Accumulation at ', a )

            ! Accumulate TRACERS into STT_TEMPO as tracer #'s 1 - NTRACE 
            ! Add parallel DO-loops (bmy, 10/13/00)
!$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, J, L, N )
            DO N = 1, NTRACE
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               STT_TEMPO2(I,J,L,N) = STT_TEMPO2(I,J,L,N) +
     &              STT(I,J,L,N) * TCVV(N) / AD(I,J,L)
            ENDDO
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            ! Only process O3 and NO for full chemistry (bmy, 10/13/00)
            IF ( NSRCX == 3 ) THEN

               ! Accumulate O3 into STT_TEMPO as tracer # NTRACE+1 
               N = NTRACE + 1

!$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, J, L )
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  CALL O3COMP(I,J,L)
                  
                  STT_TEMPO2(I,J,L,N) = STT_TEMPO2(I,J,L,N) +
     &                                  STT(I,J,L,IDTOX) * TCVV(IDTOX) *
     &                                  FO3(5)           / AD(I,J,L)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

               ! Accumulate NO into STT_TEMPO as tracer # NNPAR+2
               N = NTRACE + 2
               DO J = 1, JJPAR
               DO I = 1, IIPAR

                  ! Get local time
                  XLOCTM = GET_LOCALTIME( I ) 

                  ! Only save NO where the local time is
                  ! between HR1_NO and HR2_NO
                  IF ( XLOCTM >= HR1_NO .and. XLOCTM <= HR2_NO ) THEN
                     STT_TEMPO2(I,J,:,N) = 
     &                    STT_TEMPO2(I,J,:,N) + SAVENO(I,J,:)

                     COUNT_NO(I,J) = COUNT_NO(I,J) + 1
                  ENDIF
               ENDDO
               ENDDO
            ENDIF               ! NSRCX == 3
         ENDIF                  ! nd45
!
!*****************************************************************************
!  If it is not a new day.. 
!
!  Accumulate in the zarray.
!  For now, we have only the ND45 diag -
!*****************************************************************************
!
      ELSE
         IF ( DIAG_AREA == 45 ) THEN
            COUNT_IN_DAY = COUNT_IN_DAY + 1

            ! Echo output
            STAMP = TIMESTAMP_STRING()
            WRITE( 6, 110 ) STAMP

            ! Accumulate TRACERS into STT_TEMPO as tracer #'s 1 - NTRACE
            ! Add parallel DO-loops (bmy, 10/13/00)
!$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, J, L, N )
            DO N = 1, NTRACE
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               STT_TEMPO2(I,J,L,N) = STT_TEMPO2(I,J,L,N) +
     &              STT(I,J,L,N) * TCVV(N) / AD(I,J,L)
            ENDDO
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            ! Only process O3 and NO for full chemistry (bmy, 10/13/00)
            IF ( NSRCX == 3 ) THEN

               ! Accumulate O3 into STT_TEMPO as tracer # NNPAR+1 
               ! Add parallel DO-loops (bmy, 10/13/00)
               N = NTRACE + 1

!$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, J, L )
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  CALL O3COMP(I,J,L)
               
                  STT_TEMPO2(I,J,L,N) = STT_TEMPO2(I,J,L,N) +
     &                                  STT(I,J,L,IDTOX) * TCVV(IDTOX) *
     &                                  FO3(5)           / AD(I,J,L)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

               ! Accumulate NO into STT_TEMPO as tracer # NNPAR+2
               N = NTRACE + 2
               DO J = 1, JJPAR
               DO I = 1, IIPAR

                  ! Get local time
                  XLOCTM = GET_LOCALTIME( I )

                  ! Only save NO where the local time is
                  ! between HR1_NO and HR2_NO
                  IF ( XLOCTM >= HR1_NO .and. XLOCTM <= HR2_NO ) THEN
                     STT_TEMPO2(I,J,:,N) =
     &                    STT_TEMPO2(I,J,:,N) + SAVENO(I,J,:)
                     
                     COUNT_NO(I,J) = COUNT_NO(I,J) + 1
                  ENDIF
               ENDDO
               ENDDO
            ENDIF               ! NSRCX == 3
         ENDIF                  ! nd45

      ENDIF                     ! accumulation or new day...
!
!*****************************************************************************
!  Return to calling program
!*****************************************************************************
!
      END SUBROUTINE DIAG50
