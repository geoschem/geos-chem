! $Id: read_tindex.f,v 1.1 2003/06/30 20:26:01 bmy Exp $
      SUBROUTINE READ_TINDEX
!
!******************************************************************************
!  Subroutine READ_TINDEX reads the 'diag.dat' file, which contains
!  the tracer numbers that will be printed out for each diagnostic.
!  (bmy, 7/8/98, 6/27/02)
!
!  NOTES:
!  (1 ) READ_TINDEX is written using fixed-form F90 syntax.
!  (2 ) Uses the following subroutines from CHARPAK: TXTEXT, TXT2INUM, CSTRIP
!  (3 ) Now trap I/O errors with IOERROR.F (bmy, 12/9/99)
!  (4 ) Removed reference to "biomass.h" and "comtrid.h" (bmy, 9/11/00)
!  (5 ) Now reference NBFTRACE from "biofuel_mod.f" and NBIOTRCE from
!        "biomass_mod.f".  Also readjust format to read in only N and
!        LINE2 from the "diag.dat" file. (bmy, 4/17/01)
!  (6 ) Remove obsolete commented-out code (bmy, 4/23/01)
!  (7 ) Now reference routines CSTRIP, TXTEXT, TXT2INUM from "charpak_mod.f". 
!        Also updated comments, made cosmetic changes. (bmy, 10/15/01)
!  (8 ) Now reference IU_FILE and IOERROR from "file_mod.f".  Now use IU_FILE
!        instead of ILUN as the file unit number. (bmy, 6/27/02)
!******************************************************************************
!
      ! References to F90 modules
      USE BIOFUEL_MOD, ONLY : NBFTRACE
      USE BIOMASS_MOD, ONLY : NBIOTRCE
      USE CHARPAK_MOD, ONLY : CSTRIP, TXT2INUM, TXTEXT
      USE FILE_MOD,    ONLY : IU_FILE, IOERROR

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_DIAG"  ! MAXTRACER

      INTEGER, PARAMETER  :: ILARGE = 999999
      INTEGER             :: I, J, N, IOS, IMIN, IMAX, IDUM
      INTEGER             :: TDUMMY(MAXTRACER)
      INTEGER             :: IFLAG1, IFLAG2, ICOUNT
      INTEGER             :: COL1,   COL2,   COL3
      CHARACTER (LEN=255) :: LINE, LINE2, STRING 
      CHARACTER (LEN=255) :: WORD, SUBWORD, TMP1, TMP2

      !=================================================================
      ! READ_TINDEX starts here!
      !
      ! Initialize arrays using F90 notation
      !=================================================================
      TCOUNT  = 0
      TINDEX  = 0
      TMAX    = 0

      !=================================================================
      ! Open file and two (2) header lines
      !=================================================================
      OPEN( IU_FILE, FILE='diag.dat', STATUS='old', IOSTAT=IOS )
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_tindex:1' )

      !=================================================================
      ! Read data from file
      !
      ! TOFFSET = offset added to tracer number in the punch file 
      !           (for IDL prog's)
      !
      ! TINDEX  = index array of tracers to be printed out for 
      !           each NDxx diagnostic
      ! 
      ! TCOUNT  = number of tracers to be printed out for each 
      !           NDxx diagnostic
      ! 
      ! The data will also be sorted((via HPSORT) so that tracer 
      ! numbers will be in ascending order.  Also, so that all valid 
      ! (e.g. non-zero) tracer numbers will be contiguous entries 
      ! in the TINDEX array.
      ! 
      ! NOTE: IOS < 0 means end of file
      !       IOS > 0 means there is an I/O error -- trap it w/ IOERROR
      !=================================================================
      DO
         ! Read a line from the file
         READ ( IU_FILE, '(a)', IOSTAT=IOS ) LINE
         
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_tindex:2' )

         ! Go to next iteration if this is a header line (marked with #)
         IF ( LINE( 1:1 ) == '#' ) CYCLE

         ! Read list of tracer numbers from LINE and store in LINE2
         ! (IDUM is TOFFSET(N), which we don't use
         READ( LINE, '(i4,8x,a200)', IOSTAT=IOS ) N, LINE2
         IF ( IOS > 0 ) CALL IOERROR( IOS, -1, 'read_tindex:3' )
         
         ! Initialize some variables for each new line
         ! Set TDUMMY = ILARGE everywhere for sorting purposes (see below)
         TDUMMY = ILARGE
         ICOUNT = 0
         COL1   = 1
         COL2   = 1
         IFLAG1 = 0
         IFLAG2 = 0

         ! Compress LINE2 so there isn't any white space
         CALL CSTRIP ( LINE2 )

         ! Parse the tracer numbers, which are separated by commas
         DO WHILE ( IFLAG1 == 0 )

            ! Look for strings beteeen commas
            CALL TXTEXT ( ',', LINE2, COL1, WORD, IFLAG1 )

            ! Analyze each substring for dashes
            CALL TXTEXT ( '-', WORD, COL2, SUBWORD, IFLAG2 )
           
            ! We have found a dash!  Get the numbers on both sides of the 
            ! dash, since these the min and max of the tracer range
            IF ( IFLAG2 == 0 ) THEN
               TMP1 = TRIM( WORD(      1:COL2-1      ) )
               TMP2 = TRIM( WORD( COL2+1:LEN_TRIM( WORD ) ) )
            
            ! If we haven't found a dash, then there is only one number,
            ! so that number is both the min and max of the tracer range
            ELSE IF ( IFLAG2 == -1 ) THEN
               TMP1 = TRIM( WORD )
               TMP2 = TRIM( WORD )
            ENDIF
            
            ! Convert the text to integer numbers, flag I/O errors
            CALL TXT2INUM ( '(i255)', TRIM( TMP1 ), IMIN, IOS )
            IF ( IOS > 0 ) CALL IOERROR( IOS, -1, 'read_tindex:4' )

            CALL TXT2INUM ( '(i255)', TRIM( TMP2 ), IMAX, IOS )  
            IF ( IOS > 0 ) CALL IOERROR( IOS, -1, 'read_tindex:5' )

            ! Put the tracer numbers into the TINDEX array!
            J = 0
            DO I = MIN( IMIN, IMAX ), MAX( IMIN, IMAX )
               J = J + 1
               TDUMMY( ICOUNT + J ) = I
            ENDDO
            ICOUNT = ICOUNT + J 
         ENDDO

         ! Call HPSORT to perform a heap sort on the TDUMMY array
         ! Reset the large values to 0 afterwards
         CALL HPSORT_INT( MAXTRACER, TDUMMY )

         WHERE ( TDUMMY == ILARGE ) TDUMMY = 0
        
         ! Assign to the TCOUNT and TINDEX arrays!
         TCOUNT(N)   = COUNT( TDUMMY > 0 ) 
         TINDEX(N,:) = TDUMMY(:)
      ENDDO

      !=================================================================
      ! First check that the values of TINDEX are not out of range
      !
      ! Compute TMAX, the maximum number of tracers (not levels!!) 
      ! per diagnostic.  TMAX is used as a check to avoid subscript 
      ! range violations. Use F90 array masks.  TMAX has been 
      ! initialized to zero everywhere.
      !
      ! Use NBIOTRCE and NBFTRACE instead of PD28 and PD34.  This 
      ! ensures that the correct number of biomass and biofuel tracers 
      ! will be printed out. (bmy, 4/17/01)
      !=================================================================
 999  CONTINUE
      TMAX(  1 ) = MIN( TCOUNT(  1 ), PD01     )
      TMAX(  2 ) = MIN( TCOUNT(  2 ), PD02     )
      TMAX(  3 ) = MIN( TCOUNT(  3 ), PD03     )
      TMAX(  4 ) = MIN( TCOUNT(  4 ), PD04     )
      TMAX(  5 ) = MIN( TCOUNT(  5 ), PD05     )
      TMAX(  6 ) = MIN( TCOUNT(  6 ), PD06     )
      TMAX(  7 ) = MIN( TCOUNT(  7 ), PD07     )
      TMAX(  8 ) = MIN( TCOUNT(  8 ), PD08     )
      TMAX(  9 ) = MIN( TCOUNT(  9 ), PD09     ) 
      TMAX( 10 ) = MIN( TCOUNT( 10 ), PD10     )
      TMAX( 11 ) = MIN( TCOUNT( 11 ), PD11     )
      TMAX( 12 ) = MIN( TCOUNT( 12 ), PD12     )
      TMAX( 13 ) = MIN( TCOUNT( 13 ), PD13     )
      TMAX( 14 ) = MIN( TCOUNT( 14 ), PD14     )
      TMAX( 15 ) = MIN( TCOUNT( 15 ), PD15     )
      TMAX( 16 ) = MIN( TCOUNT( 16 ), PD16     )
      TMAX( 17 ) = MIN( TCOUNT( 17 ), PD17     )
      TMAX( 18 ) = MIN( TCOUNT( 18 ), PD18     )
      TMAX( 19 ) = MIN( TCOUNT( 19 ), PD19     )
      TMAX( 20 ) = MIN( TCOUNT( 20 ), PD20     )
      TMAX( 21 ) = MIN( TCOUNT( 21 ), PD21     )
      TMAX( 22 ) = MIN( TCOUNT( 22 ), PD22     )
      TMAX( 23 ) = MIN( TCOUNT( 23 ), PD23     )
      TMAX( 24 ) = MIN( TCOUNT( 24 ), PD24     )
      TMAX( 25 ) = MIN( TCOUNT( 25 ), PD25     )
      TMAX( 26 ) = MIN( TCOUNT( 26 ), PD26     )
      TMAX( 27 ) = MIN( TCOUNT( 27 ), PD27     )
      TMAX( 28 ) = MIN( TCOUNT( 28 ), NBIOTRCE )
      TMAX( 29 ) = MIN( TCOUNT( 29 ), PD29     )
      TMAX( 30 ) = MIN( TCOUNT( 30 ), PD30     )
      TMAX( 31 ) = MIN( TCOUNT( 31 ), PD31     )
      TMAX( 32 ) = MIN( TCOUNT( 32 ), PD32     )
      TMAX( 33 ) = MIN( TCOUNT( 33 ), PD33     )
      TMAX( 34 ) = MIN( TCOUNT( 34 ), NBFTRACE )
      TMAX( 35 ) = MIN( TCOUNT( 35 ), PD35     )   
      TMAX( 36 ) = MIN( TCOUNT( 36 ), PD36     )
      TMAX( 37 ) = MIN( TCOUNT( 37 ), PD37     )
      TMAX( 38 ) = MIN( TCOUNT( 38 ), PD38     )
      TMAX( 39 ) = MIN( TCOUNT( 39 ), PD39     )
      TMAX( 40 ) = MIN( TCOUNT( 40 ), PD40     )
      TMAX( 41 ) = MIN( TCOUNT( 41 ), PD41     )
      TMAX( 42 ) = MIN( TCOUNT( 42 ), PD42     )
      TMAX( 43 ) = MIN( TCOUNT( 43 ), PD43     )
      TMAX( 44 ) = MIN( TCOUNT( 44 ), PD44     )
      TMAX( 45 ) = MIN( TCOUNT( 45 ), PD45     )
      TMAX( 46 ) = MIN( TCOUNT( 46 ), PD46     )
      TMAX( 47 ) = MIN( TCOUNT( 47 ), PD47     )
      TMAX( 48 ) = MIN( TCOUNT( 48 ), PD48     )
      TMAX( 49 ) = MIN( TCOUNT( 49 ), PD49     )
      TMAX( 50 ) = MIN( TCOUNT( 50 ), PD50     )
      TMAX( 51 ) = MIN( TCOUNT( 51 ), PD51     )
      TMAX( 52 ) = MIN( TCOUNT( 52 ), PD52     )
      TMAX( 53 ) = MIN( TCOUNT( 53 ), PD53     )
      TMAX( 54 ) = MIN( TCOUNT( 54 ), PD54     ) 
      TMAX( 55 ) = MIN( TCOUNT( 55 ), PD55     )
      TMAX( 56 ) = MIN( TCOUNT( 56 ), PD56     )
      TMAX( 57 ) = MIN( TCOUNT( 57 ), PD57     )
      TMAX( 58 ) = MIN( TCOUNT( 58 ), PD58     )
      TMAX( 59 ) = MIN( TCOUNT( 59 ), PD59     )
      TMAX( 60 ) = MIN( TCOUNT( 60 ), PD60     )
      TMAX( 61 ) = MIN( TCOUNT( 61 ), PD61     )
      TMAX( 62 ) = MIN( TCOUNT( 62 ), PD62     )
      TMAX( 63 ) = MIN( TCOUNT( 63 ), PD63     )
      TMAX( 64 ) = MIN( TCOUNT( 64 ), PD64     )
      TMAX( 65 ) = MIN( TCOUNT( 65 ), PD65     )
      TMAX( 66 ) = MIN( TCOUNT( 66 ), PD66     ) 
      TMAX( 67 ) = MIN( TCOUNT( 67 ), PD67     ) 
      TMAX( 68 ) = MIN( TCOUNT( 68 ), PD68     ) 
      TMAX( 69 ) = MIN( TCOUNT( 69 ), PD69     ) 
      TMAX( 70 ) = MIN( TCOUNT( 70 ), PD70     ) 

      !=================================================================
      ! Echo back diagnostic numbers to the screen
      !=================================================================
      WRITE ( 6, '(a)' ) '============================================='
      WRITE ( 6, '(a)' ) 'NDxx -- tracers to be printed'

      DO N = 1, MAXDIAG
         IF ( TCOUNT(N) == 0 .or. TMAX(N) == 0 ) CYCLE
         
         WRITE( 6,'( ''ND'',i2,'': '',40i3)' ) 
     &      N, ( TINDEX(N,J), J=1,TMAX(N) )
      ENDDO
             
      ! Close file
      CLOSE( IU_FILE )

      ! Return to calling program
      END SUBROUTINE READ_TINDEX

