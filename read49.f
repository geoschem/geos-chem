! $Id: read49.f,v 1.2 2004/07/15 18:17:46 bmy Exp $
      SUBROUTINE READ49
!
!******************************************************************************
!  Subroutine READ49 reads the "timeseries.dat" file and returns values 
!  pertaining to timeseries geographical domain & diagnostics in the 
!  CMN_TIMES common block. (bey, bmy, 5/28/99, 7/15/04)
!
!  Sample "timeseries.dat" input file:
!  ===========================================================================
!  # Informations about geopraphical domain for timeseries diagnostic.
!  # Turn on the diag 49 in input.ctm.
!  # imin,imax          = indices of the limit boxes of the domain in long.
!  #                      (note : imin can be greater than imax if
!  #                      we are crossing the date line)
!  # jmin,jmax          = indices of the limit boxes of the domain in latitude.
!  # layermin, layermax = range of vertical layers.
!  # date1,date2        = first date and ending date for the archival.
!  #                      they are given as YYMMDD (idem as input.geos).
!  #                      the archival includes date1 and date2.
!  # frequency          = frequency for the archival.
!  # type of diag       = type of diag (see diag1.f for description).
!  #                      for now, only ND45 works.
!  # tracers            = list of tracers corresponding to the diagnostic.
!  #
!  # Writing format : start with a "*" for each area
!  #                  put a "=" after each line description
!  #                  put a "," between two data
!  #
!  *area   #1
!   imin,imax          =  57,13
!   jmin,jmax          =  18,38
!   layermin, layermax =  1,17
!   date1,date2        =  940901,940902
!   frequency          =  60
!   type of diag       =  45
!   tracers            =  1,2,3,4
!
!  NOTES:
!  (1 ) Despite the name, READ49 is actually used for ND49, ND50, ND51
!        timeseries diagnostics (bmy, 12/18/01)
!  (2 ) Now use subroutine IOERROR to trap I/O errors (bmy, 5/28/99)
!  (3 ) Now trap I/O errors (bmy, 4/9/99)
!  (4 ) Allow for more tracers to be included in input file (amf, 7/16/99)
!  (5 ) Modified error check for LMAX_AREA, LMIN_AREA.  
!  (6 ) Also updated comments. (bmy, 5/11/00)!
!  (7 ) Reference "time_mod.f" which contains routines TIMECHECK and CTM_TIME 
!        (bmy, 6/22/00)
!  (8 ) Now use IOS /= 0 to trap both I/O errors and EOF (bmy, 9/13/00)
!  (9 ) Added error check for IMAX, JMAX & changed HEADER to CHARACTER*100 
!        (bmy, 10/11/00)
!  (10) Removed obsolete code from 10/00 (bmy, 12/21/00)
!  (11) Remove double declaration of  LYEAR, LMONTH, LDAY, LLHR, LLMI
!        Also renamed NREAD to N_READ to avoid conflict w/ "CMN" (bmy, 7/16/01)
!  (12) Removed obsolete code (bmy, 9/4/01)
!  (13) Removed duplicate definitions for LDAY, NTOTDAY (bmy, 11/15/01)
!  (14) Now error check to make sure that IMAX_AREA, IMIN_AREA, JMAX_AREA, and
!        JMIN_AREA are set properly for windows less than global size.  Also
!        updated comments and made cosmetic changes. (yxw, bmy, 12/18/01)
!  (15) Now reference IOERROR from "file_mod.f" (bmy, 6/27/02)
!  (16) Now reference ERROR_STOP from "error_mod.f".  Updated comments and
!        made cosmetic changes. (bmy, 10/15/02)
!  (17) Now use functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f".
!        Now make I0, J0 local variables.  Removed references to TIMECHECK 
!        and CTM_TIME from the old "time_mod.f".  Now uses function GET_TAU0 
!        from "bpch2_mod.f".  Now use functions YMD_EXTRACT from "time_mod.f"
!        (bmy, 2/11/03)
!  (18) Now can read more than a 100 char string for tracer numbers. 
!        (rjp, bmy, 7/15/04)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD, ONLY : GET_TAU0
      USE ERROR_MOD, ONLY : ERROR_STOP 
      USE FILE_MOD,  ONLY : IOERROR
      USE GRID_MOD,  ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! NSRCX, etc
#     include "CMN_TIMES" ! Timeseries variables

      ! Local variables
      INTEGER            :: INFO(50), K, KK, M, N_READ, NBLANK, NINFO
      INTEGER            :: DAT1, DAT2, LNT, LNHMS, TOTO1, TOTO2, IOS
      INTEGER            :: LYEAR, LMONTH, LDAY, LLHR, LLMI, LLSE
      INTEGER            :: LNTOTDAY, LJDAY, LLMONTH, LJDATE, LJYEAR
      INTEGER            :: LNTAU,  LITAU,  LIDAY,  LNSEASON
      REAL*8             :: LTAU, LTOFDAY
      CHARACTER(LEN=1  ) :: COMMA = ','
      CHARACTER(LEN=1  ) :: BLANK = ' '
      CHARACTER(LEN=1  ) :: EQUAL = '='
      CHARACTER(LEN=4  ) :: LJMONTH
      CHARACTER(LEN=5  ) :: TNAME
      CHARACTER(LEN=20 ) :: NOTION
      CHARACTER(LEN=255) :: HEADER

      ! Now make I0, J0 local variables (bmy, 2/11/03)
      INTEGER            :: I0, J0

      !=================================================================
      ! READ49 begins here!
      !=================================================================

      ! Get nested-grid offsets (bmy, 2/11/03)
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()

      !=================================================================
      ! Read IMIN_AREA, IMAX_AREA, JMIN_AREA, JMAX_AREA,
      !      LMIN_AREA, LMAX_AREA, DAT1_AREA, DAT2_AREA,
      !      FREQ_AREA, TRAC_AREA
      !
      ! NREAD =  10 informations + number of tracers
      !=================================================================

      ! Echo info
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, '(a)' ) 'READ49 -- reading "timeseries.dat" file...' 

      ! Number of information fields (IMAX_AREA, IMIN_AREA, etc)
      NINFO = 10                 
      
      ! Open temporary file
      OPEN(67,FILE='temp.dat',STATUS='unknown', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 67, 'read49:1' )

      ! Open "timeseries dat" file
      OPEN(66,FILE='timeseries.dat',STATUS='old',FORM='formatted',
     &     IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 66, 'read49:2' )

      ! read fields from "timeseries.dat", use temp file for storage
      WRITE( 6, '(a80)' )
 99   READ( 66, '(a5)', IOSTAT=IOS ) TNAME 
      IF ( IOS < 0 ) GOTO 999
      IF ( IOS > 0 ) CALL IOERROR( IOS, 66, 'read49:3' )

      IF (TNAME .eq. '*area') THEN
         N_READ=1
 22      READ(66,'(a)', IOSTAT=IOS) HEADER
         IF ( IOS /= 0 ) CALL IOERROR( IOS, 66, 'read49:4' )

         M = 0
          DO KK=1,LEN_TRIM( HEADER ) + 1
            IF (HEADER(KK:KK).EQ.EQUAL) THEN 
              DO K=KK+1,LEN(HEADER) + 1
               IF(HEADER(K:K).EQ.COMMA.OR.K.EQ.LEN_TRIM(HEADER)+1) THEN
                     REWIND(67)
                     WRITE(67,*, IOSTAT=IOS )NOTION(1:M)
                     IF ( IOS /= 0 ) CALL IOERROR( IOS, 67, 'read49:5' )

                     REWIND(67, IOSTAT=IOS )
                     IF ( IOS /= 0 ) CALL IOERROR( IOS, 66, 'read49:6' )

                     READ(67,*, IOSTAT=IOS ) INFO(N_READ)
                     IF ( IOS /= 0 ) CALL IOERROR( IOS, 66, 'read49:7' )

                     CLOSE(67)
                     N_READ = N_READ+1
                     M = 0
                     IF(K == LEN_TRIM(HEADER)+1) GOTO 22
               ENDIF    

               IF(HEADER(K:K).NE.COMMA) THEN
                     IF(HEADER(K:K).EQ.BLANK)NBLANK=NBLANK+1
                     M = M+1
                     NOTION(M:M)=HEADER(K:K)
C                    PRINT*, NOTION(M:M)
               ENDIF    
               ENDDO  
            ENDIF    
         ENDDO           
      ENDIF  
      
      GOTO 99
 999  CONTINUE
      
      ! Copy stuff from the INFO array to variables
      IMIN_AREA  = INFO(1) 
      IMAX_AREA  = INFO(2) 
      JMIN_AREA  = INFO(3) 
      JMAX_AREA  = INFO(4) 
      LMIN_AREA  = INFO(5)
      LMAX_AREA  = INFO(6)
      DAT1       = INFO(7) 
      DAT2       = INFO(8)
      FREQ_AREA  = INFO(9)
      DIAG_AREA  = INFO(10)
      TRAC_AREA  = INFO(NINFO+1:NTRAC_AREA)

      ! Convert DAT1, DAT2 from YYMMDD to YYYYMMDD if necessary
      DAT1 = NYMD6_2_NYMD8( DAT1 )
      DAT2 = NYMD6_2_NYMD8( DAT2 )
      
      ! Compute number of tracers
      N_READ     = N_READ - 1 
      NTRAC_AREA = N_READ - NINFO

      !=================================================================
      ! Error check boundaries of geographical domain and 
      ! start and end times for the timeseries diagnostic output
      !=================================================================

      ! Error check IMIN_AREA
      IF ( IMIN_AREA+I0 < 1 .or. IMIN_AREA+I0 > IGLOB ) THEN
         CALL ERROR_STOP( 'IMIN+I0 < 1 OR IMIN+I0 > IGLOB!', 'read49.f')
      ENDIF

      ! Error check IMAX_AREA
      IF ( IMAX_AREA+I0 < 1 .or. IMAX_AREA+I0 > IGLOB ) THEN
         CALL ERROR_STOP( 'IMAX+I0 < 1 OR IMAX+I0 > IGLOB!', 'read49.f')
      ENDIF

      ! Compute longitude limits to write to disk (bey, bmy, 3/16/99)
      ! Also, if IMIN_AREA > IMAX_AREA, then we are wrapping around 
      ! the International Date Line (bmy, 3/16/99, 12/18/01)
      IF ( IMAX_AREA > IMIN_AREA ) THEN 
         NI = ( IMAX_AREA - IMIN_AREA ) + 1
      ELSE 
         NI = ( IIPAR - IMIN_AREA ) + 1 + IMAX_AREA
         WORD_WRAP = .TRUE.
         WRITE( 6, '(a)' ) 'We are wrapping over the date line!'
      ENDIF
      
      ! Error check JMIN_AREA
      IF ( JMIN_AREA+J0 < 1 .or. JMIN_AREA+J0 > JGLOB ) THEN
         CALL ERROR_STOP( 'JMIN+J0 < 1 OR JMIN+J0 > JGLOB!', 'read49.f')
      ENDIF
     
      ! Error check JMAX_AREA
      IF ( JMAX_AREA+J0 < 1 .or. JMAX_AREA+J0 > JGLOB ) THEN
         CALL ERROR_STOP( 'JMAX+J0 < 1 OR JMAX+J0 > JGLOB!', 'read49.f')
      ENDIF

      ! Compute latitude limits to write to disk (bey, bmy, 3/16/99)
      IF ( JMAX_AREA > JMIN_AREA ) THEN
         NJ = ( JMAX_AREA - JMIN_AREA ) + 1
      ELSE
         CALL ERROR_STOP( 'JMAX < JMIN!', 'read49.f' )
      ENDIF     
  
      ! Sigma level values out of bounds! (bey, bmy, 3/16/99)
      IF ( LMIN_AREA < 1 .or. LMAX_AREA > LLPAR ) THEN 
         CALL ERROR_STOP( 'LMIN < 1 OR LMIN > LLPAR!', 'read49.f' )
      ENDIF

      ! Compute sigma level limits to write to disk (bey, bmy, 3/16/99)
      ! Change criteria from LMAX_AREA > LMIN_AREA to  LMAX_AREA >= LMIN_AREA
      ! this allows us to only save one level (bmy, 5/11/00)
      IF ( LMAX_AREA  >= LMIN_AREA ) THEN  
         NL = ( LMAX_AREA - LMIN_AREA ) + 1
      ELSE
         CALL ERROR_STOP( 'LMAX < LMIN!', 'read49.f' )
      ENDIF

      !=================================================================
      ! Get the TAU values (DAT1_AREA, DAT2_AREA) corresponding to the 
      ! beginning and end of the timeseries diagnostic interval
      !=================================================================
      LNHMS = 0
      LNT   = 0

      !=================================================================
      ! Find TAU values for beginning, ending dates (bmy, 2/11/03)
      !=================================================================

      ! Split DAT1 (starting date) up into year, month, day
      CALL YMD_EXTRACT( DAT1, LYEAR, LMONTH, LDAY )

      ! Find TAU value corresponding to DAT1
      DAT1_AREA = GET_TAU0( LMONTH, LDAY, LYEAR ) 

      ! Split DAT2 (ending date) up into year, month, day
      CALL YMD_EXTRACT( DAT2, LYEAR, LMONTH, LDAY )

      ! Find TAU value corresponding to DAT2
      DAT2_AREA = GET_TAU0( LMONTH, LDAY, LYEAR )

      ! Error check DAT1_AREA, DAT2_AREA
      IF ( DAT1_AREA > DAT2_AREA ) THEN
         CALL ERROR_STOP( 'DATE2 < DATE1!', 'read49.f' )
      ENDIF

      !=================================================================
      ! Echo information to the screen
      !=================================================================
      WRITE(6,'(''# of fields read : '',i4)'  ) N_READ+1
      WRITE(6,'(''# of info fields : '',i4)'  ) NINFO
      WRITE(6,'(''Tracer numbers   : '',15i4)') TRAC_AREA(1:NTRAC_AREA)
      WRITE(6,'(''IMIN, IMAX       : '',2i4)' ) IMIN_AREA, IMAX_AREA
      WRITE(6,'(''JMIN, JMAX       : '',2i4)' ) JMIN_AREA, JMAX_AREA
      WRITE(6,'(''NI, NJ, NL       : '',3i4)' ) NI, NJ, NL
      WRITE(6,'(''DATE1, DATE2     : '',2f11.2)') DAT1_AREA, DAT2_AREA

      ! Pretty output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      
      ! Return to calling program
      END SUBROUTINE READ49
