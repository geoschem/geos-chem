C $Id: read_TOMSO3.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
      SUBROUTINE read_TOMSO3(NSKIPTOMS)
C**********************************************************************
      ! References to F90 modules
      USE DAO_MOD, ONLY : AIRVOL
      USE BPCH2_MOD

      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "CMN"
#     include "CMN_SETUP"
#     include "CMN_CO"
#     include "CMN_OH"

      INTEGER I,J,K,L,M,N,NTOMSyear,themonth
      INTEGER NTOMSstart,NTOMSend,NSKIPTOMS,TRACER_IN
      PARAMETER(NTOMSstart=1988,NTOMSend=1999)

      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=40)  :: CATEGORY_IN
      REAL*8             :: TAU0_IN
      REAL*4             :: ARRAY(IIPAR,JJPAR,1)

cbnd      REAL*8, EXTERNAL  :: GET_TAU0

C**********************************************************************
C This SR readds in TOMS total O3 column data from 1988 to 1999.
C The data is on a 4x5 grid. 
C**********************************************************************

       print*,'Read in TOMS O3 column data in SR read_TOMSAI.'
       print*,'The year is ',JYEAR,'.' 

       NTOMSyear=JYEAR-1900
       themonth=MONTH
       CATEGORY_IN = 'O3COLMAP'
c       TRACER_IN   = 2501
       TRACER_IN   = 1
       ctm4x5month(:,:)=1.

C If NSKIPTOMS =1 then use O3 column climatology instead of TOMS data.
C Skip years of no TOMS data.
         NSKIPTOMS=0
       IF(NTOMSyear.EQ.94.OR.NTOMSyear.EQ.95) NSKIPTOMS=1
       IF(NTOMSyear.EQ.93.AND.themonth.GT.4) NSKIPTOMS=1
       IF(NTOMSyear.EQ.96.AND.themonth.LT.8) NSKIPTOMS=1
C acf processed TOMS data from Jan. 1988 to Dec. 1999.
       IF(NTOMSyear.GT.99) NSKIPTOMS=1
       IF(NTOMSyear.LT.88) NSKIPTOMS=1

       IF(NSKIPTOMS.EQ.1) THEN
       print*,'There is no TOMS available for ',JYEAR,'.' 
       ENDIF

       IF(NSKIPTOMS.EQ.0) THEN
C Read in O3 columns by for specific year.
       IF(NTOMSyear.EQ.88) THEN
        FILENAME= TRIM( DATA_DIR ) // 'o3col.TOMS.1988.geos1.4x5'
       ENDIF
       IF(NTOMSyear.EQ.89) THEN
        FILENAME= TRIM( DATA_DIR ) // 'o3col.TOMS.1989.geos1.4x5'
       ENDIF
       IF(NTOMSyear.EQ.90) THEN
        FILENAME= TRIM( DATA_DIR ) // 'o3col.TOMS.1990.geos1.4x5'
       ENDIF
       IF(NTOMSyear.EQ.91) THEN
        FILENAME= TRIM( DATA_DIR ) // 'o3col.TOMS.1991.geos1.4x5'
       ENDIF
       IF(NTOMSyear.EQ.92) THEN
        FILENAME= TRIM( DATA_DIR ) // 'o3col.TOMS.1992.geos1.4x5'
       ENDIF
       IF(NTOMSyear.EQ.93) THEN
        FILENAME= TRIM( DATA_DIR ) // 'o3col.TOMS.1993.geos1.4x5'
       ENDIF
       IF(NTOMSyear.EQ.96) THEN
        FILENAME= TRIM( DATA_DIR ) // 'o3col.TOMS.1996.geos1.4x5'
       ENDIF
       IF(NTOMSyear.EQ.97) THEN
        FILENAME= TRIM( DATA_DIR ) // 'o3col.TOMS.1997.geos1.4x5'
       ENDIF
       IF(NTOMSyear.EQ.98) THEN
        FILENAME= TRIM( DATA_DIR ) // 'o3col.TOMS.1998.geos1.4x5'
       ENDIF
       IF(NTOMSyear.EQ.99) THEN
        FILENAME= TRIM( DATA_DIR ) // 'o3col.TOMS.1999.geos1.4x5'
       ENDIF

C       FILENAME=TRIM( FILENAME )
C       CATEGORY_IN=TRIM( CATEGORY_IN )

       TAU0_IN  = GET_TAU0( themonth, NTOMSyear+1900 )


          print*,TRACER_IN,TAU0_IN
          print*,FILENAME,CATEGORY_IN
       CALL READ_BPCH2( FILENAME, CATEGORY_IN, TRACER_IN, TAU0_IN,
     &                   IIPAR,    JJPAR,       1,         ARRAY )

       ctm4x5month(:,:)=ARRAY(:,:,1)

       ENDIF

       RETURN
       END

!---------------------------------------------------------------------------

      FUNCTION GET_TAU0( MONTH, YEAR ) RESULT( THIS_TAU0 )
!
!*****************************************************************************
!  Function GET_TAU0 returns the corresponding TAU0 value for a given
!  MONTH and YEAR.  This is necessary to index binary punch files.
!  (bmy, 7/20/00, 8/1/00)
!
!  Arguments as Input:
!
!===========================================================================
!  (1) MONTH (INTEGER) : Month number (1-12)
!  (2) YEAR  (INTEGER) : 4-digit year number (e.g. 1985)
!
!  NOTES: 
!  (1) 1985 is the first year of the GEOS epoch.
!  (2) Add TAU0 values for years 1985-2001 (bmy, 8/1/00)
!*****************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN)           :: MONTH
      INTEGER, INTENT(IN), OPTIONAL :: YEAR

      ! Local variables
      INTEGER                       :: TMP_YEAR
      REAL*8                        :: TAU0(12)

      ! Return value
      REAL*8                        :: THIS_TAU0

      !=================================================================
      ! GET_TAU0 begins here!
      !=================================================================

      ! Error checking 
      IF ( MONTH < 1 .or. MONTH > 12 ) THEN
         WRITE( 6, '(a)' ) 'GET_TAU0: Invalid MONTH selection!'
         WRITE( 6, '(a)' ) 'STOP in GET_TAU0_VALUES ("bpch2_mod.f")!'
      ENDIF

      ! If YEAR is not passed, default to 1985 (first year of GEOS epoch)
      IF ( PRESENT( YEAR ) ) THEN
         TMP_YEAR = YEAR
      ELSE
         TMP_YEAR = 1985
      ENDIF
      ! CASE statement for year
      SELECT CASE ( TMP_YEAR )
         CASE ( 1985 )
            TAU0(:) = (/      0d0,    744d0,   1416d0,   2160d0, 
     &                     2880d0,   3624d0,   4344d0,   5088d0, 
     &                     5832d0,   6552d0,   7296d0,   8016d0 /)

         CASE ( 1986 )
            TAU0(:) = (/   8760d0,   9504d0,  10176d0,  10920d0, 
     &                    11640d0,  12384d0,  13104d0,  13848d0, 
     &                    14592d0,  15312d0,  16056d0,  16776d0 /)

         CASE ( 1987 )
            TAU0(:) = (/  17520d0,  18264d0,  18936d0,  19680d0,
     &                    20400d0,  21144d0,  21864d0,  22608d0,
     &                    23352d0,  24072d0,  24816d0,  25536d0 /)
 
         CASE ( 1988 )
            TAU0(:) = (/  26280d0,  27024d0,  27720d0,  28464d0,
     &                    29184d0,  29928d0,  30648d0,  31392d0,
     &                    32136d0,  32856d0,  33600d0,  34320d0 /)

         CASE ( 1989 )
            TAU0(:) = (/  35064d0,  35808d0,  36480d0,  37224d0,
     &                    37944d0,  38688d0,  39408d0,  40152d0,
     &                    40896d0,  41616d0,  42360d0,  43080d0 /)

         CASE ( 1990 )
            TAU0(:) = (/  43824d0,  44568d0,  45240d0,  45984d0,
     &                    46704d0,  47448d0,  48168d0,  48912d0,
     &                    49656d0,  50376d0,  51120d0,  51840d0 /)

         CASE ( 1991 )
            TAU0(:) = (/  52584d0,  53328d0,  54000d0,  54744d0,
     &                    55464d0,  56208d0,  56928d0,  57672d0,
     &                    55464d0,  56208d0,  56928d0,  57672d0 /)

         CASE ( 1992 )
            TAU0(:) = (/  61344d0,  62088d0,  62784d0,  63528d0,
     &                    64248d0,  64992d0,  65712d0,  66456d0,
     &                    67200d0,  67920d0,  68664d0,  69384d0 /)

         CASE ( 1993 )
            TAU0(:) = (/  70128d0,  70872d0,  71544d0,  72288d0,
     &                    73008d0,  73752d0,  74472d0,  75216d0,
     &                    75960d0,  76680d0,  77424d0,  78144d0 /)

         CASE ( 1994 ) 
            TAU0(:) = (/  78888d0,  79632d0,  80304d0,  81048d0,
     &                    81768d0,  82512d0,  83232d0,  83976d0,
     &                    84720d0,  85440d0,  86184d0,  86904d0 /)

         CASE ( 1995 )
            TAU0(:) = (/  87648d0,  88392d0,  89064d0,  89808d0,
     &                    90528d0,  91272d0,  91992d0,  92736d0,
     &                    93480d0,  94200d0,  94944d0,  95664d0 /)

         CASE ( 1996 )
            TAU0(:) = (/  96408d0,  97152d0,  97848d0,  98592d0,
     &                    99312d0, 100056d0, 100776d0, 101520d0,
     &                   102264d0, 102984d0, 103728d0, 104448d0 /)

         CASE ( 1997 )
            TAU0(:) = (/ 105192d0, 105936d0, 106608d0, 107352d0,
     &                   108072d0, 108816d0, 109536d0, 110280d0,
     &                   111024d0, 111744d0, 112488d0, 113208d0 /)

         CASE ( 1998 )
            TAU0(:) = (/ 113952d0, 114696d0, 115368d0, 116112d0,
     &                   116832d0, 117576d0, 118296d0, 119040d0,
     &                   119784d0, 120504d0, 121248d0, 121968d0 /)
        CASE ( 1999 )
            TAU0(:) = (/ 122712d0, 123456d0, 124128d0, 124872d0,
     &                   125592d0, 126336d0, 127056d0, 127800d0,
     &                   128544d0, 129264d0, 130008d0, 130728d0 /)

         CASE ( 2000 )
            TAU0(:) = (/ 131472d0, 132216d0, 132912d0, 133656d0,
     &                   134376d0, 135120d0, 135840d0, 136584d0,
     &                   137328d0, 138048d0, 138792d0, 139512d0 /)

         CASE ( 2001 )
            TAU0(:) = (/ 140256d0, 141000d0, 141672d0, 142416d0,
     &                   143136d0, 143880d0, 144600d0, 145344d0,
     &                   146088d0, 146808d0, 147552d0, 148272d0 /)
 
         CASE DEFAULT
            WRITE( 6, '(a)' ) 'GET_TAU0: Invalid YEAR selection!'
            WRITE( 6, '(a)' ) 'STOP in GET_TAU0 ("bpch2_mod.f")!'
            STOP

      END SELECT

      ! Select the TAU0 value for the given month
      THIS_TAU0 = TAU0( MONTH )

      ! Return to calling program
      END FUNCTION GET_TAU0



