! $Id: jvaluein.f,v 1.1 2003/06/30 20:26:02 bmy Exp $
      SUBROUTINE JVALUEIN
C
C     reads in parameters for jvalue calculation 
C
      USE ERROR_MOD, ONLY : ERROR_STOP

      IMPLICIT NONE

C**********************************************************************
C                                                                     *
C  HARVARD TROPOSPHERIC CHEMISTRY MODULE FOR 3-D APPLICATIONS         *
C  by Larry Horowitz, Jinyou Liang, Gerry Gardner, Prof. Daniel Jacob *
C  of HARVARD UNIVERSITY    (Release V2.0)                            *
C**********************************************************************

#     include "CMN_SIZE"
#     include "comsol.h"
#     include "CMN_SETUP" ! (bmy, 4/3/02)
      CHARACTER*8 TITLE

C Standard atmosphere for air density 0-15-30-45-60N, summer.
      REAL*8 PATM(35),GATM(35),RATM(35),PSTD(NSTDL)

      INTEGER I,J,II,NTATM,NTATZ


      REAL*8 CSTAT0,CBOLTZ,AER2,CPLOG1,DLOGP,CPLOG2,AER1,
     1     HEIGHT,SCALE

      ! Added FILENAME variable (bmy, 4/3/02)
      CHARACTER(LEN=255) :: FILENAME 

      !=================================================================
      ! JVALUEIN begins here!
      !=================================================================

      ! Define file name
      FILENAME = TRIM( DATA_DIR ) // 'slowj_200203/jvalue.dat' 

      ! Echo FILENAME to stdout
      WRITE( 6, 50 ) TRIM( FILENAME )
 50   FORMAT( 'JVALUEIN: Reading ', a )

      ! Now read "jvalue.dat" from DATA_DIR/slowj_200203/ (bmy, 4/3/02)
      OPEN( 42, FILE=TRIM( FILENAME ), STATUS='UNKNOWN', ERR=160 )

C Read DATA TO CONSTANT ARRAYS
      READ(42,201) TITLE
      READ(42,201) TITLE
      READ(42,202) WL
C     
      READ(42,201) TITLE
      READ(42,202) FL
C
C Read in standard atmos. T distribution for 0-15-30-45-60N summer
C  GATM      = ACCELERATION OF GRAVITY AT GROUND LEVEL
C  RATM      = RADIUS OF THE EARTH AT 2km INTERVALS
      READ (42,204) NTATM
      DO 100 J       =   1, NTATM
         READ (42 ,205) PATM(J),GATM(J),RATM(J)
         READ (42 ,202) (TATM(I,J),I=1,NSTDL)
 100  CONTINUE
C Read in standard O3 distribution for 9-15-30-45 N four seasons. For
C 15N, 30N, one season is supposed to represent all.
      READ(42, 204) NTATZ
      DO 110 J       =   1,  NTATZ
         READ(42, 201) TITLE
         READ (42 ,202) (STDO3(I,J), I=1, 24)
C Interpolate upper level O3 exponentially
         SCALE      = STDO3(24,J)/STDO3(23,J)
         DO I     =   24, NSTDL
            STDO3(I,J) = STDO3(24,J) * (SCALE**(REAL(I)-24))
         ENDDO
 110  CONTINUE

C
C Read in aerosol parameters
      READ(42, 201) TITLE
      READ(42, 207) AERSOL
C Get AERXCT
      DO 120 I    = 1, NSTDL
         HEIGHT   = (I-1)*2.D5
         AER1 = AERSOL(1)*EXP(-HEIGHT/AERSOL(2))
         AER2 = AERSOL(3)
         IF (HEIGHT .GT. AERSOL(4)) AER2 = AER2*EXP(-
     2        (HEIGHT-AERSOL(4))/4.0D5)
         AERXCT(I) = AER1 + AER2
 120  CONTINUE

      CLOSE (42) 
      
      ! Define file name
      FILENAME = TRIM( DATA_DIR ) // 'slowj_200203/o3du.dat' 

      ! Echo FILENAME to stdout
      WRITE( 6, 50 ) TRIM( FILENAME )

      ! Now read "jvalue.dat" from DATA_DIR/slowj_200203/ (bmy, 4/3/02)
      OPEN( 42, FILE=TRIM( FILENAME ), STATUS='UNKNOWN', ERR=160 )      

      READ(42,201) TITLE
      READ(42,201) TITLE
      READ(42, 203) O3DU

      CLOSE (42)

C Convert O3DU for Dobson unit to #/cc for use in model
      DO I = 1, 11
         DO J = 1, 12
            O3DU(I,J) = O3DU(I,J)* 2.687D16
         ENDDO
      ENDDO
C Get air density for standard atmosphere
      CBOLTZ  = 1.38D-19 
      DO 150 I = 1, NTATM
C       HYDROSTATIC EQUILIBRIUM  RECALCULATE DENSITIES(DM)
         PSTD(1)=PATM(I)
C IDEAL GAS LAW:
         STDAIR(1,I) = (7.340D21/TATM(1,I)) * (PATM(I)/1013.25D0)
C WE NOW CALCULATE THE EXPONENTIAL DECREASE IN PRESSURE BASED ON THE 
C SCALE HEIGHT OF THE ATMOSPHERE H = RT/GM.
C CSTAT0 IS GM/R, WHERE G IS SURFACE GRAVITY, M IS THE MOLAR WEIGHT OF AIR,
C R IS THE GAS CONSTANT.  THE COEFFICIENT 3.416E-4 IS G0M/R, WHERE
C G0 IS 980.665 cm sec^-2.
         CSTAT0 = 3.416D-4*(GATM(I)/980.665D0)
         DO 140 II=1,NSTDL
C CPLOG2 CORRECTS GRAVITY FOR INCREASING DISTANCE FROM THE EARTH.  THE 
C FACTOR 0.5 IS JUST TO GET THE AVERAGE IN DLOGP.
            CPLOG2=0.5D0*CSTAT0*(RATM(I)/(RATM(I)+(II-1)*2.D5))**2
            IF(II.EQ.1)GOTO 130
C DLOGP IS -DZ/H, WHERE H IS THE SCALE HEIGHT OF THE ATMOSPHERE.  AN 
C AVERAGE VALUE OF H BETWEEN THE TWO LAYERS IS TAKEN.
            DLOGP=(CPLOG1/TATM(II-1,I)+CPLOG2/TATM(II,I))*(-2.D5)
C EXPONENTIAL DROPOFF OF PRESSURE, FROM HYDROSTATIC EQUILIBRIUM:
            PSTD(II)=PSTD(II-1)*EXP(DLOGP)
            STDAIR(II,I)=PSTD(II)/(CBOLTZ*TATM(II,I))
 130        CPLOG1=CPLOG2
 140     CONTINUE
 150  CONTINUE
C Get cross sections XSECT and related information
C ************************
      CALL READCS(.FALSE.)
C ************************
* 
      RETURN
C OPEN ERROR
 160  CONTINUE
      CALL ERROR_STOP( 'Open error in "jvalue.dat"', 'jvaluein.f' )

 201  FORMAT(A8)
 202  FORMAT(8E10.3)
 203  FORMAT(5X, 11F5.1)
 204  FORMAT(10X, I5)
 205  FORMAT(10X, 3E10.3)
 206  FORMAT(9F5.1)
 207  FORMAT(6E10.3)
      END
