C $Id: CO_readfields.f,v 1.1 2003/06/30 20:26:03 bmy Exp $
      SUBROUTINE CO_readfields(TEMPO,NCLIMATOLOGY,TAU_MONTH)
C
C Created by bnd/bey (12/98).
C   SR reads in DAO monthly average concentration fields.
C
        ! References to F90 modules
        USE BPCH2_MOD

      IMPLICIT NONE
#     include "CMN_SIZE"
#     include "CMN"
#     include "CMN_SETUP"
#     include "CMN_CO"
C
      CHARACTER ( LEN=5)    :: TEMPO
      CHARACTER ( LEN=8 )   :: HEADER(2),HEADER2,HEADEROH
      CHARACTER ( LEN=255 ) :: CHAROP
      INTEGER I,M,J,K,N,JJ
      INTEGER NFIELDID(NFIELDS2)
      REAL*8 TAU_MONTH
      REAL*4 ARRAY(IIPAR,JJPAR,LLPAR)
c      DATA HEADER/'IJ-24H-$','IJ-AVG-$'/
      DATA HEADER/'IJ-AVG-$','IJ-AVG-$'/
cbnd      DATA NFIELDID/1,5,6,9,18,19,21,25,4,22,23,1001/
      DATA NFIELDID/1,5,6,9,18,19,21,25,4,22,23,1/

      HEADEROH ='CHEM-L=$'

C
C NFIELDS2= number of fields (i.e., species) read
C           in SR CO_readfields.f
C
C NFIELDID
C    1  = NOx
C    5  = ALK4
C    6  = ISOP
C    9  = ACET
C   18  = PRPE
C   19  = C3H8
C   21  = C2H6
C   25  = O3
C    4  = CO
C   22  = N2O5
C   23  = HNO4
C         OH
C
      CHAROP = TRIM(TEMP_DIR)//TEMPO
C
      DO JJ=1,NFIELDS2-1
C
       IF(NFIELDID(JJ).EQ.25.AND.NCLIMATOLOGY.EQ.1) GOTO 33
C
         HEADER2=HEADER(1) 
         IF(NFIELDID(JJ).EQ.6) HEADER2=HEADER(2)
C
          print*,JJ,NFIELDID(JJ),HEADER2
C
         CALL READ_BPCH2( CHAROP,HEADER2,NFIELDID(JJ),TAU_MONTH,
     &                  IIPAR,JJPAR,LLPAR,ARRAY )
C
         BBIJ(:,:,:,JJ)=ARRAY(:,:,:) 
C
 33        CONTINUE
C
      ENDDO
C
C OH fields
C
C NOTE: bmy had to read in the monthly files and change mod.
C       OH's tracer number showed up as 1001 instead of 1.
C       Need to get him to do this for other files!
C
      DO JJ=NFIELDS2,NFIELDS2
C
          print*,JJ,NFIELDID(JJ),HEADEROH
C
         CALL READ_BPCH2( CHAROP,HEADEROH,NFIELDID(JJ),TAU_MONTH,
     &                  IIPAR,JJPAR,LLPAR,ARRAY )
C
         BBIJ(:,:,:,JJ)=ARRAY(:,:,:)
C
      ENDDO
C
      RETURN 
      END
