! $Id: rdland.f,v 1.3 2004/12/02 21:48:39 bmy Exp $
      SUBROUTINE RDLAND

      ! References to F90 modules
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE GRID_MOD,      ONLY : GET_XOFFSET, GET_YOFFSET

      IMPLICIT NONE

C
C***********************************************************************
C                                                                    ***
C     Module contains:  RDLAND, MODIN, RDLAI, FINDMON, READLAI,      ***
C                       RDDRYCF                                      ***
C                                                                    ***
C***********************************************************************
C
C     Read in land types and frkactions(times 1000) from 'vegtype.global'
C
C***********************************************************************
C     Read in gridded:
C      IREG(I,J)      - # of landtypes in grid square
C      ILAND(I,J,LDT) - Land type ID for element LDT =1, IREG(I,J)
C      IUSE(I,J,LDT)  - fraction ((per mil) of gridbox area occupied by 
C                       land type element LDT
C      FRCLND(I,J)    - Land fraction
C     Output non-gridded:
C      IJREG(IJLOOP)
C      IJLAND(IJLOOP,LDT)
C      IJUSE(IJLOOP,LDT)
C***********************************************************************
C Replaced IM with IIPAR and JM with JJPAR (bmy, 6/25/02)
C Now references ERROR_STOP from "error_mod.f" (bmy, 10/15/02)
C***********************************************************************
C

#     include "CMN_SIZE"
#     include "CMN_DEP"
#     include "CMN_VEL"

      INTEGER I,J,K,IJLOOP,IREF,JREF
      
      ! Now make I0, J0 local variables (bmy, 2/11/03)
      INTEGER :: I0, J0

      ! For filename
      CHARACTER(LEN=255) :: FILENAME
 
      !=================================================================
      ! RDLAND begins here!
      !=================================================================

      ! Get nested-grid offsets (bmy, 2/11/03)
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()

      ! File name
      FILENAME = TRIM( DATA_DIR ) // 
     &           'leaf_area_index_200202/vegtype.global'

      WRITE( 6, 50 ) TRIM( FILENAME )
 50   FORMAT( '     - RDLAND: Reading ', a )

      ! Now read "vegtype.global" from the leaf_area_index_200202/
      ! subdirectory of DATA_DIR (bmy, 4/2/02)
      OPEN( 65, FILE=TRIM( FILENAME ), STATUS='OLD',
     &          FORM='FORMATTED',      ERR=700 )

 100  READ(65,101,end=110,ERR=800) I,J,IREG(I,J),
     1     (ILAND(I,J,K),K=1,IREG(I,J)),
     2     (IUSE(I,J,K),K=1,IREG(I,J))
 101  FORMAT(20I4)
      GO TO 100
 110  CONTINUE
      CLOSE (65)
      IJLOOP = 0
      DO 500 J = 1, JJPAR
         JREF = J + J0
         DO 400 I = 1, IIPAR
            FRCLND(I,J) = 1000.
            IREF = I + I0
            IJLOOP = IJLOOP + 1
            IJREG(IJLOOP) = IREG(IREF,JREF)
            DO 300 K=1,IJREG(IJLOOP)
               IJLAND(IJLOOP,K) = ILAND(IREF,JREF,K)
               IJUSE(IJLOOP,K)  = IUSE(IREF,JREF,K)
               IF (IJLAND(IJLOOP,K) .EQ. 0 )
     x              FRCLND(I,J) = FRCLND(I,J) - IJUSE(IJLOOP,K)
 300        CONTINUE
            FRCLND(I,J) = FRCLND(I,J) / 1000.
 400     CONTINUE
 500  CONTINUE
      
      ! Return
      RETURN

      ! Trap File open error
 700  CONTINUE
      CALL ERROR_STOP( 'Error opening "vegtype.global"', 'rdland.f' )
      
      ! Trap file read error
 800  CONTINUE
      CALL ERROR_STOP( 'Error reading "vegtype.global"', 'rdland.f' )

      ! Return to calling program
      END SUBROUTINE RDLAND
