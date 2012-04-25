!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: rdlai
!
! !DESCRIPTION: Subroutine RDLAI is used for soil NOx emissions
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RDLAI( JDAY, MONTH, YEAR )
!
! !USES:
!
      USE GRID_MOD, ONLY : GET_XOFFSET, GET_YOFFSET

      USE CMN_SIZE_MOD
      USE CMN_VEL_MOD
      USE CMN_DEP_MOD

      IMPLICIT NONE
#     include "define.h"

!
! !INPUT PARAMETERS:
!
      INTEGER JDAY        ! Simulated day
      INTEGER MONTH       ! Simulated month
      INTEGER YEAR        ! Simulation year
!
! !REVISION HISTORY:
!  Y. Wang, G. Gardner, D. Jacob - Original version (release v2.1)
!  06 Oct 1999 - R. Yantosca     - Be sure to force double precision with 
!                                  the DBLE function and the "D" exponent, 
!                                  wherever necessary
!  25 Jun 2002 - R. Yantosca     - Replace IMX with IIPAR and JMX with JJPAR
!  19 Nov 2009 - M. Barkley      - Included the simulation and LAI years
!EOP
!------------------------------------------------------------------------------
!BOC
!
! ! LOCAL VARIABLES:
!
!------------------------------------------------------------------------------
! update daily the LAIs (Leaf Area Index)                             
!------------------------------------------------------------------------------
! IREG      = Number of landtypes in grid square (I,J)                
! XLAI      = Leaf Area Index of land type element K  (I,J,K)         
!             current month                                           
! XLAI2     = Leaf Area Index of land type element K  (I,J,K)         
!             following month                                         
! XYLAI     = Leaf Area Index of land type element K  (IJLOOP,K)      
!------------------------------------------------------------------------------

      INTEGER STARTDAY(13),ISAVE
      DATA STARTDAY /15,45,74,105,135,166,196,227,258,288,319,349,380/
      DATA ISAVE /0/
      SAVE ISAVE

      INTEGER IMUL
      INTEGER I,J,K,IJLOOP,MM,ITD
      INTEGER IREF,JREF

      ! Need to add I0, J0 as local variables (bmy, 2/11/03)
      INTEGER I0, J0

      ! For MODIS LAI need yearly input (mpb,2009)
      INTEGER YYYY      ! LAI Year
      
      !-----------------------------------------------------------------
      ! RDLAI begins here !
      !-----------------------------------------------------------------

      ! Get nested-grid offsets (bmy, 2/11/03)
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()

      IF (ISAVE.EQ.0) THEN
         ISAVE=1

         ! Must now include the year (mpb,2009)
         CALL FINDMON(JDAY,MONTH,YEAR,MM,YYYY,STARTDAY)
         IF (JDAY.LT.STARTDAY(1)) THEN
            IMUL=365-STARTDAY(12)+JDAY
            ITD = 31
         ELSE
            IMUL=JDAY-STARTDAY(MM)
            ITD = STARTDAY(MM+1) - STARTDAY(MM)
         END IF

         ! Must now include the year (mpb,2009)
         CALL READLAI( MM, YYYY )
         DO J=1,JJPAR
         DO I=1,IIPAR
         DO K=1,IREG(I,J)
            XLAI2(I,J,K) = (XLAI2(I,J,K)-XLAI(I,J,K))/(DBLE(ITD))
            XLAI(I,J,K)=XLAI(I,J,K)+ XLAI2(I,J,K) * DBLE(IMUL)
         END DO
         END DO
         END DO
      ELSE

         ! Must now include the year (mpb,2009)
         CALL FINDMON(JDAY,MONTH,YEAR,MM,YYYY,STARTDAY)
         IF (JDAY.EQ.STARTDAY(MM)) THEN
            ITD = STARTDAY(MM+1) - STARTDAY(MM)

            ! Must now include the year (mpb,2009)
            CALL READLAI( MM , YYYY )

            DO J=1,JJPAR
            DO I=1,IIPAR
            DO K=1,IREG(I,J)
               XLAI2(I,J,K) = (XLAI2(I,J,K)-XLAI(I,J,K))/(DBLE(ITD))
            END DO
            END DO
            END DO
         ELSE
            DO J=1,JJPAR
            DO I=1,IIPAR
            DO K=1,IREG(I,J)
               XLAI(I,J,K)=XLAI(I,J,K)+ XLAI2(I,J,K)
            END DO
            END DO
            END DO
         END IF
      END IF

      IJLOOP = 0
      DO J = 1, JJPAR
         JREF = J + J0
         DO I = 1, IIPAR
            IJLOOP = IJLOOP + 1
            DO K=1,IJREG(IJLOOP)
               IREF = I + I0
               XYLAI(IJLOOP,K)=XLAI(IREF,JREF,K)
            END DO
         END DO
      END DO

      ! Return to calling program
      END SUBROUTINE RDLAI
!EOC
