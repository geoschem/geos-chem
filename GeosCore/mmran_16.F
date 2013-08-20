! $Id: mmran_16.f,v 1.1 2009/09/16 14:06:19 bmy Exp $
      SUBROUTINE MMRAN_16( NCB,     NLON,      NLAT,   YLAT, DAY,  
     &                     MONTH,   DAY_OF_YR, CSZA,   TEMP, SFCA, 
     &                     OPTDUST, OPTAER,    MAXBLK, FMAX, ODNEW,   
     &                     KBOT,    KTOP )
!
!******************************************************************************
!  Subroutine MMRAN_16 does the maximum random cloud overlap for 1 to 6 cloud
!  blocks at a time,  and calls PHOTOJ to compute J-Values for one column.
!  (hyl, phs, bmy, 9/18/07, 11/29/07)
!
!  Arguments as Input: 
!  ============================================================================
!  Variable  Type    Dimension  Units   Description
!  --------  ----    ---------  -----   -----------
!  Those for PHOTOJ:
!
!  NLON      INT        -         -     Longitude index
!  NLAT      INT        -         -     Latitude index
!  YLAT      DBLE       -         -     Latitude
!  MONTH     INT        -         -     Month of year (1-12)
!  DAY       INT        -         -     Day of the month
!  DAY_OF_YR INT        -         -     Day of the year
!  CSZA      DBLE       -         -     Cosine of solar zenith angle 
!                                        at nlon, nlat
!  PRES      DBLE       -        [mb]   Column pressure at nlon, nlat
!  TEMP      DBLE    [LMAX]      [K]    Layer temperatures at nlon, nlat 
!  SFCA      DBLE       -         -     Surface albedo at nlon, nlat  
!  OPTDUST   DBLE    [LMAX,NDUST] -     Dust optical depths 
!                                        (for NDUST dust types)
!  OPTAER    DBLE [LMAX,NAER*NRH]  -     Aerosol optical depths
!                                          (for NAER aerosol types)
!
!  and those specifically for MMRAN:
!
!  NCB      INT         -         -    Number of cloud blocks
!  MAXBLK   INT         -         -    Dimension of FMAX, 
!  FMAX     DBLE    [MAXBLK]      -    Largest cloud fraction in block
!  ODNEW    DBLE     [LPAR]       -    In-cloud optical depth
!  KBOT     INT      [LPAR]       -    Index of bottom layer of each block
!  KTOP     INT      [LPAR]       -    Index of top layer of each block
!
! LOCAL VARIABLE:
!  OPTD     DBLE    [LPAR]        -    Layer optical depths at nlon, nlat
!  JSUM     DBLE  [LPAR,JPMAX]    -    accumulate the J-values for the column
!  
!
!  NOTES:
!  (1 ) Remove PRES as an argument, since we no longer need to pass that
!        to PHOTOJ. (bmy, 11/29/07)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      USE CMN_SIZE_MOD, ONLY : NDUST, NAER, NRH, LLPAR
      USE CMN_FJ_MOD, ONLY : JPMAX, JPPJ, JPNL, ZPJ
!      USE JV_CMN_MOD, ONLY : ZPJ

      IMPLICIT NONE

      ! Local variables
      INTEGER, INTENT(IN)    :: NCB ! Number of Cloud Blocks
      INTEGER, INTENT(IN)    :: NLON, NLAT
      REAL*8,  INTENT(IN)    :: CSZA, SFCA, YLAT
      INTEGER, INTENT(IN)    :: DAY, MONTH, DAY_OF_YR
      REAL*8,  INTENT(IN)    :: TEMP(LLPAR)
      REAL*8,  INTENT(IN)    :: OPTDUST(LLPAR,NDUST)
      REAL*8,  INTENT(IN)    :: OPTAER(LLPAR,NAER*NRH)
      INTEGER, INTENT(IN)    :: MAXBLK
      REAL*8,  INTENT(IN)    :: FMAX(MAXBLK)
      REAL*8,  INTENT(IN)    :: ODNEW(LLPAR)
      INTEGER, INTENT(IN)    :: KBOT(LLPAR)
      INTEGER, INTENT(IN)    :: KTOP(LLPAR)

      ! Local variables
      INTEGER    :: II,  JJ,  KK,  LL,  MM,  NN
      INTEGER    :: II2, JJ2,      LL2, MM2, NN2
      REAL*8     :: P1, P2, P3, P4, P5, P6
      REAL*8     :: JSUM(LLPAR,JPMAX)
      REAL*8     :: OPTD(LLPAR)


      !=================================================================
      ! MMRAN_16 begins here!
      !=================================================================

      ! Initialize J-value array
      JSUM = 0d0
 
      ! Initialize Pi
      P1=1d0
      P2=1d0
      P3=1d0 
      P4=1d0
      P5=1d0
      P6=1d0
 
      ! Define the number of loops
      II2 = 1
      JJ2 = 1
      LL2 = 1
      MM2 = 1
      NN2 = 1     
      
      IF ( NCB > 1 ) LL2 = 2         ! At least 2 block-clouds
      IF ( NCB > 2 ) MM2 = 2         ! At least 3 block-clouds
      IF ( NCB > 3 ) NN2 = 2         ! At least 4 block-clouds
      IF ( NCB > 4 ) II2 = 2         ! At least 5 block-clouds
      IF ( NCB > 5 ) JJ2 = 2         ! At least 6 block-clouds


      ! Loop over cloud blocks
      DO KK = 1, 2
      DO LL = 1, LL2
      DO MM = 1, MM2
      DO NN = 1, NN2
      DO II = 1, II2
      DO JJ = 1, JJ2

         ! Zero optical depth
         OPTD(:) = 0d0

         ! 1st cloud block
         IF ( KK == 1 ) THEN
            OPTD(KBOT(1):KTOP(1)) = 0d0
            P1                    = 1d0 - FMAX(1)
         ELSE
            OPTD(KBOT(1):KTOP(1)) = ODNEW(KBOT(1):KTOP(1))
            P1                    = FMAX(1)
         ENDIF


         ! 2nd cloud block
         IF ( NCB > 1 ) THEN
         IF ( LL == 1 ) THEN
            OPTD(KBOT(2):KTOP(2)) = 0d0
            P2                    = 1d0 - FMAX(2)
         ELSE
            OPTD(KBOT(2):KTOP(2)) = ODNEW(KBOT(2):KTOP(2))
            P2                    = FMAX(2)
         ENDIF


         ! 3rd cloud block
         IF ( NCB > 2 ) THEN
         IF ( MM == 1 ) THEN
            OPTD(KBOT(3):KTOP(3)) = 0d0
            P3                    = 1d0 - FMAX(3)
         ELSE
            OPTD(KBOT(3):KTOP(3)) = ODNEW(KBOT(3):KTOP(3))
            P3                    = FMAX(3)
         ENDIF


         ! 4th cloud block
         IF ( NCB > 3 ) THEN         
         IF ( NN == 1 ) THEN
            OPTD(KBOT(4):KTOP(4)) = 0d0
            P4                    = 1d0 - FMAX(4)
         ELSE
            OPTD(KBOT(4):KTOP(4)) = ODNEW(KBOT(4):KTOP(4))
            P4                    = FMAX(4)
         ENDIF


         ! 5th cloud block
         IF ( NCB > 4 ) THEN
         IF ( II == 1 ) THEN
            OPTD(KBOT(5):KTOP(5)) = 0d0
            P5                    = 1d0 - FMAX(5)
         ELSE
            OPTD(KBOT(5):KTOP(5)) = ODNEW(KBOT(5):KTOP(5))
            P5                    = FMAX(5)
         ENDIF


         ! 6th cloud block
         IF ( NCB > 5 ) THEN
         IF ( JJ == 1 ) THEN
            OPTD(KBOT(6):KTOP(6)) = 0d0
            P6                    = 1d0 - FMAX(6)
         ELSE
            OPTD(KBOT(6):KTOP(6)) = ODNEW(KBOT(6):KTOP(6))
            P6                    = FMAX(6)
         ENDIF
         
         ENDIF
         ENDIF
         ENDIF
         ENDIF
         ENDIF

         ! Call the photolysis routine with the OPTD as
         ! computed from the cloud overlaps
         CALL PHOTOJ( NLON, NLAT, YLAT, DAY_OF_YR, MONTH,   DAY,
     &                CSZA, TEMP, SFCA, OPTD,      OPTDUST, OPTAER )

         ! Store the J values into JSUM array
         JSUM(:,:) = JSUM(:,:) +
     &        ( P1 * P2 * P3 * P4 * P5 * P6 * ZPJ(:,:,NLON,NLAT) )
     
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      ! Update J-Values
      ZPJ(:,:,NLON,NLAT) = JSUM(:,:)


      ! Return to caller
      END SUBROUTINE MMRAN_16
