#ifdef APM
!************************************************************************
! This is the module to read surface albedo from MODIS satellite data
! There are 7 spectral bands, i.e. 
! band 1: 620-670 nm
!      2: 841-876 nm
!      3: 459-479 nm
!      4: 545-565 nm
!      5: 1230-1250 nm
!      6: 1628-1652 nm
!      7: 2105-2155 nm
! Written by
! Xiaoyan Ma 
! SUNY-Albany
! 05/2011
!************************************************************************
      MODULE APM_ALB_MOD

      implicit none

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these variables ...
      PUBLIC :: APM_ALB


      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

       contains

!************************************************************************

       SUBROUTINE APM_ALB(DAY_OF_YR,KYEAR,II,JJ,NBND,ALBD)
!
! !USES:
!
      USE APM_INIT_MOD,  ONLY : DATA_DIR_1x1
      IMPLICIT NONE

      ! Parameters

      INTEGER, INTENT(IN)    :: DAY_OF_YR   ! day of year
      INTEGER, INTENT(IN)    :: KYEAR       ! year
      INTEGER, INTENT(IN)    :: II,JJ
      INTEGER, INTENT(IN)    :: NBND 

!      INTEGER, PARAMETER   :: NBND=7  ! spectral bands from MODIS
!      INTEGER, PARAMETER   :: IM=72 
!      INTEGER, PARAMETER   :: JM=46 

      REAL               :: ALBD(II,JJ,NBND)

      CHARACTER(LEN=4)   :: NYEAR 
      CHARACTER(LEN=3)   :: NDAY
      CHARACTER(LEN=255) :: FILENAME, YPATH
      INTEGER   :: I, J, IB

      !=================================================================
      ! Read MODIS data from disk
      !=================================================================

      WRITE(NDAY,110) DAY_OF_YR
110   format(i3) 
      WRITE(NYEAR,120) KYEAR 
120   format(i4)
!
!Yu+ put MODIS_ALB under APM_data for now
!Ma+ add flexiblity for different year
      YPATH = TRIM(DATA_DIR_1x1)//'/APM_data_201906/MODIS_ALB/'//
     &        NYEAR//'/'
      print*,'YPATH',TRIM(YPATH)
      ! Now prefix the data directory
      IF(II.EQ.72.and.JJ.EQ.46)THEN
         FILENAME = 'modis_surface_albedo_4x5_'//NYEAR//NDAY
         if(filename(30:30).eq.' ') filename(30:30) = '0'
         if(filename(31:31).eq.' ') filename(31:31) = '0'
      ELSEIF(II.EQ.144.and.JJ.EQ.91)THEN
         FILENAME = 'modis_surface_albedo_2x2.5_'//NYEAR//NDAY
         if(filename(32:32).eq.' ') filename(32:32) = '0'
         if(filename(33:33).eq.' ') filename(33:33) = '0'
      ELSE
         WRITE(6,*) "APM_ALBD_mod.f: Need to check"
         STOP
      ENDIF

!      if(filename(29:29).eq.' ') filename(29:29) = '0'
!      if(filename(30:30).eq.' ') filename(30:30) = '0'
      print*,'FILENAME',NYEAR,NDAY,FILENAME

      !-----------------------------
      ! MODIS surface albedo 
      !-----------------------------

      ! Read data
       open(10,file=TRIM(YPATH)//FILENAME, status='old')
       do IB = 1, NBND
       do I = 1, II   
       read(10,100) (ALBD(I,J,IB),J=1,JJ)
c       print*,'ALBD_MOD',IB,I
c       write(*,100) (ALBD(I,J,IB),J=1,JJ)
100     format(10f10.3)
       enddo
       enddo

       close(10)

      END SUBROUTINE APM_ALB

! *****************************************************************************
      END MODULE APM_ALB_MOD
!----------------------------------------------------------------------------------
#endif
