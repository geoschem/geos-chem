!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: set_aer
!
! !DESCRIPTION: Routine SET\_AER sets aerosol/cloud types and define black 
!  carbon profile for the FAST-J photolysis scheme.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_AER( am_I_Root )
!
! !USES:
!
      USE CMN_FJ_MOD, ONLY : JPMAX, JPPJ
      USE JV_CMN_MOD

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! !REMARKS:
!  MX       Number of different types of aerosol to be considered
!  MIEDX    Index of aerosol types in jv_spec.dat - hardwire in here
! 
! 
! !REVISION HISTORY: 
!  01 Jun 1999 - O. Wild     - Initial version
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!  30 Jul 2012 - R. Yantosca - Added ProTeX headers
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer i
c
c  Initialise aerosol index
      do i=1,MX
        MIEDX(i) = 0
      enddo
c
c  Select Aerosol/Cloud types to be used - define types here
c  Each of these types must be listed in the order used by OPMIE.F
      MIEDX(1)  =  3   !  Black carbon absorber
      MIEDX(2)  = 10   !  Water Cloud (Deirmenjian 8 micron)
      MIEDX(3)  = 14   !  Irregular Ice Cloud (Mishchenko)
      MIEDX(4)  = 15   !  Mineral Dust  .15 micron    (rvm, 9/30/00)
      MIEDX(5)  = 16   !  Mineral Dust  .25 micron    (rvm, 9/30/00)
      MIEDX(6)  = 17   !  Mineral Dust  .4  micron    (rvm, 9/30/00)
      MIEDX(7)  = 18   !  Mineral Dust  .8  micron    (rvm, 9/30/00)
      MIEDX(8)  = 19   !  Mineral Dust 1.5  micron    (rvm, 9/30/00)
      MIEDX(9)  = 20   !  Mineral Dust 2.5  micron    (rvm, 9/30/00)
      MIEDX(10) = 21   !  Mineral Dust 4.0  micron    (rvm, 9/30/00)
      MIEDX(11) = 22   !  Tropospheric Sulfate, RH=0  (rvm, bmy, 2/27/02)
      MIEDX(12) = 23   !  Tropospheric Sulfate, RH=50 (rvm, bmy, 2/27/02)
      MIEDX(13) = 24   !  Tropospheric Sulfate, RH=70 (rvm, bmy, 2/27/02)
      MIEDX(14) = 25   !  Tropospheric Sulfate, RH=80 (rvm, bmy, 2/27/02)
      MIEDX(15) = 26   !  Tropospheric Sulfate, RH=90 (rvm, bmy, 2/27/02)
      MIEDX(16) = 29   !  Black Carbon,         RH=0  (rvm, bmy, 2/27/02)
      MIEDX(17) = 30   !  Black Carbon,         RH=50 (rvm, bmy, 2/27/02)
      MIEDX(18) = 31   !  Black Carbon,         RH=70 (rvm, bmy, 2/27/02)
      MIEDX(19) = 32   !  Black Carbon,         RH=80 (rvm, bmy, 2/27/02)
      MIEDX(20) = 33   !  Black Carbon,         RH=90 (rvm, bmy, 2/27/02)
      MIEDX(21) = 36   !  Organic Carbon,       RH=0  (rvm, bmy, 2/27/02)
      MIEDX(22) = 37   !  Organic Carbon,       RH=50 (rvm, bmy, 2/27/02)
      MIEDX(23) = 38   !  Organic Carbon,       RH=70 (rvm, bmy, 2/27/02)
      MIEDX(24) = 39   !  Organic Carbon,       RH=80 (rvm, bmy, 2/27/02)
      MIEDX(25) = 40   !  Organic Carbon,       RH=90 (rvm, bmy, 2/27/02)
      MIEDX(26) = 43   !  Sea Salt (accum),     RH=0  (rvm, bmy, 2/27/02)
      MIEDX(27) = 44   !  Sea Salt (accum),     RH=50 (rvm, bmy, 2/27/02)
      MIEDX(28) = 45   !  Sea Salt (accum),     RH=70 (rvm, bmy, 2/27/02)
      MIEDX(29) = 46   !  Sea Salt (accum),     RH=80 (rvm, bmy, 2/27/02)
      MIEDX(30) = 47   !  Sea Salt (accum),     RH=90 (rvm, bmy, 2/27/02)
      MIEDX(31) = 50   !  Sea Salt (coarse),    RH=0  (rvm, bmy, 2/27/02)
      MIEDX(32) = 51   !  Sea Salt (coarse),    RH=50 (rvm, bmy, 2/27/02)
      MIEDX(33) = 52   !  Sea Salt (coarse),    RH=70 (rvm, bmy, 2/27/02)
      MIEDX(34) = 53   !  Sea Salt (coarse),    RH=80 (rvm, bmy, 2/27/02)
      MIEDX(35) = 54   !  Sea Salt (coarse),    RH=90 (rvm, bmy, 2/27/02)
      
c
c  Ensure all 'MX' types are valid selections
      do i=1,MX
         IF ( am_I_Root ) THEN 
            write(6,1000) MIEDX(i),TITLEA(MIEDX(i))
         ENDIF
         if(MIEDX(i).gt.NAA.or.MIEDX(i).le.0) then
            IF ( am_I_Root ) THEN
               write(6,1200) MIEDX(i),NAA
            ENDIF
          stop
        endif
      enddo
c
c Approximate Black Carbon up to 10 km; surface 200 ng/m3  (Liousse et al)
c Scale: 1 ng/m3 = 1.0d-15 g/cm3 (1.0d-11 g/m2/cm as BREF is in cm))
c
c Simple place-holder profile
      do i=1,51
        BREF(i)=10.d0*1.0d-11
        if(i.gt.6) BREF(i)=0.d0
      enddo
c
      return
 1000 format('Using Aerosol type: ',i3,1x,a)
 1200 format('Aerosol type ',i3,' unsuitable; supplied values must be ',
     &       'between 1 and ',i3)
      END SUBROUTINE SET_AER
