! $Id: rd_prof.f,v 1.1 2009/09/16 14:06:14 bmy Exp $
      subroutine rd_prof(nj2,namfil)
!-----------------------------------------------------------------------
c  Routine to input T and O3 reference profiles
!-----------------------------------------------------------------------
#     include "define.h"

      USE CMN_FJ_MOD, ONLY : JPMAX, JPPJ
      USE JV_CMN_MOD

      IMPLICIT NONE

      integer ia, i, m, l, lat, mon, ntlats, ntmons, n216, nj2
      real*8  ofac,ofak
      character*11 namfil
c
      open(NJ2,file=namfil)
      read(NJ2,'(A)') TITLE0
      write(6,'(1X,A)') TITLE0
      read(NJ2,'(2I5)') NTLATS,NTMONS
      write(6,1000) NTLATS,NTMONS
      N216 = MIN0(216, NTLATS*NTMONS)
      do IA=1,N216
        read(NJ2,'(1X,I3,3X,I2)') LAT, MON
        M = MIN(12, MAX(1, MON))
        L = MIN(18, MAX(1, (LAT+95)/10))
        read(NJ2,'(3X,11F7.1)') (TREF(I,L,M), I=1,41)
        read(NJ2,'(3X,11F7.4)') (OREF(I,L,M), I=1,31)
      enddo
      close(NJ2)
c
c  Extend climatology to 100 km
      ofac=exp(-2.d5/ZZHT)
      do i=32,51
        ofak=ofac**(i-31)
        do m=1,ntmons
          do l=1,ntlats
            oref(i,l,m)=oref(31,l,m)*ofak
          enddo
        enddo
      enddo
      do l=1,ntlats
        do m=1,ntmons
          do i=42,51
            tref(i,l,m)=tref(41,l,m)
          enddo
        enddo
      enddo
      !### DEBUG
      print*, '### in rd_prof.F'
      print*, '### TREF: ', TREF(:,1,1)
      print*, '### OREF: ', OREF(:,1,1)
c
      return
 1000 format(1x,'Data: ',i3,' Lats x ',i2,' Months')
      end
