      MODULE FAST_JX_MOD
      IMPLICIT NONE

      PRIVATE

      PUBLIC INPHOT
      PUBLIC PHOTOJ

      CONTAINS


c-----------------------------------------------------------------------
      SUBROUTINE INPHOT( NLAYER, NREACS )
C-----------------------------------------------------------------------
C  Routine to initialise photolysis rate data, called directly from the
C  cinit routine in ASAD. Currently use it to read the JPL spectral data
C  and standard O3 and T profiles and to set the appropriate reaction index.
C-----------------------------------------------------------------------
c
c     IPH       Channel number for reading all data files
c     RAD       Radius of Earth (cm)
c     ZZHT      Effective scale height above top of atmosphere (cm)
c     DATUMX    Maximum opt.depth above which sub layers should be inserted
c     SZAMAX    Solar zenith angle cut-off, above which to skip calculation
c
C-----------------------------------------------------------------------
      ! References to F90 modules (bmy, 6/27/02)
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE FILE_MOD,     ONLY : IU_FASTJ
      USE FJ2008_MOD

      implicit none
c$$$      include 'parm_CTM.f'
c$$$      include 'cmn_metdat.f'
c$$$      include 'cmn_JVdat.f'
      integer  IPH, I, J, K
      real*8, parameter::   PI = 3.141592653589793d0
      real*8, parameter::  PI180=3.141592653589793d0/180.d0
      INTEGER, INTENT(IN) :: NLAYER,NREACS

c Use channel 8 to read files at the moment
c      IPH  = 8
c
c Defaults & constants
      RAD  = 6375.d5
      ZZHT = 5.d5
      ! # of layers to do chemistry
      JPNL  = NLAYER             

      ! # of reactions in chemistry
      JPPJ  = NREACS + 4         

      ! Error check # of layers
      IF ( JPNL > LPAR ) THEN 
         CALL ERROR_STOP( 'JPNL > LPAR!', 'inphot.f' )
      ENDIF

      ! Error check # of rxns
      IF ( JPPJ > JPMAX ) THEN
         CALL ERROR_STOP( 'JPPJ > JPMAX!', 'inphot.f' )
      ENDIF
c Read in labels of photolysis rates required   >>>>> keyed to users chem code
c   this is a tranfer map from the J's automatically calculated in fast-JX
c   onto the names and order in the users chemistry code
      CALL RD_JS(IU_FASTJ,'chemJ2008.d')

      CALL jv_index
c
c Read in fast-J X-sections (spectral data) <<<<<<<<<<<<<< new fast-JX
      CALL RD_XXX(IU_FASTJ,'spec2008.dat')
c
      !Read in aerosol optics to be used to calculate AODs
      CALL RD_AOD(IU_FASTJ,'jv_spec_aod.dat')

c Read in T & O3 climatology                    >>>> general backup clim.
      CALL RD_PROF(IU_FASTJ,'jv_atms.dat')
c
      ! Select Aerosol/Cloud types to be used
      CALL SET_AER
      
      return
      END SUBROUTINE INPHOT

      SUBROUTINE RD_JS(NJ1,NAMFIL)
c jfacta, jlabel,rnames,branch in common block
C-----------------------------------------------------------------------
c  Reread the chem_Js.dat file to map photolysis rate to reaction
c  Read in quantum yield 'jfacta' and fastj2 label 'jlabel'
C-----------------------------------------------------------------------
c
c     jfacta    Quantum yield (or multiplication factor) for photolysis
c     jlabel    Reference label identifying appropriate J-value to use
c     ipr       Photolysis reaction counter - should total 'JVN_'
c
C-----------------------------------------------------------------------
      USE FJ2008_MOD
      implicit none
c      include 'parm_CTM.f'
c      include 'cmn_metdat.f'
c      include 'cmn_JVdat.f'
c
      integer, intent(in) ::  NJ1
      character(*), intent(in) ::  NAMFIL

      integer  IPR, I, J, K
      character*120 CLINE
c
c Reread the chem_Js.dat file to map photolysis rate to reaction
c                     Read in quantum yield jfacta and fastj2 label jlabel
      IPR = 0
      open (NJ1,file=NAMFIL,status='old',form='formatted')
 10   read (NJ1,'(A)',err=20)  CLINE
      if (IPR .eq. JVN_) goto 20

      if (CLINE(2:5).eq.'9999') then
        go to 20
      elseif (CLINE(1:1).eq.'#') then
        go to 10
      elseif (CLINE(5:5).eq.'$') then
        go to 10
      else
        IPR = IPR+1
        read (CLINE(79:83),'(F5.1)') JFACTA(IPR)
        read (CLINE(86:92),'(A7)')   JLABEL(IPR)
        JFACTA(IPR) = JFACTA(IPR)/100.d0
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C  Additional code to read reaction names and branch numbers
C  (ppm, 6/98, bmy, 9/99)     
         read (cline(7:10),"(a4)") rnames(ipr)
         rnames(ipr) = trim(rnames(ipr))
         branch(ipr) = 1
         do i=1,ipr-1
            if (rnames(ipr) == rnames(i)) branch(ipr) = branch(i) + 1
         enddo
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        go to 10
      endif
 20   close(NJ1)

      NRATJ = IPR

C-----------------------------------------------------------------------
c  compare Xsections titles with J-values listed in chem code (jratd.dat)
c  map J-values needed for chemistry (chem_Js.dat) onto the fast-JX rates
c  >>>>>>>>>>>>>>>>current code revised to JPL-02 ver 8.5 (5/05)<<<<<<<<<
c          >>>this must now follow the read in of Xsects, etc<<<
C-----------------------------------------------------------------------

c---Zero / Set index arrays that map Jvalue(j) onto rates
c$$$      do J = 1,JVN_
c$$$        JIND(J) = 0
c$$$      enddo
c$$$      do J = 1,NJVAL
c$$$      do K = 1,NRATJ
c$$$        if (JLABEL(K) .eq. TITLEJ(J)) JIND(K)=J
c$$$      enddo
c$$$      enddo
c$$$
c$$$      write(6,'(a,i4,a)') ' Photochemistry Scheme with ',IPR,' J-values'
c$$$      do K=1,NRATJ
c$$$        J = JIND(K)
c$$$        if (J.eq.0) then
c$$$         write(6,'(i5,a9,f6.2,a,i4,a9)') K,JLABEL(K),JFACTA(K),
c$$$     &         ' has no mapping onto onto fast-JX'
c$$$        else
c$$$         write(6,'(i5,a9,f6.2,a,i4,a9)') K,JLABEL(K),JFACTA(K),
c$$$     &         ' mapped onto fast-JX:',J,TITLEJ(J)
c$$$        endif
c$$$      enddo  

      return
      END SUBROUTINE RD_JS

! $Id: jv_index.f,v 1.1.1.1 2008/01/25 20:30:42 bmy Exp $
      SUBROUTINE JV_INDEX
c DEFPRAT in comode.h?????
!
!******************************************************************************
!  Subroutine JV_INDEX computes the mapping between the CTM indices
!  (from "chem.dat") for J-values to the FAST-J indices (from "ratj.d")
!  for J-values.  (bmy, 10/5/98, 10/16/06)
!
!  NOTES:
!  (1 ) Assumes the ordering of a species with several branches in 
!        "ratj.d" is the same as in "chem.dat".
!  (2 ) Updated comments, cosmetic changes (bmy, 11/15/01)
!  (3 ) NAMESPEC is now NAMEGAS for SMVGEAR II.   We don't need to reference 
!        CMN anymore. Now loop from NCS = 1..NCSGAS (bdf, bmy, 4/8/03)
!  (4 ) Now reset NCS to NCSURBAN after loop (dbm, bmy, 10/16/06)
!******************************************************************************
!
      USE FJ2008_MOD
      USE COMODE_LOOP_MOD
      IMPLICIT NONE

      ! Local variables
      INTEGER            :: I, IFNC, IBRCH, N, NK
      CHARACTER (LEN=4)  :: SPECNAME

      !=================================================================
      ! JV_INDEX begins here!
      !=================================================================

      ! Zero the RINDEX array
      RINDEX(:) = 0

      ! Loop over photolysis rxns (urban chemistry only)
      DO NCS = 1, NCSGAS
      DO I   = 1, NPHOT

         !==============================================================
         ! I        = Index of photo rxns    from "globchem.dat"
         ! NK       = Absolute rxn number (adds offset to I)
         ! SPECNAME = Name of species I,     from "globchem.dat"
         ! IBRCH    = Branch # of species I, from "globchem.dat"
         !==============================================================
         NK       = NRATES(NCS) + I
         SPECNAME = NAMEGAS(IRM(1,NK,NCS)) 
         IFNC     = DEFPRAT(NK,NCS) + 0.01d0
         IBRCH    = 10d0*( DEFPRAT(NK,NCS) - IFNC ) + 0.5d0

         !==============================================================
         ! N      = Index of photolysis reactions as listed in "ratj.d"
         ! RNAMES = Name of species N,            as listed in "ratj.d" 
         ! BRANCH = Branch number of species N,   as listed in "ratj.d" 
         !  
         ! If the species names and branch numbers from both "chem.dat" 
         ! and "ratj.d" match, then store N (the "ratj.d" index) in the 
         ! Ith element of RINDEX.
         !  
         ! Thus, when looping over I (the chem.dat" indices), as is 
         ! done in FJFUNC.F, RINDEX(I) will access the correct J-value 
         ! according to the ordering in "ratj.d".
         !==============================================================
         DO N = 1, JPPJ
            IF ( SPECNAME == RNAMES(N) .and. IBRCH == BRANCH(N) ) THEN
               RINDEX(I) = N

               WRITE ( 6, 100 ) I,         SPECNAME,  IBRCH, 
     &                          RINDEX(I), RNAMES(N), BRANCH(N)
 100           FORMAT('Harvard #: ', i3, 1x, a4, ' Branch: ', i2, 
     &                ' --->  Fast-J #: ', i3, 1x, a4, ' Branch: ',i2 )
               EXIT
            ENDIF
         ENDDO
      ENDDO  
      ENDDO

      ! Reset NCS to NCSURBAN for safety's sake (bmy, 10/16/06)
      NCS = NCSURBAN

      ! Return to calling program      
      END SUBROUTINE JV_INDEX
C-----------------------------------------------------------------------

      SUBROUTINE RD_XXX(NJ1,NAMFIL)
C-----------------------------------------------------------------------
c  Read in wavelength bins, solar fluxes, Rayleigh parameters, 
c      T-dependent X-sections. 

c  >>>current code revised to JPL-02 ver 8.5 (5/05)<<<<<

C-----------------------------------------------------------------------
c     NAMFIL   Name of spectral data file (j2_spec.dat) >> j2 for fast-J2
c     NJ1      Channel number for reading data file
c
c     NJVAL    Number of species to calculate J-values for
c     NWWW     Number of wavelength bins, from 1:NWWW
c     WBIN     Boundaries of wavelength bins
c     WL       Centres of wavelength bins - 'effective wavelength'
c     FL       Solar flux incident on top of atmosphere (cm-2.s-1)
c     QRAYL    Rayleigh parameters (effective cross-section) (cm2)
c     QO2      O2 cross-sections
c     QO3      O3 cross-sections
c     Q1D      O3 => O(1D) quantum yield
c     TQQ      Temperature for supplied cross sections
c     QQQ      Supplied cross sections in each wavelength bin (cm2)
C-----------------------------------------------------------------------
      USE FJ2008_MOD
      implicit none
c$$$      include 'fj2008.h'
      integer, intent(in) :: NJ1
      character(*), intent(in) ::  NAMFIL
      character*7  lpdep(3) 
      integer  I, J, K, NK, IW, NQQQ, NWWW
      character*10, temp

c      TQQ(:,:) = 0.d0

C----------spectral data----set for new format data J-ver8.3------------------
c         note that NJVAL = # J-values, but NQQQ (>NJVAL) = # Xsects read in
c         for 2005a data, NJVAL = 62 (including a spare XXXX) and 
c              NQQQ = 64 so that 4 wavelength datasets read in for acetone
c         note NQQQ is not used outside this subroutine!

      ! SDE 09/20/12: Changed NJVAL to 64 and NQQQ to 66 to accommodate
      ! ClNO2 and OCS
      ! SDE 01/17/13: Changed NJVAL to 66 and NQQQ to 68 to accommodate
      ! gaseous H2SO4

      open (NJ1,FILE=NAMFIL,status='old',form='formatted')
      read (NJ1,100) TITLE0
      read (NJ1,101) NJVAL,NQQQ, NWWW,NW1,NW2
      if (NJVAL.gt.X_ .or. NQQQ.gt.X_) then
        write(6,201) NJVAL,X_
        stop
      endif
      write(6,'(1X,A)') TITLE0
C----J-values:  1=O2, 2=O3P,3=O3D 4=readin Xsects
      read (NJ1,102) (WL(IW),IW=1,NWWW)
      read (NJ1,102) (FL(IW),IW=1,NWWW)
      read (NJ1,102) (QRAYL(IW),IW=1,NWWW)

C---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
      read (NJ1,103) TITLEJ(1),TQQ(1,1), (QO2(IW,1),IW=1,NWWW)
      read (NJ1,103) TITLEJ2,TQQ(2,1), (QO2(IW,2),IW=1,NWWW)
      read (NJ1,103) TITLEJ3,TQQ(3,1), (QO2(IW,3),IW=1,NWWW)

      read (NJ1,103) TITLEJ(2),TQQ(1,2), (QO3(IW,1),IW=1,NWWW)
      read (NJ1,103) TITLEJ2,TQQ(2,2), (QO3(IW,2),IW=1,NWWW)
      read (NJ1,103) TITLEJ3,TQQ(3,2), (QO3(IW,3),IW=1,NWWW)

      read (NJ1,103) TITLEJ(3),TQQ(1,3), (Q1D(IW,1),IW=1,NWWW)
      read (NJ1,103) TITLEJ2,TQQ(2,3), (Q1D(IW,2),IW=1,NWWW)
      read (NJ1,103) TITLEJ3,TQQ(3,3), (Q1D(IW,3),IW=1,NWWW)

      do J = 1,3
        write(6,200) TITLEJ(J),(TQQ(I,J),I=1,3)
      enddo

C---Read remaining species:  X-sections at 2 T_s
      do J = 4,NQQQ

!        write(6,*) 'NQQQ: ', NQQQ
!        write(6,*) 'J:    ', J
        read (NJ1,103) TITLEJ(J),TQQ(1,J),(QQQ(IW,1,J),IW=1,NWWW)
        read (NJ1,103) TITLEJ2,TQQ(2,J),(QQQ(IW,2,J),IW=1,NWWW)
          write(6,200) TITLEJ(J),(TQQ(I,J),I=1,2)
      enddo

c  Reset the titles for NJVAL-1 & NJVAL to be the two acetone J_s
c   61: C3H6O  = Acet-a     (CH3CO + CH3) 
c   62: Q2-Ac  = Acet-b     (CH3 + CO + CH3)

      TITLEJ(NJVAL-1) = 'Acet-a'
      TITLEJ(NJVAL)   = 'Acet-b'

c$$$c---Pressure dependencies
c$$$      read(NJ1,104) npdep
c$$$      do k=1,npdep
c$$$        read(NJ1,105) lpdep(k),(zpdep(iw,k),iw=1,nwww)
c$$$        write(6,201)  lpdep(k),(zpdep(iw,k),iw=1,nwww)
c$$$      enddo
c$$$      read(NJ1,'(A)') TITLE0
c
c$$$c---Zero index arrays
c$$$      do j=1,jppj
c$$$        jind(j)=0
c$$$      enddo
c$$$      do j=1,NJVAL
c$$$        jpdep(j)=0
c$$$      enddo
c$$$c
C---Set mapping index
      do j=1,NJVAL
        do k=1,jppj
          if(jlabel(k).eq.titlej(j)) jind(k)=j
        enddo
c        do k=1,npdep
c          if(lpdep(k).eq.titlej(1,j)) jpdep(j)=k
c        enddo
      enddo
      do k=1,jppj
        if(jfacta(k).eq.0.d0)
     &             write(6,*) 'Not using photolysis reaction ',k
        if(jind(k).eq.0) then
          if(jfacta(k).eq.0.d0) then
            jind(k)=1
          else
            write(6,*) 'Which J-rate for photolysis reaction ',k,' ?'
            stop
          endif
        endif
      enddo

c
C---Read aerosol phase functions:
      read (NJ1,'(5x,i5,2f10.5)') JTAUMX,ATAU,ATAU0
      write(6,*)'JTAUMX'
c      read (NJ1,'(A10,I5,/)') temp,NAA
      read (NJ1,'(A10,I5)') temp,NAA
      read(NJ1, '( a )' ) title0
      write(6,*)'NAA'
      NK=5        ! Fix number of wavelengths at 4
      do j=1,NAA
        read(NJ1,110) TITLEA(j)
        do k=1,NK
          read(NJ1,*) WAA(k,j),QAA(k,j),RAA(k,j),SSA(k,j),
     &                                             (PAA(i,k,j),i=1,8)
c          write(6,*) k
        enddo
      enddo
c
      write(6,*) 'Aerosol phase functions & wavelengths'
      do J=1,NAA
        write(6,'(1x,A8,I2,A,9F8.1)')
     $                   TITLEA(J),J,'  wavel=',(WAA(K,J),K=1,NK)
        write(6,'(9x,I2,A,9F8.4)') J,'  Qext =',(QAA(K,J),K=1,NK)
      enddo
      close(NJ1)
      
  100 format(a)
  101 format(10x,5i5)
  102 format(10x,    6e10.3/(10x,6e10.3)/(10x,6e10.3))
  103 format(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))
  110 format(3x,a20)
  200 format(1x,' x-sect:',a10,3(3x,f6.2))
  201 format(' Number of x-sections supplied to Fast-J2: ',i3,/,
     &       ' Maximum number allowed (X_) only set to: ',i3,
     &       ' - increase in cmn_jv.f')

      return
      END SUBROUTINE RD_XXX
C-----------------------------------------------------------------------
      SUBROUTINE RD_PROF(NJ2,NAMFIL)
C-----------------------------------------------------------------------
c  Routine to input T and O3 reference profiles
C-----------------------------------------------------------------------
      USE FJ2008_MOD
      implicit none
c      include 'parm_CTM.f'
c      include 'cmn_metdat.f'
      integer, intent(in) ::  NJ2
      character(*), intent(in) ::  NAMFIL
c
      integer IA, I, M, L, LAT, MON, NTLATS, NTMONS, N216
      real*8  OFAC, OFAK
      character*72 TITLE3
c
      open (NJ2,file=NAMFIL,status='old',form='formatted')
      read (NJ2,'(A)') TITLE3
      read (NJ2,'(2I5)') NTLATS,NTMONS
      write(6,'(1X,A)') TITLE0
      write(6,1000) NTLATS,NTMONS
      N216  = min(216, NTLATS*NTMONS)
      do IA = 1,N216
        read (NJ2,'(1X,I3,3X,I2)') LAT, MON
        M = min(12, max(1, MON))
        L = min(18, max(1, (LAT+95)/10))
        read (NJ2,'(3X,11F7.1)') (TREF(I,L,M), I=1,41)
        read (NJ2,'(3X,11F7.4)') (OREF(I,L,M), I=1,31)
      enddo
      close (NJ2)

c  Extend climatology to 100 km
      OFAC = exp(-2.d5/5.d5)
      do I = 32,51
        OFAK = OFAC**(I-31)
        do M = 1,NTMONS
        do L = 1,NTLATS
          OREF(I,L,M) = OREF(31,L,M)*OFAK
        enddo
        enddo
      enddo
      do L = 1,NTLATS
      do M = 1,NTMONS
      do I = 42,51
        TREF(I,L,M) = TREF(41,L,M)
      enddo
      enddo
      enddo

      return
 1000 format(1x,'std atmos profiles: ',i3,' lat x ',i2,' mon')
      END SUBROUTINE RD_PROF

C-----------------------------------------------------------------------
C $Id: set_aer.f,v 1.1.1.1 2008/01/25 20:30:40 bmy Exp $
      SUBROUTINE SET_AER
C-----------------------------------------------------------------------
c  Set aerosol/cloud types and define black carbon profile
C-----------------------------------------------------------------------
c     MX       Number of different types of aerosol to be considered
c     MIEDX    Index of aerosol types in jv_spec.dat - hardwire in here
C-----------------------------------------------------------------------
      USE FJ2008_MOD
      IMPLICIT NONE

c#     include "cmn_fj.h"
c#     include "jv_cmn.h"
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
        write(6,1000) MIEDX(i),TITLEA(MIEDX(i))
        if(MIEDX(i).gt.NAA.or.MIEDX(i).le.0) then
          write(6,1200) MIEDX(i),NAA
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
     $       'between 1 and ',i3)
      END SUBROUTINE SET_AER

c<<<<<<<<<<<<<<<<<<<<<<<<end CTM-specific subroutines<<<<<<<<<<<<<<<<<<<




c<<<<<<<<<<<<<<<<<<<<<begin CTM-fastJX linking subroutines<<<<<<<<<<<<<<
      
c-----------------------------------------------------------------------
c      SUBROUTINE PHOTOJ(IDAY,ILNG,JLAT, SOLF,SZA,U0,FREFL,ZPJ)
      SUBROUTINE PHOTOJ( NLON, NLAT, YLAT, DAY_OF_YR, MONTH,   DAY,  
     &              CSZA, TEMPPROF,SFCA_IN, ODCOL_IN,  OPTDUST, OPTAER )
C-----------------------------------------------------------------------
c
c  PHOTOJ is the gateway to fast-JX calculations:
c        only access to CTM 3-D GLOBAL arrays
c        sets up the 1-D column arrays for calculating J's
c v6.1   calculates the optical properties for each wavelength
c          and passes single column, single wavelength data to OPMIE 
c        OPMIE no longer needs aerosol data, just the overall properties:
c           DTAUX() = opt depth of each layer, 
c           POMEGAX(8,) = scat phase fn, includes s-s albedo factor.
C-----------------------------------------------------------------------
      USE FJ2008_MOD
      USE ERROR_MOD
      implicit none
c$$$      include 'parm_CTM.f'
c$$$      include 'parm_MIE.f'
c$$$      include 'cmn_metdat.f'
c$$$      include 'cmn_JVdat.f'
      ! Arguments
      INTEGER, INTENT(IN)    :: DAY,   DAY_OF_YR, MONTH
      INTEGER, INTENT(IN)    :: NLAT,  NLON 
      REAL*8,  INTENT(IN)    :: YLAT,  TEMPPROF(LPAR)

c$$$      real*8, intent(inout)  ::  SZA
	real*8,	 INTENT(IN)	   :: CSZA, SFCA_IN
c$$$      integer,intent(in)  ::  IDAY,ILNG,JLAT
      REAL*8,  INTENT(INOUT) :: OPTDUST(LPAR,NDUST)   !(rvm, bmy, 9/30/00)
      REAL*8,  INTENT(INOUT) :: OPTAER(LPAR,NAER*NRH) !(rvm, bmy, 2/27/02)
c      real*8              ::  ZJ(JVL_,JVN_)   !2-D array of J's indexed to CTM chemistry!
      REAL*8 ::ZJ(LPAR,JPMAX)!from jv_cmn.h
      real*8				::  FREFL,U0,SZA_IN         !fraction of energy reflected
c-----------------------------------------------------------------------

c--------key amtospheric data needed to solve plane-parallel J---------
      real*8, dimension(L1_+1) :: ZHL
c      real*8, dimension(L1_+1) ::  
      integer,dimension(L2_+1) :: JXTRA
      
      real*8             :: DTAUX(L1_),POMEGAX(8,L1_),ODABS,ODRAY
      real*8             :: SOLF,FREFS,FREFI,FJTOP,FJBOT,FSBOT
      real*8             :: FJFLX(L_),FLXD(L1_),FLXD0
      real*8             :: AMF2(2*L1_+1,2*L1_+1)
c------------key SCATTERING arrays for clouds+aerosols------------------
c$$$      real*8             :: OPTX(5),SSAX(5),SLEGX(8,5)
c$$$      real*8             :: OD(5,L1_),SSA(5,L1_),SLEG(8,5,L1_)
      real*8             :: OD600(L1_),  PATH,DENSWP,DENS,RH
c------------key arrays AFTER solving for J's---------------------------
      real*8  FFF(W_,LPAR),VALJ(X_)
      real*8  FLXUP(W_),FLXDN(W_),DIRUP(W_),DIRDN(W_)
      real*8  VALJL(LPAR,NJVAL) !2-D array of J_s returned by JRATET

      integer  I,J,K,k0,L,M,KMIE,KW,NCLD, RATIO(W_)
      real*8   AVGF(L_),XQO3(NB),XQO2(NB),DTAUC,WAVE, TTT

c---flux/heating arrays (along with FJFLX,FLXD,FLXD0)
      real*8  FLXJ(L1_),FFX(W_,L1_),FFXNET(W_,8),FFX0,FXBOT,FABOT
	REAL*8, PARAMETER      :: PI = 3.14159265358979324D0
C-----------------------------------------------------------------------
C---new variables added by jmao
	real*8,  INTENT(INOUT) :: ODCOL_IN(LPAR) 
	real*8 xlo2,xlo3,xlray,xltau
	real*8 QXMIE(MX),XLAER(MX),SSALB(MX)
	real*8 PIRAY(NB),PIAER(MX,NB)
	integer MFIT, J1
      real*8 :: start, finish
c	real POMEGAJ(2*M_,L2_+1)
C-----------------------------------------------------------------------
c      call cpu_time(start)
      ZJ(:,:)  = 0.d0
      FFF(:,:) = 0.d0
      FREFL = 0.d0
      FREFS = 0.d0
      ZZHT = 5.d5
C-----------------------------------------------------------------------
c      call SOLARZ(UTIME,IDAY,YGRD(JLAT),XGRD(ILNG), SZA,U0,SOLF)
      ! Import the cosine of the SZA from the CTM (bmy, 9/10/99)
      U0  = CSZA
      SZA_IN = ACOS(CSZA) * ( 180.0d0 / PI )

      !-----------------------------------------------------------------
      !### If you want to set SZA = 0 degrees for testing,
      !### then uncomment the following lines (bmy, 9/13/99) 
      !U0  = 1.0d0
      !SZA = 0.0d0
      !-----------------------------------------------------------------
cccc solf calculation is from JRATE.f
      ! Return if sun is below the horizon
      IF ( SZA_IN > SZAMAX ) RETURN
c$$$      solf=1.d0-(0.034d0*cos(dble(iday-172)*2.d0*pi/365.d0))
	solf=1.d0-(0.034d0*cos(dble(DAY_OF_YR-172)*2.d0*pi/365.d0))
C-----------------------------------------------------------------------
c---  SOLF = 1.d0   ! this needs to be dropped to include 6.7% annual cycle

c      write(6,*) ' Solar Flux factor:', SOLF

c---check for dark conditions SZA > 98.0 deg => tan ht = 63 km
c                        or         99.                  80 km
      if (SZA_IN .gt. SZAMAX) goto 99
!---------------------------------------------------------------------------
! Prior to 2/13/07:
! We no longer need to pass P to SET_PROF, as we have initialized
! the PJ array in "fast_j.f" (bmy, 2/13/07)
!      CALL SET_PROF( NLON, NLAT, YLAT, MONTH, DAY, 
!     &               P,    T,    SA,   OD,    OPTDUST, OPTAER )
!---------------------------------------------------------------------------
      CALL SET_PROF( NLON, NLAT, YLAT, MONTH,   DAY, 
     &               TEMPPROF, SFCA_IN, ODCOL_IN,  OPTDUST, OPTAER )

      do L = 1,L1_
c$$$        ZHL(L) = ZH(ILNG,JLAT,L)
c$$$	ZHL(L) = ZH(NLON,NLAT,L)
	ZHL(L) = Z(L)
      enddo
        ZHL(L1_+1) = ZHL(L1_) + ZZHT

C-----------------------------------------------------------------------
c      call SPHERE(U0,RAD,ZHL,ZZHT,AMF,L1_)
      call SPHERE2 (U0,RAD,ZHL,AMF2,L1_)
C-----------------------------------------------------------------------

c$$$c---calculate the optical properties (opt-depth, single-scat-alb, phase-fn(1:8)
c$$$c---  at the 5 std wavelengths 200-300-400-600-999 nm for cloud+aerosols
c$$$      do L = 1,L1_
c$$$         OD600(L) = 0.D0
c$$$       do K=1,5
c$$$         OD(K,L)  = 0.d0
c$$$         SSA(K,L) = 0.d0
c$$$        do I=1,8
c$$$         SLEG(I,K,L) = 0.d0
c$$$        enddo
c$$$       enddo
c$$$      enddo
          do L = 1,L1_
             OD600(L) = 0.D0
          enddo
       do L = 1,L_

c$$$         RH = RELH(L)
c$$$c---cloud in layer:  if valid cloud index (now 4:13)
c$$$c---this uses data from met fields /cloud fraction: CLNNDX,CLDLWP
c$$$       NCLD = CLDNDX(ILNG,JLAT,L)
c$$$       if (NCLD .ge. 4 .and. NCLD.le.13) then
c$$$         PATH = CLDLWP(ILNG,JLAT,L)
c$$$         DENSWP = PATH / ZH(ILNG,JLAT,L)
c$$$         call OPTICD (OPTX,SSAX,SLEGX,  PATH,DENSWP,NCLD)
c$$$         do K=1,5
c$$$            OD(K,L)  = OD(K,L)  + OPTX(K)
c$$$            SSA(K,L) = SSA(K,L) + SSAX(K)*OPTX(K)
c$$$           do I=1,8
c$$$            SLEG(I,K,L)=SLEG(I,K,L) + SLEGX(I,K)*SSAX(K)*OPTX(K)
c$$$           enddo
c$$$         enddo
c$$$       endif

c$$$c---use OD of clouds (not aerosols) at 600 nm to determine added layers
c$$$       OD600(L) = OD(4,L)
         OD600(L)=ODCOL_IN(L)
        enddo

c---when combining with Rayleigh and O2-O3 abs, remember the SSA and 
c---  phase fn SLEG are weighted by OD and OD*SSA, respectively.

c---Given the aerosol+cloud OD/layer in visible (600 nm) calculate how to add 
C       additonal levels at top of clouds (now uses log spacing)
C-----------------------------------------------------------------------
      call EXTRAL(OD600,L1_,L2_,N_,JTAUMX,ATAU,ATAU0, JXTRA)
C-----------------------------------------------------------------------

c---set surface reflectance
        RFLECT = SFCA_IN
        RFLECT = max(0.d0,min(1.d0,RFLECT))

C---Loop over all wavelength bins to calc mean actinic flux AVGF(L)
      do K = NW1,NW2
        WAVE = WL(K)
C---Pick nearest Mie wavelength to get scattering properites------------
                               KMIE=1  ! use 200 nm prop for <255 nm
        if( WAVE .gt. 255.d0 ) KMIE=2  ! use 300 nm prop for 255-355 nm
        if( WAVE .gt. 355.d0 ) KMIE=3  ! use 400 nm prop for 355-500 nm
        if( WAVE .gt. 500.d0 ) KMIE=4  ! use 600 nm prop for 500-800 nm
c        if( WAVE .gt. 800.d0 ) KMIE=5  ! use 999 nm prop for 800- nm

C---For Mie code scale extinction at 1000 nm to wavelength WAVEL (QXMIE)
      do I=1,MX
        QXMIE(I) = QAA(KMIE,MIEDX(I))/QAA(5,MIEDX(I))
        SSALB(I) = SSA(KMIE,MIEDX(I))
      enddo
c

c---Combine: Rayleigh scatters & O2 & O3 absorbers to get optical properties
c---values at L1_=L_+1 are a pseudo/climatol layer above the top CTM layer (L_)

        do L = 1,L1_
c         TTT     = TTJ(L)
         TTT    = TJ(L)
         XQO3(L) = FLINT(TTT,TQQ(1,2),TQQ(2,2),TQQ(3,2)
     &                      ,QO3(K,1),QO3(K,2),QO3(K,3))

         XQO2(L) = FLINT(TTT,TQQ(1,1),TQQ(2,1),TQQ(3,1)
     &                      ,QO2(K,1),QO2(K,2),QO2(K,3))

        XLO3=DO3(L)*XQO3(L)
        XLO2=DM(L)*XQO2(L)*0.20948d0
        XLRAY=DM(L)*QRAYL(K)
        do I=1,MX
          XLAER(I)=AER(I,L)*QXMIE(I)
        enddo
        DTAUX(L)=XLO3+XLO2+XLRAY
        do I=1,MX
          DTAUX(L)=DTAUX(L)+XLAER(I)
        enddo
c  Fractional extinction for Rayleigh scattering and each aerosol type
       enddo !This is for the L=1,L1_
        MFIT = 2*M_
	  J1 = NLBATM
        do J=J1,NB
           do k0=1,MX
             POMEGAX(I,J) = SSALB(K0)*XLAER(K0)*PAA(I,KMIE,MIEDX(K0))
           enddo
          POMEGAX(1,J) = POMEGAX(1,J) + 1.0d0 * XLRAY
          POMEGAX(3,J) = POMEGAX(3,J) + 0.5d0 * XLRAY
        enddo

          ! Normalize
          DO L=1,NB
          DO I=1,8
             POMEGAX(I,L) = POMEGAX(I,L)/DTAUX(L)
          ENDDO
          ENDDO

        enddo

c      call cpu_time(start)
C-----------------------------------------------------------------------
c      call OPMIE (K, DTAUX,POMEGAX,U0,RFLECT,AMF,JXTRA,
c     &            AVGF,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0)

      call OPMIE (K, DTAUX,POMEGAX,U0,RFLECT,AMF2,JXTRA,
     &            AVGF,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0)

C-----------------------------------------------------------------------

c----direct(DIR) and diffuse(FLX) fluxes at top(UP) (solar = negative by convention)
c----     also at bottom (DN), does not include diffuse reflected flux.
c        FLXUP(K) =  FJTOP
c        DIRUP(K) = -FLXD0
c        FLXDN(K) = -FJBOT
c        DIRDN(K) = -FSBOT

        do L = 1,LPAR
          FFF(K,L) = FFF(K,L) + SOLF*FL(K)*AVGF(L)
        enddo
c          FREFI = FREFI + SOLF*FL(K)*FLXD0/WAVE
c          FREFL = FREFL + SOLF*FL(K)*FJTOP/WAVE
c          FREFS = FREFS + SOLF*FL(K)/WAVE

c---for each wavelength calculate the flux budget/heating rates:
c  FLXD(L) = direct flux deposited in layer L  [approx = MU0*(F(L+1) -F(L)]
c            but for spherical atmosphere!
c  FJFLX(L) = diffuse flux across top of layer L

c---calculate divergence of diffuse flux in each CTM layer (& t-o-a)
c---     need special fix at top and bottom: 
c---FABOT = total abs at L.B. &  FXBOT = net diffusive flux at L.B.
c        FABOT = (1.d0-RFLECT)*(FJBOT+FSBOT)
c        FXBOT = -FJBOT + RFLECT*(FJBOT+FSBOT)
c        FLXJ(1) = FJFLX(1) - FXBOT
c      do L=2,L_
c        FLXJ(L) = FJFLX(L) - FJFLX(L-1)
c      enddo
c        FLXJ(L_+1) = FJTOP - FJFLX(L_)
cc---calculate net flux deposited in each CTM layer (direct & diffuse):
c        FFX0 = 0.d0
c      do L=1,L1_
c        FFX(K,L) = FLXD(L) - FLXJ(L)
c        FFX0 = FFX0 + FFX(K,L)
c      enddo

c  NB: the radiation level ABOVE the top CTM level is included in these budgets
c      these are the flux budget/heating terms for the column:
c  FFXNET(K,1) = FLXD0        direct(solar) flux dep into atmos (spherical)
c  FFXNET(K,2) = FSBOT        direct(solar) flux dep onto LB (surface) 
c  FFXNET(K,3) = FLXD0+FSBOT  TOTAL solar into atmopshere+surface
c  FFXNET(K,4) = FJTOP        diffuse flux leaving top-of-atmos  
c  FFXNET(K,5) = FFX0         diffuse flux absorbed in atmos 
c  FFXNET(K,6) = FABOT        total (dir+dif) absorbed at LB (surface)
c       these are surface fluxes to compare direct vs. diffuse:
c  FFXNET(K,7) = FSBOT        direct flux dep onto LB (surface) - for srf diags
c  FFXNET(K,8) = FJBOT        diffuse flux dep onto LB (surface)

c      FFXNET(K,1) = FLXD0
c      FFXNET(K,2) = FSBOT
c      FFXNET(K,3) = FLXD0+FSBOT
c      FFXNET(K,4) = FJTOP
c      FFXNET(K,5) = FFX0
c      FFXNET(K,6) = FABOT
c      FFXNET(K,7) = FSBOT
c      FFXNET(K,8) = FJBOT
               


      enddo       ! end loop over wavelength K


c          FREFL = FREFL/FREFS      !calculate reflected flux (energy weighted)
c          FREFI = FREFI/FREFS

c---NB UVB = 280-320 = bins 12:15, UVA = 320-400 = bins 16:17, VIS = bin 18 (++)

C-----------------------------------------------------------------------
      call JRATET(PJ,TJ,FFF, VALJL)
C-----------------------------------------------------------------------
c       DO I = 1, JPNL
c       DO J = 1, JPPJ
c          ZPJ(I,J,NLON,NLAT) = ZJ(I,J)
c       ENDDO
c       ENDDO
c---map the J-values from fast-JX onto ASAD ones (use JIND & JFACTA)
c      do L = 1,JVL_
       do I = 1,JPNL
        do J = 1,JPPJ
          if (JIND(J).gt.0) then 
            ZJ(I,J) = VALJL(I,JIND(J))*JFACTA(J)
            ZPJ(I,J,NLON,NLAT) = ZJ(I,J)
          else
            ZJ(I,J) = 0.d0
            ZPJ(I,J,NLON,NLAT) = ZJ(I,J)
          endif
        enddo
      enddo

c---diagnostics that are NOT returned to the CTM code

C-----------------------------------------------------------------------
c      write(6,*)'fast-JX-(6.1)----PHOTOJ internal print: Atmosphere----'
c---used last called values of DTAUX and POMEGAX, should be 600 nm

c      call JP_ATM(PJ,TJ,DM,DO3,Z,ZZHT,DTAUX,POMEGAX,JXTRA)

C---PRINT SUMMARY of mean intensity, flux, heating rates:
c      write(6,*)
c      write(6,*)'fast-JX(6.1)----PHOTOJ internal print: Mean Intens----'
c      write(6,'(a,5f10.4)') 
c     & ' SUMMARY fast-JX: albedo/SZA/u0/F-incd/F-refl/',
c     &  RFLECT,SZA_IN,U0,FREFI,FREFL

c      write(6,'(a5,18i8)')   ' bin:',(K, K=NW2,NW1,-1)
c      write(6,'(a5,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1)
c      write(6,'(a)') ' ----  100000=Fsolar   MEAN INTENSITY per wvl bin'
c      do L = LPAR,1,-1
c       do K=NW1,NW2
c        RATIO(K) = (1.d5*FFF(K,L)/FL(K))
c       enddo
c        write(6,'(i3,2x,18i8)') L,(RATIO(K),K=NW2,NW1,-1)
c      enddo
      
c      write(6,*)
c      write(6,*)'fast-JX(6.1)----PHOTOJ internal print: Net Fluxes----'
c      write(6,'(a11,18i8)')   ' bin:',(K, K=NW2,NW1,-1)
c      write(6,'(a11,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1)
c      write(6,'(a11,18f8.4)') ' sol in atm',(FFXNET(K,1), K=NW2,NW1,-1)
c      write(6,'(a11,18f8.4)') ' sol at srf',(FFXNET(K,2), K=NW2,NW1,-1)
c      write(6,*) ' ---NET FLUXES--- '
c      write(6,'(a11,18f8.4)') ' sol TOTAL ',(FFXNET(K,3), K=NW2,NW1,-1)
c      write(6,'(a11,18f8.4)') ' dif outtop',(FFXNET(K,4), K=NW2,NW1,-1)
c      write(6,'(a11,18f8.4)') ' abs in atm',(FFXNET(K,5), K=NW2,NW1,-1)
c      write(6,'(a11,18f8.4)') ' abs at srf',(FFXNET(K,6), K=NW2,NW1,-1)
c      write(6,*) ' ---SRF FLUXES--- '
c      write(6,'(a11,18f8.4)') ' srf direct',(FFXNET(K,7), K=NW2,NW1,-1)
c      write(6,'(a11,18f8.4)') ' srf diffus',(FFXNET(K,8), K=NW2,NW1,-1)
c      write(6,'(2a)') '  ---NET ABS per layer:       10000=Fsolar',
c     & '  [NB: values <0 = numerical error w/clouds or SZA>90, colm OK]'
c      do L = LPAR,1,-1
c       do K=NW1,NW2
c        RATIO(K) = 1.d5*FFX(K,L)
c       enddo
c        write(6,'(i9,2x,18i8)') L,(RATIO(K),K=NW2,NW1,-1)
c      enddo
      
C-----------------------------------------------------------------------
      
   99 continue
c      call cpu_time(finish)
c      write(6,'(a,f8.5)'),'time=',finish-start
      return
      END SUBROUTINE PHOTOJ

c<<<<<<<<<<<<<<<<<<<<<<<end CTM-fastJX linking subroutines<<<<<<<<<<<<<<

! $Id: set_prof.f,v 1.1.1.1 2008/01/25 20:30:42 bmy Exp $
!----------------------------------------------------------------------------
! Prior to 2/13/07:
!      SUBROUTINE SET_PROF( NLON, NLAT, YLAT, MONTH, DAY, 
!     &                     P,    T,    SA,   ODCOL, OPTDUST, OPTAER )
!----------------------------------------------------------------------------
      SUBROUTINE SET_PROF( NLON, NLAT, YLAT,  MONTH,   DAY, 
     &                     TEMPPROF2,SFCA_IN2,ODCOL, OPTDUST, OPTAER )
!
!******************************************************************************
!  Subroutine SET_PROF sets up atmospheric profiles required by Fast-J using a
!  doubled version of the level scheme used in the CTM.  First pressure and z* 
!  altitude are defined, then O3 and T are taken from the supplied climatology
!  and integrated to the CTM levels (may be overwritten with values directly 
!  from the CTM, if desired) and then black carbon and aerosol profiles are 
!  constructed. (Oliver Wild, 4/7/99, mje, bmy, 7/14/03, 10/30/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NLON    (INTEGER) : Grid box longitude index               [unitless]
!  (2 ) NLAT    (INTEGER) : Grid box latitude index                [unitless]
!  (3 ) YLAT    (REAL*8)  : Grid box latitude                      [degrees]
!  (4 ) MONTH   (INTEGER) : Current month number                   [1-12]
!  (5 ) DAY     (INTEGER) : Current day of month                   [1-31]
!  (6 ) T       (REAL*8)  : Vertical temperature profile           [K]
!  (7 ) SA      (REAL*8)  : Surface albedo                         [unitless]
!  (8 ) ODCOL   (REAL*8)  : Vertical optical depth profile         [unitless]
!  (9 ) OPTDUST (REAL*8)  : Mineral dust opt. depths (1-D profile) [unitless]
!  (10) OPTAER  (REAL*8)  : Aerosol optical depths (1-D profile)   [unitless]
!
!  Important varables passed via "cmn_fj.h" and "jv_cmn.h"
!  ============================================================================
!  (1 ) PJ     :  Pressure at boundaries of model levels [hPa]
!  (2 ) Z      :  Altitude of boundaries of model levels [cm]
!  (3 ) ODCOL  :  Optical depth at each model level
!  (4 ) MASFAC :  Conversion factor for pressure to column density
!  (5 ) TJ     :  Temperature profile on model grid
!  (6 ) DM     :  Air column for each model level [molecules/cm2])
!  (7 ) DO3    :  Ozone column for each model level [molecules/cm2]
!  (8 ) DBC    :  Mass of Black Carbon at each model level [g/cm3]  
!  (9 ) PSTD   :  Approximate pressures of levels for supplied climatology
!
!  References:
!  ============================================================================
!  TOMS/SBUV MERGED TOTAL OZONE DATA, Version 8, Revision 3.
!  Resolution:  5 x 10 deg.
!
!  Source: http://code916.gsfc.nasa.gov/Data_services/merged/index.html
!
!  Contact person for the merged data product:
!  Stacey Hollandsworth Frith (smh@hyperion.gsfc.nasa.gov)
!
!  NOTES:
!  (1 ) Since we parallelize over columns, T, ODCOL, OPTDUST, and OPTAER
!        are 1-D vectors. In the original code from Oliver Wild, these were 
!        3-D arrays.  Also P and SA are just scalars since we just pass one 
!        surface location at a time w/in the parallel loop. (bmy, 9/13/99)
!  (2 ) Mineral dust profiles are also constructed (rvm, 06/04/00)
!  (3 ) Other aerosol profiles are also constructed (rvm, bmy, 2/27/02)
!  (4 ) Added NLON, NLAT, DAY to the arg list.  Now weight the O3 column by 
!        the observed monthly mean EP-TOMS data.  Also updated comments and 
!        added standard GEOS-CHEM documentation header. (mje, bmy, 7/13/03)
!  (5 ) We don't need to initialize the PJ array with ETAA and ETAB anymore.
!        PJ is now defined in "fast_j.f".  Updated comments. (bmy, 10/30/07)
!******************************************************************************
!
      ! References to F90 modules
      USE TOMS_MOD, ONLY : TOMS, DTOMS1, DTOMS2
      USE FJ2008_MOD
      USE ERROR_MOD
      IMPLICIT NONE

c$$$#     include "cmn_fj.h"  ! IPAR, JPAR, LPAR, JPPJ, JPNL
c$$$#     include "jv_cmn.h"  ! NDUST, NAER, PJ

      ! Argument
      INTEGER, INTENT(IN)    :: DAY, MONTH,  NLAT, NLON
      !-------------------------------------------------------------
      ! Prior to 2/13/07:
      !REAL*8,  INTENT(IN)    :: P,   T(LPAR), SA,   YLAT
      !-------------------------------------------------------------
      REAL*8,  INTENT(IN)    :: SFCA_IN2,   YLAT,  TEMPPROF2(LPAR)
      REAL*8,  INTENT(INOUT) :: ODCOL(LPAR)
      REAL*8,  INTENT(IN)    :: OPTDUST(LPAR,NDUST)  
      REAL*8,  INTENT(IN)    :: OPTAER(LPAR,NAER*NRH)  

      ! Local variables
      INTEGER                :: I, K, L, M, N
      REAL*8                 :: DLOGP,F0,T0,B0,PB,PC,XC,MASFAC,SCALEH
      REAL*8                 :: PSTD(52),OREF2(51),TREF2(51),BREF2(51)
      REAL*8                 :: PROFCOL, DAYTOMS

      !=================================================================
      ! SET_PROF begins here!
      !=================================================================

!-----------------------------------------------------------------------------
! Prior to 2/13/07:
! Now intialize PJ two levels up in "fast_j.f".  ETAA and ETAB are meaningless
! for GEOS-5, so we need to set PJ directly from the GET_PEDGE function.
! (bmy, 2/13/07)
!      ! Use the hybrid pressure formulation (bmy, 8/22/02)
!      DO I = 1, NB
!         PJ(I) = ETAA(I) + ( ETAB(I) * P )
!      ENDDO
!      PJ(NB+1) = 0.d0
!-----------------------------------------------------------------------------

      ! Set up cloud and surface properties
      CALL CLDSRF( ODCOL, SFCA_IN2)

      !=================================================================      
      ! Set up pressure levels for O3/T climatology - assume that value
      ! given for each 2 km z* level applies from 1 km below to 1 km 
      ! above, so select pressures at these boundaries. Surface level 
      ! values at 1000 mb are assumed to extend down to the actual 
      ! P(nslon,nslat).
      !=================================================================      
      PSTD(1) = MAX( PJ(1), 1000.D0 )
      PSTD(2) = 1000.D0 * 10.D0**( -1.D0/16.D0 )
      DLOGP   = 10.D0**( -2.D0/16.D0 )
      DO I = 3, 51
         PSTD(I) = PSTD(I-1) * DLOGP
      ENDDO
      PSTD(52) = 0.D0

      ! Mass factor - delta-Pressure [hPa] to delta-Column [molec/cm2]
      MASFAC = 100.D0 * 6.022D+23 / ( 28.97D0 * 9.8D0 * 10.D0 )

      ! Select appropriate monthly and latitudinal profiles
      ! Now use YLAT instead of Oliver's YDGRD(NSLAT) (bmy, 9/13/99) 
      M = MAX( 1, MIN( 12, MONTH                   ) )
      L = MAX( 1, MIN( 18, ( INT(YLAT) + 99 ) / 10 ) )

      ! Temporary arrays for climatology data
      DO I = 1, 51
	 OREF2(I) = OREF(I,L,M)
	 TREF2(I) = TREF(I,L,M)
	 BREF2(I) = BREF(I)
      ENDDO

      ! Apportion O3 and T on supplied climatology z* levels onto CTM levels 
      ! with mass (pressure) weighting, assuming constant mixing ratio and
      ! temperature half a layer on either side of the point supplied.

      DO I = 1, NB
         F0 = 0.D0
         T0 = 0.D0
         B0 = 0.D0
         DO K = 1, 51
            PC = MIN( PJ(I),   PSTD(K)   )
            PB = MAX( PJ(I+1), PSTD(K+1) )
            IF ( PC .GT. PB ) THEN
               XC = ( PC - PB ) / ( PJ(I) - PJ(I+1) )
               F0 = F0 + OREF2(K)*XC
               T0 = T0 + TREF2(K)*XC
               B0 = B0 + BREF2(K)*XC
            ENDIF
         ENDDO
         TJ(I)  = T0
         DO3(I) = 0.D0
         DO3(I) = F0 * 1.D-6
         DBC(I) = B0
      ENDDO

      !=================================================================
      ! Insert model values here to replace or supplement climatology.
      ! Note that CTM temperature is always used in x-section 
      ! calculations (see JRATET); TJ is used in actinic flux 
      ! calculation only.
      !=================================================================    
      !DO I=1,LPAR
      !   DO3(I) = MY_OZONE(i)       ! Volume Mixing Ratio
      !   TJ(I)  = T(I)              ! Kelvin
      !ENDDO
      !DO3(LPAR+1) = MY_OZONE*EXP()  ! Above top of model (or use climatology)
      !TJ(LPAR+1)  = MY_TEMP(LPAR)   ! Above top of model (or use climatology)

      !=================================================================
      ! Calculate effective altitudes using scale height at each level
      !=================================================================
      Z(1) = 0.D0
      DO I = 1, LPAR
         SCALEH = 1.3806D-19 * MASFAC * TJ(I)
         Z(I+1) = Z(I) - ( LOG( PJ(I+1) / PJ(I) ) * SCALEH )
      ENDDO

      !=================================================================
      ! Add Aerosol Column - include aerosol types here. Currently use 
      ! soot water and ice; assume black carbon x-section of 10 m2/g, 
      ! independent of wavelength; assume limiting temperature for 
      ! ice of -40 deg C.
      !=================================================================
      DO I = 1, LPAR
         AER(1,I) = DBC(I) * 10.D0 * ( Z(I+1) - Z(I) )

         ! Turn off uniform black carbon profile (rvm, bmy, 2/27/02)
         AER(1,I) = 0D0

         IF ( TEMPPROF2(I) .GT. 233.D0 ) THEN
            AER(2,I) = ODCOL(I)
            AER(3,I) = 0.D0
         ELSE
            AER(2,I) = 0.D0
            AER(3,I) = ODCOL(I)
         ENDIF   

         ! Also add in aerosol optical depth columns (rvm, bmy, 9/30/00)
         DO N = 1, NDUST
            AER(3+N,I) = OPTDUST(I,N)	
         ENDDO
        
         ! Also add in other aerosol optical depth columns (rvm, bmy, 2/27/02)
         DO N = 1, NAER*NRH
            AER(3+N+NDUST,I) = OPTAER(I,N)
         ENDDO

      ENDDO

      DO K = 1, MX
         AER(K,LPAR+1) = 0.D0
      ENDDO

      !=================================================================
      ! Calculate column quantities for FAST-J
      !=================================================================
      PROFCOL = 0d0

      DO I = 1, NB

         ! Monthly mean air Column [molec/cm2]
         DM(I)  = ( PJ(I) - PJ(I+1) ) * MASFAC

         ! Monthly mean O3 column [molec/cm2]
         DO3(I) = DO3(I) * DM(I)
    
         ! Monthly mean O3 column [DU] 
         PROFCOL = PROFCOL + ( DO3(I) / 2.69d16 )
      ENDDO
      
      !=================================================================
      ! Now weight the O3 column by the observed monthly mean TOMS.
      ! Missing data is denoted by the flag -999. (mje, bmy, 7/15/03)
      ! 
      ! TOMS/SBUV MERGED TOTAL OZONE DATA, Version 8, Revision 3.
      ! Resolution:  5 x 10 deg.
      !
      ! Methodology (bmy, 2/12/07)
      ! ----------------------------------------------------------------
      ! FAST-J comes with its own default O3 column climatology (from 
      ! McPeters 1992 & Nagatani 1991), which is stored in the input 
      ! file "jv_atms.dat".  These "FAST-J default" O3 columns are used 
      ! in the computation of the actinic flux and other optical 
      ! quantities for the FAST-J photolysis.  
      !
      ! The TOMS/SBUV O3 columns and 1/2-monthly O3 trends (contained 
      ! in the TOMS_200701 directory) are read into GEOS-Chem by routine 
      ! READ_TOMS in "toms_mod.f".  Missing values (i.e. locations where 
      ! there are no data) in the TOMS/SBUV O3 columns are defined by 
      ! the flag -999.  
      ! 
      ! After being read from disk in routine READ_TOMS, the TOMS/SBUV 
      ! O3 data are then passed to the FAST-J routine "set_prof.f".  In 
      ! "set_prof.f", a test is done to make sure that the TOMS/SBUV O3 
      ! columns and 1/2-monthly trends do not have any missing values 
      ! for (lat,lon) location for the given month.  If so, then the 
      ! TOMS/SBUV O3 column data is interpolated to the current day and 
      ! is used to weight the "FAST-J default" O3 column.  This 
      ! essentially "forces" the "FAST-J default" O3 column values to 
      ! better match the observations, as defined by TOMS/SBUV.
      !
      ! If there are no TOMS/SBUV O3 columns (and 1/2-monthly trends) 
      ! at a (lat,lon) location for given month, then FAST-J will revert 
      ! to its own "default" climatology for that location and month.  
      ! Therefore, the TOMS O3 can be thought of as an  "overlay" data 
      ! -- it is only used if it exists.
      !
      ! Note that there are no TOMS/SBUV O3 columns at the higher 
      ! latitudes.  At these latitudes, the code will revert to using 
      ! the "FAST-J default" O3 columns.
      !
      ! As of February 2007, we have TOMS/SBUV data for 1979 thru 2005.  
      ! 2006 TOMS/SBUV data is incomplete as of this writing.  For years
      ! 2006 and onward, we use 2005 TOMS O3 columns.
      !
      ! This methodology was originally adopted by Mat Evans.  Symeon 
      ! Koumoutsaris was responsible for creating the downloading and 
      ! processing the TOMS O3 data files from 1979 thru 2005 in the 
      ! TOMS_200701 directory.
      !=================================================================
      DAYTOMS = 0d0

      IF ( DAY <= 15 ) THEN 

         ! Interpolate O3 to current day (w/in first half of month)
         IF ( TOMS(NLON,NLAT)   > -999d0  .AND.
     &        DTOMS1(NLON,NLAT) > -999d0 ) THEN  
            DAYTOMS = TOMS(NLON,NLAT) + DTOMS1(NLON,NLAT) * ( DAY - 15 )
         ENDIF

      ELSE

         ! Interpolate O3 to current day (w/in 2nd half of month)
         IF ( TOMS(NLON,NLAT)   > -999d0  .AND.
     &        DTOMS2(NLON,NLAT) > -999d0 ) THEN  
            DAYTOMS = TOMS(NLON,NLAT) + DTOMS2(NLON,NLAT) * ( DAY - 15 )
         ENDIF

      ENDIF
      
      ! Scale monthly O3 profile to the daily O3 profile (if available)
      IF ( DAYTOMS > 0d0 ) THEN 
         DO I = 1, NB
            DO3(I) = DO3(I) * ( DAYTOMS / PROFCOL )
         ENDDO
      ENDIF
!### Debug
!      write (987,100) nlon,nlat,toms(nlon,nlat), profcol, daytoms,
!     $     dtoms1(nlon,nlat), dtoms2(nlon,nlat), SUM( DO3(:) / 2.69d16 )
! 100  format(i7,x,i7,x,6(f8.2,x))

      ! Return to calling program
      END SUBROUTINE SET_PROF

C $Id: CLDSRF.f,v 1.1.1.1 2008/01/25 20:30:41 bmy Exp $
      SUBROUTINE CLDSRF( ODCOL, SFCA_IN3)
C-----------------------------------------------------------------------
c  Routine to set cloud and surface properties
C-----------------------------------------------------------------------
C  Add the following input variables for CTM interface (bmy, 9/13/99)
C
C  Variable  Type    Dimensn Units   Description
C  --------  ----    ------- -----   -----------
C  ODCOL     dble     [LPAR]   -     Vertical optical depth profile
C  SA        dble       -      -     Surface Albedo
C-----------------------------------------------------------------------
c     rflect   Surface albedo (Lambertian)
c     odmax    Maximum allowed optical depth, above which they are scaled
c     odcol    Optical depth at each model level
c     odsum    Column optical depth
c     nlbatm   Level of lower photolysis boundary - usually surface
C-----------------------------------------------------------------------
      USE FJ2008_MOD
      IMPLICIT NONE
c$$$#     include "cmn_fj.h"
c$$$#     include "jv_cmn.h"

C=============== INPUT PARAMETERS ======================================
      REAL*8, INTENT(INOUT) :: ODCOL(LPAR)
      REAL*8, INTENT(IN)    :: SFCA_IN3

C=============== LOCAL VARIABLES =======================================
      integer i, j, k
      real*8  odsum, odmax, odtot
c
c Default lower photolysis boundary as bottom of level 1
      nlbatm = 1
c
c I move this part to another place.jmao 09/2008
c$$$c Set surface albedo
c$$$      RFLECT = dble(SFCA_IN3)
c$$$      RFLECT = max(0.d0,min(1.d0,RFLECT))
c
c Zero aerosol column
      do k=1,MX
        do i=1,NB
          AER(k,i) = 0.d0
        enddo
      enddo
c
c Scale optical depths as appropriate - limit column to 'odmax'
      odmax = 200.d0
      odsum =   0.d0
      do i=1,lpar
        odcol(i) = dble(odcol(i))
        odsum    = odsum + odcol(i)
      enddo
      if(odsum.gt.odmax) then
        odsum = odmax/odsum
        do i=1,lpar
          odcol(i) = odcol(i)*odsum
        enddo
        odsum = odmax
      endif
c
c  Use clear-sky conditions
c      do i=1,jpnl
c        odcol(i)=0.d0
c      enddo
c
c Set sub-division switch if appropriate
c$$$      odtot=0.d0
c$$$      jadsub(nb)=0
c$$$      jadsub(nb-1)=0
c$$$      do i=nb-1,1,-1
c$$$        k=2*i
c$$$        jadsub(k)=0
c$$$        jadsub(k-1)=0
c$$$        odtot=odtot+odcol(i)
c$$$        if(odtot.gt.0.d0.and.odcol(i).ne.0.d0.and.
c$$$     $                                     dtausub.gt.0.d0) then
c$$$          if(odtot.le.dtausub) then
c$$$            jadsub(k)=1
c$$$            jadsub(k-1)=1
c$$$          elseif(odtot.gt.dtausub) then
c$$$            jadsub(k)=1
c$$$            jadsub(k-1)=0
c$$$            do j=1,2*(i-1)
c$$$              jadsub(j)=0
c$$$            enddo
c$$$            go to 20
c$$$          endif
c$$$        endif
c$$$      enddo
c$$$ 20   continue
c
      return
      END SUBROUTINE CLDSRF

C-----------------------------------------------------------------------

c<<<<<<<<<<<<<<<<<<<<<begin core fast-J subroutines<<<<<<<<<<<<<<<<<<<<<

c-----------------------------------------------------------------------
      SUBROUTINE JRATET(PPJ,TTJ,FFF, VALJL)
c-----------------------------------------------------------------------
c in:
c        PPJ(L1_+1) = pressure profile at edges
c        TTJ(L1_) = = temperatures at mid-level
c        FFF(K=1:NW, L=1:JVL_) = mean actinic flux 
c out:
c        VALJ(JVL_,NJVAL)  JVL_ = no of levels
C-----------------------------------------------------------------------
      USE FJ2008_MOD
      USE ERROR_MOD
      USE LOGICAL_MOD, ONLY : LPRT
      implicit none

      real*8, intent(in)  ::  PPJ(L1_+1),TTJ(L1_)
      real*8, intent(in)  ::  FFF(W_,LPAR)
      real*8, intent(out) ::  VALJL(LPAR,NJVAL)

      !real*8  :: FLINT             ! external function for X-sections
      real*8  VALJ(X_)          ! temp for call J's at one L
      real*8  QO2TOT, QO3TOT, QO31DY, QO31D, QQQT, TFACT
      real*8  TT,PP,DD,TT200,TFACA,TFAC0,TFAC1,TFAC2,QQQA,QQ2,QQ1A,QQ1B
      integer J,K,L, IV
c
      do L = 1,LPAR    ! master loop over layer = L

c---need temperature and density (for some quantum yields):
c---in this case the Pressures PPJ are defined at the boundaries,
c---                Temperatures in the middle of each layer
        TT   = TTJ(L)
        PP  = (PPJ(L)+PPJ(L+1))*0.5d0
          if (L.eq.1) PP = PPJ(1)
        DD = 7.24e18*PP/TT

        do J = 1,NJVAL
          VALJ(J) = 0.d0
        enddo

        do K = NW1,NW2                    ! Using model 'T's here
           QO3TOT = FLINT(TT,TQQ(1,2),TQQ(2,2),TQQ(3,2)
     &                       ,QO3(K,1),QO3(K,2),QO3(K,3))
           QO2TOT = FLINT(TT,TQQ(1,1),TQQ(2,1),TQQ(3,1)
     &                       ,QO2(K,1),QO2(K,2),QO2(K,3))
           QO31DY = FLINT(TT,TQQ(1,3),TQQ(2,3),TQQ(3,3)
     &                       ,Q1D(K,1),Q1D(K,2),Q1D(K,3))
           QO31D  = QO31DY*QO3TOT
          VALJ(1) = VALJ(1) + QO2TOT*FFF(K,L)
          VALJ(2) = VALJ(2) + QO3TOT*FFF(K,L)
          VALJ(3) = VALJ(3) + QO31D*FFF(K,L)
        enddo

        do J = 4,NJVAL

          if (TQQ(2,J) .gt. TQQ(1,J)) then
           TFACT = max(0.d0,min(1.d0,(TT-TQQ(1,J))/(TQQ(2,J)-TQQ(1,J))))
          else
           TFACT = 0.d0
          endif

          do K = NW1,NW2
            QQQT    = QQQ(K,1,J) + (QQQ(K,2,J) - QQQ(K,1,J))*TFACT
            VALJ(J) = VALJ(J) + QQQT*FFF(K,L)
          enddo

c #52 Methylvinyl ketone   'MeVK  '     q(M) = 1/(1 + 1.67e-19*[M])
          if (TITLEJ(J).eq.'MeVK  ') then
            VALJ(J) = VALJ(J)/(1.0 + 1.67e-19*DD)
          endif
c #55 Methylethyl ketone   MEKeto     q(M) = 1/(1 + 2.0*[M/2.5e19])
          if (TITLEJ(J).eq.'MEKeto') then
            VALJ(J) = VALJ(J)/(1.0 + 0.80E-19*DD)
          endif
c #57 Methyl glyoxal       MGlyxl     q(M) = 1/(1 + 4.15*[M/2.5E19])
          if (TITLEJ(J).eq.'MGlyxl') then
            VALJ(J) = VALJ(J)/(1.0 + 1.66e-19*DD)
          endif

        enddo

      if (TITLEJ(NJVAL-1).eq.'Acet-a') then
c--------------J-ref v8.3 includes Blitz ACETONE q-yields--------------
c---Acetone is a special case:   (as per Blitz et al GRL, 2004)
c---     61 = NJVAL-1 = J1(acetone-a) ==> CH3CO + CH3
c---     62 = NJVAL   = J2(acetone-b) ==> CH3 + CO + CH3
          VALJ(NJVAL-1) = 0.d0
          VALJ(NJVAL)   = 0.d0
c---IV=NJVAL-1 = Xsect (total abs) for Acetone - pre-calc Temp interp factors
        IV    = NJVAL-1
        TFACA = (TT-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV))
        TFACA = max(0.d0, min(1.d0, TFACA))
c---IV=NJVAL = Q2 for Acetone=>(2), specifically designed for quadratic interp.
c---      but force to Q2=0 by 210K
        IV    = NJVAL
        TFAC0 = ( (TT-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV)) )**2
        if (TT .lt. TQQ(1,IV)) then
          TFAC0 = (TT - 210.d0)/(TQQ(1,IV)-210.d0)
        endif
        TFAC0 = max(0.d0, min(1.d0, TFAC0))
c---IV=NJVAL+1 = Q1A for Acetone => (1), allow full range of T = 200K-300K
        IV    = NJVAL+1
        TT200 = min(300.d0, max(200.d0, TT))
        TFAC1 = (TT200-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV))
c---IV=NJVAL+2 = Q1B for Acetone => (1)
        IV    = NJVAL+2
        TFAC2 = (TT200-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV))

c---now integrate over wavelengths
        do K = NW1,NW2
c---NJVAL-1 = Xsect (total abs) for Acetone
          IV   = NJVAL-1
          QQQA = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFACA
c---NJVAL   = Q2 for Acetone=>(2), specifically designed for quadratic interp.
          IV   = NJVAL
          QQ2  = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFAC0
          if (TT .lt. TQQ(1,IV)) then
            QQ2 = QQQ(K,1,IV)*TFAC0
          endif
c---NJVAL+1 = Q1A for Acetone => (1), allow full range of T = 200K-300K
          IV   = NJVAL+1
          QQ1A = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFAC1
c---NJVAL+2 = Q1B for Acetone => (1)   ! scaled to [M]=2.5e19
          IV   = NJVAL+2
          QQ1B = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFAC2
          QQ1B = QQ1B*4.d-20
c---J(61)
          VALJ(NJVAL-1) = VALJ(NJVAL-1)
     &         + FFF(K,L)*QQQA*(1.d0-QQ2)/(QQ1A + QQ1B*DD)
c---J(62)
          VALJ(NJVAL) = VALJ(NJVAL) + FFF(K,L)*QQQA*QQ2

        enddo    !K
c-----------end v-8.3 includes Blitz ACETONE q-yields--------------
      endif

! SDE 02/03/12: ClONO2 photolysis also has T-dependent q-yields

c----Load array of J-values in native order, need to be indexed/scaled
c    by ASAD-related code later: ZPJ(L,JJ) = VALJL(L,JIND(JJ))*JFACTA(JJ)
        do J=1,NJVAL
          VALJL(L,J) = VALJ(J)
        enddo

      enddo    ! master loop over L=1,JVL_
      return
      END SUBROUTINE JRATET


c-----------------------------------------------------------------------
      SUBROUTINE JP_ATM(PPJ,TTJ,DDJ,ZZJ,ZHL,ZZHT1,DTAUX,POMEGAX,JXTRA)
c-----------------------------------------------------------------------
      USE FJ2008_MOD
      implicit none
c$$$      include 'parm_CTM.f'
c-----------------------------------------------------------------------
c--------key amtospheric data needed to solve plane-parallel J----------
      real*8, dimension(L1_) :: TTJ,DDJ,ZZJ,ZHL
      real*8, dimension(L1_) :: PPJ 
      integer,dimension(L2_+1) :: JXTRA
      real*8                 :: ZZHT1
      real*8                 :: DTAUX(L1_),POMEGAX(8,L1_)
c-----------------------------------------------------------------------
      integer  I,J,K,L
      real*8   COLO2,COLO3,ZKM,DELZ,ZTOP

c     write(6,'(4a)') '   L z(km)     p      T   ',
c    & '    d(air)   d(O3)','  col(O2)  col(O3)     d-TAU   SS-alb',
c    & '  g(cos) CTM lyr=>'

          COLO2 = 0.d0
          COLO3 = 0.d0
          ZTOP = ZHL(L1_) + ZZHT

        do L = L1_,1,-1
          COLO2 = COLO2 + DDJ(L)*0.20948d0  
          COLO3 = COLO3 + ZZJ(L)
          DELZ = ZTOP-ZHL(L)
          ZTOP = ZHL(L)
          ZKM = ZHL(L)*1.d-5

c      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)') 
c     &      L,ZKM,PPJ(L),TTJ(L),DDJ(L)/DELZ,ZZJ(L)/DELZ,
c     &      COLO2,COLO3,DTAUX(L),POMEGAX(1,L),POMEGAX(2,L)/3.d0,
c     &      JXTRA(L+L),JXTRA(L+L-1)


        enddo

      return
      END SUBROUTINE JP_ATM



C-----------------------------------------------------------------------

c-----------------------------------------------------------------------
      real*8 FUNCTION FLINT (TINT,T1,T2,T3,F1,F2,F3)
C-----------------------------------------------------------------------
c  Three-point linear interpolation function
C-----------------------------------------------------------------------
      real*8  TINT,T1,T2,T3,F1,F2,F3
      if (TINT .le. T2)  then
        if (TINT .le. T1)  then
          FLINT = F1
        else
          FLINT = F1 + (F2 - F1)*(TINT -T1)/(T2 -T1)
        endif
      else
        if (TINT .ge. T3)  then
          FLINT = F3
        else
          FLINT = F2 + (F3 - F2)*(TINT -T2)/(T3 -T2)
        endif
      endif
      return
      END FUNCTION FLINT


c-----------------------------------------------------------------------
      subroutine SPHERE2(GMU,RAD,ZHL,AMF2,L1_)
C-----------------------------------------------------------------------
c----new v6.2: does AirMassFactors for mid-layer, needed for SZA ~ 90
c  This new AMF2 does each of the half-layers of the CTM separately,
c     whereas the original, based on the pratmo code did the whole layers
c     and thus calculated the ray-path to the CTM layre edges, NOT the middle.
c  Since fast-JX is meant to calculate the intensity at the mid-layer, the
c     solar beam at low sun (interpolated between layer edges) was incorrect.
c  This new model does make some approximations of the geometry of the layers:
c     the CTM layer is split evenly in mass (good) and in height (approx).
c
c  Calculation of spherical geometry; derive tangent heights, slant path
c  lengths and air mass factor for each layer. Not called when
c  SZA > 98 degrees.  Beyond 90 degrees, include treatment of emergent
c  beam (where tangent height is below altitude J-value desired at).
C-----------------------------------------------------------------------
c in:
c     GMU     = MU0 = cos(solar zenith angle)
c     RAD     radius of Earth mean sea level (cm)
c     ZHL(L)  height (cm) of the bottome edge of CTM level L
c     ZZHT    scale height (cm) used above top of CTM (ZHL(L_+1)
c     L1_     dimension of CTM = levels +1 (L+1 = above-CTM level)
c out:
c     AMF2(I,J) = air mass factor for CTM level I for sunlight reaching J
C-----------------------------------------------------------------------

      implicit none
      integer, intent(in) ::   L1_
      real*8, intent(in)  ::   RAD,GMU,ZHL(L1_+1)
      real*8, intent(out) ::   AMF2(2*L1_+1,2*L1_+1)

c     RZ      Distance from centre of Earth to each point (cm)
c     RQ      Square of radius ratios
c     SHADHT  Shadow height for the current SZA
c     XL      Slant path between points

      integer  I, J, K, II, L2
      real*8   XMU1,XMU2,XL,DIFF,SHADHT,RZ(L1_+1)
      real*8   RZ2(2*L1_+1),RQ2(2*L1_+1)
c
c--- must have top-of-atmos (NOT top-of-CTM) defined
c      ZHL(L1_+1) = ZHL(L1_) + ZZHT

        RZ(1) = RAD + ZHL(1)
      do II = 2,L1_+1
        RZ(II)   = RAD + ZHL(II)
      enddo

c---calculate heights for edges of split CTM-layers
      L2 = 2*L1_
      do II = 2,L2,2
        I = II/2
        RZ2(II-1) = RZ(I)
        RZ2(II) = 0.5d0*(RZ(I)+RZ(I+1))
      enddo
        RZ2(L2+1) = RZ(L1_+1)
      do II = 1,L2
        RQ2(II) = (RZ2(II)/RZ2(II+1))**2 
      enddo


c---shadow height for SZA > 90
      if (GMU .lt. 0.0d0)  then
        SHADHT = RZ2(1)/dsqrt(1.0d0-GMU**2)
      else
        SHADHT = 0.d0
      endif

c---up from the surface calculating the slant paths between each level
c---  and the level above, and deriving the appropriate Air Mass Factor
         AMF2(:,:) = 0.d0

      do 16 J = 1,2*L1_+1

c  Air Mass Factors all zero if below the tangent height
        if (RZ2(J) .lt. SHADHT) goto 16
c  Ascend from layer J calculating AMF2s
        XMU1 = abs(GMU)
        do I = J,2*L1_
          XMU2     = dsqrt(1.0d0 - RQ2(I)*(1.0d0-XMU1**2))
          XL       = RZ2(I+1)*XMU2 - RZ2(I)*XMU1
          AMF2(I,J) = XL / (RZ2(I+1)-RZ2(I))
          XMU1     = XMU2
        enddo
c--fix above top-of-atmos (L=L1_+1), must set DTAU(L1_+1)=0
          AMF2(2*L1_+1,J) = 1.d0
c
c  Twilight case - Emergent Beam, calc air mass factors below layer
        if (GMU .ge. 0.0d0) goto 16

c  Descend from layer J 
          XMU1       = abs(GMU)
         do II = J-1,1,-1
          DIFF        = RZ2(II+1)*sqrt(1.0d0-XMU1**2)-RZ2(II)
          if (II.eq.1)  DIFF = max(DIFF,0.d0)   ! filter
c  Tangent height below current level - beam passes through twice
          if (DIFF .lt. 0.0d0)  then
            XMU2      = sqrt(1.0d0 - (1.0d0-XMU1**2)/RQ2(II))
            XL        = abs(RZ2(II+1)*XMU1-RZ2(II)*XMU2)
            AMF2(II,J) = 2.d0*XL/(RZ2(II+1)-RZ2(II))
            XMU1      = XMU2
c  Lowest level intersected by emergent beam
          else
            XL        = RZ2(II+1)*XMU1*2.0d0
            AMF2(II,J) = XL/(RZ2(II+1)-RZ2(II))
            goto 16
          endif
         enddo

   16 continue
      return
      END SUBROUTINE SPHERE2
      

     

C-----------------------------------------------------------------------
      SUBROUTINE EXTRAL(DTAUX,L1X,L2X,NX,JTAUMX,ATAU,ATAU0, JXTRA)
C-----------------------------------------------------------------------
c
c    new version 6.1, add sub-layers (JXTRA) to thick cloud/aerosol layers
c    this version sets up log-spaced sub-layers of increasing thickness ATAU
c
c     DTAUX(L=1:L1X) = Optical Depth in layer L (generally 600 nm OD)
c        This can be just cloud or cloud+aerosol, it is used only to set
c        the number in levels to insert in each layer L
c        Set for log-spacing of tau levels, increasing top-down.
c
c     N.B. the TTAU, etc calculated here are NOT used elsewhere

c---The log-spacing parameters have been tested for convergence and chosen
c---  to be within 0.5% for ranges OD=1-500, rflect=0-100%, mu0=0.1-1.0
c---  use of ATAU = 1.18 and min = 0.01, gives at most +135 pts for OD=100 
c---  ATAU = 1.12 now recommended for more -accurate heating rates (not J's)
C-----------------------------------------------------------------------
c
      implicit none
      integer, intent(in) ::  JTAUMX,L1X,L2X  !index of cloud/aerosol
      integer, intent(in) ::  NX              !Mie scattering array size
      real*8,  intent(in) ::  DTAUX(L1X)      !cloud+3aerosol OD in each layer
      real*8,  intent(in) ::  ATAU,ATAU0
      integer, intent(out)::  JXTRA(L2X+1)    !number of sub-layers to be added
c
      integer JTOTL,I,L,L2 
      real*8  TTAU(L2X+1),DTAUJ, ATAU1,ATAULN,ATAUM,ATAUN1
c
C---Reinitialize arrays
      TTAU(:)  = 0.d0
      JXTRA(:) = 0
c
c---combine these edge- and mid-layer points into grid of size:
c---              L2X+1 = 2*L1X+1 = 2*L_+3
c---calculate column optical depths above each level, TTAU(1:L2X+1)
c---      note that TTAU(L2X+1)=0 and TTAU(1)=total OD
c
c---Divide thick layers to achieve better accuracy in the scattering code
c---In the original fast-J, equal sub-layers were chosen, this is wasteful
c---and this new code (ver 5.3) uses log-scale:  
c---        Each succesive layer (down) increase thickness by ATAU > 1
c---        e.g., if ATAU = 2, a layer with OD = 15 could be divided into
c---        4 sub-layers with ODs = 1 - 2 - 4 - 8
c---The key parameters are:
c---        ATAU = factor increase from one layer to the next
c---        ATAUMN = the smallest OD layer desired
c---        JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
c---These are hardwired below, can be changed, but have been tested/optimized

      ATAU1  = ATAU - 1.d0
      ATAULN = log(ATAU)
        TTAU(L2X+1)  = 0.0d0
      do L2 = L2X,1,-1
        L         = (L2+1)/2
        DTAUJ     = 0.5d0 * DTAUX(L)
        TTAU(L2)  = TTAU(L2+1) + DTAUJ
c---Now compute the number of log-spaced sub-layers to be added in
c---   the interval TTAU(L2) > TTAU(L2+1)
c---The objective is to have successive TAU-layers increasing by factor ATAU >1
c---the number of sub-layers + 1
        if (TTAU(L2) .lt. ATAU0) then
          JXTRA(L2) = 0
        else
          ATAUM    = max(ATAU0, TTAU(L2+1))
          ATAUN1 = log(TTAU(L2)/ATAUM) / ATAULN
          JXTRA(L2) = min(JTAUMX, max(0, int(ATAUN1 - 0.5d0)))
        endif
      enddo

c---check on overflow of arrays, cut off JXTRA at lower L if too many levels
      JTOTL    = L2X + 2
      do L2 = L2X,1,-1
        JTOTL  = JTOTL + JXTRA(L2)
        if (JTOTL .gt. NX/2)  then
          write(6,'(A,2I5,F9.2)') 'N_/L2_/L2-cutoff JXTRA:',NX,L2X,L2
          do L = L2,1,-1
            JXTRA(L) = 0
          enddo
          go to 10
        endif
      enddo
  10  continue

      return
      END SUBROUTINE EXTRAL


C-----------------------------------------------------------------------
      subroutine OPMIE (KW, DTAUX,POMEGAX,U0,RFLECT_IN1,AMF2,JXTRA,
     &                  FJACT,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0)
C-----------------------------------------------------------------------
C  fast-J Mie code for J_s, only uses 8-term expansion, 4-Gauss pts
c in:    
c     KW = wavelength bin # (1:18)  - used only for some diagnostics !!
c     DTAUX(1:L1_) = optical depth of each lauer
c     POMEGAX(1:8,1:L1_) = scattering phase fn (multiplied by s-s abledo)
c     U0  = cos (SZA)
c     RFLECT = Lambertian albedo of surface
c     AMF2(1:2*L1_+1,1:2*L1_+1) = air mass factor (I,L)=wt of layer-I to layer-L
c        AMF2 now does both edges and middle of CTM layers
c     JXTRA(1:L1_) = number 0:J = no. of additional levels to be inserted
c out:
c     FJACT(1:L_) = mean actinic flux(diff+direct) at std CTM levels(mid-layer)
c  (new ver 5.7 diagnostics for fluxes, deposition)  fluxes 'down' are <0
c     FJTOP = diffuse flux out top-of-atmosphere (TAU=0, above top model lauer)
c     FJBOT = diffuse flux onto surface (<0 by definition)
c     FSBOT = direct/solar flux onto surface  (<0 by definition)
c     FJFLX(1:L_) = diffuse flux across top of model layer L
C        this connects with FJBOT = FJFLX(0) & FJTOP = FJFLX(L_+1) (not dim!!)
c     FLXD(1:L_+1) = solar flux deposited in layer L (includes layer above CTM)
c        this should take into account sphericity, and is not just = mu0
c     FLXD0 = sum of solar flux deposited in atmos
c        does NOT include flux on lower surface, does NOT mean absorbed!
C-----------------------------------------------------------------------
c
c     DTAU     Local optical depth of each CTM level
c     TTAU     Optical depth of air vertically above each point (to top of atm)
c     FTAU2     Attenuation of solar beam
c     POMEGAJ  Scattering phase function
c
c---new Ver 5.3 code adds sub-layers (# = JXTRA(L2)) using ATAU as the 
c   factor increase from sub-layer to sub-layer 

C  
c
C-----------------------------------------------------------------------
      USE FJ2008_MOD
      USE ERROR_MOD
      implicit none
c$$$      include 'parm_CTM.f'
c$$$      include 'parm_MIE.f'
c$$$      include 'cmn_JVdat.f'

      integer, intent(in) ::  KW
      real*8, intent(in) ::   DTAUX(L1_),POMEGAX(8,L1_)
      real*8, intent(in) ::   AMF2(2*L1_+1,2*L1_+1)
      real*8, intent(in) ::   U0,RFLECT_IN1
      integer, intent(in) ::  JXTRA(L2_+1)
      real*8, intent(out) ::  FJACT(L_),FJTOP,FJBOT,FSBOT
      real*8, intent(out) ::  FJFLX(L_),FLXD(L1_),FLXD0
c
      integer JNDLEV(L_),JNELEV(L1_)
      integer JADDLV(L2_+1),JADDTO(L2_+1),L2LEV(L2_+1)
      integer JTOTL,I,II,J,K,L,LL,IX,JK,   L2,L2L,L22,LZ,LZZ,NDZ
      integer LZ0,LZ1,LZMID
      real*8   SUMT,SUMJ,FJACT2(L_)
c
      real*8  DTAU(L1_+1),POMEGAJ(2*M_,L2_+1),TTAU(L2_+1)
      real*8  FTAU2(L2_+1),POMEGAB(2*M_),FBTM,FTOP
      real*8  ATAUA,ATAUZ,XLTAU,TAUDN,TAUUP,DTAUJ,TAUBTM,TAUTOP
c--- variables used in mie code-----------------------------------------
      real*8   FJ(N_),POMEGA(2*M_,N_),FZ(N_),ZTAU(N_),ZREFL,ZU0,ZFLUX
      real*8   FJT,FJB, FJFLX0, FLXD2(2*L1_)
      integer  MFIT, ND
c      real*8 :: start, finish
C---Reinitialize arrays
c      call cpu_time(start)
      ZTAU(:)     = 0.d0
      FZ(:)       = 0.d0
      POMEGA(:,:) = 0.d0

C---Set up optical depth DTAU(L)
      do L = 1,L1_
        DTAU(L) = DTAUX(L)
      enddo
        DTAU(L1_+1) = 0.d0

c---Define the total scattering phase fn for each CTM layer L=1:L_+1
c---   from a DTAU-wt_d mix of aerosols, cloud & Rayleigh
C---No. of quadrature pts fixed at 4(M_), expansion of phase fn @ 8
      MFIT = 2*M_
      do L = 1,L1_
        do I = 1,MFIT
          POMEGAJ(I,L) = POMEGAX(I,L)
        enddo
      enddo

c---version 6.2 fix to do mid-layers correctly at large SZA (mp, 6/2008)

C---Calculate attenuated incident beam exp(-TTAU/U0 = DTAU * AirMassFactor)
c---      at the middle & edges of the CTM layers L=1:2*L1_+1
c---  L1_ is top-edge of CTM (ie, L=38 = 2 hPa) which has TAU > 0
c---  note that DTAU(L1_) is optical depth in the FULL CTM layer just above

        FTAU2(:) = 0.d0
        FTAU2(L2_+1) = 1.0d0

      do LL = 1,2*L1_+1
         L = (LL+1)/2
       if (AMF2(LL,LL) .gt. 0.0d0) then
           XLTAU = 0.0d0
        do II = 1,2*L1_+1
           I = (II+1)/2
           XLTAU = XLTAU + 0.5d0*DTAU(I)*AMF2(II,LL)
        enddo
        if (XLTAU .lt. 76.d0) then   ! zero out flux at 1e-33
          FTAU2(LL) = exp(-XLTAU)
        endif
       endif
      enddo


c---calculate direct solar flux deposited in each CTM half-layer: L=1:L2_
c---     use FSBOT for surface flux, cannot do layer above CTM (L_+1)
          FLXD2(:) = 0.d0
      do LL = 1,2*L1_
        if (AMF2(LL,LL) .gt. 0.d0) then 
          FLXD2(LL) = (FTAU2(LL+1) - FTAU2(LL))/AMF2(LL,LL)
        endif
      enddo
        if (AMF2(1,1) .gt. 0.d0) then 
          FSBOT = FTAU2(1)/AMF2(1,1)
        else
          FSBOT = 0.d0
        endif

      do LL = 2,2*L1_,2
         L=LL/2
         FLXD(L) = FLXD2(LL)+FLXD2(LL-1)
      enddo

c---integrate solar flux depositied in CTM layers L=1:L_, cannot do top layer
c---  note FLXD0 .ne. (1.d0 - FTAU(L_+1))/AMF(L_+1,L_+1) with spherical atmos
        FLXD0 = 0.d0
      if (AMF2(2*L1_,2*L1_) .gt. 0.d0) then
       do L=1,L1_
        FLXD0 = FLXD0 + FLXD(L)
       enddo
      endif

C------------------------------------------------------------------------
c  Take optical properties on CTM layers and convert to a photolysis
c  level grid corresponding to layer centres and boundaries. This is
c  required so that J-values can be calculated for the centre of CTM
c  layers; the index of these layers is kept in the JNDLEV array.
C------------------------------------------------------------------------
c
c---Now combine the CTM layer edges (1:L_+2) with the CTM mid-layer
c---    points (1:L_) plus 1 for the mid point of added top layer.

c---combine these edge- and mid-layer points into grid of size:
c---              L2_+1 = 2*L1_+1 = 2*L_+3
c---calculate column optical depths above each level, TTAU(1:L2_+1)
c---      note that TTAU(L2_+1)=0 and TTAU(1)=total OD
        TTAU(L2_+1) = 0.0d0
      do L2 = L2_,1,-1
        L          = (L2+1)/2
        DTAUJ      = 0.5d0 * DTAU(L)
        TTAU(L2)   = TTAU(L2+1) + DTAUJ
      enddo
c
c----solar flux incident on lower boundary & Lambertian reflect factor:
      if (FSBOT .gt. 0.d0) then
c---        FSBOT = U0*FTAU(1)
c---        ZFLUX = U0*FTAU(1)*RFLECT/(1.d0+RFLECT)
        ZFLUX = FSBOT*RFLECT_IN1/(1.d0+RFLECT_IN1)
      else
        ZFLUX = 0.d0
      endif

c  Calculate scattering properties, level centres then level boundaries
c ***be careful of order, we are overwriting/shifting the 'POMEGAJ' upward in index***
      do L2 = L2_,2,-2
        L   = L2/2
        do I = 1,MFIT
          POMEGAJ(I,L2) = POMEGAJ(I,L)
        enddo
      enddo
c---lower boundary value is set (POMEGAJ(I,1), but set upper:
        do I = 1,MFIT
          POMEGAJ(I,L2_+1) = POMEGAJ(I,L2_)
        enddo
c---now have POMEGAJ filled at even points from L2=3:L2_-1
c---use inverse interpolation for correct tau-weighted values at edges
      do L2 = 3,L2_-1,2
        TAUDN = TTAU(L2-1)-TTAU(L2)
        TAUUP = TTAU(L2)-TTAU(L2+1)
        do I = 1,MFIT
          POMEGAJ(I,L2) = (POMEGAJ(I,L2-1)*TAUDN + 
     &             POMEGAJ(I,L2+1)*TAUUP) / (TAUDN+TAUUP)
        enddo
      enddo

C---at this point FTAU2(1:L2_+1) and POMEAGJ(1:8, 1:L2_+1)
c---    where FTAU2(L2_+1) = 1.0 = top-of-atmos, FTAU2(1) = surface

C------------------------------------------------------------------------
c  Calculate cumulative total and define levels we want J-values at.
c  Sum upwards for levels, and then downwards for Mie code readjustments.
c
c     JXTRA(L2)  Number of new levels to add between (L2) and (L2+1)
c           ***JXTRA(1:L2_+1) is calculated based on the aerosol+cloud OD_s
c     JADDLV(L2)  Number of new levels actually added at each wavelength
c            where JADDLV = 0 when there is effectively no FTAU2 
c     JADDTO(L2)   Total number of new levels to add to and above level (L2)
c     JNDLEV(L) = L2 index that maps on CTM mid-layer L
c
C------------------------------------------------------------------------

c---JADDLV(L2=1:L2_) = number of levels to add between TTAU2(L2) and TTAU(L2+1)
c---    JADDLV is taken from JXTRA, which is based on visible OD.
c---    JADDTO(L2=1:L2_+1) is the cumulative number of levels to be added
c---note that JADDLV and JADDTO will change with wavelength and solar zenith

c--now try to insert additional levels for thick clouds, ONLY IF FTAU2 > 1.e-8
c-- this will cut off additional levels where the solar beam is negligible.

c---new v5.6-----keep all wavelengths the same for now
c      do L2 = 1,L2_,1
c        if (FTAU2(L2+1) .lt. 1.d-30) then
c          JADDLV(L2) = 0
c        else
c          JADDLV(L2) = JXTRA(L2)
c        endif
c      enddo
      do L2 = 1,L2_,1
        JADDLV(L2) = JXTRA(L2)
      enddo
        JADDTO(L2_+1) = 0
      do L2 = L2_,1,-1
        JADDTO(L2) = JADDTO(L2+1) + JADDLV(L2)
      enddo

c---expanded grid now included CTM edge and mid layers plus expanded 
c---    grid to allow for finer delta-tau at tops of clouds.
c---    DIM of new grid = L2_ + JADDTO(1) + 1

c---L2LEV(L2) = L2-index for old level L2 in expanded J-grid (w/JADDLV)
c     in absence of JADDLV, L2LEV(L2) = L2
        L2LEV(1)  = 1
      do L2 = 2,L2_+1
        L2LEV(L2) = L2LEV(L2-1) + 1 + JADDLV(L2-1)
      enddo

c---JNDLEV(L=1:L_) = L2-index in expanded grid for CTM mid-layer L
c---JNELEV(L=1:L_) = L2-index for top of layer L
      do L = 1,L_
        JNDLEV(L) = L2LEV(2*L)
        JNELEV(L) = L2LEV(2*L+1)
      enddo
        JNELEV(L_+1) = 0  !need to set this to top-of-atmosphere
C---------------------SET UP FOR MIE CODE-------------------------------
c
c  Transpose the ascending TTAU grid to a descending ZTAU grid.
c  Double the resolution - TTAU points become the odd points on the
c  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.
c  Odd point added at top of grid for unattenuated beam   (Z='inf')
c
c  The following mapping holds for JADDLV=0
c        Surface:   TTAU(1)    ==> ZTAU(2*L2_+1)
c        Top:       TTAU(L2_)  ==> ZTAU(3)
c        Infinity:     0.0     ==> ZTAU(1)
c        index: 2*(L2_+1-L2)+1 ==> LZ
c
c  Mie scattering code only used from surface to level L2_
C------------------------------------------------------------------------
c
C------------------------------------------------------------------------
c  Insert new levels, working downwards from the top of the atmosphere
c  to the surface (down in 'LZ', up in 'L2'). This allows ztau and pomega
c  to be incremented linearly, and the flux fz to be attenuated top-down 
c    (avoiding problems where lower level fluxes are zero).
C------------------------------------------------------------------------
c
c  Ascend through atmosphere transposing grid and adding extra points
c  remember L2=1 is surface of CTM, but last layer (LZ) in scattering code.
c  there are twice the number of layers in the LZ arrays (2*L2_ + 2*JADDTO + 1)
c    because we need to insert the intermediate layers (even LZ) for the 
c    asymmetric scattering code.


c  Transfer the L2=1:L2_+1 values (TTAU,FTAU2,POMEGAJ) onto the reverse
c    order, expanded, doubled-level scatter grid. 
c    Note that we need to deal with the expansion by JADD levels (L2L).
c      These JADDLV levels are skipped and need to be interpolated later.
c    Note that only odd LZ levels are filled, 

      NDZ = 2*L2_ + 2*JADDTO(1) + 1

c   Note that the successive sub-layers have the ratio in OD of ATAU
c      ATAUA = (ATAU - 1.d0)/ATAU     ! this is the limit for L22=>inf

      do L2 = 1,L2_+1          ! L2 = index of CTM edge- and mid-layers
        L2L = L2LEV(L2)        ! L2L = index for L2 in expanded scale(JADD)
        LZ  = NDZ + 2 - 2*L2L  ! LZ = index for L2 in scatt arrays
          ZTAU(LZ) = TTAU(L2)
          FZ(LZ)   = FTAU2(L2)
        do I=1,MFIT
          POMEGA(I,LZ) = POMEGAJ(I,L2)
        enddo
      enddo

c   Now go thru the pairs of L2 levels to see if we need JADD levels
      do L2 = 1,L2_             ! L2 = index of CTM edge- and mid-layers
        L2L = L2LEV(L2)         ! L2L = index for L2 in expanded scale(JADD)
        LZ  = NDZ + 2 - 2*L2L   ! LZ = index for L2 in scatt arrays
        L22 = L2LEV(L2+1) - L2LEV(L2) - 1   ! L22 = 0 if no added levels
       if (L22 .gt. 0) then
          TAUBTM = TTAU(L2)
          TAUTOP = TTAU(L2+1)
          FBTM   = FTAU2(L2)
          FTOP   = FTAU2(L2+1)
         do I = 1,MFIT
          POMEGAB(I) = POMEGAJ(I,L2)
         enddo

c---to fit L22 new layers between TAUBOT > TAUTOP, calculate new 1/ATAU factor
c---  such that TAU(just above TAU-btm) = ATUAZ * TAUBTM < TAUBTM

          ATAUZ = exp(-log(TAUBTM/max(TAUTOP,ATAU0))/float(L22+1))

        do L = 1,L22           ! add odd levels between L2LEV(L2) & L2LEV(L2+1)
          LZZ = LZ - 2*L       ! LZZ = index(odd) of added level in scatt arrays
          ZTAU(LZZ) = TAUBTM * ATAUZ

          ATAUA=(TAUBTM-ZTAU(LZZ))/(TAUBTM-TAUTOP) !fraction from TAUBTM=>TAUTOP
          if(IT_IS_NAN(ATAUA)) then
             CALL ERROR_STOP('ATAUA is NAN','LINE2134')
           endif
c---version 5.6 fix   (mp, 3/2007)
c---solar flux at interp-levels: use exp(TAU/U0) if U0>0.02 (89 deg), 
c---   else scale by TAU
          if (U0 .gt. 0.02d0) then
            FZ(LZZ) = FTOP * exp((TAUTOP-ZTAU(LZZ))/U0)
          else
            if (FBTM .lt. 1.d-32) then
              FZ(LZZ) = 0.d0
            else    
              FZ(LZZ) = FBTM * (FTOP/FBTM)**ATAUA
            endif
          endif

          do I = 1,MFIT
            POMEGA(I,LZZ) = POMEGAB(I) + 
     &               ATAUA*(POMEGAJ(I,L2+1)-POMEGAB(I))
          enddo
            TAUBTM       = ZTAU(LZZ)
            FBTM         = FZ(LZZ)
          do I = 1,MFIT
            POMEGAB(I) = POMEGA(I,LZZ)
          enddo
        enddo
       endif
      enddo

c   Now fill in the even points with simple interpolation in scatter arrays:
      do LZ = 2,NDZ-1,2
        ZTAU(LZ) = 0.5d0*(ZTAU(LZ-1)+ZTAU(LZ+1))
        FZ(LZ)   = sqrt(FZ(LZ-1)*FZ(LZ+1))
       do I=1,MFIT
        POMEGA(I,LZ) = 0.5d0*(POMEGA(I,LZ-1)+POMEGA(I,LZ+1))
       enddo
      enddo
      
      ND = NDZ
      ZU0 = U0
      ZREFL = RFLECT_IN1

c---PRINT diagnostics
c----now check integral of FZ over each CTM layer to ensure that it equals FLXD
c      if (KW.eq.18) then
c       write(6,'(A,3I6)') 'OPMIE levels: L,TAU,Fs,pomega',KW, ND,N_
c      do L=1,ND
c       write(6,'(i5,f15.5,1p,e15.5,0p,10f10.5)') L,
c     & ZTAU(L),FZ(L),(POMEGA(I,L),I=1,3) 
c      enddo
c
c      write(6,'(a,2i6)') '  net flux in each CTM layer',KW,L_
c      write(6,'(2i5,a8,1p,e12.5)')(L,JNELEV(L),' flxd ',FLXD(L),L=1,L1_)
c
c      endif

      if(ND .gt. N_) then
        write(6,'(a,2i9)') ' overflow of scatter arrays:',ND,N_
        stop
      endif

C-----------------------------------------------------------------------
      call MIESCT(FJ,FJT,FJB,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,ZU0,MFIT,ND)
C-----------------------------------------------------------------------

c---Move mean intensity from scatter array FJ(LZ=1:ND) 
c---              to CTM mid-level array FJACT(L=1:L_)

c---mean intensity:  4*<I> + solar at mid-layer
      do L = 1,L_
        L2L = JNDLEV(L)
        LZ  = ND+2 - 2*L2L
        FJACT(L) = 4.d0*FJ(LZ) + FZ(LZ)
      enddo

c---mean diffuse flux:  4<I*mu> (not solar) at top of layer L
c---      average (tau-wtd) the h's just above and below the L-edge
      do L = 1,L_
        L2L = JNELEV(L)
        LZ  = ND+2 - 2*L2L
c---       FJFLX(L) = 2.0d0*(FJ(LZ-1) + FJ(LZ+1))
         FJFLX0 = (ZTAU(LZ+1)-ZTAU(LZ))/(ZTAU(LZ+1)-ZTAU(LZ-1))
         FJFLX(L) = 4.d0*(FJ(LZ-1)*FJFLX0 + FJ(LZ+1)*(1.d0-FJFLX0))
      enddo

c---NB if one needs the mean intensity throughout layer L (instead of mid-pt)
c---   then average (tau-weighted) the odd-points from: NELEV(L-1) to NELEV(L)
c---NB This is NOT now used.  Errors in cloudy layers are 0-10%  (low sun)
c---   the results are stored locally in FJACT2 (actinic)
       LZ1 = ND
      do L = 1,L_
       LZMID = ND+2-2*JNDLEV(L)
       LZ0 = ND+2-2*JNELEV(L)
         SUMT = 0.d0
         SUMJ = 0.d0
       do L2 = LZ0,LZ1-2,2
         SUMT = SUMT + ZTAU(L2+2)-ZTAU(L2)
         SUMJ = SUMJ + (ZTAU(L2+2)-ZTAU(L2))*(FJ(L2)+FJ(L2+2))*0.5
       enddo
         FJACT2(L) = 4.d0*SUMJ/SUMT + FZ(LZMID)
c-----LZ are indices to the 1:ND array     LZ1 > LZMID > LZ0
c       write(6,'(4i5,3f10.3)') L, LZ0,LZMID,LZ1,
c     &      ZTAU(LZMID),FJACT(L),FJACT2(L)
         LZ1 = LZ0
      enddo

c---diffuse fluxes reflected at top, incident at bottom 
        FJTOP = FJT
        FJBOT = FJB
c     call cpu_time(finish)
c      write(6,'(a,f6.3)'),'OPMIEtime=',finish-start
      return
      END  SUBROUTINE OPMIE



c<<<<<<<<<<<<<<<<<<<<<<<end core fast-J subroutines<<<<<<<<<<<<<<<<<<<<<


c<<<<<<<<<<<<<<<<<<<<<<<begin core scattering subroutines<<<<<<<<<<<<<<<
C-----------------------------------------------------------------------
      SUBROUTINE MIESCT(FJ,FJT,FJB,
     &                    POMEGA,FZ,ZTAU,ZFLUX,ZREFL,ZU0,MFIT,ND)
C-----------------------------------------------------------------------
C   This is an adaption of the Prather radiative transfer code, (mjp, 10/95)
C     Prather, 1974, Astrophys. J. 192, 787-792.
C         Sol_n of inhomogeneous Rayleigh scattering atmosphere. 
C         (original Rayleigh w/ polarization)
C     Cochran and Trafton, 1978, Ap.J., 219, 756-762.
C         Raman scattering in the atmospheres of the major planets.
C         (first use of anisotropic code)
C     Jacob, Gottlieb and Prather, 1989, J.Geophys.Res., 94, 12975-13002.
C         Chemistry of a polluted cloudy boundary layer,
C         (documentation of extension to anisotropic scattering)
C
C    takes atmospheric structure and source terms from std J-code
C    ALSO limited to 4 Gauss points, only calculates mean field!
C
C   mean rad. field ONLY (M=1)
C-----------------------------------------------------------------------
      USE FJ2008_MOD
      implicit none
c$$$      include 'parm_MIE.f'
c--- expect parameters M_, N_ in parm_MIE.f------------------------------
c
      integer, intent(in) ::  MFIT, ND
      real*8, intent(in)  ::  POMEGA(2*M_,N_),FZ(N_),ZTAU(N_)
     &                       ,ZREFL,ZU0,ZFLUX
      real*8, intent(out) ::  FJ(N_),FJT,FJB
c
      real*8  WT(M_),EMU(M_),PM(M_,2*M_),PM0(2*M_),CMEQ1
      integer I, ID, IM, M, N
C-----------------------------------------------------------------------
C---fix scattering to 4 Gauss pts = 8-stream
      call GAUSSP (N,EMU,WT)
c---calc in OPMIE:  ZFLUX = (ZU0*FZ(ND)*ZREFL)/(1.0d0+ZREFL)
      M = 1
      do I = 1,N
        call LEGND0 (EMU(I),PM0,MFIT)
        do IM = M,MFIT
          PM(I,IM) = PM0(IM)
        enddo
      enddo
C
      CMEQ1 = 0.25d0
      call LEGND0 (-ZU0,PM0,MFIT)
      do IM=M,MFIT
        PM0(IM) = CMEQ1*PM0(IM)
      enddo
C
      call BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0
     &            ,FJT,FJB,  M,N,MFIT,ND)
C
c      do ID=1,ND,2             !!!! this is now done in MIESCT()
c        FJ(ID) = 4.0d0*FJ(ID) + FZ(ID)
c      enddo

      return
      END SUBROUTINE MIESCT


C-----------------------------------------------------------------------
      SUBROUTINE BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0
     &                 ,FJTOP,FJBOT,  M,N,MFIT,ND)
C-----------------------------------------------------------------------
C  Sets up and solves the block tri-diagonal system:  
C               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
c  This goes back to the old, dumb, fast version 5.3
C-----------------------------------------------------------------------
      USE FJ2008_MOD
      USE ERROR_MOD
      implicit none
c$$$      include 'parm_MIE.f'
c--- expect parameters M_, N_ in parm_MIE.f------------------------------

      integer, intent(in) ::  M, N, MFIT, ND
      real*8, intent(in)  ::  POMEGA(2*M_,N_),FZ(N_),ZTAU(N_)
     &                       ,WT(M_),EMU(M_),PM(M_,2*M_),PM0(2*M_)
     &                       ,ZFLUX,ZREFL
      real*8, intent(out) ::  FJ(N_),FJTOP,FJBOT
c
      real*8, dimension(M_)    :: A, C1, H
      real*8, dimension(M_,M_) :: B, AA, CC
      real*8                      DD(M_,M_,N_), RR(M_,N_)
      real*8  SUMM, FIPLUS
      integer I, J, K, ID
c
C-----------UPPER BOUNDARY ID=1
      call GEN(POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0
     &        ,B,CC,AA,A,H,C1,M,N,MFIT,ND,1)
      call MATIN4 (B)
      do I = 1,N
         RR(I,1) = 0.0d0
        do J = 1,N
          SUMM = 0.0d0
         do K = 1,N
          SUMM = SUMM - B(I,K)*CC(K,J)
         enddo
         DD(I,J,1) = SUMM
         RR(I,1) = RR(I,1) + B(I,J)*H(J)
        enddo
      enddo
C----------CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1
      do ID = 2,ND-1
        call GEN(POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0
     &          ,B,CC,AA,A,H,C1,M,N,MFIT,ND,ID)
        do I = 1,N
          do J = 1,N
          B(I,J) = B(I,J) + A(I)*DD(I,J,ID-1)
          enddo
          H(I) = H(I) - A(I)*RR(I,ID-1)
        enddo
        call MATIN4 (B)
        do I = 1,N
          RR(I,ID) = 0.0d0
          do J = 1,N
          RR(I,ID) = RR(I,ID) + B(I,J)*H(J)
          DD(I,J,ID) = - B(I,J)*C1(J)
          enddo
        enddo
      enddo
C---------FINAL DEPTH POINT: ND
      call GEN(POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0
     &        ,B,CC,AA,A,H,C1,M,N,MFIT,ND,ND)
      do I = 1,N
        do J = 1,N
          SUMM = 0.0d0
          do K = 1,N
          SUMM = SUMM + AA(I,K)*DD(K,J,ND-1)
          enddo
        B(I,J) = B(I,J) + SUMM
        H(I) = H(I) - AA(I,J)*RR(J,ND-1)
        enddo
      enddo
      call MATIN4 (B)
      do I = 1,N
        RR(I,ND) = 0.0d0
        do J = 1,N
        RR(I,ND) = RR(I,ND) + B(I,J)*H(J)
        enddo
      enddo
C-----------BACK SOLUTION
      do ID = ND-1,1,-1
       do I = 1,N
        do J = 1,N
         RR(I,ID) = RR(I,ID) + DD(I,J,ID)*RR(J,ID+1)
        enddo
       enddo
      enddo

C----------MEAN J & H
      do ID = 1,ND,2
        FJ(ID) = 0.0d0
       do I = 1,N
        FJ(ID) = FJ(ID) + RR(I,ID)*WT(I)
       enddo
      enddo
      do ID = 2,ND,2
        FJ(ID) = 0.0d0
       do I = 1,N
        FJ(ID) = FJ(ID) + RR(I,ID)*WT(I)*EMU(I)
       enddo
      enddo

c---FJTOP = scaled diffuse flux out top-of-atmosphere (limit = mu0)
c---FJBOT = scaled diffuse flux onto surface:  
        FJTOP = 0.0d0
        FJBOT = 0.0d0
       do I = 1,N
        FJTOP = FJTOP + RR(I,1)*WT(I)*EMU(I)
        FJBOT = FJBOT + RR(I,ND)*WT(I)*EMU(I)
       enddo
        FJTOP = 4.d0*FJTOP
c---         ZFLUX = reflect/(1 + reflect) * mu0 * Fsolar(lower boundary)
        FIPLUS = 4.d0*ZREFL*FJBOT/(1.0d0 + ZREFL) + ZFLUX
        FJBOT = 4.d0*FJBOT - FIPLUS

      return
      END SUBROUTINE BLKSLV


C-----------------------------------------------------------------------
      SUBROUTINE GEN(POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0
     &              ,B,CC,AA,A,H,C1,M,N,MFIT,ND,ID)
C-----------------------------------------------------------------------
C  Generates coefficient matrices for the block tri-diagonal system:
C               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
C-----------------------------------------------------------------------
      USE FJ2008_MOD
      USE ERROR_MOD
      USE LOGICAL_MOD
      implicit none
c$$$      include 'parm_MIE.f'
c--- expect parameters M_, N_ in parm_MIE.f------------------------------
      integer, intent(in) ::  M, N, MFIT, ND, ID
      real*8, intent(in)  ::  POMEGA(2*M_,N_),FZ(N_),ZTAU(N_)
     &                       ,WT(M_),EMU(M_),PM(M_,2*M_),PM0(2*M_)
     &                       ,ZFLUX,ZREFL
      real*8, intent(out) ::  B(M_,M_),AA(M_,M_),CC(M_,M_),A(M_),C1(M_)
     &                       ,H(M_)
c
      integer ID0, ID1, IM, I, J, K, MSTART
      real*8  SUM0, SUM1, SUM2, SUM3
      real*8  DELTAU, D1, D2, SURFAC
c
      real*8  S(M_,M_), W(M_,M_), U1(M_,M_), V1(M_)
C---------------------------------------------
      if (ID.eq.1 .or. ID.eq.ND) then
C---------calculate generic 2nd-order terms for boundaries
       ID0 = ID
       ID1 = ID+1
       if (ID .ge. ND) ID1 = ID-1

       do I = 1,N
          SUM0 = 0.0d0
          SUM1 = 0.0d0
          SUM2 = 0.0d0
          SUM3 = 0.0d0
        do IM = M,MFIT,2
          SUM0 = SUM0 + POMEGA(IM,ID0)*PM(I,IM)*PM0(IM)
          SUM2 = SUM2 + POMEGA(IM,ID1)*PM(I,IM)*PM0(IM)
        enddo
        do IM = M+1,MFIT,2
          SUM1 = SUM1 + POMEGA(IM,ID0)*PM(I,IM)*PM0(IM)
          SUM3 = SUM3 + POMEGA(IM,ID1)*PM(I,IM)*PM0(IM)
        enddo
         H(I) = 0.5d0*(SUM0*FZ(ID0) + SUM2*FZ(ID1))
         A(I) = 0.5d0*(SUM1*FZ(ID0) + SUM3*FZ(ID1))
        do J = 1,I
          SUM0 = 0.0d0
          SUM1 = 0.0d0
          SUM2 = 0.0d0
          SUM3 = 0.0d0
         do IM = M,MFIT,2
          SUM0 = SUM0 + POMEGA(IM,ID0)*PM(I,IM)*PM(J,IM)
          SUM2 = SUM2 + POMEGA(IM,ID1)*PM(I,IM)*PM(J,IM)
         enddo
         do IM = M+1,MFIT,2
          SUM1 = SUM1 + POMEGA(IM,ID0)*PM(I,IM)*PM(J,IM)
          SUM3 = SUM3 + POMEGA(IM,ID1)*PM(I,IM)*PM(J,IM)
         enddo
         S(I,J) = - SUM2*WT(J)
         S(J,I) = - SUM2*WT(I)
         W(I,J) = - SUM1*WT(J)
         W(J,I) = - SUM1*WT(I)
         U1(I,J) = - SUM3*WT(J)
         U1(J,I) = - SUM3*WT(I)
          SUM0 = 0.5d0*(SUM0 + SUM2)
         B(I,J) = - SUM0*WT(J)
         B(J,I) = - SUM0*WT(I)
        enddo
         S(I,I) = S(I,I) + 1.0d0
         W(I,I) = W(I,I) + 1.0d0
         U1(I,I) = U1(I,I) + 1.0d0
         B(I,I) = B(I,I) + 1.0d0
       enddo

       do I = 1,N
         SUM0 = 0.0d0
        do J = 1,N
         SUM0 = SUM0 + S(I,J)*A(J)/EMU(J)
        enddo
        C1(I) = SUM0
       enddo
       do I = 1,N
        do J = 1,N
          SUM0 = 0.0d0
          SUM2 = 0.0d0
         do K = 1,N
          SUM0 = SUM0 + S(J,K)*W(K,I)/EMU(K)
          SUM2 = SUM2 + S(J,K)*U1(K,I)/EMU(K)
         enddo
         A(J) = SUM0
         V1(J) = SUM2
        enddo
        do J = 1,N
         W(J,I) = A(J)
         U1(J,I) = V1(J)
        enddo
       enddo
       if (ID .eq. 1) then
C-------------upper boundary, 2nd-order, C-matrix is full (CC)
        DELTAU = ZTAU(2) - ZTAU(1)
        D2 = 0.25d0*DELTAU
        do I = 1,N
          D1 = EMU(I)/DELTAU
          do J = 1,N
           B(I,J) = B(I,J) + D2*W(I,J)
           CC(I,J) = D2*U1(I,J)
          enddo
          B(I,I) = B(I,I) + D1
          CC(I,I) = CC(I,I) - D1
C         H(I) = H(I) + 2.0d0*D2*C1(I) + D1*SISOTP
          H(I) = H(I) + 2.0d0*D2*C1(I)
          A(I) = 0.0d0
        enddo
       else
C-------------lower boundary, 2nd-order, A-matrix is full (AA)
        DELTAU = ZTAU(ND) - ZTAU(ND-1)
        D2 = 0.25d0*DELTAU
        SURFAC = 4.0d0*ZREFL/(1.0d0 + ZREFL)
        do I = 1,N
          D1 = EMU(I)/DELTAU
          H(I) = H(I) - 2.0d0*D2*C1(I)
           SUM0 = 0.0d0
          do J = 1,N
           SUM0 = SUM0 + W(I,J)
          enddo
           SUM0 = D1 + D2*SUM0
           SUM1 = SURFAC*SUM0
          do J = 1,N
           B(I,J) = B(I,J) + D2*W(I,J) - SUM1*EMU(J)*WT(J)
          enddo
          B(I,I) = B(I,I) + D1
          H(I) = H(I) + SUM0*ZFLUX
          do J = 1,N
           AA(I,J) = - D2*U1(I,J)
          enddo
           AA(I,I) = AA(I,I) + D1
           C1(I) = 0.0d0
        enddo
       endif
C------------intermediate points:  can be even or odd, A & C diagonal
      else
        DELTAU = ZTAU(ID+1) - ZTAU(ID-1)
        MSTART = M + MOD(ID+1,2)
        do I = 1,N
          A(I) = EMU(I)/DELTAU
          C1(I) = -A(I)
           SUM0 = 0.0d0
          do IM = MSTART,MFIT,2
           SUM0 = SUM0 + POMEGA(IM,ID)*PM(I,IM)*PM0(IM)
          enddo
          H(I) = SUM0*FZ(ID)
          do J=1,I
            SUM0 = 0.0d0
           do IM = MSTART,MFIT,2
            SUM0 = SUM0 + POMEGA(IM,ID)*PM(I,IM)*PM(J,IM)
           enddo
            B(I,J) =  - SUM0*WT(J)
            B(J,I) =  - SUM0*WT(I)
          enddo
          B(I,I) = B(I,I) + 1.0d0
        enddo
      endif
      return
      END SUBROUTINE GEN


C-----------------------------------------------------------------------
      SUBROUTINE LEGND0 (X,PL,N)
C-----------------------------------------------------------------------
C---Calculates ORDINARY Legendre fns of X (real) 
C---   from P[0] = PL(1) = 1,  P[1] = X, .... P[N-1] = PL(N)
      implicit none
      integer, intent(in) :: N
      real*8, intent(in)  :: X
      real*8, intent(out) :: PL(N)
      integer I
      real*8  DEN
C---Always does PL(2) = P[1]
        PL(1) = 1.d0
        PL(2) = X
        do I = 3,N
         DEN = (I-1)
         PL(I) = PL(I-1)*X*(2.d0-1.0/DEN) - PL(I-2)*(1.d0-1.d0/DEN)
        enddo
      return
      END SUBROUTINE LEGND0


C-----------------------------------------------------------------------
      subroutine MATIN4 (A)
C-----------------------------------------------------------------------
C  invert 4x4 matrix A(4,4) in place with L-U decomposition (mjp, old...)
      USE ERROR_MOD
      USE LOGICAL_MOD
      implicit none
      real*8, intent(inout)  ::  A(4,4)
      real*8 :: B(4,4)

      B=A
C---SETUP L AND U
      A(2,1) = A(2,1)/A(1,1)
      A(2,2) = A(2,2)-A(2,1)*A(1,2)
      A(2,3) = A(2,3)-A(2,1)*A(1,3)
      A(2,4) = A(2,4)-A(2,1)*A(1,4)
      A(3,1) = A(3,1)/A(1,1)
      A(3,2) = (A(3,2)-A(3,1)*A(1,2))/A(2,2)
      A(3,3) = A(3,3)-A(3,1)*A(1,3)-A(3,2)*A(2,3)
      A(3,4) = A(3,4)-A(3,1)*A(1,4)-A(3,2)*A(2,4)
      A(4,1) = A(4,1)/A(1,1)
      A(4,2) = (A(4,2)-A(4,1)*A(1,2))/A(2,2)
      A(4,3) = (A(4,3)-A(4,1)*A(1,3)-A(4,2)*A(2,3))/A(3,3)
      A(4,4) = A(4,4)-A(4,1)*A(1,4)-A(4,2)*A(2,4)-A(4,3)*A(3,4)
C---INVERT L
      A(4,3) = -A(4,3)
      A(4,2) = -A(4,2)-A(4,3)*A(3,2)
      A(4,1) = -A(4,1)-A(4,2)*A(2,1)-A(4,3)*A(3,1)
      A(3,2) = -A(3,2)
      A(3,1) = -A(3,1)-A(3,2)*A(2,1)
      A(2,1) = -A(2,1)
C---INVERT U
      A(4,4) = 1.d0/A(4,4)
      A(3,4) = -A(3,4)*A(4,4)/A(3,3)
      A(3,3) = 1.d0/A(3,3)
      A(2,4) = -(A(2,3)*A(3,4)+A(2,4)*A(4,4))/A(2,2)
      A(2,3) = -A(2,3)*A(3,3)/A(2,2)
      A(2,2) = 1.d0/A(2,2)
      A(1,4) = -(A(1,2)*A(2,4)+A(1,3)*A(3,4)+A(1,4)*A(4,4))/A(1,1)
      A(1,3) = -(A(1,2)*A(2,3)+A(1,3)*A(3,3))/A(1,1)
      A(1,2) = -A(1,2)*A(2,2)/A(1,1)
      A(1,1) = 1.d0/A(1,1)
C---MULTIPLY (U-INVERSE)*(L-INVERSE)
      A(1,1) = A(1,1)+A(1,2)*A(2,1)+A(1,3)*A(3,1)+A(1,4)*A(4,1)
      A(1,2) = A(1,2)+A(1,3)*A(3,2)+A(1,4)*A(4,2)
      A(1,3) = A(1,3)+A(1,4)*A(4,3)
      A(2,1) = A(2,2)*A(2,1)+A(2,3)*A(3,1)+A(2,4)*A(4,1)
      A(2,2) = A(2,2)+A(2,3)*A(3,2)+A(2,4)*A(4,2)
      A(2,3) = A(2,3)+A(2,4)*A(4,3)
      A(3,1) = A(3,3)*A(3,1)+A(3,4)*A(4,1)
      A(3,2) = A(3,3)*A(3,2)+A(3,4)*A(4,2)
      A(3,3) = A(3,3)+A(3,4)*A(4,3)
      A(4,1) = A(4,4)*A(4,1)
      A(4,2) = A(4,4)*A(4,2)
      A(4,3) = A(4,4)*A(4,3)
      return
      END SUBROUTINE MATIN4


C-----------------------------------------------------------------------
      SUBROUTINE GAUSSP (N,XPT,XWT)
C-----------------------------------------------------------------------
C  Loads in pre-set Gauss points for 4 angles from 0 to +1 in cos(theta)=mu
      implicit none
      integer, intent(out) :: N
      real*8, intent(out)  ::  XPT(*),XWT(*)
      real*8   GPT4(4),GWT4(4)
      integer  I
      data GPT4/.06943184420297d0,.33000947820757d0,.66999052179243d0,
     &          .93056815579703d0/
      data GWT4/.17392742256873d0,.32607257743127d0,.32607257743127d0,
     &          .17392742256873d0/
      N = 4
      do I = 1,N
        XPT(I) = GPT4(I)
        XWT(I) = GWT4(I)
      enddo
      return
      END SUBROUTINE GAUSSP
c<<<<<<<<<<<<<<<<<<<<<<<<<end core scattering subroutines<<<<<<<<<<<<<<<
      END MODULE FAST_JX_MOD
