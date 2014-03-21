!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: 
!
! !DESCRIPTION: Read in wavelength bins, solar fluxes, Rayleigh parameters, 
!  temperature-dependent cross sections and Rayleigh/aerosol scattering phase 
!  functions with temperature dependences.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RD_TJPL( NJ1, NAMFIL, am_I_Root )
!
! !USES:
!
      USE CMN_FJ_MOD, ONLY : JPMAX, JPPJ
      USE JV_CMN_MOD

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      INTEGER,          INTENT(IN) :: NJ1         ! Logical unit #
      CHARACTER(LEN=*), INTENT(IN) :: NAMFIL      ! File name
      LOGICAL,          INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! !INPUT PARAMETERS: 
!
!
! !REMARKS:
!     NAMFIL   Name of spectral data file (jv_spec.dat)
!     NJ1      Channel number for reading data file
!     NJVAL    Number of species to calculate J-values for
!     NWWW     Number of wavelength bins, from NW1:NW2
!     WBIN     Boundaries of wavelength bins
!     WL       Centres of wavelength bins - 'effective wavelength'
!     FL       Solar flux incident on top of atmosphere (cm-2.s-1)
!     QRAYL    Rayleigh parameters (effective cross-section) (cm2)
!     QBC      Black Carbon absorption extinct. (specific cross-sect.) (m2/g)
!     QO2      O2 cross-sections
!     QO3      O3 cross-sections
!     Q1D      O3 => O(1D) quantum yield
!     TQQ      Temperature for supplied cross sections
!     QQQ      Supplied cross sections in each wavelength bin (cm2)
!     NAA      Number of categories for scattering phase functions
!     QAA      Aerosol scattering phase functions
!     NK       Number of wavelengths at which functions supplied (set as 4)
!     WAA      Wavelengths for the NK supplied phase functions
!     PAA      Phase function: first 8 terms of expansion
!     RAA      Effective radius associated with aerosol type
!     SSA      Single scattering albedo
!     npdep    Number of pressure dependencies
!     zpdep    Pressure dependencies by wavelength bin
!     jpdep    Index of cross sections requiring pressure dependence
!     lpdep    Label for pressure dependence
!
! !REVISION HISTORY: 
!  01 Jun 1998 - P. Murti - Initial version
!  (1 ) Updated to include new pressure-dependancy function for GLYX and MGLY. 
!        (tmf, 1/7/09)
!  (2 ) Added a pressure-dependancy function selector 'pdepf'. 
!        (tmf, ccc, 1/7/09)
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  30 Jul 2012 - R. Yantosca - AddedoTeX headers
!  01 Aug 2012 - R. Yantosca - Now restore NJ1 to INTENT(IN) status
!  01 Mar 2013 - R. Yantosca - Block some extra print statements w/ am_I_Root
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
      integer i, j, k, iw, nk, nqqq, nwww
      character*7  lpdep(7)
      do J=1,NS
        do K=1,3
          TQQ(K,J) = 0.d0
        enddo
      enddo
!-------------spectral data---------------------------------------------
      open(NJ1, FILE=NAMFIL)
      read(NJ1,'(A)') TITLE0
      IF ( am_I_Root ) THEN
         write(6,'(1X,A)') TITLE0
      ENDIF
      read(NJ1,'(10X,14I5)') NJVAL,NWWW,NW1,NW2
      if(NJVAL.gt.NS) then
        IF ( am_I_Root ) write(6,300) NJVAL,NS
        stop
      endif
!------------NQQQ = no. additional J-values from X-sects (O2,O3P,O3D+NQQQ)
!- NQQQ is changed to NJVAL-1 because there are 2 dummy species at the end
! used for acetone pressure dependency only. (ccc, 4/20/09)
      NQQQ = NJVAL-1
      read(NJ1,102) (WBIN(IW),IW=1,NWWW)
      read(NJ1,102) (WBIN(IW+1),IW=1,NWWW)
      read(NJ1,102) (WL(IW),IW=1,NWWW)
      read(NJ1,102) (FL(IW),IW=1,NWWW)
      read(NJ1,102) (QRAYL(IW),IW=1,NWWW)
      read(NJ1,102) (QBC(IW),IW=1,NWWW)   !  From Loiusse et al. [JGR, 1996]
c
!---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
      do K=1,3
        read(NJ1,103) TITLEJ(K,1),TQQ(K,1), (QO2(IW,K),IW=1,NWWW)
      enddo
      do K=1,3
        read(NJ1,103) TITLEJ(K,2),TQQ(K,2), (QO3(IW,K),IW=1,NWWW)
      enddo
      do K=1,3
        read(NJ1,103) TITLEJ(K,3),TQQ(K,3), (Q1D(IW,K),IW=1,NWWW)
      enddo
      IF ( am_I_Root ) THEN
         do K=1,3
            write(6,200) titlej(1,k),(tqq(i,k),i=1,3)
         enddo
      ENDIF
c
!---Read remaining species:  X-sections at 2 T's
      do J=1,NQQQ
        read(NJ1,103) TITLEJ(1,J+3),TQQ(1,J+3),(QQQ(IW,1,J),IW=1,NWWW)
        read(NJ1,103) TITLEJ(2,J+3),TQQ(2,J+3),(QQQ(IW,2,J),IW=1,NWWW)
        IF ( am_I_Root ) THEN 
           write(6,200) titlej(1,j+3),(tqq(i,j+3),i=1,2)
        ENDIF
      enddo
      read(NJ1,'(A)') TITLE0
c
c---Pressure dependencies
      read(NJ1,104) npdep
      do k=1,npdep
         read(NJ1,105) lpdep(k), pdepf(k), (zpdep(iw,k),iw=1,nwww)
         IF ( am_I_Root ) THEN 
            write(6,201)  lpdep(k), pdepf(k), (zpdep(iw,k),iw=1,nwww)
         ENDIF

         !--------------------------------------
         ! Special treatment for MGLY pressure dependency
         ! (tmf, 11/16/06)
         !--------------------------------------
         if ( pdepf(k) .eq. 4 ) then           
            ! pass zpdep to mglypdep
            mglypdep(:,1) = zpdep(:,k)
            read(NJ1,105) lpdep(k), pdepf(k), (mglypdep(iw,2),iw=1,nwww)
            read(NJ1,105) lpdep(k), pdepf(k), (mglypdep(iw,3),iw=1,nwww)
         endif
      enddo
      read(NJ1,'(A)') TITLE0

c
c---Zero index arrays
      do j=1,jppj
        jind(j)=0
      enddo
      do j=1,NJVAL
        jpdep(j)=0
      enddo
c
!---Set mapping index
      do j=1,NJVAL
        do k=1,jppj
          if (jlabel(k).eq.titlej(1,j)) jind(k)=j
        enddo
        do k=1,npdep
          if (lpdep(k).eq.titlej(1,j)) jpdep(j)=k
        enddo
      enddo
      do k=1,jppj
        if(jfacta(k).eq.0.d0 .and. am_I_Root )
     &             write(6,*) 'Not using photolysis reaction ',k
        if(jind(k).eq.0) then
          if(jfacta(k).eq.0.d0) then
            jind(k)=1
          else
            IF ( am_I_Root ) THEN
             write(6,*) 'Which J-rate for photolysis reaction ',k,' ?'
             stop
            ENDIF
          endif
        endif
      enddo
c
!---Read aerosol phase functions:
      read(NJ1,'(A10,I5,/)') TITLE0,NAA
      NK=4        ! Fix number of wavelengths at 4
      do j=1,NAA
        read(NJ1,110) TITLEA(j)
        do k=1,NK
          read(NJ1,*) WAA(k,j),QAA(k,j),RAA(k,j),SSA(k,j),
     &                                             (PAA(i,k,j),i=1,8)
        enddo
      enddo
c
      IF ( am_I_Root ) THEN 
         write(6,*) 'Aerosol phase functions & wavelengths'
         do J=1,NAA
            write(6,'(1x,A8,I2,A,9F8.1)')
     &           TITLEA(J),J,'  wavel=',(WAA(K,J),K=1,NK)
            write(6,'(9x,I2,A,9F8.4)') J,'  Qext =',(QAA(K,J),K=1,NK)
         enddo
      ENDIF
!--------
! Modify reading and writing formats 105 & 201 for pressure dependancy 
c (ccc, 1/7/09)

  101 FORMAT(8E10.3)
  102 FORMAT(10X,7E10.3)
  103 FORMAT(A7,F3.0,7E10.3)
c 103 FORMAT(A7,F3.0,7E10.3/(10X,7E10.3))
  104 FORMAT(13x,i2)
  105 FORMAT(A7,2x,I1,7E10.3)
  110 format(3x,a20)
  200 format(1x,' x-sect:',a10,3(3x,f6.2))
  201 format(1x,' pr.dep:',a10,1x,I1,7(1pE10.3))
  300 format(' Number of x-sections supplied to Fast-J: ',i3,/,
     &       ' Maximum number allowed (NS) only set to: ',i3,
     &       ' - increase in jv_cmn.h')
      close(NJ1)
      return
      END SUBROUTINE RD_TJPL
!EOC
