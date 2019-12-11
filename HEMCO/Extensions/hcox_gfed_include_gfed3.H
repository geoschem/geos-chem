!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hcox_gfed_include_gfed3.H
!
! !DESCRIPTION: Include file with GFED emission factor data that was
!  originally contained in file GFED\_emission\_factors.txt.  We have now
!  transformed this file into hardwired F90 commands in order to avoid reading
!  an ASCII file in the ESMF environment.
!
! !REMARKS:
!  ABOUT THIS FILE:
!  ----------------
!  This file was created by script HEMCO/Extensions/Preprocess/gfed.pl.
!  This script can be executed with the following command:
!
!    cd HEMCO/Extensions/Preprocess
!    make gfed
!
!  This will regenerate this include file from the original data and
!  automatically place it in the HEMCO/Extensions directory.
!
!  White space has been removed in order to reduce the file size as much
!  as possible.  If you have to recreate this file, then it is easier to
!  generate via the Perl script than to try to hand edit the code below.
!
!  DATA:
!  -----
!  The GFED3_EMFAC array contains emission factors in kg/kgDM or kgC/kgDM
!  GFED3_EMFAC(N,1) = Agricultural Waste   Emission Factor for species N
!  GFED3_EMFAC(N,1) = Deforestation        Emission Factor for species N
!  GFED3_EMFAC(N,1) = Extratropical Forest Emission Factor for species N
!  GFED3_EMFAC(N,1) = Peat                 Emission Factor for species N
!  GFED3_EMFAC(N,1) = Savanna              Emission Factor for species N
!  GFED3_EMFAC(N,1) = Woodland             Emission Factor for species N
!
!  References:
!  ============================================================================
!  (1 ) Akagi, S.K., R.J. Yokelson, C. Wiedinmyer, M.J. Alvarado, J.S. Reid,
!       T. Karl, J.D. Crounse, and P.O. Wennberg, Emission factors for open and
!       domestic biomass burning for use in atmospheric models, Atmos. Chem.
!       Phys., 11, 4039-4072, 10.5194/acp-11-4039-2011, 2011.
!
! !REVISION HISTORY:
!  08 Aug 2014 - R. Yantosca - Initial version
!  10 Jul 2015 - R. Yantosca - Fixed minor issues in ProTeX header
!  16 Oct 2015 - M. Sulprizio- Added the following updates from E. Fischer:
!                              - Added EOH and MTPA for PAN updates;
!                              - Updated the emission factors for NMVOCs and
!                                NOx from extratropical forests, savannas and
!                                agricultural fires from Akagi et al. (2011,ACP)
!EOP
!------------------------------------------------------------------------------
!BOC

! NO
GFED3_SPEC_NAME(1)="NO"
GFED3_EMFAC(1,1)=2.29E-03_hp
GFED3_EMFAC(1,2)=2.26E-03_hp
GFED3_EMFAC(1,3)=1.12E-03_hp
GFED3_EMFAC(1,4)=2.26E-03_hp
GFED3_EMFAC(1,5)=2.12E-03_hp
GFED3_EMFAC(1,6)=2.19E-03_hp

! CO
GFED3_SPEC_NAME(2)="CO"
GFED3_EMFAC(2,1)=9.22E-02_hp
GFED3_EMFAC(2,2)=1.01E-01_hp
GFED3_EMFAC(2,3)=1.06E-01_hp
GFED3_EMFAC(2,4)=2.10E-01_hp
GFED3_EMFAC(2,5)=6.25E-02_hp
GFED3_EMFAC(2,6)=8.17E-02_hp

! ALK4
GFED3_SPEC_NAME(3)="ALK4"
GFED3_EMFAC(3,1)=1.19E-04_hp
GFED3_EMFAC(3,2)=1.92E-04_hp
GFED3_EMFAC(3,3)=3.12E-04_hp
GFED3_EMFAC(3,4)=1.92E-04_hp
GFED3_EMFAC(3,5)=1.06E-04_hp
GFED3_EMFAC(3,6)=1.53E-04_hp

! ACET
GFED3_SPEC_NAME(4)="ACET"
GFED3_EMFAC(4,1)=2.79E-04_hp
GFED3_EMFAC(4,2)=3.91E-04_hp
GFED3_EMFAC(4,3)=4.65E-04_hp
GFED3_EMFAC(4,4)=3.91E-04_hp
GFED3_EMFAC(4,5)=9.92E-05_hp
GFED3_EMFAC(4,6)=3.47E-04_hp

! MEK
GFED3_SPEC_NAME(5)="MEK"
GFED3_EMFAC(5,1)=6.89E-04_hp
GFED3_EMFAC(5,2)=3.91E-04_hp
GFED3_EMFAC(5,3)=2.31E-04_hp
GFED3_EMFAC(5,4)=3.91E-04_hp
GFED3_EMFAC(5,5)=1.93E-04_hp
GFED3_EMFAC(5,6)=2.93E-04_hp

! ALD2
GFED3_SPEC_NAME(6)="ALD2"
GFED3_EMFAC(6,1)=6.76E-04_hp
GFED3_EMFAC(6,2)=1.23E-03_hp
GFED3_EMFAC(6,3)=5.34E-04_hp
GFED3_EMFAC(6,4)=1.23E-03_hp
GFED3_EMFAC(6,5)=3.11E-04_hp
GFED3_EMFAC(6,6)=7.52E-04_hp

! PRPE
GFED3_SPEC_NAME(7)="PRPE"
GFED3_EMFAC(7,1)=1.02E-03_hp
GFED3_EMFAC(7,2)=6.54E-04_hp
GFED3_EMFAC(7,3)=1.36E-03_hp
GFED3_EMFAC(7,4)=6.54E-04_hp
GFED3_EMFAC(7,5)=7.45E-04_hp
GFED3_EMFAC(7,6)=4.85E-04_hp

! C3H8
GFED3_SPEC_NAME(8)="C3H8"
GFED3_EMFAC(8,1)=2.29E-04_hp
GFED3_EMFAC(8,2)=4.39E-04_hp
GFED3_EMFAC(8,3)=3.43E-04_hp
GFED3_EMFAC(8,4)=4.39E-04_hp
GFED3_EMFAC(8,5)=8.16E-05_hp
GFED3_EMFAC(8,6)=2.59E-04_hp

! CH2O
GFED3_SPEC_NAME(9)="CH2O"
GFED3_EMFAC(9,1)=2.08E-03_hp
GFED3_EMFAC(9,2)=2.22E-03_hp
GFED3_EMFAC(9,3)=1.92E-03_hp
GFED3_EMFAC(9,4)=2.22E-03_hp
GFED3_EMFAC(9,5)=8.66E-04_hp
GFED3_EMFAC(9,6)=1.54E-03_hp

! C2H6
GFED3_SPEC_NAME(10)="C2H6"
GFED3_EMFAC(10,1)=9.91E-04_hp
GFED3_EMFAC(10,2)=7.05E-04_hp
GFED3_EMFAC(10,3)=1.36E-03_hp
GFED3_EMFAC(10,4)=7.05E-04_hp
GFED3_EMFAC(10,5)=3.12E-04_hp
GFED3_EMFAC(10,6)=5.09E-04_hp

! SO2
GFED3_SPEC_NAME(11)="SO2"
GFED3_EMFAC(11,1)=4.00E-04_hp
GFED3_EMFAC(11,2)=7.59E-04_hp
GFED3_EMFAC(11,3)=9.99E-04_hp
GFED3_EMFAC(11,4)=7.59E-04_hp
GFED3_EMFAC(11,5)=3.37E-04_hp
GFED3_EMFAC(11,6)=5.48E-04_hp

! NH3
GFED3_SPEC_NAME(12)="NH3"
GFED3_EMFAC(12,1)=1.62E-03_hp
GFED3_EMFAC(12,2)=9.36E-04_hp
GFED3_EMFAC(12,3)=1.55E-03_hp
GFED3_EMFAC(12,4)=9.36E-04_hp
GFED3_EMFAC(12,5)=8.49E-04_hp
GFED3_EMFAC(12,6)=8.93E-04_hp

! BC
GFED3_SPEC_NAME(13)="BC"
GFED3_EMFAC(13,1)=4.21E-04_hp
GFED3_EMFAC(13,2)=5.65E-04_hp
GFED3_EMFAC(13,3)=5.62E-04_hp
GFED3_EMFAC(13,4)=5.65E-04_hp
GFED3_EMFAC(13,5)=4.61E-04_hp
GFED3_EMFAC(13,6)=5.13E-04_hp

! OC
GFED3_SPEC_NAME(14)="OC"
GFED3_EMFAC(14,1)=4.15E-03_hp
GFED3_EMFAC(14,2)=4.30E-03_hp
GFED3_EMFAC(14,3)=9.14E-03_hp
GFED3_EMFAC(14,4)=4.30E-03_hp
GFED3_EMFAC(14,5)=3.21E-03_hp
GFED3_EMFAC(14,6)=3.76E-03_hp

! GLYX
GFED3_SPEC_NAME(15)="GLYX"
GFED3_EMFAC(15,1)=9.95E-04_hp
GFED3_EMFAC(15,2)=9.95E-04_hp
GFED3_EMFAC(15,3)=9.95E-04_hp
GFED3_EMFAC(15,4)=9.95E-04_hp
GFED3_EMFAC(15,5)=9.95E-04_hp
GFED3_EMFAC(15,6)=9.95E-04_hp

! MGLY
GFED3_SPEC_NAME(16)="MGLY"
GFED3_EMFAC(16,1)=7.30E-04_hp
GFED3_EMFAC(16,2)=7.30E-04_hp
GFED3_EMFAC(16,3)=7.30E-04_hp
GFED3_EMFAC(16,4)=7.30E-04_hp
GFED3_EMFAC(16,5)=7.30E-04_hp
GFED3_EMFAC(16,6)=7.30E-04_hp

! BENZ
GFED3_SPEC_NAME(17)="BENZ"
GFED3_EMFAC(17,1)=2.83E-04_hp
GFED3_EMFAC(17,2)=3.50E-04_hp
GFED3_EMFAC(17,3)=4.81E-04_hp
GFED3_EMFAC(17,4)=3.50E-04_hp
GFED3_EMFAC(17,5)=1.84E-04_hp
GFED3_EMFAC(17,6)=3.09E-04_hp

! TOLU
GFED3_SPEC_NAME(18)="TOLU"
GFED3_EMFAC(18,1)=1.73E-04_hp
GFED3_EMFAC(18,2)=2.10E-04_hp
GFED3_EMFAC(18,3)=4.84E-04_hp
GFED3_EMFAC(18,4)=2.10E-04_hp
GFED3_EMFAC(18,5)=7.84E-05_hp
GFED3_EMFAC(18,6)=1.86E-04_hp

! XYLE
GFED3_SPEC_NAME(19)="XYLE"
GFED3_EMFAC(19,1)=1.03E-04_hp
GFED3_EMFAC(19,2)=7.88E-05_hp
GFED3_EMFAC(19,3)=2.41E-04_hp
GFED3_EMFAC(19,4)=7.88E-05_hp
GFED3_EMFAC(19,5)=3.89E-05_hp
GFED3_EMFAC(19,6)=5.89E-05_hp

! C2H4
GFED3_SPEC_NAME(20)="C2H4"
GFED3_EMFAC(20,1)=1.14E-03_hp
GFED3_EMFAC(20,2)=9.80E-04_hp
GFED3_EMFAC(20,3)=9.40E-04_hp
GFED3_EMFAC(20,4)=9.80E-04_hp
GFED3_EMFAC(20,5)=7.41E-04_hp
GFED3_EMFAC(20,6)=8.61E-04_hp

! C2H2
GFED3_SPEC_NAME(21)="C2H2"
GFED3_EMFAC(21,1)=2.00E-04_hp
GFED3_EMFAC(21,2)=1.82E-04_hp
GFED3_EMFAC(21,3)=2.38E-04_hp
GFED3_EMFAC(21,4)=1.82E-04_hp
GFED3_EMFAC(21,5)=2.92E-04_hp
GFED3_EMFAC(21,6)=2.37E-04_hp

! GLYC
GFED3_SPEC_NAME(22)="GLYC"
GFED3_EMFAC(22,1)=0.00E+00_hp
GFED3_EMFAC(22,2)=9.57E-04_hp
GFED3_EMFAC(22,3)=9.57E-04_hp
GFED3_EMFAC(22,4)=9.57E-04_hp
GFED3_EMFAC(22,5)=9.57E-04_hp
GFED3_EMFAC(22,6)=9.57E-04_hp

! HAC
GFED3_SPEC_NAME(23)="HAC"
GFED3_EMFAC(23,1)=3.77E-03_hp
GFED3_EMFAC(23,2)=0.00E+00_hp
GFED3_EMFAC(23,3)=0.00E+00_hp
GFED3_EMFAC(23,4)=0.00E+00_hp
GFED3_EMFAC(23,5)=4.49E-04_hp
GFED3_EMFAC(23,6)=0.00E+00_hp

! CO2
GFED3_SPEC_NAME(24)="CO2"
GFED3_EMFAC(24,1)=1.31E+00_hp
GFED3_EMFAC(24,2)=1.63E+00_hp
GFED3_EMFAC(24,3)=1.57E+00_hp
GFED3_EMFAC(24,4)=1.70E+00_hp
GFED3_EMFAC(24,5)=1.65E+00_hp
GFED3_EMFAC(24,6)=1.64E+00_hp

! CH4
GFED3_SPEC_NAME(25)="CH4"
GFED3_EMFAC(25,1)=8.80E-03_hp
GFED3_EMFAC(25,2)=6.60E-03_hp
GFED3_EMFAC(25,3)=4.80E-03_hp
GFED3_EMFAC(25,4)=2.08E-02_hp
GFED3_EMFAC(25,5)=2.20E-03_hp
GFED3_EMFAC(25,6)=4.40E-03_hp

! EOH
GFED3_SPEC_NAME(26)="EOH"
GFED3_EMFAC(26,1)=7.90E-06_hp
GFED3_EMFAC(26,2)=8.65E-06_hp
GFED3_EMFAC(26,3)=2.86E-05_hp
GFED3_EMFAC(26,4)=1.80E-05_hp
GFED3_EMFAC(26,5)=3.22E-06_hp
GFED3_EMFAC(26,6)=1.80E-05_hp

! MTPA=APINE+BPINE+SABIN+CARE
GFED3_SPEC_NAME(27)="MTPA"
GFED3_EMFAC(27,1)=5.98E-04_hp
GFED3_EMFAC(27,2)=0.00E+00_hp
GFED3_EMFAC(27,3)=2.72E-03_hp
GFED3_EMFAC(27,4)=0.00E+00_hp
GFED3_EMFAC(27,5)=5.98E-05_hp
GFED3_EMFAC(27,6)=0.00E+00_hp

! HG0
GFED3_SPEC_NAME(28)="Hg0"
GFED3_EMFAC(28,1)=6.60E-08_hp
GFED3_EMFAC(28,2)=7.23E-08_hp
GFED3_EMFAC(28,3)=7.56E-08_hp
GFED3_EMFAC(28,4)=1.50E-07_hp
GFED3_EMFAC(28,5)=4.48E-08_hp
GFED3_EMFAC(28,6)=5.85E-08_hp

! SOAP
GFED3_SPEC_NAME(29)="SOAP"
GFED3_EMFAC(29,1)=GFED3_EMFAC(2,1)
GFED3_EMFAC(29,2)=GFED3_EMFAC(2,2)
GFED3_EMFAC(29,3)=GFED3_EMFAC(2,3)
GFED3_EMFAC(29,4)=GFED3_EMFAC(2,4)
GFED3_EMFAC(29,5)=GFED3_EMFAC(2,5)
GFED3_EMFAC(29,6)=GFED3_EMFAC(2,6)
!EOC
