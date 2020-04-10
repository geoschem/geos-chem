!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hcox_gfed_include_gfed4.H
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
!  The Inst%GFED4_EMFAC array contains emission factors in kg/kgDM or kgC/kgDM
!  Inst%GFED4_EMFAC(N,1) = Savanna, grassland, and shrubland fires Emission Factor
!  Inst%GFED4_EMFAC(N,2) = Boreal forest fires Emission Factor
!  Inst%GFED4_EMFAC(N,3) = Temperate forest fires Forest Emission Factor
!  Inst%GFED4_EMFAC(N,4) = Deforestation (tropical) and degradation Emission Factor
!  Inst%GFED4_EMFAC(N,5) = Peat fires Emission Factor
!  Inst%GFED4_EMFAC(N,6) = Agricultural waste burning Emission Factor

! !REVISION HISTORY:
!  08 Aug 2014 - R. Yantosca - Initial version
!  10 Jul 2015 - R. Yantosca - Fixed minor issues in ProTeX header
!  08 Mar 2017 - M. Sulprizio- Added the following updates from E. Fischer:
!                              - Added EOH and MTPA for PAN updates
!  29 Mar 2018 - K. Travis   - Updated to latest GFED4 emission factors from
!                           https://www.geo.vu.nl/~gwerf/GFED/GFED4/ancill/
!                            - Added HAC, GLYX, Hg0
!EOP
!------------------------------------------------------------------------------
!BOC

! CO
GFED4_SPEC_NAME(1)="CO"
Inst%GFED4_EMFAC(1,1)=6.30E-02_hp
Inst%GFED4_EMFAC(1,2)=1.27E-01_hp
Inst%GFED4_EMFAC(1,3)=8.80E-02_hp
Inst%GFED4_EMFAC(1,4)=9.30E-02_hp
Inst%GFED4_EMFAC(1,5)=2.10E-01_hp
Inst%GFED4_EMFAC(1,6)=1.02E-01_hp

! ALK4 (Higher_Alkanes) - already in kgC
GFED4_SPEC_NAME(2)="ALK4"
Inst%GFED4_EMFAC(2,1)=5.50E-05_hp
Inst%GFED4_EMFAC(2,2)=3.49E-04_hp
Inst%GFED4_EMFAC(2,3)=2.25E-04_hp
Inst%GFED4_EMFAC(2,4)=7.20E-05_hp
Inst%GFED4_EMFAC(2,5)=7.20E-05_hp
Inst%GFED4_EMFAC(2,6)=3.40E-04_hp

! ACET (C3H6O) - convert to kgC
GFED4_SPEC_NAME(3)="ACET"
Inst%GFED4_EMFAC(3,1)=1.60E-04_hp * (3.0_hp*12.0_hp)/58.08_hp
Inst%GFED4_EMFAC(3,2)=7.50E-04_hp * (3.0_hp*12.0_hp)/58.08_hp
Inst%GFED4_EMFAC(3,3)=5.40E-04_hp * (3.0_hp*12.0_hp)/58.08_hp
Inst%GFED4_EMFAC(3,4)=6.30E-04_hp * (3.0_hp*12.0_hp)/58.08_hp
Inst%GFED4_EMFAC(3,5)=1.25E-03_hp * (3.0_hp*12.0_hp)/58.08_hp
Inst%GFED4_EMFAC(3,6)=4.50E-04_hp * (3.0_hp*12.0_hp)/58.08_hp

! MEK
GFED4_SPEC_NAME(4)="MEK"
Inst%GFED4_EMFAC(4,1)=1.81E-04_hp * (4.0_hp*12.0_hp)/72.11_hp
Inst%GFED4_EMFAC(4,2)=2.20E-04_hp * (4.0_hp*12.0_hp)/72.11_hp
Inst%GFED4_EMFAC(4,3)=1.30E-04_hp * (4.0_hp*12.0_hp)/72.11_hp
Inst%GFED4_EMFAC(4,4)=5.00E-04_hp * (4.0_hp*12.0_hp)/72.11_hp
Inst%GFED4_EMFAC(4,5)=5.00E-04_hp * (4.0_hp*12.0_hp)/72.11_hp
Inst%GFED4_EMFAC(4,6)=9.00E-04_hp * (4.0_hp*12.0_hp)/72.11_hp

! ALD2 (C2H4O)
GFED4_SPEC_NAME(5)="ALD2"
Inst%GFED4_EMFAC(5,1)=5.70E-04_hp * (2.0_hp*12.0_hp)/44.05_hp
Inst%GFED4_EMFAC(5,2)=7.70E-04_hp * (2.0_hp*12.0_hp)/44.05_hp
Inst%GFED4_EMFAC(5,3)=7.70E-04_hp * (2.0_hp*12.0_hp)/44.05_hp
Inst%GFED4_EMFAC(5,4)=1.55E-03_hp * (2.0_hp*12.0_hp)/44.05_hp
Inst%GFED4_EMFAC(5,5)=3.27E-03_hp * (2.0_hp*12.0_hp)/44.05_hp
Inst%GFED4_EMFAC(5,6)=1.24E-03_hp * (2.0_hp*12.0_hp)/44.05_hp

! PRPE (C3H6 + Higher_Alkenes)
GFED4_SPEC_NAME(6)="PRPE"
Inst%GFED4_EMFAC(6,1)=7.90E-04_hp * (3.0_hp*12.0_hp)/42.08_hp+1.33E-04_hp
Inst%GFED4_EMFAC(6,2)=1.13E-03_hp * (3.0_hp*12.0_hp)/42.08_hp+3.85E-04_hp
Inst%GFED4_EMFAC(6,3)=6.10E-04_hp * (3.0_hp*12.0_hp)/42.08_hp+3.69E-04_hp
Inst%GFED4_EMFAC(6,4)=6.40E-04_hp * (3.0_hp*12.0_hp)/42.08_hp+2.67E-04_hp
Inst%GFED4_EMFAC(6,5)=3.05E-03_hp * (3.0_hp*12.0_hp)/42.08_hp+2.27E-04_hp
Inst%GFED4_EMFAC(6,6)=6.80E-04_hp * (3.0_hp*12.0_hp)/42.08_hp+3.33E-04_hp

! C3H8
GFED4_SPEC_NAME(7)="C3H8"
Inst%GFED4_EMFAC(7,1)=1.00E-04_hp * (3.0_hp*12.0_hp)/44.1_hp
Inst%GFED4_EMFAC(7,2)=4.40E-04_hp * (3.0_hp*12.0_hp)/44.1_hp
Inst%GFED4_EMFAC(7,3)=2.20E-04_hp * (3.0_hp*12.0_hp)/44.1_hp
Inst%GFED4_EMFAC(7,4)=1.26E-04_hp * (3.0_hp*12.0_hp)/44.1_hp
Inst%GFED4_EMFAC(7,5)=1.26E-04_hp * (3.0_hp*12.0_hp)/44.1_hp
Inst%GFED4_EMFAC(7,6)=2.80E-04_hp * (3.0_hp*12.0_hp)/44.1_hp

! CH2O
GFED4_SPEC_NAME(8)="CH2O"
Inst%GFED4_EMFAC(8,1)=7.30E-04_hp
Inst%GFED4_EMFAC(8,2)=1.86E-03_hp
Inst%GFED4_EMFAC(8,3)=2.09E-03_hp
Inst%GFED4_EMFAC(8,4)=1.73E-03_hp
Inst%GFED4_EMFAC(8,5)=1.40E-03_hp
Inst%GFED4_EMFAC(8,6)=2.08E-03_hp

! C2H6
GFED4_SPEC_NAME(9)="C2H6"
Inst%GFED4_EMFAC(9,1)=6.60E-04_hp * (2.0_hp*12.0_hp)/30.07_hp
Inst%GFED4_EMFAC(9,2)=1.79E-03_hp * (2.0_hp*12.0_hp)/30.07_hp
Inst%GFED4_EMFAC(9,3)=6.30E-04_hp * (2.0_hp*12.0_hp)/30.07_hp
Inst%GFED4_EMFAC(9,4)=7.10E-04_hp * (2.0_hp*12.0_hp)/30.07_hp
Inst%GFED4_EMFAC(9,5)=7.10E-04_hp * (2.0_hp*12.0_hp)/30.07_hp
Inst%GFED4_EMFAC(9,6)=9.10E-04_hp * (2.0_hp*12.0_hp)/30.07_hp

! SO2
GFED4_SPEC_NAME(10)="SO2"
Inst%GFED4_EMFAC(10,1)=4.80E-04_hp
Inst%GFED4_EMFAC(10,2)=1.10E-03_hp
Inst%GFED4_EMFAC(10,3)=1.10E-03_hp
Inst%GFED4_EMFAC(10,4)=4.00E-04_hp
Inst%GFED4_EMFAC(10,5)=4.00E-04_hp
Inst%GFED4_EMFAC(10,6)=4.00E-04_hp

! NH3
GFED4_SPEC_NAME(11)="NH3"
Inst%GFED4_EMFAC(11,1)=5.20E-04_hp
Inst%GFED4_EMFAC(11,2)=2.72E-03_hp
Inst%GFED4_EMFAC(11,3)=8.40E-04_hp
Inst%GFED4_EMFAC(11,4)=1.33E-03_hp
Inst%GFED4_EMFAC(11,5)=1.33E-03_hp
Inst%GFED4_EMFAC(11,6)=2.17E-03_hp

! BC - already per carbon
GFED4_SPEC_NAME(12)="BC"
Inst%GFED4_EMFAC(12,1)=3.70E-04_hp
Inst%GFED4_EMFAC(12,2)=5.00E-04_hp
Inst%GFED4_EMFAC(12,3)=5.00E-04_hp
Inst%GFED4_EMFAC(12,4)=5.20E-04_hp
Inst%GFED4_EMFAC(12,5)=4.00E-05_hp
Inst%GFED4_EMFAC(12,6)=7.50E-04_hp

! OC - already per carbon
GFED4_SPEC_NAME(13)="OC"
Inst%GFED4_EMFAC(13,1)=2.62E-03_hp
Inst%GFED4_EMFAC(13,2)=9.60E-03_hp
Inst%GFED4_EMFAC(13,3)=9.60E-03_hp
Inst%GFED4_EMFAC(13,4)=4.71E-03_hp
Inst%GFED4_EMFAC(13,5)=6.02E-03_hp
Inst%GFED4_EMFAC(13,6)=2.30E-03_hp

! NO (NOx)
GFED4_SPEC_NAME(14)="NO"
Inst%GFED4_EMFAC(14,1)=3.90E-03_hp
Inst%GFED4_EMFAC(14,2)=9.00E-04_hp
Inst%GFED4_EMFAC(14,3)=1.92E-03_hp
Inst%GFED4_EMFAC(14,4)=2.55E-03_hp
Inst%GFED4_EMFAC(14,5)=1.00E-03_hp
Inst%GFED4_EMFAC(14,6)=3.11E-03_hp

! MGLY (CH3COCHO) - not traced per carbon
GFED4_SPEC_NAME(15)="MGLY"
Inst%GFED4_EMFAC(15,1)=7.30E-04_hp
Inst%GFED4_EMFAC(15,2)=7.30E-04_hp
Inst%GFED4_EMFAC(15,3)=7.30E-04_hp
Inst%GFED4_EMFAC(15,4)=7.30E-04_hp
Inst%GFED4_EMFAC(15,5)=7.30E-04_hp
Inst%GFED4_EMFAC(15,6)=7.30E-04_hp

! BENZ (C6H6)
GFED4_SPEC_NAME(16)="BENZ"
Inst%GFED4_EMFAC(16,1)=2.00E-04_hp * (6.0_hp*12.0_hp)/78.11_hp
Inst%GFED4_EMFAC(16,2)=1.11E-03_hp * (6.0_hp*12.0_hp)/78.11_hp
Inst%GFED4_EMFAC(16,3)=2.70E-04_hp * (6.0_hp*12.0_hp)/78.11_hp
Inst%GFED4_EMFAC(16,4)=3.90E-04_hp * (6.0_hp*12.0_hp)/78.11_hp
Inst%GFED4_EMFAC(16,5)=3.19E-03_hp * (6.0_hp*12.0_hp)/78.11_hp
Inst%GFED4_EMFAC(16,6)=1.50E-04_hp * (6.0_hp*12.0_hp)/78.11_hp

! TOLU (C7H8)
GFED4_SPEC_NAME(17)="TOLU"
Inst%GFED4_EMFAC(17,1)=8.00E-05_hp * (7.0_hp*12.0_hp)/92.14_hp
Inst%GFED4_EMFAC(17,2)=4.80E-04_hp * (7.0_hp*12.0_hp)/92.14_hp
Inst%GFED4_EMFAC(17,3)=1.90E-04_hp * (7.0_hp*12.0_hp)/92.14_hp
Inst%GFED4_EMFAC(17,4)=2.60E-04_hp * (7.0_hp*12.0_hp)/92.14_hp
Inst%GFED4_EMFAC(17,5)=1.55E-03_hp * (7.0_hp*12.0_hp)/92.14_hp
Inst%GFED4_EMFAC(17,6)=1.90E-04_hp * (7.0_hp*12.0_hp)/92.14_hp

! XYLE (C8H10)
GFED4_SPEC_NAME(18)="XYLE"
Inst%GFED4_EMFAC(18,1)=1.40E-05_hp * (8.0_hp*12.0_hp)/106.17_hp
Inst%GFED4_EMFAC(18,2)=1.80E-04_hp * (8.0_hp*12.0_hp)/106.17_hp
Inst%GFED4_EMFAC(18,3)=1.30E-04_hp * (8.0_hp*12.0_hp)/106.17_hp
Inst%GFED4_EMFAC(18,4)=1.10E-04_hp * (8.0_hp*12.0_hp)/106.17_hp
Inst%GFED4_EMFAC(18,5)=1.10E-04_hp * (8.0_hp*12.0_hp)/106.17_hp
Inst%GFED4_EMFAC(18,6)=1.14E-04_hp * (8.0_hp*12.0_hp)/106.17_hp

! C2H4
GFED4_SPEC_NAME(19)="C2H4"
Inst%GFED4_EMFAC(19,1)=8.20E-04_hp * (2.0_hp*12.0_hp)/28.05_hp
Inst%GFED4_EMFAC(19,2)=1.42E-03_hp * (2.0_hp*12.0_hp)/28.05_hp
Inst%GFED4_EMFAC(19,3)=1.17E-03_hp * (2.0_hp*12.0_hp)/28.05_hp
Inst%GFED4_EMFAC(19,4)=1.06E-03_hp * (2.0_hp*12.0_hp)/28.05_hp
Inst%GFED4_EMFAC(19,5)=2.57E-03_hp * (2.0_hp*12.0_hp)/28.05_hp
Inst%GFED4_EMFAC(19,6)=1.46E-03_hp * (2.0_hp*12.0_hp)/28.05_hp

! C2H2
GFED4_SPEC_NAME(20)="C2H2"
Inst%GFED4_EMFAC(20,1)=2.40E-04_hp * (2.0_hp*12.0_hp)/26.04_hp
Inst%GFED4_EMFAC(20,2)=1.80E-04_hp * (2.0_hp*12.0_hp)/26.04_hp
Inst%GFED4_EMFAC(20,3)=2.60E-04_hp * (2.0_hp*12.0_hp)/26.04_hp
Inst%GFED4_EMFAC(20,4)=4.40E-04_hp * (2.0_hp*12.0_hp)/26.04_hp
Inst%GFED4_EMFAC(20,5)=6.00E-05_hp * (2.0_hp*12.0_hp)/26.04_hp
Inst%GFED4_EMFAC(20,6)=2.70E-04_hp * (2.0_hp*12.0_hp)/26.04_hp

! GLYC (HOCH2CHO) - not carried per carbon
GFED4_SPEC_NAME(21)="GLYC"
Inst%GFED4_EMFAC(21,1)=2.50E-04_hp
Inst%GFED4_EMFAC(21,2)=8.60E-04_hp
Inst%GFED4_EMFAC(21,3)=8.60E-04_hp
Inst%GFED4_EMFAC(21,4)=7.40E-04_hp
Inst%GFED4_EMFAC(21,5)=7.40E-04_hp
Inst%GFED4_EMFAC(21,6)=7.10E-04_hp

! CO2
GFED4_SPEC_NAME(22)="CO2"
Inst%GFED4_EMFAC(22,1)=1.69E+00_hp
Inst%GFED4_EMFAC(22,2)=1.49E+00_hp
Inst%GFED4_EMFAC(22,3)=1.65E+00_hp
Inst%GFED4_EMFAC(22,4)=1.64E+00_hp
Inst%GFED4_EMFAC(22,5)=1.70E+00_hp
Inst%GFED4_EMFAC(22,6)=1.59E+00_hp

! CH4
GFED4_SPEC_NAME(23)="CH4"
Inst%GFED4_EMFAC(23,1)=1.94E-03_hp
Inst%GFED4_EMFAC(23,2)=5.96E-03_hp
Inst%GFED4_EMFAC(23,3)=3.36E-03_hp
Inst%GFED4_EMFAC(23,4)=5.07E-03_hp
Inst%GFED4_EMFAC(23,5)=2.08E-02_hp
Inst%GFED4_EMFAC(23,6)=5.82E-03_hp

! HCOOH - not carried per carbon
GFED4_SPEC_NAME(24)="HCOOH"
Inst%GFED4_EMFAC(24,1)=2.10E-04_hp
Inst%GFED4_EMFAC(24,2)=5.70E-04_hp
Inst%GFED4_EMFAC(24,3)=2.80E-04_hp
Inst%GFED4_EMFAC(24,4)=7.90E-04_hp
Inst%GFED4_EMFAC(24,5)=3.80E-04_hp
Inst%GFED4_EMFAC(24,6)=1.00E-03_hp

! DMS (C2H6S) - not carried per carbon
GFED4_SPEC_NAME(25)="DMS"
Inst%GFED4_EMFAC(25,1)=1.30E-06_hp
Inst%GFED4_EMFAC(25,2)=4.65E-06_hp
Inst%GFED4_EMFAC(25,3)=8.00E-06_hp
Inst%GFED4_EMFAC(25,4)=1.35E-06_hp
Inst%GFED4_EMFAC(25,5)=1.35E-06_hp
Inst%GFED4_EMFAC(25,6)=1.30E-06_hp

! ISOP (C5H8)
GFED4_SPEC_NAME(26)="ISOP"
Inst%GFED4_EMFAC(26,1)=3.90E-05_hp * (5.0_hp*12.0_hp)/68.12_hp
Inst%GFED4_EMFAC(26,2)=1.50E-04_hp * (5.0_hp*12.0_hp)/68.12_hp
Inst%GFED4_EMFAC(26,3)=9.90E-05_hp * (5.0_hp*12.0_hp)/68.12_hp
Inst%GFED4_EMFAC(26,4)=1.30E-04_hp * (5.0_hp*12.0_hp)/68.12_hp
Inst%GFED4_EMFAC(26,5)=1.38E-03_hp * (5.0_hp*12.0_hp)/68.12_hp
Inst%GFED4_EMFAC(26,6)=3.80E-04_hp * (5.0_hp*12.0_hp)/68.12_hp

! MTPA (C10H16) - not carried per carbon
GFED4_SPEC_NAME(27)="MTPA"
Inst%GFED4_EMFAC(27,1)=8.10E-05_hp
Inst%GFED4_EMFAC(27,2)=2.00E-03_hp
Inst%GFED4_EMFAC(27,3)=2.00E-03_hp
Inst%GFED4_EMFAC(27,4)=1.50E-04_hp
Inst%GFED4_EMFAC(27,5)=1.50E-04_hp
Inst%GFED4_EMFAC(27,6)=5.00E-06_hp

! MOH (CH3OH) - not carried per carbon
GFED4_SPEC_NAME(28)="MOH"
Inst%GFED4_EMFAC(28,1)=1.18E-03_hp
Inst%GFED4_EMFAC(28,2)=2.82E-03_hp
Inst%GFED4_EMFAC(28,3)=1.74E-03_hp
Inst%GFED4_EMFAC(28,4)=2.43E-03_hp
Inst%GFED4_EMFAC(28,5)=8.46E-03_hp
Inst%GFED4_EMFAC(28,6)=3.29E-03_hp

! EOH (C2H5OH)
GFED4_SPEC_NAME(29)="EOH"
Inst%GFED4_EMFAC(29,1)=2.40E-05_hp * (2.0_hp*12.0_hp)/46.07_hp
Inst%GFED4_EMFAC(29,2)=5.50E-05_hp * (2.0_hp*12.0_hp)/46.07_hp
Inst%GFED4_EMFAC(29,3)=1.00E-04_hp * (2.0_hp*12.0_hp)/46.07_hp
Inst%GFED4_EMFAC(29,4)=3.70E-05_hp * (2.0_hp*12.0_hp)/46.07_hp
Inst%GFED4_EMFAC(29,5)=3.70E-05_hp * (2.0_hp*12.0_hp)/46.07_hp
Inst%GFED4_EMFAC(29,6)=3.50E-05_hp * (2.0_hp*12.0_hp)/46.07_hp

! ACTA (CH3COOH) - not carried per carbon
GFED4_SPEC_NAME(30)="ACTA"
Inst%GFED4_EMFAC(30,1)=3.55E-03_hp
Inst%GFED4_EMFAC(30,2)=4.41E-03_hp
Inst%GFED4_EMFAC(30,3)=2.13E-03_hp
Inst%GFED4_EMFAC(30,4)=3.05E-03_hp
Inst%GFED4_EMFAC(30,5)=8.97E-03_hp
Inst%GFED4_EMFAC(30,6)=5.59E-03_hp

! SOAP
GFED4_SPEC_NAME(31)="SOAP"
Inst%GFED4_EMFAC(31,1)=Inst%GFED4_EMFAC(1,1) ! SAVA
Inst%GFED4_EMFAC(31,2)=Inst%GFED4_EMFAC(1,2) ! BORF
Inst%GFED4_EMFAC(31,3)=Inst%GFED4_EMFAC(1,3) ! TEMP
Inst%GFED4_EMFAC(31,4)=Inst%GFED4_EMFAC(1,4) ! DEFO
Inst%GFED4_EMFAC(31,5)=Inst%GFED4_EMFAC(1,5) ! PEAT
Inst%GFED4_EMFAC(31,6)=Inst%GFED4_EMFAC(1,6) ! AGRI

! HG0
GFED4_SPEC_NAME(32)="Hg0"
Inst%GFED4_EMFAC(32,1)=7.23E-08_hp ! SAV -- SAVA
Inst%GFED4_EMFAC(32,2)=1.50E-07_hp ! FOR -- BORF
Inst%GFED4_EMFAC(32,3)=1.50E-07_hp ! FOR -- TEMP
Inst%GFED4_EMFAC(32,4)=5.85E-08_hp ! DEF -- DEFO
Inst%GFED4_EMFAC(32,5)=7.56E-08_hp ! PET -- PEAT
Inst%GFED4_EMFAC(32,6)=4.48E-08_hp ! AGW -- AGRI

! HAC - not traced per carbon - from Akagi et al, 2011
GFED4_SPEC_NAME(33)="HAC"
Inst%GFED4_EMFAC(33,1)=4.50E-04_hp !SAV - Savannah
Inst%GFED4_EMFAC(33,2)=0.00E+00_hp !BORF - Boreal
Inst%GFED4_EMFAC(33,3)=0.00E+00_hp !TEMP - Temperate
Inst%GFED4_EMFAC(33,4)=1.13E-03_hp !DEFO - Tropical
Inst%GFED4_EMFAC(33,5)=1.92E-03_hp !PET -  Peatland
Inst%GFED4_EMFAC(33,6)=3.77E-03_hp !AGW - Crop residue

! GLYX - not traced per carbon
GFED4_SPEC_NAME(34)="GLYX"
Inst%GFED4_EMFAC(34,1)=9.95E-04_hp
Inst%GFED4_EMFAC(34,2)=9.95E-04_hp
Inst%GFED4_EMFAC(34,3)=9.95E-04_hp
Inst%GFED4_EMFAC(34,4)=9.95E-04_hp
Inst%GFED4_EMFAC(34,5)=9.95E-04_hp
Inst%GFED4_EMFAC(34,6)=9.95E-04_hp
!EOC
