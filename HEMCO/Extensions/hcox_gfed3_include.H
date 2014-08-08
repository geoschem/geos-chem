!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hcox_gfed3_include.H
!
! !DESCRIPTION: Include file with GFED3 emission factor data that was 
!  originally contained in file GFED3_emission_factors.txt.  We now hardwire
!  this so that we can avoid reading an ASCII file in the ESMF environment.
!
! !REMARKS:
!  This file was generated with the Perl script gfed3.pl, plus some
!  subsequent hand editing (to add the ProTeX header.
!
!  The GFED3_EMFAC array contains emission factors in kg/kgDM or kgC/kgDM
!  GFED3_EMFAC(N,1) = Agricultural Waste   Emission Factor for species N
!  GFED3_EMFAC(N,1) = Deforestation        Emission Factor for species N
!  GFED3_EMFAC(N,1) = Extratropical Forest Emission Factor for species N
!  GFED3_EMFAC(N,1) = Peat                 Emission Factor for species N
!  GFED3_EMFAC(N,1) = Savanna              Emission Factor for species N       
!  GFED3_EMFAC(N,1) = Woodland             Emission Factor for species N
! 
! !REVISION HISTORY: 
!  08 Aug 2014 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

   ! NO
   GFED3_SPEC_NAME(1  ) = "NO"
   GFED3_EMFAC    (1,1) = 2.29D-03
   GFED3_EMFAC    (1,2) = 2.26D-03
   GFED3_EMFAC    (1,3) = 3.41D-03
   GFED3_EMFAC    (1,4) = 2.26D-03
   GFED3_EMFAC    (1,5) = 2.12D-03
   GFED3_EMFAC    (1,6) = 2.19D-03

   ! CO
   GFED3_SPEC_NAME(2  ) = "CO"
   GFED3_EMFAC    (2,1) = 9.22D-02
   GFED3_EMFAC    (2,2) = 1.01D-01
   GFED3_EMFAC    (2,3) = 1.06D-01
   GFED3_EMFAC    (2,4) = 2.10D-01
   GFED3_EMFAC    (2,5) = 6.25D-02
   GFED3_EMFAC    (2,6) = 8.17D-02

   ! ALK4
   GFED3_SPEC_NAME(3  ) = "ALK4"
   GFED3_EMFAC    (3,1) = 1.91D-04
   GFED3_EMFAC    (3,2) = 1.92D-04
   GFED3_EMFAC    (3,3) = 2.99D-04
   GFED3_EMFAC    (3,4) = 1.92D-04
   GFED3_EMFAC    (3,5) = 1.11D-04
   GFED3_EMFAC    (3,6) = 1.53D-04

   ! ACET
   GFED3_SPEC_NAME(4  ) = "ACET"
   GFED3_EMFAC    (4,1) = 7.01D-04
   GFED3_EMFAC    (4,2) = 3.91D-04
   GFED3_EMFAC    (4,3) = 4.18D-04
   GFED3_EMFAC    (4,4) = 3.91D-04
   GFED3_EMFAC    (4,5) = 3.02D-04
   GFED3_EMFAC    (4,6) = 3.47D-04

   ! MEK
   GFED3_SPEC_NAME(5  ) = "MEK"
   GFED3_EMFAC    (5,1) = 6.89D-04
   GFED3_EMFAC    (5,2) = 3.91D-04
   GFED3_EMFAC    (5,3) = 3.87D-04
   GFED3_EMFAC    (5,4) = 3.91D-04
   GFED3_EMFAC    (5,5) = 1.93D-04
   GFED3_EMFAC    (5,6) = 2.93D-04

   ! ALD2
   GFED3_SPEC_NAME(6  ) = "ALD2"
   GFED3_EMFAC    (6,1) = 1.51D-03
   GFED3_EMFAC    (6,2) = 1.23D-03
   GFED3_EMFAC    (6,3) = 5.34D-04
   GFED3_EMFAC    (6,4) = 1.23D-03
   GFED3_EMFAC    (6,5) = 2.72D-04
   GFED3_EMFAC    (6,6) = 7.52D-04

   ! PRPE
   GFED3_SPEC_NAME(7  ) = "PRPE"
   GFED3_EMFAC    (7,1) = 4.90D-04
   GFED3_EMFAC    (7,2) = 6.54D-04
   GFED3_EMFAC    (7,3) = 4.89D-04
   GFED3_EMFAC    (7,4) = 6.54D-04
   GFED3_EMFAC    (7,5) = 3.15D-04
   GFED3_EMFAC    (7,6) = 4.85D-04

   ! C3H8
   GFED3_SPEC_NAME(8  ) = "C3H8"
   GFED3_EMFAC    (8,1) = 1.35D-04
   GFED3_EMFAC    (8,2) = 4.39D-04
   GFED3_EMFAC    (8,3) = 2.19D-04
   GFED3_EMFAC    (8,4) = 4.39D-04
   GFED3_EMFAC    (8,5) = 7.94D-05
   GFED3_EMFAC    (8,6) = 2.59D-04

   ! CH2O
   GFED3_SPEC_NAME(9  ) = "CH2O"
   GFED3_EMFAC    (9,1) = 2.12D-03
   GFED3_EMFAC    (9,2) = 2.22D-03
   GFED3_EMFAC    (9,3) = 2.24D-03
   GFED3_EMFAC    (9,4) = 2.22D-03
   GFED3_EMFAC    (9,5) = 8.66D-04
   GFED3_EMFAC    (9,6) = 1.54D-03

   ! C2H6
   GFED3_SPEC_NAME(10  ) = "C2H6"
   GFED3_EMFAC    (10,1) = 9.91D-04
   GFED3_EMFAC    (10,2) = 7.05D-04
   GFED3_EMFAC    (10,3) = 5.37D-04
   GFED3_EMFAC    (10,4) = 7.05D-04
   GFED3_EMFAC    (10,5) = 3.12D-04
   GFED3_EMFAC    (10,6) = 5.09D-04

   ! SO2
   GFED3_SPEC_NAME(11  ) = "SO2"
   GFED3_EMFAC    (11,1) = 4.00D-04
   GFED3_EMFAC    (11,2) = 7.59D-04
   GFED3_EMFAC    (11,3) = 9.99D-04
   GFED3_EMFAC    (11,4) = 7.59D-04
   GFED3_EMFAC    (11,5) = 3.37D-04
   GFED3_EMFAC    (11,6) = 5.48D-04

   ! NH3
   GFED3_SPEC_NAME(12  ) = "NH3"
   GFED3_EMFAC    (12,1) = 1.62D-03
   GFED3_EMFAC    (12,2) = 9.36D-04
   GFED3_EMFAC    (12,3) = 1.55D-03
   GFED3_EMFAC    (12,4) = 9.36D-04
   GFED3_EMFAC    (12,5) = 8.49D-04
   GFED3_EMFAC    (12,6) = 8.93D-04

   ! BC
   GFED3_SPEC_NAME(13  ) = "BC"
   GFED3_EMFAC    (13,1) = 4.21D-04
   GFED3_EMFAC    (13,2) = 5.65D-04
   GFED3_EMFAC    (13,3) = 5.62D-04
   GFED3_EMFAC    (13,4) = 5.65D-04
   GFED3_EMFAC    (13,5) = 4.61D-04
   GFED3_EMFAC    (13,6) = 5.13D-04

   ! OC
   GFED3_SPEC_NAME(14  ) = "OC"
   GFED3_EMFAC    (14,1) = 4.15D-03
   GFED3_EMFAC    (14,2) = 4.30D-03
   GFED3_EMFAC    (14,3) = 9.14D-03
   GFED3_EMFAC    (14,4) = 4.30D-03
   GFED3_EMFAC    (14,5) = 3.21D-03
   GFED3_EMFAC    (14,6) = 3.76D-03

   ! GLYX
   GFED3_SPEC_NAME(15  ) = "GLYX"
   GFED3_EMFAC    (15,1) = 9.95D-04
   GFED3_EMFAC    (15,2) = 9.95D-04
   GFED3_EMFAC    (15,3) = 9.95D-04
   GFED3_EMFAC    (15,4) = 9.95D-04
   GFED3_EMFAC    (15,5) = 9.95D-04
   GFED3_EMFAC    (15,6) = 9.95D-04

   ! MGLY
   GFED3_SPEC_NAME(16  ) = "MGLY"
   GFED3_EMFAC    (16,1) = 7.30D-04
   GFED3_EMFAC    (16,2) = 7.30D-04
   GFED3_EMFAC    (16,3) = 7.30D-04
   GFED3_EMFAC    (16,4) = 7.30D-04
   GFED3_EMFAC    (16,5) = 7.30D-04
   GFED3_EMFAC    (16,6) = 7.30D-04

   ! BENZ
   GFED3_SPEC_NAME(17  ) = "BENZ"
   GFED3_EMFAC    (17,1) = 2.83D-04
   GFED3_EMFAC    (17,2) = 3.50D-04
   GFED3_EMFAC    (17,3) = 4.81D-04
   GFED3_EMFAC    (17,4) = 3.50D-04
   GFED3_EMFAC    (17,5) = 2.68D-04
   GFED3_EMFAC    (17,6) = 3.09D-04

   ! TOLU
   GFED3_SPEC_NAME(18  ) = "TOLU"
   GFED3_EMFAC    (18,1) = 1.67D-04
   GFED3_EMFAC    (18,2) = 2.10D-04
   GFED3_EMFAC    (18,3) = 3.63D-04
   GFED3_EMFAC    (18,4) = 2.10D-04
   GFED3_EMFAC    (18,5) = 1.62D-04
   GFED3_EMFAC    (18,6) = 1.86D-04

   ! XYLE
   GFED3_SPEC_NAME(19  ) = "XYLE"
   GFED3_EMFAC    (19,1) = 1.03D-04
   GFED3_EMFAC    (19,2) = 7.88D-05
   GFED3_EMFAC    (19,3) = 1.85D-04
   GFED3_EMFAC    (19,4) = 7.88D-05
   GFED3_EMFAC    (19,5) = 3.89D-05
   GFED3_EMFAC    (19,6) = 5.89D-05

   ! C2H4
   GFED3_SPEC_NAME(20  ) = "C2H4"
   GFED3_EMFAC    (20,1) = 1.14D-03
   GFED3_EMFAC    (20,2) = 9.80D-04
   GFED3_EMFAC    (20,3) = 9.40D-04
   GFED3_EMFAC    (20,4) = 9.80D-04
   GFED3_EMFAC    (20,5) = 7.41D-04
   GFED3_EMFAC    (20,6) = 8.61D-04

   ! C2H2
   GFED3_SPEC_NAME(21  ) = "C2H2"
   GFED3_EMFAC    (21,1) = 2.00D-04
   GFED3_EMFAC    (21,2) = 1.82D-04
   GFED3_EMFAC    (21,3) = 2.38D-04
   GFED3_EMFAC    (21,4) = 1.82D-04
   GFED3_EMFAC    (21,5) = 2.92D-04
   GFED3_EMFAC    (21,6) = 2.37D-04

   ! GLYC
   GFED3_SPEC_NAME(22  ) = "GLYC"
   GFED3_EMFAC    (22,1) = 0.00E+00
   GFED3_EMFAC    (22,2) = 9.57D-04
   GFED3_EMFAC    (22,3) = 9.57D-04
   GFED3_EMFAC    (22,4) = 9.57D-04
   GFED3_EMFAC    (22,5) = 9.57D-04
   GFED3_EMFAC    (22,6) = 9.57D-04

   ! HAC
   GFED3_SPEC_NAME(23  ) = "HAC"
   GFED3_EMFAC    (23,1) = 0.00E+00
   GFED3_EMFAC    (23,2) = 0.00E+00
   GFED3_EMFAC    (23,3) = 0.00E+00
   GFED3_EMFAC    (23,4) = 0.00E+00
   GFED3_EMFAC    (23,5) = 0.00E+00
   GFED3_EMFAC    (23,6) = 0.00E+00

   ! CO2
   GFED3_SPEC_NAME(24  ) = "CO2"
   GFED3_EMFAC    (24,1) = 1.31E+00
   GFED3_EMFAC    (24,2) = 1.63E+00
   GFED3_EMFAC    (24,3) = 1.57E+00
   GFED3_EMFAC    (24,4) = 1.70E+00
   GFED3_EMFAC    (24,5) = 1.65E+00
   GFED3_EMFAC    (24,6) = 1.64E+00

   ! CH4
   GFED3_SPEC_NAME(25  ) = "CH4"
   GFED3_EMFAC    (25,1) = 8.80D-03
   GFED3_EMFAC    (25,2) = 6.60D-03
   GFED3_EMFAC    (25,3) = 4.80D-03
   GFED3_EMFAC    (25,4) = 2.08D-02
   GFED3_EMFAC    (25,5) = 2.20D-03
   GFED3_EMFAC    (25,6) = 4.40D-03
!EOC
