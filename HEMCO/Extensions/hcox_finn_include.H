!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hcox_finn_include.H
!
! !DESCRIPTION: Include file with FINN emission factor data that was 
!  originally contained in files FINN_EFratios_CO2.csv and
!  FINN_VOC_speciation.csv.  We transform these data into hardwired F90
!  assignment statements so that we can avoid reading ASCII files in the 
!  ESMF environment.
!
! !REMARKS:
!  This file was generated with script HEMCO/Extensions/Preprocess/finn.pl.
! 
! !REVISION HISTORY: 
!  11 Aug 2014 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      !---------------------------------------------------
      ! Data from FINN_EFratios_CO2.csv
      !---------------------------------------------------

      ! Species from FINN_EFratios_CO2.csv
      FINN_SPEC_NAME( 1) = "CO2"
      FINN_SPEC_NAME( 2) = "CO"
      FINN_SPEC_NAME( 3) = "CH4"
      FINN_SPEC_NAME(58) = "H2"
      FINN_SPEC_NAME( 4) = "NOx"
      FINN_SPEC_NAME( 5) = "SO2"
      FINN_SPEC_NAME( 6) = "OC"
      FINN_SPEC_NAME( 7) = "BC"
      FINN_SPEC_NAME( 8) = "NH3"
      FINN_SPEC_NAME( 9) = "NO"
      FINN_SPEC_NAME(10) = "NO2"

      ! SavannaGrasslands emission factors from FINN_EFratios_CO2.csv
      EMFAC_IN( 1,1) = 1_dp
      EMFAC_IN( 2,1) = 17.0_dp
      EMFAC_IN( 3,1) = 316.7_dp
      EMFAC_IN(58,1) = 210.1_dp
      EMFAC_IN( 4,1) = 45.5_dp
      EMFAC_IN( 5,1) = 294.8_dp
      EMFAC_IN( 6,1) = 5112.7_dp
      EMFAC_IN( 7,1) = 14.7_dp
      EMFAC_IN( 8,1) = 103.5_dp
      EMFAC_IN( 9,1) = 1254.6_dp
      EMFAC_IN(10,1) = 1553.6_dp

      ! WoodySavannah emission factors from FINN_EFratios_CO2.csv
      EMFAC_IN( 1,2) = 1_dp
      EMFAC_IN( 2,2) = 15.3_dp
      EMFAC_IN( 3,2) = 226.5_dp
      EMFAC_IN(58,2) = 398.7_dp
      EMFAC_IN( 4,2) = 80.6_dp
      EMFAC_IN( 5,2) = 353.0_dp
      EMFAC_IN( 6,2) = 3645.4_dp
      EMFAC_IN( 7,2) = 10.5_dp
      EMFAC_IN( 8,2) = 29.8_dp
      EMFAC_IN( 9,2) = 523.0_dp
      EMFAC_IN(10,2) = 829.5_dp

      ! TropicalForest emission factors from FINN_EFratios_CO2.csv
      EMFAC_IN( 1,3) = 1_dp
      EMFAC_IN( 2,3) = 11.2_dp
      EMFAC_IN( 3,3) = 118.1_dp
      EMFAC_IN(58,3) = 97.6_dp
      EMFAC_IN( 4,3) = 22.4_dp
      EMFAC_IN( 5,3) = 439.4_dp
      EMFAC_IN( 6,3) = 5978.8_dp
      EMFAC_IN( 7,3) = 7.9_dp
      EMFAC_IN( 8,3) = 71.8_dp
      EMFAC_IN( 9,3) = 478.0_dp
      EMFAC_IN(10,3) = 1231.1_dp

      ! TemperateForest emission factors from FINN_EFratios_CO2.csv
      EMFAC_IN( 1,4) = 1_dp
      EMFAC_IN( 2,4) = 10.3_dp
      EMFAC_IN( 3,4) = 127.7_dp
      EMFAC_IN(58,4) = 109.3_dp
      EMFAC_IN( 4,4) = 30.3_dp
      EMFAC_IN( 5,4) = 555.2_dp
      EMFAC_IN( 6,4) = 2594.2_dp
      EMFAC_IN( 7,4) = 4.0_dp
      EMFAC_IN( 8,4) = 65.1_dp
      EMFAC_IN( 9,4) = 437.9_dp
      EMFAC_IN(10,4) = 1015.1_dp

      ! Boreal emission factors from FINN_EFratios_CO2.csv
      EMFAC_IN( 1,5) = 1_dp
      EMFAC_IN( 2,5) = 7.5_dp
      EMFAC_IN( 3,5) = 91.1_dp
      EMFAC_IN(58,5) = 78.5_dp
      EMFAC_IN( 4,5) = 29.7_dp
      EMFAC_IN( 5,5) = 1128.1_dp
      EMFAC_IN( 6,5) = 2167.4_dp
      EMFAC_IN( 7,5) = 4.3_dp
      EMFAC_IN( 8,5) = 169.2_dp
      EMFAC_IN( 9,5) = 211.8_dp
      EMFAC_IN(10,5) = 676.9_dp

      ! Crops emission factors from FINN_EFratios_CO2.csv
      EMFAC_IN( 1,6) = 1_dp
      EMFAC_IN( 2,6) = 9.9_dp
      EMFAC_IN( 3,6) = 99.3_dp
      EMFAC_IN(58,6) = 95.3_dp
      EMFAC_IN( 4,6) = 28.1_dp
      EMFAC_IN( 5,6) = 347.5_dp
      EMFAC_IN( 6,6) = 5767.7_dp
      EMFAC_IN( 7,6) = 15.7_dp
      EMFAC_IN( 8,6) = 48.0_dp
      EMFAC_IN( 9,6) = 282.6_dp
      EMFAC_IN(10,6) = 635.8_dp

      !---------------------------------------------------
      ! Data from FINN_VOC_speciation.csv
      !---------------------------------------------------

      ! Species from FINN_VOC_speciation.csv
      FINN_SPEC_NAME(11) = "ACET"
      FINN_SPEC_NAME(12) = "ACTA"
      FINN_SPEC_NAME(13) = "ALD2"
      FINN_SPEC_NAME(14) = "ALK4"
      FINN_SPEC_NAME(15) = "APINE"
      FINN_SPEC_NAME(16) = "AROM"
      FINN_SPEC_NAME(17) = "BENZ"
      FINN_SPEC_NAME(18) = "BPINE"
      FINN_SPEC_NAME(19) = "C2H2"
      FINN_SPEC_NAME(20) = "C2H4"
      FINN_SPEC_NAME(21) = "C2H6"
      FINN_SPEC_NAME(22) = "C3H8"
      FINN_SPEC_NAME(23) = "CARENE"
      FINN_SPEC_NAME(24) = "CH2Br2"
      FINN_SPEC_NAME(25) = "CH2O"
      FINN_SPEC_NAME(26) = "CH3Br"
      FINN_SPEC_NAME(27) = "CH3CN"
      FINN_SPEC_NAME(28) = "CH3I"
      FINN_SPEC_NAME(29) = "DMS"
      FINN_SPEC_NAME(30) = "EOH"
      FINN_SPEC_NAME(31) = "ETBENZ"
      FINN_SPEC_NAME(32) = "FUR"
      FINN_SPEC_NAME(33) = "GLYC"
      FINN_SPEC_NAME(34) = "HAC"
      FINN_SPEC_NAME(35) = "HCN"
      FINN_SPEC_NAME(36) = "HCOOH"
      FINN_SPEC_NAME(37) = "HNO2"
      FINN_SPEC_NAME(38) = "IGNOR"
      FINN_SPEC_NAME(39) = "ISOP"
      FINN_SPEC_NAME(40) = "LIMO"
      FINN_SPEC_NAME(41) = "MACR"
      FINN_SPEC_NAME(42) = "MEK"
      FINN_SPEC_NAME(43) = "MGLY"
      FINN_SPEC_NAME(44) = "MNO3"
      FINN_SPEC_NAME(45) = "MOH"
      FINN_SPEC_NAME(46) = "MTPO"
      FINN_SPEC_NAME(47) = "MVK "
      FINN_SPEC_NAME(48) = "PRPE"
      FINN_SPEC_NAME(49) = "R4N2"
      FINN_SPEC_NAME(50) = "RCHO"
      FINN_SPEC_NAME(51) = "RCOOH"
      FINN_SPEC_NAME(52) = "ROH"
      FINN_SPEC_NAME(53) = "SESQ"
      FINN_SPEC_NAME(54) = "STYR"
      FINN_SPEC_NAME(55) = "TMB"
      FINN_SPEC_NAME(56) = "TOLU"
      FINN_SPEC_NAME(57) = "XYLE"

      ! SavannaGrasslands emission factors from FINN_VOC_speciation.csv
      EMFAC_IN(11,1) = 1.70E-01_dp
      EMFAC_IN(12,1) = 3.54E+00_dp
      EMFAC_IN(13,1) = 7.80E-01_dp
      EMFAC_IN(14,1) = 1.20E-01_dp
      EMFAC_IN(15,1) = 0.00E+00_dp
      EMFAC_IN(16,1) = 3.30E-01_dp
      EMFAC_IN(17,1) = 1.50E-01_dp
      EMFAC_IN(18,1) = 0.00E+00_dp
      EMFAC_IN(19,1) = 5.50E-01_dp
      EMFAC_IN(20,1) = 1.75E+00_dp
      EMFAC_IN(21,1) = 1.32E+00_dp
      EMFAC_IN(22,1) = 1.40E-01_dp
      EMFAC_IN(23,1) = 0.00E+00_dp
      EMFAC_IN(24,1) = 0.00E+00_dp
      EMFAC_IN(25,1) = 1.46E+00_dp
      EMFAC_IN(26,1) = 0.00E+00_dp
      EMFAC_IN(27,1) = 1.60E-01_dp
      EMFAC_IN(28,1) = 0.00E+00_dp
      EMFAC_IN(29,1) = 0.00E+00_dp
      EMFAC_IN(30,1) = 1.00E-02_dp
      EMFAC_IN(31,1) = 0.00E+00_dp
      EMFAC_IN(32,1) = 6.70E-01_dp
      EMFAC_IN(33,1) = 8.10E-01_dp
      EMFAC_IN(34,1) = 7.60E-01_dp
      EMFAC_IN(35,1) = 9.10E-01_dp
      EMFAC_IN(36,1) = 2.70E-01_dp
      EMFAC_IN(37,1) = 2.50E-01_dp
      EMFAC_IN(38,1) = 1.90E-01_dp
      EMFAC_IN(39,1) = 3.00E-02_dp
      EMFAC_IN(40,1) = 0.00E+00_dp
      EMFAC_IN(41,1) = 0.00E+00_dp
      EMFAC_IN(42,1) = 1.21E+00_dp
      EMFAC_IN(43,1) = 6.10E-01_dp
      EMFAC_IN(44,1) = 0.00E+00_dp
      EMFAC_IN(45,1) = 5.50E-01_dp
      EMFAC_IN(46,1) = 1.00E-02_dp
      EMFAC_IN(47,1) = 0.00E+00_dp
      EMFAC_IN(48,1) = 1.79E+00_dp
      EMFAC_IN(49,1) = 0.00E+00_dp
      EMFAC_IN(50,1) = 1.90E-01_dp
      EMFAC_IN(51,1) = 0.00E+00_dp
      EMFAC_IN(52,1) = 5.00E-02_dp
      EMFAC_IN(53,1) = 0.00E+00_dp
      EMFAC_IN(54,1) = 1.00E-02_dp
      EMFAC_IN(55,1) = 0.00E+00_dp
      EMFAC_IN(56,1) = 5.00E-02_dp
      EMFAC_IN(57,1) = 1.00E-02_dp

      ! WoodySavannah emission factors from FINN_VOC_speciation.csv
      EMFAC_IN(11,2) = 7.10E-01_dp
      EMFAC_IN(12,2) = 1.54E+00_dp
      EMFAC_IN(13,2) = 9.50E-01_dp
      EMFAC_IN(14,2) = 5.70E-01_dp
      EMFAC_IN(15,2) = 0.00E+00_dp
      EMFAC_IN(16,2) = 4.00E-01_dp
      EMFAC_IN(17,2) = 3.00E-01_dp
      EMFAC_IN(18,2) = 0.00E+00_dp
      EMFAC_IN(19,2) = 6.50E-01_dp
      EMFAC_IN(20,2) = 2.25E+00_dp
      EMFAC_IN(21,2) = 1.01E+00_dp
      EMFAC_IN(22,2) = 3.60E-01_dp
      EMFAC_IN(23,2) = 0.00E+00_dp
      EMFAC_IN(24,2) = 0.00E+00_dp
      EMFAC_IN(25,2) = 2.33E+00_dp
      EMFAC_IN(26,2) = 0.00E+00_dp
      EMFAC_IN(27,2) = 4.10E-01_dp
      EMFAC_IN(28,2) = 0.00E+00_dp
      EMFAC_IN(29,2) = 0.00E+00_dp
      EMFAC_IN(30,2) = 2.00E-02_dp
      EMFAC_IN(31,2) = 1.00E-02_dp
      EMFAC_IN(32,2) = 9.50E-01_dp
      EMFAC_IN(33,2) = 3.20E-01_dp
      EMFAC_IN(34,2) = 0.00E+00_dp
      EMFAC_IN(35,2) = 1.18E+00_dp
      EMFAC_IN(36,2) = 1.10E-01_dp
      EMFAC_IN(37,2) = 7.30E-01_dp
      EMFAC_IN(38,2) = 2.20E-01_dp
      EMFAC_IN(39,2) = 3.00E-02_dp
      EMFAC_IN(40,2) = 0.00E+00_dp
      EMFAC_IN(41,2) = 0.00E+00_dp
      EMFAC_IN(42,2) = 1.19E+00_dp
      EMFAC_IN(43,2) = 8.50E-01_dp
      EMFAC_IN(44,2) = 0.00E+00_dp
      EMFAC_IN(45,2) = 5.30E-01_dp
      EMFAC_IN(46,2) = 1.00E-02_dp
      EMFAC_IN(47,2) = 0.00E+00_dp
      EMFAC_IN(48,2) = 1.58E+00_dp
      EMFAC_IN(49,2) = 0.00E+00_dp
      EMFAC_IN(50,2) = 2.60E-01_dp
      EMFAC_IN(51,2) = 0.00E+00_dp
      EMFAC_IN(52,2) = 8.00E-02_dp
      EMFAC_IN(53,2) = 0.00E+00_dp
      EMFAC_IN(54,2) = 2.00E-02_dp
      EMFAC_IN(55,2) = 0.00E+00_dp
      EMFAC_IN(56,2) = 1.70E-01_dp
      EMFAC_IN(57,2) = 3.00E-02_dp

      ! TropicalForest emission factors from FINN_VOC_speciation.csv
      EMFAC_IN(11,3) = 3.40E-01_dp
      EMFAC_IN(12,3) = 1.60E+00_dp
      EMFAC_IN(13,3) = 1.11E+00_dp
      EMFAC_IN(14,3) = 1.20E-01_dp
      EMFAC_IN(15,3) = 0.00E+00_dp
      EMFAC_IN(16,3) = 1.70E-01_dp
      EMFAC_IN(17,3) = 1.60E-01_dp
      EMFAC_IN(18,3) = 0.00E+00_dp
      EMFAC_IN(19,3) = 5.30E-01_dp
      EMFAC_IN(20,3) = 1.19E+00_dp
      EMFAC_IN(21,3) = 7.50E-01_dp
      EMFAC_IN(22,3) = 9.00E-02_dp
      EMFAC_IN(23,3) = 0.00E+00_dp
      EMFAC_IN(24,3) = 0.00E+00_dp
      EMFAC_IN(25,3) = 1.82E+00_dp
      EMFAC_IN(26,3) = 0.00E+00_dp
      EMFAC_IN(27,3) = 3.20E-01_dp
      EMFAC_IN(28,3) = 0.00E+00_dp
      EMFAC_IN(29,3) = 0.00E+00_dp
      EMFAC_IN(30,3) = 1.00E-02_dp
      EMFAC_IN(31,3) = 1.00E-02_dp
      EMFAC_IN(32,3) = 1.44E+00_dp
      EMFAC_IN(33,3) = 1.49E+00_dp
      EMFAC_IN(34,3) = 4.80E-01_dp
      EMFAC_IN(35,3) = 4.90E-01_dp
      EMFAC_IN(36,3) = 5.30E-01_dp
      EMFAC_IN(37,3) = 7.90E-01_dp
      EMFAC_IN(38,3) = 1.90E-01_dp
      EMFAC_IN(39,3) = 6.00E-02_dp
      EMFAC_IN(40,3) = 0.00E+00_dp
      EMFAC_IN(41,3) = 5.00E-02_dp
      EMFAC_IN(42,3) = 1.45E+00_dp
      EMFAC_IN(43,3) = 3.20E-01_dp
      EMFAC_IN(44,3) = 0.00E+00_dp
      EMFAC_IN(45,3) = 6.00E-01_dp
      EMFAC_IN(46,3) = 3.00E-02_dp
      EMFAC_IN(47,3) = 1.70E-01_dp
      EMFAC_IN(48,3) = 9.10E-01_dp
      EMFAC_IN(49,3) = 0.00E+00_dp
      EMFAC_IN(50,3) = 6.00E-01_dp
      EMFAC_IN(51,3) = 0.00E+00_dp
      EMFAC_IN(52,3) = 3.00E-02_dp
      EMFAC_IN(53,3) = 0.00E+00_dp
      EMFAC_IN(54,3) = 1.00E-02_dp
      EMFAC_IN(55,3) = 0.00E+00_dp
      EMFAC_IN(56,3) = 9.00E-02_dp
      EMFAC_IN(57,3) = 3.00E-02_dp

      ! TemperateForest emission factors from FINN_VOC_speciation.csv
      EMFAC_IN(11,4) = 3.10E-01_dp
      EMFAC_IN(12,4) = 1.69E+00_dp
      EMFAC_IN(13,4) = 4.30E-01_dp
      EMFAC_IN(14,4) = 2.90E-01_dp
      EMFAC_IN(15,4) = 2.60E-01_dp
      EMFAC_IN(16,4) = 8.50E-01_dp
      EMFAC_IN(17,4) = 3.20E-01_dp
      EMFAC_IN(18,4) = 2.30E-01_dp
      EMFAC_IN(19,4) = 1.80E-01_dp
      EMFAC_IN(20,4) = 1.23E+00_dp
      EMFAC_IN(21,4) = 1.34E+00_dp
      EMFAC_IN(22,4) = 2.30E-01_dp
      EMFAC_IN(23,4) = 1.00E-02_dp
      EMFAC_IN(24,4) = 0.00E+00_dp
      EMFAC_IN(25,4) = 1.56E+00_dp
      EMFAC_IN(26,4) = 0.00E+00_dp
      EMFAC_IN(27,4) = 3.30E-01_dp
      EMFAC_IN(28,4) = 0.00E+00_dp
      EMFAC_IN(29,4) = 0.00E+00_dp
      EMFAC_IN(30,4) = 3.00E-02_dp
      EMFAC_IN(31,4) = 1.00E-02_dp
      EMFAC_IN(32,4) = 6.80E-01_dp
      EMFAC_IN(33,4) = 2.90E-01_dp
      EMFAC_IN(34,4) = 0.00E+00_dp
      EMFAC_IN(35,4) = 1.29E+00_dp
      EMFAC_IN(36,4) = 2.80E-01_dp
      EMFAC_IN(37,4) = 2.30E-01_dp
      EMFAC_IN(38,4) = 4.20E-01_dp
      EMFAC_IN(39,4) = 5.00E-02_dp
      EMFAC_IN(40,4) = 1.30E-01_dp
      EMFAC_IN(41,4) = 2.00E-02_dp
      EMFAC_IN(42,4) = 8.40E-01_dp
      EMFAC_IN(43,4) = 2.50E-01_dp
      EMFAC_IN(44,4) = 0.00E+00_dp
      EMFAC_IN(45,4) = 5.20E-01_dp
      EMFAC_IN(46,4) = 1.80E-01_dp
      EMFAC_IN(47,4) = 7.00E-02_dp
      EMFAC_IN(48,4) = 1.40E+00_dp
      EMFAC_IN(49,4) = 1.00E-02_dp
      EMFAC_IN(50,4) = 2.90E-01_dp
      EMFAC_IN(51,4) = 2.00E-02_dp
      EMFAC_IN(52,4) = 4.50E-01_dp
      EMFAC_IN(53,4) = 1.00E-02_dp
      EMFAC_IN(54,4) = 1.00E-02_dp
      EMFAC_IN(55,4) = 2.00E-02_dp
      EMFAC_IN(56,4) = 1.20E-01_dp
      EMFAC_IN(57,4) = 8.00E-02_dp

      ! Boreal emission factors from FINN_VOC_speciation.csv
      EMFAC_IN(11,5) = 3.10E-01_dp
      EMFAC_IN(12,5) = 1.76E+00_dp
      EMFAC_IN(13,5) = 5.90E-01_dp
      EMFAC_IN(14,5) = 2.50E-01_dp
      EMFAC_IN(15,5) = 2.90E-01_dp
      EMFAC_IN(16,5) = 7.80E-01_dp
      EMFAC_IN(17,5) = 3.40E-01_dp
      EMFAC_IN(18,5) = 2.60E-01_dp
      EMFAC_IN(19,5) = 1.70E-01_dp
      EMFAC_IN(20,5) = 1.22E+00_dp
      EMFAC_IN(21,5) = 1.43E+00_dp
      EMFAC_IN(22,5) = 2.40E-01_dp
      EMFAC_IN(23,5) = 0.00E+00_dp
      EMFAC_IN(24,5) = 0.00E+00_dp
      EMFAC_IN(25,5) = 1.49E+00_dp
      EMFAC_IN(26,5) = 0.00E+00_dp
      EMFAC_IN(27,5) = 3.60E-01_dp
      EMFAC_IN(28,5) = 0.00E+00_dp
      EMFAC_IN(29,5) = 0.00E+00_dp
      EMFAC_IN(30,5) = 3.00E-02_dp
      EMFAC_IN(31,5) = 1.00E-02_dp
      EMFAC_IN(32,5) = 6.10E-01_dp
      EMFAC_IN(33,5) = 3.10E-01_dp
      EMFAC_IN(34,5) = 6.80E-01_dp
      EMFAC_IN(35,5) = 1.35E+00_dp
      EMFAC_IN(36,5) = 2.90E-01_dp
      EMFAC_IN(37,5) = 6.20E-01_dp
      EMFAC_IN(38,5) = 2.60E-01_dp
      EMFAC_IN(39,5) = 5.00E-02_dp
      EMFAC_IN(40,5) = 0.00E+00_dp
      EMFAC_IN(41,5) = 2.00E-02_dp
      EMFAC_IN(42,5) = 5.50E-01_dp
      EMFAC_IN(43,5) = 2.40E-01_dp
      EMFAC_IN(44,5) = 0.00E+00_dp
      EMFAC_IN(45,5) = 5.30E-01_dp
      EMFAC_IN(46,5) = 5.40E-01_dp
      EMFAC_IN(47,5) = 7.00E-02_dp
      EMFAC_IN(48,5) = 1.14E+00_dp
      EMFAC_IN(49,5) = 1.00E-02_dp
      EMFAC_IN(50,5) = 1.50E-01_dp
      EMFAC_IN(51,5) = 0.00E+00_dp
      EMFAC_IN(52,5) = 3.00E-02_dp
      EMFAC_IN(53,5) = 0.00E+00_dp
      EMFAC_IN(54,5) = 3.00E-02_dp
      EMFAC_IN(55,5) = 2.00E-02_dp
      EMFAC_IN(56,5) = 1.30E-01_dp
      EMFAC_IN(57,5) = 4.00E-02_dp

      ! Crops emission factors from FINN_VOC_speciation.csv
      EMFAC_IN(11,6) = 2.30E-01_dp
      EMFAC_IN(12,6) = 2.76E+00_dp
      EMFAC_IN(13,6) = 8.30E-01_dp
      EMFAC_IN(14,6) = 1.60E-01_dp
      EMFAC_IN(15,6) = 0.00E+00_dp
      EMFAC_IN(16,6) = 1.80E-01_dp
      EMFAC_IN(17,6) = 6.00E-02_dp
      EMFAC_IN(18,6) = 0.00E+00_dp
      EMFAC_IN(19,6) = 3.10E-01_dp
      EMFAC_IN(20,6) = 1.54E+00_dp
      EMFAC_IN(21,6) = 9.00E-01_dp
      EMFAC_IN(22,6) = 1.90E-01_dp
      EMFAC_IN(23,6) = 0.00E+00_dp
      EMFAC_IN(24,6) = 0.00E+00_dp
      EMFAC_IN(25,6) = 2.05E+00_dp
      EMFAC_IN(26,6) = 0.00E+00_dp
      EMFAC_IN(27,6) = 1.50E-01_dp
      EMFAC_IN(28,6) = 0.00E+00_dp
      EMFAC_IN(29,6) = 0.00E+00_dp
      EMFAC_IN(30,6) = 1.00E-02_dp
      EMFAC_IN(31,6) = 2.00E-02_dp
      EMFAC_IN(32,6) = 5.90E-01_dp
      EMFAC_IN(33,6) = 9.90E-01_dp
      EMFAC_IN(34,6) = 0.00E+00_dp
      EMFAC_IN(35,6) = 3.20E-01_dp
      EMFAC_IN(36,6) = 6.30E-01_dp
      EMFAC_IN(37,6) = 1.40E-01_dp
      EMFAC_IN(38,6) = 4.10E-01_dp
      EMFAC_IN(39,6) = 1.70E-01_dp
      EMFAC_IN(40,6) = 0.00E+00_dp
      EMFAC_IN(41,6) = 0.00E+00_dp
      EMFAC_IN(42,6) = 2.47E+00_dp
      EMFAC_IN(43,6) = 3.00E-01_dp
      EMFAC_IN(44,6) = 0.00E+00_dp
      EMFAC_IN(45,6) = 7.60E-01_dp
      EMFAC_IN(46,6) = 0.00E+00_dp
      EMFAC_IN(47,6) = 0.00E+00_dp
      EMFAC_IN(48,6) = 1.05E+00_dp
      EMFAC_IN(49,6) = 0.00E+00_dp
      EMFAC_IN(50,6) = 6.80E-01_dp
      EMFAC_IN(51,6) = 0.00E+00_dp
      EMFAC_IN(52,6) = 4.00E-02_dp
      EMFAC_IN(53,6) = 0.00E+00_dp
      EMFAC_IN(54,6) = 1.00E-02_dp
      EMFAC_IN(55,6) = 0.00E+00_dp
      EMFAC_IN(56,6) = 6.00E-02_dp
      EMFAC_IN(57,6) = 4.00E-02_dp

!EOC
