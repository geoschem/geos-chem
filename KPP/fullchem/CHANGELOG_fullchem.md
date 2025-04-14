This file documents all notable changes to the GEOS-Chem `fullchem` chemistry mechanism.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

# Mechanism history

## [14.5.0] - TBD
- Added four new species ALK4N1, ALK4N2, ALK4O2, and ALK4P to address issues in ALK4 and R4N2 chemistry following Brewer et al. (2023, JGR)
- Added ALK4N1 and ALK4N2 to Ox family in KPP
- Added new monoterpene mechanism, RCOOH oxidation, chemistry for new PNs (PHAN, AROMPN, MEKPN, APAN, LIMPAN, PINPAN) and aromatic AN (RNO3) per Travis et al., 2024
- Photolysis of PPN from Horner2024 (BH)
- PPN+OH and PAN+OH based on the structure activity relationship approach (BH)

## [14.4.0] - 2024-05-30
- Bug fix: Change `3.4.e-31` to `3.4.d-31` in `NO2 + O = NO + O2` rxn
- Update rxns with zero Arrhenius `B` parameters to use function `GCARR_ac` instead of `GCARR_abc`

## [14.3.0] - 2024-02-07
### Added
- Added PH2SO2 and PSO4AQ to track production of SO4 for use in TOMAS

### Changed
- Updated rate coefficients and products in 63 reactions per JPL / IUPAC recommendations (JPL 19-5; Bates2023)
- Consolidated product branches to remove 25 reactions (Bates2023 Table S4)

### Fixed
- Fixed C and N balance in 63 reactions (Bates2023 Table S2)
- Replaced the HMS + OH -> 2SO4 + CH2O - SO2 reaction with HMS + OH + SO2 -> 2SO4 + CH2O reaction and divided the rate constant be [SO2] to improve stability

## [14.2.1] - Oct 2023
### Changed
- The `fullchem` mechanism must now be built with KPP 3.0.0 or later

## [14.2.0] - Oct 2023
### Added
- Added lumped furan chemistry following Carter2020
- Restored sink reactions for HOI, IONO, IONO2
- Added S(IV)+HOBr and S(IV)+HOCl rxns (they had been inadvertently omitted)
- Added nitrate aerosol (NIT, NITs) to Ox family in `gckpp.kpp`

### Changed
- Set `k(SALAAL+SO2)` and `k(SALCAL+SO2)` to zero if O3 < 1e10 molec/cm3

### Fixed
- Fix bugs in HOBr uptake rate calculation in `fullchem_RateLawFuncs.F90`
- Now cap `State_Het%f_Alk_SSA` and `State_Het%f_Alk_SSC` at 1.0
- Unbalanced rxn IONO = ISALA is now balanced: IONO = ISALA + HNO2
- Unbalanced rxn IONO = ISALC is now balanced: IONO = ISALC + HNO2
- Unbalanced rxn IONO2 = ISALA is now balanced: IONO2 = ISALA + HNO3
- Unbalanced rxn IONO2 = ISALC is now balanced: IONO2 = ISALC + HNO3

### Changed
  - Restored sink reactions for HOI, IONO, IONO2
  - Use `GCARR_ac` for rxns where the Arrhenius `B` parameter is zero

## [14.1.0] - Feb 2023
### Added
- MO2 + NO3 reaction (IUPAC ROO_19)

### Fixed
- Bug fix: the "a" coefficient in rxn ETO = HO2 + 2.000CH2O was changed from 9.5d-13 (incorrect) to 9.5d+13.  See geoschem/geos-chem #1274.
- Bug fixes: HOBr + SO2 and HOCl + SO2 should produce SO4 (not SO4s)

## [13.3.0] - Sep 2021
### Added
- HMS chemistry (Moch2020)
- C2H2/C2H4 chemistry (Kwon2020)
- CH3O2 + OH reaction (Bates2021a)
- Sulfur reactions for future development (commented out)
- Aromatic SOA reactions (Bates2021b)

## [13.3.0] - Nov 2022
### Changed
- Use double precision numeric constants to each numeric value (e.g. `1.0d0` instead of `1.0e0`).  This will prevent having to call the `DBLE()` functions in the rate law functions, which wastes CPU cycles, and also is a "lossy" conversion. (BMY)

## [12.8.0] - Feb 2020
### Changed
- Update isoprene chemistry from Bates2019 (KHB)

## [12.7.0] -- Dec 2019
### Added
- MENO3, ETNO3, PRNO3 chemistry from Fisher2018 (JAF)
- "OTHRO2", which is equivalent to ETO2 but is not derived from C2H6  oxidation. Necessary to prevent overestimates of ETNO3. All ETO2 reactions are duplicated except ETO2+NO->ETNO3 channel (JAF)

### Changed
- Make MOH an active species (XC, DBM)

### Removed
- Old MNO3 species (same as MENO3 but not actually used (JAF)

## [12.6.0] -- Jul 2019
### Added
- Photolysis of NITs (off by default) (TMS, PK)
- Aerosol heterogeneous uptake for NOx (CDH)

## [v11-02d] - Sep 2017
### Added
- Halogen chemistry from Sherwen2016b/Sherwen2017 (TS,JAS,SDE,LZHU)
- HOBr + S(IV) from Chen2017 (QJC)

## [v11-02c] - Jul 2017
### Added
- Isoprene SOA updates from Marais2016 (EAM,MPS)
- Fixes for carbon-creating reactions (SAS,BHH,MJE)

### Changed
- Updated isoprene and monoterpene chemistry (KRT,JAF,CCM,EAM,KHB,RHS)
- Based on Travis2016, Fisher2016, ChanMiller2017, Marais2016
- Add Bates2014 epoxide scheme
- Update isoprene nitrate chemistry following Lee2014
- Add Muller2014 fast photolysis of carbonyl nitrates
- Add HNO2 chemistry from Lee2014
- Updated product yields and rx rate for RIO2+RIO2 (Xie2013)

## [v11-02a] - Mar 2017
### Changed
- Update rate constants based on JPL 15-10 (MJE,BHH)
- See wiki.geos-chem.org/Updates_in_JPL_Publication_15-10
- PAN chemistry updates (EVF)
- Added several new NMVOCs. The extended mechanism includes ethanol,  benzene, toluene and ethylbenzene (lumped), xylenes and trimethyl  benzenes (lumped), and monoterpenes (lumped).
- Treatment of monoterpene oxidation is adopted from the RACM2 chemical  mechanism (Goliff et al., 2013), lumping terpenes with one double bond  (alpha-pinene, beta-pinene, sabinene, delta-3-carene) into one proxy.

### Fixed
- ALK4 lumping issue in R4O2 + NO reaction (BHH)

## [v11-01g] - Sep 2016
### Added
- Initial version for FlexChem (MSL,MJE,MPS,EWL)

