This file documents all notable changes to the GEOS-Chem `fullchem` chemistry mechanism.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

# Mechanism history

## [14.1.0] - Feb 2023
### Added
- MO2 + NO3 reaction (IUPAC ROO_19)

### Fixed
- Bug fix: the "a" coefficient in rxn ETO = HO2 + 2.000CH2O was changed
  from 9.5d-13 (incorrect) to 9.5d+13.  See geoschem/geos-chem #1274.
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
- Use double precision numeric constants to each numeric value (e.g.
  `1.0d0` instead of `1.0e0`).  This will prevent having to call the
  `DBLE()` functions in the rate law functions, which wastes CPU cycles,
   and also is a "lossy" conversion. (BMY)

## [12.8.0] - Feb 2020
### Changed
- Update isoprene chemistry from Bates2019 (KHB)

## [12.7.0] -- Dec 2019
### Added
- MENO3, ETNO3, PRNO3 chemistry from Fisher2018 (JAF)
- "OTHRO2", which is equivalent to ETO2 but is not derived from C2H6
  oxidation. Necessary to prevent overestimates of ETNO3. All ETO2
  reactions are duplicated except ETO2+NO->ETNO3 channel (JAF)

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
  - Added several new NMVOCs. The extended mechanism includes ethanol,
    benzene, toluene and ethylbenzene (lumped), xylenes and trimethyl
    benzenes (lumped), and monoterpenes (lumped).
  - Treatment of monoterpene oxidation is adopted from the RACM2 chemical
    mechanism (Goliff et al., 2013), lumping terpenes with one double bond
    (alpha-pinene, beta-pinene, sabinene, delta-3-carene) into one proxy.

### Fixed
- ALK4 lumping issue in R4O2 + NO reaction (BHH)

## [v11-01g] - Sep 2016
### Added
- Initial version for FlexChem (MSL,MJE,MPS,EWL)


# Developers
- BHH:  Barron Henderson;        @barronh
- BMY:  Bob Yantosca;            @yantosca
- CCM:  Christopher Chan Miller; cmiller@fas.harvard.edu
- CDH:  Christopher Holmes;      @cdholmes
- DBM:  Dylan Millet;            @dylanbm
- EAM:  Eloise Marais;           @eamarais
- ECB:  Ellie Browne             eleanor.browne@colorado.edu
- EVF:  Emily Fischer;           evf@rams.colostate.edu
- EWL:  Lizzie Lundgren;         @lizziel
- FP :  Fabien Paulot;           fabien.paulot@noaa.gov
- HOTP: Havala Pye;              pye.havala@epa.gov
- JAF:  Jenny Fisher;            @jennyfisher
- JAS:  Johan Schmidt;           johanalbrechtschmidt@gmail.com
- JMAO: Jingqiu Mao;             @jingqiumao
- JMM:  Jonathan Moch            jmoch@g.harvard.edu
- JPP:  Justin Parrella;         justin.parrella@gmail.com
- KHB:  Kelvin Bates;            @kelvinhb
- KRT:  Katie Travis;            @ktravis213
- LZHU: Lei Zhu;                 leizhu@fas.harvard.edu
- MJE:  Mat Evans;               @msl3v
- MPS:  Melissa Sulprizio;       @msulprizio
- MSL:  Michael Long;            @msl3v
- PK:   Prasad Kasibhatla;       @pkasibhatla
- QJC:  Qianjie Chen;            chenqjie@uw.edu
- RHS:  Rebecca Schwantes;       rschwant@caltech.edu
- SAS:  Sarah Safieddine;        sarahsaf@mit.edu
- SDE:  Sebastian Eastham;       @sdeastham
- TMS:  Tomas Sherwen;           @tsherwen
- XC:   Xin Chen;                @xin-chen-github
- XW:   Xuan Wang;               @xuanw0316


# References (alphabetical order)
1.  Atkinson2003: https://doi.org/10.1021/cr0206420
2.  Bates2014:      https://doi.org/10.1021/jp4107958
3.  Bates2019:      https://doi.org/10.5194/acp-19-9613-2019
4.  Bates2021a:     https://doi.org/10.1029/2020JD033439
5.  Bates2021b:     https://doi.org/10.5194/acp-2021-605
6.  Browne2011:     https://doi.org/10.5194/acp-11-4209-2011
7.  Browne2014:     https://doi.org/10.5194/acp-14-1225-2014
8.  Chen2017:       https://doi.org/10.1002/2017GL073812
9.  Crounse2012:    https://doi.org/10.1021/jp211560u
10. Eastham2014:    https://doi.org/10.1016/j.atmosenv.2014.02.001
11. Fischer2014:    https://doi.org/10.5194/acp-14-2679-2014
12. Fisher2016:     https://doi.org/10.5194/acp-16-5969-2016
13. Fisher2018:     https://doi.org/10.1029/2018JD029046
14. Fry2014:        https://doi.org/10.1021/es502204x
15. Gill2002:       https://doi.org/10.1021/jp013532, 2002.
16. Goliff2013:     https://doi.org/10.1016/j.atmosenv.2012.11.038
17. Jacobs2014:     https://doi.org/10.5194/acp-14-8933-2014
18. Jenkin2015:     https://doi.org/10.5194/acp-15-11433-2015
19. Kasibhatla2018: https://doi.org/10.5194/acp-18-11185-2018
20. IUPAC ROO_19:   https://iupac-aeris.ipsl.fr/htdocs/datasheets/pdf/ROO_19_CH3O2_NO3.pdf
21. JPL 10-6:       https://jpldataeval.jpl.nasa.gov/previous_evaluations.html
22. JPL 15-10:      https://jpldataeval.jpl.nasa.gov, 2015.
23. Kwon2020:       https://doi.org/10.1525/elementa.2021.00109
24. Lee2014:        https://doi.org/10.1021/jp4107603
25. Marais2016:     https://doi.org/10.5194/acp-16-1603-2016
26. Miller2017:     https://doi.org/10.5194/acp-2016-1042
27. Millet2015:     https://doi.org/10.5194/acp-15-6283-2015
28. Moch2020:       https;//doi.org/10.1029/2020JD032706, 2020.
29. Muller2014:     https://doi.org/10.5194/acp-14-2497-2014
30. Parrella2012:   https://doi.org/10.5194/acp-12-6723-2012
31. Paulot2009:     https://doi.org/10.5194/acp-9-1479-2009 and https://doi.org/10.1126/science.1172910
32. Peeters2010:    https://doi.org/10.1039/C0CP00811G
33. Peeters2014:    https://doi.org/10.1021/jp5033146.
34. Pye2010:        https://doi.org/10.5194/acp-10-11261-2010
35. Roberts1992:    https://doi.org/10.1002/kin.550240307
36. Sherwen2016b:   https://doi.org/10.5194/acp-16-12239-2016
37. Sherwen2017:    https://doi.org/10.1039/C7FD00026J
38. StClair2016:    https://doi.org/10.1021/acs.jpca.5b065322016
39. Travis2016:     https://doi.org/10.5194/acp-16-13561-2016
40. Wolfe2012:      https://doi.org/ 10.1039/C2CP40388A, 2012
41. Xie2013:        https://doi.org/10.5194/acp-13-8439-2013
