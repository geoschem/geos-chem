! $Id: tracerid_mod.f,v 1.3 2010/03/15 19:33:20 ccarouge Exp $
      MODULE TRACERID_MOD
!
!******************************************************************************
!  Module TRACERID_MOD contains variables which point to SMVGEAR species,
!  CTM Tracers, Biomass species, and biofuel species located within various
!  GEOS-CHEM arrays. (bmy, ccc, 11/12/02, 5/17/10)
!
!  Module Variables:
!  ============================================================================
!  (1  ) NNNTRID   (INTEGER) : Max # of tracers
!  (2  ) MMMEMBER  (INTEGER) : Max # of members per chemical family
!  (3  ) NMEMBER   (INTEGER) : Number of members per each family tracer
!  (4  ) IDTRMB    (INTEGER) : Species # for each component of family tracer
!  (5  ) IDEMIS    (INTEGER) : Emission array for SMVGEAR
!  (6  ) CTRMB     (REAL*8 ) : Array for # of moles member/moles tracer
!  (7  ) IDO3      (INTEGER) : O3     index w/in CSPEC array ("comode_mod.f")
!  (8  ) IDNO2     (INTEGER) : NO2    index w/in CSPEC array ("comode_mod.f")
!  (9  ) IDNO3     (INTEGER) : NO3    index w/in CSPEC array ("comode_mod.f")
!  (10 ) IDN2O5    (INTEGER) : N2O5   index w/in CSPEC array ("comode_mod.f")
!  (11 ) IDHNO4    (INTEGER) : HNO4   index w/in CSPEC array ("comode_mod.f") 
!  (12 ) IDOX      (INTEGER) : OX     index w/in CSPEC array ("comode_mod.f")
!  (13 ) IDNOX     (INTEGER) : NOX    index w/in CSPEC array ("comode_mod.f") 
!  (14 ) IDHC1     (INTEGER) : HC1    index w/in CSPEC array ("comode_mod.f")
!  (15 ) IDNO      (INTEGER) : NO     index w/in CSPEC array ("comode_mod.f") 
!  (16 ) IDHNO2    (INTEGER) : HNO2   index w/in CSPEC array ("comode_mod.f") 
!  (17 ) IDCO      (INTEGER) : CO     index w/in CSPEC array ("comode_mod.f") 
!  (18 ) IDPRPE    (INTEGER) : PRPE   index w/in CSPEC array ("comode_mod.f")
!  (19 ) IDISOP    (INTEGER) : ISOP   index w/in CSPEC array ("comode_mod.f")
!  (20 ) IDALK4    (INTEGER) : ALK4   index w/in CSPEC array ("comode_mod.f")
!  (21 ) IDC3H8    (INTEGER) : C3H8   index w/in CSPEC array ("comode_mod.f") 
!  (22 ) IDPAN     (INTEGER) : PAN    index w/in CSPEC array ("comode_mod.f")
!  (23 ) IDGLPAN   (INTEGER) : GLPAN  index w/in CSPEC array ("comode_mod.f")
!  (24 ) IDGPAN    (INTEGER) : GPAN   index w/in CSPEC array ("comode_mod.f")
!  (25 ) IDPMN     (INTEGER) : PMN    index w/in CSPEC array ("comode_mod.f")
!  (26 ) IDPPN     (INTEGER) : PPN    index w/in CSPEC array ("comode_mod.f")
!  (27 ) IDHNO3    (INTEGER) : HNO3   index w/in CSPEC array ("comode_mod.f")
!  (28 ) IDOH      (INTEGER) : OH     index w/in CSPEC array ("comode_mod.f")
!  (29 ) IDHO2     (INTEGER) : HO2    index w/in CSPEC array ("comode_mod.f")
!  (30 ) IDH2O2    (INTEGER) : H2O2   index w/in CSPEC array ("comode_mod.f")
!  (31 ) IDACET    (INTEGER) : ACET   index w/in CSPEC array ("comode_mod.f")
!  (32 ) IDMEK     (INTEGER) : MEK    index w/in CSPEC array ("comode_mod.f")
!  (33 ) IDALD2    (INTEGER) : ALD2   index w/in CSPEC array ("comode_mod.f")
!  (34 ) IDRCHO    (INTEGER) : RCHO   index w/in CSPEC array ("comode_mod.f") 
!  (35 ) IDMVK     (INTEGER) : MVK    index w/in CSPEC array ("comode_mod.f") 
!  (36 ) IDMACR    (INTEGER) : MACR   index w/in CSPEC array ("comode_mod.f") 
!  (37 ) IDISN2    (INTEGER) : ISN2   index w/in CSPEC array ("comode_mod.f") 
!  (38 ) IDR4N2    (INTEGER) : R4N2   index w/in CSPEC array ("comode_mod.f")
!  (39 ) IDCH2O    (INTEGER) : CH2O   index w/in CSPEC array ("comode_mod.f")
!  (40 ) IDC2H6    (INTEGER) : C2H6   index w/in CSPEC array ("comode_mod.f")
!  (41 ) IDMP      (INTEGER) : MP     index w/in CSPEC array ("comode_mod.f")
!  (42 ) IDDMS     (INTEGER) : DMS    index w/in CSPEC array ("comode_mod.f")
!  (43 ) IDSO2     (INTEGER) : SO2    index w/in CSPEC array ("comode_mod.f")
!  (44 ) IDSO4     (INTEGER) : SO4    index w/in CSPEC array ("comode_mod.f")
!  (45 ) IDMSA     (INTEGER) : MSA    index w/in CSPEC array ("comode_mod.f")
!  (46 ) IDDRYO3   (INTEGER) : DRYO3  index w/in CSPEC array ("comode_mod.f")
!  (47 ) IDDRYPAN  (INTEGER) : DRYPAN index w/in CSPEC array ("comode_mod.f")
!  (48 ) IDDRYNO2  (INTEGER) : DRYNO2 index w/in CSPEC array ("comode_mod.f")
!  (49 ) IDTNOX    (INTEGER) : NOx   index w/in STT array ("tracer_mod.f")
!  (50 ) IDTOX     (INTEGER) : Ox    index w/in STT array ("CMN_mod.F")     
!  (51 ) IDTPAN    (INTEGER) : PAN   index w/in STT array ("tracer_mod.f")
!  (52 ) IDTCO     (INTEGER) : CO    index w/in STT array ("tracer_mod.f")     
!  (53 ) IDTALK4   (INTEGER) : ALK4  index w/in STT array ("tracer_mod.f")    
!  (54 ) IDTISOP   (INTEGER) : ISOP  index w/in STT array ("tracer_mod.f")    
!  (55 ) IDTHNO3   (INTEGER) : HNO3  index w/in STT array ("tracer_mod.f")     
!  (56 ) IDTH2O2   (INTEGER) : H2O2  index w/in STT array ("tracer_mod.f")    
!  (57 ) IDTACET   (INTEGER) : ACET  index w/in STT array ("tracer_mod.f")    
!  (58 ) IDTMEK    (INTEGER) : MEK   index w/in STT array ("tracer_mod.f")    
!  (59 ) IDTALD2   (INTEGER) : ALD2  index w/in STT array ("tracer_mod.f")    
!  (60 ) IDTRCHO   (INTEGER) : RCHO  index w/in STT array ("tracer_mod.f")     
!  (61 ) IDTMVK    (INTEGER) : MVK   index w/in STT array ("tracer_mod.f")    
!  (62 ) IDTMACR   (INTEGER) : MACR  index w/in STT array ("tracer_mod.f")    
!  (63 ) IDTPMN    (INTEGER) : PMN   index w/in STT array ("tracer_mod.f")    
!  (64 ) IDTPPN    (INTEGER) : PPN   index w/in STT array ("tracer_mod.f")    
!  (65 ) IDTISN2   (INTEGER) : ISN2  index w/in STT array ("tracer_mod.f")     
!  (66 ) IDTR4N2   (INTEGER) : R4N2  index w/in STT array ("tracer_mod.f")    
!  (67 ) IDTPRPE   (INTEGER) : PRPE  index w/in STT array ("tracer_mod.f")    
!  (68 ) IDTC3H8   (INTEGER) : C3H8  index w/in STT array ("tracer_mod.f")    
!  (69 ) IDTCH2O   (INTEGER) : CH2O  index w/in STT array ("tracer_mod.f")    
!  (70 ) IDTMP     (INTEGER) : MP    index w/in STT array ("tracer_mod.f")    
!  (71 ) IDTN2O5   (INTEGER) : N2O5  index w/in STT array ("tracer_mod.f")    
!  (72 ) IDTHNO4   (INTEGER) : HNO4  index w/in STT array ("tracer_mod.f")    
!  (73 ) IDTC2H6   (INTEGER) : C2H6  index w/in STT array ("tracer_mod.f")    
!  (74 ) IDTDMS    (INTEGER) : DMS   index w/in STT array ("tracer_mod.f")     
!  (75 ) IDTSO2    (INTEGER) : SO2   index w/in STT array ("tracer_mod.f")
!  (76 ) IDTSO4    (INTEGER) : SO4   index w/in STT array ("tracer_mod.f")    
!  (77 ) IDTSO4aq  (INTEGER) : SO4aq index w/in STT array ("tracer_mod.f")    
!  (78 ) IDTSO4s   (INTEGER) : SO4s  index w/in STT array ("tracer_mod.f")    
!  (79 ) IDTMSA    (INTEGER) : MSA   index w/in STT array ("tracer_mod.f")    
!  (80 ) IDTNH3    (INTEGER) : NH3   index w/in STT array ("tracer_mod.f")    
!  (81 ) IDTNH4    (INTEGER) : NH4   index w/in STT array ("tracer_mod.f")    
!  (82 ) IDTNIT    (INTEGER) : NIT   index w/in STT array ("tracer_mod.f")    
!  (83 ) IDTNITs   (INTEGER) : NITs  index w/in STT array ("tracer_mod.f")    
!  (84 ) IDTRN     (INTEGER) : Rn    index w/in STT array ("tracer_mod.f")    
!  (85 ) IDTPB     (INTEGER) : Pb    index w/in STT array ("tracer_mod.f")    
!  (86 ) IDTBE7    (INTEGER) : Be7   index w/in STT array ("tracer_mod.f")    
!  (87 ) IDTBCPI   (INTEGER) : BCPI  index w/in STT array ("tracer_mod.f")
!  (88 ) IDTBCPO   (INTEGER) : BCPO  index w/in STT array ("tracer_mod.f")
!  (89 ) IDTOCPI   (INTEGER) : OCPI  index w/in STT array ("tracer_mod.f")
!  (90 ) IDTOCPO   (INTEGER) : OCPO  index w/in STT array ("tracer_mod.f")
!  (91 ) IDTALPH   (INTEGER) : ALPH  index w/in STT array ("tracer_mod.f")
!  (92 ) IDTLIMO   (INTEGER) : LIMO  index w/in STT array ("tracer_mod.f")
!  (93 ) IDTALCO   (INTEGER) : ALCO  index w/in STT array ("tracer_mod.f")
!  (94 ) IDTSOG1   (INTEGER) : SOG1  index w/in STT array ("tracer_mod.f")
!  (95 ) IDTSOG2   (INTEGER) : SOG2  index w/in STT array ("tracer_mod.f")
!  (96 ) IDTSOG3   (INTEGER) : SOG3  index w/in STT array ("tracer_mod.f")
!  (97 ) IDTSOA1   (INTEGER) : SOA1  index w/in STT array ("tracer_mod.f")
!  (98 ) IDTSOA2   (INTEGER) : SOA2  index w/in STT array ("tracer_mod.f")
!  (99 ) IDTSOA3   (INTEGER) : SOA3  index w/in STT array ("tracer_mod.f")
!  (100) IDTDST1   (INTEGER) : DST1  index w/in STT array ("tracer_mod.f")
!  (101) IDTDST2   (INTEGER) : DST2  index w/in STT array ("tracer_mod.f")
!  (102) IDTDST3   (INTEGER) : DST3  index w/in STT array ("tracer_mod.f")
!  (103) IDTDST4   (INTEGER) : DST4  index w/in STT array ("tracer_mod.f")
!  (104) IDTSALA   (INTEGER) : SALA  index w/in STT array ("tracer_mod.f")
!  (105) IDTSALC   (INTEGER) : SALC  index w/in STT array ("tracer_mod.f")
!  (109) IDENO     (INTEGER) : NOx   index w/in EMISRRN array ("CMN_O3")  
!  (110) IDEO3     (INTEGER) : Ox    index w/in EMISRR  array ("CMN_O3")  
!  (111) IDECO     (INTEGER) : CO    index w/in EMISRR  array ("CMN_O3")     
!  (112) IDEPRPE   (INTEGER) : PRPE  index w/in EMISRR  array ("CMN_O3")     
!  (113) IDEC3H8   (INTEGER) : C3H8  index w/in EMISRR  array ("CMN_O3")     
!  (114) IDEALK4   (INTEGER) : ALK4  index w/in EMISRR  array ("CMN_O3")     
!  (115) IDEC2H6   (INTEGER) : C2H6  index w/in EMISRR  array ("CMN_O3")    
!  (116) IDEISOP   (INTEGER) : ISOP  index w/in EMISRR  array ("CMN_O3")    
!  (117) IDEACET   (INTEGER) : ACET  index w/in EMISRR  array ("CMN_O3")     
!  (118) IDEMEK    (INTEGER) : MEK   index w/in EMISRR  array ("CMN_O3")     
!  (119) IDEALD2   (INTEGER) : ALD2  index w/in EMISRR  array ("CMN_O3")    
!  (120) IDECH2O   (INTEGER) : CH2O  index w/in EMISRR  array ("CMN_O3")     
!  (121) NEMBIOG   (INTEGER) : # of biogenic emission species for SMVGEAR
!  (122) NEMANTHRO (INTEGER) : # of anthro   emission species for SMVGEAR
!  (132) IDBFPRPE  (INTEGER) : PRPE index w/in BURNEMIS array (biofuel_mod.f)
!  (133) IDBALK4   (INTEGER) : ALD4 index w/in BURNEMIS array (biomass_mod.f)
!  (134) IDBFNOX   (INTEGER) : NOx  index w/in BIOFUEL array (biofuel_mod.f)
!  (135) IDBFCO    (INTEGER) : CO   index w/in BIOFUEL array (biofuel_mod.f)
!  (136) IDBFALK4  (INTEGER) : ALK4 index w/in BIOFUEL array (biofuel_mod.f)
!  (137) IDBFACET  (INTEGER) : ACET index w/in BIOFUEL array (biofuel_mod.f)
!  (138) IDBFMEK   (INTEGER) : MEK  index w/in BIOFUEL array (biofuel_mod.f)
!  (139) IDBFALD2  (INTEGER) : ALD2 index w/in BIOFUEL array (biofuel_mod.f)
!  (140) IDBFPRPE  (INTEGER) : PRPE index w/in BIOFUEL array (biofuel_mod.f)
!  (141) IDBFC3H8  (INTEGER) : NOx  index w/in BIOFUEL array (biofuel_mod.f)
!  (142) IDBFCH2O  (INTEGER) : NOx  index w/in BIOFUEL array (biofuel_mod.f)
!  (143) IDBFC2H6  (INTEGER) : NOx  index w/in BIOFUEL array (biofuel_mod.f)   
!  (144) IDNAP     (INTEGER) : NAP   index w/in CSPEC array ("comode_mod.f")
!  (145) IDNRO2    (INTEGER) : NRO2  index w/in CSPEC array ("comode_mod.f")
!  (146) IDLNRO2H  (INTEGER) : LRO2H index w/in CSPEC array ("comode_mod.f")
!  (147) IDLNRO2N  (INTEGER) : LRO2N index w/in CSPEC array ("comode_mod.f")
!  (148) IDTNAP    (INTEGER) : NAP   index w/in STT array ("tracer_mod.f")
!  (149) IDTMTPA   (INTEGER) : MTPA  index w/in STT array ("tracer_mod.f")
!  (150) IDTMTPO   (INTEGER) : MTPO  index w/in STT array ("tracer_mod.f")
!  (151) IDTASOAN  (INTEGER) : ASOAN index w/in STT array ("tracer_mod.f")
!  (152) IDTASOA1  (INTEGER) : ASOA1 index w/in STT array ("tracer_mod.f")
!  (153) IDTASOA2  (INTEGER) : ASOA2 index w/in STT array ("tracer_mod.f")
!  (154) IDTASOA3  (INTEGER) : ASOA3 index w/in STT array ("tracer_mod.f")
!  (155) IDTASOG1  (INTEGER) : ASOG1 index w/in STT array ("tracer_mod.f")
!  (156) IDTASOG2  (INTEGER) : ASOG2 index w/in STT array ("tracer_mod.f")
!  (157) IDTASOG3  (INTEGER) : ASOG3 index w/in STT array ("tracer_mod.f")
!  (158) IDTTSOA1  (INTEGER) : TSOA1 index w/in STT array ("tracer_mod.f")
!  (159) IDTTSOA2  (INTEGER) : TSOA2 index w/in STT array ("tracer_mod.f")
!  (160) IDTTSOA3  (INTEGER) : TSOA3 index w/in STT array ("tracer_mod.f")
!  (161) IDTTSOG1  (INTEGER) : TSOG1 index w/in STT array ("tracer_mod.f")
!  (162) IDTTSOG2  (INTEGER) : TSOG2 index w/in STT array ("tracer_mod.f")
!  (163) IDTTSOG3  (INTEGER) : TSOG3 index w/in STT array ("tracer_mod.f")
!  (164) IDTTSOA0  (INTEGER) : TSOA0 index w/in STT array ("tracer_mod.f")
!  (165) IDTTSOG0  (INTEGER) : TSOG0 index w/in STT array ("tracer_mod.f")
!  (166) IDTISOA1  (INTEGER) : ISOA1 index w/in STT array ("tracer_mod.f")
!  (167) IDTISOA2  (INTEGER) : ISOA2 index w/in STT array ("tracer_mod.f")
!  (168) IDTISOA3  (INTEGER) : ISOA3 index w/in STT array ("tracer_mod.f")
!  (169) IDTISOG1  (INTEGER) : ISOG1 index w/in STT array ("tracer_mod.f")
!  (170) IDTISOG2  (INTEGER) : ISOG2 index w/in STT array ("tracer_mod.f")
!  (171) IDTISOG3  (INTEGER) : ISOG3 index w/in STT array ("tracer_mod.f")
!  (172) IDTPOA1   (INTEGER) : POA1  index w/in STT array ("tracer_mod.f")
!  (173) IDTPOA2   (INTEGER) : POA2  index w/in STT array ("tracer_mod.f")
!  (174) IDTPOG1   (INTEGER) : POG1  index w/in STT array ("tracer_mod.f")
!  (175) IDTPOG2   (INTEGER) : POG2  index w/in STT array ("tracer_mod.f")
!  (176) IDTOPOA1  (INTEGER) : OPOA1 index w/in STT array ("tracer_mod.f")
!  (177) IDTOPOA2  (INTEGER) : OPOA2 index w/in STT array ("tracer_mod.f")
!  (178) IDTOPOG1  (INTEGER) : OPOG1 index w/in STT array ("tracer_mod.f")
!  (179) IDTOPOG2  (INTEGER) : OPOG2 index w/in STT array ("tracer_mod.f")
!  (180) IDBFNAP   (INTEGER) : FNAP  index w/in BIOFUEL array (biofuel_mod.f)
!  (181) IDBNAP    (INTEGER) : NAP   index w/in BIOFUEL array (biofuel_mod.f)
!
!  Module Routines:
!  ============================================================================
!  (1 ) TRACERID      : Defines tracer, biomass, biofuel, & anthro ID numbers
!  (2 ) SETTRACE      : Defines ID numbers for species in SMVGEAR mechanism
!  (3 ) INIT_TRACERID : Zeroes all module variables
!
!  GEOS-CHEM modules referenced by tracerid_mod.f
!  ============================================================================
!  (1 ) charpak_mod.f : Module containing string handling routines
!  (2 ) error_mod.f   : Module containing I/O error and NaN check routines
!
!  References:
!  ============================================================================
!  (1 ) Pye, H.O.T., A.W.H. Chan, M.O. Barkley, and J.H. Seinfeld, "Global
!        modeling of organic aerosol: The importance of reactive nitrogen (NOx
!        and NO3)", Atmos. Chem. Phys., Vol 10, pp 11261-11276, 2010.
!
!  NOTES:
!  (1 ) Added additional SMVGEAR species flags for DMS, SO2, SO4, MSA, so that
!        these species can be handled w/in SMVGEAR (rjp, bmy, 3/23/03)
!  (2 ) Added modifications for SMVGEAR II (bdf, bmy, 4/1/03)
!  (3 ) Added extra flags for carbon & dust tracers (rjp, tdf, bmy, 4/1/04)
!  (4 ) Added extra flags for seasalt tracers (rjp, bec, bmy, 4/20/04)
!  (5 ) Increase NNNTRID for carb+dust+seasalt tracers (bmy, 4/26/04)
!  (6 ) Increase NNNTRID & add extra flags for SOA tracers. (rjp, bmy, 7/13/04)
!  (7 ) Bug fix: reverse IDECH2O and IDEISOP (bmy, 11/15/04)
!  (8 ) Added IDTHG0, IDTHG2, IDTHGP + tagged Hg's (eck, bmy, 12/7/04)
!  (9 ) Added IDTAS, IDTAHS, IDTLET, IDTNH4aq, IDTSO4aq (cas, bmy, 12/20/04)
!  (10) Added IDTSO4s, IDTNITs
!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (12) Added functions IS_Hg0, IS_Hg2, IS_HgP.  Also now use index arrays
!        ID_Hg0, ID_Hg2, ID_HgP for tagged Hg tracers.  (cdh, bmy, 1/5/06)
!  (13) Remove IDBxxxx biomass flags; these aren't needed. (bmy, 4/5/06)
!  (14) Add IDTSOG4 and IDTSOA4 (dkh, bmy, 5/18/06)
!  (15) Minor fixes for CH3I simulation (bmy, 7/25/06)
!  (16) Add IDTH2 and IDTHD for H2/HD simulation (hup, lyj, phs, 9/18/07)
!  (17) Set IDECO=1 for Tagged CO simulation (jaf, mak, bmy, 2/14/08)
!  (18) Add IDEHNO3 to deal with ship NOx emissions (phs, 3/4/08)
!  (19) Added tracers and emissions for dicarbonyl simulation (tmf, 1/7/09)
!  (20) Increase NNNTRID for TOMAS species 30bin*(8mass+1number) from 
!       73 to 320 (win, 6/23/09)
!  (21) Added IDH2SO4 for SMVGEAR species (win, 6/23/09)
!  (22) Added TOMAS tracers IDTH2SO4, IDTNK1-IDTNK30, IDTSF1-IDTSF30, 
!        IDTSS1-IDTSS30, IDTECIL1-IDTECIL30, IDTECOB1-IDTECOB30, 
!        IDTOCIL1-IDTOCIL30, IDTOCOB1-IDTOCOB30, and IDTAW1-IDTAW30 (win, 6/23/09)
!  (23) Updates for mercury simulation (ccc, 5/17/10)
!  17 Feb 2011 - R. Yantosca - Add modifications for APM microphysics (G. Luo)
!  07 Feb 2012 - E. Corbitt  - Add mercury tagged tracers
!  13 Aug 2012 - R. Yantosca - Add routine FILL_State_Chm_NAME_IDs
!  11 Jan 2013 - M. Payer    - Add methyl peroxy nitrate (MPN) as tracer
!                              (E. Browne)
!  07 Feb 2013 - M. Payer    - Add NO, NO2, NO2, HNO2, and O3 as tracers to
!                              replace NOx and Ox
!  01 Apr 2013 - M. Payer    - Restore IDTOX for tagged Ox simulation
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations; Removed tracers
!                              for old SOA simulation (H. Pye)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!

      ! for CTM tracers

      USE CMN_SIZE_MOD, ONLY : NNNTRID=>NNPAR

      IMPLICIT NONE

      INTEGER, PARAMETER :: MMMEMBER = 15
      INTEGER            :: NMEMBER(NNNTRID) 
      INTEGER            :: IDTRMB(NNNTRID,MMMEMBER)
      INTEGER            :: IDEMIS(NNNTRID)
      REAL*8             :: CTRMB(NNNTRID,MMMEMBER)
     
      ! ID's for SMVGEAR species
      INTEGER            :: IDO3,    IDNO2,    IDNO3,   IDNO,    IDHNO2
      INTEGER            :: IDN2O5,  IDHNO4,   IDHC1
      INTEGER            :: IDCO,    IDPRPE,   IDISOP,  IDALK4,  IDC3H8
      INTEGER            :: IDPAN,   IDGLPAN,  IDGPAN,  IDPMN,   IDPPN
      INTEGER            :: IDHNO3,  IDOH,     IDHO2,   IDH2O2,  IDACET
      INTEGER            :: IDMEK,   IDALD2,   IDRCHO,  IDMVK,   IDMACR
      INTEGER            :: IDISN2,  IDR4N2,   IDCH2O,  IDC2H6,  IDMP
      INTEGER            :: IDDMS,   IDSO2,    IDSO4,   IDMSA
      INTEGER            :: IDDRYO3, IDDRYPAN, IDDRYNO2, IDSO4s         
      ! IDs of aromatics not use (hotp 8/1/09)
!      INTEGER            :: IDGLYX,  IDMGLY
!      INTEGER            :: IDBENZ,  IDTOLU,   IDXYLE,  IDMONX
      INTEGER            :: IDDRYGLYX, IDDRYMGLY
      INTEGER            :: IDC2H2, IDC2H4
      INTEGER            :: IDMBO, IDGLYC
      INTEGER            :: IDHAC
      INTEGER            :: IDAPAN, IDENPAN, IDMPAN, IDNIPAN
      INTEGER            :: IDDRYAPAN, IDDRYENPAN, IDDRYGLPAN
      INTEGER            :: IDDRYGPAN, IDDRYMPAN, IDDRYNIPAN
      INTEGER            :: IDMPN
      ! jpp, 6/5/09, adding IDDRY__ id's for dry deposition
      INTEGER            :: IDDRYHOBr, IDDRYHBr, IDDRYBrNO3

      ! jpp, 6/5/09: bromine id's
      INTEGER            :: IDBr2,   IDBr,     IDBrO,  IDHOBr, IDHBr
      INTEGER            :: IDBrNO2, IDBrNO3,  IDCHBr3, IDCH2Br2
      INTEGER            :: IDCH3Br

      ! added for aromatic (dkh, 10/06/06)  
      INTEGER            :: IDLBRO2H, IDLBRO2N  
      INTEGER            :: IDLTRO2H, IDLTRO2N  
      INTEGER            :: IDLXRO2H, IDLXRO2N

      ! SOAupdate: added gas phase NAP chemistry (hotp 6/24/09)
      ! semivolpoa5: add IVOC chemistry (hotp 3/30/09)
      INTEGER            :: IDNAP,    IDNRO2,  IDLNRO2H, IDLNRO2N

      INTEGER            :: IDH2SO4   !(win, 6/23/09)
      !sfarina - leaving IDH2SO4 defined regardless of sim.

      ! GEOS-CHEM tracer ID's
      INTEGER            :: IDTNO,   IDTNO2,   IDTNO3,  IDTHNO2, IDTO3
      INTEGER            :: IDTPAN,  IDTCO,    IDTALK4
      INTEGER            :: IDTISOP, IDTHNO3,  IDTH2O2, IDTACET, IDTMEK
      INTEGER            :: IDTALD2, IDTRCHO,  IDTMVK,  IDTMACR, IDTPMN
      INTEGER            :: IDTPPN,  IDTISN2,  IDTR4N2, IDTPRPE, IDTC3H8
      INTEGER            :: IDTCH2O, IDTMP,    IDTN2O5, IDTHNO4, IDTC2H6
      INTEGER            :: IDTDMS,  IDTSO2,   IDTSO4,  IDTMSA,  IDTNH3
      INTEGER            :: IDTNH4,  IDTNIT,   IDTRN,   IDTPB,   IDTBE7
      INTEGER            :: IDTBCPI, IDTBCPO,  IDTOCPI, IDTOCPO, IDTDST1
      INTEGER            :: IDTDST2, IDTDST3,  IDTDST4, IDTSALA, IDTSALC
      INTEGER            :: IDTHG0,  IDTHg2,  IDTHgP,  IDTAS
      INTEGER            :: IDTAHS,  IDTLET,   IDTNH4aq,IDTSO4aq,IDTSO4s
      INTEGER            :: IDTNITs
      INTEGER            :: IDTBENZ, IDTTOLU,   IDTXYLE,  IDTMONX
      INTEGER            :: IDTGLYX, IDTMGLY
      INTEGER            :: IDTSOAG, IDTSOAM
      INTEGER            :: IDTC2H2, IDTC2H4
      INTEGER            :: IDTMBO,  IDTGLYC
      INTEGER            :: IDTAPAN, IDTENPAN, IDTMPAN, IDTNIPAN
      INTEGER            :: IDTGLPAN, IDTGPAN
      INTEGER            :: IDTHAC
      INTEGER            :: IDTISOPN
      INTEGER            :: IDTPROPNN
      INTEGER            :: IDTMOBA,  IDTMMN
      INTEGER            :: IDTRIP
      INTEGER            :: IDTIEPOX
      INTEGER            :: IDTPYPAN,IDTMAP
      INTEGER            :: IDTAP
      INTEGER            :: IDTMPN

      ! For tagged Ox simulation
      INTEGER            :: IDTO3Strt

      INTEGER            :: IDTHCOOH, IDTACTA  ! (mjc, 4/13/11)

      ! For POPs simulation (clf, 2/11/2011)
      INTEGER            :: IDTPOPPOC,  IDTPOPPBC,   IDTPOPG

      ! For H2/HD simulation
      INTEGER            :: IDTH2, IDTHD ! (hup, phs, 9/18/07)

      !jpp, 6/5/09: bromine ids
      INTEGER            :: IDTBr2,   IDTBr,   IDTBrO,  IDTHOBr, IDTHBr
      INTEGER            :: IDTBrNO2, IDTBrNO3, IDTCHBr3, IDTCH2Br2
      INTEGER            :: IDTCH3Br

      ! SOAudpate (hotp 6/15/09)
      ! semivolpoa: replace OCPO and OCPI with POA (hotp 2/17/09)
      INTEGER            :: IDTPOA1,  IDTPOA2
      ! semivolpoa2: add gas phase POA (hotp 2/27/09)
      INTEGER            :: IDTPOG1,  IDTPOG2
      ! semivolpoa4opoa: add OPOA,OPOG (hotp 3/17/09)
      INTEGER            :: IDTOPOA1, IDTOPOA2, IDTOPOG1, IDTOPOG2
      ! SOAupdate: added gas phase NAP chemistry (hotp 6/24/09)
      ! semivolpoa5: add IVOC chemistry (hotp 3/30/09)
      INTEGER            :: IDTNAP

      ! SOAudpate: all aromatic and IVOC SOA lumped together (hotp 5/12/10)
      ! ASOAN: non-volatile, low-NOx SOA product
      ! ASOA1-3: semi-volatile, high NOx SOA products
      ! LUMPAROMIVOC
      INTEGER            :: IDTASOAN 
      INTEGER            :: IDTASOA1, IDTASOG1
      INTEGER            :: IDTASOA2, IDTASOG2
      INTEGER            :: IDTASOA3, IDTASOG3

      ! SOAupdate: new mtp gas phase precursors (hotp 5/20/10)
      INTEGER            :: IDTMTPA, IDTLIMO, IDTMTPO
      ! new mtp semivolatile tracers (hotp 5/20/10)
      ! monoterpenes + sesquiterpenes (TSOA)
      INTEGER            :: IDTTSOA1, IDTTSOA2, IDTTSOA3
      INTEGER            :: IDTTSOG1, IDTTSOG2, IDTTSOG3
      ! C*=0.1 product (hotp 6/12/10)
      INTEGER            :: IDTTSOA0, IDTTSOG0
      ! isoprene (ISOA)
      INTEGER            :: IDTISOA1, IDTISOA2, IDTISOA3
      INTEGER            :: IDTISOG1, IDTISOG2, IDTISOG3

#if   defined( TOMAS )
      ! TOMAS tracer ID's  !(win, 6/23/09)
      INTEGER            :: IDTH2SO4

      INTEGER ::  IDTNK1   , IDTSF1   , IDTSS1   , IDTECOB1 , IDTECIL1 
      INTEGER ::  IDTOCOB1 , IDTOCIL1 , IDTDUST1 , IDTAW1   
      INTEGER ::  IDTNK2   , IDTSF2   , IDTSS2   , IDTECOB2 , IDTECIL2 
      INTEGER ::  IDTOCOB2 , IDTOCIL2 , IDTDUST2 , IDTAW2   
      INTEGER ::  IDTNK3   , IDTSF3   , IDTSS3   , IDTECOB3 , IDTECIL3 
      INTEGER ::  IDTOCOB3 , IDTOCIL3 , IDTDUST3 , IDTAW3   
      INTEGER ::  IDTNK4   , IDTSF4   , IDTSS4   , IDTECOB4 , IDTECIL4 
      INTEGER ::  IDTOCOB4 , IDTOCIL4 , IDTDUST4 , IDTAW4   
      INTEGER ::  IDTNK5   , IDTSF5   , IDTSS5   , IDTECOB5 , IDTECIL5 
      INTEGER ::  IDTOCOB5 , IDTOCIL5 , IDTDUST5 , IDTAW5   
      INTEGER ::  IDTNK6   , IDTSF6   , IDTSS6   , IDTECOB6 , IDTECIL6 
      INTEGER ::  IDTOCOB6 , IDTOCIL6 , IDTDUST6 , IDTAW6   
      INTEGER ::  IDTNK7   , IDTSF7   , IDTSS7   , IDTECOB7 , IDTECIL7 
      INTEGER ::  IDTOCOB7 , IDTOCIL7 , IDTDUST7 , IDTAW7   
      INTEGER ::  IDTNK8   , IDTSF8   , IDTSS8   , IDTECOB8 , IDTECIL8 
      INTEGER ::  IDTOCOB8 , IDTOCIL8 , IDTDUST8 , IDTAW8   
      INTEGER ::  IDTNK9   , IDTSF9   , IDTSS9   , IDTECOB9 , IDTECIL9 
      INTEGER ::  IDTOCOB9 , IDTOCIL9 , IDTDUST9 , IDTAW9   
      INTEGER ::  IDTNK10  , IDTSF10  , IDTSS10  , IDTECOB10, IDTECIL10
      INTEGER ::  IDTOCOB10, IDTOCIL10, IDTDUST10, IDTAW10  
      INTEGER ::  IDTNK11  , IDTSF11  , IDTSS11  , IDTECOB11, IDTECIL11
      INTEGER ::  IDTOCOB11, IDTOCIL11, IDTDUST11, IDTAW11  
      INTEGER ::  IDTNK12  , IDTSF12  , IDTSS12  , IDTECOB12, IDTECIL12
      INTEGER ::  IDTOCOB12, IDTOCIL12, IDTDUST12, IDTAW12  
# if  !defined( TOMAS12 )
      INTEGER ::  IDTNK13  , IDTSF13  , IDTSS13  , IDTECOB13, IDTECIL13
      INTEGER ::  IDTOCOB13, IDTOCIL13, IDTDUST13, IDTAW13  
      INTEGER ::  IDTNK14  , IDTSF14  , IDTSS14  , IDTECOB14, IDTECIL14
      INTEGER ::  IDTOCOB14, IDTOCIL14, IDTDUST14, IDTAW14  
      INTEGER ::  IDTNK15  , IDTSF15  , IDTSS15  , IDTECOB15, IDTECIL15
      INTEGER ::  IDTOCOB15, IDTOCIL15, IDTDUST15, IDTAW15  
#  if !defined( TOMAS15 )
      INTEGER ::  IDTNK16  , IDTSF16  , IDTSS16  , IDTECOB16, IDTECIL16
      INTEGER ::  IDTOCOB16, IDTOCIL16, IDTDUST16, IDTAW16  
      INTEGER ::  IDTNK17  , IDTSF17  , IDTSS17  , IDTECOB17, IDTECIL17
      INTEGER ::  IDTOCOB17, IDTOCIL17, IDTDUST17, IDTAW17  
      INTEGER ::  IDTNK18  , IDTSF18  , IDTSS18  , IDTECOB18, IDTECIL18
      INTEGER ::  IDTOCOB18, IDTOCIL18, IDTDUST18, IDTAW18  
      INTEGER ::  IDTNK19  , IDTSF19  , IDTSS19  , IDTECOB19, IDTECIL19
      INTEGER ::  IDTOCOB19, IDTOCIL19, IDTDUST19, IDTAW19  
      INTEGER ::  IDTNK20  , IDTSF20  , IDTSS20  , IDTECOB20, IDTECIL20
      INTEGER ::  IDTOCOB20, IDTOCIL20, IDTDUST20, IDTAW20  
      INTEGER ::  IDTNK21  , IDTSF21  , IDTSS21  , IDTECOB21, IDTECIL21
      INTEGER ::  IDTOCOB21, IDTOCIL21, IDTDUST21, IDTAW21  
      INTEGER ::  IDTNK22  , IDTSF22  , IDTSS22  , IDTECOB22, IDTECIL22
      INTEGER ::  IDTOCOB22, IDTOCIL22, IDTDUST22, IDTAW22  
      INTEGER ::  IDTNK23  , IDTSF23  , IDTSS23  , IDTECOB23, IDTECIL23
      INTEGER ::  IDTOCOB23, IDTOCIL23, IDTDUST23, IDTAW23  
      INTEGER ::  IDTNK24  , IDTSF24  , IDTSS24  , IDTECOB24, IDTECIL24
      INTEGER ::  IDTOCOB24, IDTOCIL24, IDTDUST24, IDTAW24  
      INTEGER ::  IDTNK25  , IDTSF25  , IDTSS25  , IDTECOB25, IDTECIL25
      INTEGER ::  IDTOCOB25, IDTOCIL25, IDTDUST25, IDTAW25  
      INTEGER ::  IDTNK26  , IDTSF26  , IDTSS26  , IDTECOB26, IDTECIL26
      INTEGER ::  IDTOCOB26, IDTOCIL26, IDTDUST26, IDTAW26  
      INTEGER ::  IDTNK27  , IDTSF27  , IDTSS27  , IDTECOB27, IDTECIL27
      INTEGER ::  IDTOCOB27, IDTOCIL27, IDTDUST27, IDTAW27  
      INTEGER ::  IDTNK28  , IDTSF28  , IDTSS28  , IDTECOB28, IDTECIL28
      INTEGER ::  IDTOCOB28, IDTOCIL28, IDTDUST28, IDTAW28  
      INTEGER ::  IDTNK29  , IDTSF29  , IDTSS29  , IDTECOB29, IDTECIL29
      INTEGER ::  IDTOCOB29, IDTOCIL29, IDTDUST29, IDTAW29  
      INTEGER ::  IDTNK30  , IDTSF30  , IDTSS30  , IDTECOB30, IDTECIL30
      INTEGER ::  IDTOCOB30, IDTOCIL30, IDTDUST30, IDTAW30  
#  endif
# endif
# if  defined( TOMAS40 )
      INTEGER ::  IDTNK31  , IDTSF31  , IDTSS31  , IDTECOB31, IDTECIL31
      INTEGER ::  IDTOCOB31, IDTOCIL31, IDTDUST31, IDTAW31  
      INTEGER ::  IDTNK32  , IDTSF32  , IDTSS32  , IDTECOB32, IDTECIL32
      INTEGER ::  IDTOCOB32, IDTOCIL32, IDTDUST32, IDTAW32  
      INTEGER ::  IDTNK33  , IDTSF33  , IDTSS33  , IDTECOB33, IDTECIL33
      INTEGER ::  IDTOCOB33, IDTOCIL33, IDTDUST33, IDTAW33  
      INTEGER ::  IDTNK34  , IDTSF34  , IDTSS34  , IDTECOB34, IDTECIL34
      INTEGER ::  IDTOCOB34, IDTOCIL34, IDTDUST34, IDTAW34  
      INTEGER ::  IDTNK35  , IDTSF35  , IDTSS35  , IDTECOB35, IDTECIL35
      INTEGER ::  IDTOCOB35, IDTOCIL35, IDTDUST35, IDTAW35  
      INTEGER ::  IDTNK36  , IDTSF36  , IDTSS36  , IDTECOB36, IDTECIL36
      INTEGER ::  IDTOCOB36, IDTOCIL36, IDTDUST36, IDTAW36  
      INTEGER ::  IDTNK37  , IDTSF37  , IDTSS37  , IDTECOB37, IDTECIL37
      INTEGER ::  IDTOCOB37, IDTOCIL37, IDTDUST37, IDTAW37  
      INTEGER ::  IDTNK38  , IDTSF38  , IDTSS38  , IDTECOB38, IDTECIL38
      INTEGER ::  IDTOCOB38, IDTOCIL38, IDTDUST38, IDTAW38  
      INTEGER ::  IDTNK39  , IDTSF39  , IDTSS39  , IDTECOB39, IDTECIL39
      INTEGER ::  IDTOCOB39, IDTOCIL39, IDTDUST39, IDTAW39  
      INTEGER ::  IDTNK40  , IDTSF40  , IDTSS40  , IDTECOB40, IDTECIL40
      INTEGER ::  IDTOCOB40, IDTOCIL40, IDTDUST40, IDTAW40  
# endif
#endif

      ! For tagged Hg simulation
      INTEGER              :: N_Hg_CATS
      INTEGER, ALLOCATABLE :: ID_Hg0(:),  ID_Hg2(:), ID_HgP(:)
!- eds 8/31/10 ------------------------------------------------------
!      INTEGER              :: ID_Hg_tot,  ID_Hg_na,  ID_Hg_eu 
!      INTEGER              :: ID_Hg_as,   ID_Hg_rw,  ID_Hg_oc
!      INTEGER              :: ID_Hg_ln,   ID_Hg_nt

      INTEGER              :: ID_Hg_tot,  ID_Hg_can,  ID_Hg_usa
      INTEGER              :: ID_Hg_cam,  ID_Hg_sam,  ID_Hg_waf
      INTEGER              :: ID_Hg_eaf,  ID_Hg_saf,  ID_Hg_naf
      INTEGER              :: ID_Hg_eur,  ID_Hg_eeu,  ID_Hg_sov
      INTEGER              :: ID_Hg_mde,  ID_Hg_sas,  ID_Hg_eas
      INTEGER              :: ID_Hg_sea,  ID_Hg_jpn,  ID_Hg_oce
      INTEGER              :: ID_Hg_so,   ID_Hg_bb,   ID_Hg_geo
      INTEGER              :: ID_Hg_atl,  ID_Hg_nat,  ID_Hg_sat
      INTEGER              :: ID_Hg_npa,  ID_Hg_arc,  ID_Hg_ant
      INTEGER              :: ID_Hg_ocn,  ID_Hg_str
!------------------------------------------------------------------------


      ! GEOS-CHEM emission ID's
      INTEGER            :: IDENO,   IDEO3,    IDECO,   IDEPRPE, IDEC3H8
      INTEGER            :: IDEALK4, IDEC2H6,  IDEISOP, IDEACET, IDEMEK
      INTEGER            :: IDEALD2, IDECH2O,  IDEHNO3
      INTEGER            :: NEMBIOG, NEMANTHRO
      INTEGER            :: IDEBENZ,  IDETOLU,   IDEXYLE,  IDEMONX
      INTEGER            :: IDEC2H2, IDEC2H4
      INTEGER            :: IDEMBO
      INTEGER            :: IDEGLYC
      INTEGER            :: IDEGLYX, IDEMGLY
      INTEGER            :: IDEHAC
      !%%%%(bmy 4/13/11)
      INTEGER            :: IDERCHO, IDEHCOOH, IDEACTA       
      !jpp, 8/24/07: emission ID's for Br biogenic emissions in fullchem runs
      INTEGER            :: IDECHBr3, IDECH2Br2, IDEBr2      
      ! SOAupdate: added gas phase NAP chemistry (hotp 6/24/09)
      ! semivolpoa5: add IVOC chemistry (hotp 3/30/09)
      INTEGER            :: IDENAP

      ! GEOS-CHEM biofuel ID's
      INTEGER            :: IDBFNO,   IDBFCO,   IDBFALK4, IDBFACET 
      INTEGER            :: IDBFMEK,  IDBFALD2, IDBFPRPE, IDBFC3H8
      INTEGER            :: IDBFCH2O, IDBFC2H6
      INTEGER            :: IDBFBENZ, IDBFTOLU, IDBFXYLE
      INTEGER            :: IDBFC2H2, IDBFC2H4
      INTEGER            :: IDBFGLYC
      INTEGER            :: IDBFGLYX, IDBFMGLY
      INTEGER            :: IDBFHAC
      ! SOAupdate: added gas phase NAP chemistry (hotp 6/24/09)
      ! semivolpoa5: add IVOC chemistry (hotp 3/30/09)
      INTEGER            :: IDBFNAP

      !BIOMASS BURNING  !(FP,01/27/10)
      INTEGER            :: IDBNO,   IDBCO,   IDBALK4, IDBACET 
      INTEGER            :: IDBMEK,  IDBALD2, IDBPRPE, IDBC3H8
      INTEGER            :: IDBCH2O, IDBC2H6
      INTEGER            :: IDBSO2,  IDBNH3
      INTEGER            :: IDBBC,   IDBOC
      INTEGER            :: IDBCO2
!kjw
      INTEGER            :: IDBCH4
!kjw

      INTEGER            :: IDBTOLU, IDBBENZ, IDBXYLE
      ! add dicarbonyl emissions (ccc, 01/27/10)
      INTEGER            :: IDBGLYX, IDBMGLY, IDBC2H4, IDBC2H2
      INTEGER            :: IDBGLYC, IDBHAC      
      ! add biomass burning emiss too (hotp 8/24/09)
      INTEGER            :: IDBNAP

      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE TRACERID( am_I_Root )
!
!******************************************************************************
!  Subroutine TRACERID reads the "tracer.dat" file and determines which
!  tracers, emission species, biomass burning species, and biofuel burning
!  species are turned on/off. (bmy, 3/16/01, 9/18/07)
!
!  NOTES:
!  (1 ) Original code from Loretta's version of the GISS-II model.  Now we
!        loop thru the tracer names and flag tracers that way. (bmy, 11/12/02)
!  (2 ) Added extra CASEs to the CASE statement for carbon & dust tracers.
!        (rjp, tdf, bmy, 4/1/04)
!  (3 ) Added extra CASEs to the CASE statement for seasalt tracers.
!        (rjp, bec, bmy, 4/20/04)
!  (4 ) Added extra CASEs to the CASE statement for SOA tracers.
!        (rjp, bmy, 7/13/04)
!  (5 ) Now references "tracer_mod.f".  NAME is now CHAR*14. (bmy, 7/20/04)
!  (6 ) Reverse the position of IDEISOP and IDECH2O so as to keep all of the
!        anthropogenic tracers together in IDEMS (bmy, 11/15/04)
!  (7 ) Added IDTHG0, IDTHG2, IDTHGP flags (eck, bmy, 12/7/04)
!  (8 ) Added IDTAS, IDTAHS, IDTLET, IDTNH4aq, IDTSO4aq.  Now no longer need 
!        to declare IDTCO, IDBCO, IDBFCO for offline aerosol simulations. 
!        (cas, bmy, 1/26/05)
!  (9 ) Added IDTSO4s and IDTNITs (bec, bmy, 4/13/05)
!  (10) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (11) Add alternate names for tagged Hg tracers.  Also define ocean mercury 
!        flux categories.  Now references LSPLIT from "logical_mod.f".
!        (cdh, bmy, 12/15/05)
!  (12) Now remove IDBxxx biomass flags (bmy, 4/5/06)
!  (13) Now look for IDTSOG4 and IDTSOA4 (bmy, 5/18/06)
!  (14) Minor fixes for CH3I simulation (bmy, 7/25/06)
!  (15) Now define IDTH2, IDTHD (hup, lyj, phs, 9/18/07)
!  (16) To satisfy IF statement in EMISSDR for using EMFOSSIL, we need 
!        to set IDECO=1 instead of IDECO=2. (jaf, mak, bmy, 2/14/08)
!  (17) Increase NEMANTHRO from 10 to 12 and set IDEOX and IDEHNO3 (phs, 3/4/08)
!  (18) Added CASE matching for TOMAS tracers (win, 6/23/09)
!  (19) Add online definition of IDBs. (fp, hotp, 01/10)
!  (20) Add definition of IDTSOG5, IDTSOA5, IDTISOPN, PROPNN, 
!        IDTAP, IDTMOBA, IDTMMN, IDTRIP, IDTIEPOX, IDTMAP, IDTPYPAN
!        (fp, hotp, 01/10)
!  (20) Change hard-wired IDEs to dynamically defined IDEs. (fp, hotp, 01/10)
!  (21) Add IDEMS definitions for new species (fp, hotp, 01/10)
!  (22) Add writing check on IDs. (hotp, 01/10)
!  (23) Added IDTPOPG, IDTPOPPOC, and IDTPOPPBC for POPs simulation (eck, 9/20/10)
!  17 Aug 2011 - R. Yantosca - Print values of 3 extra emission IDExxx flags
!  22 Aug 2011 - R. Yantosca - Add modifications for APM microphysics (G. Luo)
!  16 Feb 2012 - M. Payer    - Account for alternate spellings of the Rn and Pb
!                              tracers in input.geos
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  07 Feb 2013 - M. Payer    - Add NO, NO2, NO2, HNO2, and O3 as tracers to
!                              replace NOx and Ox
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations; Removed tracers
!                              for old SOA simulation (H. Pye)
!******************************************************************************
!
      ! References to F90 modules
      USE CHARPAK_MOD, ONLY : TRANUC
      USE LOGICAL_MOD, ONLY : LSPLIT
      USE TRACER_MOD,  ONLY : ITS_A_C2H6_SIM,  ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,  ONLY : ITS_A_TAGCO_SIM, ITS_A_MERCURY_SIM
      USE TRACER_MOD,  ONLY : ITS_A_POPS_SIM
      USE TRACER_MOD,  ONLY : ITS_A_H2HD_SIM
      USE TRACER_MOD,  ONLY : N_TRACERS,       TRACER_NAME
      USE ERROR_MOD,   ONLY : ERROR_STOP

      USE CMN_SIZE_MOD                    ! Size parameters
      USE COMODE_LOOP_MOD                 ! IDEMS
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN)  :: am_I_Root   ! Is this the root CPU?

      ! Local variables
      INTEGER              :: N, COUNT, COUNT_Hg0, COUNT_Hg2, COUNT_HgP
      CHARACTER(LEN=14)    :: NAME

      !(fp, 6/09)
      ! COUNT_BB is the number of BB tracers for a given sim
      INTEGER              :: COUNT_BB
      
      !=================================================================
      ! TRACERID begins here!
      !
      ! NOTE: There are still some vestiges of historical baggage, we
      !       will get rid of this as time allows (bmy, 11/12/02)
      !=================================================================

      ! Array dimension check (hotp 8/1/09)
      ! NNNTRID is used in this module to dimension some arrays
      ! NNNTRID should be the same as NNPAR in CMN_SIZE
      IF ( NNNTRID .NE. NNPAR ) THEN
#if   defined( APM )
         write(*,*)'NNNTRID .NE. NNPAR!',NNNTRID,NNPAR
#else
         CALL ERROR_STOP ( 'NNNTRID .NE. NNPAR!',
     &                                     'tracerid_mod,CMN_SIZE' )
#endif
         ! If you get this error, an appropriate fix is to increase 
         ! either NNNTRID (in tracerid_mod) or NNPAR (in CMN_SIZE) so
         ! they are the same and greater than or equal to the number 
         ! of tracers in your simulation (hotp 8/1/09)
      ENDIF

      ! Zero all ID #'s and allocate Hg index arrays (if necessary)
      CALL INIT_TRACERID
      
      ! Initialize counters
      COUNT     = 0
      COUNT_Hg0 = 0
      COUNT_Hg2 = 0
      COUNT_HgP = 0

      !(fp, 6/09)
      !Initialize COUNT_BB
      COUNT_BB = 0

      !=================================================================
      ! Assign tracer, biomass, biofuel, and anthro emission ID's
      !=================================================================
      DO N = 1, N_TRACERS

         ! Save TRACER_NAME into a CHARACTER*14 variable
         NAME = TRACER_NAME(N)

         ! Convert name to uppercase
         CALL TRANUC( NAME )

         ! Find each tracer
         SELECT CASE ( TRIM( NAME ) )
        
            !------------------------
            ! Full chem tracers
            !------------------------
            CASE ( 'NO' )
               COUNT    = COUNT + 1
               IDTNO    = N
               IDBFNO   = COUNT

               !(fp)              
               COUNT_BB = COUNT_BB + 1
               IDBNO    = COUNT_BB

            CASE ( 'NO2' )
               IDTNO2   = N

            CASE ( 'NO3' )
               IDTNO3   = N

            CASE ( 'HNO2' )
               IDTHNO2  = N

            CASE ( 'O3' )
               IDTO3    = N

            CASE ( 'PAN' )
               IDTPAN   = N
               
            CASE ( 'CO' )
               COUNT    = COUNT + 1
               IDTCO    = N
               IDBFCO   = COUNT

               !(fp)             
               COUNT_BB = COUNT_BB + 1
               IDBCO   = COUNT_BB

               ! Special case: Tagged CO
               ! Set some emission flags and then exit
               ! NOTE: To satisfy IF statement in EMISSDR for using 
               ! EMFOSSIL, we need to set IDECO=1 instead of IDECO=2.
               ! (jaf, mak, bmy, 2/14/08)
               IF ( ITS_A_TAGCO_SIM() ) THEN 
                  NEMANTHRO = 1
                  IDECO     = 1
                  IDTISOP   = 1
                  EXIT
               ENDIF

            !-----------------------------------
            ! FEW ASSUMPTIONS FOR H2HD SIM:
            ! IDTH2=1, IDTHD=2, IDTCO=N('H2')
            ! H2/HD simulation requires CO...
            ! (hup, lyj, phs, 9/18/07)
            !-----------------------------------
            CASE ( 'H2' )
               COUNT    = COUNT + 1
               IDTCO    = N
               IDBFCO   = COUNT

               ! Special case: Tagged H2 (hup 4/28/2004)
               ! Set some emissions flags then exit
               IF ( ITS_A_H2HD_SIM() ) THEN
                  NEMANTHRO = 1
                  IDECO     = 2
                  IDTISOP   = 1
                  IDTH2     = 1 ! (hup 7/14/2004)

               ENDIF

            ! ... and HD
            CASE ( 'HD' )
               COUNT    = COUNT + 1
               IDTHD    = N

            CASE ( 'ALK4' )
               COUNT    = COUNT + 1
               IDTALK4  = N
               IDBFALK4 = COUNT

               !(fp)              
               COUNT_BB = COUNT_BB + 1
               IDBALK4  = COUNT_BB

            CASE ( 'ISOP' )
               IDTISOP  = N
               
            CASE ( 'HNO3' )
               IDTHNO3  = N

            CASE ( 'H2O2' )
               IDTH2O2  = N

            CASE ( 'ACET' )
               COUNT    = COUNT + 1
               IDTACET  = N
               IDBFACET = COUNT

               !(fp)              
               COUNT_BB = COUNT_BB + 1
               IDBACET  = COUNT_BB

            CASE ( 'MEK' )
               COUNT    = COUNT + 1
               IDTMEK   = N
               IDBFMEK  = COUNT

               !(fp)              
               COUNT_BB = COUNT_BB + 1
               IDBMEK   = COUNT_BB

            CASE ( 'ALD2' )
               COUNT    = COUNT + 1
               IDTALD2  = N
               IDBFALD2 = COUNT

               !(fp)              
               COUNT_BB = COUNT_BB + 1
               IDBALD2   = COUNT_BB

            CASE ( 'RCHO' )
               IDTRCHO  = N

            CASE ( 'MVK' )
               IDTMVK   = N

            CASE ( 'MACR' )
               IDTMACR  = N

            CASE ( 'PMN' )
               IDTPMN   = N

            CASE ( 'PPN' )
               IDTPPN   = N

            CASE ( 'R4N2' )
               IDTR4N2  = N

            CASE ( 'PRPE' )
               COUNT    = COUNT + 1
               IDTPRPE  = N
               IDBFPRPE = COUNT

               !(fp)              
               COUNT_BB = COUNT_BB + 1
               IDBPRPE  = COUNT_BB

            CASE ( 'C3H8' )
               COUNT    = COUNT + 1
               IDTC3H8  = N
               IDBFC3H8 = COUNT

               !(fp)              
               COUNT_BB = COUNT_BB + 1
               IDBC3H8  = COUNT_BB

            CASE ( 'CH2O' )
               COUNT    = COUNT + 1
               IDTCH2O  = N
               IDBFCH2O = COUNT

               !(fp)              
               COUNT_BB = COUNT_BB + 1
               IDBCH2O  = COUNT_BB

            CASE ( 'C2H6' )
               COUNT    = COUNT + 1
               IDTC2H6  = N
               IDBFC2H6 = COUNT
               !(fp)              
               COUNT_BB = COUNT_BB + 1
               IDBC2H6  = COUNT_BB

               ! Special case: tagged C2H6
               ! Set emission flags and then exit
               IF ( ITS_A_C2H6_SIM() ) THEN
                  NEMANTHRO = 1
                  IDEC2H6   = 1
                  EXIT
               ENDIF

            CASE ( 'N2O5' )
               IDTN2O5  = N

            CASE ( 'HNO4' )
               IDTHNO4  = N

            CASE ( 'MPN' )
               IDTMPN  = N

            CASE ( 'MP' )
               IDTMP    = N

            !--------------------------------
            ! Sulfur & nitrate aerosols
            !--------------------------------
            CASE ( 'DMS' )
               IDTDMS   = N

            CASE ( 'SO2' )
               IDTSO2   = N

               !(fp)              
               COUNT_BB = COUNT_BB + 1
               IDBSO2   = COUNT_BB

            CASE ( 'SO4' )
               IDTSO4   = N

            CASE ( 'SO4S' )
               IDTSO4s  = N

            CASE ( 'MSA' )
               IDTMSA   = N

            CASE ( 'NH3' )
               IDTNH3   = N
               !(fp)              
               COUNT_BB = COUNT_BB + 1
               IDBNH3   = COUNT_BB

            CASE ( 'NH4' )
               IDTNH4   = N

            CASE ( 'NIT' )
               IDTNIT   = N

            CASE ( 'NITS' )
               IDTNITs  = N

            !--------------------------------
            ! Crystalline & aqueous aerosols
            !--------------------------------
            CASE ( 'AS' ) 
               IDTAS    = N

            CASE ( 'AHS' ) 
               IDTAHS   = N

            CASE ( 'LET' )
               IDTLET   = N

            CASE ( 'NH4AQ' )
               IDTNH4aq = N
              
            CASE ( 'SO4AQ' )
               IDTSO4aq = N
             
            !--------------------------------
            ! Carbon & 2dy organic aerosols
            !--------------------------------
            CASE ( 'BCPI' )
               IDTBCPI  = N

               !(fp) WEIRD BUT TO BE CONSISTENT CHECK
               COUNT_BB = COUNT_BB + 1
               IDBBC    = COUNT_BB

            CASE ( 'OCPI' )
               IDTOCPI  = N
               !(fp) WEIRD BUT TO BE CONSISTENT CHECK
               COUNT_BB = COUNT_BB + 1
               IDBOC    = COUNT_BB

            CASE ( 'BCPO' )
               IDTBCPO  = N
    
            CASE ( 'OCPO' )
               IDTOCPO  = N

            ! SOAupdate (hotp 6/15/09)
            ! semivolpoa: add POA (hotp 2/17/09)
            CASE ( 'POA1' )
               IDTPOA1   = N
               COUNT_BB = COUNT_BB + 1
               IDBOC    = COUNT_BB

            CASE ( 'POA2' )
               IDTPOA2   = N

            ! SOAudpate: semivol2: add POG (hotp 2/27/09)
            CASE ( 'POG1' )
               IDTPOG1   = N

            CASE ( 'POG2' )
               IDTPOG2   = N

            ! SOAupdate: semivolpoa4opoa: add OPOA OPOG (hotp 3/17/09)
            CASE ( 'OPOA1' )
               IDTOPOA1  = N

            CASE ( 'OPOA2' )
               IDTOPOA2  = N

            CASE ( 'OPOG1' )
               IDTOPOG1  = N

            CASE ( 'OPOG2' )
               IDTOPOG2  = N

            ! SOAupdate: remove old and add new mtp (hotp 5/20/10)
            ! Gas-phase precursors (hotp 5/20/10)
            CASE ( 'MTPA' )
               IDTMTPA  = N

            CASE ( 'LIMO' )
               IDTLIMO  = N

            CASE ( 'MTPO' )
               IDTMTPO  = N

            ! new mtp gas-phase semivols (hotp 5/20/10)
            CASE ( 'TSOG1' )
               IDTTSOG1  = N

            CASE ( 'TSOG2' )
               IDTTSOG2  = N

            CASE ( 'TSOG3' )
               IDTTSOG3  = N

            CASE ( 'TSOG0' )
               IDTTSOG0  = N

            ! new mtp aerosol-phase semivols (hotp 5/20/10)
            CASE ( 'TSOA1' )
               IDTTSOA1  = N

            CASE ( 'TSOA2' )
               IDTTSOA2  = N

            CASE ( 'TSOA3' )
               IDTTSOA3  = N

            CASE ( 'TSOA0' )
               IDTTSOA0  = N

            ! new isoprene semivols
            CASE ( 'ISOG1' )
               IDTISOG1  = N

            CASE ( 'ISOG2' )
               IDTISOG2  = N

            CASE ( 'ISOG3' )
               IDTISOG3  = N

            CASE ( 'ISOA1' )
               IDTISOA1  = N

            CASE ( 'ISOA2' )
               IDTISOA2  = N

            CASE ( 'ISOA3' )
               IDTISOA3  = N

            ! replace SOG5 and SOG6 with lumped species (hotp 5/12/10)
            CASE ( 'ASOG1' )
               IDTASOG1 = N

            CASE ( 'ASOG2' )
               IDTASOG2 = N

            CASE ( 'ASOG3' )
               IDTASOG3 = N

            CASE ( 'ASOAN' )
               IDTASOAN = N

            CASE ( 'ASOA1' )
               IDTASOA1 = N

            CASE ( 'ASOA2' )
               IDTASOA2 = N

            CASE ( 'ASOA3' )
               IDTASOA3 = N
            ! end hotp SOAupdate

            !--------------------------------
            ! Mineral dust aerosols
            !--------------------------------
            CASE ( 'DST1' )
               IDTDST1  = N

            CASE ( 'DST2' )
               IDTDST2  = N
        
            CASE ( 'DST3' )
               IDTDST3  = N

            CASE ( 'DST4' )
               IDTDST4  = N

            !--------------------------------
            ! Seasalt aerosols
            !--------------------------------
            CASE ( 'SALA' )
               IDTSALA  = N

            CASE ( 'SALC' )
               IDTSALC  = N

            !--------------------------------
            ! Bromine (jpp, 6/13/07)
            !--------------------------------

            CASE ( 'BR2' )
               IDTBr2     = N
            
            CASE ( 'BR' )
               IDTBr      = N

            CASE ( 'BRO' )
               IDTBrO     = N
            
            CASE ( 'HOBR' )
               IDTHOBr    = N

            CASE ( 'HBR' )
               IDTHBr     = N

            CASE ( 'BRNO2' )
               IDTBrNO2   = N
                           
            CASE ( 'BRNO3' )
               IDTBrNO3   = N

            CASE ( 'CHBR3' )
               IDTCHBr3   = N

            CASE ( 'CH2BR2' )
               IDTCH2Br2   = N

            CASE ( 'CH3BR' )
               IDTCH3Br   = N

            !--------------------------------
            ! Dicarbonyls GLYX & MGLY
            !--------------------------------
            CASE ( 'GLYX' )
               COUNT    = COUNT + 1
               IDTGLYX  = N
               IDBFGLYX = COUNT

               COUNT_BB  = COUNT_BB + 1
               IDBGLYX   = COUNT_BB

            CASE ( 'MGLY' )
               COUNT    = COUNT + 1
               IDTMGLY  = N
               IDBFMGLY = COUNT

               COUNT_BB  = COUNT_BB + 1
               IDBMGLY   = COUNT_BB

            !--------------------------------
            ! Aromatics tracers
            !--------------------------------
            CASE ( 'BENZ' )
               COUNT    = COUNT + 1
               IDTBENZ  = N
               IDBFBENZ = COUNT

               COUNT_BB = COUNT_BB + 1
               IDBBENZ  = COUNT_BB

            CASE ( 'TOLU' )
               COUNT    = COUNT + 1
               IDTTOLU  = N
               IDBFTOLU = COUNT

               COUNT_BB = COUNT_BB + 1
               IDBTOLU  = COUNT_BB
               
            CASE ( 'XYLE' )
               COUNT    = COUNT + 1
               IDTXYLE  = N
               IDBFXYLE = COUNT

               COUNT_BB = COUNT_BB + 1
               IDBXYLE  = COUNT_BB

            ! SOAupdate: for gas phase NAP chemistry (hotp 6/24/09)
            !-------------------------------
            ! semivolpoa5: IVOC tracer (hotp 3/30/09)
            !-------------------------------
            CASE ( 'NAP' )
               COUNT    = COUNT + 1
               IDTNAP   = N
               IDBFNAP  = COUNT
               COUNT_BB = COUNT_BB + 1
               IDBNAP   = COUNT_BB

            !--------------------------------
            ! Monoterpene
            !--------------------------------
            CASE ( 'MONX' )
               IDTMONX  = N

            !--------------------------------
            ! SOA from GLYX and MGLY
            !--------------------------------
            CASE ( 'SOAG' )
               IDTSOAG  = N

            CASE ( 'SOAM' )
               IDTSOAM  = N

            !--------------------------------
            ! C2H4
            !--------------------------------
            CASE ( 'C2H4' )
               COUNT    = COUNT + 1
               IDTC2H4  = N
               IDBFC2H4 = COUNT

               COUNT_BB = COUNT_BB + 1
               IDBC2H4  = COUNT_BB
               
            !--------------------------------
            ! C2H2
            !--------------------------------
            CASE ( 'C2H2' )
               COUNT    = COUNT + 1
               IDTC2H2  = N
               IDBFC2H2 = COUNT

               COUNT_BB = COUNT_BB + 1
               IDBC2H2  = COUNT_BB
               
            !--------------------------------
            ! MBO
            !--------------------------------
            CASE ( 'MBO' )
               IDTMBO   = N

            !--------------------------------
            ! GLYC
            !--------------------------------
            CASE ( 'GLYC' )
               COUNT    = COUNT + 1
               IDTGLYC  = N
               IDBFGLYC = COUNT

               COUNT_BB  = COUNT_BB + 1
               IDBGLYC   = COUNT_BB

            !--------------------------------
            ! HAC
            !--------------------------------
            CASE ( 'HAC' )
               COUNT    = COUNT + 1
               IDTHAC   = N
               IDBFHAC  = COUNT

               COUNT_BB  = COUNT_BB + 1
               IDBHAC    = COUNT_BB

            !--------------------------------
            ! new PAN species
            !--------------------------------
            CASE ( 'APAN' )
               IDTAPAN   = N

            CASE ( 'ENPAN' )
               IDTENPAN   = N

            CASE ( 'GLPAN' )
               IDTGLPAN   = N

            CASE ( 'GPAN' )
               IDTGPAN   = N

            CASE ( 'MPAN' )
               IDTMPAN   = N

            CASE ( 'NIPAN' )
               IDTNIPAN   = N

            !
            !(fp, 6/09)
            !

            CASE ( 'ISOPN' )
               IDTISOPN  = N

            CASE ( 'PROPNN' )
               IDTPROPNN = N

            CASE ( 'AP' )
               IDTAP     = N

            CASE (  'MOBA' )
               IDTMOBA   = N

            CASE (  'MMN' )
               IDTMMN    = N

            CASE (  'RIP'  )
               IDTRIP    = N

            CASE (  'IEPOX'  )
               IDTIEPOX  = N

            CASE (  'MAP'  )
               IDTMAP    = N

            CASE (  'PYPAN'  )
               IDTPYPAN  = N

            !--------------------------------
            ! Tagged Ox Tracers
            !--------------------------------
             
            CASE ( 'O3STRT', 'O3STRAT' )
               IDTO3Strt = N

            !--------------------------------
            ! Rn-Pb-Be tracers
            !--------------------------------
            CASE ( 'RN', 'RN222' )
               IDTRN    = N
 
            CASE ( 'PB', 'PB210' )
               IDTPB    = N

            CASE ( 'BE7' )
               IDTBE7   = N

            !--------------------------------
            ! CH3I and HCN tracers
            !--------------------------------

            ! Special case: CH3I needs CO biomass/biofuel
            CASE ( 'CH3I', 'CH3IOC' )
               COUNT    = COUNT + 1
               IDTCO    = 1
               IDBFCO   = COUNT
               COUNT_BB = COUNT_BB + 1
               IDBCO    = COUNT_BB 
               NEMANTHRO= 8      ! Reset NEMANTHRO here too (bmy, 7/25/06)
               EXIT

            ! Special case: HCN needs CO biomass/biofuel
            CASE ( 'HCN' )
               COUNT    = COUNT + 1
               IDTCO    = 1
               IDBFCO   = COUNT
               COUNT_BB = COUNT_BB + 1
               IDBCO    = COUNT_BB 
               EXIT

            !--------------------------------
            ! Total & tagged mercury tracers 
            ! (eck, cdh, bmy, 12/15/05)
            !--------------------------------
            CASE ( 'HG0' )
               COUNT_Hg0         = COUNT_Hg0 + 1              
               ID_Hg_tot         = COUNT_Hg0
               IDTHg0            = N
               ID_Hg0(COUNT_Hg0) = N

               ! Special case: Hg BB emissions are scaled from CO:
               ! (ccc, 5/21/10)
               COUNT_BB          = COUNT_BB + 1
               IDBCO             = COUNT_BB 

            CASE ( 'HG2' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HGP' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

!- eds 8/31/10----------------------------------------------------
!            CASE ( 'HG0_AN_NA', 'HG0_AN' )
!               COUNT_Hg0         = COUNT_Hg0 + 1
!               ID_Hg_na          = COUNT_Hg0
!               ID_Hg0(COUNT_Hg0) = N
!
!            CASE ( 'HG0_AN_EU', 'HG0_AE' )
!               COUNT_Hg0         = COUNT_Hg0 + 1
!               ID_Hg_eu          = COUNT_Hg0
!               ID_Hg0(COUNT_Hg0) = N
!
!            CASE ( 'HG0_AN_AS', 'HG0_AA' )
!               COUNT_Hg0         = COUNT_Hg0 + 1
!               ID_Hg_as          = COUNT_Hg0
!               ID_Hg0(COUNT_Hg0) = N
!
!            CASE ( 'HG0_AN_RW', 'HG0_AR' )
!               COUNT_Hg0         = COUNT_Hg0 + 1
!               ID_Hg_rw          = COUNT_Hg0
!               ID_Hg0(COUNT_Hg0) = N
!
!            CASE ( 'HG0_OC' )
!               COUNT_Hg0         = COUNT_Hg0 + 1
!               ID_Hg_oc          = COUNT_Hg0
!               ID_Hg0(COUNT_Hg0) = N
!
!            CASE ( 'HG0_LN' )
!               COUNT_Hg0         = COUNT_Hg0 + 1
!               ID_Hg_ln          = COUNT_Hg0
!               ID_Hg0(COUNT_Hg0) = N
!
!            CASE ( 'HG0_NT' )
!               COUNT_Hg0         = COUNT_Hg0 + 1
!               ID_Hg_nt          = COUNT_Hg0
!               ID_Hg0(COUNT_Hg0) = N
!
!            CASE ( 'HG2_AN_NA', 'HG2_AN' )
!               COUNT_Hg2         = COUNT_Hg2 + 1
!               ID_Hg2(COUNT_Hg2) = N
!
!            CASE ( 'HG2_AN_EU', 'HG2_AE' )
!               COUNT_Hg2         = COUNT_Hg2 + 1
!               ID_Hg2(COUNT_Hg2) = N
!
!            CASE ( 'HG2_AN_AS', 'HG2_AA' )
!               COUNT_Hg2         = COUNT_Hg2 + 1
!               ID_Hg2(COUNT_Hg2) = N
!
!            CASE ( 'HG2_AN_RW', 'HG2_AR' )
!               COUNT_Hg2         = COUNT_Hg2 + 1
!               ID_Hg2(COUNT_Hg2) = N
!
!            CASE ( 'HG2_OC' )
!               COUNT_Hg2         = COUNT_Hg2 + 1
!               ID_Hg2(COUNT_Hg2) = N
!
!            CASE ( 'HG2_LN' )
!               COUNT_Hg2         = COUNT_Hg2 + 1
!               ID_Hg2(COUNT_Hg2) = N
!
!            CASE ( 'HG2_NT' )
!               COUNT_Hg2         = COUNT_Hg2 + 1
!               ID_Hg2(COUNT_Hg2) = N
!
!            CASE ( 'HGP_AN_NA', 'HGP_AN' )
!               COUNT_HgP         = COUNT_HgP + 1
!               ID_HgP(COUNT_HgP) = N
!
!            CASE ( 'HGP_AN_EU', 'HGP_AE' )
!               COUNT_HgP         = COUNT_HgP + 1
!               ID_HgP(COUNT_HgP) = N
!
!            CASE ( 'HGP_AN_AS', 'HGP_AA' )
!               COUNT_HgP         = COUNT_HgP + 1
!               ID_HgP(COUNT_HgP) = N
!
!            CASE ( 'HGP_AN_RW', 'HGP_AR' )
!               COUNT_HgP         = COUNT_HgP + 1
!               ID_HgP(COUNT_HgP) = N
!
!            CASE ( 'HGP_OC' )
!               COUNT_HgP         = COUNT_HgP + 1
!               ID_HgP(COUNT_HgP) = N
!
!            CASE ( 'HGP_LN' )
!               COUNT_HgP         = COUNT_HgP + 1
!               ID_HgP(COUNT_HgP) = N
!
!            CASE ( 'HGP_NT' )
!               COUNT_HgP         = COUNT_HgP + 1
!               ID_HgP(COUNT_HgP) = N

            ! New Tagged Hg Simulation (eds 8/31/10)
            CASE ( 'HG0_AN_CAN', 'HG0_CAN' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_can          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_AN_USA', 'HG0_USA' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_usa          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_AN_CAM', 'HG0_CAM' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_cam          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_AN_SAM', 'HG0_SAM' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_sam          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_AN_WAF', 'HG0_WAF' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_waf          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_AN_EAF', 'HG0_EAF' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_eaf          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_AN_SAF', 'HG0_SAF' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_saf          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_AN_NAF', 'HG0_NAF' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_naf          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_AN_EUR', 'HG0_EUR' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_eur          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_AN_EEU', 'HG0_EEU' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_eeu          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_AN_MDE', 'HG0_MDE' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_mde          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_AN_SOV', 'HG0_SOV' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_sov          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_AN_SAS', 'HG0_SAS' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_sas          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_AN_EAS', 'HG0_EAS' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_eas          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_AN_SEA', 'HG0_SEA' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_sea          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_AN_JPN', 'HG0_JPN' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_jpn          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_AN_OCE', 'HG0_OCE' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_oce          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_SO' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_so          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_BB' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_bb          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_GEO' )
               COUNT_Hg0          = COUNT_Hg0 + 1
               ID_Hg_geo          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_ATL' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_atl         = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_NAT' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_nat         = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_SAT' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_sat         = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_NPA' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_npa         = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_ARC' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_arc         = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_ANT' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_ant         = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_OCN' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_ocn         = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_STR' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_str         = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG2_AN_CAN', 'HG2_CAN' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_AN_USA', 'HG2_USA' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_AN_CAM', 'HG2_CAM' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_AN_SAM', 'HG2_SAM' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_AN_WAF', 'HG2_WAF' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_AN_EAF', 'HG2_EAF' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_AN_SAF', 'HG2_SAF' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_AN_NAF', 'HG2_NAF' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_AN_EUR', 'HG2_EUR' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_AN_EEU', 'HG2_EEU' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_AN_MDE', 'HG2_MDE' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_AN_SOV', 'HG2_SOV' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_AN_SAS', 'HG2_SAS' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_AN_EAS', 'HG2_EAS' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_AN_SEA', 'HG2_SEA' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_AN_JPN', 'HG2_JPN' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_AN_OCE', 'HG2_OCE' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_SO' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_BB' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_GEO' )
               COUNT_Hg2          = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_ATL' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_NAT' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_SAT' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_NPA' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_ARC' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_ANT' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_OCN' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_STR' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HGP_AN_CAN', 'HGP_CAN' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_AN_USA', 'HGP_USA' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_AN_CAM', 'HGP_CAM' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_AN_SAM', 'HGP_SAM' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_AN_WAF', 'HGP_WAF' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_AN_EAF', 'HGP_EAF' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_AN_SAF', 'HGP_SAF' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_AN_NAF', 'HGP_NAF' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_AN_EUR', 'HGP_EUR' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_AN_EEU', 'HGP_EEU' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_AN_MDE', 'HGP_MDE' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_AN_SOV', 'HGP_SOV' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_AN_SAS', 'HGP_SAS' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_AN_EAS', 'HGP_EAS' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_AN_SEA', 'HGP_SEA' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_AN_JPN', 'HGP_JPN' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_AN_OCE', 'HGP_OCE' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_SO' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_BB' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_GEO' )
               COUNT_HgP          = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_ATL' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_NAT' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_SAT' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_NPA' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_ARC' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_ANT' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_OCN' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_STR' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N
!--------------------------------------------------------------------------------

            !--------------------------------
            ! POPs
            ! (eck, 9/20/10); (clf, 2/11/2011)
            !--------------------------------
            CASE ( 'POPG' )
               COUNT     = COUNT+1              
               IDTPOPG         = N
               
            CASE ( 'POPPOC' )
               COUNT     = COUNT+1
               IDTPOPPOC       = N

            CASE ( 'POPPBC' )
               COUNT     = COUNT+1
               IDTPOPPBC       = N

#if   defined( TOMAS )
            !-----------------------------------
            ! TOMAS Aerosol microphysics tracers (win, 6/23/09)
            !-----------------------------------
               
            CASE ('H2SO4')
                IDTH2SO4   = N

            CASE (  'NK1')
                  IDTNK1   = N
            CASE (  'NK2')
                  IDTNK2   = N
            CASE (  'NK3')
                  IDTNK3   = N
            CASE (  'NK4')
                  IDTNK4   = N
            CASE (  'NK5')
                  IDTNK5   = N
            CASE (  'NK6')
                  IDTNK6   = N
            CASE (  'NK7')
                  IDTNK7   = N
            CASE (  'NK8')
                  IDTNK8   = N
            CASE (  'NK9')
                  IDTNK9   = N
            CASE (  'NK10')
                  IDTNK10   = N
            CASE (  'NK11')
                  IDTNK11   = N
            CASE (  'NK12')
                  IDTNK12   = N
# if  !defined( TOMAS12 )
            CASE (  'NK13')
                  IDTNK13   = N
            CASE (  'NK14')
                  IDTNK14   = N
            CASE (  'NK15')
                  IDTNK15   = N
#  if !defined( TOMAS15 )
            CASE (  'NK16')
                  IDTNK16   = N
            CASE (  'NK17')
                  IDTNK17   = N
            CASE (  'NK18')
                  IDTNK18   = N
            CASE (  'NK19')
                  IDTNK19   = N
            CASE (  'NK20')
                  IDTNK20   = N
            CASE (  'NK21')
                  IDTNK21   = N
            CASE (  'NK22')
                  IDTNK22   = N
            CASE (  'NK23')
                  IDTNK23   = N
            CASE (  'NK24')
                  IDTNK24   = N
            CASE (  'NK25')
                  IDTNK25   = N
            CASE (  'NK26')
                  IDTNK26   = N
            CASE (  'NK27')
                  IDTNK27   = N
            CASE (  'NK28')
                  IDTNK28   = N
            CASE (  'NK29')
                  IDTNK29   = N
            CASE (  'NK30')
                  IDTNK30   = N
#  endif
# endif
# if   defined( TOMAS40 )
            CASE (  'NK31')
                  IDTNK31   = N
            CASE (  'NK32')
                  IDTNK32   = N
            CASE (  'NK33')
                  IDTNK33   = N
            CASE (  'NK34')
                  IDTNK34   = N
            CASE (  'NK35')
                  IDTNK35   = N
            CASE (  'NK36')
                  IDTNK36   = N
            CASE (  'NK37')
                  IDTNK37   = N
            CASE (  'NK38')
                  IDTNK38   = N
            CASE (  'NK39')
                  IDTNK39   = N
            CASE (  'NK40')
                  IDTNK40   = N
# endif

            CASE (  'SF1')
                  IDTSF1   = N
            CASE (  'SF2')
                  IDTSF2   = N
            CASE (  'SF3')
                  IDTSF3   = N
            CASE (  'SF4')
                  IDTSF4   = N
            CASE (  'SF5')
                  IDTSF5   = N
            CASE (  'SF6')
                  IDTSF6   = N
            CASE (  'SF7')
                  IDTSF7   = N
            CASE (  'SF8')
                  IDTSF8   = N
            CASE (  'SF9')
                  IDTSF9   = N
            CASE (  'SF10')
                  IDTSF10   = N
            CASE (  'SF11')
                  IDTSF11   = N
            CASE (  'SF12')
                  IDTSF12   = N
# if  !defined( TOMAS12 )
            CASE (  'SF13')
                  IDTSF13   = N
            CASE (  'SF14')
                  IDTSF14   = N
            CASE (  'SF15')
                  IDTSF15   = N
#  if !defined( TOMAS15 )
            CASE (  'SF16')
                  IDTSF16   = N
            CASE (  'SF17')
                  IDTSF17   = N
            CASE (  'SF18')
                  IDTSF18   = N
            CASE (  'SF19')
                  IDTSF19   = N
            CASE (  'SF20')
                  IDTSF20   = N
            CASE (  'SF21')
                  IDTSF21   = N
            CASE (  'SF22')
                  IDTSF22   = N
            CASE (  'SF23')
                  IDTSF23   = N
            CASE (  'SF24')
                  IDTSF24   = N
            CASE (  'SF25')
                  IDTSF25   = N
            CASE (  'SF26')
                  IDTSF26   = N
            CASE (  'SF27')
                  IDTSF27   = N
            CASE (  'SF28')
                  IDTSF28   = N
            CASE (  'SF29')
                  IDTSF29   = N
            CASE (  'SF30')
                  IDTSF30   = N
#  endif
# endif
# if   defined( TOMAS40 )
            CASE (  'SF31')
                  IDTSF31   = N
            CASE (  'SF32')
                  IDTSF32   = N
            CASE (  'SF33')
                  IDTSF33   = N
            CASE (  'SF34')
                  IDTSF34   = N
            CASE (  'SF35')
                  IDTSF35   = N
            CASE (  'SF36')
                  IDTSF36   = N
            CASE (  'SF37')
                  IDTSF37   = N
            CASE (  'SF38')
                  IDTSF38   = N
            CASE (  'SF39')
                  IDTSF39   = N
            CASE (  'SF40')
                  IDTSF40   = N
# endif

            CASE (  'SS1')
                  IDTSS1   = N
            CASE (  'SS2')
                  IDTSS2   = N
            CASE (  'SS3')
                  IDTSS3   = N
            CASE (  'SS4')
                  IDTSS4   = N
            CASE (  'SS5')
                  IDTSS5   = N
            CASE (  'SS6')
                  IDTSS6   = N
            CASE (  'SS7')
                  IDTSS7   = N
            CASE (  'SS8')
                  IDTSS8   = N
            CASE (  'SS9')
                  IDTSS9   = N
            CASE (  'SS10')
                  IDTSS10   = N
            CASE (  'SS11')
                  IDTSS11   = N
            CASE (  'SS12')
                  IDTSS12   = N
# if  !defined( TOMAS12 )
            CASE (  'SS13')
                  IDTSS13   = N
            CASE (  'SS14')
                  IDTSS14   = N
            CASE (  'SS15')
                  IDTSS15   = N
#  if !defined( TOMAS15 )
            CASE (  'SS16')
                  IDTSS16   = N
            CASE (  'SS17')
                  IDTSS17   = N
            CASE (  'SS18')
                  IDTSS18   = N
            CASE (  'SS19')
                  IDTSS19   = N
            CASE (  'SS20')
                  IDTSS20   = N
            CASE (  'SS21')
                  IDTSS21   = N
            CASE (  'SS22')
                  IDTSS22   = N
            CASE (  'SS23')
                  IDTSS23   = N
            CASE (  'SS24')
                  IDTSS24   = N
            CASE (  'SS25')
                  IDTSS25   = N
            CASE (  'SS26')
                  IDTSS26   = N
            CASE (  'SS27')
                  IDTSS27   = N
            CASE (  'SS28')
                  IDTSS28   = N
            CASE (  'SS29')
                  IDTSS29   = N
            CASE (  'SS30')
                  IDTSS30   = N
#  endif
# endif
# if   defined( TOMAS40 )
            CASE (  'SS31')
                  IDTSS31   = N
            CASE (  'SS32')
                  IDTSS32   = N
            CASE (  'SS33')
                  IDTSS33   = N
            CASE (  'SS34')
                  IDTSS34   = N
            CASE (  'SS35')
                  IDTSS35   = N
            CASE (  'SS36')
                  IDTSS36   = N
            CASE (  'SS37')
                  IDTSS37   = N
            CASE (  'SS38')
                  IDTSS38   = N
            CASE (  'SS39')
                  IDTSS39   = N
            CASE (  'SS40')
                  IDTSS40   = N
# endif

            CASE ('ECOB1')
                IDTECOB1   = N
            CASE ('ECOB2')
                IDTECOB2   = N
            CASE ('ECOB3')
                IDTECOB3   = N
            CASE ('ECOB4')
                IDTECOB4   = N
            CASE ('ECOB5')
                IDTECOB5   = N
            CASE ('ECOB6')
                IDTECOB6   = N
            CASE ('ECOB7')
                IDTECOB7   = N
            CASE ('ECOB8')
                IDTECOB8   = N
            CASE ('ECOB9')
                IDTECOB9   = N
            CASE ('ECOB10')
                IDTECOB10   = N
            CASE ('ECOB11')
                IDTECOB11   = N
            CASE ('ECOB12')
                IDTECOB12   = N
# if  !defined( TOMAS12 )
            CASE ('ECOB13')
                IDTECOB13   = N
            CASE ('ECOB14')
                IDTECOB14   = N
            CASE ('ECOB15')
                IDTECOB15   = N
#  if !defined( TOMAS15 )
            CASE ('ECOB16')
                IDTECOB16   = N
            CASE ('ECOB17')
                IDTECOB17   = N
            CASE ('ECOB18')
                IDTECOB18   = N
            CASE ('ECOB19')
                IDTECOB19   = N
            CASE ('ECOB20')
                IDTECOB20   = N
            CASE ('ECOB21')
                IDTECOB21   = N
            CASE ('ECOB22')
                IDTECOB22   = N
            CASE ('ECOB23')
                IDTECOB23   = N
            CASE ('ECOB24')
                IDTECOB24   = N
            CASE ('ECOB25')
                IDTECOB25   = N
            CASE ('ECOB26')
                IDTECOB26   = N
            CASE ('ECOB27')
                IDTECOB27   = N
            CASE ('ECOB28')
                IDTECOB28   = N
            CASE ('ECOB29')
                IDTECOB29   = N
            CASE ('ECOB30')
                IDTECOB30   = N
#  endif
# endif
# if   defined( TOMAS40 )
            CASE ('ECOB31')
                IDTECOB31   = N
            CASE ('ECOB32')
                IDTECOB32   = N
            CASE ('ECOB33')
                IDTECOB33   = N
            CASE ('ECOB34')
                IDTECOB34   = N
            CASE ('ECOB35')
                IDTECOB35   = N
            CASE ('ECOB36')
                IDTECOB36   = N
            CASE ('ECOB37')
                IDTECOB37   = N
            CASE ('ECOB38')
                IDTECOB38   = N
            CASE ('ECOB39')
                IDTECOB39   = N
            CASE ('ECOB40')
                IDTECOB40   = N
# endif

            CASE ('ECIL1')
                IDTECIL1   = N
            CASE ('ECIL2')
                IDTECIL2   = N
            CASE ('ECIL3')
                IDTECIL3   = N
            CASE ('ECIL4')
                IDTECIL4   = N
            CASE ('ECIL5')
                IDTECIL5   = N
            CASE ('ECIL6')
                IDTECIL6   = N
            CASE ('ECIL7')
                IDTECIL7   = N
            CASE ('ECIL8')
                IDTECIL8   = N
            CASE ('ECIL9')
                IDTECIL9   = N
            CASE ('ECIL10')
                IDTECIL10   = N
            CASE ('ECIL11')
                IDTECIL11   = N
            CASE ('ECIL12')
                IDTECIL12   = N
# if  !defined( TOMAS12 )
            CASE ('ECIL13')
                IDTECIL13   = N
            CASE ('ECIL14')
                IDTECIL14   = N
            CASE ('ECIL15')
                IDTECIL15   = N
#  if !defined( TOMAS15 )
            CASE ('ECIL16')
                IDTECIL16   = N
            CASE ('ECIL17')
                IDTECIL17   = N
            CASE ('ECIL18')
                IDTECIL18   = N
            CASE ('ECIL19')
                IDTECIL19   = N
            CASE ('ECIL20')
                IDTECIL20   = N
            CASE ('ECIL21')
                IDTECIL21   = N
            CASE ('ECIL22')
                IDTECIL22   = N
            CASE ('ECIL23')
                IDTECIL23   = N
            CASE ('ECIL24')
                IDTECIL24   = N
            CASE ('ECIL25')
                IDTECIL25   = N
            CASE ('ECIL26')
                IDTECIL26   = N
            CASE ('ECIL27')
                IDTECIL27   = N
            CASE ('ECIL28')
                IDTECIL28   = N
            CASE ('ECIL29')
                IDTECIL29   = N
            CASE ('ECIL30')
                IDTECIL30   = N
#  endif
# endif
# if   defined( TOMAS40 )
            CASE ('ECIL31')
                IDTECIL31   = N
            CASE ('ECIL32')
                IDTECIL32   = N
            CASE ('ECIL33')
                IDTECIL33   = N
            CASE ('ECIL34')
                IDTECIL34   = N
            CASE ('ECIL35')
                IDTECIL35   = N
            CASE ('ECIL36')
                IDTECIL36   = N
            CASE ('ECIL37')
                IDTECIL37   = N
            CASE ('ECIL38')
                IDTECIL38   = N
            CASE ('ECIL39')
                IDTECIL39   = N
            CASE ('ECIL40')
                IDTECIL40   = N
# endif

            CASE ('OCOB1')
                IDTOCOB1   = N
            CASE ('OCOB2')
                IDTOCOB2   = N
            CASE ('OCOB3')
                IDTOCOB3   = N
            CASE ('OCOB4')
                IDTOCOB4   = N
            CASE ('OCOB5')
                IDTOCOB5   = N
            CASE ('OCOB6')
                IDTOCOB6   = N
            CASE ('OCOB7')
                IDTOCOB7   = N
            CASE ('OCOB8')
                IDTOCOB8   = N
            CASE ('OCOB9')
                IDTOCOB9   = N
            CASE ('OCOB10')
                IDTOCOB10   = N
            CASE ('OCOB11')
                IDTOCOB11   = N
            CASE ('OCOB12')
                IDTOCOB12   = N
# if  !defined( TOMAS12 )
            CASE ('OCOB13')
                IDTOCOB13   = N
            CASE ('OCOB14')
                IDTOCOB14   = N
            CASE ('OCOB15')
                IDTOCOB15   = N
#  if !defined( TOMAS15 )
            CASE ('OCOB16')
                IDTOCOB16   = N
            CASE ('OCOB17')
                IDTOCOB17   = N
            CASE ('OCOB18')
                IDTOCOB18   = N
            CASE ('OCOB19')
                IDTOCOB19   = N
            CASE ('OCOB20')
                IDTOCOB20   = N
            CASE ('OCOB21')
                IDTOCOB21   = N
            CASE ('OCOB22')
                IDTOCOB22   = N
            CASE ('OCOB23')
                IDTOCOB23   = N
            CASE ('OCOB24')
                IDTOCOB24   = N
            CASE ('OCOB25')
                IDTOCOB25   = N
            CASE ('OCOB26')
                IDTOCOB26   = N
            CASE ('OCOB27')
                IDTOCOB27   = N
            CASE ('OCOB28')
                IDTOCOB28   = N
            CASE ('OCOB29')
                IDTOCOB29   = N
            CASE ('OCOB30')
                IDTOCOB30   = N
#  endif
# endif
# if   defined( TOMAS40 )
            CASE ('OCOB31')
                IDTOCOB31   = N
            CASE ('OCOB32')
                IDTOCOB32   = N
            CASE ('OCOB33')
                IDTOCOB33   = N
            CASE ('OCOB34')
                IDTOCOB34   = N
            CASE ('OCOB35')
                IDTOCOB35   = N
            CASE ('OCOB36')
                IDTOCOB36   = N
            CASE ('OCOB37')
                IDTOCOB37   = N
            CASE ('OCOB38')
                IDTOCOB38   = N
            CASE ('OCOB39')
                IDTOCOB39   = N
            CASE ('OCOB40')
                IDTOCOB40   = N
# endif

            CASE ('OCIL1')
                IDTOCIL1   = N
            CASE ('OCIL2')
                IDTOCIL2   = N
            CASE ('OCIL3')
                IDTOCIL3   = N
            CASE ('OCIL4')
                IDTOCIL4   = N
            CASE ('OCIL5')
                IDTOCIL5   = N
            CASE ('OCIL6')
                IDTOCIL6   = N
            CASE ('OCIL7')
                IDTOCIL7   = N
            CASE ('OCIL8')
                IDTOCIL8   = N
            CASE ('OCIL9')
                IDTOCIL9   = N
            CASE ('OCIL10')
                IDTOCIL10   = N
            CASE ('OCIL11')
                IDTOCIL11   = N
            CASE ('OCIL12')
                IDTOCIL12   = N
# if  !defined( TOMAS12 )
            CASE ('OCIL13')
                IDTOCIL13   = N
            CASE ('OCIL14')
                IDTOCIL14   = N
            CASE ('OCIL15')
                IDTOCIL15   = N
#  if !defined( TOMAS15 )
            CASE ('OCIL16')
                IDTOCIL16   = N
            CASE ('OCIL17')
                IDTOCIL17   = N
            CASE ('OCIL18')
                IDTOCIL18   = N
            CASE ('OCIL19')
                IDTOCIL19   = N
            CASE ('OCIL20')
                IDTOCIL20   = N
            CASE ('OCIL21')
                IDTOCIL21   = N
            CASE ('OCIL22')
                IDTOCIL22   = N
            CASE ('OCIL23')
                IDTOCIL23   = N
            CASE ('OCIL24')
                IDTOCIL24   = N
            CASE ('OCIL25')
                IDTOCIL25   = N
            CASE ('OCIL26')
                IDTOCIL26   = N
            CASE ('OCIL27')
                IDTOCIL27   = N
            CASE ('OCIL28')
                IDTOCIL28   = N
            CASE ('OCIL29')
                IDTOCIL29   = N
            CASE ('OCIL30')
                IDTOCIL30   = N
#  endif
# endif
# if   defined( TOMAS40 )
            CASE ('OCIL31')
                IDTOCIL31   = N
            CASE ('OCIL32')
                IDTOCIL32   = N
            CASE ('OCIL33')
                IDTOCIL33   = N
            CASE ('OCIL34')
                IDTOCIL34   = N
            CASE ('OCIL35')
                IDTOCIL35   = N
            CASE ('OCIL36')
                IDTOCIL36   = N
            CASE ('OCIL37')
                IDTOCIL37   = N
            CASE ('OCIL38')
                IDTOCIL38   = N
            CASE ('OCIL39')
                IDTOCIL39   = N
            CASE ('OCIL40')
                IDTOCIL40   = N
# endif

            CASE ('DUST1')
                IDTDUST1   = N
            CASE ('DUST2')
                IDTDUST2   = N
            CASE ('DUST3')
                IDTDUST3   = N
            CASE ('DUST4')
                IDTDUST4   = N
            CASE ('DUST5')
                IDTDUST5   = N
            CASE ('DUST6')
                IDTDUST6   = N
            CASE ('DUST7')
                IDTDUST7   = N
            CASE ('DUST8')
                IDTDUST8   = N
            CASE ('DUST9')
                IDTDUST9   = N
            CASE ('DUST10')
                IDTDUST10   = N
            CASE ('DUST11')
                IDTDUST11   = N
            CASE ('DUST12')
                IDTDUST12   = N
# if  !defined( TOMAS12 )
            CASE ('DUST13')
                IDTDUST13   = N
            CASE ('DUST14')
                IDTDUST14   = N
            CASE ('DUST15')
                IDTDUST15   = N
#  if !defined( TOMAS15 )
            CASE ('DUST16')
                IDTDUST16   = N
            CASE ('DUST17')
                IDTDUST17   = N
            CASE ('DUST18')
                IDTDUST18   = N
            CASE ('DUST19')
                IDTDUST19   = N
            CASE ('DUST20')
                IDTDUST20   = N
            CASE ('DUST21')
                IDTDUST21   = N
            CASE ('DUST22')
                IDTDUST22   = N
            CASE ('DUST23')
                IDTDUST23   = N
            CASE ('DUST24')
                IDTDUST24   = N
            CASE ('DUST25')
                IDTDUST25   = N
            CASE ('DUST26')
                IDTDUST26   = N
            CASE ('DUST27')
                IDTDUST27   = N
            CASE ('DUST28')
                IDTDUST28   = N
            CASE ('DUST29')
                IDTDUST29   = N
            CASE ('DUST30')
                IDTDUST30   = N
#  endif
# endif
# if   defined( TOMAS40 )
            CASE ('DUST31')
                IDTDUST31   = N
            CASE ('DUST32')
                IDTDUST32   = N
            CASE ('DUST33')
                IDTDUST33   = N
            CASE ('DUST34')
                IDTDUST34   = N
            CASE ('DUST35')
                IDTDUST35   = N
            CASE ('DUST36')
                IDTDUST36   = N
            CASE ('DUST37')
                IDTDUST37   = N
            CASE ('DUST38')
                IDTDUST38   = N
            CASE ('DUST39')
                IDTDUST39   = N
            CASE ('DUST40')
                IDTDUST40   = N
# endif

            CASE (  'AW1')
                  IDTAW1   = N
            CASE (  'AW2')
                  IDTAW2   = N
            CASE (  'AW3')
                  IDTAW3   = N
            CASE (  'AW4')
                  IDTAW4   = N
            CASE (  'AW5')
                  IDTAW5   = N
            CASE (  'AW6')
                  IDTAW6   = N
            CASE (  'AW7')
                  IDTAW7   = N
            CASE (  'AW8')
                  IDTAW8   = N
            CASE (  'AW9')
                  IDTAW9   = N
            CASE (  'AW10')
                  IDTAW10   = N
            CASE (  'AW11')
                  IDTAW11   = N
            CASE (  'AW12')
                  IDTAW12   = N
# if  !defined( TOMAS12 )
            CASE (  'AW13')
                  IDTAW13   = N
            CASE (  'AW14')
                  IDTAW14   = N
            CASE (  'AW15')
                  IDTAW15   = N
#  if !defined( TOMAS15 )
            CASE (  'AW16')
                  IDTAW16   = N
            CASE (  'AW17')
                  IDTAW17   = N
            CASE (  'AW18')
                  IDTAW18   = N
            CASE (  'AW19')
                  IDTAW19   = N
            CASE (  'AW20')
                  IDTAW20   = N
            CASE (  'AW21')
                  IDTAW21   = N
            CASE (  'AW22')
                  IDTAW22   = N
            CASE (  'AW23')
                  IDTAW23   = N
            CASE (  'AW24')
                  IDTAW24   = N
            CASE (  'AW25')
                  IDTAW25   = N
            CASE (  'AW26')
                  IDTAW26   = N
            CASE (  'AW27')
                  IDTAW27   = N
            CASE (  'AW28')
                  IDTAW28   = N
            CASE (  'AW29')
                  IDTAW29   = N
            CASE (  'AW30')
                  IDTAW30   = N
#  endif
# endif
# if   defined( TOMAS40 )
            CASE (  'AW31')
                  IDTAW31   = N
            CASE (  'AW32')
                  IDTAW32   = N
            CASE (  'AW33')
                  IDTAW33   = N
            CASE (  'AW34')
                  IDTAW34   = N
            CASE (  'AW35')
                  IDTAW35   = N
            CASE (  'AW36')
                  IDTAW36   = N
            CASE (  'AW37')
                  IDTAW37   = N
            CASE (  'AW38')
                  IDTAW38   = N
            CASE (  'AW39')
                  IDTAW39   = N
            CASE (  'AW40')
                  IDTAW40   = N
# endif

#endif

            CASE DEFAULT
               ! Nothing

         END SELECT
      ENDDO
     
      !=================================================================
      ! Initialize the IDExxx flags, which denote tracers that
      ! are emitted.  Remove the old previous hardwiring
      ! (hotp, bmy, 8/17/11)
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         ! Initialize counters
         ! Number of anthro emissions handled in anthroems and SMVGEAR:
         NEMANTHRO = 0

         ! Number of additional emissions handled in SMVGEAR:
         NEMBIOG   = 0

         ! The first 12 tracers are in the original order and thus
         ! IDEs are the same as listed below
         IF ( IDTNO > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDENO     = NEMANTHRO
         ENDIF
         IF ( IDTCO > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDECO     = NEMANTHRO
         ENDIF
         IF ( IDTPRPE > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDEPRPE   = NEMANTHRO
         ENDIF
         IF ( IDTC3H8 > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDEC3H8   = NEMANTHRO
         ENDIF
         IF ( IDTALK4 > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDEALK4   = NEMANTHRO
         ENDIF
         IF ( IDTC2H6 > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDEC2H6   = NEMANTHRO
         ENDIF
         IF ( IDTACET > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDEACET   = NEMANTHRO
         ENDIF
         IF ( IDTMEK  > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDEMEK    = NEMANTHRO  ! NOX typo fix (hotp 8/5/09)
         ENDIF
         IF ( IDTALD2 > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDEALD2   = NEMANTHRO  ! NOX typo fix (hotp 8/5/09)
         ENDIF
         IF ( IDTCH2O > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDECH2O   = NEMANTHRO
         ENDIF
         IF ( IDTO3 > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDEO3     = NEMANTHRO
         ENDIF
         IF ( IDTHNO3 > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDEHNO3   = NEMANTHRO
         ENDIF
         ! Aromatics for SOA5 simulation
         IF ( IDTBENZ > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDEBENZ   = NEMANTHRO
         ENDIF
         IF ( IDTTOLU > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDETOLU   = NEMANTHRO
         ENDIF
         IF ( IDTXYLE > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDEXYLE   = NEMANTHRO
         ENDIF
         ! SOAupdate: Add NAP emissions (hotp 8/24/09)
         ! Used for SOA simulation
         IF ( IDTNAP  > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDENAP    = NEMANTHRO
         ENDIF
         IF ( IDTC2H4 > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDEC2H4   = NEMANTHRO
         ENDIF
         IF ( IDTC2H2 > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDEC2H2   = NEMANTHRO
         ENDIF
         
         !%%%%%%bmy add these here%%%%%%%%%%

         ! For now, comment this out, as RETRO does not have RCHO emissions
         ! (bmy, 8/17/11)
         !IF ( IDTRCHO > 0 ) THEN
         !   NEMANTHRO = NEMANTHRO + 1
         !   IDERCHO  = NEMANTHRO
         !ENDIF

         IF ( IDTHCOOH > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDEHCOOH  = NEMANTHRO
         ENDIF

         IF ( IDTACTA > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDEACTA   = NEMANTHRO
          ENDIF
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         ! Now add additional species treated in SMVGEAR that are NOT
         ! anthropogenic, increment NEMBIOG
         ! All non-anthropogenic species should be last (see anthroems.f
         ! and loop over NEMANTHRO)
         IF ( IDTISOP > 0 ) THEN
           NEMBIOG    = NEMBIOG   + 1
           IDEISOP    = NEMANTHRO + NEMBIOG
         ENDIF
         IF ( IDTGLYX > 0 ) THEN
            NEMBIOG   = NEMBIOG + 1
            IDEGLYX   = NEMANTHRO + NEMBIOG
         ENDIF
         IF ( IDTMGLY > 0 ) THEN
            NEMBIOG   = NEMBIOG + 1
            IDEMGLY   = NEMANTHRO + NEMBIOG
         ENDIF
         IF ( IDTGLYC > 0 ) THEN
            NEMBIOG   = NEMBIOG + 1
            IDEGLYC   = NEMANTHRO + NEMBIOG
         ENDIF
         IF ( IDTHAC > 0 ) THEN
            NEMBIOG   = NEMBIOG + 1
            IDEHAC    = NEMANTHRO + NEMBIOG
         ENDIF
         IF ( IDTMONX > 0 ) THEN
            NEMBIOG   = NEMBIOG + 1
            IDEMONX   = NEMANTHRO + NEMBIOG
         ENDIF
         IF ( IDTMBO > 0 ) THEN
            NEMBIOG   = NEMBIOG + 1
            IDEMBO    = NEMANTHRO + NEMBIOG
         ENDIF
         ! Add bromine species biogenic emissions (jpp, mpayer, 12/28/11)
         IF ( IDTCHBr3 > 0 ) THEN
            NEMBIOG   = NEMBIOG + 1
            IDECHBr3  = NEMANTHRO + NEMBIOG
         ENDIF
         IF ( IDTCH2Br2 > 0 ) THEN
            NEMBIOG   = NEMBIOG + 1
            IDECH2Br2 = NEMANTHRO + NEMBIOG
         ENDIF
         IF ( IDTBr2 > 0 ) THEN
            NEMBIOG   = NEMBIOG + 1
            IDEBr2    = NEMANTHRO + NEMBIOG
         ENDIF
      ENDIF
      
      !=================================================================
      ! Fill IDEMS with appropriate tracer ID #'s
      !
      ! NOTE: IDEMS is in "comode.h", maybe later split this off into
      ! an F90 module somehow.  Think about this later. (bmy, 11/12/02)
      !=================================================================
      IF ( IDENO   /= 0 ) IDEMS(IDENO  ) = IDTNO
      IF ( IDECO   /= 0 ) IDEMS(IDECO  ) = IDTCO
      IF ( IDEPRPE /= 0 ) IDEMS(IDEPRPE) = IDTPRPE
      IF ( IDEC3H8 /= 0 ) IDEMS(IDEC3H8) = IDTC3H8
      IF ( IDEALK4 /= 0 ) IDEMS(IDEALK4) = IDTALK4
      IF ( IDEC2H6 /= 0 ) IDEMS(IDEC2H6) = IDTC2H6
      IF ( IDEISOP /= 0 ) IDEMS(IDEISOP) = IDTISOP
      IF ( IDEACET /= 0 ) IDEMS(IDEACET) = IDTACET
      IF ( IDEMEK  /= 0 ) IDEMS(IDEMEK ) = IDTMEK
      IF ( IDEALD2 /= 0 ) IDEMS(IDEALD2) = IDTALD2
      IF ( IDECH2O /= 0 ) IDEMS(IDECH2O) = IDTCH2O
      IF ( IDEO3   /= 0 ) IDEMS(IDEO3  ) = IDTO3
      IF ( IDEHNO3 /= 0 ) IDEMS(IDEHNO3) = IDTHNO3  ! PHS
      IF ( IDEGLYX /= 0 ) IDEMS(IDEGLYX) = IDTGLYX
      IF ( IDEMGLY /= 0 ) IDEMS(IDEMGLY) = IDTMGLY
      IF ( IDEBENZ /= 0 ) IDEMS(IDEBENZ) = IDTBENZ
      IF ( IDETOLU /= 0 ) IDEMS(IDETOLU) = IDTTOLU
      IF ( IDEXYLE /= 0 ) IDEMS(IDEXYLE) = IDTXYLE
      IF ( IDEMONX /= 0 ) IDEMS(IDEMONX) = IDTMONX
      IF ( IDEC2H4 /= 0 ) IDEMS(IDEC2H4) = IDTC2H4
      IF ( IDEC2H2 /= 0 ) IDEMS(IDEC2H2) = IDTC2H2
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF ( IDERCHO /= 0 )  IDEMS(IDERCHO) = IDTRCHO
      IF ( IDEHCOOH /= 0 ) IDEMS(IDEHCOOH)= IDTHCOOH
      IF ( IDEACTA /= 0 )  IDEMS(IDEACTA) = IDTACTA
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF ( IDEMBO  /= 0 ) IDEMS(IDEMBO ) = IDTMBO
      IF ( IDEGLYC /= 0 ) IDEMS(IDEGLYC) = IDTGLYC
      IF ( IDEHAC  /= 0 ) IDEMS(IDEHAC ) = IDTHAC
      ! jpp:
      IF ( IDECHBr3  /= 0 ) IDEMS(IDECHBr3 ) = IDTCHBr3
      IF ( IDECH2Br2 /= 0 ) IDEMS(IDECH2Br2) = IDTCH2Br2
      IF ( IDEBr2    /= 0 ) IDEMS(IDEBr2   ) = IDTBr2
      ! SOAupdate: for gas phase NAP chemistry (hotp 6/24/09)
      ! semivolpoa5: update for IVOC chemistry (hotp 3/30/09)
      IF ( IDENAP    /= 0 ) IDEMS(IDENAP   ) = IDTNAP

      ! Echo anthro & biogenic emitted tracers
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) IDEMS ( 1:NEMANTHRO+NEMBIOG )
 100     FORMAT( /, 'TRACERID: Emitted tracers (anthro & bio) :', 27i3 )

         ! print additional information (hotp 8/4/09)
         WRITE( 6, 115 ) 'TRACERID: CHECK FULLCHEM IDEs'
         WRITE( 6, 115 ) '---------------------------------------------'
         WRITE( 6, 115 ) 'NEMANTHRO: 10 for standard sims b4 v811'
         WRITE( 6, 115 ) '          + 2 for HNO3 and Ox'
         WRITE( 6, 115 ) '          + 5 for aromatics (B,T,X,C2H2,C2H4)'
         WRITE( 6, 115 ) '          + 1 for NAP'  ! hotp
         WRITE( 6, 115 ) 'NEMANTHRO: 18 max'
         WRITE( 6, 120 ) 'THIS NEMANTHRO: ', NEMANTHRO
         WRITE( 6, 115 ) '---------------------------------------------'
         WRITE( 6, 115 ) 'NEMBIOG:    1 for isoprene'
         WRITE( 6, 115 ) '          + 6 for dicarbonyls'
         WRITE( 6, 115 ) '          + 3 for bromine'
         WRITE( 6, 115 ) 'NEMBIOG:   10 max'
         WRITE( 6, 120 ) 'THIS NEMBIOG:   ', NEMBIOG
         WRITE( 6, 115 ) '-----------------------------------------'
         WRITE( 6, 120 ) 'IDENO   ', IDENO
         WRITE( 6, 120 ) 'IDECO   ', IDECO
         WRITE( 6, 120 ) 'IDEPRPE ', IDEPRPE
         WRITE( 6, 120 ) 'IDEC3H8 ', IDEC3H8
         WRITE( 6, 120 ) 'IDEALK4 ', IDEALK4
         WRITE( 6, 120 ) 'IDEC2H6 ', IDEC2H6
         WRITE( 6, 120 ) 'IDEACET ', IDEACET
         WRITE( 6, 120 ) 'IDEMEK  ', IDEMEK
         WRITE( 6, 120 ) 'IDEALD2 ', IDEALD2
         WRITE( 6, 120 ) 'IDECH2O ', IDECH2O
         WRITE( 6, 120 ) 'IDEO3   ', IDEO3
         WRITE( 6, 120 ) 'IDEHNO3 ', IDEHNO3
         WRITE( 6, 120 ) 'IDEISOP ', IDEISOP
         WRITE( 6, 120 ) 'IDEGLYX ', IDEGLYX
         WRITE( 6, 120 ) 'IDEMGLY ', IDEMGLY
         WRITE( 6, 120 ) 'IDEBENZ ', IDEBENZ
         WRITE( 6, 120 ) 'IDETOLU ', IDETOLU
         WRITE( 6, 120 ) 'IDEXYLE ', IDEXYLE
         WRITE( 6, 120 ) 'IDEMONX ', IDEMONX
         WRITE( 6, 120 ) 'IDEC2H4 ', IDEC2H4
         WRITE( 6, 120 ) 'IDEC2H2 ', IDEC2H2
         WRITE( 6, 120 ) 'IDERCHO ', IDERCHO
         WRITE( 6, 120 ) 'IDEHCOOH', IDEHCOOH
         WRITE( 6, 120 ) 'IDEACTA ', IDEACTA
         WRITE( 6, 120 ) 'IDEMBO  ', IDEMBO
         WRITE( 6, 120 ) 'IDEGLYC ', IDEGLYC
         WRITE( 6, 120 ) 'IDEHAC  ', IDEHAC
         WRITE( 6, 120 ) 'IDECHBr3 ', IDECHBr3 
         WRITE( 6, 120 ) 'IDECH2Br2', IDECH2Br2 
         WRITE( 6, 120 ) 'IDEBr2   ', IDEBr2 
         ! Add NAP info (hotp 8/24/09)
         WRITE( 6, 120 ) 'IDENAP   ', IDENAP
         WRITE( 6, 115 ) '-----------------------------------------'
      ENDIF
 115  FORMAT( a )
 120  FORMAT( a, i3 )

      END SUBROUTINE TRACERID

!------------------------------------------------------------------------------

      SUBROUTINE SETTRACE( am_I_Root, Input_Opt, State_Chm, RC )
!
!******************************************************************************
!  Subroutine SETTRACE flags certain chemical species w/in the SMVGEAR full
!  chemistry mechanism. (lwh, jyl, gmg, djj, 1990's; bmy, 11/12/02, 10/3/05)
!
!  Arguments as Input: 
!  ============================================================================
!  (1 ) NTRACER : Number of GEOS-CHEM tracers to process
!
!  NOTES:
!  (1 ) Added comment header.
!  (2 ) Now initialize IDDMS, IDSO2, IDSO4, IDMSA.  Updated comments,
!        cosmetic changes. (rjp, bmy, 3/23/03)
!  (3 ) Currently there are only families for the troposphere, so manually 
!        set NCS = NCSURBAN.  Replace NAMESPEC w/ NAMEGAS for SMVGEAR II. 
!        (bdf, bmy, 4/23/03)
!  (4 ) Make sure IDEMIS etc doesn't go out of array bounds (bmy, 4/26/04)
!  (5 ) Removed NTRACER from the arg list, we can use N_TRACERS from 
!        "tracer_mod.f".  Now references "tracer_mod.f".  Now does not have 
!        to read the "tracer.dat" file. (bmy, 7/20/04)
!  (6 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (7 ) Added matching for IDH2SO4 for sulfuric acid (win, 6/23/09)
!  (8 ) Comment out IDBENZ, IDTOLU, IDXYLE, IDMONX, IDGLYX, IDMGLY 
!        (ccc, 01/29/10)
!  (8 ) Add IDs for Henze's aromatics (hotp, 01/10)
!  (9 ) Add debug print. (hotp, 01/10)
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  16 Oct 2012 - R. Yantosca - Renamed chem state object to State_Chm
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  13 Nov 2012 - R. Yantosca - Now pass Input Options object as an argument
!  05 Mar 2013 - R. Yantosca - Now TRIM the TRACER_NAME in the ND70 debug print
!  06 Mar 2013 - R. Yantosca - Now move call to REGISTER_SPECIES to routine
!                              GIGC_INIT_SIMULATION 
!  13 Aug 2013 - M. Sulprizio- Add modifications for IVOC chemistry (H. Pye)
!******************************************************************************
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE COMODE_LOOP_MOD
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput          
      USE GIGC_State_Chm_Mod, ONLY : ChmState          
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!     
! !INPUT/OUTPUT PARAMETERS
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State obj
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL            :: LPRT
      INTEGER            :: N_TRACERS
      INTEGER            :: I, J, T, C, STAT

      ! Arrays
      INTEGER            :: ID_EMITTED    ( Input_Opt%MAX_TRCS )
      REAL*8             :: TRACER_COEFF  ( Input_Opt%MAX_TRCS, 
     &                                      Input_Opt%MAX_MEMB )
      CHARACTER(LEN=255) :: TRACER_CONST  ( Input_Opt%MAX_TRCS, 
     &                                      Input_Opt%MAX_MEMB )
      INTEGER            :: TRACER_N_CONST( Input_Opt%MAX_TRCS )
      CHARACTER(LEN=255) :: TRACER_NAME   ( Input_Opt%MAX_TRCS )
 
      !=================================================================
      ! SETTRACE begins here!
      !=================================================================

      ! Assume success
      RC             = GIGC_SUCCESS

      ! Copy values from Input_Opt to local variables
      LPRT           = Input_Opt%LPRT
      N_TRACERS      = Input_Opt%N_TRACERS
      ID_EMITTED     = Input_Opt%ID_EMITTED
      TRACER_COEFF   = Input_Opt%TRACER_COEFF
      TRACER_CONST   = Input_Opt%TRACER_CONST 
      TRACER_N_CONST = Input_Opt%TRACER_N_CONST
      TRACER_NAME    = Input_Opt%TRACER_NAME

      ! Reset NCS = NCSURBAN, since we have defined our GEOS-CHEM
      ! mechanism in the urban slot of SMVGEAR II (bmy, 4/21/03)
      NCS = NCSURBAN

      DO I = 1, NSPEC(NCS)
         IF ( NAMEGAS(I) == 'O3'     ) IDO3     = I
         IF ( NAMEGAS(I) == 'NO2'    ) IDNO2    = I
         IF ( NAMEGAS(I) == 'NO3'    ) IDNO3    = I
         IF ( NAMEGAS(I) == 'N2O5'   ) IDN2O5   = I
         IF ( NAMEGAS(I) == 'HNO4'   ) IDHNO4   = I
         IF ( NAMEGAS(I) == 'HNO2'   ) IDHNO2   = I
         IF ( NAMEGAS(I) == 'NO'     ) IDNO     = I
         IF ( NAMEGAS(I) == 'CO'     ) IDCO     = I
         IF ( NAMEGAS(I) == 'PRPE'   ) IDPRPE   = I
         IF ( NAMEGAS(I) == 'C3H8'   ) IDC3H8   = I
         IF ( NAMEGAS(I) == 'ISOP'   ) IDISOP   = I
         IF ( NAMEGAS(I) == 'ALK4'   ) IDALK4   = I
         IF ( NAMEGAS(I) == 'PAN'    ) IDPAN    = I
         IF ( NAMEGAS(I) == 'GLPAN'  ) IDGLPAN  = I
         IF ( NAMEGAS(I) == 'GPAN'   ) IDGPAN   = I
         IF ( NAMEGAS(I) == 'PMN'    ) IDPMN    = I
         IF ( NAMEGAS(I) == 'PPN'    ) IDPPN    = I
         IF ( NAMEGAS(I) == 'HNO3'   ) IDHNO3   = I
         IF ( NAMEGAS(I) == 'OH'     ) IDOH     = I
         IF ( NAMEGAS(I) == 'HO2'    ) IDHO2    = I !(rvm, bmy, 2/27/02)
         IF ( NAMEGAS(I) == 'H2O2'   ) IDH2O2   = I
         IF ( NAMEGAS(I) == 'ACET'   ) IDACET   = I
         IF ( NAMEGAS(I) == 'MEK'    ) IDMEK    = I
         IF ( NAMEGAS(I) == 'ALD2'   ) IDALD2   = I
         IF ( NAMEGAS(I) == 'RCHO'   ) IDRCHO   = I
         IF ( NAMEGAS(I) == 'MVK'    ) IDMVK    = I
         IF ( NAMEGAS(I) == 'MACR'   ) IDMACR   = I
         IF ( NAMEGAS(I) == 'ISN2'   ) IDISN2   = I
         IF ( NAMEGAS(I) == 'R4N2'   ) IDR4N2   = I
         IF ( NAMEGAS(I) == 'CH2O'   ) IDCH2O   = I
         IF ( NAMEGAS(I) == 'C2H6'   ) IDC2H6   = I
         IF ( NAMEGAS(I) == 'DMS'    ) IDDMS    = I !(rjp, bmy, 3/23/03)
         IF ( NAMEGAS(I) == 'SO2'    ) IDSO2    = I !(rjp, bmy, 3/23/03)
         IF ( NAMEGAS(I) == 'SO4'    ) IDSO4    = I !(rjp, bmy, 3/23/03)
         IF ( NAMEGAS(I) == 'MSA'    ) IDMSA    = I !(rjp, bmy, 3/23/03)
         IF ( NAMEGAS(I) == 'DRYNO2' ) IDDRYNO2 = I
         IF ( NAMEGAS(I) == 'DRYPAN' ) IDDRYPAN = I
         IF ( NAMEGAS(I) == 'DRYO3 ' ) IDDRYO3  = I
!-Not used. (ccc, 01/27/10)
!         IF ( NAMEGAS(I) == 'BENZ'   ) IDBENZ   = I 
!         IF ( NAMEGAS(I) == 'TOLU'   ) IDTOLU   = I 
!         IF ( NAMEGAS(I) == 'XYLE'   ) IDXYLE   = I 
!         IF ( NAMEGAS(I) == 'MONX'   ) IDMONX   = I 
!         IF ( NAMEGAS(I) == 'GLYX'   ) IDGLYX   = I
!         IF ( NAMEGAS(I) == 'MGLY'   ) IDMGLY   = I
         IF ( NAMEGAS(I) == 'DRYGLYX') IDDRYGLYX = I
         IF ( NAMEGAS(I) == 'DRYMGLY') IDDRYMGLY = I
         IF ( NAMEGAS(I) == 'C2H4'   ) IDC2H4    = I
         IF ( NAMEGAS(I) == 'C2H2'   ) IDC2H2    = I
         IF ( NAMEGAS(I) == 'MBO'    ) IDMBO     = I
         IF ( NAMEGAS(I) == 'GLYC'   ) IDGLYC    = I
         IF ( NAMEGAS(I) == 'HAC'    ) IDHAC     = I
         IF ( NAMEGAS(I) == 'APAN'   ) IDAPAN    = I
         IF ( NAMEGAS(I) == 'ENPAN'  ) IDENPAN   = I
         IF ( NAMEGAS(I) == 'MPAN'   ) IDMPAN    = I
         IF ( NAMEGAS(I) == 'NIPAN'  ) IDNIPAN   = I
         IF ( NAMEGAS(I) == 'MPN'    ) IDMPN     = I
         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !jpp: for bromine
         IF ( NAMEGAS(I) == 'Br2'     ) IDBr2     = I
         IF ( NAMEGAS(I) == 'Br'      ) IDBr      = I
         IF ( NAMEGAS(I) == 'BrO'     ) IDBrO     = I
         IF ( NAMEGAS(I) == 'HOBr'    ) IDHOBr    = I
         IF ( NAMEGAS(I) == 'HBr'     ) IDHBr     = I
         IF ( NAMEGAS(I) == 'BrNO3'   ) IDBrNO3   = I
         IF ( NAMEGAS(I) == 'BrNO2'   ) IDBrNO2   = I
         IF ( NAMEGAS(I) == 'CH3Br'   ) IDCH3Br   = I
         IF ( NAMEGAS(I) == 'CH2Br2'  ) IDCH2Br2  = I
         IF ( NAMEGAS(I) == 'CHBr3'   ) IDCHBr3   = I
         IF ( NAMEGAS(I) == 'DRYHOBr' ) IDDRYHOBr = I
         IF ( NAMEGAS(I) == 'DRYHBr'  ) IDDRYHBr  = I
         IF ( NAMEGAS(I) == 'DRYBrNO3') IDDRYBrNO3= I
         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++

         IF ( NAMEGAS(I) == 'DRYAPAN' )  IDDRYAPAN  = I
         IF ( NAMEGAS(I) == 'DRYENPAN')  IDDRYENPAN = I
         IF ( NAMEGAS(I) == 'DRYGLPAN')  IDDRYGLPAN = I
         IF ( NAMEGAS(I) == 'DRYGPAN' )  IDDRYGPAN  = I
         IF ( NAMEGAS(I) == 'DRYMPAN' )  IDDRYMPAN  = I
         IF ( NAMEGAS(I) == 'DRYNIPAN')  IDDRYNIPAN = I
         ! added for aromatics (dkh, 10/06/06)  
         IF ( NAMEGAS(I) == 'LBRO2H'   ) IDLBRO2H   = I
         IF ( NAMEGAS(I) == 'LBRO2N'   ) IDLBRO2N   = I
         IF ( NAMEGAS(I) == 'LTRO2H'   ) IDLTRO2H   = I
         IF ( NAMEGAS(I) == 'LTRO2N'   ) IDLTRO2N   = I
         IF ( NAMEGAS(I) == 'LXRO2H'   ) IDLXRO2H   = I
         IF ( NAMEGAS(I) == 'LXRO2N'   ) IDLXRO2N   = I         
         ! SOAupdate: for gas phase NAP chemistry (hotp 6/24/09)
         ! semivolpoa5: add IVOC chemistry (hotp 3/30/09)
         IF ( NAMEGAS(I) == 'NAP'      ) IDNAP      = I 
         IF ( NAMEGAS(I) == 'NRO2'     ) IDNRO2     = I
         IF ( NAMEGAS(I) == 'LNRO2H'   ) IDLNRO2H   = I
         IF ( NAMEGAS(I) == 'LNRO2N'   ) IDLNRO2N   = I
      ENDDO

      !=================================================================
      ! Initialize arrays
      !=================================================================
      DO I=1, NNNTRID
         NMEMBER(I)  = 0
         IDEMIS(I)   = 0   
         DO J=1, MMMEMBER
            IDTRMB(I, J)= 0
            CTRMB(I, J)= 0.
         ENDDO
      ENDDO

      !=================================================================
      ! Save IDs for tracers (sequence in NAMESPEC.)
      !
      ! IDTRMB(T,C)  = species number for J'th component of tracer I
      ! CTRMB(T,C)+1 = coefficient of tracer constituent (e.g., each NO3
      !                molec. represents 2 units of Ox, so CTRMB=1)
      ! TRACER_N(T)   = number of component species in tracer I
      ! IDEMIS(T)    = which component of tracer I (in IDTRMB sense)
      !                receives the emissions
      ! NIDEMIS      = 0,1 -- indicates which species is emitting species.
      !                If there is only one species in tracer family and
      !                it's emitted, you still need a "1" in the spot.
      ! ljm changes: now read input from data file, tracer.dat
      !=================================================================

      ! Loop over tracers
      DO T = 1, N_TRACERS
         
         ! Number of constituents that tracer T has
         NMEMBER(T) = TRACER_N_CONST(T)

         ! Index of which tracer constituent
         ! will receive the emissions
         IF ( ID_EMITTED(T) > 0 ) THEN 
            IDEMIS(T) = ID_EMITTED(T)
         ENDIF

         ! Loop over all the species which make up the tracer
         DO C = 1, NMEMBER(T)
            
            ! Store tracer coefficient in CTRMB
            CTRMB(T,C) = TRACER_COEFF(T,C) - 1

            ! Loop over all species in "globchem.dat"
            DO J = 1, NSPEC(NCS)

               ! Special case: hydrocarbon tracers as atoms C 
               IF ( TRACER_CONST(T,C) == 'C' ) THEN

                  ! Test SMVGEAR species name against TRACER_NAME
                  IF ( NAMEGAS(J) == TRIM(TRACER_NAME(T)) ) THEN
                     IDTRMB(T,C) = J 
                  ENDIF

               ELSE 

                  ! Test SMVGEAR species name TRACER_CONST
                  IF ( NAMEGAS(J) == TRIM(TRACER_CONST(T,C)) ) THEN
                     IDTRMB(T,C) = J
                  ENDIF

               ENDIF
            ENDDO

#if !defined( ESMF_ )
            ! Print when Debug output selected (hotp 8/5/09)
            IF ( LPRT .and. am_I_Root ) THEN
               PRINT*, '###--------------------'
               PRINT*, '### T, C       : ', T, C
               PRINT*, '### NAME       : ', TRIM( TRACER_NAME(T) )
               PRINT*, '### NMEMBER    : ', NMEMBER(T)
               PRINT*, '### CONST(T,C) : ', TRIM( TRACER_CONST(T,C) )
               PRINT*, '### CTRMB(T,C) : ', CTRMB(T,C)
               PRINT*, '### IDEMIS(T)  : ', IDEMIS(T)
               PRINT*, '### IDTRMB(T,C): ', IDTRMB(T,C)
            ENDIF
#endif
         ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE SETTRACE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fill_chem_state_name_ids
!
! !DESCRIPTION: Subroutine FILL_State_Chm_NAME_IDs populates the trac_name,
!  trac_id, and smvg_id fields of an object based on the CHEMSTATE derived
!  type.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE FILL_CHEM_STATE_NAME_IDs( am_I_Root, State_Chm )
!
! !USES:
!
      ! References to F90 modules
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE TRACER_MOD,         ONLY : N_TRACERS
      USE TRACER_MOD,         ONLY : TRACER_NAME
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root    ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm    ! Obj for chemical state
! 
! !REVISION HISTORY: 
!  13 Aug 2012 - R. Yantosca - Initial version
!  16 Oct 2012 - R. Yantosca - Renamed chem state object to State_Chm
!  19 Oct 2012 - R. Yantosca - Now references gigc_state_chm_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: N

      ! Loop over defined tracers
      DO N = 1, N_TRACERS
         State_Chm%TRAC_ID(N)   = N
         State_Chm%TRAC_NAME(N) = TRIM( TRACER_NAME(N) )
      ENDDO

      END SUBROUTINE FILL_CHEM_STATE_NAME_IDs
!EOC 
!------------------------------------------------------------------------------

      FUNCTION IS_Hg0( N ) RESULT( IT_IS_Hg0 )
!
!******************************************************************************
!  Function IS_Hg0 returns TRUE if tracer N is a total or tagged Hg0 tracer.
!  (cdh, bmy, 12/15/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) N (INTEGER) : GEOS-CHEM tracer number
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: N

      ! Local variables
      LOGICAL             :: IT_IS_Hg0
      INTEGER             :: C

      !=================================================================
      ! IS_Hg0 begins here!
      !=================================================================

      ! Initialize
      IT_IS_Hg0 = .FALSE.

      ! Loop over Hg0 categories
      DO C = 1, N_Hg_CATS
         
         ! Exit with TRUE if corresponds to an Hg0 tracer
         IF ( N == ID_Hg0(C) ) THEN
            IT_IS_Hg0 = .TRUE.
            EXIT
         ENDIF

      ENDDO

      ! Return to calling program
      END FUNCTION IS_Hg0

!------------------------------------------------------------------------------

      FUNCTION IS_Hg2( N ) RESULT( IT_IS_Hg2 )
!
!******************************************************************************
!  Function IS_Hg2 returns TRUE if tracer N is a total or tagged Hg2 tracer.
!  (cdh, bmy, 12/15/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) N (INTEGER) : GEOS-CHEM tracer number
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: N

      ! Local variables
      LOGICAL             :: IT_IS_Hg2
      INTEGER             :: C

      !=================================================================
      ! IS_Hg2 begins here!
      !=================================================================

      ! Initialize
      IT_IS_Hg2 = .FALSE.

      ! Loop over Hg2 categories
      DO C = 1, N_Hg_CATS
         
         ! Exit with TRUE if corresponds to an Hg2 tracer
         IF ( N == ID_Hg2(C) ) THEN
            IT_IS_Hg2 = .TRUE.
            EXIT
         ENDIF

      ENDDO

      ! Return to calling program
      END FUNCTION IS_Hg2

!------------------------------------------------------------------------------

      FUNCTION IS_HgP( N ) RESULT( IT_IS_HgP )
!
!******************************************************************************
!  Function IS_HgP returns TRUE if tracer N is a total or tagged HgP tracer.
!  (cdh, bmy, 12/15/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) N (INTEGER) : GEOS-CHEM tracer number
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: N

      ! Local variables
      LOGICAL             :: IT_IS_HgP
      INTEGER             :: C

      !=================================================================
      ! IS_HgP begins here!
      !=================================================================

      ! Initialize
      IT_IS_HgP = .FALSE.

      ! Loop over Hg2 categories
      DO C = 1, N_Hg_CATS
         
         ! Exit with TRUE if corresponds to an HgP tracer
         IF ( N == ID_HgP(C) ) THEN
            IT_IS_HgP = .TRUE.
            EXIT
         ENDIF

      ENDDO

      ! Return to calling program
      END FUNCTION IS_HgP

!------------------------------------------------------------------------------

      FUNCTION GET_Hg0_CAT( N ) RESULT( NN )
!
!******************************************************************************
!  Function GET_Hg0_CAT the Hg0 category number given the tracer number. 
!  (eck, sas, cdh, bmy, 1/6/05)
!
!  Arguments as Input:
!  ----------------------------------------------------------------------------
!  (1 ) N (INTEGER) : GEOS-CHEM tracer number
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: N
      
      ! Function value
      INTEGER             :: NN

      !=================================================================
      ! GET_Hg0_CAT begins here!
      !=================================================================

      ! Pick the Hg0 category number from the tracer number
      IF ( N == ID_Hg0(ID_Hg_tot) ) THEN 

         ! Total
         NN = ID_Hg_tot

!- eds 8/31/10 -----------------------------------------------
!      ELSE IF ( N == ID_Hg0(ID_Hg_na) ) THEN 
!
!         ! Anthro North America
!         NN = ID_Hg_na
!
!      ELSE IF ( N == ID_Hg0(ID_Hg_eu) ) THEN 
!
!         ! Anthro Europe
!         NN = ID_Hg_eu
!
!      ELSE IF ( N == ID_Hg0(ID_Hg_as) ) THEN 
!
!         ! Anthro Asia
!         NN = ID_Hg_as
!
!      ELSE IF ( N == ID_Hg0(ID_Hg_rw) ) THEN 
!
!         ! Anthro Rest of World
!         NN = ID_Hg_rw
!
!      ELSE IF ( N == ID_Hg0(ID_Hg_oc) ) THEN 
!
!         ! Oceans
!         NN = ID_Hg_oc
!
!      ELSE IF ( N == ID_Hg0(ID_Hg_ln) ) THEN 
!
!         ! Land re-emission
!         NN = ID_Hg_ln
!
!      ELSE IF ( N == ID_Hg0(ID_Hg_nt) ) THEN 
!
!         ! Natural source
!         NN = ID_Hg_nt

      ELSE IF ( N == ID_Hg0(ID_Hg_can) ) THEN  

         ! Anthro Canada
         NN = ID_Hg_can

      ELSE IF ( N == ID_Hg0(ID_Hg_usa) ) THEN

         ! Anthro United States
         NN = ID_Hg_usa

      ELSE IF ( N == ID_Hg0(ID_Hg_cam) ) THEN

         ! Anthro Central America
         NN = ID_Hg_cam

      ELSE IF ( N == ID_Hg0(ID_Hg_sam) ) THEN

         ! Anthro South America
         NN = ID_Hg_sam

      ELSE IF ( N == ID_Hg0(ID_Hg_waf) ) THEN

         ! Anthro West Africa
         NN = ID_Hg_waf

      ELSE IF ( N == ID_Hg0(ID_Hg_eaf) ) THEN

         ! Anthro East Africa
         NN = ID_Hg_eaf

      ELSE IF ( N == ID_Hg0(ID_Hg_saf) ) THEN

         ! Anthro South Africa
         NN = ID_Hg_saf

      ELSE IF ( N == ID_Hg0(ID_Hg_naf) ) THEN

         ! Anthro North Africa
         NN = ID_Hg_naf

      ELSE IF ( N == ID_Hg0(ID_Hg_eur) ) THEN

         ! Anthro OECD Europe
         NN = ID_Hg_eur

      ELSE IF ( N == ID_Hg0(ID_Hg_eeu) ) THEN

         ! Anthro Eastern Europe
         NN = ID_Hg_eeu

      ELSE IF ( N == ID_Hg0(ID_Hg_mde) ) THEN

         ! Anthro Middle East
         NN = ID_Hg_mde

      ELSE IF ( N == ID_Hg0(ID_Hg_sov) ) THEN

         ! Anthro Former Soviet Union
         NN = ID_Hg_sov

      ELSE IF ( N == ID_Hg0(ID_Hg_sas) ) THEN

         ! Anthro South Asia
         NN = ID_Hg_sas

      ELSE IF ( N == ID_Hg0(ID_Hg_eas) ) THEN

         ! Anthro East Asia
         NN = ID_Hg_eas

      ELSE IF ( N == ID_Hg0(ID_Hg_sea) ) THEN

         ! Anthro Southeast Asia
         NN = ID_Hg_sea

      ELSE IF ( N == ID_Hg0(ID_Hg_jpn) ) THEN

         ! Anthro Japan
         NN = ID_Hg_jpn

      ELSE IF ( N == ID_Hg0(ID_Hg_oce) ) THEN

         ! Anthro Oceania
         NN = ID_Hg_oce

      ELSE IF ( N == ID_Hg0(ID_Hg_so) ) THEN

         ! Soil
         NN = ID_Hg_so

      ELSE IF ( N == ID_Hg0(ID_Hg_bb) ) THEN

         ! Biomass burning
         NN = ID_Hg_bb

      ELSE IF ( N == ID_Hg0(ID_Hg_geo) ) THEN

         ! Geogenic
         NN = ID_Hg_geo

      ELSE IF ( N == ID_Hg0(ID_Hg_atl) ) THEN

         ! Below Mixed Layer Mid-Atlantic Ocean
         NN = ID_Hg_atl

      ELSE IF ( N == ID_Hg0(ID_Hg_nat) ) THEN

         ! Below Mixed Layer North Atlantic Ocean
         NN = ID_Hg_nat

      ELSE IF ( N == ID_Hg0(ID_Hg_sat) ) THEN

         ! Below Mixed Layer South Atlantic Ocean
         NN = ID_Hg_sat

      ELSE IF ( N == ID_Hg0(ID_Hg_npa) ) THEN

         ! Below Mixed Layer North Pacific Ocean
         NN = ID_Hg_npa

      ELSE IF ( N == ID_Hg0(ID_Hg_arc) ) THEN

         ! Below Mixed Layer Arctic Ocean
         NN = ID_Hg_arc

      ELSE IF ( N == ID_Hg0(ID_Hg_ant) ) THEN

         ! Below Mixed Layer Antarctic Ocean
         NN = ID_Hg_ant

      ELSE IF ( N == ID_Hg0(ID_Hg_ocn) ) THEN

         ! Below Mixed Layer All Other Ocean
         NN = ID_Hg_ocn

      ELSE IF ( N == ID_Hg0(ID_Hg_str) ) THEN

         ! Residual Stratosphere from Spin-up
         NN = ID_Hg_str

!-------------------------------------------------------

      ELSE

         ! Invalid category
         NN = -1

      ENDIF

      ! Return to calling program
      END FUNCTION GET_Hg0_CAT

!------------------------------------------------------------------------------

      FUNCTION GET_Hg2_CAT( N ) RESULT( NN )
!
!******************************************************************************
!  Function GET_Hg2_CAT the Hg2 category number (i.e. index for DD_Hg2 and
!  WD_Hg2) given the tracer number. (eck, sas, cdh, bmy, 1/6/05)
!
!  Arguments as Input:
!  ----------------------------------------------------------------------------
!  (1 ) N (INTEGER) : GEOS-CHEM tracer number
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: N
      
      ! Function value
      INTEGER             :: NN

      !=================================================================
      ! GET_Hg2_CAT begins here!
      !=================================================================

      ! Pick the Hg2 category number from the tracer number
      IF ( N == ID_Hg2(ID_Hg_tot) ) THEN 

         ! Total
         NN = ID_Hg_tot

!- eds 8/31/10 --------------------------------------------------------------
!      ELSE IF ( N == ID_Hg2(ID_Hg_na) ) THEN 
!
!         ! Anthro North America
!         NN = ID_Hg_na
!
!      ELSE IF ( N == ID_Hg2(ID_Hg_eu) ) THEN 
!
!         ! Anthro Europe
!         NN = ID_Hg_eu
!
!      ELSE IF ( N == ID_Hg2(ID_Hg_as) ) THEN 
!
!         ! Anthro Asia
!         NN = ID_Hg_as
!
!      ELSE IF ( N == ID_Hg2(ID_Hg_rw) ) THEN 
!
!         ! Anthro Rest of World
!         NN = ID_Hg_rw
!
!      ELSE IF ( N == ID_Hg2(ID_Hg_oc) ) THEN 
!
!         ! Oceans
!         NN = ID_Hg_oc
!
!      ELSE IF ( N == ID_Hg2(ID_Hg_ln) ) THEN 
!
!         ! Land re-emission
!         NN = ID_Hg_ln
!
!      ELSE IF ( N == ID_Hg2(ID_Hg_nt) ) THEN 
!
!         ! Natural source
!         NN = ID_Hg_nt

      ! New Tagged Hg Simulation (eds 8/31/10)
      ELSE IF ( N == ID_Hg2(ID_Hg_can) ) THEN  

         ! Anthro Canada
         NN = ID_Hg_can

      ELSE IF ( N == ID_Hg2(ID_Hg_usa) ) THEN

         ! Anthro United States
         NN = ID_Hg_usa

      ELSE IF ( N == ID_Hg2(ID_Hg_cam) ) THEN

         ! Anthro Central America
         NN = ID_Hg_cam

      ELSE IF ( N == ID_Hg2(ID_Hg_sam) ) THEN

         ! Anthro South America
         NN = ID_Hg_sam

      ELSE IF ( N == ID_Hg2(ID_Hg_waf) ) THEN

         ! Anthro West Africa
         NN = ID_Hg_waf

      ELSE IF ( N == ID_Hg2(ID_Hg_eaf) ) THEN

         ! Anthro East Africa
         NN = ID_Hg_eaf

      ELSE IF ( N == ID_Hg2(ID_Hg_saf) ) THEN

         ! Anthro South Africa
         NN = ID_Hg_saf

      ELSE IF ( N == ID_Hg2(ID_Hg_naf) ) THEN

         ! Anthro North Africa
         NN = ID_Hg_naf

      ELSE IF ( N == ID_Hg2(ID_Hg_eur) ) THEN

         ! Anthro OECD Europe
         NN = ID_Hg_eur

      ELSE IF ( N == ID_Hg2(ID_Hg_eeu) ) THEN

         ! Anthro Eastern Europe
         NN = ID_Hg_eeu

      ELSE IF ( N == ID_Hg2(ID_Hg_mde) ) THEN

         ! Anthro Middle East
         NN = ID_Hg_mde

      ELSE IF ( N == ID_Hg2(ID_Hg_sov) ) THEN

         ! Anthro Former Soviet Union
         NN = ID_Hg_sov

      ELSE IF ( N == ID_Hg2(ID_Hg_sas) ) THEN

         ! Anthro South Asia
         NN = ID_Hg_sas

      ELSE IF ( N == ID_Hg2(ID_Hg_eas) ) THEN

         ! Anthro East Asia
         NN = ID_Hg_eas

      ELSE IF ( N == ID_Hg2(ID_Hg_sea) ) THEN

         ! Anthro Southeast Asia
         NN = ID_Hg_sea

      ELSE IF ( N == ID_Hg2(ID_Hg_jpn) ) THEN

         ! Anthro Japan
         NN = ID_Hg_jpn

      ELSE IF ( N == ID_Hg2(ID_Hg_oce) ) THEN

         ! Anthro Oceania
         NN = ID_Hg_oce

      ELSE IF ( N == ID_Hg2(ID_Hg_so) ) THEN

         ! Soil
         NN = ID_Hg_so

      ELSE IF ( N == ID_Hg2(ID_Hg_bb) ) THEN

         ! Biomass burning
         NN = ID_Hg_bb

      ELSE IF ( N == ID_Hg2(ID_Hg_geo) ) THEN

         ! Geogenic
         NN = ID_Hg_geo

      ELSE IF ( N == ID_Hg2(ID_Hg_atl) ) THEN

         ! Below Mixed Layer Mid-Atlantic Ocean
         NN = ID_Hg_atl

      ELSE IF ( N == ID_Hg2(ID_Hg_nat) ) THEN

         ! Below Mixed Layer North Atlantic Ocean
         NN = ID_Hg_nat

      ELSE IF ( N == ID_Hg2(ID_Hg_sat) ) THEN

         ! Below Mixed Layer South Atlantic Ocean
         NN = ID_Hg_sat

      ELSE IF ( N == ID_Hg2(ID_Hg_npa) ) THEN

         ! Below Mixed Layer North Pacific Ocean
         NN = ID_Hg_npa

      ELSE IF ( N == ID_Hg2(ID_Hg_arc) ) THEN

         ! Below Mixed Layer Arctic Ocean
         NN = ID_Hg_arc

      ELSE IF ( N == ID_Hg2(ID_Hg_ant) ) THEN

         ! Below Mixed Layer Antarctic Ocean
         NN = ID_Hg_ant

      ELSE IF ( N == ID_Hg2(ID_Hg_ocn) ) THEN

         ! Below Mixed Layer All Other Ocean
         NN = ID_Hg_ocn

      ELSE IF ( N == ID_Hg2(ID_Hg_str) ) THEN

         ! Residual Stratosphere from Spin-up
         NN = ID_Hg_str

!-------------------------------------------------------

      ELSE

         ! Invalid category
         NN = -1

      ENDIF

      ! Return to calling program
      END FUNCTION GET_Hg2_CAT

!------------------------------------------------------------------------------

      FUNCTION GET_HgP_CAT( N ) RESULT( NN )
!
!******************************************************************************
!  Function GET_HgP_CAT the HgP category number (i.e. index for DD_HgP and
!  WD_HgP) given the tracer number. (eck, sas, cdh, bmy, 1/6/05)
!
!  Arguments as Input:
!  ----------------------------------------------------------------------------
!  (1 ) N (INTEGER) : GEOS-CHEM tracer number
!
!  NOTES:
!  17 May 2010 - C. Carouge   - Added to standard version
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: N
      
      ! Function value
      INTEGER             :: NN

      !=================================================================
      ! GET_HgP_CAT begins here!
      !=================================================================

      ! Pick the HgP category number from the tracer number
      IF ( N == ID_HgP(ID_Hg_tot) ) THEN 

         ! Total
         NN = ID_Hg_tot

!- eds 8/31/10 ---------------------------------------------------------
!      ELSE IF ( N == ID_HgP(ID_Hg_na) ) THEN 
!
!         ! Anthro North America
!         NN = ID_Hg_na
!
!      ELSE IF ( N == ID_HgP(ID_Hg_eu) ) THEN 
!
!         ! Anthro Europe
!         NN = ID_Hg_eu
!
!      ELSE IF ( N == ID_HgP(ID_Hg_as) ) THEN 
!
!         ! Anthro Asia
!         NN = ID_Hg_as
!
!      ELSE IF ( N == ID_HgP(ID_Hg_rw) ) THEN 
!
!         ! Anthro Rest of World
!         NN = ID_Hg_rw
!
!      ELSE IF ( N == ID_HgP(ID_Hg_oc) ) THEN 
!
!         ! Oceans
!         NN = ID_Hg_oc
!
!      ELSE IF ( N == ID_HgP(ID_Hg_ln) ) THEN 
!
!         ! Land re-emission
!         NN = ID_Hg_ln
!
!      ELSE IF ( N == ID_HgP(ID_Hg_nt) ) THEN 
!
!         ! Natural source
!         NN = ID_Hg_nt

      ELSE IF ( N == ID_HgP(ID_Hg_can) ) THEN  

         ! Anthro Canada
         NN = ID_Hg_can

      ELSE IF ( N == ID_HgP(ID_Hg_usa) ) THEN

         ! Anthro United States
         NN = ID_Hg_usa

      ELSE IF ( N == ID_HgP(ID_Hg_cam) ) THEN

         ! Anthro Central America
         NN = ID_Hg_cam

      ELSE IF ( N == ID_HgP(ID_Hg_sam) ) THEN

         ! Anthro South America
         NN = ID_Hg_sam

      ELSE IF ( N == ID_HgP(ID_Hg_waf) ) THEN

         ! Anthro West Africa
         NN = ID_Hg_waf

      ELSE IF ( N == ID_HgP(ID_Hg_eaf) ) THEN

         ! Anthro East Africa
         NN = ID_Hg_eaf

      ELSE IF ( N == ID_HgP(ID_Hg_saf) ) THEN

         ! Anthro South Africa
         NN = ID_Hg_saf

      ELSE IF ( N == ID_HgP(ID_Hg_naf) ) THEN

         ! Anthro North Africa
         NN = ID_Hg_naf

      ELSE IF ( N == ID_HgP(ID_Hg_eur) ) THEN

         ! Anthro OECD Europe
         NN = ID_Hg_eur

      ELSE IF ( N == ID_HgP(ID_Hg_eeu) ) THEN

         ! Anthro Eastern Europe
         NN = ID_Hg_eeu

      ELSE IF ( N == ID_HgP(ID_Hg_mde) ) THEN

         ! Anthro Middle East
         NN = ID_Hg_mde

      ELSE IF ( N == ID_HgP(ID_Hg_sov) ) THEN

         ! Anthro Former Soviet Union
         NN = ID_Hg_sov

      ELSE IF ( N == ID_HgP(ID_Hg_sas) ) THEN

         ! Anthro South Asia
         NN = ID_Hg_sas

      ELSE IF ( N == ID_HgP(ID_Hg_eas) ) THEN

         ! Anthro East Asia
         NN = ID_Hg_eas

      ELSE IF ( N == ID_HgP(ID_Hg_sea) ) THEN

         ! Anthro Southeast Asia
         NN = ID_Hg_sea

      ELSE IF ( N == ID_HgP(ID_Hg_jpn) ) THEN

         ! Anthro Japan
         NN = ID_Hg_jpn

      ELSE IF ( N == ID_HgP(ID_Hg_oce) ) THEN

         ! Anthro Oceania
         NN = ID_Hg_oce

      ELSE IF ( N == ID_HgP(ID_Hg_so) ) THEN

         ! Soil
         NN = ID_Hg_so

      ELSE IF ( N == ID_HgP(ID_Hg_bb) ) THEN

         ! Biomass burning
         NN = ID_Hg_bb

      ELSE IF ( N == ID_HgP(ID_Hg_geo) ) THEN

         ! Geogenic
         NN = ID_Hg_geo

      ELSE IF ( N == ID_HgP(ID_Hg_atl) ) THEN

         ! Below Mixed Layer Mid-Atlantic Ocean
         NN = ID_Hg_atl

      ELSE IF ( N == ID_HgP(ID_Hg_nat) ) THEN

         ! Below Mixed Layer North Atlantic Ocean
         NN = ID_Hg_nat

      ELSE IF ( N == ID_HgP(ID_Hg_sat) ) THEN

         ! Below Mixed Layer South Atlantic Ocean
         NN = ID_Hg_sat

      ELSE IF ( N == ID_HgP(ID_Hg_npa) ) THEN

         ! Below Mixed Layer North Pacific Ocean
         NN = ID_Hg_npa

      ELSE IF ( N == ID_HgP(ID_Hg_arc) ) THEN

         ! Below Mixed Layer Arctic Ocean
         NN = ID_Hg_arc

      ELSE IF ( N == ID_HgP(ID_Hg_ant) ) THEN

         ! Below Mixed Layer Antarctic Ocean
         NN = ID_Hg_ant

      ELSE IF ( N == ID_HgP(ID_Hg_ocn) ) THEN

         ! Below Mixed Layer All Other Ocean
         NN = ID_Hg_ocn

      ELSE IF ( N == ID_HgP(ID_Hg_str) ) THEN

         ! Residual Stratosphere from Spin-up
         NN = ID_Hg_str

!-------------------------------------------------------

      ELSE

         ! Invalid category
         NN = -1

      ENDIF

      ! Return to calling program
      END FUNCTION GET_HgP_CAT

!------------------------------------------------------------------------------

      SUBROUTINE INIT_TRACERID
!
!******************************************************************************
!  Subroutine INIT_TRACERID zeroes module variables. (bmy, 11/12/02, 9/18/07)
!
!  NOTES:
!  (1 ) Now also zero IDDMS, IDSO2, IDSO4, IDMSA (rjp, bmy, 3/23/03)
!  (2 ) Now zero extra flags for carbon & dust tracers (rjp, tdf, bmy, 4/1/04)
!  (3 ) Now zero extra flags for seasalt tracers (rjp, bec, bmy, 4/1/04)
!  (4 ) Now zero extra flags for SOA tracers (rjp, bmy, 7/13/04)
!  (5 ) Now zero IDTHG0, IDTHG2, IDTHGP + tagged Hg's (eck, bmy, 12/7/04)
!  (6 ) Now zero IDTAS, IDTAHS, IDTLET, IDTNH4aq, IDTSO4aq (cas, bmy, 12/20/04)
!  (7 ) Now allocate ID_Hg0, ID_Hg2, ID_HgP (bmy, 12/16/05)
!  (8 ) Now zero IDTSOG4, IDTSOA4 (dkh, bmy, 5/18/06)
!  (9 ) Now zero IDTH2, IDTHD (hup, lyj, phs, 9/18/07)
!  (10) Now zero IDEHNO3 (PHS, 3/4/08)
!  (11) Zero IDTNK1-30, IDTSF1-30, IDTSS1-30, IDTAW1-30, IDTH2SO4, 
!       IDTECOC1-30, IDTECIL1-30, IDTOCOB1-30, IDTOCIL1-30 and 
!       IDTDUST1..IDTDUST30 (win, 6/23/09)
!  (12) Comment out IDBENZ, IDTOLU, IDXYLE, IDMONX, IDGLYX, IDMGLY
!        (ccc, 01/29/10)
!  (13) Now zero IDs for Henze's aromatics (dkh, 01/10)
!  (14) Now zero IDTs for new tracers in isoprene chemistry, Henze's aromatics,
!       (fp, dkh, hotp, 01/10)
!  (15) Now zero IDBs for all biomass burning emissions. (fp, hotp, 01/10)
!  07 Feb 2013 - M. Payer    - Add NO, NO2, NO2, HNO2, and O3 as tracers to
!                              replace NOx and Ox
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations (H. Pye)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE LOGICAL_MOD, ONLY : LSPLIT
      USE TRACER_MOD,  ONLY : ITS_A_MERCURY_SIM
      
      ! Local variables
      INTEGER              :: AS

      ! SMVGEAR species ID #'s
      IDO3      = 0
      IDNO2     = 0
      IDNO3     = 0 
      IDN2O5    = 0
      IDHNO4    = 0
      IDHC1     = 0
      IDNO      = 0
      IDHNO2    = 0
      IDCO      = 0
      IDPRPE    = 0
      IDISOP    = 0
      IDALK4    = 0
      IDC3H8    = 0
      IDPAN     = 0
      IDGLPAN   = 0
      IDGPAN    = 0
      IDPMN     = 0
      IDPPN     = 0
      IDHNO3    = 0
      IDOH      = 0
      IDHO2     = 0
      IDH2O2    = 0
      IDACET    = 0
      IDMEK     = 0
      IDALD2    = 0
      IDRCHO    = 0 
      IDMVK     = 0
      IDMACR    = 0
      IDISN2    = 0
      IDR4N2    = 0
      IDCH2O    = 0
      IDC2H6    = 0
      IDMP      = 0
      IDDMS     = 0
      IDSO2     = 0
      IDSO4     = 0
      IDMSA     = 0
      IDDRYO3   = 0
      IDDRYPAN  = 0
      IDDRYNO2  = 0       
!- Not used (ccc, 01/27/10)
!      IDBENZ    = 0
!      IDTOLU    = 0
!      IDXYLE    = 0
!      IDMONX    = 0
!      IDGLYX    = 0
!      IDMGLY    = 0
      IDDRYGLYX = 0
      IDDRYMGLY = 0
      IDC2H4    = 0
      IDC2H2    = 0
      IDMBO     = 0
      IDGLYC    = 0
      IDHAC     = 0
      IDMPN     = 0

      IDAPAN    = 0
      IDENPAN   = 0
      IDMPAN    = 0
      IDNIPAN   = 0

      IDDRYAPAN  = 0
      IDDRYENPAN = 0
      IDDRYGLPAN = 0
      IDDRYGPAN  = 0
      IDDRYMPAN  = 0
      IDDRYNIPAN = 0
      ! added for aroms (dkh, 10/06/06)  
      IDLBRO2H  = 0
      IDLBRO2N  = 0
      IDLTRO2H  = 0
      IDLTRO2N  = 0
      IDLXRO2H  = 0
      IDLXRO2N  = 0
      IDH2SO4   = 0 !(win, 6/23/09)
                     !(sfarina - leaving this for all sims)
      ! SOAupdate: for gas phase NAP chemistry (hotp 6/24/09)
      ! semivolpoa5: IVOCs (hotp 3/30/09)
      IDNAP     = 0
      IDNRO2    = 0
      IDLNRO2H  = 0
      IDLNRO2N  = 0

      ! +++++++++++++++++
      !jpp, 6/13/07,
      IDBr2     = 0
      IDBr      = 0
      IDBrO     = 0
      IDHBr     = 0
      IDHOBr    = 0
      IDBrNO2   = 0
      IDBrNO3   = 0
      IDCHBr3   = 0
      IDCH2Br2  = 0
      IDCH3Br   = 0
      !jpp, 2/27/08
      !ID's for drydep
      IDDRYHOBr = 0
      IDDRYHBr  = 0
      IDDRYBrNO3= 0
      ! +++++++++++++++++

      ! GEOS-CHEM Tracer ID #'s
      IDTNO     = 0
      IDTNO2    = 0
      IDTNO3    = 0
      IDTHNO2   = 0
      IDTO3     = 0
      IDTPAN    = 0
      IDTCO     = 0
      IDTH2     = 0 ! (hup, 7/14/2004)
      IDTHD     = 0 ! (jaegle, 11/07/2005)
      IDTALK4   = 0
      IDTISOP   = 0
      IDTHNO3   = 0
      IDTH2O2   = 0
      IDTACET   = 0
      IDTMEK    = 0
      IDTALD2   = 0
      IDTRCHO   = 0
      IDTMVK    = 0
      IDTMACR   = 0
      IDTPMN    = 0
      IDTPPN    = 0
      IDTISN2   = 0
      IDTR4N2   = 0
      IDTPRPE   = 0
      IDTC3H8   = 0
      IDTCH2O   = 0
      IDTC2H6   = 0
      IDTN2O5   = 0
      IDTHNO4   = 0
      IDTMP     = 0
      IDTDMS    = 0 
      IDTSO2    = 0 
      IDTSO4    = 0    
      IDTMSA    = 0 
      IDTNH3    = 0 
      IDTNH4    = 0 
      IDTNIT    = 0 
      IDTAS     = 0
      IDTAHS    = 0
      IDTNH4aq  = 0
      IDTLET    = 0
      IDTSO4aq  = 0
      IDTBCPI   = 0
      IDTOCPI   = 0
      IDTBCPO   = 0
      IDTOCPO   = 0
      IDTO3Strt = 0
      IDTDST1   = 0
      IDTDST2   = 0
      IDTDST3   = 0
      IDTDST4   = 0
      IDTSALA   = 0
      IDTSALC   = 0
      IDTRN     = 0
      IDTPB     = 0
      IDTBE7    = 0
      IDTGLYX   = 0
      IDTMGLY   = 0
      IDTBENZ   = 0
      IDTTOLU   = 0
      IDTXYLE   = 0
      IDTMONX   = 0
      IDTSOAG   = 0
      IDTSOAM   = 0
      IDTC2H4   = 0
      IDTC2H2   = 0
      IDTMBO    = 0
      IDTGLYC   = 0
      IDTHAC    = 0
      IDTAPAN   = 0
      IDTENPAN  = 0
      IDTGLPAN  = 0
      IDTGPAN   = 0
      IDTMPAN   = 0
      IDTNIPAN  = 0
      IDTMPN    = 0
      ! +++++++++++++++++
      ! jpp for bromine
      ! 6/5/09
      IDTBr2     = 0
      IDTBr      = 0
      IDTBrO     = 0
      IDTHBr     = 0
      IDTHOBr    = 0
      IDTBrNO2   = 0
      IDTBrNO3   = 0
      IDTCHBr3   = 0
      IDTCH2Br2  = 0
      IDTCH3Br   = 0
      ! +++++++++++++++++

      !(fp, 6/09)
      IDTISOPN   = 0
      IDTAP      = 0
      IDTPROPNN  = 0
      IDTMOBA    = 0
      IDTMMN     = 0
      IDTRIP     = 0
      IDTIEPOX   = 0
      IDTMAP     = 0
      IDTPYPAN   = 0

      ! SOAupdate: convert all POA/POG to 1 and 2 (hotp 10/11/09)
       !semivolpoa (hotp 8/24/09)
      IDTPOA1   = 0
      IDTPOA2   = 0
      ! SOAupdate: semivol2: add POG (hotp 2/27/09)
      IDTPOG1   = 0
      IDTPOG2   = 0
      ! SOAupdate: semivolpoa4opoa: add OPOA, OPOG (hotp 3/17/09)
      IDTOPOA1  = 0
      IDTOPOA2  = 0
      IDTOPOG1  = 0
      IDTOPOG2  = 0
      ! SOAupdate: for gas phase NAP chemsitry (hotp 6/24/09)
      ! semivolpoa5: IVOC (hotp 3/30/09)
      IDTNAP    = 0
      ! SOAudpate: new lumped aromatic and IVOC (hotp 5/12/10)
      IDTASOAN  = 0
      IDTASOA1  = 0
      IDTASOA2  = 0
      IDTASOA3  = 0
      IDTASOG1  = 0
      IDTASOG2  = 0
      IDTASOG3  = 0
      ! SOAupdate: new mtp: new terpene and isoprene IDTs (hotp 5/20/10)
      IDTMTPA   = 0
      IDTLIMO   = 0
      IDTMTPO   = 0
      IDTTSOG1  = 0
      IDTTSOG2  = 0
      IDTTSOG3  = 0
      IDTTSOG0  = 0
      IDTTSOA1  = 0
      IDTTSOA2  = 0
      IDTTSOA3  = 0
      IDTTSOA0  = 0
      IDTISOG1  = 0
      IDTISOG2  = 0
      IDTISOG3  = 0
      IDTISOA1  = 0
      IDTISOA2  = 0
      IDTISOA3  = 0
      ! End SOAupdate 

#if   defined( TOMAS )
      IDTH2SO4  = 0 !(win, 6/23/09)


                  IDTNK1  = 0
                  IDTNK2  = 0
                  IDTNK3  = 0
                  IDTNK4  = 0
                  IDTNK5  = 0
                  IDTNK6  = 0
                  IDTNK7  = 0
                  IDTNK8  = 0
                  IDTNK9  = 0
                  IDTNK10 = 0
                  IDTNK11 = 0
                  IDTNK12 = 0
# if  !defined( TOMAS12 )
                  IDTNK13 = 0
                  IDTNK14 = 0
                  IDTNK15 = 0
#  if !defined( TOMAS15 )
                  IDTNK16 = 0
                  IDTNK17 = 0
                  IDTNK18 = 0
                  IDTNK19 = 0
                  IDTNK20 = 0
                  IDTNK21 = 0
                  IDTNK22 = 0
                  IDTNK23 = 0
                  IDTNK24 = 0
                  IDTNK25 = 0
                  IDTNK26 = 0
                  IDTNK27 = 0
                  IDTNK28 = 0
                  IDTNK29 = 0
                  IDTNK30 = 0
#  endif
# endif
# if   defined( TOMAS40 )
                  IDTNK31 = 0
                  IDTNK32 = 0
                  IDTNK33 = 0
                  IDTNK34 = 0
                  IDTNK35 = 0
                  IDTNK36 = 0
                  IDTNK37 = 0
                  IDTNK38 = 0
                  IDTNK39 = 0
                  IDTNK40 = 0
# endif

                  IDTSF1  = 0
                  IDTSF2  = 0
                  IDTSF3  = 0
                  IDTSF4  = 0
                  IDTSF5  = 0
                  IDTSF6  = 0
                  IDTSF7  = 0
                  IDTSF8  = 0
                  IDTSF9  = 0
                  IDTSF10 = 0
                  IDTSF11 = 0
                  IDTSF12 = 0
# if  !defined( TOMAS12 )
                  IDTSF13 = 0
                  IDTSF14 = 0
                  IDTSF15 = 0
#  if !defined( TOMAS15 )
                  IDTSF16 = 0
                  IDTSF17 = 0
                  IDTSF18 = 0
                  IDTSF19 = 0
                  IDTSF20 = 0
                  IDTSF21 = 0
                  IDTSF22 = 0
                  IDTSF23 = 0
                  IDTSF24 = 0
                  IDTSF25 = 0
                  IDTSF26 = 0
                  IDTSF27 = 0
                  IDTSF28 = 0
                  IDTSF29 = 0
                  IDTSF30 = 0
#  endif
# endif
# if   defined( TOMAS40 )
                  IDTSF31 = 0
                  IDTSF32 = 0
                  IDTSF33 = 0
                  IDTSF34 = 0
                  IDTSF35 = 0
                  IDTSF36 = 0
                  IDTSF37 = 0
                  IDTSF38 = 0
                  IDTSF39 = 0
                  IDTSF40 = 0
# endif

                  IDTSS1  = 0
                  IDTSS2  = 0
                  IDTSS3  = 0
                  IDTSS4  = 0
                  IDTSS5  = 0
                  IDTSS6  = 0
                  IDTSS7  = 0
                  IDTSS8  = 0
                  IDTSS9  = 0
                  IDTSS10 = 0
                  IDTSS11 = 0
                  IDTSS12 = 0
# if  !defined( TOMAS12 )
                  IDTSS13 = 0
                  IDTSS14 = 0
                  IDTSS15 = 0
#  if !defined( TOMAS15 )
                  IDTSS16 = 0
                  IDTSS17 = 0
                  IDTSS18 = 0
                  IDTSS19 = 0
                  IDTSS20 = 0
                  IDTSS21 = 0
                  IDTSS22 = 0
                  IDTSS23 = 0
                  IDTSS24 = 0
                  IDTSS25 = 0
                  IDTSS26 = 0
                  IDTSS27 = 0
                  IDTSS28 = 0
                  IDTSS29 = 0
                  IDTSS30 = 0
#  endif
# endif
# if   defined( TOMAS40 )
                  IDTSS31 = 0
                  IDTSS32 = 0
                  IDTSS33 = 0
                  IDTSS34 = 0
                  IDTSS35 = 0
                  IDTSS36 = 0
                  IDTSS37 = 0
                  IDTSS38 = 0
                  IDTSS39 = 0
                  IDTSS40 = 0
# endif

                IDTECOB1  = 0
                IDTECOB2  = 0
                IDTECOB3  = 0
                IDTECOB4  = 0
                IDTECOB5  = 0
                IDTECOB6  = 0
                IDTECOB7  = 0
                IDTECOB8  = 0
                IDTECOB9  = 0
                IDTECOB10 = 0
                IDTECOB11 = 0
                IDTECOB12 = 0
# if  !defined( TOMAS12 )
                IDTECOB13 = 0
                IDTECOB14 = 0
                IDTECOB15 = 0
#  if !defined( TOMAS15 )
                IDTECOB16 = 0
                IDTECOB17 = 0
                IDTECOB18 = 0
                IDTECOB19 = 0
                IDTECOB20 = 0
                IDTECOB21 = 0
                IDTECOB22 = 0
                IDTECOB23 = 0
                IDTECOB24 = 0
                IDTECOB25 = 0
                IDTECOB26 = 0
                IDTECOB27 = 0
                IDTECOB28 = 0
                IDTECOB29 = 0
                IDTECOB30 = 0
#  endif
# endif
# if   defined( TOMAS40 )
                IDTECOB31 = 0
                IDTECOB32 = 0
                IDTECOB33 = 0
                IDTECOB34 = 0
                IDTECOB35 = 0
                IDTECOB36 = 0
                IDTECOB37 = 0
                IDTECOB38 = 0
                IDTECOB39 = 0
                IDTECOB40 = 0
# endif

                IDTECIL1  = 0
                IDTECIL2  = 0
                IDTECIL3  = 0
                IDTECIL4  = 0
                IDTECIL5  = 0
                IDTECIL6  = 0
                IDTECIL7  = 0
                IDTECIL8  = 0
                IDTECIL9  = 0
                IDTECIL10 = 0
                IDTECIL11 = 0
                IDTECIL12 = 0
# if  !defined( TOMAS12 )
                IDTECIL13 = 0
                IDTECIL14 = 0
                IDTECIL15 = 0
#  if !defined( TOMAS15 )
                IDTECIL16 = 0
                IDTECIL17 = 0
                IDTECIL18 = 0
                IDTECIL19 = 0
                IDTECIL20 = 0
                IDTECIL21 = 0
                IDTECIL22 = 0
                IDTECIL23 = 0
                IDTECIL24 = 0
                IDTECIL25 = 0
                IDTECIL26 = 0
                IDTECIL27 = 0
                IDTECIL28 = 0
                IDTECIL29 = 0
                IDTECIL30 = 0
#  endif
# endif
# if   defined( TOMAS40 )
                IDTECIL31 = 0
                IDTECIL32 = 0
                IDTECIL33 = 0
                IDTECIL34 = 0
                IDTECIL35 = 0
                IDTECIL36 = 0
                IDTECIL37 = 0
                IDTECIL38 = 0
                IDTECIL39 = 0
                IDTECIL40 = 0
# endif

                IDTOCOB1  = 0
                IDTOCOB2  = 0
                IDTOCOB3  = 0
                IDTOCOB4  = 0
                IDTOCOB5  = 0
                IDTOCOB6  = 0
                IDTOCOB7  = 0
                IDTOCOB8  = 0
                IDTOCOB9  = 0
                IDTOCOB10 = 0
                IDTOCOB11 = 0
                IDTOCOB12 = 0
# if  !defined( TOMAS12 )
                IDTOCOB13 = 0
                IDTOCOB14 = 0
                IDTOCOB15 = 0
#  if !defined( TOMAS15 )
                IDTOCOB16 = 0
                IDTOCOB17 = 0
                IDTOCOB18 = 0
                IDTOCOB19 = 0
                IDTOCOB20 = 0
                IDTOCOB21 = 0
                IDTOCOB22 = 0
                IDTOCOB23 = 0
                IDTOCOB24 = 0
                IDTOCOB25 = 0
                IDTOCOB26 = 0
                IDTOCOB27 = 0
                IDTOCOB28 = 0
                IDTOCOB29 = 0
                IDTOCOB30 = 0
#  endif
# endif
# if   defined( TOMAS40 )
                IDTOCOB31 = 0
                IDTOCOB32 = 0
                IDTOCOB33 = 0
                IDTOCOB34 = 0
                IDTOCOB35 = 0
                IDTOCOB36 = 0
                IDTOCOB37 = 0
                IDTOCOB38 = 0
                IDTOCOB39 = 0
                IDTOCOB40 = 0
# endif

                IDTOCIL1  = 0
                IDTOCIL2  = 0
                IDTOCIL3  = 0
                IDTOCIL4  = 0
                IDTOCIL5  = 0
                IDTOCIL6  = 0
                IDTOCIL7  = 0
                IDTOCIL8  = 0
                IDTOCIL9  = 0
                IDTOCIL10 = 0
                IDTOCIL11 = 0
                IDTOCIL12 = 0
# if  !defined( TOMAS12 )
                IDTOCIL13 = 0
                IDTOCIL14 = 0
                IDTOCIL15 = 0
#  if !defined( TOMAS15 )
                IDTOCIL16 = 0
                IDTOCIL17 = 0
                IDTOCIL18 = 0
                IDTOCIL19 = 0
                IDTOCIL20 = 0
                IDTOCIL21 = 0
                IDTOCIL22 = 0
                IDTOCIL23 = 0
                IDTOCIL24 = 0
                IDTOCIL25 = 0
                IDTOCIL26 = 0
                IDTOCIL27 = 0
                IDTOCIL28 = 0
                IDTOCIL29 = 0
                IDTOCIL30 = 0
#  endif
# endif
# if   defined( TOMAS40 )
                IDTOCIL31 = 0
                IDTOCIL32 = 0
                IDTOCIL33 = 0
                IDTOCIL34 = 0
                IDTOCIL35 = 0
                IDTOCIL36 = 0
                IDTOCIL37 = 0
                IDTOCIL38 = 0
                IDTOCIL39 = 0
                IDTOCIL40 = 0
# endif

                IDTDUST1  = 0
                IDTDUST2  = 0
                IDTDUST3  = 0
                IDTDUST4  = 0
                IDTDUST5  = 0
                IDTDUST6  = 0
                IDTDUST7  = 0
                IDTDUST8  = 0
                IDTDUST9  = 0
                IDTDUST10 = 0
                IDTDUST11 = 0
                IDTDUST12 = 0
# if  !defined( TOMAS12 )
                IDTDUST13 = 0
                IDTDUST14 = 0
                IDTDUST15 = 0
#  if !defined( TOMAS15 )
                IDTDUST16 = 0
                IDTDUST17 = 0
                IDTDUST18 = 0
                IDTDUST19 = 0
                IDTDUST20 = 0
                IDTDUST21 = 0
                IDTDUST22 = 0
                IDTDUST23 = 0
                IDTDUST24 = 0
                IDTDUST25 = 0
                IDTDUST26 = 0
                IDTDUST27 = 0
                IDTDUST28 = 0
                IDTDUST29 = 0
                IDTDUST30 = 0
#  endif
# endif
# if   defined( TOMAS40 )
                IDTDUST31 = 0
                IDTDUST32 = 0
                IDTDUST33 = 0
                IDTDUST34 = 0
                IDTDUST35 = 0
                IDTDUST36 = 0
                IDTDUST37 = 0
                IDTDUST38 = 0
                IDTDUST39 = 0
                IDTDUST40 = 0
# endif

                  IDTAW1  = 0
                  IDTAW2  = 0
                  IDTAW3  = 0
                  IDTAW4  = 0
                  IDTAW5  = 0
                  IDTAW6  = 0
                  IDTAW7  = 0
                  IDTAW8  = 0
                  IDTAW9  = 0
                  IDTAW10 = 0
                  IDTAW11 = 0
                  IDTAW12 = 0
# if  !defined( TOMAS12 )
                  IDTAW13 = 0
                  IDTAW14 = 0
                  IDTAW15 = 0
#  if !defined( TOMAS15 )
                  IDTAW16 = 0
                  IDTAW17 = 0
                  IDTAW18 = 0
                  IDTAW19 = 0
                  IDTAW20 = 0
                  IDTAW21 = 0
                  IDTAW22 = 0
                  IDTAW23 = 0
                  IDTAW24 = 0
                  IDTAW25 = 0
                  IDTAW26 = 0
                  IDTAW27 = 0
                  IDTAW28 = 0
                  IDTAW29 = 0
                  IDTAW30 = 0
#  endif
# endif
# if   defined( TOMAS40 )
                  IDTAW31 = 0
                  IDTAW32 = 0
                  IDTAW33 = 0
                  IDTAW34 = 0
                  IDTAW35 = 0
                  IDTAW36 = 0
                  IDTAW37 = 0
                  IDTAW38 = 0
                  IDTAW39 = 0
                  IDTAW40 = 0
# endif

#endif

      ! GEOS-CHEM Emission ID #'s
      NEMANTHRO = 0
      NEMBIOG   = 0
      IDENO     = 0
      IDEO3     = 0
      IDECO     = 0
      IDEPRPE   = 0
      IDEC3H8   = 0
      IDEALK4   = 0
      IDEC2H6   = 0
      IDEACET   = 0
      IDEMEK    = 0
      IDEALD2   = 0
      IDEISOP   = 0
      IDECH2O   = 0 
      IDEHNO3   = 0 !phs (3/4/08)
      IDEGLYX   = 0
      IDEMGLY   = 0
      IDEBENZ   = 0
      IDETOLU   = 0
      IDEXYLE   = 0
      IDEMONX   = 0
      IDEC2H4   = 0
      IDEC2H2   = 0
      IDEMBO    = 0
      IDEGLYC   = 0
      IDEHAC    = 0
      IDECHBr3  = 0 ! jpp, 8/30/07
      IDECH2Br2 = 0 ! jpp, 9/16/07
      IDEBr2    = 0 ! jpp, 3/2/10
      ! SOAudpate: for gas phase NAP chemistry (hotp 6/24/09)
      ! semivolpoa5: IVOCs (hotp 3/30/09)
      IDENAP    = 0

      ! GEOS-CHEM Biofuel ID #'s
      IDBFNO    = 0
      IDBFCO    = 0
      IDBFALK4  = 0
      IDBFACET  = 0
      IDBFMEK   = 0
      IDBFALD2  = 0
      IDBFPRPE  = 0 
      IDBFC3H8  = 0
      IDBFCH2O  = 0
      IDBFC2H6  = 0
      IDBFBENZ  = 0
      IDBFTOLU  = 0
      IDBFXYLE  = 0
      IDBFC2H2  = 0
      IDBFC2H4  = 0
      IDBFGLYC  = 0
      IDBFGLYX  = 0
      IDBFMGLY  = 0
      IDBFHAC   = 0
      ! SOAudpate: for gas phase NAP chemistry (hotp 6/24/09)
      ! semivolpoa5: IVOCs (hotp 3/30/09)
      IDBFNAP   = 0

      ! Initialize IDBs for BIOMASS BURNING (hotp 8/24/09)
      IDBNO     = 0
      IDBCO     = 0
      IDBALK4   = 0
      IDBACET   = 0
      IDBMEK    = 0
      IDBALD2   = 0
      IDBPRPE   = 0
      IDBC3H8   = 0
      IDBCH2O   = 0
      IDBC2H6   = 0
      IDBSO2    = 0
      IDBNH3    = 0
      IDBBC     = 0
      IDBOC     = 0
      IDBCO2    = 0
      IDBTOLU   = 0
      IDBBENZ   = 0
      IDBXYLE   = 0
      IDBC2H2   = 0
      IDBC2H4   = 0
      IDBHAC    = 0
      IDBGLYC   = 0
      IDBGLYX   = 0
      IDBMGLY   = 0
!kjw
      IDBCH4    = 0
!kjw
      IDBNAP    = 0 ! SOAupdate hotp 7/25/10

      !-----------------------------------
      ! Initialize tagged Hg index arrays
      !-----------------------------------
      IF ( ITS_A_MERCURY_SIM() ) THEN

         ! Initialize category flags
!- eds 8/31/10 -------------------------------------------------------
         ID_Hg_tot = 0
!         ID_Hg_na  = 0
!         ID_Hg_eu  = 0
!         ID_Hg_as  = 0
!         ID_Hg_rw  = 0
!         ID_Hg_oc  = 0
!         ID_Hg_ln  = 0
!         ID_Hg_nt  = 0

         ID_Hg_can = 0
         ID_Hg_usa = 0
         ID_Hg_cam = 0
         ID_Hg_sam = 0
         ID_Hg_waf = 0
         ID_Hg_eaf = 0
         ID_Hg_saf = 0
         ID_Hg_naf = 0
         ID_Hg_eur = 0
         ID_Hg_eeu = 0
         ID_Hg_mde = 0
         ID_Hg_sov = 0
         ID_Hg_sas = 0
         ID_Hg_eas = 0
         ID_Hg_sea = 0
         ID_Hg_jpn = 0
         ID_Hg_oce = 0
         ID_Hg_so  = 0
         ID_Hg_bb  = 0
         ID_Hg_geo = 0
         ID_Hg_atl = 0
         ID_Hg_nat = 0
         ID_Hg_sat = 0
         ID_Hg_npa = 0
         ID_Hg_arc = 0
         ID_Hg_ant = 0
         ID_Hg_ocn = 0
         ID_Hg_str = 0
!------------------------------------------------------------------------


         ! Number of Hg categories
         IF ( LSPLIT ) THEN
!- eds 8/31/10 ----------------------------------------------------------
!            N_Hg_CATS = 8
            N_Hg_CATS = 29 ! New Tagged
!------------------------------------------------------------------------
         ELSE
            N_Hg_CATS = 1
         ENDIF

         ! Index array for Hg0 tracers
         ALLOCATE( ID_Hg0( N_Hg_CATS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ID_Hg0' )
         ID_Hg0 = 0

         ! Index array for Hg2 tracers
         ALLOCATE( ID_Hg2( N_Hg_CATS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ID_Hg2' )
         ID_Hg2 = 0

         ! Index array for HgP tracers
         ALLOCATE( ID_HgP( N_Hg_CATS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ID_HgP' )
         ID_HgP = 0

      ENDIF

      ! Return to calling program
      END SUBROUTINE INIT_TRACERID

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_TRACERID
!
!******************************************************************************
!  Subroutine CLEANUP_TRACERID deallocates all module arrays (bmy, 12/16/05)
!
!  NOTES:
!******************************************************************************
!      
      !=================================================================
      ! CLEANUP_TRACERID begins here!
      !=================================================================
      IF ( ALLOCATED( ID_Hg0 ) ) DEALLOCATE( ID_Hg0 )
      IF ( ALLOCATED( ID_Hg2 ) ) DEALLOCATE( ID_Hg2 )
      IF ( ALLOCATED( ID_HgP ) ) DEALLOCATE( ID_HgP )
      
      ! Return to calling program
      END SUBROUTINE CLEANUP_TRACERID

!------------------------------------------------------------------------------

      ! End of module
      END MODULE TRACERID_MOD
