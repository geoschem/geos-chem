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
!  (50 ) IDTOX     (INTEGER) : Ox    index w/in STT array ("CMN")     
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
!  (109) IDENOX    (INTEGER) : NOx   index w/in EMISRRN array ("F77_CMN_O3")  
!  (110) IDEOX     (INTEGER) : Ox    index w/in EMISRR  array ("F77_CMN_O3")  
!  (111) IDECO     (INTEGER) : CO    index w/in EMISRR  array ("F77_CMN_O3")     
!  (112) IDEPRPE   (INTEGER) : PRPE  index w/in EMISRR  array ("F77_CMN_O3")     
!  (113) IDEC3H8   (INTEGER) : C3H8  index w/in EMISRR  array ("F77_CMN_O3")     
!  (114) IDEALK4   (INTEGER) : ALK4  index w/in EMISRR  array ("F77_CMN_O3")     
!  (115) IDEC2H6   (INTEGER) : C2H6  index w/in EMISRR  array ("F77_CMN_O3")    
!  (116) IDEISOP   (INTEGER) : ISOP  index w/in EMISRR  array ("F77_CMN_O3")    
!  (117) IDEACET   (INTEGER) : ACET  index w/in EMISRR  array ("F77_CMN_O3")     
!  (118) IDEMEK    (INTEGER) : MEK   index w/in EMISRR  array ("F77_CMN_O3")     
!  (119) IDEALD2   (INTEGER) : ALD2  index w/in EMISRR  array ("F77_CMN_O3")    
!  (120) IDECH2O   (INTEGER) : CH2O  index w/in EMISRR  array ("F77_CMN_O3")     
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
!  (20) Updates for mercury simulation (ccc, 5/17/10)
!  17 Feb 2011 - R. Yantosca - Add modifications for APM microphysics (G. Luo)
!******************************************************************************
!
      IMPLICIT NONE
#     include "define.h"

      ! for CTM tracers
!-- add tracers for new isoprene chemistry. (fp, 01/27/10)
      INTEGER, PARAMETER :: NNNTRID  = 75
!      INTEGER, PARAMETER :: MMMEMBER = 10
!      INTEGER, PARAMETER :: NNNTRID  = 100
      INTEGER, PARAMETER :: MMMEMBER = 15
      INTEGER            :: NMEMBER(NNNTRID) 
      INTEGER            :: IDTRMB(NNNTRID,MMMEMBER)
      INTEGER            :: IDEMIS(NNNTRID)
      REAL*8             :: CTRMB(NNNTRID,MMMEMBER)
     
      ! ID's for SMVGEAR species
      INTEGER            :: IDO3,    IDNO2,    IDNO3,   IDN2O5,  IDHNO4
      INTEGER            :: IDOX,    IDNOX,    IDHC1,   IDNO,    IDHNO2
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
      ! added for aromatic (dkh, 10/06/06)  
      INTEGER            :: IDLBRO2H, IDLBRO2N  
      INTEGER            :: IDLTRO2H, IDLTRO2N  
      INTEGER            :: IDLXRO2H, IDLXRO2N  

      ! GEOS-CHEM tracer ID's
      INTEGER            :: IDTNOX,  IDTOX,    IDTPAN,  IDTCO,   IDTALK4
      INTEGER            :: IDTISOP, IDTHNO3,  IDTH2O2, IDTACET, IDTMEK
      INTEGER            :: IDTALD2, IDTRCHO,  IDTMVK,  IDTMACR, IDTPMN
      INTEGER            :: IDTPPN,  IDTISN2,  IDTR4N2, IDTPRPE, IDTC3H8
      INTEGER            :: IDTCH2O, IDTMP,    IDTN2O5, IDTHNO4, IDTC2H6
      INTEGER            :: IDTDMS,  IDTSO2,   IDTSO4,  IDTMSA,  IDTNH3
      INTEGER            :: IDTNH4,  IDTNIT,   IDTRN,   IDTPB,   IDTBE7
      INTEGER            :: IDTBCPI, IDTBCPO,  IDTOCPI, IDTOCPO, IDTDST1

      INTEGER            :: IDTDST2, IDTDST3,  IDTDST4, IDTSALA, IDTSALC
      INTEGER            :: IDTALPH, IDTLIMO,  IDTALCO, IDTSOG1, IDTSOG2 
      INTEGER            :: IDTSOG3, IDTSOG4,  IDTSOA1, IDTSOA2, IDTSOA3
      INTEGER            :: IDTSOA4, IDTHG0,   IDTHg2,  IDTHgP,  IDTAS
      INTEGER            :: IDTAHS,  IDTLET,   IDTNH4aq,IDTSO4aq,IDTSO4s
      INTEGER            :: IDTNITs
      INTEGER            :: IDTBENZ,  IDTTOLU,   IDTXYLE,  IDTMONX
      INTEGER            :: IDTGLYX, IDTMGLY
      INTEGER            :: IDTSOAG, IDTSOAM
      INTEGER            :: IDTC2H2, IDTC2H4
      INTEGER            :: IDTMBO, IDTGLYC
      INTEGER            :: IDTAPAN, IDTENPAN, IDTMPAN, IDTNIPAN
      INTEGER            :: IDTGLPAN,  IDTGPAN
      INTEGER            :: IDTHAC
      INTEGER            :: IDTISOPN
      INTEGER            :: IDTPROPNN
      INTEGER            :: IDTMOBA,  IDTMMN
      INTEGER            :: IDTRIP
      INTEGER            :: IDTIEPOX
      INTEGER            :: IDTPYPAN,IDTMAP
      INTEGER            :: IDTAP

      INTEGER		 :: IDTHCOOH, IDTACTA  ! (mjc, 4/13/11)

      ! For H2/HD simulation
      INTEGER            :: IDTH2, IDTHD ! (hup, phs, 9/18/07)

      ! update for arom (dkh, 10/05/06)  
      INTEGER            :: IDTSOA5, IDTSOG5


      ! For tagged Hg simulation
      INTEGER              :: N_Hg_CATS
      INTEGER, ALLOCATABLE :: ID_Hg0(:),  ID_Hg2(:), ID_HgP(:)
      INTEGER              :: ID_Hg_tot,  ID_Hg_na,  ID_Hg_eu 
      INTEGER              :: ID_Hg_as,   ID_Hg_rw,  ID_Hg_oc
      INTEGER              :: ID_Hg_ln,   ID_Hg_nt

      ! GEOS-CHEM emission ID's
      INTEGER            :: IDENOX,  IDEOX,    IDECO,   IDEPRPE, IDEC3H8
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
      
      ! GEOS-CHEM biofuel ID's
      INTEGER            :: IDBFNOX,  IDBFCO,   IDBFALK4, IDBFACET 
      INTEGER            :: IDBFMEK,  IDBFALD2, IDBFPRPE, IDBFC3H8
      INTEGER            :: IDBFCH2O, IDBFC2H6
      INTEGER            :: IDBFBENZ,  IDBFTOLU,   IDBFXYLE
      INTEGER            :: IDBFC2H2, IDBFC2H4
      INTEGER            :: IDBFGLYC
      INTEGER            :: IDBFGLYX, IDBFMGLY
      INTEGER            :: IDBFHAC

      !BIOMASS BURNING  !(FP,01/27/10)
      INTEGER            :: IDBNOX,  IDBCO,   IDBALK4, IDBACET 
      INTEGER            :: IDBMEK,  IDBALD2, IDBPRPE, IDBC3H8
      INTEGER            :: IDBCH2O, IDBC2H6
      INTEGER            :: IDBSO2,  IDBNH3
      INTEGER            :: IDBBC,   IDBOC
      INTEGER            :: IDBCO2

      INTEGER            :: IDBTOLU, IDBBENZ, IDBXYLE
      ! add dicarbonyl emissions (ccc, 01/27/10)
      INTEGER            :: IDBGLYX, IDBMGLY, IDBC2H4, IDBC2H2
      INTEGER            :: IDBGLYC, IDBHAC      
      

      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE TRACERID
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
!  (18) Add online definition of IDBs. (fp, hotp, 01/10)
!  (19) Add definition of IDTSOG5, IDTSOA5, IDTISOPN, PROPNN, 
!        IDTAP, IDTMOBA, IDTMMN, IDTRIP, IDTIEPOX, IDTMAP, IDTPYPAN
!        (fp, hotp, 01/10)
!  (20) Change hard-wired IDEs to dynamically defined IDEs. (fp, hotp, 01/10)
!  (21) Add IDEMS definitions for new species (fp, hotp, 01/10)
!  (22) Add writing check on IDs. (hotp, 01/10)
!  17 Aug 2011 - R. Yantosca - Print values of 3 extra emission IDExxx flags
!  22 Aug 2011 - R. Yantosca - Add modifications for APM microphysics (G. Luo)
!  16 Feb 2012 - M. Payer    - Account for alternate spellings of the Rn and Pb
!                              tracers in input.geos
!******************************************************************************
!
      ! References to F90 modules
      USE CHARPAK_MOD, ONLY : TRANUC
      USE LOGICAL_MOD, ONLY : LSPLIT
      USE TRACER_MOD,  ONLY : ITS_A_C2H6_SIM,  ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,  ONLY : ITS_A_TAGCO_SIM, ITS_A_MERCURY_SIM
      USE TRACER_MOD,  ONLY : ITS_A_H2HD_SIM
      USE TRACER_MOD,  ONLY : N_TRACERS,       TRACER_NAME
      ! Error check for array dimensions (hotp 8/1/09)
      USE ERROR_MOD,   ONLY : ERROR_STOP

      USE CMN_SIZE_MOD      ! Size parameters
      USE COMODE_LOOP_MOD   ! IDEMS

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
      ! NNNTRID should be the same as NNPAR in F77_CMN_SIZE
      IF ( NNNTRID .NE. NNPAR ) THEN
#if   defined( APM )
         write(*,*)'NNNTRID .NE. NNPAR!',NNNTRID,NNPAR
#else
         CALL ERROR_STOP ( 'NNNTRID .NE. NNPAR!',
     &                                     'tracerid_mod,F77_CMN_SIZE' )
#endif
         ! If you get this error, an appropriate fix is to increase 
         ! either NNNTRID (in tracerid_mod) or NNPAR (in F77_CMN_SIZE) so
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

         ! Convert tracer name to upper case.  TCNAME is in the "CMN" header
         ! file -- we might use something better later on (bmy, 11/12/02)
         NAME = TRACER_NAME(N)
         CALL TRANUC( NAME )

         ! Find each tracer
         SELECT CASE ( TRIM( NAME ) )
        
            !------------------------
            ! Full chem tracers
            !------------------------
            CASE ( 'NOX' )
               COUNT    = COUNT + 1
               IDTNOX   = N
               IDBFNOX  = COUNT

               !(fp)              
               COUNT_BB = COUNT_BB + 1
               IDBNOX   = COUNT_BB

            CASE ( 'OX' )
               IDTOX    = N

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

            CASE ( 'ALPH' )
               IDTALPH  = N

            CASE ( 'LIMO' )
               IDTLIMO  = N

            CASE ( 'ALCO' )
               IDTALCO  = N

            CASE ( 'SOG1' )
               IDTSOG1  = N

            CASE ( 'SOG2' )
               IDTSOG2  = N

            CASE ( 'SOG3' )
               IDTSOG3  = N

            CASE ( 'SOG4' )
               IDTSOG4  = N

            !(hotp 5/25/09)
            CASE ( 'SOG5' )
               IDTSOG5  = N

            CASE ( 'SOA1' )
               IDTSOA1  = N

            CASE ( 'SOA2' )
               IDTSOA2  = N

            CASE ( 'SOA3' )
               IDTSOA3  = N

            CASE ( 'SOA4' )
               IDTSOA4  = N

            !(hotp 5/25/09)
            CASE ( 'SOA5' )
               IDTSOA5  = N

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

            CASE ( 'HG0_AN_NA', 'HG0_AN' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_na          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_AN_EU', 'HG0_AE' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_eu          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_AN_AS', 'HG0_AA' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_as          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_AN_RW', 'HG0_AR' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_rw          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_OC' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_oc          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_LN' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_ln          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG0_NT' )
               COUNT_Hg0         = COUNT_Hg0 + 1
               ID_Hg_nt          = COUNT_Hg0
               ID_Hg0(COUNT_Hg0) = N

            CASE ( 'HG2_AN_NA', 'HG2_AN' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_AN_EU', 'HG2_AE' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_AN_AS', 'HG2_AA' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_AN_RW', 'HG2_AR' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_OC' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_LN' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HG2_NT' )
               COUNT_Hg2         = COUNT_Hg2 + 1
               ID_Hg2(COUNT_Hg2) = N

            CASE ( 'HGP_AN_NA', 'HGP_AN' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_AN_EU', 'HGP_AE' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_AN_AS', 'HGP_AA' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_AN_RW', 'HGP_AR' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_OC' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_LN' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

            CASE ( 'HGP_NT' )
               COUNT_HgP         = COUNT_HgP + 1
               ID_HgP(COUNT_HgP) = N

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
         IF ( IDTNOX > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDENOX    = NEMANTHRO
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
         IF ( IDTOX > 0 ) THEN
            NEMANTHRO = NEMANTHRO + 1
            IDEOX     = NEMANTHRO
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
      ENDIF
      
      !=================================================================
      ! Fill IDEMS with appropriate tracer ID #'s
      !
      ! NOTE: IDEMS is in "comode.h", maybe later split this off into
      ! an F90 module somehow.  Think about this later. (bmy, 11/12/02)
      !=================================================================
      IF ( IDENOX  /= 0 ) IDEMS(IDENOX ) = IDTNOX
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
      IF ( IDEOX   /= 0 ) IDEMS(IDEOX  ) = IDTOX    ! PHS
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
      IF ( IDEHAC /= 0 )  IDEMS(IDEHAC ) = IDTHAC

      ! Echo anthro & biogenic emitted tracers
      WRITE( 6, 100 ) IDEMS ( 1:NEMANTHRO+NEMBIOG )
 100  FORMAT( /, 'TRACERID: Emitted tracers (anthro & bio) :', 20i3 )

      ! print additional information (hotp 8/4/09)
       WRITE( 6, 115 ) 'TRACERID: CHECK FULLCHEM IDEs'
       WRITE( 6, 115 ) '-----------------------------------------------'
       WRITE( 6, 115 ) 'NEMANTHRO: 10 for standard sims b4 v811'
       WRITE( 6, 115 ) '          + 2 for HNO3 and Ox'
       WRITE( 6, 115 ) '          + 5 for aromatics (B,T,X, C2H2, C2H4)'
       WRITE( 6, 115 ) 'NEMANTHRO: 17 max'
       WRITE( 6, 120 ) 'THIS NEMANTHRO: ', NEMANTHRO
       WRITE( 6, 115 ) '-----------------------------------------------'
       WRITE( 6, 115 ) 'NEMBIOG:    1 for isoprene'
       WRITE( 6, 115 ) '          + 6 for dicarbonyls'
       WRITE( 6, 115 ) 'NEMBIOG:    7 max'
       WRITE( 6, 120 ) 'THIS NEMBIOG:   ', NEMBIOG
       WRITE( 6, 115 ) '-----------------------------------------'
       WRITE( 6, 120 ) 'IDENOX  ', IDENOX
       WRITE( 6, 120 ) 'IDECO   ', IDECO
       WRITE( 6, 120 ) 'IDEPRPE ', IDEPRPE
       WRITE( 6, 120 ) 'IDEC3H8 ', IDEC3H8
       WRITE( 6, 120 ) 'IDEALK4 ', IDEALK4
       WRITE( 6, 120 ) 'IDEC2H6 ', IDEC2H6
       WRITE( 6, 120 ) 'IDEACET ', IDEACET
       WRITE( 6, 120 ) 'IDEMEK  ', IDEMEK
       WRITE( 6, 120 ) 'IDEALD2 ', IDEALD2
       WRITE( 6, 120 ) 'IDECH2O ', IDECH2O
       WRITE( 6, 120 ) 'IDEOX   ', IDEOX
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
       WRITE( 6, 115 ) '-----------------------------------------'
 115   FORMAT( a )
 120   FORMAT( a, i3 )

      ! Return to calling program
      END SUBROUTINE TRACERID

!------------------------------------------------------------------------------

      SUBROUTINE SETTRACE
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
!  (7 ) Comment out IDBENZ, IDTOLU, IDXYLE, IDMONX, IDGLYX, IDMGLY 
!        (ccc, 01/29/10)
!  (8 ) Add IDs for Henze's aromatics (hotp, 01/10)
!  (9 ) Add debug print. (hotp, 01/10)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,  ONLY : ERROR_STOP
      USE TRACER_MOD, ONLY : ID_EMITTED,     N_TRACERS
      USE TRACER_MOD, ONLY : TRACER_COEFF,   TRACER_CONST
      USE TRACER_MOD, ONLY : TRACER_N_CONST, TRACER_NAME
      ! for Debug (hotp 8/5/09)
      USE LOGICAL_MOD, ONLY : LPRT

      USE CMN_SIZE_MOD      ! Size parameters
      USE COMODE_LOOP_MOD   ! NAMEGAS

      ! Local variabales
      INTEGER             :: I, J, T, C
 
      !=================================================================
      ! SETTRACE begins here!
      !=================================================================

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
                  IF ( NAMEGAS(J) == TRACER_NAME(T) ) THEN
                     IDTRMB(T,C) = J 
                  ENDIF

               ELSE 

                  ! Test SMVGEAR species name TRACER_CONST
                  IF ( NAMEGAS(J) == TRACER_CONST(T,C) ) THEN
                     IDTRMB(T,C) = J
                  ENDIF

               ENDIF
            ENDDO

            ! Print when Debug output selected (hotp 8/5/09)
            IF ( LPRT ) THEN
               PRINT*, '###--------------------'
               PRINT*, '### T, C       : ', T, C
               PRINT*, '### NAME       : ', TRACER_NAME(T)
               PRINT*, '### NMEMBER    : ', NMEMBER(T)
               PRINT*, '### CONST(T,C) : ', TRACER_CONST(T,C)
               PRINT*, '### CTRMB(T,C) : ', CTRMB(T,C)
               PRINT*, '### IDEMIS(T)  : ', IDEMIS(T)
               PRINT*, '### IDTRMB(T,C): ', IDTRMB(T,C)
            ENDIF

         ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE SETTRACE
    
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
      ! GET_Hg2_CAT begins here!
      !=================================================================

      ! Pick the Hg2 category number from the tracer number
      IF ( N == ID_Hg0(ID_Hg_tot) ) THEN 

         ! Total
         NN = ID_Hg_tot

      ELSE IF ( N == ID_Hg0(ID_Hg_na) ) THEN 

         ! Anthro North America
         NN = ID_Hg_na

      ELSE IF ( N == ID_Hg0(ID_Hg_eu) ) THEN 

         ! Anthro Europe
         NN = ID_Hg_eu

      ELSE IF ( N == ID_Hg0(ID_Hg_as) ) THEN 

         ! Anthro Asia
         NN = ID_Hg_as

      ELSE IF ( N == ID_Hg0(ID_Hg_rw) ) THEN 

         ! Anthro Rest of World
         NN = ID_Hg_rw

      ELSE IF ( N == ID_Hg0(ID_Hg_oc) ) THEN 

         ! Oceans
         NN = ID_Hg_oc

      ELSE IF ( N == ID_Hg0(ID_Hg_ln) ) THEN 

         ! Land re-emission
         NN = ID_Hg_ln

      ELSE IF ( N == ID_Hg0(ID_Hg_nt) ) THEN 

         ! Natural source
         NN = ID_Hg_nt

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

      ELSE IF ( N == ID_Hg2(ID_Hg_na) ) THEN 

         ! Anthro North America
         NN = ID_Hg_na

      ELSE IF ( N == ID_Hg2(ID_Hg_eu) ) THEN 

         ! Anthro Europe
         NN = ID_Hg_eu

      ELSE IF ( N == ID_Hg2(ID_Hg_as) ) THEN 

         ! Anthro Asia
         NN = ID_Hg_as

      ELSE IF ( N == ID_Hg2(ID_Hg_rw) ) THEN 

         ! Anthro Rest of World
         NN = ID_Hg_rw

      ELSE IF ( N == ID_Hg2(ID_Hg_oc) ) THEN 

         ! Oceans
         NN = ID_Hg_oc

      ELSE IF ( N == ID_Hg2(ID_Hg_ln) ) THEN 

         ! Land re-emission
         NN = ID_Hg_ln

      ELSE IF ( N == ID_Hg2(ID_Hg_nt) ) THEN 

         ! Natural source
         NN = ID_Hg_nt

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

      ELSE IF ( N == ID_HgP(ID_Hg_na) ) THEN 

         ! Anthro North America
         NN = ID_Hg_na

      ELSE IF ( N == ID_HgP(ID_Hg_eu) ) THEN 

         ! Anthro Europe
         NN = ID_Hg_eu

      ELSE IF ( N == ID_HgP(ID_Hg_as) ) THEN 

         ! Anthro Asia
         NN = ID_Hg_as

      ELSE IF ( N == ID_HgP(ID_Hg_rw) ) THEN 

         ! Anthro Rest of World
         NN = ID_Hg_rw

      ELSE IF ( N == ID_HgP(ID_Hg_oc) ) THEN 

         ! Oceans
         NN = ID_Hg_oc

      ELSE IF ( N == ID_HgP(ID_Hg_ln) ) THEN 

         ! Land re-emission
         NN = ID_Hg_ln

      ELSE IF ( N == ID_HgP(ID_Hg_nt) ) THEN 

         ! Natural source
         NN = ID_Hg_nt

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
!  (11) Comment out IDBENZ, IDTOLU, IDXYLE, IDMONX, IDGLYX, IDMGLY
!        (ccc, 01/29/10)
!  (12) Now zero IDs for Henze's aromatics (dkh, 01/10)
!  (13) Now zero IDTs for new tracers in isoprene chemistry, Henze's aromatics,
!       (fp, dkh, hotp, 01/10)
!  (14) Now zero IDBs for all biomass burning emissions. (fp, hotp, 01/10)
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
      IDOX      = 0
      IDNOX     = 0
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

      ! GEOS-CHEM Tracer ID #'s
      IDTNOX    = 0
      IDTOX     = 0  
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


      IDTALPH   = 0
      IDTLIMO   = 0
      IDTALCO   = 0
      IDTSOG1   = 0
      IDTSOG2   = 0
      IDTSOG3   = 0
      IDTSOG4   = 0
      IDTSOA1   = 0
      IDTSOA2   = 0
      IDTSOA3   = 0
      IDTSOA4   = 0
      ! (hotp 5/25/09)
      IDTSOG5   = 0
      IDTSOA5   = 0
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
 
      ! GEOS-CHEM Emission ID #'s
      NEMANTHRO = 0
      NEMBIOG   = 0
      IDENOX    = 0
      IDEOX     = 0
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
      
      ! GEOS-CHEM Biofuel ID #'s
      IDBFNOX   = 0
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

      ! Initialize IDBs for BIOMASS BURNING (hotp 8/24/09)
      IDBNOX    = 0
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

      !-----------------------------------
      ! Initialize tagged Hg index arrays
      !-----------------------------------
      IF ( ITS_A_MERCURY_SIM() ) THEN

         ! Initialize category flags
         ID_Hg_tot = 0
         ID_Hg_na  = 0
         ID_Hg_eu  = 0
         ID_Hg_as  = 0
         ID_Hg_rw  = 0
         ID_Hg_oc  = 0
         ID_Hg_ln  = 0
         ID_Hg_nt  = 0

         ! Number of Hg categories
         IF ( LSPLIT ) THEN
            N_Hg_CATS = 8
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
