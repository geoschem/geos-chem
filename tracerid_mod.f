! $Id: tracerid_mod.f,v 1.10 2005/02/10 19:53:28 bmy Exp $
      MODULE TRACERID_MOD
!
!******************************************************************************
!  Module TRACERID_MOD contains variables which point to SMVGEAR species,
!  CTM Tracers, Biomass species, and biofuel species located within various
!  GEOS-CHEM arrays. (bmy, 11/12/02, 12/20/04)
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
!  (49 ) IDTNOX    (INTEGER) : NOx  index w/in STT array ("tracer_mod.f")
!  (50 ) IDTOX     (INTEGER) : Ox   index w/in STT array ("CMN")     
!  (51 ) IDTPAN    (INTEGER) : PAN  index w/in STT array ("tracer_mod.f")      
!  (52 ) IDTCO     (INTEGER) : CO   index w/in STT array ("tracer_mod.f")     
!  (53 ) IDTALK4   (INTEGER) : ALK4 index w/in STT array ("tracer_mod.f")    
!  (54 ) IDTISOP   (INTEGER) : ISOP index w/in STT array ("tracer_mod.f")    
!  (55 ) IDTHNO3   (INTEGER) : HNO3 index w/in STT array ("tracer_mod.f")     
!  (56 ) IDTH2O2   (INTEGER) : H2O2 index w/in STT array ("tracer_mod.f")    
!  (57 ) IDTACET   (INTEGER) : ACET index w/in STT array ("tracer_mod.f")    
!  (58 ) IDTMEK    (INTEGER) : MEK  index w/in STT array ("tracer_mod.f")    
!  (59 ) IDTALD2   (INTEGER) : ALD2 index w/in STT array ("tracer_mod.f")    
!  (60 ) IDTRCHO   (INTEGER) : RCHO index w/in STT array ("tracer_mod.f")     
!  (61 ) IDTMVK    (INTEGER) : MVK  index w/in STT array ("tracer_mod.f")    
!  (62 ) IDTMACR   (INTEGER) : MACR index w/in STT array ("tracer_mod.f")    
!  (63 ) IDTPMN    (INTEGER) : PMN  index w/in STT array ("tracer_mod.f")    
!  (64 ) IDTPPN    (INTEGER) : PPN  index w/in STT array ("tracer_mod.f")    
!  (65 ) IDTISN2   (INTEGER) : ISN2 index w/in STT array ("tracer_mod.f")     
!  (66 ) IDTR4N2   (INTEGER) : R4N2 index w/in STT array ("tracer_mod.f")    
!  (67 ) IDTPRPE   (INTEGER) : PRPE index w/in STT array ("tracer_mod.f")    
!  (68 ) IDTC3H8   (INTEGER) : C3H8 index w/in STT array ("tracer_mod.f")    
!  (69 ) IDTCH2O   (INTEGER) : CH2O index w/in STT array ("tracer_mod.f")    
!  (70 ) IDTMP     (INTEGER) : MP   index w/in STT array ("tracer_mod.f")    
!  (71 ) IDTN2O5   (INTEGER) : N2O5 index w/in STT array ("tracer_mod.f")    
!  (72 ) IDTHNO4   (INTEGER) : HNO4 index w/in STT array ("tracer_mod.f")    
!  (73 ) IDTC2H6   (INTEGER) : C2H6 index w/in STT array ("tracer_mod.f")    
!  (74 ) IDTDMS    (INTEGER) : DMS  index w/in STT array ("tracer_mod.f")     
!  (75 ) IDTSO2    (INTEGER) : SO2  index w/in STT array ("tracer_mod.f")      
!  (76 ) IDTSO4    (INTEGER) : SO4  index w/in STT array ("tracer_mod.f")    
!  (77 ) IDTMSA    (INTEGER) : MSA  index w/in STT array ("tracer_mod.f")    
!  (78 ) IDTNH3    (INTEGER) : NH3  index w/in STT array ("tracer_mod.f")    
!  (79 ) IDTNH4    (INTEGER) : NH4  index w/in STT array ("tracer_mod.f")    
!  (80 ) IDTNIT    (INTEGER) : NIT  index w/in STT array ("tracer_mod.f")    
!  (81 ) IDTRN     (INTEGER) : Rn   index w/in STT array ("tracer_mod.f")    
!  (82 ) IDTPB     (INTEGER) : Pb   index w/in STT array ("tracer_mod.f")    
!  (83 ) IDTBE7    (INTEGER) : Be7  index w/in STT array ("tracer_mod.f")    
!  (84 ) IDTBCPI   (INTEGER) : BCPI index w/in STT array ("tracer_mod.f")
!  (85 ) IDTBCPO   (INTEGER) : BCPO index w/in STT array ("tracer_mod.f")
!  (86 ) IDTOCPI   (INTEGER) : OCPI index w/in STT array ("tracer_mod.f")
!  (87 ) IDTOCPO   (INTEGER) : OCPO index w/in STT array ("tracer_mod.f")
!  (88 ) IDTALPH   (INTEGER) : ALPH index w/in STT array ("tracer_mod.f")
!  (89 ) IDTLIMO   (INTEGER) : LIMO index w/in STT array ("tracer_mod.f")
!  (90 ) IDTALCO   (INTEGER) : ALCO index w/in STT array ("tracer_mod.f")
!  (91 ) IDTSOG1   (INTEGER) : SOG1 index w/in STT array ("tracer_mod.f")
!  (92 ) IDTSOG2   (INTEGER) : SOG2 index w/in STT array ("tracer_mod.f")
!  (93 ) IDTSOG3   (INTEGER) : SOG3 index w/in STT array ("tracer_mod.f")
!  (94 ) IDTSOA1   (INTEGER) : SOA1 index w/in STT array ("tracer_mod.f")
!  (95 ) IDTSOA2   (INTEGER) : SOA2 index w/in STT array ("tracer_mod.f")
!  (96 ) IDTSOA3   (INTEGER) : SOA3 index w/in STT array ("tracer_mod.f")
!  (97 ) IDTDST1   (INTEGER) : DST1 index w/in STT array ("tracer_mod.f")
!  (98 ) IDTDST2   (INTEGER) : DST2 index w/in STT array ("tracer_mod.f")
!  (99 ) IDTDST3   (INTEGER) : DST3 index w/in STT array ("tracer_mod.f")
!  (100) IDTDST4   (INTEGER) : DST4 index w/in STT array ("tracer_mod.f")
!  (101) IDTSALA   (INTEGER) : SALA index w/in STT array ("tracer_mod.f")
!  (102) IDTSALC   (INTEGER) : SALC index w/in STT array ("tracer_mod.f")
!  (103) IDTHG0    (INTEGER) : Hg0  index w/in STT array ("tracer_mod.f")
!  (104) IDTHG2    (INTEGER) : Hg2  index w/in STT array ("tracer_mod.f")
!  (105) IDTHGP    (INTEGER) : HgP  index w/in STT array ("tracer_mod.f")
!  (106) IDENOX    (INTEGER) : NOx  index w/in EMISRRN array ("CMN_O3")  
!  (107) IDEOX     (INTEGER) : Ox   index w/in EMISRR  array ("CMN_O3")  
!  (108) IDECO     (INTEGER) : CO   index w/in EMISRR  array ("CMN_O3")     
!  (109) IDEPRPE   (INTEGER) : PRPE index w/in EMISRR  array ("CMN_O3")     
!  (110) IDEC3H8   (INTEGER) : C3H8 index w/in EMISRR  array ("CMN_O3")     
!  (111) IDEALK4   (INTEGER) : ALK4 index w/in EMISRR  array ("CMN_O3")     
!  (112) IDEC2H6   (INTEGER) : C2H6 index w/in EMISRR  array ("CMN_O3")    
!  (113) IDEISOP   (INTEGER) : ISOP index w/in EMISRR  array ("CMN_O3")    
!  (114) IDEACET   (INTEGER) : ACET index w/in EMISRR  array ("CMN_O3")     
!  (115) IDEMEK    (INTEGER) : MEK  index w/in EMISRR  array ("CMN_O3")     
!  (116) IDEALD2   (INTEGER) : ALD2 index w/in EMISRR  array ("CMN_O3")    
!  (117) IDECH2O   (INTEGER) : CH2O  index w/in EMISRR  array ("CMN_O3")     
!  (118) NEMBIOG   (INTEGER) : # of biogenic emission species for SMVGEAR
!  (119) NEMANTHRO (INTEGER) : # of anthro   emission species for SMVGEAR
!  (120) IDBNOX    (INTEGER) : NOx  index w/in BURNEMIS array (biomass_mod.f)
!  (121) IDBCO     (INTEGER) : CO   index w/in BURNEMIS array (biomass_mod.f)
!  (122) IDBC2H6   (INTEGER) : C2H6 index w/in BURNEMIS array (biomass_mod.f)
!  (123) IDBPRPE   (INTEGER) : PRPE index w/in BURNEMIS array (biomass_mod.f)
!  (124) IDBACET   (INTEGER) : ACET index w/in BURNEMIS array (biomass_mod.f)
!  (125) IDBMEK    (INTEGER) : MEK  index w/in BURNEMIS array (biomass_mod.f)
!  (126) IDBC3H8   (INTEGER) : C3H8 index w/in BURNEMIS array (biomass_mod.f)
!  (127) IDBALD2   (INTEGER) : ALD2 index w/in BURNEMIS array (biomass_mod.f)
!  (128) IDBCH2O   (INTEGER) : CH2O index w/in BURNEMIS array (biomass_mod.f)
!  (129) IDBFPRPE  (INTEGER) : PRPE index w/in BURNEMIS array (biofuel_mod.f)
!  (130) IDBALK4   (INTEGER) : ALD4 index w/in BURNEMIS array (biomass_mod.f)
!  (131) IDBFNOX   (INTEGER) : NOx  index w/in BIOFUEL array (biofuel_mod.f)
!  (132) IDBFCO    (INTEGER) : CO   index w/in BIOFUEL array (biofuel_mod.f)
!  (133) IDBFALK4  (INTEGER) : ALK4 index w/in BIOFUEL array (biofuel_mod.f)
!  (134) IDBFACET  (INTEGER) : ACET index w/in BIOFUEL array (biofuel_mod.f)
!  (135) IDBFMEK   (INTEGER) : MEK  index w/in BIOFUEL array (biofuel_mod.f)
!  (136) IDBFALD2  (INTEGER) : ALD2 index w/in BIOFUEL array (biofuel_mod.f)
!  (137) IDBFPRPE  (INTEGER) : PRPE index w/in BIOFUEL array (biofuel_mod.f)
!  (138) IDBFC3H8  (INTEGER) : NOx  index w/in BIOFUEL array (biofuel_mod.f)
!  (139) IDBFCH2O  (INTEGER) : NOx  index w/in BIOFUEL array (biofuel_mod.f)
!  (140) IDBFC2H6  (INTEGER) : NOx  index w/in BIOFUEL array (biofuel_mod.f)   
!  Plus 21 flags for tagged Hg tracers
!
!  Module Routines:
!  ============================================================================
!  (1 ) TRACERID      : Defines tracer, biomass, biofuel, & anthro ID numbers
!  (2 ) SETTRACE      : Defines ID numbers for species in SMVGEAR mechanism
!  (3 ) INIT_TRACERID : Zeroes all module variables
!
!  GEOS-CHEM modules referenced by biomass_mod.f
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
!******************************************************************************
!
      IMPLICIT NONE

      ! for CTM tracers
      INTEGER, PARAMETER :: NNNTRID  = 59
      INTEGER, PARAMETER :: MMMEMBER = 10
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
      INTEGER            :: IDDRYO3, IDDRYPAN, IDDRYNO2         

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
      INTEGER            :: IDTSOG3, IDTSOA1,  IDTSOA2, IDTSOA3, IDTHG0
      INTEGER            :: IDTHg2,  IDTHgP,   IDTAS,   IDTAHS,  IDTLET
      INTEGER            :: IDTNH4aq,IDTSO4aq

      ! Extra tracer ID's for tagged mercury simulation
      ! (we need these or else the wetdep code won't work)
      INTEGER            :: IDTHg0_NA, IDTHg0_EU, IDTHg0_AS, IDTHg0_RW
      INTEGER            :: IDTHg0_OC, IDTHg0_LN, IDTHg0_NT
      INTEGER            :: IDTHg2_NA, IDTHg2_EU, IDTHg2_AS, IDTHg2_RW
      INTEGER            :: IDTHg2_OC, IDTHg2_LN, IDTHg2_NT
      INTEGER            :: IDTHgP_NA, IDTHgP_EU, IDTHgP_AS, IDTHgP_RW
      INTEGER            :: IDTHgP_OC, IDTHgP_LN, IDTHgP_NT

      ! GEOS-CHEM emission ID's
      INTEGER            :: IDENOX,  IDEOX,    IDECO,   IDEPRPE, IDEC3H8
      INTEGER            :: IDEALK4, IDEC2H6,  IDEISOP, IDEACET, IDEMEK
      INTEGER            :: IDEALD2, IDECH2O,  NEMBIOG, NEMANTHRO

      ! GEOS-CHEM biomass burning ID's
      INTEGER            :: IDBNOX,  IDBCO,    IDBC2H6, IDBPRPE, IDBACET
      INTEGER            :: IDBMEK,  IDBC3H8,  IDBALD2, IDBCH2O, IDBALK4 

      ! GEOS-CHEM biofuel ID's
      INTEGER            :: IDBFNOX,  IDBFCO,   IDBFALK4, IDBFACET 
      INTEGER            :: IDBFMEK,  IDBFALD2, IDBFPRPE, IDBFC3H8
      INTEGER            :: IDBFCH2O, IDBFC2H6
      
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE TRACERID
!
!******************************************************************************
!  Subroutine TRACERID reads the "tracer.dat" file and determines which
!  tracers, emission species, biomass burning species, and biofuel burning
!  species are turned on/off. (bmy, 3/16/01, 1/26/05)
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
!******************************************************************************
!
      ! References to F90 modules
      USE CHARPAK_MOD, ONLY : TRANUC
      USE TRACER_MOD

#     include "CMN_SIZE"  ! Size parameters
#     include "comode.h"  ! IDEMS

      ! Local variables
      INTEGER           :: N, COUNT
      CHARACTER(LEN=14) :: NAME
      
      !=================================================================
      ! TRACERID begins here!
      !
      ! NOTE: There are still some vestiges of historical baggage, we
      !       will get rid of this as time allows (bmy, 11/12/02)
      !=================================================================

      ! Zero all ID #'s
      CALL INIT_TRACERID
      
      ! Initialize counter
      COUNT = 0

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
               IDBNOX   = COUNT
               IDBFNOX  = COUNT

            CASE ( 'OX' )
               IDTOX    = N

            CASE ( 'PAN' )
               IDTPAN   = N
               
            CASE ( 'CO' )
               COUNT    = COUNT + 1
               IDTCO    = N
               IDBCO    = COUNT
               IDBFCO   = COUNT

               ! Special case: Tagged CO
               ! Set some emission flags and then exit
               IF ( ITS_A_TAGCO_SIM() ) THEN 
                  NEMANTHRO = 1
                  IDECO     = 2
                  IDTISOP   = 1
                  EXIT
               ENDIF

            CASE ( 'ALK4' )
               COUNT    = COUNT + 1
               IDTALK4  = N
               IDBALK4  = COUNT
               IDBFALK4 = COUNT

            CASE ( 'ISOP' )
               IDTISOP  = N
               
            CASE ( 'HNO3' )
               IDTHNO3  = N

            CASE ( 'H2O2' )
               IDTH2O2  = N

            CASE ( 'ACET' )
               COUNT    = COUNT + 1
               IDTACET  = N
               IDBACET  = COUNT
               IDBFACET = COUNT

            CASE ( 'MEK' )
               COUNT    = COUNT + 1
               IDTMEK   = N
               IDBMEK   = COUNT
               IDBFMEK  = COUNT

            CASE ( 'ALD2' )
               COUNT    = COUNT + 1
               IDTALD2  = N
               IDBALD2  = COUNT
               IDBFALD2 = COUNT

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
               IDBPRPE  = COUNT
               IDBFPRPE = COUNT

            CASE ( 'C3H8' )
               COUNT    = COUNT + 1
               IDTC3H8  = N
               IDBC3H8  = COUNT
               IDBFC3H8 = COUNT

            CASE ( 'CH2O' )
               COUNT    = COUNT + 1
               IDTCH2O  = N
               IDBCH2O  = COUNT
               IDBFCH2O = COUNT

            CASE ( 'C2H6' )
               COUNT    = COUNT + 1
               IDTC2H6  = N
               IDBC2H6  = COUNT
               IDBFC2H6 = COUNT

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

            CASE ( 'SO4' )
               IDTSO4   = N

            CASE ( 'MSA' )
               IDTMSA   = N

            CASE ( 'NH3' )
               IDTNH3   = N

            CASE ( 'NH4' )
               IDTNH4   = N

            CASE ( 'NIT' )
               IDTNIT   = N

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

            CASE ( 'OCPI' )
               IDTOCPI  = N

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

            CASE ( 'SOA1' )
               IDTSOA1  = N

            CASE ( 'SOA2' )
               IDTSOA2  = N

            CASE ( 'SOA3' )
               IDTSOA3  = N

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
            ! Rn-Pb-Be tracers
            !--------------------------------
            CASE ( 'RN' )
               IDTRN    = N
 
            CASE ( 'PB' )
               IDTPB    = N

            CASE ( 'BE7' )
               IDTBE7   = N

            !--------------------------------
            ! CH3I and HCN tracers
            !--------------------------------

            ! Special case: CH3I needs CO biomass/biofuel
            CASE ( 'CH3I' )
               COUNT    = COUNT + 1
               IDTCO    = 1
               IDBCO    = COUNT
               IDBFCO   = COUNT
               EXIT

            ! Special case: HCN needs CO biomass/biofuel
            CASE ( 'HCN' )
               COUNT    = COUNT + 1
               IDTCO    = 1
               IDBCO    = COUNT
               IDBFCO   = COUNT
               EXIT

            !--------------------------------
            ! Total & tagged mercury tracers 
            ! (eck, bmy, 12/14/04)
            !--------------------------------
            CASE ( 'HG0' )
               IDTHG0   = N
 
            CASE ( 'HG2' )
               IDTHG2   = N

            CASE ( 'HGP' )
               IDTHGP   = N

            CASE ( 'HG0_AN_NA' )
               IDTHG0_NA = N

            CASE ( 'HG0_AN_EU' )
               IDTHG0_EU = N

            CASE ( 'HG0_AN_AS' )
               IDTHG0_AS = N

            CASE ( 'HG0_AN_RW' )
               IDTHG0_RW = N

            CASE ( 'HG0_OC' )
               IDTHG0_OC = N

            CASE ( 'HG0_LN' )
               IDTHG0_LN = N

            CASE ( 'HG0_NT' )
               IDTHG0_NT = N

            CASE ( 'HG2_AN_NA' )
               IDTHG2_NA = N

            CASE ( 'HG2_AN_EU' )
               IDTHG2_EU = N

            CASE ( 'HG2_AN_AS' )
               IDTHG2_AS = N

            CASE ( 'HG2_AN_RW' )
               IDTHG2_RW = N

            CASE ( 'HG2_OC' )
               IDTHG2_OC = N

            CASE ( 'HG2_LN' )
               IDTHG2_LN = N

            CASE ( 'HG2_NT' )
               IDTHG2_NT = N

            CASE ( 'HGP_AN_NA' )
               IDTHGP_NA = N

            CASE ( 'HGP_AN_EU' )
               IDTHGP_EU = N

            CASE ( 'HGP_AN_AS' )
               IDTHGP_AS = N

            CASE ( 'HGP_AN_RW' )
               IDTHGP_RW = N

            CASE ( 'HGP_OC' )
               IDTHGP_OC = N

            CASE ( 'HGP_LN' )
               IDTHGP_LN = N

            CASE ( 'HGP_NT' )
               IDTHGP_NT = N

            CASE DEFAULT
               ! Nothing

         END SELECT
      ENDDO
      
      !=================================================================
      ! SPECIAL CASE: we need to hardwire the emission flags so that
      ! they are in the same order as the old emissions code.  The 
      ! order should be: 1 4 18 19 5 21 9 10 11 20 6.  Think of a 
      ! better way to implement this later on. (bmy, 12/20/04)
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN
         NEMANTHRO = 10
         NEMBIOG   = 1
         IDENOX    = 1
         IDECO     = 2
         IDEPRPE   = 3
         IDEC3H8   = 4
         IDEALK4   = 5
         IDEC2H6   = 6
         IDEACET   = 7
         IDEMEK    = 8
         IDEALD2   = 9
         IDECH2O   = 10
         IDEISOP   = 11
      ENDIF
      
      !----------------------------------------------------------------------
      ! Prior to 1/25/05:
      !!=================================================================
      !! SPECIAL CASE: For the offline sulfate simulation, we also need
      !! to turn on the biofuel and biomass CO flags (bmy, 11/12/02)
      !!=================================================================
      !IF ( ITS_AN_AEROSOL_SIM() ) THEN
      !   IDTCO     = 1
      !   IDBCO     = 1
      !   IDBFCO    = 1
      !ENDIF
      !----------------------------------------------------------------------

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

      ! Echo anthro & biogenic emitted tracers
      WRITE( 6, 100 ) IDEMS ( 1:NEMANTHRO+NEMBIOG )
 100  FORMAT( /, 'TRACERID: Emitted tracers (anthro & bio) :', 20i3 )

      ! Return to calling program
      END SUBROUTINE TRACERID

!------------------------------------------------------------------------------

      SUBROUTINE SETTRACE
!
!******************************************************************************
!  Subroutine SETTRACE flags certain chemical species w/in the SMVGEAR full
!  chemistry mechanism. (lwh, jyl, gmg, djj, 1990's; bmy, 11/12/02, 7/20/04)
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
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP
      USE TRACER_MOD 

#     include "CMN_SIZE"  ! Size parameters
#     include "comode.h"  ! NAMEGAS

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

            !### Debug
            !PRINT*, '###--------------------'
            !PRINT*, '### T, C       : ', T, C
            !PRINT*, '### NAME       : ', TRACER_NAME(T)
            !PRINT*, '### NMEMBER    : ', NMEMBER(T)
            !PRINT*, '### CONST(T,C) : ', TRACER_CONST(T,C)
            !PRINT*, '### CTRMB(T,C) : ', CTRMB(T,C)
            !PRINT*, '### IDEMIS(T)  : ', IDEMIS(T)
            !PRINT*, '### IDTRMB(T,C): ', IDTRMB(T,C)

         ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE SETTRACE
    
!------------------------------------------------------------------------------

      SUBROUTINE INIT_TRACERID
!
!******************************************************************************
!  Subroutine INIT_TRACERID zeroes module variables. (bmy, 11/12/02, 12/20/04)
!
!  NOTES:
!  (1 ) Now also zero IDDMS, IDSO2, IDSO4, IDMSA (rjp, bmy, 3/23/03)
!  (2 ) Now zero extra flags for carbon & dust tracers (rjp, tdf, bmy, 4/1/04)
!  (3 ) Now zero extra flags for seasalt tracers (rjp, bec, bmy, 4/1/04)
!  (4 ) Now zero extra flags for SOA tracers (rjp, bmy, 7/13/04)
!  (5 ) Now zero IDTHG0, IDTHG2, IDTHGP + tagged Hg's (eck, bmy, 12/7/04)
!  (6 ) Now zero IDTAS, IDTAHS, IDTLET, IDTNH4aq, IDTSO4aq (cas, bmy, 12/20/04)
!******************************************************************************
!
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

      ! GEOS-CHEM Tracer ID #'s
      IDTNOX    = 0
      IDTOX     = 0  
      IDTPAN    = 0
      IDTCO     = 0
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
      IDTSOA1   = 0
      IDTSOA2   = 0
      IDTSOA3   = 0
      IDTDST1   = 0
      IDTDST2   = 0
      IDTDST3   = 0
      IDTDST4   = 0
      IDTSALA   = 0
      IDTSALC   = 0
      IDTRN     = 0
      IDTPB     = 0
      IDTBE7    = 0
      IDTHG0    = 0
      IDTHG2    = 0
      IDTHGP    = 0

      ! GEOS-CHEM Tagged Hg tracers
      IDTHG2_NA = 0 
      IDTHG2_EU = 0 
      IDTHG2_AS = 0 
      IDTHG2_RW = 0
      IDTHG2_OC = 0 
      IDTHG2_LN = 0 
      IDTHG2_NT = 0
      IDTHGP_NA = 0 
      IDTHGP_EU = 0 
      IDTHGP_AS = 0 
      IDTHGP_RW = 0
      IDTHGP_OC = 0 
      IDTHGP_LN = 0 
      IDTHGP_NT = 0
      
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
      
      ! GEOS-CHEM Biomass ID #'s 
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

      ! Return to calling program
      END SUBROUTINE INIT_TRACERID

!------------------------------------------------------------------------------

      END MODULE TRACERID_MOD
