!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_input_opt_mod
!
! !DESCRIPTION: Module GIGC\_INPUT\_OPT\_MOD contains the derived type
!  for GEOS-Chem options and logical switches.
!\\
!\\
! !INTERFACE:
!
MODULE GIGC_Input_Opt_Mod
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Set_GIGC_Input_Opt
  PUBLIC :: Cleanup_GIGC_Input_Opt
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for Input Options 
  !=========================================================================
  TYPE, PUBLIC :: OptInput

     !----------------------------------------
     ! SIZE PARAMETER fields
     !----------------------------------------
     INTEGER                     :: MAX_DIAG
     INTEGER                     :: MAX_TRCS
     INTEGER                     :: MAX_MEMB
     INTEGER                     :: MAX_FAMS
     INTEGER                     :: MAX_DEP


     !----------------------------------------
     ! SIMULATION MENU fields 
     !----------------------------------------
     INTEGER                     :: NYMDb              
     INTEGER                     :: NHMSb              
     INTEGER                     :: NYMDe              
     INTEGER                     :: NHMSe              
     CHARACTER(LEN=255)          :: RUN_DIR            
     CHARACTER(LEN=255)          :: IN_RST_FILE        
     LOGICAL                     :: LSVGLB             
     CHARACTER(LEN=255)          :: OUT_RST_FILE       
     CHARACTER(LEN=255)          :: DATA_DIR           
     CHARACTER(LEN=255)          :: GCAP_DIR           
     CHARACTER(LEN=255)          :: GEOS_4_DIR         
     CHARACTER(LEN=255)          :: GEOS_5_DIR         
     CHARACTER(LEN=255)          :: GEOS_FP_DIR        
     CHARACTER(LEN=255)          :: MERRA_DIR          
     CHARACTER(LEN=255)          :: DATA_DIR_1x1       
     CHARACTER(LEN=255)          :: TEMP_DIR           
     LOGICAL                     :: LUNZIP             
     LOGICAL                     :: LWAIT              
     LOGICAL                     :: LVARTROP           
     INTEGER                     :: NESTED_I0          
     INTEGER                     :: NESTED_J0          

     !----------------------------------------
     ! TRACER MENU fields
     !----------------------------------------
     INTEGER                     :: N_TRACERS          
     INTEGER,            POINTER :: ID_TRACER(:)       
     CHARACTER(LEN=255), POINTER :: TRACER_NAME(:)     
     REAL*8,             POINTER :: TRACER_MW_G(:)     
     REAL*8,             POINTER :: TRACER_MW_KG(:)    
     REAL*8,             POINTER :: TCVV(:)            
     REAL*8,             POINTER :: XNUMOL(:)          
     INTEGER,            POINTER :: TRACER_N_CONST(:)  
     CHARACTER(LEN=255), POINTER :: TRACER_CONST(:,:)  
     REAL*8,             POINTER :: TRACER_COEFF(:,:)  
     INTEGER,            POINTER :: ID_EMITTED(:)  
     LOGICAL                     :: LSPLIT
     LOGICAL                     :: ITS_A_RnPbBe_SIM
     LOGICAL                     :: ITS_A_CH3I_SIM
     LOGICAL                     :: ITS_A_FULLCHEM_SIM
     LOGICAL                     :: ITS_A_HCN_SIM
     LOGICAL                     :: ITS_A_TAGOX_SIM
     LOGICAL                     :: ITS_A_TAGCO_SIM
     LOGICAL                     :: ITS_A_C2H6_SIM
     LOGICAL                     :: ITS_A_CH4_SIM
     LOGICAL                     :: ITS_AN_AEROSOL_SIM
     LOGICAL                     :: ITS_A_MERCURY_SIM
     LOGICAL                     :: ITS_A_CO2_SIM
     LOGICAL                     :: ITS_A_H2HD_SIM
     LOGICAL                     :: ITS_A_POPS_SIM
     LOGICAL                     :: ITS_A_SPECIALTY_SIM
     LOGICAL                     :: ITS_NOT_COPARAM_OR_CH4

     !----------------------------------------
     ! AEROSOL MENU fields
     !----------------------------------------
     LOGICAL                     :: LSULF              
     LOGICAL                     :: LCRYST             
     LOGICAL                     :: LCARB              
     LOGICAL                     :: LSOA               
     LOGICAL                     :: LSVPOA
     REAL*8                      :: NAPEMISS
     REAL*8                      :: POAEMISSSCALE
     LOGICAL                     :: LDUST              
     LOGICAL                     :: LDEAD              
     LOGICAL                     :: LSSALT             
     LOGICAL                     :: LDICARB            
     REAL*8,             POINTER :: SALA_REDGE_um(:)   
     REAL*8,             POINTER :: SALC_REDGE_um(:)   

     !----------------------------------------
     ! EMISSIONS MENU fields
     !----------------------------------------
     LOGICAL                     :: LEMIS
     INTEGER                     :: TS_EMIS
     LOGICAL                     :: LANTHRO
     LOGICAL                     :: FSCALYR
     LOGICAL                     :: LEMEP 
     LOGICAL                     :: LBRAVO
     LOGICAL                     :: LEDGAR
     LOGICAL                     :: LSTREETS
     LOGICAL                     :: LCAC
     LOGICAL                     :: LNEI05
     LOGICAL                     :: LNEI08
     LOGICAL                     :: LRETRO
     LOGICAL                     :: LNEI99
     LOGICAL                     :: LICARTT
     LOGICAL                     :: LVISTAS
     LOGICAL                     :: LBIOFUEL
     LOGICAL                     :: LBIOGENIC
     LOGICAL                     :: LMEGAN
     LOGICAL                     :: LPECCA 
     LOGICAL                     :: LMEGANMONO
     REAL*8                      :: ISOP_SCALING 
     LOGICAL                     :: LBIOMASS 
     LOGICAL                     :: LBBSEA
     LOGICAL                     :: LTOMSAI
     LOGICAL                     :: LGFED2BB
     LOGICAL                     :: L8DAYBB
     LOGICAL                     :: L3HRBB
     LOGICAL                     :: LSYNOPBB 
     LOGICAL                     :: LGFED3BB
     LOGICAL                     :: LDAYBB3
     LOGICAL                     :: L3HRBB3
     LOGICAL                     :: LAEIC
     LOGICAL                     :: LLIGHTNOX
     LOGICAL                     :: LOTDLOC
     LOGICAL                     :: LSOILNOX
     CHARACTER(LEN=255)          :: SOIL_RST_FILE
     LOGICAL                     :: LFERTILIZERNOX
     REAL*8                      :: NOx_SCALING
     LOGICAL                     :: LEDGARSHIP
     LOGICAL                     :: LICOADSSHIP
     LOGICAL                     :: LEMEPSHIP
     LOGICAL                     :: LSHIPSO2
     LOGICAL                     :: LARCSHIP
     LOGICAL                     :: LCOOKE
     LOGICAL                     :: LHIST
     LOGICAL                     :: HISTYR
     LOGICAL                     :: LWARWICK_VSLS
     LOGICAL                     :: LSSABr2
     LOGICAL                     :: LFIX_PBL_BRO
     REAL*8                      :: Br_SCALING
     LOGICAL                     :: LEDGARNOx
     LOGICAL                     :: LEDGARCO
     LOGICAL                     :: LEDGARSOx
     LOGICAL                     :: LRCP
     LOGICAL                     :: LRCPSHIP
     LOGICAL                     :: LRCPAIR

     !----------------------------------------
     ! CO2 MENU fields
     !----------------------------------------
     LOGICAL                     :: LGENFF
     LOGICAL                     :: LANNFF
     LOGICAL                     :: LMONFF
     LOGICAL                     :: LCHEMCO2
     LOGICAL                     :: LSEASBB
     LOGICAL                     :: LBIODAILY
     LOGICAL                     :: LBIODIURNAL
     LOGICAL                     :: LBIONETORIG
     LOGICAL                     :: LBIONETCLIM
     LOGICAL                     :: LOCN1997
     LOGICAL                     :: LOCN2009ANN
     LOGICAL                     :: LOCN2009MON
     LOGICAL                     :: LSHIPEDG
     LOGICAL                     :: LSHIPICO
     LOGICAL                     :: LPLANE
     LOGICAL                     :: LFFBKGRD
     LOGICAL                     :: LBIOSPHTAG
     LOGICAL                     :: LFOSSILTAG
     LOGICAL                     :: LSHIPTAG
     LOGICAL                     :: LPLANETAG

     !----------------------------------------
     ! FUTURE MENU fields
     !----------------------------------------
     LOGICAL                     :: LFUTURE
     INTEGER                     :: FUTURE_YEAR
     CHARACTER(LEN=255)          :: FUTURE_SCEN

     !----------------------------------------
     ! CHEMISTRY MENU fields
     !----------------------------------------
     LOGICAL                     :: LCHEM
     LOGICAL                     :: LSCHEM
     LOGICAL                     :: LLINOZ
     INTEGER                     :: TS_CHEM
     LOGICAL                     :: LSVCSPEC
     CHARACTER(LEN=255)          :: SPEC_RST_FILE
     LOGICAL                     :: LKPP
     REAL*8                      :: GAMMA_HO2

     !----------------------------------------
     ! CHEMISTRY MENU fields
     !----------------------------------------
     LOGICAL                     :: LTRAN
     LOGICAL                     :: LFILL
     LOGICAL                     :: TPCORE_IORD
     LOGICAL                     :: TPCORE_JORD
     LOGICAL                     :: TPCORE_KORD
     INTEGER                     :: TS_DYN

     !----------------------------------------
     ! CONVECTION MENU fields
     !----------------------------------------
     LOGICAL                     :: LCONV
     LOGICAL                     :: LTURB
     LOGICAL                     :: LNLPBL
     INTEGER                     :: TS_CONV

     !----------------------------------------
     ! DEPOSITION MENU fields
     !----------------------------------------
     LOGICAL                     :: LDRYD
     LOGICAL                     :: LWETD
     LOGICAL                     :: USE_OLSON_2001

     !----------------------------------------
     ! GAMAP MENU fields
     !----------------------------------------
     CHARACTER(LEN=255)          :: GAMAP_DIAGINFO
     CHARACTER(LEN=255)          :: GAMAP_TRACERINFO

     !----------------------------------------
     ! OUTPUT MENU fields
     !----------------------------------------
     INTEGER,            POINTER :: NJDAY(:)

     !----------------------------------------
     ! DIAGNOSTIC MENU fields
     !----------------------------------------
     INTEGER                     :: ND01, LD01
     INTEGER                     :: ND02, LD02
     INTEGER                     :: ND03, LD03
     INTEGER                     :: ND04, LD04
     INTEGER                     :: ND05, LD05
     INTEGER                     :: ND06, LD06
     INTEGER                     :: ND07, LD07
     INTEGER                     :: ND08, LD08
     INTEGER                     :: ND09, LD09
     INTEGER                     :: ND10, LD10
     INTEGER                     :: ND11, LD11
     INTEGER                     :: ND12, LD12
     INTEGER                     :: ND13, LD13
     INTEGER                     :: ND14, LD14
     INTEGER                     :: ND15, LD15
     INTEGER                     :: ND16, LD16
     INTEGER                     :: ND17, LD17
     INTEGER                     :: ND18, LD18
     INTEGER                     :: ND19, LD19
     INTEGER                     :: ND20, LD20
     INTEGER                     :: ND21, LD21
     INTEGER                     :: ND22, LD22
     INTEGER                     :: ND23, LD23
     INTEGER                     :: ND24, LD24
     INTEGER                     :: ND25, LD25
     INTEGER                     :: ND26, LD26
     INTEGER                     :: ND27, LD27
     INTEGER                     :: ND28, LD28
     INTEGER                     :: ND29, LD29
     INTEGER                     :: ND30, LD30
     INTEGER                     :: ND31, LD31
     INTEGER                     :: ND32, LD32
     INTEGER                     :: ND33, LD33
     INTEGER                     :: ND34, LD34
     INTEGER                     :: ND35, LD35
     INTEGER                     :: ND36, LD36
     INTEGER                     :: ND37, LD37
     INTEGER                     :: ND38, LD38
     INTEGER                     :: ND39, LD39
     INTEGER                     :: ND40, LD40
     INTEGER                     :: ND41, LD41
     INTEGER                     :: ND42, LD42
     INTEGER                     :: ND43, LD43
     INTEGER                     :: ND44, LD44
     INTEGER                     :: ND45, LD45
     INTEGER                     :: ND46, LD46
     INTEGER                     :: ND47, LD47
     INTEGER                     :: ND48, LD48
     INTEGER                     :: ND49, LD49
     INTEGER                     :: ND50, LD50
     INTEGER                     :: ND51, LD51
     INTEGER                     :: ND52, LD52
     INTEGER                     :: ND53, LD53
     INTEGER                     :: ND54, LD54
     INTEGER                     :: ND55, LD55
     INTEGER                     :: ND56, LD56
     INTEGER                     :: ND57, LD57
     INTEGER                     :: ND58, LD58
     INTEGER                     :: ND59, LD59
     INTEGER                     :: ND60, LD60
     INTEGER                     :: ND61, LD61
     INTEGER                     :: ND62, LD62
     INTEGER                     :: ND63, LD63
     INTEGER                     :: ND64, LD64
     INTEGER                     :: ND66, LD66
     INTEGER                     :: ND67, LD67
     INTEGER                     :: ND68, LD68
     INTEGER                     :: ND69, LD69
     INTEGER                     :: ND70, LD70
     INTEGER                     :: TS_DIAG
     LOGICAL                     :: LPRT
     INTEGER,            POINTER :: TINDEX(:,:)
     INTEGER,            POINTER :: TCOUNT(:) 				  
     INTEGER,            POINTER :: TMAX(:)
     LOGICAL                     :: DO_DIAG_WRITE

     !----------------------------------------
     ! PLANEFLIGHT MENU fields
     !----------------------------------------
     LOGICAL                     :: DO_PF
     CHARACTER(LEN=255)          :: PF_IFILE
     CHARACTER(LEN=255)          :: PF_OFILE

     !----------------------------------------
     ! ND48 MENU fields
     !----------------------------------------
     LOGICAL                     :: DO_ND48
     CHARACTER(LEN=255)          :: ND48_FILE
     INTEGER                     :: ND48_FREQ
     INTEGER                     :: ND48_N_STA
     INTEGER,            POINTER :: ND48_IARR(:)
     INTEGER,            POINTER :: ND48_JARR(:)
     INTEGER,            POINTER :: ND48_LARR(:)
     INTEGER,            POINTER :: ND48_NARR(:)

     !----------------------------------------
     ! ND49 MENU fields
     !----------------------------------------
     LOGICAL                     :: DO_ND49
     CHARACTER(LEN=255)          :: ND49_FILE
     INTEGER,            POINTER :: ND49_TRACERS(:)
     INTEGER                     :: ND49_FREQ
     INTEGER                     :: ND49_IMIN
     INTEGER                     :: ND49_IMAX
     INTEGER                     :: ND49_JMIN
     INTEGER                     :: ND49_JMAX
     INTEGER                     :: ND49_LMIN
     INTEGER                     :: ND49_LMAX

     !----------------------------------------
     ! ND50 MENU fields
     !----------------------------------------
     LOGICAL                     :: DO_ND50
     CHARACTER(LEN=255)          :: ND50_FILE
     LOGICAL                     :: LND50_HDF
     INTEGER,            POINTER :: ND50_TRACERS(:)
     INTEGER                     :: ND50_IMIN
     INTEGER                     :: ND50_IMAX
     INTEGER                     :: ND50_JMIN
     INTEGER                     :: ND50_JMAX
     INTEGER                     :: ND50_LMIN
     INTEGER                     :: ND50_LMAX

     !----------------------------------------
     ! ND51 MENU fields
     !----------------------------------------
     LOGICAL                     :: DO_ND51
     CHARACTER(LEN=255)          :: ND51_FILE
     INTEGER                     :: LND51_HDF
     INTEGER,            POINTER :: ND51_TRACERS(:)
     REAL*8                      :: ND51_HR_WRITE
     REAL*8                      :: ND51_HR1
     REAL*8                      :: ND51_HR2
     INTEGER                     :: ND51_IMIN
     INTEGER                     :: ND51_IMAX
     INTEGER                     :: ND51_JMIN
     INTEGER                     :: ND51_JMAX
     INTEGER                     :: ND51_LMIN
     INTEGER                     :: ND51_LMAX

     !----------------------------------------
     ! ND51b MENU fields
     !----------------------------------------
     LOGICAL                     :: DO_ND51b
     CHARACTER(LEN=255)          :: ND51b_FILE
     INTEGER                     :: LND51b_HDF
     INTEGER,            POINTER :: ND51b_TRACERS(:)
     REAL*8                      :: ND51b_HR_WRITE
     REAL*8                      :: ND51b_HR1
     REAL*8                      :: ND51b_HR2
     INTEGER                     :: ND51b_IMIN
     INTEGER                     :: ND51b_IMAX
     INTEGER                     :: ND51b_JMIN
     INTEGER                     :: ND51b_JMAX
     INTEGER                     :: ND51b_LMIN
     INTEGER                     :: ND51b_LMAX

     !----------------------------------------
     ! ND63 MENU fields
     !----------------------------------------
     LOGICAL                     :: DO_ND63
     CHARACTER(LEN=255)          :: ND63_FILE
     INTEGER,            POINTER :: ND63_TRACERS(:)
     INTEGER                     :: ND63_FREQ
     INTEGER                     :: ND63_IMIN
     INTEGER                     :: ND63_IMAX
     INTEGER                     :: ND63_JMIN
     INTEGER                     :: ND63_JMAX

     !----------------------------------------
     ! PROD LOSS MENU fields
     !----------------------------------------
     LOGICAL                     :: DO_SAVE_PL
     LOGICAL                     :: LFAMILY
     INTEGER                     :: ND65, LD65
     LOGICAL                     :: DO_SAVE_O3
     INTEGER                     :: NFAM
     REAL*8,             POINTER :: FAM_COEF(:,:)
     CHARACTER(LEN=255), POINTER :: FAM_MEMB(:,:)
     CHARACTER(LEN=255), POINTER :: FAM_NAME(:  )
     INTEGER,            POINTER :: FAM_NMEM(:  )
     CHARACTER(LEN=255), POINTER :: FAM_TYPE(:  )

     !----------------------------------------
     ! UNIX CMDS fields
     !----------------------------------------
     CHARACTER(LEN=255)          :: BACKGROUND
     CHARACTER(LEN=255)          :: REDIRECT
     CHARACTER(LEN=255)          :: REMOVE_CMD
     CHARACTER(LEN=255)          :: SEPARATOR
     CHARACTER(LEN=255)          :: WILD_CARD
     CHARACTER(LEN=255)          :: UNZIP_CMD
     CHARACTER(LEN=255)          :: ZIP_SUFFIX

     !----------------------------------------
     ! NESTED GRID MENU fields
     !----------------------------------------
     LOGICAL                     :: LWINDO
     LOGICAL                     :: LWINDO2x25
     LOGICAL                     :: LWINDO_NA
     CHARACTER(LEN=255)          :: TPBC_DIR_NA
     LOGICAL                     :: LWINDO_EU
     CHARACTER(LEN=255)          :: TPBC_DIR_EU
     LOGICAL                     :: LWINDO_CH
     CHARACTER(LEN=255)          :: TPBC_DIR_CH
     INTEGER                     :: LWINDO_SE
     CHARACTER(LEN=255)          :: TPBC_DIR_SE
     INTEGER                     :: LWINDO_CU
     CHARACTER(LEN=255)          :: TPBC_DIR
     INTEGER                     :: NESTED_TS
     INTEGER                     :: NESTED_I1
     INTEGER                     :: NESTED_J1
     INTEGER                     :: NESTED_I2
     INTEGER                     :: NESTED_J2
     INTEGER                     :: NESTED_I0W
     INTEGER                     :: NESTED_J0W

     !----------------------------------------
     ! BENCHMARK MENU fields
     !----------------------------------------
     LOGICAL                     :: LSTDRUN
     CHARACTER(LEN=255)          :: STDRUN_INIT_FILE
     CHARACTER(LEN=255)          :: STDRUN_FINAL_FILE

     !----------------------------------------
     ! ARCHIVED OH MENU fields
     !----------------------------------------
     CHARACTER(LEN=255)          :: OH_DIR

     !----------------------------------------
     ! O3PL MENU fields
     !----------------------------------------
     CHARACTER(LEN=255)          :: O3PL_DIR

     !----------------------------------------
     ! MERCURY MENU fields
     !----------------------------------------     
     INTEGER                     :: ANTHRO_Hg_YEAR
     CHARACTER(LEN=255)          :: HG_SCENARIO
     LOGICAL                     :: USE_CHECKS
     LOGICAL                     :: LDYNOCEAN
     LOGICAL                     :: LPREINDHG
     CHARACTER(LEN=255)          :: Hg_RST_FILE
     LOGICAL                     :: LGTMM
     CHARACTER(LEN=255)          :: GTMM_RST_FILE

     !----------------------------------------
     ! CH4 MENU fields
     !----------------------------------------  
     LOGICAL                     :: LCH4BUD
     LOGICAL                     :: LGAO
     LOGICAL                     :: LCOL
     LOGICAL                     :: LLIV
     LOGICAL                     :: LWAST
     LOGICAL                     :: LBFCH4
     LOGICAL                     :: LRICE
     LOGICAL                     :: LOTANT
     LOGICAL                     :: LBMCH4
     LOGICAL                     :: LWETL
     LOGICAL                     :: LSOABS
     LOGICAL                     :: LOTNAT

     !----------------------------------------
     ! APM MENU fields
     !----------------------------------------  
     LOGICAL                     :: IFNUCL
     REAL*8                      :: FE0

     !----------------------------------------
     ! POPS MENU fields
     !----------------------------------------
     CHARACTER(LEN=3)           :: POP_TYPE
     LOGICAL                    :: CHEM_PROCESS
     CHARACTER(LEN=255)         :: POP_EMISFILE
     REAL*8                     :: POP_XMW
     REAL*8                     :: POP_KOA
     REAL*8                     :: POP_KBC
     REAL*8                     :: POP_K_POPG_OH
     REAL*8                     :: POP_K_POPP_O3A
     REAL*8                     :: POP_K_POPP_O3B
     REAL*8                     :: POP_HSTAR
     REAL*8                     :: POP_DEL_H
     REAL*8                     :: POP_DEL_Hw

     !----------------------------------------
     ! Fields for drydep and dust.  These get
     ! set in the init stage based on info 
     ! from file "input.geos". (mlong, 1/5/13)
     !----------------------------------------
     INTEGER                    :: NUMDEP
     INTEGER,           POINTER :: NDVZIND(:)
     INTEGER,           POINTER :: IDDEP(:)
     REAL*8,            POINTER :: DUSTREFF(:)
     REAL*8,            POINTER :: DUSTDEN(:)
     CHARACTER(LEN=14), POINTER :: DEPNAME(:)

     !----------------------------------------
     ! Fields for interface to GEOS-5 GCM
     !----------------------------------------
     LOGICAL                    :: haveImpRst
     INTEGER                    :: myCpu

     !----------------------------------------
     ! Fields for LINOZ strat chem
     !----------------------------------------
     INTEGER                    :: LINOZ_NLEVELS
     INTEGER                    :: LINOZ_NLAT
     INTEGER                    :: LINOZ_NMONTHS
     INTEGER                    :: LINOZ_NFIELDS
     REAL*8,            POINTER :: LINOZ_TPARM(:,:,:,:)

     !----------------------------------------
     ! Fields for overhead O3
     ! This gets set in main.F based on met
     ! field and year (mpayer, 12/13/13)
     !----------------------------------------
     LOGICAL                    :: USE_O3_FROM_MET

  END TYPE OptInput
!
! !REMARKS:
!  This will eventually replace the switches in logical_mod.F.  
!
! !REVISION HISTORY:
!  01 Nov 2012 - R. Yantosca - Initial version, based on logical_mod.F
!                              newer Olson 2001 land map & drydep inputs
!  07 Nov 2012 - R. Yantosca - Added Input_Opt%ITS_A_*_SIM fields
!  08 Nov 2012 - R. Yantosca - Added APM MENU fields
!  09 Nov 2012 - R. Yantosca - Added LD* variables for diagnostic levels
!  28 Nov 2012 - R. Yantosca - Add USE_OLSON_2001 logical flag
!  22 May 2013 - M. Payer    - Add GAMMA_HO2 variable for chemistry menu
!  26 Feb 2013 - M. Long     - Add extra fields from input.geos
!  26 Feb 2013 - M. Long     - Bug fix: timesteps are now INTEGER, not LOGICAL
!  28 Feb 2013 - R. Yantosca - Add haveImpRst field for GEOS-5 GCM interface
!  08 Mar 2013 - R. Yantosca - Add myCpu field to pass CPU # to GEOS-Chem
!  15 Mar 2013 - R. Yantosca - Add fields for LINOZ strat chemistry
!  27 Mar 2013 - R. Yantosca - Add extra fields for tagged CO2
!  27 Mar 2013 - R. Yantosca - Add extra fields for tagged EDGAR
!  29 Mar 2013 - R. Yantosca - Add DO_DIAG_WRITE field (to shut diags in MPI)
!  22 Jul 2013 - M. Sulprizio- Add extra fields for RCP emissions
!  31 Jul 2013 - M. Sulprizio- Add extra field for AEIC aircraft emissions and
!                              remove LAIRNOX field
!  13 Aug 2013 - M. Sulprizio- Add extra fields for semivolatile POA (H. Pye)
!  22 Aug 2013 - R. Yantosca - Add fields for soil NOx & species restart files
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57_DIR to GEOS_FP_DIR
!  03 Oct 2013 - M. Sulprizio- Removed obsolete LMFCT for flux correction
!  03 Oct 2013 - M. Sulprizio- Removed obsolete LAVHRRLAI and LMODISLAI
!  13 Dec 2013 - M. Sulprizio- Add USE_O3_FROM_MET logical flag
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_gigc_input_opt
!
! !DESCRIPTION: Subroutine INIT\_GIGC\_INPUT\_OPT intializes all GEOS-Chem
!  options carried in Input Options derived type object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_GIGC_Input_Opt( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE CMN_SIZE_Mod,     ONLY : NDSTBIN
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  Set the following fields of Input_Opt outside of this routine:
!  (1 ) Input_Opt%MAX_DIAG      : Max # of diagnostics
!  (2 ) Input_Opt%MAX_TRCS      : Max # of tracers
!  (3 ) Input_Opt%MAX_MEMB      : Max # of members per family tracer
!  (4 ) Input_Opt%MAX_FAMS      : Max # of P/L diagnostic families
!  (5 ) Input_Opt%MAX_DEP       : Max # of dry depositing species
!  (6 ) Input_Opt%LINOZ_NLEVELS : Number of levels    in LINOZ climatology
!  (7 ) Input_Opt%LINOZ_NLAT    : Number of latitudes in LINOZ climatology
!  (8 ) Input_Opt%LINOZ_NMONTHS : Number of months    in LINOZ climatology
!  (9 ) Input_Opt%LINOZ_NFIELDS : Number of species   in LINOZ climatology
!                                                                             .
!  We also need to implement better error checking.
!
! !REVISION HISTORY: 
!  01 Nov 2012 - R. Yantosca - Initial version
!  07 Nov 2012 - R. Yantosca - Now add size parameter fields to Input_Opt
!                              that can be set prior to calling this routine
!  09 Nov 2012 - R. Yantosca - Now zero LD* fields for diagnostic levels
!  28 Nov 2012 - R. Yantosca - Now set USE_OLSON_2001 logical flag
!  29 Nov 2012 - M. Payer    - Add Input_Opt%ITS_A_POPS_SIM
!  26 Feb 2013 - M. Long     - Add extra fields from input.geos
!  28 Feb 2013 - R. Yantosca - Add haveImpRst field for GEOS-5 GCM interface
!  08 Mar 2013 - R. Yantosca - Now initialize the myCpu field
!  15 Mar 2013 - R. Yantosca - Now initialize the LINOZ_TPARM field
!  27 Mar 2013 - R. Yantosca - Add extra fields for tagged CO2
!  27 Mar 2013 - R. Yantosca - Add extra fields for EDGAR
!  29 Mar 2013 - R. Yantosca - Add DO_DIAG_WRITE field (to shut diags in MPI)
!  22 Apr 2013 - R. Yantosca - Now dimension ND48_*ARR to 1000 so that we are
!                              consistent with the settings in diag48_mod.F
!  22 Jul 2013 - M. Sulprizio- Add extra fields for RCP emissions
!  07 Aug 2013 - M. Sulprizio- Add extra fields for SOA + SVPOA simulation
!  22 Aug 2013 - R. Yantosca - Add fields for soil NOx & species restart files
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57_DIR to GEOS_FP_DIR
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: MAX_DIAG, MAX_TRCS, MAX_MEMB, MAX_FAMS, MAX_DEP

    ! Assume success
    RC                               = GIGC_SUCCESS

    !----------------------------------------
    ! SIZE PARAMETER fields 
    !----------------------------------------
    MAX_DIAG                         = Input_Opt%MAX_DIAG
    MAX_TRCS                         = Input_Opt%MAX_TRCS
    MAX_MEMB                         = Input_Opt%MAX_MEMB
    MAX_FAMS                         = Input_Opt%MAX_FAMS
    MAX_DEP                          = Input_Opt%MAX_DEP
  
    !----------------------------------------
    ! SIMULATION MENU fields 
    !----------------------------------------
    Input_Opt%NYMDb                  = 0
    Input_Opt%NHMSb                  = 0
    Input_Opt%NYMDe                  = 0
    Input_Opt%NHMSe                  = 0
    Input_Opt%RUN_DIR                = ''
    Input_Opt%IN_RST_FILE            = ''
    Input_Opt%LSVGLB                 = .FALSE.
    Input_Opt%OUT_RST_FILE           = ''
    Input_Opt%DATA_DIR               = ''
    Input_Opt%GCAP_DIR               = ''
    Input_Opt%GEOS_4_DIR             = ''
    Input_Opt%GEOS_5_DIR             = ''
    Input_Opt%GEOS_FP_DIR            = ''
    Input_Opt%MERRA_DIR              = ''
    Input_Opt%DATA_DIR_1x1           = ''
    Input_Opt%TEMP_DIR               = ''
    Input_Opt%LUNZIP                 = .FALSE.
    Input_Opt%LWAIT                  = .FALSE.
    Input_Opt%LVARTROP               = .FALSE.
    Input_Opt%NESTED_I0              = 0
    Input_Opt%NESTED_J0              = 0
     
    !----------------------------------------
    ! TRACER MENU fields
    !----------------------------------------
    ALLOCATE( Input_Opt%ID_TRACER     ( MAX_TRCS           ), STAT=RC )
    ALLOCATE( Input_Opt%TRACER_NAME   ( MAX_TRCS           ), STAT=RC )
    ALLOCATE( Input_Opt%TRACER_MW_G   ( MAX_TRCS           ), STAT=RC )
    ALLOCATE( Input_Opt%TRACER_MW_KG  ( MAX_TRCS           ), STAT=RC )
    ALLOCATE( Input_Opt%TCVV          ( MAX_TRCS           ), STAT=RC )
    ALLOCATE( Input_Opt%XNUMOL        ( MAX_TRCS           ), STAT=RC )     
    ALLOCATE( Input_Opt%TRACER_N_CONST( MAX_TRCS           ), STAT=RC )     
    ALLOCATE( Input_Opt%TRACER_CONST  ( MAX_TRCS, MAX_MEMB ), STAT=RC )     
    ALLOCATE( Input_Opt%TRACER_COEFF  ( MAX_TRCS, MAX_MEMB ), STAT=RC )
    ALLOCATE( Input_Opt%ID_EMITTED    ( MAX_TRCS           ), STAT=RC )     

    Input_Opt%N_TRACERS              = 0
    Input_Opt%ID_TRACER              = 0
    Input_Opt%TRACER_NAME            = ''
    Input_Opt%TRACER_MW_G            = 0d0
    Input_Opt%TRACER_MW_KG           = 0d0
    Input_Opt%TCVV                   = 0d0
    Input_Opt%XNUMOL                 = 0d0
    Input_Opt%TRACER_N_CONST         = 0
    Input_Opt%TRACER_CONST           = ''  
    Input_Opt%TRACER_COEFF           = 0d0
    Input_Opt%ID_EMITTED             = 0
    Input_Opt%LSPLIT                 = .FALSE.
    Input_Opt%ITS_A_RnPbBe_SIM       = .FALSE.
    Input_Opt%ITS_A_CH3I_SIM         = .FALSE.
    Input_Opt%ITS_A_FULLCHEM_SIM     = .FALSE.
    Input_Opt%ITS_A_HCN_SIM          = .FALSE.
    Input_Opt%ITS_A_TAGOX_SIM        = .FALSE.
    Input_Opt%ITS_A_TAGCO_SIM        = .FALSE.
    Input_Opt%ITS_A_C2H6_SIM         = .FALSE.
    Input_Opt%ITS_A_CH4_SIM          = .FALSE.
    Input_Opt%ITS_AN_AEROSOL_SIM     = .FALSE.
    Input_Opt%ITS_A_MERCURY_SIM      = .FALSE.
    Input_Opt%ITS_A_CO2_SIM          = .FALSE.
    Input_Opt%ITS_A_H2HD_SIM         = .FALSE.
    Input_Opt%ITS_A_POPS_SIM         = .FALSE.
    Input_Opt%ITS_A_SPECIALTY_SIM    = .FALSE.
    Input_Opt%ITS_NOT_COPARAM_OR_CH4 = .FALSE.

    !----------------------------------------
    ! AEROSOL MENU fields
    !----------------------------------------
    ALLOCATE( Input_Opt%SALA_REDGE_um( 2 ), STAT=RC )
    ALLOCATE( Input_Opt%SALC_REDGE_um( 2 ), STAT=RC )     
    
    Input_Opt%LSULF                  = .FALSE.
    Input_Opt%LCRYST                 = .FALSE.
    Input_Opt%LCARB                  = .FALSE.
    Input_Opt%LSOA                   = .FALSE.
    Input_Opt%LSVPOA                 = .FALSE.
    Input_Opt%NAPEMISS               = 0d0
    Input_Opt%POAEMISSSCALE          = 0d0
    Input_Opt%LDUST                  = .FALSE.
    Input_Opt%LDEAD                  = .FALSE.
    Input_Opt%LSSALT                 = .FALSE.
    Input_Opt%LDICARB                = .FALSE.
    Input_Opt%SALA_REDGE_um          = 0d0
    Input_Opt%SALC_REDGE_um          = 0d0

    !----------------------------------------
    ! EMISSIONS MENU fields
    !----------------------------------------
    Input_Opt%LEMIS                  = .FALSE.
    Input_Opt%TS_EMIS                = 0
    Input_Opt%LANTHRO                = .FALSE.
    Input_Opt%FSCALYR                = .FALSE.
    Input_Opt%LEMEP                  = .FALSE.
    Input_Opt%LBRAVO                 = .FALSE.
    Input_Opt%LEDGAR                 = .FALSE.
    Input_Opt%LSTREETS               = .FALSE.
    Input_Opt%LCAC                   = .FALSE.
    Input_Opt%LNEI05                 = .FALSE.
    Input_Opt%LNEI08                 = .FALSE.
    Input_Opt%LRETRO                 = .FALSE.
    Input_Opt%LNEI99                 = .FALSE.
    Input_Opt%LICARTT                = .FALSE.
    Input_Opt%LVISTAS                = .FALSE.
    Input_Opt%LBIOFUEL               = .FALSE.
    Input_Opt%LBIOGENIC              = .FALSE.
    Input_Opt%LMEGAN                 = .FALSE.
    Input_Opt%LPECCA                 = .FALSE.
    Input_Opt%LMEGANMONO             = .FALSE.
    Input_Opt%ISOP_SCALING           = 0d0
    Input_Opt%LBIOMASS               = .FALSE.
    Input_Opt%LBBSEA                 = .FALSE.
    Input_Opt%LTOMSAI                = .FALSE.
    Input_Opt%LGFED2BB               = .FALSE.
    Input_Opt%L8DAYBB                = .FALSE.
    Input_Opt%L3HRBB                 = .FALSE.
    Input_Opt%LSYNOPBB               = .FALSE.
    Input_Opt%LGFED3BB               = .FALSE.
    Input_Opt%LDAYBB3                = .FALSE.
    Input_Opt%L3HRBB3                = .FALSE.
    Input_Opt%LAEIC                  = .FALSE.
    Input_Opt%LLIGHTNOX              = .FALSE.
    Input_Opt%LOTDLOC                = .FALSE.
    Input_Opt%LSOILNOX               = .FALSE.
    Input_Opt%SOIL_RST_FILE          = ''
    Input_Opt%LFERTILIZERNOX         = .FALSE.
    Input_Opt%NOx_SCALING            = 0d0
    Input_Opt%LEDGARSHIP             = .FALSE.
    Input_Opt%LICOADSSHIP            = .FALSE.
    Input_Opt%LEMEPSHIP              = .FALSE.
    Input_Opt%LSHIPSO2               = .FALSE.
    Input_Opt%LARCSHIP               = .FALSE.
    Input_Opt%LCOOKE                 = .FALSE.
    Input_Opt%LHIST                  = .FALSE.
    Input_Opt%HISTYR                 = .FALSE.
    Input_Opt%LWARWICK_VSLS          = .FALSE.
    Input_Opt%LSSABr2                = .FALSE.
    Input_Opt%LFIX_PBL_BRO           = .FALSE.
    Input_Opt%Br_SCALING             = 0d0    
    Input_Opt%LEDGARNOx              = .FALSE.
    Input_Opt%LEDGARCO               = .FALSE. 
    Input_Opt%LEDGARSOX              = .FALSE.
    Input_Opt%LRCP                   = .FALSE.
    Input_Opt%LRCPSHIP               = .FALSE.
    Input_Opt%LRCPAIR                = .FALSE.

    !----------------------------------------
    ! CO2 MENU fields
    !----------------------------------------
    Input_Opt%LGENFF                 = .FALSE.
    Input_Opt%LANNFF                 = .FALSE.
    Input_Opt%LMONFF                 = .FALSE.
    Input_Opt%LCHEMCO2               = .FALSE.
    Input_Opt%LSEASBB                = .FALSE.
    Input_Opt%LBIOFUEL               = .FALSE.
    Input_Opt%LBIODAILY              = .FALSE.
    Input_Opt%LBIODIURNAL            = .FALSE.
    Input_Opt%LBIONETORIG            = .FALSE.
    Input_Opt%LBIONETCLIM            = .FALSE.
    Input_Opt%LOCN1997               = .FALSE.
    Input_Opt%LOCN2009ANN            = .FALSE.
    Input_Opt%LOCN2009MON            = .FALSE.
    Input_Opt%LSHIPEDG               = .FALSE.
    Input_Opt%LSHIPICO               = .FALSE.
    Input_Opt%LPLANE                 = .FALSE.
    Input_Opt%LFFBKGRD               = .FALSE.
    Input_Opt%LBIOSPHTAG             = .FALSE.
    Input_Opt%LFOSSILTAG             = .FALSE.
    Input_Opt%LSHIPTAG               = .FALSE.
    Input_Opt%LPLANETAG              = .FALSE.

    !----------------------------------------
    ! FUTURE MENU fields
    !----------------------------------------
    Input_Opt%LFUTURE                = .FALSE.
    Input_Opt%FUTURE_YEAR            = 0
    Input_Opt%FUTURE_SCEN            = ''

    !----------------------------------------
    ! CHEMISTRY MENU fields
    !----------------------------------------
    Input_Opt%LCHEM                  = .FALSE.
    Input_Opt%LSCHEM                 = .FALSE.
    Input_Opt%LLINOZ                 = .FALSE. 
    Input_Opt%TS_CHEM                = 0
    Input_Opt%LSVCSPEC               = .FALSE. 
    Input_Opt%SPEC_RST_FILE          = ''
    Input_Opt%LKPP                   = .FALSE. 
    Input_Opt%GAMMA_HO2              = 0d0

    !----------------------------------------
    ! CHEMISTRY MENU fields
    !----------------------------------------
    Input_Opt%LTRAN                  = .FALSE.
    Input_Opt%LFILL                  = .FALSE.
    Input_Opt%TPCORE_IORD            = .FALSE.
    Input_Opt%TPCORE_JORD            = .FALSE.
    Input_Opt%TPCORE_KORD            = .FALSE.
    Input_Opt%TS_DYN                 = 0

    !----------------------------------------
    ! CONVECTION MENU fields
    !----------------------------------------
    Input_Opt%LCONV                  = .FALSE.
    Input_Opt%LTURB                  = .FALSE.
    Input_Opt%LNLPBL                 = .FALSE.
    Input_Opt%TS_CONV                = 0

    !----------------------------------------
    ! DEPOSITION MENU fields
    !----------------------------------------
    Input_Opt%LDRYD                  = .FALSE.
    Input_Opt%LWETD                  = .FALSE.
    Input_Opt%USE_OLSON_2001         = .FALSE.

    !----------------------------------------
    ! GAMAP_MENU fields
    !----------------------------------------
    Input_Opt%GAMAP_DIAGINFO         = ''
    Input_Opt%GAMAP_TRACERINFO       = ''

    !----------------------------------------
    ! OUTPUT MENU fields
    !----------------------------------------
    ALLOCATE( Input_Opt%NJDAY( 366 ), STAT=RC )

    Input_Opt%NJDAY                  = 0

    !----------------------------------------
    ! DIAGNOSTIC MENU fields
    !----------------------------------------
    Input_Opt%TS_DIAG                = 0
    ALLOCATE( Input_Opt%TINDEX( MAX_DIAG, MAX_TRCS ), STAT=RC )
    ALLOCATE( Input_Opt%TCOUNT( MAX_DIAG           ), STAT=RC )
    ALLOCATE( Input_Opt%TMAX  ( MAX_DIAG           ), STAT=RC )

    Input_Opt%ND01                   = 0
    Input_Opt%ND02                   = 0
    Input_Opt%ND03                   = 0
    Input_Opt%ND04                   = 0
    Input_Opt%ND05                   = 0
    Input_Opt%ND06                   = 0
    Input_Opt%ND07                   = 0
    Input_Opt%ND08                   = 0
    Input_Opt%ND09                   = 0
    Input_Opt%ND10                   = 0
    Input_Opt%ND11                   = 0
    Input_Opt%ND12                   = 0
    Input_Opt%ND13                   = 0
    Input_Opt%ND14                   = 0
    Input_Opt%ND15                   = 0
    Input_Opt%ND16                   = 0
    Input_Opt%ND17                   = 0
    Input_Opt%ND18                   = 0
    Input_Opt%ND19                   = 0
    Input_Opt%ND20                   = 0
    Input_Opt%ND21                   = 0
    Input_Opt%ND22                   = 0
    Input_Opt%ND23                   = 0
    Input_Opt%ND24                   = 0
    Input_Opt%ND25                   = 0
    Input_Opt%ND26                   = 0
    Input_Opt%ND27                   = 0
    Input_Opt%ND28                   = 0
    Input_Opt%ND29                   = 0
    Input_Opt%ND30                   = 0
    Input_Opt%ND31                   = 0
    Input_Opt%ND32                   = 0
    Input_Opt%ND33                   = 0
    Input_Opt%ND34                   = 0
    Input_Opt%ND35                   = 0
    Input_Opt%ND36                   = 0
    Input_Opt%ND37                   = 0
    Input_Opt%ND38                   = 0
    Input_Opt%ND39                   = 0
    Input_Opt%ND40                   = 0
    Input_Opt%ND41                   = 0
    Input_Opt%ND42                   = 0
    Input_Opt%ND43                   = 0
    Input_Opt%ND44                   = 0
    Input_Opt%ND45                   = 0
    Input_Opt%ND46                   = 0
    Input_Opt%ND47                   = 0
    Input_Opt%ND48                   = 0
    Input_Opt%ND49                   = 0
    Input_Opt%ND50                   = 0
    Input_Opt%ND51                   = 0
    Input_Opt%ND52                   = 0
    Input_Opt%ND53                   = 0
    Input_Opt%ND54                   = 0
    Input_Opt%ND55                   = 0
    Input_Opt%ND56                   = 0
    Input_Opt%ND57                   = 0
    Input_Opt%ND58                   = 0
    Input_Opt%ND59                   = 0
    Input_Opt%ND60                   = 0
    Input_Opt%ND61                   = 0
    Input_Opt%ND62                   = 0
    Input_Opt%ND63                   = 0
    Input_Opt%ND64                   = 0
    Input_Opt%ND65                   = 0
    Input_Opt%ND66                   = 0
    Input_Opt%ND67                   = 0
    Input_Opt%ND68                   = 0
    Input_Opt%ND69                   = 0
    Input_Opt%ND70                   = 0
    Input_Opt%LD01                   = 0
    Input_Opt%LD02                   = 0
    Input_Opt%LD03                   = 0
    Input_Opt%LD04                   = 0
    Input_Opt%LD05                   = 0
    Input_Opt%LD06                   = 0
    Input_Opt%LD07                   = 0
    Input_Opt%LD08                   = 0
    Input_Opt%LD09                   = 0
    Input_Opt%LD10                   = 0
    Input_Opt%LD11                   = 0
    Input_Opt%LD12                   = 0
    Input_Opt%LD13                   = 0
    Input_Opt%LD14                   = 0
    Input_Opt%LD15                   = 0
    Input_Opt%LD16                   = 0
    Input_Opt%LD17                   = 0
    Input_Opt%LD18                   = 0
    Input_Opt%LD19                   = 0
    Input_Opt%LD20                   = 0
    Input_Opt%LD21                   = 0
    Input_Opt%LD22                   = 0
    Input_Opt%LD23                   = 0
    Input_Opt%LD24                   = 0
    Input_Opt%LD25                   = 0
    Input_Opt%LD26                   = 0
    Input_Opt%LD27                   = 0
    Input_Opt%LD28                   = 0
    Input_Opt%LD29                   = 0
    Input_Opt%LD30                   = 0
    Input_Opt%LD31                   = 0
    Input_Opt%LD32                   = 0
    Input_Opt%LD33                   = 0
    Input_Opt%LD34                   = 0
    Input_Opt%LD35                   = 0
    Input_Opt%LD36                   = 0
    Input_Opt%LD37                   = 0
    Input_Opt%LD38                   = 0
    Input_Opt%LD39                   = 0
    Input_Opt%LD40                   = 0
    Input_Opt%LD41                   = 0
    Input_Opt%LD42                   = 0
    Input_Opt%LD43                   = 0
    Input_Opt%LD44                   = 0
    Input_Opt%LD45                   = 0
    Input_Opt%LD46                   = 0
    Input_Opt%LD47                   = 0
    Input_Opt%LD48                   = 0
    Input_Opt%LD49                   = 0
    Input_Opt%LD50                   = 0
    Input_Opt%LD51                   = 0
    Input_Opt%LD52                   = 0
    Input_Opt%LD53                   = 0
    Input_Opt%LD54                   = 0
    Input_Opt%LD55                   = 0
    Input_Opt%LD56                   = 0
    Input_Opt%LD57                   = 0
    Input_Opt%LD58                   = 0
    Input_Opt%LD59                   = 0
    Input_Opt%LD60                   = 0
    Input_Opt%LD61                   = 0
    Input_Opt%LD62                   = 0
    Input_Opt%LD63                   = 0
    Input_Opt%LD64                   = 0
    Input_Opt%LD65                   = 0
    Input_Opt%LD66                   = 0
    Input_Opt%LD67                   = 0
    Input_Opt%LD68                   = 0
    Input_Opt%LD69                   = 0
    Input_Opt%LD70                   = 0
    Input_Opt%LPRT                   = .FALSE.
    Input_Opt%TINDEX(:,:)            = 0
    Input_Opt%TCOUNT(:)              = 0	  
    Input_Opt%TMAX(:)	             = 0
#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
    ! Need to shut off G-C diagnostics when 
    ! connecting to an external GCM (bmy, 3/29/13)
    Input_Opt%DO_DIAG_WRITE          = .FALSE.
#else
    ! For traditional G-C runs, always write diags (bmy, 3/29/13)
    Input_Opt%DO_DIAG_WRITE          = .TRUE.
#endif   

    !----------------------------------------
    ! PLANEFLIGHT MENU fields
    !----------------------------------------
    Input_Opt%DO_PF                  = .FALSE.
    Input_Opt%PF_IFILE               = ''
    Input_Opt%PF_OFILE               = ''

    !----------------------------------------
    ! ND48 MENU fields
    !----------------------------------------
    ALLOCATE( Input_Opt%ND48_IARR( 1000 ), STAT=RC )
    ALLOCATE( Input_Opt%ND48_JARR( 1000 ), STAT=RC )
    ALLOCATE( Input_Opt%ND48_LARR( 1000 ), STAT=RC )
    ALLOCATE( Input_Opt%ND48_NARR( 1000 ), STAT=RC )

    Input_Opt%DO_ND48                = .FALSE.
    Input_Opt%ND48_FILE              = ''
    Input_Opt%ND48_FREQ              = 0
    Input_Opt%ND48_N_STA             = 0
    Input_Opt%ND48_IARR              = 0
    Input_Opt%ND48_JARR              = 0
    Input_Opt%ND48_LARR              = 0
    Input_Opt%ND48_NARR              = 0

    !----------------------------------------
    ! ND49 MENU fields
    !----------------------------------------
    ALLOCATE( Input_Opt%ND49_TRACERS( MAX_TRCS ), STAT=RC )

    Input_Opt%DO_ND49                = .FALSE.
    Input_Opt%ND49_FILE              = ''
    Input_Opt%ND49_TRACERS           = 0
    Input_Opt%ND49_FREQ              = 0
    Input_Opt%ND49_IMIN              = 0
    Input_Opt%ND49_IMAX              = 0
    Input_Opt%ND49_JMIN              = 0
    Input_Opt%ND49_JMAX              = 0
    Input_Opt%ND49_LMIN              = 0
    Input_Opt%ND49_LMAX              = 0

    !----------------------------------------
    ! ND50 MENU fields
    !----------------------------------------
    ALLOCATE( Input_Opt%ND50_TRACERS( MAX_TRCS ), STAT=RC )

    Input_Opt%DO_ND50                = .FALSE.
    Input_Opt%ND50_FILE              = ''
    Input_Opt%LND50_HDF              = .FALSE.
    Input_Opt%ND50_TRACERS           = 0
    Input_Opt%ND50_IMIN              = 0
    Input_Opt%ND50_IMAX              = 0
    Input_Opt%ND50_JMIN              = 0
    Input_Opt%ND50_JMAX              = 0
    Input_Opt%ND50_LMIN              = 0
    Input_Opt%ND50_LMAX              = 0

    !----------------------------------------
    ! ND51 MENU fields
    !----------------------------------------
    ALLOCATE( Input_Opt%ND51_TRACERS( MAX_TRCS ), STAT=RC )

    Input_Opt%DO_ND51                = .FALSE.
    Input_Opt%ND51_FILE              = ''
    Input_Opt%LND51_HDF              = .FALSE.
    Input_Opt%ND51_TRACERS           = 0
    Input_Opt%ND51_HR_WRITE          = 0d0
    Input_Opt%ND51_HR1               = 0d0
    Input_Opt%ND51_HR2               = 0d0
    Input_Opt%ND51_IMIN              = 0
    Input_Opt%ND51_IMAX              = 0
    Input_Opt%ND51_JMIN              = 0
    Input_Opt%ND51_JMAX              = 0
    Input_Opt%ND51_LMIN              = 0

    !----------------------------------------
    ! ND51b MENU fields
    !----------------------------------------
    ALLOCATE( Input_Opt%ND51b_TRACERS( MAX_TRCS ), STAT=RC )

    Input_Opt%DO_ND51b               = .FALSE.
    Input_Opt%ND51b_FILE             = ''
    Input_Opt%LND51b_HDF             = .FALSE.
    Input_Opt%ND51b_TRACERS          = 0
    Input_Opt%ND51b_HR_WRITE         = 0d0
    Input_Opt%ND51b_HR1              = 0d0
    Input_Opt%ND51b_HR2              = 0d0
    Input_Opt%ND51b_IMIN             = 0
    Input_Opt%ND51b_IMAX             = 0
    Input_Opt%ND51b_JMIN             = 0
    Input_Opt%ND51b_JMAX             = 0
    Input_Opt%ND51b_LMIN             = 0

    !----------------------------------------
    ! ND63 MENU fields
    !----------------------------------------
    ALLOCATE( Input_Opt%ND63_TRACERS( MAX_TRCS ), STAT=RC )

    Input_Opt%DO_ND63                = .FALSE.
    Input_Opt%ND63_FILE              = ''
    Input_Opt%ND63_TRACERS           = 0
    Input_Opt%ND63_FREQ              = 0
    Input_Opt%ND63_IMIN              = 0
    Input_Opt%ND63_IMAX              = 0
    Input_Opt%ND63_JMIN              = 0
    Input_Opt%ND63_JMAX              = 0

    !----------------------------------------
    ! PROD LOSS MENU fields
    !----------------------------------------
    ALLOCATE( Input_Opt%FAM_COEF( MAX_MEMB, MAX_FAMS ), STAT=RC )
    ALLOCATE( Input_Opt%FAM_MEMB( MAX_MEMB, MAX_FAMS ), STAT=RC )
    ALLOCATE( Input_Opt%FAM_NAME(           MAX_FAMS ), STAT=RC )
    ALLOCATE( Input_Opt%FAM_NMEM(           MAX_FAMS ), STAT=RC )
    ALLOCATE( Input_Opt%FAM_TYPE(           MAX_FAMS ), STAT=RC )

    Input_Opt%DO_SAVE_PL             = .FALSE.
    Input_Opt%LFAMILY                = .FALSE.
    Input_Opt%ND65                   = 0
    Input_Opt%DO_SAVE_O3             = .FALSE.
    Input_Opt%NFAM                   = 0
    Input_Opt%FAM_COEF               = 0d0
    Input_Opt%FAM_MEMB               = ''
    Input_Opt%FAM_NAME               = ''
    Input_Opt%FAM_NMEM               = 0
    Input_Opt%FAM_TYPE               = ''

    !----------------------------------------
    ! UNIX CMDS fields
    !----------------------------------------
    Input_Opt%BACKGROUND             = ''
    Input_Opt%REDIRECT               = ''  
    Input_Opt%REMOVE_CMD             = ''
    Input_Opt%SEPARATOR              = ''
    Input_Opt%WILD_CARD              = ''
    Input_Opt%UNZIP_CMD              = ''
    Input_Opt%ZIP_SUFFIX             = ''

    !----------------------------------------
    ! NESTED GRID MENU fields
    !----------------------------------------
    Input_Opt%LWINDO                 = .FALSE.
    Input_Opt%LWINDO2x25             = .FALSE.
    Input_Opt%LWINDO_NA              = .FALSE.
    Input_Opt%TPBC_DIR_NA            = ''
    Input_Opt%LWINDO_EU              = .FALSE.
    Input_Opt%TPBC_DIR_EU            = ''
    Input_Opt%LWINDO_CH              = .FALSE.
    Input_Opt%TPBC_DIR_CH            = ''
    Input_Opt%LWINDO_SE              = .FALSE.
    Input_Opt%TPBC_DIR_SE            = ''
    Input_Opt%LWINDO_CU              = 0
    Input_Opt%TPBC_DIR               = ''
    Input_Opt%NESTED_TS              = 0
    Input_Opt%NESTED_I1              = 0
    Input_Opt%NESTED_J1              = 0
    Input_Opt%NESTED_I2              = 0
    Input_Opt%NESTED_J2              = 0
    Input_Opt%NESTED_I0W             = 0
    Input_Opt%NESTED_J0W             = 0 

    !----------------------------------------
    ! BENCHMARK MENU fields
    !----------------------------------------
    Input_Opt%LSTDRUN                = .FALSE.
    Input_Opt%STDRUN_INIT_FILE       = ''
    Input_Opt%STDRUN_FINAL_FILE      =''

    !----------------------------------------
    ! ARCHIVED OH MENU fields
    !----------------------------------------
    Input_Opt%OH_DIR                 = ''

    !----------------------------------------
    ! O3PL MENU fields
    !----------------------------------------
    Input_Opt%O3PL_DIR               = ''

    !----------------------------------------
    ! MERCURY MENU fields
    !----------------------------------------     
    Input_Opt%ANTHRO_Hg_YEAR         = 0
    Input_Opt%HG_SCENARIO            = ''
    Input_Opt%USE_CHECKS             = .FALSE.
    Input_Opt%LDYNOCEAN              = .FALSE.
    Input_Opt%LPREINDHG              = .FALSE.
    Input_Opt%Hg_RST_FILE            = ''
    Input_Opt%LGTMM                  = .FALSE.
    Input_Opt%GTMM_RST_FILE          = ''

    !----------------------------------------
    ! CH4 MENU fields
    !----------------------------------------  
    Input_Opt%LCH4BUD                = .FALSE.
    Input_Opt%LGAO                   = .FALSE.
    Input_Opt%LCOL                   = .FALSE.
    Input_Opt%LLIV                   = .FALSE.
    Input_Opt%LWAST                  = .FALSE.
    Input_Opt%LBFCH4                 = .FALSE.
    Input_Opt%LRICE                  = .FALSE.
    Input_Opt%LOTANT                 = .FALSE.
    Input_Opt%LBMCH4                 = .FALSE.
    Input_Opt%LWETL                  = .FALSE.
    Input_Opt%LSOABS                 = .FALSE.
    Input_Opt%LOTNAT                 = .FALSE.

    !----------------------------------------
    ! APM MENU fields
    !----------------------------------------  
    Input_Opt%IFNUCL                 = .FALSE.
    Input_Opt%FE0                    = 0d0

    !----------------------------------------
    ! POPS MENU fields
    !----------------------------------------
    Input_Opt%POP_TYPE               = ''
    Input_Opt%CHEM_PROCESS           = .FALSE.
    Input_Opt%POP_EMISFILE           = ''
    Input_Opt%POP_XMW                = 0d0
    Input_Opt%POP_KOA                = 0d0
    Input_Opt%POP_KBC                = 0d0
    Input_Opt%POP_K_POPG_OH          = 0d0
    Input_Opt%POP_K_POPP_O3A         = 0d0
    Input_Opt%POP_K_POPP_O3B         = 0d0
    Input_Opt%POP_HSTAR              = 0d0
    Input_Opt%POP_DEL_H              = 0d0
    Input_Opt%POP_DEL_Hw             = 0d0

    !----------------------------------------
    ! Fields for DRYDEP and DUST based on
    ! input from the file "input.geos"
    !----------------------------------------
    ALLOCATE( Input_Opt%NDVZIND ( MAX_DEP ), STAT=RC ) ! Drydep
    ALLOCATE( Input_Opt%DEPNAME ( MAX_DEP ), STAT=RC ) ! Drydep
    ALLOCATE( Input_Opt%IDDEP   ( NDSTBIN ), STAT=RC ) ! Dust_mod
    ALLOCATE( Input_Opt%DUSTREFF( NDSTBIN ), STAT=RC ) ! Dust_mod
    ALLOCATE( Input_Opt%DUSTDEN ( NDSTBIN ), STAT=RC ) ! Dust_mod

    Input_Opt%NUMDEP                 = 0
    Input_Opt%NDVZIND                = 0
    Input_Opt%IDDEP                  = 0
    Input_Opt%DUSTREFF               = 0d0
    Input_Opt%DUSTDEN                = 0d0
    Input_Opt%DEPNAME                = ''

    !----------------------------------------
    ! Fields for interface to GEOS-5 GCM
    !----------------------------------------
    Input_Opt%haveImpRst             = .FALSE.
    Input_Opt%myCpu                  = -1


    !----------------------------------------
    ! Fields for LINOZ strat chem
    !----------------------------------------
    ALLOCATE( Input_Opt%LINOZ_TPARM( Input_Opt%LINOZ_NLEVELS,            &
                                     Input_Opt%LINOZ_NLAT,               &
                                     Input_Opt%LINOZ_NMONTHS,            &
                                     Input_Opt%LINOZ_NFIELDS ), STAT=RC )

    Input_Opt%LINOZ_TPARM            = 0d0

    !----------------------------------------
    ! Fields for overhead O3
    !----------------------------------------
    Input_Opt%USE_O3_FROM_MET        = .FALSE.

  END SUBROUTINE Set_GIGC_Input_Opt
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_gigc_input_opt
!
! !DESCRIPTION: Subroutine CLEANUP\_GIGC\_INPUT\_OPT deallocates all 
!  allocatable fields of the Input Options object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_GIGC_Input_Opt( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY: 
!  02 Nov 2012 - R. Yantosca - Initial version
!  07 Nov 2012 - R. Yantosca - Now deallocate fields from prod/loss menu
!  26 Feb 2013 - M. Long     - Now deallocate extra fields from input.geos
!  15 Mar 2013 - R. Yantosca - Now deallocate the LINOZ_TPARM field
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Assume success
    RC = GIGC_SUCCESS

    !======================================================================
    ! Deallocate fields of the Input Options object
    !======================================================================
    IF ( ASSOCIATED( Input_Opt%ID_TRACER ) ) THEN
       DEALLOCATE( Input_Opt%ID_TRACER  ) 
    ENDIF

    IF ( ASSOCIATED( Input_Opt%TRACER_NAME ) ) THEN
       DEALLOCATE( Input_Opt%TRACER_NAME )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%TRACER_MW_G ) ) THEN
       DEALLOCATE( Input_Opt%TRACER_MW_G )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%TRACER_MW_KG ) ) THEN
       DEALLOCATE( Input_Opt%TRACER_MW_KG   )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%TCVV ) ) THEN
       DEALLOCATE( Input_Opt%TCVV )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%XNUMOL ) ) THEN
       DEALLOCATE( Input_Opt%XNUMOL )
    ENDIF

    IF ( ASSOCIATED (Input_Opt%TRACER_N_CONST ) ) THEN
       DEALLOCATE( Input_Opt%TRACER_N_CONST )
    ENDIF

    IF ( ASSOCIATED (Input_Opt%TRACER_CONST ) ) THEN
       DEALLOCATE( Input_Opt%TRACER_CONST )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%TRACER_COEFF ) ) THEN
       DEALLOCATE( Input_Opt%TRACER_COEFF )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%ID_EMITTED ) ) THEN
       DEALLOCATE( Input_Opt%ID_EMITTED )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%SALA_REDGE_um ) ) THEN
       DEALLOCATE( Input_Opt%SALA_REDGE_um  )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%SALC_REDGE_um ) ) THEN 
       DEALLOCATE( Input_Opt%SALC_REDGE_um  )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%NJDAY ) ) THEN
       DEALLOCATE( Input_Opt%NJDAY )
    ENDIF
    
    IF ( ASSOCIATED( Input_Opt%TINDEX ) ) THEN
       DEALLOCATE( Input_Opt%TINDEX )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%TCOUNT ) ) THEN
       DEALLOCATE( Input_Opt%TCOUNT )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%TMAX ) ) THEN
       DEALLOCATE( Input_Opt%TMAX )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%ND48_IARR ) ) THEN
       DEALLOCATE( Input_Opt%ND48_IARR )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%ND48_JARR ) ) THEN
       DEALLOCATE( Input_Opt%ND48_JARR )
    ENDIF
     
    IF ( ASSOCIATED( Input_Opt%ND48_LARR ) ) THEN
       DEALLOCATE( Input_Opt%ND48_LARR )
    ENDIF
    
    IF ( ASSOCIATED( Input_Opt%ND48_NARR) ) THEN
       DEALLOCATE( Input_Opt%ND48_NARR )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%ND49_TRACERS ) ) THEN
       DEALLOCATE( Input_Opt%ND49_TRACERS )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%ND50_TRACERS ) ) THEN
       DEALLOCATE( Input_Opt%ND50_TRACERS )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%ND51_TRACERS ) ) THEN
       DEALLOCATE( Input_Opt%ND51_TRACERS )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%ND51b_TRACERS ) ) THEN
       DEALLOCATE( Input_Opt%ND51b_TRACERS )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%ND63_TRACERS ) ) THEN
       DEALLOCATE( Input_Opt%ND63_TRACERS )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%FAM_COEF ) ) THEN
       DEALLOCATE( Input_Opt%FAM_COEF )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%FAM_MEMB ) ) THEN
       DEALLOCATE( Input_Opt%FAM_MEMB )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%FAM_NAME ) ) THEN
       DEALLOCATE( Input_Opt%FAM_NAME )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%FAM_NMEM ) ) THEN
       DEALLOCATE( Input_Opt%FAM_NMEM )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%FAM_TYPE ) ) THEN
       DEALLOCATE( Input_Opt%FAM_TYPE )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%NDVZIND ) ) THEN
       DEALLOCATE( Input_Opt%NDVZIND )
    ENDIF
    
    IF ( ASSOCIATED( Input_Opt%DEPNAME ) ) THEN
       DEALLOCATE( Input_Opt%DEPNAME )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%IDDEP ) ) THEN
       DEALLOCATE( Input_Opt%IDDEP )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%DUSTREFF ) ) THEN
       DEALLOCATE( Input_Opt%DUSTREFF )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%DUSTDEN ) ) THEN
       DEALLOCATE( Input_Opt%DUSTDEN )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%LINOZ_TPARM ) ) THEN
       DEALLOCATE( Input_Opt%LINOZ_TPARM )
    ENDIF

  END SUBROUTINE Cleanup_GIGC_Input_Opt
!EOC
END MODULE GIGC_Input_Opt_Mod
