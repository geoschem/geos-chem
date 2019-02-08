!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: input_opt_mod.F90
!
! !DESCRIPTION: Module INPUT\_OPT\_MOD contains the derived type for GEOS-Chem
!  options and logical switches.
!\\
!\\
! !INTERFACE:
!
MODULE Input_Opt_Mod
!
! !USES:
!
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Set_Input_Opt
  PUBLIC :: Set_Input_Opt_Advect
  PUBLIC :: Cleanup_Input_Opt
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for Input Options 
  !=========================================================================
  TYPE, PUBLIC :: OptInput

     !----------------------------------------
     ! General Runtime & Distributed Comp Info
     !----------------------------------------
     INTEGER                     :: NPES      ! Number of MPI procs
     INTEGER                     :: myCpu     ! Local MPI process handle
     INTEGER                     :: MPICOMM   ! MPI Communicator Handle
     LOGICAL                     :: HPC       ! Is this an HPC (ESMF or otherwise) sim?
     LOGICAL                     :: RootCPU   ! Is this the root cpu?

     !----------------------------------------
     ! SIZE PARAMETER fields
     !----------------------------------------
     INTEGER                     :: Max_BPCH_Diag
     INTEGER                     :: Max_Families
     INTEGER                     :: Max_AdvectSpc
     INTEGER                     :: Max_PassiveSpc

     !----------------------------------------
     ! SIMULATION MENU fields 
     !----------------------------------------
     INTEGER                     :: NYMDb              
     INTEGER                     :: NHMSb              
     INTEGER                     :: NYMDe              
     INTEGER                     :: NHMSe              
     INTEGER                     :: SimLengthSec
     CHARACTER(LEN=255)          :: RUN_DIR            
     CHARACTER(LEN=255)          :: DATA_DIR           
     CHARACTER(LEN=255)          :: CHEM_INPUTS_DIR
     CHARACTER(LEN=255)          :: RES_DIR
     LOGICAL                     :: LCAPTROP
     REAL(fp)                    :: OZONOPAUSE
     INTEGER                     :: NESTED_I0          
     INTEGER                     :: NESTED_J0          
     CHARACTER(LEN=255)          :: HcoConfigFile

     !----------------------------------------
     ! PASSIVE SPECIES MENU fields
     !----------------------------------------
     INTEGER                     :: NPASSIVE
     INTEGER                     :: NPASSIVE_DECAY
     CHARACTER(LEN=63),  POINTER :: PASSIVE_NAME    (:)
     CHARACTER(LEN=255), POINTER :: PASSIVE_LONGNAME(:)
     INTEGER,            POINTER :: PASSIVE_ID      (:)
     REAL(fp),           POINTER :: PASSIVE_MW      (:)
     REAL(fp),           POINTER :: PASSIVE_TAU     (:)
     REAL(fp),           POINTER :: PASSIVE_INITCONC(:)
     INTEGER,            POINTER :: PASSIVE_DECAYID (:)

     !----------------------------------------
     ! ADVECTED SPECIES MENU fields
     !----------------------------------------
     INTEGER                     :: N_ADVECT
     CHARACTER(LEN=255), POINTER :: AdvectSpc_Name(:) 
     INTEGER                     :: SIM_TYPE
     CHARACTER(LEN=255)          :: SIM_NAME
     LOGICAL                     :: LSPLIT
     LOGICAL                     :: ITS_A_RnPbBe_SIM
     LOGICAL                     :: ITS_A_CH3I_SIM
     LOGICAL                     :: ITS_A_FULLCHEM_SIM
     LOGICAL                     :: ITS_A_HCN_SIM
     LOGICAL                     :: ITS_A_TAGO3_SIM
     LOGICAL                     :: ITS_A_TAGCO_SIM
     LOGICAL                     :: ITS_A_C2H6_SIM
     LOGICAL                     :: ITS_A_CH4_SIM
     LOGICAL                     :: ITS_AN_AEROSOL_SIM
     LOGICAL                     :: ITS_A_MERCURY_SIM
     LOGICAL                     :: ITS_A_CO2_SIM
     LOGICAL                     :: ITS_A_H2HD_SIM
     LOGICAL                     :: ITS_A_POPS_SIM

     !----------------------------------------
     ! AEROSOL MENU fields
     !----------------------------------------
     LOGICAL                     :: LSULF              
     LOGICAL                     :: LMETALCATSO2
     LOGICAL                     :: LCARB
     LOGICAL                     :: LBRC              
     LOGICAL                     :: LSOA
     LOGICAL                     :: LMPOA
     LOGICAL                     :: LSVPOA
     LOGICAL                     :: LOMOC
     LOGICAL                     :: LDUST              
     LOGICAL                     :: LDEAD              
     LOGICAL                     :: LSSALT             
     LOGICAL                     :: LDSTUP
     REAL(fp),           POINTER :: SALA_REDGE_um(:)   
     REAL(fp),           POINTER :: SALC_REDGE_um(:)   
     LOGICAL                     :: LGRAVSTRAT
     LOGICAL                     :: LSOLIDPSC
     LOGICAL                     :: LHOMNUCNAT
     REAL(fp)                    :: T_NAT_SUPERCOOL
     REAL(fp)                    :: P_ICE_SUPERSAT
     LOGICAL                     :: LPSCCHEM
     LOGICAL                     :: LSTRATOD
     !for BC absorption enhancement, (xnw, 8/24/15)
     LOGICAL                     :: LBCAE
     REAL(fp)                    :: BCAE_1
     REAL(fp)                    :: BCAE_2

     !----------------------------------------
     ! EMISSIONS MENU fields
     !----------------------------------------
     LOGICAL                     :: LEMIS
     INTEGER                     :: TS_EMIS
     LOGICAL                     :: LBIOFUEL
     LOGICAL                     :: LOTDLOC
     LOGICAL                     :: LSOILNOX
     LOGICAL                     :: LWARWICK_VSLS
     LOGICAL                     :: LSSABr2
     LOGICAL                     :: LFIX_PBL_BRO
     LOGICAL                     :: LCH4EMIS
     LOGICAL                     :: LCH4SBC
     LOGICAL                     :: LOCSEMIS
     LOGICAL                     :: LCFCEMIS
     LOGICAL                     :: LCLEMIS
     LOGICAL                     :: LBREMIS
     LOGICAL                     :: LN2OEMIS
     LOGICAL                     :: LBASICEMIS
     LOGICAL                     :: LSETH2O
     INTEGER                     :: CFCYEAR

     !----------------------------------------
     ! CO MENU fields
     !----------------------------------------
     LOGICAL                     :: LPCO_CH4
     LOGICAL                     :: LPCO_NMVOC

     !----------------------------------------
     ! CO2 MENU fields
     !----------------------------------------
     LOGICAL                     :: LFOSSIL
     LOGICAL                     :: LCHEMCO2
     LOGICAL                     :: LBIODIURNAL
     LOGICAL                     :: LBIONETCLIM
     LOGICAL                     :: LOCEAN
     LOGICAL                     :: LSHIP
     LOGICAL                     :: LPLANE
     LOGICAL                     :: LFFBKGRD
     LOGICAL                     :: LBIOSPHTAG
     LOGICAL                     :: LFOSSILTAG
     LOGICAL                     :: LSHIPTAG
     LOGICAL                     :: LPLANETAG

     !----------------------------------------
     ! CHEMISTRY MENU fields
     !----------------------------------------
     LOGICAL                     :: LCHEM
     LOGICAL                     :: LSCHEM
     LOGICAL                     :: LLINOZ
     LOGICAL                     :: LSYNOZ
     INTEGER                     :: TS_CHEM
     REAL(fp)                    :: GAMMA_HO2
     LOGICAL                     :: LUCX
     LOGICAL                     :: LACTIVEH2O
     LOGICAL                     :: LINITSPEC
     LOGICAL                     :: USE_ONLINE_O3
     LOGICAL                     :: USE_O3_FROM_MET
     LOGICAL                     :: USE_TOMS_O3

     !----------------------------------------
     ! RADIATION MENU fields
     !----------------------------------------
     LOGICAL                     :: LRAD
     LOGICAL                     :: LLWRAD
     LOGICAL                     :: LSWRAD
     LOGICAL,            POINTER :: LSKYRAD(:)
     INTEGER                     :: TS_RAD
     INTEGER                     :: NWVSELECT
     REAL(8),            POINTER :: WVSELECT(:)
     CHARACTER(LEN=5),   POINTER :: STRWVSELECT(:)
     INTEGER                     :: NSPECRADMENU
     INTEGER,            POINTER :: LSPECRADMENU(:)

     !----------------------------------------
     ! TRANSPORT MENU fields
     !----------------------------------------
     LOGICAL                     :: LTRAN
     LOGICAL                     :: LFILL
     INTEGER                     :: TPCORE_IORD
     INTEGER                     :: TPCORE_JORD
     INTEGER                     :: TPCORE_KORD
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
     REAL(fp)                    :: WETD_CONV_SCAL
     LOGICAL                     :: PBL_DRYDEP      

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
     CHARACTER(LEN=255)          :: HistoryInputFile
     INTEGER                     :: ND01,             LD01
     INTEGER                     :: ND02,             LD02
     INTEGER                     :: ND03,             LD03
     INTEGER                     :: ND04,             LD04
     INTEGER                     :: ND05,             LD05
     INTEGER                     :: ND06,             LD06
     INTEGER                     :: ND07,             LD07
     INTEGER                     :: ND08,             LD08
     INTEGER                     :: ND09,             LD09
     INTEGER                     :: ND10,             LD10
     INTEGER                     :: ND11,             LD11
     INTEGER                     :: ND12,             LD12
     INTEGER                     :: ND13,             LD13
     INTEGER                     :: ND14,             LD14
     INTEGER                     :: ND15,             LD15
     INTEGER                     :: ND16,             LD16
     INTEGER                     :: ND17,             LD17
     INTEGER                     :: ND18,             LD18
     INTEGER                     :: ND19,             LD19
     INTEGER                     :: ND21,             LD21
     INTEGER                     :: ND22,             LD22
     INTEGER                     :: ND24,             LD24
     INTEGER                     :: ND25,             LD25
     INTEGER                     :: ND26,             LD26
     INTEGER                     :: ND27,             LD27
     INTEGER                     :: ND28,             LD28
     INTEGER                     :: ND29,             LD29
     INTEGER                     :: ND30,             LD30
     INTEGER                     :: ND31,             LD31
     INTEGER                     :: ND32,             LD32
     INTEGER                     :: ND33,             LD33
     INTEGER                     :: ND34,             LD34
     INTEGER                     :: ND35,             LD35
     INTEGER                     :: ND36,             LD36
     INTEGER                     :: ND37,             LD37
     INTEGER                     :: ND38,             LD38
     INTEGER                     :: ND39,             LD39
     INTEGER                     :: ND41,             LD41
     INTEGER                     :: ND42,             LD42
     INTEGER                     :: ND43,             LD43
     INTEGER                     :: ND44,             LD44
     INTEGER                     :: ND45,             LD45
     INTEGER                     :: ND46,             LD46
     INTEGER                     :: ND47,             LD47
     INTEGER                     :: ND52,             LD52
     INTEGER                     :: ND53,             LD53
     INTEGER                     :: ND54,             LD54
     INTEGER                     :: ND55,             LD55
     INTEGER                     :: ND56,             LD56
     INTEGER                     :: ND57,             LD57
     INTEGER                     :: ND59,             LD59
     INTEGER                     :: ND60,             LD60
     INTEGER                     :: ND61,             LD61
     INTEGER                     :: ND62,             LD62
     INTEGER                     :: ND64,             LD64
     INTEGER                     :: ND66,             LD66
     INTEGER                     :: ND67,             LD67
     INTEGER                     :: ND68,             LD68
     INTEGER                     :: ND69,             LD69
     INTEGER                     :: ND70,             LD70
     INTEGER                     :: ND71,             LD71
     INTEGER                     :: ND72,             LD72
     INTEGER                     :: ND73,             LD73

     INTEGER                     :: TS_DIAG
     LOGICAL                     :: LPRT
     INTEGER,            POINTER :: TINDEX(:,:)
     INTEGER,            POINTER :: TCOUNT(:) 				  
     INTEGER,            POINTER :: TMAX(:)
     LOGICAL                     :: DO_DIAG_WRITE

     ! Collection ids
     INTEGER                     :: DIAG_COLLECTION
     INTEGER                     :: GC_RST_COLLECTION ! Used only for NetCDF

     !----------------------------------------
     ! PLANEFLIGHT MENU fields
     !----------------------------------------
     LOGICAL                     :: DO_PF
     CHARACTER(LEN=255)          :: PF_IFILE
     CHARACTER(LEN=255)          :: PF_OFILE

     !----------------------------------------
     ! OBSPACK MENU fields
     !----------------------------------------
     LOGICAL                     :: Do_ObsPack
     LOGICAL                     :: ObsPack_Quiet
     CHARACTER(LEN=255)          :: ObsPack_InputFile
     CHARACTER(LEN=255)          :: ObsPack_OutputFile
     INTEGER                     :: ObsPack_nSpc
     CHARACTER(LEN=255), POINTER :: ObsPack_SpcName(:) 

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
     INTEGER                     :: N_ND49
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
     INTEGER                     :: N_ND50
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
     INTEGER                     :: N_ND51
     CHARACTER(LEN=255)          :: ND51_FILE
     LOGICAL                     :: LND51_HDF
     INTEGER,            POINTER :: ND51_TRACERS(:)
     REAL(fp)                    :: ND51_HR_WRITE
     REAL(fp)                    :: ND51_HR1
     REAL(fp)                    :: ND51_HR2
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
     INTEGER                     :: N_ND51b
     CHARACTER(LEN=255)          :: ND51b_FILE
     LOGICAL                     :: LND51b_HDF
     INTEGER,            POINTER :: ND51b_TRACERS(:)
     REAL(fp)                    :: ND51b_HR_WRITE
     REAL(fp)                    :: ND51b_HR1
     REAL(fp)                    :: ND51b_HR2
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
     INTEGER                     :: N_ND63
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
     INTEGER                     :: ND65, LD65
     LOGICAL                     :: DO_SAVE_O3
     LOGICAL                     :: DO_SAVE_PCO
     INTEGER                     :: NFAM
     CHARACTER(LEN=255), POINTER :: FAM_NAME(:)
     CHARACTER(LEN=255), POINTER :: FAM_TYPE(:)

     !----------------------------------------
     ! NESTED GRID MENU fields
     !----------------------------------------
     LOGICAL                     :: ITS_A_NESTED_GRID
     LOGICAL                     :: LWINDO
     LOGICAL                     :: LWINDO2x25
     LOGICAL                     :: LWINDO_NA
     CHARACTER(LEN=255)          :: TPBC_DIR_NA
     LOGICAL                     :: LWINDO_EU
     CHARACTER(LEN=255)          :: TPBC_DIR_EU
     LOGICAL                     :: LWINDO_CH
     CHARACTER(LEN=255)          :: TPBC_DIR_CH
     LOGICAL                     :: LWINDO_AS
     CHARACTER(LEN=255)          :: TPBC_DIR_AS
     LOGICAL                     :: LWINDO_CU
     CHARACTER(LEN=255)          :: TPBC_DIR
     INTEGER                     :: NESTED_TS
     INTEGER                     :: NESTED_I1
     INTEGER                     :: NESTED_J1
     INTEGER                     :: NESTED_I2
     INTEGER                     :: NESTED_J2
     INTEGER                     :: NESTED_I0W
     INTEGER                     :: NESTED_J0W
     INTEGER                     :: NESTED_I0E
     INTEGER                     :: NESTED_J0E

     !----------------------------------------
     ! BENCHMARK MENU fields
     !----------------------------------------
     LOGICAL                     :: LSTDRUN
     CHARACTER(LEN=255)          :: STDRUN_INIT_FILE
     CHARACTER(LEN=255)          :: STDRUN_FINAL_FILE

     !----------------------------------------
     ! MERCURY MENU fields
     !----------------------------------------     
     INTEGER                     :: ANTHRO_Hg_YEAR
     CHARACTER(LEN=255)          :: HG_SCENARIO
     LOGICAL                     :: USE_CHECKS
     LOGICAL                     :: LDYNOCEAN
     LOGICAL                     :: LPREINDHG
     LOGICAL                     :: LGTMM
     CHARACTER(LEN=255)          :: GTMM_RST_FILE
     LOGICAL                     :: LARCTICRIV
     LOGICAL                     :: LKRedUV

     !----------------------------------------
     ! CH4 MENU fields
     !----------------------------------------  
     LOGICAL                     :: GOSAT_CH4_OBS
     LOGICAL                     :: TCCON_CH4_OBS

     !----------------------------------------
     ! POPS MENU fields
     !----------------------------------------
     CHARACTER(LEN=3)            :: POP_TYPE
     LOGICAL                     :: CHEM_PROCESS
     REAL(fp)                    :: POP_XMW
     REAL(fp)                    :: POP_KOA
     REAL(fp)                    :: POP_KBC
     REAL(fp)                    :: POP_K_POPG_OH
     REAL(fp)                    :: POP_K_POPP_O3A
     REAL(fp)                    :: POP_K_POPP_O3B
     REAL(fp)                    :: POP_HSTAR
     REAL(fp)                    :: POP_DEL_H
     REAL(fp)                    :: POP_DEL_Hw

     !----------------------------------------
     ! Fields for interface to GEOS-5 GCM
     !----------------------------------------
#if defined( MODEL_GEOS )
     LOGICAL                     :: haveImpRst   = .FALSE.
     LOGICAL                     :: AlwaysSetH2O = .TRUE.
     LOGICAL                     :: UseOnlineVUD = .FALSE.
     INTEGER                     :: LLFASTJX     = 601
     INTEGER                     :: NN_RxnRates             ! # of diagnosed reaction rates
     INTEGER, POINTER            :: RxnRates_IDs(:)         ! Reaction rate numbers to be diagnosed
     INTEGER                     :: NN_RxnRconst            ! # of diagnosed reaction rates
     INTEGER, POINTER            :: RxnRconst_IDs(:)        ! Reaction rate numbers to be diagnosed
     INTEGER                     :: NN_Jvals                ! # of diagnosed Jvalues 
     INTEGER, POINTER            :: Jval_IDs(:)             ! J-values to be diagnosed
     INTEGER                     :: FJX_EXTRAL_ITERMAX = 5
     LOGICAL                     :: FJX_EXTRAL_ERR     = .TRUE.
     LOGICAL                     :: KppStop            = .TRUE. ! Stop KPP if integration fails twice
#else
     LOGICAL                     :: haveImpRst
     LOGICAL                     :: AlwaysSetH2O
#endif

     !----------------------------------------
     ! Fields for LINOZ strat chem
     !----------------------------------------
     INTEGER                     :: LINOZ_NLEVELS
     INTEGER                     :: LINOZ_NLAT
     INTEGER                     :: LINOZ_NMONTHS
     INTEGER                     :: LINOZ_NFIELDS
     REAL(fp),           POINTER :: LINOZ_TPARM(:,:,:,:)

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
!  17 Sep 2013 - M. Sulprizio- Add LDSTUP flag for acid uptake on dust aerosols
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57_DIR to GEOS_FP_DIR
!  03 Oct 2013 - M. Sulprizio- Removed obsolete LMFCT for flux correction
!  03 Oct 2013 - M. Sulprizio- Removed obsolete LAVHRRLAI and LMODISLAI
!  13 Dec 2013 - M. Sulprizio- Add USE_O3_FROM_MET logical flag
!  16 Apr 2014 - M. Sulprizio- Add field for PSC restart file
!  23 Jun 2014 - R. Yantosca - Add POP_EMISDIR field for POPs simlulation
!  25 Jun 2014 - R. Yantosca - Now add Input_Opt%SIM_TYPE field
!  29 Sep 2014 - R. Yantosca - Now add Input_Opt%N_DUST_BINS field
!  03 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!  03 Dec 2014 - M. Sulprizio- Add fields for Radiation Menu
!  16 Dec 2014 - R. Yantosca - Removed JLOP, JLOP_PREV; these are in State_Chm
!  01 Apr 2015 - R. Yantosca - Add extra nested-grid fields
!  09 Apr 2015 - M. Sulprizio- Removed fields for NAPEMISS, POAEMISSSCALE,
!                              and PST_RST_FILE. These options are now handled
!                              by HEMCO.
!  11 Aug 2015 - R. Yantosca - Add MERRA2_DIR field to OptInput
!  26 Jan 2016 - E. Lundgren - Add fields for netcdf diagnostics
!  04 Feb 2016 - C. Keller   - Add LINITSPEC. Used in ESMF to initialize species
!                              concentrations from globchem.dat.
!  04 Feb 2016 - M. Sulprizio- Add Hg_CAT and Hg_CAT_FULL arrays for tagged Hg
!                              simulations
!  27 Apr 2016 - R. Yantosca - Remove Hg_Cat, Hg_Cat_Full fields
!  17 May 2016 - R. Yantosca - Remove TRACER_N_CONST, TRACER_CONST, ID_EMITTED,
!                              TRACER_COEFF
!  31 May 2016 - E. Lundgren - Remove TRACER_MW_KG, TRACER_MW_G, and XNUMOL
!  23 Jun 2016 - R. Yantosca - Remove references to APM code; it is no longer
!                              compatible with the FlexChem implementation
!  13 Jul 2016 - R. Yantosca - Remove some unused drydep fields
!  16 Aug 2016 - M. Sulprizio- Rename from gigc_input_opt_mod.F90 to
!                              input_opt_mod.F90. The "gigc" nomenclature is
!                              no longer used.
!  29 Aug 2016 - M. Sulprizio- Rename N_TRACERS to N_ADVECT and TRACER_NAME to
!                              AdvectSpc_Name to reflect that we now refer to
!                              tracers as advected species
!  20 Sep 2016 - R. Yantosca - LND51_HDF and LND51b_HDF are now declared
!                              as LOGICAL, not INTEGER.  This chokes Gfortran.
!  03 Oct 2016 - R. Yantosca - LWINDO_CU has to be LOGICAL, not INTEGER
!  16 Jun 2017 - M. Sulprizio- Remove switches for CH4 emissions since these
!                              are now handled by HEMCO
!  12 Jul 2017 - R. Yantosca - Add Input_Opt%HistoryInputFile field
!  13 Jul 2017 - E. Lundgren - Add passive species variables
!  24 Aug 2017 - M. Sulprizio- Remove obsolete options: GCAP_DIR, GEOS_4_DIR,
!                              GEOS_5_DIR, MERRA_DIR, TEMP_DIR, LUNZIP, LWAIT
!  13 Sep 2017 - M. Sulprizio- Remove USE_OLSON_2001. Olson 2001 is now the
!                              default.
!  14 Sep 2017 - M. Sulprizio- Add USE_ONLINE_O3 and USE_TOMS_O3 to options for
!                              overhead O3 in chemistry menu
!  02 Nov 2017 - R. Yantosca - Bug fix: LBIOFUEL should be LOGICAL
!  07 Nov 2017 - R. Yantosca - Remove LVARTROP; it's not needed
!  29 Dec 2017 - C. Keller   - Added LLSTRAT. Used in gc_environment_mod.F90
!  29 Dec 2017 - C. Keller   - Added AlwaysSetH2O.
!  04 Apr 2018 - E. Lundgren - Remove MAX_PASV; use # from input.geos instead
!  30 Aug 2018 - C. Keller   - Remove LLSTRAT. Only used in GEOS-5, obtained
!                              from gridded comp module directly.
!  15 Oct 2018 - E. Lundgren - Remove LFUTURECFC; no longer needed with ucx_mod updates
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
! !IROUTINE: Set_Input_Opt
!
! !DESCRIPTION: Subroutine SET\_INPUT\_OPT intializes all GEOS-Chem
!  options carried in Input Options derived type object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_Input_Opt( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
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
!  12 Sep 2013 - M. Sulprizio- Double size of IDDEP to account for dust
!                              alkalinity (tdf 04/10/08)
!  17 Sep 2013 - M. Sulprizio- Add LDSTUP flag for acid uptake on dust aerosols
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57_DIR to GEOS_FP_DIR
!  25 Jun 2014 - R. Yantosca - Now initialize Input_Opt%SIM_TYPE field
!  03 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!  05 Mar 2015 - R. Yantosca - Added RES_DIR, CHEM_INPUTS_DIR fields
!  06 Mar 2015 - R. Yantosca - Now initialize directory names with './'
!  01 Apr 2015 - R. Yantosca - Now initialize extra nested-grid fields
!  04 Mar 2016 - C. Keller   - Added WETD_CONV_SCAL, LSYNOZ, LCAPTROP, and 
!                              OZONOPAUSE. These are only used within ESMF. 
!  17 May 2016 - R. Yantosca - Remove TRACER_N_CONST, TRACER_CONST, ID_EMITTED,
!                              TRACER_COEFF
!  31 May 2016 - E. Lundgren - Remove TRACER_MW_G, TRACER_MW_KG, and XNUMOL
!  13 Jul 2016 - R. Yantosca - Remove some obsolete drydep fields
!  13 Jul 2016 - R. Yantosca - Remove ID_TRACER, NUMDEP
!  16 Mar 2017 - R. Yantosca - Remove obsolete family and drydep variables
!  17 Mar 2017 - R. Yantosca - Remove IDDEP, DUSTREFF, DUSTDEN
!  12 Jul 2017 - R. Yantosca - Initialize Input_Opt%HistoryInputFile field
!  02 Nov 2017 - R. Yantosca - LWINDO_CU should be .FALSE., not 0
!  07 Nov 2017 - R. Yantosca - Remove LVARTROP; it's not needed
!  08 Mar 2018 - R. Yantosca - Bug fix, remove reference to TINDEX here
!  06 Nov 2018 - R. Yantosca - Add error trapping for allocation statements
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=30) :: arrayId

    !----------------------------------------
    ! Initialize
    ! Set pointers to NULL for safety's sake
    !----------------------------------------
    RC                               =  GC_SUCCESS
    Input_Opt%PASSIVE_NAME           => NULL()
    Input_Opt%PASSIVE_ID             => NULL()
    Input_Opt%PASSIVE_MW             => NULL()
    Input_Opt%PASSIVE_TAU            => NULL()
    Input_Opt%PASSIVE_INITCONC       => NULL()
    Input_Opt%PASSIVE_DECAYID        => NULL()
    Input_Opt%AdvectSpc_Name         => NULL()
    Input_Opt%SALA_REDGE_um          => NULL()
    Input_Opt%SALC_REDGE_um          => NULL()
    Input_Opt%LSKYRAD                => NULL()
    Input_Opt%LSPECRADMENU           => NULL()
    Input_Opt%NJDAY                  => NULL()
    Input_Opt%TINDEX                 => NULL()
    Input_Opt%TCOUNT                 => NULL()
    Input_Opt%TMAX                   => NULL()
    Input_Opt%ND48_IARR              => NULL()
    Input_Opt%ND48_JARR              => NULL()
    Input_Opt%ND48_LARR              => NULL()
    Input_Opt%ND48_NARR              => NULL()
    Input_Opt%ND49_TRACERS           => NULL()
    Input_Opt%ND50_TRACERS           => NULL()
    Input_Opt%ND51_TRACERS           => NULL()
    Input_Opt%ND51b_TRACERS          => NULL()
    Input_Opt%ND63_TRACERS           => NULL()
    Input_Opt%FAM_NAME               => NULL()
    Input_Opt%FAM_TYPE               => NULL()
    Input_Opt%LINOZ_TPARM            => NULL()
    
    !----------------------------------------
    ! General Runtime & Distributed Comp Info
    !----------------------------------------
    Input_Opt%NPES                   = 1       ! Assume Serial Sim.
    Input_Opt%HPC                    = .false. ! Assume Serial Sim.
    Input_Opt%myCpu                  = -1
    Input_Opt%RootCPU                = .false.
    
    !----------------------------------------
    ! SIZE PARAMETER fields
    !
    ! Set to large placeholder values
    !----------------------------------------
#if   defined( RRTMG )
    Input_Opt%Max_BPCH_Diag          = 187 ! Mirror MAX_DIAG in CMN_DIAG_mod.F
#else
    Input_Opt%Max_BPCH_Diag          = 80  ! Mirror MAX_DIAG in CMN_DIAG_mod.F
#endif
    Input_Opt%Max_Families           = 250
    Input_Opt%Max_AdvectSpc          = 300
    Input_Opt%Max_PassiveSpc         = 50

    !----------------------------------------
    ! SIMULATION MENU fields 
    !----------------------------------------
    Input_Opt%NYMDb                  = 0
    Input_Opt%NHMSb                  = 0
    Input_Opt%NYMDe                  = 0
    Input_Opt%NHMSe                  = 0
    Input_Opt%SimLengthSec           = 0
    Input_Opt%RUN_DIR                = './'
    Input_Opt%DATA_DIR               = './'
    Input_Opt%RES_DIR                = './'
    Input_Opt%CHEM_INPUTS_DIR        = './'
    Input_Opt%RES_DIR                = './'
    Input_Opt%LCAPTROP               = .FALSE.
    Input_Opt%OZONOPAUSE             = -999.0 
    Input_Opt%NESTED_I0              = 0
    Input_Opt%NESTED_J0              = 0
    Input_Opt%HcoConfigFile          = ''

    !----------------------------------------
    ! ADVECTED SPECIES MENU fields
    !----------------------------------------
    arrayId = 'Input_Opt%AdvectSpc_Name'
    ALLOCATE( Input_Opt%AdvectSpc_Name( Input_Opt%Max_AdvectSpc ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    
    Input_Opt%N_ADVECT               = 0
    Input_Opt%AdvectSpc_Name         = ''
    Input_Opt%SIM_TYPE               = 0
    Input_Opt%SIM_NAME               = ''
    Input_Opt%LSPLIT                 = .FALSE.
    Input_Opt%ITS_A_RnPbBe_SIM       = .FALSE.
    Input_Opt%ITS_A_CH3I_SIM         = .FALSE.
    Input_Opt%ITS_A_FULLCHEM_SIM     = .FALSE.
    Input_Opt%ITS_A_HCN_SIM          = .FALSE.
    Input_Opt%ITS_A_TAGO3_SIM        = .FALSE.
    Input_Opt%ITS_A_TAGCO_SIM        = .FALSE.
    Input_Opt%ITS_A_C2H6_SIM         = .FALSE.
    Input_Opt%ITS_A_CH4_SIM          = .FALSE.
    Input_Opt%ITS_AN_AEROSOL_SIM     = .FALSE.
    Input_Opt%ITS_A_MERCURY_SIM      = .FALSE.
    Input_Opt%ITS_A_CO2_SIM          = .FALSE.
    Input_Opt%ITS_A_H2HD_SIM         = .FALSE.
    Input_Opt%ITS_A_POPS_SIM         = .FALSE.

    !----------------------------------------
    ! PASSIVE SPECIES MENU fields
    !----------------------------------------

    ALLOCATE( Input_Opt%PASSIVE_NAME    ( Input_Opt%Max_PassiveSpc ), STAT=RC )
    ALLOCATE( Input_Opt%PASSIVE_LONGNAME( Input_Opt%Max_PassiveSpc ), STAT=RC )
    ALLOCATE( Input_Opt%PASSIVE_ID      ( Input_Opt%Max_PassiveSpc ), STAT=RC )
    ALLOCATE( Input_Opt%PASSIVE_MW      ( Input_Opt%Max_PassiveSpc ), STAT=RC )
    ALLOCATE( Input_Opt%PASSIVE_TAU     ( Input_Opt%Max_PassiveSpc ), STAT=RC )
    ALLOCATE( Input_Opt%PASSIVE_INITCONC( Input_Opt%Max_PassiveSpc ), STAT=RC )
    ALLOCATE( Input_Opt%PASSIVE_DECAYID ( Input_Opt%Max_PassiveSpc ), STAT=RC )

    Input_Opt%NPASSIVE               = 0 
    Input_Opt%NPASSIVE_DECAY         = 0 
    Input_Opt%PASSIVE_NAME           = ''
    Input_Opt%PASSIVE_LONGNAME       = ''
    Input_Opt%PASSIVE_ID             = 0
    Input_Opt%PASSIVE_MW             = 0e+0_fp
    Input_Opt%PASSIVE_TAU            = 0e+0_fp
    Input_Opt%PASSIVE_INITCONC       = 0e+0_fp
    Input_Opt%PASSIVE_DECAYID        = 0

    !----------------------------------------
    ! AEROSOL MENU fields
    !----------------------------------------
    arrayId = 'Input_Opt%SALA_REDGE_um'
    ALLOCATE( Input_Opt%SALA_REDGE_um( 2 ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    arrayId = 'Input_Opt%SALC_REDGE_um'
    ALLOCATE( Input_Opt%SALC_REDGE_um( 2 ), STAT=RC )     
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    
    Input_Opt%LSULF                  = .FALSE.
    Input_Opt%LMETALCATSO2           = .FALSE.
    Input_Opt%LCARB                  = .FALSE.
    Input_Opt%LBRC                   = .FALSE.
    Input_Opt%LSOA                   = .FALSE.
    Input_Opt%LMPOA                  = .FALSE.
    Input_Opt%LSVPOA                 = .FALSE.
    Input_Opt%LOMOC                  = .FALSE.
    Input_Opt%LDUST                  = .FALSE.
    Input_Opt%LDEAD                  = .FALSE.
    Input_Opt%LDSTUP                 = .FALSE.
    Input_Opt%LSSALT                 = .FALSE.
    Input_Opt%SALA_REDGE_um          = 0e+0_fp
    Input_Opt%SALC_REDGE_um          = 0e+0_fp
    Input_Opt%LGRAVSTRAT             = .FALSE.
    Input_Opt%LSOLIDPSC              = .FALSE.
    Input_Opt%LHOMNUCNAT             = .FALSE.
    Input_Opt%T_NAT_SUPERCOOL        = 0e+0_fp
    Input_Opt%P_ICE_SUPERSAT         = 0e+0_fp
    Input_Opt%LPSCCHEM               = .FALSE.
    Input_Opt%LSTRATOD               = .FALSE.

    !----------------------------------------
    ! EMISSIONS MENU fields
    !----------------------------------------
    Input_Opt%LEMIS                  = .FALSE.
    Input_Opt%TS_EMIS                = 0
    Input_Opt%LSOILNOX               = .FALSE.
    Input_Opt%LWARWICK_VSLS          = .FALSE.
    Input_Opt%LSSABr2                = .FALSE.
    Input_Opt%LFIX_PBL_BRO           = .FALSE.
    Input_Opt%LCH4EMIS               = .FALSE.
    Input_Opt%LCH4SBC                = .FALSE.
    Input_Opt%LOCSEMIS               = .FALSE.
    Input_Opt%LCFCEMIS               = .FALSE.
    Input_Opt%LCLEMIS                = .FALSE.
    Input_Opt%LBREMIS                = .FALSE.
    Input_Opt%LN2OEMIS               = .FALSE.
    Input_Opt%LBASICEMIS             = .FALSE.
    Input_Opt%LSETH2O                = .FALSE.
    Input_Opt%CFCYEAR                = 0

    !----------------------------------------
    ! CO MENU fields
    !----------------------------------------
    Input_Opt%LPCO_CH4               = .FALSE.
    Input_Opt%LPCO_NMVOC             = .FALSE.

    !----------------------------------------
    ! CO2 MENU fields
    !----------------------------------------
    Input_Opt%LFOSSIL                = .FALSE.
    Input_Opt%LCHEMCO2               = .FALSE.
    Input_Opt%LBIOFUEL               = .FALSE.
    Input_Opt%LBIODIURNAL            = .FALSE.
    Input_Opt%LBIONETCLIM            = .FALSE.
    Input_Opt%LOCEAN                 = .FALSE.
    Input_Opt%LSHIP                  = .FALSE.
    Input_Opt%LPLANE                 = .FALSE.
    Input_Opt%LFFBKGRD               = .FALSE.
    Input_Opt%LBIOSPHTAG             = .FALSE.
    Input_Opt%LFOSSILTAG             = .FALSE.
    Input_Opt%LSHIPTAG               = .FALSE.
    Input_Opt%LPLANETAG              = .FALSE.

    !----------------------------------------
    ! CHEMISTRY MENU fields
    !----------------------------------------
    Input_Opt%LCHEM                  = .FALSE.
    Input_Opt%LSCHEM                 = .FALSE.
    Input_Opt%LLINOZ                 = .FALSE. 
    Input_Opt%LSYNOZ                 = .FALSE. 
    Input_Opt%TS_CHEM                = 0
    Input_Opt%GAMMA_HO2              = 0e+0_fp
    Input_Opt%LUCX                   = .FALSE.
    Input_Opt%LACTIVEH2O             = .FALSE.
    Input_Opt%LINITSPEC              = .FALSE.
    Input_Opt%USE_ONLINE_O3          = .FALSE.
    Input_Opt%USE_O3_FROM_MET        = .FALSE.
    Input_Opt%USE_TOMS_O3            = .FALSE.

    !----------------------------------------
    ! RADIATION MENU fields (for RRTMG only)
    !----------------------------------------
    arrayId = 'Input_Opt%LSKYRAD'
    ALLOCATE( Input_Opt%LSKYRAD( 2 ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    arrayId = 'Input_Opt%WVSELECT'
    ALLOCATE( Input_Opt%WVSELECT( 3 ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    arrayId = 'Input_Opt%STRWVSELECT'
    ALLOCATE( Input_Opt%STRWVSELECT( 3 ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Number of RRTMG outputs (change as necessary)
    Input_Opt%NSpecRadMenu           = 11

    arrayId = 'Input_Opt%LSPECRADMENU'
    ALLOCATE( Input_Opt%LSPECRADMENU( Input_Opt%NSpecRadMenu ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Input_Opt%LSpecRadMenu           = 0

    Input_Opt%LRAD                   = .FALSE.
    Input_Opt%LLWRAD                 = .FALSE.
    Input_Opt%LSWRAD                 = .FALSE.
    Input_Opt%LSKYRAD                = .FALSE.
    Input_Opt%TS_RAD                 = 0
    Input_Opt%NWVSELECT              = 0
    Input_Opt%WVSELECT               = 0e+0_fp
    Input_Opt%STRWVSELECT            = ''

    !----------------------------------------
    ! TRANSPORT MENU fields
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
    Input_Opt%WETD_CONV_SCAL         = 1.0_fp 
    Input_Opt%PBL_DRYDEP             = .FALSE.

    !----------------------------------------
    ! GAMAP_MENU fields
    !----------------------------------------
    Input_Opt%GAMAP_DIAGINFO         = ''
    Input_Opt%GAMAP_TRACERINFO       = ''

    !----------------------------------------
    ! OUTPUT MENU fields
    !----------------------------------------
    arrayId = 'Input_Opt%NJDAY'
    ALLOCATE( Input_Opt%NJDAY( 366 ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    Input_Opt%NJDAY                  = 0

    !----------------------------------------
    ! DIAGNOSTIC MENU fields
    !----------------------------------------
    Input_Opt%HistoryInputFile       = ''
    Input_Opt%DIAG_COLLECTION        = -999
    Input_Opt%TS_DIAG                = 0
    ALLOCATE( Input_Opt%TCOUNT( Input_Opt%Max_BPCH_Diag ), STAT=RC )
    ALLOCATE( Input_Opt%TMAX  ( Input_Opt%Max_BPCH_Diag ), STAT=RC )

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
    Input_Opt%ND21                   = 0
    Input_Opt%ND22                   = 0
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
    Input_Opt%ND41                   = 0
    Input_Opt%ND42                   = 0
    Input_Opt%ND43                   = 0
    Input_Opt%ND44                   = 0
    Input_Opt%ND45                   = 0
    Input_Opt%ND46                   = 0
    Input_Opt%ND47                   = 0
    Input_Opt%ND52                   = 0
    Input_Opt%ND53                   = 0
    Input_Opt%ND54                   = 0
    Input_Opt%ND55                   = 0
    Input_Opt%ND56                   = 0
    Input_Opt%ND57                   = 0
    Input_Opt%ND59                   = 0
    Input_Opt%ND60                   = 0
    Input_Opt%ND61                   = 0
    Input_Opt%ND62                   = 0
    Input_Opt%ND64                   = 0
    Input_Opt%ND65                   = 0
    Input_Opt%ND66                   = 0
    Input_Opt%ND67                   = 0
    Input_Opt%ND68                   = 0
    Input_Opt%ND69                   = 0
    Input_Opt%ND70                   = 0
    Input_Opt%ND71                   = 0
    Input_Opt%ND72                   = 0
    Input_Opt%ND73                   = 0
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
    Input_Opt%LD21                   = 0
    Input_Opt%LD22                   = 0
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
    Input_Opt%LD41                   = 0
    Input_Opt%LD42                   = 0
    Input_Opt%LD43                   = 0
    Input_Opt%LD44                   = 0
    Input_Opt%LD45                   = 0
    Input_Opt%LD46                   = 0
    Input_Opt%LD47                   = 0
    Input_Opt%LD52                   = 0
    Input_Opt%LD53                   = 0
    Input_Opt%LD54                   = 0
    Input_Opt%LD55                   = 0
    Input_Opt%LD56                   = 0
    Input_Opt%LD57                   = 0
    Input_Opt%LD59                   = 0
    Input_Opt%LD60                   = 0
    Input_Opt%LD61                   = 0
    Input_Opt%LD62                   = 0
    Input_Opt%LD64                   = 0
    Input_Opt%LD65                   = 0
    Input_Opt%LD66                   = 0
    Input_Opt%LD67                   = 0
    Input_Opt%LD68                   = 0
    Input_Opt%LD69                   = 0
    Input_Opt%LD70                   = 0
    Input_Opt%LD71                   = 0
    Input_Opt%LD72                   = 0
    Input_Opt%LD73                   = 0
    Input_Opt%LPRT                   = .FALSE.
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
    ! PLANEFLIGHT MENU fields
    !----------------------------------------
    ALLOCATE( Input_Opt%ObsPack_SpcName( 1000 ), STAT=RC )

    Input_Opt%Do_ObsPack             = .FALSE.
    Input_Opt%ObsPack_Quiet          = .FALSE.
    Input_Opt%ObsPack_InputFile      = ''
    Input_Opt%ObsPack_OutputFile     = ''
    Input_Opt%ObsPack_nSpc           = 0
    Input_Opt%ObsPack_SpcName        = ''

    !----------------------------------------
    ! ND48 MENU fields
    !----------------------------------------
    arrayId = 'Input_Opt%ND48_IARR'
    ALLOCATE( Input_Opt%ND48_IARR( 1000 ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    arrayId = 'Input_Opt%ND48_JARR'
    ALLOCATE( Input_Opt%ND48_JARR( 1000 ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    arrayId = 'Input_Opt%ND48_LARR'
    ALLOCATE( Input_Opt%ND48_LARR( 1000 ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    arrayId = 'Input_Opt%ND48_NARR'
    ALLOCATE( Input_Opt%ND48_NARR( 1000 ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

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
    Input_Opt%DO_ND49                = .FALSE.
    Input_Opt%N_ND49                 = 0
    Input_Opt%ND49_FILE              = ''
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
    Input_Opt%DO_ND50                = .FALSE.
    Input_Opt%N_ND50                 = 0
    Input_Opt%ND50_FILE              = ''
    Input_Opt%LND50_HDF              = .FALSE.
    Input_Opt%ND50_IMIN              = 0
    Input_Opt%ND50_IMAX              = 0
    Input_Opt%ND50_JMIN              = 0
    Input_Opt%ND50_JMAX              = 0
    Input_Opt%ND50_LMIN              = 0
    Input_Opt%ND50_LMAX              = 0

    !----------------------------------------
    ! ND51 MENU fields
    !----------------------------------------
    Input_Opt%DO_ND51                = .FALSE.
    Input_Opt%N_ND51                 = 0
    Input_Opt%ND51_FILE              = ''
    Input_Opt%LND51_HDF              = .FALSE.
    Input_Opt%ND51_HR_WRITE          = 0e+0_fp
    Input_Opt%ND51_HR1               = 0e+0_fp
    Input_Opt%ND51_HR2               = 0e+0_fp
    Input_Opt%ND51_IMIN              = 0
    Input_Opt%ND51_IMAX              = 0
    Input_Opt%ND51_JMIN              = 0
    Input_Opt%ND51_JMAX              = 0
    Input_Opt%ND51_LMIN              = 0

    !----------------------------------------
    ! ND51b MENU fields
    !----------------------------------------
    Input_Opt%DO_ND51b               = .FALSE.
    Input_Opt%N_ND51b                = 0
    Input_Opt%ND51b_FILE             = ''
    Input_Opt%LND51b_HDF             = .FALSE.
    Input_Opt%ND51b_HR_WRITE         = 0e+0_fp
    Input_Opt%ND51b_HR1              = 0e+0_fp
    Input_Opt%ND51b_HR2              = 0e+0_fp
    Input_Opt%ND51b_IMIN             = 0
    Input_Opt%ND51b_IMAX             = 0
    Input_Opt%ND51b_JMIN             = 0
    Input_Opt%ND51b_JMAX             = 0
    Input_Opt%ND51b_LMIN             = 0

    !----------------------------------------
    ! ND63 MENU fields
    !----------------------------------------
    Input_Opt%DO_ND63                = .FALSE.
    Input_Opt%N_ND63                 = 0
    Input_Opt%ND63_FILE              = ''
    Input_Opt%ND63_FREQ              = 0
    Input_Opt%ND63_IMIN              = 0
    Input_Opt%ND63_IMAX              = 0
    Input_Opt%ND63_JMIN              = 0
    Input_Opt%ND63_JMAX              = 0

    !----------------------------------------
    ! PROD LOSS MENU fields
    !---------------------------------------
    arrayId = 'Input_Opt%FAM_NAME'
    ALLOCATE( Input_Opt%FAM_NAME( Input_Opt%Max_Families ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    arrayId = 'Input_Opt%FAM_TYPE'    
    ALLOCATE( Input_Opt%FAM_TYPE( Input_Opt%Max_Families ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    
    Input_Opt%DO_SAVE_PL             = .FALSE.
    Input_Opt%ND65                   = 0
    Input_Opt%DO_SAVE_O3             = .FALSE.
    Input_Opt%DO_SAVE_PCO            = .FALSE.
    Input_Opt%NFAM                   = 0
    Input_Opt%FAM_NAME               = ''
    Input_Opt%FAM_TYPE               = ''

    !----------------------------------------
    ! NESTED GRID MENU fields
    !----------------------------------------
    Input_Opt%ITS_A_NESTED_GRID      = .FALSE.
    Input_Opt%LWINDO                 = .FALSE.
    Input_Opt%LWINDO2x25             = .FALSE.
    Input_Opt%LWINDO_NA              = .FALSE.
    Input_Opt%TPBC_DIR_NA            = './'
    Input_Opt%LWINDO_EU              = .FALSE.
    Input_Opt%TPBC_DIR_EU            = './'
    Input_Opt%LWINDO_CH              = .FALSE.
    Input_Opt%TPBC_DIR_CH            = './'
    Input_Opt%LWINDO_AS              = .FALSE.
    Input_Opt%TPBC_DIR_AS            = './'
    Input_Opt%LWINDO_CU              = .FALSE.
    Input_Opt%TPBC_DIR               = './'
    Input_Opt%NESTED_TS              = 0
    Input_Opt%NESTED_I1              = 0
    Input_Opt%NESTED_J1              = 0
    Input_Opt%NESTED_I2              = 0
    Input_Opt%NESTED_J2              = 0
    Input_Opt%NESTED_I0W             = 0
    Input_Opt%NESTED_J0W             = 0 
    Input_Opt%NESTED_I0E             = 0
    Input_Opt%NESTED_J0E             = 0 

    !----------------------------------------
    ! BENCHMARK MENU fields
    !----------------------------------------
    Input_Opt%LSTDRUN                = .FALSE.
    Input_Opt%STDRUN_INIT_FILE       = ''
    Input_Opt%STDRUN_FINAL_FILE      =''

    !----------------------------------------
    ! MERCURY MENU fields
    !----------------------------------------     
    Input_Opt%ANTHRO_Hg_YEAR         = 0
    Input_Opt%HG_SCENARIO            = ''
    Input_Opt%USE_CHECKS             = .FALSE.
    Input_Opt%LDYNOCEAN              = .FALSE.
    Input_Opt%LPREINDHG              = .FALSE.
    Input_Opt%LGTMM                  = .FALSE.
    Input_Opt%GTMM_RST_FILE          = ''

    !----------------------------------------
    ! CH4 MENU fields
    !----------------------------------------  
    Input_Opt%GOSAT_CH4_OBS          = .FALSE.
    Input_Opt%TCCON_CH4_OBS          = .FALSE.

    !----------------------------------------
    ! POPS MENU fields
    !----------------------------------------
    Input_Opt%POP_TYPE               = ''
    Input_Opt%CHEM_PROCESS           = .FALSE.
    Input_Opt%POP_XMW                = 0e+0_fp
    Input_Opt%POP_KOA                = 0e+0_fp
    Input_Opt%POP_KBC                = 0e+0_fp
    Input_Opt%POP_K_POPG_OH          = 0e+0_fp
    Input_Opt%POP_K_POPP_O3A         = 0e+0_fp
    Input_Opt%POP_K_POPP_O3B         = 0e+0_fp
    Input_Opt%POP_HSTAR              = 0e+0_fp
    Input_Opt%POP_DEL_H              = 0e+0_fp
    Input_Opt%POP_DEL_Hw             = 0e+0_fp

    !----------------------------------------
    ! Fields for interface to GEOS-5 GCM
    !----------------------------------------
#if defined( MODEL_GEOS )
!    Input_Opt%haveImpRst             = .FALSE.
!    Input_Opt%AlwaysSetH2O           = .FALSE.
!    Input_Opt%LLFASTJX               = -999
    Input_Opt%NN_RxnRates            = -999
    Input_Opt%RxnRates_IDs           => NULL()
    Input_Opt%NN_RxnRconst           = -999
    Input_Opt%RxnRconst_IDs          => NULL()
    Input_Opt%NN_Jvals               = -999
    Input_Opt%Jval_IDs               => NULL()
#else
    Input_Opt%haveImpRst             = .FALSE.
    Input_Opt%AlwaysSetH2O           = .FALSE.
#endif

    !----------------------------------------
    ! Fields for LINOZ strat chem
    !----------------------------------------
    Input_Opt%LINOZ_NLEVELS          = 25
    Input_Opt%LINOZ_NLAT             = 18
    Input_Opt%LINOZ_NMONTHS          = 12
    Input_Opt%LINOZ_NFIELDS          = 7

    arrayId = 'Input_Opt%LINOZ_TPARM'
    ALLOCATE( Input_Opt%LINOZ_TPARM( Input_Opt%LINOZ_NLEVELS,            &
                                     Input_Opt%LINOZ_NLAT,               &
                                     Input_Opt%LINOZ_NMONTHS,            &
                                     Input_Opt%LINOZ_NFIELDS ), STAT=RC )
    CALL GC_CheckVar( arrayId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    Input_Opt%LINOZ_TPARM            = 0e+0_fp

  END SUBROUTINE Set_Input_Opt
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Input_Opt_Advect
!
! !DESCRIPTION: Subroutine SET\_INPUT\_OPT\_ADVECT intializes all GEOS-Chem
!  options carried in Input Options derived type object that depend on
!  the number of advected species (Input\_Opt%N_ADVECT).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_Input_Opt_Advect( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
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
!  NOTE: These arrays are all for bpch diagnostics, and will eventually 
!  be removed from GEOS-Chem.

! !REVISION HISTORY: 
!  26 Jan 2018 - M. Sulprizio- Initial version
!  04 Apr 2018 - E. Lundgren - Renamed from Set_Input_Opt_Extra to 
!                              Set_Input_Opt_Advect
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Initialize
    RC = GC_SUCCESS

    !=======================================================================
    ! Allocate arrays
    !=======================================================================

    ALLOCATE( Input_Opt%TINDEX(Input_Opt%Max_BPCH_Diag,Input_Opt%N_ADVECT), &
              STAT=RC )
    CALL GC_CheckVar( 'Input_Opt%TINDEX', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Input_Opt%TINDEX = 0

    ALLOCATE( Input_Opt%ND49_TRACERS(Input_Opt%N_ADVECT+Input_Opt%Max_BPCH_Diag),&
              STAT=RC )
    CALL GC_CheckVar( 'Input_Opt%ND49_TRACERS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Input_Opt%ND49_TRACERS = 0

    ALLOCATE( Input_Opt%ND50_TRACERS (Input_Opt%N_ADVECT+Input_Opt%Max_BPCH_Diag),&
              STAT=RC )
    CALL GC_CheckVar( 'Input_Opt%ND50_TRACERS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Input_Opt%ND50_TRACERS = 0

    ALLOCATE( Input_Opt%ND51_TRACERS (Input_Opt%N_ADVECT+Input_Opt%Max_BPCH_Diag),&
              STAT=RC )
    CALL GC_CheckVar( 'Input_Opt%ND51_TRACERS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Input_Opt%ND51_TRACERS = 0

    ALLOCATE( Input_Opt%ND51b_TRACERS(Input_Opt%N_ADVECT+Input_Opt%Max_BPCH_Diag),&
              STAT=RC )
    CALL GC_CheckVar( 'Input_Opt%ND51b_TRACERS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Input_Opt%ND51b_TRACERS = 0

    ALLOCATE( Input_Opt%ND63_TRACERS (Input_Opt%N_ADVECT+Input_Opt%Max_BPCH_Diag),&
              STAT=RC )
    CALL GC_CheckVar( 'Input_Opt%ND63_TRACERS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Input_Opt%ND63_TRACERS = 0

  END SUBROUTINE Set_Input_Opt_Advect
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Input_Opt
!
! !DESCRIPTION: Subroutine CLEANUP\_INPUT\_OPT deallocates all 
!  allocatable fields of the Input Options object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Input_Opt( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
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
!  17 May 2016 - R. Yantosca - Remove TRACER_N_CONST, TRACER_CONST, ID_EMITTED,
!                              TRACER_COEFF
!  31 May 2016 - E. Lundgren - Remove TRACER_MW_G, TRACER_MW_KG, and XNUMOL
!  13 Jul 2016 - R. Yantosca - Remove ID_TRACER
!  16 Mar 2017 - R. Yantosca - Remove obsolete family & drydep fields
!  17 Mar 2017 - R. Yantosca - Remove IDDEP, DUSTREFF, DUSTDEN
!  06 Nov 2018 - R. Yantosca - Now trap errors at deallocation
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Assume success
    RC = GC_SUCCESS

    !======================================================================
    ! Deallocate fields of the Input Options object
    !======================================================================
    IF ( ASSOCIATED( Input_Opt%PASSIVE_NAME ) ) THEN
       DEALLOCATE( Input_Opt%PASSIVE_NAME, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%PASSIVE_NAME', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%PASSIVE_NAME => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%PASSIVE_LONGNAME ) ) THEN
       DEALLOCATE( Input_Opt%PASSIVE_LONGNAME )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%PASSIVE_ID ) ) THEN
       DEALLOCATE( Input_Opt%PASSIVE_ID, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%PASSIVE_ID', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%PASSIVE_ID => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%PASSIVE_MW ) ) THEN
       DEALLOCATE( Input_Opt%PASSIVE_MW, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%PASSIVE_MW', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%PASSIVE_MW => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%PASSIVE_TAU ) ) THEN
       DEALLOCATE( Input_Opt%PASSIVE_TAU, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%PASSIVE_TAU', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%PASSIVE_TAU => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%PASSIVE_INITCONC ) ) THEN
       DEALLOCATE( Input_Opt%PASSIVE_INITCONC, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%PASSIVE_INITCONC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%PASSIVE_INITCONC => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%PASSIVE_DECAYID ) ) THEN
       DEALLOCATE( Input_Opt%PASSIVE_DECAYID )
    ENDIF

    IF ( ASSOCIATED( Input_Opt%AdvectSpc_Name ) ) THEN
       DEALLOCATE( Input_Opt%AdvectSpc_Name, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%AdvectSpcName', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%AdvectSpc_Name => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%SALA_REDGE_um ) ) THEN
       DEALLOCATE( Input_Opt%SALA_REDGE_um, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%SALA_REDGE_um', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%SALA_REDGE_um => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%SALC_REDGE_um ) ) THEN 
       DEALLOCATE( Input_Opt%SALC_REDGE_um, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%SALC_REDGE_um', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%SALC_REDGE_um => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%NJDAY ) ) THEN
       DEALLOCATE( Input_Opt%NJDAY, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%NJDAY', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%NJDAY => NULL()
    ENDIF
    
    IF ( ASSOCIATED( Input_Opt%TINDEX ) ) THEN
       DEALLOCATE( Input_Opt%TINDEX, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%TINDEX', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%TINDEX => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%TCOUNT ) ) THEN
       DEALLOCATE( Input_Opt%TCOUNT, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%TCOUNT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%TCOUNT => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%TMAX ) ) THEN
       DEALLOCATE( Input_Opt%TMAX, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%TMAX', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%TMAX => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%ND48_IARR ) ) THEN
       DEALLOCATE( Input_Opt%ND48_IARR, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%ND48_IARR', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%ND48_IARR => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%ND48_JARR ) ) THEN
       DEALLOCATE( Input_Opt%ND48_JARR, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%ND48_JARR', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%ND48_JARR => NULL()
    ENDIF
     
    IF ( ASSOCIATED( Input_Opt%ND48_LARR ) ) THEN
       DEALLOCATE( Input_Opt%ND48_LARR, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%ND48_LARR', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%ND48_LARR => NULL()
    ENDIF
    
    IF ( ASSOCIATED( Input_Opt%ND48_NARR ) ) THEN
       DEALLOCATE( Input_Opt%ND48_NARR, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%ND48_NARR', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%ND48_NARR => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%ND49_TRACERS ) ) THEN
       DEALLOCATE( Input_Opt%ND49_TRACERS, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%ND49_TRACERS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%ND49_TRACERS => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%ND50_TRACERS ) ) THEN
       DEALLOCATE( Input_Opt%ND50_TRACERS, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%ND50_TRACERS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%ND50_TRACERS => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%ND51_TRACERS ) ) THEN
       DEALLOCATE( Input_Opt%ND51_TRACERS, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%ND51_TRACERS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%ND51_TRACERS => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%ND51b_TRACERS ) ) THEN
       DEALLOCATE( Input_Opt%ND51b_TRACERS, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%ND51b_TRACERS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%ND51b_TRACERS => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%ND63_TRACERS ) ) THEN
       DEALLOCATE( Input_Opt%ND63_TRACERS, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%ND63_TRACERS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%ND63_TRACERS => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%FAM_NAME ) ) THEN
       DEALLOCATE( Input_Opt%FAM_NAME, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%FAM_NAME', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%FAM_NAME => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%LINOZ_TPARM ) ) THEN
       DEALLOCATE( Input_Opt%LINOZ_TPARM, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%LINOZ_TPARM', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%LINOZ_TPARM => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%LSPECRADMENU ) ) THEN
       DEALLOCATE( Input_Opt%LSPECRADMENU, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%LSPECRADMENU', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%LSPECRADMENU => NULL()
    ENDIF

#if defined( MODEL_GEOS )
    !=======================================================================
    ! These fields of Input_Opt are only finalized when
    ! GEOS-Chem is coupled to the online NASA/GEOS ESM
    !=======================================================================
    IF ( ASSOCIATED( Input_Opt%RxnRconst_IDs ) ) THEN
       DEALLOCATE( Input_Opt%RxnRconst_IDs, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%RxnRconst_IDs', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%RxnRconst_IDs => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%RxnRates_IDs ) ) THEN
       DEALLOCATE( Input_Opt%RxnRates_IDs, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%RxnRates_IDs', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%RxnRates_IDs => NULL()
    ENDIF

    IF ( ASSOCIATED( Input_Opt%Jval_IDs ) ) THEN
       DEALLOCATE( Input_Opt%Jval_IDs, STAT=RC )
       CALL GC_CheckVar( 'Input_Opt%Jval_Ids', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Input_Opt%Jval_Ids => NULL()
    ENDIF
#endif

  END SUBROUTINE Cleanup_Input_Opt
!EOC
END MODULE Input_Opt_Mod
