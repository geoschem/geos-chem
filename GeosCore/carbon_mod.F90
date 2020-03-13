!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: carbon_mod.F90
!
! !DESCRIPTION: Module CARBON\_MOD contains arrays and routines for performing
!  a carbonaceous aerosol simulation.  Original code taken from Mian Chin's
!  GOCART model and modified accordingly. (rjp, bmy, 4/2/04, 6/30/10)
!\\
!\\
! !INTERFACE:
!
MODULE CARBON_MOD
!
! !USES:
!
  USE AEROSOL_MOD, ONLY : OCFPOA, OCFOPOA
  USE HCO_ERROR_MOD     ! For HEMCO error reporting
  USE PhysConstants     ! Physical constants
  USE PRECISION_MOD     ! For GEOS-Chem Precisions

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: CHEMCARBON
  PUBLIC  :: EMISSCARBON
  PUBLIC  :: CLEANUP_CARBON
  PUBLIC  :: INIT_CARBON
#ifdef TOMAS
  PUBLIC  :: EMISSCARBONTOMAS
#endif
!
! !PUBLIC DATA MEMBERS:
!
  ! SOAupdate: for branching ratio diagnostic (hotp 5/24/10)
  PUBLIC :: BETANOSAVE

  ! ORVC_SESQ needs to be public so that it can be added as a
  ! restart variable to GEOS-5. (ckeller, 10/10/17)
  PUBLIC :: ORVC_SESQ
!
! !REMARKS:
!  4 Aerosol species : Organic and Black carbon
!                    : hydrophilic (soluble) and hydrophobic of each
!                                                                             .
!  For secondary organic aerosol (SOA) simulation orginal code developed
!  by Chung and Seinfeld [2002] and Hong Liao from John Seinfeld's group
!  at Caltech was taken and further modified accordingly (rjp, bmy, 7/15/04)
!                                                                             .
!  SOAupdate: Traditional SOA simulation updated by hotp 7/2010
!    New code treats semivolatile or nonvolatile POA, aerosol from IVOCs,
!      and has updated biogenic SOA
!    For more details on the updated SOA/POA simulation, see comments
!      in SOA_CHEMISTRY, Pye and Seinfeld ACP 2010, Pye et al. in prep
!      for ACP 2010
!    Note that modifications were made throughout the code for SOAupdate
!                                                                             .
!  References:
!  ============================================================================
!  (1 ) Bond, T.C., E. Bhardwaj, R. Dong, R. Jogani, S. Jung, C. Roden, D.G.
!        Streets, and N.M. Trautmann, "Historical emissions of black and
!        organic carbon aerosol from energy-related combustion, 1850-2000",
!        Global Biogeochem. Cycles, 21, GB2018, doi:10.1029/2006GB002840, 2007.
!  (2 ) Chan, A.W.H., K.E. Kautzman, P.S. Chhabra, J.D. Surratt, M.N. Chan,
!        J.D. Crounse, A. Kurten, P.O. Wennberg, R.C. Flagan, and J.H.
!        Seinfeld, "Secondary orgainc aerosol formation from photooxidation of
!        naphthlene and alkylnaphthalenes: implications for oxidation of
!        intermediate volatility orgainc compounds (IVOCs)", Atmos. Chem. Phys,
!        Vol 9, 3049-3060, doi:10.5194/acp-9-3049-2009, 2009.
!  (3 ) Chung, S.H., and J.H. Seinfeld. "Global distribution and climate
!        forcing of carbonaceous aerosols", J. Geophys. Res., Vol 107(D19),
!        4407, doi:10.1029/2001JD001397, 2002.
!  (4 ) Grieshop, A.P., J.M. Logue, N.M. Donahue, and A.L. Robinson,
!        "Laboratory investigation of photochemical oxidation of organic
!        aerosol deom wood fires 1: Measurement and simulation of organic
!        aerosol evolution", Atmos. Chem. Phys., Vol 9, 1263-1277,
!        doi:10.5194/acp-9-1263-2009, 2009.
!  (5 ) Griffin, R.J., D.R. Cocker, R.C. Flagan, and J.H. Seinfeld, "Orgainc
!        aerosol formation from the oxidation of biogenic hydrocarbons", J.
!        Geophys. Res., 104(D3), 3555-3567, 1999.
!  (6 ) Henze, D.K., and J.H. Seinfeld, "Global secondary organic aerosol from
!        isoprene oxidation", Geophys. Res. Lett., Vol 33, L09812,
!        doi:10.1029/2006GL025976, 2006.
!  (7 ) Henze, D.K., J.H. Seinfeld, N.L. Ng, J.H. Kroll, T.-M. Fu, D.J. Jacob,
!        and C.L. Heald, "Global modeling of secondary orgainc aerosol
!        formation from aromatic hydrocarbons: high vs. low-yield pathways",
!        Atmos. Chem. Phys., Vol 8, 2405-2420, doi:10.5194/acp-8-2405-2008,
!        2008.
!  (8 ) Kroll, J.H., N.L. Ng, S.M. Murphy, R.C. Flagan, and J.H. Seinfeld,
!        "Secondary orgainc aerosol formation from isoprene photooxidation",
!        Environ. Sci. Technol, Vol 40, 1869-1877, doi:10.1021/Es0524301, 2006.
!  (9 ) Liao, H., D.K. Henze, J.H. Seinfeld, S.L Wu, and L.J. Mickley,
!        "Biogenic secondary aerosol over the United States: Comparison of
!        climatological simulations with observations, J. Geophys. Res. Vol
!        112, D06201, doi:10.1029/2006JD007813, 2007.
!  (10) Ng, N.L., P.S. Chhabra, A.W.H. Chan, J.D. Surratt, J.H. Kroll, A.J.
!        Kwan, D.C. McCabe, P.O. Wennberg, A. Sorooshian, S.M. Murphy, N.F.
!        Dalleska, R.C. Flagan, and J.H. Seinfeld, "Effect of NOx level on
!        secondary orgainc aerosol (SOA) formation from the photooxidation of
!        terpenes", Atmos. Chem. Phys., Vol 7, 5159-5174,
!        doi:10.5194/acp-7-5195-2007, 2007a.
!  (11) Ng, N.L., J.H. Kroll, A.W.H. Chan, P.S. Chhabra, R.C. Flagan, and J.H.
!        Seinfeld, "Secondary orgainc aerosol formation from m-xylene, toluele,
!        and benzene", Atmos. Chem. Phys., Vol 7, 3909-3922,
!        doi:10.5194/acp-7-3909-2007, 2007b.
!  (12) Ng, N.L., A.J. Kwan, J.D. Surratt, A.W.H. Chan, P.S. Chhabra, A.
!        Sorooshian, H.O.T. Pye, J.D. Crounse, P.O. Wennberg, R.C. Flagan, and
!        J.H. Seinfeld, "Secondary organic aerosol (SOA) formation from
!        reaction of isoprene with nitrate radicals (NO3)", Atmos. Chem. Phys.,
!        Vol 8, 4117-4140, doi:10.5194/acp-8-4117-2008, 2008.
!  (13) Pye, H.O.T., and J.H. Seinfeld, "A global perspective on aesorol from
!        low-volatility orgnaic compounds", Atmos. Chem. Phys., Vol 10, 4377-
!        4401, doi:10.5194/acp-10-4377-2010, 2010.
!  (14) Pye. H.O.T., A.W.H Chan, M.P. Barkley, and J.H. Seinfeld, "Global
!        modeling of organic aerosol: The importance of reactive nitrogen (NOx
!        and NO3)", Atmos. Chem. Phys., Vol 10, 11261-11276,
!        doi:10.5194/acp-10-11261-2010, 2010.
!  (15) Shilling, J.E., Q. Chen, S.M. King, T. Rosenoern, J.H. Kroll, D.R.
!        Worsnop, K.A. McKinney, S.T., Martin, "Particle mass yield in
!        secondary orgainc aerosol formed by the dark ozonolysis of a-pinene",
!        Atmos Chem Phys, Vol 8, 2073-2088, doi: 10.5194/acp-8-2073-2008, 2008.
!  (16) Shrivastava, M.K., E.M. Lipsky, C.O. Stanier, A.L. Robinson, "Modeling
!       semivolatile organic mass emissions from combustion systems", Environ.
!       Sci. Technol., Vol 40, 2671-2677, doi:10.1021/ES0522231, 2006.
!  (17) Zhang, J.Y., K.E.H. Hartz, S.N. Pandis, N.M. Donahue, "Secondary
!        organic aerosol formation from limonene ozonolysis: Homogeneous and
!        heterogeneous influences as a function of NOx", J. Phys. Chem. A, Vol
!        110, 11053-11063, doi:10.1021/Jp06286f, 2006.
!                                                                             .
!      Base Year is 2000. More at http://www.hiwater.org
!
! !REVISION HISTORY:
!  01 Apr 1994 - R. Park - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  ! Number of BB fires for parameterization (ramnarine 12/27/2018)
  REAL(f4), POINTER :: FIRE_NUM(:,:) => NULL()

  ! Molecules OH  per kg OH [molec/kg]
  REAL(fp), PARAMETER   :: XNUMOL_OH  = AVO / 17e-3_fp ! hard-coded MW
  REAL(fp), PARAMETER   :: CM3PERM3   = 1.e+6_fp

  ! SOAupdate:(hotp 5/20/10) new mtp
  ! one parent HC removed (only 3 instead of 4 monoterps)
  ! all monoterp and sesquiterp SOA lumped together
  ! NOX: now used to indicate high NOx (1),
  !      low NOx(2), and +NO3 (3) so MNOX is 3
  ! PROD: indicates # of volatilities/products
  INTEGER,  PARAMETER   :: MHC        = 11 ! max # HCs
  INTEGER,  PARAMETER   :: MSV        = 5  ! max # lumped semivols
  INTEGER,  PARAMETER   :: MPROD      = 4  ! max # volatility products
  INTEGER,  PARAMETER   :: MNOX       = 3  ! max # NOx levels/oxidants

  REAL(fp), PARAMETER   :: SMALLNUM   = 1e-20_fp

  ! Indicate number of parent HC based on simulation species
  ! (hotp 8/24/09)
  !INTEGER, SAVE        :: MAXSIMHC
  ! Now loop over number of semivolatiles (hotp 5/13/10)
  INTEGER,  SAVE        :: MAXSIMSV

  ! Identify parent hydrocarbon by numbers (hotp 5/12/10)
  INTEGER,  PARAMETER   :: PARENTMTPA = 1  ! bicyclic monoterpenes
  INTEGER,  PARAMETER   :: PARENTLIMO = 2  ! limonene
  INTEGER,  PARAMETER   :: PARENTMTPO = 3  ! other monoterpenes
  INTEGER,  PARAMETER   :: PARENTSESQ = 4  ! sesquiterpenes
  INTEGER,  PARAMETER   :: PARENTISOP = 5  ! isoprene
  INTEGER,  PARAMETER   :: PARENTBENZ = 6  ! aromatic benzene
  INTEGER,  PARAMETER   :: PARENTTOLU = 7  ! aromatic toluene
  INTEGER,  PARAMETER   :: PARENTXYLE = 8  ! aromatic xylene
  INTEGER,  PARAMETER   :: PARENTPOA  = 9  ! SVOCs (primary SVOCs)
  INTEGER,  PARAMETER   :: PARENTOPOA = 10 ! oxidized SVOCs (secondary SVOCs)
  INTEGER,  PARAMETER   :: PARENTNAP  = 11 ! IVOC surrogate (naphthalene)
  ! if NAP isn't last, check CHEM_NVOC

  ! NOx levels (oxidants) examined (hotp 5/13/10)
  INTEGER,  PARAMETER   :: NHIGHNOX   = 1  ! R + OH, RO2 + NO
  INTEGER,  PARAMETER   :: NLOWNOX    = 2  ! R + OH, RO2 + HO2
  INTEGER,  PARAMETER   :: NNO3RXN    = 3  ! R + NO3
  INTEGER,  PARAMETER   :: NONLYNOX   = 1  ! R + any oxidant
!
! !LOCAL VARIABLES:
!
  ! Scalars
  ! Rate constant for RO2+NO and RO2+HO2
  ! k=Aexp(B/T) like globchem.dat (hotp 5/7/10)
  REAL(fp)              :: AARO2NO,  BBRO2NO
  REAL(fp)              :: AARO2HO2, BBRO2HO2

  ! Arrays
  INTEGER               :: NPROD(MSV) !hotp 5/13/10 now MSV not MHC
  INTEGER               :: NNOX(MSV)  !hotp 5/13/10
  ! now only 4 offline oxidations (hotp 5/20/10)
  REAL(fp)              :: KO3_REF(4), KOH_REF(4), KNO3_REF(4)
  ! KOM_REF now has dims of MPROD, MSV (hotp 5/22/10)
  REAL(fp)              :: KOM_REF(MPROD,MSV)
  REAL(fp)              :: ALPHA(MNOX,MPROD,MHC)

  ! Array for mapping parent HC to semivolatiles (SV) (hotp 5/14/10)
  INTEGER               :: IDSV(MHC)

  ! Diagnostic that tracks how much parent HC reacts
  ! with each allowed reactant (hotp 5/24/10)
  REAL(fp)              :: DELTAHCSAVE(MNOX,MHC)

  REAL(fp), ALLOCATABLE :: BCCONV(:,:,:)
  REAL(fp), ALLOCATABLE :: OCCONV(:,:,:)
  REAL(fp), ALLOCATABLE :: TCOSZ(:,:)
  REAL(fp), ALLOCATABLE :: ORVC_SESQ(:,:,:)
  REAL(fp), ALLOCATABLE :: GLOB_DARO2(:,:,:,:,:) ! Diagnostic (dkh, 11/10/06)

#ifdef TOMAS
  REAL(fp), ALLOCATABLE, TARGET :: BCFF(:,:,:,:)
  REAL(fp), ALLOCATABLE, TARGET :: OCFF(:,:,:,:)
  REAL(fp), ALLOCATABLE, TARGET :: BCBF(:,:,:,:)
  REAL(fp), ALLOCATABLE, TARGET :: OCBF(:,:,:,:)
  REAL(fp), ALLOCATABLE, TARGET :: BCBB(:,:,:,:)
  REAL(fp), ALLOCATABLE, TARGET :: OCBB(:,:,:,:)

  REAL(fp), ALLOCATABLE, TARGET :: BCPI_ANTH_BULK(:,:)
  REAL(fp), ALLOCATABLE, TARGET :: BCPO_ANTH_BULK(:,:)
  REAL(fp), ALLOCATABLE, TARGET :: OCPI_ANTH_BULK(:,:)
  REAL(fp), ALLOCATABLE, TARGET :: OCPO_ANTH_BULK(:,:)

  REAL(fp), ALLOCATABLE, TARGET :: BCPI_BIOB_BULK(:,:)
  REAL(fp), ALLOCATABLE, TARGET :: BCPO_BIOB_BULK(:,:)
  REAL(fp), ALLOCATABLE, TARGET :: OCPI_BIOB_BULK(:,:)
  REAL(fp), ALLOCATABLE, TARGET :: OCPO_BIOB_BULK(:,:)

  REAL(fp), ALLOCATABLE :: TERP_ORGC(:,:)
  REAL(fp), ALLOCATABLE :: CO_ANTH(:,:)
#endif

  ! Diagnostic that tracks how much SOA is formed/evaporated
  ! (hotp 6/5/10)
  REAL(fp), ALLOCATABLE :: SPECSOAPROD(:,:,:,:,:)
  REAL(fp), ALLOCATABLE :: SPECSOAEVAP(:,:,:,:,:)

  ! semivolpoa4: diagnostic to keep track of POG reacted (hotp 3/27/09)
  REAL(fp), SAVE, ALLOCATABLE :: GLOB_POGRXN(:,:,:,:)

  ! diagnostic added for RO2 branching ratio (hotp 5/24/10)
  REAL(fp), SAVE, ALLOCATABLE :: BETANOSAVE(:,:,:)

  ! semivolpoa2: array for POA emissions (hotp 2/27/09)
  REAL(fp), SAVE, ALLOCATABLE :: POAEMISS(:,:,:,:)

  ! Array for initial OA+OG (hotp 5/17/10)
  ! Diagnostic only, dims: I,J,L,MPROD,MSV (hotp 5/22/10)
  REAL(fp), SAVE, ALLOCATABLE :: OAGINITSAVE(:,:,:,:,:)

  ! Array for change in SOG (diagnostic) (hotp 5/17/10)
  ! dims: I,J,L,MNOX,MHC
  REAL(fp), SAVE, ALLOCATABLE :: DELTASOGSAVE(:,:,:,:,:)

  ! Days per month (based on 1998)
  INTEGER :: NDAYS(12) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

  ! Pointers to fields contained in the HEMCO data structure
  ! These must all be declared as REAL(f4), aka REAL*4.
  REAL(f4), POINTER     :: O3(:,:,:)  => NULL()
  REAL(f4), POINTER     :: OH(:,:,:)  => NULL()
  REAL(f4), POINTER     :: NO3(:,:,:) => NULL()

  ! Species ID flags
  INTEGER :: id_ASOG1,   id_ASOG2,  id_ASOG3,  id_ASOA1,  id_ASOA2
  INTEGER :: id_ASOA3,   id_ASOAN,  id_AW1,    id_BCPI,   id_BCPO
  INTEGER :: id_BENZ,    id_ECIL1,  id_ECOB1,  id_HO2
  INTEGER :: id_ISOP,    id_LIMO,   id_MTPA
  INTEGER :: id_MTPO,    id_NAP,    id_NK1,    id_NH4,    id_NO
  INTEGER :: id_NO3,     id_OCIL1,  id_OCOB1,  id_O3,     id_OH
  INTEGER :: id_OCPO,    id_OCPI,   id_OPOA1,  id_OPOG1,  id_OPOA2
  INTEGER :: id_OPOG2,   id_POA1,   id_POA2,   id_POG1,   id_POG2
  INTEGER :: id_TOLU,    id_TSOA0,  id_TSOA1
  INTEGER :: id_TSOA2,   id_TSOA3,  id_TSOG0,  id_TSOG1,  id_TSOG2
  INTEGER :: id_TSOG3,   id_XYLE,   id_LBRO2N, id_LBRO2H, id_LTRO2N
  INTEGER :: id_LTRO2H,  id_LXRO2N, id_LXRO2H, id_LNRO2N, id_LNRO2H
  INTEGER :: id_LISOPOH, id_LISOPNO3
  INTEGER :: id_SOAS,    id_SOAP

#ifdef APM
  REAL(fp), ALLOCATABLE :: BCCONVNEW(:,:,:)
  REAL(fp), ALLOCATABLE :: OCCONVNEW(:,:,:)
#endif

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chemcarbon
!
! !DESCRIPTION: Subroutine CHEMCARBON is the interface between the GEOS-Chem
!  main program and the carbon aerosol chemistry routines that calculates dry
!  deposition, chemical conversion between hydrophilic and hydrophobic, and
!  SOA production.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEMCARBON( Input_Opt,  State_Chm, State_Diag, &
                         State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : DEBUG_MSG
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
    USE TIME_MOD,           ONLY : GET_TS_CHEM
#ifdef APM
    USE APM_INIT_MOD,       ONLY : APMIDS
    USE APM_INIT_MOD,       ONLY : NBCOC,CEMITBCOC1
    USE HCO_DIAGN_MOD
    USE HCO_ERROR_MOD
    USE HCO_TYPES_MOD,      ONLY : DiagnCont
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_STATE_MOD,      ONLY : HCO_GetHcoID
#endif
#ifdef BPCH_DIAG
    USE CMN_O3_MOD,         ONLY : SAVEOA
#endif
#ifdef TOMAS
    USE TOMAS_MOD,          ONLY : SOACOND, IBINS              !(win, 1/25/10)
    USE TOMAS_MOD,          ONLY : CHECKMN                     !(sfarina)
    USE PRESSURE_MOD,       ONLY : GET_PCENTER
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  01 Apr 1994 - R. Park - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVED scalars
    LOGICAL, SAVE      :: FIRSTCHEM = .TRUE.

    ! Scalars
    LOGICAL            :: prtDebug
    LOGICAL            :: IT_IS_AN_AEROSOL_SIM
    LOGICAL            :: LSOA
    LOGICAL            :: LEMIS
    REAL(fp)           :: NEWSOA
    REAL(fp)           :: DTCHEM, SOAP_LIFETIME  ! [=] seconds
    INTEGER            :: L

#ifdef TOMAS
    INTEGER            :: I, J
    REAL*4             :: BOXVOL, TEMPTMS, PRES
#endif

    ! Pointers
    REAL(fp), POINTER  :: Spc(:,:,:,:)

    ! For getting fields from HEMCO
    CHARACTER(LEN=255) :: LOC = 'CHEMCARBON (carbon_mod.F90)'

#ifdef APM
    TYPE(DiagnCont), POINTER :: DiagnCnt
    INTEGER            :: FLAG,I,J,N,IDCARBON
    REAL(fp)           :: A_M2, E_CARBON, DTSRCE
    REAL(fp)           :: EMITRATE(State_Grid%NX,State_Grid%NY)
#endif

    !=================================================================
    ! CHEMCARBON begins here!
    !=================================================================

    ! Assume success
    RC                   = GC_SUCCESS

    ! Copy fields from INPUT_OPT to local variables for use below
    LSOA                 = Input_Opt%LSOA
    LEMIS                = Input_Opt%LEMIS
    IT_IS_AN_AEROSOL_SIM = Input_Opt%ITS_AN_AEROSOL_SIM

    DTCHEM               = GET_TS_CHEM()
    !                      ( days )*(hrs/day)*(mins/hr)*(sec/min)
    SOAP_LIFETIME        = 1.00_fp * 24.0_fp * 60.0_fp * 60.0_fp

    ! Do we have to print debug output?
    prtDebug             = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Point to chemical species array [kg]
    Spc                  => State_Chm%Species

    ! First-time initialization
    IF ( FIRSTCHEM ) THEN

       ! Zero SOG4 and SOA4 (SOA from ISOP in gas & aerosol form)
       ! for offline aerosol simulations.  Eventually we should have
       ! archived isoprene oxidation fields available for offline
       ! simulations but for now we just set them to zero.
       ! (dkh, bmy, 6/1/06)
       IF ( IT_IS_AN_AEROSOL_SIM ) THEN

          ! lumped arom/IVOC (hotp 5/17/10)
          ! LUMPAROMIVOC: lump arom/IVOC not supported for offline sims
          IF ( id_ASOAN > 0 ) Spc(:,:,:,id_ASOAN) = 0.0_fp
          IF ( id_ASOA1 > 0 ) Spc(:,:,:,id_ASOA1) = 0.0_fp
          IF ( id_ASOA2 > 0 ) Spc(:,:,:,id_ASOA2) = 0.0_fp
          IF ( id_ASOA3 > 0 ) Spc(:,:,:,id_ASOA3) = 0.0_fp
          IF ( id_ASOG1 > 0 ) Spc(:,:,:,id_ASOG1) = 0.0_fp
          IF ( id_ASOG2 > 0 ) Spc(:,:,:,id_ASOG2) = 0.0_fp
          IF ( id_ASOG3 > 0 ) Spc(:,:,:,id_ASOG3) = 0.0_fp

          IF ( LSOA ) THEN
             ! Get offline oxidant fields from HEMCO (mps, 9/23/14)
             CALL HCO_GetPtr( HcoState, 'GLOBAL_OH',  OH,  RC )
             IF ( RC /= GC_SUCCESS ) &
                CALL ERROR_STOP( 'Cannot get pointer to GLOBAL_OH',  LOC)

             CALL HCO_GetPtr( HcoState, 'GLOBAL_NO3', NO3, RC )
             IF ( RC /= GC_SUCCESS ) &
                CALL ERROR_STOP( 'Cannot get pointer to GLOBAL_NO3', LOC)

             CALL HCO_GetPtr( HcoState, 'O3',         O3,  RC )
             IF ( RC /= GC_SUCCESS ) &
                CALL ERROR_STOP( 'Cannot get pointer to O3',         LOC)
          ENDIF

          ! initialize SOA Precursor and SOA for simplified SOA sims
          IF ( id_SOAP > 0 ) THEN
             Spc(:,:,:,id_SOAP) = 0e+0_fp
             Spc(:,:,:,id_SOAS) = 0e+0_fp
          ENDIF

       ENDIF

       ! Determine number of semivolatile parent HC (hotp 8/24/09)
       !MAXSIMHC = 0 ! for non-volatile sim
       ! Now use SV instead of HC (hotp 5/13/10)
       MAXSIMSV = 0
       IF ( LSOA ) THEN
          ! updated (hotp 5/20/10) new mtp
          MAXSIMSV = 3  ! mono+sesq (1) + isop (2) + aromatics (3)
          IF ( id_POA1  > 0 ) MAXSIMSV = MAXSIMSV + 1
          IF ( id_OPOA1 > 0 ) MAXSIMSV = MAXSIMSV + 1
          IF ( MAXSIMSV > MSV ) THEN
             CALL ERROR_STOP('YOUVE GOT A PROBLEM W/ SEMIVOLATILES', &
                             'carbon_mod.F90')
          ENDIF

          ! Print to log for record
          IF ( prtDebug ) THEN
             print*,'Number of SOA semivols (MAXSIMSV): ', MAXSIMSV
             print*,'This number should be 5 for semivol POA' ! hotp 5/20/10
          ENDIF
       ENDIF

       ! Reset first-time flag
       FIRSTCHEM = .FALSE.
    ENDIF

    !=================================================================
    ! Do chemistry for carbon aerosol species
    !=================================================================

    ! Chemistry for hydrophobic BC
    IF ( id_BCPO > 0 ) THEN
       CALL CHEM_BCPO( Input_Opt, State_Diag, State_Grid, &
                       Spc(:,:,:,id_BCPO),   RC          )
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### CHEMCARBON: a CHEM_BCPO' )
       ENDIF
    ENDIF

    ! Chemistry for hydrophilic BC
    IF ( id_BCPI > 0 ) THEN
       CALL CHEM_BCPI( Input_Opt, State_Diag, State_Grid, &
                       Spc(:,:,:,id_BCPI),   RC )
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### CHEMCARBON: a CHEM_BCPI' )
       ENDIF
    ENDIF

    ! Chemistry for hydrophobic OC (traditional POA only)
    IF ( id_OCPO > 0 ) THEN
       CALL CHEM_OCPO( Input_Opt, State_Diag, State_Grid, &
                       Spc(:,:,:,id_OCPO),   RC          )
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### CHEMCARBON: a CHEM_OCPO' )
       ENDIF
    ENDIF

    ! Chemistry for hydrophilic OC (traditional POA only)
    IF ( id_OCPI > 0 ) THEN
       CALL CHEM_OCPI( Input_Opt, State_Diag, State_Grid, &
                       Spc(:,:,:,id_OCPI),   RC )
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### CHEMCARBON: a CHEM_OCPI' )
       ENDIF
    ENDIF

#ifdef APM
    !=====================================================================
    ! APM Microphysics
    !=====================================================================
    CALL BCDRY_SETTLINGBIN( Input_Opt,  State_Chm, State_Diag, &
                            State_Grid, State_Met, RC )

    CALL OCDRY_SETTLINGBIN( Input_Opt,  State_Chm, State_Diag, &
                            State_Grid, State_Met, RC )

    IF ( ( id_BCPO+id_BCPI ) > 2 )THEN

       ! Get biomass BCPO diagnostics from HEMCO
       DiagnCnt => NULL()
       CALL Diagn_Get( HcoState, .FALSE., DiagnCnt, &
                       FLAG,  RC, cName='BIOMASS_BCPO',        &
                       AutoFill=1,                             &
                       COL=HcoState%Diagn%HcoDiagnIDManual )

       ! Add into EMITRATE array (or set to zero if not found)
       IF ( FLAG == HCO_SUCCESS ) THEN
          EMITRATE = DiagnCnt%Arr2D%Val
       ELSE
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, ' (a)' )'APM not found BIOMASS_BCPO'
          ENDIF
          EMITRATE = 0.0e+0_fp
       ENDIF

       ! Get anthropogenic BCPO diagnostics from HEMCO
       DiagnCnt => NULL()
       CALL Diagn_Get( HcoState, .FALSE., DiagnCnt, &
                       FLAG,  RC, cName='ANTHROPOGENIC_BCPO',  &
                       AutoFill=1,                             &
                       COL=HcoState%Diagn%HcoDiagnIDManual )


       IF ( FLAG==HCO_SUCCESS ) THEN
          !$OMP PARALLEL DO       &
          !$OMP DEFAULT( SHARED ) &
          !$OMP PRIVATE( J, I )
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX
             EMITRATE(I,J) = EMITRATE(I,J) + &
                             SUM(DiagnCnt%Arr3D%Val(I,J,1:State_Grid%NZ))

             IF( EMITRATE(I,J) > 0.0e+0_fp ) THEN
                EMITRATE(I,J) = 1.0e+0_fp - &
                                SUM(DiagnCnt%Arr3D%Val(I,J,:))/EMITRATE(I,J)
             ELSE
                EMITRATE(I,J) = 0.0e+0_fp
             ENDIF
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO
       ELSE
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, '(a)' ) 'APM not found ANTHROPOGENIC_BCPO'
          ENDIF
          EMITRATE = 0.0e+0_fp
       ENDIF

       ! Emission timestep
       DTSRCE = HcoState%TS_EMIS

       IDCARBON   = HCO_GetHcoID( 'BCPO',   HcoState )
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( L, J, I, A_M2, E_CARBON, N )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Grid box surface area [m2]
          A_M2 = HcoState%Grid%AREA_M2%Val( I, J )

          ! Get emissions [kg/m2/s] and convert to [kg/box]
          E_CARBON = HcoState%Spc(IDCARBON)%Emis%Val(I,J,L) * A_M2 * DTSRCE

          DO N = 1, NBCOC
             Spc(I,J,L,APMIDS%id_BCBIN1+N-1)= &
                Spc(I,J,L,APMIDS%id_BCBIN1+N-1)+ &
                E_CARBON*( &
                (1.0e+0_fp-EMITRATE(I,J))*CEMITBCOC1(N,1)+ &
                EMITRATE(I,J)*CEMITBCOC1(N,2))
          ENDDO
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

       IDCARBON   = HCO_GetHcoID( 'BCPI',   HcoState )
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( L, J, I, A_M2, E_CARBON, N )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Grid box surface area [m2]
          A_M2 = HcoState%Grid%AREA_M2%Val( I, J )

          ! Get emissions [kg/m2/s] and convert to [kg/box]
          E_CARBON = HcoState%Spc(IDCARBON)%Emis%Val(I,J,L) * A_M2 * DTSRCE

          DO N = 1, NBCOC
             Spc(I,J,L,APMIDS%id_BCBIN1+N-1)= &
                Spc(I,J,L,APMIDS%id_BCBIN1+N-1)+ &
                E_CARBON*( &
                (1.0e+0_fp-EMITRATE(I,J))*CEMITBCOC1(N,1)+ &
                EMITRATE(I,J)*CEMITBCOC1(N,2))
          ENDDO
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    !GanLuo add OCs
    IF ( ( id_OCPO + id_OCPI ) > 2 ) THEN
       ! Get diagnostics from HEMCO
       DiagnCnt => NULL()
       CALL Diagn_Get( HcoState, .FALSE., DiagnCnt, &
                       FLAG,  RC, cName='BIOMASS_OCPO',        &
                       AutoFill=1,                             &
                       COL=HcoState%Diagn%HcoDiagnIDManual )

       ! Add biomass OCPO diagnostic to EMITRATE (if it's found)
       IF ( FLAG==HCO_SUCCESS ) THEN
          EMITRATE = DiagnCnt%Arr2D%Val
       ELSE
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, '(a)' )'APM not found BIOMASS_OCPO'
          ENDIF
          EMITRATE = 0.0e+0_fp
       ENDIF

       DiagnCnt => NULL()
       CALL Diagn_Get( HcoState, .FALSE., DiagnCnt, &
                       FLAG,  RC, cName='ANTHROPOGENIC_OCPO',  &
                       AutoFill=1,                             &
                       COL=HcoState%Diagn%HcoDiagnIDManual )

       ! Add anthropogenic OCPO diagnostic to EMITRATE (if it's found)
       IF(FLAG==HCO_SUCCESS)THEN
          !$OMP PARALLEL DO       &
          !$OMP DEFAULT( SHARED ) &
          !$OMP PRIVATE( J, I )
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX
             EMITRATE(I,J) = EMITRATE(I,J) + &
                             SUM(DiagnCnt%Arr3D%Val(I,J,1:State_Grid%NZ))
             IF(EMITRATE(I,J) > 0.0e+0_fp )THEN
               EMITRATE(I,J) = 1.0e+0_fp- &
                               SUM(DiagnCnt%Arr3D%Val(I,J,:))/EMITRATE(I,J)
            ELSE
               EMITRATE(I,J) = 0.0e+0_fp
            ENDIF
         ENDDO
         ENDDO
         !$OMP END PARALLEL DO
      ELSE
         IF ( Input_Opt%amIRoot ) THEN
            WRITE( 6, '(a)' )'APM not found ANTHROPOGENIC_OCPO'
         ENDIF
         EMITRATE = 0.0e+0_fp
      ENDIF

      ! Emission timestep
      DTSRCE = HcoState%TS_EMIS

      IDCARBON   = HCO_GetHcoID( 'OCPO',   HcoState )
      !$OMP PARALLEL DO       &
      !$OMP DEFAULT( SHARED ) &
      !$OMP PRIVATE( L, J, I, A_M2, E_CARBON, N )
      DO L = 1, State_Grid%NZ
      DO J = 1, State_Grid%NY
      DO I = 1, State_Grid%NX

         ! Grid box surface area [m2]
         A_M2 = HcoState%Grid%AREA_M2%Val( I, J )

         ! Get emissions [kg/m2/s] and convert to [kg/box]
         E_CARBON = HcoState%Spc(IDCARBON)%Emis%Val(I,J,L) * A_M2 * DTSRCE

         DO N=1,NBCOC
            Spc(I,J,L,APMIDS%id_OCBIN1+N-1)= &
               Spc(I,J,L,APMIDS%id_OCBIN1+N-1)+ &
               E_CARBON*( &
               (1.0e+0_fp-EMITRATE(I,J))*CEMITBCOC1(N,1)+ &
               EMITRATE(I,J)*CEMITBCOC1(N,2))
         ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      IDCARBON   = HCO_GetHcoID( 'OCPI',   HcoState )
      !$OMP PARALLEL DO       &
      !$OMP DEFAULT( SHARED ) &
      !$OMP PRIVATE( L, J, I, A_M2, E_CARBON, N )
      DO L = 1, State_Grid%NZ
      DO J = 1, State_Grid%NY
      DO I = 1, State_Grid%NX

         ! Grid box surface area [m2]
         A_M2 = HcoState%Grid%AREA_M2%Val( I, J )

         ! Get emissions [kg/m2/s] and convert to [kg/box]
         E_CARBON = HcoState%Spc(IDCARBON)%Emis%Val(I,J,L) * A_M2 * DTSRCE

         DO N=1,NBCOC
            Spc(I,J,L,APMIDS%id_OCBIN1+N-1)= &
               Spc(I,J,L,APMIDS%id_OCBIN1+N-1)+ &
               E_CARBON*( &
               (1.D0-EMITRATE(I,J))*CEMITBCOC1(N,1)+ &
               EMITRATE(I,J)*CEMITBCOC1(N,2))
         ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   ENDIF

   !write(*,'(a4,10e15.7)')'Luo1', &
   ! (sum(Spc(:,:,1,id_BCPO))+sum(Spc(:,:,1,id_BCPI))), &
   !  sum(Spc(:,:,1,APMIDS%id_BCBIN1:(APMIDS%id_BCBIN1+14)))

   IF( ( id_POG1 + id_POA1 ) > 2 )THEN
      !$OMP PARALLEL DO       &
      !$OMP DEFAULT( SHARED ) &
      !$OMP PRIVATE( L, J, I )
      DO L = 1, State_Grid%NZ
      DO J = 1, State_Grid%NY
      DO I = 1, State_Grid%NX
         DO N=1,NBCOC
            Spc(I,J,L,APMIDS%id_OCBIN1+N-1)= &
               Spc(I,J,L,APMIDS%id_OCBIN1+N-1)+ &
               POAEMISS(I,J,L,1)*0.9d0*CEMITBCOC1(N,1)
         ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO
   ENDIF
#endif

#ifdef TOMAS
   CALL CHECKMN( 0, 0, 0, Input_Opt, State_Chm, State_Grid, &
                 State_Met, 'CHECKMN from chemcarbon', RC)
   ! Chemistry (aging) for size-resolved EC and OC (win, 1/25/10)
   IF ( id_ECIL1 > 0 .and. id_ECOB1 > 0 ) THEN
      CALL AGING_CARB( State_Grid, &
                       Spc(:,:,:,id_ECIL1:id_ECIL1+IBINS-1), &
                       Spc(:,:,:,id_ECOB1:id_ECOB1+IBINS-1) )
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMCARBO: AGING_CARB EC' )
      ENDIF
   ENDIF
   IF ( id_OCIL1 > 0 .and. id_OCOB1 > 0 ) THEN
      CALL AGING_CARB( State_Grid, &
                       Spc(:,:,:,id_OCIL1:id_OCIL1+IBINS-1), &
                       Spc(:,:,:,id_OCOB1:id_OCOB1+IBINS-1) )
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMCARBO: AGING_CARB OC' )
      ENDIF
   ENDIF
#endif

   IF ( id_SOAP > 0 ) THEN
      ! AGE SOAP -> SOA

#ifdef TOMAS
      CALL CHECKMN( 0, 0, 0, Input_Opt, State_Chm, State_Grid, &
                    State_Met, 'CHECKMN from chemcarbon', RC)

      !$OMP PARALLEL DO       &
      !$OMP DEFAULT( SHARED ) &
      !$OMP PRIVATE( I, J, L, NEWSOA, BOXVOL, TEMPTMS, PRES )
      DO L = 1, State_Grid%NZ
      DO J = 1, State_Grid%NY
      DO I = 1, State_Grid%NX
         NEWSOA  = Spc(I,J,L,id_SOAP) * (1.e+0_fp - DEXP(-DTCHEM/SOAP_LIFETIME))
         BOXVOL  = State_Met%AIRVOL(I,J,L) * 1.e6 !convert from m3 -> cm3
         TEMPTMS = State_Met%T(I,J,L)
         PRES    = GET_PCENTER(I,j,L)*100.0 ! in Pa
         IF ( NEWSOA > 0.0e+0_fp ) THEN
            !sfarina16: SOAP -> size Resolved TOMAS SOA
            CALL SOACOND( NEWSOA, I, J, L, BOXVOL, TEMPTMS, PRES, &
                          State_Chm, State_Grid, RC)
         ENDIF
         Spc(I,J,L,id_SOAS) = Spc(I,J,L,id_SOAS) + NEWSOA
         Spc(I,J,L,id_SOAP) = Spc(I,J,L,id_SOAP) - NEWSOA
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO
#else
      !$OMP PARALLEL DO       &
      !$OMP DEFAULT( SHARED ) &
      !$OMP PRIVATE( L, NEWSOA )
      DO L = 1, State_Grid%NZ
         !NEWSOA used in a different context than above.
         !above is absolute mass, here is a relative decay factor
         NEWSOA = DEXP(-DTCHEM/SOAP_LIFETIME)
         Spc(:,:,L,id_SOAS) = Spc(:,:,L,id_SOAS) + &
                              Spc(:,:,L,id_SOAP) * (1.0_fp - NEWSOA)
         Spc(:,:,L,id_SOAP) = Spc(:,:,L,id_SOAP) * NEWSOA
      ENDDO
      !$OMP END PARALLEL DO
#endif

      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMCARBO: SIMPLIFIED SOA')
      ENDIF
   ENDIF

   ! Free pointer
   Spc => NULL()

   !=================================================================
   ! Do chemistry for secondary organic aerosols
   !
   ! %%% NOTE: We are not planning on using the SOA mechanism   %%%
   ! %%% with the ESMF interface at this point. (bmy, 11/14/12) %%%
   !=================================================================
   IF ( LSOA ) THEN

      IF ( IT_IS_AN_AEROSOL_SIM ) THEN

         ! Compute time scaling arrays for offline OH, NO3
         ! but only if it hasn't been done in EMISSCARBON
         IF ( LSOA .and. ( .not. LEMIS ) ) THEN
            CALL OHNO3TIME( State_Grid )
            IF ( prtDebug ) THEN
               CALL DEBUG_MSG( '### CHEMCARB: a OHNO3TIME' )
            ENDIF
         ENDIF

      ENDIF

      ! Compute SOA chemistry
      ! NOTE: This is SOA production from the reversible mechanism only
      ! (tmf, 12/07/07)
      CALL SOA_CHEMISTRY( Input_Opt,  State_Chm, State_Diag, &
                          State_Grid, State_Met, RC )

      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMCARBON: a SOA_CHEM' )
      ENDIF

#ifdef BPCH_DIAG
      ! NOTE: This is only needed for ND51, ND51b (bmy, 10/4/19)
      ! Get total organic aerosol:
      CALL OASAVE( SAVEOA, Input_Opt, State_Chm, State_Grid, State_Met, RC )
#endif
   ENDIF

 END SUBROUTINE CHEMCARBON
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_bcpo
!
! !DESCRIPTION: Subroutine CHEM\_BCPO converts hydrophobic BC to hydrophilic
!  BC and calculates the dry deposition of hydrophobic BC.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE CHEM_BCPO( Input_Opt, State_Diag, State_Grid, TC, RC )
!
! !USES:
!
   USE ErrCode_Mod
   USE Input_Opt_Mod,  ONLY : OptInput
   USE State_Diag_Mod, ONLY : DgnState
   USE State_Grid_Mod, ONLY : GrdState
   USE TIME_MOD,       ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt             ! Input Options
   TYPE(GrdState), INTENT(IN)    :: State_Grid            ! Grid State
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(DgnState), INTENT(INOUT) :: State_Diag            ! Diags State
   REAL(fp),       INTENT(INOUT) :: TC(State_Grid%NX,State_Grid%NY,State_Grid%NZ) ! H-phobic BC [kg]
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC                    ! Success?
!
! !REMARKS:
!  Drydep is now applied in mixing_mod.F90.
!
! !REVISION HISTORY:
!  01 Apr 2004 - R. Park - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Scalars
   INTEGER             :: I,      J,   L
   REAL(fp)            :: DTCHEM, KBC, FREQ, TC0, CNEW, RKT
!
! !DEFINED PARAMETERS:
!
   REAL(fp), PARAMETER :: BC_LIFE = 1.15e+0_fp

   !=================================================================
   ! CHEM_BCPO begins here!
   !=================================================================

   ! Assume success
   RC        = GC_SUCCESS

   ! Initialize
   KBC       = 1.e+0_fp / ( 86400e+0_fp * BC_LIFE )
   DTCHEM    = GET_TS_CHEM()
   BCCONV    = 0e+0_fp

   !=================================================================
   ! For species with dry deposition, the loss rate of dry dep is
   ! combined in chem loss term.
   !
   ! Conversion from hydrophobic to hydrophilic:
   ! e-folding time 1.15 days
   ! ----------------------------------------
   ! Use an e-folding time of 1.15 days or a convertion rate
   ! of 1.0e-5 /sec.
   !
   ! Hydrophobic(2) --> Hydrophilic(1) ,  k  = 1.0e-5
   ! Both aerosols are dry-deposited,     kd = Dvel/DELZ (sec-1)
   !=================================================================
   !$OMP PARALLEL DO       &
   !$OMP DEFAULT( SHARED ) &
   !$OMP PRIVATE( I, J, L, TC0, FREQ, RKT, CNEW ) &
   !$OMP SCHEDULE( DYNAMIC )
   DO L = 1, State_Grid%NZ
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      ! Initial BC mass [kg]
      TC0  = TC(I,J,L)

      ! Zero drydep freq
      ! ### NOTE: Remove this later, but need to make
      ! ### sure we don't incur numerical diffs (bmy, 6/12/15)
      FREQ = 0e+0_fp

      ! Amount of BCPO left after chemistry and drydep [kg]
      RKT  = ( KBC + FREQ ) * DTCHEM
      CNEW = TC0 * EXP( -RKT )

      ! Prevent underflow condition
      IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

      ! Amount of BCPO converted to BCPI [kg/timestep]
      BCCONV(I,J,L) = ( TC0 - CNEW ) * KBC / ( KBC + FREQ )

      !==============================================================
      ! HISTORY (aka netCDF diagnostics)
      !
      ! Archive production of hydrophilic black carbon (BCPI) from
      ! hydrophobic black carbon (BCPO) [kg]
      !
      ! NOTE: Consider converting to area-independent units kg/m2/s.
      !==============================================================
      IF ( State_Diag%Archive_ProdBCPIfromBCPO ) THEN
         State_Diag%ProdBCPIfromBCPO(I,J,L) = BcConv(I,J,L)
      ENDIF

      ! Store new concentration back into species array
      TC(I,J,L) = CNEW
   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

 END SUBROUTINE CHEM_BCPO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_bcpi
!
! !DESCRIPTION: Subroutine CHEM\_BCPI calculates dry deposition of
!  hydrophilic BC.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE CHEM_BCPI( Input_Opt, State_Diag, State_Grid, TC, RC )
!
! !USES:
!
   USE ErrCode_Mod
   USE Input_Opt_Mod,  ONLY : OptInput
   USE State_Diag_Mod, ONLY : DgnState
   USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
   TYPE(OptInput), INTENT(IN)    :: Input_Opt             ! Input Options
   TYPE(GrdState), INTENT(IN)    :: State_Grid            ! Grid State
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(DgnState), INTENT(INOUT) :: State_Diag            ! Diags State
   REAL(fp),       INTENT(INOUT) :: TC(State_Grid%NX,State_Grid%NY,State_Grid%NZ) ! H-philic BC [kg]
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC                    ! Success?
!
! !REMARKS:
!  Drydep is now applied in mixing_mod.F90.
!
! !REVISION HISTORY:
!  01 Apr 2004 - R. Park - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Scalars
   INTEGER  :: I,      J,      L
   REAL(fp) :: L_FRAC, TC0, CNEW, CCV

   !=================================================================
   ! CHEM_BCPI begins here!
   !=================================================================

   ! Assume success
   RC = GC_SUCCESS

   !$OMP PARALLEL DO       &
   !$OMP DEFAULT( SHARED ) &
   !$OMP PRIVATE( I, J, L, TC0, CCV, CNEW ) &
   !$OMP SCHEDULE( DYNAMIC )
   DO L = 1, State_Grid%NZ
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      ! Initial H-philic BC [kg]
      TC0 = TC(I,J,L)

      ! H-philic BC that used to be H-phobic BC [kg]
      CCV = BCCONV(I,J,L)

      ! Add the amount of converted BCPO to BCPI
      CNEW = TC0 + CCV

      ! Prevent underflow condition
      IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

      ! Save new concentration of H-philic IC in species array
      TC(I,J,L) = CNEW

   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   !=================================================================
   ! Zero out the BCCONV array for the next iteration
   !=================================================================
   BCCONV = 0e+0_fp

 END SUBROUTINE CHEM_BCPI
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_ocpo
!
! !DESCRIPTION: Subroutine CHEM\_OCPO converts hydrophobic OC to hydrophilic
!  OC and calculates the dry deposition of hydrophobic OC.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE CHEM_OCPO( Input_Opt, State_Diag, State_Grid, TC, RC )
!
! !USES:
!
   USE ErrCode_Mod
   USE Input_Opt_Mod,  ONLY : OptInput
   USE State_Diag_Mod, ONLY : DgnState
   USE State_Grid_Mod, ONLY : GrdState
   USE TIME_MOD,       ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt             ! Input Options
   TYPE(GrdState), INTENT(IN)    :: State_Grid            ! Grid State
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(DgnState), INTENT(INOUT) :: State_Diag            ! Diags State
   REAL(fp),       INTENT(INOUT) :: TC(State_Grid%NX,State_Grid%NY,State_Grid%NZ) ! H-phobic OC [kg]
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC                    ! Success?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  01 Apr 2004 - R. Park - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Scalars
   INTEGER              :: I,      J,   L
   REAL(fp)             :: DTCHEM, KOC, TC0, CNEW, RKT, FREQ
!
! !DEFINED PARAMETERS:
!
   REAL(fp),  PARAMETER :: OC_LIFE = 1.15e+0_fp

   !=================================================================
   ! CHEM_OCPO begins here!
   !=================================================================

   ! Assume success
   RC        = GC_SUCCESS

   ! Initialize
   KOC       = 1.e+0_fp / ( 86400e+0_fp * OC_LIFE )
   DTCHEM    = GET_TS_CHEM()
   OCCONV    = 0e+0_fp

   !=================================================================
   ! For species with dry deposition, the loss rate of dry dep is
   ! combined in chem loss term.
   !
   ! Conversion from hydrophobic to hydrophilic:
   ! e-folding time 1.15 days
   ! ----------------------------------------
   ! Use an e-folding time of 1.15 days or a convertion rate
   ! of 1.0e-5 /sec.
   !    Hydrophobic --> Hydrophilic,  k  = 1.0e-5
   !    Aerosols are dry-deposited,   kd = DEPSAV (sec-1)
   !=================================================================
   !$OMP PARALLEL DO       &
   !$OMP DEFAULT( SHARED ) &
   !$OMP PRIVATE( I, J, L, TC0, FREQ, RKT, CNEW ) &
   !$OMP SCHEDULE( DYNAMIC )
   DO L = 1, State_Grid%NZ
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      ! Initial OC [kg]
      TC0  = TC(I,J,L)

      ! Zero drydep freq
      ! ### NOTE: Remove this later, but need to make
      ! ### sure we don't incur numerical diffs (bmy, 6/12/15)
      FREQ = 0e+0_fp

      ! Amount of OCPO left after chemistry and drydep [kg]
      RKT  = ( KOC + FREQ ) * DTCHEM
      CNEW = TC0 * EXP( -RKT )

      ! Prevent underflow condition
      IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

      ! Amount of OCPO converted to OCPI [kg/timestep]
      OCCONV(I,J,L) = ( TC0 - CNEW ) * KOC / ( KOC + FREQ )

      !=================================================================
      ! HISTORY (aka netCDF diagnostics)
      !
      ! Archive production of hydrophilic organic carbon (OCPI) from
      ! hydrophobic organic carbon (OCPO) [kg]
      !
      ! NOTE: Consider converting to area-independent units kg/m2/s.
      !=================================================================
      IF ( State_Diag%Archive_ProdOCPIfromOCPO ) THEN
         State_Diag%ProdOCPIfromOCPO(I,J,L) = OcConv(I,J,L)
      ENDIF

      ! Store modified OC concentration back in species array
      TC(I,J,L) = CNEW

   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

 END SUBROUTINE CHEM_OCPO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_ocpi
!
! !DESCRIPTION: Subroutine CHEM\_BCPI calculates dry deposition of
!  hydrophilic OC.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE CHEM_OCPI( Input_Opt, State_Diag, State_Grid, TC, RC )
!
! !USES:
!
   USE ErrCode_Mod
   USE Input_Opt_Mod,   ONLY : OptInput
   USE State_Diag_Mod,  ONLY : DgnState
   USE State_Grid_Mod,  ONLY : GrdState
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt             ! Input Options
   TYPE(GrdState), INTENT(IN)    :: State_Grid            ! Grid State
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(DgnState), INTENT(INOUT) :: State_Diag            ! Diags State
   REAL(fp),       INTENT(INOUT) :: TC(State_Grid%NX,State_Grid%NY,State_Grid%NZ) ! H-philic OC [kg]
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC                    ! Success?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  01 Apr 2004 - R. Park - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Scalars
   INTEGER  :: I,   J,     L
   REAL(fp) :: TC0, CNEW, CCV

   !=================================================================
   ! CHEM_OCPI begins here!
   !=================================================================

   ! Assume success
   RC = GC_SUCCESS

   !$OMP PARALLEL DO       &
   !$OMP DEFAULT( SHARED ) &
   !$OMP PRIVATE( I, J, L, TC0, CCV, CNEW ) &
   !$OMP SCHEDULE( DYNAMIC )
   DO L = 1, State_Grid%NZ
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      ! Initial H-philic OC [kg]
      TC0 = TC(I,J,L)

      ! H-philic OC that used to be H-phobic OC [kg]
      CCV = OCCONV(I,J,L)

      ! Add the amount of converted OCPO to OCPI
      CNEW = TC0 + CCV

      ! Prevent underflow condition
      IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

      ! Store modified concentration back in species array [kg]
      TC(I,J,L) = CNEW

   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   !=================================================================
   ! Zero OCCONV array for next timestep
   !=================================================================
   OCCONV = 0e+0_fp

 END SUBROUTINE CHEM_OCPI
!EOC
#ifdef TOMAS
!-------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aging_carb
!
! !DESCRIPTION: Subroutine AGING\_CARB converts the size-resolved hydrophobic
!  EC or OC to hydrophilic EC or OC with an assumed e-folding time.
!  (win, 9/11/07)
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE AGING_CARB( State_Grid, MIL, MOB )
!
! !USES:
!
   USE State_Grid_Mod, ONLY : GrdState
   USE TIME_MOD,       ONLY : GET_TS_CHEM    ! [=] second
   USE TOMAS_MOD,      ONLY : IBINS
!
! !INPUT PARAMETERS:
!
   TYPE(GrdState), INTENT(IN) :: State_Grid
!
! !INPUT/OUTPUT PARAMETERS:
!
   REAL(fp), INTENT(INOUT) :: MIL(State_Grid%NX,State_Grid%NY, &
                                  State_Grid%NZ, IBINS)
   REAL(fp), INTENT(INOUT) :: MOB(State_Grid%NX,State_Grid%NY, &
                                  State_Grid%NZ, IBINS)
!
! !REMARKS:
!  11 Sep 2007 - W. Trivitayanurak - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER                :: I, J, L, K
   REAL(fp)               :: DTCHEM
   REAL(fp),  PARAMETER   :: TAU_HYDRO = 1.15e+0_fp  ! [=]day

   !=================================================================
   ! AGING_CARB begins here!
   !=================================================================

   DTCHEM = GET_TS_CHEM() / 3600e+0_fp / 24e+0_fp    ![=] day

   DO K = 1, IBINS
      MIL(:,:,:,K) = MIL(:,:,:,K) + &
                     MOB(:,:,:,K) * (1.e+0_fp - DEXP(-DTCHEM/TAU_HYDRO))
      MOB(:,:,:,K) = MOB(:,:,:,K) * (DEXP(-DTCHEM/TAU_HYDRO))
   ENDDO

 END SUBROUTINE AGING_CARB
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: soa_chemistry
!
! !DESCRIPTION: Subroutine SOA\_CHEMISTRY performs SOA formation. This code is
!  from the Caltech group (Hong Liao, Serena Chung, et al) and was modified for
!  GEOS-CHEM. (rjp, bmy, 7/8/04, 12/21/09)
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE SOA_CHEMISTRY( Input_Opt,  State_Chm, State_Diag, &
                           State_Grid, State_Met, RC )
!
! !USES:
!
   USE ErrCode_Mod
   USE ERROR_MOD,      ONLY : DEBUG_MSG
   USE Input_Opt_Mod,  ONLY : OptInput
   USE State_Met_Mod,  ONLY : MetState
   USE State_Chm_Mod,  ONLY : ChmState
   USE State_Diag_Mod, ONLY : DgnState
   USE State_Grid_Mod, ONLY : GrdState
#ifdef APM
   USE APM_INIT_MOD,   ONLY : APMIDS
   USE APM_INIT_MOD,   ONLY : NGCOND,NSO4,NSEA,NBCOC
   USE APM_INIT_MOD,   ONLY : NCTSO4,NCTBC,NCTOC,NCTDST,NCTSEA
#endif
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
   TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
   TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
   TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  Procedure:
!  ============================================================================
!  (1 ) Read in NO3, O3, OH in CHEM_SOA
!  (2 ) Scales these fields using OHNO3TIME in sulfate_mod.f (see GET_OH)
!  (3 ) Calculate reaction rates (Serena's OCHEMPARAETER)
!  (4 ) CALCULATE DELHC
!  (5 ) get T0M gas products
!  (6 ) equilibrium calculation
!                                                                             .
!  As of 5/20/10: Havala's New formulation
!                                                                             .
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  % FOR SEMIVOLATILE POA and IVOC (aka SOA_SVPOA) simulations:      %
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  GEOS-Chem treats formation of aerosol from 11 parent hydrocarbons
!  and oxidation by OH, O3, and NO3:
!                                                                             .
!  The parent hydrocarbons are lumped into 5 semivolatile systems:
!  TSOA/G: the lumped semivolatile oxidation products of
!          monoterpene and sesquiterpene oxidation
!  ISOA/G: the lumped semivolatile oxidation products of isoprene ox
!           (REMOVED IN GEOS-Chem 12.6.0, July 2019)
!  ASOA/G: the lumped semivolatile (and nonvolatile) products of
!          benzene, toluene, xylene, and naphthalene (IVOC surrogate)
!          oxidation
!  POA/G : the lumped primary semivolatile emissions
!  OPOA/G: the lumped products of primary semivolatile oxidation
!                                                                             .
!  Parent HC      Oxidized by       Products
!  =============  ================  ==================================
!  MTPA           OH, O3, NO3       TSOA/G0-3
!  LIMO           OH, O3, NO3       TSOA/G1-3
!  MTPO           OH, O3, NO3       TSOA/G0-3
!  SESQ           OH, O3, NO3       TSOA/G1-3
!  ISOP           OH, NO3           ISOA/G1-3 (REMOVED IN 12.6.0, Jul 2019)
!  BENZ           OH, (+NO,HO2)     ASOAN, ASOA/G1-3
!  TOLU           OH, (+NO,HO2)     ASOAN, ASOA/G1-3
!  XYLE           OH, (+NO,HO2)     ASOAN, ASOA/G1-3
!  SVOC/POA       OH                POA/G1-2
!  O-SVOC/OPOA    OH                OPOA/G1-2
!  NAP            OH, (+NO,HO2)     ASOAN, ASOA/G1-3
!                                                                             .
!  Species that must be defined in input.geos (in addition to standard
!  full chem species) (34 additional):
!  TSOA1      TSOG1      ASOA1      ASOG1
!  TSOA2      TSOG2      ASOA2      ASOG2
!  TSOA3      TSOG3      ASOA3      ASOG3
!  ASOAN      TSOA0      TSOG0
!  BENZ       TOLU       XYLE       MTPA       LIMO       MTPO
!  NAP
!  POA1       POG1       POA2       POG2
!  OPOA1      OPOG1      OPOA2      OPOG2
!                                                                             .
!  The following should NOT be defined for semivol POA: OCPI, OCPO
!                                                                             .
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  % FOR NON-VOLATILE TRADITIONAL POA (aka SOA) simulations:         %
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  GEOS-Chem treats formation of aerosol from 8 parent hydrocarbons
!  and oxidation by OH, O3, and NO3:
!                                                                             .
!  Two non-volatile,traditional primary OC species exist:
!  OCPO: hydrophobic POA
!  OCPI: hydrophillic POA
!                                                                             .
!  The parent hydrocarbons are lumped into 3 semivolatile systems:
!  TSOA/G: the lumped semivolatile oxidation products of
!          monoterpene and sesquiterpene oxidation
!  ISOA/G: the lumped semivolatile oxidation products of isoprene ox
!           (REMOVED IN GEOS-Chem 12.6.0, July 2019)
!  ASOA/G: the lumped semivolatile (and nonvolatile) products of
!          benzene, toluene, and xylene oxidation
!                                                                             .
!  Parent HC      Oxidized by       Products
!  =============  ================  ==================================
!  MTPA           OH, O3, NO3       TSOA/G0-3
!  LIMO           OH, O3, NO3       TSOA/G1-3
!  MTPO           OH, O3, NO3       TSOA/G0-3
!  SESQ           OH, O3, NO3       TSOA/G1-3
!  ISOP           OH, NO3           ISOA/G1-3 (REMOVED IN 12.6.0, Jul 2019)
!  BENZ           OH, (+NO,HO2)     ASOAN, ASOA/G1-3
!  TOLU           OH, (+NO,HO2)     ASOAN, ASOA/G1-3
!  XYLE           OH, (+NO,HO2)     ASOAN, ASOA/G1-3
!                                                                             .
!  Species that must be defined in input.geos (in addition to standard
!  full chem species) (25 additional):
!  TSOA1      TSOG1      ASOA1      ASOG1
!  TSOA2      TSOG2      ASOA2      ASOG2
!  TSOA3      TSOG3      ASOA3      ASOG3
!  ASOAN      TSOA0      TSOG0
!  BENZ       TOLU       XYLE       MTPA       LIMO       MTPO
!                                                                             .
!  The following should NOT be defined for traditional POA:
!     NAP, POA/G OPOA/G
!                                                                             .
!  References (see above for full citations):
!  ===========================================================================
!  Monoterpenes and sesquiterpenes:
!    Experimental Data:
!      Griffin et al. 1999 JGR      (sesquiterps low NOx)
!      Shilling et al. 2008 ACP     (a-pinene ozonolysis for MTPO/MTPA)
!      Zhang et al. 2006 JPhysChemA (limonene ozonolysis)
!      Ng et al. 2007 ACP           (data for NOx effect on sesq aerosol)
!    Modeling:
!      Chung and Seinfeld 2002 JGR  (original formulation in GEOS-Chem)
!      Liao et al. 2007 JGR         (comparison to measurements)
!      Pye et al. in prep 2010      (new lumping scheme, NOx effect)
!  Isoprene
!      Kroll et al. 2006 ES&T       (low NOx experiments)
!      Ng et al. 2008 ACP           (isoprene + NO3 experiments)
!      Henze et al. 2006 GRL        (low NOx isoprene modeling in GEOS-Chem)
!      Pye et al. in prep 2010      (new lumping scheme and isop+no3 modeling)
!  Aromatics: benz, tolu, xyle
!      Ng et al. 2007 ACP           (experiments)
!      Henze et al. 2008 ACP        (global modeling)
!  POA/OPOA
!      Shrivastava et al. 2006 ES&T (POA experiments)
!      Grieshop et al. 2009 ACP     (POA/SVOC oxidation experiments)
!      Pye and Seinfeld 2010 ACP    (global modeling)
!  IVOC/Naphthalene
!      Chan et al. 2009 ACP         (experiments)
!      Pye and Seinfeld 2010 ACP    (global modeling)
!
! !REVISION HISTORY:
!  08 Jul 2004 - R. Park - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   LOGICAL, SAVE   :: FIRST = .TRUE.
   LOGICAL         :: prtDebug
   INTEGER         :: I,        J,        L,        N
   ! no more NOx (hotp 5/22/10)
   INTEGER         :: JHC,      IPR!,      NOX ! (dkh, 10/30/06)
   INTEGER         :: JSV ! (hotp 5/13/10)
   INTEGER         :: JSVPOA, JSVOPOA ! for diag (hotp 5/17/10)
   REAL(fp)        :: RTEMP,    VOL,      FAC,      MPOC
   REAL(fp)        :: MNEW,     MSOA_OLD, MPRODUCT, CAIR
   REAL(fp)        :: LOWER,    UPPER,    TOL,      VALUE
   REAL(fp)        :: KO3(MHC), KOH(MHC), KNO3(MHC)
   REAL(fp)        :: KOM(MPROD,MSV)
   REAL(fp)        :: GM0(MPROD,MSV), AM0(MPROD,MSV)
   REAL(fp)        :: ORG_AER(MPROD,MSV)
   REAL(fp)        :: ORG_GAS(MPROD,MSV)

   ! Parent HC reacted diag for arom/IVOC (hotp 5/17/10)
   REAL(fp)        :: DARO2_TOT_0(8)             ! added NAP (hotp 7/22/09)
   REAL(fp), SAVE  :: DARO2_TOT(8) = 0e+0_fp     ! added NAP (hotp 7/22/09)

   ! Rate constant for RO2+HO2, RO2+NO rxns (hotp 5/7/10)
   REAL(fp)          :: KRO2NO, KRO2HO2

   ! semivolpoa2: add diagnostic info for POA (hotp 3/11/09)
   ! semivolpoa4: add dimension for OPOA (hotp 3/27/09)
   REAL(fp)        :: GLOB_AM0_POA(State_Grid%NX,State_Grid%NY,State_Grid%NZ, &
                                   MNOX,MPROD,2)
   REAL(fp)        :: GLOB_AM0_POA_0(State_Grid%NX,State_Grid%NY,State_Grid%NZ,&
                                     MNOX,MPROD,2)
   REAL(fp), SAVE  :: GLOB_POA_PROD
   REAL(fp), SAVE  :: GLOB_OPOA_PROD

   ! Update IEPOX aerosol-phase species concentrations:
   INTEGER         :: SPIND
   REAL(fp)        :: fORGS,   fORGN
   REAL(fp)        :: DELDIOL, DELORGS, DELORGN

   ! Pointers
   REAL(fp), POINTER :: Spc(:,:,:,:)

   ! Debug
   Integer           :: IIDebug, JJDebug

#ifdef APM
   INTEGER           :: IFINORG
   INTEGER           :: NTEMP
#endif

   !=================================================================
   ! SOA_CHEMISTRY begins here!
   !=================================================================

   ! Point to chemical species array [kg]
   Spc          => State_Chm%Species

   ! Do we have to print debug output?
   prtDebug     = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

   ! Zero some diagnostics (hotp 5/17/10)
   GLOB_POGRXN  = 0e+0_fp
   OAGINITSAVE  = 0e+0_fp
   DELTASOGSAVE = 0e+0_fp

   ! Save initial OA+OG for diagnostic (hotp 5/17/10)
   IF ( prtDebug ) THEN
      CALL SAVE_OAGINIT( State_Chm, State_Grid, State_Met )
   ENDIF

   !this section is need to initialize the equalibirum constants and
   !yield parameters, (jje 06/22/10)
   IF ( FIRST ) THEN

      ! initialize alphas, Koms and rxn rate constants
      CALL SOA_PARA_INIT( Input_Opt )

      ! Diagnostic/debug info (hotp 5/22/10)
      IF ( prtDebug ) THEN
         WRITE(*,*) 'HC and SV IDs'
         print*, 'Monoterpenes and sesquiterpenes'
         print*, 'MTPA ', PARENTMTPA, IDSV(PARENTMTPA)
         print*, 'LIMO ', PARENTLIMO, IDSV(PARENTLIMO)
         print*, 'MTPO ', PARENTMTPO, IDSV(PARENTMTPO)
         print*, 'SESQ ', PARENTSESQ, IDSV(PARENTSESQ)
         print*, 'Isoprene'
         print*, 'ISOP ', PARENTISOP, IDSV(PARENTISOP)
         print*, 'Arom/IVOC'
         print*, 'BENZ ', PARENTBENZ, IDSV(PARENTBENZ)
         print*, 'TOLU ', PARENTTOLU, IDSV(PARENTTOLU)
         print*, 'XYLE ', PARENTXYLE, IDSV(PARENTXYLE)
         print*, 'NAP  ', PARENTNAP,  IDSV(PARENTNAP )
         print*, 'POA and OPOA'
         print*, 'POA  ', PARENTPOA,  IDSV(PARENTPOA)
         print*, 'OPOA ', PARENTOPOA, IDSV(PARENTOPOA)
         print*,'NPROD and NNOX for SV groups'
         DO JSV = 1, MSV
            print*, (JSV),NPROD(JSV),NNOX(JSV)
         ENDDO
      ENDIF

      ! Zero diagnostic (hotp 5/24/10)
      GLOB_POA_PROD  = 0e+0_fp
      GLOB_OPOA_PROD = 0e+0_fp
      DELTAHCSAVE    = 0e+0_fp
      SPECSOAPROD    = 0e+0_fp
      SPECSOAEVAP    = 0e+0_fp

      FIRST = .FALSE.
   ENDIF

   IF ( prtDebug ) THEN
      print*, ' MAX DARO2 = ', MAXLOC(GLOB_DARO2(:,:,:,1,1)), &
                               MAXVAL(GLOB_DARO2(:,:,:,1,1))
   ENDIF

   ! semivolpoa2: diagnostic info (hotp 3/11/09)
   GLOB_AM0_POA   = 0.e+0_fp
   GLOB_AM0_POA_0 = 0.e+0_fp

   ! added debug (hotp 8/24/09)
   IF ( prtDebug ) THEN
      IIDebug = Min(State_Grid%NX,20)
      JJDebug = Min(State_Grid%NY,33)
      print*, 'START SOA_CHEMISTRY'
      print*, 'species   ','global sum   ', &
              'box ',IIDEBUG,',',JJDEBUG,',2     ', &
              'box ',IIDEBUG,',',JJDEBUG,',10    '

      IF ( Input_Opt%LSVPOA ) THEN
         print*,'POA1', SUM(Spc(:,:,:,id_POA1))
         print*,'POA2', SUM(Spc(:,:,:,id_POA2))
         print*,'POG1', SUM(Spc(:,:,:,id_POG1))
         print*,'POG2', SUM(Spc(:,:,:,id_POG2))
         print*,'OPOA1',SUM(Spc(:,:,:,id_OPOA1))
         print*,'OPOA2',SUM(Spc(:,:,:,id_OPOA2))
         print*,'OPOG1',SUM(Spc(:,:,:,id_OPOG1))
         print*,'OPOG2',SUM(Spc(:,:,:,id_OPOG2))
      ENDIF

      ! semivolpoa3: debug POAEMISS
      print*, 'POAEMISS,1,POG1', &
               SUM(POAEMISS(:,:,:,1)),POAEMISS(IIDEBUG,JJDEBUG,2,1), &
                   POAEMISS(IIDEBUG,JJDEBUG,10,1)
      print*, 'POAEMISS,2,POG2', &
               SUM(POAEMISS(:,:,:,2)),POAEMISS(IIDEBUG,JJDEBUG,2,2), &
                   POAEMISS(IIDEBUG,JJDEBUG,10,2)

      ! 10/12/09 debug
      IF ( Input_Opt%LSVPOA ) THEN
         print*,'POAG tot', SUM(Spc(:,:,:,id_POA1))+ &
                            SUM(Spc(:,:,:,id_POA2))+ &
                            SUM(Spc(:,:,:,id_POG1))+ &
                            SUM(Spc(:,:,:,id_POG2))
      ENDIF
   ENDIF

   ! Locate POA and OPOA in AM0 (hotp 5/17/10)
   JSVPOA  = IDSV(PARENTPOA)
   JSVOPOA = IDSV(PARENTOPOA)

   ! For parallel do:
   ! add KRO2XXX to private (hotp 5/7/10)
   ! add JSV to private (hotp 5/13/10)
   ! remove NOX (hotp 5/22/10)
   !$OMP PARALLEL DO       &
   !$OMP DEFAULT( SHARED ) &
   !$OMP PRIVATE( I,        J,        L,     JHC,   IPR,   GM0,  AM0  ) &
#ifdef APM
   !$OMP PRIVATE( N, NTEMP) &
#endif
   !$OMP PRIVATE( VOL,      FAC,      RTEMP, KO3,   KOH,   KNO3, CAIR ) &
   !$OMP PRIVATE( MPRODUCT, MSOA_OLD, VALUE, UPPER, LOWER, MNEW, TOL  ) &
   !$OMP PRIVATE( ORG_AER,  ORG_GAS,  KOM,   MPOC                     ) &
   !$OMP PRIVATE( KRO2NO,   KRO2HO2,  JSV                             )
   DO L = 1, State_Grid%MaxChemLev
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      ! Skip non-chemistry boxes
      IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

      ! Volume of grid box [m3]
      VOL    = State_Met%AIRVOL(I,J,L)

      ! conversion factor from kg to ug/m3
      FAC    = 1.e+9_fp / VOL

      ! air conc. in kg/m3
      CAIR   = State_Met%AD(I,J,L) / VOL

      ! Temperature [K]
      RTEMP  = State_Met%T(I,J,L)

      ! Get SOA yield parameters
      ! ALPHA is a module variable now. (ccc, 2/2/10)
      ! add arguments for RO2+NO, RO2+HO2 rates (hotp 5/7/10)
      CALL SOA_PARA( RTEMP, KO3, KOH, KNO3, KOM, &
                     I,     J,   L,   KRO2NO, KRO2HO2, State_Met )

      ! Partition mass of gas & aerosol species
      ! according to 5 VOC classes & 3 oxidants
      CALL SOA_PARTITION( I, J, L, GM0, AM0, State_Chm )

      ! hotp diagnostic (3/11/09)
      GLOB_AM0_POA_0(I,J,L,1,:,:) = AM0(:,JSVPOA:JSVOPOA)

      ! Compute oxidation of hydrocarbons by O3, OH, NO3
      ! ALPHA is a module variable now (ccc, 2/2/10)
      ! semivolpoa2: emit POA into semivolatiles here (hotp 2/27/09)
      ! add RO2+NO,HO2 rate constants (hotp 5/7/10)
      CALL CHEM_NVOC( I,          J,         L,          &
                      KO3,        KOH,       KNO3,       &
                      GM0,        KRO2NO,    KRO2HO2,    &
                      Input_Opt,  State_Chm, State_Diag, &
                      State_Grid, State_Met, RC          )

      !==============================================================
      ! Equilibrium calculation between GAS (SOG) and Aerosol (SOA)
      !==============================================================

      ! Initialize other arrays to be safe  (dkh, 11/10/06)
      ! update dims (hotp 5/22/10)
      ORG_AER(:,:) = 0e+0_fp
      ORG_GAS(:,:) = 0e+0_fp

      ! Individual SOA's: convert from [kg] to [ug/m3] or [kgC] to [ugC/m3]
      DO JSV = 1, MAXSIMSV
      DO IPR = 1, NPROD(JSV)
         ORG_GAS(IPR,JSV) = GM0(IPR,JSV) * FAC
         ORG_AER(IPR,JSV) = AM0(IPR,JSV) * FAC
      ENDDO
      ENDDO

      ! semivolpoa2: include OA mass with POA (hotp 3/2/09)
      ! Check to make sure POA is defined (hotp 8/24/09)
      ! Convert from [ugC/m3] to [ug/m3]
      IF ( id_POA1 > 0 ) THEN
         JHC = PARENTPOA
         JSV = IDSV(JHC)
         DO IPR = 1, NPROD(JSV)
            ORG_GAS(IPR,JSV) = ORG_GAS(IPR,JSV) * OCFPOA(I,J)
            ORG_AER(IPR,JSV) = ORG_AER(IPR,JSV) * OCFPOA(I,J)
         ENDDO
      ENDIF

      ! semivolpoa4opoa: add OPOA mass (hotp 3/18/09)
      ! Check to make sure OPOA is defined (hotp 8/24/09)
      ! Convert from [ugC/m3] to [ug/m3]
      IF ( id_OPOA1 > 0 ) THEN
         JHC = PARENTOPOA
         JSV = IDSV(JHC)
         DO IPR = 1, NPROD(JSV)
            ORG_GAS(IPR,JSV) = ORG_GAS(IPR,JSV) * OCFOPOA(I,J)
            ORG_AER(IPR,JSV) = ORG_AER(IPR,JSV) * OCFOPOA(I,J)
         ENDDO
      ENDIF

      !-----------------------------------------------------------
      ! Compute SOG condensation onto OC aerosol
      !
      ! Primary organic aerosol concentrations [ug/m3]
      ! We carry carbon mass only in the Spc array and here
      ! multiply by 2.1 to account for non-carbon mass in the SOA.
      !
      ! Partitioning theory (Pankow, 1994) describes organic
      ! phase partitioning assuming absorption into pre-existing
      ! organic mass.  There is currently no theoretical or
      ! laboratory support for absorption of organics into
      ! inorganics.
      !
      ! Note that previous versions of the standard code
      ! (v7-04-07 through v8-02-04) did include absorption into
      ! inorganics.
      !
      ! (Colette Heald, 12/3/09)
      !-----------------------------------------------------------
#ifdef APM
      MPOC = FAC * SUM(Spc(I,J,L,APMIDS%id_OCBIN1:(APMIDS%id_OCBIN1+NBCOC-1)))
      MPOC = MPOC * 2.1d0

      IFINORG = 2
      IF(IFINORG.EQ.1) THEN !Yu+  consider inorg in the SOA partition

         IF ( APMIDS%id_SO4 > 0 .and. &
              APMIDS%id_NH4 > 0 .and. &
              APMIDS%id_NIT > 0 ) THEN
            ! Then compute SOG condensation onto SO4, NH4, NIT aerosols
            MPOC = MPOC + ( Spc(I,J,L,APMIDS%id_NH4) + &
                            Spc(I,J,L,APMIDS%id_NIT) ) * FAC

            IF(NSO4>=1)THEN
               NTEMP=APMIDS%id_SO4BIN1-1
               DO N=(NTEMP+1),(NTEMP+NSO4)
                  MPOC = MPOC + Spc(I,J,L,N) * FAC
               ENDDO
            ENDIF

            IF(NCTSO4>=1)THEN
               NTEMP=APMIDS%id_CTSO4-1
               DO N=(NTEMP+1),(NTEMP+NCTSO4)
                  MPOC = MPOC + Spc(I,J,L,N) * FAC
               ENDDO
            ENDIF

            IF(NCTBC>=1)THEN
               NTEMP=APMIDS%id_CTBC-1
               DO N=(NTEMP+1),(NTEMP+NCTBC)
                  MPOC = MPOC + Spc(I,J,L,N) * FAC
               ENDDO
            ENDIF

            IF(NCTOC>=1)THEN
               NTEMP=APMIDS%id_CTOC-1
               DO N=(NTEMP+1),(NTEMP+NCTOC)
                  MPOC = MPOC + Spc(I,J,L,N) * FAC
               ENDDO
            ENDIF

            IF(NCTSEA>=1)THEN
               NTEMP=APMIDS%id_CTSEA-1
               DO N=(NTEMP+1),(NTEMP+NCTSEA)
                  MPOC = MPOC + Spc(I,J,L,N) * FAC
               ENDDO
            ENDIF

            IF(NSEA>=1)THEN
               NTEMP=APMIDS%id_SEABIN1-1
               DO N=(NTEMP+1),(NTEMP+16) ! SALTbin16 = 1 um
                  MPOC = MPOC + Spc(I,J,L,N) * FAC
               ENDDO
            ENDIF

            !Add MSA
            MPOC = MPOC + Spc(I,J,L,APMIDS%id_MSA) * FAC

         ENDIF

      ELSEIF(IFINORG.EQ.2) THEN !Consider SV SOA partition on LV SOA

         MPOC = MPOC + FAC * (Spc(I,J,L,APMIDS%id_CTSO4) + & !MSULFLV
                Spc(I,J,L,APMIDS%id_CTBC+1)  + & !MBCLV
                Spc(I,J,L,APMIDS%id_CTOC+1)  + & !MOCLV
                Spc(I,J,L,APMIDS%id_CTDST+1) + & !MDSTLV
                Spc(I,J,L,APMIDS%id_CTSEA+1))    !MSALTLV

      ELSE

         ! Compute SOG condensation onto OC aerosol
         MPOC = ( Spc(I,J,L,id_OCPI) + Spc(I,J,L,id_OCPO) ) * FAC
         MPOC = MPOC * 2.1d0

      ENDIF
#else
      ! Now treat either traditional POA or semivolatile POA (hotp 7/25/10)
      IF ( id_OCPI > 0 .and. id_OCPO > 0 ) THEN
         MPOC = ( Spc(I,J,L,id_OCPI) + Spc(I,J,L,id_OCPO) ) * FAC
         MPOC = MPOC * OCFOPOA(I,J)
      ELSE
         ! semivolpoa2: MPOC is zero now (hotp 2/27/09)
         MPOC = 1e-30_fp
      ENDIF
#endif

      !==============================================================
      ! Solve for MNEW by solving for SOA=0
      !==============================================================
      IF ( ( MPOC / ( CAIR*1.e+9_fp ) ) <= 2.1e-18_fp ) THEN
         VALUE = 0.e+0_fp
         UPPER = 0.e+0_fp

         ! Now use SV (hotp 5/13/10)
         ! update dims (hotp 5/22/10)
         DO JSV = 1, MAXSIMSV
         DO IPR = 1, NPROD(JSV)
            VALUE = VALUE + KOM(IPR,JSV) * (ORG_GAS(IPR,JSV) + ORG_AER(IPR,JSV))
            UPPER = UPPER + ORG_GAS(IPR,JSV) + ORG_AER(IPR,JSV)
         ENDDO
         ENDDO

         IF ( VALUE <= 1.e+0_fp ) THEN
            MNEW  = 0.e+0_fp
         ELSE
            LOWER = 1.e-18_fp * ( CAIR * 1.e+9_fp )
            TOL   = 1.e-18_fp
            MNEW  = ZEROIN(LOWER,UPPER,TOL,MPOC,ORG_AER,ORG_GAS,KOM)
         ENDIF

      ELSE

         UPPER = MPOC

         ! Now use SV (hotp 5/13/10)
         ! update dims (hotp 5/22/10)
         DO JSV = 1, MAXSIMSV
         DO IPR = 1, NPROD(JSV)
            UPPER = UPPER + ORG_GAS(IPR,JSV) + ORG_AER(IPR,JSV)
         ENDDO
         ENDDO

         LOWER = MPOC
         TOL   = 1.e-9_fp*MPOC
         MNEW  = ZEROIN(LOWER,UPPER,TOL,MPOC,ORG_AER,ORG_GAS,KOM)

      ENDIF

      !==============================================================
      ! Equilibrium partitioning into new gas and aerosol
      ! concentrations for individual contributions of SOA
      !==============================================================
      IF ( MNEW > 0.e+0_fp ) THEN

         ! Use actual number of HC (hotp 8/24/09)
         ! Now use SV (hotp 5/13/10)
         ! updated dims (hotp 7/28/1)
         DO JSV = 1, MAXSIMSV
         DO IPR = 1, NPROD(JSV)
            ORG_AER(IPR,JSV) = KOM(IPR,JSV)*MNEW / &
                               (1.e+0_fp + KOM(IPR,JSV) * MNEW ) * &
                               (ORG_AER(IPR,JSV) + ORG_GAS(IPR,JSV))

            IF ( KOM(IPR,JSV).NE.0e+0_fp ) THEN
               ORG_GAS(IPR,JSV) = ORG_AER(IPR,JSV) * 1.e+8_fp / &
                                  ( KOM(IPR,JSV) * MNEW * 1.e+8_fp)
            ELSE
               ORG_GAS(IPR,JSV) = 0.e+0_fp
            ENDIF

         ENDDO
         ENDDO

         ! semivolpoa2: remove OA mass from POA (hotp 3/2/09)
         ! Check if POA defined (hotp 8/24/09)
         IF ( id_POA1 > 0 ) THEN
            JHC = PARENTPOA
            JSV = IDSV(JHC)
            DO IPR = 1, NPROD(JSV)
               ORG_GAS(IPR,JSV) = ORG_GAS(IPR,JSV) / OCFPOA(I,J)
               ORG_AER(IPR,JSV) = ORG_AER(IPR,JSV) / OCFPOA(I,J)
            ENDDO
         ENDIF

         ! semivolpoa4opoa: remove OA mass from OPOA (hotp 3/18/09)
         ! Check if OPOA defined (hotp 8/24/09)
         IF ( id_OPOA1 > 0 ) THEN
            JHC = PARENTOPOA
            JSV = IDSV(JHC)
            DO IPR = 1, NPROD(JSV)
               ORG_GAS(IPR,JSV) = ORG_GAS(IPR,JSV) / OCFOPOA(I,J)
               ORG_AER(IPR,JSV) = ORG_AER(IPR,JSV) / OCFOPOA(I,J)
            ENDDO
         ENDIF

         ! STORE PRODUCT INTO T0M
         ! Use actual number of HC for sim (hotp 8/24/09)
         ! change to SV (hotp 5/13/10)
         ! update dims (hotp 5/22/10)
         DO JSV = 1, MAXSIMSV
         DO IPR = 1, NPROD(JSV)
            GM0(IPR,JSV) = ORG_GAS(IPR,JSV) / FAC
            AM0(IPR,JSV) = ORG_AER(IPR,JSV) / FAC
         ENDDO
         ENDDO

      !==============================================================
      ! Mnew=0.e+0_fp, all SOA evaporates to the gas-phase
      !==============================================================
      ELSE

         ! Use actual number of HC for sim (hotp 8/24/09)
         ! Change to SV (hotp 5/13/10)
         DO JSV = 1, MAXSIMSV
         DO IPR = 1, NPROD(JSV)
            GM0(IPR,JSV) = GM0(IPR,JSV) + AM0(IPR,JSV)
            !AM0(IPR,JSV) = 1.D-18 * State_Met%AD(I,J,L)
            ! try this to fix MB problem (hotp 5/25/10)
            AM0(IPR,JSV) = 1.e-20_fp
         ENDDO
         ENDDO

      ENDIF

      ! enforce direct yield for low nox aromatics
      ! no longer loop (hotp 7/28/10)
      JHC = PARENTBENZ
      JSV = IDSV(JHC)
      !DO IPR = 1, NPROD(JSV)
      IPR = 4 ! HARDWIRED!!!!!!!!!
      AM0(IPR,JSV) = AM0(IPR,JSV) + GM0(IPR,JSV)
      GM0(IPR,JSV) = 0e+0_fp

      ! Lump SOA
      CALL SOA_LUMP( I, J, L, GM0, AM0, State_Chm, State_Diag )

      ! hotp diagnostic (3/11/09)
      GLOB_AM0_POA(I,J,L,1,:,:) = AM0(:,JSVPOA:JSVOPOA)

      ! Check equilibrium (hotp 5/18/10)
      IF ( prtDebug ) THEN
         ! IDSV for lumped arom/IVOC is hardwired (=3) (hotp 5/20/10)
         ! Low NOX (non-volatile) aromatic product is IPR=4
         CALL CHECK_EQLB( I, J, L, KOM, FAC, MNEW, LOWER, TOL, &
                          ORG_GAS(4,3), ORG_AER(4,3), MPOC, State_Chm )
      ENDIF

   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Debug: check mass balance (hotp 5/18/10)
   IF ( prtDebug ) THEN
      CALL CHECK_MB( Input_Opt,State_Chm, State_Grid, State_Met )
   ENDIF

   !------------------------------------------------------------------------
   !### Now only print when ND70 is turned on (bmy, 4/21/10)
   IF ( prtDebug ) THEN

      IF ( Input_Opt%LSVPOA ) THEN
         ! 10/12/09 debug
         print*,'POAG tot', SUM(Spc(:,:,:,id_POA1))+ &
                            SUM(Spc(:,:,:,id_POA2))+ &
                            SUM(Spc(:,:,:,id_POG1))+ &
                            SUM(Spc(:,:,:,id_POG2))
      ENDIF

      ! dkh print some diagnostics
      ! Parent hydrocarbon reacted diagnostic (hotp 5/17/10)
      DARO2_TOT_0(:)    = DARO2_TOT(:)

      print*, ' MAX DARO2 = ', MAXLOC(GLOB_DARO2(:,:,:,1,1)), &
                               MAXVAL(GLOB_DARO2(:,:,:,1,1))

      DARO2_TOT(1) = DARO2_TOT(1) + SUM(GLOB_DARO2(:,:,:,1,1)) / 1d9
      DARO2_TOT(2) = DARO2_TOT(2) + SUM(GLOB_DARO2(:,:,:,2,1)) / 1d9
      DARO2_TOT(3) = DARO2_TOT(3) + SUM(GLOB_DARO2(:,:,:,1,2)) / 1d9
      DARO2_TOT(4) = DARO2_TOT(4) + SUM(GLOB_DARO2(:,:,:,2,2)) / 1d9
      DARO2_TOT(5) = DARO2_TOT(5) + SUM(GLOB_DARO2(:,:,:,1,3)) / 1d9
      DARO2_TOT(6) = DARO2_TOT(6) + SUM(GLOB_DARO2(:,:,:,2,3)) / 1d9
      ! added NAP diagnostic info (hotp 7/22/09)
      ! amount of RO2 reacted in high and low NOx pathways
      DARO2_TOT(7) = DARO2_TOT(7) + SUM(GLOB_DARO2(:,:,:,1,4)) / 1d9
      DARO2_TOT(8) = DARO2_TOT(8) + SUM(GLOB_DARO2(:,:,:,2,4)) / 1d9

      ! DARO2 is not mass of parent HC, not RO2 (hotp 5/14/10)
      print*,'Accumulated parent HC reacted to RO2H,N products in Tg'
      print*,'Units are Tg of parent'
      print*, 'GLOB_DBRO2 NOX =', DARO2_TOT(1), (DARO2_TOT(1) - DARO2_TOT_0(1))
      print*, 'GLOB_DBRO2 HO2 =', DARO2_TOT(2), (DARO2_TOT(2) - DARO2_TOT_0(2))
      print*, 'GLOB_DTRO2 NOX =', DARO2_TOT(3), (DARO2_TOT(3) - DARO2_TOT_0(3))
      print*, 'GLOB_DTRO2 HO2 =', DARO2_TOT(4), (DARO2_TOT(4) - DARO2_TOT_0(4))
      print*, 'GLOB_DXRO2 NOX =', DARO2_TOT(5), (DARO2_TOT(5) - DARO2_TOT_0(5))
      print*, 'GLOB_DXRO2 HO2 =', DARO2_TOT(6), (DARO2_TOT(6) - DARO2_TOT_0(6))
      print*, 'GL_DAR NOX NAP =', DARO2_TOT(7), (DARO2_TOT(7) - DARO2_TOT_0(7))
      print*, 'GL_DAR HO2 NAP =', DARO2_TOT(8), (DARO2_TOT(8) - DARO2_TOT_0(8))
      ! end arom parent HC reacted diag (hotp 5/17/10)

      ! semivolpoa2: diagnostic info (hotp 3/11/09)
      ! initial info
      print*, 'AFTER SOA_CHEMISTRY'
      print*, 'species   ','global sum   ', &
              'box ',IIDEBUG,',',JJDEBUG,',2     ', &
              'box ',IIDEBUG,',',JJDEBUG,',10    '

      IF ( Input_Opt%LSVPOA ) THEN
         print*,'POA1 ', SUM(Spc(:,:,:,id_POA1))
         print*,'POA2 ', SUM(Spc(:,:,:,id_POA2))
         print*,'POG1 ', SUM(Spc(:,:,:,id_POG1))
         print*,'POG2 ', SUM(Spc(:,:,:,id_POG2))
         print*,'OPOA1', SUM(Spc(:,:,:,id_OPOA1))
         print*,'OPOA2', SUM(Spc(:,:,:,id_OPOA2))
         print*,'OPOG1', SUM(Spc(:,:,:,id_OPOG1))
         print*,'OPOG2', SUM(Spc(:,:,:,id_OPOG2))
      ENDIF

      ! semivolpoa4: diag for POG reacted (hotp 3/27/09)
      print*, 'POGRXN ', SUM(GLOB_POGRXN(:,:,:,:)), &
                         SUM( GLOB_POGRXN(IIDEBUG,JJDEBUG,2,:)), &
                         SUM(GLOB_POGRXN(IIDEBUG,JJDEBUG,10,:))
      ! POGRXN (hotp 10/11/09)
      print*, 'POGRXN p1', SUM(GLOB_POGRXN(:,:,:,1))
      print*, 'POGRXN p2', SUM(GLOB_POGRXN(:,:,:,2))

      ! semivolpoa3: POG1 + POG2
      print*, 'POAEMISS tot', SUM(POAEMISS(:,:,:,:)), &
               POAEMISS(IIDEBUG,JJDEBUG,2,1)+POAEMISS(IIDEBUG,JJDEBUG,2,2), &
               POAEMISS(IIDEBUG,JJDEBUG,10,1)+POAEMISS(IIDEBUG,JJDEBUG,10,2)

      ! semivolpoa4opoa: add OPOA (hotp 3/27/09)
      IF ( id_OPOA1 > 0 ) THEN ! hotp 8/24/09 ! hotp 10/11/09
         GLOB_POA_PROD  = GLOB_POA_PROD + &
                          SUM(GLOB_AM0_POA(:,:,:,:,:,1)) - &
                          SUM(GLOB_AM0_POA_0(:,:,:,:,:,1))
         GLOB_OPOA_PROD = GLOB_OPOA_PROD + &
                          SUM(GLOB_AM0_POA(:,:,:,:,:,2)) - &
                          SUM(GLOB_AM0_POA_0(:,:,:,:,:,2))

         ! semivolpoa4opoa: diagnostic info (hotp 3/27/09)
         print*, 'POA produced (cumulative tp date) ', GLOB_POA_PROD
         print*, 'POA P1',  SUM(GLOB_AM0_POA(:,:,:,1,1,1)) - &
                            SUM(GLOB_AM0_POA_0(:,:,:,1,1,1))
         print*, 'POA P2',  SUM(GLOB_AM0_POA(:,:,:,1,2,1)) - &
                            SUM(GLOB_AM0_POA_0(:,:,:,1,2,1))
         print*, 'OPOA P1', SUM(GLOB_AM0_POA(:,:,:,1,1,2)) - &
                            SUM(GLOB_AM0_POA_0(:,:,:,1,1,2))
         print*, 'OPOA P2', SUM(GLOB_AM0_POA(:,:,:,1,2,2)) - &
                            SUM(GLOB_AM0_POA_0(:,:,:,1,2,2))

         print*, 'OPOA produced (cumulative to date) ',GLOB_OPOA_PROD
         print*, 'POA products '
         print*, 'product 1', SUM(GLOB_AM0_POA(:,:,:,1,1,1)), &
                              GLOB_AM0_POA(IIDEBUG,JJDEBUG,2,1,1,1), &
                              GLOB_AM0_POA(IIDEBUG,JJDEBUG,10,1,1,1)
         print*, 'product 2', SUM(GLOB_AM0_POA(:,:,:,1,2,1)), &
                              GLOB_AM0_POA(IIDEBUG,JJDEBUG,2,1,2,1), &
                              GLOB_AM0_POA(IIDEBUG,JJDEBUG,10,1,2,1)
         print*, 'OPOA products'
         print*, 'product 1', SUM(GLOB_AM0_POA(:,:,:,1,1,2)), &
                              GLOB_AM0_POA(IIDEBUG,JJDEBUG,2,1,1,2), &
                              GLOB_AM0_POA(IIDEBUG,JJDEBUG,10,1,1,2)
         print*, 'product 2', SUM(GLOB_AM0_POA(:,:,:,1,2,2)), &
                              GLOB_AM0_POA(IIDEBUG,JJDEBUG,2,1,2,2), &
                              GLOB_AM0_POA(IIDEBUG,JJDEBUG,10,1,2,2)

         print*, 'POA to trop', SUM(Spc(:,:,1:State_Grid%MaxChemLev,id_POA1))+ &
                                SUM(Spc(:,:,1:State_Grid%MaxChemLev,id_POA2))
      ENDIF ! OPOA (hotp 8/24/09)

   ENDIF
   !------------------------------------------------------------------------

   ! Free pointer
   Spc => NULL()

 END SUBROUTINE SOA_CHEMISTRY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: soa_equil
!
! !DESCRIPTION: Function SOA\_EQUIL solves SOAeqn=0 to determine Mnew (= mass)
!  See Eqn (27) on page 70 of notes.  Originally written by Serena Chung at
!  Caltech, and modified for inclusion into GEOS-CHEM. (rjp, bmy, 7/8/04)
!\\
!\\
! !INTERFACE:
!
 FUNCTION SOA_EQUIL( MASS, MPOC, AEROSOL, GAS, KOM ) &
      RESULT( SOA_MASS )
!
! !INPUT/OUTPUT PARAMETERS:
!
   REAL(fp), INTENT(IN) :: MASS                ! Pre-existing aer mass [ug/m3]
   REAL(fp), INTENT(IN) :: MPOC                ! POA Mass [ug/m3]
   REAL(fp), INTENT(IN) :: AEROSOL(MPROD,MSV)  ! Aerosol concentration [ug/m3]
   REAL(fp), INTENT(IN) :: GAS(MPROD,MSV)      ! Gas-phase conc [ug/m3]
   REAL(fp), INTENT(IN) :: KOM(MPROD,MSV)      ! Equilibrium gas-aerosol
                                               !  partition coeff. [m3/ug]
!
! !RETURN VALUE:
!
   REAL(fp)             :: SOA_MASS
!
! !REMARKS:
!  This version does NOT assume that the gas and aerosol phases are in
!  equilibrium before chemistry; therefore, gas phase concentrations are
!  needed explicitly.  The gas and aerosol phases are assumed to be in
!  equilibrium after chemistry.
!                                                                             .
!  Note: Unlike FUNCTION SOA, this function assumes no reactions.  It only
!  considers the partitioning of existing products of VOC oxidation.
!                                                                             .
!  HC_JHC + OXID_IOXID - >
!    alpha(1,IOXID,JHC) [SOAprod_gas(1,IOXID,JHC)+SOAprod(1,IOXID,JHC)]+
!    alpha(2,IOXID,JHC) [SOAprod_gas(2,IOXID,JHC)+SOAprod(2,IOXID,JHC)]
!                                                                             .
!  SOAprod_gas(IPR,IOXID,JHC) <--> SOAprod(IPR,IOXID,JHC)
!                                           (aerosol phase)
!                                                                             .
!  w/ equilibrium partitioning:
!                                                                             .
!                                   SOAprod(IPR,IOXID,JHC)
!    SOAprod_gas(IPR,IOXID,JHC) = ------------------------
!                                     Kom(IPR,IOXID,JHC)
!
!  NOTES:
!  08 Jul 2004 - R. Park - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER  :: JHC,   IPR!,     NOX (hotp 5/22/10)
   INTEGER  :: JSV ! hotp 5/13/10
   REAL(fp) :: VALUE

   !=================================================================
   ! SOA_EQUIL begins here!
   !=================================================================

   ! Equation (39) on page 139 of notes:
   VALUE = 0.e+0_fp

   ! Use SV not HC (hotp 5/13/10)
   ! update dims (remove NOX) (hotp 5/22/10)
   DO JSV = 1, MAXSIMSV
   DO IPR = 1, NPROD(JSV)
      VALUE = VALUE + KOM(IPR,JSV)                        / &
                      ( 1.e+0_fp + KOM(IPR,JSV) * MASS  ) * &
                      ( GAS(IPR,JSV) + AEROSOL(IPR,JSV) )
   ENDDO
   ENDDO

   ! Compute SOA mass
   SOA_MASS = VALUE + ( 1.e+5_fp * MPOC ) / ( 1.e+5_fp * MASS ) - 1.0e+0_fp

 END FUNCTION SOA_EQUIL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: zeroin
!
! !DESCRIPTION: Function ZEROIN computes a zero of the function f(x) in the
!  interval ax,bx.
!\\
!\\
! !INTERFACE:
!
 FUNCTION ZEROIN(AX,BX,TOL,MPOC,AEROSOL,GAS,KOM) &
      RESULT( MNEW )
!
! !INPUT PARAMETERS:
!
   REAL(fp), INTENT(IN) :: ax
   REAL(fp), INTENT(IN) :: bx
   REAL(fp), INTENT(IN) :: tol
   REAL(fp), INTENT(IN) :: Mpoc
   REAL(fp), INTENT(IN) :: Aerosol(MPROD,MSV)
   REAL(fp), INTENT(IN) :: Gas(MPROD,MSV)
   REAL(fp), INTENT(IN) :: Kom(MPROD,MSV)
!
! !RETURN VALUE:
!
   REAL(fp)             :: MNEW
!
! !REMARKS:
! NOTE: This function may be problematic -- it uses GOTO's, which are not
! good for parallelization. (bmy, 7/8/04)
!                                                                             .
! shc I got this code from http://www.netlib.org
!                                                                             .
!      a zero of the function  f(x)  is computed in the interval ax,bx .
!                                                                             .
!  input..
!                                                                             .
!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!         the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
!         final result ( .ge. 0.0e+0_fp)
!                                                                             .
!  output..
!                                                                             .
!  zeroin abcissa approximating a zero of  f  in the interval ax,bx
!                                                                             .
!      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  without  a  check.  zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
!  is the relative machine precision.
!      this function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice - hall, inc. (1973).
!
! !REVISION HISTORY:
!  08 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   REAL(fp)             :: a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
!
!  compute eps, the relative machine precision
!
   eps = 1.0e+0_fp
10 eps = eps/2.0e+0_fp
   tol1 = 1.0e+0_fp + eps
   if (tol1 .gt. 1.0e+0_fp) go to 10
!
! initialization
!
   a  = ax
   b  = bx
   fa = SOA_equil( A, MPOC, Aerosol, GAS, Kom )
   fb = SOA_equil( B, MPOC, Aerosol, GAS, Kom )
!
! begin step
!
20 c = a
   fc = fa
   d = b - a
   e = d

30 if (ABS(fc) .ge. ABS(fb)) go to 40
   a = b
   b = c
   c = a
   fa = fb
   fb = fc
   fc = fa
!
! convergence test
!
40 tol1 = 2.0e+0_fp*eps*ABS(b) + 0.5e+0_fp*tol
   xm = 0.5e+0_fp*(c - b)
   if (ABS(xm) .le. tol1) go to 90
   if (fb .eq. 0.0e+0_fp) go to 90
!
! is bisection necessary
!
   if (ABS(e) .lt. tol1) go to 70
   if (ABS(fa) .le. ABS(fb)) go to 70
!
! is quadratic interpolation possible
!
   if (a .ne. c) go to 50
!
! linear interpolation
!
   s = fb/fa
   p = 2.0e+0_fp*xm*s
   q = 1.0e+0_fp - s
   go to 60
!
! inverse quadratic interpolation
!
50 q = fa/fc
   r = fb/fc
   s = fb/fa
   p = s*(2.0e+0_fp*xm*q*(q - r) - (b - a)*(r - 1.0e+0_fp))
   q = (q - 1.0e+0_fp)*(r - 1.0e+0_fp)*(s - 1.0e+0_fp)
!
! adjust signs
!
60 if (p .gt. 0.0e+0_fp) q = -q
   p = ABS(p)
!
! is interpolation acceptable
!
   if ((2.0e+0_fp*p) .ge. (3.0e+0_fp*xm*q - ABS(tol1*q))) go to 70
   if (p .ge. ABS(0.5e+0_fp*e*q)) go to 70

   e = d
   d = p/q
   go to 80
!
! bisection
!
70 d = xm
   e = d
!
! complete step
!
   80 a = b
   fa = fb
   if (ABS(d) .gt. tol1) b = b + d
   if (ABS(d) .le. tol1) b = b + SIGN(tol1, xm)

   fb = SOA_equil( B, MPOC, Aerosol, GAS, Kom )
   if ((fb*(fc/ABS(fc))) .gt. 0.0e+0_fp) go to 20
   go to 30
!
! done
!
90 MNEW = b

 END FUNCTION ZEROIN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rtbis
!
! !DESCRIPTION: Function RTBIS finds the root of the function SOA\_EQUIL via
!  the bisection method.  Original algorithm from "Numerical Recipes" by Press
!  et al, Cambridge UP, 1986.  Modified for inclusion into GEOS-CHEM.
!  (bmy, 7/8/04)
!\\
!\\
! !INTERFACE:
!
 FUNCTION RTBIS( X1, X2, XACC, MPOC, AEROSOL, GAS, KOM ) &
      RESULT( ROOT )
!
! !USES:
!
   USE ERROR_MOD, ONLY : ERROR_STOP
!
! !INPUT PARAMETERS:
!
   REAL(fp), INTENT(IN) :: X1                 ! Endpoint #1
   REAL(fp), INTENT(IN) :: X2                 ! Endpoint #2
   REAL(fp), INTENT(IN) :: XACC               ! Desired accuracy of solution
   REAL(fp), INTENT(IN) :: MPOC               ! POA mass [ug/m3]
   REAL(fp), INTENT(IN) :: AEROSOL(MPROD,MSV) ! Aerosol concentration [ug/m3]
   REAL(fp), INTENT(IN) :: GAS(MPROD,MSV)     ! Gas-phase concentration [ug/m3]
   REAL(fp), INTENT(IN) :: KOM(MPROD,MSV)     ! Equilibrium gas-aerosol
                                              !  partition coeff. [m3/ug]
!
! !RETURN VALUE:
!
   REAL(fp)             :: ROOT
!
! !REVISION HISTORY:
!  08 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
   INTEGER, PARAMETER :: JMAX = 100
!
! !LOCAL VARIABLES:
!
   INTEGER            :: J
   REAL(fp)           :: DX, F, FMID, XMID

   !=================================================================
   ! RTBIS begins here!
   !=================================================================

   ! Compute value of function SOA_EQUIL at endpoints
   FMID = SOA_EQUIL( X2, MPOC, AEROSOL, GAS, KOM )
   F    = SOA_EQUIL( X1, MPOC, AEROSOL, GAS, KOM )

   ! Test if we are bracketing a root
   IF ( F * FMID >= 0e+0_fp ) THEN
      CALL ERROR_STOP( 'Root must be bracketed!', 'RTBIS ("carbon_mod.F90")' )
   ENDIF

   ! Set initial root and interval
   IF ( F < 0e+0_fp ) THEN
      ROOT = X1
      DX   = X2 - X1
   ELSE
      ROOT = X2
      DX   = X1 - X2
   ENDIF

   ! Loop until max iteration count
   DO J = 1, JMAX

      ! Halve the existing interval
      DX   = DX * 0.5e+0_fp

      ! Compute midpoint of new interval
      XMID = ROOT + DX

      ! Compute value of function SOA_EQUIL at new midpoint
      FMID = SOA_EQUIL( XMID, MPOC, AEROSOL, GAS, KOM )

      ! We have found the root!
      IF ( FMID <= 0e+0_fp ) ROOT = XMID

      ! We have reached the tolerance, so return
      IF ( ABS( DX ) < XACC .OR. FMID == 0.e+0_fp ) RETURN
   ENDDO

   ! Stop with error condition
   CALL ERROR_STOP( 'Too many bisections!', 'RTBIS ("carbon_mod.F90")' )

 END FUNCTION RTBIS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: soa_para
!
! !DESCRIPTION: Subroutine SOA\_PARA gves mass-based stoichiometric
!  coefficients for semi-volatile products from the oxidation of hydrocarbons.
!  It calculates secondary organic aerosol yield parameters.  Temperature
!  effects are included.  Original code from the CALTECH group and modified for
!  inclusion to GEOS-CHEM. (rjp, bmy, 7/8/04, 6/30/08)
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE SOA_PARA( TEMP, KO3, KOH, KNO3, KOM, II, JJ, LL, &
                      KRO2NO, KRO2HO2, State_Met )
!
! !USES:
!
   USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
   INTEGER,        INTENT(IN)  :: II              ! Longitude index
   INTEGER,        INTENT(IN)  :: JJ              ! Latitude index
   INTEGER,        INTENT(IN)  :: LL              ! Altitude index
   REAL(fp),       INTENT(IN)  :: TEMP            ! Temperature [k]
   TYPE(MetState), INTENT(IN)  :: State_Met       ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
   REAL(fp),       INTENT(OUT) :: KO3(MHC)        ! Rxn rate for HC oxidation
                                                  !  by O3 [cm3/molec/s]
   REAL(fp),       INTENT(OUT) :: KOH(MHC)        ! Rxn rate for HC oxidation
                                                  !  by OH [cm3/molec/s]
   REAL(fp),       INTENT(OUT) :: KNO3(MHC)       ! Rxn rate for HC oxidation
                                                  !  by NO3 [cm3/molec/s]
   REAL(fp),       INTENT(OUT) :: KOM(MPROD,MSV)  ! Equilibrium gas-aerosol
                                                  !  partition coeff [m3/ug]

   ! RO2+NO,HO2 rate constants (hotp 5/7/10)
   REAL(fp),       INTENT(OUT) :: KRO2NO          ! RO2+NO  rate constant
   REAL(fp),       INTENT(OUT) :: KRO2HO2         ! RO2+HO2 rate constant
!
! !REMARKS:
!  References:
!  ============================================================================
!  PHOTO-OXIDATION RATE CONSTANTS OF ORGANICS come from:
!  (1 ) Atkinson, el al., Int. J. Chem.Kinet., 27: 941-955 (1995)
!  (2 ) Shu and Atkinson, JGR 100: 7275-7281 (1995)
!  (3 ) Atkinson, J. Phys. Chem. Ref. Data 26: 215-290 (1997)
!  (4 ) Some are reproduced in Table 1 of Griffin, et al., JGR 104: 3555-3567
!  (5 ) Chung and Seinfeld (2002)
!                                                                             .
!  ACTIVATION ENERGIES come from:
!  (6 ) Atkinson, R. (1994) Gas-Phase Tropospheric Chemistry of Organic
!        Compounds.  J. Phys. Chem. Ref. Data, Monograph No.2, 1-216.
!  (7 ) They are also reproduced in Tables B.9 and B.10 of Seinfeld and
!        Pandis (1988).
!                                                                             .
!  PARAMETERS FOR ISOPRENE:
!  (8 ) Kroll et al., GRL, 109, L18808 (2005)
!  (9 ) Kroll et al., Environ Sci Tech, in press (2006)
!  (10) Henze and Seinfeld, GRL, submitted (2006)
!
! !REVISION HISTORY:
!  08 Jul 2004 - R. Park - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
   ! Activation Energy/R [K] for O3, OH, NO3 (see Refs #6-7)
   REAL(fp), PARAMETER :: ACT_O3     =  732.0e+0_fp
   REAL(fp), PARAMETER :: ACT_OH     = -400.0e+0_fp
   REAL(fp), PARAMETER :: ACT_NO3    = -490.0e+0_fp

   ! Heat of vaporization (from CRC Handbook of Chemistry & Physics)
   REAL(fp), PARAMETER :: HEAT_VAPOR = 5.e+3_fp

   ! Reciprocal reference temperatures at 298K and 310K
   !REAL(fp), PARAMETER :: REF295     = 1e+0_fp / 295e+0_fp !hotp 5/21/10
   REAL(fp), PARAMETER :: REF298     = 1e+0_fp / 298e+0_fp
   !REAL(fp), PARAMETER :: REF310     = 1e+0_fp / 310e+0_fp !hotp 5/21/10
   ! semivolpoa2: reference T for POA (hotp 2/27/09)
   REAL(fp), PARAMETER :: REF300     = 1e+0_fp / 300e+0_fp
!
! !LOCAL VARIABLES:
!
   INTEGER             :: IPR,  JHC
   INTEGER             :: JSV ! (hotp 5/13/10)
   REAL(fp)            :: TMP1, TMP2, TMP3, OVER

   !=================================================================
   ! SOA_PARA begins here!
   !=================================================================

   ! move to SOA_PARA_INIT (dkh, 11/12/06)
   !! Photo-oxidation rates of O3 [cm3/molec/s] (See Refs #1-4)
   !KO3(1) = 56.15d-18
   !KO3(2) = 200.d-18
   !KO3(3) = 7707.d-18
   !KO3(4) = 422.5d-18
   !KO3(5) = ( 11600.D0 + 11700.e+0_fp ) / 2.e+0_fp * 1.D-18
   !
   !! Photo-oxidation rates of OH [cm3/molec/s] (See Refs #1-4)
   !KOH(1) = 84.4d-12
   !KOH(2) = 171.d-12
   !KOH(3) = 255.d-12
   !KOH(4) = 199.d-12
   !KOH(5) = ( 197.e+0_fp + 293.e+0_fp ) / 2.e+0_fp * 1.d-12
   !
   !! Photo-oxidation rate of NO3 [cm3/molec/s] (See Refs #1-4)
   !KNO3(1) = 6.95d-12
   !KNO3(2) = 12.2d-12
   !KNO3(3) = 88.7d-12
   !KNO3(4) = 14.7d-12
   !KNO3(5) = ( 19.e+0_fp + 35.e+0_fp ) / 2.e+0_fp * 1.d-12

   !=================================================================
   ! Temperature Adjustments of KO3, KOH, KNO3
   !=================================================================

   ! Initialize to zero (hotp 5/21/10)
   KO3  = 0e+0_fp
   KOH  = 0e+0_fp
   KNO3 = 0e+0_fp

   ! Reciprocal temperature [1/K]
   OVER = 1.0e+0_fp / TEMP

   ! Compute the exponentials once outside the DO loop
   TMP1 = EXP( ACT_O3  * ( REF298 - OVER ) )
   TMP2 = EXP( ACT_OH  * ( REF298 - OVER ) )
   TMP3 = EXP( ACT_NO3 * ( REF298 - OVER ) )

   ! Multiply photo-oxidation rates by exponential of temperature
   !(dkh, 10/08/05)
   !DO JHC = 1, 5
   DO JHC = 1, 4 ! now 4 (hotp 5/21/10)
      !KO3(JHC)  = KO3(JHC)  * TMP1
      !KOH(JHC)  = KOH(JHC)  * TMP2
      !KNO3(JHC) = KNO3(JHC) * TMP3
      KO3(JHC)  = KO3_REF(JHC)  * TMP1
      KOH(JHC)  = KOH_REF(JHC)  * TMP2
      KNO3(JHC) = KNO3_REF(JHC) * TMP3
   ENDDO

   !=================================================================
   ! Calculate KRO2NO, KRO2HO2 at TEMPERATURE (hotp 5/7/10)
   !=================================================================
   KRO2NO  = AARO2NO  * EXP( BBRO2NO  * OVER )
   KRO2HO2 = AARO2HO2 * EXP( BBRO2HO2 * OVER )

   !!=================================================================
   !! SOA YIELD PARAMETERS
   !!
   !! Aerosol yield parameters for photooxidation of biogenic organics
   !! The data (except for C7-C10 n-carbonyls, aromatics, and higher
   !! ketones are from:
   !!
   !! (7) Tables 1 and 2 of Griffin, et al., Geophys. Res. Lett.
   !!      26: (17)2721-2724 (1999)
   !!
   !! These parameters neglect contributions of the photooxidation
   !! by NO3.
   !!
   !! For the aromatics, the data are from
   !! (8) Odum, et al., Science 276: 96-99 (1997).
   !!
   !! Isoprene (dkh, bmy, 5/22/06)
   !! Unlike the other species, we consider oxidation by purely OH.
   !! CHEM_NVOC has been adjusted accordingly. There's probably
   !! significant SOA formed from NO3 oxidation, but we don't know
   !! enough to include that yet.  Data for the high NOX and low NOX
   !! parameters are given in Kroll 05 and Kroll 06, respectively.
   !! The paramters for low NOX are given in Table 1 of Henze 06.
   !!=================================================================
   !
   !! Average of ALPHA-PINENE, BETA-PINENE, SABINENE, D3-CARENE
   !RALPHA(1,1) = 0.067e+0_fp
   !RALPHA(2,1) = 0.35425e+0_fp
   !
   !! LIMONENE
   !RALPHA(1,2) = 0.239e+0_fp
   !RALPHA(2,2) = 0.363e+0_fp
   !
   !! Average of TERPINENES and TERPINOLENE
   !RALPHA(1,3) = 0.0685e+0_fp
   !RALPHA(2,3) = 0.2005e+0_fp
   !
   !! Average of MYRCENE, LINALOOL, TERPINENE-4-OL, OCIMENE
   !RALPHA(1,4) = 0.06675e+0_fp
   !RALPHA(2,4) = 0.135e+0_fp
   !
   !! Average of BETA-CARYOPHYLLENE and and ALPHA-HUMULENE
   !RALPHA(1,5) = 1.0e+0_fp
   !RALPHA(2,5) = 0.0e+0_fp
   !
   !! Using BETA-PINENE for all species for NO3 oxidation
   !! Data from Table 4 of Griffin, et al., JGR 104 (D3): 3555-3567 (1999)
   !RALPHA(3,:) = 1.e+0_fp
   !
   !! Here we define some alphas for isoprene (dkh, bmy, 5/22/06)
   !
   !! high NOX  [Kroll et al, 2005]
   !!RALPHA(1,6) = 0.264e+0_fp
   !!RALPHA(2,6) = 0.0173e+0_fp
   !!RALPHA(3,6) = 0e+0_fp
   !
   !! low NOX   [Kroll et al, 2006; Henze and Seinfeld, 2006]
   !RALPHA(1,6) = 0.232e+0_fp
   !RALPHA(2,6) = 0.0288e+0_fp
   !RALPHA(3,6) = 0e+0_fp
   !
   !!=================================================================
   !! Equilibrium gas-particle partition coefficients of
   !! semi-volatile compounds [ug-1 m**3]
   !!=================================================================
   !
   !! Average of ALPHA-PINENE, BETA-PINENE, SABINENE, D3-CARENE
   !KOM(1,1) = 0.1835e+0_fp
   !KOM(2,1) = 0.004275e+0_fp
   !
   !! LIMONENE
   !KOM(1,2) = 0.055e+0_fp
   !KOM(2,2) = 0.0053e+0_fp
   !
   !! Average of TERPINENES and TERPINOLENE
   !KOM(1,3) = 0.133e+0_fp
   !KOM(2,3) = 0.0035e+0_fp
   !
   !! Average of MYRCENE, LINALOOL, TERPINENE-4-OL, OCIMENE
   !KOM(1,4) = 0.22375e+0_fp
   !KOM(2,4) = 0.0082e+0_fp
   !
   !! Average of BETA-CARYOPHYLLENE and and ALPHA-HUMULENE
   !KOM(1,5) = ( 0.04160e+0_fp + 0.0501e+0_fp ) / 2.e+0_fp
   !KOM(2,5) = 0.0e+0_fp
   !
   !! NOT APPLICABLE -- using BETA-PINENE for all species
   !! Data from Table 4 of Griffin, et al., JGR 104 (D3): 3555-3567 (1999)
   !KOM(3,:) = 0.0163e+0_fp
   !
   !! Again, for isoprene we only consider two products,
   !! both from OH oxidation. (dkh, bmy, 5/22/06)
   !
   !! High NOX
   !!KOM(1,6) = 0.00115e+0_fp
   !!KOM(2,6) = 1.52e+0_fp
   !!KOM(3,6) = 0e+0_fp
   !
   !! Low NOX
   !KOM(1,6) = 0.00862e+0_fp
   !KOM(2,6) = 1.62e+0_fp
   !KOM(3,6) = 0e+0_fp

   !=================================================================
   ! Temperature Adjustments of KOM
   !=================================================================

   !--------------------------------------------------------
   ! Lumped semivolatiles 1-3 (hotp 5/21/10)
   !--------------------------------------------------------
   ! First 3 semivolatile systems are at Tref = 298K
   ! SV 1: MTPA,LIMO,MTPO,SESQ
   ! SV 2: ISOP
   ! SV 3: BENZ,TOLU,XYLE,(NAP)

   ! Reciprocal temperature [1/K]
   OVER = 1.0e+0_fp / TEMP

   !! Divide TEMP by 310K outside the DO loop
   !TMP1 = ( TEMP / 310.e+0_fp )
   ! Divide TEMP by 298K outside the DO loop
   TMP1 = ( TEMP / 298.e+0_fp )

   ! Compute the heat-of-vaporization exponential term outside the DO loop
   !TMP2 = EXP( HEAT_VAPOR * ( OVER - REF310 ) )
   TMP2 = EXP( HEAT_VAPOR * ( OVER - REF298 ) )

   ! Multiply KOM by the temperature and heat-of-vaporization terms
   ! now use JSV (hotp 5/21/10)
   ! update dims (hotp 5/22/10)
   DO JSV = 1, 3
   DO IPR = 1, NPROD(JSV)
      KOM(IPR,JSV) = KOM_REF(IPR,JSV) * TMP1 * TMP2
   ENDDO
   ENDDO

   !--------------------------------------------------------
   ! POA (primary semivolatiles) (hotp 5/13/10)
   !--------------------------------------------------------
   ! semivolpoa2: reference for POA is 300 K (hotp 2/27/09)
   ! Divide TEMP by 300K outside the DO loop
   TMP1 = ( TEMP / 300.e+0_fp )

   ! Compute the heat-of-vaporization exponential term outside the DO loop
   TMP2 = EXP( HEAT_VAPOR * ( OVER - REF300 ) )

   ! Multiply KOM by the temperature and heat-of-vaporization terms
   ! Adjust POA from reference of 300K
   JHC = PARENTPOA
   JSV = IDSV(JHC)
   DO IPR = 1, NPROD(JSV)
      KOM(IPR,JSV) = KOM_REF(IPR,JSV) * TMP1 * TMP2
   ENDDO

   !--------------------------------------------------------
   ! OPOA (oxidized semivolatiles) (hotp 5/13/10)
   !--------------------------------------------------------
   ! Divide TEMP by 300K outside the DO loop
   TMP1 = ( TEMP / 300.e+0_fp )

   ! Compute the heat-of-vaporization exponential term outside the DO loop
   TMP2 = EXP( HEAT_VAPOR * ( OVER - REF300 ) )

   ! Adjust OPOA KOM
   JHC = PARENTOPOA
   JSV = IDSV(JHC)
   DO IPR = 1, NPROD(JSV)
      KOM(IPR,JSV) = KOM_REF(IPR,JSV) * TMP1 * TMP2
   ENDDO

 END SUBROUTINE SOA_PARA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: soa_para_init
!
! !DESCRIPTION: Subroutine SOA\_PARA\_INIT initializes the ALPHAS and KOMS, the
!  latter at their reference temperature. It is faster to define these
!  seperately as it only needs to be done once. (dkh, 11/12/06)
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE SOA_PARA_INIT( Input_Opt )
!
! !USES:
!
   USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !REMARKS:
!  NOTE: REFT for KOM_REF depends on hydrocarbon.
!
! !REVISION HISTORY:
!  12 Nov 2006 - D. Henze - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! for debug purposes (hotp 7/22/09)
   INTEGER :: ai,bj,cl

   INTEGER :: NOX ! (hotp 5/21/10)

   !=================================================================
   ! SOA_PARA_INIT begins here!
   !=================================================================

   ! update reaction rates
   ! Still based on same underlying data (as summarized by
   ! Griffin 1999 and Chung 2002) but lumped by contribution of HC to
   ! global emissions for that category (hotp 5/21/10)
   ! K(MTPA) = K_REF(1) =
   !      0.53*K(A-PINE) + 0.25*K(B-PINE) + 0.12*K(SABI) + 0.10*K(CAR)
   ! K(LIMO) = K_REF(2)
   ! K(MTPO) = K_REF(3) =
   !      0.11*K(TERPINENE) + 0.11*K(TERPINOLENE) + 0.11*K(MYRCENE) +
   !      0.11*K(LINALOOL) + 0.11*K(terpinene-4-ol) + 0.45*K(OCIMENE)
   ! K(SESQ) = K_REF(4) = 0.5*K(B-CARYOPHYLLENE) + 0.5*K(A-HUMULENE)

   ! Photo-oxidation rates of O3 [cm3/molec/s] (See Refs #1-4)
   KO3_REF(1) =    63.668e-18_fp
   KO3_REF(2) =   200.000e-18_fp
   KO3_REF(3) =  1744.500e-18_fp
   KO3_REF(4) = 11650.000e-18_fp

   ! Photo-oxidation rates of OH [cm3/molec/s] (See Refs #1-4)
   KOH_REF(1) =    71.026e-12_fp
   KOH_REF(2) =   171.000e-12_fp
   KOH_REF(3) =   227.690e-12_fp
   KOH_REF(4) =   245.000e-12_fp

   ! Photo-oxidation rate of NO3 [cm3/molec/s] (See Refs #1-4)
   KNO3_REF(1) =    6.021e-12_fp
   KNO3_REF(2) =   12.200e-12_fp
   KNO3_REF(3) =   33.913e-12_fp
   KNO3_REF(4) =   27.000e-12_fp

   ! Rate constants for branching ratio (hotp 5/7/10)
   ! k=A*exp(B/T)
   ! Reference: Henze et al., 2008 ACP
   ! RO2+NO
   AARO2NO = 2.6e-12_fp
   BBRO2NO = 350.e+0_fp
   ! RO2+HO2
   AARO2HO2 = 1.4e-12_fp
   BBRO2HO2 = 700.e+0_fp

   !=================================================================
   ! SOA YIELD PARAMETERS
   !
   ! Aerosol yield parameters for photooxidation of biogenic organics
   ! The data (except for C7-C10 n-carbonyls, aromatics, and higher
   ! ketones are from:
   !
   ! (7) Tables 1 and 2 of Griffin, et al., Geophys. Res. Lett.
   !      26: (17)2721-2724 (1999)
   !
   ! These parameters neglect contributions of the photooxidation
   ! by NO3.
   !
   ! For the aromatics, the data are from
   ! (8) Odum, et al., Science 276: 96-99 (1997).
   !
   ! Isoprene (dkh, bmy, 5/22/06)
   ! Unlike the other species, we consider oxidation by purely OH.
   ! CHEM_NVOC has been adjusted accordingly. There's probably
   ! significant SOA formed from NO3 oxidation, but we don't know
   ! enough to include that yet.  Data for the high NOX and low NOX
   ! parameters are given in Kroll 05 and Kroll 06, respectively.
   ! The paramters for low NOX are given in Table 1 of Henze 06.
   !=================================================================

   ! SOAupdate: new yield parameterizations
   ! Initialize all ALPHAs to zero (hotp 5/12/10)
   ! ALPHAs are indexed by PARENT HYDROCARBON
   ! all monoterpenes use b-pinene + NO3 for NO3 yields
   ALPHA = 0e+0_fp

   !----------------------------
   ! MTPA
   !----------------------------
   ! MTPA based on Shilling 2008 a-pinene ozonolysis
   ! updated 6/12/10 (hotp)
   ! Product 4 has C*=0.1
   NOX = NHIGHNOX
   ALPHA(NOX,1,PARENTMTPA) = 0.0095e+0_fp
   ALPHA(NOX,2,PARENTMTPA) = 0.0900e+0_fp
   ALPHA(NOX,3,PARENTMTPA) = 0.0150e+0_fp
   ALPHA(NOX,4,PARENTMTPA) = 0.0400e+0_fp
   NOX = NLOWNOX
   ALPHA(NOX,1,PARENTMTPA) = 0.019e+0_fp
   ALPHA(NOX,2,PARENTMTPA) = 0.180e+0_fp
   ALPHA(NOX,3,PARENTMTPA) = 0.030e+0_fp
   ALPHA(NOX,4,PARENTMTPA) = 0.080e+0_fp
   NOX = NNO3RXN
   ALPHA(NOX,1,PARENTMTPA) = 0.0000e+0_fp
   ALPHA(NOX,2,PARENTMTPA) = 0.3207e+0_fp
   ALPHA(NOX,3,PARENTMTPA) = 1.0830e+0_fp

   !----------------------------
   ! LIMO
   !----------------------------
   ! Use higher LIMO yields of Zhang 2006
   ! Assumed density of 1.3 g/cm3 (hotp 6/12/10)
   NOX = NHIGHNOX
   ALPHA(NOX,1,PARENTLIMO) = 0.4743e+0_fp
   ALPHA(NOX,2,PARENTLIMO) = 0.1174e+0_fp
   ALPHA(NOX,3,PARENTLIMO) = 1.4190e+0_fp
   NOX = NLOWNOX
   ALPHA(NOX,1,PARENTLIMO) = 0.3661e+0_fp
   ALPHA(NOX,2,PARENTLIMO) = 0.3214e+0_fp
   ALPHA(NOX,3,PARENTLIMO) = 0.8168e+0_fp
   NOX = NNO3RXN
   ALPHA(NOX,1,PARENTLIMO) = 0.0000e+0_fp
   ALPHA(NOX,2,PARENTLIMO) = 0.3207e+0_fp
   ALPHA(NOX,3,PARENTLIMO) = 1.0830e+0_fp

   !----------------------------
   ! MTPO
   !----------------------------
   ! MTPO based on Shilling 2008 a-pinene ozonolysis
   ! updated 6/12/10 (hotp)
   ! Product 4 has C*=0.1
   NOX = NHIGHNOX
   ALPHA(NOX,1,PARENTMTPO) = 0.0095e+0_fp
   ALPHA(NOX,2,PARENTMTPO) = 0.0900e+0_fp
   ALPHA(NOX,3,PARENTMTPO) = 0.0150e+0_fp
   ALPHA(NOX,4,PARENTMTPO) = 0.040e+0_fp
   NOX = NLOWNOX
   ALPHA(NOX,1,PARENTMTPO) = 0.019e+0_fp
   ALPHA(NOX,2,PARENTMTPO) = 0.180e+0_fp
   ALPHA(NOX,3,PARENTMTPO) = 0.030e+0_fp
   ALPHA(NOX,4,PARENTMTPO) = 0.080e+0_fp
   NOX = NNO3RXN
   ALPHA(NOX,1,PARENTMTPO) = 0.0000e+0_fp
   ALPHA(NOX,2,PARENTMTPO) = 0.3207e+0_fp
   ALPHA(NOX,3,PARENTMTPO) = 1.0830e+0_fp

   !----------------------------
   ! SESQ
   !----------------------------
   ! update high and low NOx (hotp 6/4/2010)
   ! Griffin1999 VOC/NO>3ppbC/ppb is low NOx
   ! high NOx is double the Y for a given Mo
   NOX = NHIGHNOX
   ALPHA(NOX,1,PARENTSESQ) = 0.0005e+0_fp
   ALPHA(NOX,2,PARENTSESQ) = 1.1463e+0_fp
   ALPHA(NOX,3,PARENTSESQ) = 2.9807e+0_fp
   NOX = NLOWNOX
   ALPHA(NOX,1,PARENTSESQ) = 0.0000e+0_fp
   ALPHA(NOX,2,PARENTSESQ) = 0.5738e+0_fp
   ALPHA(NOX,3,PARENTSESQ) = 1.4893e+0_fp
   NOX = NNO3RXN
   ALPHA(NOX,1,PARENTSESQ) = 0.0000e+0_fp
   ALPHA(NOX,2,PARENTSESQ) = 0.3207e+0_fp
   ALPHA(NOX,3,PARENTSESQ) = 1.0830e+0_fp

   !----------------------------
   ! ISOP
   !----------------------------
   NOX = 1 ! low NOx/all OH rxn
   ALPHA(NOX,1,PARENTISOP) = 0.0306e+0_fp
   ALPHA(NOX,2,PARENTISOP) = 0.0000e+0_fp
   ALPHA(NOX,3,PARENTISOP) = 0.0945e+0_fp
   NOX = 2 ! NO3 rxn
   ALPHA(NOX,1,PARENTISOP) = 0.0000e+0_fp
   ALPHA(NOX,2,PARENTISOP) = 0.2171e+0_fp
   ALPHA(NOX,3,PARENTISOP) = 0.0919e+0_fp

   !----------------------------
   ! BENZ, TOLU, XYLE
   !----------------------------
   ! Replace Daven's numbers for BENZ, TOLU, XYLE with new numbers
   ! Numbers based on a 3 product fit to Ng 2007 data
   ! These numbers are for parent HC (no adjustment for ARO2)
   ! and correspond to C* of 1, 10, 100 in HIGH NOx case (hotp 5/12)

   ! HIGH NOX BENZ
   ALPHA(1,1,PARENTBENZ) = 0.0778e+0_fp
   ALPHA(1,2,PARENTBENZ) = 0.0000e+0_fp
   ALPHA(1,3,PARENTBENZ) = 0.7932e+0_fp
   ! LOW NOX BENZ (non-volatile)
   ALPHA(2,4,PARENTBENZ) = 0.37e+0_fp

   ! HIGH NOX TOLU
   ALPHA(1,1,PARENTTOLU) = 0.0315e+0_fp
   ALPHA(1,2,PARENTTOLU) = 0.0944e+0_fp
   ALPHA(1,3,PARENTTOLU) = 0.0800e+0_fp
   ! LOW NOX TOLU
   ALPHA(2,4,PARENTTOLU) = 0.30e+0_fp

   ! HIGH NOX XYLE
   ALPHA(1,1,PARENTXYLE) = 0.0250e+0_fp
   ALPHA(1,2,PARENTXYLE) = 0.0360e+0_fp
   ALPHA(1,3,PARENTXYLE) = 0.0899e+0_fp
   ! LOW NOX XYLE
   ALPHA(2,4,PARENTXYLE) = 0.36e+0_fp

   !----------------------------
   ! POA
   !----------------------------
   ! semivolpoa2: alphas for POA (hotp 2/27/09)
   ! based on Shrivastava et al. 2006 ES&T
   ! Only 2 products (wood smoke)
   ALPHA(1,1,PARENTPOA) = 0.49e+0_fp
   ALPHA(1,2,PARENTPOA) = 0.51e+0_fp
   ! No high NOx parameters
   ! semivolpoa3: add diesel/anthropogenic POA (hotp 3/13/09)
   !ALPHA(2:MNOX,1:MPROD,10) = 0e+0_fp

   !----------------------------
   ! OPOA
   !----------------------------
   ! semivolpoa4opoa: alphas for OPOA (hotp 3/18/09)
   ! remove semivolpoa3 changes (hotp 3/27/09)
   ! biomass burning
   ! (note that this is the carbon yield)
   ALPHA(1,1,PARENTOPOA) = 1.e+0_fp
   ALPHA(1,2,PARENTOPOA) = 1.e+0_fp
   ! anthropogenic
   !ALPHA(2:MNOX,1:MPROD,11) = 0e+0_fp

   !----------------------------
   ! SOA from oxidation of IVOCs
   !----------------------------
   ! NAPSOA: SOA from oxidation of IVOCs (hotp 7/22/09)
   ! Values from Chan et al. 2009 ACP (refit)
   ! ALPHAs are set up for the aromatic (NAP) as the parent HC
   ! ALPHAs must be consistent with GET_DARO2 units!
   ! HIGH NOX
   ! Ox = NO
   ALPHA(1,1,PARENTNAP) = 0.0387e+0_fp
   ALPHA(1,2,PARENTNAP) = 0.2956e+0_fp
   ALPHA(1,3,PARENTNAP) = 0.2349e+0_fp
   ! LOW NOX
   ! Ox = HO2
   ALPHA(2,4,PARENTNAP) = 0.73e+0_fp

   !=================================================================
   ! Equilibrium gas-particle partition coefficients of
   ! semi-volatile compounds [ug-1 m**3]
   !=================================================================

   ! SOAupdate: KOM for semivolatile systems
   ! Initialize to zero (hotp 5/12/10)
   ! KOM_REF are indexed by SEMIVOLATILE SPECIES (hotp 5/13/10)
   KOM_REF = 0e+0_fp

   !---------------------------------------
   ! SEMIVOLATILE 1: MTPA, LIMO, MTPO, SESQ
   ! (hotp 5/21/10)
   !---------------------------------------
   KOM_REF(1,IDSV(PARENTMTPA)) = 1.0e+0_fp/1.0e+0_fp
   KOM_REF(2,IDSV(PARENTMTPA)) = 1.0e+0_fp/10.0e+0_fp
   KOM_REF(3,IDSV(PARENTMTPA)) = 1.0e+0_fp/100.0e+0_fp
   KOM_REF(4,IDSV(PARENTMTPA)) = 1.0e+0_fp/0.1e+0_fp ! C*=0.1 hotp 6/12/10

   !---------------------------------------
   ! SEMIVOLATILE 2: ISOP
   ! (hotp 5/21/10)
   !---------------------------------------
   KOM_REF(1,IDSV(PARENTISOP)) = 1.0e+0_fp/1.0e+0_fp
   KOM_REF(2,IDSV(PARENTISOP)) = 1.0e+0_fp/10.0e+0_fp
   KOM_REF(3,IDSV(PARENTISOP)) = 1.0e+0_fp/100.0e+0_fp

   !---------------------------------------
   ! SEMIVOLATILE 3: BENZ, TOLU, XYLE, NAP
   !---------------------------------------
   ! Update aromatics to new fits (hotp 5/12/10)
   ! BENZ, TOLU, XYLE, NAP/IVOC all lumped together
   KOM_REF(1,IDSV(PARENTBENZ)) = 1.0e+0_fp/1.0e+0_fp
   KOM_REF(2,IDSV(PARENTBENZ)) = 1.0e+0_fp/10.0e+0_fp
   KOM_REF(3,IDSV(PARENTBENZ)) = 1.0e+0_fp/100.0e+0_fp
   ! Low NOX (HO2) non-volatile
   !KOM_REF(4,IDSV(PARENTBENZ)) = 1.d6
   KOM_REF(4,IDSV(PARENTBENZ)) = 1.e+10_fp ! more non-vol (hotp 5/28/10)

   !---------------------------------------
   ! SEMIVOLATILE 4: POA/SVOCs
   !---------------------------------------
   ! semivolpoa2: KOM for POA (hotp 2/27/09)
   ! based on Shrivastava et al. 2006 ES&T
   ! Only 2 products (wood smoke)
   ! Tref is 27 C = 300 K
   KOM_REF(1,IDSV(PARENTPOA)) = 1e+0_fp/1646e+0_fp
   KOM_REF(2,IDSV(PARENTPOA)) = 1e+0_fp/20e+0_fp
   ! No high NOx parameters
   ! remove semivolpoa3 changes (hotp 3/27/09)
   ! semivolpoa3: add diesel/anthropogenic POA (hotp 3/13/09)
   !KOM_REF(2:MNOX,1:MPROD,10) = 0e+0_fp

   !---------------------------------------
   ! SEMIVOLATILE 5: OPOA/O-SVOCs
   !---------------------------------------
   ! semivolpoa4opoa: OPOA parameters (hotp 3/18/09)
   ! parameters are a factor of 100 more than POA param
   KOM_REF(1,IDSV(PARENTOPOA)) = KOM_REF(1,IDSV(PARENTPOA)) * 100e+0_fp
   KOM_REF(2,IDSV(PARENTOPOA)) = KOM_REF(2,IDSV(PARENTPOA)) * 100e+0_fp

   ! debug print checks (hotp 7/22/09)
   IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) THEN
      print*, 'Semivolatile POA settings:---------------'
      print*, ' ALPHA:   ', ALPHA(1,1,9), ALPHA(1,2,9)
      ! OCFPOA and OCFOPOA are now 2D arrays
      !print*, ' POA OA/OC ratio:    ', OCFPOA(I,J)
      !print*, ' OPOA OA/OC ratio:   ', OCFOPOA(I,J)
      print*, ' LSVPOA is set to:   ', Input_Opt%LSVPOA

      print*, 'CHECK MHC, NOX, PR', MHC, MNOX, MPROD
      print*, 'CHECK MSV', MSV
      print*, '      NOX, PROD, HC/SV'

      DO ai = 1, MHC
      DO bj = 1, MNOX
      DO cl = 1, MPROD
         print*,'Alpha', bj,cl,ai
         print*, ALPHA(bj,cl,ai)
      ENDDO
      ENDDO
      ENDDO

      ! Check KOM_REF (hotp 5/13/10)
      DO ai = 1, MSV
      DO cl = 1, MPROD
         print*,'KOM_REF', cl,ai
         print*, KOM_REF(cl,ai)
      ENDDO
      ENDDO
   ENDIF

 END SUBROUTINE SOA_PARA_INIT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_nvoc
!
! !DESCRIPTION: Subroutine CHEM\_NVOC computes the oxidation of Hydrocarbon by
!  O3, OH, and NO3.  This comes from the Caltech group (Hong Liao, Serena
!  Chung, et al) and was incorporated into GEOS-CHEM. (rjp, bmy, 7/6/04,6/1/06)
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE CHEM_NVOC( I, J, L, KO3, KOH, KNO3, GM0, KNO, KHO2, &
                       Input_Opt,  State_Chm, State_Diag,       &
                       State_Grid, State_Met, RC )
!
! !USES:
!
   USE ErrCode_Mod
   USE Input_Opt_Mod,      ONLY : OptInput
   USE State_Chm_Mod,      ONLY : ChmState
   USE State_Diag_Mod,     ONLY : DgnState
   USE State_Grid_Mod,     ONLY : GrdState
   USE State_Met_Mod,      ONLY : MetState
   USE TIME_MOD,           ONLY : GET_TS_CHEM, GET_MONTH
!
! !INPUT PARAMETERS:
!
   INTEGER,        INTENT(IN)    :: I             ! Longitude index
   INTEGER,        INTENT(IN)    :: J             ! Latitude index
   INTEGER,        INTENT(IN)    :: L             ! Altitude index
   REAL(fp),       INTENT(IN)    :: KO3(MHC)      ! Rxn rate for HC oxidation
                                                  !  by O3 [cm3/molec/s]
   REAL(fp),       INTENT(IN)    :: KOH(MHC)      ! Rxn rate for HC oxidation
                                                  !  by OH [cm3/molec/s]
   REAL(fp),       INTENT(IN)    :: KNO3(MHC)     ! Rxn rate for HC oxidation
                                                  !  by NO3 [cm3/molec/s]
   ! RO2+NO, RO2+HO2 rate constants (hotp 5/7/10)
   REAL(fp),       INTENT(IN)    :: KNO           ! RO2+NO  rate constant
   REAL(fp),       INTENT(IN)    :: KHO2          ! RO2+HO2 rate constant
   TYPE(OptInput), INTENT(IN)    :: Input_Opt     ! Input Options object
   TYPE(GrdState), INTENT(IN)    :: State_Grid    ! Grid State object
   TYPE(MetState), INTENT(IN)    :: State_Met     ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
   REAL(fp),       INTENT(INOUT) :: GM0(MPROD,MSV)! Gas mass for HCs and
                                                  !  oxidation products [kg]
   TYPE(ChmState), INTENT(INOUT) :: State_Chm     ! Chemistry State object
   TYPE(DgnState), INTENT(INOUT) :: State_Diag    ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC            ! Success or failure?
!
! !REMARKS:
!  SVOCs should immediately partition upon emission
!  SVOCs also react in the gas-phase
!  If SVOCs were emitted before reactions, we wouldn't know how
!  much to put in each phase
!  H.O.T. Pye decided to emit them after the existing SVOCs
!  react in the gas-phase. Thus the order of operations is:
!    SVOC + OH in gas-phase
!    SVOC emission (added to gas-phase GM0)
!    partitioning
!    dry dep
!    wet dep
!    etc
!
! !REVISION HISTORY:
!  06 Jul 2004 - R. Park - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER  :: JHC, IPR, NOX, JSV !(hotp 5/14/10)
   INTEGER  :: MAXLOOP ! (hotp 6/7/10)
   REAL(fp) :: CHANGE(MHC), NMVOC(MHC), DELHC(MNOX)
   REAL(fp) :: OHMC, TTNO3, TTO3, DTCHEM, RK
   REAL(fp) :: OVER, DO3, DOH, DNO3

   ! for RO2+NO, RO2+HO2 branching ratio (hotp 5/7/10)
   REAL(fp) :: NOTEMP, HO2TEMP, BETANO

   ! for debug (hotp 5/10/10)
   REAL(fp) :: TEMPHC

   ! Pointers
   REAL(fp), POINTER :: Spc(:,:,:,:)

   !=================================================================
   ! CHEM_NVOC begins here!
   !=================================================================

   ! Assume success
   RC      = GC_SUCCESS

   ! Point to chemical species array [kg]
   Spc     => State_Chm%Species

   ! Chemistry timestep [s]
   DTCHEM  = GET_TS_CHEM()

   ! Get offline OH, NO3, O3 concentrations [molec/cm3]
   OHMC    = GET_OH(  I, J, L, Input_Opt, State_Chm, State_Met )
   TTNO3   = GET_NO3( I, J, L, Input_Opt, State_Chm, State_Met )
   TTO3    = GET_O3(  I, J, L, Input_Opt, State_Chm, State_Grid, State_Met )

   ! Get RO2+NO, RO2+HO2 branching ratio (hotp 5/7/10)
   NOTEMP  = GET_NO(  I, J, L, Input_Opt, State_Chm, State_Met )
   HO2TEMP = GET_HO2( I, J, L, Input_Opt, State_Chm, State_Met )

   IF ( NOTEMP .GT. 0.0 ) THEN
      BETANO  = ( KNO * NOTEMP ) / ( KNO * NOTEMP + KHO2 * HO2TEMP )
   ELSEIF ( HO2TEMP .GT. 0.0 ) THEN
      BETANO = 0.e+0_fp
   ELSE
      ! default value if no State_Chm%Species value
      BETANO = 0.5e+0_fp
   ENDIF

   ! save for diagnostic purposes (hotp 5/24/10)
   BETANOSAVE(I,J,L) = BETANO

   ! Save to State_Diag for output as netcdf diagnostic (ewl, 8/28/18)
   IF ( State_Diag%Archive_BetaNO ) THEN
      State_Diag%BetaNO(I,J,L) = BETANO
   ENDIF

   ! update for new mtp lumping (hotp 5/22/10)
   NMVOC(1) = Spc(I,J,L,id_MTPA)
   NMVOC(2) = Spc(I,J,L,id_LIMO)
   NMVOC(3) = Spc(I,J,L,id_MTPO)
   NMVOC(4) = ORVC_SESQ(I,J,L)

   ! Initialize DELHC so that the values from the previous
   ! time step are not carried over.
   DELHC(:) = 0.e+0_fp

   !=================================================================
   ! Change in NVOC concentration due to photooxidation [kg]
   !=================================================================

   ! semivolpoa2: update for POA (hotp 2/27/09)
   ! add POA emissions to GMO here (not to Spc in EMITHIGH)

   ! Only loop over parent hydrocarbons defined for a given simulation
   ! Max should be 11 for semivolatile POA/IVOC (PARENTNAP =11)
   ! Max should be 8  for nonvolatile POA/ traditional simulation
   IF ( id_POA1 > 0 ) THEN
      MAXLOOP = PARENTNAP   !11
   ELSE
      MAXLOOP = PARENTXYLE  ! 8
   ENDIF

   DO JHC = 1, MAXLOOP

      ! Initialize again for safety (hotp 5/22/10)
      DELHC = 0e+0_fp

      ! Get JSV (hotp 5/14/10)
      JSV = IDSV(JHC)

      ! update for new mtp (hotp 5/22/10)
      IF ( JHC == PARENTMTPA .or. JHC == PARENTLIMO .or. &
           JHC == PARENTMTPO .or. JHC == PARENTSESQ      ) THEN

         !------------------------------------------
         ! Oxidize parent hydrocarbon by OH, O3, NO3
         ! (unmodified from original implemenation)
         !------------------------------------------
         RK          = KO3(JHC)*TTO3 + KOH(JHC)*OHMC + KNO3(JHC)*TTNO3
         CHANGE(JHC) = NMVOC(JHC) * ( 1.e+0_fp - &
                       EXP( -RK * DTCHEM ) ) !changed to EXP (myan, 12/14)

         ! In case that the biogenic hydrocarbon is the limiting reactant
         IF ( CHANGE(JHC) >= NMVOC(JHC) ) CHANGE(JHC) = NMVOC(JHC)

         ! NMVOC concentration after oxidation reactions
         NMVOC(JHC) = NMVOC(JHC) - CHANGE(JHC)

         IF( CHANGE(JHC) > 1.e-16_fp ) THEN
            OVER  = 1.e+0_fp / RK
            DO3   = CHANGE(JHC) * KO3(JHC)  * TTO3  * OVER ![kg]
            DOH   = CHANGE(JHC) * KOH(JHC)  * OHMC  * OVER ![kg]
            DNO3  = CHANGE(JHC) * KNO3(JHC) * TTNO3 * OVER ![kg]
         ELSE
            DO3   = 0.e+0_fp
            DOH   = 0.e+0_fp
            DNO3  = 0.e+0_fp
         ENDIF

         !------------------------------------------
         ! Determine DELTAHC that corresponds to the alphas
         !------------------------------------------
         ! For HC 1-4 (hotp 5/22/10)
         NOX = NHIGHNOX ! NOX=1, high NOx photooxidation
         DELHC(NOX) = ( DO3 + DOH ) * BETANO
         NOX = NLOWNOX  ! NOX=2, low NOx photooxidation
         DELHC(NOX) = ( DO3 + DOH ) * ( 1e+0_fp - BETANO )
         NOX = NNO3RXN  ! NOX=3, NO3 oxidation
         DELHC(NOX) = ( DNO3 )

         ! debug check (updated hotp 5/26/10)
         !IF ( CHANGE(JHC) .GT. 1d-16 ) THEN
         !TEMPHC = ABS(SUM(DELHC(:))-CHANGE(JHC))
         !TEMPHC = ABS(TEMPHC/CHANGE(JHC))
         !IF ( (TEMPHC) .GE. 1d-14 ) THEN
         !   print*,'DELHC Problem in CHEM_NVOC',I,J,L,JHC
         !   print*,DELHC,CHANGE(JHC),TEMPHC
         !ENDIF
         !ENDIF

         ! Save diagnostic info for bug check (hotp 5/22/10)
         DELTASOGSAVE(I,J,L,:,JHC) = DELHC(:)

         !------------------------------------------
         ! Compute amount of semivolatile formed
         ! and add to initial SOG
         !------------------------------------------
         ! update dims and switch order (hotp 5/22/10)
         DO NOX = 1, NNOX(JSV)
         DO IPR = 1, NPROD(JSV)
            GM0(IPR,JSV) = GM0(IPR,JSV) + ALPHA(NOX,IPR,JHC) * DELHC(NOX)
         ENDDO
         ENDDO

      ELSEIF ( JHC == PARENTISOP ) THEN

         !-------------------------------
         ! SOA from ISOPRENE: Parent is oxidized in
         ! gas-phase chemsitry
         !-------------------------------

         ! Get ISOP lost to rxn with OH [kg]
         !DOH = GET_DOH( I, J, L, Input_Opt )
         ! Save as DELHC (hotp 5/22/10)
         DELHC(1) = GET_DOH( I, J, L, Input_Opt, State_Chm, State_Met )

         ! Get ISOP lost to rxn with NO3 [kgC]
         ! No longer need this for isoprene SOA simulation (eam, 02/2015):
         !DELHC(2) = GET_ISOPNO3( I, J, L, Input_Opt, State_Chm, State_Met )

         ! Save diagnostic info for bug check (hotp 5/22/10)
         ! convert from kgC to kg
         DELTASOGSAVE(I,J,L,:,JHC) = DELHC(:) * 68e+0_fp/60e+0_fp

         !------------------------------------------
         ! Compute amount of semivolatile formed
         ! and add to initial SOG (hotp 7/28/10)
         !------------------------------------------
         ! update dims (hotp 5/22/10)
         DO NOX = 1, NNOX(JSV)
         DO IPR = 1, NPROD(JSV)
            GM0(IPR,JSV) = GM0(IPR,JSV) + ALPHA(NOX,IPR,JHC) * DELHC(NOX) &
                           * 68e+0_fp / 60e+0_fp ! (dkh, 11/04/05)
         ENDDO
         ENDDO

      ! Add NAP/IVOC here (hotp 5/22/10)
      ELSEIF ( JHC == PARENTBENZ .or. JHC == PARENTTOLU .or. &
               JHC == PARENTXYLE .or. JHC == PARENTNAP  ) THEN

         !-------------------------------
         ! SOA from AROMATICS
         !-------------------------------

         ! Locate IDSV (hotp 5/14/10)
         JSV = IDSV(JHC)

         ! Determine parent hydrocarbon reacted
         ! For an online calculation, GET_DARO2 can be called
         ! with 1 for high NOx, 2 for low NOx
         ! Here, we add the two pathways together and use an
         ! offline branching ratio (BETANO) (hotp 5/22/10)
         NOX = NHIGHNOX ! NOX=1
         DELHC(NOX) = ( GET_DARO2(I,J,L,1,JHC,Input_Opt,State_Chm,State_Met) + &
                        GET_DARO2(I,J,L,2,JHC,Input_Opt,State_Chm,State_Met) ) &
                      * BETANO
         NOX = NLOWNOX  ! NOX=2
         DELHC(NOX) = ( GET_DARO2(I,J,L,1,JHC,Input_Opt,State_Chm,State_Met) + &
                        GET_DARO2(I,J,L,2,JHC,Input_Opt,State_Chm,State_Met) ) &
                      * (1e+0_fp-BETANO)

         ! Determine SOG yield and add to GM0 (hotp 5/22/10)
         DO NOX = 1, NNOX(JSV)
         DO IPR = 1, NPROD(JSV)
            GM0(IPR,JSV) = GM0(IPR,JSV) + ALPHA(NOX,IPR,JHC) * DELHC(NOX)
         ENDDO
         ENDDO

         ! Diagnostic/debug info (hotp 5/22/10)
         IF ( JHC == PARENTBENZ .or. JHC == PARENTTOLU .or. &
              JHC == PARENTXYLE ) THEN
            GLOB_DARO2(I,J,L,1:2,JHC-5) = DELHC(1:2)
         ELSE ! NAP
            GLOB_DARO2(I,J,L,1:2,4) = DELHC(1:2)
         ENDIF

         ! Total SOG production diagnostic (hotp 5/18/10)
         DELTASOGSAVE(I,J,L,:,JHC)=DELHC(:)

      ! semivolpoa2: emit POA into 2 semivolatiles here (hotp 2/27/09)
      ELSEIF ( JHC == PARENTPOA ) THEN

         ! semivolpoa4opoa: DO NOTHING NOW
         !
         ! SVOCs should immediately partition upon emission
         ! SVOCs also react in the gas-phase
         ! If SVOCs were emitted here, how would you know how
         ! much to put in each phase?
         ! hotp decided to emit them after the existing SVOCs
         ! react in the gas-phase. Thus the order of operations
         ! is:
         ! SVOC + OH in gas-phase
         ! SVOC emission (added to gas-phase GM0)
         ! partitioning
         ! dry dep
         ! wet dep
         ! etc
         !
         ! DO IPR = 1, NPROD(JHC)
         ! DO NOX =1, NNOX(JHC)
         !    ! DELHC is now emission of POA
         !    DELHC(IPR) = POAEMISS(I,J,L) ! DELHC not a function of IPR
         !    DELHC(IPR) = POAEMISS(I,J,L,NOX)
         !    GM0(NOX,IPR,JHC) = GM0(NOX,IPR,JHC)
         ! &                     + ALPHA(NOX,IPR,JHC)*DELHC(IPR)
         ! ENDDO
         ! ENDDO

      ! semivolpoa4opoa: perform OPOA production (hotp 3/18/09)
      ELSEIF ( JHC == PARENTOPOA ) THEN

         ! here we oxidize gas phase POA (POG) to OPOG by reaction with OH
         ! use constant KOH = 2e-11 for now (hotp 3/18/09)
         OHMC = GET_OH( I, J, L, Input_Opt, State_Chm, State_Met )
         RK   = 2.e-11_fp * OHMC

         ! Identify IDSV (hotp 5/14/10)
         JSV = IDSV(JHC)

         DO IPR = 1, NPROD(JSV)
         DO NOX = 1, NNOX(JSV)
            ! compute loss of POG due to conversion to OPOG
            DOH = GM0(IPR,IDSV(PARENTPOA)) * (1.e+0_fp - EXP( -RK * DTCHEM) )
            DOH = MAX( DOH, 1.e-32_fp )

            ! add OPOG mass and update GM0 (ALPHA=1)
            GM0(IPR,JSV) = GM0(IPR,JSV) + ALPHA(NOX,IPR,JHC) * DOH

            ! update POG mass
            GM0(IPR,IDSV(PARENTPOA)) = GM0(IPR,IDSV(PARENTPOA)) - DOH

            ! check (hotp 10/11/09)
            GM0(IPR,IDSV(PARENTPOA)) = MAX( GM0(IPR,IDSV(PARENTPOA)), 1e-32_fp)

            ! diagnostic information (hotp 3/28/09)
            GLOB_POGRXN(I,J,L,IPR) = DOH

            ! Total SOG production diagnostic (hotp 5/18/10)
            ! Caution: the 4th index is actually NOX, but we use
            ! IPR here
            DELTASOGSAVE(I,J,L,IPR,JHC) = DOH

         ENDDO
         ENDDO
      ENDIF
   ENDDO  ! JHC

   ! semivolpoa4opoa: emit POA last (after OPOA formation) (hotp 3/18/09)
   ! SVOC emissions are added to gas-phase GM0
   IF ( id_POA1 > 0 ) THEN
      JHC = PARENTPOA

      ! Use IDSV (hotp 5/14/10)
      JSV = IDSV(JHC)

      DO IPR = 1, NPROD(JSV)
      DO NOX = 1, NNOX(JSV)   ! update dims (hotp 5/22/10)
         ! DELHC is now emission of SVOC (POG1 + POG2)
         DELHC(IPR)   = POAEMISS(I,J,L,1) + POAEMISS(I,J,L,2)
         GM0(IPR,JSV) = GM0(IPR,JSV) + ALPHA(NOX,IPR,JHC)*DELHC(IPR)

         ! Total SOG production diagnostic (hotp 5/18/10)
         ! Caution: the 4th index is actually NOX, but we use
         ! IPR here
         DELTASOGSAVE(I,J,L,IPR,JHC) = DELHC(IPR)

      ENDDO
      ENDDO
   ENDIF

   !=================================================================
   ! Store Hydrocarbon remaining after oxidation rxn back into Spc
   !=================================================================
   ! Nothing to do for isoprene or aromatics here,  as their oxidation
   ! is treated online.
   ! The same now applies to MTPA and LIMO. As of v11-02c, their
   ! oxidation is treated online (mps, 9/7/17)
   !Spc(I,J,L,id_MTPA) = MAX( NMVOC(1), 1.e-32_fp )
   !Spc(I,J,L,id_LIMO) = MAX( NMVOC(2), 1.e-32_fp )
   Spc(I,J,L,id_MTPO) = MAX( NMVOC(3), 1.e-32_fp )
   ORVC_SESQ(I,J,L)   = MAX( NMVOC(4), 1.e-32_fp )

   ! Free pointer
   Spc => NULL()

 END SUBROUTINE CHEM_NVOC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: soa_partition
!
! !DESCRIPTION: Subroutine SOA\_PARTITION partitions the mass of gas and
!  aerosol species according to five Hydrocarbon species and three oxidants.
!  (rjp, bmy, 7/7/04, 5/22/06)
!\\
!\\
!  Revised purpose: SOA\_PARTITION assigns the mass in the chemical
!  species array to the GM0 and AM0 arrays (hotp 5/13/10)
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE SOA_PARTITION( I, J, L, GM0, AM0, State_Chm )
!
! !USES:
!
   USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
   INTEGER,  INTENT(IN)  :: I              ! Longitude index
   INTEGER,  INTENT(IN)  :: J              ! Latitude index
   INTEGER,  INTENT(IN)  :: L              ! Altitude index
!
! !OUTPUT PARAMETERS:
!
   REAL(fp), INTENT(OUT) :: GM0(MPROD,MSV) ! Gas mass for HCs and
                                           !  oxidation products [kg]
   REAL(fp), INTENT(OUT) :: AM0(MPROD,MSV) ! Aer mass for HCs and
                                           !  oxidation products [kg]
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
!
! !REMARKS:
!  NOTE: GPROD and APROD are mass ratios of individual oxidation
!        products of gas/aerosol to the sum of all.
!
! !REVISION HISTORY:
!  13 May 2010 - H.O.T. Pye - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Scalars
   INTEGER           :: JHC, IPR, NOX, JSV

   ! Pointers
   REAL(fp), POINTER :: Spc(:,:,:,:)

   !=================================================================
   ! SOA_PARTITION begins here!
   !=================================================================

   ! Point to the chemical species array [kg]
   Spc => State_Chm%Species

   ! Initialize everything to zero (hotp 5/17/10)
   GM0 = 0e+0_fp
   AM0 = 0e+0_fp

   !---------------------------------------
   ! SEMIVOLATILE 1: MTPA, LIMO, MTPO, SESQ
   ! hotp 5/21/10
   !---------------------------------------
   JHC = PARENTMTPA
   JSV = IDSV(JHC)
   ! gas phase
   GM0(1,JSV)=Spc(I,J,L,id_TSOG1) ! C* =   1
   GM0(2,JSV)=Spc(I,J,L,id_TSOG2) ! C* =  10
   GM0(3,JSV)=Spc(I,J,L,id_TSOG3) ! C* = 100
   GM0(4,JSV)=Spc(I,J,L,id_TSOG0) ! C* =   0.1
   ! aerosol phase
   AM0(1,JSV)=Spc(I,J,L,id_TSOA1)
   AM0(2,JSV)=Spc(I,J,L,id_TSOA2)
   AM0(3,JSV)=Spc(I,J,L,id_TSOA3)
   AM0(4,JSV)=Spc(I,J,L,id_TSOA0)

   !---------------------------------------------------------------------------
   ! Prior to 7/15/19:
   ! Remove isoprene from VBS (mps, 7/15/19)
   !!---------------------------------------
   !! SEMIVOLATILE 2: ISOP
   !!---------------------------------------
   !JHC = PARENTISOP
   !JSV = IDSV(JHC)
   !! gas phase
   !GM0(1,JSV)=Spc(I,J,L,id_ISOG1)
   !GM0(2,JSV)=Spc(I,J,L,id_ISOG2)
   !GM0(3,JSV)=Spc(I,J,L,id_ISOG3)
   !! aerosol phase
   !AM0(1,JSV)=Spc(I,J,L,id_ISOA1)
   !AM0(2,JSV)=Spc(I,J,L,id_ISOA2)
   !AM0(3,JSV)=Spc(I,J,L,id_ISOA3)
   !---------------------------------------------------------------------------

   !---------------------------------------
   ! SEMIVOLATILE 3: BENZ, TOLU, XYLE, NAP
   !---------------------------------------
   ! Lumped arom/IVOC/NAP semivolatiles (hotp 5/13/10)
   JHC = PARENTBENZ ! IDSV(B)=IDSV(T)=IDSV(X)=IDSV(N)
   JSV = IDSV(JHC)
   ! gas phase
   GM0(1,JSV)=Spc(I,J,L,id_ASOG1)
   GM0(2,JSV)=Spc(I,J,L,id_ASOG2)
   GM0(3,JSV)=Spc(I,J,L,id_ASOG3)
   ! aerosol phase
   AM0(1,JSV)=Spc(I,J,L,id_ASOA1)
   AM0(2,JSV)=Spc(I,J,L,id_ASOA2)
   AM0(3,JSV)=Spc(I,J,L,id_ASOA3)
   AM0(4,JSV)=Spc(I,J,L,id_ASOAN)

   !---------------------------------------
   ! SEMIVOLATILE 4: POA/SVOCs
   !---------------------------------------
   ! POA-Primary SVOCs (hotp 5/13/10)
   JHC = PARENTPOA
   JSV = IDSV(JHC)
   IF ( id_POA1 > 0 .and. id_POA2 > 0 .and. &
        id_POG1 > 0 .and. id_POG2 > 0 ) THEN
      ! gas phase
      GM0(1,JSV) = Spc(I,J,L,id_POG1)
      GM0(2,JSV) = Spc(I,J,L,id_POG2)
      ! aerosol phase
      AM0(1,JSV) = Spc(I,J,L,id_POA1)
      AM0(2,JSV) = Spc(I,J,L,id_POA2)
   ENDIF

   !---------------------------------------
   ! SEMIVOLATILE 5: OPOA/O-SVOCs
   !---------------------------------------
   ! OPOA-Oxidized SVOCs (hotp 5/13/10)
   JHC = PARENTOPOA
   JSV = IDSV(JHC)
   IF ( id_OPOA1 > 0 .and. id_OPOA2 > 0 .and. &
        id_OPOG1 > 0 .and. id_OPOG2 > 0 ) THEN
      ! gas phase
      GM0(1,JSV) = Spc(I,J,L,id_OPOG1)
      GM0(2,JSV) = Spc(I,J,L,id_OPOG2)
      ! aerosol phase
      AM0(1,JSV) = Spc(I,J,L,id_OPOA1)
      AM0(2,JSV) = Spc(I,J,L,id_OPOA2)
   ENDIF

   ! Free pointer
   Spc => NULL()

 END SUBROUTINE SOA_PARTITION
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: soa_lump
!
! !DESCRIPTION: Subroutine SOA\_LUMP returns the organic gas and aerosol back
!  to the STT array.  (rjp, bmy, 7/7/04, 2/6/07)
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE SOA_LUMP( I, J, L, GM0, AM0, State_Chm, State_Diag )
!
! !USES:
!
   USE State_Chm_Mod,      ONLY : ChmState
   USE State_Diag_Mod,     ONLY : DgnState
!
! !INPUT PARAMETERS:
!
   INTEGER,        INTENT(IN)    :: I              ! Longitude index
   INTEGER,        INTENT(IN)    :: J              ! Latitude index
   INTEGER,        INTENT(IN)    :: L              ! Altitude index
   REAL(fp),       INTENT(IN)    :: GM0(MPROD,MSV) ! Gas mass for HCs and
                                                   !  oxidation products [kg]
   REAL(fp),       INTENT(IN)    :: AM0(MPROD,MSV) ! Aer mass for HCs and
                                                   !  oxidation products [kg]
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(ChmState), INTENT(INOUT) :: State_Chm      ! Chemistry State object
   TYPE(DgnState), INTENT(INOUT) :: State_Diag     ! Diagnostics State obj
!
! !REVISION HISTORY:
!  07 Jul 2004 - R. Park - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER           :: JHC, IPR, NOX, JSV ! JSV (hotp 5/13/10)
   REAL(fp)          :: GASMASS, AERMASS
   INTEGER           :: id_SPECIES ! hotp 6/5/10
   REAL(fp)          :: AERCHANGE  ! hotp 6/5/10

   ! Pointers
   REAL(fp), POINTER :: Spc(:,:,:,:)

   !=================================================================
   ! SOA_LUMP begins here!
   !=================================================================

   ! Point to the chemical species array [kg]
   Spc => State_Chm%Species

   !=================================================================
   ! Semivolatile Group 1: monoterpenes and sesquiterpenes (hotp 5/22/10)
   !=================================================================

   ! Initialize
   GASMASS = 0e+0_fp
   AERMASS = 0e+0_fp
   JHC     = PARENTMTPA
   JSV     = IDSV(JHC)

   ! Save diagnostic info
   DO IPR = 1, NPROD(JSV) ! change JHC to JSV
      GASMASS = GASMASS + GM0(IPR,JSV)
      AERMASS = AERMASS + AM0(IPR,JSV)
   ENDDO

   !-----------------------------
   ! Transient mass bal prod/evap
   ! (hotp 6/5/10)
   !-----------------------------
   id_SPECIES = id_TSOA1
   IPR = 1
   AERCHANGE = AM0(IPR,JSV) - Spc(I,J,L,id_SPECIES)
   IF ( AERCHANGE .GT. 0e+0_fp ) THEN
      SPECSOAPROD(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAPROD(I,J,L,IPR,JSV)
   ELSE
      SPECSOAEVAP(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAEVAP(I,J,L,IPR,JSV)
   ENDIF

   id_SPECIES = id_TSOA2
   IPR = 2
   AERCHANGE = AM0(IPR,JSV) - Spc(I,J,L,id_SPECIES)
   IF ( AERCHANGE .GT. 0e+0_fp ) THEN
      SPECSOAPROD(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAPROD(I,J,L,IPR,JSV)
   ELSE
      SPECSOAEVAP(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAEVAP(I,J,L,IPR,JSV)
   ENDIF

   id_SPECIES = id_TSOA3
   IPR = 3
   AERCHANGE = AM0(IPR,JSV) - Spc(I,J,L,id_SPECIES)
   IF ( AERCHANGE .GT. 0e+0_fp ) THEN
      SPECSOAPROD(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAPROD(I,J,L,IPR,JSV)
   ELSE
      SPECSOAEVAP(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAEVAP(I,J,L,IPR,JSV)
   ENDIF

   ! Add C*=0.1 product (hotp 6/12/10)
   id_SPECIES = id_TSOA0
   IPR = 4
   AERCHANGE = AM0(IPR,JSV) - Spc(I,J,L,id_SPECIES)
   IF ( AERCHANGE .GT. 0e+0_fp ) THEN
      SPECSOAPROD(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAPROD(I,J,L,IPR,JSV)
   ELSE
      SPECSOAEVAP(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAEVAP(I,J,L,IPR,JSV)
   ENDIF

   !-----------------------------
   ! Update species [kg]
   !-----------------------------
   ! gas phase
   Spc(I,J,L,id_TSOG1) = MAX( GM0(1,JSV), 1e-32_fp )
   Spc(I,J,L,id_TSOG2) = MAX( GM0(2,JSV), 1e-32_fp )
   Spc(I,J,L,id_TSOG3) = MAX( GM0(3,JSV), 1e-32_fp )
   Spc(I,J,L,id_TSOG0) = MAX( GM0(4,JSV), 1e-32_fp )
   ! aerosol phase
   Spc(I,J,L,id_TSOA1) = MAX( AM0(1,JSV), 1e-32_fp )
   Spc(I,J,L,id_TSOA2) = MAX( AM0(2,JSV), 1e-32_fp )
   Spc(I,J,L,id_TSOA3) = MAX( AM0(3,JSV), 1e-32_fp )
   Spc(I,J,L,id_TSOA0) = MAX( AM0(4,JSV), 1e-32_fp )

   !=================================================================
   ! Semivolatile Group 2: isoprene (hotp 5/22/10)
   !=================================================================

   ! Initialize
   GASMASS = 0e+0_fp
   AERMASS = 0e+0_fp
   JHC = PARENTISOP
   JSV = IDSV(JHC)

   ! Save diagnostic info
   DO IPR = 1, NPROD(JSV) ! change JHC to JSV
      GASMASS = GASMASS + GM0(IPR,JSV)
      AERMASS = AERMASS + AM0(IPR,JSV)
   ENDDO

   !---------------------------------------------------------------------------
   ! Prior to 7/15/19:
   ! Remove isoprene from VBS (mps, 7/15/19)
   !!-----------------------------
   !! Transient mass bal prod/evap
   !! (hotp 6/5/10)
   !!-----------------------------
   !id_SPECIES = id_ISOA1
   !IPR = 1
   !AERCHANGE = AM0(IPR,JSV) - Spc(I,J,L,id_SPECIES)
   !IF ( AERCHANGE .GT. 0e+0_fp ) THEN
   !   SPECSOAPROD(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAPROD(I,J,L,IPR,JSV)
   !ELSE
   !   SPECSOAEVAP(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAEVAP(I,J,L,IPR,JSV)
   !ENDIF
   !
   !id_SPECIES = id_ISOA2
   !IPR = 2
   !AERCHANGE = AM0(IPR,JSV) - Spc(I,J,L,id_SPECIES)
   !IF ( AERCHANGE .GT. 0e+0_fp ) THEN
   !   SPECSOAPROD(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAPROD(I,J,L,IPR,JSV)
   !ELSE
   !   SPECSOAEVAP(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAEVAP(I,J,L,IPR,JSV)
   !ENDIF
   !
   !id_SPECIES = id_ISOA3
   !IPR = 3
   !AERCHANGE = AM0(IPR,JSV) - Spc(I,J,L,id_SPECIES)
   !IF ( AERCHANGE .GT. 0e+0_fp ) THEN
   !   SPECSOAPROD(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAPROD(I,J,L,IPR,JSV)
   !ELSE
   !   SPECSOAEVAP(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAEVAP(I,J,L,IPR,JSV)
   !ENDIF
   !
   !!-----------------------------
   !! Update species [kg]
   !!-----------------------------
   !! gas phase
   !Spc(I,J,L,id_ISOG1) = MAX( GM0(1,JSV), 1e-32_fp )
   !Spc(I,J,L,id_ISOG2) = MAX( GM0(2,JSV), 1e-32_fp )
   !Spc(I,J,L,id_ISOG3) = MAX( GM0(3,JSV), 1e-32_fp )
   !! aerosol phase
   !Spc(I,J,L,id_ISOA1) = MAX( AM0(1,JSV), 1e-32_fp )
   !Spc(I,J,L,id_ISOA2) = MAX( AM0(2,JSV), 1e-32_fp )
   !Spc(I,J,L,id_ISOA3) = MAX( AM0(3,JSV), 1e-32_fp )
   !---------------------------------------------------------------------------

   !=================================================================
   ! Semivolatile Group 3: benzene, toluene, xylene, naphthalene/IVOC
   ! Lump of products of 7-9 Hydrocarbon class (aromatics) (dkh, 11/11/06)
   ! Lumped aromatic/IVOC (hotp 5/13/10
   !=================================================================

   ! Initialize
   GASMASS = 0e+0_fp
   AERMASS = 0e+0_fp
   JHC = PARENTBENZ
   JSV = IDSV(JHC)

   ! Save diagnostic info
   ! This is a lumped species (hotp 5/13/10)
   DO IPR = 1, NPROD(JSV) ! change JHC to JSV
      GASMASS = GASMASS + GM0(IPR,JSV)
      AERMASS = AERMASS + AM0(IPR,JSV)
   ENDDO

   !-----------------------------
   ! Transient mass bal prod/evap
   ! (hotp 6/5/10)
   !-----------------------------
   id_SPECIES = id_ASOA1
   IPR = 1
   AERCHANGE = AM0(IPR,JSV) - Spc(I,J,L,id_SPECIES)
   IF ( AERCHANGE .GT. 0e+0_fp ) THEN
      SPECSOAPROD(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAPROD(I,J,L,IPR,JSV)
   ELSE
      SPECSOAEVAP(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAEVAP(I,J,L,IPR,JSV)
   ENDIF

   id_SPECIES = id_ASOA2
   IPR = 2
   AERCHANGE = AM0(IPR,JSV) - Spc(I,J,L,id_SPECIES)
   IF ( AERCHANGE .GT. 0e+0_fp ) THEN
      SPECSOAPROD(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAPROD(I,J,L,IPR,JSV)
   ELSE
      SPECSOAEVAP(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAEVAP(I,J,L,IPR,JSV)
   ENDIF

   id_SPECIES = id_ASOA3
   IPR = 3
   AERCHANGE = AM0(IPR,JSV) - Spc(I,J,L,id_SPECIES)
   IF ( AERCHANGE .GT. 0e+0_fp ) THEN
      SPECSOAPROD(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAPROD(I,J,L,IPR,JSV)
   ELSE
      SPECSOAEVAP(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAEVAP(I,J,L,IPR,JSV)
   ENDIF

   id_SPECIES = id_ASOAN
   IPR = 4
   AERCHANGE = AM0(IPR,JSV) - Spc(I,J,L,id_SPECIES)
   IF ( AERCHANGE .GT. 0e+0_fp ) THEN
      SPECSOAPROD(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAPROD(I,J,L,IPR,JSV)
   ELSE
      SPECSOAEVAP(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAEVAP(I,J,L,IPR,JSV)
   ENDIF

   !-----------------------------
   ! Update species [kg]
   !-----------------------------
   ! APROD and GPROD are no longer used, but GM0 and AM0
   ! need to be saved to species arrays (hotp 5/13/10)
   ! HIGH NOX  ! update dims (hotp 5/22/10)
   !NOX = NHIGHNOX
   ! gas phase
   Spc(I,J,L,id_ASOG1) = MAX( GM0(1,JSV), 1e-32_fp )
   Spc(I,J,L,id_ASOG2) = MAX( GM0(2,JSV), 1e-32_fp )
   Spc(I,J,L,id_ASOG3) = MAX( GM0(3,JSV), 1e-32_fp )
   ! aerosol phase
   Spc(I,J,L,id_ASOA1) = MAX( AM0(1,JSV), 1e-32_fp )
   Spc(I,J,L,id_ASOA2) = MAX( AM0(2,JSV), 1e-32_fp )
   Spc(I,J,L,id_ASOA3) = MAX( AM0(3,JSV), 1e-32_fp )
   ! LOW NOX (only 1 aerosol phase)
   !NOX = NLOWNOX ! store in spot 4 (hotp 5/22/10)
   Spc(I,J,L,id_ASOAN) = MAX( AM0(4,JSV), 1e-32_fp )

   !=================================================================
   ! Semivolatile 4: POA/primary SVOCs
   ! Lump of products of 10th Hydrocarbon class (POA)
   ! semivolpoa2: lump POA (hotp 2/27/09)
   !=================================================================
   IF ( id_POA1 > 0 .and. id_POA2 > 0 .and. &
        id_POG1 > 0 .and. id_POG2 > 0 ) THEN

      ! Initialize
      !JHC     = 10
      GASMASS = 0e+0_fp
      AERMASS = 0e+0_fp
      JHC     = PARENTPOA
      JSV     = IDSV(JHC)

      ! Replace JHC with JSV (hotp 5/13/10)
      DO IPR = 1, NPROD(JSV)
         GASMASS = GASMASS + GM0(IPR,JSV)
         AERMASS = AERMASS + AM0(IPR,JSV)
      ENDDO

      !---------------------------
      ! Transient SOA PROD/EVAP
      ! (hotp 6/5/10)
      !---------------------------
      id_SPECIES = id_POA1
      IPR = 1
      AERCHANGE = AM0(IPR,JSV) - Spc(I,J,L,id_SPECIES)
      IF ( AERCHANGE .GT. 0e+0_fp ) THEN
         SPECSOAPROD(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAPROD(I,J,L,IPR,JSV)
      ELSE
         SPECSOAEVAP(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAEVAP(I,J,L,IPR,JSV)
      ENDIF

      id_SPECIES = id_POA2
      IPR = 2
      AERCHANGE = AM0(IPR,JSV) - Spc(I,J,L,id_SPECIES)
      IF ( AERCHANGE .GT. 0e+0_fp ) THEN
         SPECSOAPROD(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAPROD(I,J,L,IPR,JSV)
      ELSE
         SPECSOAEVAP(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAEVAP(I,J,L,IPR,JSV)
      ENDIF

      !---------------------------
      ! Update species [kg]
      !---------------------------
      ! gas phase
      Spc(I,J,L,id_POG1) = MAX( GM0(1,JSV), 1.e-32_fp )
      Spc(I,J,L,id_POG2) = MAX( GM0(2,JSV), 1.e-32_fp )
      ! aerosol phase
      Spc(I,J,L,id_POA1) = MAX( AM0(1,JSV), 1.e-32_fp )
      Spc(I,J,L,id_POA2) = MAX( AM0(2,JSV), 1.e-32_fp )

   ENDIF ! POA

   !=================================================================
   ! Semivolatile 5: OPOA/oxidized primary SVOCs
   ! Lump of products of 11th Hydrocarbon class (OPOA)
   ! semivolpoa4opoa: lump OPOA (hotp 2/27/09)
   !=================================================================
   IF ( id_OPOA1 > 0 .and. id_OPOA2 > 0 .and. &
        id_OPOG1 > 0 .and. id_OPOG2 > 0 ) THEN

      ! Initialize
      GASMASS = 0e+0_fp
      AERMASS = 0e+0_fp
      JHC     = PARENTOPOA
      JSV     = IDSV(JHC)

      ! Save diagnostic info
      DO IPR = 1, NPROD(JSV)
         GASMASS = GASMASS + GM0(IPR,JSV)
         AERMASS = AERMASS + AM0(IPR,JSV)
      ENDDO

      !---------------------------
      ! Transient SOA PROD/EVAP
      ! (hotp 6/5/10)
      !---------------------------
      id_SPECIES = id_OPOA1
      IPR = 1
      AERCHANGE = AM0(IPR,JSV) - Spc(I,J,L,id_SPECIES)
      IF ( AERCHANGE .GT. 0e+0_fp ) THEN
         SPECSOAPROD(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAPROD(I,J,L,IPR,JSV)
      ELSE
         SPECSOAEVAP(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAEVAP(I,J,L,IPR,JSV)
      ENDIF

      id_SPECIES = id_OPOA2
      IPR = 2
      AERCHANGE = AM0(IPR,JSV) - Spc(I,J,L,id_SPECIES)
      IF ( AERCHANGE .GT. 0e+0_fp ) THEN
         SPECSOAPROD(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAPROD(I,J,L,IPR,JSV)
      ELSE
         SPECSOAEVAP(I,J,L,IPR,JSV) = AERCHANGE + SPECSOAEVAP(I,J,L,IPR,JSV)
      ENDIF

      !---------------------------
      ! Update species [kg]
      !---------------------------
      ! gas phase
      Spc(I,J,L,id_OPOG1) = MAX( GM0(1,JSV), 1.e-32_fp )
      Spc(I,J,L,id_OPOG2) = MAX( GM0(2,JSV), 1.e-32_fp )
      ! aerosol phase
      Spc(I,J,L,id_OPOA1) = MAX( AM0(1,JSV), 1.e-32_fp )
      Spc(I,J,L,id_OPOA2) = MAX( AM0(2,JSV), 1.e-32_fp )

   ENDIF ! OPOA

   ! Free pointer
   Spc => NULL()

 END SUBROUTINE SOA_LUMP
!EOC
#ifdef TOMAS
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
!! !IROUTINE: emitsgc
!
! !DESCRIPTION: Subroutine EMITSGC calculates sub-grid coagulation for the size
!  distribution of emission. (win, 10/6/07)
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE EMITSGC( Input_Opt, State_Chm, State_Grid, State_Met, &
                     EMISMASS,  CTYPE )
!
! !USES:
!
#ifdef BPCH_DIAG
   USE CMN_DIAG_MOD             ! ND59
   USE DIAG_MOD,           ONLY : AD59_ECIL,   AD59_ECOB
   USE DIAG_MOD,           ONLY : AD59_OCIL,   AD59_OCOB
   USE DIAG_MOD,           ONLY : AD59_NUMB
#endif
   USE ERROR_MOD,          ONLY : IT_IS_NAN
   USE Input_Opt_Mod,      ONLY : OptInput
   USE State_Chm_Mod,      ONLY : ChmState
   USE State_Grid_Mod,     ONLY : GrdState
   USE State_Met_Mod,      ONLY : MetState
   USE TOMAS_MOD,          ONLY : IBINS,    AVGMASS,  ICOMP,   IDIAG
   USE TOMAS_MOD,          ONLY : SRTECIL,  SRTECOB,  SRTOCIL
   USE TOMAS_MOD,          ONLY : SRTOCOB,  SRTSO4,   SRTNH4
   USE TOMAS_MOD,          ONLY : SRTH2O,   MNFIX
   USE TOMAS_MOD,          ONLY : SUBGRIDCOAG
   USE TOMAS_MOD,          ONLY : SGCTSCALE
   USE TOMAS_MOD,          ONLY : NH4BULKTOBIN
!
! !INPUT PARAMETERS:
!
   INTEGER,        INTENT(IN) :: CTYPE       ! 1 = EC and 2 = OC
   TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
   TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
   TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
   REAL(fp),       INTENT(IN) :: EMISMASS(State_Grid%NX,State_Grid%NY,IBINS)
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  06 Oct 2007 - W. Trivitayanurak - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   REAL(fp)  :: NDISTINIT(IBINS)
   REAL(fp)  :: NDISTFINAL(IBINS)
   REAL(fp)  :: MADDFINAL(IBINS)
   REAL(fp)  :: NDIST(IBINS)
   REAL(fp)  :: MDIST(IBINS,ICOMP)
   REAL(fp)  :: NDIST2(IBINS)
   REAL(fp)  :: MDIST2(IBINS,ICOMP)
   REAL*4    :: BOXVOL, TEMP, PRES
   INTEGER   :: I, J, L, K, C, PBL_MAX
   REAL(fp)  :: F_OF_PBL
   LOGICAL   :: ERRORSWITCH, PDBUG
   REAL*4    :: N0(State_Grid%NZ,IBINS)
   REAL*4    :: N1(State_Grid%NZ,IBINS)
   REAL*4    :: MIL0(State_Grid%NZ,IBINS)
   REAL*4    :: MIL1(State_Grid%NZ,IBINS)
   REAL*4    :: MOB0(State_Grid%NZ,IBINS)
   REAL*4    :: MOB1(State_Grid%NZ,IBINS)
   INTEGER   :: ii, jj, ll
   LOGICAL   :: dbg = .false.
   DATA ii, jj, ll /53, 29, 8 /

   ! Pointers
   REAL(fp), POINTER :: Spc(:,:,:,:)

   !=================================================================
   ! EMITSGC begins here!
   !=================================================================

   IF ( Input_Opt%LNLPBL ) THEN
      print *,'Currently subroutine EMITSGC does not support ', &
              'the new non-local PBL scheme!'
      stop
   ENDIF

   ! Point to the chemical species array [kg]
   Spc => State_Chm%Species

   ! Maximum extent of PBL [model levels]
   PBL_MAX = State_Met%PBL_MAX_L

   !temp debug      if( sum(emismass(ii,jj,:)) > 0e+0_fp) dbg = .true.
   if( dbg) then
      print *,'===== Entering EMITSGC ===== at',ii,jj,ll
      print *,'Nk'
      print *,Spc(ii,jj,ll,id_NK1:id_NK1+ibins-1)
      print *,'Mk'
      do k=1,icomp-idiag
         print *,'comp',k
         print *,Spc(ii,jj,ll,id_NK1+k*IBINS:id_NK1+IBINS-1+k*IBINS)
      enddo
      print *,'EMISSION'
      print *,emismass(ii,jj,:)
   endif
   !temp debug --------

   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX
      IF ( SUM( EMISMASS(I,J,:) ) == 0.e+0_fp ) GOTO 100

      DO L = 1, PBL_MAX
         ! Fraction of PBL spanned by grid box (I,J,L) [unitless]
         F_OF_PBL = State_Met%F_OF_PBL(I,J,L)

         DO K = 1, IBINS
            NDISTINIT(K) = EMISMASS(I,J,K) * F_OF_PBL / AVGMASS(K)
            NDIST(K) = Spc(I,J,L,id_NK1+K-1)
            DO C = 1, ICOMP-IDIAG
               MDIST(K,C) = Spc(I,J,L,id_NK1+IBINS*C+K-1)
               IF( IT_IS_NAN( MDIST(K,C) ) ) THEN
                  PRINT *,'+++++++ Found NaN in EMITSGC ++++++++'
                  PRINT *,'Location (I,J,L):',I,J,L,'Bin',K,'comp',C
               ENDIF
            ENDDO
            MDIST(K,SRTH2O) = Spc(I,J,L,id_AW1-1+K)
            NDISTFINAL(K) = 0e+0_fp
            MADDFINAL(K) = 0e+0_fp
         ENDDO

         IF ( SRTNH4 > 0 ) THEN
            CALL NH4BULKTOBIN( MDIST(:,SRTSO4),   &
                               Spc(I,J,L,id_NH4), &
                               MDIST(:,SRTNH4) )
         ENDIF

         ! Save initial info for diagnostic
         N0(L,:) = NDIST(:)
         IF(CTYPE == 1) THEN
            MIL0(L,:) = MDIST(:,SRTECIL)
            MOB0(L,:) = MDIST(:,SRTECOB)
         ELSE
            MIL0(L,:) = MDIST(:,SRTOCIL)
            MOB0(L,:) = MDIST(:,SRTOCOB)
         ENDIF

         BOXVOL  = State_Met%AIRVOL(I,J,L) * 1.e6 !convert from m3 -> cm3
         TEMP    = State_Met%T(I,J,L)
         PRES    = State_Met%PMID(i,j,l)*100.0 ! in Pa

         PDBUG = .FALSE.
         !temp debug
         if( dbg .and. i==ii .and. j==jj .and. l==ll ) then
            print *,'===== NDISTINIT ===== at',ii,jj,ll
            print *, NDISTINIT(:)
         endif
         !temp debug  if( dbg .and. i==ii .and. j==jj .and. l==ll ) PDBUG = .TRUE.

         CALL SUBGRIDCOAG( NDISTINIT, NDIST, MDIST, BOXVOL,TEMP, &
                           PRES, SGCTSCALE, NDISTFINAL, MADDFINAL,pdbug)
         IF ( PDBUG ) THEN
            PRINT *,'Found error in SUBGRIDCOAG at', I,J,L
            PRINT *,'Nk',Spc(I,J,L,id_NK1:id_NK1+ibins-1)
            do k=1,8
               print *,'Mk comp',k
               print *,Spc(I,J,L,id_NK1+k*IBINS:id_NK1+IBINS-1+k*IBINS)
            enddo
         ENDIF

         DO K = 1, IBINS
            NDIST(K) = NDIST(K) + NDISTFINAL(K)
            IF( CTYPE == 1 ) THEN
               MDIST(K,SRTECIL) = MDIST(K,SRTECIL) + &
                                  NDISTFINAL(K) * AVGMASS(K) * 0.2e+0_fp + &
                                  MADDFINAL(K) * 0.2e+0_fp
               MDIST(K,SRTECOB) = MDIST(K,SRTECOB) + &
                                  NDISTFINAL(K) * AVGMASS(K) * 0.8e+0_fp + &
                                  MADDFINAL(K) * 0.8e+0_fp
            ELSE
               MDIST(K,SRTOCIL) = MDIST(K,SRTOCIL) + &
                                  NDISTFINAL(K) * AVGMASS(K) * 0.5e+0_fp + &
                                  MADDFINAL(K) * 0.5e+0_fp
               MDIST(K,SRTOCOB) = MDIST(K,SRTOCOB) + &
                                  NDISTFINAL(K) * AVGMASS(K) * 0.5e+0_fp + &
                                  MADDFINAL(K) * 0.5e+0_fp
            ENDIF
         ENDDO
         !temp debug
         if( dbg .and. i==ii .and. j==jj .and. l==ll ) then
            print *,'===== After SUBGRIDCOAG ===== at',ii,jj,ll
            print *,'Nk'
            print *, NDIST(:)
            print *,'xxx___NDISTFINAL__xxx'
            print *, NDISTFINAL(:)

            print *,'Mk'
            do k=1,icomp
               print *,'comp',k
               print *,MDIST(:,k)
            enddo
         endif
         !temp debug --------

         ! Fix any inconsistencies in size dist
         DO K= 1, IBINS
            NDIST2(K) = NDIST(K)
            DO C = 1, ICOMP
               MDIST2(K,C) = MDIST(K,C)
            ENDDO
         ENDDO

         ERRORSWITCH = .FALSE.

         CALL MNFIX( NDIST2, MDIST2, ERRORSWITCH )

         IF( ERRORSWITCH ) PRINT *,'EMITSGC: MNFIX found error ', &
                                   'after SUBGRIDCOAG at ',I,J,L

         DO K = 1, IBINS
            Spc(I,J,L,id_NK1-1+K) = NDIST2(K)
            DO C = 1, ICOMP-IDIAG
               Spc(I,J,L,id_NK1+K-1+C*IBINS) = MDIST2(K,C)
            ENDDO
            Spc(I,J,L,id_AW1-1+K)  = MDIST2(K,SRTH2O)
         ENDDO

         ! Save final info for diagnostic
         N1(L,:) = NDIST2(:)
         IF(CTYPE == 1) THEN
            MIL1(L,:) = MDIST2(:,SRTECIL)
            MOB1(L,:) = MDIST2(:,SRTECOB)
         ELSE
            MIL1(L,:) = MDIST2(:,SRTOCIL)
            MOB1(L,:) = MDIST2(:,SRTOCOB)
         ENDIF

      ENDDO ! L loop

      !=======================================================================
      !  ND59 Diagnostic: Size-resolved primary emission in
      !                 [kg/box/timestep] and the corresponding
      !                  number emission [no./box/timestep]
      !=======================================================================
#ifdef BPCH_DIAG
      IF ( ND59 > 0 ) THEN
         DO L = 1, PBL_MAX
         DO K = 1, IBINS
            SELECT CASE (CTYPE)
            CASE (1)
               AD59_ECIL(I,J,1,K) = AD59_ECIL(I,J,1,K) + &
                                    ( MIL1(L,K) - MIL0(L,K) )
               AD59_ECOB(I,J,1,K) = AD59_ECOB(I,J,1,K) + &
                                    ( MOB1(L,K) - MOB0(L,K) )
            CASE (2)
               AD59_OCIL(I,J,1,K) = AD59_OCIL(I,J,1,K) + &
                                    ( MIL1(L,K) - MIL0(L,K) )
               AD59_OCOB(I,J,1,K) = AD59_OCOB(I,J,1,K) + &
                                    ( MOB1(L,K) - MOB0(L,K) )
            END SELECT
            AD59_NUMB(I,J,1,K) = AD59_NUMB(I,J,1,K) + ( N1(L,K) - N0(L,K) )
         ENDDO
         ENDDO
      ENDIF
#endif

100   CONTINUE

   ENDDO ! I loop
   ENDDO ! J loop

   ! Free pointer
   Spc => NULL()

 END SUBROUTINE EMITSGC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: scalecarb
!
! !DESCRIPTION: Function SCALECARB split the carbonaceous emission from each
!  source into the TOMAS aerosol size bins using different mass distribution
!  for fossil fuel and biomass burning.  The mass size distributions
! are different for EC and OC. (win, 9/4/07)
!\\
!\\
! !INTERFACE:
!
 FUNCTION SCALECARB( State_Grid, BULKEMIS, STYPE, CTYPE ) &
      RESULT( tomasdist )
!
! !USES:
!
   USE State_Grid_Mod, ONLY : GrdState
   USE TOMAS_MOD,      ONLY : IBINS
   USE TOMAS_MOD,      ONLY : OCSCALE30,  ECSCALE30
   USE TOMAS_MOD,      ONLY : OCSCALE100, ECSCALE100
!
! !INPUT PARAMETERS:
!
   TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
   REAL(fp),       INTENT(IN) :: BULKEMIS(State_Grid%NX,State_Grid%NY)
   INTEGER,        INTENT(IN) :: STYPE
   INTEGER,        INTENT(IN) :: CTYPE
!
! !RETURN VALUE:
!
   REAL(fp) :: tomasdist(State_Grid%NX, State_Grid%NY, IBINS)
!
! !REMARKS:
!    STYPE (source type): 1 = Fossil fuel
!                         3 = Biomass burning
!    CTYPE (carbon type): 1 = EC
!                         2 = OC
!                                                                              .
!  Array ECSCALE30 and OCSCALE100 specify how mass is distributed into bins
!  for a 30 nm number peak and a 100 nm peak.  Similary for OC size split.
!                                                                              .
!  This function is adapted from emisOCbond.f and emisBCbond.f by Jeff Pierce
!  (Jan, 2007) used in GISS GCM-II'.  Introduced to GEOS-Chem by Win T.(9/4/07)
!
!
! !REVISION HISTORY:
!  04 Sep 2007 - W. Trivitayanurak - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER               :: I, J, K
   REAL(fp)              :: scalefactor(IBINS)

   scalefactor = 0.0d0

   SELECT CASE (CTYPE)
   CASE (1)
      SELECT CASE (STYPE)
      CASE (1)
         scalefactor(:) = ECSCALE30(:)
      CASE (2)
         scalefactor(:) = ECSCALE100(:)
      CASE (3)
         scalefactor(:) = ECSCALE100(:)
      END SELECT
   CASE (2)
      SELECT CASE (STYPE)
      CASE (1)
         scalefactor(:) = OCSCALE30(:)
      CASE (2)
         scalefactor(:) = OCSCALE100(:)
      CASE (3)
         scalefactor(:) = OCSCALE100(:)
      END SELECT
   END SELECT

   DO I = 1, State_Grid%NX
   DO J = 1, State_Grid%NY
      tomasdist(I,J,:) = BULKEMIS(I,J)*scalefactor(:)
   ENDDO
   ENDDO

 END FUNCTION SCALECARB
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emisscarbon
!
! !DESCRIPTION: Subroutine EMISSCARBON is the emissions routine for the carbon
! module. All carbon emissions, incl. SESQ and SVOC, are calculated through
! HEMCO and this module simply makes sure that the SESQ and SVOC emissions (if
! defined) are properly passed to the internal arrays.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE EMISSCARBON( Input_Opt, State_Grid, State_Met, RC )
!
! !USES:
!
   USE ErrCode_Mod
   USE HCO_INTERFACE_MOD,     ONLY : HcoState, GetHcoID, GetHcoVal
   USE HCO_ERROR_MOD
   USE Input_Opt_Mod,         ONLY : OptInput
   USE State_Grid_Mod,        ONLY : GrdState
   USE State_Met_Mod,         ONLY : MetState
   USE TIME_MOD,              ONLY : GET_TS_EMIS
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput),  INTENT(IN   )  :: Input_Opt   ! Input Options object
   TYPE(GrdState),  INTENT(IN   )  :: State_Grid  ! Grid State object
   TYPE(MetState),  INTENT(IN   )  :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
   INTEGER,         INTENT(INOUT)  :: RC          ! Failure?
!
! !REMARKS:
! SVOC emissions are expected to be fully calculated by HEMCO, i.e. its
! emissions need be specified in the HEMCO configuration file. In the original
! code, the emissions were calculated by scaling anthropogenic and
! biomass burning OC emissions. The same behavior can be achieved in HEMCO
! by assigning the desired SVOC species name to the given source type, e.g.:
!
! 0 BOND\_ANTH\_POG1 Bond\_fossil.nc OC 2000/1-12/1/0 C xy kg/m2/s POG1 74 1 1
!
! All POG1 emissions (anthropogenic + biomass burning) will go into
! POAEMISS(:,:,:,1) and all POG2 emissions will go into POAEMISS(:,:,:,2). SVOC
! emissions are assigned to POG1 and POG2 in HEMCO using a ratio of 0.49:0.51.
! We no longer separate anthropogenic from biomass burning since
! this appears to have been done only for debugging purposes. Routine CHEM_NVOC
! handles passing POAEMISS to the two gas-phase semivolatile species in the
! GM0 array.
!
! IMPORTANT: The SVOC emissions scale factor should be applied through HEMCO.
! In the example above, scale factor 74 represents the scale factor POGSCAL.
! The SCALING_POG1 scale factor is applied to the GFED biomass burning
! extensions. The two scale factors should be set to the same value in the
! HEMCO configuration file. The recommended value is 1.27.
!
! !REVISION HISTORY:
!  11 Nov 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER                  :: I, J, L, PBL_MAX
   INTEGER                  :: HCOPOG1, HCOPOG2
   REAL(fp)                 :: EMIS, TMPFLX
   REAL(fp)                 :: F_OF_PBL
   LOGICAL                  :: FOUND
   INTEGER, SAVE            :: SESQID = -999

   !=================================================================
   ! EMISSCARBON begins here!
   !=================================================================

   ! Assume success
   RC = GC_SUCCESS

   ! Initialize
   POAEMISS =  0e+0_fp

   ! Check if using complex SOA scheme
   IF ( Input_Opt%LSOA ) THEN

      ! Get HEMCO ID of species SESQ
      IF ( SESQID == -999 ) THEN
         SESQID = GetHcoID( 'SESQ' )
      ENDIF
      IF ( SESQID > 0 ) THEN
         IF ( .NOT. ASSOCIATED(HcoState%Spc(SESQID)%Emis%Val) ) THEN
            SESQID = -1
         ENDIF
      ENDIF

      ! Get HEMCO ID of species POG1 and POG2
      HCOPOG1 = GetHcoID( SpcID=id_POG1 )
      IF ( HCOPOG1 > 0 ) THEN
         IF ( .NOT. ASSOCIATED(HcoState%Spc(HCOPOG1)%Emis%Val) ) THEN
            HCOPOG1 = -1
         ENDIF
      ENDIF
      HCOPOG2 = GetHcoID( SpcID=id_POG2 )
      IF ( HCOPOG2 > 0 ) THEN
         IF ( .NOT. ASSOCIATED(HcoState%Spc(HCOPOG2)%Emis%Val) ) THEN
            HCOPOG2 = -1
         ENDIF
      ENDIF

   ELSE

      ! Do not get emissions of SESQ, POG1, POG2 if complex SOA is off
      SESQID  = -1
      HCOPOG1 = -1
      HCOPOG2 = -1

   ENDIF

   ! Nothing to do if none of the species are defined
   IF ( SESQID <= 0 .AND. HCOPOG1 <= 0 .AND. HCOPOG2 <=0 ) RETURN

   ! Maximum extent of PBL [model levels]
   PBL_MAX = State_Met%PBL_MAX_L

   !$OMP PARALLEL DO       &
   !$OMP DEFAULT( SHARED ) &
   !$OMP PRIVATE( I, J, L, F_OF_PBL, TMPFLX, Emis, FOUND )
   DO L = 1, PBL_MAX
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      ! Fraction of PBL spanned by grid box (I,J,L) [unitless]
      F_OF_PBL = State_Met%F_OF_PBL(I,J,L)

      ! Add sesquiterpene emissions from HEMCO to ORVC_SESQ array.
      ! We assume all SESQ emissions are placed in surface level.
      IF ( SESQID > 0 ) THEN
         CALL GetHcoVal( SESQID, I, J, 1, FOUND, Emis=EMIS )
         IF ( FOUND ) THEN
            ! Units from HEMCO are kgC/m2/s. Convert to kgC/box here.
            TMPFLX           = Emis * GET_TS_EMIS() * State_Grid%Area_M2(I,J)
            ORVC_SESQ(I,J,L) = ORVC_SESQ(I,J,L) + ( F_OF_PBL  * TMPFLX )
         ENDIF
      ENDIF

      ! Add SVOC emissions from HEMCO to POAEMISS array.
      ! Mix entire column emissions evenly in the PBL.
      !
      ! All SVOC emissions are now assigned to the POG1 and POG2 species in
      ! HEMCO to reflect that these emissions are added to the gas-phase
      ! species. The assignment of SVOC emissions to the two gas-phase
      ! species is actually performed in routine CHEM_NVOC. We also no
      ! longer separate anthropogenic from BF and BB emissions because
      ! this appears to have been done only for debugging purposes.
      ! (mps, 1/14/16)
      IF ( HCOPOG1 > 0 ) THEN
         ! Units from HEMCO are kgC/m2/s. Convert to kgC/box here.
         TMPFLX = SUM(HcoState%Spc(HCOPOG1)%Emis%Val(I,J,:)) &
                  * GET_TS_EMIS() * State_Grid%Area_M2(I,J)
         POAEMISS(I,J,L,1) = F_OF_PBL * TMPFLX
      ENDIF
      IF ( HCOPOG2 > 0 ) THEN
         ! Units from HEMCO are kgC/m2/s. Convert to kgC/box here.
         TMPFLX = SUM(HcoState%Spc(HCOPOG2)%Emis%Val(I,J,:)) &
                  * GET_TS_EMIS() * State_Grid%Area_M2(I,J)
         POAEMISS(I,J,L,2) = F_OF_PBL * TMPFLX
      ENDIF
   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Reset SVOC emissions to zero to make sure that they are
   ! not double-counted (when doing PBL mixing)
   IF ( HCOPOG1>0) HcoState%Spc(HCOPOG1)%Emis%Val = 0.0d0

 END SUBROUTINE EMISSCARBON
!EOC
#ifdef TOMAS
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emisscarbontomas
!
! !DESCRIPTION: Subroutine emisscarbontomas scales BULK HEMCO emissions into
! TOMAS arrays. Only use for TOMAS simulations. This is essential a re-write of
! the TOMAS portions of the v9 emisscarbon (JKodros 6/2/15)
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE EMISSCARBONTOMAS( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
   USE ErrCode_Mod
   USE ERROR_MOD
   USE Input_Opt_Mod,      ONLY : OptInput
   USE State_Chm_Mod,      ONLY : ChmState
   USE State_Grid_Mod,     ONLY : GrdState
   USE State_Met_Mod,      ONLY : MetState
   USE UnitConv_Mod,       ONLY : Convert_Spc_Units
   USE PRESSURE_MOD,       ONLY : GET_PCENTER
   USE TOMAS_MOD,          ONLY : IBINS,     AVGMASS, SOACOND
   USE TOMAS_MOD,          ONLY : ICOMP,     IDIAG
   USE TOMAS_MOD,          ONLY : CHECKMN
   USE HCO_INTERFACE_MOD,  ONLY : HcoState, GetHcoDiagn
   USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr !(ramnarine 12/27/2018)
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
   TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
   TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  12 Jun 2015 - J. Kodros   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Scalars
   INTEGER                  :: I, I0, IREF, J, J0, JREF, N
   INTEGER                  :: FLAG, ERR
   REAL(fp)                 :: DTSRCE, AREA_M2
   REAL(fp)                 :: XTRA_ORG_A(State_Grid%NX,State_Grid%NY)
   REAL(fp)                 :: CO_ANTH_TOTAL
   REAL(fp)                 :: BCSRC(State_Grid%NX,State_Grid%NY,IBINS,2)
   REAL(fp)                 :: OCSRC(State_Grid%NX,State_Grid%NY,IBINS,2)
   REAL(fp)                 :: NUMBSRC(State_Grid%NX,State_Grid%NY,IBINS)
   REAL(fp)                 :: AREA(State_Grid%NX, State_Grid%NY)
   REAL(fp)                 :: SIZE_DIST(State_Grid%NX,State_Grid%NY,IBINS,4) !(ramnarine 12/27/2018)
   REAL*4                   :: BOXVOL  ! calculated from State_Met
   REAL*4                   :: TEMPTMS ! calculated from State_Met
   REAL*4                   :: PRES    ! calculated from State_Met
   REAL(fp)                 :: TMP_MASS(State_Grid%NX,State_Grid%NY,IBINS)
   REAL(fp)                 :: OC2OM = 1.8d0
   LOGICAL                  :: SGCOAG = .True.
   INTEGER                  :: L, K, EMTYPE
   INTEGER                  :: ii=53, jj=29
   CHARACTER(LEN=63)        :: OrigUnit
   LOGICAL, SAVE            :: FIRST = .TRUE. !(ramnarine 12/27/2018)
   LOGICAL                  :: FND !(ramnarine 1/2/2019)
   LOGICAL                  :: prtDebug

   ! Strings
   CHARACTER(LEN= 63)       :: DgnName
   CHARACTER(LEN=255)       :: MSG
   CHARACTER(LEN=255)       :: LOC='EMISSCARBONTOMAS (carbon_mod.F90)'

   ! Pointers
   REAL(fp),        POINTER :: emis2D(:,:)
   REAL(f4),        POINTER :: Ptr2D(:,:)
   REAL(f4),        POINTER :: Ptr3D(:,:,:)

   !=================================================================
   ! EMISSCARBONTOMAS begins here!
   !=================================================================

   ! Assume success
   RC                  = GC_SUCCESS

   ! Print debug output?
   prtDebug            = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

   ! Get nested-grid offsets
   I0                  = State_Grid%XMinOffset
   J0                  = State_Grid%YMinOffset

   ! Initialize pointers
   emis2d              => NULL()
   Ptr2D               => NULL()
   Ptr3D               => NULL()

   ! Import emissions from HEMCO (through HEMCO state)
   IF ( .NOT. ASSOCIATED(HcoState) ) THEN
      CALL ERROR_STOP ( 'HcoState not defined!', LOC )
   ENDIF

   IF ( .NOT. (id_NK1 > 0   .AND. id_ECIL1 > 0 .AND. &
               id_ECOB1 > 0 .AND. id_OCIL1 > 0 .AND. &
               id_OCOB1 > 1 ) ) THEN
      CALL ERROR_STOP ( 'TOMAS Species not defined!', LOC )
   ENDIF

   ! Emission timestep [seconds]
   DTSRCE = HcoState%TS_EMIS

   ! Grid box aarea
   AREA = HcoState%Grid%AREA_M2%Val(:,:)

   ! Convert State_Chm%Species to [kg] for TOMAS. This will be
   ! removed once TOMAS uses mixing ratio instead of mass
   ! as species units (ewl, 9/11/15)
   CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                           'kg', RC, OrigUnit=OrigUnit )
   IF ( RC /= GC_SUCCESS ) THEN
      CALL GC_Error( 'Unit conversion error', RC, &
                     'Start of EMISSCARBONTOMAS in carbon_mod.F90' )
      RETURN
   ENDIF

   ! ---------FOSSIL FUEL EMISSIONS IN 3d---------------------------
   DO EMTYPE = 1,4
      SELECT CASE (EMTYPE)
      CASE (1)
         DgnName = 'BCPI_ANTH'
         emis2d  => BCPI_ANTH_BULK
      CASE (2)
         DgnName = 'OCPI_ANTH'
         emis2d  => OCPI_ANTH_BULK
      CASE (3)
         DgnName = 'BCPO_ANTH'
         emis2d  => BCPO_ANTH_BULK
      CASE (4)
         DgnName = 'OCPO_ANTH'
         emis2d  => OCPO_ANTH_BULK
      END SELECT

      CALL GetHcoDiagn( DgnName, .FALSE., ERR, Ptr3D=Ptr3D )
      IF ( .NOT. ASSOCIATED(Ptr3D) ) THEN
         CALL HCO_WARNING( 'Not found: '//TRIM(DgnName),ERR, THISLOC=LOC )
      ELSE
         emis2d(:,:) = 0.0d0

         !flatten emissions into a 2d array for now
         !they get distributed over the pbl by emithigh or subgridcoag
         DO L = 1, State_Grid%NZ
            emis2d(:,:) = emis2d(:,:) + Ptr3D(:,:,L)
         ENDDO
         ! [kg/box/time step]
         emis2d(:,:) = emis2d(:,:) * AREA(:,:) * DTSRCE
      ENDIF
      emis2d => NULL()
      Ptr3D  => NULL()
   ENDDO

   ! ---- SCALE INTO TOMAS BINS ---------------------------
   BCFF(:,:,:,1) = SCALECARB( State_Grid, BCPI_ANTH_BULK(:,:), 1, 1 )
   BCFF(:,:,:,2) = SCALECARB( State_Grid, BCPO_ANTH_BULK(:,:), 1, 1 )
   OCFF(:,:,:,1) = SCALECARB( State_Grid, OCPI_ANTH_BULK(:,:), 1, 2 ) * OC2OM
   OCFF(:,:,:,2) = SCALECARB( State_Grid, OCPO_ANTH_BULK(:,:), 1, 2 ) * OC2OM

   !end 3d emis

   DgnName = 'BCPI_BB'
   CALL GetHcoDiagn( DgnName, .FALSE., ERR, Ptr2D=Ptr2D )
   IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
      CALL HCO_WARNING('Not found: '//TRIM(DgnName),ERR,THISLOC=LOC)
   ELSE
      BCPI_BIOB_BULK = Ptr2D(:,:)
   ENDIF
   Ptr2D => NULL()

   DgnName = 'BCPO_BB'
   CALL GetHcoDiagn( DgnName, .FALSE., ERR, Ptr2D=Ptr2D )
   IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
      CALL HCO_WARNING('Not found: '//TRIM(DgnName),ERR,THISLOC=LOC)
   ELSE
      BCPO_BIOB_BULK = Ptr2D(:,:)
   ENDIF
   Ptr2D => NULL()

   DgnName = 'OCPI_BB'
   CALL GetHcoDiagn( DgnName, .FALSE., ERR, Ptr2D=Ptr2D )
   IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
      CALL HCO_WARNING('Not found: '//TRIM(DgnName),ERR,THISLOC=LOC)
   ELSE
      OCPI_BIOB_BULK = Ptr2D(:,:)
   ENDIF
   Ptr2D => NULL()

   DgnName = 'OCPO_BB'
   CALL GetHcoDiagn( DgnName, .FALSE., ERR, Ptr2D=Ptr2D )
   IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
      CALL HCO_WARNING('Not found: '//TRIM(DgnName),ERR,THISLOC=LOC)
   ELSE
      OCPO_BIOB_BULK = Ptr2D(:,:)
   ENDIF
   Ptr2D => NULL()

   ! ---- Fire Number -------------------------------------
   ! (ramnarine 12/27/2018)
   IF ( FIRST ) THEN
      ! Get a pointer to the fire number field if in use
      CALL HCO_GetPtr( HcoState, 'FINN_DAILY_NUMBER', FIRE_NUM, RC, FOUND=FND )
      !if biomass burning subgrid coagulation is in use,
      !FIRE_NUM will not be NULL() and therefore will not be ASSOCIATED()

      !reset first time flag
      FIRST = .FALSE.
   ENDIF

   IF ( ASSOCIATED(FIRE_NUM) ) THEN
      ! ---- Calling BB subgrid coag parameterization --------
      ! (ramnarine 12/27/2018)
      SIZE_DIST = SAKAMOTO_SIZE( State_Grid, State_Met, FIRE_NUM, &
                                 OCPI_BIOB_BULK, BCPI_BIOB_BULK, &
                                 OCPO_BIOB_BULK, BCPO_BIOB_BULK, &
                                 AREA )
   ENDIF

   !2d emis
   BCPI_BIOB_BULK = BCPI_BIOB_BULK(:,:) * AREA(:,:) * DTSRCE
   BCPO_BIOB_BULK = BCPO_BIOB_BULK(:,:) * AREA(:,:) * DTSRCE
   OCPI_BIOB_BULK = OCPI_BIOB_BULK(:,:) * AREA(:,:) * DTSRCE
   OCPO_BIOB_BULK = OCPO_BIOB_BULK(:,:) * AREA(:,:) * DTSRCE

   !2d bioburn
   IF ( ASSOCIATED(FIRE_NUM) ) THEN
      DO K = 1, IBINS        !ramnarine 12/27/2018
         BCBB(:,:,K,1) = SIZE_DIST(:,:,K,1) * AREA(:, :) * DTSRCE
         BCBB(:,:,K,2) = SIZE_DIST(:,:,K,2) * AREA(:, :) * DTSRCE
         OCBB(:,:,K,1) = SIZE_DIST(:,:,K,3) * AREA(:, :) * DTSRCE * OC2OM
         OCBB(:,:,K,2) = SIZE_DIST(:,:,K,4) * AREA(:, :) * DTSRCE * OC2OM
      ENDDO
   ELSE
      BCBB(:,:,:,1) = SCALECARB(State_Grid, BCPI_BIOB_BULK(:,:), 3,1)
      BCBB(:,:,:,2) = SCALECARB(State_Grid, BCPO_BIOB_BULK(:,:), 3,1)
      OCBB(:,:,:,1) = SCALECARB(State_Grid, OCPI_BIOB_BULK(:,:), 3,2) * OC2OM
      OCBB(:,:,:,2) = SCALECARB(State_Grid, OCPO_BIOB_BULK(:,:), 3,2) * OC2OM
   ENDIF

   ! Add into BCSRC and OCSRC
   BCSRC(:,:,:,1) = BCFF(:,:,:,1) + BCBF(:,:,:,1) + BCBB(:,:,:,1)
   BCSRC(:,:,:,2) = BCFF(:,:,:,2) + BCBF(:,:,:,2) + BCBB(:,:,:,2)
   OCSRC(:,:,:,1) = OCFF(:,:,:,1) + OCBF(:,:,:,1) + OCBB(:,:,:,1)
   OCSRC(:,:,:,2) = OCFF(:,:,:,2) + OCBF(:,:,:,2) + OCBB(:,:,:,2)

   IF ( SGCOAG ) THEN
      !emitsgc uses total OC or EC mass
      ! SUM mass terms into TEMP_MASS and pass to EMITSFC
      TMP_MASS = BCSRC(:,:,:,1) + BCSRC(:,:,:,2)
      CALL EMITSGC( Input_Opt, State_Chm, State_Grid, State_Met, TMP_MASS,  1 )

      TMP_MASS = OCSRC(:,:,:,1) + OCSRC(:,:,:,2)
      CALL EMITSGC( Input_Opt, State_Chm, State_Grid, State_Met, TMP_MASS,  2 )
   ELSE
      !-----------------------------------------
      ! Add emission w/o sub-grid coagulation
      !-----------------------------------------

      ! Convert the total mass emission to number emisison [No.]
      DO K = 1, IBINS
         NUMBSRC(:,:,K) = ( BCSRC(:,:,K,1) + BCSRC(:,:,K,2) + &
                            OCSRC(:,:,K,1) + OCSRC(:,:,K,2) )/ AVGMASS(K)
      ENDDO

      CALL EMITHIGH2( Input_Opt, State_Chm, State_Grid, State_Met, &
                      BCSRC,     OCSRC,     NUMBSRC )
   ENDIF  !sgcoag

   !end anthro emissions

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   !%%%% NOTE: BIOGENIC_SOAS is defined in hcoi_gc_diagn_mod.F90,   %%%%
   !%%%% which is only called if BPCH_DIAG=y.(bmy, 8/7/18)          %%%%
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! READ IN directly emitted SOAS (sfarina / jkodros)
   Ptr2D => NULL()
   DgnName = 'BIOGENIC_SOAS'
   CALL GetHcoDiagn( DgnName, .FALSE., RC, Ptr2D=Ptr2D )
   IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
      CALL GC_Error('Not found: '//TRIM(DgnName), RC, THISLOC=LOC)
      RETURN
   ENDIF
   TERP_ORGC = Ptr2D(:,:)
   TERP_ORGC = TERP_ORGC(:,:) * AREA(:,:) * DTSRCE
   Ptr2D => NULL()

   !$OMP PARALLEL DO       &
   !$OMP DEFAULT( SHARED ) &
   !$OMP PRIVATE( I, J, BOXVOL, TEMPTMS, PRES )
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX
      CALL CHECKMN( I, J, 1, Input_Opt, State_Chm, State_Grid, &
                    State_Met, 'CHECKMN from emisscarbontomas', RC)
      IF ( TERP_ORGC(I,J) > 0.d0 ) THEN
         BOXVOL  = State_Met%AIRVOL(I,J,1) * 1.e6 !convert from m3 -> cm3
         TEMPTMS = State_Met%T(I,J,1)
         PRES    = GET_PCENTER(I,J,1)*100.0 ! in Pa
         CALL SOACOND( TERP_ORGC(I,J), I, J, 1, BOXVOL, TEMPTMS, PRES, &
                       State_Chm, State_Grid, RC )
      END IF
   END DO
   END DO
   !$OMP END PARALLEL DO

   IF ( prtDebug ) CALL DEBUG_MSG( '### EMISCARB: after SOACOND (BIOG) ' )

   ! Convert State_Chm%Species back to original units (ewl, 9/11/15)
   CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                           OrigUnit, RC )
   IF ( RC /= GC_SUCCESS ) THEN
      CALL GC_Error('Unit conversion error', RC, &
                    'End of EMISSCARBONTOMAS in carbon_mod.F90')
      RETURN
   ENDIF

 END SUBROUTINE EMISSCARBONTOMAS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emithigh2
!
! !DESCRIPTION: Subroutine EMITHIGH2 mixes species completely from the surface
!  to the PBL top. This is a copy of subroutine EMITHIGH modified to work with
!  30-bin EC and OC mass and also aerosol number.  (win, 9/4/07)
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE EMITHIGH2( Input_Opt, State_Chm, State_Grid, State_Met, &
                       BCSRC,     OCSRC,     NUMBSRC )
!
! !USES:
!
   USE Input_Opt_Mod,      ONLY : OptInput
   USE State_Chm_Mod,      ONLY : ChmState
   USE State_Grid_Mod,     ONLY : GrdState
   USE State_Met_Mod,      ONLY : MetState
   USE TOMAS_MOD,          ONLY : IBINS
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)  :: Input_Opt               ! Input Options
   TYPE(GrdState), INTENT(IN)  :: State_Grid              ! Grid State
   TYPE(MetState), INTENT(IN)  :: State_Met               ! Meteorology State
   REAL(fp),       INTENT(IN)  :: BCSRC(State_Grid%NX,State_Grid%NY,IBINS,2) ! Total BC [kg]
   REAL(fp),       INTENT(IN)  :: OCSRC(State_Grid%NX,State_Grid%NY,IBINS,2) ! Total OC [kg]
   REAL(fp),       INTENT(IN)  :: NUMBSRC(State_Grid%NX,State_Grid%NY,IBINS)
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  04 Sep 2007 - W. Trivitayanurak - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER           :: I, J, L, PBL_MAX, K
   REAL(fp)          :: F_OF_PBL

   ! Pointers
   REAL(fp), POINTER :: Spc(:,:,:,:)

   !=================================================================
   ! EMITHIGH2 begins here!
   !=================================================================

   ! Point to chemical species array [kg]
   Spc => State_Chm%Species

   ! Maximum extent of PBL [model levels]
   PBL_MAX = State_Met%PBL_MAX_L

   !=================================================================
   ! Partition emissions throughout the boundary layer
   !=================================================================
   IF ( Input_Opt%LNLPBL ) THEN
      print *,'Currently subroutine EMITHIGH2 does not support ', &
              'the new non-local PBL scheme!'
      stop
   ENDIF

   !$OMP PARALLEL DO       &
   !$OMP DEFAULT( SHARED ) &
   !$OMP PRIVATE( I, J, L, K, F_OF_PBL )
   DO L = 1, PBL_MAX
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX
   DO K = 1, IBINS

      ! Fraction of PBL spanned by grid box (I,J,L) [unitless]
      F_OF_PBL = State_Met%F_OF_PBL(I,J,L)

      ! Hydrophilic ELEMENTAL CARBON
      Spc(I,J,L,id_ECIL1-1+K) = Spc(I,J,L,id_ECIL1-1+K) + &
                                ( F_OF_PBL * BCSRC(I,J,K,1) )

      ! Hydrophobic ELEMENTAL CARBON
      Spc(I,J,L,id_ECOB1-1+K) = Spc(I,J,L,id_ECOB1-1+K) + &
                                ( F_OF_PBL * BCSRC(I,J,K,2) )

      ! Hydrophilic ORGANIC CARBON
      Spc(I,J,L,id_OCIL1-1+K) = Spc(I,J,L,id_OCIL1-1+K) + &
                                ( F_OF_PBL * OCSRC(I,J,K,1) )

      ! Hydrophobic ORGANIC CARBON
      Spc(I,J,L,id_OCOB1-1+K) = Spc(I,J,L,id_OCOB1-1+K) + &
                                ( F_OF_PBL * OCSRC(I,J,K,2) )

      ! Number corresponding to EC + OC [No.]
      Spc(I,J,L,id_NK1-1+K)   = Spc(I,J,L,id_NK1-1+K) + &
                                ( F_OF_PBL * NUMBSRC(I,J,K) )

   ENDDO
   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Free pointer
   NULLIFY( Spc )

 END SUBROUTINE EMITHIGH2
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ohno3time
!
! !DESCRIPTION: Subroutine OHNO3TIME computes the sum of cosine of the solar
!  zenith angle over a 24 hour day, as well as the total length of daylight.
!  This is needed to scale the offline OH and NO3 concentrations.
!  (rjp, bmy, 12/16/02, 1/18/05)
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE OHNO3TIME( State_Grid )
!
! !USES:
!
   USE State_Grid_Mod, ONLY : GrdState
   USE TIME_MOD,       ONLY : GET_NHMSb,   GET_ELAPSED_SEC
   USE TIME_MOD,       ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT
!
! !INPUT PARAMETERS:
!
   TYPE(GrdState), INTENT(IN) :: State_Grid
!
! !REVISION HISTORY:
!  16 Dec 2002 - R. Park - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   LOGICAL, SAVE  :: FIRST = .TRUE.
   INTEGER        :: I, J, L, N, NT, NDYSTEP
   REAL(fp)       :: A0, A1, A2, A3, B1, B2, B3
   REAL(fp)       :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
   REAL(fp)       :: SUNTMP(State_Grid%NX,State_Grid%NY)

   !=================================================================
   ! OHNO3TIME begins here!
   !=================================================================

   !  Solar declination angle (low precision formula, good enough for us):
   A0 = 0.006918
   A1 = 0.399912
   A2 = 0.006758
   A3 = 0.002697
   B1 = 0.070257
   B2 = 0.000907
   B3 = 0.000148
   R  = 2.* PI * float( GET_DAY_OF_YEAR() - 1 ) / 365.

   DEC = A0 - A1*cos(  R) + B1*sin(  R)  &
             - A2*cos(2*R) + B2*sin(2*R) &
             - A3*cos(3*R) + B3*sin(3*R)

   LHR0 = int(float( GET_NHMSb() )/10000.)

   ! Only do the following at the start of a new day
   IF ( FIRST .or. GET_GMT() < 1e-5 ) THEN

      ! Zero arrays
      TCOSZ(:,:) = 0e+0_fp

      ! NDYSTEP is # of chemistry time steps in this day
      NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 3600 / GET_TS_CHEM()

      ! NT is the elapsed time [s] since the beginning of the run
      NT = GET_ELAPSED_SEC()

      ! Loop forward through NDYSTEP "fake" timesteps for this day
      DO N = 1, NDYSTEP

         ! Zero SUNTMP array
         SUNTMP = 0e+0_fp

         ! Loop over surface grid boxes
         !$OMP PARALLEL DO       &
         !$OMP DEFAULT( SHARED ) &
         !$OMP PRIVATE( I, J, YMID_R, TIMLOC, AHR )
         DO J = 1, State_Grid%NY
         DO I = 1, State_Grid%NX

            ! Grid box latitude center [radians]
            YMID_R = State_Grid%YMid_R(I,J)

            TIMLOC = real(LHR0) + real(NT)/3600.0 + State_Grid%XMid(I,J) / 15.0

            DO WHILE (TIMLOC .lt. 0)
               TIMLOC = TIMLOC + 24.0
            ENDDO

            DO WHILE (TIMLOC .gt. 24.0)
               TIMLOC = TIMLOC - 24.0
            ENDDO

            AHR = abs(TIMLOC - 12.) * 15.0 * PI_180

            !===========================================================
            ! The cosine of the solar zenith angle (SZA) is given by:
            !
            !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR)
            !
            ! where LAT = the latitude angle,
            !       DEC = the solar declination angle,
            !       AHR = the hour angle, all in radians.
            !
            ! If SUNCOS < 0, then the sun is below the horizon, and
            ! therefore does not contribute to any solar heating.
            !===========================================================

            ! Compute Cos(SZA)
            SUNTMP(I,J) = sin(YMID_R) * sin(DEC) + &
                          cos(YMID_R) * cos(DEC) * cos(AHR)

            ! TCOSZ is the sum of SUNTMP at location (I,J)
            ! Do not include negative values of SUNTMP
            TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP(I,J), 0e+0_fp )

         ENDDO
         ENDDO
         !$OMP END PARALLEL DO

         ! Increment elapsed time [sec]
         NT = NT + GET_TS_CHEM()
      ENDDO

      ! Reset first-time flag
      FIRST = .FALSE.
   ENDIF

 END SUBROUTINE OHNO3TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_oh
!
! !DESCRIPTION: Function GET\_OH returns OH from State\_Chm%Species (for
!  coupled runs) or monthly mean OH (for offline runs).  Imposes a diurnal
!  variation on OH for offline simulations. (bmy, 7/9/04)
!\\
!\\
! !INTERFACE:
!
 FUNCTION GET_OH( I, J, L, Input_Opt, State_Chm, State_Met ) &
      RESULT( OH_MOLEC_CM3 )
!
! !USES:
!
   USE ERROR_MOD,          ONLY : ERROR_STOP
   USE Input_Opt_Mod,      ONLY : OptInput
   USE State_Chm_Mod,      ONLY : ChmState
   USE State_Met_Mod,      ONLY : MetState
   USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
   INTEGER,        INTENT(IN)  :: I           ! Longitude index
   INTEGER,        INTENT(IN)  :: J           ! Latitude index
   INTEGER,        INTENT(IN)  :: L           ! Altitude index
   TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
   TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
   TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !RETURN VALUE:
!
   REAL(fp)                    :: OH_MOLEC_CM3
!
! !REVISION HISTORY:
!  09 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   REAL(fp) :: MolecRatio ! moles C / moles species
   REAL(fp) :: OH_MW_kg   ! kg OH / mol

   !=================================================================
   ! GET_OH begins here!
   !=================================================================
   IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

      !---------------------
      ! Coupled simulation
      !---------------------

      ! OH is defined only in the chemistry grid
      IF ( State_Met%InChemGrid(I,J,L) ) THEN

         ! Get OH from State_Chm%Species [kg] and convert to [molec/cm3]
         MolecRatio = State_Chm%SpcData(id_OH)%Info%MolecRatio
         OH_MW_kg   = State_Chm%SpcData(id_OH)%Info%emMW_g * 1.e-3_fp

         OH_MOLEC_CM3 = State_Chm%Species(I,J,L,id_OH) &
                        * ( AVO / OH_MW_kg )           &
                        / ( State_Met%AIRVOL(I,J,L)    &
                        * 1e+6_fp * MolecRatio )

      ELSE

         OH_MOLEC_CM3 = 0e+0_fp

      ENDIF

   ELSE IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

      !---------------------
      ! Offline simulation
      !---------------------

      ! Test for sunlight...
      IF ( State_Met%SUNCOS(I,J) > 0e+0_fp .and. TCOSZ(I,J) > 0e+0_fp ) THEN

         ! Impose a diurnal variation on OH during the day
         OH_MOLEC_CM3 = OH(I,J,L)                                * &
                        ( State_Met%SUNCOS(I,J) / TCOSZ(I,J) )   * &
                        ( 86400e+0_fp / GET_TS_CHEM() )

         ! OH is in kg/m3 (from HEMCO), convert to molec/cm3 (mps, 9/18/14)
         OH_MOLEC_CM3 = OH_MOLEC_CM3 * XNUMOL_OH / CM3PERM3

         ! Make sure OH is not negative
         OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0e+0_fp )

      ELSE

         ! At night, OH goes to zero
         OH_MOLEC_CM3 = 0e+0_fp

      ENDIF

   ELSE

      !---------------------
      ! Invalid sim type!
      !---------------------
      CALL ERROR_STOP( 'Invalid Simulation Type!', &
                       'GET_OH ("carbon_mod.F90")' )

   ENDIF

 END FUNCTION GET_OH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_no3
!
! !DESCRIPTION: Function GET\_NO3 returns NO3 from State\_Chm%Species (for
!  coupled runs) or monthly mean OH (for offline runs).  For offline runs, the
!  concentration of NO3 is set to zero during the day. (rjp, bmy, 12/16/02,
!  7/20/04)
!\\
!\\
! !INTERFACE:
!
 FUNCTION GET_NO3( I, J, L, Input_Opt, State_Chm, State_Met ) &
      RESULT( NO3_MOLEC_CM3 )
!
! !USES:
!
   USE ERROR_MOD,          ONLY : ERROR_STOP
   USE Input_Opt_Mod,      ONLY : OptInput
   USE State_Chm_Mod,      ONLY : ChmState
   USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
   INTEGER,        INTENT(IN)  :: I           ! Longitude index
   INTEGER,        INTENT(IN)  :: J           ! Latitude index
   INTEGER,        INTENT(IN)  :: L           ! Altitude index
   TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
   TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
   TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !RETURN VALUE:
!
   REAL(fp)                    :: NO3_MOLEC_CM3
!
! !REVISION HISTORY:
!  16 Dec 2002 - R. Park - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
   REAL(fp),  PARAMETER :: XNUMOL_NO3 = AVO / 62e-3_fp ! hard-coded MW
!
! !LOCAL VARIABLES:
!
   REAL(fp)             :: MolecRatio ! moles C / moles species
   REAL(fp)             :: NO3_MW_kg  ! kg NO3 / mol
   REAL(fp)             :: BOXVL

   !=================================================================
   ! GET_NO3 begins here!
   !=================================================================
   IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

      !----------------------
      ! Fullchem simulation
      !----------------------

      ! NO3 is defined only in the chemistry grid
      IF ( State_Met%InChemGrid(I,J,L) ) THEN

         ! Get NO3 from State_Chm%Species [kg] and convert to [molec/cm3]
         MolecRatio  = State_Chm%SpcData(id_NO3)%Info%MolecRatio
         NO3_MW_kg   = State_Chm%SpcData(id_NO3)%Info%emMW_g*1.e-3_fp

         NO3_MOLEC_CM3 = State_Chm%Species(I,J,L,id_NO3) &
                         * ( AVO / NO3_MW_kg ) &
                         / ( State_Met%AIRVOL(I,J,L) &
                         * 1e+6_fp * MolecRatio )

      ELSE

         NO3_MOLEC_CM3 = 0e+0_fp

      ENDIF

   ELSE IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

      !==============================================================
      ! Offline simulation: Read monthly mean GEOS-CHEM NO3 fields
      ! in [v/v].  Convert these to [molec/cm3] as follows:
      !
      !  vol NO3   moles NO3    kg air     kg NO3/mole NO3
      !  ------- = --------- * -------- * ---------------- =  kg NO3
      !  vol air   moles air      1        kg air/mole air
      !
      ! And then we convert [kg NO3] to [molec NO3/cm3] by:
      !
      !  kg NO3   molec NO3   mole NO3     1     molec NO3
      !  ------ * --------- * -------- * ----- = ---------
      !     1     mole NO3     kg NO3     cm3       cm3
      !          ^                    ^
      !          |____________________|
      !            this is XNUMOL_NO3
      !
      ! If at nighttime, use the monthly mean NO3 concentration from
      ! the NO3 array of "global_no3_mod.f".  If during the daytime,
      ! set the NO3 concentration to zero.  We don't have to relax to
      ! the monthly mean concentration every 3 hours (as for HNO3)
      ! since NO3 has a very short lifetime. (rjp, bmy, 12/16/02)
      !==============================================================

      ! Test if daylight
      IF ( State_Met%SUNCOS(I,J) > 0e+0_fp ) THEN

         ! NO3 goes to zero during the day
         NO3_MOLEC_CM3 = 0e+0_fp

      ELSE

         ! At night: Get NO3 [v/v] and convert it to [kg]
         NO3_MOLEC_CM3 = NO3(I,J,L) * State_Met%AD(I,J,L) * ( 62e+0_fp/AIRMW )

         ! Convert NO3 from [kg] to [molec/cm3]
         BOXVL         = State_Met%AIRVOL(I,J,L) * 1e+6_fp
         NO3_MOLEC_CM3 = NO3_MOLEC_CM3 * XNUMOL_NO3 / BOXVL

      ENDIF

      ! Make sure NO3 is not negative
      NO3_MOLEC_CM3  = MAX( NO3_MOLEC_CM3, 0e+0_fp )

   ELSE

      !----------------------
      ! Invalid sim type!
      !----------------------
      CALL ERROR_STOP( 'Invalid Simulation Type!', &
                       'GET_NO3 ("carbon_mod.F90")' )

   ENDIF

 END FUNCTION GET_NO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_o3
!
! !DESCRIPTION: Function GET\_O3 returns monthly mean O3 for offline sulfate
!  aerosol simulations. (bmy, 12/16/02, 7/20/04)
!\\
!\\
! !INTERFACE:
!
 FUNCTION GET_O3( I, J, L, Input_Opt, State_Chm, State_Grid, State_Met ) &
      RESULT( O3_MOLEC_CM3 )
!
! !USES:
!
   USE ERROR_MOD,          ONLY : ERROR_STOP
   USE Input_Opt_Mod,      ONLY : OptInput
   USE State_Chm_Mod,      ONLY : ChmState
   USE State_Grid_Mod,     ONLY : GrdState
   USE State_Met_Mod,      ONLY : MetState
   USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
   INTEGER,        INTENT(IN)  :: I           ! Longitude index
   INTEGER,        INTENT(IN)  :: J           ! Latitude index
   INTEGER,        INTENT(IN)  :: L           ! Altitude index
   TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
   TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
   TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
   TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !RETURN VALUE:
!
   REAL(fp)                    :: O3_MOLEC_CM3
!
! !REVISION HISTORY:
!  16 Dec 2002 - R. Park - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
   REAL(fp),  PARAMETER :: XNUMOL_O3 = AVO / 48e-3_fp ! hard-coded MW
!
! !LOCAL VARIABLES:
!
   REAL(fp)             :: MolecRatio ! moles C / moles species
   REAL(fp)             :: O3_MW_kg   ! kg O3 / mol
   REAL(fp)             :: BOXVL

   !=================================================================
   ! GET_O3 begins here!
   !=================================================================
   IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

      !--------------------
      ! Coupled simulation
      !--------------------

      ! O3 is defined only in the chemistry grid
      IF ( State_Met%InChemGrid(I,J,L) ) THEN

         ! Get O3 from State_Chm%Species [kg] and convert to [molec/cm3]
         MolecRatio = State_Chm%SpcData(id_O3)%Info%MolecRatio
         O3_MW_kg   = State_Chm%SpcData(id_O3)%Info%emMW_g * 1.e-3_fp

         O3_MOLEC_CM3 = State_Chm%Species(I,J,L,id_O3) &
                        * ( AVO / O3_MW_kg ) &
                        / ( State_Met%AIRVOL(I,J,L) &
                        * 1e+6_fp * MolecRatio )

      ELSE

         O3_MOLEC_CM3 = 0e+0_fp

      ENDIF

   ELSE IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

      !--------------------
      ! Offline simulation
      !--------------------

      ! Get O3 [v/v] for this gridbox & month
      ! O3 is defined only in the chemistry grid
      IF ( L <= State_Grid%MaxChemLev ) THEN

         ! Get O3 [v/v] and convert it to [kg]
         O3_MOLEC_CM3 = O3(I,J,L) * State_Met%AD(I,J,L) * &
                        ( 48e+0_fp / AIRMW )            ! hard-coded MW

         ! Convert O3 from [kg] to [molec/cm3]
         BOXVL        = State_Met%AIRVOL(I,J,L) * 1e+6_fp
         O3_MOLEC_CM3 = O3_MOLEC_CM3 * XNUMOL_O3 / BOXVL

      ELSE
         O3_MOLEC_CM3 = 0e+0_fp
      ENDIF

      ! Test for sunlight...
      IF ( State_Met%SUNCOS(I,J) > 0e+0_fp .and. TCOSZ(I,J) > 0e+0_fp ) THEN

         ! Impose a diurnal variation on OH during the day
         O3_MOLEC_CM3 = O3_MOLEC_CM3                             * &
                        ( State_Met%SUNCOS(I,J) / TCOSZ(I,J) )   * &
                        ( 86400e+0_fp / GET_TS_CHEM() )

         ! Make sure OH is not negative
         O3_MOLEC_CM3 = MAX( O3_MOLEC_CM3, 0e+0_fp )

      ELSE
         O3_MOLEC_CM3 = 0e+0_fp
      ENDIF

   ELSE

      !--------------------
      ! Invalid sim type!
      !--------------------
      CALL ERROR_STOP( 'Invalid Simulation Type!', &
                       'GET_O3 ("carbon_mod.F90")' )

   ENDIF

 END FUNCTION GET_O3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_daro2
!
! !DESCRIPTION: Function GET\_DARO2 returns the amount of aromatic peroxy
!  radical that reacted with HO2 or NO during the last chemistry timestep.
!  (dkh, 11/10/06)
!\\
!\\
! !INTERFACE:
!
 FUNCTION GET_DARO2( I, J, L, NOX, JHC, Input_Opt, State_Chm, State_Met ) &
      RESULT( DARO2 )
!
! !USES:
!
   USE ERROR_MOD,          ONLY : ERROR_STOP
   USE Input_Opt_Mod,      ONLY : OptInput
   USE PhysConstants,      ONLY : AVO
   USE State_Chm_Mod,      ONLY : ChmState
   USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
   INTEGER,        INTENT(IN)  :: I           ! Longitude index
   INTEGER,        INTENT(IN)  :: J           ! Latitude index
   INTEGER,        INTENT(IN)  :: L           ! Altitude index
   INTEGER,        INTENT(IN)  :: NOX
   INTEGER,        INTENT(IN)  :: JHC
   TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
   TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
   TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !RETURN VALUE:
!
   REAL(fp)                    :: DARO2
!
! !REVISION HISTORY:
!  10 Nov 2006 - D. Henze - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER  :: id_LARO2      ! Species ID for AROM lost to HO2 or NO
   INTEGER  :: id_AROM       ! Species ID of the parent aromatic
   REAL(fp) :: ARO2CARB      ! kg C of ARO2 / kg ARO2
   REAL(fp) :: AROM_MW_kg    ! g C of AROM  / mol C
   REAL(fp) :: LARO2_MW_kg   ! g C of LARO2 / mol C
   REAL(fp) :: MolecRatio    ! moles C / moles species

   !=================================================================
   ! GET_DARO2 begins here!
   !=================================================================

   ! Initialize
   DARO2 = 0e+0_fp

   IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

      !--------------------
      ! Coupled simulation
      !--------------------

      ! Test if we are in the chemistry grid
      IF ( State_Met%InChemGrid(I,J,L) ) THEN

         ! Get information on the
         ! specified type of aromatic peroxy radical
         ! Benzene
         IF    ( JHC == PARENTBENZ ) THEN

            ! Loss due to NO2 corresponds to high NOX experiments
            ! (NOX = 1) while loss due to HO2 corresponds to
            ! low NOX experiments (NOX = 2).
            IF ( NOX == 1 ) THEN
               id_LARO2 = id_LBRO2N
            ELSEIF ( NOX == 2 ) THEN
               id_LARO2 = id_LBRO2H
            ELSE
               CALL ERROR_STOP('Bad NOX', 'GET_DARO2')
            ENDIF

            ! kg C of ARO2 / kg ARO2
            ! dkh ARMv4 (hotp 7/31/2008)
            !ARO2CARB = 0.5669e+0_fp ! = 6*12/(6*12+3*16+7)
            !ARO2CARB = 0.4528e+0_fp ! = 6*12/(6*12+5*16+7)
            ! Now use parent HC instead of RO2 (hotp 5/13/10)
            ! C6H6
            ARO2CARB = 6e+0_fp * 12e+0_fp / ( 6e+0_fp * 12e+0_fp + 6e+0_fp )

            ! Species index of the parent aromatic
            id_AROM = id_BENZ

         ! Toluene
         ELSEIF ( JHC == PARENTTOLU ) THEN

            IF ( NOX == 1 ) THEN
               id_LARO2 = id_LTRO2N
            ELSEIF ( NOX == 2 ) THEN
               id_LARO2 = id_LTRO2H
            ELSE
               CALL ERROR_STOP('Bad NOX', 'GET_DARO2')
            ENDIF

            ! kg C of ARO2 / kg ARO2
            ! dkh ARMv4 (hotp 7/31/2008)
            !ARO2CARB = 0.5874 ! = 7*12/(7*12+3*16+11)  ! This was wrong for 2 reasons
            !ARO2CARB = 0.5957e+0_fp ! = 7*12/(7*12+3*16+9) ! <-- just change 11 to 9
            !ARO2CARB = 0.48e+0_fp ! = 7*12/(7*12+5*16+11)  ! <-- just change 3*16 to 5*16
            !ARO2CARB = 0.4855e+0_fp ! = 7*12/(7*12+5*16+9)  ! <-- change both
            ! Now use parent HC instead of RO2 (hotp 5/13/10)
            ! C7H8
            ARO2CARB = 7e+0_fp * 12e+0_fp / ( 7e+0_fp * 12e+0_fp + 8e+0_fp )

            ! Species index of the parent aromatic
            id_AROM = id_TOLU

         ! XYLENE
         ELSEIF ( JHC == PARENTXYLE ) THEN

            IF ( NOX == 1 ) THEN
               id_LARO2 = id_LXRO2N
            ELSEIF ( NOX == 2 ) THEN
               id_LARO2 = id_LXRO2H
            ELSE
               CALL ERROR_STOP('Bad NOX', 'GET_DARO2')
            ENDIF

            ! kg C of ARO2 / kg ARO2
            ! dkh ARMv4 (hotp 7/31/2008)
            !ARO2CARB = 0.6194e+0_fp ! = 8*12/(8*12+3*16+11)
            !ARO2CARB = 0.5134e+0_fp ! = 8*12/(8*12+3*16+11)
            ! comments on above are bad (hotp 7/22/09)
            ! ARO2CARB for XYL is = 8*12/(8*12+5*16+11) (hotp 7/22/09)
            ! Now use parent HC instead of RO2 (hotp 5/13/10)
            ! old value based on RO2: 0.5134e+0_fp
            ! C8H10
            ARO2CARB = 8e+0_fp * 12e+0_fp / ( 8e+0_fp * 12e+0_fp + 10e+0_fp )

            ! Species index of the parent aromatic
            id_AROM = id_XYLE

         ! NAPHTHALENE (IVOC surrogate) (hotp 7/22/09)
         ELSEIF ( JHC == PARENTNAP ) THEN

            IF ( NOX == 1 ) THEN
               id_LARO2 = id_LNRO2N
            ELSEIF ( NOX == 2 ) THEN
               id_LARO2 = id_LNRO2H
            ELSE
               CALL ERROR_STOP('Bad NOX', 'GET_DARO2')
            ENDIF

            ! kg C of ARO2 / kg ARO2
            ! NOTE: for NAP, GET_DARO2 is the kg of NAP reacted,
            ! not the kg of ARO2
            ! ALPHAs are set up for this
            ARO2CARB = 12e+0_fp * 10e+0_fp / ( 12e+0_fp * 10e+0_fp + 8e+0_fp )

            ! Species index of the parent aromatic
            id_AROM = id_NAP

         ELSE

            CALL ERROR_STOP('Bad JHC', 'GET_DAR2')

         ENDIF

         !-----------------------------------------------------------
         ! Get LARO2 (AROM lost to HO2 or NO) from State_Chm%Species
         ! [kg LARO2] and convert to [kg AROM]
         !
         ! We use MolecRatio for the parent aromatic hydrocarbon,
         ! AROM, because:
         !   atom  C / mol ARO2  = atom C / mol AROM
         !-----------------------------------------------------------

         ! Now get the species coefficient from the species database
         ! instead of from Input_Opt (bmy, 5/17/16)
         MolecRatio  = State_Chm%SpcData(id_AROM)%Info%MolecRatio
         AROM_MW_kg  = State_Chm%SpcData(id_AROM)%Info%emMW_g * 1.e-3_fp
         LARO2_MW_kg = State_Chm%SpcData(id_LARO2)%Info%emMW_g * 1.e-3_fp

         DARO2 = State_Chm%Species(I,J,L,id_LARO2) * ( AVO / LARO2_MW_kg ) &
                 / ( AVO / AROM_MW_kg  ) * MolecRatio * ARO2CARB

      ELSE

         ! Otherwise set DOH=0
         DARO2 = 0e+0_fp

      ENDIF

   ELSE IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

      !--------------------
      ! Offline simulation
      !--------------------

      ! Not supported yet for
      ! offline aerosol simulations, set DOH=0
      DARO2 = 0e+0_fp

   ELSE

      !--------------------
      ! Invalid sim type!
      !--------------------
      CALL ERROR_STOP( 'Invalid simulation type!', &
                       'GET_DARO2 ("carbon_mod.F90")' )

   ENDIF

 END FUNCTION GET_DARO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_doh
!
! !DESCRIPTION: Function GET\_DOH returns the amount of isoprene [kg] that has
!  reacted with OH during the last chemistry time step. (dkh, bmy, 6/01/06)
!\\
!\\
! !INTERFACE:
!
 FUNCTION GET_DOH( I, J, L, Input_Opt, State_Chm, State_Met ) &
      RESULT( DOH )
!
! !USES:
!
   USE ERROR_MOD,          ONLY : ERROR_STOP
   USE Input_Opt_Mod,      ONLY : OptInput
   USE PhysConstants,      ONLY : AVO
   USE State_Chm_Mod,      ONLY : ChmState
   USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
   INTEGER,        INTENT(IN)  :: I           ! Longitude index
   INTEGER,        INTENT(IN)  :: J           ! Latitude index
   INTEGER,        INTENT(IN)  :: L           ! Altitude index
   TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
   TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
   TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !RETURN VALUE:
!
   REAL(fp)                    :: DOH
!
! !REVISION HISTORY:
!  01 Jun 2006 - D. Henze - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   REAL(fp) :: MolecRatio    ! moles C / moles ISOP
   REAL(fp) :: ISOP_MW_kg    ! kg C ISOP    / mol C
   REAL(fp) :: LISOPOH_MW_kg ! kg C LISOPOH / mol C

   !=================================================================
   ! GET_DOH begins here!
   !=================================================================

   IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

      !--------------------
      ! Coupled simulation
      !--------------------

      ! Test if we are in the chemistry grid
      IF ( State_Met%InChemGrid(I,J,L) ) THEN

         !-----------------------------------------------------------
         ! Get LISOPOH (ISOP lost to OH) from State_Chm%Species
         ! [kg LISOPOH] and convert to [kg C ISOP]
         !-----------------------------------------------------------
         MolecRatio    = State_Chm%SpcData(id_ISOP)%Info%MolecRatio
         ISOP_MW_kg    = State_Chm%SpcData(id_ISOP)%Info%emMW_g * 1.e-3_fp
         LISOPOH_MW_kg = State_Chm%SpcData(id_LISOPOH)%Info%emMW_g * 1.e-3_fp

         DOH = State_Chm%Species(I,J,L,id_LISOPOH) * ( AVO / LISOPOH_MW_kg ) &
               / ( AVO / ISOP_MW_kg ) * MolecRatio

      ELSE

         ! Otherwise set DOH=0
         DOH = 0e+0_fp

      ENDIF

   ELSE IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

      !--------------------
      ! Offline simulation
      !--------------------

      ! ISOP from OH not is yet supported for
      ! offline aerosol simulations, set DOH=0
      DOH = 0e+0_fp

   ELSE

      !--------------------
      ! Invalid sim type!
      !--------------------
      CALL ERROR_STOP( 'Invalid simulation type!', &
                       'GET_DOH ("carbon_mod.F90")' )

   ENDIF

 END FUNCTION GET_DOH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_eqlb
!
! !DESCRIPTION: Subroutine CHECK\_EQLB makes sure aerosols are at equilibrium
!  (checks SOA=SOG*KOM*Mo). Called inside SOA\_SVPOA\_CHEMISTRY I, J, L loop
!  after SOA\_SVPOA\_LUMP. Created by Havala Pye (5/18/10).
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE CHECK_EQLB( I,   J,   L, KOMIJL, CONVFAC, MSOACHEM, &
                        LOW, TOL, ASOANGAS, ASOANAER, OCPIOCPO, State_Chm )
!
! !USES:
!
   USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
   INTEGER,        INTENT(IN)  :: I           ! Longitude index
   INTEGER,        INTENT(IN)  :: J           ! Latitude index
   INTEGER,        INTENT(IN)  :: L           ! Altitude index
   REAL(fp),       INTENT(IN)  :: KOMIJL(MPROD,MSV) ! KOM at grid box (adj T)
   REAL(fp),       INTENT(IN)  :: CONVFAC     ! Conversion factor kg to ug/m3
   REAL(fp),       INTENT(IN)  :: OCPIOCPO    ! POA mass [ug/m3]

   ! Arguments for debugging
   REAL(fp),       INTENT(IN)  :: MSOACHEM    ! MNEW from calling prog
   REAL(fp),       INTENT(IN)  :: LOW         ! Lower bound on soln
   REAL(fp),       INTENT(IN)  :: TOL         ! Tolerance on soln
   REAL(fp),       INTENT(IN)  :: ASOANGAS    ! Gas phase ASOAN (should =0)
   REAL(fp),       INTENT(IN)  :: ASOANAER    ! Aer phase ASOAN [ug/m3]

   TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
!
! !REMARKS:
!  Note: There are some deviations from equilibrium due to the fact
!  that ASOAN is supposed to be nonvolatile, but is modeled with a KOM of
!  10^6. An adjustment is made in SOA_CHEMISTRY to force all ASOAN to
!  the aerosol phase. This was found to lead to error up to 1e-5 ug/m3
!  in Mo. This error is small, but the effects can be investigated
!  here if you're interested!
!                                                                             .
!  As of 6/2010, KOM for ASOAN was increased and the error in Mo reduced.
!
! !REVISION HISTORY:
!  18 May 2010 - H.O.T. Pye - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
   !REAL(fp), PARAMETER   :: ACCEPTERRORUG = 1d-6 ! error in ug/m3
   ! KOM_REF for non-vol is larger, so error smaller (hotp 5/28/10)
   !REAL(fp), PARAMETER   :: ACCEPTERRORUG = 1d-10 ! error in ug/m3
   ! relax error tolerance (KOM_REF for ASOAN still not perfect)
   ! hotp 6/11/10
   REAL(fp), PARAMETER :: ACCEPTERRORUG = 1e-8_fp ! error in ug/m3
!
! !LOCAL VARIABLES:}
!
   INTEGER             :: NOX, IPR, JHC, JSV
   REAL(fp)            :: MOTEMP, OATEMP, EQLBDIFF

   ! Pointers
   REAL(fp), POINTER   :: Spc(:,:,:,:)

   !=================================================================
   ! CHECK_EQLB starts here
   !=================================================================

   ! Point to chemical species array [kg]
   Spc => State_Chm%Species

   ! Calculate mass of absorbing organic medium
   MOTEMP = Spc(I,J,L,id_ASOAN) + &
            Spc(I,J,L,id_ASOA1) + &
            Spc(I,J,L,id_ASOA2) + &
            Spc(I,J,L,id_ASOA3) + &
            Spc(I,J,L,id_TSOA1) + &
            Spc(I,J,L,id_TSOA2) + &
            Spc(I,J,L,id_TSOA3) + &
            Spc(I,J,L,id_TSOA0)

   ! Add primary material as appropriate
   IF ( id_POA1 > 0 ) THEN
      MOTEMP = MOTEMP              + &
               Spc(I,J,L,id_POA1 ) * OCFPOA(I,J)  + &
               Spc(I,J,L,id_POA2 ) * OCFPOA(I,J)  + &
               Spc(I,J,L,id_OPOA1) * OCFOPOA(I,J) + &
               Spc(I,J,L,id_OPOA2) * OCFOPOA(I,J)
   ELSEIF ( id_OCPI > 0 ) THEN
      MOTEMP = MOTEMP + &
               ( Spc(I,J,L,id_OCPI) + Spc(I,J,L,id_OCPO) ) * 2.1e+0_fp
   ENDIF

   ! Convert Mo from [kg] to [ug/m3]
   MOTEMP = MOTEMP * CONVFAC

   ! Check Mo calculation
   ! Forcing ASOAN to aerosol phase causes errors in MO that will
   ! manifest themselves here
   ! Require MO to be accurate within 1d-8 ug/m3
   EQLBDIFF = ABS( MOTEMP - MSOACHEM )
   !IF ( EQLBDIFF > 1d-4 ) print*, 'CHECK_EQLB ERROR: MO disagree',
   ! KOM_REF for non-vol is larger, so tighten here (hotp 5/28/10)
   IF ( EQLBDIFF > 1e-8_fp ) print*, 'CHECK_EQLB ERROR: MO disagree', &
                                     I,J,L,MSOACHEM,MOTEMP,LOW,TOL,   &
                                     ASOANGAS, ASOANAER, OCPIOCPO

   ! quick check
   !MOTEMP = MSOACHEM

   !----------------------------------------------------
   ! Semivolatile 1: monoterpene + sesquiterpene SOA
   !----------------------------------------------------
   JHC = PARENTMTPA
   JSV = IDSV(JHC)

   ! Product 1
   IPR = 1
   ! Compute OA in kg
   OATEMP = Spc(I,J,L,id_TSOG1) * KOMIJL(IPR,JSV) * MOTEMP
   ! Compute difference in ug/m3
   EQLBDIFF = ABS( OATEMP - Spc(I,J,L,id_TSOA1) )* CONVFAC
   ! Assess error
   IF ( EQLBDIFF > ACCEPTERRORUG ) THEN
      WRITE(*,*) 'EQLB Problem PR, JSV',  IPR, JSV, ' in box ', I, J, L
   ENDIF

   ! Product 2
   IPR = 2
   ! Compute OA in kg
   OATEMP = Spc(I,J,L,id_TSOG2) * KOMIJL(IPR,JSV) * MOTEMP
   ! Compute difference in ug/m3
   EQLBDIFF = ABS( OATEMP - Spc(I,J,L,id_TSOA2) )* CONVFAC
   ! Assess error
   IF ( EQLBDIFF > ACCEPTERRORUG ) THEN
      WRITE(*,*) 'EQLB Problem PR, JSV', IPR, JSV, ' in box ', I, J, L, &
                  MSOACHEM,MOTEMP,LOW,TOL, &
                  ASOANGAS, ASOANAER, OCPIOCPO, &
                  Spc(I,J,L,id_TSOA2),Spc(I,J,L,id_TSOG2)
      WRITE(*,*) 'KOM',KOMIJL(IPR,JSV),OATEMP,CONVFAC
   ENDIF

   ! Product 3
   IPR = 3
   ! Compute OA in kg
   OATEMP = Spc(I,J,L,id_TSOG3) * KOMIJL(IPR,JSV) * MOTEMP
   ! Compute difference in ug/m3
   EQLBDIFF = ABS( OATEMP - Spc(I,J,L,id_TSOA3) )* CONVFAC
   ! Assess error
   IF ( EQLBDIFF > ACCEPTERRORUG ) THEN
      WRITE(*,*) 'EQLB Problem PR, JSV', IPR, JSV, ' in box ', I, J, L
   ENDIF

   ! Product 4, C*=0.1
   IPR = 4
   ! Compute OA in kg
   OATEMP = Spc(I,J,L,id_TSOG0) * KOMIJL(IPR,JSV) * MOTEMP
   ! Compute difference in ug/m3
   EQLBDIFF = ABS( OATEMP - Spc(I,J,L,id_TSOA0) )* CONVFAC
   ! Assess error
   IF ( EQLBDIFF > ACCEPTERRORUG ) THEN
      WRITE(*,*) 'EQLB Problem PR, JSV', IPR, JSV, ' in box ', I, J, L, &
                  MSOACHEM,MOTEMP,LOW,TOL, &
                  ASOANGAS, ASOANAER, OCPIOCPO, &
                  Spc(I,J,L,id_TSOA0),Spc(I,J,L,id_TSOG0)
      WRITE(*,*) 'KOM',KOMIJL(IPR,JSV),OATEMP,CONVFAC
   ENDIF

   !---------------------------------------------------------------------------
   ! Prior to 7/15/19:
   ! Remove isoprene from VBS (mps, 7/15/19)
   !!----------------------------------------------------
   !! Semivolatile 2: isoprene SOA
   !!----------------------------------------------------
   !JHC = PARENTISOP
   !JSV = IDSV(JHC)
   !
   !! Product 1
   !IPR = 1
   !! Compute OA in kg
   !OATEMP = Spc(I,J,L,id_ISOG1) * KOMIJL(IPR,JSV) * MOTEMP
   !! Compute difference in ug/m3
   !EQLBDIFF = ABS( OATEMP - Spc(I,J,L,id_ISOA1) )* CONVFAC
   !! Assess error
   !IF ( EQLBDIFF > ACCEPTERRORUG ) THEN
   !   WRITE(*,*) 'EQLB Problem PR, JSV', IPR, JSV, ' in box ', I, J, L
   !ENDIF
   !
   !! Product 2
   !IPR = 2
   !! Compute OA in kg
   !OATEMP = Spc(I,J,L,id_ISOG2) * KOMIJL(IPR,JSV) * MOTEMP
   !! Compute difference in ug/m3
   !EQLBDIFF = ABS( OATEMP - Spc(I,J,L,id_ISOA2) )* CONVFAC
   !! Assess error
   !IF ( EQLBDIFF > ACCEPTERRORUG ) THEN
   !   WRITE(*,*) 'EQLB Problem PR, JSV', IPR, JSV,' in box ', I, J, L
   !ENDIF
   !
   !! Product 3
   !IPR = 3
   !! Compute OA in kg
   !OATEMP = Spc(I,J,L,id_ISOG3) * KOMIJL(IPR,JSV) * MOTEMP
   !! Compute difference in ug/m3
   !EQLBDIFF = ABS( OATEMP - Spc(I,J,L,id_ISOA3) )* CONVFAC
   !! Assess error
   !IF ( EQLBDIFF > ACCEPTERRORUG ) THEN
   !   WRITE(*,*) 'EQLB Problem PR, JSV', IPR, JSV, ' in box ', I, J, L
   !ENDIF
   !---------------------------------------------------------------------------

   !----------------------------------------------------
   ! Semivolatile 3: lumped arom/IVOC: total SOA+SOG in kg
   !----------------------------------------------------
   JHC = PARENTBENZ
   JSV = IDSV(JHC)

   ! High NOx, Product 1
   !NOX = NHIGHNOX
   IPR = 1
   ! Compute OA in kg
   OATEMP = Spc(I,J,L,id_ASOG1) * KOMIJL(IPR,JSV) * MOTEMP
   ! Compute difference in ug/m3
   EQLBDIFF = ABS( OATEMP - Spc(I,J,L,id_ASOA1) )* CONVFAC
   ! Assess error
   IF ( EQLBDIFF > ACCEPTERRORUG ) THEN
      WRITE(*,*) 'EQLB Problem PR, JSV', IPR, JSV, ' in box ', I, J, L
   ENDIF

   ! High NOx, Product 2
   !NOX = NHIGHNOX
   IPR = 2
   ! Compute OA in kg
   OATEMP = Spc(I,J,L,id_ASOG2) * KOMIJL(IPR,JSV) * MOTEMP
   ! Compute difference in ug/m3
   EQLBDIFF = ABS( OATEMP - Spc(I,J,L,id_ASOA2) )* CONVFAC
   ! Assess error
   IF ( EQLBDIFF > ACCEPTERRORUG ) THEN
      WRITE(*,*) 'EQLB Problem PR, JSV', IPR, JSV, ' in box ', I, J, L
   ENDIF

   ! High NOx, Product 3
   !NOX = NHIGHNOX
   IPR = 3
   ! Compute OA in kg
   OATEMP = Spc(I,J,L,id_ASOG3) * KOMIJL(IPR,JSV) * MOTEMP
   ! Compute difference in ug/m3
   EQLBDIFF = ABS( OATEMP - Spc(I,J,L,id_ASOA3) )* CONVFAC
   ! Assess error
   IF ( EQLBDIFF > ACCEPTERRORUG ) THEN
      WRITE(*,*) 'EQLB Problem PR, JSV', IPR, JSV, ' in box ', I, J, L
   ENDIF

   ! LOW NOx, Product 1
   ! nonvolatile so don't need to check partitioning

   !----------------------------------------------------
   ! POA: total POA+POG in kgC
   !----------------------------------------------------
   IF ( id_POA1 > 0 ) THEN
      JHC = PARENTPOA
      JSV = IDSV(JHC)
      NOX = NONLYNOX

      ! Product 1
      NOX = NONLYNOX
      IPR = 1
      ! Compute OA in kg
      OATEMP = Spc(I,J,L,id_POG1) * KOMIJL(IPR,JSV) * MOTEMP
      ! Compute difference in ug/m3
      EQLBDIFF = ABS( OATEMP - Spc(I,J,L,id_POA1) )* CONVFAC
      ! Assess error
      IF ( EQLBDIFF > ACCEPTERRORUG ) THEN
         WRITE(*,*) 'EQLB Problem NOX, PR, JSV', NOX, IPR, JSV, &
                    ' in box ', I, J, L
      ENDIF

      ! Product 2
      NOX = NONLYNOX
      IPR = 2
      ! Compute OA in kg
      OATEMP = Spc(I,J,L,id_POG2) * KOMIJL(IPR,JSV) * MOTEMP
      ! Compute difference in ug/m3
      EQLBDIFF = ABS( OATEMP - Spc(I,J,L,id_POA2) )* CONVFAC
      ! Assess error
      IF ( EQLBDIFF > ACCEPTERRORUG ) THEN
         WRITE(*,*) 'EQLB Problem NOX, PR, JSV', NOX, IPR, JSV, &
                    ' in box ', I, J, L
      ENDIF

   ENDIF ! POA

   !----------------------------------------------------
   ! OPOA: total SOA+SOG in kgC
   !----------------------------------------------------
   IF ( id_OPOA1 > 0 ) THEN
      JHC = PARENTOPOA
      JSV = IDSV(JHC)
      NOX = NONLYNOX

      ! Product 1
      NOX = NONLYNOX
      IPR = 1
      ! Compute OA in kg
      OATEMP = Spc(I,J,L,id_OPOG1) * KOMIJL(IPR,JSV) * MOTEMP
      ! Compute difference in ug/m3
      EQLBDIFF = ABS( OATEMP - Spc(I,J,L,id_OPOA1) )* CONVFAC
      ! Assess error
      IF ( EQLBDIFF > ACCEPTERRORUG ) THEN
         WRITE(*,*) 'EQLB Problem NOX, PR, JSV', NOX, IPR, JSV, &
                    ' in box ', I, J, L
      ENDIF

      ! Product 2
      NOX = NONLYNOX
      IPR = 2
      ! Compute OA in kg
      OATEMP = Spc(I,J,L,id_OPOG2) * KOMIJL(IPR,JSV) * MOTEMP
      ! Compute difference in ug/m3
      EQLBDIFF = ABS( OATEMP - Spc(I,J,L,id_OPOA2) )* CONVFAC
      ! Assess error
      IF ( EQLBDIFF > ACCEPTERRORUG ) THEN
         WRITE(*,*) 'EQLB Problem NOX, PR, JSV', NOX, IPR, JSV, &
                    ' in box ', I, J, L
      ENDIF

   ENDIF ! OPOA

   ! Free pointer
   Spc => NULL()

 END SUBROUTINE CHECK_EQLB
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: save_oaginit
!
! !DESCRIPTION: Subroutine SAVE\_OAGINIT saves total SOA+SOG before
!  partitioning for diagnostic purposes. Units are the same as the STT array
!  ([kg] or [kgC per box]).
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE SAVE_OAGINIT( State_Chm, State_Grid, State_Met )
!
! !USES:
!
   USE State_Chm_Mod,      ONLY : ChmState
   USE State_Grid_Mod,     ONLY : GrdState
   USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
   TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
   TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
   TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY:
!  17 May 2010 - H.O.T. Pye - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER           :: I, J, L, NOX, JHC, JSV

   ! Pointers
   REAL(fp), POINTER :: Spc(:,:,:,:)

   !=================================================================
   ! SAVE_OAGINIT starts here
   !=================================================================

   ! Point to chemical species array [kg]
   Spc => State_Chm%Species

   !$OMP PARALLEL DO       &
   !$OMP DEFAULT( SHARED ) &
   !$OMP PRIVATE( NOX, JHC, JSV, I, J, L )
   DO L = 1, State_Grid%MaxChemLev
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

      !----------------------------------------------------
      ! Semivolatile 1: monoterpene and sesquiterpene SOA+SOG in kg
      !----------------------------------------------------
      JHC = PARENTMTPA
      JSV = IDSV(JHC)
      OAGINITSAVE(I,J,L,1,JSV) = Spc(I,J,L,id_TSOA1) + Spc(I,J,L,id_TSOG1)
      OAGINITSAVE(I,J,L,2,JSV) = Spc(I,J,L,id_TSOA2) + Spc(I,J,L,id_TSOG2)
      OAGINITSAVE(I,J,L,3,JSV) = Spc(I,J,L,id_TSOA3) + Spc(I,J,L,id_TSOG3)
      OAGINITSAVE(I,J,L,4,JSV) = Spc(I,J,L,id_TSOA0) + Spc(I,J,L,id_TSOG0)

      !------------------------------------------------------------------------
      ! Prior to 7/15/19:
      ! Remove isoprene from VBS (mps, 7/15/19
      !!----------------------------------------------------
      !! Semivolatile 2: isoprene SOA+SOG in kg
      !!----------------------------------------------------
      !JHC = PARENTISOP
      !JSV = IDSV(JHC)
      !OAGINITSAVE(I,J,L,1,JSV) = Spc(I,J,L,id_ISOA1) + Spc(I,J,L,id_ISOG1)
      !OAGINITSAVE(I,J,L,2,JSV) = Spc(I,J,L,id_ISOA2) + Spc(I,J,L,id_ISOG2)
      !OAGINITSAVE(I,J,L,3,JSV) = Spc(I,J,L,id_ISOA3) + Spc(I,J,L,id_ISOG3)
      !------------------------------------------------------------------------

      !----------------------------------------------------
      ! Semivolatile 3: lumped arom/IVOC: total SOA+SOG in kg
      !----------------------------------------------------
      JHC = PARENTBENZ
      JSV = IDSV(JHC)
      ! High NOx
      OAGINITSAVE(I,J,L,1,JSV) = Spc(I,J,L,id_ASOA1) + Spc(I,J,L,id_ASOG1)
      OAGINITSAVE(I,J,L,2,JSV) = Spc(I,J,L,id_ASOA2) + Spc(I,J,L,id_ASOG2)
      OAGINITSAVE(I,J,L,3,JSV) = Spc(I,J,L,id_ASOA3) + Spc(I,J,L,id_ASOG3)
      ! Low NOx
      OAGINITSAVE(I,J,L,4,JSV) = Spc(I,J,L,id_ASOAN)

      !----------------------------------------------------
      ! POA: total POA+POG in kgC (if semivol simulation)
      !----------------------------------------------------
      IF ( id_POA1 > 0 ) THEN
         JHC = PARENTPOA
         JSV = IDSV(JHC)
         OAGINITSAVE(I,J,L,1,JSV) = Spc(I,J,L,id_POA1) + Spc(I,J,L,id_POG1)
         OAGINITSAVE(I,J,L,2,JSV) = Spc(I,J,L,id_POA2) + Spc(I,J,L,id_POG2)
      ENDIF

      !----------------------------------------------------
      ! OPOA: total SOA+SOG in kgC (if semivol simulation)
      !----------------------------------------------------
      IF ( id_OPOA1 > 0 ) THEN
         JHC = PARENTOPOA
         JSV = IDSV(JHC)
         OAGINITSAVE(I,J,L,1,JSV) = Spc(I,J,L,id_OPOA1) + Spc(I,J,L,id_OPOG1)
         OAGINITSAVE(I,J,L,2,JSV) = Spc(I,J,L,id_OPOA2) + Spc(I,J,L,id_OPOG2)
      ENDIF

   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Free pointer
   Spc => NULL()

 END SUBROUTINE SAVE_OAGINIT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_mb
!
! !DESCRIPTION: Subroutine CHECK\_MB checks total SOA+SOG mass balance for
!  diagnostic/debugging purposes. Units are the same as the STT array ([kg] or
!  [kgC per box]). Routine also prints helpful budget info. Created by Havala
!  Pye (5/18/10).
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE CHECK_MB( Input_Opt, State_Chm, State_Grid, State_Met )
!
! !USES:
!
   USE Input_Opt_Mod,      ONLY : OptInput
   USE State_Chm_Mod,      ONLY : ChmState
   USE State_Grid_Mod,     ONLY : GrdState
   USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
   TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
   TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
   TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY:
!  18 May 2010 - H.O.T. Pye - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
   ! [kg/kg] of acceptable error, 1e-12_fp = 1e-10_fp %
   !REAL(fp), PARAMETER  :: ACCEPTERROR = 1e-12_fp
   ! more strict (hotp 5/26/10): 1e-12_fp %
   REAL(fp), PARAMETER :: ACCEPTERROR = 1e-14_fp
!
! !LOCAL VARIABLES:
!
   ! Scalars
   LOGICAL              :: prtDebug
   INTEGER              :: I, J, L, NOX, JHC, JSV, IPR
   REAL(fp)             :: TEMPDELTA(MNOX,MPROD)
   REAL(fp)             :: TEMPSOAG
   REAL(fp)             :: MBDIFF

   ! Pointers
   REAL(fp), POINTER    :: Spc(:,:,:,:)

   !=================================================================
   ! CHECK_MB starts here
   !=================================================================

   ! Point to chemical species array [kg]
   Spc => State_Chm%Species

   ! Do we have to print debug output?
   prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

   ! run in serial now (hotp 6/5/10)
   ! Make parallel again (mpayer, 9/14/11)
   !$OMP PARALLEL DO       &
   !$OMP DEFAULT( SHARED ) &
   !$OMP PRIVATE( NOX,         JHC,        JSV,    IPR ) &
   !$OMP PRIVATE( TEMPDELTA,   TEMPSOAG,   MBDIFF      ) &
   !$OMP PRIVATE( I, J, L                              )
   DO L = 1, State_Grid%MaxChemLev
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

      !----------------------------------------------------
      ! Semivolatile 1: terpene SOA+SOG in kg
      !----------------------------------------------------
      ! LIMO-MTPO bug/typo fix (hotp 5/26/10)
      TEMPDELTA = 0e+0_fp
      JHC = PARENTMTPA
      JSV = IDSV(JHC)
      DO IPR = 1, NPROD(JSV)
      DO NOX = 1, NNOX(JSV)
         TEMPDELTA(NOX,IPR) = &
            DELTASOGSAVE(I,J,L,NOX,PARENTMTPA)*ALPHA(NOX,IPR,PARENTMTPA) &
          + DELTASOGSAVE(I,J,L,NOX,PARENTLIMO)*ALPHA(NOX,IPR,PARENTLIMO) &
          + DELTASOGSAVE(I,J,L,NOX,PARENTMTPO)*ALPHA(NOX,IPR,PARENTMTPO) &
          + DELTASOGSAVE(I,J,L,NOX,PARENTSESQ)*ALPHA(NOX,IPR,PARENTSESQ)
      ENDDO
      ENDDO

      ! Product 1, C* = 1 ug/m3
      IPR      = 1
      TEMPSOAG = OAGINITSAVE(I,J,L,IPR,JSV) + SUM(TEMPDELTA(:,IPR))
      MBDIFF   = ABS( TEMPSOAG - ( Spc(I,J,L,id_TSOA1) + Spc(I,J,L,id_TSOG1) ))
      MBDIFF   = MBDIFF/TEMPSOAG ! convert to fractional error
      IF ( prtDebug .and. MBDIFF > ACCEPTERROR ) THEN
         WRITE(*,*) 'MB Problem with NOX, IPR, JSV:', NOX, IPR, JSV, &
                    'in box ', I, J, L
         print*,'CK_MB ',NOX,IPR,JSV, &
                TEMPSOAG,MBDIFF,OAGINITSAVE(I,J,L,IPR,JSV), &
                TEMPDELTA(:,IPR)
      ENDIF

      ! Product 2, C* = 10 ug/m3
      IPR      = 2
      TEMPSOAG = OAGINITSAVE(I,J,L,IPR,JSV) + SUM(TEMPDELTA(:,IPR))
      MBDIFF   = ABS( TEMPSOAG - ( Spc(I,J,L,id_TSOA2) +  Spc(I,J,L,id_TSOG2) ))
      MBDIFF   = MBDIFF/TEMPSOAG ! convert to fractional error
      IF ( prtDebug .and. MBDIFF > ACCEPTERROR ) THEN
         WRITE(*,*) 'MB Problem with NOX, IPR, JSV:', NOX, IPR, JSV, &
                    'in box ', I, J, L
         print*,'CK_MB ',NOX,IPR,JSV, &
                TEMPSOAG,MBDIFF,OAGINITSAVE(I,J,L,IPR,JSV), &
                TEMPDELTA(:,IPR)
      ENDIF

      ! Product 3, C* = 100 ug/m3
      IPR      = 3
      TEMPSOAG = OAGINITSAVE(I,J,L,IPR,JSV) + SUM(TEMPDELTA(:,IPR))
      MBDIFF   = ABS( TEMPSOAG - ( Spc(I,J,L,id_TSOA3) + Spc(I,J,L,id_TSOG3) ))
      MBDIFF   = MBDIFF/TEMPSOAG ! convert to fractional error
      IF ( prtDebug .and. MBDIFF > ACCEPTERROR ) THEN
         WRITE(*,*) 'MB Problem with NOX, IPR, JSV:', NOX, IPR, JSV, &
                    'in box ', I, J, L
         print*,'CK_MB ',NOX,IPR,JSV, &
                TEMPSOAG,MBDIFF,OAGINITSAVE(I,J,L,IPR,JSV), &
                TEMPDELTA(:,IPR)
      ENDIF

      ! Product 4, C* = 0.1 ug/m3
      IPR      = 4
      TEMPSOAG = OAGINITSAVE(I,J,L,IPR,JSV) + SUM(TEMPDELTA(:,IPR))
      MBDIFF   = ABS( TEMPSOAG - ( Spc(I,J,L,id_TSOA0) + Spc(I,J,L,id_TSOG0) ))
      MBDIFF   = MBDIFF/TEMPSOAG ! convert to fractional error
      IF ( prtDebug .and. MBDIFF > ACCEPTERROR ) THEN
         WRITE(*,*) 'MB Problem with NOX, IPR, JSV:', NOX, IPR, JSV, &
                    'in box ', I, J, L
         print*,'CK_MB ',NOX,IPR,JSV, &
                TEMPSOAG,MBDIFF,OAGINITSAVE(I,J,L,IPR,JSV), &
                TEMPDELTA(:,IPR)
      ENDIF

      !------------------------------------------------------------------------
      ! Prior to 7/15/19:
      ! Remove isoprene from VBS (mps, 7/15/19)
      !!----------------------------------------------------
      !! Semivolatile 2: isoprene SOA+SOG in kg
      !!----------------------------------------------------
      !TEMPDELTA = 0e+0_fp
      !JHC = PARENTISOP
      !JSV = IDSV(JHC)
      !DO IPR = 1, NPROD(JSV)
      !DO NOX = 1, NNOX(JSV)
      !   TEMPDELTA(NOX,IPR) =
      !   DELTASOGSAVE(I,J,L,NOX,PARENTISOP)*ALPHA(NOX,IPR,PARENTISOP)
      !ENDDO
      !ENDDO
      !
      !! Product 1, C* = 1 ug/m3
      !IPR      = 1
      !TEMPSOAG = OAGINITSAVE(I,J,L,IPR,JSV) + SUM(TEMPDELTA(:,IPR))
      !MBDIFF   = ABS( TEMPSOAG - ( Spc(I,J,L,id_ISOA1) + Spc(I,J,L,id_ISOG1) ))
      !MBDIFF   = MBDIFF/TEMPSOAG ! convert to fractional error
      !IF ( prtDebug .and. MBDIFF > ACCEPTERROR ) THEN
      !   WRITE(*,*) 'MB Problem with NOX, IPR, JSV:', NOX, IPR, JSV, &
      !              'in box ', I, J, L
      !   print*,'CK_MB ',NOX,IPR,JSV, &
      !          TEMPSOAG,MBDIFF,OAGINITSAVE(I,J,L,IPR,JSV), &
      !          TEMPDELTA(:,IPR)
      !   print*,'DELSOGSAVE NOx=1',DELTASOGSAVE(I,J,L,1,5)
      !   print*,'DELSOGSAVE NOx=2',DELTASOGSAVE(I,J,L,2,5)
      !   print*,'DELSOGSAVE NOx=3',DELTASOGSAVE(I,J,L,3,5)
      !   print*,'Spc',Spc(I,J,L,id_ISOA1),Spc(I,J,L,id_ISOG1)
      !   print*,'NNOX',NNOX(JSV)
      !   print*,'strat?',State_Met%InStratosphere(I,J,L)
      !ENDIF
      !
      !! Product 2, C* = 10 ug/m3
      !IPR      = 2
      !TEMPSOAG = OAGINITSAVE(I,J,L,IPR,JSV) + SUM(TEMPDELTA(:,IPR))
      !MBDIFF   = ABS( TEMPSOAG - ( Spc(I,J,L,id_ISOA2) + Spc(I,J,L,id_ISOG2) ))
      !MBDIFF   = MBDIFF/TEMPSOAG ! convert to fractional error
      !IF ( prtDebug .and. MBDIFF > ACCEPTERROR ) THEN
      !   WRITE(*,*) 'MB Problem with NOX, IPR, JSV:', NOX, IPR, JSV, &
      !              'in box ', I, J, L
      !   print*,'CK_MB ',NOX,IPR,JSV, &
      !          TEMPSOAG,MBDIFF,OAGINITSAVE(I,J,L,IPR,JSV), &
      !          TEMPDELTA(:,IPR)
      !   print*,'DELSOGSAVE NOx=1',DELTASOGSAVE(I,J,L,1,5)
      !   print*,'DELSOGSAVE NOx=2',DELTASOGSAVE(I,J,L,2,5)
      !   print*,'DELSOGSAVE NOx=3',DELTASOGSAVE(I,J,L,3,5)
      !   print*,'Spc',Spc(I,J,L,id_ISOA2),Spc(I,J,L,id_ISOG2)
      !   print*,'NNOX',NNOX(JSV)
      !   print*,'strat?',State_Met%InStratosphere(I,J,L)
      !
      !ENDIF
      !
      !! Product 3, C* = 100 ug/m3
      !IPR      = 3
      !TEMPSOAG = OAGINITSAVE(I,J,L,IPR,JSV) + SUM(TEMPDELTA(:,IPR))
      !MBDIFF   = ABS( TEMPSOAG - ( Spc(I,J,L,id_ISOA3) + Spc(I,J,L,id_ISOG3) ))
      !MBDIFF   = MBDIFF/TEMPSOAG ! convert to fractional error
      !IF ( prtDebug .and. MBDIFF > ACCEPTERROR ) THEN
      !   WRITE(*,*) 'MB Problem with NOX, IPR, JSV:', NOX, IPR, JSV, &
      !              'in box ', I, J, L
      !   print*,'CK_MB ',NOX,IPR,JSV, &
      !          TEMPSOAG,MBDIFF,OAGINITSAVE(I,J,L,IPR,JSV), &
      !          TEMPDELTA(:,IPR)
      !   print*,'DELSOGSAVE NOx=1',DELTASOGSAVE(I,J,L,1,5)
      !   print*,'DELSOGSAVE NOx=2',DELTASOGSAVE(I,J,L,2,5)
      !   print*,'DELSOGSAVE NOx=3',DELTASOGSAVE(I,J,L,3,5)
      !   print*,'Spc',Spc(I,J,L,id_ISOA3),Spc(I,J,L,id_ISOG3)
      !   print*,'NNOX',NNOX(JSV)
      !   print*,'strat?',State_Met%InStratosphere(I,J,L)
      !ENDIF
      !------------------------------------------------------------------------

      !----------------------------------------------------
      ! Semivolatile 3: lumped arom/IVOC: total SOA+SOG in kg
      !----------------------------------------------------
      TEMPDELTA = 0e+0_fp
      JHC = PARENTBENZ
      JSV = IDSV(JHC)
      DO IPR = 1, NPROD(JSV)
      DO NOX = 1, NNOX(JSV)
         TEMPDELTA(NOX,IPR) = &
           DELTASOGSAVE(I,J,L,NOX,PARENTBENZ)*ALPHA(NOX,IPR,PARENTBENZ) &
         + DELTASOGSAVE(I,J,L,NOX,PARENTTOLU)*ALPHA(NOX,IPR,PARENTTOLU) &
         + DELTASOGSAVE(I,J,L,NOX,PARENTXYLE)*ALPHA(NOX,IPR,PARENTXYLE) &
         + DELTASOGSAVE(I,J,L,NOX,PARENTNAP )*ALPHA(NOX,IPR,PARENTNAP )
      ENDDO
      ENDDO

      ! Low NOx
      NOX      = NLOWNOX
      IPR      = 4
      TEMPSOAG = OAGINITSAVE(I,J,L,IPR,JSV) + TEMPDELTA(NOX,IPR)
      MBDIFF   = ABS( TEMPSOAG - Spc(I,J,L,id_ASOAN) )
      MBDIFF   = MBDIFF/TEMPSOAG ! convert to fractional error
      IF ( prtDebug .and. MBDIFF > ACCEPTERROR ) THEN
         WRITE(*,*) 'MB Problem with NOX, IPR, JSV:', NOX, IPR, JSV, &
                    'in box ', I, J, L
         print*,'CK_MB ',NOX,IPR,JSV, &
                TEMPSOAG,MBDIFF,OAGINITSAVE(I,J,L,IPR,JSV), &
                TEMPDELTA(:,IPR)
         !print*,'DELSOGSAVE NOx=2',DELTASOGSAVE(I,J,L,2,6:8)
         !print*,'DELSOGSAVE NOx=2',DELTASOGSAVE(I,J,L,2,11)
         !print*,'Spc',Spc(I,J,L,id_ASOAN)
         !print*,'NNOX',NNOX(JSV)
         !print*,'strat?',State_Met%InStratosphere(I,J,L)
      ENDIF

      ! Debug print to screen
      !IF ( prtDebug ) THEN
      !   IF ( I == 37 .AND. J == 25 .AND. L == 4 ) THEN
      !      print*,'CK_MB ',NOX,IPR,JSV, &
      !             TEMPSOAG,MBDIFF,OAGINITSAVE(I,J,L,NOX,IPR,JSV)
      !       print*,'strat?',State_Met%InStratosphere(I,J,L)
      !   ENDIF
      !ENDIF

      ! HIGH NOx, Product 1
      NOX      = NHIGHNOX
      IPR      = 1
      TEMPSOAG = OAGINITSAVE(I,J,L,IPR,JSV) + TEMPDELTA(NOX,IPR)
      MBDIFF   = ABS( TEMPSOAG - ( Spc(I,J,L,id_ASOA1) + Spc(I,J,L,id_ASOG1) ))
      MBDIFF   = MBDIFF/TEMPSOAG ! convert to fractional error
      IF ( prtDebug .and. MBDIFF > ACCEPTERROR ) THEN
         WRITE(*,*) 'MB Problem with NOX, IPR, JSV:', NOX, IPR, JSV, &
                    'in box ', I, J, L
         print*,'CK_MB ',NOX,IPR,JSV, &
                TEMPSOAG,MBDIFF,OAGINITSAVE(I,J,L,IPR,JSV), &
                TEMPDELTA(:,IPR)
      ENDIF

      ! HIGH NOx, Product 2
      NOX      = NHIGHNOX
      IPR      = 2
      TEMPSOAG = OAGINITSAVE(I,J,L,IPR,JSV) + TEMPDELTA(NOX,IPR)
      MBDIFF   = ABS( TEMPSOAG - ( Spc(I,J,L,id_ASOA2) + Spc(I,J,L,id_ASOG2) ))
      MBDIFF   = MBDIFF/TEMPSOAG ! convert to fractional error
      IF ( prtDebug .and. MBDIFF > ACCEPTERROR ) THEN
         WRITE(*,*) 'MB Problem with NOX, IPR, JSV:', NOX, IPR, JSV, &
                    'in box ', I, J, L
         print*,'CK_MB ',NOX,IPR,JSV, &
                TEMPSOAG,MBDIFF,OAGINITSAVE(I,J,L,IPR,JSV), &
                TEMPDELTA(:,IPR)
      ENDIF

      ! HIGH NOx, Product 3
      NOX      = NHIGHNOX
      IPR      = 3
      TEMPSOAG = OAGINITSAVE(I,J,L,IPR,JSV) + TEMPDELTA(NOX,IPR)
      MBDIFF   = ABS( TEMPSOAG - ( Spc(I,J,L,id_ASOA3) + Spc(I,J,L,id_ASOG3) ))
      MBDIFF   = MBDIFF/TEMPSOAG ! convert to fractional error
      IF ( prtDebug .and. MBDIFF > ACCEPTERROR ) THEN
         WRITE(*,*) 'MB Problem with NOX, IPR, JSV:', NOX, IPR, JSV, &
                    'in box ', I, J, L
         print*,'CK_MB ',NOX,IPR,JSV, &
                TEMPSOAG,MBDIFF,OAGINITSAVE(I,J,L,IPR,JSV), &
                TEMPDELTA(:,IPR)
      ENDIF

      !----------------------------------------------------
      ! POA: total POA+POG in kgC
      !----------------------------------------------------
      ! Note that POA+G is both increased (due to emission)
      ! and decreased (due to conversion to OPOG)
      IF ( id_POA1 > 0 ) THEN
         TEMPDELTA = 0e+0_fp
         JHC = PARENTPOA
         JSV = IDSV(JHC)
         NOX = NONLYNOX
         DO IPR = 1, NPROD(JSV)
            TEMPDELTA(NOX,IPR) = DELTASOGSAVE(I,J,L,IPR,PARENTPOA ) * &
                                 ALPHA(NOX,IPR,PARENTPOA ) - &
                                 DELTASOGSAVE(I,J,L,IPR,PARENTOPOA) * &
                                 ALPHA(NOX,IPR,PARENTOPOA)
         ENDDO

         ! Only NOx, Product 1
         IPR      = 1
         TEMPSOAG = OAGINITSAVE(I,J,L,IPR,JSV) + TEMPDELTA(NOX,IPR)
         MBDIFF   = ABS( TEMPSOAG - ( Spc(I,J,L,id_POA1) + Spc(I,J,L,id_POG1) ))
         MBDIFF   = MBDIFF/TEMPSOAG ! convert to fractional error
         IF ( prtDebug .and. MBDIFF > ACCEPTERROR ) THEN
            WRITE(*,*) 'MB Problem with NOX, IPR, JSV:', NOX, IPR, &
                       JSV, 'in box ', I, J, L
            print*,'CK_MB ',NOX,IPR,JSV, &
                   TEMPSOAG,MBDIFF,OAGINITSAVE(I,J,L,IPR,JSV), &
                   TEMPDELTA(:,IPR)
         ENDIF

         ! Only NOx, Product 2
         IPR      = 2
         TEMPSOAG = OAGINITSAVE(I,J,L,IPR,JSV) + TEMPDELTA(NOX,IPR)
         MBDIFF   = ABS( TEMPSOAG - ( Spc(I,J,L,id_POA2) + Spc(I,J,L,id_POG2) ))
         MBDIFF   = MBDIFF/TEMPSOAG ! convert to fractional error
         IF ( prtDebug .and. MBDIFF > ACCEPTERROR ) THEN
            WRITE(*,*) 'MB Problem with NOX, IPR, JSV:', NOX, IPR, &
                       JSV, 'in box ', I, J, L
            print*,'CK_MB ',NOX,IPR,JSV, &
                   TEMPSOAG,MBDIFF,OAGINITSAVE(I,J,L,IPR,JSV), &
                   TEMPDELTA(:,IPR)
         ENDIF
      ENDIF ! POA1

      !----------------------------------------------------
      ! OPOA: total SOA+SOG in kgC
      !----------------------------------------------------
      IF ( id_OPOA1 > 0 ) THEN
         TEMPDELTA = 0e+0_fp
         JHC = PARENTOPOA
         JSV = IDSV(JHC)
         NOX = NONLYNOX
         DO IPR = 1, NPROD(JSV)
            TEMPDELTA(NOX,IPR) = DELTASOGSAVE(I,J,L,IPR,PARENTOPOA) * &
                                 ALPHA(NOX,IPR,PARENTOPOA)
         ENDDO

         ! Only NOx, Product 1
         IPR      = 1
         TEMPSOAG = OAGINITSAVE(I,J,L,IPR,JSV) + TEMPDELTA(NOX,IPR)
         MBDIFF   = ABS( TEMPSOAG - ( Spc(I,J,L,id_OPOA1) + Spc(I,J,L,id_OPOG1) ))
         MBDIFF   = MBDIFF/TEMPSOAG ! convert to fractional error
         IF ( prtDebug .and. MBDIFF > ACCEPTERROR ) THEN
            WRITE(*,*) 'MB Problem with NOX, IPR, JSV:', NOX, IPR, &
                       JSV, 'in box ', I, J, L
            print*,'CK_MB ',NOX,IPR,JSV, &
                   TEMPSOAG,MBDIFF,OAGINITSAVE(I,J,L,IPR,JSV), &
                   TEMPDELTA(:,IPR)
         ENDIF

         ! Only NOx, Product 2
         IPR      = 2
         TEMPSOAG = OAGINITSAVE(I,J,L,IPR,JSV) + TEMPDELTA(NOX,IPR)
         MBDIFF   = ABS( TEMPSOAG - ( Spc(I,J,L,id_OPOA2) + Spc(I,J,L,id_OPOG2) ))
         MBDIFF   = MBDIFF/TEMPSOAG ! convert to fractional error
         IF ( prtDebug .and. MBDIFF > ACCEPTERROR ) THEN
            WRITE(*,*) 'MB Problem with NOX, IPR, JSV:', NOX, IPR, &
                       JSV, 'in box ', I, J, L
            print*,'CK_MB ',NOX,IPR,JSV, &
                   TEMPSOAG,MBDIFF,OAGINITSAVE(I,J,L,IPR,JSV), &
                   TEMPDELTA(:,IPR)
         ENDIF
      ENDIF ! OPOA1

   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Save information in [Tg]
   DO JHC = 1, MHC
   DO NOX = 1, MNOX
      DELTAHCSAVE(NOX,JHC) = DELTAHCSAVE(NOX,JHC) + &
                             1e-9_fp * SUM(DELTASOGSAVE(:,:,:,NOX,JHC))
   ENDDO
   ENDDO

   ! Print diagnostic information to screen
   IF ( prtDebug ) THEN
      print*,'Global cumulative amount reacted in gas phase [Tg]'
      JHC = 1
      print*,'MTPA High NOx Rxn : ', DELTAHCSAVE(1,JHC)
      print*,'MTPA Low  NOx Rxn : ', DELTAHCSAVE(2,JHC)
      print*,'MTPA NO3      Rxn : ', DELTAHCSAVE(3,JHC)
      JHC = 2
      print*,'LIMO High NOx Rxn : ', DELTAHCSAVE(1,JHC)
      print*,'LIMO Low  NOx Rxn : ', DELTAHCSAVE(2,JHC)
      print*,'LIMO NO3      Rxn : ', DELTAHCSAVE(3,JHC)
      JHC = 3
      print*,'MTPO High NOx Rxn : ', DELTAHCSAVE(1,JHC)
      print*,'MTPO Low  NOx Rxn : ', DELTAHCSAVE(2,JHC)
      print*,'MTPO NO3      Rxn : ', DELTAHCSAVE(3,JHC)
      JHC = 4
      print*,'SESQ High NOx Rxn : ', DELTAHCSAVE(1,JHC)
      print*,'SESQ Low  NOx Rxn : ', DELTAHCSAVE(2,JHC)
      print*,'SESQ NO3      Rxn : ', DELTAHCSAVE(3,JHC)
      JHC = 5
      print*,'ISOP OH       Rxn : ', DELTAHCSAVE(1,JHC)
      print*,'ISOP NO3      Rxn : ', DELTAHCSAVE(2,JHC)
      JHC = 6
      print*,'BENZ High NOx Rxn : ', DELTAHCSAVE(1,JHC)
      print*,'BENZ Low  NOx Rxn : ', DELTAHCSAVE(2,JHC)
      JHC = 7
      print*,'TOLU High NOx Rxn : ', DELTAHCSAVE(1,JHC)
      print*,'TOLU Low  NOx Rxn : ', DELTAHCSAVE(2,JHC)
      JHC = 8
      print*,'XYLE High NOx Rxn : ', DELTAHCSAVE(1,JHC)
      print*,'XYLE Low  NOx Rxn : ', DELTAHCSAVE(2,JHC)
      JHC = 11
      print*,'NAP  High NOx Rxn : ', DELTAHCSAVE(1,JHC)
      print*,'NAP  Low  NOx Rxn : ', DELTAHCSAVE(2,JHC)
      JHC = 10
      print*,'POG1 OH       Rxn : ', DELTAHCSAVE(1,JHC)
      print*,'POG2 OH       Rxn : ', DELTAHCSAVE(2,JHC)

      ! Check Spc for debug purposes (hotp 6/4/10)
      !IF ( prtDebug ) THEN
      !   print*,Spc(37,25,4,:)
      !ENDIF

      ! Print diagnostic info about SOA production
      print*,'Aerosol production and evaporation (cumulative kg)'
      JSV = 1
      IPR = 1
      print*,'TSOA1 prod and evap: ', &
         SUM(SPECSOAPROD(:,:,:,IPR,JSV)),SUM(SPECSOAEVAP(:,:,:,IPR,JSV))
      IPR = 2
      print*,'TSOA2 prod and evap: ', &
         SUM(SPECSOAPROD(:,:,:,IPR,JSV)),SUM(SPECSOAEVAP(:,:,:,IPR,JSV))
      IPR = 3
      print*,'TSOA3 prod and evap: ', &
         SUM(SPECSOAPROD(:,:,:,IPR,JSV)),SUM(SPECSOAEVAP(:,:,:,IPR,JSV))
      IPR = 4
      print*,'TSOA0 prod and evap: ', &
         SUM(SPECSOAPROD(:,:,:,IPR,JSV)),SUM(SPECSOAEVAP(:,:,:,IPR,JSV))

      !------------------------------------------------------------------------
      ! Prior to 7/15/19:
      ! Remove isoprene from VBS (mps, 7/15/19)
      !JSV = 2
      !IPR = 1
      !print*,'ISOA1 prod and evap: ', &
      !   SUM(SPECSOAPROD(:,:,:,IPR,JSV)),SUM(SPECSOAEVAP(:,:,:,IPR,JSV))
      !IPR = 2
      !print*,'ISOA2 prod and evap: ', &
      !   SUM(SPECSOAPROD(:,:,:,IPR,JSV)),SUM(SPECSOAEVAP(:,:,:,IPR,JSV))
      !IPR = 3
      !print*,'ISOA3 prod and evap: ', &
      !   SUM(SPECSOAPROD(:,:,:,IPR,JSV)),SUM(SPECSOAEVAP(:,:,:,IPR,JSV))
      !------------------------------------------------------------------------

      JSV = 3
      IPR = 1
      print*,'ASOA1 prod and evap: ', &
         SUM(SPECSOAPROD(:,:,:,IPR,JSV)),SUM(SPECSOAEVAP(:,:,:,IPR,JSV))
      IPR = 2
      print*,'ASOA2 prod and evap: ', &
        SUM(SPECSOAPROD(:,:,:,IPR,JSV)),SUM(SPECSOAEVAP(:,:,:,IPR,JSV))
      IPR = 3
      print*,'ASOA3 prod and evap: ', &
         SUM(SPECSOAPROD(:,:,:,IPR,JSV)),SUM(SPECSOAEVAP(:,:,:,IPR,JSV))
      IPR = 4
      print*,'ASOAN prod and evap: ', &
         SUM(SPECSOAPROD(:,:,:,IPR,JSV)),SUM(SPECSOAEVAP(:,:,:,IPR,JSV))

      JSV = 4
      IPR = 1
      print*,'POA1  prod and evap: ', &
         SUM(SPECSOAPROD(:,:,:,IPR,JSV)),SUM(SPECSOAEVAP(:,:,:,IPR,JSV))
      IPR = 2
      print*,'POA2  prod and evap: ', &
         SUM(SPECSOAPROD(:,:,:,IPR,JSV)),SUM(SPECSOAEVAP(:,:,:,IPR,JSV))

      JSV = 5
      IPR = 1
      print*,'OPOA1 prod and evap: ', &
         SUM(SPECSOAPROD(:,:,:,IPR,JSV)),SUM(SPECSOAEVAP(:,:,:,IPR,JSV))
      IPR = 2
      print*,'OPOA2 prod and evap: ', &
         SUM(SPECSOAPROD(:,:,:,IPR,JSV)),SUM(SPECSOAEVAP(:,:,:,IPR,JSV))

   ENDIF

   ! Free pointer
   Spc => NULL()
   
 END SUBROUTINE CHECK_MB
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_no
!
! !DESCRIPTION: Function GET\_NO returns NO from State\_Chm%Species
! (for coupled runs). (hotp 5/7/2010)
!\\
!\\
! !INTERFACE:
!
 FUNCTION GET_NO( I, J, L, Input_Opt, State_Chm, State_Met ) &
      RESULT( NO_MOLEC_CM3 )
!
! !USES:
!
   USE ERROR_MOD,          ONLY : ERROR_STOP
   USE Input_Opt_Mod,      ONLY : OptInput
   USE State_Chm_Mod,      ONLY : ChmState
   USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
   INTEGER,        INTENT(IN) :: I              ! Longitude index
   INTEGER,        INTENT(IN) :: J              ! Latitude index
   INTEGER,        INTENT(IN) :: L              ! Altitude index
   TYPE(OptInput), INTENT(IN) :: Input_Opt      ! Input Options object
   TYPE(ChmState), INTENT(IN) :: State_Chm      ! Chemistry State object
   TYPE(MetState), INTENT(IN) :: State_Met      ! Meteorology State object
!
! !RETURN VALUE
!
   REAL(fp)                   :: NO_MOLEC_CM3   ! NO conc [molec/cm3]
!
! !REVISION HISTORY:
!  07 May 2010 - H.O.T. Pye - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   REAL(fp) :: MolecRatio ! moles C / moles species
   REAL(fp) :: NO_MW_kg   ! kg NO / mol

   !=================================================================
   ! GET_NO begins here!
   !=================================================================
   IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

      !---------------------
      ! Coupled simulation
      !---------------------

      ! NO is defined only in the chemistry grid
      IF ( State_Met%InChemGrid(I,J,L) ) THEN

         ! Get NO from State_Chm%Species [kg] and convert to [molec/cm3]
         MolecRatio = State_Chm%SpcData(id_NO)%Info%MolecRatio
         NO_MW_kg   = State_Chm%SpcData(id_NO)%Info%emMW_g * 1.e-3_fp

         NO_MOLEC_CM3 = State_Chm%Species(I,J,L,id_NO) &
                        * ( AVO / NO_MW_kg )           &
                        / ( State_Met%AIRVOL(I,J,L)    &
                        * 1e+6_fp * MolecRatio )

      ELSE

         NO_MOLEC_CM3 = 0e+0_fp

      ENDIF

   ELSE

      !---------------------
      ! Invalid sim type!
      !---------------------
      CALL ERROR_STOP( 'Invalid Simulation Type!', &
                       'GET_NO ("carbon_mod.F90")' )

   ENDIF

 END FUNCTION GET_NO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ho2
!
! !DESCRIPTION: Function GET\_HO2 returns HO2 from State\_Chm%Species
! (for coupled runs). Created by Havala Pye (5/7/2010).
!\\
!\\
! !INTERFACE:
!
 FUNCTION GET_HO2( I, J, L, Input_Opt, State_Chm, State_Met ) &
      RESULT( HO2_MOLEC_CM3 )
!
! !USES:
!
   USE ERROR_MOD,          ONLY : ERROR_STOP
   USE Input_Opt_Mod,      ONLY : OptInput
   USE State_Chm_Mod,      ONLY : ChmState
   USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
   INTEGER,        INTENT(IN) :: I               ! Longitude index
   INTEGER,        INTENT(IN) :: J               ! Latitude index
   INTEGER,        INTENT(IN) :: L               ! Altitude index
   TYPE(OptInput), INTENT(IN) :: Input_Opt       ! Input Options object
   TYPE(ChmState), INTENT(IN) :: State_Chm       ! Chemistry State object
   TYPE(MetState), INTENT(IN) :: State_Met       ! Meteorology State object
!
! !RETURN VALUE
!
   REAL(fp)                   :: HO2_MOLEC_CM3   ! HO2 conc [molec/cm3]
!
! !REVISION HISTORY:
!  13 May 2010 - H.O.T. Pye - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   REAL(fp) :: MolecRatio ! moles C / moles species
   REAL(fp) :: HO2_MW_kg  ! kg HO2 / mol

   !=================================================================
   ! GET_HO2 begins here!
   !=================================================================
   IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

      !---------------------
      ! Coupled simulation
      !---------------------

      ! HO2 is defined only in the chemistry grid
      IF ( State_Met%InChemGrid(I,J,L) ) THEN

         ! Get HO2 from State_Chm%Species [kg] and convert to [molec/cm3]
         MolecRatio = State_Chm%SpcData(id_HO2)%Info%MolecRatio
         HO2_MW_kg  = State_Chm%SpcData(id_HO2)%Info%emMW_g*1.e-3_fp

         HO2_MOLEC_CM3 = State_Chm%Species(I,J,L,id_HO2) &
                         * ( AVO / HO2_MW_kg ) &
                         / ( State_Met%AIRVOL(I,J,L) &
                         * 1e+6_fp * MolecRatio )

      ELSE

         HO2_MOLEC_CM3 = 0e+0_fp

      ENDIF

   ELSE

      !---------------------
      ! Invalid sim type!
      !---------------------
      CALL ERROR_STOP( 'Invalid Simulation Type!', &
                       'GET_HO2 ("carbon_mod.F90")' )

   ENDIF

 END FUNCTION GET_HO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_isopno3
!
! !DESCRIPTION: Modification of GET\_DOH that returns the amount of isoprene
!  [kgC] that has reacted with NO3 during the last chemistry time step.
!  (hotp 5/22/10)
!\\
!\\
! !INTERFACE:
!
 FUNCTION GET_ISOPNO3( I, J, L, Input_Opt, State_Chm, State_Met ) &
      RESULT( ISOPNO3 )
!
! !USES:
!
   USE ERROR_MOD,          ONLY : ERROR_STOP
   USE Input_Opt_Mod,      ONLY : OptInput
   USE PhysConstants,      ONLY : AVO
   USE State_Chm_Mod,      ONLY : ChmState
   USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
   INTEGER,        INTENT(IN) :: I           ! Longitude index
   INTEGER,        INTENT(IN) :: J           ! Latitude index
   INTEGER,        INTENT(IN) :: L           ! Altitude index
   TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
   TYPE(ChmState), INTENT(IN) :: State_Chm   ! Chemistry State object
   TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !RETURN VALUE
!
   REAL(fp)                   :: ISOPNO3     ! ISOP replaced w/ NO3 [kg C]
!
! !REVISION HISTORY:
!  22 May 2010 - H.O.T. Pye - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   REAL(fp) :: MolecRatio     ! molec C / moles ISOP
   REAL(fp) :: ISOP_MW_kg     ! kg C ISOP      / mol C
   REAL(fp) :: LISOPNO3_MW_kg ! kg C LISOPONO3 / mol Ckg
      
   !=================================================================
   ! GET_ISOPNO3 begins here!
   !=================================================================

   IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

      !--------------------
      ! Coupled simulation
      !--------------------

      ! Test if we are in the chemistry grid
      IF ( State_Met%InChemGrid(I,J,L) ) THEN

         !-----------------------------------------------------------
         ! Get ISOPNO3 (ISOP list to NO3) from State_Chm%Species
         ! [kg ISOPNO3] and convert to [kg C ISOP]
         !-----------------------------------------------------------
         MolecRatio     = State_Chm%SpcData(id_ISOP)%Info%MolecRatio
         ISOP_MW_kg     = State_Chm%SpcData(id_ISOP)%Info%emMW_g * 1.e-3_fp
         LISOPNO3_MW_kg = State_Chm%SpcData(id_LISOPNO3)%Info%emMW_g * 1.e-3_fp

         ISOPNO3 = State_Chm%Species(I,J,L,id_LISOPNO3) &
                   * ( AVO / LISOPNO3_MW_kg ) &
                   / ( AVO / ISOP_MW_kg     ) * MolecRatio

      ELSE

         ! Otherwise set ISOPNO3=0
         ISOPNO3 = 0e+0_fp
         
      ENDIF

   ELSE IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

      !--------------------
      ! Offline simulation
      !--------------------

      ! ISOP from NO3 not is yet supported for
      ! offline aerosol simulations, set DOH=0
      ISOPNO3 = 0e+0_fp

   ELSE

      !--------------------
      ! Invalid sim type!
      !--------------------
      CALL ERROR_STOP( 'Invalid simulation type!', &
                       'GET_ISOPNO3 ("carbon_mod.F90")' )

   ENDIF

 END FUNCTION GET_ISOPNO3
!EOC
#ifdef TOMAS
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: sakamoto_size
!
! !DESCRIPTION: Function SAKAMOTO\_SIZE finds the peak and standard
! deviation of a lognormal distribution parameterized to fit the number
! size distribution of the fires in the gridbox 24 hours downwind
! of the current fires using K. M. Sakamoto's parameterization from her
! 2016 paper. Then uses these parameters to seperate the size distributions
! into lognormal bins.
!\\
!\\
! !INTERFACE:
!
 FUNCTION SAKAMOTO_SIZE( State_Grid, State_Met, FIRE_NUM,  &
                         OCPIBULKEMIS, BCPIBULKEMIS,       &
                         OCPOBULKEMIS, BCPOBULKEMIS, AREA) &
      RESULT ( VALUE )
!
! !USES:
!
   USE TOMAS_MOD,      ONLY : IBINS, Xk
   USE State_Grid_Mod, ONLY : GrdState
   USE State_Met_Mod,  ONLY : MetState
   USE TOMAS_MOD,      ONLY : AVGMASS
!
! !INPUT PARAMETERS:
!
   TYPE(GrdState),  INTENT(IN)  :: State_Grid ! Grid State object
   TYPE(MetState),  INTENT(IN)  :: State_Met  ! Meteorology State object
   REAL(sp),        INTENT(IN)  :: FIRE_NUM(State_Grid%NX,State_Grid%NY)
   REAL(fp),        INTENT(IN)  :: OCPIBULKEMIS(State_Grid%NX,State_Grid%NY) ![kgm-2s-1]
   REAL(fp),        INTENT(IN)  :: BCPIBULKEMIS(State_Grid%NX,State_Grid%NY) ![kgm-2s-1]
   REAL(fp),        INTENT(IN)  :: OCPOBULKEMIS(State_Grid%NX,State_Grid%NY) ![kgm-2s-1]
   REAL(fp),        INTENT(IN)  :: BCPOBULKEMIS(State_Grid%NX,State_Grid%NY) ![kgm-2s-1]
   REAL(fp),        INTENT(IN)  :: AREA(State_Grid%NX,State_Grid%NY) ![m2]
!
! !RETURN VALUES:
!
   REAL(fp)                     :: VALUE(State_Grid%NX,State_Grid%NY,IBINS,4)
!
! !REVISION HISTORY:
!  01 Jan 2017 - E. Ramnarine - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER             :: I, J, K
   REAL(fp)            :: EMIS
   REAL(fp)            :: W10M, PBLH
   REAL(fp)            :: dM_dxdz
   REAL(fp)            :: Dpm, sig
   REAL(fp)            :: Dl, Dh, Dk
   REAL(fp)            :: NUM_FRAC(IBINS)
   REAL(fp)            :: MASS(IBINS)
   REAL(fp)            :: MASS_FRAC(IBINS)
   REAL(fp), PARAMETER :: pi=3.14159

   ! conditions
   REAL(fp), PARAMETER            :: D_0=100.0 ![nm] doing 100, 150 with sig_0=2.0
   REAL(fp), PARAMETER            :: sig_0=2.0 !doing 1.6, 2.0 with D_0=100.0
   REAL(fp), PARAMETER            :: time= 720.0 !1440.0 ![min] !720 is for a time sensitivity

   ! constants from paper
   REAL(fp), PARAMETER            :: A_d=84.58
   REAL(fp), PARAMETER            :: b_d=0.4191
   REAL(fp), PARAMETER            :: c_d=0.4870

   REAL(fp), PARAMETER            :: A_sig=0.2390
   REAL(fp), PARAMETER            :: b_sig=0.1889
   REAL(fp), PARAMETER            :: c_sig=0.3540

   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      ! Kim's parameterization
      EMIS = (1.8 * (OCPIBULKEMIS(I,J) + OCPOBULKEMIS(I,J)) + &
              BCPIBULKEMIS(I,J) + BCPOBULKEMIS(I,J)) !total emissions [kgm-2s-1]

      IF ( EMIS == 0.0 ) THEN
         DO K = 1, IBINS
            VALUE(I,J,K,1) = 0.0
            VALUE(I,J,K,2) = 0.0
            VALUE(I,J,K,3) = 0.0
            VALUE(I,J,K,4) = 0.0
         ENDDO
      ELSE
         W10M = SQRT( State_Met%U10M(I,J)**2 &
              +       State_Met%V10M(I,J)**2 ) * 60.0 ![m/min]
         PBLH = State_Met%PBLH(I,J) ![m]
         !avg emission per fire [kg/min] (note 320 is from regridding 0.25x0.25 to 4x5)
         EMIS = ( EMIS * AREA(I,J) * 60.0 ) / &
                ( max(1.0, FIRE_NUM(I,J) * 320.0) ) !no plume mixing case
         !EMIS = EMIS * AREA(I,J) * 60.0 !complete plume mixing case (no longer in use)
         dM_dxdz = EMIS / ( max(W10M, 2.0) * max(PBLH, 10.0) ) ![kgm-2]

         Dpm = D_0 + (A_D * (dM_dxdz)**b_D * (time)**c_D)
         sig = sig_0 +(A_sig * dM_dxdz**b_sig * time**c_sig * (1.2-sig_0))
         sig = max(sig, 1.2)

         !Dpm = D_0 ! no coag case
         !sig = sig_0 !no coag case

         DO K = 1, IBINS
            ! splitting into size bins (seinfeld and pandis eq 8.54)
            !Calculate diameter of this size bin
            Dl=1.0e+9*((6.0*Xk(K))/(1400.0*3.14))**0.3333
            Dh=1.0e+9*((6.0*Xk(K+1))/(1400.0*3.14))**0.3333
            Dk=sqrt(Dl*Dh)

            !Calculate number fraction
            NUM_FRAC(K) = 1.0 / ( sqrt(2.0*pi)*Dk*log(sig) ) * &
                          exp( -( (log(Dk/Dpm))**2.0 / &
                          (2.0 * (log(sig))**2.0) )) &
                          * (Dh - Dl)
            MASS(K) = NUM_FRAC(K) * AVGMASS(K)
         ENDDO

         DO K = 1, IBINS
            !Calculate mass fraction
            MASS_FRAC(K) = MASS(K) / SUM( MASS(:) )

            VALUE(I,J,K,1) = BCPIBULKEMIS(I,J) * MASS_FRAC(K)
            VALUE(I,J,K,2) = BCPOBULKEMIS(I,J) * MASS_FRAC(K)
            VALUE(I,J,K,3) = 1.8 * OCPIBULKEMIS(I,J) * MASS_FRAC(K)
            VALUE(I,J,K,4) = 1.8 * OCPOBULKEMIS(I,J) * MASS_FRAC(K)

         ENDDO

      ENDIF

   ENDDO
   ENDDO

 END FUNCTION SAKAMOTO_SIZE
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_carbon
!
! !DESCRIPTION: Subroutine INIT\_CARBON initializes all module arrays.
!  (rjp, bmy, 4/1/04, 12/19/09)
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE INIT_CARBON( Input_Opt, State_Chm, State_Diag, State_Grid, RC )
!
! !USES:
!
   USE ErrCode_Mod
   USE ERROR_MOD,          ONLY : ALLOC_ERR, ERROR_STOP
   USE Input_Opt_Mod,      ONLY : OptInput
   USE State_Chm_Mod,      ONLY : Ind_
   USE State_Chm_Mod,      ONLY : ChmState
   USE State_Diag_Mod,     ONLY : DgnState
   USE State_Grid_Mod,     ONLY : GrdState
#ifdef TOMAS
   USE TOMAS_MOD,          ONLY : IBINS
#endif
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
   TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
   TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT) :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  01 Apr 2004 - R. Park - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER  :: AS, INDICES(4), YYYYMMDD, HHMMSS, N
   REAL(fp) :: COORDS(4), TAU

   !=================================================================
   ! INIT_CARBON begins here!
   !=================================================================

   ! Assume success
   RC = GC_SUCCESS

   ! Define species IDs
   id_ASOG1    = IND_('ASOG1'   )
   id_ASOG2    = IND_('ASOG2'   )
   id_ASOG3    = IND_('ASOG3'   )
   id_ASOA1    = IND_('ASOA1'   )
   id_ASOA2    = IND_('ASOA2'   )
   id_ASOA3    = IND_('ASOA3'   )
   id_ASOAN    = IND_('ASOAN'   )
   id_AW1      = IND_('AW1'     )
   id_BCPI     = IND_('BCPI'    )
   id_BCPO     = IND_('BCPO'    )
   id_BENZ     = IND_('BENZ'    )
   id_ECIL1    = IND_('ECIL1'   )
   id_ECOB1    = IND_('ECOB1'   )
   id_HO2      = IND_('HO2'     )
   id_ISOP     = IND_('ISOP'    )
   id_LIMO     = IND_('LIMO'    )
   id_MTPA     = IND_('MTPA'    )
   id_MTPO     = IND_('MTPO'    )
   id_NAP      = IND_('NAP'     )
   id_NK1      = IND_('NK1'     )
   id_NH4      = IND_('NH4'     )
   id_NO       = IND_('NO'      )
   id_NO3      = IND_('NO3'     )
   id_OCIL1    = IND_('OCIL1'   )
   id_OCOB1    = IND_('OCOB1'   )
   id_O3       = IND_('O3'      )
   id_OH       = IND_('OH'      )
   id_OCPO     = IND_('OCPO'    )
   id_OCPI     = IND_('OCPI'    )
   id_OPOA1    = IND_('OPOA1'   )
   id_OPOG1    = IND_('OPOG1'   )
   id_OPOA2    = IND_('OPOA2'   )
   id_OPOG2    = IND_('OPOG2'   )
   id_POA1     = IND_('POA1'    )
   id_POA2     = IND_('POA2'    )
   id_POG1     = IND_('POG1'    )
   id_POG2     = IND_('POG2'    )
   id_SOAP     = IND_('SOAP'    )
   id_SOAS     = IND_('SOAS'    )
   id_TOLU     = IND_('TOLU'    )
   id_TSOA0    = IND_('TSOA0'   )
   id_TSOA1    = IND_('TSOA1'   )
   id_TSOA2    = IND_('TSOA2'   )
   id_TSOA3    = IND_('TSOA3'   )
   id_TSOG0    = IND_('TSOG0'   )
   id_TSOG1    = IND_('TSOG1'   )
   id_TSOG2    = IND_('TSOG2'   )
   id_TSOG3    = IND_('TSOG3'   )
   id_XYLE     = IND_('XYLE '   )
   id_LBRO2N   = IND_('LBRO2N'  )
   id_LBRO2H   = IND_('LBRO2H'  )
   id_LTRO2N   = IND_('LTRO2N'  )
   id_LTRO2H   = IND_('LTRO2H'  )
   id_LXRO2N   = IND_('LXRO2N'  )
   id_LXRO2H   = IND_('LXRO2H'  )
   id_LNRO2N   = IND_('LNRO2N'  )
   id_LNRO2H   = IND_('LNRO2H'  )
   id_LISOPOH  = IND_('LISOPOH' )
   id_LISOPNO3 = IND_('LISOPNO3')

   ! Some parent hydrocarbons are lumped together into 1 or more
   ! semivolatiles. Map the parent HC to lumped semivolatiles here
   ! (hotp 5/13/10)
   ! mono + sesq
   IDSV(PARENTMTPA) = 1
   IDSV(PARENTLIMO) = 1
   IDSV(PARENTMTPO) = 1
   IDSV(PARENTSESQ) = 1
   ! isoprene
   IDSV(PARENTISOP) = 2
   ! Lumped arom/IVOC
   IDSV(PARENTBENZ) = 3
   IDSV(PARENTTOLU) = 3
   IDSV(PARENTXYLE) = 3
   IDSV(PARENTNAP ) = 3
   ! More individuals
   IDSV(PARENTPOA ) = 4
   IDSV(PARENTOPOA) = 5

   ! Define number of products per semivolatile (hotp 5/14/10)
   NPROD(IDSV(PARENTMTPA)) = 4 ! 3 add C*=0.1 product (hotp 6/12/10)
   NPROD(IDSV(PARENTISOP)) = 3
   NPROD(IDSV(PARENTBENZ)) = 4
   NPROD(IDSV(PARENTPOA )) = 2
   NPROD(IDSV(PARENTOPOA)) = 2
   ! Check to make sure NPROD doesn't exceed max
   IF ( MAXVAL(NPROD(:)) > MPROD ) THEN
      CALL ERROR_STOP('Too many PRODs per SV','carbon_mod.F90')
   ENDIF

   ! Define number of NOx/Ox conditions per semivolatile
   ! (hotp 5/14/10)
   NNOX(IDSV(PARENTMTPA)) = 3 ! high NOx, low NOx, NO3
   NNOX(IDSV(PARENTISOP)) = 2 ! low  NOx, NO3
   NNOX(IDSV(PARENTBENZ)) = 2 ! high NOx, low NOx
   NNOX(IDSV(PARENTPOA )) = 1 ! just OH
   NNOX(IDSV(PARENTOPOA)) = 1 ! just OH
   ! Check to make sure NNOx doesn't exceed max
   IF ( MAXVAL(NNOX(:)) > MNOX ) THEN
      CALL ERROR_STOP('Too many NOx levels','carbon_mod.F90')
   ENDIF

   ALLOCATE( BCCONV(State_Grid%NX,State_Grid%NY,State_Grid%NZ), STAT=AS )
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCCONV' )
   BCCONV = 0e+0_fp

   ALLOCATE( OCCONV(State_Grid%NX,State_Grid%NY,State_Grid%NZ), STAT=AS )
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCCONV' )
   OCCONV = 0e+0_fp

   ! semivolpoa2: for POA emissions (hotp 2/27/09)
   ! Store POG1 and POG2 separately (mps, 1/14/16)
   ALLOCATE( POAEMISS(State_Grid%NX,State_Grid%NY,State_Grid%NZ,2),STAT=AS )
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'POAEMISS' )
   POAEMISS = 0e+0_fp

#ifdef APM
   ALLOCATE( BCCONVNEW(State_Grid%NX,State_Grid%NY,State_Grid%NZ), STAT=AS )
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCCONVNEW' )
   BCCONVNEW = 0e+0_fp
   ALLOCATE( OCCONVNEW(State_Grid%NX,State_Grid%NY,State_Grid%NZ), STAT=AS )
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCCONVNEW' )
   OCCONVNEW = 0e+0_fp
#endif

   ! JKODROS
#ifdef TOMAS
   !SFARINA the next six are introduced with the idea that
   ! emisscarbontomas needs to be restructured and these
   ! data structures eliminated
   ALLOCATE( BCFF(State_Grid%NX,State_Grid%NY,IBINS,2), STAT=AS)
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCFF' )
   BCFF = 0e+0_fp

   ALLOCATE( OCFF(State_Grid%NX,State_Grid%NY,IBINS,2), STAT=AS)
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCFF' )
   OCFF = 0e+0_fp

   ALLOCATE( BCBF(State_Grid%NX,State_Grid%NY,IBINS,2), STAT=AS)
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCBF' )
   BCBF = 0e+0_fp

   ALLOCATE( OCBF(State_Grid%NX,State_Grid%NY,IBINS,2), STAT=AS)
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCBF' )
   OCBF = 0e+0_fp

   ALLOCATE( BCBB(State_Grid%NX,State_Grid%NY,IBINS,2), STAT=AS)
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCBB' )
   BCBB = 0e+0_fp

   ALLOCATE( OCBB(State_Grid%NX,State_Grid%NY,IBINS,2), STAT=AS)
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCBB' )
   OCBB = 0e+0_fp

   ! BC
   ALLOCATE( BCPI_ANTH_BULK(State_Grid%NX,State_Grid%NY), STAT=AS)
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCPI_ANTH_BULK' )
   BCPI_ANTH_BULK = 0e+0_fp

   ALLOCATE( BCPO_ANTH_BULK(State_Grid%NX,State_Grid%NY), STAT=AS)
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCPO_ANTH_BULK' )
   BCPO_ANTH_BULK = 0e+0_fp

   ALLOCATE( BCPI_BIOB_BULK(State_Grid%NX,State_Grid%NY), STAT=AS)
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCPI_BIOB_BULK' )
   BCPI_BIOB_BULK = 0e+0_fp

   ALLOCATE( BCPO_BIOB_BULK(State_Grid%NX,State_Grid%NY), STAT=AS)
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCPO_BIOB_BULK' )
   BCPO_BIOB_BULK = 0e+0_fp

   ! OC ----------------
   ALLOCATE( OCPI_ANTH_BULK(State_Grid%NX,State_Grid%NY), STAT=AS)
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCPI_ANTH_BULK' )
   OCPI_ANTH_BULK = 0e+0_fp

   ALLOCATE( OCPO_ANTH_BULK(State_Grid%NX,State_Grid%NY), STAT=AS)
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCPO_ANTH_BULK' )
   OCPO_ANTH_BULK = 0e+0_fp

   ALLOCATE( OCPI_BIOB_BULK(State_Grid%NX,State_Grid%NY), STAT=AS)
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCPI_BIOB_BULK' )
   OCPI_BIOB_BULK = 0e+0_fp

   ALLOCATE( OCPO_BIOB_BULK(State_Grid%NX,State_Grid%NY), STAT=AS)
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCPO_BIOB_BULK' )
   OCPO_BIOB_BULK = 0e+0_fp

   !biogenic
   ALLOCATE( TERP_ORGC(State_Grid%NX,State_Grid%NY), STAT=AS )
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'TERP_ORGC' )
   TERP_ORGC = 0e+0_fp

   ! CO anth for scaling xtraSOA
   ALLOCATE( CO_ANTH(State_Grid%NX,State_Grid%NY), STAT=AS )
   IF ( AS /= 0 ) CALL ALLOC_ERR( 'CO_ANTH' )
   CO_ANTH = 0e+0_fp
#endif

   !=================================================================
   ! SOA arrays only have to be allocated if LSOA = T
   !=================================================================
   IF ( Input_Opt%LSOA ) THEN

      ALLOCATE( TCOSZ(State_Grid%NX,State_Grid%NY), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCOSZ' )
      TCOSZ = 0e+0_fp

      ALLOCATE( ORVC_SESQ(State_Grid%NX,State_Grid%NY, State_Grid%NZ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ORVC_SESQ' )
      ORVC_SESQ = 0e+0_fp

      ! diagnostic  (dkh, 11/11/06)
      ! increase last dimension by 1 to add NAP (hotp 7/22/09)
      ALLOCATE( GLOB_DARO2(State_Grid%NX,State_Grid%NY,State_Grid%NZ,2,4), &
                STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GLOB_DARO2' )
      GLOB_DARO2 = 0e+0_fp

      ! semivolpoa4: diagnostic (hotp 3/27/09)
      ALLOCATE( GLOB_POGRXN(State_Grid%NX,State_Grid%NY,State_Grid%NZ,2 ), &
                STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GLOB_POGRXN' )
      GLOB_POGRXN = 0e+0_fp

      ! Initial OA+OG diagnostic (hotp 5/17/10)
      ALLOCATE( OAGINITSAVE(State_Grid%NX,State_Grid%NY,State_Grid%NZ,MPROD,MSV), &
                STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OAGINITSAVE' )
      OAGINITSAVE = 0e+0_fp

      ! Change in OA+OG diagnostic (hotp 5/17/10)
      ALLOCATE( DELTASOGSAVE(State_Grid%NX,State_Grid%NY,State_Grid%NZ,MNOX,MHC), &
                STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DELTASOGSAVE' )
      DELTASOGSAVE = 0e+0_fp

      ! Diagnostic for NO branching ratio (hotp 5/24/10)
      ALLOCATE( BETANOSAVE(State_Grid%NX,State_Grid%NY,State_Grid%NZ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BETANOSAVE' )
      BETANOSAVE = 0e+0_fp

      ! Diagnostic (hotp 6/5/10)
      ALLOCATE( SPECSOAPROD(State_Grid%NX,State_Grid%NY,State_Grid%NZ,MPROD,MSV), &
                STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SPECSOAPROD' )
      SPECSOAPROD = 0e+0_fp

      ! Diagnostic (hotp 6/5/10)
      ALLOCATE( SPECSOAEVAP(State_Grid%NX,State_Grid%NY,State_Grid%NZ,MPROD,MSV), &
                STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SPECSOAEVAP' )
      SPECSOAEVAP = 0e+0_fp

   ENDIF

 END SUBROUTINE INIT_CARBON
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_carbon
!
! !DESCRIPTION: Subroutine CLEANUP\_CARBON deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE CLEANUP_CARBON
!
! !REVISION HISTORY:
!  01 Apr 2004 - R. Park - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

   !=================================================================
   ! CLEANUP_CARBON begins here!
   !=================================================================
   IF ( ALLOCATED( BCCONV        ) ) DEALLOCATE( BCCONV        )
   IF ( ALLOCATED( OCCONV        ) ) DEALLOCATE( OCCONV        )
   IF ( ALLOCATED( TCOSZ         ) ) DEALLOCATE( TCOSZ         )
   IF ( ALLOCATED( ORVC_SESQ     ) ) DEALLOCATE( ORVC_SESQ     )
   IF ( ALLOCATED( GLOB_DARO2    ) ) DEALLOCATE( GLOB_DARO2    )
   IF ( ALLOCATED( POAEMISS      ) ) DEALLOCATE( POAEMISS      )
   IF ( ALLOCATED( GLOB_POGRXN   ) ) DEALLOCATE( GLOB_POGRXN   )
   IF ( ALLOCATED( OAGINITSAVE   ) ) DEALLOCATE( OAGINITSAVE   )
   IF ( ALLOCATED( DELTASOGSAVE  ) ) DEALLOCATE( DELTASOGSAVE  )
   IF ( ALLOCATED( BETANOSAVE    ) ) DEALLOCATE( BETANOSAVE    )
   IF ( ALLOCATED( SPECSOAPROD   ) ) DEALLOCATE( SPECSOAPROD   )
   IF ( ALLOCATED( SPECSOAEVAP   ) ) DEALLOCATE( SPECSOAEVAP   )
#ifdef APM
   IF ( ALLOCATED( BCCONVNEW     ) ) DEALLOCATE( BCCONVNEW     )
   IF ( ALLOCATED( OCCONVNEW     ) ) DEALLOCATE( OCCONVNEW     )
#endif
#ifdef TOMAS
   IF ( ALLOCATED( BCFF           )) DEALLOCATE( BCFF          )
   IF ( ALLOCATED( OCFF           )) DEALLOCATE( OCFF          )
   IF ( ALLOCATED( BCBF           )) DEALLOCATE( BCBF          )
   IF ( ALLOCATED( OCBF           )) DEALLOCATE( OCBF          )
   IF ( ALLOCATED( BCBB           )) DEALLOCATE( BCBB          )
   IF ( ALLOCATED( OCBB           )) DEALLOCATE( OCBB          )
   IF ( ALLOCATED( BCPI_ANTH_BULK )) DEALLOCATE( BCPI_ANTH_BULK)
   IF ( ALLOCATED( BCPO_ANTH_BULK )) DEALLOCATE( BCPO_ANTH_BULK)
   IF ( ALLOCATED( BCPI_BIOB_BULK )) DEALLOCATE( BCPI_BIOB_BULK)
   IF ( ALLOCATED( BCPO_BIOB_BULK )) DEALLOCATE( BCPO_BIOB_BULK)
   IF ( ALLOCATED( OCPI_ANTH_BULK )) DEALLOCATE( OCPI_ANTH_BULK)
   IF ( ALLOCATED( OCPO_ANTH_BULK )) DEALLOCATE( OCPO_ANTH_BULK)
   IF ( ALLOCATED( OCPI_BIOB_BULK )) DEALLOCATE( OCPI_BIOB_BULK)
   IF ( ALLOCATED( OCPO_BIOB_BULK )) DEALLOCATE( OCPO_BIOB_BULK)
   IF ( ALLOCATED( TERP_ORGC      )) DEALLOCATE( TERP_ORGC     )
   IF ( ALLOCATED( CO_ANTH        )) DEALLOCATE( CO_ANTH       )
   FIRE_NUM                       => NULL() !(ramnarine 12/27/2018)
#endif

 END SUBROUTINE CLEANUP_CARBON
!EOC
#ifdef APM
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_bcponew
!
! !DESCRIPTION: Subroutine CHEM\_BCPONEW converts hydrophobic BC to hydrophilic
!  BC and calculates the dry deposition of hydrophobic BC. Modified for
!  APM simulation.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE CHEM_BCPONEW( Input_Opt, State_Grid, TC, RC )
!
! !USES:
!
   USE ErrCode_Mod
   USE Input_Opt_Mod,  ONLY : OptInput
   USE State_Grid_Mod, ONLY : GrdState
   USE APM_INIT_MOD,   ONLY : APMIDS
   USE TIME_MOD,       ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt          ! Input Options
   TYPE(GrdState), INTENT(IN)    :: State_grid         ! Grid State
!
! !INPUT/OUTPUT PARAMETERS:
!
   ! H-phobic BC [kg]
   REAL(fp),       INTENT(INOUT) :: TC(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC                 ! Success?
!
! !REVISION HISTORY:
!  16 Feb 2011 - R. Yantosca - Initial version, from G. Luo
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Scalars
   LOGICAL             :: LNLPBL
   INTEGER             :: I,       J,   L,   N_TRACERS
   REAL(fp)            :: DTCHEM, KBC,  FREQ
   REAL(fp)            :: TC0,    CNEW, RKT, BL_FRAC
!
! !DEFINED PARAMETERS:
!
   REAL(fp), PARAMETER :: BC_LIFE = 1.15D0

   !=================================================================
   ! CHEM_BCPONEW begins here!
   !=================================================================

   ! Return if BCPO isn't defined
   IF ( id_BCPO == 0 ) RETURN

   ! Assume success
   RC        = GC_SUCCESS

   ! Initialize
   KBC       = 1.0e+0_fp / ( 86400e+0_fp * BC_LIFE )
   DTCHEM    = GET_TS_CHEM()

   ! Zero BCPO -> BCPI conversion array
   BCCONVNEW  = 0e+0_fp

   !$OMP PARALLEL DO       &
   !$OMP DEFAULT( SHARED ) &
   !$OMP PRIVATE( I, J, L, TC0, FREQ, BL_FRAC, RKT, CNEW ) &
   !$OMP SCHEDULE( DYNAMIC )
   DO L = 1, State_Grid%NZ
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      ! Initial BC mass [kg]
      TC0  = TC(I,J,L)

      ! Zero drydep freq
      FREQ = 0e+0_fp

      ! Amount of BCPO left after chemistry and drydep [kg]
      RKT  = ( KBC + FREQ ) * DTCHEM
      CNEW = TC0 * EXP( -RKT )

      ! Prevent underflow condition
      IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

      ! Amount of BCPO converted to BCPI [kg/timestep]
      BCCONVNEW(I,J,L) = ( TC0 - CNEW ) * KBC / ( KBC + FREQ )

      ! Store new concentration back into tracer array
      TC(I,J,L) = CNEW
   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

 END SUBROUTINE CHEM_BCPONEW
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_bcpinew
!
! !DESCRIPTION: Subroutine CHEM\_BCPINEW calculates dry deposition of
!  hydrophilic BC. Modified for APM simulation.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE CHEM_BCPINEW( Input_Opt, State_Grid, TC, RC )
!
! !USES:
!
   USE ErrCode_Mod
   USE Input_Opt_Mod,  ONLY : OptInput
   USE State_Grid_Mod, ONLY : GrdState
   USE APM_INIT_MOD,   ONLY : APMIDS
   USE TIME_MOD,       ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt         ! Input Options
   TYPE(GrdState), INTENT(IN)    :: State_Grid        ! Grid State
!
! !INPUT/OUTPUT PARAMETERS:
!
   ! H-philic BC [kg]
   REAL(fp),       INTENT(INOUT) :: TC(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC               ! Success?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  16 Feb 2011 - R. Yantosca - Initial version, from G. Luo
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Scalars
   LOGICAL  :: LNLPBL
   INTEGER  :: I,      J,      L,       N_TRACERS
   REAL(fp) :: DTCHEM, BL_FRAC
   REAL(fp) :: TC0,    CNEW,   CCV,     FREQ

   !=================================================================
   ! CHEM_BCPINEW begins here!
   !=================================================================

   ! Return if BCPI isn't defined
   IF ( id_BCPI == 0 ) RETURN

   ! Assume success
   RC     = GC_SUCCESS

   ! Chemistry timestep [s]
   DTCHEM = GET_TS_CHEM()

   !$OMP PARALLEL DO       &
   !$OMP DEFAULT( SHARED ) &
   !$OMP PRIVATE( I, J, L, TC0, CCV, FREQ, BL_FRAC, CNEW ) &
   !$OMP SCHEDULE( DYNAMIC )
   DO L = 1, State_Grid%NZ
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      ! Initial H-philic BC [kg]
      TC0 = TC(I,J,L)

      ! H-philic BC that used to be H-phobic BC [kg]
      CCV = BCCONVNEW(I,J,L)

      ! Otherwise, omit the exponential to save on clock cycles
      CNEW = TC0 + CCV

      ! Prevent underflow condition
      IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

      ! Save new concentration of H-philic IC in tracer array
      TC(I,J,L) = CNEW

   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Zero BCPO -> BCPI conversion array
   BCCONVNEW = 0e+0_fp

 END SUBROUTINE CHEM_BCPINEW
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_ocponew
!
! !DESCRIPTION: Subroutine CHEM\_OCPONEW converts hydrophobic OC to hydrophilic
!  OC and calculates the dry deposition of hydrophobic OC. Modified for APM
!  simulation.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE CHEM_OCPONEW( Input_Opt, State_Grid, TC, RC )
!
! !USES:
!
   USE ErrCode_Mod
   USE Input_Opt_Mod,  ONLY : OptInput
   USE State_Grid_Mod, ONLY : GrdState
   USE APM_INIT_MOD,   ONLY : APMIDS
   USE TIME_MOD,       ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt          ! Input Options
   TYPE(GrdState), INTENT(IN)    :: State_Grid         ! Grid State
!
! !INPUT/OUTPUT PARAMETERS:
!
   ! H-phobic OC [kg]
   REAL(fp),       INTENT(INOUT) :: TC(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC                 ! Success?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  16 Feb 2011 - R. Yantosca - Initial version, from G. Luo
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Scalars
   LOGICAL             :: LNLPBL
   INTEGER             :: I,      J,    L,      N_TRACERS
   REAL(fp)            :: DTCHEM, KOC,  BL_FRAC
   REAL(fp)            :: TC0,    FREQ, CNEW,   RKT
!
! !DEFINED PARAMETERS:
!
   REAL(fp), PARAMETER :: OC_LIFE = 1.15e+0_fp

   !=================================================================
   ! CHEM_OCPONEW begins here!
   !=================================================================

   ! Return if OCPO isn't defined
   IF ( MAX(id_OCPO,id_POA1) == 0 ) RETURN

   ! Assume success
   RC        = GC_SUCCESS

   ! Initialize
   KOC       = 1.0e+0_fp / ( 86400e+0_fp * OC_LIFE )
   DTCHEM    = GET_TS_CHEM()

   ! Zero OCPO -> OCPI conversion array
   OCCONVNEW = 0e+0_fp

   !$OMP PARALLEL DO       &
   !$OMP DEFAULT( SHARED ) &
   !$OMP PRIVATE( I, J, L, TC0, FREQ, BL_FRAC, RKT, CNEW ) &
   !$OMP SCHEDULE( DYNAMIC )
   DO L = 1, State_Grid%NZ
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      ! Initial OC [kg]
      TC0  = TC(I,J,L)

      ! Zero drydep freq
      FREQ = 0e+0_fp

      ! Amount of OCPO left after chemistry and drydep [kg]
      RKT  = ( KOC + FREQ ) * DTCHEM
      CNEW = TC0 * EXP( -RKT )

      ! Prevent underflow condition
      IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

      ! Amount of OCPO converted to OCPI [kg/timestep]
      OCCONVNEW(I,J,L) = ( TC0 - CNEW ) * KOC / ( KOC + FREQ )

      ! Store modified OC concentration back in tracer array
      TC(I,J,L) = CNEW

   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

 END SUBROUTINE CHEM_OCPONEW
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_ocpinew
!
! !DESCRIPTION: Subroutine CHEM\_OCPINEW calculates dry deposition of
!  hydrophilic OC. Modified for APM simulation.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE CHEM_OCPINEW( Input_Opt, State_Grid, TC, RC )
!
! !USES:
!
   USE ErrCode_Mod
   USE Input_Opt_Mod,  ONLY : OptInput
   USE State_Grid_Mod, ONLY : GrdState
   USE APM_INIT_MOD,   ONLY : APMIDS
   USE TIME_MOD,       ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt          ! Input Options
   TYPE(GrdState), INTENT(IN)    :: State_Grid         ! Grid State
!
! !INPUT/OUTPUT PARAMETERS:
!
   ! H-philic OC [kg]
   REAL(fp),       INTENT(INOUT) :: TC(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC                 ! Success?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  16 Feb 2011 - R. Yantosca - Initial version, from G. Luo
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Scalars
   LOGICAL  :: LNLPBL
   INTEGER  :: I,      J,      L,   N_TRACERS
   REAL(fp) :: DTCHEM, BL_FRAC
   REAL(fp) :: TC0,    CNEW,   CCV, FREQ

   !=================================================================
   ! CHEM_OCPINEW begins here!
   !=================================================================

   ! Return if OCPI isn't defined
   IF ( MAX(id_OCPI,id_POA1) == 0 ) RETURN

   ! Assume success
   RC        = GC_SUCCESS

   ! Chemistry timestep [s]
   DTCHEM = GET_TS_CHEM()

   !$OMP PARALLEL DO       &
   !$OMP DEFAULT( SHARED ) &
   !$OMP PRIVATE( I, J, L, TC0, CCV, BL_FRAC, FREQ, CNEW ) &
   !$OMP SCHEDULE( DYNAMIC )
   DO L = 1, State_Grid%NZ
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      ! Initial H-philic OC [kg]
      TC0 = TC(I,J,L)

      ! H-philic OC that used to be H-phobic OC [kg]
      CCV = OCCONVNEW(I,J,L)

      CNEW = TC0 + CCV

      ! Prevent underflow condition
      IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

      ! Store modified concentration back in tracer array [kg]
      TC(I,J,L) = CNEW

   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Zero OCPO -> OCPI conversion array
   OCCONVNEW = 0e+0_fp

 END SUBROUTINE CHEM_OCPINEW
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dry_settlingbin
!
! !DESCRIPTION: Subroutine DRY\_SETTLINGBIN computes the dry settling of
!  aerosol tracers. Modified for APM simulation. (G. Luo)
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE BCDRY_SETTLINGBIN( Input_Opt,  State_Chm, State_Diag, &
                               State_Grid, State_Met, RC )
!
! !USES:
!
   USE ErrCode_Mod
   USE Input_Opt_Mod,  ONLY : OptInput
   USE State_Chm_Mod,  ONLY : ChmState
   USE State_Grid_Mod, ONLY : GrdState
   USE State_Diag_Mod, ONLY : DgnState
   USE State_Met_Mod,  ONLY : MetState
   USE PRESSURE_MOD,   ONLY : GET_PCENTER
   USE TIME_MOD,       ONLY : GET_TS_CHEM
   USE APM_INIT_MOD,   ONLY : APMIDS
   USE APM_INIT_MOD,   ONLY : NCTBC,NBCOC
   USE APM_INIT_MOD,   ONLY : RBCOC, DENBC
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
   TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
   TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
   TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  16 Feb 2011 - R. Yantosca - Initial version, from G. Luo
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Local variables
   INTEGER           :: I, J, L, N, K
   INTEGER           :: IDTEMP
   REAL(fp)          :: DT_SETTL, DELZ,  DELZ1
   REAL(fp)          :: REFF,     DEN,   CONST
   REAL(fp)          :: NUM,      LAMDA, FLUX
   REAL(fp)          :: AREA_CM2, TC0(State_Grid%NZ)
   REAL(fp)          :: TOT1,     TOT2

   ! Pressure in Kpa 1 mb = 100 pa = 0.1 kPa
   REAL(fp)          :: P

   ! Diameter of aerosol [um]
   REAL(fp)          :: Dp
   
   ! Pressure * DP
   REAL(fp)          :: PDp

   ! Temperature (K)
   REAL(fp)          :: TEMP

   ! Slip correction factor
   REAL(fp)          :: Slip

   ! Viscosity of air (Pa s)
   REAL(fp)          :: Visc

   ! Settling velocity of particle (m/s)
   REAL(fp)          :: VTS(State_Grid%NZ)
   REAL(fp)          :: MASS(State_Grid%NZ)
   REAL(fp)          :: OLD(State_Grid%NZ,NCTBC)

   ! Make a pointer to the tracer array
   REAL(fp), POINTER :: Spc(:,:,:,:)

   !=================================================================
   ! DRY_SETTLINGBIN begins here!
   !=================================================================

   ! Assume success
   RC        = GC_SUCCESS

   ! Point to Spc
   Spc => State_Chm%species

   ! Aerosol settling timestep [s]
   DT_SETTL = GET_TS_CHEM()

   IDTEMP = APMIDS%id_BCBIN1+NBCOC-1

   !$OMP PARALLEL DO       &
   !$OMP DEFAULT( SHARED ) &
   !$OMP PRIVATE( I, J, L, N, K, DEN, REFF, DP )       &
   !$OMP PRIVATE( CONST, VTS, TEMP, P, PDP, SLIP )     &
   !$OMP PRIVATE( MASS, OLD, VISC, TC0, DELZ, DELZ1  ) &
   !$OMP SCHEDULE( DYNAMIC )
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      DO L = 1, State_Grid%NZ
         MASS(L) = SUM(Spc(I,J,L,APMIDS%id_BCBIN1:IDTEMP))
         DO K = 1, NCTBC
            OLD(L,K) = Spc(I,J,L,(APMIDS%id_CTBC+K-1))
            Spc(I,J,L,(APMIDS%id_CTBC+K-1)) = 0.e+0_fp
         ENDDO
      ENDDO

      ! Loop over aerosol bins
      DO N = 1, NBCOC

         DO L = 1, State_Grid%NZ

            TC0(L) = Spc(I,J,L,(APMIDS%id_BCBIN1+N-1))

            IF(TC0(L)>1.e-30_fp)THEN
               ! Initialize
               DEN   = DENBC
               REFF  = RBCOC(N)
               DP    = 2e+0_fp * REFF * 1.e+6_fp ! Dp [um] = particle diameter
               CONST = 2e+0_fp * DEN * REFF**2 * G0 / 9e+0_fp

               ! Get P [kPa], T [K], and P*DP
               P    = GET_PCENTER(I,J,L) * 0.1e+0_fp
               TEMP = State_Met%T(I,J,L)
               PDP  = P * DP

               ! Slip correction factor as function of (P*dp)
               SLIP = 1e+0_fp + ( 15.60e0 + 7.0e0 * EXP(-0.059e0*PDP) ) / PDP

               ! Viscosity [Pa s] of air as a function of temp (K)
               VISC = 1.458e-6 * (TEMP)**(1.5e0) / ( TEMP + 110.4e0 )

               ! Settling velocity [m/s]
               VTS(L) = CONST * SLIP / VISC
            ELSE
               VTS(L) = 0.e+0_fp
            ENDIF

         ENDDO

         ! Method is to solve bidiagonal matrix
         ! which is implicit and first order accurate in Z
         L    = State_Grid%NZ
         IF(MASS(L)>1.e-30)THEN
            DELZ = State_Met%BXHEIGHT(I,J,L)

            Spc(I,J,L,(APMIDS%id_BCBIN1+N-1)) = &
               Spc(I,J,L,(APMIDS%id_BCBIN1+N-1)) / &
               ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )

            DO K = 1, NCTBC
               Spc(I,J,L,(APMIDS%id_CTBC+K-1)) = &
                  Spc(I,J,L,(APMIDS%id_CTBC+K-1))+ &
                  OLD(L,K)*TC0(L)/MASS(L) / &
                  ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )
            ENDDO
         ENDIF

         DO L = State_Grid%NZ-1, 1, -1
            IF((MASS(L)*MASS(L+1))>1.e-30)THEN
               DELZ  = State_Met%BXHEIGHT(I,J,L)
               DELZ1 = State_Met%BXHEIGHT(I,J,L+1)
               Spc(I,J,L,(APMIDS%id_BCBIN1+N-1)) = 1.e+0_fp / &
                  ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ ) &
                  * (Spc(I,J,L,(APMIDS%id_BCBIN1+N-1)) &
                  + DT_SETTL * VTS(L+1) / DELZ1 &
                  * TC0(L+1) )

               DO K = 1, NCTBC
                  Spc(I,J,L,(APMIDS%id_CTBC+K-1)) = &
                     Spc(I,J,L,(APMIDS%id_CTBC+K-1))+ &
                     1.e+0_fp / &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ ) &
                     * (OLD(L,K)*TC0(L)/MASS(L) &
                     + DT_SETTL * VTS(L+1) / DELZ1 &
                     * OLD(L+1,K)*TC0(L+1)/MASS(L+1) )
               ENDDO

            ELSE IF(MASS(L)>1.e-30)THEN
               DELZ  = State_Met%BXHEIGHT(I,J,L)
               DELZ1 = State_Met%BXHEIGHT(I,J,L+1)
               Spc(I,J,L,(APMIDS%id_BCBIN1+N-1)) = 1.e+0_fp / &
                  ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ ) &
                  * (Spc(I,J,L,(APMIDS%id_BCBIN1+N-1)) &
                  + DT_SETTL * VTS(L+1) / DELZ1 &
                  * TC0(L+1) )

               DO K = 1, NCTBC
                  Spc(I,J,L,(APMIDS%id_CTBC+K-1)) = &
                     Spc(I,J,L,(APMIDS%id_CTBC+K-1))+ &
                     1.e+0_fp / &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ ) &
                     * OLD(L,K)*TC0(L)/MASS(L)
               ENDDO

            ELSE IF(MASS(L+1)>1.e-30)THEN
               DELZ  = State_Met%BXHEIGHT(I,J,L)
               DELZ1 = State_Met%BXHEIGHT(I,J,L+1)
               Spc(I,J,L,(APMIDS%id_BCBIN1+N-1)) = 1.e+0_fp / &
                  ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ ) &
                  * (Spc(I,J,L,(APMIDS%id_BCBIN1+N-1)) &
                  + DT_SETTL * VTS(L+1) / DELZ1 &
                  * TC0(L+1) )

               DO K = 1, NCTBC
                  Spc(I,J,L,(APMIDS%id_CTBC+K-1)) = &
                     Spc(I,J,L,(APMIDS%id_CTBC+K-1))+ &
                     1.e+0_fp / &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ ) &
                     * DT_SETTL * VTS(L+1) / DELZ1 &
                     * OLD(L+1,K)*TC0(L+1)/MASS(L+1)
               ENDDO
            ENDIF

         ENDDO

      ENDDO

      DO L = 1, State_Grid%NZ
         DO K = 1, NCTBC
            Spc(I,J,L,(APMIDS%id_CTBC+K-1)) = &
               MAX(1.d-30,Spc(I,J,L,(APMIDS%id_CTBC+K-1)))
         ENDDO
      ENDDO

   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Clear the pointer
   NULLIFY( Spc )

 END SUBROUTINE BCDRY_SETTLINGBIN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dry_settlingbin
!
! !DESCRIPTION: Subroutine DRY\_SETTLINGBIN computes the dry settling of
!  aerosol tracers. Modified for APM simulation. (G. Luo)
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE OCDRY_SETTLINGBIN( Input_Opt,  State_Chm, State_Diag, &
                               State_Grid, State_Met, RC )
!
! !USES:
!
   USE ErrCode_Mod
   USE Input_Opt_Mod,  ONLY : OptInput
   USE State_Chm_Mod,  ONLY : ChmState
   USE State_Grid_Mod, ONLY : GrdState
   USE State_Diag_Mod, ONLY : DgnState
   USE State_Met_Mod,  ONLY : MetState
   USE PRESSURE_MOD,   ONLY : GET_PCENTER
   USE TIME_MOD,       ONLY : GET_TS_CHEM
   USE APM_INIT_MOD,   ONLY : APMIDS
   USE APM_INIT_MOD,   ONLY : NCTOC,NBCOC
   USE APM_INIT_MOD,   ONLY : RBCOC, DENOC
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
   TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
   TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
   TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  16 Feb 2011 - R. Yantosca - Initial version, from G. Luo
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Local variables
   INTEGER           :: I, J, L, N, K
   INTEGER           :: IDTEMP
   REAL(fp)          :: DT_SETTL, DELZ,  DELZ1
   REAL(fp)          :: REFF,     DEN,   CONST
   REAL(fp)          :: NUM,      LAMDA, FLUX
   REAL(fp)          :: AREA_CM2, TC0(State_Grid%NZ)
   REAL(fp)          :: TOT1,     TOT2

   ! Pressure in Kpa 1 mb = 100 pa = 0.1 kPa
   REAL(fp)          :: P

   ! Diameter of aerosol [um]
   REAL(fp)          :: Dp

   ! Pressure * DP
   REAL(fp)          :: PDp

   ! Temperature (K)
   REAL(fp)          :: TEMP

   ! Slip correction factor
   REAL(fp)          :: Slip

   ! Viscosity of air (Pa s)
   REAL(fp)          :: Visc

   ! Settling velocity of particle (m/s)
   REAL(fp)          :: VTS(State_Grid%NZ)
   REAL(fp)          :: MASS(State_Grid%NZ)
   REAL(fp)          :: OLD(State_Grid%NZ,NCTOC)

   ! Make a pointer to the tracer array
   REAL(fp), POINTER :: Spc(:,:,:,:)

   !=================================================================
   ! DRY_SETTLINGBIN begins here!
   !=================================================================

   ! Assume success
   RC        = GC_SUCCESS

   ! Point to Spc
   Spc => State_Chm%species

   ! Aerosol settling timestep [s]
   DT_SETTL = GET_TS_CHEM()

   IDTEMP = APMIDS%id_OCBIN1+NBCOC-1

   !$OMP PARALLEL DO       &
   !$OMP DEFAULT( SHARED ) &
   !$OMP PRIVATE( I, J, L, N, K, DEN, REFF, DP )       &
   !$OMP PRIVATE( CONST, VTS, TEMP, P, PDP, SLIP )     &
   !$OMP PRIVATE( MASS, OLD, VISC, TC0, DELZ, DELZ1  ) &
   !$OMP SCHEDULE( DYNAMIC )
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      DO L = 1, State_Grid%NZ
         MASS(L) = SUM(Spc(I,J,L,APMIDS%id_OCBIN1:IDTEMP))
         DO K = 1, NCTOC
            OLD(L,K) = Spc(I,J,L,(APMIDS%id_CTOC+K-1))
            Spc(I,J,L,(APMIDS%id_CTOC+K-1)) = 0.e+0_fp
         ENDDO
      ENDDO

      ! Loop over aerosol bins
      DO N = 1, NBCOC

         DO L = 1, State_Grid%NZ

            TC0(L) = Spc(I,J,L,(APMIDS%id_OCBIN1+N-1))

            IF(TC0(L)>1.e-30)THEN
               ! Initialize
               DEN   = DENOC
               REFF  = RBCOC(N)
               DP    = 2e+0_fp * REFF * 1.e+6_fp ! Dp [um] = particle diameter
               CONST = 2e+0_fp * DEN * REFF**2 * G0 / 9e+0_fp

               ! Get P [kPa], T [K], and P*DP
               P    = GET_PCENTER(I,J,L) * 0.1e+0_fp
               TEMP = State_Met%T(I,J,L)
               PDP  = P * DP

               ! Slip correction factor as function of (P*dp)
               SLIP = 1e+0_fp + ( 15.60e0 + 7.0e0 * EXP(-0.059e0*PDP) ) / PDP

               ! Viscosity [Pa s] of air as a function of temp (K)
               VISC = 1.458e-6 * (TEMP)**(1.5e0) / ( TEMP + 110.4e0 )

               ! Settling velocity [m/s]
               VTS(L) = CONST * SLIP / VISC
            ELSE
               VTS(L) = 0.e+0_fp
            ENDIF

         ENDDO

         ! Method is to solve bidiagonal matrix
         ! which is implicit and first order accurate in Z
         L    = State_Grid%NZ
         IF(MASS(L)>1.e-30)THEN
            DELZ = State_Met%BXHEIGHT(I,J,L)

            Spc(I,J,L,(APMIDS%id_OCBIN1+N-1)) = &
               Spc(I,J,L,(APMIDS%id_OCBIN1+N-1)) / &
               ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )

            DO K = 1, NCTOC
               Spc(I,J,L,(APMIDS%id_CTOC+K-1)) = &
                  Spc(I,J,L,(APMIDS%id_CTOC+K-1))+ &
                  OLD(L,K)*TC0(L)/MASS(L) / &
                  ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )
            ENDDO
         ENDIF

         DO L = State_Grid%NZ-1, 1, -1
            IF((MASS(L)*MASS(L+1))>1.e-30)THEN
               DELZ  = State_Met%BXHEIGHT(I,J,L)
               DELZ1 = State_Met%BXHEIGHT(I,J,L+1)
               Spc(I,J,L,(APMIDS%id_OCBIN1+N-1)) = 1.e+0_fp / &
                  ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ ) &
                  * (Spc(I,J,L,(APMIDS%id_OCBIN1+N-1)) &
                  + DT_SETTL * VTS(L+1) / DELZ1 &
                  * TC0(L+1) )

               DO K = 1, NCTOC
                  Spc(I,J,L,(APMIDS%id_CTOC+K-1)) = &
                     Spc(I,J,L,(APMIDS%id_CTOC+K-1))+ &
                     1.e+0_fp / &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ ) &
                     * (OLD(L,K)*TC0(L)/MASS(L) &
                     + DT_SETTL * VTS(L+1) / DELZ1 &
                     * OLD(L+1,K)*TC0(L+1)/MASS(L+1) )
               ENDDO

            ELSE IF(MASS(L)>1.e-30)THEN
               DELZ  = State_Met%BXHEIGHT(I,J,L)
               DELZ1 = State_Met%BXHEIGHT(I,J,L+1)
               Spc(I,J,L,(APMIDS%id_OCBIN1+N-1)) = 1.e+0_fp / &
                  ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ ) &
                  * (Spc(I,J,L,(APMIDS%id_OCBIN1+N-1)) &
                  + DT_SETTL * VTS(L+1) / DELZ1 &
                  * TC0(L+1) )

               DO K = 1, NCTOC
                  Spc(I,J,L,(APMIDS%id_CTOC+K-1)) = &
                     Spc(I,J,L,(APMIDS%id_CTOC+K-1))+ &
                     1.e+0_fp / &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ ) &
                     * OLD(L,K)*TC0(L)/MASS(L)
               ENDDO

            ELSE IF(MASS(L+1)>1.e-30)THEN
               DELZ  = State_Met%BXHEIGHT(I,J,L)
               DELZ1 = State_Met%BXHEIGHT(I,J,L+1)
               Spc(I,J,L,(APMIDS%id_OCBIN1+N-1)) = 1.e+0_fp / &
                  ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ ) &
                  * (Spc(I,J,L,(APMIDS%id_OCBIN1+N-1)) &
                  + DT_SETTL * VTS(L+1) / DELZ1 &
                  * TC0(L+1) )

               DO K = 1, NCTOC
                  Spc(I,J,L,(APMIDS%id_CTOC+K-1)) = &
                     Spc(I,J,L,(APMIDS%id_CTOC+K-1))+ &
                     1.e+0_fp / &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ ) &
                     * DT_SETTL * VTS(L+1) / DELZ1 &
                     * OLD(L+1,K)*TC0(L+1)/MASS(L+1)
               ENDDO
            ENDIF

         ENDDO

      ENDDO

      DO L = 1, State_Grid%NZ
         DO K = 1, NCTOC
            Spc(I,J,L,(APMIDS%id_CTOC+K-1)) = &
               MAX(1.d-30,Spc(I,J,L,(APMIDS%id_CTOC+K-1)))
         ENDDO
      ENDDO

   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Clear the pointer
   NULLIFY( Spc )

 END SUBROUTINE OCDRY_SETTLINGBIN
#endif
!EOC
      END MODULE CARBON_MOD
