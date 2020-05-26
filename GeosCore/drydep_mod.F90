!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: drydep_mod.F90
!
! !DESCRIPTION: Module DRYDEP\_MOD contains variables and routines for the
!  GEOS-Chem dry deposition scheme.
!\\
!\\
! !INTERFACE:
!
MODULE DRYDEP_MOD
!
! !USES:
!
  USE CMN_SIZE_Mod,     ONLY : NPOLY, NSURFTYPE
  USE ERROR_MOD              ! Error handling routines
#ifdef TOMAS                 
  USE TOMAS_MOD              ! For TOMAS microphysics
#endif                       
  USE PhysConstants          ! Physical constants
  USE PRECISION_MOD          ! For GEOS-Chem Precision (fp)
  USE TIME_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: CLEANUP_DRYDEP
  PUBLIC :: DO_DRYDEP
  PUBLIC :: INIT_DRYDEP
  PUBLIC :: INIT_WEIGHTSS
!
! !PUBLIC DATA MEMBERS:
!
  PUBLIC :: DEPNAME
  PUBLIC :: NUMDEP
  PUBLIC :: NTRAIND
  PUBLIC :: IDEP,   IRGSS,  IRAC, IRCLS
  PUBLIC :: IRGSO,  IRLU,   IRI,  IRCLO, DRYCOEFF
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Baldocchi, D.D., B.B. Hicks, and P. Camara, "A canopy stomatal
!        resistance model for gaseous deposition to vegetated surfaces",
!        Atmos. Environ. 21, 91-101, 1987.
!  (2 ) Brutsaert, W., "Evaporation into the Atmosphere", Reidel, 1982.
!  (3 ) Businger, J.A., et al., "Flux-profile relationships in the atmospheric
!        surface layer", J. Atmos. Sci., 28, 181-189, 1971.
!  (4 ) Dwight, H.B., "Tables of integrals and other mathematical data",
!        MacMillan, 1957.
!  (5 ) Guenther, A., and 15 others, A global model of natural volatile
!         organic compound emissions, J. Geophys. Res., 100, 8873-8892, 1995.
!  (6 ) Hicks, B.B., and P.S. Liss, "Transfer of SO2 and other reactive
!        gases across the air-sea interface", Tellus, 28, 348-354, 1976.
!  (7 ) Jacob, D.J., and S.C. Wofsy, "Budgets of reactive nitrogen,
!        hydrocarbons, and ozone over the Amazon forest during the wet season",
!        J.  Geophys. Res., 95, 16737-16754, 1990.
!  (8 ) Jacob, D.J., et al, "Deposition of ozone to tundra", J. Geophys. Res.,
!        97, 16473-16479, 1992.
!  (9 ) Levine, I.N., "Physical Chemistry, 3rd ed.", McGraw-Hill,
!        New York, 1988.
!  (10) Munger, J.W., et al, "Atmospheric deposition of reactive nitrogen
!        oxides and ozone in a temperate deciduous forest and a sub-arctic
!        woodland", J. Geophys. Res., in press, 1996.
!  (11) Walcek, C.J., R.A. Brost, J.S. Chang, and M.L. Wesely, "SO2, sulfate,
!        and HNO3 deposition velocities computed using regional landuse and
!        meteorological data", Atmos. Environ., 20, 949-964, 1986.
!  (12) Wang, Y.H., paper in preparation, 1996.
!  (13) Wesely, M.L, "Improved parameterizations for surface resistance to
!        gaseous dry deposition in regional-scale numerical models",
!        Environmental Protection Agency Report EPA/600/3-88/025,
!        Research Triangle Park (NC), 1988.
!  (14) Wesely, M. L., Parameterization of surface resistance to gaseous dry
!        deposition in regional-scale numerical models.  Atmos. Environ., 23
!        1293-1304, 1989.
!  (15) Price, H., L. Jaeglé, A. Rice, P. Quay, P.C. Novelli, R. Gammon,
!        Global Budget of Molecular Hydrogen and its Deuterium Content:
!        Constraints from Ground Station, Cruise, and Aircraft Observations,
!        submitted to J. Geophys. Res., 2007.
!  (16) Karl, T., Harley, P., Emmons, L., Thornton, B., Guenther, A., Basu, C.,
!        Turnipseed, A., and Jardine, K.: Efficient Atmospheric Cleansing of
!        Oxidized Organic Trace Gases by Vegetation, Science, 330, 816-819,
!        10.1126/science.1192534, 2010.
!  (17) JaeglÃ©, L., Shah, V.,et al (2018). Nitrogen oxides emissions, chemistry,
!        deposition,and export over the Northeast United States during the
!        WINTER aircraft campaign. J Geophys Res: Atmospheres, 123.
!        https://doi.org/10.1029/2018JD029133
!
! !REVISION HISTORY:
!  27 Jan 2003 - R. Yantosca - Moved standalone routines into this module
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  INTEGER,  PARAMETER :: NR_MAX    = 200       ! # of seasalt bins
  INTEGER,  PARAMETER :: NDRYDTYPE = 11        ! # of drydep land types
!
! PRIVATE TYPES:
!
  !========================================================================
  !  MODULE VARIABLES:
  !
  !  DRYDHNO3  : Internal flag for location of HNO3 in DEPVEL
  !  DRYDNO2   : Internal flag for location of NO2  in DEPVEL
  !  DRYDPAN   : Internal flag for location of PAN  in DEPVEL
  !  NUMDEP    : Actual number of drydep species
  !  NWATER    : Number of Olson's surface types that are water
  !  AIROSOL   : Array flags to denote aerosol drydep species
  !  IDEP      : ID #'s for dry deposition surface types
  !  IRAC      : ???       resistance for drydep land type
  !  IRCLO     : ???       resistance for drydep land type
  !  IRCLS     : ???       resistance for drydep land type
  !  IRGSO     : ???       resistance for drydep land type
  !  IRGSS     : ???       resistance for drydep land type
  !  IRI       : Internal  resistance for drydep land types
  !  IRLU      : Cuticular resistance for drydep land types
  !  IVSMAX    : ???       resistance for drydep land type
  !  IWATER    : ID #'s for Olson surface types that are water
  !  IZO       : Roughness heights for each Olson surface type
  !  NDVZIND   : Index array for ordering drydep species in DEPVEL
  !  NTRAIND   : Stores species numbers of drydep species
  !  PBLFRAC   : Array for multiplicative factor for drydep freq
  !  DRYCOEFF  : Polynomial coefficients for dry deposition
  !  HSTAR     : Henry's law constant
  !  F0        : Reactivity factor for biological oxidation
  !  XMW       : Molecular weight of drydep species [kg]
  !  A_RADI    : Radius of aerosol for size-resolved drydep [um]
  !  A_DEN     : Density of aerosol for size-res'd drydep [kg/m3]
  !  DEPNAME   : Names of dry deposition species
  !
  !  NOTE: these variables are defined in CMN_SIZE_mod.F
  !    NTYPE     : Max # of landtypes / grid box
  !    NPOLY     : Number of drydep polynomial coefficients
  !    NSURFTYPE : Number of Olson land types
  !
  !  NOTE: these grid-dependent variables are defined in State_Chm_Mod.F90
  !    DEPSAV    : Array containing dry deposition frequencies [s-1]
  !========================================================================

  ! Scalars
  INTEGER                        :: NUMDEP,   NWATER
  INTEGER                        :: DRYHg0,   DRYHg2,   DryHgP
  INTEGER                        :: id_ACET,  id_ALD2,  id_O3
  INTEGER                        :: id_MENO3, id_ETNO3
  INTEGER                        :: id_NK1
  INTEGER                        :: id_HNO3,  id_PAN,   id_IHN1

  ! Arrays for Baldocchi drydep polynomial coefficients
  REAL(fp), TARGET               :: DRYCOEFF(NPOLY    )

  ! Arrays that hold information for each of the 74 Olson land types
  INTEGER                        :: INDOLSON(NSURFTYPE )
  INTEGER                        :: IDEP    (NSURFTYPE )
  INTEGER                        :: IZO     (NSURFTYPE )
  INTEGER                        :: IWATER  (NSURFTYPE )

  ! Arrays that hold information for each of the 11 drydep land types
  INTEGER                        :: IDRYDEP (NDRYDTYPE)
  INTEGER                        :: IRAC    (NDRYDTYPE)
  INTEGER                        :: IRCLO   (NDRYDTYPE)
  INTEGER                        :: IRCLS   (NDRYDTYPE)
  INTEGER                        :: IRGSS   (NDRYDTYPE)
  INTEGER                        :: IRGSO   (NDRYDTYPE)
  INTEGER                        :: IRI     (NDRYDTYPE)
  INTEGER                        :: IRLU    (NDRYDTYPE)
  INTEGER                        :: IVSMAX  (NDRYDTYPE)

  ! Arrays that hold information about the dry-depositing species
  LOGICAL,           ALLOCATABLE :: AIROSOL (:    ) ! Is Aerosol? (T/F)
  INTEGER,           ALLOCATABLE :: NDVZIND (:    ) ! Drydep index
  INTEGER,           ALLOCATABLE :: FLAG    (:    ) ! Drydep scaling flag
  INTEGER,           ALLOCATABLE :: NTRAIND (:    ) ! Species index
  REAL(f8),          ALLOCATABLE :: HSTAR   (:    ) ! Henry's K0 [M/atm]
  REAL(f8),          ALLOCATABLE :: KOA     (:    ) ! POP's KOA
  REAL(f8),          ALLOCATABLE :: F0      (:    ) ! Reactivity factor [1]
  REAL(f8),          ALLOCATABLE :: XMW     (:    ) ! Mol wt. [kg/mol]
  REAL(f8),          ALLOCATABLE :: A_RADI  (:    ) ! Aer radius [m]
  REAL(f8),          ALLOCATABLE :: A_DEN   (:    ) ! Aer density [kg/m3]
  CHARACTER(LEN=14), ALLOCATABLE :: DEPNAME (:    ) ! Species name

  REAL(f4), POINTER :: HCO_Iodide(:,:)   => NULL()
  REAL(f4), POINTER :: HCO_Salinity(:,:) => NULL()

  ! Allocatable arrays
  REAL(f8),          ALLOCATABLE :: DMID    (:    )
  REAL(f8),          ALLOCATABLE :: SALT_V  (:    )

  !=================================================================
  ! MODULE ROUTINES -- follow below the "CONTAINS" statement
  !=================================================================
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_drydep
!
! !DESCRIPTION: Subroutine DO\_DRYDEP is the driver for the GEOS-CHEM dry
!  deposition scheme. DO\_DRYDEP calls DEPVEL to compute deposition velocities
!  [m/s], which are then converted to [cm/s].  Drydep frequencies are also
!  computed. (lwh, gmg, djj, 1989, 1994; bmy, 2/11/03, 5/25/05)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_DRYDEP( Input_Opt,  State_Chm, State_Diag, &
                        State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_ERROR_MOD
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
    USE HCO_DIAGN_MOD,      ONLY : Diagn_Update
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE Time_Mod,           ONLY : Get_Ts_Chem
    USE UnitConv_Mod
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
!  NOTE: Modeled aerosol dry deposition velocities over snow and ice surfaces
!  in the Arctic are much higher than estimated from measured values (e.g.,
!  Ibrahim et al. [1983]; Duan et al. [1988]; Nilsson and Rannik [2001]).
!  We will impose a dry deposition velocity of 0.03 cm/s for all aerosols
!  over snow and ice surfaces. (Jenny Fisher, 01 Aug 2011)
!
!  References (see full citations above):
!  ============================================================================
!  (1 ) Wesely, M. L., 1989
!  (2 ) Jacob, D.J., and S.C. Wofsy, 1990
!
! !REVISION HISTORY:
!  19 Nov 2002 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I,   J,   L,   D,   N,  NDVZ,  A
    REAL(f8)           :: DVZ, THIK
    CHARACTER(LEN=255) :: ErrMsg,  ThisLoc


    ! Arrays
    LOGICAL       :: LSNOW (State_Grid%NX,State_Grid%NY) ! Flag for snow/ice on the sfc
    REAL(f8)      :: CZ1   (State_Grid%NX,State_Grid%NY) ! Midpt ht of 1st model lev[m]
    REAL(f8)      :: TC0   (State_Grid%NX,State_Grid%NY) ! Grid box sfc temperature [K]
    REAL(f8)      :: ZH    (State_Grid%NX,State_Grid%NY) ! PBL height [m]
    REAL(f8)      :: OBK   (State_Grid%NX,State_Grid%NY) ! Monin-Obhukov Length [m]
    REAL(f8)      :: CFRAC (State_Grid%NX,State_Grid%NY) ! Column cloud frac [unitless]
    REAL(f8)      :: RADIAT(State_Grid%NX,State_Grid%NY) ! Solar radiation [W/m2]
    REAL(f8)      :: USTAR (State_Grid%NX,State_Grid%NY) ! Grid box friction vel [m/s]
    REAL(f8)      :: RHB   (State_Grid%NX,State_Grid%NY) ! Relative humidity [unitless]
    REAL(f8)      :: DVEL  (State_Grid%NX,State_Grid%NY,NUMDEP) ! Drydep velocities [m/s]
    REAL(f8)      :: PRESSU(State_Grid%NX,State_Grid%NY) ! Local surface pressure [Pa]
    REAL(f8)      :: W10   (State_Grid%NX,State_Grid%NY) ! 10m windspeed [m/s]
    REAL(f8)      :: AZO   (State_Grid%NX,State_Grid%NY)        ! Z0, per (I,J) square
    REAL(f8)      :: SUNCOS_MID(State_Grid%NX,State_Grid%NY)    ! COS(SZA) @ midpoint of the
                                                  !  current chemistry timestep

    ! Pointers
    REAL(fp), POINTER :: DEPSAV (:,:,:   )      ! Dry deposition frequencies [s-1]

    ! For ESMF, need to assign these from Input_Opt
    LOGICAL       :: PBL_DRYDEP
    LOGICAL       :: prtDebug

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    !=================================================================
    ! DO_DRYDEP begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Initialize
    SpcInfo => NULL()
    ErrMsg  = ''
    ThisLoc = ' -> at Do_DryDep  (in module GeosCore/drydep_mod.F90)'

    ! Point to columns of derived-type object fields
    DEPSAV     => State_Chm%DryDepSav

    ! Copy values from the Input Options object to local variables
    PBL_DRYDEP = Input_Opt%PBL_DRYDEP
    prtDebug   = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Get fields for oceanic O3 drydeposition
    CALL HCO_GetPtr( HcoState, 'surf_iodide',   HCO_Iodide,   RC )
    CALL HCO_GetPtr( HcoState, 'surf_salinity', HCO_Salinity, RC )

    ! Call METERO to obtain meterological fields (all 1-D arrays)
    ! Added sfc pressure as PRESSU and 10m windspeed as W10
    !  (jaegle 5/11/11, mpayer 1/10/12)
    CALL METERO( State_Grid, State_Met, CZ1,     TC0, OBK,  CFRAC, &
                 RADIAT,     AZO,       USTAR,   ZH,        LSNOW, &
                 RHB,        PRESSU,    W10,     SUNCOS_MID        )

    ! Call DEPVEL to compute dry deposition velocities [m/s]
    CALL DEPVEL( Input_Opt, State_Chm,  State_Diag, State_Grid, &
                 State_Met, RADIAT,     TC0,        SUNCOS_MID, &
                 F0,        HSTAR,      XMW,        AIROSOL,    &
                 USTAR,     CZ1,        OBK,        CFRAC,      &
                 ZH,        LSNOW,      DVEL,       AZO,        &
                 RHB,       PRESSU,     W10,        RC          )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in call to "DEPVEL!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Compute dry deposition frequencies; archive diagnostics
    !=================================================================
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, THIK, D, N, NDVZ, DVZ, SpcInfo, A )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! THIK = thickness of surface layer [m]
       THIK   = State_Met%BXHEIGHT(I,J,1)

       ! Now we calculate drydep throughout the entire PBL.
       ! Make sure that the PBL depth is greater than or equal
       ! to the thickness of the 1st layer (rjp, bmy, 7/21/03)
       ! Add option for non-local PBL mixing scheme: THIK must
       ! be the first box height. (Lin, 03/31/09)
       ! Now use PBL_DRYDEP instead of LNLPBL (ckeller, 3/5/15).
       IF (PBL_DRYDEP) THIK = MAX( ZH(I,J), THIK )

       ! Loop over drydep species
       DO D = 1, State_Chm%nDryDep

          ! GEOS-CHEM species number
          N       =  State_Chm%Map_DryDep(D)

          ! Get info about this species from the database
          SpcInfo => State_Chm%SpcData(N)%Info

          ! Get the "DryAltID" index, that is used to archive species
          ! concentrations at a user-defined altitude above the surface
          A       =  SpcInfo%DryAltID

          ! Index of drydep species in the DVEL array
          ! as passed back from subroutine DEPVEL
          NDVZ    =  NDVZIND(D)

          ! Dry deposition velocity [cm/s]
          DVZ     =  DVEL(I,J,NDVZ) * 100.e+0_f8

          ! Scale relative to specified species (krt, 3/1/15)
          IF ( FLAG(D) .eq. 1 )  THEN

             ! Scale species to HNO3
             DVZ = DVZ * sqrt(State_Chm%SpcData(id_HNO3)%Info%MW_g) &
                       / sqrt(SpcInfo%MW_g)

          ELSE IF ( FLAG(D) .eq. 2 ) THEN
             
             ! Scale species to PAN
             DVZ = DVZ * sqrt(State_Chm%SpcData(id_PAN)%Info%MW_g) &
                       / sqrt(SpcInfo%MW_g)

          ELSE IF ( FLAG(D) .eq. 3 ) THEN

             ! Scale species to ISOPN
             DVZ = DVZ * sqrt(State_Chm%SpcData(id_IHN1)%Info%MW_g) &
                       / sqrt(SpcInfo%MW_g)

          ENDIF

          !-----------------------------------------------------------
          ! Special treatment for snow vs. ice
          !-----------------------------------------------------------
          IF ( LSNOW(I,J) ) THEN

             !-------------------------------------
             ! %%% SURFACE IS SNOW OR ICE %%%
             !-------------------------------------
             IF ( SpcInfo%DD_DvzAerSnow > 0.0_fp ) THEN

                ! For most aerosol species (basically everything
                ! except sea salt and dust species), we just set
                ! the deposition velocity over snow to a fixed value.
                ! (Modification by Jenny Fisher, dated 8/1/11)
                DVZ = DBLE( SpcInfo%DD_DvzAerSnow )

             ELSE

                ! Otherwise, enforce a minimum drydep velocity over snow
                ! (cf. the GOCART model).  NOTE: In practice this will
                ! only apply to the species SO2, SO4, MSA, NH3, NH4, NIT.
                DVZ = MAX( DVZ, DBLE( SpcInfo%DD_DvzMinVal(1) ) )

             ENDIF

          ELSE

             !-------------------------------------
             ! %%% SURFACE IS NOT SNOW OR ICE %%%
             !-------------------------------------

             ! Enforce a minimum drydep velocity over land (cf. the
             ! GOCART model).  NOTE: In practice this will only apply
             ! to the species SO2, SO4, MSA, NH3, NH4, NIT.
             DVZ = MAX( DVZ, DBLE( SpcInfo%DD_DvzMinVal(2) ) )

          ENDIF

          !-----------------------------------------------------------
          ! Special treatment for ACETONE
          !-----------------------------------------------------------

          ! For ACET, we need to only do drydep over the land
          ! and not over the oceans.
          IF ( N == id_ACET ) THEN
             IF ( State_Met%IsLand(I,J) ) THEN
                DVZ = 0.1e+0_f8
             ELSE
                DVZ = 0e+0_f8
             ENDIF
          ENDIF

          !-----------------------------------------------------------
          ! Special treatment for ALD2
          !-----------------------------------------------------------

          ! For ALD2, we need to only do drydep over the land
          ! and not over the oceans.
          IF ( N == id_ALD2 ) THEN
             IF ( .not. State_Met%IsLand(I,J) ) THEN
                DVZ = 0e+0_f8
             ENDIF
          ENDIF

          !-----------------------------------------------------------
          ! Special treatment for MENO3
          !-----------------------------------------------------------

          ! For MENO3, we need to only do drydep over the land
          ! and not over the oceans.
          IF ( N == id_MENO3 ) THEN
             IF ( .not. State_Met%IsLand(I,J) ) THEN
                DVZ = 0e+0_f8
             ENDIF
          ENDIF

          !-----------------------------------------------------------
          ! Special treatment for ETNO3
          !-----------------------------------------------------------

          ! For ETNO3, we need to only do drydep over the land
          ! and not over the oceans.
          IF ( N == id_ETNO3 ) THEN
             IF ( .not. State_Met%IsLand(I,J) ) THEN
                DVZ = 0e+0_f8
             ENDIF
          ENDIF

          !-----------------------------------------------------------
          ! Compute drydep frequency and update diagnostics
          !-----------------------------------------------------------

          ! Dry deposition frequency [1/s]
          DEPSAV(I,J,D) = ( DVZ / 100.e+0_f8 ) / THIK

          ! Archive dry dep velocity [cm/s]
          IF ( State_Diag%Archive_DryDepVel ) THEN
             State_Diag%DryDepVel(I,J,D) = DVZ
          ENDIF

          ! Archive dry dep velocity [cm/s] only for those species
          ! that are requested at a given altitude (e.g. 10m)
          IF ( State_Diag%Archive_DryDepVelForALT1 ) THEN
             IF ( A > 0 ) THEN
                State_Diag%DryDepVelForALT1(I,J,A) = DVZ
             ENDIF
          ENDIF

          ! Free pointer
          SpcInfo => NULL()
       ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !### Debug
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### DO_DRYDEP: after dry dep' )
    ENDIF

    ! Nullify pointers
    NULLIFY( DEPSAV )

  END SUBROUTINE DO_DRYDEP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: OceanO3
!
! !DESCRIPTION: Function OCEANO3 calculates the dry deposition velcoity of O3
!     to the ocean using method described in Pound et.al (2019)
!     currently under discussion in ACPD. 
!     Accounts for the turbulence of the ocean surface,Iodide concentration
!     and surface temperature effects on the dry deposition velocity of
!     ozone to the ocean.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE OCEANO3( TEMPK, USTAR, DEPV, I, J )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
!     
! !INPUT PARAMETERS: 
!
    REAL(f8), INTENT(IN)         :: TEMPK ! Temperature [K]
    REAL(f8), INTENT(IN)         :: USTAR ! Fictional Velocity [m/s]
    INTEGER,  INTENT(IN)         :: I,J
    REAL(f8), INTENT(OUT)        :: DEPV  ! the new deposition vel [cm/s]
! 
! !REVISION HISTORY: 
!  21 Aug 2018 - R. Pound - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(f8) :: a,D,DelM,b,PSI,LAM,EP,USTARWater,K0,K1,Iodide

    !=================================================================
    ! OCEANO3 begins here!
    !=================================================================
      
    USTARWater = 0.0345_f8*USTAR !waterside friction velocity
    
    Iodide = HCO_Iodide(I,J)*1.0E-9_f8 !retrieve iodide from HEMCO
     
    a = Iodide*EXP((-8772.2/TEMPK)+51.5) !chemical reactivity

    D = 1.1E-6*EXP(-1896.0/TEMPK) ! diffusivity

    DelM = SQRT(D/a) ! reaction-diffusion length
      
    b = 2.0_f8/(0.4_f8*USTARWater)
   
    LAM = DelM*SQRT(a/D) ! this cancels to 1 but here for completeness
                         ! of equations

    EP = SQRT(2.0_f8*a*b*(DelM+(b*D/2.0_f8)))

    PSI = EP/SQRT(a*b**2*D)

    CALL K0K1_APROX(EP,K0,K1)

    DEPV = SQRT(a*D)*((PSI*K1*COSH(LAM)+K0*SINH(LAM))/(PSI*K1* &
           SINH(LAM)+K0*COSH(LAM)))
   
  END SUBROUTINE OCEANO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: K0K1_APROX
!
! !DESCRIPTION:  Function to estimate the modified Bessel functions of
!     the second kind order zero and one (K0,K1). Approach initially based on 
!     that described in Numerical Recipes in Fortran 90 second edition
!     (1996). This implementation is designed to be specific to the use
!     case required for calculating oceanic deposition velocity. Uses a
!     polynomial fit of each type of modified bessel function to
!     estimate the value of the function for each input. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE K0K1_APROX( input_arg, K0, K1 )
!     
! !INPUT PARAMETERS: 
!     
    REAL(f8), INTENT(IN)  :: input_arg !the value we want the soln for
    REAL(f8), INTENT(OUT) :: K0,K1 !the values of the modified bessel fncs
!     
! !REVISION HISTORY: 
!  21 Aug 2018 - R. Pound - Initial version    
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!     
! !LOCAL VARIABLES:
!    
    REAL(f8), DIMENSION(7) :: coeff !coefficients for polynomial fit
                                    ! of each bessel function
    REAL(f8)               :: I0,I1 !modified bessel functions of
                                    ! first kind order 0 and 1

    ! determine which fit method is best for the bessel functions
    IF (input_arg <= 2.0_f8) THEN
       ! begin the calculation of k0 by estimating i0
       coeff = (/1.0,3.5156229,3.0899424,1.2067492,0.2659732, &
                 0.360768e-1,0.45813e-2/)
       I0 = poly_fit((input_arg/3.75_f8)**2,coeff)
       !now we can use this estimate of i0 to calculate k0
       coeff = (/-0.57721566,0.42278420,0.23069756,0.3488590e-1, &
                  0.262698e-2,0.10750e-3,0.74e-5/)
       K0 = (-LOG(0.5_f8*input_arg)*I0)+ &
             poly_fit(0.25_f8*input_arg**2,coeff)

       !begin the calculation of k0 by estimating i1
       coeff = (/0.5,0.87890594,0.51498869,0.15084934,0.2658733e-1, &
                 0.301532e-2,0.32411e-3/)
       I1 = input_arg*poly_fit((input_arg/3.75_f8)**2,coeff)
       ! now we can use this to estimate to get a value for k1
       coeff = (/1.0,0.15443144,-0.67278579,-0.18156897, &
                -0.1919402e-1,-0.110404e-2,-0.4686e-4/)
       K1 = (LOG(0.5_f8*input_arg)*I1)+(1.0_f8/input_arg)* &
             poly_fit(0.25_f8*input_arg**2,coeff)
    ELSE !use a different approximation that doesn't need I0/I1
       coeff = (/1.25331414,-0.7832358e-1,0.2189568e-1,-0.1062446e-1, &
                 0.587872e-2,-0.251540e-2,0.53208e-3/)
       K0 = (EXP(-input_arg)/SQRT(input_arg))* &
             poly_fit((2.0_f8/input_arg),coeff)
       coeff = (/1.25331414,0.23498619,-0.3655620e-1,0.1504268e-1, &
                -0.780353e-2,0.325614e-2,-0.68245e-3/)
       K1 = (EXP(-input_arg)/SQRT(input_arg))* &
             poly_fit((2.0_f8/input_arg),coeff)
    ENDIF

  END SUBROUTINE K0K1_APROX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: poly_fit
!
! !DESCRIPTION: Function to calculate the value of a polynomial fit,
! used in the K0K1_APPROX function in estimating the values of a
! modified bessel function.
!\\
!\\
! !INTERFACE:
!
  FUNCTION poly_fit ( input, coeffs )
!     
! !INPUT PARAMETERS: 
!     
    REAL(f8), INTENT(IN)               :: input
    REAL(f8), DIMENSION(:), INTENT(IN) :: coeffs
!
! !REVISION HISTORY: 
!  21 Aug 2018 - R. Pound - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!     
! !LOCAL VARIABLES:
!     
    REAL(f8)                           :: poly_fit
    INTEGER                            :: i

    poly_fit = 0

    DO i = 1,7,1
       poly_fit = poly_fit+coeffs(i)*input**i
    ENDDO

  END FUNCTION poly_fit
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: metero
!
! !DESCRIPTION: Subroutine METERO calculates meteorological constants needed
!  for the dry deposition velocity module. (lwh, gmg, djj, 1989, 1994; bmy,
!  10/3/05)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE METERO( State_Grid, State_Met, CZ1,  TC0, OBK, CFRAC, &
                     RADIAT,     AZO,       USTR, ZH,  LSNOW,      &
                     RHB,        PRESSU,    W10,  SUNCOS_MID        )
!
! !USES:
!
    USE Calc_Met_Mod,       ONLY : GET_OBK
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,  INTENT(OUT) :: LSNOW (State_Grid%NX,State_Grid%NY)  ! Flag for denoting snow/ice
    REAL(f8), INTENT(OUT) :: CZ1   (State_Grid%NX,State_Grid%NY)  ! Midpt ht of 1st model lev [m]
    REAL(f8), INTENT(OUT) :: TC0   (State_Grid%NX,State_Grid%NY)  ! Grid box sfc temp [K]
    REAL(f8), INTENT(OUT) :: OBK   (State_Grid%NX,State_Grid%NY)  ! Monin-Obhukov length [m]
    REAL(f8), INTENT(OUT) :: CFRAC (State_Grid%NX,State_Grid%NY)  ! Column cloud fraction [unitless]
    REAL(f8), INTENT(OUT) :: RADIAT(State_Grid%NX,State_Grid%NY)  ! Solar radiation @ ground [W/m2]
    REAL(f8), INTENT(OUT) :: RHB   (State_Grid%NX,State_Grid%NY)  ! Rel humidity at sfc [unitless]
    REAL(f8), INTENT(OUT) :: USTR  (State_Grid%NX,State_Grid%NY)  ! Friction velocity [m/s]
    REAL(f8), INTENT(OUT) :: ZH    (State_Grid%NX,State_Grid%NY)  ! PBL height [m]
    REAL(f8), INTENT(OUT) :: PRESSU(State_Grid%NX,State_Grid%NY)  ! Local surface press [Pa]
    REAL(f8), INTENT(OUT) :: W10   (State_Grid%NX,State_Grid%NY)  ! 10 meter windspeed [m/s]
    REAL(f8), INTENT(OUT) :: SUNCOS_MID(State_Grid%NX,State_Grid%NY) ! COS(SZA) @ midpt of current chem timestep
    REAL(f8), INTENT(OUT) :: AZO(State_Grid%NX,State_Grid%NY) ! Roughness heights, by grid box
!
! !REMARKS:
!                                                                             .
!  References (see full citations above):
!  ============================================================================
!  (1 ) Wesely, M. L., 1989.
!  (2 ) Jacob, D.J., and S.C. Wofsy, 1990
!
! !REVISION HISTORY:
!  20 Nov 2002 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I,  J
    REAL(f8) :: THIK
    REAL(f8) :: SP
    REAL(f8) :: SFCWINDSQR

    !=================================================================
    ! METERO begins here!
    !=================================================================

    ! Loop over surface grid boxes
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, THIK, SP, SFCWINDSQR )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! THIK = thickness of layer 1 [m]
       THIK        = State_Met%BXHEIGHT(I,J,1)

       ! Midpoint height of first model level [m]
       CZ1(I,J)    = THIK / 2.0e+0_f8

       ! Local surface pressure [hPa] (mpayer, 1/10/12)
       ! Use moist air pressure for mean free path  (ewl, 3/2/15)
       SP          = State_Met%PEDGE(I,J,1)

       ! Convert from hPa to Pa for SFCPRESS
       PRESSU(I,J) = SP * 1.e+2_f8

       !==============================================================
       ! Return meterological quantities as 1-D arrays for DEPVEL
       !==============================================================

       ! Roughness height [m]
       AZO(I,J)    = State_Met%Z0(I,J)

       ! Column cloud fraction [unitless]
       CFRAC(I,J)  = State_Met%CLDFRC(I,J)

       ! Set logical LSNOW if snow and sea ice (ALBEDO > 0.4)
       LSNOW(I,J)  = ( State_Met%ALBD(I,J) > 0.4 )

       ! Monin-Obhukov length [m]
       OBK(I,J)    = GET_OBK( I, J, State_Met )

       ! Solar insolation @ ground [W/m2]
       RADIAT(I,J) = State_Met%SWGDN(I,J)

       ! Surface temperature [K]
       TC0(I,J)    = State_Met%TS(I,J)

       ! Friction velocity [m/s]
       USTR(I,J)   = State_Met%USTAR(I,J)

       ! Mixed layer depth [m]
       ZH(I,J)     = State_Met%PBL_TOP_m(I,J)

       ! Relative humidity @ surface [unitless] (bec, bmy, 4/13/05)
       !RHB(I,J)    = MIN( 0.99e+0_f8, RH(I,J,1) * 1.d-2 )
       !  changed to 98% due to vapor pressure lowering above sea water
       !  (Lewis & Schwartz, 2004)
       !  jaegle (5/11/11)
       RHB(I,J)    = MIN( 0.98e+0_f8, State_Met%RH(I,J,1) * 1.e-2_f8 )

       ! 10m windspeed [m/s] (jaegle 5/11/11)
       SFCWINDSQR  = State_Met%U10M(I,J)**2 + State_Met%V10M(I,J)**2
       W10(I,J)    = SQRT( SFCWINDSQR )

       ! Cosine of solar zenith angle at midpoint
       ! of the current chemistry timestep.
       SUNCOS_MID(I,J) = State_Met%SUNCOSmid(I,J)

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE METERO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: depvel
!
! !DESCRIPTION: Subroutine DEPVEL computes the dry deposition velocities using
!  a resistance-in-series model.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DEPVEL( Input_Opt, State_Chm, State_Diag, State_Grid, &
                     State_Met, RADIAT,    TEMP,       SUNCOS,     &
                     F0,        HSTAR,     XMW,        AIROSOL,    &
                     USTAR,     CZ1,       OBK,        CFRAC,      &
                     ZH,        LSNOW,     DVEL,       ZO,         &
                     RHB,       PRESSU,    W10,        RC          )
!
! !USES:
!
    USE CMN_SIZE_MOD,       ONLY : NTYPE
    USE Drydep_Toolbox_Mod, ONLY : BioFit
    USE ErrCode_Mod
    USE ERROR_MOD
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
#ifdef APM
    USE APM_INIT_MOD,       ONLY : APMIDS
    USE APM_INIT_MOD,       ONLY : RDRY, RSALT, RDST, DENDST
    USE APM_DRIV_MOD,       ONLY : GFTOT3D, DENWET3D, MWSIZE3D
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt      ! Input Options object
    TYPE(DgnState), INTENT(IN) :: State_Diag     ! Diagnostics state object
    TYPE(GrdState), INTENT(IN) :: State_Grid     ! Grid state object
    TYPE(MetState), INTENT(IN) :: State_Met      ! Meteorology state object

    REAL(f8), INTENT(IN) :: RADIAT (State_Grid%NX,State_Grid%NY) ! Solar radiation [W/m2]
    REAL(f8), INTENT(IN) :: TEMP   (State_Grid%NX,State_Grid%NY) ! Temperature [K]
    REAL(f8), INTENT(IN) :: SUNCOS (State_Grid%NX,State_Grid%NY) ! Cos of solar zenith angle
    LOGICAL,  INTENT(IN) :: AIROSOL(NUMDEP)      ! =T denotes aerosol species
    REAL(f8), INTENT(IN) :: F0     (NUMDEP)      ! React. factor for oxidation
                                                 !  of biological substances
    REAL(f8), INTENT(IN) :: HSTAR  (NUMDEP)      ! Henry's law constant
    REAL(f8), INTENT(IN) :: XMW    (NUMDEP)      ! Molecular weight [kg/mol]
    REAL(f8), INTENT(IN) :: USTAR  (State_Grid%NX,State_Grid%NY) ! Friction velocity [m/s]
    REAL(f8), INTENT(IN) :: CZ1    (State_Grid%NX,State_Grid%NY) ! Alt @ which Vd is compute [m]
    REAL(f8), INTENT(IN) :: OBK    (State_Grid%NX,State_Grid%NY) ! Monin-Obhukov length [m]
    REAL(f8), INTENT(IN) :: CFRAC  (State_Grid%NX,State_Grid%NY) ! Surface cloud fraction
    REAL(f8), INTENT(IN) :: ZH     (State_Grid%NX,State_Grid%NY) ! Roughness height [m]
    REAL(f8), INTENT(IN) :: RHB    (State_Grid%NX,State_Grid%NY) ! Relative humidity [%]
    REAL(f8), INTENT(IN) :: PRESSU (State_Grid%NX,State_Grid%NY) ! Surface pressure [hPa]
    REAL(f8), INTENT(IN) :: W10    (State_Grid%NX,State_Grid%NY) ! Wind speed @ 10m altitude [m/s]
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(OUT) :: RC                  ! Success or failure?
    REAL(f8), INTENT(OUT) :: DVEL(State_Grid%NX,State_Grid%NY,NUMDEP) ! Drydep velocity [m/s]
!
! !REMARKS:
!  Need as landtype input for each grid square (I,J); see CMN_DEP_mod.F
!     IREG(I,J)       - # of landtypes in grid square
!     ILAND(I,J,LDT)  - Land type ID for element LDT =1, IREG(I,J)
!                        (could be from any source - mapped to deposition
!                        surface ID in input unit 65)
!     IJUSE(I,J,LDT) - Fraction ((per mil) of gridbox area occupied by
!                       land type element LDT
!                                                                             .
!  Need as leaf area index; see CMN_DEP_mod.F
!     XLAI(I,J,LDT)  - Leaf Area Index of land type element LDT
!                                                                             .
!  Need as meteorological input for each grid square(I,J) (passed):
!     RADIAT(I,J)    - Solar radiation in W m-2
!     TEMP(I,J)      - Surface air temperature in K
!     SUNCOS(I,J)    - Cosine of solar zenith angle
!     LSNOW(I,J)     - Logical for snow and sea ice
!     RHB(I,J)       - Relative humidity at the surface
!     PRESSU(I,J)    - Local surface pressure
!     W10(I,J)       - 10m wind speed
!                                                                             .
!  Need as input for each species K (passed):
!     F0(K)          - reactivity factor for oxidation of biological substances
!     HSTAR(K)       - Henry's Law constant
!     XMW(K)         - Molecular weight (kg/mole) of species K
!                      (used to calculate molecular diffusivities)
!     AIROSOL(K)     - LOGICAL flag (T = aerosol species;
!                                    F = gas-phase species)
!                                                                             .
!  Also need to call the following subroutines to read drydep input data:
!     READ_DRYDEP_INPUTS    - (in this module) Reads in Olson land type
!                             indices, dry deposition land type indices,
!                             default roughness heights, and polynomial
!                             coefficients.  (This supersedes MODIN, RDDRYCF)
!     COMPUTE_OLSON_LANDMAP - (in olson_landmap_mod.F90).  Reads in the
!                             Olson land types at native resolution and re-bins
!                             them on-the-fly to the GEOS-Chem grid resolution.
!                             (This supersedes RDLAND)
!     "rdlai.f"             - reads Leaf Area Indices from files "lai**.global"
!                                                                             .
!  Some variables used in the subroutine (passed):
!     LRGERA(I,J)    T -> stable atmosphere; a high aerodynamic resistance
!                        (RA=1.E4 m s-1) is imposed; else RA is calculated
!     USTAR(I,J)     - Friction velocity (m s-1)
!     CZ1(I,J)       - Altitude (m) at which deposition velocity is computed
!     OBK(I,J)       - Monin-Obukhov length (m): set to 1.E5 m under neutral
!                      conditions
!     CFRAC(I,J)     - Fractional cloud cover
!     ZH(I,J)        - Mixing depth (m)
!                                                                             .
!  Some variables used in the subroutine:
!     MAXDEP         - the maximum number of species for which the dry
!                      deposition calculation is done
!     ZO(LDT)        - Roughness height (m) for specific surface type indexed
!                      by LDT
!     RSURFC(K,LDT)  - Bulk surface resistance (s m-1) for species K to
!                      surface LDT
!     C1X(K)         - Total resistance to deposition (s m-1) for species K
!                                                                             .
!  Returned:
!     DVEL(I,J,K) - Deposition velocity (m s-1) of species K
!                                                                             .
!  References:
!  ============================================================================
!     Baldocchi, D.D., B.B. Hicks, and P. Camara, A canopy stomatal
!       resistance model for gaseous deposition to vegetated surfaces,
!       Atmos. Environ. 21, 91-101, 1987.
!     Brutsaert, W., Evaporation into the Atmosphere, Reidel, 1982.
!     Businger, J.A., et al., Flux-profile relationships in the atmospheric
!       surface layer, J. Atmos. Sci., 28, 181-189, 1971.
!     Dwight, H.B., Tables of integrals and other mathematical data,
!       MacMillan, 1957.
!     Guenther, A., and 15 others, A global model of natural volatile
!       organic compound emissions, J. Geophys. Res., 100, 8873-8892, 1995.
!     Hicks, B.B., and P.S. Liss, Transfer of SO2 and other reactive
!       gases across the air-sea interface, Tellus, 28, 348-354, 1976.
!     Jacob, D.J., and S.C. Wofsy, Budgets of reactive nitrogen,
!       hydrocarbons, and ozone over the Amazon forest during the wet season,
!       J.  Geophys. Res., 95, 16737-16754, 1990.
!     Jacob, D.J., and 9 others, Deposition of ozone to tundra,
!       J. Geophys. Res., 97, 16473-16479, 1992.
!     Levine, I.N., Physical Chemistry, 3rd ed., McGraw-Hill, New York, 1988.
!     Munger, J.W., and 8 others, Atmospheric deposition of reactive
!       nitrogen oxides and ozone in a temperate deciduous forest and a
!       sub-arctic woodland, J. Geophys. Res., in press, 1996.
!     Walcek, C.J., R.A. Brost, J.S. Chang, and M.L. Wesely, SO2, sulfate, and
!       HNO3 deposition velocities computed using regional landuse and
!       meteorological data, Atmos. Environ., 20, 949-964, 1986.
!     Wang, Y.H., paper in preparation, 1996.
!     Wesely, M.L, Improved parameterizations for surface resistance to
!       gaseous dry deposition in regional-scale numerical models,
!       Environmental Protection Agency Report EPA/600/3-88/025,
!       Research Triangle Park (NC), 1988.
!     Wesely, M.L., same title, Atmos. Environ., 23, 1293-1304, 1989.
!
! !REVISION HISTORY:
!  Circa 1990 - D.J. Jacob - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(f8) :: RI(NTYPE),RLU(NTYPE),RAC(NTYPE),RGSS(NTYPE)
    REAL(f8) :: RGSO(NTYPE),RCLS(NTYPE),RCLO(NTYPE)
    REAL(f8) :: RSURFC(NUMDEP,NTYPE)
    REAL(f8) :: C1X(NUMDEP),VD(NUMDEP),VK(NUMDEP)
    REAL(f8) :: ZO(State_Grid%NX,State_Grid%NY)

    LOGICAL  :: LDEP(NUMDEP)
    LOGICAL  :: LRGERA(State_Grid%NX,State_Grid%NY)

    REAL(f8) :: VDS
    REAL(f8) :: CZ,C1,RT,XNU,RAD0,RIX,GFACT,GFACI
    REAL(f8) :: RDC,RLUXX,RGSX,RCL,DTMP1,DTMP2,DTMP3,DTMP4
    REAL(f8) :: CZH,CKUSTR,REYNO,CORR1,CORR2,Z0OBK
    REAL(f8) :: RA,RB,DUMMY1,DUMMY2,DUMMY3,DUMMY4
    REAL(f8) :: XMWH2O,DAIR,TEMPK,TEMPC
    INTEGER  :: IOLSON,II,IW
    INTEGER  :: K,LDT
    REAL(f8) :: RCLX,RIXX    !,BIOFIT

    ! for corr O3, krt,11/2017
    REAL(f8) :: RA_Alt, DUMMY2_Alt, DUMMY4_Alt, Z0OBK_Alt

#ifdef TOMAS
    ! For TOMAS aerosol (win, 7/15/09)
    INTEGER  :: BIN
    REAL(f8) :: SIZ_DIA(State_Grid%NX,State_Grid%NY,IBINS)
    REAL(f8) :: SIZ_DEN(State_Grid%NX,State_Grid%NY,IBINS)
#endif

    ! Logical for snow and sea ice
    LOGICAL  ::LSNOW(State_Grid%NX,State_Grid%NY)

    ! Loop indices (bmy, 3/29/12)
    INTEGER  :: I, J

    ! Integer for the GEOS-CHEM species number
    INTEGER  :: N_SPC
    REAL(f8) :: alpha
    REAL(f8) :: DEPVw

    ! Size of DRYCOEFF (ckeller, 05/19/14)
    INTEGER  :: NN

    ! Added these to pass particle density & number to DUST_SFCRSII
    ! so as to avoid parallelization errors with TOMAS (bmy, 1/31/14)
    REAL(f8) :: DIAM, DEN

    ! Logical switch for POPS in order to use octanol-air partitioning instead
    ! of Henry's Law for scaling of cuticular resitances (clf, 1/3/2011)
    LOGICAL  :: IS_POPS

    ! Pointers to fields in State_Met
    INTEGER,  POINTER :: IREG(:,:)
    INTEGER,  POINTER :: ILAND(:,:,:)
    INTEGER,  POINTER :: IUSE(:,:,:)
    REAL(fp), POINTER :: XLAI(:,:,:)

    ! For making sure that all inputs to BIOFIT are of the same type
    REAL(fp)          :: XLAI_FP
    REAL(fp)          :: SUNCOS_FP
    REAL(fp)          :: CFRAC_FP

#ifdef MODEL_GEOS
    ! Archive Ra?
    REAL(fp)          :: RA2M,  Z0OBK2M
    REAL(fp)          :: RA10M, Z0OBK10M
    !LOGICAL           :: WriteRa2m
    !LOGICAL           :: WriteRa10m
#endif

    ! For the species database
    INTEGER                :: SpcId
    TYPE(Species), POINTER :: SpcInfo

    ! Strings
    CHARACTER(LEN=255)  :: ThisLoc
    CHARACTER(LEN=512)  :: ErrMsg
!
! !DEFINED PARAMETERS:
!
    ! Small number
    REAL(fp), PARAMETER :: SMALL = 1.0e-10_f8

    !=================================================================
    ! DEPVEL begins here!
    !=================================================================

    ! Assume success
    RC      =  GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at DEPVEL (in module GeosCore/drydep_mod.F)'

    ! Error check: PBL_DRYDEP must not be used together with non-local
    ! PBL mixing.
    IF ( Input_Opt%LNLPBL .AND. Input_Opt%PBL_DRYDEP ) THEN
       ErrMsg = 'PBL_DRYDEP must be disabled if LNLPBL=T!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Initialize pointers
    IREG    => State_Met%IREG
    ILAND   => State_Met%ILAND
    IUSE    => State_Met%IUSE
    XLAI    => State_Met%XLAI
    SpcInfo => NULL()

    ! Is this a POPs simmulation?
    IS_POPS = Input_Opt%ITS_A_POPS_SIM  ! clf, 1/3/2011

    ! Size of drycoeff (ckeller, 05/19/14)
    NN = SIZE(DRYCOEFF)

#ifdef MODEL_GEOS
    ! Logical flag for Ra (ckeller, 12/29/17)
    State_Chm%DryDepRa2m  = 0.0_fp
    State_Chm%DryDepRa10m = 0.0_fp
    !WriteRa2m = ASSOCIATED ( State_Diag%DryDepRa2m )
    !IF ( WriteRa2m ) THEN
    !   State_Diag%DryDepRa2m = 0.0_fp
    !ENDIF
    !WriteRa10m = ASSOCIATED ( State_Diag%DryDepRa10m )
    !IF ( WriteRa10m ) THEN
    !   State_Diag%DryDepRa10m = 0.0_fp
    !ENDIF
#endif

    !***********************************************************************
    !
    !** If LDEP(K)=F, species does not deposit.
    !** Deposition is applied only to species with LDEP=T.
    DO K = 1,NUMDEP

       ! Better test for depositing species K: We need both HSTAR and XMW
       ! to be nonzero, OR the value of AIROSOL to be true.  This should
       ! avoid any futher floating point invalid issues caused by putting
       ! a zero value in a denominator. (bmy, 8/29/13)
       ! (bmy, 8/29/13)
       IF ( ( HSTAR(K) > 0e+0_f8   .and. &
              XMW  (K) > 0e+0_f8 ) .or.  AIROSOL(K) ) THEN
          LDEP(K) = .TRUE.
       ELSE
          LDEP(K) = .FALSE.
       ENDIF
    ENDDO

    ! Initialize DVEL
    DVEL = 0.0e+0_f8

    !***********************************************************************
    !*
    !*    Begin section for computing deposition velocities
    !*
#ifdef TOMAS
    !=================================================================
    ! FOR TOMAS MICROPHYSICS:
    !
    ! Calculate 30-bin aerosol size and diameter at each grid box
    ! This is done differently than 2-bin seasalt and 4-bin dust
    ! because at each box the aerosol size&density is different
    ! depending on the internally-mixed composition (win, 7/15/09)
    !=================================================================

    IF ( id_NK1 > 0 ) THEN
       CALL AERO_DIADEN( 1, Input_Opt, State_Chm, State_Grid, State_Met, &
                         SIZ_DIA, SIZ_DEN, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in call to "AERO_DIADEN"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF
#endif

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( CZ,      TEMPK,   TEMPC,   K,      VD                ) &
    !$OMP PRIVATE( LDT,     RSURFC,  C1,      XNU,    RT,      IOLSON   ) &
    !$OMP PRIVATE( II,      RI,      RLU,     RAC,    RGSS,    RGSO     ) &
    !$OMP PRIVATE( RCLS,    RCLO,    RAD0,    RIX,    GFACT,   GFACI    ) &
    !$OMP PRIVATE( RDC,     XMWH2O,  RIXX,    RLUXX,  RGSX,    RCLX     ) &
    !$OMP PRIVATE( DTMP1,   DTMP2,   DTMP3,   DTMP4,  VDS,     CZH      ) &
    !$OMP PRIVATE( CKUSTR,  REYNO,   CORR1,   CORR2,  Z0OBK,   RA       ) &
    !$OMP PRIVATE( Z0OBK_Alt, RA_Alt,DUMMY2_Alt,      DUMMY4_Alt        ) &
    !$OMP PRIVATE( DUMMY1,  DUMMY2,  DUMMY3,  DUMMY4, DAIR,    RB       ) &
    !$OMP PRIVATE( C1X,     VK,      I,       J,      IW                ) &
    !$OMP PRIVATE( DIAM,    DEN,     XLAI_FP, SUNCOS_FP,       CFRAC_FP ) &
    !$OMP PRIVATE( N_SPC,   alpha,   DEPVw                              ) &
#ifdef MODEL_GEOS
    !$OMP PRIVATE( RA2M,    Z0OBK2M, RA10M, Z0OBK10M                    ) &
#endif
#ifdef TOMAS
    !$OMP PRIVATE( BIN                                                  ) &
#endif
    !$OMP PRIVATE( SpcId,  SpcInfo                                      )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Initialize
       Z0OBK_Alt  = 0.0_f8
       DUMMY2_Alt = 0.0_f8
       DUMMY4_Alt = 0.0_f8
       RA_Alt     = 0.0_f8

       !** CZ is Altitude (m) at which deposition velocity is computed
       CZ = CZ1(I,J)

       !** TEMPK and TEMPC are surface air temperatures in K and in C
       TEMPK = TEMP(I,J)
       TEMPC = TEMP(I,J)-273.15e+0_f8

       !* Initialize variables
       DO K = 1,NUMDEP
          VD(K) = 0.0e+0_f8
          DO LDT = 1,NTYPE
             RSURFC(K,LDT) = 0.e+0_f8
          END DO
       END DO

       !** Calculate the kinematic viscosity XNU (m2 s-1) of air
       !** as a function of temperature.
       !** The kinematic viscosity is used to calculate the roughness heights
       !** over water surfaces and to diagnose whether such surfaces are
       !** aerodynamically rough or smooth using a Reynolds number criterion.
       !** The expression for the temperature dependence of XNU
       !** is from the FORTRAN code in Appendix II of Wesely [1988];
       !** I wasn't able to find an original reference but it seems benign enough.
       C1 = TEMPK/273.15e+0_f8
       XNU = 0.151e+0_f8*(C1**1.77e+0_f8)*1.0e-04_f8

       !* Compute bulk surface resistance for gases.
       !*
       !* Adjust external surface resistances for temperature;
       !* from Wesely [1989], expression given in text on p. 1296.
       !*
       !* BUG FIX!  Wesely [1989] gives RT = 1000.0*EXP(-TEMPC-4.0)
       !* so the inner parentheses are not needed (bmy, 3/4/99)
       !*        RT = 1000.0*EXP(-(TEMPC-4.0))
       RT = 1000.0e+0_f8*EXP(-TEMPC-4.0e+0_f8)
       !*
       !    Get surface resistances - loop over land types LDT
       !***********************************************************************
       !* The land types within each grid square are defined using the Olson
       !* land-type database.  Each of the Olson land types is assigned a
       !* corresponding "deposition land type" with characteristic values of
       !* surface resistance components.  There are 74 Olson land-types but only
       !* 11 deposition land-types (i.e., many of the Olson land types share the
       !* same deposition characteristics).  Surface resistance components for
       !* the "deposition land types" are from Wesely [1989] except for tropical
       !* forests [Jacob and Wofsy, 1990] and for tundra [Jacob et al., 1992].
       !* All surface resistance components are normalized to a leaf area index
       !* of unity.
       !*
       !* Olson land types, deposition land types, and surface resistance
       !* components are read from file 'drydep.table'; check that file for
       !* further details.
       !***********************************************************************

       ! Loop over the # of Olson land types in this grid box (I,J)
       DO 170 LDT = 1, IREG(I,J)

          ! If the land type is not represented in grid
          ! box  (I,J), then skip to the next land type
          IF ( IUSE(I,J,LDT) == 0 ) GOTO 170

          ! Olson land type index + 1
          IOLSON = ILAND(I,J,LDT)+1

          ! Dry deposition land type index
          II     = IDEP(IOLSON)
          !
          !** If the surface to be snow or ice;
          !** set II to 1 instead.
          !
          IF(LSNOW(I,J)) II=1

          !* Read the internal resistance RI (minimum stomatal resistance for
          !* water vapor,per unit area of leaf) from the IRI array; a '9999'
          !* value means no deposition to stomata so we impose a very large
          !* value for RI.
          RI(LDT) = DBLE(IRI(II))
          IF (RI(LDT)   .GE. 9999.e+0_f8) RI(LDT)   = 1.e+12_f8

          !** Cuticular resistances IRLU read in from 'drydep.table'
          !** are per unit area of leaf;
          !** divide them by the leaf area index to get a cuticular resistance
          !** for the bulk canopy.  If IRLU is '9999' it means there are no
          !** cuticular surfaces on which to deposit so we impose a very large
          !** value for RLU.
          IF ( IRLU(II) >= 9999 .or. XLAI(I,J,LDT) <= 0.e+0_f8 ) THEN
             RLU(LDT) = 1.e+6_f8
          ELSE
             RLU(LDT) = DBLE( IRLU(II) ) / XLAI(I,J,LDT)
             ! Additional resistance at low temperatures.
             ! Limit increase to a factor of 2.
             ! V. Shah 23 Oct 18
             ! Ref Jaegle et al. 2018
             RLU(LDT) = MIN( RLU(LDT) + RT, 2.e+0_f8 * RLU(LDT) )
          ENDIF
          !** The following are the remaining resistances for the Wesely
          !** resistance-in-series model for a surface canopy
          !** (see Atmos. Environ. paper, Fig.1).
          RAC(LDT)  = MAX(DBLE(IRAC(II)), 1.e+0_f8)
          IF (RAC(LDT)  .GE. 9999.e+0_f8) RAC(LDT)  = 1.e+12_f8
          RGSS(LDT) = MAX(DBLE(IRGSS(II)), 1.e+0_f8)
          ! Additional resistance at low temperatures.
          ! Limit increase of RGSS, RGSO, RCLS, & RCLO to a factor of 2.
          ! V. Shah 23 Oct 18
          ! Ref Jaegle et al. 2018
          RGSS(LDT) = MIN( RGSS(LDT) + RT, 2.e+0_f8 * RGSS(LDT) )
          IF (RGSS(LDT) .GE. 9999.e+0_f8) RGSS(LDT) = 1.e12_f8
          RGSO(LDT) = MAX(DBLE(IRGSO(II)) ,1.e+0_f8)
          RGSO(LDT) = MIN( RGSO(LDT) + RT, 2.e+0_f8 * RGSO(LDT) )
          IF (RGSO(LDT) .GE. 9999.e+0_f8) RGSO(LDT) = 1.e+12_f8
          RCLS(LDT) = DBLE(IRCLS(II))
          RCLS(LDT) = MIN( RCLS(LDT) + RT, 2.e+0_f8 * RCLS(LDT) )
          IF (RCLS(LDT) .GE. 9999.e+0_f8) RCLS(LDT) = 1.e+12_f8
          RCLO(LDT) = DBLE(IRCLO(II))
          RCLO(LDT) = MIN( RCLO(LDT) + RT, 2.e+0_f8 * RCLO(LDT) )
          IF (RCLO(LDT) .GE. 9999.e+0_f8) RCLO(LDT) = 1.e+12_f8
          !********************************************************************
          !*
          !*    Adjust stomatal resistances for insolation and temperature:
          !*
          !*     Temperature adjustment is from Wesely [1989], equation (3).
          !*
          !*     Light adjustment by the function BIOFIT is described by Wang
          !*     [1996]. It combines
          !*       - Local dependence of stomal resistance on the intensity I
          !*         of light impinging the leaf; this is expressed as a
          !*         mutliplicative factor I/(I+b) to the stomatal resistance
          !*         where b = 50 W m-2 (equation (7) of Baldocchi et al.[1987])
          !*       - radiative transfer of direct and diffuse radiation in the
          !*         canopy using equations (12)-(16) from Guenther et al.[1995]
          !*       - separate accounting of sunlit and shaded leaves using
          !*         equation (12) of Guenther et al. [1995]
          !*       - partitioning of the radiation at the top of the canopy into
          !*         direct and diffuse components using a parameterization to
          !*         results from an atmospheric radiative transfer model
          !*         [Wang, 1996]
          !*     The dependent variables of the function BIOFIT are the leaf
          !*     area index (XYLAI), the cosine of zenith angle (SUNCOS) and
          !*     the fractional cloud cover (CFRAC).  The factor GFACI
          !*     integrates the light dependence over the canopy depth; sp even
          !*     though RI is input per unit area of leaf it need not be scaled
          !*     by LAI to yield a bulk canopy value because that's already
          !*     done in the GFACI formulation.
          !********************************************************************
          RAD0 = RADIAT(I,J)
          RIX = RI(LDT)
          IF (RIX .GE. 9999.e+0_f8) GO TO 150
          GFACT = 100.0e+0_f8
          IF (TEMPC .GT. 0.e+0_f8 .AND. TEMPC .LT. 40.e+0_f8) &
              GFACT = 400.e+0_f8/TEMPC/(40.0e+0_f8-TEMPC)
          GFACI = 100.e+0_f8
          IF ( RAD0 > 0.e+0_f8 .and. XLAI(I,J,LDT) > 0.e+0_f8 ) THEN
             ! Now make sure all inputs to BIOFIT are flexible precision
             ! so that the code will compile properly (bmy, 12/18/14)
             XLAI_FP   = XLAI(I,J,LDT)
             SUNCOS_FP = SUNCOS(I,J)
             CFRAC_FP  = CFRAC(I,J)
             GFACI     = 1.e+0_f8 / BIOFIT( DRYCOEFF,  XLAI_FP, &
                                            SUNCOS_FP, CFRAC_FP, NN )
          ENDIF

          RIX = RIX*GFACT*GFACI

          ! Apply scaling factor to RIX when CO2 effect is turned on
          ! (ayhwong, 6/25/19)
          IF (Input_Opt%CO2_EFFECT) THEN
             RIX = RIX*Input_Opt%RS_SCALE
          ENDIF

150       CONTINUE
          !*
          !*    Compute aerodynamic resistance to lower elements in lower part
          !*    of the canopy or structure, assuming level terrain -
          !*    equation (5) of Wesely [1989].
          !*
          RDC = 100.e+0_f8*(1.0e+0_f8+1000.0e+0_f8/(RAD0+10.e+0_f8))
          !*
          !*    Loop over species; species-dependent corrections to resistances
          !*    are from equations (6)-(9) of Wesely [1989].
          !*
          DO 160  K = 1,NUMDEP

             !** exit for non-depositing species or aerosols.
             IF (.NOT. LDEP(K) .OR. AIROSOL(K)) GOTO 155

             ! Test for special treatment for O3 drydep to ocean
             N_SPC = State_Chm%Map_DryDep(K)
             IF ((N_SPC .EQ. ID_O3) .AND. (II .EQ. 11)) THEN

                IF (HCO_Salinity(I,J) .GT. 20.0_f8) THEN

                   ! Now apply the Luhar et al. [2018] equations for the
                   ! special treatment of O3 dry deposition to the ocean
                   ! surface 
                   CALL OCEANO3(State_Met%TSKIN(I,J),USTAR(I,J),DEPVw,I,J)

                   ! Now convert to the new rc value(s) can probably tidy
                   ! this up a bit
                   alpha = 10.0_f8**(-0.25-0.013*( &
                                     State_Met%TSKIN(I,J)-273.16_f8))
                   RSURFC(K,LDT) = 1.0_f8/(alpha*DEPVw)

                ELSE

                   ! It's not 'ocean' so we instead don't change it from
                   ! 'default' rc
                   RSURFC(K,LDT) = 2000.0_f8

                ENDIF

             ELSE

                !XMWH2O = 18.e-3_f8 ! Use global H2OMW (ewl, 1/6/16)
                XMWH2O = H2OMW * 1.e-3_f8
                RIXX = RIX*DIFFG(TEMPK,PRESSU(I,J),XMWH2O)/ &
                     DIFFG(TEMPK,PRESSU(I,J),XMW(K)) &
                     + 1.e+0_f8/(HSTAR(K)/3000.e+0_f8+100.e+0_f8*F0(K))
                RLUXX = 1.e+12_f8
                IF (RLU(LDT).LT.9999.e+0_f8) &
                     RLUXX = RLU(LDT)/(HSTAR(K)/1.0e+05_f8 + F0(K))

                ! If POPs simulation, scale cuticular resistances with octanol-
                ! air partition coefficient (Koa) instead of HSTAR 
                ! (clf, 1/3/2011)
                IF (IS_POPS) &
                     RLUXX = RLU(LDT)/(KOA(K)/1.0e+05_f8 + F0(K))

                !*
                !* To prevent virtually zero resistance to species with huge
                !* HSTAR, such as HNO3, a minimum value of RLUXX needs to be
                !* set. The rationality of the existence of such a minimum is
                !* demonstrated by the observed relationship between Vd(NOy-NOx)
                !* and Ustar in Munger et al.[1996];
                !* Vd(HNO3) never exceeds 2 cm s-1 in observations. The
                !* corresponding minimum resistance is 50 s m-1. This correction
                !* was introduced by J.Y. Liang on 7/9/95.
                !*
                RGSX = 1.e+0_f8/(HSTAR(K)/1.0e+05_f8/RGSS(LDT) + &
                       F0(K)/RGSO(LDT))
                RCLX = 1.e+0_f8/(HSTAR(K)/1.0e+05_f8/RCLS(LDT) + &
                       F0(K)/RCLO(LDT))
                !*
                !** Get the bulk surface resistance of the canopy, RSURFC, from
                !** the network of resistances in parallel and in series (Fig. 1
                !** of Wesely [1989])
                DTMP1=1.e+0_f8/RIXX
                DTMP2=1.e+0_f8/RLUXX
                DTMP3=1.e+0_f8/(RAC(LDT)+RGSX)
                DTMP4=1.e+0_f8/(RDC+RCLX)
                RSURFC(K,LDT) = 1.e+0_f8/(DTMP1 + DTMP2 + DTMP3 + DTMP4)

             ENDIF

             !** get surface deposition velocity for aerosols if needed;
             !** equations (15)-(17) of Walcek et al. [1986]
155          IF (.NOT. AIROSOL(K)) GOTO 160

             ! Get information about this species from the database
             SpcId   =  NTRAIND(K)
             SpcInfo => State_Chm%SpcData(SpcId)%Info

             IF ( SpcInfo%DD_AeroDryDep ) THEN

                !=====================================================
                ! Use size-resolved dry deposition calculations for
                ! seasalt aerosols.  We need to account for the
                ! hygroscopic growth of the aerosol particles.
                !=====================================================

                ! [Zhang et al., 2001]
                RSURFC(K,LDT) = AERO_SFCRSII( K,                   &
                                              II,                  &
                                              PRESSU(I,J)*1e-3_f8, &
                                              TEMPK,               &
                                              USTAR(I,J),          &
                                              RHB(I,J),            &
                                              W10(I,J),            &
                                              Input_Opt )

             ELSE IF ( SpcInfo%DD_DustDryDep ) THEN

                !=====================================================
                ! Use size-resolved dry deposition calculations for
                ! dust aerosols only.  Do not account for hygroscopic
                ! growth of the dust aerosol particles.
                !=====================================================  
#ifdef TOMAS
                !-------------------------------
                !%%% TOMAS SIMULATIONS %%%
                !-------------------------------
                IF ( TRIM(SpcInfo%Name) == 'DST1' .or. &
                     TRIM(SpcInfo%Name) == 'DST2' .or. &
                     TRIM(SpcInfo%Name) == 'DST3' .or. &
                     TRIM(SpcInfo%Name) == 'DST4' ) THEN

                   ! Particle diameter, convert [m] -> [um]
                   DIAM  = A_RADI(K) * 2.e+0_f8

                   ! Particle density [kg/m3]
                   DEN   = A_DEN(K)

                ELSE

                   ! Get current aerosol diameter and density
                   ! NOTE: In TOMAS the aerosol diameter and density
                   ! evolves with time.  We have to save these in the
                   ! DIAM and DEN variables so that we can hold these
                   ! PRIVATE for the parallel loop.
                   BIN  = NTRAIND(K) - id_NK1 + 1
                   DIAM = SIZ_DIA( I, J, BIN )     ! Diameter [m]
                   DEN  = SIZ_DEN( I, J, BIN )     ! Density  [kg/m3]

                ENDIF
#else
                !-------------------------------
                !%%% NON-TOMAS SIMULATIONS %%%
                !-------------------------------

                ! Particle diameter, convert [m] -> [um]
                DIAM  = A_RADI(K) * 2.e+0_f8

                ! Particle density [kg/m3]
                DEN   = A_DEN(K)
#endif

#ifdef APM
                IF(SpcInfo%Name(1:8)=='APMSPBIN')THEN
                   DIAM = RDRY(SpcId-APMIDS%id_SO4BIN1+1)* &
                        GFTOT3D(I,J,1,1)*2.D0
                   DEN = DENWET3D(I,J,1,1)*1.D3
                ENDIF
                IF(SpcInfo%Name(1:9)=='APMSEABIN')THEN
                   DIAM = RSALT(SpcId-APMIDS%id_SEABIN1+1)* &
                        GFTOT3D(I,J,1,2)*2.D0
                   DEN = DENWET3D(I,J,1,2)*1.D3
                ENDIF
                IF(SpcInfo%Name(1:9)=='APMDSTBIN')THEN
                   DIAM = RDST(SpcId-APMIDS%id_DSTBIN1+1)*2.D0
                   DEN = DENDST(SpcId-APMIDS%id_DSTBIN1+1)*1.D3
                ENDIF
                IF(SpcInfo%Name(1:8)=='APMLVSOA')THEN
                   DIAM = MWSIZE3D(I,J,1,1)*2.D0
                   DEN = DENWET3D(I,J,1,1)*1.D3
                ENDIF
                IF(SpcInfo%Name(1:8)=='APMCTSEA')THEN
                   DIAM = MWSIZE3D(I,J,1,2)*2.D0
                   DEN = DENWET3D(I,J,1,2)*1.D3
                ENDIF
                IF(SpcInfo%Name(1:8)=='APMCTDST')THEN
                   DIAM = MWSIZE3D(I,J,1,3)*2.D0
                   DEN = DENDST(12)*1.D3
                ENDIF
#endif

                ! [Zhang et al., 2001]
                RSURFC(K,LDT) = DUST_SFCRSII( K,                   &
                                              II,                  &
                                              PRESSU(I,J)*1e-3_f8, &
                                              TEMPK,               &
                                              USTAR(I,J),          &
                                              DIAM,                &
                                              DEN )

             ELSE

                !=====================================================
                ! Replace original code to statement 160 here: only
                ! do this for non-size-resolved species where
                ! AIROSOL(K)=T. (rjp, tdf, bec, bmy, 4/20/04)
                !=====================================================
                !! use Zhang et al for all aerosols (hotp 10/26/07)
                !VDS = 0.002D0*USTAR(I,J)
                !IF (OBK(I,J) .LT. 0.0D0) THEN
                !   VDS = VDS*(1.D0+(-300.D0/OBK(I,J))**0.6667D0)
                !ENDIF
                !
                !IF ( OBK(I,J) .EQ. 0.0D0 ) &
                !    WRITE(6,156) OBK(I,J),I,J,LDT
                ! 156 FORMAT(1X,'OBK(I,J)=',E11.2,1X,' I,J =',2I4, 1X,'LDT=',I3/)
                !CZH  = ZH(I,J)/OBK(I,J)
                !IF (CZH.LT.-30.0D0) VDS = 0.0009D0*USTAR(I,J)* &
                !   (-CZH)**0.6667D0

                RSURFC(K, LDT) = ADUST_SFCRSII(K, II, PRESSU(I,J)*1e-3_f8, &
                                            TEMPK, USTAR(I,J))

                !*
                !*    Set VDS to be less than VDSMAX (entry in input file
                !*    divided by 1.D4) VDSMAX is taken from Table 2 of Walcek
                !*    et al. [1986]. Invert to get corresponding R

                !RSURFC(K,LDT) = 1.D0/MIN(VDS, DBLE(IVSMAX(II))/1.D4)
             ENDIF

             ! Free pointer
             SpcInfo => NULL()

160       CONTINUE
             !*
170    CONTINUE
       !*
       !*    Set max and min values for bulk surface resistances
       !*
       DO 190 K = 1,NUMDEP
          IF (.NOT.LDEP(K)) GOTO 190
          DO 180 LDT = 1, IREG(I,J)
             IF ( IUSE(I,J,LDT) == 0 ) GOTO 180
             ! because of high resistance values, different rule applied for
             ! ocean ozone               
             N_SPC = State_Chm%Map_DryDep(K)
             IF ((N_SPC .EQ. ID_O3) .AND. (II .EQ. 11)) THEN
                RSURFC(K,LDT)= MAX(1.e+0_f8, MIN(RSURFC(K,LDT),999999.e+0_f8))
             ELSE
                RSURFC(K,LDT)= MAX(1.e+0_f8, MIN(RSURFC(K,LDT),9999.e+0_f8))
             ENDIF
             ! Set Rc for strong acids (HNO3,HCl,HBr) to 1 s/m
             ! V. Shah (23 Oct 18)
             ! Ref. Jaegle et al. 2018, cf. Erisman,van Pul,Ayers 1994
             IF ( HSTAR(K) .gt. 1.e+10_f8 ) THEN
                RSURFC(K,LDT)= 1.e+0_f8
             ENDIF
180       CONTINUE
190    CONTINUE
       !*
       !*    Loop through the different landuse types present in 
       !*    the grid square
       !*
       DO 500 LDT = 1, IREG(I,J)
          IF ( IUSE(I,J,LDT) == 0 ) GOTO 500
          IOLSON = ILAND(I,J,LDT) + 1

          !***** Get aerodynamic resistances Ra and Rb. ***********
          !   The aerodynamic resistance Ra is integrated from altitude z0+d up
          !   to the altitude z1 at which the dry deposition velocity is to be
          !   referenced. The integration corrects for stability using Monin-
          !   Obukhov similarity formulas from Businger et al. [1971] which
          !   apply over the range -2.5 < z/zMO < 1.5 (see their Figure 2).
          !   Under very unstable conditions when z1 > -2.5 zMO, we assume that
          !   there is no resistance to transfer in the convective column
          !   between zMO and z1. Under very stable conditions when z1 > 1.5 zMO
          !   we assume that vertical transfer in the column between zMO and z1
          !   is strongly suppressed so that the deposition velocity at altitude
          !   z1 is very low.  Under these conditions we just specify a very
          !   large Ra=1.E4 s m-1 (LRGERA = T).
          !**
          !   The Reynolds number REYNO diagnoses whether a surface is
          !   aerodynamically rough (REYNO > 1) or smooth.
          !
          !   NOTE: The criterion "REYNO > 1" was originally "REYNO > 10". See
          !   below for an explanation of why it was changed (hyl, 10/15/99)
          !
          !   Surface is rough in all cases except over water with low wind
          !   speeds. In the smooth case, vertical transport IN THE SUBLAYER
          !   near the surface is limited by molecular diffusion and is
          !   therefore very slow; we assign a large value we assign a large
          !   value of Ra + Rb to account for this effect.  [In Versions 3.2
          !   and earlier we used the formulation for Ra + Rb given in Equation
          !   (12) of Walcek et al [1986] to calculate the aerodynamic
          !   resistance over smooth surfaces.  However, that expression fails
          !   when u* is very small, as it yields negative values of Ra + Rb].
          !   (djj, hyl, bmy, 5/8/00)
          !**
          !   In the aerodynamically rough case, the expression for Ra is as
          !   given in equation (5) of Jacob et al. [1992]:
          !
          !          Ra = (1/ku*)*int(from z0 to z1) (phi(x)/z)dz
          !
          !   where x = (z-D)/zMO, z is the height above ground, and D is the
          !   displacement height which is typically 70-80% of the canopy
          !   height [Brutsaert, 1982].  We change the vertical coordinate so
          !   that z=0 at the displacement height; that's OK since for all
          !   practical applications z1 >> D.  In this manner we don't need
          !   to assume any specific value for the displacement height.
          !   Applying the variable transformation z -> x = z/zMO, the equation
          !   above becomes
          !
          !          Ra = (1/ku*)*int(from x0 to x1) (phi(x)/x)dx with x=z/zMO
          !
          !   Here phi is a stability correction function originally formulated
          !   by Businger et al. [1971] and given in eqns 5a and 5b of Jacob et
          !   al. [1992]. For unstable conditions,
          !
          !          phi(x) = a/sqrt(1-bx)  where a=0.74, b = 9
          !
          !   The analytical solution to the integral is [Dwight, 1957,
          !   integral 192.11]:
          !
          !          int(dx/(x*sqrt(1-bx))) = log(abs((sqrt(1-bx)-1)
          !                                   /(sqrt(1-bx)+1)))
          !
          !   which yields the expression for Ra used in the code for
          !   unstable conditions.  For stable conditions,
          !
          !          phi(x) = a + bx        where a=0.74, b = 4.7
          !
          !   and the analytical solution to the integral is
          !
          !          int((a/x)+b)dx = a*ln(x) + bx
          !
          !   which yields the expression of Ra used in the code for stable
          !   conditions.
          !**
          !   The formulation of RB for gases is equation (12) of Walcek et al.
          !   [1986].  The parameterization for deposition of aerosols does not
          !   include an RB term so RB for aerosols is set to zero.
          !   Modify phi(x) according to the non-local mixing scheme
          !   by Holtslag and Boville [1993] ( Lin, 07/18/08 )
          !   For unstable conditions,
          !          phi(x) = a/sqrt(1-bx)  where a=1.0, b=15.0
          !
          !   For stable conditions,
          !          phi(x) = a + bx
          !              where a=1.0, b=5.0 for 0 <= x <= 1, and
          !                    a=5.0, b=1.0 for x > 1.0
          !********************************************************

          CKUSTR = VON_KARMAN * USTAR(I,J)

          REYNO = USTAR(I,J)*ZO(I,J)/XNU

          IF ( OBK(I,J) .EQ. 0.0e+0_f8 ) &
               WRITE(6,211) OBK(I,J),I,J,LDT
211       FORMAT(1X,'OBK(I,J)=',E11.2,1X,' I,J = ',2I4,1X,'LDT=',I3/)
          CORR1 = CZ/OBK(I,J)

          ! Define Z0OBK
          Z0OBK = ZO(I,J)/OBK(I,J)
#ifdef MODEL_GEOS
          Z0OBK2M  = MAX(Z0OBK,  2e+0_fp/OBK(I,J) )
          Z0OBK10M = MAX(Z0OBK, 10e+0_fp/OBK(I,J) )
#endif
          !--------------------------------------------------------
          ! Z0OBK_Alt is Z0OBK for a user-specified height above
          ! the surface.  This is required for diagnostics.
          ! (krt, bmy, 7/10/19).
          !--------------------------------------------------------
          IF ( State_Diag%Archive_ConcAboveSfc ) THEN
             Z0OBK_Alt = Input_Opt%RA_Alt_Above_Sfc / OBK(I,J)
          ENDIF

          LRGERA(I,J) = .FALSE.
          ! Add option for non-local PBL (Lin, 03/31/09)
          IF (.NOT. Input_Opt%LNLPBL) THEN
             IF (CORR1 .GT. 0.e+0_f8) THEN
                IF (CORR1 .GT.  1.5e+0_f8) LRGERA(I,J) = .TRUE.
             ELSEIF(CORR1 .LE. 0.e+0_f8) THEN
                IF (CORR1 .LE. -2.5e+0_f8) CORR1 = -2.5e+0_f8
                CORR2 = LOG(-CORR1)
             ENDIF
          ENDIF

          IF (CKUSTR.EQ.0.0e+0_f8) THEN
             WRITE(6,212) I,J,CKUSTR,VON_KARMAN,USTAR(I,J)
212          FORMAT(1X,'I,J= ',2I4,1X,'CKUSTR=',E10.1,1X, &
                    'VON_KARMAN= ',E12.4,1X,'USTAR(I,J)= ',E12.4)
             CLOSE(98)
             STOP             ! debug
          ENDIF

          !...aerodynamically rough or smooth surface
          ! "In the classic study by Nikuradse (1933) the transition from
          ! smooth to rough was examined in pipe flow. He introduced a
          ! roughness Reynolds number Rr = U* Z0 / Nu and found the flow to
          ! be smooth for Rr < 0.13 and rough for Rr > 2.5 with a transition
          ! regime in between." (E.B. Kraus and J.A. Businger, Atmosphere-Ocean
          ! Interaction, second edition, P.144-145, 1994). 
          ! Similar statements can be found in the books: Evaporation into the
          ! atmosphere, by Wilfried Brutsaert ,P.59,89, 1982; or Seinfeld &
          ! Pandis, P.858, 1998.
          ! Here we assume a sudden transition point Rr = 1 from smooth to
          ! rough, following L. Merlivat (1978, The dependence of bulk
          ! evaporation coefficients on air-water interfacial conditions as
          ! determined by the isotopic method, J. Geophys. Res., Oceans &
          ! Atmos., 83, C6, 2977-2980). Also refer to Brutsaert's book, P.125.
          ! We used to use the criterion "REYNO > 10" for aerodynamically rough
          ! surface and now change to "REYNO > 1". (hyl, 10/15/99)
          !
          ! D. J. Jacob says to change the criterion for aerodynamically rough
          ! surface to REYNO > 0.1 (eck, djj, bmy, 11/17/05)
          IF ( REYNO < 0.1e+0_f8 ) GOTO 220

          ! Add option for non-local PBL (Lin, 03/31/09)
          IF (.NOT. Input_Opt%LNLPBL) THEN

             !...aerodynamically rough surface.
             !*
             IF (CORR1.LE.0.0e+0_f8 .AND. Z0OBK .LT. -1.e+0_f8)THEN
                !*... unstable condition; set RA to zero.
                !*    (first implemented in V. 3.2)
                RA     = 0.e+0_f8
#ifdef MODEL_GEOS
                RA2M   = 0.e+0_f8
                RA10M  = 0.e+0_f8
#endif

                !*... error trap: prevent CORR1 or Z0OBK from being
                !*... zero or close to zero (ckeller, 3/15/16)
             ELSEIF ( ABS(CORR1)<=SMALL .OR. ABS(Z0OBK)<=SMALL ) THEN
                RA = 0.e+0_f8
#ifdef MODEL_GEOS
                RA2M  = 0.e+0_f8
                RA10M = 0.e+0_f8
#endif

             ELSEIF (CORR1.LE.0.0e+0_f8 .AND. Z0OBK .GE. -1.e+0_f8) THEN
                !*... unstable conditions;
                !*... compute Ra as described above
                DUMMY1 = (1.e+0_f8 - 9e+0_f8*CORR1)**0.5e+0_f8
                DUMMY2 = (1.e+0_f8 - 9e+0_f8*Z0OBK)**0.5e+0_f8
                DUMMY3 = ABS((DUMMY1 - 1.e+0_f8)/(DUMMY1 + 1.e+0_f8))
                DUMMY4 = ABS((DUMMY2 - 1.e+0_f8)/(DUMMY2 + 1.e+0_f8))
                RA = 0.74e+0_f8* (1.e+0_f8/CKUSTR) * LOG(DUMMY3/DUMMY4)
#ifdef MODEL_GEOS
                ! 2M
                DUMMY1 = (1.e+0_f8 - 9e+0_f8*Z0OBK2M)**0.5e+0_f8
                DUMMY3 = ABS((DUMMY1 - 1.e+0_f8)/(DUMMY1 + 1.e+0_f8))
                RA2M   = 0.74e+0_f8* (1.e+0_f8/CKUSTR) * LOG(DUMMY3/DUMMY4)
                DUMMY1 = (1.e+0_f8 - 9e+0_f8*Z0OBK10M)**0.5e+0_f8
                DUMMY3 = ABS((DUMMY1 - 1.e+0_f8)/(DUMMY1 + 1.e+0_f8))
                RA10M  = 0.74e+0_f8* (1.e+0_f8/CKUSTR) * LOG(DUMMY3/DUMMY4)
#endif

             ELSEIF((CORR1.GT.0.0e+0_f8).AND.(.NOT.LRGERA(I,J)))  THEN
                !*... moderately stable conditions (z/zMO <1);
                !*... compute Ra as described above
                RA = (1e+0_f8/CKUSTR) * (.74e+0_f8*LOG(CORR1/Z0OBK) + &
                     4.7e+0_f8*(CORR1-Z0OBK))
#ifdef MODEL_GEOS
                RA2M = (1e+0_f8/CKUSTR) * (0.74_f8*LOG(Z0OBK2M/Z0OBK)+4.7_f8* &
                       (Z0OBK2M-Z0OBK))
                RA10M = (1e+0_f8/CKUSTR) * (0.74_f8*LOG(Z0OBK10M/Z0OBK)+4.7_f8*&
                        (Z0OBK10M-Z0OBK))
#endif

             ELSEIF(LRGERA(I,J)) THEN
                !*... very stable conditions
                RA     = 1.e+04_f8
#ifdef MODEL_GEOS
                RA2M   = 1.e+04_f8
                RA10M  = 1.e+04_f8
#endif
             ENDIF
             !* check that RA is positive; if RA is negative (as occassionally
             !* happened in version 3.1) send a warning message.

          ELSE

             IF (CORR1.LT.0.0e+0_f8) THEN
                !*... unstable conditions; compute Ra as described
                !*... above.
                !coef_a=1.e+0_f8
                !coef_b=15.e+0_f8
                DUMMY1 = (1.D0 - 15.e+0_f8*CORR1)**0.5e+0_f8
                DUMMY2 = (1.D0 - 15.e+0_f8*Z0OBK)**0.5e+0_f8
                DUMMY3 = ABS((DUMMY1 - 1.e+0_f8)/(DUMMY1 + 1.e+0_f8))
                DUMMY4 = ABS((DUMMY2 - 1.e+0_f8)/(DUMMY2 + 1.e+0_f8))
                RA = 1.e+0_f8 * (1.e+0_f8/CKUSTR) * LOG(DUMMY3/DUMMY4)
#ifdef MODEL_GEOS
                ! 2M
                DUMMY1 = (1.D0 - 15.e+0_f8*Z0OBK2M)**0.5e+0_f8
                DUMMY3 = ABS((DUMMY1 - 1.e+0_f8)/(DUMMY1 + 1.e+0_f8))
                RA2M = 1.e+0_f8 * (1.e+0_f8/CKUSTR)*LOG(DUMMY3/DUMMY4)
                ! 10M
                DUMMY1 = (1.D0 - 15.e+0_f8*Z0OBK10M)**0.5e+0_f8
                DUMMY3 = ABS((DUMMY1 - 1.e+0_f8)/(DUMMY1 + 1.e+0_f8))
                RA10M= 1.e+0_f8 * (1.e+0_f8/CKUSTR)*LOG(DUMMY3/DUMMY4)
#endif
                !--------------------------------------------------
                ! Compute RA at user-defined altitude above surface
                ! (krt, bmy, 7/10/19)
                !--------------------------------------------------
                IF ( State_Diag%Archive_ConcAboveSfc ) THEN
                   DUMMY2_Alt = SQRT( 1.0_f8 - 15.0_f8 * Z0OBK_Alt )
                   DUMMY4_Alt = ABS( ( DUMMY2_Alt - 1.0_f8 ) &
                                   / ( DUMMY2_Alt + 1.0_f8 ) )
                   RA_Alt     = ( 1.0_f8 / CKUSTR ) &
                                * LOG( DUMMY3 / DUMMY4_Alt )
                ENDIF

             ELSEIF((CORR1.GE.0.0e+0_f8).AND.(CORR1.LE.1.0e+0_f8)) THEN
                !coef_a=1.e+0_f8
                !coef_b=5.e+0_f8
                RA = (1D0/CKUSTR) * (1.e+0_f8*LOG(CORR1/Z0OBK) + &
                     5.e+0_f8*(CORR1-Z0OBK))
#ifdef MODEL_GEOS
                RA2M = (1D0/CKUSTR) * (1.e+0_f8*LOG(Z0OBK2M/Z0OBK)+ &
                        5.e+0_f8*(Z0OBK2M-Z0OBK))
                RA10M = (1D0/CKUSTR) * (1.e+0_f8*LOG(Z0OBK10M/Z0OBK)+ &
                        5.e+0_f8*(Z0OBK10M-Z0OBK))
#endif
                !--------------------------------------------------
                ! Compute RA at user-defined altitude above surface
                ! for diagnostic output (krt, bmy, 7/10/19)
                !--------------------------------------------------
                IF ( State_Diag%Archive_ConcAboveSfc ) THEN
                   RA_Alt = ( 1.0_fp / CKUSTR ) &
                          * ( 1.0_f8 * LOG( CORR1 / Z0OBK_Alt ) &
                            + 5.0_f8 * ( CORR1 - Z0OBK_Alt ))
                ENDIF

             ELSE ! CORR1 .GT. 1.0D0
                !coef_a=5e+0_f8
                !coef_b=1.e+0_f8
                RA = (1D0/CKUSTR) * (5.e+0_f8*LOG(CORR1/Z0OBK) + &
                     1.e+0_f8*(CORR1-Z0OBK))
#ifdef MODEL_GEOS
                RA2M = (1D0/CKUSTR) * (5.e+0_f8*LOG(Z0OBK2M/Z0OBK)+ &
                       1.e+0_f8*(Z0OBK2M-Z0OBK))
                RA10M = (1D0/CKUSTR) * (5.e+0_f8*LOG(Z0OBK10M/z0OBK)+ &
                        1.e+0_f8*(Z0OBK10M-Z0OBK))
#endif
                !--------------------------------------------------
                ! Compute RA at user-defined altitude above surface
                ! for diagnostic output (krt, bmy, 7/10/19)
                !--------------------------------------------------
                IF ( State_Diag%Archive_ConcAboveSfc ) THEN
                   RA_Alt = ( 1.0_f8 / CKUSTR ) &
                          * ( 5.0_f8 * LOG( CORR1 / Z0OBK_Alt ) &
                            + 1.0_f8 *( CORR1 - Z0OBK_Alt ) )
                ENDIF
             ENDIF

          ENDIF

          RA   = MIN(RA,1.e+4_f8)
#ifdef MODEL_GEOS
          RA2M   = MIN(RA2M,  1.e+4_f8)
          RA10M  = MIN(RA10M, 1.e+4_f8)
#endif
          ! If RA is < 0, set RA = 0 (bmy, 11/12/99)
          IF (RA .LT. 0.e+0_f8) THEN
             WRITE (6,1001) I,J,RA,CZ,ZO(I,J),OBK(I,J)
             RA = 0.0e+0_f8
          ENDIF
#ifdef MODEL_GEOS
          ! Adjust 2M Ra if needed
          IF ( RA2M  < 0.e+0_f8 ) RA2M  = 0.e+0_f8
          IF ( RA10M < 0.e+0_f8 ) RA10M = 0.e+0_f8
#endif
          !----------------------------------------------------
          ! Compute RA at user-defined altitude above surface
          ! for diagnostic output (krt, bmy, 7/10/19)
          !----------------------------------------------------
          IF ( State_Diag%Archive_ConcAboveSfc ) THEN
             RA_Alt = MIN( RA_Alt, 1.0e+4_f8 )

             ! Make sure RA_Alt is not negative
             IF ( RA_Alt <  0.0_f8 ) THEN
                IF ( Input_Opt%amIRoot ) THEN
                   WRITE (6,1001) I,J, RA_Alt, CZ, ZO(I,J) ,OBK(I,J)
                ENDIF
                RA_Alt = 0.0_f8
             ENDIF

             !-----------------------------------------------------
             ! HISTORY (aka netCDF diagnostics)
             !
             ! Store drydep resistance RA at user-defined
             ! altitude above surface into State_Diag
             !-----------------------------------------------------
             IF ( State_Diag%Archive_DryDepRaALT1 ) THEN
                State_Diag%DryDepRaALT1(I,J) = RA_Alt
             ENDIF

          ENDIF

1001      FORMAT('WARNING: RA < 0 IN SUBROUTINE DEPVEL',2I4,4(1X,E12.5))
          !* Get total resistance for deposition
          !* loop over species.
          DO 215 K = 1,NUMDEP
             IF (.NOT.LDEP(K)) GOTO 215
             !** DAIR is the thermal diffusivity of air;
             !** value of 0.2*1.E-4 m2 s-1 cited on p. 16,476 of 
             !** Jacob et al. [1992]
             DAIR = 0.2e0_f8*1.e-4_f8
             RB = (2.e+0_f8/CKUSTR)* &
                  (DAIR/DIFFG(TEMPK,PRESSU(I,J),XMW(K))) &
                  **0.667e+0_f8
             IF (AIROSOL(K)) RB=0.e+0_f8
             C1X(K) = RA + RB + RSURFC(K,LDT)
215       CONTINUE
          GOTO 240
220       CONTINUE
          !** ... aerodynamically smooth surface
          !** BUG FIX -- suppress drydep over smooth surfaces by setting Ra to
          !** a large value (1e4).  This prevents negative dry deposition
          !** velocities when u* is very small (djj, bmy, 5/8/00)
          DO 230 K = 1,NUMDEP
             IF ( LDEP(K) ) THEN
                RA     = 1.0D4
#ifdef MODEL_GEOS
                RA2M   = 1.0D4
                RA10M  = 1.0D4
#endif
                C1X(K) = RA + RSURFC(K,LDT)
             ENDIF
230       CONTINUE
240       CONTINUE
          !*
          !* IJUSE is the fraction of the grid square occupied by surface LDT
          !* in units of per mil (IJUSE=500 -> 50% of the grid square).  Add
          !* the contribution of surface type LDT to the deposition velocity;
          !* this is a loop over all surface types in the gridbox.
          !*
          DO 400 K = 1,NUMDEP
             IF (.NOT.LDEP(K)) GOTO 400
             VK(K) = VD(K)
             VD(K) = VK(K) +.001D0* DBLE( IUSE(I,J,LDT) )/C1X(K)
400       CONTINUE

#ifdef MODEL_GEOS
          !--- Eventually archive aerodynamic resistance.
          !--- Convert s m-1 to s cm-1.
          !--- Ra is set to an arbitrary large value of 1.0e4 in stable
          !--- conditions. Adjust archived Ra's downward to avoid excessive
          !--- surface concentration adjustments when using C'=(1-Ra*vd)*C. 
          !--- (ckeller, 2/2/18)
          !IF ( RA2M >= 1.0d4 ) RA2M = 0.0
          State_Chm%DryDepRa2m(I,J) = State_Chm%DryDepRa2m(I,J) + &
               0.001d0 * DBLE(IUSE(I,J,LDT)) * ( MAX(0.0d0,RA-RA2M)/100.0d0 )
          !IF ( RAKT >= 1.0d4 ) RAKT = 0.0
          State_Chm%DryDepRa10m(I,J) = State_Chm%DryDepRa10m(I,J) + &
               0.001d0 * DBLE(IUSE(I,J,LDT)) * ( MAX(0.0d0,RA-RA10M)/100.0d0 )
#endif

500    CONTINUE

       !** Load array DVEL
       DO 550 K=1,NUMDEP
          IF (.NOT.LDEP(K)) GOTO 550
          DVEL(I,J,K) = VD(K)

          ! Now check for negative deposition velocity before returning to
          ! calling program (bmy, 4/16/00)
          ! Also call CLEANUP to deallocate arrays (bmy, 10/15/02)
          IF ( DVEL(I,J,K) < 0e+0_f8 ) THEN
             !$OMP CRITICAL
             PRINT*, 'DEPVEL: Deposition velocity is negative!'
             PRINT*, 'Dep. Vel = ', DVEL(I,J,K)
             PRINT*, 'Species  = ', K
             PRINT*, 'I, J     = ', I, J
             PRINT*, 'RADIAT   = ', RADIAT(I,J)
             PRINT*, 'TEMP     = ', TEMP(I,J)
             PRINT*, 'SUNCOS   = ', SUNCOS(I,J)
             PRINT*, 'USTAR    = ', USTAR(I,J)
             PRINT*, 'CZ1      = ', CZ1(I,J)
             PRINT*, 'OBK      = ', OBK(I,J)
             PRINT*, 'CFRAC    = ', CFRAC(I,J)
             PRINT*, 'ZH       = ', ZH(I,J)
             PRINT*, 'LRGERA   = ', LRGERA(I,J)
             PRINT*, 'ZO       = ', ZO(I,J)
             PRINT*, 'STOP in depvel.f!'
             CALL CLEANUP
             STOP
             !$OMP END CRITICAL
          ENDIF

          ! Now check for IEEE NaN (not-a-number) condition before returning to
          ! calling program (bmy, 4/16/00)
          ! Also call CLEANUP to deallocate arrays (bmy, 10/15/02)
          IF ( IT_IS_NAN( DVEL(I,J,K) ) ) THEN
             !$OMP CRITICAL
             PRINT*, 'DEPVEL: Deposition velocity is NaN!'
             PRINT*, 'Dep. Vel = ', DVEL(I,J,K)
             PRINT*, 'Species  = ', K
             PRINT*, 'I, J     = ', I, J
             PRINT*, 'RADIAT   = ', RADIAT(I,J)
             PRINT*, 'TEMP     = ', TEMP(I,J)
             PRINT*, 'SUNCOS   = ', SUNCOS(I,J)
             PRINT*, 'USTAR    = ', USTAR(I,J)
             PRINT*, 'CZ1      = ', CZ1(I,J)
             PRINT*, 'OBK      = ', OBK(I,J)
             PRINT*, 'CFRAC    = ', CFRAC(I,J)
             PRINT*, 'ZH       = ', ZH(I,J)
             PRINT*, 'LRGERA   = ', LRGERA(I,J)
             PRINT*, 'ZO       = ', ZO(I,J)
             CALL CLEANUP
             STOP
             !$OMP END CRITICAL
          ENDIF
550    CONTINUE
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

#ifdef MODEL_GEOS
    ! Pass Monin Obukhov diagnostics
    IF ( ASSOCIATED ( State_Diag%MoninObukhov ) ) THEN
       State_Diag%MoninObukhov(:,:) = OBK(:,:)
    ENDIF
#endif

    ! Free pointers
    NULLIFY( IREG  )
    NULLIFY( ILAND )
    NULLIFY( IUSE  )
    NULLIFY( XLAI  )

  END SUBROUTINE DEPVEL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diffg
!
! !DESCRIPTION: Subroutine DIFFG calculates the molecular diffusivity [m2/s] in
!  air for a gas X of molecular weight XM [kg] at temperature TK [K] and
!  pressure PRESS [Pa].  (bmy, 5/16/06)
!\\
!\\
! !INTERFACE:
!
  FUNCTION DIFFG( TK, PRESS, XM ) RESULT( DIFF_G )
!
! !INPUT PARAMETERS:
!
    REAL(f8), INTENT(IN) :: TK     ! Temperature [K]
    REAL(f8), INTENT(IN) :: PRESS  ! Pressure [Pa]
    REAL(f8), INTENT(IN) :: XM     ! Molecular weight of gas [kg]
!
! !REMARKS:
!  We specify the molecular weight of air (XMAIR) and the hard-sphere molecular
!  radii of air (RADAIR) and of the diffusing gas (RADX).  The molecular
!  radius of air is given in a Table on p. 479 of Levine [1988].  The Table
!  also gives radii for some other molecules.  Rather than requesting the user
!  to supply a molecular radius we specify here a generic value of 2.E-10 m for
!  all molecules, which is good enough in terms of calculating the diffusivity
!  as long as molecule is not too big.
!
! !REVISION HISTORY:
!  16 May 2006 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(f8)             :: AIRDEN, Z, DIAM, FRPATH, SPEED, DIFF_G
!
! !DEFINED PARAMETERS:
!
    REAL(f8), PARAMETER  :: XMAIR  = 28.8e-3_f8 ! Moist air molec wt? (ewl)
    REAL(f8), PARAMETER  :: RADAIR = 1.2e-10_f8
    REAL(f8), PARAMETER  :: RADX   = 1.5e-10_f8

    !=================================================================
    ! DIFFG begins here!
    !=================================================================

    ! Air density [molec/m3]
    AIRDEN = ( PRESS * AVO ) / ( RSTARG * TK )

    ! DIAM is the collision diameter for gas X with air.
    DIAM   = RADX + RADAIR

    ! Calculate the mean free path for gas X in air:
    ! eq. 8.5 of Seinfeld [1986];
    Z      = XM  / XMAIR
    FRPATH = 1e+0_f8 /( PI * SQRT( 1e+0_f8 + Z ) * AIRDEN * ( DIAM**2 ) )

    ! Calculate average speed of gas X; eq. 15.47 of Levine [1988]
    SPEED  = SQRT( 8e+0_f8 * RSTARG * TK / ( PI * XM ) )

    ! Calculate diffusion coefficient of gas X in air;
    ! eq. 8.9 of Seinfeld [1986]
    DIFF_G = ( 3e+0_f8 * PI / 32e+0_f8 ) * ( 1e+0_f8 + Z ) * FRPATH * SPEED

  END FUNCTION DIFFG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_drydep_inputs
!
! !DESCRIPTION: Subroutine READ\_DRYDEP\_INPUTS reads inputs for the dry
!  deposition module corresponding to either the Olson 1992 (GEOS-Chem default)
!  or Olson 2001 (planned replacement for Olson 1992) land map.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_DRYDEP_INPUTS( Input_Opt, &
                                 DRYCOEFF,  IOLSON,   IDEP,   &
                                 IWATER,    NWATER,   IZO,    &
                                 IDRYDEP,   IRI,      IRLU,   &
                                 IRAC,      IRGSS,    IRGSO,  &
                                 IRCLS,     IRCLO,    IVSMAX, &
                                 RC )
!
! !USES:
!
    USE CMN_SIZE_Mod,  ONLY : NPOLY, NSURFTYPE
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput

    ! Modules for netCDF read
    USE m_netcdf_io_open
    USE m_netcdf_io_get_dimlen
    USE m_netcdf_io_read
    USE m_netcdf_io_readattr
    USE m_netcdf_io_close

#     include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    !-----------------------------------------------------------------------
    ! DRYCOEFF : Baldocchi polynomial coeffs
    ! IOLSON   : Olson land type indices (+1)
    ! IDEP     : Mapping: Olson ==> drydep ID
    ! IWATER   : Olson types that represent water
    ! NWATER   : Number of Olson types that are water
    ! IZO      : Default Z0 (routgness height) for each Olson land type
    ! IDRYDEP  : Dry deposition land type indices
    ! IRI      : RI   resistance for drydep
    ! IRLU     : RLU  resistance for drydep
    ! IRAC     : RAC  resistance for drydep
    ! IRGSS    : RGSS resistance for drydep
    ! IRGSO    : RGSO resistance for drydep
    ! IRCLS    : RCLS resistance for drydep
    ! IRCLO    : RCLO resistance for drydep
    ! IVSMAX   : Max drydep velocity (for aerosol) perr drydep land type
    !-----------------------------------------------------------------------
    REAL(fp), INTENT(OUT) :: DRYCOEFF(NPOLY    )
    INTEGER,  INTENT(OUT) :: IOLSON  (NSURFTYPE)
    INTEGER,  INTENT(OUT) :: IDEP    (NSURFTYPE)
    INTEGER,  INTENT(OUT) :: IWATER  (NSURFTYPE)
    INTEGER,  INTENT(OUT) :: NWATER
    INTEGER,  INTENT(OUT) :: IZO     (NSURFTYPE)
    INTEGER,  INTENT(OUT) :: IDRYDEP (NDRYDTYPE)
    INTEGER,  INTENT(OUT) :: IRI     (NDRYDTYPE)
    INTEGER,  INTENT(OUT) :: IRLU    (NDRYDTYPE)
    INTEGER,  INTENT(OUT) :: IRAC    (NDRYDTYPE)
    INTEGER,  INTENT(OUT) :: IRGSS   (NDRYDTYPE)
    INTEGER,  INTENT(OUT) :: IRGSO   (NDRYDTYPE)
    INTEGER,  INTENT(OUT) :: IRCLS   (NDRYDTYPE)
    INTEGER,  INTENT(OUT) :: IRCLO   (NDRYDTYPE)
    INTEGER,  INTENT(OUT) :: IVSMAX  (NDRYDTYPE)

    ! Success or failure flag
    INTEGER,  INTENT(OUT) :: RC
!
! !REMARKS:
!  Routine READ_DRYDEP_INPUTS replaces routines MODIN (which read the ASCII
!  file "drydep.table") and RDDRYCF (which read the ASCII file "drydep.coef").
!                                                                             .
!  READ_DRYDEP_INPUTS was generated from the Perl script "ncCodeRead", which
!  is part of the NcdfUtilities package (with subsequent hand-editing).
!                                                                             .
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!
! !REVISION HISTORY:
!  26 Mar 2012 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: FileExists
    INTEGER            :: fId                 ! netCDF file ID

    ! Strings
    CHARACTER(LEN=255) :: nc_dir              ! netCDF directory name
    CHARACTER(LEN=255) :: nc_file             ! netCDF file name
    CHARACTER(LEN=255) :: nc_path             ! netCDF path name
    CHARACTER(LEN=255) :: v_name              ! netCDF variable name
    CHARACTER(LEN=255) :: a_name              ! netCDF attribute name
    CHARACTER(LEN=255) :: a_val               ! netCDF attribute value
    CHARACTER(LEN=255) :: FileMsg
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    ! Arrays for netCDF start and count values
    INTEGER            :: st1d(1), ct1d(1)    ! For 1D arrays

    ! Shadow variable for reading in data at REAL*8
    REAL(f8)           :: DRYCOEFF_R8(NPOLY)

    !=================================================================
    ! In dry-run mode, print file path to dryrun log and exit.
    ! Otherwise, print file path to stdout and continue.
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at READ_DRYDEP_INPUTS (in module GeosCore/drydep_mod.F90)'

    ! Input file for Olson 2001
    nc_file = 'Olson_2001_Drydep_Inputs.nc'

    ! Open file for reading (construct full data path)
    nc_dir  = TRIM( Input_Opt%CHEM_INPUTS_DIR ) // &
              'Olson_Land_Map_201203/'
    nc_path = TRIM( nc_dir ) // TRIM( nc_file )

    ! Test if the file exists
    INQUIRE( FILE=TRIM( nc_path ), EXIST=FileExists )

    ! Test if the file exists and define an output string
    IF ( FileExists ) THEN
       FileMsg = 'READ_DRYDEP_INPUTS: Opening'
    ELSE
       FileMsg = 'READ_DRYDEP_INPUTS: REQUIRED FILE NOT FOUND'
    ENDIF

    ! Write to stdout for both regular and dry-run simulations
    IF ( Input_Opt%AmIRoot ) THEN
       WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( nc_path )
300    FORMAT( a, ' ', a )
    ENDIF

    ! For dry-run simulations, return to calling program.
    ! For regular simulations, throw an error if we can't find the file.
    IF ( Input_Opt%DryRun ) THEN
       RETURN
    ELSE
       IF ( .not. FileExists ) THEN
          WRITE( ErrMsg, 300 ) TRIM( FileMsg ), TRIM( nc_path )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! Open and read data from the netCDF file
    !=================================================================       
    CALL Ncop_Rd( fId, TRIM(nc_path) )

    !----------------------------------------
    ! VARIABLE: DRYCOEFF
    !----------------------------------------

    ! Variable name
    v_name = "DRYCOEFF"

    ! Read DRYCOEFF from file
    st1d   = (/ 1     /)
    ct1d   = (/ NPOLY /)
    CALL NcRd( DRYCOEFF_R8, fId, TRIM(v_name), st1d, ct1d )

    ! Read the DRYCOEFF:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
130    FORMAT( '%% Successfully read ',       a, ' [', a, ']' )
    ENDIF

    ! Cast from REAL*8 to flexible precision
    DRYCOEFF = DRYCOEFF_R8

    !----------------------------------------
    ! VARIABLE: IOLSON
    !----------------------------------------

    ! Variable name
    v_name = "IOLSON"

    ! Read IOLSON from file
    st1d   = (/ 1         /)
    ct1d   = (/ NSURFTYPE /)
    CALL NcRd( IOLSON, fId, TRIM(v_name), st1d, ct1d )

    ! Read the IOLSON:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF

    !----------------------------------------
    ! VARIABLE: IDEP
    !----------------------------------------

    ! Variable name
    v_name = "IDEP"

    ! Read IDEP from file
    st1d   = (/ 1         /)
    ct1d   = (/ NSURFTYPE /)
    CALL NcRd( IDEP, fId, TRIM(v_name), st1d, ct1d )

    ! Read the IDEP:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF

    !----------------------------------------
    ! VARIABLE: IWATER
    !----------------------------------------

    ! Variable name
    v_name = "IWATER"

    ! Get the # of Olson types that are water
    ! (NOTE: IWATER is an index array, dimension name = variable name)
    CALL NcGet_DimLen( fId, TRIM(v_name), NWATER )

    ! Initialize
    IWATER = 0

    ! Read IWATER from file
    ! NOTE: IWATER is declared with NNSURFTYPE, but has NWATER values
    ! The rest can be zeroed out
    st1d   = (/ 1      /)
    ct1d   = (/ NWATER /)
    CALL NcRd( IWATER(1:NWATER), fId, TRIM(v_name), st1d, ct1d )

    ! Read the IWATER:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF

    !----------------------------------------
    ! VARIABLE: IZO
    !----------------------------------------

    ! Variable name
    v_name = "IZO"

    ! Read IZO from file
    st1d   = (/ 1         /)
    ct1d   = (/ NSURFTYPE /)
    CALL NcRd( IZO, fId, TRIM(v_name), st1d, ct1d )

    ! Read the IZO:long_name attribute
    a_name = "long_name"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF

    !----------------------------------------
    ! VARIABLE: IDRYDEP
    !----------------------------------------

    ! Variable name
    v_name = "IDRYDEP"

    ! Read IDRYDEP from file
    st1d   = (/ 1         /)
    ct1d   = (/ NDRYDTYPE /)
    CALL NcRd( IDRYDEP, fId, TRIM(v_name), st1d, ct1d )

    ! Read the IDRYDEP:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF

    !----------------------------------------
    ! VARIABLE: IRI
    !----------------------------------------

    ! Variable name
    v_name = "IRI"

    ! Read IRI from file
    st1d   = (/ 1         /)
    ct1d   = (/ NDRYDTYPE /)
    CALL NcRd( IRI, fId, TRIM(v_name), st1d, ct1d )

    ! Read the IRI:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! For Olson 2001 land map, change IRI for coniferous forests
    ! to match IRI for deciduous forests (skim, mps, 2/3/14)
    IRI(3) = 200

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF

    !----------------------------------------
    ! VARIABLE: IRLU
    !----------------------------------------

    ! Variable name
    v_name = "IRLU"

    ! Read IRLU from file
    st1d   = (/ 1         /)
    ct1d   = (/ NDRYDTYPE /)
    CALL NcRd( IRLU, fId, TRIM(v_name), st1d, ct1d )

    ! Read the IRLU:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF

    !----------------------------------------
    ! VARIABLE: IRAC
    !----------------------------------------

    ! Variable name
    v_name = "IRAC"

    ! Read IRAC from file
    st1d   = (/ 1         /)
    ct1d   = (/ NDRYDTYPE /)
    CALL NcRd( IRAC, fId, TRIM(v_name), st1d, ct1d )

    ! Read the IRAC:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF

    !----------------------------------------
    ! VARIABLE: IRGSS
    !----------------------------------------

    ! Variable name
    v_name = "IRGSS"

    ! Read IRGSS from file
    st1d   = (/ 1         /)
    ct1d   = (/ NDRYDTYPE /)
    CALL NcRd( IRGSS, fId, TRIM(v_name), st1d, ct1d )

    ! Read the IRGSS:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF

    !----------------------------------------
    ! VARIABLE: IRGSO
    !----------------------------------------

    ! Variable name
    v_name = "IRGSO"

    ! Read IRGSO from file
    st1d   = (/ 1         /)
    ct1d   = (/ NDRYDTYPE /)
    CALL NcRd( IRGSO, fId, TRIM(v_name), st1d, ct1d )

    ! Read the IRGSO:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF

    !----------------------------------------
    ! VARIABLE: IRCLS
    !----------------------------------------

    ! Variable name
    v_name = "IRCLS"

    ! Read IRCLS from file
    st1d   = (/ 1         /)
    ct1d   = (/ NDRYDTYPE /)
    CALL NcRd( IRCLS, fId, TRIM(v_name), st1d, ct1d )

    ! Read the IRCLS:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF

    !----------------------------------------
    ! VARIABLE: IRCLO
    !----------------------------------------

    ! Variable name
    v_name = "IRCLO"

    ! Read IRCLO from file
    st1d   = (/ 1         /)
    ct1d   = (/ NDRYDTYPE /)
    CALL NcRd( IRCLO, fId, TRIM(v_name), st1d, ct1d )

    ! Read the IRCLO:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF

    !----------------------------------------
    ! VARIABLE: IVSMAX
    !----------------------------------------

    ! Variable name
    v_name = "IVSMAX"

    ! Read IVSMAX from file
    st1d   = (/ 1         /)
    ct1d   = (/ NDRYDTYPE /)
    CALL NcRd( IVSMAX, fId, TRIM(v_name), st1d, ct1d )

    ! Read the IVSMAX:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF

    !=================================================================
    ! Cleanup and quit
    !=================================================================

    ! Close netCDF file
    CALL NcCl( fId )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(a)' ) REPEAT( '%', 79 )
    ENDIF

  END SUBROUTINE READ_DRYDEP_INPUTS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aero_sfcrsii
!
! !DESCRIPTION: Function AERO\_SFCRSII computes the aerodynamic resistance of
!  seasalt aerosol species according to Zhang et al 2001.  We account for
!  hygroscopic growth of the seasalt aerosol particles. (rjp, tdf, bec, bmy,
!  4/1/04, 6/11/08)
!\\
!\\
! !INTERFACE:
!
  FUNCTION AERO_SFCRSII( K, II, PRESS, TEMP, USTAR, RHB, W10, Input_Opt ) &
       RESULT(RS)
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN) :: K     ! Drydep species index (range: 1-NUMDEP)
    INTEGER,        INTENT(IN) :: II    ! Surface type index of GEOS-CHEM
    REAL(f8),       INTENT(IN) :: PRESS ! Pressure [kPa] (1 mb=100 Pa=0.1 kPa)
    REAL(f8),       INTENT(IN) :: TEMP  ! Temperature [K]
    REAL(f8),       INTENT(IN) :: USTAR ! Friction velocity [m/s]
    REAL(f8),       INTENT(IN) :: RHB   ! Relative humidity (fraction)
    REAL(f8),       INTENT(IN) :: W10   ! 10 m windspeed [m/s]
    TYPE(OptInput), INTENT(IN) :: Input_Opt ! Input Options object
!
! !RETURN VALUE:
!
    REAL(f8)                   :: RS    ! Surface resistance for particles [s/m]
!
! !REMARKS:
!  Do computations internally with REAL*8 (8-byte) floating-point precision,
!  in order to avoid a loss of precision.
!
! !REVISION HISTORY:
!  01 Apr 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: N
    REAL(f8), PARAMETER   :: C1 = 0.7674e+0_f8
    REAL(f8), PARAMETER   :: C2 = 3.079e+0_f8
    REAL(f8), PARAMETER   :: C3 = 2.573e-11_f8
    REAL(f8), PARAMETER   :: C4 = -1.424e+0_f8
    REAL(f8), PARAMETER   :: BETA  = 2.e+0_f8
    REAL(f8), PARAMETER   :: E0 = 3.e+0_f8
    REAL(f8)  :: AIRVS       ! kinematic viscosity of Air (m^2/s)
    REAL(f8)  :: DP          ! Diameter of aerosol [um]
    REAL(f8)  :: PDP         ! Press * Dp
    REAL(f8)  :: CONST       ! Constant for settling velocity calculations
    REAL(f8)  :: SLIP        ! Slip correction factor
    REAL(f8)  :: VISC        ! Viscosity of air (Pa s)
    REAL(f8)  :: DIFF        ! Brownian Diffusion constant for particles (m2/s)
    REAL(f8)  :: SC, ST      ! Schmidt and Stokes number (nondim)
    REAL(f8)  :: RHBL        ! Relative humidity local

    ! replace RCM with RUM (radius in microns instead of cm) - jaegle 5/11/11
    !REAL(f8)  :: DIAM, DEN, RATIO_R, RWET, RCM
    REAL(f8)  :: DIAM, DEN, RATIO_R, RWET, RUM
    REAL(f8)  :: FAC1, FAC2
    REAL(f8)  :: EB, EIM, EIN, R1, AA, VTS
    ! New variables added (jaegle 5/11/11)
    REAL(f8)  :: SW
    REAL(f8)  :: SALT_MASS, SALT_MASS_TOTAL, VTS_WEIGHT, DMIDW ! for weighting the settling velocity
    REAL(f8)  :: D0, D1  !lower and upper bounds of sea-salt dry diameter bins
    REAL(f8)  :: DEDGE
    REAL(f8)  :: DEN1, WTP
    INTEGER   :: ID,NR
    LOGICAL, SAVE          :: FIRST = .TRUE.

    !increment of radius for integration of settling velocity (um)
    REAL(f8), PARAMETER      :: DR    = 5.e-2_f8

    ! Parameters for polynomial coefficients to derive seawater
    ! density. From Tang et al. (1997) - jaegle 5/11/11
    REAL(f8),  PARAMETER     :: A1 =  7.93e-3_f8
    REAL(f8),  PARAMETER     :: A2 = -4.28e-5_f8
    REAL(f8),  PARAMETER     :: A3 =  2.52e-6_f8
    REAL(f8),  PARAMETER     :: A4 = -2.35e-8_f8
    REAL(f8),  PARAMETER     :: EPSI = 1.0e-4_f8

    ! parameters for assumed size distribution of accumulation and coarse
    ! mode sea salt aerosols, as described in Jaegle et al. (ACP, 11, 2011)
    ! (jaegle, 5/11/11)
    ! 1) geometric dry mean diameters (microns)
    REAL(f8),  PARAMETER     ::   RG_A = 0.085e+0_f8
    REAL(f8),  PARAMETER     ::   RG_C = 0.4e+0_f8
    ! 2) sigma of the size distribution
    REAL(f8),  PARAMETER     ::   SIG_A = 1.5e+0_f8
    REAL(f8),  PARAMETER     ::   SIG_C = 1.8e+0_f8

    !=======================================================================
    !   #  LUC [Zhang et al., 2001]                GEOS-CHEM LUC (Corr. #)
    !-----------------------------------------------------------------------
    !   1 - Evergreen needleleaf trees             Snow/Ice          (12)
    !   2 - Evergreen broadleaf trees              Deciduous forest  ( 4)
    !   3 - Deciduous needleleaf trees             Coniferous forest ( 1)
    !   4 - Deciduous broadleaf trees              Agricultural land ( 7)
    !   5 - Mixed broadleaf and needleleaf trees   Shrub/grassland   (10)
    !   6 - Grass                                  Amazon forest     ( 2)
    !   7 - Crops and mixed farming                Tundra            ( 9)
    !   8 - Desert                                 Desert            ( 8)
    !   9 - Tundra                                 Wetland           (11)
    !  10 - Shrubs and interrupted woodlands       Urban             (15)
    !  11 - Wet land with plants                   Water             (14)
    !  12 - Ice cap and glacier
    !  13 - Inland water
    !  14 - Ocean
    !  15 - Urban
    !=======================================================================
    ! GEOS-CHEM LUC              1, 2, 3, 4, 5, 6, 7  8, 9,10,11
    INTEGER :: LUCINDEX(11) = (/12, 4, 1, 7,10, 2, 9, 8,11,15,14/)
    INTEGER :: LUC

    !=================================================================
    !   LUC       1,    2,    3,    4,    5,    6,    7,    8,
    !   alpha   1.0,  0.6,  1.1,  0.8,  0.8,  1.2,  1.2, 50.0,
    !   gamma  0.56, 0.58, 0.56, 0.56, 0.56, 0.54, 0.54, 0.54
    !
    !   LUC       9,   10,   11,   12,   13,   14,   15
    !   alpha  50.0,  1,3,  2.0, 50.0,100.0,100.0,  1.5
    !   gamma  0.54, 0.54, 0.54, 0.54, 0.50, 0.50, 0.56
    !=================================================================
    REAL(f8)  :: ALPHA(15) = (/   1.0e+0_f8,   0.6e+0_f8,  1.1e+0_f8, &
                                  0.8e+0_f8,   0.8e+0_f8,  1.2e+0_f8, &
                                  1.2e+0_f8,  50.0e+0_f8, 50.0e+0_f8, &
                                  1.3e+0_f8,   2.0e+0_f8, 50.0e+0_f8, &
                                100.0e+0_f8, 100.0e+0_f8,  1.5e+0_f8  /)

    REAL(f8)  :: GAMMA(15) = (/ 0.56e+0_f8, 0.58e+0_f8, 0.56e+0_f8, &
                                0.56e+0_f8, 0.56e+0_f8, 0.54e+0_f8, &
                                0.54e+0_f8, 0.54e+0_f8, 0.54e+0_f8, &
                                0.54e+0_f8, 0.54e+0_f8, 0.54e+0_f8, &
                                0.50e+0_f8, 0.50e+0_f8, 0.56e+0_f8  /)

    !...A unit is (mm) so multiply by 1.D-3 to (m)
    !   LUC       1,    2,    3,    4,    5,    6,    7,    8,
    !   SC1     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,
    !   SC2     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,
    ! A SC3     2.0,  5.0,  5.0, 10.0,  5.0,  5.0,  5.0,-999.,
    !   SC4     2.0,  5.0,  5.0, 10.0,  5.0,  5.0,  5.0,-999.,
    !   SC5     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,
    !
    !   LUC       9,   10,   11,   12,   13,   14,   15
    !   SC1   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
    !   SC2   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
    ! A SC3   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
    !   SC4   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
    !   SC5   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
    REAL(f8)  :: A(15,5)

    REAL(f8)  :: Aavg(15)

    DATA   A /  2.0e+0_f8,   5.0e+0_f8,   2.0e+0_f8,   5.0e+0_f8,  5.0e+0_f8, &
                2.0e+0_f8,   2.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
               10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
               
                2.0e+0_f8,   5.0e+0_f8,   2.0e+0_f8,   5.0e+0_f8,  5.0e+0_f8, &
                2.0e+0_f8,   2.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
               10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
               
                2.0e+0_f8,   5.0e+0_f8,   5.0e+0_f8,  10.0e+0_f8,  5.0e+0_f8, &
                5.0e+0_f8,   5.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
               10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
               
                2.0e+0_f8,   5.0e+0_f8,   5.0e+0_f8,  10.0e+0_f8,  5.0e+0_f8, &
                5.0e+0_f8,   5.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
               10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
               
                2.0e+0_f8,   5.0e+0_f8,   2.0e+0_f8,   5.0e+0_f8,  5.0e+0_f8, &
                2.0e+0_f8,   2.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
               10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8  /

    ! Annual average of A
    Aavg(:) = (A(:,1)+A(:,2)+A(:,3)+A(:,4)+A(:,5))/5.
    LUC     = LUCINDEX(II)
    AA      = Aavg(LUC) * 1.e-3_f8

    !=================================================================
    !...Ref. Zhang et al., AE 35(2001) 549-560
    !.
    !...Model theroy
    !    Vd = Vs + 1./(Ra+Rs)
    !      where Vs is the gravitational settling velocity,
    !      Ra is the aerodynamic resistance above the canopy
    !      Rs  is the surface resistance
    !    Here we calculate Rs only..
    !    Rs = 1 / (Eo*Ustar*(Eb+Eim+Ein)*R1)
    !      where Eo is an empirical constant ( = 3.)
    !      Ustar is the friction velocity
    !      Collection efficiency from
    !        Eb,  [Brownian diffusion]
    !        Eim, [Impaction]
    !        Ein, [Interception]
    !      R1 is the correction factor representing the fraction
    !         of particles that stick to the surface.
    !=======================================================================
    !      Eb is a funciont of Schmidt number, Eb = Sc^(-gamma)
    !         Sc = v/D, v (the kinematic viscosity of air)
    !                   D (particle brownian diffusivity)
    !         r usually lies between 1/2 and 2/3
    !      Eim is a function of Stokes number, St
    !          St = Vs * Ustar / (g0 * A)   for vegetated surfaces
    !          St = Vs * Ustar * Ustar / v  for smooth surface
    !          A is the characteristic radius of collectors.
    !
    !       1) Slinn (1982)
    !           Eim = 10^(-3/St)          for smooth surface
    !           Eim = St^2 / ( 1 + St^2 ) for vegetative canopies
    !       2) Peters and Eiden (1992)
    !           Eim = ( St / ( alpha + St ) )^(beta)
    !                alpha(=0.8) and beta(=2) are constants
    !       3) Giorgi (1986)
    !           Eim = St^2 / ( 400 + St^2 )     for smooth surface
    !           Eim = ( St / (0.6 + St) )^(3.2) for vegetative surface
    !       4) Davidson et al.(1982)
    !           Eim = St^3 / (St^3+0.753*St^2+2.796St-0.202) for grassland
    !       5) Zhang et al.(2001) used 2) method with alpha varying with
    !          vegetation type and beta equal to 2
    !
    !      Ein = 0.5 * ( Dp / A )^2
    !
    !      R1 (Particle rebound)  = exp(-St^0.5)
    !=================================================================
    !      Update (jaegle 5/11/2011): The above formulation of Zhang et al
    !      (2001) is valid for land surfaces and was originally based on the
    !      work of Slinn (1982). Over water surfaces, the work of reference
    !      is that of Slinn and Slinn (1980) who use the term "viscous
    !      sublayer" to refer to the thin layer extending 0.1-1mm above the
    !      water surface. Due to the proximity of the water, the RH in this
    !      layer is much higher than the ambient RH in the surface layer.
    !      According to Lewis and Schwartz (2004): "Relative humidities of
    !      99% and 100% were considered by Slinn and Slinn for the viscous
    !      sublayer, however near the ocean surface RH would be limited to
    !      near 98% because of the vapor pressure lowering of water over
    !      seawater due to the salt content". We will thus use a constant
    !      value RH=98% over all ocean boxes. This affects the growth of
    !      particles (the wet radius at RH=98% is x4 the dry radius) and thus
    !      affects all the terms depending on particle size.
    !
    !      Other updates for ocean surfaces:
    !         a)   Over ocean surfaces the formulation from Slinn & Slinn for
    !              the resistance in the viscous layer is
    !                Rs = 1 / (Cd/VON_KARMAN*U10m*(Eb+Eim)+VTS)
    !              with  Cd=(Ustar/U10m)**2, and VTS is the gravitational
    !              settling in the viscous layer. Note that the gravitational
    !              settling calculated here for the viscous layer is >> than
    !              the one calculated for the surface layer in seasalt_mod.f
    !              because of the higher RH.
    !         b)   Eim = 10^(-3/St) based on Slinn and Slinn (1980)
    !
    ! References:
    !  LEWIS and SCHWARTZ (2004), "SEA SALT AEROSOL PRODUCTION, MECHANISMS,
    !    METHODS AND MODELS" AGU monograph 152.
    !  SLINN and SLINN (1980), "PREDICTIONS FOR PARTICLE DEPOSITION ON
    !    NATURAL-WATERS" Atmos Environ (1980) vol. 14 (9) pp. 1013-1016.
    !  SLINN (1982), "PREDICTIONS FOR PARTICLE DEPOSITION TO VEGETATIVE
    !    CANOPIES" Atmos Environ (1982) vol. 16 (7) pp. 1785-1794.
    !==================================================================

    ! Number of bins for sea salt size distribution
    NR = INT((( Input_Opt%SALC_REDGE_um(2) - Input_Opt%SALA_REDGE_um(1) ) &
         / DR ) + 0.5e+0_f8 )

    ! Particle radius [cm]
    ! Bug fix: The Gerber [1985] growth should use the dry radius
    ! in micromenters and not cm. Replace RCM with RUM (jaegle 5/11/11)
    !RCM  = A_RADI(K) * 1.d2
    RUM  = A_RADI(K) * 1.e+6_f8

    ! Exponential factors used for hygroscopic growth
    ! Replace RCM with RUM (jaegle 5/11/11)
    !FAC1 = C1 * ( RCM**C2 )
    !FAC2 = C3 * ( RCM**C4 )
    FAC1 = C1 * ( RUM**C2 )
    FAC2 = C3 * ( RUM**C4 )

    ! Aerosol growth with relative humidity in radius [m]
    ! (Gerber, 1985) (bec, 12/8/04)
    ! Added safety check for LOG (phs, 6/11/08)
    RHBL    = MAX( TINY(RHB), RHB )

    ! Check whether we are over the oceans or not:
    ! Over oceans the RH in the viscous sublayer is set to 98%, following
    ! Lewis and Schwartz (2004), see discussion above (jaegle 5/11/11)
    IF (LUC == 14) THEN
       RHBL = 0.98
    ENDIF
    ! Corrected bug in Gerber formulation: use of LOG10  (jaegle 5/11/11)
    !RWET    = 0.01e+0_f8*(FAC1/(FAC2-DLOG(RHBL))+RCM**3.e+0_f8)**0.33e+0_f8
    !RWET = 1.d-6*(FAC1/(FAC2-LOG10(RHBL))+RUM**3.e+0_f8)**0.33333e+0_f8

    ! Use equation 5 in Lewis and Schwartz (2006) for sea salt growth [m]
    ! (jaegle 5/11/11)
    RWET = A_RADI(K) * (4.e+0_f8 / 3.7e+0_f8) * &
          ( (2.e+0_f8 - RHBL)/(1.e+0_f8 - RHBL) )**(1.e+0_f8/3.e+0_f8)

    ! Ratio dry over wet radii at the cubic power
    !RATIO_R = ( A_RADI(K) / RWET )**3.e+0_f8

    ! Diameter of the wet aerosol [m]
    DIAM  = RWET * 2.e+0_f8

    ! Density of the wet aerosol [kg/m3] (bec, 12/8/04)
    !DEN   = RATIO_R * A_DEN(K) + ( 1.e+0_f8 - RATIO_R ) * 1000.e+0_f8

    ! Above density calculation is chemically unsound because it ignores
    ! chemical solvation.
    ! Iteratively solve Tang et al., 1997 equation 5 to calculate density of
    ! wet aerosol (kg/m3)
    ! (bec, 6/17/10, jaegle 5/11/11)
    ! Redefine RATIO_R
    RATIO_R = A_RADI(K) / RWET

    ! Assume an initial density of 1000 kg/m3
    DEN  = 1000.e+0_f8
    DEN1 = 0.e+0_f8 !initialize (bec, 6/21/10)
    DO WHILE ( ABS( DEN1-DEN ) .gt. EPSI )
       ! First calculate weight percent of aerosol (kg_RH=0.8/kg_wet)
       WTP    = 100.e+0_f8 * A_DEN(K)/DEN * RATIO_R**3.e+0_f8
       ! Then calculate density of wet aerosol using equation 5
       ! in Tang et al., 1997 [kg/m3]
       DEN1   = ( 0.9971e+0_f8 + (A1 * WTP) + (A2 * WTP**2) + &
                (A3 * WTP**3) + (A4 * WTP**4) ) * 1000.e+0_f8

       ! Now calculate new weight percent using above density calculation
       WTP    = 100.e+0_f8 * A_DEN(K)/DEN1 * RATIO_R**3
       ! Now recalculate new wet density [kg/m3]
       DEN   = ( 0.9971e+0_f8 + (A1 * WTP) + (A2 * WTP**2) + &
               (A3 * WTP**3) + (A4 * WTP**4) ) * 1000.e+0_f8
    ENDDO

    ! Dp [um] = particle diameter
    DP    = DIAM * 1.e+6_f8

    ! Constant for settling velocity calculation
    CONST = DEN * DIAM**2 * g0 / 18.e+0_f8

    !=================================================================
    !   # air molecule number density
    !   num = P * 1d3 * 6.023d23 / (8.314 * Temp)
    !   # gas mean free path
    !   lamda = 1.d6/( 1.41421 * num * 3.141592 * (3.7d-10)**2 )
    !   # Slip correction
    !   Slip = 1. + 2. * lamda * (1.257 + 0.4 * exp( -1.1 * Dp &
    !          / (2. * lamda))) / Dp
    !=================================================================
    ! Note, Slip correction factor calculations following Seinfeld,
    ! pp464 which is thought to be more accurate but more computation
    ! required.
    !=================================================================

    ! Slip correction factor as function of (P*dp)
    PDP  = PRESS * DP
    SLIP = 1e+0_f8 + ( 15.60e+0_f8 + 7.0e+0_f8 * &
           EXP( -0.059e+0_f8 * PDP) ) / PDP

    !=================================================================
    ! Note, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
    ! which produce slip correction factore with small error
    ! compared to the above with less computation.
    !=================================================================

    ! Viscosity [Pa s] of air as a function of temp (K)
    VISC = 1.458e-6_f8 * (TEMP)**(1.5e+0_f8) / (TEMP + 110.4e+0_f8)

    ! Kinematic viscosity (Dynamic viscosity/Density)
    AIRVS= VISC / 1.2928e+0_f8

    ! Settling velocity [m/s]
    VTS  = CONST * SLIP / VISC

    ! This settling velocity is for the mid-point of the size bin.
    ! Need to integrate over the size bin, taking into account the
    ! mass distribution of sea-salt and the dependence of VTS on aerosol
    ! size. See WET_SETTLING in SEASALT_MOD.f for more details.
    ! (jaegle 5/11/11)
    SALT_MASS_TOTAL = 0e+0_f8
    VTS_WEIGHT      = 0e+0_f8
    ! Check what the min/max range of the SS size bins are
    IF ( RUM .le. Input_Opt%SALA_REDGE_um(2) ) THEN
      D0 = Input_Opt%SALA_REDGE_um(1)*2e+0_f8
      D1 = Input_Opt%SALA_REDGE_um(2)*2e+0_f8
    ELSE
      D0 = Input_Opt%SALC_REDGE_um(1)*2e+0_f8
      D1 = Input_Opt%SALC_REDGE_um(2)*2e+0_f8
    ENDIF

    DO ID = 1, NR
       ! Calculate mass of wet aerosol (Dw = wet diameter, D = dry diamter):
       ! Overall = dM/dDw = dV/dlnD * Rwet/Rdry * DEN /Rw
       IF (DMID(ID) .ge. D0 .and. DMID(ID) .le. D1 ) THEN
          DMIDW = DMID(ID) * RWET/A_RADI(K)   ! wet radius [um]
          SALT_MASS   = SALT_V(ID) * RWET/A_RADI(K) * DEN / &
                        (DMIDW*0.5e+0_f8)
          VTS_WEIGHT  = VTS_WEIGHT + &
                        SALT_MASS * VTS * (DMIDW/(RWET*1d6*2e+0_f8) )** &
                        2e+0_f8 * (2e+0_f8 * DR *  RWET/A_RADI(K))
          SALT_MASS_TOTAL = SALT_MASS_TOTAL+SALT_MASS * &
                            (2e+0_f8 * DR *  RWET/A_RADI(K))
       ENDIF
    ENDDO

    ! Final mass weighted setting velocity:
    VTS = VTS_WEIGHT/SALT_MASS_TOTAL

    ! Brownian diffusion constant for particle (m2/s)
    DIFF = BOLTZ * TEMP * SLIP / (3.e+0_f8 * PI * VISC * DIAM)

    ! Schmidt number
    SC   = AIRVS / DIFF
    EB   = 1.e+0_f8/SC**(gamma(LUC))

    ! Stokes number
    IF ( AA < 0e+0_f8 ) then
       ST   = VTS * USTAR * USTAR / ( AIRVS * g0 ) ! for smooth surface
       EIN  = 0e+0_f8
    ELSE
       ST   = VTS   * USTAR / ( g0 * AA )          ! for vegetated surfaces
       EIN  = 0.5e+0_f8 * ( DIAM / AA )**2
    ENDIF

    ! Use the formulation of Slinn and Slinn (1980) for the impaction over
    ! water surfaces (jaegle 5/11/11)
    IF (LUC == 14) THEN
       EIM  = 10.e+0_f8**( -3.e+0_f8/ ST )         ! for water surfaces
    ELSE
       EIM  = ( ST / ( ALPHA(LUC) + ST ) )**(BETA)
       EIM  = MIN( EIM, 0.6e+0_f8 )
    ENDIF

    IF (LUC == 11 .OR. LUC == 13 .OR. LUC == 14) THEN
       R1 = 1.e+0_f8
    ELSE
       R1 = EXP( -1e+0_f8 * SQRT( ST ) )
    ENDIF

    ! surface resistance for particle
    ! Use the formulation of Slinn and Slinn (1980) for the impaction over
    ! water surfaces (jaegle 5/11/11)
    IF (LUC == 14) THEN
#ifdef MODEL_GCHPCTM
       ! Include check that winds are non-zero to avoid div by 0 error
       IF ( IS_SAFE_DIV(1.e+0_f8, W10) ) THEN
          RS = 1.e+0_f8 / (USTAR**2.e+0_f8/ (W10*VON_KARMAN) * &
               (EB + EIM ) + VTS)
       ELSE
          RS = 1.e+0_f8 / (USTAR**2.e+0_f8/ (1.e-20_f8*VON_KARMAN) * &
               (EB + EIM ) + VTS)
          !PRINT *, "WARNING: 10 m winds are zero. Using small value", &
          !         "in drydep calculation of particle surface ", &
          !         "resistance to avoid divide by zero error."
       ENDIF
#else
       RS   = 1.e+0_f8 / (USTAR**2.e+0_f8/ (W10*VON_KARMAN) * &
              (EB + EIM ) + VTS)
#endif
    ELSE
       RS   = 1.e+0_f8 / (E0 * USTAR * (EB + EIM + EIN) * R1 )
    ENDIF

  END FUNCTION AERO_SFCRSII
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_weightss
!
! !DESCRIPTION: Subroutine INIT\_WEIGHTSS calculates the volume size
!  distribution of sea-salt. This only has to be done once. We assume that
!  sea-salt is the combination of a coarse mode and accumulation model
!  log-normal distribution functions. The resulting arrays are: DMID = diameter
!  of bin and SALT\_V = dV/dln(D) [in um3]. (jaegle 5/11/11)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_WEIGHTSS( Input_Opt )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt
!
! !REVISION HISTORY:
!  11 May 2011 - L. Jaegle   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: N
    REAL(f8)            :: SALT_MASS, SALT_MASS_TOTAL, VTS_WEIGHT
    REAL(f8)            :: DEDGE
    INTEGER             :: ID,NR
!
! !DEFINED PARAMETERS:
!
    ! increment of radius for integration of settling velocity (um)
    REAL(f8), PARAMETER :: DR    = 5.e-2_f8

    ! parameters for assumed size distribution of acc and coarse mode
    ! sea salt aerosols
    ! geometric dry mean diameters (microns)
    REAL(f8), PARAMETER :: RG_A  = 0.085e+0_f8
    REAL(f8), PARAMETER :: RG_C  = 0.4e+0_f8
    ! sigma of the size distribution
    REAL(f8), PARAMETER :: SIG_A = 1.5e+0_f8
    REAL(f8), PARAMETER :: SIG_C = 1.8e+0_f8

    ! Number of bins between the lowest bound of of the accumulation mode
    ! sea salt and the upper bound of the coarse mode sea salt.
    NR = INT((( Input_Opt%SALC_REDGE_um(2) - Input_Opt%SALA_REDGE_um(1) ) &
         / DR ) + 0.5e+0_f8 )

    !=================================================================
    ! Define the volume size distribution of sea-salt. This only has
    ! to be done once. We assume that sea-salt is the combination of a
    ! coarse mode and accumulation model log-normal distribution functions
    !=================================================================

    ! Lower edge of 0th bin diameter [um]
    DEDGE=Input_Opt%SALA_REDGE_um(1) * 2e+0_f8

    ! Loop over diameters
    DO ID = 1, NR

       ! Diameter of mid-point in microns
       DMID(ID)  = DEDGE + ( DR )

       ! Calculate the dry volume size distribution as the sum of two
       ! log-normal size distributions. The parameters for the size
       ! distribution are based on Reid et al. and Quinn et al.
       ! The scaling factors 13. and 0.8 for acc and coarse mode aerosols
       ! are chosen to obtain a realistic distribution
       ! SALT_V (D) = dV/dln(D) [um3]
       SALT_V(ID) = PI / 6e+0_f8* (DMID(ID)**3) * (          &
                    13e+0_f8*exp(-0.5*( LOG(DMID(ID))-       &
                    LOG(RG_A*2e+0_f8) )**2e+0_f8/            &
                              LOG(SIG_A)**2e+0_f8 )          &
                    /( sqrt(2e+0_f8 * PI) * LOG(SIG_A) )  +  &
                    0.8e+0_f8*exp(-0.5*( LOG(DMID(ID))-      &
                    LOG(RG_C*2e+0_f8) )**2e+0_f8/            &
                              LOG(SIG_C)**2e+0_f8)           &
                    /( sqrt(2e+0_f8 * PI) * LOG(SIG_C) )  )

       ! update the next edge
       DEDGE = DEDGE + DR*2e+0_f8
    ENDDO

  END SUBROUTINE INIT_WEIGHTSS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dust_sfcrsi
!
! !DESCRIPTION: Function DUST\_SFCRSI computes the aerodynamic resistance of
!  dust aerosol species according to Seinfeld et al 96.  We do not consider
!  hygroscopic growth of the dust aerosol particles. (rjp, tdf, bmy, bec,
!  4/1/04, 4/15/05)
!\\
!\\
! !INTERFACE:
!
  FUNCTION DUST_SFCRSI( K, II, PRESS, TEMP, USTAR ) RESULT( RS )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: K       ! Drydep species (range: 1-NUMDEP)
    INTEGER,  INTENT(IN) :: II      ! Surface type index of GEOS-CHEM
    REAL(f8), INTENT(IN) :: PRESS   ! Pressure [kPa]
    REAL(f8), INTENT(IN) :: TEMP    ! Temperature [K]
    REAL(f8), INTENT(IN) :: USTAR   ! Friction velocity [m/s]
!
! !RETURN VALUE:
!
    REAL(f8)             :: RS      ! Surface resistance for particles [s/m]
!
! !REVISION HISTORY:
!  01 Apr 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: N
    REAL(f8), PARAMETER   :: C1 = 0.7674e+0_f8
    REAL(f8), PARAMETER   :: C2 = 3.079e+0_f8
    REAL(f8), PARAMETER   :: C3 = 2.573e-11_f8
    REAL(f8), PARAMETER   :: C4 = -1.424e+0_f8
    REAL(f8), PARAMETER   :: BETA  = 2.e+0_f8
    REAL(f8), PARAMETER   :: E0 = 1.e+0_f8
    REAL(f8)  :: AIRVS       ! kinematic viscosity of Air (m^2/s)
    REAL(f8)  :: DP          ! Diameter of aerosol [um]
    REAL(f8)  :: PDP         ! Press * Dp
    REAL(f8)  :: CONST       ! Constant for settling velocity calculations
    REAL(f8)  :: SLIP        ! Slip correction factor
    REAL(f8)  :: VISC        ! Viscosity of air (Pa s)
    REAL(f8)  :: DIFF        ! Brownian Diffusion constant for particles (m2/s)
    REAL(f8)  :: SC, ST      ! Schmidt and Stokes number (nondim)
    REAL(f8)  :: DIAM, DEN
    REAL(f8)  :: EB, EIM, EIN, R1, AA, VTS

    !=================================================================
    ! Ref. Zhang et al., AE 35(2001) 549-560 and Seinfeld(1986)
    !
    ! Model theory
    !    Vd = Vs + 1./(Ra+Rs)
    !      where Vs is the gravitational settling velocity,
    !      Ra is the aerodynamic resistance above the canopy
    !      Rs  is the surface resistance
    !    Here we calculate Rs only..
    !    Rs = 1 / (Eo*Ustar*(Eb+Eim+Ein)*R1)
    !      where Eo is an empirical constant ( = 3.)
    !      Ustar is the friction velocity
    !      Collection efficiency from
    !        Eb,  [Brownian diffusion]
    !        Eim, [Impaction]
    !        Ein, [Interception]
    !      R1 is the correction factor representing the fraction
    !         of particles that stick to the surface.
    !=================================================================
    !      Eb is a funciont of Schmidt number, Eb = Sc^(-gamma)
    !         Sc = v/D, v (the kinematic viscosity of air)
    !                   D (particle brownian diffusivity)
    !         r usually lies between 1/2 and 2/3
    !      Eim is a function of Stokes number, St
    !          St = Vs * Ustar / (g0 * A)   for vegetated surfaces
    !          St = Vs * Ustar * Ustar / v  for smooth surface
    !          A is the characteristic radius of collectors.
    !
    !       1) Slinn (1982)
    !           Eim = 10^(-3/St)          for smooth surface
    !           Eim = St^2 / ( 1 + St^2 ) for vegetative canopies
    !       2) Peters and Eiden (1992)
    !           Eim = ( St / ( alpha + St ) )^(beta)
    !                alpha(=0.8) and beta(=2) are constants
    !       3) Giorgi (1986)
    !           Eim = St^2 / ( 400 + St^2 )     for smooth surface
    !           Eim = ( St / (0.6 + St) )^(3.2) for vegetative surface
    !       4) Davidson et al.(1982)
    !           Eim = St^3 / (St^3+0.753*St^2+2.796St-0.202) for grassland
    !       5) Zhang et al.(2001) used 2) method with alpha varying with
    !          vegetation type and beta equal to 2
    !
    !      Ein = 0.5 * ( Dp / A )^2
    !
    !      R1 (Particle rebound)  = exp(-St^0.5)
    !=================================================================

    ! Particle diameter [m]
    DIAM  = A_RADI(K) * 2.e+0_f8

    ! Particle density [kg/m3]
    DEN   = A_DEN(K)

    ! Dp [um] = particle diameter
    DP    = DIAM * 1.e+6_f8

    ! Constant for settling velocity calculation
    CONST = DEN * DIAM**2 * g0 / 18.e+0_f8

    !=================================================================
    !   # air molecule number density
    !   num = P * 1d3 * 6.023d23 / (8.314 * Temp)
    !   # gas mean free path
    !   lamda = 1.d6/( 1.41421 * num * 3.141592 * (3.7d-10)**2 )
    !   # Slip correction
    !   Slip = 1. + 2. * lamda * (1.257 + 0.4 * exp( -1.1 * Dp &
    !          / (2. * lamda))) / Dp
    !================================================================
    ! Note, Slip correction factor calculations following Seinfeld,
    ! pp464 which is thought to be more accurate but more computation
    ! required.
    !=================================================================

    ! Slip correction factor as function of (P*dp)
    PDP  = PRESS * DP
    SLIP = 1e+0_f8 + ( 15.60e+0_f8 + 7.0e+0_f8 * &
           EXP( -0.059e+0_f8 * PDP ) ) / PDP

    !=================================================================
    ! Note, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
    ! which produce slip correction factore with small error
    ! compared to the above with less computation.
    !=================================================================

    ! Viscosity [Pa s] of air as a function of temp (K)
    VISC = 1.458e-6_f8 * (TEMP)**(1.5e+0_f8) / (TEMP + 110.4e+0_f8)

    ! Kinematic viscosity (Dynamic viscosity/Density)
    AIRVS= VISC / 1.2928e+0_f8

    ! Settling velocity [m/s]
    VTS  = CONST * SLIP / VISC

    ! Brownian diffusion constant for particle (m2/s)
    DIFF = BOLTZ * TEMP * SLIP / (3.e+0_f8 * Pi * VISC * DIAM)

    ! Schmidt number and Diffusion term
    SC   = AIRVS / DIFF
    EB   = SC**(-0.666667e+0_f8)

    ! Stokes number and impaction term
    ST   = VTS * USTAR * USTAR / ( AIRVS * g0 )
    EIM  = 10.e+0_f8**(-3.e+0_f8 / ST)

    ! surface resistance for particle
    RS   = 1.e+0_f8 / ( E0 * USTAR * (EB + EIM) )

  END FUNCTION DUST_SFCRSI
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: adust_sfcrsii
!
! !DESCRIPTION: Function ADUST\_SFCRSII computes the aerodynamic resistance of
!  non-size resolved aerosol according to Zhang et al 2001.  We do not consider
!  the hygroscopic growth of the aerosol particles. (rjp, tdf, bec, bmy,
!  4/1/04, 4/15/05)
!\\
!\\
!  This routine is used for all aerosols except dust, sulfate, and seasalt
!  (hotp 7/31/09)
!\\
!\\
! !INTERFACE:
!
  FUNCTION ADUST_SFCRSII( K, II, PRESS, TEMP, USTAR ) RESULT( RS )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: K     ! Drydep species index (range: 1-NUMDEP)
    INTEGER,  INTENT(IN) :: II    ! Surface type index of GEOS-CHEM
    REAL(f8), INTENT(IN) :: PRESS ! Pressure [kPa] (1 mb = 100 Pa = 0.1 kPa)
    REAL(f8), INTENT(IN) :: TEMP  ! Temperature [K]
    REAL(f8), INTENT(IN) :: USTAR ! Friction velocity [m/s]
!
! !RETURN VALUE:
!
    REAL(f8)             :: RS    ! Surface resistance for particles [s/m]
!
! !REVISION HISTORY:
!  01 Apr 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: N
    REAL(f8), PARAMETER   :: C1 = 0.7674e+0_f8
    REAL(f8), PARAMETER   :: C2 = 3.079e+0_f8
    REAL(f8), PARAMETER   :: C3 = 2.573e-11_f8
    REAL(f8), PARAMETER   :: C4 = -1.424e+0_f8
    REAL(f8), PARAMETER   :: BETA  = 2.e+0_f8
    REAL(f8), PARAMETER   :: E0 = 3.e+0_f8
    REAL(f8)  :: AIRVS       ! kinematic viscosity of Air (m^2/s)
    REAL(f8)  :: DP          ! Diameter of aerosol [um]
    REAL(f8)  :: PDP         ! Press * Dp
    REAL(f8)  :: CONST       ! Constant for settling velocity calculations
    REAL(f8)  :: SLIP        ! Slip correction factor
    REAL(f8)  :: VISC        ! Viscosity of air (Pa s)
    REAL(f8)  :: DIFF        ! Brownian Diffusion constant for particles (m2/s)
    REAL(f8)  :: SC, ST      ! Schmidt and Stokes number (nondim)
    REAL(f8)  :: DIAM, DEN
    REAL(f8)  :: EB, EIM, EIN, R1, AA, VTS

    !=======================================================================
    !   #  LUC [Zhang et al., 2001]                GEOS-CHEM LUC (Corr. #)
    !-----------------------------------------------------------------------
    !   1 - Evergreen needleleaf trees             Snow/Ice          (12)
    !   2 - Evergreen broadleaf trees              Deciduous forest  ( 4)
    !   3 - Deciduous needleleaf trees             Coniferous forest ( 1)
    !   4 - Deciduous broadleaf trees              Agricultural land ( 7)
    !   5 - Mixed broadleaf and needleleaf trees   Shrub/grassland   (10)
    !   6 - Grass                                  Amazon forest     ( 2)
    !   7 - Crops and mixed farming                Tundra            ( 9)
    !   8 - Desert                                 Desert            ( 8)
    !   9 - Tundra                                 Wetland           (11)
    !  10 - Shrubs and interrupted woodlands       Urban             (15)
    !  11 - Wet land with plants                   Water             (14)
    !  12 - Ice cap and glacier
    !  13 - Inland water
    !  14 - Ocean
    !  15 - Urban
    !=======================================================================
    ! GEOS-CHEM LUC              1, 2, 3, 4, 5, 6, 7  8, 9,10,11
    INTEGER :: LUCINDEX(11) = (/12, 4, 1, 7,10, 2, 9, 8,11,15,14/)
    INTEGER :: LUC

    !=======================================================================
    !   LUC       1,    2,    3,    4,    5,    6,    7,    8,
    !   alpha   1.0,  0.6,  1.1,  0.8,  0.8,  1.2,  1.2, 50.0,
    !   gamma  0.56, 0.58, 0.56, 0.56, 0.56, 0.54, 0.54, 0.54
    !
    !   LUC       9,   10,   11,   12,   13,   14,   15
    !   alpha  50.0,  1,3,  2.0, 50.0,100.0,100.0,  1.5
    !   gamma  0.54, 0.54, 0.54, 0.54, 0.50, 0.50, 0.56
    !=======================================================================
    REAL(f8)  :: ALPHA(15) = (/   1.0e+0_f8,   0.6e+0_f8,  1.1e+0_f8, & 
                                  0.8e+0_f8,   0.8e+0_f8,  1.2e+0_f8, &
                                  1.2e+0_f8,  50.0e+0_f8, 50.0e+0_f8, &
                                  1.3e+0_f8,   2.0e+0_f8, 50.0e+0_f8, &
                                100.0e+0_f8, 100.0e+0_f8,  1.5e+0_f8  /)

    REAL(f8)  :: GAMMA(15) = (/ 0.56e+0_f8, 0.58e+0_f8, 0.56e+0_f8, &
                                0.56e+0_f8, 0.56e+0_f8, 0.54e+0_f8, &
                                0.54e+0_f8, 0.54e+0_f8, 0.54e+0_f8, &
                                0.54e+0_f8, 0.54e+0_f8, 0.54e+0_f8, &
                                0.50e+0_f8, 0.50e+0_f8, 0.56e+0_f8  /)

    !...A unit is (mm) so multiply by 1.D-3 to (m)
    !   LUC       1,    2,    3,    4,    5,    6,    7,    8,
    !   SC1     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,
    !   SC2     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,
    ! A SC3     2.0,  5.0,  5.0, 10.0,  5.0,  5.0,  5.0,-999.,
    !   SC4     2.0,  5.0,  5.0, 10.0,  5.0,  5.0,  5.0,-999.,
    !   SC5     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,
    !
    !   LUC       9,   10,   11,   12,   13,   14,   15
    !   SC1   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
    !   SC2   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
    ! A SC3   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
    !   SC4   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
    !   SC5   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
    REAL(f8)  :: A(15,5)

    REAL(f8)  :: Aavg(15)

    DATA   A / 2.0e+0_f8,   5.0e+0_f8,   2.0e+0_f8,   5.0e+0_f8,  5.0e+0_f8, &
               2.0e+0_f8,   2.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
              10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
              
               2.0e+0_f8,   5.0e+0_f8,   2.0e+0_f8,   5.0e+0_f8,  5.0e+0_f8, &
               2.0e+0_f8,   2.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
              10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
              
               2.0e+0_f8,   5.0e+0_f8,   5.0e+0_f8,  10.0e+0_f8,  5.0e+0_f8, &
               5.0e+0_f8,   5.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
              10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
              
               2.0e+0_f8,   5.0e+0_f8,   5.0e+0_f8,  10.0e+0_f8,  5.0e+0_f8, &
               5.0e+0_f8,   5.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
              10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
              
               2.0e+0_f8,   5.0e+0_f8,   2.0e+0_f8,   5.0e+0_f8,  5.0e+0_f8, &
               2.0e+0_f8,   2.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
              10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8  /

    ! Annual average of A
    Aavg(:) = (A(:,1)+A(:,2)+A(:,3)+A(:,4)+A(:,5))/5.
    LUC     = LUCINDEX(II)
    AA      = Aavg(LUC) * 1.e-3_f8

    !=================================================================
    !...Ref. Zhang et al., AE 35(2001) 549-560
    !.
    !...Model theroy
    !    Vd = Vs + 1./(Ra+Rs)
    !      where Vs is the gravitational settling velocity,
    !      Ra is the aerodynamic resistance above the canopy
    !      Rs  is the surface resistance
    !    Here we calculate Rs only..
    !    Rs = 1 / (Eo*Ustar*(Eb+Eim+Ein)*R1)
    !      where Eo is an empirical constant ( = 3.)
    !      Ustar is the friction velocity
    !      Collection efficiency from
    !        Eb,  [Brownian diffusion]
    !        Eim, [Impaction]
    !        Ein, [Interception]
    !      R1 is the correction factor representing the fraction
    !         of particles that stick to the surface.
    !=======================================================================
    !      Eb is a funciont of Schmidt number, Eb = Sc^(-gamma)
    !         Sc = v/D, v (the kinematic viscosity of air)
    !                   D (particle brownian diffusivity)
    !         r usually lies between 1/2 and 2/3
    !      Eim is a function of Stokes number, St
    !          St = Vs * Ustar / (g0 * A)   for vegetated surfaces
    !          St = Vs * Ustar * Ustar / v  for smooth surface
    !          A is the characteristic radius of collectors.
    !
    !       1) Slinn (1982)
    !           Eim = 10^(-3/St)          for smooth surface
    !           Eim = St^2 / ( 1 + St^2 ) for vegetative canopies
    !       2) Peters and Eiden (1992)
    !           Eim = ( St / ( alpha + St ) )^(beta)
    !                alpha(=0.8) and beta(=2) are constants
    !       3) Giorgi (1986)
    !           Eim = St^2 / ( 400 + St^2 )     for smooth surface
    !           Eim = ( St / (0.6 + St) )^(3.2) for vegetative surface
    !       4) Davidson et al.(1982)
    !           Eim = St^3 / (St^3+0.753*St^2+2.796St-0.202) for grassland
    !       5) Zhang et al.(2001) used 2) method with alpha varying with
    !          vegetation type and beta equal to 2
    !
    !      Ein = 0.5 * ( Dp / A )^2
    !
    !      R1 (Particle rebound)  = exp(-St^0.5)
    !=================================================================

    ! Particle diameter [m] hotp 10/26/07
    DIAM  = 0.5e-6_f8

    ! Particle density [kg/m3] hotp 10/26/07
    DEN   = 1500

    ! Dp [um] = particle diameter
    DP    = DIAM * 1.e+6_f8

    ! Constant for settling velocity calculation
    CONST = DEN * DIAM**2 * g0 / 18.e+0_f8

    !=================================================================
    !   # air molecule number density
    !   num = P * 1d3 * 6.023d23 / (8.314 * Temp)
    !   # gas mean free path
    !   lamda = 1.d6/( 1.41421 * num * 3.141592 * (3.7d-10)**2 )
    !   # Slip correction
    !   Slip = 1. + 2. * lamda * (1.257 + 0.4 * exp( -1.1 * Dp &
    !          / (2. * lamda))) / Dp
    !=================================================================
    ! Note, Slip correction factor calculations following Seinfeld,
    ! pp464 which is thought to be more accurate but more computation
    ! required.
    !=================================================================

    ! Slip correction factor as function of (P*dp)
    PDP  = PRESS * DP
    SLIP = 1e+0_f8 + ( 15.60e+0_f8 + 7.0e+0_f8 * &
           EXP( -0.059e+0_f8 * PDP) ) / PDP

    !=================================================================
    ! Note, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
    ! which produce slip correction factore with small error
    ! compared to the above with less computation.
    !=================================================================

    ! Viscosity [Pa s] of air as a function of temp (K)
    VISC = 1.458e-6_f8 * (TEMP)**(1.5e+0_f8) / (TEMP + 110.4e+0_f8)

    ! Kinematic viscosity (Dynamic viscosity/Density)
    AIRVS= VISC / 1.2928e+0_f8

    ! Settling velocity [m/s]
    VTS  = CONST * SLIP / VISC

    ! Brownian diffusion constant for particle (m2/s)
    DIFF = BOLTZ * TEMP * SLIP / (3.e+0_f8 * Pi * VISC * DIAM)

    ! Schmidt number
    SC   = AIRVS / DIFF
    EB   = 1.e+0_f8/SC**(gamma(LUC))

    ! Stokes number
    IF ( AA < 0e+0_f8 ) then
       ST   = VTS * USTAR * USTAR / ( AIRVS * g0 ) ! for smooth surface
       EIN  = 0e+0_f8
    ELSE
       ST   = VTS   * USTAR / ( g0 * AA )          ! for vegetated surfaces
       EIN  = 0.5e+0_f8 * ( DIAM / AA )**2
    ENDIF

    EIM  = ( ST / ( ALPHA(LUC) + ST ) )**(BETA)

    EIM  = MIN( EIM, 0.6e+0_f8 )

    IF (LUC == 11 .OR. LUC == 13 .OR. LUC == 14) THEN
       R1 = 1.e+0_f8
    ELSE
       R1 = EXP( -1e+0_f8 * SQRT( ST ) )
    ENDIF

    ! surface resistance for particle
    RS   = 1.e0_f8 / (E0 * USTAR * (EB + EIM + EIN) * R1 )

  END FUNCTION ADUST_SFCRSII
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dust_sfcrsii
!
! !DESCRIPTION: Function DUST\_SFCRSII computes the aerodynamic resistance of
!  dust aerosol species according to Zhang et al 2001.  We do not consider the
!  hygroscopic growth of the aerosol particles. (rjp, tdf, bec, bmy, 4/1/04,
!  4/15/05)
!\\
!\\
! !INTERFACE:
!

  FUNCTION DUST_SFCRSII( K, II, PRESS, TEMP, USTAR, DIAM, DEN ) &
       RESULT( RS )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: K       ! Drydep species index (range: 1-NUMDEP)
    INTEGER,  INTENT(IN) :: II      ! Surface type index of GEOS-CHEM
    REAL(f8), INTENT(IN) :: PRESS   ! Pressure [kPa]
    REAL(f8), INTENT(IN) :: TEMP    ! Temperature [K]
    REAL(f8), INTENT(IN) :: USTAR   ! Friction velocity [m/s]
    REAL(f8), INTENT(IN) :: DIAM    ! Particle diameter [m]
    REAL(f8), INTENT(IN) :: DEN     ! Particle density [kg/m3]
!
! !RETURN VALUE:
!
    REAL(f8)             :: RS      ! Surface resistance for particles [s/m]
!
! !REVISION HISTORY:
!  01 Apr 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: N
    REAL(f8), PARAMETER   :: C1 = 0.7674e+0_f8
    REAL(f8), PARAMETER   :: C2 = 3.079e+0_f8
    REAL(f8), PARAMETER   :: C3 = 2.573e-11_f8
    REAL(f8), PARAMETER   :: C4 = -1.424e+0_f8
    REAL(f8), PARAMETER   :: BETA  = 2.e+0_f8
    REAL(f8), PARAMETER   :: E0 = 3.e+0_f8
    REAL(f8)  :: AIRVS       ! kinematic viscosity of Air (m^2/s)
    REAL(f8)  :: DP          ! Diameter of aerosol [um]
    REAL(f8)  :: PDP         ! Press * Dp
    REAL(f8)  :: CONST       ! Constant for settling velocity calculations
    REAL(f8)  :: SLIP        ! Slip correction factor
    REAL(f8)  :: VISC        ! Viscosity of air (Pa s)
    REAL(f8)  :: DIFF        ! Brownian Diffusion constant for particles (m2/s)
    REAL(f8)  :: SC, ST      ! Schmidt and Stokes number (nondim)
    REAL(f8)  :: EB, EIM, EIN, R1, AA, VTS

    !=======================================================================
    !   #  LUC [Zhang et al., 2001]                GEOS-CHEM LUC (Corr. #)
    !-----------------------------------------------------------------------
    !   1 - Evergreen needleleaf trees             Snow/Ice          (12)
    !   2 - Evergreen broadleaf trees              Deciduous forest  ( 4)
    !   3 - Deciduous needleleaf trees             Coniferous forest ( 1)
    !   4 - Deciduous broadleaf trees              Agricultural land ( 7)
    !   5 - Mixed broadleaf and needleleaf trees   Shrub/grassland   (10)
    !   6 - Grass                                  Amazon forest     ( 2)
    !   7 - Crops and mixed farming                Tundra            ( 9)
    !   8 - Desert                                 Desert            ( 8)
    !   9 - Tundra                                 Wetland           (11)
    !  10 - Shrubs and interrupted woodlands       Urban             (15)
    !  11 - Wet land with plants                   Water             (14)
    !  12 - Ice cap and glacier
    !  13 - Inland water
    !  14 - Ocean
    !  15 - Urban
    !=======================================================================
    ! GEOS-CHEM LUC              1, 2, 3, 4, 5, 6, 7  8, 9,10,11
    INTEGER :: LUCINDEX(11) = (/12, 4, 1, 7,10, 2, 9, 8,11,15,14/)
    INTEGER :: LUC

    !=======================================================================
    !   LUC       1,    2,    3,    4,    5,    6,    7,    8,
    !   alpha   1.0,  0.6,  1.1,  0.8,  0.8,  1.2,  1.2, 50.0,
    !   gamma  0.56, 0.58, 0.56, 0.56, 0.56, 0.54, 0.54, 0.54
    !
    !   LUC       9,   10,   11,   12,   13,   14,   15
    !   alpha  50.0,  1,3,  2.0, 50.0,100.0,100.0,  1.5
    !   gamma  0.54, 0.54, 0.54, 0.54, 0.50, 0.50, 0.56
    !=======================================================================
    REAL(f8)  :: ALPHA(15) = (/   1.0e+0_f8,   0.6e+0_f8,  1.1e+0_f8, &
                                  0.8e+0_f8,   0.8e+0_f8,  1.2e+0_f8, &
                                  1.2e+0_f8,  50.0e+0_f8, 50.0e+0_f8, &
                                  1.3e+0_f8,   2.0e+0_f8, 50.0e+0_f8, &
                                100.0e+0_f8, 100.0e+0_f8,  1.5e+0_f8  /)

    REAL(f8)  :: GAMMA(15) = (/ 0.56e+0_f8, 0.58e+0_f8, 0.56e+0_f8, &
                                0.56e+0_f8, 0.56e+0_f8, 0.54e+0_f8, &
                                0.54e+0_f8, 0.54e+0_f8, 0.54e+0_f8, &
                                0.54e+0_f8, 0.54e+0_f8, 0.54e+0_f8, &
                                0.50e+0_f8, 0.50e+0_f8, 0.56e+0_f8  /)

    !...A unit is (mm) so multiply by 1.D-3 to (m)
    !   LUC       1,    2,    3,    4,    5,    6,    7,    8,
    !   SC1     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,
    !   SC2     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,
    ! A SC3     2.0,  5.0,  5.0, 10.0,  5.0,  5.0,  5.0,-999.,
    !   SC4     2.0,  5.0,  5.0, 10.0,  5.0,  5.0,  5.0,-999.,
    !   SC5     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,
    !
    !   LUC       9,   10,   11,   12,   13,   14,   15
    !   SC1   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
    !   SC2   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
    ! A SC3   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
    !   SC4   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
    !   SC5   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
    REAL(f8)  :: A(15,5)

    REAL(f8)  :: Aavg(15)

    DATA   A /  2.0e+0_f8,   5.0e+0_f8,   2.0e+0_f8,   5.0e+0_f8,  5.0e+0_f8, &
                2.0e+0_f8,   2.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
               10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
               
                2.0e+0_f8,   5.0e+0_f8,   2.0e+0_f8,   5.0e+0_f8,  5.0e+0_f8, &
                2.0e+0_f8,   2.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
               10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
               
                2.0e+0_f8,   5.0e+0_f8,   5.0e+0_f8,  10.0e+0_f8,  5.0e+0_f8, &
                5.0e+0_f8,   5.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
               10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
               
                2.0e+0_f8,   5.0e+0_f8,   5.0e+0_f8,  10.0e+0_f8,  5.0e+0_f8, &
                5.0e+0_f8,   5.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
               10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
               
                2.0e+0_f8,   5.0e+0_f8,   2.0e+0_f8,   5.0e+0_f8,  5.0e+0_f8, &
                2.0e+0_f8,   2.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, &
               10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8  /

    ! Annual average of A
    Aavg(:) = (A(:,1)+A(:,2)+A(:,3)+A(:,4)+A(:,5))/5.
    LUC     = LUCINDEX(II)
    AA      = Aavg(LUC) * 1.e-3_f8

    !=================================================================
    !...Ref. Zhang et al., AE 35(2001) 549-560
    !.
    !...Model theroy
    !    Vd = Vs + 1./(Ra+Rs)
    !      where Vs is the gravitational settling velocity,
    !      Ra is the aerodynamic resistance above the canopy
    !      Rs  is the surface resistance
    !    Here we calculate Rs only..
    !    Rs = 1 / (Eo*Ustar*(Eb+Eim+Ein)*R1)
    !      where Eo is an empirical constant ( = 3.)
    !      Ustar is the friction velocity
    !      Collection efficiency from
    !        Eb,  [Brownian diffusion]
    !        Eim, [Impaction]
    !        Ein, [Interception]
    !      R1 is the correction factor representing the fraction
    !         of particles that stick to the surface.
    !=======================================================================
    !      Eb is a funciont of Schmidt number, Eb = Sc^(-gamma)
    !         Sc = v/D, v (the kinematic viscosity of air)
    !                   D (particle brownian diffusivity)
    !         r usually lies between 1/2 and 2/3
    !      Eim is a function of Stokes number, St
    !          St = Vs * Ustar / (g0 * A)   for vegetated surfaces
    !          St = Vs * Ustar * Ustar / v  for smooth surface
    !          A is the characteristic radius of collectors.
    !
    !       1) Slinn (1982)
    !           Eim = 10^(-3/St)          for smooth surface
    !           Eim = St^2 / ( 1 + St^2 ) for vegetative canopies
    !       2) Peters and Eiden (1992)
    !           Eim = ( St / ( alpha + St ) )^(beta)
    !                alpha(=0.8) and beta(=2) are constants
    !       3) Giorgi (1986)
    !           Eim = St^2 / ( 400 + St^2 )     for smooth surface
    !           Eim = ( St / (0.6 + St) )^(3.2) for vegetative surface
    !       4) Davidson et al.(1982)
    !           Eim = St^3 / (St^3+0.753*St^2+2.796St-0.202) for grassland
    !       5) Zhang et al.(2001) used 2) method with alpha varying with
    !          vegetation type and beta equal to 2
    !
    !      Ein = 0.5 * ( Dp / A )^2
    !
    !      R1 (Particle rebound)  = exp(-St^0.5)
    !=================================================================

    ! Dp [um] = particle diameter
    DP    = DIAM * 1.e+6_f8

    ! Constant for settling velocity calculation
    CONST = DEN * DIAM**2 * g0 / 18.e+0_f8

    !=================================================================
    !   # air molecule number density
    !   num = P * 1d3 * 6.023d23 / (8.314 * Temp)
    !   # gas mean free path
    !   lamda = 1.d6/( 1.41421 * num * 3.141592 * (3.7d-10)**2 )
    !   # Slip correction
    !   Slip = 1. + 2. * lamda * (1.257 + 0.4 * exp( -1.1 * Dp &
    !          / (2. * lamda))) / Dp
    !=================================================================
    ! Note, Slip correction factor calculations following Seinfeld,
    ! pp464 which is thought to be more accurate but more computation
    ! required.
    !=================================================================

    ! Slip correction factor as function of (P*dp)
    PDP  = PRESS * DP
    SLIP = 1e+0_f8 + ( 15.60e+0_f8 + 7.0e+0_f8 * &
           EXP( -0.059e+0_f8 * PDP) ) / PDP

    !=================================================================
    ! Note, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
    ! which produce slip correction factore with small error
    ! compared to the above with less computation.
    !=================================================================

    ! Viscosity [Pa s] of air as a function of temp (K)
    VISC = 1.458e-6_f8 * (TEMP)**(1.5e+0_f8) / (TEMP + 110.4e+0_f8)

    ! Kinematic viscosity (Dynamic viscosity/Density)
    AIRVS= VISC / 1.2928e+0_f8

    ! Settling velocity [m/s]
    VTS  = CONST * SLIP / VISC

    ! Brownian diffusion constant for particle (m2/s)
    DIFF = BOLTZ * TEMP * SLIP / (3.e+0_f8 * Pi * VISC * DIAM)

    ! Schmidt number
    SC   = AIRVS / DIFF
    EB   = 1.e+0_f8/SC**(gamma(LUC))

    ! Stokes number
    IF ( AA < 0e+0_f8 ) then
       ST   = VTS * USTAR * USTAR / ( AIRVS * g0 ) ! for smooth surface
       EIN  = 0e+0_f8
    ELSE
       ST   = VTS   * USTAR / ( g0 * AA )          ! for vegetated surfaces
       EIN  = 0.5e+0_f8 * ( DIAM / AA )**2
    ENDIF

    EIM  = ( ST / ( ALPHA(LUC) + ST ) )**(BETA)

    EIM  = MIN( EIM, 0.6e+0_f8 )

    IF (LUC == 11 .OR. LUC == 13 .OR. LUC == 14) THEN
       R1 = 1.D0
    ELSE
       R1 = EXP( -1e+0_f8 * SQRT( ST ) )
    ENDIF

    ! surface resistance for particle
    RS   = 1.e+0_f8 / (E0 * USTAR * (EB + EIM + EIN) * R1 )

  END FUNCTION DUST_SFCRSII
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_drydep
!
! !DESCRIPTION: Subroutine INIT\_DRYDEP initializes certain variables for the
!  GEOS-CHEM dry deposition subroutines. (bmy, 11/19/02, 10/19/09)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_DRYDEP( Input_Opt,  State_Chm, State_Diag, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE Species_Mod,    ONLY : Species
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  We now know how many drydep species there are before INIT_DRYDEP is
!  called.  This allows us to get rid of MAXDEP.  NUMDEP should be
!  equal to State_Chm%nDryDep, otherwise there is an error.
!
!  Also note: we need to use the actual molecular weights instead of
!  the emitted molecular weights.  These are necessary for the Schmidt #
!  computation.
!
! !REVISION HISTORY:
!  19 Nov 2002 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL                :: LDRYD
    LOGICAL                :: IS_Hg
    INTEGER                :: N

    ! Strings
    CHARACTER(LEN=255)     :: Msg, ErrMsg, ThisLoc

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    !=================================================================
    ! INIT_DRYDEP begins here!
    !=================================================================

    ! Initialize
    RC        = GC_SUCCESS
    ErrMsg    = ''
    ThisLoc   = ' -> at Init_Drydep (in module GeosCore/drydep_mod.F)'

#ifdef MODEL_WRF
    ! If the dry deposition module has already been initialized,
    ! the arrays do not need to be allocated again, as they are only
    ! dependent on the chemistry configuration (State_Chm%nDryDep)
    !
    ! This is necessary for integrating GEOS-Chem with a variable
    ! domain model like WRF-GC, where multiple instances of GEOS-Chem
    ! run in the same CPU. (hplin, 2/16/2019)
    IF ( ALLOCATED( A_DEN ) ) RETURN
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%% NOTE: Because READ_DRYDEP_INPUTS reads info from a netCDF %%%
    !%%% file, we may have to broadcast these.  However, the file  %%%
    !%%% dimensions are not very great (10 or 20 indices each)     %%%
    !%%% (bmy, 12/11/12)                                           %%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Read drydep inputs from the netCDF file
    ! Save Olson indices in INDOLSON array, in order to avoid
    ! confusion w/ previously-assinged variable name IOLSON
    !
    ! NOTE: For dry-run simulations, print filename and exit.
    CALL READ_DRYDEP_INPUTS( Input_Opt,                   &
                             DRYCOEFF,  INDOLSON, IDEP,   &
                             IWATER,    NWATER,   IZO,    &
                             IDRYDEP,   IRI,      IRLU,   &
                             IRAC,      IRGSS,    IRGSO,  &
                             IRCLS,     IRCLO,    IVSMAX, &
                             RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Read_Drydep_Inputs"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Exit if this is a dry-run simulation
    IF ( Input_Opt%DryRun ) RETURN

    !=================================================================
    ! For regular simulations, continue to initialize drydep
    !=================================================================

    IS_Hg     = Input_Opt%ITS_A_MERCURY_SIM
    LDRYD     = Input_Opt%LDRYD
    NUMDEP    = 0
    id_ACET   = 0
    id_O3     = 0
    id_ALD2   = 0
    id_MENO3  = 0
    id_ETNO3  = 0
    id_HNO3   = IND_('HNO3'  )
    id_PAN    = IND_('PAN'   )
    id_IHN1   = IND_('IHN1'  )

    !===================================================================
    ! Arrays that hold information about dry-depositing species
    ! Only allocate these if dry deposition is activated
    !===================================================================
    IF ( State_Chm%nDryDep > 0 ) THEN

       ! Aerosol density [kg/m3]
       ALLOCATE( A_DEN( State_Chm%nDryDep ), STAT=RC )
       CALL GC_CheckVar( 'drydep_mod:A_DEN', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       A_DEN(:)   = 0e+0_f8

       ! Aerosol radius [um]
       ALLOCATE( A_RADI( State_Chm%nDryDep ), STAT=RC )
       CALL GC_CheckVar( 'drydep_mod:A_RADI', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       A_RADI = 0e+0_f8

       ! Is the species an aerosol? (T/F)
       ALLOCATE( AIROSOL( State_Chm%nDryDep ), STAT=RC )
       CALL GC_CheckVar( 'drydep_mod:AIROSOL', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       AIROSOL = .FALSE.

       ! Drydep species name
       ALLOCATE( DEPNAME( State_Chm%nDryDep ), STAT=RC )
       CALL GC_CheckVar( 'drydep_mod:DEPNAME', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       DEPNAME = ''

       ! Reactivity factor
       ALLOCATE( F0( State_Chm%nDryDep ), STAT=RC )
       CALL GC_CheckVar( 'drydep_mod:F0', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       F0 = 0e+0_f8

       ! Henry's law K0
       ALLOCATE( HSTAR( State_Chm%nDryDep ), STAT=RC )
       CALL GC_CheckVar( 'drydep_mod:HSTAR', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       HSTAR = 0e+0_f8

       ! POPs KOA
       ALLOCATE( KOA( State_Chm%nDryDep ), STAT=RC )
       CALL GC_CheckVar( 'drydep_mod:KOA', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       KOA  = 0e+0_f8

       ! Drydep species indicies
       ALLOCATE( NDVZIND( State_Chm%nDryDep ), STAT=RC )
       CALL GC_CheckVar( 'drydep_mod:NDVZIND', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       NDVZIND = 0

       ! Drydep scaling flag
       ALLOCATE( FLAG( State_Chm%nDryDep ), STAT=RC )
       CALL GC_CheckVar( 'drydep_mod:FLAG', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       FLAG = 0

       ! Species indices
       ALLOCATE( NTRAIND( State_Chm%nDryDep ), STAT=RC )
       CALL GC_CheckVar( 'drydep_mod:NTRAIND', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       NTRAIND = 0

       ! Molecular weight [kg/mol]
       ALLOCATE( XMW( State_Chm%nDryDep ), STAT=RC )
       CALL GC_CheckVar( 'drydep_mod:XMW', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       XMW = 0e+0_f8

    ENDIF

    !=================================================================
    ! First identify species that dry deposit and then initialize
    ! the dry deposition quantities accordingly:
    !
    ! Quantity  Description
    ! ----------------------------------------------------------------
    ! NUMDEP    Number of dry depositing species
    ! NTRAIND   GEOS-Chem species ID number (advected index)
    ! NDVZIND   Coresponding index in the DVEL drydep velocity array
    ! HSTAR     Henry's law solubility constant [M atm-1]
    ! F0        Reactivity (0.0 = not reactive, 1.0=very reactive)
    ! XMW       Molecular weight of species [kg mol-1]
    ! AIROSOL    = T if the species is aerosol, =F if gas
    ! A_DEN     Aerosol density [kg m-3]
    ! A_RADIUS  Aerosol radius [um]
    ! KOA       POPs KOA parameter
    !
    ! NOTES:
    ! (1) For XMW, we take the species emitted molecular weight
    !      (EmMW_g * 1e-3).  Several hydrocarbons are transported
    !      as equivalent molecules of carbon.  In this case the
    !      EmMw_g will be 12.0.
    ! (2) We have edited the molecular weights of some species to
    !      match the prior code.  Some of these definitions are
    !      inconsistent with the listed molecular weights. (e.g.
    !      ACET is transported as 3 carbons, but we give it a
    !      molecular weight of 58 instead of 12).  These will need
    !      to be researched further.
    ! (3) The deposition names of SO4s and NITs need to be in
    !      uppercase.  Therefore, we overwrite the values from
    !      the species database with SO4S, NITS.
    !=================================================================
    DO N = 1, State_Chm%nAdvect

       ! Point to the Nth species in the species database
       SpcInfo => State_Chm%SpcData(N)%Info

       ! Only proceed if the species dry deposits
       IF ( SpcInfo%Is_Drydep ) THEN

          ! Initialize dry deposition quantities
          NUMDEP            = NUMDEP + 1
          NTRAIND(NUMDEP)   = SpcInfo%ModelID
          NDVZIND(NUMDEP)   = SpcInfo%DryDepID
          DEPNAME(NUMDEP)   = TRIM( SpcInfo%Name )
          HSTAR(NUMDEP)     = DBLE( SpcInfo%DD_Hstar_old   )
          F0(NUMDEP)        = DBLE( SpcInfo%DD_F0          )
          KOA(NUMDEP)       = DBLE( SpcInfo%DD_KOA         )
          XMW(NUMDEP)       = DBLE( SpcInfo%MW_g * 1e-3_fp )
          AIROSOL(NUMDEP)   = ( .not. SpcInfo%Is_Gas       )

          ! Only copy DENSITY if it's not a missing value
          IF ( SpcInfo%Density > 0.0_fp ) THEN
             A_DEN(NUMDEP)  = DBLE( SpcInfo%Density        )
          ENDIF

          ! Only copy RADIUS if it's not a missing value
          IF ( SpcInfo%Radius > 0.0_fp ) THEN
             A_RADI(NUMDEP) = DBLE( SpcInfo%Radius         )
          ENDIF

          !-----------------------------------------------------
          ! Kludges to match behavior of older code
          !-----------------------------------------------------
          SELECT CASE ( TRIM( SpcInfo%Name ) )

          CASE( 'ACET' )
             ! Flag the species ID of ACET for use above.
             id_ACET = SpcInfo%ModelId

          CASE( 'O3' )
             ! Flag the species ID of O3 for use above
             ID_O3 = SpcInfo%ModelId

          CASE( 'ALD2' )
             ! Flag the species ID of ALD2 for use above.
             id_ALD2 = SpcInfo%ModelId

          CASE( 'MENO3' )
             ! Flag the species ID of MENO3 for use above.
             id_MENO3 = SpcInfo%ModelId

          CASE( 'ETNO3' )
             ! Flag the species ID of ETNO3 for use above.
             id_ETNO3 = SpcInfo%ModelId

          CASE( 'NITs', 'NITS' )
             ! DEPNAME for NITs has to be in all caps, for
             ! backwards compatibility with older code.
             DEPNAME(NUMDEP) = 'NITS'

          CASE( 'N2O5', 'HC187' )
             ! These species scale to the Vd of HNO3. We will
             ! explicitly compute the Vd of these species instead
             ! of assigning the Vd of HNO3 from the DVEL array.
             ! The scaling is applied in DO_DRYDEP using FLAG=1.
             !
             ! Make sure to set XMW to the MW of HNO3
             ! for the computation of Vd to work properly.
             XMW(NUMDEP)  = State_Chm%SpcData(id_HNO3)%Info%MW_g * 1e-3_fp
             FLAG(NUMDEP) = 1

          CASE(  'MPAN', 'PPN', 'R4N2' )
             ! These specied scale to the Vd of PAN.  We will
             ! explicitly compute the Vd of these species instead
             ! of assigning the Vd of PAN from the DVEL array.
             ! The scaling is applied in DO_DRYDEP using FLAG=2.
             !
             ! Make sure to set XMW to the MW of PAN
             ! for the computation of Vd to work properly.
             XMW(NUMDEP)  = State_Chm%SpcData(id_PAN)%Info%MW_g * 1e-3_fp
             FLAG(NUMDEP) = 2

          CASE( 'MONITS', 'MONITU', 'HONIT' )
             ! These species scale to the Vd of ISOPN. We will
             ! explicitly compute the Vd of these species instead
             ! of assigning the Vd of ISOPN from the DVEL array.
             ! The scaling is applied in DO_DRYDEP using FLAG=3.
             !
             ! Make sure to set XMW to the MW of ISOPN
             ! for the computation of Vd to work properly.
             XMW(NUMDEP)  = State_Chm%SpcData(id_IHN1)%Info%MW_g * 1e-3_fp
             FLAG(NUMDEP) = 3

          CASE( 'SO4s', 'SO4S' )
             ! DEPNAME for SO4s has to be in all caps, for
             ! backwards compatibility with older code
             DEPNAME(NUMDEP) = 'SO4S'

          CASE DEFAULT
             ! Do nothing

          END SELECT

       ENDIF

       ! Free pointer
       SpcInfo => NULL()
    ENDDO

    ! For TOMAS
    id_NK1   = Ind_('NK1'  )

    !=================================================================
    ! Allocate arrays
    ! add allocation for SALT_V and DMID (jaegle 5/11/11)
    !=================================================================
    ALLOCATE( SALT_V( NR_MAX ), STAT=RC )
    CALL GC_CheckVar( 'drydep_mod:SALT_V', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    SALT_V = 0e+0_f8

    ALLOCATE( DMID( NR_MAX ), STAT=RC )
    CALL GC_CheckVar( 'drydep_mod:DMID', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    DMID = 0e+0_f8

    !=================================================================
    ! Echo information to stdout
    !=================================================================
    IF ( Input_Opt%amIRoot ) THEN

       ! Line 1
       MSG = 'INIT_DRYDEP: List of dry deposition species:'
       WRITE( 6, '(/,a)' ) TRIM( MSG )

       ! Line 2
       MSG =  '  #   Name      Species DEPVEL Henry''s    React.' &
           // '   Molec.   Aerosol?'
       WRITE( 6, '(/,a)'   ) TRIM( MSG )

       ! Line 3
       MSG =  '                Number Index  Law Const  Factor' &
           // '   Weight   (T or F)'
       WRITE( 6, '(a)'   ) TRIM( MSG )

       ! Separator
       WRITE( 6, '(a)'   ) REPEAT( '-', 70 )

       ! Output
       DO N = 1, NUMDEP
          WRITE( 6, 100 ) N,          ADJUSTL( DEPNAME(N) ), &
                          NTRAIND(N), NDVZIND(N),            &
                          HSTAR(N),   F0(N),                 &
                          XMW(N),     AIROSOL(N)
       ENDDO
100    FORMAT( i3, 3x, a8, 2(3x,i3), 4x, es8.1, 2(3x,f6.3), 3x, L3 )

    ENDIF

    ! Calls INIT_WEIGHTSS to calculate the volume distribution of
    ! sea salt aerosols (jaegle 5/11/11)
    CALL INIT_WEIGHTSS( Input_Opt )

  END SUBROUTINE INIT_DRYDEP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_drydep
!
! !DESCRIPTION: Subroutine CLEANUP\_DRYDEP deallocates all module arrays.
!  (bmy, 2/27/03, 2/22/05)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_DRYDEP
!
! !REVISION HISTORY:
!  27 Feb 2003 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! CLEANUP_DRYDEP begins here!
    !=================================================================
    IF ( ALLOCATED( A_DEN    ) ) DEALLOCATE( A_DEN    )
    IF ( ALLOCATED( A_RADI   ) ) DEALLOCATE( A_RADI   )
    IF ( ALLOCATED( AIROSOL  ) ) DEALLOCATE( AIROSOL  )
    IF ( ALLOCATED( DEPNAME  ) ) DEALLOCATE( DEPNAME  )
    IF ( ALLOCATED( DMID     ) ) DEALLOCATE( DMID     )
    IF ( ALLOCATED( F0       ) ) DEALLOCATE( F0       )
    IF ( ALLOCATED( FLAG     ) ) DEALLOCATE( FLAG     )
    IF ( ALLOCATED( HSTAR    ) ) DEALLOCATE( HSTAR    )
    IF ( ALLOCATED( KOA      ) ) DEALLOCATE( KOA      )
    IF ( ALLOCATED( NDVZIND  ) ) DEALLOCATE( NDVZIND  )
    IF ( ALLOCATED( NTRAIND  ) ) DEALLOCATE( NTRAIND  )
    IF ( ALLOCATED( SALT_V   ) ) DEALLOCATE( SALT_V   )
    IF ( ALLOCATED( XMW      ) ) DEALLOCATE( XMW      )

  END SUBROUTINE CLEANUP_DRYDEP
!EOC
END MODULE DRYDEP_MOD
