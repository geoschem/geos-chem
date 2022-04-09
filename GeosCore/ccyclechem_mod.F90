!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: global_ch4_co_co2_mod.F
!
! !DESCRIPTION: Module CCYLECHEM_MOD contains variables and routines
! for simulating CH4 CO and CO2 with an online calculation of the
! chemistry between them using KPP. It was adapted directly from
! the module CH4_CO_CO2_MOD.F.
!

!\\
!\\
! !INTERFACE: 
!
      MODULE CCYCLECHEM_MOD
!
! !USES:
!
      USE HCO_ERROR_MOD       ! For HEMCO error reporting
      USE PhysConstants
      USE PRECISION_MOD       ! For GEOS-Chem Precision (fp, f4, f8)
      USE inquireMod,    ONLY : findFreeLUN

      IMPLICIT NONE
      PRIVATE

!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: CALC_DIURNAL_CH4COCO2
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: EMISS_CH4COCO2
      PUBLIC :: CHEM_CH4COCO2
      PUBLIC :: INIT_CH4COCO2
      PUBLIC :: CLEANUP_CH4COCO2
!
! !REVISION HISTORY:
!  04 Apr 2022 - M.S. Long   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PUBLIC DATA MEMBERS:
!
      REAL(fp), PARAMETER,   PUBLIC :: XNUMOL_CH4 = AVO / 16d-3 ! hard-coded MW

      ! Make CH4_EMIS now public so that it can be used by vdiff_mod.F90
      ! Methane emissions units are [kg/m2/s]
      REAL(fp),  ALLOCATABLE, PUBLIC :: CH4_EMIS_J(:,:,:)

!------------------------------------------------------------------------------
!BOC
!     
! !DEFINED PARAMETERS:
!
      !========================================================================
      ! Module Variables:
      ! BAIRDENS   : Array for air density                      [molec/cm3]
      ! BOH        : Array for OH values                        [kg/m3]
      ! BCl        : Array for Cl values                        [ppbv]
      ! COPROD     : Array for zonal mean P(CO)                 [v/v/s]
      ! CH4LOSS    : Array for monthly average CH4 loss freq    [1/s]
      ! FMOL_CH4   : Molecular weight of CH4                    [kg/mole]
      ! XNUMOL_CH4 : Molecules CH4 / kg CH4                     [molec/kg]
      ! CH4_EMIS   : Array for CH4 Emissions                    [kg/m2/s]
      !========================================================================

      REAL(fp), PARAMETER   :: XNUMOL_OH = AVO / 17e-3_fp  ! molec OH / kg OH
                                                           ! hard-coded MW
      REAL(fp), PARAMETER   :: CM3PERM3  = 1.e+6_fp
!
! !LOCAL VARIABLES:
!
      ! Diagnostic flags
      LOGICAL               :: Do_ND43

      ! Species ID flag
      INTEGER               :: id_CH4

      ! Various arrays      
      REAL(fp), ALLOCATABLE :: BAIRDENS(:,:,:)

      ! Pointers to fields in the HEMCO data structure.
      ! These need to be declared as REAL(f4), aka REAL*4.
      ! NOTE: These are globally SAVEd variables so we can
      ! nullify these in the declaration statement (bmy, 4/29/16)
      REAL(f4), POINTER     :: BOH    (:,:,:) => NULL()
      REAL(f4), POINTER     :: BCl    (:,:,:) => NULL()
      REAL(f4), POINTER     :: CH4LOSS(:,:,:) => NULL()
      REAL(f4), POINTER     :: CLUSTERS(:,:)  => NULL()

      REAL(fp)              :: TROPOCH4

! !PRIVATE TYPES:
!
      ! Arrays
      REAL(fp), ALLOCATABLE :: SUMACETCO  (:,:  )   ! P(CO) from Acetone
      REAL(fp), ALLOCATABLE :: SUMCH3OHCO (:,:  )   ! P(CO) from CH3OH
      REAL(fp), ALLOCATABLE :: SUMISOPCO  (:,:  )   ! P(CO) from Isoprene
      REAL(fp), ALLOCATABLE :: SUMMONOCO  (:,:  )   ! P(CO) from Monoterpenes
      REAL(fp), ALLOCATABLE :: TCOSZ      (:,:  )   ! Daily sum COS(SZA)
      REAL(fp), ALLOCATABLE :: CH4_OH_TROP(:,:,:)   ! CH4 loss in trop
      REAL(fp), ALLOCATABLE :: CH4_OH_STRAT(:,:,:)   ! CH4 loss in strat 

      ! Pointers to fields in the HEMCO data structure.
      ! These need to be declared REAL(f4), aka REAL*4.
      REAL(f4), POINTER     :: OH         (:,:,:) => NULL() ! Global OH
      REAL(f4), POINTER     :: GMI_PROD_CO(:,:,:) => NULL() ! Global P(CO)
      REAL(f4), POINTER     :: GMI_LOSS_CO(:,:,:) => NULL() ! Global L(CO)
      REAL(f4), POINTER     :: PCO_CH4    (:,:,:) => NULL() ! CH4 P(CO)
      REAL(f4), POINTER     :: PCO_NMVOC  (:,:,:) => NULL() ! NMVOC P(CO)
      REAL(f4), POINTER     :: SFC_CH4    (:,:  ) => NULL() ! Global sfc
CH4
!
! !DEFINED PARAMETERS:
!
      ! FMOL_CO2     - kg CO2 / mole CO2 
      REAL(fp),  PARAMETER   :: FMOL_CO2   = 44e-3_fp
      ! FMOL_C       - kg C   / mole C 
      REAL(fp),  PARAMETER   :: FMOL_C     = 12e-3_fp
      ! FMOL_CO       - kg CO   / mole CO 
      REAL(fp),  PARAMETER   :: FMOL_CO    = 28e-3_fp


      ! XNUMOL_CO2   - molecules CO2 / kg CO2 
      REAL(fp),  PARAMETER   :: XNUMOL_CO2 = AVO / FMOL_CO2
      ! XNUMOL_C     - molecules C / kg C
      REAL(fp),  PARAMETER   :: XNUMOL_C = AVO / FMOL_C
      ! XNUMOL_CO     - molecules CO / kg CO
      REAL(fp),  PARAMETER   :: XNUMOL_CO = AVO / FMOL_CO


      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emissch4
!
! !DESCRIPTION: Subroutine EMISSCH4 places emissions of CH4 [kg] into the
!  chemical species array.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_CH4COCO2( am_I_Root, Input_Opt,
     &                     State_Met, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE HCO_INTERFACE_MOD,  ONLY : HcoState, GetHcoDiagn
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REMARKS:
!  WARNING: Soil absorption has to be the 15th field in CH4_EMIS
!  Also: the ND58 diagnostics have now been removed.  We still need to
!  read the HEMCO manual diagnostics into CH4_EMIS for the analytical
!  inversion.  Therefore, we will keep EmissCh4 for the time-being
!  but only remove the bpch diagnostic.
! 
! !REVISION HISTORY:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) EMISSCH4 is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) GLOBASEAEMIS, GLOBSEAEMIS are diagnostics by jsw.
!  (4 ) Do not multiply CO emissions by 1.28 anymore (jsw, bmy, 2/12/01)
!  (5 ) Renamed input files to CH4_monthly.geos.{RES} and 
!        CH4_aseasonal.geos.{RES}. (bmy, 2/12/01)
!  (6 ) Add reference to "CMN_SETUP" for the DATA_DIR variable (bmy, 2/13/01)
!  (7 ) Removed references to "biofuel_mod.f" and "biomass_mod.f"; these
!        weren't necessary (bmy, 3/20/01)
!  (8 ) Now reference IU_FILE and IOERROR from "file_mod.f".  Now use IU_FILE
!        instead of IUNIT as the file unit #. (bmy, 6/27/02)
!  (9 ) Now reference BXHEIGHT and SUNCOS from "dao_mod.f".  Remove reference 
!        to header file "comtrid.h" -- it's not used.  Make FIRSTEMISS a local
!        SAVEd variable.  Also use MONTH from "CMN" instead of the variable
!        LMN. (bmy, 11/15/02)
!  (10) Now replace DXYP(JREF)*1d4 with routine GET_AREA_CM2 of "grid_mod.f".
!        Now use function GET_MONTH and GET_TS_EMIS from "time_mod.f". 
!        Now use functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f".
!        I0 and J0 are now local variables. (bmy, 3/27/03)
!  (11) Now reference STT from "tracer_mod.f".  Now reference DATA_DIR from
!        "directory_mod.f". (bmy, 7/20/04)
!  (12) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (13) Add non-local PBL capability (ccc, 8/31/09)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  11 Apr 2014 - R. Yantosca - Remove call to INIT_GLOBAL_CH4
!  09 Sep 2014 - C. Keller   - Implemented HEMCO
!  11 Sep 2015 - E. Lundgren - Remove area-dependency (except for global mass
!                              sums) by outputting diagnostics as kg/m2 not kg
!  11 Sep 2015 - E. Lundgren - Remove State_Chm from arg list since not used
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER                  :: I, J, N
      REAL(fp)                 :: DTSRCE, AREA_M2

      ! Strings
      CHARACTER(LEN= 63)       :: DgnName
      CHARACTER(LEN=255)       :: ErrMsg
      CHARACTER(LEN=255)       :: ThisLoc

      ! For fields from Input_Opt
      LOGICAL                  :: ITS_A_CH4COCO2_SIM
      LOGICAL                  :: LPRT
      LOGICAL                  :: prtDebug
      LOGICAL, SAVE            :: FIRST = .TRUE.

      ! Pointers
      REAL(f4),        POINTER :: Ptr2D(:,:)

      !=================================================================
      ! EMISSCH4 begins here!
      !=================================================================
      ! Nullify pointers
      Ptr2D => NULL()

      ! Assume success
      RC      = GC_SUCCESS
      ErrMsg  = ''
      ThisLoc = ' -> at EMISSCH4 (in GeosCore/global_ch4_mod.F)'

      ! Copy values from Input_Opt
      ITS_A_CH4COCO2_SIM  = Input_Opt%ITS_A_CH4COCO2_SIM
      LPRT           = Input_Opt%LPRT

      ! Do we have to print debug output?
      prtDebug       = ( LPRT .and. am_I_Root )

      ! Exit with error if we can't find the HEMCO state object
      IF ( .NOT. ASSOCIATED( HcoState ) ) THEN
         ErrMsg = 'The HcoState object is undefined!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Emission timestep
      DTSRCE = HcoState%TS_EMIS

      IF ( ITS_A_CH4COCO2_SIM .and. LPRT ) THEN
         print*,'BEGIN SUBROUTINE: EMISSCH4' 
      ENDIF

      ! =================================================================
      ! --> All emission calculations are now done through HEMCO
      ! HEMCO stores emissions of all species internally in the HEMCO
      ! state object. Here, we pass these emissions into module array 
      ! CH4_EMIS in units kg/m2/s. These values are then either added to
      ! the species array (full mixing scheme) or used later on in
      ! vdiff_mod.F90 if the non-local PBL mixing scheme is used.
      !
      ! The CH4_EMIS array is mostly used for backwards compatibility
      ! (especially the diagnostics). It is also used to ensure that 
      ! in a multi-species simulation, species 1 (total CH4) is properly
      ! defined. 
      !
      ! NOTE: the ND60 diagnostics (wetland fraction) is currently not
      ! written. To support this diagnostics, create a manual diagnostics
      ! in the wetland extension (HEMCO/Extensions/hcox_ch4wetland_mod.F90),
      ! activate it in hcoi\_gc\_diagn\_mod.F90 and import it here.
      !
      !                                              (ckeller, 9/12/2013)
      ! =================================================================
      CH4_EMIS_J(:,:,:) = 0e+0_fp

      !-------------------
      ! Oil
      !-------------------
      DgnName = 'CH4_OIL'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

      ! Trap potential errors
      IF ( RC /= HCO_SUCCESS ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF 

      ! Assign HEMCO pointer to array
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL HCO_WARNING( ErrMsg, RC, ThisLoc=ThisLoc )
      ELSE
         CH4_EMIS_J(:,:,2) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      !-------------------
      ! Gas
      !-------------------
      DgnName = 'CH4_GAS'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

      ! Trap potential errors
      IF ( RC /= HCO_SUCCESS ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF 

      ! Assign HEMCO pointer to array
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL HCO_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
      ELSE
         CH4_EMIS_J(:,:,3) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()


      !-------------------
      ! Coal 
      !-------------------
      DgnName = 'CH4_COAL'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

      ! Trap potential errors
      IF ( RC /= HCO_SUCCESS ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF 

      ! Assign HEMCO pointer to array
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL HCO_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
      ELSE
         CH4_EMIS_J(:,:,4) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      !-------------------
      ! Livestock
      !-------------------
      DgnName = 'CH4_LIVESTOCK'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

      ! Trap potential errors
      IF ( RC /= HCO_SUCCESS ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF 

      ! Assign HEMCO pointer to array
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL GC_Error( ErrMsg, RC, ThisLoc=ThisLoc )
      ELSE
         CH4_EMIS_J(:,:,5) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      !-------------------
      ! Landfills
      !-------------------
      DgnName = 'CH4_LANDFILLS'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

      ! Trap potential errors
      IF ( RC /= HCO_SUCCESS ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF 

      ! Assign HEMCO pointer to array
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL HCO_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
      ELSE
         CH4_EMIS_J(:,:,6) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      !-------------------
      ! Wastewater
      !-------------------
      DgnName = 'CH4_WASTEWATER'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

      ! Trap potential errors
      IF ( RC /= HCO_SUCCESS ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF 

      ! Assign HEMCO pointer to array
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL HCO_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
      ELSE
         CH4_EMIS_J(:,:,7) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      !-------------------
      ! Rice 
      !-------------------
      DgnName = 'CH4_RICE'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

      ! Trap potential errors
      IF ( RC /= HCO_SUCCESS ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF 

      ! Assign HEMCO pointer to array
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL HCO_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
      ELSE
         CH4_EMIS_J(:,:,8) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      !-------------------
      ! Other anthropogenic 
      !-------------------
      DgnName = 'CH4_ANTHROTHER'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

      ! Trap potential errors
      IF ( RC /= HCO_SUCCESS ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF 

      ! Assign HEMCO pointer to array
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL HCO_Warning( ErrMsg, RC, ThisLoc = ThisLoc )
      ELSE
         CH4_EMIS_J(:,:,9) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      !-------------------
      ! Biomass burning 
      !-------------------
!      DgnName = 'BIOMASS_CH4'
      DgnName = 'CH4_BIOMASS'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

      ! Trap potential errors
      IF ( RC /= HCO_SUCCESS ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF 

      ! Assign HEMCO pointer to array
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL HCO_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
      ELSE
         CH4_EMIS_J(:,:,10) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      !-------------------
      ! Wetland 
      !-------------------
      DgnName = 'CH4_WETLAND'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

      ! Trap potential errors
      IF ( RC /= HCO_SUCCESS ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF 

      ! Assign HEMCO pointer to array
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL HCO_WARNING( ErrMsg, RC, ThisLoc=ThisLoc )
      ELSE
         CH4_EMIS_J(:,:,11) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      !-------------------
      ! Global seeps
      !-------------------
      DgnName = 'CH4_SEEPS'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

      ! Trap potential errors
      IF ( RC /= HCO_SUCCESS ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF 

      ! Assign HEMCO pointer to array
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL HCO_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
      ELSE
         CH4_EMIS_J(:,:,12) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      !-------------------
      ! Lakes
      !-------------------
      DgnName = 'CH4_LAKES'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

      ! Trap potential errors
      IF ( RC /= HCO_SUCCESS ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF 

      ! Assign HEMCO pointer to array
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL HCO_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
      ELSE
         CH4_EMIS_J(:,:,13) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      !-------------------
      ! Termites
      !-------------------
      DgnName = 'CH4_TERMITES'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

      ! Trap potential errors
      IF ( RC /= HCO_SUCCESS ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF 

      ! Assign HEMCO pointer to array
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL HCO_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
      ELSE
         CH4_EMIS_J(:,:,14) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      !-------------------
      ! Soil absorption (those are negative!) 
      !-------------------
      DgnName = 'CH4_SOILABSORB'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

      ! Trap potential errors
      IF ( RC /= HCO_SUCCESS ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF 

      ! Assign HEMCO pointer to array
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
         CALL HCO_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
      ELSE
         CH4_EMIS_J(:,:,15) =  Ptr2D(:,:) * -1.0_hp
      ENDIF
      Ptr2D => NULL()

      ! =================================================================
      ! Total emission: sum of all emissions - (2*soil absorption)
      ! We have to substract soil absorption twice because it is added 
      ! to other emissions in the SUM function. (ccc, 7/23/09)
      ! =================================================================
      CH4_EMIS_J(:,:,1) = SUM(CH4_EMIS_J, 3) - (2 * CH4_EMIS_J(:,:,15))

      IF ( prtDebug ) THEN
         WRITE(*,*) 'CH4_EMIS (kg/m2/s):'
         WRITE(*,*) 'Total        : ', SUM(CH4_EMIS_J(:,:,1))
         WRITE(*,*) 'Oil          : ', SUM(CH4_EMIS_J(:,:,2))
         WRITE(*,*) 'Gas          : ', SUM(CH4_EMIS_J(:,:,3))
         WRITE(*,*) 'Coal         : ', SUM(CH4_EMIS_J(:,:,4))
         WRITE(*,*) 'Livestock    : ', SUM(CH4_EMIS_J(:,:,5))
         WRITE(*,*) 'Landfills    : ', SUM(CH4_EMIS_J(:,:,6))
         WRITE(*,*) 'Wastewater   : ', SUM(CH4_EMIS_J(:,:,7))
         WRITE(*,*) 'Rice         : ', SUM(CH4_EMIS_J(:,:,8))
         WRITE(*,*) 'Other anth   : ', SUM(CH4_EMIS_J(:,:,9))
         WRITE(*,*) 'Biomass burn : ', SUM(CH4_EMIS_J(:,:,10))
         WRITE(*,*) 'Wetlands     : ', SUM(CH4_EMIS_J(:,:,11))
         WRITE(*,*) 'Seeps        : ', SUM(CH4_EMIS_J(:,:,12))
         WRITE(*,*) 'Lakes        : ', SUM(CH4_EMIS_J(:,:,13))
         WRITE(*,*) 'Termites     : ', SUM(CH4_EMIS_J(:,:,14))
         WRITE(*,*) 'Soil absorb  : ', SUM(CH4_EMIS_J(:,:,15))
      ENDIF

      IF ( ITS_A_CH4COCO2_SIM .and. LPRT ) THEN
         print*,'END SUBROUTINE: EMISSCH4'
      ENDIF

      END SUBROUTINE EMISS_CH4COCO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chemch4
!
! !DESCRIPTION: Subroutine CHEMCH4 computes the chemical loss of CH4
!  (sources - sinks). (jsw, bnd, bmy, 6/8/00, 10/3/05)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_CH4COCO2( am_I_Root, Input_Opt,  State_Met,
     &                    State_Chm, State_Diag, RC         )
!
! !USES:
!
      USE CMN_SIZE_MOD
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD
      USE DIAG_MOD,           ONLY : AD43, AD65
      USE DIAG04_MOD,         ONLY : ND04
      USE DIAG04_MOD,         ONLY : AD04_chem
#endif
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : CHECK_VALUE
      USE GC_GRID_MOD,        ONLY : GET_YMID
      USE HCO_INTERFACE_MOD,  ONLY : HcoState, GetHcoID
      USE HCO_Error_Mod
      USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PhysConstants,      ONLY : AVO
      USE State_Chm_Mod,      ONLY : ChmState, Ind_
      USE State_Diag_Mod,     ONLY : DgnState
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM,     GET_TS_EMIS
      USE TIME_MOD,           ONLY : GET_MONTH,       GET_YEAR
      USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH, ITS_A_NEW_YEAR
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!  CH4 SOURCES
!  ============================================================================
!  (1 ) Oxidation of methane, isoprene and monoterpenes (SRCO_fromHCs).
!  (2 ) Direct emissions of CO from fossil fuel combustion, biomass 
!        burning and wood (for fuel) burning (SR SETEMIS).
!  (3 ) Emissions.
!                                                                             .
!  CH4 SINKS:
!  ============================================================================
!  (1 ) Removal of CO by OH (SR OHparam & CO_decay).
!  (2 ) CO uptake by soils (neglected).
!  (3 ) Transport of CO to stratosphere from troposphere 
!        (in dynamical subroutines).
!  (4 ) Removal by OH (Clarissa's OH--climatol_OH.f and CO_decay.f)
!  (5 ) Transport of CH4 between troposphere and stratosphere, and 
!        destruction in strat (CH4_strat.f).
!  (6 ) Removel by Cl
! 
! !REVISION HISTORY:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (6/8/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CHEMCH4 is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Updated comments (jsw, bmy, 2/12/01)
!  (4 ) LD43 is already declared in CMN_DIAG; don't redefine it (bmy, 11/15/01)
!  (5 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (6 ) Now reference AD from "dao_mod.f".  Now reference GEOS_CHEM_STOP from
!        "error_mod.f"  Now make FIRSTCHEM a local SAVEd variable.  Now 
!        reference ALBD from "dao_mod.f".  Now use MONTH and JDATE from "CMN"
!        instead of LMN and LDY. (bmy, 11/15/02)
!  (7 ) Remove NYMDb, NYMDe from the arg list.  Now use functions GET_MONTH,
!        GET_NYMDb, GET_NYMDe, GET_MONTH, GET_DAY from the new "time_mod.f"
!        (bmy, 3/27/03) 
!  (8 ) Now reference DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (9 ) Remove reference to BPCH2_MOD, it's not needed (bmy, 10/3/05)
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  23 Oct 2013 - R. Yantosca - Now pass objects to GET_GLOBAL_OH routine
!  12 Feb 2014 - K. Wecht    - Disable CH4 budget diagnostic (bracket the 
!                              code out with #ifdef blocks so it can be used)
!  23 Jul 2014 - R. Yantosca - Remove reference to obsolete CMN_mod.F
!  23 Jul 2014 - R. Yantosca - Reference ITS_IN_THE_CHEMGRID (chemgrid_mod.F)
!  24 Jul 2014 - R. Yantosca - Now compute BOXVL internally
!  17 Sep 2014 - C. Keller   - Now use HEMCO to get CH4 loss and OH conc. field.
!  06 Jan 2016 - E. Lundgren - Use global physical parameters
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  10 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!  03 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE      :: FIRSTCHEM = .TRUE.
      INTEGER            :: I, J, L

      REAL(fp)           :: PREVCH4(IIPAR, JJPAR, LLPAR)

      ! Number of days per month
      INTEGER            :: NODAYS(12) = (/ 31, 28, 31, 30, 
     &                                      31, 30, 31, 31, 
     &                                      30, 31, 30, 31 /)

      ! For fields from Input_Opt
      LOGICAL            :: LSPLIT
      LOGICAL            :: LPRT

      ! Pointers
      REAL(fp), POINTER  :: Spc(:,:,:,:)

      ! Scalars
      LOGICAL                  :: FOUND
      LOGICAL                  :: LPCO_CH4,      LPCO_NMVOC
      INTEGER                  :: nAdvect,       HcoID,       NA
      INTEGER                  :: N,             MONTH,       YEAR
      REAL(fp)                 :: ALPHA_ISOP,    DTCHEM,      GCO
      REAL(fp)                 :: STTCO,         KRATE,       CH4
      REAL(fp)                 :: CO_CH4,        CO_ISOP,     CO_MONO
      REAL(fp)                 :: CO_CH3OH,      CO_OH,       CO_ACET
      REAL(fp)                 :: CO_NMVOC,      E_CO2
      REAL(fp)                 :: CH4RATE,       DENS,        CORATE
      REAL(fp)                 :: YMID,          BOXVL,       FMOL_CO
      REAL(fp)                 :: KHI1,          KLO1,        XYRAT1
      REAL(fp)                 :: BLOG1,         FEXP1,       KHI2
      REAL(fp)                 :: KLO2,          XYRAT2,      BLOG2
      REAL(fp)                 :: FEXP2,         KCO1,        KCO2
      REAL(fp)                 :: OH_MOLEC_CM3,  FAC_DIURNAL, SUNCOS
      REAL(fp)                 :: PCO_CH4_MCM3S, PCO_NMVOC_MCM3S
      REAL(fp)                 :: kgs_to_atomsC, Emis,        DTEMIS
      REAL(fp)                 :: kgm3_to_mcm3OH,kgm3_to_mcm3sCO

      ! Strings
      CHARACTER(LEN=255)       :: ERR_VAR
      CHARACTER(LEN=255)       :: ERR_MSG
      CHARACTER(LEN=63)        :: DgnName
      CHARACTER(LEN=255)       :: ErrMsg
      CHARACTER(LEN=255)       :: ThisLoc

      ! Arrays
      INTEGER                  :: ERR_LOC(4)

      ! Pointers
      REAL(fp),        POINTER :: AD    (:,:,:  )
      REAL(fp),        POINTER :: AIRVOL(:,:,:  )
      !REAL(fp),        POINTER :: Spc   (:,:,:,:)
      REAL(fp),        POINTER :: T     (:,:,:  )

      ! SAVED scalars
      LOGICAL,   SAVE          :: FIRST = .TRUE.
      REAL(fp),  SAVE          :: A3090S, A0030S, A0030N, A3090N
      INTEGER,   SAVE          :: IDch4   = -1
      INTEGER,   SAVE          :: IDnmvoc = -1
      INTEGER,   SAVE          :: IDisop  = -1
      INTEGER,   SAVE          :: IDch3oh = -1
      INTEGER,   SAVE          :: IDmono  = -1
      INTEGER,   SAVE          :: IDacet  = -1
!
! !DEFINED PARAMETERS:
!
      ! Switch to scale yield of isoprene from NOx concentration or not
      LOGICAL,   PARAMETER     :: ALPHA_ISOP_FROM_NOX = .FALSE.
      ! Yield of CO from CH4
      REAL(fp),  PARAMETER     :: ALPHA_CH4  = 1.0_fp

      ! Yield of CO from monoterpenes
      REAL(fp),  PARAMETER     :: ALPHA_MONO = 2e-1_fp

      ! Yield of CO from acetone
      REAL(fp),  PARAMETER     :: ALPHA_ACET = 2.0_fp / 3.0_fp
      ! For CO2 conversion
      REAL(fp),        PARAMETER :: CM2PERM2 = 1.d4


      !=================================================================
      ! CHEMCH4 begins here!
      !=================================================================
      ! Assume success
      RC      = GC_SUCCESS
      ErrMsg  = ''
      ThisLoc = ' -> at CHEMCH4 (in module GeosCore/ch4coco2_mod.F)'

      ! Copy values from Input_Opt
      LSPLIT  = Input_Opt%LSPLIT
      LPRT    = Input_Opt%LPRT

      ! Point to the chemical species 
      Spc     => State_Chm%Species

      ! DTCHEM is the number of seconds per chemistry timestep
      DTCHEM    = GET_TS_CHEM()

      IF ( LPRT ) THEN
         WRITE( 6, '(a)' ) '% --- ENTERING CHEMCH4! ---'
      ENDIF

      !=================================================================
      ! (0) Calculate each box's air density [molec/cm3]
      !     Do this for saving mean OH concentrations in CH4_DECAY
      !     and in CH4_OHSAVE (kjw, 6/12/09)
      !=================================================================

      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box volume [cm3]
         BOXVL           = State_Met%AIRVOL(I,J,L) * 1e+6_fp

         ! Air density [molec/cm3]
         BAIRDENS(I,J,L) = State_Met%AD(I,J,L) * 1000e+0_fp   /
     &                     BOXVL               * AVO / AIRMW

      ENDDO
      ENDDO
      ENDDO

      !================================================================
      ! (1) Get CH4 loss rates from HEMCO. the target is automatically 
      ! updated by HEMCO (ckeller, 9/16/2014)
      !================================================================
      IF ( FIRSTCHEM ) THEN

         ! Import CH4 loss frequencies from HEMCO. The target will be 
         ! updated automatically by HEMCO (ckeller, 9/16/2014)
         CALL HCO_GetPtr( am_I_Root, HcoState, 'CH4_LOSS', CH4LOSS, RC )

         ! Trap potential errors
         IF ( RC /= HCO_SUCCESS ) THEN
            ErrMsg = 'Cannot get pointer to HEMCO field CH4_LOSS!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF 
      ENDIF

      !================================================================
      ! (2) Get OH and Cl fields from HEMCO. The targets will be
      ! updated automatically by HEMCO (ckeller, 9/16/2014)
      !================================================================
      IF ( FIRSTCHEM ) THEN

         ! Get pointer to GLOBAL_OH
         CALL HCO_GetPtr( am_I_Root, HcoState, 'GLOBAL_OH_CH4', BOH, RC)

         ! Trap potential errors
         IF ( RC /= HCO_SUCCESS ) THEN
            ErrMsg = 'Cannot get pointer to HEMCO field GLOBAL_OH_CH4!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF 

         ! Get pointer to GLOBAL_Cl
         CALL HCO_GetPtr( am_I_Root, HcoState, 'GLOBAL_Cl', BCl, RC )

         ! Trap potential errors
         IF ( RC /= HCO_SUCCESS ) THEN
            ErrMsg = 'Cannot get pointer to HEMCO field GLOBAL_Cl!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF 
      ENDIF

      ! Compute diurnal cycle for OH every day (check for new day inside
      ! subroutine) - jaf 7/10/14
      CALL CALC_DIURNAL_CH4COCO2()

      ! Check HEMCO state object
      IF ( .NOT. ASSOCIATED(HcoState) ) THEN
         ErrMsg = 'HcoState object is not associated!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !=================================================================
      ! (3.1) ND43 diagnostics...save [OH] in molecules/cm3
      ! BOH is in kg/m3 (from HEMCO), convert to molecules/cm3 
      ! (ckeller, 9/16/2014)
      !=================================================================
      IF ( Do_ND43 .or. State_Diag%Archive_OHconcAfterChem ) THEN

         ! Zero the diagnostic array for OH to avoid leftover values
         IF ( State_Diag%Archive_OHconcAfterChem ) THEN
            State_Diag%OHconcAfterChem = 0.0_f4
         ENDIF

         ! NOTE: Can parallelize this loop later
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Same as for CO add a dirunal OH cycle (bb 2019)
            ! Cosine of the solar zenith angle [unitless]
            SUNCOS = State_Met%SUNCOSmid(I,J)

            ! Scaling factor for OH diurnal cycles - zero at night
            IF ( SUNCOS > 0.0_fp .and. TCOSZ(I,J) > 0.0_fp ) THEN
               FAC_DIURNAL = ( SUNCOS    / TCOSZ(I,J)    ) *
     &                       ( 86400.0_fp / GET_TS_CHEM() )
            ELSE
               FAC_DIURNAL = 0.0_fp
            ENDIF

            IF ( State_Met%InChemGrid(I,J,L) ) THEN

#if defined( BPCH_DIAG )
               !--------------------------------------------------------
               ! ND43 (bpch) diagnostic
               ! OH concentration in [molec/cm3] after chemistry
               !--------------------------------------------------------
               IF ( Do_ND43 .and. L <= LD43 ) THEN
                  AD43(I,J,L,1) = AD43(I,J,L,1) +
     &            ( BOH(I,J,L) * XNUMOL_OH / CM3PERM3 * FAC_DIURNAL )
               ENDIF
#endif


#if defined( NC_DIAG )
               !--------------------------------------------------------
               ! HISTORY (aka netCDF diagnostics)
               ! OH concentration in [molec/cm3] after chemistry
               !--------------------------------------------------------
               IF ( State_Diag%Archive_OHconcAfterChem ) THEN
                  State_Diag%OHconcAfterChem(I,J,L) =
     &             ( BOH(I,J,L) * XNUMOL_OH / CM3PERM3 * FAC_DIURNAL )
               ENDIF
#endif

            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      !=================================================================
      ! (4) Save OH concentrations for printing of global mean [OH] and
      !     CH3CCLl3 at end of simulation.
      !=================================================================
      CALL CH4_OHSAVE_CH4COCO2( State_Met, State_Chm )

      !=================================================================
      ! (5) If multi-CH4 species, we store the CH4 total conc. to
      !     distribute the sink after the chemistry. (ccc, 2/10/09)
      !================================================================= 
      IF ( LSPLIT ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            PREVCH4(I,J,L) = Spc(I,J,L,1)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF

      !=================================================================
      ! (6) calculate rate of decay of CH4 by OH oxidation.
      !=================================================================
      CALL CH4_DECAY_CH4COCO2( am_I_Root, Input_Opt,  State_Met, 
     &                State_Chm, State_Diag, RC         )


      !=================================================================
      ! (7) calculate CH4 chemistry in layers above tropopause.
      !=================================================================
      CALL CH4_STRAT_CH4COCO2( am_I_Root, Input_Opt,  State_Met, 
     &                State_Chm, State_Diag, RC         )

      !=================================================================
      ! (8) distribute the chemistry sink from total CH4 to other CH4 
      !     species. (ccc, 2/10/09)
      !=================================================================
      IF ( LSPLIT ) THEN
         CALL CH4_DISTRIB_CH4COCO2( PREVCH4, Input_Opt, State_Chm )
      ENDIF

      !=================================================================
      ! CHEM_TAGGED_CO begins here!
      !=================================================================
      ! Assume success
      !RC        = GC_SUCCESS
      !ErrMsg    = ''
      !ThisLoc   = ' -> at CHEM_TAGGED_CO (in GeosCore/tagged_co_mod.F)'

      ! Number of advected species
      nAdvect   = State_Chm%nAdvect

      ! Get fields from Input_Opt
      !LSPLIT     = Input_Opt%LSPLIT
      LPCO_CH4   = Input_Opt%LPCO_CH4
      LPCO_NMVOC = Input_Opt%LPCO_NMVOC

      ! DTCHEM is the number of seconds per chemistry timestep
      !DTCHEM    = GET_TS_CHEM()

      ! DTEMIS is the number of seconds per emission timestep
      ! used here for conversion from HEMCO
      DTEMIS    = GET_TS_EMIS()

      ! Initialize pointers
      AD        => State_Met%AD
      AIRVOL    => State_Met%AIRVOL
      Spc       => State_Chm%Species
      T         => State_Met%T

      ! Get the molecular weight of CO [kg/mol] from the species
      ! database
      ! (All tagged species are CO, so we can use the value for species
      ! #1)
      FMOL_CO   =  State_Chm%SpcData(16)%Info%MW_g * 1.0e-3_fp

      ! Factor to convert OH from kg/m3 (from HEMCO) to molec/cm3 
      kgm3_to_mcm3OH = ( AVO / 17.0e-3_fp ) * 1.0e-6_fp

      ! Factor to convert P(CO) from kg/m3 (from HEMCO) to molec/cm3/s
      ! HEMCO uses emission timestep to convert from kg/m3/s to kg/m3;
      ! apply that here.
      kgm3_to_mcm3sCO = ( AVO / (FMOL_CO * DTEMIS) ) * 1.0e-6_fp

      ! Zero diagnostic archival arrays to make sure that we don't have
      ! any
      ! leftover values from the last timestep near the top of the
      ! chemgrid
      IF ( State_Diag%Archive_Loss ) State_Diag%Loss = 0.0_f4
      IF ( State_Diag%Archive_ProdCOfromCH4 ) THEN
         State_Diag%ProdCOfromCH4 = 0.0_f4
      ENDIF

      IF ( State_Diag%Archive_ProdCOfromNMVOC ) THEN
         State_Diag%ProdCOfromNMVOC = 0.0_f4
      ENDIF

      !=================================================================
      ! Get pointers from HEMCO for OH, P(CO), and L(CO) fields
      ! 
      ! NOTES: 
      ! (1) We only need to get the pointers on the very first
      !      timestep.  HEMCO will update the targets automatically.
      ! (2) These calls have to be placed here instead of in routine
      !      INIT_TAGGED_CO.  This is because when INIT_TAGGED_CO is
      !      called, the HEMCO_Config file has not yet been read in.
      !=================================================================
      IF ( FIRST ) THEN

         ! Get species IDs
         IDch4    = Ind_( 'COch4'   )
         IDnmvoc  = Ind_( 'COnmvoc' )
         IDisop   = Ind_( 'COisop'  )
         IDch3oh  = Ind_( 'COch3oh' )
         IDmono   = Ind_( 'COmono'  )
         IDacet   = Ind_( 'COacet'  )

         ! Get a pointer to the OH field from the HEMCO list (bmy,
         ! 3/11/15)
         CALL HCO_GetPtr( am_I_Root, HcoState, 'GLOBAL_OH_CO', OH, RC )
         IF ( RC /= HCO_SUCCESS ) THEN
            ErrMsg = 'Cannot get pointer to GLOBAL_OH_CO!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

         ! Get pointer to strat P(CO) from GMI
         CALL HCO_GetPtr( am_I_Root, HcoState,
     &                   'GMI_PROD_CO', GMI_PROD_CO, RC )
         IF ( RC /= HCO_SUCCESS ) THEN
            ErrMsg = 'Cannot get pointer to GMI_PROD_CO!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

         ! Get pointer to strat L(CO) from GMI
         CALL HCO_GetPtr( am_I_Root, HcoState,
     &                   'GMI_LOSS_CO', GMI_LOSS_CO, RC )
         IF ( RC /= HCO_SUCCESS ) THEN
            ErrMsg = 'Cannot get pointer to GMI_LOSS_CO!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

         ! Get pointer to trop P(CO) from CH4 if needed
         IF ( LPCO_CH4 ) THEN
            CALL HCO_GetPtr( am_I_Root, HcoState,
     &                      'PCO_CH4', PCO_CH4, RC )
            IF ( RC /= HCO_SUCCESS ) THEN
               ErrMsg = 'Cannot get pointer to PCO_CH4!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
         ENDIF

         IF ( LPCO_NMVOC ) THEN
            CALL HCO_GetPtr( am_I_Root, HcoState,
     &                      'PCO_NMVOC', PCO_NMVOC, RC )
            IF ( RC /= HCO_SUCCESS ) THEN
               ErrMsg = 'Cannot get pointer to PCO_NMVOC!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
         ENDIF

         ! Get pointer to surface CH4 data
         CALL HCO_GetPtr( am_I_Root, HcoState,
     &                    'NOAA_GMD_CH4', SFC_CH4, RC )
         IF ( RC /= HCO_SUCCESS ) THEN
            ErrMsg = 'Cannot get pointer to NOAA_GMD_CH4!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Read emissions from HEMCO into SUMACETCO, SUMISOPCO, SUMMONOCO
      ! arrays.  These are needed to update the chemically-produced
      ! tagged CO species: COch4, COisop, COacet, etc. (bmy, 6/1/16)
      !
      ! NOTE: HEMCO returns 3-D emissions, so we will sum in the
      ! vertical because the code expects SUMACETCO, SUMISOPCO, and
      ! SUMMONOCO to be 2-D arrays.  We may change that later.
      ! (bmy, 6/1/16)
      !
      ! Also note, skip this section if emissions are turned off,
      ! which will keep the SUM*CO arrays zeroed out. (bmy, 10/11/16)
      !=================================================================
      IF ( Input_Opt%LEMIS ) THEN

         ! Conversion factor from [kg/s] --> [atoms C]
         ! (atoms C /mole C) / (kg C /mole C) * chemistry timestep [s]
         kgs_to_atomsC = ( AVO / 12e-3_fp ) * DTCHEM

         ! SUMACETCO (convert [kgC/m2/s] to [atoms C])
         HcoId = GetHcoId( 'ACET' )
         IF ( HcoId > 0 ) THEN
            SUMACETCO = SUM( HcoState%Spc(HcoID)%Emis%Val, 3 )
!kgC/m2/s
            SUMACETCO = SUMACETCO * HcoState%Grid%AREA_M2%Val     !kgC/s
            SUMACETCO = SUMACETCO * kgs_to_atomsC                 !atoms
C
         ELSE
            ErrMsg = 'ACET not turned on in the MEGAN!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

         ! SUMISOPCO (convert [kgC/m2/s] to [atoms C])
         HcoId = GetHcoId( 'ISOP' )
         IF ( HcoId > 0 ) THEN
            SUMISOPCO = SUM( HcoState%Spc(HcoID)%Emis%Val, 3 )
!kgC/m2/s
            SUMISOPCO = SUMISOPCO * HcoState%Grid%AREA_M2%Val     !kgC/s
            SUMISOPCO = SUMISOPCO * kgs_to_atomsC                 !atoms
C
         ELSE
            ErrMsg = 'ISOP not turned on in MEGAN!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

         ! SUMMONOCO (Total monoterpene = MTPA + LIMO + MTPO)
         HcoId = GetHcoId( 'MTPA' )
         IF ( HcoId > 0 ) THEN
            ! kgC/m2/s
            SUMMONOCO = SUM( HcoState%Spc(HcoID)%Emis%Val, 3 )
         ELSE
            ErrMsg = 'MTPA not turned on in Megan_Mono !'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
         HcoId = GetHcoId( 'LIMO' )
         IF ( HcoId > 0 ) THEN
            ! kgC/m2/s
            SUMMONOCO = SUMMONOCO + SUM(HcoState%Spc(HcoID)%Emis%Val,3)
         ELSE
            ErrMsg = 'LIMO not turned on in Megan_Mono !'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
         HcoId = GetHcoId( 'MTPO' )
         IF ( HcoId > 0 ) THEN
            ! kgC/m2/s
            SUMMONOCO = SUMMONOCO + SUM(HcoState%Spc(HcoID)%Emis%Val,3)
         ELSE
            ErrMsg = 'MTPO not turned on in Megan_Mono !'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
         ! SUMMONOCO (convert [kgC/m2/s] to [atoms C])
         SUMMONOCO = SUMMONOCO * HcoState%Grid%AREA_M2%Val ! kgC/s
         SUMMONOCO = SUMMONOCO * kgs_to_atomsC ! atoms C

      ENDIF

      !=================================================================
      ! Do tagged CO chemistry -- Put everything within a large 
      ! DO-loop over all grid boxes to facilitate parallelization
      !=================================================================
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED                                                  ) 
!$OMP+PRIVATE( I,       J,            L,           N,          STTCO   )
!$OMP+PRIVATE( GCO,     DENS,         CH4RATE,     CO_CH4,     CH4     )
!$OMP+PRIVATE( KRATE,   CO_ISOP,      CO_CH3OH,    CO_MONO,    CO_ACET )
!$OMP+PRIVATE( CORATE,  CO_OH,        YMID,        ALPHA_ISOP, KHI1    )
!$OMP+PRIVATE( KLO1,    XYRAT1,       BLOG1,       FEXP1,      KHI2    )
!$OMP+PRIVATE( KLO2,    XYRAT2,       BLOG2,       FEXP2,      KCO1    )
!$OMP+PRIVATE( KCO2,    OH_MOLEC_CM3, FAC_DIURNAL, SUNCOS,     ERR_LOC ) 
!$OMP+PRIVATE( PCO_CH4_MCM3S,         PCO_NMVOC_MCM3S,         CO_NMVOC)
!$OMP+PRIVATE( ERR_VAR, ERR_MSG,      BOXVL, E_CO2                     )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Latitude of grid box
         YMID     = GET_YMID( I, J, L )

         ! Grid box volume [cm3]
         BOXVL    = State_Met%AIRVOL(I,J,L) * 1e+6_fp

         !==============================================================
         ! (0) Define useful quantities
         !==============================================================

         ! STTCO [molec CO/cm3/kg CO] converts [kg CO] --> [molec
         ! CO/cm3]
         ! kg CO/box * box/cm3 * mole/0.028 kg CO * Avog.#/mole
         STTCO = 1.0_fp  / AIRVOL(I,J,L) / 1e+6_fp /
     &           FMOL_CO * AVO

         ! GCO is CO concentration in [molec CO/cm3]
         GCO   = Spc(I,J,L,16) * STTCO

         ! DENS is the number density of air [molec air/cm3]
         DENS  = AD(I,J,L) * 1000.e+0_fp   /
     &           BOXVL     * AVO  / AIRMW

         ! Cosine of the solar zenith angle [unitless]
         SUNCOS = State_Met%SUNCOSmid(I,J)

         ! Scaling factor for diurnal cycles - zero at night
         IF ( SUNCOS > 0.0_fp .and. TCOSZ(I,J) > 0.0_fp ) THEN
            FAC_DIURNAL = ( SUNCOS    / TCOSZ(I,J)    ) *
     &                    ( 86400.0_fp / GET_TS_CHEM() )
         ELSE
            FAC_DIURNAL = 0.0_fp
         ENDIF

         ! Now impose a diurnal cycle on OH.
         ! This is done in other offline simulations but was
         ! missing from tagged CO (jaf, 3/12/14)
         !
         ! NOTE: HEMCO brings in OH in kg/m3, so we need to also
         ! apply a conversion to molec/cm3 here. (bmy, 10/12/16)
         OH_MOLEC_CM3 = ( OH(I,J,L) * kgm3_to_mcm3OH ) * FAC_DIURNAL

         ! Make sure OH is not negative
         OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0e+0_fp )

         ! Also impose diurnal cycle on P(CO) from CH4, NMVOC
         IF ( LPCO_CH4 ) THEN

            ! OLD: HEMCO brings in PCO in kg/m3, so we need to also
            ! apply a conversion to molec/cm3/s 
            ! NEW: Now use the online claculated CH4 from OH
            ! Units in kg CH4, so convert first to kg CO
            ! use AIRVOL to get kg/m3
            ! Remove the FAC_DIURNAL since it is now applied
            ! to the OH used for CH4 directly , (bb 2019)
            PCO_CH4_MCM3S = CH4_OH_TROP(I,J,L) 
     &                      * XNUMOL_CH4 / XNUMOL_CO
     &                      / AIRVOL(I,J,L)
     &                      * kgm3_to_mcm3sCO

            ! Make sure PCO_CH4 is not negative
            PCO_CH4_MCM3S = MAX( PCO_CH4_MCM3S, 0e+0_fp )


         ENDIF

         IF ( LPCO_NMVOC ) THEN

            ! HEMCO brings in PCO in kg/m3/s, so we need to also
            ! apply a conversion to molec/cm3/s 
            PCO_NMVOC_MCM3S = PCO_NMVOC(I,J,L) * kgm3_to_mcm3sCO *
     &                            FAC_DIURNAL

            ! Make sure PCO_NMVOC is not negative
            PCO_NMVOC_MCM3S = MAX( PCO_NMVOC_MCM3S, 0e+0_fp )

         ENDIF

         !==============================================================
         ! (1a) Production of CO by reaction with CH4
         !==============================================================

         ! Initialize
         CO_CH4 = 0e+0_fp

         ! Test level for stratosphere or troposphere
         IF ( State_Met%InStratosphere(I,J,L) ) THEN

            !===========================================================
            ! (1a-1) Production of CO from CH4 in the stratosphere 

            ! Call GET_PCO_LCO_STRAT to get the P(CO) rate from CH4
            CH4RATE = GMI_PROD_CO(I,J,L)
           
            ! Convert units of CH4RATE from [v/v/s] to [molec CO/cm3]
            CO_CH4  = CH4RATE * DTCHEM * DENS

         ELSE

            !===========================================================
            ! (1a-2) Production of CO from CH4 in the troposphere 
            !===========================================================

            ! CH4 concentration from HEMCO [ppbv]
            CH4 = SFC_CH4(I,J)

            ! Convert CH4 from [ppbv] to [molec CH4/cm3]
            CH4 = CH4 * 1e-9_fp * DENS

            ! Use rates saved from full chemistry
            IF ( LPCO_CH4 ) THEN

               ! Diurnal cycle & unit conversion applied above
               ! Multiply by DTCHEM to convert to [molec CO/cm3]
               CO_CH4 = PCO_CH4_MCM3S * DTCHEM

            ! Original behaviour based on latitude bands
            ELSE

               ! CH4 concentration [ppbv] for the given latitude band 
                ! (bmy, 1/2/01)
               CH4 = A3090S
               IF ( YMID >= -30.0 .and. YMID < 0.0  ) CH4 = A0030S
               IF ( YMID >=   0.0 .and. YMID < 30.0 ) CH4 = A0030N
               IF ( YMID >=  30.0                   ) CH4 = A3090N

               ! Convert CH4 from [ppbv] to [molec CH4/cm3]
               CH4 = CH4 * 1e-9_fp * DENS

               ! Calculate updated rate constant [s-1] (bnd, bmy,
               ! 1/2/01)
               KRATE = 2.45e-12_fp * EXP( -1775.e+0_fp / T(I,J,L) )

               ! Production of CO from CH4 = alpha * k * [CH4] * [OH] *
               ! dt 
               ! Units are [molec CO/cm3]
               CO_CH4 = ALPHA_CH4 * KRATE * CH4 * OH_MOLEC_CM3 * DTCHEM

            ENDIF

         ENDIF

         ! Check CO_CH4 for NaN or Infinity
         ERR_LOC = (/ I, J, L, 0 /)
         ERR_VAR = 'CO_CH4'
         ERR_MSG = 'STOP at tagged_co_mod:1'
         CALL CHECK_VALUE( CO_CH4, ERR_LOC, ERR_VAR, ERR_MSG )

         ! Use rates saved from full chemistry
         IF ( LPCO_NMVOC) THEN
            !===========================================================
            ! (1b) Production of CO from NMVOCs (all but CH4)
            !===========================================================

            ! Initialize
            CO_NMVOC = 0e+0_fp

            ! CO production only happens in the troposphere. However,
            ! this is already taken into account in the input files,
            ! which are filled with zeros above the flexchem levels.

            ! Diurnal cycle & unit conversion applied above
            ! Multiply by DTCHEM to convert to [molec CO/cm3]
            CO_NMVOC = PCO_NMVOC_MCM3S * DTCHEM

            ! Make sure it is not negative
            CO_NMVOC = MAX( CO_NMVOC, 0e+0_fp )

            ! Check CO_NMVOC for NaN or Infinity
            ERR_LOC = (/ I, J, L, 0 /)
            ERR_VAR = 'CO_NMVOC'
            ERR_MSG = 'STOP at tagged_co_mod:1'
            CALL CHECK_VALUE( CO_NMVOC, ERR_LOC, ERR_VAR, ERR_MSG )

           ! Set individual VOC contributions to zero to avoid
           ! diagnostic problems
           CO_ISOP  = 0e+0_fp
           CO_CH3OH = 0e+0_fp
           CO_ACET  = 0e+0_fp
           CO_MONO  = 0e+0_fp

         ! Otherwise use original behaviour
         ELSE
            !===========================================================
            ! (1b) Production of CO from ISOPRENE and METHANOL (CH3OH)
            !===========================================================

            ! Initialize
            CO_ISOP  = 0e+0_fp
            CO_CH3OH = 0e+0_fp

            ! Isoprene is emitted only into the surface layer 
            IF ( L == 1 ) THEN

               !========================================================
               ! Yield of CO from ISOP: 30%, from Miyoshi et al., 1994.
               ! They estimate globally 105 Tg C/yr of CO is produced 
               ! from isoprene oxidation. 
               !--------------------------------------------------------
               ! We need to scale the Isoprene flux to get the CH3OH 
               ! (methanol) flux.  Currently, the annual isoprene flux
               ! in 
               ! GEOS-CHEM is ~ 397 Tg C.
               !
               ! Daniel Jacob recommends a flux of 100 Tg/yr CO from
               ! CH3OH 
               ! oxidation based on Singh et al. 2000 [JGR 105,
               ! 3795-3805] 
               ! who estimate a global methanol source of 122 Tg yr-1,
               ! of 
               ! which most (75 Tg yr-1) is "primary biogenic".  He also 
               ! recommends for now that the CO flux from CH3OH
               ! oxidation 
               ! be scaled to monthly mean isoprene flux.
               !
               ! To get CO from METHANOL oxidation, we must therefore 
               ! multiply the ISOPRENE flux by the following scale
               ! factor:
               !  ( 100 Tg CO / 397 Tg C ) * ( 12 g C/mole / 28 g
               !  CO/mole )
               !-----------------------------------------------------------
               ! We now call GET_ALPHA_ISOP to get the yield factor of
               ! CO produced from isoprene, as a function of NOx, or
               ! as a constant. (bnd, bmy, 6/14/01)
               !=======================================================

               !--------------------------------------------------------
               ! Prior to 10/11/18:
               ! Disable this routine for now, to avoid a call
               ! to ERROR_STOP within a parallel loop.  We currently
               ! set ALPHA_ISOP to a constant value anyway.
               ! (bmy, 10/11/18)
               !! Get CO yield from isoprene
               !ALPHA_ISOP = GET_ALPHA_ISOP( .FALSE. )
               !--------------------------------------------------------

               ! Get CO yield from ISOPRENE
               ! Use a 30% yield from Miyoshi et al., 1994.
               ! They estimate globally 105 Tg C/yr of CO is produced
               ! from isoprene oxidation.
               ! ALPHA_ISOP = (0.3 molec CO/atoms C) x (5 atoms C/molec
               ! ISOP)
               ALPHA_ISOP = 1.5e+0_fp

               ! P(CO) from Isoprene Flux = ALPHA_ISOP * Flux(ISOP)
               ! Convert from [molec ISOP/box] to [molec CO/cm3]
               !
               ! Units of SUMISOPCO are [atoms C/box/time step]. 
               ! Division by 5 is necessary to convert to 
               ! [molec ISOP/box/timestep].
               !
               ! Units of ALPHA_ISOP are [molec CO/molec ISOP]
               ! Units of CO_ISOP are [molec CO/cm3]
               CO_ISOP  = SUMISOPCO(I,J)   / BOXVL
     &                  / 5.0_fp           * ALPHA_ISOP

               ! P(CO) from CH3OH is scaled to Isoprene Flux (see above)
               ! Units are [molec CO/cm3]
               CO_CH3OH = ( SUMISOPCO(I,J) / BOXVL    )
     &                  * ( 100.0_fp       / 397.0_fp )
     &                  * ( 12.0_fp        / 28.0_fp  )

               ! Zero SUMISOPCO and SUMCH3OHCO for the next emission
               ! step
               SUMISOPCO(I,J)  = 0e+0_fp
               SUMCH3OHCO(I,J) = 0e+0_fp

               ! Check CO_ISOP for NaN or Infinity
               ERR_LOC = (/ I, J, L, 0 /)
               ERR_VAR = 'CO_ISOP'
               ERR_MSG = 'STOP at tagged_co_mod:2'
               CALL CHECK_VALUE( CO_ISOP,  ERR_LOC, ERR_VAR, ERR_MSG )

               ! Check CO_CH3OH for NaN or Infinity
               ERR_VAR = 'CO_CH3OH'
               ERR_MSG = 'STOP at tagged_co_mod:3'
               CALL CHECK_VALUE( CO_CH3OH, ERR_LOC, ERR_VAR, ERR_MSG )
            ENDIF

            !===========================================================
            ! (1c) Production of CO from MONOTERPENE oxidation
            !===========================================================

            ! Initialize
            CO_MONO = 0.e+0_fp

            ! Monoterpenes are emitted only into the surface layer
            IF ( L == 1 ) THEN

               !=======================================================
               ! Assume the production of CO from monoterpenes is 
               ! instantaneous even though the lifetime of intermediate 
               ! species may be on the order of hours or days.  This 
               ! assumption will likely cause CO from monoterpene 
               ! oxidation to be too high in the box in which the 
               ! monoterpene is emitted.
               !-------------------------------------------------------
               ! The CO yield here is taken from:
               !   Hatakeyama et al. JGR, Vol. 96, p. 947-958 (1991)
               !   Vinckier et al. Fresenius Env. Bull., Vol. 7,
               !   p.361-368 
               !     (1998)
               !
               ! Hatakeyama:  "The ultimate yield of CO from the 
               !   tropospheric oxidation of terpenes (including both O3 
               !   and OH reactions) was estimated to be 20% on the
               !   carbon 
               !   number basis."  They studied ALPHA- & BETA-pinene.
               !
               ! Vinckier  :  "R(CO)=1.8+/-0.3" : 1.8/10 is about 20%.
               !--------------------------------------------------------
               ! Calculate source of CO per time step from monoterpene 
               ! flux (assume lifetime very short) using the C number
               ! basis:
               !
               !   CO [molec CO/cm3] = Flux [atoms C from MONO/box] /
               !                       Grid Box Volume [cm^-3]       *
               !                       ALPHA_MONO 
               !
               ! where ALPHA_MONO = 0.2 as explained above.
               !========================================================

               ! P(CO) from Monoterpene Flux =  alpha * Flux(Mono)
               ! Units are [molec CO/cm3]
               CO_MONO = ( SUMMONOCO(I,J) / BOXVL ) * ALPHA_MONO

               ! Zero SUMMONOCO for the next emission step
               SUMMONOCO(I,J) = 0e+0_fp

               ! Check CO_MONO for NaN or Infinity
               ERR_LOC = (/ I, J, L, 0 /)
               ERR_VAR = 'CO_MONO'
               ERR_MSG = 'STOP at tagged_co_mod:4'
               CALL CHECK_VALUE( CO_MONO, ERR_LOC, ERR_VAR, ERR_MSG )
            ENDIF

            !===========================================================
            ! (1d) Production of CO from oxidation of ACETONE 
            !
            ! ALPHA_ACET = 2/3 to get a yield for CO.  This accounts 
            ! for acetone loss from reaction with OH And photolysis.  
            ! The acetone sources taken into account are:
            ! 
            ! (a) Primary emissions of acetone from biogenic sources
            ! (b) Secondary production of acetone from monoterpene 
            !      oxidation
            ! (c) Secondary production of acetone from ALK4 and 
            !      propane oxidation
            ! (d) Direct emissions of acetone from biomass burning and 
            !      fossil fuels
            ! (e) direct emissions from ocean
            !
            ! Calculate source of CO per time step from biogenic acetone 
            ! # molec CO/cc = ALPHA * ACET Emission Rate * dt
            !===========================================================

            ! Initialize
            CO_ACET = 0.e+0_fp

            ! Biogenic acetone sources are emitted only into the surface
            ! layer
            IF ( L == 1 ) THEN

               ! Units are [molec CO/cc]
               CO_ACET = SUMACETCO(I,J) / BOXVL * ALPHA_ACET

               ! Zero SUMACETCO for the next emission step
               SUMACETCO(I,J) = 0e+0_fp

               ! Check CO_ACET for NaN or Infinity
               ERR_LOC = (/ I, J, L, 0 /)
               ERR_VAR = 'CO_ACET'
               ERR_MSG = 'STOP at tagged_co_mod:5'
               CALL CHECK_VALUE( CO_ACET, ERR_LOC, ERR_VAR, ERR_MSG )
            ENDIF

            ! Add individual NMVOC contributions together to get total
            ! NMVOC contribution
            CO_NMVOC = CO_ISOP + CO_CH3OH + CO_MONO + CO_ACET

         ENDIF !Saved rates vs surface fluxes

         !==============================================================
         ! (1e) Add production of CO into the following tagged species:
         !
         ! (a) Species #12: CO produced from CH4
         ! (a) Species #13: CO produced from NMVOC
         ! (b) Species #14: CO produced from ISOPRENE - old only
         ! (c) Species #15: CO produced from MONOTERPENES - old only
         ! (d) Species #16: CO produced from METHANOL (CH3OH) - old only
         ! (e) Species #17: CO produced from ACETONE - old only
         !
         ! %%% NOTE: If you are modifying the tagged CO simulation, 
         ! %%% and your simulation has less than 12 species, then 
         ! %%% then comment out this section.  If you don't you can
         ! %%% get an array-out-of-bounds error (bmy, 6/11/08)
         !==============================================================
         IF ( LSPLIT ) THEN
            IF ( IDch4 >= 0 )
     &         Spc(I,J,L,IDch4) = Spc(I,J,L,IDch4) + CO_CH4   / STTCO
            IF ( IDnmvoc >= 0 )
     &         Spc(I,J,L,IDnmvoc) = Spc(I,J,L,IDnmvoc) + CO_NMVOC /STTCO
            IF (.not. LPCO_NMVOC) THEN
               IF ( IDisop >= 0 )
     &            Spc(I,J,L,IDisop) = Spc(I,J,L,IDisop) + CO_ISOP /STTCO
               IF ( IDmono >= 0 )
     &            Spc(I,J,L,IDmono) = Spc(I,J,L,IDmono) + CO_MONO /STTCO
               IF ( IDch3oh >= 0 )
     &            Spc(I,J,L,IDch3oh) = Spc(I,J,L,IDch3oh)+CO_CH3OH/STTCO
               IF ( IDacet >= 0 )
     &            Spc(I,J,L,IDacet) = Spc(I,J,L,IDacet) + CO_ACET /STTCO
            ENDIF
         ENDIF

         !==============================================================
         ! (2a) Loss of CO due to chemical reaction w/ OH
         !==============================================================
         ! Select out tropospheric or stratospheric boxes
         IF ( State_Met%InStratosphere(I,J,L) ) THEN

            !===========================================================
            ! (2a-1) Stratospheric loss of CO due to chemical rxn w/ OH
            !===========================================================

            ! Get the L(CO) rate in the stratosphere in [s-1]
            CORATE  = GMI_LOSS_CO(I,J,L)

            ! CO_OH is the fraction of CO lost to OH [unitless]
            CO_OH   = CORATE * DTCHEM

            ! Check CO_OH for NaN or Infinity
            ERR_LOC = (/ I, J, L, 0 /)

            ERR_VAR = 'CO_OH'
            ERR_MSG = 'STOP at tagged_co_mod:6'
            CALL CHECK_VALUE( CO_OH, ERR_LOC, ERR_VAR, ERR_MSG )

            ! Handle strat loss by OH for regional CO species
            IF ( LSPLIT ) THEN

               ! Loop over regional CO species, but exclude CO2
               DO NA = 16, nAdvect-11

                  ! Advected species ID
                  N            = State_Chm%Map_Advect(NA)

                  !------------------------------------------------------
                  ! NOTE: The proper order should be:
                  !   (1) Calculate CO loss rate
                  !   (2) Update AD65 array
                  !   (3) Update the SPC array using the loss rate
                  ! 
                  ! Therefore, we have now moved the computation of the
                  ! ND65 diagnostic before we apply the loss to the 
                  ! tagged CO concentrations stored in the SPC array.
                  !
                  !    -- Jenny Fisher (27 Mar 2017)
                  !
                  !------------------------------------------------------

#if defined( BPCH_DIAG )
                  !-----------------------------------------------------
                  ! ND65 (bpch) diagnostic
                  !
                  ! Loss of CO by OH for "tagged" species
                  !-----------------------------------------------------

                  ! Units: [molec CO/cm3/s]
                  IF ( ND65 > 0 .and. L <= LD65 ) THEN
                     AD65(I,J,L,N) = AD65(I,J,L,N) +
     &                               ( CORATE * Spc(I,J,L,N) * STTCO)
                  ENDIF
#endif

#if defined( NC_DIAG )
                  !-----------------------------------------------------
                  ! HISTORY (aka netCDF diagnostics)
                  !
                  ! Loss of CO by OH for "tagged" species
                  !-----------------------------------------------------
                  !Units: [kg/s]
                  IF ( State_Diag%Archive_Loss ) THEN
                     State_Diag%Loss(I,J,L,N) = ( CORATE
     &                                        *   Spc(I,J,L,N) )
                  ENDIF
#endif

                  ! Loss
                  Spc(I,J,L,N) = Spc(I,J,L,N) * ( 1.0_fp - CO_OH )

                  ! Species shouldn't be less than zero
                  IF ( Spc(I,J,L,N) < 0.0_fp ) THEN
                     Spc(I,J,L,N) = 0.0_fp
                  ENDIF

                  ! Error check
                  ERR_LOC = (/ I, J, L, N /)
                  ERR_VAR = 'Spc (points to State_Chm%Species)'
                  ERR_MSG = 'STOP at tagged_co_mod:7'
                  CALL CHECK_VALUE( Spc(I,J,L,N), ERR_LOC,
     &                              ERR_VAR,      ERR_MSG )

               ENDDO
            ENDIF

            ! CO_OH above is just the fraction of CO lost by OH.  Here
            ! we multiply it by GCO (the initial value of Spc in
            ! molec/cm3)
            ! to convert it to an amount of CO lost by OH [molec/cm3]
            ! (bmy, 2/19/02)
            CO_OH = GCO * CO_OH

         ELSE

            !===========================================================
            ! (2a-2) Tropospheric loss of CO due to chemical rxn w/ OH
            !
            !  DECAY RATE
            !  The decay rate (KRATE) is calculated by:
            !
            !     No change from JPL '97 to JPL '03 (jaf, 2/27/09)
            !     OH + CO -> products (JPL '03)
            !     k = (1 + 0.6Patm) * 1.5E-13
            !
            !  KRATE has units of [ molec^2 CO / cm6 / s ]^-1, 
            !  since this is a 2-body reaction.
            !
            ! Updated rate constant from JPL 2006; now more complicated.
            ! From JPL 2006: "The  reaction between HO and CO to yield
            ! H + CO2 akes place on a potential energy surface that 
            ! contains the radical HOCO.  The yield of H and CO2 is 
            ! diminished as the pressure rises.  The loss of reactants 
            ! is thus the sum of two processes, an association to yield 
            ! HOCO and the chemical activation process yielding H and
            ! CO2." So we now need two complicated reactions. The code
            ! is more or less copied from calcrate.f as implemented by
            ! jmao (jaf, 3/4/09)
            !
            ! Update rate constant for JPL 15-10 (mps, 4/24/17):
            ! GY( A0 = 5.9e-33, B0 = 1.4e0,  A1 = 1.1e-12, B1 = -1.3e0,
            !     A2 = 1.5e-13, B2 = -0.6e0, A3 = 2.1e09,  B3 = -6.1e0 )
            !
            ! is now
            !
            ! GY( A0 = 5.9e-33, B0 = 1.,     A1 = 1.1e-12, B1 = -1.3e0,
            !     A2 = 1.5e-13, B2 = 0.,     A3 = 2.1e09,  B3 = -6.1e0 )
            !===========================================================

            ! Decay rate
            ! NOTE: This code is nearly identical to function GC_OHCO
            !  found in KPP/Standard/gckpp_Rates.F90)
            ! new JPL 2006 version (jaf, 3/4/09)
            ! KLO1 = k_0(T) from JPL Data Eval (page 2-1)
            KLO1   = 5.9e-33_fp * ( 300 / T(I,J,L) )**(1.e+0_fp)
            ! KHI1 = k_inf(T) from JPL Data Eval (page 2-1)
            KHI1   = 1.1e-12_fp * ( 300 / T(I,J,L) )**(-1.3e+0_fp)
            XYRAT1 = KLO1 * DENS / KHI1
            BLOG1  = LOG10(XYRAT1)
            FEXP1  = 1.e+0_fp / ( 1.e+0_fp + BLOG1 * BLOG1 )
            ! KCO1 = k_f([M],T) from JPL Data Eval (page 2-1)
            KCO1   = KLO1 * DENS * 0.6**FEXP1 / ( 1.e+0_fp + XYRAT1 )
            ! KLO2 = k_0(T) from JPL Data Eval (page 2-1)
            KLO2   = 1.5e-13_fp * ( 300 / T(I,J,L) )**(0.e+0_fp)
            ! KHI2 = k_inf(T) from JPL Data Eval (page 2-1)
            KHI2   = 2.1e+09_fp * ( 300 / T(I,J,L) )**(-6.1e+0_fp)
            XYRAT2 = KLO2 * DENS / KHI2
            BLOG2  = LOG10(XYRAT2)
            FEXP2  = 1.e+0_fp / ( 1.e+0_fp + BLOG2 * BLOG2 )
            ! KCO2 = k_f^ca([M],T) from JPL Data Eval (page 2-2)
            KCO2   = KLO2 * 0.6**FEXP2 / ( 1.e+0_fp + XYRAT2 )

            ! KRATE is the sum of the two.
            KRATE  = KCO1 + KCO2

            ! CO_OH = Tropospheric loss of CO by OH [molec/cm3]
            ! Now use OH_MOLEC_CM3, which includes a diurnal cycle.
            CO_OH = KRATE * GCO * OH_MOLEC_CM3 * DTCHEM

            ! Handle trop loss by OH for regional CO species
            IF ( LSPLIT ) THEN

               ! Loop over regional CO species
               DO NA = 16, nAdvect-11

                  ! Advected species ID
                  N            = State_Chm%Map_Advect(NA)

                  !-----------------------------------------------------
                  ! NOTE: The proper order should be:
                  !   (1) Calculate CO loss rate
                  !   (2) Update AD65 array
                  !   (3) Update the SPC array using the loss rate
                  ! 
                  ! Therefore, we have now moved the computation of the
                  ! ND65 diagnostic before we apply the loss to the 
                  ! tagged CO concentrations stored in the SPC array.
                  !
                  !    -- Jenny Fisher (27 Mar 2017)
                  !-----------------------------------------------------

#if defined( BPCH_DIAG )
                  !-----------------------------------------------------
                  ! ND65 (bpch) diagnostic
                  !
                  ! Loss of CO by OH for "tagged" species
                  !-----------------------------------------------------

                  ! Units: [molec CO/cm3/s]
                  IF ( ND65 > 0 .and. L <= LD65 ) THEN
                     AD65(I,J,L,N) = AD65(I,J,L,N) +
     &                   ( KRATE * OH_MOLEC_CM3 * Spc(I,J,L,N) * STTCO )
                  ENDIF
#endif

#if defined( NC_DIAG )
                  !-----------------------------------------------------
                  ! HISTORY (aka netCDF diagnostics)
                  !
                  ! Loss of CO by OH for "tagged" species
                  !-----------------------------------------------------

                  ! Units: [kg/s]
                  IF ( State_Diag%Archive_Loss ) THEN
                     State_Diag%Loss(I,J,L,N) = ( KRATE
     &                                        *   OH_MOLEC_CM3
     &                                        *   Spc(I,J,L,N) )
                  ENDIF
#endif

                  ! Use tropospheric rate constant 
                  Spc(I,J,L,N) = Spc(I,J,L,N) *
     &                 ( 1e+0_fp - KRATE * OH_MOLEC_CM3 * DTCHEM )

                  ! Error check
                  ERR_LOC = (/ I, J, L, N /)
                  ERR_VAR = 'Spc (points to State_Chm%Species)'
                  ERR_MSG = 'STOP at tagged_co_mod:8'
                  CALL CHECK_VALUE( Spc(I,J,L,N), ERR_LOC,
     &                              ERR_VAR,      ERR_MSG )
               ENDDO
            ENDIF

         ENDIF

         !==============================================================
         ! Save the total chemical production from various sources
         ! into the total CO species Spc(I,J,L,1)
         !==============================================================

         ! GCO is the total CO before chemistry was applied [molec
         ! CO/cm3]
         ! Add to GCO the sources and sinks listed above
         GCO = GCO + CO_CH4 + CO_NMVOC - CO_OH
         ! Convert net CO from [molec CO/cm3] to [kg] and store in
         ! State_Chm%Species
         Spc(I,J,L,16) = GCO / STTCO

         !==============================================================
         ! Calculate the chemical production of CO2
         ! based on the CO loss by OH
         !==============================================================

         IF ( Input_Opt%LCHEMCO2 ) THEN

            ! Production is in [molec CO/cm3], convert to [molec/cm2/s]
            ! 1 molec CO = 1 molec CO2 (bb 2019)
            E_CO2 = CO_OH    !molecCO/cm3
     &            / DTCHEM   !=>molec/cm3/s
     &            *State_Met%BXHEIGHT(I,J,L) * 100   !=>molec/cm2/s

#if defined( BPCH_DIAG )
            !==========================================================
            ! %%%%% ND04 (bpch) DIAGNOSTICS %%%%%
            ! 
            ! Save production of CO2 from CO oxidation [molec/cm2/s]
            !==========================================================
            IF ( ND04 > 0 ) THEN
               AD04_chem(I,J,L) = AD04_chem(I,J,L) + E_CO2
            ENDIF
#endif

#if defined( NC_DIAG )
            !==========================================================
            ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
            ! 
            ! Save production of CO2 from CO oxidation [kg/m2/s]
            !==========================================================
            IF ( State_Diag%Archive_ProdCO2fromCO ) THEN
               State_Diag%ProdCO2fromCO(I,J,L) = E_CO2      !molec/cm2/s 
     &                                         / XNUMOL_CO2  !=>kg/cm2/s
     &                                         * 1e4_fp      ! =>kg/m2/s

            ENDIF
#endif

            ! Convert emissions from [molec/cm2/s] to [kg/kg dry air]
            ! (ewl, 9/11/15)
            E_CO2  =  E_CO2 * DTCHEM * CM2PERM2 /
     &                ( XNUMOL_CO2 * State_Met%DELP(I,J,L)
     &                * G0_100 * ( 1.0e+0_fp
     &                - State_Met%SPHU(I,J,L) * 1.0e-3_fp ) )


            ! Spc array has been converted to kg for all species, inc.                                                        
            ! CO2 when joint simulation is used. So convert E_CO2 from                                                        
            ! kg/kg dry air to kg (jaf, bb 2019)
            E_CO2 = E_CO2 * AD(I,J,L)

            ! Add to Species #1: Total CO2 [kg/kg]
            Spc(I,J,L,29) = Spc(I,J,L,29) + E_CO2

            ! Add to Species #10: Chemical Source of CO2 [kg/kg]
            !IF ( nAdvect > 22 ) THEN
            Spc(I,J,L,38) = Spc(I,J,L,38) + E_CO2
            !ENDIF
         ENDIF
#if defined( BPCH_DIAG )
         !==============================================================
         ! ND65 (bpch) diagnostics -- Production & Loss of CO
         ! Also save P(CO) from CH3OH and MONOTERPENES (bmy, 1/2/01)
         !==============================================================
         IF ( ND65 > 0 .and. L <= LD65 ) THEN

            ! Loss of CO by OH (global) [molec CO/cm3/s] 
            N             = 1
            AD65(I,J,L,N) = AD65(I,J,L,N) + ( CO_OH / DTCHEM )

            ! Production of CO from CH4 [molec CO/cm3/s]
            N             = nAdvect + 1
            AD65(I,J,L,N) = AD65(I,J,L,N) + ( CO_CH4 / DTCHEM )

            ! Production of CO from NMVOCs [molec CO/cm3/s]
            N             = nAdvect + 2
            AD65(I,J,L,N) = AD65(I,J,L,N) + ( CO_NMVOC / DTCHEM )

            IF ( .not. LPCO_NMVOC ) THEN
               ! Production of CO from Isoprene [molec CO/cm3/s]
               N             = nAdvect + 3
               AD65(I,J,L,N) = AD65(I,J,L,N) + ( CO_ISOP / DTCHEM )

               ! Production of CO from CH3OH [molec CO/cm3/s]            
               N             = nAdvect + 4
               AD65(I,J,L,N) = AD65(I,J,L,N) + ( CO_CH3OH / DTCHEM )

               ! Production of CO from MONO [molec CO/cm3/s]
               N             = nAdvect + 5
               AD65(I,J,L,N) = AD65(I,J,L,N) + ( CO_MONO / DTCHEM )

               ! Production of CO from ACET [molec CO/cm3/s]
               N             = nAdvect + 6
               AD65(I,J,L,N) = AD65(I,J,L,N) + ( CO_ACET / DTCHEM )
            ENDIF
         ENDIF
#endif

#if defined( NC_DIAG )
         !==============================================================
         ! HISTORY (aka netCDF diagnostics)
         !
         ! Production of CO species 
         !==============================================================

         ! Units: [kg/s] Production of CO from CH4
         IF ( State_Diag%Archive_ProdCOfromCH4 ) THEN
            State_Diag%ProdCOfromCH4(I,J,L) = CO_CH4
     &                                     / STTCO / DTCHEM
         ENDIF

         ! Units: [kg/s] Production of CO from NMVOCs
         IF ( State_Diag%Archive_ProdCOfromNMVOC ) THEN
            State_Diag%ProdCOfromNMVOC(I,J,L) = CO_NMVOC
     &                                     / STTCO / DTCHEM
         ENDIF

         !==============================================================
         ! HISTORY (aka netCDF diagnostics)
         !
         ! Loss of total CO species 
         !==============================================================

         ! Units: [kg/s] 
         IF ( State_Diag%Archive_Loss ) THEN
            State_Diag%Loss(I,J,L,16) = ( CO_OH / STTCO / DTCHEM )
         ENDIF
#endif
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointers
      AD       => NULL()
      AIRVOL   => NULL()
      Spc      => NULL()
      T        => NULL()

      ! Set FIRSTCHEM to FALSE
      FIRSTCHEM = .FALSE.

      END SUBROUTINE CHEM_CH4COCO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_decay
!
! !DESCRIPTION: Subroutine CH4\_DECAY calculates the decay rate of CH4 by OH.
!  OH is the only sink for CH4 considered here. (jsw, bnd, bmy, 1/16/01,
!  7/20/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CH4_DECAY_CH4COCO2( am_I_Root, Input_Opt,  State_Met, 
     &                      State_Chm, State_Diag, RC         )
!
! !USES:
!
      USE CMN_SIZE_MOD
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD
      USE DIAG_MOD,       ONLY : AD19   	      
#endif
      USE ErrCode_Mod
      USE Input_Opt_Mod,  ONLY : OptInput
      USE State_Chm_Mod,  ONLY : ChmState
      USE State_Diag_Mod, ONLY : DgnState
      USE State_Met_Mod,  ONLY : MetState
      USE TIME_MOD,       ONLY : GET_TS_CHEM
      USE TIME_MOD,       ONLY : GET_MONTH
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REMARKS:
!  We now use function ITS_IN_THE_CHEMGRID from chemgrid_mod.F to diagnose
!  if box (I,J,L) is in the troposphere or stratosphere.
!                                                                             .
!  Monthly loss of CH4 is summed in TCH4(3)
!     TCH4(3)  = CH4 sink by OH
! 
! !REVISION HISTORY:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_DECAY is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Now use function GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (4 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  12 Feb 2014 - K. Wecht    - Disable CH4 budget diagnostic (bracket the 
!                              code out with #ifdef blocks so it can be used)
!  24 Jul 2014 - R. Yantosca - Now compute BOXVL internally
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  16 Jun 2017 - M. Sulprizio- Add loss of CH4 by Cl
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER           :: I,  J,    L
      REAL(fp)          :: DT, GCH4, Spc2GCH4
      REAL(fp)          :: KRATE, C_OH, FAC_DIURNAL, SUNCOS
      REAL(fp)          :: KRATE_Cl, C_Cl

      ! Pointers
      REAL(fp), POINTER :: Spc(:,:,:,:)

      !=================================================================
      ! CH4_DECAY begins here!
      !=================================================================
      ! Assume success
      RC = GC_SUCCESS

      ! Chemistry timestep in seconds
      DT = GET_TS_CHEM()

      ! Point to the chemical species array
      Spc => State_Chm%Species

#if defined( NC_DIAG )
      !=================================================================
      ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
      !
      ! Zero the relevant diagnostic fields of State_Diag because the 
      ! position of the tropopause changes from one timestep to the next
      !=================================================================
      IF ( State_Diag%Archive_LossCH4byClinTrop ) THEN
         State_Diag%LossCH4byClinTrop = 0.0_f4
      ENDIF

      IF ( State_Diag%Archive_LossCH4byOHinTrop ) THEN
         State_Diag%LossCH4byOHinTrop = 0.0_f4
      ENDIF
#endif

      !=================================================================
      ! Compute decay of CH4 by OH and Cl in the troposphere
      !
      ! The decay for CH4 is calculated by:
      !    OH + CH4 -> CH3 + H2O 
      !    k = 2.45E-12 exp(-1775/T)
      !
      !    This is from JPL '97.
      !    JPL '00, '06, & '11 do not revise '97 value. (jsw, kjw, ajt)
      !
      ! The decay for CH4 by Cl is calculated by:
      !    Cl + CH4 -> HCl + CH3
      !    k = 9.6E-12 exp(-1360/T)
      !
      !    This is from Kirschke et al., Nat. Geosci., 2013.
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( L, J, I, KRATE, Spc2GCH4, GCH4, C_OH )
!$OMP+PRIVATE( C_Cl, KRATE_Cl,FAC_DIURNAL, SUNCOS )
!$OMP+REDUCTION( +:TROPOCH4 )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Only consider tropospheric boxes
         IF ( State_Met%InChemGrid(I,J,L) ) THEN

            ! Calculate rate coefficients
            KRATE    = 2.45e-12_fp * EXP( -1775e+0_fp /
     &                                    State_Met%T(I,J,L))
            KRATE_Cl = 9.60e-12_fp * EXP( -1360e+0_fp /
     &                                    State_Met%T(I,J,L))

            ! Conversion from [kg/box] --> [molec/cm3]
            ! [kg CH4/box] * [box/cm3] * XNUMOL_CH4 [molec CH4/kg CH4]
            Spc2GCH4 = 1e+0_fp / State_Met%AIRVOL(I,J,L) / 1e+6_fp 
     &                 * XNUMOL_CH4 

            ! CH4 in [molec/cm3]
            GCH4 = Spc(I,J,L,1) * Spc2GCH4

            ! Cosine of the solar zenith angle [unitless]
            SUNCOS = State_Met%SUNCOSmid(I,J)

            ! Scaling factor for diurnal cycles - zero at night
            IF ( SUNCOS > 0.0_fp .and. TCOSZ(I,J) > 0.0_fp ) THEN
               FAC_DIURNAL = ( SUNCOS    / TCOSZ(I,J)    ) *
     &                       ( 86400.0_fp / DT )
            ELSE
               FAC_DIURNAL = 0.0_fp
            ENDIF

            ! OH in [molec/cm3]
            ! BOH is imported from HEMCO in units of kg/m3, convert
            !  here to molec/cm3 (ckeller, 9/16/2014)
            C_OH = BOH(I,J,L) * XNUMOL_OH / CM3PERM3 * FAC_DIURNAL

            ! Cl in [molec/cm3]
            ! BCl is imported from HEMCO in units of ppbv, convert
            !  here to molec/cm3 (mps, 6/16/2017)
            C_Cl = BCl(I,J,L) * BAIRDENS(I,J,L) * 1e-9_fp

            TROPOCH4=TROPOCH4 + GCH4 * KRATE    * C_OH * DT / Spc2GCH4
     &                        + GCH4 * KRATE_Cl * C_Cl * DT / Spc2GCH4

#if defined( BPCH_DIAG )
            !-----------------------------------------------------------
            ! %%%%% ND19 (bpch) diagnostic %%%%%
            !
            ! (1) How much CH4 (kg) is lost by reaction with OH
            ! (3) How much CH4 (kg) is lost by reaction with Cl
            !-----------------------------------------------------------
            IF ( ND19 > 0 ) THEN
               IF ( L <= LD19 ) THEN
                  AD19(I,J,L,1) = AD19(I,J,L,1) + 
     &                 ( GCH4 * KRATE    * C_OH * DT ) / Spc2GCH4

                  AD19(I,J,L,3) = AD19(I,J,L,3) + 
     &                 ( GCH4 * KRATE_Cl * C_Cl * DT ) / Spc2GCH4

               ENDIF
            ENDIF
#endif

#if defined( NC_DIAG )
            !-----------------------------------------------------------
            ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
            !
            ! Archive Loss of CH4 (kg/s) reactions with OH and Cl
            !-----------------------------------------------------------       

            ! Loss CH4 by reaction with Cl [kg/s]
            IF ( State_Diag%Archive_LossCH4byClinTrop ) THEN
               State_Diag%LossCH4byClinTrop(I,J,L) = 
     &            ( GCH4 * KRATE_Cl * C_Cl ) / Spc2GCH4
            ENDIF

            IF ( State_Diag%Archive_LossCH4byOHinTrop ) THEN
               State_Diag%LossCH4byOHinTrop(I,J,L) = 
     &            ( GCH4 * KRATE * C_OH ) / Spc2GCH4
            ENDIF
#endif

            ! Calculate new CH4 value: [CH4]=[CH4](1-k[OH]*delt) 
            GCH4 = GCH4 * ( 1e+0_fp - KRATE    * C_OH * DT )
            GCH4 = GCH4 * ( 1e+0_fp - KRATE_Cl * C_Cl * DT )

            ! Convert back from [molec/cm3] --> [kg/box]
            Spc(I,J,L,1) = GCH4 / Spc2GCH4

            !Write the loss into  an array, need for the CO chmistry 
            !Unit kg CH4 (bb 2019)
            CH4_OH_TROP(I,J,L)=( GCH4 * KRATE * C_OH * DT ) / Spc2GCH4

         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

!      print*,'% --- CHEMCH4: CH4_DECAY: TROP DECAY (Tg): ',TROPOCH4/1e9
!      print*,'Trop decay should be over 1Tg per day globally'
!      print*,'    ~ 500Tg/365d ~ 1.37/d'

      ! Free pointers
      Spc => NULL()

      END SUBROUTINE CH4_DECAY_CH4COCO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_ohsave
!
! !DESCRIPTION: Subroutine CH4\_OHSAVE archives the CH3CCl3 lifetime from the
!  OH used in the CH4 simulation. (bnd, jsw, bmy, 1/16/01, 7/20/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CH4_OHSAVE_CH4COCO2( State_Met, State_Chm )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE DIAG_OH_MOD,        ONLY : DO_DIAG_OH_CH4
      USE GC_GRID_MOD,        ONLY : GET_AREA_CM2
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM

!
! !INPUT PARAMETERS:
!
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
! 
! !REMARKS:
!  We now use function ITS_IN_THE_CHEMGRID from chemgrid_mod.F to diagnose
!  if box (I,J,L) is in the troposphere or stratosphere.
! 
! !REVISION HISTORY:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_OHSAVE is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Now call DO_DIAG_OH_CH4 to pass OH diagnostic info to the
!        "diag_oh_mod.f" (bmy, 7/20/04)
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  24 Jul 2014 - R. Yantosca - Now compute BOXVL internally
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  14 Dec 2017 - M. Sulprizio- Parallelize DO loop
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER           :: I,       J,           L
      REAL(fp)          :: MASST,   AREA_M2
      REAL(fp)          :: KCLO,    LOSS,        OHMASS  
      REAL(fp)          :: KCH4,    CH4LOSE,     CH4MASS
      REAL(fp)          :: CH4EMIS, CH4TROPMASS, BOXVL
      REAL(fp)          :: C_OH, DT, FAC_DIURNAL, SUNCOS 

      ! Pointers
      REAL(fp), POINTER :: Spc(:,:,:,:)

      !=================================================================
      ! CH4_OHSAVE begins here!
      !
      ! (1) Pass OH mass, total air mass, and  to "diag_oh_mod.f"
      ! (2) ND59: Diagnostic for CH3CCl3 calculation
      !=================================================================
      ! Chemistry timestep in seconds
      DT = GET_TS_CHEM()

      ! Point to chemical species array [kg]
      Spc => State_Chm%Species

      ! Calculate OH mass and total air mass
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, BOXVL, C_OH, OHMASS, MASST, KCLO, LOSS, KCH4 )
!$OMP+PRIVATE( CH4TROPMASS, CH4MASS, CH4LOSE, CH4EMIS, AREA_M2 )
!$OMP+PRIVATE( FAC_DIURNAL, SUNCOS )
      DO L = 1, LLPAR
      DO J = 1, JJPAR 
      DO I = 1, IIPAR 

         ! Only process tropospheric boxes (bmy, 4/17/00)
         IF ( State_Met%InChemGrid(I,J,L) ) THEN

            ! Grid box volume [cm3]
            BOXVL  = State_Met%AIRVOL(I,J,L) * 1e+6_fp

            ! Cosine of the solar zenith angle [unitless]
            SUNCOS = State_Met%SUNCOSmid(I,J)

            ! Scaling factor for diurnal cycles - zero at night
            IF ( SUNCOS > 0.0_fp .and. TCOSZ(I,J) > 0.0_fp ) THEN
               FAC_DIURNAL = ( SUNCOS    / TCOSZ(I,J)    ) *
     &                       ( 86400.0_fp / DT )
            ELSE
               FAC_DIURNAL = 0.0_fp
            ENDIF

            ! OH concentration in molec/cm3. BOH is imported
            ! from HEMCO in kg/m3 (ckeller, 9/16/2014)
            C_OH = BOH(I,J,L) * XNUMOL_OH / CM3PERM3 * FAC_DIURNAL

            ! Calculate OH mass [molec / box]
            OHMASS = C_OH * BAIRDENS(I,J,L) * BOXVL

            ! Calculate total air mass [molec / box]
            MASST  = BAIRDENS(I,J,L) * BOXVL

            ! Calculate CH3CCl3 + OH rate constant from JPL '06
            ! [cm3 / molec / s]
            KCLO = 1.64e-12_fp * EXP( -1520.e+0_fp / State_Met%T(I,J,L))

            ! Calculate Loss term [molec / box / s]
            LOSS   = KCLO            * C_OH  *
     &               BAIRDENS(I,J,L) * BOXVL


            ! Calculate CH4 + OH rate constant from JPL '06
            ! [cm3 / molec / s]
            KCH4 = 2.45e-12_fp * EXP( -1775e+0_fp / State_Met%T(I,J,L) )

            ! Calculate CH4 mass [molec / box] from [kg / box]
            CH4TROPMASS = Spc(I,J,L,1) * XNUMOL_CH4 
            CH4MASS     = Spc(I,J,L,1) * XNUMOL_CH4 

            ! Calculate loss term  [molec /box / s]
            CH4LOSE = KCH4            * C_OH  *
     &                BAIRDENS(I,J,L) * BOXVL

            ! Calculate CH4 emissions [molec / box / s]
            !   Only for surface level
            ! Grid box surface area [cm2]
            ! HEMCO update: CH4_EMIS now in kg/m2/s (ckeller, 9/12/2014)
            IF ( L .GT. 1 ) THEN 
               CH4EMIS = 0e+0_fp
            ELSE
               AREA_M2 = GET_AREA_CM2( I, J, L ) / 10000.0e+0_fp

               ! [kg/m2/s]  --> [molec/box/s]
               CH4EMIS  = CH4_EMIS_J(I,J,1)
               CH4EMIS  = CH4EMIS * AREA_M2 * XNUMOL_CH4

            ENDIF

         ELSE

            OHMASS      = 0e+0_fp
            MASST       = 0e+0_fp
            LOSS        = 0e+0_fp
            CH4LOSE     = 0e+0_fp
            CH4TROPMASS = 0e+0_fp
            CH4EMIS     = 0e+0_fp
            CH4MASS     = Spc(I,J,L,1) * XNUMOL_CH4 

         ENDIF

         ! Pass OH mass, total mass, and loss to "diag_oh_mod.f",
         ! which calculates mass-weighted mean [OH] and CH3CCl3
         ! lifetime.
         CALL DO_DIAG_OH_CH4( I, J, L, OHMASS, MASST, LOSS, 
     &        CH4LOSE, CH4TROPMASS, CH4EMIS, CH4MASS )

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      Spc => NULL()

      END SUBROUTINE CH4_OHSAVE_CH4COCO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_strat
!
! !DESCRIPTION: Subroutine CH4\_STRAT calculates uses production rates for CH4
!  to  calculate loss of CH4 in above the tropopause. (jsw, bnd, bmy, 1/16/01,
!  7/20/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CH4_STRAT_CH4COCO2( am_I_Root, Input_Opt,  State_Met, 
     &                      State_Chm, State_Diag, RC         )
!
! !USES:
!
      USE CMN_SIZE_MOD
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD
      USE DIAG_MOD,       ONLY : AD19
#endif
      USE ErrCode_Mod
      USE Input_Opt_Mod,  ONLY : OptInput
      USE State_Chm_Mod,  ONLY : ChmState
      USE State_Diag_Mod, ONLY : DgnState
      USE State_Met_Mod,  ONLY : MetState
      USE TIME_MOD,       ONLY : GET_TS_CHEM

! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
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
!  Production (mixing ratio/sec) rate provided by Dylan Jones.  
!  Only production by CH4 + OH is considered.
!                                                                             .
!  We now use function ITS_IN_THE_CHEMGRID from chemgrid_mod.F to diagnose
!  if box (I,J,L) is in the troposphere or stratosphere.
!
! !REVISION HISTORY:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_STRAT is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Removed LMN from the arg list and made it a local variable.  Now use 
!        functions GET_MONTH and GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (4 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  23 Jun 2014 - R. Yantosca - Now accept am_I_Root and RC
!  24 Jul 2014 - R. Yantosca - Now compute BOXVL internally
!  24 Mar 2015 - E. Lundgren - Now pass Input_Opt to Check_STT
!  06 Jan 2016 - E. Lundgren - Use global physical parameters
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  14 Dec 2017 - M. Sulprizio- Parallelize DO loop
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER           :: I,  J,    L
      REAL(fp)          :: DT, GCH4, Spc2GCH4, LRATE,  BOXVL

      ! Local variables for quantities from Input_Opt
      LOGICAL           :: LUCX

      ! Pointers
      REAL(fp), POINTER :: Spc(:,:,:,:)

      !=================================================================
      ! CH4_STRAT begins here!
      !=================================================================
      ! Assume success
      RC                  =  GC_SUCCESS

      ! Copy fields from INPUT_OPT
      LUCX                =  Input_Opt%LUCX

      ! Point to chemical species
      Spc                 => State_Chm%Species

      ! If unified chemistry is active, ignore all of this
      IF ( .not. LUCX ) THEN

#if defined( NC_DIAG )
         !==============================================================
         ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
         !
         ! Zero the relevant diagnostic fields of State_Diag because 
         ! the position of the tropopause changes from one timestep 
         ! to the next
         !==============================================================
         IF ( State_Diag%Archive_LossCH4inStrat ) THEN
            State_Diag%LossCH4inStrat = 0.0_f4
         ENDIF
#endif

         ! Chemistry timestep [s]
         DT  = GET_TS_CHEM()

         !=================================================================
         ! Loop over stratospheric boxes only
         !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, Spc2GCH4, GCH4, LRATE )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Only proceed if we are outside of the chemistry grid
            IF ( .not. State_Met%InChemGrid(I,J,L) ) THEN

               ! Conversion factor [kg/box] --> [molec/cm3]
               ! [kg/box] / [AIRVOL * 1e6 cm3] * [XNUMOL_CH4 molec/mole]
               Spc2GCH4 = 1e+0_fp / State_Met%AIRVOL(I,J,L) / 1e+6_fp *
     &                    XNUMOL_CH4

               ! CH4 in [molec/cm3]
               GCH4 = Spc(I,J,L,1) * Spc2GCH4

               ! Loss rate [molec/cm3/s]
               LRATE = GCH4 * CH4LOSS( I,J,L )

               ! Update Methane concentration in this grid box [molec/cm3]
               GCH4 = GCH4 - ( LRATE * DT )

               ! Convert back from [molec CH4/cm3] --> [kg/box] 
               Spc(I,J,L,1) = GCH4 / Spc2GCH4

#if defined( BPCH_DIAG )
               !------------------------------------------------------------
               ! %%%%% ND19 (bpch) diagnostic %%%%%
               !
               ! Loss of CH4 by OH above tropopause [kg]
               !------------------------------------------------------------
               IF ( ND19 > 0 ) THEN ! --> [kg/box]
                  AD19(I,J,L,2) = AD19(I,J,L,2) +
     &	                          ( LRATE * DT ) / Spc2GCH4
               ENDIF
#endif

#if defined( NC_DIAG )
               !------------------------------------------------------------
               ! %%%%%% HISTORY (aka netCDF diagnostics) %%%%%
               !
               ! Loss of CH4 by OH above tropopause [kg/s]
               !------------------------------------------------------------
               IF ( State_Diag%Archive_LossCH4inStrat ) THEN
                  State_Diag%LossCH4inStrat(I,J,L) = LRATE / Spc2GCH4
               ENDIF
#endif

            ENDIF
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF ! not LUCX (SDE 03/25/13)

      ! Free pointer
      Spc => NULL()

      END SUBROUTINE CH4_STRAT_CH4COCO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_distrib
!
! !DESCRIPTION: Subroutine CH4\_DISTRIB allocates the chemistry sink to
!  different emission species. (ccc, 10/2/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CH4_DISTRIB_CH4COCO2( PREVCH4, Input_Opt, State_Chm )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : SAFE_DIV
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
   
      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      REAL(fp)                      :: PREVCH4(IIPAR,JJPAR,LLPAR)! CH4 bef chem
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
! 
! !REVISION HISTORY:
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  25 Mar 2013 - R. Yantosca - Now accept Input_Opt, State_Chm args
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER           :: I, J, L, N, NA, nAdvect

      ! Pointers
      REAL(fp), POINTER :: Spc(:,:,:,:)

      !========================================================================
      ! CH4_DISTRIB begins here
      !========================================================================
      ! Point to chemical species array [kg]
      Spc => State_Chm%Species

      ! fix nAdvect (Xueying Yu, 12/10/2017)
      nAdvect = State_Chm%nAdvect

      ! Loop over the number of advected species
      DO NA = 2, nAdvect-24

         ! Advected species ID
         N = State_Chm%Map_Advect(NA)

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
           Spc(I,J,L,N) = SAFE_DIV(Spc(I,J,L,N),PREVCH4(I,J,L),0.e+0_fp)
     &                     * Spc(I,J,L,1)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDDO

      ! Free pointer
      Spc => NULL()

      END SUBROUTINE CH4_DISTRIB_CH4COCO2

!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_diurnal
!
! !DESCRIPTION: Subroutine CALC\_DIRUNAL computes the sume of the cosine
!  of the solar zenith angle over a 24 hour day as well as the total
!  length of daylight to scale the offline OH concentrations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_DIURNAL_CH4COCO2
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GC_GRID_MOD, ONLY : GET_YMID_R
      USE TIME_MOD, ONLY : ITS_A_NEW_DAY
      USE TIME_MOD, ONLY : GET_MINUTE,    GET_SECOND,      GET_HOUR
      USE TIME_MOD, ONLY : GET_TS_CHEM,   GET_DAY_OF_YEAR, GET_LOCALTIME
! 
! !REVISION HISTORY: 
!  12 Mar 2014 - J. Fisher   - Copied from OHNO3TIME in carbon_mod and
!                              COSSZA in dao_mod
!  12 Feb 2018 - E. Lundgren - Updates for new timestep units (sec not
!  min)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER             :: I, J, N, NDYSTEP
      INTEGER             :: SECOND,  MINUTE, TS_SUN
      REAL*8              :: GMT_MID, TIMLOC, FACTOR
      REAL*8              :: R,       AHR,    DEC
      REAL*8              :: YMID_R,  SUNTMP_MID
!
! !DEFINED PARAMETERS:
!
      ! Coefficients for solar declination angle
      REAL*8,  PARAMETER :: A0 = 0.006918d0
      REAL*8,  PARAMETER :: A1 = 0.399912d0
      REAL*8,  PARAMETER :: A2 = 0.006758d0
      REAL*8,  PARAMETER :: A3 = 0.002697d0
      REAL*8,  PARAMETER :: B1 = 0.070257d0
      REAL*8,  PARAMETER :: B2 = 0.000907d0
      REAL*8,  PARAMETER :: B3 = 0.000148d0

      !=================================================================
      ! CALC_DIURNAL begins here!
      !=================================================================
      ! Only do at the start of a new day
      IF ( FIRST .or. ITS_A_NEW_DAY() ) THEN

         ! Zero array
         TCOSZ = 0d0

         ! Get time for central chemistry timestep
         TS_SUN = GET_TS_CHEM()                     ! Chemistry interval
         SECOND = GET_SECOND()                      ! Current seconds
         MINUTE = GET_MINUTE()                      ! Current minutes
         FACTOR = ( MINUTE * 60 + SECOND ) / TS_SUN ! Multiplying factor

         ! GMT at the midpoint of the chemistry time interval for first
         ! timestep of the day
         GMT_MID  = ( DBLE( GET_HOUR()        )      )
     &            + ( DBLE( TS_SUN * FACTOR ) / 3600d0 )
     &            + ( DBLE( TS_SUN / 2      ) / 3600d0 )

         ! Solar declination angle (low precision formula):
         ! Path length of earth's orbit traversed since Jan 1 [radians]
         R = ( 2d0 * PI / 365d0 ) * FLOAT( GET_DAY_OF_YEAR() - 1 )
         DEC = A0 - A1*COS(    R) + B1*SIN(    R)
     &            - A2*COS(2d0*R) + B2*SIN(2d0*R)
     &            - A3*COS(3d0*R) + B3*SIN(3d0*R)

         ! NDYSTEP is # of chemistry time steps
         NDYSTEP = INT( 24d0 * 3600d0 / GET_TS_CHEM() )

         ! Loop forward through NDYSTEP "fake" timesteps for this day
         DO N = 1, NDYSTEP

            ! Increment GMT (hours) to midpoint of next timestep
            IF ( N > 1 ) GMT_MID = GMT_MID + TS_SUN / 3600d0

            ! Loop over surface grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, YMID_R, TIMLOC, AHR, SUNTMP_MID )
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! Zero SUNTMP_MID
               SUNTMP_MID = 0d0

               ! Grid box latitude center [radians]
               YMID_R = GET_YMID_R( I, J, 1 )

               ! Local time at box (I,J) [hours]
               TIMLOC = GET_LOCALTIME( I, J, 1, GMT=GMT_MID )

               ! Hour angle at box (I,J) [radians]
               AHR = ABS( TIMLOC - 12d0 ) * 15d00 * PI_180

               !===========================================================
               ! The cosine of the solar zenith angle (SZA) is given by:
               !
               !  cos(SZA) = sin(LAT)*sin(DEC) +
               !  cos(LAT)*cos(DEC)*cos(AHR)
               !
               ! where LAT = the latitude angle,
               !       DEC = the solar declination angle,
               !       AHR = the hour angle, all in radians.
               !
               ! If SUNCOS < 0, then the sun is below the horizon, and
               ! therefore does not contribute to any solar heating.
               !===========================================================

               ! Compute Cos(SZA)
               SUNTMP_MID = sin(YMID_R) * sin(DEC) +
     &                      cos(YMID_R) * cos(DEC) * cos(AHR)

               ! TCOSZ is the sum of SUNTMP_MID at location (I,J)
               ! Do not include negative values of SUNTMP_MID

               TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP_MID, 0d0 )

            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDDO

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! Return to calling program
      END SUBROUTINE CALC_DIURNAL_CH4COCO2

!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_global_ch4
!
! !DESCRIPTION: Subroutine INIT\_GLOBAL\_CH4 allocates and zeroes module 
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_CH4COCO2( am_I_Root, Input_Opt,
     &  State_Diag, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD
#endif
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : Ind_
      USE State_Diag_Mod,     ONLY : DgnState
      USE ERROR_MOD,          ONLY : ALLOC_ERR

!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  This routine is called from GC_INIT_EXTRA (in GeosCore/input_mod.f)
! 
! !REVISION HISTORY:
!  (1 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  12 Feb 2014 - K. Wecht    - Disable CH4 budget diagnostic (bracket the 
!                              code out with #ifdef blocks so it can be used
!  11 Apr 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC arguments
!  16 Jun 2016 - M. Sulprizio- Now define IDT_CH4 locally
!  20 Jun 2016 - R. Yantosca - Rename IDTCH4 to id_CH4 for consistency
!  20 Jun 2016 - R. Yantosca - Now stop run if id_CH4 is undefined
!  09 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Strings
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: ThisLoc

      LOGICAL, SAVE    :: IS_INIT = .FALSE.
      INTEGER          :: AS

      ! For values from Input_Opt
      LOGICAL          :: LFOSSIL
      LOGICAL          :: LCHEMCO2
      LOGICAL          :: LBIODIURNAL
      LOGICAL          :: LBIONETCLIM
      LOGICAL          :: LOCEAN
      LOGICAL          :: LSHIP
      LOGICAL          :: LPLANE
      LOGICAL          :: LFFBKGRD
      LOGICAL          :: LBIOSPHTAG,  LFOSSILTAG,  LSHIPTAG
      LOGICAL          :: LPLANETAG

      !=================================================================
      ! INIT_GLOBAL_CH4 begins here!
      !=================================================================
      ! Assume Success
      RC      = GC_SUCCESS
      ErrMsg  = ''
      ThisLoc = ' -> INIT_CH4 (in module GeosCore/global_ch4_mod.F)'

      ! Define species ID flag
      id_CH4  = Ind_('CH4')

      ! Make sure CH4 is a defined species (bmy, 6/20/16)
      IF ( id_CH4 <= 0 ) THEN
         ErrMsg = 'CH4 is an undefined species!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ALLOCATE( BAIRDENS( IIPAR, JJPAR, LLPAR ), STAT=RC )
      CALL GC_CheckVar( 'ch4coco2_mod.F:BAIRDENS', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      BAIRDENS = 0e+0_fp

      ALLOCATE( CH4_EMIS_J( IIPAR, JJPAR, 15 ), STAT=RC )
      CALL GC_CheckVar( 'ch4coco2_mod.F:CH4_EMIS', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      CH4_EMIS_J = 0e+0_fp

      ! Initialize tropoch4 (counts total decay of CH4 due to OH)
      TROPOCH4 = 0e+0_fp

      ! Set a flag to denote if bpch diagnostics are activated
#if defined( BPCH_DIAG )
      Do_ND43 = ( ND43 > 0 )
#else
      Do_ND43 = .FALSE.
#endif

      !=================================================================
      ! INIT CO begins here!
      !=================================================================

      ! Allocate SUMISOPCO -- array for CO from isoprene
      ALLOCATE( SUMISOPCO( IIPAR, JJPAR ), STAT=RC )
      CALL GC_CheckVar( 'tagged_co_mod.F:SUMISOPCO', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      SUMISOPCO = 0.0_fp

      ! Allocate SUMMONOCO -- array for CO from monoterpenes
      ALLOCATE( SUMMONOCO( IIPAR, JJPAR ), STAT=RC )
      CALL GC_CheckVar( 'tagged_co_mod.F:SUMMONOCO', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      SUMMONOCO = 0.0_fp

      ! Allocate SUMCH3OH -- array for CO from CH3OH
      ALLOCATE( SUMCH3OHCO( IIPAR, JJPAR ), STAT=RC )
      CALL GC_CheckVar( 'tagged_co_mod.F:SUMCH3OHCO', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      SUMCH3OHCO = 0.0_fp

      ! Allocate SUMACETCO -- array for CO from isoprene
      ALLOCATE( SUMACETCO( IIPAR, JJPAR ), STAT=RC )
      CALL GC_CheckVar( 'tagged_co_mod.F:SUMACETCO', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      SUMACETCO = 0.0_fp

      ! Allocate TCOSZ -- array for sum of COS(SZA)
      ALLOCATE( TCOSZ( IIPAR, JJPAR ), STAT=RC )
      CALL GC_CheckVar( 'tagged_co_mod.F:TCOSZ', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      TCOSZ = 0.0_fp

      ! Allocate CH4_OH_TROP -- array for CH4 loss by OH in tropo
      ALLOCATE( CH4_OH_TROP( IIPAR, JJPAR, LLPAR ), STAT=RC )
      CALL GC_CheckVar( 'tagged_co_mod.F:CH4OHTROP', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      CH4_OH_TROP = 0.0_fp

      ! Allocate CH4_OH_STRAT -- array for CH4 loss by OH in strat
      ALLOCATE( CH4_OH_STRAT( IIPAR, JJPAR, LLPAR ), STAT=RC )
      CALL GC_CheckVar( 'tagged_co_mod.F:CH4OHSTRAT', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      CH4_OH_STRAT = 0.0_fp

      !=================================================================
      ! INIT_CO2 begins here!
      !=================================================================

      ! Exit if we have already intialized 
      IF ( IS_INIT ) RETURN
        
      ! Copy values from Input_Opt
      LFOSSIL     = Input_Opt%LFOSSIL
      LCHEMCO2    = Input_Opt%LCHEMCO2
      LBIODIURNAL = Input_Opt%LBIODIURNAL
      LBIONETCLIM = Input_Opt%LBIONETCLIM
      LOCEAN      = Input_Opt%LOCEAN
      LSHIP       = Input_Opt%LSHIP
      LPLANE      = Input_Opt%LPLANE
      LFFBKGRD    = Input_Opt%LFFBKGRD
      LBIOSPHTAG  = Input_Opt%LBIOSPHTAG
      LFOSSILTAG  = Input_Opt%LFOSSILTAG
      LSHIPTAG    = Input_Opt%LSHIPTAG
      LPLANETAG   = Input_Opt%LPLANETAG

      ! Reset IS_INIT flag
      IS_INIT = .TRUE.

      END SUBROUTINE INIT_CH4COCO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_global_ch4
!
! !DESCRIPTION: Subroutine CLEANUP\_GLOBAL\_CH4 deallocates module arrays.
!  (bmy, 1/16/01)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_CH4COCO2( am_I_Root, RC )
!
! !USES:
!
      USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY:
!  (1 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  12 Feb 2014 - K. Wecht    - Disable CH4 budget diagnostic (bracket the 
!                              code out with #ifdef blocks so it can be used)
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! CLEANUP_GLOBAL_CH4 begins here!
      !=================================================================
      ! Initialize
      RC = GC_SUCCESS

      ! Deallocate variables
      IF ( ALLOCATED( BAIRDENS ) ) THEN
         DEALLOCATE( BAIRDENS, STAT=RC )     
         CALL GC_CheckVar( 'global_br_mod.F:BAIRDENS', 2, RC )
         RETURN
      ENDIF

      IF ( ALLOCATED( CH4_EMIS_J ) ) THEN
         DEALLOCATE( CH4_EMIS_J, STAT=RC )     
         CALL GC_CheckVar( 'global_br_mod.F:CH4_EMIS', 2, RC )
         RETURN
      ENDIF

      ! Deallocate variables
      IF ( ALLOCATED( SUMISOPCO ) ) THEN
         DEALLOCATE( SUMISOPCO, STAT=RC )
         CALL GC_CheckVar( 'tagged_co_mod.F:SUMISOPCO', 2, RC )
         RETURN
      ENDIF

      IF ( ALLOCATED( SUMMONOCO ) ) THEN
         DEALLOCATE( SUMMONOCO , STAT=RC )
         CALL GC_CheckVar( 'tagged_co_mod.F:SUMMONOCO', 2, RC )
         RETURN
      ENDIF

      IF ( ALLOCATED( SUMCH3OHCO ) ) THEN
         DEALLOCATE( SUMCH3OHCO, STAT=RC )
         CALL GC_CheckVar( 'tagged_co_mod.F:SUMCH3OHCO', 2, RC )
         RETURN
      ENDIF

      IF ( ALLOCATED( SUMACETCO ) ) THEN
         DEALLOCATE( SUMACETCO, STAT=RC )
         CALL GC_CheckVar( 'tagged_co_mod.F:SUMACETCO', 2, RC )
         RETURN
      ENDIF

      IF ( ALLOCATED( TCOSZ ) ) THEN
         DEALLOCATE( TCOSZ, STAT=RC )
         CALL GC_CheckVar( 'tagged_co_mod.F:TCOSZ', 2, RC )
         RETURN
      ENDIF

      IF ( ALLOCATED( CH4_OH_TROP ) ) THEN
         DEALLOCATE( CH4_OH_TROP, STAT=RC )
         CALL GC_CheckVar( 'tagged_co_mod.F:CH4OHTROP', 2, RC )
         RETURN
      ENDIF

      IF ( ALLOCATED( CH4_OH_STRAT ) ) THEN
         DEALLOCATE( CH4_OH_STRAT, STAT=RC )
         CALL GC_CheckVar( 'tagged_co_mod.F:CH4OHSTRAT', 2, RC )
         RETURN
      ENDIF

      ! Free pointers
      IF ( ASSOCIATED( BOH      ) ) BOH     => NULL()
      IF ( ASSOCIATED( BCl      ) ) BCl     => NULL()
      IF ( ASSOCIATED( CH4LOSS  ) ) CH4LOSS => NULL()

      ! Free pointers
      GMI_PROD_CO => NULL()
      GMI_LOSS_CO => NULL()
      OH          => NULL()
      SFC_CH4     => NULL()


      END SUBROUTINE CLEANUP_CH4COCO2
!EOC
      SUBROUTINE SETKPPVALS___( State_Chm, State_Met, I, J, L, &
           BAIRDENS, BOH, BCl )
        USE State_Chm_Mod,      ONLY : ChmState
        USE State_Met_Mod,      ONLY : MetState
        USE TIME_MOD,           ONLY : GET_TS_CHEM

        USE GCKPP_GLOBAL
        USE GCKPP_PARAMETERS

!
! !INPUT PARAMETERS:
!
        TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
        TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
        INTEGER,        INTENT(IN)  :: I, J, L
        REAL(fp),       INTENT(IN)  :: BOH ! OH
        REAL(fp),       INTENT(IN)  :: BCl ! Cl
        REAL(fp),       INTENT(IN)  :: BAIRDENS

        REAL(fp) :: DT, FAC_DIURNAL
        
        ! Chemistry timestep in seconds
        DT = GET_TS_CHEM()

! Species
        ! molec/cm3
        C(ind_CH4)     = State_Chm%Species(I,J,L,1) &
             / State_Met%AIRVOL(I,J,L) / 1e+6_fp * XNUMOL_CH4
        C(ind_CO)      = State_Chm%Species(I,J,L,Ind('CO')) &
             / State_Met%AIRVOL(I,J,L) / 1e+6_fp * XNUMOL_CO 
        C(ind_CH4_E)   = 1.E0 ! Dummy quantity.
        C(ind_NMVOC_E) = 1.E0 ! Dummy quantity.
           
        IF ( State_Met%InStratosphere(I,J,L) ) THEN
           TROPK = 0.d0
           TROP  = 0.d0 ! Toggle

        ELSE ! In the trop
           ! Temperature (K)
           TEMP            = State_Met%T(I,J,L)

           ! Cosine of the solar zenith angle [unitless]
           SUNCOS = State_Met%SUNCOSmid(I,J)
           
           ! Scaling factor for diurnal cycles - zero at night
           IF ( State_Met%SUNCOSmid(I,J) > 0.0_fp .and. &
                TCOSZ(I,J) > 0.0_fp ) THEN
              FAC_DIURNAL = ( State_Met%SUNCOSmid(I,J) / TCOSZ(I,J) ) *
              &                       ( 86400.0_fp / DT )
           ELSE
              FAC_DIURNAL = 0.0_fp
           ENDIF

           C(ind_OH_E)    = BOH * XNUMOL_OH / CM3PERM3 * FAC_DIURNAL
           C(ind_Cl_E)    = BCl * BAIRDENS(I,J,L) * 1e-9_fp
        
! Rates
           STRATK = 0.d0
           ! Troposphere
           ! CH4 + OH: This rate > 0 if in trop. Otherwise. 0
           TROP            = 1.e0
           TROP_NMVOC_RATE = PCO_NMVOC_MCM3S(I,J,L)
           TOGGLE ! LPCO_CH4 TRUE = 1, FALSE = 0
           ! Stratosphere
           STRAT_CH4_RATE  = CH4LOSS(I,J,L)
        ENDIF ! In the strat/trop

      END SUBROUTINE SETKPPVALS___
    END MODULE CCYCLECHEM_MOD
