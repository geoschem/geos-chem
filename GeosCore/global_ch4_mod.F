!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: global_ch4_mod.F
!
! !DESCRIPTION: Module GLOBAL\_CH4\_MOD contains variables and routines for
!  simulating CH4 chemistry in the troposphere.
!\\
!\\
! !INTERFACE: 
!
      MODULE GLOBAL_CH4_MOD
!
! !USES:
!
      USE HCO_ERROR_MOD       ! For HEMCO error reporting
      USE PhysConstants, ONLY : AVO, AIRMW
      USE PRECISION_MOD       ! For GEOS-Chem Precision (fp, f4, f8)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: EMISSCH4
      PUBLIC :: CHEMCH4
      PUBLIC :: INIT_GLOBAL_CH4
      PUBLIC :: CLEANUP_GLOBAL_CH4
!
! !PUBLIC DATA MEMBERS:
!
      REAL(fp), PARAMETER,   PUBLIC :: XNUMOL_CH4 = AVO / 16d-3 ! hard-coded MW

      ! Make CH4_EMIS now public so that it can be used by vdiff_mod.F90
      ! Methane emissions units are [kg/m2/s]
      REAL(fp),  ALLOCATABLE, PUBLIC :: CH4_EMIS(:,:,:)
! 
! !REVISION HISTORY: 
!  17 Jan 2001- J. Wang, B. Duncan, R. Yantosca -- Initial version
!  (1 ) Merged routines from jsw's CH4 code  into "global_ch4_mod.f" 
!        (bmy, 1/16/01)
!  (2 ) XNUMOL_CH4 and TCH4 have to be public - all other variables can
!        be made private, so as not to conflict with other common-block
!        definitions (bmy, 1/17/01)
!  (3 ) Minor fixes from jsw added (jsw, bmy, 2/17/01)
!  (4 ) Removed some F90 module references from EMISSCH4 (bmy, 3/20/01)
!  (5 ) Eliminate obsolete commented-out code (bmy, 4/20/01)
!  (6 ) Updated comments (bmy, 9/4/01)
!  (7 ) Fixes for binary punch file in READ_COPROD (bmy, 9/26/01)
!  (8 ) Removed obsolete code from READ_COPROD (bmy, 10/24/01)
!  (9 ) Minor bug fixes for compilation on ALPHA (bmy, 11/15/01)
!  (10) Eliminate obsolete code from 11/01 (bmy, 2/27/02)
!  (11) Now eliminate PS from the arg list to CH4_AVGTP (4/11/02)
!  (12) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (13) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (14) Now reference "file_mod.f".  Also removed obsolete code. (bmy, 6/27/02)
!  (15) Now references "pressure_mod.f" (bmy, 8/21/02)
!  (16) Now reference AD and T from "dao_mod.f".  Now reference "error_mod.f".
!        Remove obsolete code from various routines.  Remove reference to
!        header file "comtrid.h" -- it's not used. (bmy, 11/6/02)
!  (17) Minor bug fix in FORMAT statements (bmy, 3/23/03)
!  (18) Now references "grid_mod.f" and "time_mod.f" (bmy, 3/27/03)
!  (19) Updates to GET_GLOBAL_CH4 (bmy, 7/1/03)
!  (20) Now references "directory_mod.f", "tracer_mod.f", and "diag_oh_mod.f"
!        (bmy, 7/20/04)
!  (21) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (22) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (23) Updated CH4 simulation (kjw, cph, ccarouge, 10/1/09)
!  (24) Added modifications for MERRA (bmy, 8/13/10)
!  08 Feb 2012 - R. Yantosca - Added modifications for GEOS-5.7.x
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  12 Feb 2014 - K. Wecht    - Add modifications for 0.25 x 0.3125 NA grid
!  12 Feb 2014 - K. Wecht    - Disable CH4 budget diagnostic (bracket the 
!                              code out with #ifdef blocks so it can be used)
!  20 Jun 2014 - R. Yantosca - Remove obsolete emissions code; we now use HEMCO
!  09 Sep 2014 - C. Keller   - Implemented HEMCO for emissions. 
!  16 Sep 2014 - C. Keller   - Now also get chemical fields (OH and CH4 prod)
!                              through HEMCO. Removed READ_COPROD and 
!                              READ_CH4LOSS (now obsolete)
!  16 Sep 2014 - C. Keller   - Removed all references to Clarisa's OH clima-
!                              tology (not used anywhere)
!  14 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  04 Mar 2015 - R. Yantosca - Now use REAL(f4) for HCO_GetPtr pointer args
!  06 Jan 2016 - E. Lundgren - Use global physical parameters
!  16 Jun 2016 - M. Sulprizio- Now define IDTCH4 locally
!  20 Jun 2016 - R. Yantosca - Rename IDTCH4 to id_CH4 for consistency
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!  16 Jun 2017 - M. Sulprizio- Remove routines CH4_AVGTP, GET_CH4_ANTHRO, and
!                              GET_CH4_BIOG - they are not used anymore
!  18 Oct 2017 - M. Sulprizio- Remove CH4_BUDGET code. It did not work anymore 
!                              and was in need of rewriting. This information
!                              can be calculated using GEOS-Chem output.
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
      SUBROUTINE EMISSCH4( am_I_Root, Input_Opt,
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
      LOGICAL                  :: ITS_A_CH4_SIM
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
      ITS_A_CH4_SIM  = Input_Opt%ITS_A_CH4_SIM
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

      IF ( ITS_A_CH4_SIM .and. LPRT ) THEN
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
      CH4_EMIS(:,:,:) = 0e+0_fp

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
         CH4_EMIS(:,:,2) =  Ptr2D(:,:)
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
         CH4_EMIS(:,:,3) =  Ptr2D(:,:)
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
         CH4_EMIS(:,:,4) =  Ptr2D(:,:)
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
         CH4_EMIS(:,:,5) =  Ptr2D(:,:)
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
         CH4_EMIS(:,:,6) =  Ptr2D(:,:)
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
         CH4_EMIS(:,:,7) =  Ptr2D(:,:)
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
         CH4_EMIS(:,:,8) =  Ptr2D(:,:)
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
         CH4_EMIS(:,:,9) =  Ptr2D(:,:)
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
         CH4_EMIS(:,:,10) =  Ptr2D(:,:)
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
         CH4_EMIS(:,:,11) =  Ptr2D(:,:)
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
         CH4_EMIS(:,:,12) =  Ptr2D(:,:)
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
         CH4_EMIS(:,:,13) =  Ptr2D(:,:)
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
         CH4_EMIS(:,:,14) =  Ptr2D(:,:)
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
         CH4_EMIS(:,:,15) =  Ptr2D(:,:) * -1.0_hp
      ENDIF
      Ptr2D => NULL()

      ! =================================================================
      ! Total emission: sum of all emissions - (2*soil absorption)
      ! We have to substract soil absorption twice because it is added 
      ! to other emissions in the SUM function. (ccc, 7/23/09)
      ! =================================================================
      CH4_EMIS(:,:,1) = SUM(CH4_EMIS, 3) - (2 * CH4_EMIS(:,:,15))

      IF ( prtDebug ) THEN
         WRITE(*,*) 'CH4_EMIS (kg/m2/s):'
         WRITE(*,*) 'Total        : ', SUM(CH4_EMIS(:,:,1))
         WRITE(*,*) 'Oil          : ', SUM(CH4_EMIS(:,:,2))
         WRITE(*,*) 'Gas          : ', SUM(CH4_EMIS(:,:,3))
         WRITE(*,*) 'Coal         : ', SUM(CH4_EMIS(:,:,4))
         WRITE(*,*) 'Livestock    : ', SUM(CH4_EMIS(:,:,5))
         WRITE(*,*) 'Landfills    : ', SUM(CH4_EMIS(:,:,6))
         WRITE(*,*) 'Wastewater   : ', SUM(CH4_EMIS(:,:,7))
         WRITE(*,*) 'Rice         : ', SUM(CH4_EMIS(:,:,8))
         WRITE(*,*) 'Other anth   : ', SUM(CH4_EMIS(:,:,9))
         WRITE(*,*) 'Biomass burn : ', SUM(CH4_EMIS(:,:,10))
         WRITE(*,*) 'Wetlands     : ', SUM(CH4_EMIS(:,:,11))
         WRITE(*,*) 'Seeps        : ', SUM(CH4_EMIS(:,:,12))
         WRITE(*,*) 'Lakes        : ', SUM(CH4_EMIS(:,:,13))
         WRITE(*,*) 'Termites     : ', SUM(CH4_EMIS(:,:,14))
         WRITE(*,*) 'Soil absorb  : ', SUM(CH4_EMIS(:,:,15))
      ENDIF

      IF ( ITS_A_CH4_SIM .and. LPRT ) THEN
         print*,'END SUBROUTINE: EMISSCH4'
      ENDIF

      END SUBROUTINE EMISSCH4
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
      SUBROUTINE CHEMCH4( am_I_Root, Input_Opt,  State_Met,
     &                    State_Chm, State_Diag, RC         )
!
! !USES:
!
      USE CMN_SIZE_MOD
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD
      USE DIAG_MOD,           ONLY : AD43
#endif
      USE ErrCode_Mod
      USE HCO_INTERFACE_MOD,  ONLY : HcoState
      USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Diag_Mod,     ONLY : DgnState
      USE State_Met_Mod,      ONLY : MetState
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
      REAL(fp)           :: BOXVL  

      ! Number of days per month
      INTEGER            :: NODAYS(12) = (/ 31, 28, 31, 30, 
     &                                      31, 30, 31, 31, 
     &                                      30, 31, 30, 31 /)

      ! For fields from Input_Opt
      LOGICAL            :: LSPLIT
      LOGICAL            :: LPRT

      ! Pointers
      REAL(fp), POINTER  :: Spc(:,:,:,:)

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg
      CHARACTER(LEN=255) :: ThisLoc

      !=================================================================
      ! CHEMCH4 begins here!
      !=================================================================

      ! Assume success
      RC      = GC_SUCCESS
      ErrMsg  = ''
      ThisLoc = ' -> at CHEMCH4 (in module GeosCore/global_ch4_mod.F)'

      ! Copy values from Input_Opt
      LSPLIT  = Input_Opt%LSPLIT
      LPRT    = Input_Opt%LPRT

      ! Point to the chemical species 
      Spc     => State_Chm%Species

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
         CALL HCO_GetPtr( am_I_Root, HcoState, 'GLOBAL_OH', BOH, RC )

         ! Trap potential errors
         IF ( RC /= HCO_SUCCESS ) THEN
            ErrMsg = 'Cannot get pointer to HEMCO field GLOBAL_OH!'
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
            IF ( State_Met%InChemGrid(I,J,L) ) THEN

#if defined( BPCH_DIAG )
               !--------------------------------------------------------
               ! ND43 (bpch) diagnostic
               ! OH concentration in [molec/cm3] after chemistry
               !--------------------------------------------------------
               IF ( Do_ND43 .and. L <= LD43 ) THEN
                  AD43(I,J,L,1) = AD43(I,J,L,1) +
     &                            ( BOH(I,J,L) * XNUMOL_OH / CM3PERM3 )
               ENDIF
#endif


#if defined( NC_DIAG )
               !--------------------------------------------------------
               ! HISTORY (aka netCDF diagnostics)
               ! OH concentration in [molec/cm3] after chemistry
               !--------------------------------------------------------
               IF ( State_Diag%Archive_OHconcAfterChem ) THEN
                  State_Diag%OHconcAfterChem(I,J,L) =
     &                 ( BOH(I,J,L) * XNUMOL_OH / CM3PERM3 )
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
      CALL CH4_OHSAVE( State_Met, State_Chm )

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
      CALL CH4_DECAY( am_I_Root, Input_Opt,  State_Met, 
     &                State_Chm, State_Diag, RC         )

      !=================================================================
      ! (7) calculate CH4 chemistry in layers above tropopause.
      !=================================================================
      CALL CH4_STRAT( am_I_Root, Input_Opt,  State_Met, 
     &                State_Chm, State_Diag, RC         )

      !=================================================================
      ! (8) distribute the chemistry sink from total CH4 to other CH4 
      !     species. (ccc, 2/10/09)
      !=================================================================
      IF ( LSPLIT ) THEN
         CALL CH4_DISTRIB( PREVCH4, Input_Opt, State_Chm )
      ENDIF

      ! Free pointer
      Spc => NULL()

      ! Set FIRSTCHEM to FALSE
      FIRSTCHEM = .FALSE.

      END SUBROUTINE CHEMCH4
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
      SUBROUTINE CH4_DECAY( am_I_Root, Input_Opt,  State_Met, 
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
      REAL(fp)          :: KRATE, C_OH
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
!$OMP+PRIVATE( C_Cl, KRATE_Cl )
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

            ! OH in [molec/cm3]
            ! BOH is imported from HEMCO in units of kg/m3, convert
            !  here to molec/cm3 (ckeller, 9/16/2014)
            C_OH = BOH(I,J,L) * XNUMOL_OH / CM3PERM3

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

      END SUBROUTINE CH4_DECAY
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
      SUBROUTINE CH4_OHSAVE( State_Met, State_Chm )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE DIAG_OH_MOD,        ONLY : DO_DIAG_OH_CH4
      USE GC_GRID_MOD,        ONLY : GET_AREA_CM2
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
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
      REAL(fp)          :: C_OH 

      ! Pointers
      REAL(fp), POINTER :: Spc(:,:,:,:)

      !=================================================================
      ! CH4_OHSAVE begins here!
      !
      ! (1) Pass OH mass, total air mass, and  to "diag_oh_mod.f"
      ! (2) ND59: Diagnostic for CH3CCl3 calculation
      !=================================================================

      ! Point to chemical species array [kg]
      Spc => State_Chm%Species

      ! Calculate OH mass and total air mass
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, BOXVL, C_OH, OHMASS, MASST, KCLO, LOSS, KCH4 )
!$OMP+PRIVATE( CH4TROPMASS, CH4MASS, CH4LOSE, CH4EMIS, AREA_M2 )
      DO L = 1, LLPAR
      DO J = 1, JJPAR 
      DO I = 1, IIPAR 

         ! Only process tropospheric boxes (bmy, 4/17/00)
         IF ( State_Met%InChemGrid(I,J,L) ) THEN

            ! Grid box volume [cm3]
            BOXVL  = State_Met%AIRVOL(I,J,L) * 1e+6_fp

            ! OH concentration in molec/cm3. BOH is imported
            ! from HEMCO in kg/m3 (ckeller, 9/16/2014)
            C_OH = BOH(I,J,L) * XNUMOL_OH / CM3PERM3

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
               CH4EMIS  = CH4_EMIS(I,J,1)
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

      END SUBROUTINE CH4_OHSAVE
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
      SUBROUTINE CH4_STRAT( am_I_Root, Input_Opt,  State_Met, 
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

      END SUBROUTINE CH4_STRAT
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
      SUBROUTINE CH4_DISTRIB( PREVCH4, Input_Opt, State_Chm )
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
      DO NA = 2, nAdvect

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

      END SUBROUTINE CH4_DISTRIB
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
      SUBROUTINE INIT_GLOBAL_CH4( am_I_Root, Input_Opt, State_Diag, RC )
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
      CALL GC_CheckVar( 'global_ch4_mod.F:BAIRDENS', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      BAIRDENS = 0e+0_fp

      ALLOCATE( CH4_EMIS( IIPAR, JJPAR, 15 ), STAT=RC )
      CALL GC_CheckVar( 'global_ch4_mod.F:CH4_EMIS', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      CH4_EMIS = 0e+0_fp

      ! Initialize tropoch4 (counts total decay of CH4 due to OH)
      TROPOCH4 = 0e+0_fp

      ! Set a flag to denote if bpch diagnostics are activated
#if defined( BPCH_DIAG )
      Do_ND43 = ( ND43 > 0 )
#else
      Do_ND43 = .FALSE.
#endif

      END SUBROUTINE INIT_GLOBAL_CH4
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
      SUBROUTINE CLEANUP_GLOBAL_CH4( am_I_Root, RC )
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

      IF ( ALLOCATED( CH4_EMIS ) ) THEN
         DEALLOCATE( CH4_EMIS, STAT=RC )     
         CALL GC_CheckVar( 'global_br_mod.F:CH4_EMIS', 2, RC )
         RETURN
      ENDIF

      ! Free pointers
      IF ( ASSOCIATED( BOH      ) ) BOH     => NULL()
      IF ( ASSOCIATED( BCl      ) ) BCl     => NULL()
      IF ( ASSOCIATED( CH4LOSS  ) ) CH4LOSS => NULL()

      END SUBROUTINE CLEANUP_GLOBAL_CH4
!EOC
      END MODULE GLOBAL_CH4_MOD
