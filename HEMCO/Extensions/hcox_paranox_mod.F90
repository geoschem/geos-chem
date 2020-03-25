!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_paranox_mod.F90
!
! !DESCRIPTION: Module HCOX\_PARANOX\_MOD contains routines to
! compute ship emissions and associated concentrations of NO, NO2, HNO3
! and O3 from NO ship emission data. This follows the implementation of
! the PARANOX ship plume model in GEOS-Chem.
!\\
!\\
! This module calculates production rates of NO, NO2, HNO3, and O3, as well
! as loss rates of O3 and HNO3. All fluxes are in kg species/m2/s. The
! O3 and HNO3 loss fluxes are not converted to a deposition velocity, but
! rather saved out as mass fluxes (kg/m2/s) into diagnostics
! 'PARANOX\_O3\_DEPOSITION\_FLUX' and 'PARANOX\_HNO3\_DEPOSITION\_FLUX',
! respectively. In order to use them, they must be imported explicitly via
! routine Diagn\_Get (from module hco\_diagn\_mod.F90). This approach avoids
! problems with uncrealistically high loss rates for loss ambient air
! concentrations of O3 or HNO3.
!\\
!\\
! The PARANOx look-up-table can be provided in netCDF or ASCII (txt) format.
! The latter is particularly useful for running PARANOx in an ESMF environment,
! where 7-dimensional netCDF files are currently not supported. The input data
! format can be specified in the HEMCO configuration file (in the PARANOx
! extensions section).
! The txt-files can be generated from the previously read netCDF data using
! subroutine WRITE\_LUT\_TXTFILE.
!\\
!\\
! References:
! \begin{itemize}
! \item Vinken, G. C. M., Boersma, K. F., Jacob, D. J., and Meijer, E. W.:
! Accounting for non-linear chemistry of ship plumes in the
! GEOS-Chem global chemistry transport model, Atmos. Chem. Phys., 11,
! 11707-11722, doi:10.5194/acp-11-11707-2011, 2011.
! \end{itemize}
!
! The initial look up tables (LUT) distributed with GEOS-Chem v9-01-03
! used 7 input variables: Temperature, J(NO2), J(O1D), solar elevation angles
! at emission time and 5 hours later, and ambient concentrations of NOx
! and O3. This version was documented by  Vinken et al. (2011). Subsequently,
! we added wind speed as an input variable. We also use J(OH) rather than J(O1D)
! to index the LUT (C. Holmes,  6 May 2013)
!
! The LUTs contain 3 quantities:
!     FracNOx : The fraction of NOx emitted from ships that remains as NOx
!               after 5 hours of plume aging. mol/mol
!     OPE     : Ozone production efficiency, mol(O3)/mol(HNO3)
!               The net production of O3 per mole of ship NOx oxidized over
!               5 hours of plume aging. Can be negative!
!               Defined as OPE = [ P(O3) - L(O3) ] / P(HNO3), where each P
!               and L term is an integral over 5 hours. Net O3 production
!               in the plume is E(NOx) * (1-FracNOx) * OPE, where E(NOx) is
!               the emission rate of NOx from the ship (e.g. units: mol/s).
!     MOE     : Methane oxidation efficiency, mol(CH4)/mol(NOx)
!               The net oxidation of CH4 per mole of NOx emitted from ships
!               over 5 hours of plume aging.
!               Defined as MOE = L(CH4) / E(NOx).
!\\
!\\
! The solar elevation angles 5 hours ago are calculated using HEMCO subroutine
! HCO\_GetSUNCOS. This is the same routine that is used to calculate the solar
! zenith angles for the current time.
!\\
!\\
! !INTERFACE:
!
MODULE HCOX_ParaNOx_MOD
!
! !USES:
!
  USE HCO_Error_MOD
  USE HCO_Diagn_MOD
  USE HCO_State_MOD,  ONLY : HCO_State
  USE HCOX_State_MOD, ONLY : Ext_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCOX_ParaNOx_Run
  PUBLIC  :: HCOX_ParaNOx_Init
  PUBLIC  :: HCOX_ParaNOx_Final
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !REMARKS:
!  Adapted from the code in GeosCore/paranox_mod.F prior to GEOS-Chem v10-01.
!
! !REVISION HISTORY:
!  06 Aug 2013 - C. Keller   - Initial version
!  03 Jun 2013 - C. Holmes   - Rewritten to include wind speed in the look-up
!                              table and to take input from netCDF
!  15 Oct 2013 - C. Keller   - Now a HEMCO extension
!  06 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  06 Jun 2014 - R. Yantosca - Now indended with F90 free-format
!  25 Jun 2014 - R. Yantosca - Now pass the look-up-table filenames
!  15 Jul 2014 - C. Holmes   - Make module variables allocatable, since they
!                              are used only in full chemistry simulations.
!  22 Jul 2014 - R. Yantosca - Added shadow copy of FAST-JX function FJXFUNC
!  28 Jul 2014 - C. Keller   - Now pass J-Values through ExtState. This makes
!                              the FJXFUNC shadow copy obsolete
!  13 Aug 2014 - C. Keller   - Added manual diagnostics
!  16 Oct 2014 - C. Keller   - Now store SUNCOSmid values internally over the
!                              past 5 hours and use these values for SUNCOSmid5.
!                              This is required for standalone mode.
!  05 Feb 2015 - C. Keller   - Modified to bring in the updates from Chris
!                              Holmes (input data in netCDF format, include
!                              wind speed, calculated dry deposition freq.
!                              using whole troposheric column mass).
!  23 Feb 2015 - C. Keller   - Historic j-values can now be provided through
!                              HEMCO configuration file.
!  10 Apr 2015 - C. Keller   - Now exchange deposition fluxes via diagnostics.
!                              Keep units of kg/m2/s for loss rates.
!  20 Apr 2016 - M. Sulprizio- Get J(OH) directly from FAST-JX and remove all
!                              references to J(O1D). In FlexChem, adjustment of
!                              photolysis rates are now done in routine
!                              PHOTRATE_ADJ (found in GeosCore/fast_jx_mod.F).
!  14 Oct 2016 - C. Keller   - Now use HCO_EvalFld instead of HCO_GetPtr.
!  12 Sep 2018 - C. Keller   - Added instance wrapper
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:

  ! Number of values for each variable in the look-up table
  INTEGER, PARAMETER ::  nT=4, nJ=4, nO3=4, nNOx=5, nSEA=12, nWS=5

  ! Now place all module variables in a lderived type object (for a linked
  ! list) so that we can have one instance per node in an MPI environment.
  TYPE :: MyInst

     ! Scalars
     INTEGER               :: Instance
     INTEGER               :: ExtNr
     INTEGER               :: IDTNO
     INTEGER               :: IDTNO2
     INTEGER               :: IDTHNO3
     INTEGER               :: IDTO3
     REAL*8                :: MW_O3
     REAL*8                :: MW_NO
     REAL*8                :: MW_NO2
     REAL*8                :: MW_HNO3
     REAL*8                :: MW_AIR

     ! Arrays
     REAL(hp), POINTER     :: ShipNO(:,:,:)

     ! For SunCosMid 5hrs ago
     REAL(hp), POINTER     :: SC5(:,:)

     ! Deposition fluxes in kg/m2/s
     REAL(sp), POINTER     :: DEPO3  (:,:)
     REAL(sp), POINTER     :: DEPHNO3(:,:)

     ! Reference values of variables in the look-up tables
     REAL*4                :: Tlev(nT)
     REAL*4                :: JNO2lev(nJ)
     REAL*4                :: O3lev(nO3)
     REAL*4                :: SEA0lev(nSEA)
     REAL*4                :: SEA5lev(nSEA)
     REAL*4                :: JRATIOlev(nJ)
     REAL*4                :: NOXlev(nNOx)
     REAL*4                :: WSlev(nWS)

     ! Look-up tables currently used in GEOS-Chem (likely in v10-01)
     ! Described by Holmes et al. (2014), now includes effects of wind speed
     ! Last two digits in LUT names indicate wind speed in m/s
     REAL(sp), POINTER     :: FRACNOX_LUT02(:,:,:,:,:,:,:)
     REAL(sp), POINTER     :: FRACNOX_LUT06(:,:,:,:,:,:,:)
     REAL(sp), POINTER     :: FRACNOX_LUT10(:,:,:,:,:,:,:)
     REAL(sp), POINTER     :: FRACNOX_LUT14(:,:,:,:,:,:,:)
     REAL(sp), POINTER     :: FRACNOX_LUT18(:,:,:,:,:,:,:)
     REAL(sp), POINTER     :: OPE_LUT02    (:,:,:,:,:,:,:)
     REAL(sp), POINTER     :: OPE_LUT06    (:,:,:,:,:,:,:)
     REAL(sp), POINTER     :: OPE_LUT10    (:,:,:,:,:,:,:)
     REAL(sp), POINTER     :: OPE_LUT14    (:,:,:,:,:,:,:)
     REAL(sp), POINTER     :: OPE_LUT18    (:,:,:,:,:,:,:)
     REAL(sp), POINTER     :: MOE_LUT02    (:,:,:,:,:,:,:)
     REAL(sp), POINTER     :: MOE_LUT06    (:,:,:,:,:,:,:)
     REAL(sp), POINTER     :: MOE_LUT10    (:,:,:,:,:,:,:)
     REAL(sp), POINTER     :: MOE_LUT14    (:,:,:,:,:,:,:)
     REAL(sp), POINTER     :: MOE_LUT18    (:,:,:,:,:,:,:)
     REAL(sp), POINTER     :: DNOX_LUT02   (:,:,:,:,:,:,:)
     REAL(sp), POINTER     :: DNOX_LUT06   (:,:,:,:,:,:,:)
     REAL(sp), POINTER     :: DNOX_LUT10   (:,:,:,:,:,:,:)
     REAL(sp), POINTER     :: DNOX_LUT14   (:,:,:,:,:,:,:)
     REAL(sp), POINTER     :: DNOX_LUT18   (:,:,:,:,:,:,:)

     ! Location and type of look up table data
     CHARACTER(LEN=255)    :: LutDir
     LOGICAL               :: IsNc

     TYPE(MyInst), POINTER :: NextInst => NULL()
  END TYPE MyInst

  ! Pointer to instances
  TYPE(MyInst), POINTER    :: AllInst => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_ParaNOx_Run
!
! !DESCRIPTION: Subroutine HCOX\_ParaNOx\_Run is the driver routine to
! calculate ship NOx emissions for the current time step. Emissions in
! [kg/m2/s] are added to the emissions array of the passed
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_ParaNOx_Run( ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_Calc_Mod, ONLY : HCO_CalcEmis
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER       :: ExtState    ! External data fields
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState    ! HEMCO State object
    INTEGER,         INTENT(INOUT) :: RC          ! Success or failure?

! !REVISION HISTORY:
!  06 Aug 2013 - C. Keller   - Initial Version
!  06 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  06 Jun 2014 - R. Yantosca - Now indended with F90 free-format
!  28 Jul 2014 - C. Keller   - Now call Hco_CalcEmis instead of Hco_Run.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL               :: DefScaleEmis
    CHARACTER(LEN=255)    :: MSG
    TYPE(MyInst), POINTER :: Inst

    !=================================================================
    ! HCOX_PARANOX_RUN begins here!
    !=================================================================

    ! Return if extension disabled
    IF ( ExtState%ParaNOx <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER(HcoState%Config%Err,'HCOX_ParaNOx_Run (hcox_paranox_mod.F90)', RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get local instance
    Inst => NULL()
    CALL InstGet ( ExtState%ParaNOx, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(MSG,*) 'Cannot find ParaNOx instance Nr. ', ExtState%ParaNOx
       CALL HCO_ERROR(HcoState%Config%Err,MSG,RC)
       RETURN
    ENDIF

    ! ----------------------------------------------------------------
    ! Use HEMCO core routines to get ship NO emissions
    ! ----------------------------------------------------------------

    ! Prepare HEMCO core run (Hco_CalcEmis):
    ! --> Set tracer and category range + extension number.
    ! Note: Set species min and max to the full range of species.
    ! For the ParaNox extension, emission fields of only one species
    ! should be defined. Hco_CalcEmis will exit w/ error if this is
    ! not the case.
    HcoState%Options%SpcMin =  1
    HcoState%Options%SpcMax = -1
    HcoState%Options%CatMin =  1
    HcoState%Options%CatMax = -1
    HcoState%Options%ExtNr  = Inst%ExtNr

    ! --> Define array to write emissions into. ShipNO is reset to
    ! zero within subroutine EVOLVE_PLUME, so no need to do this
    ! here.
!    ShipNO                      = 0.0d0
    HcoState%Options%AutoFillDiagn = .FALSE.
    HcoState%Options%FillBuffer    =  .TRUE.
    HcoState%Buffer3D%Val          => Inst%ShipNO

    ! Calculate ship NO emissions and write them into the ShipNO
    ! array [kg/m2/s].
    CALL HCO_CalcEmis( HcoState, .FALSE., RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Reset settings to standard
    HcoState%Buffer3D%Val          => NULL()
    HcoState%Options%FillBuffer    = .FALSE.
    HcoState%Options%ExtNr         = 0
    HcoState%Options%AutoFillDiagn = .TRUE.

    ! Calculate production rates of NO, HNO3 and O3 based upon ship NO
    ! emissions and add these values to the respective emission
    ! arrays.
    ! Note: For O3, it is possible to get negative emissions (i.e.
    ! deposition), in which case these values will be added to the
    ! drydep array.
    CALL Evolve_Plume( ExtState, Inst%ShipNO, HcoState, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Leave w/ success
    Inst => NULL()
    CALL HCO_Leave( HcoState%Config%Err, RC )

  END SUBROUTINE HCOX_ParaNOx_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Evolve_Plume
!
! !DESCRIPTION: Subroutine EVOLVE\_PLUME performs plume dilution and chemistry
!  of ship NO emissions for every grid box and writes the resulting NO, HNO3
!  and O3 emission (production) rates into State\_Chm%NomixS.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Evolve_Plume( ExtState, ShipNoEmis, HcoState, Inst, RC )
!
! !USES:
!
    USE HCO_Types_Mod,    ONLY : DiagnCont
    USE HCO_FluxArr_mod,  ONLY : HCO_EmisAdd
    USE HCO_FluxArr_mod,  ONLY : HCO_DepvAdd
    USE HCO_Clock_Mod,    ONLY : HcoClock_First
    USE HCO_Calc_Mod,     ONLY : HCO_CheckDepv
    USE HCO_GeoTools_Mod, ONLY : HCO_GetSUNCOS
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER        :: ExtState           ! External data
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(hp),        INTENT(INOUT)  :: ShipNoEmis(:,:,:)  ! Emissions
    TYPE(HCO_State), POINTER        :: HcoState           ! HEMCO State obj
    TYPE(MyInst),    POINTER        :: Inst               ! Local instance
    INTEGER,         INTENT(INOUT)  :: RC                 ! Success or failure
!
! !REVISION HISTORY:
!  06 Aug 2013 - C. Keller   - Initial Version
!  06 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  06 Jun 2014 - R. Yantosca - Now indended with F90 free-format
!  24 Jun 2014 - R. Yantosca - Now pass LUT_FILENAME to READ_PARANOX_LUT
!  22 Jul 2014 - R. Yantosca - Comment out debug print statements
!  28 Jul 2014 - C. Keller   - Now get J-values through ExtState
!  12 Aug 2014 - R. Yantosca - READ_PARANOX_LUT is now called from Init phase
!  10 Nov 2014 - C. Keller   - Added div-zero error trap for O3 deposition.
!  25 Nov 2014 - C. Keller   - Now convert NO fluxes to HNO3 and O3 using
!                              corresponding molecular weight ratios. Safe
!                              division check for O3 deposition calculation.
!  08 May 2015 - C. Keller   - Now read/write restart variables from here to
!                              accomodate replay runs in GEOS-5.
!  25 May 2015 - C. Keller   - Now calculate SC5 via HCO_GetSUNCOS
!  29 Mar 2016 - C. Keller   - Bug fix: archive O3 deposition as positive flux.
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!  12 May 2017 - C. Keller   - Force option ScaleEmis to off.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                  :: I, J, L
    LOGICAL                  :: ERR
    LOGICAL                  :: FILLED
    LOGICAL                  :: FIRST
    LOGICAL                  :: DefScaleEmis
    REAL(hp)                 :: iFlx, TMP
    CHARACTER(LEN=255)       :: MSG
    CHARACTER(LEN=1)         :: CHAR1

    ! Arrays
    REAL(hp), TARGET         :: FLUXNO  (HcoState%NX,HcoState%NY)
    REAL(hp), TARGET         :: FLUXNO2 (HcoState%NX,HcoState%NY)
    REAL(hp), TARGET         :: FLUXHNO3(HcoState%NX,HcoState%NY)
    REAL(hp), TARGET         :: FLUXO3  (HcoState%NX,HcoState%NY)
!    REAL(hp), TARGET         :: DEPO3   (HcoState%NX,HcoState%NY)
!    REAL(hp), TARGET         :: DEPHNO3 (HcoState%NX,HcoState%NY)

    ! Pointers
    REAL(hp), POINTER        :: Arr2D(:,:)

    ! For diagnostics
    REAL(hp), TARGET         :: DIAGN   (HcoState%NX,HcoState%NY,5)
    LOGICAL, SAVE            :: DODIAGN = .FALSE.
    CHARACTER(LEN=31)        :: DiagnName
    TYPE(DiagnCont), POINTER :: TmpCnt

    ! Paranox update
    REAL(dp)                 :: SHIP_FNOx, SHIP_DNOx, SHIP_OPE, SHIP_MOE
    REAL(dp)                 :: FNO_NOx
    REAL(hp)                 :: iMass
    REAL(hp)                 :: ExpVal

!------------------------------------------------------------------------------
!### DEBUG -- COMMENT OUT FOR NOW
!    ! testing only
!    REAL*8             :: FRAC, TOTPRES, DELTPRES
!    INTEGER            :: TOP
!    integer            :: ix, jx
!    logical, parameter :: add2hemco = .true.
!------------------------------------------------------------------------------

    !=================================================================
    ! EVOLVE_PLUME begins here!
    !=================================================================

    ! Enter
    CALL HCO_ENTER(HcoState%Config%Err,'Evolve_Plume (hcox_paranox_mod.F90)', RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Leave here if none of the tracers defined
    IF ( Inst%IDTNO <= 0 .AND. Inst%IDTO3 <= 0 .AND. Inst%IDTHNO3 <= 0 ) THEN
       RC = HCO_SUCCESS
       RETURN
    ENDIF

    ! Nullify
    Arr2D  => NULL()
    TmpCnt => NULL()

    ! ------------------------------------------------------------------
    ! First call: check for diagnostics to write and fill restart values
    ! ------------------------------------------------------------------
    FIRST = HcoClock_First( HcoState%Clock, .TRUE. )

    IF ( FIRST ) THEN
       ! See if we have to write out manual diagnostics
       IF ( .NOT. DoDiagn ) THEN
          DiagnName = 'PARANOX_NOXFRAC_REMAINING'
          CALL DiagnCont_Find ( HcoState%Diagn, -1, -1, -1, -1, -1, &
                                DiagnName, 0, DoDiagn, TmpCnt )
          TmpCnt => NULL()
       ENDIF
       IF ( .NOT. DoDiagn ) THEN
          DiagnName = 'PARANOX_O3_PRODUCTION'
          CALL DiagnCont_Find ( HcoState%Diagn, -1, -1, -1, -1, -1, &
                                DiagnName, 0, DoDiagn, TmpCnt )
          TmpCnt => NULL()
       ENDIF
       IF ( .NOT. DoDiagn ) THEN
          DiagnName = 'PARANOX_NO_PRODUCTION'
          CALL DiagnCont_Find ( HcoState%Diagn, -1, -1, -1, -1, -1, &
                                DiagnName, 0, DoDiagn, TmpCnt )
          TmpCnt => NULL()
       ENDIF
       IF ( .NOT. DoDiagn ) THEN
          DiagnName = 'PARANOX_TOTAL_SHIPNOX'
          CALL DiagnCont_Find ( HcoState%Diagn, -1, -1, -1, -1, -1, &
                                DiagnName, 0, DoDiagn, TmpCnt )
          TmpCnt => NULL()
       ENDIF
       IF ( .NOT. DoDiagn ) THEN
          DiagnName = 'PARANOX_OPE'
          CALL DiagnCont_Find ( HcoState%Diagn, -1, -1, -1, -1, -1, &
                                DiagnName, 0, DoDiagn, TmpCnt )
          TmpCnt => NULL()
       ENDIF
    ENDIF

    IF ( DoDiagn ) DIAGN(:,:,:) = 0.0_hp

    ! ------------------------------------------------------------------
    ! Update SC5
    ! ------------------------------------------------------------------
    ! SC5 holds the SUNCOS values of 5 hours ago.
    CALL HCO_getSUNCOS( HcoState, Inst%SC5, -5, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Error check
    ERR = .FALSE.

    ! Init
    FLUXNO   = 0.0_hp
    FLUXNO2  = 0.0_hp
    FLUXHNO3 = 0.0_hp
    FLUXO3   = 0.0_hp

    ! Deposition fluxes
    Inst%DEPO3    = 0.0_sp
    Inst%DEPHNO3  = 0.0_sp

!------------------------------------------------------------------------
!    ! Debug
!    print*, '### In EVOLVE_PLUME:'
!    print*, '### JOH: ',  SUM   ( ExtState%JOH%Arr%Val ),  &
!                          MAXVAL( ExtState%JOH%Arr%Val )
!    print*, '### JNO2: ', SUM   ( ExtState%JNO2%Arr%Val ),  &
!                          MAXVAL( ExtState%JNO2%Arr%Val )
!    print*, '### SC5 : ', SUM   ( SC5 ), MAXVAL(SC5)
!    print*, '### EMIS: ', SUM   ( SHIPNOEMIS(:,:,1) ), MAXVAL(SHIPNOEMIS(:,:,1))
!------------------------------------------------------------------------

    ! Loop over all grid boxes
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Note: there seems to be a problem with the OMP loop in that
    ! the species concentrations (O3molec, NOmolec, NO2molec)
    ! differ slightly in a few grid boxes. Don't know exactly what
    ! is going on here, but uncomment for now! Needs more
    ! evaluation and testing.
    !
    ! Now use #if defined( 0 ) to block of this code (bmy, 6/6/14)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!$OMP PARALLEL DO                                                   &
    !!!$OMP DEFAULT( SHARED )                                             &
    !!!$OMP PRIVATE( I, J, L,   MSG, iFlx, iMass,    TMP                ) &
    !!!$OMP PRIVATE( SHIP_FNOx, SHIP_DNOx, SHIP_OPE, SHIP_MOE, FNO_NOx  ) &
    !!!$OMP SCHEDULE( DYNAMIC )
    DO J = 1, HcoState%NY
    DO I = 1, HcoState%NX

       ! Skip if no ship emissions in this grid box
       IF ( ShipNoEmis(I,J,1) == 0d0 ) CYCLE

       ! Production Efficiency for ship emiss (gvinken,mpayer,2/7/12)
       ! Updated to include effects of wind speed (cdh, 3/25/2014)
       ! Updated for HEMCO (ckeller, 02/04/2015)
       CALL PARANOX_LUT( ExtState,  HcoState, Inst, I, J, RC, &
                         SHIP_FNOx, SHIP_DNOx, SHIP_OPE, SHIP_MOE )
       IF ( RC /= HCO_SUCCESS ) THEN
          ERR = .TRUE.; EXIT
       ENDIF

!       ! for debugging only
!       if(I==3.and.J==35)then
!          write(*,*) 'PARANOX ship emissions @',I,J
!          write(*,*) 'Emis [kg/m2/s]: ', ShipNoEmis(I,J,1)
!          write(*,*) 'SHIP_FNOx: ', SHIP_FNOx
!          write(*,*) 'SHIP_DNOx: ', SHIP_DNOx
!          write(*,*) 'SHIP_OPE : ', SHIP_OPE
!          write(*,*) 'SHIP_MOE : ', SHIP_MOE
!       endif

       ! Split the ship NOx emission into NO and NO2
       ! following the ambient ratio
       ! Now check for zero ambient air concentrations. In this
       ! case, arbitrarily emit everything as NO (ckeller, 04/10/15).
       IF ( ExtState%NO%Arr%Val (I,J,1) > 0.0_hp .OR. &
            ExtState%NO2%Arr%Val(I,J,1) > 0.0_hp       ) THEN
          FNO_NOx = (ExtState%NO%Arr%Val(I,J,1)/Inst%MW_NO) / &
                  ( (ExtState%NO%Arr%Val(I,J,1)/Inst%MW_NO) + &
                    (ExtState%NO2%Arr%Val(I,J,1)/Inst%MW_NO2) )
       ELSE
          FNO_NOx = 1.0_hp
       ENDIF

       !---------------------------
       ! Calculate NO emissions
       !---------------------------
       IF ( Inst%IDTNO > 0 ) THEN

           ! Of the total ship NOx, the fraction SHIP_FNOx
           ! survives after plume dilution and chemistry.
           ! FNO_NOx is the ratio of NO / NOx.
           ! Unit: kg/m2/s
           FLUXNO(I,J) = ShipNoEmis(I,J,1) * SHIP_FNOx * FNO_NOx

       ENDIF

       !---------------------------
       ! Calculate NO2 emissions
       !---------------------------
       IF ( Inst%IDTNO2 > 0 ) THEN

           ! NO2 emissions complement NO emissions, so that total NOx
           ! emissions are preserved.
           FLUXNO2(I,J) = ShipNoEmis(I,J,1) * SHIP_FNOx * (1.d0-FNO_NOx) &
                        * ( Inst%MW_NO2 / Inst%MW_NO )
       ENDIF

       !---------------------------
       ! Calculate HNO3 emissions
       !---------------------------
       IF ( Inst%IDTHNO3 > 0 ) THEN

          ! Of the total ship NOx, the fraction 1-SHIP_FNOx-SHIP_DNOx
          ! is converted to HNO3 during plume dilution and chemistry.
          ! Unit: kg/m2/s
          FLUXHNO3(I,J) = ShipNoEmis(I,J,1) * ( 1d0-SHIP_FNOx-SHIP_DNOx ) &
                        * ( Inst%MW_HNO3 / Inst%MW_NO )
       ENDIF

       !--------------------------------------------------------------------
       ! NOy deposition (as HNO3) from the sub-grid plume.
       ! The calculated deposition flux is in kg/m2/s, which has to be
       ! converted to 1/s. The species mass is either the species mass in
       ! the first grid box or of the entire PBL column, depending on the
       ! HEMCO setting 'PBL_DRYDEP'.
       !
       ! As of 4/10/15, exchange loss rates in original units of kg/m2/s.
       ! (ckeller)
       !--------------------------------------------------------------------
       IF ( (Inst%IDTHNO3 > 0) .AND. (SHIP_DNOx > 0.0_dp) ) THEN

          ! Deposition flux in kg/m2/s.
          Inst%DEPHNO3(I,J) = ShipNoEmis(I,J,1) * SHIP_DNOx * ( Inst%MW_HNO3 / Inst%MW_NO )
!          iFlx = ShipNoEmis(I,J,1) * SHIP_DNOx * ( MW_HNO3 / MW_NO )
!
!          ! Get mass of species. This can either be the total PBL
!          ! column mass or the first layer only, depending on the
!          ! HEMCO setting.
!          iMass = ExtState%HNO3%Arr%Val(I,J,1) &
!                * ExtState%FRAC_OF_PBL%Arr%Val(I,J,1)
!          IF ( HcoState%Options%PBL_DRYDEP ) THEN
!             DO L = 1, HcoState%NZ
!                IF ( ExtState%FRAC_OF_PBL%Arr%Val(I,J,L) == 0.0_hp ) EXIT
!                iMass = iMass + ( ExtState%HNO3%Arr%Val(I,J,L) *       &
!                                  ExtState%FRAC_OF_PBL%Arr%Val(I,J,L) )
!             ENDDO
!          ENDIF
!
!          ! Calculate deposition velocity (1/s) from flux
!          ! Now avoid div-zero error (ckeller, 11/10/2014).
!          IF ( iMass > TINY(1.0_hp) ) THEN
!             TMP = ABS(iFlx) * HcoState%Grid%AREA_M2%Val(I,J)
!
!             ! Check if it's safe to do division
!             IF ( (EXPONENT(TMP)-EXPONENT(iMass)) < MAXEXPONENT(TMP) ) THEN
!                DEPHNO3(I,J) = TMP / iMass
!             ENDIF
!
!             ! Check deposition velocity
!             CALL HCO_CheckDepv( HcoState, DEPHNO3(I,J), RC )
!          ENDIF

       ENDIF

       !---------------------------
       ! Calculate O3 emissions
       !---------------------------
       IF ( Inst%IDTO3 > 0 ) THEN

          ! Of the total ship NOx, the fraction
          ! (1-FRACTION_NOX)*INT_OPE is converted to O3 during
          ! plume dilution and chemistry.
          ! Unit: kg/m2/s
          iFlx = ShipNoEmis(I,J,1) * (1d0-SHIP_FNOx) * SHIP_OPE &
               * ( Inst%MW_O3 / Inst%MW_NO )

          ! For positive fluxes, add to emission flux array
          IF ( iFlx >= 0.0_hp ) THEN
             FLUXO3(I,J) = iFlx

          ! For negative fluxes, calculate deposition velocity based
          ! on current surface O3 concentration and pass to deposition
          ! array. See comment on dry dep calculation of HNO3 above.
          ! As of 4/10/15, exchange loss rates in original units of kg/m2/s.
          ! (ckeller)
          ELSE

             ! Deposition flux in kg/m2/s.
             ! Make sure ozone deposition flux is positive (ckeller, 3/29/16).
             !DEPO3(I,J) = iFlx
             Inst%DEPO3(I,J) = ABS(iFlx)

!             ! Get mass of species. This can either be the total PBL
!             ! column mass or the first layer only, depending on the
!             ! HEMCO setting.
!             iMass = ExtState%O3%Arr%Val(I,J,1) &
!                   * ExtState%FRAC_OF_PBL%Arr%Val(I,J,1)
!             IF ( HcoState%Options%PBL_DRYDEP ) THEN
!                DO L = 1, HcoState%NZ
!                   IF ( ExtState%FRAC_OF_PBL%Arr%Val(I,J,L) == 0.0_hp ) EXIT
!                   iMass = iMass + ( ExtState%O3%Arr%Val(I,J,L) *       &
!                                     ExtState%FRAC_OF_PBL%Arr%Val(I,J,L) )
!                ENDDO
!             ENDIF
!
!             ! Calculate deposition velocity (1/s) from flux
!             ! Now avoid div-zero error (ckeller, 11/10/2014).
!             IF ( iMass > TINY(1.0_hp) ) THEN
!                TMP = ABS(iFlx) * HcoState%Grid%AREA_M2%Val(I,J)
!
!                ! Check if it's safe to do division
!                IF ( (EXPONENT(TMP)-EXPONENT(iMass)) < MAXEXPONENT(TMP) ) THEN
!                   DEPO3(I,J) = TMP / iMass
!                ENDIF
!
!                ! Check deposition velocity
!                CALL HCO_CheckDepv( HcoState, DEPO3(I,J), RC )
!             ENDIF

          ENDIF
       ENDIF

       ! Eventually write out into diagnostics array
       IF ( DoDiagn ) THEN
          DIAGN(I,J,1) = SHIP_FNOx
          DIAGN(I,J,2) = SHIP_OPE
          DIAGN(I,J,3) = FLUXO3(I,J)
          DIAGN(I,J,4) = ShipNoEmis(I,J,1)
          DIAGN(I,J,5) = FLUXNO(I,J)
       ENDIF

       ! Reset ship NO emissions to zero. Will be refilled on next
       ! emission step!
       ShipNoEmis(I,J,1) = 0.0d0

    ENDDO !I
    ENDDO !J
    !!!$OMP END PARALLEL DO

    ! Error check
    IF ( ERR ) THEN
       RC = HCO_FAIL
       RETURN
    ENDIF

!------------------------------------------------------------------------
!    ! Debug
!    print*, '### In EVOLVE_PLUME (B):'
!    print*, '### DIAG 1: ',  SUM   ( DIAGN(:,:,1) ),  &
!                             MAXVAL( DIAGN(:,:,1) )
!    print*, '### DIAG 2: ',  SUM   ( DIAGN(:,:,2) ),  &
!                             MAXVAL( DIAGN(:,:,2) )
!    print*, '### DIAG 3: ',  SUM   ( DIAGN(:,:,3) ),  &
!                             MAXVAL( DIAGN(:,:,3) )
!    print*, '### DIAG 4: ',  SUM   ( DIAGN(:,:,4) ),  &
!                             MAXVAL( DIAGN(:,:,4) )
!    print*, '### DIAG 5: ',  SUM   ( DIAGN(:,:,5) ),  &
!                             MAXVAL( DIAGN(:,:,5) )
!------------------------------------------------------------------------

    !=======================================================================
    ! PASS TO HEMCO STATE AND UPDATE DIAGNOSTICS
    !=======================================================================

    ! Turn off emission scaling. We don't want the computed fluxes to be
    ! scaled any more. If a uniform scale factor is defined for NO, it
    ! has been applied to the ship NO emissions already (ckeller, 5/11/17).
    DefScaleEmis               = HcoState%Options%ScaleEmis
    HcoState%Options%ScaleEmis = .FALSE.

    ! NO
    IF ( Inst%IDTNO > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( HcoState, FLUXNO, Inst%IDTNO, &
                         RC,       ExtNr=Inst%ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'HCO_EmisAdd error: FLUXNO', RC )
          RETURN
       ENDIF
    ENDIF

    ! NO2
    IF ( Inst%IDTNO2 > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( HcoState, FLUXNO2, Inst%IDTNO2, &
                         RC,       ExtNr=Inst%ExtNr )
    ENDIF

    ! HNO3
    IF ( Inst%IDTHNO3 > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( HcoState, FLUXHNO3, Inst%IDTHNO3, &
                         RC,       ExtNr=Inst%ExtNr )
    ENDIF

    ! O3
    IF ( Inst%IDTO3 > 0 ) THEN

       ! Add flux to emission array (kg/m2/s)
       CALL HCO_EmisAdd( HcoState, FLUXO3, Inst%IDTO3, &
                         RC,       ExtNr=Inst%ExtNr )
    ENDIF


    ! Eventually update manual diagnostics
    IF ( DoDiagn ) THEN
       DiagnName =  'PARANOX_NOXFRAC_REMAINING'
       Arr2D     => DIAGN(:,:,1)
       CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                          cName=TRIM(DiagnName), Array2D=Arr2D, RC=RC)
       IF ( RC /= HCO_SUCCESS ) RETURN
       Arr2D => NULL()

       DiagnName =  'PARANOX_OPE'
       Arr2D     => DIAGN(:,:,2)
       CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                          cName=TRIM(DiagnName), Array2D=Arr2D, RC=RC)
       IF ( RC /= HCO_SUCCESS ) RETURN
       Arr2D => NULL()

       DiagnName =  'PARANOX_O3_PRODUCTION'
       Arr2D     => DIAGN(:,:,3)
       CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                          cName=TRIM(DiagnName), Array2D=Arr2D, RC=RC)
       IF ( RC /= HCO_SUCCESS ) RETURN
       Arr2D => NULL()

       DiagnName =  'PARANOX_TOTAL_SHIPNOX'
       Arr2D     => DIAGN(:,:,4)
       CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                          cName=TRIM(DiagnName), Array2D=Arr2D, RC=RC)
       IF ( RC /= HCO_SUCCESS ) RETURN
       Arr2D => NULL()

       DiagnName =  'PARANOX_NO_PRODUCTION'
       Arr2D     => DIAGN(:,:,5)
       CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                          cName=TRIM(DiagnName), Array2D=Arr2D, RC=RC)
       IF ( RC /= HCO_SUCCESS ) RETURN
       Arr2D => NULL()
    ENDIF

    ! Reset option ScaleEmis to default value
    HcoState%Options%ScaleEmis = DefScaleEmis

    ! Return w/ success
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE Evolve_Plume
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_ParaNOx_Init
!
! !DESCRIPTION: Subroutine HcoX\_ParaNOx\_Init initializes the HEMCO
! PARANOX extension.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE HCOX_ParaNOx_Init( HcoState, ExtName, ExtState, RC )
!
! !USES:
!
   USE HCO_Chartools_Mod, ONLY : HCO_CharParse
   USE HCO_State_MOD,     ONLY : HCO_GetHcoID
   USE HCO_State_MOD,     ONLY : HCO_GetExtHcoID
   USE HCO_ExtList_Mod,   ONLY : GetExtNr
   USE HCO_ExtList_Mod,   ONLY : GetExtOpt
   USE HCO_Restart_Mod,   ONLY : HCO_RestartDefine
!   USE ParaNOx_Util_Mod,  ONLY : Read_ParaNOx_LUT
!
! !INPUT PARAMETERS:
!
   CHARACTER(LEN=*), INTENT(IN   )  :: ExtName       ! Extension name
   TYPE(Ext_State),  POINTER        :: ExtState      ! Module options
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(HCO_State),  POINTER        :: HcoState      ! HEMCO state object
   INTEGER,          INTENT(INOUT)  :: RC            ! Success or failure?
!
! !REVISION HISTORY:
!  06 Aug 2013 - C. Keller   - Initial Version
!  06 Jun 2014 - R. Yantosca - Cosmetic changes in ProTex Headers
!  06 Jun 2014 - R. Yantosca - Now indented using F90 free-format
!  13 Aug 2014 - R. Yantosca - Now read the PARANOX look-up tables here
!  14 Aug 2014 - R. Yantosca - Minor fix, read the PARANOX look-up tables
!                              after displaying text about PARANOX extension
!  16 Oct 2014 - C. Keller   - Added error check after READ_PARANOX_LUT
!  17 Oct 2014 - C. Keller   - Now parse input files via HCO_CharParse
!  17 Apr 2015 - C. Keller   - Now assign PARANOX_SUNCOS1 to SC5(:,:,1), etc.
!  25 May 2015 - C. Keller   - Now calculate SC5 via HCO_GetSUNCOS
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER                        :: ExtNr, I, NN, tmpID, nSpc
   INTEGER,           ALLOCATABLE :: HcoIDs(:)
   CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
   CHARACTER(LEN=31)              :: Dummy
   CHARACTER(LEN=31)              :: DiagnName
   CHARACTER(LEN=255)             :: MSG, LOC
   CHARACTER(LEN= 1)              :: CHAR1
   TYPE(MyInst), POINTER          :: Inst

   !========================================================================
   ! HCOX_PARANOX_INIT begins here!
   !========================================================================

   ! Assume success
   RC = HCO_SUCCESS

   ! Extension Nr.
   ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
   IF ( ExtNr <= 0 ) RETURN

   ! Enter
   CALL HCO_ENTER( HcoState%Config%Err,                                      &
                   'HCOX_ParaNOx_Init (hcox_paranox_mod.F90)', RC           )
   IF ( RC /= HCO_SUCCESS ) RETURN

   ! Create local instance
   Inst => NULL()
   CALL InstCreate( ExtNr, ExtState%ParaNOx, Inst, RC                       )
   IF ( RC /= HCO_SUCCESS ) THEN
      CALL HCO_ERROR( HcoState%Config%Err,                                   &
                      'Cannot create ParaNOx instance', RC                  )
      RETURN
   ENDIF

   !========================================================================
   ! Skip the following for GEOS-Chem dry-run or HEMCO-standalone dry-run
   !========================================================================
   IF ( .not. HcoState%Options%IsDryRun ) THEN

      !---------------------------------------------------------------------
      ! Initialize fields of Inst object for safety's sake (bmy, 10/17/18)
      !---------------------------------------------------------------------
      Inst%IDTNO         = -1
      Inst%IDTNO2        = -1
      Inst%IDTO3         = -1
      Inst%IDTHNO3       = -1
      Inst%MW_O3         =  0.0d0
      Inst%MW_NO         =  0.0d0
      Inst%MW_NO2        =  0.0d0
      Inst%MW_HNO3       =  0.0d0
      Inst%MW_AIR        =  0.0d0
      Inst%Tlev          =  0.0e0
      Inst%JNO2lev       =  0.0e0
      Inst%O3lev         =  0.0e0
      Inst%SEA0lev       =  0.0e0
      Inst%SEA5lev       =  0.0e0
      Inst%JRATIOlev     =  0.0e0
      Inst%NOXlev        =  0.0e0
      Inst%WSlev         =  0.0e0
      Inst%ShipNO        => NULL()
      Inst%SC5           => NULL()
      Inst%DEPO3         => NULL()
      Inst%DEPHNO3       => NULL()
      Inst%FRACNOX_LUT02 => NULL()
      Inst%FRACNOX_LUT06 => NULL()
      Inst%FRACNOX_LUT10 => NULL()
      Inst%FRACNOX_LUT14 => NULL()
      Inst%FRACNOX_LUT18 => NULL()
      Inst%OPE_LUT02     => NULL()
      Inst%OPE_LUT06     => NULL()
      Inst%OPE_LUT10     => NULL()
      Inst%OPE_LUT14     => NULL()
      Inst%OPE_LUT18     => NULL()
      Inst%MOE_LUT02     => NULL()
      Inst%MOE_LUT06     => NULL()
      Inst%MOE_LUT10     => NULL()
      Inst%MOE_LUT14     => NULL()
      Inst%MOE_LUT18     => NULL()
      Inst%DNOX_LUT02    => NULL()
      Inst%DNOX_LUT06    => NULL()
      Inst%DNOX_LUT10    => NULL()
      Inst%DNOX_LUT14    => NULL()
      Inst%DNOX_LUT18    => NULL()

      !------------------------------------------------------------------------
      ! Get species IDs
      !------------------------------------------------------------------------

      ! Get HEMCO species IDs
      CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Check for NO, NO2, O3, and HNO3
      DO I = 1, nSpc
         SELECT CASE ( TRIM(SpcNames(I)) )
            CASE ( "NO" )
               Inst%IDTNO = HcoIDs(I)
            CASE ( "NO2" )
               Inst%IDTNO2 = HcoIDs(I)
            CASE ( "O3" )
               Inst%IDTO3 = HcoIDs(I)
            CASE ( "HNO3" )
               Inst%IDTHNO3 = HcoIDs(I)
            CASE DEFAULT
               ! leave empty
         END SELECT
      ENDDO

      ! Get MW of all species. If species not set in the configuration
      ! file (e.g. they are not being used by PARANOx), determine MW from
      ! default values.

      ! O3
      IF ( Inst%IDTO3 <= 0 ) THEN
         tmpID = HCO_GetHcoID('O3', HcoState )
         MSG = 'O3 not produced/removed in PARANOX'
         CALL HCO_WARNING(HcoState%Config%Err, MSG, RC )
      ELSE
         tmpID = Inst%IDTO3
      ENDIF
      IF ( tmpID > 0 ) THEN
         Inst%MW_O3 = HcoState%Spc(tmpID)%MW_g
      ELSE
         MSG = 'Use default O3 molecular weight of 48g/mol'
         CALL HCO_WARNING(HcoState%Config%Err, MSG, RC )
         Inst%MW_O3 = 48.0_dp
      ENDIF

      ! NO
      IF ( Inst%IDTNO <= 0 ) THEN
         tmpID = HCO_GetHcoID('NO', HcoState )
         MSG = 'NO not produced in PARANOX'
         CALL HCO_WARNING(HcoState%Config%Err, MSG, RC )
      ELSE
         tmpID = Inst%IDTNO
      ENDIF
      IF ( tmpID > 0 ) THEN
         Inst%MW_NO = HcoState%Spc(tmpID)%MW_g
      ELSE
         MSG = 'Use default NO molecular weight of 30g/mol'
         CALL HCO_WARNING(HcoState%Config%Err, MSG, RC )
         Inst%MW_NO = 30.0_dp
      ENDIF

      ! NO2
      IF ( Inst%IDTNO2 <= 0 ) THEN
         tmpID = HCO_GetHcoID('NO2', HcoState )
         MSG = 'NO2 not produced in PARANOX'
         CALL HCO_WARNING(HcoState%Config%Err, MSG, RC )
      ELSE
         tmpID = Inst%IDTNO2
      ENDIF
      IF ( tmpID > 0 ) THEN
         Inst%MW_NO2 = HcoState%Spc(tmpID)%MW_g
      ELSE
         MSG = 'Use default NO2 molecular weight of 46g/mol'
         CALL HCO_WARNING(HcoState%Config%Err, MSG, RC )
         Inst%MW_NO2 = 46.0_dp
      ENDIF

      ! HNO3
      IF ( Inst%IDTHNO3 <= 0 ) THEN
         tmpID = HCO_GetHcoID('HNO3', HcoState )
         MSG = 'HNO3 not produced/removed in PARANOX'
         CALL HCO_WARNING(HcoState%Config%Err, MSG, RC )
      ELSE
         tmpID = Inst%IDTHNO3
      ENDIF
      IF ( tmpID > 0 ) THEN
         Inst%MW_HNO3 = HcoState%Spc(tmpID)%MW_g
      ELSE
         MSG = 'Use default HNO3 molecular weight of 63g/mol'
         CALL HCO_WARNING(HcoState%Config%Err, MSG, RC )
         Inst%MW_HNO3 = 63.0_dp
      ENDIF

      ! Verbose mode
      IF ( HcoState%amIRoot ) THEN
         MSG = 'Use ParaNOx ship emissions (extension module)'
         CALL HCO_MSG(HcoState%Config%Err,MSG, SEP1='-' )
         MSG = '    - Use the following species: (MW, emitted as HEMCO ID) '
         CALL HCO_MSG(HcoState%Config%Err,MSG )
         WRITE(MSG,"(a,F5.2,I5)") '     NO  : ', Inst%MW_NO, Inst%IDTNO
         CALL HCO_MSG(HcoState%Config%Err,MSG)
         WRITE(MSG,"(a,F5.2,I5)") '     NO2 : ', Inst%MW_NO2, Inst%IDTNO2
         CALL HCO_MSG(HcoState%Config%Err,MSG)
         WRITE(MSG,"(a,F5.2,I5)") '     O3  : ', Inst%MW_O3, Inst%IDTO3
         CALL HCO_MSG(HcoState%Config%Err,MSG)
         WRITE(MSG,"(a,F5.2,I5)") '     HNO3: ', Inst%MW_HNO3, Inst%IDTHNO3
         CALL HCO_MSG(HcoState%Config%Err,MSG)
      ENDIF

      !--------------------------------
      ! Allocate module arrays
      !--------------------------------

      ! FNOX
      ALLOCATE( Inst%FRACNOX_LUT02(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'FRACNOX_LUT02', RC )
         RETURN
      ENDIF
      Inst%FRACNOX_LUT02 = 0.0_sp

      ALLOCATE( Inst%FRACNOX_LUT06(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'FRACNOX_LUT06', RC )
         RETURN
      ENDIF
      Inst%FRACNOX_LUT06 = 0.0_sp

      ALLOCATE( Inst%FRACNOX_LUT10(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'FRACNOX_LUT10', RC )
         RETURN
      ENDIF
      Inst%FRACNOX_LUT10 = 0.0_sp

      ALLOCATE( Inst%FRACNOX_LUT14(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'FRACNOX_LUT014', RC )
         RETURN
      ENDIF
      Inst%FRACNOX_LUT14 = 0.0_sp

      ALLOCATE( Inst%FRACNOX_LUT18(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'FRACNOX_LUT18', RC )
         RETURN
      ENDIF
      Inst%FRACNOX_LUT18 = 0.0_sp

      ! OPE
      ALLOCATE( Inst%OPE_LUT02(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'OPE_LUT02', RC )
         RETURN
      ENDIF
      Inst%OPE_LUT02 = 0.0_sp

      ALLOCATE( Inst%OPE_LUT06(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'OPE_LUT06', RC )
         RETURN
      ENDIF
      Inst%OPE_LUT06 = 0.0_sp

      ALLOCATE( Inst%OPE_LUT10(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= 0 ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'OPE_LUT10', RC )
         RETURN
      ENDIF
      Inst%OPE_LUT10 = 0.0_sp

      ALLOCATE( Inst%OPE_LUT14(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= 0 ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'OPE_LUT014', RC )
         RETURN
      ENDIF
      Inst%OPE_LUT14 = 0.0_sp

      ALLOCATE( Inst%OPE_LUT18(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'OPE_LUT18', RC )
         RETURN
      ENDIF
      Inst%OPE_LUT18 = 0.0_sp

      ! MOE
      ALLOCATE( Inst%MOE_LUT02(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'MOE_LUT02', RC )
         RETURN
      ENDIF
      Inst%MOE_LUT02 = 0.0_sp

      ALLOCATE( Inst%MOE_LUT06(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'MOE_LUT06', RC )
         RETURN
      ENDIF
      Inst%MOE_LUT06 = 0.0_sp

      ALLOCATE( Inst%MOE_LUT10(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'MOE_LUT10', RC )
         RETURN
      ENDIF
      Inst%MOE_LUT10 = 0.0_sp

      ALLOCATE( Inst%MOE_LUT14(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'MOE_LUT014', RC )
         RETURN
      ENDIF
      Inst%MOE_LUT14 = 0.0_sp

      ALLOCATE( Inst%MOE_LUT18(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'MOE_LUT18', RC )
         RETURN
      ENDIF
      Inst%MOE_LUT18 = 0.0_sp

      ! DNOx
      ALLOCATE( Inst%DNOx_LUT02(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'DNOx_LUT02', RC )
         RETURN
      ENDIF
      Inst%DNOx_LUT02 = 0.0_sp

      ALLOCATE( Inst%DNOx_LUT06(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'DNOx_LUT06', RC )
         RETURN
      ENDIF
      Inst%DNOx_LUT06 = 0.0_sp

      ALLOCATE( Inst%DNOx_LUT10(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'DNOx_LUT10', RC )
         RETURN
      ENDIF
      Inst%DNOx_LUT10 = 0.0_sp

      ALLOCATE( Inst%DNOx_LUT14(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'DNOx_LUT014', RC )
         RETURN
      ENDIF
      Inst%DNOx_LUT14 = 0.0_sp

      ALLOCATE( Inst%DNOx_LUT18(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'DNOx_LUT18', RC )
         RETURN
      ENDIF
      Inst%DNOx_LUT18 = 0.0_sp

      ALLOCATE(Inst%DEPO3  (HcoState%NX,HcoState%NY),        &
               Inst%DEPHNO3(HcoState%NX,HcoState%NY), STAT=RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'Deposition arrays', RC )
         RETURN
      ENDIF
      Inst%DEPO3   = 0.0_sp
      Inst%DEPHNO3 = 0.0_sp

   !   ! O3 loss and HNO3 deposition
   !   ALLOCATE( Inst%SHIPO3LOSS(HcoState%NX,HcoState%NY), STAT=RC )
   !   IF ( RC /= HCO_SUCCESS ) THEN
   !      CALL HCO_ERROR ( HcoState%Config%Err, 'SHIPO3LOSS', RC )
   !      RETURN
   !   ENDIF
   !   Inst%SHIPO3LOSS = 0d0

   !   ALLOCATE( Inst%SHIPHNO3DEP(HcoState%NX,HcoState%NY), STAT=RC )
   !   IF ( RC /= HCO_SUCCESS ) THEN
   !        CALL HCO_ERROR ( HcoState%Config%Err, 'SHIPHNO3DEP', RC ); RETURN
   !   ENDIF
   !   Inst%SHIPHNO3DEP = 0d0

   ENDIF

   !========================================================================
   ! Initialize the PARANOX look-up tables
   !========================================================================

   ! LUT data directory
   CALL GetExtOpt( HcoState%Config, Inst%ExtNr, 'LUT source dir', &
                   OptValChar=Inst%LutDir, RC=RC)
   IF ( RC /= HCO_SUCCESS ) RETURN

   ! Call HEMCO parser to replace tokens such as $ROOT, $MET, or $RES.
   ! There shouldn't be any date token in there ($YYYY, etc.), so just
   ! provide some dummy variables here
   CALL HCO_CharParse( HcoState%Config, Inst%LutDir, -999, -1, -1, -1, -1, RC )
   IF ( RC /= HCO_SUCCESS ) RETURN

   ! Data format: ncdf (default) or txt
   Inst%IsNc = .TRUE.
   CALL GetExtOpt( HcoState%Config, Inst%ExtNr, 'LUT data format', &
                   OptValChar=Dummy, RC=RC)
   IF ( RC /= HCO_SUCCESS ) RETURN
   IF ( TRIM(Dummy) == 'txt' ) Inst%IsNc = .FALSE.

   ! Read PARANOX look-up tables from disk. This can be netCDF or txt
   ! format, as determined above.
   !
   ! NOTE: For the GEOS-Chem dry-run or HEMCO-standalone dry-run,
   ! these routines will print file paths to the dry-run log file,
   ! but will not actually read any data.
   IF ( Inst%IsNc ) THEN
      CALL READ_PARANOX_LUT_NC( HcoState, Inst, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
   ELSE
      CALL READ_PARANOX_LUT_TXT( HcoState, Inst, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
   ENDIF

   !========================================================================
   ! Exit if this is a GEOS-Chem dry-run or HEMCO-standalone dry-run
   !========================================================================
   IF ( HcoState%Options%IsDryRun ) THEN
      Inst => NULL()
      CALL HCO_LEAVE( HcoState%Config%Err,RC )
      RETURN
   ENDIF

   !========================================================================
   ! Continue initializing PARANOX for regular simulations
   !========================================================================

   !------------------------------------------------------------------------
   ! Set other module variables
   !------------------------------------------------------------------------
   ALLOCATE ( Inst%ShipNO(HcoState%NX,HcoState%NY,HcoState%NZ), STAT=RC )
   IF ( RC /= HCO_SUCCESS ) THEN
      CALL HCO_ERROR ( HcoState%Config%Err, 'ShipNO', RC )
      RETURN
   ENDIF
   Inst%ShipNO = 0.0_hp

   ! Allocate variables for SunCosMid from 5 hours ago.
   ALLOCATE ( Inst%SC5(HcoState%NX,HcoState%NY), STAT=RC )
   IF ( RC /= HCO_SUCCESS ) THEN
      CALL HCO_ERROR ( HcoState%Config%Err, 'SC5', RC )
      RETURN
   ENDIF
   Inst%SC5 = 0.0_hp

   ! Prompt warning if chemistry time step is more than 60 mins
   IF ( HcoState%TS_CHEM > 3600.0_hp ) THEN
      IF ( HcoState%amIRoot ) THEN
         MSG = ' Cannot properly store SUNCOS values ' // &
               ' because chemistry time step is more than 60 mins!'
         CALL HCO_WARNING(HcoState%Config%Err, MSG, RC )
      ENDIF
   ENDIF

   ! Molecular weight of AIR
   Inst%MW_AIR = HcoState%Phys%AIRMW

   !------------------------------------------------------------------------
   ! Define PARANOX diagnostics for the O3 and HNO3 deposition fluxes (in
   ! kg/m2/s).
   !------------------------------------------------------------------------
   CALL Diagn_Create ( HcoState,                                             &
                       cName    = 'PARANOX_O3_DEPOSITION_FLUX',              &
                       Trgt2D   = Inst%DEPO3,                                &
                       SpaceDim = 2,                                         &
                       OutUnit  = 'kg/m2/s',                                 &
                       COL = HcoState%Diagn%HcoDiagnIDManual,                &
                       RC       = RC                                        )
   IF ( RC /= HCO_SUCCESS ) RETURN

   CALL Diagn_Create ( HcoState,                                             &
                       cName    = 'PARANOX_HNO3_DEPOSITION_FLUX',            &
                       Trgt2D   = Inst%DEPHNO3,                              &
                       SpaceDim = 2,                                         &
                       OutUnit  = 'kg/m2/s',                                 &
                       COL = HcoState%Diagn%HcoDiagnIDManual,                &
                       RC       = RC                                         )
   IF ( RC /= HCO_SUCCESS ) RETURN

   !------------------------------------------------------------------------
   ! Set HEMCO extension variables
   !------------------------------------------------------------------------

   ! Met. data required by module
   ExtState%O3%DoUse          = .TRUE.
   ExtState%NO2%DoUse         = .TRUE.
   ExtState%NO%DoUse          = .TRUE.
   ExtState%AIR%DoUse         = .TRUE.
   ExtState%AIRVOL%DoUse      = .TRUE.
   ExtState%SUNCOS%DoUse      = .TRUE.
   ExtState%T2M%DoUse         = .TRUE.
   ExtState%U10M%DoUse        = .TRUE.
   ExtState%V10M%DoUse        = .TRUE.
   ExtState%FRAC_OF_PBL%DoUse = .TRUE.
   IF ( Inst%IDTHNO3 > 0 ) THEN
      ExtState%HNO3%DoUse     = .TRUE.
   ENDIF
   ExtState%JNO2%DoUse        = .TRUE.
   ExtState%JOH%DoUse         = .TRUE.

   !------------------------------------------------------------------------
   ! Leave w/ success
   !------------------------------------------------------------------------
   IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  )
   IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
   Inst => NULL()
   CALL HCO_LEAVE( HcoState%Config%Err,RC )

 END SUBROUTINE HCOX_ParaNOx_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_ParaNOx_Final
!
! !DESCRIPTION: Subroutine HcoX\_ParaNox\_Final finalizes the HEMCO
! PARANOX extension.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE HCOX_ParaNOx_Final( HcoState, ExtState, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState      ! HEMCO State obj
    TYPE(Ext_State), POINTER        :: ExtState      ! Module options
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  06 Aug 2013 - C. Keller - Initial Version
!  06 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  06 Jun 2014 - R. Yantosca - Now indended with F90 free-format
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!

   !=================================================================
   ! HCOX_PARANOX_FINAL begins here!
   !=================================================================
    CALL InstRemove( ExtState%ParaNOx )

   RC = HCO_SUCCESS

 END SUBROUTINE HCOX_ParaNOx_Final
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_paranox_lut_nc
!
! !DESCRIPTION: Subroutine READ\_PARANOX\_LUT\_NC reads look-up tables in
!  netCDF format for use in the PARANOX ship plume model (G.C.M. Vinken)
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE READ_PARANOX_LUT_NC ( HcoState, Inst, RC )
!
! !USES:
!
! !INPUT ARGUMENTS:
!
   TYPE(HCO_State), POINTER     :: HcoState    ! HEMCO State object
   TYPE(MyInst),    POINTER     :: Inst
!
! !INPUT/OUTPUT ARGUMENTS:
!
   INTEGER, INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  06 Feb 2012 - M. Payer    - Initial version modified from code provided by
!                              G.C.M. Vinken
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  03 Jun 2013 - C. Holmes   - Rewritten to include wind speed in the look-up
!                              table and to take input from netCDF
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
   INTEGER             :: IOS
   CHARACTER(LEN=255)  :: FILENAME
   CHARACTER(LEN=255)  :: MSG
   INTEGER             :: fID

   !=================================================================
   ! READ_PARANOX_LUT_NC begins here
   !=================================================================

   ! NetCDF reading of PARANOX LUT not supported in ESMF environment
#if defined(ESMF_)
   MSG = 'In ESMF, cannot read PARANOX look-up-table in netCDF ' // &
         'format. Please set `LUT data format` to `txt` in the ' // &
         'HEMCO configuration file.'
   CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, &
           THISLOC = 'READ_PARANOX_LUT_NC (hcox_paranox_mod.F90)' )
   RETURN
#else

   ! Clear FILENAME
   FILENAME = ''

   ! FILENAME format string
 101  FORMAT( a, '/ship_plume_lut_', I2.2, 'ms.nc'  )

   ! Wind speed levels correspond to the files we will read below
   Inst%WSlev = (/ 2.0e0, 6.0e0, 10.0e0, 14.0e0, 18.0e0 /)

   ! Read 2m/s LUT
   WRITE( FILENAME, 101 ) TRIM(Inst%LutDir), 2
   CALL READ_LUT_NCFILE( HcoState, TRIM( FILENAME ),                         &
        Inst%FRACNOX_LUT02, Inst%DNOx_LUT02, Inst%OPE_LUT02, Inst%MOE_LUT02, &
        Inst%Tlev, Inst%JNO2lev, Inst%O3lev, Inst%SEA0lev, Inst%SEA5lev,     &
        Inst%JRATIOlev, Inst%NOXlev, RC=RC)

   ! Read 6 m/s LUT
   WRITE( FILENAME, 101 ) TRIM(Inst%LutDir), 6
   CALL READ_LUT_NCFILE( HcoState, TRIM( FILENAME ),                         &
        Inst%FRACNOX_LUT06, Inst%DNOx_LUT06, Inst%OPE_LUT06,                 &
        Inst%MOE_LUT06, RC=RC )

   ! Read 10 m/s LUT
   WRITE( FILENAME, 101 ) TRIM(Inst%LutDir), 10
   CALL READ_LUT_NCFILE( HcoState, TRIM( FILENAME ),                         &
        Inst%FRACNOX_LUT10, Inst%DNOx_LUT10, Inst%OPE_LUT10,                 &
        Inst%MOE_LUT10, RC=RC )

   ! Read 14 m/s LUT
   WRITE( FILENAME, 101 ) TRIM(Inst%LutDir), 14
   CALL READ_LUT_NCFILE( HcoState, TRIM( FILENAME ),                         &
        Inst%FRACNOX_LUT14, Inst%DNOx_LUT14, Inst%OPE_LUT14,                 &
        Inst%MOE_LUT14, RC=RC )

   ! Read 18 m/s LUT
   WRITE( FILENAME, 101 ) TRIM(Inst%LutDir), 18
   CALL READ_LUT_NCFILE( HcoState, TRIM( FILENAME ),                         &
        Inst%FRACNOX_LUT18, Inst%DNOx_LUT18, Inst%OPE_LUT18,                 &
        Inst%MOE_LUT18, RC=RC )

!   ! To write into txt-file, uncomment the following lines
!   FILENAME = TRIM(LutDir)//'/ship_plume_lut_02ms.txt'
!   CALL WRITE_LUT_TXTFILE( HcoState, TRIM( FILENAME ), &
!        FRACNOX_LUT02, DNOx_LUT02, OPE_LUT02, MOE_LUT02, RC, &
!        Tlev, JNO2lev, O3lev, SEA0lev, SEA5lev, JRATIOlev, NOXlev )
!   IF ( RC /= HCO_SUCCESS ) RETURN
!
!   FILENAME = TRIM(LutDir)//'/ship_plume_lut_06ms.txt'
!   CALL WRITE_LUT_TXTFILE( HcoState, TRIM( FILENAME ), &
!        FRACNOX_LUT06, DNOx_LUT06, OPE_LUT06, MOE_LUT06, RC, &
!        Tlev, JNO2lev, O3lev, SEA0lev, SEA5lev, JRATIOlev, NOXlev )
!   IF ( RC /= HCO_SUCCESS ) RETURN
!
!   FILENAME = TRIM(LutDir)//'/ship_plume_lut_10ms.txt'
!   CALL WRITE_LUT_TXTFILE( HcoState, TRIM( FILENAME ), &
!        FRACNOX_LUT10, DNOx_LUT10, OPE_LUT10, MOE_LUT10, RC, &
!        Tlev, JNO2lev, O3lev, SEA0lev, SEA5lev, JRATIOlev, NOXlev )
!   IF ( RC /= HCO_SUCCESS ) RETURN
!
!   FILENAME = TRIM(LutDir)//'/ship_plume_lut_14ms.txt'
!   CALL WRITE_LUT_TXTFILE( HcoState, TRIM( FILENAME ), &
!        FRACNOX_LUT14, DNOx_LUT14, OPE_LUT14, MOE_LUT14, RC, &
!        Tlev, JNO2lev, O3lev, SEA0lev, SEA5lev, JRATIOlev, NOXlev )
!   IF ( RC /= HCO_SUCCESS ) RETURN
!
!   FILENAME = TRIM(LutDir)//'/ship_plume_lut_14ms.txt'
!   CALL WRITE_LUT_TXTFILE( HcoState, TRIM( FILENAME ), &
!        FRACNOX_LUT18, DNOx_LUT18, OPE_LUT18, MOE_LUT18, RC, &
!        Tlev, JNO2lev, O3lev, SEA0lev, SEA5lev, JRATIOlev, NOXlev )
!   IF ( RC /= HCO_SUCCESS ) RETURN

   ! Return w/ success
   RC = HCO_SUCCESS
#endif

 END SUBROUTINE READ_PARANOX_LUT_NC
!EOC
#if !defined(ESMF_)
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_lut_ncfile
!
! !DESCRIPTION: Subroutine READ\_LUT\_NCFILE reads look up tables for use in
!  the PARANOX ship plume model (C. Holmes)
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE READ_LUT_NCFILE( HcoState, FILENAME, FNOX,                       &
                             DNOx,      OPE,      MOE,      T,               &
                             JNO2,      O3,       SEA0,     SEA5,            &
                             JRATIO,    NOX,      RC                        )
!
! !USES:
!
   ! Modules for netCDF read
   USE m_netcdf_io_open
   USE m_netcdf_io_get_dimlen
   USE m_netcdf_io_read
   USE m_netcdf_io_readattr
   USE m_netcdf_io_close

#  include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
   TYPE(HCO_State), POINTER     :: HcoState    ! HEMCO State object
   CHARACTER(LEN=*),INTENT(IN)  :: FILENAME
!
! !OUTPUT PARAMETERS:
!
   REAL*4,  INTENT(OUT), DIMENSION(:,:,:,:,:,:,:) :: FNOX,OPE,MOE,DNOx
   REAL*4,  INTENT(OUT), OPTIONAL :: T(:),      JNO2(:), O3(:)
   REAL*4,  INTENT(OUT), OPTIONAL :: SEA0(:),   SEA5(:)
   REAL*4,  INTENT(OUT), OPTIONAL :: JRATIO(:), NOX(:)
   INTEGER, INTENT(OUT), OPTIONAL :: RC
!
! !REVISION HISTORY:
!  06 Feb 2012 - M. Payer    - Initial version modified from code provided by
!                              G.C.M. Vinken
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  03 Jun 2013 - C. Holmes   - Rewritten to include wind speed in the look-up
!                              table and to take input from netCDF
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

   ! Scalars
   LOGICAL             :: FileExists
   INTEGER             :: AS, IOS
   INTEGER             :: fID, HMRC

   ! arrays
   INTEGER             :: st1d(1), ct1d(1)
   INTEGER             :: st7d(7), ct7d(7)

   CHARACTER(LEN=255)  :: MSG,     FileMsg

   !=================================================================
   ! In dry-run mode, print file path to dryrun log and exit.
   ! Otherwise, print file path to the HEMCO log file and continue.
   !=================================================================

   ! Test if the file exists
   INQUIRE( FILE=TRIM( FileName ), EXIST=FileExists )

   ! Create a display string based on whether or not the file is found
   IF ( FileExists ) THEN
      FileMsg = 'HEMCO (PARANOX): Opening'
   ELSE
      FileMsg = 'HEMCO (PARANOX): REQUIRED FILE NOT FOUND'
   ENDIF

   ! Print file status to stdout and the HEMCO log
   IF ( HcoState%amIRoot ) THEN
      WRITE( 6,   300 ) TRIM( FileMsg ), TRIM( FileName )
      WRITE( MSG, 300 ) TRIM( FileMsg ), TRIM( FileName )
      CALL HCO_MSG( HcoState%Config%Err, MSG )
 300  FORMAT( a, ' ', a )
   ENDIF

   ! For dry-run simulations, return to calling program.
   ! For regular simulations, throw an error if we can't find the file.
   IF ( HcoState%Options%IsDryRun ) THEN
      RETURN
   ELSE
      IF ( .not. FileExists ) THEN
         WRITE( MSG, 300 ) TRIM( FileMsg ), TRIM( FileName )
         CALL HCO_ERROR(HcoState%Config%Err, MSG, HMRC )
         IF ( PRESENT( RC ) ) RC = HMRC
         RETURN
      ENDIF
   ENDIF

   !=================================================================
   ! READ_LUT_NCFILE begins here!
   !=================================================================

   ! Open file for reading
   CALL Ncop_Rd( fId, TRIM(FILENAME) )

   !-----------------------------------------------------------------
   ! Read reference values used to construct the Look-up table
   ! These are 1d values
   !-----------------------------------------------------------------

   ! Temperature, K
   IF ( PRESENT(T) ) THEN
      st1d = (/ 1  /)
      ct1d = (/ nT /)
      CALL NcRd( T, fId, 'T', st1d, ct1d )
   ENDIF

   ! J(NO2), 1/s
   IF ( PRESENT(JNO2) ) THEN
      st1d = (/ 1  /)
      ct1d = (/ nJ /)
      CALL NcRd( JNO2, fId, 'JNO2', st1d, ct1d )
   ENDIF

   ! Ambient O3, ppb
   IF ( PRESENT(O3) ) THEN
      st1d = (/ 1   /)
      ct1d = (/ nO3 /)
      CALL NcRd( O3, fId, 'O3', st1d, ct1d )
   ENDIF

   ! Solar elevation angle at emission time t=0, deg
   IF ( PRESENT(SEA0) ) THEN
      st1d = (/ 1    /)
      ct1d = (/ nSEA /)
      CALL NcRd( SEA0, fId, 'SEA0', st1d, ct1d )
   ENDIF

   ! Solar elevation angle at time t=5h, deg
   IF ( PRESENT(SEA5) ) THEN
      st1d = (/ 1    /)
      ct1d = (/ nSEA /)
      CALL NcRd( SEA5, fId, 'SEA5', st1d, ct1d )
   ENDIF

   ! J(OH) / J(NO2) ratio, s/s
   IF ( PRESENT(JRATIO) ) THEN
      st1d = (/ 1  /)
      ct1d = (/ nJ /)
      CALL NcRd( JRatio, fId, 'Jratio', st1d, ct1d )
   ENDIF

   ! Ambient NOx, ppt
   IF ( PRESENT(NOX) ) THEN
      st1d = (/ 1  /)
      ct1d = (/ nNOx /)
      CALL NcRd( NOX, fId, 'NOx', st1d, ct1d )
   ENDIF

   ! Define 7D variables used below
   st7d = (/ 1, 1, 1,  1,   1,   1, 1    /)
   ct7d = (/ nT,nJ,nO3,nSEA,nSEA,nJ,nNOx /)

   !-----------------------------------------------------------------
   ! Read look up table for fraction of NOx remaining for ship
   ! emissions after 5 h [unitless]
   !-----------------------------------------------------------------

   CALL NcRd( FNOx, fId, 'FNOx', st7d, ct7d )

  ! testing only
!   PRINT*, "binary_fracnox: ", Fnox(1:4,1,1,2,3,4,4)

   !-----------------------------------------------------------------
   ! Read look up table for 5-h integrated Ozone Production Efficiency
   ! for ship emissions [molec O3 produced / molec NOx lost]
   !-----------------------------------------------------------------

   CALL NcRd( OPE, fId, 'OPE', st7d, ct7d )

  ! testing only
!   PRINT*, "binary_intope: ", OPE(1:4,1,1,2,3,4,4)

   !-----------------------------------------------------------------
   ! Read look up table for 5-h integrated Methane Oxidation Efficiency
   ! for ship emissions [molec CH4 oxidized / molec NOx emitted]
   !-----------------------------------------------------------------

   CALL NcRd( MOE, fId, 'MOE', st7d, ct7d )

  ! testing only
!   PRINT*, "binary_intmoe: ", MOE(1:4,1,1,2,3,4,4)

   !-----------------------------------------------------------------
   ! Read look up table for 5-h integrated NOx deposition fraction
   ! for ship emissions [molec NOx deposited / molec NOx emitted]
   !-----------------------------------------------------------------

   CALL NcRd( DNOx, fId, 'DNOx', st7d, ct7d )

  ! testing only
!  PRINT*, "binary_depnox: ", DNOx(1:4,1,1,2,3,4,4)

   ! Close netCDF file
   CALL NcCl( fId )

 END SUBROUTINE READ_LUT_NCFILE
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_paranox_lut_txt
!
! !DESCRIPTION: Subroutine READ\_PARANOX\_LUT\_TXT reads look-up tables in
!  txt format for use in the PARANOX ship plume model (G.C.M. Vinken)
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE READ_PARANOX_LUT_TXT ( HcoState, Inst, RC )
!
! !USES:
!
! !INPUT ARGUMENTS:
!
   TYPE(HCO_State), POINTER       :: HcoState    ! HEMCO State object
   TYPE(MyInst),    POINTER       :: Inst
!
! !INPUT/OUTPUT ARGUMENTS:
!
   INTEGER, INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  05 Feb 2015 - C. Keller   - Initial version modified from code provided by
!                              G.C.M. Vinken and C. Holmes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
   INTEGER             :: IOS
   CHARACTER(LEN=255)  :: FILENAME
   CHARACTER(LEN=255)  :: MSG
   INTEGER             :: fID

   !=================================================================
   ! READ_PARANOX_LUT_TXT begins here
   !=================================================================

   ! Clear FILENAME
   FILENAME = ''

   ! FILENAME format string
 101  FORMAT( a, '/ship_plume_lut_', I2.2, 'ms.txt'  )

   ! Wind speed levels correspond to the files that we will read below
   Inst%WSlev = (/ 2.0e0, 6.0e0, 10.0e0, 14.0e0, 18.0e0 /)

   ! Read 2m/s LUT
   WRITE( FILENAME, 101 ) TRIM(Inst%LutDir), 2
   CALL READ_LUT_TXTFILE( HcoState, TRIM( FILENAME ), &
        Inst%FRACNOX_LUT02, Inst%DNOx_LUT02, Inst%OPE_LUT02, Inst%MOE_LUT02, RC, &
        Inst%Tlev, Inst%JNO2lev, Inst%O3lev, Inst%SEA0lev, Inst%SEA5lev, Inst%JRATIOlev, Inst%NOXlev )
   IF ( RC /= HCO_SUCCESS ) RETURN

   ! Read 6 m/s LUT
   WRITE( FILENAME, 101 ) TRIM(Inst%LutDir), 6
   CALL READ_LUT_TXTFILE( HcoState, TRIM( FILENAME ), &
        Inst%FRACNOX_LUT06, Inst%DNOx_LUT06, Inst%OPE_LUT06, Inst%MOE_LUT06, RC )
   IF ( RC /= HCO_SUCCESS ) RETURN

   ! Read 10 m/s LUT
   WRITE( FILENAME, 101 ) TRIM(Inst%LutDir), 10
   CALL READ_LUT_TXTFILE( HcoState, TRIM( FILENAME ), &
        Inst%FRACNOX_LUT10, Inst%DNOx_LUT10, Inst%OPE_LUT10, Inst%MOE_LUT10, RC )
   IF ( RC /= HCO_SUCCESS ) RETURN

   ! Read 14 m/s LUT
   WRITE( FILENAME, 101 ) TRIM(Inst%LutDir), 14
   CALL READ_LUT_TXTFILE( HcoState, TRIM( FILENAME ), &
        Inst%FRACNOX_LUT14, Inst%DNOx_LUT14, Inst%OPE_LUT14, Inst%MOE_LUT14, RC )
   IF ( RC /= HCO_SUCCESS ) RETURN

   ! Read 18 m/s LUT
   WRITE( FILENAME, 101 ) TRIM(Inst%LutDir), 18
   CALL READ_LUT_TXTFILE( HcoState, TRIM( FILENAME ), &
        Inst%FRACNOX_LUT18, Inst%DNOx_LUT18, Inst%OPE_LUT18, Inst%MOE_LUT18, RC )
   IF ( RC /= HCO_SUCCESS ) RETURN

   ! Return w/ success
   RC = HCO_SUCCESS

 END SUBROUTINE READ_PARANOX_LUT_TXT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_lut_txtfile
!
! !DESCRIPTION: Subroutine READ\_LUT\_TXTFILE reads look up tables for use in
!  the PARANOX ship plume model (C. Holmes)
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE READ_LUT_TXTFILE( HcoState, FILENAME, FNOX, DNOx, OPE, MOE, RC, &
                              T, JNO2, O3, SEA0, SEA5, JRATIO, NOX )
!
! !USES:
!
   USE inquireMod, ONLY : findFreeLUN
!
! !INPUT PARAMETERS:
!
   TYPE(HCO_State), POINTER     :: HcoState    ! HEMCO State object
   CHARACTER(LEN=*),INTENT(IN)  :: FILENAME
!
! !INPUT/OUTPUT PARAMETERS:
!
   INTEGER, INTENT(INOUT)       :: RC
!
! !OUTPUT PARAMETERS:
!
   REAL*4, INTENT(OUT), TARGET, DIMENSION(:,:,:,:,:,:,:) :: FNOX,OPE,MOE,DNOx
   REAL*4, INTENT(OUT), OPTIONAL :: T(:),      JNO2(:), O3(:)
   REAL*4, INTENT(OUT), OPTIONAL :: SEA0(:),   SEA5(:)
   REAL*4, INTENT(OUT), OPTIONAL :: JRATIO(:), NOX(:)
!
! !REVISION HISTORY:
!  05 Feb 2015 - C. Keller   - Initial version modified from code provided by
!                              G.C.M. Vinken and C. Holmes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Scalars
   LOGICAL             :: FileExists
   INTEGER             :: fId, IOS

!   INTEGER             :: I, I1, I2, I3, I4, I5, I6, I7
!   REAL*4, POINTER     :: TMPARR(:,:,:,:,:,:,:) => NULL()

   CHARACTER(LEN=255)  :: MSG, FileMsg
!
! !DEFINED PARAMETERS:
!
   CHARACTER(LEN=255), PARAMETER :: FMAT = "(E40.32)"
   CHARACTER(LEN=255), PARAMETER :: LOC  = &
                        'READ_LUT_TXTFILE (hcox_paranox_mod.F90)'

   !=================================================================
   ! In dry-run mode, print file path to dryrun log and exit.
   ! Otherwise, print file path to the HEMCO log file and continue.
   !=================================================================

   ! Test if the file exists
   INQUIRE ( FILE=TRIM( FileName ), EXIST=FileExists )

   ! Create a display string based on whether or not the file is found
   IF ( FileExists ) THEN
      FileMsg = 'HEMCO (PARANOX): Opening'
   ELSE
      FileMsg = 'HEMCO (PARANOX): REQUIRED FILE NOT FOUND'
   ENDIF

   ! Print file status to stdout and the HEMCO log file
   IF ( HcoState%amIRoot ) THEN
      WRITE( 6,   300 ) TRIM( FileMsg ), TRIM( FileName )
      WRITE( MSG, 300 ) TRIM( FileMsg ), TRIM( FileName )
      CALL HCO_MSG(HcoState%Config%Err,MSG)
 300  FORMAT( a, ' ', a )
   ENDIF

   ! For dry-run simulations, return to calling program.
   ! For regular simulations, throw an error if we can't find the file.
   IF ( HcoState%Options%IsDryRun ) THEN
      RETURN
   ELSE
      IF ( .not. FileExists ) THEN
         WRITE( MSG, 300 ) TRIM( FileMsg ), TRIM( FileName )
         CALL HCO_ERROR(HcoState%Config%Err, MSG, RC )
         RETURN
      ENDIF
   ENDIF

   !=================================================================
   ! READ_LUT_TXTFILE begins here
   !=================================================================

   ! Find a free file LUN
   fId = findFreeLUN()

   ! Open file for reading
   OPEN ( fID, FILE=TRIM(FILENAME), FORM="FORMATTED", IOSTAT=IOS )
   IF ( IOS /= 0 ) THEN
      CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile:1', RC, THISLOC=LOC )
      RETURN
   ENDIF

   ! Read FNOx
   READ( fId, FMT=FMAT, IOSTAT=IOS ) FNOx
   IF ( IOS /= 0 ) THEN
      CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: FNOx', RC, THISLOC=LOC )
      RETURN
   ENDIF

   ! Read OPE
   READ( fId, FMT=FMAT, IOSTAT=IOS ) OPE
   IF ( IOS /= 0 ) THEN
      CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: OPE', RC, THISLOC=LOC )
      RETURN
   ENDIF

   ! Read MOE
   READ( fId, FMT=FMAT, IOSTAT=IOS ) MOE
   IF ( IOS /= 0 ) THEN
      CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: MOE', RC, THISLOC=LOC )
      RETURN
   ENDIF

   ! Read DNOx
   READ( fId, FMT=FMAT, IOSTAT=IOS ) DNOx
   IF ( IOS /= 0 ) THEN
      CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: DNOx', RC, THISLOC=LOC )
      RETURN
   ENDIF

   ! Read optional values
   IF ( PRESENT(T) ) THEN
      READ( fId, FMT=FMAT, IOSTAT=IOS ) T
      IF ( IOS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: T', RC, THISLOC=LOC )
         RETURN
      ENDIF
   ENDIF

   IF ( PRESENT(JNO2) ) THEN
      READ( fId, FMT=FMAT, IOSTAT=IOS ) JNO2
      IF ( IOS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: JNO2', RC, THISLOC=LOC )
         RETURN
      ENDIF
   ENDIF

   IF ( PRESENT(O3) ) THEN
      READ( fId, FMT=FMAT, IOSTAT=IOS ) O3
      IF ( IOS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: O3', RC, THISLOC=LOC )
         RETURN
      ENDIF
   ENDIF

   IF ( PRESENT(SEA0) ) THEN
      READ( fId, FMT=FMAT, IOSTAT=IOS ) SEA0
      IF ( IOS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: SEA0', RC, THISLOC=LOC )
         RETURN
      ENDIF
   ENDIF

   IF ( PRESENT(SEA5) ) THEN
      READ( fId, FMT=FMAT, IOSTAT=IOS ) SEA5
      IF ( IOS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: SEA5', RC, THISLOC=LOC )
         RETURN
      ENDIF
   ENDIF

   IF ( PRESENT(JRATIO) ) THEN
      READ( fId, FMT=FMAT, IOSTAT=IOS ) JRATIO
      IF ( IOS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: JRATIO', RC, THISLOC=LOC )
         RETURN
      ENDIF
   ENDIF

   IF ( PRESENT(NOX) ) THEN
      READ( fId, FMT=FMAT, IOSTAT=IOS ) NOX
      IF ( IOS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: NOX', RC, THISLOC=LOC )
         RETURN
      ENDIF
   ENDIF

!   ! Read
!   DO I = 1,4
!
!      ! Set pointer to output array
!      SELECT CASE ( I )
!         CASE ( 1 )
!            TMPARR => FNOX
!         CASE ( 2 )
!            TMPARR => OPE
!         CASE ( 3 )
!            TMPARR => MOE
!         CASE ( 4 )
!            TMPARR => DNOx
!         CASE DEFAULT
!            CALL HCO_ERROR( HcoState%Config%Err, 'I > 4', RC, THISLOC=LOC )
!            RETURN
!      END SELECT
!
!      DO I1 = 1,nT
!      DO I2 = 1,nJ
!      DO I3 = 1,nO3
!      DO I4 = 1,nSEA
!      DO I5 = 1,nSEA
!      DO I6 = 1,nJ
!      DO I7 = 1,nNOx
!         READ( fId, FMT=FMAT, IOSTAT=IOS ) TMPARR(I1,I2,I3,I4,I5,I6,I7)
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDDO
!   ENDDO

   ! Close file
   CLOSE( fId )

   ! Return w/ success
   RC = HCO_SUCCESS

 END SUBROUTINE READ_LUT_TXTFILE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_lut_txtfile
!
! !DESCRIPTION: write\_lut\_txtfile
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE WRITE_LUT_TXTFILE( HcoState, FILENAME, FNOX, DNOx, OPE, MOE, RC, &
                              T, JNO2, O3, SEA0, SEA5, JRATIO, NOX )
!
! !USES:
!
   USE inquireMod, ONLY : findFreeLUN
!
! !INPUT PARAMETERS:
!
   TYPE(HCO_State), INTENT(INOUT)  :: HcoState          ! HEMCO state obj
   CHARACTER(LEN=*),INTENT(IN)     :: FILENAME
!
! !INPUT/OUTPUT PARAMETERS:
!
   INTEGER, INTENT(INOUT)       :: RC
!
! !OUTPUT PARAMETERS:
!
   REAL*4, INTENT(IN), TARGET, DIMENSION(:,:,:,:,:,:,:) :: FNOX,OPE,MOE,DNOx
   REAL*4, INTENT(IN), OPTIONAL :: T(:),      JNO2(:), O3(:)
   REAL*4, INTENT(IN), OPTIONAL :: SEA0(:),   SEA5(:)
   REAL*4, INTENT(IN), OPTIONAL :: JRATIO(:), NOX(:)
!
! !REVISION HISTORY:
!  05 Feb 2015 - C. Keller   - Initial version modified from code provided by
!                              G.C.M. Vinken and C. Holmes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
   INTEGER             :: fId, IOS

!   INTEGER             :: I, I1, I2, I3, I4, I5, I6, I7
!   REAL*4, POINTER     :: TMPARR(:,:,:,:,:,:,:) => NULL()

   CHARACTER(LEN=255)  :: MSG
   CHARACTER(LEN=255), PARAMETER :: FMAT = "(E40.32)"
   CHARACTER(LEN=255), PARAMETER :: LOC  = &
                        'WRITE_LUT_TXTFILE (hcox_paranox_mod.F90)'

   !=================================================================
   ! WRITE_LUT_TXTFILE begins here
   !=================================================================

   ! Echo info
   IF ( Hcostate%amIRoot ) THEN
      WRITE( MSG, 100 ) TRIM( FILENAME )
      CALL HCO_MSG(HcoState%Config%Err,MSG)
100   FORMAT( 'WRITE_LUT_TXTFILE: Writing ', a )
   ENDIF

   ! Find a free file LUN
   fId = findFreeLUN()

   ! Open file for reading
   OPEN ( fID, FILE=TRIM(FILENAME), ACTION="WRITE", FORM="FORMATTED", IOSTAT=IOS )
   IF ( IOS /= 0 ) THEN
      CALL HCO_ERROR( HcoState%Config%Err, 'write_lut_txtfile:1', RC, THISLOC=LOC )
      RETURN
   ENDIF

   ! Read FNOx
   WRITE( fId, FMT=FMAT, IOSTAT=IOS ) FNOx
   IF ( IOS /= 0 ) THEN
      CALL HCO_ERROR( HcoState%Config%Err, 'write_lut_txtfile: FNOx', RC, THISLOC=LOC )
      RETURN
   ENDIF

   ! Read OPE
   WRITE( fId, FMT=FMAT, IOSTAT=IOS ) OPE
   IF ( IOS /= 0 ) THEN
      CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: OPE', RC, THISLOC=LOC )
      RETURN
   ENDIF

   ! Read MOE
   WRITE( fId, FMT=FMAT, IOSTAT=IOS ) MOE
   IF ( IOS /= 0 ) THEN
      CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: MOE', RC, THISLOC=LOC )
      RETURN
   ENDIF

   ! Read DNOx
   WRITE( fId, FMT=FMAT, IOSTAT=IOS ) DNOx
   IF ( IOS /= 0 ) THEN
      CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: DNOx', RC, THISLOC=LOC )
      RETURN
   ENDIF

   ! Read optional values
   IF ( PRESENT(T) ) THEN
      WRITE( fId, FMT=FMAT, IOSTAT=IOS ) T
      IF ( IOS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: T', RC, THISLOC=LOC )
         RETURN
      ENDIF
   ENDIF

   IF ( PRESENT(JNO2) ) THEN
      WRITE( fId, FMT=FMAT, IOSTAT=IOS ) JNO2
      IF ( IOS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: JNO2', RC, THISLOC=LOC )
         RETURN
      ENDIF
   ENDIF

   IF ( PRESENT(O3) ) THEN
      WRITE( fId, FMT=FMAT, IOSTAT=IOS ) O3
      IF ( IOS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: O3', RC, THISLOC=LOC )
         RETURN
      ENDIF
   ENDIF

   IF ( PRESENT(SEA0) ) THEN
      WRITE( fId, FMT=FMAT, IOSTAT=IOS ) SEA0
      IF ( IOS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: SEA0', RC, THISLOC=LOC )
         RETURN
      ENDIF
   ENDIF

   IF ( PRESENT(SEA5) ) THEN
      WRITE( fId, FMT=FMAT, IOSTAT=IOS ) SEA5
      IF ( IOS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: SEA5', RC, THISLOC=LOC )
         RETURN
      ENDIF
   ENDIF

   IF ( PRESENT(JRATIO) ) THEN
      WRITE( fId, FMT=FMAT, IOSTAT=IOS ) JRATIO
      IF ( IOS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: JRATIO', RC, THISLOC=LOC )
         RETURN
      ENDIF
   ENDIF

   IF ( PRESENT(NOX) ) THEN
      WRITE( fId, FMT=FMAT, IOSTAT=IOS ) NOX
      IF ( IOS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err, 'read_lut_txtfile: NOX', RC, THISLOC=LOC )
         RETURN
      ENDIF
   ENDIF

   ! Close file
   CLOSE( fId )

   ! Return w/ success
   RC = HCO_SUCCESS

 END SUBROUTINE WRITE_LUT_TXTFILE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INTERPOL_LINWEIGHTS
!
! !DESCRIPTION:  Subroutine INTERPOL\_LINWEIGHTS finds the array elements and
!      weights for piecewise 1-D linear interpolation. The input array of NODES
!      must be in monotonic ascending order. (C. Holmes 3/27/2014)
!
!      If Y is an array containing values of a function evaluated at the points
!      given in NODES, then its interpolated value at the point VALUESIN will be
!      Y(VALUEIN) =  Y(INDICES(1)) * WEIGHTS(1) +
!                    Y(INDICES(2)) * WEIGHTS(2)
!
!      This subroutine finds indices of consecutive nodes that bracket VALUEIN and
!      weights such that
!      VALUEIN = NODES(INDICES(1))   * WEIGHTS(1)     +
!                NODES(INDICES(1)+1) * (1-WEIGHTS(1))
!
!      For convenience, the returned values of INDICES and WEIGHTS are 2-element
!      arrays, where
!          INDICES(2) = INDICES(1)+1 and
!          WEIGHTS(2) = 1 - WEIGHTS(1)
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE INTERPOL_LINWEIGHTS( NODES, VALUEIN, INDICES, WEIGHTS )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
   REAL*4,INTENT(IN)   :: NODES(:), VALUEIN
!
! !OUTPUT PARAMETERS:
!
   ! These arrays are always 2 elements each, but declaring
   ! as deferred shape avoids array temporaries
   INTEGER,INTENT(OUT) :: INDICES(:)
   REAL*4, INTENT(OUT) :: WEIGHTS(:)
!
! !REVISION HISTORY:
!  03 Jun 2013 - C. Holmes      - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
   INTEGER :: I
   REAL*8  :: VALUE

   !=================================================================
   ! INTERPOL_LINWEIGHTS begins here!
   !=================================================================

   ! If larger than largest in LUT, assign largest level values
   VALUE = MIN( VALUEIN, MAXVAL( NODES ) )

   ! If smaller, assign smallest level value
   !GanLuo+VALUE = MAX( VALUE,   MINVAL( NODES ) )
   VALUE = MAX( VALUE,   MINVAL( NODES )*1.d0 )

   ! Initialize
   INDICES = (/ 1, 1 /)

   ! Loop over interpolation nodes until we find the largest node value
   ! that is less than the desired value
   DO I=1, SIZE(NODES)
      INDICES(1) = I
      IF ( VALUE <= NODES(I+1) ) EXIT
   END DO

   ! The next node
   INDICES(2) = INDICES(1) + 1

   ! Weights for the corresponding node indices
   WEIGHTS(1) = ( NODES(I+1) - VALUE ) / ( NODES(I+1) - NODES(I) )
   WEIGHTS(2) = 1.0 - WEIGHTS(1)

 END SUBROUTINE INTERPOL_LINWEIGHTS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: paranox_lut
!
! !DESCRIPTION:  Subroutine PARANOX\_LUT returns fractional remainder of
! ship NOx (FNOx), fraction of NOx that dry deposits as NOy species (DNOx),
! ozone production efficiency (OPE), and methane oxidation
! efficiency (MOE) after 5-hrs of plume aging. Values are taken taken from a
! lookup table using piecewise linear interpolation. The look-up table is derived
! from the PARANOx gaussian plume model (Vinken et al. 2011; Holmes et al. 2014)
! (G.C.M. Vinken, KNMI, June 2010; C. Holmes June 2013)
!
! The lookup table uses 8 input variables:
!     TEMP   : model temperature, K
!     JNO2   : J(NO2) value, 1/s
!     O3     : concentration O3 in ambient air, ppb
!     SEA0   : solar elevation angle at emission time 5 hours ago, degree
!     SEA5   : solar elevation angle at this time, degree
!     JRatio : ratio J(OH)/J(NO2), unitless
!     NOx    : concentration NOx in ambient air, ppt
!     WS     : wind speed, m/s
!
! In GEOS-Chem v9-01-03 through v9-02, the effects of wind speed on FNOx and OPE
! were not included (wind speed set at 6 m/s). The JRatio also used J(O1D)
! rather than J(OH); this has only a small effect on interpolated values.
! To reproduce the behavior of these earlier versions, modify code below marked
! with ******* and call READ\_PARANOX\_LUT\_v913 in emissions\_mod.F
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE PARANOX_LUT( ExtState,  HcoState, Inst, &
                         I, J, RC,  FNOX, DNOx, OPE, MOE_OUT )
!
! !USES:
!
   USE HCO_STATE_MOD,        ONLY : HCO_State
   USE HCOX_STATE_MOD,       ONLY : Ext_State
!
! !INPUT PARAMETERS:
!
   TYPE(Ext_State), POINTER    :: ExtState
   TYPE(HCO_State), POINTER    :: HcoState
   TYPE(MyInst),    POINTER    :: Inst
   INTEGER, INTENT(IN)         :: I, J      ! Grid indices
!
! !OUTPUT PARAMETERS:
!
   REAL*8, INTENT(OUT)           :: FNOX    ! fraction of NOx remaining, mol/mol
   REAL*8, INTENT(OUT)           :: DNOX    ! fraction of NOx deposited, mol/mol
   REAL*8, INTENT(OUT)           :: OPE     ! net OPE, mol(net P(O3))/mol(P(HNO3))
   REAL*8, INTENT(OUT), OPTIONAL :: MOE_OUT ! net MOE, mol(L(CH4))/mol(E(NOx))
!
! !INPUT/OUTPUT PARAMETERS:
!
   INTEGER, INTENT(INOUT)        :: RC      ! Return code
!
! !REVISION HISTORY:
!     Jun 2010 - G.C.M. Vinken - Initial version
!  03 Jun 2013 - C. Holmes     - Heavily modified and simplified from previous
!                                LUT interpolation code by G.C.M. Vinken and
!                                M. Payer. LUT now includes wind speed.
!  04 Feb 2015 - C. Keller     - Updated for use in HEMCO.
!  24 Sep 2015 - E. Lundgren   - ExtState vars O3, NO2, and NO now in
!                                kg/kg dry air (previously kg)
!  07 Jan 2016 - E. Lundgren   - Update H2O molec wt to match GC value
!  20 Apr 2016 - M. Sulprizio  - Remove calculation of J(OH). We now get J(OH),
!                                the effective rate for O3+hv(+H2O)->OH+OH,
!                                directly from FAST-JX. In FlexChem, adjustment
!                                of the photolysis rates are done in routine
!                                PHOTRATE_ADJ (found in GeosCore/fast_jx_mod.F).
!  20 Sep 2016 - R. Yantosca   - Replace non-standard ASIND function with ASIN,
!                                and convert to degrees (divide by PI/180)
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER                    :: I1,I2,I3,I4,I5,I6,I7,I8
   REAL(sp)                   :: FNOX_TMP, DNOX_TMP, OPE_TMP, MOE_TMP
   REAL(sp)                   :: WEIGHT
   REAL(sp)                   :: JNO2, JOH, TAIR
   REAL(sp)                   :: AIR
   REAL*8                     :: MOE

   ! Interpolation variables, indices, and weights
   REAL(sp), DIMENSION(8)     :: VARS
   INTEGER,  DIMENSION(8,2)   :: INDX
   REAL(sp), DIMENSION(8,2)   :: WTS

   REAL(sp), POINTER          :: FRACNOX_LUT(:,:,:,:,:,:,:)
   REAL(sp), POINTER          :: DNOX_LUT   (:,:,:,:,:,:,:)
   REAL(sp), POINTER          :: OPE_LUT    (:,:,:,:,:,:,:)
   REAL(sp), POINTER          :: MOE_LUT    (:,:,:,:,:,:,:)

   CHARACTER(LEN=255)         :: MSG
   CHARACTER(LEN=255)         :: LOC = 'PARANOX_LUT'

   !=================================================================
   ! PARANOX_LUT begins here!
   !=================================================================

   ! Initialize for safety's sake
   FRACNOX_LUT => NULL()
   DNOX_LUT    => NULL()
   OPE_LUT     => NULL()
   MOE_LUT     => NULL()
   FNOX        = 0.0_sp
   DNOX        = 0.0_sp
   OPE         = 0.0_sp
   IF ( PRESENT( MOE_OUT ) ) THEN
      MOE_OUT  = 0.0_sp
   ENDIF

   ! Air mass [kg]
   AIR = ExtState%AIR%Arr%Val(I,J,1)

   ! Air temperature, K
   Tair = ExtState%T2M%Arr%Val(I,J)

!   ! for debugging only
!   if(I==3.and.J==35)then
!      write(*,*) 'Call PARANOX_LUT @ ',I,J
!      write(*,*) 'Tair: ', Tair
!      write(*,*) 'SUNCOSmid: ', SC5(I,J)
!   endif

   ! Check if sun is up
   IF ( ExtState%SUNCOS%Arr%Val(I,J) > 0.0_hp ) THEN

      ! J(NO2), 1/s
      JNO2 = ExtState%JNO2%Arr%Val(I,J)

!      ! J(O1D), 1/s
!      JO1D = ExtState%JO1D%Arr%Val(I,J)
!
!      ! H2O, molec/cm3. Get from specific humidity, which is in kg/kg.
!      ! NOTE: SPHU is mass H2O / mass total air so use of dry air molecular
!      ! weight is slightly inaccurate. C (ewl, 9/11/15)
!      H2O = ExtState%SPHU%Arr%Val(I,J,1) * DENS &
!          * HcoState%Phys%AIRMW / MWH2O
!
!      ! Calculate J(OH), the effective rate for O3+hv -> OH+OH,
!      ! assuming steady state for O(1D).
!      ! Rate coefficients are cm3/molec/s; concentrations are molec/cm3
!      ! This should match the O3+hv (+H2O) -> OH+OH kinetics in calcrate.F
!      JOH = JO1D *                                            &
!            1.63e-10 * EXP( 60.e0/Tair) * H2O /               &
!          ( 1.63e-10 * EXP( 60.e0/Tair) * H2O +             &
!            1.20e-10                    * DENS * 0.5000e-6  + &
!            2.15e-11 * EXP(110.e0/Tair) * DENS * 0.7808e0   + &
!            3.30e-11 * EXP( 55.e0/Tair) * DENS * 0.2095e0   )

      ! J(OH) - effective rate for O3+hv(+H2O)-> OH+OH, 1/s
      JOH = ExtState%JOH%Arr%Val(I,J)

   ELSE

      ! J-values are zero when sun is down
      JNO2 = 0e0_sp
      JOH  = 0e0_sp

   ENDIF

!   ! for debugging only
!   if(I==3.and.J==35)then
!      write(*,*) 'JNO2: ', JNO2
!      write(*,*) 'JOH : ', JOH
!   endif

   !========================================================================
   ! Load all variables into a single array
   !========================================================================

   ! Temperature, K
   VARS(1) = Tair

   ! J(NO2), 1/s
   VARS(2) = JNO2

! old
!   ! O3 concentration in ambient air, ppb
!   VARS(3) = ExtState%O3%Arr%Val(I,J,1) / AIR   &
!           * HcoState%Phys%AIRMW        / MW_O3 &
!           * 1.e9_sp
! new
   ! O3 concentration in ambient air, ppb
   ! NOTE: ExtState%O3 units are now kg/kg dry air (ewl, 9/11/15)
   VARS(3) = ExtState%O3%Arr%Val(I,J,1)         &
           * HcoState%Phys%AIRMW        / Inst%MW_O3 &
           * 1.e9_sp
! end new (ewl)

   ! Solar elevation angle, degree
   ! SEA0 = SEA when emitted from ship, 5-h before current model time
   ! SEA5 = SEA at current model time, 5-h after emission from ship
   ! Note: Since SEA = 90 - SZA, then cos(SZA) = sin(SEA) and
   ! thus SEA = arcsin( cos( SZA ) )
   !VARS(4) = ASIND( SC5(I,J) )
   !VARS(5) = ASIND( ExtState%SUNCOS%Arr%Val(I,J) )
   VARS(4) = ASIN( Inst%SC5(I,J)                ) / HcoState%Phys%PI_180
   VARS(5) = ASIN( ExtState%SUNCOS%Arr%Val(I,J) ) / HcoState%Phys%PI_180

   ! J(OH)/J(NO2), unitless
   ! Note J(OH) is the loss rate (1/s) of O3 to OH, which accounts for
   ! the temperature, pressure and water vapor dependence of these reactions
   VARS(6) = 0.0_sp
   IF ( JNO2 /= 0.0_sp ) THEN
      IF ( (EXPONENT(JOH)-EXPONENT(JNO2)) < MAXEXPONENT(JOH) ) THEN
         VARS(6) = JOH / JNO2
      ENDIF
   ENDIF

! old
!   ! NOx concetration in ambient air, ppt
!   VARS(7) = ( ( ExtState%NO%Arr%Val(I,J,1)  / AIR        &
!           *     HcoState%Phys%AIRMW         / MW_NO  )   &
!           +   ( ExtState%NO2%Arr%Val(I,J,1) / AIR        &
!           *     HcoState%Phys%AIRMW         / MW_NO2 ) ) &
!           * 1.e12_sp
! new
   ! NOx concetration in ambient air, ppt
   ! NOTE: ExtState vars NO and NO2 units are now kg/kg dry air (ewl, 9/11/15)
   VARS(7) = ( ( ExtState%NO%Arr%Val(I,J,1)               &
           *     HcoState%Phys%AIRMW         / Inst%MW_NO  )   &
           +   ( ExtState%NO2%Arr%Val(I,J,1)              &
           *     HcoState%Phys%AIRMW         / Inst%MW_NO2 ) ) &
           * 1.e12_sp
! end new (ewl)

      ! Wind speed, m/s
   VARS(8) = SQRT( ExtState%U10M%Arr%Val(I,J)**2 &
           +       ExtState%V10M%Arr%Val(I,J)**2 )

!   ! for debugging only
!   if(I==1.and.J==35)then
!      write(*,*) 'VARS(1)      : ', VARS(1)
!      write(*,*) 'VARS(2)      : ', VARS(2)
!      write(*,*) 'VARS(3)      : ', VARS(3)
!      write(*,*) 'VARS(4)      : ', VARS(4)
!      write(*,*) 'VARS(5)      : ', VARS(5)
!      write(*,*) 'VARS(6)      : ', VARS(6)
!      write(*,*) 'VARS(7)      : ', VARS(7)
!      write(*,*) 'VARS(8)      : ', VARS(8)
!      write(*,*) 'AIR          : ', ExtState%AIR%Arr%Val(I,J,1)
!      write(*,*) 'AIRMW        : ', HcoState%Phys%AIRMW
!      write(*,*) 'MWNO, NO2, O3: ', Inst%MW_NO, Inst%MW_NO2, Inst%MW_O3
!      write(*,*) 'O3conc       : ', ExtState%O3%Arr%Val(I,J,1)
!      write(*,*) 'NO,NO2 conc  : ', ExtState%NO%Arr%Val(I,J,1), &
!                                    ExtState%NO2%Arr%Val(I,J,1)
!      write(*,*) 'U, V         : ', ExtState%U10M%Arr%Val(I,J), &
!                                    ExtState%V10M%Arr%Val(I,J)
!   endif

      !*****************************************************
      ! Restoring the following lines reproduces the behavior of
      ! GEOS-Chem v9-01-03 through v9-02
      ! the LUT was indexed with the ratio J(O1D)/J(NO2) (cdh, 3/27/2014)
!      VARS(6)   = SAFE_DIV( JO1D, JNO2, 0D0 )
!      JRATIOlev = (/ 5.e-4, 0.0015, 0.0025, 0.0055 /)
      !*****************************************************

   !========================================================================
   ! Find the indices of nodes and their corresponding weights for the
   ! interpolation
   !========================================================================

   ! Temperature:
   CALL INTERPOL_LINWEIGHTS( Inst%Tlev, VARS(1), INDX(1,:), WTS(1,:) )

   ! J(NO2):
   CALL INTERPOL_LINWEIGHTS( Inst%JNO2lev, VARS(2), INDX(2,:), WTS(2,:) )

   ! [O3]:
   CALL INTERPOL_LINWEIGHTS( Inst%O3lev, VARS(3), INDX(3,:), WTS(3,:) )

   ! SEA0:
   CALL INTERPOL_LINWEIGHTS( Inst%SEA0lev, VARS(4), INDX(4,:), WTS(4,:) )

   ! SEA5:
   CALL INTERPOL_LINWEIGHTS( Inst%SEA5lev, VARS(5), INDX(5,:), WTS(5,:) )

   ! JRATIO:
   CALL INTERPOL_LINWEIGHTS( Inst%JRATIOlev, VARS(6), INDX(6,:), WTS(6,:))

   ! [NOx]:
   CALL INTERPOL_LINWEIGHTS( Inst%NOXlev, VARS(7), INDX(7,:), WTS(7,:) )

   ! Wind speed:
   CALL INTERPOL_LINWEIGHTS( Inst%WSlev, VARS(8), INDX(8,:), WTS(8,:) )

   !========================================================================
   ! Piecewise linear interpolation
   !========================================================================

   ! Initialize
   FNOX = 0.0d0
   DNOx = 0.0d0
   OPE  = 0.0d0
   MOE  = 0.0d0

   ! Loop over wind speed
   DO I8=1,2

      ! Point at the LUT for this wind speed
      ! Last two digits in fortran variable names indicate wind speed in m/s
      SELECT CASE ( NINT( Inst%WSlev(INDX(8,I8)) ) )
         CASE (  2 )
            FRACNOX_LUT => Inst%FRACNOX_LUT02
            DNOx_LUT    => Inst%DNOx_LUT02
            OPE_LUT     => Inst%OPE_LUT02
            MOE_LUT     => Inst%MOE_LUT02
         CASE (  6 )
            FRACNOX_LUT => Inst%FRACNOX_LUT06
            DNOx_LUT    => Inst%DNOx_LUT06
            OPE_LUT     => Inst%OPE_LUT06
            MOE_LUT     => Inst%MOE_LUT06
         CASE ( 10 )
            FRACNOX_LUT => Inst%FRACNOX_LUT10
            DNOx_LUT    => Inst%DNOx_LUT10
            OPE_LUT     => Inst%OPE_LUT10
            MOE_LUT     => Inst%MOE_LUT10
         CASE ( 14 )
            FRACNOX_LUT => Inst%FRACNOX_LUT14
            DNOx_LUT    => Inst%DNOx_LUT14
            OPE_LUT     => Inst%OPE_LUT14
            MOE_LUT     => Inst%MOE_LUT14
         CASE ( 18 )
            FRACNOX_LUT => Inst%FRACNOX_LUT18
            DNOx_LUT    => Inst%DNOx_LUT18
            OPE_LUT     => Inst%OPE_LUT18
            MOE_LUT     => Inst%MOE_LUT18
         CASE DEFAULT
             MSG = 'LUT error: Wind speed interpolation error!'
             CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
             RETURN
      END SELECT

      !*****************************************************
      ! Restoring the following lines reproduces the behavior of
      ! GEOS-Chem v9-01-03 through v9-02 in which wind speed
      ! effects on FNOx and OPE are neglected (cdh, 3/25/2014)
!      FRACNOX_LUT => FRACNOX_LUT_v913
!      OPE_LUT     => OPE_LUT_v913
      !*****************************************************

      ! loop over all other variables
      DO I7=1,2
      DO I6=1,2
      DO I5=1,2
      DO I4=1,2
      DO I3=1,2
      DO I2=1,2
      DO I1=1,2

         !------------------------------------------
         ! Nodes and weights used in the interpolation
         !------------------------------------------

         ! Fraction NOx from the LUT
         FNOX_TMP = FRACNOX_LUT( INDX(1,I1), INDX(2,I2), INDX(3,I3), &
                                 INDX(4,I4), INDX(5,I5),             &
                                 INDX(6,I6), INDX(7,I7) )

         DNOX_TMP = DNOx_LUT(    INDX(1,I1), INDX(2,I2), INDX(3,I3), &
                                 INDX(4,I4), INDX(5,I5),             &
                                 INDX(6,I6), INDX(7,I7) )

         ! OPE from the LUT
         OPE_TMP  = OPE_LUT(     INDX(1,I1), INDX(2,I2), INDX(3,I3), &
                                 INDX(4,I4), INDX(5,I5),             &
                                 INDX(6,I6), INDX(7,I7) )

         ! MOE from the LUT
         MOE_TMP  = MOE_LUT(     INDX(1,I1), INDX(2,I2), INDX(3,I3), &
                                 INDX(4,I4), INDX(5,I5),             &
                                 INDX(6,I6), INDX(7,I7) )

         ! Interpolation weight for this element
         WEIGHT = WTS(1,I1) * WTS(2,I2) * WTS(3,I3) * WTS(4,I4) * &
                  WTS(5,I5) * WTS(6,I6) * WTS(7,I7) * WTS(8,I8)

         !-----------------------------------
         ! Error Check
         !-----------------------------------

         !IF ENCOUNTER -999 IN THE LUT PRINT ERROR!!
         IF ( ( FNOX_TMP < 0. ) .or. ( FNOX_TMP > 1. ) ) THEN

            PRINT*, 'PARANOX_LUT: fracnox = ', FNOX_TMP
            PRINT*, 'This occured at grid box ', I, J
            PRINT*, 'Lon/Lat: ', HcoState%Grid%XMID%Val(I,J), HcoState%Grid%YMID%Val(I,J)
            PRINT*, 'SZA 5 hours ago : ', VARS(4)
            PRINT*, 'SZA at this time: ', VARS(5)
            PRINT*, 'The two SZAs should not be more than 75 deg apart!'
            PRINT*, 'If they are, your restart SZA might be wrong.'
            PRINT*, 'You can try to coldstart PARANOx by commenting'
            PRINT*, 'all PARANOX_SUNCOS entries in your HEMCO'
            PRINT*, 'configuration file.'

            !print*, I1, I2, I3, I4, I5, I6, I7, I8
            !print*, INDX(1,I1), INDX(2,I2), INDX(3,I3),  INDX(4,I4), &
            !        INDX(5,I5), INDX(6,I6), INDX(7,I7), INDX(8,I8)
            !print*, VARS

            MSG = 'LUT error: Fracnox should be between 0 and 1!'
            CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
            RETURN
         ENDIF

         !-----------------------------------
         ! Final interpolated values
         !-----------------------------------

         ! Weighted sum of FNOx from the LUT
         FNOx = FNOx + FNOX_TMP * WEIGHT

         ! Weighted sum of DNOx from the LUT
         DNOx = DNOx + DNOX_TMP * WEIGHT

         ! Weighted sum of OPE from the LUT
         OPE  = OPE + OPE_TMP * WEIGHT

         ! Weighted sum of MOE from the LUT
         MOE  = MOE + MOE_TMP * WEIGHT

      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO

      ! Free pointers
      FRACNOX_LUT => NULL()
      DNOx_LUT    => NULL()
      OPE_LUT     => NULL()
      MOE_LUT     => NULL()
   END DO

   ! Transfer MOE if optional output parameter is present
   IF ( PRESENT( MOE_OUT ) ) MOE_OUT = MOE

   ! Nullify pointers
   NULLIFY( FRACNOX_LUT )
   NULLIFY( DNOx_LUT )
   NULLIFY( OPE_LUT  )
   NULLIFY( MOE_LUT  )

   ! Return w/ success
   RC = HCO_SUCCESS

 END SUBROUTINE PARANOX_LUT
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstGet
!
! !DESCRIPTION: Subroutine InstGet returns a poiner to the desired instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstGet ( Instance, Inst, RC, PrevInst )
!
! !INPUT PARAMETERS:
!
    INTEGER                             :: Instance
    TYPE(MyInst),     POINTER           :: Inst
    INTEGER                             :: RC
    TYPE(MyInst),     POINTER, OPTIONAL :: PrevInst
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst),     POINTER    :: PrvInst

    !=================================================================
    ! InstGet begins here!
    !=================================================================

    ! Get instance. Also archive previous instance.
    PrvInst => NULL()
    Inst    => AllInst
    DO WHILE ( ASSOCIATED(Inst) )
       IF ( Inst%Instance == Instance ) EXIT
       PrvInst => Inst
       Inst    => Inst%NextInst
    END DO
    IF ( .NOT. ASSOCIATED( Inst ) ) THEN
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Pass output arguments
    IF ( PRESENT(PrevInst) ) PrevInst => PrvInst

    ! Cleanup & Return
    PrvInst => NULL()
    RC = HCO_SUCCESS

  END SUBROUTINE InstGet
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstCreate
!
! !DESCRIPTION: Subroutine InstCreate creates a new instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstCreate ( ExtNr, Instance, Inst, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,       INTENT(IN)       :: ExtNr
!
! !OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(  OUT)    :: Instance
    TYPE(MyInst),  POINTER          :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(INOUT)    :: RC
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst), POINTER          :: TmpInst
    INTEGER                        :: nnInst

    !=================================================================
    ! InstCreate begins here!
    !=================================================================

    ! ----------------------------------------------------------------
    ! Generic instance initialization
    ! ----------------------------------------------------------------

    ! Initialize
    Inst => NULL()

    ! Get number of already existing instances
    TmpInst => AllInst
    nnInst = 0
    DO WHILE ( ASSOCIATED(TmpInst) )
       nnInst  =  nnInst + 1
       TmpInst => TmpInst%NextInst
    END DO

    ! Create new instance
    ALLOCATE(Inst)
    Inst%Instance = nnInst + 1
    Inst%ExtNr    = ExtNr

    ! Attach to instance list
    Inst%NextInst => AllInst
    AllInst       => Inst

    ! Update output instance
    Instance = Inst%Instance

    ! ----------------------------------------------------------------
    ! Type specific initialization statements follow below
    ! ----------------------------------------------------------------

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE InstCreate
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstRemove
!
! !DESCRIPTION: Subroutine InstRemove creates a new instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstRemove ( Instance )
!
! !INPUT PARAMETERS:
!
    INTEGER                         :: Instance
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER                     :: RC
    TYPE(MyInst), POINTER       :: PrevInst
    TYPE(MyInst), POINTER       :: Inst

    !=================================================================
    ! InstRemove begins here!
    !=================================================================

    ! Init
    PrevInst => NULL()
    Inst     => NULL()

    ! Get instance. Also archive previous instance.
    CALL InstGet ( Instance, Inst, RC, PrevInst=PrevInst )

    ! Instance-specific deallocation
    IF ( ASSOCIATED(Inst) ) THEN

       ! Pop off instance from list
       IF ( ASSOCIATED(PrevInst) ) THEN

          IF ( ASSOCIATED( Inst%ShipNO) ) DEALLOCATE ( Inst%ShipNO )
          IF ( ASSOCIATED( Inst%SC5   ) ) DEALLOCATE ( Inst%SC5    )

          IF ( ASSOCIATED( Inst%FRACNOX_LUT02 ) ) DEALLOCATE( Inst%FRACNOX_LUT02 )
          IF ( ASSOCIATED( Inst%FRACNOX_LUT06 ) ) DEALLOCATE( Inst%FRACNOX_LUT06 )
          IF ( ASSOCIATED( Inst%FRACNOX_LUT10 ) ) DEALLOCATE( Inst%FRACNOX_LUT10 )
          IF ( ASSOCIATED( Inst%FRACNOX_LUT14 ) ) DEALLOCATE( Inst%FRACNOX_LUT14 )
          IF ( ASSOCIATED( Inst%FRACNOX_LUT18 ) ) DEALLOCATE( Inst%FRACNOX_LUT18 )

          IF ( ASSOCIATED( Inst%OPE_LUT02     ) ) DEALLOCATE( Inst%OPE_LUT02     )
          IF ( ASSOCIATED( Inst%OPE_LUT06     ) ) DEALLOCATE( Inst%OPE_LUT06     )
          IF ( ASSOCIATED( Inst%OPE_LUT10     ) ) DEALLOCATE( Inst%OPE_LUT10     )
          IF ( ASSOCIATED( Inst%OPE_LUT14     ) ) DEALLOCATE( Inst%OPE_LUT14     )
          IF ( ASSOCIATED( Inst%OPE_LUT18     ) ) DEALLOCATE( Inst%OPE_LUT18     )

          IF ( ASSOCIATED( Inst%MOE_LUT02     ) ) DEALLOCATE( Inst%MOE_LUT02     )
          IF ( ASSOCIATED( Inst%MOE_LUT06     ) ) DEALLOCATE( Inst%MOE_LUT06     )
          IF ( ASSOCIATED( Inst%MOE_LUT10     ) ) DEALLOCATE( Inst%MOE_LUT10     )
          IF ( ASSOCIATED( Inst%MOE_LUT14     ) ) DEALLOCATE( Inst%MOE_LUT14     )
          IF ( ASSOCIATED( Inst%MOE_LUT18     ) ) DEALLOCATE( Inst%MOE_LUT18     )

          IF ( ASSOCIATED( Inst%DNOx_LUT02    ) ) DEALLOCATE( Inst%DNOx_LUT02    )
          IF ( ASSOCIATED( Inst%DNOx_LUT06    ) ) DEALLOCATE( Inst%DNOx_LUT06    )
          IF ( ASSOCIATED( Inst%DNOx_LUT10    ) ) DEALLOCATE( Inst%DNOx_LUT10    )
          IF ( ASSOCIATED( Inst%DNOx_LUT14    ) ) DEALLOCATE( Inst%DNOx_LUT14    )
          IF ( ASSOCIATED( Inst%DNOx_LUT18    ) ) DEALLOCATE( Inst%DNOx_LUT18    )

          IF ( ASSOCIATED( Inst%DEPO3         ) ) DEALLOCATE( Inst%DEPO3         )
          IF ( ASSOCIATED( Inst%DEPHNO3       ) ) DEALLOCATE( Inst%DEPHNO3       )

          PrevInst%NextInst => Inst%NextInst
       ELSE
          AllInst => Inst%NextInst
       ENDIF
       DEALLOCATE(Inst)
       Inst => NULL()
    ENDIF

   END SUBROUTINE InstRemove
!EOC
END MODULE HCOX_ParaNOx_mod
