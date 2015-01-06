!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_paranox_mod.F90
!
! !DESCRIPTION: Module HCOX\_PARANOX\_MOD contains routines to 
! compute ship emissions and associated concentrations of NO, HNO3 and
! O3 from NO ship emission data.  This follows the implementation of
! the PARANOX ship plume model in GEOS-Chem.
!\\
!\\
! This is a HEMCO extension module that uses many of the HEMCO core
! utilities.
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
!  15 Oct 2013 - C. Keller   - Now a HEMCO extension
!  06 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  06 Jun 2014 - R. Yantosca - Now indended with F90 free-format
!  25 Jun 2014 - R. Yantosca - Now pass the look-up-table filenames
!  22 Jul 2014 - R. Yantosca - Added shadow copy of FAST-JX function FJXFUNC
!  28 Jul 2014 - C. Keller   - Now pass J-Values through ExtState. This makes
!                              the FJXFUNC shadow copy obsolete
!  13 Aug 2014 - C. Keller   - Added manual diagnostics
!  16 Oct 2014 - C. Keller   - Now store SUNCOSmid values internally over the
!                              past 5 hours and use these values for SUNCOSmid5.
!                              This is required for standalone mode.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:

  ! Scalars
  INTEGER                       :: ExtNr
  INTEGER                       :: IDTNO 
  INTEGER                       :: IDTHNO3
  INTEGER                       :: IDTO3
  CHARACTER(LEN=255)            :: FracNox_FILE
  CHARACTER(LEN=255)            :: IntOPE_FILE
  REAL*8                        :: MW_O3
  REAL*8                        :: MW_NO
  REAL*8                        :: MW_NO2
  REAL*8                        :: MW_HNO3
  REAL*8                        :: MW_AIR

  ! Arrays
  REAL(hp), ALLOCATABLE, TARGET :: ShipNO(:,:,:)

  ! For SunCosMid 5hrs ago
  REAL(hp), ALLOCATABLE         :: SC5(:,:,:)
  INTEGER                       :: SC5ID

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
  SUBROUTINE HCOX_ParaNOx_Run( am_I_Root, ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_Calc_Mod, ONLY : HCO_CalcEmis
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   ) :: am_I_Root   ! Are we on the root CPU?
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

    !=================================================================
    ! HCOX_PARANOX_RUN begins here!
    !=================================================================

    ! Return if extension disabled 
    IF ( .NOT. ExtState%ParaNOx ) RETURN

    ! Enter
    CALL HCO_Enter( 'HCOX_ParaNOx_Run (hcox_paranox_mod.F90)', RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

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
    HcoState%Options%ExtNr  = ExtNr 

    ! --> Define array to write emissions into. ShipNO is reset to
    ! zero within subroutine EVOLVE_PLUME, so no need to do this
    ! here.
!    ShipNO                      = 0.0d0
    HcoState%Options%AutoFillDiagn = .FALSE.
    HcoState%Options%FillBuffer    =  .TRUE.
    HcoState%Buffer3D%Val          => ShipNO 
      
    ! Calculate ship NO emissions and write them into the ShipNO
    ! array [kg/m2/s]. 
    CALL HCO_CalcEmis( am_I_Root, HcoState, .FALSE., RC )
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
    CALL Evolve_Plume( am_I_Root, ExtState, ShipNO, HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Leave w/ success
    CALL HCO_Leave( RC ) 

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
  SUBROUTINE Evolve_Plume( am_I_Root, ExtState, ShipNoEmis, HcoState, RC )  
!
! !USES:
!
    USE HCO_FluxArr_mod,  ONLY : HCO_EmisAdd
    USE HCO_FluxArr_mod,  ONLY : HCO_DepvAdd
    USE HCO_Clock_Mod,    ONLY : HcoClock_Get
    USE ParaNOx_Util_Mod, ONLY : Interpolate_LUT2
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root          ! Root CPU?
    TYPE(Ext_State), POINTER        :: ExtState           ! External data
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(hp),        INTENT(INOUT)  :: ShipNoEmis(:,:,:)  ! Emissions
    TYPE(HCO_State), POINTER        :: HcoState           ! HEMCO State obj
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, MM
    LOGICAL            :: ERR
    LOGICAL, SAVE      :: FIRST = .TRUE.
    REAL*8             :: JNO2, JO1D, TS, SUNCOSmid5, SUNCOSmid
    REAL*8             :: O3molec, NOmolec, NO2molec, AIRmolec
    REAL*4             :: FRACTION_NOx, INT_OPE
    REAL(hp)           :: iFlx, TMP
    CHARACTER(LEN=255) :: MSG

    ! Arrays
    REAL(hp), TARGET   :: FLUXNO  (HcoState%NX,HcoState%NY)
    REAL(hp), TARGET   :: FLUXHNO3(HcoState%NX,HcoState%NY)
    REAL(hp), TARGET   :: FLUXO3  (HcoState%NX,HcoState%NY)
    REAL(hp), TARGET   :: DEPO3   (HcoState%NX,HcoState%NY)

    ! Pointers
    REAL(hp), POINTER  :: Arr2D(:,:) => NULL()

    ! For diagnostics
    REAL(hp), TARGET   :: DIAGN   (HcoState%NX,HcoState%NY,4)
    LOGICAL, SAVE      :: DODIAGN = .FALSE.
    CHARACTER(LEN=31)  :: DiagnName
    TYPE(DiagnCont), POINTER :: TmpCnt => NULL()

    ! For internal SC5 array
    INTEGER            :: HH
    INTEGER, SAVE      :: lastHH = -1

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
    CALL HCO_Enter( 'Evolve_Plume (hcox_paranox_mod.F90)', RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Leave here if none of the tracers defined
    IF ( IDTNO <= 0 .AND. IDTO3 <= 0 .AND. IDTHNO3 <= 0 ) THEN
       RC = HCO_SUCCESS 
       RETURN
    ENDIF

    ! Get simulation month
    CALL HcoClock_Get( cMM=MM, cH=HH, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! On first call, see if we need to write internal diagnostics
    IF ( FIRST ) THEN
       ! See if we have to write out manual diagnostics
       IF ( .NOT. DoDiagn ) THEN
          DiagnName = 'PARANOX_NOXFRAC_REMAINING'
          CALL DiagnCont_Find ( -1, -1, -1, -1, -1, DiagnName, 0, DoDiagn, TmpCnt )
          TmpCnt => NULL()
       ENDIF
       IF ( .NOT. DoDiagn ) THEN
          DiagnName = 'PARANOX_O3_PRODUCTION'
          CALL DiagnCont_Find ( -1, -1, -1, -1, -1, DiagnName, 0, DoDiagn, TmpCnt )
          TmpCnt => NULL()
       ENDIF
       IF ( .NOT. DoDiagn ) THEN
          DiagnName = 'PARANOX_TOTAL_SHIPNOX'
          CALL DiagnCont_Find ( -1, -1, -1, -1, -1, DiagnName, 0, DoDiagn, TmpCnt )
          TmpCnt => NULL()
       ENDIF    
       IF ( .NOT. DoDiagn ) THEN
          DiagnName = 'PARANOX_OPE'
          CALL DiagnCont_Find ( -1, -1, -1, -1, -1, DiagnName, 0, DoDiagn, TmpCnt )
          TmpCnt => NULL()
       ENDIF  

       ! Also make sure that the SC5 array holds values. Initialize them to
       ! current one until we have gone through an entire 5-hour simulation cycle.
       DO I = 1,6
          SC5(:,:,I) = ExtState%SUNCOSmid%Arr%Val(:,:)
       ENDDO

       ! Not first call any more...
       FIRST = .FALSE.  
    ENDIF
    IF ( DoDiagn ) DIAGN(:,:,:) = 0.0_hp

    ! Update current active index in array SC5. This array holds the
    ! SUNCOSmid values of the past 5 runs. They become stored chronologically, 
    ! i.e. on the first time step, we archive the current SUNCOS in slice 1, on
    ! the second time step in slice 2, etc. Index SC5ID refers to the slice 
    ! that holds the SUNCOS value from 5 hours ago. It becomes moved forward
    ! every time the simulation hour changes, and the SUNCOS values from the
    ! previous time step become written into the formerly active index. This 
    ! index will again become SC5ID in 5 hours from now!
    IF ( HH /= lastHH ) THEN
   
       ! Copy SUNCOSmid from last time step from buffer into current slot
       SC5(:,:,SC5ID) = SC5(:,:,6)

       ! Archive current SUNCOSmid for future.
       SC5(:,:,6) = ExtState%SUNCOSmid%Arr%Val(:,:)

       ! Increase index by 1. Cycle back to one if we hit end of array
       SC5ID = SC5ID + 1
       IF ( SC5ID > 5 ) SC5ID = 1

    ENDIF       

    ! Error check
    ERR = .FALSE.

    ! Init
    FLUXNO   = 0.0_hp
    FLUXHNO3 = 0.0_hp
    FLUXO3   = 0.0_hp
    DEPO3    = 0.0_hp

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
#if defined( NULL )
!$OMP PARALLEL DO                                                   &
!$OMP DEFAULT( SHARED )                                             &
!$OMP PRIVATE( I, J, L, RC, iFlx, NK, SPECNAME, JNO2, JO1D, TMP   ) &
!$OMP PRIVATE( O3molec, NOmolec, NO2molec, AIRmolec, INT_OPE      ) &
!$OMP PRIVATE( FRACTION_NOx, LMAX, TS, SUNCOSmid5, SUNCOSmid, MSG ) &     
!$OMP SCHEDULE( DYNAMIC )
#endif
    DO J = 1, HcoState%NY
    DO I = 1, HcoState%NX

       ! Skip if no ship emissions in this grid box
       IF ( ShipNoEmis(I,J,1) == 0d0 ) CYCLE

       !--------------------------------------------------------------------
       ! Get J-Values for J(NO2) and J(O3)
       ! These values are provided through the ExtState object
       !--------------------------------------------------------------------
       JNO2 = ExtState%JNO2%Arr%Val(I,J)
       JO1D = ExtState%JO1D%Arr%Val(I,J)

       !--------------------------------------------------------------------
       ! Determine fraction of NOx remaining and integrated
       ! Ozone
       ! Production Efficiency for ship emiss
       ! (gvinken,mpayer,2/7/12)
       ! Uses surface-layer concentrations of O3, NO, NO2 [molec] and 
       ! air mass (molec).
       ! Note: the ExtState concentrations are in kg, not molecules. In 
       ! INTERPOLATE_LUT2, the concentrations are used to determine nox
       ! and o3 in pptv and ppbv, respectively. So we can just normalize
       ! by the respective mol. weights and will get the same results as
       ! when using [molec]. (ckeller, 08/06/13) 
       !--------------------------------------------------------------------
       O3molec    = ExtState%O3%Arr%Val(I,J,1)  / MW_O3
       NOmolec    = ExtState%NO%Arr%Val(I,J,1)  / MW_NO
       NO2molec   = ExtState%NO2%Arr%Val(I,J,1) / MW_NO2
       AIRmolec   = ExtState%AIR%Arr%Val(I,J,1) / MW_AIR
       TS         = ExtState%T2M%Arr%Val(I,J)
       SUNCOSmid5 = SC5(I,J,SC5ID)
       SUNCOSmid  = ExtState%SUNCOSmid%Arr%Val(I,J)

       CALL Interpolate_Lut2( HcoState,     I,    J,      &
                              O3molec,      NOmolec,      &
                              NO2molec,     AIRmolec,     &
                              JO1D,         JNO2,         &
                              TS,           SUNCOSmid5,   &
                              SUNCOSmid,                  &
                              FRACTION_NOx, INT_OPE, RC    )
       IF ( RC /= HCO_SUCCESS ) THEN
          ERR = .TRUE.
          EXIT
       ENDIF

       !---------------------------
       ! Calculate NO emissions
       !---------------------------
       IF ( IDTNO > 0 ) THEN
          
          ! Of the total ship NOx, the fraction FRACTION_NOX
          ! survives after plume dilution and chemistry.
          ! Unit: kg/m2/s 
          FLUXNO(I,J) = ShipNoEmis(I,J,1) * FRACTION_NOx

       ENDIF

       !---------------------------
       ! Calculate HNO3 emissions
       !---------------------------
       IF ( IDTHNO3 > 0 ) THEN

          ! Of the total ship NOx, the fraction 1-FRACTION_NOX
          ! is converted to HNO3 during plume dilution and chemistry. 
          ! Unit: kg/m2/s 
          FLUXHNO3(I,J) = ShipNoEmis(I,J,1) * ( 1d0 - FRACTION_NOx ) &
                        * ( MW_HNO3 / MW_NO )
       ENDIF

       !---------------------------
       ! Calculate O3 emissions
       !---------------------------
       IF ( IDTO3 > 0 ) THEN
          
          ! Of the total ship NOx, the fraction
          ! (1-FRACTION_NOX)*INT_OPE is converted to O3 during 
          ! plume dilution and chemistry. 
          ! Unit: kg/m2/s 
          iFlx = ShipNoEmis(I,J,1) * (1d0-FRACTION_NOx) * INT_OPE &
               * ( MW_O3 / MW_NO )

          ! For positive fluxes, add to emission flux array 
          IF ( iFlx >= 0.0_hp ) THEN
             FLUXO3(I,J) = iFlx

          ! For negative fluxes, calculate deposition velocity based
          ! on current surface O3 concentration and pass to deposition
          ! array
          ELSE

             ! Calculate deposition velocity (1/s) from flux
             ! NOTE: the calculated deposition flux is in kg/m2/s,
             ! which has to be converted to 1/s. Use here the O3 conc.
             ! [kg] of the lowest model box.
             ! Now avoid div-zero error (ckeller, 11/10/2014).
             IF ( ExtState%O3%Arr%Val(I,J,1) > 0.0_hp ) THEN
                TMP = ABS(iFlx) * HcoState%Grid%AREA_M2%Val(I,J)

                ! Check if it's safe to do division
                IF ( (EXPONENT(TMP)-EXPONENT(ExtState%O3%Arr%Val(I,J,1))) &
                     < MAXEXPONENT(TMP) ) THEN
                   DEPO3(I,J) = TMP / ExtState%O3%Arr%Val(I,J,1)
                ENDIF

                ! Sanity check: if deposition velocities are above one, 
                ! something must have gone wrong (they are on the order
                ! of <1e-9)
                IF ( DEPO3(I,J) > 1.0_hp ) THEN
                   DEPO3(I,J) = 0.0_hp
                   WRITE(MSG,*) 'O3 deposition velocity > 1., set to zero', &
                      I, J, DEPO3(I,J), ABS(iFlx), ExtState%O3%Arr%Val(I,J,1)
                   CALL HCO_WARNING(MSG, RC)
                ENDIF
             ENDIF

          ENDIF
       ENDIF

       ! Eventually write out into diagnostics array
       IF ( DoDiagn ) THEN
          DIAGN(I,J,1) = FRACTION_NOx
          DIAGN(I,J,2) = INT_OPE
          DIAGN(I,J,3) = ShipNoEmis(I,J,1) * ( 1.0d0-FRACTION_NOx) * INT_OPE 
          DIAGN(I,J,4) = ShipNoEmis(I,J,1)
       ENDIF

       ! Reset ship NO emissions to zero. Will be refilled on next
       ! emission step!
       ShipNoEmis(I,J,1) = 0.0d0

    ENDDO !I
    ENDDO !J
#if defined( NULL )
!$OMP END PARALLEL DO
#endif

    ! Error check
    IF ( ERR ) THEN
       RC = HCO_FAIL
       RETURN 
    ENDIF

    !=======================================================================
    ! PASS TO HEMCO STATE AND UPDATE DIAGNOSTICS 
    !=======================================================================

    ! NO
    IF ( IDTNO > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( HcoState, FLUXNO, IDTNO, RC)
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXNO', RC )
          RETURN 
       ENDIF

       ! Eventually update diagnostics
       IF ( Diagn_AutoFillLevelDefined(2) ) THEN
          Arr2D => FLUXNO
          CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNr, &
                             Cat=-1, Hier=-1, HcoID=IDTNO,     &
                             AutoFill=1, Array2D=Arr2D, RC=RC   )
          IF ( RC /= HCO_SUCCESS ) RETURN 
          Arr2D => NULL() 
       ENDIF
    ENDIF

    ! HNO3 
    IF ( IDTHNO3 > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( HcoState, FLUXHNO3, IDTHNO3, RC)
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXHNO3', RC )
          RETURN 
       ENDIF

       ! Eventually update diagnostics
       IF ( Diagn_AutoFillLevelDefined(2) ) THEN
          Arr2D => FLUXHNO3
          CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNr, &
                             Cat=-1, Hier=-1, HcoID=IDTHNO3,   &
                             AutoFill=1, Array2D=Arr2D, RC=RC   )
          IF ( RC /= HCO_SUCCESS ) RETURN 
          Arr2D => NULL() 
       ENDIF
    ENDIF

    ! O3 
    IF ( IDTO3 > 0 ) THEN

       ! Add flux to emission array (kg/m2/s)
       CALL HCO_EmisAdd( HcoState, FLUXO3, IDTO3, RC)
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXO3', RC )
          RETURN 
       ENDIF

       ! Eventually update diagnostics
       IF ( Diagn_AutoFillLevelDefined(2) ) THEN
          Arr2D => FLUXO3
          CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNr, &
                             Cat=-1, Hier=-1, HcoID=IDTO3,   &
                             AutoFill=1, Array2D=Arr2D, RC=RC   )
          IF ( RC /= HCO_SUCCESS ) RETURN 
          Arr2D => NULL() 
       ENDIF

       ! Add flux to emission array (1/s)
       CALL HCO_DepvAdd( HcoState, DEPO3, IDTO3, RC)
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! TODO: Add deposition diagnostics
    ENDIF


    ! Eventually update manual diagnostics
    IF ( DoDiagn ) THEN
       DiagnName =  'PARANOX_NOXFRAC_REMAINING'
       Arr2D     => DIAGN(:,:,1)
       CALL Diagn_Update( am_I_Root, HcoState,   ExtNr=ExtNr, &
                          cName=TRIM(DiagnName), Array2D=Arr2D, RC=RC)
       IF ( RC /= HCO_SUCCESS ) RETURN
       Arr2D => NULL()       
       
       DiagnName =  'PARANOX_OPE'
       Arr2D     => DIAGN(:,:,2)
       CALL Diagn_Update( am_I_Root, HcoState,   ExtNr=ExtNr, &
                          cName=TRIM(DiagnName), Array2D=Arr2D, RC=RC)
       IF ( RC /= HCO_SUCCESS ) RETURN
       Arr2D => NULL()       

       DiagnName =  'PARANOX_O3_PRODUCTION'
       Arr2D     => DIAGN(:,:,3)
       CALL Diagn_Update( am_I_Root, HcoState,   ExtNr=ExtNr, &
                          cName=TRIM(DiagnName), Array2D=Arr2D, RC=RC)
       IF ( RC /= HCO_SUCCESS ) RETURN
       Arr2D => NULL()       

       DiagnName =  'PARANOX_TOTAL_SHIPNOX'
       Arr2D     => DIAGN(:,:,4)
       CALL Diagn_Update( am_I_Root, HcoState,   ExtNr=ExtNr, &
                          cName=TRIM(DiagnName), Array2D=Arr2D, RC=RC)
       IF ( RC /= HCO_SUCCESS ) RETURN
       Arr2D => NULL()       
    ENDIF

    ! Update last hour to current one for next call.
    lastHH = HH

    ! Return w/ success
    CALL HCO_LEAVE ( RC )

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
 SUBROUTINE HCOX_ParaNOx_Init( am_I_Root, HcoState, ExtName, ExtState, RC ) 
!
! !USES:
!
   USE HCO_Chartools_Mod, ONLY : HCO_CharParse
   USE HCO_State_MOD,     ONLY : HCO_GetHcoID
   USE HCO_State_MOD,     ONLY : HCO_GetExtHcoID
   USE HCO_ExtList_Mod,   ONLY : GetExtNr
   USE HCO_ExtList_Mod,   ONLY : GetExtOpt
   USE ParaNOx_Util_Mod,  ONLY : Read_ParaNOx_LUT
!
! !INPUT PARAMETERS:
!
   LOGICAL,          INTENT(IN   )  :: am_I_Root
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER                        :: AS, I, tmpID, IDTNO2, nSpc
   INTEGER,           ALLOCATABLE :: HcoIDs(:)
   CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
   CHARACTER(LEN=255)             :: MSG, LOC

   !========================================================================
   ! HCOX_PARANOX_INIT begins here!
   !========================================================================

   ! Extension Nr.
   ExtNr = GetExtNr( TRIM(ExtName) )
   IF ( ExtNr <= 0 ) RETURN

   ! Enter
   CALL HCO_ENTER( 'HCOX_ParaNOx_Init (hcox_paranox_mod.F90)', RC )
   IF ( RC /= HCO_SUCCESS ) RETURN

   !------------------------------------------------------------------------
   ! Get species IDs
   !------------------------------------------------------------------------

   ! Get HEMCO species IDs
   CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
   IF ( RC /= HCO_SUCCESS ) RETURN

   ! Init species ID
   IDTNO   = -1
   IDTNO2  = -1
   IDTO3   = -1
   IDTHNO3 = -1

   ! Check for NO, NO2, O3, and HNO3
   DO I = 1, nSpc
      SELECT CASE ( TRIM(SpcNames(I)) )
         CASE ( "NO" )
            IDTNO = HcoIDs(I)
         CASE ( "NO2" )
            IDTNO2 = HcoIDs(I)
         CASE ( "O3" )
            IDTO3 = HcoIDs(I)
         CASE ( "HNO3" )
            IDTHNO3 = HcoIDs(I)
         CASE DEFAULT
            ! leave empty 
      END SELECT
   ENDDO

   ! NO must be defined
   IF ( IDTNO <= 0 ) THEN
      MSG = 'Species NO not defined in PARANOX - this is unrealistic!' 
      CALL HCO_ERROR ( MSG, RC )
      RETURN
   ENDIF
   MW_NO = HcoState%Spc(IDTNO)%MW_g

   ! Get MW of other three species. If species not set in the configuration
   ! file (e.g. they are not being used by PARANOx), determine MW from 
   ! default values.

   ! O3
   IF ( IDTO3 <= 0 ) THEN
      tmpID = HCO_GetHcoID('O3', HcoState )
      MSG = 'O3 not produced/removed in PARANOX'
      CALL HCO_WARNING ( MSG, RC )
   ELSE
      tmpID = IDTO3
   ENDIF
   IF ( tmpID > 0 ) THEN
      MW_O3 = HcoState%Spc(tmpID)%MW_g
   ELSE
      MSG = 'Use default O3 molecular weight of 48g/mol'
      CALL HCO_WARNING ( MSG, RC )
      MW_O3 = 48.0_dp
   ENDIF
   
   ! NO2
   IF ( IDTNO2 <= 0 ) THEN
      tmpID = HCO_GetHcoID('NO2', HcoState )
   ELSE
      tmpID = IDTNO2
   ENDIF
   IF ( tmpID > 0 ) THEN
      MW_NO2 = HcoState%Spc(tmpID)%MW_g
   ELSE
      MSG = 'Use default NO2 molecular weight of 46g/mol'
      CALL HCO_WARNING ( MSG, RC )
      MW_NO2 = 46.0_dp
   ENDIF
   
   ! HNO3
   IF ( IDTHNO3 <= 0 ) THEN
      tmpID = HCO_GetHcoID('HNO3', HcoState )
      MSG = 'HNO3 not produced/removed in PARANOX'
      CALL HCO_WARNING ( MSG, RC )
   ELSE
      tmpID = IDTHNO3
   ENDIF
   IF ( tmpID > 0 ) THEN
      MW_HNO3 = HcoState%Spc(tmpID)%MW_g
   ELSE
      MSG = 'Use default HNO3 molecular weight of 63g/mol'
      CALL HCO_WARNING ( MSG, RC )
      MW_HNO3 = 63.0_dp
   ENDIF

   ! Verbose mode
   IF ( am_I_Root ) THEN
      MSG = 'Use ParaNOx ship emissions (extension module)'
      CALL HCO_MSG( MSG, SEP1='-' )
      MSG = '    - Use the following species: (MW, emitted as HEMCO ID) ' 
      CALL HCO_MSG( MSG )
      WRITE(MSG,"(a,F5.2,I5)") '     NO  : ', MW_NO, IDTNO
      CALL HCO_MSG(MSG)
      WRITE(MSG,"(a,F5.2,I5)") '     O3  : ', MW_O3, IDTO3
      CALL HCO_MSG(MSG)
      WRITE(MSG,"(a,F5.2,I5)") '     HNO3: ', MW_HNO3, IDTHNO3
      CALL HCO_MSG(MSG)
      WRITE(MSG,"(a,F5.2,a)" ) '     NO2 : ', MW_NO2, ' - '
      CALL HCO_MSG(MSG)
   ENDIF

   !------------------------------------------------------------------------
   ! Initialize the PARANOX look-up tables
   !------------------------------------------------------------------------

   ! Fraction of NOx remaining for ship emissions
   CALL GetExtOpt ( ExtNr, 'FracNOx table', OptValChar=FRACNOX_FILE, RC=RC)
   IF ( RC /= HCO_SUCCESS ) RETURN

   ! Integrated Ozone production efficiency (OPE)
   CALL GetExtOpt ( ExtNr, 'IntOPE table', OptValChar=INTOPE_FILE, RC=RC)
   IF ( RC /= HCO_SUCCESS ) RETURN

   ! Call HEMCO parser to replace tokens such as $ROOT, $MET, or $RES.
   ! There shouldn't be any date token in there ($YYYY, etc.), so just
   ! provide some dummy variables here
   CALL HCO_CharParse( FRACNOX_FILE, -999, -1, -1, -1, RC )
   IF ( RC /= HCO_SUCCESS ) RETURN

   CALL HCO_CharParse( INTOPE_FILE, -999, -1, -1, -1, RC )
   IF ( RC /= HCO_SUCCESS ) RETURN

   ! Read PARANOX look-up tables from disk
   ! NOTE: Currently these are read from binary file, which is incompatible
   ! with the ESMF/MAPL run environment.  We are currently working on a
   ! better implementation of this, stay tuned. (bmy, 8/13/14)
   CALL READ_PARANOX_LUT( am_I_Root, FracNOx_FILE, IntOPE_FILE, RC )
   IF ( RC /= HCO_SUCCESS ) RETURN

   !------------------------------------------------------------------------
   ! Set other module variables 
   !------------------------------------------------------------------------ 
   ALLOCATE ( ShipNO(HcoState%NX,HcoState%NY,HcoState%NZ), STAT=AS )
   IF ( AS /= 0 ) THEN
      CALL HCO_ERROR ( 'ShipNO', RC )
      RETURN
   ENDIF
   ShipNO = 0.0_hp

   ! Allocate variables for SunCosMid from 5 hours ago. We internally store 
   ! the SunCosMid values from the past 5 hours in array SC5 and cycle through 
   ! that array to get the SUNCOS value from 5 hours ago. The variable SC5ID 
   ! is used to identify the slice representing the values from 5 hours ago for
   ! the given time. It is updated every time the simulation hour changes. 
   ! The sixth slice acts as a 'buffer' that holds the SUNCOS value from the 
   ! last time step.
   ! We assume explicitly that the chemistry time step is not larger that 60 
   ! mins, i.e. that PARANOX is called at least once per hour. If that's not
   ! the case, the SC5 array will hold values from further back!
   ALLOCATE ( SC5(HcoState%NX,HcoState%NY,6), STAT=AS )
   IF ( AS /= 0 ) THEN
      CALL HCO_ERROR ( 'SC5', RC )
      RETURN
   ENDIF
   SC5   = 0.0_hp
   SC5ID = 1

   ! Prompt warning if chemistry time step is more than 60 mins
   IF ( HcoState%TS_CHEM > 3600.0_hp ) THEN
      IF ( am_I_Root ) THEN
         MSG = '    Cannot properly store SUNCOSmid values ' // &
               ' because chemistry time step is more than 60 mins!'
         CALL HCO_WARNING ( MSG, RC )
      ENDIF
   ENDIF

   ! Molecular weight of AIR
   MW_AIR = HcoState%Phys%AIRMW

   ! Met. data required by module
   ExtState%O3%DoUse         = .TRUE.
   ExtState%NO2%DoUse        = .TRUE.
   ExtState%NO%DoUse         = .TRUE.
   ExtState%AIR%DoUse        = .TRUE.
   ExtState%SUNCOSmid%DoUse  = .TRUE.
   ExtState%T2M%DoUse        = .TRUE.
   IF ( IDTHNO3 > 0 ) THEN
      ExtState%HNO3%DoUse    = .TRUE.
   ENDIF
   ExtState%JNO2%DoUse       = .TRUE.
   ExtState%JO1D%DoUse       = .TRUE.

   ! Enable module
   ExtState%ParaNOx = .TRUE.

   !------------------------------------------------------------------------
   ! Leave w/ success
   !------------------------------------------------------------------------
   IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  )
   IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
   CALL HCO_LEAVE ( RC )

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
 SUBROUTINE HCOX_ParaNOx_Final()
!
! !REVISION HISTORY:
!  06 Aug 2013 - C. Keller - Initial Version
!  06 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  06 Jun 2014 - R. Yantosca - Now indended with F90 free-format
!EOP
!------------------------------------------------------------------------------
!BOC

   !=================================================================
   ! HCOX_PARANOX_FINAL begins here!
   !=================================================================

   IF ( ALLOCATED(ShipNO) ) DEALLOCATE ( ShipNO )
   IF ( ALLOCATED(SC5   ) ) DEALLOCATE ( SC5    )

 END SUBROUTINE HCOX_ParaNOx_Final
!EOC
END MODULE HCOX_ParaNOx_mod
