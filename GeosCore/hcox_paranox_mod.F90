!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
  PUBLIC :: HCOX_ParaNOx_Run
  PUBLIC :: HCOX_ParaNOx_Init
  PUBLIC :: HCOX_ParaNOx_Final
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
  REAL*8                        :: MW_AIR

  ! Arrays
  REAL(hp), ALLOCATABLE, TARGET :: ShipNO(:,:,:)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
    USE HCO_DRIVER_MOD,     ONLY : HCO_RUN
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root
    TYPE(Ext_State), POINTER        :: ExtState     ! Module options
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState   ! Hemco state 
    INTEGER,         INTENT(INOUT)  :: RC 

! !REVISION HISTORY:
!  06 Aug 2013 - C. Keller   - Initial Version
!  06 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  06 Jun 2014 - R. Yantosca - Now indended with F90 free-format
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
    CALL HCO_ENTER( 'HCOX_ParaNOx_Run (hcox_paranox_mod.F90)', RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! ----------------------------------------------------------------
    ! Use HEMCO core routines to get ship NO emissions 
    ! ----------------------------------------------------------------

    ! Prepare HEMCO core run:
    ! --> Set tracer and category range + extension number.
    ! Note: Set species min and max to the full range of species. 
    ! For the ParaNox extension, emission fields of only one species
    ! should be defined. HCO_RUN will exit w/ error if this is not 
    ! the case. 
    HcoState%Options%SpcMin =  1 
    HcoState%Options%SpcMax = -1 
    HcoState%Options%CatMin =  1 
    HcoState%Options%CatMax = -1
    HcoState%Options%ExtNr  = ExtNr 

      ! --> Define array to write emissions into. ShipNO is reset to
      ! zero within subroutine EVOLVE_PLUME, so no need to do this
      ! here.
!      ShipNO                      = 0.0d0
    HcoState%Options%AutoFillDiagn = .FALSE.
    HcoState%Options%FillBuffer    =  .TRUE.
    HcoState%Buffer3D%Val          => ShipNO 
      
    ! Calculate ship NO emissions and write them into the ShipNO
    ! array [kg/m2/s]. 
    CALL HCO_RUN ( am_I_Root, HcoState, RC )
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
    CALL EVOLVE_PLUME ( am_I_Root, ExtState, ShipNO, HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Leave w/ success
    CALL HCO_LEAVE ( RC ) 

  END SUBROUTINE HCOX_ParaNOx_Run
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Evolve_Plume 
!
! !DESCRIPTION: Subroutine EVOLVE\_PLUME performs plume
! dilution/chemistry of ship NO emissions for every grid box and writes
! the resulting NO, HNO3 and O3 emission (production) rates into
! State\_Chm%NomixS. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Evolve_Plume( am_I_Root, ExtState, ShipNoEmis, HcoState, RC )  
!
! !USES:
!
    USE ParaNOx_Util_Mod,  ONLY : INTERPOLATE_LUT2
    USE ParaNOx_Util_Mod,  ONLY : READ_PARANOX_LUT
    USE Comode_Loop_Mod,   ONLY : NCS, JPHOTRAT, NRATES
    USE Comode_loop_mod,   ONLY : NAMEGAS, IRM
    USE HCO_FluxArr_mod,   ONLY : HCO_EmisAdd, HCO_DepvAdd
    USE HCO_Clock_Mod,     ONLY : HcoClock_Get
    USE FAST_JX_MOD,       ONLY : FJXFUNC

    ! testing only
    USE PBL_MIX_MOD,       ONLY : GET_PBL_TOP_L
    USE PRESSURE_MOD,      ONLY : GET_PEDGE
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root
    TYPE(Ext_State), POINTER        :: ExtState       ! Module options
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(hp),        INTENT(INOUT)  :: ShipNoEmis(:,:,:)
    TYPE(HCO_State), POINTER        :: HcoState     ! Output obj
    INTEGER,         INTENT(INOUT)  :: RC 
!
! !REVISION HISTORY:
!  06 Aug 2013 - C. Keller   - Initial Version
!  06 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  06 Jun 2014 - R. Yantosca - Now indended with F90 free-format
!  24 Jun 2014 - R. Yantosca - Now pass LUT_FILENAME to READ_PARANOX_LUT
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L, MM
    INTEGER            :: NK, LMAX
    INTEGER, SAVE      :: SAVEMM = -999
    LOGICAL            :: ERR
    REAL*8             :: JNO2, JO1D, TS, SUNCOSmid5, SUNCOSmid
    REAL*8             :: O3molec, NOmolec, NO2molec, AIRmolec
    REAL*4             :: FRACTION_NOx, INT_OPE
    CHARACTER(LEN=8)   :: SPECNAME
    CHARACTER(LEN=255) :: MSG
    REAL(hp), TARGET   :: FLUXNO(HcoState%NX,HcoState%NY)
    REAL(hp), TARGET   :: FLUXHNO3(HcoState%NX,HcoState%NY)
    REAL(hp), TARGET   :: FLUXO3(HcoState%NX,HcoState%NY)
    REAL(hp), TARGET   :: DEPO3(HcoState%NX,HcoState%NY)
    REAL(hp)           :: iFlx
    REAL(hp), POINTER  :: Arr2D(:,:) => NULL()

    ! testing only
    REAL*8  :: FRAC, TOTPRES, DELTPRES
    INTEGER :: TOP
    integer :: ix, jx
    logical, parameter :: add2hemco = .true.

    !=================================================================
    ! EVOLVE_PLUME begins here!
    !=================================================================

    ! testing only
    ix = -1 !71
    jx = -1 !35

    ! Enter
    CALL HCO_ENTER( 'Evolve_Plume (hcox_paranox_mod.F90)', RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Leave here if none of the tracers defined
    IF ( IDTNO <= 0 .AND. IDTO3 <= 0 .AND. IDTHNO3 <= 0 ) THEN
       RC = HCO_SUCCESS 
       RETURN
    ENDIF

    ! Get simulation month
    CALL HcoClock_Get( cMM=MM, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Read look up tables every new month
    IF ( MM /= SAVEMM ) THEN
       CALL READ_PARANOX_LUT( FracNOx_FILE, IntOPE_FILE  )
       SAVEMM = MM
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
#if defined( 0 )
!$OMP PARALLEL DO                                                   &
!$OMP DEFAULT( SHARED )                                             &
!$OMP PRIVATE( I, J, L, RC, iFlx, NK, SPECNAME, JNO2, JO1D        ) &
!$OMP PRIVATE( O3molec, NOmolec, NO2molec, AIRmolec, INT_OPE      ) &
!$OMP PRIVATE( FRACTION_NOx, LMAX, TS, SUNCOSmid5, SUNCOSmid      ) &     
!$OMP SCHEDULE( DYNAMIC )
#endif
    DO J = 1, HcoState%NY
    DO I = 1, HcoState%NX

       ! Skip if no ship emissions in this grid box
       IF ( ShipNoEmis(I,J,1) == 0d0 ) CYCLE

       !--------------------------------------------------------------------
       ! Get J-Values for J(NO2) and J(O3)
       !--------------------------------------------------------------------

       ! Check if sun is up
       IF ( ExtState%SUNCOSmid%Arr%Val(I,J) > 0d0 ) THEN

          ! Loop over photolysis reactions to find NO2, O3
          LMAX = JPHOTRAT(NCS)
          DO L = 1, LMAX 

             ! Reaction number
             NK  = NRATES(NCS) + L

             ! Name of species being photolyzed
             SPECNAME = NAMEGAS(IRM(1,NK,NCS))

             ! Check if this is NO2 or O3, store values, 1/s
             SELECT CASE ( TRIM( SPECNAME ) )
                CASE ( 'NO2' )
                   JNO2 = FJXFUNC(I,J,1,L,1,SPECNAME)
                CASE ( 'O3' )
#if defined( UCX )      
                   ! IMPORTANT: Need branch *2* for O1D
                   ! Branch 1 is O3P!
                   JO1D = FJXFUNC(I,J,1,L,2,SPECNAME)
#else
                   JO1D = FJXFUNC(I,J,1,L,1,SPECNAME)
#endif
                CASE DEFAULT
             END SELECT

          ENDDO

       ELSE

          ! J-values are zero when sun is down
          JNO2 = 0d0
          JO1D = 0d0
       ENDIF

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
       TS         = ExtState%TSURFK%Arr%Val(I,J)
       SUNCOSmid5 = ExtState%SUNCOSmid5%Arr%Val(I,J)
       SUNCOSmid  = ExtState%SUNCOSmid%Arr%Val(I,J)

       CALL INTERPOLATE_LUT2( I,            J,            &
                              O3molec,      NOmolec,      &
                              NO2molec,     AIRmolec,     &
                              JO1D,         JNO2,         &
                              TS,           SUNCOSmid5,   &
                              SUNCOSmid,                  &
                              FRACTION_NOx, INT_OPE )

       ! testing only
       if ( i==ix .and. j==jx ) then
          write(*,*) 'O3molec: ', O3molec
          write(*,*) 'NOmolec: ', NOmolec
          write(*,*) 'NO2molec: ', NO2molec
          write(*,*) 'AIRmolec: ', AIRmolec
          write(*,*) 'JO1D    : ', JO1D 
          write(*,*) 'JNO2    : ', JNO2 
          write(*,*) 'SUNCOSmid5:' , SUNCOSmid5
          write(*,*) 'SUNCOSmid :' , SUNCOSmid 
          write(*,*) 'TS        :', TS
          write(*,*) 'FRACTION_NOX: ', FRACTION_NOX
          write(*,*) 'INT_OPE: ', Int_OPE
          write(*,*) 'ShipEmis [kg/m2/s]: ', ShipNoEmis(I,J,1)
       endif

       !---------------------------
       ! Calculate NO emissions
       !---------------------------
       IF ( IDTNO > 0 ) THEN
          
          ! Of the total ship NOx, the fraction FRACTION_NOX
          ! survives after plume dilution and chemistry.
          ! Unit: kg/m2/s 
          FLUXNO(I,J) = ShipNoEmis(I,J,1) * FRACTION_NOx

!            ! Add to flux array
!            if ( add2hemco ) then
!            CALL HCO_EmisAdd ( HcoState, FLUXNO, IDTNO, I, J, 1, RC ) 
!            IF ( RC /= HCO_SUCCESS ) THEN
!               ERR = .TRUE.
!               EXIT
!            ENDIF 
!            endif

          ! testing only
          if ( i==ix .and. j==jx ) then
             write(*,*) 'FLUXNO [kg/m2/s]: ', FLUXNO(I,J)
          endif

       ENDIF

       !---------------------------
       ! Calculate HNO3 emissions
       !---------------------------
       IF ( IDTHNO3 > 0 ) THEN

          ! Of the total ship NOx, the fraction 1-FRACTION_NOX
          ! is converted to HNO3 during plume dilution and chemistry. 
          ! Unit: kg/m2/s 
          FLUXHNO3(I,J) = ShipNoEmis(I,J,1) * ( 1d0 - FRACTION_NOx )

!            ! Add to flux array 
!            if ( add2hemco ) then
!            CALL HCO_EmisAdd ( HcoState, FLUXHNO3, IDTHNO3, I,J,1, RC ) 
!            IF ( RC /= HCO_SUCCESS ) THEN
!               ERR = .TRUE.
!               EXIT
!            ENDIF 
!            endif

          ! testing only
          if ( i==ix .and. j==jx ) then
             write(*,*) 'FLUXHNO3 [kg/m2/s]: ', FLUXHNO3(I,J) 
          endif

       ENDIF

       !---------------------------
       ! Calculate O3 emissions
       !---------------------------
       IF ( IDTO3 > 0 ) THEN
          
          ! Of the total ship NOx, the fraction
          ! (1-FRACTION_NOX)*INT_OPE is converted to O3 during 
          ! plume dilution and chemistry. 
          ! Unit: kg/m2/s 
          iFlx = ShipNoEmis(I,J,1) * (1d0-FRACTION_NOx) * INT_OPE

          ! testing only
          if ( i==ix .and. j==jx ) then
             write(*,*) 'FLUXO3 [kg/m2/s]: ', iFlx 
          endif

          ! For positive fluxes, add to emission flux array 
!            if ( add2hemco ) then
          IF ( iFlx >= 0d0 ) THEN
             FLUXO3(I,J) = iFlx

!               CALL HCO_EmisAdd ( HcoState, FLUXO3, IDTO3, I,J,1, RC )
!               IF ( RC /= HCO_SUCCESS ) THEN
!                  ERR = .TRUE.
!                  EXIT
!               ENDIF 

          ! For negative fluxes, calculate deposition velocity based
          ! on current surface O3 concentration and pass to deposition
          ! array
          ELSE

             ! Calculate deposition velocity (m/s) from flux
             ! NOTE: the calculated deposition flux is in kg/m2/s,
             ! which has to be converted to m/s. Use here the O3 conc.
             ! [kg/m3] of the lowest model box. 
             DEPO3(I,J) = ABS(iFlx) / ExtState%O3%Arr%Val(I,J,1) &
                          * ExtState%AIRVOL%Arr%Val(I,J,1)

             ! testing only
             if ( i==ix .and. j==jx ) then
                write(*,*) 'O3 deposition [m/s]: ', DEPO3(I,J) 
             endif

!               CALL HCO_DepvAdd ( HcoState, FLUXO3, IDTO3, I,J, RC )
!               IF ( RC /= HCO_SUCCESS ) THEN
!                  ERR = .TRUE.
!                  EXIT
!               ENDIF 

          ENDIF
!        endif  
       ENDIF

       ! Reset ship NO emissions to zero. Will be refilled on next
       ! emission step!
       ShipNoEmis(I,J,1) = 0.0d0

    ENDDO !I
    ENDDO !J
#if defined( 0 )
!$OMP END PARALLEL DO
#endif

    ! Error check
    IF ( ERR ) THEN
       RC = HCO_FAIL
       RETURN 
    ENDIF

    !=================================================================
    ! PASS TO HEMCO STATE AND UPDATE DIAGNOSTICS 
    !=================================================================
    if ( add2hemco ) then

       ! NO
       IF ( IDTNO > 0 ) THEN

          ! Add flux to emission array
          CALL HCO_EmisAdd( HcoState, FLUXNO, IDTNO, RC)
          IF ( RC /= HCO_SUCCESS ) RETURN 

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
          IF ( RC /= HCO_SUCCESS ) RETURN 

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

          ! Add flux to emission array
          CALL HCO_EmisAdd( HcoState, FLUXO3, IDTO3, RC)
          IF ( RC /= HCO_SUCCESS ) RETURN 

          ! Eventually update diagnostics
          IF ( Diagn_AutoFillLevelDefined(2) ) THEN
             Arr2D => FLUXO3
             CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNr, &
                                Cat=-1, Hier=-1, HcoID=IDTO3,   &
                                AutoFill=1, Array2D=Arr2D, RC=RC   )
             IF ( RC /= HCO_SUCCESS ) RETURN 
             Arr2D => NULL() 
          ENDIF

          ! Add flux to emission array
          CALL HCO_DepvAdd( HcoState, DEPO3, IDTO3, RC)
          IF ( RC /= HCO_SUCCESS ) RETURN 

          ! TODO: Add deposition diagnostics

       ENDIF

       ! add to tracer array directly (testing only)
    else

       DO J=1,HcoState%NY
       DO I=1,HcoState%NX

          ! Top level of the boundary layer
          ! guard for b.l. being in first level.
          TOP = FLOOR( GET_PBL_TOP_L( I, J ) )
          IF ( TOP == 0 ) TOP = 1
          
          ! Pressure thickness of entire boundary layer [hPa]
          TOTPRES = GET_PEDGE(I,J,1) - GET_PEDGE(I,J,TOP+1)

          ! Add emissions to tracer arrays. Mix over entire boundary
          ! layer
          DO L = 1, TOP

             ! Thickness of level L [mb]
             DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)
             FRAC     = DELTPRES / TOTPRES

             ! testing only
             if ( i==ix .and. j==jx ) then 
                write(*,*) 'Level ', L, ':'
                write(*,*) 'FRAC: ', FRAC
             endif

             ! Add proportion of mass to this level [kg]
             IF ( IDTNO > 0 ) THEN

                ! testing only
                if ( i==ix .and. j==jx ) then 
                   write(*,*) 'NO before [kg]: ', ExtState%NO%Arr%Val(I,J,L)
                endif

                ExtState%NO%Arr%Val(I,J,L) = &
                     ExtState%NO%Arr%Val(I,J,L) + &
                     ( FLUXNO(I,J) * HcoState%Grid%AREA_M2(I,J) * &
                     HcoState%TS_EMIS * FRAC ) 

                ! testing only
                if ( i==ix .and. j==jx ) then 
                   write(*,*) 'NO after [kg]: ', ExtState%NO%Arr%Val(I,J,L)
                endif
             ENDIF

             IF ( IDTHNO3 > 0 ) THEN
                
                ! testing only
                if ( i==ix .and. j==jx ) then 
                   write(*,*) 'HNO3 before [kg]: ', ExtState%HNO3%Arr%Val(I,J,L)
                endif

                ExtState%HNO3%Arr%Val(I,J,L) = &
                     ExtState%HNO3%Arr%Val(I,J,L) + &
                     ( FLUXHNO3(I,J) * HcoState%Grid%AREA_M2(I,J) * &
                     HcoState%TS_EMIS * FRAC ) 

                ! testing only
                if ( i==ix .and. j==jx ) then 
                   write(*,*) 'HNO3 after [kg]: ', ExtState%HNO3%Arr%Val(I,J,L)
                endif

             ENDIF

             IF ( IDTO3 > 0 ) THEN

                ! testing only
                if ( i==ix .and. j==jx ) then 
                   write(*,*) 'O3 before [kg]: ', ExtState%O3%Arr%Val(I,J,L)
                endif

                ExtState%O3%Arr%Val(I,J,L) = &
                     ExtState%O3%Arr%Val(I,J,L) + &
                     ( FLUXO3(I,J) * HcoState%Grid%AREA_M2(I,J) * &
                     HcoState%TS_EMIS * FRAC ) 

                iFlx = DEPO3(I,J) * ExtState%O3%Arr%Val(I,J,1) &
                     / ExtState%AIRVOL%Arr%Val(I,J,1)
                ExtState%O3%Arr%Val(I,J,L) = &
                     ExtState%O3%Arr%Val(I,J,L) - &
                     ( iFlx * HcoState%Grid%AREA_M2(I,J) * &
                     HcoState%TS_EMIS * FRAC ) 

                ! Make sure O3 is not negative!
                IF ( ExtState%O3%Arr%Val(I,J,L) < 0d0 ) THEN
                   ExtState%O3%Arr%Val(I,J,L) = 0d0
                ENDIF

                ! testing only
                if ( i==ix .and. j==jx ) then 
                   write(*,*) 'O3 after [kg]: ', ExtState%O3%Arr%Val(I,J,L)
                endif

             ENDIF

          ENDDO !L
      ENDDO !I
      ENDDO !J

   endif ! add2hemco

   ! Return w/ success
   CALL HCO_LEAVE ( RC )

 END SUBROUTINE Evolve_Plume
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
 SUBROUTINE HCOX_ParaNOx_Init ( am_I_Root, HcoState, ExtName, &
                                ExtState,    RC                  ) 
!
! !USES:
!
   USE HCO_State_MOD,     ONLY : HCO_GetHcoID
   USE HCOX_ExtList_Mod,  ONLY : GetExtNr
   USE HCOX_ExtList_Mod,  ONLY : GetExtHcoID
   USE HCOX_ExtList_Mod,  ONLY : GetExtOpt
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER                        :: AS, IDTNO2, nSpc
   INTEGER,           ALLOCATABLE :: HcoIDs(:)
   CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
   CHARACTER(LEN=255)             :: MSG, LOC

   !=================================================================
   ! HCOX_PARANOX_INIT begins here!
   !=================================================================

   ! Extension Nr.
   ExtNr = GetExtNr( TRIM(ExtName) )
   IF ( ExtNr <= 0 ) RETURN

   ! Enter
   CALL HCO_ENTER( 'HCOX_ParaNOx_Init (hcox_paranox_mod.F90)', RC )
   IF ( RC /= HCO_SUCCESS ) RETURN

   ! ---------------------------------------------------------------------- 
   ! Get the look-up-table filenames
   ! ---------------------------------------------------------------------- 

   ! Fraction of NOx remaining for ship emissions
   CALL GetExtOpt ( ExtNr, 'FracNOx table', OptValChar=FRACNOX_FILE, RC=RC)
   IF ( RC /= HCO_SUCCESS ) RETURN

   ! Integrated Ozone production efficiency (OPE)
   CALL GetExtOpt ( ExtNr, 'IntOPE table', OptValChar=INTOPE_FILE, RC=RC)
   IF ( RC /= HCO_SUCCESS ) RETURN

   ! ---------------------------------------------------------------------- 
   ! Get species IDs
   ! ---------------------------------------------------------------------- 

   ! Get HEMCO species IDs
   CALL GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
   IF ( RC /= HCO_SUCCESS ) RETURN
   IF ( nSpc /= 4 ) THEN
      MSG = 'Four species (NO/NO2/O3/HNO3) expected for ParaNOx!'
      CALL HCO_ERROR ( MSG, RC )
      RETURN
   ENDIF
   IDTNO   = HcoIDs(1)
   IDTNO2  = HcoIDs(2)
   IDTO3   = HcoIDs(3)
   IDTHNO3 = HcoIDs(4)

   IF ( IDTNO <= 0 ) THEN
      MSG = 'Species NO not defined - needed by ParaNOx!' 
      CALL HCO_ERROR ( MSG, RC )
      RETURN
   ELSE
      MW_NO = HcoState%Spc(IDTNO)%MW_g
   ENDIF

   IDTO3 = HCO_GetHcoID('O3',   HcoState )
   IF ( IDTO3 <= 0 ) THEN
      MSG = 'Species O3 not defined - needed by ParaNOx!' 
      CALL HCO_ERROR ( MSG, RC )
      RETURN
   ELSE
      MW_O3 = HcoState%Spc(IDTO3)%MW_g
   ENDIF

   IDTHNO3 = HCO_GetHcoID('HNO3', HcoState )
   IF ( IDTHNO3 <= 0 ) THEN
      MSG = 'Species HNO3 not defined - needed by ParaNOx!' 
      CALL HCO_WARNING ( MSG, RC )
   ENDIF

   IDTNO2  = HCO_GetHcoID('NO2',  HcoState )
   IF ( IDTNO2 <= 0 ) THEN
      MSG = 'Species NO2 not defined - needed by ParaNOx!' 
      CALL HCO_ERROR ( MSG, RC )
      RETURN
   ELSE
      MW_NO2 = HcoState%Spc(IDTNO2)%MW_g
   ENDIF

   ! Verbose mode
   MSG = 'Use ParaNOx ship emissions (extension module)'
   CALL HCO_MSG( MSG )
   MSG = '    - Use the following species: ' 
   CALL HCO_MSG( MSG )
   WRITE(MSG,*) '     NO  : ', TRIM(SpcNames(1)), IDTNO
   CALL HCO_MSG(MSG)
   WRITE(MSG,*) '     NO2 : ', TRIM(SpcNames(2)), IDTNO2
   CALL HCO_MSG(MSG)
   WRITE(MSG,*) '     O3  : ', TRIM(SpcNames(3)), IDTO3
   CALL HCO_MSG(MSG)
   WRITE(MSG,*) '     HNO3: ', TRIM(SpcNames(4)), IDTHNO3
   CALL HCO_MSG(MSG)

   ! ---------------------------------------------------------------------- 
   ! Set other module variables 
   ! ---------------------------------------------------------------------- 
   ALLOCATE ( ShipNO(HcoState%NX,HcoState%NY,HcoState%NZ), STAT=AS )
   IF ( AS /= 0 ) THEN
      CALL HCO_ERROR ( 'ShipNO', RC )
      RETURN
   ENDIF
   ShipNO = 0.0_hp

   ! Molecular weight of AIR
   MW_AIR = HcoState%Phys%AIRMW

   ! Met. data required by module
   ! TODO: Add J-Values as field
   ExtState%O3%DoUse         = .TRUE.
   ExtState%NO2%DoUse        = .TRUE.
   ExtState%NO%DoUse         = .TRUE.
   ExtState%AIR%DoUse        = .TRUE.
   ExtState%SUNCOSmid%DoUse  = .TRUE.
   ExtState%SUNCOSmid5%DoUse = .TRUE.
   ExtState%TSURFK%DoUse     = .TRUE.
   ExtState%AIRVOL%DoUse     = .TRUE.
   IF ( IDTHNO3 > 0 ) ExtState%HNO3%DoUse = .TRUE.
   !ExtState%JVAL%DoUse = .TRUE.

   ! Enable module
   ExtState%ParaNOx = .TRUE.

   ! Leave w/ success
   IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  )
   IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
   CALL HCO_LEAVE ( RC )

 END SUBROUTINE HCOX_ParaNOx_Init
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!  06 Jun 2014 - R. Yantosca - Now indended with F90 free-format!
! !NOTES: 
!EOP
!------------------------------------------------------------------------------
!BOC

   !=================================================================
   ! HCOX_PARANOX_FINAL begins here!
   !=================================================================

   IF ( ALLOCATED(ShipNO) ) DEALLOCATE ( ShipNO )

 END SUBROUTINE HCOX_ParaNOx_Final
!EOC
END MODULE HCOX_ParaNOx_mod
