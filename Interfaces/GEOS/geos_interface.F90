#include "MAPL_Generic.h"
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: geos_interface
!
! !DESCRIPTION: Module with routines and variables to interface GEOS-Chem with
!  GEOS 
!\\
!\\
! !INTERFACE:
!
MODULE GEOS_Interface
!
! !USES:
!
  ! MAPL/ESMF
  USE ESMF     
  USE MAPL_Mod 
  ! GEOS-Chem
  USE Precision_Mod
  USE ErrCode_Mod                                    ! Error numbers
  USE PHYSCONSTANTS
  USE Input_Opt_Mod,         ONLY : OptInput
  USE State_Chm_Mod,         ONLY : ChmState, Ind_   ! Chemistry State obj
  USE State_Met_Mod,         ONLY : MetState         ! Meteorology State obj
  USE State_Diag_Mod,        ONLY : DgnState         ! Diagnostics State obj
  USE State_Grid_Mod,        ONLY : GrdState         ! Grid State obj

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC   :: MetVars_For_Lightning_Init 
  PUBLIC   :: MetVars_For_Lightning_Run  
  PUBLIC   :: GEOS_Diagnostics 
  PUBLIC   :: GEOS_CalcTotOzone
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE  :: CalcColumns_
  PRIVATE  :: CalcSpeciesDiagnostics_
!
! !PRIVATE TYPES:
!
! !REVISION HISTORY:
!  01 Jul 2022 - C. Keller - initial version (refactored Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for full history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_Diagnostics 
!
! !DESCRIPTION: Wrapper routine to handle all GEOS-specific diagnostics.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS_Diagnostics( GC, IMPORT,  EXPORT, Clock, Phase, &
                                Input_Opt,  State_Met, State_Chm, &
                                State_Diag, State_Grid, RC )
!
! !USES:
!
    USE Diagnostics_Mod,    ONLY : Set_Diagnostics_EndofTimestep
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT), TARGET :: GC     ! Ref to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT)         :: Import   ! Import State
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export   ! Export State
    TYPE(ESMF_Clock),    INTENT(INOUT)         :: Clock  ! ESMF Clock object
    INTEGER,             INTENT(IN   )         :: Phase  ! Run phase (-1/1/2)
    TYPE(OptInput),      INTENT(INOUT)         :: Input_Opt
    TYPE(MetState),      INTENT(INOUT)         :: State_Met
    TYPE(ChmState),      INTENT(INOUT)         :: State_Chm
    TYPE(DgnState),      INTENT(INOUT)         :: State_Diag
    TYPE(GrdState),      INTENT(INOUT)         :: State_Grid
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(INOUT)         :: RC       ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  08 Oct 2020 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    TYPE(MAPL_MetaComp), POINTER :: STATE => NULL()
    TYPE(ESMF_Alarm)             :: ALARM
    TYPE(ESMF_STATE)             :: IntState

    LOGICAL                      :: am_I_Root
    LOGICAL                      :: IsChemTime    ! Chemistry alarm proxy
    LOGICAL                      :: IsTendTime    ! Time to calculate tendencies

    INTEGER                      :: indO3, indSpc
    INTEGER                      :: I,  J,  L, N, LB
    INTEGER                      :: IM, JM, LM

    REAL, POINTER                :: Q(:,:,:)     => NULL()
    REAL, POINTER                :: PLE(:,:,:)   => NULL()
    REAL, POINTER                :: TROPP(:,:)   => NULL()

    REAL, POINTER                :: Ptr2D(:,:)     => NULL()
    REAL, POINTER                :: Ptr3D(:,:,:)   => NULL()
    REAL(fp), POINTER            :: PTR_O3(:,:,:)  => NULL()
    REAL, POINTER                :: OX(:,:,:)      => NULL()
    REAL, POINTER                :: O3(:,:,:)      => NULL()
    REAL, POINTER                :: O3PPMV(:,:,:)  => NULL()
    REAL, POINTER                :: PTR_O3P(:,:,:) => NULL()
    REAL, POINTER                :: PTR_O1D(:,:,:) => NULL()
    REAL, PARAMETER              :: OMW = 16.0

    REAL(f4), POINTER            :: O3_MASS(:,:,:) => NULL()

    ! LFR diag
    REAL                         :: lp1, lp2      ! lightning potentials
    REAL, POINTER                :: PtrEmis(:,:)  => NULL()
    REAL, POINTER                :: LWI(:,:)      => NULL()
    REAL, POINTER                :: LFR(:,:)      => NULL()
    REAL, POINTER                :: CNV_FRC(:,:)  => NULL()

    CHARACTER(LEN=ESMF_MAXSTR)   :: OrigUnit

    INTEGER, PARAMETER           :: NRATS = 5
    CHARACTER(LEN=15), PARAMETER :: RatsNames(NRATS) = (/ 'CH4', 'N2O', 'CFC11', 'CFC12', 'HCFC22' /)

    __Iam__('GEOS_Diagnostics')

    !=======================================================================
    ! GEOS_Diagnostics starts here
    !=======================================================================

    ! Are we on the root PET?
    am_I_Root = MAPL_Am_I_Root()

    ! Get MAPL Generic State
    CALL MAPL_GetObjectFromGC(GC, STATE, __RC__)

    ! Start timers
    CALL MAPL_TimerOn(STATE, "GC_DIAGN")

    ! Get Internal state
    CALL MAPL_Get ( STATE, INTERNAL_ESMF_STATE=IntState, __RC__ )

    ! Timers
    CALL MAPL_Get(STATE, RUNALARM=ALARM, __RC__)
    IsChemTime = ESMF_AlarmIsRinging(ALARM, __RC__)
    IsTendTime = ( IsChemTime .AND. Phase /= 1 )

    CALL MAPL_GetPointer( IMPORT,     Q,     'Q', __RC__ )
    CALL MAPL_GetPointer( IMPORT,   PLE,   'PLE', __RC__ )
    CALL MAPL_GetPointer( IMPORT, TROPP, 'TROPP', __RC__ )

    ! Grid size
    IM = SIZE(Q,1); JM = SIZE(Q,2); LM = SIZE(Q,3)

    ! Move 'regular' GEOS-Chem diagnostics from gchp_chunk_mod.F90 to here to
    ! make sure that these diagnostics see any post-run updates.
    ! Diagnostics routine expects units of kg/kg dry. 
    CALL Convert_Spc_Units ( Input_Opt, State_Chm, State_Grid, State_Met, &
                             'kg/kg dry', RC, OrigUnit=OrigUnit )
    _ASSERT(RC==GC_SUCCESS, 'Error calling CONVERT_SPC_UNITS')
    CALL Set_Diagnostics_EndofTimestep( Input_Opt,  State_Chm, State_Diag, &
                                        State_Grid, State_Met, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling Set_Diagnostics_EndofTimestep')
    CALL Convert_Spc_Units ( Input_Opt, State_Chm, State_Grid, State_Met, &
                             OrigUnit, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling CONVERT_SPC_UNITS')

    !=======================================================================
    ! Dry volume mixing ratios and PM2.5 diagnostics
    !=======================================================================
    CALL CalcSpeciesDiagnostics_ ( am_I_Root, Input_Opt, State_Met, State_Chm, &
                                   State_Diag, IMPORT, EXPORT, Q, __RC__ )

    !=======================================================================
    ! Ozone diagnostics for GEOS coupling with other components. Do these
    ! via the export state directly, rather than using the State_Diag obj.
    !=======================================================================

    ! PTR_O3: kg kg-1 total air
    !CALL MAPL_GetPointer( INTSTATE, PTR_O3, 'SPC_O3', NotFoundOk=.TRUE., __RC__ )

    ! Fill ozone export states if GC is the analysis OX provider:
    !      OX: volume mixing ratio
    !      O3: mass mixing ratio
    !  O3PPMV: volume mixing ratio in ppm
    ! Get pointers to analysis OX exports
    CALL MAPL_GetPointer ( INTSTATE,    OX, 'OX'         , NotFoundOK=.TRUE., __RC__ )
    CALL MAPL_GetPointer ( EXPORT,      O3, 'GCC_O3'     , NotFoundOK=.TRUE., __RC__ )
    CALL MAPL_GetPointer ( EXPORT,  O3PPMV, 'GCC_O3PPMV' , NotFoundOK=.TRUE., __RC__ )

    IF ( ASSOCIATED(O3) .OR. ASSOCIATED(O3PPMV) .OR. ASSOCIATED(OX) ) THEN
       indO3 = -1
       DO I=1,State_Chm%nSpecies
          IF ( TRIM(State_Chm%SpcData(I)%Info%Name) == 'O3' ) THEN
             indO3 = I
             EXIT
          ENDIF
       ENDDO
       ASSERT_(indO3>0)
       PTR_O3 => State_Chm%Species(indO3)%Conc(:,:,LM:1:-1)
    ENDIF
    IF ( ASSOCIATED(O3)     ) O3     = PTR_O3
    IF ( ASSOCIATED(O3PPMV) ) O3PPMV = PTR_O3 * MAPL_AIRMW / MAPL_O3MW * 1.00E+06
    IF ( ASSOCIATED(OX) ) THEN
       OX = PTR_O3  * MAPL_AIRMW / MAPL_O3MW
       IndSpc = Ind_('O')
       ASSERT_(IndSpc>0)
       OX = OX + ( State_Chm%Species(indSpc)%Conc(:,:,LM:1:-1)*MAPL_AIRMW/State_Chm%SpcData(IndSpc)%Info%MW_g )
       IndSpc = Ind_('O1D')
       ASSERT_(IndSpc>0)
       OX = OX + ( State_Chm%Species(indSpc)%Conc(:,:,LM:1:-1)*MAPL_AIRMW/State_Chm%SpcData(IndSpc)%Info%MW_g )
    ENDIF

    !=======================================================================
    ! Ozone diagnostics handled through State_Diag object
    !=======================================================================

    ! Total ozone and total tropospheric ozone for export [dobsons]. 2.69E+20 per dobson.
    CALL GEOS_CalcTotOzone( am_I_Root, State_Met, State_Chm, State_Diag, PLE, TROPP, __RC__ )

    ! O3 mass in kg/m2
    IF ( State_Diag%Archive_O3_MASS .AND. ASSOCIATED(State_Diag%O3_MASS) ) THEN
       O3_MASS => State_Diag%O3_MASS(:,:,LM:1:-1)
       LB = LBOUND(PLE,3)
       DO L=1,LM
          O3_MASS(:,:,L)=PTR_O3(:,:,L)*(g0_100*(PLE(:,:,L+LB)-PLE(:,:,L+LB-1)))
       ENDDO
    ENDIF
    O3_MASS => NULL()

    !=======================================================================
    ! Fill RATS export states if GC is the RATS provider
    ! The tracer concentrations of the RATS export states are in mol mol-1.
    ! These fields are required for coupling with other components. Don't
    ! do this via the State_Diag object but use the EXPORT state directly.
    !=======================================================================
    ! Get pointers to RATS exports
    DO I=1,NRATS
       CALL MAPL_GetPointer ( EXPORT, Ptr3D, TRIM(RatsNames(I)), NotFoundOK=.TRUE., __RC__ )
       IF ( ASSOCIATED(Ptr3D) ) THEN
          IndSpc = Ind_(TRIM(RatsNames(I)))
          ASSERT_(IndSpc>0)
          Ptr3D = State_Chm%Species(IndSpc)%Conc(:,:,LM:1:-1) &
                * ( MAPL_AIRMW / State_Chm%SpcData(IndSpc)%Info%MW_g )
       ENDIF
    ENDDO

    !=======================================================================
    ! Total and tropospheric columns
    !=======================================================================
    CALL CalcColumns_( am_I_Root, Input_Opt, State_Chm, State_Diag, PLE, TROPP, __RC__ )

    !=======================================================================
    ! Derived met. diagnostics relevant to chemistry processes
    !=======================================================================
    IF ( Phase /= 1 ) THEN
       ! chemistry top level
       IF ( State_Diag%Archive_CHEMTOP .AND. &
            ASSOCIATED(State_Diag%CHEMTOP) ) THEN
          DO J = 1, JM
          DO I = 1, IM
             State_Diag%CHEMTOP(I,J) = LM - State_Met%ChemGridLev(I,J) + 1
          ENDDO
          ENDDO
       ENDIF

       ! chemistry tropopause
       IF ( State_Diag%Archive_CHEMTROPP .AND. &
            ASSOCIATED(State_Diag%CHEMTOP) ) THEN
          State_Diag%CHEMTOP(:,:) = State_Met%TROPP(:,:) * 100.0 ! hPa -> Pa
       ENDIF
    ENDIF

    ! convective cloud top height
    IF ( Phase /= 2 ) THEN
       IF ( State_Diag%Archive_CONVCLDTOP .AND. &
            ASSOCIATED(State_Diag%CONVCLDTOP) ) THEN
          State_Diag%CONVCLDTOP(:,:) = 0.0
          DO J = 1, JM
          DO I = 1, IM
             DO L = 1, LM
                IF ( State_Met%CMFMC(I,J,L) > 0.0d0 ) THEN
                   State_Diag%CONVCLDTOP(I,J) = REAL(LM-L+1,f4)
                   EXIT
                ENDIF
             ENDDO
          ENDDO
          ENDDO
       ENDIF
    ENDIF

    !=======================================================================
    ! Lightning potential (from GEOS lightning flash rates and convective
    ! fraction)
    !=======================================================================
    IF ( Phase /= 2 ) THEN
       ! convective cloud top height
       CALL MAPL_GetPointer( EXPORT, Ptr2D, 'LightningPotential', &
                             NotFoundOk=.TRUE., __RC__ )
       IF ( State_Diag%Archive_LGHTPOTENTIAL .AND. &
            ASSOCIATED(State_Diag%LightningPotential) ) THEN
          CALL MAPL_GetPointer( IMPORT, LFR,     'LFR_GCC', __RC__ )
          CALL MAPL_GetPointer( IMPORT, LWI,     'LWI', __RC__ )
          CALL MAPL_GetPointer( IMPORT, CNV_FRC, 'CNV_FRC', __RC__ )
          CALL MAPL_GetPointer( EXPORT, PtrEmis, 'EMIS_NO_LGHT', NotFoundOk=.TRUE., __RC__ )
          State_Diag%LightningPotential(:,:) = 0.0
          DO J = 1, JM
          DO I = 1, IM
             lp1 = 0.0
             lp2 = 0.0

             ! If there are HEMCO lightning emissions in current grid box set
             ! lightning potential accordingly
             IF ( ASSOCIATED(PtrEmis) ) THEN
                IF ( LWI(I,J) == 1 ) THEN
                   lp1 = PtrEmis(I,J) / 1.0e-11 ! Land
                ELSE
                   lp1 = PtrEmis(I,J) / 1.0e-13 ! Water/Ice
                ENDIF
                lp1 = MIN(MAX(0.25,lp1),1.00)
             ENDIF

             ! Lightning flash rate
             IF ( LFR(I,J) > 0.0 ) THEN
                IF ( LWI(I,J) == 1 ) THEN
                   lp2 = LFR(I,J) / 5.0e-07 ! Land
                ELSE
                   lp2 = LFR(I,J) / 1.0e-08 ! Water/Ice
                ENDIF
                lp2 = MIN(MAX(0.25,lp2),1.00)

             ! Convective fraction
             ELSE
                lp2 = CNV_FRC(I,J)
             ENDIF

             ! Take highest value
             State_Diag%LightningPotential(I,J) = MAX(lp1,lp2)
          ENDDO
          ENDDO
       ENDIF
       PtrEmis => NULL()
    ENDIF

    ! Start timers
    CALL MAPL_TimerOff(STATE, "GC_DIAGN")

    _RETURN(ESMF_SUCCESS)

    END SUBROUTINE GEOS_Diagnostics
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_CalcTotOzone
!
! !DESCRIPTION: GEOS_CalcTotOzone calculates total ozone for the entire
!  atmosphere and troposphere only (in dobsons) and writes them into
!  the export variables GCCTO3 and GCCTTO3, respectively. Expects O3 in the
!  internal state in kg/kg total.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS_CalcTotOzone ( am_I_Root, State_Met, State_Chm, State_Diag, PLE, TROPP, RC )
!
! !USES:
!
    USE Precision_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)     :: am_I_Root
    TYPE(MetState),   INTENT(INOUT)  :: State_Met
    TYPE(ChmState),   INTENT(INOUT)  :: State_Chm
    TYPE(DgnState),   INTENT(INOUT)  :: State_Diag
    REAL,             POINTER        :: PLE  (:,:,:)
    REAL,             POINTER        :: TROPP(:,:)
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT), OPTIONAL   :: RC
!
! !REVISION HISTORY:
!  25 Oct 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
! 
    REAL(fp), POINTER            :: O3 (:,:,:) => NULL()
    REAL(fp), POINTER            :: TO3fp(:,:) => NULL()
    REAL(f4), POINTER            :: TO3 (:,:)  => NULL()
    REAL(f4), POINTER            :: TTO3(:,:)  => NULL()

    REAL,  ALLOCATABLE           :: DUsLayerL(:,:)! Dobsons in a layer,
                                                  !  for total ozone
    REAL,  ALLOCATABLE           :: wgt(:,:)      ! Layer thickness weighting
                                                  !  for total ozone
    REAL                         :: const
    INTEGER                      :: indO3
    INTEGER                      :: IM, JM, LM, LB, L, STATUS
    CHARACTER(LEN=ESMF_MAXSTR)   :: Iam

    !=======================================================================
    ! GEOS_CalcTotOzone begins here
    !=======================================================================

    ! Traceback handle
    Iam = 'GEOS_CalcTotOzone'

    ! Check if we need to compute this
    IF ( ASSOCIATED( State_Met%TO3 )    ) TO3fp => State_Met%TO3
    IF ( State_Diag%Archive_GCCTO3 .AND. &
         ASSOCIATED(State_Diag%GCCTO3)  ) TO3   => State_Diag%GCCTO3
    IF ( State_Diag%Archive_GCCTTO3 .AND. &
         ASSOCIATED(State_Diag%GCCTTO3) ) TTO3  => State_Diag%GCCTTO3

    ! Nothing to do if neither of the arrays is associated
    IF ( .NOT. ASSOCIATED(TO3) .AND. .NOT. ASSOCIATED(TTO3) .AND. .NOT. ASSOCIATED(TO3fp) ) THEN
       RC = ESMF_SUCCESS
       RETURN
    ENDIF

    ! Get O3 from species array (kg/kg total)
    indO3 = Ind_('O3')
    O3 => State_Chm%Species(indO3)%Conc(:,:,:)

    ! Grid size
    IM = SIZE(O3,1)
    JM = SIZE(O3,2)
    LM = SIZE(O3,3)

    ! Pressure edges
    LB = LBOUND(PLE,3)

    ! Reset values
    IF ( ASSOCIATED(TO3fp ) ) TO3fp  = 0.0
    IF ( ASSOCIATED(TO3   ) ) TO3  = 0.0
    IF ( ASSOCIATED(TTO3  ) ) TTO3 = 0.0

    ! Allocate local variables
    ALLOCATE(DUsLayerL(IM,JM), STAT=STATUS)
    _VERIFY(STATUS)
    ALLOCATE(wgt(IM,JM), STAT=STATUS)
    _VERIFY(STATUS)

    ! constant
    const = 0.01 * MAPL_AVOGAD / ( MAPL_GRAV * (MAPL_AIRMW/1000.0) )
    const = const * MAPL_AIRMW / MAPL_O3MW ! convert kg/kg total to v/v total

    ! Calculate total ozone
    DO L = 1,LM
       DUsLayerL(:,:) = O3(:,:,LM-L+1) * ((PLE(:,:,L+LB)-PLE(:,:,L+LB-1))/100.0) &
                        * const / 2.69e16 / 1000.0
       IF ( ASSOCIATED(TO3fp) ) TO3fp = TO3fp+DUsLayerL
       IF ( ASSOCIATED(TO3  ) ) TO3   = TO3  +DUsLayerL
       IF ( ASSOCIATED(TTO3) ) THEN
          wgt  = MAX(0.0,MIN(1.0,(PLE(:,:,L+LB)-TROPP(:,:)) &
                 /(PLE(:,:,L+LB)-PLE(:,:,L+LB-1))))
          TTO3 = TTO3+DUsLayerL*wgt
       END IF
    END DO

    ! Cleanup
    IF ( ASSOCIATED(TO3fp) ) TO3fp => NULL()
    IF ( ASSOCIATED(TO3  ) ) TO3   => NULL()
    IF ( ASSOCIATED(TTO3 ) ) TTO3  => NULL()
    DEALLOCATE(DUsLayerL, STAT=STATUS)
    _VERIFY(STATUS)
    DEALLOCATE(wgt, STAT=STATUS)
    _VERIFY(STATUS)

    ! Successful return
    RC = ESMF_SUCCESS

  END SUBROUTINE GEOS_CalcTotOzone
!EOC



!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetVars_For_Lightning_Init
!
! !DESCRIPTION: Initialize the imports to fill the met variables needed for
!               lightning NOx computation
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MetVars_For_Lightning_Init( GC, CF, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC        ! Ref. to this GridComp
    TYPE(ESMF_Config),   INTENT(INOUT)         :: CF        ! GEOSCHEM*.rc
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)           :: RC        ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  20 Jan 2020 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    CHARACTER(LEN=31)  :: LfrSrc, CnvSrc

    __Iam__('MetVars_For_Lightning_Init')

    ! Get source for lightning fields
    CALL MetVars_For_Lightning_Run( GC, DryRun=.TRUE., CF=CF, &
                                    LfrSrc=LfrSrc, CnvSrc=CnvSrc, __RC__ )

    ! LFR import - always pass LFR and LFR_GCC. Depending on specification,
    ! also provide import from external file
    call MAPL_AddImportSpec(GC,                     &
               SHORT_NAME='LFR',                    &
               LONG_NAME ='lightning_flash_rate',   &
               UNITS     ='km-2 s-1',               &
               DIMS      = MAPL_DimsHorzOnly,       &
               VLOCATION = MAPL_VLocationNone,      &
                                              __RC__ )

    call MAPL_AddImportSpec(GC,                     &
               SHORT_NAME='LFR_GCC',                &
               LONG_NAME ='lightning_flash_rate',   &
               UNITS     ='km-2 s-1',               &
               DIMS      = MAPL_DimsHorzOnly,       &
               VLOCATION = MAPL_VLocationNone,      &
                                              __RC__ )

    IF ( (TRIM(LfrSrc)/='LFR') .AND. &
         (TRIM(LfrSrc)/='LFR_GCC')    ) THEN
       call MAPL_AddImportSpec(GC,                     &
                  SHORT_NAME=TRIM(LfrSrc),             &
                  LONG_NAME ='lightning_flash_rate',   &
                  UNITS     ='km-2 s-1',               &
                  DIMS      = MAPL_DimsHorzOnly,       &
                  VLOCATION = MAPL_VLocationNone,      &
                                                 __RC__ )
    ENDIF

    ! Import fields needed to compute convective height, depending on specification
    SELECT CASE ( TRIM(CnvSrc) )
       CASE ( 'CNV_MFC' )
          ! CNV_MFC is always imported, nothing to do here
          !CONTINUE

       CASE ( 'BYNCY' )
          call MAPL_AddImportSpec(GC,                   &
             SHORT_NAME = 'BYNCY',                      &
             LONG_NAME  ='buoyancy_of surface_parcel',  &
             UNITS      ='m s-2',                       &
             DIMS       = MAPL_DimsHorzVert,            &
             VLOCATION  = MAPL_VLocationCenter,         &
                                                  __RC__ )

       CASE DEFAULT
          call MAPL_AddImportSpec(GC,                       &
               SHORT_NAME=TRIM(CnvSrc),                     &
               LONG_NAME ='convective_cloud_top_from_file', &
               UNITS     ='m',                              &
               DIMS      = MAPL_DimsHorzOnly,               &
               VLOCATION = MAPL_VLocationNone,              &
                                                      __RC__ )

    END SELECT

    ! Also add export for CONV_DEPTH_GCC & LFR diagnostics
    call MAPL_AddExportSpec(GC,                                    &
               SHORT_NAME='GCC_CONV_DEPTH',                        &
               LONG_NAME ='Convective_depth_seen_by_GEOSCHEMchem', &
               UNITS     ='m',                                     &
               DIMS      = MAPL_DimsHorzOnly,                      &
               VLOCATION = MAPL_VLocationNone,                     &
                                                             __RC__ )

    call MAPL_AddExportSpec(GC,                                     &
               SHORT_NAME='GCC_LFR',                                &
               LONG_NAME ='Lightning_flash_rate_seen_GEOSCHEMchem', &
               UNITS     ='km-2 s-1',                               &
               DIMS      = MAPL_DimsHorzOnly,                       &
               VLOCATION = MAPL_VLocationNone,                      &
                                                              __RC__ )

    ! All done
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE MetVars_For_Lightning_Init
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MetVars_For_Lightning_Run
!
! !DESCRIPTION: Fill the State_Met variables needed for lightning NOx calculation
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MetVars_For_Lightning_Run( GC, Import, Export, State_Met, State_Grid, &
                                        DryRun, CF, LfrSrc, CnvSrc, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT)           :: GC         ! Ref. to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT), OPTIONAL :: Import     ! Import State
    TYPE(ESMF_State),    INTENT(INOUT), OPTIONAL :: Export     ! Export State
    TYPE(MetState),      INTENT(INOUT), OPTIONAL :: State_Met  ! Met. state object
    TYPE(GrdState),      INTENT(IN),    OPTIONAL :: State_Grid ! Grid state
    LOGICAL,             INTENT(IN),    OPTIONAL :: DryRun     ! Don't fill fields
    TYPE(ESMF_Config),   INTENT(INOUT), OPTIONAL :: CF         ! GEOSCHEM*.rc
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*),    INTENT(OUT), OPTIONAL   :: LfrSrc     ! Lightning flash rate source ID
    CHARACTER(LEN=*),    INTENT(OUT), OPTIONAL   :: CnvSrc     ! Convective height source ID
    INTEGER,             INTENT(OUT)             :: RC         ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  20 Jan 2020 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
  INTEGER                       :: I,J,L
  INTEGER                       :: LTOP
  LOGICAL                       :: am_I_Root
  LOGICAL                       :: Skip
  CHARACTER(LEN=31), SAVE       :: LFR_SOURCE = ""
  CHARACTER(LEN=31), SAVE       :: CNV_SOURCE = ""
  INTEGER, SAVE                 :: CNV_ID = -1
  REAL, SAVE                    :: SCAL_STRP = 1.0
  REAL, SAVE                    :: SCAL_TROP = 1.0
  REAL, SAVE                    :: SCAL_NTRP = 1.0
  REAL, SAVE                    :: SCAL_LFR  = 1.0
  REAL, POINTER                 :: Ptr2d(:,:)
  REAL, POINTER                 :: BYNCY(:,:,:)
  REAL, POINTER                 :: CNV_FRC(:,:)

  __Iam__('MetVars_For_Lightning_Run')
  am_I_Root = MAPL_Am_I_Root()

!-LFR source
  IF ( TRIM(LFR_SOURCE)=="" .OR. CNV_ID<0 ) THEN
     ASSERT_(PRESENT(CF))
     CALL ESMF_ConfigGetAttribute( CF, LFR_SOURCE,                     &
                                 Label="LIGHTNING_FLASH_RATE_SOURCE:", &
                                 Default="LFR_GCC",                    &
                                 __RC__                                 )
     CALL ESMF_ConfigGetAttribute( CF, SCAL_STRP,                      &
                                 Label="LFR_SCALING_SOUTHERN_TROP:",   &
                                 Default=1.0,                          &
                                 __RC__                                 )
     CALL ESMF_ConfigGetAttribute( CF, SCAL_TROP,                      &
                                 Label="LFR_SCALING_TROPICS:",         &
                                 Default=1.0,                          &
                                 __RC__                                 )
     CALL ESMF_ConfigGetAttribute( CF, SCAL_NTRP,                      &
                                 Label="LFR_SCALING_NORTHERN_TROP:",   &
                                 Default=1.0,                          &
                                 __RC__                                 )
     CALL ESMF_ConfigGetAttribute( CF, SCAL_LFR,                       &
                                 Label="LFR_SCALING_GLOBAL:",          &
                                 Default=1.0,                          &
                                 __RC__                                 )
     ! Verbose
     IF (am_I_Root) THEN
        WRITE(*,*) 'GEOSCHEMchem lightning flash rate source: ',TRIM(LFR_SOURCE)
        WRITE(*,*) '--> LFR scaling southern trop (<23S)    : ',SCAL_STRP
        WRITE(*,*) '--> LFR scaling tropics (23S-23N)       : ',SCAL_TROP
        WRITE(*,*) '--> LFR scaling northern trop (>23N)    : ',SCAL_NTRP
        WRITE(*,*) '--> LFR scaling global                  : ',SCAL_LFR
     ENDIF

!----Convective height source
     CALL ESMF_ConfigGetAttribute( CF, CNV_SOURCE,                         &
                                 Label="LIGHTNING_CONVECTIVE_TOP_SOURCE:", &
                                 Default="CNV_MFC",                        &
                                 __RC__                                     )
     SELECT CASE ( TRIM(CNV_SOURCE) )
        CASE ( "CNV_MFC" )
           CNV_ID = 0
        CASE ( "BYNCY" )
           CNV_ID = 1
        CASE DEFAULT
           CNV_ID = 2
     END SELECT

     ! Verbose
     IF (am_I_Root) THEN
        WRITE(*,*) 'GEOSCHEMchem lightning convective height source: ',TRIM(CNV_SOURCE)
     ENDIF

  ENDIF

!-Fill state met
  IF ( PRESENT(DryRun) ) THEN
     Skip = DryRun
  ELSE
     Skip = .FALSE.
  ENDIF
  IF ( .NOT. Skip ) THEN

!----Lightning flash rate density [km-2 s-1]
     call MAPL_GetPointer ( IMPORT, Ptr2D, TRIM(LFR_SOURCE), __RC__ )
     State_Met%FLASH_DENS = Ptr2D

     ! Rescale flash rates as specified in GEOSCHEMchem_GridComp.rc
     ! southern extratropics
     IF ( SCAL_STRP /= 1.0 ) THEN
         WHERE ( State_Grid%YMID < -23.0 )
             State_Met%FLASH_DENS = State_Met%FLASH_DENS * SCAL_STRP
         END WHERE
     ENDIF
     ! tropics
     IF ( SCAL_TROP /= 1.0 ) THEN
         WHERE ( State_Grid%YMID >= -23.0 .AND. State_Grid%YMID <= 23.0 )
             State_Met%FLASH_DENS = State_Met%FLASH_DENS * SCAL_TROP
         END WHERE
     ENDIF
     ! northern extratropics
     IF ( SCAL_NTRP /= 1.0 ) THEN
         WHERE ( State_Grid%YMID > 23.0 )
             State_Met%FLASH_DENS = State_Met%FLASH_DENS * SCAL_NTRP
         END WHERE
     ENDIF
     ! overall LFR scaling
     IF ( SCAL_LFR /= 1.0 ) THEN
        State_Met%FLASH_DENS = State_Met%FLASH_DENS * SCAL_LFR
     ENDIF

     ! Eventually add to Export
     Ptr2D => NULL()
     call MAPL_GetPointer ( EXPORT, Ptr2D, 'GCC_LFR', NotFoundOk=.TRUE., __RC__ )
     IF ( ASSOCIATED(Ptr2D) ) Ptr2D = State_Met%FLASH_DENS

!----Convective depth [m]
     SELECT CASE ( CNV_ID )
        ! Convective mass flux
        ! Get highest level with positive convective mass flux. CMFMC is  
        ! on level edges.
        CASE ( 0 )
           DO J=1,State_Grid%NY
           DO I=1,State_Grid%NX
              LTOP = 0
              DO L = State_Grid%NZ+1,2,-1
                 IF ( State_Met%CMFMC(I,J,L) > 0.0 ) THEN
                    LTOP = L-1
                    EXIT
                 ENDIF
              ENDDO
              IF ( LTOP > 0 ) THEN
                 State_Met%CONV_DEPTH(I,J) = SUM(State_Met%BXHEIGHT(I,J,1:LTOP))
              ELSE
                 State_Met%CONV_DEPTH(I,J) = 0.0
              ENDIF
           ENDDO
           ENDDO

        ! Buoyancy and convective fraction
        ! Get highest level with positive buoyancy and where convective fraction
        ! is non-zero. BYNCY is on GEOS coordinates (--> 1=top of atmosphere) and
        ! on level mid-points. LM captures the dimension of CNV_MFC, which is on
        ! level edges.
        CASE ( 1 )
           call MAPL_GetPointer ( IMPORT, BYNCY,   'BYNCY'  , __RC__ )
           call MAPL_GetPointer ( IMPORT, CNV_FRC, 'CNV_FRC', __RC__ )
           DO J=1,State_Grid%NY
           DO I=1,State_Grid%NX
              LTOP = 0
              IF ( CNV_FRC(I,J) > 0.0 ) THEN
                 DO L = 1,State_Grid%NZ
                    IF ( BYNCY(I,J,L) > 0.0 ) THEN
                       LTOP = State_Grid%NZ - L + 1
                       EXIT
                    ENDIF
                 ENDDO
              ENDIF
              IF ( LTOP > 0 ) THEN
                 State_Met%CONV_DEPTH(I,J) = SUM(State_Met%BXHEIGHT(I,J,1:LTOP))
              ELSE
                 State_Met%CONV_DEPTH(I,J) = 0.0
              ENDIF
           ENDDO
           ENDDO
           BYNCY   => NULL()
           CNV_FRC => NULL()

        ! Offline file
        CASE ( 2 )
           call MAPL_GetPointer ( IMPORT, Ptr2D, TRIM(CNV_SOURCE), __RC__ )
           State_Met%CONV_DEPTH = Ptr2D
     END SELECT

     ! Eventually add to Export
     Ptr2D => NULL()
     call MAPL_GetPointer ( EXPORT, Ptr2D, 'GCC_CONV_DEPTH', NotFoundOk=.TRUE., __RC__ )
     IF ( ASSOCIATED(Ptr2D) ) Ptr2D = State_Met%CONV_DEPTH

  ENDIF ! Skip

!-Cleanup
  IF ( PRESENT(LfrSrc) ) LfrSrc = LFR_SOURCE
  IF ( PRESENT(CnvSrc) ) CnvSrc = CNV_SOURCE
  RETURN_(ESMF_SUCCESS)

  END SUBROUTINE MetVars_For_Lightning_Run
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CalcColumns_
!
! !DESCRIPTION: CalcColumns_ calculates total and tropospheric columns for a
!  number of species.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CalcColumns_ ( am_I_Root, Input_Opt, State_Chm, State_Diag, PLE, TROPP, RC )
!
! !USES:
!
    USE State_Diag_Mod, ONLY : DgnMap
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)            :: am_I_Root
    TYPE(OptInput),   INTENT(INOUT)         :: Input_Opt
    TYPE(ChmState),   INTENT(INOUT)         :: State_Chm
    TYPE(DgnState),   INTENT(INOUT)         :: State_Diag
    REAL,             POINTER               :: PLE  (:,:,:)
    REAL,             POINTER               :: TROPP(:,:  )
    INTEGER,          INTENT(OUT)           :: RC
!
! !REVISION HISTORY:
!  25 Oct 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    REAL,  POINTER               :: ExpTOTCOL(:,:)
    REAL,  POINTER               :: ExpTRPCOL(:,:)
    REAL(fp), POINTER            :: Spc3D    (:,:,:)
    REAL,  ALLOCATABLE           :: DUsLayerL(:,:)! Dobsons in a layer, 
                                                  !  for total ozone
    REAL,  ALLOCATABLE           :: wgt(:,:)      ! Layer thickness weighting
                                                  !  for total ozone
    REAL                         :: MW, const
    INTEGER                      :: I, J, IM, JM, LM, LB, L, STATUS
    INTEGER                      :: ID, TotID, TropID
    CHARACTER(LEN=ESMF_MAXSTR)   :: Iam
    CHARACTER(LEN=15)            :: ISPEC

    ! Objects
    TYPE(DgnMap), POINTER :: mapTotCol  => NULL()
    TYPE(DgnMap), POINTER :: mapTropCol => NULL()

    !=======================================================================
    ! CalcColumns_ begins here
    !=======================================================================

    ! Traceback handle
    Iam = 'CalcColumns_'

    ! Nothing to do if not active
    IF ( .NOT. State_Diag%Archive_TotCol .AND. &
         .NOT. State_Diag%Archive_TropCol       ) THEN
       RC = ESMF_SUCCESS
       RETURN
    ENDIF

    ! Grid size
    IM = SIZE(PLE,1)
    JM = SIZE(PLE,2)
    LM = SIZE(PLE,3)-1
    LB = LBOUND(PLE,3)

    ! mapping objects
    IF ( State_Diag%Archive_TotCol  ) THEN
       mapTotCol  => State_Diag%Map_TotCol
       State_Diag%TotCol(:,:,:) = 0.0
    ENDIF
    IF ( State_Diag%Archive_TropCol ) THEN
       mapTropCol => State_Diag%Map_TropCol
       State_Diag%TropCol(:,:,:) = 0.0
    ENDIF

    ! Allocate local variables
    ALLOCATE(DUsLayerL(IM,JM), STAT=STATUS)
    _VERIFY(STATUS)
    ALLOCATE(wgt(IM,JM), STAT=STATUS)
    _VERIFY(STATUS)

    ! Check all species
    DO I = 1, State_Chm%nSpecies

       ! Check if total column and/or trop. column requested for this species
       TotID = -1
       DO J = 1,mapTotCol%nSlots
           IF ( mapTotCol%slot2id(J)==I ) THEN
              TotID = J
              EXIT
           ENDIF
       ENDDO
       TropID = -1
       DO J = 1,mapTropCol%nSlots
           IF ( mapTropCol%slot2id(J)==I ) THEN
              TropID = J
              EXIT
           ENDIF
       ENDDO
       IF ( (TotID<0) .AND. (TropID<0) ) CYCLE

       ! Species info
       ISPEC = State_Chm%SpcData(I)%Info%Name
       ID    = IND_(TRIM(ISPEC))
       MW    = State_Chm%SpcData(ID)%Info%MW_g

       ! Get species from internal state
       Spc3D => State_Chm%Species(ID)%Conc(:,:,LM:1:-1)

       ! constant 
       const = MAPL_AVOGAD / ( MAPL_GRAV * MW )

       ! Calculate total and trop. column
       DO L = 1,LM
          DUsLayerL(:,:) = Spc3D(:,:,L) * ( PLE(:,:,L+LB) &
                           - PLE(:,:,L+LB-1) ) * const
          ! rescale: molec/m2 --> molec/cm2
          ! rescale: molec/cm2 ==> 1.0e15 molec/cm2
          DUsLayerL(:,:) = DUsLayerL(:,:) / 1.0e4 / 1.0e15
          ! Add to total column
          IF ( TotID > 0 ) THEN
             State_Diag%TotCol(:,:,TotID) = State_Diag%TotCol(:,:,TotID) &
                                          + DUsLayerL(:,:)
          ENDIF
          ! Add to tropospheric column
          IF ( TropID > 0 ) THEN
             wgt = MAX(0.0,MIN(1.0,(PLE(:,:,L+LB)-TROPP(:,:)) &
                 / (PLE(:,:,L+LB)-PLE(:,:,L+LB-1))))
             State_Diag%TropCol(:,:,TropID) = State_Diag%TropCol(:,:,TropID) &
                                            + DUsLayerL(:,:)*wgt(:,:)
          END IF
       END DO
    ENDDO

    ! Cleanup
    DEALLOCATE(DUsLayerL, STAT=STATUS)
    _VERIFY(STATUS)
    DEALLOCATE(wgt, STAT=STATUS)
    _VERIFY(STATUS)

    ! Successful return
    RC = ESMF_SUCCESS

  END SUBROUTINE CalcColumns_
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CalcSpeciesDiagnostics_
!
! !DESCRIPTION: CalcSpeciesDiagnostics_ computes species' diagnostics
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CalcSpeciesDiagnostics_( am_I_Root, Input_Opt, State_Met, &
                                      State_Chm, State_Diag, IMPORT, EXPORT, &
                                      Q, RC )
!
! !USES:
!
!    USE TENDENCIES_MOD,          ONLY : Tend_Get
    USE Species_Mod,   ONLY : Species
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)            :: am_I_Root
    TYPE(OptInput),      INTENT(INOUT)         :: Input_Opt
    TYPE(MetState),      INTENT(INOUT)         :: State_Met
    TYPE(ChmState),      INTENT(INOUT)         :: State_Chm
    TYPE(DgnState),      INTENT(INOUT)         :: State_Diag
    TYPE(ESMF_State),    INTENT(INOUT)         :: Import   ! Import State
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export   ! Export State
    REAL,                POINTER               :: Q(:,:,:)
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(INOUT)         :: RC       ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  05 Dec 2017 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Objects

    ! Scalars
    INTEGER                    :: STATUS
    INTEGER                    :: I, J, N, IM, JM, LM, DryID
    INTEGER                    :: IndSpc
    LOGICAL                    :: IsBry, IsNOy,  IsCly, IsOrgCl
    LOGICAL                    :: RunMe
    CHARACTER(LEN=ESMF_MAXSTR) :: Iam           ! Gridded component name
    CHARACTER(LEN=ESMF_MAXSTR) :: FieldName, SpcName
    REAL                       :: MW
    REAL                       :: BrCoeff, ClCoeff, OrgClCoeff
    REAL(fp), POINTER          :: PtrTmp(:,:,:)
    TYPE(Species), POINTER     :: SpcInfo
    REAL(f4), POINTER          :: NOy(:,:,:) => NULL()
    REAL(f4), POINTER          :: Bry(:,:,:) => NULL()
    REAL(f4), POINTER          :: Cly(:,:,:) => NULL()
    REAL(f4), POINTER          :: OrgCl(:,:,:) => NULL()

    LOGICAL, SAVE              :: FIRST = .TRUE.

    !=======================================================================
    ! Routine starts here
    !=======================================================================

    ! Identify this routine to MAPL
    Iam = 'GCC::CalcSpeciesDiagnostics_'

    ! Grid size
    IM = SIZE(Q,1)
    JM = SIZE(Q,2)
    LM = SIZE(Q,3)

    !=======================================================================
    ! Exports in dry vol mixing ratio (v/v dry). Includes NOy. Convert from
    ! kg/kg total.
    !=======================================================================
    IF ( State_Diag%Archive_NOy .AND.       &
         ASSOCIATED(State_Diag%NOy)          ) NOy => State_Diag%NOy(:,:,LM:1:-1)
    IF ( State_Diag%Archive_Bry .AND.       &
         ASSOCIATED(State_Diag%Bry)          ) Bry => State_Diag%Bry(:,:,LM:1:-1)
    IF ( State_Diag%Archive_Cly .AND.       &
         ASSOCIATED(State_Diag%Cly)          ) Cly => State_Diag%Cly(:,:,LM:1:-1)
    IF ( State_Diag%Archive_OrganicCl .AND. &
         ASSOCIATED(State_Diag%OrganicCl)    ) OrgCl => State_Diag%OrganicCl(:,:,LM:1:-1)
    IF ( ASSOCIATED(NOy)   ) NOy(:,:,:)   = 0.0
    IF ( ASSOCIATED(Bry)   ) Bry(:,:,:)   = 0.0
    IF ( ASSOCIATED(Cly)   ) Cly(:,:,:)   = 0.0
    IF ( ASSOCIATED(OrgCl) ) OrgCl(:,:,:) = 0.0

    DO N=1,State_Chm%nSpecies
       SpcInfo   => State_Chm%SpcData(N)%Info ! Species database
       SpcName   =  TRIM(SpcInfo%Name)

       ! Need to fill at least one export?
       RunMe = .FALSE.

       ! Is this a NOy species?
       IF ( ASSOCIATED(NOy) ) THEN
          SELECT CASE ( TRIM(SpcName) )
             CASE ( 'BrNO3', 'ClNO3', 'DHDN', 'ETHLN', 'HNO2', &
                    'HNO3',  'HNO4',  'HONIT'  )
                IsNOy = .TRUE.
             CASE ( 'IONITA', 'IPMN', 'ISN1', 'ISNIOA', 'ISNIOG' )
                IsNOy = .TRUE.
             CASE ( 'ISOPNB', 'ISOPND', 'MACRN', 'MPN', 'MVKN', &
                    'N2O5',   'NIT',    'NO',    'NO2', 'NO3' )
                IsNOy = .TRUE.
             CASE ( 'NPMN', 'ONIT', 'PAN', 'PROPNN', 'R4N2' )
                IsNOy = .TRUE.
             CASE DEFAULT
                IsNOy = .FALSE.
          END SELECT
       ELSE
          IsNOy = .FALSE.
       ENDIF
       IF ( IsNOy ) RunMe = .TRUE.

       ! Is this a Bry species?
       BrCoeff = 0.0
       IF ( ASSOCIATED(Bry) ) THEN
          SELECT CASE ( TRIM(SpcName) )
             CASE ( 'Br', 'BrO', 'HOBr', 'HBr', 'BrNO2', 'BrNO3', 'BrCl', 'IBr' )
                BrCoeff = 1.0
                IsBry   = .TRUE.
             CASE ( 'Br2' )
                BrCoeff = 2.0
                IsBry   = .TRUE.
             CASE DEFAULT
                IsBry = .FALSE.
          END SELECT
       ELSE
          IsBry = .FALSE.
       ENDIF
       IF ( IsBry ) RunMe = .TRUE.

       ! Is this a Cly species?
       ClCoeff = 0.0
       IF ( ASSOCIATED(Cly) ) THEN
          SELECT CASE ( TRIM(SpcName) )
             CASE ( 'Cl', 'ClO', 'OClO', 'ClOO', 'HOCl', 'HCl', 'ClNO2', 'ClNO3', 'BrCl', 'ICl' )
                ClCoeff = 1.0
                IsCly   = .TRUE.
             CASE ( 'Cl2', 'Cl2O2' )
                ClCoeff = 2.0
                IsCly   = .TRUE.
             CASE DEFAULT
                IsCly = .FALSE.
          END SELECT
       ELSE
          IsCly = .FALSE.
       ENDIF
       IF ( IsCly ) RunMe = .TRUE.

       ! Is this an OrgCl species?
       OrgClCoeff = 0.0
       IF ( ASSOCIATED(Cly) ) THEN
          SELECT CASE ( TRIM(SpcName) )
             CASE ( 'H1211', 'CFC115', 'CH3Cl', 'HCFC142b', 'HCFC22', 'CH2ICl' )
                OrgClCoeff = 1.0
                IsOrgCl    = .TRUE.
             CASE ( 'CFC114', 'CFC12', 'HCFC141b', 'HCFC123', 'CH2Cl2' )
                OrgClCoeff = 2.0
                IsOrgCl    = .TRUE.
             CASE ( 'CFC11', 'CFC113', 'CH3CCl3', 'CHCl3' )
                OrgClCoeff = 3.0
                IsOrgCl    = .TRUE.
             CASE ( 'CCl4' )
                OrgClCoeff = 4.0
                IsOrgCl    = .TRUE.
             CASE DEFAULT
                IsOrgCl = .FALSE.
          END SELECT
       ELSE
          IsOrgCl = .FALSE.
       ENDIF
       IF ( IsOrgCl ) RunMe = .TRUE.

       ! Fill exports
       IF ( RunMe ) THEN
          !FieldName = 'SPC_'//TRIM(SpcName)
          MW = SpcInfo%MW_g
          IF ( MW < 0.0 ) THEN
             ! Get species and set MW to 1.0. This is ok because the internal
             ! state uses a MW of 1.0 for all species
             MW = 1.0
             ! Cannot add to NOy if MW is unknown because it would screw up
             ! unit conversion
             IF ( IsNOy ) THEN
                IsNOy = .FALSE.
                IF ( am_I_Root .AND. FIRST ) THEN
                   write(*,*) 'WARNING: Ignore species for NOy computation' //&
                              '  because MW is unknown: ', TRIM(SpcName)
                ENDIF
             ENDIF
          ENDIF
          PtrTmp => State_Chm%Species(N)%Conc(:,:,LM:1:-1)
          IF ( STATUS /= ESMF_SUCCESS ) THEN
             WRITE(*,*) 'Error reading ',TRIM(SpcName)
             VERIFY_(STATUS)
          ENDIF

          ! NOy concentration
          IF ( IsNOy ) NOy = NOy + PtrTmp * ( MAPL_AIRMW / MW ) / ( 1.0 - Q )

          ! Bry concentration
          IF ( IsBry ) Bry = Bry + BrCoeff * PtrTmp * ( MAPL_AIRMW / MW ) / ( 1.0 - Q )

          ! Cly concentration
          IF ( IsCly ) Cly = Cly + ClCoeff * PtrTmp * ( MAPL_AIRMW / MW ) / ( 1.0 - Q )

          ! OrgCl concentration
          IF ( IsOrgCl ) OrgCl = OrgCl + OrgClCoeff * PtrTmp * ( MAPL_AIRMW / MW ) / ( 1.0 - Q )
       ENDIF
    ENDDO

    !=======================================================================
    ! All done
    !=======================================================================

    ! Cleanup
    IF ( ASSOCIATED(NOy)   ) NOy   => NULL()
    IF ( ASSOCIATED(Bry)   ) Bry   => NULL()
    IF ( ASSOCIATED(Cly)   ) Cly   => NULL()
    IF ( ASSOCIATED(OrgCl) ) OrgCl => NULL()

    ! Successful return
    FIRST = .FALSE.
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE CalcSpeciesDiagnostics_
!EOC
END MODULE GEOS_Interface
