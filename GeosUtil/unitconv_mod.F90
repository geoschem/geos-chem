!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: unitconv_mod.F90
!
! !DESCRIPTION: Module UNITCONV\_MOD contains routines which are used to 
!  convert the units of tracer concentrations between mass mixing ratio 
!  [kg/kg air], mass per grid box [kg], molar mixing ratio [vol/vol], and 
!  molecular number density [molecules/cm3]. There are different conversion 
!  routines for dry air and total (wet) air mixing ratios.
!\\  
!\\
! !INTERFACE: 
!
MODULE UnitConv_Mod
!
! !USES:
!
  ! GEOS-Chem Modules
  USE CMN_SIZE_MOD          ! Size parameters
  USE PRECISION_MOD         ! GEOS-Chem Flexible Precision (fp)
  USE GIGC_ErrCode_Mod
                    
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  ! kg/kg dry air <-> kg/kg moist air
  PUBLIC  :: Convert_DryKgKg_to_MoistKgKg
  PUBLIC  :: Convert_MoistKgKg_to_DryKgKg

  ! v/v dry air <-> v/v moist air
  PUBLIC  :: Convert_DryVV_to_MoistVV
  PUBLIC  :: Convert_MoistVV_to_DryVV

  ! kg/kg dry air <-> v/v dry air
  PUBLIC  :: Convert_DryKgKg_to_DryVV
  PUBLIC  :: Convert_DryVV_to_DryKgKg

  ! kg/kg moist air <-> kg/grid box
  PUBLIC  :: Convert_MoistKgKg_to_KG
  PUBLIC  :: Convert_KG_to_MoistKgKg

  ! kg/kg dry air <-> kg/grid box
  PUBLIC  :: Convert_DryKgKg_to_KG
  PUBLIC  :: Convert_KG_to_DryKgKg

  ! kg/kg dry air <-> molecules/cm3
  PUBLIC  :: Convert_DryKgKg_to_MND
  PUBLIC  :: Convert_MND_to_DryKgKg
!
! !REMARKS:
!  The routines in this module are used to convert the units of tracer 
!  concentrations in various GEOS-Chem routines.
!
! !REVISION HISTORY:
!  23 Jun 2015 - E. Lundgren - Initial version
!  10 Jul 2015 - R. Yantosca - Added cosmetic changes; fixed ProTeX issues
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_drykgkg_to_moistkgkg
!
! !DESCRIPTION: Subroutine Convert\_DryKgKg\_to\_MoistKgKg converts the units 
!  of tracer mass mixing ratio from tracer mass per dry air mass [kg tracer/
!  kg dry air] to tracer mass per moist air mass [kg tracer/kg moist air]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_DryKgKg_to_MoistKgKg( am_I_Root, N_TRACERS,  &
                                           State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: N_TRACERS   ! Number of tracers
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  16 Apr 2015 - E. Lundgren - Initial version
!  10 Jul 2015 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !====================================================================
    ! Convert_DryKgKg_to_MoistKgKg begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg tracer(N)     kg dry air     
    !   ------------  *  ------------    
    !   kg dry air       kg moist air            
    !
    !   = kg tracer(N) / kg moist air
    !
    ! Therefore, with kg dry air / kg moist air defined as the
    !  complement of specific humidity (kg water vapor / kg moist air),
    !  the conversion is:
    ! 
    ! TRACERS(I,J,L,N) [kg/kg moist]
    !
    !    = TRACERS(I,J,L,N) [kg/kg dry] * ( 1 - SPHU(I,J,L) )
    !
    ! Note that State_Met%SPHU is in units of [g/kg] and so must
    ! be converted to [kg/kg]. 
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)        &
                                  * ( 1.0e+0_fp - State_Met%SPHU(I,J,L) &
                                  * 1.e-3_fp )
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Convert_DryKgKg_to_MoistKgKg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_moistkgkg_to_drykgkg
!
! !DESCRIPTION: Subroutine Convert\_MoistKgKg\_to\_DryKgKg converts the units 
!  of tracer mass mixing ratio from tracer mass per moist air mass [kg tracer/
!  kg moist air] to tracer mass per dry air mass [kg tracer/kg dry air]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_MoistKgKg_to_DryKgKg( am_I_Root, N_TRACERS, &
                                           State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: N_TRACERS   ! Number of tracers
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  16 Apr 2015 - E. Lundgren - Initial version
!  10 Jul 2015 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !====================================================================
    ! Convert_MoistKgKg_to_DryKgKg begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg tracer(N)     kg moist air     
    !   ------------  *  ------------    
    !   kg moist air     kg dry air            
    !
    !   = kg tracer(N) / kg dry air
    !
    ! Therefore, with kg dry air / kg moist air defined as the
    !  complement of specific humidity (kg water vapor / kg moist air),
    !  the conversion is:
    ! 
    ! TRACERS(I,J,L,N) [kg/kg dry]
    !
    !    = TRACERS(I,J,L,N) [kg/kg moist] / ( 1 - SPHU(I,J,L) )
    !      
    ! Note that State_Met%SPHU is in units of [g/kg] and so must
    ! be converted to [kg/kg]. 
    !                                
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)          &
                                  / ( 1.0e+0_fp - State_Met%SPHU(I,J,L) &
                                  * 1.e-3_fp )
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO
  
  END SUBROUTINE Convert_MoistKgKg_to_DryKgKg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_dryvv_to_moistvv
!
! !DESCRIPTION: Subroutine Convert\_DryVV\_to\_MoistVV converts the units of 
!  tracer volume mixing ratio from vol tracer per vol dry air [mol tracer/
!  mol dry air] to vol tracer per vol moist air [mol tracer/mol moist air]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_DryVV_to_MoistVV( am_I_Root, N_TRACERS, &
                                       State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: N_TRACERS   ! Number of tracers
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  16 Apr 2015 - E. Lundgren - Initial version
!  10 Jul 2015 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !====================================================================
    ! Convert_DryVV_to_MoistVV begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   mol tracer(N)     mol dry air     
    !   ------------  *  ------------    
    !   mol dry air       mol moist air            
    !
    !   = mol tracer(N) / mol moist air
    !
    ! Therefore, with mol dry air / mol moist air defined as the ratio
    !  of dry partial pressure to moist pressure, the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [v/v moist]
    !
    !    = TRACERS(I,J,L,N) [v/v dry] * State_Met%PMID_DRY(I,J,L) &   
    !                              / State_Met%PMID(I,J,L)
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N) &
                                  * State_Met%PMID_DRY(I,J,L)  &
                                  / State_Met%PMID(I,J,L)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Convert_DryVV_to_MoistVV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_moistvv_to_dryvv
!
! !DESCRIPTION: Subroutine Convert\_MoistVV\_to\_DryVV converts the units of 
!  tracer volume mixing ratio from vol tracer per vol moist air [mol tracer/
!  mol moist air] to mol tracer per mol dry air [mol tracer/mol dry air]. 
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE Convert_MoistVV_to_DryVV( am_I_Root, N_TRACERS, &
                                         State_Met, State_Chm, RC ) 
!
! !USES:
!
  USE GIGC_State_Met_Mod, ONLY : MetState
  USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: N_TRACERS   ! Number of tracers
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  16 Apr 2015 - E. Lundgren - Initial version
!  10 Jul 2015 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !====================================================================
    ! Convert_MoistVV_to_DryVV begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   mol tracer(N)     mol moist air     
    !   ------------  *  ------------    
    !   mol moist air     mol dry air            
    !
    !   = mol tracer(N) / mol dry air
    !
    ! Therefore, with mol moist air / mol dry air defined as the ratio
    !  of moist pressure to dry air partial pressure, the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [v/v dry]
    !
    !    = TRACERS(I,J,L,N) [v/v moist] * State_Met%PMID(I,J,L)
    !                               / State_Met%PMID_DRY(I,J,L)
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N) &
                                  * State_Met%PMID(I,J,L)      &
                                  / State_Met%PMID_DRY(I,J,L)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Convert_MoistVV_to_DryVV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_drykgkg_to_dryvv
!
! !DESCRIPTION: Subroutine Convert\_DryKgKg\_to\_DryVV converts the units of 
!  tracer concentrations from mass mixing ratio (KGKG) [kg/kg] to 
!  volume ratio (VR) [vol/vol] (same as molar ratio [mol/mol]). 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_DryKgKg_to_DryVV( am_I_Root, N_TRACERS, Input_Opt, &
                                       State_Chm, RC ) 
!
! USES: 
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: N_TRACERS   ! Number of tracers
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!  10 Jul 2015 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    
    !====================================================================
    ! Convert_DryKgKg_to_DryVV begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg tracer(N)    g dry air      mol tracer(N)    
    !   -----------  * ----------  *  -------------  
    !     kg air       mol air         g tracer(N)          
    !
    !   = mass mixing ratio * ratio of air to tracer molecular weights  
    !   
    !   = molar ratio
    !
    ! Therefore, with:
    !
    !  TCVV(N) = dry air molecular wt / tracer molecular wt 
    !     
    ! the conversion is:
    ! 
    !  Tracers(I,J,L,N) [vol/vol]
    !
    !    = Tracers(I,J,L,N) [kg/kg] * TCVV(N)
    !                   
    !====================================================================
 
    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%Tracers(I,J,L,N) = State_Chm%Tracers(I,J,L,N)  &
                                  * Input_Opt%TCVV(N)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Convert_DryKgKg_to_DryVV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_dryvv_to_drykgkg
!
! !DESCRIPTION: Subroutine Convert\_DryVV\_to\_DryKgKg converts the units of 
!  tracer concentrations from volume ratio (VR) [vol/vol] (same 
!  as molar mixing ratio [mol/mol]) to mass mixing ratio [kg/kg]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_DryVV_to_DryKgKg( am_I_Root, N_TRACERS, Input_Opt, &
                                       State_Chm, RC ) 
!
! USES: 
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: N_TRACERS   ! Number of tracers
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!  10 Jul 2015 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !====================================================================
    ! Convert_DryVV_to_DryKgKg begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   mol tracer(N)  mol air     g tracer(N)         
    !   -----------  * -------  *  -------------  
    !     mol air       g air      mol tracer(N)           
    !
    !   = volume ratio / ratio of air to tracer molecular wts  
    !   
    !   = mass mixing ratio ([g/g] is equivalent to [kg/kg])
    !
    ! Therefore, with:
    !
    !  TCVV(N) = dry air molecular wt / tracer molecular wt 
    !     
    ! the conversion is:
    ! 
    !  Tracers(I,J,L,N) [vol/vol]
    !
    !    = Tracers(I,J,L,N) [kg/kg] / TCVV(N)
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%Tracers(I,J,L,N) = State_Chm%Tracers(I,J,L,N)   &
                                  / Input_Opt%TCVV(N)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Convert_DryVV_to_DryKgKg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_drymmr_to_mnd
!
! !DESCRIPTION: Subroutine Convert\_DryKgKg\_to\_MND converts the units of 
!  tracer concentrations from dry mass mixing ratio [kg/kg dry air] to 
!  molecular number density (MND) [molecules/cm3].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_DryKgKg_to_MND( am_I_Root, N_TRACERS, State_Met, &
                                     Input_Opt, State_Chm, RC          )
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: N_TRACERS   ! Number of tracers
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!  10 Jul 2015 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !====================================================================
    ! Convert_DryKgKg_to_MND begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg tracer(N)    kg air          m3        molecules tracer(N)     
    !   -----------  * --------  *  ---------  *  ------------------      
    !   kg dry air        m3         1E6 cm3        kg tracer(N)  
    !
    !   = mass mixing ratio * dry air density * molecules / kg tracer
    !   
    !   = molecules per cm3
    !
    ! Therefore, with:
    !
    !  XNUMOL(N)       = molecules tracer / kg tracer
    !  AIRDEN(I,J,L)   = grid box dry air density [kg/m3]
    !     
    ! the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [molecules/cm3]
    !
    !    = TRACERS(I,J,L,N) [kg/kg] * XNUMOL(N) * AIRDEN(I,J,L) * 1e-6
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)        &
                                  * Input_Opt%XNUMOL(N)               &
                                  * State_Met%AIRDEN(I,J,L) * 1E-6_fp  
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Convert_DryKgKg_to_MND
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_mnd_to_drymmr
!
! !DESCRIPTION: Subroutine Convert\_MND\_to\_DryKgKg converts the units of 
!  tracer concentrations from molecular number density (MND)
!  [molecules/cm3] to dry mass mixing ratio [kg/kg dry air].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_MND_to_DryKgKg( am_I_Root, N_TRACERS, State_Met, &
                                     Input_Opt, State_Chm, RC          )
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: N_TRACERS   ! Number of tracers
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Object containing tracer concentration
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!  10 Jul 2015 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !====================================================================
    ! Convert_MND_to_DryKgKg begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   molecules tracer(N)     m3       1E6 cm3       kg tracer(N)
    !   ------------------  * ------  *  -------  * ------------------
    !          cm3            kg air        m3      molecules tracer(N)  
    !
    !   = molecules per vol / air density / (molecules / kg tracer) 
    !   
    !   = mass mixing ratio
    !
    ! Therefore, with:
    !
    !  XNUMOL(N)      = molecules tracer / kg tracer
    !  AIRDEN(I,J,L)  = grid box dry air density [kg/m3]
    !     
    ! the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [kg/kg dry air]
    !
    !   = TRACERS(I,J,L,N) [molecules/cm3] * 1E+6 
    !                          / ( XNUMOL(N) * AIRDEN(I,J,L) ) 
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)       &
                                  * 1E+6_fp / ( Input_Opt%XNUMOL(N)  &
                                  * State_Met%AIRDEN(I,J,L) ) 
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Convert_MND_to_DryKgKg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_moistmmr_to_kg
!
! !DESCRIPTION: Subroutine Convert\_MoistKgKg\_to\_Kg converts the units of 
!  tracer concentrations from moist mass mixing ratio 
!  [kg tracer/kg total (wet) air] to tracer mass per grid box [kg]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_MoistKgKg_to_Kg( am_I_Root, N_TRACERS,    &
                                      State_Met, State_Chm, RC   ) 
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: N_TRACERS   ! Number of tracers
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Object containing tracer concentration
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!  10 Jul 2015 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !=================================================================
    ! Convert_MoistKgKg_to_Kg begins here!
    !=================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !==============================================================
    !
    !  The conversion is as follows:
    !
    !   kg tracer(N)            
    !   -----------  *  kg wet air       
    !   kg wet air                   
    !
    !   = mass mixing ratio * wet air mass  
    !   
    !   = kg tracer(N)
    !
    ! Therefore, with:
    !
    !  AIRMASS(I,J,L)   = grid box total air mass [kg] (wet)
    !     
    ! the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [kg]
    !
    !    = TRACERS(I,J,L,N) [kg/kg] * AIRMASS(I,J,L)
    !                   
    !==============================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)  &
                                  * State_Met%ADMOIST(I,J,L)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Convert_MoistKgKg_to_Kg
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kg_to_moistmmr
!
! !DESCRIPTION: Subroutine Convert\_Kg\_to\_MoistKgKg converts the units of 
!  tracer concentrations from tracer mass per grid box [kg] to mass 
!  mixing ratio [kg tracer/kg total (wet) air]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_Kg_to_MoistKgKg( am_I_Root, N_TRACERS,    &
                                      State_Met, State_Chm, RC   ) 
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: N_TRACERS   ! Number of tracers
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Object containing tracer concentration
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!  10 Jul 2015 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !====================================================================
    ! Convert_Kg_to_MoistKgKg begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !                         1          
    !   kg tracer(N)  * --------------      
    !                    kg total air              
    !
    !   = kg tracer(N) / total air mass
    !   
    !   = mass mixing ratio
    !
    ! Therefore, with:
    !
    !  AIRMASS(I,J,L)    = grid box total air mass [kg] (wet)
    !     
    ! the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [kg/kg]
    !
    !    = TRACERS(I,J,L,N) [kg] / AIRMASS(I,J,L) 
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N) &
                                  / State_Met%ADMOIST(I,J,L)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Convert_Kg_to_MoistKgKg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_drymmr_to_kg
!
! !DESCRIPTION: Subroutine Convert\_DryKgKg\_to\_Kg converts the units of 
!  tracer concentrations from dry mass mixing ratio 
!  [kg tracer/kg dry air] to tracer mass per grid box [kg]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_DryKgKg_to_Kg( am_I_Root, N_TRACERS,    &
                                    State_Met, State_Chm, RC   ) 
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: N_TRACERS   ! Number of tracers
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Object containing tracer concentration
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!  10 Jul 2015 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !====================================================================
    ! Convert_DryKgKg_to_Kg begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg tracer(N)            
    !   -----------  *  kg dry air       
    !   kg dry air                   
    !
    !   = mass mixing ratio * dry air mass  
    !   
    !   = kg tracer(N)
    !
    ! Therefore, with:
    !
    !  AIRMASS(I,J,L)   = grid box dry air mass [kg]
    !     
    ! the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [kg]
    !
    !    = TRACERS(I,J,L,N) [kg/kg] * AIRMASS(I,J,L)
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N) &
                                  * State_Met%AD(I,J,L)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Convert_DryKgKg_to_Kg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kg_to_drymmr
!
! !DESCRIPTION: Subroutine Convert\_Kg\_to\_DryKgKg converts the units of 
!  tracer concentrations from tracer mass per grid box [kg] to dry mass 
!  mixing ratio [kg tracer/kg dry air].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_Kg_to_DryKgKg( am_I_Root, N_TRACERS,    &
                                    State_Met, State_Chm, RC   ) 
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: N_TRACERS   ! Number of tracers
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!  10 Jul 2015 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !====================================================================
    ! Convert_Kg_to_DryKgKg begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !                         1          
    !   kg tracer(N)  * --------------      
    !                    kg dry air              
    !
    !   = kg tracer(N) / dry air mass
    !   
    !   = mass mixing ratio
    !
    ! Therefore, with:
    !
    !  AIRMASS(I,J,L)    = grid box dry air mass [kg]
    !     
    ! the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [kg/kg]
    !
    !    = TRACERS(I,J,L,N) [kg] / AIRMASS(I,J,L) 
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N) &
                                  / State_Met%AD(I,J,L)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Convert_Kg_to_DryKgKg
!EOC
END MODULE UnitConv_Mod
