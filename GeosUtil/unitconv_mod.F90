!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: unitconv_mod.F90
!
! !DESCRIPTION: Module UNITCONV\_MOD contains routines which are used to 
!  convert the units of tracer concentrations between mass mixing ratio 
!  [kg/kg air], mass per grid box per area [kg/m2], molar mixing ratio 
!  [vol/vol], and molecular number density [molecules/cm3]. There are 
!  different conversion routines for dry air and total (wet) air mixing 
!  ratios. Conversions involving column area will be phased out for
!  grid-independent GEOS-Chem. 
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
  USE ERROR_MOD
                    
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! State_Chm%TRACERS: TOTAL <-> DRY 
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! kg/kg dry air <-> kg/kg total air
  PUBLIC  :: Convert_KgKgDry_to_KgKgTotal
  PUBLIC  :: Convert_KgKgTotal_to_KgKgDry

  ! v/v dry air <-> v/v total air (not currently used)
  PUBLIC  :: Convert_VVDry_to_VVTotal
  PUBLIC  :: Convert_VVTotal_to_VVDry

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! State_Chm%TRACERS: KG/KG <-> V/V
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! kg/kg dry air <-> v/v dry air
  PUBLIC  :: Convert_KgKgDry_to_VVDry
  PUBLIC  :: Convert_VVDry_to_KgKgDry

  ! kg/kg total air <-> v/v total air
  PUBLIC  :: Convert_KgKgTotal_to_VVTotal
  PUBLIC  :: Convert_VVTotal_to_KgKgTotal

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! State_Chm%TRACERS: MIXING RATIO <-> AREA DENSITY
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! kg/kg dry air <-> kg/m2
  PUBLIC  :: Convert_KgKgDry_to_Kgm2
  PUBLIC  :: Convert_kgm2_to_KgKgDry

  !%%%%%%%%%%%%%%%%
  ! AREA-DEPENDENT
  !%%%%%%%%%%%%%%%%

  ! kg/kg dry air <-> molec/cm3 dry air
  ! (not currently used, could be adapted to be area-independent)
  PUBLIC  :: Convert_KgKgDry_to_MND
  PUBLIC  :: Convert_MND_to_KgKgDry

  ! kg/kg moist air <-> kg/grid box (State_Chm%TRACERS)
  PUBLIC  :: Convert_KgKgTotal_to_Kg
  PUBLIC  :: Convert_Kg_to_KgKgTotal

  ! kg/kg dry air <-> kg/grid box (State_Chm%TRACERS)
  PUBLIC  :: Convert_KgKgDry_to_Kg
  PUBLIC  :: Convert_Kg_to_KgKgDry

  ! kg <-> kg/m2 (State_Chm%TRACERS)
  ! (area-dependent, eventually remove)
  PUBLIC  :: Convert_Kg_to_Kgm2
  PUBLIC  :: Convert_Kgm2_to_Kg 
!
! !REMARKS:
!  The routines in this module are used to convert the units of tracer 
!  concentrations in various GEOS-Chem routines.
!
! !REVISION HISTORY:
!  23 Jun 2015 - E. Lundgren - Initial version
!  13 Aug 2015 - E. Lundgren - Add tracer unit error handling
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kgkgdry_to_kgkgtotal
!
! !DESCRIPTION: Subroutine Convert\_KgKgDry\_to\_KgKgTotal converts the units 
!  of tracer mass mixing ratio from tracer mass per dry air mass [kg tracer/
!  kg dry air] to tracer mass per moist air mass [kg tracer/kg moist air]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_KgKgDry_to_KgKgTotal( am_I_Root, N_TRACERS,  &
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_KgKgDry_to_KgKgTotal begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'kg/kg dry' ) THEN
       MSG = 'Incorrect initial units:' // TRIM( State_Chm%Trac_Units )
       LOC = 'UNITCONV_MOD: Convert_KgKgDry_to_KgKgTotal'
       CALL GIGC_Error( MSG, RC, LOC)
    ENDIF

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

    ! Update tracer units
    State_Chm%Trac_Units = 'kg/kg total'

  END SUBROUTINE Convert_KgKgDry_to_KgKgTotal
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kgkgtotal_to_kgkgdry
!
! !DESCRIPTION: Subroutine Convert\_KgKgTotal\_to\_KgKgDry converts the units 
!  of tracer mass mixing ratio from tracer mass per moist air mass [kg tracer/
!  kg moist air] to tracer mass per dry air mass [kg tracer/kg dry air]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_KgKgTotal_to_KgKgDry( am_I_Root, N_TRACERS, &
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_KgKgTotal_to_KgKgDry begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'kg/kg total' ) THEN
       MSG = 'Incorrect initial units:' // TRIM( State_Chm%Trac_Units )
       LOC = 'UNITCONV_MOD: Convert_KgKgTotal_to_KgKgDry'
       CALL GIGC_Error( MSG, RC, LOC)
       RETURN
    ENDIF

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

    ! Update tracer units
    State_Chm%Trac_Units = 'kg/kg dry'
  
  END SUBROUTINE Convert_KgKgTotal_to_KgKgDry
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_vvdry_to_vvtotal
!
! !DESCRIPTION: Subroutine Convert\_VVDry\_to\_VVTotal converts the units of 
!  tracer volume mixing ratio from vol tracer per vol dry air [mol tracer/
!  mol dry air] to vol tracer per vol moist air [mol tracer/mol moist air]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_VVDry_to_VVTotal( am_I_Root, N_TRACERS, &
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_VVDry_to_VVTotal begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'v/v dry' ) THEN
       MSG = 'Incorrect initial units:' // TRIM( State_Chm%Trac_Units )
       LOC = 'UNITCONV_MOD: Convert_VVDry_to_VVTotal'
       CALL GIGC_Error( MSG, RC, LOC)
       RETURN
    ENDIF

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

    ! Update tracer units
    State_Chm%Trac_Units = 'v/v total'

  END SUBROUTINE Convert_VVDry_to_VVTotal
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_vvtotal_to_vvdry
!
! !DESCRIPTION: Subroutine Convert\_VVTotal\_to\_VVDry converts the units of 
!  tracer volume mixing ratio from vol tracer per vol moist air [mol tracer/
!  mol moist air] to mol tracer per mol dry air [mol tracer/mol dry air]. 
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE Convert_VVTotal_to_VVDry( am_I_Root, N_TRACERS, &
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_VVTotal_to_VVDry begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'v/v total' ) THEN
       MSG = 'Incorrect initial units:' // TRIM( State_Chm%Trac_Units )
       LOC = 'UNITCONV_MOD: Convert_VVTotal_to_VVDry'
       CALL GIGC_Error( MSG, RC, LOC)
       RETURN
    ENDIF

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

    ! Update tracer units
    State_Chm%Trac_Units = 'v/v dry'

  END SUBROUTINE Convert_VVTotal_to_VVDry
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kgkgdry_to_vvdry
!
! !DESCRIPTION: Subroutine Convert\_KgKgDry\_to\_VVDry converts the units of 
!  tracer concentrations from mass mixing ratio (KGKG) [kg/kg] to 
!  volume ratio (VR) [vol/vol] (same as molar ratio [mol/mol]). 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_KgKgDry_to_VVDry( am_I_Root, N_TRACERS, Input_Opt, &
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC
    
    !====================================================================
    ! Convert_KgKgDry_to_VVDry begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'kg/kg dry' ) THEN
       MSG = 'Incorrect initial units:' // TRIM( State_Chm%Trac_Units )
       LOC = 'UNITCONV_MOD: Convert_KgKgDry_to_VVDry'
       CALL GIGC_Error( MSG, RC, LOC)
       RETURN
    ENDIF

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

    ! Update tracer units
    State_Chm%Trac_Units = 'v/v dry'

  END SUBROUTINE Convert_KgKgDry_to_VVDry
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_vvdry_to_kgkgdry
!
! !DESCRIPTION: Subroutine Convert\_VVDry\_to\_KgKgDry converts the units of 
!  tracer concentrations from volume ratio (VR) [vol/vol] (same 
!  as molar mixing ratio [mol/mol]) to mass mixing ratio [kg/kg]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_VVDry_to_KgKgDry( am_I_Root, N_TRACERS, Input_Opt, &
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_VVDry_to_KgKgDry begins here!
    !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'v/v dry' ) THEN
       MSG = 'Incorrect initial units:' // TRIM( State_Chm%Trac_Units )
       LOC = 'UNITCONV_MOD: Convert_VVDry_to_KgKgDry'
       CALL GIGC_Error( MSG, RC, LOC)
       RETURN
    ENDIF

         !==============================================================
         !
         !  The conversion is as follows:
         !
         !   mol tracer(N)  mol dry air     g tracer(N)         
         !   -----------  * -----------  *  -------------  
         !     mol air       g dry air      mol tracer(N)           
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
         !==============================================================

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

    ! Update tracer units
    State_Chm%Trac_Units = 'kg/kg dry'

    END SUBROUTINE Convert_VVDry_to_KgKgDry
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kgkgtotal_to_vvtotal
!
! !DESCRIPTION: Subroutine Convert\_KgKgTotal\_to\_VVTotal converts the units 
!  of tracer concentrations from mass mixing ratio (KGKG) [kg/kg] to 
!  volume ratio (VR) [vol/vol] (same as molar ratio [mol/mol]). 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_KgKgTotal_to_VVTotal( am_I_Root, N_TRACERS, Input_Opt, &
                                           State_Met, State_Chm, RC ) 
!
! USES: 
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: N_TRACERS   ! Number of tracers
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
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
!  08 Sep 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC
    
    !====================================================================
    ! Convert_KgKgTotal_to_VVTotal begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'kg/kg total' ) THEN
       MSG = 'Incorrect initial units:' // TRIM( State_Chm%Trac_Units )
       LOC = 'UNITCONV_MOD: Convert_KgKgTotal_to_VVTotal'
       CALL GIGC_Error( MSG, RC, LOC)
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg tracer(N)   g total air     mol tracer(N)    
    !   -----------  * ----------  *  -------------  
    !     kg air       mol air         g tracer(N)          
    !
    !   = mass mixing ratio * ratio of air to tracer molecular weights  
    !   
    !   = molar ratio
    !
    ! Therefore, with:
    !
    !  Tracer_MW_G(N)   = tracer molecular wt 
    !  MoistMW(I,J,L) = total air molecular wt for grid box (I,J,L)  
    !
    ! the conversion is:
    ! 
    !  Tracers(I,J,L,N) [vol/vol]
    !
    !    = Tracers(I,J,L,N) [kg/kg] * MoistMW(I,J,L) / Tracer_MW(N)
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
                                    * State_Met%MoistMW(I,J,L)  &
                                    / Input_Opt%Tracer_MW_G(N)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Update tracer units
    State_Chm%Trac_Units = 'v/v total'

  END SUBROUTINE Convert_KgKgTotal_to_VVTotal
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_vvtotal_to_kgkgtotal
!
! !DESCRIPTION: Subroutine Convert\_VVTotal\_to\_KgKgTotal converts the 
!  units of tracer concentrations from volume ratio (VR) [vol/vol] (same 
!  as molar mixing ratio [mol/mol]) to mass mixing ratio [kg/kg]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_VVTotal_to_KgKgTotal( am_I_Root, N_TRACERS, Input_Opt, &
                                           State_Met, State_Chm, RC ) 
!
! USES: 
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: N_TRACERS   ! Number of tracers
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_VVTotal_to_KgKgTotal begins here!
    !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'v/v total' ) THEN
       MSG = 'Incorrect initial units:' // TRIM( State_Chm%Trac_Units )
       LOC = 'UNITCONV_MOD: Convert_VVTotal_to_KgKgTotal'
       CALL GIGC_Error( MSG, RC, LOC)
       RETURN
    ENDIF

         !==============================================================
         !
         !  The conversion is as follows:
         !
         !   mol tracer(N)  mol total air     g tracer(N)         
         !   -----------  * -------------  *  -------------  
         !     mol air       g total air      mol tracer(N)           
         !
         !   = volume ratio / ratio of air to tracer molecular wts  
         !   
         !   = mass mixing ratio ([g/g] is equivalent to [kg/kg])
         !
         ! Therefore, with:
         !
         !  Tracer_MW_G(N) = tracer molecular wt 
         !  MoistMW(I,J,L) = total air molecular wt for grid box (I,J,L)  
         !
         ! the conversion is:
         ! 
         !  Tracers(I,J,L,N) [vol/vol]
         !
         !    = Tracers(I,J,L,N) [kg/kg] * Tracer_MW(N) / MoistMW(I,J,L) 
         !                   
         !==============================================================

      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( I, J, L, N ) 
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
        State_Chm%Tracers(I,J,L,N) = State_Chm%Tracers(I,J,L,N)   &
                                    * Input_Opt%Tracer_MW_G(N)    &
                                    / State_Met%MoistMW(I,J,L)  

      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ! Update tracer units
    State_Chm%Trac_Units = 'kg/kg total'

    END SUBROUTINE Convert_VVTotal_to_KgKgTotal
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kgkgdry_to_kgm2
!
! !DESCRIPTION: Subroutine Convert\_kgkgdry\_to\_kgm2 converts the units of 
!  a 3D array from dry mass mixing ratio [kg/kg dry air] to area density 
!  [kg/m2].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_KgKgDry_to_Kgm2( am_I_Root, N_TRACERS, State_Met,  &
                                      State_Chm, RC             )
!
! !USES:
!
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE CMN_GCTM_MOD

!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: N_TRACERS   ! Number of tracers
    TYPE(MetState), INTENT(IN)    :: State_Met     ! Meteorology state object

!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC            ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  14 Aug 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC
    REAL(fp)           :: SPHU_kgkg

    !====================================================================
    ! Convert_KgKgDry_to_Kgm2 begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'kg/kg dry' ) THEN
       MSG = 'Incorrect initial units:' // TRIM( State_Chm%Trac_Units )
       LOC = 'UNITCONV_MOD: Convert_KgKgDry_to_Kgm2'
       CALL GIGC_Error( MSG, RC, LOC)
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg tracer      Delta P [hPa]   100 [Pa]     kg dry air 
    !   -----------  * ------------- * -------- * --------------     
    !   kg dry air     g [m/s2]        [hPa]      kg total air  
    !
    !   = kg / m2
    !
    ! where:
    !
    !  Delta P = edge pressure difference across level
    !  g = acceleration due to gravity
    !  kg dry air / kg total air  = 1 - specific humidity
    !     
    !====================================================================

    !$OMP PARALLEL DO        &
    !$OMP DEFAULT( SHARED  ) &
    !$OMP PRIVATE( I, J, L, SPHU_kgkg ) 
    DO N = 1, N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Convert specific humidity from [g/kg] to [kg/kg]
       SPHU_kgkg = State_Met%SPHU(I,J,L) * 1.0e-3_fp 

       ! Area-independent conversion results in +/-1e-5% precision diffs
       ! while equivalent area-dependent conversion below does not
!       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)         &
!                                    * ( ( 1.0e+0_fp - SPHU_kgkg )      &
!                                    * ( g0_100                         &
!                                    * State_Met%DELP(I,J,L) ) )          

!       ! Equivalent area-dependent conversion
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)          &
                                    * ( State_Met%AD(I,J,L)             &
                                        / State_Met%AREA_M2(I,J,1) ) 

    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Update tracer units
    State_Chm%Trac_Units = 'kg/m2'

  END SUBROUTINE Convert_KgKgDry_to_Kgm2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kgm2_to_kgkgdry
!
! !DESCRIPTION: Subroutine Convert\_Kgm2\_to\_kgkgdry converts the units of 
!  tracer concentrations from area density [kg/m2] to dry mass mixing ratio 
!  [kg/kg dry air].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_Kgm2_to_KgKgDry( am_I_Root, N_TRACERS, State_Met, &
                                      State_Chm, RC          )
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE CMN_GCTM_MOD
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
!  14 Aug 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC
    REAL(fp)           :: SPHU_kgkg

    !====================================================================
    ! Convert_Kgm2_to_KgKgDry begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'kg/m2' ) THEN
       MSG = 'Incorrect initial units:' // TRIM( State_Chm%Trac_Units )
       LOC = 'UNITCONV_MOD: Convert_Kgm2_to_KgKgDry'
       CALL GIGC_Error( MSG, RC, LOC)
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg tracer(N)   g [m/s2]        [hPa]      kg total air   
    !   -----------  * ------------- * -------- * --------------     
    !        m2        Delta P [hPa]   100 [Pa]    kg dry air
    !
    !   = kg tracer(N) / kg dry air
    !
    ! where:
    !
    !  Delta P = edge pressure difference across level
    !  g = acceleration due to gravity
    !  kg dry air / kg total air  = 1 - specific humidity
    !     
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N, SPHU_kgkg ) 
    DO N = 1, N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Convert specific humidity from [g/kg] to [kg/kg]
       SPHU_kgkg = State_Met%SPHU(I,J,L) * 1.0e-3_fp 

       ! Area-independent conversion results in +/-1e-5% precision diffs
       ! while equivalent area-dependent conversion below does not
!       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)          &
!                                    * ( 1.0e+0_fp                       &      
!                                    / ( ( 1.0e+0_fp - SPHU_kgkg )       &
!                                    * ( g0_100                          &
!                                    * State_Met%DELP(I,J,L) ) ) )          

       ! Equivalent area-dependent conversion
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)          &
                                  / ( State_Met%AD(I,J,L)               &
                                      / State_Met%AREA_M2(I,J,1) )         
                               

    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Update tracer units
    State_Chm%Trac_Units = 'kg/kg dry'

  END SUBROUTINE Convert_Kgm2_to_KgKgDry
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kgkgdry_to_mnd
!
! !DESCRIPTION: Subroutine Convert\_KgKgDry\_to\_MND converts the units of 
!  tracer concentrations from dry mass mixing ratio [kg/kg dry air] to 
!  molecular number density (MND) [molecules/cm3].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_KgKgDry_to_MND( am_I_Root, N_TRACERS, State_Met, &
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_KgKgDry_to_MND begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'kg/kg dry' ) THEN
       MSG = 'Incorrect initial units:' // TRIM( State_Chm%Trac_Units )
       LOC = 'UNITCONV_MOD: Convert_KgKgDry_to_MND'
       CALL GIGC_Error( MSG, RC, LOC)
       RETURN
    ENDIF

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

    ! Update tracer units
    State_Chm%Trac_Units = 'molec/cm3'

  END SUBROUTINE Convert_KgKgDry_to_MND
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_mnd_to_kgkgdry
!
! !DESCRIPTION: Subroutine Convert\_MND\_to\_KgKgDry converts the units of 
!  tracer concentrations from molecular number density (MND)
!  [molecules/cm3] to dry mass mixing ratio [kg/kg dry air].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_MND_to_KgKgDry( am_I_Root, N_TRACERS, State_Met, &
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_MND_to_KgKgDry begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'molec/cm3' ) THEN
       MSG = 'Incorrect initial units:' // TRIM( State_Chm%Trac_Units )
       LOC = 'UNITCONV_MOD: Convert_MND_to_KgKgDry'
       CALL GIGC_Error( MSG, RC, LOC)
       RETURN
    ENDIF

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

    ! Update tracer units
    State_Chm%Trac_Units = 'kg/kg dry'

  END SUBROUTINE Convert_MND_to_KgKgDry
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kgkgtotal_to_kg
!
! !DESCRIPTION: Subroutine Convert\_KgKgTotal\_to\_Kg converts the units of 
!  tracer concentrations from moist mass mixing ratio 
!  [kg tracer/kg total (wet) air] to tracer mass per grid box [kg]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_KgKgTotal_to_Kg( am_I_Root, N_TRACERS,    &
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !=================================================================
    ! Convert_KgKgTotal_to_Kg begins here!
    !=================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'kg/kg total' ) THEN
       MSG = 'Incorrect initial units:' // TRIM( State_Chm%Trac_Units )
       LOC = 'UNITCONV_MOD: Convert_KgKgTotal_to_Kg'
       CALL GIGC_Error( MSG, RC, LOC)
       RETURN
    ENDIF

    !==============================================================
    !
    !  The conversion is as follows:
    !
    !   kg tracer(N)            
    !   -----------  *  kg total air       
    !   kg total air                   
    !
    !   = mass mixing ratio * total air mass  
    !   
    !   = kg tracer(N)
    !
    ! Therefore, with:
    !
    !  ADMOIST(I,J,L)   = grid box total air mass [kg] (wet)
    !     
    ! the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [kg]
    !
    !    = TRACERS(I,J,L,N) [kg/kg] * ADMOIST(I,J,L)
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

    ! Update tracer units
    State_Chm%Trac_Units = 'kg'

  END SUBROUTINE Convert_KgKgTotal_to_Kg
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kg_to_kgkgtotal
!
! !DESCRIPTION: Subroutine Convert\_Kg\_to\_KgKgTotal converts the units of 
!  tracer concentrations from tracer mass per grid box [kg] to mass 
!  mixing ratio [kg tracer/kg total (wet) air]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_Kg_to_KgKgTotal( am_I_Root, N_TRACERS,    &
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_Kg_to_KgKgTotal begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'kg' ) THEN
       MSG = 'Incorrect initial units:' // TRIM( State_Chm%Trac_Units )
       LOC = 'UNITCONV_MOD: Convert_Kg_to_KgKgTotal'
       CALL GIGC_Error( MSG, RC, LOC)
       RETURN
    ENDIF

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
    !  ADMOIST(I,J,L)    = grid box total air mass [kg] (wet)
    !     
    ! the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [kg/kg]
    !
    !    = TRACERS(I,J,L,N) [kg] / ADMOIST(I,J,L) 
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

    ! Update tracer units
    State_Chm%Trac_Units = 'kg/kg total'

  END SUBROUTINE Convert_Kg_to_KgKgTotal
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kgkgdry_to_kg
!
! !DESCRIPTION: Subroutine Convert\_KgKgDry\_to\_Kg converts the units of 
!  tracer concentrations from dry mass mixing ratio 
!  [kg tracer/kg dry air] to tracer mass per grid box [kg]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_KgKgDry_to_Kg( am_I_Root, N_TRACERS,    &
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_KgKgDry_to_Kg begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'kg/kg dry' ) THEN
       MSG = 'Incorrect initial units:' // TRIM( State_Chm%Trac_Units )
       LOC = 'UNITCONV_MOD: Convert_KgKgDry_to_Kg'
       CALL GIGC_Error( MSG, RC, LOC)
       RETURN
    ENDIF

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
    !  AD(I,J,L)   = grid box dry air mass [kg]
    !     
    ! the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [kg]
    !
    !    = TRACERS(I,J,L,N) [kg/kg] * AD(I,J,L)
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

    ! Update tracer units
    State_Chm%Trac_Units = 'kg'

  END SUBROUTINE Convert_KgKgDry_to_Kg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kg_to_kgkgdry
!
! !DESCRIPTION: Subroutine Convert\_Kg\_to\_KgKgDry converts the units of 
!  tracer concentrations from tracer mass per grid box [kg] to dry mass 
!  mixing ratio [kg tracer/kg dry air].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_Kg_to_KgKgDry( am_I_Root, N_TRACERS,    &
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_Kg_to_KgKgDry begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'kg' ) THEN
       MSG = 'Incorrect initial units:' // TRIM( State_Chm%Trac_Units )
       LOC = 'UNITCONV_MOD: Convert_Kg_to_KgKgDry'
       CALL GIGC_Error( MSG, RC, LOC)
       RETURN
    ENDIF

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
    !  AD(I,J,L)    = grid box dry air mass [kg]
    !     
    ! the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [kg/kg]
    !
    !    = TRACERS(I,J,L,N) [kg] / AD(I,J,L) 
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

    ! Update tracer units
    State_Chm%Trac_Units = 'kg/kg dry'

  END SUBROUTINE Convert_Kg_to_KgKgDry 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kg_to_kgm2
!
! !DESCRIPTION: Subroutine Convert\_Kg\_to\_kgm2 converts the units of 
! mass [kg] to area density [kg/m2].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_Kg_to_Kgm2( am_I_Root, N_TRACERS, State_Met,  &
                                      State_Chm, RC             )
!
! !USES:
!
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE CMN_GCTM_MOD
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: N_TRACERS   ! Number of tracers
    TYPE(MetState), INTENT(IN)    :: State_Met     ! Meteorology state object

!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC            ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  14 Aug 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_Kg_to_Kgm2 begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !====================================================================

    !$OMP PARALLEL DO        &
    !$OMP DEFAULT( SHARED  ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)    &
                                     / State_Met%AREA_M2(I,J,1)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Convert_Kg_to_Kgm2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kgm2_to_kg
!
! !DESCRIPTION: Subroutine Convert\_Kgm2\_to\_Kg converts the units of area 
!  density [kg/m2] to mass [kg].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_Kgm2_to_Kg( am_I_Root, N_TRACERS, State_Met,  &
                                      State_Chm, RC             )
!
! !USES:
!
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE CMN_GCTM_MOD
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: N_TRACERS   ! Number of tracers
    TYPE(MetState), INTENT(IN)    :: State_Met     ! Meteorology state object

!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC            ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  14 Aug 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_Kgm2_to_Kg begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !====================================================================

    !$OMP PARALLEL DO        &
    !$OMP DEFAULT( SHARED  ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)    &
                                     * State_Met%AREA_M2(I,J,1)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Convert_Kgm2_to_Kg
!EOC

END MODULE UnitConv_Mod
