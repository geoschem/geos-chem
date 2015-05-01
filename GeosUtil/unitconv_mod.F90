!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: unitconv_mod.F90
!
! !DESCRIPTION: Module UNITCONV\_MOD contains routines which are used to 
!  convert the units of tracer concentrations between mass mixing ratio 
!  [kg/kg], mass per grid box [kg], molar mixing ratio [vol/vol], and molecular 
!  number density [molecules/cm3]. All mixing ratios are assumed per quantity 
!  of dry air except in the conversion between dry and moist mixing ratios.
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
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Convert_DryMMR_to_MoistMMR ! kg/kg dry air   -> kg/kg moist air
  PUBLIC  :: Convert_MoistMMR_to_DryMMR ! kg/kg moist air -> kg/kg dry air
  PUBLIC  :: Convert_DryVV_to_MoistVV   ! v/v dry air     -> v/v moist air
  PUBLIC  :: Convert_MoistVV_to_DryVV   ! v/v moist air   -> v/v dry air
  PUBLIC  :: Convert_MMR_to_VR          ! kg/kg           -> v/v
  PUBLIC  :: Convert_VR_to_MMR          ! v/v             -> kg/kg
  PUBLIC  :: Convert_MMR_to_KG          ! kg/kg           -> kg/grid box
  PUBLIC  :: Convert_KG_to_MMR          ! kg/grid box     -> kg/kg
  PUBLIC  :: Convert_MMR_to_MND         ! kg/kg           -> molec/cm3
  PUBLIC  :: Convert_MND_to_MMR         ! molec/cm3       -> kg/kg
!
! !REMARKS:
!  The routines in this module are used to convert the units of tracer 
!  concentrations in various GEOS-Chem routines. Tracer concentrations 
!  are stored primarily in units of mass mixing ratio [kg/kg] but some
!  subroutines  require tracer units in mass per grid box [kg], 
!  molar mixing ratio [vol/vol], or molecular number density [molecules/cm3]
!  for certain calculation. Use of mass per grid box will be phased out
!  for grid-independence.
!
!  Species concentrations are stored in units of molecular number density
!  and the unit conversion routines within this module may therefore be 
!  used to convert species concentration units as well.
!
! !REVISION HISTORY:
!  06 Jan 2015 - E. Lundgren - Initial version
!  04 Mar 2015 - E. Lundgren - Change conversions to use dry air quantities
!                              following air quantity updates.
!  16 Apr 2015 - E. Lundgren - Add dry <-> moist mixing ratio conversions
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_drymmr_to_moistmmr
!
! !DESCRIPTION: Subroutine Convert\_DryMMR\_to\_MoistMMR converts the units of 
!  tracer mass mixing ratio from tracer mass per dry air mass [kg tracer/
!  kg dry air] to tracer mass per moist air mass [kg tracer/kg moist air]. 
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE Convert_DryMMR_to_MoistMMR( am_I_Root, N_TRACERS, &
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


    ! Number of tracers
    INTEGER,        INTENT(IN)    :: N_TRACERS 

    ! Object containing meteorological state variables
    TYPE(MetState), INTENT(IN)    :: State_Met 
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Object containing tracer concentration [kg/kg]
    TYPE(ChmState), INTENT(INOUT) :: State_Chm  
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)    :: RC          ! Success or failure?
!
! !REMARKS
!
! !REVISION HISTORY: 
!  16 Apr 2015 - E. Lundgren - Initial version
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !=================================================================
    ! Convert_DryMMR_to_MoistMMR begins here!

         !==============================================================
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
         ! STT(I,J,L,N) [kg/kg moist]
         !
         !    = STT(I,J,L,N) [kg/kg dry] * ( 1 - SPHU(I,J,L) )
         !
         ! Note that State_Met%SPHU is in units of [g/kg] and so must
         ! be converted to [kg/kg]. 
         !                   
         !==============================================================

    !=================================================================

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

    ! Return to calling program
    END SUBROUTINE Convert_DryMMR_to_MoistMMR
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_moistmmr_to_drymmr
!
! !DESCRIPTION: Subroutine Convert\_MoistMMR\_to\_DryMMR converts the units of 
!  tracer mass mixing ratio from tracer mass per moist air mass [kg tracer/
!  kg moist air] to tracer mass per dry air mass [kg tracer/kg dry air]. 
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE Convert_MoistMMR_to_DryMMR( am_I_Root, N_TRACERS, &
                                           State_Met, State_Chm, RC ) 
!
! !USES:
!
  USE GIGC_State_Met_Mod, ONLY : MetState
  USE GIGC_State_Chm_Mod, ONLY : ChmState
  USE GIGC_ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?

    ! Number of tracers
    INTEGER,        INTENT(IN)    :: N_TRACERS 

    ! Object containing meteorological state variables
    TYPE(MetState), INTENT(IN)    :: State_Met 
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Object containing tracer concentration [v/v]
    TYPE(ChmState), INTENT(INOUT) :: State_Chm  
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS
!
! !REVISION HISTORY: 
!  16 Apr 2015 - E. Lundgren - Initial version
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !=================================================================
    ! Convert_MoistMMR_to_DryMMR begins here!

         !==============================================================
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
         ! STT(I,J,L,N) [kg/kg dry]
         !
         !    = STT(I,J,L,N) [kg/kg moist] / ( 1 - SPHU(I,J,L) )
         !      
         ! Note that State_Met%SPHU is in units of [g/kg] and so must
         ! be converted to [kg/kg]. 
         !                                
         !==============================================================

    !=================================================================

      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( I, J, L, N ) 
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
        State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)         &
                                  / ( 1.0e+0_fp - State_Met%SPHU(I,J,L) &
                                  * 1.e-3_fp )
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ! Return to calling program
    END SUBROUTINE Convert_MoistMMR_to_DryMMR
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
  USE GIGC_ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?

    ! Number of tracers
    INTEGER,        INTENT(IN)    :: N_TRACERS 

    ! Object containing meteorological state variables
    TYPE(MetState), INTENT(IN)    :: State_Met 
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Object containing tracer concentration [v/v]
    TYPE(ChmState), INTENT(INOUT) :: State_Chm  
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS
!
! !REVISION HISTORY: 
!  16 Apr 2015 - E. Lundgren - Initial version
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !=================================================================
    ! Convert_DryVV_to_MoistVV begins here!

         !==============================================================
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
         !  STT(I,J,L,N) [v/v moist]
         !
         !    = STT(I,J,L,N) [v/v dry] * State_Met%PMID_DRY(I,J,L) &   
         !                              / State_Met%PMID(I,J,L)
         !                   
         !==============================================================

    !=================================================================

      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( I, J, L, N ) 
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
        State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)     &
                                       * State_Met%PMID_DRY(I,J,L)  &
                                       / State_Met%PMID(I,J,L)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ! Return to calling program
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
  USE GIGC_ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?

    ! Number of tracers
    INTEGER,        INTENT(IN)    :: N_TRACERS 

    ! Object containing meteorological state variables
    TYPE(MetState), INTENT(IN)    :: State_Met 

! !INPUT/OUTPUT PARAMETERS:
!
    ! Object containing tracer concentration [v/v]
    TYPE(ChmState), INTENT(INOUT) :: State_Chm  
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS
!
! !REVISION HISTORY: 
!  16 Apr 2015 - E. Lundgren - Initial version
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !=================================================================
    ! Convert_MoistVV_to_DryVV begins here!

         !==============================================================
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
         !  STT(I,J,L,N) [v/v dry]
         !
         !    = STT(I,J,L,N) [v/v moist] * State_Met%PMID(I,J,L)
         !                               / State_Met%PMID_DRY(I,J,L)
         !                   
         !==============================================================

    !=================================================================

      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( I, J, L, N ) 
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
        State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N) &
                                       * State_Met%PMID(I,J,L)  &
                                       / State_Met%PMID_DRY(I,J,L)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ! Return to calling program
    END SUBROUTINE Convert_MoistVV_to_DryVV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_mmr_to_vmr
!
! !DESCRIPTION: Subroutine Convert\_MMR\_to\_VR converts the units of 
!  tracer concentrations from mass mixing ratio [kg tracer / kg dry
!  air] to volume mixing ratio [vol tracer /vol dry air]. 
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE Convert_MMR_to_VMR( N_TRACERS, TCVV, STT ) 
!
! !INPUT PARAMETERS: 
!
    ! Number of tracers
    INTEGER, INTENT(IN)      :: N_TRACERS 

    ! Array containing ratio of dry air to tracer molecular weights
    REAL(fp),  INTENT(IN)    :: TCVV(N_TRACERS)

! !OUTPUT PARAMETERS:
!
    ! Array containing tracer volume mixing ratio [vol/vol]
    REAL(fp),  INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,N_TRACERS)
!
! !REMARKS (edit this)
!  The volume ratio is the same as the molar ratio [mol/mol] under the 
!  same temperature and pressure conditions.  
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !=================================================================
    ! Convert_MMR_to_VMR begins here!
    !=================================================================

         !==============================================================
         !
         !  The conversion is as follows:
         !
         !   kg tracer(N)    g dry air      mol tracer(N)    
         !   -----------  * -----------  *  -------------  
         !   kg dry air     mol dry air      g tracer(N)          
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
         !  STT(I,J,L,N) [vol/vol]
         !
         !    = STT(I,J,L,N) [kg/kg] * TCVV(N)
         !                   
         !==============================================================
 
      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( I, J, L, N ) 
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
        STT(I,J,L,N) = STT(I,J,L,N) * TCVV(N)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ! Return to calling program
    END SUBROUTINE Convert_MMR_to_VMR
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_vmr_to_mmr
!
! !DESCRIPTION: Subroutine Convert\_VMR\_to\_MMR converts the units of 
!  tracer concentrations from volume mixing ratio [vol tracer /vol dry air] 
!  to mass mixing ratio [kg tracer / kg dry air]. 
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE Convert_VMR_to_MMR( N_TRACERS, TCVV, STT ) 
!
! !INPUT PARAMETERS: 
!
    ! Number of tracers
    INTEGER, INTENT(IN)      :: N_TRACERS 

    ! Array containing ratio of air to tracer molecular weights
    REAL(fp),  INTENT(IN)    :: TCVV(N_TRACERS)

! !OUTPUT PARAMETERS:
!
    ! Array containing tracer concentration [kg/kg]
    REAL(fp),  INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,N_TRACERS)
!
! !REMARKS (edit this)
!  The volume ratio is the same as the molar ratio [mol/mol] under the 
!  same temperature and pressure conditions.  
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !=================================================================
    ! Convert_VMR_to_MMR begins here!
    !=================================================================

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
         !  STT(I,J,L,N) [vol/vol]
         !
         !    = STT(I,J,L,N) [kg/kg] / TCVV(N)
         !                   
         !==============================================================

      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( I, J, L, N ) 
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
        STT(I,J,L,N) = STT(I,J,L,N) / TCVV(N)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ! Return to calling program
    END SUBROUTINE Convert_VMR_to_MMR
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_mmr_to_mnd
!
! !DESCRIPTION: Subroutine Convert\_MMR\_to\_MND converts the units of 
!  tracer concentrations from mass mixing ratio (MMR) [kg/kg] to 
!  molecular number density (MND) [molecules/cm3].  
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE Convert_MMR_to_MND( N_TRACERS, AIRDEN, XNUMOL, STT ) 
!
! !INPUT PARAMETERS: 
!
    ! Number of tracers
    INTEGER, INTENT(IN)      :: N_TRACERS 

    ! Array containing grid box air density [kg/m3]
    REAL(fp),  INTENT(IN)    :: AIRDEN(IIPAR,JJPAR,LLPAR)

    ! Array containing molecules tracer / kg tracer
    REAL(fp),  INTENT(IN)    :: XNUMOL(N_TRACERS)
!
! !OUTPUT PARAMETERS:
!
    ! Array containing tracer concentration [molecules/cm3]
    REAL(fp),  INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,N_TRACERS)
!
! !REMARKS
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !=================================================================
    ! Convert_MMR_to_MND begins here!
    !=================================================================

         !==============================================================
         !
         !  The conversion is as follows:
         !
         !   kg tracer(N)    kg air          m3        molecules tracer(N)     
         !   -----------  * --------  *  ---------  *  ------------------      
         !     kg air          m3         1E6 cm3        kg tracer(N)  
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
         !  STT(I,J,L,N) [molecules/cm3]
         !
         !    = STT(I,J,L,N) [kg/kg] * XNUMOL(N) * AIRDEN(I,J,L) * 1e-6
         !                   
         !==============================================================

      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( I, J, L, N ) 
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
        STT(I,J,L,N) = STT(I,J,L,N) * XNUMOL(N) * AIRDEN(I,J,L) * 1E-6_fp  
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ! Return to calling program
    END SUBROUTINE Convert_MMR_to_MND
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_mnd_to_mmr
!
! !DESCRIPTION: Subroutine Convert\_MND\_to\_MMR converts the units of 
!  tracer concentrations from molecular number density (MND)
!  [molecules/cm3] to mass mixing ratio (MMR) [kg/kg].  
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE Convert_MND_to_MMR( N_TRACERS, AIRDEN, XNUMOL, STT ) 
!
! !INPUT PARAMETERS: 
!
    ! Number of tracers
    INTEGER, INTENT(IN)      :: N_TRACERS 

    ! Array containing grid box dry air density [kg/m3]
    REAL(fp),  INTENT(IN)    :: AIRDEN(IIPAR,JJPAR,LLPAR)

    ! Array containing molecules tracer / kg tracer
    REAL(fp),  INTENT(IN)    :: XNUMOL(N_TRACERS)
!
! !OUTPUT PARAMETERS:
!
      ! Array containing tracer concentration [molecules/cm3]
    REAL(fp),  INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,N_TRACERS)
!
! !REMARKS
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !=================================================================
    ! Convert_MND_to_MMR begins here!
    !=================================================================

         !==============================================================
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
         !  XNUMOL(N)       = molecules tracer / kg tracer
         !  AIRDEN(I,J,L)  = grid box dry air density [kg/m3]
         !     
         ! the conversion is:
         ! 
         !  STT(I,J,L,N) [kg/kg]
         !
         !   = STT(I,J,L,N) [molecules/cm3] * 1E+6 
         !                          / ( XNUMOL(N) * AIRDEN(I,J,L) ) 
         !                   
         !==============================================================

      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( I, J, L, N ) 
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
        STT(I,J,L,N) = STT(I,J,L,N) * 1E+6_fp / ( XNUMOL(N) * AIRDEN(I,J,L) ) 
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ! Return to calling program
    END SUBROUTINE Convert_MND_to_MMR
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_mmr_to_kg
!
! !DESCRIPTION: Subroutine Convert\_MMR\_to\_KG converts the units of 
!  tracer concentrations from mass mixing ratio (MMR) 
!  [kg tracer/kg dry air] to tracer mass per grid box [kg]. 
!
!  NOTE: This will go away in the future.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_MMR_to_KG( N_TRACERS, AIRMASS, STT ) 
!
! !INPUT PARAMETERS: 
!
    ! Number of tracers
    INTEGER, INTENT(IN)      :: N_TRACERS 

    ! Array containing grid box dry air mass [kg]
    REAL(fp),  INTENT(IN)    :: AIRMASS(IIPAR,JJPAR,LLPAR)
!
! !OUTPUT PARAMETERS:
!
    ! Array containing tracer concentration [kg]
    REAL(fp),  INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,N_TRACERS)
!
! !REMARKS
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !=================================================================
    ! Convert_MMR_to_KG begins here!
    !=================================================================

         !==============================================================
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
         !  STT(I,J,L,N) [kg]
         !
         !    = STT(I,J,L,N) [kg/kg] * AIRMASS(I,J,L)
         !                   
         !==============================================================

      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( I, J, L, N ) 
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
        STT(I,J,L,N) = STT(I,J,L,N) * AIRMASS(I,J,L)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ! Return to calling program
    END SUBROUTINE Convert_MMR_to_KG
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kg_to_mmr
!
! !DESCRIPTION: Subroutine Convert\_KG\_to\_MMR converts the units of 
!  tracer concentrations from tracer mass per grid box [kg] to mass 
!  mixing ratio (MMR) [kg tracer/kg air]. 
!  
!  This will go away in the future with grid-independence.
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE Convert_KG_to_MMR( N_TRACERS, AIRMASS, STT ) 
!
! !INPUT PARAMETERS: 
!
    ! Number of tracers
    INTEGER, INTENT(IN)      :: N_TRACERS 

    ! Array containing grid box dry air mass [kg]
    REAL(fp),  INTENT(IN)    :: AIRMASS(IIPAR,JJPAR,LLPAR)

!
! !OUTPUT PARAMETERS:
!
    ! Array containing tracer concentration [kg/kg]
    REAL(fp),  INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,N_TRACERS)
!
! !REMARKS
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N

    !=================================================================
    ! Convert_KG_to_MMR begins here!
    !=================================================================

         !==============================================================
         !
         !  The conversion is as follows:
         !
         !                         1          
         !   kg tracer(N)  * --------------      
         !                      kg dry air              
         !
         !   = kg tracer(N) / dry air mass
         !   
         !   = mass mixing ratio
         !
         ! Therefore, with:
         !
         !  AIRMASS(I,J,L)    = grid box air mass [kg]
         !     
         ! the conversion is:
         ! 
         !  STT(I,J,L,N) [kg/kg]
         !
         !    = STT(I,J,L,N) [kg] / AIRMASS(I,J,L) 
         !                   
         !==============================================================

      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( I, J, L, N ) 
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
        STT(I,J,L,N) = STT(I,J,L,N) / AIRMASS(I,J,L)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ! Return to calling program
    END SUBROUTINE Convert_KG_to_MMR
!EOC

END MODULE UnitConv_Mod
