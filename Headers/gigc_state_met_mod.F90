#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_state_met_mod
!
! !DESCRIPTION: Module GIGC\_STATE\_MET\_MOD contains the derived type
!  used to define the Meteorology State object for the Grid-Independent 
!  GEOS-Chem implementation (abbreviated "GIGC").
!\\
!\\
!  This module also contains the routines that allocate and deallocate memory 
!  to the Meteorology State object.  The chemistry state object is not defined
!  in this module.  It must be be declared as variable in the top-level 
!  driver routine, and then passed to lower-level routines as an argument.
!\\
!\\
! !INTERFACE: 
!
MODULE GIGC_State_Met_Mod
!
! USES:
!
  IMPLICIT NONE
# include "define.h"
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_GIGC_State_Met
  PUBLIC :: Cleanup_GIGC_State_Met
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for Meteorology State
  !=========================================================================
  TYPE, PUBLIC :: MetState

     ! Surface fields
     REAL*8, POINTER :: ALBD    (:,:  )   ! Visible surface albedo [unitless]
     REAL*8, POINTER :: CLDFRC  (:,:  )   ! Column cloud fraction [unitless]
     REAL*8, POINTER :: FRCLND  (:,:  )   ! Olson land fraction [unitless]
     REAL*8, POINTER :: GWETTOP (:,:  )   ! Top soil moisture [unitless]
     REAL*8, POINTER :: HFLUX   (:,:  )   ! Sensible heat flux [W/m2]
     REAL*8, POINTER :: LWI     (:,:  )   ! Land/water indices [unitless]
     REAL*8, POINTER :: PARDR   (:,:  )   ! Direct  photsyn active rad [W/m2]
     REAL*8, POINTER :: PARDF   (:,:  )   ! Diffuse photsyn active rad [W/m2]
     REAL*8, POINTER :: PBLH    (:,:  )   ! PBL height [m]
     REAL*8, POINTER :: PRECCON (:,:  )   ! Conv  precip @ ground [kg/m2/s]
     REAL*8, POINTER :: PRECTOT (:,:  )   ! Total precip @ ground [kg/m2/s]
     REAL*8, POINTER :: RADSWG  (:,:  )   ! Solar radiation @ ground [W/m2]
     REAL*8, POINTER :: SST     (:,:  )   ! Sea surface temperature [K]
     REAL*8, POINTER :: SUNCOS  (:    )   ! Cosine of solar zenith angle
     REAL*8, POINTER :: TO3     (:,:  )   ! Total overhead O3 column [DU]
     REAL*8, POINTER :: TROPP   (:,:  )   ! Tropopause pressure [hPa]
     REAL*8, POINTER :: TS      (:,:  )   ! Surface temperature [K]
     REAL*8, POINTER :: U10M    (:,:  )   ! E/W wind speed @ 10m height [m/s]
     REAL*8, POINTER :: USTAR   (:,:  )   ! Friction velocity [m/s]
     REAL*8, POINTER :: UVALBEDO(:,:  )   ! UV surface albedo [unitless]
     REAL*8, POINTER :: V10M    (:,:  )   ! N/S wind speed @ 10m height [m/s]
     REAL*8, POINTER :: Z0      (:,:  )   ! Surface roughness height [m]

     ! 3-D Fields
     REAL*8, POINTER :: AD      (:,:,:)   ! Air mass [kg]
     REAL*8, POINTER :: AIRDENS (:,:,:)   ! Air density [kg/m3]
     REAL*8, POINTER :: AIRVOL  (:,:,:)   ! Grid box volume [m3]
     REAL*8, POINTER :: AREA_M2 (:,:,:)   ! Grid box surface area [cm2]
     REAL*8, POINTER :: BXHEIGHT(:,:,:)   ! Grid box height [m]
     REAL*8, POINTER :: CLDF    (:,:,:)   ! 3-D cloud fraction [unitless]
     REAL*8, POINTER :: CMFMC   (:,:,:)   ! Cloud mass flux [kg/m2/s]
     REAL*8, POINTER :: DQIDTMST(:,:,:)   ! Ice tendency, mst proc [kg/kg/s]
     REAL*8, POINTER :: DQLDTMST(:,:,:)   ! H2O tendency, mst proc [kg/kg/s]
     REAL*8, POINTER :: DQVDTMST(:,:,:)   ! Vapor tendency, mst proc [kg/kg/s]
     REAL*8, POINTER :: DTRAIN  (:,:,:)   ! Detrainment flux [kg/m2/s]
     REAL*8, POINTER :: MOISTQ  (:,:,:)   ! Tendency in sp. humidity [kg/kg/s]
     REAL*8, POINTER :: OPTD    (:,:,:)   ! Visible optical depth [unitless]
     REAL*8, POINTER :: PEDGE   (:,:,:)   ! Pressure @ level edges [Pa]
     REAL*8, POINTER :: PMID    (:,:,:)   ! Pressure @ level centers [Pa]
     REAL*8, POINTER :: DELP    (:,:,:)   ! Delta-P extent  of a grid box [mb]
     REAL*8, POINTER :: RH      (:,:,:)   ! Relative humidity [unitless]
     REAL*8, POINTER :: SPHU    (:,:,:)   ! Specific humidity [kg/kg]
     REAL*8, POINTER :: T       (:,:,:)   ! Temperature [K]
     REAL*8, POINTER :: TAUCLI  (:,:,:)   ! Opt depth of ice clouds [unitless]
     REAL*8, POINTER :: TAUCLW  (:,:,:)   ! Opt depth of H2O clouds [unitless]

  END TYPE MetState
!
! !REVISION HISTORY: 
!  19 Oct 2012 - R. Yantosca - Initial version, split off from gc_type_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gigc_state_met
!
! !DESCRIPTION: Subroutine INIT\_GIGC\_STATE\_MET allocates all fields of 
!  the Grid-Indpendent GEOS-Chem (aka "GIGC") Meteorology State object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_GIGC_State_Met( am_I_Root, IM, JM, LM, State_Met, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod                         ! Error codes
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    INTEGER,        INTENT(IN)    :: IM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: JM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: LM          ! # longitudes on this PET
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Obj for meteorology state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!  For consistency, maybe this should be moved to a different module.
!
! !REVISION HISTORY: 
!  19 Oct 2012 - R. Yantosca - Initial version, based on gc_environment_mod.F90
!  19 Oct 2012 - R. Yantosca - Now pass all dimensions as arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Assume success
    RC = GIGC_SUCCESS

    !=======================================================================
    ! Allocate 2-D Fields
    !=======================================================================
    ALLOCATE( State_Met%ALBD    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%CLDFRC  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%FRCLND  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%GWETTOP ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%HFLUX   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%LWI     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%PARDR   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%PARDF   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%PBLH    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%PRECCON ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%PRECTOT ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%RADSWG  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%SST     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%SUNCOS  ( IM* JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%TO3     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%TROPP   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%TS      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%U10M    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%USTAR   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%UVALBEDO( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%V10M    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Met%Z0      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    !=======================================================================
    ! Allocate 3-D Arrays
    !=======================================================================
    ALLOCATE( State_Met%AD      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%AIRDENS ( LM, IM, JM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%AIRVOL  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%AREA_M2 ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%BXHEIGHT( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%CLDF    ( LM, IM, JM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%CMFMC   ( IM, JM, LM+1 ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%DQIDTMST( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%DQLDTMST( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%DQVDTMST( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%DTRAIN  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%MOISTQ  ( LM, IM, JM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%OPTD    ( LM, IM, JM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%PEDGE   ( IM, JM, LM+1 ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%PMID    ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%DELP    ( LM, IM, JM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%RH      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%SPHU    ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%T       ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%TAUCLI  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
                                               
    ALLOCATE( State_Met%TAUCLW  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    
  END SUBROUTINE Init_GIGC_State_Met
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_gigc_state_met
!
! !DESCRIPTION: Subroutine CLEANUP\_GIGC\_STATE\_MET allocates all fields 
!  of the Grid-Independent GEOS-Chem (aka "GIGC") Meteorology State object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_GIGC_State_Met( am_I_Root, State_Met, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod                         ! Error codes
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Obj for meteorology state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!  For consistency, maybe this should be moved to a different module.
!
! !REVISION HISTORY: 
!  19 Oct 2012 - R. Yantosca - Initial version, based on gc_environment_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
     ! Surface fields
     IF ( ASSOCIATED( State_Met%ALBD     )  ) DEALLOCATE( State_Met%ALBD     ) 
     IF ( ASSOCIATED( State_Met%CLDFRC   )  ) DEALLOCATE( State_Met%CLDFRC   ) 
     IF ( ASSOCIATED( State_Met%FRCLND   )  ) DEALLOCATE( State_Met%FRCLND   )
     IF ( ASSOCIATED( State_Met%GWETTOP  )  ) DEALLOCATE( State_Met%GWETTOP  )
     IF ( ASSOCIATED( State_Met%HFLUX    )  ) DEALLOCATE( State_Met%HFLUX    )
     IF ( ASSOCIATED( State_Met%LWI      )  ) DEALLOCATE( State_Met%LWI      )
     IF ( ASSOCIATED( State_Met%PARDR    )  ) DEALLOCATE( State_Met%PARDR    )
     IF ( ASSOCIATED( State_Met%PARDF    )  ) DEALLOCATE( State_Met%PARDF    )
     IF ( ASSOCIATED( State_Met%PBLH     )  ) DEALLOCATE( State_Met%PBLH     )
     IF ( ASSOCIATED( State_Met%PRECCON  )  ) DEALLOCATE( State_Met%PRECCON  )
     IF ( ASSOCIATED( State_Met%PRECTOT  )  ) DEALLOCATE( State_Met%PRECTOT  )
     IF ( ASSOCIATED( State_Met%RADSWG   )  ) DEALLOCATE( State_Met%RADSWG   )
     IF ( ASSOCIATED( State_Met%SST      )  ) DEALLOCATE( State_Met%SST      )
     IF ( ASSOCIATED( State_Met%SUNCOS   )  ) DEALLOCATE( State_Met%SUNCOS   )
     IF ( ASSOCIATED( State_Met%TO3      )  ) DEALLOCATE( State_Met%TO3      )
     IF ( ASSOCIATED( State_Met%TROPP    )  ) DEALLOCATE( State_Met%TROPP    )
     IF ( ASSOCIATED( State_Met%TS       )  ) DEALLOCATE( State_Met%TS       )
     IF ( ASSOCIATED( State_Met%U10M     )  ) DEALLOCATE( State_Met%U10M     )
     IF ( ASSOCIATED( State_Met%USTAR    )  ) DEALLOCATE( State_Met%USTAR    )
     IF ( ASSOCIATED( State_Met%UVALBEDO )  ) DEALLOCATE( State_Met%V10M     )
     IF ( ASSOCIATED( State_Met%V10M     )  ) DEALLOCATE( State_Met%UVALBEDO )
     IF ( ASSOCIATED( State_Met%Z0       )  ) DEALLOCATE( State_Met%Z0       )
     IF ( ASSOCIATED( State_Met%AD       )  ) DEALLOCATE( State_Met%AD       )
     IF ( ASSOCIATED( State_Met%AIRDENS  )  ) DEALLOCATE( State_Met%AIRDENS  )
     IF ( ASSOCIATED( State_Met%AIRVOL   )  ) DEALLOCATE( State_Met%AIRVOL   )
     IF ( ASSOCIATED( State_Met%AREA_M2  )  ) DEALLOCATE( State_Met%AREA_M2  )
     IF ( ASSOCIATED( State_Met%BXHEIGHT )  ) DEALLOCATE( State_Met%BXHEIGHT )
     IF ( ASSOCIATED( State_Met%CLDF     )  ) DEALLOCATE( State_Met%CLDF     )
     IF ( ASSOCIATED( State_Met%CMFMC    )  ) DEALLOCATE( State_Met%CMFMC    )
     IF ( ASSOCIATED( State_Met%DQIDTMST )  ) DEALLOCATE( State_Met%DQIDTMST )
     IF ( ASSOCIATED( State_Met%DQLDTMST )  ) DEALLOCATE( State_Met%DQLDTMST )
     IF ( ASSOCIATED( State_Met%DQVDTMST )  ) DEALLOCATE( State_Met%DQVDTMST )
     IF ( ASSOCIATED( State_Met%DTRAIN   )  ) DEALLOCATE( State_Met%DTRAIN   )
     IF ( ASSOCIATED( State_Met%MOISTQ   )  ) DEALLOCATE( State_Met%MOISTQ   )
     IF ( ASSOCIATED( State_Met%OPTD     )  ) DEALLOCATE( State_Met%OPTD     )
     IF ( ASSOCIATED( State_Met%PEDGE    )  ) DEALLOCATE( State_Met%PEDGE    )
     IF ( ASSOCIATED( State_Met%PMID     )  ) DEALLOCATE( State_Met%PMID     )
     IF ( ASSOCIATED( State_Met%DELP     )  ) DEALLOCATE( State_Met%DELP     )
     IF ( ASSOCIATED( State_Met%RH       )  ) DEALLOCATE( State_Met%RH       )
     IF ( ASSOCIATED( State_Met%SPHU     )  ) DEALLOCATE( State_Met%SPHU     )
     IF ( ASSOCIATED( State_Met%T        )  ) DEALLOCATE( State_Met%T        )
     IF ( ASSOCIATED( State_Met%TAUCLI   )  ) DEALLOCATE( State_Met%TAUCLI   )
     IF ( ASSOCIATED( State_Met%TAUCLW   )  ) DEALLOCATE( State_Met%TAUCLW   ) 

   END SUBROUTINE Cleanup_GIGC_State_Met
!EOC
END MODULE GIGC_State_Met_Mod
#endif
