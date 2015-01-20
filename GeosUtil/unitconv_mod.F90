!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: unitconv_mod.F90
!
! !DESCRIPTION: Module UNITCONV\_MOD contains routines which are used to 
!  convert the units of tracer concentrations between mass mixing ratio 
!  [kg/kg], mass per grid box [kg], volume ratio [vol/vol], and molecular 
!  number density [molecules/cm3].
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
  USE ERROR_MOD             ! Error-handling routines
  USE PRECISION_MOD         ! GEOS-Chem Flexible Precision (fp)
                    
  IMPLICIT NONE
  PRIVATE
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Convert_MMR_to_VR     ! kg / kg       -> vol / vol
  PUBLIC  :: Convert_VR_to_MMR     ! vol / vol     -> kg / kg
  PUBLIC  :: Convert_MMR_to_KG     ! kg / kg       -> kg / grid box
  PUBLIC  :: Convert_KG_to_MMR     ! kg / grid box -> kg / kg
  PUBLIC  :: Convert_MMR_to_MND    ! kg / kg       -> molec / cm3
  PUBLIC  :: Convert_MND_to_MMR    ! molec / cm3   -> kg / kg
!
! !REMARKS:
!  The routines in this module are used to convert the units of tracer 
!  concentrations for use in other GEOS-Chem routines. Tracer concentrations 
!  are stored in units of mass mixing ratio [kg/kg] but most subroutines 
!  require tracer units in mass per grid box [kg] or in volume ratio 
!  [vol/vol] (same as molar ratio [mol/mol]). Some routines also require 
!  conversion to and from molecular number density [molecules/cm3]. 
!
!  Species concentrations are stored in units of molecular number density
!  and the unit conversion routines within this module may therefore be 
!  used to convert species concentration units, if needed.
!
! !REVISION HISTORY:
!  06 Jan 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_mmr_to_vr
!
! !DESCRIPTION: Subroutine Convert\_MMR\_to\_VR converts the units of 
!  tracer concentrations (STT) from mass mixing ratio (MMR) [kg/kg] to 
!  volume ratio (VR) [vol/vol] (same as molar ratio [mol/mol]). 
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE Convert_MMR_to_VR( N_TRACERS, TCVV, STT ) 
!
! !INPUT PARAMETERS: 
!
    ! Number of tracers
    INTEGER, INTENT(IN)      :: N_TRACERS 

    ! Array containing ratio of air to tracer molecular weights
    REAL(fp),  INTENT(IN)    :: TCVV(N_TRACERS)

! !OUTPUT PARAMETERS:
!
    ! Array containing tracer concentration [vol/vol]
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
    ! Convert_MMR_to_VR begins here!
    !=================================================================

         !==============================================================
         !
         !  The conversion is as follows:
         !
         !   kg tracer(N)    g air      mol tracer(N)    
         !   -----------  * -------  *  -------------  
         !     kg air       mol air      g tracer(N)          
         !
         !   = mass mixing ratio * ratio of air to tracer molecular weights  
         !   
         !   = molar ratio
         !
         ! Therefore, with:
         !
         !  TCVV(N) = air molecular wt / tracer molecular wt 
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
    END SUBROUTINE Convert_MMR_to_VR
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_vr_to_mmr
!
! !DESCRIPTION: Subroutine Convert\_VR\_to\_MMR converts the units of 
!  tracer concentrations (STT) from volume ratio (VR) [vol/vol] (same 
!  as molar ratio [mol/mol]) to mass mixing ratio (MMR) [kg/kg]. 
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE Convert_VR_to_MMR( N_TRACERS, TCVV, STT ) 
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
    ! Convert_VR_to_MMR begins here!
    !=================================================================

         !==============================================================
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
         !  TCVV(N) = air molecular wt / tracer molecular wt 
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
    END SUBROUTINE Convert_VR_to_MMR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_mmr_to_kg
!
! !DESCRIPTION: Subroutine Convert\_MMR\_to\_KG converts the units of 
!  tracer concentrations (STT) from mass mixing ratio (MMR) 
!  [kg tracer/kg air] to tracer mass per grid box [kg]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_MMR_to_KG( N_TRACERS, MAIRDEN, AIRVOL, STT ) 
!
! !INPUT PARAMETERS: 
!
    ! Number of tracers
    INTEGER, INTENT(IN)      :: N_TRACERS 

    ! Array containing grid box moist air density [kg/m3]
    REAL(fp),  INTENT(IN)    :: MAIRDEN(IIPAR,JJPAR,LLPAR)

    ! Array containing grid box volume [m3]
    REAL(fp),  INTENT(IN)    :: AIRVOL(IIPAR,JJPAR,LLPAR)
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
         !   kg tracer(N)   kg moist air          
         !   -----------  * ------------  *  m3       
         !     kg air            m3              
         !
         !   = mass mixing ratio * moist air density * grid box volume  
         !   
         !   = kg tracer(N)
         !
         ! Therefore, with:
         !
         !  MAIRDEN(I,J,L)   = grid box moist air density [kg/m3]
         !  AIRVOL(I,J,L)     = grid box volume [m3]
         !     
         ! the conversion is:
         ! 
         !  STT(I,J,L,N) [kg]
         !
         !    = STT(I,J,L,N) [kg/kg] * MAIRDEN(I,J,L) * AIRVOL(I,J,L)
         !                   
         !==============================================================

      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( I, J, L, N ) 
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
        STT(I,J,L,N) = STT(I,J,L,N) * AIRVOL(I,J,L) * MAIRDEN(I,J,L)
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
!  tracer concentrations (STT) from tracer mass per grid box [kg] to mass 
!  mixing ratio (MMR) [kg tracer/kg air]. 
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE Convert_KG_to_MMR( N_TRACERS, MAIRDEN, AIRVOL, STT ) 
!
! !INPUT PARAMETERS: 
!
    ! Number of tracers
    INTEGER, INTENT(IN)      :: N_TRACERS 

    ! Array containing grid box moist air density [kg/m3]
    REAL(fp),  INTENT(IN)    :: MAIRDEN(IIPAR,JJPAR,LLPAR)

    ! Array containing grid box volume [m3]
    REAL(fp),  INTENT(IN)    :: AIRVOL(IIPAR,JJPAR,LLPAR)
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
         !                          m3          
         !   kg tracer(N)  * --------------  /  m3       
         !                    kg moist air              
         !
         !   = kg tracer(N) / moist air density / grid box volume  
         !   
         !   = mass mixing ratio
         !
         ! Therefore, with:
         !
         !  MAIRDEN(I,J,L)   = grid box moist air density [kg/m3]
         !  AIRVOL(I,J,L)     = grid box volume [m3]
         !     
         ! the conversion is:
         ! 
         !  STT(I,J,L,N) [kg/kg]
         !
         !    = STT(I,J,L,N) [kg] / ( AIRVOL(I,J,L) * MAIRDEN(I,J,L) )
         !                   
         !==============================================================

      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( I, J, L, N ) 
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
        STT(I,J,L,N) = STT(I,J,L,N) / ( AIRVOL(I,J,L) * MAIRDEN(I,J,L) ) 
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ! Return to calling program
    END SUBROUTINE Convert_KG_to_MMR
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_mmr_to_mnd
!
! !DESCRIPTION: Subroutine Convert\_MMR\_to\_MND converts the units of 
!  tracer concentrations (STT) from mass mixing ratio (MMR) [kg/kg] to 
!  molecular number density (MND) [molecules/cm3].  
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE Convert_MMR_to_MND( N_TRACERS, MAIRDEN, XNUMOL, STT ) 
!
! !INPUT PARAMETERS: 
!
    ! Number of tracers
    INTEGER, INTENT(IN)      :: N_TRACERS 

    ! Array containing grid box moist air density [kg/m3]
    REAL(fp),  INTENT(IN)    :: MAIRDEN(IIPAR,JJPAR,LLPAR)

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
         !   = mass mixing ratio * moist air density * molecules / kg tracer
         !   
         !   = molecules per cm3
         !
         ! Therefore, with:
         !
         !  XNUMOL(N)        = molecules tracer / kg tracer
         !  MAIRDEN(I,J,L)   = grid box moist air density [kg/m3]
         !     
         ! the conversion is:
         ! 
         !  STT(I,J,L,N) [molecules/cm3]
         !
         !    = STT(I,J,L,N) [kg/kg] * XNUMOL(N) * MAIRDEN(I,J,L) * 1e-6
         !                   
         !==============================================================

      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( I, J, L, N ) 
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
        STT(I,J,L,N) = STT(I,J,L,N) * XNUMOL(N) * MAIRDEN(I,J,L) * 1E-6_fp  
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
!  tracer concentrations (STT) from molecular number density (MND)
!  [molecules/cm3] to mass mixing ratio (MMR) [kg/kg].  
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE Convert_MND_to_MMR( N_TRACERS, MAIRDEN, XNUMOL, STT ) 
!
! !INPUT PARAMETERS: 
!
    ! Number of tracers
    INTEGER, INTENT(IN)      :: N_TRACERS 

    ! Array containing grid box moist air density [kg/m3]
    REAL(fp),  INTENT(IN)    :: MAIRDEN(IIPAR,JJPAR,LLPAR)

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
         !  MAIRDEN(I,J,L)  = grid box moist air density [kg/m3]
         !     
         ! the conversion is:
         ! 
         !  STT(I,J,L,N) [kg/kg]
         !
         !   = STT(I,J,L,N) [molecules/cm3] * 1E+6 
         !                          / ( XNUMOL(N) * MAIRDEN(I,J,L) ) 
         !                   
         !==============================================================

      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( I, J, L, N ) 
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
        STT(I,J,L,N) = STT(I,J,L,N) * 1E+6_fp / ( XNUMOL(N) * MAIRDEN(I,J,L) ) 
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ! Return to calling program
    END SUBROUTINE Convert_MND_to_MMR
!EOC
END MODULE UnitConv_Mod
