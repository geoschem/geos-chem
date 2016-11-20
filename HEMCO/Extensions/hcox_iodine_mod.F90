!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_Iodine_mod.F90
!
! !DESCRIPTION: Module HCOX\_Iodine\_Mod contains routines to calculate
! oceanic iodine emissions (HOI and I2), following carpenter et al. (2014). 
! The emission is parameterised herein using online feilds for O3, 10 metre
! wind speed, and ocean surface iodide concentration (parameterised from 
! STT following Chance et al (2014)).
!
!\\
!\\ 
! This is a HEMCO extension module that uses many of the HEMCO core
! utilities.
!\\
!\\
! !INTERFACE: 
!
MODULE HCOX_Iodine_Mod
!
! !USES:
! 
  USE HCO_Error_Mod
  USE HCO_Diagn_Mod
  USE HCO_State_Mod,  ONLY : HCO_State
  USE HCOX_State_Mod, ONLY : Ext_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOX_Iodine_Init
  PUBLIC :: HCOX_Iodine_Run
  PUBLIC :: HCOX_Iodine_Final
!
! !PRIVATE MEMBER FUNCTIONS:
!
! N/A
!
! !REVISION HISTORY:
!  15 Mar 2013 - T. Sherwen - Initial implementation (v9-3-01)
!  15 Jul 2015 - T. Sherwen - Now a HEMCO extension module 
!EOP
!------------------------------------------------------------------------------
!
! !PRIVATE TYPES:
!
  ! Tracer IDs 
  INTEGER             :: ExtNr
  INTEGER             :: IDTI2            ! I2 model species ID
  INTEGER             :: IDTHOI           ! HOI model species ID
  LOGICAL             :: CalcI2           ! Calculate I2 oceanic emissions?
  LOGICAL             :: CalcHOI          ! Calculate HOI oceanic emissions?

!
! !DEFINED PARAMETERS:
!
   ! Molecular weight of I2 [kg/mol] 
   REAL*8,  PARAMETER   :: MWT_I2 = 2.54d-1
   ! Molecular weight of HOI [kg/mol]
   REAL*8,  PARAMETER   :: MWT_HOI = 1.44d-1

CONTAINS


!EOC
!-------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Iodine_Run 
!
! !DESCRIPTION: Subroutine HcoX\_Iodine\_Run is the driver run routine to 
! calculate ocean inorganic iodine emissions in HEMCO.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Iodine_Run( am_I_Root, ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_FluxArr_Mod,      ONLY : HCO_EmisAdd
    USE HCO_GeoTools_Mod,     ONLY : HCO_LANDTYPE
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   ) :: am_I_Root  ! root CPU?
    TYPE(HCO_State), POINTER       :: HcoState   ! Output obj
    TYPE(Ext_State), POINTER       :: ExtState   ! Module options  
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC         ! Success or failure?
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Sherwen et al. 2015
!  (2 ) Carpenter et al. 2013
!  (3 ) Chance et al. 2014
!  (4 ) Macdonal et al. 2014
!
! !REVISION HISTORY: 
!  15 Mar 2013 - T. Sherwen - Initial implementation (v9-3-01)
!  15 Jul 2015 - T. Sherwen - Now a HEMCO extension module 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: I, J
    REAL*8                 :: EMIS_HOI
    REAL*8                 :: EMIS_I2, IODIDE, O3_CONC
    REAL*8                 :: SST
    REAL*8                 :: A_M2
    REAL*8                 :: W10M
    REAL(hp), TARGET       :: FLUXHOI (HcoState%NX,HcoState%NY)
    REAL(hp), TARGET       :: FLUXI2 (HcoState%NX,HcoState%NY)

    ! Error handling
    LOGICAL                :: ERR
    CHARACTER(LEN=255)     :: MSG

    !=================================================================
    ! HCOX_Iodine_Run begins here!
    !=================================================================

    ! Return if extension disabled 
    IF ( .NOT. ExtState%Inorg_Iodine ) RETURN

    ! Enter 
    CALL HCO_ENTER ( 'HCOX_Iodine_Run (hcox_iodine_mod.F90)', RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Exit status
    ERR = .FALSE.

    ! Initialize flux arrays/variables
    FLUXHOI  = 0.0_hp
    FLUXI2   = 0.0_hp
    EMIS_HOI = 0.0_hp
    EMIS_I2  = 0.0_hp
    
    ! Loop over surface boxes 
    DO J = 1, HcoState%NY 
    DO I = 1, HcoState%NX

       ! Advance to next grid box if box is not over water
       ! further check needed for ocean, but not available
       ! ( as parameterisation of iodide based on ocean data)
       IF ( HCO_LANDTYPE( ExtState%WLI%Arr%Val(I,J), &
                          ExtState%ALBD%Arr%Val(I,J) ) /= 0 ) CYCLE

       ! Grid box surface area on simulation grid [m2]
       A_M2 = HcoState%Grid%AREA_M2%Val( I, J )

       ! Wind speed at 10 m altitude [m/s]
       W10M = SQRT( ExtState%U10M%Arr%Val(I,J)**2 &
                  + ExtState%V10M%Arr%Val(I,J)**2 ) 

       ! limit W10M to a minimium of 5 m/s to avoid overestimation of fluxes                                                               
       ! from CARPENTER et al. (2013) (per. comm.)
       IF ( W10M .LE. 5d0  ) THEN
           W10M   =  5d0
       ENDIF

       ! Sea surface temperature in Celcius 
!       SST = ExtState%TSKIN%Arr%Val(I,J) - 273.15d0
       ! Sea surface temperature in Kelvin 
       SST = ExtState%TSKIN%Arr%Val(I,J) 

       ! Empirical SST scaling factor (jaegle 5/11/11)
!       SCALE = 0.329d0 + 0.0904d0*SST - &
!               0.00717d0*SST**2d0 + 0.000207d0*SST**3d0

!       ! SST dependence of iodide - Chance et al. 2014, in press
!       IODIDE = ( (0.225d0 * ( (SST)**2d0) )  + 19d0 )  / 1d9
!       ! SST dependence of iodide - Macdonald et al. 2014
       IODIDE = 1.46d6 * EXP( (-9134d0/SST) )

       ! Get O3 concentration at the surface ( in mol/mol )
       ! ExtState%O3 is in units of kg/kg dry air
       O3_CONC = ExtState%O3%Arr%Val(I,J,1)         &
               * HcoState%Phys%AIRMW / 48.0_dp &
               * 1.e9_dp

       ! Reset to using original Gong (2003) emissions (jaegle 6/30/11)
       !SCALE = 1.0d0

       ! Eventually apply wind scaling factor. 
!       SCALE = SCALE * WindScale

       ! If I2 & emitting, use parameterisation from
       ! Carpenter et al (2013) to give emissions in nmol m-2 d-1.
       ! Then convert this to kg/m2/s
       IF ( CalcI2 ) THEN
           EMIS_I2 = ( O3_CONC * (IODIDE**1.3d0) * &
               ( ( 1.74d9 - ( 6.54d8*LOG( W10M ) )   ) )/ &
                     24d0/60d0/60d0/1d9*MWT_I2 )
!
          ! If parametsation results in negative ( W10 too high )
          ! flux set to zero
          IF ( EMIS_I2 .LT. 0d0 ) THEN
             EMIS_I2 = 0d0
          ENDIF
!
       ENDIF
!                                                                                                                                                 
       IF ( CalcHOI ) THEN
       ! If HOI & emitting, use parameterisation from
       ! Carpenter et al (2013) to give emissions in nmol m-2 d-1.
       ! Then convert this to kg/m2/s

         EMIS_HOI =  O3_CONC * &
            ( ( 4.15d5 * ( SQRT(IODIDE)/ W10M ) ) - &
            ( 20.6 / W10M ) - ( 2.36d4  * SQRT(IODIDE) ) ) / &
                      24d0/60d0/60d0/1d9*MWT_HOI   

         ! If parametsation results in negative ( W10 too high )
         ! flux set to zero
         IF ( EMIS_HOI .LT. 0d0 ) THEN
                EMIS_HOI = 0d0
         ENDIF

       ENDIF

       ! Store HOI flux in tendency array in [kg/m2/s]
       IF ( CalcHOI ) THEN 

          ! kg --> kg/m2/s
          FLUXHOI(I,J) = EMIS_HOI
	  	       

       ENDIF

       ! store I2 flux in tendency array in [kg/m2/s]
       IF ( CalcI2 ) THEN 

          ! kg --> kg/m2/s
          FLUXI2(I,J) = EMIS_I2

       ENDIF

    ENDDO !I
    ENDDO !J

    ! Check exit status 
    IF ( ERR ) THEN
       RC = HCO_FAIL
       RETURN 
    ENDIF

    !=================================================================
    ! PASS TO HEMCO STATE AND UPDATE DIAGNOSTICS 
    !=================================================================

    ! HOI 
    IF ( CalcHOI ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXHOI, IDTHOI, & 
                         RC,        ExtNr=ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXHOI', RC )
          RETURN 
       ENDIF

    ENDIF

    ! I2
    IF ( CalcI2 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXI2, IDTI2, & 
                         RC,        ExtNr=ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXI2', RC )
          RETURN 
       ENDIF

    ENDIF
      
    ! Leave w/ success
    CALL HCO_LEAVE ( RC )

  END SUBROUTINE HCOX_Iodine_Run

!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Iodine_Init
!
! !DESCRIPTION: Subroutine HcoX\_Iodine\_Init initializes all
!  extension variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Iodine_Init( am_I_Root, HcoState, ExtName, ExtState, RC ) 
!
! !USES:
!
    USE HCO_State_Mod,          ONLY : HCO_GetHcoID
    USE HCO_STATE_MOD,          ONLY : HCO_GetExtHcoID
    USE HCO_ExtList_Mod,        ONLY : GetExtNr
    USE HCO_ExtList_Mod,        ONLY : GetExtOpt
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root   ! root CPU?
    TYPE(HCO_State),  POINTER        :: HcoState    ! HEMCO state object 
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName     ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState    ! Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC          ! Return status
!
! !REVISION HISTORY:
!  15 Mar 2013 - T. Sherwen - Initial implementation (v9-3-01)
!  15 Jul 2015 - T. Sherwen - Now a HEMCO extension module adapted from hcox_seasalt_mod 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                        :: N, R, AS
    CHARACTER(LEN=255)             :: MSG
    INTEGER                        :: nSpc, minLen
    LOGICAL                        :: FOUND
    INTEGER, ALLOCATABLE           :: HcoIDs(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)

    !=================================================================
    ! HCOX_Iodine_Init begins here!
    !=================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN
 
    ! Enter 
    CALL HCO_ENTER ( 'HCOX_iodine_Init (hcox_iodine_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! ---------------------------------------------------------------------- 
    ! Get species IDs and settings 
    ! ---------------------------------------------------------------------- 
  
    ! Read settings specified in configuration file
    ! Note: the specified strings have to match those in 
    !       the config. file!
    CALL GetExtOpt ( ExtNr, 'Emit I2', OptValBool=CalcI2, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL GetExtOpt ( ExtNr, 'Emit HOI', OptValBool=CalcHOI, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    # 
    minLen = 0
    WRITE(*,*) 'tms debug 1', CalcI2, CalcHOI, minLen

    IF ( CalcI2 ) THEN
       minLen = minLen +1
    ENDIF
    IF ( CalcHOI ) THEN
       minLen = minLen +1
    ENDIF    

    WRITE(*,*) 'tms debug 2', CalcI2, CalcHOI, minLen

    ! Get HEMCO species IDs
    CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( nSpc < minLen ) THEN
       MSG = 'Not enough iodine emission species set' 
       CALL HCO_ERROR ( MSG, RC ) 
       RETURN
    ENDIF

    WRITE(*,*) 'tms debug 3', IDTI2, IDTHOI, HcoIDs
    WRITE(*,*) 'tms debug 3.1', nSpc, SpcNames
    IDTHOI = HcoIDs(1)
    IDTI2 = HcoIDs(2)
    WRITE(*,*) 'tms debug 4', IDTI2, IDTHOI, HcoIDs
    WRITE(*,*) 'tms debug 3.1', nSpc, SpcNames

    ! Final I2/HOI flag
    CalcI2 = ( CalcI2 .AND. IDTI2 > 0 )
    CalcHOI = ( CalcHOI .AND. IDTHOI > 0 )

    ! Verbose mode
    IF ( am_I_Root ) THEN
       MSG = 'Use inorganic iodine emissions (extension module)'
       CALL HCO_MSG( MSG, SEP1='-' )

       WRITE(*,*) 'tms debug 5', IDTI2, IDTHOI, HcoIDs   
       WRITE(*,*) 'tms debug 5.1',TRIM(SpcNames(1)), TRIM(SpcNames(2))

       IF ( CalcHOI ) THEN
          WRITE(MSG,*) 'HOI: ', TRIM(SpcNames(1)), IDTHOI
          CALL HCO_MSG(MSG)
       ENDIF
   
       IF ( CalcI2 ) THEN
          WRITE(MSG,*) 'I2: ', TRIM(SpcNames(2)), IDTI2
          CALL HCO_MSG(MSG)
       ENDIF
    ENDIF

    !=======================================================================
    ! Create diagnostics. The number densities of both modes are always
    ! written into a diagnostics so that they can be used by other routines
    ! and from outside of HEMCO. These two diagnostics just hold a pointer
    ! to the respective density arrays filled by the run method of this
    ! module.
    !=======================================================================
!    CALL Diagn_Create ( am_I_Root,                          &
!                        HcoState   = HcoState,              & 
!                        cName      = 'SEASALT_DENS_FINE',   &
!                        ExtNr      = ExtNr,                 &
!                        Cat        = -1,                    &
!                        Hier       = -1,                    &
!                        HcoID      = IDTSALA,               &
!                        SpaceDim   = 2,                     &
!                        OutUnit    = 'number_dens',         &
!                        AutoFill   = 0,                     &
!                        Trgt2D     = NDENS_SALA,            &
!                        COL        = HcoDiagnIDManual,      &
!                        RC         = RC                      )
!    IF ( RC /= HCO_SUCCESS ) RETURN

!   CALL Diagn_Create ( am_I_Root,                                 &
!                       cName    = 'PARANOX_HNO3_DEPOSITION_FLUX', &
!                       Trgt2D   = DEPHNO3,                        &
!                       SpaceDim = 2,                              &
!                       OutUnit  = 'kg/m2/s',                      &
!                       COL      = HcoDiagnIDManual,               &
!                       RC       = RC                               )
!   IF ( RC /= HCO_SUCCESS ) RETURN


!    CALL Diagn_Create ( am_I_Root,                          &
!                        HcoState   = HcoState,              & 
!                        cName      = 'HOI',   &
!                        ExtNr      = ExtNr,                 &
!                        Cat        = -1,                    &
!                        Hier       = -1,                    &
!                        HcoID      = IDT,                   &
!                        SpaceDim   = 2,                     &
!                        OutUnit    = 'number_dens',         &
!                        AutoFill   = 0,                     &
!                        Trgt2D     = NDENS_SALA,            &
!                        COL        = HcoDiagnIDManual,      &
!                        RC         = RC                      )
!    IF ( RC /= HCO_SUCCESS ) RETURN


    !=======================================================================
    ! Activate this module and the fields of ExtState that it uses
    !=======================================================================

    ! Activate met fields used by this module
    ExtState%WLI%DoUse   = .TRUE.
    ExtState%ALBD%DoUse  = .TRUE.
    ExtState%TSKIN%DoUse = .TRUE.
    ExtState%U10M%DoUse  = .TRUE.
    ExtState%V10M%DoUse  = .TRUE.
    ExtState%O3%DoUse    = .TRUE.
    ExtState%AIR%DoUse   = .TRUE.

    ! Enable module
    ExtState%Inorg_Iodine = .TRUE.

    ! Return w/ success
    IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  )
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
    CALL HCO_LEAVE ( RC ) 
 
  END SUBROUTINE HCOX_Iodine_Init

!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Iodine_Final 
!
! !DESCRIPTION: Subroutine HcoX\_Iodine\_Final deallocates 
!  all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Iodine_Final
!
! !REVISION HISTORY:
!  15 Mar 2013 - T. Sherwen - Initial implementation (v9-3-01)
!  15 Jul 2015 - T. Sherwen - Now a HEMCO extension module adapted from hcox_seasalt_final
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !=================================================================
    ! HCOX_Iodine_Final begins here!
    !=================================================================

    ! Cleanup module arrays
!    IF ( ALLOCATED ( HcoIDs     ) ) DEALLOCATE( HcoIDs      )
!    IF ( ALLOCATED ( SpcNames   ) ) DEALLOCATE( SpcNames    )

  END SUBROUTINE HCOX_Iodine_Final
!EOC
END MODULE HCOX_Iodine_Mod
