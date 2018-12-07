!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_gc_RnPbBe_mod.F90
!
! !DESCRIPTION: Defines the HEMCO extension for the GEOS-Chem Rn-Pb-Be 
!  specialty simulation. 
!\\
!\\
!  This extension parameterizes emissions of Rn and/or Pb based upon the
!  literature given below. The emission fields become automatically added 
!  to the HEMCO emission array of the given species. It is possible to
!  select only one of the two species (Rn or Pb) in the HEMCO configuration
!  file. This may be useful if a gridded data inventory shall be applied to
!  one of the species (through the standard HEMCO interface).
!\\
!\\
! !INTERFACE:
!
MODULE HCOX_GC_RnPbBe_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_Diagn_Mod
  USE HCO_State_Mod,  ONLY : HCO_State   ! Derived type for HEMCO state
  USE HCOX_State_Mod, ONLY : Ext_State   ! Derived type for External state

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HcoX_GC_RnPbBe_Run
  PUBLIC  :: HcoX_GC_RnPbBe_Init
  PUBLIC  :: HcoX_Gc_RnPbBe_Final
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Init_7Be_Emissions
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Liu,H., D.Jacob, I.Bey, and R.M.Yantosca, Constraints from 210Pb 
!        and 7Be on wet deposition and transport in a global three-dimensional
!        chemical tracer model driven by assimilated meteorological fields, 
!        JGR, 106, D11, 12,109-12,128, 2001.
!  (2 ) Jacob et al.,Evaluation and intercomparison of global atmospheric 
!        transport models using Rn-222 and other short-lived tracers, 
!        JGR, 1997 (102):5953-5970
!  (3 ) Dorothy Koch, JGR 101, D13, 18651, 1996.
!  (4 ) Lal, D., and B. Peters, Cosmic ray produced radioactivity on the 
!        Earth. Handbuch der Physik, 46/2, 551-612, edited by K. Sitte, 
!        Springer-Verlag, New York, 1967.
!
! !REVISION HISTORY:
!  07 Jul 2014 - R. Yantosca - Initial version
!  15 Aug 2014 - C. Keller   - Targets now in hp precision. Cosmetic changes
!  21 Aug 2014 - R. Yantosca - Add Pb as a species
!  21 Aug 2014 - R. Yantosca - Add HEMCO species indices as module variables
!  04 Sep 2014 - R. Yantosca - Remove IDTPb; Pb210 only has a chemical source
!  04 Sep 2014 - R. Yantosca - Modified for GCAP simulation
!  05 Nov 2014 - C. Keller   - Now allow Rn or Pb to be not specified.
!  07 Jan 2016 - E. Lundgren - Update Avogadro's # to NIST 2014 value
!  24 Aug 2017 - M. Sulprizio- Remove support for GCAP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Emissions indices etc.
  INTEGER                       :: ExtNr  = -1       ! HEMCO Extension number
  INTEGER                       :: IDTRn  = -1       ! Tracer index # for Rn
  INTEGER                       :: IDTBe7 = -1       ! Tracer index # for Be7

  ! For tracking Rn222 and Be7 emissions
  REAL(hp), ALLOCATABLE, TARGET :: EmissRn (:,:  )
  REAL(hp), ALLOCATABLE, TARGET :: EmissBe7(:,:,:)

  ! For Lal & Peters 7Be emissions input data
  REAL(hp), ALLOCATABLE         :: LATSOU  (:    )   ! Array for latitudes
  REAL(hp), ALLOCATABLE         :: PRESOU  (:    )   ! Array for pressures
  REAL(hp), ALLOCATABLE         :: BESOU   (:,:  )   ! Array for 7Be emissions
!
! !DEFINED PARAMETERS:
!
  ! To convert kg to atoms
  REAL*8,  PARAMETER            :: XNUMOL_Rn = ( 6.022140857d23 / 222.0d-3 )    
  REAL*8,  PARAMETER            :: XNUMOL_Be = ( 6.022140857d23 /   7.0d-3 )

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Gc_RnPbBe_run 
!
! !DESCRIPTION: Subroutine HcoX\_Gc\_RnPbBe\_Run computes emissions of 222Rn
!  and 7Be for the GEOS-Chem Rn-Pb-Be specialty simulation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Gc_RnPbBe_Run( am_I_Root, ExtState, HcoState, RC )
!
! !USES:
!
    ! HEMCO modules
    USE HCO_FluxArr_Mod, ONLY : HCO_EmisAdd
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root   ! Are we on the root CPU?
    TYPE(Ext_State),  POINTER       :: ExtState    ! Options for Rn-Pb-Be sim
    TYPE(HCO_State),  POINTER       :: HcoState    ! HEMCO state 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  This code is based on routine EMISSRnPbBe in prior versions of GEOS-Chem.
!
! !REVISION HISTORY:
!  07 Jul 2014 - R. Yantosca - Initial version
!  03 Sep 2014 - R. Yantosca - Bug fix: Prevent div-by-zero errors
!  06 Oct 2014 - C. Keller   - Now calculate pressure centers from edges.
!  29 Oct 2014 - R. Yantosca - Use latitude centers of the grid box to
!                              facilitate running in ESMF/MPI environment
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER           :: I,        J,          L,          N
    INTEGER           :: HcoID
    REAL*8            :: A_CM2,    ADD_Be,     ADD_Rn,     Rn_LAND
    REAL*8            :: Rn_WATER, DTSRCE
    REAL*8            :: Rn_TMP,     LAT,        F_LAND     
    REAL*8            :: F_WATER,  F_BELOW_70, F_BELOW_60, F_ABOVE_60
    REAL*8            :: DENOM
    REAL(hp)          :: LAT_TMP,  P_TMP,      Be_TMP

    ! Pointers
    REAL(hp), POINTER :: Arr2D(:,:  )
    REAL(hp), POINTER :: Arr3D(:,:,:)

    !=======================================================================
    ! HCOX_GC_RnPbBe_RUN begins here!
    !=======================================================================

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_GC_RnPbBe_Run (hcox_gc_RnPbBe_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Set error flag
    !ERR = .FALSE.

    ! Sanity check: return if extension not turned on
    IF ( ExtState%Gc_RnPbBe < 0 ) RETURN

    ! Emission timestep [s]
    DTSRCE = HcoState%TS_EMIS 

    ! Nullify
    Arr2D => NULL()
    Arr3D => NULL()

    !=======================================================================
    ! Compute 222Rn emissions [kg/m2/s], according to the following:
    !
    ! (1) 222Rn emission poleward of 70 degrees = 0.0 [atoms/cm2/s]
    ! 
    ! (2) For latitudes 70S-60S and 60N-70N (both land & ocean),
    !     222Rn emission is 0.005 [atoms/cm2/s]
    !
    ! (3) For latitudes between 60S and 60N, 
    !     222Rn emission is 1     [atoms/cm2/s] over land or
    !                       0.005 [atoms/cm2/s] over oceans
    !
    ! (4) For grid boxes where the surface temperature is below 
    !     0 deg Celsius, reduce 222Rn emissions by a factor of 3.
    ! 
    ! Reference: Jacob et al.,Evaluation and intercomparison of 
    !  global atmospheric transport models using Rn-222 and other 
    !  short-lived tracers, JGR, 1997 (102):5953-5970
    !=======================================================================
    IF ( IDTRn > 0 ) THEN

!$OMP PARALLEL DO                                            &
!$OMP DEFAULT( SHARED )                                       &
!$OMP PRIVATE( I,          J,          LAT,        DENOM   ) &
!$OMP PRIVATE( F_BELOW_70, F_BELOW_60, F_ABOVE_60, Rn_LAND ) &
!$OMP PRIVATE( Rn_WATER,   F_LAND,     F_WATER,    ADD_Rn  ) &
!$OMP SCHEDULE( DYNAMIC )
       DO J = 1, HcoState%Ny
       DO I = 1, HcoState%Nx
   
          ! Get ABS( latitude ) of the grid box
          LAT           = ABS( HcoState%Grid%YMID%Val( I, J ) )
   
          ! Zero for safety's sake
          F_BELOW_70    = 0d0
          F_BELOW_60    = 0d0
          F_ABOVE_60    = 0d0
   
          ! Baseline 222Rn emissions 
          ! Rn_LAND [kg/m2/s] = [1 atom 222Rn/cm2/s] / [atoms/kg] * [1d4 cm2/m2]
          Rn_LAND       = ( 1d0 / XNUMOL_Rn ) * 1d4
   
          ! Baseline 222Rn emissions over water or ice [kg]
          Rn_WATER      = Rn_LAND * 0.005d0
   
          ! Fraction of grid box that is land
          F_LAND        = ExtState%FRCLND%Arr%Val(I,J)
   
          ! Fraction of grid box that is water
          F_WATER       = 1d0 - F_LAND
   
          !--------------------
          ! 90S-70S or 70N-90N
          !--------------------
          IF ( LAT >= 70d0 ) THEN 
   
             ! 222Rn emissions are shut off poleward of 70 degrees
             ADD_Rn = 0.0d0
   
          !--------------------
          ! 70S-60S or 60N-70N 
          !--------------------
          ELSE IF ( LAT >= 60d0 ) THEN    
   
             IF ( LAT <= 70d0 ) THEN             
   
                ! If the entire grid box lies equatorward of 70 deg,
                ! then 222Rn emissions here are 0.005 [atoms/cm2/s]
                ADD_Rn = Rn_WATER
                  
             ELSE
   
                ! N-S extent of grid box [degrees]
                DENOM = HcoState%Grid%YMID%Val( I, J+1 ) &
                      - HcoState%Grid%YMID%Val( I, J   )               
   
                ! Compute the fraction of the grid box below 70 degrees
                F_BELOW_70 = ( 70.0d0 - LAT ) / DENOM
   
                ! If the grid box straddles the 70S or 70N latitude line,
                ! then only count 222Rn emissions equatorward of 70 degrees.
                ! 222Rn emissions here are 0.005 [atoms/cm2/s].
                ADD_Rn = F_BELOW_70 * Rn_WATER
                  
             ENDIF
               
          ELSE 
   
             !--------------------
             ! 70S-60S or 60N-70N
             !--------------------
             IF ( LAT > 60d0 ) THEN
                
                ! N-S extent of grid box [degrees]
                DENOM  = HcoState%Grid%YMID%Val( I, J+1 ) &
                       - HcoState%Grid%YMID%Val( I, J   )
   
                ! Fraction of grid box with ABS( lat ) below 60 degrees
                F_BELOW_60 = ( 60.0d0 - LAT ) / DENOM
   
                ! Fraction of grid box with ABS( lat ) above 60 degrees
                F_ABOVE_60 = F_BELOW_60
                
                ADD_Rn =                                                &
                         ! Consider 222Rn emissions equatorward of 
                         ! 60 degrees for both land (1.0 [atoms/cm2/s]) 
                         ! and water (0.005 [atoms/cm2/s])
                         F_BELOW_60 *                                   &
                         ( Rn_LAND  * F_LAND  ) +                       &
                         ( Rn_WATER * F_WATER ) +                       &
   
                         ! If the grid box straddles the 60 degree boundary
                         ! then also consider the emissions poleward of 60
                         ! degrees.  222Rn emissions here are 0.005 [at/cm2/s].
                         F_ABOVE_60 * Rn_WATER                    
   
   
             !--------------------
             ! 60S-60N
             !--------------------
             ELSE 
                  
                ! Consider 222Rn emissions equatorward of 60 deg for
                ! land (1.0 [atoms/cm2/s]) and water (0.005 [atoms/cm2/s])
                ADD_Rn = ( Rn_LAND * F_LAND ) + ( Rn_WATER * F_WATER )
   
             ENDIF
          ENDIF
   
          ! For boxes below freezing, reduce 222Rn emissions by 3x
          IF ( ExtState%T2M%Arr%Val(I,J) < 273.15 ) THEN
             ADD_Rn = ADD_Rn / 3d0
          ENDIF
   
          ! Save 222Rn emissions into an array [kg/m2/s]
          EmissRn(I,J) = ADD_Rn
       ENDDO
       ENDDO
!$OMP END PARALLEL DO
   
       !------------------------------------------------------------------------
       ! Add 222Rn emissions to HEMCO data structure & diagnostics
       !------------------------------------------------------------------------
   
       ! Add emissions
       Arr2D => EmissRn(:,:)
       CALL HCO_EmisAdd( am_I_Root, HcoState, Arr2D, IDTRn, RC, ExtNr=ExtNr )
       Arr2D => NULL()
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'HCO_EmisAdd error: EmisRn', RC )
          RETURN 
       ENDIF
   
    ENDIF ! IDTRn > 0

    !=======================================================================
    ! Compute 7Be emissions [kg/m2/s]    
    !
    ! Original units of 7Be emissions are [stars/g air/sec],
    ! where "stars" = # of nuclear disintegrations of cosmic rays
    !
    ! Now interpolate from 33 std levels onto GEOS-CHEM levels 
    !=======================================================================
    IF ( IDTBe7 > 0 ) THEN
!$OMP PARALLEL DO                                        &
!$OMP DEFAULT( SHARED )                                  &
!$OMP PRIVATE( I, J, L, LAT_TMP, P_TMP, Be_TMP, ADD_Be ) &
!$OMP SCHEDULE( DYNAMIC )
       DO L = 1, HcoState%Nz
       DO J = 1, HcoState%Ny
       DO I = 1, HcoState%Nx
   
          ! Get absolute value of latitude, since we will assume that 
          ! the 7Be distribution is symmetric about the equator
          LAT_TMP = ABS( HcoState%Grid%YMID%Val( I, J ) )
   
          ! Pressure at (I,J,L) [hPa]
          ! Now calculate from edge points (ckeller, 10/06/1014)
          P_TMP = ( HcoState%Grid%PEDGE%Val(I,J,L) + &
                    HcoState%Grid%PEDGE%Val(I,J,L+1) ) / 200.0_hp 
                    
          ! Interpolate 7Be [stars/g air/sec] to GEOS-Chem levels
          CALL SLQ( LATSOU, PRESOU, BESOU, 10, 33, LAT_TMP, P_TMP, Be_TMP )
   
          ! Be_TMP = [stars/g air/s] * [0.045 atom/star] * 
          !          [kg air] * [1e3 g/kg] = 7Be emissions [atoms/s]
          Be_TMP  = Be_TMP * 0.045e+0_hp * ExtState%AIR%Arr%Val(I,J,L) * 1.e+3_hp 
                     
          ! ADD_Be = [atoms/s] / [atom/kg] / [m2] = 7Be emissions [kg/m2/s]
          ADD_Be  = ( Be_TMP / XNUMOL_Be ) / HcoState%Grid%AREA_M2%Val(I,J)
   
          ! Save emissions into an array for use below
          EmissBe7(I,J,L) = ADD_Be
   
       ENDDO
       ENDDO
       ENDDO
!$OMP END PARALLEL DO
   
       !------------------------------------------------------------------------
       ! Add 7Be emissions to HEMCO data structure & diagnostics
       !------------------------------------------------------------------------
   
       ! Add emissions
       Arr3D => EmissBe7(:,:,:)
       CALL HCO_EmisAdd( am_I_Root, HcoState, Arr3D, IDTBe7, RC, ExtNr=ExtNr )
       Arr3D => NULL()
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'HCO_EmisAdd error: EmisBe7', RC )
          RETURN 
       ENDIF
   
    ENDIF !IDTBe7 > 0

    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Return w/ success
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_Gc_RnPbBe_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Gc_RnPbBe_Init 
!
! !DESCRIPTION: Subroutine HcoX\_Gc\_RnPbBe\_Init initializes the HEMCO
! GC\_Rn-Pb-Be extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Gc_RnPbBe_Init( am_I_Root, HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_ExtList_Mod, ONLY : GetExtNr
    USE HCO_State_Mod,   ONLY : HCO_GetExtHcoID
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName     ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState    ! Module options      
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState    ! Hemco state 
    INTEGER,          INTENT(INOUT)  :: RC 

! !REVISION HISTORY:
!  07 Jul 2014 - R. Yantosca - Initial version
!  21 Aug 2014 - R. Yantosca - Now define HEMCO indices as well
!  04 Sep 2014 - R. Yantosca - Activate ExtState%TROPP for GCAP simulation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                        :: N, nSpc
    CHARACTER(LEN=255)             :: MSG 

    ! Arrays
    INTEGER,           ALLOCATABLE :: HcoIDs(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)

    !=======================================================================
    ! HCOX_GC_RnPbBe_INIT begins here!
    !=======================================================================

    ! Get the extension number
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter HEMCO
    CALL HCO_ENTER( HcoState%Config%Err, 'HcoX_GC_RnPbBe_Init (hcox_gc_RnPbBe_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Set species IDs      
    CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Verbose mode
    IF ( am_I_Root ) THEN
       MSG = 'Use gc_RnPbBe emissions module (extension module)'
       CALL HCO_MSG(HcoState%Config%Err,MSG )

       MSG = 'Use the following species (Name: HcoID):'
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       DO N = 1, nSpc
          WRITE(MSG,*) TRIM(SpcNames(N)), ':', HcoIDs(N)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDDO
    ENDIF

    ! Set up tracer and HEMCO indices
    DO N = 1, nSpc
       SELECT CASE( TRIM( SpcNames(N) ) )
          CASE( 'Rn', 'Rn222', '222Rn' )
             IDTRn   = HcoIDs(N)
          CASE( 'Be', 'Be7', '7Be' )
             IDTBe7  = HcoIDs(N)
          CASE DEFAULT
             ! Do nothing
       END SELECT
    ENDDO

    ! WARNING: Rn tracer is not found!
    IF ( IDTRn <= 0 .AND. am_I_Root ) THEN
       CALL HCO_WARNING(HcoState%Config%Err, 'Cannot find 222Rn tracer in list of species!', RC )
    ENDIF
    
    ! WARNING: Be7 tracer is not found
    IF ( IDTBe7 <= 0 .AND. am_I_Root ) THEN
       CALL HCO_WARNING(HcoState%Config%Err, 'Cannot find 7Be tracer in list of species!', RC )
    ENDIF

    ! ERROR: No tracer defined
    IF ( IDTRn <= 0 .AND. IDTBe7 <= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot use RnPbBe extension: no valid species!', RC )
    ENDIF

    ! Activate met fields required by this extension
    ExtState%FRCLND%DoUse  = .TRUE. 
    ExtState%T2M%DoUse     = .TRUE. 
    ExtState%AIR%DoUse     = .TRUE. 

    ! Activate this extension
    ExtState%Gc_RnPbBe     = .TRUE.

    !=======================================================================
    ! Initialize data arrays
    !=======================================================================

    IF ( IDTRn > 0 ) THEN
       ALLOCATE( EmissRn( HcoState%Nx, HcoState%NY ), STAT=RC )
       IF ( RC /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate EmissRn', RC )
          RETURN
       ENDIF 
    ENDIF

    IF ( IDTBe7 > 0 ) THEN
       ALLOCATE( EmissBe7( HcoState%Nx, HcoState%NY, HcoState%NZ ), STAT=RC )
       IF ( RC /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate EmissBe7', RC )
          RETURN
       ENDIF 
       IF ( RC /= 0 ) RETURN

       ! Array for latitudes (Lal & Peters data)
       ALLOCATE( LATSOU( 10 ), STAT=RC )
       IF ( RC /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate LATSOU', RC )
          RETURN
       ENDIF 
   
       ! Array for pressures (Lal & Peters data)
       ALLOCATE( PRESOU( 33 ), STAT=RC )
       IF ( RC /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate PRESOU', RC )
          RETURN
       ENDIF 
   
       ! Array for 7Be emissions ( Lal & Peters data)
       ALLOCATE( BESOU( 10, 33 ), STAT=RC )
       IF ( RC /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate BESOU', RC )
          RETURN
       ENDIF 
       
       ! Initialize the 7Be emisisons data arrays
       CALL Init_7Be_Emissions()
    ENDIF

    !=======================================================================
    ! Leave w/ success
    !=======================================================================
    IF ( ALLOCATED( HcoIDs   ) ) DEALLOCATE( HcoIDs   )
    IF ( ALLOCATED( SpcNames ) ) DEALLOCATE( SpcNames )

    CALL HCO_LEAVE( HcoState%Config%Err,RC ) 

  END SUBROUTINE HCOX_Gc_RnPbBe_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Gc_RnPbBe_Final
!
! !DESCRIPTION: Subroutine HcoX\_Gc\_RnPbBe\_Final finalizes the HEMCO
!  extension for the GEOS-Chem Rn-Pb-Be specialty simulation.  All module
!  arrays will be deallocated.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Gc_RnPbBe_Final()
!
! !REVISION HISTORY:
!  13 Dec 2013 - C. Keller   - Now a HEMCO extension
!  06 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  06 Jun 2014 - R. Yantosca - Now indended with F90 free-format
!EOP
!------------------------------------------------------------------------------
!BOC

    !=======================================================================
    ! HCOX_GC_RNPBBE_FINAL begins here!
    !=======================================================================
    IF ( ALLOCATED( EmissRn  ) ) DEALLOCATE( EmissRn  ) 
    IF ( ALLOCATED( EmissBe7 ) ) DEALLOCATE( EmissBe7 )
    IF ( ALLOCATED( LATSOU   ) ) DEALLOCATE( LATSOU   )
    IF ( ALLOCATED( PRESOU   ) ) DEALLOCATE( PRESOU   )
    IF ( ALLOCATED( BESOU    ) ) DEALLOCATE( BESOU    )

  END SUBROUTINE HCOX_Gc_RnPbBe_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_7Be_Emissions
!
! !DESCRIPTION: Subroutine Init\_7Be\_Emissions initializes the 7Be emissions 
!  from Lal \& Peters on 33 pressure levels.  This data used to be read from 
!  a file, but we have now hardwired it to facilitate I/O in the ESMF 
!  environment.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_7Be_Emissions()
!
! !REMARKS:
!  (1) Reference: Lal, D., and B. Peters, Cosmic ray produced radioactivity 
!       on the Earth. Handbuch der Physik, 46/2, 551-612, edited by K. Sitte, 
!        Springer-Verlag, New York, 1967.
!                                                                             .
!  (2) In prior versions of GEOS-Chem, this routine was named READ_7BE, and
!      it read the ASCII file "7Be.Lal".   Because this data set is not placed
!      on a lat/lon grid, ESMF cannot regrid it.  To work around this, we now
!      hardwire this data in module arrays rather than read it from disk.
!                                                                             .
!  (3) Units of 7Be emissions are [stars/g air/s].  
!      Here, "stars" = # of nuclear disintegrations of cosmic rays
!                                                                             .
!  (4) Original data from Lal & Peters (1967), w/ these modifications:
!      (a) Replace data at (0hPa, 70S) following Koch 1996:
!          (i ) old value = 3000 
!          (ii) new value = 1900
!      (b) Copy data from 70S to 80S and 90S at all levels
!                                                                             .
! !REVISION HISTORY: 
!  07 Aug 2002 - H. Liu - Initial version
!  (1 ) This code was split off from routine EMISSRnPbBe below. (bmy, 8/7/02)
!  (2 ) Now reference DATA_DIR from "directory_mod.f" (bmy, 7/19/04)
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  02 Jul 2014 - R. Yantosca - Now hardwire the data instead of reading it
!                              from an ASCII file.  This facilitates ESMF I/O.
!  07 Jul 2014 - R. Yantosca - Now renamed to INIT_7Be_Emissions and added
!                              as a HEMCO extension
!  07 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!   8 Aug 2014 - R. Yantosca - Now split off into hcox_gc_RnPbBe_include.H
!  05 Nov 2014 - C. Keller   - Converted from double-precision to flexible
!                              (HEMCO) precision hp.
!  26 Feb 2015 - R. Yantosca - Now inline the code that used to be in the
!                              include file hcox_gc_RnPbBe_include.H.  This
!                              will result in faster compilation.
!  08 Jan 2016 - R. Yantosca - Change 54_hp to 54.0_hp to avoid error
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Define latitudes [degrees North]
    LATSOU      = (/     0.0_hp,     10.0_hp,     20.0_hp,    30.0_hp,    &
                        40.0_hp,     50.0_hp,     60.0_hp,    70.0_hp,    &
                        80.0_hp,     90.0_hp  /)
 
    ! Define pressures [hPa]
    PRESOU      = (/     0.0_hp,     50.0_hp,     70.0_hp,    90.0_hp,    &
                       110.0_hp,    130.0_hp,    150.0_hp,   170.0_hp,    &
                       190.0_hp,    210.0_hp,    230.0_hp,   250.0_hp,    &
                       270.0_hp,    290.0_hp,    313.0_hp,   338.0_hp,    &
                       364.0_hp,    392.0_hp,    420.0_hp,   451.0_hp,    &
                       485.0_hp,    518.0_hp,    555.0_hp,   592.0_hp,    &
                       633.0_hp,    680.0_hp,    725.0_hp,   772.0_hp,    &
                       822.0_hp,    875.0_hp,    930.0_hp,   985.0_hp,    &
                      1030.0_hp  /)

    ! Define 7Be emissions [stars/g air/s]
    ! 1 "star" = 1 nuclear disintegration via cosmic rays
    !
    ! NOTE: These statements were defined from printout of the file
    ! and need to be multiplied by 1d-5 below.
    BESOU(:,1)  = (/   150.0_hp,    156.0_hp,    188.0_hp,   285.0_hp,    &
                       500.0_hp,    910.0_hp,   1700.0_hp,  1900.0_hp,    &
                      1900.0_hp,   1900.0_hp  /)
                       
    BESOU(:,2)  = (/   280.0_hp,    310.0_hp,    390.0_hp,   590.0_hp,    &
                       880.0_hp,   1390.0_hp,   1800.0_hp,  1800.0_hp,    &
                      1800.0_hp,   1800.0_hp  /)
                       
    BESOU(:,3)  = (/   310.0_hp,    330.0_hp,    400.0_hp,   620.0_hp,    &
                       880.0_hp,   1280.0_hp,   1450.0_hp,  1450.0_hp,    &
                      1450.0_hp,   1450.0_hp  /)
                       
    BESOU(:,4)  = (/   285.0_hp,    310.0_hp,    375.0_hp,   570.0_hp,    &
                       780.0_hp,   1100.0_hp,   1180.0_hp,  1180.0_hp,    &
                      1180.0_hp,   1180.0_hp  /)
                       
    BESOU(:,5)  = (/   255.0_hp,    275.0_hp,    330.0_hp,   510.0_hp,    &
                       680.0_hp,    950.0_hp,   1000.0_hp,  1000.0_hp,    &
                      1000.0_hp,   1000.0_hp  /)
                       
    BESOU(:,6)  = (/   230.0_hp,    245.0_hp,    292.0_hp,   450.0_hp,    &
                       600.0_hp,    820.0_hp,    875.0_hp,   875.0_hp,    &
                       875.0_hp,    875.0_hp  /)
                       
    BESOU(:,7)  = (/   205.0_hp,    215.0_hp,    260.0_hp,   400.0_hp,    &
                       530.0_hp,    730.0_hp,    750.0_hp,   750.0_hp,    &
                       750.0_hp,    750.0_hp  /)
                       
    BESOU(:,8)  = (/   182.0_hp,    195.0_hp,    235.0_hp,   355.0_hp,    &
                       480.0_hp,    630.0_hp,    650.0_hp,   650.0_hp,    &
                       650.0_hp,    650.0_hp  /)
                       
    BESOU(:,9)  = (/   160.0_hp,    173.0_hp,    208.0_hp,   315.0_hp,    &
                       410.0_hp,    543.0_hp,    550.0_hp,   550.0_hp,    &
                       550.0_hp,    550.0_hp  /)
                       
    BESOU(:,10) = (/   148.0_hp,    152.0_hp,    185.0_hp,   280.0_hp,    &
                       370.0_hp,    480.0_hp,    500.0_hp,   500.0_hp,    &
                       500.0_hp,    500.0_hp  /)
                       
    BESOU(:,11) = (/   130.0_hp,    139.0_hp,    167.0_hp,   250.0_hp,    &
                       320.0_hp,    425.0_hp,    430.0_hp,   430.0_hp,    &
                       430.0_hp,    430.0_hp  /)
                       
    BESOU(:,12) = (/   116.0_hp,    123.0_hp,    148.0_hp,   215.0_hp,    &
                       285.0_hp,    365.0_hp,    375.0_hp,   375.0_hp,    &
                       375.0_hp,    375.0_hp  /)
                       
    BESOU(:,13) = (/   104.0_hp,    110.0_hp,    130.0_hp,   198.0_hp,    &
                       250.0_hp,    320.0_hp,    330.0_hp,   330.0_hp,    &
                       330.0_hp,    330.0_hp  /)
                       
    BESOU(:,14) = (/    93.0_hp,     99.0_hp,    118.0_hp,   170.0_hp,    &
                       222.0_hp,    280.0_hp,    288.0_hp,   288.0_hp,    &
                       288.0_hp,    288.0_hp  /)
                       
    BESOU(:,15) = (/    80.0_hp,     84.0_hp,    100.0_hp,   145.0_hp,    &
                       190.0_hp,    235.0_hp,    250.0_hp,   250.0_hp,    &
                       250.0_hp,    250.0_hp  /)
                       
    BESOU(:,16) = (/    72.0_hp,     74.0_hp,     88.0_hp,   129.0_hp,    &
                       168.0_hp,    210.0_hp,    218.0_hp,   218.0_hp,    &
                       218.0_hp,    218.0_hp  /)
                       
    BESOU(:,17) = (/    59.5_hp,     62.5_hp,     73.5_hp,   108.0_hp,    &
                       138.0_hp,    171.0_hp,    178.0_hp,   178.0_hp,    &
                       178.0_hp,    178.0_hp  /)
                       
    BESOU(:,18) = (/    50.0_hp,     53.0_hp,     64.0_hp,    90.0_hp,    &
                       115.0_hp,    148.0_hp,    150.0_hp,   150.0_hp,    &
                       150.0_hp,    150.0_hp  /)
                       
    BESOU(:,19) = (/    45.0_hp,     46.5_hp,     52.5_hp,    76.0_hp,    &
                        98.0_hp,    122.0_hp,    128.0_hp,   128.0_hp,    &
                       128.0_hp,    128.0_hp  /)
                       
    BESOU(:,20) = (/    36.5_hp,     37.5_hp,     45.0_hp,    61.0_hp,    &
                        77.0_hp,     98.0_hp,    102.0_hp,   102.0_hp,    &
                       102.0_hp,    102.0_hp  /)
                       
    BESOU(:,21) = (/    30.8_hp,     32.0_hp,     37.5_hp,    51.5_hp,    &
                        65.0_hp,     81.0_hp,     85.0_hp,    85.0_hp,    &
                        85.0_hp,     85.0_hp  /)

    BESOU(:,22) = (/    25.5_hp,     26.5_hp,     32.0_hp,    40.5_hp,    &
                        54.0_hp,     67.5_hp,     69.5_hp,    69.5_hp,    &
                        69.5_hp,     69.5_hp  /)

    BESOU(:,23) = (/    20.5_hp,     21.6_hp,     25.5_hp,    33.0_hp,    &
                        42.0_hp,     53.5_hp,     55.0_hp,    55.0_hp,    &
                        55.0_hp,     55.0_hp  /)

    BESOU(:,24) = (/    16.8_hp,     17.3_hp,     20.0_hp,    26.0_hp,    &
                        33.5_hp,     41.0_hp,     43.0_hp,    43.0_hp,    &
                        43.0_hp,     43.0_hp  /)

    BESOU(:,25) = (/    13.0_hp,     13.8_hp,     15.3_hp,    20.5_hp,    &
                        26.8_hp,     32.5_hp,     33.5_hp,    33.5_hp,    &
                        33.5_hp,     33.5_hp  /)

    BESOU(:,26) = (/    10.1_hp,     10.6_hp,     12.6_hp,    15.8_hp,    &
                        20.0_hp,     24.5_hp,     25.8_hp,    25.8_hp,    &
                        25.8_hp,     25.8_hp  /)

    BESOU(:,27) = (/     7.7_hp,     8.15_hp,      9.4_hp,    11.6_hp,    &
                        14.8_hp,     17.8_hp,     18.5_hp,    18.5_hp,    &
                        18.5_hp,     18.5_hp  /)
 
    BESOU(:,28) = (/     5.7_hp,     5.85_hp,     6.85_hp,    8.22_hp,    &
                        11.0_hp,     13.1_hp,     13.2_hp,    13.2_hp,    &
                        13.2_hp,     13.2_hp  /)

    BESOU(:,29) = (/     3.9_hp,      4.2_hp,     4.85_hp,     6.0_hp,    &
                         7.6_hp,      9.0_hp,      9.2_hp,     9.2_hp,    &
                         9.2_hp,      9.2_hp  /)

    BESOU(:,30) = (/     3.0_hp,     3.05_hp,     3.35_hp,     4.2_hp,    &
                         5.3_hp,      5.9_hp,     6.25_hp,    6.25_hp,    &
                        6.25_hp,     6.25_hp  /)

    BESOU(:,31) = (/    2.05_hp,      2.1_hp,     2.32_hp,     2.9_hp,    &
                         3.4_hp,      3.9_hp,      4.1_hp,     4.1_hp,    &
                         4.1_hp,      4.1_hp  /)

    BESOU(:,32) = (/    1.45_hp,     1.43_hp,     1.65_hp,    2.03_hp,    &
                         2.4_hp,     2.75_hp,     2.65_hp,    2.65_hp,    &
                        2.65_hp,     2.65_hp  /)

    BESOU(:,33) = (/    1.04_hp,     1.08_hp,     1.21_hp,     1.5_hp,    &
                        1.68_hp,      1.8_hp,      1.8_hp,     1.8_hp,    &
                         1.8_hp,      1.8_hp  /)

    ! All the numbers of BESOU need to be multiplied by 1e-5 in order to put 
    ! them into the correct data range.  NOTE: This multiplication statement
    ! needs to be preserved here in order to  ensure identical output to the
    ! prior code! (bmy, 7/7/14)
    BESOU = BESOU * 1.e-5_hp

  END SUBROUTINE Init_7Be_Emissions
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SLQ
!
! !DESCRIPTION: Subroutine SLQ is an interpolation subroutine from a 
!  Chinese reference book (says Hongyu Liu).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SLQ( X, Y, Z, N, M, U, V, W )
!
! !INPUT PARAMETERS: 
!
    INTEGER :: N        ! First dimension of Z
    INTEGER :: M        ! Second dimension of Z
    REAL(hp)  :: X(N)     ! X-axis coordinate on original grid
    REAL(hp)  :: Y(M)     ! Y-axis coordinate on original grid
    REAL(hp)  :: Z(N,M)   ! Array of data on original grid
    REAL(hp)  :: U        ! X-axis coordinate for desired interpolated value
    REAL(hp)  :: V        ! Y-axis coordinate for desired interpolated value
!
! !OUTPUT PARAMETERS:
!
    REAL(hp)  :: W        ! Interpolated value of Z array, at coords (U,V) 
!
! !REMARKS:
!  This routine was taken from the old RnPbBe_mod.F.
! 
! !REVISION HISTORY: 
!  17 Mar 1998 - H. Liu      - Initial version
!  (1 ) Added to "RnPbBe_mod.f" (bmy, 7/16/01)
!  (2 ) Removed duplicate definition of IQ.  Added comments. (bmy, 11/15/01)
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!   7 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  04 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(hp)  :: B(3), HH
    INTEGER :: NN,   IP, I, J, L, IQ, K, MM

    !=======================================================================
    ! SLQ begins here!
    !=======================================================================
    NN=3
    IF(N.LE.3) THEN
       IP=1
       NN=N
    ELSE IF (U.LE.X(2)) THEN
       IP=1
    ELSE IF (U.GE.X(N-1)) THEN
       IP=N-2
    ELSE
       I=1
       J=N
10     IF (IABS(I-J).NE.1) THEN
          L=(I+J)/2
          IF (U.LT.X(L)) THEN
             J=L
          ELSE
             I=L
          END IF
          GOTO 10
       END IF
       IF (ABS(U-X(I)).LT.ABS(U-X(J))) THEN
          IP=I-1
       ELSE
          IP=I
       END IF
    END IF
    MM=3
    IF (M.LE.3) THEN
       IQ=1
       MM=N
    ELSE IF (V.LE.Y(2)) THEN
       IQ=1
    ELSE IF (V.GE.Y(M-1)) THEN
       IQ=M-2
    ELSE
       I=1
       J=M
20     IF (IABS(J-I).NE.1) THEN
          L=(I+J)/2
          IF (V.LT.Y(L)) THEN
             J=L
          ELSE
             I=L
          END IF
          GOTO 20
       END IF
       IF (ABS(V-Y(I)).LT.ABS(V-Y(J))) THEN
          IQ=I-1
       ELSE
          IQ=I
       END IF
    END IF
    DO 50 I=1,NN
       B(I)=0.0
       DO 40 J=1,MM
          HH=Z(IP+I-1,IQ+J-1)
          DO 30 K=1,MM
             IF (K.NE.J) THEN
                HH=HH*(V-Y(IQ+K-1))/(Y(IQ+J-1)-Y(IQ+K-1))
             END IF
30        CONTINUE
          B(I)=B(I)+HH
40     CONTINUE
50  CONTINUE
    W=0.0
    DO 70 I=1,NN
       HH=B(I)
       DO 60 J=1,NN
          IF (J.NE.I) THEN
             HH=HH*(U-X(IP+J-1))/(X(IP+I-1)-X(IP+J-1))
          END IF
60     CONTINUE
        W=W+HH
70   CONTINUE

  END SUBROUTINE SLQ
!EOC
END MODULE HCOX_GC_RnPbBe_Mod
