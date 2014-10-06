!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_gc_RnPbBe_mod.F90
!
! !DESCRIPTION: Defines the HEMCO extension for the GEOS-Chem Rn-Pb-Be 
!  specialty simulation.
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
  REAL*8,  ALLOCATABLE          :: LATSOU  (:    )   ! Array for latitudes
  REAL*8,  ALLOCATABLE          :: PRESOU  (:    )   ! Array for pressures
  REAL*8,  ALLOCATABLE          :: BESOU   (:,:  )   ! Array for 7Be emissions
!
! !DEFINED PARAMETERS:
!
  ! To convert kg to atoms
  REAL*8,  PARAMETER            :: XNUMOL_Rn = ( 6.0225d23 / 222.0d-3 )    
  REAL*8,  PARAMETER            :: XNUMOL_Be = ( 6.0225d23 /   7.0d-3 )

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER           :: I,          J,          L,          N
    INTEGER           :: HcoID
    REAL*8            :: A_CM2,      ADD_Be,     ADD_Rn,     Rn_LAND
    REAL*8            :: Rn_WATER,   DTSRCE,     LAT_TMP,    P_TMP
    REAL*8            :: Be_TMP,     Rn_TMP,     LAT_S,      LAT_N
    REAL*8            :: LAT_H,      LAT_L,      F_LAND,     F_WATER
    REAL*8            :: F_BELOW_70, F_BELOW_60, F_ABOVE_60, DENOM

    ! Pointers
    REAL(hp), POINTER :: Arr2D(:,:  ) => NULL()
    REAL(hp), POINTER :: Arr3D(:,:,:) => NULL()

    !=======================================================================
    ! HCOX_GC_RnPbBe_RUN begins here!
    !=======================================================================

    ! Enter
    CALL HCO_ENTER( 'HCOX_GC_RnPbBe_Run (hcox_gc_RnPbBe_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Set error flag
    !ERR = .FALSE.

    ! Sanity check: return if extension not turned on
    IF ( .NOT. ExtState%Gc_RnPbBe ) RETURN

    ! Emission timestep [s]
    DTSRCE = HcoState%TS_EMIS 

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
!$OMP PARALLEL DO                                                      &
!$OMP DEFAULT( SHARED )                                                &
!$OMP PRIVATE( I,       J,        LAT_S,      LAT_N,      LAT_H      ) &
!$OMP PRIVATE( LAT_L,   DENOM,    F_BELOW_70, F_BELOW_60, F_ABOVE_60 ) &
!$OMP PRIVATE( Rn_LAND, Rn_WATER, F_LAND,     F_WATER,    ADD_Rn     ) &
!$OMP SCHEDULE( DYNAMIC )
    DO J = 1, HcoState%Ny
    DO I = 1, HcoState%Nx

       ! Get ABS( latitude ) at S and N edges of grid box
       LAT_S         = ABS( HcoState%Grid%YEDGE%Val( I, J   ) ) 
       LAT_N         = ABS( HcoState%Grid%YEDGE%Val( I, J+1 ) )
       LAT_H         = MAX( LAT_S, LAT_N )
       LAT_L         = MIN( LAT_S, LAT_N ) 

       ! Grid box extent, for use in denominators below
       DENOM         = ( LAT_H - LAT_L )
       
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
       IF ( LAT_L >= 70d0 ) THEN 

          ! 222Rn emissions are shut off poleward of 70 degrees
          ADD_Rn = 0.0d0

       !--------------------
       ! 70S-60S or 60N-70N 
       !--------------------
       ELSE IF ( LAT_L >= 60d0 ) THEN    

          IF ( LAT_H <= 70d0 ) THEN             

             ! If the entire grid box lies equatorward of 70 deg,
             ! then 222Rn emissions here are 0.005 [atoms/cm2/s]
             ADD_Rn = Rn_WATER
               
          ELSE
               
             ! Compute the fraction of the grid box below 70 degrees
             F_BELOW_70 = ( 70.0d0 - LAT_L ) / DENOM

             ! If the grid box straddles the 70S or 70N latitude line,
             ! then only count 222Rn emissions equatorward of 70 degrees.
             ! 222Rn emissions here are 0.005 [atoms/cm2/s].
             ADD_Rn = F_BELOW_70 * Rn_WATER
               
          ENDIF
            
       ELSE 

          !--------------------
          ! 70S-60S or 60N-70N
          !--------------------
          IF ( LAT_H > 60d0 ) THEN

             ! Fraction of grid box with ABS( lat ) below 60 degrees
             F_BELOW_60 = ( 60.0d0 - LAT_L ) / DENOM

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
       IF ( ExtState%TSURFK%Arr%Val(I,J) < 273.15 ) THEN
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
    CALL HCO_EmisAdd( HcoState, Arr2D, IDTRn, RC )
    Arr2D => NULL()
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Add diagnostics
    IF ( Diagn_AutoFillLevelDefined(2) ) THEN
       Arr2D => EmissRn(:,:)
       CALL Diagn_Update( am_I_Root,  HcoState,      ExtNr=ExtNr,  &
                          Cat=-1,     Hier=-1,       HcoID=IDTRn,  &
                          AutoFill=1, Array2D=Arr2D, RC=RC        )
       Arr2D => NULL()
       IF ( RC /= HCO_SUCCESS ) RETURN

    ENDIF

    !=======================================================================
    ! Compute 7Be emissions [kg/m2/s]    
    !
    ! Original units of 7Be emissions are [stars/g air/sec],
    ! where "stars" = # of nuclear disintegrations of cosmic rays
    !
    ! Now interpolate from 33 std levels onto GEOS-CHEM levels 
    !=======================================================================
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
       P_TMP   = ExtState%PCENTER%Arr%Val( I, J, L ) / 100.0_hp
                 
       ! Interpolate 7Be [stars/g air/sec] to GEOS-Chem levels
       CALL SLQ( LATSOU, PRESOU, BESOU, 10, 33, LAT_TMP, P_TMP, Be_TMP )

       ! Be_TMP = [stars/g air/s] * [0.045 atom/star] * 
       !          [kg air] * [1e3 g/kg] = 7Be emissions [atoms/s]
       Be_TMP  = Be_TMP * 0.045d0 * ExtState%AIR%Arr%Val(I,J,L) * 1.d3 
                  
       ! ADD_Be = [atoms/s] / [atom/kg] / [m2] = 7Be emissions [kg/m2/s]
       ADD_Be  = ( Be_TMP / XNUMOL_Be ) / HcoState%Grid%AREA_M2%Val(I,J)

#if defined( GCAP ) 
       !%%%%% FOR GCAP SIMULATION: Divide emissions flux by 3.5 to correct
       !%%%%% for the strat-trop exchange!  This replicates the prior code.
       !%%%%% (bmy, 9/4/14)
       IF ( .not. ( HcoState%Grid%PEDGE%Val(I,J,L) >          &
                    ExtState%TROPP%Arr%Val(I,J)     ) ) THEN
          ADD_Be = ADD_Be / 3.5d0
       ENDIF
#endif

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
    CALL HCO_EmisAdd( HcoState, Arr3D, IDTBe7, RC )
    Arr3D => NULL()
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Add diagnostics
    IF ( Diagn_AutoFillLevelDefined(2) ) THEN
       Arr3D => EmissBe7(:,:,:)
       CALL Diagn_Update( am_I_Root,  HcoState,      ExtNr=ExtNr,  &
                          Cat=-1,     Hier=-1,       HcoID=IDTBe7, &
                          AutoFill=1, Array3D=Arr3D, RC=RC        )
       Arr3D => NULL()
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Return w/ success
    CALL HCO_LEAVE ( RC )

  END SUBROUTINE HCOX_Gc_RnPbBe_Run
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
    LOGICAL                        :: verb
    CHARACTER(LEN=255)             :: MSG 

    ! Arrays
    INTEGER,           ALLOCATABLE :: HcoIDs(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)

    !=======================================================================
    ! HCOX_GC_RnPbBe_INIT begins here!
    !=======================================================================

    ! Get the extension number
    ExtNr = GetExtNr( TRIM( ExtName ) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter HEMCO
    CALL HCO_ENTER( 'HcoX_GC_RnPbBe_Init (hcox_gc_RnPbBe_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Verbose output?
    verb = ( am_I_Root .AND. HCO_VERBOSE_CHECK() )

    ! Set species IDs      
    CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Verbose mode
    IF ( verb ) THEN
       MSG = 'Use gc_RnPbBe emissions module (extension module)'
       CALL HCO_MSG( MSG )

       MSG = 'Use the following species (Name: HcoID):'
       CALL HCO_MSG(MSG)
       DO N = 1, nSpc
          WRITE(MSG,*) TRIM(SpcNames(N)), ':', HcoIDs(N)
          CALL HCO_MSG(MSG)
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

    ! ERROR: Rn tracer is not found!
    IF ( IDTRn <= 0 ) THEN
       CALL HCO_ERROR( 'Cannot find 222Rn tracer in list of species!', RC )
       RETURN
    ENDIF
    
    ! ERROR! Be7 tracer is not found
    IF ( IDTBe7 <= 0 ) THEN
       CALL HCO_ERROR( 'Cannot find 7Be tracer in list of species!', RC )
       RETURN
    ENDIF

    ! Activate met fields required by this extension
    ExtState%FRCLND%DoUse  = .TRUE. 
    ExtState%TSURFK%DoUse  = .TRUE. 
    ExtState%AIR%DoUse     = .TRUE. 
    ExtState%PCENTER%DoUse = .TRUE.
#if defined( GCAP ) 
    ExtState%TROPP%DoUse   = .TRUE.
#endif

    ! Activate this extension
    ExtState%Gc_RnPbBe     = .TRUE.

    !=======================================================================
    ! Initialize data arrays
    !=======================================================================

    ALLOCATE( EmissRn( HcoState%Nx, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( 'Cannot allocate EmissRn', RC )
       RETURN
    ENDIF 

    ALLOCATE( EmissBe7( HcoState%Nx, HcoState%NY, HcoState%NZ ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( 'Cannot allocate EmissBe7', RC )
       RETURN
    ENDIF 
    IF ( RC /= 0 ) RETURN

    ! Array for latitudes (Lal & Peters data)
    ALLOCATE( LATSOU( 10 ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( 'Cannot allocate LATSOU', RC )
       RETURN
    ENDIF 

    ! Array for pressures (Lal & Peters data)
    ALLOCATE( PRESOU( 33 ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( 'Cannot allocate PRESOU', RC )
       RETURN
    ENDIF 

    ! Array for 7Be emissions ( Lal & Peters data)
    ALLOCATE( BESOU( 10, 33 ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( 'Cannot allocate BESOU', RC )
       RETURN
    ENDIF 
    
    ! Initialize the 7Be emisisons data arrays
    CALL Init_7Be_Emissions()

    !=======================================================================
    ! Leave w/ success
    !=======================================================================
    IF ( ALLOCATED( HcoIDs   ) ) DEALLOCATE( HcoIDs   )
    IF ( ALLOCATED( SpcNames ) ) DEALLOCATE( SpcNames )

    CALL HCO_LEAVE ( RC ) 

  END SUBROUTINE HCOX_Gc_RnPbBe_Init
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_7Be_Emissions
!
! !DESCRIPTION: Subroutine Init\_7Be\_Emissions initializes the 7Be emissions 
!  from Lal \& Peters on 33 pressure levels.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_7Be_Emissions()
!
! !USES:
!
#include "hcox_gc_RnPbBe_include.H"
!
! !REMARKS:
!  You can update the 7Be emissions by including a new version of the
!  header file "hcox_gc_RnPbBe_include.H".
!
! !REVISION HISTORY: 
!  08 Aug 2014 - R. Yantosca - Now get code from an include file
!EOP
!------------------------------------------------------------------------------
!BOC
  END SUBROUTINE Init_7Be_Emissions
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
    REAL*8  :: X(N)     ! X-axis coordinate on original grid
    REAL*8  :: Y(M)     ! Y-axis coordinate on original grid
    REAL*8  :: Z(N,M)   ! Array of data on original grid
    REAL*8  :: U        ! X-axis coordinate for desired interpolated value
    REAL*8  :: V        ! Y-axis coordinate for desired interpolated value
!
! !OUTPUT PARAMETERS:
!
    REAL*8  :: W        ! Interpolated value of Z array, at coords (U,V) 
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*8  :: B(3), HH
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
