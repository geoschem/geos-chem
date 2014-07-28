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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Emissions indices etc.
  INTEGER                      :: ExtNr  = -1       ! HEMCO Extension number
  INTEGER                      :: IDTRn  = -1       ! Index # for Rn tracer
  INTEGER                      :: IDTBe7 = -1       ! Index # for 7Be tracer

  ! For tracking Rn222 and Be7 emissions
  REAL*8,  ALLOCATABLE, TARGET :: EmissRn (:,:  )
  REAL*8,  ALLOCATABLE, TARGET :: EmissBe7(:,:,:)

  ! For Lal & Peters 7Be emissions input data
  REAL*8,  ALLOCATABLE         :: LATSOU  (:    )   ! Array for latitudes
  REAL*8,  ALLOCATABLE         :: PRESOU  (:    )   ! Array for pressures
  REAL*8,  ALLOCATABLE         :: BESOU   (:,:  )   ! Array for 7Be emissions
!
! !DEFINED PARAMETERS:
!
  ! To convert kg to atoms
  REAL*8,  PARAMETER           :: XNUMOL_Rn = ( 6.0225d23 / 222.0d-3 )    
  REAL*8,  PARAMETER           :: XNUMOL_Be = ( 6.0225d23 /   7.0d-3 )

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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER           :: I,          J,          L,         N
    REAL*8            :: A_CM2,      ADD_Be,     ADD_Rn,    Rn_LAND
    REAL*8            :: Rn_WATER,   DTSRCE,     LAT_TMP,   P_TMP
    REAL*8            :: Be_TMP,     Rn_TMP,     LAT_S,     LAT_N
    REAL*8            :: LAT_H,      LAT_L,      F_LAND,    F_WATER
    REAL*8            :: F_BELOW_70, F_BELOW_60, F_ABOVE_60

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
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, LAT_S, LAT_N, LAT_H, LAT_L, F_BELOW_70     )
!!$OMP+PRIVATE( F_BELOW_60, F_ABOVE_60, A_CM2, Rn_LAND, Rn_WATER )
!!$OMP+PRIVATE( F_LAND, F_WATER, ADD_Rn                          )
    DO J = 1, HcoState%Ny
    DO I = 1, HcoState%Nx

       ! Get ABS( latitude ) at S and N edges of grid box
       LAT_S      = ABS( HcoState%Grid%YEDGE( I, J   ) ) 
       LAT_N      = ABS( HcoState%Grid%YEDGE( I, J+1 ) )
       LAT_H      = MAX( LAT_S, LAT_N )
       LAT_L      = MIN( LAT_S, LAT_N ) 
       
       ! Fraction of grid box w/ ABS( latitude ) less than 70 degrees
       F_BELOW_70 = ( 70.0d0 - LAT_L ) / ( LAT_H - LAT_L )

       ! Fraction of grid box w/ ABS( latitude ) less than 60 degrees
       F_BELOW_60 = ( 60.0d0 - LAT_L ) / ( LAT_H - LAT_L )

       ! Fraction of grid box w/ ABS( latitude ) greater than 60 degrees
       F_ABOVE_60 = 1d0 - F_BELOW_60

!-----------------------------------------------------------------------------
! Prior to 7/7/14:
! Now compute emissions in kg/m2/s for tracking in HEMCO (bmy, 7/7/14)
!       ! Baseline 222Rn emissions over land [kg]
!       ! Rn_LAND [kg] = [1 atom 222Rn/cm2/s] / [atoms/kg] * [s] * [cm2]
!       Rn_LAND    = 1d0 / XNUMOL_Rn * DTSRCE * A_CM2
!
!       ! Baseline 222Rn emissions over water or ice [kg]
!       Rn_WATER   = Rn_LAND * 0.005d0
!-----------------------------------------------------------------------------

       ! Baseline 222Rn emissions 
       ! Rn_LAND [kg/m2/s] = [1 atom 222Rn/cm2/s] / [atoms/kg] * [1d4 cm2/m2]
       Rn_LAND    = ( 1d0 / XNUMOL_Rn ) * 1d4

       ! Baseline 222Rn emissions over water or ice [kg]
       Rn_WATER   = Rn_LAND * 0.005d0

       ! Fraction of grid box that is land
       F_LAND     = ExtState%FRCLND%Arr%Val(I,J)

       ! Fraction of grid box that is water
       F_WATER    = 1d0 - F_LAND

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

             ADD_Rn =                                                    &
                        ! Consider 222Rn emissions equatorward of 
                        ! 60 degrees for both land (1.0 [atoms/cm2/s]) 
                        ! and water (0.005 [atoms/cm2/s])
                        F_BELOW_60 *                                     &
                        ( Rn_LAND  * F_LAND  ) +                         &
                        ( Rn_WATER * F_WATER ) +                         &

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
       IF ( ExtState%TSURFK%Arr%Val(I,J) < 273.15 ) ADD_Rn = ADD_Rn / 3d0

       ! Save 222Rn emissions into an array [kg/m2/s]
       EmissRn(I,J) = ADD_Rn
    ENDDO
    ENDDO
!!$OMP END PARALLEL DO

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
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L, LAT_TMP, P_TMP, Be_TMP, ADD_Be )
!!$OMP+SCHEDULE( DYNAMIC )
    DO L = 1, HcoState%Nz
    DO J = 1, HcoState%Ny
    DO I = 1, HcoState%Nx

       ! Get absolute value of latitude, since we will assume that 
       ! the 7Be distribution is symmetric about the equator
       LAT_TMP = ABS( HcoState%Grid%YMID( I, J ) )

       ! Pressure at (I,J,L)
       P_TMP   = ExtState%PCENTER%Arr%Val( I, J, L )
                 
       ! Interpolate 7Be [stars/g air/sec] to GEOS-CHEM levels
       CALL SLQ( LATSOU, PRESOU, BESOU, 10, 33, LAT_TMP, P_TMP, Be_TMP )

       ! Be_TMP = [stars/g air/s] * [0.045 atom/star] * 
       !          [kg air] * [1e3 g/kg] = 7Be emissions [atoms/s]
       Be_TMP  = Be_TMP * 0.045d0 * ExtState%AIR%Arr%Val(I,J,L) * 1.d3 
                  
!----------------------------------------------------------------------------
! Prior to 7/7/14:
! Convert Be7 emissions to [kg/m2/s] for tracking in HEMCO (bmy, 7/7/14)
!       ! ADD_Be = [atoms/s] * [s] / [atom/kg] = 7Be emissions [kg]
!       ADD_Be  = Be_TMP * DTSRCE / XNUMOL_Be 
!----------------------------------------------------------------------------

       ! ADD_Be = [atoms/s] / [atom/kg] / [m2] = 7Be emissions [kg/m2/s]
       ADD_Be  = ( Be_TMP / XNUMOL_Be ) / HcoState%Grid%AREA_M2(I,J)


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% TODO: This function is only really defined for GCAP, so we may need
!%%% to figure out how to test for ITS_IN_THE_STRATMESO here (bmy, 7/7/14)
!%%%       ! Correct the strat-trop exchange of 7Be
!%%%       IF ( ITS_IN_THE_STRATMESO( I, J, L, State_Met ) ) THEN
!%%%          CALL CORRECT_STE( ADD_Be )
!%%%       ENDIF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       ! Save emissions into an array for use below
       EmissBe7(I,J,L) = ADD_Be

    ENDDO
    ENDDO
    ENDDO
!!$OMP END PARALLEL DO

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
                          Cat=-1,     Hier=-1,       HcoID=IDTRn,  &
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
    USE HCO_ExtList_Mod,   ONLY : GetExtNr
    USE HCO_STATE_MOD,     ONLY : HCO_GetExtHcoID
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

    ! Set up tracer indices
    DO N = 1, nSpc
       SELECT CASE( TRIM( SpcNames(N) ) )
          CASE( 'Rn', 'Rn222', '222Rn' )
             IDTRn = N
          CASE( 'Be', 'Be7', '7Be' )
             IDTBe7 = N
          CASE DEFAULT
             ! Do nothing
       END SELECT
    ENDDO

    ! ERROR: Rn tracer is not found!
    IF ( IDTRn <= 0 ) THEN
       RC = HCO_FAIL
       CALL HCO_ERROR( 'Cannot find 222Rn tracer in list of species!', RC )
       RETURN
    ENDIF
    
    ! ERROR! Be7 tracer is not found
    IF ( IDTBe7 <= 0 ) THEN
       RC = HCO_FAIL
       CALL HCO_ERROR( 'Cannot find 7Be tracer in list of species!', RC )
       RETURN
    ENDIF

    ! Activate met fields required by this extension
    ExtState%FRCLND%DoUse  = .TRUE. 
    ExtState%TSURFK%DoUse  = .TRUE. 
    ExtState%AIR%DoUse     = .TRUE. 
    ExtState%PCENTER%DoUse = .TRUE.

    ! Activate this extension
    ExtState%Gc_RnPbBe    = .TRUE.

    !=======================================================================
    ! Initialize data arrays
    !=======================================================================

    ALLOCATE( EmissRn( HcoState%Nx, HcoState%NY ), STAT=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ALLOCATE( EmissBe7( HcoState%Nx, HcoState%NY, HcoState%NZ ), STAT=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Array for latitudes (Lal & Peters data)
    ALLOCATE( LATSOU( 10 ), STAT=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Array for pressures (Lal & Peters data)
    ALLOCATE( PRESOU( 33 ), STAT=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Array for 7Be emissions ( Lal & Peters data)
    ALLOCATE( BESOU( 10, 33 ), STAT=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    
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
! !REMARKS:
!  (1) In prior versions of GEOS-Chem, this routine was named READ_7BE, and
!      it read the ASCII file "7Be.Lal".   Because this data set is not placed
!      on a lat/lon grid, ESMF cannot regrid it.  To work around this, we now
!      hardwire this data in module arrays rather than read it from disk.
!                                                                             .
!  (2) Units of 7Be emissions are [stars/g air/s].  
!      Here, "stars" = # of nuclear disintegrations of cosmic rays
!                                                                             .
!  (3) Original data from Lal & Peters (1967), w/ these modifications:
!      (a) Replace data at (0hPa, 70S) following Koch 1996:
!          (i ) old value = 3000 
!          (ii) new value = 1900
!      (b) Copy data from 70S to 80S and 90S at all levels
!                                                                             .
!  (4) This new subroutine yielded identical results to the prior code 
!      in a difference test done by Bob Yantosca on 07 Jul 2014.
!
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
!EOP
!------------------------------------------------------------------------------
!BOC

    !=======================================================================
    ! READ_7BE begins here!
    !=======================================================================
    
    ! Define latitudes [degrees North]
    LATSOU      = (/     0d0,     10d0,     20d0,    30d0,    40d0,   &
                        50d0,     60d0,     70d0,    80d0,    90d0  /)
 
    ! Define pressures [hPa]
    PRESOU      = (/     0d0,     50d0,     70d0,    90d0,   110d0,   &
                       130d0,    150d0,    170d0,   190d0,   210d0,   &
                       230d0,    250d0,    270d0,   290d0,   313d0,   &
                       338d0,    364d0,    392d0,   420d0,   451d0,   &
                       485d0,    518d0,    555d0,   592d0,   633d0,   &
                       680d0,    725d0,    772d0,   822d0,   875d0,   &
                       930d0,    985d0,   1030d0                    /)

    ! Define 7Be emissions [stars/g air/s]
    ! 1 "star" = 1 nuclear disintegration via cosmic rays
    !
    ! NOTE: These statements were defined from printout of the file
    ! and need to be multiplied by 1d-5 below.
    BESOU(:,1)  = (/   150d0,    156d0,    188d0,   285d0,   500d0,   &
                       910d0,   1700d0,   1900d0,  1900d0,  1900d0  /)
                       
    BESOU(:,2)  = (/   280d0,    310d0,    390d0,   590d0,   880d0,   &
                      1390d0,   1800d0,   1800d0,  1800d0,  1800d0  /)
                       
    BESOU(:,3)  = (/   310d0,    330d0,    400d0,   620d0,   880d0,   &
                      1280d0,   1450d0,   1450d0,  1450d0,  1450d0  /)
                       
    BESOU(:,4)  = (/   285d0,    310d0,    375d0,   570d0,   780d0,   &
                      1100d0,   1180d0,   1180d0,  1180d0,  1180d0  /)
                       
    BESOU(:,5)  = (/   255d0,    275d0,    330d0,   510d0,   680d0,   &
                       950d0,   1000d0,   1000d0,  1000d0,  1000d0  /)
                       
    BESOU(:,6)  = (/   230d0,    245d0,    292d0,   450d0,   600d0,   &
                       820d0,    875d0,    875d0,   875d0,   875d0  /)
                       
    BESOU(:,7)  = (/   205d0,    215d0,    260d0,   400d0,   530d0,   &
                       730d0,    750d0,    750d0,   750d0,   750d0  /)
                       
    BESOU(:,8)  = (/   182d0,    195d0,    235d0,   355d0,   480d0,   &
                       630d0,    650d0,    650d0,   650d0,   650d0  /)
                       
    BESOU(:,9)  = (/   160d0,    173d0,    208d0,   315d0,   410d0,   &
                       543d0,    550d0,    550d0,   550d0,   550d0  /)
                       
    BESOU(:,10) = (/   148d0,    152d0,    185d0,   280d0,   370d0,   &
                       480d0,    500d0,    500d0,   500d0,   500d0  /)
                       
    BESOU(:,11) = (/   130d0,    139d0,    167d0,   250d0,   320d0,   &
                       425d0,    430d0,    430d0,   430d0,   430d0  /)
                       
    BESOU(:,12) = (/   116d0,    123d0,    148d0,   215d0,   285d0,   &
                       365d0,    375d0,    375d0,   375d0,   375d0  /)
                       
    BESOU(:,13) = (/   104d0,    110d0,    130d0,   198d0,   250d0,   &
                       320d0,    330d0,    330d0,   330d0,   330d0  /)
                       
    BESOU(:,14) = (/    93d0,     99d0,    118d0,   170d0,   222d0,   &
                       280d0,    288d0,    288d0,   288d0,   288d0  /)
                       
    BESOU(:,15) = (/    80d0,     84d0,    100d0,   145d0,   190d0,   &
                       235d0,    250d0,    250d0,   250d0,   250d0  /)
                       
    BESOU(:,16) = (/    72d0,   74.5d0,     88d0,   129d0,   168d0,   &
                       210d0,    218d0,    218d0,   218d0,   218d0  /)
                       
    BESOU(:,17) = (/  59.5d0,   62.5d0,   73.5d0,   108d0,   138d0,   &
                       171d0,    178d0,    178d0,   178d0,   178d0  /)
                       
    BESOU(:,18) = (/    50d0,     53d0,     64d0,    90d0,   115d0,   &
                       148d0,    150d0,    150d0,   150d0,   150d0  /)
                       
    BESOU(:,19) = (/    45d0,   46.5d0,   52.5d0,    76d0,    98d0,   &
                       122d0,    128d0,    128d0,   128d0,   128d0  /)
                       
    BESOU(:,20) = (/  36.5d0,   37.5d0,     45d0,    61d0,    77d0,   &
                        98d0,    102d0,    102d0,   102d0,   102d0  /)
                       
    BESOU(:,21) = (/  30.8d0,     32d0,   37.5d0,  51.5d0,    65d0,   &
                        81d0,     85d0,     85d0,    85d0,    85d0  /)

    BESOU(:,22) = (/  25.5d0,   26.5d0,     32d0,  40.5d0,    54d0,   &
                      67.5d0,   69.5d0,   69.5d0,  69.5d0,  69.5d0  /)

    BESOU(:,23) = (/  20.5d0,   21.6d0,   25.5d0,    33d0,    42d0,   &
                      53.5d0,     55d0,     55d0,    55d0,    55d0  /)

    BESOU(:,24) = (/  16.8d0,   17.3d0,     20d0,    26d0,  33.5d0,   &
                        41d0,     43d0,     43d0,    43d0,    43d0  /)

    BESOU(:,25) = (/    13d0,   13.8d0,   15.3d0,  20.5d0,  26.8d0,   &
              &       32.5d0,   33.5d0,   33.5d0,  33.5d0,  33.5d0  /)

    BESOU(:,26) = (/  10.1d0,   10.6d0,   12.6d0,  15.8d0,    20d0,   &
                      24.5d0,   25.8d0,   25.8d0,  25.8d0,  25.8d0  /)

    BESOU(:,27) = (/   7.7d0,   8.15d0,    9.4d0,  11.6d0,  14.8d0,   &
                      17.8d0,   18.5d0,   18.5d0,  18.5d0,  18.5d0  /)
 
    BESOU(:,28) = (/   5.7d0,   5.85d0,   6.85d0,  8.22d0,    11d0,   &
                      13.1d0,   13.2d0,   13.2d0,  13.2d0,  13.2d0  /)

    BESOU(:,29) = (/   3.9d0,    4.2d0,   4.85d0,     6d0,   7.6d0,   &
              &          9d0,    9.2d0,    9.2d0,   9.2d0,   9.2d0  /)

    BESOU(:,30) = (/     3d0,   3.05d0,   3.35d0,   4.2d0,   5.3d0,   &
                       5.9d0,   6.25d0,   6.25d0,  6.25d0,  6.25d0  /)

    BESOU(:,31) = (/  2.05d0,    2.1d0,   2.32d0,   2.9d0,   3.4d0,   &
                       3.9d0,    4.1d0,    4.1d0,   4.1d0,   4.1d0  /)

    BESOU(:,32) = (/  1.45d0,   1.43d0,   1.65d0,  2.03d0,   2.4d0,   &
                      2.75d0,   2.65d0,   2.65d0,  2.65d0,  2.65d0  /)

    BESOU(:,33) = (/  1.04d0,   1.08d0,   1.21d0,   1.5d0,  1.68d0,   &
                       1.8d0,    1.8d0,    1.8d0,   1.8d0,   1.8d0  /)

    ! All the numbers of BESOU need to be multiplied by 1e-5 in order to put 
    ! them into the correct data range.  NOTE: This multiplication statement
    ! needs to be preserved here in order to  ensure identical output to the
    ! prior code! (bmy, 7/7/14)
    BESOU = BESOU * 1d-5

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
