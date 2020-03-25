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
!  (5 ) Koch and Rind, Beryllium 10/beryllium 7 as a tracer of stratospheric
!        transport, JGR, 103, D4, 3907-3917, 1998.
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
!  07 Dec 2018 - M. Sulprizio- Add Be10, Be7Strat, and Be10Strat emissions
!  25 Jan 2019 - M. Sulprizio- Add instance wrapper
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  TYPE :: MyInst

   ! Emissions indices etc.
   INTEGER               :: Instance
   INTEGER               :: ExtNr         ! HEMCO Extension number
   INTEGER               :: IDTRn222      ! Index # for Rn222
   INTEGER               :: IDTBe7        ! Index # for Be7
   INTEGER               :: IDTBe7Strat   ! Index # for Be7Strat
   INTEGER               :: IDTBe10       ! Index # for Be10
   INTEGER               :: IDTBe10Strat  ! Index # for Be10Strat

   ! For tracking Rn222, Be7, and Be10 emissions
   REAL(hp), POINTER     :: EmissRn222    (:,:  )
   REAL(hp), POINTER     :: EmissBe7      (:,:,:)
   REAL(hp), POINTER     :: EmissBe7Strat (:,:,:)
   REAL(hp), POINTER     :: EmissBe10     (:,:,:)
   REAL(hp), POINTER     :: EmissBe10Strat(:,:,:)

   ! For Lal & Peters 7Be emissions input data
   REAL(hp), POINTER     :: LATSOU(:    ) ! Array for latitudes
   REAL(hp), POINTER     :: PRESOU(:    ) ! Array for pressures
   REAL(hp), POINTER     :: BESOU (:,:  ) ! Array for 7Be emissions

   TYPE(MyInst), POINTER :: NextInst => NULL()
  END TYPE MyInst

  ! Pointer to instances
  TYPE(MyInst), POINTER  :: AllInst => NULL()
!
! !DEFINED PARAMETERS:
!
  ! To convert kg to atoms
  REAL*8,  PARAMETER     :: XNUMOL_Rn   = ( 6.022140857d23 / 222.0d-3 )
  REAL*8,  PARAMETER     :: XNUMOL_Be7  = ( 6.022140857d23 /   7.0d-3 )
  REAL*8,  PARAMETER     :: XNUMOL_Be10 = ( 6.022140857d23 /  10.0d-3 )

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Gc_RnPbBe_run
!
! !DESCRIPTION: Subroutine HcoX\_Gc\_RnPbBe\_Run computes emissions of 222Rn,
!  7Be, and 10Be for the GEOS-Chem Rn-Pb-Be specialty simulation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Gc_RnPbBe_Run( ExtState, HcoState, RC )
!
! !USES:
!
    ! HEMCO modules
    USE HCO_FluxArr_Mod, ONLY : HCO_EmisAdd
!
! !INPUT PARAMETERS:
!
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
    REAL*8            :: A_CM2,    ADD_Rn,     Add_Be7,    Add_Be10
    REAL*8            :: Rn_LAND,  Rn_WATER,   DTSRCE
    REAL*8            :: Rn_TMP,   LAT,        F_LAND
    REAL*8            :: F_WATER,  F_BELOW_70, F_BELOW_60, F_ABOVE_60
    REAL*8            :: DENOM
    REAL(hp)          :: LAT_TMP,  P_TMP,      Be_TMP
    CHARACTER(LEN=255):: MSG

    ! Pointers
    TYPE(MyInst), POINTER :: Inst
    REAL(hp),     POINTER :: Arr2D(:,:  )
    REAL(hp),     POINTER :: Arr3D(:,:,:)

    !=======================================================================
    ! HCOX_GC_RnPbBe_RUN begins here!
    !=======================================================================

    ! Return if extension not turned on
    IF ( ExtState%GC_RnPbBe <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_GC_RnPbBe_Run (hcox_gc_RnPbBe_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Set error flag
    !ERR = .FALSE.

    ! Get instance
    Inst   => NULL()
    CALL InstGet ( ExtState%GC_RnPbBe, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(MSG,*) 'Cannot find GC_RnPbBe instance Nr. ', ExtState%GC_RnPbBe
       CALL HCO_ERROR(HcoState%Config%Err,MSG,RC)
       RETURN
    ENDIF

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
    IF ( Inst%IDTRn222 > 0 ) THEN

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
          Inst%EmissRn222(I,J) = ADD_Rn
       ENDDO
       ENDDO
!$OMP END PARALLEL DO

       !------------------------------------------------------------------------
       ! Add 222Rn emissions to HEMCO data structure & diagnostics
       !------------------------------------------------------------------------

       ! Add emissions
       Arr2D => Inst%EmissRn222(:,:)
       CALL HCO_EmisAdd( HcoState, Arr2D, Inst%IDTRn222, &
                         RC,       ExtNr=Inst%ExtNr )
       Arr2D => NULL()
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'HCO_EmisAdd error: EmissRn222', RC )
          RETURN
       ENDIF

    ENDIF ! IDTRn222 > 0

    !=======================================================================
    ! Compute 7Be and 10Be emissions [kg/m2/s]
    !
    ! Original units of 7Be and 10Be emissions are [stars/g air/sec],
    ! where "stars" = # of nuclear disintegrations of cosmic rays
    !
    ! Now interpolate from 33 std levels onto GEOS-CHEM levels
    !
    ! 7Be and 10Be have identical source distributions (Koch and Rind, 1998)
    !=======================================================================
    IF ( Inst%IDTBe7 > 0 .or. Inst%IDTBe10 > 0 ) THEN
!$OMP PARALLEL DO                                                   &
!$OMP DEFAULT( SHARED )                                             &
!$OMP PRIVATE( I, J, L, LAT_TMP, P_TMP, Be_TMP, ADD_Be7, ADD_Be10 ) &
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
          CALL SLQ( Inst%LATSOU, Inst%PRESOU, Inst%BESOU, 10, 33, &
                    LAT_TMP,     P_TMP,       Be_TMP )

          ! Be_TMP = [stars/g air/s] * [0.045 atom/star] *
          !          [kg air] * [1e3 g/kg] = 7Be/10Be emissions [atoms/s]
          Be_TMP  = Be_TMP * 0.045e+0_hp * ExtState%AIR%Arr%Val(I,J,L) * 1.e+3_hp

          ! ADD_Be = [atoms/s] / [atom/kg] / [m2] = 7Be/10Be emissions [kg/m2/s]
          ADD_Be7  = ( Be_TMP / XNUMOL_Be7  ) / HcoState%Grid%AREA_M2%Val(I,J)
          ADD_Be10 = ( Be_TMP / XNUMOL_Be10 ) / HcoState%Grid%AREA_M2%Val(I,J)

          ! Save emissions into an array for use below
          Inst%EmissBe7 (I,J,L) = ADD_Be7
          Inst%EmissBe10(I,J,L) = ADD_Be10
          IF ( L > ExtState%TropLev%Arr%Val(I,J) ) THEN
             IF ( Inst%IDTBe7Strat > 0 ) THEN
                Inst%EmissBe7Strat (I,J,L) = Add_Be7
             ENDIF
             IF ( Inst%IDTBe10Strat > 0 ) THEN 
                Inst%EmissBe10Strat(I,J,L) = Add_Be10
             ENDIF
          ENDIF

       ENDDO
       ENDDO
       ENDDO
!$OMP END PARALLEL DO

       !------------------------------------------------------------------------
       ! Add Be7 and Be10 emissions to HEMCO data structure & diagnostics
       !------------------------------------------------------------------------

       ! Add emissions
       IF ( Inst%IDTBe7 > 0 ) THEN
          Arr3D => Inst%EmissBe7(:,:,:)
          CALL HCO_EmisAdd( HcoState, Arr3D, Inst%IDTBe7, &
                            RC,       ExtNr=Inst%ExtNr )
          Arr3D => NULL()
          IF ( RC /= HCO_SUCCESS ) THEN
             CALL HCO_ERROR( HcoState%Config%Err, &
                             'HCO_EmisAdd error: EmissBe7', RC )
             RETURN
          ENDIF
       ENDIF

       ! Add emissions
       IF ( Inst%IDTBe7Strat > 0 ) THEN
          Arr3D => Inst%EmissBe7Strat(:,:,:)
          CALL HCO_EmisAdd( HcoState, Arr3D, Inst%IDTBe7Strat, &
                            RC,       ExtNr=Inst%ExtNr )
          Arr3D => NULL()
          IF ( RC /= HCO_SUCCESS ) THEN
             CALL HCO_ERROR( HcoState%Config%Err, &
                             'HCO_EmisAdd error: EmissBe7Strat', RC )
             RETURN
          ENDIF
       ENDIF

       ! Add emissions
       IF ( Inst%IDTBe10 > 0 ) THEN
          Arr3D => Inst%EmissBe10(:,:,:)
          CALL HCO_EmisAdd( HcoState, Arr3D, Inst%IDTBe10, &
                            RC,       ExtNr=Inst%ExtNr )
          Arr3D => NULL()
          IF ( RC /= HCO_SUCCESS ) THEN
             CALL HCO_ERROR( HcoState%Config%Err, &
                             'HCO_EmisAdd error: EmissBe10', RC )
             RETURN
          ENDIF
       ENDIF

       ! Add emissions
       IF ( Inst%IDTBe10Strat > 0 ) THEN
          Arr3D => Inst%EmissBe10Strat(:,:,:)
          CALL HCO_EmisAdd( HcoState, Arr3D, Inst%IDTBe10Strat, &
                            RC,       ExtNr=Inst%ExtNr )
          Arr3D => NULL()
          IF ( RC /= HCO_SUCCESS ) THEN
             CALL HCO_ERROR( HcoState%Config%Err, &
                             'HCO_EmisAdd error: EmissBe10Strat', RC )
             RETURN
          ENDIF
       ENDIF

    ENDIF !IDTBe7 > 0 or IDTBe10 > 0

    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Nullify pointers
    Inst    => NULL()

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
  SUBROUTINE HCOX_Gc_RnPbBe_Init( HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_ExtList_Mod, ONLY : GetExtNr
    USE HCO_State_Mod,   ONLY : HCO_GetExtHcoID
!
! !INPUT PARAMETERS:
!
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
    INTEGER                        :: N, nSpc, ExtNr
    CHARACTER(LEN=255)             :: MSG

    ! Arrays
    INTEGER,           ALLOCATABLE :: HcoIDs(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)

    ! Pointers
    TYPE(MyInst), POINTER          :: Inst

    !=======================================================================
    ! HCOX_GC_RnPbBe_INIT begins here!
    !=======================================================================

    ! Get the extension number
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'HcoX_GC_RnPbBe_Init (hcox_gc_RnPbBe_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Create Instance
    Inst => NULL()
    CALL InstCreate ( ExtNr, ExtState%GC_RnPbBe, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot create GC_RnPbBe instance', RC )
       RETURN
    ENDIF
    ! Also fill Inst%ExtNr
    Inst%ExtNr = ExtNr

    ! Set HEMCO species IDs
    CALL HCO_GetExtHcoID( HcoState, Inst%ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Verbose mode
    IF ( HcoState%amIRoot ) THEN
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
             Inst%IDTRn222     = HcoIDs(N)
          CASE( 'Be', 'Be7', '7Be' )
             Inst%IDTBe7       = HcoIDs(N)
          CASE( 'Be7Strat', '7BeStrat' )
             Inst%IDTBe7Strat  = HcoIDs(N)
          CASE( 'Be10', '10Be' )
             Inst%IDTBe10      = HcoIDs(N)
          CASE( 'Be10Strat', '10BeStrat' )
             Inst%IDTBe10Strat = HcoIDs(N)
          CASE DEFAULT
             ! Do nothing
       END SELECT
    ENDDO

    ! WARNING: Rn tracer is not found!
    IF ( Inst%IDTRn222 <= 0 .AND. HcoState%amIRoot ) THEN
       CALL HCO_WARNING( HcoState%Config%Err, &
                         'Cannot find Rn222 tracer in list of species!', RC )
    ENDIF

    ! WARNING: Be7 tracer is not found
    IF ( Inst%IDTBe7 <= 0 .AND. HcoState%amIRoot ) THEN
       CALL HCO_WARNING( HcoState%Config%Err, &
                         'Cannot find Be7 tracer in list of species!', RC )
    ENDIF

    ! WARNING: Be10 tracer is not found
    IF ( Inst%IDTBe10 <= 0 .AND. HcoState%amIRoot ) THEN
       CALL HCO_WARNING( HcoState%Config%Err, &
                        'Cannot find Be10 tracer in list of species!', RC )
    ENDIF

    ! ERROR: No tracer defined
    IF ( Inst%IDTRn222 <= 0 .AND. Inst%IDTBe7 <= 0 .AND. Inst%IDTBe10 <= 0) THEN
       CALL HCO_ERROR( HcoState%Config%Err, &
                       'Cannot use RnPbBe extension: no valid species!', RC )
    ENDIF

    ! Activate met fields required by this extension
    ExtState%FRCLND%DoUse  = .TRUE.
    ExtState%T2M%DoUse     = .TRUE.
    ExtState%AIR%DoUse     = .TRUE.
    ExtState%TropLev%DoUse = .TRUE.

    !=======================================================================
    ! Initialize data arrays
    !=======================================================================

    IF ( Inst%IDTRn222 > 0 ) THEN
       ALLOCATE( Inst%EmissRn222( HcoState%Nx, HcoState%NY ), STAT=RC )
       IF ( RC /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, &
                           'Cannot allocate EmissRn222', RC )
          RETURN
       ENDIF
    ENDIF

    IF ( Inst%IDTBe7 > 0 ) THEN
       ALLOCATE( Inst%EmissBe7( HcoState%Nx, HcoState%NY, HcoState%NZ ), &
                 STAT=RC )
       IF ( RC /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, &
                           'Cannot allocate EmissBe7', RC )
          RETURN
       ENDIF
       IF ( RC /= 0 ) RETURN

       ! Array for latitudes (Lal & Peters data)
       ALLOCATE( Inst%LATSOU( 10 ), STAT=RC )
       IF ( RC /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, &
                           'Cannot allocate LATSOU', RC )
          RETURN
       ENDIF

       ! Array for pressures (Lal & Peters data)
       ALLOCATE( Inst%PRESOU( 33 ), STAT=RC )
       IF ( RC /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, &
                           'Cannot allocate PRESOU', RC )
          RETURN
       ENDIF

       ! Array for 7Be emissions ( Lal & Peters data)
       ALLOCATE( Inst%BESOU( 10, 33 ), STAT=RC )
       IF ( RC /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, &
                           'Cannot allocate BESOU', RC )
          RETURN
       ENDIF

       ! Initialize the 7Be emisisons data arrays
       CALL Init_7Be_Emissions( Inst )
    ENDIF

    IF ( Inst%IDTBe7Strat > 0 ) THEN
       ALLOCATE( Inst%EmissBe7Strat( HcoState%Nx, HcoState%NY, HcoState%NZ ), &
                 STAT=RC )
       IF ( RC /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, &
                           'Cannot allocate EmissBe7Strat', RC )
          RETURN
       ENDIF
       IF ( RC /= 0 ) RETURN
    ENDIF

    IF ( Inst%IDTBe10 > 0 ) THEN
       ALLOCATE( Inst%EmissBe10( HcoState%Nx, HcoState%NY, HcoState%NZ ), &
                 STAT=RC )
       IF ( RC /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, &
                           'Cannot allocate EmissBe10', RC )
          RETURN
       ENDIF
       IF ( RC /= 0 ) RETURN
    ENDIF

    IF ( Inst%IDTBe10Strat > 0 ) THEN
       ALLOCATE( Inst%EmissBe10Strat( HcoState%Nx, HcoState%NY, HcoState%NZ ), &
                 STAT=RC )
       IF ( RC /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, &
                           'Cannot allocate EmissBe10Strat', RC )
          RETURN
       ENDIF
       IF ( RC /= 0 ) RETURN
    ENDIF

    !=======================================================================
    ! Leave w/ success
    !=======================================================================
    IF ( ALLOCATED( HcoIDs   ) ) DEALLOCATE( HcoIDs   )
    IF ( ALLOCATED( SpcNames ) ) DEALLOCATE( SpcNames )

    ! Nullify pointers
    Inst    => NULL()

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
  SUBROUTINE HCOX_Gc_RnPbBe_Final( ExtState )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options
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

    CALL InstRemove ( ExtState%GC_RnPbBe )

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
  SUBROUTINE Init_7Be_Emissions( Inst )
!
! !INPUT PARAMETERS:
!
    TYPE(MyInst),    POINTER        :: Inst      ! Instance
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
    Inst%LATSOU      = (/     0.0_hp,     10.0_hp,     20.0_hp,    30.0_hp,  &
                             40.0_hp,     50.0_hp,     60.0_hp,    70.0_hp,  &
                             80.0_hp,     90.0_hp  /)

    ! Define pressures [hPa]
    Inst%PRESOU      = (/     0.0_hp,     50.0_hp,     70.0_hp,    90.0_hp,  &
                            110.0_hp,    130.0_hp,    150.0_hp,   170.0_hp,  &
                            190.0_hp,    210.0_hp,    230.0_hp,   250.0_hp,  &
                            270.0_hp,    290.0_hp,    313.0_hp,   338.0_hp,  &
                            364.0_hp,    392.0_hp,    420.0_hp,   451.0_hp,  &
                            485.0_hp,    518.0_hp,    555.0_hp,   592.0_hp,  &
                            633.0_hp,    680.0_hp,    725.0_hp,   772.0_hp,  &
                            822.0_hp,    875.0_hp,    930.0_hp,   985.0_hp,  &
                           1030.0_hp  /)

    ! Define 7Be emissions [stars/g air/s]
    ! 1 "star" = 1 nuclear disintegration via cosmic rays
    !
    ! NOTE: These statements were defined from printout of the file
    ! and need to be multiplied by 1d-5 below.
    Inst%BESOU(:,1)  = (/   150.0_hp,    156.0_hp,    188.0_hp,   285.0_hp,  &
                            500.0_hp,    910.0_hp,   1700.0_hp,  1900.0_hp,  &
                           1900.0_hp,   1900.0_hp  /)

    Inst%BESOU(:,2)  = (/   280.0_hp,    310.0_hp,    390.0_hp,   590.0_hp,  &
                            880.0_hp,   1390.0_hp,   1800.0_hp,  1800.0_hp,  &
                           1800.0_hp,   1800.0_hp  /)

    Inst%BESOU(:,3)  = (/   310.0_hp,    330.0_hp,    400.0_hp,   620.0_hp,  &
                            880.0_hp,   1280.0_hp,   1450.0_hp,  1450.0_hp,  &
                           1450.0_hp,   1450.0_hp  /)

    Inst%BESOU(:,4)  = (/   285.0_hp,    310.0_hp,    375.0_hp,   570.0_hp,  &
                            780.0_hp,   1100.0_hp,   1180.0_hp,  1180.0_hp,  &
                           1180.0_hp,   1180.0_hp  /)

    Inst%BESOU(:,5)  = (/   255.0_hp,    275.0_hp,    330.0_hp,   510.0_hp,  &
                            680.0_hp,    950.0_hp,   1000.0_hp,  1000.0_hp,  &
                           1000.0_hp,   1000.0_hp  /)

    Inst%BESOU(:,6)  = (/   230.0_hp,    245.0_hp,    292.0_hp,   450.0_hp,  &
                            600.0_hp,    820.0_hp,    875.0_hp,   875.0_hp,  &
                            875.0_hp,    875.0_hp  /)

    Inst%BESOU(:,7)  = (/   205.0_hp,    215.0_hp,    260.0_hp,   400.0_hp,  &
                            530.0_hp,    730.0_hp,    750.0_hp,   750.0_hp,  &
                            750.0_hp,    750.0_hp  /)

    Inst%BESOU(:,8)  = (/   182.0_hp,    195.0_hp,    235.0_hp,   355.0_hp,  &
                            480.0_hp,    630.0_hp,    650.0_hp,   650.0_hp,  &
                            650.0_hp,    650.0_hp  /)

    Inst%BESOU(:,9)  = (/   160.0_hp,    173.0_hp,    208.0_hp,   315.0_hp,  &
                            410.0_hp,    543.0_hp,    550.0_hp,   550.0_hp,  &
                            550.0_hp,    550.0_hp  /)

    Inst%BESOU(:,10) = (/   148.0_hp,    152.0_hp,    185.0_hp,   280.0_hp,  &
                            370.0_hp,    480.0_hp,    500.0_hp,   500.0_hp,  &
                            500.0_hp,    500.0_hp  /)

    Inst%BESOU(:,11) = (/   130.0_hp,    139.0_hp,    167.0_hp,   250.0_hp,  &
                            320.0_hp,    425.0_hp,    430.0_hp,   430.0_hp,  &
                            430.0_hp,    430.0_hp  /)

    Inst%BESOU(:,12) = (/   116.0_hp,    123.0_hp,    148.0_hp,   215.0_hp,  &
                            285.0_hp,    365.0_hp,    375.0_hp,   375.0_hp,  &
                            375.0_hp,    375.0_hp  /)

    Inst%BESOU(:,13) = (/   104.0_hp,    110.0_hp,    130.0_hp,   198.0_hp,  &
                            250.0_hp,    320.0_hp,    330.0_hp,   330.0_hp,  &
                            330.0_hp,    330.0_hp  /)

    Inst%BESOU(:,14) = (/    93.0_hp,     99.0_hp,    118.0_hp,   170.0_hp,  &
                            222.0_hp,    280.0_hp,    288.0_hp,   288.0_hp,  &
                            288.0_hp,    288.0_hp  /)

    Inst%BESOU(:,15) = (/    80.0_hp,     84.0_hp,    100.0_hp,   145.0_hp,  &
                            190.0_hp,    235.0_hp,    250.0_hp,   250.0_hp,  &
                            250.0_hp,    250.0_hp  /)

    Inst%BESOU(:,16) = (/    72.0_hp,     74.0_hp,     88.0_hp,   129.0_hp,  &
                            168.0_hp,    210.0_hp,    218.0_hp,   218.0_hp,  &
                            218.0_hp,    218.0_hp  /)

    Inst%BESOU(:,17) = (/    59.5_hp,     62.5_hp,     73.5_hp,   108.0_hp,  &
                            138.0_hp,    171.0_hp,    178.0_hp,   178.0_hp,  &
                            178.0_hp,    178.0_hp  /)

    Inst%BESOU(:,18) = (/    50.0_hp,     53.0_hp,     64.0_hp,    90.0_hp,  &
                            115.0_hp,    148.0_hp,    150.0_hp,   150.0_hp,  &
                            150.0_hp,    150.0_hp  /)

    Inst%BESOU(:,19) = (/    45.0_hp,     46.5_hp,     52.5_hp,    76.0_hp,  &
                             98.0_hp,    122.0_hp,    128.0_hp,   128.0_hp,  &
                            128.0_hp,    128.0_hp  /)

    Inst%BESOU(:,20) = (/    36.5_hp,     37.5_hp,     45.0_hp,    61.0_hp,  &
                             77.0_hp,     98.0_hp,    102.0_hp,   102.0_hp,  &
                            102.0_hp,    102.0_hp  /)

    Inst%BESOU(:,21) = (/    30.8_hp,     32.0_hp,     37.5_hp,    51.5_hp,  &
                             65.0_hp,     81.0_hp,     85.0_hp,    85.0_hp,  &
                             85.0_hp,     85.0_hp  /)

    Inst%BESOU(:,22) = (/    25.5_hp,     26.5_hp,     32.0_hp,    40.5_hp,  &
                             54.0_hp,     67.5_hp,     69.5_hp,    69.5_hp,  &
                             69.5_hp,     69.5_hp  /)

    Inst%BESOU(:,23) = (/    20.5_hp,     21.6_hp,     25.5_hp,    33.0_hp,  &
                             42.0_hp,     53.5_hp,     55.0_hp,    55.0_hp,  &
                             55.0_hp,     55.0_hp  /)

    Inst%BESOU(:,24) = (/    16.8_hp,     17.3_hp,     20.0_hp,    26.0_hp,  &
                             33.5_hp,     41.0_hp,     43.0_hp,    43.0_hp,  &
                             43.0_hp,     43.0_hp  /)

    Inst%BESOU(:,25) = (/    13.0_hp,     13.8_hp,     15.3_hp,    20.5_hp,  &
                             26.8_hp,     32.5_hp,     33.5_hp,    33.5_hp,  &
                             33.5_hp,     33.5_hp  /)

    Inst%BESOU(:,26) = (/    10.1_hp,     10.6_hp,     12.6_hp,    15.8_hp,  &
                             20.0_hp,     24.5_hp,     25.8_hp,    25.8_hp,  &
                             25.8_hp,     25.8_hp  /)

    Inst%BESOU(:,27) = (/     7.7_hp,     8.15_hp,      9.4_hp,    11.6_hp,  &
                             14.8_hp,     17.8_hp,     18.5_hp,    18.5_hp,  &
                             18.5_hp,     18.5_hp  /)

    Inst%BESOU(:,28) = (/     5.7_hp,     5.85_hp,     6.85_hp,    8.22_hp,  &
                             11.0_hp,     13.1_hp,     13.2_hp,    13.2_hp,  &
                             13.2_hp,     13.2_hp  /)

    Inst%BESOU(:,29) = (/     3.9_hp,      4.2_hp,     4.85_hp,     6.0_hp,  &
                              7.6_hp,      9.0_hp,      9.2_hp,     9.2_hp,  &
                              9.2_hp,      9.2_hp  /)

    Inst%BESOU(:,30) = (/     3.0_hp,     3.05_hp,     3.35_hp,     4.2_hp,  &
                              5.3_hp,      5.9_hp,     6.25_hp,    6.25_hp,  &
                             6.25_hp,     6.25_hp  /)

    Inst%BESOU(:,31) = (/    2.05_hp,      2.1_hp,     2.32_hp,     2.9_hp,  &
                              3.4_hp,      3.9_hp,      4.1_hp,     4.1_hp,  &
                              4.1_hp,      4.1_hp  /)

    Inst%BESOU(:,32) = (/    1.45_hp,     1.43_hp,     1.65_hp,    2.03_hp,  &
                              2.4_hp,     2.75_hp,     2.65_hp,    2.65_hp,  &
                             2.65_hp,     2.65_hp  /)

    Inst%BESOU(:,33) = (/    1.04_hp,     1.08_hp,     1.21_hp,     1.5_hp,  &
                             1.68_hp,      1.8_hp,      1.8_hp,     1.8_hp,  &
                              1.8_hp,      1.8_hp  /)

    ! All the numbers of BESOU need to be multiplied by 1e-5 in order to put
    ! them into the correct data range.  NOTE: This multiplication statement
    ! needs to be preserved here in order to  ensure identical output to the
    ! prior code! (bmy, 7/7/14)
    Inst%BESOU = Inst%BESOU * 1.e-5_hp

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
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstGet
!
! !DESCRIPTION: Subroutine InstGet returns a poiner to the desired instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstGet ( Instance, Inst, RC, PrevInst )
!
! !INPUT PARAMETERS:
!
    INTEGER                             :: Instance
    TYPE(MyInst),     POINTER           :: Inst
    INTEGER                             :: RC
    TYPE(MyInst),     POINTER, OPTIONAL :: PrevInst
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst),     POINTER    :: PrvInst

    !=================================================================
    ! InstGet begins here!
    !=================================================================

    ! Get instance. Also archive previous instance.
    PrvInst => NULL()
    Inst    => AllInst
    DO WHILE ( ASSOCIATED(Inst) )
       IF ( Inst%Instance == Instance ) EXIT
       PrvInst => Inst
       Inst    => Inst%NextInst
    END DO
    IF ( .NOT. ASSOCIATED( Inst ) ) THEN
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Pass output arguments
    IF ( PRESENT(PrevInst) ) PrevInst => PrvInst

    ! Cleanup & Return
    PrvInst => NULL()
    RC = HCO_SUCCESS

  END SUBROUTINE InstGet
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstCreate
!
! !DESCRIPTION: Subroutine InstCreate creates a new instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstCreate ( ExtNr, Instance, Inst, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,       INTENT(IN)       :: ExtNr
!
! !OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(  OUT)    :: Instance
    TYPE(MyInst),  POINTER          :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(INOUT)    :: RC
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst), POINTER          :: TmpInst
    INTEGER                        :: nnInst

    !=================================================================
    ! InstCreate begins here!
    !=================================================================

    ! ----------------------------------------------------------------
    ! Generic instance initialization
    ! ----------------------------------------------------------------

    ! Initialize
    Inst => NULL()

    ! Get number of already existing instances
    TmpInst => AllInst
    nnInst = 0
    DO WHILE ( ASSOCIATED(TmpInst) )
       nnInst  =  nnInst + 1
       TmpInst => TmpInst%NextInst
    END DO

    ! Create new instance
    ALLOCATE(Inst)
    Inst%Instance = nnInst + 1
    Inst%ExtNr    = ExtNr

    ! Attach to instance list
    Inst%NextInst => AllInst
    AllInst       => Inst

    ! Update output instance
    Instance = Inst%Instance

    ! ----------------------------------------------------------------
    ! Type specific initialization statements follow below
    ! ----------------------------------------------------------------

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE InstCreate
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!BOP
!
! !IROUTINE: InstRemove
!
! !DESCRIPTION: Subroutine InstRemove creates a new instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstRemove ( Instance )
!
! !INPUT PARAMETERS:
!
    INTEGER                         :: Instance
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER                     :: RC
    TYPE(MyInst), POINTER       :: PrevInst
    TYPE(MyInst), POINTER       :: Inst

    !=================================================================
    ! InstRemove begins here!
    !=================================================================

    ! Init
    PrevInst => NULL()
    Inst     => NULL()

    ! Get instance. Also archive previous instance.
    CALL InstGet ( Instance, Inst, RC, PrevInst=PrevInst )

    ! Instance-specific deallocation
    IF ( ASSOCIATED(Inst) ) THEN

       ! Pop off instance from list
       IF ( ASSOCIATED(PrevInst) ) THEN

          ! Cleanup module arrays
          IF ( ASSOCIATED(Inst%EmissRn222    ) ) DEALLOCATE(Inst%EmissRn222    )
          IF ( ASSOCIATED(Inst%EmissBe7      ) ) DEALLOCATE(Inst%EmissBe7      )
          IF ( ASSOCIATED(Inst%EmissBe7Strat ) ) DEALLOCATE(Inst%EmissBe7Strat )
          IF ( ASSOCIATED(Inst%EmissBe10     ) ) DEALLOCATE(Inst%EmissBe10     )
          IF ( ASSOCIATED(Inst%EmissBe10Strat) ) DEALLOCATE(Inst%EmissBe10Strat)
          IF ( ASSOCIATED(Inst%LATSOU        ) ) DEALLOCATE(Inst%LATSOU        )
          IF ( ASSOCIATED(Inst%PRESOU        ) ) DEALLOCATE(Inst%PRESOU        )
          IF ( ASSOCIATED(Inst%BESOU         ) ) DEALLOCATE(Inst%BESOU         )

          PrevInst%NextInst => Inst%NextInst
       ELSE
          AllInst => Inst%NextInst
       ENDIF
       DEALLOCATE(Inst)
       Inst => NULL()
    ENDIF

   END SUBROUTINE InstRemove
!EOC
END MODULE HCOX_GC_RnPbBe_Mod
