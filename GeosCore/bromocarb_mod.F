! 8/13/07, jpp, bromocarb_mod.f
      MODULE BROMOCARB_MOD

!
!******************************************************************************
!  Module BROMOCARB_MOD contains variables and routines for the GEOS-CHEM
!  bromocarbon simulation
!
!  Module Variables:
!  ============================================================================
!  (1 ) ECHBr3_oc  (REAL*8 ) : CHBr3 emissions from oceans [kg/s]
!  (2 ) ECH2Br2_oc (REAL*8 ) : CH2Br2 emissions from oceans [kg/s]


!******************************************************************************
!
      IMPLICIT NONE
#     include "define.h"

      ! Make everything PRIVATE
      PRIVATE

      ! except for these routines
      PUBLIC :: EMIS_CHBr3
      PUBLIC :: EMIS_CH2Br2
      PUBLIC :: SET_CH3Br
      PUBLIC :: SET_BrO
!jp      PUBLIC :: SRC_VSLB
      PUBLIC :: SEA_SURFACE
      PUBLIC :: INIT_BROMOCARB
      PUBLIC :: CLEANUP_BROMOCARB

      PUBLIC :: Br_SCALING

      !===============================================
      ! Declare Module Variables
      !===============================================

      ! For scaling bromine emissions (mpayer, 5/15/12)
      REAL*8              :: Br_SCALING

      ! For TOT_OC_AREA:
      ! first element represents region a = between 20S and 20
      ! 2nd element represents region b = (20 to 50, NH + SH)
      REAL*8,  ALLOCATABLE :: TOT_OC_AREA(:)
      REAL*8,  ALLOCATABLE :: FRAC_IN_ab(:,:,:) ! (region, I, J)
                                                ! if region = 1, band btwn 20S and 20N
                                                ! region = 2 is band (20 to 50, NH + SH)
      REAL*8,  ALLOCATABLE :: A_M2(:) ! surface area of the box;
                                      ! only depends on latitude.
      REAL*8,  ALLOCATABLE :: Kg_CHBr3_sec(:)
      REAL*8,  ALLOCATABLE :: Kg_CH2Br2_sec(:)
      REAL*4,  ALLOCATABLE :: ARRAY(:,:,:)
      REAL*8,  ALLOCATABLE :: NPP(:,:)
      real*8 :: total_ocean_sarea
      real*8 :: molec_chbr3_sec, tot_npp
      ! Qing Liang's emissions variables
      real*8, allocatable :: chbr3_emiss(:,:), 
     &     ch2br2_emiss(:,:)


      !===============================================
      ! Useful Constants
      !===============================================
!      REAL*8,  PARAMETER   :: Gg_CHBr3_yr  = 400d0
      REAL*8,  PARAMETER   :: Gg_CHBr3_yr  = 500d0
      ! jpp, 7/15/09: FLAG, trying the quack parameterization and scaling
      !               by NPP to get spatial distribution.
!      REAL*8,  PARAMETER   :: Gg_CHBr3_yr  = 2530d0 ! trying Quack 2003 scaled by NPP
      REAL*8,  PARAMETER   :: Gg_CH2Br2_yr = 113d0
      REAL*8,  PARAMETER   :: Gg2Kg = 1.d6      ![kg/Gg]
      REAL*8,  PARAMETER   :: yr2sec = 3.1536d7 !number of seconds per year
      REAL*8,  PARAMETER   :: EMISS_F_a = 0.75d0, EMISS_F_b = 0.25d0 
            ! fraction of total VSLB emissions in region a (20S to 20N) and 
            ! region b (20 to 50, NH + SH)
      REAL*8,  PARAMETER   :: MWT_CHBr3  = 2.53d-1 ! Molecular weight of CHBr3  [kg/mol]
      REAL*8,  PARAMETER   :: MWT_CH2Br2 = 1.74d-1 ! Molecular weight of CH2Br2 [kg/mol]
      INTEGER, PARAMETER   :: NUM_REGIONS = 2    ! number of elements for TOT_OC_AREA

      !jpp: debugging
      ! needs more attention, better solution
      LOGICAL, SAVE   :: FIRST_INIT = .TRUE.

      ! for new month... can use this for calling npp
      INTEGER, SAVE             :: MONTHSAVE = 0 


      !=================================================================     
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS


!===============================================================================

      FUNCTION EMIS_CHBr3(I,J) RESULT( E_R_CHBr3 )

!===============================================================================
! jpp 8/23/07: this is the emissions driver for very short lived bromocarbon
!              (VSLB) species and the one long-lived compound, methyl bromide
!              (CH3Br)
!
! NOTES:
! 
! only ocean emissions for all bromocarbons
! plan: (1) Bromoform: 400 Gg CHBr3/yr emitted from ocean, broken into
!              latitudinal bands: 75% between 20deg south and 20deg north
!                 25% between 20deg and 50deg north and south
!              - This emission scheme follows the work of
!                Warwick et al. (2006) Global Modeling of Bromocarbons
!                   --> scheme A (eventually, should try B as well,
!                       with coastal and shelf emissions...)
!                & Yang et al. (2005) Tropospheric Bromine Chemistry
!       (2) Dibromomethane:
!              - Warwick et al. use same spatial scheme as used for
!                bromoform in scenario 3 (lat bands above...)
!                   --> schemes A & B: 113 Gg CH2Br2/yr global flux
!              - they say they require higher emissions than some previous
!                lit probably because they center emissions in the tropics,
!                yielding shorter lifetimes for bromocarbons...
!====================================================================
!
      ! References to F90 modules
      USE ERROR_MOD,         ONLY : DEBUG_MSG, is_safe_div
      USE ERROR_MOD,         ONLY : geos_chem_stop, it_is_finite
      USE LOGICAL_MOD,       ONLY : LPRT, LWARWICK_VSLS  !, LDYNOCEAN
      USE DIAG_MOD,          ONLY : AD46
      USE BPCH2_MOD,         ONLY : GET_TAU0,  READ_BPCH2
      USE TIME_MOD,          ONLY : GET_MONTH
      USE GRID_MOD,          ONLY : GET_YMID,  GET_YEDGE
      USE TRANSFER_MOD,      ONLY : transfer_2d


!jp      USE TIME_MOD,          ONLY : GET_MONTH, ITS_A_NEW_MONTH !use for CH3Br
!jp      USE TRACER_MOD,        ONLY : STT

      USE CMN_SIZE_MOD            ! Size parameters
      USE CMN_DIAG_MOD            ! ND57 -- diagnostics
      USE COMODE_LOOP_MOD         ! Avogadro's #, called 'AVG' = 6.02252d+23
      
!jp--for altix
!jp      IMPLICIT NONE

      !=================================================================
      ! Local variables
      !=================================================================

      LOGICAL, SAVE :: FIRST = .TRUE. 
!jp      LOGICAL, INTENT(IN) :: test_in
!jp      INTEGER       :: THISMONTH  !will be used for CH3Br
      INTEGER, INTENT(IN) :: I,   J
      REAL*8              :: E_R_CHBr3
      ! Northern lats (>30 degrees) seasonal cycle scaling
      real*8,  dimension(12) :: nh_scale
      integer                :: this_month

      
      ! variables for dealing with NPP
      REAL*8               :: XTAU
      integer              ::  ix, iy
      CHARACTER(LEN=255)   :: FILENAME
      logical              :: err_check

!      REAL*8, SAVE        :: test
!jp      REAL*8, PARAMETER   :: AVG = 6.02252d+23
      !=================================================================
      ! EMISSBROMOCARB begins here!
      !=================================================================

      ! set the scaling factors: jan to december
      ! see Parrella et al. [2011]
      nh_scale(1)=0.51d0; nh_scale(2) = 0.56d0; nh_scale(3) = 0.61d0
      nh_scale(4)=0.80d0; nh_scale(5) = 0.96d0; nh_scale(6) = 1.12d0
      nh_scale(7) = 1.11d0; nh_scale(8) = 1.07d0
      nh_scale(9) = 1.027d0
      nh_scale(10) = 0.94d0; nh_scale(11) = 0.78d0
      nh_scale(12) = 0.671d0

      !----------------------------------------------------------------------
      ! Initialize arrays: SEA_SURFACE is also called within INIT_BROMOCARB 
      ! to retrieve surface areas for chosen lat-band regions
      !----------------------------------------------------------------------


      !jpp: debugging
      ! this is a quick fix... need a better solution to
      ! pull out CHBr3 and CH2Br2 emissions
      IF (FIRST_INIT) THEN
         CALL INIT_BROMOCARB
         FIRST_INIT = .FALSE.
      ENDIF

      !----------------------
      ! CHBr3 Emissions
      !----------------------
      this_month = get_month()

      ! --------------------------------------------
      ! Return with 0'd emission if the VSL source
      ! has been turned off inside input.geos
      ! --------------------------------------------
      IF ( .not. LWARWICK_VSLS ) then
         E_R_CHBr3 = 0.d0
         RETURN
      ENDIF

      ! --------------------------------------------------------
      ! Calculating the CHBr3 emission rate [molecules/box/s]
      ! from Qing Liang's emissions distribution.
      ! --------------------------------------------------------
      if ( GET_YMID(I,J,1) > 30.d0 ) then
         ! use seasonal scaling
         E_R_CHBr3 = chbr3_emiss(I,J) * 
     &        nh_scale(this_month)
      else
         E_R_CHBr3 = chbr3_emiss(I,J)
      endif


      return

      !----------------------------
      ! Return to calling program
      !----------------------------

      END FUNCTION EMIS_CHBr3

!===============================================================================


!===============================================================================

      FUNCTION EMIS_CH2Br2(I,J) RESULT( E_R_CH2Br2 )

!===============================================================================
! jpp 8/23/07: this is the emissions driver for very short lived bromocarbon
!              (VSLB) species and the one long-lived compound, methyl bromide
!              (CH3Br)
!
! NOTES:
! 
! only ocean emissions for all bromocarbons
! plan: (1) Bromoform: 400 Gg CHBr3/yr emitted from ocean, broken into
!              latitudinal bands: 75% between 20deg south and 20deg north
!                 25% between 20deg and 50deg north and south
!              - This emission scheme follows the work of
!                Warwick et al. (2006) Global Modeling of Bromocarbons
!                   --> scheme A (eventually, should try B as well,
!                       with coastal and shelf emissions...)
!                & Yang et al. (2005) Tropospheric Bromine Chemistry
!       (2) Dibromomethane:
!              - Warwick et al. use same spatial scheme as used for
!                bromoform in scenario 3 (lat bands above...)
!                   --> schemes A & B: 113 Gg CH2Br2/yr global flux
!              - they say they require higher emissions than some previous
!                lit probably because they center emissions in the tropics,
!                yielding shorter lifetimes for bromocarbons...
!====================================================================
!
      ! References to F90 modules
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE LOGICAL_MOD,       ONLY : LPRT, LWARWICK_VSLS  !, LDYNOCEAN
      USE DIAG_MOD,          ONLY : AD46
!jp      USE TIME_MOD,          ONLY : GET_MONTH, ITS_A_NEW_MONTH !use for CH3Br
!jp      USE TRACER_MOD,        ONLY : STT

      USE CMN_SIZE_MOD            ! Size parameters
      USE CMN_DIAG_MOD            ! ND57 -- diagnostics
      USE COMODE_LOOP_MOD         ! Avogadro's #, called 'AVG' = 6.02252d+23
      
!jp--for altix
!jp      IMPLICIT NONE

      !=================================================================
      ! Local variables
      !=================================================================

      LOGICAL, SAVE :: FIRST = .TRUE. 
!jp      LOGICAL, INTENT(IN) :: test_in
!jp      INTEGER       :: THISMONTH  !will be used for CH3Br
      INTEGER, INTENT(IN) :: I,   J
      REAL*8              :: E_R_CH2Br2
!      REAL*8, SAVE        :: test
!jp      REAL*8, PARAMETER   :: AVG = 6.02252d+23
      !=================================================================
      ! EMISSBROMOCARB begins here!
      !=================================================================

      !----------------------------------------------------------------------
      ! Initialize arrays: SEA_SURFACE is also called within INIT_BROMOCARB 
      ! to retrieve surface areas for chosen lat-band regions
      !----------------------------------------------------------------------

      !jpp: debugging
      ! this is a quick fix... need a better solution to
      ! pull out CHBr3 and CH2Br2 emissions
      IF (FIRST_INIT) THEN
         CALL INIT_BROMOCARB
         FIRST_INIT = .FALSE.
      ENDIF


      !----------------------
      ! CHBr3 Emissions
      !----------------------

      ! --------------------------------------------
      ! Return with 0'd emission if the VSL source
      ! has been turned off inside input.geos
      ! --------------------------------------------
      IF ( .not. LWARWICK_VSLS ) then
         E_R_CH2Br2 = 0.d0
         RETURN
      ENDIF

      E_R_CH2Br2 = ch2br2_emiss(I,J)

      !----------------------------
      ! Return to calling program
      !----------------------------

      END FUNCTION EMIS_CH2Br2

!===============================================================================

!------------------------------------------------------------------------------

!jp      FUNCTION SRC_VSLB(I,J) RESULT( ECHBr3 ) ! in [molecules/box/s]
!jp!
!jp!******************************************************************************
!jp!  Subroutine SRC_VSLB is the subroutine for bromoform emissions.  
!jp!  Emissions of CHBr3 will be distributed throughout the boundary layer. 
!jp!
!jp! 
!jp!  Arguments as Input/Output:
!jp!  ============================================================================
!jp!  (1 ) (REAL*8) : Tracer concentration of CHBr3 [kg]
!jp!
!jp!  NOTES:
!jp!  (1 ) SMVGear will mix emissions evenly throughout the PBL for the fullchem()
!jp!       runs... so we don't need a mixing scheme for fullchem, (jpp, 8/24/07)
!jp!  (2 ) We do need to implement our own PBL mixing for an offline run
!jp!
!jp!  (3 ) Can't use STT either for a fullchem run, though that will be the method
!jp!       for an off-line simulation... for fullchem, we must write the emissions
!jp!       rates into the array EMISRR(I,J,N) inside the emissdr.f emission driver
!jp!       - necessary units for STT: [kg]
!jp!       - necessary units for EMISRR: [molecules/box/s]
!jp!
!jp!******************************************************************************
!jp!
!jp      ! Reference to diagnostic arrays
!jp!jp      USE DIAG03_MOD,   ONLY : AD03, ND03
!jp!jp      USE ERROR_MOD,    ONLY : ERROR_STOP
!jp!jp      USE LOGICAL_MOD,  ONLY : LSPLIT
!jp      USE PBL_MIX_MOD,  ONLY : GET_FRAC_OF_PBL, GET_PBL_MAX_L
!jp      USE TIME_MOD,     ONLY : GET_TS_EMIS !returns in minutes
!jp!jp      USE TRACER_MOD,   ONLY : STT
!jp      USE TRACERID_MOD, ONLY : IDTCHBr3
!jp      USE DIAG_MOD,     ONLY : AD57
!jp 
!jp 
!jp#     include "CMN_SIZE"     ! Size parameters
!jp#     include "CMN_DEP"      ! FRCLND == fraction of box that's land
!jp#     include "CMN_DIAG"     ! ND57 -- diagnostics
!jp#     include "comode.h"     ! Avogadro's #, called 'AVG' = 6.02252d+23
!jp 
!jp      ! Local variables
!jp      INTEGER, INTENT(IN)  :: I,     J   !,     N ,    PBL_MAX
!jp 
!jp      ! For a given box (I, J):
!jp      ! Calculating the CHBr3 emission rate [molecules/box/s]
!jp      ECHBr3(I,J) = ( EMISS_F_a*FRAC_IN_ab(1,I,J)/TOT_OC_AREA(1) + 
!jp     &     EMISS_F_b*FRAC_IN_ab(2,I,J)/TOT_OC_AREA(2) ) *
!jp     &     Kg_CHBr3_sec(1) * AVG / MWT_CHBr3
!jp 
!jp      ! For ND575 diagnostic: store emission rate in [kg/m2/s]
!jp      IF ( ND57 > 0 ) THEN
!jp         AD57(I,J) = ( ECHBr3(I,J) / A_M2 ) * ( MWT_CHBr3 / AVG )
!jp      END IF
!jp 
!jp      ! Return to calling program
!jp      END FUNCTION SRC_VSLB

!=============================================================================


!=============================================================================
      SUBROUTINE SEA_SURFACE

!-----------------------------------------------------------------------------
!   jpp, 8/16/07
!
!   The purpose of this subroutine is to calculate the total sea
!   surface area within two specified regions:
!   (1) total area between 20S and 20N
!   (2) total area between 20 and 50 degrees, North + South
!
!   These surface area values are used to set the emission fluxes
!   for each of the aforementioned regions
!
!   NOTES:                               
!-----------------------------------------------------------------------------

      !===================================
      ! Modules Called within this routine
      !===================================
      USE GRID_MOD,     ONLY : GET_AREA_M2, GET_YEDGE


      USE CMN_SIZE_MOD      !grid box loops and other size params
                            !includes jjpar, iipar...
      USE CMN_DEP_MOD       !FRCLND = returns land-fraction of given box

!jp-- for altix
!jp      IMPLICIT NONE

      !===================================
      ! Naming variables
      !===================================

      ! Local variables, adapted from RnPbBe_mod.f
      INTEGER     :: I,          J,          L,         N
      REAL*8      :: LAT_H,      LAT_L,      F_LAND,    F_WATER
      REAL*8      :: F_ABOVE_50, F_BELOW_20
      REAL*8      :: LAT_S,      LAT_N

      !===================================
      ! Initializing Variables
      !===================================

      F_ABOVE_50    = 0d0
      F_BELOW_20    = 0d0
      LAT_N         = 0d0
      LAT_S         = 0d0
      F_LAND        = 0d0
      F_WATER       = 0d0


      !===================================
      ! SUBROUTINE SEA_SURFACE BEGINS HERE
      !===================================

      ! Loop over grid boxes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Get ABS( latitude ) at S and N edges of grid box
         LAT_S      = ABS( GET_YEDGE(I,J,1)   ) 
         LAT_N      = ABS( GET_YEDGE(I,J+1,1) )
         LAT_H      = MAX( LAT_S, LAT_N )
         LAT_L      = MIN( LAT_S, LAT_N ) 

         IF ( LAT_L >= 50d0 ) then
            F_ABOVE_50 = 1d0
            F_BELOW_20 = 0d0
         ELSE IF ( (LAT_H > 50d0) .and. (LAT_L < 50d0) ) then
            F_ABOVE_50 = ( LAT_H - 50d0 ) / ( LAT_H - LAT_L )
            F_BELOW_20 = 0d0
         ELSE IF ( (LAT_H <= 50d0) .and. (LAT_L >= 20d0) ) then
            F_ABOVE_50 = 0d0
            F_BELOW_20 = 0d0
         ELSE IF ( (LAT_H > 20d0) .and. (LAT_L < 20d0) ) then
            F_ABOVE_50 = 0d0
            F_BELOW_20 = ( 20d0 - LAT_L )/ ( LAT_H - LAT_L )
         ELSE IF ( LAT_H <= 20d0 ) THEN
            F_ABOVE_50 = 0d0
            F_BELOW_20 = 1d0
         END IF

         ! Grid box surface area [m2]-- it's only a fn' of latitude
         A_M2(J)    = GET_AREA_M2( I, J, 1 )


! use this if you have problems selecting only ocean.
! they used something like this in ocean_mercury_mod.
!jp         !======================================
!jp         ! Make sure we are in an ocean box     
!jp         !======================================
!jp         IF ( ( ALBD(I,J) <= 0.4d0 ) .and.      
!jp     &        ( FRAC_L    <  0.8d0 )  THEN


            ! Fraction of grid box that is land
            F_LAND  = FRCLND(I,J)

            ! Fraction of grid box that is water
            F_WATER = 1d0 - F_LAND

            ! Find the ocean fraction for box in region a:
            ! Between 20S and 20N
            FRAC_IN_ab(1,I,J) = F_BELOW_20 * F_WATER

            ! Find the ocean fraction for box in region b:
            ! Between 20 and 50, N + S
            FRAC_IN_ab(2,I,J) = (1d0 - F_BELOW_20 - F_ABOVE_50) *
     &           F_WATER
               ! note: if we just add the fractions of box below 50
               ! and above 20, then we would need to subtract the
               ! union to avoid double counting. Trick to get around
               ! this is subtracting the total area outside of b from
               ! 1d0
            ! sum up the total areas in both reagions
            TOT_OC_AREA(1) = TOT_OC_AREA(1) + FRAC_IN_ab(1,I,J) *
     &           A_M2(J)
            TOT_OC_AREA(2) = TOT_OC_AREA(2) + FRAC_IN_ab(2,I,J) *
     &           A_M2(J)

      END DO                 !END i-loop over longitudes
      END DO                 !END j-loop over latitudes

      END SUBROUTINE SEA_SURFACE
!=============================================================================

!=============================================================================
      SUBROUTINE SET_CH3Br( N_TRACERS, TCVV, AD, STT, unit_flag )

!-----------------------------------------------------------------------------
!   jpp, 2/12/08
!
!   The purpose of this subroutine is to set CH3Br Concentrations
!   in the planetary boundary layer. Based on latitude bands
!   (1) 90-55N, (2) 55N-0, (3) 0-55S, (4) 55-90S
!
!   Values for setting pbl flux were determined by surface
!   measurements from NOAA 2006 data.
!
!  Arguments as Input:
!  ======================================================================
!  (1 ) NTRACE (INTEGER) : 
!  (2 ) TCVV   (REAL*8 ) : Array containing [Air MW / Tracer MW] for tracers
!  (3 ) AD     (REAL*8 ) : Array containing grid box air masses
!
!  Arguments as Input/Output:
!  ======================================================================
!  (4 ) STT    (REAL*8 ) : Array containing tracer conc. [kg] in this case
!
!   NOTES:
!   1) STT is converted back and forth between units of [kg] and
!      [v/v]. Placement of the call to SET_CH3Br in main.f (it's
!      with the emissions) means that it should be in [kg].
!-----------------------------------------------------------------------------

      !===================================
      ! Modules Called within this routine
      !===================================
      USE GRID_MOD,     ONLY : GET_YMID, GET_YEDGE
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_MAX_L
      USE TRACER_MOD,   ONLY : TRACER_NAME
      USE LOGICAL_MOD,  ONLY : LWARWICK_VSLS

      USE CMN_SIZE_MOD       ! grid box loops and other size params
                             ! includes jjpar, iipar...
      USE CMN_DEP_MOD        ! FRCLND = returns land-fraction of given box

!jp-- for altix
!jp      IMPLICIT NONE

      !===================================
      ! Naming variables
      !===================================

      ! Arguments
      LOGICAL, INTENT(IN)    :: unit_flag
      INTEGER, INTENT(IN)    :: N_TRACERS 
      INTEGER                :: I,  J,  L,  N
      REAL*8,  INTENT(IN)    :: TCVV(N_TRACERS)
      REAL*8,  INTENT(IN)    :: AD(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,N_TRACERS)


      REAL*8                 :: LAT_MID
      ! CH3Br values ( from pptv )
      REAL*8, PARAMETER      :: gt55N = 8.35d-12 
      REAL*8, PARAMETER      :: gt0_lt55N = 8.27d-12
      REAL*8, PARAMETER      :: lt0_gt55S = 6.94d-12
      REAL*8, PARAMETER      :: lt55S = 6.522d-12

      ! CH3Br variable
      REAL*8                 :: CH3Br_conc
      INTEGER                :: CH3Br_sel

      ! for testing
      INTEGER :: count
      LOGICAL :: FIRST_COUNT = .TRUE.

      !===================================
      ! Initializing Variables
      !===================================

      LAT_MID    = 0d0
      CH3Br_conc = 0d0

      count = 0
      ! get ID # for CH3Br
      DO N = 1, N_TRACERS
         IF(TRACER_NAME(N) == 'CH3Br') THEN
            CH3Br_sel = N
            print *, 'jpp: test worked in bromocarb_mod.f'
            print *, N
         ENDIF
      ENDDO

      !===================================
      ! SUBROUTINE SET_CH3Br BEGINS HERE
      !===================================

      ! -----------------------------------------
      ! If we aren't using bromocarbons, then
      ! set the CH3Br equal to zero.
      ! -----------------------------------------
      IF ( .not. LWARWICK_VSLS ) THEN
         STT(:,:,:,CH3Br_sel) = 0.d0
         RETURN
      ENDIF

      ! jpp: 2/12/08, continue to edit!!!
      ! jpp: 2/13/08, trying based on box midpoints

      ! Loop over grid boxes
      DO J = 1, JJPAR
         DO I = 1, IIPAR
            DO L = 1, LLPAR
               IF ( GET_FRAC_UNDER_PBLTOP( I, J, L ) > 0d0 ) THEN
                  
                  ! testing
                  count = count + 1

                  ! base lat band selection on midpoint
                  ! latitude of the box
                  LAT_MID = GET_YMID( I, J, L )

                  ! Selecting the latitude bands:
                  IF ( LAT_MID > 55d0 ) THEN
                     CH3Br_conc = gt55N
                  ELSEIF ( (LAT_MID >= 0d0) .and. (LAT_MID <= 55d0) )
     &                    THEN
                     CH3Br_conc = gt0_lt55N
                  ELSEIF ( (LAT_MID < 0d0) .and. (LAT_MID >= -55d0) )
     &                    THEN
                     CH3Br_conc = lt0_gt55S
                  ELSEIF ( LAT_MID < -55d0 ) THEN
                     CH3Br_conc = lt55S
                  ENDIF

                  !test STT for CH3Br
!                  IF ( STT(I,J,L,CH3Br_sel) == 0d0) THEN
!                     print *, 'count is'
!                     print *, count
!                  ENDIF


                  ! Make sure we're using the correct units
                  if ( unit_flag ) then
                     ! if the flag is true, then STT has been
                     ! converted from kg/box to v/v mixing ratio.
                     ! so we must supply v/v
                     STT(I,J,L,CH3Br_sel) = CH3Br_conc
                  else
                     ! Now convert the [v/v] units to [kg]
                     ! as in convert_units subroutine in dao_mod.f
                     STT(I,J,L,CH3Br_sel) = CH3Br_conc * 
     &                    AD(I,J,L) / TCVV(CH3Br_sel)
                  endif
!                 IF ( FIRST_COUNT ) THEN
!                     print *, 'New STT for CH3Br'
!                     print *, STT(I,J,L,CH3Br_sel), ' [kg] '
!                     print *, CH3Br_conc, ' [v/v] '
!                     print *, 'LAT_MID = ', LAT_MID
!                     FIRST_COUNT = .FALSE.
!                  ENDIF

               ENDIF            ! end selection of PBL boxes

               END DO           !END l-loop over altitudes

         END DO                 !END i-loop over longitudes

      END DO                    !END j-loop over latitudes
      

      END SUBROUTINE SET_CH3Br
!=============================================================================

!=============================================================================
      SUBROUTINE SET_BRO( N_TRACERS, TCVV, AD, SUNCOS, 
     &     STT, unit_flag )

!-----------------------------------------------------------------------------
!   jpp, 2/12/08
!
!   The purpose of this subroutine is to set Bro Concentrations
!   in the planetary boundary layer. Based on latitude bands
!   (1) 90-55N, (2) 55N-0, (3) 0-55S, (4) 55-90S
!
!   Values for setting pbl flux were determined by surface
!   measurements from NOAA 2006 data.
!
!  Arguments as Input:
!  ======================================================================
!  (1 ) NTRACE (INTEGER) : 
!  (2 ) TCVV   (REAL*8 ) : Array containing [Air MW / Tracer MW] for tracers
!  (3 ) AD     (REAL*8 ) : Array containing grid box air masses
!
!  Arguments as Input/Output:
!  ======================================================================
!  (4 ) STT    (REAL*8 ) : Array containing tracer conc. [kg] in this case
!
!   NOTES:
!   1) STT is converted back and forth between units of [kg] and
!      [v/v]. Placement of the call to SET_Bro in main.f (it's
!      with the emissions) means that it should be in [kg].
!-----------------------------------------------------------------------------

      !===================================
      ! Modules Called within this routine
      !===================================
      USE GRID_MOD,     ONLY : GET_YMID, GET_YEDGE
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_MAX_L
      USE TRACER_MOD,   ONLY : TRACER_NAME
      USE LOGICAL_MOD,  ONLY : LFIX_PBL_BRO
      USE DAO_MOD,      ONLY : IS_WATER

      USE CMN_SIZE_MOD       ! grid box loops and other size params
                             ! includes jjpar, iipar...
      USE CMN_DEP_MOD        ! FRCLND = returns land-fraction of given box

!jp-- for altix
!jp      IMPLICIT NONE

      !===================================
      ! I/O variables
      !===================================

      ! Arguments
      LOGICAL, INTENT(IN)    :: unit_flag
      INTEGER, INTENT(IN)    :: N_TRACERS 
      INTEGER                :: I,  J,  L,  N
      REAL*8,  INTENT(IN)    :: TCVV(N_TRACERS)
      REAL*8,  INTENT(IN)    :: SUNCOS(MAXIJ) 
      REAL*8,  INTENT(IN)    :: AD(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,N_TRACERS)

      ! --------------------------
      ! Local Variables
      ! --------------------------
      INTEGER :: IJLOOP

      ! Bro variable
      REAL*8                 :: Bro_conc
      INTEGER                :: Bro_sel

      ! for testing
      INTEGER :: count
      LOGICAL :: FIRST_COUNT = .TRUE.



      ! -----------------------------------------
      ! If we aren't using this 1pptv experiment
      ! then return without updating STT array.
      ! -----------------------------------------
      IF ( .not. LFIX_PBL_BRO ) THEN
         RETURN
      ENDIF

      !===================================
      ! Initializing Variables
      !===================================

      ! ------------------------------------------------------
      ! Set the BrO concentration to 1 pptv inside the PBL.
      ! ------------------------------------------------------
      BrO_conc = 1.0d-12

      count = 0
      ! get ID # for Bro
      DO N = 1, N_TRACERS
         IF(TRACER_NAME(N) == 'BrO') THEN
            BrO_sel = N
            print *, 'jpp: test worked in bromocarb_mod.f'
            print *, N
         ENDIF
      ENDDO

      !===================================
      ! SUBROUTINE SET_CH3Br BEGINS HERE
      !===================================

      ! do loop over latitudes
      DO J = 1, JJPAR
         lon_loop: DO I = 1, IIPAR

            ! -----------------------------------------
            ! !. Determine if we're in the marine
            !    boundary layer. If so, procede,
            !    otherwise, skip. (note, we should NOT
            !    0 the concentration... it can be
            !    contributed from other timesteps and
            !    sources.
            ! -----------------------------------------
            IF ( .not. IS_WATER(I,J) ) THEN
               CYCLE lon_loop
            ENDIF

            ! -----------------------------------------
            ! 2. Get the cosine of the SZA to determine
            !    if there's available sunlight for
            !    activation of bromine-chemistry.
            !    If so, set the PBL BrO to 1ppt.
            ! -----------------------------------------
            IJLOOP = ( (J-1) * IIPAR ) + I

            IF ( SUNCOS(IJLOOP) > 0.d0 ) THEN
               BrO_conc = 1.0d-12 ! 1pptv if daytime
            ELSE
               BrO_conc = 0.d0    ! 0 otherwise
            ENDIF

            DO L = 1, LLPAR
               IF ( GET_FRAC_UNDER_PBLTOP( I, J, L ) > 0d0 ) THEN

                  ! Make sure we're using the correct units
                  if ( unit_flag ) then
                     ! if the flag is true, then STT has been
                     ! converted from kg/box to v/v mixing ratio.
                     ! so we must supply v/v
                     STT(I,J,L,BrO_sel) = BrO_conc
                  else
                     ! Now convert the [v/v] units to [kg]
                     ! as in convert_units subroutine in dao_mod.f
                     STT(I,J,L,BrO_sel) = BrO_conc * 
     &                    AD(I,J,L) / TCVV(BrO_sel)
                  endif

               ENDIF            ! end selection of PBL boxes

            END DO              !END l-loop over altitudes

         END DO lon_loop        !END i-loop over longitudes
         
      END DO                    !END j-loop over latitudes
      
      RETURN

      END SUBROUTINE SET_BRO
!=============================================================================

!=============================================================================

      SUBROUTINE INIT_BROMOCARB
!-----------------------------------------------------------------------------
!
!  Subroutine INIT_BROMOCARB allocates and zeroes BROMOCARB 
!  module arrays, and 
!-----------------------------------------------------------------------------  
!

      ! References to F90 modules
      USE ERROR_MOD,         ONLY : ALLOC_ERR, DEBUG_MSG
!      USE GRID_MOD,          ONLY : GET_XMID, GET_YMID
      USE LOGICAL_MOD,       ONLY : LPRT           !, LDYNOCEAN
      USE GRID_MOD,          ONLY : GET_AREA_M2
      USE DIRECTORY_MOD,     ONLY : RUN_DIR

!jpp: use this after i make a logical for bromocarb run
!      USE LOGICAL_MOD, ONLY : LBRAVO


      USE CMN_SIZE_MOD            ! Size parameters
      USE COMODE_LOOP_MOD         ! Avogadro's #, called 'AVG' = 6.02252d+23

      ! Local variables
      integer :: i, j
      INTEGER           :: AS
      character(len=60) :: vsl_file, fmt, vsldir
      character(len=10) :: cnlon

      !=============================
      ! INIT_BROMOCARB begins here!
      !=============================

!jpp:
! write a logical for bromocarbon emissions to be read in 
! from the input.geos file. e.g. from bravo_mod:
!      IF ( .not. LBRAVO ) RETURN
      
      !--------------------------
      ! Allocate and zero arrays
      !--------------------------

!jp      ALLOCATE( ECHBr3( IIPAR, JJPAR ), STAT=AS )
!jp      IF ( AS /= 0 ) CALL ALLOC_ERR( 'bromocarb_emiss_CHBr3' )
!jp      ECHBr3 = 0d0
!jp 
!jp      ALLOCATE( ECH2Br2( IIPAR, JJPAR ), STAT=AS )
!jp      IF ( AS /= 0 ) CALL ALLOC_ERR( 'bromocarb_emiss_CH2Br2' )
!jp      ECH2Br2 = 0d0

      ALLOCATE( TOT_OC_AREA( NUM_REGIONS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'total regional areas allocation' )
      TOT_OC_AREA = 0d0

      ALLOCATE( FRAC_IN_ab( NUM_REGIONS, IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 
     &     'fraction of box in a or b regions' )
      FRAC_IN_ab = 0d0

      ALLOCATE( A_M2( JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'bromocarb_mod: A_M2' )
      A_M2 = 0d0

      ALLOCATE( Kg_CHBr3_sec(1), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'bromocarb_mod: CHBr3[kg/s]' )
      Kg_CHBr3_sec = 0d0

      ALLOCATE( Kg_CH2Br2_sec(1), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'bromocarb_mod: CH2Br2[kg/s]' )
      Kg_CH2Br2_sec = 0d0

      ALLOCATE( array(IIPAR, JJPAR, 1), stat=as)
      if ( as /= 0 ) call Alloc_err( 'bromocarb_mod: array' )

      ALLOCATE( npp(IIPAR, JJPAR), stat=as)
      if ( as /= 0 ) call Alloc_err( 'bromocarb_mod: npp' )

      ALLOCATE( chbr3_emiss(IIPAR, JJPAR), stat=as)
      if ( as /= 0 ) call Alloc_err( 'bromocarb_mod: chbr3_emiss' )

      ALLOCATE( ch2br2_emiss(IIPAR, JJPAR), stat=as)
      if ( as /= 0 ) call Alloc_err( 'bromocarb_mod: ch2br2_emiss' )

      !------------------------------------------
      ! Call Sea Surface to start filling arrays
      !------------------------------------------

      CALL SEA_SURFACE
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSBROMOCARB: SEA_SURFACE' )

      ! ----------------------------------------------------
      ! Now read in Qing Liang's CHBr3 and CH2Br2 Emissions
      ! **
      !  These emissions are in 2 x 2.5 (lat x lon),
      !  or a regridded version 4 x 5 (jpp, 7/31/09)
      ! ----------------------------------------------------

      ! jpp, ** VSL DIRECTORY ** hardwired
!      vsldir = '/home/jpp/qing_liang/'
      VSLDIR = TRIM(RUN_DIR)//'VSL_emissions/'

#if   defined( GRID2x25 )
      vsl_file = 'SLBromine_Emission_2x2.5.dat'

      IF ( LPRT ) CALL DEBUG_MSG( ' ### Pre-Reading VSL files' )

      open(unit=234, file = trim(adjustl(vsldir))//
     &     trim(adjustl(vsl_file)), status='OLD')

      ! set the format to read
      write(cnlon,'(i6)') IIPAR
      fmt = '('//trim(adjustl(cnlon))//'es11.4)'

      ! Read the CHBr3 emissions [kg(CHBr3)/m2/s]
      do j = 1, JJPAR
         read(234,fmt) (chbr3_emiss(i,j), i=1,iipar)
      enddo

      ! Read the CH2Br2 emissions [kg(CH2Br2)/m2/s]
      do j = 1, JJPAR
         read(234,fmt) (ch2br2_emiss(i,j), i=1,iipar)
      enddo

      close(234)


#elif defined( GRID4x5 )

      ! first read the CHBr3 emissions file
      vsl_file = 'chbr3_emission_4x5.dat'

      IF ( LPRT ) CALL DEBUG_MSG( ' ### Reading CHBr3 4x5 emissions' )

      open(unit=234, file = trim(adjustl(vsldir))//
     &     trim(adjustl(vsl_file)), status='OLD')

      ! set the format to read
      write(cnlon,'(i6)') IIPAR
      fmt = '('//trim(adjustl(cnlon))//'es11.4)'

      ! Read the CHBr3 emissions [kg(CHBr3)/m2/s]
      do j = 1, JJPAR
         read(234,fmt) (chbr3_emiss(i,j), i=1,iipar)
      enddo

      ! close the unit
      close(234)

      ! now open and read the CH2Br2 emissions file
      vsl_file = 'ch2br2_emission_4x5.dat'

      IF ( LPRT ) CALL DEBUG_MSG( ' ### Reading CH2Br2 4x5 emissions' )

      open(unit=234, file = trim(adjustl(vsldir))//
     &     trim(adjustl(vsl_file)), status='OLD')

      ! Read the CH2Br2 emissions [kg(CH2Br2)/m2/s]
      do j = 1, JJPAR
         read(234,fmt) (ch2br2_emiss(i,j), i=1,iipar)
      enddo

      close(234)


#endif




      ! Calculating the CHBr3 emission rate [molecules/box/s]
      ! from Qing Liang's emissions distribution...
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box surface area [m2]-- it's only a fn' of latitude
         A_M2(J)    = GET_AREA_M2( I, J, 1 )

            ! Conversions:
            ! kg/m2/s ---> molecules/box/second
            chbr3_emiss(i,j) = chbr3_emiss(i,j) / MWT_CHBr3
     &           * AVG * A_M2(J)

            ch2br2_emiss(i,j) = ch2br2_emiss(i,j) / MWT_CH2Br2
     &           * AVG * A_M2(J)

      ENDDO
      ENDDO

      IF ( LPRT ) CALL DEBUG_MSG( ' ### Post-Reading VSL files' )

      !---------------------------
      ! Return to calling program
      !---------------------------

      END SUBROUTINE INIT_BROMOCARB

!=============================================================================



!=============================================================================

      SUBROUTINE CLEANUP_BROMOCARB
!
!******************************************************************************
!jpp:
!  Subroutine CLEANUP_BROMOCARB deallocates all BROMOCARB 
!  module arrays
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_BROMOCARBON begins here!
      !=================================================================
!jp      IF ( ALLOCATED( ECHBr3       ) ) DEALLOCATE( ECHBr3       )
!jp      IF ( ALLOCATED( ECH2Br2      ) ) DEALLOCATE( ECH2Br2      )
      IF ( ALLOCATED( TOT_OC_AREA   ) ) DEALLOCATE( TOT_OC_AREA   )
      IF ( ALLOCATED( FRAC_IN_ab    ) ) DEALLOCATE( FRAC_IN_ab    )
      IF ( ALLOCATED( A_M2          ) ) DEALLOCATE( A_M2          )
      IF ( ALLOCATED( Kg_CHBr3_sec  ) ) DEALLOCATE( Kg_CHBr3_sec  )
      IF ( ALLOCATED( Kg_CH2Br2_sec ) ) DEALLOCATE( Kg_CH2Br2_sec )
      IF ( ALLOCATED( array     ) ) DEALLOCATE( array     )
      IF ( ALLOCATED( npp     ) ) DEALLOCATE( npp     )
      IF ( ALLOCATED( chbr3_emiss   ) ) DEALLOCATE( chbr3_emiss   )
      IF ( ALLOCATED( ch2br2_emiss  ) ) DEALLOCATE( ch2br2_emiss  )

      ! Return to calling program
      END SUBROUTINE CLEANUP_BROMOCARB

!=============================================================================


      ! End of module
      END MODULE BROMOCARB_MOD
