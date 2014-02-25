!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: tracer_mod
!
! !DESCRIPTION: Module TRACER\_MOD contains GEOS-CHEM tracer array STT plus
!  various other related quantities.  TRACER\_MOD also contains inquiry
!  functions that can be used to determine the type of GEOS-CHEM simulation.
!\\
!\\
! !INTERFACE:
!
      MODULE TRACER_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS: 
!
      PUBLIC :: ITS_A_RnPbBe_SIM
      PUBLIC :: ITS_A_CH3I_SIM
      PUBLIC :: ITS_A_FULLCHEM_SIM
      PUBLIC :: ITS_A_HCN_SIM
      PUBLIC :: ITS_A_TAGOX_SIM
      PUBLIC :: ITS_A_TAGCO_SIM
      PUBLIC :: ITS_A_C2H6_SIM
      PUBLIC :: ITS_A_CH4_SIM
      PUBLIC :: ITS_AN_AEROSOL_SIM
      PUBLIC :: ITS_A_MERCURY_SIM
      PUBLIC :: ITS_A_CO2_SIM
      PUBLIC :: ITS_A_H2HD_SIM
      PUBLIC :: ITS_A_POPS_SIM
      PUBLIC :: ITS_NOT_COPARAM_OR_CH4
      PUBLIC :: GET_SIM_NAME
      PUBLIC :: CHECK_STT
      PUBLIC :: CHECK_STT_05x0666
      PUBLIC :: CHECK_STT_025x03125
      PUBLIC :: INIT_TRACER
      PUBLIC :: CLEANUP_TRACER
!
! !PUBLIC DATA MEMBERS:
!
      !========================================================================
      ! Module Variables:
      ! SIM_TYPE       : Number denoting simulation type
      ! N_TRACERS      : Number of GEOS-CHEM tracers
      ! N_MEMBERS      : Max # of constituents a tracer can have
      ! ID_TRACER      : Array of tracer numbers
      ! ID_EMITTED     : Index of which constituent has the emissions
      ! STT            : GEOS-CHEM Tracer array [kg] 
      ! TCVV           : Molecular weight air / molecular weight tracer
      ! TRACER_COEFF   : Coefficient of each tracer constituent
      ! TRACER_MW_G    : Tracer molecular weight [g/mole]
      ! TRACER_MW_KG   : Tracer molecular weight [kg/mole]
      ! TRACER_N_CONST : Array of number of constituents per tracer
      ! TRACER_NAME    : Array of tracer names
      ! TRACER_CONST   : Array of names for tracer constituents
      ! SALA_REDGE_um  : Accum  mode seasalt radii bin edges [um]
      ! SALC_REDGE_um  : Coarse mode seasalt radii bin edges [um]
      ! XNUMOL         : Ratio of (molec/mole) / (kg/mole) = molec/kg
      ! XNUMOLAIR      : XNUMOL ratio for air
      !========================================================================

      ! Scalars
      INTEGER,           PUBLIC              :: SIM_TYPE
      INTEGER,           PUBLIC              :: N_TRACERS
#if   defined( APM )
      INTEGER,           PUBLIC              :: N_APMTRA
#endif
      ! N_MEMBERS increased from 10 to 15 (FP 8/2009)
      INTEGER,           PUBLIC, PARAMETER   :: N_MEMBERS = 15
      REAL*8,            PUBLIC, PARAMETER   :: XNUMOLAIR = 6.022d+23 / 
     &                                                      28.9644d-3

      ! Arrays
      INTEGER,           PUBLIC, ALLOCATABLE :: ID_TRACER(:)
      INTEGER,           PUBLIC, ALLOCATABLE :: ID_EMITTED(:)
      INTEGER,           PUBLIC, ALLOCATABLE :: TRACER_N_CONST(:)
      REAL*8,            PUBLIC, ALLOCATABLE :: TCVV(:)
      REAL*8,            PUBLIC, ALLOCATABLE :: TRACER_COEFF(:,:)
      REAL*8,            PUBLIC, ALLOCATABLE :: TRACER_MW_G(:)
      REAL*8,            PUBLIC, ALLOCATABLE :: TRACER_MW_KG(:)
      REAL*8,            PUBLIC, ALLOCATABLE :: XNUMOL(:)
      CHARACTER(LEN=14), PUBLIC, ALLOCATABLE :: TRACER_NAME(:)
      CHARACTER(LEN=14), PUBLIC, ALLOCATABLE :: TRACER_CONST(:,:)

      ! Define seasalt radii bin edges [um] here since these
      ! need to be used both in "seasalt_mod.f" and "drydep_mod.f"
      REAL*8,            PUBLIC              :: SALA_REDGE_um(2)
      REAL*8,            PUBLIC              :: SALC_REDGE_um(2)
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca -  Initial version
!  (1 ) Added function GET_SIM_NAME (bmy, 5/3/05)
!  (2 ) Removed ITS_A_COPARAM_SIM; the CO-OH param is obsolete (bmy, 6/24/05)
!  (3 ) Added ITS_A_CO2_SIM (pns, bmy, 7/25/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Now added XNUMOL, XNUMOLAIR as module variables (bmy, 10/25/05)
!  (6 ) Added public routine ITS_A_H2HD_SIM (phs, 9/18/07)
!  (7 ) Added public routine ITS_A_POPS_SIM (eck, 9/20/10)
!  16 Feb 2011 - R. Yantosca - Add modifications for APM microphysics (G. Luo)
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!  04 Apr 2013 - R. Yantosca - Removed STT (now in State_Chm)
!EOP
!------------------------------------------------------------------------------
!BOC
!
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_a_rnpbbe_sim
!
! !DESCRIPTION: Function ITS\_A\_RnPbBe\_SIM returns TRUE if we are doing a 
!  GEOS-CHEM Rn-Pb-Be simulation.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_RnPbBe_SIM() RESULT( VALUE )
! 
! !REVISION HISTORY:
!  15 Jul 2004 - R. Yantosca - Initial version
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_RnPbBe_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 1 )

      END FUNCTION ITS_A_RnPbBe_SIM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_a_ch3i_sim
!
! !DESCRIPTION: Function ITS\_A\_CH3I\_SIM returns TRUE if we are doing a
!  GEOS-CHEM CH3I (Methyl Iodide) simulation.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_CH3I_SIM() RESULT( VALUE )
! 
! !REVISION HISTORY: 
!  15 Jul 2004 - R. Yantosca - Initial version
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_CH3I_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 2 )

      END FUNCTION ITS_A_CH3I_SIM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_a_fullchem_sim
!
! !DESCRIPTION: Function ITS\_A\_FULLCHEM\_SIM returns TRUE if we are doing a
!  GEOS-CHEM full chemistry/aerosol simulation (i.e. via SMVGEAR).
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_FULLCHEM_SIM() RESULT( VALUE )
! 
! !REVISION HISTORY: 
!  15 Jul 2004 - R. Yantosca - Initial version
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_FULLCHEM_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 3 )

      END FUNCTION ITS_A_FULLCHEM_SIM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_a_hcn_sim
!
! !DESCRIPTION: Function ITS\_A\_HCN\_SIM returns TRUE if we are doing a
!  GEOS-CHEM HCN (Hydrogen Cyanide) simulation.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_HCN_SIM() RESULT( VALUE )
! 
! !REVISION HISTORY: 
!  15 Jul 2004 - R. Yantosca - Initial version
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_HCN_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 4 )

      END FUNCTION ITS_A_HCN_SIM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_a_tagox_sim
!
! !DESCRIPTION: Function ITS\_A\_TAGOX\_SIM returns TRUE if we are doing a
!  GEOS-CHEM tagged Ox simulation.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_TAGOX_SIM() RESULT( VALUE )
! 
! !REVISION HISTORY: 
!  15 Jul 2004 - R. Yantosca - Initial version
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_TAGOX_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 6 )

      END FUNCTION ITS_A_TAGOX_SIM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_a_tagco_sim
!
! !DESCRIPTION: Function ITS\_A\_TAGCO\_SIM returns TRUE if we are doing a
!  GEOS-CHEM tagged CO simulation.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_TAGCO_SIM() RESULT( VALUE )
! 
! !REVISION HISTORY: 
!  15 Jul 2004 - R. Yantosca - Initial version
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_TAGCO_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 7 )

      END FUNCTION ITS_A_TAGCO_SIM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_a_c2h6_sim
!
! !DESCRIPTION: Function ITS\_A\_C2H6\_SIM returns TRUE if we are doing a
!  GEOS-CHEM C2H6 (Ethane) simulation.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_C2H6_SIM() RESULT( VALUE )
! 
! !REVISION HISTORY: 
!  15 Jul 2004 - R. Yantosca - Initial version
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_RnPbBe_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 8 )

      END FUNCTION ITS_A_C2H6_SIM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_a_ch4_sim
!
! !DESCRIPTION: Function ITS\_A\_CH4\_SIM returns TRUE if we are doing a
!  GEOS-CHEM CH4 (Methane) simulation.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_CH4_SIM() RESULT( VALUE )
! 
! !REVISION HISTORY: 
!  15 Jul 2004 - R. Yantosca - Initial version
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_CH4_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 9 )

      END FUNCTION ITS_A_CH4_SIM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_an_aerosol_sim
!
! !DESCRIPTION: Function ITS\_AN\_AEROSOL\_SIM returns TRUE if we are doing a
!  GEOS-CHEM offline Sulfate/Carbon/dust/seasalt aerosol simulation.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_AN_AEROSOL_SIM() RESULT( VALUE )
! 
! !REVISION HISTORY: 
!  15 Jul 2004 - R. Yantosca - Initial version
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_AN_AEROSOL_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 10 )

      END FUNCTION ITS_AN_AEROSOL_SIM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_a_mercury_sim
!
! !DESCRIPTION: Function ITS\_A\_MERCURY\_SIM returns TRUE if we are doing a
!  GEOS-CHEM Hg0/Hg2/HgP offline mercury simulation.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_MERCURY_SIM() RESULT( VALUE )
! 
! !REVISION HISTORY: 
!  15 Jul 2004 - R. Yantosca - Initial version
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_MERCURY_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 11 )

      END FUNCTION ITS_A_MERCURY_SIM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_a_pops_sim
!
! !DESCRIPTION: Function ITS\_A\_POPS\_SIM returns TRUE if we are doing a
!  GEOS-CHEM offline POPs simulation.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_POPS_SIM() RESULT( VALUE )
! 
! !REVISION HISTORY: 
!  20 Sep 2010 - N.E. Selin  - Initial version
!  26 Nov 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_POPS_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 14 )

      END FUNCTION ITS_A_POPS_SIM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_a_co2_sim
!
! !DESCRIPTION: Function ITS\_A\_CO2\_SIM returns TRUE if we are doing a
!  GEOS-CHEM CO2 offline simulation.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_CO2_SIM() RESULT( VALUE )
! 
! !REVISION HISTORY: 
!  25 Jul 2004 - R. Yantosca - Initial version
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_CO2_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 12 )

      END FUNCTION ITS_A_CO2_SIM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_a_h2hd_sim
!
! !DESCRIPTION: Function ITS\_A\_H2HD\_SIM returns TRUE if we are doing a
!  GEOS-CHEM H2-HD simulation.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_A_H2HD_SIM() RESULT( VALUE )
! 
! !REVISION HISTORY: 
!  18 Sep 2007 - P. Le Sager - Initial version
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_H2HD_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 13 )

      END FUNCTION ITS_A_H2HD_SIM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_not_coparam_or_ch4
!
! !DESCRIPTION: Function ITS\_NOT\_COPARAM\_OR\_CH4 returns TRUE if we are
!  doing a GEOS-CHEM simulation other than CO with parameterized OH or CH4. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_NOT_COPARAM_OR_CH4() RESULT( VALUE )
!
! !REMARKS:
!  The CO-OH param (SIM_TYPE=5) is now obsolete (bmy, 6/24/05)
!
! !REVISION HISTORY: 
!  15 Jul 2004 - R. Yantosca - Initial version
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_NOT_COPARAM_OR_CH4 begins here!
      !=================================================================
      VALUE = ( SIM_TYPE /= 9 )

      END FUNCTION ITS_NOT_COPARAM_OR_CH4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_sim_name
!
! !DESCRIPTION: Function GET\_SIM\_NAME returns the name (e.g. 
!  "NOx-Ox-Hydrocarbon-Aerosol", "Tagged CO", etc.) of the GEOS-CHEM
!  simulation.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_SIM_NAME() RESULT( NAME )
!
! !RETURN VALUE:
!
      CHARACTER(LEN=40) :: NAME
!
! !REVISION HISTORY:
!  03 May 2005 - R. Yantosca - Initial version
!  (1 ) The CO-OH simulation has been removed (bmy, 6/24/05)
!  (2 ) Added CASE blocks for CO2 and H2/HD simulations (bmy, 9/18/07)
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! GET_SIM_NAME begins here!
      !=================================================================

      ! Pick proper name for each simulation type
      SELECT CASE( SIM_TYPE )
         CASE( 1 ) 
            NAME = 'Rn-Pb-Be'
         CASE( 2 ) 
            NAME = 'CH3I'
         CASE( 3 ) 
            NAME = 'NOx-Ox-Hydrocarbon-Aerosol'
         CASE( 4 )
            NAME = 'HCN'
         CASE( 5 )
            NAME = ''
         CASE( 6 )
            NAME = 'Tagged Ox'
         CASE( 7 )
            NAME = 'Tagged CO'
         CASE( 8 ) 
            NAME = 'Tagged C2H6'
         CASE( 9 )
            NAME = 'CH4'
         CASE( 10 ) 
            NAME = 'Offline Aerosol'
         CASE( 11 ) 
            NAME = 'Mercury'
         CASE( 12 )
            NAME = 'CO2'
         CASE( 13 )
            NAME = 'H2 and HD'
         CASE( 14 )
            NAME = 'POPs'
         CASE DEFAULT
            NAME = 'UNKNOWN'
       END SELECT

      END FUNCTION GET_SIM_NAME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_stt
!
! !DESCRIPTION: Subroutine CHECK\_STT checks the STT tracer array for negative
!  values, NaN values, or Infinity values.  If any of these are found, the code
!  will stop with an error message.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_STT( State_Chm, LOCATION )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE ERROR_MOD,          ONLY : IT_IS_NAN
      USE ERROR_MOD,          ONLY : IT_IS_FINITE
      USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN) :: LOCATION
!
! !INPUT PARAMETERS: 
!
      TYPE(ChmState),   INTENT(IN) :: State_Chm   ! Chemistry State object
! 
! !REVISION HISTORY:
!  (1 ) CHECK_STT uses the interfaces defined above -- these will do the
!        proper error checking for either SGI or DEC/Compaq platforms.
!        (bmy, 3/8/01)
!  (2 ) Now call GEOS_CHEM_STOP to shutdown safely.  Now use logicals LNAN,
!        LNEG, LINF to flag if we have error conditions, and then stop the
!        run outside of the parallel DO-loop. (bmy, 11/27/02)
!  (3 ) Bug fix in FORMAT statement: replace missing commas (bmy, 3/23/03)
!  (4 ) Moved from "error_mod.f" to "tracer_mod.f" (bmy, 7/15/04)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL                      :: LNEG, LNAN, LINF
      INTEGER                      :: I,    J,    L,   N, N_TOT_TRC
      
      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHECK_STT begins here!
      !=================================================================

      ! Pick the proper # of tracers
#if   defined( APM )
      N_TOT_TRC = N_TRACERS + N_APMTRA   ! G-C advected tracers + APM tracers
#else
      N_TOT_TRC = N_TRACERS              ! G-C advected tracers only
#endif

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Initialize
      LNEG = .FALSE.
      LNAN = .FALSE.
      LINF = .FALSE.

      ! Loop over grid boxes
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, N_TOT_TRC
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !---------------------------
         ! Check for Negatives
         !---------------------------
!clb fill nega
         STT(I,J,L,N) = ABS(STT(I,J,L,N))
         IF ( STT(I,J,L,N) < 0d0 ) THEN 
!$OMP CRITICAL
            LNEG = .TRUE.
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
            PRINT*, STT(I,J,L,N)
!$OMP END CRITICAL

         !---------------------------
         ! Check for NaN's
         !---------------------------
         ELSE IF ( IT_IS_NAN( STT(I,J,L,N) ) ) THEN
!$OMP CRITICAL
            LNAN = .TRUE.
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
!$OMP END CRITICAL

         !----------------------------
         ! Check STT's for Infinities
         !----------------------------
         ELSE IF ( .not. IT_IS_FINITE( STT(I,J,L,N) ) ) THEN
!$OMP CRITICAL
            LINF = .TRUE.
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
!$OMP END CRITICAL            

         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      !=================================================================
      ! Stop the run if any of LNEG, LNAN, LINF is true
      !=================================================================
      IF ( LNEG .or. LNAN .or. LINF ) THEN
         WRITE( 6, 120 ) TRIM( LOCATION ), LNEG, LNAN, LINF
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! FORMAT statements
      !=================================================================
 100  FORMAT( 'CHECK_STT: STT(',i3,',',i3,',',i3,',',i3,') = ', f13.6 )
 120  FORMAT( 'CHECK_STT: STOP at ', a , 3L2 )

      END SUBROUTINE CHECK_STT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_stt_05x0666
!
! !DESCRIPTION: Subroutine CHECK\_STT\_05x0666 checks the STT tracer array for
!  negative values, NaN values, or Infinity values.  If any of these are found,
!  the STT array will be set to a specified value.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_STT_05x0666( State_Chm, LOCATION )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : IT_IS_NAN
      USE ERROR_MOD,          ONLY : IT_IS_FINITE
      USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN)    :: LOCATION
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState),   INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  05 Mar 2012 - M. Payer    - Initial version based on CHECK_STT and updates
!                              for nested grid by Yuxuan Wang.
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                      :: I,    J,    L,   N
      
      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHECK_STT_05x0666 begins here!
      !=================================================================

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Loop over grid boxes
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !---------------------------
         ! Check for Negatives
         !---------------------------
         IF ( STT(I,J,L,N) < 0d0 ) THEN 
!$OMP CRITICAL
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
            PRINT*, 'Neg STT ' // TRIM( LOCATION ) //
     &              '. SET STT TO BE ZERO.'
            STT(I,J,L,N) = 0d0
!$OMP END CRITICAL

         !---------------------------
         ! Check for NaN's
         !---------------------------
         ELSE IF ( IT_IS_NAN( STT(I,J,L,N) ) ) THEN 
!$OMP CRITICAL
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
            PRINT*, 'NaN STT ' // TRIM( LOCATION ) //
     &              '. SET STT TO BE LOWER LEVEL.'
            STT(I,J,L,N) = STT(I,J,L-1,N)
!$OMP END CRITICAL

         !----------------------------
         ! Check STT's for Infinities
         !----------------------------
         ELSE IF ( .not. IT_IS_FINITE( STT(I,J,L,N) ) ) THEN
!$OMP CRITICAL
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
            PRINT*, 'Inf STT ' // TRIM( LOCATION ) //
     &              '. SET STT TO BE LOWER LEVEL.'
            STT(I,J,L,N) =  STT(I,J,L-1,N)
!$OMP END CRITICAL

         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

 100  FORMAT( ' STT(',i3,',',i3,',',i3,',',i3,') = ', f13.6 )

      END SUBROUTINE CHECK_STT_05x0666
!EOC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_stt_025x03125
!
! !DESCRIPTION: Subroutine CHECK\_STT\_025x03125 checks the STT tracer 
!  array for negative values, NaN values, or Infinity values.  If any of 
!  these are found, the STT array will be set to a specified value.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_STT_025x03125( State_Chm, LOCATION )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE ERROR_MOD,          ONLY : IT_IS_NAN
      USE ERROR_MOD,          ONLY : IT_IS_FINITE
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN)    :: LOCATION
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState),   INTENT(INOUT) :: State_Chm   ! Chemistry State object
! 
! !REVISION HISTORY:
!  05 Mar 2012 - M. Payer    - Initial version based on CHECK_STT and updates
!                              for nested grid by Yuxuan Wang.
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!  07 Jun 2013 - R. Yantosca - Now pass State_Chm object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER         :: I, J, L, N
      
      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHECK_STT_05x0666 begins here!
      !=================================================================

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Loop over grid boxes
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !---------------------------
         ! Check for Negatives
         !---------------------------
!clb         STT (I,J,L,N)=ABS(STT(I,J,L,N))
         IF ( STT(I,J,L,N) < 0d0 ) THEN
!$OMP CRITICAL
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
            PRINT*, 'Neg STT ' // TRIM( LOCATION ) //
     &              '. SET STT TO BE ZERO.'
            STT(I,J,L,N) = 0d0
!$OMP END CRITICAL

         !---------------------------
         ! Check for NaN's
         !---------------------------
         ELSE IF ( IT_IS_NAN( STT(I,J,L,N) ) ) THEN
!$OMP CRITICAL
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
            PRINT*, 'NaN STT ' // TRIM( LOCATION ) //
     &              '. SET STT TO BE LOWER LEVEL.'
            STT(I,J,L,N) = STT(I,J,L-1,N)
!$OMP END CRITICAL

         !----------------------------
         ! Check STT's for Infinities
         !----------------------------
         ELSE IF ( .not. IT_IS_FINITE( STT(I,J,L,N) ) ) THEN
!$OMP CRITICAL
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
            PRINT*, 'Inf STT ' // TRIM( LOCATION ) //
     &              '. SET STT TO BE LOWER LEVEL.'
            STT(I,J,L,N) =  STT(I,J,L-1,N)
!$OMP END CRITICAL

         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

 100  FORMAT( ' STT(',i3,',',i3,',',i3,',',i3,') = ', f13.6 )

      ! Free pointer memory
      NULLIFY( STT )

      END SUBROUTINE CHECK_STT_025x03125
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_tracer
!
! !DESCRIPTION: Subroutine INIT\_TRACER initializes all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_TRACER( am_I_Root, Input_Opt, RC )
!
! !USES:
!
#if   defined( APM ) 
      USE APM_INIT_MOD,       ONLY : APM_NTRACERS
      USE APM_INIT_MOD,       ONLY : LAPM
#endif
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure
! 
! !REVISION HISTORY:
!  15 Jul 2004 - R. Yantosca - Initial version
!  25 Oct 2005 - R. Yantosca - Now allocate XNUMOL
!  16 Feb 2011 - R. Yantosca - Add modifications for APM from G. Luo
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!  26 Feb 2013 - R. Yantosca - Now pass am_I_Root, Input_Opt, RC as arguments
!  07 Mar 2013 - R. Yantosca - Now assign RC=GIGC_SUCCESS outside of APM block
!  04 Apr 2013 - R. Yantosca - Removed STT (now in State_Chm)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS, N_TOT_TRC

#if   defined( APM ) 
      INTEGER :: N_APMTRAC1

      !=================================================================
      ! We are using APM microphysics!
      !=================================================================

      ! Force LAPM=T if we compile with APM=yes makefile option
      ! NOTE: LAPM is somewhat superfluous, since we now can use the
      ! C-preprocessor switch to activate the code for APM simulations
      ! (bmy, 2/17/11)
      LAPM      = .TRUE.

      ! Add additional tracers for APM
      CALL APM_NTRACERS( Input_Opt%N_TRACERS, N_APMTRAC1 )
      N_APMTRA  = N_APMTRAC1  
 
      ! N_TOT_TRC is the total # of G-C advected tracers + APM tracers
      N_TOT_TRC = Input_Opt%N_TRACERS + N_APMTRA 
      
#else

      !=================================================================
      ! We are not using APM microphysics
      !=================================================================

      ! N_TOT_TRC is the total # of G-C advected tracers only
      N_TOT_TRC = Input_Opt%N_TRACERS

#endif

      ! Assume success
      RC        = GIGC_SUCCESS

      !=================================================================
      ! Allocate arrays 
      !=================================================================
      ALLOCATE( ID_TRACER( Input_Opt%N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ID_TRACER' )
      ID_TRACER = 0

      ALLOCATE( ID_EMITTED( Input_Opt%N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ID_EMITTED' )
      ID_EMITTED = 0

      ALLOCATE( TCVV( N_TOT_TRC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCVV' )
      TCVV = 0d0

      ALLOCATE( TRACER_NAME( Input_Opt%N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_NAME' )
      TRACER_NAME = ''

      ALLOCATE( TRACER_MW_G( Input_Opt%N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_MW_G' )
      TRACER_MW_G = 0d0

      ALLOCATE( TRACER_MW_KG( Input_Opt%N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_MW_KG' )
      TRACER_MW_KG = 0d0

      ALLOCATE( TRACER_COEFF( Input_Opt%N_TRACERS, N_MEMBERS ), STAT=AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_COEFF' )
      TRACER_COEFF = 0d0

      ALLOCATE( TRACER_CONST( Input_Opt%N_TRACERS, N_MEMBERS ), STAT=AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_CONST' )
      TRACER_CONST = ''

      ALLOCATE( TRACER_N_CONST( Input_Opt%N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_N_CONST' )
      TRACER_N_CONST = 0

      ALLOCATE( XNUMOL( N_TOT_TRC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'XNUMOL' )
      XNUMOL = 0

      ! Print a debug msg
      IF ( Input_Opt%LPRT .and. am_I_Root ) THEN
         CALL Debug_Msg( '### a INIT_TRACER' )
      ENDIF

      END SUBROUTINE INIT_TRACER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_tracer
!
! !DESCRIPTION: Subroutine CLEANUP\_TRACER deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_TRACER
! 
! !REVISION HISTORY:
!  15 Jul 2004 - R. Yantosca - Initial version
!  25 Oct 2005 - R. Yantosca - Now deallocates XNUMOL
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!  04 Apr 2013 - R. Yantosca - Removed STT (now in State_Chm)
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! CLEANUP_TRACER begins here!
      !=================================================================
      IF ( ALLOCATED( ID_TRACER      ) ) DEALLOCATE( ID_TRACER      )
      IF ( ALLOCATED( ID_EMITTED     ) ) DEALLOCATE( ID_EMITTED     )
      IF ( ALLOCATED( TCVV           ) ) DEALLOCATE( TCVV           )
      IF ( ALLOCATED( TRACER_NAME    ) ) DEALLOCATE( TRACER_NAME    )
      IF ( ALLOCATED( TRACER_COEFF   ) ) DEALLOCATE( TRACER_COEFF   )
      IF ( ALLOCATED( TRACER_CONST   ) ) DEALLOCATE( TRACER_CONST   )
      IF ( ALLOCATED( TRACER_N_CONST ) ) DEALLOCATE( TRACER_N_CONST )
      IF ( ALLOCATED( TRACER_MW_G    ) ) DEALLOCATE( TRACER_MW_G    )
      IF ( ALLOCATED( TRACER_MW_KG   ) ) DEALLOCATE( TRACER_MW_KG   )
      IF ( ALLOCATED( XNUMOL         ) ) DEALLOCATE( XNUMOL         )

      END SUBROUTINE CLEANUP_TRACER
!EOC
      END MODULE TRACER_MOD
