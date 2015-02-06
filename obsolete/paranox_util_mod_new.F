!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: paranox_mod
!
! !DESCRIPTION: Module PARANOX\_MOD contains subroutines for reading and 
! interpolating look up tables necessary for the PARANOX (PARAmeterization
! of emitted NOX) ship plume model developed by G.C.M. Vinken.
!
! The initial initial look up tables (LUT) distributed with GEOS-Chem v9-01-03
! used 7 input variables: Temperature, J(NO2), J(O1D), solar elevation angles
! at emission time and 5 hours later, and ambient concentrations of NOx
! and O3. This version was documented by  Vinken et al. (2011). Subsequently, 
! we added wind speed as an input variable. We also use J(OH) rather than J(O1D)
! to index the LUT (C. Holmes,  6 May 2013)
!
! The LUTs contain 3 quantities: 
!     FracNOx : The fraction of NOx emitted from ships that remains as NOx 
!               after 5 hours of plume aging. mol/mol
!     OPE     : Ozone production efficiency, mol(O3)/mol(HNO3)
!               The net production of O3 per mole of ship NOx oxidized over
!               5 hours of plume aging. Can be negative! 
!               Defined as OPE = [ P(O3) - L(O3) ] / P(HNO3), where each P
!               and L term is an integral over 5 hours. Net O3 production
!               in the plume is E(NOx) * (1-FracNOx) * OPE, where E(NOx) is
!               the emission rate of NOx from the ship (e.g. units: mol/s).
!     MOE     : Methane oxidation efficiency, mol(CH4)/mol(NOx)
!               The net oxidation of CH4 per mole of NOx emitted from ships
!               over 5 hours of plume aging.
!               Defined as MOE = L(CH4) / E(NOx).
! 
!\\
!\\
! !INTERFACE: 
!
      MODULE PARANOX_MOD
!
! !USES:
!
      USE HCO_ERROR_MOD
      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: READ_PARANOX_LUT
      PUBLIC  :: PARANOX_LUT
      PUBLIC  :: CLEANUP_PARANOX
!      PUBLIC  :: SHIPO3LOSS, SHIPHNO3DEP
! Old PARANOX version from GEOS-Chem v9-01-03 through v9-02, cdh 3/27/2014
!      PUBLIC  :: READ_PARANOX_LUT_v913 
!
! !MODULE VARIABLES
!
      ! Number of values for each variable in the look-up table
      INTEGER, PARAMETER ::  nT=4, nJ=4, nO3=4, nNOx=5, nSEA=12, nWS=5

      ! Reference values of variables in the look-up tables
      REAL*4  :: Tlev(nT),      JNO2lev(nJ),   O3lev(nO3)
      REAL*4  :: SEA0lev(nSEA), SEA5lev(nSEA)
      REAL*4  :: JRATIOlev(nJ), NOXlev(nNOx)
      REAL*4  :: WSlev(nWS)

      ! Look-up tables currently used in GEOS-Chem (likely in v10-01)
      ! Described by Holmes et al. (2014), now includes effects of wind speed
      ! Last two digits in LUT names indicate wind speed in m/s
      REAL*4, TARGET, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE ::
     &     FRACNOX_LUT02, FRACNOX_LUT06, FRACNOX_LUT10, 
     &     FRACNOX_LUT14, FRACNOX_LUT18,
     &     OPE_LUT02,     OPE_LUT06,     OPE_LUT10, 
     &     OPE_LUT14,     OPE_LUT18,
     &     MOE_LUT02,     MOE_LUT06,     MOE_LUT10, 
     &     MOE_LUT14,     MOE_LUT18,
     &     DNOX_LUT02,    DNOX_LUT06,    DNOX_LUT10, 
     &     DNOX_LUT14,    DNOX_LUT18
      
!      REAL*8, ALLOCATABLE :: SHIPO3LOSS(:,:), SHIPHNO3DEP(:,:)

      ! LEGACY VARIABLES
      ! Look-up tables used in GEOS-Chem v9-01-03 through v9-02
      ! Originally described by Vinken et al. (2011)
      ! These are no longer needed, but currently preserved
!      REAL*4, TARGET, ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:)  ::   
!     &     FRACNOX_LUT_V913, OPE_LUT_V913

!
! !REMARKS
!  References:
!  ============================================================================
!  (1 ) Vinken, G.C.M., Boersma, K.F., Jacob, D.J., and Meijer, E.W.:
!       Accounting for non-linear chemistry of ship plumes in the GEOS-Chem
!       global chemistry transport model, Atmos. Chem. Phys., 11, 11707-11722,
!       doi:10.5194/acp-11-11707-2011, 2011.
!  (2 ) Holmes, C.D., M.J. Prather, G.C.M. Vinken: The climate impact of
!       ship NOx emissions: uncertainties due to plume chemistry, 
!       Atmos. Chem. Phys. Discuss., 14, 3427-3458, 
!       doi:10.5194/acpd-14-3427-2014, 2014.
!
! !REVISION HISTORY:
!  06 Feb 2012 - M. Payer    - Initial version
!  01 Mar 2012 - R. Yantosca - Use updated GET_LOCALTIME from time_mod.F
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  03 Jun 2013 - C. Holmes   - Rewritten to include wind speed in the look-up
!                              table and to take input from netCDF
!  15 Jul 2014 - C. Holmes   - Make module variables allocatable, since they are 
!                              used only in full chemistry simulations.
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
! !IROUTINE: read_paranox_lut
!
! !DESCRIPTION: Subroutine READ\_PARANOX\_LUT reads look-up tables for use in
!  the PARANOX ship plume model (G.C.M. Vinken)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_PARANOX_LUT
!
! !USES:
!
!
! !REVISION HISTORY:
!  06 Feb 2012 - M. Payer    - Initial version modified from code provided by
!                              G.C.M. Vinken
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  03 Jun 2013 - C. Holmes   - Rewritten to include wind speed in the look-up 
!                              table and to take input from netCDF
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER             :: IOS
      CHARACTER(LEN=255)  :: FILENAME
      INTEGER             :: fID
      LOGICAL, SAVE       :: FIRST = .TRUE.

      !=================================================================
      ! READ_PARANOX_LUT begins here
      !=================================================================
      
      ! Initialize module variables on first use
      IF (FIRST) THEN
         CALL INIT_PARANOX
         FIRST=.FALSE.
      ENDIF

      ! Clear FILENAME
      FILENAME = ''

      ! FILENAME format string
 101  FORMAT( '/home/ckeller/geosChem/diesdas/Ship_NOx/',
     &     'ship_plume_lut_', I2.2, 'ms.nc' )
      
      ! Read 2m/s LUT
      WRITE( FILENAME, 101 ) 2
      CALL READ_LUT_NCFILE( TRIM( FILENAME ), 
     &     FRACNOX_LUT02, DNOx_LUT02, OPE_LUT02, MOE_LUT02, 
     &     Tlev, JNO2lev, O3lev, SEA0lev, SEA5lev, JRATIOlev, NOXlev )

      ! Read 6 m/s LUT
      WRITE( FILENAME, 101 ) 6
      CALL READ_LUT_NCFILE( TRIM( FILENAME ), 
     &     FRACNOX_LUT06, DNOx_LUT06, OPE_LUT06, MOE_LUT06 )

      ! Read 10 m/s LUT
      WRITE( FILENAME, 101 ) 10
      CALL READ_LUT_NCFILE( TRIM( FILENAME ), 
     &     FRACNOX_LUT10, DNOx_LUT10, OPE_LUT10, MOE_LUT10 )

      ! Read 14 m/s LUT
      WRITE( FILENAME, 101 ) 14
      CALL READ_LUT_NCFILE( TRIM( FILENAME ), 
     &     FRACNOX_LUT14, DNOx_LUT14, OPE_LUT14, MOE_LUT14 )

      ! Read 18 m/s LUT
      WRITE( FILENAME, 101 ) 18
      CALL READ_LUT_NCFILE( TRIM( FILENAME ), 
     &     FRACNOX_LUT18, DNOx_LUT18, OPE_LUT18, MOE_LUT18 )

      ! Wind speed levels correspond to the files that we just read
      WSlev = (/ 2, 6, 10, 14, 18 /)
      
      END SUBROUTINE READ_PARANOX_LUT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_lut_ncfile
!
! !DESCRIPTION: Subroutine READ\_LUT\_NCFILE reads look up tables for use in
!  the PARANOX ship plume model (C. Holmes)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_LUT_NCFILE( FILENAME, FNOX, DNOx, OPE, MOE,  
     &     T, JNO2, O3, SEA0, SEA5, JRATIO, NOX )
!
! !USES:
!
      ! Modules for netCDF read
      USE m_netcdf_io_open
      USE m_netcdf_io_get_dimlen
      USE m_netcdf_io_read
      USE m_netcdf_io_readattr
      USE m_netcdf_io_close
      
#     include "netcdf.inc"
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*),INTENT(IN)  :: FILENAME
!
! !OUTPUT PARAMETERS: 
!
      REAL*4, INTENT(OUT), DIMENSION(:,:,:,:,:,:,:) :: FNOX,OPE,MOE,DNOx
      REAL*4, INTENT(OUT), OPTIONAL :: T(:),      JNO2(:), O3(:) 
      REAL*4, INTENT(OUT), OPTIONAL :: SEA0(:),   SEA5(:)
      REAL*4, INTENT(OUT), OPTIONAL :: JRATIO(:), NOX(:)
!
! !REVISION HISTORY:
!  06 Feb 2012 - M. Payer    - Initial version modified from code provided by
!                              G.C.M. Vinken
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  03 Jun 2013 - C. Holmes   - Rewritten to include wind speed in the look-up 
!                              table and to take input from netCDF
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER             :: IOS
      INTEGER             :: fID

      !=================================================================
      ! READ_LUT_NCFILE begins here
      !=================================================================

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_PARANOX_LUT: Reading ', a )

      ! Open file for reading
      CALL Ncop_Rd( fId, FILENAME )

      !-----------------------------------------------------------------
      ! Read reference values used to construct the Look-up table
      !-----------------------------------------------------------------

      ! Temperature, K
      IF ( PRESENT(T) ) 
     &     CALL NcRd( T,      fId, 'T',      (/1/), shape(T) )
      
      ! J(NO2), 1/s
      IF ( PRESENT(JNO2) )
     &     CALL NcRd( JNO2,   fId, 'JNO2',   (/1/), shape(JNO2) )
      
      ! Ambient O3, ppb
      IF ( PRESENT(O3) )
     &     CALL NcRd( O3,     fId, 'O3',     (/1/), shape(O3) )
      
      ! Solar elevation angle at emission time t=0, deg
      IF ( PRESENT(SEA0) ) 
     &     CALL NcRd( SEA0,   fId, 'SEA0',   (/1/), shape(SEA0) )

      ! Solar elevation angle at time t=5h, deg
      IF ( PRESENT(SEA5) )
     &     CALL NcRd( SEA5,   fId, 'SEA5',   (/1/), shape(SEA5) )

      ! J(OH) / J(NO2) ratio, s/s
      IF ( PRESENT(JRATIO) ) 
     &     CALL NcRd( JRatio, fId, 'Jratio', (/1/), shape(JRatio) )

      ! Ambient NOx, ppt
      IF ( PRESENT(NOX) )
     &     CALL NcRd( NOX,    fId, 'NOx',    (/1/), shape(NOx) )

      !-----------------------------------------------------------------
      ! Read look up table for fraction of NOx remaining for ship
      ! emissions after 5 h [unitless]
      !-----------------------------------------------------------------

      CALL NcRd( FNOx, fId, 'FNOx', (/1,1,1,1,1,1,1/), 
     &     shape(FRACNOX_LUT02) )

!      PRINT*, "binary_fracnox: ", Fnox(1:4,1,1,2,3,4,4)

      !-----------------------------------------------------------------
      ! Read look up table for 5-h integrated Ozone Production Efficiency
      ! for ship emissions [molec O3 produced / molec NOx lost]
      !-----------------------------------------------------------------

      CALL NcRd( OPE, fId, 'OPE', (/1,1,1,1,1,1,1/), 
     &     shape(OPE_LUT02) )

!      PRINT*, "binary_intope: ", OPE(1:4,1,1,2,3,4,4)

      !-----------------------------------------------------------------
      ! Read look up table for 5-h integrated Methane Oxidation Efficiency
      ! for ship emissions [molec CH4 oxidized / molec NOx emitted]
      !-----------------------------------------------------------------

      CALL NcRd( MOE, fId, 'MOE', (/1,1,1,1,1,1,1/), 
     &     shape(MOE_LUT02) )

!      PRINT*, "binary_intmoe: ", MOE(1:4,1,1,2,3,4,4)

      !-----------------------------------------------------------------
      ! Read look up table for 5-h integrated NOx deposition fraction
      ! for ship emissions [molec NOx deposited / molec NOx emitted]
      !-----------------------------------------------------------------

      CALL NcRd( DNOx, fId, 'DNOx', (/1,1,1,1,1,1,1/), 
     &     shape(DNOx) )

!      PRINT*, "binary_depnox: ", DNOx(1:4,1,1,2,3,4,4)

      ! Close netCDF file
      CALL NcCl( fId )

      END SUBROUTINE READ_LUT_NCFILE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INTERPOL_LINWEIGHTS
!
! !DESCRIPTION:  Subroutine INTERPOL\_LINWEIGHTS finds the array elements and 
!      weights for piecewise 1-D linear interpolation. The input array of NODES
!      must be in monotonic ascending order. (C. Holmes 3/27/2014)
!
!      If Y is an array containing values of a function evaluated at the points 
!      given in NODES, then its interpolated value at the point VALUESIN will be
!      Y(VALUEIN) =  Y(INDICES(1)) * WEIGHTS(1) + 
!                    Y(INDICES(2)) * WEIGHTS(2)
!  
!      This subroutine finds indices of consecutive nodes that bracket VALUEIN and
!      weights such that
!      VALUEIN = NODES(INDICES(1))   * WEIGHTS(1)     + 
!                NODES(INDICES(1)+1) * (1-WEIGHTS(1))
!
!      For convenience, the returned values of INDICES and WEIGHTS are 2-element 
!      arrays, where
!          INDICES(2) = INDICES(1)+1 and
!          WEIGHTS(2) = 1 - WEIGHTS(1)
!     
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INTERPOL_LINWEIGHTS( NODES, VALUEIN, INDICES, WEIGHTS )
!
! !USES:
!
!
! !INPUT PARAMETERS: 
!
      REAL*4,INTENT(IN)   :: NODES(:), VALUEIN
!
! !OUTPUT PARAMETERS: 
!
      ! These arrays are always 2 elements each, but declaring 
      ! as deferred shape avoids array temporaries
      INTEGER,INTENT(OUT) :: INDICES(:) 
      REAL*4, INTENT(OUT) :: WEIGHTS(:) 

!
! !REVISION HISTORY:
!  03 Jun 2013 - C. Holmes      - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER :: I
      REAL*8  :: VALUE

      !=================================================================
      ! INTERPOL_LINWEIGHTS begins here!
      !=================================================================

      ! If larger than largest in LUT, assign largest level values
      VALUE = MIN( VALUEIN, MAXVAL( NODES ) )

      ! If smaller, assign smallest level value
      VALUE = MAX( VALUE,   MINVAL( NODES ) )

      ! Initialize
      INDICES = (/ 1, 1 /)

      ! Loop over interpolation nodes until we find the largest node value
      ! that is less than the desired value
      DO I=1, SIZE(NODES)
         INDICES(1) = I 
         IF ( VALUE <= NODES(I+1) ) EXIT
      END DO
      
      ! The next node
      INDICES(2) = INDICES(1) + 1

      ! Weights for the corresponding node indices
      WEIGHTS(1) = ( NODES(I+1) - VALUE ) / ( NODES(I+1) - NODES(I) )
      WEIGHTS(2) = 1.0 - WEIGHTS(1)

      END SUBROUTINE INTERPOL_LINWEIGHTS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: paranox_lut
!
! !DESCRIPTION:  Subroutine PARANOX\_LUT returns fractional remainder of
! ship NOx (FNOx), fraction of NOx that dry deposits as NOy species (DNOx), 
! ozone production efficiency (OPE), and methane oxidation 
! efficiency (MOE) after 5-hrs of plume aging. Values are taken taken from a 
! lookup table using piecewise linear interpolation. The look-up table is derived
! from the PARANOx gaussian plume model (Vinken et al. 2011; Holmes et al. 2014)
! (G.C.M. Vinken, KNMI, June 2010; C. Holmes June 2013)
!
! The lookup table uses 8 input variables:
!     TEMP   : model temperature, K
!     JNO2   : J(NO2) value, 1/s
!     O3     : concentration O3 in ambient air, ppb
!     SEA0   : solar elevation angle at emission time 5 hours ago, degree
!     SEA5   : solar elevation angle at this time, degree
!     JRatio : ratio J(OH)/J(NO2), unitless
!     NOx    : concentration NOx in ambient air, ppt
!     WS     : wind speed, m/s
!
! In GEOS-Chem v9-01-03 through v9-02, the effects of wind speed on FNOx and OPE
! were not included (wind speed set at 6 m/s). The JRatio also used J(O1D) rather 
! than J(OH); this has only a small effect on interpolated values.
! To reproduce the behavior of these earlier versions, modify code below marked 
! with ******* and call READ_PARANOX_LUT_v913 in emissions_mod.F
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE PARANOX_LUT(  am_I_Root, ExtState,  HcoState, 
     &                        I, J, FNOX, DNOx, OPE, RC, MOE_OUT ) 
!
! !USES:
!
!      USE TRACERID_MOD,       ONLY : IDO3,   IDNO,  IDNO2
!      USE TIME_MOD,           ONLY : GET_LOCALTIME
!      USE ERROR_MOD,          ONLY : ERROR_STOP, SAFE_DIV
!      USE COMODE_MOD,         ONLY : CSPEC, JLOP
!      USE COMODE_LOOP_MOD,    ONLY : NCS, NAMEGAS, JPHOTRAT, NRATES
!      USE COMODE_LOOP_MOD,    ONLY : IRM, IH2O
!      USE GIGC_State_Met_Mod, ONLY : MetState
!      USE CMN_SIZE_MOD
      USE HCO_STATE_MOD,        ONLY : HCO_State
      USE HCOX_STATE_MOD,       ONLY : Ext_State
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN)         :: am_I_Root
      TYPE(Ext_State), POINTER    :: ExtState
      TYPE(HCO_State), POINTER    :: HcoState
      INTEGER, INTENT(IN)         :: I, J        ! Grid indices
!
! !OUTPUT PARAMETERS: 
!
      REAL*8, INTENT(OUT)           :: FNOX ! fraction of NOx remaining, mol/mol
      REAL*8, INTENT(OUT)           :: DNOX ! fraction of NOx deposited, mol/mol
      REAL*8, INTENT(OUT)           :: OPE  ! net OPE, mol(net P(O3))/mol(P(HNO3))
      REAL*8, INTENT(OUT), OPTIONAL :: MOE_OUT ! net MOE, mol(L(CH4))/mol(E(NOx))
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER, INTENT(INOUT)        :: RC   ! Return code
!
! !REVISION HISTORY:
!     Jun 2010 - G.C.M. Vinken - Initial version
!  03 Jun 2013 - C. Holmes     - Heavily modified and simplified from previous LUT 
!                                interpolation code by G.C.M. Vinken and 
!                                M. Payer. LUT now includes wind speed.
!  04 Feb 2015 - C. Keller     - Updated for use in HEMCO.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!      INTEGER                    :: NK, JLOOP
      INTEGER                    :: I1,I2,I3,I4,I5,I6,I7,I8
      REAL*4                     :: FNOX_TMP, DNOX_TMP, OPE_TMP, MOE_TMP
      REAL*4                     :: WEIGHT, MOE
      REAL*8                     :: DENS, JNO2, JO1D, JOH, TAIR
      REAL(dp)                   :: H2O
      REAL(hp)                   :: AIR
!      CHARACTER*8                :: SPECNAME
!      REAL*8, EXTERNAL           :: FJFUNC

      ! Interpolation variables, indices, and weights
      REAL*4, DIMENSION(8)       :: VARS
      INTEGER,DIMENSION(8,2)     :: INDX
      REAL*4, DIMENSION(8,2)     :: WTS

      REAL*4,POINTER, DIMENSION(:,:,:,:,:,:,:) :: FRACNOX_LUT => NULL(),
     &     DNOX_LUT => NULL(), OPE_LUT  => NULL(),   MOE_LUT => NULL() 
    
      CHARACTER(LEN=255)         :: MSG
      CHARACTER(LEN=255)         :: LOC = 'PARANOX_LUT' 

      REAL(dp), PARAMETER        :: MWH2O = 18.0_dp

!      LOGICAL, SAVE              :: FIRST=.TRUE.
!      INTEGER, SAVE              :: Ind_JNO2, Ind_JO1D

      !=================================================================
      ! PARANOX_LUT begins here!
      !=================================================================

!      ! Index for CSPEC and AIRDENS arrays
!      JLOOP   = JLOP(I,J,1)

      ! Air mass, kg
      AIR  = ExtState%AIR%Arr%Val(I,J,1)

      ! Air density, molec/cm3
      !DENS = State_Met%AIRDEN(1,I,J) / 28.97d-3 * 6.02d23 / 1d6
      DENS = AIR                            * 1.d3
     &     / ExtState%AIRVOL%Arr%Val(I,J,1) / HcoState%Phys%AIRMW 
     &     * HcoState%Phys%Avgdr            / 1.d6

      ! Air temperature, K
      Tair = ExtState%T2M%Arr%Val(I,J)

!      !----------------------------------
!      ! Get J-Values for J(NO2) and J(OH)
!      !----------------------------------
!
!      IF (FIRST) THEN
!
!         ! Loop over photolysis reactions to find NO2, O3 photolysis reactions
!         ! Note: This assumes that the reaction numbers
!         ! do not change during a simulation
!         DO I1 = 1, JPHOTRAT(NCS)
!                        
!            ! Reaction number
!            NK  = NRATES(NCS) + I1
!                        
!            ! Name of species being photolyzed
!            SPECNAME = NAMEGAS(IRM(1,NK,NCS))
!
!            ! Check if this is NO2 or O3, 
!            ! store their index numbers
!            SELECT CASE ( TRIM( SPECNAME ) )
!            CASE ( 'NO2' )
!               Ind_JNO2 = I1
!            CASE ( 'O3' )
!               Ind_JO1D = I1
!            CASE DEFAULT
!            END SELECT
!                        
!         ENDDO
!
!         FIRST = .FALSE.
!
!      ENDIF

      ! Check if sun is up
      IF ( State_Met%SUNCOSmid(I,J) > 0d0 ) THEN

         ! J(NO2), 1/s
         !JNO2 = FJFUNC(I,J,1,Ind_JNO2,1,SPECNAME)
         JNO2 = ExtState%JNO2%Arr%Val(I,J)

         ! J(O1D), 1/s
         !JO1D = FJFUNC(I,J,1,Ind_JO1D,1,SPECNAME)
         JO1D = ExtState%JO1D%Arr%Val(I,J)

         ! H2O, molec/cm3. Get from specific humidity, which is in kg/kg.
         H2O = ExtState%SPHU%Arr%Val(I,J,1) * 1.d3
     &       * ExtState%AIR%Arr%Val(I,J,1)  / MWH2O
     &       * HcoState%Phys%Avgdr          / 1d6
         
         ! Calculate J(OH), the effective rate for O3+hv -> OH+OH,
         ! assuming steady state for O(1D).
         ! Rate coefficients are cm3/molec/s; concentrations are molec/cm3
         ! This should match the O3+hv (+H2O) -> OH+OH kinetics in calcrate.F
         JOH = JO1D *
!     $          1.63d-10 * EXP( 60.d0/Tair) * CSPEC(JLOOP,IH2O) /
!     $        ( 1.63d-10 * EXP( 60.d0/Tair) * CSPEC(JLOOP,IH2O) +
     $          1.63d-10 * EXP( 60.d0/Tair) * H2O /
     $        ( 1.63d-10 * EXP( 60.d0/Tair) * H2O +
     $          1.20d-10                    * DENS * 0.5000d-6  +
     $          2.15d-11 * EXP(110.d0/Tair) * DENS * 0.7808d0   +
     $          3.30d-11 * EXP( 55.d0/Tair) * DENS * 0.2095d0   )

      ELSE

         ! J-values are zero when sun is down
         JNO2 = 0d0
         JO1D = 0d0
         JOH  = 0d0

      ENDIF

      !========================================================================
      ! Load all variables into a single array
      !========================================================================

      ! Temperature, K
      VARS(1) = Tair

      ! J(NO2), 1/s
      VARS(2) = JNO2                           
      
      ! O3 concentration in ambient air, ppb
      VARS(3) = ExtState%O3%Arr%Val(I,J,1) / AIR
     &        * HcoState%Phys%AIRMW        / MW_O3
     &        * 1.E9
!      VARS(3) = CSPEC(JLOOP,IDO3) / DENS * 1.E9

      ! Solar elevation angle, degree
      ! SEA0 = SEA when emitted from ship, 5-h before current model time
      ! SEA5 = SEA at current model time, 5-h after emission from ship
      ! Note: Since SEA = 90 - SZA, then cos(SZA) = sin(SEA) and 
      ! thus SEA = arcsin( cos( SZA ) )
!      VARS(4) = ASIND( State_Met%SUNCOSmid5(I,J) )
!      VARS(5) = ASIND( State_Met%SUNCOSmid(I,J)  )
      VARS(4) = ASIND( SC5(I,J,SC5ID) )
      VARS(5) = ExtState%SUNCOSmid%Arr%Val(I,J)

      ! J(OH)/J(NO2), unitless
      ! Note J(OH) is the loss rate (1/s) of O3 to OH, which accounts for 
      ! the temperature, pressure and water vapor dependence of these reactions
      VARS(6) = SAFE_DIV( JOH, JNO2, 0D0 )

      ! NOx concetration in ambient air, ppt
!      VARS(7) = ( CSPEC(JLOOP,IDNO) + CSPEC(JLOOP,IDNO2) ) /DENS * 1.E12
      VARS(7) = ( ( ExtState%NO%Arr%Val(I,J,1)  / AIR
     &        *     HcoState%Phys%AIRMW         / MW_NO  )
     &        +   ( ExtState%NO2%Arr%Val(I,J,1) / AIR
     &        *     HcoState%Phys%AIRMW         / MW_NO2 ) )
     &        * 1.E12

      ! Wind speed, m/s
!      VARS(8) = SQRT( State_Met%U(I,J,1)**2 + State_Met%V(I,J,1)**2 )         
      VARS(8) = SQRT( ExtState%U10M%Arr%Val(I,J)**2 
     &        +       ExtState%V10M%Arr%Val(I,J)**2 )

      !*****************************************************
      ! Restoring the following lines reproduces the behavior of 
      ! GEOS-Chem v9-01-03 through v9-02
      ! the LUT was indexed with the ratio J(O1D)/J(NO2) (cdh, 3/27/2014)
!      VARS(6)   = SAFE_DIV( JO1D, JNO2, 0D0 )
!      JRATIOlev = (/ 5.e-4, 0.0015, 0.0025, 0.0055 /)
      !*****************************************************

      !========================================================================
      ! Find the indices of nodes and their corresponding weights for the 
      ! interpolation
      !========================================================================

      ! Temperature:
      CALL INTERPOL_LINWEIGHTS( Tlev, VARS(1), INDX(1,:), WTS(1,:) )
      
      ! J(NO2):
      CALL INTERPOL_LINWEIGHTS( JNO2lev, VARS(2), INDX(2,:), WTS(2,:) )

      ! [O3]:
      CALL INTERPOL_LINWEIGHTS( O3lev, VARS(3), INDX(3,:), WTS(3,:) )

      ! SEA0:
      CALL INTERPOL_LINWEIGHTS( SEA0lev, VARS(4), INDX(4,:), WTS(4,:) )

      ! SEA5:
      CALL INTERPOL_LINWEIGHTS( SEA5lev, VARS(5), INDX(5,:), WTS(5,:) )

      ! JRATIO:
      CALL INTERPOL_LINWEIGHTS( JRATIOlev, VARS(6), INDX(6,:), WTS(6,:))

      ! [NOx]:
      CALL INTERPOL_LINWEIGHTS( NOXlev, VARS(7), INDX(7,:), WTS(7,:) )

      ! Wind speed:
      CALL INTERPOL_LINWEIGHTS( WSlev, VARS(8), INDX(8,:), WTS(8,:) )

      !========================================================================
      ! Piecewise linear interpolation
      !========================================================================
      
      ! Initialize
      FNOX = 0.0
      DNOx = 0.0
      OPE  = 0.0
      MOE  = 0.0

      ! Loop over wind speed
      DO I8=1,2
         
         ! Point at the LUT for this wind speed
         ! Last two digits in fortran variable names indicate wind speed in m/s
         SELECT CASE ( NINT( WSlev(INDX(8,I8)) ) )
            CASE (  2 )
               FRACNOX_LUT => FRACNOX_LUT02
               DNOx_LUT    => DNOx_LUT02
               OPE_LUT     => OPE_LUT02
               MOE_LUT     => MOE_LUT02
            CASE (  6 )
               FRACNOX_LUT => FRACNOX_LUT06
               DNOx_LUT    => DNOx_LUT06
               OPE_LUT     => OPE_LUT06
               MOE_LUT     => MOE_LUT06
            CASE ( 10 )
               FRACNOX_LUT => FRACNOX_LUT10
               DNOx_LUT    => DNOx_LUT10
               OPE_LUT     => OPE_LUT10
               MOE_LUT     => MOE_LUT10
            CASE ( 14 )
               FRACNOX_LUT => FRACNOX_LUT14
               DNOx_LUT    => DNOx_LUT14
               OPE_LUT     => OPE_LUT14
               MOE_LUT     => MOE_LUT14
            CASE ( 18 ) 
               FRACNOX_LUT => FRACNOX_LUT18
               DNOx_LUT    => DNOx_LUT18
               OPE_LUT     => OPE_LUT18
               MOE_LUT     => MOE_LUT18
            CASE DEFAULT
                MSG = 'LUT error: Wind speed interpolation error!'
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN
         END SELECT

         !*****************************************************
         ! Restoring the following lines reproduces the behavior of 
         ! GEOS-Chem v9-01-03 through v9-02 in which wind speed
         ! effects on FNOx and OPE are neglected (cdh, 3/25/2014)
!         FRACNOX_LUT => FRACNOX_LUT_v913
!         OPE_LUT     => OPE_LUT_v913
         !*****************************************************

         ! loop over all other variables
         DO I7=1,2 
         DO I6=1,2
         DO I5=1,2
         DO I4=1,2
         DO I3=1,2
         DO I2=1,2
         DO I1=1,2

            !------------------------------------------
            ! Nodes and weights used in the interpolation
            !------------------------------------------

            ! Fraction NOx from the LUT
            FNOX_TMP = FRACNOX_LUT( INDX(1,I1), INDX(2,I2), INDX(3,I3), 
     &                              INDX(4,I4), INDX(5,I5), 
     &                              INDX(6,I6), INDX(7,I7) )

            DNOX_TMP = DNOx_LUT(    INDX(1,I1), INDX(2,I2), INDX(3,I3), 
     &                              INDX(4,I4), INDX(5,I5), 
     &                              INDX(6,I6), INDX(7,I7) )
            
            ! OPE from the LUT
            OPE_TMP  = OPE_LUT(     INDX(1,I1), INDX(2,I2), INDX(3,I3), 
     &                              INDX(4,I4), INDX(5,I5), 
     &                              INDX(6,I6), INDX(7,I7) )

            ! MOE from the LUT
            MOE_TMP  = MOE_LUT(     INDX(1,I1), INDX(2,I2), INDX(3,I3), 
     &                              INDX(4,I4), INDX(5,I5), 
     &                              INDX(6,I6), INDX(7,I7) )
            
            ! Interpolation weight for this element
            WEIGHT = WTS(1,I1) * WTS(2,I2) * WTS(3,I3) * WTS(4,I4) * 
     &               WTS(5,I5) * WTS(6,I6) * WTS(7,I7) * WTS(8,I8)

            !-----------------------------------
            ! Error Check
            !-----------------------------------

            !IF ENCOUNTER -999 IN THE LUT PRINT ERROR!!       
            IF ( ( FNOX_TMP < 0. ) .or. ( FNOX_TMP > 1. ) ) THEN
               
               PRINT*, 'PARANOX_LUT: fracnox = ,', fnox
            
               print*, I1, I2, I3, I4, I5, I6, I7, I8
               print*, INDX(1,I1), INDX(2,I2), INDX(3,I3),  INDX(4,I4), 
     &                 INDX(5,I5), INDX(6,I6), INDX(7,I7), INDX(8,I8)
               print*, VARS

               MSG = 'LUT error: Fracnox should be between 0 and 1!'
               CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
            ENDIF

            !-----------------------------------
            ! Final interpolated values
            !-----------------------------------

            ! Weighted sum of FNOx from the LUT
            FNOx = FNOx + FNOX_TMP * WEIGHT

            ! Weighted sum of DNOx from the LUT
            DNOx = DNOx + DNOX_TMP * WEIGHT

            ! Weighted sum of OPE from the LUT
            OPE  = OPE + OPE_TMP * WEIGHT

            ! Weighted sum of MOE from the LUT
            MOE  = MOE + MOE_TMP * WEIGHT

         END DO
         END DO
         END DO
         END DO
         END DO
         END DO
         END DO

      END DO

      ! Transfer MOE if optional output parameter is present
      IF ( PRESENT( MOE_OUT ) ) MOE_OUT = MOE

      ! Nullify pointers
      NULLIFY( FRACNOX_LUT )
      NULLIFY( DNOx_LUT )
      NULLIFY( OPE_LUT  )
      NULLIFY( MOE_LUT  )

      END SUBROUTINE PARANOX_LUT
!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: read_paranox_lut_v913
!!
!! !DESCRIPTION: Subroutine READ\_PARANOX\_LUT reads look up tables for use in
!!  the PARANOX ship plume model (G.C.M. Vinken)
!!\\
!!\\
!! !INTERFACE:
!!
!      SUBROUTINE READ_PARANOX_LUT_v913
!!
!! !USES:
!!
!      USE DIRECTORY_MOD, ONLY : DATA_DIR_1x1
!      USE FILE_MOD,      ONLY : IOERROR
!!
!! !REVISION HISTORY:
!!  06 Feb 2012 - M. Payer    - Initial version modified from code provided by
!!                              G.C.M. Vinken
!!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!      INTEGER             :: IOS
!      CHARACTER(LEN=255)  :: FILENAME
!      INTEGER             :: IU_FILE
!
!      !=================================================================
!      ! READ_PARANOX_LUT_v913 begins here
!      !=================================================================
!
!      !=================================================================
!      ! Read look up table for fraction of NOx remaining for ship
!      ! emissions [unitless]
!      !=================================================================
!
!      ! File name
!      FILENAME = TRIM( DATA_DIR_1x1 ) // 'PARANOX_201202/' //
!     &           'FracNOx_binary_5hrs_20gs.dat'
!
!      ! Echo info
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( 'READ_PARANOX_LUT_v913: Reading ', a )
!
!      ! Find a free file LUN
!      IU_FILE = findFreeLUN()
!
!      ! Open file to read  
!      OPEN( IU_FILE, FILE=TRIM( FILENAME ), FORM="binary", 
!     &      IOSTAT=IOS )
!      IF ( IOS /= 0 ) 
!     &     CALL IOERROR( IOS, IU_FILE, 'read_paranox_lut_v913:1' )
!
!      ! Read file
!      READ( IU_FILE, IOSTAT=IOS ) FRACNOX_LUT_v913
!      IF ( IOS /= 0 )
!     &     CALL IOERROR( IOS, IU_FILE, 'read_paranox_lut_v913:2' )
!
!      !PRINT*, "binary_fracnox: ", FRACNOX(1:4,1,1,2,3,4,4)      
!
!      ! Close file
!      CLOSE( IU_FILE )
!
!      !=================================================================
!      ! Read look up table for integrated Ozone Production Efficiency
!      ! for ship emissions [molec O3 produced / molec NOx lost]
!      !=================================================================
!
!      ! File name
!      FILENAME = TRIM( DATA_DIR_1x1 ) // 'PARANOX_201202/' //
!     &           'IntOPE_binary_5hrs_20gs.dat'
!
!      ! Echo info
!      WRITE( 6, 101 ) TRIM( FILENAME )
! 101  FORMAT( 'READ_PARANOX_LUT_v913: Reading ', a )
!
!      ! Find a free file LUN
!      IU_FILE = findFreeLUN()
!
!      ! Open file to read
!      OPEN( IU_FILE, FILE=TRIM( FILENAME ), FORM="BINARY" )
!      IF ( IOS /= 0 ) 
!     &     CALL IOERROR( IOS, IU_FILE, 'read_paranox_lut_v913:3' )
!
!      ! Read file
!      READ( IU_FILE, IOSTAT=IOS ) OPE_LUT_v913
!      IF ( IOS /= 0 ) 
!     &     CALL IOERROR( IOS, IU_FILE, 'read_paranox_lut_v913:4' )
!
!      PRINT*, "binary_intope: ", OPE_LUT_v913(1:4,1,1,2,3,4,4)
!
!      ! Close file
!      CLOSE( IU_FILE )
!
!      END SUBROUTINE READ_PARANOX_LUT_v913
!!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_paranox
!
! !DESCRIPTION: Initializes module variables for paranox_mod.F
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_PARANOX
!
! !REVISION HISTORY:
!  16 July 2014 - C. Holmes - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: AS
      CHARACTER(LEN=255) :: LOC = 'INIT_PARANOX'

      !=================================================================
      ! INIT_PARANOX begins here
      !=================================================================

      !--------------------------------
      ! FNOX
      !--------------------------------
      ALLOCATE( FRACNOX_LUT02(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRACNOX_LUT02' )
      FRACNOX_LUT02 = 0d0      

      ALLOCATE( FRACNOX_LUT06(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRACNOX_LUT06' )
      FRACNOX_LUT06 = 0d0      

      ALLOCATE( FRACNOX_LUT10(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRACNOX_LUT10' )
      FRACNOX_LUT10 = 0d0      

      ALLOCATE( FRACNOX_LUT14(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRACNOX_LUT014' )
      FRACNOX_LUT14 = 0d0      

      ALLOCATE( FRACNOX_LUT18(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRACNOX_LUT18' )
      FRACNOX_LUT18 = 0d0      


      !--------------------------------
      ! OPE
      !--------------------------------
      ALLOCATE( OPE_LUT02(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OPE_LUT02' )
      OPE_LUT02 = 0d0      

      ALLOCATE( OPE_LUT06(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OPE_LUT06' )
      OPE_LUT06 = 0d0      

      ALLOCATE( OPE_LUT10(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OPE_LUT10' )
      OPE_LUT10 = 0d0      

      ALLOCATE( OPE_LUT14(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OPE_LUT014' )
      OPE_LUT14 = 0d0      

      ALLOCATE( OPE_LUT18(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OPE_LUT18' )
      OPE_LUT18 = 0d0      


      !--------------------------------
      ! MOE
      !--------------------------------
      ALLOCATE( MOE_LUT02(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MOE_LUT02' )
      MOE_LUT02 = 0d0      

      ALLOCATE( MOE_LUT06(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MOE_LUT06' )
      MOE_LUT06 = 0d0      

      ALLOCATE( MOE_LUT10(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MOE_LUT10' )
      MOE_LUT10 = 0d0      

      ALLOCATE( MOE_LUT14(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MOE_LUT014' )
      MOE_LUT14 = 0d0      

      ALLOCATE( MOE_LUT18(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MOE_LUT18' )
      MOE_LUT18 = 0d0      


      !--------------------------------
      ! DNOx
      !--------------------------------
      ALLOCATE( DNOx_LUT02(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DNOx_LUT02' )
      DNOx_LUT02 = 0d0      

      ALLOCATE( DNOx_LUT06(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DNOx_LUT06' )
      DNOx_LUT06 = 0d0      

      ALLOCATE( DNOx_LUT10(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DNOx_LUT10' )
      DNOx_LUT10 = 0d0      

      ALLOCATE( DNOx_LUT14(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DNOx_LUT014' )
      DNOx_LUT14 = 0d0      

      ALLOCATE( DNOx_LUT18(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DNOx_LUT18' )
      DNOx_LUT18 = 0d0      

      
!      !--------------------------------
!      ! O3 loss and HNO3 deposition
!      !--------------------------------
!      ALLOCATE( SHIPO3LOSS( IIPAR, JJPAR ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SHIPO3LOSS' )
!      SHIPO3LOSS = 0d0      

!      ALLOCATE( SHIPHNO3DEP( IIPAR, JJPAR ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SHIPHNO3DEP' )
!      SHIPHNO3DEP = 0d0      


!      !--------------------------------
!      ! Legacy variables, can be deleted with READ_PARANOX_LUT_V913 is removed
!      !--------------------------------
!      ALLOCATE( FRACNOX_LUT_V913(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRACNOX_LUT_v913' )
!      FRACNOX_LUT_v913 = 0d0      
!
!      ALLOCATE( OPE_LUT_V913(nT,nJ,nO3,nSEA,nSEA,nJ,nNOx), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OPE_LUT_v913' )
!      OPE_LUT_v913 = 0d0      
 
      
      END SUBROUTINE INIT_PARANOX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_paranox
!
! !DESCRIPTION: Deallocates module variables for paranox_mod.F
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_PARANOX
!
! !USES:
!
!
! !REVISION HISTORY:
!  16 July 2014 - C. Holmes - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

!      IF ( ALLOCATED( SHIPHNO3DEP   ) ) DEALLOCATE( SHIPHNO3DEP   )
!      IF ( ALLOCATED( SHIPO3LOSS    ) ) DEALLOCATE( SHIPO3LOSS    )

      IF ( ALLOCATED( FRACNOX_LUT02 ) ) DEALLOCATE( FRACNOX_LUT02 )
      IF ( ALLOCATED( FRACNOX_LUT06 ) ) DEALLOCATE( FRACNOX_LUT06 )
      IF ( ALLOCATED( FRACNOX_LUT10 ) ) DEALLOCATE( FRACNOX_LUT10 )
      IF ( ALLOCATED( FRACNOX_LUT14 ) ) DEALLOCATE( FRACNOX_LUT14 )
      IF ( ALLOCATED( FRACNOX_LUT18 ) ) DEALLOCATE( FRACNOX_LUT18 )

      IF ( ALLOCATED( OPE_LUT02     ) ) DEALLOCATE( OPE_LUT02     )
      IF ( ALLOCATED( OPE_LUT06     ) ) DEALLOCATE( OPE_LUT06     )
      IF ( ALLOCATED( OPE_LUT10     ) ) DEALLOCATE( OPE_LUT10     )
      IF ( ALLOCATED( OPE_LUT14     ) ) DEALLOCATE( OPE_LUT14     )
      IF ( ALLOCATED( OPE_LUT18     ) ) DEALLOCATE( OPE_LUT18     )

      IF ( ALLOCATED( MOE_LUT02     ) ) DEALLOCATE( MOE_LUT02     )
      IF ( ALLOCATED( MOE_LUT06     ) ) DEALLOCATE( MOE_LUT06     )
      IF ( ALLOCATED( MOE_LUT10     ) ) DEALLOCATE( MOE_LUT10     )
      IF ( ALLOCATED( MOE_LUT14     ) ) DEALLOCATE( MOE_LUT14     )
      IF ( ALLOCATED( MOE_LUT18     ) ) DEALLOCATE( MOE_LUT18     )

      IF ( ALLOCATED( DNOx_LUT02    ) ) DEALLOCATE( DNOx_LUT02    )
      IF ( ALLOCATED( DNOx_LUT06    ) ) DEALLOCATE( DNOx_LUT06    )
      IF ( ALLOCATED( DNOx_LUT10    ) ) DEALLOCATE( DNOx_LUT10    )
      IF ( ALLOCATED( DNOx_LUT14    ) ) DEALLOCATE( DNOx_LUT14    )
      IF ( ALLOCATED( DNOx_LUT18    ) ) DEALLOCATE( DNOx_LUT18    )

      ! Legacy variables
      IF ( ALLOCATED( FRACNOX_LUT_v913 ) ) DEALLOCATE( FRACNOX_LUT_v913)
      IF ( ALLOCATED( OPE_LUT_v913     ) ) DEALLOCATE( OPE_LUT_v913    )

      END SUBROUTINE CLEANUP_PARANOX
!EOC
      END MODULE PARANOX_MOD
