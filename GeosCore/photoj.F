!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: photoj
!
! !DESCRIPTION: Subroutine PHOTOJ is the driver routine for the FAST-J 
!  photolysis package.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE PHOTOJ( NLON,   NLAT,    YLAT,    DAY_OF_YR, 
     &                   MONTH,  DAY,     CSZA,    T,    
     &                   SA,     OD,      OPTDUST, OPTAER,
     &                   O3COL,  am_I_Root                    )
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : NDUST, NAER, NRH,  LLPAR
      USE CMN_FJ_MOD,   ONLY : JPMAX, JPPJ, JPNL, ZPJ
      USE JV_CMN_MOD,   ONLY : ZJ,    U0,   SZA,  SZAMAX

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)    :: NLON                   ! Grid box lon index
      INTEGER, INTENT(IN)    :: NLAT                   ! Grid box lat index 
      REAL*8,  INTENT(IN)    :: YLAT                   ! Latitude [degrees]
      INTEGER, INTENT(IN)    :: DAY_OF_YR              ! Day of year
      INTEGER, INTENT(IN)    :: MONTH                  ! Current month
      INTEGER, INTENT(IN)    :: DAY                    ! Day of month
      REAL*8,  INTENT(IN)    :: CSZA                   ! Cosine(SZA) [unitless]
      REAL*8,  INTENT(IN)    :: T(LLPAR)               ! Temperature [K]
      REAL*8,  INTENT(IN)    :: SA                     ! UV albedo [unitless]
      REAL*8,  INTENT(IN)    :: OD(LLPAR)              ! Visible OD [unitless]
      REAL*8,  INTENT(IN)    :: O3COL                  ! Overhead O3 column [DU]
      LOGICAL, INTENT(IN)    :: am_I_Root              ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(INOUT) :: OPTDUST(LLPAR,NDUST)   ! Dust OD [unitless]
      REAL*8,  INTENT(INOUT) :: OPTAER(LLPAR,NAER*NRH) ! Aerosol OD [unitless]
!
! !AUTHOR:
!  Oliver Wild & Michael Prather
!
! !REMARKS:
!  New FAST J-Value code, troposphere only (mjprather 6/96); uses special 
!  wavelength quadrature spectral data (jv_spec.dat) that includes only 
!  289 nm - 800 nm  (later a single 205 nm add-on); uses special compact Mie 
!  code based on Feautrier/Auer/Prather vers. 
!                                                                             .
!  Important variables from other modules:
!  (1) ZJ   : Column array for J-values 
!  (2) ZPJ  : Global array for J-values (passed to SMVGEAR)
!  (3) JPNL : # of GEOS-CHEM layers in which to compute J-values
!  (4) JPPJ : # of photolysis rxns for FAST-J
!                                                                             .
!  NOTE: The value of PI listed here is slightly different than the value
!  in CMN_GCTM_mod.F.  The last digit is 4, whereas in CMN_GCTM_mod.F, the
!  last digit is 3.  Keep for now during testing of grid-independent code, but
!  this may be something to revisit at a later data. (bmy, 3/6/12)
!
! !REVISION HISTORY: 
!  01 Jun 1996 - M. Prather & O. Wild - Initial version
!  (1 ) Renamed NSLON to NLON and NSLAT to NLAT.  Now add DAY_OF_YR 
!        (formerly IDAY) and DAY to the arg list.  Swap places in arg list 
!        of SA and OD.  Now pass NLON, NLAT, DAY_OF_YR and DAY to "set_prof.f".
!        Added standard documentation header; cosmetic changes. (bmy, 7/15/03)
!  (2 ) We don't need to pass "P" via the arg list (bmy, 2/13/07)
!  06 Mar 2012 - R. Yantosca - Now pass O3COL via the arg list
!  06 Mar 2012 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I, J
      REAL*8, PARAMETER :: PI = 3.14159265358979324D0

      !=================================================================
      ! PHOTOJ begins here!
      !=================================================================

      ! Zero ZJ (column J-value array) and ZPJ (global J-value array)
      DO I = 1, JPNL
      DO J = 1, JPPJ
         ZJ(I,J)            = 0.D0
         ZPJ(I,J,NLON,NLAT) = 0.D0
      ENDDO
      ENDDO

      ! Import the cosine of the SZA from the CTM (bmy, 9/10/99)
      U0  = CSZA
      SZA = ACOS(CSZA) * ( 180.0d0 / PI )

      !-----------------------------------------------------------------
      !### If you want to set SZA = 0 degrees for testing,
      !### then uncomment the following lines (bmy, 9/13/99) 
      !U0  = 1.0d0
      !SZA = 0.0d0
      !-----------------------------------------------------------------

      ! Return if sun is below the horizon
      IF ( SZA > SZAMAX ) RETURN

      ! Set up Air, O3, BC profiles on GEOS-CHEM vertical levels
      CALL SET_PROF( NLON, NLAT, YLAT, MONTH,   DAY, 
     &               T,    SA,   OD,   OPTDUST, OPTAER, O3COL )

      ! Compute actinic flux at each GEOS-CHEM vertical level
      CALL JVALUE( SA, am_I_Root )
      
      ! Calculate J-values for all species
      CALL JRATET( T, DAY_OF_YR )

      ! ZJ is the J-value array for this column only
      ! Store in ZPJ (global array) for passing to SMVGEAR
      DO I = 1, JPNL
      DO J = 1, JPPJ
         ZPJ(I,J,NLON,NLAT) = ZJ(I,J)
      ENDDO
      ENDDO

      END SUBROUTINE PHOTOJ
!EOC
