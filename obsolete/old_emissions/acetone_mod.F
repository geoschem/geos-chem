!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: acetone_mod
!
! !DESCRIPTION: Module ACETONE\_MOD contains subroutines to emit the biogenic 
!  flux of acetone into the full chemistry simulation.
!\\
!\\
! !INTERFACE:
!
      MODULE ACETONE_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CLEANUP_ACETONE   
      PUBLIC :: EMISS_BIOACET
      PUBLIC :: OCEAN_SOURCE_ACET
      PUBLIC :: OCEAN_SINK_ACET
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Jacob, D.J., B.D. Field, E. Jin, I. Bey, Q. Li, J.A. Logan, and 
!        R.M. Yantosca, "Atmospheric budget of acetone", Geophys. Res. Lett., 
!        107(D11), 4100, 2002. 
!  (2 ) Nightingale et al [2000a], J. Geophys. Res, 14, 373-387
!  (3 ) Nightingale et al [2000b], Geophys. Res. Lett, 27, 2117-2120
!
! !REVISION HISTORY:
!  NOTES:
!  18 Sep 2001 - B. Field, R. Yantosca - Initial version  
!  (1 ) Added changes from bdf and updated comments (bmy, 9/5/01)
!  (2 ) Updated comments (bmy, 9/12/01)
!  (3 ) Removed VERBOSE flag and all "print-to-log-file" diagnostics.  The
!        ND11 diagnostic produces the same totals. (bdf, bmy, 9/18/01)
!  (4 ) Now cal GET_TAU0 w/ 3 arguments instead of 2.  Also minor bug
!        fix in READ_RESP (bmy, 11/15/01)
!  (5 ) Implemented fix for ocean source/sink from Mat Evans.  Also deleted 
!        obsolete code from 11/01. (bmy, 11/26/01)
!  (6 ) Eliminated more obsolete code from 11/01 (bmy, 2/27/02)
!  (7 ) Removed duplicate variable definitions (bmy, 3/20/02)
!  (8 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (9 ) Bug fix: Now apply true exponential loss in OCEAN_SINK_ACET, instead
!        of just the 1st order approximation. (bdf, bmy, 7/11/02)
!  (10) Scale the ocean source of acetone for GEOS-3 meteorology in order to
!        match the total listed in Jacob et al 2002. (bdf, bmy, 9/16/02)
!  (11) Now references "error_mod.f" (bmy, 10/15/02)
!  (12) Minor modifications to READ_JO1D, READ_RESP (bmy, 3/14/03)
!  (13) Add surface area scale factor for ocean source for 1x1 nested
!        grids.  (yxw, bmy, 5/16/03)
!  (14) Scale ACET ocean source to Jacob et al 2002 for GEOS-4, and now
!        account for surface area ratio for all GEOS grids. (bmy, 3/15/04)
!  (15) Now references "directory_mod.f" (bmy, 7/19/04)
!  (16) Now can read data from GEOS and GCAP grids.  Also now use Nightingale
!        et al 2000b formulation for piston velocity KL. (swu, bmy, 8/16/05)
!  (17) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (18) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (19) Updates for nested EU and NA grids (amv, bmy, 12/18/09)
!  (20) Updates for GEOS-4 1 x 1.25 grid (lok, bmy, 1/13/10)
!  13 Aug 2010 - R. Yantosca - Add modifications for MERRA (treat like GEOS-5)
!  04 Nov 2010 - R. Yantosca - Added ProTeX headers
!  06 Dec 2011 - E. Fischer  - Updated Ocean exchange, MEGAN biogenic emiss.
!                              Removed obsolete code.
!  19 Mar 2012 - M. Payer    - Removed obsolete subroutines READ_JO1D and
!                              READ_RESP (E. Fischer)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Avogadro's number
      REAL*8,  PARAMETER   :: AVO      = 6.022d23 

      ! Molecules C / kg C 
      REAL*8,  PARAMETER   :: XNUMOL_C = AVO / 12d-3 

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ocean_source_acet
!
! !DESCRIPTION: Subroutine OCEAN\_SOURCE\_ACET specifies the ocean source of 
!  acetone.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OCEAN_SOURCE_ACET( I, J, ACETONE, State_Met )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : CHECK_VALUE
      USE DIAG_MOD,           ONLY : AD11
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE TIME_MOD,           ONLY : GET_TS_EMIS

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND11
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)    :: I           ! GEOS-Chem longitude index
      INTEGER,        INTENT(IN)    :: J           ! GEOS-Chem latitude index
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*8,         INTENT(INOUT) :: ACETONE     ! Acetone emiss [atoms C/s]
!
! !REMARKS:

! 
! !REVISION HISTORY: 
!  14 Sep 2001 - B. Field    - Initial version
!  (1 ) Now compute u = SQRT( U10M^2 + V10M^2 ) as SQRT( SFCWINDSQR(I,J) ).
!        This is necessary since U10M and V10M are missing for 1996, and
!        need to be computed from UWND and VWND.  (bmy, 9/5/01)
!  (2 ) Bug fixes: multiply kg by 360000 and use exponent to the -0.5 power
!        in the expression for Kl.  Also update value of the OCEAN_SCALE
!        factor to 3.63e11.  Also updated comments. (bdf, bmy, 9/5/01)
!  (3 ) Bug fix: ACETONE has units of [atoms C/box/s], to match those of
!        EMISRR.  This involves an extra division by DTSRCE. (bmy, 9/14/01)
!  (4 ) Removed diagnostic variable OCEAN_SOURCE (bmy, 9/18/01)
!  (5 ) JO1D(IREF,JREF) is now JO1D(I,J).  Bug fix: Zero the ocean source
!        of acetone in grid boxes that are covered by less than 50% ocean.  
!        Bug fix: make sure -5 <= TC <= 30, in order to prevent the power
!        series for Schmidt # from going negative.  Also eliminate IREF,
!        JREF, we don't need them anymore. (mje, rvm, bmy, 11/26/01)
!  (6 ) Eliminated obsolete code from 11/01 (bmy, 2/27/02)
!  (7 ) Scale the ocean source of acetone for GEOS-3 meteorology in order to
!        match the total listed in Jacob et al 2002. (bdf, bmy, 9/16/02)
!  (8 ) Now use function GET_AREA_CM2 of "grid_mod.f" to return the
!        grid box area in cm2.  Use function GET_TS_EMIS from "time_mod.f".
!        Remove reference to CMN header file. (bmy, 2/11/03)
!  (9 ) Apply surface area scale factor for 1x1 nested grids, in order to
!        make the total ocean source the same as for 4x5. (yxw, bmy, 5/16/03)
!  (10) Scale the ocean source to Jacob et al 2002 for GEOS-4.  Also account
!        for surface area ratio for all GEOS grids. (bmy, 3/15/04)
!  (11) Added space in #ifdef block for GEOS-4 x 1x125 grid (bmy, 12/1/04)
!  (12) Now use Nightingale et al 2000b formulation for piston velocity KL.
!        (swu, bmy, 8/16/05)
!  (13) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (14) Adjust SCALE_FACTOR for 0.5 x 0.667 grid (dan, bmy, 11/6/08)
!  (15) Additional scale factors for NESTED_NA and NESTED_EU calculated and 
!        included (amv, bmy, 12/18/09)
!  (16) Added scale factor for GEOS-4 1 x 1.25 grid (lok, bmy, 1/13/10)
!  13 Aug 2010 - R. Yantosca - Treat MERRA in the same way as GEOS-5
!  04 Nov 2010 - R. Yantosca - Added ProTeX headers
!  04 Nov 2010 - R. Yantosca - Cleaned up #if statements for clarity
!  06 Dec 2011 - E. Fischer  - Updated ocean source and sink terms to be 
!                              different than Jacob et al. [2002]. Ocean mixed
!                              layer is now set to a constant concentration of
!                              acetone (15 nM). Fluxes are now calculated
!                              using a direct application of the standard
!                              two-film model described by Liss and Slater
!                              [1974]. The fluxes are calculated using an
!                              updated Henry's law coefficient and transfer
!                              velocities have been updated following Johnson
!                              [2010]. The model now reproduces aircraft
!                              measurements over the remote oceans well.
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8                 :: KG, U, TC, SC, kl, KKL, HSTAR, KL600
      REAL*8                 :: DTSRCE, OCEAN_ACET, AREA_CM2, FOCEAN
      REAL*8,  EXTERNAL      :: SFCWINDSQR
! evf(6/20/2011)
      REAL*8                 :: USTAR, CD, Cl
!
! !DEFINED PARAMETERS:
!
      ! ALPHA scale factor in kg*s^2/cm
      REAL*8,  PARAMETER     :: OCEAN_SCALE = 3.63d11     

      ! Coefficients for fitting the Schmidt number for acetone
      REAL*8,  PARAMETER     :: A0 =  3287.687d0
      REAL*8,  PARAMETER     :: A1 = -136.2176d0
      REAL*8,  PARAMETER     :: A2 =  2.20642d0
      REAL*8,  PARAMETER     :: A3 = -0.01410642d0


      !=================================================================
      ! OCEAN_SOURCE_ACET begins here!
      !=================================================================

      ! Emission timestep in seconds
      DTSRCE   = GET_TS_EMIS() * 60d0 

      ! Fraction of (I,J) that is ocean
      FOCEAN   = 1d0 - State_Met%FRCLND(I,J)    

      ! Area of grid box (I,J) in cm^2
      AREA_CM2 = GET_AREA_CM2( I, J, 1 )    

      !=================================================================
      ! Compute ocean source by Henry's law
      ! Now make sure only to compute the ocean source if more than
      ! 50% of the box is covered by water (mje, bdf, bmy, 11/26/01)
      !=================================================================
      IF ( FOCEAN > 0.5d0 ) THEN 

         ! Updated Henry's Law to 27 M/atm following Benkelberg et al. [1995],
         ! Johnson [2010], and Zhou and Mopper [1990](evf, 5/11/11)
         HSTAR = (660d0)*EXP(-5090d0*(1./298d0 - 1./State_Met%TS(I,J)))
  
         ! Want Henry constant exressed as
         ! concentration in air/concentration in water.
         ! Take reciprocal.(evf, 5/11/11) 
 
         HSTAR = 1d0/HSTAR

         ! Now HENCONST = dimensionless H 
         ! [mass Acetone/volume air]/[mass Acetone/volume H2O]

         ! Magnitude of resultant wind [m/s]
         ! SFCWINDSQR(I,J) is needed since this will compute the square
         ! of the surface wind correctly for all GEOS models (bmy, 9/5/01)
         U        = SQRT( SFCWINDSQR( I, J,
     &                                State_Met%U10M,
     &                                State_Met%V10M ) )

         ! TC is temperature in Celsius
         ! Also make sure -5 <= TC <= 30 (mje, rvm, bmy, 11/26/01)
         TC       = State_Met%TS(I,J) - 273.15d0                           
         TC       = MIN( MAX( TC, -5d0 ), 30d0 )

         ! SC is Schmidt # for acetone [unitless]
         SC       = A0 + A1*TC + A2*TC**2 + A3*TC**3 

         ! KL is conductance for mass transfer in liquid phase 
         ! (Nightingale et al 2000b), which has units of [cm/hr]
         KL       = ( 0.24d0*U*U + 0.061d0*U ) * SQRT( 600d0/Sc )  

         ! KG is conductance for mass transfer in gas phase (Asher 1997)
         ! Multiply KG by 360000 to convert from [m/s] to [cm/hr]

         ! Updated KG to the Johnson [2010] parameterization (evf, 5/13/2011)
         ! USTAR = friction velocity (U* in Johnson [2010]
         USTAR = SQRT(6.1d-4 + U*6.3d-5)*U
         ! CD = drag coefficient
         CD = (USTAR/U)**2

         !KG is airside transfer velocity (Johnson 2010)
         ! Multiply KG by 360000 to convert from [m/s] to [cm/hr]
         KG = 1d-3 + (USTAR/(13.3*SC**(1/2) + CD**(-1/2)-5+LOG(SC)/0.8))
         KG = KG * 360000d0

         ! KKL is the air-to-sea transfer velocity (Liss and Slater 1974)
         ! Multiply KKL by 3600 to convert from [cm/hr] to [cm/s]
         KKL      = 1d0 / ( 1d0/KL + 1d0/( HSTAR * KG ) ) 
         KKL      = KKL / 3600d0

         ! Turn off ocean uptake over snow/ice (ALBEDO > 0.4)
         IF ( State_Met%ALBD(I,J) > 0.4d0 ) KKL = 0d0
      
         !===============================================================
         !(evf, 5/11/11)
         ! Remove photochemical acetone source hypothesized by Jacob.
         ! et al. [2002]. 

         ! Assume a constant seawater concentration
         ! of Cl = 15 nM (Available measurements of acetone in seawater:
         ! Williams et al., (2004) GRL, VOL. 31, L23SO6, 
         ! doi:10.1029/2004GL020012
         ! Marandino et al., (2005) GRL, VOL 32, L15806,
         ! doi:10.1029/2005GL02385
         ! Kameyama et al., (2010) Marine Chemistry, VOL 122, 59-73
         ! Zhou and Mopper (1997) Marine Chemistry, VOL 56, 201-213

         ! convert Cl to kg acetone/cm3 (evf, 5/11/11)
         Cl = 15.0d-9*58.08d0/(1000.0d0*1000.0d0)

         !correct for the fraction of the grid cell that is ocean
         !and compute the flux ( kg/cm2/s)(evf, 5/11/11)
         OCEAN_ACET = Cl * KKL * FOCEAN 
          
         ! Convert to kg Acetone / box / step (evf, 5/11/11)
         OCEAN_ACET  = OCEAN_ACET * DTSRCE * AREA_CM2

         ! Convert to kg C / box / step (evf, 5/11/11)
         OCEAN_ACET  = OCEAN_ACET * 36d0/58.08d0
        

      ELSE

         ! If there is less than 50% water in the grid box, zero 
         ! the ocean source from acetone (mje, rvm, bmy, 11/26/01)
         OCEAN_ACET = 0d0
         
      ENDIF

      ! Add ocean source to total biogenic source in [atoms C/box/s]
      ACETONE = ACETONE + ( OCEAN_ACET * XNUMOL_C / DTSRCE )
            
      !=================================================================
      ! ND11 diag -- save ACETONE from the ocean in [atoms C/cm2/s]
      !=================================================================
      IF ( ND11 > 0 ) THEN
         AD11(I,J,4) = AD11(I,J,4) + ( OCEAN_ACET * XNUMOL_C ) / 
     &                               ( AREA_CM2   * DTSRCE   )
      ENDIF

      END SUBROUTINE OCEAN_SOURCE_ACET
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ocean_sink_acet
!
! !DESCRIPTION: Subroutine OCEAN\_SINK\_ACET applies the ocean sink to global
!  acetone concentrations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OCEAN_SINK_ACET( ACETONE, State_Met )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE CMN_DIAG_MOD
      USE DIAG_MOD,           ONLY : AD11
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL*8,         INTENT(INOUT) :: ACETONE(IIPAR,JJPAR) ! Acet mass [kg C]
! 
! !REVISION HISTORY: 
!  14 Sep 2001 - B. Field    - Initial version
!  (1 ) Remove references to CMN_UV10M and CMN_LWI -- these are now
!        obsolete in GEOS-CHEM versions 4.18 and higher (bmy, 9/5/01)
!  (2 ) Now compute u = SQRT( U10M^2 + V10M^2 ) as SQRT( SFCWINDSQR(I,J) ).
!        This is necessary since U10M and V10M are missing for 1996, and
!        need to be computed from UWND and VWND.  (bmy, 8/2/01)
!  (3 ) Now declare OCEANSINK_SCALE = 0.15 as a parameter.  This is the
!        optimized value of BETA from Emily Jin's analysis.  Also updated
!        comments. (bdf, bmy, 9/5/01)
!  (4 ) Updated comments.  Also parallellized DO loops. (bmy, 9/14/01)
!  (5 ) Removed diagnostic variable OCEAN_LOSS (bmy, 9/18/01)
!  (6 ) Bug fix: Zero the ocean sink of acetone in grid boxes where there
!        is less than 50% of ocean, and where there is ice on the surface.
!        Bug fix: Make sure -5 <= TC <= 30, in order to prevent the power
!        series for Schmidt # from going negative.  Also eliminate IREF,
!        JREF, we don't need them. (mje, rvm, bmy, 11/26/01)
!  (7 ) Eliminated obsolete code from 11/01 (bmy, 2/27/02)
!  (8 ) Bug fix: now use true exponential for loss instead of just 1st
!        order term.  Also added PRE_ACET variable to save previous acetone
!        mass for diagnostic, before applying loss.  (bdf, bmy, 7/11/02)
!  (9 ) Now use function GET_AREA_CM2 of "grid_mod.f" to return the
!        grid box area in cm2.  Now use function GET_TS_CHEM from
!        "time_mod.f".  Remove reference to CMN header file. (bmy, 2/11/03)
!  (12) Now use Nightingale et al 2000b formulation for piston velocity KL.
!        (swu, bmy, 8/16/05)
!  04 Nov 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  02 Aug 2012 - R. Yantosca - Add error trap for acetone under DEVEL tag
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I, IREF, J, JREF
      REAL*8            :: KH298, DHR, KH, U, TC, SC, KL, KG 
      REAL*8            :: KKL, CG, F, T1L, H, KL600, FLUX, HSTAR
      REAL*8            :: USTAR, CD
      REAL*8            :: AREA_CM2, DTCHEM, FOCEAN, OCEAN_ACET
      REAL*8            :: PRE_ACET
!
! !EXTERNAL FUNCTIONS:
!
      REAL*8, EXTERNAL  :: BOXVL, SFCWINDSQR
!
! !DEFINED PARAMETERS:
!
      ! Optimized value of BETA for ocean sink found by
      ! Emily Jin from inverse modeling analysis
      REAL*8, PARAMETER :: OCEANSINK_SCALE = 0.15d0

      ! Coefficients for fitting the Schmidt number for acetone
      REAL*8, PARAMETER :: A0 =  3287.687d0
      REAL*8, PARAMETER :: A1 = -136.2176d0
      REAL*8, PARAMETER :: A2 =  2.20642d0
      REAL*8, PARAMETER :: A3 = -0.01410642d0

      !=================================================================
      ! OCEAN_SINK_ACET begins here! 
      !
      ! Compute acetone lost to ocean sink and subtract from STT
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Loop over grid boxes
      DO J = 1, JJPAR

         DO I = 1, IIPAR

            ! Grid box area in cm2
            AREA_CM2 = GET_AREA_CM2( I, J, 1 )

            ! Fraction of grid box that is ocean
            FOCEAN = 1d0 - State_Met%FRCLND(I,J)

            !===========================================================
            ! Only compute ocean sink if there is more than 50% ocean
            ! in the grid box, and if it is not ice (albedo > 0.4)
            ! (mje, rvm, bmy, 11/26/01)
            !===========================================================
            IF ( FOCEAN > 0.5d0 .and.
     &           State_Met%ALBD(I,J) <= 0.4d0 ) THEN

               ! Updated Henry's Law to 27 M/atm following Benkelberg et al.
               ! [1995], Johnson [2010], and Zhou and Mopper [1990]
               ! (evf, 5/11/11)
               HSTAR = (660d0)*EXP(-5090d0*(1./298d0 -
     &                 1./State_Met%TS(I,J)))
  
               ! Want Henry constant exressed as
               ! concentration in air/concentration in water.
               ! Take reciprocal.(evf, 5/11/11)  
               HSTAR = 1d0/HSTAR

               ! Now HENCONST = dimensionless H 
               ! [mass Acetone/volume air]/[mass Acetone/volume H2O]       

               ! Magnitude of surface wind [m/s]
               ! SFCWINDSQR(I,J) is needed since this will compute the 
               ! square of the surface wind correctly for all GEOS models 
               U     = SQRT( SFCWINDSQR( I, J,
     &                                   State_Met%U10M,
     &                                   State_Met%V10M ) )

            
               ! TC is temperature in Celsius
               ! Bug fix: make sure -5 <= TC <= 30 (mje, rvm, bmy, 11/26/01) 
               TC    = State_Met%TS(I,J) - 273.15d0  
               TC    = MIN( MAX( TC, -5d0 ), 30d0 )

               ! SC is Schmidt # for acetone [unitless]
               SC    = A0 + A1*TC + A2*TC**2 + A3*TC**3 

               ! KL is conductance for mass transfer in liquid phase 
               ! (Nightingale et al 2000b), which has units of [cm/hr]
               KL    = ( 0.24d0*U*U + 0.061d0*U ) * SQRT( 600d0/Sc )  

               ! KG is conductance for mass transfer in gas phase (Asher 1997)
               ! Multiply KG by 360000 to convert from [m/s] to [cm/hr]
               ! Updated KG to the Johnson [2010] parameterization (evf, 5/13/2011)
               ! USTAR = friction velocity (U* in Johnson [2010]
               USTAR  = SQRT(6.1d-4 + U*6.3d-5)*U
               ! CD = drag coefficient
               CD = (USTAR/U)**2

               ! KG is airside transfer velocity (Johnson 2010)
               ! Multiply KG by 360000 to convert from [m/s] to [cm/hr]
               KG = 1d-3 + (USTAR/(13.3*SC**(1/2) + CD**(-1/2)-5+LOG(SC) 
     &         /0.8))
               KG = KG * 360000d0      

               ! KKL is the air-to-sea transfer velocity (Liss and Slater 1974)
               ! Multiply KKL by 3600 to convert from [cm/hr] to [cm/s]
               KKL   = 1d0 / ( 1d0/KL + 1d0/( HSTAR * KG ) ) 
               KKL   = KKL / 3600d0       

               ! CG is the gas concentration of acetone [kg C/cm3]
               CG    = ACETONE(I,J) / BOXVL(I,J,1,State_Met)

               ! FLUX is the air-to-sea flux of acetone in [kg C/cm2/s].
               FLUX  = KKL * CG / HSTAR    
               !========================================================
               ! Ocean loss of acetone consists of the following terms:
               !
               ! (1) FLUX, the air-to-sea flux of acetone in [kg C/cm2/s]
               !
               ! (2) AREA_CM2, the grid box surface area
               !
               ! (3) DTCHEM, the number of seconds per chemistry timestep
               !
               ! (4) FOCEAN, the fraction of the grid box that is ocean.
               !
               ! The units of the resultant ocean loss (in OCEAN_ACET) are:
               !
               !     kg C  | AREA_CM2 cm2 |  DTCHEM s           kg C
               !  ---------+--------------+--------------- = ------------
               !   cm2 * s |   grid box   | chem timestep     box * step
               !========================================================
               OCEAN_ACET = ( FLUX * AREA_CM2 * DTCHEM * FOCEAN )
               
            ELSE

               ! If there is less than 50% water in the grid box, or  
               ! if there is ice on the ocean, then zero the ocean sink
               ! for acetone (mje, rvm, bmy, 11/26/01)
               OCEAN_ACET = 0d0

            ENDIF

            ! Save mass of acetone in tmp variable for diagnostic
            PRE_ACET = ACETONE(I,J) 

            ! Apply exponential loss to acetone mass
#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
            ! Eric Nielsen added this error trap .  We can eventually
            ! bring this into the std G-C but we need to benchmark first.
            IF ( ABS( ACETONE(I,J) ) .gt. 1d-99 )
     &      ACETONE(I,J) = ACETONE(I,J) * EXP(-OCEAN_ACET/ACETONE(I,J))

#else
            ACETONE(I,J) = ACETONE(I,J) * EXP(-OCEAN_ACET/ACETONE(I,J))
#endif

            !===========================================================
            ! Diagnostics: save ACETONE lost to ocean in [atoms C/cm2/s]
            !===========================================================
            IF ( ND11 > 0 ) THEN
               AD11(I,J,5) = AD11(I,J,5) + 
     &              ( ( PRE_ACET - ACETONE(I,J) ) * XNUMOL_C ) / 
     &              ( AREA_CM2 * DTCHEM   )
            ENDIF
         ENDDO
      ENDDO

      END SUBROUTINE OCEAN_SINK_ACET
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_bioacet
!
! !DESCRIPTION: Subroutine EMISS\_BIOACET computes the biogenic emissions of 
!  ACETONE from monoterpenes, isoprene, methyl butenol, dry leaf matter, and 
!  grasslands.
!\\
!\\
! !INTERFACE:
!
!evf, edits to use MEGAN biogenic acetone emissions (5/25/2011)
      SUBROUTINE EMISS_BIOACET( I,    J,    TMMP,  EMMO, SUNCOS, Q_DIR,
     &                          Q_DIFF, XNUMOL_C, EMIS, EMMB, GRASS,
     &                          ACETONE )
!
! !USES:
!
      USE DIAG_MOD, ONLY : AD11
      USE GRID_MOD, ONLY : GET_AREA_M2
      USE GRID_MOD, ONLY : GET_XMID
      USE GRID_MOD, ONLY : GET_YMID
      USE TIME_MOD, ONLY : GET_TS_EMIS
 !(evf, 5/25/2011)
      USE MEGAN_MOD, ONLY : GET_EMACET_MEGAN

      USE CMN_SIZE_MOD                  ! Size parameters
      USE CMN_DIAG_MOD                  ! ND11
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)    :: I         ! Grid box longitude index
      INTEGER, INTENT(IN)    :: J         ! Grid box latitude index
      REAL*8,  INTENT(IN)    :: TMMP      ! Local Surface Air temperature [K]
      REAL*8,  INTENT(IN)    :: EMMO      ! Monoterpene emission [atoms C]
      REAL*8,  INTENT(IN)    :: EMIS      ! Isoprene emission [atoms C]
      REAL*8,  INTENT(IN)    :: EMMB      ! Methylbutenol emission  [atoms C]
      REAL*8,  INTENT(IN)    :: GRASS     ! Isoprene from grasslands [atoms C]
!evf, edits to use MEGAN biogenic acetone emissions (5/25/2011)
      REAL*8,  INTENT(IN)    :: SUNCOS    ! Cosine of Solar Zenith Angle
      REAL*8,  INTENT(IN)    :: Q_DIR     ! Flux of direct PAR above canopy
      REAL*8,  INTENT(IN)    :: Q_DIFF    ! Flux of diffuse PAR above canopy
      REAL*8,  INTENT(IN)    :: XNUMOL_C  ! Number of atoms C / kg C

!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(INOUT) :: ACETONE   ! Biogenic acetone [atoms C/s]
! 
! !REVISION HISTORY: 
!  18 Sep 2001 - B. Field    - Initial version
!  (1 ) Now pass acetone array (e.g. from STT) thru the argument list, since
!        this avoids dependence on IDTACET in this program (bmy, 8/1/01)
!  (2 ) Updated scale factors (bdf, bmy, 9/5/01)
!  (3 ) Updated comments (bmy, 9/14/01)
!  (4 ) Removed diagnostic variables: MONOTERPENES, ISOPRENE, ISOP_TOTAL,
!        MONO_TOTAL, NA_TOT, RESP_TOT, GRASS_TOT.  These have now been 
!        supplanted by the ND11 acetone source diagnostic. (bdf, bmy, 9/18/01)
!  (5 ) XRESP(I+I0,J+J0) is now XRESP(I,J) (bmy, 11/26/01)
!  (6 ) Eliminated obsolete code from 11/01 (bmy, 2/27/02)
!  (7 ) Removed duplicate definitions of EMMB and GRASS (bmy, 3/20/02)
!  (8 ) Now use functions from "grid_mod.f" to get surface area, lon, and
!        lat of grid box (I,J).   Use function GET_AREA_M2 to get the grid
!        box surface area in m2, then convert to cm2.   Now use function
!        GET_TS_EMIS from "time_mod.f".  Remove reference to CMN header
!        file. (bmy, 2/11/03)
!  04 Nov 2010 - R. Yantosca - Added ProTeX headers
!  06 Dec 2011 - E. Fischer  - Direct biogenic emissions of acetone from
!                              metabolism and decay are now calculated using
!                              the MEGAN biogenic emission model [Guenther et
!                              al., 2006]. The code assumes the fraction of
!                              emissions that are light-independent is 0.20,
!                              and the temperature response factor (beta) is
!                              0.10 as recommended by Alex Guenther. The
!                              dependency on GEIA has been removed.
!  08 Dec 2011 - M. Payer    - Remove use of obsolete CMN_MONOT_MOD
!  01 Mar 2012 - R. Yantosca - Use new grid routines from grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8                 :: EMMO_MOL,    YIELD_MO   
      REAL*8                 :: ACET_MOL,    ACET_MB    
      REAL*8                 :: ACET_MO,     DTSRCE,      AREA_CM2 
      REAL*8                 :: YIELD_ISOP,  MB_MOL 
      REAL*8                 :: ACET_C,      ACET_ISOP,   YIELD_GRASS 
      REAL*8                 :: ACET_GRASS,  ACETSCAL,    AREA_M2
      REAL*8                 :: X,           Y

      ! Scale factors for a posteriori
      REAL*8                 :: MONO_SCALE, DIRECT_SCALE, MB_SCALE 
      REAL*8                 :: DP_SCALE,   GRASS_SCALE
!evf, edits to use MEGAN biogenic acetone emissions (5/25/2011)
      REAL*8                 :: EMISS_MEGAN

      !=================================================================
      ! EMISS_BIOACET begins here!
      !
      ! The yield for acetone from biogenic sources comes from 
      ! experimental yields from monoterpenes, (Reissell et. al. 1999), 
      ! methyl butenol (Alvarado et. al. 1999) and emissions for 
      ! monoterpenes, methyl butenol, and acetone from Guenther et. al. 
      ! 1999.  Guenther's emissions are for North America, and have 
      ! been scaled to the entire globe
      !=================================================================

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

      ! Grid box areas in [m2] and [cm2]
      AREA_M2  = GET_AREA_M2( I, J, 1 )
      AREA_CM2 = AREA_M2 * 1d4

      !=================================================================
      ! (1) BIOGENIC EMISSIONS OF ACETONE FROM MONOTERPENES
      !
      ! Monoterpenes has same # molecules/kg of carbon as isoprene
      ! The yield for monoterpenes is .12 mol/mol from Reisell et.al. 
      ! 1999 (this does not includes direct acetone emissions)
      !=================================================================

      ! Convert [atoms C/box/step] to [molec MONOTERPENE/box/step]
      ! There are 10 C atoms per molecule of MONOTERPENE
      EMMO_MOL   = EMMO / 10d0      

      ! Apply yield from monoterpenes to get [molec ACET/box/step]
      YIELD_MO   = 0.116d0
      ACET_MOL   = EMMO_MOL * YIELD_MO

      ! Convert acetoneemissions back into [atoms C/box/step] 
      ACET_MO    = ACET_MOL * 3.d0

      ! Scale to a posteriori source from Jacob et al 2001 (bdf, 9/5/01)
      MONO_SCALE = 0.89d0
      ACET_MO    = ACET_MO * MONO_SCALE

      ! Convert monoterpene yield to [atoms C/box/s] and 
      ! add to the total biogenic acetone emissions
      ACETONE    = ACETONE + ( ACET_MO / DTSRCE )

      ! Diagnostics -- save ACETONE from MONOTERPENES in [atoms C/cm2/s]
      IF ( ND11 > 0 ) THEN
         AD11(I,J,1) = AD11(I,J,1) + 
     &                 ( ACET_MO / ( AREA_CM2 * DTSRCE ) )
      ENDIF
 
      !=================================================================
      ! (2) BIOGENIC ACETONE FROM METHYL BUTENOL -- NORTH AMERICA
      !
      ! Methyl Butenol (a.k.a. MBO) produces acetone with a molar yield 
      ! of 0.6 [Alvarado (1999)].  The biogenic source of MBO is thought 
      ! to be restricted to North America.  According to Guenther (1999) 
      ! North america emits 3.2Tg-C of MBO, producing 1.15 Tg-C of 
      ! Acetone in North America.
      !=================================================================
      ACET_MB = 0D0

      ! Lon and lat of grid box (I,J) in degrees
      X = GET_XMID( I, J, 1 )
      Y = GET_YMID( I, J, 1 )

      ! Methyl butenol is emitted only in North America, where
      ! ( -167.5 <= lon <= -52.5 ) and ( 16.0 <= lat <= 72.0 ) 
      IF ( ( X >= -167.5d0 .and. X <= -52.5d0 ) .AND. 
     &     ( Y >=   16.0d0 .and. Y <=  72.0d0 ) ) THEN

         ! Convert from [atoms C/box/step] to [molec MBO/box/step] 
         ! There are 5 C atoms per molecule MBO
         MB_MOL   = EMMB / 5.d0      

         ! Apply yield from MBO to get [molec ACET/box/step]
         MB_SCALE = 0.6d0                
         ACET_MOL = MB_MOL * MB_SCALE

         ! Convert from [molec ACET/box/step] to [atoms C/box/step]
         ! There are 3 C atoms per acetone molecule
         ACET_MB  = ACET_MOL * 3.d0      

         ! Scale to a posteriori source from Jacob et al 2001 (bdf, 9/5/01)
         MB_SCALE = 0.76d0 
         ACET_MB  = ACET_MB * MB_SCALE

         ! Convert MBO yield to [atoms C/box/s] and add 
         ! to the total biogenic acetone emissions
         ACETONE  = ACETONE + ( ACET_MB / DTSRCE )
            
         ! Diagnostics -- save ACETONE from MBO in [atoms C/cm2/s]
         IF ( ND11 > 0 ) THEN
            AD11(I,J,2) = AD11(I,J,2) + 
     &                   ( ACET_MB / ( AREA_CM2 * DTSRCE ) )
         ENDIF
      ENDIF

      !=================================================================
      ! (3) BIOGENIC ACETONE -- DIRECT EMISSION 
      ! evf, removed obsolete code, replaced with MEGAN acetone
      ! emissions (5/25/2011) Direct Emission now includes emission
      ! from grasses and emission from dry leaf matter
      ! 
      !=================================================================

      !evf, edits to use MEGAN biogenic acetone emissions (5/25/2011)
      ! Acetone Emissions from MEGAN (atoms C / box / step)         
      EMISS_MEGAN = GET_EMACET_MEGAN(I, J, SUNCOS, TMMP, 
     &              Q_DIR, Q_DIFF, XNUMOL_C )
     
      !convert MEGAN Acetone Emissions to atoms C/box/s (evf, 5/25/2011)
      EMISS_MEGAN = EMISS_MEGAN / DTSRCE

      ! Add MEGAN Acetone Emissions to total biogenic acetone emissions (evf, 5/25/2011)
      ACETONE     = ACETONE +  EMISS_MEGAN

      ! Save MEGAN Acetone Emissions to DIRECT EMISSIONS diag [atoms C/cm2/s]
      IF ( ND11 > 0 ) THEN
         AD11(I,J,3) = AD11(I,J,3) + 
     &                 ( EMISS_MEGAN / ( AREA_CM2 ) )
      ENDIF

      END SUBROUTINE EMISS_BIOACET
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_acetone
!
! !DESCRIPTION: Subroutine CLEANUP\_ACETONE deallocates module arrays
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_ACETONE
! 
! !REVISION HISTORY: 
!  14 Sep 2001 - R. Yantosca - Initial version
!  04 Nov 2010 - R. Yantosca - Added ProTeX headers
!  19 Mar 2012 - M. Payer    - Removed deallocation of JO1D and XRESP. They
!                              are obsolete after removal of READ_JO1D and
!                              READ_RESP routines.
!EOP
!------------------------------------------------------------------------------
!BOC

      END SUBROUTINE CLEANUP_ACETONE   
!EOC
      END MODULE ACETONE_MOD
