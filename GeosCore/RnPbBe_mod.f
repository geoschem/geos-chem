! $Id: RnPbBe_mod.f,v 1.1 2009/09/16 14:06:46 bmy Exp $
      MODULE RnPbBe_MOD
!
!******************************************************************************
!  Module RnPbBe_MOD contains variables and routines used for the 
!  222Rn-210Pb-7Be simulation. (hyl, swu, bmy, 6/14/01, 8/4/06)
!
!  Module Variables:
!  ============================================================================
!  (1 ) LATSOU           : Array holding 10 latitudes for 7Be emissions
!  (2 ) PRESOU           : Array holding 33 pressure levels for 7Be emissions
!  (3 ) BESOU            : Array holding 7Be emissions for 10 lat x 33 prs levs
!  (4 ) XNUMOL_Rn        : Atoms 222Rn per kg 222Rn
!  (5 ) XNUMOL_Pb        : Atoms 210Pb per kg 210Pb
!  (6 ) XNUMOL_Be        : Atoms   7Be per kg   7Be
!
!  Module Procedures:
!  ============================================================================
!  (1 ) READ_7BE         : Reads Lal & Peters 7Be emissions from a file
!  (2 ) CORRECT_STE      : Corrects S-T exchange for 210Pb and 7Be
!  (3 ) EMISSRnPbBe      : Adds emissions of Rn, 210Pb, 7Be, to tracer array  
!  (4 ) CHEMRnPbBe       : Performs radioactive decay for Rn, 210Pb, 7Be
!  (5 ) SLQ              : Interpolation subroutine (cf. Numerical Recpies)
!
!  GEOS-CHEM modules referenced by RnPbBe_mod.f
!  ============================================================================
!  (1 ) dao_mod.f        : Module w/ arrays for DAO met fields
!  (2 ) diag_mod.f       : Module w/ GEOS-CHEM diagnostic arrays
!  (3 ) directory_mod.f  : Module w/ GEOS-CHEM data & met field dires
!  (4 ) file_mod.f       : Module w/ file unit numbers and error checks
!  (5 ) logical_mod.f    : Module w/ GEOS-CHEM logical switches
!  (6 ) tracer_mod.f     : Module w/ GEOS-CHEM tracer array STT etc.
!  (7 ) tropopause_mod.f : Module w/ routines to read in ann mean tropopause
!
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
!  NOTES:
!  (1 ) Added existing routines to this module (bmy, 6/14/01)
!  (2 ) Updated comments (bmy, 9/4/01)
!  (3 ) Eliminate AVGF; redimensioned XTRA2 (bmy, 9/25/01)
!  (4 ) Replace references to PW(I,J) with P(I,J) (bmy, 10/3/01)
!  (5 ) Remove obsolete code from 9/01 and 10/01 (bmy, 10/23/01)
!  (6 ) Removed duplicate variable declarations (bmy, 11/15/01)
!  (7 ) Now read files from DATA_DIR/RnPbBe_200203/ directory.  
!        Also updated comments. (bmy, 3/29/02)
!  (8 ) Incorporated latest changes from Hongyu Liu.  Also split off the
!        code to read in the 7Be emissions into a separate routine. 
!        Add parallel DO-loops in several places.  Cleaned up DRYFLXRnPbBe,
!        and now make sure ND44 accurately represents the drydep fluxes
!        of 210Pb and 7Be. (hyl, bmy, 8/7/02)
!  (9 ) Now reference AD from "dao_mod.f".  Now references "error_mod.f".
!        Moved routine DRYFLXRnPbBe into "drydep_mod.f".  (bmy, 1/27/03)
!  (10) Now references the new "time_mod.f" (bmy, 2/11/03)
!  (11) Bug fix in EMISSRnPbBe -- take abs( lat) for 7Be emiss. (bmy, 6/10/03)
!  (12) Bug fix in EMISSRnPbBe -- shut off 222Rn emissions in polar regions
!        (swu, bmy, 10/28/03)
!  (13) Now references "directory_mod.f", "logical_mod.f", and "tracer_mod.f"
!        (bmy, 7/20/04)
!  (14) Now modified for GCAP and GEOS-5 met fields (swu, bmy, 5/24/05)
!  (15) Now references "tropopause_mod.f"
!  (16) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!******************************************************************************
!
      IMPLICIT NONE 

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "RnPbBe_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE
       
      ! ... except these routines
      PUBLIC :: EMISSRnPbBe 
      PUBLIC :: CHEMRnPbBe

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      REAL*8            :: LATSOU(10), PRESOU(33), BESOU(10,33)
      REAL*8, PARAMETER :: XNUMOL_Rn = ( 6.0225d23 / 222.0d-3 )    
      REAL*8, PARAMETER :: XNUMOL_Pb = ( 6.0225d23 / 210.0d-3 )    
      REAL*8, PARAMETER :: XNUMOL_Be = ( 6.0225d23 /   7.0d-3 )

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
     
!------------------------------------------------------------------------------

      SUBROUTINE READ_7BE
!
!******************************************************************************
!  Subroutine READ_7BE reads the 7Be emissions from Lal & Peters on 33 
!  pressure levels.  This only needs to be done on the very first timestep.  
!  (hyl, bmy, 8/7/02, 7/19/04)
!
!  NOTES:
!  (1 ) This code was split off from routine EMISSRnPbBe below. (bmy, 8/7/02)
!  (2 ) Now reference DATA_DIR from "directory_mod.f" (bmy, 7/19/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IU_FILE, IOERROR

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      INTEGER            :: IOS, J, L
      CHARACTER(LEN=255) :: FILENAME

      !==============================================================
      ! READ_7BE begins here!
      !
      ! Units of 7Be emissions are [stars/g air/s].  
      ! Here, "stars" = # of nuclear disintegrations of cosmic rays
      !==============================================================

      ! Define the file name
      FILENAME = TRIM( DATA_DIR ) // 'RnPbBe_200203/7Be.Lal'

      ! Open the 7Be file
      OPEN( IU_FILE,      FILE=TRIM( FILENAME ), 
     &      STATUS='OLD', IOSTAT=IOS )
      
      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'emissRnPbBe:1' )
 
      ! Read latitudes in southern hemisphere
      READ ( IU_FILE, '(13X,F5.0,7F8.0)', IOSTAT=IOS ) 
     &     ( LATSOU(J), J=1,8 )

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'emissRnPbBe:2' )

      ! Add latitudes for 80S and 90S
      LATSOU(9)  = 80d0
      LATSOU(10) = 90d0

      ! For 33 levels read the pressure and the Be concentration
      ! at each of the above-defined southern latitudes
      DO L = 1, 33
         READ( IU_FILE, '(F5.0,8X,8F8.2)', IOSTAT=IOS ) 
     &        PRESOU(L), ( BESOU(J,L), J=1,8 )

         ! Error check
         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_FILE, 'emissRnPbBe:3' )
         ENDIF
      ENDDO

      ! Overwrite 70S at the top (as recommended by Koch 1996)
      BESOU(8,1) = 1900d0
      
      ! Copy value from 70S into 80S and 90S at all levels
      DO L = 1, 33
         BESOU(9,L)  = BESOU(8,L)
         BESOU(10,L) = BESOU(8,L)
      ENDDO

      ! All the numbers in the file need to be multiplied by 1e-5
      ! in order to put them into the correct data range. 
      BESOU = BESOU * 1d-5

      ! Close the file
      CLOSE( IU_FILE )

      ! Return to calling program
      END SUBROUTINE READ_7BE
               
!------------------------------------------------------------------------------

      SUBROUTINE CORRECT_STE( EMISSION )
!
!******************************************************************************
!  Subroutine CORRECT_STE reduces the emission of 210Pb and/or 7Be in the
!  stratosphere, to correct for too fast STE in the GEOS-CHEM model.
!  (hyl, bmy, 8/7/02, 8/4/06)
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1  ) EMISSION (REAL*8)  : Emissions to be corrected [kg]
!
!  NOTES:
!  (1 ) Now updated for GCAP met fields (swu, bmy, 5/24/05)
!  (2 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!******************************************************************************
!
#     include "define.h"    ! Switches

      ! Arguments
      REAL*8, INTENT(INOUT) :: EMISSION
      
      !=================================================================
      ! CORRECT_STE begins here!
      !
      ! Correction factors were computed by Hongyu Liu (hyl, 8/6/02)
      !=================================================================
#if   defined( GEOS_3 )
      EMISSION = EMISSION / 3.5d0 

#elif defined( GEOS_4 )
      !EMISSION = 0d0           ! to be determined later

#elif defined( GEOS_5 )
      !EMISSION = 0d0           ! to be determined later

#elif defined( GCAP )
      EMISSION = EMISSION / 3.5d0

#endif
      
      ! Return to calling program
      END SUBROUTINE CORRECT_STE

!------------------------------------------------------------------------------

      SUBROUTINE EMISSRnPbBe
!
!******************************************************************************
!  Subroutine EMISSRnPbBe emits 222Rn and 7Be into the tracer array STT.
!  (hyl, bey, bmy, 5/28/99, 10/28/03)
!
!  NOTES:
!  (1 ) Also added Hongyu's code for emission of Be7 (bmy, 3/22/99)
!  (2 ) Now trap I/O errors with subroutine IOERROR (bmy, 5/28/99)
!  (3 ) Eliminate obsolete code and ND63 diagnostic (bmy, 4/12/00)
!  (4 ) Now reference TS from "dao_mod.f" instead of from common block
!        header file "CMN_TS". (bmy, 6/23/00)
!  (5 ) Cosmetic changes (bmy, 7/12/00)
!  (6 ) Now use IOS /= 0 criterion to trap both I/O errors and EOF 
!        condition. (bmy, 9/13/00)
!  (7 ) Added to module "RnPbBe_mod.f".  Also updated comments and made
!        cosmetic changes. (bmy, 6/14/01)
!  (8 ) Replace PW(I,J) with P(I,J) (bmy, 10/3/01)
!  (9 ) Now reference DATA_DIR from "CMN_SETUP".  Added FILENAME variable.
!        Now read "7Be.Lal" file from DATA_DIR/RnPbBe_200203/ directory.
!        (bmy, 3/29/02)
!  (10) Add diagnostics for Rn/Be emissions.  Also cleaned up some old code
!        and added parallel DO-loops.  Correct for S-T exchange for 7Be
!        emissions. Updated comments, cosmetic changes. (hyl, 8/6/02)
!  (11) Now reference routine GET_PCENTER from "pressure_mod.f", which
!        returns the correct "floating" pressure. (dsa, bdf, bmy, 8/20/02)
!  (12) Now reference AD from "dao_mod.f".  Now make FIRSTEMISS a local SAVEd 
!        variable instead of an argument.  (bmy, 1/27/03)
!  (13) Now use routine GET_YMID from "grid_mod.f" instead of common block
!        variable YLMID.  Now replace DXYP(JREF)*1d4 with routine GET_AREA_CM2
!        of "grid_mod.f".  Now use routine GET_TS_EMIS from time_mod.
!        (bmy, 2/11/03)
!  (14) Bug fix: take the absolute value of latitude -- this was a bug when
!        implementing the GET_YMID function from v5-04. (bmy, 6/10/03)
!  (15) Now reference GET_YEDGE from "grid_mod.f".  
!  (16) Bug fix: the Rn emission in antarctic area in the original code would 
!        lead to enormously hight Rn concentrations there, esp. after boundary
!        layer mixing.  Now apply different emissions over land and water,
!        and also shut off emissions poleward of 70 deg. (swu, bmy, 10/28/03)
!  (17) Now reference LEMIS from "logical_mod.f".  Now reference STT and
!        N_TRACERS from "tracer_mod.f" (bmy, 7/20/04)
!  (18) Remove reference to CMN; it's obsolete.  Now use inquiry functions
!        from "tropopause_mod.f" to diagnose strat boxes. (bmy, 8/15/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,        ONLY : AD, TS
      USE DIAG_MOD,       ONLY : AD01 
      USE GRID_MOD,       ONLY : GET_AREA_CM2, GET_YMID, GET_YEDGE 
      USE LOGICAL_MOD,    ONLY : LEMIS
      USE TIME_MOD,       ONLY : GET_TS_EMIS
      USE TRACER_MOD,     ONLY : STT, N_TRACERS
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_STRAT
      USE PRESSURE_MOD,   ONLY : GET_PCENTER

#     include "CMN_SIZE" ! Size parameters
#     include "CMN_DIAG" ! ND02 
#     include "CMN_DEP"  ! FRCLND

      ! Local variables
      LOGICAL, SAVE      :: FIRSTEMISS = .TRUE.
      INTEGER            :: I,          J,          L,         N
      REAL*8             :: A_CM2,      ADD_Be,     ADD_Rn,    Rn_LAND
      REAL*8             :: Rn_WATER,   DTSRCE,     LAT_TMP,   P_TMP
      REAL*8             :: Be_TMP,     Rn_TMP,     LAT_S,     LAT_N
      REAL*8             :: LAT_H,      LAT_L,      F_LAND,    F_WATER
      REAL*8             :: F_BELOW_70, F_BELOW_60, F_ABOVE_60

      !=================================================================
      ! EMISSRnPbBe begins here!
      !=================================================================

      ! Return if we are not doing emissions!
      IF ( .not. LEMIS ) RETURN   

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

      !=================================================================
      ! Add 222Rn emissions into tracer #1 according to the following:
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
      !=================================================================
      
      ! Loop over latitudes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, LAT_S, LAT_N, LAT_H, LAT_L, F_BELOW_70     )
!$OMP+PRIVATE( F_BELOW_60, F_ABOVE_60, A_CM2, Rn_LAND, Rn_WATER )
!$OMP+PRIVATE( F_LAND, F_WATER, ADD_Rn                          )
      DO J = 1, JJPAR
         
         ! Get ABS( latitude ) at S and N edges of grid box
         LAT_S      = ABS( GET_YEDGE(J)   ) 
         LAT_N      = ABS( GET_YEDGE(J+1) )
         LAT_H      = MAX( LAT_S, LAT_N )
         LAT_L      = MIN( LAT_S, LAT_N ) 

         ! Fraction of grid box w/ ABS( latitude ) less than 70 degrees
         F_BELOW_70 = ( 70.0d0 - LAT_L ) / ( LAT_H - LAT_L )

         ! Fraction of grid box w/ ABS( latitude ) less than 60 degrees
         F_BELOW_60 = ( 60.0d0 - LAT_L ) / ( LAT_H - LAT_L )

         ! Fraction of grid box w/ ABS( latitude ) greater than 60 degrees
         F_ABOVE_60 = 1d0 - F_BELOW_60

         ! Grid box surface area [cm2]
         A_CM2      = GET_AREA_CM2( J )

         ! Baseline 222Rn emissions over land [kg]
         ! Rn_LAND [kg] = [1 atom 222Rn/cm2/s] / [atoms/kg] * [s] * [cm2]
         Rn_LAND    = 1d0 / XNUMOL_Rn * DTSRCE * A_CM2 

         ! Baseline 222Rn emissions over water or ice [kg]
         Rn_WATER   = Rn_LAND * 0.005d0

      ! Loop over longitudes
      DO I = 1, IIPAR

         ! Fraction of grid box that is land
         F_LAND  = FRCLND(I,J)

         ! Fraction of grid box that is water
         F_WATER = 1d0 - F_LAND

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

               ADD_Rn = 
                        ! Consider 222Rn emissions equatorward of 
                        ! 60 degrees for both land (1.0 [atoms/cm2/s]) 
                        ! and water (0.005 [atoms/cm2/s])
     &                  F_BELOW_60 * 
     &                  ( Rn_LAND  * F_LAND  ) + 
     &                  ( Rn_WATER * F_WATER ) +

                        ! If the grid box straddles the 60 degree boundary
                        ! then also consider the emissions poleward of 60
                        ! degrees.  222Rn emissions here are 0.005 [at/cm2/s].
     &                  F_ABOVE_60 * Rn_WATER


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
         IF ( TS(I,J) < 273.15 ) ADD_Rn = ADD_Rn / 3d0

         ! Save 222Rn into STT array [kg]
         STT(I,J,1,1) = STT(I,J,1,1) + ADD_Rn

         ! ND01 diag: 222Rn emission [kg/s] 
         IF ( ND01 > 0 ) THEN
            AD01(I,J,1,1) = AD01(I,J,1,1) + ( ADD_Rn / DTSRCE )
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Add 7Be emissions into tracer #3 (if necessary)
      !
      ! Original units of 7Be emissions are [stars/g air/sec],
      ! where "stars" = # of nuclear disintegrations of cosmic rays
      !=================================================================
      IF ( N_TRACERS >= 3 ) THEN

         ! Read 7Be emissions on the first timestep only
         IF ( FIRSTEMISS ) CALL READ_7BE

         !==============================================================
         ! Now interpolate from 33 std levels onto GEOS-CHEM levels 
         !==============================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, LAT_TMP, P_TMP, Be_TMP, ADD_Be )
!$OMP+SCHEDULE( DYNAMIC )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Get absolute value of latitude, since we will assume that 
            ! the 7Be distribution is symmetric about the equator
            LAT_TMP = ABS( GET_YMID( J ) )

            ! Pressure at (I,J,L) -- need to change for fvDAS!
            P_TMP   = GET_PCENTER( I, J, L )
                 
            ! Interpolate 7Be [stars/g air/sec] to GEOS-CHEM levels
            CALL SLQ( LATSOU,PRESOU,BESOU,10,33,LAT_TMP,P_TMP,Be_TMP)

            ! Be_TMP = [stars/g air/s] * [0.045 atom/star] * 
            !          [kg air] * [1e3 g/kg] = 7Be emissions [atoms/s]
            Be_TMP  = Be_TMP * 0.045d0 * AD(I,J,L) * 1.d3 
                  
            ! ADD_Be = [atoms/s] * [s] / [atom/kg] = 7Be emissions [kg]
            ADD_Be  = Be_TMP * DTSRCE / XNUMOL_Be 

            ! Correct the strat-trop exchange of 7Be
            IF ( ITS_IN_THE_STRAT( I, J, L ) ) THEN
               CALL CORRECT_STE( ADD_Be )
            ENDIF

            ! Add 7Be into STT tracer array [kg]
            STT(I,J,L,3) = STT(I,J,L,3) + ADD_Be

            ! ND01 diag: 7Be emission [kg/s]
            IF ( ND01 > 0 ) THEN
               AD01(I,J,L,3) = AD01(I,J,L,3) + ( ADD_Be / DTSRCE )
            ENDIF
         ENDDO 
         ENDDO 
         ENDDO
!$OMP END PARALLEL DO
      ENDIF
     
      ! Reset FIRSTEMISS
      FIRSTEMISS = .FALSE.

      ! Return to calling program
      END SUBROUTINE EMISSRnPbBe     

!------------------------------------------------------------------------------
      
      SUBROUTINE CHEMRnPbBe
!
!******************************************************************************
!  Subroutine CHEMRnPbBe performs loss chemistry on 222Rn, 210Pb, and 7Be.
!  (hyl, amf, bey, bmy, 10/13/99, 8/15/05)
!
!  NOTES:
!  (1 ) Now use F90 syntax (bmy, hyl, 3/22/99)
!  (2 ) Add FIRSTCHEM as an argument.  Only compute the exponential terms
!        when FIRSTCHEM = .TRUE., and save the values for later use
!        (bmy, 3/24/99)
!  (3 ) Cosmetic changes (bmy, 10/13/99)
!  (4 ) Eliminate obsolete code and ND63 diagnostic (bmy, 4/12/00)
!  (5 ) Cosmetic changes (bmy, 7/12/00)
!  (6 ) Added to module "RnPbBe_mod.f".  Also updated comments 
!        and made cosmetic changes. (bmy, 6/14/01)
!  (7 ) Add diagnostics for Rn/Be emissions.  Also cleaned up some old code
!        and added parallel DO-loops.  Updated comments. (hyl, 8/6/02)
!  (8 ) Now make FIRSTCHEM a local SAVEd variable.  (bmy, 1/27/03)
!  (9 ) Now use function GET_TS_CHEM from "time_mod.f" (bmy, 2/11/03)
!  (10) Now references STT and N_TRACERS from "tracer_mod.f" (bmy, 7/20/04)
!  (11) Remove reference to CMN; it's obsolete.  Now use inquiry functions 
!        from "tropopause_mod.f" to diagnose strat boxes. (bmy, 8/15/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,       ONLY : AD01, AD02
      USE TIME_MOD,       ONLY : GET_TS_CHEM
      USE TRACER_MOD,     ONLY : STT, N_TRACERS
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_STRAT

#     include "CMN_SIZE" ! Size parameters
#     include "CMN_DIAG" ! ND01, ND02

      ! Local variables
      LOGICAL, SAVE     :: FIRSTCHEM = .TRUE.
      INTEGER           :: I, J, L, N            
      REAL*8            :: ADD_Pb, Be_LOST ,DTCHEM, Pb_LOST 
      REAL*8            :: Rn_LOST(IIPAR,JJPAR,LLPAR)

      ! Static variables
      REAL*8, SAVE      :: EXP_Rn, EXP_Pb, EXP_Be

      ! Ratio of molecular weights of 210Pb/222Rn
      REAL*8, PARAMETER :: Pb_Rn_RATIO = 210d0 / 222d0

      !=================================================================
      ! CHEMRnPbBe begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Pre-compute exponential terms only on first timestep
      IF ( FIRSTCHEM ) THEN 
         
         ! Fraction of (222Rn, 210Pb, 7Be) left after radioactive decay
         EXP_Rn = EXP( -DTCHEM * 2.097d-6  )
         EXP_Pb = EXP( -DTCHEM * 9.725d-10 ) 
         EXP_Be = EXP( -DTCHEM * 1.506d-7  )

         ! Reset FIRSTCHEM flag
         FIRSTCHEM = .FALSE.
      ENDIF

      !=================================================================
      ! Radioactive decay of 222Rn (tracer #1)
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Rn_LOST = amount of 222Rn lost to decay [kg]
         Rn_LOST(I,J,L) = STT(I,J,L,1) * ( 1d0 - EXP_Rn )

         ! ND02 diag: 222Rn lost to decay [kg/s]
         IF ( ND02 > 0 ) THEN
            AD02(I,J,L,1) = AD02(I,J,L,1) + ( Rn_LOST(I,J,L) / DTCHEM )
         ENDIF

         ! Subtract Rn_LOST from STT [kg]
         STT(I,J,L,1) = STT(I,J,L,1) - Rn_LOST(I,J,L)
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Radioactive decay of 210Pb (tracer #2)
      !=================================================================
      IF ( N_TRACERS >= 2 ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, ADD_Pb, Pb_LOST )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
           
            ! ADD_Pb = Amount of 210Pb gained by decay from 222Rn [kg]
            ADD_Pb = Rn_LOST(I,J,L) * Pb_Rn_RATIO 

            ! Correct strat-trop exchange of 210Pb in stratosphere
            IF ( ITS_IN_THE_STRAT( I, J, L ) ) THEN
               CALL CORRECT_STE( ADD_Pb )
            ENDIF

            ! ND01 diag: 210Pb emission from 222Rn decay [kg/s]
            IF ( ND01 > 0 ) THEN
               AD01(I,J,L,2) = AD01(I,J,L,2) + ( ADD_Pb / DTCHEM )
            ENDIF

            ! Add 210Pb gained by decay from 222Rn into STT [kg]
            STT(I,J,L,2) = STT(I,J,L,2) + ADD_Pb          

            ! Amount of 210Pb lost to radioactive decay [kg]
            ! NOTE: we've already added in the 210Pb gained from 222Rn
            Pb_LOST = STT(I,J,L,2) * ( 1d0 - EXP_Pb )

            ! ND02 diag: 210Pb lost to decay [kg/s]
            IF ( ND02 > 0 ) THEN
               AD02(I,J,L,2) = AD02(I,J,L,2) + ( Pb_LOST / DTCHEM )
            ENDIF

            ! Subtract 210Pb lost to decay from STT [kg]
            STT(I,J,L,2) = STT(I,J,L,2) - Pb_LOST
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !=================================================================
      ! Radioactive decay of 7Be (tracer #3)
      !=================================================================
      IF ( N_TRACERS >= 3 ) THEN 

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, Be_LOST )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Amount of 7Be lost to decay [kg]
            Be_LOST = STT(I,J,L,3) * ( 1d0 - EXP_Be )

            ! ND02 diag: 7Be lost to decay [kg/s]
            IF ( ND02 > 0 ) THEN
               AD02(I,J,L,3) = AD02(I,J,L,3) + ( Be_LOST / DTCHEM )
            ENDIF

            ! Subtract amount of 7Be lost to decay from STT [kg]
            STT(I,J,L,3) = STT(I,J,L,3) - Be_LOST
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Return to calling program  
      END SUBROUTINE CHEMRnPbBe

!------------------------------------------------------------------------------

      SUBROUTINE SLQ( X, Y, Z, N, M, U, V, W )
!
!******************************************************************************
!  Subroutine SLQ is an interpolation subroutine from a Chinese 
!  reference book (says Hongyu). (hyl, bmy, 3/17/98, 11/15/01)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) X (REAL*8) : X-axis coordinate on original grid
!  (2 ) Y (REAL*8) : Y-axis coordinate on original grid
!  (3 ) Z (REAL*8) : Array of data on original grid
!  (4 ) N (REAL*8) : First dimension of Z
!  (5 ) M (REAL*8) : Second dimension of Z
!  (6 ) U (REAL*8) : X-axis coordinate for desired interpolated value
!  (7 ) V (REAL*8) : Y-axis coordinate for desired interpolated value
!
!  Arguments as Output:
!  ============================================================================
!  (8 ) W (REAL*8) : Interpolated value of Z array, at coordinates (U,V) 
!
!  NOTES:
!  (1 ) Added to "RnPbBe_mod.f" (bmy, 7/16/01)
!  (2 ) Removed duplicate definition of IQ.  Added comments. (bmy, 11/15/01)
!******************************************************************************
!
      ! Arguments
      INTEGER   :: N, M 
      REAL*8    :: X, Y, Z, U, V, W, B, HH
      DIMENSION :: X(N), Y(M), Z(N,M), B(3)

      ! Local variables
      INTEGER   NN, IP, I, J, L, IQ, K, MM

      !=================================================================
      ! SLQ begins here!
      !=================================================================
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
 10      IF (IABS(I-J).NE.1) THEN
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
 20      IF (IABS(J-I).NE.1) THEN
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
 30         CONTINUE
            B(I)=B(I)+HH
 40      CONTINUE
 50   CONTINUE
      W=0.0
      DO 70 I=1,NN
         HH=B(I)
         DO 60 J=1,NN
            IF (J.NE.I) THEN
               HH=HH*(U-X(IP+J-1))/(X(IP+I-1)-X(IP+J-1))
            END IF
 60      CONTINUE
         W=W+HH
 70   CONTINUE

      ! Return to calling program
      END SUBROUTINE SLQ

!------------------------------------------------------------------------------

      END MODULE RnPbBe_MOD


