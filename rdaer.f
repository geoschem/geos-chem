! $Id: rdaer.f,v 1.1 2003/06/30 20:26:01 bmy Exp $
      SUBROUTINE RDAER( THISMONTH, THISYEAR, SO4_NH4_NIT )
!
!******************************************************************************
!  Subroutine RDAER reads global aerosol concentrations as determined by
!  Mian Chin.  Calculates optical depth at each level for "set_prof.f".
!  Also calculates surface area for heterogeneous chemistry.  It uses aerosol
!  parameters in FAST-J input file "jv_spec.dat" for these calculations.
!  (rvm, bmy, 11/04/01, 3/27/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISMONTH (INTEGER) : Number of the current month (1-12)
!  (2 ) THISYEAR  (INTEGER) : 4-digit year value (e.g. 1997, 2002)
!
!  NOTES:
!  (1 ) At the point in which "rdaer.f" is called, ABSHUM is actually
!        absolute humidity and not relative humidity (rvm, bmy, 2/28/02)
!  (2 ) Now force double-precision arithmetic by using the "D" exponent.
!        (bmy, 2/28/02)
!  (3 ) At present aerosol growth is capped at 90% RH.  The data
!        in jv_spec.dat could be used to allow a particle to grow to
!        99% RH if desired. (rvm, 3/15/02)
!  (4 ) Bug fix: TEMP2 needs to be sized (IIPAR,JJPAR,LLPAR) (bmy, 5/30/02)
!  (5 ) Now reference BXHEIGHT from "dao_mod.f".  Also references ERROR_STOP
!        from "error_mod.f".  Delete local declaration of TIME, since that
!        is also declared w/in comode.h -- this causes compile-time errors
!        on the ALPHA platform. (gcc, bmy, 11/6/02)
!  (6 ) Now use the online SO4, NH4, NIT aerosol, taken from the STT array, 
!        and passed via SO4_NH4_NIT argument if sulfate chemistry is turned on.
!        Otherwise, read monthly mean sulfate from disk.  (rjp, bmy, 3/23/03)
!  (7 ) Now call READ_BPCH2 with QUIET=.TRUE., which prevents info from being
!        printed to stdout.  Also made cosmetic changes. (bmy, 3/27/03)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE COMODE_MOD,   ONLY : ABSHUM, ERADIUS, IXSAVE, 
     &                         IYSAVE, IZSAVE,  TAREA 
      USE DAO_MOD,      ONLY : BXHEIGHT
      USE DIAG_MOD,     ONLY : AD21
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE TRANSFER_MOD, ONLY : TRANSFER_3D

      IMPLICIT NONE

#     include "cmn_fj.h"   ! LPAR, CMN_SIZE
#     include "jv_cmn.h"   ! ODAER, QAA, RAA
#     include "CMN_DIAG"   ! ND21, LD21
#     include "CMN_SETUP"  ! DATA_DIR
#     include "comode.h"   ! NTLOOP

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH, THISYEAR
      REAL*8,  INTENT(IN) :: SO4_NH4_NIT(IIPAR,JJPAR,LLPAR)

      ! Local variables
      CHARACTER(LEN=255)  :: FILENAME
      INTEGER             :: I, J, L, N, R, JLOOP, IRH, IRHN
      INTEGER, SAVE       :: MONTH_LAST = -999
      REAL*4              :: TEMP(IGLOB,JGLOB,LGLOB)
      REAL*8              :: TEMP2(IIPAR,JJPAR,LLPAR)
      REAL*8              :: MSDENS(NAER), XTAU, DRYAREA

      ! Mass of hydrophobic aerosol from Mian Chin
      REAL*8, SAVE        :: DAERSL(IIPAR,JJPAR,LLPAR,2)       

      ! Mass of hydrophilic aerosol from Mian Chin
      REAL*8, SAVE        :: WAERSL(IIPAR,JJPAR,LLPAR,NAER)    

      ! Fraction of aerosol from H2O
      REAL*8		  :: FWET      

      ! Effective radius at RH bins read in from "jv_spec.dat"
      REAL*8		  :: RW(NRH)	

      ! Effective radius at RH after interpolation
      REAL*8		  :: REFF       

      ! Q at different RH bins read in from "jv_spec.dat"
      REAL*8		  :: QW(NRH)	
 
      ! Used to interpolate between sizes
      REAL*8		  :: FRAC       
 
      ! Change in Q (extinction efficiency)
      REAL*8		  :: SCALEQ     

      ! Change in Radius with RH
      REAL*8		  :: SCALER     

      ! Chnge in Optical Depth vs RH
      REAL*8		  :: SCALEOD(IIPAR,JJPAR,LLPAR,NRH) 

      ! Change in Vol vs RH 
      REAL*8		  :: SCALEVOL(IIPAR,JJPAR,LLPAR)  

      ! Relative Humidities
      REAL*8,  SAVE       :: RH(NRH)   = (/0d0,0.5d0,0.7d0,0.8d0,0.9d0/)

      ! Index to aerosol types in jv_spec.dat
      ! The following are ordered according to the mass densities below
      INTEGER, SAVE	  :: IND(NAER) = (/22, 29, 36, 43, 50/)

      !=================================================================
      ! RDAER begins here!
      !
      ! Read aerosol data from the binary punch file during the first 
      ! chemistry timestep and, after that, at the start of each month.
      !=================================================================
      IF ( THISMONTH /= MONTH_LAST ) THEN   
         
         ! Save the current month
         MONTH_LAST = THISMONTH

#if   defined( GEOS_STRAT )

         ! Select proper dust file name
         ! Get TAU0 value used to index the punch file
         SELECT CASE ( THISYEAR )
           
            ! GEOS-STRAT -- 1996 dust fields from P. Ginoux
            ! Use 1996 fields as a proxy for December 1995
            CASE ( 1995, 1996 )
               FILENAME = TRIM( DATA_DIR )    // 
     &              'aerosol_200106/aerosol.' //
     &              GET_NAME_EXT() // '.'     // 
     &              GET_RES_EXT()  // '.1996'

               XTAU = GET_TAU0( THISMONTH, 1, 1996 )

            ! GEOS-STRAT -- 1997 dust fields from P. Ginoux
            CASE ( 1997 )
               FILENAME = TRIM( DATA_DIR )    // 
     &              'aerosol_200106/aerosol.' //
     &              GET_NAME_EXT() // '.'     // 
     &              GET_RES_EXT()  // '.1997'
               
               XTAU = GET_TAU0( THISMONTH, 1, 1997 )
           
            ! 1995, 1996, 1997 are the only valid GEOS-STRAT years
            CASE DEFAULT
               CALL ERROR_STOP( 'Invalid GEOS-STRAT year!', 'rdaer.f' )

         END SELECT

#elif defined( GEOS_1 )

         ! Filename for GEOS-1
         FILENAME = TRIM( DATA_DIR ) // 'aerosol_200106/aerosol.' // 
     &              GET_NAME_EXT()   // '.' // GET_RES_EXT()

         ! Use the "generic" year 1990
         XTAU = GET_TAU0( THISMONTH, 1, 1990 )

#else 

         ! Filename for GEOS-3 or GEOS-4
         FILENAME = TRIM( DATA_DIR ) // 'aerosol_200106/aerosol.' // 
     &              GET_NAME_EXT()   // '.' // GET_RES_EXT()

         ! Use the "generic" year 1996
         XTAU = GET_TAU0( THISMONTH, 1, 1996 )

#endif

         ! Print filename
         WRITE( 6, 110 ) TRIM( FILENAME )
 110     FORMAT( '     - RDAER: Reading ', a )

         !==============================================================
         ! Read aerosol concentrations [kg/m3] for each type from the 
         ! binary punch file.  Also merge the coarse sea salt aerosols 
         ! into a combined bin rather than carrying them separately.
         !==============================================================

         !-------------------
         ! Sulfate
         !-------------------
         IF ( LSULF ) THEN 

            ! For sulfate chemistry, save lumped SO4, NH4, NIT into
            ! WAERSL for both FAST-J and hetchem (rjp, bmy, 3/23/03)
            WAERSL(:,:,:,1) = SO4_NH4_NIT

         ELSE
            
            ! For sulfate chemistry turned off, then read monthly
            ! mean sulfate aerosol from disk (rjp, bmy, 3/23/03)
            CALL READ_BPCH2( FILENAME, 'ARSL-L=$', 1,     XTAU,
     &                       IGLOB,     JGLOB,     LGLOB, TEMP )

            CALL TRANSFER_3D( TEMP, WAERSL(:,:,:,1) )
            
         ENDIF

         !-------------------
         ! Hydrophobic BC
         !-------------------
         CALL READ_BPCH2( FILENAME, 'ARSL-L=$', 2,     
     &                    XTAU,      IGLOB,     JGLOB,     
     &                    LGLOB,     TEMP,      QUIET=.TRUE. )

         CALL TRANSFER_3D( TEMP, DAERSL(:,:,:,1) )

         !-------------------
         ! Hydrophilic BC
         !-------------------
         CALL READ_BPCH2( FILENAME, 'ARSL-L=$', 3,     
     &                    XTAU,      IGLOB,     JGLOB,     
     &                    LGLOB,     TEMP,      QUIET=.TRUE. )

         CALL TRANSFER_3D( TEMP, WAERSL(:,:,:,2) )

         !-------------------
         ! Hydrophobic OC
         !-------------------
         CALL READ_BPCH2( FILENAME, 'ARSL-L=$', 4,     
     &                    XTAU,      IGLOB,     JGLOB,     
     &                    LGLOB,     TEMP,      QUIET=.TRUE. )

         CALL TRANSFER_3D( TEMP, DAERSL(:,:,:,2) )

         !-------------------
         ! Hydrophilic OC
         !-------------------
         CALL READ_BPCH2( FILENAME, 'ARSL-L=$', 5,     
     &                    XTAU,      IGLOB,     JGLOB,     
     &                    LGLOB,     TEMP,      QUIET=.TRUE. )

         CALL TRANSFER_3D( TEMP, WAERSL(:,:,:,3) )

         !-------------------
         ! Sea Salt (accum)
         !-------------------
         CALL READ_BPCH2( FILENAME, 'ARSL-L=$', 6,     
     &                    XTAU,      IGLOB,     JGLOB,     
     &                    LGLOB,     TEMP,      QUIET=.TRUE. )

         CALL TRANSFER_3D( TEMP, WAERSL(:,:,:,4) )

         !-------------------
         ! Sea Salt (coarse)
         !-------------------
         CALL READ_BPCH2( FILENAME, 'ARSL-L=$', 7,     
     &                    XTAU,      IGLOB,     JGLOB,     
     &                    LGLOB,     TEMP,      QUIET=.TRUE. )

         CALL TRANSFER_3D( TEMP, WAERSL(:,:,:,5) )

         !-------------------
         ! Sea Salt (coarse)
         !-------------------
         CALL READ_BPCH2( FILENAME, 'ARSL-L=$', 8,     
     &                    XTAU,      IGLOB,     JGLOB,     
     &                    LGLOB,     TEMP,      QUIET=.TRUE. )

         CALL TRANSFER_3D( TEMP, TEMP2 )

         ! Accumulate into one size bin
         WAERSL(:,:,:,5) = WAERSL(:,:,:,5) + TEMP2 

         !-------------------
         ! Sea Salt (coarse)
         !-------------------
         CALL READ_BPCH2( FILENAME, 'ARSL-L=$', 9,     
     &                    XTAU,      IGLOB,     JGLOB,     
     &                    LGLOB,     TEMP,      QUIET=.TRUE. )

         CALL TRANSFER_3D( TEMP, TEMP2 )

         ! Accumulate into one size bin
         WAERSL(:,:,:,5) = WAERSL(:,:,:,5) + TEMP2 

      ENDIF 

      !=================================================================
      ! Calculate optical depth and surface area at each timestep
      ! to account for the change in relative humidity
      !
      ! For the optical depth calculation, this involves carrying the 
      ! optical depth at each RH as separate aerosols since OPMIE.f 
      ! treats the phase functions and single scattering albedos 
      ! separately. (An alternative would be to rewrite OPMIE.f)
      !
      ! Scaling is sufficient for the surface area calculation
      !=================================================================
      MSDENS(1) = 1700.0    !SO4
      MSDENS(2) = 1000.0    !BC 
      MSDENS(3) = 1800.0    !OC 
      MSDENS(4) = 2200.0    !SS (accum)
      MSDENS(5) = 2200.0    !SS (coarse)

      ! Loop over types of aerosol
      DO N = 1, NAER

         ! Zero array
         SCALEOD(:,:,:,:) = 0d0
         
         !==============================================================
         ! Determine aerosol growth rates from the relative 
         ! humidity in each box
         !
         ! The optical depth scales with the radius and Q alone
         ! since SCALEDENS cancels as follows
         ! 
         !    SCALER 	= RW / RDRY
         !    SCALEDENS = DENSWET / DENSDRY
         !    SCALEM 	= SCALEDENS * SCALER**3
         !    SCALEOD 	= (SCALEQ * SCALEM) / (SCALEDENS * SCALER)
         !          	= SCALEQ * SCALER**2
         !
         ! Cap aerosol values at 90% relative humidity since
         ! aerosol growth at that point becomes highly nonlinear and 
         ! relative humidities above this value essentially mean
         ! there is a cloud in that grid box
         !
         ! Q is the extinction efficiency
         !
         ! Each grid box (I,J,L) will fall into one of the RH bins, 
         ! since each grid box will have a different RH value.  So,
         ! for SCALEOD(I,J,L,:), only one of the IRH bins will contain
         ! nonzero data, while the other IRH bins will all be zero.
         !==============================================================
 
         ! Loop over relative humidity bins
         DO R = 1, NRH

            ! Wet radius in "jv_spec.dat"
            RW(R) = RAA(4,IND(N)+R-1)	

            ! Wet frac of aerosol 
            FWET  = (RW(R)**3 - RW(1)**3) / RW(R)**3 

            ! Extinction efficiency Q for each RH bin
            QW(R) = QAA(4,IND(N)+R-1)*FWET + QAA(4,IND(N))*(1.d0-FWET)
         ENDDO

         ! Loop over SMVGEAR grid boxes
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, IRH, JLOOP, SCALEQ, SCALER, REFF, FRAC )
!$OMP+SCHEDULE( DYNAMIC )
         DO JLOOP = 1, NTLOOP

            ! Get 3-D grid box indices
            I = IXSAVE(JLOOP)
            J = IYSAVE(JLOOP)
            L = IZSAVE(JLOOP)

            ! Sort into relative humidity bins
            IF (      ABSHUM(JLOOP) <= RH(2) ) THEN   
               IRH = 1
            ELSE IF ( ABSHUM(JLOOP) <= RH(3) ) THEN
               IRH = 2
            ELSE IF ( ABSHUM(JLOOP) <= RH(4) ) THEN
               IRH = 3
            ELSE IF ( ABSHUM(JLOOP) <= RH(5) ) THEN
               IRH = 4
            ELSE 
               IRH = 5
            ENDIF

            ! For the NRHth bin, we don't have to interpolate
            ! For the other bins, we have to interpolate 
            IF ( IRH == NRH ) THEN
               SCALEQ = QW(NRH) / QW(1)  !QW(1) is dry extinction eff.
               REFF   = RW(NRH) 

            ELSE                

               ! Interpolate between different RH
               FRAC = (ABSHUM(JLOOP)-RH(IRH)) / (RH(IRH+1)-RH(IRH))
               IF ( FRAC > 1.0d0 ) FRAC = 1.0d0
               
               SCALEQ = (FRAC*QW(IRH+1) + (1.d0-FRAC)*QW(IRH)) / QW(1)
               REFF   = FRAC*RW(IRH+1)  + (1.d0-FRAC)*RW(IRH)

            ENDIF

            SCALER                 = REFF / RW(1)
            SCALEOD(I,J,L,IRH)     = SCALEQ * SCALER * SCALER
            SCALEVOL(I,J,L)        = SCALER**3
            ERADIUS(JLOOP,NDUST+N) = 1.0D-4 * REFF

            !==============================================================
            ! ND21 Diagnostic: 
            !
            ! Computed here:
            ! --------------
            ! #7  Hygroscopic growth of SO4                [unitless]
            ! #10 Hygroscopic growth of Black Carbon       [unitless]
            ! #13 Hygroscopic growth of Organic Carbon     [unitless]
            ! #16 Hygroscopic growth of Sea Salt (accum)   [unitless]
            ! #19 Hygroscopic growth of Sea Salt (coarse)  [unitless]
            !==============================================================
            IF ( ND21 > 0 .and. L <= LD21 ) THEN
               AD21(I,J,L,4+3*N) = AD21(I,J,L,4+3*N) +SCALEOD(I,J,L,IRH)
            ENDIF

         ENDDO
!$OMP END PARALLEL DO
      
         !==============================================================
         ! Convert concentration [kg/m3] to optical depth [unitless].
         !
         ! ODAER = ( 0.75 * BXHEIGHT * AERSL * QAA ) / 
         !         ( MSDENS * RAA * 1e-6 )
         ! (see Tegen and Lacis, JGR, 1996, 19237-19244, eq. 1)
         !
         ! Units ==> AERSL    [ kg/m3    ]
         !           MSDENS   [ kg/m3    ]
         !           RAA      [ um       ]  
         !           BXHEIGHT [ m        ]
         !           QAA      [ unitless ]
         !           ODAER    [ unitless ]
         !
         ! NOTES: 
         ! (1 ) Do the calculation at QAA(4,:) (i.e. 999 nm).          
         ! (2 ) RAA is the 'effective radius', Hansen and Travis, 1974
         ! (3 ) Report at the more relevant QAA(2,:) (i.e. 400 nm)   
         !       Although SCALEOD would be slightly different at 400nm 
         !       than at 1000nm as done here, FAST-J does currently 
         !       allow one to provide different input optical depths at 
         !       different wavelengths.  Therefore the reported value at
         !       determined with QAA(2,:) is as used in FAST-J. 
         ! (4 ) Now use explicit indices in parallel DO-loops, since
         !       some compilers may not like array masks in parallel
         !       regions (bmy, 2/28/02)
         !==============================================================
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, R, IRHN ) 
!$OMP+SCHEDULE( DYNAMIC )
         DO R = 1, NRH
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Bin for aerosol type and relative humidity
            IRHN = ( (N-1) * NRH ) + R

            ! Save aerosol optical depth for each combination 
            ! of aerosol type and relative humidity into ODAER, 
            ! which will get passed to the FAST-J routines
            ODAER(I,J,L,IRHN) = SCALEOD(I,J,L,R) 
     &                        * 0.75d0 * BXHEIGHT(I,J,L) 
     &                        * WAERSL(I,J,L,N) * QAA(4,IND(N)) / 
     &                        ( MSDENS(N) * RAA(4,IND(N)) * 1.0D-6 )

         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         !==============================================================
         !  Calculate Aerosol Surface Area
         !
         !  Units ==> AERSL    [ kg aerosol m^-3 air ]
         !            MSDENS   [ kg aerosol m^-3 aerosol ]
         !            ERADIUS  [ cm      ]
         !            TAREA    [ cm^2 dry aerosol/cm^3 air ]
         !
         !  Note: first find volume of aerosol (cm^3 arsl/cm^3 air), then
         !        multiply by 3/radius to convert to surface area in cm^2
         !
         !  Wet Volume = AERSL * SCALER**3 / MSDENS
         !  Wet Surface Area = 3 * (Wet Volume) / ERADIUS 
         !
         !  Effective radius for surface area and optical depths 
         !  are identical.
         !==============================================================
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, JLOOP ) 
!$OMP+SCHEDULE( DYNAMIC )
         DO JLOOP = 1, NTLOOP

            ! Get 3-D grid box indices
            I = IXSAVE(JLOOP)
            J = IYSAVE(JLOOP)
            L = IZSAVE(JLOOP)
          
            ! Store aerosol surface areas in TAREA, and be sure
            ! to list them following the dust surface areas
            TAREA(JLOOP,N+NDUST) = 3.D0                    * 
     &                             WAERSL(I,J,L,N)         *  
     &                             SCALEVOL(I,J,L)         / 
     &                            ( ERADIUS(JLOOP,NDUST+N) * 
     &                              MSDENS(N) )  
         ENDDO
!$OMP END PARALLEL DO

      ENDDO  !Loop over NAER

      !### Debug
      !CALL DEBUG_MSG( '### RDAER: after loop over NAER' )

      !==============================================================
      ! Account for hydrophobic aerosols (BC and OC), N=2 and N=3
      !==============================================================
      DO N = 2, 3

         ! Index for combination of aerosol type and RH
         IRHN = ( (N-1) * NRH ) + 1

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L ) 
!$OMP+SCHEDULE( DYNAMIC )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Aerosol optical depth
            ODAER(I,J,L,IRHN) = ODAER(I,J,L,IRHN) + 
     &                          0.75d0            * BXHEIGHT(I,J,L) * 
     &                          DAERSL(I,J,L,N-1) * QAA(4,IND(N))   / 
     &                          ( MSDENS(N) * RAA(4,IND(N)) * 1.0D-6 )
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! Effective radius
         REFF = 1.0D-4 * RAA(4,IND(N))

         ! Loop over grid boxes
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, JLOOP, DRYAREA ) 
!$OMP+SCHEDULE( DYNAMIC )
         DO JLOOP = 1, NTLOOP

            ! Get 3-D grid box indices
            I = IXSAVE(JLOOP)
            J = IYSAVE(JLOOP)
            L = IZSAVE(JLOOP)

            ! Dry surface area
            DRYAREA = 3.D0 * DAERSL(I,J,L,N-1) / ( REFF * MSDENS(N) )  

            ! Add surface area to TAREA array
            TAREA(JLOOP,N+NDUST) = TAREA(JLOOP,N+NDUST) + DRYAREA

            ! Define a new effective radius that accounts 
            ! for the hydrophobic aerosol 
            ERADIUS(JLOOP,NDUST+N) = ( ERADIUS(JLOOP,NDUST+N) * 
     &                                  TAREA(JLOOP,N+NDUST)  +
     &                                  REFF * DRYAREA)       / 
     &                               ( TAREA(JLOOP,N+NDUST) + DRYAREA )

         ENDDO
!$OMP END PARALLEL DO

      ENDDO
       
      !==============================================================
      ! ND21 Diagnostic: Aerosol OD's, Growth Rates, Surface Areas
      !
      ! Computed in other routines:
      ! ---------------------------------
      ! #1: Cloud optical depths (1000 nm) --> from "optdepth_mod.f"       
      ! #2: Max Overlap Cld Frac           --> from "optdepth_mod.f" 
      ! #3: Random Overlap Cld Frac        --> from "optdepth_mod.f" 
      ! #4: Dust optical depths (400 nm)   --> from "rdust.f"
      ! #5: Dust surface areas             --> from "rdust.f"
      !
      ! Computed previously in "rdaer.f":
      ! ---------------------------------
      ! #7  Hygroscopic growth of SO4                [unitless]
      ! #10 Hygroscopic growth of Black Carbon       [unitless]
      ! #13 Hygroscopic growth of Organic Carbon     [unitless]
      ! #16 Hygroscopic growth of Sea Salt (accum)   [unitless]
      ! #19 Hygroscopic growth of Sea Salt (coarse)  [unitless]
      !
      ! Computed here:
      ! ---------------------------------
      ! #6  Sulfate Optical Depth (400 nm)           [unitless]
      ! #8  Sulfate Surface Area                     [cm2/cm3 ]
      ! #9  Black Carbon Optical Depth (400 nm)      [unitless]
      ! #11 Black Carbon Surface Area                [cm2/cm3 ]
      ! #12 Organic Carbon Optical Depth (400 nm)    [unitless]
      ! #14 Organic Carbon Surface Area              [cm2/cm3 ]
      ! #15 Sea Salt (accum) Opt Depth (400 nm)      [unitless]
      ! #17 Sea Salt (accum) Surface Area            [cm2/cm3 ]
      ! #18 Sea Salt (coarse) Opt Depth(400 nm)      [unitless]
      ! #20 Sea Salt (coarse) Surface Area           [cm2/cm3 ]
      !
      ! NOTE: The cloud optical depths are actually recorded at
      !       1000 nm, but vary little with wavelength.
      !==============================================================
      IF ( ND21 > 0 ) THEN

         ! Loop over aerosol types
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, IRHN, J, JLOOP, L, N, R ) 
!$OMP+SCHEDULE( DYNAMIC )
         DO N = 1, NAER

            !------------------------------------
            ! Optical Depths 
            ! Scale of optical depths w/ RH 
            !------------------------------------
            DO R = 1, NRH
            DO L = 1, LD21
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! Index for type of aerosol and RH value
               IRHN = ( (N-1) * NRH ) + R

               ! Optical Depths
               AD21(I,J,L,3+3*N) = AD21(I,J,L,3+3*N) + 
     &                             ODAER(I,J,L,IRHN) * 
     &                             QAA(2,IND(N)) /QAA(4,IND(N))

            ENDDO
            ENDDO
            ENDDO
            ENDDO

            !------------------------------------
            ! Surface areas
            !------------------------------------
            DO JLOOP = 1, NTLOOP

               ! Get 3-D grid box indices
               I = IXSAVE(JLOOP)
               J = IYSAVE(JLOOP)
               L = IZSAVE(JLOOP)

               ! Add aerosol surface areas 
               IF ( L <= LD21 ) THEN
                  AD21(I,J,L,5+3*N) = AD21(I,J,L,5+3*N) + 
     &                                TAREA(JLOOP,N+NDUST)
               ENDIF 
            ENDDO

         ENDDO
!$OMP END PARALLEL DO

      ENDIF 

      !=================================================================
      ! To turn off the radiative effects of different aerososl
      ! uncomment the following lines
      !=================================================================
      !DO R = 1,NRH
      !  ODAER(:,:,:,R)       = 0.d0  !sulfate
      !  ODAER(:,:,:,R+NRH)   = 0.d0  !BC
      !  ODAER(:,:,:,R+2*NRH) = 0.d0  !OC
      !  ODAER(:,:,:,R+3*NRH) = 0.d0  !SS(accum)
      !  ODAER(:,:,:,R+4*NRH) = 0.d0  !SS(coarse)
      !ENDDO

      !=================================================================
      ! To turn off heterogeneous chemistry on different aerosols
      ! uncomment the following lines
      !=================================================================
      !TAREA(:,NDUST+1) = 0.d0	!Sulfate
      !TAREA(:,NDUST+2) = 0.d0	!BC 
      !TAREA(:,NDUST+3) = 0.d0	!OC 
      !TAREA(:,NDUST+4) = 0.d0	!SS (accum)
      !TAREA(:,NDUST+5) = 0.d0	!SS (coarse)

      ! Return to calling program
      END SUBROUTINE RDAER
