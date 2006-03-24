! $Id: emfossil.f,v 1.9 2006/03/24 20:22:46 bmy Exp $
      SUBROUTINE EMFOSSIL( I, J, N, NN, IREF, JREF, JSCEN )
!
!******************************************************************************
!  Subroutine EMFOSSIL emits fossil fuels into the EMISRR and EMISRRN 
!  arrays, which are then passed to SMVGEAR. (bmy, 4/19/99, 2/1/06)
!
!  Arguments as input:
!  ========================================================================
!  (1-2) I, J       : longitude and latitude indices
!  (3-4) N, NN      : Emission index and tracer index
!  (5-6) IREF, JREF : Offset indices I+I0 and J+J0
!  (7  ) JSCEN      : Day index (Sat=1, Sun=2, Weekday=3)
!
!  NOTES:
!  (1 ) Uses the correct seasonal NOx and multi-level NOx (anthroems.f)
!  (2 ) Uses anthro scale factors for years since 1985 (from anthroems.f)
!  (3 ) Scales emissions based on weekday/weekend (emf_scale.f)
!  (4 ) Preserves old sensitivity study cases (emf_scale.f, emissdr.f)
!  (5 ) Scales emissions based on time of day (emfossil.f)
!  (6 ) Get rid of all GISS and PLUMES code (bmy, 4/19/99)
!  (7 ) Now use F90 syntax for declarations, etc. (bmy, 4/19/99)
!  (8 ) Now use allocatable arrays for ND29 and ND36 diagnostics.
!        Also made minor cosmetic changes & updated comments. (bmy, 3/16/00)
!  (9 ) Eliminate obsolete code and ND63 diagnostic (bmy, 4/12/00)
!  (10) Enhance anthropogenic CO emission by 8%, to account for CO production
!        from oxidation of anthropogenic VOC's (bnd, bmy, 1/2/01)
!  (11) Comment out scaling by 1.08 for anthro CO (bmy, 2/12/01)
!  (12) Eliminate obsolete commented-out code (bmy, 4/20/01)
!  (13) Now use 2% as the enhancment factor for CO instead of 1.08,
!        according to new jal numbers (bmy, 4/26/01)
!  (14) Now references "tracerid_mod.f" (bmy, 11/6/02)
!  (15) Now replaced DXYP(JREF)*1d4 with GET_AREA_CM2(J).  Now use function
!        GET_TS_EMIS() from "time_mod.f" (bmy, 2/11/03)
!  (16) Now can overwrite existing emissions with EPA/NEI data over the 
!        continental USA if LNEI99=T.  Now reference LNEI99 from F90
!        module "logical_mod.f".  Now reference GET_EPA_ANTHRO and
!        GET_USA_MASK from "epa_nei_mod.f". (rch, rjp, bmy, 11/5/04)
!  (17) Now references GET_DAY_OF_WEEK from "time_mod.f" to correctly figure
!        out if this is a weekday or weekend. (bmy, 7/6/05)
!  (18) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (19) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (20) Now apply EMEP European emissions if necessary.  Remove reference
!        to CMN, it's now obsolete. (bdf, bmy, 11/1/05)
!  (21) Rewrite IF statements to avoid seg fault errors when LEMEP and LNEI99 
!        are turned off. (bmy, 2/1/06)
!******************************************************************************
!          
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD29, AD32_an,   AD36
      USE EMEP_MOD,     ONLY : GET_EMEP_ANTHRO, GET_EUROPE_MASK
      USE EPA_NEI_MOD,  ONLY : GET_EPA_ANTHRO,  GET_USA_MASK
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,  ONLY : LNEI99,          LEMEP
      USE TIME_MOD,     ONLY : GET_TS_EMIS,     GET_DAY_OF_WEEK
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDENOX, IDEOX,   IDTCO
  
      IMPLICIT NONE

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! Diagnostic switches & arrays
#     include "CMN_O3"       ! EMISR, EMISRR, etc...
#     include "comode.h"     ! IHOUR

      ! Arguments
      INTEGER, INTENT(IN)   :: I, J, N, NN, IREF, JREF, JSCEN

      ! Local Variables & external functions
      LOGICAL               :: WEEKDAY
      INTEGER               :: L, LL, K, DOW
      REAL*8                :: TODX, DTSRCE, AREA_CM2           
      REAL*8                :: EMX(NOXLEVELS)  
      REAL*8                :: XEMISR
      REAL*8                :: XEMISRN(NOXLEVELS)
      REAL*8                :: EPA_NEI, EMEP

      !=================================================================
      ! EMFOSSIL begins here!
      !=================================================================

      ! Emission timestep [s]
      DTSRCE   = GET_TS_EMIS() * 60d0

      ! Surface area of grid box
      AREA_CM2 = GET_AREA_CM2( J )

      ! Flag for weekday or weekend for NEI emissions
      DOW      = GET_DAY_OF_WEEK()
      WEEKDAY  = ( DOW > 0 .and. DOW < 6 )

      !=================================================================
      ! Call EMF_SCALE to do the following:
      ! (1) Save original values of EMISR, EMISRN
      ! (2) If LFFNOX=F, turn off NOx, Ox emissions
      ! (3) Scale emissions to weekend/weekday usage
      !=================================================================
      CALL EMF_SCALE( I, J, N, NN, IREF, JREF, JSCEN, XEMISR, XEMISRN )       

      !=================================================================    
      ! ADD ANTHROPOGENIC EMISSIONS TO TRACER TOTALS
      ! NOTE APPROPRIATE TIME-OF-DAY FACTOR (TOD) MUST BE
      ! ESTABLISHED FOR EACH TRACER; 
      ! WITH IHOUR = 1-6 (1 = 10pm-2am)
      ! and tracer index distinguishing NOx-HC- BIO
      !
      ! NOx only: account for all NOx levels (LL=1,NOXLEVELS)
      !=================================================================    
      IF ( N == IDENOX ) THEN

         ! Initialize work variables
         EMX(:)  = 0d0

         ! Time of day factor for NOx
         TODX = TODN(IHOUR)
         
         ! Loop over all of the emission levels for NOx (e.g. surface, 100m)
         DO LL = 1, NOXLEVELS
            EMX(LL)  = TODX * EMISRN(IREF,JREF,LL)

            !-----------------------------------------------------------
            ! Get NOx from EPA/NEI inventory over the USA 
            !-----------------------------------------------------------

            ! If we are using EPA/NEI emissions
            IF ( LNEI99 ) THEN

               ! If we are over the USA ...               
               IF ( GET_USA_MASK( I, J ) > 0d0 ) THEN 
            
                  ! Get EPA emissions for NOx (and apply time-of-day factor)
                  EPA_NEI = GET_EPA_ANTHRO( I, J, NN, WEEKDAY )
                  EPA_NEI = EPA_NEI * TODX
            
                  IF ( LL == 1 ) THEN
                  
                     ! Replace GEIA with EPA/NEI emissions at surface
                     EMX(LL) = EPA_NEI * 
     &                         ( DTSRCE * AREA_CM2 ) / XNUMOL(NN) 
            
                  ELSE 
            
                     ! Zero GEIA emissions in the 2nd level 
                     ! where the EPA/NEI emissions are nonzero
                     EMX(LL) = 0d0                   
                  
                  ENDIF
               ENDIF
            ENDIF

            !-----------------------------------------------------------
            ! Get NOx from EMEP inventory over Europe 
            !-----------------------------------------------------------

            ! If we are using EMEP ...
            IF ( LEMEP ) THEN

               ! If we are over the European region ...
               IF ( GET_EUROPE_MASK( I, J ) > 0d0 ) THEN 

                  ! Get EMEP emissions for NOx (and apply time-of-day factor)
                  EMEP = GET_EMEP_ANTHRO( I, J, NN )
                  EMEP = EMEP * TODX

                  IF ( LL == 1 ) THEN
            
                     ! Replace GEIA with EMEP emissions at surface
                     EMX(LL) = EMEP * ( DTSRCE * AREA_CM2 ) / XNUMOL(NN) 
                  
                  ELSE 
            
                     ! Zero GEIA emissions in the 2nd level 
                     ! where the EMEP emissions are nonzero
                     EMX(LL) = 0d0                   
                  
                  ENDIF
               ENDIF
            ENDIF

            !-----------------------------------------------------------
            ! Store in EMISRRN array and archive diagnostics
            !-----------------------------------------------------------

            ! EMISRRN [molec/box/s] is referenced by LL
            EMISRRN(I,J,LL) = EMISRRN(I,J,LL) +
     &                       ( EMX(LL) * XNUMOL(NN) / DTSRCE )

            ! ND32 = save anthro NOx for levels L=1,NOXEXTENT [molec/cm2/s]
            IF ( ND32 > 0 ) THEN
               AD32_an(I,J,LL) = AD32_an(I,J,LL) + 
     &              ( EMX(LL) * XNUMOL(NN) / ( DTSRCE * AREA_CM2 ) )   
            ENDIF

            ! ND36 = save anthro emissions in [molec/cm2]
            ! and then convert to [molec/cm2/s] in DIAG3.F
            IF ( ND36 > 0 ) THEN
               AD36(I,J,N) = AD36(I,J,N) +
     &            ( EMX(LL) * XNUMOL(NN) / AREA_CM2 )
            ENDIF    
         ENDDO     

      !=================================================================    
      ! All other emitted tracers except NOx! 
      !=================================================================    
      ELSE

         ! Initialize work variables
         EMX(:)  = 0d0

         ! Use appropriate scale factor for time of day
         IF ( N == IDEOX ) THEN
            TODX = TODN(IHOUR)
         ELSE
            TODX = TODH(IHOUR)
         ENDIF

         EMX(1)        = TODX * EMISR(IREF,JREF,N) 

         ! Enhance CO production by 2%, to account for CO
         ! production from anthropogenic VOC's (bnd, bmy, 4/26/01)
         IF ( NN == IDTCO ) EMX(1) = EMX(1) * 1.02d0

         !--------------------------------------------------------------
         ! Get CO & Hydrocarbons from EPA/NEI inventory over the USA 
         !--------------------------------------------------------------

         ! If we are using EPA/NEI99 emissions ...
         IF ( LNEI99 ) THEN

            ! If we are over the USA ...
            IF ( GET_USA_MASK( I, J ) > 0d0 ) THEN
         
               ! Get EPA/NEI emissions (and apply time-of-day factor)
               EPA_NEI = GET_EPA_ANTHRO( I, J, NN, WEEKDAY )
               EPA_NEI = EPA_NEI * TODX
         
               ! Convert from molec/cm2/s to kg/box/timestep in order
               ! to be in the proper units for EMISRR array
               EMX(1)  = EPA_NEI * ( DTSRCE * AREA_CM2 ) / XNUMOL(NN) 

            ENDIF
         ENDIF

         !--------------------------------------------------------------
         ! Get CO & Hydrocarbons from EMEP inventory over Europe
         !--------------------------------------------------------------

         ! If we are using EMEP emissions ...
         IF ( LEMEP ) THEN

            ! If we are over the European region ...
            IF ( GET_EUROPE_MASK( I, J ) > 0d0 ) THEN
         
               ! Get EMEP emissions 
               EMEP = GET_EMEP_ANTHRO( I, J, NN )
         
               ! -1 indicates tracer NN does not have EMEP emissions
               IF ( EMEP > 0d0 ) THEN

                  ! Apply time-of-day factor
                  EMEP   = EMEP * TODX

                  ! Convert from molec/cm2/s to kg/box/timestep in order
                  ! to be in the proper units for EMISRR array
                  EMX(1) = EMEP * ( DTSRCE * AREA_CM2 ) / XNUMOL(NN) 

               ENDIF
            ENDIF
         ENDIF

         !--------------------------------------------------------------
         ! Store in EMISRR array and archive diagnostics
         !--------------------------------------------------------------
         EMISRR(I,J,N) = EMX(1) * XNUMOL(NN) / DTSRCE

         ! ND29 = CO source diagnostic... 
         ! store as [molec/cm2/s] in AD29(:,:,1)
         IF ( ND29 > 0 .and. NN == IDTCO ) THEN
            AD29(I,J,1) = AD29(I,J,1) +
     &           ( EMX(1) * XNUMOL(NN) / ( DTSRCE * AREA_CM2 ) ) 
         ENDIF

         ! ND36 = Anthro source diagnostic...store as [molec/cm2] 
         ! and convert to [molec/cm2/s] in DIAG3.F
         IF ( ND36 > 0 ) THEN
            AD36(I,J,N) = AD36(I,J,N) + 
     &           ( EMX(1) * XNUMOL(NN) / AREA_CM2 ) 
         ENDIF    
      ENDIF

      !=================================================================
      ! Restore EMISR, EMISRN to original values
      !=================================================================
      IF ( N == IDENOX ) THEN
         EMISRN(IREF,JREF,1:NOXLEVELS) = XEMISRN(1:NOXLEVELS)
      ELSE
         EMISR(IREF,JREF,N) = XEMISR
      ENDIF  

      ! Return to calling program
      END SUBROUTINE EMFOSSIL
