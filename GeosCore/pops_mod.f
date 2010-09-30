!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: pops_mod.f
!
! !DESCRIPTION: This module contains variables and routines for the 
!  GEOS-Chem POPs simulation. 
!\\
!\\
! !INTERFACE: 
!
      MODULE POPS_MOD
! 
! !USES:
!
      IMPLICIT NONE
! Make everything Private ...
      PRIVATE
!
! !PUBLIC TYPES:
!

! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CHEMPOPS
!      PUBLIC :: INIT_POPS

!
! !PUBLIC DATA MEMBERS:
!

! !REVISION HISTORY:
!  20 September 2010 N.E. Selin - Initial Version
!
! !REMARKS:
! Under construction
!
!EOP

      ! Parameters
      REAL*8,  PARAMETER   :: SMALLNUM = 1D-20
      ! Arrays
      REAL*8,  ALLOCATABLE :: TCOSZ(:,:)
      REAL*8,  ALLOCATABLE :: TTDAY(:,:)
      REAL*8,  ALLOCATABLE :: ZERO_DVEL(:,:)
      REAL*8,  ALLOCATABLE :: COSZM(:,:)
     
      CONTAINS
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CHEMPOPS
!
! !DESCRIPTION: This routine is the driver routine for POPs chemistry 
!  (eck, 9/20/10)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEMPOPS
!
! !INPUT PARAMETERS: 
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
 
      ! References to F90 modules
!      USE DIAG_MOD,      ONLY : AD44
      USE DRYDEP_MOD,    ONLY : DEPSAV
      USE ERROR_MOD,     ONLY : DEBUG_MSG
      USE GLOBAL_OH_MOD, ONLY : GET_GLOBAL_OH
      USE PBL_MIX_MOD,   ONLY : GET_PBL_MAX_L
      USE LOGICAL_MOD,   ONLY : LPRT, LGTMM, LNLPBL !CDH added LNLPBL
      USE TIME_MOD,      ONLY : GET_MONTH, ITS_A_NEW_MONTH
      USE TRACER_MOD,    ONLY : N_TRACERS
      USE TRACERID_MOD,  ONLY : N_HG_CATS
      USE DRYDEP_MOD,    ONLY : DRYPOPG, DRYPOPP

#     include "CMN_SIZE"      ! Size parameters

!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version
!
! !REMARKS:
! (1) Based initially on CHEMMERCURY from MERCURY_MOD (eck, 9/20/10)
!
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Local variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: I, J, L, MONTH, N, PBL_MAX

      REAL*8                 :: K_DRYD2(IIPAR,JJPAR)

      !=================================================================
      ! CHEMPOPS begins here!
      !
      ! Read monthly mean OH fields
      !=================================================================
      IF ( ITS_A_NEW_MONTH() ) THEN 

         ! Get the current month
         MONTH = GET_MONTH()

         ! Read monthly mean OH and O3 from disk
         CALL GET_GLOBAL_OH( MONTH )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMPOPS: a GET_GLOBAL_OH' )

      ENDIF
      
      !=================================================================
      ! Perform chemistry on POPS TRACERS
      !=================================================================
      
      ! Compute diurnal scaling for OH
      CALL OHNO3TIME
      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMPOPS: a OHNO3TIME' )

      !-------------------------
      ! GAS AND PARTICLE PHASE chemistry
      !-------------------------
      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMPOPS: b CHEM_GASPART' )
      
      ! Add option for non-local PBL (cdh, 08/27/09)
      IF ( LNLPBL ) THEN

         ! Dry deposition occurs with PBL mixing,
         ! pass zero deposition frequency
         CALL CHEM_POPGP( ZERO_DVEL, ZERO_DVEL)
         
      ELSE

         IF ( DRYPOPG > 0 .and. DRYPOPP > 0 ) THEN
         
            ! Dry deposition active for both POP-Gas and POP-Particle; 
            ! pass drydep frequency to CHEM_POPGP (NOTE: DEPSAV has units 1/s)
            CALL CHEM_POPGP( K_DRYD2, DEPSAV(:,:,DRYPOPP) )
            
         ELSEIF (DRYPOPG > 0 .and. DRYPOPP .le. 0) THEN

            ! Only POPG dry deposition is active
            CALL CHEM_POPGP( K_DRYD2, ZERO_DVEL) 
            
         ELSEIF (DRYPOPG <= 0 .and. DRYPOPP > 0) THEN

            ! Only POPP dry deposition is active
            CALL CHEM_POPGP( ZERO_DVEL , DEPSAV(:,:,DRYPOPP))
            
         ELSE

            ! No dry deposition, pass zero deposition frequency
            CALL CHEM_POPGP( ZERO_DVEL , ZERO_DVEL)

         ENDIF

      ENDIF      

      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMPOPS: a CHEM_GASPART' )
   
    
      ! Return to calling program
      END SUBROUTINE CHEMPOPS

!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CHEM_POPGP
!
! !DESCRIPTION: This routine does chemistry for POPs gas and particles
!  (eck, 9/20/10)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_POPGP (V_DEP_G, V_DEP_P)

      USE TRACER_MOD,   ONLY : STT,    XNUMOL
      USE TRACERID_MOD, ONLY: IDTPOPG, IDTPOPP
      USE DIAG53_MOD,   ONLY: AD53_PG_PP
      USE TIME_MOD,     ONLY: GET_TS_CHEM
      USE DIAG_MOD,     ONLY : AD44
      USE PBL_MIX_MOD,  ONLY: GET_FRAC_UNDER_PBLTOP
#     include "CMN_SIZE" !Size parameters
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN)    :: V_DEP_G(IIPAR,JJPAR)
      REAL*8, INTENT(IN)    :: V_DEP_P(IIPAR,JJPAR)

!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version
!
! !REMARKS:
! (1) Based initially on CHEM_HG0_HG2 from MERCURY_MOD (eck, 9/20/10)
!
!EOP
!------------------------------------------------------------------------------
!BOC


!Here, put some code that reacts popg into popp and deposits both
! save out dry dep tracers
! save into stt

      END SUBROUTINE CHEM_POPGP 

!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  EMISSPOPS
!
! !DESCRIPTION: This routine is the driver routine for POPs emissions
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISSPOPS
!
! !INPUT PARAMETERS: 
!

!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version
!
! !REMARKS:
! (1) Based initially on EMISSMERCURY from MERCURY_MOD (eck, 9/20/10)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      ! References to F90 modules
      USE ERROR_MOD,         ONLY : DEBUG_MSG, ERROR_STOP
      USE LOGICAL_MOD,       ONLY : LPRT, LNLPBL !CDH added LNLPBL
      USE TIME_MOD,          ONLY : GET_MONTH, ITS_A_NEW_MONTH
      USE TRACER_MOD,        ONLY : STT
      USE VDIFF_PRE_MOD,     ONLY : EMIS_SAVE !cdh for LNLPBL
      USE GRID_MOD,          ONLY : GET_XMID, GET_YMID
      ! Reference to diagnostic arrays
      USE DIAG53_MOD,   ONLY : AD53, ND53, AD53_PG_PP
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_OF_PBL, GET_PBL_MAX_L
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE TRACERID_MOD, ONLY : IDTPOPG, IDTPOPP
      
#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DEP"      ! FRCLND

      ! Local variables
      INTEGER               :: I,      J,    L,        N,    PBL_MAX
      REAL*8                :: DTSRCE, E_POP, F_OF_PBL, T_POP
      LOGICAL, SAVE          :: FIRST = .TRUE.
      !=================================================================
      ! EMISSPOPS begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Read anthro emissions from disk
         !CALL POPS_READYR
         !To be written when we read from disk

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! Emission timestep [s]
      DTSRCE  = GET_TS_EMIS() * 60d0

      ! Maximum extent of the PBL [model levels]
      PBL_MAX = GET_PBL_MAX_L() 

     
      ! If we are using the non-local PBL mixing,
      ! we need to initialize the EMIS_SAVE array (cdh, 08/27/09)
      IF (LNLPBL) EMIS_SAVE = 0d0

      ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N,  T_POP, F_OF_PBL, E_POP)
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !Here, save the total from the anthro emissions array
         !into the T_POP variable (for now emit one unit)
         T_POP=1d0
 
         !==============================================================
         ! Partition POP throughout PBL; store into STT [kg]
         ! Now make sure STT does not underflow (cdh, bmy, 4/6/06; eck 9/20/10)
         !==============================================================

         ! Loop up to max PBL level
         DO L = 1, PBL_MAX

            ! Fraction of box (I,J,L) w/in the PBL [unitless]
            F_OF_PBL        = GET_FRAC_OF_PBL( I, J, L )

            !-----------------
            ! POP GAS (no particulate emission)
            !-----------------
            N               = IDTPOPG
            E_POP            = F_OF_PBL * T_POP * DTSRCE
            CALL EMITPOP( I, J, L, N, E_POP )
         ENDDO



      ENDDO
      ENDDO
!$OMP END PARALLEL DO

         !==============================================================
         ! ND53 diagnostic: Total POP emissions [kg]
         !==============================================================
         IF ( ND53 > 0 ) THEN
            AD53(I,J,1) = AD53(I,J,1) + ( T_POP       * DTSRCE )
     
         ENDIF


      ! Return to calling program
      END SUBROUTINE EMISSPOPS

!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  EMITPOP
!
! !DESCRIPTION: This routine directs emission either to STT directly or to EMIS_SAVE
!  for use by the non-local PBL mixing. This is a programming convenience.
!  (cdh, 08/27/09, modified for pops by eck, 9/20/10)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMITPOP( I, J, L, ID, E_POP )
!! 
! !USES:
      ! Reference to diagnostic arrays
      USE TRACER_MOD,   ONLY : STT
      USE LOGICAL_MOD,  ONLY : LNLPBL
      USE VDIFF_PRE_MOD,ONLY : EMIS_SAVE
! !INPUT PARAMETERS: 
      INTEGER, INTENT(IN)   :: I, J, L, ID
      REAL*8,  INTENT(IN)   :: E_POP
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version
!
! !REMARKS:
! (1) Based initially on EMITHG from MERCURY_MOD (eck, 9/20/10)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! EMITPOP begins here!
      !=================================================================

      ! Save emissions for non-local PBL mixing or emit directly.
      ! Make sure that emitted mass is non-negative
      ! This is here only for consistency with old code which warned of
      ! underflow error (cdh, 08/27/09, modified for POPs 9/20/10)
      IF (LNLPBL) THEN
         EMIS_SAVE(I,J,ID) = EMIS_SAVE(I,J,ID) + MAX( E_POP, 0D0 )
      ELSE
         STT(I,J,L,ID) = STT(I,J,L,ID) + MAX( E_POP, 0D0 )
      ENDIF

      END SUBROUTINE EMITPOP
!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  OHNO3TIME
!
! !DESCRIPTION: Subroutine OHNO3TIME computes the sum of cosine of the solar zenith
!  angle over a 24 hour day, as well as the total length of daylight. 
!  This is needed to scale the offline OH and NO3 concentrations.
!  (rjp, bmy, 12/16/02, 12/8/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OHNO3TIME
!! 
! !USES:
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XMID,    GET_YMID_R
      USE TIME_MOD, ONLY : GET_NHMSb,   GET_ELAPSED_SEC
      USE TIME_MOD, ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT


#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_GCTM"  ! Physical constants
! !INPUT PARAMETERS: 
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version for POPS_MOD
!
! !REMARKS:
!  (1 ) Copy code from COSSZA directly for now, so that we don't get NaN
!        values.  Figure this out later (rjp, bmy, 1/10/03)
!  (2 ) Now replace XMID(I) with routine GET_XMID from "grid_mod.f".  
!        Now replace RLAT(J) with routine GET_YMID_R from "grid_mod.f". 
!        Removed NTIME, NHMSb from the arg list.  Now use GET_NHMSb,
!        GET_ELAPSED_SEC, GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT from 
!        "time_mod.f". (bmy, 3/27/03)
!  (3 ) Now store the peak SUNCOS value for each surface grid box (I,J) in 
!        the COSZM array. (rjp, bmy, 3/30/04)
!  (4 ) Also added parallel loop over grid boxes (eck, bmy, 12/8/04)
!  (5 ) copied from mercury_mod by eck (9/20/10)
!******************************************************************************
!
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Local variables
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER             :: I, IJLOOP, J, L, N, NT, NDYSTEP
      REAL*8              :: A0, A1, A2, A3, B1, B2, B3
      REAL*8              :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
      REAL*8              :: SUNTMP(MAXIJ)
      
      !=================================================================
      ! OHNO3TIME begins here!
      !=================================================================

      !  Solar declination angle (low precision formula, good enough for us):
      A0 = 0.006918
      A1 = 0.399912
      A2 = 0.006758
      A3 = 0.002697
      B1 = 0.070257
      B2 = 0.000907
      B3 = 0.000148
      R  = 2.* PI * float( GET_DAY_OF_YEAR() - 1 ) / 365.

      DEC = A0 - A1*cos(  R) + B1*sin(  R)
     &         - A2*cos(2*R) + B2*sin(2*R)
     &         - A3*cos(3*R) + B3*sin(3*R)

      LHR0 = int(float( GET_NHMSb() )/10000.)

      ! Only do the following at the start of a new day
      IF ( FIRST .or. GET_GMT() < 1e-5 ) THEN 
      
         ! Zero arrays
         TTDAY(:,:) = 0d0
         TCOSZ(:,:) = 0d0
         COSZM(:,:) = 0d0

         ! NDYSTEP is # of chemistry time steps in this day
         NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 60 / GET_TS_CHEM()         

         ! NT is the elapsed time [s] since the beginning of the run
         NT = GET_ELAPSED_SEC()

         ! Loop forward through NDYSTEP "fake" timesteps for this day 
         DO N = 1, NDYSTEP
            
            ! Zero SUNTMP array
            SUNTMP(:) = 0d0

            ! Loop over surface grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, YMID_R, IJLOOP, TIMLOC, AHR )
            DO J = 1, JJPAR

               ! Grid box latitude center [radians]
               YMID_R = GET_YMID_R( J )

            DO I = 1, IIPAR

               ! Increment IJLOOP
               IJLOOP = ( (J-1) * IIPAR ) + I
               TIMLOC = real(LHR0) + real(NT)/3600.0 + GET_XMID(I)/15.0
         
               DO WHILE (TIMLOC .lt. 0)
                  TIMLOC = TIMLOC + 24.0
               ENDDO

               DO WHILE (TIMLOC .gt. 24.0)
                  TIMLOC = TIMLOC - 24.0
               ENDDO

               AHR = abs(TIMLOC - 12.) * 15.0 * PI_180

            !===========================================================
            ! The cosine of the solar zenith angle (SZA) is given by:
            !     
            !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR) 
            !                   
            ! where LAT = the latitude angle, 
            !       DEC = the solar declination angle,  
            !       AHR = the hour angle, all in radians. 
            !
            ! If SUNCOS < 0, then the sun is below the horizon, and 
            ! therefore does not contribute to any solar heating.  
            !===========================================================

               ! Compute Cos(SZA)
               SUNTMP(IJLOOP) = sin(YMID_R) * sin(DEC) +
     &                          cos(YMID_R) * cos(DEC) * cos(AHR)

               ! TCOSZ is the sum of SUNTMP at location (I,J)
               ! Do not include negative values of SUNTMP
               TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP(IJLOOP), 0d0 )

               ! COSZM is the peak value of SUMTMP during a day at (I,J)
               ! (rjp, bmy, 3/30/04)
               COSZM(I,J) = MAX( COSZM(I,J), SUNTMP(IJLOOP) )

               ! TTDAY is the total daylight time at location (I,J)
               IF ( SUNTMP(IJLOOP) > 0d0 ) THEN
                  TTDAY(I,J) = TTDAY(I,J) + DBLE( GET_TS_CHEM() )
               ENDIF
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            ! Increment elapsed time [sec]
            NT = NT + ( GET_TS_CHEM() * 60 )             
         ENDDO

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! Return to calling program
      END SUBROUTINE OHNO3TIME

!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  INIT_POPS
!
! !DESCRIPTION: Subroutine INIT_POPS allocates and zeroes all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_POPS
!
      ! References to F90 modules
      USE DRYDEP_MOD,   ONLY : DEPNAME,   NUMDEP
      USE ERROR_MOD,    ONLY : ALLOC_ERR, ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LSPLIT,    LDRYD,     LNLPBL
      USE TRACER_MOD,   ONLY : N_TRACERS
c$$$
#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! ND44

! !INPUT PARAMETERS: 
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version
!
! !REMARKS:
! (1) Based initially on INIT_MERCURY from MERCURY_MOD (eck, 9/20/10)
!
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Local variables
      LOGICAL, SAVE         :: IS_INIT = .FALSE. 
      INTEGER               :: AS, N
      !=================================================================
      ! INIT_POPS begins here!
      !=================================================================

      ! Return if we have already allocated arrays
      IF ( IS_INIT ) RETURN

      !=================================================================
      ! Allocate arrays
      !=================================================================
      ALLOCATE( COSZM( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COSZM' )
      COSZM = 0d0

      ALLOCATE( TCOSZ( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCOSZ' )
      TCOSZ = 0d0

      ALLOCATE( TTDAY( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TTDAY' )
      TTDAY = 0d0

      ! Allocate ZERO_DVEL if we use non-local PBL mixing or
      ! if drydep is turned off 
      IF ( LNLPBL .OR. (.not. LDRYD) ) THEN
         ALLOCATE( ZERO_DVEL( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ZERO_DVEL' )
         ZERO_DVEL = 0d0
      ENDIF

      !=================================================================
      ! Done
      !=================================================================

      ! Reset IS_INIT, since we have already allocated arrays
      IS_INIT = .TRUE.
      
      ! Return to calling program
      END SUBROUTINE INIT_POPS

!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CLEANUP_POPS
!
! !DESCRIPTION: Subroutine CLEANUP_POPS deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_POPS
!

! !INPUT PARAMETERS: 
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version
!
! !REMARKS:
! (1) Based initially on INIT_MERCURY from MERCURY_MOD (eck, 9/20/10)
!
!EOP
!------------------------------------------------------------------------------
!BOC

      IF ( ALLOCATED( COSZM    ) ) DEALLOCATE( COSZM   )     
      IF ( ALLOCATED( TCOSZ    ) ) DEALLOCATE( TCOSZ   )
      IF ( ALLOCATED( TTDAY    ) ) DEALLOCATE( TTDAY   )
      IF ( ALLOCATED( ZERO_DVEL) ) DEALLOCATE( ZERO_DVEL )

      END SUBROUTINE CLEANUP_POPS
!EOC
!------------------------------------------------------------------------------
      END MODULE POPS_MOD
