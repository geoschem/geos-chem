!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: linoz_mod
!
! !DESCRIPTION: Module LINOZ\_MOD contains routines to perform the Linoz 
!  stratospheric ozone chemistry.
!\\
!\\
! !INTERFACE:
!
      MODULE LINOZ_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PRIVATE DATA MEMBERS:
!
      REAL*8, ALLOCATABLE :: TLSTT(:,:,:)
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_LINOZ
      PUBLIC  :: DO_LINOZ    
      PUBLIC  :: INIT_LINOZ
      PUBLIC  :: LINOZ_READ
!
! !PRIVATE MEMBER FUNCTIONS:
! 
      PRIVATE :: LINOZ_CHEM3
      PRIVATE :: LINOZ_STRATL
      PRIVATE :: LINOZ_STRT2M
      PRIVATE :: LINOZ_SOMLFQ
      PRIVATE :: LINOZ_INTPL
!
! !REMARKS:
!  LINOZ Climatology:
!  ============================================================================
!  The LINOZ stratospheric chemistry tables for ozone consist of:
!                                                                             .
!  7 tables, each a function of:
!    12 months,
!    18 latitudes (-85 to 85 in 10 deg. increments)
!    25 altitudes ( z*=10-58 km in 2 km increments)
!                                                                             .
!  The 7 data fields are:
!    1- ozone (Logan climatology), v/v
!    2- Temperature climatology, K
!    3- Column ozone climatology, Logan ozone integrated above box, DU
!    4- ozone (P-L) for climatological ozone, v/v/s
!    5- d(P-L) / dO3, 1/s
!    6- d(P-L) / dT, v/v/s/K
!    7- d(P-L) / d(column O3), v/v/s/DU
!                                                                             .
!  Implementation notes:
!  ============================================================================
!  Dylan Jones (dbj@atmosp.physics.utoronto.ca) wrote:
!                                                                             .
!    Testing this code [in v8-02-04] was more difficult that I thought. 
!    I began by trying to compare the output of v8-02-04 with our previous 
!    runs with v8-02-01.  I accounted for the changes in the transport_mod.f 
!    and I tried to undo the changes in when the diagnostics are archived in 
!    v8-02-04, but I was still getting large differences between v8-02-04 
!    and v8-02-01. I finally gave up on this since I may have made a mistake
!    in reverting to the old way of doing the diagnostics in v8-02-04.  In 
!    the end I took the new linoz code from v8-02-04 and used it in v8-02-01. 
!    I ran two GEOS-5 full chemistry simulations for 2007 and the output 
!    were consistent over the full year.  
!                                                                             .
!    I think that it is safe to release [Linoz in v8-02-04].  However, we 
!    should acknowledge that it was [only] tested in v8-02-01, since I was 
!    not able to assess the quality of the output in v8-02-04. 
!                                                                             .
!  Bob Yantosca (yantosca@seas.harvard.edu) wrote:
!                                                                             .
!     We have also modified the code for use within the GEOS-5 GCM.  We now
!     declare the TPARM array as part of the Input_Opt object.  The LINOZ
!     climatology ASCII file is now read on the root CPU and MPI-broadcasted
!     to the non-root CPUs.  Also, the INIT_LINOZ routine is now called
!     not on the first chemistry timestep but rather in the initialization
!     phase at the start of the run. (bmy, 3/18/13)

! REVISION HISTORY:
!  23 Mar 2000 - P. Cameron-Smith    - Initial version adapted heavily 
!                                      from McLinden's original file.
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!  28 May 2009 - D. Jones            - Further modifications
!  18 Nov 2009 - D. Jones            - Further modifications
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  15 Mar 2013 - R. Yantosca - Now use fields from Input_Opt and made other
!                              modifications for GIGC interface to GEOS-5 GCM
!  18 Mar 2013 - R. Yantosca - Comment out STRAT_INIT, it's not called
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
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
! !IROUTINE: do_linoz
!
! !DESCRIPTION: Subroutine DO\_LINOZ is the main driver for the Linoz
!  stratospheric Ozone chemistry package.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_LINOZ( am_I_Root, Input_Opt, 
     &                     State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_MONTH
      USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!  10 Jun 2012 - L. Murray - Move call to DO_LINOZ from transport code to 
!                chemistry code, so the use of DT_TS_CHEM is now correct.
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  18 Mar 2013 - R. Yantosca - Now accept Input_Opt, RC as arguments
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! SAVEd scalars
      INTEGER, SAVE :: LASTMONTH = -99

      ! Non-SAVEd sclars
      INTEGER       :: MONTH
      REAL*8        :: NSCHEM

      !=================================================================
      ! DO_LINOZ begins here!
      !=================================================================

      ! Assume success
      RC    = GIGC_SUCCESS

      ! Current month
      MONTH = GET_MONTH()

      ! if new month, get new parameters?
      IF ( MONTH /= LASTMONTH ) THEN
         CALL LINOZ_STRATL( am_I_Root, Input_Opt, RC )
         LASTMONTH =  MONTH
      ENDIF

      ! Linoz needs time step in seconds
      NSCHEM = GET_TS_CHEM() * 60D0 

      ! Call the Linoz chemistry
      CALL LINOZ_CHEM3( NSCHEM,    am_I_Root, Input_Opt, 
     &                  State_Met, State_Chm, RC )

      END SUBROUTINE DO_LINOZ    
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: linoz_chem3
!
! !DESCRIPTION: Subroutine LINOZ\_CHEM3 applies linearized chemistry based on 
!  tables from PRATMO model using climatological T, O3, time of year
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LINOZ_CHEM3( DTCHEM,    am_I_Root, Input_Opt, 
     &                        State_Met, State_Chm, RC         )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2   
      USE TRACERID_MOD                          
      USE TROPOPAUSE_MOD,     ONLY : GET_TPAUSE_LEVEL 
      USE TROPOPAUSE_MOD,     ONLY : GET_MAX_TPAUSE_LEVEL
      USE PRESSURE_MOD,       ONLY : GET_PEDGE
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
!
! !INPUT PARAMETERS:
!
      REAL*8,         INTENT(IN)    :: DTCHEM      ! Time step [seconds]
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  Replace fields from tracer_mod.f with fields from Input_Opt.  When we
!  use GEOS-Chem within the GEOS-5 GCM, the fields within Input_Opt will
!  be read on the root CPU and MPI-broadcasted to all other CPUs.
!
! !REVISION HISTORY: 
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!  18 Nov 2009 - D. Jones            - For now, set tagged stratospheric 
!                                      tracer to total O3 in the overworld 
!                                      to avoid issues with spin ups
!  08 Feb 2010 - R. Yantosca         - Deleted obsolete local variables
!  22 Oct 2010 - R. Yantosca         - Added OMP parallel loop
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  14 Mar 2013 - M. Payer    - Replace Ox with O3 as part of removal of NOx-Ox
!                              partitioning
!  18 Mar 2013 - R. Yantosca - Now accept Input_Opt, RC as arguments
!  19 Mar 2013 - R. Yantosca - Now copy Input_Opt%TCVV(1:N_TRACERS)
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!  06 Nov 2013 - R. Yantosca - Now activate the parallel loop.  Unit testing
!                              revealed that this no longer causes errors.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER :: IM,      JM,         LM,   N_TRACERS
      INTEGER :: I,       J,          L,    N
      INTEGER :: K,       M,          LBOT, L_OVERWRLD
      INTEGER :: NTRACER, NUM_TRACER, LPOS, ITRC
      REAL*8  :: CLIMO3,  CLIMPML,    DCO3, DERO3, DERTMP
      REAL*8  :: DERCO3,  DMASS,      DTMP, SSO3

      ! Arrays
      REAL*8  :: DCOLO3  (IIPAR,JJPAR,LLPAR                    )
      REAL*8  :: COLO3   (IIPAR,JJPAR,LLPAR                    )
      REAL*8  :: OUT_DATA(IIPAR,JJPAR,LLPAR                    )
      REAL*8  :: TCVV    (                  Input_Opt%N_TRACERS)

      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)
      
      !=================================================================
      ! LINOZ_CHEM3 begins here!
      !=================================================================

      ! Assume success
      RC         = GIGC_SUCCESS

      ! Dimensions
      IM         = IIPAR
      JM         = JJPAR
      LM         = LLPAR
      N_TRACERS  = Input_Opt%N_TRACERS
      TCVV       = Input_Opt%TCVV(1:N_TRACERS)

      ! Initialize GEOS-Chem tracer array [kg] from State_Chm object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ) 
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! When we are connecting to the GEOS-5 GCM, we don't define
      ! L_OVERWRLD.  We cannot know the maximum extent of the
      ! tropopause in the GCM; we instead have to diagnose it at
      ! every timestep by comparing the pressure at a grid box to
      ! the tropopause pressure.  (bmy, 3/18/13)
      !-----------------------------------------------------------------
#else
      !-----------------------------------------------------------------
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! Current practice in the std GEOS-Chem is to call function
      ! GET_MAX_TPAUSE_LEVEL (from tropopause_mod.F) to return the
      ! highest possible tropopause level.  This is only required
      ! when doing a tagged Ox simulation. (bmy, 3/18/13)
      !-----------------------------------------------------------------
      L_OVERWRLD = GET_MAX_TPAUSE_LEVEL()
#endif

      !=================================================================
      ! Select the proper tracer number to store O3 into, depending on
      ! whether this is a full chemistry run or a tagged Ox run.
      ! If tagged Ox, tracer 2 should be the stratospheric tracer.  (dbj)
      !=================================================================
      IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
         NUM_TRACER = 1
      ELSE
         IF ( Input_Opt%ITS_A_TAGOX_SIM ) THEN
            IF ( N_TRACERS > 1 ) THEN
               NUM_TRACER = 2
            ELSE
               NUM_TRACER = 1
            ENDIF
         ELSE
            ! All other simulations don't use O3...print error message
            IF ( am_I_Root ) THEN
               WRITE( 6, '(a)' ) 'This simulation does not use O3!!'
               WRITE( 6, '(a)' ) 'STOP in linoz_chem3.f!'
            ENDIF
            STOP
         ENDIF
      ENDIF
      
      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) 
 100     FORMAT( '     - LINOZ_CHEM3: Doing LINOZ' )
      ENDIF

      !=================================================================
      ! Perform stratospheric chemistry
      !=================================================================

      ! **** note dbj: check STT(I,J,20:LLPAR,NTRACER) = with trop level
      ! ****         : check DMASS
      DO ITRC = 1, NUM_TRACER

         IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
            NTRACER = IDTO3
         ELSE
            NTRACER = ITRC
         ENDIF

         ! Start at top layer and continue to lowest layer for strat. chem
         OUT_DATA = 0d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,       J,     LBOT,   LPOS,   L    )
!$OMP+PRIVATE( CLIMPML, DERO3, CLIMO3, DERCO3, DCO3 )
!$OMP+PRIVATE( DERTMP,  DTMP,  SSO3,   DMASS        )
         DO J = 1, JM
         DO I = 1, IM

            LBOT = GET_TPAUSE_LEVEL( I, J, State_Met ) + 1
            LPOS = 1
            DO WHILE (GET_PEDGE(I,J,LPOS+1) .GE. 0.3D0)
               LPOS = LPOS +1
            ENDDO
            LPOS = LPOS-1

#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
            !-----------------------------------------------------------
            !  %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
            !
            ! When we are connecting to the GEOS-5 GCM, we don't define
            ! L_OVERWRLD.  We cannot know the maximum extent of the
            ! tropopause in the GCM; we instead have to diagnose it at
            ! every timestep by comparing the pressure at a grid box to
            ! the tropopause pressure.  (bmy, 3/18/13)
            !-----------------------------------------------------------
#else
            !-----------------------------------------------------------
            !           %%%%% TRADITIONAL GEOS-Chem %%%%%
            !
            ! Current practice in the std GEOS-Chem is to call function
            ! GET_MAX_TPAUSE_LEVEL (from tropopause_mod.F) to return the
            ! highest possible tropopause level.  This is only required
            ! when doing a tagged Ox simulation. (bmy, 3/18/13)
            !-----------------------------------------------------------

            ! dbj: for now, set tagged stratospheric tracer to total 
            ! O3 in the overworld to avoid issues with spin ups
            IF ( Input_Opt%ITS_A_TAGOX_SIM ) THEN
              STT(I,J,(L_OVERWRLD+1):LLPAR,NTRACER) =
     &                STT(I,J,(L_OVERWRLD+1):LLPAR,1)
            ENDIF
#endif

            ! Loop over levels
            DO L = LM, LBOT, -1

               ! Skip if tracer is negative
               IF ( STT(I,J,L,NTRACER) .LE. 0.D0 ) CYCLE

               ! calculate ozone column above box (and save)
               ! dcolo3 = ozone column (in DU) in given layer
               ! colo3 =  ozone column above layer + half of
               ! column in layer

               ! bdf stt is in v/v, make conversion to DU
               if (l.eq.lm) then !top model layer
                  dcolo3(i,j,l) = (stt(i,j,l,NTRACER) *
     &                 State_Met%AD(I,J,L) /
     &                 TCVV(NTRACER))      / GET_AREA_CM2(I,J,L) *
     &                 6.022d23 / (28.97/TCVV(NTRACER)*1d-3) / 2.687d16
                  colo3(i,j,l) = dcolo3(i,j,l)*0.5
               else
                  dcolo3(i,j,l) = (stt(i,j,l,NTRACER) * 
     &                 State_Met%AD(I,J,L) / 
     &                 TCVV(NTRACER))      / GET_AREA_CM2(I,J,L) *
     &                 6.022d23 / (28.97/TCVV(NTRACER)*1d-3) / 2.687d16
                  colo3(i,j,l) = colo3(i,j,l+1) +
     &                 (dcolo3(i,j,l)+dcolo3(i,j,l+1))*0.5
               endif
               out_data(i,j,l) = colo3(i,j,l)

               ! ++++++ climatological P-L:   ++++++          
               climpml=tlstt(j,l,4)      ! Climatological P-L = (P-L)^o

               ! ++++++ local ozone feedback: ++++++ 
               dero3=tlstt(j,l,5)               ! Derivative w.r.t. O3. 
                                                !  dero3=-1/(time constant)
               IF (dero3.EQ.0) CYCLE            ! Skip Linoz if lifetime 
                                                !  is infinite.
               climo3=tlstt(j,l,1)              ! Climatological O3 = f^o
               derco3=tlstt(j,l,7)              ! Derivative w.r.t. Column O3
               dco3=(colo3(i,j,l)-tlstt(j,l,3)) ! deviation from o3 climatology.
               ! ++++++ temperature feedback: ++++++ 
               dertmp=tlstt(j,l,6)              ! Derivative w.r.t. Temperature
               dtmp=(State_Met%T(I,J,L) -       ! Deviation in Temperature 
     &                tlstt(j,l,2))             !  from climatology.

               ! ++++++ calculate steady-state ozone: ++++++ 
               sso3=climo3 - (climpml+dtmp*dertmp+dco3*derco3)/dero3 

               ! ++++++ change in ozone mass due to chemistry: ++++++
               !ssO3 = f^*
               dmass=(sso3-STT(I,J,L,NTRACER))*(1.0-exp(dero3*dtchem))

               ! ++++++ update ozone mass ++++++ 
               ! LINOX valid only up to 58 km, so do not use above 0.3 hPa
               ! dbj: impose exponential fall off of mixing ratio 
               ! between 0.3 and 0.01 hPa (with fall off of a scale height)
               IF (GET_PEDGE(I,J,L) .LE. 0.3D0) THEN
                  STT(I,J,L,NTRACER) = 
     &                (GET_PCENTER(I,J,L)/GET_PCENTER(I,J,LPOS-1)) * 
     &                STT(I,J,LPOS-1,NTRACER)
               ELSE
                  STT(I,J,L,NTRACER) = STT(I,J,L,NTRACER)+DMASS 
               ENDIF

            ENDDO       ! loop over L

         ENDDO          ! loop over I

      ENDDO             ! loop pver J
!$OMP END PARALLEL DO

      !write our calculated column o3 maximum
      !write(6,*) 'max of columns= ',maxval(out_data)

      ENDDO

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE LINOZ_CHEM3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: linoz_stratl
!
! !DESCRIPTION: Subroutine LINOZ\_STRATL performs a monthly fixup of chemistry 
!  parameters for the Linoz stratospheric ozone chemistry.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LINOZ_STRATL( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GRID_MOD,           ONLY : GET_YMID 
      USE TIME_MOD,           ONLY : GET_MONTH
      USE PRESSURE_MOD
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  Replace size fields NLAT_LINOZ etc. with fields from Input_Opt.  When we
!  use GEOS-Chem within the GEOS-5 GCM, the fields within Input_Opt will
!  be read on the root CPU and MPI-broadcasted to all other CPUs.
!                                                                             .
!  The LINOZ climatology array is Input_Opt%LINOZ_TPARM(25,18,12,N),
!  which has the following dimensions
!    * 25 layers from 58 km to 10 km by 2 km intervals
!    * 18 latitudes (85S, 75S, ...85N)
!    * 12 months
!    *  N fields (currently N=7)
!
! !REVISION HISTORY: 
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  18 Mar 2013 - R. Yantosca - Accept Input_Opt, RC arguments.  Use fields
!                              from Input_Opt to facilitate GC/GEOS-5 interface
!  18 Mar 2013 - R. Yantosca - Cosmetic changes, updated comments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER           :: I, IM, IM1, IM2
      INTEGER           :: J, JM, JJ,  JXXX
      INTEGER           :: K, L,  LR,  MONTH, N

      ! Dimension sizes
      INTEGER           :: NFIELDS_LINOZ
      INTEGER           :: NLAT_LINOZ
      INTEGER           :: NLEVELS_LINOZ
      INTEGER           :: NMONTHS_LINOZ

      ! Arrays
      INTEGER           :: JLATMD(JJPAR                  )
      REAL*8            :: P0L   (LLPAR+1                )
      REAL*8            :: STRT0L(LLPAR+1                )
      REAL*8            :: STRT1L(LLPAR+1                )
      REAL*8            :: STRT2L(LLPAR+1                )
      REAL*8            :: STRTX (Input_Opt%LINOZ_NLEVELS)
      REAL*8            :: YSTRT (Input_Opt%LINOZ_NLAT   )

      ! Month names
      CHARACTER(LEN=3)  :: CMONTH(12) = (/'jan', 'feb', 'mar', 'apr',
     &                                    'may', 'jun', 'jul', 'aug',
     &                                    'sep', 'oct', 'nov', 'dec'/)
!
! !DEFINED PARAMETERS:
!
      REAL*8, PARAMETER :: PSF = 1010D0   ! Surface pressure [hPa]

      !=================================================================
      ! LINOZ_STRATL begins here!
      !=================================================================

      ! Initialization
      IM             = IIPAR
      JM             = JJPAR
      MONTH          = GET_MONTH()
      NFIELDS_LINOZ  = Input_Opt%LINOZ_NFIELDS
      NLAT_LINOZ     = Input_Opt%LINOZ_NLAT
      NLEVELS_LINOZ  = Input_Opt%LINOZ_NLEVELS
      NMONTHS_LINOZ  = Input_Opt%LINOZ_NMONTHS

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, '(a)' ) REPEAT( '#', 79 )
         WRITE( 6,  50   ) CMONTH(MONTH)
         WRITE( 6, '(a)' ) REPEAT( '#', 79 )
 50      FORMAT( '# Interpolating Linoz fields for ', a )
      ENDIF

!------------------------------------------------------------------------------
! THE FOLLOWING LINEAR INTERPOLATION OPTION IS NOT USED, BUT CAN BE RESTORED
! BY UNCOMMENTING THE APPROPRIATE LINES OF CODE BELOW (bmy, 3/18/13)
!
!      ! get weights for month interpolation
!      do i=1,nmonths_linoz
!        jdmc(i) = jdofm(i+1) - (jdofm(i+1)-jdofm(i))/2
!      enddo
!
!      im1=0
!      do i=1,nmonths_linoz
!        if (jdmc(i).lt.jday) then
!          im1=i
!        endif
!      enddo
!      if (im1.eq.0) then
!        im1=nmonths_linoz
!        im2=1
!        wm1=(jdmc(im2)-jday)*1.0/(jdmc(im2)-(jdmc(im1)-365.0))
!      elseif (im1.eq.nmonths_linoz) then
!        im2=1
!        wm1=(jdmc(im2)+365.0-jday)/(jdmc(im2)+365.0-jdmc(im1))
!      else
!        im2=im1+1
!        wm1=(jdmc(im2)-jday)*1.0/(jdmc(im2)-jdmc(im1))
!      endif
!      wm2=1.0-wm1
!
!!      write(6,*)iday,jday,' weights: ',wm1,wm2
!!      write(6,*)'months: ',im1,im2,month
!!      write(6,*)'between: ',jdmc(im1),jdmc(im2)
!------------------------------------------------------------------------------


      ! YSTRT(J) are the latitudes (-85, -75, -65, ... +65, +75, +85)
      ! contained in the LINOZ climatology
      YSTRT(1) = -85.d0   
      do J = 2, NLAT_LINOZ
        YSTRT(J) = YSTRT(J-1) + 10.d0
      enddo

      ! JLATMD(J) is the nearest-neighbor LINOZ data column corresponding
      ! to each GEOS-Chem latitude index J. (dbj, 6/25/03)
      DO J = 1, JJPAR
         JXXX      = INT( 0.1d0 * GET_YMID(1,J,1) + 10.d0 )
         JLATMD(J) = MIN( 18,     MAX( 1, JXXX )          )
      ENDDO

      ! P0L are the pressure at the level edges of each GEOS-Chem
      ! grid box, assuming a surface pressure of 1010 hPa (dbj, 6/25/03)
      DO L = 1, LLPAR+1
         P0L(L) = GET_AP(LLPAR+2-L) + ( GET_BP(LLPAR+2-L) * PSF )
      ENDDO

      !=================================================================
      ! Lookup data in the LINOZ climatology 
      !=================================================================

      ! Loop over the # of fields in the LINOZ climatology
      DO N = 1, NFIELDS_LINOZ

!------------------------------------------------------------------------------
! INTERPOLATION BETWEEN LATITUDES IS CURRENTLY NOT USED, BUT CAN BE RESTORED
! BY UNCOMMENTING THE APPROPRIATE LINES OF CODE BELOW (bmy, 3/18/13)
!
! ***** Interpolation between latitudes is not currently used {PJC} *****
!c----- interpolating along latitude, from TPAR2 to STRTXY
!        do K = 1,nlevels_linoz
!          do J = 1,nlat_linoz
!c            TPAR2(K,J) = TPARM(K,J,MONTH,N)
!            TPAR2(K,J) = TPARM(K,J,im1,N)
!          enddo
!        enddo
!        call LINOZ_INTPL(nlevels_linoz,NLAT_LINOZ,JPAR,JM,YSTRT,YDGRD,
!     &       TPAR2,STRTXY1)
!        do K = 1,nlevels_linoz
!          do J = 1,nlat_linoz
!            TPAR2(K,J) = TPARM(K,J,im2,N)
!          enddo
!        enddo
!        call LINOZ_INTPL(nlevels_linoz,NLAT_LINOZ,JPAR,JM,YSTRT,YDGRD,
!     &       TPAR2,STRTXY2)
!------------------------------------------------------------------------------

         ! Loop over GEOS-Chem latitudes
         DO J = 1, JM

            ! JJ is the index of the nearest LINOZ data column
            ! corresponding to GEOS-Chem latitude index J
            JJ = JLATMD(J)

            ! Loop over the # of levels in the LINOZ climatology
            DO K = 1, NLEVELS_LINOZ

               ! STRTX(K) is the column of data from the LINOZ climatology
               ! for the given month, latitude (JJ), level (K), and species
               ! (N).  One of the following interpolation options may be used.

!------------------------------------------------------------------------------
! THE FOLLOWING OPTIONS ARE CURRENTLY NOT USED, BUT CAN BE RESTORED BY
! UNCOMMENTING THE APPROPRIATE LINES OF CODE BELOW (bmy, 3/18/13)
!
!               ! linearly interpolate in latitude and month
!               STRTX(K) = STRTXY1(K,J)*wm1 + STRTXY2(K,J)*wm2
!
!               ! linearly interpolate in latitude, single month
!               STRTX(K) = STRTXY2(K,J)
!
!               ! nearest latitude, linearly interpolate in month
!               !STRTX(K) = Input_Opt%LINOZ_TPARM(K,JJ,im1,N)*wm1
!     &                   + Input_Opt%LINOZ_TPARM(K,JJ,im2,N)*wm2
!------------------------------------------------------------------------------

               ! Nearest neighbor, no interpolation
               STRTX(K) = Input_Opt%LINOZ_TPARM(K,JJ,MONTH,N)
            ENDDO

            ! *PJC* Interpolate and calculate moments of column distribution 
            CALL LINOZ_STRT2M( am_I_Root, Input_Opt, LLPAR,
     &                         STRTX,     P0L,       STRT0L,
     &                         STRT1L,    STRT2L,    RC      )

            ! Store loss freq/yields & moments in TLSTT/SWT/SWW 
            ! for exact CTM layers LM down
            ! Order reversed from C.McLinden version {PJC}
            DO LR = 1,LLPAR  
               TLSTT(J,LR,N) = STRT0L(LLPAR+1-LR)
            ENDDO
         ENDDO     ! loop over J
      ENDDO      ! loop over N

      END SUBROUTINE LINOZ_STRATL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: linoz_strt2m
!
! !DESCRIPTION: Subroutine LINOZ\_STRT2M interpolates quantities from the
!  LINOZ vertical grid to the GEOS-Chem vertical grid.  It also computes
!  the 1st \& 2nd moments of the distribution.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LINOZ_STRT2M( am_I_Root, Input_Opt, NSTRT,
     &                         STRTX,     P0L,       STRT0L,
     &                         STRT1L,    STRT2L,    RC      )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      !--------------------------------------------------------------------
      ! am_I_Root : Are we on the root CPU?
      ! Input_Opt : Input Options object
      ! NSTRT     : # of levels in the GEOS-Chem grid (= LLPAR)
      ! STRTX     : Quantity on the LINOZ vertical grid
      !             (i.e. fields #1-7 of the LINOZ climatology)
      ! P0L       : Pressure edges on the GEOS-Chem grid
      !--------------------------------------------------------------------
      LOGICAL,        INTENT(IN)  :: am_I_Root
      TYPE(OptInput), INTENT(IN)  :: Input_Opt
      INTEGER,        INTENT(IN)  :: NSTRT
      REAL*8,         INTENT(IN)  :: STRTX(Input_Opt%LINOZ_NLEVELS)
      REAL*8,         INTENT(IN)  :: P0L(LLPAR+1)
!
! !OUTPUT PARAMETERS:
!
      !--------------------------------------------------------------------
      ! STRT0L    : 0th moment of distribution, on GEOS-Chem grid edges
      ! STRT1L    : 1st moment of distribution, on GEOS-Chem grid edges
      ! STRT2L    : 2nd moment of distribution, on GEOS-Chem grid edges
      ! RC        : Success or failure?
      !--------------------------------------------------------------------
      REAL*8,         INTENT(OUT) :: STRT0L(LLPAR+1)
      REAL*8,         INTENT(OUT) :: STRT1L(LLPAR+1)
      REAL*8,         INTENT(OUT) :: STRT2L(LLPAR+1)
      INTEGER,        INTENT(OUT) :: RC
!
! !REMARKS:
!  Comments from Chris McLinden to Peter Cameron-Smith:
!  ===========================================================================
!  CALL SOMLFQ(P1,P2,F0,F1,F2,PS,F,NL)
!  - P1,P2 are the pressure EDGES for the CTM layer onto which the
!    coefficients will be mapped. [P1>P2 I believe {PJC}]
!  - F0,F1,F2 are the CTM layer vertical moments determined in SOMLFQ
!  - PS are the pressure layer edges of the original [ie Linox] grid
!  - F is the column of coefficients (on the original grid); note
!    F is flipped relative to STRTX and since the coefficients begin
!    at z*=10, F(1)=F(2)=...=F(5)=0
!  - NL is 30; size of F()
!                                                                             .
!   The box model calculations were performed at z*=10km, 12km, ... and
!   so these would represent the centres with the corresponding edges at
!   9,11km ; 11,13km; ...
!   PS() represents the edges (although PS(1) is set to 1000mb).
!   The first few values are:
!     PS(1)=1000
!     PS(2)=874.947105    (note PS(2) is not quite 1000 exp(-1/16) as the
!     PS(3)=656.117767     the average pressure is used - not the pressure
!     PS(4)=492.018914     at the average z*)
!     PS(5)=368.96213
!     PS(6)=276.68257
!     PS(7)=207.48266
!     ...
!     PS(30)=0.276682568
!     PS(31)=0.0
!                                                                             .
!     F(1) spans PS(1)-PS(2)
!     F(2) spans PS(2)-PS(3)
!     ...
!     F(30) spans PS(30)-PS(31)

! !REVISION HISTORY:
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!  18 Mar 2013 - R. Yantosca - Now pass am_I_root, Input_Opt, RC arguments
!  18 Mar 2013 - R. Yantosca - Rearrange argument list so that inputs are
!                              all listed before outputs
!  18 Mar 2013 - R. Yantosca - Now make NL a local varable and set its
!                              value from Input_Opt%LINOZ_NLEVELS
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER :: NL,  NCBOX, L,      K,    NX
      REAL*8  :: P1,  P2,    F0,     F1
      REAL*8  :: F2,  XPSD,  XPSLM1, XPSL

      ! Arrays
      REAL*8  :: PS(Input_Opt%LINOZ_NLEVELS+6)  ! Old code: PS(NL+1)
      REAL*8  :: F (Input_Opt%LINOZ_NLEVELS+5)  ! Old code: PS(NL)

      !=================================================================
      ! Initialization
      !=================================================================

      ! Assume success
      RC         = GIGC_SUCCESS

      ! Now make NX a local variable, since we already pass
      ! its value via the Input_Opt object (bmy, 3/13/18)
      NX         = Input_Opt%LINOZ_NLEVELS

      ! Now make NL a local variable instead of a parameter, because
      ! we need to construct its value usign Input_Opt (bmy, 3/18/13)
      NL         = Input_Opt%LINOZ_NLEVELS + 5

      !=================================================================
      ! Set up std z* atmosphere: p = 1000 * 10**(-z*/16 km)
      !
      ! Assume that stratospheric chemical parameters always start at
      ! 52 km (N=27).  Scan downward from 52 km to 14 km (NX=20) by
      ! 2 km.
      !
      ! 58 km (N=30) scan downward from 58 km to 10 km (NX=25) by 2 km
      ! intervals, constant >58km
      !
      !  N.B. F(@30km) assumed to be constant from 29-31 km (by mass)
      !=================================================================

      XPSD       = 10.D0 **(-0.125D0)
      XPSLM1     = 1000.D0
      PS(1)      = 1000.D0
      DO L = 2,NL
         XPSL    = XPSLM1 *XPSD
         PS(L)   = 0.5D0 *(XPSLM1 +XPSL)
         XPSLM1  = XPSL
      ENDDO
      PS(NL+1)   = 0.D0
      DO L = 1,NL-NX
         F(L)     = 0.D0
      ENDDO

      ! K=1 is at the top of atmosphere
      DO K = 1,NX
        F(NL+1-K)= STRTX(K) !STRTX has increasing preasure. {PJC}
      ENDDO

      DO K = 1,NSTRT
         P1       = P0L(K+1)
         P2       = P0L(K)
         CALL LINOZ_SOMLFQ(P1,P2,F0,F1,F2,PS,F,NL)
         STRT0L(K)= F0
         STRT1L(K)= F1
         STRT2L(K)= F2
      ENDDO

      END SUBROUTINE LINOZ_STRT2M
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: linoz_somlfq
!
! !DESCRIPTION: subroutine LINOZ\_SOMLFQ calculates loss freq moments from a 
!  set of loss frequencies at std z*, given a CTM model interval pressure
!  range: P1 > P2 (decreasing up)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LINOZ_SOMLFQ(P1,P2,F0,F1,F2,PS,F,NL)
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)  :: NL
      REAL*8,  INTENT(IN)  :: F(NL)
      REAL*8,  INTENT(IN)  :: PS(NL+1)
      REAL*8,  INTENT(IN)  :: P1
      REAL*8,  INTENT(IN)  :: P2
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: F0
      REAL*8,  INTENT(OUT) :: F1
      REAL*8,  INTENT(OUT) :: F2
!
! REMARKS:
! The pressure levels BETWEEN z* values are:
!      PS(i) > PS(i+1) bounds z*(i)
!                                                                             .
! NL:  z* levels, ==> PS(NL+1) = 0  (extrapolate chemical loss to top)
!      Z1 = 16.D0*LOG10(1000.D0/P1)
!      Z2 = 16.D0*LOG10(1000.D0/P2)
!                                                                             .
! The MOMENTS for a square-wave or 'bar': F(x)=f0  b<=x<=c, =0.0 else
!      S0 =   f0 (x)                      [from x=b to x=c]
!      S1 = 3 f0 (x^2 - x)                [from x=b to x=c]
!      S2 = 5 f0 (2x^3 - 3x^2 + x)        [from x=b to x=c]
!
! !REVISION HISTORY: 
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!  19 Mar 2013 - R. Yantosca - P1, P2 are now declared as INTENT(IN)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer  I
      real*8   XB,XC,PC,PB,THIRD,sgnf0

      F0     = 0.D0
      F1     = 0.D0
      F2     = 0.D0
      DO I = 1,NL
         PC   = MIN(P1,PS(I))
         PB   = MAX(P2,PS(I+1))
         IF (PC .GT. PB)  THEN

            ! have condition:  P1>=PC > PB>=P2, 0<=XB < XC<=1
            XC = (PC-P2)/(P1-P2)
            XB = (PB-P2)/(P1-P2)

            ! assume that the loss freq, F, is constant over interval [XB,XC],
            ! F0: (c-b),   
            ! F1: 6((c2-c)-(b2-b)),  
            ! F2: 5((2c3-3c2+c)-(2b3-3b2+b))
            ! calculate its contribution to the moments in the interval [0,1]
            F0 = F0 +F(I) *(XC -XB)
            F1 = F1 +F(I) *3.D0 *((XC *XC -XC) - (XB *XB -XB))
            F2 = F2 +F(I) *5.D0 *
     &           ((XC+XC-1.D0)*(XC*XC -XC) - (XB+XB-1.D0)*(XB*XB -XB))
         ENDIF
      ENDDO

      ! RESTRAIN moments: force monotonicity & positive at min end pt

      ! cam: tables can be + or -
      if (f0.ne.0.0) then
         sgnf0=f0 / abs(f0)
      else
         sgnf0=1.0
      endif
      f0=abs(f0)
      
      !F0 = MAX(F0, 0.D0)
      THIRD = 1.D0/3.D0
      IF (F2 .GT. 0.D0)  THEN


         ! do not allow reversal of curvature: F2 > 0
         F2   = MIN(F2, ABS(F1)*THIRD, 5.D-1*F0)
         IF (F1 .LT .0.D0)  THEN
            F1 = MAX(-(F0+F2), F1)
         ELSE
            F1 = MIN(+(F0+F2), F1)
         ENDIF
      ELSE

         ! F2 < 0 = curved down at ends, allow if F1 < F0
         F1  = MIN(F0,MAX(-F0,F1))
         F2  = MAX(F2,(ABS(F1)-F0),(-ABS(F1)*THIRD))
      ENDIF

      ! cam: apply sign
      f0=sgnf0 * f0
      f1=sgnf0 * f1
      f2=sgnf0 * f2

      END SUBROUTINE LINOZ_SOMLFQ
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: linoz_read
!
! !DESCRIPTION: Subroutine LINOZ\_READ reads the input data file for the
!  Linoz stratospheric ozone chemistry.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LINOZ_READ( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE FILE_MOD,           ONLY : IOERROR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE InquireMod,         ONLY : findFreeLun
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  LINOZ_READ is called from "main.f" at the start of the simulation.
!  LINOZ_READ will also call INIT_LINOZ to initialize the arrays.
!
! !REVISION HISTORY: 
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!  16 Oct 2009 - R. Yantosca         - Now use IU_FILE instead of IU_LINOZ
!  16 Oct 2009 - R. Yantosca         - Read file from DATA_DIR_1x1
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  15 Mar 2013 - R. Yantosca - Now call INIT_LINOZ from GIGC_Init_Extra
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      ! Scalars
      INTEGER            :: IU_FILE, K,     J
      INTEGER            :: M,       NTBLS, IOS
      REAL*8             :: TMAX,    TMIN

      ! Dimension sizes
      INTEGER            :: NFIELDS_LINOZ
      INTEGER            :: NLAT_LINOZ
      INTEGER            :: NLEVELS_LINOZ
      INTEGER            :: NMONTHS_LINOZ

      ! Strings
      CHARACTER(LEN=80)  :: HEADING, TITL1
      CHARACTER(LEN=255) :: ERRMSG
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! LINOZ_READ begins here!
      !=================================================================

      ! Assume success
      RC             = GIGC_SUCCESS

      ! Get fields from Input_Opt
      NFIELDS_LINOZ  = Input_Opt%LINOZ_NFIELDS
      NLAT_LINOZ     = Input_Opt%LINOZ_NLAT
      NLEVELS_LINOZ  = Input_Opt%LINOZ_NLEVELS
      NMONTHS_LINOZ  = Input_Opt%LINOZ_NMONTHS

      ! Filename
      FILENAME       = TRIM( Input_Opt%DATA_DIR_1x1 ) //
     &                 'Linoz_200910/Linoz_March2007.dat'

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - LINOZ_READ: Reading ', a )
      ENDIF

      ! Find a free file LUN
      IU_FILE = findFreeLUN()

      ! Define error message
      ERRMSG  = 'LINOZ_READ (in GeosCore/linoz_mod.F)'

      ! new std z*=2km levels from model:  z*=10,12,...(25*2)+8 km
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', 
     &               FORM='FORMATTED',      IOSTAT=IOS )

      ! Stop on error
      IF ( IOS /= 0 ) THEN
         CALL IOERROR( IOS, IU_FILE, ERRMSG )
      ENDIF

      ! Read header
      READ ( IU_FILE, '(a)' ) HEADING
      IF ( am_I_Root ) THEN
         WRITE(6,*) TRIM( HEADING )
      ENDIF

      ! Loop over # of fields in the LINOZ climatology
      DO NTBLS = 1, NFIELDS_LINOZ

         ! Zero min & max values
         TMIN = +1.d30
         TMAX = -1.d30

         ! Skip header line
         READ( IU_FILE, '(a)' ) TITL1

         ! Loop over # of months in the LINOZ climatology
         do M = 1, NMONTHS_LINOZ

            ! Loop over # of latitudes in the LINOZ climatology
            do J = 1, NLAT_LINOZ

               ! Read data into Input_Opt%LINOZ_TPARM
               READ( IU_FILE, '(20X,6E11.4/(8E11.4))', IOSTAT=IOS )
     &               ( Input_Opt%LINOZ_TPARM(K,J,M,NTBLS),
     &                 K=NLEVELS_LINOZ,1,-1                 )

               ! Stop on error
               IF ( IOS > 0 ) THEN
                  CALL IOERROR( IOS, IU_FILE, ERRMSG )
               ENDIF

               ! Loop over # of levels in the LINOZ climatology
               ! and compute the overall min & max
               do K = 1, NLEVELS_LINOZ
                  TMAX = MAX( TMAX, Input_Opt%LINOZ_TPARM(K,J,M,ntbls) )
                  TMIN = MIN( TMIN, Input_Opt%LINOZ_TPARM(K,J,M,ntbls) )
               enddo
            enddo
         enddo

         ! Write overall min & max
         write (6,912) TITL1,TMIN,TMAX
 912     FORMAT('  Linoz Data:  ',a80,1p,2e10.3)

      enddo

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, '(a)' ) '$$ Finished Reading Linoz Data $$'
         WRITE( 6, '(a)' )
      ENDIF

      ! Close the files
      CLOSE( IU_FILE ) 

      END SUBROUTINE LINOZ_READ
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: linoz_intpl
!
! !DESCRIPTION: Subroutine LINOZ\_INTPL does some kind of interpolation.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LINOZ_INTPL(KE,IE,ND,NE,XI,XN,YI,YN)
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)  :: KE
      INTEGER, INTENT(IN)  :: IE
      INTEGER, INTENT(IN)  :: ND
      INTEGER, INTENT(IN)  :: NE
      REAL*8,  INTENT(IN)  :: XI(IE) 
      REAL*8,  INTENT(IN)  :: XN(ND)
      REAL*8,  INTENT(IN)  :: YI(KE,IE)
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: YN(KE,ND)

!
! !REVISION HISTORY: 
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer  I,II,J,K
      real*8   CNST1,CNST2

      ! k=height; i=lat
      J       = 2
      do I = 1,NE
        if (XN(I) .gt. XI(1        ))  then
          if (XN(I) .lt. XI(IE))  then
            CNST1     = (XI(J) - XN(I)) / (XI(J) - XI(J-1))
            CNST2     = (XN(I) - XI(J-1)) / (XI(J) - XI(J-1))
            do K = 1,KE
              YN(K,I) = CNST1 * YI(K,J-1) + CNST2 * YI(K,J)
            enddo
            II    = min(I+1,NE)
            if (XN(II) .gt. XI(J))  J = min(IE,J+1)
          else
            do K = 1 ,KE
              YN(K,I) = YI(K,IE)
            enddo
          endif
        else
          do K = 1,KE
            YN(K,I)   = YI(K,1)
          enddo
        endif
        !write(6,*)i,(yn(k,i),k=1,ke)
      enddo

      END SUBROUTINE LINOZ_INTPL
!EOC
!----------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_linoz
!
! !DESCRIPTION: Subroutine INIT\_LINOZ allocates and zeroes the module arrays 
!  used in the Linoz stratospheric ozone algorithm.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_LINOZ( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD 
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  16 Oct 2009 - R. Yantosca - Initial version
!  18 Mar 2013 - R. Yantosca - Remove TPARM array since that is now carried
!                              within the Input_Opt object.
!  18 Mar 2013 - R. Yantosca - Accept am_I_Root, Input_Opt, RC arguments
!  14 Mar 2013 - M. Payer    - Replace Ox with O3 for full-chemistry simulation
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC = GIGC_SUCCESS

      ! Allocate TLSTT array
      ALLOCATE( TLSTT( JJPAR, 
     &                 LLPAR, 
     &                 Input_Opt%LINOZ_NFIELDS ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'TPARM' )
      TLSTT = 0d0

      END SUBROUTINE INIT_LINOZ
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_linoz
!
! !DESCRIPTION: Subroutine CLEANUP\_LINOZ deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_LINOZ
!
! !REVISION HISTORY: 
!  16 Oct 2009 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Deallocate arrays
      IF ( ALLOCATED( TLSTT ) ) DEALLOCATE( TLSTT )
      
      END SUBROUTINE CLEANUP_LINOZ
!EOC
      END MODULE LINOZ_MOD
