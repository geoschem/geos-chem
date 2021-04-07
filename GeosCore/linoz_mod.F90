!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: linoz_mod.F90
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
  USE ERROR_MOD, ONLY : DEBUG_MSG  ! Routine for debug output
  USE PRECISION_MOD                ! For GEOS-Chem Precision (fp, f4, f8)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: DO_LINOZ
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
!    runs with v8-02-01.  I accounted for the changes in the transport_mod.F90 
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
!  See https://github.com/geoschem/geos-chem for complete history
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
  SUBROUTINE DO_LINOZ( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_MONTH
    USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    REAL(fp)        :: NSCHEM

    !=================================================================
    ! DO_LINOZ begins here!
    !=================================================================

    ! Assume success
    RC    = GC_SUCCESS

    ! Current month
    MONTH = GET_MONTH()

    ! if new month, get new parameters?
    IF ( MONTH /= LASTMONTH ) THEN
       CALL LINOZ_STRATL( Input_Opt, State_Chm, State_Grid, RC )
       LASTMONTH =  MONTH
    ENDIF

    ! Linoz needs time step in seconds
    NSCHEM = GET_TS_CHEM()

    ! Call the Linoz chemistry
    CALL LINOZ_CHEM3( NSCHEM, Input_Opt, State_Chm, State_Grid, State_Met, RC )

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
  SUBROUTINE LINOZ_CHEM3( DTCHEM, Input_Opt, State_Chm, State_Grid, &
                          State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PhysConstants,      ONLY : AIRMW, AVO
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    REAL(fp),       INTENT(IN)    :: DTCHEM      ! Time step [seconds]
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
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
! !REVISION HISTORY:
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER  :: IM,      JM,         LM,   nAdvect
    INTEGER  :: I,       J,          L,    N,         NA
    INTEGER  :: K,       M,          LBOT, L_OVERWRLD
    INTEGER  :: NTRACER, NUM_TRACER, LPOS, ITRC
    REAL(fp) :: CLIMO3,  CLIMPML,    DCO3, DERO3,     DERTMP
    REAL(fp) :: DERCO3,  DMASS,      DTMP, SSO3

    ! Arrays
    REAL(fp) :: DCOLO3  (State_Grid%NX,State_Grid%NY,State_Grid%NZ )
    REAL(fp) :: COLO3   (State_Grid%NX,State_Grid%NY,State_Grid%NZ )
    REAL(fp) :: OUT_DATA(State_Grid%NX,State_Grid%NY,State_Grid%NZ )

    ! We need to define local arrays to hold corresponding values
    ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
    REAL(fp), POINTER :: Spc(:,:,:,:)

    ! Define local arrays to hold the previously global TLSTT array
    ! from the State_Chm object. (hplin, 1/22/19)
    REAL(fp), POINTER :: TLSTT (:,:,:,:)   ! IM, JM, LM, Input_Opt%LINOZ_NFIELDS

    ! Print debug output?
    LOGICAL       :: prtDebug

    ! SAVED scalars
    LOGICAL, SAVE :: FIRST = .TRUE.
    INTEGER, SAVE :: id_O3

    !=================================================================
    ! LINOZ_CHEM3 begins here!
    !=================================================================

    ! Print debug output?
    prtDebug   = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Assume success
    RC         = GC_SUCCESS

    ! Dimensions
    IM         = State_Grid%NX
    JM         = State_Grid%NY
    LM         = State_Grid%NZ
    nAdvect    = State_Chm%nAdvect
    L_OVERWRLD = 0

    ! Point to chemical species array [v/v dry air]
    Spc        => State_Chm%Species

    ! Point to the TLSTT array
    TLSTT      => State_Chm%TLSTT

    ! Look up the species ID of O3 only on the first call
    IF ( FIRST ) THEN
       id_O3   = Ind_('O3')
       FIRST   = .FALSE.
    ENDIF

    ! For the tagged O3 simulation only, get the highest
    ! level of the chemistry grid in the column at (I,J)
    IF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN
       L_OVERWRLD = MAXVAL( State_Met%ChemGridLev )
    ENDIF

    !=================================================================
    ! Select the proper tracer number to store O3 into, depending on
    ! whether this is a full chemistry run or a tagged O3 run.
    ! If tagged O3, tracer 2 should be the stratospheric tracer.  (dbj)
    !=================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
       NUM_TRACER = 1
    ELSE
       IF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN
          IF ( nAdvect > 1 ) THEN
             NUM_TRACER = 2
          ELSE
             NUM_TRACER = 1
          ENDIF
       ELSE
          ! All other simulations don't use O3...print error message
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, '(a)' ) 'This simulation does not use O3!!'
             WRITE( 6, '(a)' ) 'STOP in linoz_chem3!'
          ENDIF
          STOP
       ENDIF
    ENDIF

    ! Echo info
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 )
100    FORMAT( '     - LINOZ_CHEM3: Doing LINOZ' )
    ENDIF

    !=================================================================
    ! Perform stratospheric chemistry
    !=================================================================

    ! **** note dbj: check Spc(I,J,20:State_Grid%NZ,NTRACER) = with trop level
    ! ****         : check DMASS
    DO ITRC = 1, NUM_TRACER

       IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
          NTRACER = id_O3
       ELSE
          NTRACER = ITRC
       ENDIF

       ! Start at top layer and continue to lowest layer for strat. chem
       OUT_DATA = 0e+0_fp

       !IF ( preDebug ) THEN
       !   PRINT*, '### NUM_TRACER: ', NUM_TRACER
       !   PRINT*, '### ITRC: ', ITRC, 'JM: ', JM , 'IM', IM
       !   CALL FLUSH(6)
       !END IF

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I,       J,     LBOT,   LPOS,   L      ) &
       !$OMP PRIVATE( CLIMPML, DERO3, CLIMO3, DERCO3, DCO3   ) &
       !$OMP PRIVATE( DERTMP,  DTMP,  SSO3,   DMASS )
       DO J = 1, JM
       DO I = 1, IM

          !IF ( prtDebug ) THEN
          !   PRINT*, '### I: ', I, 'J: ', J
          !   CALL FLUSH(6)
          !END IF

          LBOT = State_Met%ChemGridLev(I,J) + 1
          LPOS = 1
          DO WHILE (State_Met%PEDGE(I,J,LPOS+1) .GE. 0.3e+0_f8)
             LPOS = LPOS +1
          ENDDO
          LPOS = LPOS-1

          !IF ( prtDebug ) CALL DEBUG_MSG('DONE GET_PEDGE')

#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
          !-----------------------------------------------------------
          !       %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
          !
          ! When we are connecting to the GEOS-5 GCM, we don't define
          ! L_OVERWRLD.  We cannot know the maximum extent of the
          ! tropopause in the GCM; we instead have to diagnose it at
          ! every timestep by comparing the pressure at a grid box to
          ! the tropopause pressure.  (bmy, 3/18/13)
          !-----------------------------------------------------------
#else
          !-----------------------------------------------------------
          !       %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
          !-----------------------------------------------------------

          ! dbj: for now, set tagged stratospheric tracer to total
          ! O3 in the overworld to avoid issues with spin ups
          IF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN
             Spc(I,J,(L_OVERWRLD+1):State_Grid%NZ,NTRACER) = &
                  Spc(I,J,(L_OVERWRLD+1):State_Grid%NZ,1)
          ENDIF
#endif
          !IF ( prtDebug ) CALL DEBUG_MSG( '### LINOZ_CHEM3: at LM, LBOT')
          !IF ( prtDebug ) CALL DEBUG_MSG('DONE TAGO3')

          ! Loop over levels
          DO L = LM, LBOT, -1

             !IF ( prtDebug ) THEN
             !   PRINT*, '### Spc: ', Spc(I,J,L,NTRACER)
             !   CALL FLUSH(6)
             !ENDIF

             ! Skip if tracer is negative
             ! SDE 2016-04-06: Changed from LE to LT
             IF ( Spc(I,J,L,NTRACER) .LT. 0.e+0_fp ) CYCLE

             ! calculate ozone column above box (and save)
             ! dcolo3 = ozone column (in DU) in given layer
             ! colo3 =  ozone column above layer + half of
             ! column in layer

             ! bdf Spc is in v/v, make conversion to DU
             if (l.eq.lm) then !top model layer
                dcolo3(i,j,l) = (Spc(i,j,l,NTRACER) *             &
                      State_Met%AD(I,J,L) / ( AIRMW               &
                      / State_Chm%SpcData(NTRACER)%Info%emMW_g )) &
                      / ( State_Grid%Area_M2(I,J) * 1e+4_fp )     &
                      * AVO / (AIRMW / ( AIRMW                    &
                      / State_Chm%SpcData(NTRACER)%Info%emMW_g )  &
                      *1e-3_fp) / 2.687e+16_fp
                colo3(i,j,l) = dcolo3(i,j,l)*0.5
             else
                dcolo3(i,j,l) = (Spc(i,j,l,NTRACER) *             &
                      State_Met%AD(I,J,L) / ( AIRMW               &
                      / State_Chm%SpcData(NTRACER)%Info%emMW_g )) &
                      / ( State_Grid%Area_M2(I,J) * 1e+4_fp )     &
                      * AVO / (AIRMW / ( AIRMW                    &
                      / State_Chm%SpcData(NTRACER)%Info%emMW_g )  &
                      *1e-3_fp) / 2.687e+16_fp
                colo3(i,j,l) = colo3(i,j,l+1) +                   &
                     (dcolo3(i,j,l)+dcolo3(i,j,l+1))*0.5
             endif
             out_data(i,j,l) = colo3(i,j,l)

             ! ++++++ climatological P-L:   ++++++
             climpml=tlstt(i,j,l,4)      ! Climatological P-L = (P-L)^o

             ! ++++++ local ozone feedback: ++++++
             dero3=tlstt(i,j,l,5)               ! Derivative w.r.t. O3. 
                                                !  dero3=-1/(time constant)
             IF (dero3.EQ.0) CYCLE              ! Skip Linoz if lifetime
                                                !  is infinite.
             climo3=tlstt(i,j,l,1)              ! Climatological O3 = f^o
             derco3=tlstt(i,j,l,7)              ! Derivative w.r.t. Column O3
             dco3=(colo3(i,j,l)-tlstt(i,j,l,3)) ! deviation from o3 climatology.
             ! ++++++ temperature feedback: ++++++
             dertmp=tlstt(i,j,l,6)              ! Derivative w.r.t. Temperature
             dtmp=(State_Met%T(I,J,L) - &       ! Deviation in Temperature
                      tlstt(i,j,l,2))           !  from climatology.

             ! ++++++ calculate steady-state ozone: ++++++
             sso3=climo3 - (climpml+dtmp*dertmp+dco3*derco3)/dero3

             ! ++++++ change in ozone mass due to chemistry: ++++++
             !ssO3 = f^*
             dmass=(sso3-Spc(I,J,L,NTRACER))*(1.0-exp(dero3*dtchem))

             ! ++++++ update ozone mass ++++++
             ! LINOX valid only up to 58 km, so do not use above 0.3 hPa
             ! dbj: impose exponential fall off of mixing ratio
             ! between 0.3 and 0.01 hPa (with fall off of a scale height)
             IF (State_Met%PEDGE(I,J,L) .LE. 0.3e+0_f8) THEN
                Spc(I,J,L,NTRACER) = &
                     (State_Met%PMID(I,J,L)/State_Met%PMID(I,J,LPOS-1))  &
                     * Spc(I,J,LPOS-1,NTRACER)
             ELSE
                Spc(I,J,L,NTRACER) = Spc(I,J,L,NTRACER)+DMASS
             ENDIF

          ENDDO       ! loop over L

       ENDDO          ! loop over I
       ENDDO          ! loop pver J
       !$OMP END PARALLEL DO

       !write our calculated column o3 maximum
       !write(6,*) 'max of columns= ',maxval(out_data)

    ENDDO
    IF ( prtDebug ) CALL DEBUG_MSG('DONE LINOZ_CHEM3')

    ! Free pointer
    Spc => NULL()

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
  SUBROUTINE LINOZ_STRATL( Input_Opt, State_Chm, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE PRESSURE_MOD
    USE State_Grid_Mod,     ONLY : GrdState
    USE TIME_MOD,           ONLY : GET_MONTH
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    REAL(fp)            :: P0L   (State_Grid%NZ+1        )
    REAL(fp)            :: STRT0L(State_Grid%NZ+1        )
    REAL(fp)            :: STRT1L(State_Grid%NZ+1        )
    REAL(fp)            :: STRT2L(State_Grid%NZ+1        )
    REAL(fp)            :: STRTX (Input_Opt%LINOZ_NLEVELS)
    REAL(fp)            :: YSTRT (Input_Opt%LINOZ_NLAT   )

    REAL(fp), POINTER   :: TLSTT (:,:,:,:) ! IM, JM, LM, Input_Opt%LINOZ_NFIELDS

    ! Month names
    CHARACTER(LEN=3)  :: CMONTH(12) = (/'jan', 'feb', 'mar', 'apr', &
                                        'may', 'jun', 'jul', 'aug', &
                                        'sep', 'oct', 'nov', 'dec'/)
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER :: PSF = 1010e+0_fp   ! Surface pressure [hPa]

    !=================================================================
    ! LINOZ_STRATL begins here!
    !=================================================================

    ! Initialization
    IM             = State_Grid%NX
    JM             = State_Grid%NY
    MONTH          = GET_MONTH()
    NFIELDS_LINOZ  = Input_Opt%LINOZ_NFIELDS
    NLAT_LINOZ     = Input_Opt%LINOZ_NLAT
    NLEVELS_LINOZ  = Input_Opt%LINOZ_NLEVELS
    NMONTHS_LINOZ  = Input_Opt%LINOZ_NMONTHS

    ! Point to the TLSTT array in state_chm_mod.F90
    TLSTT          => State_Chm%TLSTT

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(a)' ) REPEAT( '#', 79 )
       WRITE( 6,  50   ) CMONTH(MONTH)
       WRITE( 6, '(a)' ) REPEAT( '#', 79 )
50     FORMAT( '# Interpolating Linoz fields for ', a )
    ENDIF

    !--------------------------------------------------------------------------
    ! THE FOLLOWING LINEAR INTERPOLATION OPTION IS NOT USED, BUT CAN BE RESTORED
    ! BY UNCOMMENTING THE APPROPRIATE LINES OF CODE BELOW (bmy, 3/18/13)
    !
    !! get weights for month interpolation
    !do i=1,nmonths_linoz
    !  jdmc(i) = jdofm(i+1) - (jdofm(i+1)-jdofm(i))/2
    !enddo
    !
    !im1=0
    !do i=1,nmonths_linoz
    !  if (jdmc(i).lt.jday) then
    !    im1=i
    !  endif
    !enddo
    !if (im1.eq.0) then
    !  im1=nmonths_linoz
    !  im2=1
    !  wm1=(jdmc(im2)-jday)*1.0/(jdmc(im2)-(jdmc(im1)-365.0))
    !elseif (im1.eq.nmonths_linoz) then
    !  im2=1
    !  wm1=(jdmc(im2)+365.0-jday)/(jdmc(im2)+365.0-jdmc(im1))
    !else
    !  im2=im1+1
    !  wm1=(jdmc(im2)-jday)*1.0/(jdmc(im2)-jdmc(im1))
    !endif
    !wm2=1.0-wm1
    !
    !!write(6,*)iday,jday,' weights: ',wm1,wm2
    !!write(6,*)'months: ',im1,im2,month
    !!write(6,*)'between: ',jdmc(im1),jdmc(im2)
    !--------------------------------------------------------------------------

    ! YSTRT(J) are the latitudes (-85, -75, -65, ... +65, +75, +85)
    ! contained in the LINOZ climatology
    YSTRT(1) = -85.e+0_fp
    do J = 2, NLAT_LINOZ
       YSTRT(J) = YSTRT(J-1) + 10.e+0_fp
    enddo

    ! TLSTT is now 4D. Need to calculate latitude indeces inside of loop.
    ! (ckeller, 10/19/15)
    !! JLATMD(J) is the nearest-neighbor LINOZ data column corresponding
    !! to each GEOS-Chem latitude index J. (dbj, 6/25/03)
    !DO J = 1, State_Grid%NY
    !   JXXX      = INT( 0.1e+0_fp * State_Grid%YMid(1,J) + 10.e+0_fp )
    !   JLATMD(J) = MIN( 18,     MAX( 1, JXXX )          )
    !ENDDO

    ! P0L are the pressure at the level edges of each GEOS-Chem
    ! grid box, assuming a surface pressure of 1010 hPa (dbj, 6/25/03)
    DO L = 1, State_Grid%NZ+1
       P0L(L) = GET_AP(State_Grid%NZ+2-L) + &
              ( GET_BP(State_Grid%NZ+2-L) * PSF )
    ENDDO

    !=================================================================
    ! Lookup data in the LINOZ climatology
    !=================================================================

    ! Loop over the # of fields in the LINOZ climatology
    DO N = 1, NFIELDS_LINOZ

       !----------------------------------------------------------------------
       ! INTERPOLATION BETWEEN LATITUDES IS CURRENTLY NOT USED, BUT CAN BE
       ! RESTORED BY UNCOMMENTING THE APPROPRIATE LINES OF CODE BELOW
       ! (bmy, 3/18/13)
       !
       ! ***** Interpolation between latitudes is not currently used {PJC} ***
       !!----- interpolating along latitude, from TPAR2 to STRTXY
       !do K = 1, nlevels_linoz
       !do J = 1, nlat_linoz
       !   TPAR2(K,J) = TPARM(K,J,MONTH,N)
       !   TPAR2(K,J) = TPARM(K,J,im1,N)
       !enddo
       !enddo
       !call LINOZ_INTPL(nlevels_linoz,NLAT_LINOZ,JPAR,JM,YSTRT,YDGRD, &
       !                 TPAR2,STRTXY1)
       !do K = 1, nlevels_linoz
       !do J = 1, nlat_linoz
       !   TPAR2(K,J) = TPARM(K,J,im2,N)
       !enddo
       !enddo
       !call LINOZ_INTPL(nlevels_linoz,NLAT_LINOZ,JPAR,JM,YSTRT,YDGRD, &
       !     TPAR2,STRTXY2)
       !----------------------------------------------------------------------

       ! Loop over GEOS-Chem latitudes
       DO J = 1, JM
       DO I = 1, IM

          ! JJ is the index of the nearest LINOZ data column
          ! corresponding to GEOS-Chem latitude index J
          ! Now explicitly calculate for every grid box to account for
          ! curvilinear grids. (ckeller, 10/19/15)
          ! JJ = JLATMD(J)
          JXXX = INT( 0.1e+0_fp * State_Grid%YMid(I,J) + 10.e+0_fp )
          JJ   = MIN( 18, MAX( 1, JXXX ) )

          ! Loop over the # of levels in the LINOZ climatology
          DO K = 1, NLEVELS_LINOZ

             ! STRTX(K) is the column of data from the LINOZ climatology
             ! for the given month, latitude (JJ), level (K), and species
             ! (N).  One of the following interpolation options may be used.

             !-----------------------------------------------------------------
             ! THE FOLLOWING OPTIONS ARE CURRENTLY NOT USED, BUT CAN BE RESTORED
             ! BY UNCOMMENTING THE APPROPRIATE LINES OF CODE BELOW (bmy,3/18/13)
             !
             !! linearly interpolate in latitude and month
             !STRTX(K) = STRTXY1(K,J)*wm1 + STRTXY2(K,J)*wm2
             !
             !! linearly interpolate in latitude, single month
             !STRTX(K) = STRTXY2(K,J)
             !
             !! nearest latitude, linearly interpolate in month
             !!STRTX(K) = Input_Opt%LINOZ_TPARM(K,JJ,im1,N)*wm1 &
             !          + Input_Opt%LINOZ_TPARM(K,JJ,im2,N)*wm2
             !-----------------------------------------------------------------

             ! Nearest neighbor, no interpolation
             STRTX(K) = Input_Opt%LINOZ_TPARM(K,JJ,MONTH,N)
          ENDDO

          ! *PJC* Interpolate and calculate moments of column distribution
          CALL LINOZ_STRT2M( Input_Opt,     State_Grid,        &
                             State_Grid%NZ, STRTX,     P0L,    &
                             STRT0L,        STRT1L,    STRT2L, RC )

          ! Store loss freq/yields & moments in TLSTT/SWT/SWW
          ! for exact CTM layers LM down
          ! Order reversed from C.McLinden version {PJC}
          DO LR = 1,State_Grid%NZ
             TLSTT(I,J,LR,N) = STRT0L(State_Grid%NZ+1-LR)
          ENDDO

       ENDDO   ! loop over I
       ENDDO   ! loop over J

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
  SUBROUTINE LINOZ_STRT2M( Input_Opt, State_Grid,         &
                           NSTRT,     STRTX,      P0L,    &
                           STRT0L,    STRT1L,     STRT2L, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Grid_Mod,     ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt    ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid   ! Grid State object
    INTEGER,        INTENT(IN)  :: NSTRT        ! # of levels( = State_Grid%NZ)
    REAL(fp),       INTENT(IN)  :: P0L(State_Grid%NZ+1) ! Pressure edges

    ! Quantity on the LINOZ vertical grid (i.e. fields #1-7 of the LINOZ clim.
    REAL(fp),       INTENT(IN)  :: STRTX(Input_Opt%LINOZ_NLEVELS)
!
! !OUTPUT PARAMETERS:
!
    ! 0th moment of distribution, on GEOS-Chem grid edges
    REAL(fp),       INTENT(OUT) :: STRT0L(State_Grid%NZ+1)

    ! 1st moment of distribution, on GEOS-Chem grid edges
    REAL(fp),       INTENT(OUT) :: STRT1L(State_Grid%NZ+1)

    ! 2nd moment of distribution, on GEOS-Chem grid edges
    REAL(fp),       INTENT(OUT) :: STRT2L(State_Grid%NZ+1)

    ! Success or failure?
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER  :: NL,  NCBOX, L,      K,    NX
    REAL(fp) :: P1,  P2,    F0,     F1
    REAL(fp) :: F2,  XPSD,  XPSLM1, XPSL

    ! Arrays
    REAL(fp) :: PS(Input_Opt%LINOZ_NLEVELS+6)  ! Old code: PS(NL+1)
    REAL(fp) :: F (Input_Opt%LINOZ_NLEVELS+5)  ! Old code: PS(NL)

    !=================================================================
    ! Initialization
    !=================================================================

    ! Assume success
    RC         = GC_SUCCESS

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

    XPSD       = 10.e+0_fp **(-0.125e+0_fp)
    XPSLM1     = 1000.e+0_fp
    PS(1)      = 1000.e+0_fp
    DO L = 2,NL
       XPSL    = XPSLM1 *XPSD
       PS(L)   = 0.5e+0_fp *(XPSLM1 +XPSL)
       XPSLM1  = XPSL
    ENDDO
    PS(NL+1)   = 0.e+0_fp
    DO L = 1,NL-NX
       F(L)     = 0.e+0_fp
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
    INTEGER,   INTENT(IN)  :: NL
    REAL(fp),  INTENT(IN)  :: F(NL)
    REAL(fp),  INTENT(IN)  :: PS(NL+1)
    REAL(fp),  INTENT(IN)  :: P1
    REAL(fp),  INTENT(IN)  :: P2
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),  INTENT(OUT) :: F0
    REAL(fp),  INTENT(OUT) :: F1
    REAL(fp),  INTENT(OUT) :: F2
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
    real(fp)   XB,XC,PC,PB,THIRD,sgnf0

    F0     = 0.e+0_fp
    F1     = 0.e+0_fp
    F2     = 0.e+0_fp
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
          F1 = F1 +F(I) *3.e+0_fp *((XC *XC -XC) - (XB *XB -XB))
          F2 = F2 +F(I) *5.e+0_fp * &
               ((XC+XC-1.e+0_fp)*(XC*XC -XC) - &
                (XB+XB-1.e+0_fp)*(XB*XB -XB))
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
    THIRD = 1.e+0_fp/3.e+0_fp
    IF (F2 .GT. 0.e+0_fp)  THEN

       ! do not allow reversal of curvature: F2 > 0
       F2   = MIN(F2, ABS(F1)*THIRD, 5.e-1_fp*F0)
       IF (F1 .LT. 0.e+0_fp)  THEN
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
  SUBROUTINE LINOZ_READ( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE InquireMod,    ONLY : findFreeLun
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
!  LINOZ_READ is called from "main.F90" at the start of the simulation.
!  LINOZ_READ will also call INIT_LINOZ to initialize the arrays.
!
! !REVISION HISTORY:
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: FileExists
    INTEGER            :: IU_FILE, K,     J
    INTEGER            :: M,       NTBLS, IOS
    REAL(fp)           :: TMAX,    TMIN

    ! Dimension sizes
    INTEGER            :: NFIELDS_LINOZ
    INTEGER            :: NLAT_LINOZ
    INTEGER            :: NLEVELS_LINOZ
    INTEGER            :: NMONTHS_LINOZ

    ! Strings
    CHARACTER(LEN=80)  :: HEADING, TITL1
    CHARACTER(LEN=255) :: ERRMSG
    CHARACTER(LEN=255) :: FILENAME
    CHARACTER(LEN=255) :: FileMsg
    CHARACTER(LEN=255) :: ThisLoc

    !=================================================================
    ! LINOZ_READ begins here!
    !=================================================================

    ! Assume success
    RC             = GC_SUCCESS
    ErrMsg         = ''
    ThisLoc        = ' -> at LINOZ_READ (in GeosCore/linoz_mod.F90)'

    ! Get fields from Input_Opt
    NFIELDS_LINOZ  = Input_Opt%LINOZ_NFIELDS
    NLAT_LINOZ     = Input_Opt%LINOZ_NLAT
    NLEVELS_LINOZ  = Input_Opt%LINOZ_NLEVELS
    NMONTHS_LINOZ  = Input_Opt%LINOZ_NMONTHS

    ! Filename
    FILENAME       = TRIM( Input_Opt%CHEM_INPUTS_DIR ) // &
                     'Linoz_200910/Linoz_March2007.dat'

    !=================================================================
    ! In dry-run mode, print file path to dryrun log and exit.
    ! Otherwise, print file path to stdout and continue.
    !=================================================================

    ! Test if the file exists
    INQUIRE( FILE=TRIM( FileName ), EXIST=FileExists )

    ! Test if the file exists and define an output string
    IF ( FileExists ) THEN
       FileMsg = 'LINOZ (LINOZ_READ): Opening'
    ELSE
       FileMsg = 'LINOZ (LINOZ_READ): REQUIRED FILE NOT FOUND'
    ENDIF

    ! Write message to stdout for both regular and dry-run simulations
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( FileName )
300    FORMAT( a, ' ', a )
    ENDIF

    ! For dry-run simulations, return to calling program.
    ! For regular simulations, throw an error if we can't find the file.
    IF ( Input_Opt%DryRun ) THEN
       RETURN
    ELSE
       IF ( .not. FileExists ) THEN
          WRITE( ErrMsg, 300 ) TRIM( FileMsg ), TRIM( FileName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! Read climatological data from the LINOZ lookup tables
    !=================================================================

    ! Find a free file LUN
    IU_FILE = findFreeLUN()

    ! Define error message
    ERRMSG  = 'LINOZ_READ (in GeosCore/linoz_mod.F90)'

    ! new std z*=2km levels from model:  z*=10,12,...(25*2)+8 km
    OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', &
          FORM='FORMATTED',      IOSTAT=IOS )

    ! Return if there was an error opening the file
    IF ( IOS /= 0 ) THEN
       WRITE( FileMsg, 300 ) TRIM( FileMsg ), TRIM( FileName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Read header
    READ ( IU_FILE, '(a)' ) HEADING
    IF ( Input_Opt%amIRoot ) THEN
       WRITE(6,*) TRIM( HEADING )
    ENDIF

    ! Loop over # of fields in the LINOZ climatology
    DO NTBLS = 1, NFIELDS_LINOZ

       ! Zero min & max values
       TMIN = +1.e+30_fp
       TMAX = -1.e+30_fp

       ! Skip header line
       READ( IU_FILE, '(a)' ) TITL1

       ! Loop over # of months in the LINOZ climatology
       do M = 1, NMONTHS_LINOZ

          ! Loop over # of latitudes in the LINOZ climatology
          do J = 1, NLAT_LINOZ

             ! Read data into Input_Opt%LINOZ_TPARM
             READ( IU_FILE, '(20X,6E11.4/(8E11.4))', IOSTAT=IOS ) &
                   ( Input_Opt%LINOZ_TPARM(K,J,M,NTBLS), &
                     K=NLEVELS_LINOZ,1,-1                 )

             ! Stop on error
             IF ( IOS > 0 ) THEN
                ErrMsg = 'Error reading ' // TRIM( FileName )
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
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
       IF ( Input_Opt%amIRoot ) write (6,912) TITL1,TMIN,TMAX
912    FORMAT('  Linoz Data:  ',a80,1p,2e10.3)

    enddo

    ! Echo info
    IF ( Input_Opt%amIRoot ) THEN
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
    INTEGER,  INTENT(IN)  :: KE
    INTEGER,  INTENT(IN)  :: IE
    INTEGER,  INTENT(IN)  :: ND
    INTEGER,  INTENT(IN)  :: NE
    REAL(fp), INTENT(IN)  :: XI(IE)
    REAL(fp), INTENT(IN)  :: XN(ND)
    REAL(fp), INTENT(IN)  :: YI(KE,IE)
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: YN(KE,ND)
!
! !REVISION HISTORY:
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer  I,II,J,K
    real(fp)   CNST1,CNST2

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
END MODULE LINOZ_MOD
