!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: modis_lai_mod.F90
!
! !DESCRIPTION: Module MODIS\_LAI\_MOD reads the MODIS LAI and CHLR data at 
!  native resolution and then regrids them to the GEOS-Chem resolution on the
!  fly. CHLR data is only read if marine organic aerosol tracers are enabled. 
!
! !INTERFACE: 
!
MODULE Modis_Lai_Mod
!
! !USES:
!
  USE CMN_SIZE_Mod                                ! Size parameters
  USE Error_Mod                                   ! Error checking routines
  USE PRECISION_MOD                               ! For GEOS-Chem Precision (fp)
  USE Mapping_Mod                                 ! Mapping weights & areas
  USE Time_Mod                                    ! EXPAND_DATE

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Init_Modis_Lai
  PUBLIC  :: Read_Modis_Lai
  PUBLIC  :: Compute_Modis_Lai
  PUBLIC  :: Find_Lai_Month
  PUBLIC  :: Cleanup_Modis_Lai
!
! !PRIVATE MEMBER FUNCTIONS:
!
#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
  PUBLIC  :: Compute_XLAI_GCHP
#endif
  PUBLIC  :: Read_Modis
  PUBLIC  :: Compute_Modis
  PRIVATE :: RoundOff
!
! !REMARKS:
!  (1) MODIS LAI data resolution is the same as the Olson land map. The Olson
!      2001 landmap is the default, therefore, you will use MODIS LAI data at
!      0.25 x 0.25 resolution. 
!      Note that MODIS CHLR data is only available at 0.25 x 0.25 resolution.
!  (2) In HEMCO, MEGAN uses 'offline' MODIS LAI (State_Met%MODISLAI) which
!      is computed in this module.
!      in an ESMF environment, in which case State_Met%LAI is used instead. 
!  (3) MODIS LAI arrays and where they are used in GEOS-Chem:
!       (a) State_Met%XLAI     --> dry deposition routine DEPVEL
!       (b) State_Met%MODISLAI --> MEGAN (if using standard GC); Hg(0)
!                                  emissions; several diagnostics
!
! !REVISION HISTORY:
!  03 Apr 2012 - R. Yantosca - Initial version
!  05 Apr 2012 - R. Yantosca - Added descriptive comments
!  09 Apr 2012 - R. Yantosca - Fixed error in ROUNDOFF function that caused
!                              numbers to be rounded up incorrectly.
!  09 Apr 2012 - R. Yantosca - Changed variables to REAL(fp)
!  09 Apr 2012 - R. Yantosca - Now set MODIS_START and MODIS_END depending
!                              on which version of MODIS LAI we are using
!  13 Dec 2012 - R. Yantosca - Remove reference to obsolete CMN_DEP_mod.F;
!                              XLAI, XLAI2 now are carried in State_Met
!  23 Jun 2014 - R. Yantosca - Removed references to logical_mod.F
!  09 Oct 2014 - C. Keller   - Removed GC_LAI_PM, GC_LAI_CM, GC_LAI_NM and
!                              MODIS_LAI_PM.
!  17 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  07 Jul 2015 - E. Lundgren - Now also read and compute MODIS chlorophyll-a 
!                              (B. Gantt, M. Johnson). Use separate end years.
!  18 Oct 2016 - E. Lundgren - Move module vars GC_LAI and GC_CHLR to State_Met;
!                              rename as MODISLAI and MODISCHLR; make
!                              State_Met vars XLAI2 and XCHLR2 module vars
!  22 Jan 2019 - H.P. Lin    - Move XLAI2 and XCHLR2 to State_Met for compatibility
!                              with multi-domain WRF-GC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Scalars
  INTEGER :: I_MODIS             ! # of longitudes, MODIS grid
  INTEGER :: J_MODIS             ! # of latitudes,  MODIS grid
  INTEGER :: MODIS_START         ! First year of MODIS data (both LAI and CHLR) 
  INTEGER :: LAI_END             ! Last  year of MODIS LAI data
  INTEGER :: CHLR_END            ! Last  year of MODIS CHLR data
                                              
  ! MODIS grid arrays                                    
  REAL*4,  ALLOCATABLE, TARGET :: MODIS_LAI   (:,:) ! Daily LAI on MODIS grid
  REAL*4,  ALLOCATABLE, TARGET :: MODIS_LAI_CM(:,:) ! MODIS LAI for current mo
  REAL*4,  ALLOCATABLE, TARGET :: MODIS_LAI_NM(:,:) ! MODIS LAI for next month 
  REAL*4,  ALLOCATABLE, TARGET :: MODIS_CHLR(:,:)   ! Daily CHLR on MODIS grid
  REAL*4,  ALLOCATABLE, TARGET :: MODIS_CHLR_CM(:,:)! MODIS CHLR for current mo
  REAL*4,  ALLOCATABLE, TARGET :: MODIS_CHLR_NM(:,:)! MODIS CHLR for next month

  ! specify midmonth day for year 2000
  INTEGER, PARAMETER  :: startDay(13) = (/  15,  45,  74, 105,      &
                                            135, 166, 196, 227,     &
                                            258, 288, 319, 349, 380/)

  ! specify number of digits of precision for LAI and CHLR
  INTEGER, PARAMETER  :: numRoundLAI  = 1 ! Data precision for LAI
  INTEGER, PARAMETER  :: numRoundCHLR = 3 ! Data precision for CHLR

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: compute_modis_lai
!
! !DESCRIPTION: Subroutine COMPUTE\_MODIS\_LAI is the wrapper routine to
!  compute the daily MODIS leaf area indices for GEOS-Chem directly from the 
!  native grid resolution  (0.25 x 0.25 or 0.5 x 0.5). If marine organic
!  aerosol tracers are used, then daily MODIS chlorophyll is also computed.
!  This routine is not used in High Performance GEOS-Chem (GCHP).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Modis_Lai( am_I_Root,    Input_Opt, State_Met,  &
                                doy,          mm,        mapping,    &
                                wasModisRead, RC                    )
!
! !USES:
!
    USE State_Met_Mod, ONLY : MetState
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput 
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN)  :: am_I_Root     ! Are we on the root CPU?
    TYPE(OptInput),  INTENT(IN)  :: Input_Opt     ! Input Options object
    INTEGER,         INTENT(IN)  :: doy           ! Day of year
    INTEGER,         INTENT(IN)  :: mm            ! Month for LAI data
    TYPE(MapWeight), POINTER     :: mapping(:,:)  ! "fine" -> "coarse" grid map
    LOGICAL,         INTENT(IN)  :: wasModisRead  ! Was LAI data just read in?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState),  INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT) :: RC            ! Success or failure?
!
! !REVISION HISTORY: 
!  03 Apr 2012 - R. Yantosca - Initial version
!  05 Apr 2012 - R. Yantosca - Renamed arg "doMonthly" to "wasModisRead"
!  09 Apr 2012 - R. Yantosca - Changed variables to REAL(fp)
!  09 Apr 2012 - R. Yantosca - Now follows same algorithm as rdlai.F for
!                              populating XLAI array
!  09 Apr 2012 - R. Yantosca - Remove refs to CMN_VEL_mod.F and XYLAI array;
!                              these are now obsolete
!  17 Apr 2012 - R. Yantosca - Now rename "map" object to "mapping" to avoid
!                              name confusion w/ an F90 intrinsic function
!  13 Dec 2012 - R. Yantosca - Add am_I_Root, State_Met, RC arguments
!  13 Dec 2012 - R. Yantosca - XLAI, XLAI2 are now carried in State_Met
!                              instead of in obsolete Headers/CMN_DEP_mod.F
!  23 Jun 2014 - R. Yantosca - Now accept Input_Opt via the arg list
!  09 Oct 2014 - C. Keller   - Removed GC_LAI_PM, GC_LAI_CM, GC_LAI_NM and
!                              MODIS_LAI_PM.
!  08 Jul 2015 - E. Lundgren - Now compute LAI and CHLR data. Abstracted 
!                              compute code to new routine Compute_MODIS 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL             :: ComputeLAI       ! T = compute LAI, F = compute CHLR

    ! Assume success
    RC = GC_SUCCESS

    ! Always compute LAI
    ComputeLAI = .true.
    CALL Compute_Modis( am_I_Root,  ComputeLAI, Input_Opt,    State_Met,  &
                        doy,  mm,   mapping,    wasModisRead, RC )

    ! Only compute CHLR if organic marine aerosols are tracers
    IF ( Input_Opt%LMPOA ) THEN

       ComputeLAI = .false.
       CALL Compute_Modis( am_I_Root,  ComputeLAI, Input_Opt,    State_Met, &
                           doy,  mm,   mapping,    wasModisRead, RC )
    ENDIF

  END SUBROUTINE Compute_Modis_LAI
!EOC
#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: compute_xlai_gchp
!
! !DESCRIPTION: Subroutine COMPUTE\_XLAI\_GCHP computes MODIS-based leaf
!  area indices (LAI) per land type and grid cell. This computation uses 
!  offline 0.25x0.25 MODIS LAI/CHLR and Olson landmap data regridded to 
!  the cubed sphere. Variables set include State\_Met%XLAI.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_XLAI_GCHP( am_I_Root, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN)  :: am_I_Root     ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState),  INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT) :: RC            ! Success or failure?
!
! !REMARKS:
!
!
! !REVISION HISTORY: 
!  18 Oct 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, S, T
    CHARACTER(len=2)   :: landStr
    REAL(fp)           :: landFrac
    REAL(fp), POINTER  :: XModisPtr(:,:,:) 
    REAL(fp), POINTER  :: XPtr(:,:,:)  

    !======================================================================
    ! Initialize
    !======================================================================

    ! Assume success
    RC                = GC_SUCCESS

    ! Initialize LAI variables
    State_Met%XLAI = 0.0_fp
    State_Met%MODISLAI = 0.0_fp

    ! Loop over all grid cells
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Loop over all surface types present in this grid cell
       DO S = 1, State_Met%IREG(I,J)
       
          ! Set current surface type index
          T = State_Met%ILAND(I,J,S) + 1

          ! Get fraction of cell with this surface type by retrieving it
          ! from the land type fraction calculated by ExtData
          landFrac = State_Met%LandTypeFrac(I,J,T)

          ! Set XLAI to average LAI for this surface type ( as calculated 
          ! by ExtData using zeros for coverage by other surface types ) 
          ! divided by the fractional coverage of this surface type.
          ! The resultant XLAI is the average LAI for only the area 
          ! with the current surface type, and therefore is larger than 
          ! XLAI_NATIVE when other surface types exist within the cell.
          ! 
          ! NOTE: Unlike XLAI_NATIVE and LandTypeFrac, the 3rd 
          ! dimension indexes of XLAI are NOT surface types 1-73! Instead,
          ! It is the surface indexes ILAND-1 and therefore contains 
          ! zeros beyond the number of surface types present in the cell 
          ! (IREG). This is for backwards compatibility with GC classic
          ! legacy drydep code.
          IF ( landFrac .gt. 1.e-9_fp ) THEN
             State_Met%XLAI(I,J,S) = State_Met%XLAI_NATIVE(I,J,T) / landFrac
          ENDIF
       
       ENDDO

       ! Calculate average LAI for this grid cell across all land types
       State_Met%MODISLAI(I,J) = SUM( State_Met%XLAI_NATIVE(I,J,:) )

    ENDDO
    ENDDO

  END SUBROUTINE Compute_XLAI_GCHP
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: compute_modis
!
! !DESCRIPTION: Subroutine COMPUTE\_MODIS computes either the daily MODIS 
!  leaf area indices or the daily chlorophyll for GEOS-Chem directly from 
!  the native grid resolution (0.25 x 0.25 or 0.5 x 0.5).  The XLAI array 
!  (used in the legacy soil NOx and dry deposition routines) are populated 
!  accordingly. The XYLAI array is now obsolete and has been replaced by XLAI.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Modis( am_I_Root,  ComputeLAI,                      &
                            Input_Opt,  State_Met,    doy,               & 
                            mm,         mapping,      wasModisRead,  RC )

!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN)  :: am_I_Root     ! Are we on the root CPU?
    LOGICAL,         INTENT(IN)  :: ComputeLAI 
    TYPE(OptInput),  INTENT(IN)  :: Input_Opt     ! Input Options object
    INTEGER,         INTENT(IN)  :: doy           ! Day of year
    INTEGER,         INTENT(IN)  :: mm            ! Month for LAI data
    TYPE(MapWeight), POINTER     :: mapping(:,:)  ! "fine" -> "coarse" grid map
    LOGICAL,         INTENT(IN)  :: wasModisRead  ! Was LAI data just read in?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState),  INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT) :: RC            ! Success or failure?
!
! !REMARKS:
!  Uses same algorithm as RDISOLAI in the existing lai_mod.F.
!
! !REVISION HISTORY: 
!  07 Jul 2015 - E. Lundgren - Initial version, contains old Compute_Modis_Lai
!                              code plus modifications to compute CHLR
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE :: FIRST = .TRUE.

    ! Scalars
    INTEGER   :: I,      J,     IMUL,     ITD
    INTEGER   :: C,      II,    JJ,       type,  K,      numRound
    REAL(fp)  :: mapWt,  area,  sumArea,  DMON,  DITD,   DIMUL

    ! Arrays
    REAL(fp)  :: tempArea   (0:NSURFTYPE-1)
    REAL(fp)  :: tempModis  (0:NSURFTYPE-1)
    REAL(fp)  :: tempModisCm(0:NSURFTYPE-1)
    REAL(fp)  :: tempModisNm(0:NSURFTYPE-1)

    ! Pointers
    REAL(fp), POINTER  :: GC_PTR(:,:)
    REAL*4,   POINTER  :: MODIS_PTR(:,:)
    REAL*4,   POINTER  :: MODIS_PTR_CM(:,:)
    REAL*4,   POINTER  :: MODIS_PTR_NM(:,:)
    REAL(fp), POINTER  :: XTMP(:,:,:)
    REAL(fp), POINTER  :: XTMP2(:,:,:)

    !======================================================================
    ! Interpolate the data on the MODIS grid to current day
    ! Use same algorithm as in routines RDISOLAI (in lai_mod.F)
    !======================================================================
    
    ! Assume success
    RC                = GC_SUCCESS

    ! Assign pointers and precision based on whether computing LAI or CHLR
    IF ( ComputeLAI ) THEN
       GC_PTR            => State_Met%MODISLAI  ! Daily computed LAI per
                                                ! GC grid cell
       MODIS_PTR         => MODIS_LAI
       MODIS_PTR_CM      => MODIS_LAI_CM
       MODIS_PTR_NM      => MODIS_LAI_NM
       XTMP              => State_Met%XLAI  ! LAI/land type, GCgrid, this month
       XTMP2             => State_Met%XLAI2 ! LAI/land type, GCgrid, next month
       numRound          =  numRoundLAI
    ELSE
       GC_PTR            => State_Met%MODISCHLR ! Daily computed CHLR per
                                                ! GC grid cell
       MODIS_PTR         => MODIS_CHLR 
       MODIS_PTR_CM      => MODIS_CHLR_CM
       MODIS_PTR_NM      => MODIS_CHLR_NM
       XTMP              => State_Met%XCHLR ! CHLR/land type, GCgrid, this month
       XTMP2             => State_Met%XCHLR2! CHLR/land type, GCgrid, next month
       numRound          =  numRoundCHLR 
    ENDIF

    ! IMUL is days since midmonth
    ! ITD  is days between midmonths
    IF ( doy < startDay(1) ) THEN
       IMUL           = 365 + doy - startDay(12) 
       ITD            = 31
    ELSE
       IMUL           = doy            - startDay(mm)
       ITD            = startDay(mm+1) - startDay(mm)
    ENDIF

!    ! Archive the days between midmonths in the LAI data
!    DAYS_BTW_MON      = ITD

    ! Cast ITD, IMUL to REAL(fp)
    DITD              = DBLE( ITD  )
    DIMUL             = DBLE( IMUL )

    ! Fraction of the LAI month that we are in
    DMON              = REAL( IMUL ) / REAL( ITD ) 
       
    ! Interpolate to daily values on the MODIS grid
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) & 
    !$OMP PRIVATE( I, J   )
    DO J = 1, J_MODIS
    DO I = 1, I_MODIS
       MODIS_PTR(I,J) = MODIS_PTR_CM(I,J)  &
                      + ( ( MODIS_PTR_NM(I,J) - MODIS_PTR_CM(I,J) ) * DMON )
    ENDDO
    ENDDO 
    !$OMP END PARALLEL DO

    !======================================================================
    ! Bin data from the "fine" MODIS grid to the "coarse" GEOS-Chem grid.
    ! Populate arrays for backwards-compatibility w/ existing routines
    ! Use same algorithm as in routine rdlai.F
    !======================================================================
    !$OMP PARALLEL DO                                                 &
    !$OMP DEFAULT( SHARED                                           ) &
    !$OMP PRIVATE( I,           J,           tempArea, tempModis    ) &
    !$OMP PRIVATE( tempModisCm, tempModisNm, sumArea                ) &
    !$OMP PRIVATE( C,           II,          JJ,       type         ) & 
    !$OMP PRIVATE( area,        K                                   )      
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Initialize
       tempArea             = 0e+0_fp
       tempModis            = 0e+0_fp
       tempModisCm          = 0e+0_fp
       tempModisNm          = 0e+0_fp
       sumArea              = mapping(I,J)%sumarea
       GC_PTR(I,J)          = 0e+0_fp

       !-------------------------------------------------------------------
       ! Sum up the leaf area indices or chlorophyl-a from all of the the 
       ! "fine" grid boxes (II,JJ) that are located within "coarse" grid 
       ! box (I,J)
       !-------------------------------------------------------------------
       DO C = 1, mapping(I,J)%count

          ! Extract fields from MAP object
          II                = mapping(I,J)%II(C)
          JJ                = mapping(I,J)%JJ(C)
          type              = mapping(I,J)%olson(C)
          area              = mapping(I,J)%area(C)

          ! Sum of areas corresponding to each Olson
          ! for "coarse" GEOS-Chem grid box (I,J)
          tempArea(type)    = tempArea(type)  + area 

          ! Compute the total leaf area or chlorophyll-a in "coarse" 
          ! GEOS-Chem grid box (I,J) corresponding to each Olson land type
          tempModisCm(type) = tempModisCm(type) +                &
                                     ( MODIS_PTR_CM(II,JJ) * area )
          tempModisNm(type) = tempModisNm(type) +                &
                                     ( MODIS_PTR_NM(II,JJ) * area )

          ! Compute the total leaf area or chlorophyll-a in "coarse" 
          ! GEOS-Chem grid box (I,J), irrespective of Olson land type
          GC_PTR(I,J) = GC_PTR(I,J) + ( MODIS_PTR(II,JJ) * area )

       ENDDO

       !-------------------------------------------------------------------
       ! Compute the resultant (i.e. for all land types) daily-interpolated 
       ! values for the "coarse" GEOS-Chem grid box (I,J).  DAILY_LAI is
       ! a replacement for the ISOLAI array from "lai_mod.F".
       !-------------------------------------------------------------------

       ! Convert leaf area [cm2 leaf] to LAI [cm2 leaf/cm2 grid box],
       ! or chlorophyll-a [mg/m] to [mg/m/cm2] for grid box (I,J), 
       ! irrespective of Olson land type. 
       GC_PTR(I,J) = GC_PTR(I,J) / sumArea

       !-------------------------------------------------------------------
       ! Compute the LAI or CHLR for each Olson land type at GEOS-Chem grid 
       ! box (I,J).  These will be used to populate the XLAI & XYLAI arrays.
       !-------------------------------------------------------------------
       DO C = 0, NSURFTYPE-1
          
          ! Skip land types that are not in "coarse" grid box (I,J)
          IF ( tempArea(C) > 0e+0_fp ) THEN
          
             ! Ordering for ILAND, IUSE, XLAI, XYLAI etc arrays
             K = mapping(I,J)%ordOlson(C)

             ! Convert leaf area [cm2 leaf] to LAI [cm2 leaf/cm2 grid box]
             tempModisCm(C) = tempModisCm(C) / tempArea(C)
             tempModisNm(C) = tempModisNm(C) / tempArea(C)

             ! Round off to digits of precision to mimic the fact that
             ! the LAI in the lai*.global files only had one decimal point,
             ! Round CHLR to three digits of precision.
             tempModisCm(C) = RoundOff( tempModisCm(C), numRound )
             tempModisNm(C) = RoundOff( tempModisNm(C), numRound )

             ! This IF statement mimics the algorithm in the obsolete
             ! routine rdlai.F.  We need to keep the same algorithm
             ! for backwards compatibility.
             IF ( FIRST ) THEN 

                !----------------------------------------------------------
                ! %%%%% START OF SIMULATION, FIRST DATA READ %%%%%
                !
                ! Follow original algorithm in the old rdland.F. Example for
                ! LAI is:
                ! (1) XLAI  gets read in as the current month of LAI data
                ! (2) XLAI2 gets read in as the next    month of LAI data
                ! (3) XLAI2 is recomputed as the Delta-LAI this month
                !     ([ next month - this month ] / # of days in month)
                ! (4) XLAI is incremented by the amount 
                !      ( Delta-LAI ) * # of days since start of month
                ! 
                ! NOTE: As of Dec 2012, XLAI and XLAI2 are now part of 
                ! the Meteorology State object (bmy, 12/13/12)
                !----------------------------------------------------------
                XTMP2(I,J,K) = ( tempModisNm(C) - tempModisCm(C) ) / DITD       
                XTMP(I,J,K) = tempModisCm(C) + ( XTMP2(I,J,K) * DIMUL ) 

             ELSE

                IF ( wasModisRead ) THEN

                   !----------------------------------------------------------
                   ! %%%% SUBSEQUENT DATA READ @ START OF NEW LAI/CHLR MONTH %
                   !
                   ! Follow original algorithm in the old rdland.F
                   ! (1) XLAI (or XCHLR) gets read in as current month of data
                   ! (2) XLAI2 (or XCHLR2) is computed as the Delta of this 
                   !     month (i.e. [next month - this month ] / # of days
                   !----------------------------------------------------------
                   XTMP2(I,J,K) = ( tempModisNm(C) - tempModisCm(C) ) / DITD
                   XTMP(I,J,K) = tempModisCm(C)       
 
                ELSE
                
                   !----------------------------------------------------------
                   ! %%%%% ALL OTHER TIMES OF THE MONTH (NO DATA READS %%%%%
                   !
                   ! Follow original algorithm in the old rdland.F
                   ! (1) Increment LAI or CHLR by this month's Delta-LAI
                   !----------------------------------------------------------
                   XTMP(I,J,K) = XTMP(I,J,K) + XTMP2(I,J,K)

                ENDIF
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Save
    FIRST = .FALSE.

    ! Cleanup pointers
    GC_PTR       => NULL()
    MODIS_PTR    => NULL()
    MODIS_PTR_CM => NULL()
    MODIS_PTR_NM => NULL()
    XTMP         => NULL()
    XTMP2        => NULL()

  END SUBROUTINE Compute_Modis
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_modis_lai
!
! !DESCRIPTION: Subroutine READ\_MODIS\_LAI is the wrapper routine to read 
!  the MODIS LAI from disk in netCDF format for the current month, and for 
!  next month. If enabled, MODIS CHLR is also read in the same way as LAI.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Read_Modis_Lai( am_I_Root, Input_Opt, yyyy, mm, wasModisRead, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput 
   
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: am_I_Root     ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)  :: Input_Opt     ! Input Options object
    INTEGER,        INTENT(IN)  :: yyyy          ! Year for LAI data
    INTEGER,        INTENT(IN)  :: mm            ! Month for LAI data
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,        INTENT(OUT) :: wasModisRead  ! Was LAI data just read in?
    INTEGER,        INTENT(OUT) :: RC            ! Success or failure?   
!
! !REVISION HISTORY:
!  03 Apr 2012 - R. Yantosca - Initial version
!  05 Apr 2012 - R. Yantosca - Renamed arg "doMonthly" to "wasModisRead"
!  05 Jun 2013 - R. Yantosca - Bug fix, use "mm" for current month index
!  20 Jun 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!  09 Oct 2014 - C. Keller   - MODIS_LAI_PM not needed anymore.
!  05 Mar 2015 - R. Yantosca - Now read data w/r/t ExtData/CHEM_INPUTS
!  07 Jul 2015 - E. Lundgren - Generalized to read either LAI or CHLR
!  08 Jul 2015 - E. Lundgren - Now read LAI and CHLR data. Abstracted 
!                              read code to new routine Read_Modis.
!  13 Sep 2017 - M. Sulprizio- Remove Input_Opt%USE_OLSON_2001. Olson 2001 is
!                              now the default.
!  15 Mar 2018 - M. Sulprizio- Move nc_dir definition to here from Read_Modis.
!                              We now use updates MODIS LAI files stored in
!                              MODIS_LAI_201707/, while CHLR files are still
!                              read from MODIS_LAI_201204/.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)   :: nc_dir          ! netCDF directory name
    CHARACTER(LEN=255)   :: nc_tmpl         ! netCDF file name template 
    LOGICAL              :: ReadLAI         ! T = read LAI, F = read CHLR

    ! Assume success
    RC = GC_SUCCESS

    ! Construct file path from directory & file name (use LAI directory)
    nc_dir  = TRIM( Input_Opt%CHEM_INPUTS_DIR ) // 'MODIS_LAI_201707/'

    ! Set filename template for LAI
    nc_tmpl = 'For_Olson_2001/MODIS.LAIv.V5.generic.025x025.YYYY.nc'
    nc_tmpl = TRIM( nc_dir ) // TRIM( nc_tmpl )

    ! Always read LAI file
    ReadLAI = .true.
    CALL Read_Modis( am_I_Root, ReadLAI, nc_tmpl,      Input_Opt,  &
                     yyyy,      mm,      wasModisRead, RC         )

    ! Read CHLR only if organic marine aerosols tracers are turned on
    IF ( Input_Opt%LMPOA ) THEN

       ! Construct file path from directory & file name (use LAI directory)
       nc_dir  = TRIM( Input_Opt%CHEM_INPUTS_DIR ) // 'MODIS_LAI_201204/'

       ! Set filename template for CHLR
       nc_tmpl = 'For_Olson_2001/MODIS.CHLRv.V5.generic.025x025.YYYY.nc'
       nc_tmpl = TRIM( nc_dir ) // TRIM( nc_tmpl )

       ! Read CHLR file
       ReadLAI = .false.
       CALL Read_Modis(  am_I_Root, ReadLAI, nc_tmpl,      Input_Opt,  &
                         yyyy,      mm,      wasModisRead, RC         )
    ENDIF

  END SUBROUTINE Read_Modis_Lai
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_modis
!
! !DESCRIPTION: Subroutine READ\_MODIS reads the MODIS LAI or CHLR from disk
!  (in netCDF format) for the current month, and for next month.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Read_Modis( am_I_Root, ReadLAI, nc_tmpl,      Input_Opt,  &
                         yyyy,      mm,      wasModisRead, RC         )
!
! !USES:
!
    USE m_netcdf_io_open                         ! netCDF file open
    USE m_netcdf_io_read                         ! netCDF read
    USE m_netcdf_io_readattr                     ! netCDF attribute reads
    USE m_netcdf_io_close                        ! netCDF file close
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput 
   
#   include "netcdf.inc"                         ! netCDF settings & parameters
!
! !INPUT PARAMETERS:
!
    LOGICAL,            INTENT(IN)  :: am_I_Root ! Are we on the root CPU?
    LOGICAL,            INTENT(IN)  :: ReadLAI   ! T for LAI, F to read CHLR
    TYPE(OptInput),     INTENT(IN)  :: Input_Opt ! Input Options object
    CHARACTER(LEN=255), INTENT(IN)  :: nc_tmpl
    INTEGER,            INTENT(IN)  :: yyyy      ! Year for LAI data
    INTEGER,            INTENT(IN)  :: mm        ! Month for LAI data
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,        INTENT(OUT) :: wasModisRead  ! Was data just read in?
    INTEGER,        INTENT(OUT) :: RC            ! Success or failure?   
!
! !REVISION HISTORY:
!  07 Jul 2015 - E. Lundgren - Initial version, containing legacy Read_Modis_Lai
!                              code plus modifications to read CHLR
!  13 Sep 2017 - M. Sulprizio- Remove Input_Opt%USE_OLSON_2001. Olson 2001 is
!                              now the default.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: fId                ! netCDF file ID
    INTEGER              :: Pyyyy              ! Previous LAI year
    INTEGER              :: Nyyyy              ! Next     LAI year
    INTEGER              :: Pmm                ! Previous LAI month
    INTEGER              :: Nmm                ! Next     LAI month
    INTEGER              :: yyyymmdd           ! Date variable
    INTEGER              :: MODIS_END          ! Last year or MODIS data
                         
    ! Character strings  
    CHARACTER(LEN=255)   :: nc_file            ! netCDF file name
    CHARACTER(LEN=255)   :: v_name             ! netCDF variable name 
    CHARACTER(LEN=255)   :: a_name             ! netCDF attribute name
    CHARACTER(LEN=255)   :: a_val              ! netCDF attribute value

    ! Arrays for netCDF  start and count values
    INTEGER              :: st1d(1), ct1d(1)   ! For 1D arrays    
    INTEGER              :: st3d(3), ct3d(3)   ! For 3D arrays 

    ! SAVED variable
    INTEGER, SAVE        :: mmLastLAI  = -1
    INTEGER, SAVE        :: mmLastCHLR = -1

    ! Pointers
    REAL*4, POINTER      :: MODIS_PTR_CM(:,:) => NULL()
    REAL*4, POINTER      :: MODIS_PTR_NM(:,:) => NULL()

    !======================================================================
    ! Test if it is time to read data
    !======================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Assign pointers etc. based on whether reading LAI or CHLR
    IF ( ReadLAI ) THEN
       MODIS_PTR_CM => MODIS_LAI_CM
       MODIS_PTR_NM => MODIS_LAI_NM
       MODIS_END    =  LAI_END
    ELSE
       MODIS_PTR_CM => MODIS_CHLR_CM
       MODIS_PTR_NM => MODIS_CHLR_NM
       MODIS_END    =  CHLR_END
    ENDIF

    ! If we enter a new month, then read MODIS from disk
    ! Otherwise, just exit, since it is not time to read data yet
    IF ( ReadLAI .and. mm /= mmLastLAI ) THEN
       mmLastLAI    = mm
       wasModisRead = .TRUE.
    ELSEIF ( ( .not. ReadLAI ) .and. mm /= mmLastCHLR ) THEN
       mmLastCHLR   = mm
       wasModisRead = .TRUE.
    ELSE
       wasModisRead = .FALSE.
       RETURN
    ENDIF

    ! Save for next iteration

    !======================================================================
    ! Read current month's LAI or CHLR
    !======================================================================

    ! Test if yyyy is w/in the valid range of MODIS data
    IF ( yyyy >= MODIS_START .and. yyyy <= MODIS_END ) THEN

       ! Here, yyyy lies w/in the MODIS data timespan
       nc_file  = nc_tmpl
       yyyymmdd = yyyy*10000 + mm*100 + 01
       CALL Expand_Date( nc_file, yyyymmdd, 000000 )

    ELSE

       ! Here, yyyy lies outside the MODIS data timespan,
       ! so we have to read data from a different year.
       IF ( yyyy > MODIS_END ) THEN
          yyyymmdd = MODIS_END*10000   + mm*100 + 01   ! Use final year
       ELSE IF ( yyyy < MODIS_START ) THEN             !
          yyyymmdd = MODIS_START*10000 + mm*100 + 01   ! Use 1st year
       ENDIF

       ! Expand date tokens in filename
       nc_file  = nc_tmpl
       CALL Expand_Date( nc_file, yyyymmdd, 000000 )

    ENDIF

    ! Open file for read
    CALL Ncop_Rd( fId, TRIM(nc_file) )  

    ! Echo info to stdout
    WRITE( 6, 100 ) REPEAT( '%', 79 )
    WRITE( 6, 110 ) TRIM(nc_file)

    ! Variable name
    v_name = "MODIS"
    
    ! Read OLSON from file
    st3d   = (/ 1,       1,       mm /)
    ct3d   = (/ I_MODIS, J_MODIS, 1  /)
    CALL NcRd( MODIS_PTR_CM, fId, TRIM(v_name), st3d, ct3d )

    ! Read the OLSON:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
    
    ! Echo info to stdout
    WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val), mm

    ! Close netCDF file
    CALL NcCl( fId )
    
    ! Echo info to stdout
    WRITE( 6, 140 )
 
    !======================================================================
    ! Read next month's LAI or CHLR
    !======================================================================

    ! Next LAI month
    Nmm = mm + 1

    ! Year corresponding to next LAI month (readjust for December)
    IF ( Nmm == 13 ) THEN
       Nmm   = 1
       Nyyyy = yyyy + 1
    ELSE
       Nyyyy = yyyy
    ENDIF

    ! Test if Nyyy is w/in the valid range of MODIS data
    IF ( Nyyyy >= MODIS_START .and. Nyyyy <= MODIS_END ) THEN

       ! Here, Nyyyy lies w/in the MODIS data timespan
       nc_file  = nc_tmpl
       yyyymmdd = Nyyyy*10000 + Nmm*100 + 01
       CALL Expand_Date( nc_file, yyyymmdd, 000000 )

    ELSE

       ! Here, yyyy lies outside the MODIS data timespan,
       ! so we have to read data from a different year.
       IF ( Nyyyy > MODIS_END ) THEN
          yyyymmdd = MODIS_END*10000   + Nmm*100 + 01  ! Use final year
       ELSE IF ( Nyyyy < MODIS_START ) THEN            !
          yyyymmdd = MODIS_START*10000 + Nmm*100 + 01  ! Use 1st year
       ENDIF

       ! Expand date tokens in filename
       nc_file  = nc_tmpl
       CALL Expand_Date( nc_file, yyyymmdd, 000000 )
          
    ENDIF

    ! Open file for read
    CALL Ncop_Rd( fId, TRIM(nc_file) )
     
    ! Echo info to stdout
    WRITE( 6, 100 ) REPEAT( '%', 79 )
    WRITE( 6, 110 ) TRIM(nc_file)

    ! Variable name
    v_name = "MODIS"
    
    ! Read OLSON from file
    st3d   = (/ 1,       1,       Nmm /)
    ct3d   = (/ I_MODIS, J_MODIS, 1   /)
    CALL NcRd( MODIS_PTR_NM, fId, TRIM(v_name), st3d, ct3d )

    ! Read the OLSON:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
    
    ! Echo info to stdout
    WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val), Nmm

    ! Close netCDF file
    CALL NcCl( fId )
    
    ! Echo info to stdout
    WRITE( 6, 140 )
    WRITE( 6, 100 ) REPEAT( '%', 79 )

    ! FORMAT statements
100 FORMAT( a                                                                 )
110 FORMAT( '%% Opening file  : ',         a                                  )
120 FORMAT( '%%  in directory : ',         a, / , '%%'                        )
130 FORMAT( '%% Successfully read ',       a, ' [', a, '] for month = ', i2.2 )
140 FORMAT( '%% Successfully closed file!'                                    )

    ! Cleanup pointers
    MODIS_PTR_CM => NULL()
    MODIS_PTR_NM => NULL()

  END SUBROUTINE Read_Modis
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: find_lai_month
!
! !DESCRIPTION: Function FIND\_LAI\_MONTH returns the corresponding LAI 
!  month and year for the current calendar date.  Note that the LAI data
!  starts at mid-month.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Find_Lai_Month( doy, month, year, mm, yyyy )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN)  :: doy       ! Current day of year
    INTEGER, INTENT(IN)  :: month     ! Current month
    INTEGER, INTENT(IN)  :: year      ! Current year
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: mm        ! Output month for LAI data
    INTEGER, INTENT(OUT) :: yyyy      ! Output year  for LAI data
!
! !REVISION HISTORY: 
!  05 Jan 1994 - Y. H. Wang, G.M. Gardner, D. Jacob - Initial version
!  (1 ) Updated comments, cosmetic changes (bmy, 4/4/03)
!  (2 ) Add the current simulation year as input & the current LAI as output.
!       This is necessary for reading in MODIS LAI (mpb,2009).
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  03 Apr 2012 - R. Yantosca - Renamed to FIND_LAI_MONTH; made PUBLIC
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( doy < startDay(1) ) THEN

       ! LAI month & year in December of the preceding calendar year
       mm   = 12
       yyyy = year - 1       

    ELSE IF ( doy < startDay(month) ) THEN

       ! LAI month is the preceding month of this year
       mm   = month - 1
       yyyy = year           

    ELSE

       ! LAI month is the the current month of this year
       mm   = month
       yyyy = year           

    ENDIF

  END SUBROUTINE Find_Lai_Month
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: RoundOff
!
! !DESCRIPTION: Rounds a number X to N decimal places of precision.
!\\
!\\
! !INTERFACE:
!
  FUNCTION RoundOff( X, N ) RESULT( Y )
!
! !INPUT PARAMETERS:
! 
    REAL(fp),  INTENT(IN) :: X   ! Number to be rounded
    INTEGER, INTENT(IN)   :: N   ! Number of decimal places to keep
!
! !RETURN VALUE:
!
    REAL(fp)              :: Y   ! Number rounded to N decimal places
!
! !REMARKS:
!  The algorithm to round X to N decimal places is as follows:
!  (1) Multiply X by 10**(N+1)
!  (2) If X < 0, then add -5 to X; otherwise add 5 to X
!  (3) Take the integer part of X
!  (4) Divide X by 10**(N+1)
!  (5) Truncate X to N decimal places: INT( X * 10**N ) / 10**N
!                                                                             .
!  Rounding algorithm from: Hultquist, P.F, "Numerical Methods for Engineers 
!   and Computer Scientists", Benjamin/Cummings, Menlo Park CA, 1988, p. 20.
!                                                                             .
!  Truncation algorithm from: http://en.wikipedia.org/wiki/Truncation
!                                                                             .
!  The two algorithms have been merged together for efficiency.
!
! !REVISION HISTORY:
!  06 Apr 2012 - R. Yantosca - Initial version
!  09 Apr 2012 - R. Yantosca - Changed all variables & arguments to REAL(fp)
!EOP
!------------------------------------------------------------------------------
!BOC
!f
! !LOCAL VARIABLES
!
    REAL(fp) :: TEN_TO_THE_N                   ! Term for 10**N
    REAL(fp) :: TEN_TO_THE_Np1                 ! Term for 10**(N+1)

    ! Pre-compute exponential terms
    TEN_TO_THE_N   = 10e+0_fp**N
    TEN_TO_THE_Np1 = 10e+0_fp**(N+1)
    
    ! Steps (1) through (4) above
    Y = INT( ( X * TEN_TO_THE_Np1 ) + SIGN( 5e+0_fp, X ) ) / TEN_TO_THE_Np1

    ! Step (5) above
    Y = INT( Y * TEN_TO_THE_N ) / TEN_TO_THE_N
  
  END FUNCTION RoundOff
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_modis
!
! !DESCRIPTION: Subroutine INIT\_MODIS\_LAI initializes and allocates
!  all module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Modis_Lai( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Diag_Mod,     ONLY : DgnState
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  03 Apr 2012 - R. Yantosca - Initial version
!  03 Feb 2014 - M. Sulprizio- Force last year of MODIS data to 2008. There is
!                              a large difference in the 2009 file that still
!                              needs to be investigated (skim, 1/29/14)
!  23 Jun 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC 
!  08 Jul 2015 - E. Lundgren - New end years to match files from M. Johnson
!  13 Sep 2017 - M. Sulprizio- Remove Input_Opt%USE_OLSON_2001. Olson 2001 is
!                              now the default.
!  15 Mar 2018 - M. Sulprizio- Update MODIS LAI end year to 2011.
!  07 Aug 2018 - H.P. Lin    - Now accepts State_Chm, State_Diag to unify input
!  23 Jan 2019 - H.P. Lin    - Do not allocate MODIS grid grid vars if already done
!EOP
!------------------------------------------------------------------------------
!BOC
    ! If already initialized, do not allocate the "fine" MODIS grid
    ! (which is fixed and not changing between domain calls)
    IF ( ALLOCATED( MODIS_LAI ) ) THEN
      RETURN
    ENDIF

    !======================================================================
    ! Allocate arrays on the "fine" MODIS grid grid
    !======================================================================
    I_MODIS     = 1440             ! For Olson 2001, use MODIS LAI
    J_MODIS     = 720              ! on the 0.25 x 0.25 native grid
    MODIS_START = 2005             ! First year of MODIS data  
    LAI_END     = 2011             ! Last  year of MODIS data
    CHLR_END    = 2011             ! Last  year of MODIS CHLR data

    ALLOCATE( MODIS_LAI( I_MODIS, J_MODIS ), STAT=RC ) 
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'MODIS_LAI' )
    MODIS_LAI = 0e+0_fp

    ALLOCATE( MODIS_LAI_CM( I_MODIS, J_MODIS ), STAT=RC ) 
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'MODIS_LAI_CM' )
    MODIS_LAI_CM = 0e+0_fp

    ALLOCATE( MODIS_LAI_NM( I_MODIS, J_MODIS ), STAT=RC ) 
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'MODIS_LAI_NM' )
    MODIS_LAI_NM = 0e+0_fp

    ALLOCATE( MODIS_CHLR( I_MODIS, J_MODIS ), STAT=RC ) 
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'MODIS_CHLR' )
    MODIS_LAI = 0e+0_fp

    ALLOCATE( MODIS_CHLR_CM( I_MODIS, J_MODIS ), STAT=RC ) 
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'MODIS_CHLR_CM' )
    MODIS_LAI_CM = 0e+0_fp

    ALLOCATE( MODIS_CHLR_NM( I_MODIS, J_MODIS ), STAT=RC ) 
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'MODIS_CHLR_NM' )
    MODIS_LAI_NM = 0e+0_fp

  END SUBROUTINE Init_Modis_Lai
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_modis_lai
!
! !DESCRIPTION: Subroutine CLEANUP\_MODIS\_LAI deallocates all
!  previously-allocated module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Modis_Lai
!
! !REVISION HISTORY:'
!  03 Apr 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( ALLOCATED( MODIS_LAI     ) ) DEALLOCATE( MODIS_LAI     )
    IF ( ALLOCATED( MODIS_LAI_CM  ) ) DEALLOCATE( MODIS_LAI_CM  )
    IF ( ALLOCATED( MODIS_LAI_NM  ) ) DEALLOCATE( MODIS_LAI_NM  )
    IF ( ALLOCATED( MODIS_CHLR    ) ) DEALLOCATE( MODIS_CHLR    )
    IF ( ALLOCATED( MODIS_CHLR_CM ) ) DEALLOCATE( MODIS_CHLR_CM )
    IF ( ALLOCATED( MODIS_CHLR_NM ) ) DEALLOCATE( MODIS_CHLR_NM )

  END SUBROUTINE Cleanup_Modis_Lai
!EOC
END MODULE Modis_Lai_Mod
