#ifdef BPCH_DIAG
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag51b_mod.F90
!
! !DESCRIPTION: Module DIAG51b\_MOD contains variables and routines to
!  generate save timeseries data where the local time is between two
!  user-defined limits. This facilitates comparisons with morning or
!  afternoon-passing satellites such as GOME.
!\\
!\\
! !INTERFACE:
!
MODULE DIAG51b_MOD
!
! !USES:
!
  USE PhysConstants, ONLY : AIRMW
  USE PRECISION_MOD       ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: CLEANUP_DIAG51b
  PUBLIC  :: DIAG51b
  PUBLIC  :: INIT_DIAG51b
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: ACCUMULATE_DIAG51b
  PRIVATE :: GET_LOCAL_TIME
  PRIVATE :: ITS_TIME_FOR_WRITE_DIAG51b
  PRIVATE :: WRITE_DIAG51b
!
! !REMARKS:
!  NOTE by Melissa Sulprizio, 26 May 2015
!  ----------------------------------------------------------------------------
!  The emission options in the timeseries diagnostics were removed
!  from v10-01 since HEMCO now handles the emission diagnostics. To
!  utilize the diagnostics capability from HEMCO and output hourly
!  isoprene emissions, you can follow these steps:
!
!    1. At the top of your HEMCO_Config.rc file, set DiagnFreq to Hourly
!       and add a line for DiagnFile:
!
!          DiagnPrefix: HEMCO_Diagnostics
!          DiagnFreq: Hourly
!          DiagnFile: DiagnFile.rc<
!
!    2.  Create a new text file in your run directory named DiagnFile.rc
!        and list the emission fields that you would like to be saved out.
!        For example:
!
!          # Name        Spec ExtNr Cat Hier Dim OutUnit
!          ISOP_BIOG     ISOP 108    1   1   2   kg/m2/s
!
!        NOTE: The ExtNr, Cat, Hier, and Dim values listed above were
!        obtained from the MEGAN entries in the HEMCO_Config.rc file.
!
!    3. You can then run GEOS-Chem as usual. HEMCO will write out the
!       specified diagnostics in a netCDF file named
!       HEMCO_Diagnostics.YYYYMMDDHHmm.nc. I recommend running a short
!       1-day simulation to make sure the diagnostic output is what
!       you expect.
!
!  For more details on the HEMCO diagnostics, please see this post
!  in the HEMCO User’s Guide:
!
!       http://wiki.geos-chem.org/The_HEMCO_User%27s_Guide#Diagnostics
!
!  ND51b tracer numbers:
!  ============================================================================
!  1 - nAdvect   : GEOS-CHEM advected species               [v/v        ]
!  501           : OH concentration                         [molec/cm3  ]
!  502           : NOy concentration                        [v/v        ]
!  503           : Relative Humidity                        [%          ]
!  504           : 3-D Cloud fractions                      [unitless   ]
!  505           : Column optical depths                    [unitless   ]
!  506           : Cloud top heights                        [hPa        ]
!  507           : Air density                              [molec/cm3  ]
!  508           : Total seasalt tracer concentration       [unitless   ]
!  509           : PBL heights                              [m          ]
!  510           : PBL heights                              [levels     ]
!  511           : Grid box heights                         [m          ]
!  512           : PEDGE-$ (Pressure @ level edges          [hPa        ]
!  513           : Sea level pressure                       [hPa        ]
!  514           : Zonal wind (a.k.a. U-wind)               [m/s        ]
!  515           : Meridional wind (a.k.a. V-wind)          [m/s        ]
!  516           : Temperature                              [K          ]
!  517           : Sulfate aerosol optical depth            [unitless   ]
!  518           : Black carbon aerosol optical depth       [unitless   ]
!  519           : Organic carbon aerosol optical depth     [unitless   ]
!  520           : Accumulation mode seasalt optical depth  [unitless   ]
!  521           : Coarse mode seasalt optical depth        [unitless   ]
!  522           : Total dust optical depth                 [unitless   ]
!  523-529       : Size resolved dust optical depth         [unitless   ]
!  530           : PAR direct                               [hPa        ]
!  531           : PAR diffuse                              [hPa        ]
!  532           : Daily LAI                                [hPa        ]
!  533           : Temperature at 2m                        [K          ]
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!

  !=================================================================
  ! MODULE VARIABLES
  ! GOOD             : Array denoting grid boxes w/in LT limits
  ! GOOD_CT          : # of "good" times per grid box
  ! GOOD_CT_CHEM     : # of "good" chemistry timesteps
  ! COUNT_CHEM3D     : Counter for 3D chemistry boxes
  ! I0               : Offset between global & nested grid
  ! J0               : Offset between global & nested grid
  ! IOFF             : Longitude offset
  ! JOFF             : Latitude offset
  ! LOFF             : Altitude offset
  ! ND51b_NI         : Number of longitudes in DIAG51b region
  ! ND51b_NJ         : Number of latitudes  in DIAG51b region
  ! ND51b_NL         : Number of levels     in DIAG51b region
  ! Q                : Accumulator array for various quantities
  ! TAU0             : Starting TAU used to index the bpch file
  ! TAU1             : Ending TAU used to index the bpch file
  ! HALFPOLAR        : Used for bpch file output
  ! CENTER180        : Used for bpch file output
  ! LONRES           : Used for bpch file output
  ! LATRES           : Used for bpch file output
  ! MODELNAME        : Used for bpch file output
  ! RESERVED         : Used for bpch file output
  !=================================================================
  ! Scalars
  INTEGER              :: IOFF,           JOFF,     LOFF
  INTEGER              :: I0,             J0
  INTEGER              :: ND51b_NI,       ND51b_NJ, ND51b_NL
  INTEGER              :: HALFPOLAR
  INTEGER, PARAMETER   :: CENTER180=1
  REAL*4               :: LONRES,         LATRES
  REAL(f8)             :: TAU0,           TAU1
  CHARACTER(LEN=20)    :: MODELNAME
  CHARACTER(LEN=40)    :: RESERVED = ''
  CHARACTER(LEN=80)    :: TITLE

  ! LUN for ND51b output file
  INTEGER              :: IU_ND51b

  ! Arrays
  INTEGER, ALLOCATABLE :: GOOD(:)
  INTEGER, ALLOCATABLE :: GOOD_CT(:)
  INTEGER, ALLOCATABLE :: COUNT_CHEM3D(:,:,:)
  REAL(fp),ALLOCATABLE :: Q(:,:,:,:)

  !=================================================================
  ! MODULE ROUTINES -- follow below the "CONTAINS" statement
  !=================================================================
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diag51b
!
! !DESCRIPTION:  Subroutine DIAG51b generates time series (averages from !
!  10am - 12pm LT or 1pm - 4pm LT) for the US grid area.  Output is to
!  binary punch files or HDF5 files.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DIAG51b( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
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
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp) :: TAU_W

    !=================================================================
    ! DIAG51b begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Construct array of where local times are between HR1, HR2
    CALL GET_LOCAL_TIME( Input_Opt, State_Grid )

    ! Accumulate data in the Q array
    CALL ACCUMULATE_DIAG51b( Input_Opt, State_Chm, State_Grid, State_Met, RC )

    ! Write data to disk at the proper time
    IF ( ITS_TIME_FOR_WRITE_DIAG51b( Input_Opt, TAU_W ) .and. &
         Input_Opt%DO_DIAG_WRITE ) THEN
       CALL WRITE_DIAG51b( Input_Opt, State_Chm, State_Grid, TAU_W, RC)
    ENDIF

  END SUBROUTINE DIAG51b
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_local_time
!
! !DESCRIPTION: Subroutine GET\_LOCAL\_TIME computes the local time and
!  returns an array of points where the local time is between two user-defined
!  limits.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_LOCAL_TIME( Input_Opt, State_Grid )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
    USE TIME_MOD,       ONLY : GET_LOCALTIME
    USE TIME_MOD,       ONLY : GET_TS_DYN
!
! !INPUT ARGUMENTS
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt  ! Input Options
    TYPE(GrdState), INTENT(IN) :: State_Grid ! Grid State object
!
! !REMARKS:
!  For now use GET_LOCALTIME( I, 1, 1 ) which will be independent of J and L
!  for a pure cartesian grid.  This may need to be revisited once G-C is
!  interfaced into a GCM.
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I
    REAL(fp) :: LT, TS_DYN

    !=================================================================
    ! GET_LOCAL_TIME begins here!
    !=================================================================
    TS_DYN = GET_TS_DYN() / 3600e+0_fp

    DO I = 1, State_Grid%NX

       ! Get local time
       LT = GET_LOCALTIME( I, 1, 1, State_Grid ) - TS_DYN
       IF ( LT < 0  ) LT = LT + 24e+0_fp

       ! GOOD indicates which boxes have local times between HR1 and HR2
       IF ( LT >= Input_Opt%ND51b_HR1 .and. &
            LT <= Input_Opt%ND51b_HR2 ) THEN
          GOOD(I) = 1
       ENDIF
    ENDDO

  END SUBROUTINE GET_LOCAL_TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: accumulate_diag51b
!
! !DESCRIPTION: Subroutine ACCUMULATE\_DIAG51b accumulates tracers into the
!  Q array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ACCUMULATE_DIAG51b( Input_Opt, State_Chm, &
                                 State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE CMN_FJX_MOD,        ONLY : ODAER, ODMDUST
    USE CMN_FJX_MOD,        ONLY : IWVSELECT, ACOEF_WV, BCOEF_WV
    USE CMN_O3_MOD               ! SAVEOH
    USE CMN_SIZE_MOD,       ONLY : NRH, NDUST
    USE PhysConstants            ! SCALE_HEIGHT, XNUMOLAIR
    USE TIME_MOD,           ONLY : GET_ELAPSED_SEC,  GET_TS_CHEM
    USE TIME_MOD,           ONLY : TIMESTAMP_STRING, GET_TS_DYN
    USE TIME_MOD,           ONLY : GET_TS_DIAG,      GET_TS_EMIS
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
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL, SAVE     :: FIRST = .TRUE.
    LOGICAL, SAVE     :: IS_FULLCHEM, IS_SEASALT
    LOGICAL, SAVE     :: IS_CLDTOPS,  IS_NOy,    IS_OPTD, IS_SLP
    INTEGER, SAVE     :: id_HNO3,     id_HNO4,   id_N2O5, id_NO
    INTEGER, SAVE     :: id_PAN,      id_MPAN,   id_PPN,  id_O3
    INTEGER, SAVE     :: id_R4N2,     id_SALA,   id_SALC, id_NO2
    INTEGER, SAVE     :: id_OH

    ! Scalars
    LOGICAL           :: IS_CHEM,     IS_DIAG,   IS_EMIS
    LOGICAL           :: LINTERP
    INTEGER           :: H, I, J, K, L, M, N, ISPC, NA, nAdvect
    INTEGER           :: PBLINT,  R, X, Y, W, XMIN
    INTEGER           :: IWV
    REAL(fp)          :: C1, C2, PBLDEC, TEMPBL, TMP, SCALEAODnm
    CHARACTER(LEN=16) :: STAMP

    ! Aerosol types (rvm, aad, bmy, 7/20/04)
    INTEGER           :: IND(6) = (/ 22, 29, 36, 43, 50, 15 /)

    ! Pointers
    REAL(fp), POINTER :: Spc(:,:,:,:)

    !=================================================================
    ! ACCUMULATE_DIAG51b begins here!
    !=================================================================

    ! Assume success
    RC        =  GC_SUCCESS

    ! Number of advected species
    nAdvect   =  State_Chm%nAdvect

    ! First-time setup
    IF ( FIRST ) THEN

       ! Define species ID flags
       id_HNO3 = Ind_('HNO3')
       id_HNO4 = Ind_('HNO4')
       id_N2O5 = Ind_('N2O5')
       id_NO   = Ind_('NO')
       id_PAN  = Ind_('PAN')
       id_MPAN = Ind_('MPAN')
       id_PPN  = Ind_('PPN')
       id_O3   = Ind_('O3')
       id_R4N2 = Ind_('R4N2')
       id_SALA = Ind_('SALA')
       id_SALC = Ind_('SALC')
       id_NO2  = Ind_('NO2')
       id_OH   = Ind_('OH')

       ! Set logical flags on first call
       IS_OPTD     = ASSOCIATED( State_Met%OPTD    )
       IS_CLDTOPS  = ASSOCIATED( State_Met%CLDTOPS )
       IS_SLP      = ASSOCIATED( State_Met%SLP     )
       IS_FULLCHEM = Input_Opt%ITS_A_FULLCHEM_SIM
       IS_SEASALT  = ( id_SALA > 0 .and. id_SALC > 0 )
       IS_NOy      = ( IS_FULLCHEM .and. id_NO   > 0 .and. &
                       id_NO2  > 0 .and. id_PAN  > 0 .and. &
                       id_HNO3 > 0 .and. id_MPAN > 0 .and. &
                       id_PPN  > 0 .and. id_R4N2 > 0 .and. &
                       id_N2O5 > 0 .and. id_HNO4 > 0 )

       ! Reset first-time flag
       FIRST       = .FALSE.
    ENDIF

    ! Force accumulation on every "heartbeat" timestep (bmy, 10/4/18)
    IS_DIAG = .TRUE.

    ! Echo info
    STAMP = TIMESTAMP_STRING()
    WRITE( 6, 100 ) STAMP
100 FORMAT( '     - DIAG51b: Accumulation at ', a )

    !=================================================================
    ! Archive tracers into accumulating array Q
    !=================================================================

    ! Archive counter array of good points
    IF ( IS_DIAG ) THEN
       DO X = 1, ND51b_NI
          I          = GET_I( X, State_Grid )
          GOOD_CT(X) = GOOD_CT(X) + GOOD(I)
       ENDDO
    ENDIF

    ! Also increment 3-D counter for boxes in the chemistry grid
    IF ( IS_FULLCHEM ) THEN

       ! Loop over levels
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( X, Y, K, I, J, L ) &
       !$OMP SCHEDULE( DYNAMIC )
       DO K = 1, ND51b_NL
          L = LOFF + K

       ! Loop over latitudes
       DO Y = 1, ND51b_NJ
          J = JOFF + Y

       ! Loop over longitudes
       DO X = 1, ND51b_NI
          I = GET_I( X, State_Grid )

          ! Only increment if we are in the
          IF (  State_Met%InChemGrid(I,J,L) ) THEN
             COUNT_CHEM3D(X,Y,K) = COUNT_CHEM3D(X,Y,K) + GOOD(I)
          ENDIF

       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    !------------------------
    ! Accumulate quantities
    !------------------------
    IF ( IS_DIAG ) THEN
       !Determine if optical properties need interpolating
       !the LUT wavelengths in IWVSELECT will match if no interpolation
       !is needed. (output is only for the first requested wavelength)
       IF(IWVSELECT(1,1).EQ.IWVSELECT(2,1)) THEN
          LINTERP=.FALSE.
       ELSE
          LINTERP=.TRUE.
       ENDIF

       ! Initialize pointers
       Spc => State_Chm%Species

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( W, N, X, Y, K, I, J, L, TMP, H, R, SCALEAODnm ) &
       !$OMP SCHEDULE( DYNAMIC )
       DO W = 1, Input_Opt%N_ND51b

          ! ND51b Tracer number
          N = Input_Opt%ND51b_TRACERS(W)

          ! Loop over levels
          DO K = 1, ND51b_NL
             L = LOFF + K

          ! Loop over latitudes
          DO Y = 1, ND51b_NJ
             J = JOFF + Y

          ! Loop over longitudes
          DO X = 1, ND51b_NI
             I = GET_I( X, State_Grid )

             ! Archive by simulation
             IF ( N <= nAdvect ) THEN

                !--------------------------------------
                ! GEOS-CHEM advected species [v/v]
                !--------------------------------------

                ! Archive afternoon points
                Q(X,Y,K,W) = Q(X,Y,K,W) + &
                           ( Spc(I,J,L,N) * ( AIRMW &
                           / State_Chm%SpcData(N)%Info%emMW_g ) &
                             * GOOD(I) )

             ELSE IF ( N == 502 .and. IS_NOy ) THEN

                !--------------------------------------
                ! NOy [v/v]
                !--------------------------------------

                ! Temp variable for accumulation
                TMP = 0e+0_fp

                ! NO
                TMP = TMP + ( ( AIRMW &
                          / State_Chm%SpcData(id_NO)%Info%emMW_g ) &
                          * GOOD(I) * Spc(I,J,L,id_NO)    )

                ! NO2
                TMP = TMP + ( ( AIRMW &
                          / State_Chm%SpcData(id_NO2)%Info%emMW_g ) &
                          * GOOD(I) * Spc(I,J,L,id_NO2)   )
                ! PAN
                TMP = TMP + ( ( AIRMW &
                          / State_Chm%SpcData(id_PAN)%Info%emMW_g ) &
                          * GOOD(I) * Spc(I,J,L,id_PAN)   )

                ! HNO3
                TMP = TMP + ( ( AIRMW &
                          / State_Chm%SpcData(id_HNO3)%Info%emMW_g ) &
                          * GOOD(I) * Spc(I,J,L,id_HNO3)  )

                ! MPAN
                TMP = TMP + ( ( AIRMW &
                          / State_Chm%SpcData(id_MPAN)%Info%emMW_g ) &
                          * GOOD(I) * Spc(I,J,L,id_MPAN)   )

                ! PPN
                TMP = TMP + ( ( AIRMW &
                          / State_Chm%SpcData(id_PPN)%Info%emMW_g ) &
                          * GOOD(I) * Spc(I,J,L,id_PPN)   )

                ! R4N2
                TMP = TMP + ( ( AIRMW &
                          / State_Chm%SpcData(id_R4N2)%Info%emMW_g ) &
                          * GOOD(I) * Spc(I,J,L,id_R4N2)  )

                ! N2O5
                TMP = TMP + ( 2e+0_fp * ( AIRMW &
                          / State_Chm%SpcData(id_N2O5)%Info%emMW_g ) &
                          * GOOD(I) * Spc(I,J,L,id_N2O5)  )

                ! HNO4
                TMP = TMP + ( ( AIRMW &
                          / State_Chm%SpcData(id_HNO4)%Info%emMW_g ) &
                          * GOOD(I) * Spc(I,J,L,id_HNO4)  )

                ! Save afternoon points
                Q(X,Y,K,W) = Q(X,Y,K,W) + TMP

             ELSE IF ( N == 503 ) THEN

                !---------------------------------------
                ! RELATIVE HUMIDITY [%]
                !---------------------------------------
                Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%RH(I,J,L) * GOOD(I) )

             ELSE IF ( N == 504 ) THEN

                !--------------------------------------
                ! 3-D CLOUD FRACTIONS [unitless]
                !--------------------------------------
                Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%CLDF(I,J,L) * GOOD(I) )

             ELSE IF ( N == 505 .and. IS_OPTD ) THEN

                !--------------------------------------
                ! COLUMN OPTICAL DEPTHS [unitless]
                !--------------------------------------
                Q(X,Y,1,W) = Q(X,Y,1,W) + ( State_Met%OPTD(I,J,L) * GOOD(I) )

             ELSE IF ( N == 506 .and. IS_CLDTOPS ) THEN

                !--------------------------------------
                ! CLOUD TOP HEIGHTS [mb]
                !--------------------------------------
                IF ( K == 1 ) THEN
                   TMP      = State_Met%PEDGE(I,J,State_Met%CLDTOPS(I,J))
                   Q(X,Y,K,W) = Q(X,Y,K,W) + ( TMP * GOOD(I) )
                ENDIF

             ELSE IF ( N == 507 ) THEN

                !--------------------------------------
                ! AIR DENSITY [molec/cm3]
                !--------------------------------------
                Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%AIRDEN(I,J,L) * &
                             XNUMOLAIR  * 1d-6 * GOOD(I) )

             ELSE IF ( N == 508 .and. IS_SEASALT ) THEN

                !--------------------------------------
                ! TOTAL SEASALT TRACER [v/v]
                !--------------------------------------
                Q(X,Y,K,W) = Q(X,Y,K,W) + &
                           ( Spc(I,J,L,id_SALA)   + &
                             Spc(I,J,L,id_SALC) ) * &
                             ( AIRMW / &
                               State_Chm%SpcData(id_SALA)%Info%emMW_g ) &
                               * GOOD(I)

             ELSE IF ( N == 509 ) THEN

                !--------------------------------------
                ! PBL HEIGHTS [m]
                !--------------------------------------
                IF ( K == 1 ) THEN
                   Q(X,Y,K,W) = Q(X,Y,K,W) + &
                              ( State_Met%PBL_TOP_m(I,J) * GOOD(I) )
                ENDIF

             ELSE IF ( N == 510 ) THEN

                !--------------------------------------
                ! PBL HEIGHTS [layers]
                !--------------------------------------
                IF ( K == 1 ) THEN
                   Q(X,Y,K,W) = Q(X,Y,K,W) + &
                              ( State_Met%PBL_TOP_L(I,J) * GOOD(I) )
                ENDIF

             ELSE IF ( N == 511 ) THEN

                !--------------------------------------
                ! GRID BOX HEIGHTS [m]
                !--------------------------------------
                Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%BXHEIGHT(I,J,L) * &
                             GOOD(I) )

             ELSE IF ( N == 512 ) THEN

                !--------------------------------------
                ! PEDGE-$ (prs @ level edges) [hPa]
                !--------------------------------------
                Q(X,Y,K,W) = Q(X,Y,K,W) + &
                           ( State_Met%PEDGE(I,J,K) * GOOD(I) )

             ELSE IF ( N == 513 .and. IS_SLP ) THEN

                !--------------------------------------
                ! SEA LEVEL PRESSURE [hPa]
                !--------------------------------------
                IF ( K == 1 ) THEN
                   Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%SLP(I,J) * GOOD(I) )
                ENDIF

             ELSE IF ( N == 514 ) THEN

                !--------------------------------------
                ! ZONAL (U) WIND [M/S]
                !--------------------------------------
                Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%U(I,J,L) * GOOD(I) )

             ELSE IF ( N == 515 ) THEN

                !--------------------------------------
                ! MERIDIONAL (V) WIND [M/S]
                !--------------------------------------
                Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%V(I,J,L) * GOOD(I) )

             ELSE IF ( N == 516 ) THEN

                !--------------------------------------
                ! TEMPERATURE [K]
                !--------------------------------------
                Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%T(I,J,L) * GOOD(I) )

             ELSE IF ( N == 517 ) THEN

                !--------------------------------------
                ! SULFATE AOD [unitless]
                ! for wavelengths set in Radiation Menu
                !
                ! NOTE: Only archive at chem timestep
                !--------------------------------------
                ISPC = 1 !sulfate
                IF ( .not. LINTERP ) THEN
                   ! Accumulate
                   Q(X,Y,K,W) = Q(X,Y,K,W) + &
                        (ODAER(I,J,L,IWVSELECT(1,1),ISPC) * GOOD(X))
                ELSE
                   ! Interpolated using angstrom exponent between
                   ! Closest available wavelengths
                   ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                   !catch any zero values before interpolation
                   IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND. &
                       (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                      Q(X,Y,K,W) = Q(X,Y,K,W) + GOOD(X)* &
                           (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)** &
                           (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/ &
                           ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                   ENDIF
                ENDIF

             ELSE IF ( N == 518 ) THEN

                !--------------------------------------
                ! BLACK CARBON AOD [unitless]
                ! for wavelengths set in Radiation Menu
                !
                ! NOTE: Only archive at chem timestep
                !--------------------------------------
                ISPC = 2 !BC
                IF ( .not. LINTERP ) THEN
                   ! Accumulate
                   Q(X,Y,K,W) = Q(X,Y,K,W) + &
                        (ODAER(I,J,L,IWVSELECT(1,1),ISPC) * GOOD(X))
                ELSE
                   ! Interpolated using angstrom exponent between
                   ! Closest available wavelengths
                   ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                   !catch any zero values before interpolation
                   IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND. &
                       (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                      Q(X,Y,K,W) = Q(X,Y,K,W) + GOOD(X)* &
                           (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)** &
                           (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/ &
                           ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                   ENDIF
                ENDIF

             ELSE IF ( N == 519 ) THEN

                !--------------------------------------
                ! ORG CARBON AOD [unitless]
                ! for wavelengths set in Radiation Menu
                !
                ! NOTE: Only archive at chem timestep
                !--------------------------------------
                ISPC = 3 !OC
                IF ( .not. LINTERP ) THEN
                   ! Accumulate
                   Q(X,Y,K,W) = Q(X,Y,K,W) + &
                        (ODAER(I,J,L,IWVSELECT(1,1),ISPC) * GOOD(X))
                ELSE
                   ! Interpolated using angstrom exponent between
                   ! Closest available wavelengths
                   ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                   !catch any zero values before interpolation
                   IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND. &
                       (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                      Q(X,Y,K,W) = Q(X,Y,K,W) + GOOD(X)* &
                           (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)** &
                           (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/ &
                           ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                   ENDIF
                ENDIF

             ELSE IF ( N == 520 ) THEN

                !--------------------------------------
                ! ACCUM SEASALT AOD [unitless]
                ! for wavelengths set in Radiation Menu
                !
                ! NOTE: Only archive at chem timestep
                !--------------------------------------
                ISPC = 4 !SSa
                IF ( .not. LINTERP ) THEN
                   ! Accumulate
                   Q(X,Y,K,W) = Q(X,Y,K,W) + &
                        (ODAER(I,J,L,IWVSELECT(1,1),ISPC) * GOOD(X))
                ELSE
                   ! Interpolated using angstrom exponent between
                   ! Closest available wavelengths
                   ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                   !catch any zero values before interpolation
                   IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND. &
                       (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                      Q(X,Y,K,W) = Q(X,Y,K,W) + GOOD(X)* &
                           (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)** &
                           (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/ &
                           ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                   ENDIF
                ENDIF

             ELSE IF ( N == 521 ) THEN

                !--------------------------------------
                ! COARSE SEASALT AOD [unitless]
                ! for wavelengths set in Radiation Menu
                !
                ! NOTE: Only archive at chem timestep
                !--------------------------------------
                ISPC = 5 !SSc
                IF ( .not. LINTERP ) THEN
                   ! Accumulate
                   Q(X,Y,K,W) = Q(X,Y,K,W) + &
                        (ODAER(I,J,L,IWVSELECT(1,1),ISPC) * GOOD(X))
                ELSE
                   ! Interpolated using angstrom exponent between
                   ! Closest available wavelengths
                   ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                   !catch any zero values before interpolation
                   IF ((ODAER(I,J,L,IWVSELECT(2,1),ISPC).GT.0).AND. &
                       (ODAER(I,J,L,IWVSELECT(1,1),ISPC).GT.0)) THEN
                      Q(X,Y,K,W) = Q(X,Y,K,W) + GOOD(X)* &
                           (ODAER(I,J,L,IWVSELECT(2,1),ISPC)*ACOEF_WV(1)** &
                           (BCOEF_WV(1)*LOG(ODAER(I,J,L,IWVSELECT(1,1),ISPC)/ &
                           ODAER(I,J,L,IWVSELECT(2,1),ISPC))))
                   ENDIF
                ENDIF

             ELSE IF ( N == 522 ) THEN

                !--------------------------------------
                ! TOTAL DUST OPTD [unitless]
                ! for wavelengths set in Radiation Menu
                !
                ! NOTE: Only archive at chem timestep
                !--------------------------------------
                DO R = 1, NDUST

                   IF ( .not. LINTERP ) THEN
                      Q(X,Y,K,W) = Q(X,Y,K,W) + GOOD(X) * &
                           ODMDUST(I,J,L,IWVSELECT(1,1),R)
                   ELSE
                      ! Interpolated using angstrom exponent between
                      ! Closest available wavelengths
                      ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                      ! AOD sometimes zero (if Q zero), must catch this
                      IF ((ODMDUST(I,J,L,IWVSELECT(1,1),R).GT.0).AND. &
                          (ODMDUST(I,J,L,IWVSELECT(2,1),R).GT.0)) THEN
                         Q(X,Y,K,W) = Q(X,Y,K,W) + GOOD(X) * &
                              (ODMDUST(I,J,L,IWVSELECT(2,1),R)* &
                              ACOEF_WV(1)**(BCOEF_WV(1)*LOG( &
                              ODMDUST(I,J,L,IWVSELECT(1,1),R)/ &
                              ODMDUST(I,J,L,IWVSELECT(2,1),R))))
                      ENDIF
                   ENDIF

                ENDDO

             ELSE IF ( ( N >= 523 ) .and. ( N <= 529) ) THEN

                !--------------------------------------
                ! Dust BINS 1-7 optical depth [unitless]
                ! for wavelengths set in Radiation Menu
                !
                ! NOTE: Only archive at chem timestep
                !--------------------------------------
                R = N - 522

                ! Accumulate
                IF ( .not. LINTERP ) THEN
                   Q(X,Y,K,W) = Q(X,Y,K,W) + GOOD(X) * &
                                ODMDUST(I,J,L,IWVSELECT(1,1),R)
                ELSE
                   ! Interpolated using angstrom exponent between
                   ! Closest available wavelengths
                   ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                   ! AOD sometimes zero (if Q zero), must catch this
                   IF ((ODMDUST(I,J,L,IWVSELECT(1,1),R).GT.0).AND. &
                       (ODMDUST(I,J,L,IWVSELECT(2,1),R).GT.0)) THEN
                      Q(X,Y,K,W) = Q(X,Y,K,W) + GOOD(X) * &
                           (ODMDUST(I,J,L,IWVSELECT(2,1),R)* &
                           ACOEF_WV(1)**(BCOEF_WV(1)*LOG( &
                           ODMDUST(I,J,L,IWVSELECT(1,1),R)/ &
                           ODMDUST(I,J,L,IWVSELECT(2,1),R))))
                   ENDIF
                ENDIF

             ELSE IF ( N == 530 ) THEN

                !--------------------------------------
                ! PAR DR [W/m2] (mpb,2009)
                !--------------------------------------

                IF ( K == 1 ) THEN
                   Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%PARDR(I,J) * GOOD(I) )
                ENDIF

             ELSE IF ( N == 531 ) THEN

                !--------------------------------------
                ! PAR DF [W/m2] (mpb,2009)
                !--------------------------------------

                IF ( K == 1 ) THEN
                   Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%PARDF(I,J) * GOOD(I) )
                ENDIF

             ELSE IF ( N == 532 ) THEN

                !--------------------------------------
                ! DAILY LAI  [cm2/cm2] (mpb,2009)
                !--------------------------------------

                IF ( K == 1 ) THEN
                   Q(X,Y,K,W) = Q(X,Y,K,W) + &
                              ( State_Met%MODISLAI(I,J) * GOOD(I) )
                ENDIF

             ELSE IF ( N == 533 ) THEN

                !--------------------------------------
                ! T at 2m [K] (mpb,2009)
                !--------------------------------------

                IF ( K == 1 ) THEN
                   Q(X,Y,K,W) = Q(X,Y,K,W) + ( State_Met%TS(I,J) * GOOD(I) )
                ENDIF

             ELSE

                ! Skip other tracers
                CYCLE

             ENDIF
          ENDDO
          ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL DO

       GOOD(:) = 0

       ! Free pointers
       Spc => NULL()

    ENDIF

  END SUBROUTINE ACCUMULATE_DIAG51b
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_write_diag51b
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_WRITE\_DIAG51b returns TRUE if it's
!  time to write the ND51b bpch file to disk.  We test the time at the next
!  dynamic timestep so that we can write to disk properly.
!\\
!\\
! !INTERFACE:
!
  FUNCTION ITS_TIME_FOR_WRITE_DIAG51b( Input_Opt, TAU_W ) &
       RESULT( ITS_TIME )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
    USE ERROR_MOD,     ONLY : GEOS_CHEM_STOP
    USE TIME_MOD,      ONLY : GET_HOUR
    USE TIME_MOD,      ONLY : GET_MINUTE
    USE TIME_MOD,      ONLY : GET_SECOND
    USE TIME_MOD,      ONLY : GET_TAU
    USE TIME_MOD,      ONLY : GET_TAUb
    USE TIME_MOD,      ONLY : GET_TAUe
    USE TIME_MOD,      ONLY : GET_TS_DYN
    USE TIME_MOD,      ONLY : GET_TS_DIAG
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: TAU_W   ! TAU at time of disk write
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL  :: ITS_TIME

    REAL(fp) :: TAU, HOUR, DYN, TS_DIAG

    !=================================================================
    ! ITS_TIME_FOR_WRITE_DIAG51b begins here!
    !=================================================================

    ! Initialize
    ITS_TIME = .FALSE.

    ! Add a check for the time to save. Must be a multiple of TS_DIAG
    ! (ccc, 7/21/09)
    TS_DIAG = ( GET_TS_DIAG() / 3600e+0_fp )

    ! NOTE: Change from equality to greater than an a small number,
    ! because 20min or 40min timesteps are irrational numbers in hours
    IF ( MOD(Input_Opt%ND51b_HR_WRITE, TS_DIAG) > 1e-5_fp ) THEN
       WRITE( 6, 100 ) Input_Opt%ND51b_HR_WRITE, TS_DIAG
100    FORMAT( 'The ND51b output frequency must be a multiple ' &
               'of the largest time step:', f9.2, f9.2 )
       CALL GEOS_CHEM_STOP
    ENDIF

    ! Current TAU, Hour, and Dynamic Timestep [hrs]
    TAU      = GET_TAU()
    HOUR     = ( GET_SECOND() / 3600e+0_fp ) + &
               ( GET_MINUTE() / 60e+0_fp ) + GET_HOUR()
    DYN      = ( GET_TS_DYN() / 3600e+0_fp )

    ! If first timestep, return FALSE
    IF ( TAU == GET_TAUb() ) RETURN

    ! If the next dyn timestep is the hour of day
    ! when we have to save to disk, return TRUE
    IF ( MOD( HOUR, 24e+0_fp ) == Input_Opt%ND51b_HR_WRITE ) THEN
       ITS_TIME = .TRUE.
       TAU_W    = TAU + DYN
       RETURN
    ENDIF

    ! If the next dyn timestep is the
    ! end of the run, return TRUE
    IF ( TAU == GET_TAUe() ) THEN
       ITS_TIME = .TRUE.
       TAU_W    = TAU + DYN
       RETURN
    ENDIF

  END FUNCTION ITS_TIME_FOR_WRITE_DIAG51b
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_diag51b
!
! !DESCRIPTION: Subroutine WRITE\_DIAG51b computes the time-average of
!  quantities between local time limits ND51b\_HR1 and ND51b\_HR2 and writes
!  them to a bpch file or HDF5 file.  Arrays and counters are also zeroed
!  for the next diagnostic interval.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE WRITE_DIAG51b( Input_Opt, State_Chm, State_Grid, TAU_W, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE BPCH2_MOD,          ONLY : BPCH2
    USE BPCH2_MOD,          ONLY : OPEN_BPCH2_FOR_WRITE
    USE ERROR_MOD,          ONLY : ALLOC_ERR
    USE inquireMod,         ONLY : findFreeLUN
    USE TIME_MOD,           ONLY : EXPAND_DATE
    USE TIME_MOD,           ONLY : GET_NYMD_DIAG
    USE TIME_MOD,           ONLY : GET_NHMS
    USE TIME_MOD,           ONLY : GET_TAU
    USE TIME_MOD,           ONLY : TIMESTAMP_STRING
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    REAL(fp),       INTENT(IN)    :: TAU_W       ! TAU value at time of write
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
!  Arguments as Input:
!  ============================================================================
!  (1 ) TAU_W (REAL(fp)) : TAU value at time of writing to disk
!
!  NOTES:
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: I,   J,  L,  W, N, GMNL, GMTRC
    INTEGER             :: IOS, X, Y, K, NHMS, nAdvect
    CHARACTER(LEN=16)   :: STAMP
    CHARACTER(LEN=40)   :: CATEGORY
    CHARACTER(LEN=40)   :: UNIT
    CHARACTER(LEN=255)  :: FILENAME

    !=================================================================
    ! WRITE_DIAG51b begins here!
    !=================================================================

    ! Assume success
    RC       =  GC_SUCCESS

    ! Copy values from Input_Opt
    nAdvect  = State_Chm%nAdvect

    ! Find a free file LUN
    IU_ND51b = findFreeLUN()

    ! Replace date tokens in FILENAME
    FILENAME = Input_Opt%ND51b_FILE

    ! Change to get the good timestamp: day that was run and not next
    ! day if saved at midnight
    NHMS = GET_NHMS()
    IF ( NHMS == 0 ) NHMS = 240000

    CALL EXPAND_DATE( FILENAME, GET_NYMD_DIAG(), NHMS )

    ! Echo info
    WRITE( 6, 100 ) TRIM( FILENAME )
100 FORMAT( '     - DIAG51b: Opening file ', a, ' on unit ', i4 )

    ! Open output file
    CALL OPEN_BPCH2_FOR_WRITE( IU_ND51b, FILENAME, TITLE )

    ! Set ENDING TAU for this bpch write
    TAU1 = TAU_W

    !=================================================================
    ! Compute time-average of tracers between local time limits
    !=================================================================

    ! Echo info
    STAMP = TIMESTAMP_STRING()
    WRITE( 6, 110 ) STAMP
110 FORMAT( '     - DIAG51b: Saving to disk at ', a )

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( X, Y, K, W )
    DO W = 1, Input_Opt%N_ND51b

       ! Loop over grid boxes
       DO K = 1, ND51b_NL
       DO Y = 1, ND51b_NJ
       DO X = 1, ND51b_NI

          SELECT CASE( Input_Opt%ND51b_TRACERS(W) )

          CASE( 502 )
             !--------------------------------------------------------
             ! Avoid div by zero for tracers which are archived each
             ! chem timestep and only available in the troposphere
             !--------------------------------------------------------
             IF ( COUNT_CHEM3D(X,Y,K) > 0 ) THEN
                Q(X,Y,K,W) = Q(X,Y,K,W) / COUNT_CHEM3D(X,Y,K)
             ELSE
                Q(X,Y,K,W) = 0e+0_fp
             ENDIF

          CASE DEFAULT

             !--------------------------------------------------------
             ! Avoid division by zero for all other tracers
             !--------------------------------------------------------
             IF ( GOOD_CT(X) > 0 ) THEN
                Q(X,Y,K,W) = Q(X,Y,K,W) / GOOD_CT(X)
             ELSE
                Q(X,Y,K,W) = 0e+0_fp
             ENDIF

          END SELECT

       ENDDO
       ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !=================================================================
    ! Write each tracer from "timeseries.dat" to the timeseries file
    !=================================================================
    DO W = 1, Input_Opt%N_ND51b

       ! ND51b tracer number
       N = Input_Opt%ND51b_TRACERS(W)

       ! Save by simulation
       IF ( N <= nAdvect ) THEN

          !---------------------
          ! GEOS-CHEM species
          !---------------------
          CATEGORY = 'IJ-AVG-$'
          UNIT     = ''              ! Let GAMAP pick unit
          GMNL     = ND51b_NL
          GMTRC    = N

       ELSE IF ( N == 501 ) THEN

          !---------------------
          ! OH
          !---------------------
          CATEGORY  = 'CHEM-L=$'
          UNIT      = 'molec/cm3'
          GMNL      = ND51b_NL
          GMTRC     = 1

       ELSE IF ( N == 502 ) THEN

          !---------------------
          ! NOy
          !---------------------
          CATEGORY = 'TIME-SER'
          UNIT     = ''              ! Let GAMAP pick unit
          GMNL     = ND51b_NL
          GMTRC    = 2

       ELSE IF ( N == 503 ) THEN

          !---------------------
          ! Relative humidity
          !---------------------
          CATEGORY = 'TIME-SER'
          UNIT     = '%'
          GMNL     = ND51b_NL
          GMTRC    = 3

       ELSE IF ( N == 504 ) THEN

          !---------------------
          ! 3-D Cloud fractions
          !---------------------
          CATEGORY = 'TIME-SER'
          UNIT     = 'unitless'
          GMNL     = ND51b_NL
          GMTRC    = 4

       ELSE IF ( N == 505 ) THEN

          !---------------------
          ! Column opt depths
          !---------------------
          CATEGORY = 'TIME-SER'
          UNIT     = 'unitless'
          GMNL     = 1
          GMTRC    = 5

       ELSE IF ( N == 506 ) THEN

          !---------------------
          ! Cloud top heights
          !---------------------
          CATEGORY = 'TIME-SER'
          UNIT     = 'hPa'
          GMNL     = 1
          GMTRC    = 6

       ELSE IF ( N == 507 ) THEN

          !---------------------
          ! Air Density
          !---------------------
          CATEGORY = 'TIME-SER'
          UNIT     = 'molec/cm3'
          GMNL     = ND51b_NL
          GMTRC    = 7

       ELSE IF ( N == 508 ) THEN

          !---------------------
          ! Total seasalt
          !---------------------
          CATEGORY = 'TIME-SER'
          UNIT     = ''              ! Let GAMAP pick unit
          GMNL     = ND51b_NL
          GMTRC    = 8

       ELSE IF ( N == 509 ) THEN

          !---------------------
          ! PBL Height [m]
          !---------------------
          CATEGORY = 'PBLDEPTH'
          UNIT     = 'm'
          GMNL     = 1
          GMTRC    = 1

       ELSE IF ( N == 510 ) THEN

          !---------------------
          ! PBL Height [levels]
          !---------------------
          CATEGORY = 'PBLDEPTH'
          UNIT     = 'levels'
          GMNL     = 1
          GMTRC    = 2

       ELSE IF ( N == 511 ) THEN

          !---------------------
          ! Grid box heights
          !---------------------
          CATEGORY = 'BXHGHT-$'
          UNIT     = 'm'
          GMNL     = ND51b_NL
          GMTRC    = 1

       ELSE IF ( N == 512 ) THEN

          !---------------------
          ! PEDGE-$
          !---------------------
          CATEGORY = 'PEDGE-$'
          UNIT     = 'hPa'
          GMNL     = ND51b_NL
          GMTRC    = 1

       ELSE IF ( N == 513 ) THEN

          !---------------------
          ! Sea level prs
          !---------------------
          CATEGORY = 'DAO-FLDS'
          UNIT     = 'hPa'
          GMNL     = 1
          GMTRC    = 18

       ELSE IF ( N == 514 ) THEN

          !---------------------
          ! U-wind
          !---------------------
          CATEGORY = 'DAO-3D-$'
          UNIT     = 'm/s'
          GMNL     = ND51b_NL
          GMTRC    = 1

       ELSE IF ( N == 515 ) THEN

          !---------------------
          ! V-wind
          !---------------------
          CATEGORY = 'DAO-3D-$'
          UNIT     = 'm/s'
          GMNL     = ND51b_NL
          GMTRC    = 2

       ELSE IF ( N == 516 ) THEN

          !---------------------
          ! Temperature
          !---------------------
          CATEGORY = 'DAO-3D-$'
          UNIT     = 'K'
          GMNL     = ND51b_NL
          GMTRC    = 3


       ELSE IF ( N == 517 ) THEN

          !---------------------
          ! Sulfate AOD
          !---------------------
          CATEGORY = 'OD-MAP-$'
          UNIT     = 'unitless'
          GMNL     = ND51b_NL
          GMTRC    = 6

       ELSE IF ( N == 518 ) THEN

          !---------------------
          ! Black Carbon AOD
          !---------------------
          CATEGORY = 'OD-MAP-$'
          UNIT     = 'unitless'
          GMNL     = ND51b_NL
          GMTRC    = 9

       ELSE IF ( N == 519 ) THEN

          !---------------------
          ! Organic Carbon AOD
          !---------------------
          CATEGORY = 'OD-MAP-$'
          UNIT     = 'unitless'
          GMNL     = ND51b_NL
          GMTRC    = 12

       ELSE IF ( N == 520 ) THEN

          !---------------------
          ! SS Accum AOD
          !---------------------
          CATEGORY = 'OD-MAP-$'
          UNIT     = 'unitless'
          GMNL     = ND51b_NL
          GMTRC    = 15

       ELSE IF ( N == 521 ) THEN

          !---------------------
          ! SS Coarse AOD
          !---------------------
          CATEGORY = 'OD-MAP-$'
          UNIT     = 'unitless'
          GMNL     = ND51b_NL
          GMTRC    = 18

       ELSE IF ( N == 522 ) THEN

          !---------------------
          ! Total dust OD
          !---------------------
          CATEGORY = 'OD-MAP-$'
          UNIT     = 'unitless'
          GMNL     = ND51b_NL
          GMTRC    = 4

       ELSE IF ( ( N >= 523 ) .and. ( N <= 529) ) THEN

          !---------------------
          ! dust OD (bins 1-7)
          !---------------------
          CATEGORY = 'OD-MAP-$'
          UNIT     = 'unitless'
          GMNL     = ND51b_NL
          GMTRC    = 21+(N-523)

       ELSE IF ( N == 530 ) THEN

          !---------------------
          ! PARDR [W/m2]
          ! (mpb,2009)
          !---------------------
          CATEGORY = 'DAO-FLDS'
          UNIT     = 'W/m2'
          GMNL     = ND51b_NL
          GMTRC    = 20

       ELSE IF ( N == 531 ) THEN

          !---------------------
          ! PARDF [W/m2]
          ! (mpb,2009)
          !---------------------
          CATEGORY = 'DAO-FLDS'
          UNIT     = 'W/m2'
          GMNL     = ND51b_NL
          GMTRC    = 21

       ELSE IF ( N == 532 ) THEN

          !---------------------
          ! DAILY LAI [W/m2]
          ! (mpb,2009)
          !---------------------
          CATEGORY = 'TIME-SER'
          UNIT     = 'm2/m2'
          GMNL     = ND51b_NL
          GMTRC    = 9

       ELSE IF ( N == 533 ) THEN

          !---------------------
          ! T at 2m
          ! (mpb,2008)
          !---------------------
          CATEGORY = 'DAO-FLDS'
          UNIT     = 'K'
          GMNL     = ND51b_NL
          GMTRC    = 5

       ELSE

          ! Otherwise skip
          CYCLE

       ENDIF

       !------------------------
       ! Save to bpch file
       !------------------------
       CALL BPCH2( IU_ND51b,     MODELNAME,    LONRES, &
                   LATRES,       HALFPOLAR,    CENTER180, &
                   CATEGORY,     GMTRC,        UNIT, &
                   TAU0,         TAU1,         RESERVED, &
                   ND51b_NI,     ND51b_NJ,     GMNL, &
                   Input_Opt%ND51b_IMIN+I0, &
                   Input_Opt%ND51b_JMIN+J0, &
                   Input_Opt%ND51b_LMIN, &
                   REAL( Q(1:ND51b_NI, 1:ND51b_NJ, 1:GMNL, W),4))

    ENDDO

    ! Echo info
    WRITE( 6, 120 ) TRIM( FILENAME )
120 FORMAT( '     - DIAG51b: Closing file ', a )

    ! Close file
    CLOSE( IU_ND51b )

    !=================================================================
    ! Re-initialize quantities for next diagnostic cycle
    !=================================================================

    ! Echo info
    STAMP = TIMESTAMP_STRING()
    WRITE( 6, 130 ) STAMP
130 FORMAT( '     - DIAG51b: Zeroing arrays at ', a )

    ! Set STARTING TAU for the next bpch write
    TAU0 = TAU_W

    ! Zero accumulating array for tracer
    Q            = 0e+0_fp

    ! Zero counter arrays
    COUNT_CHEM3D = 0e+0_fp
    GOOD_CT      = 0e+0_fp

  END SUBROUTINE WRITE_DIAG51b
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_i
!
! !DESCRIPTION: Function GET\_I returns the absolute longitude index (I),
!  given the relative longitude index (X).
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_I( X, State_Grid ) RESULT( I )
!
! !USES:
!
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN) :: X          ! Relative longitude index
    TYPE(GrdState), INTENT(IN) :: State_Grid ! Grid State object
!
! !RETURN VALUE:
!
    INTEGER             :: I   ! Absolute longitude index
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! GET_I begins here!
    !=================================================================

    ! Add the offset to X to get I
    I = IOFF + X

    ! Handle wrapping around the date line, if necessary
    IF ( I > State_Grid%NX ) I = I - State_Grid%NX

  END FUNCTION GET_I
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_diag51b
!
! !DESCRIPTION: Subroutine INIT\_DIAG51b allocates and zeroes all module
!  arrays.  It also gets values for module variables from "input\_mod.F90".
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_DIAG51b( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
    USE BPCH2_MOD,      ONLY : GET_MODELNAME
    USE BPCH2_MOD,      ONLY : GET_HALFPOLAR
    USE ErrCode_Mod
    USE ERROR_MOD,      ONLY : ALLOC_ERR
    USE ERROR_MOD,      ONLY : ERROR_STOP
    USE TIME_MOD,       ONLY : GET_TAUb
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: AS
    CHARACTER(LEN=255) :: LOCATION

    !=================================================================
    ! INIT_DIAG51b begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Exit if ND51b is turned off, or if it's a dry-run
    IF ( ( .not. Input_Opt%DO_ND51b ) .or. Input_Opt%DryRun ) RETURN

    ! Initialize
    LOCATION = 'INIT_DIAG51b ("diag51b_mod.F90")'

    !=================================================================
    ! Error check longitude, latitude, altitude limits
    !=================================================================

    ! Get grid offsets
    I0 = State_Grid%XMinOffset
    J0 = State_Grid%YMinOffset

    !-----------
    ! Longitude
    !-----------

    ! Error check ND51b_IMIN
    IF ( Input_Opt%ND51b_IMIN < 1 .or. &
         Input_Opt%ND51b_IMIN > State_Grid%NX ) THEN
       CALL ERROR_STOP( 'Bad ND51b_IMIN value!', LOCATION )
    ENDIF

    ! Error check ND51b_IMAX
    IF ( Input_Opt%ND51b_IMAX < 1 .or. &
         Input_Opt%ND51b_IMAX > State_Grid%NX ) THEN
       CALL ERROR_STOP( 'Bad ND51b_IMAX value!', LOCATION )
    ENDIF

    ! Compute longitude limits to write to disk
    ! Also handle wrapping around the date line
    IF ( Input_Opt%ND51b_IMAX >= Input_Opt%ND51b_IMIN ) THEN
       ND51b_NI = ( Input_Opt%ND51b_IMAX - Input_Opt%ND51b_IMIN ) + 1
    ELSE
       ND51b_NI = ( State_Grid%NX - Input_Opt%ND51b_IMIN ) + 1 + &
                    Input_Opt%ND51b_IMAX
       WRITE( 6, '(a)' ) 'We are wrapping over the date line!'
    ENDIF

    ! Make sure that ND50_NI <= IIPAR
    IF ( ND51b_NI > State_Grid%NX ) THEN
       CALL ERROR_STOP( 'Too many longitudes!', LOCATION )
    ENDIF

    !-----------
    ! Latitude
    !-----------

    ! Error check JMIN_AREA
    IF ( Input_Opt%ND51b_JMIN < 1 .or. &
         Input_Opt%ND51b_JMIN > State_Grid%NY ) THEN
       CALL ERROR_STOP( 'Bad ND51b_JMIN value!', LOCATION )
    ENDIF

    ! Error check JMAX_AREA
    IF ( Input_Opt%ND51b_JMAX < 1 .or. &
         Input_Opt%ND51b_JMAX > State_Grid%NY ) THEN
       CALL ERROR_STOP( 'Bad ND51b_JMAX value!', LOCATION )
    ENDIF

    ! Compute latitude limits to write to disk (bey, bmy, 3/16/99)
    IF ( Input_Opt%ND51b_JMAX >= Input_Opt%ND51b_JMIN ) THEN
       ND51b_NJ = ( Input_Opt%ND51b_JMAX - Input_Opt%ND51b_JMIN ) + 1
    ELSE
       CALL ERROR_STOP( 'ND51b_JMAX < ND51b_JMIN!', LOCATION )
    ENDIF

    !-----------
    ! Altitude
    !-----------

    ! Error check ND51b_LMIN, ND51b_LMAX
    IF ( Input_Opt%ND51b_LMIN < 1 .or. &
         Input_Opt%ND51b_LMAX > State_Grid%NZ ) THEN
       CALL ERROR_STOP( 'Bad ND51b altitude values!', LOCATION )
    ENDIF

    ! # of levels to save in ND51b timeseries
    IF ( Input_Opt%ND51b_LMAX >= Input_Opt%ND51b_LMIN ) THEN
       ND51b_NL = ( Input_Opt%ND51b_LMAX - Input_Opt%ND51b_LMIN ) + 1
    ELSE
       CALL ERROR_STOP( 'ND51b_LMAX < ND51b_LMIN!', LOCATION )
    ENDIF

    !-----------
    ! Offsets
    !-----------
    IOFF      = Input_Opt%ND51b_IMIN - 1
    JOFF      = Input_Opt%ND51b_JMIN - 1
    LOFF      = Input_Opt%ND51b_LMIN - 1

    !-----------
    ! For bpch
    !-----------
    TAU0      = GET_TAUb()
    TITLE     = 'GEOS-CHEM DIAG51b time series'
    LONRES    = State_Grid%DX
    LATRES    = State_Grid%DY
    MODELNAME = GET_MODELNAME( Input_Opt, State_Grid )
    HALFPOLAR = GET_HALFPOLAR()

    !=================================================================
    ! Allocate arrays
    !=================================================================

    ! Array denoting where LT is between HR1 and HR2
    ALLOCATE( GOOD( State_Grid%NX ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'GOOD' )
    GOOD = 0

    ! Counter of "good" times per day at each grid box
    ALLOCATE( GOOD_CT( ND51b_NI ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'GOOD_CT' )
    GOOD_CT = 0

    ! Accumulating array
    ALLOCATE( Q( ND51b_NI, ND51b_NJ, ND51b_NL, Input_Opt%N_ND51b ), STAT=AS)
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'Q' )
    Q = 0e+0_fp

    ! Accumulating array
    ALLOCATE( COUNT_CHEM3D( ND51b_NI, ND51b_NJ, ND51b_NL ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'COUNT_CHEM3D' )
    COUNT_CHEM3D = 0

  END SUBROUTINE INIT_DIAG51b
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_diag51b
!
! !DESCRIPTION: Subroutine CLEANUP\_DIAG51b deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_DIAG51b()
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! CLEANUP_DIAG51b begins here!
    !=================================================================
    IF ( ALLOCATED( COUNT_CHEM3D ) ) DEALLOCATE( COUNT_CHEM3D )
    IF ( ALLOCATED( GOOD         ) ) DEALLOCATE( GOOD         )
    IF ( ALLOCATED( GOOD_CT      ) ) DEALLOCATE( GOOD_CT      )
    IF ( ALLOCATED( Q            ) ) DEALLOCATE( Q            )
    
  END SUBROUTINE CLEANUP_DIAG51b
!EOC
END MODULE DIAG51b_MOD
#endif
