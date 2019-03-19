!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_seasalt_mod.F90
!
! !DESCRIPTION: Module HCOX\_SeaSalt\_Mod contains routines to calculate
! sea salt aerosol emissions, following the implementation in GEOS-Chem.
! Emission number densities of the fine and coarse mode sea salt aerosols
! are written into diagnostic containers `SEASALT\_DENS\_FINE` and
! `SEASALT\_DENS\_COARSE`, respectively.
!\\
!\\ 
! This is a HEMCO extension module that uses many of the HEMCO core
! utilities.
!\\
!\\
! !INTERFACE: 
!
MODULE HCOX_SeaSalt_Mod
!
! !USES:
! 
  USE HCO_Error_Mod
  USE HCO_Diagn_Mod
  USE HCO_State_Mod,  ONLY : HCO_State
  USE HCOX_State_Mod, ONLY : Ext_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOX_SeaSalt_Init
  PUBLIC :: HCOX_SeaSalt_Run
  PUBLIC :: HCOX_SeaSalt_Final
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: EMIT_SSABr2
!
! !REVISION HISTORY:
!  15 Dec 2013 - C. Keller   - Now a HEMCO extension module 
!  09 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  09 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  09 Jul 2015 - E. Lundgren - Add marine organoc aerosols (B.Gantt, M.Johnson)
!EOP
!------------------------------------------------------------------------------
!
! !PRIVATE TYPES:
!
  TYPE :: MyInst
   ! Tracer IDs 
   INTEGER                :: Instance
   INTEGER                :: ExtNr

   ! Tracer IDs 
   INTEGER             :: ExtNrSS           ! Extension number for seasalt
   INTEGER             :: ExtNrMPOA         ! Extension number for marine POA
   INTEGER             :: IDTSALA           ! Fine aerosol model species ID
   INTEGER             :: IDTSALC           ! Coarse aerosol model species ID  
   INTEGER             :: IDTMOPO           ! marine organic aerosol - phobic
   INTEGER             :: IDTMOPI           ! marine organic aerosol - philic
   INTEGER             :: IDTBr2            ! Br2 model species ID
   INTEGER             :: IDTBrSALA         ! Br- in accum. sea salt aerosol
   INTEGER             :: IDTBrSALC         ! Br- in coarse sea salt aerosol
   LOGICAL             :: CalcBr2           ! Calculate Br2 SSA emissions?
   LOGICAL             :: CalcBrSalt        ! Calculate Br- content?

   ! Scale factors
   REAL*8              :: Br2Scale          ! Br2 scale factor 
   REAL*8              :: BrContent         ! Ratio of Br- to dry SSA (mass)
   REAL*8              :: WindScale         ! Wind adjustment factor

   ! Module variables
   INTEGER              :: NSALT             ! # of seasalt tracers
   INTEGER, POINTER     :: NR(:)             ! Size bin information
   REAL*8,  POINTER     :: SRRC  (:,:)
   REAL*8,  POINTER     :: SRRC_N(:,:)
   REAL*8,  POINTER     :: RREDGE(:,:)
   REAL*8,  POINTER     :: RRMID (:,:)
   REAL*8,  POINTER     :: SS_DEN(:)         ! densities 

   ! Number densities
   REAL(sp), POINTER   :: NDENS_SALA(:,:) => NULL()
   REAL(sp), POINTER   :: NDENS_SALC(:,:) => NULL()
   REAL(sp), POINTER   :: NDENS_MOPO(:,:) => NULL() 
   REAL(sp), POINTER   :: NDENS_MOPI(:,:) => NULL() 

   TYPE(MyInst), POINTER  :: NextInst => NULL()
  END TYPE MyInst

  ! Pointer to instances
  TYPE(MyInst), POINTER   :: AllInst => NULL()
!
! !DEFINED PARAMETERS:
!
  INTEGER, PARAMETER  :: NR_MAX = 200  ! max. # of bins per mode

  ! Increment of radius for Emission integration (um)
  REAL*8, PARAMETER   :: DR    = 5.d-2
  REAL*8, PARAMETER   :: BETHA = 2.d0

  ! this needs adding to input file switches (Offline seasalt)
  LOGICAL, PARAMETER  :: OFFSEASALT = .TRUE.

CONTAINS
!EOC
!-------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_SeaSalt_Run 
!
! !DESCRIPTION: Subroutine HcoX\_SeaSalt\_Run is the driver run routine to 
! calculate SeaSalt emissions in HEMCO.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_SeaSalt_Run( am_I_Root, ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_FluxArr_Mod,      ONLY : HCO_EmisAdd
    USE HCO_CLOCK_MOD,        ONLY : HcoClock_Get
    USE HCO_GeoTools_Mod,     ONLY : HCO_LANDTYPE
    USE Input_Opt_Mod,        ONLY : OptInput
    USE HCO_INTERP_MOD,       ONLY : REGRID_MAPA2A
    USE CMN_SIZE_MOD                        ! Size parameters
    USE TIME_MOD                            ! Date & time routines
    ! NcdfUtil modules for netCDF I/O
    USE m_netcdf_io_open                    ! netCDF open
    USE m_netcdf_io_get_dimlen              ! netCDF dimension queries
    USE m_netcdf_io_read                    ! netCDF data reads
    USE m_netcdf_io_close                   ! netCDF close

    TYPE(OptInput), POINTER    :: Input_Opt  ! Input Options object      
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   ) :: am_I_Root  ! root CPU?
    TYPE(HCO_State), POINTER       :: HcoState   ! Output obj
    TYPE(Ext_State), POINTER       :: ExtState  ! Module options  
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC         ! Success or failure?
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Chin, M., P. Ginoux, S. Kinne, B. Holben, B. Duncan, R. Martin,
!        J. Logan, A. Higurashi, and T. Nakajima, "Tropospheric aerosol
!        optical thickness from the GOCART model and comparisons with
!        satellite and sunphotometers measurements", J. Atmos Sci., 2001.
!  (2 ) Gong, S., L. Barrie, and J.-P. Blanchet, "Modeling sea-salt
!        aerosols in the atmosphere. 1. Model development", J. Geophys. Res.,
!        v. 102, 3805-3818, 1997.
!  (3 ) Gong, S. L., "A parameterization of sea-salt aerosol source function
!        for sub- and super-micron particles", Global Biogeochem.  Cy., 17(4), 
!        1097, doi:10.1029/2003GB002079, 2003.
!  (4 ) Jaegle, L., P.K. Quinn, T.S. Bates, B. Alexander, J.-T. Lin, "Global
!        distribution of sea salt aerosols: New constraints from in situ and 
!        remote sensing observations", Atmos. Chem. Phys., 11, 3137-3157, 
!        doi:10.5194/acp-11-3137-2011.
!
! !REVISION HISTORY: 
!  (1 ) Now references SALA_RREDGE_um and SALC_RREDGE_um from "tracer_mod.f"
!        (bmy, 7/20/04)
!  (2 ) Now references GET_FRAC_OF_PBL and GET_PBL_TOP_L from "pbl_mix_mod.f".
!        Removed reference to header file CMN.  Removed reference to 
!        "pressure_mod.f".  (bmy, 2/22/05)
!  (3 ) Now also compute alkalinity and number density of SeaSalt emissions.
!        (bec, bmy, 4/13/05)
!  (4 ) Now references XNUMOL & XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (5 ) The source function is for wet aerosol radius (RH=80%, with a radius
!        twice the size of dry aerosols) so BETHA should be set to 2 
!        instead of 1.  Also now use LOG10 instead of LOG in the expressions
!        for the SeaSalt base source, since we need the logarithm to the base
!        10. (jaegle, bec, bmy, 11/23/09)
!  (6 ) Update to use the Gong (2003) source function (jaegle 5/11/11)
!  (7 ) Apply an empirical sea surface temperature dependence to Gong (2003)
!       (jaegle 5/11/11)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  15 Dec 2013 - C. Keller   - Now a HEMCO extension 
!  09 Jul 2015 - E. Lundgren - Add marine organic aerosols (B.Gantt, M.Johnson)
!  19 Oct 2015 - C. Keller   - Now pass I and J index to EMIT_SSABr2 to support
!                              curvilinear grids.
!  22 Oct 2015 - E. Lundgren - Bug fix: include CHLR in OMP PRIVATE statement
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(MyInst), POINTER  :: Inst
    INTEGER                :: I, J, N, R
    REAL*8                 :: SALT, SALT_N, SSA_BR2, CHLR
    REAL*8                 :: A_M2
    REAL*8                 :: W10M
    REAL                   :: FLUX
    REAL(hp), TARGET       :: FLUXSALA  (HcoState%NX,HcoState%NY)
    REAL(hp), TARGET       :: FLUXSALC  (HcoState%NX,HcoState%NY)
    REAL(hp), TARGET       :: FLUXBr2   (HcoState%NX,HcoState%NY)
    REAL(hp), TARGET       :: FLUXBrSalA(HcoState%NX,HcoState%NY)
    REAL(hp), TARGET       :: FLUXBrSalC(HcoState%NX,HcoState%NY)
    REAL(hp), TARGET       :: FLUXMOPO  (HcoState%NX,HcoState%NY)
    REAL(hp), TARGET       :: FLUXMOPI  (HcoState%NX,HcoState%NY)

    ! New variables (jaegle 5/11/11)
    REAL*8                 :: SST, SCALE
    ! jpp, 3/2/10
    REAL*8                 :: BR2_NR, SALT_NR 
    ! B. Gantt, M. Johnson (7,9/15)
    REAL*8                 :: OMSS1, OMSS2

    ! Error handling
    LOGICAL                :: ERR
    CHARACTER(LEN=255)     :: MSG
    
    ! OFFLINE SEASALT PARAMETERS
    INTEGER                :: fA1          ! netCDF file ID
    INTEGER                :: X, Y         ! netCDF file dimensions
    INTEGER                :: YYYYMMDD, HHMMSS, cYr,cMt,cDy,cHr
    INTEGER                :: II, JJ, TIME_INDEX
    CHARACTER(LEN=255)     :: dir,nc_file,v_name
    CHARACTER(LEN=4)       :: Yrs, Mts, Dys, hrs
    CHARACTER(LEN=2)       :: STRBIN
!    INTEGER                :: st1d(1), ct1d(1)    ! Start + count, for 1D arrays
    INTEGER                :: st3d(3), ct3d(3)    ! Start + count, for 3D arrays
    REAL*4                 :: Q(IIPAR,JJPAR)     ! Temporary data arrray
    !=================================================================
    ! HCOX_SeaSalt_Run begins here!
    !=================================================================

    ! Return if extension disabled 
    IF ( ExtState%SeaSalt <= 0 ) RETURN

    ! Enter 
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_SeaSalt_Run (hcox_seasalt_mod.F90)', RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Exit status
    ERR = .FALSE.

    ! Get instance
    Inst   => NULL()
    CALL InstGet ( ExtState%SeaSalt, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN 
       WRITE(MSG,*) 'Cannot find SeaSalt instance Nr. ', ExtState%SeaSalt
       CALL HCO_ERROR(HcoState%Config%Err,MSG,RC)
       RETURN
    ENDIF

    ! Init values
    FLUXSALA   = 0.0_hp
    FLUXSALC   = 0.0_hp
    FLUXBr2    = 0.0_hp
    FLUXBrSalA = 0.0_hp
    FLUXBrSalC = 0.0_hp
    FLUXMOPO   = 0.0_hp
    FLUXMOPI   = 0.0_hp

   ! We have a flag for offline seasalt
   IF (OFFSEASALT .EQ. .FALSE.) THEN      

    !=================================================================
    ! Emission is integrated over a given size range for each bin
    !=================================================================
!$OMP PARALLEL DO                                                      &
!$OMP DEFAULT( SHARED )                                                &
!$OMP PRIVATE( I, J, A_M2, W10M, SST, SCALE, SSA_BR2, N              ) &
!$OMP PRIVATE( SALT, SALT_N, R, SALT_NR, BR2_NR, RC                  ) & 
!$OMP PRIVATE( OMSS1, OMSS2, CHLR                                    ) & 
!$OMP SCHEDULE( DYNAMIC )

    ! Loop over surface boxes 
    DO J = 1, HcoState%NY 
    DO I = 1, HcoState%NX

       ! Grid box surface area on simulation grid [m2]
       A_M2 = HcoState%Grid%AREA_M2%Val( I, J )

       ! Advance to next grid box if it's not over water
       IF ( HCO_LANDTYPE( ExtState%WLI%Arr%Val(I,J), &
                          ExtState%ALBD%Arr%Val(I,J) ) /= 0 ) CYCLE

       ! Wind speed at 10 m altitude [m/s]
       W10M = SQRT( ExtState%U10M%Arr%Val(I,J)**2 &
                  + ExtState%V10M%Arr%Val(I,J)**2 ) 

       ! Sea surface temperature in Celcius (jaegle 5/11/11)
       SST = ExtState%TSKIN%Arr%Val(I,J) - 273.15d0

       ! Limit SST to 0-30C range
       SST = MAX( SST , 0d0 )  ! limit to  0C
       SST = MIN( SST , 30d0 ) ! limit to 30C

       ! Empirical SST scaling factor (jaegle 5/11/11)
       SCALE = 0.329d0 + 0.0904d0*SST - &
               0.00717d0*SST**2d0 + 0.000207d0*SST**3d0

       ! Reset to using original Gong (2003) emissions (jaegle 6/30/11)
       !SCALE = 1.0d0

       ! Eventually apply wind scaling factor. 
       SCALE = SCALE * Inst%WindScale

       ! Reset sea salt aerosol emissions of Br2, which are integrated
       ! over accumulation and coarse mode
       SSA_BR2 = 0d0
    
       ! Do for accumulation and coarse mode, and Marine POA if enabled
       DO N = 1,Inst%NSALT  

          ! Reset values for SALT and SALT_N
          SALT   = 0d0
          SALT_N = 0d0

          ! Loop over size bins
          DO R = 1, Inst%NR(N)

             ! Coarse and accumulation modes
             IF ( N .LT. 3 ) THEN 
             
                ! Update SeaSalt source into SALT [kg]
                SALT   = SALT +                                   &
                         ( SCALE * Inst%SRRC(R,N) * A_M2 * W10M**3.41d0 )
   
                ! Update SeaSalt source into SALT_N [#] 
                ! (bec, bmy, 4/13/05)
                SALT_N = SALT_N +                               &
                         ( SCALE * Inst%SRRC_N(R,N) * A_M2 * W10M**3.41d0 )
   
                ! --------------------------------------------------
                ! jpp, 3/2/10: Accounting for the bromine emissions
                ! now. Store mass flux [kg] of Br2 based on how much
                ! aerosol there is emitted in this box.
                ! --------------------------------------------------
                IF ( Inst%CalcBr2 ) THEN
   
                   ! jpp, 3/3/10: since the SALT arrays are integrations
                   !              I cannot use them for each independent
                   !              radius... that's why I'm getting too
                   !              much bromine. So I'm making a tmp
                   !              array to store only the current
                   !              Dry Radius bin. [kg]
                   SALT_NR = ( Inst%SRRC(R,N) * A_M2 * W10M**3.41d0 )
                   CALL EMIT_SSABr2( am_I_Root, ExtState,  HcoState, Inst, &
                                     I, J,      Inst%RRMID(R,N), SALT_NR,  &
                                     BR2_NR,    RC                    )
                   IF ( RC /= HCO_SUCCESS ) THEN
                      ERR = .TRUE.
                      EXIT
                   ENDIF
                   SSA_Br2 = SSA_Br2 + BR2_NR
                ENDIF

             ENDIF 

             ! Marine organic aerosols (M. Johnson, B. Gantt)
             IF ( N .EQ. 3 ) THEN 

                ! Get MODIS Chlorophyll-a
                CHLR = ExtState%CHLR%Arr%Val(I,J)

                ! Calculate organic mass fraction of SSA
                OMSS1 = 1.0 / ( 1.0 + EXP( -2.63 * 3.0 * CHLR         &
                        + 0.18 * 3.0 * W10M ) )

                OMSS2 = ( OMSS1 ) / (1.0 + 0.03                       &
                        * EXP( 6.81 * ( Inst%RRMID(R,N) * 2.0 ) ) )        &
                        + 0.03 * ( OMSS1 )

                ! Update seasalt source into SALT [kg]
                SALT  = SALT + 6.0 * ( ( Inst%SRRC(R,N) * SCALE * A_M2     &
                               * W10M**3.41d0 * OMSS2 )               &
                               * ( 1.0 / ( 2.2 / ( 1.0 - OMSS2        &
                               * (1.0 - 2200.0 / 1000.0 ) ) ) ) )

                SALT_N = SALT_N +  6.0 * ( Inst%SRRC_N(R,N) * SCALE * A_M2 &
                                   * W10M**3.41d0 * OMSS2 )

             ENDIF

          ENDDO !R

          ! ----------------------------------------------------------------
          ! Pass sea salt emissions do emission array [kg/m2/s]
          ! ----------------------------------------------------------------
          ! kg --> kg/m2/s
          IF     ( N == 1 ) THEN
             FLUXSALA(I,J) = SALT / A_M2 / HcoState%TS_EMIS 
          ELSEIF ( N == 2 ) THEN    
             FLUXSALC(I,J) = SALT / A_M2 / HcoState%TS_EMIS 
          ELSEIF ( N == 3 ) THEN    
             FLUXMOPO(I,J) = SALT / A_M2 / HcoState%TS_EMIS 
          ELSEIF ( N == 4 ) THEN    
             FLUXMOPI(I,J) = SALT / A_M2 / HcoState%TS_EMIS 
          ENDIF

          ! ----------------------------------------------------------------
          ! Write out number density for diagnostics [#]
          ! ----------------------------------------------------------------
          IF     ( N == 1 ) THEN
             Inst%NDENS_SALA(I,J) = SALT_N 
          ELSEIF ( N == 2 ) THEN
             Inst%NDENS_SALC(I,J) = SALT_N 
          ELSEIF ( N == 3 ) THEN
             Inst%NDENS_MOPO(I,J) = SALT_N 
          ELSEIF ( N == 4 ) THEN
             Inst%NDENS_MOPI(I,J) = SALT_N 
          ENDIF

!=============================================================================
!            ! Alkalinity [kg] (bec, bmy, 4/13/05)
!            ALK_EMIS(I,J,L,N) = SALT
!
!            ! Number density [#/m3] (bec, bmy, 4/13/05)
!            N_DENS(I,J,L,N)   = SALT_N / State_Met%AIRVOL(I,J,L)

!            ! ND08 diagnostic: sea salt emissions [kg]
!            IF ( ND08 > 0 ) THEN
!               AD08(I,J,N) = AD08(I,J,N) + SALT
!            ENDIF
!=============================================================================

       ENDDO !N

       ! Store Br2 flux in tendency array in [kg/m2/s]
       IF ( Inst%CalcBr2 ) THEN 

          ! kg --> kg/m2/s
          FLUXBR2(I,J) = SSA_Br2 / A_M2 / HcoState%TS_EMIS
       ENDIF

!=============================================================================
!          ! Add to diagnostics
!          ! ND08 diagnostic: sea salt emissions [kg]
!          IF ( ND08 > 0 ) THEN
!             AD08(I,J,N) = AD08(I,J,N) + SALT(I,J)
!          ENDIF
!
!          IF ( ND46 > 0 ) THEN
!             ! store the emission in the AD46 Biogenic Emissions
!             ! diagnostic array [kg/m2/s]
!             AD46(I,J,16) = AD46(I,J,16) + ( SSA_Br2 / A_M2 / DTEMIS )
!          ENDIF
!=============================================================================

    ENDDO !I
    ENDDO !J
!$OMP END PARALLEL DO

    ! Check exit status 
    IF ( ERR ) THEN
       RC = HCO_FAIL
       RETURN 
    ENDIF
    

   ELSE 
     !Read the offline emissions for the current hour
     !Format to FLUX to match what the following code expects
     ! Get emissions time step
     CALL HcoClock_Get( am_I_Root, HcoState%Clock, cYYYY=cYr, cMM=cMt, cDD=cDy, cH=cHr, RC=RC )
     YYYYMMDD  = cYr*10000 + cMt*100 + cDy
     HHMMSS    = cHr*10000

     ! Write datetime
     WRITE( Yrs, '(i4.4)' ) cYr
     WRITE( Mts, '(i2.2)' ) cMt
     WRITE( Dys, '(i2.2)' ) cDy
     WRITE( hrs, '(i2.2)' ) cHr
     
     ! replace time & date tokens in the file name
     dir     = 'offline_seasalt/' // 'GEOS_FP.2x2.5/'// trim(Yrs) // '/' // trim(Mts) // '/'

     ! Replace time & date tokens in the file name
     nc_file = 'GEOSFP.2x25.' // trim(Yrs) // trim(Mts) &
              // trim(Dys) // trim(hrs) // '00.nc'

     ! Construct complete file path
     nc_file = TRIM( '/as/data/geos/ExtData/HEMCO/OFFLINE_NATURAL_EMISSIONS/' ) //  &
               TRIM( dir ) // TRIM( nc_file )
     ! Open netCDF file
     CALL NcOp_Rd( fA1, TRIM( nc_file ) )

     ! Read the dimensions from the netCDF file
     CALL NcGet_DimLen( fA1, 'lon',   X )
     CALL NcGet_DimLen( fA1, 'lat',   Y )
     !======================================================================
     ! Read data from the netCDF file
     !======================================================================
     ! Find the proper time-slice to read from disk
     time_index = ( HHMMSS / 10000 ) + 1

     ! netCDF start & count indices
     st3d      = (/ 1,     1,     1          /)
     ct3d      = (/ IIPAR, JJPAR, 1          /)

     ! Read emissions

      CALL NcRd( Q, fA1, TRIM("SALA_TOTAL"), st3d, ct3d )
      !======================================================================
      ! Regrid the NC data if necessary (might be high res)
      !======================================================================
!      CALL REGRID_MAPA2A ( am_I_Root, HcoState, Q, LonE, LatE, CLct, RC )
      DO II=1,IIPAR
      DO JJ=1,JJPAR
       FLUXSALA(II,JJ) = Q(II,JJ) 
      ENDDO
      ENDDO 

      CALL NcRd( Q, fA1, TRIM("SALC_TOTAL"), st3d, ct3d )
!      CALL REGRID_MAPA2A ( am_I_Root, HcoState, Q, LonE, LatE, CLct, RC )
      DO II=1,IIPAR
      DO JJ=1,JJPAR
       FLUXSALC(II,JJ) = Q(II,JJ) 
      ENDDO
      ENDDO 

      CALL NcRd( Q, fA1, TRIM("Br2_TOTAL"), st3d, ct3d )
!      CALL REGRID_MAPA2A ( am_I_Root, HcoState, Q, LonE, LatE, CLct, RC )
      DO II=1,IIPAR
      DO JJ=1,JJPAR
       FLUXBr2(II,JJ) = Q(II,JJ) 
      ENDDO
      ENDDO 

      CALL NcRd( Q, fA1, TRIM("BrSALA_TOTAL"), st3d, ct3d )
!      CALL REGRID_MAPA2A ( am_I_Root, HcoState, Q, LonE, LatE, CLct, RC )
      DO II=1,IIPAR
      DO JJ=1,JJPAR
       FLUXBrSALA(II,JJ) = Q(II,JJ) 
      ENDDO
      ENDDO 

      CALL NcRd( Q, fA1, TRIM("BrSALC_TOTAL"), st3d, ct3d )
!      CALL REGRID_MAPA2A ( am_I_Root, HcoState, Q, LonE, LatE, CLct, RC )
      DO II=1,IIPAR
      DO JJ=1,JJPAR
       FLUXBrSALC(II,JJ) = Q(II,JJ) 
      ENDDO
      ENDDO 

    CALL NcCl( fA1 )
   ENDIF !OFFSEASALT


    !=================================================================
    ! PASS TO HEMCO STATE AND UPDATE DIAGNOSTICS 
    !=================================================================

    ! SALA 
    IF ( Inst%IDTSALA > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXSALA, Inst%IDTSALA, &
                         RC,        ExtNr=Inst%ExtNrSS )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'HCO_EmisAdd error: FLUXSALA', RC )
          RETURN 
       ENDIF
    ENDIF

    ! SALC 
    IF ( Inst%IDTSALC > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXSALC, Inst%IDTSALC, & 
                         RC,        ExtNr=Inst%ExtNrSS )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'HCO_EmisAdd error: FLUXSALC', RC )
          RETURN 
       ENDIF

    ENDIF

    ! BR2 
    IF ( Inst%CalcBr2 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXBr2, Inst%IDTBr2, & 
                         RC,        ExtNr=Inst%ExtNrSS )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'HCO_EmisAdd error: FLUXBr2', RC )
          RETURN 
       ENDIF

    ENDIF

    ! Bromine incorporated into sea salt 
    IF ( Inst%CalcBrSalt ) THEN

       ! Scale BrSalX emissions to SalX
       FluxBrSalA = Inst%BrContent * FluxSalA
       FluxBrSalC = Inst%BrContent * FluxSalC

       ! Add flux to emission array
       CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXBrSalA, Inst%IDTBrSalA, & 
                         RC,        ExtNr=Inst%ExtNrSS )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXBrSalA', RC )
          RETURN 
       ENDIF

       ! Add flux to emission array
       CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXBrSalC, Inst%IDTBrSalC, & 
                         RC,        ExtNr=Inst%ExtNrSS )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXBrSalC', RC )
          RETURN 
       ENDIF

    ENDIF

    ! MOPO 
    IF ( Inst%IDTMOPO > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXMOPO, Inst%IDTMOPO, & 
                         RC,        ExtNr=Inst%ExtNrMPOA )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'HCO_EmisAdd error: FLUXMOPO', RC )
          RETURN 
       ENDIF

    ENDIF

    ! MOPI
    IF ( Inst%IDTMOPI > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXMOPI, Inst%IDTMOPI, & 
                         RC,        ExtNr=Inst%ExtNrMPOA )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'HCO_EmisAdd error: FLUXMOPI', RC )
          RETURN 
       ENDIF

    ENDIF
      
    ! Cleanup
    Inst => NULL()

    ! Leave w/ success
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_SeaSalt_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_SeaSalt_Init
!
! !DESCRIPTION: Subroutine HcoX\_SeaSalt\_Init initializes all
!  extension variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_SeaSalt_Init( am_I_Root, HcoState, ExtName, ExtState, RC ) 
!
! !USES:
!
    USE HCO_State_Mod,          ONLY : HCO_GetHcoID
    USE HCO_STATE_MOD,          ONLY : HCO_GetExtHcoID
    USE HCO_ExtList_Mod,        ONLY : GetExtNr
    USE HCO_ExtList_Mod,        ONLY : GetExtOpt
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root   ! root CPU?
    TYPE(HCO_State),  POINTER        :: HcoState    ! HEMCO state object 
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName     ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState    ! Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC          ! Return status
!
! !REVISION HISTORY:
!  15 Dec 2013 - C. Keller   - Initial version
!  07 Oct 2014 - C. Keller   - Allow wind scale factor be set in config file
!  09 Jul 2015 - E. Lundgren - Add marine organic aerosols (B.Gantt, M.Johnson)
!  21 Sep 2016 - R. Yantosca - Bug fix: don't initialize SS_DEN before
!                              it is allocated.  This causes a segfault.
!  29 Dec 2017 - C. Keller   - Bug fix: define index location of BrSALA and
!                              BrSALC based upon # of species ID entries.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                        :: ExtNrSS
    INTEGER                        :: N, R, AS
    REAL*8                         :: A, B, R0, R1
    REAL*8                         :: CONST_N
    CHARACTER(LEN=255)             :: MSG
    INTEGER                        :: nSpcSS, nSpcMPOA, minLen
    REAL*8                         :: SALA_REDGE_um(2), SALC_REDGE_um(2)
    REAL(dp)                       :: tmpScale
    LOGICAL                        :: FOUND
    INTEGER, ALLOCATABLE           :: HcoIDsSS(:)
    INTEGER, ALLOCATABLE           :: HcoIDsMPOA(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNamesSS(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNamesMPOA(:)
    TYPE(MyInst), POINTER          :: Inst

    !=================================================================
    ! HCOX_SeaSalt_Init begins here!
    !=================================================================

    ! Extension number for seasalt
    ExtNrSS = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNrSS <= 0 ) RETURN

    ! Enter 
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_SeaSalt_Init (hcox_seasalt_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Create Instance
    Inst => NULL()
    CALL InstCreate ( ExtNrSS, ExtState%SeaSalt, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot create SeaSalt instance', RC )
       RETURN
    ENDIF
    ! Also fill ExtNrSS - this is the same as the parent ExtNr 
    Inst%ExtNrSS = ExtNrSS

    ! Check for marine organic aerosols option
    Inst%ExtNrMPOA = GetExtNr( HcoState%Config%ExtList, 'MarinePOA' )
 
    ! ---------------------------------------------------------------------- 
    ! Get species IDs and settings 
    ! ---------------------------------------------------------------------- 
  
    ! Read settings specified in configuration file
    ! Note: the specified strings have to match those in 
    !       the config. file!
    CALL GetExtOpt( HcoState%Config, Inst%ExtNrSS, 'Emit Br2', &
                     OptValBool=Inst%CalcBr2, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    IF ( Inst%CalcBr2 ) THEN
       minLen = 3
       CALL GetExtOpt( HcoState%Config, Inst%ExtNrSS, 'Br2 scaling', &
                       OptValDp=Inst%Br2Scale, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ELSE
       minLen   = 2
       Inst%IDTBr2   = -1
       Inst%Br2Scale = 1.0d0
    ENDIF

    Call GetExtOpt ( HcoState%Config, Inst%ExtNrSS, 'Model sea salt Br-', &
    	 	     OptValBool=Inst%CalcBrSalt, RC=RC )
    IF ( Inst%CalcBrSalt ) THEN
       minLen = minLen+2
       CALL GetExtOpt( HcoState%Config, Inst%ExtNrSS, 'Br- mass ratio', &
       	    OptValDp=Inst%BrContent, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ELSE
       Inst%IDTBrSALA = -1
       Inst%IDTBrSALC = -1
       Inst%BrContent = 0.0d0
    ENDIF

    ! Get HEMCO species IDs
    CALL HCO_GetExtHcoID( HcoState,   Inst%ExtNrSS, HcoIDsSS,     &
                          SpcNamesSS, nSpcSS,  RC           )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( nSpcSS < minLen ) THEN
       MSG = 'Not enough sea salt emission species set' 
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC ) 
       RETURN
    ENDIF
    Inst%IDTSALA = HcoIDsSS(1) 
    Inst%IDTSALC = HcoIDsSS(2)
    IF ( Inst%CalcBr2 ) Inst%IDTBR2 = HcoIDsSS(3)
    IF ( Inst%CalcBrSalt ) Inst%IDTBrSALA = HcoIDsSS(nSpcSS-1)
    IF ( Inst%CalcBrSalt ) Inst%IDTBrSALC = HcoIDsSS(nSpcSS)
    
    ! Get the marine organic aerosol species defined for MarinePOA option
    IF ( Inst%ExtNrMPOA > 0 ) THEN
       CALL HCO_GetExtHcoID( HcoState,     Inst%ExtNrMPOA, HcoIDsMPOA,  &
                             SpcNamesMPOA, nSpcMPOA,  RC          )
       IF ( RC /= HCO_SUCCESS ) RETURN
       Inst%IDTMOPO = HcoIDsMPOA(1)
       Inst%IDTMOPI = HcoIDsMPOA(2)
    ENDIF

    ! Get aerosol radius'
    SALA_REDGE_um(:) = 0.0d0
    SALC_REDGE_um(:) = 0.0d0
    CALL GetExtOpt( HcoState%Config, Inst%ExtNrSS, 'SALA lower radius', &
                    OptValDp=SALA_REDGE_um(1), RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL GetExtOpt( HcoState%Config, Inst%ExtNrSS, 'SALA upper radius', &
                    OptValDp=SALA_REDGE_um(2), RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL GetExtOpt( HcoState%Config, Inst%ExtNrSS, 'SALC lower radius', &
                    OptValDp=SALC_REDGE_um(1), RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL GetExtOpt( HcoState%Config, Inst%ExtNrSS, 'SALC upper radius', &
                    OptValDp=SALC_REDGE_um(2), RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Final Br2 flag
    Inst%CalcBr2 = ( Inst%CalcBr2 .AND. Inst%IDTBR2 > 0 )

    ! Final BrSalt flag
    Inst%CalcBrSalt = ( Inst%CalcBrSalt .and. Inst%IDTBrSALA > 0 .and. Inst%IDTBrSALC > 0 )

    ! The source function calculated with GEOS-4 2x2.5 wind speeds
    ! is too high compared to GEOS-5 at the same resolution. The 10m
    ! winds in GEOS-4 are too rapid. To correct this, apply a global
    ! scaling factor of 0.72 (jaegle 5/11/11)
    ! Now check first if this factor is specified in configuration file
    CALL GetExtOpt( HcoState%Config, Inst%ExtNrSS, 'Wind scale factor', & 
                    OptValDp=tmpScale, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN   
       tmpScale = 1.0d0
    ENDIF
    Inst%WindScale = tmpScale

    ! Verbose mode
    IF ( am_I_Root ) THEN
       MSG = 'Use sea salt aerosol emissions (extension module)'
       CALL HCO_MSG(HcoState%Config%Err,MSG, SEP1='-' )
 
       IF ( Inst%ExtNrMPOA > 0 ) THEN
          MSG = 'Use marine organic aerosols option'
          CALL HCO_MSG(HcoState%Config%Err,MSG, SEP1='-' )
       ENDIF 

       WRITE(MSG,*) 'Accumulation aerosol: ', TRIM(SpcNamesSS(1)),  &
                    ':', Inst%IDTSALA 
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) ' - size range       : ', SALA_REDGE_um
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) 'Coarse aerosol      : ', TRIM(SpcNamesSS(2)),  &
                     ':', Inst%IDTSALC
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) ' - size range       : ', SALA_REDGE_um
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) ' - wind scale factor: ', Inst%WindScale 
       CALL HCO_MSG(HcoState%Config%Err,MSG)
   
       IF ( Inst%CalcBr2 ) THEN
          WRITE(MSG,*) 'Br2: ', TRIM(SpcNamesSS(3)), Inst%IDTBr2
          CALL HCO_MSG(HcoState%Config%Err,MSG)
          WRITE(MSG,*) 'Br2 scale factor: ', Inst%Br2Scale
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

       IF ( Inst%CalcBrSalt ) THEN
          WRITE(MSG,*) 'BrSALA: ', TRIM(SpcNamesSS(nSpcSS-1)), Inst%IDTBrSALA
          CALL HCO_MSG(HcoState%Config%Err,MSG)
          WRITE(MSG,*) 'BrSALC: ', TRIM(SpcNamesSS(nSpcSS)), Inst%IDTBrSALC
          CALL HCO_MSG(HcoState%Config%Err,MSG)
          WRITE(MSG,*) 'Br- mass content: ', Inst%BrContent
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

       IF ( Inst%ExtNrMPOA > 0 ) THEN
          WRITE(MSG,*) 'Hydrophobic marine organic aerosol: ',        &
                       TRIM(SpcNamesMPOA(1)), ':', Inst%IDTMOPO 
          CALL HCO_MSG(HcoState%Config%Err,MSG)

          WRITE(MSG,*) 'Hydrophilic marine organic aerosol: ',        &
                       TRIM(SpcNamesMPOA(2)), ':', Inst%IDTMOPI 
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF
    ENDIF

    ! ---------------------------------------------------------------------- 
    ! Allocate module and subroutine arrays
    ! ---------------------------------------------------------------------- 

    ! Number of tracers dependent on MarinePOA (ewl, 7/9/15)
    IF ( Inst%ExtNrMPOA > 0 ) THEN
       Inst%NSALT = 4
    ELSE
       Inst%NSALT = 2
    ENDIF

    ALLOCATE ( Inst%NR  ( Inst%NSALT ), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate NR', RC )
       RETURN
    ENDIF
    Inst%NR = 0

    ALLOCATE ( Inst%SS_DEN  ( Inst%NSALT ), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate SS_DEN', RC )
       RETURN
    ENDIF
    Inst%SS_DEN = 2200.d0

    ALLOCATE ( Inst%SRRC   ( NR_MAX,   Inst%NSALT ), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate SRRC', RC )
       RETURN
    ENDIF
    Inst%SRRC = 0d0
    ALLOCATE ( Inst%SRRC_N ( NR_MAX,   Inst%NSALT ), STAT=AS ) 
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate SRRC_N', RC )
       RETURN
    ENDIF
    Inst%SRRC_N = 0d0
    ALLOCATE ( Inst%RREDGE ( 0:NR_MAX, Inst%NSALT ), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate RREDGE', RC )
       RETURN
    ENDIF
    Inst%RREDGE = 0d0
    ALLOCATE ( Inst%RRMID  ( NR_MAX,   Inst%NSALT ), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate RRMID', RC )
       RETURN
    ENDIF
    Inst%RRMID = 0d0

    ALLOCATE ( Inst%NDENS_SALA( HcoState%NX, HcoState%NY), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate NDENS_SALA', RC )
       RETURN
    ENDIF
    Inst%NDENS_SALA = 0.0_sp

    ALLOCATE ( Inst%NDENS_SALC( HcoState%NX, HcoState%NY), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate NDENS_SALC', RC )
       RETURN
    ENDIF
    Inst%NDENS_SALC = 0.0_sp

    IF ( Inst%ExtNrMPOA > 0 ) THEN 
   
       ! Allocate density of phobic marine organic aerosols
       ALLOCATE ( Inst%NDENS_MOPO( HcoState%NX, HcoState%NY), STAT=AS )
       IF ( AS/=0 ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate NDENS_MOPO', RC )
          RETURN
       ENDIF
       Inst%NDENS_MOPO = 0.0_sp
   
       ! Allocate density of philic marine organic aerosols
       ALLOCATE ( Inst%NDENS_MOPI( HcoState%NX, HcoState%NY), STAT=AS )
       IF ( AS/=0 ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate NDENS_MOPI', RC )
          RETURN
       ENDIF
       Inst%NDENS_MOPI = 0.0_sp

    ENDIF

    !=================================================================
    ! Define edges and midpoints of each incremental radius bin
    !=================================================================

    ! Constant [volume * time * other stuff??] 
    !CONST   = 4d0/3d0 * PI * DR * DTEMIS * 1.d-18 * 1.373d0

    !CONST_N = DTEMIS * DR * 1.373d0
    !  Constant for converting from [#/m2/s/um] to [#/m2]
    CONST_N = HcoState%TS_EMIS * (DR * BETHA)
 
    ! Do for accumulation, fine mode, and marine organics (if enabled)
    DO N = 1,Inst%NSALT

       ! Lower and upper limit of size bin N [um]
       ! Note that these are dry size bins. In order to
       ! get wet (RH=80%) sizes, we need to multiply by
       ! BETHA.

       ! Accumulation mode
       IF ( N==1 ) THEN
          R0 = SALA_REDGE_um(1) 
          R1 = SALA_REDGE_um(2)
          
       ! Coarse mode
       ELSEIF ( N==2 ) THEN 
          R0 = SALC_REDGE_um(1) 
          R1 = SALC_REDGE_um(2)
       
       ! Marine phobic (mj, bg, 7/9/15)
       ELSEIF ( N==3 ) THEN 
          R0 = SALA_REDGE_um(1) 
          R1 = SALA_REDGE_um(2)
          
       ! Marine philic (mj, bg, 7/9/15) 
       ELSEIF ( N==4 ) THEN 
          R0 = SALC_REDGE_um(1) 
          R1 = SALC_REDGE_um(2)
       ENDIF

       ! Number of radius size bins
       Inst%NR(N) = INT( ( ( R1 - R0 ) / DR ) + 0.5d0 ) 

       ! Error check
       IF ( Inst%NR(N) > NR_MAX ) THEN
          MSG = 'Too many bins'
          CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
          RETURN
       ENDIF

       ! Lower edge of 0th bin
       Inst%RREDGE(0,N) = R0
      
       ! Loop over the # of radius bins
       DO R = 1, Inst%NR(N)

          ! Midpoint of IRth bin
          Inst%RRMID(R,N)  = Inst%RREDGE(R-1,N) + ( DR / 2d0 )

          ! Upper edge of IRth bin
          Inst%RREDGE(R,N) = Inst%RREDGE(R-1,N) + DR 

          ! Sea salt base source [#/m2]. Note that the Gong formulation
          ! is for r80 (radius at 80% RH), so we need to multiply RRMID
          ! by the scaling factor BETHA=2.
          A           = 4.7*(1.+30.*(BETHA*Inst%RRMID(R,N)))             &
                       **(-0.017*(BETHA*Inst%RRMID(R,N))**(-1.44))
          B           = (0.433d0-LOG10(BETHA*Inst%RRMID(R,N))) / 0.433d0
          Inst%SRRC_N(R,N) = CONST_N * 1.373                            &
                      * (1.d0/(BETHA*Inst%RRMID(R,N))**(A))            &
                      * (1.d0+0.057d0*(BETHA*Inst%RRMID(R,N))**3.45d0) &
                      * 10d0**(1.607d0*EXP(-(B**2)))

          ! Sea salt base source [kg/m2]: multiply the number of particles
          ! by the dry volume multiplied by the dry density of sea-salt.
          Inst%SRRC(R,N) = Inst%SRRC_N(R,N) * 4d0/3d0 * HcoState%Phys%PI * 1.d-18 &
                         * Inst%SS_DEN( N ) * (Inst%RRMID(R,N))**3

          !-----------------------------------------------------------
          ! IMPORTANT NOTE!
          !
          ! In mathematics, "LOG" means "log10". 
          ! In Fortran,     "LOG" means "ln" (and LOG10 is "log10").
          !
          ! The following equations require log to the base 10, so 
          ! we need to use the Fortran function LOG10 instead of LOG. 
          ! (jaegle, bmy, 11/23/09)
          !-----------------------------------------------------------

!          ! Old Monahan et al. (1986) formulation
!          ! Sea salt base source [kg/m2]
!          CONST_N = DTEMIS * (DR * BETHA)
!          SRRC(R,N)  = CONST * SS_DEN( N )
!     &         * ( 1.d0 + 0.057d0*( BETHA * RRMID(R,N) )**1.05d0 )
!     &         * 10d0**( 1.19d0*
!     &           EXP(-((0.38d0-LOG10(BETHA*RRMID(R,N)))/0.65d0)**2))
!     &         / BETHA**2

!          ! Sea salt base source [#/m2] (bec, bmy, 4/13/05)
!          SRRC_N(R,N) = CONST_N * (1.d0/RRMID(R,N)**3)
!     &         * (1.d0+0.057d0*(BETHA*RRMID(R,N))**1.05d0)
!     &         * 10d0**(1.19d0*EXP(-((0.38d0-LOG10(BETHA*RRMID(R,N)))
!     &        /0.65d0)**2))/ BETHA**2

!### Debug
!###           WRITE( 6, 100 ) R,RREDGE(R-1,N),RRMID(R,N),RREDGE(R,N),SRRC(R,N)
!### 100        FORMAT( 'IR, R0, RRMID, R1: ', i3, 3f11.4,2x,es13.6 )
       ENDDO !R
    ENDDO !N

    !=======================================================================
    ! Create diagnostics. The number densities of both modes are always
    ! written into a diagnostics so that they can be used by other routines
    ! and from outside of HEMCO. These diagnostics just hold a pointer
    ! to the respective density arrays filled by the run method of this
    ! module.
    !=======================================================================
    CALL Diagn_Create ( am_I_Root,                          &
                        HcoState   = HcoState,              & 
                        cName      = 'SEASALT_DENS_FINE',   &
                        ExtNr      = Inst%ExtNrSS,          &
                        Cat        = -1,                    &
                        Hier       = -1,                    &
                        HcoID      = Inst%IDTSALA,          &
                        SpaceDim   = 2,                     &
                        OutUnit    = 'number_dens',         &
                        AutoFill   = 0,                     &
                        Trgt2D     = Inst%NDENS_SALA,       &
                        COL = HcoState%Diagn%HcoDiagnIDManual,      &
                        RC         = RC                      )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create ( am_I_Root,                          & 
                        HcoState   = HcoState,              & 
                        cName      = 'SEASALT_DENS_COARSE', &
                        ExtNr      = Inst%ExtNrSS,          &
                        Cat        = -1,                    &
                        Hier       = -1,                    &
                        HcoID      = Inst%IDTSALC,          &
                        SpaceDim   = 2,                     &
                        OutUnit    = 'number_dens',         &
                        AutoFill   = 0,                     &
                        Trgt2D     = Inst%NDENS_SALC,       &
                        COL = HcoState%Diagn%HcoDiagnIDManual, &
                        RC         = RC                      )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Create marine density diagnostics only if marine POA enabled
    IF ( Inst%ExtNrMPOA > 0 ) THEN

       CALL Diagn_Create ( am_I_Root,                          & 
                           HcoState   = HcoState,              & 
                           cName      = 'SEASALT_DENS_PHOBIC', &
                           ExtNr      = Inst%ExtNrMPOA,        &
                           Cat        = -1,                    &
                           Hier       = -1,                    &
                           HcoID      = Inst%IDTMOPO,          &
                           SpaceDim   = 2,                     &
                           OutUnit    = 'number_dens',         &
                           AutoFill   = 0,                     &
                           Trgt2D     = Inst%NDENS_MOPO,       &
                           COL = HcoState%Diagn%HcoDiagnIDManual, &
                           RC         = RC                      )
       IF ( RC /= HCO_SUCCESS ) RETURN
   
       CALL Diagn_Create ( am_I_Root,                          & 
                           HcoState   = HcoState,              & 
                           cName      = 'SEASALT_DENS_PHILIC', &
                           ExtNr      = Inst%ExtNrMPOA,        &
                           Cat        = -1,                    &
                           Hier       = -1,                    &
                           HcoID      = Inst%IDTMOPI,          &
                           SpaceDim   = 2,                     &
                           OutUnit    = 'number_dens',         &
                           AutoFill   = 0,                     &
                           Trgt2D     = Inst%NDENS_MOPI,       &
                           COL = HcoState%Diagn%HcoDiagnIDManual, &
                           RC         = RC                      )
       IF ( RC /= HCO_SUCCESS ) RETURN

    ENDIF

    !=======================================================================
    ! Activate this module and the fields of ExtState that it uses
    !=======================================================================

    ! Activate met fields used by this module
    ExtState%WLI%DoUse   = .TRUE.
    ExtState%ALBD%DoUse  = .TRUE.
    ExtState%TSKIN%DoUse = .TRUE.
    ExtState%U10M%DoUse  = .TRUE.
    ExtState%V10M%DoUse  = .TRUE.
    IF ( Inst%ExtNrMPOA > 0 ) THEN
       ExtState%CHLR%DoUse  = .TRUE.
    ENDIF

    ! Enable module
    !ExtState%SeaSalt = .TRUE.

    ! Return w/ success
    IF ( ALLOCATED(HcoIDsSS    ) ) DEALLOCATE(HcoIDsSS    )
    IF ( ALLOCATED(HcoIDsMPOA  ) ) DEALLOCATE(HcoIDsMPOA  )
    IF ( ALLOCATED(SpcNamesSS  ) ) DEALLOCATE(SpcNamesSS  )
    IF ( ALLOCATED(SpcNamesMPOA) ) DEALLOCATE(SpcNamesMPOA)

    CALL HCO_LEAVE( HcoState%Config%Err,RC ) 
 
  END SUBROUTINE HCOX_SeaSalt_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_SeaSalt_Final 
!
! !DESCRIPTION: Subroutine HcoX\_SeaSalt\_Final deallocates 
!  all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_SeaSalt_Final ( ExtState )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options      
!
! !REVISION HISTORY:
!  15 Dec 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !=================================================================
    ! HCOX_SeaSalt_Final begins here!
    !=================================================================

!    ! Cleanup module arrays
!    IF ( ALLOCATED ( NR         ) ) DEALLOCATE( NR         )    
!    IF ( ALLOCATED ( SS_DEN     ) ) DEALLOCATE( SS_DEN     )    
!    IF ( ALLOCATED ( SRRC       ) ) DEALLOCATE( SRRC       )
!    IF ( ALLOCATED ( SRRC_N     ) ) DEALLOCATE( SRRC_N     )
!    IF ( ALLOCATED ( RREDGE     ) ) DEALLOCATE( RREDGE     )
!    IF ( ALLOCATED ( RRMID      ) ) DEALLOCATE( RRMID      )
!
!    IF ( ASSOCIATED( NDENS_SALA ) ) DEALLOCATE( NDENS_SALA )
!    IF ( ASSOCIATED( NDENS_SALC ) ) DEALLOCATE( NDENS_SALC )    
!    IF ( ASSOCIATED( NDENS_MOPO ) ) DEALLOCATE( NDENS_MOPO )    
!    IF ( ASSOCIATED( NDENS_MOPI ) ) DEALLOCATE( NDENS_MOPI )
    CALL InstRemove ( ExtState%SeaSalt )

  END SUBROUTINE HCOX_SeaSalt_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Emit_SsaBr2
!
! !DESCRIPTION: Subroutine Emit\_SsaBr2 calculates aerosol emissions
!  of Br2.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Emit_SsaBr2( am_I_Root,  ExtState, HcoState, Inst, &
                          ilon, ilat, rmid, p_kgsalt, br2_emiss_kg, RC )
!
! !USE:
!
    USE HCO_Clock_Mod, ONLY : HcoClock_Get
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root ! root CPU?
    TYPE(Ext_State), POINTER        :: ExtState  ! Module options  
    TYPE(HCO_State), POINTER        :: HcoState  ! Output obj
    TYPE(MyInst),    POINTER        :: Inst      ! Instance
    INTEGER,         INTENT(IN)     :: ilon      ! Grid longitude index
    INTEGER,         INTENT(IN)     :: ilat      ! Grid latitude index
    REAL*8,          INTENT(IN)     :: rmid      ! Dry radius of aerosol
    REAL*8,          INTENT(IN)     :: p_kgsalt  ! SeaSalt aerosol production [kgNaCl]
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)  :: RC           ! Success or failure?
    REAL*8,          INTENT(OUT)    :: br2_emiss_kg ! Br2 emissions [kg NaCl]
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1)  Parrella, J. P., Jacob, D. J., Liang, Q., Zhang, Y., Mickley, L. J.,
!        Miller, B., Evans, M. J., Yang, X., Pyle, J. A., Theys, N., and Van
!        Roozendael, M.: Tropospheric bromine chemistry: implications for
!        present and pre-industrial ozone and mercury, Atmos. Chem. Phys., 12,
!        6723-6740, doi:10.5194/acp-12-6723-2012, 2012.
!  (2 ) Yang, X., Cox, R. A., Warwick, N. J., Pyle, J. A., Carver, G. D.,
!        O'Connor, F. M., and Savage, N. H.: Tropospheric bromine chemistry and
!        its impacts on ozone: A model study, J. Geophys. Res., 110, D23311,
!        doi:10.1029/2005JD006244, 2005.
!  (2 ) Yang, X., Pyle, J. A., and Cox, R. A.: Sea salt aerosol production and
!        bromine release: Role of snow on sea ice, Geophys. Res. Lett., 35,
!        L16815, doi:10.1029/2008GL034536, 2008.
!
! !REVISION HISTORY:
!  02 Mar 2010 - J. Parrella - Initial version
!  22 May 2012 - M. Payer    - Added ProTeX headers
!  08 Aug 2012 - M. Payer    - Modified for size-dependent depletion factors
!                              from Yang et al. (2008)
!  07 Aug 2013 - C. Keller   - Moved to SeaSalt_mod.F 
!  15 Dec 2013 - C. Keller   - Now a HEMCO extension 
!  19 Oct 2015 - C. Keller   - Now use lon and lat index to work on curvilinear
!                              grids.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    REAL*4, PARAMETER :: dfmax=0.7
    REAL*4, PARAMETER :: dfmin=0.1
    REAL*4, PARAMETER :: Ra=0.00223     ! Ratio of Br/NaCl [g/g]

    ! Use size-dependent sea salt bromine depletion factors from
    ! Table 1 of Yang et al. (2008)
    REAL*8, PARAMETER :: dmid_ref(10) = (/  0.2d0,   0.4d0,  & 
                                            0.8d0,   1.0d0,  &
                                            1.25d0,  1.5d0,  &
                                            2.0d0,   4.0d0,  &
                                            5.0d0,   10.0d0   /) 
    REAL*8, PARAMETER :: df_size(10)  = (/ -3.82d0, -2.54d0, &
                                            0.0d0,   0.23d0, &
                                            0.38d0,  0.37d0, &
                                            0.31d0,  0.21d0, &
                                            0.16d0,  0.11d0   /)
!
! !LOCAL VARIABLES:
!

    INTEGER :: month, IDF
    REAL*8  :: DF
    REAL*8  :: dmid         ! Dry diameter of aerosol [um]
    REAL*8  :: seasonal     ! Seasonal depletion factor

    !=================================================================
    ! EMIT_SSABr2 begins here!
    !=================================================================

    ! Dry diameter of aerosol [um]
    dmid = rmid * 2

    ! only do calculation if we're inside the
    ! range of aerosol sizes observed to be
    ! depeleted in bromide.
    IF ( (dmid < 0.2) .or. (dmid > 10.0) ) THEN
       br2_emiss_kg = 0.d0
       RC = HCO_SUCCESS
       RETURN
    ENDIF

    ! store the month
    CALL HcoClock_Get( am_I_Root, HcoState%Clock, cMM=month, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! --------------------------------------------
    ! 1. Calculate Depletion Factor DF, based on:
    !    (a) sea salt diameter (b) month and (c) latitude.
    !
    ! following Yang et al. 2005, 2008
    ! --------------------------------------------

    ! Sort into diameter bins
    IF (      dmid <= 0.4  ) THEN
       IDF = 1
    ELSE IF ( dmid <= 0.8  ) THEN
       IDF = 2
    ELSE IF ( dmid <= 1.0  ) THEN
       IDF = 3
    ELSE IF ( dmid <= 1.25 ) THEN
       IDF = 4
    ELSE IF ( dmid <= 1.5  ) THEN
       IDF = 5
    ELSE IF ( dmid <= 2.0  ) THEN
       IDF = 6
    ELSE IF ( dmid <= 4.0  ) THEN
       IDF = 7
    ELSE IF ( dmid <= 5.0  ) THEN
       IDF = 8
    ELSE 
       IDF = 9
    ENDIF

    ! Interpolate between sea salt diameters
    DF = df_size(IDF) + ( dmid            - dmid_ref(IDF) ) / &
                        ( dmid_ref(IDF+1) - dmid_ref(IDF) ) * &
                        ( df_size(IDF+1)  - df_size(IDF)  )


    ! Apply seasonality to latitudes south of 30S
    IF ( HcoState%Grid%YMID%Val(ilon,ilat) < -30.0 ) THEN
       ! Divide by mean value 0.4 = (dfmax+dfmin)/2 to keep
       ! seasonal dependence along with size dependence
       seasonal = ( dfmax + (dfmin - dfmax) / 2.d0 *                  &
                  ( sin( HcoState%Phys%PI*(month/6.d0 - 0.5) ) + 1 )) &
                  / 0.4
    ELSE
       ! no seasonal dependence for the NH
       seasonal = 1.d0
    ENDIF
    DF = DF * seasonal

    ! --------------------------------------------
    ! Now return the emissions for Br2 given the
    ! Sea-salt mass production.
    ! --------------------------------------------
    ! divide by 2 for stoichiometry of Br- to Br2
    br2_emiss_kg = p_kgsalt * Ra * DF / 2.0d0

    ! Apply Br scaling
    br2_emiss_kg = br2_emiss_kg * Inst%Br2Scale

    ! RETURN w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE Emit_SsaBr2
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstGet 
!
! !DESCRIPTION: Subroutine InstGet returns a poiner to the desired instance. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstGet ( Instance, Inst, RC, PrevInst ) 
!
! !INPUT PARAMETERS:
!
    INTEGER                             :: Instance
    TYPE(MyInst),     POINTER           :: Inst
    INTEGER                             :: RC
    TYPE(MyInst),     POINTER, OPTIONAL :: PrevInst
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst),     POINTER    :: PrvInst

    !=================================================================
    ! InstGet begins here!
    !=================================================================
 
    ! Get instance. Also archive previous instance.
    PrvInst => NULL() 
    Inst    => AllInst
    DO WHILE ( ASSOCIATED(Inst) ) 
       IF ( Inst%Instance == Instance ) EXIT
       PrvInst => Inst
       Inst    => Inst%NextInst
    END DO
    IF ( .NOT. ASSOCIATED( Inst ) ) THEN
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Pass output arguments
    IF ( PRESENT(PrevInst) ) PrevInst => PrvInst

    ! Cleanup & Return
    PrvInst => NULL()
    RC = HCO_SUCCESS

  END SUBROUTINE InstGet 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstCreate 
!
! !DESCRIPTION: Subroutine InstCreate creates a new instance. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstCreate ( ExtNr, Instance, Inst, RC ) 
!
! !INPUT PARAMETERS:
!
    INTEGER,       INTENT(IN)       :: ExtNr
!
! !OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(  OUT)    :: Instance
    TYPE(MyInst),  POINTER          :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(INOUT)    :: RC 
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst), POINTER          :: TmpInst
    INTEGER                        :: nnInst

    !=================================================================
    ! InstCreate begins here!
    !=================================================================

    ! ----------------------------------------------------------------
    ! Generic instance initialization 
    ! ----------------------------------------------------------------
    ! Initialize
    Inst => NULL()

    ! Get number of already existing instances
    TmpInst => AllInst
    nnInst = 0
    DO WHILE ( ASSOCIATED(TmpInst) )
       nnInst  =  nnInst + 1
       TmpInst => TmpInst%NextInst
    END DO

    ! Create new instance
    ALLOCATE(Inst)
    Inst%Instance = nnInst + 1
    Inst%ExtNr    = ExtNr 

    ! Init values
    Inst%ExtNrSS       = -1
    Inst%ExtNrMPOA     = -1
    Inst%IDTSALA       = -1
    Inst%IDTSALC       = -1
    Inst%IDTMOPI       = -1
    Inst%IDTMOPO       = -1
    Inst%IDTBr2        = -1
    Inst%IDTBrSALA     = -1
    Inst%IDTBrSALC     = -1
    Inst%CalcBr2       = .FALSE. 
    Inst%CalcBrSalt    = .FALSE. 
    Inst%Br2Scale      = 1.0 
    Inst%BrContent     = 1.0 
    Inst%WindScale     = 1.0 

    ! Attach to instance list
    Inst%NextInst => AllInst
    AllInst       => Inst

    ! Update output instance
    Instance = Inst%Instance

    ! ----------------------------------------------------------------
    ! Type specific initialization statements follow below
    ! ----------------------------------------------------------------

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE InstCreate
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstRemove 
!
! !DESCRIPTION: Subroutine InstRemove creates a new instance. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstRemove ( Instance ) 
!
! !INPUT PARAMETERS:
!
    INTEGER                         :: Instance 
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER                     :: RC
    TYPE(MyInst), POINTER       :: PrevInst
    TYPE(MyInst), POINTER       :: Inst

    !=================================================================
    ! InstRemove begins here!
    !=================================================================

    ! Init 
    PrevInst => NULL()
    Inst     => NULL()
    
    ! Get instance. Also archive previous instance.
    CALL InstGet ( Instance, Inst, RC, PrevInst=PrevInst )

    ! Instance-specific deallocation
    IF ( ASSOCIATED(Inst) ) THEN 
   
       ! Pop off instance from list
       IF ( ASSOCIATED(PrevInst) ) THEN

          ! Cleanup module arrays
          IF ( ASSOCIATED( Inst%NR         ) ) DEALLOCATE( Inst%NR         )    
          IF ( ASSOCIATED( Inst%SS_DEN     ) ) DEALLOCATE( Inst%SS_DEN     )    
          IF ( ASSOCIATED( Inst%SRRC       ) ) DEALLOCATE( Inst%SRRC       )
          IF ( ASSOCIATED( Inst%SRRC_N     ) ) DEALLOCATE( Inst%SRRC_N     )
          IF ( ASSOCIATED( Inst%RREDGE     ) ) DEALLOCATE( Inst%RREDGE     )
          IF ( ASSOCIATED( Inst%RRMID      ) ) DEALLOCATE( Inst%RRMID      )

          IF ( ASSOCIATED( Inst%NDENS_SALA ) ) DEALLOCATE( Inst%NDENS_SALA )
          IF ( ASSOCIATED( Inst%NDENS_SALC ) ) DEALLOCATE( Inst%NDENS_SALC )    
          IF ( ASSOCIATED( Inst%NDENS_MOPO ) ) DEALLOCATE( Inst%NDENS_MOPO )    
          IF ( ASSOCIATED( Inst%NDENS_MOPI ) ) DEALLOCATE( Inst%NDENS_MOPI )

          PrevInst%NextInst => Inst%NextInst
       ELSE
          AllInst => Inst%NextInst
       ENDIF
       DEALLOCATE(Inst)
       Inst => NULL() 
    ENDIF
   
   END SUBROUTINE InstRemove
!EOC
END MODULE HCOX_SeaSalt_Mod
