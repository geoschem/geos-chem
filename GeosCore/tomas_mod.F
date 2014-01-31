#if defined ( TOMAS )
! $Id: tomas_mod.f,v 1.2 2010/03/15 19:33:19 ccarouge Exp $
      MODULE TOMAS_MOD
!
!******************************************************************************
!  Module TOMAS_MOD contains variable specific to the TOMAS aerosol 
!  microphysics simulation, e.g. number of species, number of size-bins, mass
!  per particle bin boundaries and arrays used inside the microphysics 
!  subroutine. (win, 7/9/07)
!
!  NOTE:
!  This module also contains what used to be in sizecode.COM header file 
!  containing common blocks for TOMAS aerosol microphysics algorithm 
!  originally implemented in GISS GCM-II' by Peter Adams. Below is the original
!  comment from sizecode.COM
!
!     This header file includes all the variables used by the
!     size-resolved aerosol microphysics code incorporated into
!     the GISS GCM II' by Peter Adams.  The microphysics algorithm
!     conserves aerosol number and mass using the schemes developed
!     by Graham Feingold and others.
!
!  Tzivion, S., Feingold, G., and Levin, Z., An Efficient
!   Numerical Solution to the Stochastic Collection Equation,
!   J. Atmos. Sci., 44, 3139-3149, 1987.
!  Feingold, G., Tzivion, S., and Levin, Z., Evolution of
!   Raindrop Spectra. Part I: Solution to the Stochastic
!   Collection/Breakup Equation Using the Method of Moments,
!   J. Atmos. Sci., 45, 3387-3399, 1988.
!  Tzivion, S., Feingold, G., and Levin, Z., The Evolution of
!   Raindrop Spectra. Part II: Collisional Collection/Breakup
!   and Evaporation in a Rainshaft, J. Atmos. Sci., 46, 3312-
!   3327, 1989.
!  Feingold, G., Levin, Z., and Tzivion, S., The Evolution of
!   Raindrop Spectra. Part III: Downdraft Generation in an
!   Axisymmetrical Rainshaft Model, J. Atmos. Sci., 48, 315-
!   330, 1991.
!
!  The algorithms described in these papers have been extended
!  to include multicomponent aerosols and modified for a moving
!  sectional approach.  Using this approach, the boundaries
!  between size bins are defined in terms of dry aerosol mass
!  such that the actual sizes of the sections move as water
!  is added to or lost from the aerosol.
!
!  All of the subroutines needed for this aerosol microphysics
!  algorithm use only their own internal variables or the ones
!  listed here.  GISS GCM II' variables are not used (a driver
!  subroutine performs the necessary swapping between the GCM
!  and the microphysics code).  The microphysics code is,
!  therefore, completely modular.
!
!  Module Variables:
!  ============================================================================
!  (1 ) IBINS                  : Number of size bins 
!  (2 ) ICOMP                  : Number of aerosol mass species + 1 for number
!  (3 ) Nk                     : Aerosol number internal array
!  (4 ) Mk                     : Aerosol mass internal array
!  (5 ) Gc                     : Condensing gas 
!  (6 ) Nkd                    : Aerosol number diagnostic array
!  (7 ) Mkd                    : Aerosol mass diagnostic array
!  (8 ) Gcd                    : Condensing gas diagnostic array
!  (9 ) Xk                     : Size bin boundary in dry mass per particle
!  (10) MOLWT                  : Aerosol molecular weight
!  (11) SRTSO4                 : ID of sulfate
!  (12) SRTNACL                : ID of sea-salt
!  (13) SRTH2O                 : ID of aerosol water
!  (14) SRTECOB                : ID of hydrophobic EC
!  (15) SRTECIL                : ID of hydrophilic EC
!  (16) SRTOCOB                : ID of hydrophobic OC
!  (17) SRTOCIL                : ID of hydrophilic OC
!  (18) dust??                 : ID of internally mixed dust (future work)
!  (19) dust??                 : ID of externally mixed dust (future work)
!  (20) BOXVOL                 : Grid box volume  (cm3)
!  (21) BOXMASS                : Grid box air mass (kg)
!  (22) TEMPTMS                : Temperature (K) of each grid box 
!  (23) PRES                   : Pressure (Pa) in grid box 
!  (24) RHTOMAS                : Relative humidity (0-1) 
!  (25) BINACT1                : Activated bin as a function of composition
!  (26) FRACTION1              : Activated fraction as a fcn of composition
!  (27) IDIAG                  : Number of diagnostic tracer (NH4 and H2O)
!
!  Module Routines:
!  ============================================================================
!  (1 ) DO_TOMAS               : Driver subroutine to call microphysics and 
!                                dry deposition
!  (2 ) AEROPHYS               : Main microphysics calling individual processes
!  (3 ) AQOXID                 : Aqueous oxidation 
!  (4 ) MULTICOAG              : Coagulation
!  (5 ) NUCLEATION             : Nucleation
!  (6 ) SO4COND                : H2SO4 condensation
!  (7 ) TMCOND                 : Condensation calculation
!  (8 ) AERODIAG               : Save changes to diagnostic arrays
!  (9 ) INIT_TOMAS             : Initialize TOMAS microphysics variables  
!  (10) EZWATEREQM             : Update aerosol water with current RH 
!  (11) AERO_DIADEN            : Calculate the current wet diameter and density
!  (12) CHECKMN                : Wrapper for call MNFIX for error checking
!  (13) MNFIX                  : Check for error (out-of-bounds, negative, etc.)
!  (14) CLEANUP_TOMAS          : Clean up allocated array to free memory.
!  (15) STORENM                : Store value for diagnostic purposes
!  (16) SOACOND                : Simplified condensation for SOA 
!  (17) READBINACT             : Read data file to create activating bin LUT
!  (18) READFRACTION           : Read data file to create scavenging fraction LUT
!  (19) GETFRACTION            : Use lookup table and determine scavenging fraction
!  (20) GETACTBIN              : Use lookup table and determine activating bin 
!  (21) EZWATEREQM2            : Like EZWATEREQM but for calling from outside TOMAS
!  (22) SUBGRIDCOAG            : Calculate emission w/ number reduced by subgrid coag
!
!  Functions:
!  ============================================================================
!  (1 ) AERODENS               : Return the current aerosol density
!  (2 ) DMDT_INT               : Return droplet growth analytical solution 
!  (3 ) GASDIFF                : Return gas diffusion constant
!  (4 ) GETDP                  : Calculate the current aerosol diameter
!  (5 ) SPINUP                 : Return TRUE or FALSE if the current time is 
!                                after the spin-up period
!
!  NOTES:
!  (1 ) Add BINACT1 and FRACTION1 for compositionally-resolved scavenging after
!        OC and EC are included in the simulation (win, 9/10/07)
!  23 Apr 2013 - R. Yantosca - Now use findFreeLUN to get file LUN's
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars
      INTEGER                        :: ICOMP,   IDIAG         

      ! Note that there is a parameter declared in CMN_SIZE called TOMASBIN --> which defines how many size bins (win, 7/7/09)
#if   defined ( TOMAS12 )
      INTEGER,           PARAMETER   :: IBINS = 12
#elif defined ( TOMAS15 )
      INTEGER,           PARAMETER   :: IBINS = 15
#elif defined ( TOMAS40 )
      INTEGER,           PARAMETER   :: IBINS = 40
#else
      INTEGER,           PARAMETER   :: IBINS = 30
#endif

      INTEGER,           PARAMETER   :: SRTSO4  = 1
      INTEGER,           PARAMETER   :: SRTNACL = 2
      INTEGER,           PARAMETER   :: SRTECIL = 3
      INTEGER,           PARAMETER   :: SRTECOB = 4
      INTEGER,           PARAMETER   :: SRTOCIL = 5
      INTEGER,           PARAMETER   :: SRTOCOB = 6
      INTEGER,           PARAMETER   :: SRTDUST = 7
      INTEGER,           PARAMETER   :: SRTNH4  = 8
      INTEGER,           PARAMETER   :: SRTH2O  = 9

      INTEGER,           PARAMETER   :: ICOMPHARD = 9
     


      ! Arrays
      REAL*8,            SAVE        :: Xk(IBINS+1)
      REAL*4,   save,    ALLOCATABLE :: MOLWT(:)
      INTEGER,           SAVE        :: BINACT1(101,101,101)
      REAL*8,            SAVE        :: FRACTION1(101,101,101)
      INTEGER,           SAVE        :: BINACT2(101,101,101)
      REAL*8,            SAVE        :: FRACTION2(101,101,101)
      REAL*8                            AVGMASS(IBINS) ! Average mass per particle
                                                       ! mid-range of size bin
                                                       ! [kg/no.]
      REAL*8                         ::  cosmic_ions(72,46,9) !careful, this is not GIGC safe!
                                                       ! [ion pairs kg^-1 s^-1]

!      DATA AVGMASS/ 1.4142d-21, 
!     &              2.8284d-21, 5.6569d-21, 1.1314d-20, 2.2627d-20,
!     &  4.5255d-20, 9.0510d-20, 1.8102d-19, 3.6204d-19, 7.2408d-19,
!     &  1.4482d-18, 2.8963d-18, 5.7926d-18, 1.1585d-17, 2.3170d-17,
!     &  4.6341d-17, 9.2682d-17, 1.8536d-16, 3.7073d-16, 7.4146d-16,
!     &  1.4829d-15, 2.9658d-15, 5.9316d-15, 1.1863d-14, 2.3727d-14,
!     &  4.7453d-14, 9.4906d-14, 1.8981d-13, 3.7963d-13, 7.5925d-13 /

      INTEGER                        :: bin_nuc = 1, tern_nuc = 0  !Switches for nucleation type.
      INTEGER                        :: act_nuc = 1 ! in BL
      INTEGER                        :: ion_nuc = 0 ! 1 for modgil, 2 for Yu
      INTEGER                        :: absall  = 1 ! 1 for soa absorb to all spec
!napari nucleation tuned by factor of 1.0D-5

      REAL*8                         :: soaareafrac=1.0d0 ! fraction of SOA that goes
                                                          ! to SA rather than into mass
      

      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DO_TOMAS( am_I_Root, Input_Opt, 
     &                     State_Met, State_Chm, RC )
!
!******************************************************************************
!  Subroutine DO_TOMAS is the driver subroutine that calls the appropriate
!  aerosol microphysics and dry deposition subroutines. (win, 7/23/07)
!
!  NOTES: 
!  26 Apr 2013 - R. Yantosca - Remove LTOMAS flag, that is now obsolete.
!  31 May 2013 - R. Yantosca - Now pass State_Chm object via the arg list
!  31 May 2013 - R. Yantosca - Remove LDRYD; replace w/ Input_Opt%LDRYD
!  31 May 2013 - R. Yantosca - Pass State_Chm to AEROPHYS, AERO_DRYDEP
!******************************************************************************
!
      ! References to F90 modules
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState

      ! Arguments
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?

      !=================================================================
      ! DO_TOMAS begins here
      !=================================================================

      ! Assume success
      RC    = GIGC_SUCCESS

      ! Do TOMAS aerosol microphysics
      CALL AEROPHYS( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

#if defined ( DEBUG ) 
      print *, 'call checkmn in tomas_mod:222'
#endif
      CALL CHECKMN( 0,         0,         0,     Input_Opt, 
     &              State_Met, State_Chm, 'kg', 'Before Aerodrydep')

      ! Do dry deposition
      IF ( Input_Opt%LDRYD ) THEN
         CALL AERO_DRYDEP
     &        ( am_I_Root, Input_Opt, State_Met, State_Chm,RC )
      ENDIF

#if defined ( DEBUG ) 
      print *, 'call checkmn in tomas_mod:229'
#endif
      CALL CHECKMN( 0,         0,         0,    
     &              Input_Opt, State_Met, State_Chm, 
     &              'kg',      'Before exiting DO_TOMAS' )

      ! Return to calling program 
      END SUBROUTINE DO_TOMAS

!------------------------------------------------------------------------------

      SUBROUTINE AEROPHYS( am_I_Root, Input_Opt, 
     &                     State_Met, State_Chm, RC )
!
!******************************************************************************
!  Subroutine AEROPHYS does aerosol microphysics, including nucleation, 
!  coagulation, and condensation.
!
!  NOTES:
!  (1 ) Major change in nucleation scheme from binary w/ critical concentration 
!       method to calculating nucleation rate with the option to do binary 
!       (Vehkamaki et al. 2002) or ternary nucleation (Napari et al, 2002).
!       Also add total N and total S mass conservation check.
!       (win, 9/30/08)
!  (2 ) Use H2SO4RATE array which is accessible only if ND65 family is setup
!       for SO4 - refer to diag_pl_mod.f (win, 9/30/08)
!  (3 ) Reference to NFAMILIES from comode.h  --> should change in the future.
!       Prefer not to use .h file (win, 9/30/08)
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  28 Jan 2014 - R. Yantosca - Avoid array temporaries in call to CHECK_VALUE
!  29 Jan 2014 - R. Yantosca - Avoid array temporaries in call to NH4BULKTOBIN
!  30 Jan 2014 - R. Yantosca - Add missing variables to !$OMP+PRIVATE statement
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,          ONLY : ERROR_STOP, IT_IS_NAN, check_value
      USE DIAG_MOD,           ONLY : AD61,  AD61_INST
      USE DIAG_PL_MOD,        ONLY : H2SO4RATE, FAM_NAME
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE TRACERID_MOD,       ONLY : IDTNK1, IDTH2SO4
      USE TRACERID_MOD,       ONLY : IDTAW1, IDTSF1,  IDTSO4
      USE TRACERID_MOD,       ONLY : IDTNH3, IDTNH4
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_STRAT
     
      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND60, ND61,  LD61
      USE COMODE_LOOP_MOD          ! NFAMILIES

      ! Arguments
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?

      ! Local variables
      INTEGER           :: I, J, L, N, JC, K !counters
      INTEGER           :: MPNUM          !microphysical process id #
      REAL*4            :: ADT !aerosol microphysics time step (seconds)
      REAL*8            :: QSAT         !used in RH calculation
      INTEGER           :: TRACNUM
      REAL*8            :: FRAC

      ! Arguments for CHECK_VALUE; avoids array temporaries (bmy, 1/28/14)
      CHARACTER(LEN=255) :: ERR_VAR
      CHARACTER(LEN=255) :: ERR_MSG
      INTEGER            :: ERR_IND(4)

!---------
!sfarina - move definitions of these from module common
!to within each sub and pass them around for openmp
      REAL*8                         :: Nk(IBINS), Nkd(IBINS)
      REAL*8                         :: Mk(IBINS, ICOMP)
      REAL*8                         :: Mkd(IBINS,ICOMP)
      REAL*8                         :: Gc(ICOMP - 1)
      REAL*8                         :: Gcd(ICOMP - 1)

      REAL*4                         :: BOXVOL,  BOXMASS, TEMPTMS
      REAL*4                         :: PRES,    RHTOMAS

      REAL*8                         ::  surf_area     ! aerosol surface area [micon^2 cm^-3]
      REAL*8                         ::  ionrate       ! ion pair formation 
                                                       ! rate [ion pairs cm^-3 s^-1]
!---------

      REAL*8            :: Nkout(ibins), Mkout(ibins,icomp)
      REAL*8            :: Gcout(icomp-1)
      REAL*8            :: Nknuc(ibins), Mknuc(ibins,icomp)
      REAL*8            :: Nkcond(ibins), Mkcond(ibins,icomp)
      REAL*8            :: fn  ! nucleation rate of clusters cm-3 s-1
      REAL*8            :: fn1 ! formation rate of particles to first size bin cm-3 s-1
      REAL*8            :: nucrate(JJPAR,LLPAR), nucrate1(JJPAR,LLPAR)
      REAL*8            :: tot_n_1, tot_n_1a, tot_n_2, tot_n_i ! used for nitrogen mass checks
      REAL*8            :: tot_s_1, tot_s_1a, tot_s_2 ! used for sulfur mass checks
      REAL*8            :: h2so4rate_o ! H2SO4rate for the specific grid cell
      REAL*8            :: TOT_Mk, TOT_nk  ! for checking mass and number

      REAL*8            :: transfer(ibins)
      LOGICAL           :: PRINTNEG  !<step4.0-temp> (win, 3/24/05)
      logical           :: ERRORSWITCH  !<step4.2> To see where mnfix found negative value (win, 9/12/05)
      logical           :: errspot   !<step4.4> To see where so4cond found errors (win, 9/21/05)
      logical           :: printdebug !<step4.3> Print out for debugging (win, 9/16/05)
      logical           :: COND, COAG, NUCL !<step5.1> switch for each process (win 4/8/06)
      integer :: iob, job,lob !Just declare in case I want to debug (4/8/06)
      data iob, job, lob /  20  ,       33    ,       2 /
      real*8           :: NH3_to_NH4, CEPS
      parameter ( CEPS=1.d-17 )

      real*8 igR
      parameter (igR=8.314) !Ideal gas constant J/mol.K

      !The following are constants used in calculating rel. humidity
      real*8 axcons, bxcons, bytf, tf  !for RH calculation
      parameter(axcons=1.8094520287589733,
     &          bxcons=0.0021672473136556273,
     &           bytf=0.0036608580560606877, tf=273.16)
      !lhe and lhs are the latent heats of evaporation and sublimation
      
      logical, save     :: firsttime = .true.
      integer           :: num_iter

      real*8      cplevels(9) ! cosmic ray pressure levels (for interp)
      data        cplevels  / 959., 894., 786., 634., 468.,
     &                        321., 201., 103., 27. /

      integer     lev
      real*8      weight

      double precision soil_ions(9) ! ion pairs cm-3 s-1 from radioactive elements in soil
      data             soil_ions / 5.,3.3,1.8,0.7,0.,0.,0.,0.,0./

      ! For fields from Input_Opt
      LOGICAL :: LPRT

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! AEROPHYS begins here
      !=================================================================
      
      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      STT                      => State_Chm%Tracers

      ! Get logical values from Input_Opt
      LPRT = ( Input_Opt%LPRT .and. am_I_Root )

      ! Initialize debugging and error-signal switches
      printneg = .false.
      errorswitch = .false.
      PRINTDEBUG = .false.
      ERRSPOT = .FALSE.

      ! Initialize switches for each microphysical process
      COND = .TRUE.
      COAG = .TRUE.
      NUCL = .TRUE.

      ! Initialize nucleation rate arrays
      nucrate(:,:)  = 0.d0
      nucrate1(:,:) = 0.d0

      ! Make sure there is access to the H2SO4 production rate array
      ! H2SO4RATE, which saves the H2SO4 production rate for EACH chemistry
      ! timestep from the ND65 fake-family tracer.  The family has to be set
      ! with at leaste one family with the family name PSO4 with SO4 as its 
      ! one member. (win, 9/30/08)     
      if (firsttime) then
         DO N = 1, NFAMILIES
            ! If family name 'PSO4' is found, then skip error-stop
            IF ( FAM_NAME(N) == 'PSO4') GOTO 1            
         ENDDO
         ! Family name 'PSO4' not found... exit with error message
         write(*,*)'-----------------------------------------------'
         write(*,*)' Need to setup ND65 family PSO4 with SO4 as '
         write(*,*)' a member to have H2SO4RATE array '
         write(*,*)'  ... need H2SO4RATE for nucl & cond in TOMAS'
         write(*,*)'-----------------------------------------------'
         CALL ERROR_STOP('AEROPHYS','Enter microphys')
 1       CONTINUE            

         write(*,*) 'AEROPHYS: This run uses coupled condensation-',
     &        'nucleation scheme with pseudo-steady state H2SO4'
         if(tern_nuc == 1) then
            write(*,*)'  Nucleation: Ternary (Napari ',
     &           'et al. 2002) and Binary (Vehkamaki et al. 2002)'
         else
            write(*,*)'  Nucleation: Binary (Vehkamaki et al. 2002)'
         endif
         
         firsttime = .false.
      endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% NOTE: THIS PARALLEL LOOP MAY BE ABLE TO BE REVERSED TO L-J-I
!%%% WHICH IS MUCH MORE EFFICIENT (bmy, 1/28/14)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
!$OMP+PRIVATE( ADT, PRES, TEMPTMS, BOXMASS, RHTOMAS, BOXVOL )
!$OMP+PRIVATE( printneg, ionrate, lev, weight, GC, N, NK, JC )
!$OMP+PRIVATE( MK, H2SO4rate_o, tot_n_1, k, tot_s_1, MPNUM )
!$OMP+PRIVATE( Nkd, Mkd, TOT_NK, TOT_MK, TRANSFER )
!$OMP+PRIVATE( Nkout,Mkout,Gcout,fn,fn1 )
!$OMP+PRIVATE( num_iter,Nknuc,Mknuc,Nkcond )
!$OMP+PRIVATE( Mkcond, ERRORSWITCH, tot_s_1a, tot_n_1a )
!$OMP+PRIVATE( ERR_VAR, ERR_MSG, ERR_IND )
!$OMP+PRIVATE( TRACNUM, NH3_TO_NH4, SURF_AREA )
!$OMP+SCHEDULE( DYNAMIC )
      DO I = 1, IIPAR
      DO J = 1, JJPAR
      DO L = 1, LLPAR

         ! Reset the AD61_INST array used for tracking instantaneous 
         ! certain rates.  As of now, tracking nucleation (win, 10/7/08)
         AD61_INST(I,J,L,:) = 0e0

         ! Skip stratospheric boxes
         IF ( ITS_IN_THE_STRAT( I, J, L, State_Met ) ) CYCLE

cvbn         write(890,89)I,J,L,STT(i,j,l,IDTH2SO4)


         !if(printdebug)          print *,'+++++',I,J,L,'inside Aerophys'
 
         ! Get info on this grid box
         ADT     = 3600.
         PRES    = GET_PCENTER(i,j,l)*100.0 ! in Pa       
         TEMPTMS = State_Met%T(I,J,L)
         BOXMASS = State_Met%AD(I,J,L)
         RHTOMAS = State_Met%RH(I,J,L)/ 1.e2 
         IF ( RHTOMAS > 0.99 ) RHTOMAS = 0.99
         BOXVOL  = State_Met%AIRVOL(I,J,L) * 1.e6 !convert from m3 -> cm3 
         
         printneg = .FALSE.


         ! determine ion rate
c         ionrate = 10. ! set as constant now !!jrp

#if  defined( GRID4x5)
         if( (pres / 100.) .lt. cplevels(9) ) then
            ionrate = cosmic_ions(i,j,9) * boxmass / boxvol
            if(State_Met%FRCLND(I,J) .gt. 0.2) then
               ionrate = ionrate + soil_ions(9)
            endif
         elseif((pres/100.) .gt. cplevels(1)) then
            ionrate = cosmic_ions(i,j,1) * boxmass / boxvol
            if( State_Met%FRCLND(I,J) .gt. 0.2 ) then
               ionrate = ionrate + soil_ions(1)
            endif
         else
            lev=2
            do while (pres / 100. .lt. cplevels(lev))
               lev=lev+1
            enddo
            weight=( cplevels( lev - 1 ) - pres / 100. ) /
     &             ( cplevels( lev - 1 ) - cplevels(lev) )
            ionrate=( cosmic_ions(i, j, lev    ) * weight +
     &                cosmic_ions(i, j, lev - 1) * (1.d0 - weight) )
     &              * boxmass / boxvol
            if( State_Met%FRCLND(I,J) .gt. 0.2) then
               ionrate=ionrate + ( soil_ions( lev   ) * weight +
     &                             soil_ions( lev-1 ) * (1.d0-weight) )
            endif
         endif
#else
         ionrate = 0.D0
#endif
         if(ionrate .le. 1.501) ionrate = 1.501

c         print*,'i',i,'j',j,'l',l,'ionrate',ionrate

         ! Initialize all condensible gas values to zero  
         ! Gc(srtso4) will remain zero until within cond_nuc where the
         ! pseudo steady state H2SO4 concentration will be put in this place.
         DO JC=1, ICOMP-1
            Gc(JC) = 0.d0
         ENDDO

         ! Swap STT into Nk, Mk, Gc arrays
         DO N = 1, IBINS
            NK(N) = STT(I,J,L,IDTNK1-1+N)
            DO JC = 1, ICOMP-IDIAG
               MK(N,JC) = STT(I,J,L,IDTNK1-1+N+JC*IBINS)

               IF( IT_IS_NAN( MK(N,JC) ) ) THEN
                  PRINT *,'+++++++ Found NaN in AEROPHYS ++++++++'
                  PRINT *,'Location (I,J,L):',I,J,L,'Bin',N,'comp',JC
               ENDIF

            ENDDO
            MK(N,SRTH2O) = STT(I,J,L,IDTAW1-1+N)
         ENDDO

         ! Get NH4 mass from the bulk mass and scale to bin with sulfate
         IF ( SRTNH4 > 0 ) THEN

!------------------------------------------------------------------------------
! Prior to 1/29/14:
! Rewrite subroutine call to avoid array temporaries
!            CALL NH4BULKTOBIN( MK(1:IBINS,SRTSO4), STT(I,J,L,IDTNH4),
!     &                         TRANSFER(1:IBINS)  )
!------------------------------------------------------------------------------
            CALL NH4BULKTOBIN( MK(:,SRTSO4), STT(I,J,L,IDTNH4), 
     &                         TRANSFER )
            MK(1:IBINS,SRTNH4) = TRANSFER(1:IBINS)
            Gc(SRTNH4) = STT(I,J,L,IDTNH3) 
         ENDIF

         
         ! Give it the pseudo-steady state value instead later (win,9/30/08)
         !GC(SRTSO4) = STT(I,J,L,IDTH2SO4)

         H2SO4rate_o = H2SO4RATE(I,J,L)  ! [kg s-1]


                                ! nitrogen and sulfur mass checks
                                ! get the total mass of N
         tot_n_1 = Gc(srtnh4)*14.d0/17.d0
         do k=1,ibins
            tot_n_1 = tot_n_1 + Mk(k,srtnh4)*14.d0/18.d0
         enddo
         
                                ! get the total mass of S
         tot_s_1 = H2SO4rate_o*adt*32.d0/98.d0
         do k=1,ibins
            tot_s_1 = tot_s_1 + Mk(k,srtso4)*32.d0/96.d0
         enddo

         
         if (printdebug.and.i==iob .and. j==job .and. l==lob ) then
            CALL DEBUGPRINT( Nk, Mk, I, J, L,
     &        'Begin aerophys' )
            print *,'H2SO4RATE ',H2SO4rate_o
         endif
        
         !*********************
         ! Aerosol dynamics
         !*********************

         !Do water eqm at appropriate times
         CALL EZWATEREQM( MK, RHTOMAS )

         !Fix any inconsistencies in M/N distribution (because of advection)
         CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd) 

         if(printdebug .and. 
     &        i==iob.and.j==job.and.l==lob) ERRORSWITCH =.TRUE. 

#if defined ( DEBUG ) 
!      print *, 'mnfix in tomas_mod:533'
#endif
         CALL MNFIX( NK, MK, ERRORSWITCH )
         IF ( ERRORSWITCH ) THEN
            PRINT *,'Aerophys: MNFIX found error at',I,J,L
            CALL ERROR_STOP('AEROPHYS-MNFIX (1)','Enter microphys')
         ENDIF

         MPNUM = 5    
         IF ( ND60 > 0 ) THEN
            CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL )
         ENDIF

         IF ( printdebug.and.i==iob .and. j==job .and. l==lob ) THEN
            CALL DEBUGPRINT( Nk, Mk, I, J, L,
     &           'After mnfix before cond/nucl' )
         ENDIF

         ! Before doing any cond/nucl/coag, check if there's any aerosol in the current box
         TOT_NK = 0.D0
         TOT_MK = 0.d0
         do k = 1, ibins
            TOT_NK = TOT_NK + Nk(K)
            do jc=1, icomp-idiag
               TOT_MK = TOT_MK + Mk(k,jc)
            enddo
         enddo
         
         if(TOT_NK .lt. 1.d-10) then
            if( .NOT. SPINUP(5.0)) then
               print *,'No aerosol in box ',I,J,L,'-->SKIP'
            endif
            CYCLE
         endif

         !-------------------------------------
         ! Condensation and nucleation (coupled)
         !-------------------------------------
         IF ( COND .AND. NUCL ) THEN

            if(printdebug .and. 
     &           i==iob.and.j==job.and.l==lob) ERRORSWITCH =.TRUE. 

         CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd) 

!            print*,'Before COND_NUC Gc(srtso4)=',Gc(srtso4)

            CALL COND_NUC(Nk,Mk,Gc,Nkout,Mkout,Gcout,fn,fn1,
     &             H2SO4rate_o,adt,num_iter,Nknuc,Mknuc,Nkcond,Mkcond, 
     &             ionrate, surf_area, BOXVOL, BOXMASS, TEMPTMS, PRES,
     &             RHTOMAS, ERRORSWITCH, l)

!            print*,'After COND_NUC Gcout(srtso4)=',Gcout(srtso4)
            
            IF ( ERRORSWITCH ) THEN
               PRINT *,'Aerophys: found error at',I,J,L
               CALL ERROR_STOP('AEROPHYS','After cond_nuc')
            ENDIF

           ERR_VAR = 'Gcout'
           ERR_MSG = 'After COND_NUC'
            ! check for NaN and Inf (win, 10/4/08)
            do jc = 1, icomp-1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% Prior to 1/28/14:
!%%% Eliminate array temporary in call to CHECK_VALUE.  Now use ERR_IND
!%%% array instead of a constructor list (/I,J,L,0/).  (bmy, 1/28/14)
!%%%              call check_value(Gcout(jc),(/I,J,L,0/),ERR_VAR,
!%%%    &              ERR_MSG)
!%%%!               call check_value(Gcout(jc),(/I,J,L,0/),'Gcout',
!%%%!     &              'After COND_NUC')
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               ERR_IND(1) = I
               ERR_IND(2) = J
               ERR_IND(3) = L
               ERR_IND(4) = 0
               call check_value( Gcout(jc), ERR_IND, ERR_VAR, ERR_MSG )

!               if( IT_IS_FINITE(Gcout(jc))) then
!                  print *,'xxxxxxxxx Found Inf in Gcout xxxxxxxxxxxxxx'
!                  print *,'Location ',I,J,L, 'comp',jc
!                  call debugprint( Nkout, Mkout, i,j,l,'After COND_NUC')
!                  stop
!               endif
            enddo

                               !get nucleation diagnostic
            DO N = 1, IBINS
               NK(N) = NKnuc(N)
               DO JC = 1, ICOMP
                  MK(N,JC) = MKnuc(N,JC)
               ENDDO
            ENDDO

            MPNUM = 3
            IF ( ND60 > 0 ) THEN
               CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL )
            ENDIF

            MPNUM = 7
            IF ( ND61 > 0 )  THEN
               CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL )
            ENDIF

            IF ( printdebug.and.i==iob .and. j==job .and. l==lob )  THEN
               CALL DEBUGPRINT( Nk, Mk, I, J, L,'After nucleation' )
            ENDIF
                               !get condensation diagnostic
            DO N = 1, IBINS
               NK(N) = NKcond(N)
               DO JC = 1, ICOMP
                  MK(N,JC) = MKcond(N,JC)
               ENDDO
            ENDDO

            Gc(srtnh4)=Gcout(srtnh4)
            Gc(srtso4)=Gcout(srtso4)

            MPNUM = 1
            IF ( ND60 > 0 ) THEN
               CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL )
            ENDIF

            IF ( printdebug.and.i==iob .and. j==job .and. l==lob ) THEN
               CALL DEBUGPRINT( Nk, Mk, I, J, L,'After condensation' )
            ENDIF

            nucrate(j,l)=nucrate(j,l)+fn
            nucrate1(j,l)=nucrate1(j,l)+fn1
            
            ! Write nucleation rate to diagnostric ND61 (win, 10/6/08)
            IF ( ND61 > 0 ) THEN
               IF ( L <= LD61 ) 
     &              AD61(I,J,L,2) = AD61(I,J,L,2) + fn

               ! Tracks nucleation rates instantaneously for planeflight
               AD61_INST(I,J,L,2) = fn
               
            ENDIF
               

            DO N = 1, IBINS
               NK(N) = NKout(N)
               DO JC = 1, ICOMP
                  MK(N,JC) = MKout(N,JC)
               ENDDO
            ENDDO

         
         
         ENDIF

                                ! nitrogen and sulfur mass checks
                                ! get the total mass of N
         tot_n_1a = Gc(srtnh4)*14.d0/17.d0
         do k=1,ibins
            tot_n_1a = tot_n_1a + Mk(k,srtnh4)*14.d0/18.d0
         enddo
               
                                ! get the total mass of S
         tot_s_1a = 0.d0
         do k=1,ibins
            tot_s_1a = tot_s_1a + Mk(k,srtso4)*32.d0/96.d0
         enddo

         CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd) 
#if defined ( DEBUG ) 
!      print *, 'mnfix in tomas_mod:677'
#endif
         CALL MNFIX( Nk, Mk, ERRORSWITCH )
         IF ( ERRORSWITCH ) THEN
            PRINT *,'Aerophys: MNFIX found error at',I,J,L
            IF( .not. SPINUP(14.0) ) THEN
               CALL ERROR_STOP('AEROPHYS-MNFIX (2)','After cond/nucl')
            ELSE
               PRINT *,'Let error go during spin up'
            ENDIF
         ENDIF

         MPNUM = 5    
         IF ( ND60 > 0 ) THEN
            CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL )
         ENDIF

         !-----------------------------
         ! Coagulation 
         !-----------------------------

            if(printdebug .and. 
     &           i==iob.and.j==job.and.l==lob) ERRORSWITCH =.TRUE. 
C         if (i==iob .and. j==job .and. l==lob ) 
C     &      CALL DEBUGPRINT( Nk, Mk, I, J, L,'Before coagulation' )
         
         IF( COAG )  THEN
            CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd) 
            CALL MULTICOAG( ADT, Nk, Mk, BOXVOL, PRES, TEMPTMS,
     &             errorswitch )

            if ( errorswitch )
     &       CALL DEBUGPRINT( Nk, Mk, I, J, L,'After coagulation' )
C            if (i==iob .and. j==job .and. l==lob ) 
C     &       CALL DEBUGPRINT( Nk, Mk, I, J, L,'After coagulation' )

            MPNUM = 2
            IF ( ND60 > 0 ) THEN 
               CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL )
            ENDIF

         !Fix any inconsistency after coagulation (win, 4/18/06)
         CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd) 
            if(printdebug .and. i==iob.and.j==job.and.l==lob) 
     &           ERRORSWITCH=.true. !4/18/06 win

#if defined ( DEBUG ) 
!      print *, 'mnfix in tomas_mod:719'
#endif
            CALL MNFIX( NK, MK, ERRORSWITCH )

            IF ( ERRORSWITCH ) THEN
               PRINT *,'MNFIX found error at',I,J,L
               IF( .not. SPINUP(14.0) ) THEN
                  CALL ERROR_STOP('AEROPHYS-MNFIX (3)',
     &                 'After COAGULATION'  )
               ELSE
                  PRINT *,'Let error go during spin up'
               ENDIF
            ENDIF

            MPNUM = 5
            IF ( ND60 > 0 ) THEN
               CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL )
            ENDIF
         ENDIF  ! Coagulation


         ! Do water eqm at appropriate times
         CALL EZNH3EQM( Gc, Mk )
         CALL EZWATEREQM ( MK, RHTOMAS )

         !****************************
         ! End of aerosol dynamics
         !****************************

         !Fix any inconsistencies in M/N distribution (because of advection)
         CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd) 

         ! Make sure anything that leaves AEROPHYS is free of any error
         ! This MNFIX call could be temporary (?) or just leave it here and
         ! monitor if the error fixed is significantly large meaning some 
         ! serious problem needs to be investigated 
         if(printdebug .and. i==iob.and.j==job.and.l==lob) 
     &      ERRORSWITCH =.true. 

#if defined ( DEBUG ) 
!      print *, 'mnfix in tomas_mod:758'
#endif
         CALL MNFIX(NK,MK,ERRORSWITCH) 
         IF ( ERRORSWITCH ) THEN
            PRINT *,'End of Aerophys: MNFIX found error at',I,J,L
            IF( .not. SPINUP(14.0) ) THEN
               CALL ERROR_STOP('AEROPHYS-MNFIX (4)',
     &                         'End of microphysics')
            ELSE
               PRINT *,'Let error go during spin up'
            ENDIF
         ENDIF

         ! Accumulate changes by mnfix to diagnostic (win, 9/8/05)
         MPNUM = 5
         IF ( ND60 > 0 ) THEN
            CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL )
         ENDIF

         ! Swap Nk, Mk, and Gc arrays back to STT
         DO N = 1, IBINS
            TRACNUM = IDTNK1 - 1 + N
            STT(I,J,L,TRACNUM) = NK(N)
            DO JC = 1, ICOMP-IDIAG
               TRACNUM = IDTNK1 - 1 + N + IBINS * JC
               STT(I,J,L,TRACNUM) = MK(N,JC)
            ENDDO
            STT(I,J,L,IDTAW1-1+N) = MK(N,SRTH2O)
         ENDDO
         STT(I,J,L,IDTH2SO4) = GC(SRTSO4)
         
         ! print to file to check mass conserv
c         write(*,77) I,J,L, STT(I,J,L,IDTNH3), 
c     &        STT(I,J,L,IDTNH3)-GC(SRTNH4)

         ! Calculate NH3 gas lost to aerosol phase as NH4 
         NH3_to_NH4 = STT(I,J,L,IDTNH3)-GC(SRTNH4)

         ! Update the bulk NH4 aerosol tracer 
         if ( NH3_to_NH4 > 0d0 ) 
     &        STT(I,J,L,IDTNH4) = STT(I,J,L,IDTNH4) + 
     &                            NH3_to_NH4/17.d0*18.d0

         ! Update NH3 gas tracer (win, 10/6/08)
         ! plus tiny amount CEPS in case zero causes some problem
         STT(I,J,L,IDTNH3)   = GC(SRTNH4) + CEPS !MUST CHECK THIS!! (win,9/26/08)


cvbn         write(889,89)I,J,L,STT(i,j,l,IDTH2SO4)
 89      format(3I3,'STT(IDTH2SO4) kg', E13.5)

      !End of loop over grid cells
      ENDDO                     !L loop
      ENDDO                     !J loop
      ENDDO                     !I loop

!$OMP END PARALLEL DO

!      WRITE(777,*) '---------------------------'
 77   FORMAT(3I4, '  STT(IDTNH3),'E13.5,'  Used', E13.5 ) 

      IF ( COND .and. LPRT ) PRINT *,'### AEROPHYS: SO4 CONDENSATION' 
      IF ( COAG .and. LPRT ) PRINT *,'### AEROPHYS: COAGULATION' 
      IF ( NUCL .and. LPRT ) PRINT *,'### AEROPHYS: NUCLEATION' 
      
      ! Free pointer memory
      NULLIFY( STT )

      END SUBROUTINE AEROPHYS

!------------------------------------------------------------------------------

      SUBROUTINE COND_NUC(Nki,Mki,Gci,Nkf,Mkf,Gcf,fnavg,fn1avg,
     &      H2SO4rate,dti,num_iter,Nknuc,Mknuc,Nkcond,Mkcond,
     &      ionrate, surf_area, BOXVOL, BOXMASS, TEMPTMS, PRES, 
     &      RHTOMAS, errswitch, lev)
!
!******************************************************************************
C     This subroutine calculates the change in the aerosol size distribution
C     due to so4 condensation and binary/ternary nucleation during the
C     overal microphysics timestep.
C     WRITTEN BY Jeff Pierce, May 2007 for GISS GCM-II'
!     Put in GEOS-Chem by Win T. 9/30/08
!
C     Initial values of
C     =================

C     Nki(ibins) - number of particles per size bin in grid cell
C     Nnuci - number of nucleation size particles per size bin in grid cell
C     Mnuci - mass of given species in nucleation pseudo-bin (kg/grid cell)
C     Mki(ibins, icomp) - mass of a given species per size bin/grid cell
C     Gci(icomp-1) - amount (kg/grid cell) of all species present in the
C                    gas phase except water
C     H2SO4rate - rate of H2SO4 chemical production [kg s^-1]
C     dt - total model time step to be taken (s)

C     OUTPUTS

C     Nkf, Mkf, Gcf - same as above, but final values
C     Nknuc, Mknuc - same as above, final values from just nucleation
C     Nkcond, Mkcond - same as above, but final values from just condensation
C     fn, fn1
!
! !REVISION HISTORY:
!  28 Jan 2014 - R. Yantosca - Avoid div-by-zero errors
!*****************************************************************************
!
      IMPLICIT NONE

      REAL*8 Nki(ibins), Mki(ibins, icomp), Gci(icomp-1)
      REAL*8 Nkf(ibins), Mkf(ibins, icomp), Gcf(icomp-1)
      REAL*8 Nknuc(ibins), Mknuc(ibins, icomp)
      REAL*8 Nkcond(ibins),Mkcond(ibins,icomp)
      double precision H2SO4rate
      real             dti
      double precision fnavg        ! nucleation rate of clusters cm-3 s-1
      double precision fn1avg       ! formation rate of particles to first size bin cm-3 s-1
      REAL*4           BOXVOL, BOXMASS, TEMPTMS, RHTOMAS, PRES
      logical          errswitch    ! signal for error 
      integer          lev          ! layer of the model
      REAL*8   surf_area
      REAL*8   ionrate

C-----VARIABLE DECLARATIONS---------------------------------------------

      double precision dt
      integer          i,j,k,c      ! counters
      double precision fn           ! nucleation rate of clusters cm-3 s-1
      double precision fn1          ! formation rate of particles to first size bin cm-3 s-1
      double precision pi, R        ! pi and gas constant (J/mol K)
      double precision CSi,CSa      ! intial and average condensation sinks
      double precision CS1,CS2      ! guesses for condensation sink [s^-1]
      double precision CStest       ! guess for condensation sink
      REAL*8           Nk1(ibins), Mk1(ibins, icomp), Gc1(icomp-1)
      REAL*8           Nk2(ibins), Mk2(ibins, icomp), Gc2(icomp-1)
      REAL*8           Nk3(ibins), Mk3(ibins, icomp), Gc3(icomp-1)
      logical          nflg         ! returned from nucleation, says whether nucleation occurred or not
      double precision mcond,mcond1 ! mass to condense [kg]
      double precision tol          ! tolerance
      double precision eps          ! small number
      double precision sinkfrac(ibins) ! fraction of condensation sink coming from bin k
      double precision totmass      ! the total mass of H2SO4 generated during the timestep
      double precision tmass
      double precision CSch         ! fractional change in condensation sink
      double precision CSch_tol     ! tolerance in change in condensation sink
      double precision addt         ! adaptive timestep time
      double precision time_rem     ! time remaining
      integer          num_iter     ! number of iteration
      double precision sumH2SO4     ! used for finding average H2SO4 conc over timestep
      integer          iter         ! number of iteration
      double precision rnuc         ! critical radius [nm]
      double precision gasConc      ! gas concentration [kg]
      double precision mass_change  ! change in mass during nucleation.f
      double precision total_nh4_1,total_nh4_2
      double precision min_tstep    ! minimum timestep [s]
      integer          nuc_bin      ! the nucleation bin
      double precision sumfn, sumfn1 ! used for getting average nucleation rates
      logical          tempvar,  pdbg
      real*8           tnumb


C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.141592654, R=8.314) !pi and gas constant (J/mol K)
      parameter(eps=1E-40)
      parameter(CSch_tol=0.01)
      parameter(min_tstep=1.0d0)

C-----CODE--------------------------------------------------------------

      pdbg = errswitch ! transfer the signal to print debug from outside
      errswitch = .false. ! flag error to outide to terminate program. Initialize wiht .false.
!      pdbg = .true. ! transfer the signal to print debug from outside

      dt = dble(dti)

C Initialize values of Nkf, Mkf, Gcf, and time
      do j=1,icomp-1
         Gc1(j)=Gci(j)
!         Gcf(j)=Gci(j) !sfarina
      enddo
      do k=1,ibins
         Nk1(k)=Nki(k)
         Nknuc(k)=Nki(k)
         Nkcond(k)=Nki(k)
         do j=1,icomp
            Mk1(k,j)=Mki(k,j)
            Mknuc(k,j)=Mki(k,j)
            Mkcond(k,j)=Mki(k,j)
         enddo
      enddo

C     Get initial condensation sink
      CS1 = 0.d0
      call getCondSink(Nk1,Mk1,srtso4,CS1,sinkfrac,surf_area,
     &               BOXVOL,TEMPTMS, PRES)  
      if( pdbg) print*,'CS1', CS1
c      CS1 = max(CS1,eps)

C     Get initial H2SO4 concentration guess (assuming no nucleation)
C     Make sure that H2SO4 concentration doesn't exceed the amount generated
C     during that timestep (this will happen when the condensation sink is very low)

C     get the steady state H2SO4 concentration
      call getH2SO4conc(Nk1, Mk1, H2SO4rate, CS1, Gc1(srtnh4),
     &                  gasConc, ionrate, surf_area, 
     &                  BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, lev)
      if( pdbg) print*,'gasConc',gasConc
      Gc1(srtso4) = gasConc
      addt = min_tstep
c      addt = 3600.d0
      totmass = H2SO4rate*addt*96.d0/98.d0

      tempvar = pdbg

C     Get change size distribution due to nucleation with initial guess
c      call nucleation(Nk1,Mk1,Gc1,Nnuc1,Mnuc1,totmass,addt,Nk2,Mk2,Gc2,
c     &     Nnuc2,Mnuc2,nflg)    
!      call nucleation(Nk1,Mk1,Gc1,Nk2,Mk2,Gc2,fn,fn1,totmass,nuc_bin,
!     &     addt, PDBG, lev)
      call nucleation(Nk1,Mk1,Gc1,Nk2,Mk2,Gc2,fn,fn1,totmass,nuc_bin,
     &     addt, ionrate, surf_area, BOXVOL, BOXMASS, TEMPTMS,
     &     PRES, RHTOMAS, PDBG, lev)
      
      if(pdbg) then
         print*,'COND_NUC: Found an error at nucleation --> TERMINATE'
         errswitch = .true.
         return
      endif
      pdbg = tempvar !put the print debug switch back to pdbg
      if(pdbg) call debugprint(Nk2, Mk2, 0,0,0,'After nucleation[1]')

c      print*,'after nucleation'
c      print*,'Nnuc1',Nnuc1
c      print*,'Nnuc2',Nnuc2
c      print*,'Mnuc1',Mnuc1
c      print*,'Mnuc2',Mnuc2

      mass_change = 0.d0

      do k=1,ibins
         mass_change = mass_change + (Mk2(k,srtso4)-Mk1(k,srtso4))
      enddo
      if( pdbg)  print*,'mass_change',mass_change

      mcond = totmass-mass_change ! mass of h2so4 to condense

      if( pdbg) print*,'after nucleation'
      if( pdbg)  print*,'totmass',totmass,'mass_change1',mass_change,
     &     'mcond',mcond
      if( pdbg)  print*,'cs1',CS1, Gc1(srtso4)

      if (mcond.lt.0.d0)then
         tmass = 0.d0
         do k=1,ibins
            do j=1,icomp-idiag
               tmass = tmass + Mk2(k,j)
            enddo
         enddo
c     if (abs(mcond).gt.tmass*1.0D-8) then
         if (abs(mcond).gt.totmass*1.0d-8) then
            if (-mcond.lt.Mk2(nuc_bin,srtso4)) then
c               if (CS1.gt.1.0D-5)then
c                  print*,'budget fudge 1 in cond_nuc.f'
c               endif
               tmass = 0.d0
               do j=1,icomp-idiag
                  tmass = tmass + Mk2(nuc_bin,j)
               enddo
               Nk2(nuc_bin) = Nk2(nuc_bin)*(tmass+mcond)/tmass
               Mk2(nuc_bin,srtso4) = Mk2(nuc_bin,srtso4) + mcond
               mcond = 0.d0
            else
               print*,'budget fudge 2 in cond_nuc.f'
               do k=2,ibins
                  Nk2(k) = Nk1(k)
                  Mk2(k,srtso4) = Mk1(k,srtso4)
               enddo
               Nk2(1) = Nk1(1)+totmass/sqrt(xk(1)*xk(2))
               Mk2(1,srtso4) = Mk1(1,srtso4) + totmass
               mcond = 0.d0        
c               print*,'mcond < 0 in cond_nuc', mcond, totmass
c               stop
            endif
         else
            mcond = 0.d0
         endif
      endif
      
c      if (mcond.lt.0.d0)then
c         print*,'mcond < 0 in cond_nuc', mcond
c         stop
c      endif
      tmass = 0.d0
      do k=1,ibins
         do j=1,icomp-idiag
            tmass = tmass + Mk2(k,j)
         enddo
      enddo
      if( pdbg)  print*, 'mcond',mcond,'tmass',tmass,'nuc',Nk2(1)-Nk1(1)
      tempvar = pdbg
      
C     Get guess for condensation
      call ezcond(Nk2,Mk2,mcond,srtso4,Nk3,Mk3,surf_area,
     &      BOXVOL, TEMPTMS, PRES, pdbg )

      if(pdbg) then
         print*,'COND_NUC: Found an error at EZCOND --> TERMINATE'
         errswitch = .true.
         return
      endif
      pdbg = tempvar
      if(pdbg) call debugprint(Nk3, Mk3, 0,0,0,'After EZCOND[1]')
c      print*,'after ezcond',Nk2,Nk3
Cjrp      mcond1 = 0.d0
Cjrp      do k=1,ibins
Cjrp         do j=1,icomp
Cjrp            mcond1 = mcond1 + (Mk3(k,j)-Mk2(k,j))
Cjrp         enddo
Cjrp      enddo
c      print*,'mcond',mcond,'mcond1',mcond1

      Gc3(srtnh4) = Gc1(srtnh4)   

      call eznh3eqm(Gc3,Mk3)
      call ezwatereqm(Mk3, RHTOMAS)

      ! check to see how much condensation sink changed
      call getCondSink(Nk3,Mk3,srtso4,CS2,sinkfrac,surf_area,
     &               BOXVOL,TEMPTMS, PRES)  
      CSch = abs(CS2 - CS1)/CS1    
       
c      if (CSch.gt.CSch_tol) then ! condensation sink didn't change much use whole timesteps
         ! get starting adaptive timestep to not allow condensationk sink
         ! to change that much
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% Prior to 1/28/14:
!%%% BUG FIX: Do not allow a division by zero to occur.  Set ADDT to zero
!%%% if CSCH is zero.  ADDT is reset to MIN_TSTEP below so that should be OK.
!%%% (bmy, 1/28/14)
!%%%      addt = addt*CSch_tol/CSch/2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Avoid div-by-zero (bmy, 1/28/14)
         IF ( ABS( CSch ) > 0d0 ) THEN
            addt = addt*CSch_tol/CSch/2d0
         ELSE
            addt = 0d0
         ENDIF
         addt = min(addt,dt)
         addt = max(addt,min_tstep)

                
         time_rem = dt ! time remaining
         if( pdbg)    print*,'addt',addt,time_rem
         num_iter = 0
         sumH2SO4=0.d0
         sumfn = 0.d0
         sumfn1 = 0.d0
         ! do adaptive timesteps
         do while (time_rem .gt. 0.d0)
            num_iter = num_iter + 1
            if( pdbg)    print*, 'iter', num_iter, ' addt', addt,
     &           'time_rem', time_rem
C     get the steady state H2SO4 concentration
            if (num_iter.gt.1)then ! no need to recalculate for first step
      call getH2SO4conc(Nk1, Mk1, H2SO4rate, CS1, Gc1(srtnh4),
     &                  gasConc, ionrate, surf_area,
     &                  BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, lev)
               Gc1(srtso4) = gasConc
            endif
            if( pdbg)    print*,'gasConc',gasConc

            sumH2SO4 = sumH2SO4 + Gc1(srtso4)*addt
            totmass = H2SO4rate*addt*96.d0/98.d0
c            call nucleation(Nk1,Mk1,Gc1,Nnuc1,Mnuc1,totmass,addt,Nk2,
c     &           Mk2,Gc2,Nnuc2,Mnuc2,nflg,lev) 

            !Debug to see what goes in nucleation (win, 10/3/08)
            if(pdbg) then 
               print*,'Temperature',TEMPTMS,'RH',RHTOMAS
               print*,'H2SO4',Gc1(srtso4)/boxvol*1000.d0/98.d0*6.022d23
               print*,'NH3ppt',Gc1(srtnh4)/17.d0/(boxmass/29.d0)*1d12
            endif


            tempvar = pdbg
            call nucleation(Nk1,Mk1,Gc1,Nk2,Mk2,Gc2,fn,fn1,totmass,
     &           nuc_bin,addt, ionrate, surf_area, BOXVOL, BOXMASS,
     &           TEMPTMS, PRES, RHTOMAS, PDBG, lev)

            if(pdbg) then
               print*,'COND_NUC: Error at nucleation[2] --> TERMINATE'
               errswitch=.true.
               return
            endif
            pdbg = tempvar
            if(pdbg) call debugprint(Nk2, Mk2, 0,0,0,
     &                               'After nucleation[2]')
c            print*,'after nucleation iter'
            sumfn = sumfn + fn*addt
            sumfn1 = sumfn1 + fn1*addt
            
c            total_nh4_1 = Mnuc1(srtnh4)
c            total_nh4_2 = Mnuc2(srtnh4)
c            do i=1,ibins
c               total_nh4_1 = total_nh4_1 + Mk1(i,srtnh4)
c               total_nh4_2 = total_nh4_2 + Mk2(i,srtnh4)
c            enddo
c            print*,'total_nh4',total_nh4_1,total_nh4_2

            mass_change = 0.d0

            do k=1,ibins
               mass_change = mass_change + (Mk2(k,srtso4)-Mk1(k,srtso4))
            enddo 
            if( pdbg)    print*,'mass_change2',mass_change

            mcond = totmass-mass_change ! mass of h2so4 to condense

c            print*,'after nucleation'
c            print*,'totmass',totmass,'mass_change',mass_change,
c     &           'mcond',mcond

c      print*,'2 mass_change',mass_change,mcond,totmass
c      print*,'2 cs1',CS1, Gc1(srtso4)

            if (mcond.lt.0.d0)then
               tmass = 0.d0
               do k=1,ibins
                  do j=1,icomp-idiag
                     tmass = tmass + Mk2(k,j)
                  enddo
               enddo
c     if (abs(mcond).gt.tmass*1.0D-8) then
               if (abs(mcond).gt.totmass*1.0D-8) then
                  if (-mcond.lt.Mk2(nuc_bin,srtso4)) then
c                     if (CS1.gt.1.0D-5)then
c                        print*,'budget fudge 1 in cond_nuc.f'
c                     endif
                     tmass = 0.d0
                     do j=1,icomp-idiag
                        tmass = tmass + Mk2(nuc_bin,j)
                     enddo
                     Nk2(nuc_bin) = Nk2(nuc_bin)*(tmass+mcond)/tmass
                     Mk2(nuc_bin,srtso4) = Mk2(nuc_bin,srtso4) + mcond
                     mcond = 0.d0
                  else
                     print*,'budget fudge 2 in cond_nuc.f'
                     do k=2,ibins
                        Nk2(k) = Nk1(k)
                        Mk2(k,srtso4) = Mk1(k,srtso4)
                     enddo
                     Nk2(1) = Nk1(1)+totmass/sqrt(xk(1)*xk(2))
                     Mk2(1,srtso4) = Mk1(1,srtso4) + totmass
                     print*,'mcond < 0 in cond_nuc', mcond, totmass
                     mcond = 0.d0 
! should I stop or not?? (win, 10/4/08) 
                    !stop
                     ! change from stop here to stop outside with more info (win, 10/4/08)
                     print*,'COND_NUC: --> TERMINATE'
!10/4/08                     errswitch = .true.
!10/4/08                     return
                  endif
               else
                  mcond = 0.d0
               endif
            endif

            do k=1,ibins
               Nknuc(k) = Nknuc(k)+Nk2(k)-Nk1(k)
               do j=1,icomp-idiag
                  Mknuc(k,j)=Mknuc(k,j)+Mk2(k,j)-Mk1(k,j)
               enddo
            enddo
            
            
c            Gc2(srtnh4) = Gc1(srtnh4)
c            call eznh3eqm(Gc2,Mk2,Mnuc2)
c            call ezwatereqm(Mk2,Mnuc2)

c            call getCondSink(Nk2,Mk2,Nnuc2,Mnuc2,srtso4,CStest,sinkfrac)  

            ! Before entering ezcond, check if there's enough aerosol to 
            ! condense onto. After several iteration in the case with high
            ! H2SO4 amount but little existing aerosol and also lack the conditions
            ! for nucleation, the whole size distribution is grown out of our 
            ! tracked size bins, so let's exit the loop if there is no aerosol
            ! to condense onto anymore. (win, 10/4/08)
            tmass = 0.d0
            tnumb = 0.d0
            do k=1,ibins
               tnumb = tnumb + Nk2(k)
               do j=1,icomp-idiag
                  tmass = tmass + Mk2(k,j)
               enddo
            enddo
            
            if( (tmass+mcond)/tnumb  > Xk(ibins) ) then
               if( .not. SPINUP(10.0) ) then
                  print*,'Not enough aerosol for condensation!'
                  print*,'  Exiting COND_NUC iteration with '
                  print*,time_rem,'sec remaining time'
               endif

               Gc3(srtnh4)=Gc2(srtnh4)
               do k=1,ibins
                  Nk3(k)=Nk2(k)
                  do j=1,icomp
                     Mk3(k,j)=Mk2(k,j)
                  enddo
               enddo    
               goto 100
            endif     
            
            tempvar = pdbg

            call ezcond(Nk2,Mk2,mcond,srtso4,Nk3,Mk3,surf_area,
     &      BOXVOL, TEMPTMS, PRES, pdbg)
            do k=1,ibins
               Nkcond(k) = Nkcond(k)+Nk3(k)-Nk2(k)
               do j=1,icomp-idiag
                  Mkcond(k,j)=Mkcond(k,j)+Mk3(k,j)-Mk2(k,j)
               enddo
            enddo
            Gc3(srtnh4) = Gc1(srtnh4)

            if(pdbg) then
               print*,'COND_NUC: Error at EZCOND[2] --> TERMINATE'
               errswitch=.true.
               return
            endif
            pdbg = tempvar

            if(pdbg) call debugprint(Nk3, Mk3, 0,0,0,'After EZCOND[2]')

            if( pdbg)    print*,'after ezcond iter'
            call eznh3eqm(Gc3,Mk3)
            call ezwatereqm(Mk3, RHTOMAS)
            
                                ! check to see how much condensation sink changed
            call getCondSink(Nk3,Mk3,srtso4,CS2,sinkfrac,surf_area,
     &               BOXVOL,TEMPTMS, PRES)  

            time_rem = time_rem - addt
            if (time_rem .gt. 0.d0) then
               CSch = abs(CS2 - CS1)/CS1
Cjrp               if (CSch.lt.0.d0) then
Cjrp                  print*,''
Cjrp                  print*,'CSch LESS THAN ZERO!!!!!', CS1,CStest,CS2
Cjrp                  print*,'Nnuc',Nnuc1,Nnuc2
Cjrp                  print*,''
Cjrp
Cjrp                  addt = min(addt,time_rem)
Cjrp               else

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% Prior to 1/28/14:
!%%% BUG FIX: Do not allow a division by zero to occur.  Set ADDT to zero
!%%% if CSCH is zero.  ADDT is reset to MIN_TSTEP below so that should be OK.
!%%% (bmy, 1/28/14)
!%%%               addt = min(addt*CSch_tol/CSch,addt*1.5d0) ! allow adaptive timestep to change
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               ! Allow adaptive timestep to change
               ! Avoid div-by-zero error
               IF ( ABS( CSch ) > 0d0 ) THEN 
                  addt = min(addt*CSch_tol/CSch,addt*1.5d0) 
               ELSE
                  addt = 0d0
               ENDIF

               ! allow adaptive timestep to change again
               addt = min(addt,time_rem) 
               addt = max(addt,min_tstep)
Cjrp               endif
               if( pdbg)     print*,'CS1',CS1,'CS2',CS2
               CS1 = CS2
               Gc1(srtnh4)=Gc3(srtnh4)
               do k=1,ibins
                  Nk1(k)=Nk3(k)
                  do j=1,icomp
                     Mk1(k,j)=Mk3(k,j)
                  enddo
               enddo         
            endif
         enddo ! while loop

 100     continue  

         Gcf(srtso4)=sumH2SO4/dt
         fnavg = sumfn/dt
         fn1avg = sumfn1/dt
         if( pdbg)    print*,'AVERAGE GAS CONC',Gcf(srtso4)

Cjrp      else
Cjrp         num_iter = 1
Cjrp         Gcf(srtso4)=Gc1(srtso4)
Cjrp      endif
      
         if( pdbg)    print*, 'cond_nuc num_iter =', num_iter
!T0M(1,1,1,3) = double(num_iter) ! store iterations here

         if(pdbg) call debugprint(Nk3, Mk3, 0,0,0,'End of COND_NUC')

      do k=1,ibins
         Nkf(k)=Nk3(k)
         do j=1,icomp
            Mkf(k,j)=Mk3(k,j)
         enddo
      enddo      
      Gcf(srtnh4)=Gc3(srtnh4)

      return
      END SUBROUTINE COND_NUC

!------------------------------------------------------------------------------

      SUBROUTINE getCondSink(Nko, Mko, spec, CS, sinkfrac, surf_area,
     &       BOXVOL, TEMPTMS, PRES)
!
!******************************************************************************
C     This subroutine calculates the condensation sink (first order loss
C     rate of condensing gases) from the aerosol size distribution.
C     WRITTEN BY Jeff Pierce, May 2007 for GISS GCM-II' 
!     Put in GEOS-Chem by Win T. (9/30/08)
!
C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Nk(ibins) - number of particles per size bin in grid cell
C     Nnuc - number of particles per size bin in grid cell
C     Mnuc - mass of given species in nucleation pseudo-bin (kg/grid cell)
C     Mk(ibins, icomp) - mass of a given species per size bin/grid cell
C     spec - number of the species we are finding the condensation sink for

C-----OUTPUTS-----------------------------------------------------------

C     CS - condensation sink [s^-1]
C     sinkfrac(ibins) - fraction of condensation sink from a bin
!
!******************************************************************************

      IMPLICIT NONE

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision Nko(ibins), Mko(ibins, icomp)
      double precision CS, sinkfrac(ibins)
      REAL*4, INTENT(IN)       :: BOXVOL, TEMPTMS, PRES
      integer spec
      REAL*8, INTENT(OUT)     :: surf_area

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,k,c           ! counters
      double precision pi, R    ! pi and gas constant (J/mol K)
      double precision mu                  !viscosity of air (kg/m s)
      double precision mfp                 !mean free path of air molecule (m)
      double precision l_ab                !mean free path of h2so4 molecule (m)
      real Di       !diffusivity of gas in air (m2/s)
      double precision Neps     !tolerance for number
      real density  !density [kg m^-3]
      double precision mp       !mass per particle [kg]
      double precision Dpk(ibins) !diameter of particle [m]
      double precision Kn       !Knudson number
      double precision beta(ibins) !non-continuum correction factor
      double precision Mktot    !total mass in bin [kg]
      double precision c_a      !average speed of a, h2so4 molecule

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.141592654, R=8.314) !pi and gas constant (J/mol K)
      parameter(Neps=1.0d10)
      double precision alpha(icomp) ! accomodation coef  
!      data alpha/0.65,0.,0.,0.,0.,0.,0.,0.,0./
      real Sv(icomp)         !parameter used for estimating diffusivity
!      data Sv /42.88,42.88,42.88,42.88,42.88,42.88,42.88,
!     &         42.88,42.88/

C-----CODE--------------------------------------------------------------

      ! have to find a better way to simply assign contants to these array
      ! The problem is I declare the array with ICOMP - its value will be 
      ! determined at time of run, so I can't use DATA statement
      DO J=1,ICOMP
!         IF ( J == SRTSO4 ) THEN
            alpha(J) = 0.65
!         ELSE
!            alpha(J) = 0.
!         ENDIF
         Sv(J) = 42.88
      ENDDO


C     get some parameters

!      mu=2.5277e-7 * TEMPTMS**0.75302
!      mfp=2.0*mu / ( pres*sqrt( 8.0 * 0.6589 / (pi*R*TEMPTMS) ) )  !S&P eqn 8.6

!      mfp=2.0*mu / ( pres*sqrt( 8.0 * 0.0289 / (pi*R*TEMPTMS) ) )  !S&P eqn 8.6

      Di=gasdiff(TEMPTMS,pres,98.0,Sv(spec))

      c_a  = sqrt(8.0 * TEMPTMS * R / 0.098)
      l_ab = 2.0 * Di / c_a

C     get size dependent values
      do k=1,ibins
         if (Nko(k) .gt. Neps) then
            Mktot=0.d0
            do j=1,icomp
                  Mktot=Mktot+Mko(k,j)
            enddo
Ckpc  Density should be changed due to more species involed.     
            density=aerodens(Mko(k,srtso4),0.d0,
     &           Mko(k,srtnh4),Mko(k,srtnacl),Mko(k,srtecil),
     &           Mko(k,srtecob),Mko(k,srtocil),Mko(k,srtocob),
     &           Mko(k,srtdust),Mko(k,srth2o)) !assume bisulfate                     
            mp=Mktot/Nko(k)
         else
            !nothing in this bin - set to "typical value"
            density=1500.
#if  defined( TOMAS12 ) || defined( TOMAS15 )
            mp=sqrt(xk(k)*xk(k+1))
#else
            mp=1.4*xk(k)
#endif
         endif
         Dpk(k)  = ( (mp/density)*(6./pi) )**(0.333)
!        Kn      = 2.0 * mfp  / Dpk(k)                             !S&P eqn 11.35 (text)
         Kn      = 2.0 * l_ab / Dpk(k)                             !S&Pv2 chapter 12 - Kn for Dahneke correction factor
         beta(k) = ( 1.+Kn )  / ( 1.+2.*Kn*(1.+Kn)/alpha(spec) )    !S&P eqn 11.35
      enddo      

C     get condensation sink
      CS = 0.d0
      surf_area = 0.d0
      do k=1,ibins
         CS = CS + Dpk(k)*Nko(k)*beta(k)
         surf_area = surf_area+Nko(k)*pi*(Dpk(k)*1.0D6)**2
      enddo
      do k=1,ibins
         sinkfrac(k) = Dpk(k)*Nko(k)*beta(k)/CS
      enddo
      CS = 2.d0*pi*dble(Di)*CS/(dble(boxvol)*1D-6)
      surf_area = surf_area/(dble(boxvol))

      return
      end subroutine getcondsink

!------------------------------------------------------------------------------

      SUBROUTINE getH2SO4conc(Nk, Mk, H2SO4rate, CS, NH3conc, gasConc, 
     &    ionrate, surf_area, BOXVOL, BOXMASS,
     &    TEMPTMS, PRES, RHTOMAS, lev)
!
!******************************************************************************
C     This subroutine uses newtons method to solve for the steady state 
C     H2SO4 concentration when nucleation is occuring.

C     It solves for H2SO4 in 0 = P - CS*H2SO4 - M(H2SO4)

C     where P is the production rate of H2SO4, CS is the condensation sink
C     and M(H2SO4) is the loss of mass towards making new particles.
C     WRITTEN BY Jeff Pierce, May 2007 for GISS GCM-II'
!     Put in GEOS-CHEM by Win T. (9/30/08)
!C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     H2SO4rate - H2SO4 generation rate [kg box-1 s-1]
C     CS - condensation sink [s-1]
C     NH3conc - ammonium in box [kg box-1]
Cxxxxx     prev - logical flag saying if a previous guess should be used or not
Cxxxxx     gasConc_prev - the previous guess [kg/box] (not used if prev is false)

C-----OUTPUTS-----------------------------------------------------------

C     gasConc - gas H2SO4 [kg/box]
!*****************************************************************************
!
      ! Reference to F90 modules
      USE ERROR_MOD,      ONLY : ERROR_STOP, IT_IS_NAN

      IMPLICIT NONE
C-----ARGUMENT DECLARATIONS---------------------------------------------

      REAL*8                         :: Nk(IBINS)
      REAL*8                         :: Mk(IBINS, ICOMP)
      double precision H2SO4rate
      double precision CS
      double precision NH3conc
      double precision gasConc
      logical          prev
      double precision gasConc_prev
      REAL*4, INTENT(IN)       :: BOXVOL,  BOXMASS, TEMPTMS
      REAL*4, INTENT(IN)       :: PRES,    RHTOMAS
      integer          lev

      REAL*8                   :: ionrate
      REAL*8                   :: surf_area

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,k,c           ! counters
      double precision fn, rnuc ! nucleation rate [# cm-3 s-1] and critical radius [nm]
      double precision mnuc, mnuc1 ! mass of nucleated particle [kg]
      double precision fn1, rnuc1 ! nucleation rate [# cm-3 s-1] and critical radius [nm]
      double precision res      ! d[H2SO4]/dt, need to find the solution where res = 0
      double precision massnuc     ! mass being removed by nucleation [kg s-1 box-1]
      double precision gasConc1 ! perturbed gasConc
      double precision gasConc_hi, gasConc_lo
      double precision res1     ! perturbed res
      double precision res_new  ! new guess for res
      double precision dresdgasConc ! derivative for newtons method
      double precision Gci(icomp-1)      !array to carry gas concentrations
      logical nflg              !says if nucleation occured
      double precision H2SO4min !minimum H2SO4 concentration in parameterizations (molec/cm3)
      double precision pi
      integer iter,iter1
      double precision CSeps    ! low limit for CS
      double precision max_H2SO4conc !maximum H2SO4 concentration in parameterizations (kg/box)
      double precision nh3ppt   !ammonia concentration in ppt
 
C     VARIABLE COMMENTS...

C-----EXTERNAL FUNCTIONS------------------------------------------------


C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.141592654)
!      parameter(H2SO4min=1.D4) !molecules cm-3
      parameter(CSeps=1.0d-20)
        
C-----CODE--------------------------------------------------------------

      do i=1,icomp-1
         Gci(i)=0.d0
      enddo
      Gci(srtnh4)=NH3conc 
      
      ! make sure CS doesn't equal zero
c     CS = max(CS,CSeps)
      
      ! some specific stuff for napari vs. vehk
      if (ion_nuc.eq.1) then
         H2SO4min=1.0D5
      elseif (ion_nuc.eq.2) then
         H2SO4min=5.0D5
      else
         H2SO4min=1.0D4
      endif

      if ((bin_nuc.eq.1).or.(tern_nuc.eq.1).or.(ion_nuc.le.2))then
         nh3ppt = Gci(srtnh4)/17.d0/(boxmass/29.d0)*1d12*
     &             PRES/101325.*273./TEMPTMS ! corrected for pressure (because this should be concentration)
         if (ion_nuc.eq.1)then
            max_H2SO4conc=1.0D8*boxvol/1000.d0*98.d0/6.022d23
         elseif (ion_nuc.eq.2)then
            max_H2SO4conc=5.0D8*boxvol/1000.d0*98.d0/6.022d23
         elseif ((nh3ppt.gt.1.0d0).and.(tern_nuc.eq.1))then
            max_H2SO4conc=1.0D9*boxvol/1000.d0*98.d0/6.022d23
         elseif (bin_nuc.eq.1)then
            max_H2SO4conc=1.0D11*boxvol/1000.d0*98.d0/6.022d23
         else
            max_H2SO4conc = 1.0D100
         endif
      else
         max_H2SO4conc = 1.0D100
      endif

      
C     Checks for when condensation sink is very small
      if (CS.gt.CSeps) then
         gasConc = H2SO4rate/CS
      else
         if((bin_nuc.gt.0).or.(tern_nuc.gt.0).or.
     &                           (ion_nuc.gt.0))then
            gasConc = max_H2SO4conc
         else
            print*,'condesation sink too small in getH2SO4conc.f'
            STOP
         endif
      endif

      gasConc = min(gasConc,max_H2SO4conc)
      Gci(srtso4) = gasConc
      call getNucRate(Nk, Mk, Gci,fn,mnuc,nflg,ionrate, surf_area,
     &       BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, lev)

      if (fn.gt.0.d0) then      ! nucleation occured
         gasConc_lo = H2SO4min*boxvol/(1000.d0/98.d0*6.022d23) !convert to kg/box
         
C     Test to see if gasConc_lo gives a res < 0 (this means ANY nucleation is too high)
         Gci(srtso4) = gasConc_lo*1.000001d0
         call getNucRate(Nk,Mk,Gci,fn1,mnuc1,nflg,ionrate,surf_area,
     &        BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, lev)
         if (fn1.gt.0.d0) then
            massnuc = mnuc1*fn1*boxvol*98.d0/96.d0
c     massnuc = 4.d0/3.d0*pi*(rnuc1*1.d-9)**3*1350.*fn1*boxvol*
c     massnuc = 4.d0/3.d0*pi*(rnuc1*1.d-9)**3*1800.*fn1*boxvol*
c     &           98.d0/96.d0
C     jrp            print*,'res',res
C     jrp            print*,'H2SO4rate',H2SO4rate
C     jrp            print*,'CS*gasConc_lo',CS*gasConc_lo
C     jrp            print*,'mnuc',mnuc
            res = H2SO4rate - CS*gasConc_lo - massnuc
            if (res.lt.0.d0) then ! any nucleation too high
!               if (.not. spinup(14.0)) print*,'nucleation cuttoff'
               gasConc = gasConc_lo*1.000001 ! have nucleation occur and fix mass balance after
               return
            endif
         endif
         
         gasConc_hi = gasConc   ! we know this must be the upper limit (since no nucleation)
                                !take density of nucleated particle to be 1350 kg/m3
         massnuc = mnuc*fn*boxvol*98.d0/96.d0
c         print*,'H2SO4rate',H2SO4rate,'CS*gasConc',CS*gasConc,
c     &        'mnuc',mnuc
         res = H2SO4rate - CS*gasConc - massnuc
         
                                ! check to make sure that we can get solution
         if (res.gt.H2SO4rate*1.d-10) then
            print*,'gas production rate too high in getH2SO4conc.f'
            print*,H2SO4rate,CS,gasConc,massnuc,res
            return
c     STOP
         endif
         
         iter = 0
C     jrp         print*, 'iter',iter
C     jrp         print*,'gasConc_lo',gasConc_lo,'gasConc_hi',gasConc_hi
C     jrp         print*,'res',res
         do while ((abs(res/H2SO4rate).gt.1.D-4).and.(iter.lt.40))
            iter = iter+1
            if (res .lt. 0.d0) then ! H2SO4 concentration too high, must reduce
               gasConc_hi = gasConc ! old guess is new upper bound
            elseif (res .gt. 0.d0) then ! H2SO4 concentration too low, must increase
               gasConc_lo = gasConc ! old guess is new lower bound
            endif
c            print*, 'iter',iter
c            print*,'gasConc_lo',gasConc_lo,'gasConc_hi',gasConc_hi
            gasConc = sqrt(gasConc_hi*gasConc_lo) ! take new guess as logmean
            Gci(srtso4) = gasConc
            call getNucRate(Nk, Mk,Gci,fn,mnuc,nflg,ionrate,surf_area,
     &            BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, lev)
            massnuc = mnuc*fn*boxvol*98.d0/96.d0
            res = H2SO4rate - CS*gasConc - massnuc
c            print*,'res',res
c            print*,'H2SO4rate',H2SO4rate,'CS',CS,'gasConc',gasConc
            if (iter.eq.40.and.CS.gt.1.0D-4)then
               print*,'getH2SO4conc iter break'
               print*,'H2SO4rate',H2SO4rate,'CS',CS
               print*,'gasConc',gasConc,'massnuc',massnuc
               print*,'max_H2SO4conc',max_H2SO4conc
               print*,'fn',fn
               print*,'res/H2SO4rate',res/H2SO4rate
            endif
         enddo
         
c     print*,'IN getH2SO4conc'
c     print*,'fn',fn
c     print*,'H2SO4rate',H2SO4rate
c     print*,'massnuc',massnuc,'CS*gasConc',CS*gasConc
         
      else                      ! nucleation didn't occur
      endif
      
      return
      end SUBROUTINE GETH2SO4CONC

!------------------------------------------------------------------------------

      SUBROUTINE getNucRate(Nk, Mk, Gci,fn,mnuc,nflg, ionrate,surf_area,
     &            BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, lev)
!
!******************************************************************************
C     This subroutine calls the Vehkamaki 2002 and Napari 2002 nucleation
C     parameterizations and gets the binary and ternary nucleation rates.
C     WRITTEN BY Jeff Pierce, April 2007 for GISS GCM-II'
!     Put in GEOS-Chem by win T. 9/30/08
!
C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Gci(icomp-1) - amount (kg/grid cell) of all species present in the
C                    gas phase except water

C-----OUTPUTS-----------------------------------------------------------

C     fn - nucleation rate [# cm-3 s-1]
C     rnuc - radius of nuclei [nm]
C     nflg - says if nucleation happend
!*****************************************************************************
!
      ! Reference to F90 modules
      USE ERROR_MOD,      ONLY : ERROR_STOP, IT_IS_NAN

      IMPLICIT NONE
C-----ARGUMENT DECLARATIONS---------------------------------------------

      REAL*8, INTENT(INOUT)          :: Nk(IBINS)
      REAL*8, INTENT(INOUT)          :: Mk(IBINS, ICOMP)
      REAL*4, INTENT(IN)       :: BOXVOL,  BOXMASS, TEMPTMS
      REAL*4, INTENT(IN)       :: PRES,    RHTOMAS
      REAL*8, INTENT(IN)       :: Gci(icomp-1)

      REAL*8                         :: surf_area
      REAL*8                         :: ionrate

      integer j,i,k
      double precision fn       ! nucleation rate to first bin cm-3 s-1
      double precision mnuc     !mass of nucleating particle [kg]
      logical nflg
      integer lev

C-----VARIABLE DECLARATIONS---------------------------------------------

      double precision nh3ppt   ! gas phase ammonia in pptv
      double precision h2so4    ! gas phase h2so4 in molec cc-1
      double precision gtime    ! time to grow to first size bin [s]
      double precision ltc, ltc1, ltc2 ! coagulation loss rates [s-1]
      double precision Mktot    ! total mass in bin
      double precision neps
      double precision meps
      double precision density  ! density of particle [kg/m3]
      double precision pi
      double precision frac     ! fraction of particles growing into first size bin
      double precision d1,d2    ! diameters of particles [m]
      double precision mp       ! mass of particle [kg]
      double precision mold     ! saved mass in first bin
      double precision rnuc     ! critical nucleation radius [nm]
      double precision sinkfrac(ibins) ! fraction of loss to different size bins
      double precision nadd     ! number to add
      double precision CS       ! kerminan condensation sink [m-2]
      double precision Dpmean   ! the number wet mean diameter of the existing aerosol
      double precision Dp1      ! the wet diameter of bin 1
      double precision dens1    ! density in bin 1 [kg m-3]
      double precision GR       ! growth rate [nm hr-1]
      double precision gamma,eta ! used in kerminen 2004 parameterzation
      double precision drymass,wetmass,WR
      double precision fn_c     ! barrierless nucleation rate
      double precision h1,h2,h3,h4,h5,h6
      double precision dum1,dum2,dum3,dum4   ! dummy variables
      double precision rhin,tempin ! rel hum in

      real*8    mydummy
C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter (neps=1E8, meps=1E-8)
      parameter (pi=3.14159)

C-----CODE--------------------------------------------------------------

      h2so4 = Gci(srtso4)/boxvol*1000.d0/98.d0*6.022d23
      nh3ppt = Gci(srtnh4)/17.d0/(boxmass/29.d0)*1d12*
     &             PRES/101325.*273./TEMPTMS ! corrected for pressure (because this should be concentration)

      fn = 0.d0
      rnuc = 0.d0

c      print*,'h2so4',h2so4,'nh3ppt',nh3ppt

C     if requirements for nucleation are met, call nucleation subroutines
C     and get the nucleation rate and critical cluster size
      if (h2so4.gt.1.d4) then   
         if ((nh3ppt.gt.0.1).and.(tern_nuc.eq.1)) then 
c            print*, 'napari'
            call napa_nucl(TEMPTMS,RHTOMAS,h2so4,nh3ppt,fn,rnuc) !ternary nuc 
            if (ion_nuc.eq.1.and.ionrate.ge.1.d0) then
               call ion_nucl(h2so4,surf_area,TEMPTMS,ionrate,RHTOMAS,
     &              h1,h2,h3,h4,h5,h6)
            else
               h1=0.d0
            endif
            if (h1.gt.fn)then
               fn=h1
               rnuc=h5
            endif
            nflg=.true.
         elseif (bin_nuc.eq.1) then
c            print*, 'vehk'
            call vehk_nucl(TEMPTMS,RHTOMAS,h2so4,fn,rnuc) !binary nuc
            if ((ion_nuc.eq.1).and.(ionrate.ge.1.d0)) then
               call ion_nucl(h2so4,surf_area,TEMPTMS,ionrate,RHTOMAS,
     &              h1,h2,h3,h4,h5,h6)
            else
               h1=0.d0
            endif
            if (h1.gt.fn)then
               fn=h1
               rnuc=h5
            endif
            if (fn.gt.1.0d-6)then
               nflg=.true.
            else
               fn = 0.d0
               nflg=.false.
            endif
         elseif ((ion_nuc.eq.1).and.(ionrate.ge.1.d0)) then
            call ion_nucl(h2so4,surf_area,TEMPTMS,ionrate,RHTOMAS,
     &              h1,h2,h3,h4,h5,h6)
            fn=h1
            rnuc=h5
            nflg=.true.
         elseif(ion_nuc.eq.2) then
            ! Yu Ion nucleation
c            ! first we need to calculate the available surface area
c            surf_area = 0.D0
c            do k=1, ibins
c               if (Nki(k) .gt. Neps) then
c                  Mktot=0.d0
c                  do j=1,icomp
c                     Mktot=Mktot+Mki(k,j)
c                  enddo
c                  mp=Mktot/Nki(k)
c                  density=aerodens(Mki(k,srtso4),0.d0,
c     &                    Mki(k,srtnh4),0.d0,Mki(k,srth2o))  ! assume bisulfate
c                  ! diameter = ((mass/density)*(6/pi))**(1/3)
c                  d2 = 1.D6*((mp/density)*(6.D0/pi))**(1.D0/3.D0) ! (micrometers)
c                  ! surface area per particle = pi*diameter**2
c                  surf_area = surf_area + 1.D-6*(Nki(k)/boxvol)*
c     &                        pi*(d2**2.D0) ! (um2 cm-2)
c               endif
c            enddo
            rhin=dble(RHTOMAS*100.d0)
            tempin=dble(TEMPTMS)
            call YUJIMN(h2so4, rhin, tempin, ionrate, surf_area,
     &               fn, dum1, rnuc, dum2)
            nflg=.true.
         else
            nflg=.false.
         endif         
         if((act_nuc.eq.1).and.(lev.le.7))then
            call bl_nucl(h2so4,fn,rnuc)
            nflg=.true.
         endif
         call cf_nucl(TEMPTMS,RHTOMAS,h2so4,nh3ppt,fn_c) ! use barrierless nucleation as a max for ternary
         fn = min(fn,fn_c)
      else
         nflg=.false.
      endif

      if (fn.gt.0.d0) then
         call getCondSink_kerm(Nk,Mk,CS,Dpmean,Dp1,dens1,
     &         BOXVOL, TEMPTMS, PRES)
         d1 = rnuc*2.d0*1D-9
         drymass = 0.d0
         do j=1,icomp-idiag
            drymass = drymass + Mk(1,j)
         enddo
         wetmass = 0.d0
         do j=1,icomp
            wetmass = wetmass + Mk(1,j)
         enddo
         !prior 10/15/08 
         !WR = wetmass/drymass
         
         ! prevent division by zero (win, 10/15/08)
         if( drymass == 0.d0 ) then
            WR = 1.d0 
         else 
            WR = wetmass/drymass
         endif

cc         print*,'[getnucrate] Gci',Gci
cc         print*,'WR',WR, 'drymass',drymass, 'wetmass',wetmass
         call getGrowthTime(d1,Dp1,Gci(srtso4)*WR,TEMPTMS,
     &        boxvol,dens1,gtime)
         GR = (Dp1-d1)*1D9/gtime*3600.d0 ! growth rate, nm hr-1
         
         gamma = 0.23d0*(d1*1.0d9)**(0.2d0)*(Dp1*1.0d9/3.d0)**0.075d0*
     &        (Dpmean*1.0d9/150.d0)**0.048d0*(dens1*1.0d-3)**
     &        (-0.33d0)*(TEMPTMS/293.d0) ! equation 5 in kerminen
         eta = gamma*CS/GR
c         print*,'fn1',fn
         if (Dp1.gt.d1)then
            fn = fn*exp(eta/(Dp1*1.0D9)-eta/(d1*1.0D9))
         endif
c         print*,'fn2',fn
         if( IT_IS_NAN( fn ) ) then
            print*, '---------------->>> Found NAN in GETNUCRATE'
            print*,'fn',fn
            print*,'eta',eta, 'Dp1',Dp1,'d1',d1
            print*,'gamma',gamma,'CS',CS,'GR',GR,'gtime',gtime
            call ERROR_STOP('Found NaN in fn','getnucrate')
         endif
            
         mnuc = sqrt(xk(1)*xk(2))
      endif
      
      return
      end SUBROUTINE GETNUCRATE
      
!------------------------------------------------------------------------------

      SUBROUTINE VEHK_NUCL (tempi,rhi,cnai,fn,rnuc)
!
!*****************************************************************************
!  Subroutine vehk_nucl calculates the binary nucleation rate and radius of the 
!  critical nucleation cluster using the parameterization of...
!
c     Vehkamaki, H., M. Kulmala, I. Napari, K. E. J. Lehtinen, C. Timmreck, 
C     M. Noppel, and A. Laaksonen. "An Improved Parameterization for Sulfuric 
C     Acid-Water Nucleation Rates for Tropospheric and Stratospheric Conditions." 
C     Journal of Geophysical Research-Atmospheres 107, no. D22 (2002).
!
!  WRITTEN BY Jeff Pierce, April 2007 for GISS GCM-II'
!  Introduce to GEOS-Chem by Win Trivitayanurak Sep 29,2008
!
!*****************************************************************************
!
      ! Arguments
      real*4, intent(in)   :: tempi ! temperature of air [K]
      real*4, intent(in)   :: rhi ! relative humidity of air as a fraction
      real*8, intent(in)   :: cnai ! concentration of gas phase sulfuric acid [molec cm-3]
      real*8, intent(out)  :: fn ! nucleation rate [cm-3 s-1]
      real*8, intent(out)  :: rnuc ! critical cluster radius [nm]

      ! Local variables
      REAL*8  :: fb0(10),fb1(10),fb2(10),fb3(10),fb4(10),fb(10)
      REAL*8  :: gb0(10),gb1(10),gb2(10),gb3(10),gb4(10),gb(10) ! set parameters
      REAL*8  :: temp                 ! temperature of air [K]
      REAL*8  :: rh                   ! relative humidity of air as a fraction
      REAL*8  :: cna                  ! concentration of gas phase sulfuric acid [molec cm-3]
      REAL*8  :: xstar                ! mole fraction sulfuric acid in cluster
      REAL*8  :: ntot                 ! total number of molecules in cluster
      integer :: i                    ! counter

      ! ADJUSTABLE PARAMETERS

c     Nucleation Rate Coefficients
c
      data fb0 /0.14309, 0.117489, -0.215554, -3.58856, 1.14598,
     $          2.15855, 1.6241, 9.71682, -1.05611, -0.148712        /
      data fb1 /2.21956, 0.462532, -0.0810269, 0.049508, -0.600796,
     $       0.0808121, -0.0160106, -0.115048, 0.00903378, 0.00283508/
      data fb2 /-0.0273911, -0.0118059, 0.00143581, -0.00021382, 
     $       0.00864245, -0.000407382, 0.0000377124, 0.000157098,
     $       -0.0000198417, -9.24619d-6 /
      data fb3 /0.0000722811, 0.0000404196, -4.7758d-6, 3.10801d-7,
     $       -0.0000228947, -4.01957d-7, 3.21794d-8, 4.00914d-7,
     $       2.46048d-8, 5.00427d-9 /
      data fb4 /5.91822, 15.7963, -2.91297, -0.0293333, -8.44985,
     $       0.721326, -0.0113255, 0.71186, -0.0579087, -0.0127081  / 

c     Coefficients of total number of molecules in cluster     
c
      data gb0 /-0.00295413, -0.00205064, 0.00322308, 0.0474323,
     $         -0.0125211, -0.038546, -0.0183749, -0.0619974,
     $         0.0121827, 0.000320184 /
      data gb1 /-0.0976834, -0.00758504, 0.000852637, -0.000625104,
     $         0.00580655, -0.000672316, 0.000172072, 0.000906958,
     $         -0.00010665, -0.0000174762 /      
      data gb2 /0.00102485, 0.000192654, -0.0000154757, 2.65066d-6,
     $         -0.000101674, 2.60288d-6, -3.71766d-7, -9.11728d-7,
     $         2.5346d-7, 6.06504d-8 /
      data gb3 /-2.18646d-6, -6.7043d-7, 5.66661d-8, -3.67471d-9,
     $         2.88195d-7, 1.19416d-8, -5.14875d-10, -5.36796d-9,
     $         -3.63519d-10, -1.42177d-11 /
      data gb4 /-0.101717, -0.255774, 0.0338444, -0.000267251,
     $         0.0942243, -0.00851515, 0.00026866, -0.00774234,
     $         0.000610065, 0.000135751 /

      !=================================================================
      ! VEHK_NUCL begins here!
      !=================================================================
      temp=dble(tempi)
      rh=dble(rhi)
      cna=cnai

c     Respect the limits of the parameterization
      if (cna .lt. 1.d4) then ! limit sulf acid conc
         fn = 0.
         rnuc = 1.
c         print*,'cna < 1D4', cna
         goto 10
      endif
      if (cna .gt. 1.0d11) cna=1.0e11 ! limit sulfuric acid conc  
      if (temp .lt. 230.15) temp=230.15 ! limit temp
      if (temp .gt. 305.15) temp=305.15 ! limit temp
      if (rh .lt. 1d-4) rh=1d-4 ! limit rh
      if (rh .gt. 1.) rh=1. ! limit rh
c
c     Mole fraction of sulfuric acid
      xstar=0.740997-0.00266379*temp-0.00349998*log(cna)
     &   +0.0000504022*temp*log(cna)+0.00201048*log(rh)
     &   -0.000183289*temp*log(rh)+0.00157407*(log(rh))**2.
     &   -0.0000179059*temp*(log(rh))**2.
     &   +0.000184403*(log(rh))**3.
     &   -1.50345d-6*temp*(log(rh))**3.
c 
c     Nucleation rate coefficients 
      do i=1, 10
         fb(i) = fb0(i)+fb1(i)*temp+fb2(i)*temp**2.
     &        +fb3(i)*temp**3.+fb4(i)/xstar
      enddo
c
c     Nucleation rate (1/cm3-s)
      fn = exp(fb(1)+fb(2)*log(rh)+fb(3)*(log(rh))**2.
     &    +fb(4)*(log(rh))**3.+fb(5)*log(cna)
     &    +fb(6)*log(rh)*log(cna)+fb(7)*(log(rh))**2.*log(cna)
     &    +fb(8)*(log(cna))**2.+fb(9)*log(rh)*(log(cna))**2.
     &    +fb(10)*(log(cna))**3.)

c      print*,'in vehk_nuc, fn',fn
c      print*,'cna',cna,'rh',rh,'temp',temp
c      print*,'xstar',xstar
c
c   Cap at 10^6 particles/s, limit for parameterization
      if (fn.gt.1.0d6) then
         fn=1.0d6
      endif
c
c     Coefficients of total number of molecules in cluster 
      do i=1, 10
         gb(i) = gb0(i)+gb1(i)*temp+gb2(i)*temp**2.
     &        +gb3(i)*temp**3.+gb4(i)/xstar
      enddo
c     Total number of molecules in cluster
      ntot=exp(gb(1)+gb(2)*log(rh)+gb(3)*(log(rh))**2.
     &    +gb(4)*(log(rh))**3.+gb(5)*log(cna)
     &    +gb(6)*log(rh)*log(cna)+gb(7)*log(rh)**2.*log(cna)
     &    +gb(8)*(log(cna))**2.+gb(9)*log(rh)*(log(cna))**2.
     &    +gb(10)*(log(cna))**3.)

c     cluster radius
      rnuc=exp(-1.6524245+0.42316402*xstar+0.3346648*log(ntot)) ! [nm]

 10   return
      end SUBROUTINE VEHK_NUCL

!------------------------------------------------------------------------------

C     **************************************************
C     *  bl_nucl                                     *
C     **************************************************

C     WRITTEN BY Jeff Pierce, April 2007

C     This subroutine calculates a simple binary nucleation rate of 1 nm
C     particles.

C      j_1nm = A * [H2SO4]

      SUBROUTINE bl_nucl(cnai,fn,rnuc)

      IMPLICIT NONE

C-----INPUTS------------------------------------------------------------

      double precision cnai                 ! concentration of gas phase sulfuric acid [molec cm-3]

C-----OUTPUTS-----------------------------------------------------------

      double precision fn                   ! nucleation rate [cm-3 s-1]
      double precision rnuc                 ! critical cluster radius [nm]

C-----INCLUDE FILES-----------------------------------------------------

C-----ARGUMENT DECLARATIONS---------------------------------------------

C-----VARIABLE DECLARATIONS---------------------------------------------

      double precision cna                  ! concentration of gas phase sulfuric acid [molec cm-3]
      double precision A   ! prefactor... empirical

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(A=2.0D-6)

C-----CODE--------------------------------------------------------------

      cna=cnai

      fn=A*cna
      rnuc=0.5d0 ! particle diameter of 1 nm

 10   return
      end SUBROUTINE bl_nucl

!------------------------------------------------------------------------------

      SUBROUTINE napa_nucl(tempi,rhi,cnai,nh3ppti,fn,rnuc)
!
!*****************************************************************************
!  Subroutine NAPA_NUCL calculates the ternary nucleation rate and radius of the 
!  critical nucleation cluster using the parameterization of...
!
c     Napari, I., M. Noppel, H. Vehkamaki, and M. Kulmala. "Parametrization of 
c     Ternary Nucleation Rates for H2so4-Nh3-H2o Vapors." Journal of Geophysical 
c     Research-Atmospheres 107, no. D19 (2002).
!
C     WRITTEN BY Jeff Pierce, April 2007 for GISS GCM-II'
!     Introduce to GEOS-Chem by Win Trivitayanurak Sep 29, 2008
!
!  NOTES:
!*****************************************************************************
!
      ! Arguments
      real*4, intent(in) :: tempi ! temperature of air [K]
      real*4, intent(in) :: rhi ! relative humidity of air as a fraction
      real*8, intent(in) :: cnai ! concentration of gas phase sulfuric acid [molec cm-3]
      real*8, intent(in) :: nh3ppti ! concentration of gas phase ammonia
      real*8, intent(out):: fn  ! nucleation rate [cm-3 s-1]
      real*8, intent(out):: rnuc ! critical cluster radius [nm]

      ! Local variables
      real*8    ::  aa0(20),a1(20),a2(20),a3(20),fa(20) ! set parameters
      real*8    ::  fnl                  ! natural log of nucleation rate
      real*8    ::  temp                 ! temperature of air [K]
      real*8    ::  rh                   ! relative humidity of air as a fraction
      real*8    ::  cna                  ! concentration of gas phase sulfuric acid [molec cm-3]
      real*8    ::  nh3ppt               ! concentration of gas phase ammonia
      integer   ::  i                 ! counter

      ! Adjustable parameters
      data aa0 /-0.355297, 3.13735, 19.0359, 1.07605, 6.0916,
     $         0.31176, -0.0200738, 0.165536,
     $         6.52645, 3.68024, -0.066514, 0.65874,
     $         0.0599321, -0.732731, 0.728429, 41.3016,
     $         -0.160336, 8.57868, 0.0530167, -2.32736        /

      data a1 /-33.8449, -0.772861, -0.170957, 1.48932, -1.25378,
     $         1.64009, -0.752115, 3.26623, -0.258002, -0.204098,
     $         -7.82382, 0.190542, 5.96475, -0.0184179, 3.64736,
     $         -0.35752, 0.00889881, -0.112358, -1.98815, 0.0234646/
     
      data a2 /0.34536, 0.00561204, 0.000479808, -0.00796052,
     $         0.00939836, -0.00343852, 0.00525813, -0.0489703,
     $         0.00143456, 0.00106259, 0.0122938, -0.00165718,
     $         -0.0362432, 0.000147186, -0.027422, 0.000904383,
     $         -5.39514d-05, 0.000472626, 0.0157827, -0.000076519/
     
      data a3 /-0.000824007, -9.74576d-06, -4.14699d-07, 7.61229d-06,
     $         -1.74927d-05, -1.09753d-05, -8.98038d-06, 0.000146967,
     $         -2.02036d-06, -1.2656d-06, 6.18554d-05, 3.41744d-06,
     $         4.93337d-05, -2.37711d-07, 4.93478d-05, -5.73788d-07, 
     $         8.39522d-08, -6.48365d-07, -2.93564d-05, 8.0459d-08   /


      !=================================================================
      ! NAPA_NUCL begins here!
      !=================================================================
      temp=dble(tempi)
      rh=dble(rhi)
      cna=cnai
      nh3ppt=nh3ppti

c     Napari's parameterization is only valid within limited area
      if ((cna .lt. 1.d4).or.(nh3ppt.lt.0.1)) then ! limit sulf acid and nh3 conc
         fn = 0.
         rnuc = 1
         goto 10
      endif  
      if (cna .gt. 1.0d9) cna=1.0d9 ! limit sulfuric acid conc
      if (nh3ppt .gt. 100.) nh3ppt=100. ! limit temp  
      if (temp .lt. 240.) temp=240. ! limit temp
      if (temp .gt. 300.) temp=300. ! limit temp
      if (rh .lt. 0.05) rh=0.05 ! limit rh 
      if (rh .gt. 0.95) rh=0.95 ! limit rh

      do i=1,20
         fa(i)=aa0(i)+a1(i)*temp+a2(i)*temp**2.+a3(i)*temp**3.
      enddo

      fnl=-84.7551+fa(1)/log(cna)+fa(2)*log(cna)+fa(3)*(log(cna))**2.
     &  +fa(4)*log(nh3ppt)+fa(5)*(log(nh3ppt))**2.+fa(6)*rh
     &  +fa(7)*log(rh)+fa(8)*log(nh3ppt)/log(cna)+fa(9)*log(nh3ppt)
     &  *log(cna)+fa(10)*rh*log(cna)+fa(11)*rh/log(cna)
     &  +fa(12)*rh
     &  *log(nh3ppt)+fa(13)*log(rh)/log(cna)+fa(14)*log(rh)
     &  *log(nh3ppt)+fa(15)*(log(nh3ppt))**2./log(cna)+fa(16)*log(cna)
     &  *(log(nh3ppt))**2.+fa(17)*(log(cna))**2.*log(nh3ppt)
     &  +fa(18)*rh
     &  *(log(nh3ppt))**2.+fa(19)*rh*log(nh3ppt)/log(cna)+fa(20)
     &  *(log(cna))**2.*(log(nh3ppt))**2.
c
c
      fn=exp(fnl)

      ! Try scaling down the rate by 1e-5 to see how the param is 
      ! doing on the false positive nucleation (win, 12/18/08)
      !sensitivity simulation, change scaling factor down to 1e-4
      fn = fn * 1.e-5

c   Cap at 10^6 particles/cm3-s, limit for parameterization
      if (fn.gt.1.0d6) then
        fn=1.0d6
        fnl=log(fn)
      endif

      rnuc=0.141027-0.00122625*fnl-7.82211d-6*fnl**2.
     &     -0.00156727*temp-0.00003076*temp*fnl
     &     +0.0000108375*temp**2.

 10   return
      end subroutine napa_nucl

!------------------------------------------------------------------------------

      SUBROUTINE getCondSink_kerm(Nko,Mko,CS,Dpmean,Dp1,dens1,
     &         BOXVOL, TEMPTMS, PRES)
!
!*****************************************************************************
!  Subroutine GETCONDSINK_KERM calculates the condensation sink (first order 
!  loss rate of condensing gases) from the aerosol size distribution.
!
!  This is the cond sink in kerminen et al 2004 Parameterization for 
!  new particle formation AS&T Eqn 6.
!
!  Written by Jeff Pierce, May 2007 for GISS GCM-II'
!  Introduced to GEOS-Chem by Win Trivitayanurak, Sep 29, 2008
!
!  NOTES:
C     Nk(ibins) - number of particles per size bin in grid cell
C     Nnuc - number of particles per size bin in grid cell
C     Mnuc - mass of given species in nucleation pseudo-bin (kg/grid cell)
C     Mk(ibins, icomp) - mass of a given species per size bin/grid cell
C     spec - number of the species we are finding the condensation sink for
C     CS - condensation sink [s^-1]
C     sinkfrac(ibins) - fraction of condensation sink from a bin
!
!*****************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN)        :: Nko(ibins), Mko(ibins, icomp)
      REAL*8, INTENT(OUT)       :: CS       ! 
      REAL*8, INTENT(OUT)       :: Dpmean   ! the number mean diameter [m]
      REAL*8, INTENT(OUT)       :: Dp1      ! the size of the first size bin [m]
      REAL*8, INTENT(OUT)       :: dens1    ! the density of the first size bin [kg/m3]
      REAL*4, INTENT(IN)        :: BOXVOL, TEMPTMS, PRES 

      ! Local variables
      integer        :: i,j,k,c           ! counters
      REAL*8         :: pi, R    ! pi and gas constant (J/mol K)
      REAL*8         :: mu                  !viscosity of air (kg/m s)
      REAL*8         :: mfp                 !mean free path of air molecule (m)
      REAL*4         :: Di       !diffusivity of gas in air (m2/s)
      REAL*8         :: Neps     !tolerance for number
      REAL*4         :: density  !density [kg m^-3]
      REAL*8         :: mp       !mass per particle [kg]
      REAL*8         :: Dpk(ibins) !diameter of particle [m]
      REAL*8         :: Kn       !Knudson number
      REAL*8         :: beta(ibins) !non-continuum correction factor
      REAL*8         :: Mktot    !total mass in bin [kg]
      REAL*8         :: Dtot,Ntot ! used on getting the number mean diameter

      parameter(pi=3.141592654, R=8.314) !pi and gas constant (J/mol K)
      parameter(Neps=1.0d10)


      !=================================================================
      ! GETCONDSINK_KERM  begins here!
      !=================================================================

C     get some parameters
      mu=2.5277e-7*TEMPTMS**0.75302
!      mfp=2.0*mu / ( pres*sqrt( 8.0 * 0.6589 / (pi*R*TEMPTMS) ) )  !S&P eqn 8.6
      mfp=2.0*mu / ( pres*sqrt( 8.0 * 0.0289 / (pi*R*TEMPTMS) ) )  !S&P eqn 8.6
c      Di=gasdiff(temp,pres,98.0,Sv(srtso4))
c      print*,'Di',Di

C     get size dependent values
      CS = 0.d0
      Ntot = 0.d0
      Dtot = 0.d0
      do k=1,ibins
         if (Nko(k) .gt. Neps) then
            Mktot=0.d0
            do j=1,icomp
               Mktot=Mktot+Mko(k,j)
            enddo
C     kpc  Density should be changed due to more species involed.
            density=aerodens(Mko(k,srtso4),0.d0,
     &           Mko(k,srtnh4),Mko(k,srtnacl),Mko(k,srtecil),
     &           Mko(k,srtecob),Mko(k,srtocil),Mko(k,srtocob),
     &           Mko(k,srtdust),Mko(k,srth2o))
            mp=Mktot/Nko(k)
         else
                                !nothing in this bin - set to "typical value"
            density=1500.
            mp=1.4*xk(k)
         endif
         Dpk(k)=((mp/density)*(6./pi))**(0.333)
         Kn=2.0*mfp/Dpk(k)      !S&P eqn 11.35 (text)
         CS=CS+0.5d0*(Dpk(k)*Nko(k)/(dble(boxvol)*1.0D-6)*(1+Kn))/
     &        (1.d0+0.377d0*Kn+1.33d0*Kn*(1+Kn))
         Ntot = Ntot + Nko(k)
         Dtot = Dtot + Nko(k)*Dpk(k)
         if (k.eq.1)then
            Dp1=Dpk(k)
            dens1 = density
         endif
      enddo      
      
      if (Ntot.gt.1D15)then
         Dpmean = Dtot/Ntot
      else
         Dpmean = 150.d0
      endif
      
      return
      END SUBROUTINE GETCONDSINK_KERM

!------------------------------------------------------------------------------

      SUBROUTINE getGrowthTime (d1,d2,h2so4,temp,boxvol,density,gtime)
!
!*****************************************************************************
C     This subroutine calculates the time it takes for a particle to grow
C     from one size to the next by condensation of sulfuric acid (and
C     associated NH3 and water) onto particles.
!
C     This subroutine assumes that the growth happens entirely in the kinetic
C     regine such that the dDp/dt is not size dependent.  The time for growth 
C     to the first size bin may then be approximated by the time for growth via
C     sulfuric acid (not including nh4 and water) to the size of the first size bin
C     (not including nh4 and water).
!     WRITTEN BY Jeff Pierce, April 2007 for GISS GCM-II'
!     Introduce to GEOS-Chem by Win Trivitayanurak (win, 9/29/08)
!
C     d1: intial diameter [m]
C     d2: final diameter [m]
c     h2so4: h2so4 ammount [kg]
c     temp: temperature [K]
c     boxvol: box volume [cm3]
!
C     gtime: the time it takes the particle to grow to first size bin [s]
!
!*****************************************************************************
!
      ! Reference to F90 modules
      USE ERROR_MOD,      ONLY : ERROR_STOP, IT_IS_NAN

      ! Arguments
      REAL*8, INTENT(IN)     ::  d1,d2    ! initial and final diameters [m]
      REAL*8, INTENT(IN)     ::  h2so4    ! h2so4 amount [kg]
      real*4, INTENT(IN)     ::  temp     ! temperature [K]
      real*4, INTENT(IN)     ::  boxvol  ! box volume [cm3]
      REAL*8, INTENT(IN)     ::  density  ! density of particles in first bin [kg/m3]
      REAL*8, INTENT(OUT)    ::  gtime    ! the time it will take the particle to grow 
                                          ! to first size bin [s]

      ! Local variables
      REAL*8     ::  pi, R, MW
      REAL*8     ::  csulf    ! concentration of sulf acid [kmol/m3]
      REAL*8     ::  mspeed   ! mean speed of molecules [m/s]
      REAL*8     ::  alpha    ! accomidation coef

      parameter(pi=3.141592654d0, R=8.314d0) !pi and gas constant (J/mol K)
      parameter(MW=98.d0) ! density [kg/m3], mol wgt sulf [kg/kmol]
      parameter(alpha=0.65)

      !=================================================================
      ! GETGROWTHTIME begins here!
      !=================================================================
cc      print *,'h2so4',h2so4,'MW',MW,'boxvol',boxvol,dble(boxvol)

      csulf = h2so4/MW/(dble(boxvol)*1d-6) ! SA conc. [kmol/m3]
      mspeed = sqrt(8.d0*R*dble(temp)*1000.d0/(pi*MW))

C     Kinetic regime expression (S&P 11.25) solved for T
      gtime = (d2-d1)/(4.d0*MW/density*mspeed*alpha*csulf)

      if ( IT_IS_NAN(gtime) ) then
Cjrp
      print*,'IN GET GROWTH TIME'
Cjrp
      print*,'d1',d1,'d2',d2
Cjrp
      print*,'h2so4',h2so4
Cjrp
      print*,'boxvol',boxvol
Cjrp
      print*,'csulf',csulf,'mspeed',mspeed
Cjrp
      print*,'density',density,'gtime',gtime
      call ERROR_STOP('Found NaN in fn','getnucrate')

      endif

      RETURN
      END SUBROUTINE GETGROWTHTIME

!------------------------------------------------------------------------------

      SUBROUTINE NUCLEATION (Nki,Mki,Gci,Nkf,Mkf,Gcf,fn,fn1,totsulf,
     &              nuc_bin,dt,ionrate, surf_area, BOXVOL, BOXMASS,
     &              TEMPTMS, PRES, RHTOMAS, pdbg,lev)
!
!******************************************************************************
C     This subroutine calls the Vehkamaki 2002 and Napari 2002 nucleation
C     parameterizations and gets the binary and ternary nucleation rates.
C     The number of particles added to the first size bin is calculated
C     by comparing the growth rate of the new particles to the coagulation
C     sink.
!     WRITTEN BY Jeff Pierce, April 2007 for GISS GCM-II'
!     Introduce to GEOS-Chem by Win Trivitayanurak (win, 9/30/08)
!
C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Nki(ibins) - number of particles per size bin in grid cell
C     Mki(ibins, icomp) - mass of a given species per size bin/grid cell
C     Gci(icomp-1) - amount (kg/grid cell) of all species present in the
C                    gas phase except water
C     dt - total model time step to be taken (s)

C-----OUTPUTS-----------------------------------------------------------

C     Nkf, Mkf, Gcf - same as above, but final values
C     fn, fn1
!
!*****************************************************************************
!
      ! Reference to F90 modules
      USE ERROR_MOD,      ONLY : ERROR_STOP, IT_IS_NAN

      IMPLICIT NONE

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer j,i,k
      double precision Nki(ibins), Mki(ibins, icomp), Gci(icomp-1)
      double precision Nkf(ibins), Mkf(ibins, icomp), Gcf(icomp-1)
      double precision totsulf
      integer nuc_bin
      double precision dt
      double precision fn       ! nucleation rate of clusters cm-3 s-1
      double precision fn1      ! formation rate of particles to first size bin cm-3 s-1
      REAL*4, INTENT(IN)       :: BOXVOL,  BOXMASS, TEMPTMS
      REAL*4, INTENT(IN)       :: PRES,    RHTOMAS
      LOGICAL  PDBG             ! Signal print for debug
      integer lev ! layer of model

      REAL*8                     ionrate
      REAL*8                     surf_area
C-----VARIABLE DECLARATIONS---------------------------------------------

      double precision nh3ppt   ! gas phase ammonia in pptv
      double precision h2so4    ! gas phase h2so4 in molec cc-1
      double precision rnuc     ! critical nucleation radius [nm]
      double precision gtime    ! time to grow to first size bin [s]
      double precision ltc, ltc1, ltc2 ! coagulation loss rates [s-1]
      double precision Mktot    ! total mass in bin
      double precision neps
      double precision meps
      double precision density  ! density of particle [kg/m3]
      double precision pi
      double precision frac     ! fraction of particles growing into first size bin
      double precision d1,d2    ! diameters of particles [m]
      double precision mp       ! mass of particle [kg]
      double precision mold     ! saved mass in first bin
      double precision mnuc     !mass of nucleation
      double precision sinkfrac(ibins) ! fraction of loss to different size bins
      double precision nadd     ! number to add
      double precision CS       ! kerminan condensation sink [m-2]
      double precision Dpmean   ! the number wet mean diameter of the existing aerosol
      double precision Dp1      ! the wet diameter of bin 1
      double precision dens1    ! density in bin 1 [kg m-3]
      double precision GR       ! growth rate [nm hr-1]
      double precision gamma,eta ! used in kerminen 2004 parameterzation
      double precision drymass,wetmass,WR
      double precision fn_c     ! barrierless nucleation rate
      double precision h1,h2,h3,h4,h5,h6
      double precision dum1,dum2,dum3,dum4   ! dummy variables
      double precision rhin,tempin ! rel hum in

      LOGICAL ERRORSWITCH

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter (neps=1E8, meps=1E-8)
      parameter (pi=3.14159)

C-----CODE--------------------------------------------------------------

      errorswitch = .false.

      h2so4 = Gci(srtso4)/boxvol*1000.d0/98.d0*6.022d23
      nh3ppt = Gci(srtnh4)/17.d0/(boxmass/29.d0)*1d12*
     &             PRES/101325.*273./TEMPTMS ! corrected for pressure (because this should be concentration)

      fn = 0.d0
      fn1 = 0.d0
      rnuc = 0.d0
      gtime = 0.d0

C     if requirements for nucleation are met, call nucleation subroutines
C     and get the nucleation rate and critical cluster size
      if (h2so4.gt.1.d4) then
         if (nh3ppt.gt.0.1.and.tern_nuc.eq.1) then
            call napa_nucl(TEMPTMS,RHTOMAS,h2so4,nh3ppt,fn,rnuc) !ternary nuc
            if (ion_nuc.eq.1.and.ionrate.ge.1.d0) then
               call ion_nucl(h2so4,surf_area,TEMPTMS,ionrate,RHTOMAS,
     &              h1,h2,h3,h4,h5,h6)
            else
               h1=0.d0
            endif
            if (h1.gt.fn)then
               fn=h1
               rnuc=h5
            endif
         elseif (bin_nuc.eq.1) then
            call vehk_nucl(TEMPTMS,RHTOMAS,h2so4,fn,rnuc) !binary nuc
            if ((ion_nuc.eq.1).and.(ionrate.ge.1.d0)) then
               call ion_nucl(h2so4,surf_area,TEMPTMS,ionrate,RHTOMAS,
     &              h1,h2,h3,h4,h5,h6)
            else
               h1=0.d0
            endif
            if (h1.gt.fn)then
               fn=h1
               rnuc=h5
            endif
            if (fn.lt.1.0d-6)then
               fn = 0.d0
            endif
         elseif ((ion_nuc.eq.1).and.(ionrate.ge.1.d0)) then
            call ion_nucl(h2so4,surf_area,TEMPTMS,ionrate,RHTOMAS,,
     &              h1,h2,h3,h4,h5,h6)
            fn=h1
            rnuc=h5
         elseif(ion_nuc.eq.2) then
            ! Yu Ion nucleation
c            ! first we need to calculate the available surface area
c            surf_area = 0.D0
c            do k=1, ibins
c               if (Nki(k) .gt. Neps) then
c                  Mktot=0.d0
c                  do j=1,icomp
c                     Mktot=Mktot+Mki(k,j)
c                  enddo
c                  mp=Mktot/Nki(k)
c                  density=aerodens(Mki(k,srtso4),0.d0,
c     &                    Mki(k,srtnh4),0.d0,Mki(k,srth2o))  ! assume bisulfate
c                  ! diameter = ((mass/density)*(6/pi))**(1/3)
c                  d2 = 1.D6*((mp/density)*(6.D0/pi))**(1.D0/3.D0) ! (micrometers)
c                  ! surface area per particle = pi*diameter**2
c                  surf_area = surf_area + 1.D-6*(Nki(k)/boxvol)*
c     &                        pi*(d2**2.D0) ! (um2 cm-2)
c               endif
c            enddo
            rhin=dble(RHTOMAS*100.d0)
            tempin=dble(TEMPTMS)
            
            call YUJIMN(h2so4, rhin, tempin, ionrate, surf_area,
     &               fn, dum1, rnuc, dum2)
         endif
         if((act_nuc.eq.1).and.(lev.le.7))then
            call bl_nucl(h2so4,fn,rnuc)
         endif
         call cf_nucl(TEMPTMS,RHTOMAS,h2so4,nh3ppt,fn_c) ! use barrierless nucleation as a max
         fn = min(fn,fn_c)
c         if (fn.gt.1.0)then
c            print*, 'fn',fn
c            print*, 'Yu Yes!'
c            print*, 'ionrate',ionrate
c            print*, 'surf_area',surf_area
c         endif
      endif


      if (pdbg) then
         if( bin_nuc == 1 ) then
            print *, 'BINARY cluster form rate : fn',fn
         else
            print *, 'TERNARY cluster form rate: fn',fn
         endif
      endif

C     if nucleation occured, see how many particles grow to join the first size
C     section
      if (fn.gt.0.d0) then

         if(pdbg) print*,'Nki',Nki
         if(pdbg) print*,'Mki',Mki

         call getCondSink_kerm(Nki,Mki,CS,Dpmean,Dp1,dens1,
     &            BOXVOL, TEMPTMS,PRES)
         
         if(pdbg) print*,'CS',CS,'Dpmean',Dpmean,'Dp1',Dp1,'dens1',dens1


         d1 = rnuc*2.d0*1D-9
         drymass = 0.d0
         do j=1,icomp-idiag
            drymass = drymass + Mki(1,j)
         enddo
         wetmass = 0.d0
         do j=1,icomp
            wetmass = wetmass + Mki(1,j)
         enddo
         
         ! to prevent division by zero (win, 10/1/08)
         if(drymass == 0.d0) then
            WR = 1.d0
         else
            WR = wetmass/drymass
         endif
 
         if(pdbg) print*,'rnuc',rnuc,'WR',WR
         if(pdbg) print*,'d1',d1,'Gci(srtso4)',Gci(srtso4),
     &        'TEMP',temptms,'boxvol',boxvol
         
         if( IT_IS_NAN( Gci(srtso4) )) then 
            print*,'rnuc',rnuc,'WR',WR
            print*,'d1',d1,'Gci(srtso4)',Gci(srtso4)
            call ERROR_STOP('Found NaN in Gci','nucleation')
         endif
cc         print*,'[nucleation] Gci',Gci
         call getGrowthTime(d1,Dp1,Gci(srtso4)*WR,TEMPTMS,
     &        boxvol,dens1,gtime)
         if (pdbg) print*,'gtime',gtime

         GR = (Dp1-d1)*1D9/gtime*3600.d0 ! growth rate, nm hr-1
         
         gamma = 0.23d0*(d1*1.0d9)**(0.2d0)*(Dp1*1.0d9/3.d0)**0.075d0*
     &        (Dpmean*1.0d9/150.d0)**0.048d0*(dens1*1.0d-3)**
     &        (-0.33d0)*(TEMPTMS/293.d0) ! equation 5 in kerminen
         eta = gamma*CS/GR

         if (Dp1.gt.d1)then
            fn1 = fn*exp(eta/(Dp1*1.0D9)-eta/(d1*1.0D9))
         else
            fn1 = fn
         endif

         if (pdbg) print*,'eta',eta,'Dp1',Dp1,'d1',d1,'fn1',fn1

         mnuc = sqrt(xk(1)*xk(2))
 
         nadd = fn1
         
         nuc_bin = 1
         
         mold = Mki(nuc_bin,srtso4)
         Mkf(nuc_bin,srtso4) = Mki(nuc_bin,srtso4)+nadd*mnuc*
     &        boxvol*dt
         Nkf(nuc_bin) = Nki(nuc_bin)+nadd*boxvol*dt

         Gcf(srtso4) = Gci(srtso4) ! - (Mkf(nuc_bin,srtso4)-mold)
         Gcf(srtnh4) = Gci(srtnh4)

         if (pdbg) then
            print*, 'nadd',nadd
            print *,'Mass add to bin',nuc_bin,'=',
     &        nadd*mnuc*boxvol*dt
            print *,'Number added',nadd*boxvol*dt
            print *,'Gcf(srtso4)',Gcf(srtso4)
            print *,'Gcf(srtnh4)',Gcf(srtnh4)
         endif
         
         do k=1,ibins
            if (k .ne. nuc_bin)then
               Nkf(k) = Nki(k)
               do i=1,icomp
                  Mkf(k,i) = Mki(k,i)
               enddo
            else
               do i=1,icomp
                  if (i.ne.srtso4) then
                     Mkf(k,i) = Mki(k,i)
                  endif
               enddo
            endif
         enddo
         
         do k=1,ibins
            if (Nkf(k).lt.1.d0) then
               Nkf(k) = 0.d0
               do j=1,icomp
                  Mkf(k,j) = 0.d0
               enddo
            endif
         enddo
#if defined ( DEBUG ) 
!      print *, 'mnfix in tomas_mod:2679'
#endif
         call mnfix(Nkf,Mkf, ERRORSWITCH)
         pdbg = errorswitch ! carry the error signal from mnfix to outside
         if (errorswitch) print*,'NUCLEATION: Error after mnfix'

C     there is a chance that Gcf will go less than zero because we are artificially growing
C     particles into the first size bin.  don't let it go less than zero.

      else
         
         do k=1,ibins
            Nkf(k) = Nki(k)
            do i=1,icomp
               Mkf(k,i) = Mki(k,i)
            enddo
         enddo
         
      endif
      
      pdbg = errorswitch        ! carry the error signal from mnfix to outside
      
      
      RETURN
      END SUBROUTINE NUCLEATION

!------------------------------------------------------------------------------

      SUBROUTINE EZCOND (Nki,Mki,mcondi,spec,Nkf,Mkf,surf_area,
     &      BOXVOL, TEMPTMS, PRES, errswitch)
!
!******************************************************************************
C     This subroutine takes a given amount of mass and condenses it
C     across the bins accordingly.  
C     WRITTEN BY Jeff Pierce, May 2007 for GISS GCM-II'
!     Put in GEOS-Chem by Win T. 9/30/08
!
C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Nki(ibins) - number of particles per size bin in grid cell
C     Mki(ibins, icomp) - mass of a given species per size bin/grid cell [kg]
C     mcond - mass of species to condense [kg/grid cell]
C     spec - the number of the species to condense

C-----OUTPUTS-----------------------------------------------------------

C     Nkf, Mkf - same as above, but final values

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------


C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision Nki(ibins), Mki(ibins, icomp)
      double precision Nkf(ibins), Mkf(ibins, icomp)
      double precision mcondi
      REAL*8           surf_area
      integer spec
      REAL*4, INTENT(IN)       :: BOXVOL, TEMPTMS, PRES
      LOGICAL ERRSWITCH   ! signal error to outside

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,k,c           ! counters
      double precision mcond
      double precision pi, R    ! pi and gas constant (J/mol K)
      double precision CS       ! condensation sink [s^-1]
      double precision sinkfrac(ibins+1) ! fraction of CS in size bin
      double precision Nk1(ibins), Mk1(ibins, icomp)
      double precision Nk2(ibins), Mk2(ibins, icomp)
      double precision madd     ! mass to add to each bin [kg]
      double precision maddp(ibins)    ! mass to add per particle [kg]
      double precision mconds ! mass to add per step [kg]
      integer nsteps            ! number of condensation steps necessary
      integer my_floor, my_ceil       ! the floor and ceiling (temporary)
      double precision eps     ! small number
      double precision tdt      !the value 2/3
      double precision mpo,mpw  !dry and "wet" mass of particle
      double precision WR       !wet ratio
      double precision tau(ibins) !driving force for condensation
      double precision totsinkfrac ! total sink fraction not including nuc bin
      double precision CSeps    ! lower limit for condensation sink
      double precision tot_m,tot_s    !total mass, total sulfate mass
      double precision ratio    ! used in mass correction
      double precision fracch(ibins,icomp)
      double precision totch

      double precision tot_i,tot_f,tot_fa ! used for conservation of mass check
      LOGICAL          PDBG,  ERRORSWITCH
      real*8     zeros(ibins)

C     VARIABLE COMMENTS...

C-----EXTERNAL FUNCTIONS------------------------------------------------


C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.141592654, R=8.314) !pi and gas constant (J/mol K)
      parameter(eps=1.d-40)
      parameter(CSeps=1.d-20)

C-----CODE--------------------------------------------------------------

        pdbg = errswitch ! take the signal for print debug from outside
        errswitch = .false. !signal to terminate with error. Initialize with .false.

      tdt=2.d0/3.d0

      mcond=mcondi

      ! initialize variables
      do k=1,ibins
         Nk1(k)=Nki(k)
         do j=1,icomp
            Mk1(k,j)=Mki(k,j)
         enddo
      enddo


#if defined ( DEBUG ) 
!      print *, 'mnfix in tomas_mod:2804'
#endif
      call mnfix(Nk1,Mk1, errorswitch)
      if(errorswitch) then
         print *, 'EZCOND: MNFIX[1] found error --> TERMINATE'
         errswitch=.true.
         return
      endif

      ! get the sink fractions
      call getCondSink(Nk1,Mk1,spec,CS,sinkfrac,surf_area, ! set Nnuc to zero for this calc
     &               BOXVOL,TEMPTMS, PRES)  

      ! make sure that condensation sink isn't too small
      if (CS.lt.CSeps) then     ! just make particles in first bin
         Mkf(1,spec) = Mk1(1,spec) + mcond
         Nkf(1) = Nk1(1) + mcond/sqrt(xk(1)*xk(2))
         do j=1,icomp
            if (icomp.ne.spec) then
               Mkf(1,j) = Mk1(1,j)
            endif
         enddo
         do k=2,ibins
            Nkf(k) = Nk1(k)
            do j=1,icomp
               Mkf(k,j) = Mk1(k,j)
            enddo
         enddo
         return
      endif

      if (pdbg) then
         print*,'CS',CS
         print*,'sinkfrac',sinkfrac
         print*,'mcond',mcond
      endif

      ! determine how much mass to add to each size bin
      ! also determine how many condensation steps we need
      totsinkfrac = 0.d0
      do k=1,ibins
      totsinkfrac = totsinkfrac + sinkfrac(k) ! get sink frac total not including nuc bin
      enddo
      nsteps = 1
      do k=1,ibins
         if (sinkfrac(k).lt.1.0D-20)then
            madd = 0.d0
         else
            madd = mcond*sinkfrac(k)/totsinkfrac
         endif
         mpo=0.0
         do j=1,icomp-idiag
            mpo=mpo + Mk1(k,j)
         enddo
         if(mpo == 0.0 ) then  ! prevent division by zero (win, 10/16/08)
            my_floor = 0
         else
            my_floor = int(madd*0.00001/mpo)
         endif
          my_ceil = my_floor + 1
         nsteps = max(nsteps,my_ceil) ! don't let the mass increase by more than 10%
      enddo

      if(pdbg) print*,'nsteps',nsteps

      ! mass to condense each step
      mconds = mcond/nsteps

      ! do steps of condensation
      do i=1,nsteps
         if (i.ne.1) then
            call getCondSink(Nk1,Mk1,spec,CS,sinkfrac,surf_area,      ! set Nnuc to zero for this calculation
     &               BOXVOL,TEMPTMS, PRES)  
            totsinkfrac = 0.d0
            do k=1,ibins
             totsinkfrac = totsinkfrac + sinkfrac(k) ! get sink frac total not including nuc bin
            enddo
         endif      
         
         tot_m=0.d0
         tot_s=0.d0
         do k=1,ibins
            do j=1,icomp-idiag
               tot_m = tot_m + Mk1(k,j)
               if (j.eq.srtso4) then
                  tot_s = tot_s + Mk1(k,j)
               endif
            enddo
         enddo

         if (pdbg) print *,'tot_s ',tot_s,' tot_m ',tot_m
         
         ! change criteria to bigger amount (win, 9/30/08)
         if (mcond.gt.tot_m*5.0D-2) then
!prior to 9/30/08
!         if (mcond.gt.tot_m*1.0D-3) then
            if (pdbg) print *,'Entering TMCOND '

            do k=1,ibins
               mpo=0.0
               mpw=0.0
                                !WIN'S CODE MODIFICATION 6/19/06
                                !THIS MUST CHANGED WITH THE NEW dmdt_int.f
               do j=1,icomp-idiag
                  mpo = mpo+Mk1(k,j) !accumulate dry mass
               enddo
               do j=1,icomp
                  mpw = mpw+Mk1(k,j) ! have wet mass include amso4
               enddo
               if( mpo > 0.0 ) then    ! prevent division by zero (win, 10/16/08)
                  WR = mpw/mpo  !WR = wet ratio = total mass/dry mass
               else
                  WR = 1.0
               endif
               if (Nk1(k) .gt. 0.d0) then
                  !Change maddp(k) from mass/no. to be just mass (win,10/3/08)
                  ! this is because in tmcond here, the moxd argument takes
                  ! mass to add for each bin array, not mass/no. array.
                  maddp(k) = mconds*sinkfrac(k)/totsinkfrac
                  !Prior to 10/3/08 (win)
                  !maddp(k) = mconds*sinkfrac(k)/totsinkfrac/Nk1(k)
                  mpw=mpw/Nk1(k)

                  if(pdbg) print*,'mpw',mpw,'maddp',maddp(k),'WR',WR
                  !Change the maddp(k) to accordingly -- adding the /Nk1(k) (win, 10/3/08)
                  tau(k)=1.5d0*((mpw+maddp(k)/Nk1(k)*WR)**tdt-mpw**tdt) 
                  ! Prior to 10/3/08 (win)
                  !tau(k)=1.5d0*((mpw+maddp(k)*WR)**tdt-mpw**tdt) !added WR to moxid term (win, 5/15/06)
c     tau(k)=0.d0
c     maddp(k)=0.d0
               else
                                !nothing in this bin - set tau to zero
                  tau(k)=0.d0
                  maddp(k) = 0.d0
               endif
            enddo
c     print*,'tau',tau
#if defined ( DEBUG ) 
!      print *, 'mnfix in tomas_mod:2942'
#endif
            call mnfix(Nk1,Mk1, errorswitch)
            if (errorswitch) then
               print *, 'EZCOND: MNFIX[2] found error --> TERMINATE'
               errswitch=.true.
               return
            endif
                               ! do condensation
            errorswitch = pdbg
!prior to 9/30/08 from Jeff's version
            call tmcond(tau,xk,Mk1,Nk1,Mk2,Nk2,spec,errorswitch,maddp)

! For SO4 condensation, the last argument should be zeroes (win, 9/30/08)
!            zeros(:) = 0.d0
!            call tmcond(tau,xk,Mk1,Nk1,Mk2,Nk2,spec,errorswitch,zeros)


            if( errorswitch) then
               errswitch=.true.
               print *,'EZCOND: error after TMCOND --> TERMINATE'
               return
            endif
            errorswitch = pdbg


c     call tmcond(tau,xk,Mk1,Nk1,Mk2,Nk2,spec)
C     jrp         totch=0.0
C     jrp         do k=1,ibins
C     jrp            do j=1,icomp
C     jrp               fracch(k,j)=(Mk2(k,j)-Mk1(k,j))
C     jrp               totch = totch + (Mk2(k,j)-Mk1(k,j))
C     jrp            enddo
C     jrp         enddo
c     print*,'fracch',fracch,'totch',totch
            

         elseif (mcond.gt.tot_s*1.0D-12) then
            if (pdbg) print *,'Small mcond: distrib w/ sinkfrac '
            if (pdbg) print *, 'maddp(bin) to add to SO4'
            do k=1,ibins
               if (Nk1(k) .gt. 0.d0) then
                  maddp(k) = mconds*sinkfrac(k)/totsinkfrac
               else
                  maddp(k) = 0.d0
               endif
               if(pdbg) print *, maddp(k)
               Mk2(k,srtso4)=Mk1(k,srtso4)+maddp(k)
               do j=1,icomp
                  if (j.ne.srtso4) then
                     Mk2(k,j)=Mk1(k,j)
                  endif
               enddo
               Nk2(k)=Nk1(k)
            enddo
            if(pdbg) errorswitch = .true.
#if defined ( DEBUG ) 
!      print *, 'mnfix in tomas_mod:2999'
#endif
            call mnfix(Nk2,Mk2, errorswitch)
            if(errorswitch) then
               print *, 'EZCOND: MNFIX[3] found error --> TERMINATE'
               errswitch=.true.
               return
            endif
         else ! do nothing
            if (pdbg) print *,'Very small mcond: do nothing!'
            mcond = 0.d0
            do k=1,ibins
               Nk2(k)=Nk1(k)
               do j=1,icomp
                  Mk2(k,j)=Mk1(k,j)
               enddo
            enddo
         endif
         if (i.ne.nsteps)then
            do k=1,ibins
               Nk1(k)=Nk2(k)
               do j=1,icomp
                  Mk1(k,j)=Mk2(k,j)
               enddo
            enddo            
         endif

      enddo

      do k=1,ibins
         Nkf(k)=Nk2(k)
         do j=1,icomp
            Mkf(k,j)=Mk2(k,j)
         enddo
      enddo

      ! check for conservation of mass
      tot_i = 0.d0
      tot_fa = mcond
      tot_f = 0.d0
      do k=1,ibins
         tot_i=tot_i+Mki(k,srtso4)
         tot_f=tot_f+Mkf(k,srtso4)
         tot_fa=tot_fa+Mki(k,srtso4)
      enddo
      
      if(pdbg) then
         print *,'Check conserv of mass after mcond is distrib'
         print *,' Initial total so4 ',tot_i
         print *,' Final total so4   ',tot_f
         print *,'Percent error=',abs((mcond-(tot_f-tot_i))/mcond)*1e2
      endif

      if (mcond.gt.0.d0.and.
     &    abs((mcond-(tot_f-tot_i))/mcond).gt.0.d0) then
         IF(mcond > 1.d-8 .and. tot_i > 5.d-2)  THEN  !Add a check to check error if mcond is significant (win, 10/2/08)

          IF (abs((mcond-(tot_f-tot_i))/mcond).lt.1.d0 .OR. 
     &            spinup(31.0) ) THEN
          !Prior to 10/2/08 (win)   .. original was Jeff's fix  
          !  ! do correction of mass
          !  ratio = (tot_f-tot_i)/mcond 
          !  if(pdbg) print *,'Mk at mass correction '
          !  if(pdbg) print *,'  ratio',ratio
          !  do k=1,ibins
          !     Mkf(k,srtso4)=Mki(k,srtso4)+
!     &              (Mkf(k,srtso4)-Mki(k,srtso4))/ratio
          !     if(pdbg) print *,Mkf(k,srtso4)
          !  enddo
            
            ! Do mass correction (win, 10/2/08)
            ratio = (tot_i+mcond)/tot_f
            if(pdbg) print *,'Mk at mass correction apply ratio= ',ratio
            do k=1,ibins
               Mkf(k,srtso4)=Mkf(k,srtso4) * ratio
               if(pdbg) print *,Mkf(k,srtso4)
            enddo

            if(pdbg) errorswitch=.true.
#if defined ( DEBUG ) 
!      print *, 'mnfix in tomas_mod:3079'
#endif
            call mnfix(Nkf,Mkf, errorswitch)
            if(errorswitch) then
               print *, 'EZCOND: MNFIX[4] found error --> TERMINATE'
               errswitch=.true.
               return
            endif
          else
            print*,'ERROR in ezcond'
            print*,'Condensation error',(mcond-(tot_f-tot_i))/mcond
            print*,'mcond',mcond,'change',tot_f-tot_i
            print*,'tot_i',tot_i,'tot_fa',tot_fa,'tot_f',tot_f
            print*,'Nki',Nki
            print*,'Nkf',Nkf
            print*,'Mki',Mki
            print*,'Mkf',Mkf
            !Prior to 10/2/08 (win)
!            STOP
            ! Send error signal to outside and terminate with more info (win, 10/2/08)
!!as of 10/27/08, try comment out this signal to stop the run (win, 10/27/08)
!! the problem is that maybe or mostly the mass conservation is ruined becuase of the
!! fudging inside mnfix.
!            ERRSWITCH=.TRUE.
!            RETURN
          ENDIF
         ENDIF
      endif

Cjrp      if (abs(tot_f-tot_fa)/tot_i.gt.1.0D-8)then
Cjrp         print*,'No S conservation in ezcond.f'
Cjrp         print*,'initial',tot_fa
Cjrp         print*,'final',tot_f
Cjrp         print*,'mcond',mcond,'change',tot_f-tot_i
Cjrp         print*,'ERROR',(mcond-(tot_f-tot_i))/mcond
Cjrp      endif

      ! check for conservation of mass
      tot_i = 0.d0
      tot_f = 0.d0
      do k=1,ibins
         tot_i=tot_i+Mki(k,srtnh4)
         tot_f=tot_f+Mkf(k,srtnh4)
      enddo
      if (.not. spinup(14.0)) then
         if (abs(tot_f-tot_i)/tot_i.gt.1.0D-8)then
            if ( tot_i > 1.0d-20 ) then
               print*,'No N conservation in ezcond.f'
               print*,'initial',tot_i
               print*,'final  ',tot_f
            endif
         endif
      endif


      return
      end SUBROUTINE EZCOND

!------------------------------------------------------------------------------

      SUBROUTINE AQOXID ( MOXID, KMIN, I, J, L, State_Met, State_Chm )
!
!******************************************************************************
!  Subroutine AQOXID takes an amount of SO4 produced via in-cloud oxidation and
!  condenses it onto an existing aerosol size distribution.  It assumes that 
!  only particles larger than the critical activation diameter activate and 
!  that all of these have grown to roughly the same size.  Therefore, the mass 
!  of SO4 produced by oxidation is partitioned to the various size bins
!  according to the number of particles in that size bin.  Values of tau are 
!  calculated for each size bin accordingly and the TMCOND subroutine is called
!  to update Nk and Mk. (win, 7/23/07)
!  Originally written by Peter Adams for TOMAS in GISS GCM-II', June 2000
!  Modified by Win Trivitayanurak (win@cmu.edu), Oct 3, 2005
!
!  NOTES: 
!  (1 ) Currently do aqueous oxidation pretending that it's dry, so aerosol 
!       water is just zero (win, 7/23/07)
!  29 Jan 2014 - R. Yantosca - Avoid array temporaries in call to NH4BULKTOBIN
!*****************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TRACERID_MOD,       ONLY : IDTNK1, IDTAW1, IDTNH4
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_STRAT

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND60

      ! Arguments
      REAL*8                        :: MOXID
      INTEGER                       :: KMIN, I, J, L
      TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm    ! Chemistry State object

      ! Local variables
      INTEGER,     PARAMETER :: K_MIN = 4
      REAL*8                 :: Nact, Mact, MPO, AQTAU(IBINS)
      REAL*8                 :: Nko(IBINS), Mko(IBINS, ICOMP)
      REAL*8                 :: Nkf(IBINS), Mkf(IBINS, ICOMP)
      REAL*8,      PARAMETER :: TDT = 2.D0 / 3.D0
      REAL*8                 :: M_OXID(IBINS)
      INTEGER                :: K, MPNUM, JC, TRACNUM
      INTEGER                :: NKID
      LOGICAL                :: PDBG

      REAL*8                 :: Nk(IBINS), Nkd(IBINS)
      REAL*8                 :: Mk(IBINS, ICOMP)
      REAL*8                 :: Mkd(IBINS,ICOMP)
      REAL*8                 :: Gc(ICOMP - 1)
      REAL*8                 :: Gcd(ICOMP - 1)
      REAL*4                 :: BOXVOL
      REAL*4                 :: thresh

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! AQOXID begins here
      !=================================================================

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      STT => State_Chm%Tracers

      PDBG = .FALSE.            !For print debugging 
!debug     IF ( I == 46 .AND. J == 59 .AND. L == 9) PDBG = .TRUE.

      BOXVOL  = State_Met%AIRVOL(I,J,L) * 1.e6 !convert from m3 -> cm3 
      ! Update aerosol water from the current RH
      DO K = 1, IBINS
         CALL EZWATEREQM2( I, J, L, K, State_Met, State_Chm )
      ENDDO


#if   defined ( TOMAS12 ) || defined ( TOMAS15 )
      thresh = 4.0
#else
      thresh = 1.0
#endif

      
      ! Swap GEOSCHEM variables into aerosol algorithm variables
      DO K = 1, IBINS
         NKID = IDTNK1 - 1 + K
         NK(K) = STT(I,J,L,NKID)
         DO JC = 1, ICOMP-IDIAG
            MK(K,JC) = STT(I,J,L,NKID+JC*IBINS)
         ENDDO
         MK(K,SRTH2O) = STT(I,J,L,IDTAW1-1+K)
      ENDDO
      !sfarina - initialize Gc to ensure storenm doesn't go NaN on us.
      DO JC=1, ICOMP-1
         Gc(JC) = 0.d0
      ENDDO

      ! Take the bulk NH4 and allocate to size-resolved NH4
!-----------------------------------------------------------------------------
! Prior to 1/29/14:
! Rewrite subroutine call to avoid array temporaries (bmy, 1/29/14)
!      IF ( SRTNH4 > 0 ) 
!     &     CALL NH4BULKTOBIN( MK(1:IBINS,SRTSO4), STT(I,J,L,IDTNH4), 
!    &                        MK(1:IBINS,SRTNH4) )
!-----------------------------------------------------------------------------
      IF ( SRTNH4 > 0 ) THEN
         CALL NH4BULKTOBIN( MK(:,SRTSO4), STT(I,J,L,IDTNH4), 
     &                      MK(:,SRTNH4) )
      ENDIF

      ! Fix any inconsistencies in M/N distribution 
         CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd) 
!debug      IF ( I == 46 .AND. J == 59 .AND. L == 9) PDBG = .TRUE.
#if defined ( DEBUG ) 
!      print *, 'mnfix in tomas_mod:3225'
#endif
      CALL MNFIX( NK, MK, PDBG )
      IF( PDBG ) THEN
         PRINT *,'AQOXID: MNFIX found error at',I,J,L
         CALL ERROR_STOP('Found bad error in MNFIX',
     &                   'Beginning AQOXID after MNFIX' )
      ENDIF
      MPNUM = 5
      IF ( ND60 > 0 ) 
     &      CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL )
         CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd) 

!debug      IF ( I == 46 .AND. J == 59 .AND. L == 9) 
!     &     call debugprint(Nk,Mk,I,J,L,'AQOXID after MNFIX_1')


      ! Calculate which particles activate
 10   CONTINUE ! Continue here if KMIN has to be lowered
      Nact = 0.d0
      Mact = 0.d0
      DO K = KMIN, IBINS
         Nact = Nact + Nk(k)
         DO JC = 1, ICOMP-IDIAG  !accumulate dry mass exclude NH4
            Mact = Mact + Mk(K,JC)
         ENDDO
      ENDDO
      
      ! No particles to condense on, then just exit AQOXID
      IF ( Nact == 0d0 ) GOTO 20

      ! If condensing mass is too large for the alloted portion of NK
      ! then lower KMIN 
      IF ( ( Mact + MOXID )/ Nact > XK(IBINS+1) / thresh ) THEN
      IF ( KMIN > K_MIN ) THEN 
         KMIN = KMIN - 1
         GOTO 10
      ELSE
         ! If there is really not enough number to condense onto when lower
         ! KMIN to the threshold K_MIN (set to 4), then 
         !  IF current time is within first 2 weeks from initialization
         !    (spin-up), then skip and exit 
         !  IF current time is after the first 2 weeks, then terminate
         !    with an error message.
         IF ( .not. SPINUP(14.0) ) THEN
            !WRITE(*,*) 'Location: ',I,J,L
            !WRITE(*,*) 'Kmin/Nact: ',KMIN,NACT
            !WRITE(*,*) 'MOXID/Mact: ',MOXID,Mact
!            DO K = 1, IBINS
!               WRITE(*,*) 'K, N, MSO4, MH2O: ',K,Nk(k),
!     &              MK(K,SRTSO4),MK(K,SRTH2O)
!            ENDDO
            IF ( MOXID > 5D0 .and.
     &           ITS_IN_THE_STRAT( I, J, L, State_Met ) ) THEN
               CALL ERROR_STOP( 'Too few number for condensing mass',
     &                          'AQOXID:1'                           )
            ELSE
               WRITE(*,*) 'AQOXID WARNING: SO4 mass is being discarded'
               GOTO 20
            ENDIF
         ELSE
            IF ( PDBG ) print *,'AQOXID: Discard mass (spin-up)'
            GOTO 20
         ENDIF
      ENDIF
      ENDIF
               
      ! Calculate Tau (driving force) for each size bin
      MOXID = MOXID/ Nact       !Moxid becomes kg SO4 per activated particle 
                                !NOTE: NOT using kg of H2SO4
      DO K = 1, IBINS
      IF ( K < KMIN ) THEN
         !too small to activate - no sulfate for this particle
         AQTAU(K) = 0.d0
         M_OXID(K) = 0.d0
      ELSE
         !activated particle - calculate appropriate tau
         MPO=0.d0
         DO JC = 1, ICOMP-IDIAG
            MPO = MPO + Mk(K,JC) !accumulate dry mass
         ENDDO
         M_OXID(K) = MOXID * NK(K)
         
         IF (Nk(K) > 0.d0) THEN
            ! Calculate Tau
            MPO = MPO / Nk(K)
            AQTAU(K) = 1.5d0 * ( ( ( MPO + MOXID) ** TDT ) - 
     &                           (   MPO          ** TDT )    )

            ! Error checking for negative Tau
            IF ( AQTAU(K) < 0.d0 ) THEN
            IF ( ABS(AQTAU(K)) < 1.d0 ) THEN
               AQTAU(K)=1.d-50  !0.d0  !try change to tiny number instead of 0d0 (win, 5/28/06)
            ELSE
               PRINT *,' ######### aqoxid.f:  NEGATIVE TAU'
               PRINT *,'Error at',i,j,l,'bin',k
               PRINT *,'aqtau(k)',aqtau(k)
               CALL ERROR_STOP( 'Negative Tau','AQOXID:2' )
            ENDIF
            ENDIF

         ELSE
            ! Nothing in this bin - set tau to zero
            AQTAU(K) = 0.d0
         ENDIF                  ! Nk>0
      ENDIF                     ! K<kmin
      ENDDO                     ! Loop ibins
               
      ! Call condensation algorithm

      ! Swap into Nko, Mko
      Mko(:,:) = 0D0
      DO K = 1, IBINS
         Nko(K) = Nk(K)
         DO JC = 1, ICOMP-IDIAG ! Now do aqoxid "dry" (win, 7/23/07)
            Mko(K,JC) = Mk(K,JC)
         ENDDO

      ENDDO
!debug      IF ( I == 46 .AND. J == 59 .AND. L == 9) PDBG = .TRUE.

      CALL TMCOND( AQTAU, XK, Mko, Nko, Mkf, Nkf, SRTSO4, PDBG, M_OXID )
      IF(.not.SPINUP(60.) .and.  PDBG ) THEN
         write(116,*) 'Error at',i,j,l
      ELSE
         PDBG = .false.
      ENDIF

      ! Swap out of Nkf, Mkf
      DO K = 1, IBINS
         Nk(k)=Nkf(k)
         DO JC = 1, ICOMP-IDIAG
            Mk(K,JC) = Mkf(K,JC)
         ENDDO
      ENDDO

 20   CONTINUE ! Continue here if the process is skipped

      ! Save changes to diagnostic
      MPNUM = 4
      IF ( ND60 > 0 ) 
     &      CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL )

      ! Fix any inconsistencies in M/N distribution 
      CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd) 
#if defined ( DEBUG ) 
!      print *, 'mnfix in tomas_mod:3371'
#endif
      CALL MNFIX( NK, MK, PDBG )
      IF( PDBG ) THEN
         PRINT *,'AQOXID: MNFIX found error at',I,J,L
         CALL ERROR_STOP('Found bad error in MNFIX',
     &                   'End of AQOXID after MNFIX' )
      ENDIF
      MPNUM = 5
      IF ( ND60 > 0 ) 
     &      CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL )

      ! Swap Nk and Mk arrays back to STT 
      DO K = 1, IBINS
         TRACNUM = IDTNK1 - 1 + K
         STT(I,J,L,TRACNUM) = Nk(K)
         DO JC = 1, ICOMP-IDIAG
            TRACNUM = IDTNK1 - 1 + K + IBINS*JC
            STT(I,J,L,TRACNUM) = Mk(K,JC)
         ENDDO
         STT(I,J,L,IDTAW1-1+K) = Mk(K,SRTH2O)
      ENDDO
           

      ! Free pointer memory
      NULLIFY( STT )
               
      ! Return to calling subroutine
      END SUBROUTINE AQOXID

!------------------------------------------------------------------------------

      SUBROUTINE SOACOND
     &    ( MSOA, I, J, L, BOXVOL, TEMPTMS, PRES, State_Chm )
!
!******************************************************************************
!  Subroutine SOACOND takes the SOA calculated via 10% yeild assumption and 
!  condense onto existing aerosol size distribution in a similar manner as in
!  aqoxid.  The difference is that SOA condensational driving force is a 
!  function of the amount of soluble mass existing in each bin, unlike aqoxid
!  where driving force depends on activated number (proportional to surface) of
!  each bin. (win, 9/25/07)
!
!  NOTES:
!  ( 1) There has been several versions of SOAcond using different assumptions
!       about what to condense onto, e.g. soluble mass (org+inorg), only 
!       organic (OCIL), or surface area.  Now rewrite for convenient switch 
!       of different absorptive surrogate, e.g. organic mass, surface area. 
!       (win, 3/5/08)
!  ( 2) Now use MABS as absorbing media (win, 3/5/08)
!
!******************************************************************************
!
      ! Reference to F90 modules
      USE CMN_SIZE_MOD
      USE CMN_DIAG_MOD
      USE ERROR_MOD,          ONLY : ERROR_STOP, IT_IS_NAN
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE TRACERID_MOD,       ONLY : IDTNK1, IDTAW1, IDTNH4

      ! Arguments
      REAL*8                        :: MSOA
      INTEGER,INTENT(IN)            :: I, J, L
      REAL*4, INTENT(IN)            :: BOXVOL, TEMPTMS, PRES
      TYPE(ChmState), INTENT(INOUT) :: State_Chm
      REAL*8                        :: surf_area

      ! Local variables
      REAL*8                 :: MPO, OCTAU(IBINS), ntot, mtot
      REAL*8                 :: Nko(IBINS), Mko(IBINS, ICOMP)
      REAL*8                 :: Nkf(IBINS), Mkf(IBINS, ICOMP)
      REAL*8,      PARAMETER :: TDT = 2.D0 / 3.D0
      REAL*8                 :: MEDTOT, MED(IBINS), MABS(IBINS)
      REAL*8                 :: MKTOT(IBINS), DENSITY, PI
      REAL*8                 :: M_NH4
      REAL*8                 :: CS, sinkfrac(IBINS)
      REAL*8                 :: partfrac(IBINS),avgfrac(IBINS)
      INTEGER                :: K, MPNUM, JC, TRACNUM
      LOGICAL                :: PDBG
      PARAMETER   (PI = 3.141592654D0)

!sfarina
      REAL*8                         :: Nk(IBINS), Nkd(IBINS)
      REAL*8                         :: Mk(IBINS, ICOMP)
      REAL*8                         :: Mkd(IBINS,ICOMP)
      REAL*8                         :: Gc(ICOMP - 1)
      REAL*8                         :: Gcd(ICOMP - 1)
      REAL*4                         :: thresh

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER   :: STT(:,:,:,:)

      !=================================================================
      ! SOACOND begins here
      !=================================================================

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers
      
      pdbg = .false.

#if   defined ( TOMAS12 ) || defined ( TOMAS15 )
      thresh = 1.0
#else
      thresh = 1.0
#endif

      ! Swap GEOSCHEM variables into TOMAS variables
      DO K = 1, IBINS
         TRACNUM = IDTNK1 - 1 + K
         NK(K) = STT(I,J,L,TRACNUM)
         DO JC = 1, ICOMP-IDIAG  ! do I need aerosol water here?
            TRACNUM = IDTNK1 - 1 + K + IBINS*JC
            MK(K,JC) = STT(I,J,L,TRACNUM)
            IF( IT_IS_NAN( MK(K,JC) ) ) THEN
               PRINT *,'+++++++ Found NaN in SOACOND ++++++++'
               PRINT *,'Location (I,J,L):',I,J,L,'Bin',K,'comp',JC
            ENDIF
         ENDDO
         MK(K,SRTH2O) = STT(I,J,L,IDTAW1-1+K)
      ENDDO

      ! Take the bulk NH4 and allocate to size-resolved NH4
!-----------------------------------------------------------------------------
! Prior to 1/29/14:
! Rewrite subroutine call to avoid array temporaries (bmy, 1/29/14)
!      IF ( SRTNH4 > 0 ) 
!     &     CALL NH4BULKTOBIN( MK(1:IBINS,SRTSO4), STT(I,J,L,IDTNH4), 
!     &                        MK(1:IBINS,SRTNH4) )
!-----------------------------------------------------------------------------
      IF ( SRTNH4 > 0 ) THEN
         CALL NH4BULKTOBIN( MK(:,SRTSO4), STT(I,J,L,IDTNH4), 
     &                      MK(:,SRTNH4) )
      ENDIF

      CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd) 

      ! Establish an 30-bin array and accculate the total
      ! of the absorbing media.  The choices can be:
      ! organic mass, surface area, organic+inorganic. (win, 3/5/08)

      MEDTOT = 0.d0
      MED = 0.d0
      mtot = 0.d0
      MKTOT(:) = 0.d0
      ! Accumulate the total absorbing media
      DO K = 1, IBINS

         DO JC = 1, ICOMP
            Mktot(k) = Mktot(k) + Mk(k,jc)
         ENDDO
         mtot = mtot + Mktot(k)
         
         if (absall.eq.1)then ! partition to all mass
            MED(K) = Mktot(k)
            MEDTOT = MEDTOT + Mktot(k)
         else
            MED(K) = Mk(k,srtocil) ! partition to just hydrophilic organic
            MEDTOT = MEDTOT + Mk(k,srtocil)
         endif   

      ENDDO
      
      ! Fraction to each bin for mass partitioning
      do k = 1,IBINS
         partfrac(k) = MED(K) / MEDTOT ! MSOA (kg SOA) become (kg SOA per 
                                   ! total absorbing media)
      enddo

      ! Fraction to each bin for surface condensation
      call getCondSink(Nk,Mk,srtocil,CS,sinkfrac,surf_area,
     &               BOXVOL,TEMPTMS, PRES)  

      do k = 1,IBINS
         avgfrac(k)=soaareafrac*sinkfrac(k)+
     &              (1.d0-soaareafrac)*partfrac(k)
      enddo 

      !temporary
      ntot = 0.d0
      do k = 1, ibins
         ntot = ntot + Nk(k)
      enddo

      IF ( ( Mtot + MSOA ) / Ntot > XK(IBINS+1) / thresh ) THEN
         IF ( .not. SPINUP(14.0) ) THEN
            WRITE(*,*) 'Location: ',I,J,L
            WRITE(*,*) 'Mtot_&_Ntot: ',Mtot, Ntot
            IF ( MSOA > 5d0 ) CALL ERROR_STOP('Too few no. for SOAcond',
     &                                        'SOACOND:1')
         ELSE
            WRITE(*,*) 'SOACOND WARNING: SOA mass is being discarded'
            GOTO 30
         ENDIF
      ENDIF

      DO K = 1, IBINS
         MPO = 0.d0
         DO JC = 1, ICOMP-IDIAG
            MPO = MPO + MK(K,JC)  ! Accumulate dry mass
         ENDDO
         MABS(K) = MSOA * avgfrac(K)
         
         IF ( Nk(K) > 0.d0 ) THEN
            MPO = MPO / Nk(K)
            OCTAU(K) = 1.5d0 * ( ( ( MPO + MABS(K)/Nk(K) ) ** TDT ) - 
     &                           (   MPO                   ** TDT )   )

            ! Error checking for negative Tau
            IF ( OCTAU(K) < 0.d0 ) THEN
            IF ( ABS(OCTAU(K)) < 1.d0 ) THEN
               OCTAU(K)=1.d-50  !0.d0  !try change to tiny number instead of 0d0 (win, 5/28/06)
            ELSE
               PRINT *,' ######### Subroutine SOACOND:  NEGATIVE TAU'
               PRINT *,'Error at',i,j,l,'bin',k
               PRINT *,'octau(k)',octau(k)
               CALL ERROR_STOP( 'Negative Tau','SOACOND:2' )
            ENDIF
            ENDIF
           
         ELSE
            OCTAU(K) = 0.d0
         ENDIF
      ENDDO

      ! Call condensation algorithm
      ! Swap into Nko, Mko
      Mko(:,:) = 0.d0
      DO K = 1, IBINS
         Nko(K) = Nk(K)  
         DO JC = 1, ICOMP-IDIAG    ! Now do SOA condensation "dry" 
            Mko(K,JC) = Mk(K,JC)  ! dry mass excl. nh4
         ENDDO
      ENDDO
!debug      if(i==24.and.j==13)       pdbg = .true.
      CALL TMCOND( OCTAU, XK, Mko, Nko, Mkf, Nkf, SRTOCIL, PDBG, MABS )

      IF( PDBG ) THEN
!         print 12, I,J,L
 12      FORMAT( 'Error in SOAcond at ', 3I4 )
         if( .not. SPINUP(60.) )write(116,*) 'Error in SOACOND at',i,j,l
      ELSE
         PDBG = .false.
      ENDIF

      ! Swap out of Nkf, Mkf
      DO K = 1, IBINS
         Nk(k)=Nkf(k)
         DO JC = 1, ICOMP-IDIAG
            Mk(K,JC) = Mkf(K,JC)
         ENDDO
      ENDDO

 30   CONTINUE

      ! Save changes to diagnostic
      MPNUM = 6
      IF ( ND60 > 0 ) 
     &      CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL )


      ! Fix any inconsistencies in M/N distribution

      ! Swap Nk and Mk arrays back to STT array
      DO K = 1, IBINS
         TRACNUM = IDTNK1 - 1 + K
         STT(I,J,L,TRACNUM) = Nk(K)
         DO JC = 1, ICOMP-IDIAG
            TRACNUM = IDTNK1 - 1 + K + IBINS*JC
            STT(I,J,L,TRACNUM) = Mk(K,JC)
         ENDDO
         STT(I,J,L,IDTAW1-1+K) = Mk(K,SRTH2O)
      ENDDO

      ! Free pointer memory
      NULLIFY( STT )

      ! Return to calling subroutine
      END SUBROUTINE SOACOND
      
!------------------------------------------------------------------------------

       SUBROUTINE MULTICOAG( DT, Nk, Mk, BOXVOL, PRES, TEMPTMS,
     &             PDBG )
      
!
!******************************************************************************
!  Subroutine MULTICOAG performs coagulation on the aerosol size distribution
!  defined by Nk and Mk (number and mass).  See "An Efficient Numerical
!  Solution to the Stochastic Collection Equation", S. Tzivion, G. Feingold,
!  and Z. Levin, J Atmos Sci, 44, no 21, 3139-3149, 1987.  Unless otherwise 
!  noted, all equation references refer to this paper.  Some equations are 
!  taken from "Atmospheric Chemistry and Physics: From Air Pollution to Climate 
!  Change" by Seinfeld and Pandis (S&P).  
!
!   This routine uses a "moving sectional" approach in which the
!     aerosol size bins are defined in terms of dry aerosol mass.
!     Addition or loss of water, therefore, does not affect which bin
!     a particle falls into.  As a result, this routine does not
!     change Mk(water), although water masses are needed to compute
!     particle sizes and, therefore, coagulation coefficients.  Aerosol
!     water masses in each size bin will need to be updated later
!     (in another routine) to reflect changes that result from
!     coagulation.
!  Original implementation in GISS GCM-II' by Peter Adams, June 1999
!  Modified to allow for multicomponent aerosols, February 2000
!  Introduced to GEOS-CHEM by Win Trivitayanurak (win@cmu.edu) July 2007
!  
!  Arguments:
!  DT   (REAL*4)   : Time step (s)
!  PDBG (LOGICAL)  : For signalling print debug
!
!  Some key variables
!  kij represents the coagulation coefficient (cm3/s) normalized by the
!      volume of the GCM grid cell (boxvol, cm3) such that its units are (s-1)
!  dNdt and dMdt are the rates of change of Nk and Mk.  xk contains
!     the mass boundaries of the size bins.  xbar is the average mass
!     of a given size bin (it varies with time in this algorithm).  phi
!     and eff are defined in the reference, equations 13a and b.
!
!  NOTES:
!  ( 1) Add call argument in AERODENS for EC, OC, and dust. (win, 9/3/07) 
!******************************************************************************
!
      ! Arguments
      REAL*4,    INTENT(IN)     :: DT  
      REAL*4,    INTENT(IN)     :: PRES  
      REAL*4,    INTENT(IN)     :: TEMPTMS  
      REAL*4,    INTENT(IN)     :: BOXVOL  
      REAL*8,    INTENT(INOUT)  :: Nk(IBINS)
      REAL*8,    INTENT(INOUT)  :: Mk(IBINS, ICOMP)
      LOGICAL,   INTENT(INOUT)  :: PDBG

      ! Local variables
      INTEGER     :: K, J, I, JJ
      REAL*8      :: dNdt(ibins), dMdt(ibins,icomp-idiag)
      REAL*8      :: xbar(ibins), phi(ibins), eff(ibins)
      REAL*4      :: kij(ibins,ibins)
      REAL*4      :: Dpk(ibins)             !diameter (m) of particles in bin k
      REAL*4      :: Dk(ibins)              !Diffusivity (m2/s) of bin k particles
      REAL*4      :: ck(ibins)              !Mean velocity (m/2) of bin k particles
      REAL*4      :: olddiff                !used to iterate to find diffusivity
      REAL*4      :: density                !density (kg/m3) of particles
      REAL*4      :: mu                     !viscosity of air (kg/m s)
      REAL*4      :: mfp                    !mean free path of air molecule (m)
      REAL*4      :: Kn                     !Knudsen number of particle
      REAL*8      :: mp                     !particle mass (kg)
      REAL*4      :: beta                   !correction for coagulation coeff.
!      real*8, external ::   aerodens  !<tmp> try change to double precision (win, 1/4/06)

      !temporary summation variables
      REAL*8       :: k1m(icomp-idiag),k1mx(icomp-idiag)
      REAL*8       :: k1mx2(icomp-idiag)
      REAL*8       :: k1mtot,k1mxtot
      REAL*8       :: sk2mtot, sk2mxtot
      REAL*8       :: sk2m(icomp-idiag), sk2mx(icomp-idiag)
      REAL*8       :: sk2mx2(icomp-idiag)
      REAL*8       :: High_in
      REAL*8       :: mtotal, mktot

      REAL*4      ::zeta                      !see reference, eqn 6
      REAL*4      ::tlimit, dtlimit, itlimit  !fractional change in M/N allowed in one time step
      REAL*4      ::dts  !internal time step (<dt for stability)
      REAL*4      ::tsum !time so far
      REAL*8      :: Neps !minimum value for Nk
cdbg
      character*12 limit        !description of what limits time step

      REAL*8       :: mi, mf   !initial and final masses

#if  defined( TOMAS12 ) || defined( TOMAS15 )
      parameter(zeta=1.28125 , dtlimit=0.25, itlimit=10.)
#else
      parameter(zeta=1.0625, dtlimit=0.25, itlimit=10.)
#endif
      REAL*4      ::pi, kB  !kB is Boltzmann constant (J/K)
      REAL*4      ::R       !gas constant (J/ mol K)
      parameter (pi=3.141592654, kB=1.38e-23, R=8.314, Neps=1.0e-3)

      REAL*8      :: M_NH4  

      LOGICAL     :: ERRSPOT

!sfarina

 1    format(16E15.3)

      !=================================================================
      ! MULTICOAG begins here!
      !=================================================================
      tsum = 0.0

C If any Nk are zero, then set them to a small value to avoid division by zero
      do k=1,ibins
         if (Nk(k) .lt. Neps) then
            Nk(k)=Neps
#if  defined( TOMAS12 ) || defined( TOMAS15 )
            Mk(k,srtso4)=Neps*sqrt( xk(k)*xk(k+1) ) !make the added particles SO4
#else
            Mk(k,srtso4)=Neps*1.4d0*xk(k) !make the added particles SO4
#endif
         endif
      enddo

C Calculate air viscosity and mean free path

      mu=2.5277e-7*temptms**0.75302
      mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*R*temptms)))  !S&P eqn 8.6

      !<temp> 
!      write(6,*)'+++ Nk(1:30)    =',Nk(1:30)
!      write(6,*)'+++ Mk(1:30,SO4)=',Mk(1:30,srtso4)
!      write(6,*)'+++ Mk(1:30,H2O)=',Mk(1:30,srth2o)
      if (pdbg) call debugprint(Nk,Mk,0,0,0,'Inside MULTICOAG')
C Calculate particle sizes and diffusivities
      do k=1,ibins

         IF ( SRTNH4 > 0 ) THEN 
            M_NH4 = Mk(k,SRTNH4) 
         ELSE
            M_NH4 = 0.1875d0*Mk(k,srtso4)  !assume bisulfate
         ENDIF
!tmp         write(6,*)'+++ multicoag:  Mk(',k,'srtso4)=',Mk(k,srtso4)
         density=aerodens(Mk(k,srtso4),0.d0, M_NH4,
     &           Mk(k,srtnacl), Mk(k,srtecil), Mk(k,srtecob),
     &           Mk(k,srtocil), Mk(k,srtocob), Mk(k,srtdust),
     &           Mk(k,srth2o))     !use Mk for sea salt mass(win, 4/18/06)
         !Update mp calculation to include all species (win, 4/18/06)

!prior to 9/26/08 (win)
!         Mktot=0.1875d0*Mk(k,srtso4) !start with NH4 mass

         Mktot = M_NH4         ! start with ammonium (win, 9/26/08)
         Mktot = Mktot + Mk(k,srth2o) ! then include water

         do j=1, icomp-idiag
            Mktot=Mktot+Mk(k,j)
         enddo
         mp=Mktot/Nk(k)
         Dpk(k)=((mp/density)*(6./pi))**(0.333)
         Kn=2.0*mfp/Dpk(k)                            !S&P Table 12.1
         Dk(k)=kB*temptms/(3.0*pi*mu*Dpk(k))             !S&P Table 12.1
     &   *((5.0+4.0*Kn+6.0*Kn**2+18.0*Kn**3)/(5.0-Kn+(8.0+pi)*Kn**2))
         ck(k)=sqrt(8.0*kB*temptms/(pi*mp))              !S&P Table 12.1
      enddo

C Calculate coagulation coefficients

      do i=1,ibins
         do j=1,ibins
            Kn=4.0*(Dk(i)+Dk(j))          
     &        /(sqrt(ck(i)**2+ck(j)**2)*(Dpk(i)+Dpk(j))) !S&P eqn 12.51
            beta=(1.0+Kn)/(1.0+2.0*Kn*(1.0+Kn))          !S&P eqn 12.50
            !This is S&P eqn 12.46 with non-continuum correction, beta
            kij(i,j)=2.0*pi*(Dpk(i)+Dpk(j))*(Dk(i)+Dk(j))*beta
            kij(i,j)=kij(i,j)*1.0e6/boxvol  !normalize by grid cell volume
         enddo
      enddo


 10   continue     !repeat process here if multiple time steps are needed

      if(pdbg) print*,'In the time steps loop +++++++++++++'
      

C Calculate xbar, phi and eff
#if  defined( TOMAS12 ) || defined( TOMAS15 )
      do k=1,ibins

         xbar(k)=0.0
         do j=1,icomp-idiag
            xbar(k)=xbar(k)+Mk(k,j)/Nk(k)            !eqn 8b
         enddo
       if(k.lt.ibins-1)then !from 1 to 10 bins

         eff(k)=2./9.*Nk(k)/xk(k)
     *        *(4.-xbar(k)/xk(k)) !eqn 4 in tzivion 1999
         phi(k)=2./9.*Nk(k)/xk(k)
     *        *(xbar(k)/xk(k)-1.) !eqn 4 in tzivion 1999

         !Constraints in equation 15
         if (xbar(k) .lt. xk(k)) then
            eff(k)=2./3.*Nk(k)/xk(k)
            phi(k)=0.0

         else if (xbar(k) .gt. xk(k+1)) then
            phi(k)=2./3.*Nk(k)/xk(k)
            eff(k)=0.0
         endif
       else                      ! from 11 bins to 12 bins
            eff(k)=2./31./31.*Nk(k)/xk(k)
     *           *(32.-xbar(k)/xk(k)) !eqn 4 in tzivion 1999
            phi(k)=2./31./31.*Nk(k)/xk(k)
     *           *(xbar(k)/xk(k)-1.) !eqn 4 in tzivion 1999

         !Constraints in equation 15
         if (xbar(k) .lt. xk(k)) then
            eff(k)=2./31.*Nk(k)/xk(k)
            phi(k)=0.0

         else if (xbar(k) .gt. xk(k+1)) then
            phi(k)=2./31.*Nk(k)/xk(k)
            eff(k)=0.0
         endif
       endif

      enddo

#else
      do k=1,ibins

         xbar(k)=0.0
         do j=1,icomp-idiag
            xbar(k)=xbar(k)+Mk(k,j)/Nk(k)            !eqn 8b
         enddo

         eff(k)=2.*Nk(k)/xk(k)*(2.-xbar(k)/xk(k))    !eqn 13a
         phi(k)=2.*Nk(k)/xk(k)*(xbar(k)/xk(k)-1.)    !eqn 13b
         
         !Constraints in equation 15
         if (xbar(k) .lt. xk(k)) then
            eff(k)=2.*Nk(k)/xk(k)
            phi(k)=0.0
         else if (xbar(k) .gt. xk(k+1)) then
            phi(k)=2.*Nk(k)/xk(k)
            eff(k)=0.0
         endif
      enddo
#endif

C Necessary initializations
         sk2mtot=0.0
         sk2mxtot=0.0
         do j=1,icomp-idiag
            sk2m(j)=0.0
            sk2mx(j)=0.0
            sk2mx2(j)=0.0
         enddo

C Calculate rates of change for Nk and Mk

      do k=1,ibins

         !Initialize to zero
         do j=1,icomp-idiag
            k1m(j)=0.0
            k1mx(j)=0.0
            k1mx2(j)=0.0
         enddo
         High_in=0.0
         k1mtot=0.0
         k1mxtot=0.0

         !Calculate sums
#if  defined( TOMAS12 ) || defined( TOMAS15 )
         do j=1,icomp-idiag
            if (k .gt. 1.and.k.lt.ibins) then
               do i=1,k-1
                  k1m(j)=k1m(j)+kij(k,i)*Mk(i,j)
                  k1mx(j)=k1mx(j)+kij(k,i)*Mk(i,j)*xbar(i)*zeta
                  k1mx2(j)=k1mx2(j)+kij(k,i)*Mk(i,j)*xbar(i)**2.
     *                 *zeta**3.
               enddo
            elseif(k.eq.ibins)then
                  k1m(j)= sk2m(j)+kij(k,k-1)*Mk(k-1,j)
                  k1mx(j)=sk2mx(j)+kij(k,k-1)*Mk(k-1,j)*xbar(k-1)*4.754
                  k1mx2(j)=sk2mx2(j)+kij(k,k-1)*Mk(k-1,j)*xbar(k-1)**2.
     *                 *107.4365
            endif
            k1mtot=k1mtot+k1m(j)
            k1mxtot=k1mxtot+k1mx(j)
         enddo
#else
         do j=1,icomp-idiag
            if (k .gt. 1) then
               do i=1,k-1
                  k1m(j)=k1m(j)+kij(k,i)*Mk(i,j)
                  k1mx(j)=k1mx(j)+kij(k,i)*Mk(i,j)*xbar(i)
                  k1mx2(j)=k1mx2(j)+kij(k,i)*Mk(i,j)*xbar(i)**2
               enddo
            endif
            k1mtot=k1mtot+k1m(j)
            k1mxtot=k1mxtot+k1mx(j)
         enddo
#endif

         if (k .lt. ibins) then
            do i=k+1,ibins
               High_in=High_in+Nk(i)*kij(k,i)
            enddo
         endif

         !Calculate rates of change
#if  defined( TOMAS12 ) || defined( TOMAS15 )
         if(k.lt.ibins-1)then

            dNdt(k)= -Nk(k)*High_in-kij(k,k)*Nk(k)**2.*1.125
     *       -(phi(k)*k1mtot+(eff(k)-phi(k))/6./xk(k)*k1mxtot)
     *       -kij(k,k)*(phi(k)/3.*xbar(k)*Nk(k)+(eff(k)-phi(k))/18.
     *           /xk(k)*zeta*xbar(k)*xbar(k)*Nk(k))

            if (k .gt. 1) then
cyhl Nk*low_in changes to -0.5*Kij*Nk**2.
               dNdt(k)=dNdt(k)+0.625*kij(k-1,k-1)*Nk(k-1)**2
     *          +(phi(k-1)*sk2mtot+(eff(k-1)-phi(k-1))/6./xk(k-1)
     *          *sk2mxtot)
     *          +kij(k-1,k-1)*(phi(k-1)/3.*xbar(k-1)*Nk(k-1)+(eff(k-1)
     *          -phi(k-1))/18./xk(k-1)*zeta*xbar(k-1)*xbar(k-1)
     *          *Nk(k-1))

            endif

            do j=1,icomp-idiag

               dMdt(k,j)= Nk(k)*k1m(j)-Mk(k,j)*High_in ! !term5,term6
     *            -(phi(k)*xk(k+1)*k1m(j)+
     *            (eff(k)+2.*phi(k))/6.*k1mx(j)
     *            +(eff(k)-phi(k))/6./xk(k)*k1mx2(j)) ! term3
     *            - kij(k,k)*Nk(k)*Mk(k,j)/3. ! I assume 1/2Nk and 2/3Mk for half bin
     *            - kij(k,k)*(phi(k)*xk(k+1)*Mk(k,j)/3.
     *            +(eff(k)+2.*phi(k))/6.*zeta*xbar(k)*Mk(k,j)/3.
     *            +(eff(k)-phi(k))/6./xk(k)*zeta**3.*xbar(k)**2.
     *            *Mk(k,j)/3.)

cyhl  Term9(-kij(k,k)*Nk(k)*Mk(k,j)) is cancled out by term6 (k)  

               if (k .gt. 1) then
                  dMdt(k,j)=dMdt(k,j)
     *               +(phi(k-1)*xk(k)*sk2m(j)+(eff(k-1)
     *               +2.*phi(k-1))/6.*sk2mx(j)
     *               +(eff(k-1)-phi(k-1))/6./xk(k-1)*sk2mx2(j)) !term1
     *               +kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j)/3.
     *               +kij(k-1,k-1)*(phi(k-1)*xk(k)*Mk(k-1,j)/3.
     *               +(eff(k-1)+2.*phi(k-1))/6.*zeta
     *               *xbar(k-1)*Mk(k-1,j)/3.+(eff(k-1)-phi(k-1))/6.
     *               /xk(k-1)*zeta**3.*xbar(k-1)**2.*Mk(k-1,j)/3.)
               endif
            enddo
         else if (k.eq.ibins-1)then

            dNdt(k)=0.625*kij(k-1,k-1)*Nk(k-1)**2
     *         +(phi(k-1)*sk2mtot+(eff(k-1)-phi(k-1))/6./xk(k-1)
     *         *sk2mxtot)
     *         +kij(k-1,k-1)*xbar(k-1)*Nk(k-1)/3.*(phi(k-1)
     *         +(eff(k-1)
     *         -phi(k-1))/6./xk(k-1)*zeta*xbar(k-1))

cyhl updated the following
            dNdt(k)=dNdt(k)-Nk(k)*High_in-kij(k,k)*Nk(k)**2.*1.02
     *         -(phi(k)*k1mtot+(eff(k)-phi(k))/62./xk(k)*k1mxtot)
     *         -kij(k,k)*xbar(k)*Nk(k)*0.484*(phi(k)+(eff(k)
     *         -phi(k))/62./xk(k)*4.754*xbar(k))

cyhl I am not sure how it bring 0.5*kij(k-1,k-1)*Nk(k-1)**2 here. But
cyhl It results in much closer result as 30 bins. Apr.27.08

            do j=1,icomp-idiag
               dMdt(k,j)=
     *            +(phi(k-1)*xk(k)*sk2m(j)+(eff(k-1)
     *            +2.*phi(k-1))/6.*sk2mx(j)
     *            +(eff(k-1)-phi(k-1))/6./xk(k-1)*sk2mx2(j)) !term1
     *            +kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j)/3.
     *            +kij(k-1,k-1)*(phi(k-1)*xk(k)*Mk(k-1,j)/3.
     *            +(eff(k-1)+2.*phi(k-1))/6.*zeta
     *            *xbar(k-1)*Mk(k-1,j)/3.+(eff(k-1)-phi(k-1))/6.
     *            /xk(k-1)*zeta**3.*xbar(k-1)**2.*Mk(k-1,j)/3.)

cyhl updated the following
               dMdt(k,j)= dMdt(k,j)+Nk(k)*k1m(j)-Mk(k,j)*High_in ! !term5,term6
     *            -(phi(k)*xk(k+1)*k1m(j)+(eff(k)/62.+0.484*phi(k))
     *            *k1mx(j)+(eff(k)-phi(k))/62./xk(k)*k1mx2(j)) ! term3
     *            -kij(k,k)*Nk(k)*Mk(k,j)*0.103226 ! I assume 1/2Nk and 2/3Mk for half bin
     *            -kij(k,k)*Mk(k,j)*0.484*(phi(k)*xk(k+1)+(eff(k)/62.
     *            +0.484*phi(k))*4.754*xbar(k)
     *            +(eff(k)-phi(k))/62./xk(k)*107.4365*xbar(k)**2.)
            enddo

         else if (k.eq.ibins)then
            dNdt(k)=-Nk(k)*High_in-kij(k,k)*Nk(k)**2.*1.103226
     *         -(phi(k)*k1mtot+(eff(k)-phi(k))/62./xk(k)*k1mxtot)
     *         -kij(k,k)*0.484*xbar(k)*Nk(k)*(phi(k)+(eff(k)-phi(k))
     *         /62./xk(k)*4.754*xbar(k))
     *         +0.52*kij(k-1,k-1)*Nk(k-1)**2
     *         +(phi(k-1)*sk2mtot+(eff(k-1)-phi(k-1))/62./xk(k-1)
     *         *sk2mxtot)
     *         +kij(k-1,k-1)*xbar(k-1)*Nk(k-1)*0.484*(phi(k-1)
     *         +(eff(k-1)
     *         -phi(k-1))/62./xk(k-1)*4.754*xbar(k-1))

            do j=1,icomp-idiag
               dMdt(k,j)= Nk(k)*k1m(j)-Mk(k,j)*High_in ! !term5,term6
     *            -(phi(k)*xk(k+1)*k1m(j)+(eff(k)/62.+0.484
     *            *phi(k))*k1mx(j)
     *            +(eff(k)-phi(k))/62./xk(k)*k1mx2(j)) ! term3
     *            -kij(k,k)*Nk(k)*Mk(k,j)*0.103226 ! I assume 1/2Nk and 2/3Mk for half bin
     *            -kij(k,k)*Mk(k,j)*0.484*(phi(k)*xk(k+1)+(eff(k)/62.
     *            +0.484*phi(k))*4.754*xbar(k)
     *            +(eff(k)-phi(k))/62./xk(k)*107.4365*xbar(k)**2.)
     *            +(phi(k-1)*xk(k)*sk2m(j)+(eff(k-1)/62.+0.484
     *            *phi(k-1))
     *            *sk2mx(j)+(eff(k-1)-phi(k-1))/62./xk(k-1)*sk2mx2(j)) !term1
     *            +kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j)*0.103226
     *            +kij(k-1,k-1)*Mk(k-1,j)*0.484*(phi(k-1)*xk(k)
     *            +(eff(k-1)/62.+0.484*phi(k-1))*4.754*xbar(k-1)
     *            +(eff(k-1)-phi(k-1))/62./xk(k-1)*107.4365
     *            *xbar(k-1)**2.)
            enddo
         endif

#else
         dNdt(k)= 
     &           -kij(k,k)*Nk(k)**2
     &           -phi(k)*k1mtot
     &           -zeta*(eff(k)-phi(k))/(2*xk(k))*k1mxtot
     &           -Nk(k)*High_in
         if (k .gt. 1) then
         dNdt(k)=dNdt(k)+
     &           0.5*kij(k-1,k-1)*Nk(k-1)**2
     &           +phi(k-1)*sk2mtot
     &           +zeta*(eff(k-1)-phi(k-1))/(2*xk(k-1))*sk2mxtot
         endif

         do j=1,icomp-idiag
            dMdt(k,j)= 
     &           +Nk(k)*k1m(j)
     &           -kij(k,k)*Nk(k)*Mk(k,j)
     &           -Mk(k,j)*High_in
     &           -phi(k)*xk(k+1)*k1m(j)
     &           -0.5*zeta*eff(k)*k1mx(j)
     &           +zeta**3*(phi(k)-eff(k))/(2*xk(k))*k1mx2(j)
            if (k .gt. 1) then
               dMdt(k,j)=dMdt(k,j)+
     &           kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j)
     &           +phi(k-1)*xk(k)*sk2m(j)
     &           +0.5*zeta*eff(k-1)*sk2mx(j)
     &           -zeta**3*(phi(k-1)-eff(k-1))/(2*xk(k-1))*sk2mx2(j)
            endif
cdbg            if (j. eq. srtso4) then
cdbg               if (k. gt. 1) then
cdbg                  write(*,1) Nk(k)*k1m(j), kij(k,k)*Nk(k)*Mk(k,j),
cdbg     &               Mk(k,j)*in, phi(k)*xk(k+1)*k1m(j),
cdbg     &               0.5*zeta*eff(k)*k1mx(j),
cdbg     &               zeta**3*(phi(k)-eff(k))/(2*xk(k))*k1mx2(j),
cdbg     &               kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j),
cdbg     &               phi(k-1)*xk(k)*sk2m(j),
cdbg     &               0.5*zeta*eff(k-1)*sk2mx(j),
cdbg     &               zeta**3*(phi(k-1)-eff(k-1))/(2*xk(k-1))*sk2mx2(j)
cdbg               else
cdbg                  write(*,1) Nk(k)*k1m(j), kij(k,k)*Nk(k)*Mk(k,j),
cdbg     &               Mk(k,j)*in, phi(k)*xk(k+1)*k1m(j),
cdbg     &               0.5*zeta*eff(k)*k1mx(j),
cdbg     &               zeta**3*(phi(k)-eff(k))/(2*xk(k))*k1mx2(j)
cdbg               endif
cdbg            endif
         enddo
#endif

cdbg
         if(pdbg) write(*,*) 'k,dNdt,dMdt: ', k, dNdt(k), dMdt(k,srtso4)

         !Save the summations that are needed for the next size bin
         sk2mtot=k1mtot
         sk2mxtot=k1mxtot
         do j=1,icomp-idiag
            sk2m(j)=k1m(j)
            sk2mx(j)=k1mx(j)
            sk2mx2(j)=k1mx2(j)
         enddo

      enddo  !end of main k loop

C Update Nk and Mk according to rates of change and time step

      !If any Mkj are zero, add a small amount to achieve finite
      !time steps
      do k=1,ibins
         do j=1,icomp-idiag
            if (Mk(k,j) .eq. 0.d0) then
               !add a small amount of mass
               mtotal=0.d0
               do jj=1,icomp-idiag
                  mtotal=mtotal+Mk(k,jj)
               enddo
               Mk(k,j)=1.d-10*mtotal
            endif
         enddo
      enddo

      call mnfix(NK, MK, PDBG)

      !Choose time step
      dts=dt-tsum      !try to take entire remaining time step
      limit='comp'
      do k=1,ibins
         if(pdbg) print*,'At bin ',k
         if (Nk(k) .gt. Neps) then
            !limit rates of change for this bin
            if (dNdt(k) .lt. 0.0) tlimit=dtlimit
            if (dNdt(k) .gt. 0.0) tlimit=itlimit
            if (abs(dNdt(k)*dts) .gt. Nk(k)*tlimit) then 
               dts=Nk(k)*tlimit/abs(dNdt(k))
               if(pdbg) print*,'tlimit',tlimit,'Nk(',k,')',Nk(k), 
     &              'dNdt',dNdt(k), ' == dts ',dts
               limit='number'
               if(pdbg) write(limit(8:9),'(I2)') k
               if(pdbg) write(*,*) Nk(k), dNdt(k)
            endif
            do j=1,icomp-idiag
               !limit rates of change x(win, 4/22/06)
               if (dMdt(k,j) .lt. 0.0) tlimit=dtlimit
               if (dMdt(k,j) .gt. 0.0) tlimit=itlimit
               if (abs(dMdt(k,j)*dts) .gt. Mk(k,j)*tlimit) then 
                  mtotal=0.d0
                  do jj=1,icomp-idiag
                     mtotal=mtotal+Mk(k,jj)
                  enddo
                  !only use this criteria if this species is significant
                  if ((Mk(k,j)/mtotal) .gt. 1.d-5) then
                    dts=Mk(k,j)*tlimit/abs(dMdt(k,j))
                    if(pdbg) print*,'tlimit',tlimit,'Mk(',k,j,')',
     &                Mk(k,j),'dMdt',dMdt(k,j), ' == dts ',dts
                  else
                    if (dMdt(k,j) .lt. 0.0) then
                       !set dmdt to 0 to avoid very small mk going negative
                       dMdt(k,j)=0.0
                       if(pdbg) print*,' dMdt(k,j) < 0 '
                    endif
                  endif
                  limit='mass'
                  if(pdbg) write(limit(6:7),'(I2)') k
                  if(pdbg) write(limit(9:9),'(I1)') j
                  if(pdbg) write(*,*) Mk(k,j), dMdt(k,j)
               endif
            enddo
         else
            !nothing in this bin - don't let it affect time step
            Nk(k)=Neps
#if  defined( TOMAS12 ) || defined( TOMAS15 )
            Mk(k,srtso4)=Neps*sqrt(xk(k)*xk(k+1)) !make the added particles SO4
#else
            Mk(k,srtso4)=Neps*1.4d0*xk(k) !make the added particles SO4
#endif
            !make sure mass/number don't go negative
            if (dNdt(k) .lt. 0.0) dNdt(k)=0.0
            if (pdbg) print*,' dNdt(k) < 0 '
            do j=1,icomp-idiag
               if (dMdt(k,j) .lt. 0.0) dMdt(k,j)=0.0
            enddo
         endif
      enddo  !loop bin

       if (pdbg .and. dts .lt. 1. ) then
         write(*,*), dts, 'dts < 1. in multicoag'
       endif

       if (dts .eq. 0.) then
          write(*,*) 'time step is 0 in multicoag - inf/nan/tiny error'
C          pause
         do k = 1,ibins
            print *, 'dNdt(k)', dNdt(k)
            print *, 'dMdt(k,j)'
            do j = 1,icomp - idiag
               print *, dMdt(k,j)
            end do
         end do
      
          call debugprint(nk, mk, 0,0,0,
     &         'MULTICOAG before terminate: dts=0')
          PDBG = .true.
          return
C          stop
C       go to 20
       endif

      !Change Nk and Mk
cdbg
      if(pdbg) write(*,*) 'tsum=',tsum+dts,' ',limit
      do k=1,ibins
         Nk(k)=Nk(k)+dNdt(k)*dts
         do j=1,icomp-idiag
            Mk(k,j)=Mk(k,j)+dMdt(k,j)*dts
         enddo
      enddo

      !Update time and repeat process if necessary
      tsum=tsum+dts
      if (tsum .lt. dt) then
C         print*,'tsum',tsum, 'less than 3600. loop again'
         goto 10
      endif

      
      RETURN
      END SUBROUTINE MULTICOAG

!------------------------------------------------------------------------------

      SUBROUTINE SO4COND(Nki,Mki,Gci,Nkf,Mkf,Gcf,dt,
     &       BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, errspot)
!
!******************************************************************************
!  Subroutine SO4COND determines the condensational driving force for mass 
!  transfer of sulfate between gas and aerosol phases.  It then calls a mass-
!  and number-conserving algorithm for condensation (/evaporation) of aerosol.
!
!     An adaptive time step is used to prevent numerical difficulties.
!     To account for the changing gas phase concentration of sulfuric
!     acid, its decrease during a condensational time step is well-
!     approximated by an exponential decay with a constant, sK (Hz).
!     sK is calculated from the mass and number distribution of the
!     aerosol.  Not only does this approach accurately take into account
!     the changing sulfuric acid concentration, it is also used to
!     predict (and limit) the final sulfuric acid concentration.  This
!     approach is more accurate and faster (allows longer condensational
!     time steps) than assuming a constant supersaturation of sulfate.
!
!  INPUTS: Initial values of
!     Nki(ibins) - number of particles per size bin in grid cell
!     Mki(ibins, icomp) - mass of a given species per size bin/grid cell
!     Gci(icomp-1) - amount (kg/grid cell) of all species present in the
!                    gas phase except water
!     dt - total model time step to be taken (s)
!  OUTPUTS: 
!     Nkf, Mkf, Gcf - same as above, but final values
!
!  Written by Peter Adams, June 2000, based on thermocond.f
!  Introduced to GEOS-CHEM by Win Trivitayanurak (win@cmu.edu) July 2007
!
!  NOTES: 
!  (1 ) Need clean up later (win, 8/1/07)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,      ONLY : ERROR_STOP

      ! Arguments
      REAL*8      :: Nki(ibins), Mki(ibins, icomp), Gci(icomp-1)
      REAL*8      :: Nkf(ibins), Mkf(ibins, icomp), Gcf(icomp-1)
      REAL*4      :: dt
      LOGICAL     :: errspot  

      REAL*4, INTENT(IN)       :: BOXVOL,  BOXMASS, TEMPTMS
      REAL*4, INTENT(IN)       :: PRES, RHTOMAS

C-----VARIABLE DECLARATIONS---------------------------------------------

      REAL*8      :: dp(ibins, icomp-1)  !Driving force for condensation (Pa)
      REAL*8      :: tau(ibins)          !condensation parameter (see cond.f)
      REAL*8      :: atau(ibins, icomp)  !same as tau, but all species
      REAL*8      :: atauc(ibins, icomp) !same as atau, but for const dp
      REAL*4      :: time                !amount of time (s) that has been simulated
      REAL*4      :: cdt                 !internal, adaptive time step
      REAL*4      :: mu                  !viscosity of air (kg/m s)
      REAL*4      :: mfp                 !mean free path of air molecule (m)
      REAL*4      :: Kn                  !Knudsen number of particle
      REAL*4      :: Dpk(ibins)          !diameter (m) of particles in bin k
      REAL*4      :: density             !density (kg/m3) of particles
      INTEGER     :: j,k,jj,kk        !counters
      REAL*8      :: tj(icomp-1), tk(ibins)  !factors used for calculating tau
      REAL*8      :: sK                  !exponential decay const for H2SO4(g)
      REAL*8      :: pi, R            !constants
      REAL*8      :: zeta13             !from Eqn B30 of Tzivion et al. (1989)
      REAL*4      ::Di                  !diffusivity of gas in air (m2/s)
      REAL*4      ::gmw(icomp-1)        !molecular weight of condensing gas
      REAL*4      ::Sv(icomp-1)         !parameter used for estimating diffusivity
      REAL*4      ::alpha(icomp-1)      !accomodation coefficients
      REAL*4      ::beta                !correction for non-continuum
      REAL*8      :: mp         !particle mass (kg)
      REAL*8      :: Nko(ibins), Mko(ibins, icomp), Gco(icomp-1) !output of cond routine
      REAL*8      :: mi, mf  !initial and final aerosol masses (updates Gc)
      REAL*8      :: tr      ! used to calculate time step limits
      REAL*8      :: mc, ttr
      REAL*8      :: Neps     !value below which Nk is insignificant
      REAL*8      :: cthresh  !determines minimum gas conc. for cond.
cdbg      character*12 limit        !description of what limits time step
      REAL*8      :: tdt      !the value 2/3
      REAL*8      :: Ntotf, Ntoto, dNerr  !used to track number cons.
cdbg      integer numcalls          !number of times cond routine is called
      REAL*8      :: Mktot        ! total mass (win, 4/14/06)
      REAL*8      :: zeros(IBINS)

      LOGICAL      :: negvalue  ! negative check variable
      LOGICAL      :: printdebug ! signal received from aerophys to print values for debug (win, 4/8/06)
      LOGICAL      :: tempvar    ! just a temporary variable (win, 4/12/06)

!      REAL*8, EXTERNAL :: AERODENS
!      REAL, EXTERNAL ::  GASDIFF

      PARAMETER(PI=3.141592654, R=8.314) !pi and gas constant (J/mol K)
      PARAMETER(Neps=1.0d10, zeta13=0.98483, cthresh=1.d-16)

      !=================================================================
      ! SO4COND begins here!
      !=================================================================

      negvalue = .false.
      printdebug  = .false.
      tempvar  = .false.

      ! Set some constants
      ! Note: Could have declare this using DATA statement but don't want to
      !       keep modifying when changing the multi-component mass species
      DO J = 1, ICOMP-1
         IF( J == SRTSO4 ) THEN
            gmw(J)  = 98.
            Sv(J)   = 42.88 
            alpha(J)= 0.65
            !alpha from U. Poschl et al., J. Phys. Chem. A, 102, 10082-10089, 1998
         ELSE IF( J == SRTNACL ) THEN
            gmw(J)  = 0.
            Sv(J)   = 42.88  !use 42.88 for all components following Jeff Pierce's code (win,9/26/08)
            alpha(J)= 0.
         ELSE IF( J == SRTECOB .or. J == SRTECIL  
     &           .or. J == SRTOCOB .or. J == SRTOCIL ) THEN
            gmw(J)  = 12.          ! check these values with Jeff again (win, 8/22/07)
            Sv(J)   = 42.88
            alpha(J)= 0. 
         ELSE IF( J == SRTDUST ) THEN 
            gmw(J)  = 0.
            Sv(J)   = 42.88
            alpha(J)= 0.
         ELSE IF( J == SRTNH4 ) THEN 
            gmw(J)  = 0.
            Sv(J)   = 42.88
            alpha(J)= 0.
         ELSE 
            PRINT *, 'Modify SO4cond for the new species'
            CALL ERROR_STOP('SO4COND','Need values for Gmw, Sv, alpha')
         ENDIF
      ENDDO
            
cdbg      numcalls=0
      printdebug = errspot !taking the signal to printdebug from aerophys (win, 4/8/06)
      errspot = .true. !<step4.4> Flag for showing error location outside this subroutine (win,9/21/05)

      dNerr=0.0
      tdt=2.d0/3.d0

      ! Initialize values of Nkf, Mkf, Gcf, and time
      !--------------------------------------------------
      TIME = 0.0                !subroutine exits when time=dt
      DO J = 1, ICOMP-1
         GCF(J) = GCI(J)
      ENDDO
      DO K = 1, IBINS
         NKF(K) = NKI(K)
         DO J = 1, ICOMP
            MKF(K,J) = MKI(K,J)
         ENDDO
      ENDDO

      !Leave everything the same if nothing to condense
      IF ( GCI(SRTSO4) < CTHRESH * BOXMASS ) GOTO 100

      IF ( PRINTDEBUG ) PRINT*,'COND NOW: H2SO4=',Gci(srtso4)

      ! Repeat from this point if multiple internal time steps are needed
      !------------------------------------------------------------------
 10   CONTINUE

C Call thermodynamics to get dp forcings for volatile species
      do k=1,ibins
         do j = 1, icomp-1
!bug fix         do j=1,icomp
            dp(k,j)=0.0
         enddo
      enddo

C Set dp for nonvolatile species
      do k=1,ibins
         !<step4.5> correct the MW of Gcf(srtso4) to be 98. (win, 10/13/05)
         dp(k,srtso4)=(Gcf(srtso4)/98.)/(boxmass/28.9)*pres
      enddo

C Calculate tj and tk factors needed to calculate tau values
      mu=2.5277e-7*temptms**0.75302
      mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*R*temptms)))  !S&P eqn 8.6
      do j=1,icomp-1
         Di=gasdiff(temptms,pres,gmw(j),Sv(j))
         tj(j)=2.*pi*Di*molwt(j)*1.0d-3/(R*temptms)
      enddo
      sK=0.0d0
      do k=1,ibins
         if (Nkf(k) .gt. Neps) then
!            density=aerodens(Mkf(k,srtso4),0.d0,
!     &            0.1875d0*Mkf(k,srtso4),Mkf(k,srtnacl), 
!     &            Mk(k,srtecil), Mk(k,srtecob), 
!     &            Mk(k,srtocil), Mk(k,srtocob),
!     &            Mk(k,srtdust),
!     &            Mkf(k,srth2o)) 
            density=aerodens(Mkf(k,srtso4),0.d0,
     &            0.1875d0*Mkf(k,srtso4),Mkf(k,srtnacl), 
     &            Mkf(k,srtecil), Mkf(k,srtecob), 
     &            Mkf(k,srtocil), Mkf(k,srtocob),
     &            Mkf(k,srtdust),
     &            Mkf(k,srth2o)) 
             !factor of 1.2 assumes ammonium bisulfate
             !(NH4)H has MW of 19 which is = 0.2*96
             !So the Mass of ammonium bisulfate = 1.2*mass sulfate
!win, 4/14/06             mp=(1.2*Mkf(k,srtso4)+Mkf(k,srth2o))/Nkf(k)
            !Need to include new mass species in mp (win, 4/14/06)
            !Add 0.1875x first for ammonium, and then add 1.0x in the loop (win, 4/14/06)
            Mktot=0.1875d0*Mkf(k,srtso4)
            do j=1, icomp
                  Mktot=Mktot+Mkf(k,j)
            enddo
            mp=Mktot/Nkf(k)
            
         else
            !nothing in this bin - set to "typical value"
            density=1500.
#if  defined( TOMAS12 ) || defined( TOMAS15 )
            mp=sqrt(xk(k)*xk(k+1))
#else
            mp=1.4*xk(k)
#endif
         endif
         Dpk(k)=((mp/density)*(6./pi))**(0.333)
         Kn=2.0*mfp/Dpk(k)                             !S&P eqn 11.35 (text)
         beta=(1.+Kn)/(1.+2.*Kn*(1.+Kn)/alpha(srtso4)) !S&P eqn 11.35
         tk(k)=(6./(pi*density))**(1./3.)*beta
         if (Nkf(k) .gt. 0.0) then
            do kk=1,icomp
               sK=sK+tk(k)*Nkf(k)*
!original     &         ((Mkf(k,srtso4)+Mkf(k,srth2o))
!original     &         /Nkf(k))**(1.d0/3.d0)
     &         (Mkf(k,kk)/Nkf(k))**(1.d0/3.d0)  !<step5.1> (win, 4/14/06)
            enddo
         endif
      enddo  !bin loop
      sK=sK*zeta13*tj(srtso4)*R*temptms/(molwt(srtso4)*1.d-3)
     &      /(boxvol*1.d-6)

C Choose appropriate time step

      !Try to take full time step
      cdt=dt-time
cdbg      limit='complete'

      !Not more than 15 minutes
      if (cdt .gt. 900.) then
         cdt=900.
cdbg         limit='max'
      endif

 20   continue   !if time step is shortened, repeat from here

      !Calculate tau values for all species/bins
      do k=1,ibins
         do j=1,icomp
            atauc(k,j)=0.d0
            atau(k,j)=0.d0
         enddo
         !debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if(printdebug)then
         write(*,*)'+++ k loop at',k !<temp>
         write(*,*)'+++ tj(srtso4)', tj(srtso4) !<temp>
         write(*,*)'+++ dp(k,srtso4)', dp(k,srtso4) !<temp>
         write(*,*)'+++ tk(k)',tk(k) !<temp>
         write(*,*)'+++ cdt',cdt !<temp>
         endif
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         atauc(k,srtso4)=tj(srtso4)*tk(k)*dp(k,srtso4)*cdt

         if (sK .gt. 0.d0) then
            atau(k,srtso4)=tj(srtso4)*R*temptms/(molwt(srtso4)*1.d-3)
     &                   /(boxvol*1.d-6)*tk(k)*Gcf(srtso4)/sK
     &                   *(1.d0-exp(-1.d0*sK*cdt))
         else
            !nothing to condense onto
            atau(k,srtso4)=0.d0
         endif
      enddo

      !debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (printdebug) then
         do j=1,icomp
         print *,'atauc(1:ibins,comp) at comp',j
         print *,atauc(1:ibins,j)
         print *,'atau(1:ibins,comp) at comp',j
         print *,atau(1:ibins,j)
         enddo
      endif
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      !The following sections limit the condensation time step
      !when necessary.  tr is a factor that describes by
      !how much to reduce the time step.
!      tr=1.0
      tr=1.d0  !make sure tr is double precision (win, 3/20/05)

      !Make sure masses of individual species don't change too much
      do j=1,icomp-1
         do k=1,ibins
            if (Nkf(k) .gt. Neps) then
               mc=0.d0
               do jj=1,icomp
                  mc=mc+Mkf(k,jj)/Nkf(k)
               enddo
               if (mc/xk(k) .gt. 1.0d-3) then
                  !species has significant mass in particle - limit change
                  if (abs(atau(k,j))/mc**(2.d0/3.d0) .gt. 0.1) then
                     ttr=abs(atau(k,j))/mc**(2.d0/3.d0)/5.d-2
                     if (ttr. gt. tr) then 
                        tr=ttr
cdbg                        limit='amass'
cdbg                        write(limit(7:11),'(I2,X,I2)') k,j
                     endif
                  endif
               else
                  !species is new to particle - set max time step
                if ((cdt/tr .gt. 1.d-1) .and. (atau(k,j).gt. 0.d0)) then
!                     tr=cdt/0.1 
                     tr=cdt/1.d-1 !Make sure tr is double precision (win,3/20/05)
cdbg                     limit='nspec'
cdbg                     write(limit(7:11),'(I2,X,I2)') k,j
                endif
               endif
            endif
         enddo
         !Make sure gas phase concentrations don't change too much
         if (exp(-1.d0*sK*cdt) .lt. 2.5d-1) then
            ttr=-2.d0*cdt*sK/log(2.5d-1)
            if (ttr .gt. tr) then 
               tr=ttr
cdbg               limit='gphas'
cdbg               write(limit(7:8),'(I2)') j
            endif
         endif
      enddo

      !Never shorten timestep by less than half
!      if (tr .gt. 1.0) tr=max(tr,2.0)
      if (tr .gt. 1.d0) tr=max(tr,2.d0) !make sure tr is double precision (win,3/20/05)

      !Repeat for shorter time step if necessary
!      if (tr .gt. 1.0) then
      if (tr .gt. 1.d0) then  !make sure tr is double precision (win,3/20/05)
         cdt=cdt/tr
         goto 20
      endif

C Call condensation subroutine to do mass transfer

      do j=1,icomp-1  !Loop over all aerosol components

         !debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if(printdebug) print *,'Call condensation at comp',j
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         !Swap tau values for this species into array for cond
         do k=1,ibins
            tau(k)=atau(k,j)
         enddo

cdbg         write(*,*) 'so4cond - time = ', time, ' ',limit
cdbg         if (j .eq. srtso4) then
cdbg            do k=1,ibins
cdbg               write(*,'(I3,4E12.4)') 
cdbg     &           k,sK,cdt,atauc(k,srtso4),atau(k,srtso4)
cdbg            enddo
cdbg         endif

C original         call mnfix(Nkf,Mkf)
         if (printdebug) negvalue=.true. !signal received to printdebug (win, 4/8/06)
         call mnfix(Nkf,Mkf,negvalue) !<step5.1> bug fix call argument (win, 4/15/06) !<step4.2> Add call argument to carry tell where mnfix found
                                ! the negative value (win, 9/12/05)
         !-------------------------------------------------------------------
         ! Prior to 1/25/10:
         ! Bug fix: You can't use == to compare logicals, you have to either
         ! use IF ( negvalue .eqv. .TRUE. ) or IF ( negvalue ) (bmy, 1/25/10)
         !if(negvalue==.true.)STOP 'MNFIX terminate' !(win, 9/12/05)
         !-------------------------------------------------------------------
         if ( negvalue ) STOP 'MNFIX terminate' !(win, 9/12/05)

         !Call condensation routine
         Ntotf=0.d0  !Force double precision (win, 4/20/06)
         do k=1,ibins
            Ntotf=Ntotf+Nkf(k)
         enddo

         !<step5.1> Skip tmcond call if there is absolutely no particle (win, 4/20/06)
         if(Ntotf.gt.0d0) then

            !debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if(printdebug) print *,'=== Entering TMCOND ==='
            tempvar = printdebug
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            zeros(:) = 0.d0
            call tmcond(tau,xk,Mkf,Nkf,Mko,Nko,j,printdebug,zeros)
cdbg         numcalls=numcalls+1
            errspot = printdebug !receive the error signal from inside tmcond (win,4/12/06)
            printdebug = tempvar !printdebug gets the originally assigned value (win, 4/12/06)
!8/2/07            if(errspot) goto 100 !Exit so4cond right away when found error from tmcond. (win, 4/13/06)

            !Check for number conservation
            Ntoto=0.0
            do k=1,ibins
               Ntoto=Ntoto+Nko(k)
            enddo
cdbg         write(*,*) 'Time=', time
cdbg         write(*,*) 'Ntoto=', Ntoto
cdbg         write(*,*) 'Ntotf=', Ntotf
            dNerr=dNerr+Ntotf-Ntoto
            if (abs(dNerr/Ntoto) .gt. 1.e-4) then
               write(*,*) 'ERROR in so4cond: Number not conserved'
               write(*,*) 'time=',time
               write(*,*) Ntoto, Ntotf
               write(*,*) (Nkf(k),k=1,ibins)
               errspot = .true. !<step4.4> This flag will trigger printing of location with error (win, 9/21/05)
            endif

         else !(win, 4/20/06)
            if(printdebug) print *,'so4cond: Nk=0 -> skip tmcond'
            do k=1,ibins
               nko(k) = 0d0
               do jj=1,icomp-1
                  Mko(k,jj) = 0d0
               enddo
            enddo
         endif !(win, 4/20/06)

         if(printdebug) print *,'Initial gas conc:',Gcf(j)  !<temp> (win, 4/11/06)

         !Update gas phase concentration
         mi=0.0
         mf=0.0
         do k=1,ibins
            mi=mi+Mkf(k,j)
            mf=mf+Mko(k,j)
         enddo
         Gcf(j)=Gcf(j)+(mi-mf)*gmw(j)/molwt(j)

         if(printdebug) print *,'Updated gas conc:',Gcf(j) !<temp> (win, 4/11/06)

         !Swap into Nkf, Mkf
         do k=1,ibins
            Nkf(k)=Nko(k)
            do jj=1,icomp-1
               Mkf(k,jj)=Mko(k,jj)
            enddo
         enddo

         !Update water concentrations
         call ezwatereqm(Mkf, RHTOMAS)

      enddo

C Update time
      time=time+cdt
!dbg      write(*,*) 'so4cond - time = ', time, ' ',limit
cdbg      write(*,*) 'H2SO4(g)= ', Gcf(srtso4)
      if (Gcf(srtso4) .lt. 0.0) then
         if (abs(Gcf(srtso4)) .gt. 1.d-5) then
            !Gcf is substantially less than zero - this is a problem
            write(*,*) 'ERROR in so4cond: H2SO4(g) < 0'
            write(*,*) 'time=',time
            write(*,*) 'Gcf()=',Gcf(srtso4)
!4/11/06            STOP 
!Let the run STOP outside so4cond so I can know where the run died (win, 4/11/06)
            errspot=.true. !win, 4/11/06
         else
            !Gcf is negligibly less than zero - probably roundoff error
            Gcf(srtso4)=0.0
         endif
      endif

C Repeat process if necessary
      if (time .lt. dt) goto 10

cdbg      write(*,*) 'Cond routine called ',numcalls,' times'
cdbg      write(*,*) 'Number cons. error was ', dNerr

 100  continue   !skip to here if there is no gas phase to condense

      RETURN
      END SUBROUTINE SO4COND

!------------------------------------------------------------------------------

      SUBROUTINE TMCOND(TAU,X,AMKD,ANKD,AMK,ANK,CSPECIES,pdbug,moxd)
!
!******************************************************************************
!  Subroutine TMCOND do condensation calculation.
!  Original code from Peter Adams
!  Modified for GEOS-CHEM by Win Trivitayaurak (win@cmu.edu)
!
!
!  Original Note:
c CONDENSATION
c Based on Tzivion, Feingold, Levin, JAS 1989 and 
c Stevens, Feingold, Cotton, JAS 1996
c--------------------------------------------------------------------

c TAU(k) ......... Forcing for diffusion = (2/3)*CPT*ETA_BAR*DELTA_T
c X(K) ........ Array of bin limits in mass space
c AMKD(K,J) ... Input array of mass moments
c ANKD(K) ..... Input array of number moments
c AMK(K,J) .... Output array of mass moments
c ANK(K) ...... Output array of number moments
c CSPECIES .... Index of chemical species that is condensing
c
c The supersaturation is calculated outside of the routine and assumed
c to be constant at its average value over the timestep.
c 
c The method has three basic components:
c (1) first a top hat representation of the distribution is construced
c     in each bin and these are translated according to the analytic
c     solutions
c (2) The translated tophats are then remapped to bins.  Here if a 
c     top hat entirely or in part lies below the lowest bin it is 
c     not counted.
c 

Cpja Additional notes (Peter Adams)

C     I have changed the routine to handle multicomponent aerosols.  The
C     arrays of mass moments are now two dimensional (size and species).
C     Only a single component (CSPECIES) is allowed to condense during
C     a given call to this routine.  Multicomponent condensation/evaporation
C     is accomplished via multiple calls.  Variables YLC and YUC are
C     similar to YL and YU except that they refer to the mass of the 
C     condensing species, rather than total aerosol mass.

C     I have removed ventilation variables (VSW/VNTF) from the subroutine
C     call.  They still exist internally within this subroutine, but
C     are initialized such that they do nothing.

C     I have created a new variable, AMKDRY, which is the total mass in
C     a size bin (sum of all chemical components excluding water).  I
C     have also created WR, which is the ratio of total wet mass to 
C     total dry mass in a size bin.

C     AMKC(k,j) is the total amount of mass after condensation of species
C     j in particles that BEGAN in bin k.  It is used as a diagnostic
C     for tracking down numerical errors.

Cpja End of my additional notes

!<step5.1> Add a call argument pdbug to signal print values for debugging (win, 4/10/06)
!<step5.3> Add PROC argument to signal gas phase condensation or aq.oxidation (win, 7/13/06)
!<step5.3> Remove PROC (win, 7/17/06)
!<step6.1> Update code for SOA condensation allowing mass conservation fix 
!           similar to when aqoxid calls tmcond. Now madd(ibins) is not just for 
!           SO4 mass but works for SOA mass too. (win, 9/27/07)
!         -Change the argument moxd from one element to be an array moxd(IBINS)
!          and update the way the correction to conserve mass is done (win, 3/5/08)
!******************************************************************************

      ! Local variables
      INTEGER      ::  L,I,J,K,IMN,CSPECIES
      REAL*8       :: DN,DM,DYI,TAU(ibins),XL,XU,YL,YLC,YU,YUC
      REAL*8       :: TEPS,NEPS,EX2,ZERO
      REAL*8       :: XI,XX,XP,YM,WTH,W1,W2,WW,AVG
      REAL*8       :: VSW,VNTF(ibins)
      REAL*8       :: TAU_L, maxtau
      REAL*8       :: X(ibins+1),AMKD(ibins,icomp),ANKD(ibins)
      REAL*8       :: AMK(ibins,icomp),ANK(ibins)
      REAL*8       :: AMKDRY(ibins), WR(ibins), AMKWET(ibins)
      REAL*8       :: AMKDRYSOL(ibins)
      PARAMETER (TEPS=1.0d-40,NEPS=1.0d-20)
      PARAMETER (EX2=2.d0/3.d0,ZERO=0.0d0)
      PARAMETER       NEPS2=1.0d-10
      LOGICAL      :: pdbug !(win, 4/10/06)
      LOGICAL      :: errspot !(win, 4/12/06)
!x      REAL*8       :: moxd !moxid/Nact (win, 5/25/06)
      REAL*8       :: moxd(IBINS) ! condensing mass distributed to size bins
                                  ! according to the selected absorbing media (win, 3/5/08)
      REAL*8       :: c1, c2 !correction factor (win, 5/25/06)
      REAL*8       :: madd(ibins) !condensing mass to be added by aqoxid 
                                  !or SOAcond. For error fixing (win, 9/27/07)
      REAL*8       :: xadd(ibins) !mass per particle to be added by aqoxid
                                  ! or SOAcond. For error fixing (win, 9/27/07)
      REAL*8       :: macc !accumulating the condensing mass (win, 7/24/06)
      REAL*8       :: delt1,delt2 !the delta = mass not conserved (win, 7/24/06)
      REAL*8       :: dummy, xtra,maddtot ! for mass conserv fixing (win, 9/27/07)
      integer      :: kk !counter (wint, 7/24/06)
      REAL*8       :: AMKD_tot
 
      !=================================================================
      ! TMCOND begins here!
      !=================================================================
 
 3    format(I4,200E20.11)

!<step4.5> This first check cause the error of 'number not conserved'
! though only with the small amounts because when ANKD(k) = 0.d0 from start,
! the original check just give it a value NEPS = 1.d-20, and then undergo 
! tmcond calculation.   I'm changing the check to if ANKD(k)= 0.d0,  
! then keep it that way and make the following calculations skip when 
! ANKD(k) is zero (win, 10/18/05)

C If any ANKD are zero, set them to a small value to avoid division by zero
!      do k=1,ibins
!         if (ANKD(k) .lt. NEPS) then
!            ANKD(k)=NEPS
!            AMKD(k,srtso4)=NEPS*1.4*xk(k) !make the added particles SO4
!         endif
!      enddo

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !<step5.1> Add print for debugging (win, 4/10/06)
      if (pdbug) then
         call debugprint(ANKD, AMKD, 0,0,0,'Entering TMCOND')
!         print *, 'TMCOND:entering*************************'
!         print *,'Nk(1:30)'
!         print *, ANKD(1:30)
!         print *,'Mk(1:30,comp)'
!         do j=1,icomp
!         print *,'comp',j
!         print *, AMKD(1:30,j)
!         enddo
      endif
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      errspot = .false. !initialize error signal as false (Win, 4/12/06)

Cpja Sometimes, after repeated condensation calls, the average bin mass
Cpja can be just above the bin boundary - in that case, transfer a some
Cpja to the next highest bin
      do k=1,ibins-1
         if ( ANKD(k) .lt. NEPS2) goto 300 !<step4.5> (win, 10/18/05)
         ! Modify the check to include all dry mass (win, 10/3/08)
         AMKD_tot = 0.d0
         do kk=1,icomp-idiag
            AMKD_tot = AMKD_tot + AMKD(k,kk)
         enddo
         if ((AMKD_tot)/ANKD(k).gt.xk(k+1)) then
         !Prior to 10/3/08 (win)
         !if ((AMKD(k,srtso4))/ANKD(k).gt.xk(k+1)) then
            do j=1,icomp-idiag
               AMKD(k+1,j)=AMKD(k+1,j)+0.1d0*AMKD(k,j)
               AMKD(k,j)=AMKD(k,j)*0.9d0
            enddo
            ANKD(k+1)=ANKD(k+1)+0.1d0*ANKD(k)
            ANKD(k)=ANKD(k)*0.9d0
           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           !<step5.1> Add print for debugging (win, 4/10/06)
            if (pdbug) then
            print *, 'Modified at checkpoint1: BIN',k
            print *,'ANKD(k)',ANKD(k),'ANKD(k+1)',ANKD(k+1)
            print *,'Mk(k,comp)       Mk(k+1,comp)'
            do j=1,icomp
               print *,'comp',j
               print *, AMKD(k,j), AMKD(k+1,j)
            enddo
            endif
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         endif
 300     continue   !<step4.5> If aerosol number is zero (win, 10/18/05)

      enddo

Cpja Initialize ventilation variables so they don't do anything
      VSW=0.0d0
      DO L=1,ibins
         VNTF(L)=0.0d0
      ENDDO

Cpja Initialize AMKDRY and WR
      DO L=1,ibins
         AMKDRY(L)=0.d0
         AMKWET(L)=0.d0
         AMKDRYSOL(L) = 0.d0
         DO J=1,icomp-idiag     ! dry mass excl. nh4 (win, 9/26/08)
            AMKDRY(L)=AMKDRY(L)+AMKD(L,J)
            ! Accumulate the absorbing media (win, 3/5/08)
            IF ( J == SRTOCIL  ) 
     &           AMKDRYSOL(L) = AMKDRYSOL(L) + AMKD(L,J)
        ENDDO
        DO J=1,ICOMP
           AMKWET(L) = AMKWET(L) + AMKD(L,J)
        ENDDO
         if (AMKDRY(L) .gt. 0.d0)    !<step4.5> In case there is no mass, then just skip (win, 10/18/05)
     &        WR(L)= AMKWET(L) / AMKDRY(L)
!original        WR(L)=(AMKDRY(L)+AMKD(L,icomp))/AMKDRY(L)
      ENDDO

      !debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(pdbug)then
         print*,'AMKDRY(1:ibins)'
         print *,AMKDRY(1:ibins)
         print *,'WR(1:ibins)'
         print *,WR(1:ibins)
      endif
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cpja Initialize X() array of particle masses based on xk()
      DO L=1,ibins
         X(L)=xk(L)
      ENDDO

c
c Only solve when significant forcing is available
c
      maxtau=0.0d0
      do l=1,ibins
         maxtau=max(maxtau,abs(TAU(l)))
      enddo

      !debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(pdbug) then
      print*,'tau(1:ibins)'
      print *,tau(1:ibins)
      endif
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IF(ABS(maxtau).LT.TEPS)THEN
         DO L=1,ibins
            DO J=1,icomp
               AMK(L,J)=AMKD(L,J)
            ENDDO
            ANK(L)=ANKD(L)
         ENDDO
      ELSE
         !<step5.3> Try to fix the error of mass conservation
         ! during aqueous oxidation. Too little mass is used up
         ! (win, 7/24/06)
         IF ( MAXVAL(MOXD(:)) >  0d0 ) THEN
            IF( PDBUG ) PRINT *,'Mass_to_add_by_aqoxid_or_SOAcond'
            maddtot = 0d0
            DO L = 1, IBINS
               IF(TAU(L) >  0d0 ) THEN
                  MADD(L) = MOXD(L)
                  XADD(L) = MOXD(L) / ANKD(L)
!x                  IF( CSPECIES == SRTSO4 ) THEN
!x                     MADD(L) = MOXD * ANKD(L)  ! absolute condensing mass 
!x                     XADD(L) = MOXD            ! mass per particle 
!x                  ELSE IF ( CSPECIES == SRTOCIL ) THEN
!x                     MADD(L) = MOXD * AMKDRYSOL(L)
!x                     XADD(L) = MADD(L) / ANKD(L)
!x                  ELSE
!x                     PRINT *,'TMCOND ERROR : mass fixing not supported'
!x                  ENDIF
               ELSE
                  MADD(L) = 0D0
                  XADD(L) = 0D0
               ENDIF
               IF ( PDBUG ) PRINT *,L,madd(L), xadd(L)
               maddtot = maddtot + madd(L)
            ENDDO
         ENDIF

         DO L=1,ibins
            DO J=1,icomp
               AMK(L,J)=0.d0
            ENDDO
            ANK(L)=0.d0
         ENDDO
         WW=0.5d0
c        IF(TAU.LT.0.)WW=.5d0
c
c identify tophats and do lagrangian growth
c
         DO L=1,ibins
            IF(ANKD(L) .LT. NEPS2)GOTO 200 !skip if Number is effectively zero

            !if tau is zero, leave everything in same bin
            IF (TAU(L) .EQ. 0.) THEN
               ANK(L)=ANK(L)+ANKD(L)
               DO J=1,icomp
                  AMK(L,J)=AMK(L,J)+AMKD(L,J)
               ENDDO
            ENDIF
            IF (TAU(L) .EQ. 0.) GOTO 200

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !<step5.1> Add print for debugging (win, 4/10/06)
            if (pdbug) then
               print *, 'Identify_tophat_and_grow-BIN',L
               print *,'Starting_Nk(1:ibins)'
               print *, ANK(1:ibins)
               print *,'Starting_Mk(1:ibins,comp)'
               do j=1,icomp-1
                  print *,'comp',j
                  print *, AMK(1:ibins,j)
               enddo
            endif
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cpja Limiting AVG, the average particle size to lie within the size
Cpja bounds causes particles to grow or shrink arbitrarily and is
Cpja wreacking havoc with choosing condensational timesteps and
Cpja conserving mass.  I have turned them off.
c            AVG=MAX(X(L),MIN(X(L+1),AMKDRY(L)/(NEPS+ANKD(L)))) 
!try bring the above line back, win 4/10/06
!win 4/10/06            


            AVG=AMKDRY(L)/ANKD(L)
            XX=X(L)/AVG


#if  defined( TOMAS12 ) || defined( TOMAS15 )
            if(l.lt.ibins-1)then ! bin quadrupuling
               XI=.5d0 + XX*(2.5d0 - 2.0d0*XX)
              !XI<1 means the AVG falls out of bin bounds
               if (XI .LT. 1.d0) then
                  !W1 will have sqrt of negative number
                  write(*,*)'ERROR: tmcond - XI<1 for bin: ',L
                  write(*,*)'AVG is ',AVG
                  write(*,*)'Nk is ', ANKD(k)
                  write(*,*)'Mk are ', (AMKD(k,j),j=1,icomp)
                  write(*,*)'Initial N and M are: ',ANKD(L),AMKDRY(L)
                  errspot = .true.
                  RETURN
               endif
               W1 =SQRT(12.d0*(XI-1.d0))*AVG/4.0d0 ! cyhl 4.0=xk(k+1)/xk(k)
               W2 =(MIN(X(L+1)-AVG,AVG-X(L)))*2.0d0
            else ! final 2 bins mass*32
               XI=.5d0 + XX*(16.5d0 - 16.0d0*XX)
               if (XI .LT. 1.d0) then
                  !W1 will have sqrt of negative number
                  write(*,*)'ERROR: tmcond - XI<1 for bin: ',L
                  write(*,*)'lower limit is',X(L)
                  write(*,*)'AVG is ',AVG
                  write(*,*)'Nk is ', ANKD(L)
                  write(*,*)'Mk are ', (AMKD(k,j),j=1,icomp)
                  write(*,*)'Initial N and M are: ',ANKD(L),AMKDRY(L)
                  errspot = .true.
                  RETURN
               endif
               W1 =SQRT(12.d0*(XI-1.d0))*AVG/32.0d0 ! cyhl 32.0=xk(k+1)/xk(k)
               W2 =(MIN(X(L+1)-AVG,AVG-X(L)))*2.0d0
            endif
#else
            XI=.5d0 + XX*(1.5d0 - XX)
            !XI<1 means the AVG falls out of bin bounds

            if (XI .LT. 1.d0) then
               !W1 will have sqrt of negative number
               write(*,*)'ERROR: tmcond - XI<1 for bin: ',L
               write(*,*)'AVG is ',AVG
               write(*,*)'Nk is ', ANKD(L)
               write(*,*)'Mk are ', (AMKD(L,j),j=1,icomp)
               write(*,*)'Initial N and M are: ',ANKD(L),AMKDRY(L)
               errspot = .true.
               RETURN
            endif
            W1 =SQRT(12.d0*(XI-1.d0))*AVG
            W2 =MIN(X(L+1)-AVG,AVG-X(L))
#endif

            WTH=W1*WW+W2*(1.d0-WW)
            IF(WTH.GT.1.) then
               write(*,*)'WTH>1 in cond.f, bin #',L
               errspot = .true.
               RETURN
            ENDIF

            XU=AVG+WTH*.5d0
            XL=AVG-WTH*.5d0
c Ventilation added bin-by-bin
            TAU_L=TAU(l)*MAX(1.d0,VNTF(L)*VSW)
            IF(TAU_L/TAU(l).GT. 6.) THEN
               PRINT *,'TAU..>6.',TAU(l),TAU_L,VSW,L
            ENDIF
            IF(TAU_L.GT.TAU(l)) THEN 
               PRINT *,'TAU...',TAU(l),TAU_L,VSW,L
            ENDIF
! prior to 5/25/06 (win)
!            YU=DMDT_INT(XU,TAU_L,WR(L))
!            YUC=XU*AMKD(L,CSPECIES)/AMKDRY(L)+YU-XU
!            IF (YU .GT. X(ibins+1) ) THEN
!               YUC=YUC*X(ibins+1)/YU
!               YU=X(ibins+1)
!            ENDIF
!            YL=DMDT_INT(XL,TAU_L,WR(L)) 
!            YLC=XL*AMKD(L,CSPECIES)/AMKDRY(L)+YL-XL
!add new correction factor to YU and YL (win, 5/25/06)
            YU=DMDT_INT(XU,TAU_L,WR(L)) 
            YL=DMDT_INT(XL,TAU_L,WR(L)) 

            ! change to check MOXD of current bin (win, 10/3/08)
            IF( MOXD(L) == 0d0) THEN
            !Prior to 10/3/08 (win)
            !IF( MAXVAL(MOXD(:)) == 0D0 ) THEN
               C1=1.d0          !for so4cond call, without correction factor.
            ELSE
               C1 = XADD(L)*2.d0/(YU+YL-XU-XL)
            ENDIF
            C2 = C1 - ( C1 - 1.d0 ) * ( XU + XL )/( YU + YL )
          !prior to 10/2/08 (win)  
            YU = YU * C2
            YL = YL * C2
            ! Run into a problem that YU < XU creating YUC<0
            ! So let's limit the application of C2 to only if
            ! it does not result in YU < XU and YL < XL (win, 10/2/08)
ccc            IF(TAU_L > 0.d0) YU = max( YU*C2, XU )
ccc            IF(TAU_L > 0.d0) YL = max( YL*C2, XL )
 
!end part for fudging to get higher AVG 

            YUC=XU*AMKD(L,CSPECIES)/AMKDRY(L)+YU-XU
            IF (YU .GT. X(ibins+1) ) THEN
!               IF(.not.SPINUP(60.)) write(116,*)
!     &              'YU > Xk(30+1) ++++++++++++' !debug (win, 7/17/06)
               YUC=YUC*X(ibins+1)/YU
               YU=X(ibins+1)
!               errspot=.true.  !just try temp (win, 7/30/07)
            ENDIF
            YLC=XL*AMKD(L,CSPECIES)/AMKDRY(L)+YL-XL
            DYI=1.d0/(YU-YL)
            
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !<step5.2> Debug why there is extra mass added when called
            ! by aqoxid. (win, 5/10/06)
            if (pdbug) then
               print *, 'XU',XU,'YU',YU,'YUC',YUC,'c2',c2
               print *, 'XL',XL,'YL',YL,'YLC',YLC
            endif
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            !deal with tiny negative (win, 5/28/06)
            if(YUC.lt.0d0 .or. YLC.lt.0d0)then
               if(YLC.lt.0d0) YLC=0d0
               if(YUC.lt.0d0) then
                  YUC = 0d0
                  YLC = 0d0
               endif
               if(pdbug) print *,'Fudge negative YUC, YLC to zero'
            endif
c
c deal with portion of distribution that lies below lowest gridpoint
c
            IF(YL.LT.X(1))THEN
            
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !<step5.2> Debug step-by-step (win, 5/10/06)
            if (pdbug) print *,'YL<X(1)_Just_condensing_to_current_bin'
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
cpja Instead of the following, I will just add all new condensed
cpja mass to the same size bin
c               if ((YL/XL-1.d0) .LT. 1.d-3) then
c                  !insignificant growth - leave alone
c                  ANK(L)=ANK(L)+ANKD(L)
c                  DO J=1,icomp-1
c                     AMK(L,J)=AMK(L,J)+AMKD(L,J)
c                  ENDDO
c                  GOTO 200
c               else
c                  !subtract out lower portion
c                  write(*,*)'ERROR in cond - low portion subtracted'
c                  write(*,*) 'Nk,Mk: ',ANKD(L),AMKD(L,1),AMKD(L,2)
c                  write(*,*) 'TAU: ', TAU_L
c                  write(*,*) 'XL, YL, YLC: ',XL,YL,YLC
c                  write(*,*) 'XU, YU, YUC: ',XU,YU,YUC
c                  ANKD(L)=ANKD(L)*MAX(ZERO,(YU-X(1)))*DYI
c                  YL=X(1)
c                  YLC=X(1)*AMKD(1,CSPECIES)/AMKDRY(1)
c                  DYI=1.d0/(YU-YL)
c               endif
               ANK(L)=ANK(L)+ANKD(L)
               do j=1,icomp
                  if (J.EQ.CSPECIES) then
                     AMK(L,J)=AMK(L,J)+(YUC+YLC)*.5d0*ANKD(L)
                  else
                     AMK(L,J)=AMK(L,J)+AMKD(L,J)
                  endif
               enddo
               GOTO 200
            ENDIF
            IF(YU.LT.X(1))GOTO 200
c
c Begin remapping (start search at present location if condensation)
c
            IMN=1
            IF(TAU(l).GT.0.)IMN=L
            DO I=IMN,ibins
               !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               !<step5.2> Debug step-by-step (win, 5/10/06)
               if(pdbug) print *,'Now_remapping_in_bin',I
               !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               IF(YL.LT.X(I+1))THEN 
                  ![1] lower bound of new tophat in the current I bin
                  IF(YU.LE.X(I+1))THEN
                     ![2] upper bound of new tophat also in the current I bin
                     DN=ANKD(L)      ! DN = number from the bin L being remapped
                     do j=1,icomp
                        DM=AMKD(L,J)   
                        IF (J.EQ.CSPECIES) THEN
                           !Add mass from new tophat to the existing mass of bin I
                           AMK(I,J)=(YUC+YLC)*.5d0*DN+AMK(I,J)  
                           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           !<step5.2> Debug step-by-step (win, 5/10/06)
                           if (pdbug) then
                           print *,'CASE_1:_New_Tophat_in_a_single_bin'
                           print *,'SO4_from_tophat=',(YUC+YLC)*.5d0*DN
                           endif
                           !<step5.3> Check mass conservation (win, 7/24/06)
                           if(MAXVAL(moxd(:)).gt.0d0)then
                           delt1 = (YUC+YLC)*.5d0*DN-AMKD(L,J)-madd(L)
                           if( abs(delt1)/madd(L).gt.1d-6 .and.
     &                             madd(L).gt.1d-4)then
                              ! Just print out this for debugging 
                              IF(.not.SPINUP(60.) .and. pdbug ) then
!                                 write(116,*)'CASE1_mass_conserv_fix'
                                 write(116,13) L, madd(L), delt1 
 13                              FORMAT('CASE_1 Bin ',I2,' moxid ', 
     &                                   E13.5,' delta ',E13.5 )
!                                 errspot=.true. !just try temp (win, 7/30/07)
                              ENDIF
                              AMK(I,J) = AMK(I,J)-delt1 !fix the error
                           endif
                           endif
!                           endif !temporary 
                          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                           
                        ELSE
                           !For non-condensing, migrate the mass to bin I
                           AMK(I,J)=AMK(I,J)+DM 
                           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           !<step5.2> Debug step-by-step (win, 5/10/06)
                           if (pdbug) then
                           !print *,' Migrating_mass(',j,')',DM   !use this debugging line if there are more than seasalt+so4
                           print *,'Migrating_mass',DM
                           endif
                           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                           
                        ENDIF
                     enddo
                     !Add number of old bin to ANK (which is blank for the first loop of bin I)
                     ANK(I)=ANK(I)+DN
                     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     !<step5.2> (win, 5/10/06)
                     if(pdbug) print*,'Migrating_number',DN
                     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  ELSE  
                    ![3] upper bound of new tophat grow beyond the upper bound of bin I
                     DN=ANKD(L)*(X(I+1)-YL)*DYI !DN= proportion of the number from tophat that still stays in the bin I
                     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     !<step5.2> (win, 5/10/06)
                     if ( pdbug) then
                     print*,'Case_2:_Tophat_cross_bin_boundary'
                     print *,'Number_that_remain_in_low_bin',DN
                     endif
                     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                     !<step5.3> For fixing mass conserv problem (win, 7/24/06)
                     macc=0d0

                    do j=1,icomp
                       !DM= proporation of the mass that is still in bin I
                        DM=AMKD(L,J)*(X(I+1)-YL)*DYI 
                        IF (J.EQ.CSPECIES) THEN
                           !XP= what would have grown to be X(I+1) 
                           XP=DMDT_INT(X(I+1),-1.0d0*TAU_L,WR(L)) 
                           YM=XP*AMKD(L,J)/AMKDRY(L)+X(I+1)-XP
                           !add the condensing mass to the existing sulfate of bin I
                           AMK(I,J)=DN*(YM+YLC)*0.5d0+AMK(I,J)
                           !<step5.3>Accumulating the condensing mass for error check (win, 7/24/06)
                           macc = macc + DN*(YM+YLC)*0.5d0
                           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           !<step5.2> (win, 5/10/06)
                           if(pdbug)then
                           print *,'XP',XP,'YM',YM
                           print *,'Cond_TophatLowEnd',DN*(YM+YLC)*0.5d0
                           endif
                           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        ELSE
                           !Add DM to AMK (which is blank for the first loop of bin I)
                           AMK(I,J)=AMK(I,J)+DM
                           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           if(pdbug) print*,'Other___in_low_end',DM
                           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        ENDIF
                     enddo
                     ANK(I)=ANK(I)+DN ! Add DN number to ANK (which is blank for the first loop of bin I)
                     ! Remapping loop from bin I+1 to bin30
                     DO K=I+1,ibins   
                     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     if(pdbug) print *,'Spreading_to_bin',K
                     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        IF(YU.LE.X(K+1))GOTO 100 
                        ![4] Found the bin where the high end of the tophat is in --> do the final loop
                        
                        ![5.1] This part for distributing to the bins in between 
                        !      the original and the furthest bin that growing occurs


                        !Use width of bin K to proportionate number from old bin wrt. to the top hat (YU-YL)                        
                        DN=ANKD(L)*(X(K+1)-X(K))*DYI  

                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if(pdbug) then
                           print *,'Number_migrated',DN
                        endif
                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        do j=1,icomp
                           !Proportion of old-bin mass that falls in this current bin K
                           DM=AMKD(L,J)*(X(K+1)-X(K))*DYI 
                           IF (J.EQ.CSPECIES) THEN
                              XP=DMDT_INT(X(K),-1.0d0*TAU_L,WR(L)) !what would have grown to be X(k)
                              YM=XP*AMKD(L,J)/AMKDRY(L)+X(K)-XP !what would have grown to be X(k) but just for sulfate
                              AMK(K,J)=DN*1.5d0*YM+AMK(K,J)    ! A factor of 1.5 is from averaging (YM+2*YM)
                              !<step5.3> Accumulating condensing mass for error check (win, 7/24/06)
                              macc = macc+DN*1.5d0*YM
                              !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                              !<step5.2> (win, 5/10/06)
                              if(pdbug)then
                                 print *,'XP',XP,'YM',YM
                                 print *,'Cond_mass_spread',DN*1.5d0*YM
                              endif
                              !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           ELSE
                              AMK(K,J)=AMK(K,J)+DM    !Add migrating mass of non-condensing species
                              !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                              if(pdbug) print*,'No-cond_mass_migrate',DM
                              !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           ENDIF
                        enddo
                        ANK(K)=ANK(K)+DN  !Add migrating number to the exising number of bin K
                     ENDDO
                     !This STOP is for when there's excessive growth over bin30
                     STOP 'Trying to put stuff in bin ibins+1'
                    
 100                 CONTINUE
                     ![5.2] Final section that the tophat grows to.
                     DN=ANKD(L)*(YU-X(K))*DYI  
                     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     if(pdbug) then
                        print *,'Found_right_edge_for_tophat'
                        print *,'Number_migrated',DN
                     endif
                     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     do j=1,icomp
                        DM=AMKD(L,J)*(YU-X(K))*DYI  ! proportion of old mass that gets to this furthest bin.
                        IF (J.EQ.CSPECIES) THEN
                           XP=DMDT_INT(X(K),-1.0d0*TAU_L,WR(L))   !what would have grown to be X(k)
                           YM=XP*AMKD(L,J)/AMKDRY(L)+X(K)-XP !=XP for just sulfate
                           AMK(K,J)=DN*(YUC+YM)*0.5d0+AMK(K,J) !add condensing mass to existing sulfate of bin K
                           !<step5.3>Accumulating condensing mass for error check (win, 7/24/06)
                           macc = macc+DN*(YUC+YM)*0.5d0
                           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           !<step5.2> (win, 5/10/06)
                           if(pdbug)then
                           print *,'XP',XP,'YM',YM
                           print *,'Cond_mass_spread_final',
     &                          DN*(YM+YUC)*0.5d0
                           endif
                           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        ELSE
                           AMK(K,J)=AMK(K,J)+DM  !This adds the migrating mass to the exising mass of non-condensing species
                           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           if(pdbug) print*,'No-cond_mass_migrated',DM
                           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        ENDIF
                     enddo
                     ANK(K)=ANK(K)+DN   !This adds the migrating number to the existing number of bin K

                     !<step5.3> Check mass conservation (win, 7/24/06)
                     if(MAXVAL(moxd(:)).gt.0d0)then
                        delt2 = 0d0
                        delt2 = macc-AMKD(L,CSPECIES)-madd(L)
                        if(abs(delt2)/ madd(L) > 1e-6)then
                           if( madd(L) > 10.d0 .and. 
     &                          abs(delt2)/ madd(L) > 15d-2 ) then
!                              print *,'TMCOND ERROR: mass condensation',
!     &                          'discrep >15% during aqoxid or SOAcond'
                              IF(.not.SPINUP(60.))  THEN
 14                              FORMAT('CASE_2 Bin',I2,' moxid',F7.1,
     &                                  ' delta',F7.1 )
                                 write(116,14) L, madd(L),delt2
!                                 write(116,*)'CASE_2_mass_not_conserve'
!                                 write(116,*)'For_bin',L,'moxid',madd(L)
!     &                                ,'delta',delt2
                              ENDIF
                              errspot=.true. !just try temp (win, 7/30/07)
                           endif !significant mass add (10 kg) - then print error.
                           !<step5.3> Fix the problem of mass not conserved
                           !in case of aqueous oxidation by find the missing mass
                           !and spread them equally into the bins that the final 
                           !tophat has grown to. (win, 7/24/06)
                           xtra  = 0d0
                           dummy = 0d0
                           do kk = I,K
!                              AMK(kk,CSPECIES) = AMK(kk,CSPECIES)-
!     &                                           delt2/(K-I+1) 
                              dummy = AMK(kk,CSPECIES) - 
     &                               ( delt2/(K-I+1) + xtra )
                              if(dummy < 0.d0 )then
                                 xtra = xtra + delt2/(K-I+1)
                              else                                 
                                 AMK(kk,CSPECIES) = dummy
                                 xtra = 0.d0
                              endif
                           enddo
                        endif   !error>treshold
                     endif      !moxd>0

                  ENDIF  !YU.LE.X(I+1)
                  GOTO 200
               ELSE    !YL > X(I+1)
                  IF(I == IBINS .and.(madd(L)/maddtot)> 1.5d-1) THEN
 11                  FORMAT( 'Tophat>Xk(31) at bin ',I3,' loosing ',
     &                     E13.5,' kg = ',F5.1,'%')
                     if(MAXVAL(moxd(:)) > 0d0) then
                        print 11, L, madd(L),(madd(L)/maddtot)*1.d2
                        write(116,11) L, madd(L),(madd(L)/maddtot)*1.d2
                        write(117,*) madd(L)  !for accumulating mass loss
!                     PRINT *,'Tophat > Xk(31): growth over bin30,Loss%'
!                     if(moxd >0d0)print *,madd(L),(madd(L)/maddtot)*1.d2
!jrp                        errspot = .true.
                     endif
                  ENDIF
               ENDIF   !YL.LT.X(I+1)
            ENDDO !I loop
 200        CONTINUE
         ENDDO    !L loop
      ENDIF

      !Signal error out to so4cond so the run can stop in aerophys and show i,j,l (win, 4/12/06)
      pdbug = errspot

      ! Return to call subroutine
      RETURN
      END SUBROUTINE TMCOND

!------------------------------------------------------------------------------
      
      SUBROUTINE AERODIAG( PTYPE, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL )
!
!******************************************************************************
!  Subroutine AERODIAG saves changes to the appropriate diagnostic arrays
!  (win, 7/23/07)
!
!  Variables
!  (1 ) PTYPE  (INTEGER) : Number assigned to each dianostic
!  (2 ) I,J,L  (INTEGER) : Location of grid box
!
!  NOTES: 
!  (1 ) Add error checking IT_IS_NAN (win, 7/31/07)
!  (2 ) Add new diagnostic for TMS-SOA (win, 9/25/07)
!  (3 ) Add new diagnostic TOMAS-3D to save rates of a selected process
!       in 3-D (I,J,L).  At this point just do one for NK1 nucleation (win, 5/21/08)
!  (4 ) Modified AD61 to use with 2 tracers and adjust dimension accordingly
!       (win, 10/06/08)
!*****************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,    ONLY :  AD60_COND, AD60_COAG, AD60_NUCL
      USE DIAG_MOD,    ONLY :  AD60_AQOX, AD60_ERROR, AD60_SOA
      USE DIAG_MOD,    ONLY :  AD61,      AD61_INST  
      USE ERROR_MOD,   ONLY :  IT_IS_NAN
      USE TIME_MOD,    ONLY : GET_TS_CHEM

      USE CMN_SIZE_MOD      ! Size parameters
      USE CMN_DIAG_MOD      ! ND60

      ! Arguments
      INTEGER               :: PTYPE,    I,      J,      L 

      ! Local variables
      INTEGER               :: K, JS
      REAL*4                :: ADXX(IBINS*(ICOMP-IDIAG+1))
      REAL*4,    SAVE       :: ACCUN, ACCUM(2)
      LOGICAL,   SAVE       :: FIRST = .TRUE.

      real*4                :: tempsum
      REAL*4                :: DTCHEM

!sfarina
      REAL*8                :: Nk(IBINS), Nkd(IBINS)
      REAL*8                :: Mk(IBINS, ICOMP)
      REAL*8                :: Mkd(IBINS,ICOMP)
      REAL*4,INTENT(IN)     :: BOXVOL
      !=================================================================
      ! AERODIAG begins here!
      !=================================================================

      ! PTYPE = 7 is for ND61  --- NOW use for Nucleation at tracer NK1 
      !  Note: This is created to look at 3-D rate for a selected process
      !        Right now (5/21/08) I created this to watch NUCLEATION rate
      !        We can't afford to save all 30-bin and all mass component
      !        in all (I,J,L), thus this is created. (win, 5/21/08)
      IF ( PTYPE == 7 ) THEN 
         IF ( L <= LD61 ) THEN
            DTCHEM = GET_TS_CHEM() * 60d0   ! chemistry time step in sec
            AD61(I,J,L,1) = AD61(I,J,L,1)  + 
     &           ( NK(1) - NKD(1) )/ DTCHEM / BOXVOL  ! no./cm3/sec
            AD61_INST(I,J,L,1) =  ( NK(1) - NKD(1) ) /DTCHEM / BOXVOL ! no./cm3/sec

            !IF(i==39 .and. j==29 ) then 
             !  if ( AD61_INST(I,J,L) .gt. 1e18)  write(6,*) '*********', 
!     &           'AD61_INST(',I,J,L,')', AD61_INST(I,J,L)
            !endif
         ENDIF
      ELSE ! PTYPE = 1-6 is for ND60

      ADXX(:) = 0d0
      IF ( FIRST ) THEN
         ACCUN = 0e0
         ACCUM(:) = 0e0
         FIRST = .FALSE.
      ENDIF

      ! Debug: check error fixed accumulated at each step
!      IF ( I == 1 .and. J == 1 .and. L == 1 .and. PTYPE == 2) then
!         print *, 'Accumulated diagnostic for ND60 #',PTYPE,' at',i,j,l
!         print *, '   Number :',ACCUN
!         print *, '   Sulf   :',ACCUM(1)
!         print *, '   NaCl   :',ACCUM(2)
!      ENDIF

      IF ( L <= LD60 ) THEN

      SELECT CASE ( PTYPE )
         
      CASE ( 1 )                ! Condensation diagnostic
         ADXX(:) = AD60_COND(1,J,L,:)
         
      CASE ( 2 )                ! Coagulation diagnostic
         ADXX(:) = AD60_COAG(1,J,L,:)
         
      CASE ( 3 )                ! Nucleation diagnostic
         ADXX(:) = AD60_NUCL(1,J,L,:)
         
      CASE ( 4 )                ! Aqueous oxidation diagnostic
         ADXX(:) = AD60_AQOX(1,J,L,:)
         
      CASE ( 5 )                ! Error fudging diagnostic
         ADXX(:) = AD60_ERROR(1,J,L,:)
         
      CASE ( 6 )                ! SOA condensation diagnostic
         ADXX(:) = AD60_SOA(1,J,L,:)

      END SELECT

!delete_this_after_dbg_(win,8/1/07)
!      if(i==6.and.j==34.and.l==28) then 
!         print *, I,J,L,'--------------------'
!         print *,'Nk         Nkd'
!         do k=1,ibins
!         print *,NK(k), nkd(k)
!         enddo
!         print *,'Mk(:,1)    Mkd(:,1)'
!         do k=1,ibins
!         print *,mk(k,1), mkd(k,1)
!         enddo
!         print *,'Mk(:,2)    Mkd(:,2)'
!         do k=1,ibins
!         print *,mk(k,2), mkd(k,2)
!         enddo
!      endif
!delete to this line ---------------------         

      ! Change of aerosol number    
      DO K = 1, IBINS
         ADXX(K) =  ADXX(K) + NK(K) - NKD(K)
!         IF ( PTYPE == 2 ) ACCUN = ACCUN + NK(K) - NKD(K)
      ENDDO
      IF ( IT_IS_NAN(ACCUN)) print *,'AERODIAG: Nan',I,J,L

      ! Change of aerosol mass
      DO JS = 1, ICOMP-IDIAG
         tempsum = 0e0
      DO K = 1, IBINS
         ADXX(JS*IBINS+K) = ADXX(JS*IBINS+K) + MK(K,JS) - MKD(K,JS)
!         tempsum = tempsum + MK(K,JS) - MKD(K,JS)
!         IF (PTYPE == 2 ) ACCUM(JS) = ACCUM(JS) + MK(K,JS) - MKD(K,JS)
      ENDDO
      !temp-----------
!      if(ptype == 6 .or. ptype ==4) then
!        print *, 'Component',JS,'sumchange',tempsum
!        print *, '-----------------------'
!      endif
      !temp-----------
      ENDDO
!      IF ( IT_IS_NAN(ACCUM(1))) print *,'ADIAG: Nan',I,J,L

      ! Put the updated values back into the diagnostic arrays
      SELECT CASE ( PTYPE )      
      CASE ( 1 )                ! Condensation diagnostic
         AD60_COND(1,J,L,:) = ADXX(:)
         
      CASE ( 2 )                ! Coagulation diagnostic
         AD60_COAG(1,J,L,:) = ADXX(:)
         
      CASE ( 3 )                ! Nucleation diagnostic
         AD60_NUCL(1,J,L,:) = ADXX(:)
         
      CASE ( 4 )                ! Aqueous oxidation diagnostic
         AD60_AQOX(1,J,L,:) = ADXX(:)
         
      CASE ( 5 )                ! Error fudging diagnostic
          AD60_ERROR(1,J,L,:) = ADXX(:)
         
      CASE ( 6 )                ! SOA condensation diagnostic
          AD60_SOA(1,J,L,:) = ADXX(:)
      END SELECT

      ENDIF
      ! Debug: check error fixed accumulated at each step
!      IF ( I == 3 .and. J == 41 .and. L == 30 .and. PTYPE == 5) then
!         print *, 'Accumulated diagnostic for ND60 #',PTYPE,' at',i,j,l
!         print *, '   Number :',ACCUN
!         print *, '   Sulf   :',ACCUM(1)
!         print *, '   NaCl   :',ACCUM(2)
!         print *, ' Nk',Nk(:)
!         print *, ' Nkd',Nkd(:)
!         print *, ' Mk',Mk(:,1)
!         print *, ' Mkd',Mkd(:,1)
!         print *, ' Mk',Mk(:,2)
!         print *, ' Mkd',Mkd(:,2)

!      ENDIF

      ! Debug: check error fixed accumulated at each step
      IF ( I == IIPAR .and. J == JJPAR .and. L == LLPAR 
     &     .and. PTYPE == 2 ) then
!         print *, 'Accumulated diagnostic for ND60 #',PTYPE,' at',i,j,l
         print *, ' Accumulated Coagulation'
         print *, '   Number :',ACCUN
         print *, '   Sulf   :',ACCUM(1)
         print *, '   NaCl   :',ACCUM(2)
      ENDIF

      ENDIF ! If (PTYPE == 7)

      ! Return to calling subroutine
      END SUBROUTINE AERODIAG

!------------------------------------------------------------------------------

      SUBROUTINE INIT_TOMAS( am_I_Root, Input_Opt, RC )

!
!******************************************************************************
!  Subroutine INIT_TOMAS intializes variables for TOMAS microphysics based on
!  switches from input.geos, e.g. what aerosol species to simulate.(win, 7/9/07)
!
!  NOTES:
!   (1 ) Add calculation for bin boundaries array Xk using Mo = lower mass 
!         bound for first size bin (kg) (win, 7/18/07)
!   (2 ) IDIAG is introduced.  IDIAG=2 for NH4 and H2O (win, 9/26/08)
!  23 Apr 2013 - R. Yantosca - Now get unit # with findFreeLUN
!  30 Jan 2014 - R. Yantosca - Now read big files from shared disk space
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,          ONLY : ALLOC_ERR, ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE inquireMod,         ONLY : findFreeLUN
!------------------------------------------------------------------------------
! Prior to 1/30/14:
! Now use Input_Opt for the data directory (bmy, 1/30/14)
!      USE DIRECTORY_MOD,      ONLY: DATA_DIR
!------------------------------------------------------------------------------
      USE TRACERID_MOD,       ONLY : IDTSF1, IDTSS1, IDTECIL1, IDTECOB1
      USE TRACERID_MOD,       ONLY : IDTOCIL1, IDTOCOB1, IDTDUST1
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !LOCAL VARIABLES:
!
      INTEGER              :: AS, K, I,J,L, dum1, dum2, dum3, LUN
      REAL*8               :: Mo
!------------------------------------------------------------------------------
! Prior to 1/30/14:
! Make all character strings 255 chars wide (bmy, 1/30/14)
!      character*80         :: filename
!      character*80         :: fname(4)
!------------------------------------------------------------------------------
      CHARACTER(LEN=255)   :: filename
      CHARACTER(LEN=255)   :: fname(4)     
      CHARACTER(LEN=255)   :: DATA_DIR

      !=================================================================
      ! INIT_TOMAS begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Now read large TOMAS input files from a common disk directory
      ! (bmy, 1/30/14)
      DATA_DIR = TRIM( Input_Opt%DATA_DIR_1x1 ) // 'TOMAS_201402/'

#if    defined( TOMAS40 )
      Mo = 1.0d-21*2.d0**(-10)     
#elif  defined( TOMAS15 )
      Mo = 1.0d-21*4.d0**(-3)     
#else
      Mo = 1.0d-21
#endif
      
      ICOMP = 0 
      IDIAG = 0
      K = 0
c      IF( LSULF30 )  THEN 
      IF (IDTSF1 > 0) THEN
         ICOMP = ICOMP + 1
c        SRTSO4 = ICOMP
      ENDIF
c      IF( LSALT30 )  THEN
      IF ( IDTSS1 > 0 ) THEN
         ICOMP = ICOMP + 1
c         SRTNACL = ICOMP
      ENDIF
C      IF( LCARB30 )  THEN
      IF ( IDTECIL1 > 0 .AND. IDTECOB1 > 0 .AND. 
     &     IDTOCIL1 > 0 .AND. IDTOCOB1 > 0 ) THEN
         ICOMP = ICOMP + 1
c         SRTECIL = ICOMP
         ICOMP = ICOMP + 1
c         SRTECOB = ICOMP
         ICOMP = ICOMP + 1
c         SRTOCIL = ICOMP
         ICOMP = ICOMP + 1
c         SRTOCOB = ICOMP
      ENDIF
c      IF( LDUST30 )  THEN
      IF ( IDTDUST1 > 0 ) THEN
         ICOMP = ICOMP + 1
c         SRTDUST = ICOMP
      ENDIF
      
      ! Have to add one more for aerosol water
      IF( ICOMP > 1 ) THEN 
         ICOMP = ICOMP + 1
         IDIAG = IDIAG + 1
c         SRTNH4 = ICOMP

         ICOMP = ICOMP + 1
         IDIAG = IDIAG + 1
c         SRTH2O = ICOMP
      ENDIF
      print *, 'In init_TOMAS, ICOMP = ', ICOMP
      print *, 'In init_TOMAS, IBINS = ', IBINS

      !=================================================================
      ! Allocate arrays 
      !=================================================================

      ALLOCATE( MOLWT( ICOMP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MOLWT [TOMAS]' )
      MOLWT(:) = 0d0

      !=================================================================
      ! Calculate aerosol size bin boundaries (dry mass / particle)
      !=================================================================

#if  defined( TOMAS12 ) || defined( TOMAS15 )
      DO K = 1, IBINS + 1
          if(k.lt.ibins)then
             xk(k)=Mo * 4.d0**(k-1) !mass quadrupling
          else
             xk(k)=xk(k-1) * 32.d0
          endif
      ENDDO
#else
      DO K = 1, IBINS + 1
         Xk( k ) = Mo * 2.d0 ** ( K-1 )
      ENDDO
#endif


      DO K = 1, IBINS
         AVGMASS( k ) = sqrt(Xk(k)*Xk(k+1))
      ENDDO

      DO J = 1, ICOMP
         IF ( J == SRTSO4 ) THEN
            MOLWT(J) = 98.0
         ELSE IF ( J == SRTNACL ) THEN
            MOLWT(J) = 58.5 
         ELSE IF ( J == SRTH2O ) THEN
            MOLWT(J) = 18.0
         ELSE IF ( J == SRTECIL ) THEN 
            MOLWT(J) = 12.0
         ELSE IF ( J == SRTECOB ) THEN 
            MOLWT(J) = 12.0
         ELSE IF ( J == SRTOCIL ) THEN 
            MOLWT(J) = 12.0
         ELSE IF ( J == SRTOCOB ) THEN 
            MOLWT(J) = 12.0
         ELSE IF ( J == SRTDUST ) THEN
            MOLWT(J) = 100.0
         ELSE IF ( J == SRTNH4 ) THEN
            MOLWT(J) = 18.0
         ELSE
            PRINT *,'INIT_TOMAS ERROR: Modify code for more species!!'
            CALL ERROR_STOP('INIT_TOMAS','Modify code for new species')
         ENDIF            
      ENDDO

      !=================================================================
      ! Create a look-up table for activating bin and scavenging 
      ! fraction as a function of chemical composition.
      !=================================================================
      IF ( IBINS == 12 .OR. IBINS == 15 ) THEN
!-----------------------------------------------------------------------------
! Prior to 1/30/14:
! Now read big files from a shared data directory (bmy, 1/30/14)
!         fname(1) = 'binact02_12.dat'
!         fname(2) = 'binact10_12.dat'
!         fname(3) = 'fraction02_12.dat'
!         fname(4) = 'fraction10_12.dat'
!-----------------------------------------------------------------------------
         fname(1) = TRIM( DATA_DIR ) // 'binact02_12.dat'
         fname(2) = TRIM( DATA_DIR ) // 'binact10_12.dat'
         fname(3) = TRIM( DATA_DIR ) // 'fraction02_12.dat'
         fname(4) = TRIM( DATA_DIR ) // 'fraction10_12.dat'

      ELSE IF ( IBINS == 30 .OR. IBINS == 40 ) THEN
!-----------------------------------------------------------------------------
! Prior to 1/30/14:
! Now read big files from a shared data directory (bmy, 1/30/14)
!         fname(1) = 'binact02.dat'
!         fname(2) = 'binact10.dat'
!         fname(3) = 'fraction02.dat'
!         fname(4) = 'fraction10.dat'
!-----------------------------------------------------------------------------
         fname(1) = TRIM( DATA_DIR ) // 'binact02.dat'
         fname(2) = TRIM( DATA_DIR ) // 'binact10.dat'
         fname(3) = TRIM( DATA_DIR ) // 'fraction02.dat'
         fname(4) = TRIM( DATA_DIR ) // 'fraction10.dat'
      END IF

      CALL READBINACT  ( fname(1), BINACT1   )
      CALL READBINACT  ( fname(2), BINACT2   )
      CALL READFRACTION( fname(3), FRACTION1 )
      CALL READFRACTION( fname(4), FRACTION2 )

      !initialize yu lookup table
      call READJIMN5D( am_I_Root, Input_Opt, RC )  ! yu nucleation inputs

      ! Find a free file LUN
      LUN = findFreeLUN()

      ! Read cosmic ray ion input file
      FILENAME = TRIM( DATA_DIR ) // 'IonPairs1GV.4x5'
      OPEN( unit=LUN, FILE=TRIM( FILENAME ), FORM='FORMATTED',
     &      STATUS='OLD' )

      DO L=1,9
      DO J=1,46
      DO I=1,72
         READ( LUN ,'(I5,I5,I5,E10.3)') dum1,dum2,dum3,
     &                                  cosmic_ions(I,J,L)
         if (I.eq.50.and.J.eq.20.and.L.eq.5)then
            print*,'ion test',cosmic_ions(I,J,L)
         endif
      ENDDO
      ENDDO
      ENDDO
      
      ! Close file
      CLOSE( LUN )

      ! Return to calling program
      END SUBROUTINE INIT_TOMAS

!------------------------------------------------------------------------------

      SUBROUTINE READBINACT( INFILE, BINACT )

      USE inquireMod, ONLY : findFreeLun

      CHARACTER*80 INFILE
      INTEGER INNUM, II, JJ, KK, LUN
      INTEGER BINACT(101,101,101)

      ! Find a free file LUN
      INNUM = findFreeLun()

 1    FORMAT(I2)
      OPEN(UNIT=INNUM,FILE=INFILE,FORM='FORMATTED',STATUS='OLD')
      DO II=1,101
      DO JJ=1,101
      DO KK=1,101
         READ(INNUM,1) BINACT(KK,JJ,II)
         IF (BINACT(KK,JJ,II).eq.0) BINACT(KK,JJ,II)=IBINS + 1
      ENDDO
      ENDDO
      ENDDO
      CLOSE(INNUM)
      RETURN
      END SUBROUTINE READBINACT
      
!------------------------------------------------------------------------------

      SUBROUTINE READFRACTION(INFILE,FRACTION)

      USE inquireMod, ONLY : findFreeLun

      CHARACTER*80 INFILE
      INTEGER INNUM, II, JJ, KK
      REAL*8 FRACTION(101,101,101)

      ! Find a free file LUN
      INNUM = findFreeLun()

 1    FORMAT(F6.5)
      OPEN(UNIT=INNUM,FILE=INFILE,FORM='FORMATTED',STATUS='OLD')
      DO II=1,101
      DO JJ=1,101
      DO KK=1,101
         READ(INNUM,1) FRACTION(KK,JJ,II)
         IF (FRACTION(KK,JJ,II).GT.1.) FRACTION(KK,JJ,II)=0.
      ENDDO
      ENDDO
      ENDDO
      CLOSE(INNUM)
      RETURN
      END SUBROUTINE READFRACTION

!------------------------------------------------------------------------------

      SUBROUTINE GETFRACTION
     &   ( I, J, L, N, LS, FRACTION, SOLFRAC, STT )
!
!******************************************************************************
!  Subroutine GETFRACTION calculate the mass fraction of each soluble component 
!  i.e. SO4, sea-salt, hydrophilic OC to use as inputs for a lookup table of
!  activating bin and scavenging fraction. (win, 9/10/07)
!
!  Variable comments:
!  I, J, L          : Grid box index
!  N                : Tracer ID 
!  LS               : True = Large-scale (stratiform) precip, False = convectiv
!  FRACTION         : Scavenging fraction of the given grid box
!  SOLFRAC          : Soluble mass fraction of the aerosol popultion of the 
!                     given grid box
!
!  NOTES:
!******************************************************************************

      ! Reference to F90 modules
c      USE LOGICAL_MOD, ONLY  : LSULF30, LSALT30
c      USE LOGICAL_MOD, ONLY  : LCARB30, LDUST30
      USE TRACERID_MOD, ONLY : IDTECIL1, IDTOCIL1, IDTOCOB1, IDTECOB1
      USE TRACERID_MOD, ONLY : IDTNK1,   IDTSF1,   IDTSS1,   IDTDUST1
      USE CMN_SIZE_MOD      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: I, J, L, N
      LOGICAL, INTENT(IN)    :: LS
      REAL*8,  INTENT(OUT)   :: FRACTION, SOLFRAC
      REAL*8,  INTENT(INOUT) :: STT(:,:,:,:)            ! Tracer array

      ! Local variables
      REAL*4                 ::  MECIL, MOCIL, MOCOB, MSO4, MNACL, MTOT
      REAL*4                 ::  MECOB, MDUST
      REAL*4                 ::  XOCIL, XSO4, XNACL
      INTEGER                ::  ISO4, INACL, IOCIL
      INTEGER                ::  GETBINACT
      INTEGER                ::  BIN
      INTEGER                ::  OFFSET

      !=================================================================
      ! GETFRACTION begins here
      !=================================================================

      BIN = N - IDTNK1 + 1
      IF ( BIN > IBINS ) THEN
         BIN = MOD( BIN, IBINS )
         IF ( BIN == 0 ) BIN = IBINS
      ENDIF
      
      MECIL = 0.E0
      MOCIL = 0.E0
      MOCOB = 0.E0
      MSO4  = 0.E0
      MNACL = 0.E0
      MDUST = 0.E0

      IF ( IDTECIL1 > 0 .AND.IDTOCIL1 > 0 .AND. IDTOCOB1 > 0 ) THEN
         MECIL = STT(I,J,L,IDTECIL1-1+BIN)
         MOCIL = STT(I,J,L,IDTOCIL1-1+BIN)
         MOCOB = STT(I,J,L,IDTOCOB1-1+BIN)
      ENDIF
      IF ( IDTDUST1 > 0 ) MDUST = STT(I,J,L,IDTDUST1-1+BIN)
      IF ( IDTSF1 > 0 ) MSO4  = STT(I,J,L,IDTSF1-1+BIN) * 1.2 !account for ammonium sulfate
      IF ( IDTSS1 > 0 ) MNACL = STT(I,J,L,IDTSS1-1+BIN)
      MTOT  = MECIL + MOCIL + MOCOB + MSO4 + MNACL + MDUST + 1.e-20
      XOCIL = MOCIL / MTOT
      XSO4  = MSO4  / MTOT
      XNACL = MNACL / MTOT
      ISO4  = MIN(101, INT(XSO4*100)+1)
      INACL = MIN(101, INT(XNACL*100)+1)
      IOCIL = MIN(101, INT(XOCIL*100)+1)

      !==========================================================
      ! subroutine was written considering bin 1 is 10nm
      ! in TOMAS-40, bin 1 is 1nm and bin 11 is 10nm
      !==========================================================
#if    defined( TOMAS40 )
      OFFSET = 10
#elif  defined( TOMAS15 )
      OFFSET = 3
#else
      OFFSET = 0
#endif

      IF ( LS ) THEN 
         GETBINACT = BINACT1(ISO4, INACL, IOCIL) + OFFSET
      ELSE 
         GETBINACT = BINACT2(ISO4, INACL, IOCIL) + OFFSET
      ENDIF

      if((GETBINACT.lt.0).or.(GETBINACT.gt.50))then
         print*,'BINACT ERROR GETBINACT=',GETBINACT
         stop
      endif

!      print*,'N, BINACT = ',N,GETBINACT
      
      IF ( GETBINACT > BIN ) THEN
         FRACTION = 0. !NOT ACTIVATED
      ELSE IF ( GETBINACT == BIN ) THEN
         IF ( LS ) THEN 
            FRACTION = FRACTION1(ISO4, INACL, IOCIL ) !PARTLY ACTIVATED
         ELSE
            FRACTION = FRACTION2(ISO4, INACL, IOCIL ) !PARTLY ACTIVATED
         ENDIF
      ELSE
         FRACTION = 1. !ALL ACTIVATED
      ENDIF

      ! Calculate the soluble fraction of mass
      MECOB = 0.E0
      IF ( IDTECOB1 > 0 ) MECOB = STT(I,J,L,IDTECOB1-1+BIN)
      SOLFRAC = MTOT / ( MTOT + MECOB ) 

      ! Return to calling program
      RETURN
      END SUBROUTINE GETFRACTION
      
!------------------------------------------------------------------------------

      SUBROUTINE GETACTBIN ( I, J, L, N, LS, BINACT, State_Chm )

!
!******************************************************************************
!  Subroutine GETFRACTION calculate the mass fraction of each soluble component 
!  i.e. SO4, sea-salt, hydrophilic OC to use as inputs for a lookup table of
!  activating bin and scavenging fraction. (win, 9/10/07)
!
!  Variable comments:
!  I, J, L          : Grid box index
!  N                : Tracer ID 
!  LS               : True = Large-scale (stratiform) precip, False = convectiv
!  FRACTION         : Scavenging fraction of the given grid box
!  SOLFRAC          : Soluble mass fraction of the aerosol popultion of the 
!                     given grid box
!
!  NOTES:
!******************************************************************************

      ! Reference to F90 modules
      USE CMN_SIZE_MOD
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE TRACERID_MOD, ONLY : IDTECIL1, IDTOCIL1, IDTOCOB1, IDTECOB1
      USE TRACERID_MOD, ONLY : IDTNK1,   IDTSF1,   IDTSS1 ,  IDTDUST1

      ! Arguments
      INTEGER,        INTENT(IN)    :: I, J, L, N
      LOGICAL,        INTENT(IN)    :: LS
      INTEGER,        INTENT(OUT)   :: BINACT
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object

      ! Local variables
      REAL*4                 ::  MECIL, MOCIL, MOCOB, MSO4, MNACL, MTOT
      REAL*4                 ::  MECOB, MDUST
      REAL*4                 ::  XOCIL, XSO4, XNACL
      INTEGER                ::  ISO4, INACL, IOCIL
      INTEGER                ::  BIN
      INTEGER                ::  OFFSET

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER   :: STT(:,:,:,:)

      !=================================================================
      ! GETACTBIN begins here
      !=================================================================

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers


      BIN = N - IDTNK1 + 1
      IF ( BIN > IBINS ) THEN
         BIN = MOD( BIN, IBINS )
         IF ( BIN == 0 ) BIN = IBINS
      ENDIF

      MECIL = 0.E0
      MOCIL = 0.E0
      MOCOB = 0.E0
      MSO4  = 0.E0
      MNACL = 0.E0

      IF ( IDTECIL1 > 0 .AND.IDTOCIL1 > 0 .AND. IDTOCOB1 > 0 ) THEN
c      IF (LCARB30) THEN
         MECIL = STT(I,J,L,IDTECIL1-1+BIN)
         MOCIL = STT(I,J,L,IDTOCIL1-1+BIN)
         MOCOB = STT(I,J,L,IDTOCOB1-1+BIN)
      ENDIF
       IF ( IDTDUST1 > 0 ) MDUST = STT(I,J,L,IDTDUST1-1+BIN)
c       IF (LDUST30) MDUST = STT(I,J,L,IDTDUST1-1+BIN)
      MSO4  = STT(I,J,L,IDTSF1-1+BIN) * 1.2 !account for ammonium sulfate
      MNACL = STT(I,J,L,IDTSS1-1+BIN)
      
      MTOT  = MECIL + MOCIL + MOCOB + MSO4 + MNACL + MDUST + 1.e-20
      XOCIL = MOCIL / MTOT
      XSO4  = MSO4 / MTOT
      XNACL = MNACL / MTOT
      ISO4  = MIN(101, INT(XSO4*100)+1)
      INACL = MIN(101, INT(XNACL*100)+1)
      IOCIL = MIN(101, INT(XOCIL*100)+1)

      !==========================================================
      ! subroutine was written considering bin 1 is 10nm
      ! in TOMAS-40, bin 1 is 1nm and bin 11 is 10nm
      !==========================================================
#if   defined( TOMAS40 )
      OFFSET = 10
#elif  defined( TOMAS15 )
      OFFSET = 3
#else
      OFFSET = 0
#endif

      IF ( LS ) THEN 
         BINACT = BINACT1(ISO4, INACL, IOCIL) + OFFSET
      ELSE 
         BINACT = BINACT2(ISO4, INACL, IOCIL) + OFFSET
      ENDIF

      ! Free pointer
      NULLIFY( STT )

      RETURN
      END SUBROUTINE GETACTBIN

!------------------------------------------------------------------------------

      SUBROUTINE EZWATEREQM( Mke, RHTOMAS )
!
!******************************************************************************
!     WRITTEN BY Peter Adams, March 2000
!
!     This routine uses the current RH to calculate how much water is 
!     in equilibrium with the aerosol.  Aerosol water concentrations 
!     are assumed to be in equilibrium at all times and the array of 
!     concentrations is updated accordingly.
!
!     Introduced to GEOS-CHEM by Win Trivitayanurak. May 8, 2006.
!     This file is replacing the old ezwatereqm.f that was not compatible
!     with multicomponent aerosols.  The new ezwatereqm use external 
!     functions to do ISORROPIA-result curve fitting for each aerosol 
!     component.
!     WARNING :
!      *** Watch out for the new aerosol species added in the future!
!
!     VARIABLE COMMENTS...
!
!     This version of the routine works for sulfate and sea salt
!     particles.  They are assumed to be externally mixed and their
!     associated water is added up to get total aerosol water.
!     wr is the ratio of wet mass to dry mass of a particle.  Instead
!     of calling a thermodynamic equilibrium code, this routine uses a
!     simple curve fits to estimate wr based on the current humidity.
!     The curve fit is based on ISORROPIA results for ammonium bisulfate
!     at 273 K and sea salt at 273 K.
!
!  NOTE:
!  ( 1) Add OC mass in the calculation (win, 9/3/07)
!******************************************************************************
!
      ! Arguments
      REAL*8,INTENT(INOUT) :: Mke(IBINS,ICOMP)
      REAL*4,INTENT(IN)    :: RHTOMAS

      ! Local variables
      INTEGER           :: k, j
      REAL*8            :: so4mass, naclmass, ocilmass
      REAL*8            :: wrso4, wrnacl, wrocil
      REAL*8            :: rhe

      !========================================================================
      ! EZWATEREQM begins here!
      !========================================================================
      
      rhe=100.d0*rhtomas
      if (rhe .gt. 99.d0) rhe=99.d0
      if (rhe .lt. 1.d0) rhe=1.d0

      do k=1,ibins

         so4mass=Mke(k,srtso4)*1.2  !1.2 converts kg so4 to kg nh4hso4
         wrso4=waterso4(rhe)

         ! Add condition for srtnacl in case of running so4 only. (win, 5/8/06)
         if (srtnacl.gt.0) then 
            naclmass=Mke(k,srtnacl) !already as kg nacl - no conv necessary
            wrnacl=waternacl(rhe)
         else
            naclmass = 0.d0
            wrnacl = 1.d0
         endif

         if (srtocil.gt.0) then 
            ocilmass=Mke(k,srtocil) !already as kg ocil - no conv necessary
            wrocil=waterocil(rhe)
         else
            ocilmass = 0.d0
            wrocil = 1.d0
         endif

         Mke(k,srth2o)=so4mass*(wrso4-1.d0)+naclmass*(wrnacl-1.d0)
     &                 +ocilmass*(wrocil-1.d0)

      enddo

      RETURN
      END SUBROUTINE EZWATEREQM

!------------------------------------------------------------------------------

      SUBROUTINE EZWATEREQM2( I, J, L, BIN, State_Met, State_Chm ) 
!
!******************************************************************************
!  Subroutine EZWATEREQM2 is just like EZWATEREQM but access directly to STT
!  array unlike EZWATEREQM that needs the array Mke to be passed in and out.  
!  This subroutine is for calling from outside microphysics module.
!  (win, 8/6/07)
!
!  NOTES:
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!******************************************************************************
!
      ! Reference to F90 modules
      USE CMN_SIZE_MOD             ! Size parameters
      USE GIGC_State_Chm_Mod, ONLY : ChmState    
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TRACERID_MOD,       ONLY : IDTNK1, IDTSF1, IDTSS1
      USE TRACERID_MOD,       ONLY : IDTAW1, IDTOCIL1

      ! Arguments
      INTEGER,  INTENT(IN)          :: I, J, L, BIN
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object

      ! Local variables
      REAL*8               :: RHE
      REAL*8               :: SO4MASS, NACLMASS, OCILMASS
      REAL*8               :: WRSO4, WRNACL, WROCIL

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER   :: STT(:,:,:,:)

      !========================================================================
      ! EZWATEREQM begins here!
      !========================================================================

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      rhe = State_Met%RH(i,j,l)             !RH [=] percent
      
      if (rhe .gt. 99.) rhe=99.
      if (rhe .lt. 1.) rhe=1.
      
      so4mass=STT(i,j,l,IDTSF1-1+bin)*1.2 !1.2 converts kg so4 to kg nh4hso4
      wrso4=waterso4(rhe)       !use external function
      
      ! Add condition for srtnacl in case of running so4 only. (win, 5/8/06)
      if (IDTSS1.gt.0) then 
         naclmass=STT(i,j,l,IDTSS1-1+bin) !already as kg nacl - no conv necessary
         wrnacl=waternacl(rhe)  !use external function
      else
         naclmass = 0.d0
         wrnacl = 1.d0
      endif
      
      if (IDTOCIL1 > 0) then 
         ocilmass=STT(i,j,l,IDTOCIL1-1+bin)  !already as kg ocil - no conv necessary
         wrocil=waterocil(rhe)
      else
         ocilmass = 0.d0
         wrocil = 1.d0
      endif

      STT(i,j,l,IDTAW1-1+bin)= so4mass*(wrso4-1.d0) +
     &                            naclmass*(wrnacl-1.d0)
     &                 + ocilmass*(wrocil-1.d0)

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE EZWATEREQM2   
     
!------------------------------------------------------------------------------

      SUBROUTINE EZNH3EQM( Gce, Mke ) 
!
!******************************************************************************
!  Subroutine EZNH3REQM2 puts ammonia to the particle phase until 
!     there is 2 moles of ammonium per mole of sulfate and the remainder
!     of ammonia is left in the gas phase.
!     (win, 9/30/08)
!******************************************************************************
!
      IMPLICIT NONE

      ! Arguments
!      REAL*8,  INTENT(INOUT)  :: Gce(icomp)
      REAL*8,  INTENT(INOUT)  :: Gce(icomp - 1) !sfarina - fixed incorrect definition of Gc array
      REAL*8,  INTENT(INOUT)  :: Mke(ibins,icomp)

      ! Local variables
      integer       ::  k
      REAL*8        :: tot_nh3  !total kmoles of ammonia
      REAL*8        :: tot_so4  !total kmoles of so4
      REAL*8        :: sfrac    !fraction of sulfate that is in that bin

      !========================================================================
      ! EZNH3EQM begins here!
      !========================================================================

      ! get the total number of kmol nh3
      tot_nh3 = Gce(srtnh4)/17.d0
      do k=1,ibins
         tot_nh3 = tot_nh3 + Mke(k,srtnh4)/18.d0
      enddo

      ! get the total number of kmol so4
      tot_so4 = 0.d0
      do k=1,ibins
         tot_so4 = tot_so4 + Mke(k,srtso4)/96.d0
      enddo

      ! see if there is free ammonia
      if (tot_nh3/2.d0.lt.tot_so4)then  ! no free ammonia
         Gce(srtnh4) = 0.d0 ! no gas phase ammonia
         do k=1,ibins
            sfrac = Mke(k,srtso4)/96.d0/tot_so4
            Mke(k,srtnh4) = sfrac*tot_nh3*18.d0 ! put the ammonia where the sulfate is
#if defined ( DEBUG )
            if ( Mke(k,srtnh4) < 0.0 ) then
               print *,'negative gas phase ammonia in eznh3eqm!!'
               print *,'bin  ', k
            endif
#endif
         enddo
      else ! free ammonia
         do k=1,ibins
            Mke(k,srtnh4) = Mke(k,srtso4)/96.d0*2.d0*18.d0 ! fill the particle phase
#if defined ( DEBUG )
            if ( Mke(k,srtnh4) < 0.0 ) then
               print *,'negative gas phase ammonia in eznh3eqm!!'
               print *,'bin  ', k
            endif
#endif
         enddo
         Gce(srtnh4) = (tot_nh3 - tot_so4*2.d0)*17.d0 ! put whats left over in the gas phase
#if defined ( DEBUG )
         if ( Gce(srtnh4) < 0.0 ) then
            print *,'negative gas phase ammonia in eznh3eqm!!'
         endif
#endif

      endif

      RETURN
      END SUBROUTINE EZNH3EQM      
      
!------------------------------------------------------------------------------

      SUBROUTINE AERO_DIADEN( DIA,       DENSITY,   LEV, 
     &                        Input_Opt, State_Met, State_Chm )

!******************************************************************************
!  Subroutine AERO_DIADEN calculate the diameter and density by calling external
!  functions GETDP and AERODENS respectively. (win, 7/19/07)
!  Note: This subroutine is created for supplying diameter and density for 
!        dry dep velocity calculation in DEPVEL.  Did not want to add much 
!        to DEPVEL.
!******************************************************************************

      USE CMN_SIZE_MOD             ! Size parameters
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TRACERID_MOD,       ONLY : IDTNK1,   IDTSF1,   IDTSS1
      USE TRACERID_MOD,       ONLY : IDTDUST1, IDTECIL1, IDTECOB1
      USE TRACERID_MOD,       ONLY : IDTOCIL1, IDTOCOB1

      IMPLICIT NONE

      ! Arguments 
      REAL*8,         INTENT(OUT)   :: DIA(MAXIJ, IBINS)
      REAL*8,         INTENT(OUT)   :: DENSITY(MAXIJ, IBINS)
      INTEGER,        INTENT(IN)    :: LEV  
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      
      ! Local variables
      INTEGER           :: I,J, BIN, JC, IJLOOP, TRACID, WID
!      REAL*8, EXTERNAL  :: AERODENS 
      REAL*8            :: AMASS(ICOMP)
      REAL*8            :: MSO4, MNACL, MH2O
      REAL*8            :: MECIL, MECOB, MOCIL, MOCOB, MDUST

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER   :: STT(:,:,:,:)

      !========================================================================
      ! AERO_DIADEN begins here!
      !========================================================================
      
      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Initialize mass
      MSO4  = 0d0
      MNACL = 0d0
      MH2O  = 0d0
      MECIL = 0d0
      MECOB = 0d0
      MOCIL = 0d0
      MOCOB = 0d0
      MDUST = 0d0

      CALL CHECKMN( 0,         0,         0, 
     &              Input_Opt, State_Met, State_Chm, 
     &              'kg',      'AERO_DIADEN called from DEPVEL' )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( MECIL, MECOB, MOCIL, MOCOB, MDUST )
!$OMP+PRIVATE( BIN, IJLOOP, I, J, TRACID, WID, MH2O, MSO4, MNACL )
!$OMP+SCHEDULE( DYNAMIC )
      DO IJLOOP = 1, MAXIJ
      DO BIN = 1, IBINS

         ! Calculate the I and J-index from IJLOOP index
         I = MOD( IJLOOP, IIPAR )
         IF ( I==0 ) I = 72
         J = ( ( IJLOOP - I ) / IIPAR ) + 1

         TRACID = IDTNK1 + BIN - 1
c         print *,"TRACID=",TRACID,"IDTNK1=",IDTNK1, "BIN=", BIN
         WID    = IDTNK1 + (ICOMP - 1)*IBINS - 1 + BIN  !(fixed WID to 281-310. dmw 10/3/09)
c         print *, "wid=", WID, "ICOMP=", ICOMP, "IBINS=", IBINS
         ! Get the diameter from an external function
         DIA(IJLOOP,BIN) = GETDP( I,      J,         LEV, 
     &                            TRACID, State_Met, State_Chm )
         
         ! Prepare the mass to call external function for density
         MH2O = STT(I,J,1,WID         )
c         IF ( LSULF30 ) MSO4  = STT(I,J,LEV,IDTSF1-1+BIN)
         IF ( IDTSF1 > 0 ) MSO4  = STT(I,J,LEV,IDTSF1-1+BIN)
C         IF ( LSALT30 ) MNACL = STT(I,J,LEV,IDTSS1-1+BIN)
         IF ( IDTSS1 > 0 ) MNACL = STT(I,J,LEV,IDTSS1-1+BIN)
         IF ( IDTECIL1 > 0 .AND.IDTECOB1 > 0 .AND. 
     &        IDTOCIL1 > 0 .AND. IDTOCOB1 > 0 ) THEN
C         IF ( LCARB30 ) THEN 
            MECIL = STT(I,J,LEV,IDTECIL1-1+BIN)
            MECOB = STT(I,J,LEV,IDTECOB1-1+BIN)
            MOCIL = STT(I,J,LEV,IDTOCIL1-1+BIN)
            MOCOB = STT(I,J,LEV,IDTOCOB1-1+BIN)
         ENDIF
C         IF ( LDUST30 ) MDUST = STT(I,J,LEV,IDTDUST1-1+BIN)
        IF ( IDTDUST1 > 0 ) MDUST = STT(I,J,LEV,IDTDUST1-1+BIN)

        ! Get density from external function
         DENSITY(IJLOOP,BIN) = AERODENS( MSO4, 0.d0, 1.875D-1*MSO4, 
     &                                   MNACL, MECIL, MECOB, 
     &                                   MOCIL, MOCOB, MDUST, MH2O  )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      RETURN
      END SUBROUTINE AERO_DIADEN
      
!-----------------------------------------------------------------------------
      
      SUBROUTINE CHECKMN
     &  ( II, JJ, LL, Input_Opt, State_Met, State_Chm, UNITS, LOCATION )
!      
!*****************************************************************************
!  Subroutine CHECKMN use the subroutine MNFIX to check for error in the 
!  aerosol mass and number inconsistencies. (win, 7/24/07)
!  Variables
!  I, J, L    (INTEGER )  : Grid box
!  UNITS       (CHARACTER) : Current unit of STT tracer array
!  LOCATION   (CHARACTER) : Location of calling program / subroutine
!
!  NOTES:
!  (1 ) Now just access STT directly and avoid having to pass values in-out
!       (Win, 7/24/07)
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!*****************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,            ONLY : CONVERT_UNITS      
      USE ERROR_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TRACERID_MOD,       ONLY : IDTNK1,  IDTAW1, IDTNH4, IDTNH3

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND60

      ! Arguments
      INTEGER                         :: II, JJ, LL
      CHARACTER(LEN=*), INTENT(IN)    :: UNITS, LOCATION
      TYPE(OptInput),   INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState),   INTENT(IN)    :: State_Met   ! Meteorology State object
      TYPE(ChmState),   INTENT(INOUT) :: State_Chm   ! Chemistry State object

      ! Local variables
      INTEGER                :: I, J, L
      INTEGER                :: I1, I2, J1, J2, L1, L2
      INTEGER                :: K, JC, NKID, TRACNUM, MPNUM
      LOGICAL                :: ERRORSWITCH
                               ! Make ERRORSWITCH = .TRUE. to get full print 
                               ! for debugging

!sfarina
      REAL*8                         :: Nk(IBINS), Nkd(IBINS)
      REAL*8                         :: Mk(IBINS, ICOMP)
      REAL*8                         :: Mkd(IBINS,ICOMP)
      REAL*8                         :: Gc(ICOMP - 1)
      REAL*8                         :: Gcd(ICOMP - 1)
      REAL*4                         :: BOXVOL
      REAL*8                         :: XFER(IBINS)
    
      ! For values from Input_Opt
      LOGICAL :: LPRT
      INTEGER :: N_TRACERS
      REAL*8  :: XNUMOL(Input_Opt%N_TRACERS)
      REAL*8  :: TCVV(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER   :: STT(:,:,:,:)

      !=================================================================
      ! CHECKMN begins here!
      !=================================================================

      ! Copy values from Input_Opt
      LPRT      = Input_Opt%LPRT
      N_TRACERS = Input_Opt%N_TRACERS
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)
      TCVV      = Input_Opt%TCVV(1:N_TRACERS)

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT       => State_Chm%Tracers

      ERRORSWITCH = .FALSE.

      ! We want to work with STT in [kg]
      IF ( TRIM( UNITS ) == 'V/V' .or. TRIM( UNITS ) == 'v/v' )
      ! Convert STT from v/v --> kg
     &     CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, State_Met%AD, STT ) 

      ! Check throughout all grid boxes
      IF ( II == 0 .and. JJ == 0 .and. LL == 0 ) THEN 
         I1 = 1
         I2 = IIPAR
         J1 = 1
         J2 = JJPAR
         L1 = 1
         L2 = LLPAR
      ELSE ! Check at a single grid
         I1 = II
         I2 = II
         J1 = JJ
         J2 = JJ
         L1 = LL
         L2 = LL
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
!$OMP+PRIVATE( Nk, Nkd, Mk, Mkd, K, TRACNUM, JC, MPNUM, BOXVOL )
!$OMP+PRIVATE( GC, GCd, ERRORSWITCH )
!$OMP+SCHEDULE( DYNAMIC )

!%%%%+DEFAULT( PRIVATE )
!%%%%+SHARED( STT, State_Met, IDTNK1, ICOMP, IDIAG, IDTAW1 )
      DO I = I1, I2
      DO J = J1, J2
      DO L = L1, L2
         
         BOXVOL  = State_Met%AIRVOL(I,J,L) * 1.e6 !convert from m3 -> cm3 

         ! Swap GEOSCHEM variables into aerosol algorithm variables
         DO K = 1, IBINS
            TRACNUM = IDTNK1 - 1 + K
            ! Check for nan
            IF ( IT_IS_NAN( STT(I,J,L,TRACNUM) ) ) 
     &           print *, 'Found NaN at',I, J, L,'Tracer',TRACNUM
            NK(K) = STT(I,J,L,TRACNUM)
            DO JC = 1, ICOMP-IDIAG
               TRACNUM = IDTNK1 - 1 + K + IBINS*JC
               IF ( IT_IS_NAN( STT(I,J,L,TRACNUM) ) ) 
     &              print *, 'Found NaN at',I, J, L,'Tracer',TRACNUM
               MK(K,JC) = STT(I,J,L,TRACNUM)
            ENDDO
            MK(K,SRTH2O) = STT(I,J,L,IDTAW1-1+K)
         ENDDO

         DO JC = 1, ICOMP - 1
            Gc(JC) = 0.0   !sfarina - init Gc so storenm doesn't go NaN on us.
         END DO

         ! Get NH4 mass from the bulk mass and scale to bin with sulfate
         IF ( SRTNH4 > 0 ) THEN
!-----------------------------------------------------------------------------
! Prior to 1/29/14:
! Rewrite subroutine call to avoid array temporaries (bmy, 1/29/14)
!            CALL NH4BULKTOBIN( MK(1:IBINS,SRTSO4), STT(I,J,L,IDTNH4),
!     &                         XFER(1:IBINS)  )
!-----------------------------------------------------------------------------
            CALL NH4BULKTOBIN( MK(:,SRTSO4), STT(I,J,L,IDTNH4), XFER  )
            MK(1:IBINS,SRTNH4) = XFER(1:IBINS)
            Gc(SRTNH4) = STT(I,J,L,IDTNH3) 
         ENDIF

!         if ( i==26 .and. j==57 .and. l==13 ) 
!     &        call debugprint(Nk,Mk,i,j,l,'In CHECKMN')

         CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd) 
         !if(i==47.and.j==10.and.l==7) ERRORSWITCH = .TRUE.
C         if(i==22.and.j==33.and.l==2) 
C     &      call debugprint( NK, MK, i,j,l,LOCATION)
!%%%% CRITICAL
         CALL MNFIX( NK, MK, ERRORSWITCH )
         IF ( ERRORSWITCH ) THEN
            PRINT *, 'CHECKMN is going to terminate at grid',I,J,L
!            IF( .not. SPINUP(14.0) ) THEN
               CALL ERROR_STOP( 'MNFIX found error', LOCATION )
!            ELSE
!               PRINT *,'Let error go during spin up'
!            ENDIF
         ENDIF
         
         ! Save the error fixing to diagnostic AERO-FIX
         MPNUM = 5
         IF ( ND60 > 0 )
     &      CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL )
         

         ! Swap Nk and Mk arrays back to STT 
         DO K = 1, IBINS
            TRACNUM = IDTNK1 - 1 + K
            STT(I,J,L,TRACNUM) = Nk(K)
            DO JC = 1, ICOMP-IDIAG
               TRACNUM = IDTNK1 - 1 + K + IBINS*JC
               STT(I,J,L,TRACNUM) = Mk(K,JC)
            ENDDO
            STT(I,J,L,IDTAW1-1+K) = MK(K,SRTH2O)
         ENDDO
      
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      ! In the input unit is [v/v], then convert [kg] back to [v/v]
      IF ( TRIM( UNITS ) == 'V/V' .or. TRIM( UNITS ) == 'v/v' ) 
       ! Convert STT from kg --> v/v
     &     CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, State_Met%AD, STT ) 
        
      IF ( LPRT ) WRITE(6,*)' #### CHECKMN: finish at ',LOCATION

      ! Free pointer
      NULLIFY( STT )

      ! Return to calling program
      END SUBROUTINE CHECKMN    

!-----------------------------------------------------------------------------
      
      SUBROUTINE MNFIX ( NK, MK, ERRORSWITCH )
!
!*****************************************************************************
!  Subroutine MNFIX examines the mass and number distrubution and determine if
!  any bins have an average mass outside their normal range.  This can happen 
!  because some process, e.g. advection, seems to treat the mass and number
!  tracers inconsistently.  If any bins are out of range, I shift some mass 
!  and number to a new bin in a way that conserves both. 
!  (win, 7/23/07)
!  ORiginally written by Peter Adams, September 2000
!  Modified for GEOS-CHEM by Win Trivitayanurak (win@cmu.edu)
!
!  NOTES: 
!  (1 ) Getting rid of eps which was added to the denominator to 
!       avoid a division by zero when calculating avg mass/number
!       and calculating multicomponent mass species fj=MK(k,j)/drymass.
!  (2 ) Introduce this code from a previous version of GEOSCHEM -- so this 
!       may need cleaning up later since this version is full of debugging
!       print statements (Win, 7/23/07)
!  (3 ) Bug fix - now when shifting mass or number, it can't try to go beyond
!       the max or min bin. (win, 8/1/07)
!*****************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP, IT_IS_NAN

      ! Arguments
      REAL*8,  INTENT(INOUT) :: NK(IBINS),  MK(IBINS, ICOMP)
      LOGICAL, INTENT(INOUT) :: ERRORSWITCH

      ! Local variables
      integer             :: K,J,KK !counters
      integer             :: NEWBIN !bin number into which mass is shifted
      REAL*8              :: XOLD, XNEW !average masses of old and new bins
      REAL*8              :: DRYMASS !dry mass of in a bin
      REAL*8              :: AVG !average dry mass of particles in bin
      REAL*8              :: NUM_INITIAL !number of particles initially in problem bin
      REAL*8              :: NSHIFT  !number to shift to new bin
      REAL*8              :: MSHIFT !mass to shift to new bin
      REAL*8              :: FJ !fraction of mass that is component j
      REAL*8              :: save1,save2,save3,save4,save5
!      REAL*8,   PARAMETER :: EPS2= 1.D-32 !small number for Mk 
!      REAL*8,   PARAMETER :: EPS = 1.D-2 !small number for Nk 
      REAL*8,   PARAMETER :: EPS  = 1.D-20 !small number for Nk
      REAL*8,   PARAMETER :: EPS2 = 1.D-28 !small number for Mk
      REAL*8,   PARAMETER :: TINY = 1.D-36 !small number
      REAL*8,   PARAMETER :: VTINY= 1.D-50 !very small number

      LOGICAL             :: FIXERROR
      LOGICAL             :: PRT
      REAL*8              :: TOTMAS, TOTNUM !for print debug 

      !=================================================================
      ! MNFIX begins here!
      !=================================================================
!%%%%+PRIVATE(K, J, KK, NEWBIN, XOLD, XNEW, DRYMASS, AVG)
!%%%%+PRIVATE(NUM_INITIAL, NSHIFT, MSHIFT, FJ )
!%%%%+PRIVATE( save1, save2, save3, save4, save5)
!%%%%+PRIVATE(FIXERROR, PRT, TOTMAS, TOTNUM)
      FIXERROR = .TRUE.
      PRT = .FALSE.
      PRT = ERRORSWITCH !just carrying a signal to print out value at the observed box - since mnfix does not have any information about I,J,L location. (Win, 9/27/05)
      ERRORSWITCH = .FALSE.

!      xk(1)=xk(2)/2.d0  ! jrp for some reason xk(1) is changing?!
      save1=xk(1)

      ! Check for any incoming negative values or NaN
      !--------------------------------------------------------------------------
      DO K = 1, IBINS
        IF ( IT_IS_NAN(NK(K)) ) THEN
           PRINT *,'11 Found Nan in Nk at bin',K
           ERRORSWITCH = .TRUE.
           print *,'11 MNFIX(0): Found NaN in Nk(,',k,')'
           GOTO 300
        ENDIF
        DO J = 1, ICOMP
           IF ( IT_IS_NAN(MK(K,J)) ) THEN
              PRINT *,'11 Found Nan in Mk at bin',K,'component',J
              ERRORSWITCH = .TRUE.
              GOTO 300
           ENDIF
        ENDDO
        IF ( NK(K) < 0D0 ) THEN
           IF ( PRT ) THEN 
              PRINT *,'MNFIX[0]: FOUND NEGATIVE N'
              PRINT *, 'Bin, N', K, NK(K)
           ENDIF
           IF ( ABS(NK(K)) < 1D0 .and. FIXERROR ) THEN
              NK(K) = 0D0
              IF ( PRT ) PRINT *,'Negative N > -1.0 Reset to zero'
           ELSE
              ERRORSWITCH = .TRUE.
              print *,'MNFIX(0): Found negative Nk(',k,') <-1d0'
              GOTO 300          !exit mnfix if found negative error (win, 4/18/06)
           ENDIF
        ENDIF
        IF ( IT_IS_NAN(NK(K)) ) THEN 
           PRINT *,'Found Nan in Nk at bin',K
           ERRORSWITCH = .TRUE.
           print *,'MNFIX(0): Found NaN in Nk(,',k,')'
           GOTO 300
        ENDIF
        DO J = 1, ICOMP
           IF ( MK(K,J) < 0D0 ) THEN
              IF ( PRT ) THEN 
                 PRINT *,'MNFIX[0]: FOUND NEGATIVE M'
                 PRINT *,'Bin, Comp, Mk', K, J, MK(K,J)
              ENDIF
              IF( ABS(MK(K,J)) < 1D-5 .and. FIXERROR ) THEN
                 MK(K,J) = 0D0
                 IF ( PRT ) PRINT *,'Negative M > -1.d-5 Reset to zero'
              ELSE
                 ERRORSWITCH =.TRUE.
                 print *,'MNFIX(0): Found negative Mk(',k,',comp',j,')'
                 GOTO 300       !exit mnfix if found negative error (win, 4/18/06)
              ENDIF
           ENDIF
           IF ( IT_IS_NAN(MK(K,J)) ) THEN 
              PRINT *,'Found Nan in Mk at bin',K,'component',J
              ERRORSWITCH = .TRUE.
              GOTO 300
           ENDIF
        ENDDO                   !icomp
      ENDDO                     !ibins
      save2=xk(1)

      ! JRP check for neg numbers
!      DO K = 1,IBINS
!        IF (NK(K) < 0.d0) THEN
!           print*,'1 NK < 0 in MNFIX',K,NK(K)
!        ENDIF
!        DO J=1,ICOMP
!           IF (MK(K,J) < 0.d0) THEN
!              print*,'1 MK < 0 in MNFIX',K,J,MK(K,J)
!           ENDIF
!        ENDDO
!        IF ( IT_IS_NAN(NK(K)) ) THEN
!           PRINT *,'11 Found Nan in Nk at bin',K
!           ERRORSWITCH = .TRUE.
!           print *,'11 MNFIX(0): Found NaN in Nk(,',k,')'
!           GOTO 300
!        ENDIF
!        DO J = 1, ICOMP
!           IF ( IT_IS_NAN(MK(K,J)) ) THEN
!              PRINT *,'11 Found Nan in Mk at bin',K,'component',J
!              ERRORSWITCH = .TRUE.
!              GOTO 300
!           ENDIF
!        ENDDO
!      ENDDO

      ! Check if both number and mass are zero, if yes then exit mnfix.
      !----------------------------------------------------------------
      TOTNUM = 0d0 
      TOTMAS = 0D0
      DO K = 1,IBINS
         TOTNUM = TOTNUM + NK(K)
         DO J=1,ICOMP-IDIAG
            TOTMAS = TOTMAS + MK(K,J)
         ENDDO
      ENDDO
      IF ( TOTNUM == 0D0 .AND. TOTMAS == 0D0 ) THEN
          IF ( PRT ) PRINT *,'MNFIX: Nk=Mk=0. Exit now'
         GOTO 300
      ENDIF

      ! If number is tiny ( < EPS) then set it to zero
!      DO K = 1,IBINS
!         IF ( NK(K) <= EPS ) THEN
!            NK(K) = 0D0
!            DO J= 1, ICOMP-1
!               MK(K,J) = 0d0
!            ENDDO               !STOP  !original (win, 9/1/05)

!         ENDIF
!      ENDDO

      ! If N is tiny and M is tiny, set both to zeroes
      !--------------------------------------------------------
      DO K = 1, IBINS
        IF ( IT_IS_NAN(NK(K)) ) THEN
           PRINT *,'22 Found Nan in Nk at bin',K
           ERRORSWITCH = .TRUE.
           print *,'22 MNFIX(0): Found NaN in Nk(,',k,')'
           GOTO 300
        ENDIF
        DO J = 1, ICOMP
           IF ( IT_IS_NAN(MK(K,J)) ) THEN
              PRINT *,'22 Found Nan in Mk at bin',K,'component',J
              ERRORSWITCH = .TRUE.
              GOTO 300
           ENDIF
        ENDDO
       IF ( NK(K) <= EPS .AND. NK(K)>= 0D0 ) THEN
         !print*,'1111'
         !print*,k,EPS,xk(K),xk(K+1)
         !print*,'word up'
         NK(K) = EPS
         !NK(K) = 0.d0
         !DO J = 1, ICOMP-IDIAG
         DO J = 1, ICOMP
            if (J .eq. 1) then
               !MK(K,J) = EPS*sqrt(xk(K)*xk(K+1))
               MK(K,J) = EPS*AVGMASS(k)
               !MK(K,J) = 0.d0
            else
               MK(K,J) = VTINY
            endif
         enddo
         !print*,'allbins',MK(:,1)
       ENDIF ! If tiny number
         TOTMAS = SUM(MK(K,1:ICOMP-IDIAG))
         if (TOTMAS.lt.eps2) then
            !print*,'2222'
            NK(K) = EPS
            !NK(K) = 0.d0
            DO J = 1, ICOMP
            !DO J = 1, ICOMP-IDIAG
               if (J .eq. 1) then
                  !MK(K,J) = EPS*sqrt(xk(K)*xk(K+1))
                  MK(K,J) = EPS*AVGMASS(k)
                  !MK(K,J) = 0.d0
               else
                  MK(K,J) = VTINY
               endif
            enddo
         endif
      ENDDO
      save3=xk(1)

      ! JRP check for neg numbers
      DO K = 1,IBINS
        IF (NK(K) < 0.d0) THEN
           print*,'2 NK < 0 in MNFIX',K,NK(K)
        ENDIF
        DO J=1,ICOMP
           IF (MK(K,J) < 0.d0) THEN
              print*,'2 MK < 0 in MNFIX',K,J,MK(K,J)
           ENDIF
        ENDDO
        IF ( IT_IS_NAN(NK(K)) ) THEN
           PRINT *,'2 Found Nan in Nk at bin',K
           ERRORSWITCH = .TRUE.
           print *,'2 MNFIX(0): Found NaN in Nk(,',k,')'
           GOTO 300
        ENDIF
        DO J = 1, ICOMP
           IF ( IT_IS_NAN(MK(K,J)) ) THEN
              PRINT *,'2 Found Nan in Mk at bin',K,'component',J
              ERRORSWITCH = .TRUE.
              GOTO 300
           ENDIF
        ENDDO
      ENDDO

      ! Check to see if any bins are completely out of bounds for min or max bin
      !-------------------------------------------------------------------------
      DO K = 1, IBINS
         DRYMASS = 0.d0
         DO J = 1, ICOMP-IDIAG
            DRYMASS = DRYMASS + MK(K,J)
         ENDDO

         IF ( NK(k) == 0d0 ) THEN
            !AVG = SQRT( xk(K)* xk(K+1) )
            AVG = SQRT( AVGMASS(k) )
         ELSE
            AVG = DRYMASS/ NK(K)
         ENDIF

         IF ( AVG >  xk(IBINS+1) ) THEN
            IF ( PRT ) PRINT *, 'MNFIX [1]: AVG > Xk(ibins+1) at bin',K
            IF ( FIXERROR ) THEN 
               !out of bin range - remove some mass
               MSHIFT = NK(k)* xk(IBINS+1)/ 1.2
               DO J= 1, ICOMP
                  MK(K,J) = MK(K,J)* MSHIFT/ (DRYMASS+EPS2)
               ENDDO
            ELSE
               ERRORSWITCH = .TRUE.
               print *,'MNFIX(1): AVG>Xk(ibins+1) at bin',K
               GOTO 300
            ENDIF
         ENDIF
         IF ( AVG < xk(1)) THEN
            IF( PRT ) PRINT *,'MNFIX [2]: AVG < Xk(1)'
            IF( FIXERROR ) THEN 
               !out of bin range - remove some number
               NK(K) = DRYMASS/ ( xk(1)* 1.2 )
            ELSE
               ERRORSWITCH = .TRUE.
               print *,'MNFIX(1): AVG < Xk(1) at bin',K
               GOTO 300
            ENDIF
         ENDIF
      ENDDO

      ! JRP check for neg numbers
      DO K = 1,IBINS
        IF (NK(K) < 0.d0) THEN
           print*,'3 NK < 0 in MNFIX',K,NK(K)
        ENDIF
        DO J=1,ICOMP
           IF (MK(K,J) < 0.d0) THEN
              print*,'3 MK < 0 in MNFIX',K,J,MK(K,J)
           ENDIF
        ENDDO
        IF ( IT_IS_NAN(NK(K)) ) THEN
           PRINT *,'3 Found Nan in Nk at bin',K
           ERRORSWITCH = .TRUE.
           print *,'3 MNFIX(0): Found NaN in Nk(,',k,')'
           GOTO 300
        ENDIF
        DO J = 1, ICOMP
           IF ( IT_IS_NAN(MK(K,J)) ) THEN
              PRINT *,'3 Found Nan in Mk at bin',K,'component',J
              ERRORSWITCH = .TRUE.
              GOTO 300
           ENDIF
        ENDDO
      ENDDO
      save4=xk(1)

!      if (PRT) then !<step5.1-temp>
!         print *,'After_Check2 ---------------------'
!         do k=1,ibins
!            totmas = sum(MK(k,1:icomp-1))
!            print *, totmas,NK(k), totmas/NK(k)
!         enddo
!      endif

      !print*,1,NK(1),NK(2)
      !print*,1,MK(1,:)
      !print*,1,MK(2,:)

      ! Check to see if any bins are out of bounds
      !-------------------------------------------------------------------
      DO K = 1, IBINS
!         if (PRT) print *,'Now at bin',k !<step4.4>tmp (win, 9/28/05)

         DRYMASS = 0.d0
         DO J = 1, ICOMP-IDIAG
            DRYMASS = DRYMASS + MK(K,J)
         ENDDO

         IF ( NK(K) == 0d0 ) THEN
            !AVG = SQRT(xk(K)*xk(K+1)) !set to mid-range value
            AVG = AVGMASS(k) !set to mid-range value
         ELSE
            AVG = DRYMASS/NK(K)
         ENDIF

!         if (PRT) then     !<step5.1-temp>
!            print *,'After_Check3---------------------'
!            totmas = sum(MK(k,1:icomp-1))
!            print *, totmas,NK(k), totmas/NK(k)
!         endif

         ! If over boundary of the current bin
         IF ( AVG >  xk(K+1) ) THEN
           IF ( PRT ) PRINT *, 'MNFIX [3]: AVG>Xk(',K+1,')'
           IF ( PRT ) CALL DEBUGPRINT(NK,MK,0,0,0,'inside MNFIX')
           IF ( FIXERROR ) THEN 
            !Average mass is too high - shift to higher bin
            !KK = K + 1 ! jrp, this was causing errors
!      ERRORSWITCH=.TRUE.
            KK = K
            XNEW = xk(KK+1)/ 1.1
            if ( PRT ) PRINT *, 'k',k,'AVG',AVG,' XNEW ',XNEW
 100        IF ( XNEW <= AVG ) THEN
               IF ( KK < IBINS ) THEN
                  KK = KK + 1
                  XNEW = xk(KK+1)/ 1.1
                  if (PRT) PRINT *, '..move up to bin ',KK,' XNEW ',XNEW
                  GOTO 100
               ELSE
                  ! Already reach highest bin - must remove some mass (win, 8/1/07)
                  ! Updated by jrp 3/1/2012
                  MSHIFT = NK(k)* xk(k+1)/ 1.1
                  if( PRT ) PRINT*,' Mass being discarded: '
                  DO J= 1, ICOMP
                     !if (PRT) 
!                     print*,'Removing mass in MNFIX',
!     &                    MSHIFT, DRYMASS
                     MK(K,J) = MK(K,J)* MSHIFT/ (DRYMASS)
                  ENDDO 
                  ! and recalculate dry mass (win, 8/1/07)
                  DRYMASS = 0.d0 ! jrp fix 2/29/12
                  DO J = 1, ICOMP-IDIAG
                     DRYMASS = DRYMASS + MK(K,J)
                  ENDDO
                  GOTO 111
               ENDIF         
            ENDIF

            if(PRT)print*,'Old NK',NK(k),'Old DRYMASS',DRYMASS,'bin',k

            !XOLD = SQRT( xk(K)* xk(K+1) )
            XOLD = AVGMASS(k)
            NUM_INITIAL = NK(K)
            NSHIFT = ( DRYMASS - XOLD * NUM_INITIAL )/ ( XNEW - XOLD )
            MSHIFT = XNEW * NSHIFT
            NK(K) = NK(K) - NSHIFT
            NK(KK) =NK(KK) + NSHIFT

            if(prt) then
               print*,'NSHIFT',NSHIFT, 'MSHIFT',MSHIFT
               print*,'New NK',k,NK(k),' Nk(kk)',kk,NK(kk)
               print*,'Total mass bin',k,sum(MK(k,1:icomp-idiag))
               print*,'SO4 mass bin  ',k,(MK(k,srtso4))
               print*,'Total mass bin',kk,sum(MK(kk,1:icomp-idiag))
               print*,'SO4 mass bin  ' ,kk,(MK(kk,srtso4))
            endif            

            DO J = 1, ICOMP-IDIAG
              FJ = MK(K,J)/ DRYMASS 
              MK(K,J) = XOLD * NK(K) * FJ
              MK(KK,J) = MK(KK,J) + MSHIFT * FJ
            ENDDO

            if(prt) then
               print*,'After shift mass'
               print*,'Total mass bin',k,sum(MK(k,1:icomp-idiag))
               print*,'SO4 mass bin  ',k,(MK(k,srtso4))
               print*,'Total mass bin',kk,sum(MK(kk,1:icomp-idiag))
               print*,'SO4 mass bin  ',kk,(MK(kk,srtso4))
            endif            



           ELSE
            ERRORSWITCH = .TRUE.
            PRINT *, 'MNFIX(3) : AVG>Xk(',K+1,')'
            GOTO 300
           ENDIF    ! Fixerror             
         ENDIF       ! AVG > Xk(k+1)

!         if (PRT) then     !<step5.1-temp>
!            print *,'After_Check4---------------------'
!            totmas = sum(MK(k,1:icomp-1))
!            print *, totmas,NK(k), totmas/NK(k)
!         endif

         ! If under boundary of the current bin
 111     IF ( AVG <  xk(K) ) THEN
            IF ( PRT ) PRINT *,'MNFIX [4]: AVG<Xk(',K,')'
            IF ( FIXERROR ) THEN
               !average mass is too low - shift number to lower bin
               !KK = K - 1 ! jrp potential for errors here
               KK = K
               XNEW = xk(KK)* 1.1
 200           IF ( XNEW >= AVG ) THEN
                  IF ( KK > 1 ) THEN 
                     KK = KK - 1
                     XNEW = xk(KK)* 1.1
                     GOTO 200
                  ELSE
                     ! Already reach lowest bin - must remove some number (win, 8/1/07)
                     NK(K) = DRYMASS/ ( xk(1)* 1.2 )
                     GOTO 222
                  ENDIF                     
               ENDIF
               !XOLD = SQRT(xk(K)* xk(K+1))
               XOLD = AVGMASS(k)
               NUM_INITIAL = NK(K)
               NSHIFT = NUM_INITIAL - DRYMASS/XOLD !(win, 10/20/08)
               !Prior to 10/20/08 (win)
               !NSHIFT = (DRYMASS - XOLD * NUMBER)/ ( XNEW - XOLD )
               MSHIFT = XNEW * NSHIFT
               NK(K) = NK(K) - NSHIFT
               NK(KK) = NK(KK) + NSHIFT
               DO J=1,ICOMP
                  FJ = MK(K,J)/ DRYMASS   
                  MK(K,J) = XOLD * NK(K) * FJ
                  MK(KK,J) = MK(KK,J) + MSHIFT * FJ
               ENDDO

            ELSE
               ERRORSWITCH = .TRUE.
               PRINT *, 'MNFIX(4): AVG < Xk(',k,')'
               GOTO 300
            ENDIF       
 222      ENDIF
         
!         if (PRT) then     !<step5.1-temp>
!            print *,'After_Check5---------------------'
!            totmas = sum(MK(k,1:icomp-1))
!            print *, totmas,NK(k), totmas/NK(k)
!         endif
!c         if (PRT) print *,MK(k,1),NK(k), MK(k,1)/NK(k),'Check5'!<step4.4>tmp (win, 9/28/05)


      ENDDO ! loop bin
      save5=xk(1)

      !print*,2,NK(1),NK(2)
      !print*,2,MK(1,:)
      !print*,2,MK(2,:)

      ! JRP check for neg numbers
      DO K = 1,IBINS
        IF (NK(K) < 0.d0) THEN
           print*,'4 NK < 0 in MNFIX',K,NK(K)
        ENDIF
        DO J=1,ICOMP
           IF (MK(K,J) < 0.d0) THEN
              print*,'4 MK < 0 in MNFIX',K,J,MK(K,J)
              print*,'saved xk1s',save1,save2,save3,save4,save5
              print*,'xk',xk
           ENDIF
        ENDDO
        IF ( IT_IS_NAN(NK(K)) ) THEN
           PRINT *,'4 Found Nan in Nk at bin',K
           ERRORSWITCH = .TRUE.
           print *,'4 MNFIX(0): Found NaN in Nk(,',k,')'
           GOTO 300
        ENDIF
        DO J = 1, ICOMP
           IF ( IT_IS_NAN(MK(K,J)) ) THEN
              PRINT *,'4 Found Nan in Mk at bin',K,'component',J
              ERRORSWITCH = .TRUE.
              GOTO 300
           ENDIF
        ENDDO
      ENDDO

!      if (PRT) then !<step5.1-temp>
      ! Catch any small negative values resulting from fixing
      !--------------------------------------------------------------------------
      DO K = 1, IBINS
        IF ( NK(K) < 0D0 ) THEN
           IF ( PRT ) THEN 
              PRINT *,'MNFIX[5]: FOUND NEGATIVE N'
              PRINT *, 'Bin, N', K, NK(K)
           ENDIF
           IF ( ABS(NK(K)) < 1D0 .and. FIXERROR ) THEN
              NK(K) = 0D0
              IF ( PRT ) PRINT *,'Negative N > -1.0 Reset to zero'
           ELSE
              ERRORSWITCH = .TRUE.
              PRINT *, 'MNFIX(5): Negative N after fixing at bin',k
              GOTO 300          !exit mnfix if found negative error (win, 4/18/06)
           ENDIF
        ENDIF
        DO J = 1, ICOMP
           IF ( MK(K,J) < 0D0 ) THEN
              IF ( PRT ) THEN 
                 PRINT *,'MNFIX[6]: FOUND NEGATIVE M'
                 PRINT *,'Bin, Comp, Mk', K, J, MK(K,J)
              ENDIF
              IF( ABS(MK(K,J)) < 1D-5 .and. FIXERROR ) THEN
                 MK(K,J) = 0D0
                 IF ( PRT ) PRINT *,'Negative M > -1.d-5 Reset to zero'
              ELSE
                 ERRORSWITCH =.TRUE.
                 PRINT *, 'MNFIX(6): Negative M after fixing at bin',k
                 GOTO 300       !exit mnfix if found negative error (win, 4/18/06)
              ENDIF
           ENDIF
        ENDDO                   !icomp
      ENDDO                     !ibins

      ! JRP check for neg numbers
      DO K = 1,IBINS
        IF (NK(K) < 0.d0) THEN
           print*,'5 NK < 0 in MNFIX',K,NK(K)
        ENDIF
        DO J=1,ICOMP
           IF (MK(K,J) < 0.d0) THEN
              print*,'5 MK < 0 in MNFIX',K,J,MK(K,J)
           ENDIF
        ENDDO
        IF ( IT_IS_NAN(NK(K)) ) THEN
           PRINT *,'5 Found Nan in Nk at bin',K
           ERRORSWITCH = .TRUE.
           print *,'5 MNFIX(0): Found NaN in Nk(,',k,')'
           GOTO 300
        ENDIF
        DO J = 1, ICOMP
           IF ( IT_IS_NAN(MK(K,J)) ) THEN
              PRINT *,'5 Found Nan in Mk at bin',K,'component',J
              ERRORSWITCH = .TRUE.
              GOTO 300
           ENDIF
        ENDDO
      ENDDO

      ! Check any last inconsistent M=0 or N=0 
      !--------------------------------------------------------
      DO K = 1, IBINS
         DRYMASS = 0.d0
         DO J = 1, ICOMP-IDIAG
            DRYMASS = DRYMASS + MK(K,J)
         ENDDO
         IF ( NK(K) /= 0D0 .AND. DRYMASS == 0D0 .or.
     &        NK(K) == 0D0 .AND. DRYMASS /= 0D0     ) THEN   
            PRINT *, '5.5 set nk, mk to ZERO for all bins'
            DO J = 1, ICOMP
               MK(K,J)=0.D0
               NK(K) = 0.D0
            ENDDO
            MK(K,ICOMP) = 0.d0 !Set aerosol water to zero too
         ENDIF                  ! If tiny number
      ENDDO

      ! JRP check for neg numbers
      DO K = 1,IBINS
        IF (NK(K) < 0.d0) THEN
           print*,'6 NK < 0 in MNFIX',K,NK(K)
           STOP
        ENDIF
        DO J=1,ICOMP
           IF (MK(K,J) < 0.d0) THEN
              print*,'6 MK < 0 in MNFIX',K,J,MK(K,J)
              STOP
           ENDIF
        ENDDO
      ENDDO

 300  CONTINUE
          
      IF (ERRORSWITCH) THEN
555   FORMAT (3E15.5E2)
       WRITE(6,*)'END OF MNFIX ( WHERE? )'
       WRITE(6,*)'DRYMAS-excl-NH4  NK      DRYMASS/NK'
       DO K = 1,IBINS
            TOTMAS = SUM(MK(K,1:ICOMP-1))
!            PRINT *, TOTMAS,NK(K), TOTMAS/NK(K)
       ENDDO
        WRITE(6,555)
     &          TOTMAS, NK(K),
     &          TOTMAS/ NK(K) 
        print*,'-----------'
        call debugprint( NK, MK, 0,0,0,'End of MNFIX')

!       write(*,*)'Nk'
!       write(*,*) NK(1:30)
!       write(*,*)'Mk(srtso4)'
!       write(*,*) MK(1:30,srtso4)
!       write(*,*)'Mk(srth2o)'
!       write(*,*) MK(1:30,srth2o)
       !STOP 'Negative Nk or Mk at after mnfix'  !comment out this to make it stop outside mnfix so that I can print out the i,j,l (location) of the error (win, 9/1/05)
      ENDIF

       ! Return to call subroutine
      RETURN
      END SUBROUTINE MNFIX

!-----------------------------------------------------------------------------

      SUBROUTINE SUBGRIDCOAG( NDISTINIT, NDIST, MDIST, BOXVOLUME,TEMP,
     &                        PRES,TSCALE, NDISTFINAL, MADDFINAL, pdbug)
!
!*****************************************************************************
!  Subroutine SUBGRIDCOAG determine how much of each size of freshly emitted 
!  aerosol will be scavenged by coagulation prior to being completely mixed in 
!  the gridbox and will give the new emissions size distribution along with 
!  where the mass of coagulated particles should be added.
!  Written by Jeff Pierce, December, 2006
!  Implement in GEOS-Chem by Win Trivitayanurak, 10/4/07
!
!  Varibles comments:
!  ==========================================================================
!  ndistinit(nbins)   : the number of particles being added to the gridbox 
!                       before subgrid coag
!  ndist(nbins)       : the number of particles in the box
!  mdist(nbins,icomp) : the mass of each component in the box. (kg)
!  boxvolume          : volume of box in cm3
!  tscale             : the scale time for mixing (s)
!  ndistfinal(nbins)  : the number of particles being added to the gridbox 
!                       after subgrid coag
!  maddfinal(nbins)   : the mass that should be added to each bin due to 
!                       coagulation (kg)
!
!
!  NOTES:
!*****************************************************************************

      ! Arguments
      REAL*8 ndistinit(ibins) 
      REAL*8 ndist(ibins) 
      REAL*8 mdist(ibins,icomp)
      REAL*4 boxvolume, temp , PRES
      REAL*4 tscale
      REAL*8 ndistfinal(ibins)
      REAL*8 maddfinal(ibins)
      logical pdbug  ! for pringing out during debugging

      ! Local variables
      REAL*8 mp                     ! mass of the particle (kg)
      REAL*4 density                !density (kg/m3) of particles
      REAL*8 fracdiaml(ibins,ibins) ! fraction of coagulation that occurs with each bin larger
      REAL*8 kcoag(ibins) ! the coagulation rate for the particles in each bin (s^-1)
!      REAL*4 aerodens
!      external aerodens

      REAL*4 mu                     !viscosity of air (kg/m s)
      REAL*4 mfp                    !mean free path of air molecule (m)
      REAL*4 Kn                     !Knudsen number of particle
      REAL*4 beta                   !correction for coagulation coeff.
      REAL*8 Mktot      !total mass of aerosol
      REAL*4 kij(ibins,ibins)
      REAL*4 Dpk(ibins)             !diameter (m) of particles in bin k
      REAL*4 Dk(ibins)              !Diffusivity (m2/s) of bin k particles
      REAL*4 ck(ibins)              !Mean velocity (m/2) of bin k particles
      REAL*4 neps
      REAL*4 meps
      INTEGER I, J, K, KK
      LOGICAL ERRORSWITCH

      ! Adjustable parameters
      real*4 pi, kB               !kB is Boltzmann constant (J/K)
      real*4 R       !gas constant (J/ mol K)
      parameter (pi=3.141592654, kB=1.38e-23, R=8.314)
      parameter (neps=1E8, meps=1E-8)

      !=================================================================
      ! SUBGRIDCOAG begins here!
      !=================================================================

      if (pdbug) call debugprint(Ndist,Mdist,0,0,0,
     &                           'SUBDGRIDCOAG: Entering')


      !Before going in to calculation, check and fix Nk and Mk
      ERRORSWITCH = .FALSE.
#if defined ( DEBUG ) 
!      print *, 'mnfix in tomas_mod:7495'
#endif
      CALL MNFIX( NDIST, MDIST, ERRORSWITCH )
      IF ( ERRORSWITCH ) THEN
         PRINT *,'SUBGRIDCOAG: MNFIX found error: Entering SUBGRIDCOAG'
         PDBUG = .TRUE.
         GOTO 11
      ENDIF
         
      if (pdbug) call debugprint(Ndist,Mdist,0,0,0,
     &                           'SUBDGRIDCOAG: after MNFIX_1')

      mu=2.5277e-7*temp**0.75302
      mfp=2.0*mu / ( pres*sqrt( 8.0 * 0.0289 / (pi*R*temp) ) )  !S&P eqn 9.6
C Calculate particle sizes and diffusivities
      do k=1,ibins
         Mktot=0.2*mdist(k,srtso4)
         do j=1, icomp
            Mktot=Mktot+mdist(k,j)
         enddo
         if (Mktot.gt.meps)then
            density=aerodens(mdist(k,srtso4),0d0,
     &           0.1875d0*mdist(k,srtso4),mdist(k,srtnacl),
     &           mdist(k,srtecil),mdist(k,srtecob),
     &           mdist(k,srtocil),mdist(k,srtocob),mdist(k,srtdust),
     &           mdist(k,srth2o)) !assume bisulfate
         else
            density = 1400.
         endif
         if(ndist(k).gt.neps .and. Mktot.gt.meps)then
            mp=Mktot/ndist(k)
         else
            mp=sqrt(xk(k)*xk(k+1))
         endif
         Dpk(k)=((mp/density)*(6./pi))**(0.333)
         Kn=2.0*mfp/Dpk(k)                            !S&P Table 12.1
         Dk(k)=kB*temp/(3.0*pi*mu*Dpk(k))             !S&P Table 12.1
     &   *((5.0+4.0*Kn+6.0*Kn**2+18.0*Kn**3)/(5.0-Kn+(8.0+pi)*Kn**2))
         ck(k)=sqrt(8.0*kB*temp/(pi*mp))              !S&P Table 12.1
      enddo

C Calculate coagulation coefficients

      do i=1,ibins
         do j=1,ibins
            Kn=4.0*(Dk(i)+Dk(j))          
     &        /(sqrt(ck(i)**2+ck(j)**2)*(Dpk(i)+Dpk(j))) !S&P eqn 12.51
            beta=(1.0+Kn)/(1.0+2.0*Kn*(1.0+Kn))          !S&P eqn 12.50
            !This is S&P eqn 12.46 with non-continuum correction, beta
            kij(i,j)=2.0*pi*(Dpk(i)+Dpk(j))*(Dk(i)+Dk(j))*beta
            kij(i,j)=kij(i,j)*1.0e6/boxvolume  !normalize by grid cell volume
         enddo
      enddo

C     get the first order loss rate
      kcoag(ibins)=0.0
!debug
      if(pdbug) print *,'Bin  KCOAG'
      do k=1,ibins-1
         kcoag(k)=0.0
         do kk=k+1,ibins
            kcoag(k)=kcoag(k)+kij(k,kk)*ndist(kk)
         enddo
!debug
         if(pdbug) print *, k, kcoag(k)
      enddo

C     get the fraction of the coagulation that occurs from each bin larger
      do k=1,ibins
         do kk=1,ibins
            fracdiaml(k,kk)=0.
         enddo
      enddo
      do k=1,ibins-1
!debug
         if(pdbug) print *, 'Bin k', k
!debug
         if(pdbug) print *, 'Bin kk   fracdiaml(k,kk)'         
         do kk=k+1,ibins
            if (kcoag(k).gt.0.d0)then
               fracdiaml(k,kk)=kij(k,kk)*ndist(kk)/kcoag(k)
            else
               fracdiaml(k,kk)=0
            endif
!debug
            if(pdbug) print *, kk, fracdiaml(k,kk)
         enddo
      enddo

C     determine the number of new particles left after coagulation
      do k=1,ibins
         ndistfinal(k)=ndistinit(k)*exp(-kcoag(k)*tscale)
      enddo

C     determine the mass added to each bin coagulation
      do k=1,ibins
         maddfinal(k)=0.
      enddo
      do k=1,ibins-1
         do kk=k+1,ibins
            maddfinal(kk)=maddfinal(kk) + (ndistinit(k)-ndistfinal(k))*
     &           fracdiaml(k,kk)*sqrt(xk(k)*xk(k+1))
         enddo
      enddo

      pdbug = .FALSE.

 11   continue
      return
      
      END SUBROUTINE SUBGRIDCOAG
      
!-----------------------------------------------------------------------------
      
      SUBROUTINE STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd )
!
!*****************************************************************************
!  Subroutine STORENM stores values of Nk and Mk into Nkd and Mkd for 
!  diagnostic purposes.  Also do gas phase concentrations. (win, 7/23/07)
!*****************************************************************************

      ! Local variables
      INTEGER             :: J, K

      REAL*8,INTENT(IN)    :: Nk(IBINS)
      REAL*8,INTENT(IN)    :: Mk(IBINS, ICOMP)
      REAL*8,INTENT(IN)    :: Gc(ICOMP - 1)
      REAL*8,INTENT(OUT)   :: Nkd(IBINS)
      REAL*8,INTENT(OUT)   :: Mkd(IBINS,ICOMP)
      REAL*8,INTENT(OUT)   :: Gcd(ICOMP - 1)

!sfarina 

      DO J= 1, ICOMP-1
         Gcd(J)=Gc(J)
      ENDDO
      DO K = 1, IBINS
         Nkd(K)=Nk(K)
         DO J= 1, ICOMP
            Mkd(K,J)=Mk(K,J)
         ENDDO
      ENDDO

      ! Return to calling subroutine
      RETURN
      END SUBROUTINE STORENM

!----------------------------------------------------------------------------

      SUBROUTINE DEBUGPRINT( NK, MK, I,J,L, LOCATION)
!
!*****************************************************************************
!  Subroutine DEBUGPRINT print out the Nk and Mk values for error checking
!  (win, 9/30/08)
!*****************************************************************************
!


      ! Arguments
      REAL*8,           INTENT(IN) :: Nk(IBINS), MK(IBINS,ICOMP)
      INTEGER,          INTENT(IN) :: I,J,L
      CHARACTER(LEN=*), INTENT(IN) :: LOCATION

      INTEGER  :: JC, k

      WRITE(*,*) LOCATION, I, J, L
!      write(6,*) 'Nk(1:30)'
!      write(6,*) Nk(1:30)
!      do jc=1,icomp
!         write(6,*) 'Mk(1:30) comp:',jc
!         write(6,*) Mk(1:30,jc)
!      enddo
      write(*,111) 'Bin        Num         SO4        NaCl        ',
     &             'ECIL        ECOB        OCIL       OCOB       ',
     &             'Dust        NH4         Water  '
      DO K = 1, IBINS
         write(*,110) k,Nk(k), Mk(k, 1:icomp)
      ENDDO
      write(*,*) ' '
 110  FORMAT ( I2, 10E12.5 )
 111  FORMAT (a,a,a)

      END SUBROUTINE DEBUGPRINT

!-----------------------------------------------------------------------------
      
      SUBROUTINE NH4BULKTOBIN( MSULF, NH4B, MAMMO )
!
!*****************************************************************************
!  Subroutine NH4BULKTOBIN takes the bulk ammonium aerosol from GEOS-Chem and
!  fraction it to each bin according to sulfate mole fraction in each bin
!  Written by Win Trivityanurak, Sep 26, 2008
!
!  Make sure that we work with mass or mass conc.
!*****************************************************************************
!

      ! Arguments
      REAL*8,  INTENT(IN)   :: MSULF(IBINS)  ! size-resolved sulfate [kg]
      REAL*8,  INTENT(IN)   :: NH4B          ! Bulk NH4 mass [kg]
      REAL*8,  INTENT(OUT)  :: MAMMO(IBINS)  ! size-resolved NH4 [kg]

      ! Local variables
      INTEGER               :: K
      REAL*8                :: TOTMASS, NH4TEMP

      !=================================================================
      ! NH4BULKTOBIN begins here
      !=================================================================

      MAMMO(:) = 0.d0

      ! Sum the total sulfate
      TOTMASS = 0.d0
      DO K = 1, IBINS
         TOTMASS = TOTMASS + MSULF(K)
      ENDDO

      IF ( TOTMASS .eq. 0.d0 ) RETURN

      ! Limit the amount of NH4 entering TOMAS calculation
      ! if it is very NH4-rich, just limit the amount to balance
      ! existing 30-bin-summed SO4 assuming (NH4)2SO4 in such case
      !  (NH4)2 mass = (SO4)mass / 96. * 2. * 18. = 0.375*(SO4)mass
      ! (win, 9/28/08)
      NH4TEMP = NH4B
      IF ( NH4B/TOTMASS > 0.375d0 )   !make sure we use mass ratio
     &     NH4TEMP = 0.375d0 * TOTMASS


      ! Calculate ammonium aerosol scale to each bin
      DO K = 1, IBINS
         MAMMO(K) = MSULF(K) / TOTMASS * NH4TEMP
      ENDDO


!      write(777,*) NH4B/TOTMASS
      
      RETURN
      END SUBROUTINE NH4BULKTOBIN
      
!------------------------------------------------------------------------------
      
      FUNCTION AERODENS( MSO4, MNO3, MNH4, MNACL, MECIL, MECOB, MOCIL, 
     &                   MOCOB, MDUST, MH2O )  RESULT( VALUE )
!
!******************************************************************************
!  Function AERODENS calculates the density (kg/m3) of a sulfate-nitrate-
!  ammonium-nacl-OC-EC-dust-water mixture.  Inorganic mass (sulfate-nitrate-
!  ammonium-nacl-water) is assumed to be internally mixed.  Then the density
!  of inorg and EC, OC, and dust is combined weighted by mass.  
!  WRITTEN BY Peter Adams, May 1999 in GISS GCM-II' and extened to include 
!  carbonaceous aerosol in Jan, 2002.
!
!  NOTES:
!  (1 ) Add error check to prevent division by zero (win, 3/13/08) 
!******************************************************************************
!
      ! Arguments
      REAL*8,  INTENT(IN)  ::  MSO4, MNO3, MNH4, MNACL, MH2O
      REAL*8,  INTENT(IN)  ::  MECIL, MECOB, MOCIL, MOCOB, MDUST

      ! Return value
      REAL*8                  :: VALUE
      
      ! Local variables
      real*8                  :: IDENSITY, DEC, DOC, DDUST, MTOT
      parameter(dec=2200., doc=1400., ddust=2650.)

      !=================================================================
      ! AERODENS begins here!
      !=================================================================
      
      IDENSITY = INODENS( MSO4, MNO3, MNH4, MNACL, MH2O )
      MTOT = MSO4+MNO3+MNH4+MNACL+MH2O+MECIL+MECOB+MOCIL+MDUST+MOCOB
      IF ( MTOT > 0.d0 ) THEN 
         VALUE = ( IDENSITY*(MSO4+MNO3+MNH4+MNACL+MH2O) +
     &             DEC*(MECIL+MECOB) + DOC*(MOCIL+MOCOB)+
     &             DDUST*MDUST                            )/MTOT
      ELSE
         VALUE = 1400.
      ENDIF

      END FUNCTION AERODENS

!------------------------------------------------------------------------------

      FUNCTION INODENS( MSO4_, MNO3_, MNH4_, MNACL_, MH2O_ ) 
     &     RESULT( VALUE )
!
!******************************************************************************
!  Function INODENS calculates the density (kg/m3) of a sulfate-nitrate-
!  ammonium-nacl-water mixture that is assumed to be internally mixed.  
!  WRITTEN BY Peter Adams, May 1999 in GISS GCM-II' 
!  Introduced to GEOS-CHEM by Win Trivitayanurak (win@cmu.edu) 8/6/07 first
!  as AERODENS, then change to INODENS on 9/3/07
!
!     mso4, mno3, mnh4, mh2o, mnacl - These are the masses of each aerosol
!     component.  Since the density is an intensive property,
!     these may be input in a variety of units (ug/m3, mass/cell, etc.).
!
!-----Literature cited--------------------------------------------------
!     I. N. Tang and H. R. Munkelwitz, Water activities, densities, and
!       refractive indices of aqueous sulfates and sodium nitrate droplets
!       of atmospheric importance, JGR, 99, 18,801-18,808, 1994
!     Ignatius N. Tang, Chemical and size effects of hygroscopic aerosols
!       on light scattering coefficients, JGR, 101, 19,245-19,250, 1996
!     Ignatius N. Tang, Thermodynamic and optical properties of mixed-salt
!       aerosols of atmospheric importance, JGR, 102, 1883-1893, 1997
!  NOTES:
!  ( 1) The function originally is AERODENS (in sulfate & sea-salt version)
!       got changed the name to INODENS when OC/EC aerosols are added. AERODENS
!       then call for this INODENS for calculation of inorganic aerosol
!       density (win, 9/3/07)
!******************************************************************************
!
      ! Arguments
      REAL*8,  INTENT(IN)  ::  MSO4_, MNO3_, MNH4_, MNACL_, MH2O_

      ! Return value
      REAL*8                  :: VALUE

      ! Local variables
      real*8 MSO4, MNO3, MNH4, MNACL, MH2O
c      real*8 so4temp, no3temp, nh4temp, nacltemp, h2otemp  
      real*8 mwso4, mwno3, mwnh4, mwnacl, mwh2o            !molecular weights
      real*8 ntot, mtot                          !total number of moles, mass
      real*8 nso4, nno3, nnh4, nnacl, nh2o       !moles of each species
      real*8 xso4, xno3, xnh4, xnacl, xh2o       !mole fractions
      real*8 rso4, rno3, rnh4, rnacl, rh2o       !partial molar refractions
      real*8 ran, rs0, rs1, rs15, rs2       !same, but for solute species
      real*8 asr                            !ammonium/sulfate molar ratio
      real*8 nan, ns0, ns1, ns15, ns2, nss  !moles of dry solutes (nss = sea salt)
      real*8 xan, xs0, xs1, xs15, xs2, xss  !mass % of dry solutes - Tang (1997) eq. 10
      real*8 dan, ds0, ds1, ds15, ds2, dss  !binary solution densities - Tang (1997) eq. 10
      real*8 mwan, mws0, mws1, mws15, mws2  !molecular weights
      real*8 yan, ys0, ys1, ys15, ys2, yss  !mole fractions of dry solutes
      real*8 yh2o
      real*8 d                              !mixture density
      real*8 xtot

C     In the lines above, "an" refers to ammonium nitrate, "s0" to 
C     sulfuric acid, "s1" to ammonium bisulfate, and "s2" to ammonium sulfate.
C     "nacl" or "ss" is sea salt.

      parameter(mwso4=96.d0, mwno3=62.d0, mwnh4=18.d0, mwh2o=18.d0, 
     &          mwnacl=58.45d0)
      parameter(mwan=mwnh4+mwno3, mws0=mwso4+2.d0, mws1=mwso4+1.d0+mwnh4
     &          ,mws2=2.d0*mwnh4+mwso4)

      !=================================================================
      ! INODENS begins here!
      !=================================================================

C Pass initial component masses to local variables
      mso4=mso4_
      mno3=mno3_
      mnh4=mnh4_
      mnacl=mnacl_
      mh2o=mh2o_

c      so4temp=mso4
c      no3temp=mno3
c      nh4temp=mnh4
c      h2otemp=mh2o
c      nacltemp=mnacl
      
      !<step4.7> if the aerosol mass is zero - then just return the 
      !typical density = 1500 kg/m3 (win, 1/4/06)
      if (mso4 .eq. 0.d0 .and. mno3 .eq.0.d0 .and. mnh4.eq.0.d0 
     &     .and. mnacl .eq. 0.d0 ) then
         VALUE = 1500.d0 !kg/m3
         goto 10
      endif

C Calculate mole fractions
      mtot = mso4+mno3+mnh4+mnacl+mh2o
      nso4 = mso4/mwso4
      nno3 = mno3/mwno3
      nnh4 = mnh4/mwnh4
      nnacl = mnacl/mwnacl
      nh2o = mh2o/mwh2o
      ntot = nso4+nno3+nnh4+nnacl+nh2o
      xso4 = nso4/ntot
      xno3 = nno3/ntot
      xnh4 = nnh4/ntot
      xnacl = nnacl/ntot
      xh2o = nh2o/ntot

C If there are more moles of nitrate than ammonium, treat unneutralized
C HNO3 as H2SO4
      if (nno3 .gt. nnh4) then
         !make the switch
         nso4=nso4+(nno3-nnh4)
         nno3=nnh4
         mso4=nso4*mwso4
         mno3=nno3*mwno3

         !recalculate quantities
         mtot = mso4+mno3+mnh4+mnacl+mh2o
         nso4 = mso4/mwso4
         nno3 = mno3/mwno3
         nnh4 = mnh4/mwnh4
         nnacl = mnacl/mwnacl
         nh2o = mh2o/mwh2o
         ntot = nso4+nno3+nnh4+nnacl+nh2o
         xso4 = nso4/ntot
         xno3 = nno3/ntot
         xnh4 = nnh4/ntot
         xnacl = nnacl/ntot
         xh2o = nh2o/ntot

      endif

C Calculate the mixture density
C Assume that nitrate exists as ammonium nitrate and that other ammonium
C contributes to neutralizing sulfate
      nan=nno3
      if (nnh4 .gt. nno3) then
         !extra ammonium
         asr=(nnh4-nno3)/nso4
      else
         !less ammonium than nitrate - all sulfate is sulfuric acid
         asr=0.d0
      endif
      if (asr .ge. 2.d0) asr=2.d0
      if (asr .ge. 1.d0) then
         !assume NH4HSO4 and (NH4)2(SO4) mixture
         !NH4HSO4
         ns1=nso4*(2.d0-asr)
         !(NH4)2SO4
         ns2=nso4*(asr-1.d0)
         ns0=0.d0
      else
         !assume H2SO4 and NH4HSO4 mixture
         !NH4HSO4
         ns1=nso4*asr
         !H2SO4
         ns0=nso4*(1.d0-asr)
         ns2=0.d0
      endif

      !Calculate weight percent of solutes
      xan=nan*mwan/mtot*100.d0
      xs0=ns0*mws0/mtot*100.d0
      xs1=ns1*mws1/mtot*100.d0
      xs2=ns2*mws2/mtot*100.d0
      xnacl=nnacl*mwnacl/mtot*100.d0
      xtot=xan+xs0+xs1+xs2+xnacl

      !Calculate binary mixture densities (Tang, eqn 9)
      dan=0.9971d0 +4.05d-3*xtot +9.0d-6*xtot**2.d0
      ds0=0.9971d0 +7.367d-3*xtot -4.934d-5*xtot**2.d0 
     &     +1.754d-6*xtot**3.d0 - 1.104d-8*xtot**4.d0
      ds1=0.9971d0 +5.87d-3*xtot -1.89d-6*xtot**2.d0 
     &     +1.763d-7*xtot**3.d0
      ds2=0.9971d0 +5.92d-3*xtot -5.036d-6*xtot**2.d0 
     &     +1.024d-8*xtot**3.d0
      dss=0.9971d0 +7.41d-3*xtot -3.741d-5*xtot**2.d0 
     &     +2.252d-6*xtot**3.d0   -2.06d-8*xtot**4.d0

      !Convert x's (weight percent of solutes) to fraction of dry solute (scale to 1)
      xtot=xan+xs0+xs1+xs2+xnacl
      xan=xan/xtot
      xs0=xs0/xtot
      xs1=xs1/xtot
      xs2=xs2/xtot
      xnacl=xnacl/xtot

      !Calculate mixture density
      d=1./(xan/dan+xs0/ds0+xs1/ds1+xs2/ds2+xnacl/dss)  !Tang, eq. 10

      if ((d .gt. 2.d0) .or. (d .lt. 0.997d0)) then
         write(*,*) 'ERROR in aerodens'
         write(*,*) mso4,mno3,mnh4,mnacl,mh2o
         print *, 'xtot',xtot
         print *, 'xs1',xs1, 'ns1',ns1,'mtot',mtot,'asr',asr
         write(*,*) 'density(g/cm3)',d
         STOP
      endif

C Restore masses passed
c      mso4=so4temp
c      mno3=no3temp
c      mnh4=nh4temp
c      mnacl=nacltemp
c      mh2o=h2otemp

C Return the density
      VALUE = 1000.d0*d    !Convert g/cm3 to kg/m3

      !<step4.7> negative value check (win, 1/4/06)
      if ( VALUE < 0d0 ) then
         print *, 'ERROR :: aerodens - negative', VALUE
         STOP
      endif
 
 10   CONTINUE

      END FUNCTION INODENS

!------------------------------------------------------------------------------

      FUNCTION DMDT_INT ( M0, TAU, WR ) RESULT( VALUE )
!
!******************************************************************************
!  Function DMDT_INT apply the analytic solution to the droplet growth equation
!  in mass space for a given scale length which mimics the inclusion of gas 
!  kinetic effects.
!  (win, 7/23/07)
!  Originally written by Peter Adams
!  Modified for GEOS-CHEM by Win Trivitayanurak (win@cmu.edu)
!  
!  Variables
!  (1)  M0  initial mass
!  (2)  L0  length scale
!  (3)  Tau forcing from vapor field
!
!  Original note:
Cpja I have changed the length scale.  Non-continuum effects are
Cpja assumed to be taken into account in choice of tau (in so4cond
Cpja subroutine).

Cpja I have also added another argument to the function call, WR.  This
Cpja is the ratio of wet mass to dry mass of the particle.  I use this
Cpja information to calculate the amount of growth of the wet particle,
Cpja but then return the resulting dry mass.  This is the appropriate
Cpja way to implement the condensation algorithm in a moving sectional
Cpja framework.

!Reference: Stevens et al. 1996, Elements of the Microphysical Structure
!           of Numerically Simulated Nonprecipitating Stratocumulus,
!           J. Atmos. Sci., 53(7),980-1006. 
! This calculates a solution for m(t+dt) using eqn.(A3) from the reference


      ! Arguments
      REAL*8,   INTENT(IN)  ::  M0,  TAU,  WR

      ! Return value
      REAL*8                :: VALUE

      ! Local variables
      REAL*8                ::  X,  L0,  C,  ZERO,  MH2O
      PARAMETER (C=2.d0/3.d0,L0=0.0d0,ZERO=0.0d0)

      !=================================================================
      ! DMDT_INT begins here!
      !=================================================================
 
      MH2O = ( WR - 1.d0 ) * M0
      X = ( ( M0 + MH2O ) ** C + L0 )
      X = MAX( ZERO, SQRT(MAX(ZERO,C*TAU+X))-L0 )

!<step5.3> Do aqueous oxidation dry - so no need to select process (win, 7/14/06)
      !<step5.3> For so4cond condensation, use constant water amount.
      ! For aqueous oxidation, use constant wet ratio. (win, 7/13/06)
!prior to 10/2/08
!      VALUE = X * X * X - MH2O
!!         DMDT_INT = X*X*X/WR    !<step5.2> change calculation to keep WR constant after condensation/evap (win, 5/14/06)

      !<step6.3> bring back the previously reverted back (win, 10/2/08)
      VALUE = X*X*X/WR 
Cpja Perform some numerical checks on dmdt_int
      IF ((TAU > 0.0) .and. (VALUE < M0)) VALUE = M0
      IF ((TAU < 0.0) .and. (VALUE > M0)) VALUE = M0


      ! Return to calling program
      END FUNCTION DMDT_INT
 
!------------------------------------------------------------------------------

      FUNCTION GASDIFF( TEMP, PRES, MW, SV ) RESULT( VALUE )
!
!******************************************************************************
!  Function GASDIFF returns the diffusion constant of a species in air (m2/s).
!  It uses the method of Fuller, Schettler, and Giddings as described in 
!  Perry's Handbook for Chemical Engineers.
!  WRITTEN BY Peter Adams, May 2000
!******************************************************************************

      ! Arguments
      real temp, pres  !temperature (K) and pressure (Pa) of air
      real mw          !molecular weight (g/mol) of diffusing species
      real Sv          !sum of atomic diffusion volumes of diffusing species
      real VALUE      

      ! Local variables
       real mwair, Svair   !same as above, but for air
       real mwf, Svf
       parameter(mwair=28.9, Svair=20.1)

      !========================================================================
      ! GASDIFF begins here!
      !========================================================================

       mwf=sqrt((mw+mwair)/(mw*mwair))
       Svf=(Sv**(1./3.)+Svair**(1./3.))**2.
       VALUE =1.0e-7*temp**1.75*mwf/pres*1.0e5/Svf
      

       END FUNCTION GASDIFF

!-----------------------------------------------------------------------------

      FUNCTION GETDP( I, J, L, N, State_Met, State_Chm ) RESULT( VALUE )
!
!*****************************************************************************
!  Function GETDP calculate multi-component aerosol diameter
!
!  Originally written by Peter Adams in GISS GCM-II'
!  Use in GEOS-CHEM v5-07-08 and later by Win Trivitayanurak (win@cmu.edu)
!
!  Arguments as Input : 
!  ===========================================================================
!  (1 ) I  (INTEGER) : Grid location
!  (2 ) J  (INTEGER) : Grid location
!  (3 ) L  (INTEGER) : Grid location
!  (4 ) N  (INTEGER) : Tracer index
!
!  NOTES: 
!  (1 ) Each component mass has a unit of kg/grid box
!  (2 ) Total mass assumed to be ammonium bisulfate including water
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!*****************************************************************************

      USE ERROR_MOD,          ONLY : ERROR_STOP, IT_IS_NAN
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TRACERID_MOD,       ONLY : IDTNK1,  IDTAW1

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND60

      implicit none

      ! Arguments
      INTEGER,   INTENT(IN)         :: I, J, L, N
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object

      ! Return value
      REAL*8                :: VALUE

      ! Local variables
      INTEGER               :: NUMBIN,  ID,   JC
      REAL*8                :: MSO4, MNO3, MNH4, MH2O, MNACL 
      REAL*8                :: MECIL, MECOB, MOCIL, MOCOB, MDUST
      REAL*8                :: DENSITY !density (kg/m3) of current size bin
      REAL*8                :: TOTALMASS !(kg)
      REAL*8                :: MCONC, NCONC  

!      real*8, external :: aerodens

      REAL*8                ::  pi
      parameter (pi=3.141592654d0)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER   :: STT(:,:,:,:)

      !=================================================================
      ! GETDP begins here!
      !=================================================================
      
      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      !-------------------------------------------------------------
      ! Calculate bin that we're working with
      !-------------------------------------------------------------
      NUMBIN = MOD(N-IDTNK1+1,IBINS)
      IF (NUMBIN==0) NUMBIN = IBINS
      ID = IDTNK1-1+NUMBIN   !ID = tracer ID of number at current bin

      !-------------------------------------------------------------
      ! Calculate aerosol water in case it has not been initialized elsewhere 
      !-------------------------------------------------------------
      CALL EZWATEREQM2( I, J, L, NUMBIN, State_Met, State_Chm )

      !-------------------------------------------------------------
      ! Check negative STT
      !-------------------------------------------------------------
            ! Significance limit in concentration unit
            ! Treshold for mass concentration 1d-4 ug/m3 = 1d-13 kg/m3
            ! Treshold for numb concentration 1d-1 #/cm3 = 1d5 kg/m3 (fake kg = #)
      IF( STT(I,J,L,ID) == 0d0 ) GOTO 10

      IF( STT(I,J,L,ID) < 0d0 ) THEN
         NCONC = ABS( STT(I,J,L,ID) )/ State_Met%AIRVOL(I,J,L)/1D6
         IF ( nconc <= 1d-1 ) THEN
            STT(I,J,L,ID) = 0d0
         ELSE
            PRINT *,'#### GETDP: negative NK at',I,J,L,'bin',NUMBIN
            PRINT *,'Tracer',N,'STT=',STT(I,J,L,ID)
            CALL ERROR_STOP('Negative NK', 'GETDP:1')
         ENDIF
      ENDIF
      IF(IT_IS_NAN(STT(I,J,L,ID))) PRINT *,'+++++++++ Found Nan in ' ,
     &     'GETDP at (I,J,L)',I,J,L,'Bin',NUMBIN,': Nk'
      DO JC = 1, ICOMP-IDIAG 
         IF( STT(I,J,L,ID+JC*IBINS) < 0d0 ) THEN
            MCONC = ( ABS(STT(I,J,L,ID+JC*IBINS)) * 1.D9 /
     &              State_Met%AIRVOL(I,J,L) )
            IF ( MCONC <= 1.d-4 ) THEN
               STT(I,J,L,ID+JC*IBINS) = 0d0
            ELSE
               PRINT *,'#### GETDP: negative mass at',I,J,L,'bin',NUMBIN
               PRINT *,'Tracer',N,'STT=',STT(I,J,L,ID+JC*IBINS)
               CALL ERROR_STOP('Negative mass','GETDP:2')
            ENDIF  
         ENDIF
         IF(IT_IS_NAN(STT(I,J,L,ID+JC*IBINS))) PRINT *,'+++++++++ ',
     &     'Found Nan in GETDP at (I,J,L)',I,J,L,'Bin',NUMBIN,'comp',JC
      ENDDO

      !-------------------------------------------------------------
      ! Begin calculation of diameter
      !-------------------------------------------------------------

      ! Totalmass is the total mass per particle (including water and ammonia)
      ! The factor of 0.1875 is the proportion of nh4 to make the particle
      ! ammonium bisulfate
      MSO4=0.d0
      MNACL=0.d0
      MH2O=0.d0
      MECIL=0.d0
      MECOB=0.d0
      MOCIL=0.d0
      MOCOB=0.d0
      MDUST = 0.d0

      ! Get aerosol masses from GEOS-CHEM's STT array
      DO JC = 1, ICOMP-IDIAG
         IF( JC == SRTSO4  ) MSO4  = STT(I,J,L,ID+JC*IBINS)
         IF( JC == SRTNACL ) MNACL = STT(I,J,L,ID+JC*IBINS)
         IF( JC == SRTECIL ) MECIL = STT(I,J,L,ID+JC*IBINS)
         IF( JC == SRTECOB ) MECOB = STT(I,J,L,ID+JC*IBINS)
         IF( JC == SRTOCIL ) MOCIL = STT(I,J,L,ID+JC*IBINS)
         IF( JC == SRTOCOB ) MOCOB = STT(I,J,L,ID+JC*IBINS)
         IF( JC == SRTDUST ) MDUST = STT(I,J,L,ID+JC*IBINS)
      ENDDO
      MH2O  = STT(I,J,L,IDTAW1-1+NUMBIN)

!dbg      print *,'mh2o',mh2o,'at',i,j,l

      MNO3 = 0.d0
      MNH4 = MSO4 * 1.875d-1
      TOTALMASS = ( MSO4 + MNO3 + MNH4 + MNACL + MH2O + 
     &              MECIL + MECOB + MOCIL + MOCOB + MDUST)/
     &             STT(I,J,L,ID)
      DENSITY = AERODENS( MSO4, MNO3, MNH4, MNACL, MECIL, MECOB, 
     &                    MOCIL, MOCOB, MDUST, MH2O)

      VALUE = ( TOTALMASS* 6.d0/ DENSITY/ PI )**(1.d0/3.d0) !getdp [=] meter

      GOTO 20

      !if number and mass is zero - calculate dp based on the density=1500 kg/m3
 10   CONTINUE
      TOTALMASS = 1.414d0 * Xk(NUMBIN)  ! Mid-bin mass per particle
      VALUE = ( TOTALMASS* 6.d0/ 1500.d0/ PI )**(1.d0/3.d0) !getdp [=] meter

 20   CONTINUE

      IF( IT_IS_NAN( VALUE )) 
     &     CALL ERROR_STOP('Result is NaN', 'GETDP:3')

      ! Free pointer
      NULLIFY( STT )

      ! Return to calling subroutine
      END FUNCTION GETDP

!-----------------------------------------------------------------------------

      FUNCTION SPINUP( DAYS ) RESULT( VALUE )
!
!*****************************************************************************
!  Function SPINUP retuns .TRUE. or .FALSE. whether or not the current time
!  in the run have passed the spin-up period.  This would be used to determine
!  if certain errors should be fixed and let slipped or to stop a run with
!  an error message.  (win, 8/2/07)
!  ====> Be cautious that TIMEBEGIN should be changed according to 
!         whatever your spin-up beginning time is
!  Example of TIMEBEGIN (in julian time)
!         2001/07/01 = 144600.0
!         2000/11/01 = 138792.0
!          
!  NOTES:
!  (1 ) Allow choices of spin-up period depending on where the check condition
!       is used, e.g. 2 weeks, 2 days. (win, 8/2/07)
!*****************************************************************************
!
      ! References to F90 modules
      USE TIME_MOD,     ONLY : GET_TAU , GET_TAUb 

      ! Arguments
      REAL*4,    INTENT(IN) :: DAYS   ! Spin-up duration (day)

      ! Return value
      LOGICAL               :: VALUE

      ! Local variables
      REAL*4                 :: TIMENOW, TIMEBEGIN, TIMEINIT, HOURS

      !========================================================================
      ! SPINUP begins here!
      !========================================================================
      
      TIMENOW   = GET_TAU()   ! Current time in the run (Julian time) (hrs)
      TIMEBEGIN = GET_TAUb()  ! Begin time of this run (hrs)
      TIMEINIT  = 141000. !2/1/2001    ! Start time for spin-up (hrs)
      HOURS = DAYS * 24.0     ! Period allow error to pass (hrs)
      
      ! Criteria to let error go or to terminate the run
!      IF ( TIMENOW > MIN( TIMEBEGIN, TIMEINIT ) + HOURS  ) THEN
      IF ( TIMENOW > TIMEBEGIN + HOURS  ) THEN
         VALUE = .FALSE.
      ELSE
         VALUE = .TRUE.
      ENDIF

      ! Return to calling subroutine
      END FUNCTION SPINUP

!-----------------------------------------------------------------------------

      FUNCTION STRATSCAV( DP ) RESULT( VALUE )
!
!*****************************************************************************
!  Function STRATSCAV is basically a lookup table to get the below-cloud 
!     scavenging rate (per mm of rainfall) as a function of particle
!     diameter.  The data are taken from Dana, M. T., and
!     J. M. Hales, Statistical Aspects of the Washout of Polydisperse
!     Aerosols, Atmos. Environ., 10, 45-50, 1976.  I am using the
!     monodisperse aerosol curve from Figure 2 which assumes a
!     lognormal distribution of rain drops with Rg=0.02 cm and a
!     sigma of 1.86, values typical of a frontal rain spectrum
!     (stratiform clouds).
!     WRITTEN BY Peter Adams, January 2001
!     Intoduced to GEOS-Chem by Win Trivitayanurak, 8/6/07
!
!*****************************************************************************

      ! Arguments 
      REAL*8,   INTENT(IN)  :: DP  !particle diameter (m)
      
      ! Return value
      REAL*4                  :: VALUE

      ! Local variables
      integer numpts  !number of points in lookup table
      real dpdat      !particle diameter in lookup table (m)
      real scdat      !scavenging rate in lookup table (mm-1)
      integer n1, n2  !indices of nearest data points
      parameter(numpts=37)
      dimension dpdat(numpts), scdat(numpts)

      data dpdat/ 2.0E-09, 4.0E-09, 6.0E-09, 8.0E-09, 1.0E-08,
     &            1.2E-08, 1.4E-08, 1.6E-08, 1.8E-08, 2.0E-08,
     &            4.0E-08, 6.0E-08, 8.0E-08, 1.0E-07, 1.2E-07,
     &            1.4E-07, 1.6E-07, 1.8E-07, 2.0E-07, 4.0E-07,
     &            6.0E-07, 8.0E-07, 1.0E-06, 1.2E-06, 1.4E-06,
     &            1.6E-06, 1.8E-06, 2.0E-06, 4.0E-06, 6.0E-06,
     &            8.0E-06, 1.0E-05, 1.2E-05, 1.4E-05, 1.6E-05,
     &            1.8E-05, 2.0E-05/

      data scdat/ 6.99E-02, 2.61E-02, 1.46E-02, 9.67E-03, 7.07E-03,
     &            5.52E-03, 4.53E-03, 3.87E-03, 3.42E-03, 3.10E-03,
     &            1.46E-03, 1.08E-03, 9.75E-04, 9.77E-04, 1.03E-03,
     &            1.11E-03, 1.21E-03, 1.33E-03, 1.45E-03, 3.09E-03,
     &            4.86E-03, 7.24E-03, 1.02E-02, 1.36E-02, 1.76E-02,
     &            2.21E-02, 2.70E-02, 3.24E-02, 4.86E-01, 8.36E-01,
     &            1.14E+00, 1.39E+00, 1.59E+00, 1.75E+00, 1.85E+00, 
     &            1.91E+00, 1.91E+00/

      !=================================================================
      ! STRATSCAV begins here!
      !=================================================================
      
C If particle diameter is in bounds, interpolate to find value
      if ((dp .gt. dpdat(1)) .and. (dp .lt. dpdat(numpts))) then
         !loop over lookup table points to find nearest values
         n1=1
         do while (dp .gt. dpdat(n1+1))
            n1=n1+1
         enddo
         n2=n1+1
         VALUE=scdat(n1)+(scdat(n2)-scdat(n1))
     &             *(dp-dpdat(n1))/(dpdat(n2)-dpdat(n1))
      endif

C If particle diameter is out of bounds, return reasonable value
      if (dp .gt. dpdat(numpts)) VALUE=2.0
      if (dp .lt. dpdat(1))      VALUE=7.0e-2

      ! Return to calling subroutine
      END FUNCTION STRATSCAV

!-----------------------------------------------------------------------------

      FUNCTION WATERNACL( RHE ) RESULT( VALUE )
!
!*****************************************************************************
!  Function WATERNACL uses the current RH to calculate how much water is 
!  in equilibrium with the seasalt.  Aerosol water concentrations are assumed
!  to be in equilibrium at all times and the array of concentrations is 
!  updated accordingly.
!  WRITTEN BY Peter Adams, November 2001
!  Introduced to GEOS-CHEM by Win Trivitayanurak. 8/6/07
!
!     VARIABLE COMMENTS...
!     waternacl is the ratio of wet mass to dry mass of a particle.  Instead
!     of calling a thermodynamic equilibrium code, this routine uses a
!     simple curve fit to estimate waternacl based on the current humidity.
!     The curve fit is based on ISORROPIA results for sodium sulfate
!     at 273 K.
!
!  NOTES:
!*****************************************************************************
!
      ! Arguments
      REAL*8                  :: RHE ! Relative humidity (0-100 scale)

      ! Return value
      REAL*8                  :: VALUE

      !=================================================================
      ! WATERNACL begins here!
      !=================================================================

      if (rhe .gt. 99.) rhe=99.
      if (rhe .lt. 1.) rhe=1.

         if (rhe .gt. 90.) then
            VALUE=5.1667642e-2*rhe**3-14.153121*rhe**2
     &               +1292.8377*rhe-3.9373536e4
         else
         if (rhe .gt. 80.) then
            VALUE=
     &      1.0629e-3*rhe**3-0.25281*rhe**2+20.171*rhe-5.3558e2
         else
         if (rhe .gt. 50.) then
            VALUE=
     &      4.2967e-5*rhe**3-7.3654e-3*rhe**2+.46312*rhe-7.5731
         else
         if (rhe .gt. 20.) then
            VALUE=
     &      2.9443e-5*rhe**3-2.4739e-3*rhe**2+7.3430e-2*rhe+1.3727
         else
            VALUE=1.17
         endif
         endif
         endif
         endif

         !check for error
         if (VALUE .gt. 45.) then
            write(*,*) 'ERROR in waternacl'
            write(*,*) rhe,VALUE
            STOP
         endif

      ! Return to calling subroutine
      END FUNCTION WATERNACL

!-----------------------------------------------------------------------------

      FUNCTION WATEROCIL( RHE ) RESULT( VALUE )
!
!*****************************************************************************
!  Function WATEROCIL uses the current RH to calculate how much water is 
!  in equilibrium with the hydrophillic OA.  Aerosol water concentrations
!  are assumed to be in equilibrium at all times and the array of 
!  concentrations is updated accordingly.
!  MODIFIED BY YUNHA LEE, AUG, 2006
!  Bring to GEOS-CHEM by Win Trivitayanurak 9/3/07
!
!     waterocil is the ratio of wet mass to dry mass of a particle.  Instead
!     of calling a thermodynamic equilibrium code, this routine uses a
!     simple curve fit to estimate waterocil based on the current humidity.
!     The curve fit is based on observations of Dick et al. JGR D1 1471-1479
!
!  NOTES:
!*****************************************************************************
!
      ! Arguments
      REAL*8                  :: RHE ! Relative humidity (0-100 scale)

      ! Return value
      REAL*8                  :: VALUE

      ! Local variables
      REAL*8                  :: a, b, c, d, e, f, prefactor, activcoef
      parameter(a=1.0034, b=0.1614, c=1.1693,d=-3.1,e=6.0)

      !=================================================================
      ! WATEROCIL begins here!
      !=================================================================

      if (rhe .gt. 99.) rhe=99.
      if (rhe .lt. 1.) rhe=1.

      if (rhe .gt. 85.) then
         VALUE =d+e*(rhe/100) 
cyhl Growth factor above RH 85% is not available, so it assumes linear growth 
cyhl at above 85%.  
      else
         VALUE =a+b*(rhe/100)+c*(rhe/100)**2. 
cyhl This eq is based on the extrapolation curve obtained from  
cyhl Dick et al 2000 figure 5.(High organic,density=1400g/cm3)
      endif
      
         !check for error
      if (VALUE .gt. 10.) then
         write(*,*) 'ERROR in waterocil'
         write(*,*) rhe, value
         STOP
      endif

      RETURN
      END FUNCTION WATEROCIL

!-----------------------------------------------------------------------------

      FUNCTION WATERSO4( RHE ) RESULT( VALUE )
!
!*****************************************************************************
!  Function WATERSO4 uses the current RH to calculate how much water is in 
!  equilibrium with the sulfate.  Aerosol water concentrations are assumed to 
!  be in equilibrium at all times and the array of concentrations is updated 
!  accordingly.
!     Introduced to GEOS-CHEM by Win Trivitayanurak. 8/6/07
!     Adaptation of ezwatereqm used in size-resolved sulfate only sim
!     November, 2001
!     ezwatereqm WRITTEN BY Peter Adams, March 2000
!
!     VARIABLE COMMENTS...
!
!     waterso4 is the ratio of wet mass to dry mass of a particle.  Instead
!     of calling a thermodynamic equilibrium code, this routine uses a
!     simple curve fit to estimate wr based on the current humidity.
!     The curve fit is based on ISORROPIA results for ammonium bisulfate
!     at 273 K.
!
!  NOTES:
!*****************************************************************************
!
      ! Arguments
      REAL*8                  :: RHE ! Relative humidity (0-100 scale)

      ! Return value
      REAL*8                  :: VALUE

      !=================================================================
      ! WATERSO4 begins here!
      !=================================================================

      if (rhe .gt. 99.) rhe=99.
      if (rhe .lt. 1.) rhe=1.

      if (rhe .gt. 96.) then
         value=
     &        0.7540688*rhe**3-218.5647*rhe**2+21118.19*rhe-6.801999e5
      else
         if (rhe .gt. 91.) then
            value=8.517e-2*rhe**2 -15.388*rhe +698.25
         else
            if (rhe .gt. 81.) then
               value=8.2696e-3*rhe**2 -1.3076*rhe +53.697
            else
               if (rhe .gt. 61.) then
                  value=9.3562e-4*rhe**2 -0.10427*rhe +4.3155
               else
                  if (rhe .gt. 41.) then
                     value=1.9149e-4*rhe**2 -8.8619e-3*rhe +1.2535
                  else
                     value=5.1337e-5*rhe**2 +2.6266e-3*rhe +1.0149
                  endif
               endif
            endif
         endif
      endif

      !check for error
      if (value .gt. 30.) then
         write(*,*) 'ERROR in waterso4'
         write(*,*) rhe,value
         STOP
      endif

      ! Return to calling subroutine
      END FUNCTION WATERSO4

!-----------------------------------------------------------------------------

      SUBROUTINE CLEANUP_TOMAS
!
!******************************************************************************
!  Subroutine CLEANUP_TOMAS deallocates all module arrays 
!  (win, 7/9/07)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_TRACER begins here!
      !=================================================================
      IF ( ALLOCATED( MOLWT       ) ) DEALLOCATE( MOLWT )

      ! Return to calling program
      END SUBROUTINE CLEANUP_TOMAS 

!-----------------------------------------------------------------------------

      ! End of module
      END MODULE TOMAS_MOD
#endif
