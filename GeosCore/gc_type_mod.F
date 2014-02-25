!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gc_type_mod.f
!
! !DESCRIPTION: Module GC\_TYPE\_MOD contains derived type definitions for 
!  GEOS-Chem.  These definitions are used to create objects for:
!  \begin{itemize}
!  \item ID flags for chemical species
!  \item ID flags for advected tracers
!  \item Coefficients \& other quantities that translate between chemical
!        species and advected tracers
!  \item ID flags for advected tracers that dry deposit
!  \item Logical flags for selecting the various GEOS-Chem options
!  \item GEOS-Chem columnized meteorological fields and related quantities
!  \item Error traceback output
!  \end{itemize}
!
! !INTERFACE:
!
      MODULE GC_TYPE_MOD
!
! USES:
!
      IMPLICIT NONE
      PUBLIC
!
! !PUBLIC TYPES:
! 
      !-----------------------------------------------------------------
      ! ID flags for advected tracers
      ! NOTE: This is a holdover from the column code and is presently
      ! only used in convection_mod.F (bmy, 10/22/12)
      !-----------------------------------------------------------------
      TYPE :: ID_TRAC
         INTEGER          :: NOx,     Ox,      PAN,     CO,      ALK4  
         INTEGER          :: ISOP,    HNO3,    H2O2,    ACET,    MEK   
         INTEGER          :: ALD2,    RCHO,    MVK,     MACR,    PMN
         INTEGER          :: PPN,     R4N2,    PRPE,    C3H8,    CH2O
         INTEGER          :: C2H6,    N2O5,    HNO4,    MP,      DMS
         INTEGER          :: SO2,     SO4,     SO4s,    MSA,     NH3
         INTEGER          :: NH4,     NIT,     NITs,    BCPI,    BCPO
         INTEGER          :: OCPI,    OCPO,    ALPH,    LIMO,    ALCO
         INTEGER          :: SOG1,    SOG2,    SOG3,    SOG4,    SOA1
         INTEGER          :: SOA2,    SOA3,    SOA4,    DST1,    DST2
         INTEGER          :: DST3,    DST4,    SALA,    SALC
         INTEGER          :: Hg0,     Hg2,     HgP
      END TYPE ID_TRAC

      !-----------------------------------------------------------------
      ! Coefficients & arrays that link species & tracers
      !-----------------------------------------------------------------
      TYPE :: SPEC_2_TRAC
         REAL*8,  POINTER :: SPEC_COEF(:,:)
         INTEGER, POINTER :: SPEC_ID(:,:)
         INTEGER, POINTER :: SPEC_EMITTED(:)
         INTEGER, POINTER :: SPEC_PER_TRAC(:)
         REAL*8,  POINTER :: TRAC_COEF(:)
         REAL*8,  POINTER :: MOLWT_KG(:)
         REAL*8,  POINTER :: XNUMOL(:)
      END TYPE SPEC_2_TRAC

      !-----------------------------------------------------------------
      ! Logical flags that turn various options on/off
      !-----------------------------------------------------------------
      TYPE :: GC_OPTIONS
         LOGICAL          :: USE_ANTHRO
         LOGICAL          :: USE_ANTHRO_BRAVO  
         LOGICAL          :: USE_ANTHRO_CAC       
         LOGICAL          :: USE_ANTHRO_EDGAR     
         LOGICAL          :: USE_ANTHRO_EPA 
         LOGICAL          :: USE_ANTHRO_VISTAS    
         LOGICAL          :: USE_ANTHRO_EMEP      
         LOGICAL          :: USE_BIOGENIC
         LOGICAL          :: USE_BIOMASS   
         LOGICAL          :: USE_BIOMASS_GFED2
         LOGICAL          :: USE_CARBON_AEROSOLS
         LOGICAL          :: USE_CHEMISTRY
         LOGICAL          :: USE_CONVECTION
         LOGICAL          :: USE_DEAD_DUST
         LOGICAL          :: USE_DEBUG_PRINT
         LOGICAL          :: USE_DRYDEP
         LOGICAL          :: USE_DUST_AEROSOLS
         LOGICAL          :: USE_EMISSIONS
         LOGICAL          :: USE_NOx_AIRCRAFT
         LOGICAL          :: USE_NOx_LIGHTNING
         LOGICAL          :: USE_NOx_SOIL
         LOGICAL          :: USE_PBL_MIXING
         LOGICAL          :: USE_SEC_ORG_AEROSOLS
         LOGICAL          :: USE_SHIP_ARCTAS   
         LOGICAL          :: USE_SEASALT_AEROSOLS
         LOGICAL          :: USE_SULFATE_AEROSOLS
         LOGICAL          :: USE_WETDEP
         LOGICAL          :: USE_Hg
         LOGICAL          :: USE_Hg_DYNOCEAN
         LOGICAL          :: USE_DIAG14
         LOGICAL          :: USE_DIAG38
      END TYPE GC_OPTIONS

      !-----------------------------------------------------------------
      ! Time & date values
      !-----------------------------------------------------------------
      TYPE :: GC_TIME
         INTEGER          :: YEAR         ! Current year (YYYY)
         INTEGER          :: MONTH        ! Current month (1-12)
         INTEGER          :: DAY          ! Current day (1-31)
         INTEGER          :: DOY          ! Day of year (0-365/366)
         INTEGER          :: HOUR         ! Current hour (0-23)
         INTEGER          :: MINUTE       ! Current minute (0-59)
         LOGICAL          :: FIRST_TIME   ! Is it the first timestep?
         REAL*8           :: T_ELAPSED    ! Elapsed simulation time [min]
         REAL*8           :: TS_DYN       ! Dynamic timestep [min]
         REAL*8           :: TS_CHEM      ! Chemistry timestep [min]
      END TYPE GC_TIME

      !-----------------------------------------------------------------
      ! Geographic location
      !-----------------------------------------------------------------
      TYPE :: GC_GEOLOC
         REAL*8           :: LON          ! Longitude [degrees]
         REAL*8           :: LAT          ! Latitude  [degrees]
         REAL*8           :: LOCALTIME    ! Local solar time [hrs]
      END TYPE GC_GEOLOC

      !-----------------------------------------------------------------
      ! Dimension information
      !-----------------------------------------------------------------
      TYPE :: GC_DIMS  
         INTEGER          :: L_COLUMN     ! # of boxes in the vertical column
         INTEGER          :: J_LAT        ! # of lat boxes
         INTEGER          :: I_LON        ! # of lon boxes
         INTEGER          :: N_AER        ! # of aerosol tracers
         INTEGER          :: N_DRYDEP     ! # of dry deposited tracers
         INTEGER          :: N_DUST       ! # of dust tracers
         INTEGER          :: N_JV         ! # of J-value reactions
         INTEGER          :: N_MEMBERS    ! Max # of species per family tracer
         INTEGER          :: N_RH         ! # of RH bins for photolysis
         INTEGER          :: N_REACTIONS  ! # of chemical reactons
         INTEGER          :: N_SOA_HC     ! # of SOA HC classes and products
         INTEGER          :: N_SOA_PROD   !   (dimensions of GPROD/APROD)
         INTEGER          :: N_SPECIES    ! # of chemical species
         INTEGER          :: N_TRACERS    ! # of advected tracers
         INTEGER          :: N_WETDEP     ! # of wet deposited tracers
      END TYPE GC_DIMS

      !-----------------------------------------------------------------
      ! Derived type to pass ID information to GEOS-Chem routines
      !-----------------------------------------------------------------
      TYPE :: GC_IDENT
         INTEGER            :: PET          ! # of the CPU we are on
         INTEGER            :: STDOUT_LUN   ! LUN for stdout redirect
         CHARACTER(LEN=255) :: STDOUT_FILE  ! Filename for stdout redirect
         INTEGER            :: LEV          ! Stack pointer
         CHARACTER(LEN=40)  :: I_AM(20)     ! Stack for routine names
         CHARACTER(LEN=999) :: ERRMSG       ! Error message to display
         LOGICAL            :: VERBOSE      ! Should we print out debug info?
      END TYPE GC_IDENT
!
! !REVISION HISTORY: 
!  24 Mar 2009 - R. Yantosca - Initial version
!  21 Apr 2009 - R. Yantosca - Renamed from "id_type_mod.f" to
!                              "gc_type_mod.f"; added type GC_OPTIONS
!  05 Jun 2009 - R. Yantosca - added LISOPOH to ID_SPEC
!                            - added MOLWT, XNUMOL to SPEC_2_TRAC
!  08 Jul 2009 - R. Yantosca - Added USE_CONVECTION, USE_PBL_MIXING,
!                              USE_WETDEP, USE_EMISSIONS flags to GC_OPTIONS
!  24 Aug 2009 - R. Yantosca - Added GC_TIME type
!  30 Oct 2008 - R. Yantosca - Added GC_GEOLOC type
!  05 Nov 2009 - R. Yantosca - Added GC_MET_LOCAL type
!  13 Apr 2010 - R. Yantosca - Added N_MEMBERS to GC_DIMS type
!  16 Apr 2010 - R. Yantosca - Added all chemical species to ID_SPEC
!  16 Apr 2010 - R. Yantosca - Added N_SPECIES, N_JV to GC_DIMS
!  23 Apr 2010 - R. Yantosca - Added GC_IDENT object
!  29 Apr 2010 - R. Yantosca - Added TO3, FRCLND to GC_MET_LOCAL object
!  30 Apr 2010 - R. Yantosca - Change character lengths in GC_IDENT
!  13 May 2010 - R. Yantosca - Updated comments  
!  17 May 2010 - R. Yantosca - Added L_COLUMN to DIMINFO
!  02 Jun 2010 - R. Yantosca - Added VERBOSE to GC_IDENT
!  02 Jun 2010 - R. Yantosca - Added TAUCLI, TAUCLW to GC_MET_LOCAL
!  02 Jun 2010 - R. Yantosca - Added DQ*DTMST fields to GC_MET_LOCAL
!  02 Jun 2010 - R. Yantosca - Moved type SCOX_1d here from schem_mod.F
!  22 Oct 2012 - R. Yantosca - Comment out obsolete type definitions
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
      END MODULE GC_TYPE_MOD
!EOC
