#include "MAPL_Generic.h"

!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GEOSCHEMchem_GridCompMod
!
! !DESCRIPTION: GEOSCHEMchem_GridComp is an ESMF gridded component 
!  implementing the GEOS-Chem chemical processes.  
!\\
!\\
! !INTERFACE:
!
MODULE GEOSCHEMchem_GridCompMod
!
! !USES:
!
  USE ESMF_Mod                                       ! ESMF library
  USE MAPL_Mod                                       ! MAPL library
  USE GC_Type_Mod                                    ! Derived type definitions
  USE GC_Column_Mod                                  ! GEOS-Chem column code
  USE Schem_Mod, ONLY : SCOX_1D                      ! Typedef for strat Ox

  IMPLICIT NONE
  PRIVATE

# include "smv_dimension.h"                          ! Max levels & tracers
# include "smv_errcode.h"                            ! GEOS-Chem error codes
# include "smv_physconst.h"                          ! Physical constants
# include "lai_land_info.h"                          ! Sizes for land/LAI data
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC                            :: SetServices   ! Sets ESMF entry points
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE                           :: Initialize_      ! Init method
  PRIVATE                           :: Run_             ! Run method  
  PRIVATE                           :: Finalize_        ! Finalize method
  PRIVATE                           :: Extract_         ! Get values from ESMF
  PRIVATE                           :: Setup_GeoLoc_    ! Init GEOLOC object
  PRIVATE                           :: Setup_Options_   ! Init OPTIONS object
  PRIVATE                           :: Setup_DimInfo_   ! Init DIMINFO object
  PRIVATE                           :: Setup_TrcName_   ! Init TRCNAME object
  PRIVATE                           :: Error_Trap_      ! Handle error state
  PRIVATE                           :: Roundoff         ! Truncates a number
  PRIVATE                           :: Print_Mean_OH    ! Mean OH lifetime
  PRIVATE                           :: GlobalSum        ! Sums across CPUs
!
! !PRIVATE TYPES:
!
  ! Legacy state
  TYPE GEOSCHEM_State
     PRIVATE
     TYPE(ESMF_Config)              :: myCF          ! Private ESMF Config obj
  END TYPE GEOSCHEM_State

  ! Hook for the ESMF
  TYPE GEOSCHEM_Wrap
     TYPE (GEOSCHEM_State), POINTER :: PTR => null() ! Ptr to GEOSCHEM_State
  END TYPE GEOSCHEM_Wrap

  ! Objects for GEOS-Chem column chemistry code
  TYPE(GC_Dims)                     :: DimInfo       ! Dimension information
  TYPE(GC_Options)                  :: Options       ! Logical flags
  TYPE(ID_Dryd)                     :: IdDrydep      ! ID's for drydep species
  TYPE(ID_Trac)                     :: IdTracers     ! ID's for tracers
  TYPE(ID_Spec)                     :: IdSpecies     ! ID's for chem species
  TYPE(ID_Wetd)                     :: IdWetdep      ! ID's for wetdep species
  TYPE(Spec_2_Trac)                 :: Coef          ! Species <-> tracer map  

  ! Scalars
  INTEGER                           :: logLun        ! LUN for stdout logfile
!
! !DEFINED PARAMETERS:
!
  ! Scale factor to prevent underflow in mean OH diagnostic
  REAL*8, PARAMETER                 :: OH_SCALE = 1d20
!
! !REMARKS:
!  Developed for GEOS-5 release Fortuna 2.0 and later.
!
! !REVISION HISTORY:
!  06 Dec 2009 - A. da Silva - Created the GEOSCHEM skeleton.
!  07 Apr 2010 - R. Yantosca - Updated comments, cosmetic changes 
!  08 Apr 2010 - R. Yantosca - Updated the Extract_ method
!  13 Apr 2010 - R. Yantosca - Now call down to the GC_COLUMN_INIT method
!  02 Jul 2010 - R. Yantosca - Add Print_Mean_OH, GlobalSum methods
!  02 Jul 2010 - R. Yantosca - Add OH_SCALE parameter for mean OH lifetime
!  14 Jul 2010 - R. Yantosca - Force lons & lats to 4 places of precision
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices
!
! !DESCRIPTION: The SetServices routine does the following:
!
! \begin{itemize}
! \item Defines the Initialize method for the GEOSCHEMchem gridded component
! \item Defines the Run method for the GEOSCHEMchem gridded component
! \item Defines the Finalize method for the GEOSCHEMchem gridded component
! \item Attaches an internal state (which holds a private ESMF Config object)
!       to the GEOSCHEMchem gridded component.
! \end{itemize}
!
! !INTERFACE:
!
  SUBROUTINE SetServices( GC, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp),  INTENT(INOUT) :: GC       ! Ref to this GridComp
!
! !OUTPUT PARAMETERS:
!
    INTEGER,              INTENT(OUT)   :: RC         ! Return code
!
! !REMARKS:
!  ESMF can only attach one Config object per Gridded Component.  The 
!  Config object that is defined from the "MAPL.rc" resource file is 
!  directly attached to the GEOSCHEMchem gridded component.  
!                                                                             .
!  To attach the Config object defined from the "GEOSCHEMchem_GridComp.rc" 
!  resource file, we must first create a derived type with a pointer to
!  the Config object, then attach that to the gridded component as an
!  "internal state" (also called "legacy state").
!
! !REVISION HISTORY:
!  06 Dec 2009 - A. da Silva - Initial version
!  07 Apr 2010 - R. Yantosca - Updated comments, cosmetic changes 
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    TYPE(GEOSCHEM_State), POINTER       :: myState    ! Legacy state
    TYPE(GEOSCHEM_Wrap)                 :: wrap       ! Wrapper for myState
    CHARACTER(LEN=ESMF_MAXSTR)          :: compName   ! Gridded Component name

    !=======================================================================
    ! Initialization
    !=======================================================================

    __Iam__('SetServices')

    ! Get my name and set-up traceback handle
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )
    Iam = TRIM( compName ) // '::' // TRIM( Iam )

    !=======================================================================
    ! Wrap internal state for storing in this gridded component
    ! Rename this to a "legacy state"
    !=======================================================================
    ALLOCATE( myState, stat=STATUS )
    VERIFY_(STATUS)
    wrap%ptr => myState

    !=======================================================================
    ! Define an ESMF Config object from the Resource file and set it 
    ! as an "internal state" of the GEOSCHEMchem gridded component
    !=======================================================================
    myState%myCF = ESMF_ConfigCreate(__RC__)
    call ESMF_ConfigLoadFile( myState%myCF, 'GEOSCHEMchem_GridComp.rc', __RC__)

    !=======================================================================
    !                 %%% ESMF Functional Services %%%
    !=======================================================================

    ! Set the Initialize, Run, Finalize entry points
    CALL MAPL_GridCompSetEntryPoint( GC, ESMF_SETINIT,  Initialize_, __RC__ )
    CALL MAPL_GridCompSetEntryPoint( GC, ESMF_SETRUN,   Run_,        __RC__ )
    CALL MAPL_GridCompSetEntryPoint( GC, ESMF_SETFINAL, Finalize_,   __RC__ )
        
    ! Store internal state with Config object in the gridded component
    CALL ESMF_UserCompSetInternalState( GC, 'GEOSCHEM_State', wrap, STATUS )
    VERIFY_(STATUS)
  
    !=======================================================================
    !                    %%% MAPL Data Services %%%
    !=======================================================================
!EOC
!BOP
!
! !IMPORT STATE:
!
#   include "GEOSCHEMchem_ImportSpec___.h"
!
! !INTERNAL STATE:
!
#   include "GEOSCHEMchem_InternalSpec___.h"
!
! !EXTERNAL STATE:
!
#   include "GEOSCHEMchem_ExportSpec___.h"
!EOP
!BOC

    ! Generic Set Services
    CALL MAPL_GenericSetServices( GC, __RC__ )

    !=======================================================================
    ! All done
    !=======================================================================
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE SetServices
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize_ 
!
! !DESCRIPTION: Initialize_ is the initialize method of the GEOSCHEMchem
!  gridded component.  This is a simple ESMF/MAPL wrapper which calls down
!  to the Initialize method of the GEOS-Chem column chemistry code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Initialize_( GC, Import, Export, Clock, RC )
!
! !USES:
! 
    USE Charpak_Mod, ONLY : StrRepl
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC          ! Ref. to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT) :: Import      ! Import state
    TYPE(ESMF_State),    INTENT(INOUT) :: Export      ! Export state
    TYPE(ESMF_Clock),    INTENT(INOUT) :: Clock       ! ESMF clock object
!                                                      
! !OUTPUT PARAMETERS:                                  
!                                                      
    INTEGER,             INTENT(OUT)   :: RC          ! Error return code
!
! !REVISION HISTORY:
!  06 Dec 2009 - A. da Silva - Initial version
!  08 Apr 2010 - R. Yantosca - Now uses the updated Extract_ method.
!  14 Apr 2010 - R. Yantosca - Activated call to GC_COLUMN_INIT
!  15 Apr 2010 - R. Yantosca - Add extra error checks for dimensions
!  23 Apr 2010 - R. Yantosca - Now pass IDENT obj to GC_COLUMN_INIT routine
!  30 Apr 2010 - R. Yantosca - Now use 5 digits for PET
!  02 Jun 2010 - R. Yantosca - Now set Ident%VERBOSE to FALSE
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Objects
    TYPE(ESMF_Grid)            :: Grid        ! ESMF Grid object
    TYPE(ESMF_Config)          :: MaplCF      ! ESMF Config obj (MAPL.rc)
    TYPE(ESMF_Config)          :: GeosCF      ! ESMF Config obj (GEOSCHEM*.rc) 
    TYPE(GC_Ident)             :: Ident       ! ID information 
                                      
    ! Scalars                                  
    INTEGER                    :: error       ! GEOS-Chem error code
    INTEGER                    :: myPet       ! # of the PET we are on 
    INTEGER                    :: LM          ! # of vertical levels
    REAL                       :: tsChem      ! Chemistry timestep [s]
    REAL                       :: tsDyn       ! Dynamic timestep [s]
    CHARACTER(LEN=5)           :: petStr      ! String for PET #
    CHARACTER(LEN=ESMF_MAXSTR) :: compName    ! Name of gridded component
    CHARACTER(LEN=ESMF_MAXSTR) :: logFile     ! File for stdout redirect
     
    ! Arrays
    CHARACTER(LEN=ESMF_MAXSTR) :: trcName(MAX_TRACERS) ! Tracer names

    !=======================================================================
    ! Initialization
    !=======================================================================

    __Iam__('Initialize_')

    ! Traceback info
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )
    Iam = trim( compName ) // '::' // trim( Iam )

    ! Initialize MAPL Generic
    CALL MAPL_GenericInitialize( GC, Import, Export, Clock, __RC__ )

    ! Get various parameters from the ESMF/MAPL framework
    CALL Extract_( GC, Clock,                   &
                   Grid     = Grid,             &  ! ESMF Grid object
                   MaplCF   = MaplCF,           &  ! MAPL.rc Config object
                   GeosCF   = GeosCF,           &  ! GEOSCHEM*.rc Config object
                   LM       = DimInfo%L_COLUMN, &  ! # of vertical levels
                   localPet = myPet,            &  ! PET # we are on now 
                   tsChem   = tsChem,           &  ! Chemistry timestep
                   __RC__ )

    ! Define variables for passing down to GEOS-Chem column chemistry code
    CALL Setup_DimInfo_( compName, GeosCF, DimInfo, __RC__ )
    CALL Setup_Options_( compName, GeosCF, Options, __RC__ )
    CALL Setup_TrcName_( compName, GeosCF, trcName, __RC__ )

    ! Stop if we have the wrong # of dimensions
    ASSERT_( DimInfo%L_COLUMN   == MAX_COLUMN  )   ! # of vertical levels
    ASSERT_( DimInfo%N_TRACERS  >  0           )   ! # of tracers: lower bound
    ASSERT_( DimInfo%N_TRACERS  <= MAX_TRACERS )   ! # of tracers: upper bound
    ASSERT_( N_OLSON_LOCAL      == 15          )   ! # of lai/leaf types/box
    ASSERT_( DimInfo%N_AER      >  0           )   ! # of aerosol bins
    ASSERT_( DimInfo%N_DUST     >  0           )   ! # of dust bins
    ASSERT_( DimInfo%N_MEMBERS  >  0           )   ! # of species/tracer
    ASSERT_( DimInfo%N_RH       >  0           )   ! # of RH bins
    ASSERT_( DimInfo%N_SOA_HC   >  0           )   ! # of GPROD/APROD HC types
    ASSERT_( DimInfo%N_SOA_PROD >  0           )   ! # of GPROD/APROD products
    
    ! Name of logfile for stdout redirect
    CALL ESMF_ConfigGetAttribute( GeosCF, Ident%STDOUT_FILE,   &
                                  Label   = "STDOUT_LOGFILE:", &
                                  Default = "PET%%%%%.init",    &
                                   __RC__ )

    ! Name of log LUN # for stdout redirect
    CALL ESMF_ConfigGetAttribute( GeosCF, logLun,              &
                                  Label   = "STDOUT_LOGLUN:",  &
                                  Default = 700,               &
                                   __RC__ )

    ! Fill the remaining fields of the IDENT object
    Ident%STDOUT_LUN = logLun
    Ident%PET        = myPet
    Ident%I_AM(1)    = TRIM( compName )
    Ident%I_AM(2)    = 'Initialize_'
    Ident%LEV        = 2
    Ident%ERRMSG     = ''
    Ident%VERBOSE    = .FALSE.

    !=======================================================================
    ! Open a log file on each PET where stdout will be redirected
    !=======================================================================

    ! Replace tokens w/ PET # in the filename
    logFile = Ident%STDOUT_FILE
    WRITE( petStr, '(i5.5)' ) myPet
    CALL StrRepl( logFile, '%%%%%', petStr )

    ! Open file for stdout redirect
    OPEN ( logLun, FILE=LOGFILE, STATUS='UNKNOWN' )

    ! Add descriptive header text
    WRITE( logLun, '(a)'   ) REPEAT( '#', 79 )
    WRITE( logLun, 100     ) TRIM( logFile ), TRIM( Iam ), myPet
    WRITE( logLun, '(a,/)' ) REPEAT( '#', 79 )

    !=======================================================================
    ! Initialize GEOS-Chem
    !=======================================================================

    ! Call the Initialize routine of the GEOS-Chem column chemistry
    CALL GC_Column_Init( IDENT       = Ident,                        &
                         OPTIONS     = Options,                      &
                         L_COLUMN    = DimInfo%L_COLUMN,             &
                         N_TRACERS   = DimInfo%N_TRACERS,            &   
                         N_MEMBERS   = DimInfo%N_MEMBERS,            &
                         TS_CHEM     = DBLE( tsChem ),               &
                         TRACER_NAME = trcName(1:DimInfo%N_TRACERS), &
                         COEF        = Coef,                         &
                         ID_SPECIES  = IdSpecies,                    & 
                         ID_TRACERS  = IdTracers,                    &
                         ID_DRYDEP   = IdDrydep,                     &
                         ID_WETDEP   = IdWetdep,                     &
                         N_DRYDEP    = DimInfo%N_DRYDEP,             &
                         N_WETDEP    = DimInfo%N_WETDEP,             &    
                         N_SPECIES   = DimInfo%N_SPECIES,            &
                         N_REACTIONS = DimInfo%N_REACTIONS,          &
                         N_JV        = DimInfo%N_JV,                 &
                         RC          = error )

    ! Trap the error from GEOS-Chem
    IF ( error /= SMV_SUCCESS ) THEN 
       CALL Error_Trap_( Ident, error, __RC__ )
    ENDIF

    ! Check # of drydep species
    IF ( Options%USE_DRYDEP ) THEN 
       ASSERT_( DimInfo%N_DRYDEP > 0 ) 
    ENDIF

    ! Check # of wetdep species
    IF ( Options%USE_WETDEP ) THEN
       ASSERT_( DimInfo%N_WETDEP > 0 )
    ENDIF

    ! Check # of chemical rxns & species
    IF ( Options%USE_CHEMISTRY ) THEN
       ASSERT_( DimInfo%N_JV        >  0           ) 
       ASSERT_( DimInfo%N_REACTIONS >  0           ) 
       ASSERT_( DimInfo%N_SPECIES   >  0           )
       ASSERT_( DimInfo%N_SPECIES   <= MAX_SPECIES )
    ENDIF

    !=======================================================================
    ! All done
    !=======================================================================

    ! Write a header before the timestepping begins
    WRITE( logLun, '(/,a)' ) REPEAT( '#', 79 )
    WRITE( logLun, 200     ) TRIM( compName ) // '::Run_', myPet
    WRITE( logLun, '(a,/)' ) REPEAT( '#', 79 )

    ! Successful return
    RETURN_(ESMF_SUCCESS)

    ! Formats
100 FORMAT( '### ',                                           / &
            '### ', a ,                                       / &
            '### ', a, '  |  Initialization on PET # ', i5.5, / &
            '### ' )
200 FORMAT( '### ',                                           / &
            '### ', a, '  |  Execution on PET # ',      i5.5, / &
            '###' )

  END SUBROUTINE Initialize_
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Run_
!
! !DESCRIPTION: Run_ is the run method of the GEOSCHEMchem gridded component.  
!  GC is a simple ESMF/MAPL wrapper which calls down to the Run method of 
!  the GEOS-Chem column chemistry code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Run_( GC, Import, Export, Clock, RC )
!
! !USES:
!
#   include "GEOSCHEMchem_DeclarePointer___.h"        ! Ptr decls to states
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC          ! Ref to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT) :: Import      ! Import State
    TYPE(ESMF_State),    INTENT(INOUT) :: Export      ! Export State
    TYPE(ESMF_Clock),    INTENT(INOUT) :: Clock       ! ESMF Clock object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC          ! Error return code
!
! !REVISION HISTORY:
!  06 Dec 2009 - A. da Silva - Initial version
!  08 Apr 2010 - R. Yantosca - Now uses the updated Extract_ method
!  09 Apr 2010 - R. Yantosca - Initialize Timing, GeoLoc objects
!  16 Apr 2010 - R. Yantosca - Now move the array assignments before & after
!                              the call to GC_COLUMN_RUN into separate
!                              include files, for clarity
!  30 Apr 2010 - R. Yantosca - Now use 5 digits for PET
!  02 Jun 2010 - R. Yantosca - Now use IDENT%VERBOSE to trigger debug output
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!  
    ! Objects
    TYPE(ESMF_Grid)              :: Grid          ! ESMF Grid object
    TYPE(ESMF_Config)            :: MaplCF        ! Config (MAPL.rc)
    TYPE(ESMF_Config)            :: GeosCF        ! Config (GEOSCHEM*.rc)
    TYPE(GC_Time)                :: Timing        ! G-C obj for time values
    TYPE(GC_GeoLoc)              :: GeoLoc        ! G-C obj for location
    TYPE(GC_Ident)               :: Ident         ! G-C obj for ident info
    TYPE(GC_Met_1d)              :: Met           ! G-C obj for met fields
    TYPE(SCOX_1d)                :: Schem         ! G-C obj for strat Ox
                                                  
    ! Scalars                                     
    INTEGER                      :: IM            ! # of lons   on this PET
    INTEGER                      :: JM            ! # of lats   on this PET
    INTEGER                      :: LM            ! # of levels on this PET
    INTEGER                      :: error         ! G-C error return code
    INTEGER(ESMF_KIND_I8)        :: advCount      ! # of clock advances
    INTEGER                      :: nymd          ! YYYY/MM/DD date
    INTEGER                      :: nhms          ! hh:mm:ss time
    INTEGER                      :: myPet         ! PET # we are on now
    INTEGER                      :: nPets         ! Total # of PETs
    INTEGER                      :: I, J, L       ! Loop indices
    INTEGER                      :: N_TRC         ! Shadow var: # of tracers
    INTEGER                      :: year          ! Current year    
    INTEGER                      :: month         ! Current month
    INTEGER                      :: day           ! Current day
    INTEGER                      :: dayOfYr       ! Current day of year
    INTEGER                      :: hour          ! Current hour
    INTEGER                      :: minute        ! Current minute
    REAL                         :: UTC           ! Universal time
    REAL                         :: tsChem        ! Chem timestep [min]
    REAL                         :: tsDyn         ! Dynamic timestep [min]
    REAL                         :: hElapsed      ! Elapsed time [hours]
    REAL*8                       :: lonDeg        ! Longitude [degrees]
    REAL*8                       :: latDeg        ! Latitude [degrees]
    REAL                         :: locTime       ! Local time [hours]
    REAL*8                       :: P1, P2        ! Pressure variables
    CHARACTER(LEN=4)             :: petStr        ! String for PET #
    CHARACTER(LEN=ESMF_MAXSTR)   :: logFile       ! File for stdout redirect
    CHARACTER(LEN=ESMF_MAXSTR)   :: compName      ! Gridded Component name
                                                  
    ! Arrays                                            
    REAL(ESMF_KIND_R4),  POINTER :: lonCtr(:,:)   ! Lon centers [radians]
    REAL(ESMF_KIND_R4),  POINTER :: latCtr(:,:)   ! Lat centers [radians]
    REAL(ESMF_KIND_R4),  POINTER :: latEdg(:,:)   ! Lat centers [radians]

    ! Local variables for passing values from States to GEOS-Chem code
    INTEGER                      :: iReg_1d
    INTEGER                      :: iLand_1d   ( N_OLSON_LOCAL           )
    INTEGER                      :: iUse_1d    ( N_OLSON_LOCAL           )
    REAL*8                       :: lai_1d     ( N_OLSON_LOCAL           )
    REAL*8,              SAVE    :: airMassDiag( MAX_COLUMN              )
    REAL*8,              SAVE    :: ohMassDiag ( MAX_COLUMN              )
    REAL*8                       :: h2o2s_1d   ( MAX_COLUMN              )
    REAL*8                       :: so2s_1d    ( MAX_COLUMN              )
    REAL*8                       :: terp_1d    ( MAX_COLUMN              )
    REAL*8                       :: sesq_1d    ( MAX_COLUMN              )
    REAL*8                       :: gProd_1d   ( MAX_COLUMN, 3, 6        )   
    REAL*8                       :: aProd_1d   ( MAX_COLUMN, 3, 6        )   
    REAL*8                       :: alkEmis_1d ( MAX_COLUMN, 2           )
    REAL*8                       :: nDens_1d   ( MAX_COLUMN, 2           )
    REAL*8                       :: cspec_1d   ( MAX_COLUMN, MAX_SPECIES )
    REAL*8                       :: emiss_1d   ( MAX_COLUMN, MAX_TRACERS ) 
    REAL*8                       :: tracer_1d  ( MAX_COLUMN, MAX_TRACERS ) 
    REAL*8                       :: wdLossDiag ( MAX_COLUMN, MAX_TRACERS )
    REAL*8                       :: ddVelDiag  (             MAX_TRACERS )
    REAL*8                       :: ddFluxDiag (             MAX_TRACERS )
    REAL*8                       :: ddFreqDiag (             MAX_TRACERS )
    REAL*8                       :: To3_1d
    REAL*8                       :: FrcLnd_1d
!
! !DEFINED PARAMETERS:
! 
    ! Hardwire these for now
    INTEGER,           PARAMETER :: idEmission(13) = (/ &
                                     1,4,18,19,5,21,9,10,11,20,2,7,6 /)
        
    !=======================================================================
    ! Initialization
    !=======================================================================

    __Iam__('Run_')

    ! Traceback info
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )
    Iam = TRIM( compName ) // '::' // TRIM( Iam )

    ! Get pointers to fields in import, internal, and export states
#   include "GEOSCHEMchem_GetPointer___.h"

    ! Get various parameters from the ESMF/MAPL framework
    CALL Extract_( GC, Clock,              &
                   GRID      = Grid,       &  ! ESMF Grid obj
                   MAPLCF    = MaplCF,     &  ! ESMF Config obj (MAPL*.rc) 
                   GEOSCF    = GeosCF,     &  ! ESMF Config obj (GEOSCHEM*.rc)
                   TSCHEM    = tsChem,     &  ! Chemistry timestep [min]
                   TSDYN     = tsDyn,      &  ! Dynamic timestep [min]
                   YEAR      = year,       &  ! Current year
                   MONTH     = month,      &  ! Current month
                   DAY       = day,        &  ! Current day
                   DAYOFYR   = dayOfYr,    &  ! Current day of year
                   HOUR      = hour,       &  ! Current hour
                   MINUTE    = minute,     &  ! Current minute
                   HELAPSED  = hElapsed,   &  ! Elapsed hours
                   ADVCOUNT  = advCount,   &  ! # of times clock has advanced
                   IM        = IM,         &  ! # of longitudes on this PET
                   JM        = JM,         &  ! # of latitudes  on this PET
                   LM        = LM,         &  ! # of levels     on this pET
                   UTC       = utc,        &  ! Universal time [hours]
                   LOCALPET  = myPet,      &  ! # of the PET we are on now
                   PETCOUNT  = nPets,      &  ! Total # of PETs
                   LONCTR    = lonCtr,     &  ! Array of lon centers [radians]
                   LATCTR    = latCtr,     &  ! Array of lat centers [radians]
                   LATEDG    = latEdg,     &  ! Array of lat edges   [radians]
                   __RC__ )

    ! Make sure we have proper dimension sizes
    ASSERT_( IM >  0                )
    ASSERT_( JM >  0                )
    ASSERT_( LM == MAX_COLUMN       )
    ASSERT_( LM == DimInfo%L_COLUMN )

    ! Populate the TIMING object for GEOS-Chem column code
    Timing%Year       = year     
    Timing%Month      = month    
    Timing%Day        = day
    Timing%Doy        = dayOfYr
    Timing%Hour       = hour
    Timing%Minute     = minute
    Timing%First_Time = ( advCount == 0 )
    Timing%T_Elapsed  = DBLE( hElapsed ) * 60d0  ! min
    Timing%Ts_Dyn     = DBLE( tsDyn    )         ! min
    Timing%Ts_Chem    = DBLE( tsChem   )         ! min

    !=======================================================================
    ! Print timing etc. info to the log file outside of the (I,J) loop
    !=======================================================================

    ! Write time quantities
    WRITE( logLun, 100 ) year, month, day, hour, minute, hElapsed

    !=======================================================================
    ! Loop over all grid boxes on this PET
    !=======================================================================
    DO J = 1, JM
    DO I = 1, IM

       ! Initialize fields of the IDENT object for each (I,J) location
       Ident%STDOUT_FILE = ''
       Ident%STDOUT_LUN  = logLun
       Ident%PET         = myPet
       Ident%I_AM(1)     = TRIM( compName )
       Ident%I_AM(2)     = 'Run_'
       Ident%LEV         = 2
       Ident%ERRMSG      = ''

       !--------------------------------------------------------------------
       ! Set up geographic coordinates for column chemistry code
       !--------------------------------------------------------------------
       CALL Setup_GeoLoc_( lonCtr(I,J), latCtr(I,J), UTC, GeoLoc )

!#############################################################################
!### DEBUG: Set flag for writing debug output if we are at our test box
       IF ( INT( GeoLoc%LON ) == -80 .and. INT( GeoLoc%LAT ) == 42 ) THEN
          Ident%VERBOSE = .TRUE.
       ELSE
          Ident%VERBOSE = .FALSE.
       ENDIF
!#############################################################################

       !-------------------------------------------------------------------
       ! Call GEOS-Chem column chemistry code
       !-------------------------------------------------------------------

       ! Include file containing pre-Run method array assignments
#      include "Includes_Before_Run.h"
       
       ! Shadow the # of tracers in a convenience variable
       N_TRC = DimInfo%N_TRACERS

       ! Run the GEOS-Chem column chemistry code
       CALL GC_Column_Run( COEF          = Coef,                      &
                           DIMINFO       = DimInfo,                   &
                           GEOLOC        = GeoLoc,                    &
                           ID_DRYDEP     = IdDrydep,                  &
                           ID_EMISSION   = idEmission,                &
                           ID_SPECIES    = IdSpecies,                 &
                           ID_TRACERS    = IdTracers,                 &
                           ID_WETDEP     = IdWetdep,                  &
                           IDENT         = Ident,                     &
                           MET_1d        = Met,                       &
                           OPTIONS       = Options,                   &
                           TIMING        = Timing,                    &
                           ALK_EMIS_1d   = alkEmis_1d,                &
                           N_DENS_1d     = nDens_1d,                  &
                           IREG_1d       = iReg_1d,                   &
                           ILAND_1d      = iLand_1d,                  &
                           IUSE_1d       = iUse_1d,                   &
                           LAI_1d        = lai_1d,                    &
                           ORVC_TERP     = terp_1d,                   & 
                           ORVC_SESQ     = sesq_1d,                   & 
                           GPROD         = gProd_1d,                  &
                           APROD         = aProd_1d,                  &
                           OXIDANTS_1d   = Schem,                     &
                           CSPEC_FULL_1d = cspec_1d,                  &
                           H2O2s_1d      = h2o2s_1d,                  & 
                           SO2s_1d       = so2s_1d,                   &
                           EMISSION_1d   = emiss_1d   ( :, 1:N_TRC ), &
                           TRACER_1d     = tracer_1d  ( :, 1:N_TRC ), &
                           DD_VEL_1d     = ddVelDiag  (    1:N_TRC ), &
                           DD_FREQ_1d    = ddFreqDiag (    1:N_TRC ), &
                           DD_FLUX_1d    = ddFluxDiag (    1:N_TRC ), &
                           WETD_LOSS_1d  = wdLossDiag ( :, 1:N_TRC ), &
                           AIR_MASS_1d   = airMassDiag,               & 
                           OH_MASS_1d    = ohMassDiag,                & 
                           RC            = error )

       ! Trap the error from GEOS-Chem
       IF ( error /= SMV_SUCCESS ) THEN 
          CALL Error_Trap_( Ident, error, __RC__ )
       ENDIF

       ! Include file w/ post-Run method array assignments
#      include "Includes_After_Run.h"

    ENDDO
    ENDDO 

    !=======================================================================
    ! All done
    !=======================================================================

    ! Free pointers
    NULLIFY( lonCtr, latCtr )

    ! Successful return
    RETURN_(ESMF_SUCCESS)

    ! Formats
100 FORMAT( '---> DATE: ', i4.4, '/', i2.2, '/', i2.2,      &
            '  GMT: ', i2.2, ':', i2.2, '  X-HRS: ', f11.3 )
110 FORMAT( 'Box (',i3,',',i3,') on PET ', i3, ' has coords: ', 2f7.2, &
               ' LocT = ', f9.4 )

  END SUBROUTINE Run_
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finalize_
!
! !DESCRIPTION: Finalize_ is the finalize method of the GEOSCHEMchem gridded 
!  component.  GC is a simple ESMF/MAPL wrapper which calls down to the
!  Finalize method of the GEOS-Chem column chemistry code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Finalize_( GC, Import, Export, Clock, RC )
!
! !USES:
!
#   include "GEOSCHEMchem_DeclarePointer___.h"       ! Ptr decls to states
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC         ! Ref. to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT) :: Import     ! Import State
    TYPE(ESMF_State),    INTENT(INOUT) :: Export     ! Export State
    TYPE(ESMF_Clock),    INTENT(INOUT) :: Clock      ! ESMF Clock object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC         !  Error return code
!
! !REVISION HISTORY:
!  01 Dec 2009 - A. Da Silva - Initial version
!  08 Apr 2010 - R. Yantosca - Now finalize myState%CF and myState
!  15 Apr 2010 - R. Yantosca - Activate call to GC_COLUMN_FINAL
!  30 Apr 2010 - R. Yantosca - Now use 5 digits for PET
!  02 Jun 2010 - R. Yantosca - Now set Ident%VERBOSE to FALSE
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Objects
    TYPE(ESMF_Grid)              :: Grid          ! ESMF Grid object
    TYPE(ESMF_Config)            :: MaplCF        ! Config (MAPL.rc)
    TYPE(ESMF_Config)            :: GeosCF        ! Config (GEOSCHEM*.rc)
    TYPE(GC_GeoLoc)              :: GeoLoc        ! G-C obj for location
    TYPE(GC_Ident)               :: Ident         ! G-C obj for ident info

    ! Scalars
    CHARACTER(LEN=ESMF_MAXSTR)   :: compName      ! Gridded component name
    INTEGER                      :: error         ! GEOS-Chem error code
    INTEGER                      :: myPet         ! # of PET we are on now
    INTEGER                      :: I,  J         ! Loop indices
    INTEGER                      :: IM, JM, LM    ! Loop bounds
    REAL                         :: UTC           ! UTC time [hours]
    
    !### Debug: for column mean OH lifetime
    REAL*8                       :: col_OHMASS
    REAL*8                       :: col_MASS
    REAL*8                       :: col_OHCONC

    ! Arrays                                            
    REAL(ESMF_KIND_R4),  POINTER :: lonCtr(:,:)   ! Lon centers [radians]
    REAL(ESMF_KIND_R4),  POINTER :: latCtr(:,:)   ! Lat centers [radians]

    !=======================================================================
    ! Initialization
    !=======================================================================

    __Iam__('Finalize_')

    ! Traceback info
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )
    Iam = TRIM( compName ) // '::' // TRIM( Iam )

    ! Get pointers to fields in import, internal, and export states
#   include "GEOSCHEMchem_GetPointer___.h"

    ! Get various parameters from the ESMF/MAPL framework
    CALL Extract_( GC, Clock,          &
                   GRID     = Grid,    &    ! ESMF Grid object
                   MAPLCF   = MaplCF,  &    ! ESMF Config obj (MAPL.rc)
                   GEOSCF   = GeosCF,  &    ! ESMF Config obj (GEOSCHEM*.rc)
                   IM       = IM,      &    ! # of longitudes on this PET
                   JM       = JM,      &    ! # of latitudes  on this PET
                   LM       = LM,      &    ! # of levels     on this PET
                   UTC      = utc,     &    ! Universal time [hours]
                   LONCTR   = lonCtr,  &    ! Longitude centers
                   LATCTR   = latCtr,  &    ! Latitude centers
                   LOCALPET = myPet,   &    ! PET # we are on now 
                   __RC__ )

    ! Fill the remaining fields of the IDENT object
    Ident%STDOUT_LUN = logLun
    Ident%PET        = myPet
    Ident%I_AM(1)    = TRIM( compName )
    Ident%I_AM(2)    = 'Finalize_'
    Ident%LEV        = 2
    Ident%ERRMSG     = ''
    Ident%VERBOSE    = .FALSE.

    !=======================================================================
    ! Print end-of-simulation output
    !=======================================================================

    ! Print the mean OH lifetime
    CALL Print_Mean_OH( GC, logLun, D_AIR_MASS, D_OH_MASS, __RC__ )

    !=======================================================================
    ! Finalize the Gridded Component
    !=======================================================================

    ! Call the FINALIZE method of the GEOS-Chem column chemistry code
    CALL GC_Column_Final( Ident, Coef, error )

    ! Trap the error from GEOS-Chem
    IF ( error /= SMV_SUCCESS ) THEN
       CALL Error_Trap_( Ident, error, __RC__ ) 
    ENDIF

    ! Finalize MAPL Generic
    CALL MAPL_GenericFinalize( GC, Import, Export, Clock, __RC__ )

    !=======================================================================
    ! All done
    !=======================================================================

    ! Reset fields in the IDENT object
    Ident%ERRMSG = ''
    Ident%LEV    = 2

    ! Add descriptive footer text
    WRITE( logLun, '(/,a)' ) REPEAT( '#', 79 )
    WRITE( logLun, 100     ) TRIM( Iam ), myPet
    WRITE( logLun, '(a)'   ) REPEAT( '#', 79 )

    ! Formats
100 FORMAT( '###',                                     /, &
            '### ', a, '  |  Cleanup on PET # ', i5.5, /  &
            '###' )

    ! Close file
    CLOSE( logLun )

    ! Successful return
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE Finalize_
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Extract_
!
! !DESCRIPTION: GC routine extracts several common quantities from the 
!  ESMF/MAPL environment so that they can be later passed down to the 
!  GEOS-Chem column chemistry code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Extract_( GC,       Clock,    Grid,   MaplCF,   GeosCF,   &
                       localPet, petCount, IM,     JM,       LM,       &
                       lonCtr,   latCtr,   latEdg, advCount, nymd,     &
                       nhms,     year,     month,  day,      dayOfYr,  &
                       hour,     minute,   second, utc,      hElapsed, &
                       tsChem,   tsDyn,    RC )
!
! !INPUT PARAMETERS:
!
    TYPE(ESMF_Clock),    INTENT(IN)            :: Clock       ! ESMF clock obj 
!                                                             
! !INPUT/OUTPUT PARAMETERS:                                   
!                                                             
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC          ! GC grid comp
!                                                             
! !OUTPUT PARAMETERS:                                         
!                                                             
    ! ESMF Grid object                                        
    TYPE(ESMF_Grid),     INTENT(OUT), OPTIONAL :: Grid        ! ESMF Grid obj
                                                              
    ! ESMF Config objects                                     
    TYPE(ESMF_Config),   INTENT(OUT), OPTIONAL :: MaplCF      ! MAPL.rc
    TYPE(ESMF_Config),   INTENT(OUT), OPTIONAL :: GeosCF      ! GEOSCHEM*.rc
                                                              
    ! ESMF virtual machine parameters                         
    INTEGER,             INTENT(OUT), OPTIONAL :: localPet    ! This PET
    INTEGER,             INTENT(OUT), OPTIONAL :: petCount    ! Total # of PETs
                                                              
    ! Horizontal grid coordinates                             
    INTEGER,             INTENT(OUT), OPTIONAL :: IM          ! # of lons
    INTEGER,             INTENT(OUT), OPTIONAL :: JM          ! # of lats
    INTEGER,             INTENT(OUT), OPTIONAL :: LM          ! # of levels
    REAL(ESMF_KIND_R4),  POINTER,     OPTIONAL :: lonCtr(:,:) ! Lon ctrs  [rad]
    REAL(ESMF_KIND_R4),  POINTER,     OPTIONAL :: latCtr(:,:) ! Lat ctrs  [rad]
    REAL(ESMF_KIND_R4),  POINTER,     OPTIONAL :: latEdg(:,:) ! Lat edges [rad]
                                                              
    ! Time variables     
   INTEGER(ESMF_KIND_I8),INTENT(OUT), OPTIONAL :: advCount    ! # of clock advs
    INTEGER,             INTENT(OUT), OPTIONAL :: nymd        ! YYYY/MM/DD
    INTEGER,             INTENT(OUT), OPTIONAL :: nhms        ! hh:mm:ss
    INTEGER,             INTENT(OUT), OPTIONAL :: year        
    INTEGER,             INTENT(OUT), OPTIONAL :: month       
    INTEGER,             INTENT(OUT), OPTIONAL :: day         
    INTEGER,             INTENT(OUT), OPTIONAL :: dayOfYr     
    INTEGER,             INTENT(OUT), OPTIONAL :: hour        
    INTEGER,             INTENT(OUT), OPTIONAL :: minute      
    INTEGER,             INTENT(OUT), OPTIONAL :: second      
    REAL,                INTENT(OUT), OPTIONAL :: utc         ! UTC time [hrs]
    REAL,                INTENT(OUT), OPTIONAL :: hElapsed    ! Elapsed hours
                                                              
    ! Timestep variables [min]                                
    REAL,                INTENT(OUT), OPTIONAL :: tsChem      ! Chemistry
    REAL,                INTENT(OUT), OPTIONAL :: tsDyn       ! Dynamics
    
    ! Return code 
    INTEGER,             INTENT(OUT), OPTIONAL :: RC          ! 0 = all is well
!
! !REMARKS:
!  If you need to obtain a quantity not returend by this routine, you can
!  manually extract it from the MaplCF or GeosCF configuration objects.
!
! !REVISION HISTORY:
!  01 Dec 2009 - A. Da Silva - Initial version
!  07 Apr 2010 - R. Yantosca - Added ProTeX headers
!  08 Apr 2010 - R. Yantosca - Make all outputs optional
!  08 Apr 2010 - R. Yantosca - Added outputs for localPet, petCount
!  08 Apr 2010 - R. Yantosca - Added outputs for individual time values
!                              as well as elapsed time (hours)
!  13 Apr 2010 - R. Yantosca - Now take tsDyn from the MAPL "RUN_DT:" setting
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
! 
    ! Objects
    TYPE(ESMF_Time)               :: startTime      ! ESMF start time obj
    TYPE(ESMF_Time)               :: currTime       ! ESMF current time obj
    TYPE(ESMF_TimeInterval)       :: elapsedTime    ! ESMF elapsed time obj
    TYPE(ESMF_VM)                 :: VM             ! ESMF VM object
    TYPE(GEOSCHEM_State), POINTER :: myState        ! Legacy state
    TYPE(GEOSCHEM_Wrap)           :: wrap           ! Wrapper for myState
    TYPE(MAPL_MetaComp),  POINTER :: metaComp       ! MAPL MetaComp object
    TYPE(ESMF_Array)              :: myArray

    ! Scalars
    CHARACTER(len=ESMF_MAXSTR)    :: compName       ! Gridded component name
    INTEGER(ESMF_KIND_I8)         :: count          ! # of clock advances
    INTEGER                       :: dims(3)        ! Array for dimensions
    INTEGER                       :: doy            ! Day of year (0-365/366)
    INTEGER                       :: yyyy, mm, dd   ! Year, month, day
    INTEGER                       :: h,    m,  s    ! Hour, minute, seconds
    REAL                          :: elapsedHours   ! Elapsed hours of run

    !=======================================================================
    ! Initialization
    !=======================================================================

    __Iam__('Extract_')

    ! Get my name and set-up traceback handle
    CALL ESMF_GridCompGet( GC, name=compName, vm=VM, __RC__ )
    Iam = TRIM( compName ) // '::' // TRIM( Iam )

    ! Get the internal state which holds the private Config object
    CALL ESMF_UserCompGetInternalState( GC, 'GEOSCHEM_State', wrap, STATUS )
    VERIFY_(STATUS)
    myState => wrap%ptr

    ! Assume successful return
    IF ( PRESENT( RC ) ) RC = ESMF_SUCCESS

    ! Zero variables
    dims = 0

    !=======================================================================
    ! Extract information from ESMF Config objects
    !=======================================================================

    ! Get the Config object based on "MAPL.rc"
    CALL ESMF_GridCompGet( GC, Config=MaplCF, __RC__ )
    
    ! Get the Config object based on "GEOSCHEMchem_GridComp.rc"
    GeosCF = myState%myCF

    ! Dynamic timestep
    IF ( PRESENT( tsDyn ) ) THEN
       CALL ESMF_ConfigGetAttribute( MaplCF, tsDyn,                       &
                                     Label="RUN_DT:",             __RC__ )
       tsDyn = tsDyn / 60
    ENDIF

    ! Chemistry timestep
    IF ( PRESENT( tsChem ) ) THEN
       CALL ESMF_ConfigGetAttribute( GeosCF, tsChem,                      &
                                     Label="CHEMISTRY_TIMESTEP:", __RC__ )
    ENDIF

    !=======================================================================
    ! Extract time/date information
    !=======================================================================
    
    ! Get the ESMF time object
    CALL ESMF_ClockGet( Clock,                    &
                        startTime    = startTime, &
                        currTime     = currTime,  &
                        advanceCount = count,     &
                         __RC__ )

    ! Get individual fields from the time object
    CALL ESMF_TimeGet( currTime, yy=yyyy, mm=mm, dd=dd, dayOfYear=doy, &
                                 h=h,     m=m,   s=s,   __RC__ )

    ! Save fields for return
    IF ( PRESENT( nymd     ) ) CALL MAPL_PackTime( nymd, yyyy, mm, dd )
    IF ( PRESENT( nhms     ) ) CALL MAPL_PackTime( nhms, h,    m,  s  )
    IF ( PRESENT( advCount ) ) advCount = count
    IF ( PRESENT( year     ) ) year     = yyyy
    IF ( PRESENT( month    ) ) month    = mm
    IF ( PRESENT( day      ) ) day      = dd
    IF ( PRESENT( dayOfYr  ) ) dayOfYr  = doy
    IF ( PRESENT( hour     ) ) hour     = h
    IF ( PRESENT( minute   ) ) minute   = m
    IF ( PRESENT( second   ) ) second   = s
    IF ( PRESENT( utc      ) ) utc      = ( DBLE( h )        ) + & 
                                          ( DBLE( m )/60d0   ) + &
                                          ( DBLE( s )/3600d0 )
 
    ! Compute elapsed time since start of simulation
    elapsedTime = currTime - startTime

    ! Get time fields from the elapsedTime object
    CALL ESMF_TimeIntervalGet( elapsedTime, h=h, m=m, s=s, __RC__ )

    ! Convert to decimal hours
    elapsedHours = DBLE( h ) + ( DBLE( m )/60d0 ) + ( DBLE( s )/3600d0 )
    
    ! Save fields for return
    IF ( PRESENT( hElapsed ) ) hElapsed = elapsedHours

    !=======================================================================
    ! Extract grid information
    !=======================================================================
    IF ( PRESENT( Grid ) ) THEN
    
       ! Get the ESMF grid attached to this gridded component
       CALL ESMF_GridCompGet( GC, grid=Grid, __RC__ )

       ! Get the dimensions on this PET
       CALL ESMF_GridGet( Grid,                                        &
                          localDE            = 0,                      &
                          staggerloc         = ESMF_STAGGERLOC_CENTER, &
                          computationalCount = dims,                   &
                          __RC__)
    ENDIF

    ! Save fields for return
    IF ( PRESENT( IM ) ) IM = dims(1)
    IF ( PRESENT( JM ) ) JM = dims(2)
    IF ( PRESENT( LM ) ) LM = dims(3)

    ! Get horizontal coordinate variables
    CALL MAPL_GetObjectFromGC( GC, metaComp, __RC__ )

    ! Longitude values on this PET
    IF ( PRESENT( lonCtr ) ) THEN
       CALL MAPL_Get( metaComp, lons=lonCtr, __RC__ )
    ENDIF

    ! Latitude values on this PET
    IF ( PRESENT( latCtr ) ) THEN
       CALL MAPL_Get( metaComp, lats=latCtr, __RC__ )
    ENDIF 

    ! Latitude edges on this PET
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%% NOTE: For now just hardwire latEdg = latCtr %%%
    !%%% This is only needed for the area_m2 array!  %%%
    !%%% Deal w/ this later. (bmy, 4/12/10)          %%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IF ( PRESENT( latEdg ) ) THEN
       CALL MAPL_Get( metaComp, lats=latEdg, __RC__ )
    ENDIF 

    !=======================================================================
    ! Extract information from ESMF VM object
    !=======================================================================

    ! Index of the PET we are on now
    IF ( PRESENT( localPet ) ) THEN
       CALL ESMF_VmGet( VM, localPet=localPet, __RC__ )
    ENDIF

    ! Total # of PETs used by this gridded component
    IF ( PRESENT( petCount ) ) THEN
       CALL ESMF_VmGet( VM, petCount=petCount, __RC__ )
    ENDIF

    !=======================================================================
    ! All done
    !=======================================================================
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE Extract_
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Setup_GeoLoc_
!
! !DESCRIPTION: This routine initializes the GEOLOC object, which passes 
!  the longitude, latitude, and local time to the GEOS-Chem column chemistry
!  code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Setup_GeoLoc_( lonRad, latRad, UTC, GeoLoc )
!
! !INPUT PARAMETERS:
!
    REAL,            INTENT(IN)    :: lonRad   ! Longitude [radians]
    REAL,            INTENT(IN)    :: latRad   ! Latitude  [radians]
    REAL,            INTENT(IN)    :: UTC      ! UTC       [decimal hours]
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(GC_GeoLoc), INTENT(INOUT) :: GeoLoc   ! Geographic location object
!
! !REVISION HISTORY:
!  13 Apr 2010 - R. Yantosca - Initial version
!  14 Jul 2010 - R. Yantosca - Now round off lonDeg, latDeg to 4 places
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    REAL*8                         :: lonDeg   ! Longitude  [degrees]
    REAL*8                         :: latDeg   ! Latitude   [degrees]
    REAL*8                         :: locTime  ! Local time [decimal hours]

    ! Convert radians to degrees 
    lonDeg           = lonRad * ( 180d0 / MAPL_PI )
    latDeg           = latRad * ( 180d0 / MAPL_PI )

    ! NOTE: Due to limited precision, sometimes lonDeg and latDeg are 
    ! sometimes returned as e.g. 79.9999 instead of 80.0000.  Apply an 
    ! algorithm to force exactitude to 4 decimal places.
    lonDeg           = RoundOff( lonDeg, 4 )
    latDeg           = RoundOff( latDeg, 4 )

    ! Compute Local time
    locTime          = UTC - ( lonDeg / 15d0 )
    if ( locTime < 0.0 ) locTime = locTime + 24.0d0
    
    ! Populate the GEOLOC object
    GeoLoc%Lon       = lonDeg
    GeoLoc%Lat       = latDeg
    GeoLoc%LocalTime = locTime

  END SUBROUTINE Setup_GeoLoc_
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Setup_Options_
!
! !DESCRIPTION: This routine sets up the OPTIONS object, which passes the
!  various option on/off switches down to the GEOS-Chem column chemistry code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Setup_Options_( compName, GeosCF, Options, RC )
!
! !INPUT PARAMETERS:
! 
    CHARACTER(LEN=*),  INTENT(IN)    :: compName  ! Name of Gridded Component
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_Config), INTENT(INOUT) :: GeosCF    ! Config obj for GEOSCHEM*.rc
    TYPE(GC_Options),  INTENT(INOUT) :: Options   ! Obj w/ G-C on/off switches
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  13 Apr 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!

    !=======================================================================
    ! Initialization
    !=======================================================================

    __Iam__('Setup_Options_')

    ! Traceback info
    Iam = TRIM( compName ) // '::' // TRIM( Iam )

    !=======================================================================
    ! Read from config file
    !=======================================================================

    !--------------------------
    ! Emissions Options
    !--------------------------
    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_EMISSIONS,             &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_EMISSIONS:",        __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_ANTHRO,                &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_ANTHRO:",           __RC__ )
                                                                             
    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_ANTHRO_BRAVO,          &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_ANTHRO_BRAVO:",     __RC__ )
                                                                            
    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_ANTHRO_CAC,            &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_ANTHRO_CAC:",       __RC__ )
                                                                             
    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_ANTHRO_EDGAR,          &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_ANTHRO_EDGAR:",     __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_ANTHRO_EPA,            &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_ANTHRO_EPA:",       __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_ANTHRO_VISTAS,         &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_ANTHRO_VISTAS:",    __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_ANTHRO_EPA,            &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_ANTHRO_EPA:",       __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_ANTHRO_EMEP,           &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_ANTHRO_EMEP:",      __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_BIOGENIC,              &
                                  Default = .FALSE.,                         &
                                  Label   ="USE_BIOGENIC:",          __RC__ )
                                                                             
    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_BIOMASS,               &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_BIOMASS:",          __RC__ )
                                                                             
    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_BIOMASS_GFED2,         &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_BIOMASS_GFED2:",    __RC__ )
 
    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_DEAD_DUST,             &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_DEAD_DUST:",        __RC__ )
                                                                             
    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_NOx_AIRCRAFT,          &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_NOx_AIRCRAFT:",     __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_NOX_LIGHTNING,         &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_NOx_LIGHTNING:",    __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_NOx_SOIL,              &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_NOx_SOIL:",         __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_SHIP_ARCTAS,           &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_SHIP_ARCTAS:",      __RC__ )

    !--------------------------
    ! Aerosol options
    !--------------------------
    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_CARBON_AEROSOLS,       &
                                  Default =.FALSE.,                          &
                                  Label   = "USE_CARBON_AEROSOLS:",  __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_DUST_AEROSOLS,         &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_DUST_AEROSOLS:",    __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_SEC_ORG_AEROSOLS,      &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_SEC_ORG_AEROSOLS:", __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_SEASALT_AEROSOLS,      &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_SEASALT_AEROSOLS:", __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_SULFATE_AEROSOLS,      &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_SULFATE_AEROSOLS:", __RC__ )

    !--------------------------
    ! Operations switches
    !--------------------------
    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_CHEMISTRY,             &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_CHEMISTRY:",        __RC__ )
        
    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_DRYDEP,                &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_DRYDEP:",           __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_CONVECTION,            &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_CONVECTION:",       __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_PBL_MIXING,            &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_PBL_MIXING:",       __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_WETDEP,                &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_WETDEP:",           __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, Options%USE_DEBUG_PRINT,           &
                                  Default = .FALSE.,                         &
                                  Label   = "USE_DEBUG_PRINT:",      __RC__ )

  END SUBROUTINE Setup_Options_
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Setup_DimInfo_
!
! !DESCRIPTION: This routine sets up the DIMINFO object, which passes down the
!  the following quantities to the GEOS-Chem column chemistry code.
!\\
!\\
!  The following members of DIMINFO are read from the configuration file
!  by this routine:
!
! \begin{itemize}
! \item N_AER: Number of aerosol species
! \item N_DUST: Number of Number of aerosol species
! \item N_RH: Number of relative humidity bins
! \item N_SOA_HC: Number of SOA HC classes (used to dimension GPROD/APROD)
! \item N_SOA_PROD: Number of SOA HC products (used to dimension GPROD/APROD)
! \item N_TRACERS: Number of advected tracers
! \end{itemize}
!
! wheareas the following members of the DIMINFO object are defined by the 
! Initialize method of the GEOS-Chem column code:
!
! \begin{itemize}
! \item N_DRYDEP: Number of dry deposition species
! \item N_JV: Number of J-value reactions
! \item N_SPECIES: Number of gas-phase species (active+inactive)
! \item N_REACTIONS: Number of reactions (gas-phase + photolysis)
! \end{itemize}
!
! !INTERFACE:
!
  SUBROUTINE Setup_DimInfo_( compName, GeosCF, DimInfo, RC )
!
! !INPUT PARAMETERS:
! 
    CHARACTER(LEN=*),  INTENT(IN)    :: compName  ! Name of Gridded Component
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_Config), INTENT(INOUT) :: GeosCF    ! Config obj for GEOSCHEM*.rc
    TYPE(GC_Dims),     INTENT(INOUT) :: DimInfo   ! Obj w/ dimension info
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC        ! Error return code
!
! !REVISION HISTORY:
!  13 Apr 2010 - R. Yantosca - Initial version
!  16 Apr 2010 - R. Yantosca - Now read N_SOA_HC and N_SOA_PROD
!EOP
!------------------------------------------------------------------------------
!BOC
    !=======================================================================
    ! Initialization
    !=======================================================================

    __Iam__('Setup_DimInfo_')

    ! Traceback info
    Iam = TRIM( compName ) // '::' // TRIM( Iam )

    !=======================================================================
    ! Read from config file
    !=======================================================================

    CALL ESMF_ConfigGetAttribute( GeosCF, DimInfo%N_AER,            &
                                  Default = 5,                      &
                                  Label   = "N_AER:",       __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, DimInfo%N_DUST,           &
                                  Default = 7,                      &
                                  Label   = "N_DUST:",      __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, DimInfo%N_MEMBERS,        &
                                  Default = 10,                     &
                                  Label   = "N_MEMBERS:",   __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, DimInfo%N_RH,             &
                                  Default = 5,                      &
                                  Label   = "N_RH:",        __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, DimInfo%N_SOA_HC,         &
                                  Default = 6,                      &
                                  Label   = "N_SOA_HC:",    __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, DimInfo%N_SOA_PROD,       &
                                  Default = 3,                     &
                                  Label   = "N_SOA_PROD:",  __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, DimInfo%N_TRACERS,        &
                                  Default = 54,                     &
                                  Label   = "N_TRACERS:",   __RC__ )

  END SUBROUTINE Setup_DimInfo_
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Setup_TrcName_
!
! !DESCRIPTION: This routine sets up the TRCNAME array, which holds the
!  names of all advected tracers in GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Setup_TrcName_( compName, GeosCF, trcName, RC )
!
! !USES:
!
    USE Charpak_Mod, ONLY : StrSplit
!
! !INPUT PARAMETERS:
! 
    CHARACTER(LEN=*),  INTENT(IN)    :: compName  ! Name of Gridded Component
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_Config), INTENT(INOUT) :: GeosCF    ! Config obj for GEOSCHEM*.rc
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*),  INTENT(OUT)   :: trcName(MAX_TRACERS)  ! Tracer names
    INTEGER,           INTENT(OUT)   :: RC                    ! Return code
!
! !REMARKS:
!  ESMF_ConfigGetAttribute cannot return a character array.  Therefore we
!  read in comma-separated strings with tags "TRACER_NAME1:", "TRACER_NAME2:",
!  "TRACER_NAME3:", "TRACER_NAME4:" and then concatenate them into a long
!  string.  Routine StrSplit from charpak_mod.F is then used to split the
!  strings into a substring array.
!
! !REVISION HISTORY:
!  13 Apr 2010 - R. Yantosca - Initial version
!  14 Apr 2010 - R. Yantosca - Now read comma-separated strings and then
!                              concatenate them into a big string variable
!                              for splitting. 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES: 
!
    INTEGER                          :: N, nSubStr
    CHARACTER(LEN=ESMF_MAXSTR)       :: str1, str2, str3, str4
    CHARACTER(LEN=500)               :: str,  subStr(255)

    !=======================================================================
    ! Initialization
    !=======================================================================

    __Iam__('Setup_TrcName_')

    ! Traceback info
    Iam = TRIM( compName ) // '::' // TRIM( Iam )

    !=======================================================================
    ! Read from config file
    !=======================================================================

    ! Get the list of names of advected tracers
    CALL ESMF_ConfigGetAttribute( GeosCF, str1,                      &
                                  Default = "",                      &
                                  Label   = "TRACER_NAME1:", __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, str2,                      &
                                  Default = "",                      &
                                  Label   = "TRACER_NAME2:", __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, str3,                      &
                                  Default = "",                      &
                                  Label   = "TRACER_NAME3:", __RC__ )

    CALL ESMF_ConfigGetAttribute( GeosCF, str4,                      &
                                  Default = "",                      &
                                  Label   = "TRACER_NAME4:", __RC__ )

    ! Concatenate the strings together
    str = TRIM( str1 ) // TRIM( str2 ) // TRIM( str3 ) // TRIM( str4 )

    ! Split by spaces
    CALL StrSplit( str, ',', subStr, nSubStr )
    
    ! Copy into TRCNAME array
    DO N = 1, nSubStr
       trcName(N) = TRIM( subStr(N) )
    ENDDO

  END SUBROUTINE Setup_TrcName_
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Error_Trap_
!
! !DESCRIPTION: This routine stops the run and prints the name of the
!  offending routine if an error condition is returned.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Error_Trap_( Ident, error, RC )
!
! !INPUT PARAMETERS: 
!
    TYPE(GC_IDENT),   INTENT(IN)  :: Ident      ! Obj w/ info from ESMF
    INTEGER,          INTENT(IN)  :: error      ! Error code from GEOS-Chem
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC         ! Return error code
!
! !REVISION HISTORY: 
!  22 Jun 2009 - R. Yantosca - Initial version
!  15 Jul 2009 - R. Yantosca - Updated for drydep, wetdep, PBL mixing
!  24 Aug 2009 - R. Yantosca - Updated for emissions reader etc. routines
!  03 Nov 2009 - R. Yantosca - Now trap error in the GC_INTERFACE
!  03 Nov 2009 - R. Yantosca - Cosmetic changes
!  14 Dec 2009 - R. Yantosca - Now trap errors in the GC_COLUMN routines
!  14 Apr 2010 - R. Yantosca - Adapted for use in GEOSCHEMchem_GridCompMod.F90
!  29 Apr 2010 - R. Yantosca - Now print error traceback info from IDENT
!  30 Apr 2010 - R. Yantosca - Now use 5 digits for PET
!  06 May 2010 - R. Yantosca - Remove redundant error codes
!  03 Jun 2010 - R. Yantosca - Remove more redundant error codes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N

    !=======================================================================
    ! Initialization
    !=======================================================================

     __Iam__('Error_Trap_')

    ! Traceback info
    Iam = TRIM( Ident%I_AM(1) ) // '::' // TRIM( Iam )

    !=======================================================================
    ! Error trap
    !=======================================================================

    IF ( error == SMV_FAILURE ) THEN
       
       !--------------------------------------------------
       ! Print error traceback information
       !--------------------------------------------------
   
       ! Begin header
       WRITE( logLun, '(/,a)' ) REPEAT( '=', 79 )
       WRITE( logLun,  20     ) Ident%PET
       WRITE( logLun,  21     ) TRIM( Ident%ERRMSG )
       WRITE( logLun, '(/,a) ') 'Error traceback:' 

       ! Write the calling sequence of routines
       DO N = Ident%LEV, 1, -1
          WRITE( loglun, 22   ) N, TRIM( Ident%I_AM(N) )
       ENDDO

       ! End header
       WRITE( logLun, '(a,/)' ) REPEAT( '=', 79 )
       CLOSE( logLun          )

       ! Return w/ failure
       RETURN_(ESMF_FAILURE)

    ELSE IF ( error == SMV_FAIL_GINOUX      ) THEN
       WRITE( logLun, 10 ) 'FAILURE in SRC_FUNC_GINOUX',         Ident%PET
       RETURN_(ESMF_FAILURE)
    ELSE IF ( error == SMV_FAIL_EMDUSTBOX   ) THEN
       WRITE( logLun, 10 ) 'FAILURE in EMISSDUST_BOX',           Ident%PET
       RETURN_(ESMF_FAILURE)
    ELSE

       ! This is now deprecated
       PRINT*, 'RC', error
       WRITE( logLun, 10 ) 'UNKNOWN FAILURE',                    Ident%PET
       RETURN_(ESMF_FAILURE)

    ENDIF

    ! FORMAT string
 10 FORMAT( a, ' PET=', i5.5                             )
 20 FORMAT( 'GEOS-Chem error encountered on PET: ', i5.5 )
 21 FORMAT( 'Message: ', a                               )
 22 FORMAT( i4, ' : ', a                                 )


  END SUBROUTINE Error_Trap_
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
    REAL*8,  INTENT(IN) :: X   ! Number to be rounded
    INTEGER, INTENT(IN) :: N   ! Number of decimal places to keep
!
! !RETURN VALUE:
!
    REAL*8              :: Y   ! Number rounded to N decimal places
!
! !REMARKS:
!  The algorithm to round X to N decimal places is as follows:
!  (1) Multiply X by 10**(N+1)
!  (2) If X < 0, then add -5 to X; otherwise add 5 to X
!  (3) Round X to nearest integer
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
!  14 Jul 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Round and truncate X to N decimal places
    Y = INT( NINT( X*(10d0**(N+1)) + SIGN( 5d0, X ) ) / 10d0 ) / (10d0**N)

  END FUNCTION RoundOff
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: print_mean_oh
!
! !DESCRIPTION: Subroutine Print\_Mean\_OH prints the average mass-weighted OH 
!  concentration at the end of a simulation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_Mean_OH( GC, logLun, airMass, ohMass, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC             ! Ref to this GridComp

    INTEGER,             INTENT(IN)    :: logLun         ! LUN for stdout print
    REAL,     POINTER,   INTENT(IN)    :: airMass(:,:,:) ! Air mass [molec air]
    REAL,     POINTER,   INTENT(IN)    :: ohMass (:,:,:) ! Mass-weighted OH
                                                         !  [molec OH/cm3 * 
                                                         !   molec air]
!
! !INPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC             ! Return code
!
! !REMARKS:
!  The AIRMASS and OHMASS variables are used to pass the data values from the 
!  ESMF state, which is REAL*4.  We need to make sure that we do not store into
!  the ESMF state any values that exceed 1e38, which is the maximum allowable 
!  value.  Therefore, AIRMASS and OHMASS will store the values divided by the
!  scale factor OH_SCALE = 1d20 in order to prevent this overflow situation.
!
! !REVISION HISTORY: 
!  01 Jul 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    CHARACTER(LEN=ESMF_MAXSTR) :: compName      ! Name of gridded component
    REAL*8                     :: SUM_OHMASS   
    REAL*8                     :: SUM_MASS   
    REAL*8                     :: OHCONC
 
    !=======================================================================
    ! Initialization
    !=======================================================================

     __Iam__('Print_Mean_OH')

    ! Traceback info
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )
    Iam = TRIM( compName ) // '::' // TRIM( Iam )
 
    !=======================================================================
    ! Print mean OH values
    !=======================================================================
    
    ! Total Mass-weighted OH [molec OH/cm3] * [molec air]
    SUM_OHMASS = GlobalSum( GC, ohMass,  __RC__ )

    ! Atmospheric air mass [molec air]
    SUM_MASS   = GlobalSum( GC, airMass, __RC__ )

    ! Restore proper values by applying the OH scale factor
    ! (This is necessary in order avoid overflow)
    SUM_OHMASS = SUM_OHMASS * OH_SCALE
    SUM_MASS   = SUM_MASS   * OH_SCALE

    ! Avoid divide-by-zero errors 
    IF ( SUM_MASS > 0d0 ) THEN 
            
       ! Divide OH by [molec air] and report as [1e5 molec/cm3]
       OHCONC = ( SUM_OHMASS / SUM_MASS ) / 1d5
         
       ! Write value to log file
       WRITE( logLun, '(/,a)' ) REPEAT( '=', 79 ) 
       WRITE( logLun, *       ) 'Mass-Weighted OH Concentration'
       WRITE( logLun, *       ) 'Mean OH = ', OHCONC, ' [1e5 molec/cm3]' 
       WRITE( logLun, '(  a)' ) REPEAT( '=', 79 ) 

    ELSE

       ! Write error msg if SUM_MASS is zero
       WRITE( logLun, '(/,a)' ) REPEAT( '=', 79 ) 
       WRITE( logLun, '(  a)' ) 'Could not print mass-weighted OH!'
       WRITE( logLun, '(  a)' ) 'Atmospheric air mass is zero!'
       WRITE( logLun, '(  a)' ) REPEAT( '=', 79 ) 
       
    ENDIF

  END SUBROUTINE Print_Mean_OH
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GlobalSum
!
! !DESCRIPTION: Subroutine GlobalSum prints the sum (or minimum, or maximum)
!  of an array across all PETs.  Calls the ESMF_VMAllFullReduce function
!  to do the array reduction.  The default is to compute the array sum
!  unless MINIMUM=.TRUE. or MAXIMUM=.TRUE.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GlobalSum( GC, dataPtr, minimum, maximum, RC ) RESULT( value )
!
! !INPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC              ! Gridcomp name
    REAL,    POINTER,    INTENT(IN)    :: dataPtr(:,:,:)  ! Data array
    LOGICAL, OPTIONAL,   INTENT(IN)    :: minimum         ! Compute minimum?
    LOGICAL, OPTIONAL,   INTENT(IN)    :: maximum         ! Compute maximum?
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC              ! Return code
!
! !RETURN VALUE:
!
    REAL                               :: value           ! Sum, max, or min
! 
! !REVISION HISTORY: 
!  01 Jul 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Objects
    TYPE(ESMF_VM)              :: VM            ! ESMF VM object

    ! Scalars
    INTEGER                    :: nSize         ! Size of 3-D array
    CHARACTER(LEN=ESMF_MAXSTR) :: compName      ! Gridded component name

    ! Arrays 
    REAL, ALLOCATABLE          :: data1d(:)     ! 1-D array for reduction

    !=======================================================================
    ! Initialization
    !=======================================================================

     __Iam__('GlobalSum')

    ! Traceback info
    CALL ESMF_GridCompGet( GC, name=compName, vm=VM, __RC__ )
    Iam = TRIM( compName ) // '::' // TRIM( Iam )

    !=======================================================================
    ! Do the reduction operation across all CPU's: sum, max, or min
    !=======================================================================

    ! Create a 1-D vector
    nSize = SIZE( dataPtr )
    ALLOCATE( data1d( nSize ), STAT=STATUS )
    VERIFY_(STATUS)

    ! Rearrange into a 1-D array
    data1d = RESHAPE( dataPtr, (/1/) )
    
    ! Compute the sum over all PETS
    IF ( PRESENT( maximum ) ) THEN
       CALL ESMF_VMAllFullReduce( VM, data1d, value, nSize, ESMF_MAX, __RC__ )
    ELSE IF ( PRESENT( minimum ) ) THEN
       CALL ESMF_VMAllFullReduce( VM, data1d, value, nSize, ESMF_MIN, __RC__ )
    ELSE 
       CALL ESMF_VMAllFullReduce( VM, data1d, value, nSize, ESMF_SUM, __RC__ )
    ENDIF

    ! Deallocate temporary array
    IF( ALLOCATED( data1D ) ) DEALLOCATE( data1d )

  END FUNCTION GlobalSum

 END MODULE GEOSCHEMchem_GridCompMod
