#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_environment_mod
!
! !DESCRIPTION: Module GIGC\_ENVIRONMENT\_MOD establishes the runtime 
!  environment for the Grid-Independent GEOS-Chem (aka "GIGC") model.  It is 
!  designed to receive model parameter and geophysical environment information 
!  and allocate memory based upon it.
!\\
!\\
!  It provides routines to do the following:
!
! \begin{itemize}
! \item Allocate geo-spatial arrays
! \item Initialize met. field derived type.
! \item Initialize CHEM, PHYS, and EMISSIONS states
! \end{itemize}

!  NOTE: This is mostly for testing the grid-independent code in the current 
!  GEOS-Chem.  Many of these inputs will come from the GEOS-5 interface. 
!  It will remain in DEVEL state for some time.
!\\
!\\
! !INTERFACE: 
!
MODULE GIGC_Environment_Mod
!
! !USES
!        
  IMPLICIT NONE
# include "define.h"

  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: GIGC_Allocate_All
  PUBLIC  :: GIGC_Init_All
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Get_nSchm_nSchmBry
!
! !REMARKS:
!  For consistency, we should probably move the met state initialization
!  to the same module where the met state derived type is contained.
!
! !REVISION HISTORY:
!  26 Jan 2012 - M. Long     - Created module file
!  13 Aug 2012 - R. Yantosca - Added ProTeX headers
!  19 Oct 2012 - R. Yantosca - Removed routine INIT_LOCAL_MET, this is now
!                              handled in Headers/gigc_state_met_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed to gigc_environment_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC        
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_allocate_all
!
! !DESCRIPTION: Subroutine GIGC\_ALLOCATE\_ALL allocates all LAT/LON 
!  ALLOCATABLE arrays for global use by the GEOS-Chem either as a standalone 
!  program or module.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Allocate_All( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE CMN_DEP_MOD,       ONLY : SET_CMN_DEP_MOD
    USE CMN_NOX_MOD,       ONLY : SET_CMN_NOX_MOD
    USE CMN_O3_MOD,        ONLY : SET_CMN_O3_MOD
    USE CMN_MOD,           ONLY : SET_CMN_MOD
    USE CMN_FJ_MOD,        ONLY : SET_CMN_FJ_MOD
    USE CMN_SIZE_MOD,      ONLY : SET_CMN_SIZE_MOD
    USE CMN_DIAG_MOD,      ONLY : SET_CMN_DIAG_MOD
    USE COMODE_LOOP_MOD,   ONLY : SET_COMODE_LOOP_MOD
    USE COMMSOIL_MOD,      ONLY : SET_COMMSOIL_MOD
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod
    USE JV_CMN_MOD,        ONLY : SET_JV_CMN_MOD
    USE VDIFF_PRE_MOD,     ONLY : SET_VDIFF_PRE_MOD

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
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
!  Need to add better error checking and exit upon failure.
!
! !REVISION HISTORY: 
!  26 Jan 2012 - M. Long     - Initial version
!  13 Aug 2012 - R. Yantosca - Added ProTeX headers
!  17 Oct 2012 - R. Yantosca - Add am_I_Root, RC as arguments
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Allocate_All
!  30 Oct 2012 - R. Yantosca - Now pass am_I_Root, RC to SET_COMMSOIL_MOD
!  01 Nov 2012 - R. Yantosca - Now zero the fields of the Input Options object
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Initialize fields of the Input Options object
    CALL Set_GIGC_Input_Opt( am_I_Root, Input_Opt, RC )

    ! Allocate module fields with the locally-determined 
    ! longitude (IIPAR) and latitude (JJPAR) dimensions
    CALL SET_CMN_SIZE_MOD
    CALL SET_CMN_DEP_MOD
    CALL SET_CMN_DIAG_MOD
    CALL SET_CMN_NOX_MOD
    CALL SET_CMN_O3_MOD
    CALL SET_CMN_MOD
    CALL SET_CMN_FJ_MOD
    CALL SET_COMMSOIL_MOD   ( am_I_Root, RC )
    CALL SET_COMODE_LOOP_MOD( am_I_Root, RC )
    CALL SET_JV_CMN_MOD
    
    CALL SET_VDIFF_PRE_MOD
          
  END SUBROUTINE GIGC_Allocate_All
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_init_all
!
! !DESCRIPTION: Subroutine GIGC\_INIT\_ALL initializes the top-level data 
!  structures that are either passed to/from GC or between GC components 
!  (emis->transport->chem->etc)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Init_All( State_Met, State_Chm, am_I_Root, RC ) 
!
! !USES:
!
    USE CMN_Size_Mod,       ONLY : IIPAR, JJPAR, LLPAR, NBIOMAX
    USE Comode_Loop_Mod,    ONLY : IGAS
    USE GIGC_ErrCode_Mod
    USE GIGC_State_Chm_Mod
    USE GIGC_State_Met_Mod
    USE Logical_Mod,        ONLY : LSCHEM
    USE Tracer_Mod,         ONLY : N_TRACERS
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  Need to add better error checking, currently we just return upon error.
!
! !REVISION HISTORY: 
!  26 Jan 2012 - M. Long     - Initial version
!  13 Aug 2012 - R. Yantosca - Added ProTeX headers
!  16 Oct 2012 - R. Yantosca - Renamed LOCAL_MET argument to State_Met
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE  argument to State_Chm
!  16 Oct 2012 - R. Yantosca - Call Init_Chemistry_State (in gc_type2_mod.F90,
!                              which was renamed from INIT_CHEMSTATE)
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_errcode_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference IGAS in Headers/comode_loop_mod.F
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Init_All
!  26 Oct 2012 - R. Yantosca - Now call Get_nSchm, nSchmBry to find out the
!                              number of strat chem species and Bry species
!  01 Nov 2012 - R. Yantosca - Now use LSCHEM from logical_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    INTEGER :: nSchm, nSchmBry

    !=======================================================================
    ! Initialize object for met fields
    !=======================================================================
    CALL Init_GIGC_State_Met( am_I_Root  = am_I_Root,   &
                              IM         = IIPAR,       &
                              JM         = JJPAR,       &
                              LM         = LLPAR,       &
                              State_Met  = State_Met,   &
                              RC         = RC          )

    ! Return upon error
    IF ( RC /= GIGC_SUCCESS ) RETURN

    !=======================================================================
    ! Initialize object for chemical state
    !=======================================================================

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ) 

    !-----------------------------------------------------------------------
    !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
    !
    ! At present, we have not included stratospheric chemistry in when
    ! connecting to the ESMF interface.  We will have to include all strat
    ! chem species via the ESMF import state.  For now just set nsChm
    ! and nSchmBry to default values. (bmy, 11/1/12)
    !-----------------------------------------------------------------------
    nSchm    = 1
    nSchmBry = 1

#else

    !-----------------------------------------------------------------------
    !   %%%%% TESTING GIGC INTERFACE FROM EXISTING GEOS-CHEM %%%%%
    !
    ! We can test the grid-independent implementation of stratospheric
    ! chemistry when compiling the traditional GEOS-Chem with DEVEL=yes.
    ! In this case, we need to pre-compute the # of strat chem species
    ! (nSchm) and the number of bromine species (nSchmBry) so that we can
    ! allocate the corresponding fields of the chemistry state.
    ! (bmy, 11/1/12)
    !-----------------------------------------------------------------------
    IF ( LSCHEM ) THEN

       ! Strat chem is turned on, find out the # of stratospheric 
       ! chemistry species for which we need to read rates from disk.
       ! NOTE: Bromine species are handled specially.
       CALL Get_nSchm_nSchmBry( am_I_Root  = am_I_Root,  &  ! Root CPU (Y/N)?
                                nSchm      = nSchm,      &  ! # strat chem spec
                                nSchmBry   = nSchmBry,   &  ! # strat chem spec
                                RC         = RC         )   ! Success or failure

    ELSE

       ! Strat chem is turned off
       nSchm    = 0
       nSchmBry = 0

    ENDIF

#endif

    ! Initialize chemistry state
    CALL Init_GIGC_State_Chm(  am_I_Root  = am_I_Root,   &  ! Root CPU (Y/N)?
                               IM         = IIPAR,       &  ! # of lons
                               JM         = JJPAR,       &  ! # of lats
                               LM         = LLPAR,       &  ! # of levels
                               nTracers   = N_TRACERS,   &  ! # of tracers
                               nBioMax    = NBIOMAX,     &  ! # biomass species
                               nSpecies   = IGAS,        &  ! # chemical species 
                               nSchm      = nSchm,       &  ! # strat chem spec
                               nSchmBry   = nSchmBry,    &  ! # bromine species
                               State_Chm  = State_Chm,   &  ! Chemistry State
                               RC         = RC          )   ! Success or failure
    
    ! Return upon error
    IF ( RC /= GIGC_SUCCESS ) RETURN

  END SUBROUTINE GIGC_Init_All
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_nSchm_nSchmBry
!
! !DESCRIPTION: Subroutine Get\_nSchm\_nSchmBry finds out the # of 
!  stratospheric chemistry tracers and bromine tracers so that we can
!  allocate the various Schm_* fields in the Chemistry State object.
!\\
!\\
! !INTERFACE:
!      
  SUBROUTINE Get_nSchm_nSchmBry( am_I_Root, nSchm, nSchmBry, RC )
!
! !USES:
!
    USE CMN_SIZE_MOD
    USE GIGC_ErrCode_Mod
    USE LOGICAL_MOD,       ONLY : LLINOZ
    USE TRACER_MOD,        ONLY : N_TRACERS, TRACER_NAME, STT
    USE TIME_MOD,          ONLY : GET_TAU, GET_NYMD, GET_NHMS, GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: nSchm       ! # of strat chem species
    INTEGER, INTENT(OUT) :: nSchmBry    ! # of strat chem Bry species
    INTEGER, INTENT(OUT) :: RC          ! Success or failure
! 
! !REVISION HISTORY:
!  01 Feb 2011 - L. Murray   - Initial version
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  26 Oct 2012 - R. Yantosca - Now pass Chemistry State object for GIGC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: GMI_name, GC_Name
    CHARACTER(LEN=255) :: GMI_TrName(120)
    INTEGER            :: AS, N, NN

    !=================================================================
    ! INIT_STRAT_CHEM begins here!
    !=================================================================

    ! Assume success
    RC         = GIGC_SUCCESS

    ! Initialize counters, initial times, mapping arrays
    nsChm      = 0
    nSchmBry   = 0

    ! List of available tracers with archived monthly climatological
    ! production rates, loss frequencies, and mixing ratios from the 
    ! GMI Combo model (tracer names here are as used in GMI).
    GMI_TrName = (/ &
             'A3O2',     'ACET',   'ACTA',   'ALD2',    'ALK4',  'ATO2', &
             'B3O2',       'Br',   'BrCl',    'BrO',  'BrONO2',  'C2H6', &
             'C3H8',     'CCl4', 'CF2Br2', 'CF2Cl2', 'CF2ClBr', 'CF3Br', &
           'CFC113',   'CFC114', 'CFC115',  'CFCl3',    'CH2O', 'CH3Br', &
          'CH3CCl3',    'CH3Cl',    'CH4',     'CO',      'Cl',   'Cl2', &
            'Cl2O2',      'ClO', 'ClONO2',    'EOH',    'ETO2',   'ETP', &
             'GCO3',     'GLYC',   'GLYX',     'GP',    'GPAN',     'H', &
               'H2',    'H2402',    'H2O',   'H2O2',     'HAC',   'HBr', &
         'HCFC141b', 'HCFC142b', 'HCFC22',  'HCOOH',     'HCl',  'HNO2', &
             'HNO3',     'HNO4',    'HO2',   'HOBr',    'HOCl',  'IALD', &
             'IAO2',      'IAP',   'INO2',   'INPN',    'ISN1',  'ISNP', &
             'ISOP',      'KO2',   'MACR',   'MAN2',    'MAO3',  'MAOP', &
              'MAP',     'MCO3',    'MEK',   'MGLY',     'MO2',   'MOH', &
               'MP',     'MRO2',    'MRP',    'MVK',    'MVN2',     'N', &
              'N2O',     'N2O5',     'NO',    'NO2',     'NO3',   'NOx', &
                'O',      'O1D',     'O3',   'OClO',      'OH',    'Ox', &
              'PAN',      'PMN',    'PO2',     'PP',     'PPN',  'PRN1', &
             'PRPE',     'PRPN',   'R4N1',   'R4N2',    'R4O2',   'R4P', &
             'RA3P',     'RB3P',   'RCHO',   'RCO3',   'RCOOH',  'RIO1', &
             'RIO2',      'RIP',    'ROH',     'RP',    'VRO2',   'VRP'    /)

    !=====================================================================
    ! Determine the number of stratospheric species & bromine species
    ! defined as GEOS-Chem advected tracers
    !=====================================================================

    ! Loop over all GEOS-Chem advected tracers
    DO N = 1, N_TRACERS

       ! GEOS-Chem advected tracer name
       GC_Name = TRACER_NAME(N) 
    
       !---------------------------------------------------------------
       ! For now, guarantee that GMI prod/loss rates are not used for  
       ! any bromine species
       !---------------------------------------------------------------
       IF ( ( TRIM( GC_Name ) ==      'Br' )   .or. &
            ( TRIM( GC_Name ) ==    'BrCl' )   .or. &
            ( TRIM( GC_Name ) ==     'BrO' )   .or. &
            ( TRIM( GC_Name ) ==  'BrONO2' )   .or. &
            ( TRIM( GC_Name ) ==  'CF2Br2' )   .or. &
            ( TRIM( GC_Name ) == 'CF2ClBr' )   .or. &
            ( TRIM( GC_Name ) ==   'CF3Br' )   .or. &
            ( TRIM( GC_Name ) ==   'CH3Br' )   .or. &
            ( TRIM( GC_Name ) ==     'HBr' )   .or. &
            ( TRIM( GC_Name ) ==    'HOBr' ) ) THEN

          ! Increment # of Bromine tracers
          nSchmBry = nSchmBry + 1

          ! Skip to next GEOS-Chem tracer
          CYCLE
       ENDIF

       ! Loop over all possible stratospheric species
       DO NN = 1, SIZE( GMI_TrName )     

          ! Stratospheric species name according to GMI
          GMI_Name = TRIM( GMI_TrName(NN) )

          ! Some species names don't exactly match GEOS-Chem names
          !IF ( TRIM(GMI_TrName(NN)) .eq. 'BrONO2' ) GMI_Name = 'BrNO3'
          
          !---------------------------------------------------------------
          ! Increment nSchm for each match 
          !---------------------------------------------------------------
          IF ( TRIM( GC_Name ) == TRIM( GMI_Name ) ) THEN
                
             IF ( LLINOZ .and. ( TRIM( GC_Name ) == 'Ox' ) ) THEN
                IF ( am_I_Root ) THEN
                   WRITE( 6, '(a)' ) TRIM( GC_Name ) // ' (via Linoz)'
                ENDIF
             ELSE IF ( TRIM( GC_Name ) == 'Ox' ) THEN
                IF ( am_I_Root ) THEN
                   WRITE( 6, '(a)' ) TRIM( GC_Name ) // ' (via Synoz)'
                ENDIF
             ELSE
                IF ( am_I_Root ) THEN
                   WRITE( 6, '(a)' ) TRIM( GC_Name ) //' (via GMI rates)'
                ENDIF
             ENDIF
             
             nSchm = nSchm + 1
             EXIT
          ENDIF
          
       ENDDO
    ENDDO

  END SUBROUTINE Get_nSchm_nSchmBry
!EOC
END MODULE GIGC_Environment_Mod
#endif
