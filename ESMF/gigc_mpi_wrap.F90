#if defined( ESMF_ )
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_mpi_wrap
!
! !DESCRIPTION: Module GIGC\_MPI\_WRAP is the module containing MPI-based
!  routines used for the ESMF interface to the Grid-Independent
!  GEOS-Chem (aka "GIGC").
!\\
!\\
! !INTERFACE:
!
MODULE GIGC_Mpi_Wrap
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GIGC_Input_Bcast
  PUBLIC :: GIGC_Idt_Bcast
  PUBLIC :: GIGC_Reader_Bcast
  PUBLIC :: GIGC_Readchem_Bcast
  PUBLIC :: GIGC_Bcast_Char 
  PUBLIC :: GIGC_Bcast_Int
  PUBLIC :: GIGC_Bcast_Real8
  PUBLIC :: mpiComm
!
! !REMARKS:
!  These routines are needed to broadcast values read in from ASCII input
!  (which are read only on the root CPU) to all other CPUs.
!                                                                             .
!  NOTE: If you add values to a derived type object (e.g. Input_Opt), then
!  you will also have to add the proper calls so that the extra fields will
!  get broadcasted to other CPUs.
!                                                                             .
!  NOTE: The SMVGEAR init functions READER and READCHEM touch several
!  variables and arrays in a very convoluted manner.  It may be difficult
!  to try to call these on the root CPU and then broadcast to all other
!  CPUs.  For now just call READER and READCHEM on all CPUs (bmy, 3/7/13)
!
! !REVISION HISTORY:
!  03 Jan 2013 - M. Long     - Initial version
!  07 Mar 2013 - R. Yantosca - Added more ProTeX headers + comments
!EOP
!------------------------------------------------------------------------------
!BOC
  INTEGER, SAVE :: mpiComm

CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_input_bcast
!
! !DESCRIPTION: Routine GIGC\_INPUT\_BCAST broadcasts contents of the Input_Opt
!  read in from the input.geos\_\_\_.rc input file from the Root process to all
!  processes in the MPI_COMM_WORLD.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Gigc_Input_Bcast( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE CMN_SIZE_Mod
    USE Drydep_Mod
    USE GIGC_Input_Opt_Mod,  ONLY : OptInput
    USE GIGC_Errcode_Mod,    ONLY : GIGC_SUCCESS
    USE M_MPIF
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
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  In lieu of using MPI_Type_Struct in order to broadcast a derived type
!  object, which would be simpler, the choice was made to individually
!  parse out individual elements of the Input_Opt object. In this way
!  it is more straight forward and easily interpreted by future GIGC
!  users without much MPI experience. The intent of this routine should
!  be self-evident without reference to MPI dicumentation.
!  MSL - 04 Jan 2013
!
! !REVISION HISTORY:
!  04 Jan 2013 - M. Long     - Initial version
!  28 Feb 2013 - R. Yantosca - Now MPI BCast the Input_Opt%haveImpRst field
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
! mpi_bcast (buffer, count, datatype, root, comm, ier)

    RC = GIGC_SUCCESS

    !----------------------------------------
    ! SIZE PARAMETER fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%MAX_DIAG, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%MAX_TRCS, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%MAX_MEMB, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%MAX_FAMS, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%MAX_DEP,  1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! SIMULATION MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%NYMDb, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%NHMSb, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%NYMDe, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%NHMSe, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%RUN_DIR, len(INPUT_OPT%RUN_DIR), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%IN_RST_FILE, len(INPUT_OPT%IN_RST_FILE), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LSVGLB, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%OUT_RST_FILE, len(INPUT_OPT%OUT_RST_FILE), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%DATA_DIR, len(INPUT_OPT%DATA_DIR), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%GCAP_DIR, len(INPUT_OPT%GCAP_DIR), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%GEOS_4_DIR, len(INPUT_OPT%GEOS_4_DIR), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%GEOS_5_DIR, len(INPUT_OPT%GEOS_5_DIR), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%GEOS_57_DIR, len(INPUT_OPT%GEOS_57_DIR), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%MERRA_DIR, len(INPUT_OPT%MERRA_DIR), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%DATA_DIR_1x1, len(INPUT_OPT%DATA_DIR_1x1), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TEMP_DIR, len(INPUT_OPT%TEMP_DIR), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LUNZIP, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LWAIT, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LVARTROP, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%NESTED_I0, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%NESTED_J0, 1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! TRACER MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%N_TRACERS, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ID_TRACER(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TRACER_NAME(:), (255)*INPUT_OPT%MAX_TRCS, mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TRACER_MW_G, INPUT_OPT%MAX_TRCS, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TRACER_MW_KG, INPUT_OPT%MAX_TRCS, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TCVV, INPUT_OPT%MAX_TRCS, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%XNUMOL, INPUT_OPT%MAX_TRCS, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TRACER_N_CONST(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TRACER_CONST, len(INPUT_OPT%TRACER_CONST)*INPUT_OPT%MAX_TRCS*INPUT_OPT%MAX_MEMB, mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TRACER_COEFF, INPUT_OPT%MAX_TRCS*INPUT_OPT%MAX_MEMB, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ID_EMITTED(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LSPLIT, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ITS_A_RnPbBe_SIM, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ITS_A_CH3I_SIM, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ITS_A_FULLCHEM_SIM, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ITS_A_HCN_SIM, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ITS_A_TAGOX_SIM, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ITS_A_TAGCO_SIM, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ITS_A_C2H6_SIM, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ITS_A_CH4_SIM, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ITS_AN_AEROSOL_SIM, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ITS_A_MERCURY_SIM, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ITS_A_CO2_SIM, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ITS_A_H2HD_SIM, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ITS_A_POPS_SIM, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ITS_NOT_COPARAM_OR_CH4, 1, mpi_logical, 0, mpiComm, RC )

    !----------------------------------------
    ! AEROSOL MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%LSULF, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LCRYST, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LCARB, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LSOA, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LDUST, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LDEAD, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LSSALT, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LDICARB, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%SALA_REDGE_um(:), 2, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%SALC_REDGE_um(:), 2, mpi_real8, 0, mpiComm, RC )

    !----------------------------------------
    ! EMISSIONS MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%LEMIS, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TS_EMIS, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LANTHRO, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%FSCALYR, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LEMEP, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LBRAVO, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LEDGAR, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LSTREETS, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LCAC, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LNEI05, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LRETRO, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LNEI99, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LICARTT, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LVISTAS, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LBIOFUEL, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LBIOGENIC, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LMEGAN, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LPECCA, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LMEGANMONO, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ISOP_SCALING, 1, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LBIOMASS, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LBBSEA, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LTOMSAI, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LGFED2BB, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%L8DAYBB, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%L3HRBB, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LSYNOPBB, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LGFED3BB, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LDAYBB3, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%L3HRBB3, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LAIRNOX, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LLIGHTNOX, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LOTDLOC, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LSOILNOX, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LFERTILIZERNOX, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%NOx_SCALING, 1, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LEDGARSHIP, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LICOADSSHIP, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LEMEPSHIP, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LSHIPSO2, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LARCSHIP, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LCOOKE, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LAVHRRLAI, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LMODISLAI, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LHIST, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%HISTYR, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LWARWICK_VSLS, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LSSABr2, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LFIX_PBL_BRO, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%Br_SCALING, 1, mpi_real8, 0, mpiComm, RC )

    !----------------------------------------
    ! CO2 MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%LGENFF, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LANNFF, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LMONFF, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LCHEMCO2, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LSEASBB, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LBIODAILY, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LBIODIURNAL, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LBIONETORIG, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LBIONETCLIM, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LOCN1997, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LOCN2009ANN, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LOCN2009MON, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LSHIPEDG, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LSHIPICO, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LPLANE, 1, mpi_logical, 0, mpiComm, RC )

    !----------------------------------------
    ! FUTURE MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%LFUTURE, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%FUTURE_YEAR, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%FUTURE_SCEN, len(INPUT_OPT%FUTURE_SCEN), mpi_character, 0, mpiComm, RC )

    !----------------------------------------
    ! CHEMISTRY MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%LCHEM, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LSCHEM, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LLINOZ, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TS_CHEM, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LSVCSPEC, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LKPP, 1, mpi_logical, 0, mpiComm, RC )

    !----------------------------------------
    ! CHEMISTRY MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%LTRAN, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LMFCT, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LFILL, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TPCORE_IORD, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TPCORE_JORD, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TPCORE_KORD, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TS_DYN, 1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! CONVECTION MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%LCONV, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LTURB, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LNLPBL, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TS_CONV, 1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! DEPOSITION MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%LDRYD, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LWETD, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%USE_OLSON_2001, 1, mpi_logical, 0, mpiComm, RC )

    !----------------------------------------
    ! GAMAP MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%GAMAP_DIAGINFO, len(INPUT_OPT%GAMAP_DIAGINFO), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%GAMAP_TRACERINFO, len(INPUT_OPT%GAMAP_TRACERINFO), mpi_character, 0, mpiComm, RC )

    !----------------------------------------
    ! OUTPUT MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%NJDAY(:), 366, mpi_integer, 0, mpiComm, RC )
    
    !----------------------------------------
    ! DIAGNOSTIC MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%ND01, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD01, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND02, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD02, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND03, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD03, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND04, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD04, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND05, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD05, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND06, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD06, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND07, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD07, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND08, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD08, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND09, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD09, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND10, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD10, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND11, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD11, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND12, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD12, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND13, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD13, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND14, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD14, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND15, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD15, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND16, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD16, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND17, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD17, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND18, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD18, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND19, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD19, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND20, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD20, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND21, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD21, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND22, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD22, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND23, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD23, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND24, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD24, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND25, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD25, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND26, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD26, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND27, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD27, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND28, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD28, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND29, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD29, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND30, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD30, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND31, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD31, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND32, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD32, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND33, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD33, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND34, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD34, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND35, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD35, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND36, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD36, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND37, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD37, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND38, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD38, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND39, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD39, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND40, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD40, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND41, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD41, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND42, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD42, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND43, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD43, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND44, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD44, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND45, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD45, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND46, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD46, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND47, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD47, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND48, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD48, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND49, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD49, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND50, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD50, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD51, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND52, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD52, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND53, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD53, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND54, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD54, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND55, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD55, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND56, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD56, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND57, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD57, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND58, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD58, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND59, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD59, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND60, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD60, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND61, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD61, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND62, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD62, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND63, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD63, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND64, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD64, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND66, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD66, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND67, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD67, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND68, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD68, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND69, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD69, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND70, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD70, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LPRT, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TINDEX(:,:), INPUT_OPT%MAX_DIAG*INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TCOUNT(:)  , INPUT_OPT%MAX_DIAG, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TMAX(:)    , INPUT_OPT%MAX_DIAG, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! PLANEFLIGHT MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%DO_PF, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%PF_IFILE, len(INPUT_OPT%PF_IFILE), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%PF_OFILE, len(INPUT_OPT%PF_OFILE), mpi_character, 0, mpiComm, RC )

    !----------------------------------------
    ! ND48 MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%DO_ND48, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND48_FILE, len(INPUT_OPT%ND48_FILE), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND48_FREQ, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND48_N_STA, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND48_IARR(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND48_JARR(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND48_LARR(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND48_NARR(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! ND49 MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%DO_ND49, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND49_FILE, len(INPUT_OPT%ND49_FILE), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND49_TRACERS(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND49_FREQ, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND49_IMIN, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND49_IMAX, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND49_JMIN, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND49_JMAX, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND49_LMIN, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND49_LMAX, 1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! ND50 MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%DO_ND50, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND50_FILE, len(INPUT_OPT%ND50_FILE), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LND50_HDF, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND50_TRACERS(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND50_IMIN, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND50_IMAX, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND50_JMIN, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND50_JMAX, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND50_LMIN, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND50_LMAX, 1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! ND51 MENU fields
     !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%DO_ND51, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51_FILE, len(INPUT_OPT%ND51_FILE), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LND51_HDF, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51_TRACERS(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51_HR_WRITE, 1, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51_HR1, 1, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51_HR2, 1, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51_IMIN, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51_IMAX, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51_JMIN, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51_JMAX, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51_LMIN, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51_LMAX, 1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! ND51b MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%DO_ND51b, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51b_FILE, len(INPUT_OPT%ND51b_FILE), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LND51b_HDF, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51b_TRACERS(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51b_HR_WRITE, 1, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51b_HR1, 1, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51b_HR2, 1, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51b_IMIN, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51b_IMAX, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51b_JMIN, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51b_JMAX, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51b_LMIN, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND51b_LMAX, 1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! ND63 MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%DO_ND63, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND63_FILE, len(INPUT_OPT%ND63_FILE), mpi_character, 0, mpiComm, RC )
!    CALL MPI_Bcast( INPUT_OPT%ND63_TRACERS(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND63_FREQ, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND63_IMIN, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND63_IMAX, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND63_JMIN, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND63_JMAX, 1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! PROD LOSS MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%DO_SAVE_PL, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LFAMILY, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ND65, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LD65, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%DO_SAVE_O3, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%NFAM, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%FAM_COEF(:,:), INPUT_OPT%MAX_MEMB*INPUT_OPT%MAX_FAMS, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%FAM_MEMB, len(INPUT_OPT%FAM_MEMB)*INPUT_OPT%MAX_MEMB*INPUT_OPT%MAX_FAMS, mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%FAM_NAME, len(INPUT_OPT%FAM_NAME)*INPUT_OPT%MAX_FAMS, mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%FAM_NMEM, INPUT_OPT%MAX_FAMS, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%FAM_TYPE, len(INPUT_OPT%FAM_TYPE)*INPUT_OPT%MAX_FAMS, mpi_character, 0, mpiComm, RC )

    !----------------------------------------
    ! UNIX CMDS fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%BACKGROUND, len(INPUT_OPT%BACKGROUND), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%REDIRECT, len(INPUT_OPT%REDIRECT), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%REMOVE_CMD, len(INPUT_OPT%REMOVE_CMD), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%SEPARATOR, len(INPUT_OPT%SEPARATOR), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%WILD_CARD, len(INPUT_OPT%WILD_CARD), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%UNZIP_CMD, len(INPUT_OPT%UNZIP_CMD), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%ZIP_SUFFIX, len(INPUT_OPT%ZIP_SUFFIX), mpi_character, 0, mpiComm, RC )

    !----------------------------------------
    ! NESTED GRID MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%LWINDO, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LWINDO2x25, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LWINDO_NA, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TPBC_DIR_NA, len(INPUT_OPT%TPBC_DIR_NA), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LWINDO_EU, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TPBC_DIR_EU, len(INPUT_OPT%TPBC_DIR_EU), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LWINDO_CH, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TPBC_DIR_CH, len(INPUT_OPT%TPBC_DIR_CH), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LWINDO_SE, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TPBC_DIR_SE, len(INPUT_OPT%TPBC_DIR_SE), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LWINDO_CU, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%TPBC_DIR, len(INPUT_OPT%TPBC_DIR), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%NESTED_TS, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%NESTED_I1, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%NESTED_J1, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%NESTED_I2, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%NESTED_J2, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%NESTED_I0W, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%NESTED_J0W, 1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! BENCHMARK MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%LSTDRUN, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%STDRUN_INIT_FILE, len(INPUT_OPT%STDRUN_INIT_FILE), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%STDRUN_FINAL_FILE, len(INPUT_OPT%STDRUN_FINAL_FILE), mpi_character, 0, mpiComm, RC )

    !----------------------------------------
    ! ARCHIVED OH MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%OH_DIR, len(INPUT_OPT%OH_DIR), mpi_character, 0, mpiComm, RC )

    !----------------------------------------
    ! O3PL MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%O3PL_DIR, len(INPUT_OPT%O3PL_DIR), mpi_character, 0, mpiComm, RC )

    !----------------------------------------
    ! MERCURY MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%ANTHRO_Hg_YEAR, 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%HG_SCENARIO, len(INPUT_OPT%Hg_SCENARIO), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%USE_CHECKS, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LDYNOCEAN, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LPREINDHG, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%Hg_RST_FILE, len(INPUT_OPT%Hg_RST_FILE), mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LGTMM, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%GTMM_RST_FILE, len(INPUT_OPT%GTMM_RST_FILE), mpi_character, 0, mpiComm, RC )

    !----------------------------------------
    ! CH4 MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%LCH4BUD, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LGAO, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LCOL, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LLIV, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LWAST, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LBFCH4, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LRICE, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LOTANT, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LBMCH4, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LWETL, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LSOABS, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%LOTNAT, 1, mpi_logical, 0, mpiComm, RC )

    !----------------------------------------
    ! APM MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%IFNUCL, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%FE0, 1, mpi_real8, 0, mpiComm, RC )

    !----------------------------------------
    ! POPS MENU fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%POP_TYPE, 3, mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%CHEM_PROCESS, 1, mpi_logical, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%POP_EMISFILE, 255, mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%POP_XMW, 1, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%POP_KOA, 1, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%POP_KBC, 1, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%POP_K_POPG_OH, 1, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%POP_K_POPP_O3A, 1, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%POP_K_POPP_O3B, 1, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%POP_HSTAR, 1, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%POP_DEL_H, 1, mpi_real8, 0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%POP_DEL_Hw, 1, mpi_real8, 0, mpiComm, RC )

    !----------------------------------------
    ! DRYDEP and DUST fields from input.geos
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%NUMDEP,   1,         mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%NDVZIND,  Input_Opt%MAX_DEP,    mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%IDDEP,    NDSTBIN,   mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%DUSTREFF, NDSTBIN,   mpi_real8,     0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%DUSTDEN,  NDSTBIN,   mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( INPUT_OPT%DEPNAME,  14*Input_Opt%MAX_DEP, mpi_character, 0, mpiComm, RC )

    !----------------------------------------
    ! GEOS-5 GCM INTERFACE fields
    !----------------------------------------
    CALL MPI_Bcast( INPUT_OPT%haveImpRst, 1, mpi_logical, 0, mpiComm, RC )

  END SUBROUTINE GIGC_Input_Bcast
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_idt_bcast
!
! !DESCRIPTION: Routine GIGC\_IDT_BCAST broadcasts the tracer flags (IDTxxxx)
!  and species flags (IDxxxx), etc., that are computed by the routines in
!  GeosCore/tracerid_mod.F.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Idt_Bcast( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE Drydep_Mod
    USE GIGC_Errcode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE Tracerid_Mod
    USE M_MPIF
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  04 Jan 2013 - M. Long     - Initial version
!  07 Mar 2013 - R. Yantosca - Added more ProTeX headers + comments
!   7 Mar 2013 - R. Yantosca - Reordered for clarity + cosmetic changes

!EOP
!------------------------------------------------------------------------------
!BOC

    !=======================================================================
    ! MPI Broadcast the IDxxxx chemical species flags
    !=======================================================================
    CALL MPI_Bcast( IDO3        , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDNO2       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDNO3       , 1, mpi_integer, 0, mpiComm, RC ) 
    CALL MPI_Bcast( IDN2O5      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDHNO4      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDOX        , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDNOX       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDHC1       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDNO        , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDHNO2      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDCO        , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDPRPE      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDISOP      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDALK4      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDC3H8      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDPAN       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDGLPAN     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDGPAN      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDPMN       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDPPN       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDHNO3      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDOH        , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDHO2       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDH2O2      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDACET      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDMEK       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDALD2      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDRCHO      , 1, mpi_integer, 0, mpiComm, RC ) 
    CALL MPI_Bcast( IDMVK       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDMACR      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDISN2      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDR4N2      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDCH2O      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDC2H6      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDMP        , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDDMS       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDSO2       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDSO4       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDMSA       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDDRYO3     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDDRYPAN    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDDRYNO2    , 1, mpi_integer, 0, mpiComm, RC )       
!- Not used (ccc, 01/27/10)     
!    CALL MPI_Bcast( IDBENZ      , 1, mpi_integer, 0, mpiComm, RC )
!    CALL MPI_Bcast( IDTOLU      , 1, mpi_integer, 0, mpiComm, RC )
!    CALL MPI_Bcast( IDXYLE      , 1, mpi_integer, 0, mpiComm, RC )
!    CALL MPI_Bcast( IDMONX      , 1, mpi_integer, 0, mpiComm, RC )
!    CALL MPI_Bcast( IDGLYX      , 1, mpi_integer, 0, mpiComm, RC )
!    CALL MPI_Bcast( IDMGLY      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDDRYGLYX   , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDDRYMGLY   , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDC2H4      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDC2H2      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDMBO       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDGLYC      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDHAC       , 1, mpi_integer, 0, mpiComm, RC )
                               
    CALL MPI_Bcast( IDAPAN      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDENPAN     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDMPAN      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDNIPAN     , 1, mpi_integer, 0, mpiComm, RC )

    CALL MPI_Bcast( IDDRYAPAN   , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDDRYENPAN  , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDDRYGLPAN  , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDDRYGPAN   , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDDRYMPAN   , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDDRYNIPAN  , 1, mpi_integer, 0, mpiComm, RC )

    ! Aromatics
    CALL MPI_Bcast( IDLBRO2H    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDLBRO2N    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDLTRO2H    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDLTRO2N    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDLXRO2H    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDLXRO2N    , 1, mpi_integer, 0, mpiComm, RC )

    ! Bromine
    CALL MPI_Bcast( IDBr2       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBr        , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBrO       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDHBr       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDHOBr      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBrNO2     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBrNO3     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDCHBr3     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDCH2Br2    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDCH3Br     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDDRYHOBr   , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDDRYHBr    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDDRYBrNO3  , 1, mpi_integer, 0, mpiComm, RC )

    !=======================================================================
    ! MPI Broadcast the IDTxxxx tracer ID flags
    !=======================================================================

    ! Standard simulation
    CALL MPI_Bcast( IDTNOX      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTOX       , 1, mpi_integer, 0, mpiComm, RC )  
    CALL MPI_Bcast( IDTPAN      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTCO       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTALK4     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTISOP     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTHNO3     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTH2O2     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTACET     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTMEK      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTALD2     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTRCHO     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTMVK      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTMACR     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTPMN      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTPPN      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTISN2     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTR4N2     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTPRPE     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTC3H8     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTCH2O     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTC2H6     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTN2O5     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTHNO4     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTMP       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTDMS      , 1, mpi_integer, 0, mpiComm, RC ) 
    CALL MPI_Bcast( IDTSO2      , 1, mpi_integer, 0, mpiComm, RC ) 
    CALL MPI_Bcast( IDTMSA      , 1, mpi_integer, 0, mpiComm, RC ) 
    CALL MPI_Bcast( IDTNH3      , 1, mpi_integer, 0, mpiComm, RC ) 
    CALL MPI_Bcast( IDTSO4      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTNH4      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTNIT      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTBCPI     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTBCPO     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTOCPO     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTOCPI     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTDST1     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTDST2     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTDST3     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTDST4     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTSALA     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTSALC     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTBr2      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTBr       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTBrO      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTHBr      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTHOBr     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTBrNO2    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTBrNO3    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTCHBr3    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTCH2Br2   , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTCH3Br    , 1, mpi_integer, 0, mpiComm, RC )

    ! SOA mechanism
    CALL MPI_Bcast( IDTALPH     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTLIMO     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTALCO     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTSOA1     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTSOA2     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTSOA3     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTSOA4     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTSOA5     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTSOAG     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTSOAM     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTSOG1     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTSOG2     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTSOG3     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTSOG4     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTSOG5     , 1, mpi_integer, 0, mpiComm, RC )

    ! Crystalline tracers (not used)
    CALL MPI_Bcast( IDTAS       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTAHS      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTNH4aq    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTLET      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTSO4aq    , 1, mpi_integer, 0, mpiComm, RC )

    ! SOA tracers
    CALL MPI_Bcast( IDTSOA5     , 1, mpi_integer, 0, mpiComm, RC )

    ! Dicarbonyls simulation
    CALL MPI_Bcast( IDTGLYX     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTMGLY     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTBENZ     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTTOLU     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTXYLE     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTMONX     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTC2H4     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTC2H2     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTMBO      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTGLYC     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTHAC      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTAPAN     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTENPAN    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTGLPAN    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTGPAN     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTMPAN     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTNIPAN    , 1, mpi_integer, 0, mpiComm, RC )

    ! Fabien Paulot ISOP mechanism
    CALL MPI_Bcast( IDTISOPN    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTAP       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTPROPNN   , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTMOBA     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTMMN      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTRIP      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTIEPOX    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTMAP      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTPYPAN    , 1, mpi_integer, 0, mpiComm, RC )

    ! H2/HD simulation
    CALL MPI_Bcast( IDTH2       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTHD       , 1, mpi_integer, 0, mpiComm, RC )

    ! Rn/Pb/Be simulation
    CALL MPI_Bcast( IDTRN       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTPB       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDTBE7      , 1, mpi_integer, 0, mpiComm, RC )

    ! Tagged Ox simulation
    CALL MPI_Bcast( IDTOxStrt   , 1, mpi_integer, 0, mpiComm, RC )

    !=======================================================================
    ! MPI Broadcast the GEOS-CHEM Emission ID flags
    !=======================================================================
    CALL MPI_Bcast( NEMANTHRO   , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( NEMBIOG     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDENOX      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEOX       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDECO       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEPRPE     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEC3H8     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEALK4     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEC2H6     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEACET     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEMEK      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEALD2     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEISOP     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDECH2O     , 1, mpi_integer, 0, mpiComm, RC ) 
    CALL MPI_Bcast( IDEHNO3     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEGLYX     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEMGLY     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEBENZ     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDETOLU     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEXYLE     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEMONX     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEC2H4     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEC2H2     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEMBO      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEGLYC     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEHAC      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDECHBr3    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDECH2Br2   , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDEBr2      , 1, mpi_integer, 0, mpiComm, RC )

    !=======================================================================
    ! MPI Broadcast the IDBFxxxx flags for BIOFUELS
    !=======================================================================
    CALL MPI_Bcast( IDBFNOX     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBFCO      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBFALK4    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBFACET    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBFMEK     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBFALD2    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBFPRPE    , 1, mpi_integer, 0, mpiComm, RC ) 
    CALL MPI_Bcast( IDBFC3H8    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBFCH2O    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBFC2H6    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBFBENZ    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBFTOLU    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBFXYLE    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBFC2H2    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBFC2H4    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBFGLYC    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBFGLYX    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBFMGLY    , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBFHAC     , 1, mpi_integer, 0, mpiComm, RC )

    !=======================================================================
    ! MPI Broadcast the IDBxxxx flags for BIOMASS BURNING
    !=======================================================================
    CALL MPI_Bcast( IDBNOX      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBCO       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBALK4     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBACET     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBMEK      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBALD2     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBPRPE     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBC3H8     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBCH2O     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBC2H6     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBSO2      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBNH3      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBBC       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBOC       , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBCO2      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBTOLU     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBBENZ     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBXYLE     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBC2H2     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBC2H4     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBHAC      , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBGLYC     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBGLYX     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBMGLY     , 1, mpi_integer, 0, mpiComm, RC )
    CALL MPI_Bcast( IDBCH4      , 1, mpi_integer, 0, mpiComm, RC )

    !=======================================================================
    ! MPI Broadcast the tagged Hg index arrays
    !=======================================================================
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

         ! Initialize category flags
!- eds 8/31/10 -------------------------------------------------------
       CALL MPI_Bcast( ID_Hg_tot, 1, mpi_integer, 0, mpiComm, RC )
!       CALL MPI_Bcast( ID_Hg_na, 1, mpi_integer, 0, mpiComm, RC )
!       CALL MPI_Bcast( ID_Hg_eu, 1, mpi_integer, 0, mpiComm, RC )
!       CALL MPI_Bcast( ID_Hg_as, 1, mpi_integer, 0, mpiComm, RC )
!       CALL MPI_Bcast( ID_Hg_rw, 1, mpi_integer, 0, mpiComm, RC )
!       CALL MPI_Bcast( ID_Hg_oc, 1, mpi_integer, 0, mpiComm, RC )
!       CALL MPI_Bcast( ID_Hg_ln, 1, mpi_integer, 0, mpiComm, RC )
!       CALL MPI_Bcast( ID_Hg_nt, 1, mpi_integer, 0, mpiComm, RC )

       CALL MPI_Bcast( ID_Hg_can, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_usa, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_cam, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_sam, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_waf, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_eaf, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_saf, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_naf, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_eur, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_eeu, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_mde, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_sov, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_sas, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_eas, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_sea, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_jpn, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_oce, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_so,  1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_bb,  1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_geo, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_atl, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_nat, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_sat, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_npa, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_arc, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_ant, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_ocn, 1, mpi_integer, 0, mpiComm, RC )
       CALL MPI_Bcast( ID_Hg_str, 1, mpi_integer, 0, mpiComm, RC )
    ENDIF

    ! Return success
    RC = GIGC_SUCCESS

  END SUBROUTINE GIGC_Idt_Bcast
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_reader_bcast
!
! !DESCRIPTION: Routine GIGC\_READER\_BCAST performs an MPI broadcast for 
!  all of the namelist data that is read from the "mglob.dat" file
!  by routine GeosCore/reader.F.  This is one of the SMVGEAR input files.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Reader_Bcast( RC )
!
! !USES:
!
    USE Comode_Loop_Mod
    USE GIGC_Errcode_Mod
    USE M_MPIF
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC   ! Success or failure
!
! !REMARKS:
!  NOTE: The READER (GeosCore/reader.F) subroutine is a tangled mess.  It not 
!  only reads values from mglob.dat but it also sets up other arrays used 
!  elsewhere in the SMVGEAR code.  It may just be simpler to run READER on
!  all CPUs so as not to have to worry about doing the MPI broadcast properly.
!  (bmy, 3/7/13)
!
! !REVISION HISTORY:
!  04 Jan 2013 - M. Long     - Initial version
!  07 Mar 2013 - R. Yantosca - Added ProTex header
!EOP
!------------------------------------------------------------------------------
!BOC

    !=======================================================================
    ! MPI Broadcast fields of /CTLFLG/ namelist
    !=======================================================================
    CALL MPI_Bcast( IFSOLVE     , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( ITESTGEAR   , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( IFURBAN     , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( IFTROP      , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( IFSTRAT     , 1,          mpi_integer,   0, mpiComm, RC )

    !=======================================================================   
    ! MPI Broadcast fields of /CTLDIM/ namelist
    !=======================================================================
    CALL MPI_Bcast( KULOOP      , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( LYOUT       , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( LXOUT       , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( LZOUT       , 1,          mpi_integer,   0, mpiComm, RC )

    !=======================================================================   
    ! MPI Broadcast fields of /CTLTIM/ namelist
    !=======================================================================
    CALL MPI_Bcast( CHEMINTV    , 1,          mpi_real8,     0, mpiComm, RC )

    !=======================================================================   
    ! MPI Broadcast fields of /CTLPRT/ namelist
    !=======================================================================
    CALL MPI_Bcast( IPRATES     , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( IPREADER    , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( IOSPEC      , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( IOREAC      , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( APGASA      , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( APGASB      , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( APGASC      , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( APGASD      , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( APGASE      , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( APGASF      , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( APGASG      , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( APGASH      , 1,          mpi_integer,   0, mpiComm, RC )

    !=======================================================================   
    ! MPI Broadcast fields of /CLGEAR/ namelist
    !=======================================================================
    CALL MPI_Bcast( IFREORD     , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( FRACDEC     , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( PLOURB      , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( PLOTROP     , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( HMAXNIT     , 1,          mpi_integer,   0, mpiComm, RC )

    !=======================================================================   
    ! MPI Broadcast other relevant fields
    !=======================================================================
    CALL MPI_Bcast( NMASBAL     , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( NAMEMB      , IMASBAL*14, mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( MLOP        , ILAT*ILONG, mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( NREAD       , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( NTLOOP      , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( NVERT       , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( NLAYER      , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( KGLC        , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( IGLOBCHEM   , 1,          mpi_integer,   0, mpicomm, RC )
    CALL MPI_Bcast( NCSALL      , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( NCSTRST     , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( NCSURBAN    , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( NCSTROP     , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( NCSSTRAT    , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( NCSGAS      , 1,          mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( NTLOOPNCS   , ICS,        mpi_integer,   0, mpiComm, RC )
    !CALL MPI_Bcast( NOUSE       , ICS,        mpi_integer,   0, mpiComm, RC )
    !CALL MPI_Bcast( NRATCUR     , ICS,        mpi_integer,   0, mpiComm, RC )

    ! Return success
    RC = GIGC_SUCCESS

  END SUBROUTINE GIGC_Reader_Bcast
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_readchem_bcast
!
! !DESCRIPTION: Routine GIGC\_READCHEM\_BCAST performs an MPI broadcast for 
!  data that is read from the SMVGEAR "globchem.dat" input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_ReadChem_Bcast( C1,    CSTRAT, CTROPL, CTROPS, CURBAN, &
                                  ININT, IORD,   NCOF,   RC              )
!
! !USES:
!
    USE Comode_Loop_Mod
    USE GIGC_Errcode_Mod
    USE M_MPIF
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL*8,  INTENT(INOUT) :: C1, CSTRAT, CTROPL, CTROPS, CURBAN
    INTEGER, INTENT(INOUT) :: ININT(10), IORD, NCOF
!
! !OUPTUT PARAMETERS:
!
    INTEGER, INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  NOTE: The READCHEM (GeosCore/readchem.F) subroutine is a tangled mess.  
!  It not only reads values from mglob.dat but it also sets up other arrays 
!  used elsewhere in the SMVGEAR code.  It may just be simpler to run READER 
!  on all CPUs so as not to have to worry about doing the MPI broadcast 
!  properly. (bmy, 3/7/13)
!
! !REVISION HISTORY:
!  04 Jan 2013 - M. Long     - Initial version
!  07 Mar 2013 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Globchem.dat
    CALL MPI_Bcast( C1,     1,        mpi_real8,     0, mpiComm, RC ) ! IN
    CALL MPI_Bcast( CSTRAT, 1,        mpi_real8,     0, mpiComm, RC ) ! IN
    CALL MPI_Bcast( CTROPL, 1,        mpi_real8,     0, mpiComm, RC ) ! IN
    CALL MPI_Bcast( CTROPS, 1,        mpi_real8,     0, mpiComm, RC ) ! IN
    CALL MPI_Bcast( CURBAN, 1,        mpi_real8,     0, mpiComm, RC ) ! IN
    CALL MPI_Bcast( ININT,  10,       mpi_integer,   0, mpiComm, RC ) ! IN
    CALL MPI_Bcast( DINP,   4,        mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( IORD,   1,        mpi_integer,   0, mpiComm, RC ) ! IN
    CALL MPI_Bcast( NCOF,   1,        mpi_integer,   0, mpiComm, RC ) ! IN
    CALL MPI_Bcast( SPECL,  MXCOF*2,  mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( JST,    14,       mpi_character, 0, mpiComm, RC ) 
    CALL MPI_Bcast( ARRT,   MXCOF,    mpi_real8,     0, mpiComm, RC )
    CALL MPI_Bcast( BRRT,   MXCOF,    mpi_real8,     0, mpiComm, RC )
    CALL MPI_Bcast( KCRRT,  MXCOF,    mpi_integer,   0, mpiComm, RC )
    CALL MPI_Bcast( RINP,   IMISC*14, mpi_character, 0, mpiComm, RC )
    CALL MPI_Bcast( PINP,   IMISC,    mpi_real8,     0, mpiComm, RC )
    CALL MPI_Bcast( XINP,   IMISC*14, mpi_character, 0, mpiComm, RC )
    
    !     Chemga.dat
    !      READ(I,*)
    !      READ(I,610) ASTKCF
    !      READ(I,*)
    !      READ(I,610)
    !      READ(I,*)
    
    RC = GIGC_SUCCESS
    
  END SUBROUTINE GIGC_Readchem_Bcast
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_bcast_char
!
! !DESCRIPTION: Wrapper routine to do an MPI broadcast operation on
!  a CHARACTER string variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Bcast_Char( VAL, SIZE, RC )
!
! !USES:
!
    USE GIGC_Errcode_Mod
    USE M_MPIF
!
! !INPUT PARAMETERS:
!
    INTEGER,      INTENT(IN)    :: SIZE     ! # of characters
!
! !INPUT/OUTPUT PARAMETERS:
!
    CHARACTER(*), INTENT(INOUT) :: VAL(:)   ! Character string

!
! !OUTPUT PARAMETERS:
!
    INTEGER,      INTENT(OUT)   :: RC       ! Success or failure?
!
! !REMARKS:
!  Mostly experimental.
!
! !REVISION HISTORY:
!  04 Jan 2013 - M. Long     - Initial version
!  07 Mar 2013 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
    
    ! Assume success
    RC = GIGC_SUCCESS

    ! Do MPI broadcast
    CALL MPI_Bcast( VAL, SIZE, mpi_character, 0, mpiComm, RC )
    
  END SUBROUTINE GIGC_Bcast_Char
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_bcast_int
!
! !DESCRIPTION: Wrapper routine to do an MPI broadcast operation on
!  an INTEGER variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Bcast_Int( VAL, SIZE, RC )
!
! !USES:
!
    USE GIGC_Errcode_Mod
    USE M_MPIF
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)    :: SIZE        ! Size of variable to be broadcasted
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT) :: VAL(SIZE)   ! Variable to be broadcasted
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  Mostly experimental
!
! !REVISION HISTORY:
!  04 Jan 2013 - M. Long     - Initial version
!  07 Mar 2013 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC

    RC = GIGC_SUCCESS
    
    CALL MPI_Bcast( VAL, SIZE, mpi_integer, 0, mpiComm, RC )
    
  END SUBROUTINE GIGC_Bcast_Int
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_bcast_real8
!
! !DESCRIPTION: Wrapper routine to do an MPI broadcast operation on
!  an REAL*8 variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Bcast_Real8( VAL, SIZE, RC )
!
! !USES:
!
    USE GIGC_Errcode_Mod
    USE M_MPIF
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)    :: SIZE        ! Size of variable to be broadcast
! 
! !INPUT/OUTPUT PARAMETERS:
!
    REAL*8,  INTENT(INOUT) :: VAL(SIZE)   ! Variable to be broadcast
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  Mostly experimental
!
! !REVISION HISTORY:
!  04 Jan 2013 - M. Long     - Initial version
!  07 Mar 2013 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
    
    RC = GIGC_SUCCESS
    
    CALL MPI_Bcast( VAL, SIZE, mpi_real8, 0, mpiComm, RC )
    
  END SUBROUTINE GIGC_Bcast_Real8
!EOC
END MODULE GIGC_Mpi_Wrap
#endif
