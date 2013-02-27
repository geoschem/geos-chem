#if defined( ESMF_ )
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc\_mpi\_wrap
!
! !DESCRIPTION: Module GIGC\_MPI\_WRAP is the module containing MPI-based
!  routines used for the ESMF interface to the Grid-Independent
!  GEOS-Chem (aka "GIGC").
!\\
!\\
! !INTERFACE:
!
MODULE GIGC_MPI_WRAP
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GIGC_INPUT_BCAST
  PUBLIC :: GIGC_IDT_BCAST
  PUBLIC :: GIGC_READER_BCAST
  PUBLIC :: GIGC_READCHEM_BCAST
  PUBLIC :: GIGC_BCAST_CHAR, GIGC_BCAST_INT, GIGC_BCAST_REAL8
  PUBLIC :: mpiComm
!
! !REVISION HISTORY:
!  03 Jan 2013 - M. Long     - Initial version
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
  SUBROUTINE GIGC_INPUT_BCAST( Input_Opt, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod,  ONLY : OptInput
    USE GIGC_Errcode_Mod,    ONLY : GIGC_SUCCESS
    USE DRYDEP_MOD
    USE CMN_SIZE_MOD
    USE M_MPIF
!
! !INPUT PARAMETERS:
!
!!    INTEGER,        INTENT(IN)    :: mpicomm     ! Mpi communication handle
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
    call mpi_bcast( INPUT_OPT%MAX_DIAG, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%MAX_TRCS, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%MAX_MEMB, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%MAX_FAMS, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%MAX_DEP,  1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! SIMULATION MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%NYMDb, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%NHMSb, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%NYMDe, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%NHMSe, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%RUN_DIR, len(INPUT_OPT%RUN_DIR), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%IN_RST_FILE, len(INPUT_OPT%IN_RST_FILE), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LSVGLB, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%OUT_RST_FILE, len(INPUT_OPT%OUT_RST_FILE), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%DATA_DIR, len(INPUT_OPT%DATA_DIR), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%GCAP_DIR, len(INPUT_OPT%GCAP_DIR), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%GEOS_4_DIR, len(INPUT_OPT%GEOS_4_DIR), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%GEOS_5_DIR, len(INPUT_OPT%GEOS_5_DIR), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%GEOS_57_DIR, len(INPUT_OPT%GEOS_57_DIR), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%MERRA_DIR, len(INPUT_OPT%MERRA_DIR), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%DATA_DIR_1x1, len(INPUT_OPT%DATA_DIR_1x1), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TEMP_DIR, len(INPUT_OPT%TEMP_DIR), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LUNZIP, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LWAIT, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LVARTROP, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%NESTED_I0, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%NESTED_J0, 1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! TRACER MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%N_TRACERS, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ID_TRACER(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TRACER_NAME(:), (255)*INPUT_OPT%MAX_TRCS, mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TRACER_MW_G, INPUT_OPT%MAX_TRCS, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TRACER_MW_KG, INPUT_OPT%MAX_TRCS, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TCVV, INPUT_OPT%MAX_TRCS, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%XNUMOL, INPUT_OPT%MAX_TRCS, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TRACER_N_CONST(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TRACER_CONST, len(INPUT_OPT%TRACER_CONST)*INPUT_OPT%MAX_TRCS*INPUT_OPT%MAX_MEMB, mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TRACER_COEFF, INPUT_OPT%MAX_TRCS*INPUT_OPT%MAX_MEMB, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ID_EMITTED(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LSPLIT, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ITS_A_RnPbBe_SIM, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ITS_A_CH3I_SIM, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ITS_A_FULLCHEM_SIM, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ITS_A_HCN_SIM, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ITS_A_TAGOX_SIM, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ITS_A_TAGCO_SIM, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ITS_A_C2H6_SIM, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ITS_A_CH4_SIM, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ITS_AN_AEROSOL_SIM, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ITS_A_MERCURY_SIM, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ITS_A_CO2_SIM, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ITS_A_H2HD_SIM, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ITS_A_POPS_SIM, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ITS_NOT_COPARAM_OR_CH4, 1, mpi_logical, 0, mpiComm, RC )

    !----------------------------------------
    ! AEROSOL MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%LSULF, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LCRYST, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LCARB, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LSOA, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LDUST, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LDEAD, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LSSALT, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LDICARB, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%SALA_REDGE_um(:), 2, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%SALC_REDGE_um(:), 2, mpi_real8, 0, mpiComm, RC )

    !----------------------------------------
    ! EMISSIONS MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%LEMIS, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TS_EMIS, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LANTHRO, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%FSCALYR, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LEMEP, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LBRAVO, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LEDGAR, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LSTREETS, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LCAC, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LNEI05, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LRETRO, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LNEI99, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LICARTT, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LVISTAS, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LBIOFUEL, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LBIOGENIC, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LMEGAN, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LPECCA, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LMEGANMONO, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ISOP_SCALING, 1, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LBIOMASS, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LBBSEA, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LTOMSAI, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LGFED2BB, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%L8DAYBB, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%L3HRBB, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LSYNOPBB, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LGFED3BB, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LDAYBB3, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%L3HRBB3, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LAIRNOX, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LLIGHTNOX, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LOTDLOC, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LSOILNOX, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LFERTILIZERNOX, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%NOx_SCALING, 1, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LEDGARSHIP, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LICOADSSHIP, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LEMEPSHIP, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LSHIPSO2, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LARCSHIP, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LCOOKE, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LAVHRRLAI, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LMODISLAI, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LHIST, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%HISTYR, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LWARWICK_VSLS, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LSSABr2, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LFIX_PBL_BRO, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%Br_SCALING, 1, mpi_real8, 0, mpiComm, RC )

    !----------------------------------------
    ! CO2 MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%LGENFF, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LANNFF, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LMONFF, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LCHEMCO2, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LSEASBB, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LBIODAILY, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LBIODIURNAL, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LBIONETORIG, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LBIONETCLIM, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LOCN1997, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LOCN2009ANN, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LOCN2009MON, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LSHIPEDG, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LSHIPICO, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LPLANE, 1, mpi_logical, 0, mpiComm, RC )

    !----------------------------------------
    ! FUTURE MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%LFUTURE, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%FUTURE_YEAR, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%FUTURE_SCEN, len(INPUT_OPT%FUTURE_SCEN), mpi_character, 0, mpiComm, RC )

    !----------------------------------------
    ! CHEMISTRY MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%LCHEM, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LSCHEM, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LLINOZ, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TS_CHEM, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LSVCSPEC, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LKPP, 1, mpi_logical, 0, mpiComm, RC )

    !----------------------------------------
    ! CHEMISTRY MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%LTRAN, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LMFCT, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LFILL, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TPCORE_IORD, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TPCORE_JORD, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TPCORE_KORD, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TS_DYN, 1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! CONVECTION MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%LCONV, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LTURB, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LNLPBL, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TS_CONV, 1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! DEPOSITION MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%LDRYD, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LWETD, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%USE_OLSON_2001, 1, mpi_logical, 0, mpiComm, RC )

    !----------------------------------------
    ! GAMAP MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%GAMAP_DIAGINFO, len(INPUT_OPT%GAMAP_DIAGINFO), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%GAMAP_TRACERINFO, len(INPUT_OPT%GAMAP_TRACERINFO), mpi_character, 0, mpiComm, RC )

    !----------------------------------------
    ! OUTPUT MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%NJDAY(:), 366, mpi_integer, 0, mpiComm, RC )
    
    !----------------------------------------
    ! DIAGNOSTIC MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%ND01, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD01, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND02, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD02, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND03, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD03, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND04, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD04, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND05, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD05, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND06, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD06, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND07, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD07, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND08, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD08, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND09, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD09, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND10, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD10, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND11, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD11, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND12, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD12, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND13, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD13, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND14, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD14, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND15, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD15, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND16, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD16, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND17, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD17, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND18, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD18, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND19, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD19, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND20, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD20, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND21, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD21, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND22, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD22, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND23, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD23, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND24, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD24, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND25, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD25, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND26, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD26, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND27, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD27, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND28, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD28, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND29, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD29, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND30, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD30, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND31, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD31, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND32, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD32, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND33, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD33, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND34, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD34, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND35, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD35, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND36, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD36, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND37, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD37, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND38, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD38, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND39, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD39, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND40, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD40, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND41, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD41, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND42, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD42, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND43, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD43, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND44, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD44, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND45, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD45, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND46, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD46, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND47, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD47, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND48, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD48, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND49, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD49, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND50, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD50, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD51, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND52, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD52, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND53, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD53, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND54, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD54, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND55, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD55, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND56, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD56, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND57, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD57, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND58, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD58, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND59, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD59, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND60, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD60, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND61, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD61, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND62, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD62, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND63, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD63, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND64, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD64, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND66, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD66, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND67, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD67, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND68, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD68, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND69, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD69, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND70, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD70, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LPRT, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TINDEX(:,:), INPUT_OPT%MAX_DIAG*INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TCOUNT(:)  , INPUT_OPT%MAX_DIAG, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TMAX(:)    , INPUT_OPT%MAX_DIAG, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! PLANEFLIGHT MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%DO_PF, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%PF_IFILE, len(INPUT_OPT%PF_IFILE), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%PF_OFILE, len(INPUT_OPT%PF_OFILE), mpi_character, 0, mpiComm, RC )

    !----------------------------------------
    ! ND48 MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%DO_ND48, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND48_FILE, len(INPUT_OPT%ND48_FILE), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND48_FREQ, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND48_N_STA, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND48_IARR(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND48_JARR(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND48_LARR(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND48_NARR(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! ND49 MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%DO_ND49, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND49_FILE, len(INPUT_OPT%ND49_FILE), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND49_TRACERS(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND49_FREQ, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND49_IMIN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND49_IMAX, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND49_JMIN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND49_JMAX, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND49_LMIN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND49_LMAX, 1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! ND50 MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%DO_ND50, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND50_FILE, len(INPUT_OPT%ND50_FILE), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LND50_HDF, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND50_TRACERS(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND50_IMIN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND50_IMAX, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND50_JMIN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND50_JMAX, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND50_LMIN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND50_LMAX, 1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! ND51 MENU fields
     !----------------------------------------
    call mpi_bcast( INPUT_OPT%DO_ND51, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51_FILE, len(INPUT_OPT%ND51_FILE), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LND51_HDF, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51_TRACERS(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51_HR_WRITE, 1, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51_HR1, 1, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51_HR2, 1, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51_IMIN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51_IMAX, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51_JMIN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51_JMAX, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51_LMIN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51_LMAX, 1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! ND51b MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%DO_ND51b, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51b_FILE, len(INPUT_OPT%ND51b_FILE), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LND51b_HDF, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51b_TRACERS(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51b_HR_WRITE, 1, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51b_HR1, 1, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51b_HR2, 1, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51b_IMIN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51b_IMAX, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51b_JMIN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51b_JMAX, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51b_LMIN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND51b_LMAX, 1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! ND63 MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%DO_ND63, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND63_FILE, len(INPUT_OPT%ND63_FILE), mpi_character, 0, mpiComm, RC )
!    call mpi_bcast( INPUT_OPT%ND63_TRACERS(:), INPUT_OPT%MAX_TRCS, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND63_FREQ, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND63_IMIN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND63_IMAX, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND63_JMIN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND63_JMAX, 1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! PROD LOSS MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%DO_SAVE_PL, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LFAMILY, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ND65, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LD65, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%DO_SAVE_O3, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%NFAM, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%FAM_COEF(:,:), INPUT_OPT%MAX_MEMB*INPUT_OPT%MAX_FAMS, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%FAM_MEMB, len(INPUT_OPT%FAM_MEMB)*INPUT_OPT%MAX_MEMB*INPUT_OPT%MAX_FAMS, mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%FAM_NAME, len(INPUT_OPT%FAM_NAME)*INPUT_OPT%MAX_FAMS, mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%FAM_NMEM, INPUT_OPT%MAX_FAMS, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%FAM_TYPE, len(INPUT_OPT%FAM_TYPE)*INPUT_OPT%MAX_FAMS, mpi_character, 0, mpiComm, RC )

    !----------------------------------------
    ! UNIX CMDS fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%BACKGROUND, len(INPUT_OPT%BACKGROUND), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%REDIRECT, len(INPUT_OPT%REDIRECT), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%REMOVE_CMD, len(INPUT_OPT%REMOVE_CMD), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%SEPARATOR, len(INPUT_OPT%SEPARATOR), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%WILD_CARD, len(INPUT_OPT%WILD_CARD), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%UNZIP_CMD, len(INPUT_OPT%UNZIP_CMD), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%ZIP_SUFFIX, len(INPUT_OPT%ZIP_SUFFIX), mpi_character, 0, mpiComm, RC )

    !----------------------------------------
    ! NESTED GRID MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%LWINDO, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LWINDO2x25, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LWINDO_NA, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TPBC_DIR_NA, len(INPUT_OPT%TPBC_DIR_NA), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LWINDO_EU, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TPBC_DIR_EU, len(INPUT_OPT%TPBC_DIR_EU), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LWINDO_CH, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TPBC_DIR_CH, len(INPUT_OPT%TPBC_DIR_CH), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LWINDO_SE, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TPBC_DIR_SE, len(INPUT_OPT%TPBC_DIR_SE), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LWINDO_CU, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%TPBC_DIR, len(INPUT_OPT%TPBC_DIR), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%NESTED_TS, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%NESTED_I1, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%NESTED_J1, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%NESTED_I2, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%NESTED_J2, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%NESTED_I0W, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%NESTED_J0W, 1, mpi_integer, 0, mpiComm, RC )

    !----------------------------------------
    ! BENCHMARK MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%LSTDRUN, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%STDRUN_INIT_FILE, len(INPUT_OPT%STDRUN_INIT_FILE), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%STDRUN_FINAL_FILE, len(INPUT_OPT%STDRUN_FINAL_FILE), mpi_character, 0, mpiComm, RC )

    !----------------------------------------
    ! ARCHIVED OH MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%OH_DIR, len(INPUT_OPT%OH_DIR), mpi_character, 0, mpiComm, RC )

    !----------------------------------------
    ! O3PL MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%O3PL_DIR, len(INPUT_OPT%O3PL_DIR), mpi_character, 0, mpiComm, RC )

    !----------------------------------------
    ! MERCURY MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%ANTHRO_Hg_YEAR, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%HG_SCENARIO, len(INPUT_OPT%Hg_SCENARIO), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%USE_CHECKS, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LDYNOCEAN, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LPREINDHG, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%Hg_RST_FILE, len(INPUT_OPT%Hg_RST_FILE), mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LGTMM, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%GTMM_RST_FILE, len(INPUT_OPT%GTMM_RST_FILE), mpi_character, 0, mpiComm, RC )

    !----------------------------------------
    ! CH4 MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%LCH4BUD, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LGAO, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LCOL, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LLIV, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LWAST, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LBFCH4, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LRICE, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LOTANT, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LBMCH4, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LWETL, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LSOABS, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%LOTNAT, 1, mpi_logical, 0, mpiComm, RC )

    !----------------------------------------
    ! APM MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%IFNUCL, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%FE0, 1, mpi_real8, 0, mpiComm, RC )

    !----------------------------------------
    ! POPS MENU fields
    !----------------------------------------
    call mpi_bcast( INPUT_OPT%POP_TYPE, 3, mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%CHEM_PROCESS, 1, mpi_logical, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%POP_EMISFILE, 255, mpi_character, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%POP_XMW, 1, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%POP_KOA, 1, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%POP_KBC, 1, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%POP_K_POPG_OH, 1, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%POP_K_POPP_O3A, 1, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%POP_K_POPP_O3B, 1, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%POP_HSTAR, 1, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%POP_DEL_H, 1, mpi_real8, 0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%POP_DEL_Hw, 1, mpi_real8, 0, mpiComm, RC )

    ! MSL
    call mpi_bcast( INPUT_OPT%NUMDEP,   1,         mpi_integer,   0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%NDVZIND,  Input_Opt%MAX_DEP,    mpi_integer,   0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%IDDEP,    NDSTBIN,   mpi_integer,   0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%DUSTREFF, NDSTBIN,   mpi_real8,     0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%DUSTDEN,  NDSTBIN,   mpi_integer,   0, mpiComm, RC )
    call mpi_bcast( INPUT_OPT%DEPNAME,  14*Input_Opt%MAX_DEP, mpi_character, 0, mpiComm, RC )

  END SUBROUTINE GIGC_INPUT_BCAST
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
    SUBROUTINE GIGC_IDT_BCAST( Input_Opt, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod,  ONLY : OptInput
    USE GIGC_Errcode_Mod,    ONLY : GIGC_SUCCESS
    USE TRACERID_MOD
    USE DRYDEP_MOD
    USE M_MPIF
!
! !INPUT PARAMETERS:
!
!!    INTEGER,        INTENT(IN)    :: mpicomm     ! Mpi communication handle
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  04 Jan 2013 - M. Long     - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

    ! SMVGEAR species ID #'s

    call mpi_bcast( IDTSO4, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTNH4, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTNIT, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTBCPI, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTBCPO, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTOCPO, 1, mpi_integer, 0, mpiComm, RC )

    call mpi_bcast( IDTOCPI, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTSOA1, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTSOA2, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTSOA3, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTSOA4, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTSOA5, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTSOAG, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTSOAM, 1, mpi_integer, 0, mpiComm, RC )

    call mpi_bcast( IDTDST1, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTDST2, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTDST3, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTDST4, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTSALA, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTSALC, 1, mpi_integer, 0, mpiComm, RC )

    call mpi_bcast( IDO3    , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDNO2   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDNO3   , 1, mpi_integer, 0, mpiComm, RC ) 
    call mpi_bcast( IDN2O5  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDHNO4  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDOX    , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDNOX   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDHC1   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDNO    , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDHNO2  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDCO    , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDPRPE  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDISOP  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDALK4  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDC3H8  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDPAN   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDGLPAN , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDGPAN  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDPMN   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDPPN   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDHNO3  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDOH    , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDHO2   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDH2O2  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDACET  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDMEK   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDALD2  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDRCHO  , 1, mpi_integer, 0, mpiComm, RC ) 
    call mpi_bcast( IDMVK   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDMACR  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDISN2  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDR4N2  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDCH2O  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDC2H6  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDMP    , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDDMS   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDSO2   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDSO4   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDMSA   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDDRYO3 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDDRYPAN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDDRYNO2, 1, mpi_integer, 0, mpiComm, RC )       
!- Not used (ccc, 01/27/10)
!    call mpi_bcast( IDBENZ  , 1, mpi_integer, 0, mpiComm, RC )
!    call mpi_bcast( IDTOLU  , 1, mpi_integer, 0, mpiComm, RC )
!    call mpi_bcast( IDXYLE  , 1, mpi_integer, 0, mpiComm, RC )
!    call mpi_bcast( IDMONX  , 1, mpi_integer, 0, mpiComm, RC )
!    call mpi_bcast( IDGLYX  , 1, mpi_integer, 0, mpiComm, RC )
!    call mpi_bcast( IDMGLY  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDDRYGLYX, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDDRYMGLY, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDC2H4  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDC2H2  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDMBO   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDGLYC  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDHAC   , 1, mpi_integer, 0, mpiComm, RC )

    call mpi_bcast( IDAPAN  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDENPAN , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDMPAN  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDNIPAN , 1, mpi_integer, 0, mpiComm, RC )

    call mpi_bcast( IDDRYAPAN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDDRYENPAN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDDRYGLPAN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDDRYGPAN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDDRYMPAN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDDRYNIPAN, 1, mpi_integer, 0, mpiComm, RC )
      ! added for aroms (dkh, 10/06/06)  
    call mpi_bcast( IDLBRO2H, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDLBRO2N, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDLTRO2H, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDLTRO2N, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDLXRO2H, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDLXRO2N, 1, mpi_integer, 0, mpiComm, RC )

      ! +++++++++++++++++
      !jpp, 6/13/07,
    call mpi_bcast( IDBr2   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBr    , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBrO   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDHBr   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDHOBr  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBrNO2 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBrNO3 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDCHBr3 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDCH2Br2, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDCH3Br , 1, mpi_integer, 0, mpiComm, RC )
      !jpp, 2/27/08
      !ID's for drydep
    call mpi_bcast( IDDRYHOBr, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDDRYHBr, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDDRYBrNO3, 1, mpi_integer, 0, mpiComm, RC )
      ! +++++++++++++++++

      ! GEOS-CHEM Tracer ID #'s
    call mpi_bcast( IDTNOX  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTOX   , 1, mpi_integer, 0, mpiComm, RC )  
    call mpi_bcast( IDTPAN  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTCO   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTH2   , 1, mpi_integer, 0, mpiComm, RC ) ! (hup, 7/14/2004)
    call mpi_bcast( IDTHD   , 1, mpi_integer, 0, mpiComm, RC ) ! (jaegle, 11/07/2005)
    call mpi_bcast( IDTALK4 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTISOP , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTHNO3 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTH2O2 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTACET , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTMEK  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTALD2 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTRCHO , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTMVK  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTMACR , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTPMN  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTPPN  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTISN2 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTR4N2 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTPRPE , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTC3H8 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTCH2O , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTC2H6 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTN2O5 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTHNO4 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTMP   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTDMS  , 1, mpi_integer, 0, mpiComm, RC ) 
    call mpi_bcast( IDTSO2  , 1, mpi_integer, 0, mpiComm, RC ) 
    call mpi_bcast( IDTMSA  , 1, mpi_integer, 0, mpiComm, RC ) 
    call mpi_bcast( IDTNH3  , 1, mpi_integer, 0, mpiComm, RC ) 
    call mpi_bcast( IDTAS   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTAHS  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTNH4aq, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTLET  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTSO4aq, 1, mpi_integer, 0, mpiComm, RC )

    call mpi_bcast( IDTOxStrt, 1, mpi_integer, 0, mpiComm, RC )
      
    call mpi_bcast( IDTALPH , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTLIMO , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTALCO , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTSOG1 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTSOG2 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTSOG3 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTSOG4 , 1, mpi_integer, 0, mpiComm, RC )
      ! (hotp 5/25/09)
    call mpi_bcast( IDTSOG5 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTSOA5 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTRN   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTPB   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTBE7  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTGLYX , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTMGLY , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTBENZ , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTTOLU , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTXYLE , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTMONX , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTC2H4 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTC2H2 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTMBO  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTGLYC , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTHAC  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTAPAN , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTENPAN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTGLPAN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTGPAN , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTMPAN , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTNIPAN, 1, mpi_integer, 0, mpiComm, RC )
      ! +++++++++++++++++
      ! jpp for bromine
      ! 6/5/09
    call mpi_bcast( IDTBr2   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTBr    , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTBrO   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTHBr   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTHOBr  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTBrNO2 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTBrNO3 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTCHBr3 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTCH2Br2, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTCH3Br , 1, mpi_integer, 0, mpiComm, RC )
      ! +++++++++++++++++


      !(fp, 6/09)
    call mpi_bcast( IDTISOPN , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTAP    , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTPROPNN, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTMOBA  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTMMN   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTRIP   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTIEPOX , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTMAP   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDTPYPAN , 1, mpi_integer, 0, mpiComm, RC )
 
      ! GEOS-CHEM Emission ID #'s
    call mpi_bcast( NEMANTHRO, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( NEMBIOG , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDENOX  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDEOX   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDECO   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDEPRPE , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDEC3H8 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDEALK4 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDEC2H6 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDEACET , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDEMEK  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDEALD2 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDEISOP , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDECH2O , 1, mpi_integer, 0, mpiComm, RC ) 
    call mpi_bcast( IDEHNO3 , 1, mpi_integer, 0, mpiComm, RC ) !phs (3/4/08)
    call mpi_bcast( IDEGLYX , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDEMGLY , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDEBENZ , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDETOLU , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDEXYLE , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDEMONX , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDEC2H4 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDEC2H2 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDEMBO  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDEGLYC , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDEHAC  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDECHBr3, 1, mpi_integer, 0, mpiComm, RC ) ! jpp, 8/30/07
    call mpi_bcast( IDECH2Br2,1, mpi_integer, 0, mpiComm, RC ) ! jpp, 9/16/07
    call mpi_bcast( IDEBr2  , 1, mpi_integer, 0, mpiComm, RC ) ! jpp, 3/2/10

      ! GEOS-CHEM Biofuel ID #'s
    call mpi_bcast( IDBFNOX , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBFCO  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBFALK4, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBFACET, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBFMEK , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBFALD2, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBFPRPE, 1, mpi_integer, 0, mpiComm, RC ) 
    call mpi_bcast( IDBFC3H8, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBFCH2O, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBFC2H6, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBFBENZ, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBFTOLU, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBFXYLE, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBFC2H2, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBFC2H4, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBFGLYC, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBFGLYX, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBFMGLY, 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBFHAC , 1, mpi_integer, 0, mpiComm, RC )

      ! Initialize IDBs for BIOMASS BURNING (hotp 8/24/09)
    call mpi_bcast( IDBNOX  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBCO   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBALK4 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBACET , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBMEK  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBALD2 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBPRPE , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBC3H8 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBCH2O , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBC2H6 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBSO2  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBNH3  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBBC   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBOC   , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBCO2  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBTOLU , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBBENZ , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBXYLE , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBC2H2 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBC2H4 , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBHAC  , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBGLYC , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBGLYX , 1, mpi_integer, 0, mpiComm, RC )
    call mpi_bcast( IDBMGLY , 1, mpi_integer, 0, mpiComm, RC )
!kjw
    call mpi_bcast( IDBCH4  , 1, mpi_integer, 0, mpiComm, RC )
!kjw

      !-----------------------------------
      ! Initialize tagged Hg index arrays
      !-----------------------------------
      IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

         ! Initialize category flags
!- eds 8/31/10 -------------------------------------------------------
       call mpi_bcast( ID_Hg_tot, 1, mpi_integer, 0, mpiComm, RC )
!       call mpi_bcast( ID_Hg_na, 1, mpi_integer, 0, mpiComm, RC )
!       call mpi_bcast( ID_Hg_eu, 1, mpi_integer, 0, mpiComm, RC )
!       call mpi_bcast( ID_Hg_as, 1, mpi_integer, 0, mpiComm, RC )
!       call mpi_bcast( ID_Hg_rw, 1, mpi_integer, 0, mpiComm, RC )
!       call mpi_bcast( ID_Hg_oc, 1, mpi_integer, 0, mpiComm, RC )
!       call mpi_bcast( ID_Hg_ln, 1, mpi_integer, 0, mpiComm, RC )
!       call mpi_bcast( ID_Hg_nt, 1, mpi_integer, 0, mpiComm, RC )

       call mpi_bcast( ID_Hg_can, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_usa, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_cam, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_sam, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_waf, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_eaf, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_saf, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_naf, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_eur, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_eeu, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_mde, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_sov, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_sas, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_eas, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_sea, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_jpn, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_oce, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_so,  1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_bb,  1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_geo, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_atl, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_nat, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_sat, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_npa, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_arc, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_ant, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_ocn, 1, mpi_integer, 0, mpiComm, RC )
       call mpi_bcast( ID_Hg_str, 1, mpi_integer, 0, mpiComm, RC )
    ENDIF
!------------------------------------------------------------------------
    RC = GIGC_SUCCESS

  END SUBROUTINE GIGC_IDT_BCAST

   SUBROUTINE GIGC_READER_BCAST( RC )

     USE GIGC_Errcode_Mod,    ONLY : GIGC_SUCCESS
     USE M_MPIF
     USE COMODE_LOOP_MOD

     INTEGER,        INTENT(OUT)   :: RC          ! Success or failure

     call mpi_bcast( IFSOLVE,   1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( IFURBAN,   1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( IFTROP,    1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( IFSTRAT,   1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( KULOOP,    1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( LYOUT,     1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( LXOUT,     1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( LZOUT,     1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( NREAD,     1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( IPRATES,   1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( IPREADER,  1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( IOSPEC,    1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( IOREAC,    1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( APGASA,    1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( APGASB,    1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( APGASC,    1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( APGASD,    1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( APGASE,    1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( APGASF,    1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( APGASG,    1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( APGASH,    1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( IFREORD,   1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( FRACDEC,   1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( PLOURB,    1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( PLOTROP,   1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( HMAXNIT,   1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( NTLOOP,    1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( NVERT,     1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( NLAYER,    1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( LXOUT,     1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( LYOUT,     1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( LZOUT,     1, mpi_integer, 0, mpiComm, RC )
     call mpi_bcast( KGLC,      1, mpi_integer, 0, mpiComm, RC )

     RC = GIGC_SUCCESS

   END SUBROUTINE GIGC_READER_BCAST

   SUBROUTINE GIGC_READCHEM_BCAST( C1,    CSTRAT, CTROPL, CTROPS, CURBAN, &
                                   ININT, IORD,   NCOF,   RC              )

     USE GIGC_Errcode_Mod,    ONLY : GIGC_SUCCESS
     USE M_MPIF
     USE COMODE_LOOP_MOD

     REAL*8,  INTENT(INOUT) :: C1, CSTRAT, CTROPL, CTROPS, CURBAN
     INTEGER, INTENT(INOUT) :: ININT(10), IORD, NCOF

     INTEGER,        INTENT(OUT)   :: RC          ! Success or failure

!     Globchem.dat
     call mpi_bcast( C1,     1,        mpi_real8,     0, mpiComm, RC ) ! IN
     call mpi_bcast( CSTRAT, 1,        mpi_real8,     0, mpiComm, RC ) ! IN
     call mpi_bcast( CTROPL, 1,        mpi_real8,     0, mpiComm, RC ) ! IN
     call mpi_bcast( CTROPS, 1,        mpi_real8,     0, mpiComm, RC ) ! IN
     call mpi_bcast( CURBAN, 1,        mpi_real8,     0, mpiComm, RC ) ! IN
     call mpi_bcast( ININT,  10,       mpi_integer,   0, mpiComm, RC ) ! IN
     call mpi_bcast( DINP,   4,        mpi_character, 0, mpiComm, RC )
     call mpi_bcast( IORD,   1,        mpi_integer,   0, mpiComm, RC ) ! IN
     call mpi_bcast( NCOF,   1,        mpi_integer,   0, mpiComm, RC ) ! IN
     call mpi_bcast( SPECL,  MXCOF*2,  mpi_character, 0, mpiComm, RC )
     call mpi_bcast( JST,    14,       mpi_character, 0, mpiComm, RC ) 
     call mpi_bcast( ARRT,   MXCOF,    mpi_real8,     0, mpiComm, RC )
     call mpi_bcast( BRRT,   MXCOF,    mpi_real8,     0, mpiComm, RC )
     call mpi_bcast( KCRRT,  MXCOF,    mpi_integer,   0, mpiComm, RC )
     call mpi_bcast( RINP,   IMISC*14, mpi_character, 0, mpiComm, RC )
     call mpi_bcast( PINP,   IMISC,    mpi_real8,     0, mpiComm, RC )
     call mpi_bcast( XINP,   IMISC*14, mpi_character, 0, mpiComm, RC )

!     Chemga.dat
!      READ(I,*)
!      READ(I,610) ASTKCF
!      READ(I,*)
!      READ(I,610)
!      READ(I,*)

     RC = GIGC_SUCCESS

   END SUBROUTINE GIGC_READCHEM_BCAST
   
   SUBROUTINE GIGC_BCAST_CHAR( VAL, SIZE, RC )
     USE GIGC_Errcode_Mod, ONLY : GIGC_SUCCESS
     USE M_MPIF

     CHARACTER(*), INTENT(INOUT) :: VAL(:)
     INTEGER,      INTENT(IN)    :: SIZE
     INTEGER                     :: RC

     RC = GIGC_SUCCESS

     call mpi_bcast( VAL, SIZE, mpi_character, 0, mpiComm, RC )

   END SUBROUTINE GIGC_BCAST_CHAR

   SUBROUTINE GIGC_BCAST_INT( VAL, SIZE, RC )
     USE GIGC_Errcode_Mod, ONLY : GIGC_SUCCESS
     USE M_MPIF

     INTEGER, INTENT(IN)         :: SIZE
     INTEGER, INTENT(INOUT)      :: VAL(SIZE)
     INTEGER                     :: RC

     RC = GIGC_SUCCESS

     call mpi_bcast( VAL, SIZE, mpi_integer, 0, mpiComm, RC )

   END SUBROUTINE GIGC_BCAST_INT

   SUBROUTINE GIGC_BCAST_REAL8( VAL, SIZE, RC )
     USE GIGC_Errcode_Mod, ONLY : GIGC_SUCCESS
     USE M_MPIF

     INTEGER, INTENT(IN)         :: SIZE
     REAL*8, INTENT(INOUT)       :: VAL(SIZE)
     INTEGER                     :: RC

     RC = GIGC_SUCCESS

     call mpi_bcast( VAL, SIZE, mpi_real8, 0, mpiComm, RC )

   END SUBROUTINE GIGC_BCAST_REAL8

!EOC
END MODULE GIGC_MPI_WRAP
#endif
