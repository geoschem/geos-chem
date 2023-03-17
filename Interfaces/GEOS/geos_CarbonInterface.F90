#include "MAPL_Generic.h"
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: geos_CarbonInterface
!
! !DESCRIPTION: Module to handle carbon procedures specific to GEOS, such as
!  importing CO2 from GOCART or producing CO from CO2 photolysis. 
!\\
!\\
! !INTERFACE:
!
MODULE GEOS_CarbonInterface
!
! !USES:
!
  ! MAPL/ESMF
  USE ESMF     
  USE MAPL_Mod 
  USE PHYSCONSTANTS
  USE ESMF_CFIOFileMOD
  USE MAPL_CFIOMOD

  ! GEOS-Chem
  USE Precision_Mod
  USE ErrCode_Mod                                    ! Error numbers
  USE Input_Opt_Mod,         ONLY : OptInput
  USE State_Chm_Mod,         ONLY : ChmState         ! Chemistry State obj
  USE State_Met_Mod,         ONLY : MetState         ! Meteorology State obj
  USE State_Diag_Mod,        ONLY : DgnState         ! Diagnostics State obj
  USE State_Grid_Mod,        ONLY : GrdState         ! Grid State obj

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC   :: GEOS_CarbonSetServices
  PUBLIC   :: GEOS_CarbonInit
  PUBLIC   :: GEOS_CarbonRunPhoto
  PUBLIC   :: GEOS_CarbonGetConc 

! !REVISION HISTORY:
!  12 Jan 2023 - C. Keller - initial version (from StratChem/Carbon_GridComp)
!  See https://github.com/geoschem/geos-chem for full history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_CarbonSetServices
!
! !DESCRIPTION: Set the necessary services 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS_CarbonSetServices( GC, CF, RC ) 
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT), TARGET :: GC     ! Ref to this GridComp
    TYPE(ESMF_Config),   INTENT(INOUT)         :: CF        ! GEOSCHEM*.rc
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(INOUT)         :: RC       ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  12 Jan 2023 - C. Keller   - Initial version (refactored Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    CHARACTER(LEN=*), PARAMETER  :: myname = 'GEOS_CarbonSetServices'
    CHARACTER(LEN=*), PARAMETER  :: Iam = myname    
    CHARACTER(LEN=ESMF_MAXSTR)   :: ImpCO2name
    INTEGER                      :: DoIt
    INTEGER                      :: STATUS

    ! If enabled, create import field 
    CALL ESMF_ConfigGetAttribute( CF, DoIt, Label="Import_CO2_from_GOCART:", Default=0, __RC__ )
    IF ( DoIt == 1 ) THEN
       CALL ESMF_ConfigGetAttribute( CF, ImpCO2name, Label="GOCART_CO2_FieldName:", Default="GOCART_CO2", __RC__ )
       call MAPL_AddImportSpec(GC,                      &
            SHORT_NAME         = TRIM(ImpCO2name),      &
            LONG_NAME          = 'CO2_mixing_ratio',    &
            UNITS              = 'v/v_total_air',       &   ! correct?!
            DIMS               = MAPL_DimsHorzVert,     &
            VLOCATION          = MAPL_VLocationCenter,  &
            RC=STATUS  )
       _VERIFY(STATUS)
    ENDIF

    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE GEOS_CarbonSetServices
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_CarbonInit
!
! !DESCRIPTION: Initialization routine 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS_CarbonInit( GC, CF, State_Chm, State_Grid, RC ) 
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT), TARGET :: GC     ! Ref to this GridComp
    TYPE(ESMF_Config),   INTENT(INOUT)         :: CF        ! GEOSCHEM*.rc
    TYPE(ChmState),      INTENT(INOUT)         :: State_Chm
    TYPE(GrdState),      INTENT(INOUT)         :: State_Grid
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(INOUT)         :: RC       ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  12 Jan 2023 - C. Keller   - Initial version (refactored Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    CHARACTER(LEN=*), PARAMETER  :: myname = 'GEOS_CarbonInit'
    CHARACTER(LEN=*), PARAMETER  :: Iam = myname    
    INTEGER                      :: DoIt, km
    INTEGER                      :: STATUS
    CHARACTER(LEN=ESMF_MAXSTR)   :: fnphoto
    CHARACTER(LEN=ESMF_MAXSTR)   :: ImpCO2name

    !=======================================================================
    ! GEOS_CarbonInit starts here
    !=======================================================================

    ! Import CO2 from GOCART 
    ! ----------------------
    CALL ESMF_ConfigGetAttribute( CF, DoIt, Label="Import_CO2_from_GOCART:", Default=0, __RC__ )
    State_Chm%CO2fromGOCART = ( DoIt == 1 )
    IF ( State_Chm%CO2fromGOCART ) THEN
       CALL ESMF_ConfigGetAttribute( CF, ImpCO2name, Label="GOCART_CO2_FieldName:", Default="GOCART_CO2", __RC__ )
       State_Chm%ImpCO2name = TRIM(ImpCO2name)
       IF ( MAPL_am_I_Root() ) WRITE(*,*) 'Will get CO2 from import field '//TRIM(ImpCO2name)
    ENDIF

    ! CO2 photolysis
    ! --------------

    ! Initialize
    State_Chm%numphoto = 0
    km = State_Grid%NZ

    ! Check if we want to do this
    CALL ESMF_ConfigGetAttribute( CF, DoIt, Label="CO_production_from_CO2_photolysis:", Default=0, __RC__ )
    IF ( DoIt == 1 ) THEN
       CALL ESMF_ConfigGetAttribute( CF, fnphoto, Label="CO2photolysisFile:", Default="please/provide//file.nc", __RC__ )
       IF ( MAPL_am_I_Root() ) WRITE(*,*) 'CO2 photolysis enabled, read photolysis tables: '//TRIM(fnphoto)
       CALL readPhotTables(trim(fnphoto), RC)
       VERIFY_(RC)
       State_Chm%numphoto = 55
    ENDIF

    _RETURN(ESMF_SUCCESS)

CONTAINS

   SUBROUTINE readPhotTables(fileName, rc)

   IMPLICIT NONE

!  Read tables for photolysis in GOCART ... from a NetCDF file
!
!  Input parameters:
!
   CHARACTER(LEN=*), INTENT(IN) :: fileName
!
!  Output parameters:
!
   INTEGER, INTENT(OUT) :: rc
!
!  Restrictions:
!  ASSERT that the number of pressure layers in the dataset equals km.
!
!  REVISION HISTORY:
!  Nielsen     11 May 2012: First crack.
!  Weir        29 Jan 2021: Pilferd from StratChem
!  Keller      12 Jan 2023: Making it worse by also putting it into GCC 
!-----------------------------------------------------------------------

  CHARACTER(LEN=ESMF_MAXSTR) :: Iam = "GCC::readPhotTables"

  TYPE(ESMF_VM) :: vm

  INTEGER :: comm, info, unit, status
  INTEGER :: dimid, i, n

  INTEGER :: length

  INTEGER, PARAMETER :: nD = 7
  CHARACTER(LEN=ESMF_MAXSTR) :: dimName(nD)= (/"nsza  ", "numO3 ", "layers", &
                                               "nlam  ", "nts   ", "nxdo  ", "aqsize" /)

  INTEGER, PARAMETER :: nV = 7
  CHARACTER(LEN=ESMF_MAXSTR) :: varName(nV)= (/"sza    ", &
                        "lambda ", "O3TAB  ",  "SDAT   ", &
                        "O2JDAT ", "XTAB   ",  "CH2O_AQ" /)
  rc = 0

! Grab the virtual machine
! ------------------------
  CALL ESMF_VMGetCurrent(vm, RC=status)
  VERIFY_(status)

  CALL ESMF_VMGet(vm, MPICOMMUNICATOR=comm, rc=status)
  VERIFY_(status)

#ifdef H5_HAVE_PARALLEL

  CALL MPI_Info_create(info, status)
  VERIFY_(status)
  CALL MPI_Info_set(info, "romio_cb_read", "automatic", status)
  VERIFY_(status)

#ifdef NETCDF_NEED_NF_MPIIO
  status = NF_OPEN_PAR(TRIM(fileName), IOR(NF_NOWRITE,NF_MPIIO), comm, info, unit)
#else
  status = NF_OPEN_PAR(TRIM(fileName), NF_NOWRITE, comm, info, unit)
#endif

#else

  IF(MAPL_AM_I_ROOT(vm)) THEN
   status = NF_OPEN(TRIM(fileName), NF_NOWRITE, unit)

#endif

   IF(status /= NF_NOERR) THEN
    PRINT *,'Error opening file ',TRIM(fileName), status
    PRINT *, NF_STRERROR(status)
    VERIFY_(status)
   END IF

   DO i = 1,nD

    status = NF_INQ_DIMID(unit, TRIM(dimName(i)), dimid)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error inquiring dimension ID for ", TRIM(dimName(i)), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

    status = NF_INQ_DIMLEN(unit, dimid, n)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error inquiring  dimension length for ", TRIM(dimName(i)), status
     PRINT *, NF_STRERROR(status)
    END IF

    SELECT CASE (i)
     CASE (1)
      State_Chm%nsza = n
     CASE (2)
      State_Chm%numO3 = n
     CASE (3)
      ASSERT_(n == km)
     CASE (4)
      State_Chm%nlam = n
     CASE (5)
      State_Chm%nts = n
     CASE (6)
      State_Chm%nxdo = n
     CASE (7)
      State_Chm%aqsize = n
     CASE DEFAULT
    END SELECT

   END DO

#ifndef H5_HAVE_PARALLEL

  END IF ! MAPL_AM_I_ROOT

  CALL MAPL_CommsBcast(vm, State_Chm%nsza, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, State_Chm%numO3, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, State_Chm%nlam, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, State_Chm%nts, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, State_Chm%nxdo, 1, 0, RC=status)
  VERIFY_(status)
  CALL MAPL_CommsBcast(vm, State_Chm%aqSize, 1, 0, RC=status)
  VERIFY_(status)

#endif

  ALLOCATE(State_Chm%sdat(State_Chm%nsza,State_Chm%numo3,km,State_Chm%nlam), STAT=status)
  VERIFY_(status)
  ALLOCATE(State_Chm%o2jdat(State_Chm%nsza,State_Chm%numo3,km), STAT=status)
  VERIFY_(status)
  ALLOCATE(State_Chm%o3_tab(State_Chm%numo3,km), STAT=status)
  VERIFY_(status)
  ALLOCATE(State_Chm%xtab(State_Chm%nlam,State_Chm%nxdo,State_Chm%nts), STAT=status)
  VERIFY_(status)
  ALLOCATE(State_Chm%sza_tab(State_Chm%nsza), STAT=status)
  VERIFY_(status)
  ALLOCATE(State_Chm%CH2O_aq(State_Chm%aqSize), STAT=status)
  VERIFY_(status)
  ALLOCATE(State_Chm%rlam(State_Chm%nlam), STAT=status)
  VERIFY_(status)

#ifndef H5_HAVE_PARALLEL

  IF(MAPL_AM_I_ROOT()) THEN

#endif

   DO i = 1,nV

    status = NF_INQ_VARID(unit, TRIM(varName(i)), n)
    IF(status /= NF_NOERR) THEN
     PRINT *,"Error getting varid for ", TRIM(varName(i)), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

    SELECT CASE (i)
     CASE (1)
      status = NF_GET_VAR_REAL(unit, n, State_Chm%sza_tab)
     CASE (2)
      status = NF_GET_VAR_REAL(unit, n, State_Chm%rlam)
     CASE (3)
      status = NF_GET_VAR_REAL(unit, n, State_Chm%o3_tab)
     CASE (4)
      status = NF_GET_VAR_REAL(unit, n, State_Chm%sdat)
     CASE (5)
      status = NF_GET_VAR_REAL(unit, n, State_Chm%o2jdat)
     CASE (6)
      status = NF_GET_VAR_REAL(unit, n, State_Chm%xtab)
     CASE (7)
      status = NF_GET_VAR_REAL(unit, n, State_Chm%CH2O_aq)
     CASE DEFAULT
    END SELECT

    IF(status /= NF_NOERR) THEN
     PRINT *,"Error getting values for ", TRIM(varName(i)), status
     PRINT *, NF_STRERROR(status)
     VERIFY_(status)
    END IF

   END DO

#ifdef H5_HAVE_PARALLEL

   CALL MPI_Info_free(info, status)
   VERIFY_(status)

#else

  END IF ! MAPL_AM_I_ROOT

  length = SIZE(State_Chm%sza_tab)
  CALL MPI_Bcast(State_Chm%sza_tab, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  length = SIZE(State_Chm%rlam)
  CALL MPI_Bcast(State_Chm%rlam, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  length = SIZE(State_Chm%o3_tab)
  CALL MPI_Bcast(State_Chm%o3_tab, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  length = SIZE(State_Chm%sdat)
  CALL MPI_Bcast(State_Chm%sdat, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  length = SIZE(State_Chm%o2jdat)
  CALL MPI_Bcast(State_Chm%o2jdat, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  length = SIZE(State_Chm%xtab)
  CALL MPI_Bcast(State_Chm%xtab, length, MPI_REAL, 0, comm, status)
  VERIFY_(status)

  CALL MAPL_CommsBcast(vm, State_Chm%CH2O_aq, State_Chm%aqsize, 0, RC=status)
  VERIFY_(status)

#endif

  status = NF_CLOSE(unit)
  VERIFY_(status)

  RETURN
  END SUBROUTINE readPhotTables

  END SUBROUTINE GEOS_CarbonInit
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_CarbonGetConc
!
! !DESCRIPTION: Gets the concentrations from the import state 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS_CarbonGetConc( Import,    Input_Opt,  State_Chm,  &
                                 State_Met, State_Diag, State_Grid, RC )
!
! !USES:
!
  USE UnitConv_Mod,          ONLY : Convert_Spc_Units
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_State),    INTENT(INOUT)         :: Import     ! Import State
    TYPE(OptInput),      INTENT(IN)            :: Input_Opt  ! Input Options object
    TYPE(ChmState),      INTENT(INOUT)         :: State_Chm  ! Chemistry object 
    TYPE(MetState),      INTENT(INOUT)         :: State_Met  ! Met object
    TYPE(DgnState),      INTENT(INOUT)         :: State_Diag ! Diagnostics object 
    TYPE(GrdState),      INTENT(INOUT)         :: State_Grid ! Grid object 
    INTEGER,             INTENT(INOUT)         :: RC         ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  12 Jan 2023 - C. Keller   - Initial version (refactored Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    CHARACTER(LEN=*), PARAMETER  :: myname = 'GEOS_CarbonGetConc'
    CHARACTER(LEN=*), PARAMETER  :: Iam = myname    
    CHARACTER(LEN=63)            :: OrigUnit
    INTEGER                      :: I, LM, indCO2, STATUS
    REAL, POINTER                :: CO2(:,:,:) => null()
    REAL, PARAMETER              :: MWCO2 = 44.01 ! everybody knows this

    !=======================================================================
    ! GEOS_CarbonGetConc starts here
    !=======================================================================

    IF ( State_Chm%CO2fromGOCART ) THEN

       ! Make sure concentrations are in kg/kg total (this should already be the case) 
       CALL Convert_Spc_Units( Input_Opt,         State_Chm,     State_Grid, &
                               State_Met,         'kg/kg total', RC,        &
                               OrigUnit=OrigUnit                             )
       ASSERT_(RC==GC_SUCCESS)

       ! Get index
       indCO2  = -1
       DO I = 1, State_Chm%nSpecies
          IF ( TRIM(State_Chm%SpcData(I)%Info%Name) == "CO2"  ) THEN
             indCO2 = I 
             EXIT
          ENDIF
       ENDDO
       ASSERT_(indCO2 > 0  )

       ! Get CO2 field via import. This is expected in v/v total!! 
       CALL MAPL_GetPointer ( Import, CO2, TRIM(State_Chm%ImpCO2name), __RC__ )

       ! Pass to GEOS-Chem, flip in vertical and convert v/v to kg/kg
       LM = State_Grid%NZ
       State_Chm%Species(indCO2)%Conc(:,:,:) = CO2(:,:,LM:1:-1) * ( MWCO2 / MAPL_AIRMW )

       ! testing only
       !!!IF ( MAPL_am_I_Root() ) WRITE(*,*) 'Got CO2 from GOCART: ',SUM(CO2)

       ! Convert species back to original units
       CALL Convert_Spc_Units( Input_Opt, State_Chm,  State_Grid, State_Met, &
                               OrigUnit,  RC )
       ASSERT_( RC == GC_SUCCESS )

    ENDIF

    _RETURN(ESMF_SUCCESS)
    END SUBROUTINE GEOS_CarbonGetConc

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_CarbonRunPhoto
!
! !DESCRIPTION: Calculates the photolysis rates (following StratChem) 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS_CarbonRunPhoto( Input_Opt,  State_Chm,  State_Met,  &
                                  State_Diag, State_Grid, RC )
!
! !USES:
!
  USE TIME_MOD,              ONLY : GET_TS_CHEM
  USE UnitConv_Mod,          ONLY : Convert_Spc_Units
  USE ERROR_MOD,             ONLY : SAFE_DIV
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)            :: Input_Opt  ! Input Options object
    TYPE(ChmState),      INTENT(INOUT)         :: State_Chm  ! Chemistry object 
    TYPE(MetState),      INTENT(INOUT)         :: State_Met  ! Met object
    TYPE(DgnState),      INTENT(INOUT)         :: State_Diag ! Diagnostics object 
    TYPE(GrdState),      INTENT(INOUT)         :: State_Grid ! Grid object 
    INTEGER,             INTENT(INOUT)         :: RC         ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  12 Jan 2023 - C. Keller   - Initial version (refactored Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    CHARACTER(LEN=*), PARAMETER  :: myname = 'GEOS_CarbonRunPhoto'
    CHARACTER(LEN=*), PARAMETER  :: Iam = myname    

    CHARACTER(LEN=63)            :: OrigUnit
    REAL, ALLOCATABLE            :: aj(:)
    INTEGER                      :: I, J, L, LM, STATUS
    INTEGER                      :: indCO, indCO2, indO3
    REAL                         :: sza, o3col, press, temp
    REAL                         :: photj, dCOphot, tsChem 
    REAL(fp)                     :: COinit, COpost, CO2conc 

    !=======================================================================
    ! GEOS_CarbonRunPhoto starts here
    !=======================================================================

    ! Only do this if active
    IF ( State_Chm%numphoto > 0) THEN

       ! Convert to molec/cm3 units are molec/cm3 
       CALL Convert_Spc_Units( Input_Opt,         State_Chm,    State_Grid, &
                               State_Met,         'molec/cm3'  , RC,        &
                               OrigUnit=OrigUnit                             )
       ASSERT_(RC==GC_SUCCESS)

       ! Chemistry time step in secods
       tsChem  = GET_TS_CHEM()

       ! Get index for species: need CO, O3, and CO2 
       indCO  = -1
       indO3  = -1
       indCO2 = -1
       DO I = 1, State_Chm%nSpecies
          IF ( TRIM(State_Chm%SpcData(I)%Info%Name) == "CO"  ) indCO  = I 
          IF ( TRIM(State_Chm%SpcData(I)%Info%Name) == "O3"  ) indO3  = I 
          IF ( TRIM(State_Chm%SpcData(I)%Info%Name) == "CO2" ) indCO2 = I 
       ENDDO
       ASSERT_(indCO  > 0  )
       ASSERT_(indO3  > 0  )
       ASSERT_(indCO2 > 0  )

       ! Allocate local arrays
       allocate(aj(State_Chm%numphoto), STAT=RC)
       VERIFY_(RC)

       ! Loop over entire atmosphere
       LM = State_Grid%NZ
       DO L = 1, LM 
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Solar Zenith Angle (radians)
          SZA = 0.
          IF ( State_Met%SUNCOSmid(I,J)<=1.0 ) SZA = ACOS(State_Met%SUNCOSmid(I,J))

          ! Overhead ozone (molec cm-2)
          O3col = SUM(   State_Chm%Species(indO3)%Conc(I,J,L:LM) &
                       * State_Met%BXHEIGHT(I,J,L:LM)*100.0       )

          ! Pressure (hPa) and Temperature (K) 
          Press = State_Met%PMID_DRY(I,J,L)
          Temp  = State_Met%T(I,J,L)

          ! Calculate photolysis rates (s-1) a la StratChem 
          call jcalc4(L, SZA, O3col, Press, Temp, aj, State_Chm)
          photJ = aj(12)

          ! Get CO species concentration (molec cm-3)
          COinit  = State_Chm%Species(indCO)%Conc(I,J,L)

          ! Get CO2 (molec cm-3)
          CO2conc = State_Chm%Species(indCO2)%Conc(I,J,L)

          ! production rate
          dCOphot = photJ*CO2conc

          ! Update CO concentration 
          COpost = COinit + tsChem*dCOphot

          ! Add back to concentration array 
          State_Chm%Species(indCO)%Conc(I,J,L) = COpost

          ! Add to diagnostics if requested 
          IF ( State_Diag%Archive_CO2photrate ) THEN
             State_Diag%CO2photrate(I,J,L) = photJ 
          ENDIF

          IF ( State_Diag%Archive_COincCO2phot ) THEN
             State_Diag%COincCO2phot(I,J,L) = SAFE_DIV( COpost, COinit, 1.0_fp, 1.0_fp ) - 1.0_fp
          ENDIF
       ENDDO
       ENDDO
       ENDDO

       ! Convert species back to original units
       CALL Convert_Spc_Units( Input_Opt, State_Chm,  State_Grid, State_Met, &
                               OrigUnit,  RC )
       ASSERT_( RC == GC_SUCCESS )

       ! Cleanup
       IF ( ALLOCATED(aj) ) DEALLOCATE(aj)
    ENDIF

    _RETURN(ESMF_SUCCESS)

CONTAINS

   SUBROUTINE interp_s(k,sza,o3column,s,jo2,State_Chm)
! ----------------------------------------------------------------------------
! NAME:
!   interp_s
!
! PURPOSE:
!   Interpolate S values for each wavelength in table to specified O3
!   column and zenith angle
!
! INPUTS:
!   k         Current layer number
!   szaRad    Solar zenith angle [radians]
!   o3column  Overhead o3 column value [cm^{-2}]
!   State_Chm      The GOCART::CO grid component, which contains
!     sza_tab Solar zenith angle table
!     o3_tab  Overhead O3 values table
!     sdat    Radiative source function 
!     o2jdat  Table of J(O2) values
!
! OUTPUTS:
!   s         S value for each wavelength at current k, interpolated to
!               the given o3column and sza
!   jo2       J(O2) values interpolated as above
!
! 
! PROCEDURE:
!   Bi-linear interpolation, for sza > 94 s=0, for O3 out of range use min/max
!
! MODIFICATION HISTORY: 
!   25 Aug 1993  Kawa
!   10 Jul 1996  Kawa    For 28 levels and to handle J(O2) separately
!   11 May 2012  Nielsen Accomodation for GEOS-5 FV cubed release
!   30 Jan 2021  Weir    Copied from StratChem
! ----------------------------------------------------------------------------

   IMPLICIT NONE

   TYPE(ChmState), INTENT(IN) :: State_Chm   ! Grid Component

   INTEGER, INTENT(IN) :: k
   REAL, INTENT(IN) :: sza, o3column
   REAL, INTENT(OUT) :: s(State_Chm%nlam), jo2

   INTEGER :: ijj, ik, ikk, ikkm, il, is
   REAL :: omt, omu, t, u
   REAL, PARAMETER :: PI = 3.14159265

! For each input solar zenith angle, find the first element of State_Chm%sza_tab that 
! is greater.  Use this element and previous one to determine the interpolated value.
! -----------------------------------------------------------------------------------
   DO is = 1,State_Chm%nsza
      ijj = is
      IF(State_Chm%sza_tab(is) > sza) EXIT
   ENDDO

! Zenith angle test       
! -----------------
   IF(sza > State_Chm%sza_tab(State_Chm%nsza)) THEN
!     Cell is dark, set s and jo2=0        
!     -----------------------------
      s(1:State_Chm%nlam) = 0.
      jo2 = 0.
   ELSE
!     Cell is illuminated     
!     -------------------
      t = (sza-State_Chm%sza_tab(ijj-1))/(State_Chm%sza_tab(ijj)-State_Chm%sza_tab(ijj-1))
      omt = 1.-t

! For each overhead O3 column, find the first element in State_Chm%o3_tab that is
! greater. Use this element and previous one to determine the interpolated value.
! -------------------------------------------------------------------------------
      DO is = 1,State_Chm%numo3
         ikk = is
         IF(State_Chm%o3_tab(is,k) > o3column) EXIT
      ENDDO

      ikkm = ikk-1
      IF(ikk > 1 .AND. o3column <= State_Chm%o3_tab(State_Chm%numo3,k)) THEN
         u = (o3column-State_Chm%o3_tab(ikkm,k))/(State_Chm%o3_tab(ikk,k)-State_Chm%o3_tab(ikkm,k))
         omu = 1.-u

! Do bilinear interpolation for each wavelength.
! ----------------------------------------------
         DO il = 1,State_Chm%nlam
            s(il) = omt*omu*State_Chm%sdat(ijj-1,ikkm,k,il)+t*omu*State_Chm%sdat(ijj,ikkm,k,il)+ &
                    t*u*State_Chm%sdat(ijj,ikk,k,il)+omt*u*State_Chm%sdat(ijj-1,ikk,k,il)
         ENDDO
         jo2 = omt*omu*State_Chm%o2jdat(ijj-1,ikkm,k)+t*omu*State_Chm%o2jdat(ijj,ikkm,k)+ &
               t*u*State_Chm%o2jdat(ijj,ikk,k)+omt*u*State_Chm%o2jdat(ijj-1,ikk,k)

! Extrapolate ahead of table
! --------------------------
      ELSE IF (ikk == 1) THEN
         DO il = 1,State_Chm%nlam
            s(il) = omt*State_Chm%sdat(ijj-1,1,k,il)+t*State_Chm%sdat(ijj,1,k,il)
         ENDDO
         jo2 = omt*State_Chm%o2jdat(ijj-1,1,k)+t*State_Chm%o2jdat(ijj,1,k)

! Extrapolate beyond table
! ------------------------
      ELSE
         DO il = 1,State_Chm%nlam
            s(il) = omt*State_Chm%sdat(ijj-1,State_Chm%numo3,k,il)+t*State_Chm%sdat(ijj,State_Chm%numo3,k,il)
         END DO
         jo2 = omt*State_Chm%o2jdat(ijj-1,State_Chm%numo3,k)+t*State_Chm%o2jdat(ijj,State_Chm%numo3,k)
      ENDIF
   ENDIF

   RETURN
   END SUBROUTINE interp_s

   SUBROUTINE jcalc4(k,szan,o3column,press,kel,aj,State_Chm)
! ---------------------------------------------------------------------------------
! NAME: jcalc4
! PURPOSE:
!   Calculate photolysis rates
! INPUT:
!   k         Current layer number
!   levels    Number of layers
!   szan      Solar zenith angle (radians)
!   o3column  Overhead O3 values
!   press     Mid-layer pressure (hPa)
!   kel       Mid-layer temperature (K)
! OUTPUT:
!   aj        Array of photolysis rates
! RESTRICTIONS:
!   Currently set up for 23-J set (see var State_Chm%nxdo)
! REQUIRED ROUTINES:
!   interp_s
! MODIFICATION HISTORY: 
!   26 Aug 1993 Kawa    Created
!   23 Nov 1993 Kawa    Remade xtab to do multiplication by solar flux beforehand 
!                        and removed inputs.
!   25 Feb 1994         Add 3 additional Js, incl N2O
!   18 Sep 1995         Add 2 additional Js, up to 22, and do CH2O special
!   13 May 1996 Crum    Removed fossils, move toward Fortran 90
!   10 Jul 1996         Modified to handle J(O2) separately and use 28 levels
!    1 Apr 2009 Nielsen GEOS-5 form with standardized SC_GridComp interface.
!    1 Jun 2009 Nielsen Updated to JPL 2006
!   12 Dec 2010 Nielsen Updated to JPL 2010 following Luke Oman's testing.
!   11 May 2012 Nielsen Accomodation for GEOS-5 FV cubed release
!    3 Jun 2015 Liang   Updated to the new 50-slot table with addition of halons,
!                       HCFCs, and 5 VSLSs
!                       numphoto is now updated to 52
!   30 Jan 2021 Weir    Copied from StratChem
!
! WARNING: Photolysis reaction rate numbers 38-42 are calculated in MESO_PHOT.
! ---------------------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, PARAMETER :: DBL = KIND(0.00D+00)

   TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Grid Component

   INTEGER, INTENT(IN) :: k
   REAL, INTENT(IN) :: szan, o3column, press, kel
!  REAL(KIND=DBL), INTENT(OUT) :: aj(State_Chm%numphoto)
!  bweir: demoted to single
   REAL, INTENT(OUT) :: aj(State_Chm%numphoto)

   INTEGER :: ilam,indt,ix

   REAL :: alpha300, alphat, jo2, rjs(State_Chm%nxdo), q1, q2, r1mq1
   REAL :: s(State_Chm%nlam), sx(2,State_Chm%nlam), tfac, wvl

! Start with a clean slate
! ------------------------
   aj(1:State_Chm%numphoto) = 0.

! Interpolate radiative flux function values to model conditions
! --------------------------------------------------------------
   CALL interp_s(k,szan,o3column,s,jo2,State_Chm)
   indt = kel-148.5
   indt = MAX(1,indt)
   indt = MIN(indt,200)

! cakelle2: comment stuff not needed for CO2 photolysis...

!!!! Preliminaries for CH2O quantum yield dependence on m, T, wavelength
!!!! -------------------------------------------------------------------
!!!   tfac = (kel-80.0)/80.0
!!!
!!!   DO ilam=1,State_Chm%nlam
!!!      ZeroS: IF(s(ilam) == 0.) THEN
!!!         sx(1,ilam) = 0.00
!!!         sx(2,ilam) = 0.00
!!!      ELSE
!!!
!!!         wvl = State_Chm%rlam(ilam)*0.10
!!!
!!!         IF(wvl < 250.00) THEN
!!!            q1 = 0.24
!!!         ELSE IF(wvl >= 339.00) THEN
!!!            q1 = 0.00
!!!         ELSE
!!!            q1 = State_Chm%CH2O_aq(1) + State_Chm%CH2O_aq(2)*wvl         + &
!!!                                   State_Chm%CH2O_aq(3)*wvl*wvl     + &
!!!                                   State_Chm%CH2O_aq(4)*wvl*wvl*wvl + &
!!!                                   State_Chm%CH2O_aq(5)*wvl*wvl*wvl*wvl
!!!         ENDIF
!!!
!!!         r1mq1 = 1./(1.-q1)
!!!
!!!         IF(wvl < 330.00) THEN
!!!            q2 = State_Chm%xtab(ilam,22,indt)
!!!         ELSE IF(wvl > 360.00) THEN
!!!            q2 = 0.00
!!!         ELSE
!!!            alpha300 = 1.00E-03*(1./State_Chm%xtab(ilam,22,1)-r1mq1)
!!!            alphat = alpha300*(1.+0.05*(wvl-329.)*((300.-kel)/80.))
!!!            q2 = 1.00/(r1mq1+alphat*press)
!!!         ENDIF
!!!
!!!         IF(wvl .LT. 250.00) q2=0.5
!!!
!!!         sx(2,ilam) = s(ilam)*State_Chm%xtab(ilam,21,indt)*q2
!!!         sx(1,ilam) = s(ilam)*State_Chm%xtab(ilam,21,indt)*q1
!!!      ENDIF ZeroS
!!!   ENDDO
!!!
!!!! J(BrONO2) through J(OCLO)
!!!! -------------------------
!!!   DO ix=1,14
!!!      rjs(ix) = 0.
!!!
!!!      DO ilam=1,State_Chm%nlam
!!!         rjs(ix) = rjs(ix)+s(ilam)*State_Chm%xtab(ilam,ix,indt)
!!!      ENDDO
!!!   ENDDO
!!!
!!!! J(O2)
!!!! -----
!!!   rjs(15) = jo2
!!!
!!!! J(O3_O1D) through J(N2O)
!!!! ------------------------
!!!   DO ix=16,20
!!!      rjs(ix) = 0.
!!!
!!!      DO ilam=1,State_Chm%nlam
!!!         rjs(ix) = rjs(ix)+s(ilam)*State_Chm%xtab(ilam,ix,indt)
!!!      ENDDO
!!!   ENDDO
!!!
!!!! J(CH2O)
!!!! -------
!!!   rjs(21) = 0.
!!!   rjs(22) = 0.
!!!   DO ilam=1,State_Chm%nlam
!!!      rjs(21) = rjs(21)+sx(1,ilam)
!!!      rjs(22) = rjs(22)+sx(2,ilam)
!!!   ENDDO

! J(CO2 -> CO + O) through xH1211
! -------------------------------
!!!   DO ix=23,State_Chm%nxdo
      ix=23
      rjs(ix) = 0.

      DO ilam=1,State_Chm%nlam
         rjs(ix) = rjs(ix)+s(ilam)*State_Chm%xtab(ilam,ix,indt)
      ENDDO
!!!   ENDDO

! ---------------------------------------------------------------
! Order photolysis rates to match order in full chemistry model.  
! Sort rjs into CTM photolysis rate array, aj.  Order of rjs:
!
!  1-J(BrONO2)
!  2-J(BrO)
!  3-J(Cl2O2)
!  4-J(ClONO2)
!  5-J(H2O2)
!  6-J(HCl)
!  7-J(HNO3)
!  8-J(HO2NO2)
!  9-J(HOCl)
! 10-J(N2O5)
! 11-J(NO2)
! 12-J(NO3_NO)
! 13-J(NO3_NO2)
! 14-J(OClO)
! 15-J(O2)
! 16-J(O3_O1D)
! 17-J(O3_3P)
! 18-J(HOBr)
! 19-J(CH3OOH)
! 20-J(N2O)
! 21-J(CH2O_HCO)
! 22-J(CH2O_CO)
! 23-J(CO2 -> CO + O)
! 24-xCFC-11
! 25-xCFC-12
! 26-xCCl4
! 27-xCH3CCl3
! 28-xHCFC-22
! 29-xCFC-113
! 30-xCH3Cl
! 31-xCH3Br
! 32-xH1301
! 33-xH1211 
! 34-xH1202
! 35-xH2402
! 36-xCHBr3
! 37-xCH2Br2
! 38-xCH2ClBr
! 39-xCHClBr2
! 40-xCHCl2Br
! 41-xHCFC-141b
! 42-xHCFC-142b
! 43-xCFC-114 
! 44-xCFC-115
! 45-xOCS
! 46-
! 47-
! 48-
! 49-
! 50-
! ---------------------------------------------------------------
! ---------------------------------------------------------------
! Solar cycle goes here when ready  
!     aj( 1) = rjs(15)*State_Chm%s_cycle(3,State_Chm%iscyr)
! ----------------------------------------------------------------
!!!   aj( 1) = rjs(15)
!!!   aj( 2) = rjs(16)
!!!   aj( 3) = rjs(17)
!!!! H2O
!!!! ---
!!!   aj( 4) = 0.
!!!   aj( 5) = rjs(13)
!!!   aj( 6) = rjs(7)
!!!   aj( 7) = rjs(11)
!!!   aj( 8) = rjs(5)
!!!   aj( 9) = rjs(10)
!!!   aj(10) = rjs(21)
!!!   aj(11) = rjs(22)
   aj(12) = rjs(23)
!!!   aj(13) = rjs(19)
!!!   aj(14) = rjs(20)
!!!   aj(15) = rjs(4)
!!!   aj(16) = 0.
!!!   aj(17) = rjs(12)
!!!   aj(18) = rjs(6)
!!!   aj(19) = 0.
!!!
!!!! CH3Br(20) H1301(21) H12_24(22)
!!!! ------------------------------
!!!   aj(20) = rjs(31)
!!!   aj(21) = rjs(32)
!!!   aj(22) = rjs(33)
!!!   aj(23) = rjs(9)
!!!   aj(24) = rjs(8)
!!!   aj(25) = rjs(18)
!!!   aj(26) = 0.
!!!   aj(27) = rjs(2)
!!!   aj(28) = rjs(1)
!!!
!!!! F11(29) F12(30) CCl4(31) CHCCl3(32) HCFC(33) F113(34) CH3Cl(35)
!!!! ---------------------------------------------------------------
!!!   aj(29) = rjs(24)
!!!   aj(30) = rjs(25)
!!!   aj(31) = rjs(26)
!!!   aj(32) = rjs(27)
!!!   aj(33) = rjs(28)
!!!   aj(34) = rjs(29)
!!!   aj(35) = rjs(30)
!!!   aj(36) = rjs(3)
!!!   aj(37) = rjs(14)
!!!
!!!! ------------------------------------------
!!!! WARNING: Photolysis reaction rate
!!!! numbers 38-42 are calculated in MESO_PHOT.
!!!! ------------------------------------------
!!!! Add aj(43) which is J(Cl2O2) for partitioning but not Ox loss 
!!!! which is aj(36). In lookup table J(Cl2O2) is J*qy where qy is 0.8 
!!!! so multiply by 1.25 to equal J and used in part.F and partest.F
!!!
!!!   aj(43) = rjs(3)*1.25
!!!
!!!! QingLiang -- 06/03/2015
!!!! CHBr3(44) CH2Br2(45) CH2BrCl(46) CHBrCl2(47) CHBr2Cl(48)
!!!   aj(44) = rjs(36)
!!!   aj(45) = rjs(37)
!!!   aj(46) = rjs(38)
!!!   aj(47) = rjs(39)
!!!   aj(48) = rjs(40)
!!!
!!!! QingLiang -- 06/03/2015
!!!! Add two new halons: H-1202 (49) H2402 (50) 
!!!! and two new HCFCs: HCFC-141b (51) HCFC-142b (52) 
!!!   aj(49) = rjs(34)
!!!   aj(50) = rjs(35)
!!!   aj(51) = rjs(41)
!!!   aj(52) = rjs(42)
!!!
!!!! QingLiang -- 02/05/2016
!!!! Add CFC-114 and CFC-115
!!!! Add OCS for GOCART module
!!!   aj(53) = rjs(43)
!!!   aj(54) = rjs(44)
!!!   aj(55) = rjs(45)
!!!!  aj(53) = rjs(34)
!!!!  aj(54) = rjs(34)
!!!!  aj(55) = rjs(34)

   RETURN
   END SUBROUTINE jcalc4

  END SUBROUTINE GEOS_CarbonRunPhoto
!EOC
END MODULE GEOS_CarbonInterface
