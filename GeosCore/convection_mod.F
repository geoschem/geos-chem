!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: convection_mod.F90
!
! !DESCRIPTION: Module CONVECTION\_MOD contains routines which select the 
!  proper convection code for GEOS-3, GEOS-4, GEOS-5, MERRA, or GCAP met 
!  field data sets. 
!\\
!\\
! !INTERFACE: 
!
      MODULE CONVECTION_MOD
! 
! !USES:
!
      USE ErrCode_Mod
      USE PhysConstants    ! Physical constants
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: DO_CONVECTION
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: DO_GEOS4_CONVECT
      PRIVATE :: DO_GCAP_CONVECT
      PRIVATE :: NFCLDMX
      PRIVATE :: DO_MERRA_CONVECTION
!
! !REVISION HISTORY:
!  27 Jan 2004 - R. Yantosca - Initial version
!  (1 ) Contains new updates for GEOS-4/fvDAS convection.  Also now references
!        "error_mod.f".  Now make F in routine NFCLDMX a 4-D array to avoid
!        memory problems on the Altix. (bmy, 1/27/04)
!  (2 ) Bug fix: Now pass NTRACE elements of TCVV to FVDAS_CONVECT in routine 
!        DO_CONVECTION (bmy, 2/23/04)  
!  (3 ) Now references "logical_mod.f" and "tracer_mod.f" (bmy, 7/20/04)
!  (4 ) Now also references "ocean_mercury_mod.f" and "tracerid_mod.f" 
!        (sas, bmy, 1/19/05)
!  (5 ) Now added routines DO_GEOS4_CONVECT and DO_GCAP_CONVECT by breaking 
!        off code from DO_CONVECTION, in order to implement GCAP convection
!        in a much cleaner way. (swu, bmy, 5/25/05)
!  (6 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (7 ) Shut off scavenging in shallow convection for GCAP (swu, bmy, 11/1/05)
!  (8 ) Modified for tagged Hg simulation (cdh, bmy, 1/6/06)
!  (9 ) Bug fix: now only call ADD_Hg2_WD if LDYNOCEAN=T (phs, 2/8/07)
!  (10) Fix for GEOS-5 met fields in routine NFCLDMX (swu, 8/15/07)
!  (11) Resize DTCSUM array in NFCLDMX to save memory (bmy, 1/31/08)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Treat MERRA in the same way as for GEOS-5
!  29 Sep 2010 - R. Yantosca - Added modifications for MERRA
!  05 Oct 2010 - R. Yantosca - Added ND14 and ND38 diagnostics to 
!                              DO_MERRA_CONVECTION routine
!  16 Aug 2011 - J. Fisher   - Minor bug fixes in DO_MERRA_CONVECTION
!  15 Feb 2011 - R. Yantosca - Add modifications for APM from G. Luo
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  02 Mar 2012 - R. Yantosca - Now reference the new grid_mod.F90
!  22 Oct 2012 - R. Yantosca - Now reference Headers/gigc_errcode_mod.F90
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  07 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  23 Jun 2015 - E. Lundgren - Convert tracer units from v/v dry air to
!                              kg/kg total air for convection
!  23 Jun 2016 - R. Yantosca - Remove references to APM code; it is no longer
!                              compatible with the FlexChem implementation
!  06 Jul 2016 - E. Lundgren - Spc units are now kg/kg dry air in convection
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_convection
!
! !DESCRIPTION: Subroutine DO\_CONVECTION calls the appropriate convection 
!  driver program for different met field data sets.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_CONVECTION( am_I_Root, Input_Opt,
     &                          State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : CONVFLUP
      USE ERROR_MOD,          ONLY : GC_Error
      USE UnitConv_Mod
#if defined( BPCH_DIAG )
      USE DIAG_MOD,           ONLY : AD38
#endif
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE GC_GRID_MOD,        ONLY : GET_AREA_M2
      USE Input_Opt_Mod,      ONLY : OptInput
      USE Species_Mod,        ONLY : Species
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_DYN
      USE TIME_MOD,           ONLY : GET_TS_CONV
      USE WETSCAV_MOD,        ONLY : COMPUTE_F
      USE WETSCAV_MOD,        ONLY : H2O2s
      USE WETSCAV_MOD,        ONLY : SO2s
#if defined( NC_DIAG )
      USE HCO_ERROR_MOD
      USE HCO_INTERFACE_MOD
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE HCO_DIAGN_MOD,      ONLY : Diagn_Update

#endif
#if defined( USE_TEND )
      USE TENDENCIES_MOD
      USE State_Chm_Mod,      ONLY : Ind_
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!  25 May 2005 - S. Wu       - Initial version
!  08 Feb 2007 - R. Yantosca - Now reference "CMN_SIZE".  Now references 
!                              CLDMAS, CMFMC, DTRAIN from "dao_mod.f" so that 
!                              we can pass either GEOS-5 or GEOS-3 meteorology
!                              to NFCLDMX. 
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Treat MERRA in the same way as for GEOS-5
!  29 Sep 2010 - R. Yantosca - Now call DO_MERRA_CONVECTION for MERRA met
!  05 Oct 2010 - R. Yantosca - Now attach diagnostics to MERRA conv routine
!  06 Oct 2010 - R. Yantosca - Parallelized call to DO_MERRA_CONVECTION
!  15 Oct 2010 - H. Amos     - Now get BXHEIGHT, T from dao_mod.f
!  15 Oct 2010 - R. Yantosca - Now get LDYNOCEAN from logical_mod.f
!  15 Oct 2010 - R. Yantosca - Now get ITS_A_MERCURY_SIM from tracer_mod.f
!  15 Oct 2010 - R. Yantosca - Now get IDTHg2, IDTHgP from tracerid_mod.f
!  15 Oct 2010 - R. Yantosca - Now get H2O2s, SO2s from wetscav_mod.f
!  15 Oct 2010 - H. Amos     - Now pass BXHEIGHT, T, to DO_MERRA_CONVECTION
!  15 Oct 2010 - R. Yantosca - Now pass H2O2s, SO2s to DO_MERRA_CONVECTION
!  15 Feb 2011 - R. Yantosca - Add modifications for APM from G. Luo
!  29 Aug 2011 - R. Yantosca - Bug fix: reposition #if defined(APM) statement
!  09 Feb 2012 - R. Yantosca - For GEOS-5.7, PFICU and PFLCU fields are
!                              defined on level edges.  Pass the top edge
!                              of each level to DO_MERRA_CONVECTION
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  21 Jun 2012 - R. Yantosca - Now use poiners to pass array slices to routines
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  04 Feb 2013 - S. Kim      - Bug fix: H2O2s, SO2s, STT are not in State_Met
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  31 May 2013 - R. Yantosca - Now pass objects to NFCLDMX
!  03 Jun 2013 - R. Yantosca - Bug fix: pass State_Chm to DO_MERRA_CONVECTION
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  17 Apr 2014 - R. Yantosca - Speed up MERRA and GEOS-FP convection by
!                              adding !$OMP+SCHEDULE( DYNAMIC )
!  17 Apr 2014 - R. Yantosca - Test if (ND14>0) and (ND38>0) outside OMP loop
!  18 Apr 2014 - R. Yantosca - Now use a pointer to pass a slice of the
!                              F array to subroutine COMPUTE_F
!  18 Apr 2014 - R. Yantosca - Clean up call to DO_MERRA_CONVECTION, remove
!                              stuff that was leftover from the column code
!  18 Apr 2014 - R. Yantosca - Now use proper # of tracers for APM in the
!                              call to DO_MERRA_CONVECTION (MERRA & GEOS-FP)
!  25 Jun 2014 - R. Yantosca - Now pass Input_Opt to COMPUTE_F
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE with State_Met%PEDGE.
!                              Remove dependency on pressure_mod
!  28 Apr 2015 - E. Lundgren - Consolidate NFCLDMX arguments to remove 
!                              passing State_Met array slice
!  15 Jun 2015 - E. Lundgren - Now use kg/kg total air as tracer units
!                              (previously v/v)
!  10 Aug 2015 - E. Lundgren - Incoming tracer units are now kg/kg dry air
!  11 Aug 2015 - R. Yantosca - Added support for MERRA2 data
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  06 Jul 2016 - R. Yantosca - Use State_Chm%Map_WetDep to get species ID's
!  07 Jul 2016 - R. Yantosca - Now dimension DIAG14 and DIAG38 for the #
!                              of soluble species (State_Chm%nWetDep)
!  07 Jul 2016 - E. Lundgren - Now use spc kg/kg dry instead of kg/kg total
!   8 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER           :: N, NA, nAdvect, NW
#if defined( USE_TEND )
      REAL(fp)          :: DT_TEND
#endif

      !=================================================================
      ! MERRA or GEOS-FP met fields definitions
      !=================================================================
#if defined( MERRA ) || defined( GEOS_FP ) || defined( MERRA2 )
      
      ! Scalars          
      INTEGER           :: I, J, L, NN, TS_DYN
      REAL(fp)          :: AREA_M2, DT
      LOGICAL           :: DO_ND14
      LOGICAL           :: DO_ND38
      
      ! Arrays        
      REAL(fp)          :: PEDGE  (          LLPAR+1                  )
      INTEGER           :: ISOL   (                  State_Chm%nAdvect)
      REAL(fp)          :: DIAG14 (            LLPAR,State_Chm%nWetDep)
      REAL(fp)          :: DIAG38 (            LLPAR,State_Chm%nWetDep)
      REAL(fp)          :: F      (            LLPAR,State_Chm%nAdvect)
      REAL(fp), TARGET  :: FSOL   (IIPAR,JJPAR,LLPAR,State_Chm%nAdvect)

      ! Pointers
      REAL(fp), POINTER :: p_FSOL (:,:,:)
#endif

#if defined ( NC_DIAG )
      ! For netCDF diagnostics
      INTEGER                :: ND, D, HCRC
      CHARACTER(LEN=60)      :: DiagnName
      REAL(fp), ALLOCATABLE  :: BEFORE(:,:,:,:)
      REAL(fp), POINTER      :: TMP3D(:,:,:)
      REAL(fp), POINTER      :: TMP2D(:,:  )
      REAL(fp)               :: DTIME

      ! Objects
      TYPE(Species), POINTER :: SpcInfo
#endif

      !=================================================================
      ! DO_CONVECT begins here!
      !=================================================================
#if defined ( NETCDF )
      ! Initialize pointers
      Tmp2d => NULL()
      Tmp3d => NULL()
#endif

#if defined( USE_TEND )
      !=================================================================
      ! Archive species concentrations for tendencies (ckeller,7/15/2015)
      !=================================================================
      CALL TEND_STAGE1( am_I_Root, Input_Opt, State_Met, 
     &                  State_Chm, 'CONV', RC )
#endif
 
      ! Assume success
      RC = GC_SUCCESS

      ! Number of advected species
      nAdvect = State_Chm%nAdvect

      ! For netCDF diagnostics, archive old species concentration. At the 
      ! moment, the new ND38 diagnostics only returns the total column 
      ! loss due to convection. For these values, we don't need to 
      ! write the diagnostics within the convection schemes but can
      ! simply take the difference between the original and updated 
      ! species field column.
#if defined( NC_DIAG )
      IF ( Input_Opt%ND38 > 0 ) THEN
         ALLOCATE( BEFORE(IIPAR,JJPAR,LLPAR,State_Chm%nSpecies) )
         BEFORE = State_Chm%Species
      ENDIF
#endif

      !=================================================================
      ! Call met field specific convection routines 
      !=================================================================

#if   defined( GCAP ) 

      !=================================================================
      ! GCAP met fields
      !=================================================================

      ! Assume success
      RC =  GC_SUCCESS

      ! Call GCAP driver routine
      CALL DO_GCAP_CONVECT
     &   ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

#elif defined( GEOS_4 )

      !=================================================================
      ! GEOS-4 met fields
      !=================================================================

      ! Assume success
      RC =  GC_SUCCESS

      ! Call GEOS-4 driver routine
      CALL DO_GEOS4_CONVECT
     &   ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

#elif defined( GEOS_5 )

      !=================================================================
      ! GEOS-5 met fields
      !=================================================================

      ! Assume success
      RC =  GC_SUCCESS

      ! Call the S-J Lin convection routine for GEOS-5
      CALL NFCLDMX( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

#elif defined( MERRA ) || defined( GEOS_FP ) || defined( MERRA2 )

      !=================================================================
      ! GEOS-FP, MERRA, or MERRA-2 met fields
      !=================================================================

      ! Initialize
      RC        =  GC_SUCCESS                   ! Assume success
      DO_ND14   = ( Input_Opt%ND14 > 0 )          ! Should we use ND14 diat?
      DO_ND38   = ( Input_Opt%ND38 > 0 )          ! Should we use ND38 diag?
      TS_DYN    = GET_TS_DYN()                    ! Dynamic timestep [min]
      DT        = DBLE( TS_DYN )                  ! Dynamic timestep [min]
      FSOL      = 0e+0_fp                         ! Zero the FSOL array

      ! Loop over advected species
      DO NA = 1, nAdvect
         
         ! Species ID
         N = State_Chm%Map_Advect(NA)

         ! Now point to a 3D slice of the FSOL array
         p_FSOL => FSOL(:,:,:,NA)

         ! Fraction of soluble species
         CALL COMPUTE_F( am_I_Root, N, p_FSOL, ISOL(NA), 
     &                   Input_Opt, State_Met, State_Chm, RC ) 

         ! Free pointer memory
         p_FSOL => NULL()

      ENDDO

      !-----------------------------------------------------------------
      ! Do convection column by column
      !-----------------------------------------------------------------
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( J,      I,      AREA_M2, L, F  ) 
!$OMP_PRIVATE( DIAG14, DIAG38, RC,      N, NA,    NW )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip buffer zone for nested grids (lzh, 4/1/15)
         IF ( Input_Opt%ITS_A_NESTED_GRID ) THEN
            IF ( J <=         Input_Opt%NESTED_J0W  ) CYCLE
            IF ( J >  JJPAR - Input_Opt%NESTED_J0E  ) CYCLE
            IF ( I <=         Input_Opt%NESTED_I0W  ) CYCLE
            IF ( I >  IIPAR - Input_Opt%NESTED_I0E  ) CYCLE
         ENDIF

         ! Grid box surface area [m2]
         AREA_M2     =  GET_AREA_M2( I, J, 1 ) 

         ! NOTE: For some reason the code chokes when we try
         ! to use a pointer, so we'll use a 2D array here.
         F           =  FSOL(I,J,:,:)
         
         !--------------------------
         ! Do the cloud convection
         !--------------------------
         CALL DO_MERRA_CONVECTION( am_I_Root  = am_I_Root,
     &                             Input_Opt  = Input_Opt,
     &                             State_Met  = State_Met,
     &                             State_Chm  = State_Chm,
     &                             I          = I,
     &                             J          = J,
     &                             AREA_M2    = AREA_M2,
     &                             F          = F,
     &                             TS_DYN     = DT, 
     &                             USE_DIAG14 = DO_ND14,
     &                             DIAG14     = DIAG14,
     &                             USE_DIAG38 = DO_ND38,
     &                             DIAG38     = DIAG38,
     &                             RC         = RC         )

         ! Stop if error is encountered
         IF ( RC /= GC_SUCCESS ) THEN
            CALL GEOS_CHEM_STOP
         ENDIF

         !--------------------------
         ! ND14 diagnostic [kg/s]
         !--------------------------
         IF ( DO_ND14 ) THEN

            ! Loop over only the soluble species
            DO NW = 1, State_Chm%nWetDep

               ! Get the species ID From the wetdep ID
               N = State_Chm%Map_Advect(NA)

               ! Archive diagnostic
               DO L = 1, LD14
                  CONVFLUP(I,J,L,N) = CONVFLUP(I,J,L,N) + DIAG14(L,NW)
               ENDDO
            ENDDO
         ENDIF

#if defined( BPCH_DIAG )
         !--------------------------
         ! ND38 diagnostic [kg/s]
         !--------------------------
         IF ( DO_ND38 ) THEN

            ! Loop over only the soluble species
            DO NW = 1, State_Chm%nWetDep
            DO L  = 1, LD38
               AD38(I,J,L,NW) = AD38(I,J,L,NW) + DIAG38(L,NW)
            ENDDO
            ENDDO
         ENDIF
#endif
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

#endif
#if defined( NC_DIAG )
      !-----------------------------------------------------------------
      ! NetCDF diagnostics (ND38)
      ! Write out the column mass difference before and after 
      ! convection.
      ! Note: The ND38 diagnostics is defined as a 2D diagnostics
      ! (--> diagnostics_mod.F90).
      !-----------------------------------------------------------------
      IF ( Input_Opt%ND38 > 0 ) THEN

         ! Assume success
         HCRC = HCO_SUCCESS

         ! Time step
         DTIME = GET_TS_DYN() * 60e+0_fp

         ! Allocate temporary pointer
         ALLOCATE(TMP3D(IIPAR,JJPAR,LLPAR), TMP2D(IIPAR,JJPAR) )
        
         ! Loop over the number of wet-depositing species
         DO ND = 1, State_Chm%nWetDep

            ! Get the species ID from the wetdep ID
            D = State_Chm%Map_WetDep(ND)

            ! Get information for this species from the species database
            SpcInfo => State_Chm%SpcData(D)%Info

            ! Is this a species selected for this diagnostics type?
            IF ( ANY( Input_Opt%TINDEX(38,:) == D ) ) THEN

               ! Loss rate [kg/kg total air per seconnd]
               TMP3D = ( BEFORE(:,:,:,D) - 
     &                 State_Chm%Species(:,:,:,D)) / DTIME 

               ! Convert to [kg/s]
               TMP3D = TMP3D * State_Met%DELP_DRY(I,J,L) * 
     &                 G0_100 * AREA_M2

               ! Sum over the column
               TMP2D = SUM(TMP3D,DIM=3) 

               ! Update diagnostics based on name
               DiagnName = 'WETSCAV_CONVLOSS_' 
     &                     // TRIM( SpcInfo%Name )

               CALL Diagn_Update ( am_I_Root, HcoState, 
     &                             cName   = TRIM( DiagnName ),
     &                             Array2D = TMP2D,
     &                             COL     = Input_Opt%DIAG_COLLECTION,
     &                             PosOnly = .TRUE.,
     &                             RC      = HCRC )
               IF ( HCRC /= HCO_SUCCESS ) THEN
                  CALL ERROR_STOP ( 'Error updating wetscav conv loss '
     &               // 'diagnostic: ' // TRIM( DiagnName ), 
     &               'DO_CONVECTION (convection_mod.F)' )
               ENDIF
            ENDIF

            ! Free pointer
            SpcInfo => NULL()

         ENDDO

         ! Cleanup
         DEALLOCATE( BEFORE, TMP3D, TMP2D )
      ENDIF
#endif

#if defined( USE_TEND )
      !=================================================================
      ! Calculate tendencies and write to diagnostics 
      ! (ckeller,7/15/2015)
      !=================================================================
      DT_TEND = GET_TS_CONV() * 60.0_fp
      CALL TEND_STAGE2( am_I_Root, Input_Opt, State_Met, 
     &                  State_Chm, 'CONV', DT_TEND, RC )
#endif

      END SUBROUTINE DO_CONVECTION
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_geos4_convect
!
! !DESCRIPTION: Subroutine DO\_GEOS4\_CONVECT is a wrapper for the 
!  GEOS-4/fvDAS convection code.  This was broken off from the old 
!  DO\_CONVECTION routine above.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_GEOS4_CONVECT( am_I_Root, Input_Opt,
     &                             State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD37
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE FVDAS_CONVECT_MOD,  ONLY : INIT_FVDAS_CONVECT, FVDAS_CONVECT
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CONV
      USE WETSCAV_MOD,        ONLY : COMPUTE_F
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
! !REMARKS:
!  As of July 2016, we assume that all of the advected species are listed
!  first in the species database.  This is the easiest way to pass a slab
!  to the TPCORE routine.  This may change in the future. (bmy, 7/13/16)
!                                                                             .
!  %%%% GEOS-4 IS DEPRECATED AND MAY BE REMOVED SOON %%%%
!
! !REVISION HISTORY: 
!  25 May 2005 - S. Wu       - Initial version
!  (1 ) Now use array masks to flip arrays vertically in call to FVDAS_CONVECT
!        (bmy, 5/25/05)
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (3 ) Add a check to set negative values in STT to TINY (ccc, 4/15/09)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  15 Feb 2011 - R. Yantosca - Add modifications for APM from G. Luo
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  14 Apr 2014 - R. Yantosca - Remove array temporary in call to FVDAS_CONVECT
!  18 Apr 2014 - R. Yantosca - Now use a pointer to pass a slice of the
!                              F array to subroutine COMPUTE_F
!  23 Jun 2014 - R. Yantosca - Now pass Input_Opt to CONVTRAN
!  25 Jun 2014 - R. Yantosca - Now pass Input_Opt to COMPUTE_F
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE differences with 
!                              State_Met%DELP. Remove dependency on 
!                              pressure_mod
!  23 Jun 2015 - E. Lundgren - Now use kg/kg total air as tracer units
!                              (previously v/v)
!  22 Apr 2016 - R. Yantosca - Now pass State_Chm to FVDAS_CONVECT_MOD
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  06 Jul 2016 - E. Lundgren - Now use kg/kg dry air as spc units, requiring
!                              use of DELP_DRY instead of DELP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! SAVEd scalars
      LOGICAL, SAVE     :: FIRST = .TRUE.

      ! Scalars
      LOGICAL           :: LPRT
      INTEGER           :: I, ISOL, J, L, L2, N, NSTEP, NA, nAdvect
      INTEGER           :: CONVDT
      REAL(fp)          :: P1, P2, TDT

      ! Arrays
      REAL(fp)          :: RPDEL   (IIPAR,JJPAR,LLPAR                  )
      REAL(fp)          :: DP      (IIPAR,JJPAR,LLPAR                  )
      REAL(fp), TARGET  :: F       (IIPAR,JJPAR,LLPAR,State_Chm%nAdvect)
      INTEGER           :: INDEXSOL(                  State_Chm%nAdvect)

      ! Pointers
      REAL(fp), POINTER :: Spc     (:,:,:,:) ! Ptr to advected species
      REAL(fp), POINTER :: p_HKETA (:,:,:  ) ! Ptr to HKETA  (flipped in vert)
      REAL(fp), POINTER :: p_HKBETA(:,:,:  ) ! Ptr to HKBETA (flipped in vert)
      REAL(fp), POINTER :: p_ZMMU  (:,:,:  ) ! Ptr to ZMMU   (flipped in vert)
      REAL(fp), POINTER :: p_ZMMD  (:,:,:  ) ! Ptr to ZMMD   (flipped in vert)
      REAL(fp), POINTER :: p_ZMEU  (:,:,:  ) ! Ptr to ZMEU   (flipped in vert)
      REAL(fp), POINTER :: p_Spc   (:,:,:,:) ! Ptr to Spc    (flipped in vert)
      REAL(fp), POINTER :: p_F     (:,:,:,:) ! Ptr to F      (flipped in vert)
      REAL(fp), POINTER :: p_FSOL  (:,:,:  ) ! Ptr to F      (for COMPUTE_F)

      !=================================================================
      ! DO_GEOS4_CONVECT begins here!
      !=================================================================
      
      ! Assume success
      RC        =  GC_SUCCESS

      ! Initialize
      LPRT      = Input_Opt%LPRT
      nAdvect   = State_Chm%nAdvect

      ! Convection timestep [s]
      CONVDT = GET_TS_CONV() * 60e+0_fp 
       
      ! NSTEP is the # of internal convection timesteps.  According to
      ! notes in the old convection code, 300s works well. (swu, 12/12/03)
      NSTEP  = CONVDT / 300    
      NSTEP  = MAX( NSTEP, 1 ) 

      ! TIMESTEP*2; will be divided by 2 before passing to CONVTRAN 
      TDT    = DBLE( CONVDT ) * 2.0e+0_fp / DBLE( NSTEP )

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_FVDAS_CONVECT( am_I_Root, State_Met )
         FIRST = .FALSE.
      ENDIF

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_G4_CONV: a INIT_FV' )

      !=================================================================
      ! Before calling convection, compute the fraction of insoluble
      ! species (Finsoluble) lost in updrafts.  Finsoluble = 1-Fsoluble.
      !=================================================================

      ! Loop over only the advected species
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, NA, ISOL, p_FSOL )
!$OMP+SCHEDULE( DYNAMIC )
      DO NA = 1, nAdvect

         ! Point to array slice of F
         p_FSOL => F(:,:,:,NA)

         ! Get fraction of species scavenged and the soluble species 
         ! index (ISOL). For non-soluble species, F=0 and ISOL=0.
         CALL COMPUTE_F( am_I_Root, NA,        p_FSOL,    ISOL, 
     &                   Input_Opt, State_Met, State_Chm, RC ) 
         
         ! Free pointer memory
         p_FSOL => NULL()

         ! Store ISOL in an array for later use
         INDEXSOL(NA) = ISOL

         ! Loop over grid boxes
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! ND37 diagnostic: store fraction of species 
            ! lost in moist convective updrafts ("MC-FRC-$")
            IF ( ND37 > 0 .and. ISOL > 0 .and. L <= LD37 ) THEN
               AD37(I,J,L,ISOL) = AD37(I,J,L,ISOL) + F(I,J,L,NA) 
            ENDIF

            ! GEOS-4 convection routines need the insoluble fraction
            F(I,J,L,NA) = 1e+0_fp - F(I,J,L,NA)
         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
       
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_G4_CONV: a COMPUTE_F' )

      !=================================================================
      ! Compute pressure thickness arrays DP and RPDEL
      ! These arrays are indexed from atm top --> surface
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, L2, P1, P2 )
      DO L = 1, LLPAR

         ! L2 runs from the atm top down to the surface
         L2 = LLPAR - L + 1

         ! Loop over surface grid boxes
         DO J = 1, JJPAR
         DO I = 1, IIPAR
               
            ! DP = Pressure difference between top & bottom edges [Pa]
            ! Now use DELP_DRY, equivalent to delta dry pressure across
            ! the level [hPa] (ewl, 7/6/15)
            DP(I,J,L2) = State_Met%DELP_DRY(I,J,L) * 100.0e+0_fp

            ! RPDEL = reciprocal of DP [1/hPa]
            RPDEL(I,J,L2) = 100.0e+0_fp / DP(I,J,L2) 
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_G4_CONV: a DP, RPDEL' )
   
      !=================================================================
      ! Flip arrays in the vertical and call FVDAS_CONVECT
      !=================================================================

      ! Point to chemical species array  [kg/kg total air]
      !
      ! NOTE: For now, so as to avoid having to rewrite the internals
      ! of the GEOS-4 convection routines, just point to 1:nAdvect entries 
      ! of State_Chm%Species.  This is OK for now, as of July 2016, all of
      ! the advected species are listed first.  This may change in the
      ! future, but we'll worry about that later.  The units of p_SPC
      ! will be converted to [kg/kg moist air] below. (bmy, 7/13/16)
      Spc      => State_Chm%Species( :, :, :,          1:nAdvect )   

      ! Point to the chemical species array [kg/kg moist air]
      ! flipping the arrays in the vertical (bmy, 7/18/16)
      p_Spc    => State_Chm%Species( :, :, LLPAR:1:-1, 1:nAdvect )

      ! Point to various fields and flip in vertical.  This avoids 
      ! array temporaries in the call to FVDAS_CONVECT (bmy, 4/14/14)
      p_HKETA  => State_Met%HKETA  ( :, :, LLPAR:1:-1            )
      p_HKBETA => State_Met%HKBETA ( :, :, LLPAR:1:-1            )
      p_ZMMU   => State_Met%ZMMU   ( :, :, LLPAR:1:-1            )
      p_ZMMD   => State_Met%ZMMD   ( :, :, LLPAR:1:-1            )
      p_ZMEU   => State_Met%ZMEU   ( :, :, LLPAR:1:-1            )
      p_F      => F                ( :, :, LLPAR:1:-1, :         )

      ! Call the fvDAS convection routines (originally from NCAR!)
      CALL FVDAS_CONVECT( TDT,   
     &                    nAdvect, 
     &                    p_Spc,
     &                    RPDEL,         
     &                    p_HKETA,
     &                    p_HKBETA,
     &                    p_ZMMU,
     &                    p_ZMMD,
     &                    p_ZMEU,
     &                    DP,     
     &                    NSTEP,    
     &                    p_F,
     &                    INDEXSOL,
     &                    Input_Opt,
     &                    State_Met,
     &                    State_Chm )

      !==================================================================
      ! Add a check to set negative values in Spc to TINY (ccc, 4/15/09)
      ! Only loop over the advected species (bmy, 7/18/16)
      !==================================================================
!$OMP PARALLEL DO        
!$OMP+DEFAULT( SHARED   )
!$OMP+PRIVATE( NA, N, L, J, I )
      DO NA = 1, nAdvect

         ! Get the species ID from the advected species ID
         N = State_Chm%Map_Advect(NA)
         
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            Spc(I,J,L,N) = MAX(Spc(I,J,L,N),TINY(1e+0_fp))
         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointers
      NULLIFY( Spc,     p_Spc,    p_F                    )
      NULLIFY( p_HKETA, p_HKBETA, p_ZMMU, p_ZMMD, p_ZMEU )

      !### Debug! 
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_G4_CONV: a FVDAS_CONVECT' )

      END SUBROUTINE DO_GEOS4_CONVECT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_gcap_convect
!
! !DESCRIPTION: Subroutine DO\_GCAP\_CONVECT is a wrapper for the GCAP 
!  convection code.  This was broken off from the old DO\_CONVECTION routine 
!  above.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_GCAP_CONVECT( am_I_Root, Input_Opt,
     &                            State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD37
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE GCAP_CONVECT_MOD,   ONLY : GCAP_CONVECT
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CONV
      USE WETSCAV_MOD,        ONLY : COMPUTE_F
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!  25 May 2005 - S. Wu       - Initial version
!  (1 ) Now use array masks to flip arrays vertically in call to GCAP_CONVECT
!        (bmy, 5/25/05)
!  (2 ) Shut off scavenging in shallow convection for GCAP below 700 hPa
!        (swu, bmy, 11/1/05)
!  (3 ) Add a check to set negative values in STT to TINY (ccc, 4/15/09)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  15 Feb 2011 - R. Yantosca - Add modifications for APM from G. Luo
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  04 Nov 2013 - M. Sulprizio- Now use pointer variables to flip met fields in
!                              the vertical
!  15 Apr 2014 - R. Yantosca - Remove array temporaries in call to GCAP_CONVECT
!  18 Apr 2014 - R. Yantosca - Now use a pointer to pass a slice of the
!                              F array to subroutine COMPUTE_F
!  25 Jun 2014 - R. Yantosca - Now pass Input_Opt to COMPUTE_F
!  26 Feb 2015 - E. Lundgren - Replace GET_PCENTER with State_Met%PMID.
!                              Replace PEDGE difference with State_Met%DELP.
!                              Remove dependency on pressure_mod.
!  23 Jun 2015 - E. Lundgren - Now use kg/kg total air as tracer units
!                              (previously v/v)
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  06 Jul 2016 - E. Lundgren - Now use kg/kg dry air as spc units, requiring
!                              use of DELP_DRY instead of DELP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE      :: FIRST = .TRUE.
      INTEGER            :: I, ISOL, J, L, L2, N, NSTEP, NA, nAdvect
      INTEGER            :: CONVDT    
      REAL(fp)           :: DP( IIPAR, JJPAR, LLPAR )
      REAL(fp)           :: P1, P2, TDT   
      REAL(fp)           :: GAINMASS
      INTEGER            :: N_TOT_TRC      
      INTEGER            :: INDEXSOL(State_Chm%nAdvect)
      REAL(fp),  TARGET  :: F(IIPAR,JJPAR,LLPAR,State_Chm%nAdvect)

      ! For fields from Input_Opt
      LOGICAL            :: LPRT

      ! Pointers
      REAL(fp),  POINTER :: Spc        (:,:,:,:) ! Ptr to advected species
      REAL(fp),  POINTER :: p_UPDE     (:,:,:  ) ! Ptr to UPDE     met field
      REAL(fp),  POINTER :: p_DNDE     (:,:,:  ) ! Ptr to DNDE     met field
      REAL(fp),  POINTER :: p_ENTRAIN  (:,:,:  ) ! Ptr to ENTRAIN  met field
      REAL(fp),  POINTER :: p_DETRAINE (:,:,:  ) ! Ptr to DETRAINE met field
      REAL(fp),  POINTER :: p_UPDN     (:,:,:  ) ! Ptr to UPDN     met field
      REAL(fp),  POINTER :: p_DNDN     (:,:,:  ) ! Ptr to DNDN     met field
      REAL(fp),  POINTER :: p_DETRAINN (:,:,:  ) ! Ptr to DETRAINN met field
      REAL(fp),  POINTER :: p_Spc      (:,:,:,:) ! Ptr to Spc (flipped)
      REAL(fp),  POINTER :: p_F        (:,:,:,:) ! Ptr to F   (flipped)
      REAL(fp),  POINTER :: p_FSOL     (:,:,:  ) ! Ptr to F   (for COMPUTE_F)

      !=================================================================
      ! DO_GCAP_CONVECT begins here!
      !=================================================================
      
      ! Initialize
      ! Assume success
      RC       = GC_SUCCESS
      LPRT     = Input_Opt%LPRT
      nAdvect  = State_Chm%nAdvect
      GAINMASS = 0e+0_fp

      ! Convection timestep [s]
      CONVDT   = GET_TS_CONV() * 60e+0_fp 
       
      ! NSTEP is the # of internal convection timesteps.  According to
      ! notes in the old convection code, 300s works well. (swu, 12/12/03)
      NSTEP    = CONVDT / 300    
      NSTEP    = MAX( NSTEP, 1 ) 

      ! TIMESTEP*2; will be divided by 2 before passing to CONVTRAN 
      TDT      = DBLE( CONVDT ) * 2.0e+0_fp / DBLE( NSTEP )
         
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_GCAP_CONV: a INIT_FV' )

      !=================================================================
      ! Before calling convection, compute the fraction of insoluble
      ! species (Finsoluble) lost in updrafts.  Finsoluble = 1-Fsoluble.
      !=================================================================

      ! Loop over only the advected species
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, NA, ISOL, p_FSOL )
!$OMP+SCHEDULE( DYNAMIC )
      DO NA = 1, nAdvect

         ! Point to 3D slice of F
         p_FSOL => F(:,:,:,NA)

         ! Get fraction of species scavenged and the soluble species 
         ! index (ISOL). For non-soluble species, F=0 and ISOL=0.
         CALL COMPUTE_F( am_I_Root, NA,        p_FSOL,    ISOL, 
     &                   Input_Opt, State_Met, State_Chm, RC    ) 
         
         ! Free pointer memory
         p_FSOL => NULL()

         ! Store ISOL in an array for later use
         INDEXSOL(NA) = ISOL

         ! Loop over grid boxes
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Shut off scavenging in shallow convection for GCAP
            ! (swu, bmy, 11/1/05)
            IF ( State_Met%PMID( I, J, L ) > 700e+0_fp ) 
     &            F(I,J,L,NA) = 0e+0_fp

            ! ND37 diagnostic: store fraction of species 
            ! lost in moist convective updrafts ("MC-FRC-$")
            IF ( ND37 > 0 .and. ISOL > 0 .and. L <= LD37 ) THEN
               AD37(I,J,L,ISOL) = AD37(I,J,L,ISOL) + F(I,J,L,NA) 
            ENDIF

            ! GEOS-4 convection routines need the insoluble fraction
            F(I,J,L,NA) = 1e+0_fp - F(I,J,L,NA)
         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
       
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_GCAP_CONV: a COMPUTE_F' )

      !=================================================================
      ! Compute pressure thickness arrays DP and RPDEL
      ! These arrays are indexed from atm top --> surface
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, L2, P1, P2 )
      DO L = 1, LLPAR

         ! L2 runs from the atm top down to the surface
         L2 = LLPAR - L + 1

         ! Loop over surface grid boxes
         DO J = 1, JJPAR
         DO I = 1, IIPAR
               
             ! DP = Pressure difference between top & bottom edges [Pa]
             ! Now use DELP_DRY, equivalent to the delta dry pressure 
             ! across the level [hPa] (ewl, 7/6/15)
             DP(I,J,L2) = State_Met%DELP_DRY(I,J,L) * 100.0e+0_fp
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_GCAP_CONV: a DP, RPDEL' )
   
      !=================================================================
      ! Flip arrays in the vertical and call GCAP_CONVECT
      !=================================================================

      ! Point to chemical species array  [kg/kg total air]
      !
      ! NOTE: For now, so as to avoid having to rewrite the internals
      ! of the GEOS-4 convection routines, just point to 1:nAdvect entries 
      ! of State_Chm%Species.  This is OK for now, as of July 2016, all of
      ! the advected species are listed first.  This may change in the
      ! future, but we'll worry about that later.  The units of p_SPC
      ! will be converted to [kg/kg moist air] below. (bmy, 7/13/16)
      Spc        => State_Chm%Species ( :, :, :,          1:nAdvect )   

      ! Point to the chemical species array [kg/kg moist air]
      ! flipping the arrays in the vertical (bmy, 7/18/16)
      p_Spc      => State_Chm%Species ( :, :, LLPAR:1:-1, 1:nAdvect )

      ! Point to met fields and flip in vertical
      p_UPDE     => State_Met%UPDE    ( :, :, LLPAR:1:-1            )
      p_DNDE     => State_Met%DNDE    ( :, :, LLPAR:1:-1            )
      p_ENTRAIN  => State_Met%ENTRAIN ( :, :, LLPAR:1:-1            )
      p_DETRAINE => State_Met%DETRAINE( :, :, LLPAR:1:-1            )
      p_UPDN     => State_Met%UPDN    ( :, :, LLPAR:1:-1            )
      p_DNDN     => State_Met%DNDN    ( :, :, LLPAR:1:-1            )
      p_DETRAINN => State_Met%DETRAINN( :, :, LLPAR:1:-1            )
      p_F        => F                 ( :, :, LLPAR:1:-1, :         )

      ! Call the GCAP convection routines 
      CALL GCAP_CONVECT( TDT, 
     &                   p_Spc,
     &                   nAdvect,
     &                   DP, ! I think this is the correct way (bmy, 5/25/05)
     &                   NSTEP, 
     &                   p_F,
     &                   INDEXSOL, 
     &                   p_UPDE,
     &                   p_DNDE,
     &                   p_ENTRAIN,
     &                   p_DETRAINE,
     &                   p_UPDN,
     &                   p_DNDN,
     &                   p_DETRAINN )

      !===================================================================
      ! Add a check to set negative values in Spc to TINY (ccc, 4/15/09)
      ! Only loop over the advected species (bmy, 7/5/16)
      !===================================================================
!$OMP PARALLEL DO        
!$OMP+DEFAULT( SHARED   )
!$OMP+PRIVATE( NA, N, L, J, I )
      DO NA = 1, nAdvect

         ! Get the species ID from the advected species ID
         N = State_Chm%Map_Advect(NA)

         ! Loop over grid boxes
         DO L  = 1, LLPAR
         DO J  = 1, JJPAR
         DO I  = 1, IIPAR
            Spc(I,J,L,N) = MAX(Spc(I,J,L,N),TINY(1e+0_fp))
         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( Spc,    p_Spc,  p_F                    )
      NULLIFY( p_UPDE, p_DNDE, p_ENTRAIN,  p_DETRAINE )
      NULLIFY( p_UPDN, p_DNDN, p_DETRAINN             )

      !### Debug! 
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_GCAP_CONV: a GCAP_CONVECT' )

      END SUBROUTINE DO_GCAP_CONVECT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: nfcldmx
!
! !DESCRIPTION: Subroutine NFCLDMX is S-J Lin's cumulus transport module for 
!  3D GSFC-CTM, modified for the GEOS-Chem model.  The "NF" stands for "no 
!  flipping", and denotes that you don't have to flip the species array Q in 
!  the main program before passing it to NFCLDMX.
!\\
!\\
!  NOTE: NFCLDMX can be used with GEOS-1, GEOS-STRAT, and GEOS-3 met fields.
!  For GEOS-4/fVDAS, you must use the routines in "fvdas\_convect\_mod.f"
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NFCLDMX( am_I_Root, Input_Opt, State_Met, 
     &                    State_Chm, RC )
!
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DEPO_MERCURY_MOD,   ONLY : ADD_Hg2_WD
      USE DEPO_MERCURY_MOD,   ONLY : ADD_HgP_WD
      USE DEPO_MERCURY_MOD,   ONLY : ADD_Hg2_SNOWPACK
      USE DIAG_MOD,           ONLY : AD37
#if defined( BPCH_DIAG )
      USE DIAG_MOD,           ONLY : AD38
#endif
      USE DIAG_MOD,           ONLY : CONVFLUP
      USE ErrCode_Mod
      USE GC_GRID_MOD,        ONLY : GET_AREA_M2
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PRESSURE_MOD,       ONLY : GET_BP
      USE Species_Mod,        ONLY : Species
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CONV
      USE WETSCAV_MOD,        ONLY : COMPUTE_F
      USE ERROR_MOD,          ONLY : ALLOC_ERR

      IMPLICIT NONE

!
! !INPUT PARAMETERS: 
!
      ! Are we on the root CPU
      LOGICAL,        INTENT(IN)    :: am_I_Root

      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure

! !REMARKS:
!  (1) The "NF" stands for "no flipping", and denotes that you don't have to 
!  flip the species array Q in the main program before passing it to NFCLDMX.
!  (bmy, 2/12/97, 1/31/08)
!
!  (2) This version has been customized to work with GEOS-5 met fields.
!
!  Reference:
!  ============================================================================
!  Lin, SJ.  "Description of the parameterization of cumulus transport
!     in the 3D Goddard Chemistry Transport Model, NASA/GSFC, 1996.
!
!  Vertical indexing:
!  ============================================================================
!  The indexing of the vertical sigma levels has been changed from 
!  SJ-Lin's original code:
!
!                 Old Method          New Method
!                  (SJ Lin)    
!
!               ------------------------------------- Top of Atm.
!                  k = 1               k = NLAY
!               ===================================== Max Extent 
!                  k = 2               k = NLAY-1      of Clouds
!               -------------------------------------
!              
!                   ...                 ...             
!
!               -------------------------------------
!                  k = NLAY-3          k = 4
!               -------------------------------------
!                  k = NLAY-2          k = 3
!               ------------------------------------- Cloud base
!                  k = NLAY-1          k = 2
!               -    -    -    -    -    -    -    - 
!                  k = NLAY            k = 1
!               ===================================== Ground
!
!      which means that:
!
!                 Old Method                      New Method
!                  (SJ Lin)
!
!            k-1      ^                      k+1      ^
!            ---------|---------             ---------|---------
!                     |                               |
!                  CMFMC(k)                       CMFMC(k)
!             
!                                 becomes
!            k     DTRAIN(k),                k     DTRAIN(k),       
!                 QC(k), Q(k)                     QC(k), Q(k)   
!          
!                     ^                               ^
!            ---------|---------             ---------|---------
!                     |                               |   
!            k+1   CMFMC(k+1)                k-1   CMFMC(k-1)
!
!
!      i.e., the lowest level    used to be  NLAY  but is now  1
!            the level below k   used to be  k+1   but is now  k-1.
!            the level above k   used to be  k-1   but is now  k+1
!            the top of the atm. used to be  1     but is now  NLAY.
!
!  The old method required that the vertical dimensions of the CMFMC, DTRAIN, 
!  and Q arrays had to be flipped before and after calling CLDMX.  Also, 
!  diagnostic arrays generated within CLDMX also had to be flipped.  The new 
!  indexing eliminates this requirement (and also saves on array operations).  
!
!  Major Modifications:
!  ============================================================================
!  Original Author:   Shian-Jiann Lin, Code 910.3, NASA/GSFC
!  Original Release:  12 February 1997
!                     Version 3, Detrainment and Entrainment are considered.
!                     The algorithm reduces to that of version 2 if Dtrn = 0.
!                                                                             .
!  Modified By:       Bob Yantosca, for Harvard Atmospheric Sciences
!  Modified Release:  27 January 1998
!                     Version 3.11, contains features of V.3 but also 
!                     scavenges soluble tracer in wet convective updrafts.
!                                                                             .
!                     28 April 1998
!                     Version 3.12, now includes mass flux diagnostic
!                                                                             .
!                     11 November 1999
!                     Added mass-flux diagnostics
!                                                                             .
!                     04 January 2000
!                     Updated scavenging constant AS2
!                                                                             .
!                     14 March 2000
!                     Added new wet scavenging code and diagnostics
!                     based on the GMI algorithm
!                                                                             .
!                     02 May 2000
!                     Added parallel loop over species! 
! 
! !REVISION HISTORY: 
!  12 Feb 1997 - M. Prather  - Initial version
!  (1 ) NFCLDMX is written in Fixed-Form Fortran 90.
!  (2 ) Added TCVV to the argument list.  Also cleaned up argument
!        and local variable declarations. (bey, bmy, 11/10/99)
!  (3 ) AD38 and CONVFLUP are now declared allocatable in "diag_mod.f". 
!        (bmy, 11/29/99)
!  (4 ) Bug fix for tagged CO tracer run (bey, bmy, 1/4/00)
!  (5 ) Add new routines for computing scavenging coefficients,
!        as well as adding the AD37 diagnostic array. (bmy, 3/14/00)
!  (6 ) Updated comments (bmy, 10/2/01)
!  (7 ) Now print a header to stdout on the first call, to confirm that 
!        NFCLDMX has been called (bmy, 4/15/02)
!  (8 ) Remove PZ from the arg list -- it isn't used! (bmy, 8/22/02)
!  (9 ) Fixed ND38 diagnostic so that it now reports correctly (must divide
!        by DNS).  Updatec comments, cosmetic changes. (bmy, 1/27/03)
!  (10) Bug fix: remove duplicate K from PRIVATE declaration (bmy, 3/23/03)
!  (11) Now removed all arguments except NC, TCVV, Q from the arg list -- the
!        other arguments can be supplied via F90 modules.  Now references
!        "dao_mod.f", "grid_mod.f", "pressure_mod.f", and "time_mod.f".
!        (bmy, 3/27/03)
!  (12) Bundled into "convection_mod.f" (bmy, 6/26/03)
!  (13) Make sure K does not go out of bounds in ND38 diagnostic.  Now make 
!        F a 4-D array in order to avoid memory problems on the Altix.  
!        (bmy, 1/27/04)
!  (14) Now references both "ocean_mercury_mod.f" and "tracerid_mod.f".
!        Now call ADD_Hg2_WD from "ocean_mercury_mod.f" to pass the amt of Hg2
!        lost by wet scavenging (sas, bmy, 1/19/05)
!  (15) Now references IS_Hg2 from "tracerid_mod.f".  Now pass tracer # IC
!        to ADD_Hg2_WD. (cdh, bmy, 1/6/06)
!  (16) Bug fix: now only call ADD_Hg2_WD if LDYNOCEAN=T (phs, 2/8/07)
!  (17) Now make CLDMAS, DTRN as arguments, so that we can pass either
!        GEOS-3 or GEOS-3 met data.  Redimension DTCSUM with NC instead of 
!        NNPAR.  In many cases, NC is less than NNPAR and this will help to 
!        save memory especially when running at 2x25 or greater resolution 
!        (bmy, 1/31/08)
!  (18) Add a check to set negative values in Q to TINY (ccc, 4/15/09)
!  (19) Updates for mercury simulation (ccc, 5/17/10)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  29 May 2013 - R. Yantosca - Now set TINY = 1d-60 only for TOMAS code
!  31 May 2013 - R. Yantosca - Now pass State_Chm and then have Q point to 
!                              State_Chm%Tracers.  This is for TOMAS.
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  04 Feb 2014 - R. Yantosca - Bug fix for TOMAS: call COMPUTE_F In its own
!                              separate parallel loop.  Also save the values
!                              of ISOL in the ISOL_SAVE array so that we can
!                              pass them to the parallel tracer loop.
!  18 Apr 2014 - R. Yantosca - Now use a pointer to pass a slice of the
!                              F array to subroutine COMPUTE_F
!  25 Jun 2014 - R. Yantosca - Now pass Input_Opt to COMPUTE_F
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE dinnerences with DELP and
!                              remove dependency on pressure_mod
!  20 Apr 2015 - E. Lundgren - Use DELP*100/g instead of AD/area for BMASS
!                              to keep definition as grid box moist mass/area
!  28 Apr 2015 - E. Lundgren - Remove CLDMAS and DTRM as arguments since 
!                              now included in State_Met
!  09 Jun 2015 - R. Yantosca - Now deposit Hg2, HgP to snowpack regardless of
!                              whether the dynamic ocean is used
!  15 Jun 2015 - E. Lundgren - Now use kg/kg total air as tracer units
!                              (previously v/v)
!  22 Apr 2016 - R. Yantosca - Now get Is_Hg2 and Is_HgP from species database
!  22 Apr 2016 - R. Yantosca - Now pass ThisSpc to ADD_HG2_SNOWPACK
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!  06 Jul 2016 - E. Lundgren - Now use kg/kg dry air as spc units, requiring
!                              use of DELP_DRY instead of DELP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE            :: FIRST = .TRUE.
      INTEGER                  :: I, J, K, KTOP, L, N, NDT, NA
      INTEGER                  :: IC, ISTEP, JUMP, JS, JN, NS
      INTEGER                  :: IMR, JNP, NLAY, AS, Hg_Cat
      REAL(fp)                 :: SDT, CMOUT, ENTRN, DQ, AREA_M2
      REAL(fp)                 :: T0, T1, T2, T3, T4, TSUM, DELQ
      REAL(fp)                 :: DTCSUM(IIPAR,JJPAR,LLPAR,
     &                                   State_Chm%nAdvect )

      REAL(fp), ALLOCATABLE, SAVE :: DSIG(:)

      ! F is the fraction of species lost to wet scavenging in updrafts
      REAL(fp), TARGET         :: F(IIPAR,JJPAR,LLPAR,State_Chm%nAdvect)

      ! Local Work arrays
      REAL(fp)                 :: BMASS( IIPAR, JJPAR, LLPAR )
      REAL(fp)                 :: QB( IIPAR, JJPAR )
      REAL(fp)                 :: MB( IIPAR, JJPAR )
      REAL(fp)                 :: QC( IIPAR, JJPAR ) 

      ! TINY = a very small number
#if defined( TOMAS )
      REAL(fp), PARAMETER      :: TINY  = 1e-60_fp ! Set to 1d-60 for TOMAS code
#else
      REAL(fp), PARAMETER      :: TINY  = 1e-14_fp ! Leave as 1d-14 for std G-C
#endif
      REAL(fp), PARAMETER      :: TINY2 = 1e-30_fp 

      ! ISOL is the wet deposition ID from the species database
      INTEGER                  :: ISOL

      ! QC_PRES and QC_SCAV are the amounts of species 
      ! preserved against and lost to wet scavenging
      REAL(fp)                 :: QC_PRES, QC_SCAV 

      ! DNS is the double precision value for NS
      REAL(fp)                 :: DNS
     
      ! Amt of Hg2, HgP scavenged out of the column (sas, bmy, 1/19/05)
      REAL(fp)                 :: WET_Hg2
      REAL(fp)                 :: WET_HgP 

      ! For fields from Input_Opt
      LOGICAL                  :: ITS_A_MERCURY_SIM
      LOGICAL                  :: LDYNOCEAN
      LOGICAL                  :: LGTMM
      INTEGER                  :: NC

      ! Pointers
      REAL(fp),      POINTER   :: Q        (:,:,:,:)
      REAL(fp),      POINTER   :: CLDMAS   (:,:,:  )
      REAL(fp),      POINTER   :: DTRN     (:,:,:  )
      REAL(fp),      POINTER   :: p_FSOL   (:,:,:  )
      TYPE(Species), POINTER   :: SpcInfo           

      !=================================================================
      ! NFCLDMX begins here!
      !=================================================================

      ! Assume success
      RC                = GC_SUCCESS

      ! Copy values from Input_Opt
      ITS_A_MERCURY_SIM = Input_Opt%ITS_A_MERCURY_SIM
      LDYNOCEAN         = Input_Opt%LDYNOCEAN
      LGTMM             = Input_Opt%LGTMM
      NC                = State_Chm%nAdvect

      ! Point to the chemical species array [kg/kg total air]
      Q                 => State_Chm%Species

      ! Initialize cloud mass flux at upper edge of each level [kg/m2/s]
      ! from Meteorology State object (ewl, 4/28/15)
      CLDMAS            => State_Met%CMFMC(:,:,2:LLPAR+1)

      ! Initialize detrainment mass flux [kg/m2/s] from Meterology State
      ! object (ewl, 4/28/15)
      DTRN              => State_Met%DTRAIN

      ! Initialize pointer to species database
      SpcInfo           => NULL()

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Echo info
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'N F C L D M X  -- by S-J Lin'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

#if defined( GEOS_5 ) || defined( GEOS_FP ) || defined( MERRA )
         ! Do nothing, these are all hybrid grids
#else
         ! NOTE: We don't need to do this for GEOS-5 (bmy, 6/27/07)
         ! DSIG is the sigma-level thickness (NOTE: this assumes that
         ! we are using a pure-sigma grid.  Use new routine for fvDAS.)
         ALLOCATE( DSIG(LLPAR), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'DSIG_CONVECTION' )

         DO L = 1, LLPAR        
            DSIG(L) = GET_BP(L) - GET_BP(L+1)
         ENDDO
#endif

         ! Reset first time flag
         FIRST = .FALSE.
      ENDIF
      
      ! Define dimensions
      IMR  = IIPAR
      JNP  = JJPAR
      NLAY = LLPAR

      ! Convection timestep [s]
      NDT  = GET_TS_CONV() * 60

      !=================================================================
      ! Define active convective region, from J = JS(outh) to 
      ! J = JN(orth), and to level K = KTOP. 
      !
      ! Polar regions are too cold to have moist convection.
      ! (Dry convection should be done elsewhere.)
      !
      ! We initialize the ND14 diagnostic each time we start a new 
      ! time step loop.  Only initialize DTCSUM array if the ND14 
      ! diagnostic is turned on.  This saves quite a bit of time. 
      ! (bmy, 12/15/99)       
      !=================================================================
      IF ( ND14 > 0 ) DTCSUM = 0e+0_fp

      KTOP = NLAY - 1
      JUMP = (JNP-1) / 20
      JS   = 1 + JUMP
      JN   = JNP - JS + 1

      !=================================================================
      ! Internal time step for convective mixing is 300 sec.
      ! Doug Rotman (LLNL) says that 450 sec works just as well.       
      !=================================================================
      NS  = NDT / 300
      NS  = MAX(NS,1)
      SDT = FLOAT(NDT) / FLOAT(NS)
      DNS = DBLE( NS )

!=============================================================================
!  BMASS has units of kg/m^2 and is equivalent to mass (I,J,L) / area (I,J)
!
!   Ps - Pt (mb)| P2 - P1 | 100 Pa |  s^2  | 1  |  1 kg        kg
!  -------------+---------+--------+-------+----+--------  =  -----
!               | Ps - Pt |   mb   | 9.8 m | Pa | m^2 s^2      m^2
!
!  This is done to keep BMASS in the same units as CLDMAS * SDT
!
!  We can parallelize over levels here.  The only quantities that need to 
!  be held local are the loop counters (I, IC, J, JREF, K). (bmy, 5/2/00)
!
!  Now use routine GET_AREA_M2 from "grid_mod.f" to get surface area of
!  grid boxes in m2. (bmy, 2/4/03)
!
!  Now using delta dry pressure * 100 / g for BMASS calculation. 
!  Area-dependent AD / AREA_M2 is no longer used. Note that BMASS is 
!  the dry air kg/m2 as estimated using the dry delta pressure derived 
!  from surface dry pressure, itself derived from GMAO surface wet pressure.
!  (ewl, 7/6/15)
!
!=============================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_M2, K )
!$OMP+SCHEDULE( DYNAMIC )
      DO K = 1, NLAY
      DO J = 1, JJPAR
      DO I = 1, IMR
         BMASS(I,J,K) = State_Met%DELP_DRY(I,J,K) * G0_100
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! (1)  S c a v e n g i n g   i n   C l o u d   U p d r a f t s
      !
      ! Call COMPUTE_F to compute the fraction of species scavenged
      ! in convective cloud updrafts.  COMPUTE_F works for both full
      ! chemistry (NSRCX == 3) and Rn-Pb-Be chemistry (NSRCX == 1)
      ! simulations.  It is best to compute the fraction of species
      ! scavenged at this point, outside the internal time step loop.
      ! This will avoid having to repeat the entire calculation of F
      ! for NS times.
      !
      ! ISOL, which is returned from COMPUTE_F, is the species index
      ! used for diagnostic arrays AD37 and AD38.  ISOL = 0 for all
      ! non-soluble species.
      !
      ! NOTE: We now move this into a separate parallel loop so as to
      ! ensure that TOMAS simulations result in identical results
      ! when parallelization is turned on compared to when it is
      ! turned off. (bmy, sfarina, 2/4/14)
      !=================================================================
      
      ! Loop over the advected species
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED   )
!$OMP+PRIVATE( NA, IC, ISOL, I, J, K, p_FSOL ) 
!$OMP+SCHEDULE( DYNAMIC )
      DO NA = 1, NC
         
         ! Get the species ID from the advected species ID
         IC = State_Chm%Map_Advect(NA)

         ! Point to 3D slice of F
         p_FSOL => F(:,:,:,NA)

         ! Compute the fraction of species IC scavenged in updrafts
         ! NOTE: ISOL is the wetdep Id from the species DB (bmy, 7/7/16)
         CALL COMPUTE_F( am_I_Root, IC, p_FSOL, ISOL, 
     &                   Input_Opt, State_Met, State_Chm, RC ) 

         ! Free pointer memory
         p_FSOL => NULL()

         ! ND37 diagnostic -- store F only for soluble species
         IF ( ND37 > 0 .and. ISOL > 0 ) THEN
            DO K = 1, LD37
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               AD37(I,J,K,ISOL) = AD37(I,J,K,ISOL) + F(I,J,K,NA) 
            ENDDO
            ENDDO
            ENDDO
         ENDIF
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! (2)  A d v e c t e d   S p e c i e s   L o o p 
      !
      ! We now parallelize over species, since species are independent 
      ! of each other.  The parallel loop only takes effect if you 
      ! compile with the f90 "-mp" switch.  Otherwise the compiler will 
      ! interpret the parallel-processing directives as comments, and 
      ! the loop will execute on a single thread.
      !
      ! The following types of quantities must be held local for 
      ! parallelization:
      ! (1) Loop counters ( I, IC, ISTEP, J, K )
      ! (2) Scalars that are assigned values inside the species loop: 
      !     ( CMOUT, DELQ, ENTRN, ISOL, QC_PRES, etc. )
      ! (3) Arrays independent of species ( F, MB, QB, QC )
      !=================================================================

      ! Loop over only the advected species
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( CMOUT, DELQ, ENTRN, I, IC, ISTEP, ISOL, J, AREA_M2, K  ) 
!$OMP+PRIVATE( MB, QB, QC, QC_PRES, QC_SCAV, T0, T1, T2, T3, T4, TSUM )
!$OMP+PRIVATE( WET_Hg2, WET_HgP, SpcInfo, Hg_Cat, NA                  )
!$OMP+SCHEDULE( DYNAMIC )
      DO NA = 1, NC
         
         ! Get the species ID from the advected species ID
         IC      =  State_Chm%Map_Advect(NA)

         ! Point to the Species Database entry for this species
         SpcInfo => State_Chm%SpcData(IC)%Info

         ! Also get wetdep ID that corresponds to IC
         ! If not a soluble species, then ISOL = -999
         ISOL    =  SpcInfo%WetDepID

         !==============================================================
         ! (3)  I n t e r n a l   T i m e   S t e p   L o o p
         !
         ! The internal time step is currently 300 seconds.  
         !==============================================================
         DO ISTEP = 1, NS

            !===========================================================
            ! (4)  B e l o w   C l o u d   B a s e  (K < 3)
            !
            ! Loop over longitude and latitude (I,J), and consider what
            ! is below the cloud base (i.e. below the 2nd sigma level).  
            !===========================================================
            DO J = JS, JN
            DO I = 1, IMR

               !========================================================
               ! (4.1) If Cloud Mass Flux exists at (I,J,2), 
               !       then compute QB.
               !
               ! QB is "weighted average" mixing ratio below the cloud 
               ! base.  QB is used to compute QC, which is the mixing 
               ! ratio of the air that moved in cumulus transport up 
               ! to the next level.  MB is the dry mass of air below 
               ! the cloud base.
               !========================================================

               IF ( CLDMAS(I,J,2) .gt. TINY ) THEN 
#if   defined( GEOS_5 )

                  ! Need to replace DSIG w/ the difference 
                  ! of pressure edges for GEOS-5 (bmy, 6/27/07)
                  ! Now use DELP_DRY, equivalent to delta dry pressure 
                  ! across the level [hPa] (ewl, 7/6/15)
                  QB(I,J) = 
     &                      ( Q(I,J,1,IC) * State_Met%DELP_DRY(I,J,1) + 
     &                      Q(I,J,2,IC) * State_Met%DELP_DRY(I,J,2) ) / 
     &                      ( State_Met%DELP_DRY(I,J,1) 
     &                      + State_Met%DELP_DRY(I,J,2) )

#else
                  ! for GEOS-3
                  QB(I,J) = 
     &               ( Q(I,J,1,IC) * DSIG(1)   + 
     &                 Q(I,J,2,IC) * DSIG(2) ) / ( DSIG(1) + DSIG(2) )

#endif

                  ! Dry mass of air below the cloud base
                  MB(I,J) = BMASS(I,J,1) + BMASS(I,J,2)

!=============================================================================
!         Dry mass of species below cloud base (i.e. MB  * QB ) +   
!         Subsidence into cloud base from above (i.e. SDT * C(2) * Q(3) ) 
!  QC =  -----------------------------------------------------------------
!             Dry air mass below cloud base (i.e. MB + C(2)*Q(3) )
!=============================================================================
                  QC(I,J) = 
     &               ( MB(I,J)       * QB(I,J)       + 
     &                 CLDMAS(I,J,2) * Q(I,J,3,IC)   * SDT ) /
     &               ( MB(I,J)       + CLDMAS(I,J,2) * SDT )
     
                  !=====================================================
                  ! DQ = QB - QC is the dry mass to be transported 
                  ! out of the cloud base.  Changes below cloud base 
                  ! are proportional to the background mass.
                  ! 
                  ! Subtract DQ from Q(*,*,K=1,*) and from Q(*,*,K=2,*), 
                  ! but do not make Q(*,*,K=1,*) or Q(*,*,K=2,*) < 0.
                  !=====================================================
                  Q(I,J,2,IC) = QC(I,J)
                  Q(I,J,1,IC) = QC(I,J)

               !========================================================
               ! If there is no Cloud mass flux, set QC = Q(K=3) 
               ! at this I,J location
               !========================================================
               ELSE
                  QC(I,J) = Q(I,J,3,IC)
               ENDIF
            ENDDO !I
            ENDDO !J

            !===========================================================
            ! (5)  A b o v e   C l o u d   B a s e
            !
            ! Loop over sigma levels K and longitude and latitude (I,J) 
            !
            ! Recall that K+1 is the level BELOW level K, and that K-1 
            ! is the level ABOVE level K.  This is the way that SJ Lin 
            ! indexes the sigma levels.
            !===========================================================
            DO K = 3,  KTOP
            DO J = JS, JN
            DO I = 1,  IMR

               ! Grid box surface area [m2]
               AREA_M2 = GET_AREA_M2( I, J, K )

               ! Initialize
               WET_Hg2 = 0.0_fp
               WET_HgP = 0.0_fp
               
!==============================================================================
!  (5.1)  M a s s   B a l a n c e   i n   C l o u d  ===>  QC(I,J)
!
!  QC_PRES = QC(I,J,K-1) * AP = amt of Qc preserved against wet scavenging
!  QC_SCAV = QC(I,J,K-1) * AS = amt of Qc lost to wet scavenging
!  CMOUT   = air mass flowing out of cloud at level K 
!  ENTRN   = Entrainment: air mass flowing into cloud at level K   
!
!  If ENTRN > 0 then compute the new value of QC(I,J):
!
!                CLDMAS(K-1) * QC_PRES  +  ENTRN(K) * Q(K)
!    QC(I,J) =  -------------------------------------------
!                      CLDMAS(I,J,K) + DTRN(I,J,K)
!
!            =   species mass coming in from below      (i.e. level K-1) + 
!                species mass coming in from this level (i.e. level K)
!               -----------------------------------------------------------
!                             dry mass coming into cloud
!
!  Otherwise, preserve the previous value of QC(I,J).  This will ensure
!  that TERM1 - TERM2 is not a negative quantity (see below).
!  
!  Entrainment must be >= 0 (since we cannot have a negative flux of air
!  into the cloud).  This condition is strong enough to ensure that
!  CMOUT > 0 and will prevent floating-point exception.
!==============================================================================

               IF ( CLDMAS(I,J,K-1) .gt. TINY ) THEN
                  CMOUT    = CLDMAS(I,J,K) + DTRN(I,J,K)
                  ENTRN    = CMOUT         - CLDMAS(I,J,K-1)
!-----------------------------------------------------------------------------
! Prior to 7/7/16:
! F loops over the # of advected species, so use NA index (bmy, 7/7/16)
!                  QC_PRES  = QC(I,J) * ( 1e+0_fp - F(I,J,K,IC) )
!-----------------------------------------------------------------------------
                  QC_PRES  = QC(I,J) * ( 1e+0_fp - F(I,J,K,NA) )
!-----------------------------------------------------------------------------
! Prior to 7/7/16:
! F loops over the # of advected species, so use NA index (bmy, 7/7/16)
!                  QC_SCAV  = QC(I,J) * F(I,J,K,IC)
!-----------------------------------------------------------------------------
                  QC_SCAV  = QC(I,J) * F(I,J,K,NA)

!#if defined ( DEBUG )
!       IF ( IC == 129 ) THEN
!          print *, I, J, K, QC_PRES, QC_SCAV, QC(I,J), F(I,J,K,NA)
!       ENDIF
!#endif
                  
                  IF ( ENTRN .ge. 0 ) THEN
                     QC(I,J) = ( CLDMAS(I,J,K-1) * QC_PRES       + 
     &                           ENTRN           * Q(I,J,K,IC) ) / 
     &                           CMOUT
                  ENDIF   

!==============================================================================
!  (5.2)  M a s s   B a l a n c e   i n   L e v e l  ===>  Q(I,J,K,IC)
!
!  The cumulus transport above the cloud base is done as follows:
!     C_k-1  = cloud air mass flux from level k-1 to level k
!     C_k    = cloud air mass flux from level k   to level k+1
!     QC_k-1 = mixing ratio of species INSIDE CLOUD at level k-1
!     QC_k   = mixing ratio of species INSIDE CLOUD at level k
!     Q_k    = mixing ratio of species in level k
!     Q_k+1  = mixing ratio of species in level k+1
!                                
!                       |                    |
!      k+1     ^        |       Cloud        |3)      C_k * Q_k+1
!              |        |         ^          |            |
!      --------|--------+---------|----------+------------|--------
!              |        |         |          |            V
!      k     C_k        |2)   C_k * QC_k     |
!                       |                    |
!                       |                    |
!                       |         ^          |4)    C_k-1 * Q_k           
!              ^        |         |          |            |
!      --------|--------+---------|----------+------------|----------
!              |        |         |          |            |
!      k-1   C_k-1      |1) C_k-1 * QC_k-1   |            V
!                       |         * AP       |
!
!  There are 4 terms that contribute to mass flow in and out of level k:
!
!  1) C_k-1 * QC_PRES = species convected from level k-1 to level k 
!  2) C_k   * QC_k    = species convected from level k   to level k+1 
!  3) C_k   * Q_k+1   = species subsiding from level k+1 to level k 
!  4) C_k-1 * Q_k     = species subsiding from level k   to level k-1 
!
!  Therefore the change in species concentration is given by
!     DELQ = (Term 1) - (Term 2) + (Term 3) - (Term 4)
!
!  and Q(I,J,K,IC) = Q(I,J,K,IC) + DELQ.  
!
!  The term T0 is the amount of species that is scavenged out of the box.
!  Compute that term here for the ND38 diagnostic below. (bmy, 1/27/03)
!==============================================================================
                  T0   =  CLDMAS(I,J,K-1) * QC_SCAV   
                  T1   =  CLDMAS(I,J,K-1) * QC_PRES
                  T2   = -CLDMAS(I,J,K  ) * QC(I,J       )
                  T3   =  CLDMAS(I,J,K  ) * Q (I,J,K+1,IC)
                  T4   = -CLDMAS(I,J,K-1) * Q (I,J,K,  IC)

                  TSUM = T1 + T2 + T3 + T4 

                  DELQ = ( SDT / BMASS(I,J,K) ) * TSUM

                  ! If DELQ > Q then do not make Q negative!!!
                  IF ( Q(I,J,K,IC) + DELQ < 0 ) THEN
                     DELQ = -Q(I,J,K,IC)
                  ENDIF

                  Q(I,J,K,IC) = Q(I,J,K,IC) + DELQ
                                
                  !=====================================================
                  ! ND14 Diagnostic: Upward mass flux due to wet 
                  ! convection.  The diagnostic ND14 works only for 
                  ! levels >= than 3.  This is ok for now since I 
                  ! want to have a closed budget with a box starting 
                  ! from the ground.
                  !
                  ! DTCSUM(I,J,K,IC) is the flux (kg/box/sec) in 
                  ! the box (I,J), for the species IC going out of 
                  ! the top of the layer K to the layer above (K+1) 
                  ! (bey, 11/10/99).      
                  !=====================================================
                  IF ( ND14 > 0 ) THEN
                     DTCSUM(I,J,K,IC) = DTCSUM(I,J,K,IC) + (-T2-T3)  
     &                                * AREA_M2 
                  ENDIF

#if defined( BPCH_DIAG )
                  !=====================================================
                  ! ND38 Diagnostic: loss of soluble species to wet
                  ! scavenging in cloud updrafts [kg/s].  We must 
                  ! divide by DNS, the # of internal timesteps.
                  !=====================================================
                  IF ( ( ND38 > 0 .or. LGTMM )
     &                 .and. ISOL > 0 .and. K <= LD38 ) THEN

                     AD38(I,J,K,ISOL) = AD38(I,J,K,ISOL) 
     &                                + T0 * AREA_M2 / DNS  

                  ENDIF
#endif
 
                  !=====================================================
                  ! Pass the amount of Hg2 and HgP lost in wet 
                  ! scavenging [kg] to "ocean_mercury_mod.f" via 
                  ! ADD_Hg2_WET and ADD_HgP_WET. We must also divide 
                  ! by DNS, the # of internal timesteps. 
                  ! (sas, bmy, eck, eds, 1/19/05, 1/6/06, 7/30/08)
                  !
                  ! IMPORTANT NOTE: WET_Hg2 and WET_HgP are defined
                  ! on each 3-D (I,J,K) iteration.  But in routines
                  ! ADD_Hg2_WD and ADD_Hg2_SNOWPACK, these are saved to 
                  ! diagnostic arrays that have only 2-D (I,J) spatial 
                  ! size.  This can cause slight numerical differences 
                  ! when OpenMP parallelization is turned on, because 
                  ! one CPU is trying to write to the diagnostic arrays
                  ! simultaneously.  To avoid this situation, we call 
                  ! ADD_Hg2_WD and ADD_Hg2_SNOWPACK from within an 
                  ! !$OMP CRITICAL block, which ensures that only one 
                  ! CPU at a time can write to the diagnostic arrays.
                  ! (bmy, 4/20/16)
                  !=====================================================
                  IF ( ITS_A_MERCURY_SIM ) THEN
!$OMP CRITICAL
                     ! Is it a Hg2 species?
                     IF ( SpcInfo%Is_Hg2 ) THEN

                        ! Wet scavenged Hg(II) in [kg/s], converted
                        ! to [kg] by multiplying by NDT
                        WET_Hg2 = T0 * AREA_M2 / DNS * NDT
                        
                        ! Category # for this Hg2 species
                        Hg_Cat  = SpcInfo%Hg_Cat

                        ! Pass to "ocean_mercury_mod.f"
                        CALL ADD_Hg2_WD      ( I, J, Hg_Cat, WET_Hg2   )
                        CALL ADD_Hg2_SNOWPACK( I, J, Hg_Cat, WET_Hg2,
     &                                                       State_Met )
                     ENDIF

                     ! Is it a HgP species?
                     IF ( SpcInfo%Is_HgP ) THEN 

                        ! Wet scavenged Hg(P) in [kg/s], converted 
                        ! to [kg] by multiplying by NDT
                        WET_HgP = T0 * AREA_M2 / DNS * NDT

                        ! Category # for this Hg2 species
                        Hg_Cat  = SpcInfo%Hg_Cat

                        ! Pass to "ocean_mercury_mod.f"
                        CALL ADD_HgP_WD      ( I, J, Hg_Cat, WET_HgP   ) 
!------------------------------------------------------------------------------
! NOTE: This call to ADD_HG2_SNOWPACK is known to cause very small (numerical
! noise) differences in slot #4 of the PL-HG2-$ diagnostic when comparing
! a GEOS-5 Hg simulation done on a single processor with one done on multi-
! processors.  (In other words, the unit test result is "yellow").  The 
! magnitude of the differences are very small, and on the order of numerical 
! noise.
!
! We hypothesize that this issue may have two causes: (1) updating the AD03 
! array (which is 2-D) from within a 3-D parallel loop, and (2) possible 
! numerical roundoff when casting the REAL*8 value to REAL*4.  
!
! The results of the simulation in the species restart files are identical,
! so technically, the GEOS-5 Hg simulations pass the unit test. 
!
! It is probably not possible to avoid the numerical differences in slot #4
! of the PL-HG2-$ diagnostic without redesiging the diagnostics.  For the 
! time being, we will leave the code as-is but just note that differences
! will occur in this diagnostic.
! 
!    -- Bob Yantosca (22 Apr 2016)
! 
                        CALL ADD_Hg2_SNOWPACK( I, J, Hg_Cat, WET_HgP,
     &                                                       State_Met )
!------------------------------------------------------------------------------
                     ENDIF
!$OMP END CRITICAL
                  ENDIF

               !========================================================
               ! No cloud transport if cloud mass flux < TINY; 
               ! Change Qc to q
               !========================================================
               ELSE                     
                  QC(I,J) = Q(I,J,K,IC)

#if   defined( GEOS_5 ) 
                     !--------------------------------------------------
                     ! FIX FOR GEOS-5 MET FIELDS!
                     !
                     ! Bug fix for the cloud base layer, which is not 
                     ! necessarily in the boundary layer, and for 
                     ! GEOS-5, there could be "secondary convection 
                     ! plumes - one in the PBL and another one not.
                     !
                     ! NOTE: T2 and T3 are the same terms as described
                     ! in the above section.
                     !
                     ! (swu, 08/13/2007)
                     !--------------------------------------------------
                     IF ( CLDMAS(I,J,K) > TINY ) THEN 

                        ! Species convected from K -> K+1 
                        T2   = -CLDMAS(I,J,K  ) * QC(I,J)

                        ! Species subsiding from K+1 -> K 
                        T3   =  CLDMAS(I,J,K  ) * Q (I,J,K+1,IC)

                        ! Change in species concentration
                        DELQ = ( SDT / BMASS(I,J,K) ) * (T2 + T3)

                        ! If DELQ > Q then do not make Q negative!!!
                        IF ( Q(I,J,K,IC) + DELQ < 0.0e+0_fp ) THEN 
                            DELQ = -Q(I,J,K,IC)
                        ENDIF
  
                        ! Add change in species to Q array
                        Q(I,J,K,IC) = Q(I,J,K,IC) + DELQ

                     ENDIF
#endif
                  ENDIF
               ENDDO        !I
            ENDDO           !J
            ENDDO           !K
         ENDDO              !NSTEP

         ! Free pointer
         SpcInfo => NULL()
      ENDDO                 !NA (also IC)
!$OMP END PARALLEL DO

      !=================================================================
      ! ND14 Diagnostic: Store into the CONVFLUP array.  
      ! Also divide by the number of internal timesteps (DNS).
      !================================================================= 
      IF ( ND14 > 0 ) THEN

         ! Loop over only the advected species
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, NA, IC, J, K ) 
!$OMP+SCHEDULE( DYNAMIC )
         DO NA = 1,  NC

            ! Species ID
            IC = State_Chm%Map_Advect(NA)

            ! Archive diagnostic
            DO K  = 3,  KTOP
            DO J  = JS, JN
            DO I  = 1,  IMR
               CONVFLUP(I,J,K,IC) = CONVFLUP(I,J,K,IC) +
     &                              ( DTCSUM(I,J,K,IC) / DNS )
            ENDDO
            ENDDO
            ENDDO
         ENDDO 
!$OMP END PARALLEL DO
      ENDIF

      ! Add a check to set negative values in Q to TINY2
      ! Loop over only the advected species
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED   )
!$OMP+PRIVATE( NA, IC, L, J, I )
      DO NA = 1,NC
         
         ! Species ID
         IC = State_Chm%Map_Advect(NA)

         DO L = 1,LLPAR
         DO J = 1,JJPAR
         DO I = 1,IIPAR
            Q(I,J,L,IC) = MAX(Q(I,J,L,IC),TINY2)
         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer memory
      NULLIFY( Q, CLDMAS, DTRN )

      END SUBROUTINE NFCLDMX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_merra_convection
!
! !DESCRIPTION: Subroutine DO\_MERRA\_CONVECTION (formerly called NFCLDMX)
!  is S-J Lin's cumulus transport module for 3D GSFC-CTM, modified for the 
!  GEOS-Chem model.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_MERRA_CONVECTION( am_I_Root,  
     &                                Input_Opt,
     &                                State_Met,  
     &                                State_Chm,
     &                                I,          
     &                                J,
     &                                AREA_M2,    
     &                                F,          
     &                                TS_DYN,   
     &                                USE_DIAG14, 
     &                                DIAG14,
     &                                USE_DIAG38, 
     &                                DIAG38,
     &                                RC          )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE DEPO_MERCURY_MOD,   ONLY : ADD_Hg2_SNOWPACK
      USE DEPO_MERCURY_MOD,   ONLY : ADD_Hg2_WD
      USE DEPO_MERCURY_MOD,   ONLY : ADD_HgP_WD
      USE ERROR_MOD,          ONLY : IT_IS_NAN
      USE ERROR_MOD,          ONLY : IT_IS_FINITE
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP  ! hma Nov 3, debug
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE Species_Mod,        ONLY : Species
      USE WETSCAV_MOD,        ONLY : H2O2s_3D => H2O2s ! [v/v]
      USE WETSCAV_MOD,        ONLY : SO2s_3D  => SO2s  ! [v/v]
      USE WETSCAV_MOD,        ONLY : WASHOUT
      USE WETSCAV_MOD,        ONLY : LS_K_RAIN
      USE WETSCAV_MOD,        ONLY : LS_F_PRIME
!
! !INPUT PARAMETERS:
!  
      LOGICAL,        INTENT(IN)    :: am_I_Root  ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met  ! Meteorology State object
      INTEGER,        INTENT(IN)    :: I, J       ! Lon & lat indices
      REAL(fp),       INTENT(IN)    :: AREA_M2    ! Surface area [m2]
      REAL(fp),       INTENT(IN)    :: F(:,:)     ! Fraction of soluble species
                                                  !  for updraft scavenging 
                                                  !  [unitless].  Computed by
                                                  !  routine  COMPUTE_F.
      REAL(fp),       INTENT(IN)    :: TS_DYN     ! Dynamic timestep [min]
      LOGICAL,        INTENT(IN)    :: USE_DIAG14 ! Archive DIAG14?
      LOGICAL,        INTENT(IN)    :: USE_DIAG38 ! Archive DIAG38?
!                                                  
! !INPUT/OUTPUT PARAMETERS:                        
!        
      TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
!                                                  
! !OUTPUT PARAMETERS:                              
!                     
      REAL(fp),       INTENT(OUT)   :: DIAG14(:,:) ! Array for ND14 diagnostic
      REAL(fp),       INTENT(OUT)   :: DIAG38(:,:) ! Array for ND38 diagnostic
      INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!  Reference:
!  ============================================================================
!  Lin, SJ.  "Description of the parameterization of cumulus transport
!     in the 3D Goddard Chemistry Transport Model, NASA/GSFC, 1996.
!                                                                             .
!  Unit conversion for BMASS:
!
!      Ps - Pt (mb)| P2 - P1 | 100 Pa |  s^2  | 1  |  1 kg        kg
!     -------------+---------+--------+-------+----+--------  =  -----
!                  | Ps - Pt |   mb   | 9.8 m | Pa | m^2 s^2      m^2
!
!                                                                             .
!  NOTE: We are passing I & J down to this routine so that it can call the
!  proper code from "mercury_mod.f".  Normally, we wouldn't pass I & J as
!  arguments to columnized code.  This prevents rewriting the mercury_mod.f
!  routines ADD_Hg2_
!
! !REVISION HISTORY:
!  15 Jul 2009 - R. Yantosca - Columnized and cleaned up.
!                            - CLDMAS renamed to CMFMC and DTRN renamed
!                              to DTRAIN for consistency w/ GEOS-5.
!  17 Jul 2009 - R. Yantosca - Now do unit conversion of Q array from
!                              [kg] --> [v/v] and vice versa internally
!  14 Dec 2009 - R. Yantosca - Now remove internal unit conversion, since
!                              Q now comes in as [mol/mol] (=[v/v]) from the
!                              calling routine.
!  14 Dec 2009 - R. Yantosca - Remove COEF from the argument list
!  06 May 2010 - R. Yantosca - Now add IDENT via the argument list
!  29 Sep 2010 - R. Yantosca - Modified for MERRA met fields
!  05 Oct 2010 - R. Yantosca - Now pass COEF via the argument list
!  05 Oct 2010 - R. Yantosca - Attach ND14 and ND38 diagnostics
!  15 Oct 2010 - H. Amos     - Added BXHEIGHT and T as arguments
!  15 Oct 2010 - R. Yantosca - Added I, J, H2O2s and SO2s as arguments
!  15 Oct 2010 - H. Amos     - Added scavenging below cloud base
!  06 Apr 2011 - M.Fu, H.Amos- Bug fix: make sure washout adheres to the same
!                              algorithm as in the wet deposition code.
!  27 Jul 2011 - R. Yantosca - Declare CLDBASE as INTEGER to avoid PGI errors
!  16 Aug 2011 - J. Fisher   - Bug fix: use IS_Hg2() and IS_HgP to test if 
!                              a tracer is Hg2 or HgP (for tagged species)
!  16 Aug 2011 - J. Fisher   - Now use WETLOSS instead of T0_SUM in the ND38 
!                              diagnostic below the cloud.  Using T0_SUM leads 
!                              us to over-count the tracer scavenged out of 
!                              the column.
!  22 Oct 2012 - R. Yantosca - Now reference Headers/gigc_errcode_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  31 May 2013 - R. Yantosca - Now pass State_Chm to WASHOUT
!  05 Sep 2013 - R. Yantosca - Bug fix: DT is apparently undefined, but still
!                              passed to WASHOUT.  Use SDT instead.  This
!                              avoids a floating-point error.
!  18 Apr 2014 - R. Yantosca - Now point to 3-D arrays internally
!  18 Apr 2014 - R. Yantosca - Now also pass N_TRACERS (to facilitate APM)
!  18 Apr 2014 - R. Yantosca - Remove code that we don't need anymore
!  04 Feb 2015 - M. Sulprizio- Fix calculation of WETLOSS for non-aerosol 
!                              tracers below the cloud base (C. Friedman)
!  20 Apr 2015 - E. Lundgren - Use DELP*100/g instead of AD/area as grid box
!                              moist mass per area and remove AD from routine
!  20 May 2015 - M. Sulprizio- Apply bug fixes provided by Viral Shah:
!                              -- Remove F(K,IC) > 0 condition that prevents
!                                 archiving of deposited mass in DIAG38
!                              -- Add statement that subtracts the wet deposited
!                                 amount from the atmospheric mass
!                              -- Fix inconsistency in units when T0_SUM is used
!  04 Jun 2015 - E. Lundgren - Adapt Viral Shah bug fixes to moist mixing ratio
!  09 Jun 2015 - R. Yantosca - Now deposit Hg2, HgP to snowpack regardless of
!                              whether the dynamic ocean is used
!  15 Jun 2015 - E. Lundgren - Now use kg/kg total air as tracer units not v/v
!  22 Jun 2015 - E. Lundgren - Move QB_NUM calculation to within timestep loop
!  12 Aug 2015 - R. Yantosca - Treat MERRA2 in same way as we do for GEOS-FP
!  14 Sep 2015 - E. Lundgren - Apply bug fixes provided by Viral Shah:
!                              -- Prevent ALPHA > 1 in washout of aerosols 
!                              -- Add tracer GAINED to Q before WETLOSS
!                                 calculation in aerosol washout
!  22 Apr 2016 - R. Yantosca - Now get Is_Hg2 & Is_HgP from species database
!  25 Apr 2016 - R. Yantosca - Now pass Hg category # to ADD_Hg2_* functions
!  28 Apr 2016 - R. Yantosca - Rewrite Is_Hg block to avoid unassociated
!                              pointer seg faults
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!  05 Jul 2016 - R. Yantosca - Now replace N_TRACERS argument with the
!                              State_Chm%nAdvect field
!  06 Jul 2016 - E. Lundgren - Now use kg/kg dry air as spc units, requiring
!                              use of DELP_DRY instead of DELP and PEDGE_DRY
!                              for avg mixing ratio 
!  07 Jul 2016 - R. Yantosca - Bug fix: F and ISOL need to be indexed with
!                              the advected species index NA
!  07 Jul 2016 - R. Yantosca - DIAG14 and DIAG38 now are dimensioned with
!                              of size State_Chm%nWetDep instead of N_TRACERS
!  20 Sep 2016 - R. Yantosca - Rewrote IF test to avoid Gfortran compiler error
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL(fp), PARAMETER    :: TINYNUM = 1e-14_fp
!                            
! !LOCAL VARIABLES:          
!                            
      ! Scalars              
      LOGICAL                :: AER,         IS_Hg 
      INTEGER                :: IC,          ISTEP,     K        
      INTEGER                :: KTOP,        NC,        NDT
      INTEGER                :: NLAY,        NS,        CLDBASE
      INTEGER                :: Hg_Cat,      NA,        nAdvect, NW
      REAL(fp)               :: CMFMC_BELOW, ALPHA,     ALPHA2
      REAL(fp)               :: CMOUT,       DELQ,      DQ
      REAL(fp)               :: DNS,         ENTRN,     QC
      REAL(fp)               :: QC_PRES,     QC_SCAV,   SDT
      REAL(fp)               :: T0,          T0_SUM,    T1    
      REAL(fp)               :: T2,          T3,        T4 
      REAL(fp)               :: TSUM,        LOST,      GAINED  
      REAL(fp)               :: WETLOSS,     MASS_WASH, MASS_NOWASH
      REAL(fp)               :: QDOWN,       DT,        F_WASHOUT
      REAL(fp)               :: K_RAIN,      WASHFRAC,  WET_Hg2
      REAL(fp)               :: WET_HgP,     MB,        QB
      REAL(fp)               :: QB_NUM,      DELP_DRY_NUM

                             
      ! Arrays               
      REAL(fp)               :: BMASS    (LLPAR    )    
      REAL(fp)               :: PDOWN    (LLPAR    )

      ! Pointers
      REAL(fp),      POINTER :: BXHEIGHT (:        )
      REAL(fp),      POINTER :: CMFMC    (:        )
      REAL(fp),      POINTER :: DQRCU    (:        )
      REAL(fp),      POINTER :: DTRAIN   (:        )
      REAL(fp),      POINTER :: PFICU    (:        )
      REAL(fp),      POINTER :: PFLCU    (:        )
      REAL(fp),      POINTER :: REEVAPCN (:        )
      REAL(fp),      POINTER :: DELP_DRY (:        )
      REAL(fp),      POINTER :: T        (:        )
      REAL(fp),      POINTER :: H2O2s    (:        )
      REAL(fp),      POINTER :: SO2s     (:        )
      REAL(fp),      POINTER :: Q        (:,:      )
      TYPE(Species), POINTER :: SpcInfo

      !========================================================================
      ! (0)  I n i t i a l i z a t i o n
      !========================================================================
      
      ! Point to columns of derived-type object fields
      BXHEIGHT => State_Met%BXHEIGHT(I,J,:        ) ! Box height [m]
      CMFMC    => State_Met%CMFMC   (I,J,2:LLPAR+1) ! Cloud mass flux [kg/m2/s] 
      DQRCU    => State_Met%DQRCU   (I,J,:        ) ! Precip production rate:
      DTRAIN   => State_Met%DTRAIN  (I,J,:        ) ! Detrainment flux [kg/m2/s]
      REEVAPCN => State_Met%REEVAPCN(I,J,:        ) ! Evap of precip'ing conv.
      DELP_DRY => State_Met%DELP_DRY(I,J,:        ) ! Edge dry P diff [hPa]
      T        => State_Met%T       (I,J,:        ) ! Air temperature [K]
      H2O2s    => H2O2s_3D          (I,J,:        ) ! H2O2s from sulfate_mod
      SO2s     => SO2s_3D           (I,J,:        ) ! SO2s from sulfate_mod
      Q        => State_Chm%Species (I,J,:,:      ) ! Chemical species 
                                                    ! [mol/mol dry air]
      SpcInfo  => NULL()                            ! Species database entry
                                                
#if   defined( GEOS_FP ) || defined( MERRA2 )

      !%%%%% FOR GEOS-FP/MERRA2, PFICU and PFLCU are on level edges %%%%%
      PFICU    => State_Met%PFICU   (I,J,2:LLPAR+1) ! Dwnwd flx of conv 
                                                    !  ice precip [kg/m2/s]
      PFLCU    => State_Met%PFLCU   (I,J,2:LLPAR+1) ! Dwnwd flux of conv   
                                                    !  liquid precip [kg/m2/s]
#elif defined( MERRA   )                              

      !%%%%% FOR MERRA, PFICU and PFLCU are on level centers %%%%%
      PFICU    => State_Met%PFICU   (I,J,:        )
      PFLCU    => State_Met%PFLCU   (I,J,:        )

#endif
        
      ! # of levels and # of species
      NLAY     = LLPAR
      NC       = State_Chm%nAdvect

      ! Top level for convection
      KTOP     = NLAY - 1

      ! Convection timestep [s]
      NDT      = TS_DYN * 60e+0_fp

      ! Internal time step for convective mixing is 300 sec.
      ! Doug Rotman (LLNL) says that 450 sec works just as well.       
      NS       = NDT / 300                ! Num internal timesteps (int)
      NS       = MAX( NS, 1 )             ! Set lower bound to 1
      DNS      = DBLE( NS )               ! Num internal timesteps (real)
      SDT      = DBLE( NDT ) / DBLE( NS ) ! seconds in internal timestep
      
      !-----------------------------------------------------------------
      ! Determine location of the cloud base, which is the level where
      ! we start to have non-zero convective precipitation formation
      !-----------------------------------------------------------------

      ! Minimum value of cloud base is the surface level
      CLDBASE = 1

      ! Find the cloud base
      DO K = 1, NLAY
         IF ( DQRCU(K) > 0e+0_fp ) THEN
            CLDBASE = K
            EXIT
         ENDIF
      ENDDO

      !-----------------------------------------------------------------
      ! Compute PDOWN and BMASS
      !-----------------------------------------------------------------
      DO K = 1, NLAY

         ! PDOWN is the convective precipitation leaving each
         ! box [cm3 H2O/cm2 air/s]. This accounts for the 
         ! contribution from both liquid & ice precip.
         ! PFLCU and PFICU are converted from kg/m2/s to m3/m2/s
         ! using water and ice densities, respectively.
         ! m3/m2/s becomes cm3/cm2/s using a factor of 100.
         PDOWN(K) = ( ( PFLCU(K) / 1000e+0_fp ) 
     &            +   ( PFICU(K) /  917e+0_fp ) ) * 100e+0_fp

         ! BMASS is the dry air mass per unit area for the grid box
         ! bounded by level K and K+1 [kg/m2] 
         ! BMASS is equivalent to deltaP (dry) * 100 / g
         ! This is done to keep BMASS in the same units as CMFMC * SDT
         BMASS(K) = DELP_DRY(K) * G0_100

      ENDDO

      !-----------------------------------------------------------------
      ! Compute MB, the mass per unit area of dry air below the cloud
      ! base [kg/m2]. Calculate MB by looping over levels below the 
      ! cloud base.
      !-----------------------------------------------------------------
      MB = 0e+0_fp 
      DO K = 1, CLDBASE-1
         MB = MB + BMASS(K)
      ENDDO

      ! Is this a Hg simulation?
      IS_Hg = Input_Opt%ITS_A_MERCURY_SIM

      !========================================================================
      ! (1)  A d v e c t e d   S p e c i e s   L o o p 
      !========================================================================

      ! Loop over only the advected species
      DO NA = 1, NC

         ! Get the species ID (modelID) from the advected species ID
         IC       =  State_Chm%Map_Advect(NA)

         ! Look up the corresponding entry in the species database
         SpcInfo  => State_Chm%SpcData(IC)%Info 

         ! Also get the corresponding wetdep ID
         NW       =  SpcInfo%WetDepId

         ! Zero the ND14 and ND38 diagnostic arrays 
         ! These are only sized for the # of soluble species
         IF ( NW > 0 ) THEN
            DIAG14(:,NW) = 0.0_fp 
            DIAG38(:,NW) = 0.0_fp
         ENDIF

         !=====================================================================
         ! (2)  I n t e r n a l   T i m e   S t e p   L o o p
         !=====================================================================
         DO ISTEP = 1, NS

            ! Initialize
            QC     = 0e+0_fp    ! [kg species/kg dry air]                  
            T0_SUM = 0e+0_fp    ! [kg species/m2/timestep]

            !----------------------------------------------------------
            ! B e l o w   C l o u d   B a s e   (K < CLDBASE) 
            !
            ! QB is the "weighted avg" mixing ratio below the cloud
            ! base [kg/kg dry air].
            ! QC is the mixing ratio of the air that moved in cumulus 
            ! transport up to the next level [kg/kg dry air]. 
            ! MB is the dry mass of air below the cloud base per 
            ! unit area [kg/m2] (see calculation before loop).
            !-----------------------------------------------------------

            ! We need to make this a nested IF statement so that we don't
            ! get an out-of-bounds error when CLDBASE=1 (bmy, 11/18/10)
            IF ( CLDBASE > 1 ) THEN
               
               IF ( CMFMC(CLDBASE-1) > TINYNUM ) THEN
 
                  !-----------------------------------------------------
                  ! %%% Non-negligible Cloud mass flux %%% 
                  !-----------------------------------------------------

                  ! Calculate QB_NUM, the numerator for QB. QB is the 
                  ! weighted average mixing ratio below the cloud base.
                  ! QB_NUM is equal to the grid box species concentrations 
                  ! [kg/kg dry air] weighted by the adjacent level pressure 
                  ! differences and summed over all levels up to just 
                  ! below the cloud base (ewl, 6/22/15)           
                  QB_NUM  = 0e+0_fp  
                  DELP_DRY_NUM = 0e+0_fp

                  DO K  = 1, CLDBASE-1
                     QB_NUM = QB_NUM + Q(K,IC) * DELP_DRY(K)
                     DELP_DRY_NUM = DELP_DRY_NUM + DELP_DRY(K)
                  ENDDO

                  ! Compute QB, the weighted avg mixing ratio below
                  ! the cloud base [kg/kg dry air]
                  QB = QB_NUM / DELP_DRY_NUM

                  ! Compute QC, the mixing ratio of the air that moved
                  ! in cumulus transport up to the next level [kg/kg]
                  !
                  !        Dry mass of species below cloud base  +   
                  !        Subsidence into cloud base from above  
                  ! QC =  --------------------------------------------
                  !            Dry air mass below cloud base 
                  !
                  QC = ( MB*QB + CMFMC(CLDBASE-1) * 
     &                           Q(CLDBASE,IC)    * SDT  ) /
     &                 ( MB    + CMFMC(CLDBASE-1) * SDT  ) 

                  ! Copy QC to all levels of the species array Q
                  ! that are below the cloud base level [kg/kg]
                  Q(1:CLDBASE-1,IC) = QC

               ELSE                  

                  !-----------------------------------------------------
                  ! %%% Negligible cloud mass flux %%% 
                  !-----------------------------------------------------

                  ! When CMFMC is negligible, then set QC to the species
                  ! concentration at the cloud base level [kg/kg]
                  QC = Q(CLDBASE,IC)
   
               ENDIF

            ELSE
               
               !-----------------------------------------------------
               ! If the cloud base happens at level 1, then just
               ! set QC to the species concentration at the surface 
               ! level [kg/kg]
               !-----------------------------------------------------
               QC = Q(CLDBASE,IC) 

            ENDIF

            !==================================================================
            ! (3)  A b o v e   C l o u d   B a s e
            !==================================================================
            DO K = CLDBASE, KTOP
               
               ! Initialize
               ALPHA   = 0e+0_fp
               ALPHA2  = 0e+0_fp
               CMOUT   = 0e+0_fp
               ENTRN   = 0e+0_fp
               QC_PRES = 0e+0_fp
               QC_SCAV = 0e+0_fp
               
               ! CMFMC_BELOW is the air mass [kg/m2/s] coming into the
               ! grid box (K) from the box immediately below (K-1).
               IF ( K == 1 ) THEN
                  CMFMC_BELOW = 0e+0_fp
               ELSE
                  CMFMC_BELOW = CMFMC(K-1)
               ENDIF

               ! If we have a nonzero air mass flux coming from 
               ! grid box (K-1) into (K) ...
               IF ( CMFMC_BELOW > TINYNUM ) THEN

                  !------------------------------------------------------------
                  ! (3.1)  M a s s   B a l a n c e   i n   C l o u d
                  !
                  ! F(K,NA) = fraction of species IC in level K that is 
                  !           available for wet-scavenging by cloud updrafts.  
                  !
                  ! If ENTRN > 0 then compute the new value of QC:
                  !
                  !      species mass from below      (i.e. level K-1) + 
                  !      species mass from this level (i.e. level K)
                  !  = -----------------------------------------------------
                  !             dry mass coming into cloud
                  !
                  ! Otherwise, preserve the previous value of QC.  This will 
                  ! ensure that TERM1 - TERM2 is not a negative quantity (see 
                  ! below).
                  !  
                  ! Entrainment must be >= 0 (since we cannot have a negative 
                  ! flux of air into the cloud).  This condition is strong 
                  ! enough to ensure that CMOUT > 0 and will prevent floating-
                  ! point exception.
                  !------------------------------------------------------------

                  ! Air mass flowing out of cloud at grid box (K) [kg/m2/s]
                  CMOUT   = CMFMC(K) + DTRAIN(K)

                  ! Air mass flowing into cloud at grid box (K) [kg/m2/s]
                  ENTRN   = CMOUT - CMFMC_BELOW

                  ! Amount of QC preserved against scavenging [kg/kg]
!                  QC_PRES = QC * ( 1e+0_fp - F(K,IC) )
                  QC_PRES = QC * ( 1e+0_fp - F(K,NA) )

                  ! Amount of QC lost to scavenging [kg/kg]
                  ! QC_SCAV = 0 for non-soluble species
!                  QC_SCAV = QC * F(K,IC)
                  QC_SCAV = QC * F(K,NA)

                  ! - - - - - - - - FOR SOLUBLE SPECIES ONLY - - - - - - - - - 
                  IF ( QC_SCAV > 0e+0_fp ) THEN 

                     ! The fraction ALPHA is the fraction of raindrops that 
                     ! will re-evaporate soluble species while falling from 
                     ! grid box K+1 down to grid box K.  Avoid div-by-zero.

                     ! Initialize 
                     ALPHA = 0e+0_fp
                     
                     IF ( PDOWN(K+1)  > TINYNUM ) THEN 

                        ! %%%% CASE 1 %%%%
                        ! Partial re-evaporation. Less precip is leaving 
                        ! the grid box then entered from above.
                        IF ( PDOWN(K+1) > PDOWN(K) .AND.
     &                       PDOWN(K)   > TINYNUM        ) THEN

                           ! Define ALPHA, the fraction of raindrops that 
                           ! re-evaporate when falling from grid box
                           ! (I,J,L+1) to (I,J,L):
                           ! NOTE:
                           !   REEVAPCN is in units of [kg/kg/s]
                           !   Now use BMASS [kg/m2] instead of AD/area to
                           !   remove area dependency
                           !   PDOWN is in units of [cm3/cm2/s]
                           !   Factor of 10 in denom for unit conversion
                           !     1000 kg/m3 * 0.01 m/cm = 10 kg/m2/cm
                           ALPHA = REEVAPCN(K) * BMASS(K) 
     &                             / ( PDOWN(K+1)  * 10e+0_fp )

                           ! Restrict ALPHA to be less than 1 
                           ! (>1 is unphysical)  (hma, 24-Dec-2010)  
                           IF ( ALPHA > 1e+0_fp ) THEN 
                              ALPHA = 1e+0_fp 
                           ENDIF

                           ! We assume that 1/2 of the soluble species w/in 
                           ! the raindrops actually gets resuspended into 
                           ! the atmosphere
                           ALPHA2   = ALPHA * 0.5e+0_fp

                        ENDIF

                        ! %%%% CASE 2 %%%%
                        ! Total re-evaporation. Precip entered from above, 
                        ! but no precip is leaving grid box (ALPHA = 2 so 
                        ! that  ALPHA2 = 1)
                        IF ( PDOWN(K) < TINYNUM ) THEN
                           ALPHA2 = 1e+0_fp
                        ENDIF
                        
                     ENDIF
                    
                     ! The resuspension takes 1/2 the amount of the scavenged
                     ! aerosol (QC_SCAV) and adds that back to QC_PRES
                     QC_PRES  = QC_PRES + ( ALPHA2 * QC_SCAV )
                    
                     ! ... then we decrement QC_SCAV accordingly
                     QC_SCAV  = QC_SCAV * ( 1e+0_fp    - ALPHA2     )

                  ENDIF
                  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

                  ! Update QC taking entrainment into account [kg/kg]
                  ! Prevent div by zero condition
                  IF ( ENTRN >= 0e+0_fp .and. CMOUT > 0e+0_fp ) THEN
                     QC   = ( CMFMC_BELOW * QC_PRES   + 
     &                       ENTRN       * Q(K,IC) ) / CMOUT
                  ENDIF

                  !------------------------------------------------------------
                  ! (3.2)  M a s s   B a l a n c e   i n   L e v e l  ==> Q
                  !
                  ! Terminology:
                  !
                  !  C_k-1   = cloud air mass flux from level k-1 to level k
                  !  C_k     = cloud air mass flux from level k   to level k+1
                  !  QC_k-1  = mixing ratio of species INSIDE CLOUD at level k-1
                  !  QC_k    = mixing ratio of species INSIDE CLOUD at level k
                  !  Q_k     = mixing ratio of species in level k
                  !  Q_k+1   = mixing ratio of species in level k+1
                  ! 
                  ! For convenience we denote:
                  !
                  !  QC_SCAV = Amount of species wet-scavenged in updrafts
                  !          = QC_k-1 * F(k,IC)    [kg/kg]
                  !
                  !  QC_PRES = Amount of species preserved against
                  !            wet-scavenging in updrafts [kg/kg]
                  !          = QC_k-1 * ( 1 - F(k,IC) )   
                  !
                  ! Where F(k,IC) is the fraction of species IC in level k
                  ! that is available for wet-scavenging by cloud updrafts.
                  ! F(k,IC) is computed by routine COMPUTE_UPDRAFT_FSOL
                  ! and passed to this routine as an argument.
                  !
                  ! The cumulus transport above the cloud base is done as 
                  ! follows:
                  !  
                  !                 ||///////////////////||
                  !                 ||//// C L O U D ////||
                  !                 ||                   ||
                  !   k+1     ^     ||         ^         ||3)   C_k * Q_k+1
                  !           |     ||         |         ||         |
                  !   --------|-----++---------|---------++---------|--------
                  !           |     ||         |         ||         |
                  !   k      C_k    ||2)   C_k * QC_k    ||         V
                  !                 ||                   ||
                  !                 ||                   ||
                  !           ^     ||         ^         ||4)   C_k-1 * Q_k 
                  !           |     ||         |         ||         |
                  !   --------|-----++---------|---------++---------|--------
                  !           |     ||         |         ||         |
                  !   k-1   C_k-1   ||1) C_k-1 * QC_k-1  ||         V
                  !                 ||         * (1 - F) ||
                  !                 ||                   ||
                  !                 ||//// C L O U D ////||
                  !                 ||///////////////////||
                  !
                  ! There are 4 terms that contribute to mass flow in 
                  ! and out of level k:
                  !
                  ! 1) C_k-1 * QC_PRES = species convected from k-1 to k 
                  ! 2) C_k   * QC_k    = species convected from k   to k+1 
                  ! 3) C_k   * Q_k+1   = species subsiding from k+1 to k 
                  ! 4) C_k-1 * Q_k     = species subsiding from k   to k-1 
                  !
                  ! Therefore the change in species concentration is given by
                  !
                  !    DELQ = (Term 1) - (Term 2) + (Term 3) - (Term 4)
                  !
                  ! and Q(K,IC) = Q(K,IC) + DELQ.  
                  !
                  ! The term T0 is the amount of species that is scavenged 
                  ! out of the box.
                  !
                  ! Units of T0, T1, T2, T3, T4, and TSUM are
                  ! [kg/m2/s * kg species / kg dry air]
                  !------------------------------------------------------------
                  T0      =  CMFMC_BELOW * QC_SCAV      
                  T1      =  CMFMC_BELOW * QC_PRES       
                  T2      = -CMFMC(K  )  * QC
                  T3      =  CMFMC(K  )  * Q(K+1,IC)
                  T4      = -CMFMC_BELOW * Q(K,  IC)
                 
                  TSUM    = T1 + T2 + T3 + T4          
                  
                  DELQ    = ( SDT / BMASS(K) ) * TSUM    ! change in [kg/kg]

                  ! If DELQ > Q then do not make Q negative!!!
                  IF ( Q(K,IC) + DELQ < 0 ) THEN
                     DELQ = -Q(K,IC)
                  ENDIF

                  ! Increment the species array [kg/kg]
                  Q(K,IC) = Q(K,IC) + DELQ

                  ! check for infinity
                  IF ( .not. IT_IS_FINITE( Q(K,IC) ) ) THEN
                     CALL GEOS_CHEM_STOP
                  ENDIF

                  ! Return if we encounter NaN
                  IF ( IT_IS_NAN( Q(K,IC) ) ) THEN 
                     WRITE( 6, 250 )
                     WRITE( 6, 255 ) K, IC, Q(K,IC) 
 250                 FORMAT( 'NaN encountered in DO_MERRA_CONVECTION!' )
 255                 FORMAT( 'K, IC, Q(K,IC): ', 2i4, 1x, es13.6 )
                     RC = GC_FAILURE
                     RETURN
                  ENDIF

                  ! Pass T0_SUM in units of [kg species/m2/timestep].
                  ! Converting kg dry air to kg species requires use
                  ! of the molecular weight of air including moisture
                  ! (ewl, 6/5/15)
                  T0_SUM = T0_SUM + T0 * SDT

                  !------------------------------------------------------------
                  ! (3.3)  N D 1 4   D i a g n o s t i c
                  !
                  ! Archive upward mass flux due to wet convection  
                  ! [kg/sec] in the box (I,J), for the species IC going out 
                  ! of the top of the layer K to the layer above (K+1)  
                  ! (bey, 11/10/99). We must divide by DNS, the # of internal 
                  ! timesteps so that the sum represents the average loss rate 
                  ! across all internal timesteps.
                  !------------------------------------------------------------
                  IF ( USE_DIAG14 ) THEN
                     DIAG14(K,NW) = 
     &               DIAG14(K,NW) + ( -T2-T3 ) * AREA_M2 / DNS   
                  ENDIF

                  !------------------------------------------------------------
                  ! (3.4)  N D 3 8   D i a g n o s t i c
                  !
                  ! Archive the loss of soluble species to wet scavenging in 
                  ! cloud updrafts [kg/s].  We must divide by DNS, the # of 
                  ! internal timesteps so that the sum represents the
                  ! average loss rate across all internal timesteps.
                  !------------------------------------------------------------
                  IF ( USE_DIAG38 .and. NW > 0 ) THEN
                     DIAG38(K,NW) = 
     &               DIAG38(K,NW) + ( T0 * AREA_M2 / DNS )

                     ! check for infinity (added by hma, 20101117)
                     IF ( .not. IT_IS_FINITE( DIAG38(K,NW) ) ) THEN
                        PRINT*, 'DIAG38 IS INFINITY at K = ', K
                        CALL GEOS_CHEM_STOP
                     ENDIF
                  ENDIF
                  
               ELSE

                  !------------------------------------------------------------
                  ! (3.5)  N o   C l o u d   M a s s   F l u x   B e l o w 
                  !------------------------------------------------------------

                  ! If there is no cloud mass flux coming from below, set
                  ! QC to the species concentration at this level [kg/kg]
                  QC = Q(K,IC)
                  
                  ! Bug fix for the cloud base layer, which is not necessarily
                  ! in the boundary layer, and for MERRA, there could be 
                  ! "secondary convection" plumes - one in the PBL and another 
                  ! one not.  NOTE: T2 and T3 are the same terms as described 
                  ! in the above section.  (swu, 08/13/2007)
                  IF ( CMFMC(K) > TINYNUM ) THEN 

                     ! Species convected from K -> K+1 
                     ! [kg/m2/s * kg species/kg dry air]
                     T2   = -CMFMC(K) * QC

                     ! Species subsiding from K+1 -> K [kg/m2/s]
                     ! [kg/m2/s * kg species/kg dry air]
                     T3   =  CMFMC(K) * Q(K+1,IC)

                     ! Change in species concentration [kg/kg]
                     DELQ = ( SDT / BMASS(K) ) * (T2 + T3)

                     ! If DELQ > Q then do not make Q negative!!!
                     IF ( Q(K,IC) + DELQ < 0.0e+0_fp ) THEN 
                        DELQ = -Q(K,IC)
                     ENDIF
  
                     ! Add change in species to Q array [kg/kg]
                     Q(K,IC) = Q(K,IC) + DELQ

                  ENDIF
               ENDIF
            ENDDO     ! End of loop over levels above cloud base
            
            !==================================================================
            ! (4)  B e l o w   C l o u d   B a s e
            !==================================================================
            DO K = CLDBASE-1, 1, -1

               ! Initialize
               QDOWN       = 0e+0_fp
               F_WASHOUT   = 0e+0_fp
               WASHFRAC    = 0e+0_fp
               ALPHA       = 0e+0_fp
               ALPHA2      = 0e+0_fp
               GAINED      = 0e+0_fp
               WETLOSS     = 0e+0_fp
               LOST        = 0e+0_fp
               MASS_WASH   = 0e+0_fp
               MASS_NOWASH = 0e+0_fp
               K_RAIN      = 0e+0_fp
               
               ! Check if...
               ! (1) there is precip coming into box (I,J,K) from (I,J,K+1)
               ! (2) there is re-evaporation happening in grid box (I,J,K)
               ! (3) there is species to re-evaporate
               IF ( PDOWN(K+1)  > 0 .and.   
     &              REEVAPCN(K) > 0 .and. 
     &              T0_SUM      > 0        ) THEN

                  ! Compute F_WASHOUT, the fraction of grid box (I,J,L)
                  ! experiencing washout. First, convert units of PDOWN, 
                  ! the downward flux of precip leaving grid box (K+1) 
                  ! from [cm3 H20/cm2 area/s] to [cm3 H20/cm3 air/s]
                  ! by dividing by box height in cm
                  QDOWN = PDOWN(K+1) / ( BXHEIGHT(K+1) * 100e+0_fp  )

                  ! Compute K_RAIN and F_WASHOUT based on the flux of precip 
                  ! leaving grid box (K+1). 
                  K_RAIN   = LS_K_RAIN( QDOWN )
                  F_WASHOUT= LS_F_PRIME( QDOWN, K_RAIN ) 

                  ! Call WASHOUT to compute the fraction of species lost
                  ! to washout in grid box (I,J,K)
                  CALL WASHOUT( am_I_Root, I, J, K, IC, BXHEIGHT(K), 
     &                          T(K),      QDOWN,    SDT,   
     &                          F_WASHOUT, H2O2s(K), SO2s(K), 
     &                          WASHFRAC,  AER,      Input_Opt,
     &                          State_Met, State_Chm, RC         )

                  ! Check if the species is an aerosol or not
                  IF ( AER ) THEN
 
                     !---------------------------------------------------------
                     ! Washout of aerosol species
                     ! This is modeled as a kinetic process
                     !---------------------------------------------------------

                     ! Define ALPHA, the fraction of raindrops that 
                     ! re-evaporate when falling from (I,J,L+1) to (I,J,L)
                     ! NOTE:
                     !   REEVAPCN is in units of [kg/kg/s]
                     !   Now use BMASS [kg/m2] instead of AD/area to
                     !   remove area dependency
                     !   PDOWN is in units of [cm3/cm2/s]
                     !   Factor of 10 in denom for unit conversion
                     !     1000 kg/m3 * 0.01 m/cm = 10 kg/m2/cm
                      
                     ! %%%% CASE 1 %%%%
                     ! Partial re-evaporation. Less precip is leaving
                     ! the grid box then entered from above (V. Shah, 9/14/15)
                     IF ( PDOWN(K+1) > PDOWN(K) .AND.
     &                    PDOWN(K)   > TINYNUM        ) THEN

                        ! Define ALPHA, the fraction of raindrops that
                        ! re-evaporate when falling from grid box
                        ! (I,J,L+1) to (I,J,L)
                        ALPHA = REEVAPCN(K) * BMASS(K)
     &                          / ( PDOWN(K+1) * 10e+0_fp  )   

                        ! For safety
                        ALPHA = MIN( ALPHA, 1e+0_fp ) 

                        ! ALPHA2 is the fraction of the rained-out aerosols
                        ! that gets resuspended in grid box (I,J,L)
                        ALPHA2  = 0.5e+0_fp * ALPHA

                     ENDIF

                     ! %%%% CASE 2 %%%%
                     ! Total re-evaporation. Precip entered from above,
                     ! but no precip is leaving grid box (ALPHA = 2 so
                     ! that  ALPHA2 = 1) (V. Shah, 9/14/15)
                     IF ( PDOWN(K) < TINYNUM ) THEN
                        ALPHA2 = 1e+0_fp
                     ENDIF

                     ! GAINED is the rained out aerosol coming down from 
                     ! grid box (I,J,L+1) that will evaporate and re-enter 
                     ! the atmosphere in the gas phase in grid box (I,J,L)
                     ! [kg species/m2/timestep]
                     GAINED  = T0_SUM * ALPHA2

                     ! Amount of aerosol lost to washout in grid box [kg/m2]
                     ! (V. Shah, 9/14/15)
                     WETLOSS = ( Q(K,IC) * BMASS(K) + GAINED ) *
     &                         WASHFRAC - GAINED

                     ! LOST is the rained out aerosol coming down from
                     ! grid box (I,J,L+1) that will remain in the liquid
                     ! phase in grid box (I,J,L) and will NOT re-evaporate
                     ! [kg/m2/timestep]
                     LOST    = T0_SUM - GAINED

                     ! Update species concentration (V. Shah, mps, 5/20/15)
                     ! [kg/kg]
                     Q(K,IC) = Q(K,IC) - WETLOSS / BMASS(K) 

                     ! Update T0_SUM, the total amount of scavenged
                     ! species that will be passed to the grid box below
                     ! [kg/m2/timestep]
                     T0_SUM = T0_SUM + WETLOSS  

                  ELSE

                     !---------------------------------------------------------
                     ! Washout of non-aerosol species
                     ! This is modeled as an equilibrium process
                     !---------------------------------------------------------

                     ! MASS_NOWASH is the amount of non-aerosol species in 
                     ! grid box (I,J,L) that is NOT available for washout.
                     ! Calculate in units of [kg/kg]
                     MASS_NOWASH = ( 1e+0_fp - F_WASHOUT ) * Q(K,IC)

                     ! MASS_WASH is the total amount of non-aerosol species
                     ! that is available for washout in grid box (I,J,L).
                     ! It consists of the mass in the precipitating
                     ! part of box (I,J,L), plus the previously rained-out
                     ! species coming down from grid box (I,J,L+1).
                     ! (Eq. 15, Jacob et al, 2000)
                     ! Units are [kg species/m2/timestep]
                     MASS_WASH = ( F_WASHOUT * Q(K,IC) ) * BMASS(K)
     &                            + T0_SUM

                     ! WETLOSS is the amount of species mass in 
                     ! grid box (I,J,L) that is lost to washout.
                     ! (Eq. 16, Jacob et al, 2000)
                     ! [kg species/m2/timestep]
                     WETLOSS     = MASS_WASH * WASHFRAC - T0_SUM
                     
                     ! The species left in grid box (I,J,L) is what was
                     ! originally in the non-precipitating fraction 
                     ! of the box, plus MASS_WASH, less WETLOSS. 
                     ! [kg/kg]
                     Q(K,IC) = Q(K,IC) - WETLOSS / BMASS(K) 
                     
                     ! Update T0_SUM, the total scavenged species
                     ! that will be passed to the grid box below
                     ! [kg species/m2/timestep]
                     T0_SUM      = T0_SUM + WETLOSS

                  ENDIF

                  !------------------------------------------------------------
                  ! N D 1 4   D i a g n o s t i c
                  !
                  ! Archive upward mass flux due to wet convection.  
                  ! [kg/sec] in the box (I,J), for the species IC going 
                  ! out of the top of the layer K to the layer above (K+1)  
                  ! (bey, 11/10/99). We must divide by DNS, the # of internal 
                  ! timesteps so that the sum represents the average loss 
                  ! rate across all internal timesteps.
                  !------------------------------------------------------------
                  IF ( USE_DIAG14 .and. NW > 0 ) THEN
                     DIAG14(K,NW) = 
     &               DIAG14(K,NW) + ( -T2-T3 ) * AREA_M2 / DNS 
                  ENDIF

                  !------------------------------------------------------------
                  !  N D 3 8   D i a g n o s t i c
                  !
                  ! Archive the loss of soluble species to wet scavenging in 
                  ! cloud updrafts [kg/s].  We must divide by NDT, the # of 
                  ! seconds in the convective timestep, equal to DNS * SDT, 
                  ! in order to make diag38 represent the average loss rate 
                  ! across all internal timesteps. Note that the units of
                  ! WETLOSS are [kg/m2/timestep].
                  !------------------------------------------------------------
                  !%%% NOTE: SHOULD TEST FOR NW > 0 BUT IF WE DO THAT WE
                  !%%% NO LONGER GET IDENTICAL RESULTS WITH THE REF CODE.
                  !%%% LOOK INTO THIS LATER.  (bmy, 7/7/16)
                  IF ( USE_DIAG38 .and. F(K,NA) > 0.0_fp ) THEN
                     DIAG38(K,NW) = 
     &               DIAG38(K,NW) + ( WETLOSS * AREA_M2 / NDT )
                  ENDIF
                
                  ! CHECK for infinity (added by hma, 20101117)
                  IF ( .not. IT_IS_FINITE( DIAG38(K,NW) ) ) THEN
                     PRINT*, 'DIAG38 IS INFINITY at K = ', K
                     CALL GEOS_CHEM_STOP
                  ENDIF
               ENDIF
            ENDDO     ! End of loop over levels below cloud base  

            !==================================================================
            ! (5)  M e r c u r y   O c e a n   M o d e l   A r c h i v a l
            !
            ! Pass the amount of Hg2 and HgP lost in wet  scavenging [kg] to 
            ! "ocean_mercury_mod.f" via ADD_Hg2_WET and ADD_HgP_WET.   We must 
            ! also divide  by DNS, the # of internal timesteps. 
            ! (sas, bmy, eck, eds, 1/19/05, 1/6/06, 7/30/08)
            !
            ! NOTE: Reorder 
            !==================================================================
            IF ( Is_Hg ) THEN

               !--------------------------------------
               ! Hg2 
               !--------------------------------------
               IF ( SpcInfo%IS_Hg2 ) THEN

                  ! Wet scavenged Hg(II) in [kg]
                  WET_Hg2 = ( T0_SUM * AREA_M2 )

                  ! Category # for this Hg2 species
                  Hg_Cat  = SpcInfo%Hg_Cat

                  ! Pass to "ocean_mercury_mod.f"
                  CALL ADD_Hg2_WD      ( I, J, Hg_Cat, WET_Hg2   )
                  CALL ADD_Hg2_SNOWPACK( I, J, Hg_Cat, WET_Hg2, 
     &                                                 State_Met )
               ENDIF

               !--------------------------------------
               ! HgP
               !--------------------------------------
               IF ( SpcInfo%Is_HgP ) THEN

                  ! Wet scavenged Hg(P) in [kg]
                  WET_HgP = ( T0_SUM * AREA_M2 )

                  ! Category # for this Hg2 species
                  Hg_Cat  = SpcInfo%Hg_Cat

                  ! Pass to "ocean_mercury_mod.f"
                  CALL ADD_HgP_WD      ( I, J, Hg_Cat, WET_HgP   )
                  CALL ADD_Hg2_SNOWPACK( I, J, Hg_Cat, WET_HgP, 
     &                                                 State_Met )
               ENDIF
            ENDIF
         ENDDO               ! End internal timestep loop

         ! Free pointer
         SpcInfo => NULL()
      ENDDO                  ! End loop over advected species

      !================================================================
      ! Succesful return!
      !================================================================
      
      ! Nullify pointers
      NULLIFY( BXHEIGHT )
      NULLIFY( CMFMC    )
      NULLIFY( DQRCU    )
      NULLIFY( DTRAIN   )
      NULLIFY( PFICU    )
      NULLIFY( PFLCU    )
      NULLIFY( REEVAPCN )
      NULLIFY( DELP_DRY )
      NULLIFY( T        )
      NULLIFY( H2O2s    )
      NULLIFY( SO2s     )
      NULLIFY( Q        )

      ! Set error code to success
      RC                      = GC_SUCCESS

      END SUBROUTINE DO_MERRA_CONVECTION
!EOC
      END MODULE CONVECTION_MOD
