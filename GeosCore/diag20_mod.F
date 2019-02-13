!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag20_mod.F
!
! !DESCRIPTION: Module DIAG20\_MOD contains variables and routines which are 
!  used to compute the production and loss of O3 for use in the tagged O3
!  simulation
!\\
!\\
! !INTERFACE:
!
      MODULE DIAG20_MOD
!
! !USES:
!
      USE PRECISION_MOD

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      REAL(fp), ALLOCATABLE, PUBLIC :: POx(:,:,:) ! Ox prod [molec/cm3/s]
      REAL(fp), ALLOCATABLE, PUBLIC :: LOx(:,:,:) ! Ox loss [molec/cm3/s]
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: DIAG20
      PUBLIC  :: CLEANUP_DIAG20
      PUBLIC  :: INIT_DIAG20
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: ITS_TIME_FOR_WRITE20
      PRIVATE :: WRITE20
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Add TAUe as a module variable.  Bug fixes: Make sure WRITE20 uses the 
!        global FILENAME, and also write to disk on the last timestep before
!        the end of the simulation. (bmy, 11/15/04)
!  (2 ) Added routine ITS_TIME_FOR_WRITE20 (bmy, 3/3/05)
!  (3 ) Added functions GET_NFAM, GET_FAM_MWT, GET_FAM_NAME (bmy, 5/2/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (6 ) Bug fix in DIAG20 (phs, 1/22/07)
!  (7 ) Now use LD65 as the vertical dimension instead of LLTROP or LLTROP_FIX
!        in DO_DIAG_PL, DIAG20, and WRITE20 (phs, bmy, 12/4/07)
!  (8 ) Now make COUNT a 3-D array (phs, 11/18/08)
!  (9 ) Minor fix in DIAG20 (dbj, bmy, 10/26/09)
!  (10) Make public FAM_NAME and H2SO4RATE (win, 1/25/10)  
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers
!  06 Aug 2012 - R. Yantosca - Make IU_ND20 a local module variable
!  12 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  07 Jul 2015 - R. Yantosca - Added fixes for several minor issues
!  09 Aug 2016 - M. Sulprizio- Moved code specific to ND20 from obsolete
!                              diag_pl_mod.F to new module diag20_mod.F
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! Scalars
      INTEGER                        :: IU_ND20
      INTEGER                        :: NFAM
      INTEGER                        :: YYYYMMDD
      REAL(f8)                       :: TAUb, TAUe, TAU0, TAU1
      CHARACTER(LEN=255)             :: FILENAME

      ! Arrays
      INTEGER,           ALLOCATABLE :: COUNT(:,:,:  )
      REAL(fp),          ALLOCATABLE :: PL24H(:,:,:,:)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diag20
!
! !DESCRIPTION: Subroutine DIAG20 computes production and loss rates of O3, 
!  and then calls subroutine WRITE20 to save the these rates to disk.  By 
!  saving the production and loss rates from a full-chemistry run,
!  a user can use these archived rates to perform a quick O3 chemistry
!  run at a later time.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DIAG20( am_I_Root, Input_Opt, State_Chm, State_Met,
     &                   RC )
!
! !USES:
!
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD
#endif
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PhysConstants,      ONLY : AVO
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : EXPAND_DATE,   GET_NYMD
      USE TIME_MOD,           ONLY : GET_TAU,       GET_TAUb 
      USE TIME_MOD,           ONLY : ITS_A_NEW_DAY, TIMESTAMP_STRING
      USE TIME_MOD,           ONLY : GET_LOCALTIME
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
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
!  If you want to archive the prod and loss of ozone from one of the full-
!  chemistry simulations (e.g. tropchem, UCX, SOA), then the PROD & LOSS MENU
!  in your input.geos file should look similar to this:
!                                                                             .
!      %%% PROD & LOSS MENU %%%:
!      Turn on P/L (ND65) diag?: T
!      # of levels for ND65    : 47
!      Save O3 P/L (ND20)?     : T
!                                                                             .
!  POx and LOx values for ND20 are obtained from KPP in flexchem_mod.F.
!
!  ##########################################################################
!  #####    NOTE: BINARY PUNCH INPUT IS BEING PHASED OUT.  THIS DATA    #####
!  #####    WILL EVENTUALLY BE WRITTEN TO netCDF FILES FOR HEMCO!       #####
!  #####       -- Bob Yantosca (05 Mar 2015)                            #####
!  ##########################################################################
!
! !REVISION HISTORY:
!  09 Jun 1999 - I. Bey      - Initial version
!  (1 ) Now bundled into "diag20_mod.f" (bmy, 7/20/04)
!  (2 ) Now also write to disk when it is the last timestep before the end of 
!        the run.  Now references GET_TAUE from "time_mod.f". (bmy, 11/15/04)
!  (3 ) Now call function ITS_TIME_FOR_WRITE20 to determine if the next
!        chemistry timestep is the start of a new day.  Remove reference
!        to GET_TAUe and GET_TS_CHEM.  Now archive P(Ox) and L(Ox) first
!        and then test if we have to save the file to disk. (bmy, 3/3/05)
!  (4 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (5 ) Now use LLTROP_FIX instead of LLTROP (phs, 1/22/07)
!  (6 ) Now use LD65 instead of LLTROP_FIX (phs, bmy, 12/4/07)
!  (7 ) Now take care of boxes that switch b/w stratospheric and tropospheric
!        regimes (phs, 11/17/08)
!  (8 ) Bug fix: Now just zero arrays w/o loop indices (dbj, bmy, 10/26/09)
!  15 Sep 2010 - R. Yantosca - Added ProTeX headers 
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  14 Mar 2013 - M. Payer    - Replace Ox with O3 as part of removal of
!                              NOx-Ox partitioning
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object + RC via the arg list
!  04 Apr 2013 - R. Yantosca - Now pass the Input_Opt object
!  06 Jul 2015 - R. Yantosca - Zero P_Ox and L_Ox variables at start of loop
!  06 Jul 2015 - R. Yantosca - Now skip computations if we are not in the
!                              chemgrid (where JLOOP == 0)
!  08 Jul 2015 - R. Yantosca - Now save POx as kg/m3/s and LOx as 1/m3/s,
!                              for compatibility with HEMCO
!  31 May 2016 - E. Lundgren - Use molec wt from species database rather than
!                              Input_Opt%XNUMOL 
!  09 Aug 2016 - M. Sulprizio- Move routine from diag_pl_mod.F to diag20_mod.F;
!                              Obtain prod/loss rates from KPP via the POx and
!                              LOx arrays set in flex_chemdr.F
!  10 Aug 2016 - M. Sulprizio- Remove State_Chm%Tracers from LOx calculation and
!                              replace with State_chm%Species. Convert from
!                              molec/cm3 to molec using BOXVL.
!  10 Apr 2018 - R. Yantosca - Now divide LOx by the mass of Ox, not O3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! SAVEd scalars
      LOGICAL,  SAVE      :: FIRST = .TRUE.
      REAL(fp), SAVE      :: O3_MW_g
      INTEGER,  SAVE      :: id_O3,      id_NO2,    id_NO3,    id_PAN
      INTEGER,  SAVE      :: id_NPMN,    id_IPMN,   id_PPN,    id_HNO4
      INTEGER,  SAVE      :: id_N2O5,    id_HNO3,   id_BrO,    id_HOBr
      INTEGER,  SAVE      :: id_BrNO2,   id_BrNO3,  id_MPN,    id_ETHLN
      INTEGER,  SAVE      :: id_ISN1,    id_ISOPNB, id_ISOPND, id_MACRN
      INTEGER,  SAVE      :: id_MVKN,    id_PROPNN, id_R4N2,   id_INPN
      INTEGER,  SAVE      :: id_ISNP,    id_INO2,   id_ISNOOA, id_ISNOOB
      INTEGER,  SAVE      :: id_ISNOHOO, id_MAN2,   id_PRN1,   id_PRPN
      INTEGER,  SAVE      :: id_R4N1,    id_PMNN,   id_MACRNO2,id_ClO
      INTEGER,  SAVE      :: id_HOCl,    id_ClNO2,  id_ClNO3,  id_Cl2O2
      INTEGER,  SAVE      :: id_OClO,    id_O,      id_O1D,    id_IO
      INTEGER,  SAVE      :: id_HOI,     id_IONO,   id_IONO2,  id_OIO
      INTEGER,  SAVE      :: id_I2O2,    id_I2O3,   id_I2O4

      ! Scalars
      LOGICAL             :: DO_WRITE
      INTEGER             :: I, J, L, N
      REAL(fp)            :: P_Ox, L_Ox, Ox_Mass
      REAL(fp)            :: LT
      REAL(fp)            :: BOXVL

      ! Strings
      CHARACTER(LEN=16)   :: STAMP 
!
! !DEFINED PARAMETERS:
!
      ! Local time limits (optional)
      REAL(fp), PARAMETER :: LT_START = 12.0e+0_fp
      REAL(fp), PARAMETER :: LT_END   = 16.0e+0_fp

      !=================================================================
      ! DIAG20 begins here!
      !=================================================================

      ! Assume success
      RC = GC_SUCCESS

#if defined( BPCH_DIAG )

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Get the species ID's
         id_O3      = Ind_('O3'     )
         id_NO2     = Ind_('NO2'    )
         id_NO3     = Ind_('NO3'    )
         id_PAN     = Ind_('PAN'    )
         id_NPMN    = Ind_('NPMN'   )
         id_IPMN    = Ind_('IPMN'   )
         id_PPN     = Ind_('PPN'    )
         id_HNO4    = Ind_('HNO4'   )
         id_N2O5    = Ind_('N2O5'   )
         id_HNO3    = Ind_('HNO3'   )
         id_BrO     = Ind_('BrO'    )
         id_HOBr    = Ind_('HOBr'   )
         id_BrNO2   = Ind_('BrNO2'  )
         id_BrNO3   = Ind_('BrNO3'  )
         id_MPN     = Ind_('MPN'    )
         id_ETHLN   = Ind_('ETHLN'  )
         id_ISN1    = Ind_('ISN1'   )
         id_ISOPNB  = Ind_('ISOPNB' )
         id_ISOPND  = Ind_('ISOPND' )
         id_MACRN   = Ind_('MACRN'  )
         id_MVKN    = Ind_('MVKN'   )
         id_PROPNN  = Ind_('PROPNN' )
         id_R4N2    = Ind_('R4N2'   )
         id_INPN    = Ind_('INPN'   )
         id_ISNP    = Ind_('ISNP'   )
         id_INO2    = Ind_('INO2'   )
         id_ISNOOA  = Ind_('ISNOOA' )
         id_ISNOOB  = Ind_('ISNOOB' )
         id_ISNOHOO = Ind_('ISNOHOO')
         id_MAN2    = Ind_('MAN2'   )
         id_PRN1    = Ind_('PRN1'   )
         id_PRPN    = Ind_('PRPN'   )
         id_R4N1    = Ind_('R4N1'   )
         id_PMNN    = Ind_('PMNN'   )
         id_MACRNO2 = Ind_('MACRNO2')
         id_ClO     = Ind_('ClO'    )
         id_HOCl    = Ind_('HOCl'   )
         id_ClNO2   = Ind_('ClNO2'  )
         id_ClNO3   = Ind_('ClNO3'  )
         id_Cl2O2   = Ind_('Cl2O2'  )
         id_OClO    = Ind_('OClO'   )
         id_O       = Ind_('O'      )
         id_O1D     = Ind_('O1D'    )
         id_IO      = Ind_('IO'     )
         id_HOI     = Ind_('HOI'    )
         id_IONO    = Ind_('IONO'   )
         id_IONO2   = Ind_('IONO2'  )
         id_OIO     = Ind_('OIO'    )
         id_I2O2    = Ind_('I2O2'   )
         id_I2O3    = Ind_('I2O3'   ) 
         id_I2O4    = Ind_('I2O4'   )     

         ! Error check
         IF ( id_O3 <= 0 ) THEN
            CALL ERROR_STOP( 'Ind_("O3") <= 0 !', 
     &                       'DIAG20 ("diag20_mod.f")' )
         ENDIF

         ! Get ozone molecular weight from species database
         O3_MW_g  = State_Chm%SpcData(id_O3)%Info%emMW_g ! g/mol

         ! Starting time of run
         TAUb     = GET_TAUb()

         ! Get time of run at 1st timestep
         TAU0     = TAUb

         ! Reset first-time flag
         FIRST    = .FALSE.

      ENDIF

      !=================================================================
      ! Archive P(Ox) and L(Ox) over the course of an entire day
      !=================================================================

      ! Echo info
      STAMP = TIMESTAMP_STRING()
      IF ( am_I_Root ) WRITE( 6, 120 ) STAMP
 120  FORMAT( '     - DIAG20: Archiving P(Ox) & L(Ox) at ', a )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, P_Ox, L_Ox, LT, BOXVL, Ox_Mass )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Zero P(Ox), L(Ox) variables for safety's sake
         P_Ox    = 0.0_fp
         L_Ox    = 0.0_fp
         Ox_Mass = 0.0_fp

         ! If this is a grid box where we have done chemistry, then ...
         IF ( State_Met%InChemGrid(I,J,L) ) THEN

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%% NOTE: If you want to restrict output to a given local
            !%%% time range, then uncomment these lines of code.
            !%%% (bmy, 7/9/15)
            !
            ! Get the local time at box (I,J,L)
            !LT = GET_LOCALTIME( I, J, L )
            !
            ! Skip processing if we are outside of the desired 
            ! local time range (bmy, 7/9/15)
            !IF ( LT < LT_START .or. LT > LT_END ) CYCLE
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ! Increment counter of valid grid boxes
            COUNT(I,J,L) = COUNT(I,J,L) + 1

            !-----------------------------------------------------------
            ! Production
            !-----------------------------------------------------------

            ! Convert P(Ox) from [molec/cm3/s] to [kg/cm3/s]
            P_Ox           = POx(I,J,L) / AVO * 1.e-3_fp *
     &                       O3_MW_g

            ! Convert to [kg/m3/s], for HEMCO (bmy, 7/7/15)
            P_Ox           = P_Ox * 1e+6_fp

            ! Store P(Ox) [kg/m3/s] in PL24H array
            PL24H(I,J,L,1) = PL24H(I,J,L,1) + P_Ox

            !-----------------------------------------------------------
            ! Loss
            !-----------------------------------------------------------

            ! Grid box volume [cm3]
            BOXVL          = State_Met%AIRVOL(I,J,L) * 1e+6_fp

            ! Compute the mass of Ox following the family definition
            ! in KPP/Tropchem/gckpp.kpp.  At this point, Ox_mass is 
            ! in units of [molec/cm3].
            !
            ! NOTE: If you update the Ox family definition, you will 
            ! also need to update this statement below.  
            !
            ! ALSO NOTE: We know that this is horribly kludged, but in 
            ! the interest of preparing v11-02 for release, we have 
            ! written this in the brute-force manner.  At a later time 
            ! we can examine how to make this more elegant.  This might 
            ! require modifying the KPP solver code, so we'll worry 
            ! about it later. (bmy, 4/10/18)
            !
            Ox_Mass        = State_Chm%Species(I,J,L,id_O3     ) 
     &                     + State_Chm%Species(I,J,L,id_NO2    ) 
     &                     + State_Chm%Species(I,J,L,id_NO3    )*2.0_fp
     &                     + State_Chm%Species(I,J,L,id_PAN    ) 
     &                     + State_Chm%Species(I,J,L,id_NPMN   ) 
     &                     + State_Chm%Species(I,J,L,id_IPMN   ) 
     &                     + State_Chm%Species(I,J,L,id_PPN    ) 
     &                     + State_Chm%Species(I,J,L,id_HNO4   ) 
     &                     + State_Chm%Species(I,J,L,id_N2O5   )*3.0_fp
     &                     + State_Chm%Species(I,J,L,id_HNO3   ) 
     &                     + State_Chm%Species(I,J,L,id_BrO    ) 
     &                     + State_Chm%Species(I,J,L,id_HOBr   ) 
     &                     + State_Chm%Species(I,J,L,id_BrNO2  ) 
     &                     + State_Chm%Species(I,J,L,id_BrNO3  )*2.0_fp
     &                     + State_Chm%Species(I,J,L,id_MPN    ) 
     &                     + State_Chm%Species(I,J,L,id_ETHLN  ) 
     &                     + State_Chm%Species(I,J,L,id_ISN1   ) 
     &                     + State_Chm%Species(I,J,L,id_ISOPNB ) 
     &                     + State_Chm%Species(I,J,L,id_ISOPND ) 
     &                     + State_Chm%Species(I,J,L,id_MACRN  ) 
     &                     + State_Chm%Species(I,J,L,id_MVKN   ) 
     &                     + State_Chm%Species(I,J,L,id_PROPNN ) 
     &                     + State_Chm%Species(I,J,L,id_R4N2   ) 
     &                     + State_Chm%Species(I,J,L,id_INPN   ) 
     &                     + State_Chm%Species(I,J,L,id_ISNP   ) 
     &                     + State_Chm%Species(I,J,L,id_INO2   ) 
     &                     + State_Chm%Species(I,J,L,id_ISNOOA ) 
     &                     + State_Chm%Species(I,J,L,id_ISNOOB ) 
     &                     + State_Chm%Species(I,J,L,id_ISNOHOO) 
     &                     + State_Chm%Species(I,J,L,id_MAN2   ) 
     &                     + State_Chm%Species(I,J,L,id_PRN1   ) 
     &                     + State_Chm%Species(I,J,L,id_PRPN   ) 
     &                     + State_Chm%Species(I,J,L,id_R4N1   ) 
     &                     + State_Chm%Species(I,J,L,id_PMNN   ) 
     &                     + State_Chm%Species(I,J,L,id_MACRNO2) 
     &                     + State_Chm%Species(I,J,L,id_ClO    ) 
     &                     + State_Chm%Species(I,J,L,id_HOCl   ) 
     &                     + State_Chm%Species(I,J,L,id_ClNO2  ) 
     &                     + State_Chm%Species(I,J,L,id_ClNO3  )*2.0_fp
     &                     + State_Chm%Species(I,J,L,id_Cl2O2  )*2.0_fp
     &                     + State_Chm%Species(I,J,L,id_OClO   )*2.0_fp
     &                     + State_Chm%Species(I,J,L,id_IO     ) 
     &                     + State_Chm%Species(I,J,L,id_HOI    ) 
     &                     + State_Chm%Species(I,J,L,id_IONO   ) 
     &                     + State_Chm%Species(I,J,L,id_IONO2  )*2.0_fp
     &                     + State_Chm%Species(I,J,L,id_OIO    )*2.0_fp
     &                     + State_Chm%Species(I,J,L,id_I2O2   )*2.0_fp
     &                     + State_Chm%Species(I,J,L,id_I2O3   )*3.0_fp
     &                     + State_Chm%Species(I,J,L,id_I2O4   )*4.0_fp

            ! The "Standard" KPP mechanism has all of the same species in
            ! its Ox family definition, plus O and O1D.  Add these too.
            IF ( Input_Opt%LUCX ) THEN
               Ox_Mass     = Ox_Mass
     &                     + State_Chm%Species(I,J,L,id_O      ) 
     &                     + State_Chm%Species(I,J,L,id_O1D    ) 
            ENDIF

            ! Convert Ox_Mass from [molec/cm3] to [molec]
            Ox_Mass        = Ox_Mass * BOXVL

            ! Divide L(Ox) [molec/cm3/s] by Ox mass [molec] 
            ! and then convert to [1/m3/s] for HEMCO
            IF ( Ox_Mass > 0.0_fp ) THEN
               L_Ox        = ( LOx(I,J,L) / Ox_Mass ) *  1.0e+6_fp
            ENDIF

            ! Store L(Ox) [1/m3/s] in PL24H array
            PL24H(I,J,L,2) = PL24H(I,J,L,2) + L_Ox
           
         ENDIF
 
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Write data to disk and zero counters for next timestep
      !=================================================================

      ! Check to see if the next chemistry timestep is the start of a
      ! new day, or the last chemistry timestep before the end-of-run.
      ! If so, then we need to write to disk. (bmy, 7/6/15)
      IF ( ITS_TIME_FOR_WRITE20( TAU1 ) ) THEN

         ! Compute average daily values
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, 2
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            IF ( COUNT(I,J,L) /= 0 ) THEN
               PL24H(I,J,L,N) = PL24H(I,J,L,N) / COUNT(I,J,L)
            ENDIF
         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! Get YYYYMMDD date for this day
         YYYYMMDD = GET_NYMD()        

         ! Replace YYYYMMDD in filename w/ the actual date
         FILENAME = 'rate.YYYYMMDD'
         CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )

         ! Echo info
         IF ( am_I_Root ) WRITE( 6, 110 ) TRIM( FILENAME )
 110     FORMAT( '     - DIAG20: Writing ', a )

         ! Write P(Ox) and L(Ox) to disk
         CALL WRITE20()

         ! Reset variables for the next diagnostic interval
         COUNT = 0
         PL24H = 0e+0_fp
         TAU0  = TAU1

      ENDIF
#endif

      END SUBROUTINE DIAG20
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write20
!
! !DESCRIPTION: Subroutine WRITE20 saves production and loss rates to disk, 
!  where they will be later read by subroutine CHEMO3. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE20()
!
! !USES:
!
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD      ! LD65
      USE BPCH2_MOD,   ONLY : BPCH2,         GET_HALFPOLAR
      USE BPCH2_MOD,   ONLY : GET_MODELNAME, OPEN_BPCH2_FOR_WRITE
#endif
      USE CMN_SIZE_MOD      ! Size parameters
      USE GC_GRID_MOD, ONLY : GET_XOFFSET,   GET_YOFFSET
      USE inquireMod,  ONLY : findFreeLUN
! 
! !REMARKS:
!  WRITE20 assumes that ND65 (P-L diagnostics) have been turned on.
!
!  ##########################################################################
!  #####    NOTE: BINARY PUNCH INPUT IS BEING PHASED OUT.  THIS DATA    #####
!  #####    WILL EVENTUALLY BE WRITTEN TO netCDF FILES FOR HEMCO!       #####
!  #####       -- Bob Yantosca (05 Mar 2015)                            #####
!  ##########################################################################
!
! !REVISION HISTORY:
!  09 Jun 1999 - I. Bey      - Initial version
!  (1 ) Now bundled into "diag20_mod.f" (bmy, 7/20/04)
!  (2 ) Bug fix: remove declaration of FILENAME which masked the global
!        declaration (bmy, 11/15/04)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now only write up to LD65 levels (phs, bmy, 12/4/07)
!  15 Sep 2010 - R. Yantosca - Added ProTeX headers 
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  06 Jul 2015 - R. Yantosca - Now do not use parallel loops for casting
!  06 Jul 2015 - R. Yantosca - Bug fix: restore missing JFIRST assignment
!  08 Jul 2015 - R. Yantosca - Now save out as per m3 instead of per cm3
!  09 Aug 2016 - M. Sulprizio- Move routine from diag_pl_mod.F to diag20_mod.F
!  10 Apr 2018 - R. Yantosca - Now save out levels 1:LLPAR instead of 1:LD65
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I, J, L, N, IOS
      INTEGER            :: IFIRST, JFIRST, LFIRST
      INTEGER            :: HALFPOLAR 
      INTEGER, PARAMETER :: CENTER180 = 1 
      REAL(f4)           :: LONRES, LATRES
      REAL(f4)           :: ARRAY(IIPAR,JJPAR,LLPAR)
      CHARACTER(LEN=20)  :: MODELNAME
      CHARACTER(LEN=40)  :: CATEGORY
      CHARACTER(LEN=40)  :: UNIT
      CHARACTER(LEN=40)  :: RESERVED
      CHARACTER(LEN=80)  :: TITLE

      !=================================================================
      ! WRITE20 begins here!
      !=================================================================

#if defined( BPCH_DIAG )

      ! Find a free file LUN
      IU_ND20   = findFreeLUN()

      ! Define various parameters for the BPCH file
      TITLE     = 'GEOS-CHEM archived P(Ox) and L(Ox) rates for Tag O3'
      CATEGORY  = 'PORL-L=$'
      RESERVED  = ''
      LONRES    = DISIZE
      LATRES    = DJSIZE
      MODELNAME = GET_MODELNAME()
      HALFPOLAR = GET_HALFPOLAR()
      IFIRST    = 1 + GET_XOFFSET( GLOBAL=.TRUE. )
      JFIRST    = 1 + GET_YOFFSET( GLOBAL=.TRUE. )
      LFIRST    = 1

      ! Open BPCH file for writing
      CALL OPEN_BPCH2_FOR_WRITE( IU_ND20, FILENAME, TITLE )

      !=================================================================
      ! Save P(Ox) to disk
      !=================================================================

      ! Zero for safety's sake
      ARRAY = 0e0

      ! Cast P(Ox) to REAL*4
      ARRAY = PL24H(:,:,:,1)

      ! Now save out as kg/m3/s, compatible with HEMCO (bmy, 7/8/15)
      UNIT = 'kg/m3/s'

      ! Save P(Ox) to BPCH file
      CALL BPCH2( IU_ND20,   MODELNAME, LONRES,   LATRES,    
     &            HALFPOLAR, CENTER180, CATEGORY, 1,           
     &            UNIT,      TAU0,      TAU1,     RESERVED,  
     &            IIPAR,     JJPAR,     LLPAR,    IFIRST,
     &            JFIRST,    LFIRST,    ARRAY               )

      !=================================================================
      ! Save L(Ox) to disk
      !=================================================================

      ! Zero for safety's sake
      ARRAY = 0e0

      ! Cast L(Ox) to REAL*4
      ARRAY = PL24H(:,:,:,2)

      ! Now save out as 1/m3/s, compatible with HEMCO (bmy, 7/8/15)
      UNIT = '1/m3/s'

      ! Save L(Ox) to BPCH file
      CALL BPCH2( IU_ND20,   MODELNAME, LONRES,   LATRES,    
     &            HALFPOLAR, CENTER180, CATEGORY, 2,           
     &            UNIT,      TAU0,      TAU1,     RESERVED,  
     &            IIPAR,     JJPAR,     LLPAR,    IFIRST,
     &            JFIRST,    LFIRST,    ARRAY               )

      ! Close BPCH file
      CLOSE( IU_ND20 )
#endif

      END SUBROUTINE WRITE20
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_write20
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_WRITE20 returns TRUE if it's time to 
!  write the ND20 ozone P/L rate file to disk.  We test the time at the next 
!  chemistry timestep so that we can write to disk properly. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_WRITE20( TAU_W ) RESULT( ITS_TIME )
!
! !USES:
!
      USE TIME_MOD, ONLY : GET_HOUR, GET_MINUTE, GET_SECOND,  GET_TAU
      USE TIME_MOD, ONLY : GET_TAUb, GET_TAUe,   GET_TS_CHEM, GET_TS_DYN
!
! !INPUT PARAMETERS: 
!
      REAL(f8), INTENT(OUT) :: TAU_W      ! TAU value @ time of writing to disk
!
! !RETURN VALUE:
!
      LOGICAL               :: ITS_TIME   ! =T if its time to write to disk
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  15 Sep 2010 - R. Yantosca - Added ProTeX headers 
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  06 Jul 2015 - R. Yantosca - Now use TAU+CHEM >= TAUe to test for end of run
!  09 Aug 2016 - M. Sulprizio- Move routine from diag_pl_mod.F to diag20_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: TAU, HOUR, CHEM, DYN

      !=================================================================
      ! ITS_TIME_FOR_WRITE20 begins here!
      !=================================================================

      ! Initialize
      ITS_TIME = .FALSE.

      ! Current TAU, Hour, and Dynamic Timestep [hrs]
      TAU      = GET_TAU()
      HOUR     = ( GET_SECOND()  / 3600e+0_fp ) + 
     &           ( GET_MINUTE()  / 60e+0_fp ) + GET_HOUR()
      CHEM     = ( GET_TS_CHEM() / 3600e+0_fp )
      DYN      = ( GET_TS_DYN()  / 3600e+0_fp )

      ! If first timestep, return FALSE
      IF ( TAU == GET_TAUb() ) RETURN

      ! If the next chemistry timestep is the hour of day
      ! when we have to save to disk, return TRUE
      IF ( HOUR + CHEM >= 24e+0_fp ) THEN
         ITS_TIME = .TRUE.
         TAU_W    = TAU + CHEM
         RETURN
      ENDIF

      ! If the next chem timestep is the end of the run, return TRUE.
      ! This allows us to do short debuggging runs (bmy, 7/6/15)
      IF ( TAU + CHEM >= GET_TAUe() ) THEN
         ITS_TIME = .TRUE.
         TAU_W    = TAU + CHEM
         RETURN
      ENDIF

      END FUNCTION ITS_TIME_FOR_WRITE20
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_diag20
!
! !DESCRIPTION: Subroutine INIT\_DIAG20 takes values read from the GEOS-Chem
!  input file and saves to module variables w/in "diag20\_mod.f" 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_DIAG20( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD
#endif
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS: 
!
      LOGICAL,           INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput),    INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,           INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY:
!  09 Aug 2016 - M. Sulprizio- Initial version
!  10 Apr 2018 - R. Yantosca - Allocate arrays 1:LLPAR instead of 1:LD65

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! INIT_DIAG20 begins here!
      !=================================================================

      ! Assume success
      RC = GC_SUCCESS

#if defined( BPCH_DIAG )
      ALLOCATE( COUNT( IIPAR, JJPAR, LLPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'COUNT' )
      COUNT = 0

      ALLOCATE( PL24H( IIPAR, JJPAR, LLPAR, 2 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'PL24H' )
      PL24H = 0e+0_fp

      ALLOCATE( POx( IIPAR, JJPAR, LLPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'POx' )
      POx = 0e+0_fp

      ALLOCATE( LOx( IIPAR, JJPAR, LLPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'LOx' )
      LOx = 0e+0_fp
#endif

      END SUBROUTINE INIT_DIAG20
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_diag20
!
! !DESCRIPTION: Subroutine CLEANUP\_DIAG20 deallocates all module arrays. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_DIAG20
! 
! !REVISION HISTORY: 
!  09 Aug 2016 - M. Sulprizio- Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_DIAG20 begins here!
      !=================================================================
#if defined( BPCH_DIAG )
      IF ( ALLOCATED( COUNT     ) ) DEALLOCATE( COUNT     )
      IF ( ALLOCATED( PL24H     ) ) DEALLOCATE( PL24H     )
      IF ( ALLOCATED( POx       ) ) DEALLOCATE( POx       )
      IF ( ALLOCATED( LOx       ) ) DEALLOCATE( LOx       )
#endif

      END SUBROUTINE CLEANUP_DIAG20
!EOC
      END MODULE DIAG20_MOD
