!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag42_mod.F
!
! !DESCRIPTION: Module DIAG42\_MOD contains arrays and routines for archiving 
!  the ND42 diagnostic -- secondary organic aerosols [ug/m3]. 
!\\
!\\
! !INTERFACE: 
!
      MODULE DIAG42_MOD
!
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !DEFINED PARAMETERS:
!
      ! Maximum number of output:
      INTEGER, PUBLIC, PARAMETER   :: PD42 = 25
!
! !PUBLIC DATA MEMBERS:
!
      INTEGER, PUBLIC              :: ND42            ! ND42 on/off flag
      INTEGER, PUBLIC              :: LD42            ! # of levels for ND42

      ! Arrays
      REAL*4,  PUBLIC, ALLOCATABLE :: AD42(:,:,:,:)   ! Array for SOA [ug/m3]
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC :: DIAG42
      PUBLIC :: ZERO_DIAG42
      PUBLIC :: WRITE_DIAG42
      PUBLIC :: INIT_DIAG42
      PUBLIC :: CLEANUP_DIAG42
!
! !REVISION HISTORY:
!  22 May 2006 - D. Henze, R. Yantosca - Initial version
!  (1 ) Replace TINY(1d0) with 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!  (2 ) Now use ratio of 2.1 instead of 1.4 for SOA4 (dkh, bmy, 3/29/07)
!  (3 ) Add diagnostics for SOAG and SOAM (tmf, 1/7/09)
!  (4 ) Increase PD42 to 24. (fp, hotp, 2/3/10)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations (H. Pye)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  10 Nov 2014 - M. Yannetti - Added PRECISION_MOD
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
! !IROUTINE: diag42
!
! !DESCRIPTION: Subroutine DIAG42 archives SOA concentrations [ug/m3] 
!  for the ND42 diagnostic.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DIAG42( am_I_Root, Input_Opt, State_Met, State_Chm, RC)
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
#if defined( BPCH_DIAG )
      USE State_Chm_Mod,      ONLY : Ind_
      USE AEROSOL_MOD,        ONLY : BCPI, BCPO, OCPI, OCPO, OCPISOA
      USE AEROSOL_MOD,        ONLY : SALA, SALC, SOILDUST
      USE AEROSOL_MOD,        ONLY : SO4,  NH4,  NIT
      USE AEROSOL_MOD,        ONLY : TSOA, ISOA, ASOA, OPOA
      USE AEROSOL_MOD,        ONLY : PM25
      USE AEROSOL_MOD,        ONLY : SOAGX,  SOAMG, ISOAAQ
      USE AEROSOL_MOD,        ONLY : OCFPOA, OCFOPOA
      USE CARBON_MOD,         ONLY : BETANOSAVE
      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! NDxx flags
      USE CMN_O3_MOD,         ONLY : SAVEOA
      USE ERROR_MOD
      USE UnitConv_Mod,       ONLY : Convert_Spc_Units
      USE PhysConstants,      ONLY : ATM, MWCARB
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS: 
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
! 
! !REVISION HISTORY: 
!  22 May 2006 - D. Henze, R. Yantosca - Initial version
!  (1 ) Now use ratio of 2.1 instead of 1.4 for SOA4 (dkh, bmy, 3/29/07)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations (H. Pye)
!  26 Feb 2015 - E. Lundgren - Remove dependency on pressure_mod (not used)
!  25 Mar 2015 - E. Lundgren - Change tracer units from kg to kg/kg
!  06 Jan 2016 - E. Lundgren - Use global physical parameter ATM
!  16 Jun 2016 - K. Yu       - Now define species ID's with the Ind_ function
!  17 Jun 2016 - R. Yantosca - Now only define species ID's on the first call
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  11 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!  01 Sep 2016 - M. Sulprizio- Add PM2.5 concentrations to index 18
!  18 Nov 2016 - M. Sulprizio- Move code for calculating concentrations to
!                              AEROSOL_CONC (in aerosol_mod.F) so it is in a
!                              single location
!  28 Sep 2017 - E. Lundgren - Simplify unit conversions with wrapper routine
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( BPCH_DIAG )

      INTEGER             :: I, J, L
      REAL(fp)            :: FACTOR
      CHARACTER(LEN=63)   :: OrigUnit
      CHARACTER(LEN=255)  :: MSG
     
      ! Species ID flags
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER, SAVE       :: id_POA1,   id_POA2
      INTEGER, SAVE       :: id_SOAIE,  id_SOAME, id_INDIOL
      INTEGER, SAVE       :: id_SOAGX,  id_SOAMG
      INTEGER, SAVE       :: id_LVOCOA, id_ISN1OA
      LOGICAL, SAVE       :: IS_BETANOSAVE
      LOGICAL, SAVE       :: IS_POA
      LOGICAL, SAVE       :: IS_ISOPOA

      ! Define number of carbon atoms in each irreversible isoprene
      ! SOA tracer species. Named according to the parent HC (same 
      ! number of carbons):
      REAL(fp), PARAMETER :: NCIMAE   = 4e+0_fp
      REAL(fp), PARAMETER :: NCIEPOX  = 5e+0_fp
      REAL(fp), PARAMETER :: NCINDIOL = NCIEPOX
      REAL(fp), PARAMETER :: NCGLYX   = 2e+0_fp
      REAL(fp), PARAMETER :: NCGLYC   = NCGLYX
      REAL(fp), PARAMETER :: NCMGLY   = 3e+0_fp
      REAL(fp), PARAMETER :: NCLVOC   = NCIEPOX
      REAL(fp), PARAMETER :: NCISN1   = NCIEPOX

      ! Pointers      
      REAL(fp), POINTER   :: Spc(:,:,:,:)
      REAL(fp), POINTER   :: AIRVOL(:,:,:)
#endif

      !================================================================= 
      ! DIAG42 begins here! 
      !================================================================= 

      ! Initialize
      RC = GC_SUCCESS

#if defined( BPCH_DIAG )

      ! First-time setup
      IF ( FIRST ) THEN

         ! Define species ID flags
         id_POA1       = Ind_( 'POA1'   )
         id_POA2       = Ind_( 'POA2'   )
         id_SOAIE      = Ind_( 'SOAIE'  )
         id_SOAME      = Ind_( 'SOAME'  )
         id_INDIOL     = Ind_( 'INDIOL' )
         id_SOAGX      = Ind_( 'SOAGX'  )
         id_SOAMG      = Ind_( 'SOAMG'  )
         id_LVOCOA     = Ind_( 'LVOCOA' )
         id_ISN1OA     = Ind_( 'ISN1OA' )
         IS_BETANOSAVE = ALLOCATED( BETANOSAVE )

         ! Define logical flags
         IS_POA    = ( id_POA1   > 0 .and. id_POA2   > 0 )
         IS_ISOPOA = ( id_SOAIE  > 0 .and. id_SOAME  > 0 .and. 
     &                 id_INDIOL > 0 .and. id_SOAGX  > 0 .and. 
     &                 id_SOAMG  > 0 .and. id_LVOCOA > 0 .and. 
     &                 id_ISN1OA > 0 )

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! Conversion factor from kg/m3 --> ug/m3
      FACTOR = 1e+9_fp

      ! Convert species units to [kg] for this routine (ewl, 8/12/15)
      CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, 
     &                        State_Chm, 'kg', RC, OrigUnit=OrigUnit )
      IF ( RC /= GC_SUCCESS ) THEN
         MSG = 'Unit conversion error!'
         CALL GC_Error( MSG, RC, 'Start of diag42 in diag42_mod.F')
         RETURN
      ENDIF 

      ! Initialize pointers
      Spc    => State_Chm%Species
      AIRVOL => State_Met%AIRVOL

      ! Loop over grid boxes     
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )  
      DO L = 1, LD42  
      DO J = 1, JJPAR 
      DO I = 1, IIPAR

         !--------------------------------------------------------------
         ! TSOA [ug/m3]
         ! (terpene SOA)
         !--------------------------------------------------------------
         AD42(I,J,L,1) = AD42(I,J,L,1) + TSOA(I,J,L) * FACTOR

         !--------------------------------------------------------------
         ! ISOA [ug/m3]
         ! (isoprene SOA)
         !--------------------------------------------------------------
         ! Semivolatile isoprene OA + mechanistic isoprene OA
         AD42(I,J,L,2) = AD42(I,J,L,2) +
     &                   ( ISOA(I,J,L) + ISOAAQ(I,J,L) ) * FACTOR

         !--------------------------------------------------------------
         ! ASOA [ug/m3]
         ! (benz, tolu, xyle, + NAP/IVOC SOA)
         !--------------------------------------------------------------
         AD42(I,J,L,3) = AD42(I,J,L,3) + ASOA(I,J,L) * FACTOR

         !--------------------------------------------------------------
         ! POA [ug/m3]
         !--------------------------------------------------------------
         IF ( IS_POA ) THEN
            AD42(I,J,L,4) = AD42(I,J,L,4) + OCPO(I,J,L) * FACTOR 
         ELSE
            AD42(I,J,L,4) = AD42(I,J,L,4) +
     &                       ( OCPI(I,J,L) + OCPO(I,J,L) ) * FACTOR 
         ENDIF

         !--------------------------------------------------------------
         ! OPOA [ug/m3]
         !--------------------------------------------------------------
         AD42(I,J,L,5) = AD42(I,J,L,5) + OPOA(I,J,L) * FACTOR

         !-----------------------------------------------------------
         ! Sum of all organic aerosol [ug/m3]
         !-----------------------------------------------------------
         AD42(I,J,L,6) = AD42(I,J,L,6) +
     &                   ( TSOA(I,J,L) + ISOA(I,J,L) + ASOA(I,J,L) +
     &                     OCPO(I,J,L) + OCPI(I,J,L) + OPOA(I,J,L) +
     &                     ISOAAQ(I,J,L) ) * FACTOR

         !-----------------------------------------------------------
         ! Sum of all organic carbon [ugC/m3]
         !-----------------------------------------------------------
         IF ( IS_POA ) THEN
            IF ( OCFOPOA(I,J) > 0 .and. OCFPOA(I,J) > 0 ) THEN
               AD42(I,J,L,7) = AD42(I,J,L,7) +
     &                     ( ( TSOA(I,J,L) + ISOA(I,J,L) + ASOA(I,J,L) +
     &                         OCPI(I,J,L) + OPOA(I,J,L) )/OCFOPOA(I,J)+
     &                         OCPO(I,J,L) / OCFPOA(I,J) ) * FACTOR
            ENDIF
         ELSE
            IF ( OCFOPOA(I,J) > 0 ) THEN
               AD42(I,J,L,7) = AD42(I,J,L,7)  +
     &                     ( ( TSOA(I,J,L) + ISOA(I,J,L) + ASOA(I,J,L) +
     &                         OCPO(I,J,L) + OCPI(I,J,L) + OPOA(I,J,L) )
     &                       / OCFOPOA(I,J) ) * FACTOR
            ENDIF
         ENDIF

         ! Add mechanistic isoprene OA
         IF ( IS_ISOPOA ) THEN

            ! Convert from kg to ugC/m3
            AD42(I,J,L,7) = AD42(I,J,L,7) + (
     &                Spc(I,J,L,id_SOAIE) * ( NCIEPOX * MWCARB /
     &              ( State_Chm%SpcData(id_SOAIE)%Info%MW_g * 1e-3_fp ))
     &              + Spc(I,J,L,id_SOAME) * ( NCIMAE  * MWCARB /
     &              ( State_Chm%SpcData(id_SOAME)%Info%MW_g * 1e-3_fp ))
     &              + Spc(I,J,L,id_INDIOL)* ( NCINDIOL* MWCARB / 
     &              ( State_Chm%SpcData(id_INDIOL)%Info%MW_g* 1e-3_fp ))
     &              + Spc(I,J,L,id_SOAGX) * ( NCGLYX  * MWCARB /
     &              ( State_Chm%SpcData(id_SOAGX)%Info%MW_g * 1e-3_fp ))
     &              + Spc(I,J,L,id_SOAMG) * ( NCMGLY  * MWCARB /
     &              ( State_Chm%SpcData(id_SOAMG)%Info%MW_g * 1e-3_fp ))
     &              + Spc(I,J,L,id_LVOCOA)* ( NCLVOC  * MWCARB /
     &              ( State_Chm%SpcData(id_LVOCOA)%Info%MW_g* 1e-3_fp ))
     &              + Spc(I,J,L,id_ISN1OA)* ( NCISN1  * MWCARB /
     &              ( State_Chm%SpcData(id_ISN1OA)%Info%MW_g* 1e-3_fp ))
     &              ) / AIRVOL(I,J,L) * FACTOR

         ENDIF

         !--------------------------------------------------------------
         ! Sum of biogenic aerosol [ug/m3]
         !--------------------------------------------------------------
         AD42(I,J,L,8) = AD42(I,J,L,8) +
     &                   ( TSOA(I,J,L) + ISOA(I,J,L) + ISOAAQ(I,J,L) ) *
     &                   FACTOR
          
         !--------------------------------------------------------------
         ! NO branching ratio [dimless]
         !--------------------------------------------------------------
         IF ( IS_BETANOSAVE ) THEN
            ! will have zero or junk values if not in troposphere
            AD42(I,J,L,9) = AD42(I,J,L,9) + BETANOSAVE(I,J,L)
         ENDIF

         !--------------------------------------------------------------
         ! POA [ugC/m3]
         !--------------------------------------------------------------
         IF ( IS_POA ) THEN
            IF ( OCFPOA(I,J) > 0 ) THEN
               AD42(I,J,L,10) = AD42(I,J,L,10) +
     &                        ( OCPO(I,J,L) / OCFPOA(I,J) ) * FACTOR
            ENDIF
         ELSE
            IF ( OCFOPOA(I,J) > 0 ) THEN
               AD42(I,J,L,10) = AD42(I,J,L,10) + 
     &                      ( ( OCPI(I,J,L) + OCPO(I,J,L) ) /
     &                          OCFOPOA(I,J) ) * FACTOR
            ENDIF
         ENDIF

         !--------------------------------------------------------------
         ! OPOA [ugC/m3]
         !--------------------------------------------------------------
         IF ( OCFOPOA(I,J) > 0 ) THEN
            AD42(I,J,L,11) = AD42(I,J,L,11) +
     &                       OPOA(I,J,L) / OCFOPOA(I,J) * FACTOR
         ENDIF

         !--------------------------------------------------------------
         ! Additional aerosol tracers
         !--------------------------------------------------------------

         ! OC [ugC/m3]
         IF ( OCFOPOA(I,J) > 0 ) THEN
            AD42(I,J,L,12) = AD42(I,J,L,12) + 
     &                   ( ( OCPI(I,J,L) + OCPO(I,J,L) ) / 
     &                       OCFOPOA(I,J) ) * FACTOR
         ENDIF

         ! BC [ugC/m3]
         AD42(I,J,L,13) = AD42(I,J,L,13) + 
     &                    ( BCPI(I,J,L) + BCPO(I,J,L) ) * FACTOR

         ! SO4 [ug/m3]
         AD42(I,J,L,14) = AD42(I,J,L,14) + SO4(I,J,L) * FACTOR

         ! NH4 [ug/m3]
         AD42(I,J,L,15) = AD42(I,J,L,15) + NH4(I,J,L) * FACTOR

         ! NIT [ug/m3]
         AD42(I,J,L,16) = AD42(I,J,L,16) + NIT(I,J,L) * FACTOR

         ! Sea salt [ug/m3]
         AD42(I,J,L,17) = AD42(I,J,L,17) + 
     &                    ( SALA(I,J,L) + SALC(I,J,L) ) * FACTOR

         !--------------------------------------------------------------
         ! PM2.5 [ug/m3]
         !--------------------------------------------------------------
         AD42(I,J,L,18) = AD42(I,J,L,18) + PM25(I,J,L) * FACTOR

         !--------------------------------------------------------------
         ! Additional diagnostics for SOAGX, SOAMG (tmf, 12/8/07) 
         !
         ! Assume SOAGX mass = GLYX mass, SOAMG mass = MGLY mass
         !--------------------------------------------------------------

         ! SOAGX [ug/m3]
         AD42(I,J,L,19) = AD42(I,J,L,19) + SOAGX(I,J,L) * FACTOR

         ! SOAMG [ug/m3]
         AD42(I,J,L,20) = AD42(I,J,L,20) + SOAMG(I,J,L) * FACTOR

         !--------------------------------------------------------------
         ! Isoprene SOA obtained from irreversible uptake of isoprene
         ! SOA precursors IEPOX and IMAE (eam,  2014)
         !--------------------------------------------------------------
         IF ( IS_ISOPOA ) THEN

            ! IEPOX-OA [ug/m3]:
            AD42(I,J,L,21) = AD42(I,J,L,21) + Spc(I,J,L,id_SOAIE )
     &                       * FACTOR / AIRVOL(I,J,L) 

            ! IMAE-OA [ug/m3]:
            AD42(I,J,L,22) = AD42(I,J,L,22) + Spc(I,J,L,id_SOAME )
     &                       * FACTOR / AIRVOL(I,J,L) 

            ! INDIOL [ug/m3]:
            AD42(I,J,L,23) = AD42(I,J,L,23) + Spc(I,J,L,id_INDIOL)
     &                       * FACTOR / AIRVOL(I,J,L)

            ! LVOCOA [ug/m3]:
            AD42(I,J,L,24) = AD42(I,J,L,24) + Spc(I,J,L,id_LVOCOA)
     &                       * FACTOR / AIRVOL(I,J,L)

            ! ISN1OA [ug/m3]:
            AD42(I,J,L,25) = AD42(I,J,L,25) + Spc(I,J,L,id_ISN1OA)
     &                       * FACTOR / AIRVOL(I,J,L)

         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Convert back to original unit (ewl, 8/12/15)
      CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, 
     &                        State_Chm, OrigUnit, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         MSG = 'Unit conversion error!'
         CALL GC_Error( MSG, RC, 'End of diag42 in diag42_mod.F')
         RETURN
      ENDIF 

      ! Free pointers
      Spc    => NULL()
      AIRVOL => NULL()
#endif

      END SUBROUTINE DIAG42
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: zero_diag42
!
! !DESCRIPTION: Subroutine ZERO\_DIAG42 zeroes all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ZERO_DIAG42
! 
! !REVISION HISTORY: 
!  22 May 2006 - D. Henze, R. Yantosca - Initial version
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! ZERO_DIAG42 begins here!
      !=================================================================
#if defined( BPCH_DIAG )

      ! Exit if ND42 is turned off
      IF ( ND42 == 0 ) RETURN

      ! Zero arrays
      AD42(:,:,:,:) = 0e0
#endif

      END SUBROUTINE ZERO_DIAG42
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_diag42
!
! !DESCRIPTION: Subroutine WRITE\_DIAG42 writes the ND42 diagnostic arrays 
!  to the binary punch file at the proper time.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE_DIAG42( Input_Opt )
!
! !USES:
!
      USE Input_Opt_Mod,      ONLY : OptInput
#if defined( BPCH_DIAG )
      USE BPCH2_MOD,          ONLY : BPCH2
      USE BPCH2_MOD,          ONLY : GET_MODELNAME
      USE BPCH2_MOD,          ONLY : GET_HALFPOLAR
      USE CMN_DIAG_MOD             ! TINDEX
      USE FILE_MOD,           ONLY : IU_BPCH
      USE CMN_SIZE_MOD             ! Size parameters
      USE GC_GRID_MOD,        ONLY : GET_XOFFSET
      USE GC_GRID_MOD,        ONLY : GET_YOFFSET
      USE TIME_MOD,           ONLY : GET_CT_DIAG
      USE TIME_MOD,           ONLY : GET_DIAGb
      USE TIME_MOD,           ONLY : GET_DIAGe
#endif

!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
! 
! !REVISION HISTORY: 
!  22 May 2006 - D. Henze, R. Yantosca - Initial version
!  (1 ) Replace TINY(1d0) with 1d-32 to avoid problems  on SUN 4100 platform
!        (bmy, 9/5/06)
!  (2 ) Use TS_DIAG for scaling instead of TS_DYN. (ccc, 8/18/09)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations (H. Pye)
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( BPCH_DIAG )
      INTEGER           :: CENTER180, HALFPOLAR
      INTEGER           :: L,         M,         N
      INTEGER           :: IFIRST,    JFIRST,    LFIRST        
      REAL*4            :: LONRES,    LATRES
      REAL*4            :: ARRAY(IIPAR,JJPAR,LLPAR)
      !REAL(fp)         :: SCALE(IIPAR,JJPAR)
      REAL(fp)          :: SCALE
      REAL(f8)          :: DIAGb,     DIAGe
      CHARACTER(LEN=20) :: MODELNAME 
      CHARACTER(LEN=40) :: CATEGORY
      CHARACTER(LEN=40) :: RESERVED
      CHARACTER(LEN=40) :: UNIT

      !=================================================================
      ! WRITE_DIAG42 begins here!
      !=================================================================

      ! Exit if ND42 is turned off
      IF ( ND42 == 0 ) RETURN

      ! Initialize
      CENTER180 = 1
      DIAGb     = GET_DIAGb()
      DIAGe     = GET_DIAGe()
      HALFPOLAR = GET_HALFPOLAR()
      IFIRST    = GET_XOFFSET( GLOBAL=.TRUE. ) + 1
      JFIRST    = GET_YOFFSET( GLOBAL=.TRUE. ) + 1
      LATRES    = DJSIZE
      LFIRST    = 1
      LONRES    = DISIZE
      MODELNAME = GET_MODELNAME()
      RESERVED  = ''
      SCALE     = DBLE( GET_CT_DIAG() ) + TINY( 1e0 )

      !=================================================================
      ! Write data to the bpch file
      !=================================================================

      ! Loop over ND42 diagnostic tracers
      DO M = 1, TMAX(42)

         ! Define quantities
         N        = TINDEX(42,M)
         CATEGORY = 'IJ-SOA-$'

         ! Pick proper unit
         SELECT CASE ( N )
            ! SOAupdate: update units (hotp 5/24/10)
            CASE( 7, 10, 11, 12, 13 )
               UNIT = 'ug C/m3'
            CASE( 9 )
               UNIT = 'dimless'
            CASE DEFAULT
               UNIT = 'ug/m3'
         END SELECT

         ! Apply scale factor
         DO L = 1, LD42
            ARRAY(:,:,L) = AD42(:,:,L,N) / SCALE
         ENDDO

         ! Write data to disk
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LD42,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD42) )
      ENDDO
#endif

      END SUBROUTINE WRITE_DIAG42
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_diag42
!
! !DESCRIPTION: Subroutine INIT\_DIAG42 allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_DIAG42( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
#if defined( BPCH_DIAG )

      USE CMN_SIZE_MOD   
      USE ERROR_MOD,          ONLY : ALLOC_ERR
#endif
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
! !REVISION HISTORY: 
!  22 May 2006 - D. Henze, R. Yantosca - Initial version
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  23 Jun 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! INIT_DIAG42 begins here!
      !=================================================================

      ! Initialize
      RC = GC_SUCCESS

#if defined( BPCH_DIAG )
     
      ! This is no longer necessary (bmy, 12/13/17)
      !! Turn off ND42 if SOA tracers are not used
      !IF ( .not. Input_Opt%LSOA ) THEN
      !   ND42 = 0
      !   RETURN
      !ENDIF

      ! Exit if ND42 is turned off
      IF ( ND42 == 0 ) RETURN

      ! Number of levels to save for this diagnostic
      LD42 = MIN( ND42, LLPAR )

      ! 4-D array ("IJ-SOA-$")
      ALLOCATE( AD42( IIPAR, JJPAR, LD42, PD42 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'AD42' )
      AD42 = 0.0_fp

      ! Zero arrays
      CALL ZERO_DIAG42
#endif

      END SUBROUTINE INIT_DIAG42
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_diag42
!
! !DESCRIPTION: Subroutine CLEANUP\_DIAG42 deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_DIAG42
! 
! !REVISION HISTORY: 
!  22 May 2006 - D. Henze, R. Yantosca - Initial version
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_DIAG42 begins here!
      !=================================================================
#if defined( BPCH_DIAG )
      IF ( ALLOCATED( AD42 ) ) DEALLOCATE( AD42 ) 
#endif

      END SUBROUTINE CLEANUP_DIAG42
!EOC
      END MODULE DIAG42_MOD
