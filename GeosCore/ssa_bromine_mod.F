!------------------------------------------------------------------------------
!             University of Minnesota Atmospheric Chemistry Group
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: ssa_bromine_mod
!
! !DESCRIPTION: Module SSA\_BROMINE\_MOD contains variables and routines for 
!  emissions of Br2.
!\\
!\\
! !INTERFACE:
!
      MODULE SSA_BROMINE_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: EMISS_SSA_BROMINE
      PUBLIC :: EMIT_Br2
!
! !REVISION HISTORY:
!  02 Mar 2010 - J. Parrella - Initial version
!  22 May 2012 - M. Payer    - Added ProTeX headers 
!EOP
!------------------------------------------------------------------------------
!BOC
!
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_ssa_bromine
!
! !DESCRIPTION: Subroutine EMISS\_SSA\_BROMINE calculates aerosol emissions
!  of Br2.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_SSA_BROMINE( ilat, rmid, p_kgsalt, br2_emiss_kg )
!
! !USES:
!
      USE TIME_MOD, ONLY : GET_MONTH
      USE GRID_MOD, ONLY : GET_YMID
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: ilat         ! Grid latitude index
      REAL*8,  INTENT(IN) :: rmid         ! Dry radius of aerosol
      REAL*8,  INTENT(IN) :: p_kgsalt     ! Seasalt aerosol production [kgNaCl]
!
! !OUTPUT PARAMETERS:
!
      REAL*8, INTENT(OUT) :: br2_emiss_kg ! Br2 emissions [kg NaCl]
!
! !REVISION HISTORY:
!  02 Mar 2010 - J. Parrella - Initial version
!  22 May 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL*4, PARAMETER :: dfmax=0.7
      REAL*4, PARAMETER :: dfmin=0.1
      REAL*4, PARAMETER :: Ra=0.00223     ! Ratio of Br/NaCl [g/g]
      REAL*8, PARAMETER :: pi = 3.14159265358979323846d0
!
! !LOCAL VARIABLES:
!

      INTEGER :: i, j, k, month
      REAL*8  :: DF

      !=================================================================
      ! EMISS_SSA_BROMINE begins here!
      !=================================================================

      ! only do calculation if we're inside the
      ! range of aerosol sizes observed to be
      ! depeleted in bromide.
      IF ( (rmid < 1.0) .or. (rmid > 10.0) ) THEN
         br2_emiss_kg = 0.d0
         RETURN
      ENDIF

      ! store the month
      month = GET_MONTH()

      ! --------------------------------------------
      ! 1. Calculate Depletion Factor DF, based on:
      !    (a) month and (b) latitude.
      !
      ! following Yang et al. 2005
      ! --------------------------------------------
      ! note: this selection should work fine for
      !      our grid sizes. 30S will be a model
      !      box edge (in 4x5 or 2x2.5)... so
      !      the midpoint will work fine as
      !      an indicator. >= is unnecessary... but ok.
      IF ( GET_YMID(1,ilat,1) >= -30.0 ) THEN
         ! no seasonal dependence for this range
         ! that covers the NH.
         DF = 0.5d0

!         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         ! jpp, 4/11/2010: testing the Cape Verde DF values from
!         ! Muller et al. 2010 in the appropriate ranges
!         ! ** JPP, FLAG: TESTING
!         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         IF ( (GET_YMID(ilat) > 10.0d0) .and.
!     &        (GET_YMID(ilat) < 20.0d0) ) THEN
!            DF = 0.8d0
!         ENDIF

      ELSE
         DF = dfmax + (dfmin - dfmax) / 2.d0 *
     &        ( sin( pi*(month/6.d0 - 0.5) ) + 1 )
      ENDIF

      ! --------------------------------------------
      ! Now return the emissions for Br2 given the
      ! Sea-salt mass production.
      ! --------------------------------------------
      ! divide by 2 for stoichiometry of Br- to Br2
      br2_emiss_kg = p_kgsalt * Ra * DF / 2.0d0

      RETURN

      END SUBROUTINE EMISS_SSA_BROMINE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emit_br2
!
! !DESCRIPTION: Subroutine EMIT\_BR2 takes the mass flux of Br2 [kg] emitted
!  from sea-salt and distributes it through the the boundary layer.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMIT_BR2(SSA_Br2)
! 
! !USES:
!
      USE BROMOCARB_MOD,  ONLY : Br_SCALING
      USE COMODE_MOD,     ONLY : JLOP,      REMIS,   VOLUME
      USE GRID_MOD,       ONLY : GET_AREA_M2
      USE LOGICAL_MOD,    ONLY : LVARTROP, LSSABr2
      USE PBL_MIX_MOD,    ONLY : GET_PBL_TOP_L
      USE PRESSURE_MOD,   ONLY : GET_PEDGE
      USE TRACERID_MOD,   ONLY : CTRMB,     IDEMIS,  IDEBr2
      USE TROPOPAUSE_MOD, ONLY : GET_TPAUSE_LEVEL
      USE TIME_MOD,       ONLY : GET_TS_EMIS
      USE LOGICAL_MOD,    ONLY : LNLPBL ! (Lin, 03/31/09) nonlocal PBL
      USE DIAG_MOD,       ONLY : AD46

      USE CMN_SIZE_MOD         ! Size parameters
      USE COMODE_LOOP_MOD      ! IDEMS, NEMIS, AVG(avagadro's #)
      USE CMN_DIAG_MOD         ! Diagnostic integers...
      USE CMN_O3_MOD           ! for EMISRR array
!
! !INPUT PARAMETERS:
!
      REAL*8, INTENT(INOUT)   :: SSA_Br2(IIPAR, JJPAR)
!
! !REVISION HISTORY:
!  02 Mar 2010 - J. Parrella - Initial version
!  22 May 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL*8, PARAMETER       :: mwt_br2 = 0.160d0 !kg/mole
!
! !LOCAL VARIABLES:
!
      INTEGER                 :: I, J,  JLOOP, JLOOP1, LTROP
      INTEGER                 :: L, LL, N, NN,  NBB, NBF, TOP
      REAL*8                  :: DTEMIS, AREA_M2

      ! testing
      REAL*8                  :: total_br2


      !=================================================================
      ! EMIT_BR2 begins here!
      !=================================================================

      ! Emission timestep [s]
      DTEMIS = GET_TS_EMIS() * 60d0

      ! ---------------------------------------------
      ! Debug checking... how does the total mass of
      ! Br2 emissions shape up?
      ! ---------------------------------------------
      IF (LSSABr2) THEN
         ! Note: the divide by 2 is taken care of in above function...
         ! so SSA_Br2 is actually # of Br2 emitted, not Br.
         total_br2 = sum( SSA_Br2(:,:) ) / DTEMIS * 3.1556926d7 * 1.0d-9
      ELSE
         total_br2 = 0.d0
      ENDIF

!jp      write(6, '(a)') '-------------------------------------------'
!jp      write(6, '(a)') 'jpp - total sea-salt Br2 emitted [Tg/yr]:'
!jp      write(6, '(1es12.4)') total_br2
!jp      write(6, '(a)') '-------------------------------------------'
!jp
!jp      print*, 'jpp: beginning EMIT_Br2'
!jp      print*, 'nlat =', nlat, '; nlon =',nlong
!jp      call flush(6)

      ! Now convert the total emission of SSA_Br2
      ! from a total emission over the emission timestep [kg/box] 
      ! to an emission rate [#/box/s].
      SSA_Br2(:,:) = SSA_Br2(:,:) / mwt_br2 / DTEMIS * AVG

      ! jpp, testing sensitivity to sea salt bromine emissions
!      SSA_Br2(:,:) = SSA_Br2(:,:) * 100.d0

      ! -----------------------------------------------
      ! If the sea-salt Br2 emissions logical is
      ! turned off in the input.geos file, then
      ! zero the emissions... REMIS has already been
      ! zero'd for initialization... so just return.
      ! -----------------------------------------------
      IF ( .not. LSSABr2 ) THEN
         AD46(:,:,16) = 0.d0
         EMISRR(:,:,IDEBr2) = 0.d0
         RETURN
      ENDIF

      ! Loop over Lat and Long boxes
      DO J = 1, NLAT
      DO I = 1, NLONG

         ! store the surface area of the box
         AREA_M2 = GET_AREA_M2( I, J, 1 )

         ! store the emission for use inside SMVGEAR
         ! in [molecules/box/s]
         EMISRR(I,J,IDEBr2) = SSA_Br2(I,J)

         ! Scale emissions
         EMISRR(I,J,IDEBr2) = EMISRR(I,J,IDEBr2) * Br_SCALING

         IF ( ND46 > 0 ) THEN
            ! store the emission in the AD46 Biogenic Emissions
            ! diagnostic array [kg/m2/s]
            AD46(I,J,16) = AD46(I,J,16) + ( EMISRR(I,J,IDEBr2) /
     &                     AREA_M2 )    * ( MWT_BR2 / AVG )
         ENDIF

      ENDDO
      ENDDO

      print*, 'jpp: through EMIT_Br2'
      call flush(6)

      RETURN
      
      END SUBROUTINE EMIT_Br2
!EOC
      END MODULE SSA_BROMINE_MOD
