      MODULE SSA_BROMINE_MOD


      IMPLICIT NONE


      PRIVATE

      PUBLIC :: EMISS_SSA_BROMINE
      PUBLIC :: EMIT_Br2



      contains


      SUBROUTINE EMISS_SSA_BROMINE(ilat, rmid, p_kgsalt, br2_emiss_kg)

      ! ----------------------------------------------------------------
      ! This subroutine is meant to calculate the given
      ! aerosol emissions of Br2 provided the following
      ! inputs:
      ! 1. ilat,  integer  :: latitude index for given model box
      ! 2. rmid,  real8    :: the dry radius of aerosol (midpoint of bin)
      ! 3. p_kgsalt, real8 :: production [kg NaCl] for sea-salt aerosol.
      !
      ! output:
      ! br2_emiss_kg       :: Br2 emissions in units of p_kgsalt
      ! ----------------------------------------------------------------

      USE TIME_MOD, ONLY : GET_MONTH
      USE GRID_MOD, ONLY : GET_YMID

      ! -----------------
      ! Input variables
      ! -----------------
      real*8,  intent(in) :: rmid, p_kgsalt
      integer, intent(in) :: ilat ! grid latitude index
      
      ! -----------------
      ! Output variable
      ! -----------------
      real*8, intent(out) :: br2_emiss_kg

      ! -----------------
      ! Local variables
      ! -----------------
      integer :: i, j, k, month
      real*8  :: DF

      ! -----------------
      ! parameters
      ! -----------------
      ! Ra = the ratio of Br/NaCl [g/g]
      real*4, parameter :: dfmax=0.7, dfmin=0.1, Ra=0.00223
      real*8, parameter :: pi = 3.14159265358979323846d0

      ! -----------------
      ! Begin here
      ! -----------------      

      ! only do calculation if we're inside the
      ! range of aerosol sizes observed to be
      ! depeleted in bromide.
      if ( (rmid < 1.0) .or. (rmid > 10.0) ) then
         br2_emiss_kg = 0.d0
         RETURN
      endif


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
      if ( get_ymid(ilat) >= -30.0 ) then
         ! no seasonal dependence for this range
         ! that covers the NH.
         DF = 0.5d0

!         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         ! jpp, 4/11/2010: testing the Cape Verde DF values from
!         ! Muller et al. 2010 in the appropriate ranges
!         ! ** JPP, FLAG: TESTING
!         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         if ( (get_ymid(ilat) > 10.0d0) .and.
!     &        (get_ymid(ilat) < 20.0d0) ) then
!            DF = 0.8d0
!         endif

      else
         DF = dfmax + (dfmin - dfmax) / 2.d0 *
     &        ( sin( pi*(month/6.d0 - 0.5) ) + 1 )
      endif


      ! --------------------------------------------
      ! Now return the emissions for Br2 given the
      ! Sea-salt mass production.
      ! --------------------------------------------
      br2_emiss_kg = p_kgsalt * Ra * DF / 2.0d0 ! divide by 2 for stoichiometry of Br- to Br2


      Return

      End subroutine emiss_ssa_bromine

! ----------------------------------------------------------
      subroutine emit_br2(SSA_Br2)

! **********************************************************
! This subroutine takes the mass flux of Br2 [kg]
! emitted from sea-salt and distributes it through the
! the boundary layer.
! **********************************************************
!

      USE COMODE_MOD,        ONLY : JLOP,      REMIS,   VOLUME
      USE GRID_MOD,          ONLY : GET_AREA_M2
      USE LOGICAL_MOD,       ONLY : LVARTROP, LSSABr2
      USE PBL_MIX_MOD,       ONLY : GET_PBL_TOP_L
      USE PRESSURE_MOD,      ONLY : GET_PEDGE
      USE TRACERID_MOD,      ONLY : CTRMB,     IDEMIS,  IDEBr2
      USE TROPOPAUSE_MOD,    ONLY : GET_TPAUSE_LEVEL
      USE TIME_MOD,          ONLY : GET_TS_EMIS
      USE LOGICAL_MOD,       ONLY : LNLPBL ! (Lin, 03/31/09) nonlocal PBL
      USE DIAG_MOD,          ONLY : AD46

#     include "CMN_SIZE"  ! Size parameters
#     include "comode.h"  ! IDEMS, NEMIS, AVG(avagadro's #)
#     include "CMN_DIAG"  ! Diagnostic integers...
#     include "CMN_O3"    ! for EMISRR array

      ! Input variables
      REAL*8, intent(inout)  :: SSA_Br2(IIPAR, JJPAR)

      ! Local variables
      INTEGER                :: I, J,  JLOOP, JLOOP1, LTROP
      INTEGER                :: L, LL, N, NN,  NBB, NBF, TOP
      REAL*8                 :: DTEMIS, AREA_M2

      ! parameters
      REAL*8, parameter      :: mwt_br2 = 0.160d0 !kg/mole

      ! testing
      real*8 :: total_br2

      ! ----------------------------------------------------
      ! Begin Br2 Emissions Loading
      ! ----------------------------------------------------

      ! Emission timestep [s]
      DTEMIS = GET_TS_EMIS() * 60d0

      ! ---------------------------------------------
      ! Debug checking... how does the total mass of
      ! Br2 emissions shape up?
      ! ---------------------------------------------
      if (LSSABr2) then
         total_br2 = sum( SSA_Br2(:,:) ) / DTEMIS * 3.1556926d7
     &        * 1.0d-9  ! note: the divide by 2 is taken care of in above function... so SSA_Br2 is actually # of Br2 emitted, not Br.
      else
         total_br2 = 0.d0
      endif

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
      SSA_Br2(:,:) = SSA_Br2(:,:) / mwt_br2 / DTEMIS *
     *     AVG

      ! jpp, testing sensitivity to sea salt bromine emissions
!      SSA_Br2(:,:) = SSA_Br2(:,:) * 100.d0

      ! -----------------------------------------------
      ! If the sea-salt Br2 emissions logical is
      ! turned off in the input.geos file, then
      ! zero the emissions... REMIS has already been
      ! zero'd for initialization... so just return.
      ! -----------------------------------------------
      IF ( .not. LSSABr2 ) THEN
         AD46(:,:,9) = 0.d0
         EMISRR(:,:,IDEBr2) = 0.d0
         RETURN
      ENDIF


      ! Loop over Lat and Long boxes
      DO J = 1, NLAT

         ! store the surface area of the box
         AREA_M2 = GET_AREA_M2( J )

         DO I = 1, NLONG

            ! store the emission for use inside SMVGEAR
            ! in [molecules/box/s]
            EMISRR(I,J,IDEBr2) = SSA_Br2(I,J)

            if ( ND46 > 0 ) then
               ! store the emission in the AD46 Biogenic Emissions
               ! diagnostic array [kg/m2/s]
               AD46(I,J,9) = AD46(I,J,9) + ( EMISRR(I,J,IDEBr2) /
     &              AREA_M2 ) * ( MWT_BR2 / AVG )
            endif

         ENDDO

      ENDDO


      print*, 'jpp: through EMIT_Br2'
      call flush(6)

      Return
      
      END SUBROUTINE EMIT_Br2


      ! -----------------------------

      END MODULE SSA_BROMINE_MOD
