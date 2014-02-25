!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
! !REMARKS:
!  References:
!  ============================================================================
!  (1)  Parrella, J. P., Jacob, D. J., Liang, Q., Zhang, Y., Mickley, L. J.,
!        Miller, B., Evans, M. J., Yang, X., Pyle, J. A., Theys, N., and Van
!        Roozendael, M.: Tropospheric bromine chemistry: implications for
!        present and pre-industrial ozone and mercury, Atmos. Chem. Phys., 12,
!        6723-6740, doi:10.5194/acp-12-6723-2012, 2012.
!  (2 ) Yang, X., Cox, R. A., Warwick, N. J., Pyle, J. A., Carver, G. D.,
!        O’Connor, F. M., and Savage, N. H.: Tropospheric bromine chemistry and
!        its impacts on ozone: A model study, J. Geophys. Res., 110, D23311,
!        doi:10.1029/2005JD006244, 2005.
!  (2 ) Yang, X., Pyle, J. A., and Cox, R. A.: Sea salt aerosol production and
!        bromine release: Role of snow on sea ice, Geophys. Res. Lett., 35,
!        L16815, doi:10.1029/2008GL034536, 2008.
!
! !REVISION HISTORY:
!  02 Mar 2010 - J. Parrella - Initial version
!  22 May 2012 - M. Payer    - Added ProTeX headers
!  08 Aug 2012 - M. Payer    - Modified for size-dependent depletion factors
!                              from Yang et al. (2008)
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

      INTEGER :: month, IDF
      REAL*8  :: DF
      REAL*8  :: df_size(10)  ! Size-dependent depletion factors
      REAL*8  :: dmid_ref(10) ! Reference diameters for depletion factors
      REAL*8  :: dmid         ! Dry diameter of aerosol [um]
      REAL*8  :: seasonal     ! Seasonal depletion factor

      !=================================================================
      ! EMISS_SSA_BROMINE begins here!
      !=================================================================

      ! Dry diameter of aerosol [um]
      dmid = rmid * 2

      ! only do calculation if we're inside the
      ! range of aerosol sizes observed to be
      ! depeleted in bromide.
      IF ( (dmid < 0.2) .or. (dmid > 10.0) ) THEN
         br2_emiss_kg = 0.d0
         RETURN
      ENDIF

      ! store the month
      month = GET_MONTH()

      ! --------------------------------------------
      ! 1. Calculate Depletion Factor DF, based on:
      !    (a) sea salt diameter (b) month and (c) latitude.
      !
      ! following Yang et al. 2005, 2008
      ! --------------------------------------------

      ! Use size-dependent sea salt bromine depletion factors from
      ! Table 1 of Yang et al. (2008)
      dmid_ref = (/  0.2d0,   0.4d0,  0.8d0,  1.0d0,  1.25d0,
     &               1.5d0,   2.0d0,  4.0d0,  5.0d0,  10.0d0  /) 
      df_size  = (/ -3.82d0, -2.54d0, 0.0d0,  0.23d0, 0.38d0, 
     &               0.37d0,  0.31d0, 0.21d0, 0.16d0, 0.11d0  /)

      ! Sort into diameter bins
      IF (      dmid <= 0.4  ) THEN
         IDF = 1
      ELSE IF ( dmid <= 0.8  ) THEN
         IDF = 2
      ELSE IF ( dmid <= 1.0  ) THEN
         IDF = 3
      ELSE IF ( dmid <= 1.25 ) THEN
         IDF = 4
      ELSE IF ( dmid <= 1.5  ) THEN
         IDF = 5
      ELSE IF ( dmid <= 2.0  ) THEN
         IDF = 6
      ELSE IF ( dmid <= 4.0  ) THEN
         IDF = 7
      ELSE IF ( dmid <= 5.0  ) THEN
         IDF = 8
      ELSE 
         IDF = 9
      ENDIF

      ! Interpolate between sea salt diameters
      DF = df_size(IDF) + ( dmid            - dmid_ref(IDF) ) / 
     &                    ( dmid_ref(IDF+1) - dmid_ref(IDF) ) *
     &                    ( df_size(IDF+1)  - df_size(IDF)  )


      ! Apply seasonality to latitudes south of 30S
      IF ( GET_YMID(1, ilat, 1) < -30.0 ) THEN
         ! Divide by mean value 0.4 = (dfmax+dfmin)/2 to keep
         ! seasonal dependence along with size dependence
         seasonal = ( dfmax + (dfmin - dfmax) / 2.d0 *
     &              ( sin( pi*(month/6.d0 - 0.5) ) + 1 ))
     &              / 0.4
      ELSE
         ! no seasonal dependence for the NH
         seasonal = 1.d0
      ENDIF
      DF = DF * seasonal

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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      USE GRID_MOD,       ONLY : GET_AREA_M2
      USE LOGICAL_MOD,    ONLY : LSSABr2
      USE TRACERID_MOD,   ONLY : IDEBr2
      USE TIME_MOD,       ONLY : GET_TS_EMIS
      USE DIAG_MOD,       ONLY : AD46

      USE CMN_SIZE_MOD         ! Size parameters
      USE COMODE_LOOP_MOD      ! AVG(avagadro's #)
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
!  27 Aug 2012 - M. Payer    - Now parallelize DO loop
!  18 Dec 2012 - M. Payer    - Replace NLONG and NLAT with IIPAR and JJPAR
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
      INTEGER                 :: I, J
      REAL*8                  :: DTEMIS
      REAL*8                  :: AREA_M2

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
         AD46(:,:,20) = 0.d0
         EMISRR(:,:,IDEBr2) = 0.d0
         RETURN
      ENDIF

      ! Loop over Lat and Long boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_M2 )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

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
            AD46(I,J,20) = AD46(I,J,20) + ( EMISRR(I,J,IDEBr2) /
     &                     AREA_M2 )    * ( MWT_BR2 / AVG )
         ENDIF

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      
      END SUBROUTINE EMIT_Br2
!EOC
      END MODULE SSA_BROMINE_MOD
