!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_prof_o3
!
! !DESCRIPTION: Subroutine SET\_PROF\_O3 sets up atmospheric profiles required
!  by RRTMG in the stratosphere  using a doubled version of the level scheme
!  used in the CTM.  First pressure and z* altitude are defined, then O3 and T
!  are taken from the supplied climatology and integrated to the CTM levels
!  (may be overwritten with values directly from the CTM, if desired).
!  This is a stripped down version of SET\_PROF; it does O3 only.
!\\
!\\
! !INTERFACE:
!
SUBROUTINE SET_PROF_O3( YLAT,     MONTH,   DAY,      T_CTM,     &
                        P_CTM,    O3_CTM,  O3_TOMS,  T_CLIM,    &
                        O3_CLIM,  Z_CLIM,  AIR_CLIM, Input_Opt, &
                        State_Grid )

!
! !USES:
!
  USE CMN_FJX_MOD
  USE Input_Opt_Mod,      ONLY : OptInput
  USE PhysConstants            ! Physical constants
  USE PRECISION_MOD            ! For GEOS-Chem Precision (fp)
  USE State_Grid_Mod,     ONLY : GrdState

  IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
  REAL(fp), INTENT(IN)       :: YLAT              ! Latitude (degrees)
  INTEGER,  INTENT(IN)       :: MONTH             ! Month
  INTEGER,  INTENT(IN)       :: DAY               ! Day *of month*
  REAL(fp), INTENT(IN)       :: T_CTM(L1_)        ! CTM temperatures (K)
  REAL(fp), INTENT(IN)       :: O3_TOMS           ! O3 column (DU)
  REAL(fp), INTENT(IN)       :: P_CTM(L1_+1)      ! CTM edge pressures (hPa)
  REAL(fp), INTENT(IN)       :: O3_CTM(L1_)       ! CTM ozone (molec/cm3)
  TYPE(OptInput), INTENT(IN) :: Input_Opt         ! Input options
  TYPE(GrdState), INTENT(IN) :: State_Grid        ! Grid State object
!
! !OUTPUT VARIABLES:
!
  REAL(fp), INTENT(OUT)      :: T_CLIM(L1_)       ! Clim. temperatures (K)
  REAL(fp), INTENT(OUT)      :: Z_CLIM(L1_+1)     ! Edge altitudes (cm)
  REAL(fp), INTENT(OUT)      :: O3_CLIM(L1_)      ! O3 column depth (#/cm2)
  REAL(fp), INTENT(OUT)      :: AIR_CLIM(L1_)     ! O3 column depth (#/cm2)
!
! !AUTHOR:
!  Oliver Wild & Michael Prather
!
! !REMARKS:
!
! !REVISION HISTORY:
!  01 Jun 1996 - M. Prather & O. Wild - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  INTEGER                  :: I, K, L, M, N, LCTM
  REAL(fp)                 :: DLOGP,F0,T0,B0,PB,PC,XC,MASFAC,SCALEH
  REAL(fp)                 :: PSTD(52),OREF2(51),TREF2(51)
  REAL(fp)                 :: PROFCOL
  REAL(fp), PARAMETER      :: ODMAX = 200.0e+0_fp

  ! Local variables for quantities from Input_Opt
  LOGICAL :: USE_ONLINE_O3

  !=================================================================
  ! SET_PROF_O3 begins here!
  !=================================================================

  ! Copy fields from INPUT_OPT
  USE_ONLINE_O3   = Input_Opt%USE_ONLINE_O3

  !=================================================================
  ! Set up pressure levels for O3/T climatology - assume that value
  ! given for each 2 km z* level applies from 1 km below to 1 km
  ! above, so select pressures at these boundaries. Surface level
  ! values at 1000 mb are assumed to extend down to the actual
  ! surface pressure for this lat/lon.
  !=================================================================
  PSTD(1)  = MAX(P_CTM(1),1000.e+0_fp)
  PSTD(2)  = 1000.e+0_fp * 10.e+0_fp ** (-1.e+0_fp/16.e+0_fp)
  DLOGP    = 10.e+0_fp**(-2.e+0_fp/16.e+0_fp)
  DO I=3,51
     PSTD(I) = PSTD(I-1) * DLOGP
  ENDDO
  PSTD(52) = 0.e+0_fp

  ! Mass factor - delta-Pressure [hPa] to delta-Column [molec/cm2]
  MASFAC = 100.e+0_fp * AVO / ( AIRMW * g0 * 10.e+0_fp )

  ! Select appropriate monthly and latitudinal profiles
  ! Now use YLAT instead of Oliver's YDGRD(NSLAT) (bmy, 9/13/99)
  M = MAX( 1, MIN( 12, MONTH                   ) )
  L = MAX( 1, MIN( 18, ( INT(YLAT) + 99 ) / 10 ) )

  ! Temporary arrays for climatology data
  DO I = 1, 51
     OREF2(I) = OREF(I,L,M)
     TREF2(I) = TREF(I,L,M)
  ENDDO

  ! Apportion O3 and T on supplied climatology z* levels onto CTM levels
  ! with mass (pressure) weighting, assuming constant mixing ratio and
  ! temperature half a layer on either side of the point supplied.

  DO I = 1, L1_
     F0 = 0.e+0_fp
     T0 = 0.e+0_fp
     DO K = 1, 51
        PC = MIN( P_CTM(I),   PSTD(K)   )
        PB = MAX( P_CTM(I+1), PSTD(K+1) )
        IF ( PC .GT. PB ) THEN
           XC = ( PC - PB ) / ( P_CTM(I) - P_CTM(I+1) )
           F0 = F0 + OREF2(K)*XC
           T0 = T0 + TREF2(K)*XC
        ENDIF
     ENDDO
     T_CLIM(I)  = T0
     O3_CLIM(I) = F0 * 1.e-6_fp
  ENDDO

  !=================================================================
  ! Calculate effective altitudes using scale height at each level
  !=================================================================
  Z_CLIM(1) = 0.e+0_fp
  DO I = 1, L_
     SCALEH = BOLTZ * 1e+4_fp * MASFAC * T_CLIM(I)
     Z_CLIM(I+1) = Z_CLIM(I) - ( LOG( P_CTM(I+1) / P_CTM(I) ) * SCALEH )
  ENDDO
  Z_CLIM(L1_+1)=Z_CLIM(L1_) + ZZHT

  !=================================================================
  ! Calculate column quantities for for RRTMG (which will only use
  ! tropospheric values)
  !=================================================================
  PROFCOL = 0e+0_fp

  DO I = 1, L1_

     ! Monthly mean air Column [molec/cm2]
     AIR_CLIM(I)  = ( P_CTM(I) - P_CTM(I+1) ) * MASFAC

     ! Monthly mean O3 column [molec/cm2]
     O3_CLIM(I) = O3_CLIM(I) * AIR_CLIM(I)

     ! Monthly mean O3 column [DU]
     PROFCOL = PROFCOL + ( O3_CLIM(I) / 2.69e+16_fp )
  ENDDO

  !! Top values are special (do not exist in CTM data)
  !AIR_CLIM(L1_)     = P_CTM(L1_) * MASFAC
  !O3_CLIM(L1_) = O3_CLIM(L1_) * AIR_CLIM(L1_)

  !=================================================================
  ! Now weight the O3 column by the observed monthly mean TOMS.
  ! Missing data is denoted by the flag -999. (mje, bmy, 7/15/03)
  !
  ! TOMS/SBUV MERGED TOTAL OZONE DATA, Version 8, Revision 3.
  ! Resolution:  5 x 10 deg.
  !
  ! Methodology (bmy, 2/12/07)
  ! ----------------------------------------------------------------
  ! FAST-J comes with its own default O3 column climatology (from
  ! McPeters 1992 & Nagatani 1991), which is stored in the input
  ! file "jv_atms.dat".  These "FAST-J default" O3 columns are used
  ! in the computation of the actinic flux and other optical
  ! quantities for the FAST-J photolysis.
  !
  ! The TOMS/SBUV O3 columns and 1/2-monthly O3 trends (contained
  ! in the TOMS_200701 directory) are read into GEOS-Chem by routine
  ! READ_TOMS in "toms_mod.f".  Missing values (i.e. locations where
  ! there are no data) in the TOMS/SBUV O3 columns are defined by
  ! the flag -999.
  !
  ! After being read from disk in routine READ_TOMS, the TOMS/SBUV
  ! O3 data are then passed to the FAST-J routine "set_prof.f".  In
  ! "set_prof.f", a test is done to make sure that the TOMS/SBUV O3
  ! columns and 1/2-monthly trends do not have any missing values
  ! for (lat,lon) location for the given month.  If so, then the
  ! TOMS/SBUV O3 column data is interpolated to the current day and
  ! is used to weight the "FAST-J default" O3 column.  This
  ! essentially "forces" the "FAST-J default" O3 column values to
  ! better match the observations, as defined by TOMS/SBUV.
  !
  ! If there are no TOMS/SBUV O3 columns (and 1/2-monthly trends)
  ! at a (lat,lon) location for given month, then FAST-J will revert
  ! to its own "default" climatology for that location and month.
  ! Therefore, the TOMS O3 can be thought of as an  "overlay" data
  ! -- it is only used if it exists.
  !
  ! Note that there are no TOMS/SBUV O3 columns at the higher
  ! latitudes.  At these latitudes, the code will revert to using
  ! the "FAST-J default" O3 columns.
  !
  ! As of February 2007, we have TOMS/SBUV data for 1979 thru 2005.
  ! 2006 TOMS/SBUV data is incomplete as of this writing.  For years
  ! 2006 and onward, we use 2005 TOMS O3 columns.
  !
  ! This methodology was originally adopted by Mat Evans.  Symeon
  ! Koumoutsaris was responsible for creating the downloading and
  ! processing the TOMS O3 data files from 1979 thru 2005 in the
  ! TOMS_200701 directory.
  !=================================================================

  ! Updated with UCX
  ! Since we now have stratospheric ozone calculated online, use
  ! this instead of archived profiles for all chemistry-grid cells
  ! The variable O3_CTM is obtained from State_Met%Species, and will be 0
  ! outside the chemgrid (in which case we use climatology)

  ! Scale monthly O3 profile to the daily O3 profile (if available)
  DO I = 1, L1_

     ! Use online O3 values in the chemistry grid if selected
     IF ( (USE_ONLINE_O3)              .and. &
          (I <= State_Grid%MaxChemLev) .and. &
          (O3_CTM(I) > 0e+0_fp) ) THEN

        ! Convert from molec/cm3 to molec/cm2
        O3_CLIM(I) = O3_CTM(I) * (Z_CLIM(I+1)-Z_CLIM(I))

        ! Otherwise, use O3 values from the met fields or TOMS/SBUV
     ELSEIF (O3_TOMS > 0e+0_fp) THEN

        O3_CLIM(I) = O3_CLIM(I) * ( O3_TOMS / PROFCOL )

     ENDIF

  ENDDO

END SUBROUTINE SET_PROF_O3
!EOC
