! $Id: set_prof.f,v 1.2 2003/07/21 15:09:27 bmy Exp $
      SUBROUTINE SET_PROF( NLON, NLAT, YLAT, MONTH, DAY, 
     &                     P,    T,    SA,   ODCOL, OPTDUST, OPTAER )
!
!******************************************************************************
!  Subroutine SET_PROF sets up atmospheric profiles required by Fast-J using a
!  doubled version of the level scheme used in the CTM.  First pressure and z* 
!  altitude are defined, then O3 and T are taken from the supplied climatology
!  and integrated to the CTM levels (may be overwritten with values directly 
!  from the CTM, if desired) and then black carbon and aerosol profiles are 
!  constructed. (Oliver Wild, 4/7/99, mje, bmy, 7/14/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NLON    (INTEGER) : Grid box longitude index               [unitless]
!  (2 ) NLAT    (INTEGER) : Grid box latitude index                [unitless]
!  (3 ) YLAT    (REAL*8)  : Grid box latitude                      [degrees]
!  (4 ) MONTH   (INTEGER) : Current month number                   [1-12]
!  (5 ) DAY     (INTEGER) : Current day of month                   [1-31]
!  (6 ) P       (REAL*8)  : Surface Pressure - Model top pressure  [hPa]
!  (7 ) T       (REAL*8)  : Vertical temperature profile           [K]
!  (8 ) SA      (REAL*8)  : Surface albedo                         [unitless]
!  (9 ) ODCOL   (REAL*8)  : Vertical optical depth profile         [unitless]
!  (10) OPTDUST (REAL*8)  : Mineral dust opt. depths (1-D profile) [unitless]
!  (11) OPTAER  (REAL*8)  : Aerosol optical depths (1-D profile)   [unitless]
!
!  Important varables passed via "cmn_fj.h" and "jv_cmn.h"
!  ============================================================================
!  (1 ) PJ     :  Pressure at boundaries of model levels [hPa]
!  (2 ) Z      :  Altitude of boundaries of model levels [cm]
!  (3 ) ODCOL  :  Optical depth at each model level
!  (4 ) MASFAC :  Conversion factor for pressure to column density
!  (5 ) TJ     :  Temperature profile on model grid
!  (6 ) DM     :  Air column for each model level [molecules/cm2])
!  (7 ) DO3    :  Ozone column for each model level [molecules/cm2]
!  (8 ) DBC    :  Mass of Black Carbon at each model level [g/cm3]  
!  (9 ) PSTD   :  Approximate pressures of levels for supplied climatology
!
!  NOTES:
!  (1 ) Since we parallelize over columns, T, ODCOL, OPTDUST, and OPTAER
!        are 1-D vectors. In the original code from Oliver Wild, these were 
!        3-D arrays.  Also P and SA are just scalars since we just pass one 
!        surface location at a time w/in the parallel loop. (bmy, 9/13/99)
!  (2 ) Mineral dust profiles are also constructed (rvm, 06/04/00)
!  (3 ) Other aerosol profiles are also constructed (rvm, bmy, 2/27/02)
!  (4 ) Added NLON, NLAT, DAY to the arg list.  Now weight the O3 column by 
!        the observed monthly mean EP-TOMS data.  Also updated comments and 
!        added standard GEOS-CHEM documentation header. (mje, bmy, 7/13/03)
!******************************************************************************
!
      ! References to F90 modules
      USE TOMS_MOD, ONLY : TOMS, DTOMS1, DTOMS2

      IMPLICIT NONE

#     include "cmn_fj.h"  ! IPAR, JPAR, LPAR, JPPJ, JPNL
#     include "jv_cmn.h"  ! NDUST, NAER

      ! Argument
      INTEGER, INTENT(IN)    :: DAY, MONTH,   NLAT, NLON
      REAL*8,  INTENT(IN)    :: P,   T(LPAR), SA,   YLAT
      REAL*8,  INTENT(INOUT) :: ODCOL(LPAR)
      REAL*8,  INTENT(IN)    :: OPTDUST(LPAR,NDUST)  
      REAL*8,  INTENT(IN)    :: OPTAER(LPAR,NAER*NRH)  

      ! Local variables
      INTEGER                :: I, K, L, M, N
      REAL*8                 :: DLOGP,F0,T0,B0,PB,PC,XC,MASFAC,SCALEH
      REAL*8                 :: PSTD(52),OREF2(51),TREF2(51),BREF2(51)
      REAL*8                 :: PROFCOL, DAYTOMS

      !=================================================================
      ! SET_PROF begins here!
      !=================================================================

      ! Now use the hybrid pressure formulation (bmy, 8/22/02)
      DO I = 1, NB
         PJ(I) = ETAA(I) + ( ETAB(I) * P )
      ENDDO
      PJ(NB+1) = 0.d0

      ! Set up cloud and surface properties
      CALL CLDSRF( ODCOL, SA )

      !=================================================================      
      ! Set up pressure levels for O3/T climatology - assume that value
      ! given for each 2 km z* level applies from 1 km below to 1 km 
      ! above, so select pressures at these boundaries. Surface level 
      ! values at 1000 mb are assumed to extend down to the actual 
      ! P(nslon,nslat).
      !=================================================================      
      PSTD(1) = MAX( PJ(1), 1000.D0 )
      PSTD(2) = 1000.D0 * 10.D0**( -1.D0/16.D0 )
      DLOGP   = 10.D0**( -2.D0/16.D0 )
      DO I = 3, 51
         PSTD(I) = PSTD(I-1) * DLOGP
      ENDDO
      PSTD(52) = 0.D0

      ! Mass factor - delta-Pressure [hPa] to delta-Column [molec/cm2]
      MASFAC = 100.D0 * 6.022D+23 / ( 28.97D0 * 9.8D0 * 10.D0 )

      ! Select appropriate monthly and latitudinal profiles
      ! Now use YLAT instead of Oliver's YDGRD(NSLAT) (bmy, 9/13/99) 
      M = MAX( 1, MIN( 12, MONTH                   ) )
      L = MAX( 1, MIN( 18, ( INT(YLAT) + 99 ) / 10 ) )

      ! Temporary arrays for climatology data
      DO I = 1, 51
	 OREF2(I) = OREF(I,L,M)
	 TREF2(I) = TREF(I,L,M)
	 BREF2(I) = BREF(I)
      ENDDO

      ! Apportion O3 and T on supplied climatology z* levels onto CTM levels 
      ! with mass (pressure) weighting, assuming constant mixing ratio and
      ! temperature half a layer on either side of the point supplied.
      DO I = 1, NB
         F0 = 0.D0
         T0 = 0.D0
         B0 = 0.D0
         DO K = 1, 51
            PC = MIN( PJ(I),   PSTD(K)   )
            PB = MAX( PJ(I+1), PSTD(K+1) )
            IF ( PC .GT. PB ) THEN
               XC = ( PC - PB ) / ( PJ(I) - PJ(I+1) )
               F0 = F0 + OREF2(K)*XC
               T0 = T0 + TREF2(K)*XC
               B0 = B0 + BREF2(K)*XC
            ENDIF
         ENDDO
         TJ(I)  = T0
         DO3(I) = F0 * 1.D-6
         DBC(I) = B0
      ENDDO
      
      !=================================================================
      ! Insert model values here to replace or supplement climatology.
      ! Note that CTM temperature is always used in x-section 
      ! calculations (see JRATET); TJ is used in actinic flux 
      ! calculation only.
      !=================================================================    
      !DO I=1,LPAR
      !   DO3(I) = MY_OZONE(i)       ! Volume Mixing Ratio
      !   TJ(I)  = T(I)              ! Kelvin
      !ENDDO
      !DO3(LPAR+1) = MY_OZONE*EXP()  ! Above top of model (or use climatology)
      !TJ(LPAR+1)  = MY_TEMP(LPAR)   ! Above top of model (or use climatology)

      !=================================================================
      ! Calculate effective altitudes using scale height at each level
      !=================================================================
      Z(1) = 0.D0
      DO I = 1, LPAR
         SCALEH = 1.3806D-19 * MASFAC * TJ(I)
         Z(I+1) = Z(I) - ( LOG( PJ(I+1) / PJ(I) ) * SCALEH )
      ENDDO

      !=================================================================
      ! Add Aerosol Column - include aerosol types here. Currently use 
      ! soot water and ice; assume black carbon x-section of 10 m2/g, 
      ! independent of wavelength; assume limiting temperature for 
      ! ice of -40 deg C.
      !=================================================================
      DO I = 1, LPAR
         AER(1,I) = DBC(I) * 10.D0 * ( Z(I+1) - Z(I) )

         ! Turn off uniform black carbon profile (rvm, bmy, 2/27/02)
         AER(1,I) = 0D0

         IF ( T(I) .GT. 233.D0 ) THEN
            AER(2,I) = ODCOL(I)
            AER(3,I) = 0.D0
         ELSE
            AER(2,I) = 0.D0
            AER(3,I) = ODCOL(I)
         ENDIF   

         ! Also add in aerosol optical depth columns (rvm, bmy, 9/30/00)
         DO N = 1, NDUST
            AER(3+N,I) = OPTDUST(I,N)	
         ENDDO
        
         ! Also add in other aerosol optical depth columns (rvm, bmy, 2/27/02)
         DO N = 1, NAER*NRH
            AER(3+N+NDUST,I) = OPTAER(I,N)
         ENDDO

      ENDDO

      DO K = 1, MX
         AER(K,LPAR+1) = 0.D0
      ENDDO

      !=================================================================
      ! Calculate column quantities for FAST-J
      !=================================================================
      PROFCOL = 0d0

      DO I = 1, NB

         ! Monthly mean air Column [molec/cm2]
         DM(I)  = ( PJ(I) - PJ(I+1) ) * MASFAC

         ! Monthly mean O3 column [molec/cm2]
         DO3(I) = DO3(I) * DM(I)

         ! Monthly mean O3 column [DU] 
         PROFCOL = PROFCOL + ( DO3(I) / 2.69d16 )
      ENDDO

      !=================================================================
      ! Now weight the O3 column by the observed monthly mean TOMS.
      ! Missing data is denoted by the flag -999. (mje, bmy, 7/15/03)
      !=================================================================
      DAYTOMS = 0d0

      IF ( DAY <= 15 ) THEN 

         ! Interpolate O3 to current day (w/in first half of month)
         IF ( TOMS(NLON,NLAT)   > -999d0  .AND.
     &        DTOMS1(NLON,NLAT) > -999d0 ) THEN  
            DAYTOMS = TOMS(NLON,NLAT) + DTOMS1(NLON,NLAT) * ( DAY - 15 )
         ENDIF

      ELSE

         ! Interpolate O3 to current day (w/in 2nd half of month)
         IF ( TOMS(NLON,NLAT)   > -999d0  .AND.
     &        DTOMS2(NLON,NLAT) > -999d0 ) THEN  
            DAYTOMS = TOMS(NLON,NLAT) + DTOMS2(NLON,NLAT) * ( DAY - 15 )
         ENDIF

      ENDIF
      
      ! Scale monthly O3 profile to the daily O3 profile (if available)
      IF ( DAYTOMS > 0d0 ) THEN 
         DO I = 1, NB
            DO3(I) = DO3(I) * ( DAYTOMS / PROFCOL )
         ENDDO
      ENDIF
      
!### Debug
!      write (987,100) nlon,nlat,toms(nlon,nlat), profcol, daytoms,
!     $     dtoms1(nlon,nlat), dtoms2(nlon,nlat), SUM( DO3(:) / 2.69d16 )
! 100  format(i7,x,i7,x,6(f8.2,x))

      ! Return to calling program
      END SUBROUTINE SET_PROF
