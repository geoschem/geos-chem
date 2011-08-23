! $Id: set_prof.f,v 1.1 2009/09/16 14:06:10 bmy Exp $
      SUBROUTINE SET_PROF( NLON, NLAT, YLAT,  MONTH,   DAY, 
     &                     T,    SA,   ODCOL, OPTDUST, OPTAER )
!
!******************************************************************************
!  Subroutine SET_PROF sets up atmospheric profiles required by Fast-J using a
!  doubled version of the level scheme used in the CTM.  First pressure and z* 
!  altitude are defined, then O3 and T are taken from the supplied climatology
!  and integrated to the CTM levels (may be overwritten with values directly 
!  from the CTM, if desired) and then black carbon and aerosol profiles are 
!  constructed. (Oliver Wild, 4/7/99, mje, bmy, 7/14/03, 10/30/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NLON    (INTEGER) : Grid box longitude index               [unitless]
!  (2 ) NLAT    (INTEGER) : Grid box latitude index                [unitless]
!  (3 ) YLAT    (REAL*8)  : Grid box latitude                      [degrees]
!  (4 ) MONTH   (INTEGER) : Current month number                   [1-12]
!  (5 ) DAY     (INTEGER) : Current day of month                   [1-31]
!  (6 ) T       (REAL*8)  : Vertical temperature profile           [K]
!  (7 ) SA      (REAL*8)  : Surface albedo                         [unitless]
!  (8 ) ODCOL   (REAL*8)  : Vertical optical depth profile         [unitless]
!  (9 ) OPTDUST (REAL*8)  : Mineral dust opt. depths (1-D profile) [unitless]
!  (10) OPTAER  (REAL*8)  : Aerosol optical depths (1-D profile)   [unitless]
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
!  References:
!  ============================================================================
!  TOMS/SBUV MERGED TOTAL OZONE DATA, Version 8, Revision 3.
!  Resolution:  5 x 10 deg.
!
!  Source: http://code916.gsfc.nasa.gov/Data_services/merged/index.html
!
!  Contact person for the merged data product:
!  Stacey Hollandsworth Frith (smh@hyperion.gsfc.nasa.gov)
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
!  (5 ) We don't need to initialize the PJ array with ETAA and ETAB anymore.
!        PJ is now defined in "fast_j.f".  Updated comments. (bmy, 10/30/07)
!  (6 ) Modified to use GEOS-5 O3 columns when TOMS/SBUV data don't exist,
!        i.e. after 2008. (ccc, 7/13/09)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,  ONLY : TO3
      USE TOMS_MOD, ONLY : TOMS, DTOMS1, DTOMS2
      USE TIME_MOD, ONLY : GET_YEAR

#     include "define.h"

      USE CMN_FJ_MOD, ONLY : IPAR, JPAR, LPAR, JPMAX, JPPJ  ! IPAR, JPAR, LPAR, JPPJ, JPNL
      USE JV_CMN_MOD            !, ONLY : NDUST, NAER, PJ

      IMPLICIT NONE

      ! Argument
      INTEGER, INTENT(IN)    :: DAY, MONTH,  NLAT, NLON
      REAL*8,  INTENT(IN)    :: SA,   YLAT,  T(LPAR)
      REAL*8,  INTENT(INOUT) :: ODCOL(LPAR)
      REAL*8,  INTENT(IN)    :: OPTDUST(LPAR,NDUST)  
      REAL*8,  INTENT(IN)    :: OPTAER(LPAR,NAER*NRH)  

      ! Local variables
      INTEGER                :: I, K, L, M, N
      INTEGER                :: YEAR
      REAL*8                 :: DLOGP,F0,T0,B0,PB,PC,XC,MASFAC,SCALEH
      REAL*8                 :: PSTD(52),OREF2(51),TREF2(51),BREF2(51)
      REAL*8                 :: PROFCOL, DAYTOMS
      ! TTOMS is to determine if we are using TOMS or GEOS-5 O3 columns.
      ! (ccc, 7/13/09)
      LOGICAL                :: TTOMS = .TRUE.
      LOGICAL,SAVE           :: FIRST = .TRUE.

      !=================================================================
      ! SET_PROF begins here!
      !=================================================================

      ! Get the year of the run to know if using TOMS or GEOS-5 O3 
      ! columns. (ccc, 07/13/09)
      YEAR = GET_YEAR()

#if defined(GEOS_5)
      IF ( YEAR > 2008 ) THEN
         TTOMS = .FALSE.
      ENDIF
#endif

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
      DAYTOMS = 0d0

      IF ( TTOMS ) THEN
         IF ( DAY <= 15 ) THEN 

            ! Interpolate O3 to current day (w/in first half of month)
            IF ( TOMS(NLON,NLAT)   > -999d0  .AND.
     &           DTOMS1(NLON,NLAT) > -999d0 ) THEN  
               DAYTOMS = TOMS(NLON,NLAT) + DTOMS1(NLON,NLAT) * 
     &                   ( DAY - 15 )
            ENDIF

         ELSE

            ! Interpolate O3 to current day (w/in 2nd half of month)
            IF ( TOMS(NLON,NLAT)   > -999d0  .AND.
     &           DTOMS2(NLON,NLAT) > -999d0 ) THEN  
               DAYTOMS = TOMS(NLON,NLAT) + DTOMS2(NLON,NLAT) * 
     &                   ( DAY - 15 )
            ENDIF

         ENDIF

      ELSE

         ! TOMS/SBUV data don't exist so we use GEOS_5 total O3 column, 
         ! TO3. (ccc, 7/13/09)
         IF ( FIRST ) THEN
            WRITE(6,*) 'TOMS/SBUV data are not available for this year'
            WRITE(6,*) 'We use the total O3 column given by GEOS-5.'
            WRITE(6,*) 'The validation by GMAO of these O3 columns is ' 
     &                 // 'in progress.'
            FIRST = .FALSE.
         ENDIF

         DAYTOMS = TO3(NLON, NLAT)

      ENDIF

      ! Scale monthly O3 profile to the daily O3 profile (if available)
      IF ( DAYTOMS > 0d0 ) THEN 
         DO I = 1, NB
            DO3(I) = DO3(I) * ( DAYTOMS / PROFCOL )
         ENDDO
      ENDIF
      
!### Debug
!      write (*,100) nlon,nlat,toms(nlon,nlat), profcol, daytoms,
!     $     dtoms1(nlon,nlat), dtoms2(nlon,nlat), SUM( DO3(:) / 2.69d16 )
! 100  format('DO3 debug ',i7,x,i7,x,6(f8.2,x))

      ! Return to calling program
      END SUBROUTINE SET_PROF
