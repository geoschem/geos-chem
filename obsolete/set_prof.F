!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_prof
!
! !DESCRIPTION: Subroutine SET\_PROF sets up atmospheric profiles required by 
!  Fast-J using a doubled version of the level scheme used in the CTM.  First 
!  pressure and z* altitude are defined, then O3 and T are taken from the 
!  supplied climatology and integrated to the CTM levels (may be overwritten 
!  with values directly from the CTM, if desired) and then black carbon and 
!  aerosol profiles are constructed.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_PROF( NLON,    NLAT,   YLAT, MONTH,   
     &                     DAY,     T,      SA,   ODCOL, 
     &                     OPTDUST, OPTAER, O3COL )

!
! !USES:
!
      USE CMN_SIZE_MOD,  ONLY : LLPAR
      USE CMN_FJ_MOD,    ONLY : JPMAX, JPPJ  
      USE JV_CMN_MOD            

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)    :: NLON                   ! Grid box lon index
      INTEGER, INTENT(IN)    :: NLAT                   ! Grid box lat index 
      REAL*8,  INTENT(IN)    :: YLAT                   ! Latitude [degrees]
      INTEGER, INTENT(IN)    :: MONTH                  ! Current month
      INTEGER, INTENT(IN)    :: DAY                    ! Day of month
      REAL*8,  INTENT(IN)    :: T(LLPAR)               ! Temperature [K]
      REAL*8,  INTENT(IN)    :: SA                     ! UV albedo [unitless]
      REAL*8,  INTENT(IN)    :: OPTDUST(LLPAR,NDUST)   ! Dust OD [unitless]
      REAL*8,  INTENT(IN)    :: OPTAER(LLPAR,NAER*NRH) ! Aerosol OD [unitless]
      REAL*8,  INTENT(IN)    :: O3COL                  ! Overhd O3 column [DU]
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL*8,  INTENT(INOUT) :: ODCOL(LLPAR)           ! Visible OD [unitless]
!
! !AUTHOR:
!  Oliver Wild & Michael Prather
!
! !REMARKS:
!  References:
!  ============================================================================
!  TOMS/SBUV MERGED TOTAL OZONE DATA, Version 8, Revision 3.
!  Resolution:  5 x 10 deg.
!                                                                             .
!  Source: http://code916.gsfc.nasa.gov/Data_services/merged/index.html
!                                                                             .
!  Contact person for the merged data product:
!  Stacey Hollandsworth Frith (smh@hyperion.gsfc.nasa.gov)
!                                                                             .
!  Important module variables:
!  ============================================================================
!  (1 ) PJ     :  Pressure at boundaries of model levels [hPa]
!  (2 ) Z      :  Altitude of boundaries of model levels [cm]
!  (4 ) MASFAC :  Conversion factor for pressure to column density
!  (5 ) TJ     :  Temperature profile on model grid
!  (6 ) DM     :  Air column for each model level [molecules/cm2])
!  (7 ) DO3    :  Ozone column for each model level [molecules/cm2]
!  (8 ) DBC    :  Mass of Black Carbon at each model level [g/cm3]  
!  (9 ) PSTD   :  Approximate pressures of levels for supplied climatology
!
! !REVISION HISTORY: 
!  01 Jun 1996 - M. Prather & O. Wild - Initial version
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
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  02 Aug 2012 - R. Yantosca - Use online temperature when connecting to GCM
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, K, L, M, N
      REAL*8  :: DLOGP,F0,T0,B0,PB,PC,XC,MASFAC,SCALEH
      REAL*8  :: PSTD(52),OREF2(51),TREF2(51),BREF2(51)
      REAL*8  :: PROFCOL, DAYTOMS

      !=================================================================
      ! SET_PROF begins here!
      !=================================================================

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
#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !%%% For external GCM, take temperature (and maybe eventually ozone)
      DO I=1,LLPAR
      !   DO3(I) = MY_OZONE(i)        ! Volume Mixing Ratio
         TJ(I)  = T(I)                ! Kelvin
      ENDDO
      !DO3(LPAR+1) = MY_OZONE*EXP()   ! Above top of model (or use climatology)
      !TJ(LPAR+1)  = MY_TEMP(LPAR)    ! Above top of model (or use climatology)
#else
      !%%% For traditional G-C we comment all these out
      !DO I=1,LLPAR
      !   DO3(I) = MY_OZONE(i)        ! Volume Mixing Ratio
      !   TJ(I)  = T(I)               ! Kelvin
      !ENDDO
      !DO3(LLPAR+1) = MY_OZONE*EXP()  ! Above top of model (or use climatology)
      !TJ(LLPAR+1)  = MY_TEMP(LLPAR)  ! Above top of model (or use climatology)
#endif

      !=================================================================
      ! Calculate effective altitudes using scale height at each level
      !=================================================================
      Z(1) = 0.D0
      DO I = 1, LLPAR
         SCALEH = 1.3806D-19 * MASFAC * TJ(I)
         Z(I+1) = Z(I) - ( LOG( PJ(I+1) / PJ(I) ) * SCALEH )
      ENDDO

      !=================================================================
      ! Add Aerosol Column - include aerosol types here. Currently use 
      ! soot water and ice; assume black carbon x-section of 10 m2/g, 
      ! independent of wavelength; assume limiting temperature for 
      ! ice of -40 deg C.
      !=================================================================
      DO I = 1, LLPAR
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
         AER(K,LLPAR+1) = 0.D0
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
      ! TOMS/SBUV MERGED TOTAL OZONE DATA, Version 8, Revision 5.
      ! Resolution:  5 x 10 deg.
      !
      ! Methodology (bmy, 3/7/12)
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
      ! As of March 2012, we have TOMS/SBUV data for 1979 thru 2008,
      ! which are stored in the TOMS_200906 data directory.  We will
      ! shortly update to the latest TOMS product.
      !=================================================================

      ! Scale monthly O3 profile to the daily O3 profile (if available)
      IF ( O3COL > 0d0 ) THEN 
         DO I = 1, NB
            DO3(I) = DO3(I) * ( O3COL / PROFCOL )
         ENDDO
      ENDIF

      END SUBROUTINE SET_PROF
!EOC
