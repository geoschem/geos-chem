!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: fast_jx_mod
!
! !DESCRIPTION: Module FAST\_JX\_MOD contains routines and variables for 
!  calculating photolysis rates using the Fast-JX scheme (Prather et al).
!  Current implementation is version 6.2.
!\\
!\\
! !INTERFACE: 
!
      MODULE FAST_JX_MOD
!
! !USES:
!
      USE CMN_FJX_MOD

      IMPLICIT NONE
#     include "define.h"

      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: FAST_J
      PUBLIC  :: FJFUNC
      PUBLIC  :: INPHOT
      PUBLIC  :: JRATET
      PUBLIC  :: OPMIE
      PUBLIC  :: PHOTOJ
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: BLKSLV
      PRIVATE :: CLDSRF
      PRIVATE :: FLINT
      PRIVATE :: GAUSSP
      PRIVATE :: GEN
      PRIVATE :: JV_INDEX
      PRIVATE :: JVALUE
      PRIVATE :: LEGND0
      PRIVATE :: MATIN4
      PRIVATE :: MIESCT
      PRIVATE :: MMRAN_16
      PRIVATE :: RD_AOD
      PRIVATE :: RD_JS
      PRIVATE :: RD_XXX
      PRIVATE :: READ_JV_ATMS_DAT
      PRIVATE :: SET_AER
      PRIVATE :: SET_PROF
      PRIVATE :: SPHERE
      PRIVATE :: XSEC1D
      PRIVATE :: XSECO2
      PRIVATE :: XSECO3
!
! !REVISION HISTORY:
!  27 Mar 2013 - S. D. Eastham - Initial version (based on original GEOS-
!  Chem files and instructions from J. Mao)               
!  
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      CONTAINS      
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: blkslv
!
! !DESCRIPTION: Subroutine BLKSLV solves the block tri-diagonal system
!  $A(I)\timesX(I-1) + B(I)\timesX(I) + C(I)\timesX(I+1) = H(I)$.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE BLKSLV
!
! !USES:
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  27 Mar 2013 - S. D. Eastham - Copied from GEOS-Chem v9-01-03
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, K, ID
      REAL*8  :: SUM

      !=================================================================
      ! BLKSLV begins here!
      !=================================================================

!-----------UPPER BOUNDARY ID=1
      CALL GEN(1)
      CALL MATIN4 (B)
      DO I=1,N
         RR(I,1) = 0.0d0
        DO J=1,N
          SUM = 0.0d0
         DO K=1,N
          SUM = SUM - B(I,K)*CC(K,J)
         ENDDO
         DD(I,J,1) = SUM
         RR(I,1) = RR(I,1) + B(I,J)*H(J)
        ENDDO
      ENDDO
!----------CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1
      DO ID=2,ND-1
        CALL GEN(ID)
        DO I=1,N
          DO J=1,N
          B(I,J) = B(I,J) + A(I)*DD(I,J,ID-1)
          ENDDO
          H(I) = H(I) - A(I)*RR(I,ID-1)
        ENDDO
        CALL MATIN4 (B)
        DO I=1,N
          RR(I,ID) = 0.0d0
          DO J=1,N
          RR(I,ID) = RR(I,ID) + B(I,J)*H(J)
          DD(I,J,ID) = - B(I,J)*C1(J)
          ENDDO
        ENDDO
      ENDDO
!---------FINAL DEPTH POINT: ND
      CALL GEN(ND)
      DO I=1,N
        DO J=1,N
          SUM = 0.0d0
          DO K=1,N
          SUM = SUM + AA(I,K)*DD(K,J,ND-1)
          ENDDO
        B(I,J) = B(I,J) + SUM
        H(I) = H(I) - AA(I,J)*RR(J,ND-1)
        ENDDO
      ENDDO
      CALL MATIN4 (B)
      DO I=1,N
        RR(I,ND) = 0.0d0
        DO J=1,N
        RR(I,ND) = RR(I,ND) + B(I,J)*H(J)
        ENDDO
      ENDDO
!-----------BACK SOLUTION
      DO ID=ND-1,1,-1
       DO I=1,N
        DO J=1,N
         RR(I,ID) = RR(I,ID) + DD(I,J,ID)*RR(J,ID+1)
        ENDDO
       ENDDO
      ENDDO
!----------MEAN J & H
      DO ID=1,ND,2
        FJ(ID) = 0.0d0
       DO I=1,N
        FJ(ID) = FJ(ID) + RR(I,ID)*WT(I)
       ENDDO
      ENDDO
      DO ID=2,ND,2
        FJ(ID) = 0.0d0
       DO I=1,N
        FJ(ID) = FJ(ID) + RR(I,ID)*WT(I)*EMU(I)
       ENDDO
      ENDDO
! Output fluxes for testing purposes
c      CALL CH_FLUX
c
      RETURN
      END SUBROUTINE BLKSLV
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cldsrf
!
! !DESCRIPTION: Subroutine CLDSRF sets cloud and surface properties.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLDSRF( ODCOL, SA )
!
! !USES:
!
! !INPUT PARAMETERS:
!
      REAL*8, INTENT(INOUT) :: ODCOL(LLPAR)
      REAL*8, INTENT(IN)    :: SA

! !REMARKS:
! 
! !REVISION HISTORY: 
!  27 Mar 2013 - S. D. Eastham - Copied from GEOS-Chem v9-01-03
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, K
      REAL*8  :: ODSUM, ODMAX, ODTOT

      !=================================================================
      ! CLDSRF begins here!
      !=================================================================

      ! Default lower photolysis boundary as bottom of level 1
      nlbatm = 1
 
      ! Set surface albedo
      RFLECT = dble(SA)
      RFLECT = max(0.d0,min(1.d0,RFLECT))

      ! Zero aerosol column
      do k=1,MX
        do i=1,NB
          AER(k,i) = 0.d0
        enddo
      enddo

      ! Scale optical depths as appropriate - limit column to 'odmax'
      odmax = 200.d0
      odsum =   0.d0
      do i=1,LLPAR
        odcol(i) = dble(odcol(i))
        odsum    = odsum + odcol(i)
      enddo
      if(odsum.gt.odmax) then
        odsum = odmax/odsum
        do i=1,LLPAR
          odcol(i) = odcol(i)*odsum
        enddo
        odsum = odmax
      endif
!
!      ! Use clear-sky conditions
!      do i=1,jpnl
!        odcol(i)=0.d0
!      enddo

      ! Set sub-division switch if appropriate
      odtot=0.d0
      jadsub(nb)=0
      jadsub(nb-1)=0
      do i=nb-1,1,-1
        k=2*i
        jadsub(k)=0
        jadsub(k-1)=0
        odtot=odtot+odcol(i)
        if(odtot.gt.0.d0.and.odcol(i).ne.0.d0.and.
     $                                     dtausub.gt.0.d0) then
          if(odtot.le.dtausub) then
            jadsub(k)=1
            jadsub(k-1)=1
          elseif(odtot.gt.dtausub) then
            jadsub(k)=1
            jadsub(k-1)=0
            do j=1,2*(i-1)
              jadsub(j)=0
            enddo
            go to 20
          endif
        endif
      enddo
 20   continue
c
      return
      END SUBROUTINE CLDSRF
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: fast_j.f
!
! !DESCRIPTION: Subroutine FAST\_J loops over longitude and latitude, and
!  calls PHOTOJ to compute J-Values for each column at every chemistry 
!  time-step.  
!\\
!\\
!  References:
!
! \begin{enumerate}
! \item H. Liu, J.H. Crawford, R.B. Pierce, P. Norris, S.E. Platnick, G. Chen,
!       J.A. Logan, R.M. Yantosca, M.J. Evans, C. Kittaka, Y. Feng, and 
!       X. Tie, \emph{Radiative effect of clouds on tropospheric chemistry 
!       in a global three-dimensional chemical transport model}, 
!       \underline{J. Geophys. Res.}, \textbf{111}, D20303, 
!       doi:10.1029/2005JD006403, 2006. 
!       http://research.nianet.org/~hyl/publications/liu2006\_cloud1.abs.html
! \end{enumerate}
!
! !INTERFACE:
!
#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      SUBROUTINE FAST_J( SUNCOS, OD, ALBD, PEDGE, am_I_Root )  
#else
      SUBROUTINE FAST_J( SUNCOS, OD, ALBD, am_I_Root )  
#endif
!
! !USES:
!
#     include "define.h"

      USE DAO_MOD,      ONLY : T, CLDF
      USE ERROR_MOD,    ONLY : ERROR_STOP, ALLOC_ERR
      USE GRID_MOD,     ONLY : GET_YMID
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TIME_MOD,     ONLY : GET_MONTH, GET_DAY, GET_DAY_OF_YEAR
      USE TIME_MOD,     ONLY : GET_TAU,   GET_YEAR
      USE TOMS_MOD,     ONLY : GET_OVERHEAD_O3

      USE CMN_SIZE_MOD, ONLY : IIPAR,   JJPAR, LLPAR
      USE CMN_SIZE_MOD, ONLY : NDUST,   MAXIJ, NAER,  NRH

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      ! Is this the root CPU?
      LOGICAL, INTENT(IN) :: am_I_Root

      ! Cosine of solar zenith angle [unitless]
      REAL*8,  INTENT(IN) :: SUNCOS(MAXIJ)   

      ! Cloud optical depth [unitless]
      REAL*8,  INTENT(IN) :: OD(LLPAR,IIPAR,JJPAR) 

      ! UV albedo [unitless]
      REAL*8,  INTENT(IN) :: ALBD(IIPAR,JJPAR)    

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      ! Pressure at gridbox interfaces
      REAL*8,  INTENT(IN) :: PEDGE(IIPAR, JJPAR, LLPAR+1)
#endif
!
! !REMARKS:
!
!  Parameter to choose cloud overlap algorithm:
!  ============================================================================
!  (1 ) OVERLAP (INTEGER) : 1 - Linear Approximation (used up to v7-04-12)
!                           2 - Approximate Random Overlap (default)
!                           3 - Maximum Random Overlap (computation intensive)
!
! !REVISION HISTORY:
!  01 Apr 1998 - P. Murti, R. Martin, R. Yantosca - Initial version
!  (1 ) Call this routine EACH chemistry time-step, before solver.
!  (2 ) This routine must know IMAX, JMAX, LMAX. 
!  (3 ) Now use new !$OMP compiler directives for parallelization (bmy, 5/2/00)
!  (4 ) Now reference "cmn_fj.h" and "jv_cmn.h" for the aerosol
!        optical depths (bmy, 10/2/00)
!  (5 ) Add OPTDUST as a local variable -- make OPTDUST private for
!        the parallel DO-loop, since it stores 1 column of aerosol optical
!        depth for each dust type (bmy, rvm, 10/2/00)
!  (6 ) For now, LPAR in "cmn_fj.h" = LGLOB in "CMN_SIZE".  Therefore we 
!        assume that we are always doing global runs. (bmy, 10/2/00)
!  (7 ) Removed obsolete code from 10/2/00 (bmy, 12/21/00)
!  (8 ) Replace {IJL}GLOB w/ IIPAR,JJPAR,LLPAR everywhere.  Also YLMID(NLAT)
!        needs to be referenced by YLMID(NLAT+J0). (bmy, 9/26/01)
!  (9 ) Remove obsolete code from 9/01.  Updated comments. (bmy, 10/24/01)
!  (10) Add OPTAER as a local variable, make it private for the parallel
!        DO loop, since it stores 1 column of aerosol optical depths for each
!        aerosol type.  Pass OPTAER to PHOTOJ via the argument list.  Declare
!        OPTAER as PRIVATE for the parallel DO-loop. (rvm, bmy, 2/27/02)
!  (11) Now reference GET_PEDGE from "pressure_mod.f", which returns the
!        correct "floating" pressure. (dsa, bdf, bmy, 8/20/02)
!  (12) Now reference T from "dao_mod.f" (bmy, 9/23/02)
!  (13) Now uses routine GET_YMID from "grid_mod.f" to compute grid box 
!        latitude.  Now make IDAY, MONTH local variables.  Now use function 
!        GET_DAY_OF_YEAR from "time_mod.f".  Bug fix: now IDAY (as passed to
!        photoj.f) is day of year rather than cumulative days since Jan 1, 
!        1985. (bmy, 2/11/03)
!  (14) Now reference routine GET_YEAR from "time_mod.f".  Added LASTMONTH
!        as a SAVEd variable.  Now call READ_TOMSO3 from "toms_mod.f" at the
!        beginning of a new month (or the first timestep) to read TOMS O3
!        columns which will be used by "set_prof.f".  Now also reference
!        routine GET_DAY from "time_mod.f".  Rename IDAY to DAY_OF_YR. Pass 
!        day of month to PHOTOJ.  Updated comments, cosmetic changes.
!        (bmy, 7/17/03)
!  (15) Bug fix: PRES needs to be the true surface pressure for GEOS-4, but
!        PS-PTOP for all prior GEOS models.  (bmy, 2/6/04)
!  (16) Now account for cloud overlap (Maximum-Random Overlap and Random 
!        Overlap) in each column (hyl, phs, bmy, 9/18/07)
!  (17) Now initialize the PJ array here, instead of two layers below in
!        "set_prof.f".  Now no longer pass PRES to "photoj.f". (bmy, 11/29/07)
!  (18) Now switch to approx. random overlap option (hyl, phs, bmy, 10/7/08)
!  (19) Now can handle GEOS-5 reprocessed met data with OPTDEPTH being
!        in-cloud optical depths. (bmy, hyl, 10/24/08)
!  (10) Remove references to IN_CLOUD_OD (bmy, 10/15/09)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Treat MERRA in the same way as GEOS-5
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!  06 Mar 2012 - R. Yantosca - Now call GET_OVERHEAD_O3 to get the total
!                              overhead O3 column for FAST-J
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  10 Aug 2012 - R. Yantosca - Replace IPAR, JPAR, LPAR w/ IIPAR, JJPAR, LLPAR
!  27 Mar 2013 - S. D. Eastham - Rolled in FAST_JX_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER, SAVE      :: LASTMONTH = -1
      INTEGER            :: NLON, NLAT, DAY,  MONTH, DAY_OF_YR, L
      REAL*8             :: CSZA, PRES, SFCA, YLAT,  O3COL
      REAL*8             :: TEMP(LLPAR), OPTD(LLPAR)
      REAL*8             :: OPTDUST(LLPAR,NDUST)
      REAL*8             :: OPTAER(LLPAR,NAER*NRH)

      ! Local variables for cloud overlap (hyl, phs)
      INTEGER           :: NUMB, KK, I
      INTEGER            :: INDIC(LLPAR+1)
!#if defined( DEVEL )
!      INTEGER, ALLOCATABLE :: INDGEN(:)
!      INTEGER              :: AS
!#else
      INTEGER            :: INDGEN(LLPAR+1)! = (/ (i,i=1,LLPAR+1) /)
!#endif
      INTEGER            :: KBOT(LLPAR)
      INTEGER            :: KTOP(LLPAR)
      INTEGER            :: INDICATOR(LLPAR+2)
      REAL*8             :: FMAX(LLPAR)  ! maximum cloud fraction 
                                            !  in a block, size can be to 
                                            !  FIX(LLPAR)+1
      REAL*8             :: CLDF1D(LLPAR)
      REAL*8             :: ODNEW(LLPAR)

      ! NOTE: Switch from linear approximation (OVERLAP=1) to approximate
      ! random overlap (OVERLAP=2) because we have re-processed the GEOS-5
      ! met data such that OPTDEPTH, TAUCLI, and TAUCLW are now the in-cloud
      ! optical depths rather than the grid-box optical depths. 
      ! (hyl, phs, bmy, 10/7/08)
      INTEGER, PARAMETER :: OVERLAP = 2
      
      LOGICAL, SAVE      :: FIRST = .true.

      !=================================================================
      ! FAST_J begins here!
      !=================================================================

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      IF ( FIRST ) THEN

!         ALLOCATE ( INDGEN(LLPAR+1), STAT=AS)
!         IF ( AS /=0 ) CALL ALLOC_ERR( 'INDGEN' )

         DO i = 1,LLPAR+1
            INDGEN(i) = i       !(/(i,i=1,LLPAR+1)/)
         ENDDO
         FIRST = .FALSE.

      ENDIF
#endif

      ! Get day of year (0-365 or 0-366)
      DAY_OF_YR = GET_DAY_OF_YEAR()

      ! Get current month
      MONTH = GET_MONTH()

      ! Get day of month
      DAY       = GET_DAY()

      !=================================================================
      ! For each (NLON,NLAT) location, call subroutine PHOTOJ (in a 
      ! parallel loop to compute J-values for the entire column.  
      ! J-values will be stored in the common-block variable ZPJ, and 
      ! will be later accessed via function FJFUNC. 
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( NLON,  NLAT,   YLAT,  CSZA,      OPTAER  )
!$OMP+PRIVATE( PRES,  TEMP,   OPTD,  SFCA,      OPTDUST )
!$OMP+PRIVATE( FMAX,  CLDF1D, KK,    NUMB,      L       )
!$OMP+PRIVATE( KBOT,  KTOP,   ODNEW, INDICATOR, INDIC   )
!$OMP+PRIVATE( O3COL                                    )
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over latitudes
      DO NLAT = 1, JJPAR

         ! Loop over longitudes
         DO NLON = 1, IIPAR

            ! Grid box latitude [degrees]
            YLAT         = GET_YMID( NLON, NLAT, 1 )

            ! Cosine of solar zenith angle [unitless] at (NLON,NLAT)
            CSZA         = SUNCOS( (NLAT-1)*IIPAR + NLON ) 

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
            ! Define the PJ array here, using GCM pressures
            DO L = 1, LLPAR+1
               PJ(L)     = PEDGE( NLON, NLAT, L )
            ENDDO
#else
            ! Define the PJ array here, using Ap and Bp from GEOS-Chem
            DO L = 1, NB
               PJ(L)     = GET_PEDGE( NLON, NLAT, L )
            ENDDO
#endif

            ! Top edge of PJ is top of atmosphere (bmy, 2/13/07)
            PJ(NB+1)     = 0d0

            ! Temperature profile [K] at (NLON,NLAT)
            TEMP         = T(NLON,NLAT,1:LLPAR)

            ! Surface albedo [unitless] at (NLON,NLAT)
            SFCA         = ALBD(NLON,NLAT)

            ! Overhead ozone column [DU] at (NLON, NLAT)
            O3COL        = GET_OVERHEAD_O3( NLON, NLAT )

            ! Aerosol OD profile [unitless] at (NLON,NLAT)
            OPTAER(:,:)  = ODAER(NLON,NLAT,:,:)

            ! Mineral dust OD profile [unitless] at (NLON,NLAT)
            OPTDUST(:,:) = ODMDUST(NLON,NLAT,:,:)

            ! Cloud OD profile [unitless] at (NLON,NLAT)
            OPTD         = OD(1:LLPAR,NLON,NLAT)

            !-----------------------------------------------------------
            !### If you want to exclude aerosol OD, mineral dust OD,
            !### or cloud OD, then uncomment the following lines:
            !OPTAER  = 0d0
            !OPTDUST = 0d0
            !OPTD    = 0d0
            !-----------------------------------------------------------

            !===========================================================
            ! CLOUD OVERLAP : LINEAR ASSUMPTION 
            ! Directly use OPTDEPTH = TAUCLD * CLDTOT
            ! 
            ! NOTE: Use this option if you want to compare to results
            !       from GEOS-Chem v7-04-12 and prior versions.
            !===========================================================
            IF ( OVERLAP == 1 ) then

#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_57 )

               ! Column cloud fraction (not less than zero)
               CLDF1D = CLDF(1:LLPAR,NLON,NLAT)
               WHERE ( CLDF1D < 0d0 ) CLDF1D = 0d0
               
               ! NOTE: for the reprocessed GEOS-5 met fields (i.e. with
               ! optical depth & cloud fractions regridded with RegridTau)
               ! OPTD is the in-cloud optical depth.  At this point it has
               ! NOT been multiplied by cloud fraction yet.  Therefore,
               ! we can just apply the linear overlap formula as written 
               ! above (i.e. multiply by cloud fraction). (hyl, bmy, 10/24/08)
               OPTD = OPTD * CLDF1D
#endif

               ! Call FAST-J routines to compute J-values
               CALL PHOTOJ( NLON,  NLAT,     YLAT,    DAY_OF_YR,  
     &                      MONTH, DAY,      CSZA,    TEMP,    
     &                      SFCA,  OPTD,     OPTDUST, OPTAER, 
     &                      O3COL, am_I_Root                     )

            !===========================================================
            ! CLOUD OVERLAP : APPROXIMATE RANDOM OVERLAP
            ! Use OPTDEPTH = TAUCLD * CLDTOT**1.5
            !===========================================================
            ELSE IF ( OVERLAP == 2 ) THEN

               ! Column cloud fraction (not less than zero)
               CLDF1D = CLDF(1:LLPAR,NLON,NLAT)
               WHERE ( CLDF1D < 0d0 ) CLDF1D = 0d0
               
#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_57 )

               ! NOTE: for the reprocessed GEOS-5 met fields (i.e. with
               ! optical depth & cloud fractions regridded with RegridTau)
               ! OPTD is the in-cloud optical depth.  At this point it has
               ! NOT been multiplied by cloud fraction yet.  Therefore,
               ! we can just apply the approximate random overlap formula
               ! as written above (i.e. multiply by cloud fraction^1.5).
               ! (hyl, bmy, 10/24/08)
               OPTD = OPTD * ( CLDF1D )**1.5d0
               
#else
               ! Otherwise, OPTD is the grid-box optical depth and has 
               ! already been multiplied by  the cloud fraction.  Therefore 
               ! we only need to multiply by the square root of the cloud 
               ! fraction here for the approximate random overlap option. 
               ! (hyl, bmy, 10/24/08)
               OPTD = OPTD * SQRT( CLDF1D )

#endif

               ! Call FAST-J routines to compute J-values
               CALL PHOTOJ( NLON,  NLAT,     YLAT,    DAY_OF_YR,  
     &                      MONTH, DAY,      CSZA,    TEMP,  
     &                      SFCA,  OPTD,     OPTDUST, OPTAER,
     &                      O3COL, am_I_Root                    )

            !===========================================================
            ! CLOUD OVERLAP : MAXIMUM RANDOM OVERLAP
            !
            ! The Maximum-Random Overlap (MRAN) scheme assumes that 
            ! clouds in adjacent layers are maximally overlapped to 
            ! form a cloud block and that blocks of clouds separated by 
            ! clear layers are randomly overlapped.  A vertical profile 
            ! of fractional cloudiness is converted into a series of 
            ! column configurations with corresponding fractions 
            ! (see Liu et al., JGR 2006; hyl,3/3/04). 
            !
            ! For more details about cloud overlap assumptions and 
            ! their effect on photolysis frequencies and key oxidants 
            ! in the troposphere, refer to the following articles:
            ! 
            ! (1) Liu, H., et al., Radiative effect of clouds on 
            !      tropospheric chemistry in a global three-dimensional 
            !      chemical transport model, J. Geophys. Res., vol.111, 
            !      D20303, doi:10.1029/2005JD006403, 2006.
            ! (2) Tie, X., et al., Effect of clouds on photolysis and 
            !      oxidants in the troposphere, J. Geophys. Res., 
            !      108(D20), 4642, doi:10.1029/2003JD003659, 2003.
            ! (3) Feng, Y., et al., Effects of cloud overlap in 
            !      photochemical models, J. Geophys. Res., 109, 
            !      D04310, doi:10.1029/2003JD004040, 2004.
            ! (4) Stubenrauch, C.J., et al., Implementation of subgrid 
            !      cloud vertical structure inside a GCM and its effect 
            !      on the radiation budget, J. Clim., 10, 273-287, 1997.
            !-----------------------------------------------------------
            ! MMRAN needs IN-CLOUD optical depth (ODNEW) as input 
            ! Use cloud fraction, instead of OPTD, to form cloud blocks
            ! (hyl,06/19/04)
            !===========================================================
            ELSE IF ( OVERLAP == 3 ) THEN

               ! Initialize
               FMAX(:)   = 0d0  ! max cloud fraction in each cloud block
               ODNEW(:)  = 0d0  ! in-cloud optical depth
               CLDF1D    = CLDF(1:LLPAR,NLON,NLAT)
               INDICATOR = 0

               ! set small negative CLDF or OPTD to zero. 
               ! Set indicator vector.
               WHERE ( CLDF1D <= 0d0 ) 
                  CLDF1D               = 0d0
                  OPTD                 = 0D0
               ELSEWHERE
                  INDICATOR(2:LLPAR+1) = 1
               ENDWHERE

               ! Prevent negative opt depth
               WHERE ( OPTD < 0D0 ) OPTD   = 0D0

               !--------------------------------------------------------
               ! Generate cloud blocks & get their Bottom and Top levels
               !--------------------------------------------------------
               INDICATOR = CSHIFT(INDICATOR, 1) - INDICATOR
               INDIC     = INDICATOR(1:LLPAR+1)

               ! Number of cloud block
               NUMB      = COUNT( INDIC == 1 ) 
               
               ! Bottom layer of each block
               KBOT(1:NUMB) = PACK(INDGEN, (INDIC == 1 ) ) 

               ! Top layer of each block
               KTOP(1:NUMB) = PACK(INDGEN, (INDIC == -1) ) - 1 
             
               !--------------------------------------------------------
               ! For each cloud block, get Max Cloud Fractions, and 
               ! in-cloud optical depth vertical distribution.
               !--------------------------------------------------------
               DO KK = 1, NUMB

                  ! Max cloud fraction
                  FMAX(KK) = MAXVAL( CLDF1D(KBOT(KK):KTOP(KK)) )

#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_57 )

                  ! NOTE: for the reprocessed GEOS-5 met fields (i.e. with
                  ! optical depth & cloud fractions regridded with RegridTau)
                  ! OPTD is the in-cloud optical depth.  At this point it has
                  ! NOT been multiplied by cloud fraction yet.  Therefore,
                  ! we can just set ODNEW = OPTD. (bmy, hyl, 10/24/08)

                  ! ODNEW is adjusted in-cloud OD vertical distrib.
                  ODNEW(KBOT(KK):KTOP(KK)) = OPTD(KBOT(KK):KTOP(KK))

#else

                  ! Otherwise, OPTD is the grid-box optical depth.  
                  ! Therefore, we must divide out by the cloud fraction
                  ! and thus set ODNEW = OPTD / FMAX. (bmy, hyl, 10/24/08)

                  ! ODNEW is adjusted in-cloud OD vertical distrib.
                  ODNEW(KBOT(KK):KTOP(KK)) = OPTD(KBOT(KK):KTOP(KK)) / 
     &                                       FMAX(KK)

#endif
               ENDDO
            
               !--------------------------------------------------------
               ! Apply Max RANdom if 1-6 clouds blocks, else use linear 
               !--------------------------------------------------------
               SELECT CASE( NUMB ) 
           
                  CASE( 0,7: )
                     CALL PHOTOJ( NLON,  NLAT,     YLAT,    DAY_OF_YR, 
     &                            MONTH, DAY,      CSZA,    TEMP,  
     &                            SFCA,  OPTD,     OPTDUST, OPTAER,
     &                            O3COL, am_I_Root                    )

                  CASE( 1:6 ) 
                     CALL MMRAN_16( NUMB,  NLON,  NLAT,      YLAT,   
     &                              DAY,   MONTH, DAY_OF_YR, CSZA,    
     &                              TEMP,  SFCA,  OPTDUST,   OPTAER, 
     &                              LLPAR, FMAX,  ODNEW,     KBOT,   
     &                              KTOP,  O3COL, am_I_Root )

               END SELECT
            ENDIF 
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !-----------------------------------------------------------
      ! END OF SUBROUTINE FAST-J
      !-----------------------------------------------------------
      END SUBROUTINE FAST_J
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fjfunc
!
! !DESCRIPTION: Function FJFUNC supplies J-values to the chemistry solver.  
!\\
!\\
! !INTERFACE:
!
      REAL*8 FUNCTION FJFUNC(I,J,L,NREAC,BRCH,NAME)
!
! !USES:
!
      USE ERROR_MOD, ONLY : ERROR_STOP
!
! !INPUT PARAMETERS: 
!
      INTEGER,           INTENT(IN) :: I,J,L,NREAC,BRCH
      CHARACTER (LEN=4), INTENT(IN) :: NAME
!
! !REVISION HISTORY: 
!  (1  ) "cmn_fj.h" also includes "CMN_SIZE" and "define.h".
!  (2  ) J-values are stored in array "ZPJ" from "cmn_fj.h".
!  (3  ) Now references ERROR_STOP from "error_mod.f".  Updated comments,
!         and made some cosmetic changes. (bmy, 10/15/02)
!  27 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: N
      
      !=================================================================
      ! FJFUNC begins here!
      !
      ! If your compiler has subscript-range checking (-C or 
      ! -check_bounds) then it is recommended to use this option to 
      ! test for the validity of (I,J,L), since repeated IF statements 
      ! are computationally expensive.
      !
      ! If your compiler does not have subscript-range checking, then 
      ! uncomment the following lines to do a manual test for the 
      ! validity of (I,J,L).  
      !=================================================================
      !IF ( I > IPAR .OR. J > JPAR .OR. L > JPNL ) THEN 
      !   STOP 'invalid grid-box # in call to fjfunc - check fjfunc.f'
      !ENDIF

      !=================================================================
      ! RINDEX converts the J-value index as read from "chem.dat" to 
      ! the J-value index as read from "ratj.d". (bmy, 10/5/98)
      !
      ! Make sure that we have taken the proper reaction! 
      !=================================================================
      N = RINDEX(NREAC)

      IF ( N > JPPJ ) THEN
         WRITE(6,*) 'RXN for ',name,', branch ',brch,' not found!'
         CALL ERROR_STOP( 'Check FJFUNC', 'FAST_JX_MOD.F' )
      ENDIF

      !=================================================================
      ! Return the appropriate J-value as the value of the function 
      !=================================================================
      FJFUNC = ZPJ(L,N,I,J)

      ! Return to calling program
      END FUNCTION FJFUNC
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: flint
!
! !DESCRIPTION: Function FLINT performs three-point linear interpolation. 
!\\
!\\
! !INTERFACE:
!
      REAL*8 FUNCTION FLINT (TINT,T1,T2,T3,F1,F2,F3)
!
! !USES:
!
! !INPUT PARAMETERS:
!
      REAL*8, INTENT(IN) :: TINT, T1, T2, T3, F1, F2, F3

! !REVISION HISTORY: 
!  27 Mar 2013 - S. D. Eastham - Copied from GEOS-Chem v9-01-03
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! FLINT begins here!
      !=================================================================

      IF (TINT .LE. T2)  THEN
        IF (TINT .LE. T1)  THEN
          FLINT  = F1
        ELSE
          FLINT = F1 + (F2 - F1)*(TINT -T1)/(T2 -T1)
        ENDIF
      ELSE
        IF (TINT .GE. T3)  THEN
          FLINT  = F3
        ELSE
          FLINT = F2 + (F3 - F2)*(TINT -T2)/(T3 -T2)
        ENDIF
      ENDIF
      return
      END FUNCTION FLINT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gaussp
!
! !DESCRIPTION: Subroutine GAUSSP loads in pre-set Gauss points for 4 angles
!  from 0 to +1 in cos(theta)=mu.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GAUSSP (N,XPT,XWT)
!
! !USES:
!
! !INPUT PARAMETERS:
!
      INTEGER               :: N
      REAL*8                :: XPT(N),XWT(N)

! !REVISION HISTORY: 
!  27 Mar 2013 - S. D. Eastham - Copied from GEOS-Chem v9-01-03
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I
      REAL*8             :: GPT4(4),GWT4(4)

      !=================================================================
      ! GAUSSP begins here!
      !=================================================================
      
      DATA GPT4/.06943184420297D0,.33000947820757D0,.66999052179243D0,
     G          .93056815579703D0/
      DATA GWT4/.17392742256873D0,.32607257743127D0,.32607257743127D0,
     W          .17392742256873D0/
      N = 4
      DO I=1,N
        XPT(I) = GPT4(I)
        XWT(I) = GWT4(I)
      ENDDO
      RETURN
      END SUBROUTINE GAUSSP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gen
!
! !DESCRIPTION: Subroutine GEN generates coefficient matrices for the block 
!  tri-diagonal system described in BLKSLV.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GEN(ID)
!
! !INPUT PARAMETERS:
!
      INTEGER :: ID
!
! !REVISION HISTORY: 
!  27 Mar 2013 - S. D. Eastham - Copied from GEOS-Chem v9-01-03
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: ID0, ID1, IM, I, J, K, MSTART
      REAL*8   :: SUM0, SUM1, SUM2, SUM3
      REAL*8   :: DELTAU, D1, D2, SURFAC
 
      !=================================================================
      ! GEN begins here!
      !=================================================================

!---------------------------------------------
      IF(ID.EQ.1 .OR. ID.EQ.ND) THEN
!---------calculate generic 2nd-order terms for boundaries
       ID0 = ID
       ID1 = ID+1
       IF(ID.GE.ND) ID1 = ID-1
       DO 10 I=1,N
          SUM0 = 0.0d0
          SUM1 = 0.0d0
          SUM2 = 0.0d0
          SUM3 = 0.0d0
        DO IM=M,MFIT,2
          SUM0 = SUM0 + POMEGA(IM,ID0)*PM(I,IM)*PM0(IM)
          SUM2 = SUM2 + POMEGA(IM,ID1)*PM(I,IM)*PM0(IM)
        ENDDO
        DO IM=M+1,MFIT,2
          SUM1 = SUM1 + POMEGA(IM,ID0)*PM(I,IM)*PM0(IM)
          SUM3 = SUM3 + POMEGA(IM,ID1)*PM(I,IM)*PM0(IM)
        ENDDO
         H(I) = 0.5d0*(SUM0*FZ(ID0) + SUM2*FZ(ID1))
         A(I) = 0.5d0*(SUM1*FZ(ID0) + SUM3*FZ(ID1))
        DO J=1,I
          SUM0 = 0.0d0
          SUM1 = 0.0d0
          SUM2 = 0.0d0
          SUM3 = 0.0d0
         DO IM=M,MFIT,2
          SUM0 = SUM0 + POMEGA(IM,ID0)*PM(I,IM)*PM(J,IM)
          SUM2 = SUM2 + POMEGA(IM,ID1)*PM(I,IM)*PM(J,IM)
         ENDDO
         DO IM=M+1,MFIT,2
          SUM1 = SUM1 + POMEGA(IM,ID0)*PM(I,IM)*PM(J,IM)
          SUM3 = SUM3 + POMEGA(IM,ID1)*PM(I,IM)*PM(J,IM)
         ENDDO
         S(I,J) = - SUM2*WT(J)
         S(J,I) = - SUM2*WT(I)
         W(I,J) = - SUM1*WT(J)
         W(J,I) = - SUM1*WT(I)
         U1(I,J) = - SUM3*WT(J)
         U1(J,I) = - SUM3*WT(I)
          SUM0 = 0.5d0*(SUM0 + SUM2)
         B(I,J) = - SUM0*WT(J)
         B(J,I) = - SUM0*WT(I)
        ENDDO
         S(I,I) = S(I,I) + 1.0d0
         W(I,I) = W(I,I) + 1.0d0
         U1(I,I) = U1(I,I) + 1.0d0
         B(I,I) = B(I,I) + 1.0d0
   10  CONTINUE
       DO I=1,N
         SUM0 = 0.0d0
        DO J=1,N
         SUM0 = SUM0 + S(I,J)*A(J)/EMU(J)
        ENDDO
        C1(I) = SUM0
       ENDDO
       DO I=1,N
        DO J=1,N
          SUM0 = 0.0d0
          SUM2 = 0.0d0
         DO K=1,N
          SUM0 = SUM0 + S(J,K)*W(K,I)/EMU(K)
          SUM2 = SUM2 + S(J,K)*U1(K,I)/EMU(K)
         ENDDO
         A(J) = SUM0
         V1(J) = SUM2
        ENDDO
        DO J=1,N
         W(J,I) = A(J)
         U1(J,I) = V1(J)
        ENDDO
       ENDDO
       IF (ID.EQ.1) THEN
!-------------upper boundary, 2nd-order, C-matrix is full (CC)
        DELTAU = ZTAU(2) - ZTAU(1)
        D2 = 0.25d0*DELTAU
        DO I=1,N
          D1 = EMU(I)/DELTAU
          DO J=1,N
           B(I,J) = B(I,J) + D2*W(I,J)
           CC(I,J) = D2*U1(I,J)
          ENDDO
          B(I,I) = B(I,I) + D1
          CC(I,I) = CC(I,I) - D1
!         H(I) = H(I) + 2.0d0*D2*C1(I) + D1*SISOTP
          H(I) = H(I) + 2.0d0*D2*C1(I)
          A(I) = 0.0d0
        ENDDO
       ELSE
!-------------lower boundary, 2nd-order, A-matrix is full (AA)
        DELTAU = ZTAU(ND) - ZTAU(ND-1)
        D2 = 0.25d0*DELTAU
        SURFAC = 4.0d0*ZREFL/(1.0d0 + ZREFL)
        DO I=1,N
          D1 = EMU(I)/DELTAU
          H(I) = H(I) - 2.0d0*D2*C1(I)
           SUM0 = 0.0d0
          DO J=1,N
           SUM0 = SUM0 + W(I,J)
          ENDDO
           SUM0 = D1 + D2*SUM0
           SUM1 = SURFAC*SUM0
          DO J=1,N
           B(I,J) = B(I,J) + D2*W(I,J) - SUM1*EMU(J)*WT(J)
          ENDDO
          B(I,I) = B(I,I) + D1
          H(I) = H(I) + SUM0*ZFLUX
          DO J=1,N
           AA(I,J) = - D2*U1(I,J)
          ENDDO
           AA(I,I) = AA(I,I) + D1
           C1(I) = 0.0d0
        ENDDO
       ENDIF
!------------intermediate points:  can be even or odd, A & C diagonal
      ELSE
        DELTAU = ZTAU(ID+1) - ZTAU(ID-1)
        MSTART = M + MOD(ID+1,2)
        DO I=1,N
          A(I) = EMU(I)/DELTAU
          C1(I) = -A(I)
           SUM0 = 0.0d0
          DO IM=MSTART,MFIT,2
           SUM0 = SUM0 + POMEGA(IM,ID)*PM(I,IM)*PM0(IM)
          ENDDO
          H(I) = SUM0*FZ(ID)
          DO J=1,I
            SUM0 = 0.0d0
           DO IM=MSTART,MFIT,2
            SUM0 = SUM0 + POMEGA(IM,ID)*PM(I,IM)*PM(J,IM)
           ENDDO
            B(I,J) =  - SUM0*WT(J)
            B(J,I) =  - SUM0*WT(I)
          ENDDO
          B(I,I) = B(I,I) + 1.0d0
        ENDDO
      ENDIF
      RETURN

      END SUBROUTINE GEN
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: inphot
!
! !DESCRIPTION: Subroutine INPHOT initializes values for Fast-JX photolysis,
!  including JPL spectral data (e.g. cross sections, quantum yields), standard 
!  O3 and T profiles, and translation indices between GEOS-Chem and Fast-JX
!  species names.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INPHOT( NLAYER, NREACS, am_I_Root )
!
! !USES:
!
#     include "define.h"

      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE CMN_SIZE_MOD, ONLY : LLPAR
      USE inquireMod,   ONLY : findFreeLUN

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: NLAYER      ! # of layers for FAST-JX photolysis
      INTEGER, INTENT(IN) :: NREACS      ! # of FAST-JX photolysis reactions
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
! 
! !REVISION HISTORY: 
!  01 Apr 1999 - O. Wild     - Initial version
!  (1 ) Remove PTOP from the arg list, since it is now a 
!        parameter in "CMN_SIZE" (bmy, 2/10/00).
!  (2 ) Remove SIGE from the argument list, since we are now using
!        a hybrid pressure specification.  Now define ETAA and ETAB
!        for use in "set_prof.f". (bmy, 8/23/02)
!  (3 ) Now reference ERROR_STOP from "error_mod.f".  Updated comments and
!        made cosmetic changes (bmy, 10/15/02)
!  (4 ) Remove IPH -- now use IU_FASTJ directly (bmy, 4/8/03)
!  (5 ) Removed ETAA and ETAB arrays.  We now compute PJ directly from the 
!        GET_PEDGE routine.  Also remove reference to "pressure_mod.f".  
!        Updated comments. (bmy, 10/30/07)
!  (6 ) Read jv_spec_aod.dat file for AOD diagnostics. (clh, bmy, 5/10/10)
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  30 Jul 2012 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now add LUN as a local variable
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!  27 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
      INTEGER            :: LUN
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! INPHOT begins here!
      !=================================================================

      ! # of layers to do chemistry
      JPNL  = NLAYER             

      ! # of reactions in chemistry
      JPPJ  = NREACS + 4         

      ! Error check # of layers
      IF ( JPNL > LLPAR ) THEN 
         CALL ERROR_STOP( 'JPNL > LLPAR!', 'inphot.f' )
      ENDIF

      ! Error check # of rxns
      IF ( JPPJ > JPMAX ) THEN
         CALL ERROR_STOP( 'JPPJ > JPMAX!', 'inphot.f' )
      ENDIF

      !-----------------------------------------------------------------
      ! Read in the names of FAST-J photolysis species from "ratj.d"
      !-----------------------------------------------------------------

      ! Find a free file unit
      LUN      = findFreeLUN()

      ! Name of file to read
      FILENAME = 'ratj.d'
!      FILENAME = 'chemJ2008.d'

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      ! For the ESMF interface to GEOS-5 GCM, we need to append *___.rc
      FILENAME = TRIM( FILENAME ) // '___.rc'
#endif

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( FILENAME ), LUN
 100     FORMAT( 'GEOSCHEM::INPHOT: Reading ', a ', on unit ', i4 )
      ENDIF

      ! Read J_value names 
      CALL RD_JS( LUN, TRIM( FILENAME ), am_I_Root )

      !-----------------------------------------------------------------
      ! Call JV_INDEX to translate between GEOS-Chem species 
      ! nomenclature and Fast-J species nomenclature (bmy, 9/13/99)
      !-----------------------------------------------------------------
      CALL JV_INDEX( am_I_Root )

      !-----------------------------------------------------------------
      ! Read in JPL spectral data set (e.g. X-sections, quantum yields)
      !-----------------------------------------------------------------

      ! Find a free file unit
      LUN      = findFreeLUN()

      ! Name of file to read
      FILENAME = 'jv_spec.dat'
!      FILENAME = 'spec2008.dat'

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      ! For the ESMF interface to GEOS-5 GCM, we need to append *___.rc
      FILENAME = TRIM( FILENAME ) // '___.rc'
#endif

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( FILENAME ), LUN
      ENDIF
      
      ! Read JPL data from disk
      CALL RD_XXX( LUN, TRIM( FILENAME ), am_I_Root )

      !-----------------------------------------------------------------
      ! Read in aerosol optics to be used to calculate AODs (clh)
      !-----------------------------------------------------------------

      ! Find a free file unit
      LUN      = findFreeLUN()

      ! Name of file to read
      FILENAME = 'jv_spec_aod.dat'
#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      ! For the ESMF interface to GEOS-5 GCM, we need to append *___.rc
      FILENAME = TRIM( FILENAME ) // '___.rc'
#endif

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( FILENAME ), LUN
      ENDIF

      ! Read aerosol optics from disk
      CALL RD_AOD( LUN, TRIM( FILENAME ), am_I_Root )

      !-----------------------------------------------------------------
      ! Read in T & O3 climatology (cf. Nagatani/92 and McPeters/91)
      !-----------------------------------------------------------------
      CALL READ_JV_ATMS_DAT( am_I_Root )

      !-----------------------------------------------------------------
      ! Select Aerosol/Cloud types to be used
      !-----------------------------------------------------------------
      CALL SET_AER( am_I_Root )

      END SUBROUTINE INPHOT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: jratet
!
! !DESCRIPTION: Subroutine JRATET calculates and prints J-values. Note that 
!  the loop in this routine only covers the jpnl levels actually needed by 
!  the CTM.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE JRATET( T, IDAY )
!
! !USES:
!
      USE FJX_ACET_MOD

      USE CMN_SIZE_MOD, ONLY : LLPAR

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN) :: T(LLPAR)   ! Temperature [K]
      INTEGER, INTENT(IN) :: IDAY       ! Day of year (0-365 or 0-366)
!
! !REMARKS:
!  FFF    Actinic flux at each level for each wavelength bin
!  QQQ    Cross sections for species (read in in RD_TJPL)
!  SOLF   Solar distance factor, for scaling; normally given by:
!                   1.0-(0.034*cos(real(iday-172)*2.0*pi/365.))
!  TQQ    Temperatures at which QQQ cross sections supplied
!
! !REVISION HISTORY: 
!         1997 - O. Wild     - Initial version
!  (1 ) Added a pressure-dependancy function selector 'pdepf' 
!        in 'jv_spec.dat'. (tmf, 1/7/09)
!  (2 ) Added pressure dependency for MGLY. (tmf, 1/7/09)
!  (3 ) Updated pressure dependency algorithm for ACET. (tmf, 1/7/09)
!  (4 ) Added pressure dependancy for MeCOVi, EtCOMe, MeCOCHO. Rewritten 
!        pressure dependancy for Acetone according to FAST-JX v6.4.
!        See more detailed documentation for Acetone in fjx_acet_mod.f.
!        (ccc, 4/20/09)
!  25 Aug 2011 - R. Yantosca - Rewrite IF statement to prevent PF from
!                              never being initialized.
!  31 Jul 2012 - R. Yantosca - Added ProTeX headers
!  10 Aug 2012 - R. Yantosca - Replace LPAR with LLPAR
!  27 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!     Add Pressure dependancy function selector PF. (tmf, 1/7/09) 
      integer i, j, k, l, PF
      real*8 qptemp

!     For new pressure-dependency algorithm: (tmf, 1/7/09) 
      real*8 xp, xa, xb, xc

!     For new pressure dependency algo. for acetone
!     All variables "*_F" are results from external functions from
!     fjx_acet_mod.f (ccc, 4/20/09)
      real*8 TFACA
      real*8 TFAC0
      real*8 TFAC1, TFAC2
      real*8 QQQA , QQ1A , QQ1B
      real*8 QQ2

      real*8 qo2tot, qo3tot, qo31d, qo33p, qqqt
      !real*8 xseco2, xseco3, xsec1d, solf, tfact
      real*8 solf, tfact

!     Parameters for Solar distance compensation
      real*8  PI, TWOPI
      PARAMETER (PI=3.14159265358979324D0,TWOPI=2.*PI)

!     Physical constants
      REAL*8  Na, R
      PARAMETER (Na=6.02217d23, R=8.3143d0)

!     Scale actinic flux (FFF) by Solar distance factor (SOLF)
      solf=1.d0-(0.034d0*cos(dble(iday-172)*2.d0*pi/365.d0))
!----------------------------------------------------------------------
! If you want to set SOLF = 1.0 for testing, uncomment the next line
!      SOLF = 1d0
!----------------------------------------------------------------------
!
      do I=1,jpnl
       VALJ(1) = 0.d0
       VALJ(2) = 0.d0
       VALJ(3) = 0.d0
       do K=NW1,NW2                       ! Using model 'T's here
         QO2TOT= XSECO2(K,dble(T(I))) 
         VALJ(1) = VALJ(1) + QO2TOT*FFF(K,I)
         QO3TOT= XSECO3(K,dble(T(I)))
         QO31D = XSEC1D(K,dble(T(I)))*QO3TOT
         QO33P = QO3TOT - QO31D
         VALJ(2) = VALJ(2) + QO33P*FFF(K,I)
         VALJ(3) = VALJ(3) + QO31D*FFF(K,I)
       enddo
!------Calculate remaining J-values with T-dep X-sections 
       do J=4,NJVAL
         VALJ(J) = 0.d0
         TFACT = 0.d0
         L = jpdep(J)

!        To choose different forms of pres. dependancy. (ccc, 4/20/09)
         if ( L.ne.0 ) then
            PF = pdepf(L)
         else
            PF = -1
         endif

         if(TQQ(2,J).gt.TQQ(1,J)) TFACT = max(0.d0,min(1.d0,
     $        (T(I)-TQQ(1,J))/(TQQ(2,J)-TQQ(1,J)) ))

!        FAST_JX introduces a new pres. dependancy for acetone (ccc, 4/20/09)
!        Special calculations for the temperature interpolation factors
         if ( PF.eq.2 ) then
            TFACA=TFACA_F(dble(T(I)), J      )
            TFAC0=TFAC0_F(dble(T(I)), J+1    )
            TFAC1=TFAC_F (dble(T(I)), NJVAL+1)
            TFAC2=TFAC_F (dble(T(I)), NJVAL+2)
         else if ( PF.eq.3 ) then
            TFACA=TFACA_F(dble(T(I)), J-1    )
            TFAC0=TFAC0_F(dble(T(I)), J      )
         endif

         do K=NW1,NW2
           QQQT = QQQ(K,1,J-3) + (QQQ(K,2,J-3) - QQQ(K,1,J-3))*TFACT
           if(L.eq.0) then
             VALJ(J) = VALJ(J) + QQQT*FFF(K,I)
           else

              ! Select pressure dependancy function (tmf, 1/31/06)
              if (PF .eq. 1) then
!----------------------------------------------------------------------
! Prior to 9/17/99
! Original form for acetaldehyde P-dep -- believed to be incorrect (pjc)
!             VALJ(J) = VALJ(J) + QQQT*FFF(K,I)*
!     $                   (1.d0+zpdep(K,L)*(pj(i)+pj(i+1))*0.5d0)
!----------------------------------------------------------------------
! Essentially the change is the replacement of the factor
!
!   (1 + a P)     with               1
!                           ---------------------
!                             (1 + b density)
!
! where a and b are constants, P is pressure, and density is the 
! density of air in molec-cm(-3)   (pjc, 9/17/99)
!----------------------------------------------------------------------
              VALJ(J)=VALJ(J)+QQQT*FFF(K,I)/(1 + 
     $                 (zpdep(K,L)*Na*1d-6 /(R*T(I))) * 
     $                 (pj(i)+pj(i+1))*0.5d0*1d2)

             else if ( PF .eq. 4 ) then
!-----------------------------------------------------------------------
! For MGLY
!       y = a + ( b * exp(-p/c) )
!    where y is the ratio between Omega(p) / Omega(p=0);
!          x is the atmospheric pressure [Pa]
!          a,b,c are MGLYPDEP(:,1), MGLYPDEP(:,2), MGLYPDEP(:,3)
!-----------------------------------------------------------------------
                 xp = (pj(i)+pj(i+1))*0.5d0*1.d2   ! pressure [Pa]
                 xa = mglypdep( K, 1 )
                 xb = mglypdep( K, 2 )
                 xc = mglypdep( K, 3 )
                 qptemp = 1.0d0

                 if ( abs( xc ) .ge. 1.d-10 ) then
                    qptemp = xa + ( xb * exp(-xp/xc) )
                 endif

                 VALJ(J) = VALJ(J) + QQQT*FFF(K,I)*qptemp

              else if ( PF.eq.2 ) then
!             Acetone pressure dependency from FAST-JX (ccc, 4/20/09)
!             J1(acetone-a) ==> CH3CO + CH3
!             Special values for Xsect
                 QQQA = QQ1_F (TFACA, J      , K            )
                 QQ2  = QQ2_F (TFAC0, J+1    , K, dble(T(I)))
                 QQ1A = QQ1_F (TFAC1, NJVAL+1, K            )
                 QQ1B = QQ1_F (TFAC2, NJVAL+2, K            ) * 4.d-20

                 VALJ(J) = VALJ(J) + FFF(K,L)*QQQA *
     &            (1.d0-QQ2)/(QQ1A + (QQ1B*Na*1d-6 /(R*T(I))) * 
     $            (pj(i)+pj(i+1))*0.5d0*1d2)
              else if ( PF.eq.3 ) then
!             Second acetone pressure dependency from FAST-JX (ccc, 4/20/09)
!             J2(acetone-b) ==> CH3 + CO + CH3
!             Special values for Xsect
                 QQQA = QQ1_F (TFACA, J-1    , K            )
                 QQ2  = QQ2_F (TFAC0, J      , K, dble(T(I)))

                 VALJ(J) = VALJ(J) + FFF(K,L)*QQQA*QQ2
              endif
           endif
         enddo
       enddo
       do j=1,jppj
         zj(i,j)=VALJ(jind(j))*jfacta(j)*SOLF
       enddo
cc       write(6,'(I5,1P,7E10.3/(5X,7E10.3))') I, (VALJ(J), J=1,NJVAL)
      enddo
      return
      END SUBROUTINE JRATET
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: jv_index
!
! !DESCRIPTION: Subroutine JV\_INDEX computes the mapping between the CTM 
!  indices (from "globchem.dat") for J-values to the FAST-J indices 
!  (from "ratj.d") for J-values.  (bmy, 10/5/98, 4/27/10)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE JV_INDEX( am_I_Root )
!
! !USES:
!
#     include "define.h"

      USE CMN_SIZE_MOD
      USE COMODE_LOOP_MOD  ! SMVGEAR II arrays

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
! 
! !REVISION HISTORY:
!  05 Oct 1998 - R. Yantosca - Initial version
!  (1 ) Assumes the ordering of a species with several branches in 
!        "ratj.d" is the same as in "chem.dat".
!  (2 ) Updated comments, cosmetic changes (bmy, 11/15/01)
!  (3 ) NAMESPEC is now NAMEGAS for SMVGEAR II.   We don't need to reference 
!        CMN anymore. Now loop from NCS = 1..NCSGAS (bdf, bmy, 4/8/03)
!  (4 ) Now reset NCS to NCSURBAN after loop (dbm, bmy, 10/16/06)
!  (5 ) Increase species name length (fp, 6/09)
!  (6 ) Change "Harvard #" to "GEOS-Chem #" in output (bmy, 4/27/10)
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  30 Jul 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local variables
      INTEGER            :: I, IFNC, IBRCH, N, NK

      ! increase length of name from 4 to 7 (FP 6/2009)
      !CHARACTER (LEN=4)  :: SPECNAME
      CHARACTER (LEN=7)  :: SPECNAME

      !=================================================================
      ! JV_INDEX begins here!
      !=================================================================

      ! Zero the RINDEX array
      RINDEX(:) = 0

      ! Loop over photolysis rxns (urban chemistry only)
      DO NCS = 1, NCSGAS
      DO I   = 1, NPHOT

         !==============================================================
         ! I        = Index of photo rxns    from "globchem.dat"
         ! NK       = Absolute rxn number (adds offset to I)
         ! SPECNAME = Name of species I,     from "globchem.dat"
         ! IBRCH    = Branch # of species I, from "globchem.dat"
         !==============================================================
         NK       = NRATES(NCS) + I
         SPECNAME = NAMEGAS(IRM(1,NK,NCS)) 
         IFNC     = DEFPRAT(NK,NCS) + 0.01d0
         IBRCH    = 10d0*( DEFPRAT(NK,NCS) - IFNC ) + 0.5d0

         !==============================================================
         ! N      = Index of photolysis reactions as listed in "ratj.d"
         ! RNAMES = Name of species N,            as listed in "ratj.d" 
         ! BRANCH = Branch number of species N,   as listed in "ratj.d" 
         !  
         ! If the species names and branch numbers from both "chem.dat" 
         ! and "ratj.d" match, then store N (the "ratj.d" index) in the 
         ! Ith element of RINDEX.
         !  
         ! Thus, when looping over I (the chem.dat" indices), as is 
         ! done in FJFUNC.F, RINDEX(I) will access the correct J-value 
         ! according to the ordering in "ratj.d".
         !==============================================================
         DO N = 1, JPPJ
            IF ( SPECNAME == RNAMES(N) .and. IBRCH == BRANCH(N) ) THEN
               RINDEX(I) = N

               ! change formatting for longer name (FP 6/2009)
               IF ( am_I_Root ) THEN
                  WRITE ( 6, 100 ) I,         SPECNAME,  IBRCH, 
     &                             RINDEX(I), RNAMES(N), BRANCH(N)
               ENDIF
 100           FORMAT('GEOS-Chem #: ', i3, 1x, a7, ' Branch: ', i2, 
     &                ' --->  Fast-J #: ', i3, 1x, a7, ' Branch: ',i2 )
               EXIT
            ENDIF
         ENDDO
      ENDDO  
      ENDDO

      ! Reset NCS to NCSURBAN for safety's sake (bmy, 10/16/06)
      NCS = NCSURBAN

      ! Return to calling program      
      END SUBROUTINE JV_INDEX
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: jvalue 
!
! !DESCRIPTION: subroutine JVALUE calculates the actinic flux at each level 
!  for the current solar zenith angle.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE JVALUE( SA, am_I_Root )
!
! !USES:
!
#     include "define.h"

      USE CMN_SIZE_MOD, ONLY : LLPAR

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN) :: SA          ! Surface albedo [unitless]
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! !REMARKS:
!-----------------------------------------------------------------------
!        quit when SZA > 98.0 deg ==> tangent height = 63 km
!             or         99.                           80 km
!-----------------------------------------------------------------------
!  AVGF   Attenuation of beam at each level for each wavelength
!  FFF    Actinic flux at each desired level
!  WAVE   Effective wavelength of each wavelength bin
!  XQO2   Absorption cross-section of O2
!  XQO3   Absorption cross-section of O3
! 
! !REVISION HISTORY: 
!         1997 - O. Wild     - Initial version
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  31 Jul 2012 - R. Yantosca - Added ProTeX headers
!  10 Aug 2012 - R. Yantosca - Replace LPAR with LLPAR
!  27 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer j, k
!      real*8  wave, xseco3, xseco2
      real*8  wave
      real*8  AVGF(LLPAR),XQO3(NB),XQO2(NB)
!
      do J=1,jpnl
        do K=NW1,NW2
          FFF(K,J) = 0.d0
        enddo
      enddo
c
c---SZA check
!      IF(  am_I_Root ) THEN
!         write(6,1000) SZA, RFLECT, (OD(nslon,nslat,j),j=1,LLPAR)
!      ENDIF
      if(SZA.gt.szamax) GOTO 99
c
!---Calculate spherical weighting functions
      CALL SPHERE
c
!---Loop over all wavelength bins
      do K=NW1,NW2
        WAVE = WL(K)
        do J=1,NB
          XQO3(J) = XSECO3(K,dble(TJ(J)))
        enddo
        do J=1,NB
          XQO2(J) = XSECO2(K,dble(TJ(J)))
        enddo
!-----------------------------------------
        CALL OPMIE(K,WAVE,XQO2,XQO3,AVGF)
!-----------------------------------------
        do J=1,jpnl
          FFF(K,J) = FFF(K,J) + FL(K)*AVGF(J)
        enddo
      enddo
c
   99 continue
 1000 format('  SZA=',f6.1,' Reflectvty=',f6.3,' OD=',10(1pe10.3))
      return
      END SUBROUTINE JVALUE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: legnd0
!
! !DESCRIPTION: Subroutine LEGND0 calculates ordinary legendre functions 
!  of X (real).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LEGND0 (X,PL,N)
!
! !INPUT PARAMETERS:
!
      INTEGER             :: N
      REAL*8              :: X, PL(N)
! 
! !REVISION HISTORY: 
!  27 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I
      REAL*8             :: DEN

      !=================================================================
      ! LEGND0 begins here!
      !=================================================================

!---Always does PL(2) = P[1]
        PL(1) = 1.D0
        PL(2) = X
        DO I=3,N
         DEN = (I-1)
         PL(I) = PL(I-1)*X*(2.d0-1.D0/DEN) - PL(I-2)*(1.d0-1.D0/DEN)
        ENDDO
      RETURN

      END SUBROUTINE LEGND0
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: mmran_16
!
! !DESCRIPTION: Subroutine MMRAN_16 does the maximum random cloud overlap
!  for 1 to 6 cloud blocks at a time, and calls PHOTOJ to compute J-values
!  for one column. (hyl, phs, bmy, 9/18/07, 11/29/07)
!
!  Arguments as Input: 
!  ============================================================================
!  Variable  Type    Dimension  Units   Description
!  --------  ----    ---------  -----   -----------
!  Those for PHOTOJ:
!
!  NLON      INT        -         -     Longitude index
!  NLAT      INT        -         -     Latitude index
!  YLAT      DBLE       -         -     Latitude
!  MONTH     INT        -         -     Month of year (1-12)
!  DAY       INT        -         -     Day of the month
!  DAY_OF_YR INT        -         -     Day of the year
!  CSZA      DBLE       -         -     Cosine of solar zenith angle 
!                                        at nlon, nlat
!  PRES      DBLE       -        [mb]   Column pressure at nlon, nlat
!  TEMP      DBLE    [LMAX]      [K]    Layer temperatures at nlon, nlat 
!  SFCA      DBLE       -         -     Surface albedo at nlon, nlat  
!  OPTDUST   DBLE    [LMAX,NDUST] -     Dust optical depths 
!                                        (for NDUST dust types)
!  OPTAER    DBLE [LMAX,NAER*NRH]  -     Aerosol optical depths
!                                          (for NAER aerosol types)
!
!  and those specifically for MMRAN:
!
!  NCB      INT         -         -    Number of cloud blocks
!  MAXBLK   INT         -         -    Dimension of FMAX, 
!  FMAX     DBLE    [MAXBLK]      -    Largest cloud fraction in block
!  ODNEW    DBLE     [LPAR]       -    In-cloud optical depth
!  KBOT     INT      [LPAR]       -    Index of bottom layer of each block
!  KTOP     INT      [LPAR]       -    Index of top layer of each block
!
! LOCAL VARIABLE:
!  OPTD     DBLE    [LPAR]        -    Layer optical depths at nlon, nlat
!  JSUM     DBLE  [LPAR,JPMAX]    -    accumulate the J-values for the column
!  ============================================================================
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE MMRAN_16( NCB,     NLON,      NLAT,   YLAT, DAY,  
     &                     MONTH,   DAY_OF_YR, CSZA,   TEMP, SFCA, 
     &          SAFE_OPTDUST, SAFE_OPTAER,    MAXBLK, FMAX, ODNEW,   
     &                     KBOT,    KTOP,      O3COL,  am_I_Root )
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : NDUST, NAER, NRH, LLPAR
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)    :: NCB ! Number of Cloud Blocks
      INTEGER, INTENT(IN)    :: NLON, NLAT
      REAL*8,  INTENT(IN)    :: CSZA, SFCA, YLAT
      INTEGER, INTENT(IN)    :: DAY, MONTH, DAY_OF_YR
      REAL*8,  INTENT(IN)    :: TEMP(LLPAR)
      REAL*8,  INTENT(IN)    :: SAFE_OPTDUST(LLPAR,NDUST)
      REAL*8,  INTENT(IN)    :: SAFE_OPTAER(LLPAR,NAER*NRH)
!      REAL*8,  INTENT(INOUT) :: OPTDUST(LLPAR,NDUST)
!      REAL*8,  INTENT(INOUT) :: OPTAER(LLPAR,NAER*NRH)
      INTEGER, INTENT(IN)    :: MAXBLK
      REAL*8,  INTENT(IN)    :: FMAX(MAXBLK)
      REAL*8,  INTENT(IN)    :: ODNEW(LLPAR)
      INTEGER, INTENT(IN)    :: KBOT(LLPAR)
      INTEGER, INTENT(IN)    :: KTOP(LLPAR)
      REAL*8,  INTENT(IN)    :: O3COL
      LOGICAL, INTENT(IN)    :: am_I_Root
! 
! !REVISION HISTORY: 
!  (1 ) Remove PRES as an argument, since we no longer need to pass that
!        to PHOTOJ. (bmy, 11/29/07)
!  27 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER    :: II,  JJ,  KK,  LL,  MM,  NN
      INTEGER    :: II2, JJ2,      LL2, MM2, NN2
      REAL*8     :: P1, P2, P3, P4, P5, P6
      REAL*8     :: JSUM(LLPAR,JPMAX)
      REAL*8     :: OPTD(LLPAR)
      REAL*8     :: OPTDUST(LLPAR,NDUST)
      REAL*8     :: OPTAER(LLPAR,NAER*NRH)

      !=================================================================
      ! MMRAN_16 begins here!
      !=================================================================

      ! Copy inputs across
      OPTDUST = SAFE_OPTDUST
      OPTAER  = SAFE_OPTAER

      ! Initialize J-value array
      JSUM = 0d0
 
      ! Initialize Pi
      P1=1d0
      P2=1d0
      P3=1d0 
      P4=1d0
      P5=1d0
      P6=1d0
 
      ! Define the number of loops
      II2 = 1
      JJ2 = 1
      LL2 = 1
      MM2 = 1
      NN2 = 1     
      
      IF ( NCB > 1 ) LL2 = 2         ! At least 2 block-clouds
      IF ( NCB > 2 ) MM2 = 2         ! At least 3 block-clouds
      IF ( NCB > 3 ) NN2 = 2         ! At least 4 block-clouds
      IF ( NCB > 4 ) II2 = 2         ! At least 5 block-clouds
      IF ( NCB > 5 ) JJ2 = 2         ! At least 6 block-clouds


      ! Loop over cloud blocks
      DO KK = 1, 2
      DO LL = 1, LL2
      DO MM = 1, MM2
      DO NN = 1, NN2
      DO II = 1, II2
      DO JJ = 1, JJ2

         ! Zero optical depth
         OPTD(:) = 0d0

         ! 1st cloud block
         IF ( KK == 1 ) THEN
            OPTD(KBOT(1):KTOP(1)) = 0d0
            P1                    = 1d0 - FMAX(1)
         ELSE
            OPTD(KBOT(1):KTOP(1)) = ODNEW(KBOT(1):KTOP(1))
            P1                    = FMAX(1)
         ENDIF


         ! 2nd cloud block
         IF ( NCB > 1 ) THEN
         IF ( LL == 1 ) THEN
            OPTD(KBOT(2):KTOP(2)) = 0d0
            P2                    = 1d0 - FMAX(2)
         ELSE
            OPTD(KBOT(2):KTOP(2)) = ODNEW(KBOT(2):KTOP(2))
            P2                    = FMAX(2)
         ENDIF


         ! 3rd cloud block
         IF ( NCB > 2 ) THEN
         IF ( MM == 1 ) THEN
            OPTD(KBOT(3):KTOP(3)) = 0d0
            P3                    = 1d0 - FMAX(3)
         ELSE
            OPTD(KBOT(3):KTOP(3)) = ODNEW(KBOT(3):KTOP(3))
            P3                    = FMAX(3)
         ENDIF


         ! 4th cloud block
         IF ( NCB > 3 ) THEN         
         IF ( NN == 1 ) THEN
            OPTD(KBOT(4):KTOP(4)) = 0d0
            P4                    = 1d0 - FMAX(4)
         ELSE
            OPTD(KBOT(4):KTOP(4)) = ODNEW(KBOT(4):KTOP(4))
            P4                    = FMAX(4)
         ENDIF


         ! 5th cloud block
         IF ( NCB > 4 ) THEN
         IF ( II == 1 ) THEN
            OPTD(KBOT(5):KTOP(5)) = 0d0
            P5                    = 1d0 - FMAX(5)
         ELSE
            OPTD(KBOT(5):KTOP(5)) = ODNEW(KBOT(5):KTOP(5))
            P5                    = FMAX(5)
         ENDIF


         ! 6th cloud block
         IF ( NCB > 5 ) THEN
         IF ( JJ == 1 ) THEN
            OPTD(KBOT(6):KTOP(6)) = 0d0
            P6                    = 1d0 - FMAX(6)
         ELSE
            OPTD(KBOT(6):KTOP(6)) = ODNEW(KBOT(6):KTOP(6))
            P6                    = FMAX(6)
         ENDIF
         
         ENDIF
         ENDIF
         ENDIF
         ENDIF
         ENDIF

         ! Call the photolysis routine with the OPTD as
         ! computed from the cloud overlaps
         CALL PHOTOJ( NLON, NLAT, YLAT, DAY_OF_YR, MONTH,   DAY,
     &                CSZA, TEMP, SFCA, OPTD,      OPTDUST, OPTAER,
     &                O3COL, am_I_Root )

         ! Store the J values into JSUM array
         JSUM(:,:) = JSUM(:,:) +
     &        ( P1 * P2 * P3 * P4 * P5 * P6 * ZPJ(:,:,NLON,NLAT) )
     
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      ! Update J-Values
      ZPJ(:,:,NLON,NLAT) = JSUM(:,:)


      ! Return to caller
      END SUBROUTINE MMRAN_16
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: matin4
!
! !DESCRIPTION: Subroutine MATIN4 used L-U decomposition to invert a
!  4x4 input matrix A.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE MATIN4 (A)
!
! !INPUT PARAMETERS:
!
      REAL*8, INTENT(INOUT) :: A(4,4)
! 
! !REVISION HISTORY: 
!  10 Feb 2005 - R. Yantosca - Initial version
!  09 Sep 2010 - R. Yantosca - Added ProTeX headers
!  25 Mar 2013 - S. D Eastham - Adapted from GET\_MIN\_TPAUSE\_LEVEL
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! MATIN4 begins here!
      !=================================================================

!---SETUP L AND U
      A(2,1) = A(2,1)/A(1,1)
      A(2,2) = A(2,2)-A(2,1)*A(1,2)
      A(2,3) = A(2,3)-A(2,1)*A(1,3)
      A(2,4) = A(2,4)-A(2,1)*A(1,4)
      A(3,1) = A(3,1)/A(1,1)
      A(3,2) = (A(3,2)-A(3,1)*A(1,2))/A(2,2)
      A(3,3) = A(3,3)-A(3,1)*A(1,3)-A(3,2)*A(2,3)
      A(3,4) = A(3,4)-A(3,1)*A(1,4)-A(3,2)*A(2,4)
      A(4,1) = A(4,1)/A(1,1)
      A(4,2) = (A(4,2)-A(4,1)*A(1,2))/A(2,2)
      A(4,3) = (A(4,3)-A(4,1)*A(1,3)-A(4,2)*A(2,3))/A(3,3)
      A(4,4) = A(4,4)-A(4,1)*A(1,4)-A(4,2)*A(2,4)-A(4,3)*A(3,4)
!---INVERT L
      A(4,3) = -A(4,3)
      A(4,2) = -A(4,2)-A(4,3)*A(3,2)
      A(4,1) = -A(4,1)-A(4,2)*A(2,1)-A(4,3)*A(3,1)
      A(3,2) = -A(3,2)
      A(3,1) = -A(3,1)-A(3,2)*A(2,1)
      A(2,1) = -A(2,1)
!---INVERT U
      A(4,4) = 1.D0/A(4,4)
      A(3,4) = -A(3,4)*A(4,4)/A(3,3)
      A(3,3) = 1.D0/A(3,3)
      A(2,4) = -(A(2,3)*A(3,4)+A(2,4)*A(4,4))/A(2,2)
      A(2,3) = -A(2,3)*A(3,3)/A(2,2)
      A(2,2) = 1.D0/A(2,2)
      A(1,4) = -(A(1,2)*A(2,4)+A(1,3)*A(3,4)+A(1,4)*A(4,4))/A(1,1)
      A(1,3) = -(A(1,2)*A(2,3)+A(1,3)*A(3,3))/A(1,1)
      A(1,2) = -A(1,2)*A(2,2)/A(1,1)
      A(1,1) = 1.D0/A(1,1)
!---MULTIPLY (U-INVERSE)*(L-INVERSE)
      A(1,1) = A(1,1)+A(1,2)*A(2,1)+A(1,3)*A(3,1)+A(1,4)*A(4,1)
      A(1,2) = A(1,2)+A(1,3)*A(3,2)+A(1,4)*A(4,2)
      A(1,3) = A(1,3)+A(1,4)*A(4,3)
      A(2,1) = A(2,2)*A(2,1)+A(2,3)*A(3,1)+A(2,4)*A(4,1)
      A(2,2) = A(2,2)+A(2,3)*A(3,2)+A(2,4)*A(4,2)
      A(2,3) = A(2,3)+A(2,4)*A(4,3)
      A(3,1) = A(3,3)*A(3,1)+A(3,4)*A(4,1)
      A(3,2) = A(3,3)*A(3,2)+A(3,4)*A(4,2)
      A(3,3) = A(3,3)+A(3,4)*A(4,3)
      A(4,1) = A(4,4)*A(4,1)
      A(4,2) = A(4,4)*A(4,2)
      A(4,3) = A(4,4)*A(4,3)
      RETURN

      END SUBROUTINE MATIN4
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: miesct
!
! !DESCRIPTION: Subroutine MIESCT is an adaptation of the Prather radiative
!  transfer code.
!     Prather, 1974, Astrophys. J., 192, 787-792. 
!        Sol'n of inhomogeneous Rayleigh scattering atmosphere.
!        (original Rayleigh w/ polarization)
!     Cochran and Trafton, 1978, Ap. J., 219, 756-762
!        Raman scattering in the atmospheres of the major planets.
!        (first use of anisotropic code)
!     Jacob, Gottlieb and Prather, 1989, J. Geophys. Res., 94, 12975-13002
!        Chemistry of a polluted cloudy boundary layer.
!        (documentation of extension to anisotropic scattering)
!  Takes atmospheric structure and source terms from standard J-code.
!  Limited to 4 Gauss points, only calculates mean rad. field.
!  Initialize variables FIXED/UNUSED in this special version:
!  FTOP = 1.0 = astrophysical flux (unit of pi) at SZA, -ZU0, use for scaling
!  FBOT = 0.0 = external isotropic flux on lower boundary 
!  SISOTP = 0.0 = Specific Intensity of isotropic radiation incident from top
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE MIESCT
!
! !USES:
!
! !REVISION HISTORY:
!  27 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER   :: I, ID, IM
      REAL*8    :: CMEQ1

      !=================================================================
      !MIESCT begins here!
      !=================================================================

!-----------------------------------------------------------------------
!---fix scattering to 4 Gauss pts = 8-stream
      CALL GAUSSP (N,EMU,WT)
!---solve eqn of R.T. only for first-order M=1
!      ZFLUX = (ZU0*FZ(ND)*ZREFL+FBOT)/(1.0d0+ZREFL)
      ZFLUX = (ZU0*FZ(ND)*ZREFL)/(1.0d0+ZREFL)
      M=1
      DO I=1,N
        CALL LEGND0 (EMU(I),PM0,MFIT)
        DO IM=M,MFIT
          PM(I,IM) = PM0(IM)
        ENDDO
      ENDDO
!
      CMEQ1 = 0.25D0
      CALL LEGND0 (-ZU0,PM0,MFIT)
      DO IM=M,MFIT
        PM0(IM) = CMEQ1*PM0(IM)
      ENDDO
!
      CALL BLKSLV
!
      DO ID=1,ND,2
        FJ(ID) = 4.0d0*FJ(ID) + FZ(ID)
      ENDDO

      RETURN

      END SUBROUTINE MIESCT
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: 
!
! !DESCRIPTION: NEW Mie code for J's, only uses 8-term expansion, 4-Gauss pts
!  Currently allow up to NP aerosol phase functions (at all altitudes) to
!  be associated with optical depth AER(1:NC) = aerosol opt.depth @ 1000 nm
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OPMIE( KW, WAVEL, XQO2, XQO3, FMEAN )
!
! !USES:
!
#     include "define.h"

      IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER KW
      REAL*8  WAVEL
      REAL*8  XQO2(NB)
      REAL*8  XQO3(NB)
      REAL*8  FMEAN(LLPAR)
!
! !REMARKS:
!  Pick Mie-wavelength with phase function and Qext:
!
!  01 RAYLE = Rayleigh phase
!  02 ISOTR = isotropic
!  03 ABSRB = fully absorbing 'soot', wavelength indep.
!  04 S_Bkg = backgrnd stratospheric sulfate (n=1.46,log-norm:r=.09um/sigma=.6)
!  05 S_Vol = volcanic stratospheric sulfate (n=1.46,log-norm:r=.08um/sigma=.8)
!  06 W_H01 = water haze (H1/Deirm.) (n=1.335, gamma:  r-mode=0.1um /alpha=2)
!  07 W_H04 = water haze (H1/Deirm.) (n=1.335, gamma:  r-mode=0.4um /alpha=2)
!  08 W_C02 = water cloud (C1/Deirm.) (n=1.335, gamma:  r-mode=2.0um /alpha=6)
!  09 W_C04 = water cloud (C1/Deirm.) (n=1.335, gamma:  r-mode=4.0um /alpha=6)
!  10 W_C08 = water cloud (C1/Deirm.) (n=1.335, gamma:  r-mode=8.0um /alpha=6)
!  11 W_C13 = water cloud (C1/Deirm.) (n=1.335, gamma:  r-mode=13.3um /alpha=6)
!  12 W_L06 = water cloud (Lacis) (n=1.335, r-mode=5.5um / alpha=11/3)
!  13 Ice-H = hexagonal ice cloud (Mishchenko)
!  14 Ice-I = irregular ice cloud (Mishchenko)
!                                                                             .
!  Choice of aerosol index MIEDX is made in SET_AER; optical depths are
!  apportioned to the AER array in SET_PROF
!                                                                             .
!-----------------------------------------------------------------------
!  FUNCTION RAYLAY(WAVE)---RAYLEIGH CROSS-SECTION for wave > 170 nm
!       WSQI = 1.E6/(WAVE*WAVE)
!       REFRM1 = 1.0E-6*(64.328+29498.1/(146.-WSQI)+255.4/(41.-WSQI))
!       RAYLAY = 5.40E-21*(REFRM1*WSQI)**2
!-----------------------------------------------------------------------
!                                                                             .
!     DTAUX    Local optical depth of each CTM level
!     PIRAY    Contribution of Rayleigh scattering to extinction
!     PIAER    Contribution of Aerosol scattering to extinction
!     TTAU     Optical depth of air vertically above each point (to top of atm)
!     FTAU     Attenuation of solar beam
!     POMEGA   Scattering phase function
!     FMEAN    Mean actinic flux at desired levels! 
! 
! !REVISION HISTORY: 
!  01 Oct 1995 - R. Yantosca - Initial version
!  31 Jul 2012 - R. Yantosca - Added ProTeX headers
!  27 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer jndlev(LLPAR),jaddlv(nc),jaddto(nc+1)
      integer km,i,j,k,l,ix,j1
      real*8 QXMIE(MX),XLAER(MX),SSALB(MX)
      real*8 xlo2,xlo3,xlray,xltau,zk,taudn,tauup,zk2
      real*8 POMEGAJ(2*M__,NC+1)
      real*8 DTAUX(NB),PIRAY(NB),PIAER(MX,NB),TTAU(NC+1),FTAU(NC+1)
      real*8 ftaulog,dttau,dpomega(2*M__)
      real*8 ftaulog2,dttau2,dpomega2(2*M__)

      ! For KLUDGE to fix the # of added levels (phs, 7/1/08)
      INTEGER :: loc(1)
c
!---Pick nearest Mie wavelength, no interpolation--------------
                              KM=1
      if( WAVEL .gt. 355.d0 ) KM=2
      if( WAVEL .gt. 500.d0 ) KM=3
!     if( WAVEL .gt. 800.d0 ) KM=4  !drop the 1000 nm wavelength
c
!---For Mie code scale extinction at 1000 nm to wavelength WAVEL (QXMIE)
      do I=1,MX
        QXMIE(I) = QAA(KM,MIEDX(I))/QAA(4,MIEDX(I))
        SSALB(I) = SSA(KM,MIEDX(I))
      enddo
c
!---Reinitialize arrays
      do j=1,nc+1
        ttau(j)=0.d0
        ftau(j)=0.d0
      enddo
c
!---Set up total optical depth over each CTM level, DTAUX
      J1 = NLBATM
      do J=J1,NB
        XLO3=DO3(J)*XQO3(J)
        XLO2=DM(J)*XQO2(J)*0.20948d0
        XLRAY=DM(J)*QRAYL(KW)
c  Zero absorption for testing purposes
c        call NOABS(XLO3,XLO2,XLRAY,AER(1,j),RFLECT)
        do I=1,MX
          XLAER(I)=AER(I,J)*QXMIE(I)
        enddo
c  Total optical depth from all elements
        DTAUX(J)=XLO3+XLO2+XLRAY
        do I=1,MX
          DTAUX(J)=DTAUX(J)+XLAER(I)
        enddo
c  Fractional extinction for Rayleigh scattering and each aerosol type
        PIRAY(J)=XLRAY/DTAUX(J)
        do I=1,MX
          PIAER(I,J)=SSALB(I)*XLAER(I)/DTAUX(J)
        enddo
      enddo
c
!---Define the scattering phase fn. with mix of Rayleigh(1) & Mie(MIEDX)
!   No. of quadrature pts fixed at 4 (M__), expansion of phase fn @ 8
      N = M__
      MFIT = 2*M__
      do j=j1,NB
        do i=1,MFIT
          pomegaj(i,j) = PIRAY(J)*PAA(I,KM,1)
          do k=1,MX
            pomegaj(i,j) = pomegaj(i,j) + PIAER(K,J)*PAA(I,KM,MIEDX(K))
          enddo
        enddo
      enddo
c
!---Calculate attenuated incident beam EXP(-TTAU/U0) and flux on surface
      do J=J1,NB
        if(AMF(J,J).gt.0.0D0) then
          XLTAU=0.0D0
          do I=1,NB
            XLTAU=XLTAU + DTAUX(I)*AMF(I,J)
          enddo
          if(XLTAU.gt.450.d0) then   ! for compilers with no underflow trapping
            FTAU(j)=0.d0
          else
            FTAU(J)=DEXP(-XLTAU)
          endif
        else
          FTAU(J)=0.0D0
        endif
      enddo
      if(U0.gt.0.D0) then
        ZFLUX = U0*FTAU(J1)*RFLECT/(1.d0+RFLECT)
      else
        ZFLUX = 0.d0
      endif
c
!------------------------------------------------------------------------
c  Take optical properties on CTM layers and convert to a photolysis
c  level grid corresponding to layer centres and boundaries. This is
c  required so that J-values can be calculated for the centre of CTM
c  layers; the index of these layers is kept in the jndlev array.
!------------------------------------------------------------------------
c
c  Set lower boundary and levels to calculate J-values at 
      J1=2*J1-1
      do j=1,LLPAR
        jndlev(j)=2*j
      enddo
c
c  Calculate column optical depths above each level, TTAU
      TTAU(NC+1)=0.0D0
      do J=NC,J1,-1
        I=(J+1)/2
        TTAU(J)=TTAU(J+1) + 0.5d0*DTAUX(I)
        jaddlv(j)=int(0.5d0*DTAUX(I)/dtaumax)
c  Subdivide cloud-top levels if required
! NOTE: Don't add more than DTAUSUB-1 (=9) sublevels (phs)
        if(jadsub(j).gt.0) then
          jadsub(j)=min(jaddlv(j)+1,nint(dtausub))*(nint(dsubdiv)-1) 
          jaddlv(j)=jaddlv(j)+jadsub(j)
        endif
      enddo
c
c  Calculate attenuated beam, FTAU, level boundaries then level centres
      FTAU(NC+1)=1.0D0
      do J=NC-1,J1,-2
        I=(J+1)/2
        FTAU(J)=FTAU(I)
      enddo
      do J=NC,J1,-2
        FTAU(J)=sqrt(FTAU(J+1)*FTAU(J-1))
      enddo
c
c  Calculate scattering properties, level centres then level boundaries
c  using an inverse interpolation to give correctly-weighted values
      do j=NC,J1,-2
        do i=1,MFIT
          pomegaj(i,j) = pomegaj(i,j/2)
        enddo
      enddo
      do j=J1+2,nc,2
        taudn = ttau(j-1)-ttau(j)
        tauup = ttau(j)-ttau(j+1)
        do i=1,MFIT
          pomegaj(i,j) = (pomegaj(i,j-1)*taudn + 
     $                    pomegaj(i,j+1)*tauup) / (taudn+tauup)
        enddo
      enddo
c  Define lower and upper boundaries
      do i=1,MFIT
        pomegaj(i,J1)   = pomegaj(i,J1+1)
        pomegaj(i,nc+1) = pomegaj(i,nc)
      enddo
c
!------------------------------------------------------------------------
c  Calculate cumulative total and define levels we want J-values at.
c  Sum upwards for levels, and then downwards for Mie code readjustments.
c
c     jaddlv(i)   Number of new levels to add between (i) and (i+1)
c     jaddto(i)   Total number of new levels to add to and above level (i)
c     jndlev(j)   Level needed for J-value for CTM layer (j)
c
!------------------------------------------------------------------------
c
c  Reinitialize level arrays
      do j=1,nc+1
        jaddto(j)=0
      enddo
c
      jaddto(J1)=jaddlv(J1)
      do j=J1+1,nc
        jaddto(j)=jaddto(j-1)+jaddlv(j)
      enddo

!==============================================================================
! KLUDGE TO LIMIT THE NUMBER OF ADDED LEVELS (phs, 7/1/08)
!
! PART 1: We need to replace the .gt. with .ge in this IF test
!
      if((jaddto(nc)+nc).GE.nl) then
         write(6,*)  jaddto(nc)+nc, 'NL',NL
!         write(6,1500)  jaddto(nc)+nc, 'NL',NL
!
! PART 2: We just trim the largest JADDLV until the condition is satisfied 
!         instead of simply stopping.  Remove the STOP statement.
!
         ! trim
         do while( (SUM( jaddlv(J1:nc) ) + NC) >= NL )
            loc=maxloc(jaddlv)
            jaddlv(loc(1))=jaddlv(loc(1))-1
         enddo

         ! then refill JADDTO
         jaddto(J1)=jaddlv(J1)
         do j=J1+1,nc
            jaddto(j)=jaddto(j-1)+jaddlv(j)
         enddo
         
!        ! Debug: double check
!        write(6,*) jaddto(nc)+nc
!        if((jaddto(nc)+nc).gt.nl) 
!     &      write(6,*)'OPMIE kludge: trap not working'
!==============================================================================
      endif

c     write(6,1300) jndlev
c     write(6,1300) jaddto
      do i=1,LLPAR
         jndlev(i)=jndlev(i)+jaddto(jndlev(i)-1)
      enddo

      ! this is just a transposition of the jaddto vector (phs)
      jaddto(nc)=jaddlv(nc)
      do j=nc-1,J1,-1
         jaddto(j)=jaddto(j+1)+jaddlv(j)
      enddo
c     write(6,1300) jndlev
c     write(6,1300) jaddto
c
!---------------------SET UP FOR MIE CODE-------------------------------
c
c  Transpose the ascending TTAU grid to a descending ZTAU grid.
c  Double the resolution - TTAU points become the odd points on the
c  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.
c  Odd point added at top of grid for unattenuated beam   (Z='inf')
c  
c        Surface:   TTAU(1)   now use ZTAU(2*NC+1)
c        Top:       TTAU(NC)  now use ZTAU(3)
c        Infinity:            now use ZTAU(1)
c
c  Mie scattering code only used from surface to level NC
!------------------------------------------------------------------------
!
c  Initialise all Fast-J optical property arrays
      do k=1,N__
        do i=1,MFIT
          pomega(i,k) = 0.d0
        enddo
        ztau(k) = 0.d0
        fz(k)   = 0.d0
      enddo
c
c  Ascend through atmosphere transposing grid and adding extra points
      do j=J1,nc+1
        k = 2*(nc+1-j)+2*jaddto(j)+1
        ztau(k)= ttau(j)
        fz(k)  = ftau(j)
        do i=1,MFIT
          pomega(i,k) = pomegaj(i,j)
        enddo
      enddo
c
c  Check profiles if desired
c      ND = 2*(NC+jaddto(J1)-J1)  + 3
c      if(kw.eq.1) call CH_PROF
c
!------------------------------------------------------------------------
c    Insert new levels, working downwards from the top of the atmosphere
c  to the surface (down in 'j', up in 'k'). This allows ztau and pomega
c  to be incremented linearly (in a +ve sense), and the flux fz to be
c  attenuated top-down (avoiding problems where lower level fluxes are
c  zero).
c
c    zk        fractional increment in level
c    dttau     change in ttau per increment    (linear, positive)
c    dpomega   change in pomega per increment  (linear)
c    ftaulog   change in ftau per increment    (exponential, normally < 1)
c
!------------------------------------------------------------------------
c
      do j=nc,J1,-1
          zk = 0.5d0/(1.d0+dble(jaddlv(j)-jadsub(j)))
          dttau = (ttau(j)-ttau(j+1))*zk
          do i=1,MFIT
            dpomega(i) = (pomegaj(i,j)-pomegaj(i,j+1))*zk
          enddo
c  Filter attenuation factor - set minimum at 1.0d-05
          if(ftau(j+1).eq.0.d0) then
            ftaulog=0.d0
          else
            ftaulog = ftau(j)/ftau(j+1)
            if(ftaulog.lt.1.d-150) then
              ftaulog=1.0d-05
            else
              ftaulog=exp(log(ftaulog)*zk)
            endif
          endif
          k = 2*(nc-j+jaddto(j)-jaddlv(j))+1   !  k at level j+1
          l = 0
c  Additional subdivision of first level if required
          if(jadsub(j).ne.0) then
            l=jadsub(j)/nint(dsubdiv-1)
            zk2=1.d0/dsubdiv
            dttau2=dttau*zk2
            ftaulog2=ftaulog**zk2
            do i=1,MFIT
              dpomega2(i)=dpomega(i)*zk2
            enddo
            do ix=1,2*(jadsub(j)+l)
              ztau(k+1) = ztau(k) + dttau2
              fz(k+1) = fz(k)*ftaulog2
              do i=1,MFIT
                pomega(i,k+1) = pomega(i,k) + dpomega2(i)
              enddo
              k = k+1
            enddo
          endif
          l = 2*(jaddlv(j)-jadsub(j)-l)+1
c
c  Add values at all intermediate levels
          do ix=1,l
            ztau(k+1) = ztau(k) + dttau
            fz(k+1) = fz(k)*ftaulog
            do i=1,MFIT
              pomega(i,k+1) = pomega(i,k) + dpomega(i)
            enddo
            k = k+1
          enddo
c
c  Alternate method to attenuate fluxes, fz, using 2nd-order finite
c  difference scheme - just need to comment in section below
c          ix = 2*(jaddlv(j)-jadsub(j))+1
c          if(l.le.0) then
c            l=k-ix-1
c          else
c            l=k-ix
c          endif
c          call efold(ftau(j+1),ftau(j),ix+1,fz(l))
c          if(jadsub(j).ne.0) then
c            k = 2*(nc-j+jaddto(j)-jaddlv(j))+1 !  k at level j+1
c            ix=2*(jadsub(j)+(jadsub(j)/nint(dsubdiv-1)))
c            call efold(ftau(j+1),fz(k+ix),ix,fz(k))
c          endif
c
      enddo
c
!---Update total number of levels and check doesn't exceed N__
      ND = 2*(NC+jaddto(J1)-J1)  + 3
!==============================================================================
! KLUDGE TO LIMIT THE NUMBER OF ADDED LEVELS (phs, 7/1/08)
!
! PART 3: Test to make sure that we haven't added more levels than the
!         dimension of the common block (i.e. ND <= N__).
!         
!         NOTE: this test should always be passed now that .ge. is 
!         used instead of .gt. in PART 1.
!
      if(nd.gt.N__) then
         write(6,1500) ND, 'N__',N__
         stop
      endif
!==============================================================================
c
!---Add boundary/ground layer to ensure no negative J's caused by
!---too large a TTAU-step in the 2nd-order lower b.c.
      ZTAU(ND+1) = ZTAU(ND)*1.000005d0
      ZTAU(ND+2) = ZTAU(ND)*1.000010d0
      zk=max(abs(U0),0.01d0)
      zk=dexp(-ZTAU(ND)*5.d-6/zk)
      FZ(ND+1) = FZ(ND)*zk
      FZ(ND+2) = FZ(ND+1)*zk
      do I=1,MFIT
        POMEGA(I,ND+1)   = POMEGA(I,ND)
        POMEGA(I,ND+2)   = POMEGA(I,ND)
      enddo
      ND = ND+2
c
      ZU0 = U0
      ZREFL = RFLECT
c
!-----------------------------------------
      CALL MIESCT
!-----------------------------------------
c  Accumulate attenuation for selected levels
      l=2*(NC+jaddto(J1))+3
      do j=1,LLPAR
        k=l-(2*jndlev(j))
        if(k.gt.ND-2) then
          FMEAN(j) = 0.d0
        else
          FMEAN(j) = FJ(k)
        endif
      enddo
c
      return
 1000 format(1x,i3,3(2x,1pe10.4),1x,i3)
 1300 format(1x,50(i3))
 1500 format(' Too many levels in photolysis code: need ',i5,' but ',a,
     $       ' dimensioned as ',i5)
      END SUBROUTINE OPMIE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: photoj
!
! !DESCRIPTION: Subroutine PHOTOJ is the driver routine for the FAST-J 
!  photolysis package.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE PHOTOJ( NLON,   NLAT,    YLAT,    DAY_OF_YR, 
     &                   MONTH,  DAY,     CSZA,    T,    
     &                   SA,     SAFE_OD, OPTDUST, OPTAER,
     &                   O3COL,  am_I_Root                    )
!
! !USES:
!
#     include "define.h"

      USE CMN_SIZE_MOD, ONLY : NDUST, NAER, NRH,  LLPAR

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)    :: NLON                   ! Grid box lon index
      INTEGER, INTENT(IN)    :: NLAT                   ! Grid box lat index 
      REAL*8,  INTENT(IN)    :: YLAT                   ! Latitude [degrees]
      INTEGER, INTENT(IN)    :: DAY_OF_YR              ! Day of year
      INTEGER, INTENT(IN)    :: MONTH                  ! Current month
      INTEGER, INTENT(IN)    :: DAY                    ! Day of month
      REAL*8,  INTENT(IN)    :: CSZA                   ! Cosine(SZA) [unitless]
      REAL*8,  INTENT(IN)    :: T(LLPAR)               ! Temperature [K]
      REAL*8,  INTENT(IN)    :: SA                     ! UV albedo [unitless]
      REAL*8,  INTENT(IN)    :: SAFE_OD(LLPAR)         ! Visible OD [unitless]
!      REAL*8,  INTENT(INOUT) :: OD(LLPAR)              ! Visible OD [unitless]
      REAL*8,  INTENT(IN)    :: O3COL                  ! Overhead O3 column [DU]
      LOGICAL, INTENT(IN)    :: am_I_Root              ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(INOUT) :: OPTDUST(LLPAR,NDUST)   ! Dust OD [unitless]
      REAL*8,  INTENT(INOUT) :: OPTAER(LLPAR,NAER*NRH) ! Aerosol OD [unitless]
!
! !AUTHOR:
!  Oliver Wild & Michael Prather
!
! !REMARKS:
!  New FAST J-Value code, troposphere only (mjprather 6/96); uses special 
!  wavelength quadrature spectral data (jv_spec.dat) that includes only 
!  289 nm - 800 nm  (later a single 205 nm add-on); uses special compact Mie 
!  code based on Feautrier/Auer/Prather vers. 
!                                                                             .
!  Important variables from other modules:
!  (1) ZJ   : Column array for J-values 
!  (2) ZPJ  : Global array for J-values (passed to SMVGEAR)
!  (3) JPNL : # of GEOS-CHEM layers in which to compute J-values
!  (4) JPPJ : # of photolysis rxns for FAST-J
!                                                                             .
!  NOTE: The value of PI listed here is slightly different than the value
!  in CMN_GCTM_mod.F.  The last digit is 4, whereas in CMN_GCTM_mod.F, the
!  last digit is 3.  Keep for now during testing of grid-independent code, but
!  this may be something to revisit at a later data. (bmy, 3/6/12)
!
! !REVISION HISTORY: 
!  01 Jun 1996 - M. Prather & O. Wild - Initial version
!  (1 ) Renamed NSLON to NLON and NSLAT to NLAT.  Now add DAY_OF_YR 
!        (formerly IDAY) and DAY to the arg list.  Swap places in arg list 
!        of SA and OD.  Now pass NLON, NLAT, DAY_OF_YR and DAY to "set_prof.f".
!        Added standard documentation header; cosmetic changes. (bmy, 7/15/03)
!  (2 ) We don't need to pass "P" via the arg list (bmy, 2/13/07)
!  06 Mar 2012 - R. Yantosca - Now pass O3COL via the arg list
!  06 Mar 2012 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  27 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I, J
      REAL*8, PARAMETER :: PI = 3.14159265358979324D0
      REAL*8            :: OD(LLPAR)

      !=================================================================
      ! PHOTOJ begins here!
      !=================================================================

      ! Copy over inputs
      OD = SAFE_OD

      ! Zero ZJ (column J-value array) and ZPJ (global J-value array)
      DO I = 1, JPNL
      DO J = 1, JPPJ
         ZJ(I,J)            = 0.D0
         ZPJ(I,J,NLON,NLAT) = 0.D0
      ENDDO
      ENDDO

      ! Import the cosine of the SZA from the CTM (bmy, 9/10/99)
      U0  = CSZA
      SZA = ACOS(CSZA) * ( 180.0d0 / PI )

      !-----------------------------------------------------------------
      !### If you want to set SZA = 0 degrees for testing,
      !### then uncomment the following lines (bmy, 9/13/99) 
      !U0  = 1.0d0
      !SZA = 0.0d0
      !-----------------------------------------------------------------

      ! Return if sun is below the horizon
      IF ( SZA > SZAMAX ) RETURN

      ! Set up Air, O3, BC profiles on GEOS-CHEM vertical levels
      CALL SET_PROF( NLON, NLAT, YLAT, MONTH,   DAY, 
     &               T,    SA,   OD,   OPTDUST, OPTAER, O3COL )

      ! Compute actinic flux at each GEOS-CHEM vertical level
      CALL JVALUE( SA, am_I_Root )
      
      ! Calculate J-values for all species
      CALL JRATET( T, DAY_OF_YR )

      ! ZJ is the J-value array for this column only
      ! Store in ZPJ (global array) for passing to SMVGEAR
      DO I = 1, JPNL
      DO J = 1, JPPJ
         ZPJ(I,J,NLON,NLAT) = ZJ(I,J)
      ENDDO
      ENDDO

      END SUBROUTINE PHOTOJ
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_aod
!
! !DESCRIPTION: Subroutine RD\_AOD reads aerosol phase functions that are
!  used to scale diagnostic output to an arbitrary wavelengh.  This
!  facilitates comparing with satellite observations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RD_AOD( NJ1, NAMFIL, am_I_Root )
!
! !USES:
!
#     include "define.h"

      USE ERROR_MOD,  ONLY : ERROR_STOP
      USE FILE_MOD,   ONLY : IOERROR

      IMPLICIT NONE

!
! !INPUT PARAMETERS:
!     
      INTEGER,          INTENT(IN) :: NJ1         ! Unit # of file to open
      CHARACTER(LEN=*), INTENT(IN) :: NAMFIL      ! Name of file to open
      LOGICAL,          INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
!
! !REMARKS:
!  The jv_spec_aod.dat file contains the optical properties for aerosols 
!  at a single wavelength to be used in the online calculation of the aerosol 
!  optical depth diagnostics.  The default properties are provided at 550 nm. 
!  These properties have been calculated using the same size and optical 
!  properties as the jv_spec.dat file used for the FAST-J photolysis 
!  calculations.  The user can exchange this set of properties with those at 
!  another wavelength.  We recommend that the wavelength used be included 
!  in the first line of the header for traceability (this line is output to 
!  the GEOS-Chem log file during run time). A complete set of optical 
!  properties from 250-2000 nm for aerosols is available at:  
!  ftp://ftp.as.harvard.edu/geos-chem/data/aerosol_optics/hi_spectral_res
!                                                                             .
!     -- Colette L. Heald, 05/10/10)
!
!  Important variables:
!                                                                             .
!     NAMFIL       Name of spectral data file (jv_spec_aod.dat)
!     NJ1          Channel number for reading data file
!     NAA2         Number of categories for scattering phase functions
!     QAA_AOD      Aerosol scattering phase functions
!     WAA_AOD      Wavelengths for the NK supplied phase functions
!     PAA_AOD      Phase function: first 8 terms of expansion
!     RAA_AOD      Effective radius associated with aerosol type
!     SSA_AOD      Single scattering albedo
!
! !REVISION HISTORY: 
!  10 May 2010 - C. Heald    - Initial version
!  06 Aug 2010 - C. Carouge  - Add an error check when opening the file
!  01 Aug 2012 - R. Yantosca - Now restore NJ1 to INTENT(IN) status
!  28 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
      INTEGER :: I, J, K, NAA2
      INTEGER :: IOS

      !================================================================
      ! RD_AOD begins here!
      !================================================================

      ! open file
      OPEN( NJ1, FILE=TRIM( NAMFIL ), STATUS='OLD', IOSTAT=IOS )

      ! Error check
      IF ( IOS /= 0 ) THEN
         IF ( am_I_Root ) THEN
            WRITE(6,100) trim(NAMFIL)
 100        FORMAT('Error opening filename=', a )
            CALL FLUSH(6)
         ENDIF
         CALL IOERROR( IOS, NJ1, 'RD_AOD:1')
      ENDIF

      
      ! Read header lines
      READ( NJ1,'(A)' ) TITLE0
      IF  ( am_I_Root ) WRITE( 6, '(1X,A)' ) TITLE0
      READ( NJ1,'(A)' ) TITLE0

      ! Read aerosol phase functions (one wavelength only):
      READ( NJ1,'(A10,I5,/)' ) TITLE0,NAA2
      DO j = 15, NAA
         READ(NJ1,110) TITLEA(j)
 110     FORMAT( 3x, a20 )
         IF ( am_I_Root ) WRITE(6,*) TITLEA(j)
         READ(NJ1,*) WAA_AOD(j),QAA_AOD(j),RAA_AOD(j),SSA_AOD(j),
     &               (PAA_AOD(i,j),i=1,8)
      ENDDO

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, '(a)' ) 'Aerosol Qext for AOD calculations'
         DO J=15,NAA
            WRITE( 6, * ) TITLEA(J),J,'  Qext =',(QAA_AOD(J))
         ENDDO
      ENDIF

      ! Close file
      CLOSE( NJ1 )

      END SUBROUTINE RD_AOD
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: rd_js
!
! !DESCRIPTION: Rereads the "ratj.d" file to map photolysis rate to reaction
!  Read in quantum yield 'jfacta' and fastj labels 'jlabel'
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RD_JS( NJ1, NAMFIL, am_I_Root )
!
! !USES:
!
#     include "define.h"

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      INTEGER,          INTENT(INOUT) :: nj1         ! Logical unit #
      CHARACTER(LEN=*), INTENT(IN)    :: namfil      ! File name
      LOGICAL,          INTENT(IN)    :: am_I_Root   ! Is this the root CPU? 
!
! !INPUT PARAMETERS: 
!
!
! !REMARKS:
!     jfacta    Quantum yield (or multiplication factor) for photolysis
!     jlabel    Reference label identifying appropriate J-value to use
!     ipr       Photolysis reaction counter - should total 'jppj'
! 
! !REVISION HISTORY: 
!  01 Jun 1998 - P. Murti    - Initial version
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  30 Jul 2012 - R. Yantosca - Add reference to findFreeLUN under DEVEL tag
!   3 Aug 2012 - R. Yantosca - Restore NJ1 to INTENT(IN) status
!  27 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer       ipr, i
      character*120 cline
c
c Reread the ratj.d file to map photolysis rate to reaction
c                     Read in quantum yield jfacta and fastj label jlabel
      ipr=0
      ! open file
      open(nj1,file=namfil,status='old',form='formatted')
 10   read(nj1,'(a)',err=20) cline
      if(cline(2:5).eq.'9999') then
         go to 20
      elseif(cline(1:1).eq.'#') then
         go to 10
      elseif(cline(5:5).eq.'$') then
         go to 10
      else
         ipr=ipr+1
         read(cline(79:83),'(f5.1)') jfacta(ipr)
         read(cline(86:92),'(a7)')   jlabel(ipr)
         jfacta(ipr)=jfacta(ipr)/100.d0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Additional code to read reaction names and branch numbers
!  (ppm, 6/98, bmy, 9/99)     
!         read (cline(7:10),"(a4)") rnames(ipr)
! FP increase length of name (6/2009)
         read (cline(7:13),"(a7)") rnames(ipr)
         rnames(ipr) = trim(rnames(ipr))
         branch(ipr) = 1
         do i=1,ipr-1
            if (rnames(ipr) == rnames(i)) branch(ipr) = branch(i) + 1
         enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         go to 10
      endif
 20   close(nj1)
      if(ipr.ne.jppj) then
         IF ( am_I_Root ) write(6,1000) ipr,jppj
         stop
      endif
c
c Print details to standard output
      IF ( am_I_Root ) THEN
         write(6,1100) ipr
         write(6,1200) (i, jlabel(i), jfacta(i),i=1,ipr)
      ENDIF
c
      return
 1000 format(' Error: ',i3,' photolysis labels but ',i3,' reactions')
 1100 format(' Fast-J Photolysis Scheme: considering ',i2,' reactions')
 1200 format(3x,10(3(i2,': ',a7,' (Q.Y. ',f5.3,') '),/,3x))
      END SUBROUTINE RD_JS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_jv_atms_dat
!
! !DESCRIPTION: Reads the default T and O3 profiles for FAST-J photolysis.
!  This replaces the obsolete rd\_prof.F routine, which read from the ASCII
!  file "jv\_atms.dat".
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_JV_ATMS_DAT( am_I_Root )
!
! !USES:
!
      ! Modules for netCDF read
      USE m_netcdf_io_open
      USE m_netcdf_io_read
      USE m_netcdf_io_readattr
      USE m_netcdf_io_close
   
      ! GEOS-Chem modules
      USE DIRECTORY_MOD, ONLY : DATA_DIR_1x1    ! Data directory 
         
      IMPLICIT NONE

# include "netcdf.inc"

! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root          ! Is this the root CPU?
!
! !REMARKS:
!  This file was automatically generated by the Perl scripts in the
!  NcdfUtilities package (which ships w/ GEOS-Chem) and was subsequently
!  hand-edited.
!
! !REVISION HISTORY:
!  19 Apr 2012 - R. Yantosca - Initial version
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  27 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! netCDF file handle
      INTEGER            :: fId                 ! netCDF file ID
    
      ! Character strings
      CHARACTER(LEN=255) :: nc_dir              ! netCDF directory name
      CHARACTER(LEN=255) :: nc_file             ! netCDF file name
      CHARACTER(LEN=255) :: nc_path             ! netCDF path name
      CHARACTER(LEN=255) :: v_name              ! netCDF variable name 
      CHARACTER(LEN=255) :: a_name              ! netCDF attribute name
      CHARACTER(LEN=255) :: a_val               ! netCDF attribute value
    
      ! Arrays for netCDF start and count values
      INTEGER            :: st3d(3), ct3d(3)    ! For 3D arrays 
    
      !=========================================================================
      ! Open and read data from the netCDF file
      !=========================================================================
    
      ! Directory and file names
      nc_dir  = TRIM( DATA_DIR_1x1 ) // 'FastJ_201204/'      
      nc_file = 'fastj.jv_atms_dat.nc'  
      nc_path = TRIM( nc_dir ) // TRIM( nc_file )
    
      ! Open netCDF file
      CALL Ncop_Rd( fId, TRIM(nc_path) )
      
      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) REPEAT( '%', 79 )
         WRITE( 6, 110 ) TRIM(nc_file)
         WRITE( 6, 120 ) TRIM(nc_dir)
      ENDIF
    
      !----------------------------------------
      ! VARIABLE: T
      !----------------------------------------
      
      ! Variable name
      v_name = "T"
      
      ! Read T from file
      st3d   = (/  1,  1,  1 /)
      ct3d   = (/ 51, 18, 12 /)
      CALL NcRd( TREF, fId, TRIM(v_name), st3d, ct3d )
      
      ! Read the T:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
      
      ! Echo info to stdout
      IF ( am_I_Root ) THEN 
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF
    
      !----------------------------------------
      ! VARIABLE: O3
      !----------------------------------------
      
      ! Variable name
      v_name = "O3"
      
      ! Read O3 from file
      st3d   = (/  1,  1,  1 /)
      ct3d   = (/ 51, 18, 12 /)
      CALL NcRd( OREF, fId, TRIM(v_name), st3d, ct3d )
      
      ! Read the O3:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
          
      ! Echo info to stdout
      IF ( am_I_Root ) THEN 
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF
    
      !=================================================================
      ! Cleanup and quit
      !=================================================================
    
      ! Close netCDF file
      CALL NcCl( fId )
    
      ! Echo info to stdout
      IF ( am_I_Root ) THEN 
         WRITE( 6, 140 )
         WRITE( 6, 100 ) REPEAT( '%', 79 )
      ENDIF
    
      ! FORMAT statements
 100  FORMAT( a                                              )
 110  FORMAT( '%% Opening file  : ',         a               )
 120  FORMAT( '%%  in directory : ',         a, / , '%%'     )
 130  FORMAT( '%% Successfully read ',       a, ' [', a, ']' )
 140  FORMAT( '%% Successfully closed file!'                 )
      
      END SUBROUTINE READ_JV_ATMS_DAT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: rd_xxx
!
! !DESCRIPTION: Read in wavelength bins, solar fluxes, Rayleigh parameters, 
!  temperature-dependent cross sections and Rayleigh/aerosol scattering phase 
!  functions with temperature dependences.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RD_XXX( NJ1, NAMFIL, am_I_Root )
!
! !USES:
!
#     include "define.h"

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      INTEGER,          INTENT(IN) :: NJ1         ! Logical unit #
      CHARACTER(LEN=*), INTENT(IN) :: NAMFIL      ! File name
      LOGICAL,          INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! !INPUT PARAMETERS: 
!
!
! !REMARKS:
!     NAMFIL   Name of spectral data file (jv_spec.dat)
!     NJ1      Channel number for reading data file
!     NJVAL    Number of species to calculate J-values for
!     NWWW     Number of wavelength bins, from NW1:NW2
!     WBIN     Boundaries of wavelength bins
!     WL       Centres of wavelength bins - 'effective wavelength'
!     FL       Solar flux incident on top of atmosphere (cm-2.s-1)
!     QRAYL    Rayleigh parameters (effective cross-section) (cm2)
!     QBC      Black Carbon absorption extinct. (specific cross-sect.) (m2/g)
!     QO2      O2 cross-sections
!     QO3      O3 cross-sections
!     Q1D      O3 => O(1D) quantum yield
!     TQQ      Temperature for supplied cross sections
!     QQQ      Supplied cross sections in each wavelength bin (cm2)
!     NAA      Number of categories for scattering phase functions
!     QAA      Aerosol scattering phase functions
!     NK       Number of wavelengths at which functions supplied (set as 4)
!     WAA      Wavelengths for the NK supplied phase functions
!     PAA      Phase function: first 8 terms of expansion
!     RAA      Effective radius associated with aerosol type
!     SSA      Single scattering albedo
!     npdep    Number of pressure dependencies
!     zpdep    Pressure dependencies by wavelength bin
!     jpdep    Index of cross sections requiring pressure dependence
!     lpdep    Label for pressure dependence
!
! !REVISION HISTORY: 
!  01 Jun 1998 - P. Murti - Initial version
!  (1 ) Updated to include new pressure-dependancy function for GLYX and MGLY. 
!        (tmf, 1/7/09)
!  (2 ) Added a pressure-dependancy function selector 'pdepf'. 
!        (tmf, ccc, 1/7/09)
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  30 Jul 2012 - R. Yantosca - Added ProTeX headers
!  01 Aug 2012 - R. Yantosca - Now restore NJ1 to INTENT(IN) status
!  27 Mar 2013 - S. D. Eastham - Renamed from RD_TJPL to RD_XXX and rolled
!                                into FAST_JX_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
      integer i, j, k, iw, nk, nqqq, nwww
      character*7  lpdep(7)
      do J=1,NS
        do K=1,3
          TQQ(K,J) = 0.d0
        enddo
      enddo
!-------------spectral data---------------------------------------------
      open(NJ1, FILE=NAMFIL)
      read(NJ1,'(A)') TITLE0
      write(6,'(1X,A)') TITLE0
      read(NJ1,'(10X,14I5)') NJVAL,NWWW,NW1,NW2
      if(NJVAL.gt.NS) then
        IF ( am_I_Root ) write(6,300) NJVAL,NS
        stop
      endif
!------------NQQQ = no. additional J-values from X-sects (O2,O3P,O3D+NQQQ)
!- NQQQ is changed to NJVAL-1 because there are 2 dummy species at the end
! used for acetone pressure dependency only. (ccc, 4/20/09)
      NQQQ = NJVAL-1
      read(NJ1,102) (WBIN(IW),IW=1,NWWW)
      read(NJ1,102) (WBIN(IW+1),IW=1,NWWW)
      read(NJ1,102) (WL(IW),IW=1,NWWW)
      read(NJ1,102) (FL(IW),IW=1,NWWW)
      read(NJ1,102) (QRAYL(IW),IW=1,NWWW)
      read(NJ1,102) (QBC(IW),IW=1,NWWW)   !  From Loiusse et al. [JGR, 1996]
c
!---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
      do K=1,3
        read(NJ1,103) TITLEJ(K,1),TQQ(K,1), (QO2(IW,K),IW=1,NWWW)
      enddo
      do K=1,3
        read(NJ1,103) TITLEJ(K,2),TQQ(K,2), (QO3(IW,K),IW=1,NWWW)
      enddo
      do K=1,3
        read(NJ1,103) TITLEJ(K,3),TQQ(K,3), (Q1D(IW,K),IW=1,NWWW)
      enddo
      do K=1,3
        write(6,200) titlej(1,k),(tqq(i,k),i=1,3)
      enddo
c
!---Read remaining species:  X-sections at 2 T's
      do J=1,NQQQ
        read(NJ1,103) TITLEJ(1,J+3),TQQ(1,J+3),(QQQ(IW,1,J),IW=1,NWWW)
        read(NJ1,103) TITLEJ(2,J+3),TQQ(2,J+3),(QQQ(IW,2,J),IW=1,NWWW)
        IF ( am_I_Root ) THEN 
           write(6,200) titlej(1,j+3),(tqq(i,j+3),i=1,2)
        ENDIF
      enddo
      read(NJ1,'(A)') TITLE0
c
c---Pressure dependencies
      read(NJ1,104) npdep
      do k=1,npdep
         read(NJ1,105) lpdep(k), pdepf(k), (zpdep(iw,k),iw=1,nwww)
         IF ( am_I_Root ) THEN 
            write(6,201)  lpdep(k), pdepf(k), (zpdep(iw,k),iw=1,nwww)
         ENDIF

         !--------------------------------------
         ! Special treatment for MGLY pressure dependency
         ! (tmf, 11/16/06)
         !--------------------------------------
         if ( pdepf(k) .eq. 4 ) then           
            ! pass zpdep to mglypdep
            mglypdep(:,1) = zpdep(:,k)
            read(NJ1,105) lpdep(k), pdepf(k), (mglypdep(iw,2),iw=1,nwww)
            read(NJ1,105) lpdep(k), pdepf(k), (mglypdep(iw,3),iw=1,nwww)
         endif
      enddo
      read(NJ1,'(A)') TITLE0

c
c---Zero index arrays
      do j=1,jppj
        jind(j)=0
      enddo
      do j=1,NJVAL
        jpdep(j)=0
      enddo
c
!---Set mapping index
      do j=1,NJVAL
        do k=1,jppj
          if (jlabel(k).eq.titlej(1,j)) jind(k)=j
        enddo
        do k=1,npdep
          if (lpdep(k).eq.titlej(1,j)) jpdep(j)=k
        enddo
      enddo
      do k=1,jppj
        if(jfacta(k).eq.0.d0)
     &             write(6,*) 'Not using photolysis reaction ',k
        if(jind(k).eq.0) then
          if(jfacta(k).eq.0.d0) then
            jind(k)=1
          else
            write(6,*) 'Which J-rate for photolysis reaction ',k,' ?'
            stop
          endif
        endif
      enddo
c
!---Read aerosol phase functions:
      read(NJ1,'(A10,I5,/)') TITLE0,NAA
      NK=4        ! Fix number of wavelengths at 4
      do j=1,NAA
        read(NJ1,110) TITLEA(j)
        do k=1,NK
          read(NJ1,*) WAA(k,j),QAA(k,j),RAA(k,j),SSA(k,j),
     &                                             (PAA(i,k,j),i=1,8)
        enddo
      enddo
c
      IF ( am_I_Root ) THEN 
         write(6,*) 'Aerosol phase functions & wavelengths'
         do J=1,NAA
            write(6,'(1x,A8,I2,A,9F8.1)')
     &           TITLEA(J),J,'  wavel=',(WAA(K,J),K=1,NK)
            write(6,'(9x,I2,A,9F8.4)') J,'  Qext =',(QAA(K,J),K=1,NK)
         enddo
      ENDIF
!--------
! Modify reading and writing formats 105 & 201 for pressure dependancy 
c (ccc, 1/7/09)

  101 FORMAT(8E10.3)
  102 FORMAT(10X,7E10.3)
  103 FORMAT(A7,F3.0,7E10.3)
c 103 FORMAT(A7,F3.0,7E10.3/(10X,7E10.3))
  104 FORMAT(13x,i2)
  105 FORMAT(A7,2x,I1,7E10.3)
  110 format(3x,a20)
  200 format(1x,' x-sect:',a10,3(3x,f6.2))
  201 format(1x,' pr.dep:',a10,1x,I1,7(1pE10.3))
  300 format(' Number of x-sections supplied to Fast-J: ',i3,/,
     &       ' Maximum number allowed (NS) only set to: ',i3,
     &       ' - increase in jv_cmn.h')
      close(NJ1)
      return
      END SUBROUTINE RD_XXX
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: set_aer
!
! !DESCRIPTION: Routine SET\_AER sets aerosol/cloud types and define black 
!  carbon profile for the FAST-J photolysis scheme.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_AER( am_I_Root )
!
! !USES:
!
#     include "define.h"

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! !REMARKS:
!  MX       Number of different types of aerosol to be considered
!  MIEDX    Index of aerosol types in jv_spec.dat - hardwire in here
! 
! 
! !REVISION HISTORY: 
!  01 Jun 1999 - O. Wild     - Initial version
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!  30 Jul 2012 - R. Yantosca - Added ProTeX headers
!  27 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer i
c
c  Initialise aerosol index
      do i=1,MX
        MIEDX(i) = 0
      enddo
c
c  Select Aerosol/Cloud types to be used - define types here
c  Each of these types must be listed in the order used by OPMIE.F
      MIEDX(1)  =  3   !  Black carbon absorber
      MIEDX(2)  = 10   !  Water Cloud (Deirmenjian 8 micron)
      MIEDX(3)  = 14   !  Irregular Ice Cloud (Mishchenko)
      MIEDX(4)  = 15   !  Mineral Dust  .15 micron    (rvm, 9/30/00)
      MIEDX(5)  = 16   !  Mineral Dust  .25 micron    (rvm, 9/30/00)
      MIEDX(6)  = 17   !  Mineral Dust  .4  micron    (rvm, 9/30/00)
      MIEDX(7)  = 18   !  Mineral Dust  .8  micron    (rvm, 9/30/00)
      MIEDX(8)  = 19   !  Mineral Dust 1.5  micron    (rvm, 9/30/00)
      MIEDX(9)  = 20   !  Mineral Dust 2.5  micron    (rvm, 9/30/00)
      MIEDX(10) = 21   !  Mineral Dust 4.0  micron    (rvm, 9/30/00)
      MIEDX(11) = 22   !  Tropospheric Sulfate, RH=0  (rvm, bmy, 2/27/02)
      MIEDX(12) = 23   !  Tropospheric Sulfate, RH=50 (rvm, bmy, 2/27/02)
      MIEDX(13) = 24   !  Tropospheric Sulfate, RH=70 (rvm, bmy, 2/27/02)
      MIEDX(14) = 25   !  Tropospheric Sulfate, RH=80 (rvm, bmy, 2/27/02)
      MIEDX(15) = 26   !  Tropospheric Sulfate, RH=90 (rvm, bmy, 2/27/02)
      MIEDX(16) = 29   !  Black Carbon,         RH=0  (rvm, bmy, 2/27/02)
      MIEDX(17) = 30   !  Black Carbon,         RH=50 (rvm, bmy, 2/27/02)
      MIEDX(18) = 31   !  Black Carbon,         RH=70 (rvm, bmy, 2/27/02)
      MIEDX(19) = 32   !  Black Carbon,         RH=80 (rvm, bmy, 2/27/02)
      MIEDX(20) = 33   !  Black Carbon,         RH=90 (rvm, bmy, 2/27/02)
      MIEDX(21) = 36   !  Organic Carbon,       RH=0  (rvm, bmy, 2/27/02)
      MIEDX(22) = 37   !  Organic Carbon,       RH=50 (rvm, bmy, 2/27/02)
      MIEDX(23) = 38   !  Organic Carbon,       RH=70 (rvm, bmy, 2/27/02)
      MIEDX(24) = 39   !  Organic Carbon,       RH=80 (rvm, bmy, 2/27/02)
      MIEDX(25) = 40   !  Organic Carbon,       RH=90 (rvm, bmy, 2/27/02)
      MIEDX(26) = 43   !  Sea Salt (accum),     RH=0  (rvm, bmy, 2/27/02)
      MIEDX(27) = 44   !  Sea Salt (accum),     RH=50 (rvm, bmy, 2/27/02)
      MIEDX(28) = 45   !  Sea Salt (accum),     RH=70 (rvm, bmy, 2/27/02)
      MIEDX(29) = 46   !  Sea Salt (accum),     RH=80 (rvm, bmy, 2/27/02)
      MIEDX(30) = 47   !  Sea Salt (accum),     RH=90 (rvm, bmy, 2/27/02)
      MIEDX(31) = 50   !  Sea Salt (coarse),    RH=0  (rvm, bmy, 2/27/02)
      MIEDX(32) = 51   !  Sea Salt (coarse),    RH=50 (rvm, bmy, 2/27/02)
      MIEDX(33) = 52   !  Sea Salt (coarse),    RH=70 (rvm, bmy, 2/27/02)
      MIEDX(34) = 53   !  Sea Salt (coarse),    RH=80 (rvm, bmy, 2/27/02)
      MIEDX(35) = 54   !  Sea Salt (coarse),    RH=90 (rvm, bmy, 2/27/02)
      
c
c  Ensure all 'MX' types are valid selections
      do i=1,MX
         IF ( am_I_Root ) THEN 
            write(6,1000) MIEDX(i),TITLEA(MIEDX(i))
         ENDIF
         if(MIEDX(i).gt.NAA.or.MIEDX(i).le.0) then
            IF ( am_I_Root ) THEN
               write(6,1200) MIEDX(i),NAA
            ENDIF
          stop
        endif
      enddo
c
c Approximate Black Carbon up to 10 km; surface 200 ng/m3  (Liousse et al)
c Scale: 1 ng/m3 = 1.0d-15 g/cm3 (1.0d-11 g/m2/cm as BREF is in cm))
c
c Simple place-holder profile
      do i=1,51
        BREF(i)=10.d0*1.0d-11
        if(i.gt.6) BREF(i)=0.d0
      enddo
c
      return
 1000 format('Using Aerosol type: ',i3,1x,a)
 1200 format('Aerosol type ',i3,' unsuitable; supplied values must be ',
     &       'between 1 and ',i3)
      END SUBROUTINE SET_AER
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
#     include "define.h"

      USE CMN_SIZE_MOD,  ONLY : LLPAR

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
!  27 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD
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
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sphere
!
! !DESCRIPTION: Subroutine SPHERE calculates spherical geometry; derives 
!  tangent heights, slant path lengths and air mass factor for each layer. 
!  Not called when SZA > 98 degrees.  Beyond 90 degrees, include treatment 
!  of emergent beam (where tangent height is below altitude J-value desired at).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SPHERE
!
! !USES:
!
#     include "define.h"

      USE CMN_SIZE_MOD, ONLY : LLPAR

      IMPLICIT NONE
!
! !REMARKS:
!  GMU     MU, cos(solar zenith angle)
!  RZ      Distance from centre of Earth to each point (cm)
!  RQ      Square of radius ratios
!  TANHT   Tangent height for the current SZA
!  XL      Slant path between points
!  AMF     Air mass factor for slab between level and level above
! 
! 
! !REVISION HISTORY: 
!         1997 - O. Wild     - Initial version
!  31 Jul 2012 - R. Yantosca - Added ProTeX headers
!  10 Aug 2012 - R. Yantosca - Replace LPAR with LLPAR
!  27 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer i, j, k, ii
      real*8 airmas, gmu, xmu1, xmu2, xl, diff
      REAL*8 Ux,H,RZ(NB),RQ(NB),ZBYR
c
c  Inlined air mass factor function for top of atmosphere
      AIRMAS(Ux,H) = (1.0d0+H)/SQRT(Ux*Ux+2.0d0*H*(1.0d0-
     $         0.6817d0*EXP(-57.3d0*ABS(Ux)/SQRT(1.0d0+5500.d0*H))/
     $                                             (1.0d0+0.625d0*H)))
c
      GMU = U0
      RZ(1)=RAD+Z(1)
      ZBYR = ZZHT/RAD
      DO 2 II=2,NB
        RZ(II) = RAD + Z(II)
        RQ(II-1) = (RZ(II-1)/RZ(II))**2
    2 CONTINUE
      IF (GMU.LT.0.0D0) THEN
        TANHT = RZ(nlbatm)/DSQRT(1.0D0-GMU**2)
      ELSE
        TANHT = RZ(nlbatm)
      ENDIF
c
c  Go up from the surface calculating the slant paths between each level
c  and the level above, and deriving the appropriate Air Mass Factor
      DO 16 J=1,NB
        DO K=1,NB
          AMF(K,J)=0.D0
        ENDDO
c
c  Air Mass Factors all zero if below the tangent height
        IF (RZ(J).LT.TANHT) GOTO 16
c  Ascend from layer J calculating AMFs
        XMU1=ABS(GMU)
        DO 12 I=J,LLPAR
          XMU2=DSQRT(1.0D0-RQ(I)*(1.0D0-XMU1**2))
          XL=RZ(I+1)*XMU2-RZ(I)*XMU1
          AMF(I,J)=XL/(RZ(I+1)-RZ(I))
          XMU1=XMU2
   12   CONTINUE
c  Use function and scale height to provide AMF above top of model
        AMF(NB,J)=AIRMAS(XMU1,ZBYR)
c
c  Twilight case - Emergent Beam
        IF (GMU.GE.0.0D0) GOTO 16
        XMU1=ABS(GMU)
c  Descend from layer J
        DO 14 II=J-1,1,-1
          DIFF=RZ(II+1)*DSQRT(1.0D0-XMU1**2)-RZ(II)
          if(II.eq.1) DIFF=max(DIFF,0.d0)   ! filter
c  Tangent height below current level - beam passes through twice
          IF (DIFF.LT.0.0D0) THEN
            XMU2=DSQRT(1.0D0-(1.0D0-XMU1**2)/RQ(II))
            XL=ABS(RZ(II+1)*XMU1-RZ(II)*XMU2)
            AMF(II,J)=2.d0*XL/(RZ(II+1)-RZ(II))
            XMU1=XMU2
c  Lowest level intersected by emergent beam
          ELSE
            XL=RZ(II+1)*XMU1*2.0D0
c            WTING=DIFF/(RZ(II+1)-RZ(II))
c            AMF(II,J)=(1.0D0-WTING)*2.D0**XL/(RZ(II+1)-RZ(II))
            AMF(II,J)=XL/(RZ(II+1)-RZ(II))
            GOTO 16
          ENDIF
   14   CONTINUE
c
   16 CONTINUE
      RETURN
      END SUBROUTINE SPHERE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: xsec1d
!
! !DESCRIPTION: Function XSEC1D calculates quantum yields for the photolytic
!  decomposition of O3 into O2 and O(1D), interpolated across 3 temperatures.
!\\
!\\
! !INTERFACE:
!
      REAL*8 FUNCTION XSEC1D(K,TTT)
!
! !USES:
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: K
      REAL*8,  INTENT(IN) :: TTT
!
! !REVISION HISTORY: 
!  27 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      !=================================================================
      ! XSEC1D begins here!
      !=================================================================

      XSEC1D =
     &  FLINT(TTT,TQQ(1,3),TQQ(2,3),TQQ(3,3),Q1D(K,1),Q1D(K,2),Q1D(K,3))

      END FUNCTION XSEC1D
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: xseco2
!
! !DESCRIPTION: Function XSECO2 calculates the photolytic cross-section of  
!  O2 across 3 temperatures (ignoring S-R bands).
!\\
!\\
! !INTERFACE:
!
      REAL*8 FUNCTION XSECO2(K,TTT)
!
! !USES:
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: K
      REAL*8,  INTENT(IN) :: TTT
!
! !REVISION HISTORY: 
!  27 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      !=================================================================
      ! XSECO2 begins here!
      !=================================================================

      XSECO2 =
     &  FLINT(TTT,TQQ(1,1),TQQ(2,1),TQQ(3,1),QO2(K,1),QO2(K,2),QO2(K,3))

      END FUNCTION XSECO2
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: xseco3
!
! !DESCRIPTION: Function XSECO3 calculates the photolytic cross-section of  
!  O3 across 3 temperatures.
!\\
!\\
! !INTERFACE:
!
      REAL*8 FUNCTION XSECO3(K,TTT)
!
! !USES:
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: K
      REAL*8,  INTENT(IN) :: TTT
!
! !REVISION HISTORY: 
!  27 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      !=================================================================
      ! XSECO3 begins here!
      !=================================================================

      XSECO3 =
     &  FLINT(TTT,TQQ(1,2),TQQ(2,2),TQQ(3,2),QO3(K,1),QO3(K,2),QO3(K,3))

      END FUNCTION XSECO3
!EOC
      END MODULE FAST_JX_MOD
