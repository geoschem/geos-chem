!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: regrid_a2a_mod.F90
!
! !DESCRIPTION: Module REGRID\_A2A\_MOD uses an algorithm adapted from 
!  MAP\_A2A code to regrid from one horizontal grid to another.
!\\
!\\
! !INTERFACE: 
!
MODULE Regrid_A2A_Mod
!
! !USES:
!
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Do_Regrid_A2A
  PUBLIC  :: Map_A2A
  PUBLIC  :: Init_Map_A2A
  PUBLIC  :: Cleanup_Map_A2A

  ! Map_A2A overloads these routines
  INTERFACE Map_A2A
    MODULE PROCEDURE Map_A2A_R8R8 
    MODULE PROCEDURE Map_A2A_R4R8
    MODULE PROCEDURE Map_A2A_R4R4
    MODULE PROCEDURE Map_A2A_R8R4
  END INTERFACE Map_A2A
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Read_Input_Grid
  PRIVATE :: Map_A2A_R8R8
  PRIVATE :: Map_A2A_R4R4
  PRIVATE :: Map_A2A_R4R8
  PRIVATE :: Map_A2A_R8R4
  PRIVATE :: Ymap_R8R8
  PRIVATE :: Ymap_R4R8
  PRIVATE :: Ymap_R4R4
  PRIVATE :: Ymap_R8R4
  PRIVATE :: Xmap_R8R8
  PRIVATE :: Xmap_R4R4
  PRIVATE :: Xmap_R4R8
  PRIVATE :: Xmap_R8R4
!
! !REVISION HISTORY:
!  13 Mar 2012 - M. Cooper   - Initial version
!  03 Apr 2012 - M. Payer    - Now use functions GET_AREA_CM2(I,J,L), 
!                              GET_YEDGE(I,J,L) and GET_YSIN(I,J,L) from the
!                              new grid_mod.F90
!  22 May 2012 - L. Murray   - Implemented several bug fixes
!  23 Aug 2012 - R. Yantosca - Add capability for starting from hi-res grids
!                              (generic 0.5x0.5, generic 0.25x0.25, etc.)
!  23 Aug 2012 - R. Yantosca - Add subroutine READ_INPUT_GRID, which reads the
!                              grid parameters (lon & lat edges) w/ netCDF
!  27 Aug 2012 - R. Yantosca - Now parallelize key DO loops
!  19 May 2014 - C. Keller   - MAP_A2A now accepts single and double precision
!                              input/output.
!  14 Jul 2014 - R. Yantosca - Now save IIPAR, JJPAR, OUTLON, OUTSIN, OUTAREA 
!                              as module variables.  This helps us remove a
!                              dependency for the HEMCO emissions package.
!                              input/output.
!  02 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!  11 Feb 2015 - C. Keller   - Add capability for regridding local grids onto
!                              global grids. To do so, xmap now only operates
!                              within the longitude range spanned by the input
!                              domain.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!

  !---------------------------------------------------------------------------
  ! These are now kept locally, to "shadow" variables from other parts of
  ! GEOS-Chem.  This avoids depending on GEOS-Chem code within the core
  ! HEMCO modules. (bmy, 7/14/14)
  !---------------------------------------------------------------------------
  CHARACTER(LEN=255)  :: NC_DIR        ! Directory w/ netCDF files
  INTEGER             :: IIPAR         ! # of longitudes (x-dimension) in grid
  INTEGER             :: JJPAR         ! # of latitudes  (y-dimension) in grid
  REAL(fp), ALLOCATABLE :: OUTLON (:  )  ! Longitude on output grid
  REAL(fp), ALLOCATABLE :: OUTSIN (:  )  ! Sines of latitudes on output grid
  REAL(fp), ALLOCATABLE :: OUTAREA(:,:)  ! Surface areas on output grid
!
! !DEFINED PARAMETERS:
!
  !---------------------------------------------------------------------------
  ! These were taken from CMN_GCTM_mod.F90.  This helps us to avoid depending 
  ! on GEOS-Chem modules in the core HEMCO modules.  (bmy, 7/14/14)
  !---------------------------------------------------------------------------
  REAL(fp), PARAMETER :: PI =   3.14159265358979323e+0_fp   ! Pi
  REAL(fp), PARAMETER :: Re =   6.375d6                 ! Earth radius [m]

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_regrid_a2a
!
! !DESCRIPTION: Subroutine DO\_REGRID\_A2A regrids 2-D data in the
!  horizontal direction.  This is a wrapper for the MAP\_A2A routine.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_REGRID_A2A( FILENAME, IM,      JM,              &
                            INGRID,   OUTGRID, IS_MASS, netCDF )
!
! !INPUT PARAMETERS:
!
    ! Name of file with lon and lat edge information on the INPUT GRID
    CHARACTER(LEN=*), INTENT(IN)    :: FILENAME

    ! Number of lon centers and lat centers on the INPUT GRID
    INTEGER,          INTENT(IN)    :: IM 
    INTEGER,          INTENT(IN)    :: JM

    ! Data array on the input grid
    REAL(fp),           INTENT(IN)    :: INGRID(IM,JM)

    ! IS_MASS=0 if data is units of concentration (molec/cm2/s, unitless, etc.)
    ! IS_MASS=1 if data is units of mass (kg/yr, etc.); we will need to convert
    !           INGRID to per unit area
    INTEGER,          INTENT(IN)    :: IS_MASS

    ! Read from netCDF file?  (needed for debugging, will disappear later)
    LOGICAL, OPTIONAL,INTENT(IN)    :: netCDF  
!
! !OUTPUT PARAMETERS:
!
    ! Data array on the OUTPUT GRID
    REAL(fp),           INTENT(OUT)   :: OUTGRID(IIPAR,JJPAR) 
!
! !REMARKS:
!  The netCDF optional argument is now obsolete, because we now always read
!  the grid definitions from netCDF files instead of ASCII.  Keep it for
!  the time being in order to avoid having to change many lines of code
!  everywhere.
!
! !REVISION HISTORY:
!  13 Mar 2012 - M. Cooper   - Initial version
!  22 May 2012 - L. Murray   - Bug fix: INSIN should be allocated w/ JM+1.
!  22 May 2012 - R. Yantosca - Updated comments, cosmetic changes
!  25 May 2012 - R. Yantosca - Bug fix: declare the INGRID argument as
!                              INTENT(IN) to preserve the values of INGRID
!                              in the calling routine
!  06 Aug 2012 - R. Yantosca - Now make IU_REGRID a local variable
!  06 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  23 Aug 2012 - R. Yantosca - Now use f10.4 format for hi-res grids
!  23 Aug 2012 - R. Yantosca - Now can read grid info from netCDF files
!  27 Aug 2012 - R. Yantosca - Add parallel DO loops
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS to describe parameter
!                              more clearly
!  15 Jul 2014 - R. Yantosca - Now use global module variables
!  15 Jul 2014 - R. Yantosca - Remove reading from ASCII input files
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER           :: I,        J
    INTEGER           :: IOS,      M
    INTEGER           :: IU_REGRID
    REAL(fp)            :: INAREA,   RLAT
    CHARACTER(LEN=15) :: HEADER1
    CHARACTER(LEN=20) :: FMT_LAT,  FMT_LON, FMT_LEN
    LOGICAL           :: USE_NETCDF

    ! Arrays
    REAL(fp)            :: INLON  (IM+1    )  ! Lon edges        on INPUT GRID
    REAL(fp)            :: INSIN  (    JM+1)  ! SIN( lat edges ) on INPUT GRID
    REAL(fp)            :: IN_GRID(IM, JM  )  ! Shadow variable for INGRID

    !======================================================================
    ! Initialization
    !======================================================================

    ! Read the grid specifications from a netCDF file
    CALL READ_INPUT_GRID( IM, JM, FILENAME, INLON, INSIN )

    !======================================================================
    ! Regridding
    !======================================================================

    ! Copy the input argument INGRID to a local shadow variable,
    ! so that we can preserve the value of INGRID in the calling routine
    IN_GRID = INGRID

    ! Convert input to per area units if necessary
    IF ( IS_MASS == 1 ) THEN

       !$OMP PARALLEL DO                   &
       !$OMP DEFAULT( SHARED             ) &
       !$OMP PRIVATE( I, J, RLAT, INAREA )
       DO J = 1, JM
          RLAT   = INSIN(J+1) - INSIN(J)
          INAREA = (2e+0_fp * PI * Re * RLAT * 1e+4_fp * Re) / DBLE(IM)
          DO I = 1, IM
             IN_GRID(I,J) = IN_GRID(I,J) / INAREA
          ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

    ! Call MAP_A2A to do the regridding
    CALL MAP_A2A( IM,    JM,    INLON,  INSIN,  IN_GRID,        &
                  IIPAR, JJPAR, OUTLON, OUTSIN, OUTGRID, 0, 0 )

    ! Convert back from "per area" if necessary
    IF ( IS_MASS == 1 ) THEN

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J   )
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          OUTGRID(I,J) = OUTGRID(I,J) * OUTAREA(I,J)
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

  END SUBROUTINE Do_Regrid_A2A
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Map_A2A_r8r8
!
! !DESCRIPTION: Subroutine MAP\_A2A\_R8R8 is a horizontal arbitrary grid to 
!  arbitrary grid conservative high-order mapping regridding routine by S-J 
!  Lin.  Both the input data and output data have REAL(fp) precision.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Map_A2A_r8r8( im, jm, lon1, sin1, q1, &
                           in, jn, lon2, sin2, q2, ig, iv, missval)
!
! !INPUT PARAMETERS:
!
    ! Longitude and Latitude dimensions of INPUT grid
    INTEGER, INTENT(IN)  :: im, jm

    ! Longitude and Latitude dimensions of OUTPUT grid
    INTEGER, INTENT(IN)  :: in, jn

    ! IG=0: pole to pole; 
    ! IG=1 J=1 is half-dy north of south pole
    INTEGER, INTENT(IN)  :: ig

    ! IV=0: Regrid scalar quantity
    ! IV=1: Regrid vector quantity
    INTEGER, INTENT(IN)  :: iv

    ! Longitude edges (degrees) of INPUT and OUTPUT grids
    REAL*8,  INTENT(IN)  :: lon1(im+1), lon2(in+1)

    ! Sine of Latitude Edges (radians) of INPUT and OUTPUT grids
    REAL*8,  INTENT(IN)  :: sin1(jm+1), sin2(jn+1)

    ! Quantity on INPUT grid
    REAL*8,  INTENT(IN)  :: q1(im,jm)

!
! !OUTPUT PARAMETERS:
!
    ! Regridded quantity on OUTPUT grid
    REAL*8,  INTENT(OUT) :: q2(in,jn)
!
! !OPTIONAL ARGUMENTS
!
    REAL*8,  INTENT(IN), OPTIONAL :: missval
!
! !REMARKS:
!  This routine is overloaded by the MAP_A2A interface.
!
! !REVISION HISTORY:
!  (1) Original subroutine by S-J Lin.  Converted to F90 freeform format
!      and inserted into "Geos3RegridModule" by Bob Yantosca (9/21/00)
!  (2) Added F90 type declarations to be consistent w/ TypeModule.f90.
!      Also updated comments. (bmy, 9/21/00)
!  21 Sep 2000 - R. Yantosca - Initial version
!  27 Aug 2012 - R. Yantosca - Add parallel DO loops
!  02 Mar 2015 - C. Keller   - Added optional argument missval
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: i,j,k
    REAL*8  :: qtmp(in,jm)

    ! Init
    IF ( PRESENT(missval) ) THEN
       qtmp = missval
       q2   = missval
    ELSE
       qtmp = 0.0d0
       q2   = 0.0d0
    ENDIF

    !===================================================================
    ! E-W regridding
    !===================================================================    
    IF ( im         == in         .and. &
         lon1(1)    == lon2(1)    .and. &
         lon1(im+1) == lon2(in+1)        ) THEN

       ! Don't call XMAP if both grids have the same # of longitudes
       ! but save the input data in the QTMP array
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J )
       DO j=1,jm-ig
       DO i=1,im
          qtmp(i,j+ig) = q1(i,j+ig)
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ELSE

       ! Otherwise, call XMAP to regrid in the E-W direction
       CALL xmap_r8r8(im, jm-ig, lon1, q1(1,1+ig),in, lon2, qtmp(1,1+ig) )

    ENDIF
    
    !===================================================================
    ! N-S regridding
    !===================================================================    
    IF ( jm         == jn         .and. & 
         sin1(1)    == sin2(1)    .and. &
         sin1(jm+1) == sin2(jn+1)        ) THEN

       ! Don't call XMAP if both grids have the same # of longitudes,
       ! but assign the value of QTMP to the output Q2 array
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J )      
       DO j=1,jm-ig
       DO i=1,in
          q2(i,j+ig) = qtmp(i,j+ig)
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ELSE

       ! Otherwise, call YMAP to regrid in the N-S direction
       CALL ymap_r8r8(in, jm, sin1, qtmp(1,1+ig), jn, sin2, q2(1,1+ig), ig, iv)

    ENDIF

  END SUBROUTINE Map_A2A_r8r8
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Map_A2A_R4R4
!
! !DESCRIPTION: Subroutine MAP\_A2A\_R4R4 is a horizontal arbitrary grid 
!  to arbitrary grid conservative high-order mapping regridding routine 
!  by S-J Lin.  Both the input and output data have REAL*4 precision.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Map_A2A_r4r4( im, jm, lon1, sin1, q1, &
                           in, jn, lon2, sin2, q2, ig, iv, missval)
!
! !INPUT PARAMETERS:
!
    ! Longitude and Latitude dimensions of INPUT grid
    INTEGER, INTENT(IN)  :: im, jm

    ! Longitude and Latitude dimensions of OUTPUT grid
    INTEGER, INTENT(IN)  :: in, jn

    ! IG=0: pole to pole; 
    ! IG=1 J=1 is half-dy north of south pole
    INTEGER, INTENT(IN)  :: ig

    ! IV=0: Regrid scalar quantity
    ! IV=1: Regrid vector quantity
    INTEGER, INTENT(IN)  :: iv

    ! Longitude edges (degrees) of INPUT and OUTPUT grids
    REAL*4,  INTENT(IN)  :: lon1(im+1), lon2(in+1)

    ! Sine of Latitude Edges (radians) of INPUT and OUTPUT grids
    REAL*4,  INTENT(IN)  :: sin1(jm+1), sin2(jn+1)

    ! Quantity on INPUT grid
    REAL*4,  INTENT(IN)  :: q1(im,jm)
!
! !OUTPUT PARAMETERS:
!
    ! Regridded quantity on OUTPUT grid
    REAL*4,  INTENT(OUT) :: q2(in,jn)
!
! !OPTIONAL ARGUMENTS
!
    REAL*4,  INTENT(IN), OPTIONAL :: missval
!
! !REMARKS:
!  This routine is overloaded by the MAP_A2A interface.
!
! !REVISION HISTORY:
!  (1) Original subroutine by S-J Lin.  Converted to F90 freeform format
!      and inserted into "Geos3RegridModule" by Bob Yantosca (9/21/00)
!  (2) Added F90 type declarations to be consistent w/ TypeModule.f90.
!      Also updated comments. (bmy, 9/21/00)
!  21 Sep 2000 - R. Yantosca - Initial version
!  27 Aug 2012 - R. Yantosca - Add parallel DO loops
!  02 Mar 2015 - C. Keller   - Added optional argument missval
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: i,j,k
    REAL*4  :: qtmp(in,jm)

    ! Init
    IF ( PRESENT(missval) ) THEN
       qtmp = missval
       q2   = missval
    ELSE
       qtmp = 0.0
       q2   = 0.0
    ENDIF

    !===================================================================
    ! E-W regridding
    !===================================================================    
    IF ( im         == in         .and. &
         lon1(1)    == lon2(1)    .and. &
         lon1(im+1) == lon2(in+1)        ) THEN

       ! Don't call XMAP if both grids have the same # of longitudes
       ! but save the input data in the QTMP array
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J )
       DO j=1,jm-ig
       DO i=1,im
          qtmp(i,j+ig) = q1(i,j+ig)
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ELSE

       ! Otherwise, call XMAP to regrid in the E-W direction
       CALL xmap_r4r4(im, jm-ig, lon1, q1(1,1+ig),in, lon2, qtmp(1,1+ig) )

    ENDIF

    !===================================================================
    ! N-S regridding
    !===================================================================    
    IF ( jm         == jn         .and. & 
         sin1(1)    == sin2(1)    .and. &
         sin1(jm+1) == sin2(jn+1)        ) THEN

       ! Don't call XMAP if both grids have the same # of longitudes,
       ! but assign the value of QTMP to the output Q2 array
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J )      
       DO j=1,jm-ig
       DO i=1,in
          q2(i,j+ig) = qtmp(i,j+ig)
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ELSE

       ! Otherwise, call YMAP to regrid in the N-S direction
       CALL ymap_r4r4(in, jm, sin1, qtmp(1,1+ig), jn, sin2, q2(1,1+ig), ig, iv)

    ENDIF

  END SUBROUTINE Map_A2A_r4r4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Map_A2A_r4r8
!
! !DESCRIPTION: Subroutine MAP\_A2A\_R4R8 is a horizontal arbitrary grid to 
!  arbitrary grid conservative high-order mapping regridding routine by
!  S-J Lin.  The input data has REAL*4 precision, but the output argument
!  has REAL(fp) precision.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Map_A2A_r4r8( im, jm, lon1, sin1, q1, &
                           in, jn, lon2, sin2, q2, ig, iv, missval)
!
! !INPUT PARAMETERS:
!
    ! Longitude and Latitude dimensions of INPUT grid
    INTEGER, INTENT(IN)  :: im, jm

    ! Longitude and Latitude dimensions of OUTPUT grid
    INTEGER, INTENT(IN)  :: in, jn

    ! IG=0: pole to pole; 
    ! IG=1 J=1 is half-dy north of south pole
    INTEGER, INTENT(IN)  :: ig

    ! IV=0: Regrid scalar quantity
    ! IV=1: Regrid vector quantity
    INTEGER, INTENT(IN)  :: iv

    ! Longitude edges (degrees) of INPUT and OUTPUT grids
    REAL*4,  INTENT(IN)  :: lon1(im+1), lon2(in+1)

    ! Sine of Latitude Edges (radians) of INPUT and OUTPUT grids
    REAL*4,  INTENT(IN)  :: sin1(jm+1), sin2(jn+1)

    ! Quantity on INPUT grid
    REAL*4,  INTENT(IN)  :: q1(im,jm)
!
! !OUTPUT PARAMETERS:
!
    ! Regridded quantity on OUTPUT grid
    REAL*8,  INTENT(OUT) :: q2(in,jn)
!
! !OPTIONAL ARGUMENTS
!
    REAL*8,  INTENT(IN), OPTIONAL :: missval
!
! !REMARKS:
!  This routine is overloaded by the MAP_A2A interface.
!
! !REVISION HISTORY:
!  (1) Original subroutine by S-J Lin.  Converted to F90 freeform format
!      and inserted into "Geos3RegridModule" by Bob Yantosca (9/21/00)
!  (2) Added F90 type declarations to be consistent w/ TypeModule.f90.
!      Also updated comments. (bmy, 9/21/00)
!  21 Sep 2000 - R. Yantosca - Initial version
!  27 Aug 2012 - R. Yantosca - Add parallel DO loops
!  02 Mar 2015 - C. Keller   - Added optional argument missval
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: i,j,k
    REAL*8  :: qtmp(in,jm)

    ! Init
    IF ( PRESENT(missval) ) THEN
       qtmp = missval
       q2   = missval
    ELSE
       qtmp = 0.0d0
       q2   = 0.0d0
    ENDIF

    !===================================================================
    ! E-W regridding
    !===================================================================    
    IF ( im         == in         .and. &
         lon1(1)    == lon2(1)    .and. &
         lon1(im+1) == lon2(in+1)        ) THEN

       ! Don't call XMAP if both grids have the same # of longitudes
       ! but save the input data in the QTMP array
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J )
       DO j=1,jm-ig
       DO i=1,im
          qtmp(i,j+ig) = q1(i,j+ig)
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ELSE

       ! Otherwise, call XMAP to regrid in the E-W direction
       CALL xmap_r4r8(im, jm-ig, lon1, q1(1,1+ig),in, lon2, qtmp(1,1+ig) )

    ENDIF
    
    !===================================================================
    ! N-S regridding
    !===================================================================    
    IF ( jm         == jn         .and. & 
         sin1(1)    == sin2(1)    .and. &
         sin1(jm+1) == sin2(jn+1)        ) THEN

       ! Don't call XMAP if both grids have the same # of longitudes,
       ! but assign the value of QTMP to the output Q2 array
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J )      
       DO j=1,jm-ig
       DO i=1,in
          q2(i,j+ig) = qtmp(i,j+ig)
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ELSE

       ! Otherwise, call YMAP to regrid in the N-S direction
       CALL ymap_r4r8(in, jm, sin1, qtmp(1,1+ig), jn, sin2, q2(1,1+ig), ig, iv)

    ENDIF

  END SUBROUTINE Map_A2A_r4r8
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Map_A2A_r8r4
!
! !DESCRIPTION: Subroutine MAP\_A2A\_R8R4 is a horizontal arbitrary grid to 
!  arbitrary grid conservative high-order mapping regridding routine by
!  S-J Lin.  The input data has REAL*8 precision, but the output argument
!  has REAL*4 precision.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Map_A2A_r8r4( im, jm, lon1, sin1, q1, &
                           in, jn, lon2, sin2, q2, ig, iv, missval)
!
! !INPUT PARAMETERS:
!
    ! Longitude and Latitude dimensions of INPUT grid
    INTEGER, INTENT(IN)  :: im, jm

    ! Longitude and Latitude dimensions of OUTPUT grid
    INTEGER, INTENT(IN)  :: in, jn

    ! IG=0: pole to pole; 
    ! IG=1 J=1 is half-dy north of south pole
    INTEGER, INTENT(IN)  :: ig

    ! IV=0: Regrid scalar quantity
    ! IV=1: Regrid vector quantity
    INTEGER, INTENT(IN)  :: iv

    ! Longitude edges (degrees) of INPUT and OUTPUT grids
    REAL*4,  INTENT(IN)  :: lon1(im+1), lon2(in+1)

    ! Sine of Latitude Edges (radians) of INPUT and OUTPUT grids
    REAL*4,  INTENT(IN)  :: sin1(jm+1), sin2(jn+1)

    ! Quantity on INPUT grid
    REAL*8,  INTENT(IN)  :: q1(im,jm)
!
! !OUTPUT PARAMETERS:
!
    ! Regridded quantity on OUTPUT grid
    REAL*4,  INTENT(OUT) :: q2(in,jn)
!
! !OPTIONAL ARGUMENTS
!
    REAL*4,  INTENT(IN), OPTIONAL :: missval
!
! !REMARKS:
!  This routine is overloaded by the MAP_A2A interface.
!
! !REVISION HISTORY:
!  (1) Original subroutine by S-J Lin.  Converted to F90 freeform format
!      and inserted into "Geos3RegridModule" by Bob Yantosca (9/21/00)
!  (2) Added F90 type declarations to be consistent w/ TypeModule.f90.
!      Also updated comments. (bmy, 9/21/00)
!  21 Sep 2000 - R. Yantosca - Initial version
!  27 Aug 2012 - R. Yantosca - Add parallel DO loops
!  02 Mar 2015 - C. Keller   - Added optional argument missval
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: i,j,k
    REAL*4  :: qtmp(in,jm)

    ! Init
    IF ( PRESENT(missval) ) THEN
       qtmp = missval
       q2   = missval
    ELSE
       qtmp = 0.0
       q2   = 0.0
    ENDIF

    !===================================================================
    ! E-W regridding
    !===================================================================    
    IF ( im         == in         .and. &
         lon1(1)    == lon2(1)    .and. &
         lon1(im+1) == lon2(in+1)        ) THEN

       ! Don't call XMAP if both grids have the same # of longitudes
       ! but save the input data in the QTMP array
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J )
       DO j=1,jm-ig
       DO i=1,im
          qtmp(i,j+ig) = q1(i,j+ig)
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ELSE

       ! Otherwise, call XMAP to regrid in the E-W direction
       CALL xmap_r8r4(im, jm-ig, lon1, q1(1,1+ig),in, lon2, qtmp(1,1+ig) )

    ENDIF
    
    !===================================================================
    ! N-S regridding
    !===================================================================    
    IF ( jm         == jn         .and. & 
         sin1(1)    == sin2(1)    .and. &
         sin1(jm+1) == sin2(jn+1)        ) THEN

       ! Don't call XMAP if both grids have the same # of longitudes,
       ! but assign the value of QTMP to the output Q2 array
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J )      
       DO j=1,jm-ig
       DO i=1,in
          q2(i,j+ig) = qtmp(i,j+ig)
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ELSE

       ! Otherwise, call YMAP to regrid in the N-S direction
       CALL ymap_r4r4(in, jm, sin1, qtmp(1,1+ig), jn, sin2, q2(1,1+ig), ig, iv)

    ENDIF

  END SUBROUTINE Map_A2A_r8r4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ymap_r8r8
!
! !DESCRIPTION: Routine to perform area preserving mapping in N-S from an 
!  arbitrary resolution to another.  Both the input and output arguments
!  have REAL(fp) precision.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ymap_r8r8(im, jm, sin1, q1, jn, sin2, q2, ig, iv)
!
! !INPUT PARAMETERS:
!
    ! original E-W dimension
    INTEGER, INTENT(IN)  :: im            
    
    ! original N-S dimension
    INTEGER, INTENT(IN)  :: jm            

    ! Target N-S dimension
    INTEGER, INTENT(IN)  :: jn           
    
    ! IG=0: scalars from SP to NP (D-grid v-wind is also IG=0)
    ! IG=1: D-grid u-wind
    INTEGER, INTENT(IN)  :: ig            
  
    ! IV=0: scalar; 
    ! IV=1: vector
    INTEGER, INTENT(IN)  :: iv            
  
    ! Original southern edge of the cell sin(lat1)  
    REAL*8,  INTENT(IN)  :: sin1(jm+1-ig) 
    
    ! Original data at center of the cell
    REAL*8,  INTENT(IN)  :: q1(im,jm)      
    
    ! Target cell's southern edge sin(lat2)
    REAL*8,  INTENT(IN)  :: sin2(jn+1-ig) 
!
! !OUTPUT PARAMETERS:
!
    ! Mapped data at the target resolution
    REAL*8,  INTENT(OUT) :: q2(im,jn)     
!
! !REMARKS:
!
!   sin1 (1) = -1 must be south pole; sin1(jm+1)=1 must be N pole.
!
!   sin1(1) < sin1(2) < sin1(3) < ... < sin1(jm) < sin1(jm+1)
!   sin2(1) < sin2(2) < sin2(3) < ... < sin2(jn) < sin2(jn+1)!
!
! !AUTHOR:
!   Developer: Prasad Kasibhatla
!   March 6, 2012
!
! !REVISION HISTORY
!  06 Mar 2012 - P. Kasibhatla - Initial version
!  27 Aug 2012 - R. Yantosca   - Added parallel DO loops
!  27 Aug 2012 - R. Yantosca   - Change REAL*4 variables to REAL(fp) to better
!                                ensure numerical stability
!  31 Mar 2014 - C. Keller     - Initialize qsum to zero to avoid undefined 
!                                values in nested grids
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: i, j0, m, mm, j
    REAL*8               :: dy1(jm)
    REAL*8               :: dy
    REAL*8               :: qsum, sum
    
    ! YMAP begins here!
    do j=1,jm-ig
       dy1(j) = sin1(j+1) - sin1(j)
    enddo

    !===============================================================
    ! Area preserving mapping
    !===============================================================
    
    !$OMP PARALLEL DO                          &
    !$OMP DEFAULT( SHARED                    ) &
    !$OMP PRIVATE( I, J0, J, M, QSUM, MM, DY )
    do 1000 i=1,im
       qsum = 0.0d0
       j0 = 1
       do 555 j=1,jn-ig
       do 100 m=j0,jm-ig
             
          !=========================================================
          ! locate the southern edge: sin2(i)
          !=========================================================
          if(sin2(j) .ge. sin1(m) .and. sin2(j) .le. sin1(m+1)) then
             
             if(sin2(j+1) .le. sin1(m+1)) then
                
                ! entire new cell is within the original cell
                q2(i,j)=q1(i,m)
                j0 = m
                goto 555
             else
                
                ! South most fractional area
                qsum=(sin1(m+1)-sin2(j))*q1(i,m)
                
                do mm=m+1,jm-ig
                   
                   ! locate the northern edge: sin2(j+1)
                   if(sin2(j+1) .gt. sin1(mm+1) ) then
                      
                      ! Whole layer
                      qsum = qsum + dy1(mm)*q1(i,mm)
                   else
                      
                      ! North most fractional area
                      dy = sin2(j+1)-sin1(mm)
                      qsum=qsum+dy*q1(i,mm)
                      j0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100    continue
123    q2(i,j) = qsum / ( sin2(j+1) - sin2(j) )
555    continue
1000 continue
     !$OMP END PARALLEL DO

     !===================================================================
     ! Final processing for poles
     !===================================================================
     if ( ig .eq. 0 .and. iv .eq. 0 ) then
         
        ! South pole
        sum = 0.e+0_fp
        do i=1,im
           sum = sum + q2(i,1)
        enddo

        sum = sum / DBLE( im )
        do i=1,im
           q2(i,1) = sum
        enddo

        ! North pole:
        sum = 0.e+0_fp
        do i=1,im
           sum = sum + q2(i,jn)
        enddo

        sum = sum / DBLE( im )
        do i=1,im
           q2(i,jn) = sum
        enddo

     endif

   END SUBROUTINE ymap_r8r8
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ymap_r4r8
!
! !DESCRIPTION: Routine to perform area preserving mapping in N-S from an 
!  arbitrary resolution to another.  The input argument has REAL*4 precision
!  but the output argument has REAL(fp) precision.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ymap_r4r8(im, jm, sin1, q1, jn, sin2, q2, ig, iv)
!
! !INPUT PARAMETERS:
!

    ! original E-W dimension
    INTEGER, INTENT(IN)  :: im            
    
    ! original N-S dimension
    INTEGER, INTENT(IN)  :: jm            

    ! Target N-S dimension
    INTEGER, INTENT(IN)  :: jn           
    
    ! IG=0: scalars from SP to NP (D-grid v-wind is also IG=0)
    ! IG=1: D-grid u-wind
    INTEGER, INTENT(IN)  :: ig            
  
    ! IV=0: scalar; 
    ! IV=1: vector
    INTEGER, INTENT(IN)  :: iv            
  
    ! Original southern edge of the cell sin(lat1)  
    REAL*4,  INTENT(IN)  :: sin1(jm+1-ig) 
    
    ! Original data at center of the cell
    REAL*8,  INTENT(IN)  :: q1(im,jm)      
    
    ! Target cell's southern edge sin(lat2)
    REAL*4,  INTENT(IN)  :: sin2(jn+1-ig) 
!
! !OUTPUT PARAMETERS:
!
    ! Mapped data at the target resolution
    REAL*8,  INTENT(OUT) :: q2(im,jn)     
!
! !REMARKS:
!
!   sin1 (1) = -1 must be south pole; sin1(jm+1)=1 must be N pole.
!
!   sin1(1) < sin1(2) < sin1(3) < ... < sin1(jm) < sin1(jm+1)
!   sin2(1) < sin2(2) < sin2(3) < ... < sin2(jn) < sin2(jn+1)!
!
! !AUTHOR:
!   Developer: Prasad Kasibhatla
!   March 6, 2012
!
! !REVISION HISTORY
!  06 Mar 2012 - P. Kasibhatla - Initial version
!  27 Aug 2012 - R. Yantosca   - Added parallel DO loops
!  27 Aug 2012 - R. Yantosca   - Change REAL*4 variables to REAL(fp) to better
!                                ensure numerical stability
!  31 Mar 2014 - C. Keller     - Initialize qsum to zero to avoid undefined 
!                                values in nested grids
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: i, j0, m, mm, j
    REAL*8               :: dy1(jm)
    REAL*8               :: dy
    REAL*8               :: qsum, sum
    
    ! YMAP begins here!
    do j=1,jm-ig
       dy1(j) = sin1(j+1) - sin1(j)
    enddo

    !===============================================================
    ! Area preserving mapping
    !===============================================================
    
    !$OMP PARALLEL DO                          &
    !$OMP DEFAULT( SHARED                    ) &
    !$OMP PRIVATE( I, J0, J, M, QSUM, MM, DY )
    do 1000 i=1,im
       qsum = 0.0d0
       j0 = 1
       do 555 j=1,jn-ig
       do 100 m=j0,jm-ig
             
          !=========================================================
          ! locate the southern edge: sin2(i)
          !=========================================================
          if(sin2(j) .ge. sin1(m) .and. sin2(j) .le. sin1(m+1)) then
             
             if(sin2(j+1) .le. sin1(m+1)) then
                
                ! entire new cell is within the original cell
                q2(i,j)=q1(i,m)
                j0 = m
                goto 555
             else
                
                ! South most fractional area
                qsum=(sin1(m+1)-sin2(j))*q1(i,m)
                
                do mm=m+1,jm-ig
                   
                   ! locate the northern edge: sin2(j+1)
                   if(sin2(j+1) .gt. sin1(mm+1) ) then
                      
                      ! Whole layer
                      qsum = qsum + dy1(mm)*q1(i,mm)
                   else
                      
                      ! North most fractional area
                      dy = sin2(j+1)-sin1(mm)
                      qsum=qsum+dy*q1(i,mm)
                      j0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100    continue
123    q2(i,j) = qsum / ( sin2(j+1) - sin2(j) )
555    continue
1000 continue
     !$OMP END PARALLEL DO

     !===================================================================
     ! Final processing for poles
     !===================================================================
     if ( ig .eq. 0 .and. iv .eq. 0 ) then
         
        ! South pole
        sum = 0.e+0_fp
        do i=1,im
           sum = sum + q2(i,1)
        enddo

        sum = sum / DBLE( im )
        do i=1,im
           q2(i,1) = sum
        enddo

        ! North pole:
        sum = 0.e+0_fp
        do i=1,im
           sum = sum + q2(i,jn)
        enddo

        sum = sum / DBLE( im )
        do i=1,im
           q2(i,jn) = sum
        enddo

     endif

   END SUBROUTINE ymap_r4r8
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ymap_r8r4
!
! !DESCRIPTION: Routine to perform area preserving mapping in N-S from an 
!  arbitrary resolution to another.  The input argument has REAL*8 precision
!  but the output argument has REAL*4 precision.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ymap_r8r4(im, jm, sin1, q1, jn, sin2, q2, ig, iv)
!
! !INPUT PARAMETERS:
!

    ! original E-W dimension
    INTEGER, INTENT(IN)  :: im            
    
    ! original N-S dimension
    INTEGER, INTENT(IN)  :: jm            

    ! Target N-S dimension
    INTEGER, INTENT(IN)  :: jn           
    
    ! IG=0: scalars from SP to NP (D-grid v-wind is also IG=0)
    ! IG=1: D-grid u-wind
    INTEGER, INTENT(IN)  :: ig            
  
    ! IV=0: scalar; 
    ! IV=1: vector
    INTEGER, INTENT(IN)  :: iv            
  
    ! Original southern edge of the cell sin(lat1)  
    REAL*4,  INTENT(IN)  :: sin1(jm+1-ig) 
    
    ! Original data at center of the cell
    REAL*8,  INTENT(IN)  :: q1(im,jm)      
    
    ! Target cell's southern edge sin(lat2)
    REAL*4,  INTENT(IN)  :: sin2(jn+1-ig) 
!
! !OUTPUT PARAMETERS:
!
    ! Mapped data at the target resolution
    REAL*4,  INTENT(OUT) :: q2(im,jn)     
!
! !REMARKS:
!
!   sin1 (1) = -1 must be south pole; sin1(jm+1)=1 must be N pole.
!
!   sin1(1) < sin1(2) < sin1(3) < ... < sin1(jm) < sin1(jm+1)
!   sin2(1) < sin2(2) < sin2(3) < ... < sin2(jn) < sin2(jn+1)!
!
! !AUTHOR:
!   Developer: Prasad Kasibhatla
!   March 6, 2012
!
! !REVISION HISTORY
!  06 Mar 2012 - P. Kasibhatla - Initial version
!  27 Aug 2012 - R. Yantosca   - Added parallel DO loops
!  27 Aug 2012 - R. Yantosca   - Change REAL*4 variables to REAL(fp) to better
!                                ensure numerical stability
!  31 Mar 2014 - C. Keller     - Initialize qsum to zero to avoid undefined 
!                                values in nested grids
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: i, j0, m, mm, j
    REAL*8               :: dy1(jm)
    REAL*8               :: dy
    REAL*8               :: qsum, sum
    
    ! YMAP begins here!
    do j=1,jm-ig
       dy1(j) = sin1(j+1) - sin1(j)
    enddo

    !===============================================================
    ! Area preserving mapping
    !===============================================================
    
    !$OMP PARALLEL DO                          &
    !$OMP DEFAULT( SHARED                    ) &
    !$OMP PRIVATE( I, J0, J, M, QSUM, MM, DY )
    do 1000 i=1,im
       qsum = 0.0d0
       j0 = 1
       do 555 j=1,jn-ig
       do 100 m=j0,jm-ig
             
          !=========================================================
          ! locate the southern edge: sin2(i)
          !=========================================================
          if(sin2(j) .ge. sin1(m) .and. sin2(j) .le. sin1(m+1)) then
             
             if(sin2(j+1) .le. sin1(m+1)) then
                
                ! entire new cell is within the original cell
                q2(i,j)=q1(i,m)
                j0 = m
                goto 555
             else
                
                ! South most fractional area
                qsum=(sin1(m+1)-sin2(j))*q1(i,m)
                
                do mm=m+1,jm-ig
                   
                   ! locate the northern edge: sin2(j+1)
                   if(sin2(j+1) .gt. sin1(mm+1) ) then
                      
                      ! Whole layer
                      qsum = qsum + dy1(mm)*q1(i,mm)
                   else
                      
                      ! North most fractional area
                      dy = sin2(j+1)-sin1(mm)
                      qsum=qsum+dy*q1(i,mm)
                      j0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100    continue
123    q2(i,j) = qsum / ( sin2(j+1) - sin2(j) )
555    continue
1000 continue
     !$OMP END PARALLEL DO

     !===================================================================
     ! Final processing for poles
     !===================================================================
     if ( ig .eq. 0 .and. iv .eq. 0 ) then
         
        ! South pole
        sum = 0.0_f4
        do i=1,im
           sum = sum + q2(i,1)
        enddo

        sum = sum / DBLE( im )
        do i=1,im
           q2(i,1) = sum
        enddo

        ! North pole:
        sum = 0.0_f4
        do i=1,im
           sum = sum + q2(i,jn)
        enddo

        sum = sum / DBLE( im )
        do i=1,im
           q2(i,jn) = sum
        enddo

     endif

   END SUBROUTINE ymap_r8r4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ymap_r4r4
!
! !DESCRIPTION: Routine to perform area preserving mapping in N-S from an 
!  arbitrary resolution to another.  Both the input and output arguments
!  have REAL(fp) precision.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ymap_r4r4(im, jm, sin1, q1, jn, sin2, q2, ig, iv)
!
! !INPUT PARAMETERS:
!

    ! original E-W dimension
    INTEGER, INTENT(IN)  :: im            
    
    ! original N-S dimension
    INTEGER, INTENT(IN)  :: jm            

    ! Target N-S dimension
    INTEGER, INTENT(IN)  :: jn           
    
    ! IG=0: scalars from SP to NP (D-grid v-wind is also IG=0)
    ! IG=1: D-grid u-wind
    INTEGER, INTENT(IN)  :: ig            
  
    ! IV=0: scalar; 
    ! IV=1: vector
    INTEGER, INTENT(IN)  :: iv            
  
    ! Original southern edge of the cell sin(lat1)  
    REAL*4,  INTENT(IN)  :: sin1(jm+1-ig) 
    
    ! Original data at center of the cell
    REAL*4,  INTENT(IN)  :: q1(im,jm)      
    
    ! Target cell's southern edge sin(lat2)
    REAL*4,  INTENT(IN)  :: sin2(jn+1-ig) 
!
! !OUTPUT PARAMETERS:
!
    ! Mapped data at the target resolution
    REAL*4,  INTENT(OUT) :: q2(im,jn)     
!
! !REMARKS:
!
!   sin1 (1) = -1 must be south pole; sin1(jm+1)=1 must be N pole.
!
!   sin1(1) < sin1(2) < sin1(3) < ... < sin1(jm) < sin1(jm+1)
!   sin2(1) < sin2(2) < sin2(3) < ... < sin2(jn) < sin2(jn+1)!
!
! !AUTHOR:
!   Developer: Prasad Kasibhatla
!   March 6, 2012
!
! !REVISION HISTORY
!  06 Mar 2012 - P. Kasibhatla - Initial version
!  27 Aug 2012 - R. Yantosca   - Added parallel DO loops
!  27 Aug 2012 - R. Yantosca   - Change REAL*4 variables to REAL(fp) to better
!                                ensure numerical stability
!  31 Mar 2014 - C. Keller     - Initialize qsum to zero to avoid undefined 
!                                values in nested grids
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: i, j0, m, mm, j
    REAL*4               :: dy1(jm)
    REAL*4               :: dy
    REAL*4               :: qsum, sum
    
    ! YMAP begins here!
    do j=1,jm-ig
       dy1(j) = sin1(j+1) - sin1(j)
    enddo

    !===============================================================
    ! Area preserving mapping
    !===============================================================
    
    !$OMP PARALLEL DO                          &
    !$OMP DEFAULT( SHARED                    ) &
    !$OMP PRIVATE( I, J0, J, M, QSUM, MM, DY )
    do 1000 i=1,im
       qsum = 0.0
       j0 = 1
       do 555 j=1,jn-ig
       do 100 m=j0,jm-ig
             
          !=========================================================
          ! locate the southern edge: sin2(i)
          !=========================================================
          if(sin2(j) .ge. sin1(m) .and. sin2(j) .le. sin1(m+1)) then
             
             if(sin2(j+1) .le. sin1(m+1)) then
                
                ! entire new cell is within the original cell
                q2(i,j)=q1(i,m)
                j0 = m
                goto 555
             else
                
                ! South most fractional area
                qsum=(sin1(m+1)-sin2(j))*q1(i,m)
 
                do mm=m+1,jm-ig
                   
                   ! locate the northern edge: sin2(j+1)
                   if(sin2(j+1) .gt. sin1(mm+1) ) then
                      
                      ! Whole layer
                      qsum = qsum + dy1(mm)*q1(i,mm)
                   else
                      
                      ! North most fractional area
                      dy = sin2(j+1)-sin1(mm)
                      qsum=qsum+dy*q1(i,mm)
                      j0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100    continue
123    q2(i,j) = qsum / ( sin2(j+1) - sin2(j) )
555    continue
1000 continue
     !$OMP END PARALLEL DO

     !===================================================================
     ! Final processing for poles
     !===================================================================
     if ( ig .eq. 0 .and. iv .eq. 0 ) then
         
        ! South pole
        sum = 0.e+0_fp
        do i=1,im
           sum = sum + q2(i,1)
        enddo

        sum = sum / REAL( im, 4 )
        do i=1,im
           q2(i,1) = sum
        enddo

        ! North pole:
        sum = 0.e+0_fp
        do i=1,im
           sum = sum + q2(i,jn)
        enddo

        sum = sum / REAL( im, 4 )
        do i=1,im
           q2(i,jn) = sum
        enddo
     endif

   END SUBROUTINE ymap_r4r4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: xmap_r8r8
!
! !DESCRIPTION: Routine to perform area preserving mapping in E-W from an 
!  arbitrary resolution to another.  Both the input and output arguments
!  have REAL(fp) precision.
!\\
!\\
!  Periodic domain will be assumed, i.e., the eastern wall bounding cell
!  im is lon1(im+1) = lon1(1); Note the equal sign is true geographysically.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE xmap_r8r8(im, jm, lon1, q1, iin, ilon2, iq2)
!
! !INPUT PARAMETERS:
!
    ! Original E-W dimension
    INTEGER, INTENT(IN)  :: im           

    ! Target E-W dimension
    INTEGER, INTENT(IN)  :: iin
  
    ! Original N-S dimension
    INTEGER, INTENT(IN)  :: jm           
  
    ! Original western edge of the cell
    REAL*8,  INTENT(IN)  :: lon1(im+1)   
  
    ! Original data at center of the cell
    REAL*8,  INTENT(IN)  :: q1(im,jm)    
  
    ! Target cell's western edge
    REAL*8,  INTENT(IN), TARGET  :: ilon2(iin+1)   
!
! !OUTPUT PARAMETERS:
!
    ! Mapped data at the target resolution
    REAL*8,  INTENT(OUT), TARGET :: iq2(iin,jm)    
!
! !REMARKS:
!   lon1(1) < lon1(2) < lon1(3) < ... < lon1(im) < lon1(im+1)
!   lon2(1) < lon2(2) < lon2(3) < ... < lon2(in) < lon2(in+1)
!
! !AUTHOR:
!   Developer: Prasad Kasibhatla
!   March 6, 2012
!
! !REVISION HISTORY
!  06 Mar 2012 - P. Kasibhatla - Initial version
!  27 Aug 2012 - R. Yantosca   - Added parallel DO loops
!  27 Aug 2012 - R. Yantosca   - Change REAL*4 variables to REAL(fp) to better
!                                ensure numerical stability
!  15 May 2015 - C. Keller     - Now initialize qtmp to zero, and set q2 pointer
!                                to valid range n1:(n2-1). Do not initialize q2
!                                to zero after pointer assignment. This seems to
!                                cause problems with some compilers. 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: i1, i2, i, i0, m, mm, j
    REAL*8               :: qtmp(-im:im+im)
    REAL*8               :: x1(-im:im+im+1)
    REAL*8               :: dx1(-im:im+im)
    REAL*8               :: dx
    REAL*8               :: qsum
    LOGICAL              :: found

    ! Update
    INTEGER              :: n1, n2
    INTEGER              :: in
    REAL*8, POINTER      :: lon2(:) => NULL()
    REAL*8, POINTER      :: q2(:,:) => NULL()
    REAL*8               :: minlon, maxlon 
    REAL*8               :: lon1s(im+1)

    ! XMAP begins here!
    do i=1,im+1
       x1(i) = lon1(i)
    enddo
  
    do i=1,im
       dx1(i) = x1(i+1) - x1(i)
    enddo

    !===================================================================
    ! define minimum and maximum longitude on output grid
    ! to be used. Remapping will be restricted to this
    ! domain. This procedure allows remapping of nested
    ! domains onto larger (e.g. global) domains. 
    ! ckeller, 2/11/15).
    !===================================================================
    minlon = minval(lon1)
    maxlon = maxval(lon1)

    ! check for values > 180.0
    if(maxlon > 180.0d0) then
       lon1s = lon1
       do while(maxlon > 180.0d0)
          WHERE(lon1s > 180.0d0) lon1s = lon1s - 360.0d0
          minlon = minval(lon1s)
          maxlon = maxval(lon1s)
       enddo
    endif    

    ! maxlon must represent the easter edge of the grid:
    maxlon = maxlon + ( lon1(im+1)-lon1(im) )

    ! Reduce input grid
    n1 = 1
    n2 = iin+1
    do i=1,iin+1
       if ( ilon2(i) < minlon ) n1 = i
       if ( ilon2(iin+2-i) > maxlon ) n2 = iin+2-i
    enddo
    in = n2 - n1
    lon2 => ilon2(n1:n2)
    q2   => iq2(n1:(n2-1),:)
 
    !===================================================================
    ! check to see if ghosting is necessary
    ! Western edge:
    !===================================================================
    found = .false.
    i1 = 1
    do while ( .not. found )
       if( lon2(1) .ge. x1(i1) ) then
          found = .true.
       else
          i1 = i1 - 1
          if (i1 .lt. -im) then
             write(6,*) 'Failed in Xmap_R8R8 (regrid_a2a_mod.F90)'
             stop
          else
             x1(i1) = x1(i1+1) - dx1(im+i1)
             dx1(i1) = dx1(im+i1)
          endif
       endif
    enddo
    
    !===================================================================
    ! Eastern edge:
    !===================================================================
    found = .false.
    i2 = im+1
    do while ( .not. found )
       if( lon2(in+1) .le. x1(i2) ) then
          found = .true.
       else
          i2 = i2 + 1
          if (i2 .gt. 2*im) then
             write(6,*) 'Failed in Xmap_R8R8 (regrid_a2a_mod.F90)'
             stop
          else
             dx1(i2-1) = dx1(i2-1-im)
             x1(i2) = x1(i2-1) + dx1(i2-1)
          endif
       endif
    enddo

    !$OMP PARALLEL DO                                &
    !$OMP DEFAULT( SHARED                          ) &
    !$OMP PRIVATE( J, QTMP, I, I0, M, QSUM, MM, DX )
    do 1000 j=1,jm
       
       !=================================================================
       ! Area preserving mapping
       !================================================================
      
       qtmp(:) = 0.0d0 
       qtmp(0)=q1(im,j)
       do i=1,im
          qtmp(i)=q1(i,j)
       enddo
       qtmp(im+1)=q1(1,j)

       ! check to see if ghosting is necessary
       ! Western edge
       if ( i1 .le. 0 ) then
          do i=i1,0
             qtmp(i) = qtmp(im+i)
          enddo
       endif
       
       ! Eastern edge:
       if ( i2 .gt. im+1 ) then
          do i=im+1,i2-1
             qtmp(i) = qtmp(i-im)
          enddo
       endif
        
       i0 = i1

       do 555 i=1,in
       do 100 m=i0,i2-1

          !=============================================================  
          ! locate the western edge: lon2(i)
          !=============================================================  
          if(lon2(i) .ge. x1(m) .and. lon2(i) .le. x1(m+1)) then
             
             if(lon2(i+1) .le. x1(m+1)) then
                
                ! entire new grid is within the original grid
                q2(i,j)=qtmp(m)
                i0 = m
                goto 555
             else
  
                ! Left most fractional area
                qsum=(x1(m+1)-lon2(i))*qtmp(m)
                do mm=m+1,i2-1
                   
                   ! locate the eastern edge: lon2(i+1)
                   if(lon2(i+1) .gt. x1(mm+1) ) then
                      
                      ! Whole layer
                      qsum = qsum + dx1(mm)*qtmp(mm)
                      
                   else
                      ! Right most fractional area
                      dx = lon2(i+1)-x1(mm)
                      qsum=qsum+dx*qtmp(mm)
                      i0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100    continue
123    q2(i,j) = qsum / ( lon2(i+1) - lon2(i) )
555    continue
1000 continue
     !$OMP END PARALLEL DO

    ! Cleanup
    lon2 => NULL() 
    q2   => NULL() 
    
  END SUBROUTINE xmap_r8r8
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: xmap_r4r4
!
! !DESCRIPTION: Routine to perform area preserving mapping in E-W from an 
!  arbitrary resolution to another.  Both the input and output arguments
!  have REAL*4 precision.
!\\
!\\
!  Periodic domain will be assumed, i.e., the eastern wall bounding cell
!  im is lon1(im+1) = lon1(1); Note the equal sign is true geographysically.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE xmap_r4r4(im, jm, lon1, q1, iin, ilon2, iq2)
!
! !INPUT PARAMETERS:
!
    ! Original E-W dimension
    INTEGER, INTENT(IN)  :: im           

    ! Target E-W dimension
    INTEGER, INTENT(IN)  :: iin
  
    ! Original N-S dimension
    INTEGER, INTENT(IN)  :: jm           
  
    ! Original western edge of the cell
    REAL*4,  INTENT(IN)  :: lon1(im+1)   
  
    ! Original data at center of the cell
    REAL*4,  INTENT(IN)  :: q1(im,jm)    
  
    ! Target cell's western edge
    REAL*4,  INTENT(IN), TARGET  :: ilon2(iin+1)   
!
! !OUTPUT PARAMETERS:
!
    ! Mapped data at the target resolution
    REAL*4,  INTENT(OUT), TARGET :: iq2(iin,jm)    
!
! !REMARKS:
!   lon1(1) < lon1(2) < lon1(3) < ... < lon1(im) < lon1(im+1)
!   lon2(1) < lon2(2) < lon2(3) < ... < lon2(in) < lon2(in+1)
!
! !AUTHOR:
!   Developer: Prasad Kasibhatla
!   March 6, 2012
!
! !REVISION HISTORY
!  06 Mar 2012 - P. Kasibhatla - Initial version
!  27 Aug 2012 - R. Yantosca   - Added parallel DO loops
!  27 Aug 2012 - R. Yantosca   - Change REAL*4 variables to REAL(fp) to better
!                                ensure numerical stability
!  15 May 2015 - C. Keller     - Now initialize qtmp to zero, and set q2 pointer
!                                to valid range n1:(n2-1). Do not initialize q2
!                                to zero after pointer assignment. This seems to
!                                cause problems with some compilers. 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: i1, i2, i, i0, m, mm, j
    REAL*4               :: qtmp(-im:im+im)
    REAL*4               :: x1(-im:im+im+1)
    REAL*4               :: dx1(-im:im+im)
    REAL*4               :: dx
    REAL*4               :: qsum
    LOGICAL              :: found

    ! Update
    INTEGER              :: n1, n2
    INTEGER              :: in
    REAL*4, POINTER      :: lon2(:) => NULL()
    REAL*4, POINTER      :: q2(:,:) => NULL()
    REAL*4               :: minlon, maxlon 
    REAL*4               :: lon1s(im+1)

    ! XMAP begins here!
    do i=1,im+1
       x1(i) = lon1(i)
    enddo
  
    do i=1,im
       dx1(i) = x1(i+1) - x1(i)
    enddo

    !===================================================================
    ! define minimum and maximum longitude on output grid
    ! to be used. Remapping will be restricted to this
    ! domain. This procedure allows remapping of nested
    ! domains onto larger (e.g. global) domains. 
    ! ckeller, 2/11/15).
    !===================================================================
    minlon = minval(lon1)
    maxlon = maxval(lon1)

    ! check for values > 180.0
    if(maxlon > 180.0) then
       lon1s = lon1
       do while(maxlon > 180.0)
          WHERE(lon1s > 180.0) lon1s = lon1s - 360.0
          minlon = minval(lon1s)
          maxlon = maxval(lon1s)
       enddo
    endif    

    ! maxlon must represent the easter edge of the grid:
    maxlon = maxlon + ( lon1(im+1)-lon1(im) )

    ! Reduce output grid
    n1 = 1
    n2 = iin+1
    do i=1,iin+1
       if ( ilon2(i) < minlon ) n1 = i
       if ( ilon2(iin+2-i) > maxlon ) n2 = iin+2-i
    enddo
    in = n2 - n1
    lon2 => ilon2(n1:n2)
    q2   => iq2(n1:(n2-1),:)

    ! shadow variables to selected range
 
    !===================================================================
    ! check to see if ghosting is necessary
    ! Western edge:
    !===================================================================
    found = .false.
    i1 = 1
    do while ( .not. found )
       if( lon2(1) .ge. x1(i1) ) then
          found = .true.
       else
          i1 = i1 - 1
          if (i1 .lt. -im) then
             write(6,*) 'Failed in Xmap_R4R4 (regrid_a2a_mod.F90)'
             stop
          else
             x1(i1) = x1(i1+1) - dx1(im+i1)
             dx1(i1) = dx1(im+i1)
          endif
       endif
    enddo
    
    !===================================================================
    ! Eastern edge:
    !===================================================================
    found = .false.
    i2 = im+1
    do while ( .not. found )
       if( lon2(in+1) .le. x1(i2) ) then
          found = .true.
       else
          i2 = i2 + 1
          if (i2 .gt. 2*im) then
             write(6,*) 'Failed in Xmap_R4R4 (regrid_a2a_mod.F90)'
             stop
          else
             dx1(i2-1) = dx1(i2-1-im)
             x1(i2) = x1(i2-1) + dx1(i2-1)
          endif
       endif
    enddo

    !$OMP PARALLEL DO                                &
    !$OMP DEFAULT( SHARED                          ) &
    !$OMP PRIVATE( J, QTMP, I, I0, M, QSUM, MM, DX )
    do 1000 j=1,jm
       
       !=================================================================
       ! Area preserving mapping
       !================================================================
       
       qtmp(:) = 0.0
       qtmp(0)=q1(im,j)
       do i=1,im
          qtmp(i)=q1(i,j)
       enddo
       qtmp(im+1)=q1(1,j)

       ! check to see if ghosting is necessary
       ! Western edge
       if ( i1 .le. 0 ) then
          do i=i1,0
             qtmp(i) = qtmp(im+i)
          enddo
       endif
       
       ! Eastern edge:
       if ( i2 .gt. im+1 ) then
          do i=im+1,i2-1
             qtmp(i) = qtmp(i-im)
          enddo
       endif
        
       i0 = i1

       do 555 i=1,in
       do 100 m=i0,i2-1

          !=============================================================  
          ! locate the western edge: lon2(i)
          !=============================================================  
          if(lon2(i) .ge. x1(m) .and. lon2(i) .le. x1(m+1)) then
             
             if(lon2(i+1) .le. x1(m+1)) then
                
                ! entire new grid is within the original grid
                q2(i,j)=qtmp(m)
                i0 = m
                goto 555
             else
  
                ! Left most fractional area
                qsum=(x1(m+1)-lon2(i))*qtmp(m)
                do mm=m+1,i2-1
                   
                   ! locate the eastern edge: lon2(i+1)
                   if(lon2(i+1) .gt. x1(mm+1) ) then
                      
                      ! Whole layer
                      qsum = qsum + dx1(mm)*qtmp(mm)
                      
                   else
                      ! Right most fractional area
                      dx = lon2(i+1)-x1(mm)
                      qsum=qsum+dx*qtmp(mm)
                      i0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100    continue
123    q2(i,j) = qsum / ( lon2(i+1) - lon2(i) )
555    continue
1000 continue
     !$OMP END PARALLEL DO

     ! Cleanup
     lon2 => NULL()
     q2   => NULL()

  END SUBROUTINE xmap_r4r4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: xmap_r4r8
!
! !DESCRIPTION: Routine to perform area preserving mapping in E-W from an 
!  arbitrary resolution to another.  The input argument has REAL*4 precision
!  but the output argument has REAL(fp) precision.
!\\
!\\
!  Periodic domain will be assumed, i.e., the eastern wall bounding cell
!  im is lon1(im+1) = lon1(1); Note the equal sign is true geographysically.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE xmap_r4r8(im, jm, lon1, q1, iin, ilon2, iq2)
!
! !INPUT PARAMETERS:
!
    ! Original E-W dimension
    INTEGER, INTENT(IN)  :: im           

    ! Target E-W dimension
    INTEGER, INTENT(IN)  :: iin
  
    ! Original N-S dimension
    INTEGER, INTENT(IN)  :: jm           
  
    ! Original western edge of the cell
    REAL*4,  INTENT(IN)  :: lon1(im+1)   
  
    ! Original data at center of the cell
    REAL*4,  INTENT(IN)  :: q1(im,jm)    
  
    ! Target cell's western edge
    REAL*4,  INTENT(IN), TARGET  :: ilon2(iin+1)   
!
! !OUTPUT PARAMETERS:
!
    ! Mapped data at the target resolution
    REAL*8,  INTENT(OUT), TARGET :: iq2(iin,jm)    
!
! !REMARKS:
!   lon1(1) < lon1(2) < lon1(3) < ... < lon1(im) < lon1(im+1)
!   lon2(1) < lon2(2) < lon2(3) < ... < lon2(in) < lon2(in+1)
!
! !AUTHOR:
!   Developer: Prasad Kasibhatla
!   March 6, 2012
!
! !REVISION HISTORY
!  06 Mar 2012 - P. Kasibhatla - Initial version
!  27 Aug 2012 - R. Yantosca   - Added parallel DO loops
!  27 Aug 2012 - R. Yantosca   - Change REAL*4 variables to REAL(fp) to better
!                                ensure numerical stability
!  15 May 2015 - C. Keller     - Now initialize qtmp to zero, and set q2 pointer
!                                to valid range n1:(n2-1). Do not initialize q2
!                                to zero after pointer assignment. This seems to
!                                cause problems with some compilers. 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: i1, i2, i, i0, m, mm, j
    REAL*8               :: qtmp(-im:im+im)
    REAL*8               :: x1(-im:im+im+1)
    REAL*8               :: dx1(-im:im+im)
    REAL*8               :: dx
    REAL*8               :: qsum
    LOGICAL              :: found

    ! Update
    INTEGER              :: n1, n2
    INTEGER              :: in
    REAL*4, POINTER      :: lon2(:) => NULL()
    REAL*8, POINTER      :: q2(:,:) => NULL()
    REAL*4               :: minlon, maxlon 
    REAL*4               :: lon1s(im+1)

    ! XMAP begins here!
    do i=1,im+1
       x1(i) = lon1(i)
    enddo
  
    do i=1,im
       dx1(i) = x1(i+1) - x1(i)
    enddo
    
    !===================================================================
    ! define minimum and maximum longitude on output grid
    ! to be used. Remapping will be restricted to this
    ! domain. This procedure allows remapping of nested
    ! domains onto larger (e.g. global) domains. 
    ! ckeller, 2/11/15).
    !===================================================================
    minlon = minval(lon1)
    maxlon = maxval(lon1)

    ! check for values > 180.0
    if(maxlon > 180.0) then
       lon1s = lon1
       do while(maxlon > 180.0)
          WHERE(lon1s > 180.0) lon1s = lon1s - 360.0
          minlon = minval(lon1s)
          maxlon = maxval(lon1s)
       enddo
    endif    

    ! maxlon must represent the easter edge of the grid:
    maxlon = maxlon + ( lon1(im+1)-lon1(im) )

    ! Reduce input grid
    n1 = 1
    n2 = iin+1
    do i=1,iin+1
       if ( ilon2(i) < minlon ) n1 = i
       if ( ilon2(iin+2-i) > maxlon ) n2 = iin+2-i
    enddo
    in = n2 - n1
    lon2 => ilon2(n1:n2)
    q2   => iq2(n1:(n2-1),:)
    
    !===================================================================
    ! check to see if ghosting is necessary
    ! Western edge:
    !===================================================================
    found = .false.
    i1 = 1
    do while ( .not. found )
       if( lon2(1) .ge. x1(i1) ) then
          found = .true.
       else
          i1 = i1 - 1
          if (i1 .lt. -im) then
             write(6,*) 'Failed in Xmap_R4R8 (regrid_a2a_mod.F90)'
             stop
          else
             x1(i1) = x1(i1+1) - dx1(im+i1)
             dx1(i1) = dx1(im+i1)
          endif
       endif
    enddo
    
    !===================================================================
    ! Eastern edge:
    !===================================================================
    found = .false.
    i2 = im+1
    do while ( .not. found )
       if( lon2(in+1) .le. x1(i2) ) then
          found = .true.
       else
          i2 = i2 + 1
          if (i2 .gt. 2*im) then
             write(6,*) 'Failed in Xmap_R4R8 (regrid_a2a_mod.F90)'
             stop
          else
             dx1(i2-1) = dx1(i2-1-im)
             x1(i2) = x1(i2-1) + dx1(i2-1)
          endif
       endif
    enddo

    !$OMP PARALLEL DO                                &
    !$OMP DEFAULT( SHARED                          ) &
    !$OMP PRIVATE( J, QTMP, I, I0, M, QSUM, MM, DX )
    do 1000 j=1,jm
       
       !=================================================================
       ! Area preserving mapping
       !================================================================
       
       qtmp(:) = 0.0d0
       qtmp(0)=q1(im,j)
       do i=1,im
          qtmp(i)=q1(i,j)
       enddo
       qtmp(im+1)=q1(1,j)

       ! check to see if ghosting is necessary
       ! Western edge
       if ( i1 .le. 0 ) then
          do i=i1,0
             qtmp(i) = qtmp(im+i)
          enddo
       endif
       
       ! Eastern edge:
       if ( i2 .gt. im+1 ) then
          do i=im+1,i2-1
             qtmp(i) = qtmp(i-im)
          enddo
       endif
        
       i0 = i1

       do 555 i=1,in
       do 100 m=i0,i2-1

          !=============================================================  
          ! locate the western edge: lon2(i)
          !=============================================================  
          if(lon2(i) .ge. x1(m) .and. lon2(i) .le. x1(m+1)) then
             
             if(lon2(i+1) .le. x1(m+1)) then
                
                ! entire new grid is within the original grid
                q2(i,j)=qtmp(m)
                i0 = m
                goto 555
             else
  
                ! Left most fractional area
                qsum=(x1(m+1)-lon2(i))*qtmp(m)
                do mm=m+1,i2-1
                   
                   ! locate the eastern edge: lon2(i+1)
                   if(lon2(i+1) .gt. x1(mm+1) ) then
                      
                      ! Whole layer
                      qsum = qsum + dx1(mm)*qtmp(mm)
                      
                   else
                      ! Right most fractional area
                      dx = lon2(i+1)-x1(mm)
                      qsum=qsum+dx*qtmp(mm)
                      i0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100    continue
123    q2(i,j) = qsum / ( lon2(i+1) - lon2(i) )
555    continue
1000 continue
     !$OMP END PARALLEL DO

    ! Cleanup
    lon2 => NULL() 
    q2   => NULL() 
    
  END SUBROUTINE xmap_r4r8
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: xmap_r8r4
!
! !DESCRIPTION: Routine to perform area preserving mapping in E-W from an 
!  arbitrary resolution to another.  The input argument has REAL*8 precision
!  but the output argument has REAL*4 precision.
!\\
!\\
!  Periodic domain will be assumed, i.e., the eastern wall bounding cell
!  im is lon1(im+1) = lon1(1); Note the equal sign is true geographysically.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE xmap_r8r4(im, jm, lon1, q1, iin, ilon2, iq2)
!
! !INPUT PARAMETERS:
!
    ! Original E-W dimension
    INTEGER, INTENT(IN)  :: im           

    ! Target E-W dimension
    INTEGER, INTENT(IN)  :: iin           
  
    ! Original N-S dimension
    INTEGER, INTENT(IN)  :: jm           
  
    ! Original western edge of the cell
    REAL*4,  INTENT(IN)  :: lon1(im+1)   
  
    ! Original data at center of the cell
    REAL*8,  INTENT(IN)  :: q1(im,jm)    
  
    ! Target cell's western edge
    REAL*4,  INTENT(IN), TARGET  :: ilon2(iin+1)   
!
! !OUTPUT PARAMETERS:
!
    ! Mapped data at the target resolution
    REAL*4,  INTENT(OUT), TARGET :: iq2(iin,jm)    
!
! !REMARKS:
!   lon1(1) < lon1(2) < lon1(3) < ... < lon1(im) < lon1(im+1)
!   lon2(1) < lon2(2) < lon2(3) < ... < lon2(in) < lon2(in+1)
!
! !AUTHOR:
!   Developer: Prasad Kasibhatla
!   March 6, 2012
!
! !REVISION HISTORY
!  06 Mar 2012 - P. Kasibhatla - Initial version
!  27 Aug 2012 - R. Yantosca   - Added parallel DO loops
!  27 Aug 2012 - R. Yantosca   - Change REAL*4 variables to REAL(fp) to better
!                                ensure numerical stability
!  15 May 2015 - C. Keller     - Now initialize qtmp to zero, and set q2 pointer
!                                to valid range n1:(n2-1). Do not initialize q2
!                                to zero after pointer assignment. This seems to
!                                cause problems with some compilers. 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: i1, i2, i, i0, m, mm, j
    REAL*4               :: qtmp(-im:im+im)
    REAL*4               :: x1(-im:im+im+1)
    REAL*4               :: dx1(-im:im+im)
    REAL*4               :: dx
    REAL*4               :: qsum
    LOGICAL              :: found

    ! Update
    INTEGER              :: n1, n2
    INTEGER              :: in
    REAL*4, POINTER      :: lon2(:) => NULL()
    REAL*4, POINTER      :: q2(:,:) => NULL()
    REAL*4               :: minlon, maxlon 
    REAL*4               :: lon1s(im+1)

    ! XMAP begins here!
    do i=1,im+1
       x1(i) = lon1(i)
    enddo
  
    do i=1,im
       dx1(i) = x1(i+1) - x1(i)
    enddo
    
    !===================================================================
    ! define minimum and maximum longitude on output grid
    ! to be used. Remapping will be restricted to this
    ! domain. This procedure allows remapping of nested
    ! domains onto larger (e.g. global) domains. 
    ! ckeller, 2/11/15).
    !===================================================================
    minlon = minval(lon1)
    maxlon = maxval(lon1)

    ! check for values > 180.0
    if(maxlon > 180.0) then
       lon1s = lon1
       do while(maxlon > 180.0)
          WHERE(lon1s > 180.0) lon1s = lon1s - 360.0
          minlon = minval(lon1s)
          maxlon = maxval(lon1s)
       enddo
    endif    

    ! maxlon must represent the easter edge of the grid:
    maxlon = maxlon + ( lon1(im+1)-lon1(im) )

    ! Reduce input grid
    n1 = 1
    n2 = iin+1
    do i=1,iin+1
       if ( ilon2(i) < minlon ) n1 = i
       if ( ilon2(iin+2-i) > maxlon ) n2 = iin+2-i
    enddo
    in = n2 - n1
    lon2 => ilon2(n1:n2)
    q2   => iq2(n1:(n2-1),:)
    
    !===================================================================
    ! check to see if ghosting is necessary
    ! Western edge:
    !===================================================================
    found = .false.
    i1 = 1
    do while ( .not. found )
       if( lon2(1) .ge. x1(i1) ) then
          found = .true.
       else
          i1 = i1 - 1
          if (i1 .lt. -im) then
             write(6,*) 'Failed in Xmap_R4R8 (regrid_a2a_mod.F90)'
             stop
          else
             x1(i1) = x1(i1+1) - dx1(im+i1)
             dx1(i1) = dx1(im+i1)
          endif
       endif
    enddo
    
    !===================================================================
    ! Eastern edge:
    !===================================================================
    found = .false.
    i2 = im+1
    do while ( .not. found )
       if( lon2(in+1) .le. x1(i2) ) then
          found = .true.
       else
          i2 = i2 + 1
          if (i2 .gt. 2*im) then
             write(6,*) 'Failed in Xmap_R4R8 (regrid_a2a_mod.F90)'
             stop
          else
             dx1(i2-1) = dx1(i2-1-im)
             x1(i2) = x1(i2-1) + dx1(i2-1)
          endif
       endif
    enddo

    !$OMP PARALLEL DO                                &
    !$OMP DEFAULT( SHARED                          ) &
    !$OMP PRIVATE( J, QTMP, I, I0, M, QSUM, MM, DX )
    do 1000 j=1,jm
       
       !=================================================================
       ! Area preserving mapping
       !================================================================
      
       qtmp(:) = 0.0
       qtmp(0)=q1(im,j)
       do i=1,im
          qtmp(i)=q1(i,j)
       enddo
       qtmp(im+1)=q1(1,j)

       ! check to see if ghosting is necessary
       ! Western edge
       if ( i1 .le. 0 ) then
          do i=i1,0
             qtmp(i) = qtmp(im+i)
          enddo
       endif
       
       ! Eastern edge:
       if ( i2 .gt. im+1 ) then
          do i=im+1,i2-1
             qtmp(i) = qtmp(i-im)
          enddo
       endif
        
       i0 = i1

       do 555 i=1,in
       do 100 m=i0,i2-1

          !=============================================================  
          ! locate the western edge: lon2(i)
          !=============================================================  
          if(lon2(i) .ge. x1(m) .and. lon2(i) .le. x1(m+1)) then
             
             if(lon2(i+1) .le. x1(m+1)) then
                
                ! entire new grid is within the original grid
                q2(i,j)=qtmp(m)
                i0 = m
                goto 555
             else
  
                ! Left most fractional area
                qsum=(x1(m+1)-lon2(i))*qtmp(m)
                do mm=m+1,i2-1
                   
                   ! locate the eastern edge: lon2(i+1)
                   if(lon2(i+1) .gt. x1(mm+1) ) then
                      
                      ! Whole layer
                      qsum = qsum + dx1(mm)*qtmp(mm)
                      
                   else
                      ! Right most fractional area
                      dx = lon2(i+1)-x1(mm)
                      qsum=qsum+dx*qtmp(mm)
                      i0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100    continue
123    q2(i,j) = qsum / ( lon2(i+1) - lon2(i) )
555    continue
1000 continue
     !$OMP END PARALLEL DO

    ! Cleanup
    lon2 => NULL() 
    q2   => NULL() 
    
  END SUBROUTINE xmap_r8r4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_input_grid
!
! !DESCRIPTION: Routine to read variables and attributes from a netCDF
!  file.  This routine was automatically generated by the Perl script
!  NcdfUtilities/perl/ncCodeRead.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Read_Input_Grid( IM, JM, fileName, lon_edges, lat_sines )
!
! !USES:
!
    ! Modules for netCDF read
    USE m_netcdf_io_open
    USE m_netcdf_io_get_dimlen
    USE m_netcdf_io_read
    USE m_netcdf_io_readattr
    USE m_netcdf_io_close

    IMPLICIT NONE

#   include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: IM                ! # of longitudes
    INTEGER,          INTENT(IN)  :: JM                ! # of latitudes
    CHARACTER(LEN=*), INTENT(IN)  :: fileName          ! File w/ grid info
!
! !OUTPUT PARAMETERS:
!   
    REAL(fp),           INTENT(OUT) :: lon_edges(IM+1)   ! Lon edges [degrees]
    REAL(fp),           INTENT(OUT) :: lat_sines(JM+1)   ! SIN( latitude edges )
!
! !REMARKS:
!  Created with the ncCodeRead script of the NcdfUtilities package,
!  with subsequent hand-editing.
!
! !REVISION HISTORY:
!  23 Aug 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: fId                          ! netCDF file ID

    ! Arrays
    INTEGER            :: st1d(1), ct1d(1)             ! netCDF start & count

    !======================================================================
    ! Read data from file
    !======================================================================

    ! Open file for reading
    CALL Ncop_Rd( fId, TRIM( fileName ) )

    ! Read lon_edges from file
    st1d = (/ 1    /)
    ct1d = (/ IM+1 /)
    CALL NcRd( lon_edges, fId,  "lon_edges", st1d, ct1d )
        
    ! Read lat_sines from file
    st1d = (/ 1    /)
    ct1d = (/ JM+1 /)
    CALL NcRd( lat_sines, fId,  "lat_sines", st1d, ct1d )

    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE Read_Input_Grid
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Map_A2A
!
! !DESCRIPTION: Subroutine Init\_Map\_A2A initializes all module variables.
!  This allows us to keep "shadow" copies of variables that are defined
!  elsewhere in GEOS-Chem.  This also helps us from having dependencies to
!  GEOS-Chem modules in the HEMCO core modules.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Map_A2A( NX, NY, LONS, SINES, AREAS, DIR )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: NX             ! # of longitudes
    INTEGER,          INTENT(IN) :: NY             ! # of latitudes
    REAL(fp),           INTENT(IN) :: LONS (NX+1 )   ! Longitudes
    REAL(fp),           INTENT(IN) :: SINES(NY+1 )   ! Sines of latitudes
    REAL(fp),           INTENT(IN) :: AREAS(NX,NY)   ! Surface areas [m2]
    CHARACTER(LEN=*), INTENT(IN) :: DIR            ! Dir for netCDF files w/ 
                                                   !  grid definitions
!
! !REVISION HISTORY:
!  14 Jul 2014 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: AS

    !------------------------------------------
    ! Allocate module variables
    !------------------------------------------
    IF ( .not. ALLOCATED( OUTLON ) ) THEN
       ALLOCATE( OUTLON( NX+1 ), STAT=AS )
       IF ( AS /= 0 ) THEN 
          PRINT*, '### Could not allocate OUTLON (regrid_a2a_mod.F90)'
          STOP
       ENDIF
    ENDIF

    IF ( .not. ALLOCATED( OUTSIN ) ) THEN
       ALLOCATE( OUTSIN( NY+1 ), STAT=AS ) 
       IF ( AS /= 0 ) THEN
          PRINT*, '### Could not allocate OUTSIN (regrid_a2a_mod.F90)'
          STOP
       ENDIF
    ENDIF

    IF ( .not. ALLOCATED( OUTAREA ) ) THEN 
       ALLOCATE( OUTAREA( NX, NY ), STAT=AS ) 
       IF ( AS /= 0 ) THEN
          PRINT*, '### Could not allocate OUTAREA (regrid_a2a_mod.F90)'
          STOP
       ENDIF
    ENDIF

    !------------------------------------------
    ! Store values in local shadow variables
    !------------------------------------------
    IIPAR   = NX
    JJPAR   = NY
    OUTLON  = LONS
    OUTSIN  = SINES
    OUTAREA = AREAS
    NC_DIR  = DIR

  END SUBROUTINE Init_Map_A2A
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Map_A2A
!
! !DESCRIPTION: Subroutine Cleanup\_Map\_A2A deallocates all module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Map_A2A()
!
! !REVISION HISTORY:
!  14 Jul 2014 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Cleanup module variables
    IF ( ALLOCATED( OUTLON  ) ) DEALLOCATE( OUTLON  )
    IF ( ALLOCATED( OUTSIN  ) ) DEALLOCATE( OUTSIN  )
    IF ( ALLOCATED( OUTAREA ) ) DEALLOCATE( OUTAREA )

  END SUBROUTINE Cleanup_Map_A2A
!EOC
END MODULE Regrid_A2A_Mod
