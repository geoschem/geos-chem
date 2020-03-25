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
!  See https://github.com/geoschem/geos-chem for complete history
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
  CHARACTER(LEN=255)    :: NC_DIR        ! Directory w/ netCDF files
  INTEGER               :: OUTNX         ! # of longitudes (x-dimension) in grid
  INTEGER               :: OUTNY         ! # of latitudes  (y-dimension) in grid
  REAL(fp), ALLOCATABLE :: OUTLON (:  )  ! Longitude on output grid
  REAL(fp), ALLOCATABLE :: OUTSIN (:  )  ! Sines of latitudes on output grid
  REAL(fp), ALLOCATABLE :: OUTAREA(:,:)  ! Surface areas on output grid
!
! !DEFINED PARAMETERS:
!
  !---------------------------------------------------------------------------
  ! These were taken from CMN_GCTM_mod.F90.  This helps us to avoid depending
  ! on GEOS-Chem modules in the core HEMCO modules.  (bmy, 7/14/14)
  ! NOTE: CMN_GCTM_mod.F90 is now physconstants.F90 (ewl, 1/8/2016)
  !---------------------------------------------------------------------------
  REAL(fp), PARAMETER :: PI =   3.14159265358979323e+0_fp   ! Pi
  REAL(fp), PARAMETER :: Re =   6.375d6                 ! Earth radius [m]

  !---------------------------------------------------------------------------
  ! Tiny numbers for single and double precision. These are being used for
  ! skipping missing values. miss_r4 and miss_r8 are the default missing values
  ! for single and double precision, respectively. (ckeller, 4/8/2017)
  !---------------------------------------------------------------------------
  REAL*4, PARAMETER   :: tiny_r4 = 1.0e-20
  REAL*4, PARAMETER   :: miss_r4 = 0.0e0
  REAL*8, PARAMETER   :: tiny_r8 = 1.0d-40
  REAL*8, PARAMETER   :: miss_r8 = 0.0d0

CONTAINS
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
!  See https://github.com/geoschem/geos-chem for complete history
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
       qtmp = miss_r8
       q2   = miss_r8
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
       CALL xmap_r8r8(im, jm-ig, lon1, q1(1,1+ig),in, lon2, qtmp(1,1+ig), &
                      missval=missval )

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
       CALL ymap_r8r8(in, jm, sin1, qtmp(1,1+ig), jn, sin2, q2(1,1+ig), ig, iv, &
                      missval=missval )

    ENDIF

  END SUBROUTINE Map_A2A_r8r8
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Map_A2A_r4r4
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
!  See https://github.com/geoschem/geos-chem for complete history
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
       qtmp = miss_r4
       q2   = miss_r4
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
       CALL xmap_r4r4(im, jm-ig, lon1, q1(1,1+ig),in, lon2, qtmp(1,1+ig), &
                      missval=missval )

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
       CALL ymap_r4r4(in, jm, sin1, qtmp(1,1+ig), jn, sin2, q2(1,1+ig), ig, iv, &
                      missval=missval)

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
    REAL*4,  INTENT(IN), OPTIONAL :: missval
!
! !REMARKS:
!  This routine is overloaded by the MAP_A2A interface.
!
! !REVISION HISTORY:
!  (1) Original subroutine by S-J Lin.  Converted to F90 freeform format
!      and inserted into "Geos3RegridModule" by Bob Yantosca (9/21/00)
!  See https://github.com/geoschem/geos-chem for complete history
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
       qtmp = real(missval,8)
       q2   = real(missval,8)
    ELSE
       qtmp = miss_r8
       q2   = miss_r8
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
       CALL xmap_r4r8(im, jm-ig, lon1, q1(1,1+ig),in, lon2, qtmp(1,1+ig), &
                      missval=missval )

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
       CALL ymap_r4r8(in, jm, sin1, qtmp(1,1+ig), jn, sin2, q2(1,1+ig), ig, iv, &
                      missval=missval )

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
    REAL*8,  INTENT(IN), OPTIONAL :: missval
!
! !REMARKS:
!  This routine is overloaded by the MAP_A2A interface.
!
! !REVISION HISTORY:
!  (1) Original subroutine by S-J Lin.  Converted to F90 freeform format
!      and inserted into "Geos3RegridModule" by Bob Yantosca (9/21/00)
!  See https://github.com/geoschem/geos-chem for complete history
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
       qtmp = real(missval,4)
       q2   = real(missval,4)
    ELSE
       qtmp = miss_r4
       q2   = miss_r4
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
       CALL xmap_r8r4(im, jm-ig, lon1, q1(1,1+ig),in, lon2, qtmp(1,1+ig), &
                      missval=missval )

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
       CALL ymap_r4r4(in, jm, sin1, qtmp(1,1+ig), jn, sin2, q2(1,1+ig), ig, iv, &
                      missval=real(missval,4))

    ENDIF

  END SUBROUTINE Map_A2A_r8r4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ymap_r8r8
!
! !DESCRIPTION: Routine to perform area preserving mapping in N-S from an
!  arbitrary resolution to another.  Both the input and output arguments
!  have REAL(fp) precision.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ymap_r8r8(im, jm, sin1, q1, jn, sin2, q2, ig, iv, missval )
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
! !OPTIONAL INPUT PARAMETERS:
!
    REAL*8,  INTENT(IN), OPTIONAL :: missval
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    REAL*8               :: dlat, nlon, miss

    ! YMAP begins here!
    do j=1,jm-ig
       dy1(j) = sin1(j+1) - sin1(j)
    enddo

    ! Missing value
    miss = miss_r8
    if ( present(missval) ) miss=missval

    !===============================================================
    ! Area preserving mapping
    !===============================================================

    !$OMP PARALLEL DO                                &
    !$OMP DEFAULT( SHARED                          ) &
    !$OMP PRIVATE( I, J0, J, M, QSUM, DLAT, MM, DY )
    do 1000 i=1,im
       qsum = 0.0d0
       dlat = 0.0d0
       j0 = 1
       do 555 j=1,jn-ig
       do 100 m=j0,jm-ig

          !=========================================================
          ! locate the southern edge: sin2(i)
          !=========================================================
          if(sin2(j) .ge. sin1(m) .and. sin2(j) .le. sin1(m+1)) then

             if(sin2(j+1) .le. sin1(m+1)) then

                ! entire new cell is within the original cell
                if( abs(q1(i,m)-miss)>tiny_r8 ) q2(i,j)=q1(i,m)
                j0 = m
                goto 555
             else

                ! South most fractional area
                if( abs(q1(i,m)-miss)>tiny_r8 ) then
                   dlat= sin1(m+1)-sin2(j)
                   qsum=(sin1(m+1)-sin2(j))*q1(i,m)
                endif

                do mm=m+1,jm-ig

                   ! locate the northern edge: sin2(j+1)
                   if(sin2(j+1) .gt. sin1(mm+1) ) then

                      ! Whole layer
                      if( abs(q1(i,mm)-miss)>tiny_r8 ) then
                         dlat = dlat + dy1(mm)
                         qsum = qsum + dy1(mm)*q1(i,mm)
                      endif
                   else

                      ! North most fractional area
                      dy = sin2(j+1)-sin1(mm)
                      if ( abs(q1(i,mm)-miss)>tiny_r8 ) then
                         qsum=qsum+dy*q1(i,mm)
                         dlat=dlat+dy
                      endif
                      j0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100    continue
!123    q2(i,j) = qsum / ( sin2(j+1) - sin2(j) )
123    if ( dlat /= 0.0d0 ) q2(i,j) = qsum / dlat
555    continue
1000 continue
     !$OMP END PARALLEL DO

     !===================================================================
     ! Final processing for poles
     !===================================================================
     if ( ig .eq. 0 .and. iv .eq. 0 ) then

        ! South pole
        sum = 0.e+0_fp
        nlon= 0.0d0
        do i=1,im
           if(abs(q2(i,1)-miss)>tiny_r8 ) then
              sum = sum + q2(i,1)
              nlon= nlon + 1.0d0
           endif
        enddo

        if ( nlon > 0.0d0 ) sum = sum / nlon
        do i=1,im
           q2(i,1) = sum
        enddo

        ! North pole:
        sum = 0.e+0_fp
        nlon= 0.0d0
        do i=1,im
           if( abs(q2(i,jn)-miss)>tiny_r8 ) then
              sum = sum + q2(i,jn)
              nlon= nlon + 1.0d0
           endif
        enddo

        if ( nlon > 0.0d0 ) sum = sum / DBLE( im )
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
! !IROUTINE: Ymap_r4r8
!
! !DESCRIPTION: Routine to perform area preserving mapping in N-S from an
!  arbitrary resolution to another.  The input argument has REAL*4 precision
!  but the output argument has REAL(fp) precision.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ymap_r4r8(im, jm, sin1, q1, jn, sin2, q2, ig, iv, missval)
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
! !OPTIONAL INPUT PARAMETERS:
!
    REAL*4,  INTENT(IN), OPTIONAL :: missval
!
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: i, j0, m, mm, j
    REAL*8               :: dy1(jm)
    REAL*8               :: dy
    REAL*8               :: qsum, dlat, nlon, sum
    REAL*4               :: miss

    ! YMAP begins here!
    do j=1,jm-ig
       dy1(j) = sin1(j+1) - sin1(j)
    enddo

    ! Missing value
    miss = miss_r4
    if ( present(missval) ) miss=missval

    !===============================================================
    ! Area preserving mapping
    !===============================================================

    !$OMP PARALLEL DO                                &
    !$OMP DEFAULT( SHARED                          ) &
    !$OMP PRIVATE( I, J0, J, M, QSUM, DLAT, MM, DY )
    do 1000 i=1,im
       qsum = 0.0d0
       dlat = 0.0d0
       j0 = 1
       do 555 j=1,jn-ig
       do 100 m=j0,jm-ig

          !=========================================================
          ! locate the southern edge: sin2(i)
          !=========================================================
          if(sin2(j) .ge. sin1(m) .and. sin2(j) .le. sin1(m+1)) then

             if(sin2(j+1) .le. sin1(m+1)) then

                ! entire new cell is within the original cell
                if ( abs(q1(i,m)-miss)>tiny_r4 ) q2(i,j)=q1(i,m)
                j0 = m
                goto 555
             else

                ! South most fractional area
                if( abs(q1(i,m)-miss)>tiny_r4 ) then
                   dlat= sin1(m+1)-sin2(j)
                   qsum=(sin1(m+1)-sin2(j))*q1(i,m)
                endif

                do mm=m+1,jm-ig

                   ! locate the northern edge: sin2(j+1)
                   if(sin2(j+1) .gt. sin1(mm+1) ) then

                      ! Whole layer
                      if( abs(q1(i,mm)-miss)>tiny_r4 ) then
                         qsum = qsum + dy1(mm)*q1(i,mm)
                         dlat = dlat + dy1(mm)
                      endif
                   else

                      ! North most fractional area
                      if( abs(q1(i,mm)-miss)>tiny_r4 ) then
                         dy = sin2(j+1)-sin1(mm)
                         qsum=qsum+dy*q1(i,mm)
                         dlat=dlat+dy
                      endif
                      j0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100    continue
123    if ( dlat /= 0.0d0 ) q2(i,j) = qsum / dlat
555    continue
1000 continue
     !$OMP END PARALLEL DO

     !===================================================================
     ! Final processing for poles
     !===================================================================
     if ( ig .eq. 0 .and. iv .eq. 0 ) then

        ! South pole
        sum = 0.e+0_fp
        nlon= 0.0d0
        do i=1,im
           if( abs(q2(i,1)-miss)>tiny_r4 ) then
              sum = sum + q2(i,1)
              nlon = nlon + 1.0d0
           endif
        enddo

        if ( nlon > 0.0d0 ) sum = sum / nlon
        do i=1,im
           q2(i,1) = sum
        enddo

        ! North pole:
        sum = 0.e+0_fp
        nlon = 0.0d0
        do i=1,im
           if( abs(q2(i,jn)-miss)>tiny_r4 ) then
              sum = sum + q2(i,jn)
              nlon = nlon + 1.0d0
           endif
        enddo

        if ( nlon > 0.0d0 ) sum = sum / nlon
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
! !IROUTINE: Ymap_r8r4
!
! !DESCRIPTION: Routine to perform area preserving mapping in N-S from an
!  arbitrary resolution to another.  The input argument has REAL*8 precision
!  but the output argument has REAL*4 precision.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ymap_r8r4(im, jm, sin1, q1, jn, sin2, q2, ig, iv, missval)
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
! !OPTIONAL INPUT PARAMETERS:
!
    REAL*8,  INTENT(IN), OPTIONAL :: missval
!
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: i, j0, m, mm, j
    REAL*8               :: dy1(jm)
    REAL*8               :: dy
    REAL*8               :: qsum, sum, dlat
    REAL*8               :: miss
    REAL*4               :: nlon

    ! YMAP begins here!
    do j=1,jm-ig
       dy1(j) = sin1(j+1) - sin1(j)
    enddo

    ! Missing value
    miss = miss_r8
    if ( present(missval) ) miss=missval

    !===============================================================
    ! Area preserving mapping
    !===============================================================

    !$OMP PARALLEL DO                                &
    !$OMP DEFAULT( SHARED                          ) &
    !$OMP PRIVATE( I, J0, J, M, QSUM, DLAT, MM, DY )
    do 1000 i=1,im
       qsum = 0.0d0
       dlat = 0.0d0
       j0 = 1
       do 555 j=1,jn-ig
       do 100 m=j0,jm-ig

          !=========================================================
          ! locate the southern edge: sin2(i)
          !=========================================================
          if(sin2(j) .ge. sin1(m) .and. sin2(j) .le. sin1(m+1)) then

             if(sin2(j+1) .le. sin1(m+1)) then

                ! entire new cell is within the original cell
                if( abs(q1(i,m)-miss)>tiny_r8 ) q2(i,j)=q1(i,m)
                j0 = m
                goto 555
             else

                ! South most fractional area
                if( abs(q1(i,m)-miss)>tiny_r8 ) then
                   dlat= sin1(m+1)-sin2(j)
                   qsum=(sin1(m+1)-sin2(j))*q1(i,m)
                endif

                do mm=m+1,jm-ig

                   ! locate the northern edge: sin2(j+1)
                   if(sin2(j+1) .gt. sin1(mm+1) ) then

                      ! Whole layer
                      if( abs(q1(i,mm)-miss)>tiny_r8 ) then
                         qsum = qsum + dy1(mm)*q1(i,mm)
                         dlat = dlat + dy1(mm)
                      endif
                   else

                      ! North most fractional area
                      dy = sin2(j+1)-sin1(mm)
                      if( abs(q1(i,mm)-miss)>tiny_r8 ) then
                         qsum=qsum+dy*q1(i,mm)
                         dlat=dlat+dy
                      endif
                      j0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100    continue
123    if ( dlat /= 0.0d0 ) q2(i,j) = qsum / dlat
555    continue
1000 continue
     !$OMP END PARALLEL DO

     !===================================================================
     ! Final processing for poles
     !===================================================================
     if ( ig .eq. 0 .and. iv .eq. 0 ) then

        ! South pole
        sum = 0.0_f4
        nlon= 0.0
        do i=1,im
           if( abs(q2(i,1)-miss)>tiny_r8 ) then
              sum = sum + q2(i,1)
              nlon= nlon + 1.0
           endif
        enddo

        if ( nlon > 0.0 ) sum = sum / nlon
        do i=1,im
           q2(i,1) = sum
        enddo

        ! North pole:
        sum = 0.0_f4
        nlon= 0.
        do i=1,im
           if( abs(q2(i,jn)-miss)>tiny_r8 ) then
              sum = sum + q2(i,jn)
              nlon= nlon + 1.0
           endif
        enddo

        if ( nlon > 0.0 ) sum = sum / nlon
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
! !IROUTINE: Ymap_r4r4
!
! !DESCRIPTION: Routine to perform area preserving mapping in N-S from an
!  arbitrary resolution to another.  Both the input and output arguments
!  have REAL(fp) precision.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ymap_r4r4(im, jm, sin1, q1, jn, sin2, q2, ig, iv, missval)
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
! !OPTIONAL INPUT PARAMETERS:
!
    ! Missing value
    REAL*4,  INTENT(IN), OPTIONAL  :: missval
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    REAL*4               :: dlat, nlon, miss

    ! YMAP begins here!
    do j=1,jm-ig
       dy1(j) = sin1(j+1) - sin1(j)
    enddo

    ! missing value
    miss = miss_r4
    if ( present(missval) ) miss = missval

    !===============================================================
    ! Area preserving mapping
    !===============================================================

    !$OMP PARALLEL DO                                &
    !$OMP DEFAULT( SHARED                          ) &
    !$OMP PRIVATE( I, J0, J, M, QSUM, DLAT, MM, DY )
    do 1000 i=1,im
       qsum = 0.0
       dlat = 0.0
       j0 = 1
       do 555 j=1,jn-ig
       do 100 m=j0,jm-ig

          !=========================================================
          ! locate the southern edge: sin2(i)
          !=========================================================
          if(sin2(j) .ge. sin1(m) .and. sin2(j) .le. sin1(m+1)) then

             if(sin2(j+1) .le. sin1(m+1)) then

                ! entire new cell is within the original cell
                if( abs(q1(i,m)-miss)>tiny_r4 ) q2(i,j)=q1(i,m)
                j0 = m
                goto 555
             else

                ! South most fractional area
                if( abs(q1(i,m)-miss)>tiny_r4 ) then
                   dlat=sin1(m+1)-sin2(j)
                   qsum=dlat*q1(i,m)
                endif

                do mm=m+1,jm-ig

                   ! locate the northern edge: sin2(j+1)
                   if(sin2(j+1) .gt. sin1(mm+1) ) then

                      ! Whole layer
                      if( abs(q1(i,mm)-miss)>tiny_r4 ) then
                         qsum = qsum + dy1(mm)*q1(i,mm)
                         dlat = dlat + dy1(mm)
                      endif
                   else

                      ! North most fractional area
                      if( abs(q1(i,mm)-miss)>tiny_r4 ) then
                         dy = sin2(j+1)-sin1(mm)
                         qsum=qsum+dy*q1(i,mm)
                         dlat=dlat+dy
                      endif
                      j0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100    continue
!123    q2(i,j) = qsum / ( sin2(j+1) - sin2(j) )
123    if ( dlat /= 0.0 ) q2(i,j) = qsum / dlat
555    continue
1000 continue
     !$OMP END PARALLEL DO

     !===================================================================
     ! Final processing for poles
     !===================================================================
     if ( ig .eq. 0 .and. iv .eq. 0 ) then

        ! South pole
        sum  = 0.e+0_fp
        nlon = 0.0
        do i=1,im
           if( abs(q2(i,1)-miss)>tiny_r4 ) then
              sum  = sum + q2(i,1)
              nlon = nlon + 1.0
           endif
        enddo

        if ( nlon > 0.0 ) sum = sum / nlon
        !sum = sum / REAL( im, 4 )
        do i=1,im
           q2(i,1) = sum
        enddo

        ! North pole:
        sum = 0.e+0_fp
        nlon= 0.0
        do i=1,im
           if( abs(q2(i,jn)-miss)>tiny_r4 ) then
              sum  = sum + q2(i,jn)
              nlon = nlon + 1.0
           endif
        enddo

        !sum = sum / REAL( im, 4 )
        if ( nlon > 0.0 ) sum = sum / nlon
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
! !IROUTINE: Xmap_r8r8
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
  SUBROUTINE xmap_r8r8(im, jm, lon1, q1, iin, ilon2, iq2, missval)
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
! !OPTIONAL INPUT PARAMETERS:
!
    ! Missing value
    REAL*8,  INTENT(IN), OPTIONAL  :: missval
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    REAL*8               :: qsum, dlon
    LOGICAL              :: found

    ! Update
    INTEGER              :: n1, n2
    INTEGER              :: in
    REAL*8, POINTER      :: lon2(:)
    REAL*8, POINTER      :: q2(:,:)
    REAL*8               :: minlon, maxlon
    REAL*8               :: lon1s(im+1)

    ! Ghost correction
    Logical              :: isGlobal
    Real*8               :: xSpan

    ! Missing value
    Real*8               :: miss

    ! Initialize pointers
    lon2 => NULL()
    q2   => NULL()

    ! missing value
    miss = miss_r8
    if ( present(missval) ) miss = missval

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

    ! if there is no overlap between original grid and output grid
    ! reduced will be zero and missing values should be returned
    if ( in .eq. 0 ) then
       iq2 = missval
       lon2 => NULL()
       q2 => NULL()
       return
    endif

    ! Periodic BC only valid if the variable is "global"
    xSpan = x1(im+1)-x1(1)
    isGlobal = ((xSpan.ge.355.0).and.(xSpan.le.365.0))

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

    !$OMP PARALLEL DO                                      &
    !$OMP DEFAULT( SHARED                                ) &
    !$OMP PRIVATE( J, QTMP, I, I0, M, QSUM, DLON, MM, DX )
    do 1000 j=1,jm

       !=================================================================
       ! Area preserving mapping
       !================================================================

       qtmp(:) = 0.0d0
       do i=1,im
          qtmp(i)=q1(i,j)
       enddo

       ! SDE 2017-01-07
       ! Only have shadow regions if we are on a global grid. Otherwise, we
       ! should keep the zero boundary conditions.
       If (isGlobal) Then
          qtmp(0)=q1(im,j)
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
       End If

       i0 = i1

       do 555 i=1,in
       do 100 m=i0,i2-1

          !=============================================================
          ! locate the western edge: lon2(i)
          !=============================================================
          if(lon2(i) .ge. x1(m) .and. lon2(i) .le. x1(m+1)) then

             if(lon2(i+1) .le. x1(m+1)) then

                ! entire new grid is within the original grid
                if( abs(qtmp(m)-miss)>tiny_r8 ) q2(i,j)=qtmp(m)
                i0 = m
                goto 555
             else

                ! Left most fractional area
                if( abs(qtmp(m)-miss)>tiny_r8 ) then
                   qsum=(x1(m+1)-lon2(i))*qtmp(m)
                   dlon= x1(m+1)-lon2(i)
                else
                   qsum = 0.0d0
                   dlon = 0.0d0
                endif
                do mm=m+1,i2-1

                   ! locate the eastern edge: lon2(i+1)
                   if(lon2(i+1) .gt. x1(mm+1) ) then

                      ! Whole layer
                      if( abs(qtmp(mm)-miss)>tiny_r8 ) then
                         qsum = qsum + dx1(mm)*qtmp(mm)
                         dlon = dlon + dx1(mm)
                      endif
                   else
                      ! Right most fractional area
                      if( abs(qtmp(mm)-miss)>tiny_r8 ) then
                         dx = lon2(i+1)-x1(mm)
                         qsum=qsum+dx*qtmp(mm)
                         dlon=dlon+dx
                      endif
                      i0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100    continue
123    if ( dlon /= 0.0d0 ) q2(i,j) = qsum / dlon
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
! !IROUTINE: Xmap_r4r4
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
  SUBROUTINE xmap_r4r4(im, jm, lon1, q1, iin, ilon2, iq2, missval)
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
! !OPTIONAL INPUT PARAMETERS:
!
    ! Missing value
    REAL*4,  INTENT(IN), OPTIONAL  :: missval
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    REAL*4               :: qsum, dlon
    LOGICAL              :: found

    ! Update
    INTEGER              :: n1, n2
    INTEGER              :: in
    REAL*4, POINTER      :: lon2(:)
    REAL*4, POINTER      :: q2(:,:)
    REAL*4               :: minlon, maxlon
    REAL*4               :: lon1s(im+1)

    ! Ghost correction
    Logical              :: isGlobal
    Real*4               :: xSpan

    ! Missing value
    REAL*4               :: miss

    ! Initialize
    lon2 => NULL()
    q2   => NULL()

    ! Missing value
    miss = miss_r4
    if ( present(missval) ) miss = missval

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
    ! ckeller, (2/11/15).
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

    ! if there is no overlap between original grid and output grid
    ! reduced will be zero and missing values should be returned
    if ( in .eq. 0 ) then
       iq2 = missval
       lon2 => NULL()
       q2 => NULL()
       return
    endif

    ! shadow variables to selected range
    ! Periodic BC only valid if the variable is "global"
    xSpan = x1(im+1)-x1(1)
    isGlobal = ((xSpan.ge.355.0).and.(xSpan.le.365.0))

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

    !$OMP PARALLEL DO                                      &
    !$OMP DEFAULT( SHARED                                ) &
    !$OMP PRIVATE( J, QTMP, I, I0, M, QSUM, DLON, MM, DX )
    do 1000 j=1,jm

       !=================================================================
       ! Area preserving mapping
       !================================================================

       qtmp(:) = 0.0
       do i=1,im
          qtmp(i)=q1(i,j)
       enddo

       ! SDE 2017-01-07
       ! Only have shadow regions if we are on a global grid. Otherwise, we
       ! should keep the zero boundary conditions.
       If (isGlobal) Then
          qtmp(0)=q1(im,j)
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
       End If

       i0 = i1

       do 555 i=1,in
       do 100 m=i0,i2-1

          !=============================================================
          ! locate the western edge: lon2(i)
          !=============================================================
          if(lon2(i) .ge. x1(m) .and. lon2(i) .le. x1(m+1)) then

             if(lon2(i+1) .le. x1(m+1)) then

                ! entire new grid is within the original grid
                if ( abs(qtmp(m)-miss)>tiny_r4 ) q2(i,j)=qtmp(m)
                i0 = m
                goto 555
             else

                ! Left most fractional area
                if( abs(qtmp(m)-miss)>tiny_r4 ) then
                   dlon=x1(m+1)-lon2(i)
                   qsum=dlon*qtmp(m)
                else
                   dlon=0.0
                   qsum=0.0
                endif
                do mm=m+1,i2-1

                   ! locate the eastern edge: lon2(i+1)
                   if(lon2(i+1) .gt. x1(mm+1) ) then

                      ! Whole layer
                      if( abs(qtmp(mm)-miss)>tiny_r4 ) then
                         qsum = qsum + dx1(mm)*qtmp(mm)
                         dlon = dlon + dx1(mm)
                      endif

                   else
                      ! Right most fractional area
                      if( abs(qtmp(mm)-miss)>tiny_r4 ) then
                         dx = lon2(i+1)-x1(mm)
                         qsum=qsum+dx*qtmp(mm)
                         dlon=dlon+dx
                      endif
                      i0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100    continue
123    if( dlon > 0.0 ) q2(i,j) = qsum / dlon
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
! !IROUTINE: Xmap_r4r8
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
  SUBROUTINE xmap_r4r8(im, jm, lon1, q1, iin, ilon2, iq2, missval)
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
! !OPTIONAL INPUT PARAMETERS:
!
    REAL*4,  INTENT(IN), OPTIONAL :: missval
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    REAL*8               :: qsum, dlon
    LOGICAL              :: found

    ! Update
    INTEGER              :: n1, n2
    INTEGER              :: in
    REAL*4, POINTER      :: lon2(:)
    REAL*8, POINTER      :: q2(:,:)
    REAL*4               :: minlon, maxlon
    REAL*4               :: lon1s(im+1)

    ! Ghost correction
    Logical              :: isGlobal
    Real*8               :: xSpan

    ! Missing value
    REAL*4               :: miss

    ! Initialize
    lon2 => NULL()
    q2   => NULL()

    ! Missing value
    miss = miss_r4
    if ( present(missval) ) miss = missval

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

    ! Periodic BC only valid if the variable is "global"
    xSpan = x1(im+1)-x1(1)
    isGlobal = ((xSpan.ge.355.0).and.(xSpan.le.365.0))

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

    !$OMP PARALLEL DO                                      &
    !$OMP DEFAULT( SHARED                                ) &
    !$OMP PRIVATE( J, QTMP, I, I0, M, QSUM, DLON, MM, DX )
    do 1000 j=1,jm

       !=================================================================
       ! Area preserving mapping
       !================================================================

       qtmp(:) = 0.0d0
       do i=1,im
          qtmp(i)=q1(i,j)
       enddo

       ! SDE 2017-01-07
       ! Only have shadow regions if we are on a global grid. Otherwise, we
       ! should keep the zero boundary conditions.
       If (isGlobal) Then
          qtmp(0)=q1(im,j)
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
       End If

       i0 = i1

       do 555 i=1,in
       do 100 m=i0,i2-1

          !=============================================================
          ! locate the western edge: lon2(i)
          !=============================================================
          if(lon2(i) .ge. x1(m) .and. lon2(i) .le. x1(m+1)) then

             if(lon2(i+1) .le. x1(m+1)) then

                ! entire new grid is within the original grid
                if( abs(qtmp(m)-miss)>tiny_r4 ) q2(i,j)=qtmp(m)
                i0 = m
                goto 555
             else

                ! Left most fractional area
                if( abs(qtmp(m)-miss)>tiny_r4 ) then
                   qsum=(x1(m+1)-lon2(i))*qtmp(m)
                   dlon= x1(m+1)-lon2(i)
                else
                   qsum=0.0d0
                   dlon=0.0d0
                endif
                do mm=m+1,i2-1

                   ! locate the eastern edge: lon2(i+1)
                   if(lon2(i+1) .gt. x1(mm+1) ) then

                      ! Whole layer
                      if( abs(qtmp(mm)-miss)>tiny_r4 ) then
                         qsum = qsum + dx1(mm)*qtmp(mm)
                         dlon = dlon + dx1(mm)
                      endif
                   else
                      ! Right most fractional area
                      if( abs(qtmp(mm)-miss)>tiny_r4 ) then
                         dx = lon2(i+1)-x1(mm)
                         qsum=qsum+dx*qtmp(mm)
                         dlon=dlon+dx
                      endif
                      i0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100    continue
123    if ( dlon /= 0.0d0 ) q2(i,j) = qsum / dlon
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
! !IROUTINE: Xmap_r8r4
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
  SUBROUTINE xmap_r8r4(im, jm, lon1, q1, iin, ilon2, iq2, missval)
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
! !OPTIONAL INPUT PARAMETERS:
!
    REAL*8,  INTENT(IN), OPTIONAL :: missval
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    REAL*4               :: qsum, dlon
    LOGICAL              :: found

    ! Update
    INTEGER              :: n1, n2
    INTEGER              :: in
    REAL*4, POINTER      :: lon2(:)
    REAL*4, POINTER      :: q2(:,:)
    REAL*4               :: minlon, maxlon
    REAL*4               :: lon1s(im+1)

    ! Ghost correction
    Logical              :: isGlobal
    Real*4               :: xSpan

    ! Missing value
    REAL*8               :: miss

    ! Initialize
    lon2 => NULL()
    q2   => NULL()

    ! Missing value
    miss = miss_r8
    if ( present(missval) ) miss = missval

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

    ! Periodic BC only valid if the variable is "global"
    xSpan = x1(im+1)-x1(1)
    isGlobal = ((xSpan.ge.355.0).and.(xSpan.le.365.0))

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

    !$OMP PARALLEL DO                                      &
    !$OMP DEFAULT( SHARED                                ) &
    !$OMP PRIVATE( J, QTMP, I, I0, M, QSUM, DLON, MM, DX )
    do 1000 j=1,jm

       !=================================================================
       ! Area preserving mapping
       !================================================================

       qtmp(:) = 0.0
       do i=1,im
          qtmp(i)=q1(i,j)
       enddo

       ! SDE 2017-01-07
       ! Only have shadow regions if we are on a global grid. Otherwise, we
       ! should keep the zero boundary conditions.
       If (isGlobal) Then
          qtmp(0)=q1(im,j)
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
       End If

       i0 = i1

       do 555 i=1,in
       do 100 m=i0,i2-1

          !=============================================================
          ! locate the western edge: lon2(i)
          !=============================================================
          if(lon2(i) .ge. x1(m) .and. lon2(i) .le. x1(m+1)) then

             if(lon2(i+1) .le. x1(m+1)) then

                ! entire new grid is within the original grid
                if( abs(qtmp(m)-miss)>tiny_r8 ) q2(i,j)=qtmp(m)
                i0 = m
                goto 555
             else

                ! Left most fractional area
                if( abs(qtmp(m)-miss)>tiny_r8 ) then
                   qsum=(x1(m+1)-lon2(i))*qtmp(m)
                   dlon= x1(m+1)-lon2(i)
                else
                   qsum=0.0
                   dlon=0.0
                endif
                do mm=m+1,i2-1

                   ! locate the eastern edge: lon2(i+1)
                   if(lon2(i+1) .gt. x1(mm+1) ) then

                      ! Whole layer
                      if( abs(qtmp(mm)-miss)>tiny_r8 ) then
                         qsum = qsum + dx1(mm)*qtmp(mm)
                         dlon = dlon + dx1(mm)
                      endif
                   else
                      ! Right most fractional area
                      if( abs(qtmp(m)-miss)>tiny_r8 ) then
                         dx = lon2(i+1)-x1(mm)
                         qsum=qsum+dx*qtmp(mm)
                         dlon=dlon+dx
                      endif
                      i0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100    continue
123    if( dlon /= 0.0 ) q2(i,j) = qsum / dlon
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
    REAL(fp),         INTENT(IN) :: LONS (NX+1 )   ! Longitudes
    REAL(fp),         INTENT(IN) :: SINES(NY+1 )   ! Sines of latitudes
    REAL(fp),         INTENT(IN) :: AREAS(NX,NY)   ! Surface areas [m2]
    CHARACTER(LEN=*), INTENT(IN) :: DIR            ! Dir for netCDF files w/
                                                   !  grid definitions
!
! !REVISION HISTORY:
!  14 Jul 2014 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
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
    OUTNX   = NX
    OUTNY   = NY
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
!  See https://github.com/geoschem/geos-chem for complete history
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
