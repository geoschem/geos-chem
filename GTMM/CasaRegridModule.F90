!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CasaRegridModule
!
! !DESCRIPTION: Module CasaRegridModule contains arrays and variables used to 
!  regrid the GEOS-5 data from 1 x 1 Generic to 2 x 2.5, and 4 x 5 geos grids.
!\\
!\\
! !INTERFACE: 
!
MODULE CasaRegridModule
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: regrid1x1to4x5
  PUBLIC :: regrid4x5to1x1
  PUBLIC :: regrid1x1to2x25
  PUBLIC :: regrid2x25to1x1
  PUBLIC :: CasaRegridInit
!
! !PUBLIC DATA MEMBERS:
!
  PUBLIC :: I1x1,      J1x1,      L1x1    
  PUBLIC :: I2x25,     J2x25,     L2x25   
  PUBLIC :: I4x5,      J4x5,      L4x5    
!
! !DEFINED PARAMETERS:
!
  INTEGER, PARAMETER  :: I1x1     = 360,  J1x1     = 180,  L1x1     = 72    
  INTEGER, PARAMETER  :: I2x25    = 144,  J2x25    = 91,   L2x25    = 72
  INTEGER, PARAMETER  :: I4x5     = 72,   J4x5     = 46,   L4x5     = 72
!
! !PRIVATE TYPES:
!
  ! Degrees to Radians
  REAL*8,  PARAMETER  :: D2R = 3.141592658979323d0 / 180d0

  !---------------------
  ! 1 x 1 grid
  !---------------------

  ! Lon edges
  REAL*8 :: xedge_1x1( I1x1 + 1 ) 

  ! Lat edges
  REAL*8 :: yedge_1x1( J1x1 + 1 )

  ! Sine of latitude
  REAL*8 :: sine_1x1( J1x1 + 1 )

  !---------------------
  ! 2 x 2.5 grid
  !---------------------

  ! Longitude edges
  REAL*8 :: xedge_2x25( I2x25 + 1 ) 

  ! Latitude edges
  REAL*8 :: yedge_2x25( J2x25 + 1 )

  ! Latitude edges
  REAL*8 :: sine_2x25( J2x25 + 1 )

  !---------------------
  ! 4 x 5 grid
  !---------------------

  ! Longitude edges
  REAL*8 :: xedge_4x5( I4x5 + 1 )

  ! Latitude edges
  REAL*8 :: yedge_4x5( J4x5 + 1 )

  ! Latitude edges
  REAL*8 :: sine_4x5( J4x5 + 1 )
!
! !REMARKS:
!  CasaRegridModule uses the regridding software "MAP_A2A" from S-J Lin.  
!  This is area-preserving mapping.  For example, if you have a quantity 
!  such as kg/m2/s or W/m2, MAP_A2A will multiply by the area on the
!  input grid, then regrid, and divide by the area on the output grid,
!  such that the total quantity is preserved.  
!
! !REVISION HISTORY:
!  14 Jan 2008 - R. Yantosca - Initial version
!  (1 ) Modify regriddGeos5To* routines so that if all values are zero,
!        then we just fill the output data array with zeros and return.
!        This ought to speed up program execution. (bmy, 11/14/06)
!EOP
!------------------------------------------------------------------------------
!BOC
  CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: regrid4x5to1x1
!
! !DESCRIPTION: Subroutine regrid4x5to1x1 is a wrapper for MAP\_A2A, which is 
!  called to regrid from the GEOS-5 4x5 grid to the GENERIC 1x1 grid. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE regrid4x5to1x1( iv, q1, q2 )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN)  :: iv
    REAL*8,  INTENT(IN)  :: q1( I4x5, J4x5 )
!
! !OUTPUT PARAMETERS:
!
    REAL*8,  INTENT(OUT) :: q2( I1x1, J1x1 )
!
! !REVISION HISTORY: 
!  08 Nov 2006 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !---------------------------------
    ! regridGeos5to1x1 begins here!
    !---------------------------------

    ! If all elements of q1 are zero, then set q2=0 and return
    IF ( ALL( q1 == 0d0 ) ) THEN
       q2 = 0d0
       RETURN
    ENDIF

    ! Call MAP_A2A to do the horizontal regridding
    CALL map_a2a( I4x5, J4x5, xedge_4x5, sine_4x5, q1, &
                  I1x1, J1x1, xedge_1x1, sine_1x1, q2, 0, iv )

  END SUBROUTINE regrid4x5to1x1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: regrid1x1to4x5
!
! !DESCRIPTION: Subroutine regrid1x1to4x5 is a wrapper for MAP\_A2A, which is 
!  called to regrid from the GENERIC 1x1 grid to the GEOS-5 4x5 grid.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE regrid1x1to4x5( iv, q1, q2 )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)  :: iv
    REAL*8,  INTENT(IN)  :: q1( I1x1, J1x1 )
!
! !OUTPUT PARAMETERS:
!
    REAL*8,  INTENT(OUT) :: q2( I4x5, J4x5 )
!
! !REVISION HISTORY:
!  08 Nov 2006 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !---------------------------------
    ! regridGeos5to1x1 begins here!
    !---------------------------------

    ! If all elements of q1 are zero, then set q2=0 and return
    IF ( ALL( q1 == 0d0 ) ) THEN
       q2 = 0d0
       RETURN
    ENDIF

    ! Call MAP_A2A to do the horizontal regridding
    CALL map_a2a( I1x1, J1x1, xedge_1x1, sine_1x1, q1, &
                  I4x5, J4x5, xedge_4x5, sine_4x5, q2, 0, iv )

  END SUBROUTINE regrid1x1to4x5
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: regrid2x25to1x1
!
! !DESCRIPTION: Subroutine regrid2x25to1x1 is a wrapper for MAP\_A2A, which is 
!  called to regrid from the GENERIC 1x1 grid to the GEOS 2 x 2.5 grid.
!  (bmy, 11/8/06)
!\\
!\\
! !INTERFACE
!
  SUBROUTINE regrid2x25to1x1( iv, q1, q2 )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)  :: iv
    REAL*8,  INTENT(IN)  :: q1( I2x25, J2x25 )
!
! !OUTPUT PARAMETERS:
!
    REAL*8,  INTENT(OUT) :: q2( I1x1,  J1x1  )
!
! !REVISION HISTORY:
!  08 Nov 2006 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !---------------------------------
    ! regridGeos5to1x125 begins here!
    !---------------------------------

    ! If all elements of q1 are zero, then set q2=0 and return
    IF ( ALL( q1 == 0d0 ) ) THEN
       q2 = 0d0
       RETURN
    ENDIF

    ! Call MAP_A2A to do the horizontal regridding
    CALL map_a2a( I2x25, J2x25, xedge_2x25, sine_2x25, q1, &
                  I1x1,  J1x1,  xedge_1x1,  sine_1x1,  q2, 0, iv )

  END SUBROUTINE regrid2x25to1x1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: regrid1x1to2x25
!
! !DESCRIPTION: Subroutine regridGeos5to2x25 is a wrapper for MAP\_A2A, which 
!  regrids from the GEOS-5 1x1 grid to the GEOS 2 x 2.5 grid.  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE regrid1x1to2x25( iv, q1, q2 )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)  :: iv
    REAL*8,  INTENT(IN)  :: q1( I1x1,  J1x1  )
!
! !OUTPUT PARAMETERS:
!
    REAL*8,  INTENT(OUT) :: q2( I2x25, J2x25 )
!
! !REVISION HISTORY:
!  08 Nov 2006 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !---------------------------------
    ! regridGeos5to2x25 begins here!
    !---------------------------------

    ! If all elements of q1 are zero, then set q2=0 and return
    IF ( ALL( q1 == 0d0 ) ) THEN
       q2 = 0d0
       RETURN
    ENDIF

    ! Call MAP_A2A to do the horizontal regridding
    CALL map_a2a( I1x1,  J1x1,  xedge_1x1,  sine_1x1,  q1, &
                  I2x25, J2x25, xedge_2x25, sine_2x25, q2, 0, iv )

  END SUBROUTINE regrid1x1to2x25
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_a2a
!
! !DESCRIPTION: Subroutine MAP\_A2A is a orizontal arbitrary grid to arbitrary 
!  grid conservative high-order mapping regridding routine by S-J Lin.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE map_a2a( im, jm, lon1, sin1, q1, &
                      in, jn, lon2, sin2, q2, ig, iv)
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)  :: im, jm, in, jn, ig, iv
    REAL*8,  INTENT(IN)  :: lon1(im+1), lon2(in+1)
    REAL*8,  INTENT(IN)  :: sin1(jm+1), sin2(jn+1)
    REAL*8,  INTENT(IN)  :: q1(im,jm)
!
! !OUTPUT PARAMETERS:
!
    REAL*8,  INTENT(OUT) :: q2(in,jn)
!
!  !REVISION HISTORY:
!  (1) Original subroutine by S-J Lin.  Converted to F90 freeform format
!      and inserted into "Geos3RegridModule" by Bob Yantosca (9/21/00)
!  (2) Added F90 type declarations to be consistent w/ TypeModule.f90.
!      Also updated comments. (bmy, 9/21/00)
!  21 Sep 2000 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: i,j,k
    REAL*8               :: qtmp(in,jm)

    !===================================================================
    ! MAP_A2A begins here!
    !
    ! Mapping in the E-W direction
    ! If both grids have the same longitude dimension, don't call XMAP
    !===================================================================    
    IF ( im .eq. in ) THEN
       DO j=1,jm-ig
       DO i=1,im
          qtmp(i,j+ig) = q1(i,j+ig)
       ENDDO
       ENDDO
    ELSE
       CALL xmap(im, jm-ig, lon1, q1(1,1+ig),in, lon2, qtmp(1,1+ig) )
    ENDIF
    
    !===================================================================
    ! Mapping in the N-S direction
    ! If both grids have the same latitude dimension, don't call YMAP 
    !===================================================================    
    IF ( jm .eq. jn ) THEN
       DO j=1,jm-ig
       DO i=1,in
          q2(i,j+ig) = qtmp(i,j+ig)
       ENDDO
       ENDDO
    ELSE
       CALL ymap(in, jm, sin1, qtmp(1,1+ig), jn, sin2, q2(1,1+ig), ig, iv)
    ENDIF

  END SUBROUTINE map_a2a
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ymap
!
! !DESCRIPTION: Routine to perform area preserving mapping in N-S from an 
!  arbitrary resolution to another.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ymap(im, jm, sin1, q1, jn, sin2, q2, ig, iv)
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)  :: im            ! original E-W dimension
    INTEGER, INTENT(IN)  :: jm            ! original N-S dimension
    INTEGER, INTENT(IN)  :: jn            ! Target N-S dimension
    INTEGER, INTENT(IN)  :: ig            ! ig=0: scalars from SP to NP
                                          ! D-grid v-wind is also ig 0
                                          ! ig=1: D-grid u-wind
    INTEGER, INTENT(IN)  :: iv            ! iv=0 scalar; iv=1: vector
    REAL*8,  INTENT(IN)  :: sin1(jm+1-ig) ! original southern edge of 
                                          !  the cell sin(lat1)  
    REAL*8,  INTENT(IN)  :: q1(im,jm)     ! original data at center of 
                                          !  the cell
    REAL*8,  INTENT(IN)  :: sin2(jn+1-ig) ! Target cell's southern edge
                                          !  sin(lat2)
!
! !OUTPUT PARAMETERS:
!
    REAL*8,  INTENT(OUT) :: q2(im,jn)     ! Mapped data at the 
                                          !  target resolution
!
! !REMARKS:
!  sin1 (1) = -1 must be south pole; sin1(jm+1)=1 must be N pole.
!                                                                             .
!  sin1(1) < sin1(2) < sin1(3) < ... < sin1(jm) < sin1(jm+1)
!  sin2(1) < sin2(2) < sin2(3) < ... < sin2(jn) < sin2(jn+1)
!
! !AUTHOR:
!  S.-J. Lin
!  First version: piece-wise constant mapping
!  Apr 1, 2000
!  Last modified:
!
! !REVISION HISTORY:
!  21 Sep 2000 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: i, j0, m, mm, j
    REAL*8               :: al(im,jm), ar(im,jm), a6(im,jm), dy1(jm)
    REAL*8,  PARAMETER   :: r3 = 1./3., r23 = 2./3. 
    REAL*8               :: pl, pr, qsum, esl, dy, sum
    
    ! YMAP begins here!
    do j=1,jm-ig
       dy1(j) = sin1(j+1) - sin1(j)
    enddo

    !===============================================================
    ! Area preserving mapping
    !===============================================================

    ! Construct subgrid PP distribution
    call ppm_lat(im, jm, ig, q1, al, ar, a6, 3, iv)
    
    do 1000 i=1,im
       j0 = 1
       do 555 j=1,jn-ig
       do 100 m=j0,jm-ig

          !=========================================================
          ! locate the southern edge: sin2(i)
          !=========================================================
          if(sin2(j) .ge. sin1(m) .and. sin2(j) .le. sin1(m+1)) then
             pl = (sin2(j)-sin1(m)) / dy1(m)
             
             if(sin2(j+1) .le. sin1(m+1)) then
                
                ! entire new cell is within the original cell
                pr = (sin2(j+1)-sin1(m)) / dy1(m)
                q2(i,j) = al(i,m) + 0.5d0*(a6(i,m)+ar(i,m)-al(i,m)) &
               &                    *(pr+pl)-a6(i,m)*r3*(pr*(pr+pl)+pl**2)
                j0 = m
                goto 555
             else

                ! South most fractional area
                qsum = (sin1(m+1)-sin2(j))*(al(i,m)+0.5d0*(a6(i,m)+ &
                &              ar(i,m)-al(i,m))*(1.+pl)-a6(i,m)*  &
                &               (r3*(1.+pl*(1.+pl))))

                do mm=m+1,jm-ig

                   ! locate the eastern edge: sin2(j+1)
                   if(sin2(j+1) .gt. sin1(mm+1) ) then

                      ! Whole layer
                      qsum = qsum + dy1(mm)*q1(i,mm)
                   else

                      ! North most fractional area
                      dy = sin2(j+1)-sin1(mm)
                      esl = dy / dy1(mm)
                      qsum = qsum + dy*(al(i,mm)+0.5d0*esl* &
                     &       (ar(i,mm)-al(i,mm)+a6(i,mm)*(1.-r23*esl)))
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

    !===================================================================
    ! Final processing for poles
    !===================================================================
    if ( ig .eq. 0 .and. iv .eq. 0 ) then
        
       ! South pole
       sum = 0.
       do i=1,im
          sum = sum + q2(i,1)
       enddo
        
       sum = sum / float(im)
       do i=1,im
          q2(i,1) = sum
       enddo
        
       ! North pole:
       sum = 0.
       do i=1,im
          sum = sum + q2(i,jn)
       enddo
        
       sum = sum / float(im)
       do i=1,im
          q2(i,jn) = sum
       enddo

    endif

  END SUBROUTINE YMAP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ppm_lat
!
! !DESCRIPTION: Subroutine PPM\_LAT is called by YMAP.  Written by S-J Lin, and
!  converted to F90 freeform format by Bob Yantosca.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ppm_lat(im, jm, ig, q, al, ar, a6, jord, iv)
!
! INPUT PARAMETERS:
!
    INTEGER           :: im, jm          !  Dimensions
    INTEGER           :: ig
    INTEGER           :: jord
    INTEGER           :: iv              ! iv=0 scalar
                                         ! iv=1 vector
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL*8            :: q(im,jm-ig)
!
! !OUTPUT PARAMETERS:
!
    REAL*8            :: al(im,jm-ig)
    REAL*8            :: ar(im,jm-ig)
    REAL*8            :: a6(im,jm-ig)
!
! !REVISION HISTORY:
!  21 Sep 2000 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*8            :: dm(im,jm-ig)
    REAL*8, PARAMETER :: r3 = 1./3. 
    INTEGER           :: i, j, im2, iop, jm1
    REAL*8            :: tmp, qmax, qmin, qop
    
    ! PPM_LAT begins here
    ! Compute dm: linear slope
    do j=2,jm-1-ig
       do i=1,im
          dm(i,j) = 0.25d0*(q(i,j+1) - q(i,j-1))
          qmax = max(q(i,j-1),q(i,j),q(i,j+1)) - q(i,j)
          qmin = q(i,j) - min(q(i,j-1),q(i,j),q(i,j+1))
          dm(i,j) = sign(min(abs(dm(i,j)),qmin,qmax),dm(i,j))
       enddo
    enddo

    im2 = im/2
    jm1 = jm - 1

    ! Poles:
    if (iv .eq. 1 ) then

       !===============================================================
       ! u-wind (ig=1)
       ! v-wind (ig=0)
       !===============================================================

       ! SP
       do i=1,im
          if( i .le. im2) then
             qop = -q(i+im2,2-ig)
          else
             qop = -q(i-im2,2-ig)
          endif
          tmp = 0.25d0*(q(i,2) - qop)
          qmax = max(q(i,2),q(i,1), qop) - q(i,1)
          qmin = q(i,1) - min(q(i,2),q(i,1), qop)
          dm(i,1) = sign(min(abs(tmp),qmax,qmin),tmp)
       enddo
       
       ! NP
       do i=1,im
          if( i .le. im2) then
             qop = -q(i+im2,jm1)
          else
             qop = -q(i-im2,jm1)
          endif
          tmp = 0.25d0*(qop - q(i,jm1-ig))
          qmax = max(qop,q(i,jm-ig), q(i,jm1-ig)) - q(i,jm-ig)
          qmin = q(i,jm-ig) - min(qop,q(i,jm-ig), q(i,jm1-ig))
          dm(i,jm-ig) = sign(min(abs(tmp),qmax,qmin),tmp)
       enddo
    else
        
       !===============================================================
       ! Scalar:
       ! This code segment currently works only if ig=0
       !===============================================================

       ! SP
       do i=1,im2
          tmp = 0.25d0*(q(i,2)-q(i+im2,2))
          qmax = max(q(i,2),q(i,1), q(i+im2,2)) - q(i,1)
          qmin = q(i,1) - min(q(i,2),q(i,1), q(i+im2,2))
          dm(i,1) = sign(min(abs(tmp),qmax,qmin),tmp)
       enddo
        
       do i=im2+1,im
          dm(i, 1) =  - dm(i-im2, 1)
       enddo

       ! NP
       do i=1,im2
          tmp = 0.25d0*(q(i+im2,jm1)-q(i,jm1))
          qmax = max(q(i+im2,jm1),q(i,jm), q(i,jm1)) - q(i,jm)
          qmin = q(i,jm) - min(q(i+im2,jm1),q(i,jm), q(i,jm1))
          dm(i,jm) = sign(min(abs(tmp),qmax,qmin),tmp)
       enddo

       do i=im2+1,im
          dm(i,jm) =  - dm(i-im2,jm)
       enddo
     endif
      
     do j=2,jm-ig
        do i=1,im
           al(i,j) = 0.5d0*(q(i,j-1)+q(i,j)) + r3*(dm(i,j-1) - dm(i,j))
        enddo
     enddo
      
     do j=1,jm-1-ig
        do i=1,im
           ar(i,j) = al(i,j+1)
        enddo
     enddo
     
     if ( iv .eq. 1 ) then
        
        if ( ig .eq. 0 ) then

           !============================================================
           ! Vector: ig=0
           !============================================================
           do i=1,im2
              al(i,    1) = -al(i+im2,2)
              al(i+im2,1) = -al(i,    2)
           enddo
           
           do i=1,im2
              ar(i,    jm) = -ar(i+im2,jm1)
              ar(i+im2,jm) = -ar(i,    jm1)
           enddo
        else

           !============================================================
           ! ig=1 : SP
           !============================================================
           do i=1,im
              if( i .le. im2) then
                 iop = i+im2
              else
                 iop = i-im2
              endif
              al(i,1) = 0.5d0*(q(i,1)-q(iop,1)) - r3*(dm(iop,1) + dm(i,1))
           enddo

           !============================================================
           ! NP
           !============================================================
           do i=1,im
              if( i .le. im2) then
                 iop = i+im2
              else
                 iop = i-im2
              endif
              ar(i,jm1) = 0.5d0*(q(i,jm1)-q(iop,jm1)) - &
             &                 r3*(dm(iop,jm1) + dm(i,jm1))
            enddo
        endif
     else

        ! Scalar (works for ig=0 only):
        do i=1,im2
           al(i,    1) = al(i+im2,2)
           al(i+im2,1) = al(i,    2)
        enddo
        
        do i=1,im2
           ar(i,    jm) = ar(i+im2,jm1)
           ar(i+im2,jm) = ar(i,    jm1)
        enddo
     endif
      
     do j=1,jm-ig
        do i=1,im
           a6(i,j) = 3d0*(q(i,j)+q(i,j) - (al(i,j)+ar(i,j)))
        enddo
        call lmppm(dm(1,j), a6(1,j), ar(1,j), al(1,j),  q(1,j), im, jord-3)
     enddo

  END SUBROUTINE ppm_lat
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: xmap
!
! !DESCRIPTION: Routine to perform area preserving mapping in E-W from an 
!  arbitrary resolution to another.
!  Periodic domain will be assumed, i.e., the eastern wall bounding cell
!  im is lon1(im+1) = lon1(1); Note the equal sign is true geographysically.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE xmap(im, jm, lon1, q1, in, lon2, q2)
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)  :: im           ! original E-W dimension
    INTEGER, INTENT(IN)  :: in           ! Target E-W dimension
    INTEGER, INTENT(IN)  :: jm           ! original N-S dimension
    REAL*8,  INTENT(IN)  :: lon1(im+1)   ! original western edge of 
                                         !  the cell
    REAL*8,  INTENT(IN)  :: q1(im,jm)    ! original data at center of 
                                         !  the cell
    REAL*8,  INTENT(IN)  :: lon2(in+1)   ! Target cell's western edge
!
! !OUTPUT PARAMETERS:
!
    REAL*8,  INTENT(OUT) :: q2(in,jm)    ! Mapped data at the 
                                         !  target resolution
!
! !REVISION HISTORY:
!  21 Sep 2000 - R. Yantosca - Initial version

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: i1, i2, i, i0, m, mm, j

    REAL*8               :: qtmp(-im:im+im)
    REAL*8               :: al(-im:im+im)
    REAL*8               :: ar(-im:im+im)
    REAL*8               :: a6(-im:im+im)
    REAL*8               :: x1(-im:im+im+1)
    REAL*8               :: dx1(-im:im+im)
    REAL*8,  PARAMETER   :: r3 = 1./3., r23 = 2./3. 
    REAL*8               :: pl, pr, qsum, esl, dx
    INTEGER              :: iord = 3
    LOGICAL              :: found

    ! XMAP begins here!
    do i=1,im+1
       x1(i) = lon1(i)
    enddo

    do i=1,im
       dx1(i) = x1(i+1) - x1(i)
    enddo
    
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
             write(6,*) 'failed in xmap'
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
             write(6,*) 'failed in xmap'
             stop
          else
             dx1(i2-1) = dx1(i2-1-im)
             x1(i2) = x1(i2-1) + dx1(i2-1)
          endif
       endif
    enddo
    
    !write(6,*) 'i1,i2=',i1,i2

    do 1000 j=1,jm

       !=================================================================
       ! Area preserving mapping
       !================================================================

       ! Construct subgrid PP distribution
       call ppm_cycle(im, q1(1,j), al(1), ar(1), a6(1), qtmp(0), iord)
       
       ! check to see if ghosting is necessary
       ! Western edge
       if ( i1 .le. 0 ) then
          do i=i1,0
             qtmp(i) = qtmp(im+i)
             al(i) = al(im+i)
             ar(i) = ar(im+i)
             a6(i) = a6(im+i)
          enddo
       endif
       
       ! Eastern edge:
       if ( i2 .gt. im+1 ) then
          do i=im+1,i2-1
             qtmp(i) = qtmp(i-im)
             al(i) =   al(i-im)
             ar(i) =   ar(i-im)
             a6(i) =   a6(i-im)
          enddo
       endif
        
       i0 = i1
        
       do 555 i=1,in
       do 100 m=i0,i2-1

          !=============================================================  
          ! locate the western edge: lon2(i)
          !=============================================================  
          if(lon2(i) .ge. x1(m) .and. lon2(i) .le. x1(m+1)) then
             pl = (lon2(i)-x1(m)) / dx1(m)
             
             if(lon2(i+1) .le. x1(m+1)) then
                
                ! entire new grid is within the original grid
                pr = (lon2(i+1)-x1(m)) / dx1(m)
                q2(i,j) = al(m) + 0.5d0*(a6(m)+ar(m)-al(m)) &
               &                  *(pr+pl)-a6(m)*r3*(pr*(pr+pl)+pl**2)
                i0 = m
                goto 555
             else

                ! Left most fractional area
                qsum = (x1(m+1)-lon2(i))*(al(m)+0.5d0*(a6(m)+ &
               &              ar(m)-al(m))*(1.+pl)-a6(m)*   &
               &               (r3*(1.+pl*(1.+pl))))

                do mm=m+1,i2-1

                   ! locate the eastern edge: lon2(i+1)
                   if(lon2(i+1) .gt. x1(mm+1) ) then

                      ! Whole layer
                      qsum = qsum + dx1(mm)*qtmp(mm)

                   else
                      ! Right most fractional area
                      dx = lon2(i+1)-x1(mm)
                      esl = dx / dx1(mm)
                      qsum = qsum + dx*(al(mm)+0.5d0*esl* &
                     &              (ar(mm)-al(mm)+a6(mm)*(1.-r23*esl)))
                      i0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100       continue
123       q2(i,j) = qsum / ( lon2(i+1) - lon2(i) )
555    continue
1000 continue

   END SUBROUTINE xmap
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ppm_cycle
!
! !DESCRIPTION: PPM\_CYCLE is called by XMAP
!\\
!\\
! !INTERFACE:
!
   subroutine ppm_cycle(im, q, al, ar, a6, p, iord)
!
! !INPUT PARAMETERS:
! 
     INTEGER, INTENT(IN)  :: im, iord
     REAL*8,  INTENT(IN)  :: q(1)
! 
! !OUTPUT PARAMETERS:
!
     REAL*8,  INTENT(OUT) :: al(1), ar(1), a6(1), p(0:im+1)
!
! !REVISION HISTORY:
!  21 Sep 2000 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
     REAL*8               :: dm(0:im), tmp, qmax, qmin
     INTEGER              :: i, lmt
     REAL*8,  PARAMETER   :: r3 = 1./3. 

     ! PPM_CYCLE begins here!
     p(0) = q(im)
     do i=1,im
        p(i) = q(i)
     enddo
     p(im+1) = q(1)

     ! 2nd order slope
     do i=1,im
        tmp = 0.25d0*(p(i+1) - p(i-1))
        qmax = max(p(i-1), p(i), p(i+1)) - p(i)
        qmin = p(i) - min(p(i-1), p(i), p(i+1))
        dm(i) = sign(min(abs(tmp),qmax,qmin), tmp)
     enddo
     dm(0) = dm(im)

     do i=1,im
        al(i) = 0.5d0*(p(i-1)+p(i)) + (dm(i-1) - dm(i))*r3
     enddo

     do i=1,im-1
        ar(i) = al(i+1)
     enddo
     ar(im) = al(1)

     if(iord .le. 6) then
        do i=1,im
           a6(i) = 3d0*(p(i)+p(i)  - (al(i)+ar(i)))
        enddo
        lmt = iord - 3
        if(lmt.le.2) call lmppm(dm(1),a6(1),ar(1),al(1),p(1),im,lmt)
     else
        call huynh(im, ar(1), al(1), p(1), a6(1), dm(1))
     endif

   END SUBROUTINE ppm_cycle
!EOC
!-----------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lmppm
!
! !DESCRIPTION: Subroutine LMPPM is called by PPM\_CYCLE.
!\\
!\\ 
! !INTERFACE:
!
   SUBROUTINE lmppm(dm, a6, ar, al, p, im, lmt)
!
! !INPUT PARAMETERS:
!
     INTEGER           :: im, lmt
     REAL*8            :: p(im),dm(im)
!
! !INPUT/OUTPUT PARAMETERS:
!
     REAL*8            :: a6(im),ar(im),al(im)
!
! !REMARKS:
!  LMT = 0: full monotonicity
!  LMT = 1: semi-monotonic constraint (no undershoot)
!  LMT = 2: positive-definite constraint
!
! !REVISION HISTORY:
!  21 Sep 2000 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
     INTEGER           :: i

     REAL*8            :: da1, da2, fmin, a6da
     REAL*8, PARAMETER :: r12 = 1.d0/12.d0 

     ! LMPPM begins here!
     if(lmt.eq.0) then

        ! Full constraint
        do 100 i=1,im
           if(dm(i) .eq. 0.) then
              ar(i) = p(i)
              al(i) = p(i)
              a6(i) = 0.
           else
              da1  = ar(i) - al(i)
              da2  = da1**2
              a6da = a6(i)*da1
              if(a6da .lt. -da2) then
                 a6(i) = 3d0*(al(i)-p(i))
                 ar(i) = al(i) - a6(i)
              elseif(a6da .gt. da2) then
                 a6(i) = 3d0*(ar(i)-p(i))
                 al(i) = ar(i) - a6(i)
              endif
           endif
100     continue

     elseif(lmt.eq.1) then

        ! Semi-monotonic constraint
        do 150 i=1,im
           if(abs(ar(i)-al(i)) .ge. -a6(i)) go to 150
           if(p(i).lt.ar(i) .and. p(i).lt.al(i)) then
              ar(i) = p(i)
              al(i) = p(i)
              a6(i) = 0d0
           elseif(ar(i) .gt. al(i)) then
              a6(i) = 3d0*(al(i)-p(i))
              ar(i) = al(i) - a6(i)
           else
              a6(i) = 3d0*(ar(i)-p(i))
              al(i) = ar(i) - a6(i)
           endif
150     continue
           
     elseif(lmt.eq.2) then

        ! Positive definite constraint
        do 250 i=1,im
           if(abs(ar(i)-al(i)) .ge. -a6(i)) go to 250
           fmin = p(i) + 0.25d0*(ar(i)-al(i))**2/a6(i) + a6(i)*r12
           if(fmin.ge.0d0) go to 250
           if(p(i).lt.ar(i) .and. p(i).lt.al(i)) then
              ar(i) = p(i)
              al(i) = p(i)
              a6(i) = 0d0
           elseif(ar(i) .gt. al(i)) then
              a6(i) = 3d0*(al(i)-p(i))
              ar(i) = al(i) - a6(i)
           else
              a6(i) = 3d0*(ar(i)-p(i))
              al(i) = ar(i) - a6(i)
           endif
250     continue
     endif

   END SUBROUTINE lmppm
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: huynh
!
! !DESCRIPTION: Subroutine HUYNH enforces Huynh's 2nd constraint in 1D 
!  periodic domain
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE huynh(im, ar, al, p, d2, d1)
!
! !INPUT PARAMETERS:
!
    INTEGER :: im
    REAL*8  :: p(im)
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL*8  :: ar(im), al(im), d2(im), d1(im)
!
! !REVISION HISTORY:
!  21 Sep 2000 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: i
    REAL*8  :: pmp, lac, pmin, pmax

    !===================================================================
    ! HUYNH begins here!
    ! Compute d1 and d2
    !===================================================================
    d1(1) = p(1) - p(im)
    do i=2,im
       d1(i) = p(i) - p(i-1)
    enddo
    
    do i=1,im-1
       d2(i) = d1(i+1) - d1(i)
    enddo
    d2(im) = d1(1) - d1(im)

    !===================================================================
    ! Constraint for AR
    ! i = 1
    !===================================================================
    pmp   = p(1) + 2.0d0 * d1(1)
    lac   = p(1) + 0.5d0 * (d1(1)+d2(im)) + d2(im) 
    pmin  = min(p(1), pmp, lac)
    pmax  = max(p(1), pmp, lac)
    ar(1) = min(pmax, max(ar(1), pmin))
    
    do i=2, im
       pmp   = p(i) + 2.0d0*d1(i)
       lac   = p(i) + 0.5d0*(d1(i)+d2(i-1)) + d2(i-1)
       pmin  = min(p(i), pmp, lac)
       pmax  = max(p(i), pmp, lac)
       ar(i) = min(pmax, max(ar(i), pmin))
    enddo
     
    !==================================================================
    ! Constraint for AL
    !==================================================================
    do i=1, im-1
       pmp   = p(i) - 2.0d0*d1(i+1)
       lac   = p(i) + 0.5d0*(d2(i+1)-d1(i+1)) + d2(i+1)
       pmin  = min(p(i), pmp, lac)
       pmax  = max(p(i), pmp, lac)
       al(i) = min(pmax, max(al(i), pmin))
    enddo

    !==================================================================
    ! i=im
    !==================================================================
    i = im
    pmp    = p(im) - 2.0d0*d1(1)
    lac    = p(im) + 0.5d0*(d2(1)-d1(1)) + d2(1)
    pmin   = min(p(im), pmp, lac)
    pmax   = max(p(im), pmp, lac)
    al(im) = min(pmax, max(al(im), pmin))

    !==================================================================
    ! compute A6 (d2)
    !==================================================================
    do i=1, im
       d2(i) = 3d0*(p(i)+p(i)  - (al(i)+ar(i)))
    enddo

  END SUBROUTINE huynh
!EOC
!-----------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CasaRegridInit
!
! !DESCRIPTION: Subroutine CasaRegridInit initializes the longitude and 
!  latitude edge arrays for 0.5 x 0.666, 1 x 1.25, 2 x 2.5, and 4 x 5 grids.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CasaRegridInit
!
! !REMARKS:
!  Computation is done in REAL*8 and then casted to REAL*4 in order
!  to get correct values for the high-resolution grids. 
!
! !REVISION HISTORY:
!  09 Nov 2006- R. Yantosca - Initial version
!EOP
!-----------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J
    REAL*8  :: DI, DJ

    !-----------------------
    ! 1 x 1 GENERIC GRID
    ! edged on -180, -90
    !-----------------------

    ! Size of box
    DI = 1.00d0
    DJ = 1.00d0

    ! Lon edges
    DO I = 0, I1x1
       xedge_1x1(I+1) = -180d0 + ( DI * I )
    ENDDO

    ! Lat edges
    DO J = 0, J1x1
       yedge_1x1(J+1) =  -90d0 + ( DJ * J )
    ENDDO

    ! Reset poles
    yedge_1x1(1)      = -90d0
    yedge_1x1(J1x1+1) = +90d0

    ! Sine of latitude edges
    DO J = 1, J1x1+1
       sine_1x1(J)  = SIN( yedge_1x1(J) * D2R )
    ENDDO

    !-----------------------
    ! 2 x 2.5 GEOS Grid
    ! centered on -180
    !-----------------------

    ! Size of box
    DI = 2.5d0
    DJ = 2.0d0

    ! Lon edges
    DO I = 0, I2x25
       xedge_2x25(I+1)    = -180d0 - DI/2d0 + ( DI * I )
    ENDDO

    ! Lat edges
    DO J = 0, J2x25
       yedge_2x25(J+1)    =  -90d0 - DJ/2d0 + ( DJ * J )
    ENDDO

    ! Reset poles
    yedge_2x25(1)         = -90d0
    yedge_2x25(J2x25+1)   = +90d0

    ! Sine of latitude edges
    DO J = 1, J2x25+1
       sine_2x25(J)  = SIN( yedge_2x25(J) * D2R )
    ENDDO

    !-----------------------
    ! 4 x 5 GEOS Grid
    ! centerd on -180
    !-----------------------

    ! Size of box
    DI = 5d0
    DJ = 4d0

    ! Lon edges
    DO I = 0, I4x5
       xedge_4x5(I+1) = -180d0 - DI/2d0 + ( DI * I )
    ENDDO

    ! Lat edges and sine
    DO J = 0, J4x5
       yedge_4x5(J+1) =  -90d0 - DJ/2d0 + ( DJ * J )
    ENDDO

    ! Reset poles
    yedge_4x5(1)      = -90d0
    yedge_4x5(J4x5+1) = +90d0

    ! Sine of latitude edges
    DO J = 1, J4x5+1
       sine_4x5(J)  = SIN( yedge_4x5(J) * D2R )
    ENDDO

  END SUBROUTINE CasaRegridInit
!EOC
END MODULE CasaRegridModule

