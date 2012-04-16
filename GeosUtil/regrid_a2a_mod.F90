!------------------------------------------------------------------------------
!          Matthew Cooper - Dalhousie University                              !
!------------------------------------------------------------------------------
! !MODULE: regrid_a2a_mod.F90
!
! !DESCRIPTION: REGRID\_A2A\_MOD uses an algorithm adapted from MAP\_A2A code 
!  to regrid from one horizonatal grid to another
!\\
!\\
! !INTERFACE: 
!
MODULE REGRID_A2A_MOD
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE :: XMAP
  PRIVATE :: YMAP
  PRIVATE :: MAP_A2A
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: DO_REGRID_A2A

!------------------------------------------------------------------------------

  CONTAINS

  SUBROUTINE DO_REGRID_A2A(FILENAME,IM,JM,INGRID,OUTGRID,PERAREA)
!
!******************************************************************************
!  Subroutine DO_REGRID_A2A regrids 2-D data
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME            : Name of file with longitude and latitude edge info
!  (2 ) IM      (INTEGER  ) : Length of longitude edge vector
!  (3 ) JM      (INTEGER  ) : Length of latitude edge vector
!  (4 ) INGRID  (REAL*8   ) : Input grid
!  (5 ) PERAREA (INTEGER  ) : =1 if need to convert ingrid to per area
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) OUTGRID (REAL*8   ) : Output data array
!******************************************************************************
    ! Uses
    USE GRID_MOD, ONLY : GET_XEDGE
    USE GRID_MOD, ONLY : GET_YSIN
    USE GRID_MOD, ONLY : GET_AREA_CM2
    USE FILE_MOD, ONLY : IOERROR, IU_REGRID
    USE CMN_SIZE_MOD
    USE CMN_GCTM_MOD

    ! Arguments
    INTEGER                        :: IM
    INTEGER                        :: JM
    INTEGER                        :: PERAREA
    REAL*8,           ALLOCATABLE  :: INGRID(:,:)
    REAL*8,           ALLOCATABLE  :: OUTGRID(:,:)


    CHARACTER(LEN=*)             :: FILENAME

    ! Local variables
    !INTEGER                      :: IG, IV
    INTEGER                       :: I, J
    INTEGER, PARAMETER            :: IN = IIPAR 
    INTEGER, PARAMETER            :: JN = JJPAR
    REAL*8, ALLOCATABLE           :: INSIN(:), INLON(:)
    REAL*8                        :: LON2(IN+1), SIN2(JN+1)
    CHARACTER(LEN=15)             :: HEADER1
    INTEGER                       :: IOS, M
    REAL*8                        :: INAREA,RLAT
    CHARACTER(LEN=20)             :: FMT_LAT, FMT_LON, FMT_LEN
 


    !======================================================================
    ! DO_REGRID_A2A BEGINS HERE
    !======================================================================

    ALLOCATE(INSIN(JM))
    ALLOCATE(INLON(IM+1))

    ! Get latitude and longitude edges for current model resolution
    DO I = 1, IIPAR+1
      LON2 (I) = GET_XEDGE(I)
    ENDDO
    DO J = 1, JJPAR+1
      SIN2 (J) = GET_YSIN(J)
    ENDDO  
  

      ! Open emission factor file (ASCII format)
      print*,'Filename'
      print*,TRIM(FILENAME)
      call flush(6)
      OPEN( IU_REGRID, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_REGRID, 'latlonread' )

      ! Read input grid lat/lon edges
        WRITE(FMT_LEN,*) IM+1
        FMT_LON='(' // TRIM ( FMT_LEN ) // 'F9.3)'
        WRITE(FMT_LEN,*) JM
        FMT_LAT='(' // TRIM ( FMT_LEN ) // 'F15.10)'

        READ( IU_REGRID, '(A15)',IOSTAT=IOS ) HEADER1
        READ( IU_REGRID,FMT_LON,IOSTAT=IOS )( INLON(M), M=1,IM+1 )
        READ( IU_REGRID,FMT_LAT,IOSTAT=IOS )( INSIN(M), M=1,JM )

      ! Close file
      CLOSE( IU_REGRID )

      !Convert input to per area units if necessary
       IF(PERAREA == 1) THEN
       DO J = 1, JM-1
          RLAT = INSIN(J+1) - INSIN(J)
          INAREA = 2d0*PI*Re*RLAT*1d4*Re/DBLE(IM)
          DO I = 1, IM
             INGRID(I,J)=INGRID(I,J)/INAREA
          ENDDO
       ENDDO
       ENDIF

    CALL MAP_A2A( IM, JM-1, INLON, INSIN, INGRID, &
   &       IN, JN, LON2, SIN2, OUTGRID, 0, 0)


       !Convert back from "per area" if necessary
       IF( PERAREA==1 ) THEN
       DO I=1, IN
          DO J=1, JN
             OUTGRID(I,J)=OUTGRID(I,J)*GET_AREA_CM2(J)
          ENDDO
       ENDDO
       ENDIF

     DEALLOCATE(INSIN)
     DEALLOCATE(INLON)

    END SUBROUTINE DO_REGRID_A2A

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!  (1 ) INLON   (REAL*8   ) : Longitude edges of input grid
!  (2 ) INSIN   (REAL*8   ) : Sine of input grid latitude edges
!  (3 ) INGRID  (REAL*8   ) : Data array to be regridded

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
!                   Prasad Kasibhatla - Duke University                       !
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: i, j0, m, mm, j
      REAL*8               :: dy1(jm)
      REAL*8               :: dy
      REAL*4               :: qsum, sum
    
    ! YMAP begins here!
      do j=1,jm-ig
         dy1(j) = sin1(j+1) - sin1(j)
      enddo

      !===============================================================
      ! Area preserving mapping
      !===============================================================

      do 1000 i=1,im
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
100      continue
123      q2(i,j) = qsum / ( sin2(j+1) - sin2(j) )
555      continue
1000  continue

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
!                   Prasad Kasibhatla - Duke University                       !
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
      ! Original E-W dimension
      INTEGER, INTENT(IN)  :: im           

      ! Target E-W dimension
      INTEGER, INTENT(IN)  :: in           
  
      ! Original N-S dimension
      INTEGER, INTENT(IN)  :: jm           
  
      ! Original western edge of the cell
      REAL*8,  INTENT(IN)  :: lon1(im+1)   
  
      ! Original data at center of the cell
      REAL*8,  INTENT(IN)  :: q1(im,jm)    
  
      ! Target cell's western edge
      REAL*8,  INTENT(IN)  :: lon2(in+1)   
!
! !OUTPUT PARAMETERS:
!
    ! Mapped data at the target resolution
      REAL*8,  INTENT(OUT) :: q2(in,jm)    
!
! !REMARKS:
!   lon1(1) < lon1(2) < lon1(3) < ... < lon1(im) < lon1(im+1)
!   lon2(1) < lon2(2) < lon2(3) < ... < lon2(in) < lon2(in+1)
!
! !AUTHOR:
!   Developer: Prasad Kasibhatla
!   March 6, 2012
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
      REAL*4               :: qsum
      LOGICAL              :: found
!

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
      
      do 1000 j=1,jm

         !=================================================================
         ! Area preserving mapping
         !================================================================
  
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
100         continue
123         q2(i,j) = qsum / ( lon2(i+1) - lon2(i) )
555      continue
1000  continue

   END SUBROUTINE xmap

END MODULE REGRID_A2A_MOD
