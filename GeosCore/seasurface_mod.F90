!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: seasurface_mod.F90
!
! !DESCRIPTION: Module Sea\_Surface\_Mod is under development
!  to enhance the treatment of the ocean surface in GEOS-Chem 
!  with respect to dry deposition. This over time will expand
!  to a more comprehensive treatment of the ocean surface by
!  the model
!\\
!\\
! !INTERFACE: 
!
MODULE Sea_Surface_Mod

  USE CMN_SIZE_MOD
  USE ERROR_MOD                         ! Error checking routines
  USE GC_GRID_MOD                       ! Horizontal grid definition
  USE MAPPING_MOD                       ! Mapping weights & areas
  USE PhysConstants                     ! Physical constants
  USE PRECISION_MOD                     ! For GEOS-Chem Precision

!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: SeaSurface_Iodide

CONTAINS
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_olson_landmap
!
! !DESCRIPTION: Subroutine SEASURFACE\_IODIDE returns the  
! sea surface iodide concentration from the input file for
! a given month in an (I,J) dimension array
!\\
!\\
! !INTERFACE:
  SUBROUTINE SeaSurface_Iodide(am_I_Root,State_Met,Input_Opt)

    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
    USE m_netcdf_io_open
    USE m_netcdf_io_read
    USE m_netcdf_io_readattr
    USE m_netcdf_io_close

    IMPLICIT NONE

#   include "netcdf.inc"

    LOGICAL, INTENT(IN)           :: am_I_Root
    !INTEGER, INTENT(IN)          :: Month
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    TYPE(MetState), INTENT(INOUT) :: State_Met
    !REAL(kind=f8), POINTER, INTENT(OUT) :: Iodide_Conc(IIPAR,JJPAR,1:12)
    
    !local variables
    INTEGER                    :: IIodide,JIodide,ITime,fID, &
         I,J,SumBoxs,III,II,JJ,IG
    REAL(kind=f8)              :: D_LON,D_LAT,SumIodide,dxdy4,dxdy
    REAL*4                     :: xedgeC_w,xedgeC_e,yedgeC_s,yedgeC_n
    REAL*4                     :: xedge_w,xedge_e,yedge_s,yedge_n,mapWt
    REAL(kind=f8), DIMENSION(:), ALLOCATABLE :: Lat,Lon,Time
    REAL(kind=f8), DIMENSION(:,:,:), ALLOCATABLE :: Iodide_Map
    CHARACTER(LEN=255)         :: nc_path
    CHARACTER(LEN=255) :: v_name             ! netCDF variable name 
    CHARACTER(LEN=255) :: a_name             ! netCDF attribute name
    CHARACTER(LEN=255) :: a_val              ! netCDF attribute value
     
    ! Arrays for netCDF start and count values
    INTEGER            :: st1d(1), ct1d(1)   ! For 1D arrays    
    INTEGER            :: st3d(3), ct3d(3)   ! For 3D arrays
    INTEGER            :: Month

    ! Arrays on the Iodide concentration map NATIVE GRID
    INTEGER, DIMENSION(:), ALLOCATABLE     :: indLon, shiftLon
    REAL(kind=f8), DIMENSION(:), ALLOCATABLE :: lonedge, latedge

    !initialise
    !Iodide_Conc = 0.0_f8
    nc_path = TRIM(Input_Opt%IodideFile)

    IIodide = 2880
    JIodide = 1441
    ITime   = 12

    D_LON = 0.125_f8
    D_LAT = 0.125_f8

    ALLOCATE(Lat(JIodide))
    ALLOCATE(Lon(IIodide))
    ALLOCATE(Time(1:ITime))
    ALLOCATE(indLON(IIodide))
    ALLOCATE(shiftLon(IIodide))
    ALLOCATE(lonedge(IIodide+1))
    ALLOCATE(latedge(JIodide+1))
    ALLOCATE(Iodide_Map(IIodide,JIodide,ITime))

    !open file for reading
    CALL Ncop_Rd( fId, TRIM(nc_path) )

    ! Echo info to stdout
    IF ( am_I_Root ) THEN
       WRITE( 6, 100 ) REPEAT( '%', 79 )
       WRITE( 6, 110 ) TRIM(nc_path)
    ENDIF

    !----------------------------------------
    ! VARIABLE: lon
    !----------------------------------------
     
    ! Variable name
    v_name = "lon"
    
    ! Read lon from file
    st1d   = (/ 1       /)
    ct1d   = (/ IIodide /)
    CALL NcRd( Lon, fId, TRIM(v_name), st1d, ct1d )
    ! Echo info to stdout
    IF ( am_I_Root ) THEN
       WRITE( 6, 130 ) TRIM(v_name)!, TRIM(a_val) 
    ENDIF

    !----------------------------------------
    ! VARIABLE: lat
    !----------------------------------------
    
    ! Variable name
    v_name = "lat"
    
    ! Read lat from file
    st1d   = (/ 1       /)
    ct1d   = (/ JIodide /)
    CALL NcRd( Lat, fId, TRIM(v_name), st1d, ct1d )
    ! Echo info to stdout
    IF ( am_I_Root ) THEN
       WRITE( 6, 130 ) TRIM(v_name)!, TRIM(a_val) 
    ENDIF

    !----------------------------------------
    ! VARIABLE: time
    !----------------------------------------
    
    ! Variable name
    v_name = "time"
    
    ! Read lat from file
    st1d   = (/ 1       /)
    ct1d   = (/ ITime /)
    CALL NcRd( Time, fId, TRIM(v_name), st1d, ct1d )
    ! Echo info to stdout
    IF ( am_I_Root ) THEN
       WRITE( 6, 130 ) TRIM(v_name)!, TRIM(a_val) 
    ENDIF

    !----------------------------------------
    ! VARIABLE: iodide
    !----------------------------------------
    
    ! Variable name
    v_name = "Ensemble Monthly mean"
    
    ! Read lat from file
    st3d   = (/ 1,       1,       1     /)
    ct3d   = (/ IIodide, JIodide, ITime /)
    CALL NcRd( Iodide_Map, fId, TRIM(v_name), st3d, ct3d )
    ! Echo info to stdout
    IF ( am_I_Root ) THEN
       WRITE( 6, 130 ) TRIM(v_name)!, TRIM(a_val) 
    ENDIF

    !=================================================================
    ! Cleanup and close
    !=================================================================
    
    ! Close netCDF file
    CALL NcCl( fId )
    
    ! Echo info to stdout
    IF ( am_I_Root ) THEN
       WRITE( 6, 140 )
       WRITE( 6, 100 ) REPEAT( '%', 79 )
    ENDIF

    ! FORMAT statements
100 FORMAT( a                                              )
110 FORMAT( '%% Opening file  : ',         a               )
120 FORMAT( '%%  in directory : ',         a, / , '%%'     )
130 FORMAT( '%% Successfully read ',       a, ' [', a, ']' )
140 FORMAT( '%% Successfully closed file!'                 )

    ! Be lazy, construct lon edges from lon centers
    DO I = 1, IIodide
       lonedge(I)      = DBLE( lon(I) ) - ( D_LON * 0.5e+0_f8 )
       indLon(I)       = I
    ENDDO
    lonedge(IIodide+1) = lonedge(IIodide) + D_LON
    
    ! Be lazy, construct lat edges from lat centers
    DO J = 1, JIodide
       latedge(J)      = DBLE( lat(J) ) - ( D_LAT * 0.5e+0_f8 )
    ENDDO
    latedge(JIodide+1) = latedge(JIodide) + D_LAT
    
    ! Shift longitudes by 2 degrees to the west for date-line handling
    shiftLon           = CSHIFT( indLon, -20 )

    !$OMP PARALLEL DO                                                  &
    !$OMP DEFAULT( SHARED )                                            &
    !$OMP PRIVATE( I,        J,         xedgeC_w, yedgeC_s, xedgeC_e ) &
    !$OMP PRIVATE( yedgeC_n, dxdy4,     SumBoxs,  JJ,       III      ) &
    !$OMP PRIVATE( dxdy,     mapWt,     II,       xedge_w,  yedge_s  ) &
    !$OMP PRIVATE( xedge_e,  yedge_n,   SumIodide,IG                 ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO Month = 1, 12
       print*, 'averaging iodide for month ', Month 
       DO J = 1, JJPAR
          DO I = 1, IIPAR

             IG = I + I_LO - 1

             ! Edges of this GEOS-CHEM GRID box
             xedgeC_w  = GET_XEDGE( I,   J,   1 )          ! W edge
             yedgeC_s  = GET_YEDGE( I,   J,   1 )          ! S edge
             xedgeC_e  = GET_XEDGE( I+1, J,   1 )          ! E edge
             yedgeC_n  = GET_YEDGE( I,   J+1, 1 )          ! N edge
             
             ! "Area" of the GEOS-CHEM GRID box in degrees (DLON * DLAT)
             dxdy4     = ( xedgeC_e - xedgeC_w ) * ( yedgeC_n - yedgeC_s )
             
             SumIodide = 0.0_f8
             SumBoxs = 0

             ! Loop over latitudes on the NATIVE GRID
             DO JJ  = 1, JIodide

                ! Latitude edges of this NATIVE GRID box
                yedge_s    = latedge(JJ  )                ! S edge
                yedge_n    = latedge(JJ+1)                ! N edge
             
                DO III = 1, IIodide
                   ! Initialize
                   dxdy = 0.0_f8
                   mapWt = 0.0_f8

                   IF ( IG == 1 ) THEN
                      II = shiftLon(III)
                   ELSE
                      II = indLon(III)
                   ENDIF
                
                   ! Edges of this NATIVE GRID box
                   xedge_w    = lonedge(II  )                ! W edge
                   xedge_e    = lonedge(II+1)                ! E edge

                   IF ( IG == 1 .and. II >= shiftLon(1) )  THEN
                      xedge_w = xedge_w - 360.0_f8
                      xedge_e = xedge_e - 360.0_f8
                   ENDIF

                   ! "Area" of the NATIVE GRID BOX in degrees (DLON * DLAT)
                   dxdy = ( xedge_e - xedge_w )*( yedge_n - yedge_s )
                   
                   CALL GET_MAP_WT( xedge_w, xedge_e, xedgeC_w, xedgeC_e,  &
                        yedge_s, yedge_n, yedgeC_s, yedgeC_n,  &
                        mapWt                                 )

                   IF ( mapWt <= 0e0 .or. mapWt > 1e0 ) CYCLE

                   SumBoxs = SumBoxs + 1
                   SumIodide = SumIodide + Iodide_Map(II,JJ,Month)
                END DO
             ENDDO

             State_Met%Iodide_Conc(I,J,Month) = SumIodide/REAL(SumBoxs,f8)

          END DO
       END DO
    END DO
    !$OMP END PARALLEL DO
    print*, 'Completed averaging for Iodide'

    IF ( ALLOCATED( Lon   ) ) DEALLOCATE( Lon   )
    IF ( ALLOCATED( Lat   ) ) DEALLOCATE( Lat   )
    IF ( ALLOCATED( Time  ) ) DEALLOCATE( Time  )
    IF ( ALLOCATED( Iodide_Map ) ) DEALLOCATE( Iodide_Map  )
    IF ( ALLOCATED( indLON ) ) DEALLOCATE( indLON )
    IF ( ALLOCATED( shiftLon ) ) DEALLOCATE( shiftLon )
    IF ( ALLOCATED( lonedge ) ) DEALLOCATE( lonedge )
    IF ( ALLOCATED( latedge ) ) DEALLOCATE( latedge )


  END SUBROUTINE SeaSurface_Iodide


END MODULE Sea_Surface_Mod
