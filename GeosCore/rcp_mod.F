      MODULE RCP_MOD

        IMPLICIT NONE

        PRIVATE
        PUBLIC :: GET_NETCDF_VAR_1D, GET_NETCDF_VAR_3D, 
     &            GET_NETCDF_VAR_4D,
     &            READ_AND_REGRID_SURFACE,
     &            READ_AND_REGRID_VOLUME
        
        CHARACTER(LEN=*),PARAMETER :: fmtRcpFile = 
     &       '(I4,"/IPCC_emissions_",A5,"_",3(A,"_"),'//
     &       '"_0.5x0.5_v1_01_03_2010.nc")'

        CONTAINS

        !---------------------------------------------------------------

        SUBROUTINE EXPAND_RCP_VARNAME( varname, longname )

           CHARACTER(LEN=*), INTENT(IN)  :: varname
           CHARACTER(LEN=25),INTENT(OUT) :: longname

           SELECT CASE (varname)
              CASE ('emiss_ene')
                 longname='energy production'
              CASE ('emiss_dom')
                 longname='residential combustion'
              CASE ('emiss_ind')
                 longname='industrial combustion'
              CASE ('emiss_tra')
                 longname='land transport'
              CASE ('emiss_wst')
                 longname='waste'
              CASE ('emiss_agr')
                 longname='agriculture'
              CASE ('emiss_awb')
                 longname='ag waste burning'
              CASE ('emiss_slv')
                 longname='solvent use'
              CASE ('emiss_shp')
                 longname='shipping'
              CASE ('grassfire')
                 longname='grass fires'
              CASE ('forestfire')
                 longname='forest fires'
              CASE ('emiss_air')
                 longname='aircraft'
              CASE ('emiss_nrd')
                 longname='land non-road'
              CASE DEFAULT
                 longname=''
              END SELECT


        END SUBROUTINE EXPAND_RCP_VARNAME


        !---------------------------------------------------------------
        

        SUBROUTINE READ_AND_REGRID_1SURFACE( filename, varname, 
     &        outLonEdge, outLatEdge, nMoments, outEmission )

          ! Include
          USE NETCDF

          ! Arguments
          CHARACTER(LEN=*),INTENT(IN) :: filename, varname
          REAL*8,DIMENSION(:),INTENT(IN) :: outLonEdge, outLatEdge
          INTEGER,INTENT(IN) :: nMoments
          REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(OUT) ::
     &         outEmission

          ! Local variables
          REAL*8, ALLOCATABLE, DIMENSION(:) :: inLon, inLat, inTime
          REAL*8, ALLOCATABLE, DIMENSION(:) :: inLonEdge, inLatEdge
          REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: inEmission
          REAL*8, ALLOCATABLE, DIMENSION(:,:) :: boxArea
          
          INTEGER            :: status, ncid, nDims, nVars, nAtts
          INTEGER            :: varID,  len
          INTEGER            :: nLon,   nLat, nLev, nTime
          INTEGER            :: nOutLon,      nOutLat
          INTEGER            :: I, J, T

          CHARACTER(LEN=25)  :: longname
          CHARACTER(LEN=25)  :: units !CDH caution using fixed length
          
          REAL*8        :: totIn, totOut
          REAL*8        :: molWt 

          ! Constants
          real*8, parameter ::  A0 = 6371000.d0  !radius of earth, m
          real*8, parameter ::  CPI    = 3.141592653589793d0
          real*8, parameter ::  CPI180 = CPI/180.D0
          real*8, parameter ::  !days per month, no leapyear
     &         dpm(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)

          !-----------------------
          ! Begin Subroutine
          !-----------------------
          
          !-----------------------
          ! Read data from netCDF
          !-----------------------
            
          ! Latitude and longitude
          CALL GET_NETCDF_VAR_1D( filename, 'lon',  inLon  )
          CALL GET_NETCDF_VAR_1D( filename, 'lat',  inLat  )
          CALL GET_NETCDF_VAR_1D( filename, 'time', inTime )
         
          ! Dimensions of inputs
          nLon  = SIZE( inLon  )
          nLat  = SIZE( inLat  )
          nTime = SIZE( inTime )

          ALLOCATE( inEmission( nLon, nLat, nTime ) ) 
          inEmission = 0d0

          ! Emission variable, expected units: kg/m2/s for each month
          CALL GET_NETCDF_VAR_3D( filename, varname, inEmission )
          
          ! Check the units for varname
          CALL GET_NETCDF_ATT_CHAR( filename, varname, 'units', units )

          ! Make sure units conform to expectations, otherwise convert
          SELECT CASE ( TRIM( units ) )
             CASE ( 'kg m-2 s-1', 'kg m-2 sec-1', 'kg/m2/s' )
                !Do nothing, this is what we expect
             CASE DEFAULT
             WRITE(6,'(4(A,X))') 'Unrecognized units ',
     &            TRIM(filename), TRIM(varname), TRIM(units)
             STOP
          END SELECT

          ! Check the molecular weight for varname
          CALL GET_NETCDF_ATT_FLOAT( filename, varname, 
     &         'molecular_weight', molWt )

          !--------------------------
          ! Create output array
          !--------------------------

          ! Size of output array
          nOutLon = SIZE( outLonEdge ) - 1
          nOutLat = SIZE( outLatEdge ) - 1

          ! Allocate output array
          if ( ALLOCATED( outEmission ) ) DEALLOCATE( outEmission )
          ALLOCATE( outEmission( nOutLon, nOutLat, nMoments, nTime ) )
          outEmission = 0d0


          !--------------------------
          ! Processing and regridding
          !--------------------------

          ! Edges of input grid
          ALLOCATE( inLonEdge(nLon+1) ) 
          ALLOCATE( inLatEdge(nLat+1) )
          inLonEdge=0d0
          inLatEdge=0d0
          
          ! Calculate edges, assumes that longitudes increase monotonically
          DO I=2, nLon
             IF ( inLon(I-1) > inLon(I) ) 
     &            STOP 'longitudes must increase monotonically!'
             inLonEdge(I) = ( inLon(I-1) + inLon(I) ) / 2d0
          ENDDO
          inLonEdge(1)      = ( inLon(1) + inLon(nLon) - 360d0 ) / 2d0
          inLonEdge(nLon+1) = inLonEdge(1) + 360d0

          DO J=2, nLat
             IF ( inLat(J-1) > inLat(J) ) 
     &            STOP 'latitudes must increase monotonically!'
             inLatEdge(J) = ( inLat(J-1) + inLat(J) ) / 2d0
          ENDDO
          inLatEdge(1)      = -90d0
          inLatEdge(nLat+1) = 90d0


          ! Calculate grid box areas, m^2
          ALLOCATE( boxArea(nLon,nLat) )
          boxArea = 0d0
          DO J=1, nLat
          DO I=1, nLon
             boxArea(I,J) = A0*A0 * CPI180*(inLonEdge(I+1)-inLonEdge(I)) 
     &         * (sin(CPI180*inLatEdge(J+1)) - sin(CPI180*inLatEdge(J)))
          ENDDO
          ENDDO

          ! Convert kg/m2/s -> kg/s for each month
          DO T=1, nTime
             inEmission(:,:,T) = inEmission(:,:,T) * boxArea(:,:)             
          ENDDO

          ! Convert kg/s -> kg/y for each month
          inEmission = inEmission * (365d0*86400d0)

          ! total output emissions
          totOut = 0d0

          ! Loop over time
          DO T=1, nTime

             ! Regrid emissions to CTM grid, kg/y for each month
             CALL E_GRID( inEmission(:,:,T), 
     &            inLonEdge, inLatEdge, nLon, nLat, 
     &            outEmission(:,:,:,T), 
     &            outLonEdge, outLatEdge, nOutLon, nOutLat, 
     &            nMoments )
             
             ! Calculate total annual emissions, 
             totOut = totOut + 
     &            sum( outEmission(:,:,1,T) ) / 365d0 * dpm(T)

          ENDDO

          ! If there is only one time, then we calculate annualized emission
          IF ( nTime == 1 ) THEN
             totOut = totOut / dpm(1) * 365d0
          ENDIF

          ! Print annual emissions
          CALL EXPAND_RCP_VARNAME( varname, longname )
          print '(3X,A25,F8.2," Tg")', longname, totOut/1d9

          ! UNITS for GEOS-Chem
          ! Convert kg/y -> molec/cm2/s for each month
          FORALL (T=1:nTime, I=1:nMoments) 
             outEmission(:,:,I,T) = outEmission(:,:,I,T) / 
     &         (boxArea(:,:) * 1D4) / (molWt/1D3) * 6.02D23 / 
     &         (365D0 *86400D0)
          END FORALL

          inEmission = inEmission / (molWt/1D3) * 6.02D23 / 1D4


          !-----------------------
          ! Cleanup
          !-----------------------

          ! Deallocate all local variables
          IF ( ALLOCATED(inLon) ) DEALLOCATE(inLon)
          IF ( ALLOCATED(inLat) ) DEALLOCATE(inLat)
          IF ( ALLOCATED(inTime) ) DEALLOCATE(inTime)

          IF ( ALLOCATED(inLonEdge) ) DEALLOCATE(inLonEdge)
          IF ( ALLOCATED(inLatEdge) ) DEALLOCATE(inLatEdge)

          IF ( ALLOCATED(boxArea) ) DEALLOCATE(boxArea)
          IF ( ALLOCATED(inEmission) ) DEALLOCATE(inEmission)


        END SUBROUTINE READ_AND_REGRID_1SURFACE

        !---------------------------------------------------------------
        

        SUBROUTINE READ_AND_REGRID_1VOLUME( filename, varname, 
     &       outLonEdge, outLatEdge, outLevEdge, nMoments, outEmission )

          ! Include
          USE NETCDF

          ! Arguments
          CHARACTER(LEN=*),INTENT(IN) :: filename, varname
          REAL*8,DIMENSION(:),INTENT(IN) :: outLonEdge, outLatEdge
          REAL*8,DIMENSION(:),INTENT(IN) :: outLevEdge
          INTEGER,INTENT(IN) :: nMoments
          REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:,:),INTENT(OUT) ::
     &         outEmission

          ! Local variables
          REAL*8, ALLOCATABLE, DIMENSION(:) :: inLon,  inLat,  inTime
          REAL*8, ALLOCATABLE, DIMENSION(:) :: inLev,  inLevZ, inLevP
          REAL*8, ALLOCATABLE, DIMENSION(:) :: inLonEdge,  inLatEdge
          REAL*8, ALLOCATABLE, DIMENSION(:) :: inLevEdgeZ, inLevEdgeP
          REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: inEmission
          REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: outEmissionInLev
          REAL*8, ALLOCATABLE, DIMENSION(:,:) :: boxArea, Weights
          
          INTEGER            :: status, ncid, nDims, nVars, nAtts
          INTEGER            :: varID,  len
          INTEGER            :: nLon,   nLat, nLev, nTime
          INTEGER            :: nOutLon,      nOutLat,    nOutLev
          INTEGER            :: I, J, L, LL, M, T

          CHARACTER(LEN=25)  :: longname
          CHARACTER(LEN=25)  :: units !CDH caution using fixed length
          
          REAL*8             :: totIn, totOut, press, molWt

          LOGICAL            :: FOUND

          ! Constants
          real*8, parameter ::  A0 = 6371000.d0
          real*8, parameter ::  CPI    = 3.141592653589793d0
          real*8, parameter ::  CPI180 = CPI/180.D0
          real*8, parameter ::  !days per month, no leapyear
     &         dpm(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)


          real*8::t1,t2 !cdh temporary
          !-----------------------
          ! Begin Subroutine
          !-----------------------
          
          !-----------------------
          ! Read dimensions from netCDF
          !-----------------------
            
          ! Latitude and longitude
          CALL GET_NETCDF_VAR_1D( filename, 'lon',  inLon  )
          CALL GET_NETCDF_VAR_1D( filename, 'lat',  inLat  )
          CALL GET_NETCDF_VAR_1D( filename, 'lev',  inLev  )
          CALL GET_NETCDF_VAR_1D( filename, 'time', inTime )
         
          ! Dimensions of inputs
          nLon  = SIZE( inLon  )
          nLat  = SIZE( inLat  )
          nLev  = SIZE( inLev  )
          nTime = SIZE( inTime )


          !--------------------------
          ! Convert level altitudes to pressure, if necessary
          !--------------------------

          ! Edges of levels in geometric altitude (m) and pressure (hPa)
          ALLOCATE( inLevZ(nLev) )
          ALLOCATE( inLevP(nLev) )
          ALLOCATE( inLevEdgeZ(nLev+1) )
          ALLOCATE( inLevEdgeP(nLev+1) )
          inLevZ = 0d0 
          inLevP = 0d0 
          inLevEdgeZ = 0d0 
          inLevEdgeP = 0d0 

          ! Get the units for levels
          CALL GET_NETCDF_ATT_CHAR( filename, 'lev', 'units', units )


          SELECT CASE ( TRIM( units ) )
             CASE ( 'km' )
                ! Level centers in m
                ! Convert km -> m
                inLevZ = inLev * 1d3

                ! Calculate level edges in m
                DO L=2, nLev
                   IF ( inLevZ(L-1) > inLevZ(L) ) 
     &                  STOP 'altitude must decrease monotonically!'
                   inLevEdgeZ(L) = ( inLevZ(L-1) + inLevZ(L) ) / 2d0
                ENDDO
                inLevEdgeZ(1)      = 0d0
                inLevEdgeZ(nLev+1) = 2d0*inLevZ(nLev) - inLevEdgeZ(nLev)

                ! Calculate level centers in hPa
                CALL ALT2PRESSURE( inLevZ/1d3, inLevP )
                ! Calculate level edges in hPa
                CALL ALT2PRESSURE( inLevEdgeZ/1d3, inLevEdgeP )

             CASE ( 'hPa' )
                ! Do nothing, these are the desired units
                STOP 'Must add code to calculate level edges'
             CASE DEFAULT   
                WRITE(6,'(4(A,X),I)') 'Unrecognized units ',
     &            TRIM(filename), 'lev', TRIM(units),len_trim(units)
                STOP
          END SELECT

          ! Check the molecular weight for varname
          CALL GET_NETCDF_ATT_FLOAT( filename, varname, 
     &         'molecular_weight', molWt )

          !-----------------------
          ! Read volume data from netCDF
          !-----------------------

          ALLOCATE( inEmission( nLon, nLat, nLev, nTime ) ) 
          inEmission = 0d0

          ! Emission variable, expected units: kg/m3/s for each month
          CALL GET_NETCDF_VAR_4D( filename, varname, inEmission )
          
          ! Check the units for varname
          CALL GET_NETCDF_ATT_CHAR( filename, varname, 'units', units )

          ! Make sure units conform to expectations, otherwise convert
          ! Include other conversions here if necessary so that 
          ! final units are kg/m2/s
          SELECT CASE ( TRIM( units ) )
             CASE ( 'kg m-3 s-1', 'kg(NO)/m3/s' )
                ! Convert kg/m3/s -> kg/m2/s for each level
                DO L=1, nLev
                   inEmission(:,:,L,:) = inEmission(:,:,L,:) * 
     &                  (inLevEdgeZ(L+1)-inLevEdgeZ(L))
                ENDDO

             CASE DEFAULT
                WRITE(6,'(4(A,X))') 'Unrecognized units ',
     &            TRIM(filename), TRIM(varname), TRIM(units)
          END SELECT

          !--------------------------
          ! Create output array
          !--------------------------

          ! Size of output array
          nOutLon = SIZE( outLonEdge ) - 1
          nOutLat = SIZE( outLatEdge ) - 1
          nOutLev = SIZE( outLevEdge ) - 1

          ! Allocate output array
          if ( ALLOCATED( outEmission ) ) DEALLOCATE( outEmission )
          ALLOCATE( outEmission( nOutLon, nOutLat, nOutLev, 
     &         nMoments, nTime ) )
          outEmission = 0d0


          !--------------------------
          ! Calculate grid areas
          !--------------------------
          
          ! Edges of input grid
          ALLOCATE( inLonEdge(nLon+1) ) 
          ALLOCATE( inLatEdge(nLat+1) )
          inLonEdge=0d0
          inLatEdge=0d0
          
          ! Calculate edges, assumes that longitudes increase monotonically
          DO I=2, nLon
             IF ( inLon(I-1) > inLon(I) ) 
     &            STOP 'longitudes must increase monotonically!'
             inLonEdge(I) = ( inLon(I-1) + inLon(I) ) / 2d0
          ENDDO
          inLonEdge(1)      = ( inLon(1) + inLon(nLon) - 360d0 ) / 2d0
          inLonEdge(nLon+1) = inLonEdge(1) + 360d0

          DO J=2, nLat
             IF ( inLat(J-1) > inLat(J) ) 
     &            STOP 'latitudes must increase monotonically!'
             inLatEdge(J) = ( inLat(J-1) + inLat(J) ) / 2d0
          ENDDO
          inLatEdge(1)      = -90d0
          inLatEdge(nLat+1) = 90d0

          ! Calculate grid box areas, m^2
          ALLOCATE( boxArea(nLon,nLat) )
          boxArea = 0d0

          FORALL (I=1:nLon, J=1:nLat)
             boxArea(I,J) = A0*A0 * CPI180*(inLonEdge(I+1)-inLonEdge(I)) 
     &         * (sin(CPI180*inLatEdge(J+1)) - sin(CPI180*inLatEdge(J)))
          END FORALL

          !--------------------------
          ! Horizontal regridding
          !--------------------------

          ! Array with horizontal resolution of output grid
          ! and vertical resolution of input grid
          ALLOCATE( outEmissionInLev( nOutLon,  nOutLat, nLev, 
     &                                nMoments, nTime ) )
          outEmissionInLev = 0d0

          ! Convert kg/m2/s -> kg/s for each month
          FORALL (T=1:nTime, L=1:nLev )
             inEmission(:,:,L,T) = inEmission(:,:,L,T) *
     &         boxArea(:,:)              
          END FORALL

          ! Convert kg/s -> kg/y for each month
          inEmission = inEmission * (365d0*86400d0)

          ! total output emissions
          totOut = 0d0

          ! Loop over time, altitudes
          DO T=1, nTime
             DO L=1, nLev

                ! Regrid emissions to CTM grid, kg/y for each month
                CALL E_GRID( inEmission(:,:,L,T), 
     &               inLonEdge, inLatEdge, nLon, nLat, 
     &               outEmissionInLev(:,:,L,:,T), 
     &               outLonEdge, outLatEdge, nOutLon, nOutLat, 
     &               nMoments )
             
                ! Calculate total annual emissions, kg
                totOut = totOut + 
     &               sum( outEmission(:,:,L,1,T) ) / 365d0 * dpm(T)
             ENDDO

          ENDDO
          !--------------------------
          ! Vertical regridding
          !--------------------------

          ALLOCATE( Weights( nLev, nOutLev ) )
          Weights = 0d0

          ! Calculate regridding weights, based on pressure edges
          CALL Regrid_Column_Weights( inLevEdgeP, outLevEdge, Weights )
          

          ! Regrid column by column
          ! This forall statement generates a nonfatal compiler warning, but
          ! the code is correct. The error is: All active combinations of
          ! index-names are not used within the variable being defined
          ! (i.e., leftside) of this assignment-stmt.          
          forall( I=1:nOutLon, J=1:nOutLat, L=1:nLev, LL=1:nOutLev, 
     &         M=1:nMoments, T=1:nTime )
             outEmission(I,J,LL,M,T) = outEmission(I,J,LL,M,T) +
     &         outEmissionInLev(I,J,L,M,T) * Weights(L,LL)
          end forall


          if ( Allocated( Weights ) ) Deallocate( Weights )

          !-----------------------------------------------------------------
          ! OLD METHOD using nearest neighbor interpolation
          !-----------------------------------------------------------------
          !! Locate the level of the output (CTM) grid in which emissions occur
          !! This is based on pressure levels
          !! kg/y
          !DO L=1, nLev
          !   FOUND  = .FALSE.
          !   
          !   ! Loop over levels of CTM grid
          !   DO LL=1, nOutLev
          !      IF ( (inLevP(L) <= outLevEdge(LL)) .AND. 
     &    !           (inLevP(L) > outLevEdge(LL+1)) ) THEN
          !
          !         ! Emission on input level L occur on CTM level LL
          !         outEmission(:,:,LL,:,:) = 
     &    !              outEmission(:,:,LL,:,:) + 
     &    !              outEmissionInLev(:,:,L,:,:)
          !         FOUND = .TRUE.
          !         EXIT
          !     ENDIF
          !   ENDDO
          !   
          !   ! Error check if the correct level was found
          !   IF (FOUND == .FALSE.) THEN
          !      PRINT*, 'COULD NOT FIND LEVEL FOR EMISSIONS', L
          !      STOP
          !   ENDIF
          !ENDDO

          !--------------------------
          ! Calculate emissions total on new vertical grid
          !--------------------------

          totOut = 0d0

          DO T=1, nTime

             ! Total, kg
             totOut = totOut + dpm(T) / (365d0) *
     &            sum( outEmission(:,:,:,1,T) )

          ENDDO
          
          ! Print annual emissions
          CALL EXPAND_RCP_VARNAME( varname, longname )
!          print '(3X,A25,F8.2," Tg")', longname, totOut/1d9

          ! UNITS for GEOS-Chem
          ! Convert kg/y -> molec/s for each month
          outEmission = outEmission / (molWt/1D3) * 6.02D23 / 
     &         (365D0 *86400D0)

          !-----------------------
          ! Cleanup
          !-----------------------

          ! Deallocate all local variables
          IF ( ALLOCATED(inLon) ) DEALLOCATE(inLon)
          IF ( ALLOCATED(inLat) ) DEALLOCATE(inLat)
          IF ( ALLOCATED(inLev) ) DEALLOCATE(inLev)
          IF ( ALLOCATED(inLevZ) ) DEALLOCATE(inLevZ)
          IF ( ALLOCATED(inLevP) ) DEALLOCATE(inLevP)
          IF ( ALLOCATED(inTime) ) DEALLOCATE(inTime)

          IF ( ALLOCATED(inLonEdge) )  DEALLOCATE(inLonEdge)
          IF ( ALLOCATED(inLatEdge) )  DEALLOCATE(inLatEdge)
          IF ( ALLOCATED(inLevEdgeP) ) DEALLOCATE(inLevEdgeP)
          IF ( ALLOCATED(inLevEdgeZ) ) DEALLOCATE(inLevEdgeZ)

          IF ( ALLOCATED(boxArea) ) DEALLOCATE(boxArea)
          IF ( ALLOCATED(inEmission) ) DEALLOCATE(inEmission)
          IF ( ALLOCATED(outEmissionInLev) ) 
     &         DEALLOCATE(outEmissionInLev)


        END SUBROUTINE READ_AND_REGRID_1VOLUME

        !-------------------------------------------------------------------

        SUBROUTINE ADD_ALLOC_ARRAYS_4D( arrayInOut, array2 )
 
          !****************************************************
          ! Subroutine adds the contents of two array and
          ! returns the sum as the first argument. This
          ! subroutine accepts an unallocated array as the first
          ! argument and initializes it to zero. 
          !
          !****************************************************
          ! INPUT ARGUMENTS:
          !    arrayInOut
          !    array2
          !
          ! OUTPUT ARGUMENTS:
          !    arrayInOut
          !****************************************************
       
           ! Arguments
           REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(INOUT) ::
     &       arrayInOut
           REAL*8, DIMENSION(:,:,:,:), INTENT(IN) :: array2

           ! Local variables
           INTEGER :: sz(4)

           
           ! Allocate output array if necessary, set it to zero
           IF (.NOT. ALLOCATED( arrayInOut ) ) THEN
              sz = SHAPE( array2 )
              ALLOCATE(arrayInOut(sz(1),sz(2),sz(3),sz(4)))
              arrayInOut = 0d0
           ENDIF

           ! Add second array to first array
           arrayInOut = arrayInOut + array2
           

        END SUBROUTINE ADD_ALLOC_ARRAYS_4D


        !-------------------------------------------------------------------

        SUBROUTINE ADD_ALLOC_ARRAYS_5D( arrayInOut, array2 )
 
          !****************************************************
          ! Subroutine adds the contents of two array and
          ! returns the sum as the first argument. This
          ! subroutine accepts an unallocated array as the first
          ! argument and initializes it to zero. 
          !
          !****************************************************
          ! INPUT ARGUMENTS:
          !    arrayInOut
          !    array2
          !
          ! OUTPUT ARGUMENTS:
          !    arrayInOut
          !****************************************************
       
           ! Arguments
           REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE, INTENT(INOUT) ::
     &       arrayInOut
           REAL*8, DIMENSION(:,:,:,:,:), INTENT(IN) :: array2

           ! Local variables
           INTEGER :: sz(5)

           
           ! Allocate output array if necessary, set it to zero
           IF (.NOT. ALLOCATED( arrayInOut ) ) THEN
              sz = SHAPE( array2 )
              ALLOCATE(arrayInOut(sz(1),sz(2),sz(3),sz(4),sz(5)))
              arrayInOut = 0d0
           ENDIF

           ! Add second array to first array
           arrayInOut = arrayInOut + array2
           

        END SUBROUTINE ADD_ALLOC_ARRAYS_5D


        !-------------------------------------------------------------------

        SUBROUTINE READ_AND_REGRID_SURFACE( filename, 
     &        varname1, varname2, varname3, varname4,  
     &        varname5, varname6, varname7, varname8, 
     &        outLonEdge, outLatEdge, nMoments,       
     &        outEmission )
          
          ! Include
          USE NETCDF
          
          ! Arguments
          CHARACTER(LEN=*),INTENT(IN) :: filename
          CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: 
     &         varname1, varname2, varname3, varname4, 
     &         varname5, varname6, varname7, varname8
          REAL*8,DIMENSION(:),INTENT(IN) :: outLonEdge, outLatEdge
          INTEGER,INTENT(IN) :: nMoments
          REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(OUT) :: 
     &         outEmission
          
          ! Local variables
          REAL*8, ALLOCATABLE, DIMENSION(:) :: inLon, inLat, inTime
          REAL*8, ALLOCATABLE, DIMENSION(:) :: inLonEdge, inLatEdge
          REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: emis1, emis2, 
     &          emis3, emis4, emis5, emis6, emis7, emis8
          
          INTEGER            :: status, ncid, nDims, nVars, nAtts
          INTEGER            :: varID,  len
          INTEGER            :: nLon,   nLat, nLev, nTime
          INTEGER            :: nOutLon,      nOutLat
          INTEGER            :: T
          
          !-----------------------
          ! Begin Subroutine
          !-----------------------
            
          IF ( PRESENT(varname1) ) THEN
             CALL READ_AND_REGRID_1SURFACE( filename, varname1, 
     &             outLonEdge, outLatEdge, nMoments, emis1 )
             CALL ADD_ALLOC_ARRAYS_4D( outEmission, emis1 )
          ENDIF

          IF ( PRESENT(varname2) ) THEN
             CALL READ_AND_REGRID_1SURFACE( filename, varname2, 
     &             outLonEdge, outLatEdge, nMoments, emis2 )
             CALL ADD_ALLOC_ARRAYS_4D( outEmission, emis2 )
          ENDIF

          IF ( PRESENT(varname3) ) THEN
             CALL READ_AND_REGRID_1SURFACE( filename, varname3, 
     &             outLonEdge, outLatEdge, nMoments, emis3 )
             CALL ADD_ALLOC_ARRAYS_4D( outEmission, emis3 )
          ENDIF

          IF ( PRESENT(varname4) ) THEN
             CALL READ_AND_REGRID_1SURFACE( filename, varname4, 
     &             outLonEdge, outLatEdge, nMoments, emis4 )
             CALL ADD_ALLOC_ARRAYS_4D( outEmission, emis4 )
          ENDIF

          IF ( PRESENT(varname5) ) THEN
             CALL READ_AND_REGRID_1SURFACE( filename, varname5, 
     &             outLonEdge, outLatEdge, nMoments, emis5 )
             CALL ADD_ALLOC_ARRAYS_4D( outEmission, emis5 )
          ENDIF

          IF ( PRESENT(varname6) ) THEN
             CALL READ_AND_REGRID_1SURFACE( filename, varname6, 
     &             outLonEdge, outLatEdge, nMoments, emis6 )
             CALL ADD_ALLOC_ARRAYS_4D( outEmission, emis6 )
          ENDIF

          IF ( PRESENT(varname7) ) THEN
             CALL READ_AND_REGRID_1SURFACE( filename, varname7, 
     &             outLonEdge, outLatEdge, nMoments, emis7 )
             CALL ADD_ALLOC_ARRAYS_4D( outEmission, emis7 )
          ENDIF

          IF ( PRESENT(varname8) ) THEN
             CALL READ_AND_REGRID_1SURFACE( filename, varname8, 
     &             outLonEdge, outLatEdge, nMoments, emis8 )
             CALL ADD_ALLOC_ARRAYS_4D( outEmission, emis8 )
          ENDIF
            
          !-----------------------
          ! Cleanup
          !-----------------------

          ! Deallocate all local variables
          IF ( ALLOCATED(inLon) ) DEALLOCATE(inLon)
          IF ( ALLOCATED(inLat) ) DEALLOCATE(inLat)

          IF ( ALLOCATED(inLonEdge) ) DEALLOCATE(inLonEdge)
          IF ( ALLOCATED(inLatEdge) ) DEALLOCATE(inLatEdge)

          IF ( ALLOCATED(emis1) ) DEALLOCATE(emis1)
          IF ( ALLOCATED(emis2) ) DEALLOCATE(emis2)
          IF ( ALLOCATED(emis3) ) DEALLOCATE(emis3)
          IF ( ALLOCATED(emis4) ) DEALLOCATE(emis4)
          IF ( ALLOCATED(emis5) ) DEALLOCATE(emis5)
          IF ( ALLOCATED(emis6) ) DEALLOCATE(emis6)
          IF ( ALLOCATED(emis7) ) DEALLOCATE(emis7)
          IF ( ALLOCATED(emis8) ) DEALLOCATE(emis8)

          RETURN
        END SUBROUTINE READ_AND_REGRID_SURFACE
       
        !-------------------------------------------------------------------

        SUBROUTINE READ_AND_REGRID_VOLUME( filename, 
     &        varname1, varname2, varname3, varname4,  
     &        outLonEdge, outLatEdge, outLevEdge, nMoments,       
     &        outEmission )
          
          ! Include
          USE NETCDF
          
          ! Arguments
          CHARACTER(LEN=*),INTENT(IN) :: filename
          CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: 
     &         varname1, varname2, varname3, varname4
          REAL*8,DIMENSION(:),INTENT(IN) :: outLonEdge, outLatEdge 
          REAL*8,DIMENSION(:),INTENT(IN) :: outLevEdge
          INTEGER,INTENT(IN) :: nMoments
          REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:,:),INTENT(OUT) :: 
     &         outEmission
          
          ! Local variables
          REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: emis1, emis2, 
     &          emis3, emis4 
          
          !-----------------------
          ! Begin Subroutine
          !-----------------------
            
          IF ( PRESENT(varname1) ) THEN
             CALL READ_AND_REGRID_1VOLUME( filename, varname1, 
     &             outLonEdge, outLatEdge, outLevEdge, nMoments, emis1 )
             CALL ADD_ALLOC_ARRAYS_5D( outEmission, emis1 )
          ENDIF

          IF ( PRESENT(varname2) ) THEN
!             CALL READ_AND_REGRID_1VOLUME( filename, varname2, 
!     &             outLonEdge, outLatEdge, nMoments, emis2 )
             CALL ADD_ALLOC_ARRAYS_5D( outEmission, emis2 )
          ENDIF

          IF ( PRESENT(varname3) ) THEN
!             CALL READ_AND_REGRID_1VOLUME( filename, varname3, 
!     &             outLonEdge, outLatEdge, nMoments, emis3 )
             CALL ADD_ALLOC_ARRAYS_5D( outEmission, emis3 )
          ENDIF

          IF ( PRESENT(varname4) ) THEN
!             CALL READ_AND_REGRID_1VOLUME( filename, varname4, 
!     &             outLonEdge, outLatEdge, nMoments, emis4 )
             CALL ADD_ALLOC_ARRAYS_5D( outEmission, emis4 )
          ENDIF

          !-----------------------
          ! Cleanup
          !-----------------------

          ! Deallocate all local variables
          IF ( ALLOCATED(emis1) ) DEALLOCATE(emis1)
          IF ( ALLOCATED(emis2) ) DEALLOCATE(emis2)
          IF ( ALLOCATED(emis3) ) DEALLOCATE(emis3)
          IF ( ALLOCATED(emis4) ) DEALLOCATE(emis4)

          RETURN
        END SUBROUTINE READ_AND_REGRID_VOLUME
       


        !-------------------------------------------------------------------

        SUBROUTINE GET_NETCDF_ATT_CHAR( filename, varname, attname, 
     &     data )

          !****************************************************
          ! Subroutine retrieves attribute data from a netCDF file
          ! The variable and attribute names must be known in advance, but
          ! the size of the variable is unnecessary.
          ! The subroutine allocates the variable and returns 
          ! the value from file.
          !
          !****************************************************
          ! INPUT ARGUMENTS:
          !    filename
          !    varname
          !    attname
          !
          ! OUTPUT ARGUMENTS:
          !    data
          !****************************************************

          USE NETCDF

          ! Arguments
          CHARACTER(LEN=*),INTENT(IN)   :: filename, varname, attname
          CHARACTER(LEN=*),INTENT(OUT)  :: data
          
          ! Local variables
          INTEGER :: status, ncID, varID, dimLen

          CHARACTER(LEN=*),PARAMETER :: subr='get_netcdf_att_char: '

          ! Open file
          status = NF90_OPEN( TRIM(filename), NF90_NoWrite, ncid )
          IF( status /= NF90_NOERR ) 
     &         CALL HANDLE_ERR_LOC(status,subr // 'open')
          
          ! Locate variable
          status = NF90_INQ_VARID( ncid, TRIM(varname), varID )
          IF( status /= NF90_NOERR ) 
     &         CALL HANDLE_ERR_LOC(status,subr // 'inq_varid')

          ! Read the attribute
          status = NF90_GET_ATT( ncid, varID, attname, data )
          IF( status /= NF90_NOERR ) 
     &         CALL HANDLE_ERR_LOC(status,subr // 'get_att')

          ! Trim character string
          data = TRIM( data )

          ! Close
          status = NF90_CLOSE( ncid )
          IF( status /= NF90_NOERR ) 
     &         CALL HANDLE_ERR_LOC(status, subr // 'close')


        END SUBROUTINE GET_NETCDF_ATT_CHAR

        !-------------------------------------------------------------------

        SUBROUTINE GET_NETCDF_ATT_FLOAT( filename, varname, attname, 
     &     data )

          !****************************************************
          ! Subroutine retrieves attribute data from a netCDF file
          ! The variable and attribute names must be known in advance, but
          ! the size of the variable is unnecessary.
          ! The subroutine allocates the variable and returns 
          ! the value from file.
          !
          !****************************************************
          ! INPUT ARGUMENTS:
          !    filename
          !    varname
          !    attname
          !
          ! OUTPUT ARGUMENTS:
          !    data
          !****************************************************

          USE NETCDF

          ! Arguments
          CHARACTER(LEN=*),INTENT(IN)   :: filename, varname, attname
          REAL*8,INTENT(OUT)            :: data
          
          ! Local variables
          INTEGER :: status, ncID, varID, dimLen

          CHARACTER(LEN=*),PARAMETER :: subr='get_netcdf_att_float: '

          ! Open file
          status = NF90_OPEN( TRIM(filename), NF90_NoWrite, ncid )
          IF( status /= NF90_NOERR ) 
     &         CALL HANDLE_ERR_LOC(status,subr // 'open')
          
          ! Locate variable
          status = NF90_INQ_VARID( ncid, TRIM(varname), varID )
          IF( status /= NF90_NOERR ) 
     &         CALL HANDLE_ERR_LOC(status,subr // 'inq_varid')

          ! Read the attribute
          status = NF90_GET_ATT( ncid, varID, attname, data )
          IF( status /= NF90_NOERR ) 
     &         CALL HANDLE_ERR_LOC(status,subr // 'get_att')

          ! Close
          status = NF90_CLOSE( ncid )
          IF( status /= NF90_NOERR ) 
     &         CALL HANDLE_ERR_LOC(status, subr // 'close')


        END SUBROUTINE GET_NETCDF_ATT_FLOAT

        !-------------------------------------------------------------------

        SUBROUTINE GET_NETCDF_VAR_1D( filename, varname, data )

          !****************************************************
          ! Subroutine  retrieves 1-D data from a netCDF file
          ! The variable name must be known in advance, but
          ! the size of the variable is unnecessary.
          ! The subroutine allocates the variable and returns 
          ! the value from file.
          !
          !****************************************************
          ! INPUT ARGUMENTS:
          !    filename
          !    varname
          !
          ! OUTPUT ARGUMENTS:
          !    data
          !****************************************************

          USE NETCDF

          ! Arguments
          CHARACTER(LEN=*),INTENT(IN)   :: filename, varname
          REAL*8,ALLOCATABLE,INTENT(OUT):: data(:)
          
          ! Local variables
          INTEGER :: status, ncID, varID, dimLen
          INTEGER :: dimIDs(1)

          CHARACTER(LEN=*),PARAMETER :: subr='get_netcdf_var_1d: '

          ! Open file
          status = NF90_OPEN( TRIM(filename), NF90_NoWrite, ncid )
          IF( status /= NF90_NOERR ) 
     &         CALL HANDLE_ERR_LOC(status,subr // 'open')
          
          ! Locate variable
          status = NF90_INQ_VARID( ncid, TRIM(varname), varID )
          IF( status /= NF90_NOERR ) 
     &         CALL HANDLE_ERR_LOC(status,subr // 'inq_varid')

          ! Variable information
          status = NF90_INQUIRE_VARIABLE( ncid, varID, dimIDs=dimIDs )
          IF( status /= NF90_NOERR ) 
     &         CALL HANDLE_ERR_LOC(status,subr // 'variable')
          
          ! Dimension of the variable
          status = NF90_INQUIRE_DIMENSION( ncid, dimIDs(1), len=dimLen )
          IF( status /= NF90_NOERR ) 
     &         CALL HANDLE_ERR_LOC(status,subr // 'dimension 1')

          ! Allocate the variable
          ALLOCATE( data(dimLen), STAT=status )
          IF( status < 0 ) STOP 'data allocation failed'

          ! Read the variable
          status = NF90_GET_VAR( ncid, varID, data )
          IF( status /= NF90_NOERR ) 
     &         CALL HANDLE_ERR_LOC(status,subr // 'get_var')

          ! Close
          status = NF90_CLOSE( ncid )
          IF( status /= NF90_NOERR ) 
     &         CALL HANDLE_ERR_LOC(status, subr // 'close')


        END SUBROUTINE GET_NETCDF_VAR_1D

        !-------------------------------------------------------------------

        SUBROUTINE GET_NETCDF_VAR_3D( filename, varname, data )

          USE NETCDF

          ! Arguments
          CHARACTER(LEN=*),INTENT(IN)   :: filename, varname
          REAL*8,ALLOCATABLE,INTENT(OUT):: data(:,:,:)
          
          ! Local variables
          INTEGER :: status, ncID, varID
          INTEGER :: dimIDs(3), dimLen(3)

          CHARACTER(LEN=*),PARAMETER :: subr='get_netcdf_var_3d: '

          ! Open file
          status = NF90_OPEN( TRIM(filename), NF90_NoWrite, ncid )
          IF( status /= NF90_NOERR ) 
     &          CALL HANDLE_ERR_LOC(status,subr // 'open')
          
          ! Locate variable
          status = NF90_INQ_VARID( ncid, TRIM(varname), varID )
          IF( status /= NF90_NOERR ) 
     &          CALL HANDLE_ERR_LOC(status,subr // 'inq_varid')

          ! Variable information
          status = NF90_INQUIRE_VARIABLE( ncid, varID, dimIDs=dimIDs )
          IF( status /= NF90_NOERR ) 
     &          CALL HANDLE_ERR_LOC(status,subr // 'variable')
          
          ! 1st Dimension of the variable
          status = NF90_INQUIRE_DIMENSION(ncid,dimIDs(1),len=dimLen(1))
          IF( status /= NF90_NOERR ) 
     &          CALL HANDLE_ERR_LOC(status,subr // 'dimension 1')

          ! 2nd Dimension of the variable
          status = NF90_INQUIRE_DIMENSION(ncid,dimIDs(2),len=dimLen(2))
          IF( status /= NF90_NOERR ) 
     &          CALL HANDLE_ERR_LOC(status,subr // 'dimension 2')

          ! 3rd Dimension of the variable
          status = NF90_INQUIRE_DIMENSION(ncid,dimIDs(3),len=dimLen(3))
          IF( status /= NF90_NOERR ) 
     &          CALL HANDLE_ERR_LOC(status,subr // 'dimension 2')

          ! Allocate the variable
          ALLOCATE( data(dimLen(1),dimLen(2),dimLen(3)), STAT=status )
          IF( status < 0 ) STOP 'data allocation failed'

          ! Read the variable
          status = NF90_GET_VAR( ncid, varID, data )
          IF( status /= NF90_NOERR ) 
     &          CALL HANDLE_ERR_LOC(status, subr // 'get_var')

          ! Close
          status = NF90_CLOSE( ncid )
          IF( status /= NF90_NOERR ) 
     &          CALL HANDLE_ERR_LOC(status, subr // 'close')

          

        END SUBROUTINE GET_NETCDF_VAR_3D

        !-------------------------------------------------------------------

        SUBROUTINE GET_NETCDF_VAR_4D( filename, varname, data )

          USE NETCDF

          ! Arguments
          CHARACTER(LEN=*),INTENT(IN)   :: filename, varname
          REAL*8,ALLOCATABLE,INTENT(OUT):: data(:,:,:,:)
          
          ! Local variables
          INTEGER :: status, ncID, varID
          INTEGER :: dimIDs(4), dimLen(4)

          CHARACTER(LEN=*),PARAMETER :: subr='get_netcdf_var_3d: '

          ! Open file
          status = NF90_OPEN( TRIM(filename), NF90_NoWrite, ncid )
          IF( status /= NF90_NOERR ) 
     &          CALL HANDLE_ERR_LOC(status,subr // 'open')
          
          ! Locate variable
          status = NF90_INQ_VARID( ncid, TRIM(varname), varID )
          IF( status /= NF90_NOERR ) 
     &          CALL HANDLE_ERR_LOC(status,subr // 'inq_varid')

          ! Variable information
          status = NF90_INQUIRE_VARIABLE( ncid, varID, dimIDs=dimIDs )
          IF( status /= NF90_NOERR ) 
     &          CALL HANDLE_ERR_LOC(status,subr // 'variable')
          
          ! 1st Dimension of the variable
          status = NF90_INQUIRE_DIMENSION(ncid,dimIDs(1),len=dimLen(1))
          IF( status /= NF90_NOERR ) 
     &          CALL HANDLE_ERR_LOC(status,subr // 'dimension 1')

          ! 2nd Dimension of the variable
          status = NF90_INQUIRE_DIMENSION(ncid,dimIDs(2),len=dimLen(2))      if ( total( emiss_tot_shp )  gt 0d0 ) then $

          IF( status /= NF90_NOERR ) 
     &          CALL HANDLE_ERR_LOC(status,subr // 'dimension 2')

          ! 3rd Dimension of the variable
          status = NF90_INQUIRE_DIMENSION(ncid,dimIDs(3),len=dimLen(3))
          IF( status /= NF90_NOERR ) 
     &          CALL HANDLE_ERR_LOC(status,subr // 'dimension 2')

          ! 4th Dimension of the variable
          status = NF90_INQUIRE_DIMENSION(ncid,dimIDs(4),len=dimLen(4))
          IF( status /= NF90_NOERR ) 
     &          CALL HANDLE_ERR_LOC(status,subr // 'dimension 2')

          ! Allocate the variable
          ALLOCATE( data(dimLen(1),dimLen(2),dimLen(3),dimLen(4)), 
     &         STAT=status )
          IF( status < 0 ) STOP 'data allocation failed'

          ! Read the variable
          status = NF90_GET_VAR( ncid, varID, data )
          IF( status /= NF90_NOERR ) 
     &          CALL HANDLE_ERR_LOC(status, subr // 'get_var')

          ! Close
          status = NF90_CLOSE( ncid )
          IF( status /= NF90_NOERR ) 
     &          CALL HANDLE_ERR_LOC(status, subr // 'close')

          
        END SUBROUTINE GET_NETCDF_VAR_4D

        !-------------------------------------------------------------------

        SUBROUTINE handle_err_loc(status,location)

          use netcdf
          integer, intent(in)   :: status  !Error status
          character(len=*),intent(in) :: location
          if(status/=nf90_noerr)then
             print *, trim(nf90_strerror(status))
             print *, 'Occurred at ',location
             stop "stopped"
          endif
        END SUBROUTINE handle_err_loc


      !-----------------------------------------------------------------------
      ! The E_GRID subroutine is copied verbatim from p-grid.f in order
      ! to use the regridding routine within this module. The only change is
      ! to use STOP instead of subroutine EXITC.
      ! CDH 10/21/2010
      subroutine E_GRID(EBOX,XBEDG,YBEDG,IG,JG, EDXY,XDEDG,YDEDG,IM,JM, 
     &      NEDXY)
!-----------------------------------------------------------------------
!  generic regridder of 
!        emissions EBOX(1:IG,1:JG) 
!        on regular grid XBEDG(1:IG+1),YBEDG(1:JG+1) 
!                                
!  onto CTM grid   EDXY(1:IM,1:JM)
!        EBOX = emissions or absolute quantity (kg or kg/yr, not kg/m2)
!  IF original data is %-coverage or kg/m2, THEN must be pre-mult by area
!
!        X's Y's are the grid edges in degrees
!              note that the Longitude (X) edges can jump by 360 deg
!        XEDG(1,IM+1) is created = CTM longitude-grid edges (radians, monotonic)
!
!        EDXY(1:IM, 1:JM, 1:6) = emissions partitioned into box (I,J) + moments
!           EDXY(:,:,1) = total emissions
!           EDXY(:,:,2) = Sx  moment of emissions
!           EDXY(:,:,3) = Sxx moment of emissions
!           EDXY(:,:,4) = Sy  moment of emissions
!           EDXY(:,:,5) = Syy moment of emissions
!           EDXY(:,:,6) = Sxy moment of emissions

      implicit none
      integer, intent(in) :: IG,JG, IM,JM, NEDXY
      real*8,  intent(in) :: EBOX(IG,JG),XBEDG(IG+1),YBEDG(JG+1)
      real*8,  intent(in) :: XDEDG(IM+1),YDEDG(JM+1)
      real*8, intent(out) :: EDXY(IM,JM,NEDXY)

      real*8   Y1,Y2,SY1,SY2,DYBOX,YG,  YB1,YB2,SYB1,SYB2
      real*8   B_GYA(JG,JM), B_GYY(JG,JM)
      real*8   B_GXA(IG,IM), B_GXX(IG,IM)
      real*8   KX0,KXX,KY0,KYY,KXY, EMKIKJ
      real*8   SYG(JM+1), XEDG(IM+1)
      real*8   X1,X2,DXBOX,XG,  XB1,XB2
      integer  J,J1,J2,I,I1,I2, KI,KJ, I2M

      real*8, parameter ::  CPI    = 3.141592653589793d0
      real*8, parameter ::  C2PI   = 2.d0*CPI
      real*8, parameter ::  CPI180 = CPI/180.D0
!-----------------------------------------------------------------------

          XEDG(1) = CPI180*XDEDG(1)
      do I = 2,IM+1
        if (XDEDG(I) .gt. 0.d0)  then
          XEDG(I) = CPI180*XDEDG(I)
        else
          XEDG(I) = CPI180*(XDEDG(I)+360.d0)
        endif
      enddo
!     do I = 1,361
!       if (XBEDGE(I) .lt. 0.d0) then
!         XBEDG(I) = XBEDGE(I)+360.d0
!       else
!         XBEDG(I) = XBEDGE(I)
!       endif
!     enddo

          B_GXA(:,:) = 0.0d0
          B_GXX(:,:) = 0.0d0
          B_GYA(:,:) = 0.0d0
          B_GYY(:,:) = 0.0d0

!---for relative area use sin(LATITUDE) x LONGITUDE (degrees)
      do J = 1,JM+1
        SYG(J) = sin(YDEDG(J)*CPI180)
      enddo

!---setup the LATITUDE mapping and weights for the emission boxes
      do KJ = 1,JG
        YB1 = YBEDG(KJ)
        YB2 = YBEDG(KJ+1)
         SYB1 = sin(CPI180*YB1)
         SYB2 = sin(CPI180*YB2)
        do J = 2,JM
          if (YDEDG(J) .gt. YB1) then
             J1 = J-1
             goto 12
          endif
        enddo
             J1 = JM
   12   continue
        do J = J1+1,JM+1
          if (YDEDG(J) .ge. YB2) then
             J2 = J-1
             goto 14
          endif
        enddo
             J2 = JM
   14   continue
!---the emission box KJ=1:JG is contained in CTM grid boxes J1 to J2
!---    save the indices and fractions and weightings
        do J = J1,J2
          Y1 = max(YB1, YDEDG(J))
          Y2 = min(YB2, YDEDG(J+1))
           SY1 = sin(CPI180*Y1)
           SY2 = sin(CPI180*Y2)
           DYBOX = (SY2-SY1)/(SYB2-SYB1)     ! fractional area of emission box
!    mean value of Y in grid box
           YG = (0.5d0*(SY1+SY2)-SYG(J))/(SYG(J+1)-SYG(J))

          B_GYA(KJ,J) = DYBOX
          B_GYY(KJ,J) = YG

        enddo
      enddo

!---now setup the LONGITUDE mapping and weights for the emission boxes
!---    require XEDG(I=1:IM) to be monotonic, 
!---    but XBEDG(KI=1:IG) may have shift at dateline
      do KI = 1,IG
        XB1 = XBEDG(KI) * CPI180
        XB2 = XBEDG(KI+1) * CPI180
        if (XB1 .gt. XB2) XB2 = XB2 + C2PI
        if (XB1 .gt. XB2) STOP ' error#1 emission X-grid'
        if (XB1 .lt. XEDG(1)) then
          XB1 = XB1 + C2PI
          XB2 = XB2 + C2PI
          if (XB1 .lt. XEDG(1)) then
            write(6,'(I4,6F9.3)')KI,XBEDG(KI)*CPI180,XBEDG(KI+1)*CPI180, 
     &             XB1,XB2,XEDG(1)
            STOP ' error#2 emission X-grid'
          endif
        endif
        if (XB1 .gt. XEDG(IM+1)) then
          XB1 = XB1 - C2PI
          XB2 = XB2 - C2PI
          if (XB1 .gt. XEDG(IM+1)) STOP ' error#3 in em X-grid'
        endif
        do I = 2,IM
          if (XEDG(I) .gt. XB1) then
             I1 = I-1
             goto 16
          endif
        enddo
             I1 = IM
   16   continue
        do I = I1+1,IM+1
          if (XEDG(I) .ge. XB2) then
             I2 = I-1
             goto 18
          endif
        enddo
!---allow for XB2 to extend (wrap) beyond XEDG(IM+1)
        do I = 2,IM
          if (XEDG(I)+C2PI .ge. XB2) then
             I2 = IM + I-1
             goto 18
          endif
        enddo
             I2 = IM +IM
   18   continue
!---the emission box K=1:IG is contained in CTM grid boxes I1 to I2
!---    but the I2 may be > IM to wrap around the Longitude circle.
!---    save the indices and fractions and weightings
          I2M = min(I2, IM)
        do I = I1,I2M
          X1 = max(XB1, XEDG(I))
          X2 = min(XB2, XEDG(I+1))
          DXBOX = (X2-X1)/(XB2-XB1)         ! fractional area of emission box
          XG = (0.5d0*(X1+X2)-XEDG(I))/(XEDG(I+1)-XEDG(I)) ! mean X value 

          B_GXA(KI,I) = DXBOX
          B_GXX(KI,I) = XG

        enddo
!---while XB1 is constrained ot lie within 1:IM, XB2 may wrap around:
        IF (I2 .gt. IM) then
          I2M = I2 - IM
        do I = 1,I2M
          X1 = max(XB1, XEDG(I)+C2PI)
          X2 = min(XB2, XEDG(I+1)+C2PI)
          DXBOX = (X2-X1)/(XB2-XB1)
          XG = (0.5d0*(X1+X2)-XEDG(I)-C2PI)/(XEDG(I+1)-XEDG(I))

          B_GXA(KI,I) = DXBOX
          B_GXX(KI,I) = XG

        enddo
        endif
      enddo

!---sum contributions (& moments) to CTM grid (I,J) from emission boxes (KI,KJ)
!---SOM XY moments:
!         k00 = 1.d0
!         kX0 = XG - 0.5d0
!         kXX = XG*(XG - 1.d0) + C16TH
!         kY0 = YG - 0.5d0
!         kYY = YG*(YG - 1.d0) + C16TH
!         kXY = (XG - 0.5d0)*(YG - 0.5d0)
!--- for Kronecker-delat approx:  f(x,y) = DXBOX*DYBOX*EMI-BOX * k__(XG,YG)
!         So = INTEG[ f(x,y) ] x=0:1, y=0:1
!         Sx = 6 * INTEG[ f(x,y) * kX(x) ]
!         Sxy = 36 * INTEG[ f(x,y) * kXY(x,y) ]
!         Sxx = 30 * INTEG[ f(x,y) * kXX(x) ]
!--- adjust K_s to get So,Sx,Sxx,Sy,Syy,Sxy  (1:6)
!         K00 = 1.d0
!         KX0 = 6.d0*XG - 3.0d0
!         KXX = 30.d0*XG*(XG - 1.d0) + 5.d0
!         KY0 = 6.d0*YG - 3.0d0
!         KYY = 30.d0*YG*(YG - 1.d0) + 5.d0
!         KXY = KX0 * KY0

      if (NEDXY .eq. 6)  then
        do J = 1,JM
         do I = 1,IM
          EDXY(I,J,:) = 0.d0
           do KJ = 1,JG
            if (B_GYA(KJ,J) .gt. 1.d-5) then
             do KI = 1,IG
              if (B_GXA(KI,I) .gt. 1.d-5) then
                EMKIKJ = EBOX(KI,KJ) * B_GXA(KI,I) * B_GYA(KJ,J)
                  XG = B_GXX(KI,I)
                  YG = B_GYY(KJ,J)
                    KX0 = 6.d0*XG - 3.0d0
                    KXX = 30.d0*XG*(XG - 1.d0) + 5.d0
                    KY0 = 6.d0*YG - 3.0d0
                    KYY = 30.d0*YG*(YG - 1.d0) + 5.d0
                    KXY = KX0 * KY0
                EDXY(I,J,1) = EDXY(I,J,1) + EMKIKJ
                EDXY(I,J,2) = EDXY(I,J,2) + EMKIKJ*KX0
                EDXY(I,J,3) = EDXY(I,J,3) + EMKIKJ*KXX
                EDXY(I,J,4) = EDXY(I,J,4) + EMKIKJ*KY0
                EDXY(I,J,5) = EDXY(I,J,5) + EMKIKJ*KYY
                EDXY(I,J,6) = EDXY(I,J,6) + EMKIKJ*KXY
              endif
             enddo 
            endif
           enddo
         enddo
        enddo

!---1st variable - EDXY(I,J,1) has real units of EBOX (kg or m2)
!---2nd-6th variables are the XY-moments and dimensionless.
        do J = 1,JM
         do I = 1,IM
           if (EDXY(I,J,1) .gt. 1.d-30) then
             EDXY(I,J,2) = EDXY(I,J,2)/EDXY(I,J,1)
             EDXY(I,J,3) = EDXY(I,J,3)/EDXY(I,J,1)
             EDXY(I,J,4) = EDXY(I,J,4)/EDXY(I,J,1)
             EDXY(I,J,5) = EDXY(I,J,5)/EDXY(I,J,1)
             EDXY(I,J,6) = EDXY(I,J,6)/EDXY(I,J,1)
           endif
         enddo
        enddo
      else
        do J = 1,JM
        do I = 1,IM
          EDXY(I,J,1) = 0.d0
          do KJ = 1,JG
            if (B_GYA(KJ,J) .gt. 1.d-5) then
              do KI = 1,IG
                if (B_GXA(KI,I) .gt. 1.d-5) then
                  EMKIKJ = EBOX(KI,KJ) * B_GXA(KI,I) * B_GYA(KJ,J)
                  EDXY(I,J,1) = EDXY(I,J,1) + EMKIKJ
                endif
              enddo 
            endif
          enddo
        enddo
        enddo
      endif

!---sum check for complete overlaps (done in calling program)
!       SUMIJ = 0.d0       !!!need to declare sums
!      do J = 1,JM
!       do I = 1,IM
!        SUMIJ = SUMIJ + EDXY(I,J,1)
!       enddo
!      enddo
!        SUMK  = 0.d0
!      do KJ = 1,JG
!       do KI = 1,IG
!        SUMK = SUMK + EBOX(KI,KJ)
!       enddo 
!      enddo
!      if ( abs(SUMIJ-SUMK) .gt. 1.d-6*SUMK) then
!         write (6,'(a,1p,e12.5,a,e12.5)')
!     &     ' emissions mis-mapping: ',SUMK,' to ',SUMIJ
!      endif

      return
      end subroutine

      PURE SUBROUTINE Atmosphere(alt, sigma, delta, theta)
       !   -------------------------------------------------------------------------
       ! PURPOSE - Compute the properties of the 1976 standard atmosphere to 86 km.
       ! AUTHOR - Ralph Carmichael, Public Domain Aeronautical Software
       ! NOTE - If alt > 86, the values returned will not be correct, but they will
       !   not be too far removed from the correct values for density.
       !   The reference document does not use the terms pressure and temperature
       !   above 86 km.
       IMPLICIT NONE
       !============================================================================
       !     A R G U M E N T S                                                     |
       !============================================================================
       REAL*8,INTENT(IN)::  alt        ! geometric altitude, km.
       REAL*8,INTENT(OUT):: sigma      ! density/sea-level standard density
       REAL*8,INTENT(OUT):: delta      ! pressure/sea-level standard pressure
       REAL*8,INTENT(OUT):: theta      ! temperature/sea-level standard temperature
       !============================================================================
       !     L O C A L   C O N S T A N T S                                         |
       !============================================================================
       REAL*8,PARAMETER:: REARTH = 6369.0                 ! radius of the Earth (km)
       REAL*8,PARAMETER:: GMR = 34.163195                     ! hydrostatic constant
       INTEGER,PARAMETER:: NTAB=8       ! number of entries in the defining tables
       !============================================================================
       !     L O C A L   V A R I A B L E S                                         |
       !============================================================================
       INTEGER:: i,j,k                                                  ! counters
       REAL:: h                                       ! geopotential altitude (km)
       REAL:: tgrad, tbase      ! temperature gradient and base temp of this layer
       REAL:: tlocal                                           ! local temperature
       REAL:: deltah                             ! height above base of this layer
       !============================================================================
       !     L O C A L   A R R A Y S   ( 1 9 7 6   S T D.  A T M O S P H E R E )   |
       !============================================================================
       REAL*8,DIMENSION(NTAB),PARAMETER:: htab= 
     &       (/0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852/)
       REAL*8,DIMENSION(NTAB),PARAMETER:: ttab= 
     &      (/288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 
     &      214.65, 186.946/)
       REAL*8,DIMENSION(NTAB),PARAMETER:: ptab= 
     &       (/1.0, 2.233611E-1, 5.403295E-2, 8.5666784E-3, 
     &        1.0945601E-3, 6.6063531E-4, 3.9046834E-5, 3.68501E-6/)
       REAL*8,DIMENSION(NTAB),PARAMETER:: gtab= 
     &       (/-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0/)
       !----------------------------------------------------------------------------
       h=alt*REARTH/(alt+REARTH)      ! convert geometric to geopotential altitude
       
       i=1
       j=NTAB                                       ! setting up for binary search
       DO
          k=(i+j)/2                                              ! integer division
          IF (h < htab(k)) THEN
             j=k
          ELSE
             i=k
          END IF
          IF (j <= i+1) EXIT
       END DO
       
       tgrad=gtab(i)                                     ! i will be in 1...NTAB-1
       tbase=ttab(i)
       deltah=h-htab(i)
       tlocal=tbase+tgrad*deltah
       theta=tlocal/ttab(1)                                    ! temperature ratio
       
       IF (tgrad == 0.0) THEN                                     ! pressure ratio
          delta=ptab(i)*EXP(-GMR*deltah/tbase)
       ELSE
          delta=ptab(i)*(tbase/tlocal)**(GMR/tgrad)
       END IF
       
       sigma=delta/theta        ! density ratio

       RETURN
      END Subroutine Atmosphere 


      !--------------------------------------------------------------------


      ELEMENTAL SUBROUTINE ALT2PRESSURE( ALT, PRESS ) 
      
         !==========================================
         ! Wrapper routine for Atmosphere. 
         ! Calculates pressure (hPa) for a given altitude (km)
         ! using the 1976 US Standard Atmosphere
         ! (cdh, 01/26/2011)
         !==========================================

         !Arguments
         REAL*8, INTENT(IN)  :: ALT !km
         REAL*8, INTENT(OUT) :: PRESS !hPa
        
         ! Local variables
         REAL*8, PARAMETER   :: P0 = 1013.25 !hPa 
         REAL*8              :: S, P, T

         !===================================

         CALL ATMOSPHERE( ALT, S, P, T )
         PRESS = P * P0

      END SUBROUTINE
!-----------------------------------------------------------------------
!+
! NAME:
!        REGRID_COLUMN_WEIGHTS, REGRID_COLUMN_FAST
!
! PURPOSE:
!        Vertically regrids a column quantity in such a way
!        as to preserve the total mass. We may regrid many columns
!        by calculating the weights only once, then applying them repeatedly
!
! CATEGORY:
!        Regridding
!
! CALLING SEQUENCE:
!        CALL REGRID_COLUMN_WEIGHTS( PEDGE1, PEDGE2, WEIGHTS )
!        CALL REGRID_COLUMN_FAST( DATA1, DATA2, WEIGHTS )
!
! INPUTS:
!        DATA1 -> Column vector containing data on the original
!             grid.  DATA1 must be a mass-like quantity that does
!             not have any vertical dependence (e.g. molecules,
!             g, kg, kg/m2, molec/cm2, etc.) 
!
!        PEDGE1 -> Pressures [hPa] at the edges of each layer on 
!             the old vertical grid.  PEDGE1 will thus have one 
!             more element element than DATA1 (since DATA1 is 
!             specified on the midpoints of each layer).
!
!        PEDGE2 -> Pressures [hPa] at the edges of each layer on 
!             the new vertical grid.  PEDGE2 will thus have one 
!             more element element than DATA2 (since DATA2 is 
!             specified on the midpoints of each layer).
!
! OUTPUTS:
!        DATA2 -> Column vector containing data on the new vertical 
!             grid.  The column sum of DATA2 will equal that of
!             DATA1. DATA2 will also be a mass-like quantity.
!
!        WEIGHTS -> 2D array (I,J) with regridding weights, consisting of the
!             fraction of each input grid box I contained within the 
!             output grid box J
!
!
! SUBROUTINES:
!        None
!
! REQUIREMENTS:
!        None
!
! NOTES:
!        Converted from IDL code (GAMAP2 routine regrid_column.pro) by
!        C.D. Holmes (1/28/2011)
!
! EXAMPLE:
!-----------------------------------------------------------------------

      SUBROUTINE Regrid_Column_Weights( PEdge1, PEdge2, Weights )

      IMPLICIT NONE

      ! Arguments
      REAL*8,DIMENSION(:),INTENT(IN) :: PEDGE1, PEDGE2
      REAL*8,DIMENSION(:,:),INTENT(OUT):: Weights
      
      ! Local variables
      INTEGER           :: LM1, LM2, L, K
      LOGICAL           :: FIRST, VALID
      LOGICAL,PARAMETER :: CHECK=.false., DEBUG=.false.

      !========================================
      ! Begin subroutine
      !========================================
      
      !====================================================================
      ! Initialize variables
      !====================================================================
      LM1      = size( PEdge1 ) - 1
      LM2      = size( PEdge2 ) - 1
      Weights = 0d0
      First    = .true.
      Valid    = .false.      

      ! Check if array Weights has correct dimensions LM1xLM2

      !### Debug output
      if ( Debug ) then 
         print*, '### LM1, LM2: ', LM1,  LM2
         print*, '### Input Pressure Edges'
         print*, Pedge1
      endif

      !====================================================================
      ! Determine fraction of each INPUT box 
      ! which contributes to each OUTPUT box
      !====================================================================
 
      ! Loop over INPUT layers
      do L=1, LM1
       
         ! Reset VALID flag
         Valid = .false.
  

         ! If the thickness of this pressure level is zero, then this 
         ! means that this pressure level lies below the surface 
         ! pressure (due to topography), as set up in the calling
         ! program.  Therefore, skip to the next INPUT level.
         ! This also helps avoid divide by zero errors. (bmy, 8/6/01)
         IF ( ( PEdge1(L) - Pedge1(L+1) ) < 1e-5 ) CYCLE
 
         ! Loop over OUTPUT layers
         do K=1, LM2 
                    
            !==============================================================
            ! No contribution if:
            ! -------------------
            ! Bottom of OUTPUT layer above Top    of INPUT layer  OR
            ! Top    of OUTPUT layer below Bottom of INPUT layer
            !==============================================================
            if ( ( PEdge2(K)   < PEdge1(L+1))   .OR. 
     &           ( PEdge2(K+1) > PEdge1(L)  )   ) CYCLE
 
            !==============================================================
            ! Contribution if: 
            ! ----------------
            ! Entire INPUT layer in OUTPUT layer
            !==============================================================
            if ( (PEdge2(K)   >= PEdge1(L))     .AND. 
     &           (PEdge2(K+1) <= PEdge1(L+1)) ) then 
               
               Weights(L,K) = 1d0

               !### Debug output
               IF (DEBUG)
     &              Print'(2i12,6x,f13.6,i2)', 
     &              L, K, Weights(L,K), 1
 
               ! Indicate a valid contribution from L to K
               Valid = .true.
 
               ! Go to next K iteration
               CYCLE
            endif
           
            !==============================================================
            ! Contribution if: 
            ! ----------------
            ! Top of OUTPUT layer in INPUT layer
            !==============================================================
            if ( (PEdge2(K+1) <= PEdge1(L))   .AND. 
     &           (PEdge2(K)   >= PEdge1(L)) ) THEN 
 
               Weights(L,K) = ( PEdge1(L) - PEdge2(K+1) ) / 
     &              ( PEdge1(L) - PEdge1(L+1) ) 
 
               !### Debug output
               IF (DEBUG)
     &              Print'(2i12,6x,f13.6,i2)', 
     &              L, K, Weights(L,K), 2
 
               ! Indicate a valid contribution from L to K
               Valid = .true.
 
               ! Go to next K iteration
               CYCLE
            endif

            !==============================================================
            ! Contribution if: 
            ! ----------------
            ! Entire OUTPUT layer in INPUT layer
            !==============================================================
            if ( (PEdge2(K)   <= PEdge1(L))     .AND. 
     &           (PEdge2(K+1) >= PEdge1(L+1)) ) then
 
               Weights(L,K) = ( PEdge2(K) - PEdge2(K+1) ) / 
     &              ( PEdge1(L) - PEdge1(L+1) )
 
               ! Also add the to the first OUTPUT layer the fraction
               ! of the first INPUT layer that is below sigma = 1.0
               ! This is a condition that can be found in GEOS-3 data.
               if ( ( First                 )   .AND.  
     &              ( K == 1                )   .AND.   
     &              ( PEdge1(L) > PEdge2(1) ) ) then
 
                  Weights(L,K) = Weights(L,K) +              
     &                 ( PEdge1(L) - PEdge2(1)   ) / 
     &                 ( PEdge1(L) - PEdge1(L+1) )                
 
                  ! We only need to do this once...
                  First = .false.
               endif
 
               !### Debug output
               IF (DEBUG)
     &              Print'(2i12,6x,f13.6,i2)', 
     &              L, K, Weights(L,K), 3
 
               ! Indicate a valid contribution from L to K
               Valid = .true.
 
               ! Go to next K iteration
               CYCLE
            endif
            
            !==============================================================
            ! Contribution if: 
            ! ----------------
            ! Bottom of OUTPUT layer in INPUT layer
            !==============================================================
            if ( (PEdge2(K)   >= PEdge1(L+1))   .AND.  
     &           (PEdge2(K+1) <= PEdge1(L+1)) ) then
            
               Weights(L,K) = ( PEdge2(K) - PEdge1(L+1) ) / 
     &              ( PEdge1(L) - PEdge1(L+1) )
            
               ! Also add the to the first OUTPUT layer the fraction
               ! of the first INPUT layer that is below sigma = 1.0
               ! This is a condition that can be found in GEOS-3 data.
               if ( ( First                 )   .AND. 
     &              ( K == 0                )   .AND.  
     &                 ( PEdge1(L) > PEdge2(1) ) ) then 

                  Weights(L,K) = Weights(L,K) +               
     &                 ( PEdge1(L) - PEdge2(1)   ) / 
     &                 ( PEdge1(L) - PEdge1(L+1) )                
                  
                  ! We only need to do this once...
                  First = .false.
               endif
 
               !### Debug output
               IF (DEBUG)
     &              Print'(2i12,6x,f13.6,i2)', 
     &              L, K, Weights(L,K), 4
 
               ! Indicate a valid contribution from L to K
               Valid = .true.
 
               ! Go to next K iteration
               CYCLE
            endif
            

         enddo

         !=================================================================
         ! Consistency Check:
         ! ------------------
         ! If SUM( WEIGHTS(L,:) ) does not = 1, there is a problem.
         ! Test those INPUT layers (L) which make a contribution to 
         ! OUTPUT layers (K) for this criterion.
         !
         ! NOTE: This will be skipped if /NO_CHECK is set (bmy, 3/14/02)
         !=================================================================
         if ( Valid .AND. Check ) then 
            if ( Abs( 1e0 - sum( Weights(L,:) ) ) >= 1e-4 ) THEN
               print*, 'Weights does not add to 1;;'
               print '(''L, SUM( WEIGHTS(L,:) ): '', i4, 1x, f13.7 )', 
     &              L, sum( Weights(L,:) ) 
               stop
            endif
         endif


      enddo


      END SUBROUTINE Regrid_Column_Weights

      !------------------------------------------------------------

      SUBROUTINE Regrid_Column_Fast( Data1, Data2, Weights )
      
      IMPLICIT NONE

      ! Arguments
      REAL*8,DIMENSION(:)  ,INTENT(IN) :: DATA1
      REAL*8,DIMENSION(:,:),INTENT(IN) :: Weights
      REAL*8,DIMENSION(:),  INTENT(OUT):: DATA2
      
      ! Local variables
      INTEGER           :: LM1, LM2, L, K
      LOGICAL,PARAMETER :: CHECK=.false., DEBUG=.false.

      !========================================
      ! Begin subroutine
      !========================================
      
      !====================================================================
      ! Initialize variables
      !====================================================================
      LM1      = size( Weights, DIM=1 ) 
      LM2      = size( Weights, DIM=2 ) 

      !====================================================================
      ! Compute "new" data -- multiply "old" data by fraction of
      ! "old" data residing in the "new" layer
      !====================================================================
      do K=1, LM2
         do L=1, LM1
            Data2(K) = Data2(K) + ( Data1(L) * Weights(L,K) )
         enddo
      enddo
 
      !### Debug output
      if ( Debug ) then 
!         print*, '### Output Data:'
!         print*, Data2
!         print*, '### Output Data: column sum'
!         print*, Sum( Data2 )
!         print*, '### Output Data: ratio sum(new)/sum(old):'
!         print*, Sum( Data2 ) / Sum( Data1 )
      endif
   
      END SUBROUTINE Regrid_Column_Fast

      !------------------------------------------------------------


      END MODULE RCP_MOD
