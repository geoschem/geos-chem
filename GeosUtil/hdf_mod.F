!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hdf_mod
!
! !DESCRIPTION: Module HDF\_MOD contains routines to write data to 
!  HDF5 files.
!\\
!\\
! !INTERFACE: 
!
      MODULE HDF_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CLEANUP_HDF
      PUBLIC :: CLOSE_HDF 
      PUBLIC :: INIT_HDF
      PUBLIC :: OPEN_HDF
      PUBLIC :: WRITE_HDF 
!
! !PUBLIC DATA MEMBERS:
! 
      PUBLIC :: HDFCATEGORY
      PUBLIC :: HDFDESCRIPT
      PUBLIC :: HDFNAME
      PUBLIC :: HDFFNAME
      PUBLIC :: HDFMOLC
      PUBLIC :: HDFMWT
      PUBLIC :: HDFSCALE
      PUBLIC :: HDFUNIT

      CHARACTER(LEN=40), ALLOCATABLE :: HDFCATEGORY(:)
      CHARACTER(LEN=40), ALLOCATABLE :: HDFDESCRIPT(:)
      INTEGER                        :: MAXDIAG
      INTEGER                        :: MAXTRACER
      INTEGER                        :: MAXCAT
      INTEGER,           ALLOCATABLE :: HDFMOLC(:,:)
      REAL*4,            ALLOCATABLE :: HDFMWT(:,:)
      REAL*4,            ALLOCATABLE :: HDFSCALE(:,:)
      CHARACTER(LEN=40), ALLOCATABLE :: HDFNAME(:,:)
      CHARACTER(LEN=40), ALLOCATABLE :: HDFFNAME(:,:)
      CHARACTER(LEN=40), ALLOCATABLE :: HDFUNIT(:,:)
!
! !REMARKS:
!  If you have the HDF5 library installed on your system, then you can 
!  compile GEOS-Chem with the option:
!                                                                             .
!      make HDF5=yes
!                                                                             .
!  which will activate the HDF5-specific code in this module to enable file
!  I/O to HDF5 format.  You must also specify the HDF5 include and library 
!  paths in the Makefile_header.mk.  
!                                                                             .
!  The default is not to activate the HDF5-specific code.
!
! !REVISION HISTORY:
!  19 Nov 2009 - A. van Donkelaar - Initial Version
!  21 Dec 2009 - R. Yantosca      - Modified to block out HDF5-specific code
!                                   so that users w/o HDF5 library can still
!                                   compile & run GEOS-Chem
!  21 Dec 2009 - R. Yantosca      - Updated comments
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
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
! !IROUTINE: open_hdf
!
! !DESCRIPTION: Subroutine OPEN\_HDF creates and opens an hdf file for output.
!\\
!\\
! !INTERFACE:
!
#if   defined( USE_HDF5 )
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%    NOTE: Subroutine is used only when USE_HDF5 is defined!    %%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE OPEN_HDF( IU_HDF, FILENAME, IMAX, IMIN, 
     &                     JMAX,   JMIN,     NI,   NJ    )
!
! !USES:
!
      USE HDF5
      USE GRID_MOD, ONLY : GET_XMID
      USE GRID_MOD, ONLY : GET_YMID

      USE CMN_SIZE_MOD    ! Size parameters
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN)    :: FILENAME   ! File name to open
      INTEGER,          INTENT(IN)    :: IMIN       ! Min longitude index
      INTEGER,          INTENT(IN)    :: IMAX       ! Max longitude index
      INTEGER,          INTENT(IN)    :: JMIN       ! Min latitude index
      INTEGER,          INTENT(IN)    :: JMAX       ! Max latitude index
      INTEGER,          INTENT(IN)    :: NI         ! # of longitudes
      INTEGER,          INTENT(IN)    :: NJ         ! # of latitudes
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER(HID_T),   INTENT(INOUT) :: IU_HDF     ! HDF5 file identifier
!
! !REVISION HISTORY: 
!  Nov 20 2009 - A. van Donkelaar - Initial Version
!  21 Dec 2009 - R. Yantosca      - Modified to block out HDF5-specific code
!                                   so that users w/o HDF5 library can still
!                                   compile & run GEOS-Chem
!  21 Dec 2009 - R. Yantosca      - Updated comments
!  01 Mar 2012 - R. Yantosca      - Now use GET_XMID(I,J,L), GET_YMID(I,J,L)
!                                   from the new grid_mod.F90
!  17 Oct 2012 - M. Payer         - Bug fix: Reorder indices for GET_YMID,
!                                   otherwise LAT variable in HDF files will
!                                   be wrong (S. Philip)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER(HID_T)                 :: dset_id    ! Dataset identifier
      INTEGER(HID_T)                 :: dspace_id  ! Dataspace identifier
      INTEGER(HSIZE_T), DIMENSION(1) :: dims       ! Dataset dimensions
      INTEGER                        :: rank = 1   ! Dataset rank
      CHARACTER(LEN=60)              :: dsetname   ! Dataset name
      INTEGER                        :: error      ! Error flag
      REAL*4                         :: LAT(NJ)
      REAL*4                         :: LON(NI)
      INTEGER                        :: I, J

      ! Initialize FORTRAN interface.
      CALL h5open_f( error )

      ! Create a new file using default properties.
      CALL h5fcreate_f( TRIM(FILENAME), H5F_ACC_TRUNC_F, IU_HDF, error )

      ! Write LAT/LON to new file
      ! Does not wrap around date line
      DO I = 0, NI-1
         LON(I+1) = GET_XMID(I+IMIN,1,1)
      ENDDO

      DO J = 0, NJ-1
         LAT(J+1) = GET_YMID(1,J+JMIN,1)
      ENDDO

      DO I = 1,2

         IF (I .eq. 1) THEN
            dims = NI
            dsetname = 'LON'
         ELSEIF (I .eq. 2) THEN
            dims = NJ
            dsetname = 'LAT'
         ENDIF

         ! Create the dataspace.
         CALL h5screate_simple_f(rank, dims, dspace_id, error)

         ! Create the dataset with default properties.
         CALL h5dcreate_f(IU_HDF, dsetname, H5T_NATIVE_REAL, 
     &                 dspace_id, dset_id, error)

         ! Write the dataset.
         IF ( I .eq. 1 ) THEN
            CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, LON, dims, error)
         ELSEIF ( I .eq. 2 ) THEN
            CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, LAT, dims, error)
         ENDIF

         ! End access to the dataset and release resources used by it.
         CALL h5dclose_f(dset_id, error)

         ! Terminate access to the data space.
         CALL h5sclose_f(dspace_id, error)

      ENDDO

      END SUBROUTINE OPEN_HDF

#else

      SUBROUTINE OPEN_HDF
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%    NOTE: OPEN_HDF will be a stub if USE_HDF5 is undefined!    %%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END SUBROUTINE OPEN_HDF

#endif
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: close_hdf
!
! !DESCRIPTION: Subroutine CLOSE\_HDF closes an HDF5 file that is already open.
!\\
!\\
! !INTERFACE:
!
#if   defined( USE_HDF5 )
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%    NOTE: Subroutine is used only when USE_HDF5 is defined!    %%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE CLOSE_HDF( IU_HDF )
!
! !USES:
!
      USE HDF5
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER(HID_T), INTENT(INOUT) :: IU_HDF  ! HDF5 File identifier
!
! !REVISION HISTORY:
!  20 Nov 2009 - A. van Donkelaar - Initial Version
!  21 Dec 2009 - R. Yantosca      - Modified to block out HDF5-specific code
!                                   so that users w/o HDF5 library can still
!                                   compile & run GEOS-Chem
!  21 Dec 2009 - R. Yantosca      - Updated comments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: error  ! Error flag

      ! Initialize FORTRAN interface.
      CALL h5open_f ( error )

      ! Terminate access to the file.
      CALL h5fclose_f( IU_HDF, error )

      END SUBROUTINE CLOSE_HDF

#else

      SUBROUTINE CLOSE_HDF
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%   NOTE: CLOSE_HDF will be a stub if USE_HDF5 is undefined!    %%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END SUBROUTINE CLOSE_HDF

#endif
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_hdf
!
! !DESCRIPTION: Subroutine WRITE\_HDF writes data to an open HDF5 file.
!\\
!\\
! !INTERFACE:
!
#if   defined( USE_HDF5 )
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%    NOTE: Subroutine is used only when USE_HDF5 is defined!    %%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE WRITE_HDF( IU_HDF,  N,      NDCATEGORY, NDTRACER, 
     &                      NDUNIT,  TAU0,   TAU1,       RESERVED,
     &                      NI,      NJ,     NL,         IFIRST,
     &                      JFIRST,  LFIRST, ARRAY )
!
! !USES:
!
      USE HDF5
!
! !INPUT PARAMETERS:
!
      INTEGER,           INTENT(IN) :: IU_HDF            ! HDF file unit #
      INTEGER,           INTENT(IN) :: N                 ! Actual tracer #
      CHARACTER(LEN=40), INTENT(IN) :: NDCATEGORY        ! Diagnostic category
      INTEGER,           INTENT(IN) :: NDTRACER          ! Tracer # for file
      CHARACTER(LEN=40), INTENT(IN) :: NDUNIT            ! Units of data
      REAL*8,            INTENT(IN) :: TAU0              ! TAU at start & end
      REAL*8,            INTENT(IN) :: TAU1              !  of diag interval
      CHARACTER(LEN=40), INTENT(IN) :: RESERVED          ! Descriptive string
      INTEGER,           INTENT(IN) :: NI                ! # of longitudes
      INTEGER,           INTENT(IN) :: NJ                ! # of latitudes
      INTEGER,           INTENT(IN) :: NL                ! # of levels
      INTEGER,           INTENT(IN) :: IFIRST            ! Index of 1st lon
      INTEGER,           INTENT(IN) :: JFIRST            ! Index of 1st lat
      INTEGER,           INTENT(IN) :: LFIRST            ! Index of 1st lev
      REAL*4,            INTENT(IN) :: ARRAY(NI,NJ,NL)   ! Data array
!
! !REVISION HISTORY:
!  20 Nov 2009 - A. van Donkelaar - Initial Version
!  21 Dec 2009 - R. Yantosca      - Modified to block out HDF5-specific code
!                                   so that users w/o HDF5 library can still
!                                   compile & run GEOS-Chem
!  21 Dec 2009 - R. Yantosca      - Updated comments
!  17 Oct 2012 - M. Payer         - Add PBLDEPTH to fix bad HDF category error
!                                   (S. Philip, A. van Donkelaar)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: error             ! Error flag
      INTEGER(HID_T)    :: dset_id           ! Dataset identifier
      INTEGER(HID_T)    :: dspace_id         ! Dataspace identifier
      INTEGER(HSIZE_T)  :: dims(3)           ! Dataset dimensions
      INTEGER           :: rank = 3          ! Dataset rank
      CHARACTER(LEN=60) :: dsetname          ! Dataset name
      INTEGER(HID_T)    :: attr_id           ! Attribute identifier
      INTEGER(HID_T)    :: aspace_id         ! Attribute Dataspace identifier
      INTEGER(HID_T)    :: atype_id          ! Attribute Dataspace identifier
      INTEGER(HSIZE_T)  :: adims(1) = (/1/)  ! Attribute dimension
      INTEGER           :: arank = 1         ! Attribute rank
      INTEGER(SIZE_T)   :: attrlen           ! Length of the attribute string
      CHARACTER(LEN=20) :: aname             ! Attibute Name
      CHARACTER(LEN=40) :: attr_data         ! Attribute Data
      REAL*8            :: attr_data_num     ! Attribute Data

      INTEGER           :: I, GAMAPI, GAMAPI2

      ! Initialize FORTRAN interface.
      CALL h5open_f (error)

      dims(1) = NI
      dims(2) = NJ
      dims(3) = NL

      ! Create the dataspace.
      CALL h5screate_simple_f(rank, dims, dspace_id, error)

      GAMAPI = 0
      DO I = 1, MAXDIAG
         IF ( HDFCATEGORY(I) .EQ. NDCATEGORY ) THEN
            GAMAPI = I
         ENDIF
      ENDDO


      IF ( GAMAPI .eq. 0 ) THEN
         PRINT*, 'BAD HDF CATEGORY:', NDCATEGORY
         CALL FLUSH(6)
      ENDIF

      IF ( NDCATEGORY .eq. 'IJ-AVG-$' ) THEN
         GAMAPI2 = 45
      ELSEIF ( NDCATEGORY .eq. 'BXHGHT-$' ) THEN
         GAMAPI2 = 68
      ELSEIF ( NDCATEGORY .eq. 'DAO-3D-$' ) THEN
         GAMAPI2 = 66
      ELSEIF ( NDCATEGORY .eq. 'PBLDEPTH' ) THEN
         GAMAPI2 = 32
         HDFNAME(1,GAMAPI2) = 'PBL-HGHT-M'
         HDFFNAME(1,GAMAPI2) = 'Planetary Boundary Layer Height [m]'
         HDFUNIT(1,GAMAPI2) = 'm'
         HDFMOLC(1,GAMAPI2) = 1
         HDFMWT(1,GAMAPI2) = 0
         HDFSCALE(1,GAMAPI2) = 1

         HDFNAME(2,GAMAPI2) = 'PBL-HGHT-L'
         HDFFNAME(2,GAMAPI2) = 'Planetary Boundary Layer Height [lvl]'
         HDFUNIT(2,GAMAPI2) = 'level'
         HDFMOLC(2,GAMAPI2) = 1
         HDFMWT(2,GAMAPI2) = 0
         HDFSCALE(2,GAMAPI2) = 1
      ELSEIF ( NDCATEGORY .eq. 'OD-MAP-$' ) THEN
         GAMAPI2 = 21
      ELSEIF ( NDCATEGORY .eq. 'PEDGE-$' ) THEN
         GAMAPI2 = 31
         HDFNAME(NDTRACER,GAMAPI2) = 'PEDGE'
         HDFFNAME(NDTRACER,GAMAPI2) = 'Box Edge Pressure'
         HDFUNIT(NDTRACER,GAMAPI2) = 'hPa'
         HDFMOLC(NDTRACER,GAMAPI2) = 1
         HDFMWT(NDTRACER,GAMAPI2) = 0
         HDFSCALE(NDTRACER,GAMAPI2) = 1
      ELSEIF ( NDCATEGORY .eq. 'TIME-SER' ) THEN
         GAMAPI2 = 48
      ELSE
         GAMAPI2 = 0
      ENDIF

      IF (GAMAPI2 .eq. 0) THEN
         PRINT*, 'BAD HDF CATEGORY 2:', NDCATEGORY
         CALL FLUSH(6)
      ENDIF

      WRITE(dsetname,*) TRIM(NDCATEGORY), '::',
     &                  TRIM(HDFNAME(NDTRACER,GAMAPI2))

      ! remove space at beginning of dataset name
      ! don't understand why this is needed (amv)
      dsetname=dsetname(2:60)

      ! Create the dataset with default properties.
      CALL h5dcreate_f(IU_HDF, dsetname, H5T_NATIVE_REAL, dspace_id,
     &                 dset_id, error)

      ! Write the dataset.
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, ARRAY(:,:,:) *  
     &                HDFSCALE(NDTRACER,GAMAPI2), dims, error)

      DO I = 1,7 !1,7
         IF ( I .eq. 1 ) THEN
            attrlen = 40
            adims = 1
            aname = 'CATEGORY DESCRIPTION'
            attr_data = HDFDESCRIPT(GAMAPI)
         ELSEIF ( I .eq. 2 ) THEN
            attrlen = 40
            adims = 1
            aname = 'TRACER DESCRIPTION'
            attr_data = HDFFNAME(NDTRACER,GAMAPI2)
         ELSEIF ( I .eq. 3 ) THEN
            attrlen = 40
            adims = 1
            aname = 'UNITS'
            attr_data = HDFUNIT(NDTRACER,GAMAPI2)
         ELSEIF ( I .eq. 4 ) THEN
            attrlen = 8
            adims = 1
            aname = 'SCALE'
            attr_data_num = HDFSCALE(NDTRACER,GAMAPI2)
         ELSEIF ( I .eq. 5 ) THEN
            attrlen = 8
            adims = 1
            aname = 'MWT'
            attr_data_num = HDFMWT(NDTRACER,GAMAPI2)
         ELSEIF ( I .eq. 6 ) THEN
            attrlen = 8
            adims = 1
            aname = 'MOLC'
            attr_data_num = HDFMOLC(NDTRACER,GAMAPI2)
         ELSEIF ( I .eq. 7 ) THEN
            attrlen = 8
            adims = 1
            aname = 'GMTRACER'
            attr_data_num = NDTRACER
         ENDIF

         ! Create scalar data space for the attribute.
         CALL h5screate_simple_f(arank, adims, aspace_id, error)

         ! Create datatype for the attribute.
         IF (I < 4) THEN
            CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
         ELSE
            CALL h5tcopy_f(H5T_NATIVE_DOUBLE, atype_id, error)
         ENDIF

         CALL h5tset_size_f(atype_id, attrlen, error)

         ! Create dataset attribute.
         CALL h5acreate_f(dset_id, aname, atype_id, aspace_id, 
     &                 attr_id, error)

         ! Write the attribute data.
         IF (I < 4) THEN
            CALL h5awrite_f(attr_id, atype_id, attr_data, 
     &                      adims, error)
         ELSE
            CALL h5awrite_f(attr_id, atype_id, attr_data_num, 
     &                      adims, error)
         ENDIF

         ! Close the attribute.
         CALL h5aclose_f(attr_id, error)

      ENDDO

      ! End access to the dataset and release resources used by it.
      CALL h5dclose_f(dset_id, error)

      ! Terminate access to the data space.
      CALL h5sclose_f(dspace_id, error)

      END SUBROUTINE WRITE_HDF

#else

      SUBROUTINE WRITE_HDF
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%   NOTE: WRITE_HDF will be a stub if USE_HDF5 is undefined!    %%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END SUBROUTINE WRITE_HDF

#endif
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_hdf
!
! !DESCRIPTION: Subroutine INIT\_HDF allocates all module variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_HDF( GMMAXCAT, GMMAXTRACER, GMMAXDIAG )
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: GMMAXCAT
      INTEGER, INTENT(IN) :: GMMAXTRACER
      INTEGER, INTENT(IN) :: GMMAXDIAG
!
! !REVISION HISTORY:
!  23 Nov 2009 - A. van Donkelaar - Initial Version
!  21 Dec 2009 - R. Yantosca      - Modified to block out HDF5-specific code
!                                   so that users w/o HDF5 library can still
!                                   compile & run GEOS-Chem
!  21 Dec 2009 - R. Yantosca      - Updated comments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER     ::   AS

      MAXCAT    = GMMAXCAT
      MAXTRACER = GMMAXTRACER
      MAXDIAG   = GMMAXDIAG

      ALLOCATE( HDFCATEGORY( MAXCAT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HDFCATEGORY' )
      HDFCATEGORY = ''

      ALLOCATE( HDFDESCRIPT( MAXCAT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HDFDESCRIPT' )
      HDFDESCRIPT = ''

      ALLOCATE( HDFNAME( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HDFNAME' )
      HDFNAME = ''

      ALLOCATE( HDFFNAME( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HDFFNAME' )
      HDFFNAME = ''

      ALLOCATE( HDFUNIT( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HDFUNIT' )
      HDFUNIT = ''

      ALLOCATE( HDFMOLC( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HDFMOLC' )
      HDFMOLC = 0

      ALLOCATE( HDFMWT( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HDFMWT' )
      HDFMWT = 0.0

      ALLOCATE( HDFSCALE( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HDFSCALE' )
      HDFSCALE = 0.0

      END SUBROUTINE INIT_HDF
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_hdf
!
! !DESCRIPTION: Subroutine CLEANUP\_HDF deallocates all module variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_HDF
!
! !REVISION HISTORY:
!  21 Dec 2009 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

      IF ( ALLOCATED( HDFCATEGORY ) ) DEALLOCATE( HDFCATEGORY )
      IF ( ALLOCATED( HDFDESCRIPT ) ) DEALLOCATE( HDFDESCRIPT )
      IF ( ALLOCATED( HDFNAME     ) ) DEALLOCATE( HDFNAME     )
      IF ( ALLOCATED( HDFFNAME    ) ) DEALLOCATE( HDFFNAME    )
      IF ( ALLOCATED( HDFUNIT     ) ) DEALLOCATE( HDFUNIT     )
      IF ( ALLOCATED( HDFMOLC     ) ) DEALLOCATE( HDFMOLC     )
      IF ( ALLOCATED( HDFMWT      ) ) DEALLOCATE( HDFMWT      )
      IF ( ALLOCATED( HDFSCALE    ) ) DEALLOCATE( HDFSCALE    )

      END SUBROUTINE CLEANUP_HDF
!EOC
      END MODULE HDF_MOD


