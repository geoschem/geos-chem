program try_compile
    use netcdf

    implicit none

    character (len = *), parameter :: FILE_NAME = "simple_xy.nc"
    integer, parameter :: NDIMS = 2
    integer, parameter :: NX = 6, NY = 12
    integer :: ncid, varid, dimids(NDIMS)
    integer :: x_dimid, y_dimid
    integer :: data_out(NY, NX)
    integer :: x, y
#ifndef NO_OMP
    integer nthreads, tid, OMP_GET_THREAD_NUM

    write(*,*) 'About to start multiple threads'
    ! Try OpenMP
    !$OMP PARALLEL PRIVATE(nthreads, tid)
    TID = OMP_GET_THREAD_NUM()
    write(*,*) 'Hello from thread ', tid
    !$OMP END PARALLEL
#endif

    ! Try NetCDF-F
    do x = 1, NX
       do y = 1, NY
          data_out(y, x) = (x - 1) * NY + (y - 1)
       end do
    end do
    write(*,*) 'Creating NetCDF file'
    call check( nf90_create(FILE_NAME, NF90_CLOBBER, ncid) )
    call check( nf90_def_dim(ncid, "x", NX, x_dimid) )
    call check( nf90_def_dim(ncid, "y", NY, y_dimid) )

    ! The dimids array is used to pass the IDs of the dimensions of
    ! the variables. Note that in fortran arrays are stored in
    ! column-major format.
    dimids =  (/ y_dimid, x_dimid /)

    call check( nf90_def_var(ncid, "data", NF90_INT, dimids, varid) )

    call check( nf90_enddef(ncid) )
    call check( nf90_put_var(ncid, varid, data_out) )
    call check( nf90_close(ncid) )
    write(*,*) 'Finished creating the NetCDF file'

  contains
    subroutine check(status)
      integer, intent ( in) :: status

      if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop "Stopped"
      end if
    end subroutine check

end program
