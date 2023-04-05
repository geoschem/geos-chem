# NcdfUtil: NetCDF Utility routines for GEOS-Chem

This folder contains netCDF utiliity routines for GEOS-Chem.

## Contents

- `CMakeLists.txt`: CMake build file
- `charpak_mod.F90`: Copy of `Headers/charpak_mod.F90`, used locally.
- `julday_mod.F90`: Copy of `Headers/julday_mod.F90`, used locally.
- `m_do_err_out.F90`: Error handling module
- `m_netcdf_io_checks.F90`: Error checking routines
- `m_netcdf_io_close.F90`: Routines to close netCDF files
- `m_netcdf_io_create.F90`: Routines to create netCDF files
- `m_netcdf_io_define.F90`: Routines to define netCDF variables
- `m_netcdf_io_get_dimlen.F90`: Reoutines
- `m_netcdf_io_handle_err.F90`: Error checking routines
- `m_netcdf_io_open.F90`: Routines for opening netCDF files
- `m_netcdf_io_readattr.F90`: Routines for reading netCDF attributes
- `m_netcdf_io_read.F90`: Routines for reading data to a netCDF file
- `m_netcdf_io_write.F90`: Routines for writing data to a netCDF file
- `ncdf_mod.F90`: Convenience routines for netCDF handling
- `TestNcdfUtil.F90`: Test program 

## Scripts

We have now moved netCDF utility scripts (such as `isCoards` and `nc_chunk.pl`) to a separate Github repository.  You may download them from https://github.com/geoschem/netcdf-scripts.



