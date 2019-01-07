add_library(NcdfUtil STATIC
	m_do_err_out.F90
	m_netcdf_io_checks.F90
	m_netcdf_io_close.F90
	m_netcdf_io_create.F90
	m_netcdf_io_define.F90
	m_netcdf_io_get_dimlen.F90
	m_netcdf_io_handle_err.F90
	m_netcdf_io_open.F90
	m_netcdf_io_read.F90
	m_netcdf_io_readattr.F90
	m_netcdf_io_write.F90
	ncdf_mod.F90
)
target_link_libraries(NcdfUtil 
	PUBLIC GeosUtil ${NETCDF_F90_LIBRARIES}
)
target_include_directories(NcdfUtil
	PUBLIC ${NETCDF_F90_INCLUDE_DIR}
)
