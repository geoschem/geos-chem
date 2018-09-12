add_static_target(NcdfUtil
    SOURCES
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
    
    PUBLIC_LINK
        GeosUtil
        ${NETCDF_F90_LIBRARIES}

    PUBLIC_INCLUDE
        ${NETCDF_F90_INCLUDE_DIR}
)