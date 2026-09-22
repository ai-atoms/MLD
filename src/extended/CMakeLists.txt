enable_language( Fortran )

if("$ENV{MLD_HDF5}" STREQUAL "ON")

if(NOT DEFINED HDF5_DIR)
  include(ExternalProject)

  message( STATUS "MPI_FC : " ${MPI_Fortran_COMPILER} )
  message( STATUS "MPI_CC : " ${MPI_C_COMPILER} )

  if ((Fortran_COMPILER_NAME STREQUAL "ifx") OR (Fortran_COMPILER_NAME STREQUAL "ifort"))
    set(OMP_ROOT_LIB ${OMP_ROOT}/lib64)
  elseif (Fortran_COMPILER_NAME STREQUAL "gfortran")
    set(OMP_ROOT_LIB ${OMP_ROOT}/lib)
  else()
    message( FATAL_ERROR "*** HDF5: No optimized Fortran compiler flags are known, we just try NOTHING!" )
  endif()

  message( STATUS "OMP_ROOT_LIB : " ${OMP_ROOT_LIB} )
  message( STATUS "OMP_ROOT : " ${OMP_ROOT} )

  set(hdf5_build_dir "$ENV{MLD_BUIDIR}/hdf5")
  list(
    APPEND
    hdf5_options
    "-DHDF5_BUILD_CPP_LIB=OFF"  
    "-DHDF5_ENABLE_USING_MEMCHECKER=ON"
    "-DHDF5_ENABLE_DIRECT_VFD=OFF"
    "-DHDF5_ENABLE_Z_LIB_SUPPORT=OFF"
    "-DHDF5_ENABLE_SZIP_SUPPORT=OFF"
    "-DHDF5_ENABLE_SZIP_ENCODING=OFF"
    "-DHDF5_ENABLE_PARALLEL=ON"
    "-DHDF5_BUILD_FORTRAN=ON")
  ExternalProject_Add(
    hdf5_lib
    PREFIX ${hdf5_build_dir}
    URL ${CMAKE_CURRENT_SOURCE_DIR}/hdf5-1.14.3.tar.gz
    CMAKE_ARGS  "${hdf5_options}"
    "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
    "-DCMAKE_INSTALL_PREFIX=${MLD_BUILD_DIR}"
    "-DCMAKE_BUILD_TYPE=Release"
    "-DMPI_C_COMPILER=${MPI_C_COMPILER}"
    "-DMPI_Fortran_COMPILER=${MPI_Fortran_COMPILER}"
    # "-DCMAKE_C_FLAGS=-I${OMP_ROOT}/include -L${OMP_ROOT}/lib64 -lmpi"
    # "-DCMAKE_Fortran_FLAGS=-I${OMP_ROOT}/include -L${OMP_ROOT}/lib64 -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi"
    "-DCMAKE_C_FLAGS=-I${OMP_ROOT}/include -L${OMP_ROOT_LIB} -lmpi"
    "-DCMAKE_Fortran_FLAGS=-I${OMP_ROOT}/include -L${OMP_ROOT_LIB} -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi"
    LOG_DOWNLOAD ON
    LOG_UPDATE ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
    LOG_INSTALL ON)

  add_dependencies(MILADY hdf5_lib)
  target_include_directories(MILADY PRIVATE ${MLD_BUILD_DIR}/include)
  set(MLD_HDF_LIB ${MLD_BUILD_DIR}/lib/libhdf5_fortran.so  ${MLD_BUILD_DIR}/lib/libhdf5_hl_fortran.so PARENT_SCOPE)
endif()

endif()

# order no matter.
target_sources( MILADY  PRIVATE
  soap_turbo.f90
  soap_turbo_angular.f90
  soap_turbo_compress.f90
  soap_turbo_functions.f90
  soap_turbo_radial.f90
  npy.f90
  print_matrix.f90
  nr.f90
)

