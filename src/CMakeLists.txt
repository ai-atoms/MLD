# order matter
if("$ENV{MLD_NDM}" STREQUAL "ON")
  add_subdirectory(ndm)
endif()

add_subdirectory(extended)
add_subdirectory(minextended)
add_subdirectory(mld)
add_subdirectory(scatools)

if("$ENV{MLD_NDM}" STREQUAL "ON")
  add_executable(milady_main.exe ndm/main.F90)
else()
  add_executable(milady_main.exe mld/main.F90)
endif()
# no need add_dependencies(milady_main.exe MILADY)

message( STATUS "*** MLD_MKL_LIB : ${MLD_MKL_LIB}")
message( STATUS "*** MLD_HDF_LIB : ${MLD_HDF_LIB}")

#target_link_libraries( milady_main.exe PRIVATE  MILADY ${MLD_MKL_LIB} )
#target_link_libraries( milady_main.exe PRIVATE MPI::MPI_Fortran MILADY ${MLD_HDF_LIB} )
if ("$ENV{MLD_HDF5}" STREQUAL "ON")
  target_link_libraries( milady_main.exe PRIVATE  MILADY ${MLD_MKL_LIB} ${MLD_HDF_LIB} )
else()
  target_link_libraries( milady_main.exe PRIVATE  MILADY ${MLD_MKL_LIB} )
endif()

target_include_directories(milady_main.exe PRIVATE ${MLD_BUILD_DIR}/include)
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set_target_properties(milady_main.exe PROPERTIES LINKER_LANGUAGE Fortran)

