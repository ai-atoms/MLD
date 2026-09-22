################################################################################################
# thanks to
# https://stackoverflow.com/questions/52712416/cmake-fortran-module-directory-to-be-used-with-add-library
# https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/languages/fortran/ForFortranExample
# https://gitlab.com/CLIUtils/modern-cmake/-/blob/master/examples/extended-project/CMakeLists.txt
# https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html
#
# cmake3 -S<path-to-source> -B<path-to-build>
#
# all MLD_xxx vars are for EZ cmake resume debug with something as
# cmake3 -S${SRCDIR} -B${BUIDIR} -DMLD_OPT_TRACE=ON | grep ".. MLD_"
################################################################################################

CMAKE_minimum_required( VERSION 3.17 )

project( MILADY
         VERSION 2.0.0
         DESCRIPTION "Milady project with CMake"
         HOMEPAGE_URL "https://codev-tuleap.intra.cea.fr/plugins/git/milady"
         LANGUAGES Fortran C )

set( PROJECT_AUTHOR "Cosmin Marinica" )

# get MILADY cmake macro etc

# list( APPEND CMAKE_MODULE_PATH "$ENV{ONEAPI_ROOT}/mkl/2022.0.1/lib/cmake/mkl")
# list( APPEND CMAKE_MODULE_PATH "/volatile2/intel/oneapi/mkl/latest/lib/cmake/mkl")
list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_files" )

message( STATUS "CMAKE_MODULE_PATH : ${CMAKE_MODULE_PATH}" )

################################################
# utilities
################################################

include( MiladyNewMacros )
include( CTest )

################################################
# user options
# begin with "MLD_"
################################################

option( MLD_OPT_TRACE "Debug trace mode for cmake" OFF )
option( MLD_OPT_COMPILE_ONLY_CTESTS "Compile ONLY ctests Milady" OFF )
option( MLD_OPT_COMPILE_DOC "Compile doc Milady" OFF )

# set ENV{MPI_HOME} to get it
#set(MPI_C_COMPILER "$ENV{MPI_CC}" CACHE PATH "MPI C Compiler")
#set(MPI_Fortran_COMPILER "$ENV{MPI_FC}" CACHE PATH "MPI Fortran Compiler")
#set(CMAKE_PREFIX_PATH "/usr/local/iopenmpi" ${CMAKE_PREFIX_PATH})
#set(MPI_C_COMPILER $ENV{OMP_INSDIR}/bin/mpicc)
#set(MPI_Fortran_COMPILER $ENV{OMP_INSDIR}/bin/mpif90)
#set(MPI_CXX_COMPILER $ENV{OMP_INSDIR}/bin/mpicxx)
#
#set(MPI_C_INCLUDE_PATH $ENV{OMP_INSDIR}/include)
#set(MPI_C_LIBRARIES $ENV{OMP_INSDIR}/lib64/libmpi.so)
#set(MPI_Fortran_LIBRARIES $ENV{OMP_INSDIR}/lib64/libmpi_mpifh.so)
#set(MPI_CXX_LIBRARIES $ENV{OMP_INSDIR}/lib64/libmpi_cxx.so)

set(MPI_C_COMPILER "$ENV{MPI_CC}" CACHE PATH "MPI C Compiler")
if(DEFINED ENV{MPI_FC})
set(MPI_Fortran_COMPILER "$ENV{MPI_FC}" CACHE PATH "MPI Fortran Compiler")
endif() 
if(DEFINED ENV{FC})
  set(CMAKE_Fortran_COMPILER $ENV{FC})
endif()
set(MPI_CXX_COMPILER "$ENV{MPI_CXX}" CACHE PATH "MPI CXX Compiler")

# Hint for find_package(MPI) to prefer MPI_HOME
if(DEFINED ENV{MPI_HOME})
  list(PREPEND CMAKE_PREFIX_PATH "$ENV{MPI_HOME}")
  list(PREPEND CMAKE_LIBRARY_PATH "$ENV{MPI_HOME}/lib64")
  list(PREPEND CMAKE_LIBRARY_PATH "$ENV{MPI_HOME}/lib")
  # Explicitly set MPI library paths to prevent mixing Intel MPI and OpenMPI
  set(MPI_C_LIBRARIES "$ENV{MPI_HOME}/lib64/libmpi.so" CACHE PATH "MPI C Libraries")
  set(MPI_Fortran_LIBRARIES "$ENV{MPI_HOME}/lib64/libmpi_mpifh.so;$ENV{MPI_HOME}/lib64/libmpi_usempif08.so" CACHE PATH "MPI Fortran Libraries")
  set(MPI_CXX_LIBRARIES "$ENV{MPI_HOME}/lib64/libmpi.so" CACHE PATH "MPI CXX Libraries")
endif()

find_package(MPI REQUIRED)

# Configure MKL to use the correct MPI flavor based on which MPI we're using
# This MUST be set before find_package(MKL)
# Check the MPI compiler wrapper name from environment or detected paths
message(STATUS "*** DEBUG MPI_Fortran_COMPILER: ${MPI_Fortran_COMPILER}")
get_filename_component(MPI_Fortran_COMPILER_NAME "${MPI_Fortran_COMPILER}" NAME)
message(STATUS "*** MPI_Fortran_COMPILER_NAME: ${MPI_Fortran_COMPILER_NAME}")

if(MPI_Fortran_COMPILER_NAME MATCHES "mpif90" OR MPI_Fortran_COMPILER_NAME MATCHES "mpifort")
  # Using OpenMPI - set MKL_MPI before find_package
  set(MKL_MPI "openmpi")
  message(STATUS "*** Setting MKL_MPI to 'openmpi' for OpenMPI compatibility")
elseif(MPI_Fortran_COMPILER_NAME MATCHES "mpiifort" OR MPI_Fortran_COMPILER_NAME MATCHES "mpiifx")
  # Using Intel MPI
  set(MKL_MPI "intelmpi")
  message(STATUS "*** Setting MKL_MPI to 'intelmpi' for Intel MPI compatibility")
else()
  message(STATUS "*** WARNING: Could not determine MPI type from ${MPI_Fortran_COMPILER_NAME}")
endif()

# find_package(MKL REQUIRED) # because not found on irene...
find_package(MKL)
include_directories( SYSTEM ${MKL_INCLUDE} )

if( DEFINED HDF5_DIR )
  message( STATUS "*** COMPILE WITH USER DEFINED HDF5" )
  set( HDF5_USE_STATIC_LIBRARIES ON )
  # https://cmake.org/cmake/help/latest/module/FindHDF5.html
  # find_package(HDF5 COMPONENTS Fortran HL REQUIRED) # as Fortran not compiled in matix3
  find_package(HDF5 COMPONENTS Fortran REQUIRED)
  if (NOT HDF5_FOUND)
    message( FATAL_ERROR "HDF5 not found" )
  endif()
  include_directories( SYSTEM ${HDF5_INCLUDE_DIRS} )
  message( STATUS "HDF5_DIR : " ${HDF5_DIR} )
  message( STATUS "HDF5_INCLUDE_DIRS : " ${HDF5_INCLUDE_DIRS} )
  message( STATUS "HDF5_LIBRARIES : " ${HDF5_Fortran_LIBRARIES} )
endif()

# If you want yo use Intel MPI with Intel Fortran compiler, the wrapper is mpiifort.
# If you want to use an other MPI (such as Open MPI) with use mpi
# then you likely need to rebuild the MPI library because Fortran modules
# are not interoperable between vendors,
# and sometimes not even different compiler versions of the same vendor.
# I strongly discourage you to use ifort and manually set the include and library path.


if( MLD_OPT_COMPILE_DOC )
  find_package( Doxygen )  # TODO Sphinx instead
  if (Doxygen_FOUND)
    # TODO to make doc ? add_subdirectory( docs )
  else()
    message( STATUS "Doxygen not found, not building docs" )
  endif()
endif()

# before all cmake default installation
get_filename_component( default_prefix ".." ABSOLUTE )

# where make install do
set( MLD_DEFAULT_PREFIX ${default_prefix} )

################################################
# trace user options
# begin with "MLD_"
################################################

message( STATUS "*** cmake trace : " ${MLD_OPT_TRACE} )
message( STATUS "*** compile ONLY ctests milady : " ${MLD_OPT_COMPILE_ONLY_CTESTS} )

set( MLD_BUILD_DIR $ENV{MLD_BUIDIR} ) 
if (NOT MLD_BUILD_DIR)
  set( MLD_BUILD_DIR ${MLD_DEFAULT_PREFIX}/build CACHE STRING
       "Choose the build dir, default is : ${MLD_DEFAULT_PREFIX}/build."
       FORCE )
endif()

set( MLD_INSDIR $ENV{MLD_INSDIR} ) 
if (NOT MLD_INSDIR)
  set( MLD_INSDIR ${default_prefix}/mld_install CACHE STRING
       "Choose the build dir, default is: ${MLD_DEFAULT_PREFIX}/mld_install."
        FORCE )
endif()


# set( CMAKE_BUILD_TYPE Debug )
# set( CMAKE_BUILD_TYPE Release )

set( MKL_ROOT $ENV{MKL_ROOT} ) # TODO verify if not defined if(DEFINED ENV{<name>})
message( STATUS "*** CMAKE MKL_ROOT : " ${MKL_ROOT} )

set( OMP_ROOT $ENV{OMP_ROOT} ) # TODO verify if not defined if(DEFINED ENV{<name>})
message( STATUS "*** !!!!!! MPI_HOME : " $ENV{MPI_HOME} )
message( STATUS "*** !!!!!! MPI_DETERMINE_Fortran_CAPABILITIES : " ${MPI_DETERMINE_Fortran_CAPABILITIES} )
message( STATUS "*** !!!!!! CMAKE OMP_ROOT : " ${OMP_ROOT} )
message( STATUS "*** !!!!!! MPI_Fortran_MODULE_DIR : " ${MPI_Fortran_MODULE_DIR} )
message( STATUS "*** !!!!!! MPI_Fortran_HAVE_F08_MODULE : " ${MPI_Fortran_HAVE_F08_MODULE} )
message( STATUS "*** !!!!!! MPI_Fortran_HAVE_F90_MODULE: " ${MPI_Fortran_HAVE_F90_MODULE} )
message( STATUS "*** !!!!!! MPI_Fortran_VERSION : " ${MPI_Fortran_VERSION} )
get_filename_component(MPI_Fortran_COMPILER_NAME ${MPI_Fortran_COMPILER} NAME)
#message( STATUS "*** OLLA 111 : " ${Fortran_COMPILER_NAME} )
#message( STATUS "*** OLLA 222 : " ${MPI_Fortran_COMPILER} )
################################################
# Other options
# as CMAKE usage needs
################################################

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set( CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE )
endif()

set( CMAKE_INSTALL_PREFIX ${MLD_INSDIR} CACHE STRING
     "Choose the installation directory; by default it installs in the MLD_INSTALL_DIR."
     FORCE )


################################################
# set paths mod executables
################################################
set( CMAKE_Fortran_MODULE_DIRECTORY ${MLD_BUILD_DIR}/mod )
include_directories(${MLD_BUILD_DIR}/mod)
set( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin )

# FFLAGS depend on the compiler
get_filename_component( Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME )

# avoiding if elif else endif.
message( STATUS "*** CMAKE_BUILD_TYPE : " ${CMAKE_BUILD_TYPE} )
message( STATUS "*** CMAKE_INSTALL_PREFIX : " ${CMAKE_INSTALL_PREFIX} )
message( STATUS "*** CMAKE_Fortran_COMPILER full path : " ${CMAKE_Fortran_COMPILER} )
message( STATUS "*** Fortran_COMPILER_NAME : " ${Fortran_COMPILER_NAME} )
message( STATUS "*** FC                    : " $ENV{FC} )
message( STATUS "*** CC                    : " $ENV{CC} )
message( STATUS "*** CXX                   : " $ENV{CXX} )
message( STATUS "*** MPI_FC                : " $ENV{MPI_FC} )
message( STATUS "*** MPI_CC                : " $ENV{MPI_CC} )
message( STATUS "*** MPI_CXX               : " $ENV{MPI_CXX} )
message( STATUS "*** MKL : " ${MKL_ROOT} )


if ( ((MPI_Fortran_COMPILER_NAME STREQUAL "mpiifort") OR (MPI_Fortran_COMPILER_NAME STREQUAL "mpiifx"))  AND 
         ((Fortran_COMPILER_NAME STREQUAL "ifx") OR (Fortran_COMPILER_NAME STREQUAL "ifort")) )

  set( CMAKE_COMPILE_MODE "INTEL" )
  if (Fortran_COMPILER_NAME STREQUAL "ifx")
    message( STATUS "***********************************\n   (INTEL) ifx   + intelmpi (" ${CMAKE_COMPILE_MODE} ")\n   ************************************")
  endif()

  if (Fortran_COMPILER_NAME STREQUAL "ifort")
   message( STATUS "***********************************\n   (INTEL) ifort + intelmpi (" ${CMAKE_COMPILE_MODE} ")\n   ************************************")
  endif()
  
  add_compile_definitions(IFORT)
  set( TMP_MKL "${MKL_ROOT}/lib/intel64" )  
  set( TMP_OMP ${MPI_Fortran_MODULE_DIR})  

  set( TTMP "${MKL_ROOT}/lib/intel64" )  
  # OPTIMIZED VERSION: 
  # Added -fp-model strict to fix numerical issues with -O2 optimization in angular_functions.f90
  # Also disable unsafe math optimizations and preserve NaN/Inf handling
  set( CMAKE_Fortran_FLAGS_RELEASE " -O3 -xHost -diag-disable=10448  -traceback -L${TMP_MKL} -L${TMP_OMP} -L${TTMP}" )
  #set( CMAKE_Fortran_FLAGS_RELEASE " -O2 -no-vec -fp-model strict -fpe-all=0 -no-ftz -no-fast-transcendentals -prec-div -prec-sqrt -xHost -diag-disable=10448  -traceback -L${TMP_MKL} -L${TMP_OMP} -L${TTMP}" )
  # Alternative: use O1 globally if file-specific optimization doesn't work
  #set( CMAKE_Fortran_FLAGS_RELEASE " -O1 -no-vec -fp-model strict -fpe-all=0 -no-ftz -xHost -diag-disable=10448  -traceback -L${TMP_MKL} -L${TMP_OMP} -L${TTMP}" )
  #set( CMAKE_Fortran_FLAGS_RELEASE " -O0 -g -traceback  -fstack-security-check -fsanitize=address -fno-omit-frame-pointer -check bounds  -check pointers -check all -init=snan,arrays -ftrapuv -fpe-all=0 -ftz -warn all -warn nousage -heap-arrays -L${TMP_MKL} -L${TMP_OMP}" )

  set( CMAKE_CXX_FLAGS_RELEASE " -O2   -L${TMP_MKL} -L${TMP_OMP}" )
  # MAGICAL DETECT OVERFLOW :
  #set( CMAKE_Fortran_FLAGS_RELEASE "   -O0 -debug all -check all -warn all  -diag-disable 8889 -fpe-all=0 -no-ftz -g -traceback -L${TMP_MKL} -L${TMP_OMP}" )
  # set( CMAKE_Fortran_FLAGS_RELEASE " -O0 -debug all -g -fpe0 -fpe-all=0 -no-ftz -traceback " )

  #set( CMAKE_Fortran_FLAGS_DEBUG " -O0 -debug all -g -fpe0 -traceback -C " )
  include_directories( SYSTEM  ${MPI_Fortran_MODULE_DIR} )

  # Special compilation settings for numerically sensitive files (Intel compilers only)
  # Compile angular_functions.f90 with reduced optimization to prevent numerical issues
  set_source_files_properties(src/mld/angular_functions.f90 PROPERTIES
    COMPILE_FLAGS " -O1  -diag-disable=10448"
  )

  # set it for 'gf' as 'mpifort-gnu' and 'intel' as 'mpiifort-intel'
  list( APPEND MLD_MKL_LIB
    mkl_blas95_lp64
    mkl_scalapack_lp64
    mkl_intel_lp64
    mkl_sequential
    mkl_core
    mkl_blacs_intelmpi_lp64
    pthread
    m
    dl
  )


elseif ( (MPI_Fortran_COMPILER_NAME STREQUAL "mpif90")  AND 
         ((Fortran_COMPILER_NAME STREQUAL "ifx") OR (Fortran_COMPILER_NAME STREQUAL "ifort")) )

  set( CMAKE_COMPILE_MODE "MIX" )
  if (Fortran_COMPILER_NAME STREQUAL "ifx")
   message( STATUS "***********************************\n   (MIX) ifx   + openmpi (" ${CMAKE_COMPILE_MODE} ")\n   ************************************")
  endif()

  if (Fortran_COMPILER_NAME STREQUAL "ifort")
  message( STATUS "***********************************\n   (MIX) ifort + openmpi (" ${CMAKE_COMPILE_MODE} ")\n   ************************************")
  endif()

  add_compile_definitions(IFORT)

  set( TMP_MKL "${MKL_ROOT}/lib/intel64" )  
  set( TMP_OMP ${MPI_Fortran_MODULE_DIR})  
  message( STATUS "*** CMAKE TMP_OMP : " ${TMP_OMP} )


    # Release: optimized, no checks
    #set( CMAKE_Fortran_FLAGS_RELEASE " -O3 -xHost -qopt-zmm-usage=high -traceback -diag-disable=10448 -L${TMP_MKL} -L${TMP_OMP}" )

  set( CMAKE_Fortran_FLAGS "" )
  # Separate flags for ifx vs ifort
  if (Fortran_COMPILER_NAME STREQUAL "ifx")
    # Release: optimized, no checks
    set( CMAKE_Fortran_FLAGS_RELEASE " -O3 -xHost -qopt-zmm-usage=high -traceback -diag-disable=10448 -L${TMP_OMP}" )
    # Debug: bounds check only (explicit list, no uninit → no MSan)
    # Clear any base CMAKE_Fortran_FLAGS to avoid stray -check
    # Debug: space-separated -check (no '='), avoid uninit (prevents MSan)
    # Use multiple -check tokens (robust for ifx)
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback -fpe-all=0 -check bounds -check format -check output_conversion -diag-disable=10448 -L${TMP_OMP}")
  else()
    # ifort: -C is shorthand for all checks, safe (no MSan)
    set( CMAKE_Fortran_FLAGS_RELEASE " -O3 -xHost -qopt-zmm-usage=high -traceback -diag-disable=10448 -L${TMP_OMP}" )
    set( CMAKE_Fortran_FLAGS_DEBUG " -O0 -g -traceback -C -diag-disable=10448 -L${TMP_OMP}" )
  endif()


  #set( CMAKE_Fortran_FLAGS_RELEASE " -O3 -xHost -qopt-zmm-usage=high -traceback  -diag-disable=10448    -L${TMP_MKL} -L${TMP_OMP}" )
  #set( CMAKE_Fortran_FLAGS_RELEASE " -C -O0 -fno-sanitize=memory -g -traceback  -diag-disable=10448    -L${TMP_MKL} -L${TMP_OMP}" )
  #set( CMAKE_Fortran_FLAGS_RELEASE " -O3  -traceback  -diag-disable=10448 -diag-disable 8889  -warn all  -warn nousage -L${TMP_MKL} -L${TMP_OMP}" )
  #set( CMAKE_Fortran_FLAGS_RELEASE " -O2  -traceback  -diag-disable=10448    -L${TMP_MKL} -L${TMP_OMP}" )
  #set( CMAKE_Fortran_FLAGS_RELEASE " -O2 -fprotect-parens -heap-arrays -g -traceback  -L${TMP_MKL} -L${TMP_OMP}" )
  #set( CMAKE_Fortran_FLAGS_RELEASE " -O3  -traceback -L${TMP_MKL} -L${TMP_OMP}" )
  # MAGICAL DETECT OVERFLOW ... # -C is check bounds
  #set( CMAKE_Fortran_FLAGS_RELEASE "  -C -O0 -g -traceback -check uninit -warn all -warn nousage -diag-disable 8889 -diag-disable=10448  -fpe-all=0  -traceback -L${TMP_MKL} -L${TMP_OMP}")
  #set( CMAKE_Fortran_FLAGS_RELEASE "  -C -O0 -g -traceback -check uninit -warn all -warn nousage -diag-disable 8889 -diag-disable=10448  -fpe-all=0  -L${TMP_MKL} -L${TMP_OMP}")
  # set( CMAKE_Fortran_FLAGS_RELEASE " -O2 -fpe0 -traceback -L${TMP_MKL} -L${TMP_OMP}" )

  #set( CMAKE_Fortran_FLAGS_DEBUG " -O0 -debug all -g -fpe0 -traceback -C -L${TMP_MKL} -L${TMP_OMP} " )
  set_source_files_properties(src/mld/angular_functions.f90 PROPERTIES
    COMPILE_FLAGS " -O1  -diag-disable=10448"
  )
  include_directories( SYSTEM  ${MPI_Fortran_MODULE_DIR} )  # get *.mod ifort


  add_compile_definitions(MPI_STATUSES_IGNORE=1)


  list( APPEND MLD_MKL_LIB
    ${MKL_ROOT}/lib/intel64/libmkl_blas95_lp64.a
    ${MKL_ROOT}/lib/intel64/libmkl_scalapack_lp64.a
    "-Wl,--start-group"
    ${MKL_ROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a
    ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.a
    ${MKL_ROOT}/lib/intel64/libmkl_sequential.a
    ${MKL_ROOT}/lib/intel64/libmkl_core.a
    "-Wl,--end-group"
    pthread
    iomp5
    m
    dl
  )


  # TODO mpi_mpifh

elseif ( (MPI_Fortran_COMPILER_NAME STREQUAL "mpif90")  AND ( (Fortran_COMPILER_NAME STREQUAL "gfortran")) )
  # gfortran gnu 10 mini 220203
  set( CMAKE_COMPILE_MODE "GNU" )
  message( STATUS "***********************************\n   (GNU) gfortran + openmpi (" ${CMAKE_COMPILE_MODE} ")\n   ***********************************" )
  add_compile_definitions(GNU)

  # Find all .f files in your project source directory and set specific properties
  file(GLOB_RECURSE Fortran77_SOURCES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.f")
  foreach(_file ${Fortran77_SOURCES})
    set_source_files_properties(${_file} PROPERTIES
    COMPILE_FLAGS  " -O3 -ffixed-form -L${TMP_MKL} -L${TMP_OMP} "
   )
  endforeach()

  set( TMP_MKL "${MKL_ROOT}/lib/intel64" )  
  set( TMP_OMP ${MPI_Fortran_MODULE_DIR})  
  message( STATUS "*** CMAKE TMP_OMP : " ${TMP_OMP} )
  #set(CMAKE_Fortran_MODULE_DIRECTORY ${MPI_Fortran_MODULE_DIR}/mod )
  include_directories( SYSTEM ${MPI_Fortran_MODULE_DIR} )

  #set (CMAKE_Fortran_FLAGS_RELEASE "  -funroll-all-loops -fno-f2c -O3 -ffree-form  -ffree-line-length-512 -L${TMP_MKL} -L${TMP_OMP} -J${MLD_BUILD_DIR}/mod")
  set (CMAKE_Fortran_FLAGS_RELEASE "  -funroll-all-loops -fno-f2c -O3 -ffree-form  -ffree-line-length-512 -L${TMP_MKL} -L${TMP_OMP} ")
  # set (CMAKE_Fortran_FLAGS_RELEASE " -pg -g  -fbacktrace -O0 -ffree-form  -ffree-line-length-512 -L${TMP_MKL} -L${TMP_OMP} ")
  set (CMAKE_Fortran_FLAGS_DEBUG "  -pg -g -fbacktrace -fcheck=all -finit-real=snan -ffpe-trap=zero,overflow,invalid -O0 -ffree-line-length-512 -L${TMP_MKL} -L${TMP_OMP} ")
  # -pg for profilling then gprof 
  # set (CMAKE_Fortran_FLAGS_RELEASE "  -pg -g  -fbacktrace -ffpe-trap=denormal -fcheck=all -finit-real=snan -ffpe-trap=zero,overflow,invalid -Wall -O0 -Wno-compare-reals -Wextra -ffree-form  -ffree-line-length-512 -L${TMP_MKL} -L${TMP_OMP} ")
  # set (CMAKE_Fortran_FLAGS_RELEASE "  -pg -g -fbounds-check -fbacktrace -ffpe-trap=denormal -fcheck=all -ffpe-trap=zero,overflow,invalid -Wall -Wno-compare-reals -Wextra -ffree-form  -ffree-line-length-512 -L${TMP_MKL} -L${TMP_OMP}")
  #set (CMAKE_Fortran_FLAGS_RELEASE "  -pg -g  -fbacktrace -ffpe-trap=denormal -fcheck=all -ffpe-trap=zero,overflow,invalid -Wall -Wno-compare-reals -Wextra -ffree-form  -ffree-line-length-512 -L${TMP_MKL} -L${TMP_OMP} -J${MLD_BUILD_DIR}/mod")
  #set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -ffree-form  -ffree-line-length-512 ")

  # set it for 'gf' as 'mpifort-gnu' and 'intel' as 'mpiifort-intel'
  list( APPEND MLD_MKL_LIB
    mkl_blas95_lp64
    mkl_scalapack_lp64
    mkl_gf_lp64
    mkl_sequential
    mkl_core
    mkl_blacs_openmpi_lp64
    pthread
    m
    dl
  )

else()
  # others general else as houston problem
  message( FATAL_ERROR "*** No optimized Fortran compiler flags are known, we just try NOTHING!" )
  # could be
  # set( CMAKE_Fortran_FLAGS_RELEASE "-O2" )
  # set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g" )
endif()

list( APPEND MLD_HDF_LIB ${HDF5_Fortran_LIBRARIES} )

################################################
# preprocessor definitions used in milady 220210
# cf in code '#if(PHONDY || PARAPH || MAB || ML || PARAML)'
# au cas ou, pour une cible specifique,
# vous pouvez effectuer les operations suivantes :
# target_compile_definitions(my_target PRIVATE FOO=1 BAR=1)
################################################

# useless from cosmin 220210
# add_compile_definitions(PHONDY="")
# add_compile_definitions(PARAPH="")
# add_compile_definitions(MAB="")
# add_compile_definitions(PARAM=1)

# Check if the environment variable MLD_NDM is set to ON
if("$ENV{MLD_NDM}" STREQUAL "ON")
  # Add the NDM compile definition
  add_compile_definitions(MLD_NDM=1)
endif()

# Check if the environment variable MLD_NDM is set to ON
if("$ENV{MLD_HDF5}" STREQUAL "ON")
  # Add the HDF5 compile definition
  add_compile_definitions(MLD_HDF5=1)
endif()

add_compile_definitions(ML=1)
add_compile_definitions(PARA=1)


message( STATUS "*** !!!!! CMAKE_Fortran_FLAGS_RELEASE : " ${CMAKE_Fortran_FLAGS_RELEASE} )
message( STATUS "*** !!!!! CMAKE_Fortran_FLAGS_DEBUG   : " ${CMAKE_Fortran_FLAGS_DEBUG} )
message( STATUS "*** !!!!! Module build dir            : " ${MLD_BUILD_DIR}/mod )





if( NOT MLD_OPT_COMPILE_ONLY_CTESTS )
  add_library( MILADY STATIC ${test_SOURCES} )
  add_subdirectory( src )
  install( TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
endif()

#TODOhdf5  add_subdirectory( tests )

################################################
# DOXYGEN Documentation
################################################
#option(BUILD_DOC "Build documentation" ON)
#include( doxygen_doc )

# optionally get trace of all cmake main variables in log
if( MLD_OPT_TRACE )
  printCmakeTrace()
endif()
