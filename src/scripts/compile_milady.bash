#!/usr/bin/env bash

# milady compilation etc, as bash functions
# cd .../MILADY/src/scripts && source ./compile_milady.bash

source ${MLD_SRCDIR}/src/scripts/utils_milady.bash
# set -x for script debug



unset -f f_set_mld_buildir
#$# function f_set_mld_buildir {
#$#   # append ${1} as "intel" "gnu" "mix" to $MLD_BUIDIR
#$#   # (storing original $MLD_BUIDIR in $MLD_BUIDIR_ORIG if not done)
#$#   unset MLD_BUIDIR_ORIG
#$#   [ -z "${MLD_BUIDIR_ORIG}" ] && export MLD_BUIDIR_ORIG=${MLD_BUIDIR}
#$#   export MLD_BUIDIR=${MLD_BUIDIR_ORIG}"_"${1}  # example from ../mld_build to  .../mld_build_mix
#$# }

function f_set_mld_buildir {
  # Check and store the original $MLD_BUIDIR once
  if [ -z "${MLD_BUIDIR_ORIG}" ]; then
    export MLD_BUIDIR_ORIG=${MLD_BUIDIR}
  fi

  # Construct the new directory name with the desired suffix
  local new_dir="${MLD_BUIDIR_ORIG}_${1}"

  # Check if MLD_BUIDIR already ends with the desired suffix
  if [[ "${MLD_BUIDIR}" != *"${1}" ]]; then
    export MLD_BUIDIR=${new_dir}
  fi
  # Now, MLD_BUIDIR will only change if it doesn't already end with the specified suffix
}


function f_clean_dir {
  # clean dir contents, no remove dir if existing, warning if $1 is file
  in_green "Clean directory ${1}"
  if [ ! -z "${1}" ]; then  # avoid 'rm /*'
    [ -d "${1}" ] && rm -rf ${1}/*
    [ -f "${1}" ] && in_red "expected directory ${1} to clean, but is a file"
    [ ! -e "${1}" ] && in_green "inexisting dir ${1}, clean is superfluous"
  else
    in_red 'f_clean_dir expected directory name ($1) empty, fix it.'
  fi
}

function f_clean_milady {
  # manual clean MILADY source directory
  [ -z "${MLD_SRCDIR}" ] && in_red 'problem undefined ${MLD_SRCDIR}'
  if [ ! -z "${MLD_SRCDIR}" ]; then

    # in case of (eventually) cmake pollution in MILADY source directory
    rm -rf ${MLD_SRCDIR}/mod
    rm -rf ${MLD_SRCDIR}/lib
    rm -rf ${MLD_SRCDIR}/build
    rm -rf ${MLD_SRCDIR}/CMakeFiles
    rm -f  ${MLD_SRCDIR}/CMakeCache.txt ${MLD_SRCDIR}/cmake_install.cmake
    rm -f  ${MLD_SRCDIR}/Makefile ${MLD_SRCDIR}/install_manifest.txt
    rm -f  ${MLD_SRCDIR}/libMILADY.a

    # classical clean all case
    [ ! -z "${MLD_INSDIR}" ] && rm -rf ${MLD_INSDIR} || in_red 'problem undefined ${MLD_INSDIR}'
    if [ ! -z "${MLD_BUIDIR_ORIG}" ]; then  # stay bug TODO build/mod
      f_clean_dir ${MLD_BUIDIR}
    fi
    if [ ! -z "${MLD_BUIDIR}" ]; then
      f_clean_dir ${MLD_BUIDIR}
    else
      in_red 'problem undefined ${MLD_BUIDIR}'
    fi
  fi
}

# Add after f_clean_milady
function f_deep_clean_objs {
  in_green "Deep clean object/ module files under source tree"
  find "${MLD_SRCDIR}" -type f \( -name '*.o' -o -name '*.mod' -o -name '*.smod' \) -delete
}

function f_cmake_flags_fortran {
  local mode=${1:-Release}
  # Common libs
  local mkl="-L${MKL_ROOT}/lib/intel64"
  local omp="-L${OMP_INSDIR}/lib"
  if [[ "${FC}" == *ifx ]]; then
    if [[ "${mode}" == "Debug" ]]; then
      echo "-O0 -g -traceback -C -check:no-uninit -fno-sanitize=memory -diag-disable=10448 ${mkl} ${omp}"
    else
      echo "-O3 -xHost -qopt-zmm-usage=high -traceback -diag-disable=10448 -fno-sanitize=memory ${mkl} ${omp}"
    fi
  else
    # ifort / gfortran
    if [[ "${mode}" == "Debug" ]]; then
      echo "-O0 -g -traceback -C -diag-disable=10448 ${mkl} ${omp}"
    else
      echo "-O3 -xHost -traceback -diag-disable=10448 ${mkl} ${omp}"
    fi
  fi
}


function f_explore_milady {
  # list libraries dependencies of executables, for verification
  in_red '\ndirectory ${MLD_SRCDIR}' ${MLD_SRCDIR}
  ls -alt ${MLD_SRCDIR}
  in_red '\ndirectory ${MLD_BUIDIR}' ${MLD_BUIDIR}
  ls -alt ${MLD_BUIDIR}
  tree -d ${MLD_BUIDIR}
  in_red '\ndirectory ${MLD_INSDIR}' ${MLD_INSDIR}
  tree ${MLD_INSDIR}
  f_ldd_exe
}

function f_ldd_exe {
  # $1 could be -v as verbose
  tmp=$(find ${MLD_BUIDIR} -name "*.exe")
  for i in $tmp; do
    in_green '\nldd '${1} ${i}
    ldd $1 $i
  done
}

function f_end_mld_compile {
  cd ${MLD_BUIDIR}
  in_green 'End build in directory '${MLD_BUIDIR}
  ls -alt ${MLD_BUIDIR}
  in_green '\nNow you could type:\n
  make -j'${NPROC}'\n
  make install\n
  f_ctest_milady\n
  f_explore_milady\n
  '
}

function f_compile_milady {
  # do the MILADY cmake compile job with remove all previous build
  # if 'f_compile_milady ON' then cmake trace is activated
  envs | grep 'MLD_'

  # standart entry point for cmake find hdf5
  # [ -z "${HDF5_ROOT}" ] && in_red 'problem undefined ${HDF5_ROOT} (try ${MLD_ROODIR}/hdf_install)'

  [ -z "${MLD_SRCDIR}" ] && in_red 'problem undefined ${MLD_SRCDIR}'
  [ -z "${MLD_BUIDIR}" ] && in_red 'problem undefined ${MLD_BUIDIR}'
  [ -z "${MLD_INSDIR}" ] && in_red 'problem undefined ${MLD_INSDIR}'
  [ -z "${FC}" ] && in_red 'problem undefined ${FC}'
  [ -z "${CC}" ] && in_red 'problem undefined ${CC}'
  [ -z "${CXX}" ] && in_red 'problem undefined ${CXX}'
  local BUILD_TYPE=${CMAKE_BUILD_TYPE:-Release}
  TRAC=${1:-OFF} # OFF by default as first parameter of function
  NPROC=${NPROC:-$(f_nproc)}
  # envs -g HDF_
  # envs -g OMP_
  if [ ! -z "${MLD_SRCDIR}" ]; then
    in_green "BUILD_TYPE=${BUILD_TYPE} TRACE=${TRAC}"
    f_clean_milady
    # Deep clean objects/modules in source tree (avoid MSan-instrumented leftovers)
    find "${MLD_SRCDIR}" -type f \( -name '*.o' -o -name '*.mod' -o -name '*.smod' \) -delete
    in_green FC=${FC}; in_green CC=${CC}; in_green CXX=${CXX}
    cmd="$(f_cmake3) -S${MLD_SRCDIR} -B${MLD_BUIDIR} \
         -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
         -DCMAKE_Fortran_COMPILER=${FC} \
         -DCMAKE_C_COMPILER=${CC} \
         -DCMAKE_CXX_COMPILER=${CXX} \
         -DMLD_OPT_TRACE=${TRAC}"
    #cmd="$(f_cmake3) -S${MLD_SRCDIR} -B${MLD_BUIDIR}  -DCMAKE_Fortran_COMPILER=${FC}  -DMLD_OPT_TRACE=${TRAC}"

    in_green ${PWD}
    in_green ${cmd}
    ${cmd} && f_end_mld_compile || in_red 'problem cmake command'
    in_green "END f_compile_milady"
  fi
}

function f_compile_milady_intel {
  in_green "compile milady with intel mpiifort and intel mpi"
  #old export FC=$(f_which mpiifort)
  #old export CC=$(f_which icc)
  #export FC=$(f_which ifort)
  export FC=$(f_which ifx)
  export CC=$(f_which mpiicx)
  export CXX=$(f_which mpiicpx)
  export MPI_FC=$(f_which mpiifx)
  #export MPI_FC=$(f_which mpiifx)
  export MPI_CC=$(f_which mpiicx)
  export MPI_CXX=$(f_which mpiicpx)
  f_set_mld_buildir intel
  f_compile_milady
}


function f_compile_milady_gnu {
  in_red "compile milady with gnu mpifort and openmpi\n
  this is tricky as link to oneApi mkl libraries for gnu gf"
  export FC=$(f_which gfortran)
  [ "${FC}" = "UNKNOWN" ] && f_help_compile_gnu
  export CC=$(f_which gcc)
  export CXX=$(f_which g++)
  export MPI_FC=$(f_which mpif90)
  export MPI_CC=$(f_which mpicc)
  export MPI_CXX=$(f_which mpic++)
  f_set_mld_buildir gnu
  export MPI_HOME=${OMP_INSDIR}
  export OMP_ROOT=${OMP_INSDIR} 
  export MKLROOT=${MKL_ROOT}
  [ -z "${MKL_ROOT}" ] && in_red 'problem undefined ${MKL_ROOT}' || f_compile_milady
}

function f_compile_milady_mix_ifx {
  in_green "compile milady (mix) with ifx + openmpi"
  [ -z "${MKL_ROOT}" ] && in_red 'problem undefined ${MKL_ROOT}'
  [ -z "${OMP_INSDIR}" ] && in_red 'problem undefined ${OMP_INSDIR}'
  export FC=$(f_which ifx)
  export MPI_FC=$(f_which mpif90)
  export CC=$(f_which gcc)
  export MPI_CC=$(f_which mpicc)
  export MPI_CXX=$(f_which mpic++)
  export CXX=$(f_which g++)
  f_set_mld_buildir mix
  export MPI_HOME=${OMP_INSDIR}
  export OMP_ROOT=${OMP_INSDIR}    # /usr/local/iopenmpi.4.1.2/lib64/
  # Set build type (Release or Debug)
  export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:-Release}
  in_red 'TODO have to set other way ${OMP_ROOT} ? ' ${OMP_ROOT}
  in_green "type: f_compile_milady"
  f_compile_milady # TODO uncomment
}


function f_compile_milady_mix_ifort {
  in_green "compile milady with intel ifort and intel-compiled openmpi"
  [ -z "${MKL_ROOT}" ] && in_red 'problem undefined ${MKL_ROOT}'
  [ -z "${OMP_INSDIR}" ] && in_red 'problem undefined ${OMP_INSDIR}'
  export FC=$(f_which ifort)
  #export FC=$(f_which ifx)
  export MPI_FC=$(f_which mpif90)
  #export CC=$(f_which icc)
  export CC=$(f_which gcc)
  export MPI_CC=$(f_which mpicc)
  export MPI_CXX=$(f_which mpic++)
  export CXX=$(f_which g++)
  f_set_mld_buildir mix
  export MPI_HOME=${OMP_INSDIR}
  export OMP_ROOT=${OMP_INSDIR}    # /usr/local/iopenmpi.4.1.2/lib64/
  in_red 'TODO have to set other way ${OMP_ROOT} ? ' ${OMP_ROOT}
  in_green "type: f_compile_milady"
  f_compile_milady # TODO uncomment
}


function f_make {
  [ -z "${MLD_BUIDIR}" ] && in_red 'problem undefined ${MLD_BUIDIR}'
  cd ${MLD_BUIDIR}
  make -j10
}


function f_help_ctest_milady {
 in_red "
  TODO 220405\n
  This help is for user local computers, NOT irene (where modules load do stuff)\n
  You have to set PATH and LD_LIBRARY_PATH before launch tests.\n
  this is tricky, may be TODO a user-defined module load file.\n
  In case of using openmpi intel_compiled (mix installation, or else gnu and intel)\n
  may be (as example)\n"

 echo -e '
  # For user local computers, NOT irene (where modules load do stuff)
  export MLD_ROODIR=/home/catA/${USER}/MLD
  # OAP_ as OneAPi
  export OAP_ROODIR=/volatile2/catA/${USER}/oneAPI
  # OMP_ as OpenMPi
  export OMP_INSDIR=${MLD_ROODIR}/omp_install_intel
  # set LD_LIBRARY_PATH
  export LD_LIBRARY_PATH_ORIG=${LD_LIBRARY_PATH}
  export LD_LIBRARY_PATH=${OAP_ROODIR}/oneapi/mkl/2022.0.1/lib/intel64:${OMP_INSDIR}/lib:${LD_LIBRARY_PATH}
  export LD_LIBRARY_PATH=${OAP_ROODIR}/oneapi/itac/2021.5.0/bin/rtlib:${LD_LIBRARY_PATH}
  export LD_LIBRARY_PATH=${OAP_ROODIR}/oneapi/compiler/2022.0.1/linux/compiler/lib/intel64_lin:${LD_LIBRARY_PATH}
  # set mpirun PATH
  export PATH_ORIG=${PATH}
  export PATH=${OMP_INSDIR}/bin:${PATH}
  envs -e LD_LIBRARY_PATH
  envs -e PATH
  conda env list
  conda activate yourPythonEnv  # python3 with numpy, matplotlib
  f_ctest_milady
  f_python_tests_milady
  '
}

function f_ctest_milady {
  # launch ctest milady
  # could type 'f_ctest -V' as verbose
  f_help_ctest_milady
  f_setenv_milady
  envs -e LD_LIBRARY_PATH
  envs -e PATH
  envs -g oneapi/mkl
  olddir=$(pwd)
  # could type $1 as f_ctest_milady -V as verbose or f_ctest_milady '-I 11,13'
  if [ -z "${MLD_BUIDIR}" ]; then
    in_red 'problem undefined ${MLD_BUIDIR}'
  else
    cd ${MLD_BUIDIR} || in_red 'problem inexisting directory ${MKL_BUIDIR} '${MLD_BUIDIR}
    ctest ${1}
  fi
  cd ${olddir}
}

function f_activate_python3_marinica {
  # cosmin anaconda3 installation to get python3 with numpy etc
  export MARINICA_ANACONDA="/ccc/work/cont002/den/marinica/anaconda3"
  export PATH="${MARINICA_ANACONDA}/bin:${PATH}"
  source ${WAM_HOME}/conda_activate.bash
  in_green "*** python anaconda (marinica) ***"
  which python
  # conda info --envs
  # list installed packages (as examples)
  # ls ${MARINICA_ANACONDA}/lib/python3.7/site-packages | grep Qt
  # conda list | grep numpy
  python -c "import numpy" || in_red 'problem needs python3 with numpy, at least, try conda info --envs'
}

function f_python_tests_milady {
  # launch python tests milady
  # could do 'ptes_begin=7 ptes_end=9 f_python_tests_milady' for only test 7 to 9
  export mld_plotFailed=${1:-ON}
  python3 -c "import numpy" || in_red 'Needs python3 with numpy, at least, try conda info --envs'
  [ -z "${MLD_TESDIR}" ] && in_red 'problem undefined ${MLD_TESDIR}'
  if [ -z "${MLD_SRCDIR}" ]; then
    in_red 'problem undefined ${MLD_SRCDIR}, fix it'
    in_red 'may be you have to call f_setenv_milady'
  else
    olddir=$(pwd)
    cd ${MLD_TESDIR} || in_red 'problem inexisting ${MLD_TESDIR}'
    # /volatile2/wambeke/ttmp/Tests  # ...where you put the tests
    in_green MLD_TESDIR=${MLD_TESDIR}
    # export MLD_TESPY=${MLD_SRCDIR}/tests/small_tests.py # obsolete
    export MLD_TESPY=${MLD_SRCDIR}/miladypy/AllTestLauncherMILADYPY.sh
    in_green "execute tests with ${MLD_TESPY}"
    # python3 ${MLD_TESPY} clean
    # python3 ${MLD_TESPY} run
    ${MLD_TESPY}
    cd ${olddir}
  fi
}


function f_help_end_compile {
  in_green '\nNow you could type:'
  [ -z "${MODULES_COLLECTION_TARGET}" ] && MODULES_COLLECTION_TARGET="UNKNOWN"
  if [ ${MODULES_COLLECTION_TARGET} = "irene" ]; then
    in_green '  f_compile_milady_irene_mix_ifx or'
    in_green '  (TODO) f_compile_milady_irene_intel or'
    in_green '  (TODO) f_compile_milady_irene_gnu'
  else
    in_green '  f_setenv_milady'
    in_green '  f_compile_milady_mix_ifort or'
    in_green '  f_compile_milady_mix_ifx or'
    in_green '  f_compile_milady_intel or'
    in_green '  f_compile_milady_gnu'
  fi

  in_green '\nAnd after build, for tests, you could type:'
  in_green '  f_ctest_milady and f_python_tests_milady\n'
}

source utils_milady.bash
f_help_end_compile
