#!/usr/bin/env bash

# openmpi compilation etc, as bash functions
# cd .../MILADY/scripts && source ./compilation_openmpi.bash
# see https://www.hdfgroup.org/downloads/hdf5/source-code

function f_setenv_hdf5_wambeke {
  # as van wambeke example of hdf5 compilation directories: all under $HDF_ROODIR
  if [ -z "${HDF_SETENV}" ]; then  # do only one time
    # export HDF_ROODIR=/volatile2/${USER}/MLD  # HDF_ROODIR idem MLD_ROODIR by choice
    export HDF_ROODIR=${MLD_ROODIR} # HDF_ROODIR idem MLD_ROODIR by choice
    # export HDF_SRCDIR=${HDF_ROODIR}/CMake-hdf5-1.12.1/hdf5-1.12.1
    export HDF_SRCDIR=${HDF_ROODIR}/hdf5-1.12.1
    export HDF_BUIDIR=${HDF_ROODIR}/hdf_build
    export HDF_INSDIR=${HDF_ROODIR}/hdf_install
    export HDF_SETENV=ON
  fi
}

function f_setenv_hdf5 {
  # for now default is as wambeke
  # could define some f_setenv for ${USER} and ${HOSTNAME} etc
  unset HDF_SETENV
  [ "${USER}" = "wambeke" ] && f_setenv_hdf5_wambeke
  [ "${USER}" = "marinica" ] && f_setenv_hdf5_marinica
  [ -z "${HDF_SETENV}" ] && in_red 'problem unknown ${USER}, may be fix your preferences in f_setenv_hdf5_'${USER}
  [ -z "${HDF_SETENV}" ] && f_setenv_hdf5_wambeke
  envs | grep 'HDF_'
}

function f_wget_hdf5_1_12_1 {
  if [ ! -d "${HDF_SRCDIR}" ]; then
    in_red "download hdf5 1_12_1 in "${HDF_ROODIR}
    cd ${HDF_ROODIR}
    in_green "see https://www.hdfgroup.org/downloads/hdf5/source-code"
    wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.1/src/hdf5-1.12.1.tar.gz
    in_red 'detar ...take a time, be patient...'
    tar -xf ./hdf5-1.12.1.tar.gz
  fi
}

function f_explore_hdf5 {
  # list libraries dependencies of executables, for verification
  in_red '\ndirectory ${HDF_SRCDIR}' ${HDF_SRCDIR}
  ls -alt ${HDF_SRCDIR}
  in_red '\ndirectory ${HDF_INSDIR}' ${HDF_INSDIR}
  tree -d ${HDF_INSDIR}
  in_red '\ndirectory ${HDF_INSDIR}/bin' ${HDF_INSDIR}/bin
  ls -alt ${HDF_INSDIR}/bin
}

function f_compile_hdf5 {
  [ -z "${HDF_SRCDIR}" ] && in_red 'problem undefined ${HDF_SRCDIR}'
  [ -z "${HDF_BUIDIR}" ] && in_red 'problem undefined ${HDF_BUIDIR}'
  [ -z "${HDF_INSDIR}" ] && in_red 'problem undefined ${HDF_INSDIR}'
  [ -z "${FC}" ] && in_red 'problem undefined ${FC}'
  [ -z "${CC}" ] && in_red 'problem undefined ${CC}'
  NPROC=${NPROC:-$(f_nproc)}
  f_wget_hdf5_1_12_1
  if [ ! -z "${HDF_SRCDIR}" ]; then
    cd ${HDF_ROODIR}
    rm -rf ${HDF_BUIDIR}
    mkdir ${HDF_BUIDIR}
    cd ${HDF_BUIDIR}
    cmd="$(f_cmake3) -S${HDF_SRCDIR} -B${HDF_BUIDIR} -DCMAKE_INSTALL_PREFIX=${HDF_INSDIR} -DHDF5_BUILD_FORTRAN=ON"
    in_green ${cmd}
    in_red '...take a time, be patient...'
    ${cmd}
    # ccmake # for info
    # HDF5_BUILD_FORTRAN ON
    # HDF5_ENABLE_PARALLEL ? FOR FUTURE
    in_red 'make -j'${NPROC}' ...take a time, be patient...'
    make -j${NPROC}
    in_red '...ctest hdf5 unconditionally... take a time, be patient...'
    ctest -j${NPROC}
    make install
    in_green "END f_compile_hdf5"
  fi
}

function f_set_hdf_buildir {
  # append ${1} as "intel" "gnu" to $HDF_BUIDIR
  # (storing original $HDF_BUIDIR in $HDF_BUIDIR_ORIG if not done)
  [ -z "${HDF_BUIDIR_ORIG}" ] && export HDF_BUIDIR_ORIG=${HDF_BUIDIR}
  export HDF_BUIDIR=${HDF_BUIDIR_ORIG}"_"${1}  # example from ../hdf_build to  .../hdf_build_gnu
}

function f_compile_hdf5_intel {
  in_green "compile hdf5 with intel ifort"
  export FC=$(f_which ifort)
  export F77=$(f_which ifort)
  export CC=$(f_which icc)
  export CXX=$(f_which icpc)
  f_set_hdf_buildir intel
  f_compile_hdf5
}

function f_compile_hdf5_gnu {
  in_green "compile hdf5 with gnu gfortran"
  # may be gnu 9 from system
  # module avail gcc
  # module load gcc/9.3.0
  # module avail openmpi
  # module load openmpi/gcc_9.3.0/4.0.1
  export FC=$(f_which gfortran)
  export F77=$(f_which gfortran)
  export CC=$(f_which cc)
  export CXX=$(f_which g++)
  f_set_hdf_buildir gnu
  f_compile_hdf5
}


source utils_milady.bash
in_green 'Now you could type:\n
  f_setenv_hdf5\n
  f_compile_hdf5_intel\n
  f_compile_hdf5_gnu\n
'
