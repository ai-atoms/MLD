#!/usr/bin/env bash

# openmpi compilation etc, as bash functions
# cd .../MILADY/scripts && source ./compilation_openmpi.bash
# see https://www.open-mpi.org/faq/?category=building#easy-build

source utils_milady.bash

function f_setenv_openmpi_wambeke {
  # as van wambeke example of openmpi compilation directories: all under $OMP_ROODIR
  if [ -z "${OMP_SETENV}" ]; then  # do only one time
    # export OMP_ROODIR=/volatile2/${USER}/MLD  # OMP_ROODIR idem MLD_ROODIR by choice
    # export OMP_ROODIR=/home/catA/${USER}/MLD
    export OMP_ROODIR=${MLD_ROODIR}
    # export OMP_SRCDIR=${OMP_ROODIR}/openmpi-4.1.2
    export OMP_SRCDIR=${OMP_ROODIR}/openmpi-3.1.6
    # useless autotools export OMP_BUIDIR=${OMP_ROODIR}/omp_build
    export OMP_INSDIR=${OMP_ROODIR}/omp_install
    export OMP_SETENV=ON
  fi
}

function f_setenv_openmpi_412_wambeke {
  # as van wambeke example of openmpi compilation directories: all under $OMP_ROODIR
  export OMP_ROODIR=${MLD_ROODIR}
  export OMP_SRCDIR=${OMP_ROODIR}/omp/openmpi-4.1.2
  export OMP_INSDIR=${OMP_ROODIR}/omp/install_412
  export OMP_SETENV=ON
  envs | grep 'OMP_'
}



function f_setenv_openmpi_marinica {
  # as marinica example of openmpi compilation directories
  if [ -z "${MLD_SETENV}" ]; then  # do only one time
    in_red "f_setenv_openmpi_marinica TODO"
  fi
}

function f_setenv_openmpi {
  # for now default is as wambeke
  # could define some f_setenv for ${USER} and ${HOSTNAME} etc
  unset OMP_SETENV
  [ "${USER}" = "wambeke" ] && f_setenv_openmpi_wambeke
  [ "${USER}" = "marinica" ] && f_setenv_openmpi_marinica
  [ -z "${OMP_SETENV}" ] && in_red 'problem unknown ${USER}, may be fix your preferences in f_setenv_openmpi_'${USER}
  [ -z "${OMP_SETENV}" ] && f_setenv_openmpi_wambeke
  envs | grep 'OMP_'
}

function f_wget_openmpi_412 {
  if [ ! -d "${OMP_SRCDIR}" ]; then
    in_red "download openmpi 412 in "${OMP_ROODIR}
    cd ${OMP_ROODIR}
    wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.2.tar.gz
    in_red 'detar ...take a time, be patient...'
    tar -xf ./openmpi-4.1.2.tar.gz
  fi
}

function f_wget_openmpi_316 {
  if [ ! -d "${OMP_SRCDIR}" ]; then
    in_red "download openmpi 316 in "${OMP_ROODIR}
    cd ${OMP_ROODIR}
    wget https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.6.tar.gz
    in_red 'detar ...take a time, be patient...'
    tar -xf ./openmpi-3.1.6.tar.gz
  fi
}

function f_explore_openmpi {
  # list libraries dependencies of executables, for verification
  in_red '\ndirectory ${OMP_SRCDIR}' ${OMP_SRCDIR}
  ls -alt ${OMP_SRCDIR}
  in_red '\ndirectory ${OMP_INSDIR}' ${OMP_INSDIR}
  tree -d ${OMP_INSDIR}
  in_red '\ndirectory ${OMP_INSDIR}/bin' ${OMP_INSDIR}/bin
  ls -alt ${OMP_INSDIR}/bin
}

function f_compile_openmpi {
  [ -z "${OMP_SRCDIR}" ] && in_red 'problem undefined ${OMP_SRCDIR}'
  # useless autotools [ -z "${OMP_BUIDIR}" ] && in_red 'problem undefined ${OMP_BUIDIR}'
  [ -z "${OMP_INSDIR}" ] && in_red 'problem undefined ${OMP_INSDIR}'
  [ -z "${FC}" ] && in_red 'problem undefined ${FC}'
  [ -z "${CC}" ] && in_red 'problem undefined ${CC}'
  f_wget_openmpi_316
  NPROC=4  # no more ... ${NPROC:-$(f_nproc)}
  if [ ! -z "${OMP_SRCDIR}" ]; then
    cd ${OMP_SRCDIR}
    in_green 'FC='${FC}
    in_green 'F77='${F77}
    in_green 'CC='${CC}
    in_green 'CXX='${CXX}
    in_green 'OMP_INSDIR='${OMP_INSDIR}
    in_red 'configure ...take a time, be patient...'
    # ./configure CC=icc CXX=icpc F77=ifort FC=ifort --prefix=${OMP_INSDIR}
    ./configure CC=${CC} CXX=${CXX} F77=${F77} FC=${FC} --prefix=${OMP_INSDIR}
    in_red 'make -j'${NPROC}' ...take a time, be patient...'
    make -j${NPROC} all
    make install
    in_green "END f_compile_openmpi" ${OMP_INSDIR}
  fi
}

function f_set_omp_instdir {
  # append ${1} as "intel" "gnu" to $OMP_INSDIR
  # (storing original $OMP_INSDIR in $OMP_INSDIR_ORIG if not done)
  [ -z "${OMP_INSDIR_ORIG}" ] && export OMP_INSDIR_ORIG=${OMP_INSDIR}
  export OMP_INSDIR=${OMP_INSDIR_ORIG}"_"${1}  # example from ../omp_install to  .../omp_install_gnu
}

function f_compile_openmpi_intel {
  in_green "compile openmpi with intel ifort"
  export FC=$(f_which ifort)
  export F77=$(f_which ifort)
  export CC=$(f_which icc)
  export CXX=$(f_which icpc)
  f_set_omp_instdir intel
  f_compile_openmpi
}

function f_compile_openmpi_gnu {
  in_green "compile openmpi with gnu gfortran, verify environment without oneAPI"
  export FC=$(f_which gfortran)
  export F77=$(f_which gfortran)
  export CC=$(f_which cc)
  export CXX=$(f_which g++)
  f_set_omp_instdir gnu
  envs -g oneAPI
  in_red "verify environment without oneAPI please"
  f_compile_openmpi
}

in_green "Now you could type:\n
  f_setenv_openmpi\n
  f_compile_openmpi_intel\n
  f_compile_openmpi_gnu\n
"
