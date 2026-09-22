#!/usr/bin/env bash

# common milady bash utilities, as bash functions
# cd .../MILADY/src/scripts && source src/scripts/utils_milady.bash


function in_red () {
  # Red bold echo
  echo -e "\e[31m\e[1m"$@"\e[0m"
}

function in_green () {
  # Green bold echo
  echo -e "\e[32m\e[1m"$@"\e[0m"
}

function f_which {
  # which returning 'UNKNOWN' if not found
  # avoid stderr message if not found
  tmp=$(which $1 2> /dev/null) && echo $tmp || echo "UNKNOWN"
}


function f_nproc {
  # policy is no use all procs to freeze desktop.
  if [ ${MODULES_COLLECTION_TARGET} = "irene" ]; then
    echo 10  # faster compilation on irene, TODO ... or more ?
  else
    tmp=$(nproc) || tmp=8
    case "$tmp" in
      '1')
        echo "1" ;;
      '2')
        echo "2" ;;
      '3')
        echo "3" ;;
      '4')
        echo "4" ;;  
      *)
        # echo -e ${tmp}'/2' | bc  # /2 for cores/cpu ?
        echo -e ${tmp}'-2' | bc  # -2 for no freeze desktop ;;
    esac
  fi
}

function f_cmake3 {
  tmp=$(\cmake --version)    # "\" avoid alias, caseof
  case "$tmp" in
    *'version 2'*)
      echo cmake3 ;;
    *'version 3'*)
      echo cmake ;;
    *'version 4'*)
      echo cmake ;;  
    *)
      echo "UNKNOWN cmake" ;;
  esac
}
