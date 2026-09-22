#-------------------------
# MiladyNewMacros.cmake
#-------------------------


# function to print all current variables contents
function(printCmakeTrace)
  get_cmake_property(_variableNames VARIABLES)
  foreach(_variableName ${_variableNames})
    message(STATUS ".. ${_variableName}=${${_variableName}}")
  endforeach()
endfunction()




