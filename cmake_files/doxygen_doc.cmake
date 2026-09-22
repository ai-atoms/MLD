set(DOC_OUTPUTDIR ${MLD_BUILD_DIR}/doc)
set(DOC_LANGUAGE English)
set(QUIET_DOXYGEN YES)

get_target_property(SOURCE_LIST MILADY SOURCES)
#get_target_property(SOURCE_MAIN milady_main.exe SOURCES)
#list(APPEND SOURCE_LIST ${SOURCE_MAIN})
string(REPLACE ";" " " SOURCE_LIST_2 "${SOURCE_LIST}")

find_package(Doxygen)

if (BUILD_DOC)
if (DOXYGEN_FOUND)
    set(DOXYGEN_IN ${PROJECT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${PROJECT_SOURCE_DIR}/Doxyfile)
    
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message(STATUS "Doxygen build started")

    add_custom_target(doc_doxygen ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT} 2> ${DOC_OUTPUTDIR}/err.log
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Generating documentation with Doxygen"
        VERBATIM
    )
else (DOXYGEN_FOUND)
    message(STATUS "Doxygen needs to be installed to generate the Doxygen documentation")
endif (DOXYGEN_FOUND)
endif(BUILD_DOC)
