project(xperience)

cmake_minimum_required(VERSION 2.6.2)

set(INCLUDE_INSTALL_DIR
    "${CMAKE_INSTALL_PREFIX}/include/xperience"
    CACHE INTERNAL
    "Header file install directory"
  )

set(LIBRARY_INSTALL_DIR
    "${CMAKE_INSTALL_PREFIX}/lib"
    CACHE INTERNAL
    "Library file install directory"
  )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

include_directories(src)
add_subdirectory(src lib)
add_subdirectory(test test/bin)

find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)

