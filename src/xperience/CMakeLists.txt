file(GLOB xperience_directory_files "*")

foreach(f ${xperience_directory_files})
  if(f MATCHES "\\.hpp")
    list(APPEND xperience_directory_files_to_install ${f})
  endif()
endforeach(f ${xperience_directory_files})

list(APPEND xperience_directory_files_to_install "core")

install(FILES
  ${xperience_directory_files_to_install}
  DESTINATION ${INCLUDE_INSTALL_DIR}
  )

add_library(xperience Scope.cpp Context.cpp Experiment.cpp)

install(TARGETS xperience DESTINATION ${LIBRARY_INSTALL_DIR})

target_link_libraries(xperience config++)