if(MSVC)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR})
else()
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../")
endif()

add_subdirectory(igsv)

add_executable(IGSV_bin main.cpp)

if(TARGET igsv::gui)
  # gui mode
  target_link_libraries(IGSV_bin PUBLIC igsv::core igsv::gui)
  target_compile_definitions(IGSV_bin PUBLIC IGSV_WITH_GUI)
else()
  # command-line mode
  target_link_libraries(IGSV_bin PUBLIC igsv::core)
endif()