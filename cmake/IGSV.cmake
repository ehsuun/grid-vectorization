set(CMAKE_CXX_STANDARD 14) # note: we cannot use c++17 with the current version of gmm
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(OpenGL_GL_PREFERENCE GLVND)

#### compiler-specific options
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  # using Clang
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "./")
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  # force colored output with ninja:
  # https://medium.com/@alasher/colored-c-compiler-output-with-ninja-clang-gcc-10bfe7f2b949
  add_compile_options(-fcolor-diagnostics)

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "./")
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  # force colored output with ninja, see the above link
  add_compile_options(-fdiagnostics-color=always)

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR})
  set(Boost_USE_STATIC_LIBS   ON)
  set(CMAKE_CXX_FLAGS "/Zc:__cplusplus")
  set(CMAKE_CXX_FLAGS "/EHsc")
  set(MSVC_RUNTIME "dynamic")
  add_compile_definitions(_USE_MATH_DEFINES)
  add_compile_definitions(_ENABLE_EXTENDED_ALIGNED_STORAGE)
  add_compile_definitions(WIN32)
endif()
