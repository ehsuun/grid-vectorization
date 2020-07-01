cmake_minimum_required(VERSION 3.1)
project(ALGLIB)
add_library(ALGLIB STATIC
	src/alglibinternal.cpp
	src/alglibmisc.cpp
	src/ap.cpp
	src/dataanalysis.cpp
	src/diffequations.cpp
	src/fasttransforms.cpp
	src/integration.cpp
	src/interpolation.cpp
	src/linalg.cpp
	src/optimization.cpp
	src/solvers.cpp
	src/specialfunctions.cpp
	src/statistics.cpp)
target_include_directories(ALGLIB PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/..")
set_target_properties(ALGLIB PROPERTIES POSITION_INDEPENDENT_CODE ON)
