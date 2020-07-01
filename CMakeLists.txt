cmake_minimum_required (VERSION 3.9)
if (WIN32)
    cmake_minimum_required (VERSION 3.13)
endif()

project(IGSV)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(IGSV_EXT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext)
set(IGSV_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

include(IGSV)

# -------------------------------------------------------------------
# IGSV options
# -------------------------------------------------------------------
option(IGSV_COMPILE_GUI "Compile the IGSV Gui (based on libigl viewer)" ON)

# -------------------------------------------------------------------
# Boost (required)
# -------------------------------------------------------------------
find_package(Boost 1.48 REQUIRED program_options graph filesystem)

# -------------------------------------------------------------------
# OpenCV (required)
# -------------------------------------------------------------------
find_package(OpenCV REQUIRED)

# -------------------------------------------------------------------
# OpenMP (optional) [UNUSED]
# -------------------------------------------------------------------
# note for MacOS with clang:
# if OpenMP is not found, try specifying the location before executing cmake, e.g.:
#   >> export LDFLAGS="-L/usr/local/opt/llvm/lib"
#   >> cmake .. -DCMAKE_BUILD_TYPE=Release -GNinja -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++
#find_package(OpenMP) # check OpenMP_CXX_FOUND

# -------------------------------------------------------------------
# add other dependencies: Eigen, libigl, CoMISo, SimpleSVG, ALGLIB
# -------------------------------------------------------------------
add_subdirectory(ext)

# -------------------------------------------------------------------
# add IGSV sources
# -------------------------------------------------------------------
add_subdirectory(src)
