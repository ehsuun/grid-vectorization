# This is a minimalist cmake file designed to be used
# with libigl. If you need CoMISo for other purposes
# please use the CMakeLists_original.cmake which
# supports all the features of CoMISo

# modified for IGSV by Tibor Stanko

cmake_minimum_required (VERSION 3.9)
if (WIN32)
  cmake_minimum_required (VERSION 3.13)
endif()

project(CoMISo)

add_library(CoMISo STATIC
        Base/Debug/DebCallStack.cc
        Base/Debug/DebConfig.cc
        Base/Debug/DebFile.cc
        Base/Debug/DebIndexMeshOut.cc
        Base/Debug/DebStream.cc
        Base/Progress/ProgressNode.cc
        Base/Utils/BaseError.cc
        Base/Utils/Environment.cc
        Base/Utils/FileOutput.cc
        Base/Utils/IOutputStream.cc
        Base/Utils/OStringStream.cc
        Base/Utils/StopWatch.cc
        Solver/ConstrainedSolver.cc
        Solver/EigenLDLTSolver.cc
        Solver/GMM_Tools.cc
        Solver/IterativeSolverT.cc
        Solver/MISolver.cc
        NSolver/BoundConstraint.cc
        NSolver/COMISOSolver.cc
        NSolver/FiniteElementProblem.cc
        NSolver/LinearConstraint.cc
        NSolver/NProblemInterface.cc
        NSolver/NPDerivativeChecker.cc
        Utils/CoMISoError.cc
)

target_include_directories(CoMISo PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/"
  "${CMAKE_CURRENT_SOURCE_DIR}/.."
  "${CMAKE_CURRENT_SOURCE_DIR}/../eigen"
  "${CMAKE_CURRENT_SOURCE_DIR}/ext/gmm-5.4/include")


target_compile_definitions(CoMISo PUBLIC -DINCLUDE_TEMPLATES)

# Generate position independent code
set_target_properties(CoMISo PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(WIN32)
  target_compile_definitions(CoMISo PRIVATE -D_USE_MATH_DEFINES -DNOMINMAX)
  target_compile_options(CoMISo PUBLIC /D_SCL_SECURE_NO_DEPRECATE)
  add_library(CoMISo_blas STATIC IMPORTED)
  set_target_properties(CoMISo_blas PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
      IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/ext/OpenBLAS-v0.2.14-Win64-int64/lib/libopenblas.dll.a.lib")
  target_link_libraries(CoMISo PRIVATE CoMISo_blas)
endif()

if(APPLE)
  find_library(accelerate_library Accelerate)
  target_link_libraries(CoMISo PUBLIC ${accelerate_library})
elseif(UNIX)
  find_package(BLAS REQUIRED)
  target_link_libraries(CoMISo PUBLIC ${BLAS_LIBRARIES})
endif()
