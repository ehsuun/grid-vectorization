# -------------------------------------
# Eigen (required)
# -------------------------------------
set(EIGEN_DIR "${IGSV_EXT_DIR}/eigen")
add_library(eigen INTERFACE)
add_library(Eigen3::Eigen ALIAS eigen)
target_include_directories(eigen INTERFACE
                           $<BUILD_INTERFACE:${EIGEN_DIR}>
                           $<INSTALL_INTERFACE:include>)

# -------------------------------------
# libigl (required)
# -------------------------------------
# disable libigl as static lib to avoid compiling everything
option(LIBIGL_USE_STATIC_LIBRARY            "Use libigl as static lib"  OFF)
option(LIBIGL_WITH_CGAL                     "Use CGAL"                  ON)
option(LIBIGL_WITH_TRIANGLE                 "Use Triangle"              ON)
if(IGSV_COMPILE_GUI)
    option(LIBIGL_WITH_OPENGL               "Use OpenGL"                ON)
    option(LIBIGL_WITH_OPENGL_GLFW          "Use GLFW"                  ON)
    option(LIBIGL_WITH_OPENGL_GLFW_IMGUI    "Use ImGui"                 ON)
    option(LIBIGL_WITH_PNG                  "Use PNG"                   ON)
endif()
find_package(LIBIGL REQUIRED QUIET)

# -------------------------------------
# CoMISo (required)
# -------------------------------------
set(COMISO_DIR "${IGSV_EXT_DIR}/CoMISo")
add_subdirectory(${COMISO_DIR})

# -------------------------------------
# SimpleSVG (required)
# -------------------------------------
add_library(SimpleSVG STATIC simple_svg/simple_svg_1.0.0.cpp)

# -------------------------------------
# ALGLIB (required)
# -------------------------------------
set(ALGLIB_DIR "${IGSV_EXT_DIR}/alglib")
add_subdirectory(${ALGLIB_DIR})