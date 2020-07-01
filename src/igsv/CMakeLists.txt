# -----------------------------------------------------------------------------
# IGSV : core library
# -----------------------------------------------------------------------------
set(IGSV_COMMON_SRC
        common/CustomColors.cpp
        common/Options.cpp
        common/VectorizationData.cpp
        common/print_timings.cpp)

set(IGSV_PARAMETRIZATION_SRC
        parametrization/compute_constraints.cpp
        parametrization/compute_narrow_band.cpp
        parametrization/compute_scale.cpp
        parametrization/compute_snapping_weights.cpp
        parametrization/connect_triangles.cpp
        parametrization/custom_comb_frame_field.cpp
        parametrization/find_nearest.cpp
        parametrization/FrameFieldGenerator.cpp
        parametrization/get_narrow_band_mesh.cpp
        parametrization/IGMSolver.cpp
        parametrization/pix_coords_from_xy.cpp
        parametrization/process_uv.cpp
        parametrization/read_sketch.cpp
        parametrization/simplify_cut_graph.cpp
        parametrization/trace_streamlines.cpp
        parametrization/triangulate.cpp)

set(IGSV_EXTRACTION_SRC
        extraction/assign_pixels_to_qedges.cpp
        extraction/assign_samples_to_chains.cpp
        extraction/connected_components.cpp
        extraction/determine_chain_adjacency.cpp
        extraction/exact_predicates.cpp
        extraction/export_svg.cpp
        extraction/extract_qedges.cpp
        extraction/extract_qports.cpp
        extraction/extract_qvertices.cpp
        extraction/fit_chains.cpp
        extraction/get_dir_as_complex.cpp
        extraction/init_base_labels.cpp
        extraction/init_chains.cpp
        extraction/parametrize_chains.cpp
        extraction/remove_short_chains.cpp
        extraction/store_orientation.cpp
        extraction_entities/Bezier.cpp
        extraction_entities/LineSegment.cpp)

add_library(IGSV_core STATIC
        ${IGSV_COMMON_SRC}
        ${IGSV_PARAMETRIZATION_SRC}
        ${IGSV_EXTRACTION_SRC})

target_include_directories(IGSV_core PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/ext>
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/ext/libigl/external>)

target_compile_definitions(IGSV_core PUBLIC
        CGAL_HEADER_ONLY)

set(IGSV_CORE_EXT_LIBS
        ${OpenCV_LIBS}
        ALGLIB
        Boost::boost
        Boost::graph
        Boost::filesystem
        Boost::program_options
        CoMISo
        igl::cgal
        igl::core
        igl::triangle
        SimpleSVG)

#### OpenMP is not used ####
#if (OpenMP_CXX_FOUND)
#    set(IGSV_CORE_EXT_LIBS
#            ${IGSV_CORE_EXT_LIBS}
#            OpenMP::OpenMP_CXX)
#endif()

target_link_libraries(IGSV_core PUBLIC
        ${IGSV_CORE_EXT_LIBS})

add_library(igsv::core ALIAS IGSV_core)

# -----------------------------------------------------------------------------
# IGSV : gui
# -----------------------------------------------------------------------------
if(IGSV_COMPILE_GUI)
    #### old gui, used for prototyping
#    add_library(
#            IGSV_gui STATIC
#            gui/BaseGui.cpp
#            gui/Gui.cpp
#            gui/Gui__store_entities.cpp)
#    target_link_libraries(
#            IGSV_gui PUBLIC
#            IGSV_core
#            igl::opengl
#            igl::opengl_glfw
#            igl::opengl_glfw_imgui
#            igl::png)

    #### new gui, simple
    add_library(
            IGSV_gui2 STATIC
            gui2/Gui.cpp)
    target_link_libraries(
            IGSV_gui2 PUBLIC
            ${OpenCV_LIBS}
            igl::opengl
            igl::opengl_glfw
            igl::opengl_glfw_imgui)
    target_include_directories(
            IGSV_gui2 PUBLIC
            $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    )
    target_compile_definitions(
            IGSV_gui2 PUBLIC IGSV_GUI2_DIR="${IGSV_SRC_DIR}/igsv/gui2")

    #### which gui to use?
    # add_library(igsv::gui ALIAS IGSV_gui)
    add_library(igsv::gui ALIAS IGSV_gui2)
endif()