#include <igsv/gui/Gui.h>

#include <igsv/common/CustomColors.h>
#include <igsv/common/VectorizationData.h>

#include <boost/filesystem.hpp>

#include <GLFW/glfw3.h>

#include <imgui/imgui.h>

#include <igl/colormap.h>
#include <igl/file_dialog_open.h>
#include <igl/project.h>
#include <igl/unproject.h>

namespace IGSV {

  // ===========================================================================

  /*
        keyboard shortcuts
        ----------------------------------------
        1     [gui]         next image texture
       ⇧1     [gui]         previous image texture
        2     [gui]         toggle uv texture
        3     [v]           Sobel tangents
        4     [e]           (black pixel--nearest vertex) connections
        5     [v]           frame field crosses
        6     [e]           frame field streamlines
       ⇧6     [p]           frame field streamlines
        7     [p]           frame field tangent direction
        8     [e+p]         frame field: period jumps and singularities
        9     [l]           frame field: frame angles
        0
        =     [l]           frame field: # of traced black pixels (for labels)
        A     [l]           pixel labels
        B     [l]           face labels
        C     [l]           curve labels
        D     [p]           non-manifold vertices in the initial narrow band mesh
        E     [keyshort]    compute extraction
        F
        G     [e]           integer-grid map
        H     [p]           q-vertices
        I     [e]           frame field: custom cut graph
        J     [e]           extraction labels
       ⇧J     [e]           (black pixel--sample) connections
        K     [e]           extraction: q-edges
       ⇧K     [e]           extraction: q-ports
        L     [gui]         show mesh wireframe (igl Viewer)
        M     [e]           narrow band boundary
        N     [l]           pixel labels
        O     [gui]         increase line width
       ⇧O     [gui]         increase point size
        P     [gui]         decrease line width
       ⇧P     [gui]         decrease point size
        Q     [keyshort]    compute frame field
        R     [keyshort]    serialize
        S     [gui]         brush points
        T
        U     [e+p]         snapping connections
        V     [l]           vertex labels
        W     [keyshort]    compute parametrization
        X     [e]           fitting: spline curves
       ⇧X     [e]           fitting: spline polygons
        Y
        Z     [e]           fitting: Bezier curves
       ⇧Z     [e]           fitting: Bezier polygons
        ;     [gui]         show vertex ids (igl Viewer)
        :     [gui]         show face ids (igl Viewer)
  */

  IGSVGui::IGSVGui(VectorizationData& data)
      : BaseGui(data),

        E__bezierCurves("Bezier curves", 'Z'),
        E__bezierPolygons("Bezier polygons", 'Z', GLFW_MOD_SHIFT),
        E__blackPixNearestVertex("nearest neighbors", '4'),
        E__customCuts("custom cuts", 'I'),
        E__lineMapSegments("ex: igm edges", 'G'),
        E__narrowBandBoundary("nb bd", 'M'),
        E__periodJumps("period jumps", '8'),
        E__pixSampleConnections("pix-sample connections", 'J', GLFW_MOD_SHIFT),
        E__qedges("q-edges", 'K'),
        E__qports("q-ports", 'K', GLFW_MOD_SHIFT),
        E__snappingConnections("snapping connections", 'U'),
        E__splineCurves("Spline curves", 'X'),
        E__splinePolygons("Spline polygons", 'X', GLFW_MOD_SHIFT),
        E__streamlines("field streamlines", '6'),

        L__corner("corner", 'A'),
        L__curves("curves", 'C'),
        L__faces("faces", 'B'),
        L__ffAngles("ff angles", '9'),
        L__nTracedBlackPix("ff streamline counts", '='),
        L__pixels("pixels", 'N'),
        L__qexData("ex data", 'J'),
        L__vertex("vertex", 'V'),

        P__BRUSH("BRUSH", 'S'),
        P__isNonManifold("non-manifold?", 'D'),
        P__qvertices("q-vertices", 'H'),
        P__singularities("singularities", '8'),
        P__snappingConstraints("snapping constraints", 'U', GLFW_MOD_SHIFT),
        P__streamlines("field streamlines [P]", '6', GLFW_MOD_SHIFT),
        P__tangentDir("tangent dir", '7'),

        V__ffVector0("crosses u", '5'),
        V__ffVector1("crosses v", '5'),
        V__SobelTangents("Sobel tangents", '3'),
        V__SobelNormals("Sobel normals", ',') {

    custom_window_wide         = 150.f;
    custom_window_wpos         = 0.f;
    custom_window_is_collapsed = false;
    viewer_window_wide         = 150.f;
    viewer_window_wpos         = custom_window_wide;
    viewer_window_is_collapsed = true;

    add_gui_edges(&E__bezierCurves);
    add_gui_edges(&E__bezierPolygons);
    add_gui_edges(&E__blackPixNearestVertex);
    add_gui_edges(&E__customCuts);
    add_gui_edges(&E__lineMapSegments);
    add_gui_edges(&E__narrowBandBoundary);
    add_gui_edges(&E__periodJumps);
    add_gui_edges(&E__pixSampleConnections);
    add_gui_edges(&E__qedges);
    add_gui_edges(&E__qports);
    add_gui_edges(&E__snappingConnections);
    add_gui_edges(&E__splineCurves);
    add_gui_edges(&E__splinePolygons);
    add_gui_edges(&E__streamlines);

    add_gui_fcolors(&F__ffAngles);
    add_gui_fcolors(&F__isCombingOk);
    add_gui_fcolors(&F__narrowBand);
    add_gui_fcolors(&F__uvFlipped);
    selected_fcolor = 0;

    add_gui_labels(&L__corner);
    add_gui_labels(&L__curves);
    add_gui_labels(&L__faces);
    add_gui_labels(&L__ffAngles);
    add_gui_labels(&L__nTracedBlackPix);
    add_gui_labels(&L__pixels);
    add_gui_labels(&L__qexData);
    add_gui_labels(&L__vertex);

    add_gui_points(&P__BRUSH);
    add_gui_points(&P__isNonManifold);
    add_gui_points(&P__qvertices);
    add_gui_points(&P__singularities);
    add_gui_points(&P__snappingConstraints);
    add_gui_points(&P__streamlines);
    add_gui_points(&P__tangentDir);

    add_gui_vectors(&V__ffVector0);
    add_gui_vectors(&V__ffVector1);
    add_gui_vectors(&V__SobelTangents);
    add_gui_vectors(&V__SobelNormals);

    P__BRUSH.clr = CustomColors::Magenta;
    P__BRUSH.pos.resize(0, 3);

    P__tangentDir.show();

    for (int k = 0; k < 3; k++) {
      streamline_color1[k] = CustomColors::Red[k];
      streamline_color2[k] = CustomColors::DeepSkyBlue[k];
    }
  }

  // ===========================================================================

  void IGSVGui::init(igl::opengl::glfw::Viewer* _viewer) {
    BaseGui::init(_viewer);
    viewer->data().show_texture2 = false;
  }

  // ===========================================================================

  bool IGSVGui::key_down(int key, int modifiers) {
    if (ImGui::GetIO().WantCaptureKeyboard)
      return true;

    if (modifiers & GLFW_MOD_SHIFT) {

    } else {
      switch (key) {

        // case GLFW_KEY_Q:
        //   update_frame_field();
        //   return true;
        //
        // case GLFW_KEY_W:
        //   update_parametrization();
        //   return true;
        //
        // case GLFW_KEY_E:
        //   update_extraction();
        //   return true;
        //
        // case GLFW_KEY_R:
        //   serialize_data();
        //   return true;

      default:;
      }
    }

    return BaseGui::key_down(key, modifiers);
  }

  // ===========================================================================

  void IGSVGui::update_narrow_band() {
    if (!data.compute_narrow_band())
      close_window();
    store_entities();
    redraw();
  }

  // ===========================================================================

  void IGSVGui::update_mesh() {
    data.clear_extraction();

    if (!data.init_mesh(/* recompute = */ true))
      close_window();

    if (!data.init_frame_field(/* recompute = */ true))
      close_window();

    data.UV.resize(0, 2);
    data.FUV.resize(0, 3);

    // if (!data.init_parametrization(/* recompute = */ true))
    //   close_window();

    store_entities();
    redraw();
  }

  // ===========================================================================

  void IGSVGui::update_frame_field() {

    if (!data.init_frame_field(/* recompute = */ true))
      close_window();

    store_entities();
    redraw();
  }

  // ===========================================================================

  void IGSVGui::update_parametrization() {

    if (!data.init_parametrization(/* recompute = */ true))
      close_window();

    store_entities();
    redraw();
  }

  // ===========================================================================

  void IGSVGui::update_extraction() {
    data.clear_extraction();

    if (!data.init_extraction())
      close_window();

    store_entities();
    redraw();
  }

  // ===========================================================================

  void IGSVGui::draw_viewer_menu() {
    BaseGui::draw_viewer_menu();
    if (ImGui::Checkbox("show_unused_qedges", &show_unused_qedges) ||                                             //
        ImGui::Checkbox("show_bisector_field", &show_bisector_field) ||                                           //
        ImGui::Combo("##qedge color", &selected_qedge_color, QEdgeColorLabels, IM_ARRAYSIZE(QEdgeColorLabels)) || //
        ImGui::Checkbox("uv mesh?", &render_uv_mesh)) {
      store_entities();
      redraw();
    }
  }

  // ===========================================================================

  void IGSVGui::serialize_data() {
    igl::serialize(data, "state", data.opts.in.filename_input + ".state"); // also createst the dirs if needed
  }

  // ===========================================================================

  bool IGSVGui::mouse_down(int button, int modifier) {
    if (BaseGui::mouse_down(button, modifier))
      return true;
    if (selected_texture == TextureType::I__MASK) {
      if (button == GLFW_MOUSE_BUTTON_RIGHT || button == GLFW_MOUSE_BUTTON_MIDDLE) {
        brush_mode = BrushMode::White;
        return true;
      }
      if (button == GLFW_MOUSE_BUTTON_LEFT) {
        brush_mode = BrushMode::Black;
        return true;
      }
    }
    return false;
  }

  // ===========================================================================

  bool IGSVGui::mouse_up(int button, int modifier) {
    if (brush_mode != BrushMode::Inactive) {
      brush_mode = BrushMode::Inactive;
      return true;
    }
    return BaseGui::mouse_up(button, modifier);
  }

  // ===========================================================================

  bool IGSVGui::mouse_move(int mouse_x, int mouse_y) {
    if (brush_mode != BrushMode::Inactive) {
      // project the origin to get z-value
      const float mouse_z                      = (igl::project(Eigen::Vector3f(0., 0., 0.), //
                                          viewer->core.view,           //
                                          viewer->core.proj,           //
                                          viewer->core.viewport))[2];
      const Eigen::Vector3f point_sketch_space = igl::unproject( //
          Eigen::Vector3f(                                       //
              (float)mouse_x,                                    //
              viewer->core.viewport[3] - (float)mouse_y,         //
              mouse_z),                                          //
          viewer->core.view,                                     //
          viewer->core.proj,                                     //
          viewer->core.viewport                                  //
      );
      P__BRUSH.pos.conservativeResize(P__BRUSH.pos.rows() + 1, 3);
      P__BRUSH.pos.bottomRows<1>() << point_sketch_space.cast<double>().transpose();
      cv::circle(data.detail_mask,                                                // mask
                 cv::Point(point_sketch_space(0),                                 // center.x
                           (float)data.detail_mask.rows - point_sketch_space(1)), // center.y
                 data.sw_avg,                                                     // radius
                 cv::Scalar(brush_mode == BrushMode::White ? 255 : 0),            // color
                 cv::FILLED,                                                      // thickness
                 cv::LINE_AA                                                      // line type
      );
      return redraw();
    }
    return BaseGui::mouse_move(mouse_x, mouse_y);
  }

  // ===========================================================================

  void IGSVGui::draw_custom_window() {

    float menu_width = custom_window_wide * menu_scaling();
    ImGui::SetNextWindowPos(ImVec2(custom_window_wpos, 0.0f), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(ImVec2(menu_width, -1.0f), ImVec2(menu_width, -1.0f));
    ImGui::SetNextWindowCollapsed(custom_window_is_collapsed, ImGuiSetCond_FirstUseEver);

    bool _viewer_menu_visible = true;
    ImGui::Begin("Options", &_viewer_menu_visible,
                 ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_AlwaysAutoResize);
    ImGui::PushItemWidth(0.9f * ImGui::GetWindowWidth());

    if (ImGui::SliderInt("##threshold", (int*)(&data.opts.in.threshold), 10, 250, "in|threshold= %.0f")) {
      if (!data.init_sketch())
        close_window();
      store_entities();
      redraw();
    }

    if (ImGui::SliderFloat("##nb_radius", &data.opts.tri.narrow_band_radius, 0.0, 2.0, "tri|nb_radius= %.2f"))
      update_narrow_band();

    ImGui::SliderFloat("##tanratiothresh", &data.opts.ff.tangent_ratio_threshold, 0.0f, 0.5f, "ff|labels_thresh= %.4f");
    ImGui::SliderFloat("##tanratiolen", &data.opts.ff.tangent_ratio_streamlen, 0.0f, 100.0f, "ff|labels_length= %.1f");
    ImGui::SliderFloat("##wsnap", &data.opts.uv.wsnap, 0.0f, 1000.0f, "uv|w_snap= %.4f", 3.0f);
    ImGui::SliderFloat("##mask_factor", &data.opts.in.mask_factor, 0.0f, 10.0f, "uv|mask_factor= %.1f", 2.0f);
    ImGui::SliderFloat("##scale", &data.opts.uv.scale_multiplier, 0.1, 25.0, "uv|scale= %.2f");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (ImGui::Button("Reset Options", ImVec2(-1, 0))) {
      data.clear();
      if (!data.init_sketch())
        close_window();
      update_mesh();
      data.opts.set_default();
    }

    if (ImGui::Button("Reset Mask", ImVec2(-1, 0))) {
      data.detail_mask.setTo(128);
      redraw();
    }

    if (ImGui::Button("Select File", ImVec2(-1, 0))) {
      std::string filename = igl::file_dialog_open();
      if (!filename.empty()) {
        boost::filesystem::path p(filename);
        // data.opts.in.dataname = p.stem().string();
        // data.opts.in.subdir   = p.parent_path().stem().string();
        // data.opts.in.fulldir  = data.opts.in.folder + "/" + data.opts.in.subdir + "/" + data.opts.in.dataname;
        data.clear();
        data.init_sketch();
        do_align_camera_center = true;
        store_entities();
        redraw();
      }
    }

    if (ImGui::Button("Mesh", ImVec2(-1, 0))) {
      data.opts.tri.estimate_max_area = true;
      update_mesh();
    }

    if (ImGui::Button("Trace", ImVec2(-1, 0))) {
      if (!data.trace_streamlines())
        close_window();
      store_entities();
      redraw();
    }

    if (ImGui::Button("Orientation field", ImVec2(-1, 0)))
      update_frame_field();

    if (ImGui::Button("Parametrize", ImVec2(-1, 0)))
      update_parametrization();

    if (ImGui::Button("Extract", ImVec2(-1, 0)))
      update_extraction();

    if (ImGui::Button("Serialize", ImVec2(-1, 0)))
      serialize_data();

    // if (ImGui::Button("Export info", ImVec2(-1, 0)))
    //   data.export_info();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ImGui::PopItemWidth();
    ImGui::End();
  }

  // ===========================================================================

}; // namespace IGSV
