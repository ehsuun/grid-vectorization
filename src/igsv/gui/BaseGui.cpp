#include <igsv/gui/BaseGui.h>

#include <igsv/common/CustomColors.h>
#include <igsv/common/VectorizationData.h>

#include <GLFW/glfw3.h>
#include <igl/parula.h>
#include <imgui/imgui.h>

// texture matrix
using TxMat = Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic>;

namespace IGSV {

  std::string GuiEntity::key_name() const {
    std::stringstream ss;
    if (key == -1) {
      ss << "     " << name;
    } else {
      ss << "[";
      if (modifier & GLFW_MOD_SHIFT) {
        ss << "s";
      } else if (modifier & GLFW_MOD_ALT) {
        ss << "a";
      } else if (modifier & GLFW_MOD_CONTROL) {
        ss << "c";
      } else {
        ss << " ";
      }
      ss << char(key) << "] " << name;
    }
    return ss.str();
  }

  void BaseGui::close_window() { glfwSetWindowShouldClose(viewer->window, GLFW_TRUE); }

  BaseGui::BaseGui(VectorizationData& data) : data(data) {
    F_white.clr = CustomColors::White;
    add_gui_fcolors(&F_white);
    selected_fcolor = 0;
  }

  void BaseGui::init(igl::opengl::glfw::Viewer* _viewer) {
    igl::opengl::glfw::imgui::ImGuiMenu::init(_viewer);
    viewer->data().line_width    = 1.0f;
    viewer->data().point_size    = 10.0f;
    viewer->data().show_lines    = false; // wireframe
    viewer->data().show_overlay  = true;
    viewer->data().show_texture  = true;
    viewer->data().show_texture2 = true;
    viewer->core.background_color << 0.8f, 0.8f, 0.8f, 1.0f;
    store_entities();
    if (!redraw())
      close_window();
  }

  void BaseGui::draw_labels(const igl::opengl::ViewerData& vdata) {
    ImGuiMenu::draw_labels(vdata);
    for (const auto l : gui_labels) {
      if (!l->is_shown)
        continue;
      const int nc = l->clr.rows();
      const int np = l->pos.rows();
      const int ns = l->str.size();
      if (np != ns) {
        DEBUG_PRINT_WARNING(" in BaseGui::draw_labels() : pos.rows()=%d does not match str.size()=%d", np, ns);
      } else if (nc != 1 && nc != ns) {
        DEBUG_PRINT_WARNING(" in BaseGui::draw_labels() : clr.rows()=%d does not match str.size()=%d", nc, ns);
      } else {
        for (int i = 0; i < ns; ++i) {
          draw_text(l->pos.row(i).transpose(), Eigen::Vector3d(0, 0, 1e-4), l->str[i],
                    nc == 1 ? l->clr : l->clr.row(i));
        }
        if (l->show_points) {
          viewer->data().add_points(l->pos, l->clr);
        }
      }
    }
  }

  Eigen::VectorXd band_to_vertex_weights(const Eigen::VectorXd& band_w,  //
                                         const Eigen::VectorXi& nearest, //
                                         const int nv,                   //
                                         bool normalize = true) {
    Eigen::VectorXd w(nv);
    w.fill(0.0);
    for (int k = 0; k < nearest.rows(); ++k) {
      w(nearest(k)) = band_w(k);
    }
    if (normalize) {
      if (band_w.maxCoeff() > 1e-4) {
        w.array() /= band_w.maxCoeff();
      }
    }
    return w;
  }

  void BaseGui::set_texture_image() {
    viewer->data().show_texture = true;
    switch (selected_texture) {

    case TextureType::NONE: {
      TxMat White;
      White.resize(1, 1);
      White.fill(255);
      viewer->data().set_texture(White, White, White);
      viewer->data().set_uv(data.UV_generated, data.F);
    } break;

    case TextureType::I__COLOR:
      set_texture_rgb(data.input);
      break;

    case TextureType::I__MASK:
      set_texture_grey(data.detail_mask);
      break;

    case TextureType::I__BINARY:
      set_texture_grey(data.bw);
      break;

    case TextureType::I__GREY_BLURRED:
      set_texture_grey(data.grey_blurred);
      break;

    case TextureType::V__SNAPPING_WEIGHTS:
      set_texture_distance(data.V_snappingWeights);
      break;

    case TextureType::I__INV_DIST_TRANSFORM_01:
      set_texture_grey(data.idt_01, true);
      break;

    case TextureType::I__SW:
      set_texture_grey(data.sw, true);
      break;

    default:
      DEBUG_PRINT_WARNING("invalid texture type = %d", selected_texture);
    }
  }

  bool BaseGui::redraw() {
    // save and load settings
    const auto lw  = viewer->data().line_width;
    const auto tx  = viewer->data().show_texture;
    const auto tx2 = viewer->data().show_texture2;
    const auto ln  = viewer->data().show_lines;
    viewer->data().clear();
    viewer->data().line_width   = lw;
    viewer->data().show_texture = tx;
    viewer->data().show_lines   = ln;
    viewer->core.set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_NO_ROTATION);
    // set mesh and align if needed
    if (data.UV.rows() > 0)
      viewer->data().set_uv2(data.UV, data.FUV);
    if (render_uv_mesh)
      viewer->data().set_mesh(data.UV, data.FUV);
    else
      viewer->data().set_mesh(data.V, data.F);
    viewer->data().set_colors(gui_fcolors[selected_fcolor % gui_fcolors.size()]->clr);
    // add curves
    for (auto const c : gui_curves) {
      if (c->is_shown && c->pos.rows() > 1)
        viewer->data().add_edges(c->pos.topRows(c->pos.rows() - 1),    //
                                 c->pos.bottomRows(c->pos.rows() - 1), //
                                 c->clr.topRows(c->pos.rows() - 1));
    }
    // add edges
    for (auto const e : gui_edges) {
      if (e->is_shown && e->pos0.rows() > 0)
        viewer->data().add_edges(e->pos0, e->pos1, e->clr);
    }
    // add points
    for (auto const p : gui_points) {
      if (p->is_shown && p->pos.rows() > 0)
        viewer->data().add_points(p->pos, p->clr);
    }
    // add vectors
    for (auto const v : gui_vectors) {
      if (v->is_shown && v->pos.rows() > 0)
        viewer->data().add_edges(v->pos, v->pos + vector_scale * v->vec, v->clr);
    }
    if (do_align_camera_center) {
      if (viewer->data().V.rows() > 0) {
        viewer->core.align_camera_center(viewer->data().V, viewer->data().F);
        viewer->core.camera_base_zoom *= 3.0f;
        set_texture_grid();
        do_align_camera_center = false;
      }
    }
    // set texture
    set_texture_image();
    return true;
  }

  bool BaseGui::set_texture_rgb(const cv::Mat& I) {
    const int h = I.rows;
    const int w = I.cols;
    if (h == 0 || w == 0) {
      DEBUG_PRINT_WARNING("cannot load rgb texture from empty image");
      return false;
    }
    TxMat R, G, B;
    R.setZero(w, h);
    G.setZero(w, h);
    B.setZero(w, h);
    for (int y = 0; y < h; y++) {
      for (int x = 0; x < w; x++) {
        R(x, y) = I.at<cv::Vec3b>(h - 1 - y, x)[2];
        G(x, y) = I.at<cv::Vec3b>(h - 1 - y, x)[1];
        B(x, y) = I.at<cv::Vec3b>(h - 1 - y, x)[0];
      }
    }
    viewer->data().set_texture(R, G, B);
    viewer->data().set_uv(data.UV_generated, data.F);
    return true;
  }

  bool BaseGui::set_texture_grey(const cv::Mat& I, bool normalize,   //
                                 unsigned char lo, unsigned char hi, //
                                 unsigned char clamp_lo, unsigned char clamp_hi) {
    const int h = I.rows;
    const int w = I.cols;
    if (h == 0 || w == 0) {
      DEBUG_PRINT_WARNING("cannot load grey texture from empty image");
      return false;
    }
    cv::Mat I_clone = I.clone();
    if (normalize)
      cv::normalize(I_clone, I_clone, lo, hi, cv::NORM_MINMAX, CV_8UC1);
    I_clone.setTo(clamp_lo, I_clone < clamp_lo);
    I_clone.setTo(clamp_hi, I_clone > clamp_hi);
    // if (normalize) cv::normalize(I_clone, I_clone, lo, hi, cv::NORM_MINMAX, CV_8UC1);
    TxMat J;
    J.setZero(w, h);
    for (int y = 0; y < h; y++)
      for (int x = 0; x < w; x++)
        J(x, y) = I_clone.at<unsigned char>(h - 1 - y, x);
    viewer->data().set_texture(J, J, J);
    viewer->data().set_uv(data.UV_generated, data.F);
    return true;
  }

  bool BaseGui::set_texture_distance(const Eigen::VectorXd& w) {
    if (w.rows() == 0) {
      DEBUG_PRINT_WARNING("cannot show empty distance texture");
      return false;
    }
    TxMat R, G, B;
    R.setConstant(1, 256, 255);
    G.setConstant(1, 256, 255);
    B.setConstant(1, 256, 255);
    Eigen::RowVector3d color;
    for (int i = 0; i < 256; ++i) {
      igl::parula(1. / 255. * (double)(i), color[0], color[1], color[2]);
      R(0, i) = std::round(color[0] * 255);
      G(0, i) = std::round(color[1] * 255);
      B(0, i) = std::round(color[2] * 255);
    }
    // generate texture coords
    Eigen::MatrixXd UV_1D = Eigen::MatrixXd::Zero(data.V.rows(), 2);
    UV_1D.col(1)          = 0.98 * w.array() + 0.01;

    // set texture and coords
    viewer->data().set_texture(R, G, B);
    viewer->data().set_uv(UV_1D, data.F);
    return true;
  }

  bool BaseGui::set_texture_grid() {
    const auto colorU = CustomColors::U();
    const auto colorV = CustomColors::V();

    const unsigned char ru = std::round(colorU(0) * 255), rv = std::round(colorV(0) * 255),
                        gu = std::round(colorU(1) * 255), gv = std::round(colorV(1) * 255),
                        bu = std::round(colorU(2) * 255), bv = std::round(colorV(2) * 255);
    // generate texture
    const int gtex_size   = 1024;
    const int strip_width = std::pow(2, grid_texture_log_width);
    TxMat R, G, B, A;
    R.setConstant(gtex_size, gtex_size, 255);
    G.setConstant(gtex_size, gtex_size, 255);
    B.setConstant(gtex_size, gtex_size, 255);
    A.setConstant(gtex_size, gtex_size, 0);
    // left
    R.leftCols(strip_width).fill(ru);
    G.leftCols(strip_width).fill(gu);
    B.leftCols(strip_width).fill(bu);
    A.leftCols(strip_width).fill(255);
    // right
    R.rightCols(strip_width).fill(ru);
    G.rightCols(strip_width).fill(gu);
    B.rightCols(strip_width).fill(bu);
    A.rightCols(strip_width).fill(255);
    // top
    R.topRows(strip_width).fill(rv);
    G.topRows(strip_width).fill(gv);
    B.topRows(strip_width).fill(bv);
    A.topRows(strip_width).fill(255);
    // bottom
    R.bottomRows(strip_width).fill(rv);
    G.bottomRows(strip_width).fill(gv);
    B.bottomRows(strip_width).fill(bv);
    A.bottomRows(strip_width).fill(255);
    viewer->data().set_texture2(R, G, B, A);
    return true;
  }

  bool BaseGui::key_down(int key, int modifiers) {

    if (modifiers & GLFW_MOD_SHIFT) {

      switch (key) {

      case GLFW_KEY_1:
        selected_texture += (selected_texture == 0) ? NUMBER_OF_TEXTURES - 1 : -1;
        return redraw();

      case GLFW_KEY_Q:
        selected_fcolor--;
        if (selected_fcolor < 0)
          selected_fcolor += gui_fcolors.size();
        return redraw();

      case GLFW_KEY_W:
        selected_fcolor++;
        selected_fcolor %= gui_fcolors.size();
        return redraw();

      case GLFW_KEY_O:
        viewer->data().point_size -= (viewer->data().point_size > 3.0f) ? 2.0f : 0.0f;
        return true; // decrease point size

      case GLFW_KEY_P:
        viewer->data().point_size += 2.0f;
        return true; // increase point size
      }

    } else {

      switch (key) {

      case GLFW_KEY_GRAVE_ACCENT:
        viewer->data().show_overlay = !viewer->data().show_overlay;
        return true;

      case GLFW_KEY_1:
        selected_texture = (selected_texture == NUMBER_OF_TEXTURES - 1) ? 0 : selected_texture + 1;
        return redraw();

      case GLFW_KEY_2:
        viewer->data().show_texture2 = !viewer->data().show_texture2;
        if (viewer->data().show_texture2) {
          set_texture_grid();
        }
        return true;

      case GLFW_KEY_O:
        viewer->data().line_width *= 0.5f;
        return true; // decrease line width

      case GLFW_KEY_P:
        viewer->data().line_width *= 2.0f;
        return true; // increase line width
        // default: return igl::opengl::glfw::imgui::ImGuiMenu::key_down(key, modifiers);
      }
    }

    for (auto const& x : gui_points) {
      if (x->key == key && ((x->modifier == 0 && modifiers == 0) || (x->modifier & modifiers)))
        x->toggle();
    }

    for (auto const& x : gui_vectors) {
      if (x->key == key && ((x->modifier == 0 && modifiers == 0) || (x->modifier & modifiers)))
        x->toggle();
    }

    for (auto const& x : gui_edges) {
      if (x->key == key && ((x->modifier == 0 && modifiers == 0) || (x->modifier & modifiers)))
        x->toggle();
    }

    for (auto const& x : gui_labels) {
      if (x->key == key && ((x->modifier == 0 && modifiers == 0) || (x->modifier & modifiers)))
        x->toggle();
    }

    return redraw();
  }

  void BaseGui::draw_viewer_window() {
    float menu_width = viewer_window_wide * menu_scaling();
    ImGui::SetNextWindowPos(ImVec2(viewer_window_wpos, 0.0f), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiSetCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(ImVec2(menu_width, -1.0f), ImVec2(menu_width, -1.0f));
    ImGui::SetNextWindowCollapsed(viewer_window_is_collapsed, ImGuiSetCond_FirstUseEver);

    bool _viewer_menu_visible = true;
    ImGui::Begin("Viewer", &_viewer_menu_visible, ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_AlwaysAutoResize);
    ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.8f);
    if (callback_draw_viewer_menu) {
      callback_draw_viewer_menu();
    } else {
      draw_viewer_menu();
    }
    ImGui::PopItemWidth();
    ImGui::End();
  }

  void BaseGui::draw_viewer_menu() {
    // ImGui::ShowDemoWindow();
    if (ImGui::CollapsingHeader("General", ImGuiTreeNodeFlags_DefaultOpen)) {


      ImGui::InputFloat("#cam:zoom", &(viewer->core.camera_zoom), 0.001f, 100.0f, "%.3f");
      ImGui::InputFloat("#cam:trans:x", &(viewer->core.camera_translation[0]), -1000.f, 1000.0f, "%.1f");
      ImGui::InputFloat("#cam:trans:y", &(viewer->core.camera_translation[1]), -1000.f, 1000.0f, "%.1f");

      // ImGui::SliderFloat("#cam:zoom", &(viewer->core.camera_zoom), 0.0001f, 100.0f, "zoom = %0.4f");
      // ImGui::SliderFloat("##cam:trans:x", &(viewer->core.camera_translation[0]), -1000.0f, 1000.0f, "tx = %0.2f");
      // ImGui::SliderFloat("##cam:trans:y", &(viewer->core.camera_translation[1]), -1000.0f, 1000.0f, "ty = %0.2f");
      // ImGui::DragFloat("cam bzoom", &(viewer->core.camera_base_zoom), 0.05f, 0.1f, 20.0f);
      // ImGui::DragFloat("cam btrans x", &(viewer->core.camera_base_translation[0]), 0.0f, 1000.0f, 20.0f);
      // ImGui::DragFloat("cam btrans y", &(viewer->core.camera_base_translation[1]), 0.0f, 1000.0f, 20.0f);
      // ImGui::DragFloat("cam btrans z", &(viewer->core.camera_base_translation[2]), 0.0f, 1000.0f, 20.0f);
      // ImGui::DragFloat("cam trans z", &(viewer->core.camera_translation[2]), 0.0f, 1000.0f, 20.0f);

      ImGui::Checkbox("show overlay", &viewer->data().show_overlay);
      if (ImGui::Checkbox("texture : grid", &viewer->data().show_texture2)) {
        if (viewer->data().show_texture2) {
          set_texture_grid();
        }
      }
      if (ImGui::Combo("##texture type", &selected_texture, TextureLabels, IM_ARRAYSIZE(TextureLabels))) {
        redraw();
      }
      ImGui::SliderFloat("##line width", &viewer->data().line_width, 0.01f, 100.0f, "line w = %.2f", 2.0f);
      ImGui::SliderFloat("##point size", &viewer->data().point_size, 0.01f, 50.0f, "pt size = %.1f");
      if (ImGui::SliderFloat("##vector scale", &vector_scale, 0.01f, 10.0f, "vec scale = %.1f")) {
        redraw();
      }
      if (ImGui::SliderInt("##grid density", &grid_texture_log_width, 0, 8, "grid = %.0f")) {
        set_texture_grid();
      }
    }
    if (ImGui::CollapsingHeader("Points")) {
      for (auto const& x : gui_points) {
        if (ImGui::Checkbox(x->key_name().c_str(), &(x->is_shown)))
          redraw();
      }
    }
    if (ImGui::CollapsingHeader("Vectors")) {
      for (auto const& x : gui_vectors) {
        if (ImGui::Checkbox(x->key_name().c_str(), &(x->is_shown)))
          redraw();
      }
    }
    if (ImGui::CollapsingHeader("Edges", ImGuiTreeNodeFlags_DefaultOpen)) {
      for (auto const& x : gui_edges) {
        if (ImGui::Checkbox(x->key_name().c_str(), &(x->is_shown)))
          redraw();
      }
    }
    if (ImGui::CollapsingHeader("Labels")) {
      for (auto const& x : gui_labels) {
        if (ImGui::Checkbox(x->key_name().c_str(), &(x->is_shown)))
          redraw();
      }
    }
    if (ImGui::InputInt("selected face", &selected_face)) {
      store_entities();
      redraw();
    }
  }

  /*static*/ void BaseGui::show_help_marker(const char* desc) // from imgui_demo.cpp
  {
    ImGui::SameLine(); // Always on the same line!
    ImGui::TextDisabled("(?)");
    if (ImGui::IsItemHovered()) {
      ImGui::BeginTooltip();
      ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
      ImGui::TextUnformatted(desc);
      ImGui::PopTextWrapPos();
      ImGui::EndTooltip();
    }
  }

} // namespace IGSV
