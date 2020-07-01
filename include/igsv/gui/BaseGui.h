//=============================================================================
//
//  CLASS : BaseGui
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include "Entities.h"

#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <opencv2/opencv.hpp>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== FORWARD DECLARATIONS ===================================================

  class VectorizationData;

  //== CLASS DEFINITION =======================================================

  class BaseGui : public igl::opengl::glfw::imgui::ImGuiMenu {

  public:
    BaseGui(VectorizationData&);
    virtual void draw_labels(const igl::opengl::ViewerData&) override;
    virtual void draw_viewer_window() override;
    virtual void draw_viewer_menu() override;
    virtual void init(igl::opengl::glfw::Viewer* _viewer) override;
    virtual bool key_down(int key, int modifiers) override;
    virtual bool redraw();
    virtual void store_entities() {}

    void set_render_uv_mesh(bool val) { render_uv_mesh = val; }

  protected:
    VectorizationData& data;

    void add_gui_curve(GuiCurve* c) { gui_curves.push_back(c); }
    void add_gui_edges(GuiEdges* e) { gui_edges.push_back(e); }
    void add_gui_labels(GuiLabels* l) { gui_labels.push_back(l); }
    void add_gui_points(GuiPoints* p) { gui_points.push_back(p); }
    void add_gui_vectors(GuiVectors* v) { gui_vectors.push_back(v); }
    void add_gui_fcolors(GuiFColors* c) { gui_fcolors.push_back(c); }
    void close_window();

    static void show_help_marker(const char* desc);

    //// gui settings
    float custom_window_wide        = 0.f;
    float custom_window_wpos        = 0.f;
    bool custom_window_is_collapsed = false;
    float viewer_window_wide        = 150.f;
    float viewer_window_wpos        = 0.f;
    bool viewer_window_is_collapsed = false;

    //// render settings
    bool do_align_camera_center = true;
    bool render_uv_mesh         = false;
    int grid_texture_log_width  = 5;
    int selected_face           = -1;
    int selected_fcolor         = 0;
    int selected_texture        = TextureType::I__COLOR; // '~' '1'
    float vector_scale          = 1.f;

    GuiFColors F_white; // default face colors (all white)

    enum TextureType {
      NONE = 0,
      I__COLOR,
      I__MASK,
      I__BINARY,
      I__GREY_BLURRED,
      V__SNAPPING_WEIGHTS,
      I__INV_DIST_TRANSFORM_01,
      I__SW,
      NUMBER_OF_TEXTURES
    };

  private:
    //// texture labels
    const char* TextureLabels[NUMBER_OF_TEXTURES] = {
      "no texture",           //
      "I: input",             //
      "I: mask",              //
      "I: binarized",         //
      "I: grey blurred",      //
      "V: snapping w",        //
      "I: inv DT 01",         //
      "I: SW",                //
    };

    //// data
    std::vector<GuiCurve*> gui_curves;
    std::vector<GuiEdges*> gui_edges;
    std::vector<GuiLabels*> gui_labels;
    std::vector<GuiPoints*> gui_points;
    std::vector<GuiVectors*> gui_vectors;
    std::vector<GuiFColors*> gui_fcolors;

    //// textures
    bool set_texture_rgb(const cv::Mat& I);
    bool set_texture_grey(const cv::Mat& I, bool normalize = false, unsigned char lo = 0, unsigned char hi = 255,
                          unsigned char clamp_lo = 0, unsigned char clamp_hi = 255);
    bool set_texture_distance(const Eigen::VectorXd&);
    void set_texture_image(); // wrapper

    bool set_texture_grid();

  }; // class Gui
};   // namespace IGSV
