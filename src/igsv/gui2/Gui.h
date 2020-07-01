//
// Created by Tibor Stanko on 15/04/2020.
//

#ifndef IGSV_GUI2_H
#define IGSV_GUI2_H

#include <igl/opengl/gl.h>

#include "drawable.h"

#include <igsv/extraction_entities/Bezier.h>
#include <igsv/extraction_entities/Chain.h>
#include <igsv/extraction_entities/QEdge.h>
#include <igsv/extraction_entities/QVertex.h>

#include <Eigen/Core>
#include <GLFW/glfw3.h>
#include <map>
#include <opencv2/opencv.hpp>

namespace IGSV {

  using VectorXb = Eigen::Matrix<bool, Eigen::Dynamic, 1>;

  class Gui2 {

  public:
    Gui2(cv::Mat& input, //
         cv::Mat& mask,  //
         // cv::Mat& grid,                            //
         float& narrow_band_radius, //
         float& scale_multiplier,   //
         float& mask_factor,        //
         Eigen::MatrixXd& V,        //
         Eigen::MatrixXd& UV,       //
         Eigen::MatrixXd& BC,       //
         Eigen::MatrixXd& X0,       //
         Eigen::MatrixXd& X1,       //
         Eigen::MatrixXi& F,        //
         Eigen::MatrixXi& FUV,      //
         VectorXb& F_isNarrow,      //
         std::vector<QVertex>& QV,  //
         std::vector<QEdge>& QE,    //
         std::vector<Chain>& CH     //
         )                          //
        : _input(input),            //
          _mask(mask),              //
          // _grid(grid),                             //
          _narrow_band_radius(narrow_band_radius), //
          _scale_multiplier(scale_multiplier),     //
          _mask_factor(mask_factor),               //
          _V(V),                                   //
          _UV(UV),                                 //
          _BC(BC),                                 //
          _X0(X0),                                 //
          _X1(X1),                                 //
          _F(F),                                   //
          _FUV(FUV),                               //
          _F_isNarrow(F_isNarrow),                 //
          _QV(QV),                                 //
          _QE(QE),                                 //
          _CH(CH)                                  //
    {
      initialize();
    }

    void add_points(const std::string& _name, const IGSV::drawable::points& _points);
    void add_edges(const std::string& _name, const IGSV::drawable::edges& _edges);
    void add_triangles(const std::string& _name, const IGSV::drawable::triangles& _triangles);

    // void add_texture(const char* _name, cv::Mat* _image, bool _do_normalize = false);

    void update_structures();

    void render();

    std::function<bool()> _callback;
    bool _callback_success = true;

  private:
    void initialize();

    void init_grid_texture();

    bool init_window();
    void init_textures();
    void init_shader_uniforms();
    void init_vertex_arrays();

    void fill_vertex_arrays__triangles();
    void fill_vertex_arrays__edges();
    void fill_vertex_arrays__points();

    void mouse_events();
    void keyboard_events();
    void recompute_matrices();
    void render_custom_menu();

    void unproject(float x_win, float y_win, float& x_obj, float& y_obj);

  private:
    cv::Mat& _input;
    cv::Mat& _mask;
    cv::Mat _grid;

    float& _narrow_band_radius;
    float& _scale_multiplier;
    float& _mask_factor;

    Eigen::MatrixXd& _V;
    Eigen::MatrixXd& _UV;
    Eigen::MatrixXd& _BC;
    Eigen::MatrixXd& _X0;
    Eigen::MatrixXd& _X1;
    Eigen::MatrixXi& _F;
    Eigen::MatrixXi& _FUV;
    VectorXb& _F_isNarrow;
    std::vector<QVertex>& _QV;
    std::vector<QEdge>& _QE;
    std::vector<Chain>& _CH;

    // std::vector<cv::Mat*> _textures;
    // std::vector<const char*> _texture_labels;
    // std::vector<bool> _texture_normalize;

    std::map<std::string, IGSV::drawable::points> _render_points;
    std::map<std::string, IGSV::drawable::edges> _render_edges;
    std::map<std::string, IGSV::drawable::triangles> _render_triangles;

    // ## textures (image, mask, grid, ...)
    GLuint textureID[3]{};

    // ## shader: triangles
    int n_triangles = 0;
    GLuint programID{};
    // narrow band mesh
    GLuint vao{};
    GLuint vbo_positions{};
    GLuint vbo_colors{};
    GLuint vbo_uvcoords{};
    GLuint vbo_triangles{};
    // base mesh
    GLuint vao_base{};
    GLuint vbo_base_positions{};
    GLuint vbo_base_colors{};
    GLuint vbo_base_uvcoords{};
    GLuint vbo_base_triangles{};
    // locations
    GLint locModel{};
    GLint locView{};
    GLint locProj{};
    GLint locImageTexture{};
    GLint locMaskTexture{};
    GLint locGridTexture{};
    GLint locEnableImageTexture{};
    GLint locEnableMaskTexture{};
    GLint locEnableGridTexture{};

    // ## shader: lines
    int n_edges = 0;
    GLuint programID_lines{};
    GLuint vao_lines{};
    GLuint vbo_lines_positions{};
    GLuint vbo_lines_colors{};
    GLuint vbo_lines_edges{};
    GLint locModel_lines{};
    GLint locView_lines{};
    GLint locProj_lines{};
    GLint locLineThickness_lines{};

    // ## shader : points
    int n_points = 0;
    GLuint programID_points{};
    GLuint vao_points{};
    GLuint vbo_points_positions{};
    GLuint vbo_points_colors{};
    GLint locModel_points{};
    GLint locView_points{};
    GLint locProj_points{};
    GLint locPointSize_points{};

    // ## glfw window
    GLFWwindow* window{};
    int window_w = 1280;
    int window_h = 720;

    const Eigen::Vector4f bgColor = Eigen::Vector4f(0.2f, 0.2f, 0.2f, 1.0f);

    // ## MVP matrices
    Eigen::Matrix4f model, view, proj;
    const float nearVal = 1.0;
    const float farVal  = 100.0;

    // ## model translation
    float tx = 0.0f;
    float ty = 0.0f;

    // ## model zoom
    float zoom = 0.01f;

    // ## brush stroke
    int nBrushStrokePoints = 0;
    float x0{}, y0{}; // world coords of last two
    float x1{}, y1{}; // brush stroke points

    // ## gui options
    int brushSize           = 5;
    float lineThickness     = 1.0f;
    float pointSize         = 10.0f;
    bool enableImageTexture = true;
    bool enableMaskEditing  = false;
    bool isMaskBeingEdited  = false;
    bool showWireframe      = false;
    int selectedTexture     = 0;
    int gridThickness       = 5;
  };
} // namespace IGSV

#endif // IGSV_GUI2_H
