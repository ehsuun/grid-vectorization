//
// Created by Tibor Stanko on 15/04/2020.
//

#include "Gui.h"
#include "colors.h"

#include <Eigen/Dense>

#include <imgui/imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace IGSV {

  void Gui2::add_points(const std::string& _name, const IGSV::drawable::points& _points) { //
    bool isEnabled = true;
    if (_render_points.find(_name) != _render_points.end())
      isEnabled = _render_points[_name].enabled;
    else if (_name == "verts")
      isEnabled = false; // verts are hidden by default
    _render_points[_name]         = _points;
    _render_points[_name].enabled = isEnabled;
  }

  void Gui2::add_edges(const std::string& _name, const IGSV::drawable::edges& _edges) { //
    bool isEnabled = true;
    if (_render_edges.find(_name) != _render_edges.end())
      isEnabled = _render_edges[_name].enabled;
    else if (_name == "field") // frame field is hidden by default
      isEnabled = false;
    else if (_name == "edges") // edges is hidden by default
      isEnabled = false;

    _render_edges[_name]         = _edges;
    _render_edges[_name].enabled = isEnabled;
  }

  void Gui2::add_triangles(const std::string& _name, const IGSV::drawable::triangles& _triangles) { //
    bool isEnabled = true;
    if (_render_triangles.find(_name) != _render_triangles.end())
      isEnabled = _render_triangles[_name].enabled;
    _render_triangles[_name]         = _triangles;
    _render_triangles[_name].enabled = isEnabled;
  }

  // void Gui2::add_texture(const char* _name, cv::Mat* _image, bool _do_normalize) {
  //   _textures.push_back(_image);
  //   _texture_labels.push_back(_name);
  //   _texture_normalize.push_back(_do_normalize);
  // }

  //=============================================================================

  bool set_texture_from_image(GLuint texID, const cv::Mat& I, bool normalize = false) {

    const int w = I.cols;
    const int h = I.rows;

    cv::Mat temp;
    temp = I.clone();
    if (normalize)
      cv::normalize(temp, temp, 0, 255, cv::NORM_MINMAX, CV_8UC1);
    temp.convertTo(temp, CV_8U);

    std::vector<unsigned char> data(w * h * 4);
    int i = 0;
    unsigned char val;

    switch (temp.channels()) {

    case 1:
      for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++) {
          val = temp.at<unsigned char>(h - 1 - y, x);

          data[i * 4 + 0] = val; // R
          data[i * 4 + 1] = val; // G
          data[i * 4 + 2] = val; // B
          data[i * 4 + 3] = 255;
          i++;
        }
      break;

    case 3:
      for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++) {
          data[i * 4 + 0] = temp.at<cv::Vec3b>(h - 1 - y, x)[2]; // R
          data[i * 4 + 1] = temp.at<cv::Vec3b>(h - 1 - y, x)[1]; // G
          data[i * 4 + 2] = temp.at<cv::Vec3b>(h - 1 - y, x)[0]; // B
          data[i * 4 + 3] = 255;
          i++;
        }
      break;

    case 4:
      for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++) {
          data[i * 4 + 0] = temp.at<cv::Vec4b>(h - 1 - y, x)[2]; // R
          data[i * 4 + 1] = temp.at<cv::Vec4b>(h - 1 - y, x)[1]; // G
          data[i * 4 + 2] = temp.at<cv::Vec4b>(h - 1 - y, x)[0]; // B
          data[i * 4 + 3] = temp.at<cv::Vec4b>(h - 1 - y, x)[3]; // A
          i++;
        }
      break;

    default:
      std::cout << "warning: unexpected # of channels (" << temp.channels() << ") in texture" << std::endl;
      return false;
    }

    glBindTexture(GL_TEXTURE_2D, texID);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, &data[0]);
    return true;
  }

  //=============================================================================

  void Gui2::unproject(float x_win, float y_win, float& x_obj, float& y_obj) {
    Eigen::Vector4f tmp(2.0f * x_win / (float)this->window_w - 1.0f, //
                        2.0f * y_win / (float)this->window_h - 1.0f, //
                        0.0f,                                        //
                        1.0f);

    // flip y coord
    tmp(1) = -tmp(1);
    // apply inverse of MVP
    tmp = (model * view * proj).transpose().inverse() * tmp;
    // convert from homogeneous coords
    x_obj = tmp(0) / tmp(3);
    y_obj = tmp(1) / tmp(3);
    // flip y coord back
    y_obj = (float)this->_input.rows - y_obj;
  }

  //=============================================================================

  void Gui2::render_custom_menu() {
    const ImVec4 color      = ImVec4(0.1f, 0.1f, 0.1f, 1.0f);
    const ImU32 cursorColor = ImColor(color);

    ImGui::Begin("IGSV controls");

    // if (ImGui::Combo("##texture type", &selectedTexture, _texture_labels.data(), _texture_labels.size()))
    //   set_texture_from_image(this->textureID[0], *_textures[selectedTexture], _texture_normalize[selectedTexture]);

    if (ImGui::SliderInt("##gridThickness", &this->gridThickness, 1, 10, "grid thickness: %d")) {
      init_grid_texture();
    }

    ImGui::SliderFloat("##lineThickness", &this->lineThickness, 0.0f, 10.0f, "line thickness: %0.2f");
    ImGui::SliderFloat("##pointSize", &this->pointSize, 0.0f, 100.0f, "point size: %0.2f");
    ImGui::SliderInt("##brush_size", &this->brushSize, 1, 50, "brush size: %dpx");
    ImGui::Checkbox("show image", &this->enableImageTexture);
    ImGui::Checkbox("show wireframe", &this->showWireframe);
    ImGui::Checkbox("edit mask", &this->enableMaskEditing);

    ImGui::SliderFloat("##_narrow_band_radius", &this->_narrow_band_radius, 0.01f, 10.0f, "band = %0.2f");
    ImGui::SliderFloat("##_scale_multiplier", &this->_scale_multiplier, 1.0f, 10.0f, "scale = %0.1f");
    // ImGui::SliderFloat("##_mask_factor", &this->_mask_factor, 0.25f, 4.0f, "maskf = %0.2f");

    if (ImGui::Button("recompute parametrization")) {
      if (_callback) {
        _callback_success = _callback();
        if (_callback_success) {
          this->update_structures();
          this->fill_vertex_arrays__edges();
          this->fill_vertex_arrays__points();
          this->fill_vertex_arrays__triangles();
        }
      }
    }

    for (auto& pair : _render_triangles)
      if (ImGui::Checkbox((pair.first).c_str(), &pair.second.enabled))
        this->fill_vertex_arrays__triangles();

    for (auto& pair : _render_edges)
      if (ImGui::Checkbox((pair.first).c_str(), &pair.second.enabled))
        this->fill_vertex_arrays__edges();

    for (auto& pair : _render_points)
      if (ImGui::Checkbox((pair.first).c_str(), &pair.second.enabled))
        this->fill_vertex_arrays__points();

    ImGui::Text("%.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
    const ImVec2 mouse_pos = ImGui::GetMousePos();

    static bool showDemoWindow = false;
    ImGui::Checkbox("show ImGUi demo window", &showDemoWindow);

    if (showDemoWindow)
      ImGui::ShowDemoWindow();

    //// Note: in older versions of ImGui, GetForegroundDrawList was called GetForegroundDrawList
    // ImGui::GetForegroundDrawList()->AddCircle(mouse_pos, (float)this->brushSize * 2.f, cursorColor, 30, 2.f);
    if (enableMaskEditing)
      ImGui::GetOverlayDrawList()->AddCircle(mouse_pos, (float)this->brushSize * 2.f, cursorColor, 30, 2.f);

    ImGui::End();
  }

  //=============================================================================

  std::string read_src_from_file(const char* filename) {
    std::ifstream fileStream;
    fileStream.open(filename, std::ios::in);
    if (!fileStream.is_open()) {
      throw std::runtime_error(std::string("Failed to open file: ") + filename);
    }
    std::stringstream buffer;
    buffer << fileStream.rdbuf();
    buffer << "\0";
    fileStream.close();
    return buffer.str();
  }

  //=============================================================================

  bool compile_shader(const char* shader_src, //
                      GLuint& shader,         //
                      GLenum shader_type) {

    shader = glCreateShader(shader_type);
    glShaderSource(shader, 1, &shader_src, nullptr);
    glCompileShader(shader);

    GLint isCompiled = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &isCompiled);
    if (isCompiled == GL_FALSE) {
      GLint maxLength = 0;
      glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &maxLength);

      // The maxLength includes the NULL character
      std::vector<GLchar> infoLog(maxLength);
      glGetShaderInfoLog(shader, maxLength, &maxLength, &infoLog[0]);
      for (const auto& c : infoLog)
        std::cout << c;

      // Exit with failure.
      glDeleteShader(shader); // Don't leak the shader.
      return false;
    }
    // Shader compilation is successful.
    return true;
  }

  //=============================================================================

  bool link_shader_program(GLuint& program,                //
                           GLuint vertexShader,            //
                           GLuint fragmentShader,          //
                           bool useGeometryShader = false, //
                           GLuint geometryShader  = 0) {

    program = glCreateProgram();

    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);
    if (useGeometryShader)
      glAttachShader(program, geometryShader);

    glLinkProgram(program);

    GLint isLinked = 0;
    glGetProgramiv(program, GL_LINK_STATUS, (int*)&isLinked);
    if (isLinked == GL_FALSE) {
      GLint maxLength = 0;
      glGetProgramiv(program, GL_INFO_LOG_LENGTH, &maxLength);

      // The maxLength includes the NULL character
      std::vector<GLchar> infoLog(maxLength);
      glGetProgramInfoLog(program, maxLength, &maxLength, &infoLog[0]);
      for (const auto& c : infoLog)
        std::cout << c;

      // We don't need the program anymore.
      glDeleteProgram(program);
      glDeleteShader(vertexShader);
      glDeleteShader(fragmentShader);
      if (useGeometryShader)
        glDeleteShader(geometryShader);

      return false;
    }

    // Always detach shaders after a successful link.
    glDetachShader(program, vertexShader);
    glDetachShader(program, fragmentShader);
    if (useGeometryShader)
      glDetachShader(program, geometryShader);

    return true;
  }

  //=============================================================================

  bool compile_shader_program(GLuint& program_id,                 //
                              const char* vertexShaderFilename,   //
                              const char* fragmentShaderFilename, //
                              const char* geometryShaderFilename = "") {

    GLuint vertexShader;
    GLuint fragmentShader;
    GLuint geometryShader = 0;

    std::string vertexShaderSrc = read_src_from_file(vertexShaderFilename);
    if (compile_shader(vertexShaderSrc.c_str(), vertexShader, GL_VERTEX_SHADER)) {
      // std::cout << "vert shader: ok\n";
    } else {
      std::cerr << "vert shader: error! aborting\n";
      return false;
    }

    std::string fragmentShaderSrc = read_src_from_file(fragmentShaderFilename);
    if (compile_shader(fragmentShaderSrc.c_str(), fragmentShader, GL_FRAGMENT_SHADER)) {
      // std::cout << "frag shader: ok\n";
    } else {
      std::cerr << "frag shader: error! aborting\n";
      return false;
    }

    bool useGeometryShader = geometryShaderFilename && geometryShaderFilename[0];
    if (useGeometryShader) {
      std::string geometryShaderSrc = read_src_from_file(geometryShaderFilename);
      if (compile_shader(geometryShaderSrc.c_str(), geometryShader, GL_GEOMETRY_SHADER)) {
        // std::cout << "geom shader: ok\n";
      } else {
        std::cerr << "geom shader: error! aborting\n";
        return false;
      }
    }

    if (link_shader_program(program_id, vertexShader, fragmentShader, useGeometryShader, geometryShader)) {
      // std::cout << "shader program: ok\n";
      return true;
    }

    std::cerr << "shader program: error! aborting\n";
    return false;
  }

  //=============================================================================

  static void glfw_error_callback(int error, const char* description) {
    std::cout << "Glfw error " << error << " : " << description << std::endl;
  }

  //=============================================================================

  void Gui2::init_grid_texture() {
    _grid = cv::Mat(100, 100, CV_8UC3);
    _grid.setTo(cv::Vec3b(255, 255, 255));
    cv::copyMakeBorder(_grid, _grid,                       //
                       gridThickness, gridThickness, 0, 0, //
                       CV_HAL_BORDER_CONSTANT,             //
                       cv::Vec3b(255, 185, 0));
    cv::copyMakeBorder(_grid, _grid,                       //
                       0, 0, gridThickness, gridThickness, //
                       CV_HAL_BORDER_CONSTANT,             //
                       cv::Vec3b(0, 0, 255));
    set_texture_from_image(this->textureID[2], _grid);
  }

  //=============================================================================

  bool Gui2::init_window() {

    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
      return false;

    glfwWindowHint(GLFW_SAMPLES, 8);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // 3.2+ only
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);           // Required on Mac
#endif
    // Create window with graphics context
    window = glfwCreateWindow(this->window_w, this->window_h, "Simple Things", nullptr, nullptr);
    if (!window)
      return false;
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
      std::cerr << "error: failed to initialize OpenGL loader!" << std::endl;
      return false;
    }

    //// debug info
    printf("OpenGL Version %d.%d loaded\n", GLVersion.major, GLVersion.minor);
    // int major, minor, rev;
    // major = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MAJOR);
    // minor = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MINOR);
    // rev   = glfwGetWindowAttrib(window, GLFW_CONTEXT_REVISION);
    // printf("OpenGL version received: %d.%d.%d\n", major, minor, rev);
    // printf("Supported OpenGL is %s\n", (const char*)glGetString(GL_VERSION));
    // printf("Supported GLSL is %s\n", (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION));

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    (void)io;
    io.Fonts->AddFontFromFileTTF(IGSV_GUI2_DIR "/fonts/inconsolata.otf", 20);

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();

    // Setup Platform/Renderer bindings
    const char* glsl_version = "#version 330";
    ImGui_ImplGlfw_InitForOpenGL(this->window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

    return true;
  }

  //=============================================================================

  void Gui2::init_textures() {
    glGenTextures(3, &this->textureID[0]);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, this->textureID[0]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, this->textureID[1]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, this->textureID[2]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  }

  //=============================================================================

  void Gui2::init_shader_uniforms() {
    this->locModel = glGetUniformLocation(this->programID, "model");
    this->locView  = glGetUniformLocation(this->programID, "view");
    this->locProj  = glGetUniformLocation(this->programID, "proj");

    this->locImageTexture = glGetUniformLocation(this->programID, "imageTexture");
    this->locMaskTexture  = glGetUniformLocation(this->programID, "maskTexture");
    this->locGridTexture  = glGetUniformLocation(this->programID, "gridTexture");

    this->locEnableImageTexture = glGetUniformLocation(this->programID, "enableImageTexture");
    this->locEnableMaskTexture  = glGetUniformLocation(this->programID, "enableMaskTexture");
    this->locEnableGridTexture  = glGetUniformLocation(this->programID, "enableGridTexture");

    this->locModel_lines         = glGetUniformLocation(this->programID_lines, "model");
    this->locView_lines          = glGetUniformLocation(this->programID_lines, "view");
    this->locProj_lines          = glGetUniformLocation(this->programID_lines, "proj");
    this->locLineThickness_lines = glGetUniformLocation(this->programID_lines, "lineThickness");

    this->locModel_points     = glGetUniformLocation(this->programID_points, "model");
    this->locView_points      = glGetUniformLocation(this->programID_points, "view");
    this->locProj_points      = glGetUniformLocation(this->programID_points, "proj");
    this->locPointSize_points = glGetUniformLocation(this->programID_points, "pointSize");
  }

  //=============================================================================

  void Gui2::init_vertex_arrays() {

    // ----- triangles -----

    glGenVertexArrays(1, &this->vao);
    glBindVertexArray(this->vao);
    glGenBuffers(1, &this->vbo_positions);
    glGenBuffers(1, &this->vbo_colors);
    glGenBuffers(1, &this->vbo_uvcoords);
    glGenBuffers(1, &this->vbo_triangles);

    glGenVertexArrays(1, &this->vao_base);
    glBindVertexArray(this->vao_base);
    glGenBuffers(1, &this->vbo_base_positions);
    glGenBuffers(1, &this->vbo_base_colors);
    glGenBuffers(1, &this->vbo_base_uvcoords);
    glGenBuffers(1, &this->vbo_base_triangles);

    // ----- edges -----

    glGenVertexArrays(1, &this->vao_lines);
    glGenBuffers(1, &this->vbo_lines_positions);
    glGenBuffers(1, &this->vbo_lines_colors);
    glGenBuffers(1, &this->vbo_lines_edges);

    // ----- points -----

    glGenVertexArrays(1, &this->vao_points);
    glGenBuffers(1, &this->vbo_points_positions);
    glGenBuffers(1, &this->vbo_points_colors);
  }

  //=============================================================================

  void Gui2::fill_vertex_arrays__triangles() {

    n_triangles           = 0;
    int n_triangle_points = 0;

    for (const auto& _pair : _render_triangles)
      if (_pair.second.enabled) {
        n_triangle_points += _pair.second.xy.size();
        n_triangles += _pair.second.idx.size();
      }

    if (n_triangles == 0)
      return;

    IGSV::drawable::triangles _triangles(n_triangle_points, n_triangles);

    int v = 0;
    int u = 0;
    int c = 0;
    int e = 0;
    int k = 0;

    for (const auto& _pair : _render_triangles)
      if (_pair.second.enabled) {

        const auto& _t = _pair.second;

        for (const auto& xy : _t.xy)
          _triangles.xy[v++] = xy;

        for (const auto& rgb : _t.rgb)
          _triangles.rgb[c++] = rgb;

        for (const auto& uv : _t.uv)
          _triangles.uv[u++] = uv;

        for (const auto& idx : _t.idx) {
          _triangles.idx[e][0] = idx[0] + k;
          _triangles.idx[e][1] = idx[1] + k;
          _triangles.idx[e][2] = idx[2] + k;
          e++;
        }

        k += _t.xy.size();
      }

    glBindVertexArray(this->vao);
    glBindBuffer(GL_ARRAY_BUFFER, this->vbo_positions);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 2 * n_triangle_points, _triangles.xy.data(), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, this->vbo_colors);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * n_triangle_points, _triangles.rgb.data(), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, this->vbo_uvcoords);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 2 * n_triangle_points, _triangles.uv.data(), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->vbo_triangles);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned) * 3 * n_triangles, _triangles.idx.data(), GL_DYNAMIC_DRAW);

    // ----- Render : base mesh (two triangles) -----
    static const int n_base_vertices  = 4;
    static const int n_base_triangles = 2;
    IGSV::drawable::triangles _base(n_base_vertices, n_base_triangles);
    _base.xy[0]  = { 0.0f, 0.0f };
    _base.xy[1]  = { (float)_input.cols, 0.0f };
    _base.xy[2]  = { (float)_input.cols, (float)_input.rows };
    _base.xy[3]  = { 0.0f, (float)_input.rows };
    _base.uv[0]  = { 0.0f, 0.0f };
    _base.uv[1]  = { 1.0f, 0.0f };
    _base.uv[2]  = { 1.0f, 1.0f };
    _base.uv[3]  = { 0.0f, 1.0f };
    _base.rgb[0] = { 1.0f, 1.0f, 1.0f };
    _base.rgb[1] = { 1.0f, 1.0f, 1.0f };
    _base.rgb[2] = { 1.0f, 1.0f, 1.0f };
    _base.rgb[3] = { 1.0f, 1.0f, 1.0f };
    _base.idx[0] = { 0, 1, 2 };
    _base.idx[1] = { 2, 3, 0 };
    glBindVertexArray(this->vao_base);
    glBindBuffer(GL_ARRAY_BUFFER, this->vbo_base_positions);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 2 * n_base_vertices, _base.xy.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, this->vbo_base_colors);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * n_base_vertices, _base.rgb.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, this->vbo_base_uvcoords);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 2 * n_base_vertices, _base.uv.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->vbo_base_triangles);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned) * 3 * n_base_triangles, _base.idx.data(), GL_STATIC_DRAW);
  }

  //=============================================================================

  void Gui2::fill_vertex_arrays__edges() {

    n_edges           = 0;
    int n_edge_points = 0;

    for (const auto& _pair : _render_edges)
      if (_pair.second.enabled) {
        n_edge_points += _pair.second.xy.size();
        n_edges += _pair.second.idx.size();
      }

    if (n_edges == 0)
      return;

    IGSV::drawable::edges _edges(n_edge_points, n_edges);

    int v = 0;
    int c = 0;
    int e = 0;
    int k = 0;

    for (const auto& _pair : _render_edges)
      if (_pair.second.enabled) {

        const auto& _e = _pair.second;

        for (const auto& xy : _e.xy)
          _edges.xy[v++] = xy;

        for (const auto& rgb : _e.rgb)
          _edges.rgb[c++] = rgb;

        for (const auto& idx : _e.idx) {
          _edges.idx[e][0] = idx[0] + k;
          _edges.idx[e][1] = idx[1] + k;
          e++;
        }

        k += _e.xy.size();
      }

    glBindVertexArray(this->vao_lines);
    glBindBuffer(GL_ARRAY_BUFFER, this->vbo_lines_positions);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 2 * n_edge_points, _edges.xy.data(), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, this->vbo_lines_colors);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * n_edge_points, _edges.rgb.data(), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->vbo_lines_edges);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned) * 2 * n_edges, _edges.idx.data(), GL_DYNAMIC_DRAW);
  }

  //=============================================================================

  void Gui2::fill_vertex_arrays__points() {

    n_points = 0;

    for (const auto& _pair : _render_points)
      if (_pair.second.enabled)
        n_points += _pair.second.xy.size();

    if (n_points == 0)
      return;

    IGSV::drawable::points _points(n_points);

    int v = 0;
    int c = 0;
    int k = 0;

    for (const auto& _pair : _render_points)
      if (_pair.second.enabled) {

        const auto& _p = _pair.second;

        for (const auto& xy : _p.xy)
          _points.xy[v++] = xy;

        for (const auto& rgb : _p.rgb)
          _points.rgb[c++] = rgb;

        k += _p.xy.size();
      }

    glBindVertexArray(this->vao_points);
    glBindBuffer(GL_ARRAY_BUFFER, this->vbo_points_positions);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 2 * n_points, _points.xy.data(), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, this->vbo_points_colors);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * n_points, _points.rgb.data(), GL_DYNAMIC_DRAW);
  }

  //=============================================================================

  void Gui2::initialize() {

    this->tx = (float)_input.cols * -0.5f;
    this->ty = (float)_input.rows * -0.5f;

    if (!this->init_window())
      throw std::runtime_error(std::string("Failed to create glfw window"));

    if (!compile_shader_program(this->programID,                      //
                                IGSV_GUI2_DIR "/glsl/mesh.vert.glsl", //
                                IGSV_GUI2_DIR "/glsl/mesh.frag.glsl"))
      throw std::runtime_error(std::string("Failed to compile mesh shader program"));

    if (!compile_shader_program(this->programID_lines,                 //
                                IGSV_GUI2_DIR "/glsl/lines.vert.glsl", //
                                IGSV_GUI2_DIR "/glsl/lines.frag.glsl", //
                                IGSV_GUI2_DIR "/glsl/lines.geom.glsl"))
      throw std::runtime_error(std::string("Failed to compile lines shader program"));

    if (!compile_shader_program(this->programID_points,                 //
                                IGSV_GUI2_DIR "/glsl/points.vert.glsl", //
                                IGSV_GUI2_DIR "/glsl/points.frag.glsl"))
      throw std::runtime_error(std::string("Failed to compile points shader program"));

    glEnable(GL_PROGRAM_POINT_SIZE);

    this->init_textures();
    set_texture_from_image(this->textureID[0], _input);
    set_texture_from_image(this->textureID[1], _mask);
    init_grid_texture();

    this->init_shader_uniforms();
    this->init_vertex_arrays();
  }

  //=============================================================================

  void Gui2::mouse_events() {

    if (!ImGui::GetIO().WantCaptureMouse) {

      if (ImGui::IsMouseDragging(0)) {
        this->tx += ImGui::GetIO().MouseDelta.x * 0.005f / this->zoom;
        this->ty -= ImGui::GetIO().MouseDelta.y * 0.005f / this->zoom;
      }

      if (ImGui::IsMouseReleased(/*ImGuiMouseButton_Right*/ 1)) {
        isMaskBeingEdited = false;
      }

      if (ImGui::IsMouseClicked(/*ImGuiMouseButton_Right*/ 1) && enableMaskEditing) {
        isMaskBeingEdited  = true;
        nBrushStrokePoints = 0;
      }

      if (isMaskBeingEdited) {

        nBrushStrokePoints++;
        // const auto& brushColor = ImGui::GetIO().KeyCtrl ? IGSV::colors::transparent
        //                                                 : (ImGui::GetIO().KeyShift //
        //                                                    ? IGSV::colors::orange
        //                                                    : IGSV::colors::violet);
        const auto& brushColor = ImGui::GetIO().KeyCtrl ? IGSV::colors::transparent
                                                        : (ImGui::GetIO().KeyShift //
                                                               ? IGSV::colors::black
                                                               : IGSV::colors::white);
        int brushSizeReal = std::round<int>((float)brushSize * 0.01f / this->zoom);
        if (brushSizeReal < 1)
          brushSizeReal = 1;

        this->unproject(ImGui::GetMousePos().x, ImGui::GetMousePos().y, this->x1, this->y1);

        if (nBrushStrokePoints > 1)
          cv::line(_mask,                         //
                   cv::Point(this->x0, this->y0), //
                   cv::Point(this->x1, this->y1), //
                   brushColor,                    //
                   brushSizeReal,                 //
                   cv::LINE_AA, 0);
        else
          cv::circle(_mask,                         //
                     cv::Point(this->x1, this->y1), //
                     brushSizeReal,                 //
                     brushColor,                    //
                     cv::FILLED,                    //
                     cv::LINE_AA);

        this->x0 = this->x1;
        this->y0 = this->y1;

        set_texture_from_image(this->textureID[1], _mask);
      }

      this->zoom *= std::pow(1.025, ImGui::GetIO().MouseWheel);
    }
  }

  //=============================================================================

  void Gui2::keyboard_events() {
    if (!ImGui::GetIO().WantCaptureKeyboard) {
      if (ImGui::IsKeyPressed(GLFW_KEY_ESCAPE)) {
        glfwSetWindowShouldClose(this->window, GLFW_TRUE);
      }
    }
  }

  //=============================================================================

  void Gui2::recompute_matrices() {
    proj.setIdentity();
    proj(0, 0) = (float)window_h / (float)window_w;
    proj(1, 1) = 1.0f;
    proj(2, 2) = -2.f / (farVal - nearVal);
    proj(2, 3) = -(farVal + nearVal) / (farVal - nearVal);

    view.setIdentity();
    view.topLeftCorner<3, 3>() *= zoom;

    model.setIdentity();
    model(3, 0) = tx;
    model(3, 1) = ty;
    // model(0, 3) = tx;
    // model(1, 3) = ty;
  }

  //=============================================================================

  void Gui2::render() {

    this->fill_vertex_arrays__triangles();
    this->fill_vertex_arrays__edges();
    this->fill_vertex_arrays__points();

    int framebuffer_w, framebuffer_h;
    while (!glfwWindowShouldClose(this->window) && _callback_success) {

      //----- process events -----

      glfwPollEvents();
      this->mouse_events();
      this->keyboard_events();
      this->recompute_matrices();

      //----- init new frame -----

      ImGui_ImplOpenGL3_NewFrame();
      ImGui_ImplGlfw_NewFrame();
      ImGui::NewFrame();
      this->render_custom_menu();
      ImGui::Render();

      //----- init opengl -----

      glfwGetWindowSize(this->window, &this->window_w, &this->window_h);
      glfwGetFramebufferSize(this->window, &framebuffer_w, &framebuffer_h);
      glViewport(0, 0, framebuffer_w, framebuffer_h);
      glClearColor(this->bgColor[0], //
                   this->bgColor[1], //
                   this->bgColor[2], //
                   this->bgColor[3]);
      glClear(GL_COLOR_BUFFER_BIT);

      // std::cout << std::endl;

      //----- draw triangles -----

      {

        // std::cout << "#T = " << n_triangles << std::endl;

        glUseProgram(this->programID);

        glUniformMatrix4fv(this->locModel, 1, GL_TRUE, this->model.data());
        glUniformMatrix4fv(this->locView, 1, GL_TRUE, this->view.data());
        glUniformMatrix4fv(this->locProj, 1, GL_TRUE, this->proj.data());

        glUniform1i(this->locImageTexture, 0); // GL_TEXTURE0
        glUniform1i(this->locMaskTexture, 1);  // GL_TEXTURE1
        glUniform1i(this->locGridTexture, 2);  // GL_TEXTURE2

        glEnableVertexAttribArray(0);
        glEnableVertexAttribArray(1);
        glEnableVertexAttribArray(2);

        //// base mesh
        glUniform1i(this->locEnableImageTexture, this->enableImageTexture);
        glUniform1i(this->locEnableMaskTexture, this->enableMaskEditing);
        glUniform1i(this->locEnableGridTexture, 0);

        glBindBuffer(GL_ARRAY_BUFFER, this->vbo_base_positions);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);

        glBindBuffer(GL_ARRAY_BUFFER, this->vbo_base_colors);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);

        glBindBuffer(GL_ARRAY_BUFFER, this->vbo_base_uvcoords);
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);

        glBindTexture(GL_TEXTURE_2D, this->textureID[0]); // input
        glBindTexture(GL_TEXTURE_2D, this->textureID[1]); // mask
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->vbo_base_triangles);

        glDrawElements(GL_TRIANGLES, 2 * 3, GL_UNSIGNED_INT, (void*)0);

        // narrow band mesh
        glUniform1i(this->locEnableImageTexture, 0);
        glUniform1i(this->locEnableMaskTexture, 0);
        glUniform1i(this->locEnableGridTexture, 1);

        glBindBuffer(GL_ARRAY_BUFFER, this->vbo_positions);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);

        glBindBuffer(GL_ARRAY_BUFFER, this->vbo_colors);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);

        glBindBuffer(GL_ARRAY_BUFFER, this->vbo_uvcoords);
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);

        glBindTexture(GL_TEXTURE_2D, this->textureID[2]); // grid
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->vbo_triangles);

        if (showWireframe) {
          glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
          glDrawElements(GL_TRIANGLES, n_triangles * 3, GL_UNSIGNED_INT, (void*)0);
          glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        } else {
          glEnable(GL_BLEND);
          glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
          glDrawElements(GL_TRIANGLES, n_triangles * 3, GL_UNSIGNED_INT, (void*)0);
          glDisable(GL_BLEND);
        }

        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);
        glDisableVertexAttribArray(2);
      }

      //----- draw edges -----

      if (n_edges > 0) {

        // std::cout << "#E = " << n_edges << std::endl;

        glUseProgram(this->programID_lines);

        glUniformMatrix4fv(this->locModel_lines, 1, GL_TRUE, this->model.data());
        glUniformMatrix4fv(this->locView_lines, 1, GL_TRUE, this->view.data());
        glUniformMatrix4fv(this->locProj_lines, 1, GL_TRUE, this->proj.data());

        glUniform1f(this->locLineThickness_lines, this->lineThickness);

        glEnableVertexAttribArray(0);
        glEnableVertexAttribArray(1);

        glBindBuffer(GL_ARRAY_BUFFER, this->vbo_lines_positions);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);

        glBindBuffer(GL_ARRAY_BUFFER, this->vbo_lines_colors);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->vbo_lines_edges);
        glDrawElements(GL_LINES, n_edges * 2, GL_UNSIGNED_INT, (void*)0);

        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);
      }

      //----- draw points -----

      if (n_points > 0) {

        // std::cout << "#P = " << n_points << std::endl;

        glUseProgram(this->programID_points);

        glUniformMatrix4fv(this->locModel_points, 1, GL_TRUE, this->model.data());
        glUniformMatrix4fv(this->locView_points, 1, GL_TRUE, this->view.data());
        glUniformMatrix4fv(this->locProj_points, 1, GL_TRUE, this->proj.data());

        glUniform1f(this->locPointSize_points, this->pointSize);

        glEnableVertexAttribArray(0);
        glEnableVertexAttribArray(1);

        glBindBuffer(GL_ARRAY_BUFFER, this->vbo_points_positions);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);

        glBindBuffer(GL_ARRAY_BUFFER, this->vbo_points_colors);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);

        // glDrawElements(GL_POINTS, n_points, GL_UNSIGNED_INT, (void*)0);
        glDrawArrays(GL_POINTS, 0, n_points);

        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);
      }

      //----- draw ImGui stuff -----

      ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

      glfwSwapBuffers(this->window);
    }

    // cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glDeleteTextures(3, &this->textureID[0]);

    glDeleteBuffers(1, &this->vbo_positions);
    glDeleteBuffers(1, &this->vbo_colors);
    glDeleteBuffers(1, &this->vbo_uvcoords);
    glDeleteBuffers(1, &this->vbo_triangles);
    glDeleteVertexArrays(1, &this->vao);
    glDeleteBuffers(1, &this->vbo_base_positions);
    glDeleteBuffers(1, &this->vbo_base_colors);
    glDeleteBuffers(1, &this->vbo_base_uvcoords);
    glDeleteBuffers(1, &this->vbo_base_triangles);
    glDeleteVertexArrays(1, &this->vao_base);
    glDeleteProgram(this->programID);

    glDeleteBuffers(1, &this->vbo_lines_positions);
    glDeleteBuffers(1, &this->vbo_lines_colors);
    glDeleteBuffers(1, &this->vbo_lines_edges);
    glDeleteVertexArrays(1, &this->vao_lines);
    glDeleteProgram(this->programID_lines);

    glDeleteBuffers(1, &this->vbo_points_positions);
    glDeleteBuffers(1, &this->vbo_points_colors);
    glDeleteVertexArrays(1, &this->vao_points);
    glDeleteProgram(this->programID_points);

    glfwDestroyWindow(this->window);
    glfwTerminate();
  }

  //=============================================================================

  void Gui2::update_structures() {

    // ---- Render ---- narrow band mesh and frame field

    const size_t nf_narrow = _F_isNarrow.count();

    IGSV::drawable::triangles _mesh(nf_narrow * 3u, nf_narrow);
    IGSV::drawable::edges _field(nf_narrow * 4u, nf_narrow * 2u);

    int f = 0;
    for (int f0 = 0; f0 < _F.rows(); ++f0)
      if (_F_isNarrow(f0)) {

        // narrow band triangle

        // store vertices
        for (int k = 0; k < 3; ++k) {
          const int vi    = 3 * f + k;
          _mesh.idx[f][k] = vi;
          _mesh.rgb[vi]   = { 0.9f, 0.6f, 0.3f };
          for (int d = 0; d < 2; ++d) {
            _mesh.xy[vi][d] = _V(_F(f0, k), d);
            _mesh.uv[vi][d] = _UV(_FUV(f0, k), d);
          }
        }

        // store frame field vectors
        _field.xy[4 * f + 0][0] = _BC(f0, 0);
        _field.xy[4 * f + 0][1] = _BC(f0, 1);
        _field.xy[4 * f + 2][0] = _BC(f0, 0);
        _field.xy[4 * f + 2][1] = _BC(f0, 1);
        _field.xy[4 * f + 1][0] = _BC(f0, 0) + _X0(f0, 0);
        _field.xy[4 * f + 1][1] = _BC(f0, 1) + _X0(f0, 1);
        _field.xy[4 * f + 3][0] = _BC(f0, 0) + _X1(f0, 0);
        _field.xy[4 * f + 3][1] = _BC(f0, 1) + _X1(f0, 1);

        _field.rgb[4 * f + 0][0] = _field.rgb[4 * f + 1][0] = 1.0;
        _field.rgb[4 * f + 0][1] = _field.rgb[4 * f + 1][1] = 0.0;
        _field.rgb[4 * f + 0][2] = _field.rgb[4 * f + 1][2] = 0.0;

        _field.rgb[4 * f + 2][0] = _field.rgb[4 * f + 3][0] = 0.0;
        _field.rgb[4 * f + 2][1] = _field.rgb[4 * f + 3][1] = 0.0;
        _field.rgb[4 * f + 2][2] = _field.rgb[4 * f + 3][2] = 1.0;

        _field.idx[2 * f + 0][0] = 4 * f + 0;
        _field.idx[2 * f + 0][1] = 4 * f + 1;
        _field.idx[2 * f + 1][0] = 4 * f + 2;
        _field.idx[2 * f + 1][1] = 4 * f + 3;

        f++;
      }

    // ---- Render ---- graph vertices and edges
    IGSV::drawable::points _qverts(_QV.size());
    IGSV::drawable::edges _qedges(_QV.size(), _QE.size());
    {
      int i = 0;
      for (const auto& qv : _QV) {
        _qverts.xy[i][0]  = (float)qv.x;
        _qverts.xy[i][1]  = (float)qv.y;
        _qverts.rgb[i][0] = 0.4f;
        _qverts.rgb[i][1] = 0.2f;
        _qverts.rgb[i][2] = 0.1f;

        _qedges.xy[i][0]  = (float)qv.x;
        _qedges.xy[i][1]  = (float)qv.y;
        _qedges.rgb[i][0] = 0.6f;
        _qedges.rgb[i][1] = 0.3f;
        _qedges.rgb[i][2] = 0.2f;

        i++;
      }

      i = 0;
      for (const auto& qe : _QE) {
        _qedges.idx[i][0] = qe.qvi0;
        _qedges.idx[i][1] = qe.qvi1;
        i++;
      }
    }

    // ----- Render : Bezier curves -----
    IGSV::drawable::edges bezier_curves;
    {
      const int n_points_per_segment = 100;

      double t0, t1, t, dt, c0, c1, c2, c3, d0, d1, d2, d3, x, y;

      bool start_new_chain;

      // int row = 0;
      // int cui = 0;

      for (const auto& _chain : _CH) {

        start_new_chain = true;

        for (int i = 0; i < _chain.CubicSpline.rows(); ++i) {

          t0 = _chain.CubicSpline(i, 0);
          t1 = _chain.CubicSpline(i, 1);

          c0 = _chain.CubicSpline(i, 2);
          c1 = _chain.CubicSpline(i, 3);
          c2 = _chain.CubicSpline(i, 4);
          c3 = _chain.CubicSpline(i, 5);

          d0 = _chain.CubicSpline(i, 6);
          d1 = _chain.CubicSpline(i, 7);
          d2 = _chain.CubicSpline(i, 8);
          d3 = _chain.CubicSpline(i, 9);

          dt = (t1 - t0) / ((double)n_points_per_segment - 1.);

          t = 0;
          x = c0 + c1 * t + c2 * t * t + c3 * t * t * t;
          y = d0 + d1 * t + d2 * t * t + d3 * t * t * t;

          for (int k = 0; k < n_points_per_segment - 1; ++k) {

            if (start_new_chain) {
              bezier_curves.xy.push_back({ (float)x, (float)y });
              bezier_curves.rgb.push_back({ 0.2f, 1.0f, 0.4f });
              start_new_chain = false;
            }

            t += dt;
            x = c0 + c1 * t + c2 * t * t + c3 * t * t * t;
            y = d0 + d1 * t + d2 * t * t + d3 * t * t * t;

            bezier_curves.xy.push_back({ (float)x, (float)y });
            // if (k == n_points_per_segment - 2)
            //   bezier_curves.rgb.push_back({ 0.2f, 1.0f, 0.4f });
            // else
            bezier_curves.rgb.push_back({ 0.2f, 1.0f, 0.4f });
            bezier_curves.idx.push_back({ (unsigned)bezier_curves.xy.size() - 2, //
                                          (unsigned)bezier_curves.xy.size() - 1 });
          }
        }
      }
    }

    this->add_triangles("mesh", _mesh);
    this->add_edges("field", _field);
    this->add_points("verts", _qverts);
    this->add_edges("edges", _qedges);
    this->add_edges("curves", bezier_curves);
  }

  //=============================================================================

} // namespace IGSV
