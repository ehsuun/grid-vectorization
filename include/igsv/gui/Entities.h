//=============================================================================
//
//  CLASS : GuiFColors
//          GuiEntity
//          GuiPoints
//          GuiLabels
//          GuiEdges
//          GuiCurve
//          GuiVectors
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <Eigen/Core>
#include <string>
#include <vector>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== CLASS DEFINITION =======================================================

  class GuiFColors {
  public:
    Eigen::MatrixXd clr;
  };

  //== CLASS DEFINITION =======================================================

  class GuiEntity { // abstract entity
  public:
    GuiEntity(std::string name, unsigned char key, int modifier)
        : is_shown(false), name(name), key(key), modifier(modifier) {
      clr.resize(0, 3);
    }
    bool is_shown;
    std::string name;
    int key;
    int modifier;
    Eigen::MatrixXd clr;

    void show() { is_shown = true; }
    void hide() { is_shown = false; }
    void toggle() { is_shown = !is_shown; }

    std::string key_name() const;
  };

  //== CLASS DEFINITION =======================================================

  class GuiPoints : public GuiEntity {
  public:
    GuiPoints(std::string name, unsigned char key = -1, int modifier = 0) : GuiEntity(name, key, modifier) {
      pos.resize(0, 3);
    }
    Eigen::MatrixXd pos;
  };

  //== CLASS DEFINITION =======================================================

  class GuiLabels : public GuiPoints {
  public:
    GuiLabels(std::string name, unsigned char key = -1, int modifier = 0) : GuiPoints(name, key, modifier) {
      str.clear();
    }
    std::vector<std::string> str;
    bool show_points = false;
  };

  //== CLASS DEFINITION =======================================================

  class GuiEdges : public GuiEntity {
  public:
    GuiEdges(std::string name, unsigned char key = -1, int modifier = 0) : GuiEntity(name, key, modifier) {
      pos0.resize(0, 3);
      pos1.resize(0, 3);
    }
    Eigen::MatrixXd pos0, pos1;
  };

  class GuiCurve : public GuiEntity {
  public:
    GuiCurve(std::string name, unsigned char key = -1, int modifier = 0) : GuiEntity(name, key, modifier) {
      pos.resize(0, 3);
    }
    Eigen::MatrixXd pos;
  };

  //== CLASS DEFINITION =======================================================

  class GuiVectors : public GuiEntity {
  public:
    GuiVectors(std::string name, unsigned char key = -1, int modifier = 0) : GuiEntity(name, key, modifier) {
      pos.resize(0, 3);
      vec.resize(0, 3);
    }
    Eigen::MatrixXd pos, vec;
  };

}; // namespace IGSV
