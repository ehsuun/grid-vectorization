//=============================================================================
//
//  CLASS : LineSegment
//          LineMap
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <CGAL/Simple_cartesian.h>
#include <vector>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== CLASS DEFINITION ====================================================

#define MAX_TRACING_ITER 1e6

  class LineSegment {
  public:
    LineSegment(const Eigen::MatrixXd& _V,   //
                const Eigen::MatrixXi& _F,   //
                const Eigen::MatrixXd& _UV,  //
                const Eigen::MatrixXi& _FUV, //
                int f,                       //
                int iso,                     //
                int dir,                     //
                int id);

    const int id;
    int f;
    int iso;
    int dir;
    int compid  = -1;
    int numpts  = 0;
    bool traced = false;

    const std::string str() const;
    std::vector<CGAL::Simple_cartesian<double>::Point_2> xyz;

  private:
    std::vector<double> bc;
  };

  //== CLASS DEFINITION ====================================================

  class LineMap {
  public:
    void clear();
    bool compute(const Eigen::MatrixXd& _V,   //
                 const Eigen::MatrixXi& _F,   //
                 const Eigen::MatrixXd& _UV,  //
                 const Eigen::MatrixXi& _FUV, //
                 const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow);

    std::vector<LineSegment> segmentsU, segmentsV;
    std::vector<std::vector<int>> listU, listV;
  };

} // namespace IGSV
