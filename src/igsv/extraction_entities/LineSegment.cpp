#include <igsv/extraction_entities/LineSegment.h>

#include <igsv/common/CornerType.h>
#include <igsv/common/defs.h>

#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>

#include <sstream>

using IK        = CGAL::Simple_cartesian<double>;                         // inexact predicates
using EK        = CGAL::Simple_cartesian<CGAL::Quotient<CGAL::MP_Float>>; // exact predicates
using IK_to_EK  = CGAL::Cartesian_converter<IK, EK>;
using EK_to_IK  = CGAL::Cartesian_converter<EK, IK>;
using TriCoords = CGAL::Barycentric_coordinates::Triangle_coordinates_2<IK>;

namespace IGSV {

  // ===========================================================================

  const std::string LineSegment::str() const {
    std::stringstream ss;
    ss << "#" << id << " [f=" << f << " | dir=" << dir << " | iso=" << iso << "]";
    return ss.str();
  }

  // ===========================================================================

  LineSegment::LineSegment(const Eigen::MatrixXd& _V,   //
                           const Eigen::MatrixXi& _F,   //
                           const Eigen::MatrixXd& _UV,  //
                           const Eigen::MatrixXi& _FUV, //
                           int f,                       //
                           int iso,                     //
                           int dir,                     //
                           int id)
      : f(f), iso(iso), dir(dir), id(id) {

    // 0.0 construct the triangle
    const int c0 = _FUV(f, 0);
    const int c1 = _FUV(f, 1);
    const int c2 = _FUV(f, 2);

    const EK::Triangle_2 tri_uv_EK(          //
        EK::Point_2(_UV(c0, 0), _UV(c0, 1)), //
        EK::Point_2(_UV(c1, 0), _UV(c1, 1)), //
        EK::Point_2(_UV(c2, 0), _UV(c2, 1)));

    const IK::Triangle_2 tri_uv_IK(          //
        IK::Point_2(_UV(c0, 0), _UV(c0, 1)), //
        IK::Point_2(_UV(c1, 0), _UV(c1, 1)), //
        IK::Point_2(_UV(c2, 0), _UV(c2, 1)));

    const int v0 = _F(f, 0);
    const int v1 = _F(f, 1);
    const int v2 = _F(f, 2);

    const IK::Triangle_2 tri_xy_IK(        //
        IK::Point_2(_V(v0, 0), _V(v0, 1)), //
        IK::Point_2(_V(v1, 0), _V(v1, 1)), //
        IK::Point_2(_V(v2, 0), _V(v2, 1)));

    // 0.1 construct a line given by the equation a*x + b*y + c = 0
    double a = 0, b = 0, c = (double)(-iso);
    if (dir == CORNER_TYPE_INTEGER_V) {
      b = 1; // y = iso  >>>  0*x + 1*y - iso = 0
    } else {
      a = 1; // x = iso  >>>  1*x + 0*y - iso = 0
    }
    const EK::Line_2 line(a, b, c);

    // 0.2 intersect the triangle and the line
    // result can be : o) empty, i) point, ii) segment
    // here, auto := CGAL::cpp11::result_of<EK::Intersect_2(EK::Triangle_2, EK::Line_2)>::type
    auto result = intersection(tri_uv_EK, line);

    // 0.3 get barycentric coords
    TriCoords tri_coords(tri_uv_IK[0], tri_uv_IK[1], tri_uv_IK[2]);
    // tri_coords.print_information(); // for debugging
    if (result) {
      // result is non-empty
      if (const EK::Segment_2* s_EK = boost::get<EK::Segment_2>(&*result)) {
        // result is a segment
        numpts = 2;
        EK_to_IK to_inexact;
        IK::Segment_2 s_IK = to_inexact(*s_EK);
        bc.reserve(6);
        tri_coords(s_IK[0], std::inserter(bc, bc.end()));
        tri_coords(s_IK[1], std::inserter(bc, bc.end()));
      } else {
        // result is a point
        const EK::Point_2* p_EK = boost::get<EK::Point_2>(&*result);

        numpts = 1;
        EK_to_IK to_inexact;
        IK::Point_2 p_IK = to_inexact(*p_EK);
        bc.reserve(3);
        tri_coords(p_IK, std::inserter(bc, bc.end()));
      }
    } else {
      // result is empty
      DEBUG_PRINT_INFO("f=%d : EMPTY intersection", f);
      numpts = 0;
    }

    for (int i = 0; i < numpts; ++i) {
      double px = 0, py = 0;
      for (int k = 0; k < 3; ++k) {
        px += bc[3 * i + k] * tri_xy_IK[k].x();
        py += bc[3 * i + k] * tri_xy_IK[k].y();
      }
      xyz.push_back(IK::Point_2(px, py));
    }
  }

  // ===========================================================================

  void LineMap::clear() {
    listU.clear();
    listV.clear();
    segmentsU.clear();
    segmentsV.clear();
  }

  // ===========================================================================

  bool LineMap::compute(const Eigen::MatrixXd& _V,   //
                        const Eigen::MatrixXi& _F,   //
                        const Eigen::MatrixXd& _UV,  //
                        const Eigen::MatrixXi& _FUV, //
                        const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow) {

    clear();

    listU.resize(_F.rows());
    listV.resize(_F.rows());

    // stats, number of intersections that are:
    //    0 : empty
    //    1 : point
    //    2 : segment
    int num_intersections[3] = { 0, 0, 0 };

    // Tuv : rows are uv coords of triangle's vertices
    Eigen::Matrix<double, 3, 2> Tuv;
    double u_lower, u_upper, v_lower, v_upper;

    int id = 0;
    for (unsigned f = 0; f < _F.rows(); f++) {

      listU[f].clear();
      listV[f].clear();

      if (!_F_isNarrow(f))
        continue;

      // get bounds
      for (int k = 0; k < 3; k++)
        Tuv.row(k) << _UV.row(_FUV(f, k));
      u_lower = std::ceil(Tuv.col(0).minCoeff());
      u_upper = std::floor(Tuv.col(0).maxCoeff()) + 1;
      v_lower = std::ceil(Tuv.col(1).minCoeff());
      v_upper = std::floor(Tuv.col(1).maxCoeff()) + 1;

      // loop over u
      for (int u = u_lower; u < u_upper; u++) {
        segmentsU.push_back(LineSegment(_V, _F, _UV, _FUV, f, u, 1, id++));
        listU[f].push_back(segmentsU.size() - 1);
        num_intersections[segmentsU.back().numpts]++;
      }

      // loop over v
      for (int v = v_lower; v < v_upper; v++) {
        segmentsV.push_back(LineSegment(_V, _F, _UV, _FUV, f, v, 2, id++));
        listV[f].push_back(segmentsV.size() - 1);
        num_intersections[segmentsV.back().numpts]++;
      }

    } // end loop over faces

    // std::cout << "  # of intersections:" << std::endl
    //           << "    " << num_intersections[0] << " empty" << std::endl
    //           << "    " << num_intersections[1] << " point" << std::endl
    //           << "    " << num_intersections[2] << " segment" << std::endl;

    return true;
  }

  // ===========================================================================

} // namespace IGSV
