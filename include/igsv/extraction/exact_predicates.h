//=============================================================================
//
//  wrapper for exact predicates from CGAL :
//  -- ORIENT2D(a, b, c)         :=  orientation of triangle (a, b, c)
//                                  [NEGATIVE, ZERO, POSITIVE]
//  -- IS_CW(a, b, c)            :=  ORIENT2D(a, b, c) == NEGATIVE
//  -- IS_CCW(a, b, c)           :=  ORIENT2D(a, b, c) == POSITIVE
//  -- IS_COLLINEAR(a, b, c)     :=  ORIENT2D(a, b, c) == ZERO
//  -- POINTS_INTO(a, b, c, d)   :=  IS_CCW(a, b, a + d) & IS_CCW(a, a + d, c)
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <complex>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== TYPEDEFS =================================================================

  using EP_IC_K  = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Point_2  = EP_IC_K::Point_2;
  using Vector_2 = EP_IC_K::Vector_2;

  using complexd = std::complex<double>;
  using complexi = std::complex<int>;

  //== FUNCTION DEFINITIONS ===================================================

  CGAL::Orientation ORIENT2D(const Point_2& a, const Point_2& b, const Point_2& c);
  CGAL::Orientation ORIENT2D(const complexd& a, const complexd& b, const complexd& c);

  bool IS_CW(const Point_2& a, const Point_2& b, const Point_2& c);
  bool IS_CW(const complexd& a, const complexd& b, const complexd& c);

  bool IS_CCW(const Point_2& a, const Point_2& b, const Point_2& c);
  bool IS_CCW(const complexd& a, const complexd& b, const complexd& c);

  bool IS_COLLINEAR(const Point_2& a, const Point_2& b, const Point_2& c);
  bool IS_COLLINEAR(const Vector_2& u, const Vector_2& v);
  bool IS_COLLINEAR(const complexd& a, const complexd& b, const complexd& c);
  bool IS_COLLINEAR(const complexd& u, const complexd& v);
  bool IS_COLLINEAR(const complexd& u, const complexi& v);

  bool POINTS_INTO(const Vector_2& d, const Point_2& a, const Point_2& b, const Point_2& c);
  bool POINTS_INTO(const complexi& d, const complexd& a, const complexd& b, const complexd& c);

} // namespace IGSV
