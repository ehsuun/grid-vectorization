#include <igsv/extraction/exact_predicates.h>

namespace IGSV {

  // ===========================================================================

  CGAL::Orientation ORIENT2D(const Point_2& a, const Point_2& b, const Point_2& c) {
    return CGAL::orientation(a, b, c);
  }

  CGAL::Orientation ORIENT2D(const complexd& a, const complexd& b, const complexd& c) {
    return ORIENT2D(Point_2(a.real(), a.imag()), //
                    Point_2(b.real(), b.imag()), //
                    Point_2(c.real(), c.imag()));
  }

  // ===========================================================================

  bool IS_CW(const Point_2& a, const Point_2& b, const Point_2& c) {
    return ORIENT2D(a, b, c) == CGAL::Orientation::NEGATIVE;
  }

  bool IS_CW(const complexd& a, const complexd& b, const complexd& c) {
    return IS_CW(Point_2(a.real(), a.imag()), //
                 Point_2(b.real(), b.imag()), //
                 Point_2(c.real(), c.imag()));
  }

  // ===========================================================================

  bool IS_CCW(const Point_2& a, const Point_2& b, const Point_2& c) {
    return ORIENT2D(a, b, c) == CGAL::Orientation::POSITIVE;
  }

  bool IS_CCW(const complexd& a, const complexd& b, const complexd& c) {
    return IS_CCW(Point_2(a.real(), a.imag()), //
                  Point_2(b.real(), b.imag()), //
                  Point_2(c.real(), c.imag()));
  }

  // ===========================================================================

  bool IS_COLLINEAR(const Point_2& a, const Point_2& b, const Point_2& c) {
    return ORIENT2D(a, b, c) == CGAL::Orientation::ZERO;
  }

  bool IS_COLLINEAR(const Vector_2& u, const Vector_2& v) {
    return ORIENT2D(Point_2(0, 0), Point_2(u.x(), u.y()), Point_2(v.x(), v.y())) == CGAL::Orientation::ZERO;
  }

  bool IS_COLLINEAR(const complexd& a, const complexd& b, const complexd& c) {
    return IS_COLLINEAR(Point_2(a.real(), a.imag()), //
                        Point_2(b.real(), b.imag()), //
                        Point_2(c.real(), c.imag()));
  }

  bool IS_COLLINEAR(const complexd& u, const complexd& v) {
    return IS_COLLINEAR(Vector_2(u.real(), u.imag()), //
                        Vector_2(v.real(), v.imag()));
  }

  bool IS_COLLINEAR(const complexd& u, const complexi& v) {
    return IS_COLLINEAR(Vector_2(u.real(), u.imag()), //
                        Vector_2(v.real(), v.imag()));
  }

  // ===========================================================================

  bool POINTS_INTO(const Vector_2& d, const Point_2& a, const Point_2& b, const Point_2& c) {
    return IS_CCW(a, b, a + d) & IS_CCW(a, a + d, c);
  }

  bool POINTS_INTO(const complexi& d, const complexd& a, const complexd& b, const complexd& c) {
    return POINTS_INTO(Vector_2(d.real(), d.imag()), //
                       Point_2(a.real(), a.imag()),  //
                       Point_2(b.real(), b.imag()),  //
                       Point_2(c.real(), c.imag()));
  }

  // ===========================================================================

} // namespace IGSV
