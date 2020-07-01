#pragma once

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

namespace IGSV {

  using EP_IC_K    = CGAL::Exact_predicates_inexact_constructions_kernel;
  using uVb        = CGAL::Triangulation_vertex_base_with_info_2<unsigned, EP_IC_K>;
  using uFb        = CGAL::Triangulation_face_base_with_info_2<unsigned, EP_IC_K>;
  using IndexedTds = CGAL::Triangulation_data_structure_2<uVb, uFb>;
  using IndexedDT  = CGAL::Delaunay_triangulation_2<EP_IC_K, IndexedTds>;

} // namespace IGSV
