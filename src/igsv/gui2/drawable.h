//
// Created by Tibor Stanko on 13/04/2020.
//

#ifndef IGSV_DRAWABLE_H
#define IGSV_DRAWABLE_H

#include <array>
#include <vector>
namespace IGSV {
  namespace drawable {

    struct abstract {
      abstract() {}
      abstract(int n_pts) : xy(n_pts), rgb(n_pts) {}

      bool enabled = true;
      std::vector<std::array<float, 2>> xy;
      std::vector<std::array<float, 3>> rgb;
    };

    struct points : public abstract {
      points() : abstract() {}
      points(int n_pts) : abstract(n_pts) {}
    };

    struct edges : public abstract {
      edges() : abstract() {}
      edges(int n_pts, int n_edges) : abstract(n_pts), idx(n_edges) {}
      std::vector<std::array<unsigned, 2>> idx;
    };

    struct triangles : public abstract {
      triangles() : abstract() {}
      triangles(int n_pts, int n_triangles) : abstract(n_pts), uv(n_pts), idx(n_triangles) {}
      std::vector<std::array<float, 2>> uv;
      std::vector<std::array<unsigned, 3>> idx;
    };

  } // namespace drawable
} // namespace IGSV

#endif // IGSV_DRAWABLE_H
