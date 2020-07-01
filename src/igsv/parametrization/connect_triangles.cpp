#include <Eigen/Core>
#include <opencv2/opencv.hpp>
#include <vector>

#include <igsv/parametrization/_cgal_delaunay_2.h> // using IndexedDT

#include <igsv/parametrization/connect_triangles.h>

namespace IGSV {

  const Eigen::RowVector3d barycentric_interpolation(const Eigen::MatrixXd& _V, //
                                                     const Eigen::MatrixXi& _F, //
                                                     const Eigen::Vector3d& bc, //
                                                     int f) {
    Eigen::RowVector3d x(0.0, 0.0, 0.0);
    for (int k = 0; k < 3; k++)
      x.array() += bc(k) * _V.row(_F(f, k)).array();
    return x;
  }

  // ===========================================================================

  void connect_triangles(const Eigen::MatrixXd& _V,  //
                         const Eigen::MatrixXi& _F,  //
                         const Eigen::MatrixXi& _TT, //
                         const IndexedDT& _tri,      //
                         const int f0,               //
                         const int f1,               //
                         std::vector<int>& faceList, //
                         std::vector<int>& edgeList) {
    // if no points are provided, set them as barycenters perturbed by a small
    // amount of noise the perturbation is done on the barycentric coords to make
    // sure we stay inside the triangle noise_scale : has to be strictly smaller
    // than 1/3
    const double noise_scale = 0.1; // 0 means no noise is added
    Eigen::Vector3d noise0   = noise_scale * Eigen::Vector3d::Random();
    Eigen::Vector3d noise1   = noise_scale * Eigen::Vector3d::Random();
    noise0.array() -= noise0.mean(); // center, so that Σ t_i = 1
    noise1.array() -= noise1.mean(); // center, so that Σ t_i = 1
    const Eigen::Vector3d bc0 = Eigen::Vector3d(1. / 3., 1. / 3., 1. / 3.) + noise0;
    const Eigen::Vector3d bc1 = Eigen::Vector3d(1. / 3., 1. / 3., 1. / 3.) + noise1;

    const Eigen::RowVector3d p0 = barycentric_interpolation(_V, _F, bc0, f0);
    const Eigen::RowVector3d p1 = barycentric_interpolation(_V, _F, bc1, f1);
    // call the method with all inputs
    connect_triangles(_TT, _tri, f0, f1, faceList, edgeList, p0(0), p0(1), p1(0), p1(1));
  }

  // ===========================================================================

  void connect_triangles(const Eigen::MatrixXi& _TT, //
                         const IndexedDT& _tri,      //
                         const int f0,               //
                         const int f1,               //
                         std::vector<int>& faceList, //
                         std::vector<int>& edgeList, //
                         const double p0_x,          //
                         const double p0_y,          //
                         const double p1_x,          //
                         const double p1_y,          //
                         const IndexedDT::Face_handle fh0) {
    faceList.clear();
    edgeList.clear();

    // if f0 and f1 are the same face, just return it
    if (f0 == f1) {
      faceList.push_back(f0);
      return;
    }

    // if f0 and f1 are adjacent, just return them + the index of the common edge
    for (int k0 = 0; k0 < 3; k0++)
      if (_TT(f0, k0) == f1) {
        faceList.push_back(f0);
        faceList.push_back(f1);
        edgeList.push_back(k0);
        return;
      }

    // otherwise, split at the midpoint
    const double p01_x = 0.5 * (p0_x + p1_x);
    const double p01_y = 0.5 * (p0_y + p1_y);

    // locate the midpoint in the Delaunay triangulation
    IndexedDT::Face_handle fh01 = _tri.locate(EP_IC_K::Point_2(p01_x, p01_y), fh0);
    // get the index of the face in F
    const int f01 = fh01->info();

    // divide and conquer :
    // - init two pairs of lists
    std::vector<int> faceList_L, faceList_R;
    std::vector<int> edgeList_L, edgeList_R;

    // - get the result for both sides
    connect_triangles(_TT, _tri, f0, f01, faceList_L, edgeList_L, p0_x, p0_y, p01_x, p01_y, fh01);
    connect_triangles(_TT, _tri, f01, f1, faceList_R, edgeList_R, p01_x, p01_y, p1_x, p1_y, fh01);

    // - join the two sides
    faceList.insert(faceList.end(), faceList_L.begin(), faceList_L.end());
    faceList.insert(faceList.end(), faceList_R.begin() + 1,
                    faceList_R.end()); // +1 to avoid repeating the middle face
    edgeList.insert(edgeList.end(), edgeList_L.begin(), edgeList_L.end());
    edgeList.insert(edgeList.end(), edgeList_R.begin(), edgeList_R.end());
  }

  // ===========================================================================

  // void connect_triangles(const int f0,               //
  //                        const int f1,               //
  //                        std::vector<int>& faceList, //
  //                        std::vector<int>& edgeList, //
  //                        const std::complex<double>& p0,         //
  //                        const std::complex<double>& p1) {
  //   Eigen::Vector3d p0_vec(p0.real(), p0.imag(), 0.0);
  //   Eigen::Vector3d p1_vec(p1.real(), p1.imag(), 0.0);
  //   connect_triangles(f0, f1, faceList, edgeList, p0_vec, p1_vec);
  // }

  // ===========================================================================

} // namespace IGSV
