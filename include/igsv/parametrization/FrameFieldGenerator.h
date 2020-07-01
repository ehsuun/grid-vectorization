#pragma once

//== INCLUDES =================================================================

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <complex>
#include <vector>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== CLASSES ================================================================

  class FrameFieldGenerator {

  public:
    FrameFieldGenerator(const Eigen::MatrixXd& V,                 //
                        const Eigen::MatrixXi& F,                 //
                        const Eigen::VectorXcd& C_tangents,       //
                        const Eigen::VectorXd& C_tangentWeights,  //
                        const Eigen::VectorXi& C_nearestVertices, //
                        Eigen::VectorXd& X0_arg_vtx,              //
                        Eigen::VectorXd& X0_mag_vtx,              //
                        Eigen::VectorXd& X1_arg_vtx,              //
                        Eigen::VectorXd& X1_mag_vtx,              //
                        Eigen::VectorXd& X0_arg,                  //
                        Eigen::VectorXd& X0_mag,                  //
                        Eigen::VectorXd& X1_arg,                  //
                        Eigen::VectorXd& X1_mag)
        : _V(V),
          _F(F),
          _C_tangents(C_tangents),
          _C_tangentWeights(C_tangentWeights),
          _C_nearestVertices(C_nearestVertices),
          _X0_arg_vtx(X0_arg_vtx),
          _X0_mag_vtx(X0_mag_vtx),
          _X1_arg_vtx(X1_arg_vtx),
          _X1_mag_vtx(X1_mag_vtx),
          _X0_arg(X0_arg),
          _X0_mag(X0_mag),
          _X1_arg(X1_arg),
          _X1_mag(X1_mag) {
      init();
    }

    void init();
    bool compute();

    // clang-format off
    void set_wsmooth      (double w) { _wsmooth = w; }
    void set_wtangent     (double w) { _wtangent = w; }
    void set_wortho       (double w) { _wortho = w; }
    void set_weighted_reg (bool x)   { _use_weighted_reg = x; }
    // clang-format on

  private:
    static void extract_roots(const std::complex<double>& c0, //
                              const std::complex<double>& c2, //
                              double& X0_arg,                 //
                              double& X0_mag,                 //
                              double& X1_arg,                 //
                              double& X1_mag                  //
    );

    bool build_energy_smooth(Eigen::SparseMatrix<double>& A, //
                             Eigen::VectorXd& B              //
    );

    static bool build_energy_align(const Eigen::VectorXcd& tau,    //
                                   const Eigen::VectorXd& weights, //
                                   const Eigen::VectorXi& nearest, //
                                   int nv,                         //
                                   Eigen::SparseMatrix<double>& A, //
                                   Eigen::VectorXd& B,             //
                                   bool grad_equal_zero = true     //
    );

    bool solve(Eigen::VectorXd& X);
    bool solve(Eigen::SparseMatrix<double>& A_smooth,  //
              Eigen::SparseMatrix<double>& A_tangent, //
              Eigen::SparseMatrix<double>& A_regular, //
              Eigen::VectorXd& B_smooth,              //
              Eigen::VectorXd& B_tangent,             //
              Eigen::VectorXd& B_regular,             //
              Eigen::VectorXd& X_                     //
    );

    bool extract(const Eigen::VectorXd& X);

  private:
    // input : mesh
    const Eigen::MatrixXd& _V; // vertices
    const Eigen::MatrixXi& _F; // faces

    // input : constraints
    const Eigen::VectorXcd& _C_tangents;       // tangent constraints (Sobel)
    const Eigen::VectorXd& _C_tangentWeights;  // tangent weights (Sobel)
    const Eigen::VectorXi& _C_nearestVertices; // indices of nearest vertex in the mesh (per-constraint)

    std::vector<std::vector<int>> V_constraints; // for each vertex: index of constraints for which it is nearests

    // output : per-vertex frame field
    Eigen::VectorXd& _X0_arg_vtx;
    Eigen::VectorXd& _X0_mag_vtx;
    Eigen::VectorXd& _X1_arg_vtx;
    Eigen::VectorXd& _X1_mag_vtx;

    // output : per-face frame field (interpolated from vertices)
    Eigen::VectorXd& _X0_arg;
    Eigen::VectorXd& _X0_mag;
    Eigen::VectorXd& _X1_arg;
    Eigen::VectorXd& _X1_mag;

    // settings : weights
    double _wsmooth  = 50.0;
    double _wtangent = 1.0;
    double _wortho   = 0.1;

    // settings : use weighted regularization?
    bool _use_weighted_reg = false;

    // constant setting : min angle
    const double _minangle = 20.0;

  public:
    // stats to show
    // -- number of "stiffening" iterations
    int _n_iters = 0;
    // -- number of (vertex) frames with small angles
    int _n_degenerate = 0;
  };

}; // namespace IGSV
