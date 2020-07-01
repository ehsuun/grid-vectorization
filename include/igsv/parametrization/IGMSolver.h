//=============================================================================
//
//  CLASS IGMSolver
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <CoMISo/Config/config.hh>
#include <CoMISo/NSolver/COMISOSolver.hh>
#include <CoMISo/NSolver/FiniteElementProblem.hh>
#include <CoMISo/NSolver/LinearConstraint.hh>
#include <Eigen/Core>
#include <vector>

#include "../common/CornerType.h"
#include "IGM_Elements.h"

//== FORWARD DECLARATIONS =====================================================

//== NAMESPACES ===============================================================

namespace IGSV {

  //== TYPEDEFS ===============================================================

  using VarTypePair = std::pair<unsigned int, COMISO::VariableType>;
  using CornerPair  = std::pair<unsigned int, unsigned int>;
  using CornerList  = std::vector<CornerPair>;
  using UVertexList = std::vector<int>;

  //== CLASS DEFINITION =======================================================

  class IGMSolver {
  public:
    IGMSolver(const Eigen::MatrixXd& V,                                         //
              const Eigen::MatrixXd& XU,                                        //
              const Eigen::MatrixXd& XV,                                        //
              const Eigen::MatrixXi& F,                                         //
              const Eigen::MatrixXi& TT,                                        //
              const Eigen::MatrixXi& TTi,                                       //
              const std::vector<std::vector<int>>& VT,                          //
              const std::vector<std::vector<int>>& VTi,                         //
              const Eigen::Matrix<bool, Eigen::Dynamic, 1>& V_isBlack,          //
              const Eigen::Matrix<bool, Eigen::Dynamic, 1>& V_isNarrowBoundary, //
              const Eigen::VectorXd& V_snappingWeights,                         //
              const Eigen::VectorXi& V_singularityIdx,                          //
              const Eigen::MatrixXi& E_periodJumps,                             //
              const Eigen::Matrix<bool, Eigen::Dynamic, 3>& E_isNarrow,         //
              const Eigen::Matrix<bool, Eigen::Dynamic, 3>& E_isNarrowBoundary, //
              const Eigen::Matrix<bool, Eigen::Dynamic, 3>& E_customCut,        //
              const Eigen::Matrix<bool, Eigen::Dynamic, 1>& F_isNarrow,         //
              const Eigen::VectorXi& F_labels,                                  //
              const Eigen::VectorXd& F_uvScale,                                 //
              Eigen::MatrixXd& UV,                                              //
              Eigen::MatrixXi& FUV,                                             //
              Eigen::MatrixXi& FK_snapping_constraints,                         //
              Eigen::MatrixXi& I_snapping_connections)
        : _V(V), //
          _XU(XU),
          _XV(XV),
          _F(F),
          _TT(TT),
          _TTi(TTi),
          _VT(VT),
          _VTi(VTi),
          _V_isBlack(V_isBlack),
          _V_isNarrowBoundary(V_isNarrowBoundary),
          _V_snappingWeights(V_snappingWeights),
          _V_singularityIdx(V_singularityIdx),
          _E_periodJumps(E_periodJumps),
          _E_isNarrow(E_isNarrow),
          _E_isNarrowBoundary(E_isNarrowBoundary),
          _E_customCut(E_customCut),
          _F_isNarrow(F_isNarrow),
          _F_labels(F_labels),
          _F_uvScale(F_uvScale),
          _UV(UV),
          _FUV(FUV),
          _FK_snapping_constraints(FK_snapping_constraints),
          _I_snapping_connections(I_snapping_connections) {}

    bool solve(double& solve_time);
    void set_snapping_weight(double w) { _wsnap = w; }

  private:
    bool global_indexing();
    bool get_variable_counts();
    bool setup_poisson_elements(COMISO::FiniteElementSet<PoissonElement>&);
    bool setup_snapping_elements(COMISO::FiniteElementSet<SnappingElement>&, std::vector<VarTypePair>&,
                                 std::vector<COMISO::LinearConstraint>&);
    bool setup_transition_functions(std::vector<VarTypePair>&, std::vector<COMISO::LinearConstraint>&);
    bool add_singularity_constraints(std::vector<VarTypePair>&);

    void test_objective_function();

    // variable counts
    size_t n_uvertices        = 0;
    size_t n_real             = 0;
    size_t n_integers_transfn = 0;
    size_t n_integers_aux     = 0;
    size_t n_integers         = 0;
    size_t n_complete         = 0;

    // weights
    double _wsnap                    = 25.0;  // snapping
    const double _wreg               = 1e-6;  // regularization
    const double _anisotropy_alpha   = 1.0;   // alpha for anisotropic norm
    const bool _use_sing_constraints = false; // use singularity constraints?

    // input : mesh, frame field, adjacency
    const Eigen::MatrixXd& _V;                 // vertices
    const Eigen::MatrixXd& _XU;                // combed frame field
    const Eigen::MatrixXd& _XV;                // ...
    const Eigen::MatrixXi& _F;                 // faces
    const Eigen::MatrixXi& _TT;                // triangle-triangle adjacency
    const Eigen::MatrixXi& _TTi;               // ...
    const std::vector<std::vector<int>>& _VT;  // vertex-triangle  adjacency
    const std::vector<std::vector<int>>& _VTi; // ...

    // input : vertex properties
    const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _V_isNarrowBoundary; // bool : is vertex NB boundary?
    const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _V_isBlack;          // bool : is vertex black?
    const Eigen::VectorXd& _V_snappingWeights;                         // double : per-vertex snapping weights
    const Eigen::VectorXi& _V_singularityIdx;                          // int : per-vertex singularity index

    // input : edge properties
    const Eigen::MatrixXi& _E_periodJumps;                             // int : per-edge period jumps
    const Eigen::Matrix<bool, Eigen::Dynamic, 3>& _E_isNarrow;         // bool : is edge in the narrow band?
    const Eigen::Matrix<bool, Eigen::Dynamic, 3>& _E_isNarrowBoundary; // bool : is edge on the bd of narrow band?
    const Eigen::Matrix<bool, Eigen::Dynamic, 3>& _E_customCut;        // bool : is edge a cut edge?

    // input : face properties
    const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow; // bool : is face in the narrow band?
    const Eigen::VectorXi& _F_labels;                          // int : per-triangle tangent labels
    const Eigen::VectorXd& _F_uvScale;                         // double : per-face uv scale

    // output : parametrization
    Eigen::MatrixXd& _UV;
    Eigen::MatrixXi& _FUV;

    // for optimization
    std::vector<UVertexList> _vertex_uvertices; // list of lists of uvertices for each vertex
    std::vector<CornerList> _uvertex_corners;   // list of lists of corners for each uvertex

    Eigen::VectorXi _U_snapping_constraints;
    Eigen::VectorXd _U_snapping_weights;

    // for visualisation
    Eigen::MatrixXi& _FK_snapping_constraints;
    Eigen::MatrixXi& _I_snapping_connections;
    // col 0    : face index (0, ..., nf-1)
    // col 1, 2 : corner indices (0,1,2)
    // col 3    : type (u=0 or v=1)
  };

  //===========================================================================
}; // namespace IGSV
