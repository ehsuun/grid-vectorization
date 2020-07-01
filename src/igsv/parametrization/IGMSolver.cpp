// #include <CoMISo/NSolver/NPDerivativeChecker.hh>
#include <CoMISo/Utils/StopWatch.hh>
#include <igsv/common/CornerType.h>
#include <igsv/common/defs.h>
#include <igsv/common/print_timings.h>
#include <igsv/parametrization/IGMSolver.h>
#include <igsv/parametrization/IGM_Elements.h>

// ===========================================================================
// === CLASS : custom FEM problem with a constant Hessian
// ===========================================================================
class IGM_FEMProblem : public COMISO::FiniteElementProblem {
public:
  explicit IGM_FEMProblem(const unsigned int _n) : COMISO::FiniteElementProblem(_n) {}
  ~IGM_FEMProblem() {}
  bool constant_hessian() const override { return true; }
};

namespace IGSV {

  // ===========================================================================

  bool IGMSolver::solve(double& solve_time) {

    log_message((boost::format("w_snap: %.2f") % _wsnap).str());
    log_message((boost::format("w_reg : %.1e") % _wreg).str());

    log_message("setup ...");

    // -------------------------------------------------------------------------
    // Prepare constraints
    // -------------------------------------------------------------------------
    DEBUG_PRINT_INFO("prepare constraints");
    std::vector<VarTypePair> discrete_constraints;
    std::vector<COMISO::LinearConstraint> constraints;

    // -------------------------------------------------------------------------
    // Setup global indexing
    // -------------------------------------------------------------------------
    DEBUG_PRINT_INFO("compute global indexing");
    if (!global_indexing())
      return false;

    // -------------------------------------------------------------------------
    // Get variable counts
    // -------------------------------------------------------------------------
    DEBUG_PRINT_INFO("get variable counts");
    if (!get_variable_counts())
      return false;

    DEBUG_PRINT_INFO(
        "n_complete           = %zd\n"
        "n_real               = %zd\n"
        "n_integers           = %zd\n"
        "n_integers_transfn   = %zd\n"
        "n_integers_aux       = %zd\n",
        n_complete, n_real, n_integers, n_integers_transfn, n_integers_aux);

    // -------------------------------------------------------------------------
    // Setup finite elements
    // -------------------------------------------------------------------------
    DEBUG_PRINT_INFO("setup Poisson elements");
    COMISO::FiniteElementSet<PoissonElement> fe_poisson;
    if (!setup_poisson_elements(fe_poisson))
      return false;

    COMISO::FiniteElementSet<SnappingElement> fe_snapping;
    if (_wsnap > 0) {
      DEBUG_PRINT_INFO("setup snapping elements");
      if (!setup_snapping_elements(fe_snapping, discrete_constraints, constraints))
        return false;
    }

    // -------------------------------------------------------------------------
    // Setup transition functions
    // -------------------------------------------------------------------------
    DEBUG_PRINT_INFO("setup transition functions");
    if (!setup_transition_functions(discrete_constraints, constraints))
      return false;

    // -------------------------------------------------------------------------
    // Setup transition functions
    // -------------------------------------------------------------------------
    if (_use_sing_constraints) {
      DEBUG_PRINT_INFO("setup singularity constraints");
      if (!add_singularity_constraints(discrete_constraints))
        return false;
    }

    // -------------------------------------------------------------------------
    // Setup constraints to fix the translation DoFs (1 per connected component)
    // -------------------------------------------------------------------------
    // if (!setup_translation_constraints(constraints))
    //   return false;

    // -------------------------------------------------------------------------
    // Setup FEM problem
    // -------------------------------------------------------------------------
    DEBUG_PRINT_INFO("setup FEM, n=%zd", n_complete);
    IGM_FEMProblem fe_problem(n_complete);
    fe_problem.add_set(&fe_poisson);
    fe_problem.add_set(&fe_snapping);

    // -------------------------------------------------------------------------
    // (Optional) Check derivatives
    // -------------------------------------------------------------------------
    // std::cout << "  (optional for debugging) check derivatives of problem" << std::endl;
    // IGM_FEMProblem fe_problem__DEBUG(n_complete);
    // IGM_FEMProblem fe_problem__DEBUG(n_real);
    // fe_problem__DEBUG.add_set(&fe_poisson);
    // fe_problem__DEBUG.add_set(&fe_snapping);
    // COMISO::NPDerivativeChecker npd;
    // npd.check_all(&fe_problem__DEBUG);

    // -------------------------------------------------------------------------
    // Setup the MI solver
    // -------------------------------------------------------------------------
    // DEBUG_PRINT_INFO("setup solver");
    log_message("solve ...");

    COMISO::COMISOSolver comiso;

    comiso.solver().misolver().set_noise(0); // verbosity: 0 - quiet, 1 - more noise, 2 - even more, 100 - all noise
    comiso.solver().misolver().set_stats(false);
    comiso.solver().misolver().set_multiple_rounding(false);
    comiso.solver().misolver().set_inital_full(true);
    comiso.solver().misolver().set_final_full(true);
    comiso.solver().misolver().set_iter_full(true);
    comiso.solver().misolver().set_local_iters(20);
    comiso.solver().misolver().set_cg_iters(20);
    comiso.solver().misolver().set_cg_error(1e-8);

    std::vector<COMISO::NConstraintInterface*> constraint_pointers;
    constraint_pointers.reserve(constraints.size());
    for (int i = 0; i < constraints.size(); ++i)
      constraint_pointers.push_back(&(constraints[i]));

    COMISO::StopWatch sw;
    sw.start();
    comiso.solve(&fe_problem, constraint_pointers, discrete_constraints, _wreg, false, false);
    solve_time = sw.stop() / 1000.0;

    // -------------------------------------------------------------------------
    // Map coords to the matrix
    // -------------------------------------------------------------------------
    // DEBUG_PRINT_INFO("map coords");
    _UV.resize(n_uvertices + 1, 2); // the last row is for vertices outside the narrow band
    _UV.fill(0.5);
    for (int i = 0; i < n_uvertices; ++i)
      _UV.row(i) << fe_problem.x()[2 * i], fe_problem.x()[2 * i + 1];
    _UV.row(n_uvertices).fill(0.5);

    // -------------------------------------------------------------------------
    // Sanity checks
    // -------------------------------------------------------------------------
    // energies
    // std::cout << "========================================" << std::endl
    //           << "  Final energies:" << std::endl
    //           << "  e_poisson     = " << fe_poisson.eval_f(fe_problem.x().data()) << std::endl
    //           << "  e_snapping    = " << fe_snapping.eval_f(fe_problem.x().data()) //
    //           << "    (w = " << _wsnap << ")" << std::endl
    //           << "  --------------------------" << std::endl
    //           << "  e_Total       = " << fe_problem.eval_f(fe_problem.x().data()) << std::endl
    //           << "========================================" << std::endl;
    //
    // snapping integer variables
    // int n_aux_int_zero = 0;
    // for (int i = n_real + n_integers_transfn; i < n_complete; ++i)
    //   if (std::abs(fe_problem.x()[i]) < 1e-4)
    //     n_aux_int_zero++;
    // if (n_aux_int_zero > 0)
    //   std::cerr << "WARNING: " << n_aux_int_zero << " of " << n_integers_aux //
    //             << " aux variables are zero" << std::endl
    //             << "========================================" << std::endl;

    return true;
  }

  // ===========================================================================

  bool IGMSolver::global_indexing() {

    n_uvertices = 0;
    _FUV.resize(_F.rows(), 3);
    _FUV.fill(-1);

    _vertex_uvertices.clear();
    _vertex_uvertices.resize(_V.rows());

    int f0, k0, f1, k1, n_iter;
    const int MAX_ITER = 100;

    // loop over corners
    for (int f = 0; f < _F.rows(); ++f)
      if (_F_isNarrow(f))
        for (int k = 0; k < 3; ++k)
          if (_FUV(f, k) == -1) {

            // assign an index to the corner
            _FUV(f, k) = n_uvertices;
            _vertex_uvertices[_F(f, k)].push_back(n_uvertices);

            // -------------------------------------------------------------------------
            // vertex v : loop around adjacent edges (CW)
            // -------------------------------------------------------------------------
            f0     = f;
            k0     = k;
            n_iter = 0;
            while (n_iter++ < MAX_ITER /* safety check */) {
              // stop if there is a cut
              if (_E_customCut(f0, k0))
                break;

              // edge -> opposite -> next
              f1 = _TT(f0, k0);
              k1 = (_TTi(f0, k0) + 1) % 3;
              f0 = f1;
              k0 = k1;

              // stop if leaving the narrow band or the corner is already indexed
              if (!_F_isNarrow(f0) || _FUV(f0, k0) != -1)
                break;

              // otherwise, assign the same index
              _FUV(f0, k0) = n_uvertices;
            }

            // -------------------------------------------------------------------------
            // vertex v : loop around adjacent edges (CCW)
            // -------------------------------------------------------------------------
            f0     = f;
            k0     = k;
            n_iter = 0;
            while (n_iter++ < MAX_ITER /* safety check */) {

              // stop if there is a cut
              if (_E_customCut(f0, (k0 + 2) % 3))
                break;

              // edge -> previous -> opposite
              f1 = _TT(f0, (k0 + 2) % 3); // (x+2) mod 3 = (x-1) mod 3
              k1 = _TTi(f0, (k0 + 2) % 3);
              f0 = f1;
              k0 = k1;

              // stop if leaving the narrow band or the corner is already indexed
              if (!_F_isNarrow(f0) || _FUV(f0, k0) != -1)
                break;

              // otherwise, assign the same index
              _FUV(f0, k0) = n_uvertices;
            }

            // move on
            n_uvertices++;
          }

    // fill unused corners
    for (int f = 0; f < _F.rows(); ++f)
      for (int k = 0; k < 3; ++k)
        if (_FUV(f, k) == -1)
          _FUV(f, k) = n_uvertices;

    // for each u-vertex : store a list of corners
    _uvertex_corners.resize(n_uvertices);
    for (int f = 0; f < _F.rows(); ++f)
      for (int k = 0; k < 3; ++k)
        if (_FUV(f, k) < n_uvertices)
          _uvertex_corners[_FUV(f, k)].push_back(std::make_pair(f, k));

    // *** sanity check -- print ***
    // {
    //   int v = 0;
    //   for (const auto& list : _vertex_uvertices) {
    //     if (!list.empty()) {
    //       std::cout << "v#" << v << " : ";
    //       for (auto i : list)
    //         std::cout << i << " ";
    //       std::cout << std::endl;
    //     }
    //     v++;
    //   }
    // }

    return true;
  }

  // ===========================================================================

  bool IGMSolver::get_variable_counts() {

    _U_snapping_constraints.setZero(n_uvertices, 1);
    _U_snapping_weights.setZero(n_uvertices, 1);
    _FK_snapping_constraints.setZero(_F.rows(), 3);

    // number of integer variables
    // -- two per cut edge
    n_integers_transfn = 0;
    for (int f = 0; f < _F.rows(); ++f)
      for (int k = 0; k < 3; ++k)
        if (_E_isNarrow(f, k) & _E_customCut(f, k))
          n_integers_transfn++; // each edge is counted twice, so this gives two integers per edge, as required

    // -- one per black vertex with a label
    for (int f = 0; f < _F.rows(); ++f)
      if (_F_isNarrow(f) && _F_labels(f) != CORNER_TYPE_NO_INTEGERS)
        for (int k = 0; k < 3; ++k)
          if (_V_isBlack(_F(f, k))) {
            _U_snapping_constraints(_FUV(f, k)) |= _F_labels(f);
            _U_snapping_weights(_FUV(f, k)) = _V_snappingWeights(_F(f, k));

            _FK_snapping_constraints(f, k) |= _F_labels(f);
          }

    n_integers_aux = 0;
    if (_wsnap > 0) {
      for (int u = 0; u < n_uvertices; ++u) {
        if (_U_snapping_constraints(u) & CORNER_TYPE_INTEGER_U)
          n_integers_aux++;
        if (_U_snapping_constraints(u) & CORNER_TYPE_INTEGER_V)
          n_integers_aux++;
      }
    }

    n_integers = n_integers_transfn + n_integers_aux;

    // total number of variables
    n_real     = 2 * n_uvertices;
    n_complete = n_real + n_integers;

    return true;
  }

  // ===========================================================================

  bool IGMSolver::setup_poisson_elements(COMISO::FiniteElementSet<PoissonElement>& _fe_poisson) {
    PoissonElement::VecI pvi;
    PoissonElement::VecC pvc;
    Eigen::Vector3d p0, p1, p2, fu, fv;

    for (int f = 0; f < _F.rows(); ++f) {

      // skip if triangle not in the narrow band
      if (!_F_isNarrow(f))
        continue;

      // vertex indices
      const int v0 = _F(f, 0);
      const int v1 = _F(f, 1);
      const int v2 = _F(f, 2);

      // corner indices
      const int g0 = _FUV(f, 0);
      const int g1 = _FUV(f, 1);
      const int g2 = _FUV(f, 2);
      pvi << 2 * g0, 2 * g1, 2 * g2, 2 * g0 + 1, 2 * g1 + 1, 2 * g2 + 1;

      // triangle data
      p0 << _V.row(v0).transpose();
      p1 << _V.row(v1).transpose();
      p2 << _V.row(v2).transpose();
      fu << _XU.row(f).transpose();
      fv << _XV.row(f).transpose();
      PoissonElement::constants_from_triangle(p0, p1, p2, fu, fv, //
                                              _F_uvScale(f),      //
                                              _F_uvScale(f),      //
                                              _F_labels(f),       //
                                              _anisotropy_alpha,  //
                                              pvc);
      _fe_poisson.instances().add_element(pvi, pvc);

      // *** sanity check ***
      // {
      //   PoissonElement test_element;
      //   PoissonElement::VecV test_x;
      //   PoissonElement::VecC test_c;
      //   const Eigen::Vector3d test_fu(1., 0., 0.);
      //   const Eigen::Vector3d test_fv(0., 1., 0.);
      //   PoissonElement::constants_from_triangle(p0, p1, p2, test_fu, test_fv, 1., 1., 1., 0, test_c);
      //   test_x[0]           = p0[0]; // p0_u
      //   test_x[1]           = p1[0]; // p1_u
      //   test_x[2]           = p2[0]; // p2_u
      //   test_x[3]           = p0[1]; // p0_v
      //   test_x[4]           = p1[1]; // p1_v
      //   test_x[5]           = p2[1]; // p2_v
      //   const double energy = test_element.eval_f(test_x, test_c);
      //   if (energy > 1e-8)
      //     std::cerr << "WARNING : face # " << f << ", identity test should lead be zero energy : " << energy
      //               << std::endl;
      // }
    }
    return true;
  }

  // ===========================================================================

  bool IGMSolver::setup_snapping_elements(COMISO::FiniteElementSet<SnappingElement>& _fe_snapping, //
                                          std::vector<VarTypePair>& _discrete_constraints,         //
                                          std::vector<COMISO::LinearConstraint>& _constraints) {

    // i) ADD ELEMENTS for soft snapping
    //    while introducing auxiliary integer variables
    SnappingElement::VecI svi;
    SnappingElement::VecC svc;
    svc[0] = 0.0;

    const double wsnap_alpha_sqrt = _wsnap * std::sqrt(_anisotropy_alpha);

    // initialize the index of an aux integer variable
    int ivar_idx = n_real + n_integers_transfn;

    // // for each u-vertex, store the corresponding auxiliary integer variable (if any)
    Eigen::VectorXi aux_indices_U, aux_indices_V;
    aux_indices_U.resize(n_uvertices, 1);
    aux_indices_V.resize(n_uvertices, 1);
    aux_indices_U.fill(-1);
    aux_indices_V.fill(-1);

    for (int u = 0; u < n_uvertices; ++u) {

      svc[0] = wsnap_alpha_sqrt * _U_snapping_weights(u);

      if (_U_snapping_constraints(u) & CORNER_TYPE_INTEGER_U) {
        svi[0] = 2 * u + 1;
        svi[1] = ivar_idx++;
        _fe_snapping.instances().add_element(svi, svc);
        _discrete_constraints.emplace_back(svi[1], COMISO::Integer);
        aux_indices_U(u) = svi[1];
      }

      if (_U_snapping_constraints(u) & CORNER_TYPE_INTEGER_V) {
        svi[0] = 2 * u;
        svi[1] = ivar_idx++;
        _fe_snapping.instances().add_element(svi, svc);
        _discrete_constraints.emplace_back(svi[1], COMISO::Integer);
        aux_indices_V(u) = svi[1];
      }
    }

    // ii) ADD CONSTRAINTS between auxiliary integer variables
    COMISO::LinearConstraint::SVectorNC coeffs(n_complete);

#define SNAPPING_SAME_INTEGER_THRESHOLD 0.5

    const size_t old_n_constraints = _constraints.size();

    int n_threshold_exceeded = 0; // keep track of how many constraints were NOT added due to threshold

    _I_snapping_connections.resize(_F.rows() * 6, 4); // at most 6 connections per face
    int n_snap_connections = 0;

    // loop over triangles : add constraints between auxiliary variables
    for (int f = 0; f < _F.rows(); ++f)
      if (_F_isNarrow(f)) {

        // how many u-vertices in this triangle are snapped to integers?
        int n_constrained = 0;
        for (int k = 0; k < 3; ++k)
          if (_U_snapping_constraints(_FUV(f, k)) > 0)
            n_constrained++;

        // if at most 1 is snapped, skip -- no constraints can be added
        if (n_constrained < 2)
          continue;

        // get indices of u-vertices in this triangle
        const int u0 = _FUV(f, 0);
        const int u1 = _FUV(f, 1);
        const int u2 = _FUV(f, 2);

        // which u-vertices are snapped in the u dir?
        const bool u0_u = _U_snapping_constraints(u0) & CORNER_TYPE_INTEGER_U;
        const bool u1_u = _U_snapping_constraints(u1) & CORNER_TYPE_INTEGER_U;
        const bool u2_u = _U_snapping_constraints(u2) & CORNER_TYPE_INTEGER_U;

        // which u-vertices are snapped in the v dir?
        const bool u0_v = _U_snapping_constraints(u0) & CORNER_TYPE_INTEGER_V;
        const bool u1_v = _U_snapping_constraints(u1) & CORNER_TYPE_INTEGER_V;
        const bool u2_v = _U_snapping_constraints(u2) & CORNER_TYPE_INTEGER_V;

        // ---------------------------------------------------------------------
        // Transform the triangle to uv space
        // ---------------------------------------------------------------------
        using Matrix23d = Eigen::Matrix<double, 2, 3>;
        // triangle coords in the xy space, in cols
        Matrix23d T_xy;
        for (int k = 0; k < 3; ++k)
          T_xy.col(k) << _V(_F(f, k), 0), _V(_F(f, k), 1); // only store xy coords, z=0
        // transformation : canonical frame [e0,e1] to [Xu,Xv]
        Eigen::Matrix2d M;
        M.col(0) << _XU(f, 0), _XU(f, 1);
        M.col(1) << _XV(f, 0), _XV(f, 1);
        // linear transformation to the uv space : use inverse() directly, ok for a small matrix
        Matrix23d T_uv = M.inverse() * T_xy;

        // isotropic scaling
        T_uv.array() *= _F_uvScale(f);

        // ---------------------------------------------------------------------
        // Add constraints in the u dir
        // -- for each pair, check if they are close enough
        // -- u coord is in row 0 of T_uv
        // ---------------------------------------------------------------------

        if (u0_u & u1_u) {
          if (std::abs(T_uv(0, 0) - T_uv(0, 1)) < SNAPPING_SAME_INTEGER_THRESHOLD) {
            coeffs.setZero();
            coeffs.coeffRef(aux_indices_U(u0)) = +1;
            coeffs.coeffRef(aux_indices_U(u1)) = -1;
            _constraints.emplace_back(coeffs, 0., COMISO::NConstraintInterface::NC_EQUAL);

            _I_snapping_connections.row(n_snap_connections++) << f, 0, 1, 0;

          } else
            n_threshold_exceeded++;
        }

        if (u0_u & u2_u) {
          if (std::abs(T_uv(0, 0) - T_uv(0, 2)) < SNAPPING_SAME_INTEGER_THRESHOLD) {
            coeffs.setZero();
            coeffs.coeffRef(aux_indices_U(u0)) = +1;
            coeffs.coeffRef(aux_indices_U(u2)) = -1;
            _constraints.emplace_back(coeffs, 0., COMISO::NConstraintInterface::NC_EQUAL);

            _I_snapping_connections.row(n_snap_connections++) << f, 0, 2, 0;

          } else
            n_threshold_exceeded++;
        }

        if (u1_u & u2_u) {
          if (std::abs(T_uv(0, 1) - T_uv(0, 2)) < SNAPPING_SAME_INTEGER_THRESHOLD) {
            coeffs.setZero();
            coeffs.coeffRef(aux_indices_U(u1)) = +1;
            coeffs.coeffRef(aux_indices_U(u2)) = -1;
            _constraints.emplace_back(coeffs, 0., COMISO::NConstraintInterface::NC_EQUAL);

            _I_snapping_connections.row(n_snap_connections++) << f, 1, 2, 0;

          } else
            n_threshold_exceeded++;
        }

        // ---------------------------------------------------------------------
        // Add constraints in the v dir
        // -- for each pair, check if they are close enough
        // -- v coord is in row 1 of T_uv
        // ---------------------------------------------------------------------

        if (u0_v & u1_v) {
          if (std::abs(T_uv(1, 0) - T_uv(1, 1)) < SNAPPING_SAME_INTEGER_THRESHOLD) {
            coeffs.setZero();
            coeffs.coeffRef(aux_indices_V(u0)) = +1;
            coeffs.coeffRef(aux_indices_V(u1)) = -1;
            _constraints.emplace_back(coeffs, 0., COMISO::NConstraintInterface::NC_EQUAL);

            _I_snapping_connections.row(n_snap_connections++) << f, 0, 1, 1;

          } else
            n_threshold_exceeded++;
        }

        if (u0_v & u2_v) {
          if (std::abs(T_uv(1, 0) - T_uv(1, 2)) < SNAPPING_SAME_INTEGER_THRESHOLD) {
            coeffs.setZero();
            coeffs.coeffRef(aux_indices_V(u0)) = +1;
            coeffs.coeffRef(aux_indices_V(u2)) = -1;
            _constraints.emplace_back(coeffs, 0., COMISO::NConstraintInterface::NC_EQUAL);

            _I_snapping_connections.row(n_snap_connections++) << f, 0, 2, 1;

          } else
            n_threshold_exceeded++;
        }

        if (u1_v & u2_v) {
          if (std::abs(T_uv(1, 1) - T_uv(1, 2)) < SNAPPING_SAME_INTEGER_THRESHOLD) {
            coeffs.setZero();
            coeffs.coeffRef(aux_indices_V(u1)) = +1;
            coeffs.coeffRef(aux_indices_V(u2)) = -1;
            _constraints.emplace_back(coeffs, 0., COMISO::NConstraintInterface::NC_EQUAL);

            _I_snapping_connections.row(n_snap_connections++) << f, 1, 2, 1;

          } else
            n_threshold_exceeded++;
        }
      }

    _I_snapping_connections.conservativeResize(n_snap_connections, 4);

    // if (PRINT_DEBUGINFO) {
    //   std::cout << "added " << _constraints.size() - old_n_constraints << " constrained aux int pairs" << std::endl;
    //   std::cout << "threshold exceeded for " << n_threshold_exceeded << " pairs" << std::endl;
    // }

    return true;
  }

  // ===========================================================================

  bool IGMSolver::setup_transition_functions(std::vector<VarTypePair>& discrete_constraints,
                                             std::vector<COMISO::LinearConstraint>& constraints) {

    COMISO::LinearConstraint::SVectorNC coeffs(n_complete);

    int ivar_idx = n_real;

    // loop over edges
    Eigen::Matrix<bool, Eigen::Dynamic, 3> processed;
    processed.resize(_F.rows(), 3);
    processed.fill(false);
    for (int f0 = 0; f0 < _F.rows(); ++f0)
      for (int k0 = 0; k0 < 3; ++k0)
        if (_E_isNarrow(f0, k0) & _E_customCut(f0, k0) && !processed(f0, k0)) {

          // opposite edge
          const int f1 = _TT(f0, k0);
          const int k1 = _TTi(f0, k0);

          // skip if boundary
          if (f1 == -1)
            continue;

          processed(f0, k0) = true;
          processed(f1, k1) = true;

          // corner indices in f0
          const int p0 = _FUV(f0, k0);
          const int q0 = _FUV(f0, (k0 + 1) % 3); // next vertex
          // corner indices in f1
          const int q1 = _FUV(f1, k1);
          const int p1 = _FUV(f1, (k1 + 1) % 3);
          // make sure they are different!
          if (p0 == p1 && q0 == q1) {
            std::cerr << "WARNING in " << __FUNCTION__ << " : indexing error!\n";
            std::cerr << "    p0 = " << p0 << std::endl;
            std::cerr << "    p1 = " << p1 << std::endl;
            std::cerr << "    q0 = " << q0 << std::endl;
            std::cerr << "    q1 = " << q1 << std::endl;
            return false;
          }

          //// indices of real vars
          const int up0 = 2 * p0;
          const int up1 = 2 * p1;
          const int uq0 = 2 * q0;
          const int uq1 = 2 * q1;
          const int vp0 = 2 * p0 + 1;
          const int vp1 = 2 * p1 + 1;
          const int vq0 = 2 * q0 + 1;
          const int vq1 = 2 * q1 + 1;

          // indices of integer vars
          const int tu = ivar_idx++;
          const int tv = ivar_idx++;
          discrete_constraints.emplace_back(tu, COMISO::Integer);
          discrete_constraints.emplace_back(tv, COMISO::Integer);

          // rotation angle (alpha) from the period jump
          double ca(0); // cos(alpha)
          double sa(0); // sin(alpha)
          switch (_E_periodJumps(f0, k0)) {
          case 0:
            ca = 1;
            break;
          case 1:
            sa = 1;
            break;
          case 2:
            ca = -1;
            break;
          case 3:
            sa = -1;
            break;
          default:
            std::cerr << "WARNING : Unexpected period jump : f= " << f0 << ", k= " << k0
                      << ", pj= " << _E_periodJumps(f0, k0) << std::endl;
            return false;
          }

          // (1) : ca*up0 - sa*vp0 + tu = up1
          coeffs.setZero();
          coeffs.coeffRef(up0) = ca;
          coeffs.coeffRef(vp0) = -sa;
          coeffs.coeffRef(tu)  = 1.;
          coeffs.coeffRef(up1) = -1.;
          constraints.emplace_back(coeffs, 0., COMISO::NConstraintInterface::NC_EQUAL);

          // // (2) : ca*uq0 - sa*vq0 + tu = uq1
          coeffs.setZero();
          coeffs.coeffRef(uq0) = ca;
          coeffs.coeffRef(vq0) = -sa;
          coeffs.coeffRef(tu)  = 1.;
          coeffs.coeffRef(uq1) = -1.;
          constraints.emplace_back(coeffs, 0., COMISO::NConstraintInterface::NC_EQUAL);

          // // (3) : sa*up0 + ca*vp0 + tv = vp1
          coeffs.setZero();
          coeffs.coeffRef(up0) = sa;
          coeffs.coeffRef(vp0) = ca;
          coeffs.coeffRef(tv)  = 1.;
          coeffs.coeffRef(vp1) = -1.;
          constraints.emplace_back(coeffs, 0., COMISO::NConstraintInterface::NC_EQUAL);

          // // (4) : sa*uq0 + ca*vq0 + tv = vq1
          coeffs.setZero();
          coeffs.coeffRef(uq0) = sa;
          coeffs.coeffRef(vq0) = ca;
          coeffs.coeffRef(tv)  = 1.;
          coeffs.coeffRef(vq1) = -1.;
          constraints.emplace_back(coeffs, 0., COMISO::NConstraintInterface::NC_EQUAL);
        }

    return true;
  }

  // ===========================================================================

  bool IGMSolver::add_singularity_constraints(std::vector<VarTypePair>& discrete_constraints) {
    int n_added = 0;
    for (int v = 0; v < _V.rows(); ++v)
      if ((_V_singularityIdx(v) > 0) && !_V_isNarrowBoundary(v))
        for (int u : _vertex_uvertices[v]) {
          discrete_constraints.emplace_back(2 * u, COMISO::Integer);
          discrete_constraints.emplace_back(2 * u + 1, COMISO::Integer);
          n_added++;
        }

    if (PRINT_DEBUGINFO)
      std::cout << "added " << n_added << " singularity constraints" << std::endl;

    return true;
  }

  // ===========================================================================

  void IGMSolver::test_objective_function() {

    COMISO::FiniteElementSet<PoissonElement>::VecC c;
    COMISO::FiniteElementSet<PoissonElement>::VecV x;

    PoissonElement ffe;

    const Eigen::Vector3d p0(-1., 0., 0.);
    const Eigen::Vector3d p1(+1., 0., 0.);
    const Eigen::Vector3d p2(0., +2., 0.);

    const Eigen::Vector3d fu(+1., 0., 0.);
    const Eigen::Vector3d fv(0., +1., 0.);

    PoissonElement::constants_from_triangle(p0, p1, p2, fu, fv, 1., 1., 0, 1., c);

    // u coords
    x[0] = p0[0];
    x[1] = p1[0];
    x[2] = p2[0];
    // v coords
    x[3] = p0[1];
    x[4] = p1[1];
    x[5] = p2[1];

    double d = ffe.eval_f(x, c);

    if (PRINT_DEBUGINFO)
      std::cerr << "identity transformation should lead to zero energy: " << d << std::endl;

    const Eigen::Vector2d offset(-1, 3);
    // shift u coords
    x[0] += offset[0];
    x[1] += offset[0];
    x[2] += offset[0];
    // shift v coords
    x[3] += offset[1];
    x[4] += offset[1];
    x[5] += offset[1];

    d = ffe.eval_f(x, c);

    if (PRINT_DEBUGINFO)
      std::cerr << "translated identity transformation should lead to zero energy: " << d << std::endl;
  }

  // ===========================================================================

  //// UNUSED : replaced with global regularization
  // bool IGMSolver::setup_translation_constraints(VectorizationData& data, //
  //                                    std::vector<COMISO::LinearConstraint>& _constraints) {
  //
  //   COMISO::LinearConstraint::SVectorNC coeffs(n_complete);
  //
  //   // get connected components
  //   std::vector<std::vector<int>> conncomp;
  //   connected_components(data.dual_connections, conncomp, _F.rows());
  //   int cid = 0;
  //
  //   // add one constraint per (valid) connected component
  //   for (const auto& comp : conncomp) {
  //
  //     // move on if there is only a single face outside the narrow band
  //     if (comp.size() == 1)
  //       if (!_F_isNarrow(comp[0]))
  //         continue;
  //
  //     // fix a vertex in the first face
  //     const int f    = comp[0];
  //     const int gidx = _FUV(f, 0);
  //
  //     // fix the u coord
  //     coeffs.setZero();
  //     coeffs.coeffRef(2 * gidx) = 1.;
  //     _constraints.push_back(COMISO::LinearConstraint(coeffs, 0., COMISO::NConstraintInterface::NC_EQUAL));
  //
  //     // fix the v coord
  //     coeffs.setZero();
  //     coeffs.coeffRef(2 * gidx + 1) = 1.;
  //     _constraints.push_back(COMISO::LinearConstraint(coeffs, 0., COMISO::NConstraintInterface::NC_EQUAL));
  //
  //     std::cout << "Fix translational DoFs: add a constraint for vertex " << gidx << std::endl;
  //   }
  //
  //   return true;
  // }

} // namespace IGSV
