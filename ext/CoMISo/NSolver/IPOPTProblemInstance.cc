//=============================================================================
//
//  CLASS IPOPTSolverLean - IMPLEMENTATION
//
//=============================================================================

//== INCLUDES =================================================================

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_IPOPT_AVAILABLE
//=============================================================================

#include "CoMISo/Utils/CoMISoError.hh"
#include <CoMISo/Utils/gmm.hh>

#include <Base/Debug/DebTime.hh>

#include "NProblemGmmInterface.hh"
#include "NProblemInterface.hh"
#include "NConstraintInterface.hh"
#include "BoundConstraint.hh"
#include "IPOPTCallbackParameters.hh"

#include "IPOPTProblemInstance.hh"

//== NAMESPACES ===============================================================

namespace COMISO {



//== IMPLEMENTATION PROBLEM INSTANCE===========================================

void
IPOPTProblemInstance::
split_constraints(const std::vector<NConstraintInterface*>& _constraints)
{
  DEB_enter_func;
  // split user-provided constraints into general-constraints and bound-constraints
  constraints_      .clear();       constraints_.reserve(_constraints.size());
  bound_constraints_.clear(); bound_constraints_.reserve(_constraints.size());

  for(unsigned int i=0; i<_constraints.size(); ++i)
  {
    BoundConstraint* bnd_ptr = dynamic_cast<BoundConstraint*>(_constraints[i]);

    if(bnd_ptr)
      bound_constraints_.push_back(bnd_ptr);
    else
      constraints_.push_back(_constraints[i]);
  }
}


//-----------------------------------------------------------------------------


void
IPOPTProblemInstance::
analyze_special_properties(const NProblemInterface* _problem, const std::vector<NConstraintInterface*>& _constraints)
{
  hessian_constant_ = true;
  jac_c_constant_   = true;
  jac_d_constant_   = true;

  if(!_problem->constant_hessian())
    hessian_constant_ = false;

  for(unsigned int i=0; i<_constraints.size(); ++i)
  {
    if(!_constraints[i]->constant_hessian())
      hessian_constant_ = false;

    if(!_constraints[i]->constant_gradient())
    {
      if(_constraints[i]->constraint_type() == NConstraintInterface::NC_EQUAL)
        jac_c_constant_ = false;
      else
        jac_d_constant_ = false;
    }

    // nothing else to check?
    if(!hessian_constant_ && !jac_c_constant_ && !jac_d_constant_)
      break;
  }

  //hessian of Lagrangian is only constant, if all hessians of the constraints are zero (due to lambda multipliers)
  if(!jac_c_constant_ || !jac_d_constant_)
    hessian_constant_ = false;
}


//-----------------------------------------------------------------------------


bool IPOPTProblemInstance::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  DEB_enter_func;
  // number of variables
  n = static_cast<Index>(problem_->n_unknowns());

  // number of constraints
  m = static_cast<Index>(constraints_.size());

  // get non-zeros of hessian of lagrangian and jacobi of constraints
  nnz_jac_g = 0;
  nnz_h_lag = 0;

  // get nonzero structure
  std::vector<double> x(n);
  problem_->initial_x(P(x));

  // nonzeros in the jacobian of C_ and the hessian of the lagrangian
  SMatrixNP HP;
  SVectorNC g;
  SMatrixNC H;

  if (!hessian_approximation_)
  {
    problem_->eval_hessian(P(x), HP);

    // get nonzero structure of hessian of problem
    for(int i=0; i<HP.outerSize(); ++i)
      for (SMatrixNP::InnerIterator it(HP,i); it; ++it)
        if(it.row() >= it.col())
          ++nnz_h_lag;
  }

  // get nonzero structure of constraints
  for( int i=0; i<m; ++i)
  {
    constraints_[i]->eval_gradient(P(x),g);

    nnz_jac_g += static_cast<Index>(g.nonZeros());

    if(!hessian_approximation_)
    {
      // count lower triangular elements
      constraints_[i]->eval_hessian (P(x),H);

      SMatrixNC::iterator m_it = H.begin();
      for(; m_it != H.end(); ++m_it)
        if( m_it.row() >= m_it.col())
          ++nnz_h_lag;
    }
  }

  // We use the standard fortran index style for row/col entries
  index_style = C_STYLE;

  return true;
}


//-----------------------------------------------------------------------------


bool IPOPTProblemInstance::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
  DEB_enter_func;
  // check dimensions
  DEB_warning_if(( n != (Index)problem_->n_unknowns() ), 1,
    "IPOPT #unknowns != n " << n << problem_->n_unknowns() );
  DEB_warning_if(( m != (Index)constraints_.size() ), 1,
    "Warning: IPOPT #constraints != m " << m << constraints_.size() );


  // first clear all variable bounds
  for( int i=0; i<n; ++i)
  {
    // x_l[i] = Ipopt::nlp_lower_bound_inf;
    // x_u[i] = Ipopt::nlp_upper_bound_inf;

    x_l[i] = -1.0e19;
    x_u[i] =  1.0e19;
  }

  // iterate over bound constraints and set them
  for(unsigned int i=0; i<bound_constraints_.size(); ++i)
  {
    if((Index)(bound_constraints_[i]->idx()) < n)
    {
      switch(bound_constraints_[i]->constraint_type())
      {
      case NConstraintInterface::NC_LESS_EQUAL:
      {
        x_u[bound_constraints_[i]->idx()] = bound_constraints_[i]->bound();
      }break;

      case NConstraintInterface::NC_GREATER_EQUAL:
      {
        x_l[bound_constraints_[i]->idx()] = bound_constraints_[i]->bound();
      }break;

      case NConstraintInterface::NC_EQUAL:
      {
        x_l[bound_constraints_[i]->idx()] = bound_constraints_[i]->bound();
        x_u[bound_constraints_[i]->idx()] = bound_constraints_[i]->bound();
      }break;
      }
    }
    else
      DEB_warning(2, "invalid bound constraint in IPOPTSolverLean!!!")
  }

  // set bounds for constraints
  for( int i=0; i<m; ++i)
  {
    // enum ConstraintType {NC_EQUAL, NC_LESS_EQUAL, NC_GREATER_EQUAL};
    switch(constraints_[i]->constraint_type())
    {
      case NConstraintInterface::NC_EQUAL         : g_u[i] = 0.0   ; g_l[i] =  0.0   ; break;
      case NConstraintInterface::NC_LESS_EQUAL    : g_u[i] = 0.0   ; g_l[i] = -1.0e19; break;
      case NConstraintInterface::NC_GREATER_EQUAL : g_u[i] = 1.0e19; g_l[i] =  0.0   ; break;
      default                                     : g_u[i] = 1.0e19; g_l[i] = -1.0e19; break;
    }
  }

  return true;
}


//-----------------------------------------------------------------------------


bool IPOPTProblemInstance::get_starting_point(Index n, bool init_x, Number* x,
                               bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,
                               Number* lambda)
{
  DEB_enter_func;
  // get initial value of problem instance
  problem_->initial_x(x);

  return true;
}


//-----------------------------------------------------------------------------


bool IPOPTProblemInstance::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  DEB_enter_func;
  // return the value of the objective function
  obj_value = problem_->eval_f(x);
  return true;
}


//-----------------------------------------------------------------------------


bool IPOPTProblemInstance::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  DEB_enter_func;
  problem_->eval_gradient(x, grad_f);

  return true;
}


//-----------------------------------------------------------------------------


bool IPOPTProblemInstance::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  DEB_enter_func;
  // evaluate all constraint functions
  for( int i=0; i<m; ++i)
    g[i] = constraints_[i]->eval_constraint(x);

  return true;
}


//-----------------------------------------------------------------------------


bool IPOPTProblemInstance::eval_jac_g(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values)
{
  DEB_enter_func;
  if (values == NULL)
  {
    // get x for evaluation (arbitrary position should be ok)
    std::vector<double> x_rnd(problem_->n_unknowns(), 0.0);

    int gi = 0;
    SVectorNC g;
    for( int i=0; i<m; ++i)
    {
      constraints_[i]->eval_gradient(&(x_rnd[0]), g);
      SVectorNC::InnerIterator v_it(g);
      for( ; v_it; ++v_it)
      {
        iRow[gi] = i;
        jCol[gi] = v_it.index();
        ++gi;
      }
    }
  }
  else
  {
    // return the values of the jacobian of the constraints

    // return the structure of the jacobian of the constraints
    // global index
    int gi = 0;
    SVectorNC g;

    for( int i=0; i<m; ++i)
    {
      constraints_[i]->eval_gradient(x, g);

      SVectorNC::InnerIterator v_it(g);

      for( ; v_it; ++v_it)
      {
        values[gi] = v_it.value();
        ++gi;
      }
    }

    DEB_warning_if((gi != nele_jac), 1,
      "number of non-zeros in Jacobian of C is incorrect: "
                << gi << " vs " << nele_jac)
  }

  return true;
}


//-----------------------------------------------------------------------------


bool IPOPTProblemInstance::eval_h(Index n, const Number* x, bool new_x,
                   Number obj_factor, Index m, const Number* lambda,
                   bool new_lambda, Index nele_hess, Index* iRow,
                   Index* jCol, Number* values)
{
  DEB_enter_func;
  if (values == NULL)
  {
    // return structure

    // get x for evaluation (arbitrary position should be ok)
    std::vector<double> x_rnd(problem_->n_unknowns(), 0.0);

     // global index
     int gi = 0;
     // get hessian of problem
     SMatrixNP HP;
     problem_->eval_hessian(&(x_rnd[0]), HP);

     for(int i=0; i<HP.outerSize(); ++i)
       for (SMatrixNP::InnerIterator it(HP,i); it; ++it)
       {
         // store lower triangular part only
         if(it.row() >= it.col())
         {
           //         it.value();
           iRow[gi] = static_cast<Index>(it.row());
           jCol[gi] = static_cast<Index>(it.col());
           ++gi;
         }
       }

    // Hessians of Constraints
    for(unsigned int j=0; j<constraints_.size(); ++j)
    {
      SMatrixNC H;
      constraints_[j]->eval_hessian(&(x_rnd[0]), H);

      SMatrixNC::iterator m_it  = H.begin();
      SMatrixNC::iterator m_end = H.end();

      for(; m_it != m_end; ++m_it)
      {
        // store lower triangular part only
        if( m_it.row() >= m_it.col())
        {
          iRow[gi] = m_it.row();
          jCol[gi] = m_it.col();
          ++gi;
        }
      }
    }

    // error check
    DEB_warning_if(( gi != nele_hess), 1,
      "number of non-zeros in Hessian of Lagrangian is incorrect while indexing: "
                << gi << " vs " << nele_hess )
  }
  else
  {
    // return values.

    // global index
    int gi = 0;
    // get hessian of problem
    SMatrixNP HP;
    problem_->eval_hessian(x, HP);

    for(int i=0; i<HP.outerSize(); ++i)
      for (SMatrixNP::InnerIterator it(HP,i); it; ++it)
      {
        // store lower triangular part only
        if(it.row() >= it.col())
        {
          values[gi] = obj_factor*it.value();
          ++gi;
        }
      }

    // Hessians of Constraints
    for(unsigned int j=0; j<constraints_.size(); ++j)
    {
      SMatrixNC H;
      constraints_[j]->eval_hessian(x, H);

      SMatrixNC::iterator m_it  = H.begin();
      SMatrixNC::iterator m_end = H.end();

      for(; m_it != m_end; ++m_it)
      {
        // store lower triangular part only
        if( m_it.row() >= m_it.col())
        {
          values[gi] = lambda[j]*(*m_it);
          ++gi;
        }
      }
    }

    // error check
    DEB_warning_if(( gi != nele_hess), 1,
      "number of non-zeros in Hessian of Lagrangian is incorrect2: "
                << gi << " vs " << nele_hess )
  }
  return true;
}


//-----------------------------------------------------------------------------


//inline double _QNT(double x)
//{
//	return double(float(x));
//}

//double _QNT(const double x)
//    {
//    // clear the 12 least significant mantissa bits to reduce noise
//    //const double fact = pow(2., 41);
//
//	const double fact = pow(2., 37);
//    int i;
//    double m = frexp(x, &i);
//    m *= fact;
//    int sgn_x = m < 0 ? -1 : 1;
//    m = sgn_x * floor(fabs(m));
//    m /= fact;
//    double xq = ldexp(m, i);
//    return xq;
//    }

double _QNT(const double x) { return x; }

void IPOPTProblemInstance::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* lambda,
                              Number obj_value,
                              const IpoptData* ip_data,
                              IpoptCalculatedQuantities* ip_cq)
{
  DEB_enter_func;

  // problem knows what to do
  problem_->store_result(x);

  if(store_solution_)
  {
    x_.resize(n);
    for( Index i=0; i<n; ++i)
      x_[i] = x[i];
  }


//  DEB_out(1, "Quantanizing the IPOPT solution\n");
//  std::vector<double> x_qnt(n);
//    for( Index i=0; i<n; ++i)
//      x_qnt[i] = _QNT(x[i]);
//
//
//
//  // problem knows what to do
//  problem_->store_result(&x_qnt[0]);
//
//  if(store_solution_)
//  {
//    x_.resize(n);
//    for( Index i=0; i<n; ++i)
//      x_[i] = x_qnt[i];
//  }
}


//-----------------------------------------------------------------------------
void
IPOPTProblemInstance::
set_callback_function
(std::function<bool(const IPOPTCallbackParameters &)> func)
{
  intermediate_callback_ = func;
}

bool IPOPTProblemInstance::intermediate_callback(
  Ipopt::AlgorithmMode mode,
  Index iter, Number obj_value,
  Number inf_pr, Number inf_du,
  Number mu, Number d_norm,
  Number regularization_size,
  Number alpha_du, Number alpha_pr,
  Index ls_trials,
  const IpoptData* ip_data,
  IpoptCalculatedQuantities* ip_cq
)
{
  PROGRESS_TICK;

  if(intermediate_callback_) {
    IPOPTCallbackParameters callbackParameters {
      mode,
      iter, obj_value,
      inf_pr, inf_du,
      mu, d_norm,
      regularization_size,
      alpha_du, alpha_pr,
      ls_trials,
      ip_data,
      ip_cq
    };

    return intermediate_callback_(callbackParameters);
  }

  return true;
}


//-----------------------------------------------------------------------------


bool IPOPTProblemInstance::hessian_constant() const
{
  return hessian_constant_;
}


//-----------------------------------------------------------------------------


bool IPOPTProblemInstance::jac_c_constant() const
{
  return jac_c_constant_;
}


//-----------------------------------------------------------------------------


bool IPOPTProblemInstance::jac_d_constant() const
{
  return jac_d_constant_;
}


//== IMPLEMENTATION PROBLEM INSTANCE==========================================================


bool IPOPTProblemInstanceGmm::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  DEB_enter_func;
  // number of variables
  n = problem_->n_unknowns();

  // number of constraints
  m = static_cast<Index>(constraints_.size());

  // get nonzero structure
  std::vector<double> x(n);
  problem_->initial_x(&(x[0]));
  // ToDo: perturb x

  // nonzeros in the jacobian of C_ and the hessian of the lagrangian
  SMatrixNP HP;
  SVectorNC g;
  SMatrixNC H;
  problem_->eval_hessian(&(x[0]), HP);
  nnz_jac_g = 0;
  nnz_h_lag = 0;

  // clear old data
  jac_g_iRow_.clear();
  jac_g_jCol_.clear();
  h_lag_iRow_.clear();
  h_lag_jCol_.clear();

  // get non-zero structure of initial hessian
  // iterate over rows
  for( int i=0; i<n; ++i)
  {
    SVectorNP& ri = HP.row(i);

    SVectorNP_citer v_it  = gmm::vect_const_begin(ri);
    SVectorNP_citer v_end = gmm::vect_const_end  (ri);

    for(; v_it != v_end; ++v_it)
    {
      // store lower triangular part only
      if( i >= (int)v_it.index())
      {
        h_lag_iRow_.push_back(i);
        h_lag_jCol_.push_back(static_cast<int>(v_it.index()));
        ++nnz_h_lag;
      }
    }
  }


  // get nonzero structure of constraints
  for( int i=0; i<m; ++i)
  {
    constraints_[i]->eval_gradient(&(x[0]),g);
    constraints_[i]->eval_hessian (&(x[0]),H);

    // iterate over sparse vector
    SVectorNC::InnerIterator v_it(g);
    for(; v_it; ++v_it)
    {
      jac_g_iRow_.push_back(i);
      jac_g_jCol_.push_back(v_it.index());
      ++nnz_jac_g;
    }

    // iterate over superSparseMatrix
    SMatrixNC::iterator m_it  = H.begin();
    SMatrixNC::iterator m_end = H.end();
    for(; m_it != m_end; ++m_it)
      if( m_it.row() >= m_it.col())
      {
        h_lag_iRow_.push_back(m_it.row());
        h_lag_jCol_.push_back(m_it.col());
        ++nnz_h_lag;
      }
  }

  // store for error checking...
  nnz_jac_g_ = nnz_jac_g;
  nnz_h_lag_ = nnz_h_lag;

  // We use the standard fortran index style for row/col entries
  index_style = C_STYLE;

  return true;
}


//-----------------------------------------------------------------------------


bool IPOPTProblemInstanceGmm::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
  DEB_enter_func;
  // first clear all variable bounds
  for( int i=0; i<n; ++i)
  {
    // x_l[i] = Ipopt::nlp_lower_bound_inf;
    // x_u[i] = Ipopt::nlp_upper_bound_inf;

    x_l[i] = -1.0e19;
    x_u[i] =  1.0e19;
  }

  // set bounds for constraints
  for( int i=0; i<m; ++i)
  {
    // enum ConstraintType {NC_EQUAL, NC_LESS_EQUAL, NC_GREATER_EQUAL};
    switch(constraints_[i]->constraint_type())
    {
      case NConstraintInterface::NC_EQUAL         : g_u[i] = 0.0   ; g_l[i] =  0.0   ; break;
      case NConstraintInterface::NC_LESS_EQUAL    : g_u[i] = 0.0   ; g_l[i] = -1.0e19; break;
      case NConstraintInterface::NC_GREATER_EQUAL : g_u[i] = 1.0e19; g_l[i] =  0.0   ; break;
      default                                     : g_u[i] = 1.0e19; g_l[i] = -1.0e19; break;
    }
  }

  return true;
}


//-----------------------------------------------------------------------------


bool IPOPTProblemInstanceGmm::get_starting_point(Index n, bool init_x, Number* x,
                               bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,
                               Number* lambda)
{
  DEB_enter_func;
  // get initial value of problem instance
  problem_->initial_x(x);

  return true;
}


//-----------------------------------------------------------------------------


bool IPOPTProblemInstanceGmm::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  DEB_enter_func;
  // return the value of the objective function
  obj_value = problem_->eval_f(x);
  return true;
}


//-----------------------------------------------------------------------------


bool IPOPTProblemInstanceGmm::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  DEB_enter_func;
  problem_->eval_gradient(x, grad_f);

  return true;
}


//-----------------------------------------------------------------------------


bool IPOPTProblemInstanceGmm::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  DEB_enter_func;
  // evaluate all constraint functions
  for( int i=0; i<m; ++i)
    g[i] = constraints_[i]->eval_constraint(x);

  return true;
}


//-----------------------------------------------------------------------------


bool IPOPTProblemInstanceGmm::eval_jac_g(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values)
{
  DEB_enter_func;
  if (values == NULL)
  {
    // return the (cached) structure of the jacobian of the constraints
    gmm::copy(jac_g_iRow_, VectorPTi(iRow, jac_g_iRow_.size()));
    gmm::copy(jac_g_jCol_, VectorPTi(jCol, jac_g_jCol_.size()));
  }
  else
  {
    // return the values of the jacobian of the constraints

    // return the structure of the jacobian of the constraints
    // global index
    int gi = 0;
    SVectorNC g;

    for( int i=0; i<m; ++i)
    {
      constraints_[i]->eval_gradient(x, g);

      SVectorNC::InnerIterator v_it(g);

      for( ; v_it; ++v_it)
      {
        if(gi < nele_jac)
          values[gi] = v_it.value();
        ++gi;
      }
    }

    DEB_warning_if(( gi != nele_jac), 1,
      "number of non-zeros in Jacobian of C is incorrect: "
       << gi << " vs " << nele_jac)
  }

  return true;
}


//-----------------------------------------------------------------------------


bool IPOPTProblemInstanceGmm::eval_h(Index n, const Number* x, bool new_x,
                   Number obj_factor, Index m, const Number* lambda,
                   bool new_lambda, Index nele_hess, Index* iRow,
                   Index* jCol, Number* values)
{
  DEB_enter_func;
  if (values == NULL)
  {
    // return the (cached) structure of the hessian
    gmm::copy(h_lag_iRow_, VectorPTi(iRow, h_lag_iRow_.size()));
    gmm::copy(h_lag_jCol_, VectorPTi(jCol, h_lag_jCol_.size()));
  }
  else
  {
    // return values.

    // global index
    int gi = 0;

    // get hessian of problem
    problem_->eval_hessian(x, HP_);

    for( int i=0; i<n; ++i)
    {
      SVectorNP& ri = HP_.row(i);

      SVectorNP_citer v_it  = gmm::vect_const_begin(ri);
      SVectorNP_citer v_end = gmm::vect_const_end  (ri);

      for(; v_it != v_end; ++v_it)
      {
        // store lower triangular part only
        if( i >= (int)v_it.index())
        {
          if( gi < nele_hess)
            values[gi] = obj_factor*(*v_it);
          ++gi;
        }
      }
    }

    // Hessians of Constraints
    for(unsigned int j=0; j<constraints_.size(); ++j)
    {
      SMatrixNC H;
      constraints_[j]->eval_hessian(x, H);

      SMatrixNC::iterator m_it  = H.begin();
      SMatrixNC::iterator m_end = H.end();

      for(; m_it != m_end; ++m_it)
      {
        // store lower triangular part only
        if( m_it.row() >= m_it.col())
        {
          if( gi < nele_hess)
            values[gi] = lambda[j]*(*m_it);
          ++gi;
        }
      }
    }

    // error check
    DEB_warning_if(( gi != nele_hess), 1,
      "number of non-zeros in Hessian of Lagrangian is incorrect: "
        << gi << " vs " << nele_hess);
  }
  return true;
}


//-----------------------------------------------------------------------------


void IPOPTProblemInstanceGmm::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* lambda,
                              Number obj_value,
                              const IpoptData* ip_data,
                              IpoptCalculatedQuantities* ip_cq)
{
  DEB_enter_func;
  // problem knows what to do
  problem_->store_result(x);
}

void
IPOPTProblemInstanceGmm::
set_callback_function
(std::function<bool(const IPOPTCallbackParameters &)> func)
{
  intermediate_callback_ = func;
}

bool IPOPTProblemInstanceGmm::intermediate_callback(
  Ipopt::AlgorithmMode mode,
  Index iter, Number obj_value,
  Number inf_pr, Number inf_du,
  Number mu, Number d_norm,
  Number regularization_size,
  Number alpha_du, Number alpha_pr,
  Index ls_trials,
  const IpoptData* ip_data,
  IpoptCalculatedQuantities* ip_cq
)
{
  PROGRESS_TICK;
  if(intermediate_callback_) {
    IPOPTCallbackParameters callbackParameters {
      mode,
      iter, obj_value,
      inf_pr, inf_du,
      mu, d_norm,
      regularization_size,
      alpha_du, alpha_pr,
      ls_trials,
      ip_data,
      ip_cq
    };

    return intermediate_callback_(callbackParameters);
  }

  return true;
}


//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // COMISO_IPOPT_AVAILABLE
//=============================================================================
