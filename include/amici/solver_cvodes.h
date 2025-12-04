#ifndef AMICI_SOLVER_CVODES_h
#define AMICI_SOLVER_CVODES_h

#include "amici/defines.h"
#include "amici/solver.h"
#include "amici/vector.h"

namespace amici {
class ExpData;
class ReturnData;
class Model_ODE;
class CVodeSolver;
} // namespace amici

// for serialization friend in Solver
namespace boost {
namespace serialization {
template <class Archive>
void serialize(Archive& ar, amici::CVodeSolver& s, unsigned int version);
}
} // namespace boost

namespace amici {

/**
 * @brief The CVodeSolver class is a wrapper around the SUNDIALS CVODES solver.
 */

class CVodeSolver : public Solver {
  public:
    using Solver::Solver;

    ~CVodeSolver() override = default;

    /**
     * @brief Clone this instance
     *
     * @return The clone
     */
    Solver* clone() const override;

    std::string get_class_name() const override { return "CVodeSolver"; };

    void reinit(
        realtype t0, AmiVector const& yy0, AmiVector const& yp0
    ) const override;

    void sens_reinit(
        AmiVectorArray const& yyS0, AmiVectorArray const& ypS0
    ) const override;

    void sens_toggle_off() const override;

    void reinit_b(
        int which, realtype tB0, AmiVector const& yyB0, AmiVector const& ypB0
    ) const override;

    void reinit_quad_b(int which, AmiVector const& yQB0) const override;

    int solve(realtype tout, int itask) const override;

    int solve_f(realtype tout, int itask, int* ncheckPtr) const override;

    void solve_b(realtype tBout, int itaskB) const override;

    void get_dky(realtype t, int k) const override;

    void get_sens_dky(realtype t, int k) const override;

    void get_quad_dky_b(realtype t, int k, int which) const override;

    void get_dky_b(realtype t, int k, int which) const override;

    void get_root_info(int* rootsfound) const override;

    void set_stop_time(realtype tstop) const override;

    void turn_off_root_finding() const override;

    Model const* get_model() const override;

#if !defined(EXHALE_DOXYGEN_SHOULD_SKIP_THIS)
    using Solver::set_linear_solver;

    using Solver::set_linear_solver_b;
#endif
    void set_linear_solver() const override;

    void set_linear_solver_b(int which) const override;

    void set_non_linear_solver() const override;

    void set_non_linear_solver_sens() const override;

    void set_non_linear_solver_b(int which) const override;

  protected:
    void calc_ic(realtype tout1) const override;

    void calc_ic_b(int which, realtype tout1) const override;

    void get_b(int which) const override;

    void get_sens() const override;

    void get_quad_b(int which) const override;

    void get_quad(realtype& t) const override;

    void get_quad_dky(realtype t, int k) const override;

    void reinit_post_process_f(realtype tnext) const override;

    void reinit_post_process_b(realtype tnext) const override;

    /**
     * @brief Post-processing of the solver memory after a discontinuity
     *
     * @param cv_mem pointer to CVODES solver memory object
     * @param t pointer to integration time
     * @param yout  new state vector
     * @param tout  anticipated next integration timepoint.
     */
    void reInit_post_process(
        void* cv_mem, realtype* t, AmiVector* yout, realtype tout
    ) const;

    void allocate_solver() const override;

    void set_ss_tolerances(double rtol, double atol) const override;

    void set_sens_ss_tolerances(double rtol, double const* atol) const override;

    void set_sens_err_con(bool error_corr) const override;

    void set_quad_err_con_b(int which, bool flag) const override;

    void set_quad_err_con(bool flag) const override;

    void set_user_data() const override;

    void set_user_data_b(int which) const override;

    void set_max_num_steps(long int mxsteps) const override;

    void set_stab_lim_det(int stldet) const override;

    void set_stab_lim_det_b(int which, int stldet) const override;

    void set_id(Model const* model) const override;

    void set_suppress_alg(bool flag) const override;

    /**
     * @brief reset_state reset the CVODES solver to restart integration after a
     * rhs discontinuity.
     *
     * @param cv_mem pointer to CVODES solver memory object
     * @param y0 new state vector
     */
    void reset_state(void* cv_mem, const_N_Vector y0) const;

    void set_sens_params(
        realtype const* p, realtype const* pbar, int const* plist
    ) const override;

    void adj_init() const override;

    void quad_init(AmiVector const& xQ0) const override;

    void allocate_solver_b(int* which) const override;

    void set_ss_tolerances_b(
        int which, realtype relTolB, realtype absTolB
    ) const override;

    void quad_ss_tolerances_b(
        int which, realtype reltolQB, realtype abstolQB
    ) const override;

    void quad_ss_tolerances(realtype reltolQ, realtype abstolQ) const override;

    void set_max_num_steps_b(int which, long int mxstepsB) const override;

    void diag() const override;

    void diag_b(int which) const override;

    void get_num_steps(void const* ami_mem, long int* numsteps) const override;

    void get_num_rhs_evals(
        void const* ami_mem, long int* numrhsevals
    ) const override;

    void get_num_err_test_fails(
        void const* ami_mem, long int* numerrtestfails
    ) const override;

    void get_num_non_lin_solv_conv_fails(
        void const* ami_mem, long int* numnonlinsolvconvfails
    ) const override;

    void get_last_order(void const* ami_mem, int* order) const override;

    void* get_adj_b_mem(void* ami_mem, int which) const override;

    /**
     * @brief Serialize amici::CVodeSolver to boost archive
     *
     * @param ar Archive
     * @param s Solver instance to serialize
     */
    template <class Archive>
    friend void boost::serialization::
        serialize(Archive& ar, CVodeSolver& s, unsigned int /*version*/);

    /**
     * @brief Equality operator
     *
     * @param a
     * @param b
     * @return Whether a and b are equal
     */
    friend bool operator==(CVodeSolver const& a, CVodeSolver const& b);

    void
    init(realtype t0, AmiVector const& x0, AmiVector const& dx0) const override;

    void init_steady_state(
        realtype t0, AmiVector const& x0, AmiVector const& dx0
    ) const override;

    void sens_init_1(
        AmiVectorArray const& sx0, AmiVectorArray const& sdx0
    ) const override;

    void b_init(
        int which, realtype tf, AmiVector const& xB0, AmiVector const& dxB0
    ) const override;

    void qb_init(int which, AmiVector const& xQB0) const override;

    void root_init(int ne) const override;

    void set_dense_jac_fn() const override;

    void set_sparse_jac_fn() const override;

    void set_band_jac_fn() const override;

    void set_jac_times_vec_fn() const override;

    void set_dense_jac_fn_b(int which) const override;

    void set_sparse_jac_fn_b(int which) const override;

    void set_band_jac_fn_b(int which) const override;

    void set_jac_times_vec_fn_b(int which) const override;

    void set_sparse_jac_fn_ss() const override;

    void apply_max_nonlin_iters() const override;

    void apply_max_conv_fails() const override;

    void apply_constraints() const override;

    void apply_max_step_size() const override;
};

} // namespace amici

#endif /* AMICI_SOLVER_CVODES_h */
