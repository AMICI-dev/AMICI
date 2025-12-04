#ifndef AMICI_SOLVER_IDAS_h
#define AMICI_SOLVER_IDAS_h

#include "amici/solver.h"

namespace amici {
class ExpData;
class ReturnData;
class Model_DAE;
class IDASolver;
} // namespace amici

// for serialization friend in Solver
namespace boost {
namespace serialization {
template <class Archive>
void serialize(Archive& ar, amici::IDASolver& s, unsigned int version);
}
} // namespace boost

namespace amici {

/**
 * @brief The IDASolver class is a wrapper around the SUNDIALS IDAS solver.
 */
class IDASolver : public Solver {
  public:
    using Solver::Solver;

    ~IDASolver() override = default;

    /**
     * @brief Clone this instance
     *
     * @return The clone
     */
    Solver* clone() const override;

    std::string get_class_name() const override { return "IDASolver"; };

    void reinit_post_process_f(realtype tnext) const override;

    void reinit_post_process_b(realtype tnext) const override;

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

    void quad_ss_tolerances_b(
        int which, realtype reltolQB, realtype abstolQB
    ) const override;

    void quad_ss_tolerances(realtype reltolQ, realtype abstolQ) const override;

    int solve(realtype tout, int itask) const override;

    int solve_f(realtype tout, int itask, int* ncheckPtr) const override;

    void solve_b(realtype tBout, int itaskB) const override;

    void get_root_info(int* rootsfound) const override;

    void get_dky(realtype t, int k) const override;

    void get_sens() const override;

    void get_sens_dky(realtype t, int k) const override;

    void get_b(int which) const override;

    void get_dky_b(realtype t, int k, int which) const override;

    void get_quad_b(int which) const override;

    void get_quad_dky_b(realtype t, int k, int which) const override;

    void get_quad(realtype& t) const override;

    void get_quad_dky(realtype t, int k) const override;

    void calc_ic(realtype tout1) const override;

    void calc_ic_b(int which, realtype tout1) const override;

    void set_stop_time(realtype tstop) const override;

    void turn_off_root_finding() const override;

    Model const* get_model() const override;

    void set_linear_solver() const override;

    void set_linear_solver_b(int which) const override;

    void set_non_linear_solver() const override;

    void set_non_linear_solver_sens() const override;

    void set_non_linear_solver_b(int which) const override;

  protected:
    /**
     * @brief Postprocessing of the solver memory after a discontinuity
     *
     * @param ida_mem pointer to IDAS solver memory object
     * @param t pointer to integration time
     * @param yout new state vector
     * @param ypout new state derivative vector
     * @param tout anticipated next integration timepoint.
     */
    void reinit_post_process(
        void* ida_mem, realtype* t, AmiVector* yout, AmiVector* ypout,
        realtype tout
    ) const;

    void allocate_solver() const override;

    void set_ss_tolerances(realtype rtol, realtype atol) const override;

    void
    set_sens_ss_tolerances(realtype rtol, realtype const* atol) const override;

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
     * @brief reset_state reset the IDAS solver to restart integration after a
     * rhs discontinuity.
     *
     * @param ida_mem pointer to IDAS solver memory object
     * @param yy0 new state vector
     * @param yp0 new state derivative vector
     */
    void
    reset_state(void* ida_mem, const_N_Vector yy0, const_N_Vector yp0) const;

    void set_sens_params(
        realtype const* p, realtype const* pbar, int const* plist
    ) const override;

    void adj_init() const override;

    void quad_init(AmiVector const& xQ0) const override;

    void allocate_solver_b(int* which) const override;

    void set_max_num_steps_b(int which, long int mxstepsB) const override;

    void set_ss_tolerances_b(
        int which, realtype relTolB, realtype absTolB
    ) const override;

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

#endif /* AMICI_SOLVER_IDAS_h */
