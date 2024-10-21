#include "amici/model_state.h"

namespace amici {

ModelStateDerived::ModelStateDerived(ModelDimensions const& dim)
    : sunctx_(sundials::Context())
    , J_(dim.nx_solver, dim.nx_solver, dim.nnz, CSC_MAT, sunctx_)
    , JB_(dim.nx_solver, dim.nx_solver, dim.nnz, CSC_MAT, sunctx_)
    , dxdotdw_(dim.nx_solver, dim.nw, dim.ndxdotdw, CSC_MAT, sunctx_)
    , dx_rdatadx_solver(
          dim.nx_rdata, dim.nx_solver, dim.ndxrdatadxsolver, CSC_MAT, sunctx_
      )
    , dx_rdatadtcl(
          dim.nx_rdata, dim.nx_rdata - dim.nx_solver, dim.ndxrdatadtcl, CSC_MAT,
          sunctx_
      )
    , dtotal_cldx_rdata(
          dim.nx_rdata - dim.nx_solver, dim.nx_rdata, dim.ndtotal_cldx_rdata,
          CSC_MAT, sunctx_
      )
    , dxdotdp(0, 0, sunctx_)
    , w_(dim.nw)
    , x_rdata_(dim.nx_rdata, 0.0)
    , sx_rdata_(dim.nx_rdata, 0.0)
    , x_pos_tmp_(dim.nx_solver, sunctx_) {
    // If Matlab wrapped: dxdotdp is a full AmiVector,
    // if Python wrapped: dxdotdp_explicit and dxdotdp_implicit are CSC matrices
    if (dim.pythonGenerated) {

        dwdw_ = SUNMatrixWrapper(dim.nw, dim.nw, dim.ndwdw, CSC_MAT, sunctx_);
        // size dynamically adapted for dwdx_ and dwdp_
        dwdx_ = SUNMatrixWrapper(dim.nw, dim.nx_solver, 0, CSC_MAT, sunctx_);
        dwdp_ = SUNMatrixWrapper(dim.nw, dim.np, 0, CSC_MAT, sunctx_);

        for (int irec = 0; irec <= dim.w_recursion_depth; ++irec) {
            /* for the first element we know the exact size, while for all
               others we guess the size*/
            dwdp_hierarchical_.emplace_back(SUNMatrixWrapper(
                dim.nw, dim.np, irec * dim.ndwdw + dim.ndwdp, CSC_MAT, sunctx_
            ));
            dwdx_hierarchical_.emplace_back(SUNMatrixWrapper(
                dim.nw, dim.nx_solver, irec * dim.ndwdw + dim.ndwdx, CSC_MAT,
                sunctx_
            ));
        }
        assert(
            gsl::narrow<int>(dwdp_hierarchical_.size())
            == dim.w_recursion_depth + 1
        );
        assert(
            gsl::narrow<int>(dwdx_hierarchical_.size())
            == dim.w_recursion_depth + 1
        );

        dxdotdp_explicit = SUNMatrixWrapper(
            dim.nx_solver, dim.np, dim.ndxdotdp_explicit, CSC_MAT, sunctx_
        );
        // guess size, will be dynamically reallocated
        dxdotdp_implicit = SUNMatrixWrapper(
            dim.nx_solver, dim.np, dim.ndwdp + dim.ndxdotdw, CSC_MAT, sunctx_
        );
        dxdotdx_explicit = SUNMatrixWrapper(
            dim.nx_solver, dim.nx_solver, dim.ndxdotdx_explicit, CSC_MAT,
            sunctx_
        );
        // guess size, will be dynamically reallocated
        dxdotdx_implicit = SUNMatrixWrapper(
            dim.nx_solver, dim.nx_solver, dim.ndwdx + dim.ndxdotdw, CSC_MAT,
            sunctx_
        );
        // dynamically allocate on first call
        dxdotdp_full
            = SUNMatrixWrapper(dim.nx_solver, dim.np, 0, CSC_MAT, sunctx_);

        for (int iytrue = 0; iytrue < dim.nytrue; ++iytrue)
            dJydy_.emplace_back(SUNMatrixWrapper(
                dim.nJ, dim.ny, dim.ndJydy.at(iytrue), CSC_MAT, sunctx_
            ));

        dJydy_dense_ = SUNMatrixWrapper(dim.nJ, dim.ny, sunctx_);
    } else {
        dwdx_ = SUNMatrixWrapper(
            dim.nw, dim.nx_solver, dim.ndwdx, CSC_MAT, sunctx_
        );
        dwdp_ = SUNMatrixWrapper(dim.nw, dim.np, dim.ndwdp, CSC_MAT, sunctx_);
        dJydy_matlab_
            = std::vector<realtype>(dim.nJ * dim.nytrue * dim.ny, 0.0);
    }
}

} // namespace amici
