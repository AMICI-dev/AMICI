#ifndef AMICI_SUNDIALS_LINSOL_WRAPPER_H
#define AMICI_SUNDIALS_LINSOL_WRAPPER_H

#include <amici/exception.h>
#include <sundials/sundials_linearsolver.h> // SUNLinearSolver
#include <sunlinsol/sunlinsol_band.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_klu.h>
#include <sunlinsol/sunlinsol_pcg.h>
#include <sunlinsol/sunlinsol_spbcgs.h>
#include <sunlinsol/sunlinsol_spfgmr.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunlinsol/sunlinsol_sptfqmr.h>

namespace amici {

/**
 * @brief A RAII wrapper for SUNLinearSolver structs.
 *
 * More members to be added as required.
 * Start accessing through SUNLinSolWrapper::linsol().
 */
class SUNLinSolWrapper {
public:
    /**
     * @brief SUNLinSolWrapper from existing SUNLinearSolver
     * @param linsol_
     */
    explicit SUNLinSolWrapper(SUNLinearSolver linsol_);

    virtual ~SUNLinSolWrapper();

    /**
     * @brief Copy constructor
     * @param other
     */
    SUNLinSolWrapper(const SUNLinSolWrapper& other) = delete;

    /**
     * @brief Move constructor
     * @param other
     */
    SUNLinSolWrapper(SUNLinSolWrapper&& other) noexcept;

    /**
     * @brief Copy assignment
     * @param other
     * @return
     */
    SUNLinSolWrapper& operator=(const SUNLinSolWrapper& other);

    /**
     * @brief Move assignment
     * @param other
     * @return
     */
    SUNLinSolWrapper& operator=(SUNLinSolWrapper&& other) noexcept;

    /**
     * @brief Get the wrapped SUNLinSol
     * @return SlsMat
     */
    SUNLinearSolver get() const;

    /**
     * @brief getType
     * @return
     */
    SUNLinearSolver_Type getType() const;

protected:
    /** the wrapper solver */
    SUNLinearSolver linsol = nullptr;
};


class SUNLinSolBand: public SUNLinSolWrapper {
    /**
     * @brief SUNLinSolBand
     * @param y
     * @param A
     */
    SUNLinSolBand(N_Vector y, SUNMatrix A)
        : SUNLinSolWrapper (SUNLinSol_Band(y, A))
    {
        if(!linsol)
            throw AmiException("Failed to create solver.");
    }
};


class SUNLinSolDense: public SUNLinSolWrapper {
    /**
     * @brief SUNLinSolDense
     * @param y
     * @param A
     */
    SUNLinSolDense(N_Vector y, SUNMatrix A)
        : SUNLinSolWrapper (SUNLinSol_Dense(y, A))
    {
        if(!linsol)
            throw AmiException("Failed to create solver.");
    }
};


class SUNLinSolKLU : public SUNLinSolWrapper {
    /**
     * @brief SUNLinSolKLU
     * @param y
     * @param A
     */
    SUNLinSolKLU(N_Vector y, SUNMatrix A)
        : SUNLinSolWrapper(SUNLinSol_KLU(y, A))
    {
        if(!linsol)
            throw AmiException("Failed to create solver.");
    }
};


class SUNLinSolPCG: public SUNLinSolWrapper {
    /**
     * @brief SUNLinSolPCG
     * @param y
     * @param pretype
     * @param maxl
     */
    SUNLinSolPCG(N_Vector y, int pretype, int maxl)
        : SUNLinSolWrapper(SUNLinSol_PCG(y, pretype, maxl))
    {
        if(!linsol)
            throw AmiException("Failed to create solver.");
    }
};


class SUNLinSolSPBCGS : public SUNLinSolWrapper {
    SUNLinSolSPBCGS(N_Vector y, int pretype, int maxl)
        : SUNLinSolWrapper(SUNLinSol_SPBCGS(y, pretype, maxl))
    {
        if(!linsol)
            throw AmiException("Failed to create solver.");
    }
};


class SUNLinSolSPFGMR: public SUNLinSolWrapper {
    /**
     * @brief SUNLinSolSPFGMR
     * @param y
     * @param pretype
     * @param maxl
     */
    SUNLinSolSPFGMR(N_Vector y, int pretype, int maxl)
        : SUNLinSolWrapper(SUNLinSol_SPFGMR(y, pretype, maxl))
    {
        if(!linsol)
            throw AmiException("Failed to create solver.");
    }
};


class SUNLinSolSPGMR: public SUNLinSolWrapper {
    /**
     * @brief SUNLinSolSPGMR
     * @param y
     * @param pretype
     * @param maxl
     */
    SUNLinSolSPGMR(N_Vector y, int pretype, int maxl)
        : SUNLinSolWrapper(SUNLinSol_SPGMR(y, pretype, maxl))
    {
        if(!linsol)
            throw AmiException("Failed to create solver.");
    }
};


class SUNLinSolSPTFQMR: public SUNLinSolWrapper {
    /**
     * @brief SUNLinSolSPTFQMR
     * @param y
     * @param pretype
     * @param maxl
     */
    SUNLinSolSPTFQMR(N_Vector y, int pretype, int maxl)
        : SUNLinSolWrapper(SUNLinSol_SPTFQMR(y, pretype, maxl))
    {
        if(!linsol)
            throw AmiException("Failed to create solver.");
    }
};

} // namespace amici
#endif // AMICI_SUNDIALS_LINSOL_WRAPPER_H

