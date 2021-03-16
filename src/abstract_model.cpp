#include "amici/abstract_model.h"

namespace amici {

std::string
AbstractModel::getAmiciVersion() const
{
    throw AmiException("Version not set during code generation");
}

std::string
AbstractModel::getAmiciCommit() const
{
    throw AmiException("Commit not set during code generation");
}

void
AbstractModel::fx0(realtype* /*x0*/,
                   const realtype /*t*/,
                   const realtype* /*p*/,
                   const realtype* /*k*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

bool
AbstractModel::isFixedParameterStateReinitializationAllowed() const
{
    return false;
}

void
AbstractModel::fx0_fixedParameters(realtype* /*x0*/,
                                   const realtype /*t*/,
                                   const realtype* /*p*/,
                                   const realtype* /*k*/,
                                   gsl::span<const int> /*reinitialization_state_idxs*/)
{
    // no-op default implementation
}

void
AbstractModel::fsx0_fixedParameters(realtype* /*sx0*/,
                                    const realtype /*t*/,
                                    const realtype* /*x0*/,
                                    const realtype* /*p*/,
                                    const realtype* /*k*/,
                                    const int /*ip*/,
                                    gsl::span<const int> /*reinitialization_state_idxs*/)
{
    // no-op default implementation
}

void
AbstractModel::fsx0(realtype* /*sx0*/,
                    const realtype /*t*/,
                    const realtype* /*x0*/,
                    const realtype* /*p*/,
                    const realtype* /*k*/,
                    const int /*ip*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdx0(AmiVector& /*x0*/, AmiVector& /*dx0*/)
{
    // no-op default implementation
}

void
AbstractModel::fstau(realtype* /*stau*/,
                     const realtype /*t*/,
                     const realtype* /*x*/,
                     const realtype* /*p*/,
                     const realtype* /*k*/,
                     const realtype* /*h*/,
                     const realtype* /*sx*/,
                     const int /*ip*/,
                     const int /*ie*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fy(realtype* /*y*/,
                  const realtype /*t*/,
                  const realtype* /*x*/,
                  const realtype* /*p*/,
                  const realtype* /*k*/,
                  const realtype* /*h*/,
                  const realtype* /*w*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdydp(realtype* /*dydp*/,
                     const realtype /*t*/,
                     const realtype* /*x*/,
                     const realtype* /*p*/,
                     const realtype* /*k*/,
                     const realtype* /*h*/,
                     const int /*ip*/,
                     const realtype* /*w*/,
                     const realtype* /*dwdp*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdydx(realtype* /*dydx*/,
                     const realtype /*t*/,
                     const realtype* /*x*/,
                     const realtype* /*p*/,
                     const realtype* /*k*/,
                     const realtype* /*h*/,
                     const realtype* /*w*/,
                     const realtype* /*dwdx*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fz(realtype* /*z*/,
                  const int /*ie*/,
                  const realtype /*t*/,
                  const realtype* /*x*/,
                  const realtype* /*p*/,
                  const realtype* /*k*/,
                  const realtype* /*h*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fsz(realtype* /*sz*/,
                   const int /*ie*/,
                   const realtype /*t*/,
                   const realtype* /*x*/,
                   const realtype* /*p*/,
                   const realtype* /*k*/,
                   const realtype* /*h*/,
                   const realtype* /*sx*/,
                   const int /*ip*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::frz(realtype* /*rz*/,
                   const int /*ie*/,
                   const realtype /*t*/,
                   const realtype* /*x*/,
                   const realtype* /*p*/,
                   const realtype* /*k*/,
                   const realtype* /*h*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fsrz(realtype* /*srz*/,
                    const int /*ie*/,
                    const realtype /*t*/,
                    const realtype* /*x*/,
                    const realtype* /*p*/,
                    const realtype* /*k*/,
                    const realtype* /*h*/,
                    const realtype* /*sx*/,
                    const int /*ip*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdzdp(realtype* /*dzdp*/,
                     const int /*ie*/,
                     const realtype /*t*/,
                     const realtype* /*x*/,
                     const realtype* /*p*/,
                     const realtype* /*k*/,
                     const realtype* /*h*/,
                     const int /*ip*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdzdx(realtype* /*dzdx*/,
                     const int /*ie*/,
                     const realtype /*t*/,
                     const realtype* /*x*/,
                     const realtype* /*p*/,
                     const realtype* /*k*/,
                     const realtype* /*h*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdrzdp(realtype* /*drzdp*/,
                      const int /*ie*/,
                      const realtype /*t*/,
                      const realtype* /*x*/,
                      const realtype* /*p*/,
                      const realtype* /*k*/,
                      const realtype* /*h*/,
                      const int /*ip*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdrzdx(realtype* /*drzdx*/,
                      const int /*ie*/,
                      const realtype /*t*/,
                      const realtype* /*x*/,
                      const realtype* /*p*/,
                      const realtype* /*k*/,
                      const realtype* /*h*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdeltax(realtype* /*deltax*/,
                       const realtype /*t*/,
                       const realtype* /*x*/,
                       const realtype* /*p*/,
                       const realtype* /*k*/,
                       const realtype* /*h*/,
                       const int /*ie*/,
                       const realtype* /*xdot*/,
                       const realtype* /*xdot_old*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdeltasx(realtype* /*deltasx*/,
                        const realtype /*t*/,
                        const realtype* /*x*/,
                        const realtype* /*p*/,
                        const realtype* /*k*/,
                        const realtype* /*h*/,
                        const realtype* /*w*/,
                        const int /*ip*/,
                        const int /*ie*/,
                        const realtype* /*xdot*/,
                        const realtype* /*xdot_old*/,
                        const realtype* /*sx*/,
                        const realtype* /*stau*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdeltaxB(realtype* /*deltaxB*/,
                        const realtype /*t*/,
                        const realtype* /*x*/,
                        const realtype* /*p*/,
                        const realtype* /*k*/,
                        const realtype* /*h*/,
                        const int /*ie*/,
                        const realtype* /*xdot*/,
                        const realtype* /*xdot_old*/,
                        const realtype* /*xB*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdeltaqB(realtype* /*deltaqB*/,
                        const realtype /*t*/,
                        const realtype* /*x*/,
                        const realtype* /*p*/,
                        const realtype* /*k*/,
                        const realtype* /*h*/,
                        const int /*ip*/,
                        const int /*ie*/,
                        const realtype* /*xdot*/,
                        const realtype* /*xdot_old*/,
                        const realtype* /*xB*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fsigmay(realtype* /*sigmay*/,
                       const realtype /*t*/,
                       const realtype* /*p*/,
                       const realtype* /*k*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdsigmaydp(realtype* /*dsigmaydp*/,
                          const realtype /*t*/,
                          const realtype* /*p*/,
                          const realtype* /*k*/,
                          const int /*ip*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fsigmaz(realtype* /*sigmaz*/,
                       const realtype /*t*/,
                       const realtype* /*p*/,
                       const realtype* /*k*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdsigmazdp(realtype* /*dsigmazdp*/,
                          const realtype /*t*/,
                          const realtype* /*p*/,
                          const realtype* /*k*/,
                          const int /*ip*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fJy(realtype* /*nllh*/,
                   const int /*iy*/,
                   const realtype* /*p*/,
                   const realtype* /*k*/,
                   const realtype* /*y*/,
                   const realtype* /*sigmay*/,
                   const realtype* /*my*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fJz(realtype* /*nllh*/,
                   const int /*iz*/,
                   const realtype* /*p*/,
                   const realtype* /*k*/,
                   const realtype* /*z*/,
                   const realtype* /*sigmaz*/,
                   const realtype* /*mz*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fJrz(realtype* /*nllh*/,
                    const int /*iz*/,
                    const realtype* /*p*/,
                    const realtype* /*k*/,
                    const realtype* /*z*/,
                    const realtype* /*sigmaz*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdJydy(realtype* /*dJydy*/,
                      const int /*iy*/,
                      const realtype* /*p*/,
                      const realtype* /*k*/,
                      const realtype* /*y*/,
                      const realtype* /*sigmay*/,
                      const realtype* /*my*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdJydy_colptrs(SUNMatrixWrapper &/*indexptrs*/,
                              int /*index*/) {
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdJydy_rowvals(SUNMatrixWrapper & /*indexptrs*/,
                              int /*index*/) {
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdJydsigma(realtype* /*dJydsigma*/,
                          const int /*iy*/,
                          const realtype* /*p*/,
                          const realtype* /*k*/,
                          const realtype* /*y*/,
                          const realtype* /*sigmay*/,
                          const realtype* /*my*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdJzdz(realtype* /*dJzdz*/,
                      const int /*iz*/,
                      const realtype* /*p*/,
                      const realtype* /*k*/,
                      const realtype* /*z*/,
                      const realtype* /*sigmaz*/,
                      const realtype* /*mz*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdJzdsigma(realtype* /*dJzdsigma*/,
                          const int /*iz*/,
                          const realtype* /*p*/,
                          const realtype* /*k*/,
                          const realtype* /*z*/,
                          const realtype* /*sigmaz*/,
                          const realtype* /*mz*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdJrzdz(realtype* /*dJrzdz*/,
                       const int /*iz*/,
                       const realtype* /*p*/,
                       const realtype* /*k*/,
                       const realtype* /*rz*/,
                       const realtype* /*sigmaz*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdJrzdsigma(realtype* /*dJrzdsigma*/,
                           const int /*iz*/,
                           const realtype* /*p*/,
                           const realtype* /*k*/,
                           const realtype* /*rz*/,
                           const realtype* /*sigmaz*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fw(realtype* /*w*/,
                  const realtype /*t*/,
                  const realtype* /*x*/,
                  const realtype* /*p*/,
                  const realtype* /*k*/,
                  const realtype* /*h*/,
                  const realtype* /*tcl*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdwdp(realtype* /*dwdp*/,
                     const realtype /*t*/,
                     const realtype* /*x*/,
                     const realtype* /*p*/,
                     const realtype* /*k*/,
                     const realtype* /*h*/,
                     const realtype* /*w*/,
                     const realtype* /*tcl*/,
                     const realtype* /*stcl*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdwdp_colptrs(SUNMatrixWrapper &/*dwdp*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdwdp_rowvals(SUNMatrixWrapper &/*dwdp*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdwdp(realtype* /*dwdp*/,
                     const realtype /*t*/,
                     const realtype* /*x*/,
                     const realtype* /*p*/,
                     const realtype* /*k*/,
                     const realtype* /*h*/,
                     const realtype* /*w*/,
                     const realtype* /*tcl*/,
                     const realtype* /*stcl*/,
                     const int /*ip*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdwdx(realtype* /*dwdx*/,
                     const realtype /*t*/,
                     const realtype* /*x*/,
                     const realtype* /*p*/,
                     const realtype* /*k*/,
                     const realtype* /*h*/,
                     const realtype* /*w*/,
                     const realtype* /*tcl*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdwdx_colptrs(SUNMatrixWrapper &/*dwdx*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdwdx_rowvals(SUNMatrixWrapper &/*dwdx*/)
{
    throw AmiException("Requested functionality is not supported as %s is "
                       "not implemented for this model!",
                       __func__);
}

void
AbstractModel::fdwdw(realtype */*dwdw*/,
                     realtype /*t*/,
                     const realtype */*x*/,
                     const realtype */*p*/,
                     const realtype */*k*/,
                     const realtype */*h*/,
                     const realtype */*w*/,
                     const realtype */*tcl*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__);
}

void AbstractModel::fdwdw_colptrs(SUNMatrixWrapper &/*dwdw*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__);
}

void AbstractModel::fdwdw_rowvals(SUNMatrixWrapper &/*dwdw*/) {
    throw AmiException("Requested functionality is not supported as %s "
                       "is not implemented for this model!",
                       __func__);
}

} // namespace amici
