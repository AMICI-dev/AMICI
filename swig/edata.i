%module edata

// Add necessary symbols to generated header
%{
#include "amici/edata.h"
using namespace amici;
%}

%ignore ConditionContext;

// ExpData.__repr__
%pythoncode %{
def _edata_repr(self: "ExpData"):
    n_data_y = sum(
        self.isSetObservedData(it, iy)
        for it in range(self.nt()) for
        iy in range(self.nytrue())
    )
    n_sigma_y = sum(
        self.isSetObservedDataStdDev(it, iy)
        for it in range(self.nt())
        for iy in range(self.nytrue())
    )
    n_data_z = sum(
        self.isSetObservedData(ie, iz)
        for ie in range(self.nmaxevent())
        for iz in range(self.nztrue())
    )
    n_sigma_z = sum(
        self.isSetObservedDataStdDev(ie, iz)
        for ie in range(self.nmaxevent())
        for iz in range(self.nztrue())
    )
    return "\n".join([
        self.this.__repr__()[:-1],
        f"  {self.nt()}x{self.nytrue()} time-resolved datapoints",
        f"    ({n_data_y}/{self.nt()*self.nytrue()} measurements & {n_sigma_y}/{self.nt()*self.nytrue()} sigmas set)",
        f"  {self.nmaxevent()}x{self.nztrue()} event-resolved datapoints",
        f"    ({n_data_z}/{self.nmaxevent()*self.nztrue()} measurements & {n_sigma_z}/{self.nmaxevent()*self.nztrue()} sigmas set)",
        ">"
    ])
%}
%extend amici::ExpData {
%pythoncode %{
def __repr__(self):
    return _edata_repr(self)
%}
};
%extend std::unique_ptr<amici::ExpData> {
%pythoncode %{
def __repr__(self):
    return _edata_repr(self)
%}
};


// Process symbols in header
%include "amici/edata.h"
