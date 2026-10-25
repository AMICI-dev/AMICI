`ReturnData.sx0` is no longer populated with meaningless values when
pre-equilibration is run with adjoint sensitivities together with
adjoint sensitivity analysis of the main simulation. These values were
never actually computed in that case and were not used for the
gradient; `rdata.sx0` is now empty (`None` in Python), consistent with
other fields that are unavailable for the chosen settings (#1184).
