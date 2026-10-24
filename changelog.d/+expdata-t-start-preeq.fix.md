Fixed `ExpData.t_start_preeq` being ignored: pre-equilibration started at
`Model.t0_preeq()`, or at `ExpData.t_start` if that was NaN (#3455).
Now, pre-equilibration starts at `ExpData.t_start_preeq`, or at
`ExpData.t_start` if that is NaN. This follows the rule that all settings
defining the simulation periods (`t_start`, `t_presim`, the
pre-equilibration and presimulation settings, and reinitialization) come
from `ExpData` only; `Model.t0_preeq()` was the only exception, and
pre-equilibration only happens with an `ExpData` anyway. For now, the
deprecated `Model.t0_preeq()` is still used if `ExpData.t_start_preeq` is
NaN, with a warning. `get_model_for_preeq` follows the same rules.
This changes results of PEtab v2 problems with pre-equilibration whose main
simulation does not start at `t = 0`: `ReturnData.preeq_t` (absolute time)
shifts accordingly, and states change if the model depends explicitly on
time. In the PEtab v2 simulator, a user-set `Model.set_t0_preeq()` is now
overridden by PEtab's pre-equilibration start time `t = 0`.
