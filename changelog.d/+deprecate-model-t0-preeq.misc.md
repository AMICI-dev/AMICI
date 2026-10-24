`Model.set_t0_preeq()` and `Model.t0_preeq()` are deprecated; set
`ExpData.t_start_preeq` instead. Pre-equilibration only happens with an
`ExpData`, so its start time belongs there (#3455).
