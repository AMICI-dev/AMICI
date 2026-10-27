Fixed heap corruption and crashes for DAE models with both events with
explicit trigger times (e.g., `time >= 1`) and events tracked by the solver
(#3449). The IDAS root function callback wrote all root function values into
the solver's buffer, which only holds those tracked by the solver.
