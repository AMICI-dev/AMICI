Fixed (un)pickling of `Model` instances with a custom
`reinitialization_state_idxs` subset: the unpickled model reinitialized all
fixed-parameter-dependent states instead, or unpickling failed if state
reinitialization had been disabled at import
(`allow_reinit_fixpar_initcond=False`) (#3457).
