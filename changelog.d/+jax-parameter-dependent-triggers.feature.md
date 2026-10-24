The JAX backend now accepts state-updating events whose trigger depends
on parameters or on static expressions, not only on constants and time.
Previously, any such model was rejected with "The JAX backend does not
support event assignments with implicit triggers", even though only
*state*-dependent trigger times are actually unsupported.
