Integers outside the 32-bit integer range in model expressions are now
printed as floating-point literals in the generated C++ and JAX code.
Previously, integers outside the 64-bit range were silently truncated by the
C++ compiler (with only a warning), and the JAX backend could fail with an
`OverflowError` (#3405).
