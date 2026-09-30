Events whose trigger combines multiple conditions via `And`/`Or` (as
generated, e.g., for PEtab v2 experiment/period-start events) now have
their explicit trigger time(s) resolved where possible, enabling the
"skip numerical root-finding" optimization for such events. Such
triggers are represented internally as `Min`/`Max` expressions, which
`sympy.solve` cannot solve for `t`; previously, this silently disabled
that optimization, always falling back to the slower root-finding path
(#3126).
