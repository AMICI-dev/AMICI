Model entity IDs no longer have to avoid names used by AMICI or its
generated code, for either backend: previously, model
import or compilation could fail — or silently generate incorrect code,
with no indication of the actual cause — whenever a model entity's ID
collided with something the generated code already used. That was a C++
keyword or standard-library macro (e.g. `NULL`, `EOF`) for the SUNDIALS
backend; a Python keyword, one of the generated module's imports
(`jnp`, `safe_log`, ...) or one of its methods' own arguments (`self`,
`tcl`, `my`, ...) for the JAX backend; AMICI's fixed array-parameter
names (`x`, `p`, `k`, `h`, `w`, `y`) for both; or an
AMICI-internally-derived symbol name (e.g. a reaction's own flux symbol,
a conservation-law total, an observable's measurement/sigma symbol, or a
sparse-Jacobian placeholder like the `dwdx`/`dwdp`/`dxdotdx_explicit`
entries named `d{row}_d{col}`). In particular, single-letter entity names
used to be a common source of such failures.
All of these are now handled transparently: each backend's code printer
mangles every identifier it prints, so a name that's unsafe in generated
code stays an implementation detail of code generation. Only `t` is still
renamed at the model level (it denotes simulation time symbolically, not
just by name), with the original ID restored everywhere it's reported
(state/parameter/observable/expression IDs, PEtab mapping, ...) — this
internal renaming previously used an `amici_` prefix that could itself
leak into reported IDs; that prefix is now purely an implementation
detail and never visible to users. Generated locals are also never
aliased across two different quantities that merely happen to print the
same name, and generated C++ files no longer rely on unscoped `#define`
macros for entity names either.

As a side effect, the local variable names used internally in generated
model code have changed (e.g. `STAT` is now printed as `STAT_`). This
does not affect the public API — state/parameter/observable IDs are
unaffected — but if you inspect or post-process AMICI-generated C++ or
JAX source directly, expect different local variable names (#2226,
#2461, #3237, #3240, #3243, #3246, #3255).
