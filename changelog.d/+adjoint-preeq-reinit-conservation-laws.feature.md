Adjoint pre-equilibration with reinitialization of non-constant states is now
supported for models with conservation laws, as long as no reinitialized state
shares a conservation law with other states (#1156). Previously, this raised
an error for any model with conservation laws.
