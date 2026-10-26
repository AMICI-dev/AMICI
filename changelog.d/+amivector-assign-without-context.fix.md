Fixed a use-after-free when assigning an `AmiVector` without SUNDIALS
context, e.g., an empty one, to an `AmiVector` with an `N_Vector`.
