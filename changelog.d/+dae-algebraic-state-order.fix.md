For SBML models with algebraic rules, the order of the algebraic states
no longer depends on Python's string hashing (`PYTHONHASHSEED`), which
made the generated code non-reproducible and could change the state order
of the same model between imports. The algebraic states determined by
each algebraic rule are now sorted by ID.
