Fixed the observable IDs created by
`amici.importers.petab.flatten_timepoint_specific_output_overrides` for
numeric overrides, which depended on how libpetab writes measurement tables.
Numbers that differ only beyond the 15th significant digit got the same ID,
and negative numbers or numbers in exponential notation resulted in invalid
IDs. Numbers are now written in their shortest exact representation, with
`p` for the decimal point and `m` for minus signs, e.g., `obs1__1p5__m2` for
the overrides `1.5` and `-2.0` of `obs1`. An error is raised if different
overrides would still result in the same observable ID.
