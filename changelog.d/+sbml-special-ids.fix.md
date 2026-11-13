SBML models with entities whose IDs are `INF`, `NaN`, `rem`, `avogadro`,
`exponentiale`, `log10`, `pi`, `true`, `false`, `True`, or `False` can now
be imported. Previously, import failed. In `MeasurementChannel` formulas and
sigmas for SBML import, the SBML L3 formula names `INF`, `NaN`, `rem`,
`avogadro`, `exponentiale`, and `log10` no longer have a special meaning.
SBML elements with ID `time` remain unsupported, but now result in a
clearer error message.
