SBML models with entities, and measurement channels, whose IDs are `INF`,
`NaN`, `rem`, `avogadro`, `exponentiale`, `log10`, `pi`, `true`, or `false`
can now be imported. Previously, import failed. SBML entities or
measurement channels with ID `time`, `True`, or `False` remain unsupported,
but now result in a clearer error message.
