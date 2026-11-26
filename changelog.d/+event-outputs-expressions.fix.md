Fixed the generated code of event outputs (and of the root function values
reported for events that did not trigger) that refer to expressions, e.g.,
targets of SBML assignment rules. Such models previously failed to compile
(#3558). Event outputs that depend on splines now raise a
`NotImplementedError` during import instead of producing code that does not
compile.
