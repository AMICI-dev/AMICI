Model extensions are now built with the `cmake` and `ninja` of the current
Python environment also if that environment is not activated, i.e., if its
scripts directory is not on `PATH`. Previously, the build failed with
`Required command 'ninja' not found` or used another `cmake` from `PATH`
([#3505](https://github.com/AMICI-dev/AMICI/issues/3505)).
