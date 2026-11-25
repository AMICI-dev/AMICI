Model builds now use `cmake` and `ninja` of the Python environment of the
running interpreter, also if that environment is not activated. Previously,
they were only looked up on `PATH`, which failed with "Required command
'ninja' not found" or silently used other versions (#3505).
