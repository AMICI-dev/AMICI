Event assignments to entities that are not state variables now raise an
error during model import instead of being silently ignored. PySB import
with events now disables conservation laws with a warning instead of
failing.
