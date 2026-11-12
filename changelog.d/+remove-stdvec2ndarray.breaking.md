Removed `stdVec2ndarray` and `SwigPyIterator.value()`. Both segfaulted on
every call, since their wrappers released the GIL while creating Python
objects ([#3509](https://github.com/AMICI-dev/AMICI/issues/3509)).
Use `numpy.array(vec)` (and `.reshape(...)`) to convert a `DoubleVector` or
`IntVector` to a `numpy.ndarray`, and `next(it)` to advance an iterator.
