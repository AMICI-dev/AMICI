`ReturnDataView` no longer forms a reference cycle with its `xr` accessor.
`ReturnDataView`s, and the `ReturnData` they hold, are now freed as soon as
they are no longer referenced, instead of only by Python's cyclic garbage
collector. Before, dropped simulation results could accumulate across many
`run_simulations` calls, e.g. during parameter estimation, until the next
full garbage collection.
