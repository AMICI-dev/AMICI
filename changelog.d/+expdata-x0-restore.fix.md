Fixed `ExpData.x0` remaining set as the model's initial state after
`run_simulation` if the model had no custom initial state (#3452).
Subsequent simulations without `ExpData.x0` used the initial state of the
earlier `ExpData`.
