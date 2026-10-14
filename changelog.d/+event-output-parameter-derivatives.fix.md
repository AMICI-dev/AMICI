Fixed forward sensitivities of event outputs and of the root function
values of events that do not trigger, if they depend explicitly on more
than one parameter, and the corresponding contributions to `sllh`
(#3361). For events that did not trigger until the last output time, the
contribution of the root function value to `sllh` used the parameter
derivative of the event output instead of the root function value.
