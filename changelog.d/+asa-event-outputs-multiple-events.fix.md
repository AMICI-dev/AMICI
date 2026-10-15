Fixed adjoint sensitivities for measurements of event outputs or root function
values of models with multiple events. The state derivatives of the event
likelihood were stored per occurrence index only, so the contributions of
different events with the same occurrence index overwrote each other.
