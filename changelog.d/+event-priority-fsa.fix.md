Fixed incorrect forward sensitivities for models with SBML events that
have explicit priorities and can trigger simultaneously.
Whenever a lower-priority event's assignment depended on a state that a
higher-priority event, triggering at the same instant, had just updated,
the sensitivity update used a stale, pre-cascade state instead of that
already-updated state. The state trajectory itself was unaffected.
