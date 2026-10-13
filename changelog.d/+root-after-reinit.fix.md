Fixed spurious "root after reinitialization" simulation failures
after event handling: CVODES/IDAS could re-detect an event immediately
after it was just handled and treat that as a fatal simulation error,
even though the crossing had already been processed. A re-detection of
exactly the just-processed root is now ignored, and any other,
genuinely new root (e.g. from a second, near-simultaneous event) is
handled normally instead of aborting the simulation. This also covers
discontinuities less than one solver step apart, which previously made
the simulation fail with `reInitPostProcess failed with error code 1`
(#2861, #3266).
