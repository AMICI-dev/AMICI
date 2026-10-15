Fixed the computation of consistent initial conditions for DAE models
failing if the next output time is within roundoff distance of the
current time or not finite (#3353). This also applies to the consistent
states after events.
