Fixed the computation of consistent initial conditions for DAE models
failing if the next output time is within roundoff distance of the
current time (#3353).
