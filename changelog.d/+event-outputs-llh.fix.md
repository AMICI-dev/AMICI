Fixed the likelihood of models with more than one event output: only the
first event output of the model was included, the evaluation could write
past the end of a buffer (heap corruption, crashes on Windows), and outputs
of other events were included at each event.
