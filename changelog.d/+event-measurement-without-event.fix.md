For event measurements of events that did not occur until the last output
time, only the root function value enters the likelihood, consistent with the
likelihood sensitivities. Before, the likelihood also included the event
output at the last output time, which was not accounted for in the
sensitivities. This changes `ReturnData.llh` for such data (#3386).
