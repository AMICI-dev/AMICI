function fun = am_stepfun(t,tstart,vstart,tend,vend)
% am_stepfun is the amici implementation of the step function 
%
% Parameters:
%  t: input variable @type sym
%  tstart: input variable value at which the step starts
%  vstart: value during the step
%  tend: input variable value at which the step end
%  vend: value after the step
%
% Return values:
%  fun: 0 before tstart, vstart between tstart and tend and vend after tend
fun = heaviside(t-tstart)*vstart - heaviside(t-tend)*(vstart-vend);
end