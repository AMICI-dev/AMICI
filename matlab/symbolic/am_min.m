function fun = am_min(a, b)
% am_min is the amici implementation of the symbolic min function
%
% Parameters:
%  a: first input parameter @type sym
%  b: second input parameter @type sym
%
% Return values:
%  fun: minimum of a and b
fun = -am_max(-a,-b);
end