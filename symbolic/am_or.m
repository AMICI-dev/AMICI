function fun = am_or(a,b)
% am_or is the amici implementation of the symbolic or function
%
% Parameters:
%  a: first input parameter @type sym
%  b: second input parameter @type sym
%
% Return values:
%  fun: logical value, negative for false, positive for true
fun = am_max(a,b);
end