function fun = am_and(a,b)
% am_and is the amici implementation of the symbolic and function
%
% Parameters:
%  a: first input parameter @type sym
%  b: second input parameter @type sym
%
% Return values:
%  fun: logical value, negative for false, positive for true
fun = am_min(a,b);
end