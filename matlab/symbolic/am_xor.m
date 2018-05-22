function fun = am_xor(a,b)
% am_xor is the amici implementation of the symbolic exclusive or function
%
% Parameters:
%  a: first input parameter @type sym
%  b: second input parameter @type sym
%
% Return values:
%  fun: logical value, negative for false, positive for true

fun = am_and(am_or(a,b),-am_and(a,b));
end