function fun = am_max(a, b)
% am_max is the amici implementation of the symbolic max function
%
% Parameters:
%  a: first input parameter @type sym
%  b: second input parameter @type sym
%
% Return values:
%  fun: maximum of a and b
fun = sym(['am_max(' char(sym(a)) ',' char(sym(b)) ',0.0)']);
end
