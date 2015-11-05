function fun = am_lt(a,b)
% syms x y
% f = symfun(sym('cw_lt(x,y)'),[x y]);
% fun = f(a,b);
fun = 1-heaviside(a-b);
end