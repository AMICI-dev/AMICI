function fun = am_ge(a,b)
% syms x y
% f = symfun(sym('cw_ge(x,y)'),[x y]);
% fun = f(a,b);
fun = heaviside(a-b);
end