function fun = am_or(a,b)
% syms x y
% f = symfun(sym('cw_or(x,y)'),[x y]);
% fun = f(a,b);
fun = 1-am_and(a,b);
end