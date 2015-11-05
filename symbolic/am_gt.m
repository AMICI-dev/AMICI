function fun = am_gt(a,b)
% syms x y
% f = symfun(sym('cw_gt(x,y)'),[x y]);
% fun = f(a,b);
fun = 1-heaviside(b-a);
end