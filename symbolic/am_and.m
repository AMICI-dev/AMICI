function fun = am_and(a,b)
% syms x y
% f = symfun(sym('cw_and(x,y)'),[x y]);
% fun = f(a,b);
fun = heaviside((a+b)/2-1);
end