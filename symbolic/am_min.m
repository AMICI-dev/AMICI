function fun = am_min(a, b)
% syms x y
% f = symfun(sym('am_min(x,y)'),[x y]);
% fun = f(a,b);
fun = sym(['am_min(' char(sym(a)) ',' char(sym(b)) ',0.0)']);
% fun = b + heaviside(b-a)*(a-b);
end