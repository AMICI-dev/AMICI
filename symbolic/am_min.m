function fun = am_min(a, b)
syms x y
% f = symfun(sym('am_min(x,y)'),[x y]);
% fun = f(a,b);
fun = b + heaviside(b-a)*(a-b);
end