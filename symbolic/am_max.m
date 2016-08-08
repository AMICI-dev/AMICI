function fun = am_max(a, b)
% syms x y
% f = symfun(sym('am_max(x,y)'),[x y]);
% fun = f(a,b);
fun = sym(['am_max(' char(vpa(a)) ',' char(vpa(b)) ',0.0)']);
% fun = b + heaviside(a-b)*(a-b);
end
