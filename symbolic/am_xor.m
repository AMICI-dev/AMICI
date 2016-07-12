function fun = am_xor(a,b)
% syms x y
% f = symfun(sym('cw_or(x,y)'),[x y]);
% fun = f(a,b);
fun = am_and(am_or(a,b),-am_and(a,b));
end