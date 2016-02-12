function fun = am_stepfun(t,tstart,vstart,tend,vend)
% syms x y
% f = symfun(sym('am_min(x,y)'),[x y]);
% fun = f(a,b);
fun = heaviside(t-tstart)*vstart - heaviside(t-tend)*(vstart-vend);
end