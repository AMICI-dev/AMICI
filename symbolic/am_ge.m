function fun = am_ge(varargin)
% syms x y
% f = symfun(sym('cw_ge(x,y)'),[x y]);
% fun = f(a,b);
a=varargin{1};
b=varargin{2};
fun = a-b;
if(nargin>2)
    fun = am_and(a-b,am_ge(varargin{2:end}));
end
end