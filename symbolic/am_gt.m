function fun = am_gt(varargin)
% syms x y
% f = symfun(sym('cw_gt(x,y)'),[x y]);
% fun = f(a,b);
a=varargin{1};
b=varargin{2};
fun = a-b;
if(nargin>2)
    error('N-ary gt is currently not supported!')
end
end