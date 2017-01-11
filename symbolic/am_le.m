function fun = am_le(varargin)
% syms x y
% f = symfun(sym('cw_le(x,y)'),[x y]);
% fun = f(a,b);
a=varargin{1};
b=varargin{2};
fun = b-a;
if(nargin>2)
    error('N-ary le is currently not supported!')
end
end