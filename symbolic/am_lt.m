function fun = am_lt(varargin)
% syms x y
% f = symfun(sym('cw_lt(x,y)'),[x y]);
% fun = f(a,b);
a=varargin{1};
b=varargin{2};
fun = b-a;
if(nargin>2)
    error('N-ary lt is currently not supported!')
end
end