function fun = am_lt(varargin)
% am_lt is the amici implementation of the n-ary mathml lessthan function
% this is an n-ary function, for more than 2 input parameters it will check
% whether and(varargin{1} < varargin{2},varargin{2} < varargin{3},...)
%
% Parameters:
%  varargin: chain of input parameters @type sym
%
% Return values:
%  fun: a < b logical value, negative for false, positive for true
a=varargin{1};
b=varargin{2};
fun = b-a;
if(nargin>2)
    fun = am_and(b-a,am_lt(varargin{2:end}));
end
end