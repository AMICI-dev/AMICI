function fun = am_piecewise( piece,condition,default )
% am_piecewise is the amici implementation of the mathml piecewise function
%
% Parameters:
%  piece: value if condition is true
%  condition: logical value
%  default: value if condition is false
%
% Return values:
%  fun: return value, piece if condition is true, default if not
    fun = am_if(condition,piece,default);
end

