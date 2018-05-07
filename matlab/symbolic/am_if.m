function fun = am_if(condition, truepart, falsepart)
% am_if is the amici implementation of the symbolic if function
%
% Parameters:
%  condition: logical value @type sym
%  truepart: value if condition is true @type sym
%  falsepart: value if condition is false @type sym
%
% Return values:
%  fun: if condition is true truepart, else falsepart
if(islogical(condition))
    if(condition)
        fun = truepart;
    else
        fun = falsepart;
    end
else
    if(logical(condition~=0))
        fun = falsepart + heaviside(condition)*(truepart-falsepart);
    else
        fun = falsepart;
    end
end
end