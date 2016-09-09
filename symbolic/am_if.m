function fun = am_if(condition, truepart, falsepart)
% syms x y z
% f = symfun(sym('am_if(x,y,z)'),[x y z]);
% fun = f(condition, truepart, falsepart);
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