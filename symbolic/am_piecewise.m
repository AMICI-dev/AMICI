function fun = am_piecewise( piece,condition,default )
    fun = piece*heaviside(condition) + default*heaviside(-condition);
end

