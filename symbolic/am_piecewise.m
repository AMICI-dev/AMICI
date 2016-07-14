function fun = am_piecewise( piece,condition,default )
    fun = piece*heaviside(sym(condition)) + default*heaviside(-sym(condition));
end

