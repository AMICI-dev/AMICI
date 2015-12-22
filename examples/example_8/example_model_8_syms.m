function model = example_model_8_syms()
    
    model.param = 'lin';
    model.recompile = true;
    
    syms a b c d
    
    p = [a b c d];
    
    syms u0 v0 tI I0
    
    k = [v0 I0];
    
    syms v u
    
    x = v;
    
    syms I t
    
    I = I0;

    f(1) = -a*x^2 + x^2*dirac(v-1/2);
    
    y = v;
    
    x0 = v0;
    
    sym('t');
    event(1) = amievent(v-1/2,x^2,t);
    
    model.sym.p = p;
    model.sym.k = k;
    model.sym.x = x;
    model.sym.y = y;
    model.sym.f = f;
    model.sym.x0 = x0;
    model.event = event;
    
    
end