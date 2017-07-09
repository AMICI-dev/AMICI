function model = neuron_syms()
    
    model.param = 'log10';
    
    syms a b c d 
    
    p = [a b c d];
    
    syms v0 I0
    
    k = [v0,I0];
    
    syms v u
    
    x = [v u];
    
    syms I t
    
    I = I0;
    
    f(1) = 0.04*v^2 + 5*v + 140 - u + I ;
    f(2) = a*(b*v - u);

    y(1) = v;
    
    
    x0 = [v0,b*v0];
    
    event = amievent(v-30,[-c-v,d],t);
    
    model.sym.p = p;
    model.sym.k = k;
    model.sym.x = x;
    model.sym.y = y;
    model.sym.f = f;
    model.event = event;
    model.sym.x0 = x0;
    
    
end