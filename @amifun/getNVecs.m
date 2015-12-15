function this = getNVecs(this,model)
    % getfunargs populates the nvecs property with the names of the
    % N_Vector elements which are required in the execution of the function 
    % (if applicable)
    %
    % Parameters:
    %  model: model definition object @type amimodel
    %
    % Return values:
    %  this: updated function definition object @type amifun
    %
    
    % here we need to differentiate between the ODE and the DAE case. For
    % ODEs we want to use the staggered method for which the sensitivity
    % right hand side has a different function call (the parameter index is
    % passed). For other functions the loop over the parameter index is
    % realised internally and the vector of sensitivity right hand sides is
    % passed. These differences are realised by differentiating between sx
    % and *sx etc.
    if(strcmp(model.wtype,'iw'))
        x = {'x','dx'};
        sx = {'*sx','*sdx'};
        xB = {'xB','dxB'};
        sxdot = {'*sxdot'};
    else
        x = {'x'};
        sx = {'sx'};
        xB = {'xB'};
        sxdot = {'sxdot'};
    end
    psx = {'*sx'};
    xdot = {'xdot'};
    xBdot = {'xBdot'};
    qBdot = {'qBdot'};
    x0 = {'x0'};
    dx0 = {'dx0'};
    sx0 = {'*sx0'};
    sdx0 = {'sdx0'};
    v = {'v'};
    vB = {'vB'};
    Jv = {'Jv'};
    JvB = {'JvB'};
    
    
    switch(this.funstr)
        case 'xdot'
            this.nvecs = [x,xdot];
        case 'xBdot'
            this.nvecs = [x,xB,xBdot];
        case 'qBdot'
            this.nvecs = [x,xB,qBdot];
        case 'x0'
            this.nvecs = x0;
        case 'dx0'
            this.nvecs = [x0,dx0];
        case 'Jv'
            this.nvecs = [x,v,Jv];
        case 'JvB'
            this.nvecs = [x,xB,vB,JvB];
        case 'JBand'
            this.nvecs = x;
        case 'J'
            this.nvecs = x;
        case 'JSparse'
            this.nvecs = x;
        case 'JBandB'
            this.nvecs = [x,xB];
        case 'JB'
            this.nvecs = [x,xB];
        case 'JSparseB'
            this.nvecs = [x,xB];
        case 'sxdot'
            this.nvecs = [x,sx,sxdot];
        case 'sx0'
            if(strcmp(model.wtype,'iw'))
                this.nvecs = [sx0,x];
            else
                this.nvecs = sx0;
            end
        case 'sdx0'
            if(strcmp(model.wtype,'iw'))
                this.nvecs = [sdx0,x];
            else
                this.nvecs = sdx0;
            end
        case 'y'
            this.nvecs = x(1);
        case 'sy'
            this.nvecs = [x(1),psx];
        case 'z'
            this.nvecs = x(1);
        case 'sz'
            this.nvecs = [x(1),psx];
        case 'dydp'
            this.nvecs = x(1);
        case 'dydx'
            this.nvecs = x(1);
        case 'dzdp'
            this.nvecs = x(1);
        case 'dzdx'
            this.nvecs = x(1);
        case 'deltax'
            this.nvecs = x(1);
        case 'deltaxB'
            this.nvecs = [x(1),xB(1)];
        case 'deltaqB'
            this.nvecs = [x(1),xB(1),qBdot];
        case 'deltasx'
            this.nvecs = [x(1),psx];
        case 'dxdotdp'
            this.nvecs = [x(1),psx];
        case 'sigma_y'
            this.nvecs = {};
        case 'dsigma_ydp'
            this.nvecs = {};
        case 'sigma_z'
            this.nvecs = {};
        case 'dsigma_zdp'
            this.nvecs = {};
        case 'stau'
            this.nvecs = [x(1),psx];
        case 'root'
            this.nvecs = [x(1)];
        otherwise
            % do nothing
    end
    
end