function generateM(this, amimodelo2)
    % generateM generates the matlab wrapper for the compiled C files.
    %
    % Parameters:
    %  amimodelo2: this struct must contain all necessary symbolic
    %  definitions for second order sensivities @type amimodel
    %
    % Return values:
    %  this: model definition object @type amimodel
    
    nx = this.nx;
    ny = this.ny;
    np = this.np;
    nk = this.nk;
    nz = this.nz;
    if(~isempty(amimodelo2))
        nztrue = amimodelo2.nztrue;
        nxtrue = amimodelo2.nxtrue;
        nytrue = amimodelo2.nytrue;
        o2flag = amimodelo2.o2flag;
    else
        nztrue = nz;
        nxtrue = nx;
        nytrue = ny;
        o2flag = this.o2flag;
    end
    
    
    
    if(o2flag)
        nxtrue = amimodelo2.nxtrue;
        nytrue = amimodelo2.nytrue;
    end
    
    fid = fopen(fullfile(this.wrap_path,'models',this.modelname,['simulate_',this.modelname,'.m']),'w');
    fprintf(fid,['%% simulate_' this.modelname '.m is the matlab interface to the cvodes mex\n'...
        '%%   which simulates the ordinary differential equation and respective\n'...
        '%%   sensitivities according to user specifications.\n'...
        '%%\n'...
        '%% USAGE:\n'...
        '%% ======\n'...
        '%% [...] = simulate_' this.modelname '(tout,theta)\n'...
        '%% [...] = simulate_' this.modelname '(tout,theta,kappa,data,options)\n'...
        '%% [status,tout,x,y,sx,sy] = simulate_' this.modelname '(...)\n'...
        '%%\n'...
        '%% INPUTS:\n'...
        '%% =======\n'...
        '%% tout ... 1 dimensional vector of timepoints at which a solution to the ODE is desired\n'...
        '%% theta ... 1 dimensional parameter vector of parameters for which sensitivities are desired.\n'...
        '%%           this corresponds to the specification in model.sym.p\n'...
        '%% kappa ... 1 dimensional parameter vector of parameters for which sensitivities are not desired.\n'...
        '%%           this corresponds to the specification in model.sym.k\n'...
        '%% data ... struct containing the following fields. Can have the following fields '...
        '%%     Y ... 2 dimensional matrix containing data.\n'...
        '%%           columns must correspond to observables and rows to time-points\n'...
        '%%     Sigma_Y ... 2 dimensional matrix containing standard deviation of data.\n'...
        '%%           columns must correspond to observables and rows to time-points\n'...
        '%%     T ... (optional) 2 dimensional matrix containing events.\n'...
        '%%           columns must correspond to event-types and rows to possible event-times\n'...
        '%%     Sigma_T ... (optional) 2 dimensional matrix containing standard deviation of events.\n'...
        '%%           columns must correspond to event-types and rows to possible event-times\n'...
        '%% options ... additional options to pass to the cvodes solver. Refer to the cvodes guide for more documentation.\n'...
        '%%    .cvodes_atol ... absolute tolerance for the solver. default is specified in the user-provided syms function.\n'...
        '%%    .cvodes_rtol ... relative tolerance for the solver. default is specified in the user-provided syms function.\n'...
        '%%    .cvodes_maxsteps    ... maximal number of integration steps. default is specified in the user-provided syms function.\n'...
        '%%    .tstart    ... start of integration. for all timepoints before this, values will be set to initial value.\n'...
        '%%    .sens_ind ... 1 dimensional vector of indexes for which sensitivities must be computed.\n'...
        '%%           default value is 1:length(theta).\n'...
        '%%    .sx0 ... user-provided sensitivity initialisation. this should be a matrix of dimension [#states x #parameters].\n'...
        '%%        default is sensitivity initialisation based on the derivative of the state initialisation.'...
        '%%    .lmm    ... linear multistep method for forward problem.\n'...
        '%%        1: Adams-Bashford\n'...
        '%%        2: BDF (DEFAULT)\n'...
        '%%    .iter    ... iteration method for linear multistep.\n'...
        '%%        1: Functional\n'...
        '%%        2: Newton (DEFAULT)\n'...
        '%%    .linsol   ... linear solver module.\n'...
        '%%        direct solvers:\n'...
        '%%        1: Dense (DEFAULT)\n'...
        '%%        2: Band (not implented)\n'...
        '%%        3: LAPACK Dense (not implented)\n'...
        '%%        4: LAPACK Band  (not implented)\n'...
        '%%        5: Diag (not implented)\n'...
        '%%        implicit krylov solvers:\n'...
        '%%        6: SPGMR\n'...
        '%%        7: SPBCG\n'...
        '%%        8: SPTFQMR\n'...
        '%%        sparse solvers:\n'...
        '%%        9: KLU\n'...
        '%%    .stldet   ... flag for stability limit detection. this should be turned on for stiff problems.\n'...
        '%%        0: OFF\n'...
        '%%        1: ON (DEFAULT)\n'...
        '%%    .qPositiveX   ... vector of 0 or 1 of same dimension as state vector. 1 enforces positivity of states.\n'...
        '%%    .sensi_meth   ... method for sensitivity analysis.\n'...
        '%%        ''forward'': forward sensitivity analysis (DEFAULT)\n'...
        '%%        ''adjoint'': adjoint sensitivity analysis \n'...
        '%%        ''ss'': steady state sensitivity analysis \n'...
        '%%    .adjoint   ... flag for adjoint sensitivity analysis.\n'...
        '%%        true: on \n'...
        '%%        false: off (DEFAULT)\n'...
        '%%    .ism   ... only available for sensi_meth == 1. Method for computation of forward sensitivities.\n'...
        '%%        1: Simultaneous (DEFAULT)\n'...
        '%%        2: Staggered\n'...
        '%%        3: Staggered1\n'...
        '%%    .Nd   ... only available for sensi_meth == 2. Number of Interpolation nodes for forward solution. \n'...
        '%%              Default is 1000. \n'...
        '%%    .interpType   ... only available for sensi_meth == 2. Interpolation method for forward solution.\n'...
        '%%        1: Hermite (DEFAULT for problems without discontinuities)\n'...
        '%%        2: Polynomial (DEFAULT for problems with discontinuities)\n'...
        '%%    .data_model   ... noise model for data.\n'...
        '%%        1: Normal (DEFAULT)\n'...
        '%%        2: Lognormal \n'...
        '%%    .event_model   ... noise model for events.\n'...
        '%%        1: Normal (DEFAULT)\n'...
        '%%    .ordering   ... online state reordering.\n'...
        '%%        0: AMD reordering\n'...
        '%%        1: COLAMD reordering (default)\n'...
        '%%        2: natural reordering\n'...
        '%%\n'...
        '%% Outputs:\n'...
        '%% ========\n'...
        '%% sol.status ... flag for status of integration. generally status<0 for failed integration\n'...
        '%% sol.t ... vector at which the solution was computed\n'...
        '%% sol.llh ... likelihood value\n'...
        '%% sol.chi2 ... chi2 value\n'...
        '%% sol.sllh ... gradient of likelihood\n'...
        '%% sol.s2llh ... hessian of likelihood\n'...
        '%% sol.x ... time-resolved state vector\n'...
        '%% sol.y ... time-resolved output vector\n'...
        '%% sol.sx ... time-resolved state sensitivity vector\n'...
        '%% sol.sy ... time-resolved output sensitivity vector\n'...
        '%% sol.z event output\n'...
        '%% sol.sz sensitivity of event output\n'...
        ]);
    fprintf(fid,['function varargout = simulate_' this.modelname '(varargin)\n\n']);
    fprintf(fid,['%% DO NOT CHANGE ANYTHING IN THIS FILE UNLESS YOU ARE VERY SURE ABOUT WHAT YOU ARE DOING\n']);
    fprintf(fid,['%% MANUAL CHANGES TO THIS FILE CAN RESULT IN WRONG SOLUTIONS AND CRASHING OF MATLAB\n']);
    fprintf(fid,['if(nargin<2)\n']);
    fprintf(fid,['    error(''Not enough input arguments.'');\n']);
    fprintf(fid,['else\n']);
    fprintf(fid,['    tout=varargin{1};\n']);
    fprintf(fid,['    phi=varargin{2};\n']);
    fprintf(fid,['end\n']);
    
    fprintf(fid,['if(nargin>=3)\n']);
    fprintf(fid,['    kappa=varargin{3};\n']);
    fprintf(fid,['else\n']);
    fprintf(fid,['    kappa=[];\n']);
    fprintf(fid,['end\n']);
    switch(this.param)
        case 'log'
            fprintf(fid,'theta = exp(phi(:));\n\n');
        case 'log10'
            fprintf(fid,'theta = 10.^(phi(:));\n\n');
        case 'lin'
            fprintf(fid,'theta = phi(:);\n\n');
        otherwise
            disp('No valid parametrisation chosen! Valid options are "log","log10" and "lin". Using linear parametrisation (default)!')
            fprintf(fid,'theta = phi(:);\n\n');
    end
    fprintf(fid,'\n');
    if(nk==0)
        fprintf(fid,'if(nargin==2)\n');
        fprintf(fid,'    kappa = [];\n');
        fprintf(fid,'end\n');
    end
    
    fprintf(fid,['if(length(theta)<' num2str(np) ')\n']);
    fprintf(fid,'    error(''provided parameter vector is too short'');\n');
    fprintf(fid,'end\n');
    fprintf(fid,['if(length(kappa)<' num2str(nk) ')\n']);
    fprintf(fid,'    error(''provided constant vector is too short'');\n');
    fprintf(fid,'end\n');
    fprintf(fid,'\n');
    
    fprintf(fid,'\n');
    fprintf(fid,'pbar = ones(size(theta));\n');
    fprintf(fid,'pbar(pbar==0) = 1;\n');
    fprintf(fid,'xscale = [];\n');
    fprintf(fid,['if(nargin>=5)\n']);
    fprintf(fid,['    options_ami = amioption(varargin{5});\n']);
    fprintf(fid,['else\n']);
    fprintf(fid,['    options_ami = amioption();\n']);
    fprintf(fid,['end\n']);
    fprintf(fid,['if(isempty(options_ami.sens_ind))\n']);
    fprintf(fid,['    options_ami.sens_ind = 1:' num2str(np) ';\n']);
    fprintf(fid,['end\n']);
    fprintf(fid,['options_ami.id = transpose([' num2str(transpose(double(this.id))) ']);\n\n']);
    fprintf(fid,['options_ami.z2event = [' num2str(transpose(this.z2event)) ']; %% MUST NOT CHANGE THIS VALUE\n']);
    
    if(o2flag == 2)
        fprintf(fid,['if(nargin>=6)\n']);
        fprintf(fid,['    v = varargin{6};\n']);
        switch(this.param)
            case 'log'
                fprintf(fid,'    v = v(:).*theta(options_ami.sens_ind);\n');
            case 'log10'
                fprintf(fid,'    v = v(:).*theta(options_ami.sens_ind)*log(10);\n');
        end
        fprintf(fid,['    kappa = [kappa(:);v(:)];\n']);
        fprintf(fid,['else\n']);
        fprintf(fid,['    if(options_ami.sensi==2)\n']);
        fprintf(fid,['        error(''6th argument (multiplication vector is missing'');\n']);
        fprintf(fid,['    end\n']);
        fprintf(fid,['end\n']);
    end
    
    if(o2flag)
        fprintf(fid,['if(nargout>1)\n']);
        fprintf(fid,['    if(nargout>6)\n']);
        fprintf(fid,['        options_ami.sensi = 2;\n']);
        fprintf(fid,['        options_ami.sensi_meth = ''forward'';\n']);
        fprintf(fid,['    elseif(nargout>4)\n']);
        fprintf(fid,['        options_ami.sensi = 1;\n']);
        fprintf(fid,['        options_ami.sensi_meth = ''forward'';\n']);
        fprintf(fid,['    else\n']);
        fprintf(fid,['        options_ami.sensi = 0;\n']);
        fprintf(fid,['    end\n']);
        fprintf(fid,['end\n']);
    else
        fprintf(fid,['if(nargout>1)\n']);
        fprintf(fid,['    if(nargout>4)\n']);
        fprintf(fid,['        options_ami.sensi = 1;\n']);
        fprintf(fid,['        options_ami.sensi_meth = ''forward'';\n']);
        fprintf(fid,['    else\n']);
        fprintf(fid,['        options_ami.sensi = 0;\n']);
        fprintf(fid,['    end\n']);
        fprintf(fid,['end\n']);
    end
    fprintf(fid,['if(options_ami.ss>0)\n']);
    fprintf(fid,['    if(options_ami.sensi>1)\n']);
    fprintf(fid,['        error(''Computation of steady state sensitivity only possible for first order sensitivities'');\n']);
    fprintf(fid,['    end\n']);
    fprintf(fid,['    options_ami.sensi = 0;\n']);
    fprintf(fid,['end\n']);
    if(~this.forward)
        fprintf(fid,['if(options_ami.sensi>0)\n']);
        fprintf(fid,['    if(options_ami.sensi_meth == 1)\n']);
        fprintf(fid,['        error(''forward sensitivities are disabled as necessary routines were not compiled'');\n']);
        fprintf(fid,['    end\n']);
        fprintf(fid,['end\n']);
    end
    if(~this.adjoint)
        fprintf(fid,['if(options_ami.sensi>0)\n']);
        fprintf(fid,['    if(options_ami.sensi_meth == 2)\n']);
        fprintf(fid,['        error(''adjoint sensitivities are disabled as necessary routines were not compiled'');\n']);
        fprintf(fid,['    end\n']);
        fprintf(fid,['end\n']);
    end
    fprintf(fid,['np = length(options_ami.sens_ind); %% MUST NOT CHANGE THIS VALUE\n']);
    fprintf(fid,['if(np == 0)\n']);
    fprintf(fid,['    options_ami.sensi = 0;\n']);
    fprintf(fid,['end\n']);
    fprintf(fid,['if(isempty(options_ami.qpositivex))\n']);
    fprintf(fid,['    options_ami.qpositivex = zeros(' num2str(this.nx) ',1);\n']);
    fprintf(fid,['else\n']);
    fprintf(fid,['    if(numel(options_ami.qpositivex)==' num2str(this.nx) ')\n']);
    fprintf(fid,['        options_ami.qpositivex = options_ami.qpositivex(:);\n']);
    fprintf(fid,['    else\n']);
    fprintf(fid,['        error(''Number of elements in options_ami.qpositivex does not match number of states ' num2str(this.nx) ''');\n']);
    fprintf(fid,['    end\n']);
    fprintf(fid,['end\n']);
    fprintf(fid,'plist = options_ami.sens_ind-1;\n');
    fprintf(fid,['if(nargin>=4)\n']);
    fprintf(fid,['    if(isempty(varargin{4}));\n']);
    fprintf(fid,['        data=amidata(length(tout),' num2str(this.ny) ',' num2str(this.nz) ',options_ami.nmaxevent,length(kappa));\n']);
    fprintf(fid,['    else\n']);
    fprintf(fid,['        data=amidata(varargin{4});\n']);
    fprintf(fid,['    end\n']);
    fprintf(fid,['else\n']);
    fprintf(fid,['    data=amidata(length(tout),' num2str(this.ny) ',' num2str(this.nz) ',options_ami.nmaxevent,length(kappa));\n']);
    fprintf(fid,['end\n']);
    fprintf(fid,['if(data.ne>0);\n']);
    fprintf(fid,['    options_ami.nmaxevent = data.ne;\n']);
    fprintf(fid,['else\n']);
    fprintf(fid,['    data.ne = options_ami.nmaxevent;\n']);
    fprintf(fid,['end\n']);
    fprintf(fid,['if(isempty(kappa))\n']);
    fprintf(fid,['    kappa = data.condition;\n']);
    fprintf(fid,['end\n']);
    fprintf(fid,['if(max(options_ami.sens_ind)>' num2str(np) ')\n']);
    fprintf(fid,['    error(''Sensitivity index exceeds parameter dimension!'')\n']);
    fprintf(fid,['end\n']);
    
    switch(this.param)
        case 'log'
            fprintf(fid,['if(~isempty(options_ami.sx0))\n']);
            fprintf(fid,['    if(size(options_ami.sx0,2)~=np)\n']);
            fprintf(fid,['        error(''Number of rows in sx0 field does not agree with number of model parameters!'');\n']);
            fprintf(fid,['    end\n']);
            fprintf(fid,['    options_ami.sx0 = bsxfun(@times,options_ami.sx0,1./permute(theta(options_ami.sens_ind),[2,1]));\n']);
            fprintf(fid,['end\n']);
        case 'log10'
            fprintf(fid,['if(~isempty(options_ami.sx0))\n']);
            fprintf(fid,['    if(size(options_ami.sx0,2)~=np)\n']);
            fprintf(fid,['        error(''Number of rows in sx0 field does not agree with number of model parameters!'');\n']);
            fprintf(fid,['    end\n']);
            fprintf(fid,['    options_ami.sx0 = bsxfun(@times,options_ami.sx0,1./(permute(theta(options_ami.sens_ind),[2,1])*log(10)));\n']);
            fprintf(fid,['end\n']);
        otherwise
            fprintf(fid,['if(~isempty(options_ami.sx0))\n']);
            fprintf(fid,['    if(size(options_ami.sx0,2)~=np)\n']);
            fprintf(fid,['        error(''Number of rows in sx0 field does not agree with number of model parameters!'');\n']);
            fprintf(fid,['    end\n']);
            fprintf(fid,['    options_ami.sx0 = options_ami.sx0;\n']);
            fprintf(fid,['end\n']);
    end
    
    
    if(o2flag)
        fprintf(fid,'if(options_ami.sensi<2)\n');
        fprintf(fid,['sol = ami_' this.modelname '(tout,theta(1:' num2str(np) '),kappa(1:' num2str(nk) '),options_ami,plist,pbar,xscale,data);\n']);
        fprintf(fid,'else\n');
        switch(o2flag)
            case 1
                fprintf(fid,['sol = ami_' this.modelname '_o2(tout,theta(1:' num2str(np) '),kappa(1:' num2str(nk) '),options_ami,plist,pbar,xscale,data);\n']);
            case 2
                fprintf(fid,['sol = ami_' this.modelname '_o2vec(tout,theta(1:' num2str(np) '),kappa(1:' num2str(amimodelo2.nk) '),options_ami,plist,pbar,xscale,data);\n']);
        end
        fprintf(fid,'end\n');
    else
        fprintf(fid,['sol = ami_' this.modelname '(tout,theta(1:' num2str(np) '),kappa(1:' num2str(nk) '),options_ami,plist,pbar,xscale,data);\n']);
    end
    fprintf(fid,'if(options_ami.sensi==1)\n');
    switch(this.param)
        case 'log'
            fprintf(fid,['    sol.sllh = sol.llhS.*theta(options_ami.sens_ind);\n']);
            fprintf(fid,['    sol.sx = bsxfun(@times,sol.xS,permute(theta(options_ami.sens_ind),[3,2,1]));\n']);
            fprintf(fid,['    sol.sy = bsxfun(@times,sol.yS,permute(theta(options_ami.sens_ind),[3,2,1]));\n']);
            fprintf(fid,['    sol.sz = bsxfun(@times,sol.zS,permute(theta(options_ami.sens_ind),[3,2,1]));\n']);
        case 'log10'
            fprintf(fid,['    sol.sllh = sol.llhS.*theta(options_ami.sens_ind)*log(10);\n']);
            fprintf(fid,['    sol.sx = bsxfun(@times,sol.xS,permute(theta(options_ami.sens_ind),[3,2,1])*log(10));\n']);
            fprintf(fid,['    sol.sy = bsxfun(@times,sol.yS,permute(theta(options_ami.sens_ind),[3,2,1])*log(10));\n']);
            fprintf(fid,['    sol.sz = bsxfun(@times,sol.zS,permute(theta(options_ami.sens_ind),[3,2,1])*log(10));\n']);
        case 'lin'
            fprintf(fid,['    sol.sllh = sol.llhS;\n']);
            fprintf(fid,'    sol.sx = sol.xS;\n');
            fprintf(fid,'    sol.sy = sol.yS;\n');
            fprintf(fid,'    sol.sz = sol.zS;\n');
        otherwise
            fprintf(fid,['    sol.sllh = sol.llhS;\n']);
            fprintf(fid,'    sol.sx = sol.xS;\n');
            fprintf(fid,'    sol.sy = sol.yS;\n');
            fprintf(fid,'    sol.sz = sol.zS;\n');
    end
    fprintf(fid,'end\n');
    if(o2flag)
        fprintf(fid,'if(options_ami.sensi == 2)\n');
        fprintf(fid,['    sx = sol.xS(:,1:' num2str(nxtrue) ',:);\n']);
        fprintf(fid,['    sy = sol.yS(:,1:' num2str(nytrue) ',:);\n']);
        fprintf(fid,['    for iz = 1:' num2str(nztrue) '\n']);
        fprintf(fid,['        sz(:,iz,:) = sol.zS(:,2*iz-1,:);\n']);
        fprintf(fid,['    end\n']);
        switch(o2flag)
            case 1
                fprintf(fid,['    s2x = reshape(sol.xS(:,' num2str(nxtrue+1) ':end,:),length(tout),' num2str(nxtrue) ',length(theta(options_ami.sens_ind)),length(theta(options_ami.sens_ind)));\n']);
                fprintf(fid,['    s2y = reshape(sol.yS(:,' num2str(nytrue+1) ':end,:),length(tout),' num2str(nytrue) ',length(theta(options_ami.sens_ind)),length(theta(options_ami.sens_ind)));\n']);
            case 2
                fprintf(fid,['    s2x = sol.xS(:,' num2str(nxtrue+1) ':end,:);\n']);
                fprintf(fid,['    s2y = sol.yS(:,' num2str(nytrue+1) ':end,:);\n']);
        end
        fprintf(fid,['    for iz = 1:' num2str(nztrue) '\n']);
        switch(o2flag)
            case 1
                fprintf(fid,['        s2z(:,iz,:,:) = reshape(sol.zS(:,((iz-1)*(length(theta(options_ami.sens_ind)+1))+2):((iz-1)*(length(theta(options_ami.sens_ind)+1))+length(theta(options_ami.sens_ind))+1),:),options_ami.nmaxevent,1,length(theta(options_ami.sens_ind)),length(theta(options_ami.sens_ind)));\n']);
            case 2
                fprintf(fid,['        s2z= [];\n']);
                %TBD
        end
        fprintf(fid,['    end\n']);
        fprintf(fid,['    sol.x = sol.x(:,1:' num2str(nxtrue) ');\n']);
        fprintf(fid,['    sol.y = sol.y(:,1:' num2str(nytrue) ');\n']);
        fprintf(fid,['    sol.z = sol.z(:,1:' num2str(nztrue) ');\n']);
        switch(this.param)
            case 'log'
                fprintf(fid,['    sol.sx = bsxfun(@times,sx,permute(theta(options_ami.sens_ind),[3,2,1]));\n']);
                fprintf(fid,['    sol.sy = bsxfun(@times,sy,permute(theta(options_ami.sens_ind),[3,2,1]));\n']);
                if(nztrue>0)
                    fprintf(fid,['    sol.sz = bsxfun(@times,sz,permute(theta(options_ami.sens_ind),[3,2,1]));\n']);
                end
                switch(o2flag)
                    case 1
                        fprintf(fid,['    sol.s2x = bsxfun(@times,s2x,permute(theta(options_ami.sens_ind)*transpose(theta(options_ami.sens_ind)),[4,3,2,1])) + bsxfun(@times,sx,permute(diag(theta(options_ami.sens_ind).*ones(length(theta(options_ami.sens_ind)),1)),[4,3,2,1]));\n']);
                        fprintf(fid,['    sol.s2y = bsxfun(@times,s2y,permute(theta(options_ami.sens_ind)*transpose(theta(options_ami.sens_ind)),[4,3,2,1])) + bsxfun(@times,sy,permute(diag(theta(options_ami.sens_ind).*ones(length(theta(options_ami.sens_ind)),1)),[4,3,2,1]));\n']);
                        %fprintf(fid,['    sol.s2z = bsxfun(@times,s2z,permute(theta(options_ami.sens_ind)*transpose(theta(options_ami.sens_ind)),[4,3,2,1])) + bsxfun(@times,sz,permute(diag(theta(options_ami.sens_ind).*ones(length(theta(options_ami.sens_ind)),1)),[4,3,2,1]));\n']);
                    case 2
                        fprintf(fid,['    sol.s2x = bsxfun(@times,s2x,permute(theta(options_ami.sens_ind),[3,2,1])) + bsxfun(@times,sx,permute(v,[3,2,1]));\n']);
                        fprintf(fid,['    sol.s2y = bsxfun(@times,s2y,permute(theta(options_ami.sens_ind),[3,2,1])) + bsxfun(@times,sy,permute(v,[3,2,1]));\n']);
                        %fprintf(fid,['    sol.s2z = bsxfun(@times,s2z,permute(theta(options_ami.sens_ind),[3,2,1])) + bsxfun(@times,sz,permute(theta(options_ami.sens_ind),[3,2,1]));\n']);
                end
            case 'log10'
                fprintf(fid,['    sol.sx = bsxfun(@times,sx,permute(theta(options_ami.sens_ind),[3,2,1])*log(10));\n']);
                fprintf(fid,['    sol.sy = bsxfun(@times,sy,permute(theta(options_ami.sens_ind),[3,2,1])*log(10));\n']);
                if(nztrue>0)
                    fprintf(fid,['    sol.sz = bsxfun(@times,sz,permute(theta(options_ami.sens_ind),[3,2,1])*log(10));\n']);
                end
                switch(o2flag)
                    case 1
                        fprintf(fid,['    sol.s2x = bsxfun(@times,s2x,permute(theta(options_ami.sens_ind)*transpose(theta(options_ami.sens_ind))*(log(10)^2),[4,3,2,1])) + bsxfun(@times,sx,permute(diag(log(10)^2*theta(options_ami.sens_ind).*ones(length(theta(options_ami.sens_ind)),1)),[4,3,2,1]));\n']);
                        fprintf(fid,['    sol.s2y = bsxfun(@times,s2y,permute(theta(options_ami.sens_ind)*transpose(theta(options_ami.sens_ind))*(log(10)^2),[4,3,2,1])) + bsxfun(@times,sy,permute(diag(log(10)^2*theta(options_ami.sens_ind).*ones(length(theta(options_ami.sens_ind)),1)),[4,3,2,1]));\n']);
                        %fprintf(fid,['    sol.s2z = bsxfun(@times,s2z,permute(theta(options_ami.sens_ind)*transpose(theta(options_ami.sens_ind))*(log(10)^2),[4,3,2,1])) + bsxfun(@times,sz,permute(diag(log(10)^2*theta(options_ami.sens_ind).*ones(length(theta(options_ami.sens_ind)),1)),[4,3,2,1]));\n']);
                    case 2
                        fprintf(fid,['    sol.s2x = bsxfun(@times,s2x,permute(theta(options_ami.sens_ind),[3,2,1])*log(10)) + bsxfun(@times,sx,permute(v,[3,2,1])*log(10));\n']);
                        fprintf(fid,['    sol.s2y = bsxfun(@times,s2y,permute(theta(options_ami.sens_ind),[3,2,1])*log(10)) + bsxfun(@times,sy,permute(v,[3,2,1])*log(10));\n']);
                        %fprintf(fid,['    sol.s2z = bsxfun(@times,s2z,permute(theta(options_ami.sens_ind),[3,2,1])*log(10)) + bsxfun(@times,sz,permute(theta(options_ami.sens_ind),[3,2,1])*log(10)*log(10));\n']);
                end
            case 'lin'
                fprintf(fid,'    sol.sx = sx;\n');
                fprintf(fid,'    sol.s2x = s2x;\n');
                fprintf(fid,'    sol.sy = sy;\n');
                fprintf(fid,'    sol.s2y = s2y;\n');
                if(nztrue>0)
                    fprintf(fid,'    sol.sz = sz;\n');
                    fprintf(fid,'    sol.s2z = s2z;\n');
                end
            otherwise
                fprintf(fid,'    sol.sx = sx;\n');
                fprintf(fid,'    sol.s2x = s2x;\n');
                fprintf(fid,'    sol.sy = sy;\n');
                fprintf(fid,'    sol.s2y = s2y;\n');
                if(nztrue>0)
                    fprintf(fid,'    sol.sz = sz;\n');
                    fprintf(fid,'    sol.s2z = s2z;\n');
                end
        end
        fprintf(fid,'end\n');
    end
    
    fprintf(fid,['if(options_ami.sensi_meth == 3)\n']);
    switch(this.param)
        case 'log'
            fprintf(fid,['    sol.dxdotdp = bsxfun(@times,sol.dxdotdp,permute(theta(options_ami.sens_ind),[2,1]));\n']);
            fprintf(fid,['    sol.dydp = bsxfun(@times,sol.dydp,permute(theta(options_ami.sens_ind),[2,1]));\n']);
        case 'log10'
            fprintf(fid,['    sol.dxdotdp = bsxfun(@times,sol.dxdotdp,permute(theta(options_ami.sens_ind),[2,1])*log(10));\n']);
            fprintf(fid,['    sol.dydp = bsxfun(@times,sol.dydp,permute(theta(options_ami.sens_ind),[2,1])*log(10));\n']);
        otherwise
    end
    
    fprintf(fid,['    sol.sx = -sol.J\\sol.dxdotdp;\n']);
    fprintf(fid,['    sol.sy = sol.dydx*sol.sx + sol.dydp;\n']);
    
    
    fprintf(fid,['end\n']);
    fprintf(fid,['if(nargout>1)\n']);
    fprintf(fid,['    varargout{1} = sol.status;\n']);
    fprintf(fid,['    varargout{2} = sol.t;\n']);
    fprintf(fid,['    varargout{3} = sol.x;\n']);
    fprintf(fid,['    varargout{4} = sol.y;\n']);
    fprintf(fid,['    if(nargout>4)\n']);
    fprintf(fid,['        varargout{5} = sol.sx;\n']);
    fprintf(fid,['        varargout{6} = sol.sy;\n']);
    if(o2flag)
        fprintf(fid,['        if(nargout>6)\n']);
        fprintf(fid,['            varargout{7} = sol.s2x;\n']);
        fprintf(fid,['            varargout{8} = sol.s2y;\n']);
        fprintf(fid,['        end\n']);
    end
    fprintf(fid,['    end\n']);
    fprintf(fid,['else\n']);
    fprintf(fid,['    varargout{1} = sol;\n']);
    fprintf(fid,['end\n']);
    fprintf(fid,'end\n');
    
    fclose(fid);
end

