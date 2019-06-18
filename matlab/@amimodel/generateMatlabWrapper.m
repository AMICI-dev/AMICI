function generateMatlabWrapper(nx, ny, np, nk, nz, o2flag, amimodelo2, wrapperFilename, modelname, ...
    pscale, forward, adjoint)
    % generateMatlabWrapper generates the matlab wrapper for the compiled C files.
    %
    % Parameters:
    %  nx: number of states
    %  ny: number of observables
    %  np: number of parameters
    %  nk: number of fixed parameters
    %  nz: number of events
    %  o2flag: o2flag
    %  amimodelo2: this struct must contain all necessary symbolic
    %  definitions for second order sensivities @type amimodel
    %  wrapperFilename: output filename
    %  modelname: name of the model
    %  pscale: default parameter scaling
    %  forward: has forward sensitivity equations
    %  adjoint: has adjoint sensitivity equations
    %
    % Return values:
    %  void
    amiciRootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));

    if(~isempty(amimodelo2))
        nztrue = amimodelo2.nztrue;
        nxtrue = amimodelo2.nxtrue;
        nytrue = amimodelo2.nytrue;
        o2flag = amimodelo2.o2flag;
    else
        nztrue = nz;
        nxtrue = nx;
        nytrue = ny;
        o2flag = o2flag;
    end
    
    [commit_hash,branch,url] = getCommitHash(amiciRootDir);
    if(isempty(commit_hash))
        commit_hash = '#';
    end
    
    if(o2flag)
        nxtrue = amimodelo2.nxtrue;
        nytrue = amimodelo2.nytrue;
    end
    
    
    %% Generation of the simulation file
    
    fid = fopen(wrapperFilename,'w');
    fprintf(fid,['%% simulate_' modelname '.m is the matlab interface to the cvodes mex\n'...
        '%%   which simulates the ordinary differential equation and respective\n'...
        '%%   sensitivities according to user specifications.\n'...
        '%%   this routine was generated using AMICI commit ' commit_hash ' in branch ' branch ' in repo ' url '.\n'...
        '%%\n'...
        '%% USAGE:\n'...
        '%% ======\n'...
        '%% [...] = simulate_' modelname '(tout,theta)\n'...
        '%% [...] = simulate_' modelname '(tout,theta,kappa,data,options)\n'...
        '%% [status,tout,x,y,sx,sy] = simulate_' modelname '(...)\n'...
        '%%\n'...
        '%% INPUTS:\n'...
        '%% =======\n'...
        '%% tout ... 1 dimensional vector of timepoints at which a solution to the ODE is desired\n'...
        '%% theta ... 1 dimensional parameter vector of parameters for which sensitivities are desired.\n'...
        '%%           this corresponds to the specification in model.sym.p\n'...
        '%% kappa ... 1 dimensional parameter vector of parameters for which sensitivities are not desired.\n'...
        '%%           this corresponds to the specification in model.sym.k\n'...
        '%% data ... struct containing the following fields:\n'...
        '%%     Y ... 2 dimensional matrix containing data.\n'...
        '%%           columns must correspond to observables and rows to time-points\n'...
        '%%     Sigma_Y ... 2 dimensional matrix containing standard deviation of data.\n'...
        '%%           columns must correspond to observables and rows to time-points\n'...
        '%%     T ... (optional) 2 dimensional matrix containing events.\n'...
        '%%           columns must correspond to event-types and rows to possible event-times\n'...
        '%%     Sigma_T ... (optional) 2 dimensional matrix containing standard deviation of events.\n'...
        '%%           columns must correspond to event-types and rows to possible event-times\n'...
        '%% options ... additional options to pass to the cvodes solver. Refer to the cvodes guide for more documentation.\n'...
        '%%    .atol ... absolute tolerance for the solver. default is specified in the user-provided syms function.\n'...
        '%%        default value is 1e-16\n'...
        '%%    .rtol ... relative tolerance for the solver. default is specified in the user-provided syms function.\n'...
        '%%        default value is 1e-8\n'...
        '%%    .maxsteps    ... maximal number of integration steps. default is specified in the user-provided syms function.\n'...
        '%%        default value is 1e4\n'...
        '%%    .tstart    ... start of integration. for all timepoints before this, values will be set to initial value.\n'...
        '%%        default value is 0\n'...
        '%%    .sens_ind ... 1 dimensional vector of indexes for which sensitivities must be computed.\n'...
        '%%        default value is 1:length(theta).\n'...
        '%%    .x0 ... user-provided state initialisation. This should be a vector of dimension [#states, 1].\n'...
        '%%        default is state initialisation based on the model definition.\n'...
        '%%    .sx0 ... user-provided sensitivity initialisation. this should be a matrix of dimension [#states x #parameters].\n'...
        '%%        default is sensitivity initialisation based on the derivative of the state initialisation.\n'...
        '%%    .lmm    ... linear multistep method for forward problem.\n'...
        '%%        1: Adams-Bashford\n'...
        '%%        2: BDF (DEFAULT)\n'...
        '%%    .iter    ... iteration method for linear multistep.\n'...
        '%%        1: Functional\n'...
        '%%        2: Newton (DEFAULT)\n'...
        '%%    .linsol   ... linear solver module.\n'...
        '%%        direct solvers:\n'...
        '%%        1: Dense\n'...
        '%%        2: Band (not implemented)\n'...
        '%%        3: LAPACK Dense (not implemented)\n'...
        '%%        4: LAPACK Band  (not implemented)\n'...
        '%%        5: Diag (not implemented)\n'...
        '%%        implicit krylov solvers:\n'...
        '%%        6: SPGMR\n'...
        '%%        7: SPBCG\n'...
        '%%        8: SPTFQMR\n'...
        '%%        sparse solvers:\n'...
        '%%        9: KLU (DEFAULT)\n'...
        '%%    .stldet   ... flag for stability limit detection. this should be turned on for stiff problems.\n'...
        '%%        0: OFF\n'...
        '%%        1: ON (DEFAULT)\n'...
        '%%    .sensi   ... sensitivity order.\n'...
        '%%        0: OFF (DEFAULT)\n'...
        '%%        1: first\n'...
        '%%        2: second\n'...
        '%%    .sensi_meth   ... method for sensitivity analysis.\n'...
        '%%        0: no sensitivity analysis\n'...
        '%%        1 or ''forward'': forward sensitivity analysis (DEFAULT)\n'...
        '%%        2 or ''adjoint'': adjoint sensitivity analysis \n'...
        '%%        3 or ''ss'': defined but not documented \n'...
        '%%    .adjoint   ... flag for adjoint sensitivity analysis.\n'...
        '%%        true: on \n'...
        '%%        false: off (DEFAULT)\n'...
        '%%        NO LONGER USED: Replaced by sensi_meth\n'...
        '%%    .ism   ... only available for sensi_meth == 1. Method for computation of forward sensitivities.\n'...
        '%%        1: Simultaneous (DEFAULT)\n'...
        '%%        2: Staggered\n'...
        '%%        3: Staggered1\n'...
        '%%    .Nd   ... only available for sensi_meth == 2. Number of Interpolation nodes for forward solution. \n'...
        '%%        default is 1000. \n'...
        '%%    .interpType   ... only available for sensi_meth == 2. Interpolation method for forward solution.\n'...
        '%%        1: Hermite (DEFAULT for problems without discontinuities)\n'...
        '%%        2: Polynomial (DEFAULT for problems with discontinuities)\n'...
        '%%    .ordering   ... online state reordering.\n'...
        '%%        0: AMD reordering (default)\n'...
        '%%        1: COLAMD reordering\n'...
        '%%        2: natural reordering\n'...
        '%%    .newton_maxsteps   ... maximum newton steps\n'...
        '%%        default value is 40\n'...
        '%%        a value of 0 will disable the newton solver\n'...
        '%%    .newton_maxlinsteps   ... maximum linear steps\n'...
        '%%        default value is 100\n'...
        '%%    .newton_preeq   ... preequilibration of system via newton solver\n'...
        '%%        default value is false\n'...
        '%%    .pscale   ... parameter scaling\n'...
        '%%        []: (DEFAULT) use value specified in the model (fallback: ''lin'')\n'...
        '%%        0 or ''lin'': linear\n'...
        '%%        1 or ''log'': natural log (base e)\n'...
        '%%        2 or ''log10'': log to the base 10\n'...
        '%%\n'...
        '%% Outputs:\n'...
        '%% ========\n'...
        '%% sol.status ... flag for status of integration. generally status<0 for failed integration\n'...
        '%% sol.t ... vector at which the solution was computed\n'...
        '%% sol.llh ... likelihood value\n'...
        '%% sol.chi2 ... chi2 value\n'...
        '%% sol.sllh ... gradient of likelihood\n'...
        '%% sol.s2llh ... hessian or hessian-vector-product of likelihood\n'...
        '%% sol.x ... time-resolved state vector\n'...
        '%% sol.y ... time-resolved output vector\n'...
        '%% sol.sx ... time-resolved state sensitivity vector\n'...
        '%% sol.sy ... time-resolved output sensitivity vector\n'...
        '%% sol.z ... event output\n'...
        '%% sol.sz ... sensitivity of event output\n'...
        ]);
    
    
    fprintf(fid,['function varargout = simulate_' modelname '(varargin)\n\n']);
    fprintf(fid,'%% DO NOT CHANGE ANYTHING IN THIS FILE UNLESS YOU ARE VERY SURE ABOUT WHAT YOU ARE DOING\n');
    fprintf(fid,'%% MANUAL CHANGES TO THIS FILE CAN RESULT IN WRONG SOLUTIONS AND CRASHING OF MATLAB\n');
    fprintf(fid,'if(nargin<2)\n');
    fprintf(fid,'    error(''Not enough input arguments.'');\n');
    fprintf(fid,'else\n');
    fprintf(fid,'    tout=varargin{1};\n');
    fprintf(fid,'    theta=varargin{2};\n');
    fprintf(fid,'end\n');
    
    fprintf(fid,'if(nargin>=3)\n');
    fprintf(fid,'    kappa=varargin{3};\n');
    fprintf(fid,'else\n');
    fprintf(fid,'    kappa=[];\n');
    fprintf(fid,'end\n');
    fprintf(fid,'\n');
    if(nk==0)
        fprintf(fid,'if(nargin==2)\n');
        fprintf(fid,'    kappa = [];\n');
        fprintf(fid,'end\n');
    end
    
    fprintf(fid,['if(length(theta)<' num2str(np) ')\n']);
    fprintf(fid,'    error(''provided parameter vector is too short'');\n');
    fprintf(fid,'end\n');
    fprintf(fid,'\n');
    
    fprintf(fid,'\n');
    fprintf(fid,'xscale = [];\n');
    fprintf(fid,'if(nargin>=5)\n');
    fprintf(fid,'    if(isa(varargin{5},''amioption''))\n');
    fprintf(fid,'        options_ami = varargin{5};\n');
    fprintf(fid,'    else\n');
    fprintf(fid,'        options_ami = amioption(varargin{5});\n');
    fprintf(fid,'    end\n');
    fprintf(fid,'else\n');
    fprintf(fid,'    options_ami = amioption();\n');
    fprintf(fid,'end\n');
    fprintf(fid,'if(isempty(options_ami.sens_ind))\n');
    fprintf(fid,['    options_ami.sens_ind = 1:' num2str(np) ';\n']);
    fprintf(fid,['end\n']);
    if(o2flag == 0)
        fprintf(fid,'if(options_ami.sensi>1)\n');
        fprintf(fid,'    error(''Second order sensitivities were requested but not computed'');\n');
        fprintf(fid,'end\n');
    end
    fprintf(fid,'\n');
    
    
    fprintf(fid,'if(isempty(options_ami.pscale))\n');
    fprintf(fid,['    options_ami.pscale = ''' pscale ''' ;\n']);
    fprintf(fid,'end\n');
    
    if(o2flag == 2)
        fprintf(fid,'if(nargin>=6)\n');
        fprintf(fid,'    chainRuleFactor = getChainRuleFactors(options_ami.pscale, theta, options_ami.sens_ind);\n');
        fprintf(fid,'    v = varargin{6};\n');
        fprintf(fid,'    v = v(:).*chainRuleFactor(:);\n');
        fprintf(fid,'else\n');
        fprintf(fid,'    if(options_ami.sensi==2)\n');
        fprintf(fid,'        error(''6th argument (multiplication vector is missing'');\n');
        fprintf(fid,'    end\n');
        fprintf(fid,'end\n');
    end
    
    if(o2flag)
        fprintf(fid,'if(nargout>1)\n');
        fprintf(fid,'    if(nargout>6)\n');
        fprintf(fid,'        options_ami.sensi = 2;\n');
        fprintf(fid,'        options_ami.sensi_meth = ''forward'';\n');
        fprintf(fid,'    elseif(nargout>4)\n');
        fprintf(fid,'        options_ami.sensi = 1;\n');
        fprintf(fid,'        options_ami.sensi_meth = ''forward'';\n');
        fprintf(fid,'    else\n');
        fprintf(fid,'        options_ami.sensi = 0;\n');
        fprintf(fid,'    end\n');
        fprintf(fid,'end\n');
    else
        fprintf(fid,'if(nargout>1)\n');
        fprintf(fid,'    if(nargout>4)\n');
        fprintf(fid,'        options_ami.sensi = 1;\n');
        fprintf(fid,'        options_ami.sensi_meth = ''forward'';\n');
        fprintf(fid,'    else\n');
        fprintf(fid,'        options_ami.sensi = 0;\n');
        fprintf(fid,'    end\n');
        fprintf(fid,'end\n');
    end
    if(~forward)
        fprintf(fid,'if(options_ami.sensi>0)\n');
        fprintf(fid,'    if(options_ami.sensi_meth == 1)\n');
        fprintf(fid,'        error(''forward sensitivities are disabled as necessary routines were not compiled'');\n');
        fprintf(fid,'    end\n');
        fprintf(fid,'end\n');
    end
    if(~adjoint)
        fprintf(fid,'if(options_ami.sensi>0)\n');
        fprintf(fid,'    if(options_ami.sensi_meth == 2)\n');
        fprintf(fid,'        error(''adjoint sensitivities are disabled as necessary routines were not compiled'');\n');
        fprintf(fid,'    end\n');
        fprintf(fid,'end\n');
    end
    fprintf(fid,['nplist = length(options_ami.sens_ind); %% MUST NOT CHANGE THIS VALUE\n']);
    fprintf(fid,['if(nplist == 0)\n']);
    fprintf(fid,['    options_ami.sensi = 0;\n']);
    fprintf(fid,['end\n']);
    if(o2flag)
        fprintf(fid,['if(options_ami.sensi > 1)\n']);
        fprintf(fid,['    nxfull = ' num2str(amimodelo2.nx) ';\n']);
        fprintf(fid,['else\n']);
        fprintf(fid,['    nxfull = ' num2str(nxtrue) ';\n']);
        fprintf(fid,['end\n']);
    else
        fprintf(fid,['nxfull = ' num2str(nx) ';\n']);
    end
    fprintf(fid,'plist = options_ami.sens_ind-1;\n');
    fprintf(fid,['if(nargin>=4)\n']);
    fprintf(fid,['    if(isempty(varargin{4}))\n']);
    fprintf(fid,['        data=[];\n']);
    fprintf(fid,['    else\n']);
    fprintf(fid,['        if(isa(varargin{4},''amidata''))\n']);
    fprintf(fid,['             data=varargin{4};\n']);
    fprintf(fid,['        else\n']);
    fprintf(fid,['            data=amidata(varargin{4});\n']);
    fprintf(fid,['        end\n']);
    fprintf(fid,['        if(data.ne>0)\n']);
    fprintf(fid,['            options_ami.nmaxevent = data.ne;\n']);
    fprintf(fid,['        else\n']);
    fprintf(fid,['            data.ne = options_ami.nmaxevent;\n']);
    fprintf(fid,['        end\n']);
    fprintf(fid,['        if(isempty(kappa))\n']);
    fprintf(fid,['            kappa = data.condition;\n']);
    fprintf(fid,['        end\n']);
    fprintf(fid,['        if(isempty(tout))\n']);
    fprintf(fid,['            tout = data.t;\n']);
    fprintf(fid,['        end\n']);
    fprintf(fid,['    end\n']);
    fprintf(fid,['else\n']);
    fprintf(fid,['    data=[];\n']);
    fprintf(fid,['end\n']);
    fprintf(fid,['if(~all(tout==sort(tout)))\n']);
    fprintf(fid,['    error(''Provided time vector is not monotonically increasing!'');\n']);
    fprintf(fid,['end\n']);
    fprintf(fid,['if(max(options_ami.sens_ind)>' num2str(np) ')\n']);
    fprintf(fid,['    error(''Sensitivity index exceeds parameter dimension!'')\n']);
    fprintf(fid,['end\n']);
    fprintf(fid,['if(length(kappa)<' num2str(nk) ')\n']);
    fprintf(fid,'    error(''provided condition vector is too short'');\n');
    fprintf(fid,'end\n');
    
    if(o2flag == 2)
        fprintf(fid,'if(nargin>=6)\n');
        fprintf(fid,'    kappa = [kappa(:);v(:)];\n');
        fprintf(fid,'end\n');
    end
    
    fprintf(fid,'init = struct();\n');
    fprintf(fid,'if(~isempty(options_ami.x0))\n');
    fprintf(fid,'    if(size(options_ami.x0,2)~=1)\n');
    fprintf(fid,'        error(''x0 field must be a column vector!'');\n');
    fprintf(fid,'    end\n');
    fprintf(fid,'    if(size(options_ami.x0,1)~=nxfull)\n');
    fprintf(fid,'        error(''Number of rows in x0 field does not agree with number of states!'');\n');
    fprintf(fid,'    end\n');
    fprintf(fid,'    init.x0 = options_ami.x0;\n');
    fprintf(fid,'end\n');
    fprintf(fid,'if(~isempty(options_ami.sx0))\n');
    fprintf(fid,'    if(size(options_ami.sx0,2)~=nplist)\n');
    fprintf(fid,'        error(''Number of columns in sx0 field does not agree with number of model parameters!'');\n');
    fprintf(fid,'    end\n');
    fprintf(fid,'    if(size(options_ami.sx0,1)~=nxfull)\n');
    fprintf(fid,'        error(''Number of rows in sx0 field does not agree with number of states!'');\n');
    fprintf(fid,'    end\n');
    fprintf(fid,'end\n');
    
    
    if(o2flag)
        fprintf(fid,'if(options_ami.sensi<2)\n');
        fprintf(fid,['    sol = ami_' modelname '(tout,theta(1:' num2str(np) '),kappa(1:' num2str(nk) '),options_ami,plist,xscale,init,data);\n']);
        fprintf(fid,'else\n');
        switch(o2flag)
            case 1
                fprintf(fid,['    sol = ami_' modelname '_o2(tout,theta(1:' num2str(np) '),kappa(1:' num2str(nk) '),options_ami,plist,xscale,init,data);\n']);
            case 2
                fprintf(fid,['    sol = ami_' modelname '_o2vec(tout,theta(1:' num2str(np) '),kappa(1:' num2str(amimodelo2.nk) '),options_ami,plist,xscale,init,data);\n']);
        end
        fprintf(fid,'end\n');
    else
        fprintf(fid,['sol = ami_' modelname '(tout,theta(1:' num2str(np) '),kappa(1:' num2str(nk) '),options_ami,plist,xscale,init,data);\n']);
    end
    
    if(o2flag)
        fprintf(fid,'if(options_ami.sensi == 2)\n');
        fprintf(fid, '    if(~(options_ami.sensi_meth==2))\n');
        
        fprintf(fid,['        sz = zeros(size(sol.z,1),' num2str(nztrue) ',length(theta(options_ami.sens_ind)));\n']);
        fprintf(fid,['        ssigmaz = zeros(size(sol.z,1),' num2str(nztrue) ',length(theta(options_ami.sens_ind)));\n']);
        fprintf(fid,['        srz = zeros(size(sol.z,1),' num2str(nztrue) ',length(theta(options_ami.sens_ind)));\n']);
        fprintf(fid,['        for iz = 1:' num2str(nztrue) '\n']);
        fprintf(fid,['            sz(:,iz,:) = sol.sz(:,2*iz-1,:);\n']);
        fprintf(fid,['            ssigmaz(:,iz,:) = sol.ssigmaz(:,2*iz-1,:);\n']);
        fprintf(fid,['            srz(:,iz,:) = sol.srz(:,2*iz-1,:);\n']);
        fprintf(fid,['        end\n']);
        switch(o2flag)
            case 1
                fprintf(fid,['        sol.s2x = reshape(sol.sx(:,' num2str(nxtrue+1) ':end,:),length(tout),' num2str(nxtrue) ',' num2str(np) ',length(options_ami.sens_ind));\n']);
                fprintf(fid,['        sol.s2y = reshape(sol.sy(:,' num2str(nytrue+1) ':end,:),length(tout),' num2str(nytrue) ',' num2str(np) ',length(options_ami.sens_ind));\n']);
                fprintf(fid,['        sol.s2sigmay = reshape(sol.ssigmay(:,' num2str(nytrue+1) ':end,:),length(tout),' num2str(nytrue) ',' num2str(np) ',length(options_ami.sens_ind));\n']);
            case 2
                fprintf(fid,['        sol.s2x = sol.sx(:,' num2str(nxtrue+1) ':end,:);\n']);
                fprintf(fid,['        sol.s2y = sol.sy(:,' num2str(nytrue+1) ':end,:);\n']);
                fprintf(fid,['        sol.s2sigmay = sol.ssigmay(:,' num2str(nytrue+1) ':end,:);\n']);
        end
        switch(o2flag)
            case 1
                fprintf(fid,['        s2z = zeros(size(sol.z,1),' num2str(nztrue) ',' num2str(np) ',length(options_ami.sens_ind));\n']);
                fprintf(fid,['        s2sigmaz = zeros(size(sol.z,1),' num2str(nztrue) ',' num2str(np) ',length(options_ami.sens_ind));\n']);
                fprintf(fid,['        s2rz = zeros(size(sol.z,1),' num2str(nztrue) ',' num2str(np) ',length(options_ami.sens_ind));\n']);
            case 2
                fprintf(fid,['        s2z = zeros(size(sol.z,1),' num2str(nztrue) ',length(theta(options_ami.sens_ind)));\n']);
                fprintf(fid,['        s2sigmaz = zeros(size(sol.z,1),' num2str(nztrue) ',length(theta(options_ami.sens_ind)));\n']);
                fprintf(fid,['        s2rz = zeros(size(sol.z,1),' num2str(nztrue) ',length(theta(options_ami.sens_ind)));\n']);
        end
        fprintf(fid,['        for iz = 1:' num2str(nztrue) '\n']);
        switch(o2flag)
            case 1
                fprintf(fid,['            sol.s2z(:,iz,:,:) = reshape(sol.sz(:,((iz-1)*(' num2str(np) '+1)+2):((iz-1)*(' num2str(np) '+1)+' num2str(np) '+1),:),options_ami.nmaxevent,1,' num2str(np) ',length(options_ami.sens_ind));\n']);
                fprintf(fid,['            sol.s2sigmaz(:,iz,:,:) = reshape(sol.ssigmaz(:,((iz-1)*(' num2str(np) '+1)+2):((iz-1)*(' num2str(np) '+1)+' num2str(np) '+1),:),options_ami.nmaxevent,1,' num2str(np) ',length(options_ami.sens_ind));\n']);
                fprintf(fid,['            sol.s2rz(:,iz,:,:) = reshape(sol.srz(:,((iz-1)*(' num2str(np) '+1)+2):((iz-1)*(' num2str(np) '+1)+' num2str(np) '+1),:),options_ami.nmaxevent,1,' num2str(np) ',length(options_ami.sens_ind));\n']);
            case 2
                fprintf(fid,['            sol.s2z(:,iz,:) = reshape(sol.sz(:,2*(iz-1)+2,:),options_ami.nmaxevent,1,length(theta(options_ami.sens_ind)));\n']);
                fprintf(fid,['            sol.s2sigmaz(:,iz,:) = reshape(sol.ssigmaz(:,2*(iz-1)+2,:),options_ami.nmaxevent,1,length(theta(options_ami.sens_ind)));\n']);
                fprintf(fid,['            sol.s2rz(:,iz,:) = reshape(sol.srz(:,2*(iz-1)+2,:),options_ami.nmaxevent,1,length(theta(options_ami.sens_ind)));\n']);
        end
        fprintf(fid,'        end\n');
        fprintf(fid,['        sol.sx = sol.sx(:,1:' num2str(nxtrue) ',:);\n']);
        fprintf(fid,['        sol.sy = sol.sy(:,1:' num2str(nytrue) ',:);\n']);
        fprintf(fid,['        sol.ssigmay = sol.ssigmay(:,1:' num2str(nytrue) ',:);\n']);
        fprintf(fid,['        if(iz>0)\n']);
        fprintf(fid,['            sol.sz = sz;\n']);
        fprintf(fid,['            sol.ssigmaz = ssigmaz;\n']);
        fprintf(fid,['            sol.srz = srz;\n']);
        fprintf(fid,'         end\n');
        fprintf(fid,'    end\n');
        fprintf(fid,['    sol.x = sol.x(:,1:' num2str(nxtrue) ');\n']);
        fprintf(fid,['    sol.y = sol.y(:,1:' num2str(nytrue) ');\n']);
        fprintf(fid,['    sol.sigmay = sol.sigmay(:,1:' num2str(nytrue) ');\n']);
        fprintf(fid,['    sol.z = sol.z(:,1:' num2str(nztrue) ');\n']);
        fprintf(fid,['    sol.rz = sol.rz(:,1:' num2str(nztrue) ');\n']);
        fprintf(fid,['    sol.sigmaz = sol.sigmaz(:,1:' num2str(nztrue) ');\n']);
        fprintf(fid,'end\n');
    end
    
    fprintf(fid,'if(nargout>1)\n');
    fprintf(fid,'    varargout{1} = sol.status;\n');
    fprintf(fid,'    varargout{2} = sol.t;\n');
    fprintf(fid,'    varargout{3} = sol.x;\n');
    fprintf(fid,'    varargout{4} = sol.y;\n');
    fprintf(fid,'    if(nargout>4)\n');
    fprintf(fid,'        varargout{5} = sol.sx;\n');
    fprintf(fid,'        varargout{6} = sol.sy;\n');
    if(o2flag)
        fprintf(fid,'        if(nargout>6)\n');
        fprintf(fid,'            varargout{7} = sol.s2x;\n');
        fprintf(fid,'            varargout{8} = sol.s2y;\n');
        fprintf(fid,'        end\n');
    end
    fprintf(fid,'    end\n');
    fprintf(fid,'else\n');
    fprintf(fid,'    varargout{1} = sol;\n');
    fprintf(fid,'end\n');
    
    fprintf(fid,'function chainRuleFactors = getChainRuleFactors(pscale, theta, sens_ind)\n');
    fprintf(fid,'    if(length(pscale) == 1 && length(sens_ind) ~= length(pscale))\n');
    fprintf(fid,'        chainRuleFactors = arrayfun(@(x, ip) getChainRuleFactor(x, theta(ip)), repmat(pscale, 1, length(sens_ind)), sens_ind);\n');
    fprintf(fid,'    else\n');
    fprintf(fid,'        chainRuleFactors = arrayfun(@(x, ip) getChainRuleFactor(x, theta(ip)), pscale, sens_ind);\n');
    fprintf(fid,'    end\n');
    fprintf(fid,'end\n');
    fprintf(fid,'\n');

    fprintf(fid,'function chainRuleFactor = getChainRuleFactor(pscale, parameterValue)\n');
    fprintf(fid,'    switch (pscale)\n');
    fprintf(fid,'        case 1\n');
    fprintf(fid,'            chainRuleFactor = exp(parameterValue);\n');
    fprintf(fid,'        case 2\n');
    fprintf(fid,'            chainRuleFactor = 10.^parameterValue*log(10);\n');
    fprintf(fid,'        otherwise\n');
    fprintf(fid,'            chainRuleFactor = 1.0;\n');
    fprintf(fid,'    end\n');
    fprintf(fid,'end\n');
    fprintf(fid,'\n');

    fprintf(fid,'end\n');
    
    fclose(fid);
    
end

