% simulate_villab2.m is the matlab interface to the cvodes mex
%   which simulates the ordinary differential equation and respective
%   sensitivities according to user specifications.
%
% USAGE:
% ======
% [...] = simulate_villab2(tout,theta)
% [...] = simulate_villab2(tout,theta,kappa,data,options)
% [status,tout,x,y,sx,sy] = simulate_villab2(...)
%
% INPUTS:
% =======
% tout ... 1 dimensional vector of timepoints at which a solution to the ODE is desired
% theta ... 1 dimensional parameter vector of parameters for which sensitivities are desired.
%           this corresponds to the specification in model.sym.p
% kappa ... 1 dimensional parameter vector of parameters for which sensitivities are not desired.
%           this corresponds to the specification in model.sym.k
% data ... struct containing the following fields. Can have the following fields %     Y ... 2 dimensional matrix containing data.
%           columns must correspond to observables and rows to time-points
%     Sigma_Y ... 2 dimensional matrix containing standard deviation of data.
%           columns must correspond to observables and rows to time-points
%     T ... (optional) 2 dimensional matrix containing events.
%           columns must correspond to event-types and rows to possible event-times
%     Sigma_T ... (optional) 2 dimensional matrix containing standard deviation of events.
%           columns must correspond to event-types and rows to possible event-times
% options ... additional options to pass to the cvodes solver. Refer to the cvodes guide for more documentation.
%    .cvodes_atol ... absolute tolerance for the solver. default is specified in the user-provided syms function.
%    .cvodes_rtol ... relative tolerance for the solver. default is specified in the user-provided syms function.
%    .cvodes_maxsteps    ... maximal number of integration steps. default is specified in the user-provided syms function.
%    .tstart    ... start of integration. for all timepoints before this, values will be set to initial value.
%    .sens_ind ... 1 dimensional vector of indexes for which sensitivities must be computed.
%           default value is 1:length(theta).
%    .sx0 ... user-provided sensitivity initialisation. this should be a matrix of dimension [#states x #parameters].
%        default is sensitivity initialisation based on the derivative of the state initialisation.%    .lmm    ... linear multistep method for forward problem.
%        1: Adams-Bashford
%        2: BDF (DEFAULT)
%    .iter    ... iteration method for linear multistep.
%        1: Functional
%        2: Newton (DEFAULT)
%    .linsol   ... linear solver module.
%        direct solvers:
%        1: Dense (DEFAULT)
%        2: Band (not implented)
%        3: LAPACK Dense (not implented)
%        4: LAPACK Band  (not implented)
%        5: Diag (not implented)
%        implicit krylov solvers:
%        6: SPGMR
%        7: SPBCG
%        8: SPTFQMR
%        sparse solvers:
%        9: KLU
%    .stldet   ... flag for stability limit detection. this should be turned on for stiff problems.
%        0: OFF
%        1: ON (DEFAULT)
%    .qPositiveX   ... vector of 0 or 1 of same dimension as state vector. 1 enforces positivity of states.
%    .sensi_meth   ... method for sensitivity analysis.
%        'forward': forward sensitivity analysis (DEFAULT)
%        'adjoint': adjoint sensitivity analysis 
%        'ss': steady state sensitivity analysis 
%    .adjoint   ... flag for adjoint sensitivity analysis.
%        true: on 
%        false: off (DEFAULT)
%    .ism   ... only available for sensi_meth == 1. Method for computation of forward sensitivities.
%        1: Simultaneous (DEFAULT)
%        2: Staggered
%        3: Staggered1
%    .Nd   ... only available for sensi_meth == 2. Number of Interpolation nodes for forward solution. 
%              Default is 1000. 
%    .interpType   ... only available for sensi_meth == 2. Interpolation method for forward solution.
%        1: Hermite (DEFAULT for problems without discontinuities)
%        2: Polynomial (DEFAULT for problems with discontinuities)
%    .data_model   ... noise model for data.
%        1: Normal (DEFAULT)
%        2: Lognormal 
%    .event_model   ... noise model for events.
%        1: Normal (DEFAULT)
%    .ordering   ... online state reordering.
%        0: AMD reordering
%        1: COLAMD reordering (default)
%        2: natural reordering
%
% Outputs:
% ========
% sol.status ... flag for status of integration. generally status<0 for failed integration
% sol.tout ... vector at which the solution was computed
% sol.llh ... likelihood value
% sol.chi2 ... chi2 value
% sol.sllh ... gradient of likelihood
% sol.s2llh ... hessian of likelihood
% sol.x ... time-resolved state vector
% sol.y ... time-resolved output vector
% sol.sx ... time-resolved state sensitivity vector
% sol.sy ... time-resolved output sensitivity vector
% sol.xdot time-resolved right-hand side of differential equation
% sol.z event output
% sol.sz sensitivity of event output
function varargout = simulate_villab2(varargin)

% DO NOT CHANGE ANYTHING IN THIS FILE UNLESS YOU ARE VERY SURE ABOUT WHAT YOU ARE DOING
% MANUAL CHANGES TO THIS FILE CAN RESULT IN WRONG SOLUTIONS AND CRASHING OF MATLAB
if(nargin<2)
    error('Not enough input arguments.');
else
    tout=varargin{1};
    phi=varargin{2};
end
if(nargin>=3)
    kappa=varargin{3};
else
    kappa=[];
end
theta = 10.^(phi(:));


if(nargin==2)
    kappa = [];
end
if(length(theta)<116)
    error('provided parameter vector is too short');
end
if(length(kappa)<0)
    error('provided constant vector is too short');
end

options_ami.atol = 1e-08;
options_ami.rtol = 1e-08;
options_ami.maxsteps = 10000;
options_ami.sens_ind = 1:116;
options_ami.id = transpose([0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]);

options_ami.ne = 0; % MUST NOT CHANGE THIS VALUE
options_ami.tstart = 0;
options_ami.lmm = 2;
options_ami.iter = 2;
options_ami.linsol = 9;
options_ami.stldet = 1;
options_ami.Nd = 1000;
options_ami.interpType = 1;
options_ami.lmmB = 2;
options_ami.iterB = 2;
options_ami.ism = 1;
options_ami.sensi_meth = 'forward';

options_ami.sensi = 0;

options_ami.nmaxevent = 10;
options_ami.z2event = []; % MUST NOT CHANGE THIS VALUE
options_ami.ubw = 15; % MUST NOT CHANGE THIS VALUE
options_ami.lbw = 15; % MUST NOT CHANGE THIS VALUE

options_ami.data_model = 1;
options_ami.event_model = 1;

options_ami.ordering = 1;

options_ami.ss = 0;

sol.status = 0;
sol.llh = 0;
sol.chi2 = 0;
sol.t = tout;
sol.numsteps = zeros(length(tout),1);
sol.numrhsevals = zeros(length(tout),1);
sol.order = zeros(length(tout),1);
sol.numstepsS = zeros(length(tout),1);
sol.numrhsevalsS = zeros(length(tout),1);

pbar = ones(size(theta));
pbar(pbar==0) = 1;
xscale = [];
if(nargin>=5)
    options_ami = am_setdefault(varargin{5},options_ami);
else
end
sol.z = NaN(options_ami.nmaxevent,0);
if(nargout>1)
    if(nargout>4)
        options_ami.sensi = 1;
        options_ami.sensi_meth = 'forward';
    else
        options_ami.sensi = 0;
    end
end
if(ischar(options_ami.sensi_meth))
    if(strcmp(options_ami.sensi_meth,'forward'))
        options_ami.sensi_meth = 1;
    elseif(strcmp(options_ami.sensi_meth,'adjoint'))
        options_ami.sensi_meth = 2;
    elseif(strcmp(options_ami.sensi_meth,'ss'))
        options_ami.sensi_meth = 3;
        options_ami.sensi = 0;
    else
        error('Invalid choice of options.sensi_meth. Must be either ''forward'',''adjoint'' or ''ss''');
    end
else
    error('Invalid choice of options.sensi_meth. Must be either ''forward'',''adjoint'' or ''ss''');
end
if(options_ami.ss>0)
    if(options_ami.sensi>1)
        error('Computation of steady state sensitivity only possible for first order sensitivities');
    end
    options_ami.sensi = 0;
end
options_ami.np = length(options_ami.sens_ind); % MUST NOT CHANGE THIS VALUE
if(options_ami.np == 0)
    options_ami.sensi = 0;
end
options_ami.nx = 18; % MUST NOT CHANGE THIS VALUE
options_ami.ny = 9; % MUST NOT CHANGE THIS VALUE
options_ami.nz = 0; % MUST NOT CHANGE THIS VALUE
options_ami.nnz = 90; % MUST NOT CHANGE THIS VALUE
sol.x = NaN(length(tout),18);
sol.y = NaN(length(tout),9);
sol.xdot = NaN(1,18);
sol.J = NaN(18,18);
sol.dydx = NaN(9,18);
sol.dydp = NaN(9,options_ami.np);
sol.dxdotdp = NaN(18,options_ami.np);
plist = options_ami.sens_ind-1;
if(nargin>=4)
    if(~isempty(varargin{4}))
        data=varargin{4};
    else
        data.Y=NaN(length(tout),options_ami.ny);
        data.Sigma_Y=-ones(length(tout),options_ami.ny);
    end
else
    data.Y=NaN(length(tout),options_ami.ny);
    data.Sigma_Y= NaN(length(tout),options_ami.ny);
end
if(isfield(data,'T'))
    options_ami.nmaxevent = size(data.T,1);
else
    data.Z = NaN(options_ami.nmaxevent,options_ami.nz);
    data.Sigma_Z = NaN(options_ami.nmaxevent,options_ami.nz);
end
if(options_ami.sensi>0)
    sol.llhS = zeros(length(options_ami.sens_ind),1);
    sol.xS = zeros(length(tout),18,length(options_ami.sens_ind));
    sol.yS = zeros(length(tout),9,length(options_ami.sens_ind));
    sol.zS =  NaN(options_ami.nmaxevent,0,length(options_ami.sens_ind));
end
if(max(options_ami.sens_ind)>116)
    error('Sensitivity index exceeds parameter dimension!')
end
if(isfield(options_ami,'sx0'))
    if(size(options_ami.sx0,2)~=options_ami.np)
        error('Number of rows in sx0 field does not agree with number of model parameters!');
    end
    options_ami.sx0 = bsxfun(@times,options_ami.sx0,1./(permute(theta(options_ami.sens_ind),[2,1])*log(10)));
end
ami_villab2(sol,tout,theta(1:116),kappa(1:0),options_ami,plist,pbar,xscale,data);
if(options_ami.sensi==1)
    sol.sllh = sol.llhS.*theta(options_ami.sens_ind)*log(10);
    sol.sx = bsxfun(@times,sol.xS,permute(theta(options_ami.sens_ind),[3,2,1])*log(10));
    sol.sy = bsxfun(@times,sol.yS,permute(theta(options_ami.sens_ind),[3,2,1])*log(10));
    sol.sz = bsxfun(@times,sol.zS,permute(theta(options_ami.sens_ind),[3,2,1])*log(10));
    sol = rmfield(sol,'llhS');
    sol = rmfield(sol,'xS');
    sol = rmfield(sol,'yS');
    sol = rmfield(sol,'zS');
end
if(options_ami.sensi_meth == 3)
    sol.dxdotdp = bsxfun(@times,sol.dxdotdp,permute(theta(options_ami.sens_ind),[2,1])*log(10));
    sol.dydp = bsxfun(@times,sol.dydp,permute(theta(options_ami.sens_ind),[2,1])*log(10));
    sol.sx = -sol.J\sol.dxdotdp;
    sol.sy = sol.dydx*sol.sx + sol.dydp;
end
if(nargout>1)
    varargout{1} = sol.status;
    varargout{2} = sol.t;
    varargout{3} = sol.x;
    varargout{4} = sol.y;
    if(nargout>4)
        varargout{5} = sol.sx;
        varargout{6} = sol.sy;
    end
else
    varargout{1} = sol;
end
end
