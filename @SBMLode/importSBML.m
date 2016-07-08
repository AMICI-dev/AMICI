function importSBML(this,modelname)
% importSBML parses information from the SBML definition and populates
% the SBMLode object from this information.
%
% Parameters:
%  modelname: target name of the model
%
% Return values:
%

try
    model = TranslateSBML([modelname '.sbml']);
catch
    model = TranslateSBML([modelname '.xml']);
end

%% COMPARTMENTS

fprintf('loading compartments ...\n')

% initialise
compartments_sym = sym({model.compartment.id});
this.compartment = compartments_sym;

% set initial assignments
if(isfield(model,'initialAssignment'))
    initAssignemnts_sym = sym({model.initialAssignment.symbol});
    initAssignemnts_math = sym({model.initialAssignment.math});
else
    initAssignemnts_sym = sym([]);
    initAssignemnts_math = sym([]);
end
setInitialAssignment(this,model,'compartment',initAssignemnts_sym,initAssignemnts_math)

% extract compartment sizes, default to 1
this.compartment = subs(this.compartment,compartments_sym(logical([model.compartment.isSetSize])),sym([model.compartment.size]));

% set remaining ones to default value 1
this.compartment = subs(this.compartment,compartments_sym,sym(ones(size(compartments_sym))));

%% RULES

fprintf('applying rules ...\n')

rule_types = {model.rule.typecode};
if(any(strcmp(rule_types,'SBML_ALGEBRAIC_RULE')))
    %DAE part TBD
    error('DAEs currently not supported!');
end

rulevars = sym({model.rule.variable});
rulemath = sym({model.rule.formula});
% remove rate rules
rulevars = rulevars(not(strcmp({model.rule.typecode},'SBML_RATE_RULE')));
rulemath = rulemath(not(strcmp({model.rule.typecode},'SBML_RATE_RULE')));
repeat_idx = ismember(rulevars,symvar(rulemath));
while(any(repeat_idx))
    rulemath= subs(rulemath,rulevars,rulemath);
    repeat_idx = ismember(rulevars,symvar(rulemath));
end

applyRule(this,model,'compartment',rulevars,rulemath)

%% SPECIES

fprintf('loading species ...\n')

% extract species
species_sym = sym({model.species.id});
species_sym = species_sym(:);
this.state = species_sym;

nx = length(this.state);

% extract corresponding volumes
compartments = sym({model.species.compartment});
this.volume = subs(compartments(:),sym({model.compartment.id}),this.compartment);

initConcentration = [model.species.initialConcentration];
initConcentration = initConcentration(:);
initAmount = [model.species.initialAmount];
initAmount = initAmount(:);
this.initState = species_sym;

% set initial assignments
setInitialAssignment(this,model,'initState',initAssignemnts_sym,initAssignemnts_math)


% remove conditions species (boundary condition + no initialisation)
if(~isempty(this.state))
    cond_idx = and(logical([model.species.boundaryCondition]),transpose(logical(this.state==this.initState)));
    bound_idx = and(logical([model.species.boundaryCondition]),transpose(logical(this.state~=this.initState)));
    condition_sym = this.state(cond_idx);
    conditions = this.state(cond_idx);
    boundary_sym = this.state(bound_idx);
    boundaries = this.initState(bound_idx);
else
    conditions = this.state;
    cond_idx = [];
    bound_idx = [];
    condition_sym = sym([]);
    boundary_sym = sym([]);
    boundaries = this.state;
end

nk = length(conditions);

applyRule(this,model,'condition',rulevars,rulemath)

% set initial concentrations
concentration_idx = logical([model.species.isSetInitialConcentration]);
this.initState = subs(this.initState,species_sym(concentration_idx),sym(initConcentration(concentration_idx)));

% set initial amounts
amount_idx = logical([model.species.isSetInitialAmount]);
this.initState = subs(this.initState,species_sym(amount_idx),sym(initAmount(amount_idx))./this.volume(amount_idx));

% apply rules
applyRule(this,model,'initState',rulevars,rulemath)

this.knom = double(this.initState(cond_idx));

% remove constant species
const_idx = logical([model.species.constant]) & not(cond_idx);
constant_sym = this.state(const_idx);
constants = this.initState(const_idx);

%% PARAMETERS

fprintf('loading parameters ...\n')

% extract this.param
parameter_sym = sym({model.parameter.id});
parameter_val = transpose([model.parameter.value]);
parameter_sym = parameter_sym(:);
this.param = parameter_sym;

% set initial assignments
setInitialAssignment(this,model,'parameter',initAssignemnts_sym,initAssignemnts_math)
applyRule(this,model,'param',rulevars,rulemath)

np = length(this.param);


%% REACTIONS

fprintf('parsing reactions ...\n')

nr = length(model.reaction);

kLaw = cellfun(@(x) x.math,{model.reaction.kineticLaw},'UniformOutput',false);

if(any(cell2mat(strfind(kLaw,'factorial'))))
    error('Sorry, factorial functions are not supported at the moment!')
end

if(any(cell2mat(strfind(kLaw,'ceil'))))
    error('Sorry, ceil functions are not supported at the moment!')
end

this.flux = sym(kLaw);
this.flux = this.flux(:);
% add local parameters to global parameters, make them global by
% extending them by the reaction id
species_idx = transpose(sym(1:nx));
if(length({model.reaction.id})>0)
    try
        tmp = cellfun(@(x,y) sym(cellfun(@(x) [x '_' y],{x.parameter.id},'UniformOutput',false)),{model.reaction.kineticLaw},{model.reaction.id},'UniformOutput',false);
        plocal = transpose([tmp{:}]);
        tmp = cellfun(@(x) cellfun(@double,{x.parameter.value}),{model.reaction.kineticLaw},'UniformOutput',false);
        pvallocal = transpose([tmp{:}]);
        % replace local parameters by globalized ones
        tmp = cellfun(@(x,y,z) subs(x,sym({y.parameter.id}),sym(cellfun(@(x) [x '_' z],{y.parameter.id},'UniformOutput',false))),transpose(num2cell(this.flux)),{model.reaction.kineticLaw},{model.reaction.id},'UniformOutput',false);
        this.flux = [tmp{:}];
        this.flux = this.flux(:);
        
    catch
        tmp = cellfun(@(x,y) sym(cellfun(@(x) [x '_' y],{x.localParameter.id},'UniformOutput',false)),{model.reaction.kineticLaw},{model.reaction.id},'UniformOutput',false);
        plocal = transpose([tmp{:}]);
        tmp = cellfun(@(x) cellfun(@double,{x.localParameter.value}),{model.reaction.kineticLaw},'UniformOutput',false);
        pvallocal = transpose([tmp{:}]);
        % replace local parameters by globalized ones
        tmp = cellfun(@(x,y,z) subs(x,sym({y.localParameter.id}),sym(cellfun(@(x) [x '_' z],{y.localParameter.id},'UniformOutput',false))),transpose(num2cell(this.flux)),{model.reaction.kineticLaw},{model.reaction.id},'UniformOutput',false);
        this.flux = [tmp{:}];
        this.flux = this.flux(:);
        
    end
    
    parameter_sym = [parameter_sym;plocal];
    parameter_val = [parameter_val;pvallocal];
    
    reactants = cellfun(@(x) {x.species},{model.reaction.reactant},'UniformOutput',false);
    % species index of the reactant
    reactant_sidx = double(subs(cat(2,reactants{:}),species_sym,species_idx));
    % reaction index
    tmp = cumsum(cell2mat(cellfun(@(x) [ones(1,1),zeros(1,max(length(x)-1,0))],reactants,'UniformOutput',false)));
    wreact = cell2mat(cellfun(@(x) [ones(1,length(x)),zeros(1,isempty(x))],reactants,'UniformOutput',false));
    reactant_ridx = tmp(logical(wreact));
    products = cellfun(@(x) {x.species},{model.reaction.product},'UniformOutput',false);
    % species index of the product
    product_sidx = double(subs(cat(2,products{:}),species_sym,species_idx));
    % reaction index
    tmp = cumsum(cell2mat(cellfun(@(x) [ones(1,1),zeros(1,max(length(x)-1,0))],products,'UniformOutput',false)));
    wprod = cell2mat(cellfun(@(x) [ones(1,length(x)),zeros(1,isempty(x))],products,'UniformOutput',false));
    product_ridx = tmp(logical(wprod));
    if(model.SBML_level>=3)
        reactant_stochiometry = cellfun(@(x) {x.stoichiometry},{model.reaction.reactant},'UniformOutput',false);
        product_stochiometry = cellfun(@(x) {x.stoichiometry},{model.reaction.product},'UniformOutput',false);
    else
        % addition is necessary due to 1x0 struct that is returned by libSBML which is not properly handled by MATLAB,
        % the concatenation is necessary because MATLAB does reallly weird things when you apply arrayfuns to fields of
        % 1x0 structs. Let's pray this works and these bugs are fixed in the future
        math_expr = @(y) sym(y.math);
        symbolic_expr = @(x) num2cell([sym(arrayfun(@(y) math_expr(y),[x.stoichiometryMath],'UniformOutput',false)),sym(zeros(1,isempty([x.stoichiometryMath])*length(x)))] + sym(arrayfun(@(y) y.stoichiometry,x))*isempty([x.stoichiometryMath]));
        reactant_stochiometry = cellfun(@(x) symbolic_expr(x),{model.reaction.reactant},'UniformOutput',false);
        product_stochiometry = cellfun(@(x) symbolic_expr(x),{model.reaction.product},'UniformOutput',false);
    end
    eS = sym(zeros(nx,nr));
    pS = sym(zeros(nx,nr));
    tmp = cat(2,reactant_stochiometry{:});
    eS(sub2ind(size(eS),reactant_sidx,reactant_ridx)) = [tmp{:}];
    tmp = cat(2,product_stochiometry{:});
    pS(sub2ind(size(eS),product_sidx,product_ridx)) = [tmp{:}];
    
    this.stochiometry = - eS + pS;
    
    this.xdot = this.stochiometry*this.flux;
else
    this.xdot = sym(zeros(size(this.state)));
end

%% RATE RULES

for irule = 1:length(model.rule)
    if(strcmp(model.rule(irule).typecode,'SBML_RATE_RULE'))
        this.xdot(find(this.state == sym(model.rule(irule).variable))) = sym(model.rule(irule).formula);
    end
end
%% CONVERSION FACTORS/VOLUMES

fprintf('converting to concentrations ...\n')

%extract model conversion factor
if(isfield(model,'conversionFactor'))
    if(strcmp(model.conversionFactor,''))
        conversionfactor = ones(nx,1);
    else
        conversionfactor = model.conversionFactor*ones(nx,1);
    end
else
    conversionfactor = ones(nx,1);
end

if(isfield(model.species,'conversionFactor'))
    if(any(not(strcmp({model.species.conversionFactor},''))))
        conversionfactor(~isnan({model.species.conversionFactor})) = model.species.conversionFactor(~isnan(model.species.conversionFactor));
    end
end

this.xdot = conversionfactor.*subs(this.xdot,this.state,this.state.*this.volume)./this.volume;

%% EVENTS

fprintf('loading events ...\n')


try
    tmp = cellfun(@(x) sym(x),{model.event.trigger},'UniformOutput',false);
    this.trigger = [tmp{:}];
catch
    tmp = cellfun(@(x) sym(x.math),{model.event.trigger},'UniformOutput',false);
    this.trigger = [tmp{:}];
end
this.trigger = this.trigger(:);
this.trigger = subs(this.trigger,sym('ge'),sym('am_ge'));
this.trigger = subs(this.trigger,sym('gt'),sym('am_gt'));
this.trigger = subs(this.trigger,sym('le'),sym('am_le'));
this.trigger = subs(this.trigger,sym('lt'),sym('am_lt'));

this.bolus = sym(zeros([length(this.state),length(this.trigger)]));
if(length(this.trigger)>0)
    for ievent = 1:length(this.trigger)
        tmp = cellfun(@(x) {x.variable},{model.event(ievent).eventAssignment},'UniformOutput',false);
        assignments = sym(cat(2,tmp{:}));
        
        tmp = cellfun(@(x) {x.math},{model.event(ievent).eventAssignment},'UniformOutput',false);
        assignments_math = sym(cat(2,tmp{:}));
        
        for iassign = 1:length(assignments)
            state_assign_idx = find(assignments(iassign)==species_sym);
            param_assign_idx = find(assignments(iassign)==parameter_sym);
            cond_assign_idx = find(assignments(iassign)==condition_sym);
            bound_assign_idx = find(assignments(iassign)==boundary_sym);
            
            if(np>0 && ~isempty(param_assign_idx))
                this.param(assignments_param_idx) = this.param(assignments_param_idx)*heaviside(this.trigger(ievent)) + assignments_math(iassign)*heaviside(-this.trigger(ievent));
            end
            
            if(nk>0 && ~isempty(assignments_cond_tidx))
                conditions(assignments_cond_idx) = conditions(assignments_cond_idx)*heaviside(this.trigger(ievent)) + assignments_math(iassign)*heaviside(-this.trigger(ievent));
            end
            
            if(length(boundaries)>0 && ~isempty(bound_assign_idx))
                boundaries(assignments_bound_idx) = conditions(assignments_bound_idx)*heaviside(this.trigger(ievent)) + assignments_math(iassign)*heaviside(-this.trigger(ievent));
            end
            
            if(length(this.state)>0 && ~isempty(state_assign_idx))
                
                this.bolus(state_assign_idx,ievent) = -this.state(state_assign_idx);
                addToBolus = sym(zeros(size(this.bolus(:,ievent))));
                addToBolus(state_assign_idx) = assignments_math(iassign);
                
                this.bolus(:,ievent) = this.bolus(:,ievent) + addToBolus;
            end
        end
    end
else
    addToBolus = sym([]);
end

%% FUNCTIONS

fprintf('loading functions ...\n')

tmp = cellfun(@(x) x(8:end-1),{model.functionDefinition.math},'UniformOutput',false);
lambdas = cellfun(@(x)argScan(x),tmp);

if(~isempty(lambdas))
    this.funmath = cellfun(@(x) x{end},lambdas,'UniformOutput',false);
    % temp replacement for any user defined function
    tmpfun = cellfun(@(x) ['fun_' num2str(x)],num2cell(1:length(model.functionDefinition)),'UniformOutput',false);
    this.funmath = strrep(this.funmath,{model.functionDefinition.id},tmpfun);
    % replace helper functions
    this.funmath = strrep(this.funmath,'ge(','am_ge(');
    this.funmath = strrep(this.funmath,'gt(','am_gt(');
    this.funmath = strrep(this.funmath,'le(','am_le(');
    this.funmath = strrep(this.funmath,'lt(','am_lt(');
    this.funmath = strrep(this.funmath,'piecewise(','am_piecewise(');
    this.funmath = strrep(this.funmath,'max(','am_max(');
    this.funmath = strrep(this.funmath,'min(','am_min(');
    
    this.funmath = strrep(this.funmath,tmpfun,{model.functionDefinition.id});
    this.funarg = cellfun(@(x,y) [y '(' strjoin(transpose(x(1:end-1)),',') ')'],lambdas,{model.functionDefinition.id},'UniformOutput',false);
    
    % make functions available in this file
    
    for ifun = 1:length(this.funmath)
        token = regexp(this.funarg(ifun),'\(([0-9\w\,]*)\)','tokens');
        start = regexp(this.funarg(ifun),'\(([0-9\w\,]*)\)');
        eval([this.funarg{ifun}(1:(start{1}-1)) ' = @(' token{1}{1}{1} ')' this.funmath{ifun} ';']);
    end
end

%% CLEAN-UP

fprintf('cleaning up ...\n')

% remove constant/condition states
this.state(any([cond_idx;const_idx;bound_idx])) = [];
this.volume(any([cond_idx;const_idx;bound_idx])) = [];
this.initState(any([cond_idx;const_idx;bound_idx])) = [];
this.xdot(any([cond_idx;const_idx;bound_idx])) = [];
this.bolus(any([cond_idx;const_idx;bound_idx]),:) = [];

% apply rules to dynamics
for irule = 1:length(rulevars)
    eval(['syms ' strjoin(arrayfun(@char,rulevars(irule),'UniformOutput',false),' ') ' ' strjoin(arrayfun(@char,symvar(sym(rulemath(irule))),'UniformOutput',false),' ')]);
    eval(['rule = ' char(rulemath(irule)) ';']);
    rule_idx = find(rulevars(irule)==this.state);
    this.xdot(rule_idx) = jacobian(rule,this.state)*this.xdot;
    this.initState(rule_idx) = rulemath(irule);
    if(~isempty(this.bolus))
        this.bolus(rule_idx,:) = jacobian(rule,this.state)*this.bolus;
    end
end

state_vars = [symvar(this.xdot),symvar(this.initState)];
event_vars = [symvar(addToBolus),symvar(this.trigger)];

% substitute with actual expressions
makeSubs(this,parameter_sym(1:np),this.param);
makeSubs(this,condition_sym,conditions);
makeSubs(this,constant_sym,constants);
makeSubs(this,boundary_sym,boundaries);
makeSubs(this,compartments_sym,this.compartment);

isUsedParam = or(ismember(parameter_sym,event_vars),ismember(parameter_sym,state_vars));
isPartOfRule = and(ismember(parameter_sym,symvar(sym({model.rule.formula}))),ismember(parameter_sym,symvar(sym({model.rule.variable}))));
isRuleVar = ismember(parameter_sym,sym({model.rule.variable}));
hasAssignment = ismember(parameter_sym,initAssignemnts_sym);
this.parameter = parameter_sym(and(not(isRuleVar),not(isPartOfRule)));
this.pnom = parameter_val(and(not(isRuleVar),not(isPartOfRule)));

this.condition = condition_sym;
obs_idx = all([isRuleVar,not(isPartOfRule),not(isUsedParam),not(hasAssignment)],2);
this.observable = this.param(obs_idx);
this.observable_name = parameter_sym(obs_idx);

equal_idx = logical(this.observable==this.observable_name); % this is the unused stuff
this.observable(equal_idx) = [];
this.observable_name(equal_idx) = [];

this.observable = subs(this.observable,parameter_sym(1:np),this.param);
this.observable = subs(this.observable,condition_sym,conditions);
this.observable = subs(this.observable,constant_sym,constants);
this.observable = subs(this.observable,boundary_sym,boundaries);
this.observable = subs(this.observable,compartments_sym,this.compartment);
this.time_symbol = model.time_symbol;

prohibited = sym({'null','beta'});
alt_prohib = sym({'null_sym','beta_sym'});
this.parameter = subs(this.parameter,prohibited,alt_prohib);
this.state = subs(this.state,prohibited,alt_prohib);
this.condition = subs(this.condition,prohibited,alt_prohib);
end

function setInitialAssignment(this,~,field,initAssignemnts_sym,initAssignemnts_math)
this.(field) = subs(this.(field),initAssignemnts_sym,initAssignemnts_math);
end


function applyRule(this,~,field,rulevars,rulemath)
this.(field) = subs(this.(field),rulevars,rulemath);
end

function args = argScan(str)
brl = computeBracketLevel(str);

l1_idx = find(brl==0); % store indexes

% find commas but only in lowest bracket level
c_idx = strfind(str(brl==0),',');
str(l1_idx(c_idx)) = '@'; % replace by something that should never occur in equations
args = textscan(str,'%s','Whitespace','@');
end

function makeSubs(this,old,new)
this.xdot = subs(this.xdot,old,new);
this.trigger = subs(this.trigger,old,new);
this.bolus = subs(this.bolus,old,new);
this.initState = subs(this.initState,old,new);
end
