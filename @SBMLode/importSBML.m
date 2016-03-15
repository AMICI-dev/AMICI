function importSBML(this,modelname)
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
    this.compartment = subs(this.compartment,compartments_sym([model.compartment.isSetSize]),sym([model.compartment.size]));
    
    % set remaining ones to default value 1
    this.compartment = subs(this.compartment,compartments_sym,sym(ones(size(compartments_sym))));
    
    rulevars = sym({model.rule.variable});
    rulemath = sym({model.rule.formula});
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
    cond_idx = and(logical([model.species.boundaryCondition]),transpose(logical(this.state==this.initState)));
    bound_idx = and(logical([model.species.boundaryCondition]),transpose(logical(this.state~=this.initState)));
    condition_sym = this.state(cond_idx);
    conditions = this.state(cond_idx);
    boundary_sym = this.state(bound_idx);
    boundaries = this.initState(bound_idx);

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
    
    %% RULES
    
    fprintf('applying rules ...\n')
    
    rule_types = {model.rule.typecode};
    if(any(strcmp(rule_types,'SBML_ALGEBRAIC_RULE')))
        %DAE part TBD
        error('DAEs currently not supported!');
    end
    
    if(any(strcmp(rule_types,'SBML_RATE_RULE')))
        %custom reate laws TBD
        error('custom rate laws currently not supported!');
    end
    
    %% REACTIONS
    
    fprintf('parsing reactions ...\n')
    
    nr = length(model.reaction);
    
    this.flux = sym(cellfun(@(x) x.math,{model.reaction.kineticLaw},'UniformOutput',false));
    this.flux = this.flux(:);
    % add local parameters to global parameters, make them global by
    % extending them by the reaction id
    tmp = cellfun(@(x,y) sym(cellfun(@(x) [x '_' y],{x.parameter.id},'UniformOutput',false)),{model.reaction.kineticLaw},{model.reaction.id},'UniformOutput',false);
    plocal = transpose([tmp{:}]);
    tmp = cellfun(@(x) cellfun(@double,{x.parameter.value}),{model.reaction.kineticLaw},'UniformOutput',false);
    pvallocal = transpose([tmp{:}]);
    
    % replace local parameters by globalized ones
    tmp = cellfun(@(x,y,z) subs(x,sym({y.parameter.id}),sym(cellfun(@(x) [x '_' z],{y.parameter.id},'UniformOutput',false))),transpose(num2cell(this.flux)),{model.reaction.kineticLaw},{model.reaction.id},'UniformOutput',false);
    this.flux = [tmp{:}];
    this.flux = this.flux(:);
    
    parameter_sym = [parameter_sym;plocal];
    parameter_val = [parameter_val;pvallocal];
   
    % convert to macroscopic rates
    
    species_idx = transpose(sym(1:nx));
    reactants = cellfun(@(x) {x.species},{model.reaction.reactant},'UniformOutput',false);
    % species index of the reactant
    reactant_sidx = double(subs(cat(2,reactants{:}),species_sym,species_idx));
    % reaction index
    reactant_ridx = cumsum(cell2mat(cellfun(@(x) [ones(1,min(length(x),1)),zeros(1,max(length(x)-1,0))],reactants,'UniformOutput',false)));
    products = cellfun(@(x) {x.species},{model.reaction.product},'UniformOutput',false);
    % species index of the product
    product_sidx = double(subs(cat(2,products{:}),species_sym,species_idx));
    % reaction index
    product_ridx = cumsum(cell2mat(cellfun(@(x) [ones(1,min(length(x),1)),zeros(1,max(length(x)-1,0))],products,'UniformOutput',false)));
    reactant_stochiometry = cellfun(@(x) {x.stoichiometry},{model.reaction.reactant},'UniformOutput',false);
    product_stochiometry = cellfun(@(x) {x.stoichiometry},{model.reaction.product},'UniformOutput',false);
    eS = sym(zeros(nx,nr));
    pS = sym(zeros(nx,nr));
    tmp = cat(2,reactant_stochiometry{:});
    eS(sub2ind(size(eS),reactant_sidx,reactant_ridx)) = [tmp{:}];
    tmp = cat(2,product_stochiometry{:});
    pS(sub2ind(size(eS),product_sidx,product_ridx)) = [tmp{:}];
    
    this.stochiometry = - eS + pS;
    
    this.xdot = this.stochiometry*this.flux;
    
    %% CONVERSION FACTORS/VOLUMES
    
    fprintf('converting to concentrations ...\n')
    
    %extract model conversion factor
    if(isfield(model,'conversionFactor'))
        conversionfactor = model.conversionFactor*ones(nx,1);
    else
        conversionfactor = ones(nx,1);
    end
    
    if(isfield(model.species,'conversionFactor'))
        conversionfactor(~isnan(model.spevies.conversionFactor)) = model.species.conversionFactor(~isnan(model.species.conversionFactor));
    end
    
    this.xdot = conversionfactor.*subs(this.xdot,this.state,this.state.*this.volume)./this.volume;
    
    %% EVENTS
    
    fprintf('loading events ...\n')
    try
        this.trigger = sym({model.event.trigger});
    catch
        this.trigger = sym({model.event.trigger.math});
    end
    this.trigger = this.trigger(:);
    this.trigger = subs(this.trigger,sym('ge'),sym('am_ge'));
    this.trigger = subs(this.trigger,sym('gt'),sym('am_gt'));
    this.trigger = subs(this.trigger,sym('le'),sym('am_le'));
    this.trigger = subs(this.trigger,sym('lt'),sym('am_lt'));
    
    this.bolus = repmat(-this.state,[1,length(this.trigger)]);
    
    tmp = cellfun(@(x) {x.variable},{model.event.eventAssignment},'UniformOutput',false);
    assignments = sym(cat(2,tmp{:}));
    assignments_tidx = cumsum(cell2mat(cellfun(@(x) [ones(1,min(length(x),1)),zeros(1,max(length(x)-1,0))],tmp,'UniformOutput',false)));
    
    tmp = cellfun(@(x) {x.math},{model.event.eventAssignment},'UniformOutput',false);
    assignments_math = sym(cat(2,tmp{:}));
    
    state_assign_idx = ismember(assignments,species_sym(not(any([cond_idx;const_idx;bound_idx]))));
    param_assign_idx = ismember(assignments,parameter_sym);
    cond_assign_idx = ismember(assignments,condition_sym);
    bound_assign_idx = ismember(assignments,boundary_sym);
    
    if(np>0)
    parameter_idx = transpose(sym(1:np));
    assignments_param = assignments(param_assign_idx);
    assignments_param_pidx = double(subs(assignments_param,parameter_sym(1:np),parameter_idx));
    assignments_param_tidx = assignments_tidx(param_assign_idx);
    
    this.param(assignments_param_pidx) = this.param(assignments_param_pidx).*heaviside(this.trigger(assignments_param_tidx));
    end
    
    if(nk>0)
    condition_idx = transpose(sym(1:nk));
    assignments_cond = assignments(cond_assign_idx);
    assignments_cond_kidx = double(subs(assignments_cond,condition_sym,condition_idx));
    assignments_cond_tidx = assignments_tidx(cond_assign_idx);
    
    conditions(assignments_cond_kidx) = conditions(assignments_cond_kidx).*heaviside(this.trigger(assignments_cond_tidx));
    end
    
    if(length(boundaries)>0)
    boundary_idx = transpose(sym(1:length(boundaries)));
    assignments_bound = assignments(bound_assign_idx);
    assignments_bound_bidx = double(subs(assignments_bound,boundary_sym,boundary_idx));
    assignments_bound_tidx = assignments_tidx(bound_assign_idx);
    
    boundaries(assignments_bound_bidx) = conditions(assignments_bound_bidx).*heaviside(this.trigger(assignments_bound_tidx));
    end
    
    assignments_state = assignments(state_assign_idx);
    assignments_state_sidx = double(subs(assignments_state,species_sym,species_idx));
    assignments_state_tidx = assignments_tidx(state_assign_idx);
    
    addToBolus = sym(zeros(size(this.bolus)));
    addToBolus(sub2ind(size(addToBolus),assignments_state_sidx,assignments_state_tidx)) = assignments_math(state_assign_idx);
    
    nobolus_idx = not(any(addToBolus~=0));
    this.bolus(:,nobolus_idx) = [];
    this.trigger(nobolus_idx) = [];
    addToBolus(:,nobolus_idx) = [];
    this.bolus = this.bolus + addToBolus;
    
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
    end
    
    %% CLEAN-UP
    
    fprintf('cleaning up ...\n')
    
    % remove constant/condition states
    this.state(any([cond_idx;const_idx;bound_idx])) = [];
    this.volume(any([cond_idx;const_idx;bound_idx])) = [];
    this.initState(any([cond_idx;const_idx;bound_idx])) = [];
    this.xdot(any([cond_idx;const_idx;bound_idx])) = [];
    this.bolus(any([cond_idx;const_idx;bound_idx]),:) = [];
    
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
