function makeEvents( this )
% makeEvents extracts discontiniuties from the model right hand side
% and converts them into events
%
% Parameters:
%
% Return values:
%  void

nevent = length(this.event);
nx = length(this.sym.x);

% extract trigger and bolus, also check dimensions here
for ievent = 1:nevent
    trigger{ievent} = this.event(ievent).trigger;
    if(length(this.event(ievent).bolus) == nx)
        bolus{ievent} = this.event(ievent).bolus(:);
    elseif(numel(this.event(ievent).bolus) == 1)
        bolus{ievent} = this.event(ievent).bolus*ones(nx,1);
    else
        error(['Bolus of event ' num2str(ievent) ' is neither a scalar nor does it match the state dimension!']);
    end
    z{ievent} = transpose(this.event(ievent).z(:));
end

ievent = nevent;
nx = length(this.sym.x);
triggerstr = '';
triggers = {};

% collect all the triggers
for ix = 1:nx
    tmp_str = char(this.sym.xdot(ix));
    % find discontinuity inducing functions
    idx_start = [strfind(tmp_str,'heaviside') + 9,strfind(tmp_str,'dirac') + 5]; % find
    if(~isempty(idx_start))
        % compute bracket level
        brl = cumsum(tmp_str == '(') - cumsum(tmp_str == ')');
        for iocc = 1:length(idx_start) % loop over occurances
            % extract argument
            idx_end = find(brl(idx_start(iocc):end)-brl(idx_start(iocc))==-1,1,'first');
            arg = tmp_str((idx_start(iocc)+1):(idx_start(iocc)+idx_end-2));
            triggers{end+1} = arg;
            triggers{end+1} = ['-(' arg ')']; % for dirac we need both +to- and -to+ transitions
        end
    end
end

% select the unique ones
utriggers = unique(arrayfun(@char,betterSym(triggers),'UniformOutput',false));
for itrigger = 1:length(utriggers)
    ievent = ievent + 1;
    trigger{ievent} = betterSym(utriggers{itrigger});
    bolus{ievent} = sym(zeros(nx,1));
    z{ievent} = sym.empty([0,0]);
end

% update number of events
if(exist('trigger','var'))
    nevent = length(trigger);
else
    nevent = 0;
end

if(nevent>0)
    % initialis tmp_bolus
    tmp_bolus = cell([nevent,1]);
    for ievent = 1:nevent
        tmp_bolus{ievent} = sym(zeros([nx,1]));
    end
    syms polydirac

    % initialise hflag
    hflags = zeros([nx,nevent]);

    % heaviside
    event_dependency = zeros(nevent);
    for ievent = 1:nevent
        symchar = char(trigger{ievent});
        if(~isempty(strfind(symchar,'heaviside')))
            for jevent = 1:nevent
                % remove the heaviside function and replace by h
                % variable which is update on event occurrence in the
                % solver
                triggerchar = char(trigger{jevent});
                str_arg_h = ['heaviside(' triggerchar ')' ];
                event_dependency(ievent,jevent) = event_dependency(ievent,jevent) + ~isempty(strfind(symchar,str_arg_h));
                mtriggerchar = char(-trigger{jevent});
                str_arg_h = ['heaviside(' mtriggerchar ')' ];
                event_dependency(ievent,jevent) = event_dependency(ievent,jevent) + ~isempty(strfind(symchar,str_arg_h));
            end
        end
    end

    % check for loops
    if(any(any(event_dependency^(size(event_dependency,1)))))
        error('Found loop in trigger dependency. This can lead to the simulation getting stuck and is thus currently not supported. Please check your model definition!')
    end

    P = 1:size(event_dependency,1);

    % make matrix upper triangular, this is to ensure that we dont end
    % up with partially replaced trigger functions that we no longer recognise
    while(~isempty(find(triu(event_dependency(P,P))-event_dependency(P,P))))
        for ip = 1:length(P)
            jp = ip+find(event_dependency(P(ip+1:end),P(ip)),1,'first');
            while ~isempty(jp)
                pp = P(ip);
                P(ip) = P(jp);
                P(jp) = pp;
                jp = ip+find(event_dependency(P(ip+1:end),P(ip)),1,'first');
            end
        end
    end
    trigger = trigger(P);
    bolus = bolus(P);
    z = z(P);



    for ix = 1:nx
        symchar = char(this.sym.xdot(ix));
        symvariable = this.sym.xdot(ix);
        if(strfind(symchar,'dirac'))
            for ievent = 1:nevent
                % remove the dirac function and replace by adding
                % respective bolus
                symvariable = subs(symvariable,dirac(trigger{ievent}),sym('polydirac'));
                % extract coefficient
                [cfp,t] = coeffs(symvariable,polydirac);
                if(any(double(t==sym('polydirac'))))
                    tmp_bolus{ievent}(ix) = tmp_bolus{ievent}(ix) + cfp(logical(t==sym('polydirac')));
                    idx_mirror = find(cellfun(@(x) double(x==-trigger{ievent}),trigger));
                    if(~isempty(idx_mirror))
                        tmp_bolus{idx_mirror}(ix) = tmp_bolus{idx_mirror}(ix) + cfp(logical(t==sym('polydirac'))); % for dirac we need both +to- and -to+ transitions
                    end
                end
                % remove dirac
                symvariable = subs(symvariable,sym('polydirac'),sym('0'));
            end
        end
        if(strfind(symchar,'heaviside'))

            for ievent = 1:nevent
                % remove the heaviside function and replace by h
                % variable which is updated upon event occurrence in the
                % solver

                % h variables only change for one sign change but heaviside
                % needs updating for both, thus we should
                symvariable = subs(symvariable,heaviside( trigger{ievent}),betterSym(['h_' num2str(ievent-1)']));
                symvariable = subs(symvariable,heaviside(-trigger{ievent}),betterSym(['(1-h_' num2str(ievent-1) ')']));
                % set hflag

                % we can check whether dividing cfp(2) by
                % trigger{ievent} reduced the length of the symbolic
                % expression. If it does, this suggests that
                % trigger{ievent} is a factor of cfp(2), which will be
                if(or(...
                    ismember(sym(['h_' num2str(ievent-1)']),symvar(symvariable)),...
                    ismember(sym(['h_' num2str(find(-trigger{ievent}==trigger)-1)']),symvar(symvariable))...
                    ))
                    hflags(ix,ievent) = 1;
                else
                    hflags(ix,ievent) = 0;
                end
            end
        end
        if(strfind(char(symvariable),'heaviside'))
            warning(['Missed heaviside function in state ' num2str(ix) '. AMICI will continue and produce a running model, but this should be fixed!'])
        end
        % update xdot
        this.sym.xdot(ix) = symvariable;
    end

    % loop until we no longer found any dynamic heaviside functions in the triggers in the previous loop
    nheavy = 1;
    while nheavy>0
        nheavy = 0;
        for ievent = 1:nevent
            symchar = char(trigger{ievent});
            for jevent = 1:nevent
                % remove the heaviside function and replace by h
                % variable which is update on event occurrence in the
                % solver
                triggerchar = char(trigger{jevent});
                str_arg_h = ['heaviside(' triggerchar ')' ];
                nheavy = nheavy + length(strfind(symchar,str_arg_h));
                symchar = strrep(symchar,str_arg_h,['h_' num2str(jevent-1)]);
                mtriggerchar = char(-trigger{jevent});
                str_arg_h = ['heaviside(' mtriggerchar ')' ];
                nheavy = nheavy + length(strfind(symchar,str_arg_h));
                symchar = strrep(symchar,str_arg_h,['(1-h_' num2str(jevent-1) ')']);
                % set hflag
            end
            trigger{ievent} = betterSym(symchar);
        end
    end

    % compute dtriggerdt and constant trigger functions
    for ievent = 1:nevent
        dtriggerdt(ievent) = diff(trigger{ievent},sym('t')) + jacobian(trigger{ievent},this.sym.x)*this.sym.xdot(:);
    end
    triggeridx = logical(dtriggerdt~=0);

    % multiply by the dtriggerdt factor, this should stay here as we
    % want the xdot to be cleaned of any dirac functions
    ievent = 1;
    for ievent = 1:nevent
        if(not(triggeridx(ievent)))
            if(any(bolus{ievent}~=0))
                error(['Event ' num2str(ievent) ' has a constant trigger function but non-zero bolus.' ...
                    ' Please check your model definition!'])
            end
            if(~isempty(z{ievent}))
                error(['Event ' num2str(ievent) ' has a constant trigger function but non-empty output.' ...
                    ' Please check your model definition!'])
            end
        else
            bolus{ievent} = bolus{ievent} + tmp_bolus{ievent}/abs(dtriggerdt(ievent));
            ievent = ievent+1;
        end
    end

    % update hflags according to bolus
    for ievent = 1:nevent
        if(any(double(bolus{ievent}~=0)))
            for ix = 1:nx
                if(~hflags(ix,ievent))
                    for j = transpose(find(double(bolus{ievent}~=0)))
                        if(ismember(this.sym.x(j),symvar(this.sym.xdot(ix))))
                            hflags(ix,ievent) = 1;
                        end
                    end
                end
            end
        end
    end

    this.event = amievent.empty();

    % update events
    for ievent = 1:nevent
        this.event(ievent) = amievent(trigger{ievent},bolus{ievent}(:),z{ievent});
        % do not add a (:) after z{ievent} this will transform an
        % [ empty sym ] into Empty sym: 0-by-1 which will lead to a
        % zero entry if we apply [this.event.z]
        this.event(ievent) = this.event(ievent).setHflag(hflags(:,ievent));
    end
end

if(~isfield(this.sym,'sigma_z'))
    this.sym.sigma_z = sym(ones(length([this.event.z]),1));
end
if(numel(this.sym.sigma_z) == 1)
    this.sym.sigma_z = this.sym.sigma_z*sym(ones(length([this.event.z]),1));
end

if(~isfield(this.sym,'Jz'))
    this.sym.Jz = sym(zeros(length([this.event.z]),1));
    for iz = 1:length([this.event.z])
        this.sym.Jz(iz) = betterSym(['0.5*log(2*pi*sigma_z_' num2str(iz-1) '^2) + 0.5*((z_' num2str(iz-1) '-mz_' num2str(iz-1) ')/sigma_z_' num2str(iz-1) ')^2']);
    end
end

z = sym(arrayfun(@(iz) ['z_' num2str(iz-1)],1:length([this.event.z]),'UniformOutput',false));
var_z = sym(arrayfun(@(iz) ['var_z_' num2str(iz-1)],1:length([this.event.z]),'UniformOutput',false));
sigma_z = sym(arrayfun(@(iz) ['sigma_z_' num2str(iz-1)],1:length([this.event.z]),'UniformOutput',false));
var_sigma_z = sym(arrayfun(@(iz) ['sigma_z_' num2str(iz-1)],1:length([this.event.z]),'UniformOutput',false));

this.sym.Jz = subs(this.sym.Jz,z,var_z);
this.sym.Jz = subs(this.sym.Jz,sigma_z,var_sigma_z);

rz = sym(arrayfun(@(iz) ['rz_' num2str(iz-1)],1:length([this.event.z]),'UniformOutput',false));
mz = sym(arrayfun(@(iz) ['mz_' num2str(iz-1)],1:length([this.event.z]),'UniformOutput',false));
var_rz = sym(arrayfun(@(iz) ['var_rz_' num2str(iz-1)],1:length([this.event.z]),'UniformOutput',false));

if(~isfield(this.sym,'Jrz'))
    this.sym.Jrz = sym(zeros(size(this.sym.Jz)));
    for iz = 1:length([this.event.z])
        tmp = subs(this.sym.Jz(iz,:),var_z,var_rz);
        this.sym.Jrz(iz,:) = subs(tmp,mz,sym(zeros(size(mz))));
    end
end

this.sym.Jrz = subs(this.sym.Jrz,rz,var_rz);

end
