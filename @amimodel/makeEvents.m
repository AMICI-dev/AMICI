function makeEvents( this )
% makeEvents extracts discontiniuties from the model right hand side
% and converts them into events

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
            triggers{end+1} = sym(arg);
        end
    end
end

% select the unique ones
abstriggers = cellfun(@(x) char(abs(x)),triggers,'UniformOutput',false);
utriggers = unique(abstriggers);
for itrigger = 1:length(utriggers)
    ievent = ievent + 1;
    % find the original one
    % transform to char once to get the right ordering
    trigger{ievent} = sym(char(triggers{find(strcmp(utriggers{itrigger},abstriggers),1)}));
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
    
    for ievent = 1:nevent
        dtriggerdt(ievent) = diff(trigger{ievent},sym('t')) + jacobian(trigger{ievent},this.sym.x)*this.sym.xdot(:);
    end
    triggeridx = logical(dtriggerdt~=0);
    
    
    event_dependency = zeros(nevent);
    for ievent = 1:nevent
        symchar = char(trigger{ievent});
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
                symvariable = subs(symvariable,heaviside( trigger{ievent}),sym(['h_'        num2str(ievent-1)]));
                symvariable = subs(symvariable,heaviside(-trigger{ievent}),sym(['(1-h_' num2str(ievent-1) ')']));
                % set hflag
                
                % we can check whether dividing cfp(2) by
                % trigger{ievent} reduced the length of the symbolic
                % expression. If it does, this suggests that
                % trigger{ievent} is a factor of cfp(2), which will be
                % the case for min/max functions. in that case we do
                % not need a hflag as there is no discontinuity in the
                % right hand side. This is not a perfect fix, in the
                % long run one should maybe go back to the old syntax for
                % am_max and am_min?
                try
                    [cfp,tfp] = coeffs(symvariable,sym(['h_' num2str(ievent-1) ]));
                    if(any(double(tfp==sym(['h_' num2str(ievent-1)]))))
                        if(length(char(cfp(logical(tfp==sym(['h_' num2str(ievent-1)])))/trigger{ievent}))<length(char(cfp(logical(tfp==sym(['h_' num2str(ievent-1)]))))))
                            hflags(ix,ievent) = 0;
                        else
                            hflags(ix,ievent) = 1;
                        end
                    else
                        hflags(ix,ievent) = 0;
                    end
                catch
                    hflags(ix,ievent) = 1;
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
            trigger{ievent} = sym(symchar);
        end
    end
    
    % update to prevent dirac functions.
    for ievent = 1:nevent
        dtriggerdt(ievent) = diff(trigger{ievent},sym('t')) + jacobian(trigger{ievent},this.sym.x)*this.sym.xdot(:);
    end
    
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
    this.sym.Jz = sym(0);
    for iz = 1:length([this.event.z])
        this.sym.Jz = this.sym.Jz + sym(['log(2*pi*sigma_z_' num2str(iz) '^2) + ((z_' num2str(iz) '-mz_' num2str(iz) ')/sigma_z_' num2str(iz) ')^2']);
    end
end

end

