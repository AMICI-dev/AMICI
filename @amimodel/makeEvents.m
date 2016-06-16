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
                if(ievent>0)
                    %check whether we already had this trigger function
                    ediscf = zeros(length(trigger));
                    for id = 1:length(trigger)
                        ediscf(id) = isequaln(abs(trigger{id}),abs(sym(arg)));
                    end
                    if(~any(ediscf))
                        ievent = ievent + 1;
                        trigger{ievent} = sym(arg);
                        bolus{ievent} = sym(zeros(nx,1));
                        z{ievent} = sym.empty([0,0]);
                    else
                    end
                else
                    ievent = ievent + 1;
                    trigger{ievent} = sym(arg);
                    bolus{ievent} = sym(zeros(nx,1));
                    z{ievent} = sym.empty([0,0]);
                end
            end
        end
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
        
        for ix = 1:nx
            symchar = char(this.sym.xdot(ix));
            if(strfind(symchar,'dirac'))
                for ievent = 1:nevent
                    % remove the dirac function and replace by adding
                    % respective bolus
                    triggerchar = char(trigger{ievent});
                    str_arg_d = ['dirac(' triggerchar ')' ];
                    % replace event specific functions by "polydirac" variable
                    symchar = strrep(symchar,str_arg_d,'polydirac');
                    % extract coefficient
                    [cfp,t] = coeffs(sym(symchar),polydirac);
                    if(any(double(t==sym('polydirac'))))
                        tmp_bolus{ievent}(ix) = tmp_bolus{ievent}(ix) + cfp(logical(t==sym('polydirac')));
                    end
                    % remove dirac
                    symchar = strrep(symchar,'polydirac','0');
                end
            end
            if(strfind(symchar,'heaviside'))
                
                for ievent = 1:nevent
                    % check whether we have a dynamic event
                    dtriggerdt = diff(trigger{ievent},'t') + jacobian(trigger{ievent},this.sym.x)*this.sym.xdot(:);
                    if(dtriggerdt~=0)
                        % remove the heaviside function and replace by h
                        % variable which is update on event occurrence in the
                        % solver
                        triggerchar = char(trigger{ievent});
                        str_arg_h = ['heaviside(' triggerchar ')' ];
                        symchar = strrep(symchar,str_arg_h,['h_' num2str(ievent-1)]);
                        mtriggerchar = char(-trigger{ievent});
                        str_arg_h = ['heaviside(' mtriggerchar ')' ];
                        symchar = strrep(symchar,str_arg_h,['(1-h_' num2str(ievent-1) ')']);
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
                            [cfp,tfp] = coeffs(sym(symchar),sym(['h_' num2str(ievent-1) ]));
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
            end
            % update xdot
            this.sym.xdot(ix) = sym(symchar);
        end
        % multiply by the dtriggerdt factor, this should stay here as we
        % want the xdot to be cleaned of any dirac functions
        ievent = 1;
        while ievent <= length(trigger)
            dtriggerdt = diff(trigger{ievent},'t') + jacobian(trigger{ievent},this.sym.x)*this.sym.xdot(:);
            if(dtriggerdt == 0)
                trigger(ievent) = [];
                if(any(bolus{ievent}~=0))
                    warning(['Event ' num2str(ievent) ' has a constant trigger function but non-zero bolus.' ...
                        ' This event will be ommited in the model. Please check your model definition!'])
                end
                bolus(ievent) = [];
                if(~isempty(z{ievent}))
                    warning(['Event ' num2str(ievent) ' has a constant trigger function but non-empty output.' ...
                        ' This event will be ommited in the model. Please check your model definition!'])
                end
                z(ievent) = [];
            else
                bolus{ievent} = bolus{ievent} + tmp_bolus{ievent}/abs(dtriggerdt);
                ievent = ievent+1;
            end
        end 
        nevent = length(trigger);
        
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
            this.sym.Jz = this.sym.Jz + sym(['log(2*pi*sdz_' num2str(iz) '^2) + ((z_' num2str(iz) '-mz_' num2str(iz) ')/sdz_' num2str(iz) ')^2']);
        end
    end
    
end

