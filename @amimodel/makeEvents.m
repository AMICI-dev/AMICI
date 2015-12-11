function [ this ] = makeEvents( this )
    % getEvents extracts discontiniuties from the model right hand side
    % and converts them into events
    %
    % Return values:
    %  this: updated model definition object @type amimodel
    
    nevent = length(this.event);
    
    % extract trigger and bolus
    for ievent = 1:nevent
        trigger{ievent} = this.event(ievent).trigger;
        bolus{ievent} = this.event(ievent).bolus;
        z{ievent} = this.event(ievent).z;
    end
    
    ievent = nevent;
    nx = length(this.sym.xdot);
    
    for ix = 1:nx
        tmp_str = char(this.sym.xdot(ix));
        % find discontinuity inducing functions
        idx_start = [strfind(tmp_str,'heaviside') + 9,strfind(tmp_str,'dirac') + 5]; % find
        if(~isempty(idx_start))
            % compute bracket level
            brl = cumsum(tmp_str == '(') - cumsum(tmp_str == ')');
            for iocc = 1:length(idx_start) % loop over occurances
                % extract argument
                idx_end = find(brl(idx_start(iocc):end)-brl(idx_start(iocc)),'first',1);
                arg = tmp_str((idx_start(iocc)+1):(idx_end-1));
                if(ievent>0)
                    %check whether we already had this trigger function
                    ediscf = zeros(size(trigger{ievent}));
                    for id = 1:length(trigger{ievent})
                        ediscf(id) = isequaln(abs(trigger{ievent}),abs(sym(arg)));
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
                    z{ievent} = sym.empty([0,1]);
                end
            end
        end
    end
    
    % update number of events
    nevent = length(trigger);
    
    if(nevent>0)        
        % initialis tmp_bolus
        tmp_bolus = cell([nevent,1]);
        for ievent = 1:nevent
            tmp_bolus{ievent} = sym(zeros([nx,1]));
        end
        syms polydirac
        for ix = 1:nx
            symchar = char(this.sym.xdot(ix));
            if(strfind(symchar,'dirac'))
                for ievent = 1:nevent
                    triggerchar = char(trigger{ievent});
                    str_arg_d = ['dirac(' triggerchar ')' ];
                    % replace event specific functions by "polydirac" variable
                    symchar = strrep(symchar,str_arg_d,'polydirac');
                    % extract coefficient
                    cfp = coeffs(sym(symchar),polydirac);
                    % add coefficient to bolus
                    tmp_bolus{ievent}(ix) = tmp_bolus{ievent}(ix) + cfp;
                    % remove dirac
                    symchar = strrep(symchar,'polydirac',0);
                end
                % update xdot
                this.sym.xdot(ix) = sym(symchar);
            end
        end
        % multiply by the dtriggerdt factor, this should stay here as we
        % want the xdot to be cleaned of any dirac functions
        for ievent = 1:nevent
            dtriggerdt = diff(trigger{ievent},'t') + jacobian(trigger{ievent},model.sym.x)*model.sym.xdot;
            bolus{ievent} = bolus{ievent} + tmp_bolus{ievent}/abs(dtriggerdt);
        end 
        
        % update events
        for ievent = 1:nevent
            this.event(ievent) = amievent(trigger{ievent},bolus{ievent},z{ievent});
        end
        
    end
    
end

