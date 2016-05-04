%
% @file amievent
% @brief definition of amievent class
%
classdef amievent
    % AMIEVENT defines events which later on will be transformed into appropriate
    % C code
    
    properties ( GetAccess = 'public', SetAccess = 'private' )
        % the trigger function activates the event on every zero crossing @type symbolic
        trigger = sym.empty();
        % the bolus function defines the change in states that is applied on every event occurence @type symbolic
        bolus = sym.empty();
        % output function for the event @type symbolic
        z = sym.empty();
        % flag indicating that a heaviside function is present, this helps
        % to speed up symbolic computations
        hflag = logical.empty();
    end
    
    methods
        function AE = amievent(trigger,bolus,z)
            % amievent constructs an amievent object from the provided input.
            %
            % Parameters:
            %  trigger: trigger function, the event will be triggered on
            %      at all roots of this function
            %  bolus: the bolus that will be added to all states on every
            %      occurence of the event
            %  z: the event output that will be reported on every occurence
            %      of the event
            %
            % Return values:
            %  AE: amievent object
            %
            if(~isa(trigger,'sym'))
                if(isa(trigger,'double'))
                    AE.trigger = sym(trigger);
                    warning('Constant trigger function will never trigger. Please check the event definition.')
                else
                    error('trigger function must be a symbolic expression')
                end
            else
                AE.trigger = trigger;
            end
            if(numel(AE.trigger)>1)
                error('The trigger function must be scalar.')
            end
            
            if(~isa(bolus,'sym'))
                if(isa(bolus,'double'))
                    AE.bolus = sym(bolus(:));
                else
                    error('bolus function must be a symbolic expression')
                end
            else
                AE.bolus = bolus;
            end
            
            if(~isa(z,'sym'))
                if(isa(z,'double'))
                    if(~isempty(z))
                        AE.z = sym(z);
                        warning('Constant outputs are not informative. Please check the event definition.')
                    else
                        AE.z = sym(zeros(0,1));
                    end
                else
                    error('output function must be a symbolic expression')
                end      
            else
                AE.z = z;
            end
        end
        
        this = setHflag(this,hflag);
    end
end

