%
% @file amievent
% @brief definition of amievent class
%
classdef amievent
    % the amievent class defines the prototype for all events which later
    % on will be transformed into C code
    
    properties ( GetAccess = 'public', SetAccess = 'public' )
        % the trigger function activates the event on every zero crossing @type symbolic
        trigger;
        % the bolus function defines the change in states that is applied on every event occurence @type symbolic
        bolus;
        % output function for the event @type symbolic
        z;
        
    end
    
    methods
        function AE = amievent(trigger,bolus,z)
            if(~isa(trigger,'sym'))
                error('trigger function must be a symbolic expression')
            else
                AE.trigger = trigger;
            end
            
            if(~isa(bolus,'sym'))
                error('bolus function must be a symbolic expression')
            else
                AE.bolus = bolus;
            end
            
            if(~isa(z,'sym'))
                error('bolus function must be a symbolic expression')
            else
                AE.z = z;
            end
        end
    end
end

