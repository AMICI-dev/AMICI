%
% @file amievent
% @brief definition of amievent class
%
classdef amievent
    % the amievent class defines the prototype for all events which later
    % on will be transformed into C code
    
    properties ( GetAccess = 'public', SetAccess = 'private' )
        % the trigger function activates the event on every zero crossing @type symbolic
        trigger = sym([]);
        % the bolus function defines the change in states that is applied on every event occurence @type symbolic
        bolus = sym([]);
        % output function for the event @type symbolic
        z = sym([]);
        % flag indicating that a heaviside function is present, this helps
        % to speed up symbolic computations
        hflag = logical([]);
    end
    
    methods
        function AE = amievent(trigger,bolus,z)
            % constructor of the amievent class. this function constructs an event object based on the provided
            % trigger function, bolus function and output function
            %
            % Parameters:
            %  trigger: trigger fuction, the roots of this function define
            %  the occurence of the event @type symbolic
            %  bolus: bolus fuction, this function defines the change in 
            %  the states on event occurences @type symbolic
            %  z: output function, this expression is evaluated on event
            %  occurences and returned by the simulation function @type
            %  symbolic
            % 
            % Return values:
            %  AM: model definition object
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

