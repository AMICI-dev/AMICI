%
% @file amioption
% @brief definition of amioption class
%
classdef amioption < matlab.mixin.SetGet
    %AMIOPTION provides an option container to pass simulation parameters to the
    %simulation routine.
    
    properties
        atol = 1e-16;
        rtol = 1e-8;
        maxsteps = 1e4;
        sens_ind@double;
        qpositivex@double;
        tstart = 0;
        lmm = 2;
        iter = 2;
        linsol = 9;
        stldet = true;
        interpType = 1;
        lmmB = 2;
        iterB = 2;
        ism = 1;
        sensi_meth = 1;
        sensi = 0;
        nmaxevent = 10;
        ordering = 1;
        ss = 0;
        sx0@double;
    end
    
    properties (Hidden)
        z2event@double;
        id@double;
    end
    
    methods
        function obj = amioption(varargin)
            %amioptions Construct a new amioptions object
            %
            %   OPTS = amioption() creates a set of options with each option set to its
            %   default value.
            %
            %   OPTS = amioption(PARAM, VAL, ...) creates a set of options with the named
            %   parameters altered with the specified values.
            %
            %   OPTS = amioption(OLDOPTS, PARAM, VAL, ...) creates a copy of OLDOPTS with
            %   the named parameters altered with the specified value
            %
            %   Note to see the parameters, check the
            %   documentation page for amioptions
            
            % adapted from SolverOptions
            
            if nargin > 0 
                
                % Deal with the case where the first input to the
                % constructor is a amioptions/struct object.
                if isa(varargin{1},'amioption')
                    if strcmp(class(varargin{1}),class(obj))
                        obj = varargin{1};
                    else
                        % Get the properties from options object passed
                        % into the constructor.
                        thisProps = properties(obj);
                        % Set the common properties. Note that we
                        % investigated first finding the properties that
                        % are common to both objects and just looping over
                        % those. We found that in most cases this was no
                        % quicker than just looping over the properties of
                        % the object passed in.
                        for i = 1:length(thisProps)
                            try %#ok
                                % Try to set one of the properties of the
                                % old object in the new one.
                                obj.(thisProps{i}) = varargin{1}.(thisProps{i});
                            end
                        end
                    end
                    firstInputObj = true;
                elseif isstruct(varargin{1})
                    fieldlist = fieldnames(varargin{1});
                    for ifield = 1:length(fieldlist)
                        obj.(fieldlist{ifield}) = varargin{1}.(fieldlist{ifield});
                    end
                    firstInputObj = true;
                elseif isempty(varargin{1})
                    firstInputObj = true;
                else
                    firstInputObj = false;
                end
                
                % Extract the options that the caller of the constructor
                % wants to set.
                if firstInputObj
                    pvPairs = varargin(2:end);
                else
                    pvPairs = varargin;
                end
                
                % Loop through each param-value pair and just try to set
                % the option. When the option has been fully specified with
                % the correct case, this is fast. The catch clause deals
                % with partial matches or errors.
                haveCreatedInputParser = false;
                for i = 1:2:length(pvPairs)
                    try
                        obj.(pvPairs{i}) = pvPairs{i+1};
                    catch ME %#ok
                        
                        % Create the input parser if we haven't already. We
                        % do it here to avoid creating it if possible, as
                        % it is slow to set up.
                        if ~haveCreatedInputParser
                            ip = inputParser;
                            % Structures are currently not supported as
                            % an input to optimoptions. Setting the
                            % StructExpand property of the input parser to
                            % false, forces the parser to treat the
                            % structure as a single input and not a set of
                            % param-value pairs.
                            ip.StructExpand =  false;
                            % Get list of option names
                            allOptionNames = properties(obj);
                            for j = 1:length(allOptionNames)
                                % Just specify an empty default as we already have the
                                % defaults in the options object.
                                ip.addParameter(allOptionNames{j}, []);
                            end
                            haveCreatedInputParser = true;
                        end
                        
                        % Get the p-v pair to parse.
                        thisPair = pvPairs(i:min(i+1, length(pvPairs)));
                        ip.parse(thisPair{:});
                        
                        % Determine the option that was specified in p-v pairs.
                        % These options will now be matched even if only partially
                        % specified (by 13a). Now set the specified value in the
                        % options object.
                        optionSet = setdiff(allOptionNames, ip.UsingDefaults);
                        obj.(optionSet{1}) = ip.Results.(optionSet{1});
                    end
                end
            end
            
            
        end
        
        function set.sensi_meth(this,value)
            if(ischar(value))
                switch(value)
                    case 'forward'
                        this.sensi_meth = 1;
                    case 'adjoint' 
                        this.sensi_meth = 2;
                    case 'ss'
                        this.sensi_meth = 3;
                    otherwise
                        error('Unknown sensitivity method. Must be either ''forward'',''adjoint'' or ''ss''!');
                end
            else
                assert(isnumeric(value),'The option sensi_meth must have a numeric value!')
                assert(floor(value)==value,'The option sensi_meth must be an integer!')
                assert(value<=4,'Only 1, 2, 3 are valid options for sensi_meth!')
                assert(value>=0,'Only 1, 2, 3 are valid options for sensi_meth!')
                this.sensi_meth = value;
            end
        end
        
        function set.sensi(this,value)
            assert(isnumeric(value),'The option sensi must have a numeric value!')
            assert(floor(value)==value,'The option sensi must be an integer!')
            assert(value<=4,'Only 0, 1, 2 are valid options for sensi!')
            this.sensi = value;
        end
    end
    
end

