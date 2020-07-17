%
% @file amidata
% @brief definition of amidata class
%
classdef amidata < handle
    % AMIDATA provides a data container to pass experimental data to the
    % simulation routine for likelihood computation. 
    % when any of the properties are updated, the class automatically 
    % checks consistency of dimension and updates related properties and 
    % initialises them with NaNs
    
    properties
        % number of timepoints
        nt=0;
        % number of observables
        ny=0;
        % number of event observables
        nz=0;
        % number of events
        ne=0;
        % number of conditions/constants
        nk=0;
        % timepoints of observations
        t = double.empty();
        % observations
        Y = double.empty();
        % standard deviation of observations
        Sigma_Y = double.empty();
        % event observations
        Z = double.empty();
        % standard deviation of event observations
        Sigma_Z = double.empty();
        % experimental condition
        condition = double.empty();
        % experimental condition for preequilibration
        conditionPreequilibration = double.empty();
        % reinitialize states based on fixed parameters after preeq.?
        reinitializeStates = false;
    end
    
    methods
        function D = amidata(varargin)
            % amidata creates an amidata container for experimental data 
            % with specified dimensions amidata.
            % 
            % AMIDATA(amidata) creates a copy of the input container
            %
            % AMIDATA(struct) tries to creates an amidata container from the
            %                 input struct. the struct should have the following
            %                 fields:
            %                 t [nt,1]
            %                 Y [nt,ny]
            %                 Sigma_Y [nt,ny]
            %                 Z [ne,nz]
            %                 Sigma_Z [ne,nz]
            %                 condition [nk,1]
            %                 conditionPreequilibration [nk,1]
            %                 if some fields are missing the function will try
            %                 to initialise them with NaNs with consistent
            %                 dimensions
            %
            % AMIDATA(nt,ny,nz,ne,nk) constructs an empty data container with
            %                 in the provided dimensions intialised with NaNs
            %
            % 
            %
            % Parameters:
            %    varargin:
            %
            % Return values:
            %
            
            
            % initialisation via struct
            if isa(varargin{1},'amidata')
                    if strcmp(class(varargin{1}),class(D))
                        D = varargin{1};
                    else
                        thisProps = properties(D);
                        for i = 1:length(thisProps)
                            try %#ok
                                % Try to set one of the properties of the
                                % old object in the new one.
                                D.(thisProps{i}) = varargin{1}.(thisProps{i});
                            end
                        end
                    end
            elseif(isstruct(varargin{1}))
                if(isfield(varargin{1},'t'))
                    D.nt = numel(varargin{1}.t);
                    D.t = varargin{1}.t(:);
                elseif(isfield(varargin{1},'Y'))
                    D.nt = size(varargin{1}.Y,1);
                else
                    error('Cannot construct valid amidata object from input struct. There is no field D.t or D.Y in the input struct!');
                end
                if(isfield(varargin{1},'Y'))
                    D.ny = size(varargin{1}.Y,2);
                    D.Y = varargin{1}.Y;
                else
                    D.ny = 0;
                end
                if(isfield(varargin{1},'Z'))
                    D.nz = size(varargin{1}.Z,2);
                    D.ne = size(varargin{1}.Z,1);
                    D.Z = varargin{1}.Z;
                else
                    D.nz = 0;
                end
                if(isfield(varargin{1},'Sigma_Y'))
                    D.Sigma_Y = varargin{1}.Sigma_Y;
                end
                if(isfield(varargin{1},'Sigma_Z'))
                    D.Sigma_Z = varargin{1}.Sigma_Z;
                end   
                if(isfield(varargin{1},'condition'))
                    D.nk = numel(varargin{1}.condition);
                    D.condition = varargin{1}.condition;
                else
                    D.nk = 0;
                end
                if(isfield(varargin{1},'conditionPreequilibration'))
                    assert(D.nk == numel(varargin{1}.conditionPreequilibration));
                    D.conditionPreequilibration = varargin{1}.conditionPreequilibration;
                end
                if(isfield(varargin{1},'reinitializeStates'))
                    if islogical(varargin{1}.reinitializeStates)
                        D.reinitializeStates = varargin{1}.reinitializeStates;
                    elseif isnumeric(varargin{1}.reinitializeStates)
                        D.reinitializeStates = logical(varargin{1}.reinitializeStates);
                    else
                        error('Assignment error: Value for field reinitializeStates must be logical.');
                    end
                end    
            elseif(nargin == 5)
                D.nt = varargin{1};
                D.ny = varargin{2};
                D.nz = varargin{3};
                D.ne = varargin{4};
                D.nk = varargin{5};
            end

        end
        
        function set.nt(this,nt)
            this.nt = nt;
            this.t = 1:nt;
            this.Y = NaN;
            this.Sigma_Y = NaN;
        end
        
        function set.ny(this,ny)
            this.ny = ny;
            this.Y = NaN;
            this.Sigma_Y = NaN;
        end
        
        function set.nz(this,nz)
            this.nz = nz;
            this.Z = NaN;
            this.Sigma_Z = NaN;
        end
        
        function set.ne(this,ne)
            this.ne = ne;
            this.Z = NaN;
            this.Sigma_Z = NaN;
        end
        
        function set.nk(this,nk)
            this.nk = nk;
            this.condition = NaN(nk,1);
        end
        
        function set.t(this,value)
            assert(isnumeric(value),'AMICI:amimodel:t:numeric','t must have a numeric value!')
            assert(ismatrix(value),'AMICI:amimodel:t:ndims','t must be a two dimensional matrix!')
            assert(numel(value)==this.nt,'AMICI:amimodel:t:ndims',['t must have ' num2str(this.nt) ' (D.nt) elements!'])
            this.t = double(value(:));
        end
        
        function set.condition(this,value)
            assert(isnumeric(value),'AMICI:amimodel:condition:numeric','condition must have a numeric value!')
            assert(ismatrix(value),'AMICI:amimodel:condition:ndims','condition must be a two dimensional matrix!')
            assert(numel(value)==this.nk,'AMICI:amimodel:condition:ndims',['condition must have ' num2str(this.nk) ' (D.nk) elements!'])
            this.condition = double(value(:));
        end
        
        function set.conditionPreequilibration(this,value)
            assert(isnumeric(value),'AMICI:amimodel:condition:numeric','condition must have a numeric value!')
            assert(ismatrix(value),'AMICI:amimodel:condition:ndims','condition must be a two dimensional matrix!')
            assert(numel(value)==this.nk,'AMICI:amimodel:condition:ndims',['condition must have ' num2str(this.nk) ' (D.nk) elements!'])
            this.conditionPreequilibration = double(value(:));
        end
        
        function set.Y(this,value)
            assert(ismatrix(value),'AMICI:amimodel:Y:ndims','Y must be a two dimensional matrix!')
            assert(all(all(or(isnumeric(value),isnan(value)))),'AMICI:amimodel:Y:numeric','Y must have a numeric value!')
            
            if(all(size(value)==[this.nt this.ny]))
                this.Y = double(value);
            elseif(all(size(value)==[this.nt 1]))
                this.Y = repmat(double(value),[1,this.ny]);
            elseif(all(size(value)==[1 this.ny]))
                this.Y = repmat(double(value),[this.nt,1]);
            elseif(all(size(value)==[1 1]))
                this.Y = repmat(double(value),[this.nt,this.ny]);
            else
                error('AMICI:amimodel:Y:size',['Y must have size [' num2str(this.nt) ',' num2str(this.ny) '] ([D.nt,D.ny])!'])
            end
        end
        
        function set.Sigma_Y(this,value)
            assert(ismatrix(value),'AMICI:amimodel:Sigma_Y:ndims','Sigma_Y must be a two dimensional matrix!')
            assert(all(all(or(isnumeric(value),isnan(value)))),'AMICI:amimodel:Sigma_Y:numeric','Sigma_Y must have a numeric value!')
            if(all(size(value)==[this.nt this.ny]))
                this.Sigma_Y = double(value);
            elseif(all(size(value)==[this.nt 1]))
                this.Sigma_Y = repmat(double(value),[1,this.ny]);
            elseif(all(size(value)==[1 this.ny]))
                this.Sigma_Y = repmat(double(value),[this.nt,1]);
            elseif(all(size(value)==[1 1]))
                this.Sigma_Y = repmat(double(value),[this.nt,this.ny]);
            else
                error('AMICI:amimodel:Sigma_Y:size',['Sigma_Y must have size [' num2str(this.nt) ',' num2str(this.ny) '] ([D.nt,D.ny])!'])
            end
        end
        
        function set.Z(this,value)
            assert(ismatrix(value),'AMICI:amimodel:Z:ndims','Z must be a two dimensional matrix!')
            assert(all(all(or(isnumeric(value),isnan(value)))),'AMICI:amimodel:Z:numeric','Z must have a numeric value!')
            if(all(size(value)==[this.ne this.nz]))
                this.Z = double(value);
            elseif(all(size(value)==[this.ne 1]))
                this.Z = repmat(double(value),[1,this.nz]);
            elseif(all(size(value)==[1 this.nz]))
                this.Z = repmat(double(value),[this.ne,1]);
            elseif(all(size(value)==[1 1]))
                this.Z = repmat(double(value),[this.ne,this.nz]);
            else
                error('AMICI:amimodel:Z:size',['Z must have size [' num2str(this.ne) ',' num2str(this.nz) '] ([D.ne,D.nz])!'])
            end
        end
        
        function set.Sigma_Z(this,value)
            assert(ismatrix(value),'AMICI:amimodel:Sigma_Z:ndims','Sigma_Z must be a two dimensional matrix!')
            assert(all(all(or(isnumeric(value),isnan(value)))),'AMICI:amimodel:Sigma_Z:numeric','Sigma_Z must have a numeric value!')
            if(all(size(value)==[this.ne this.nz]))
                this.Sigma_Z = double(value);
            elseif(all(size(value)==[this.ne 1]))
                this.Sigma_Z = repmat(double(value),[1,this.nz]);
            elseif(all(size(value)==[1 this.nz]))
                this.Sigma_Z = repmat(double(value),[this.ne,1]);
            elseif(all(size(value)==[1 1]))
                this.Sigma_Z = repmat(double(value),[this.ne,this.nz]);
            else
                error('AMICI:amimodel:Sigma_Z:size',['Sigma_Z must have size [' num2str(this.ne) ',' num2str(this.nz) '] ([D.ne,D.nz])!'])
            end
        end
    end
    
end

