%
% @file amidata
% @brief definition of amidata class
%
classdef amidata < handle
    %AMIDATA provides a data container to pass experimental data to the
    %simulation routine for likelihood computation. 
    
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
        t@double;
        % observations
        Y@double;
        % standard deviation of observations
        Sigma_Y@double;
        % event observations
        Z@double;
        % standard deviation of event observations
        Sigma_Z@double;
    end
    
    methods
        function D = amidata(nt,ny,nz,ne,nk)
            D.nt = nt;
            D.ny = ny;
            D.nz = nz;
            D.ne = ne;
            D.nk = nk;
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
        
        function set.t(this,value)
            assert(isnumeric(value),'AMICI:amimodel:t:numeric','t must have a numeric value!')
            assert(ismatrix(value),'AMICI:amimodel:t:ndims','t must be a two dimensional matrix!')
            assert(numel(value)==this.nt,'AMICI:amimodel:t:ndims',['t must have ' num2str(this.nt) ' (D.nt) elements!'])
            this.t = double(value(:));
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

