function [ this ] = getSensiFlag(this)
    % getSensiFlag populates the sensiflag property
    %
    % Parameters:
    %
    % Return values:
    %  this: updated function definition object @type amifun
    
    switch(this.funstr)
        case 'dxdotdp'
            this.sensiflag = true;
            
        case 'sx0'
            this.sensiflag = true;  
            
        case 'sdx0'
            this.sensiflag = true;
            
            
        case 'dydp'
            this.sensiflag = true;
          
        case 'sy'
            this.sensiflag = true;
            
        case 'qBdot'
            this.sensiflag = true;
            
        case 'dsigma_ydp'
            this.sensiflag = true;
            
        case 'dsigma_zdp'
            this.sensiflag = true;
            
        case 'drzdp'
            this.sensiflag = true;
            
        case 'ddeltaxdp'
            this.sensiflag = true;
            
        case 'deltasx'
            this.sensiflag = true;
        
        case 'deltaqB'
            this.sensiflag = true;
            
        case 'dzdp'
            this.sensiflag = true;
            
        case 'sz'
            this.sensiflag = true;
            
        case 'dtaudp'
            this.sensiflag = true;
            
        case 'stau'
            this.sensiflag = true;
            
        case 'sroot'
            this.sensiflag = true;
            
        case 'srz'
            this.sensiflag = true;
            
        case 'sx'
            this.sensiflag = true;
            
        case 'sdx'
            this.sensiflag = true;

        otherwise
            this.sensiflag = false;
    end
end

