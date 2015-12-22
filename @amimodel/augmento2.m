function [this] = augmento2(this)
    % augmento2 augments the system equation to also include equations for
    % sensitivity equation. This will enable us to compute second order
    % sensitivities in a forward-adjoint or forward-forward apporach later on.
    %
    % Parameters:
    %
    % Return values:
    %  this: augmented system which contains symbolic definition of the
    %  original system and its sensitivities @type amimodel
    
    syms Sx Sdot Sy S0
   
    
    this.nxtrue = length(this.sym.x); % number of states
    this.nytrue = length(this.sym.y); % number of observables
    this.nztrue = length([this.event.z]); % number of observables
    
    % augment states
    Sx = sym(zeros(length(this.sym.x),length(this.sym.p)));
    for j = 1:length(this.sym.x)
        for k = 1:length(this.sym.p)
            eval(['syms S' num2str(j) '_' num2str(k)]);
            eval(['Sx(j,k) = S' num2str(j) '_' num2str(k) ';']);
        end
    end
    Sdot = jacobian(this.sym.xdot,this.sym.x)*Sx+jacobian(this.sym.xdot,this.sym.p);
    
    % augment output
    Sy = jacobian(this.sym.y,this.sym.x)*Sx+jacobian(this.sym.y,this.sym.p);
    
    for ievent = 1:length(this.event);
        Sz = jacobian(this.event(ievent).z,this.sym.x)*Sx+jacobian(this.event(ievent).z,this.sym.p);
        this.event(ievent).z = [this.event(ievent).z;reshape(Sz,[numel(Sz),1])];
    end
    
    S0 = jacobian(this.sym.x0,this.sym.p);
    
    this.sym.x = [this.sym.x;reshape(Sx,[numel(Sx),1])];
    this.sym.xdot = [this.sym.xdot;reshape(Sdot,[numel(Sdot),1])];
    this.sym.f = this.sym.xdot;
    this.sym.y = [this.sym.y;reshape(Sy,[numel(Sy),1])];
    this.sym.x0 = [this.sym.x0;reshape(S0,[numel(S0),1])];
    
    this.modelname = [this.modelname '_o2'];
    
    if(~exist(fullfile(this.wrap_path,'models'),'dir'))
        mkdir(fullfile(this.wrap_path,'models'));
        mkdir(fullfile(this.wrap_path,'models',this.modelname));
    else
        if(~exist(fullfile(this.wrap_path,'models',this.modelname),'dir'))
            mkdir(fullfile(this.wrap_path,'models',this.modelname))
        end
    end
    
end


