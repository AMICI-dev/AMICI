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
   
    
    this.nxtrue = length(this.sym.x); %! number of states
    this.nytrue = length(this.sym.y); %! number of observables
    
    
    Sx = sym(zeros(length(this.sym.x),length(this.sym.p)));
    
    for j = 1:length(this.sym.x)
        for k = 1:length(this.sym.p)
            eval(['syms S' num2str(j) '_' num2str(k)]);
            eval(['Sx(j,k) = S' num2str(j) '_' num2str(k) ';']);
        end
    end
    
    Sdot = jacobian(this.sym.xdot,this.sym.x)*Sx+jacobian(this.sym.xdot,this.sym.p);
    
    Sy = jacobian(this.sym.y,this.sym.x)*Sx+jacobian(this.sym.y,this.sym.p);
    
    if(size(this.sym.y,2)>size(this.sym.y,1))
        this.sym.y = transpose(this.sym.y);
    end
    if(size(this.sym.x,2)>size(this.sym.x,1))
        this.sym.x = transpose(this.sym.x);
    end
    if(size(this.sym.xdot,2)>size(this.sym.xdot,1))
        this.sym.xdot = transpose(this.sym.xdot);
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


