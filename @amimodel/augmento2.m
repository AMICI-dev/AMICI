function [modelo2] = augmento2(this)
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
   
    augmodel.nxtrue = length(this.sym.x); % number of states
    augmodel.nytrue = length(this.sym.y); % number of observables
    if(this.nevent>0)
        augmodel.nztrue = length([this.event.z]); % number of observables
    else
        augmodel.nztrue = 0;
    end
    np = this.np;
    
    % augment states
    Sx = sym(zeros(length(this.sym.x),np));
    for j = 1:length(this.sym.x)
        for k = 1:np
            eval(['syms S' num2str(j) '_' num2str(k)]);
            eval(['Sx(j,k) = S' num2str(j) '_' num2str(k) ';']);
        end
    end
    Sdot = jacobian(this.sym.xdot,this.sym.x)*Sx+jacobian(this.sym.xdot,this.sym.p);
    
    % augment output
    Sy = jacobian(this.sym.y,this.sym.x)*Sx+jacobian(this.sym.y,this.sym.p);
   
    % generate deltasx
    this.getFun([],'deltasx');
    for ievent = 1:this.nevent;
        if(numel(this.event(ievent).z)>0)
            Sz = jacobian(this.event(ievent).z,this.sym.x)*Sx+jacobian(this.event(ievent).z,this.sym.p);
            znew = [this.event(ievent).z,reshape(Sz,[1,numel(Sz)])];
        else
            znew = this.event(ievent).z;
        end
        tmp=subs(this.fun.deltasx.sym(:,:,ievent),this.fun.xdot.strsym_old,this.fun.xdot.sym);
        tmp=subs(tmp,this.fun.xdot.strsym,subs(this.fun.xdot.sym,this.fun.x.sym,this.fun.x.sym+this.event(ievent).bolus));
        tmp=subs(subs(tmp,this.fun.stau.strsym,this.fun.stau.sym),this.fun.sx.sym, Sx);
        bolusnew = [this.event(ievent).bolus;reshape(tmp,[numel(Sx),1])];
        % replace sx by augmented x
        for ip = 1:np
            bolusnew(this.nxtrue*ip+(1:this.nxtrue)) = mysubs(bolusnew(this.nxtrue*ip+(1:this.nxtrue)), this.fun.sx.sym(:,ip),Sx(:,ip));
        end
        augmodel.event(ievent) = amievent(this.event(ievent).trigger,bolusnew,znew);
    end
    
    % augment likelihood
    this.getFun([],'dsigma_ydp');
    this.getFun([],'y');
    this.getFun([],'dydp');
    SJy = jacobian(this.sym.Jy,this.sym.p) ...
        + jacobian(this.sym.Jy,this.fun.sigma_y.strsym)*this.fun.dsigma_ydp.sym ...
        + jacobian(this.sym.Jy,this.fun.y.strsym)*Sy;
    this.getFun([],'dsigma_zdp');
    this.getFun([],'z');
    this.getFun([],'dzdp');   
    SJz = jacobian(this.sym.Jz,this.sym.p);
    if(~isempty(this.fun.sigma_z.strsym))
        SJz = SJz + jacobian(this.sym.Jz,this.fun.sigma_z.strsym)*this.fun.dsigma_zdp.sym ...
              + jacobian(this.sym.Jz,this.fun.z.strsym)*Sz;   
    end
    
    % augment sigmas
    this.getFun([],'sigma_y');
    this.getFun([],'sigma_z');
    
    S0 = jacobian(this.sym.x0,this.sym.p);
    
    augmodel.sym.x = [this.sym.x;reshape(Sx,[numel(Sx),1])];
    augmodel.sym.xdot = [this.sym.xdot;reshape(Sdot,[numel(Sdot),1])];
    augmodel.sym.f = augmodel.sym.xdot;
    augmodel.sym.y = [this.sym.y;reshape(Sy,[numel(Sy),1])];
    augmodel.sym.x0 = [this.sym.x0;reshape(S0,[numel(S0),1])];
    augmodel.sym.Jy = [this.sym.Jy,SJy];
    augmodel.sym.Jz = [this.sym.Jz,SJz];
    augmodel.sym.p = this.sym.p;
    augmodel.sym.k = this.sym.k;
    augmodel.sym.sigma_y = [this.sym.sigma_y(:)', reshape(transpose(this.fun.dsigma_ydp.sym), [1,numel(this.fun.dsigma_ydp.sym)])];
    augmodel.sym.sigma_z = [this.sym.sigma_z(:)', reshape(transpose(this.fun.dsigma_zdp.sym), [1,numel(this.fun.dsigma_zdp.sym)])];
    
    modelo2 = amimodel(augmodel,[this.modelname '_o2']);
    modelo2.o2flag = 1;
    modelo2.debug = this.debug;
    modelo2.forward = this.forward;
    modelo2.adjoint = this.adjoint;
end


function out = mysubs(in, old, new)
    % mysubs is a wrapper for ther subs matlab function
    %
    % Parameters:
    %  in: symbolic expression in which to replace @type symbolic
    %  old: expression to be replaced @type symbolic
    %  new: replacement expression @type symbolic
    %
    % Return values:
    %  out: symbolic expression with replacement @type symbolic
    if(~isnumeric(in) && ~isempty(old) && ~isempty(symvar(in)))
        matVer = ver('MATLAB');
        if(str2double(matVer.Version)>=8.1)
            out = subs(in, old(:), new(:));
        else
            out = subs(in, old(:), new(:), 0);
        end
    else
        out = in;
    end
end
