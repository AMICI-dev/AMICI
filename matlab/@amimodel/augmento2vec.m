function [modelo2vec] = augmento2vec(this)
    % augmento2 augments the system equation to also include equations for
    % sensitivity equation. This will enable us to compute second order
    % sensitivities in a forward-adjoint or forward-forward apporach later on.
    %
    % Parameters:
    %
    % Return values:
    %  modelo2vec: augmented system which contains symbolic definition of the
    %  original system and its sensitivities @type amimodel
    
    syms Sx Sdot Sy S0
    

   
    augmodel.nxtrue = length(this.sym.x); % number of states
    augmodel.nytrue = length(this.sym.y); % number of observables
    augmodel.nztrue = this.nz;
    augmodel.coptim = this.coptim;
    augmodel.debug = this.debug;
    
    % multiplication vector (extension of kappa
    vecs = cell([length(this.sym.p),1]);
    for ivec = 1:length(this.sym.p)
        vecs{ivec} = sprintf('k_%i', length(this.sym.k) + ivec-1);
    end
    vec = sym(vecs);
    
    if(this.nevent>0)
        augmodel.nztrue = length([this.event.z]); % number of observables
    else
        augmodel.nztrue = 0;
    end
    np = this.np;
    
    % augment states
    sv = sym('sv',[length(this.sym.x),1]);
    Sdot = jacobian(this.sym.xdot,this.sym.x)*sv+jacobian(this.sym.xdot,this.sym.p)*vec;
    
    % augment output
    Sy = jacobian(this.sym.y,this.sym.x)*sv+jacobian(this.sym.y,this.sym.p)*vec;
   
    % generate deltasx
    this.getFun([],'deltasx');
    for ievent = 1:this.nevent;
        if(numel(this.event(ievent).z)>0)
            Sz = jacobian(this.event(ievent).z,this.sym.x)*sv+jacobian(this.event(ievent).z,this.sym.p)*vec;
            znew = [this.event(ievent).z,Sz];
        else
            znew = this.event(ievent).z;
        end
        tmp=subs(this.fun.deltasx.sym(:,:,ievent),this.fun.xdot.strsym_old,this.fun.xdot.sym);
        tmp=subs(tmp,this.fun.xdot.strsym,subs(this.fun.xdot.sym,this.fun.x.sym,this.fun.x.sym+this.event(ievent).bolus));
        tmp=subs(subs(tmp,this.fun.stau.strsym,this.fun.stau.sym(ievent,:)),this.fun.sx.sym*vec,sv);
        bolusnew = [this.event(ievent).bolus;tmp*vec];
        % replace sx by augmented x
        bolusnew(this.nxtrue+(1:this.nxtrue)) = mysubs(bolusnew(this.nxtrue+(1:this.nxtrue)), this.fun.sx.sym(:,1),sv);

        hflagold = this.event(ievent).hflag;
        augmodel.event(ievent) = amievent(this.event(ievent).trigger,bolusnew,znew);
        augmodel.event(ievent) = augmodel.event(ievent).setHflag([hflagold;zeros([numel(sv),1])]);
    end
    
    % augment likelihood
    this.getFun([],'dsigma_ydp');
    this.getFun([],'y');
    this.getFun([],'dydp');
    this.getFun([],'Jy');
    tmp = arrayfun(@(x) sym(['var_y_' num2str(x)]),0:(augmodel.nytrue*2-1),'UniformOutput',false);
    aug_y_strsym = transpose([tmp{(augmodel.nytrue+1):end}]); % the update Jy must not contain any
    % references to x (otherwise we have to deal with sx in sJy), thus
    % we replace sy with the corresponding augmented y that we are about to
    % create
    SJy = (jacobian(this.fun.Jy.sym,this.sym.p) ...
        + jacobian(this.fun.Jy.sym,this.fun.sigma_y.strsym)*this.fun.dsigma_ydp.sym) ...
        * vec + jacobian(this.fun.Jy.sym,this.fun.y.strsym)*aug_y_strsym;
    this.getFun([],'dsigma_zdp');
    
    this.getFun([],'dzdp');
    this.getFun([],'Jz');
    SJz = jacobian(this.fun.Jz.sym,this.sym.p);
    if(~isempty(this.fun.sigma_z.strsym))
        SJz = SJz + jacobian(this.fun.Jz.sym,this.fun.sigma_z.strsym)*this.fun.dsigma_zdp.sym ...
              + jacobian(this.fun.Jz.sym,this.fun.z.strsym)*Sz;
    end

    % augment sigmas
    this.getFun([],'sigma_y');
    this.getFun([],'sigma_z');
    
    S0 = jacobian(this.sym.x0,this.sym.p)*vec;
    
    augmodel.sym.x = [this.sym.x;sv];
    augmodel.sym.xdot = [this.sym.xdot;Sdot];
    augmodel.sym.f = augmodel.sym.xdot;
    augmodel.sym.y = [this.sym.y;Sy];
    augmodel.sym.x0 = [this.sym.x0;S0];
    augmodel.sym.Jy = [this.sym.Jy,SJy];
    augmodel.sym.Jz = [this.sym.Jz,SJz*vec];
    augmodel.sym.k = [this.sym.k;vec];
    augmodel.sym.p = this.sym.p;
    augmodel.sym.sigma_y = [this.sym.sigma_y, transpose(this.fun.dsigma_ydp.sym * vec)];
    augmodel.sym.sigma_z = [this.sym.sigma_z, transpose(this.fun.dsigma_zdp.sym * vec)];
    
    modelo2vec = amimodel(augmodel,[this.modelname '_o2vec']);
    modelo2vec.o2flag = 2;
    modelo2vec.debug = this.debug;
    modelo2vec.forward = this.forward;
    modelo2vec.adjoint = this.adjoint;
    modelo2vec.param = this.param;
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
