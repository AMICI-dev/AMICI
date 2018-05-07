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
    this.getFun([],'sz');
    this.getFun([],'srz');
    for ievent = 1:this.nevent;
        if(numel(this.event(ievent).z)>0)
            Sz = this.fun.sz.sym(this.z2event==ievent,:);
            for ip = 1:np
                Sz(:,ip) = subs(Sz(:,ip),this.fun.sx.sym(:,ip),Sx(:,ip));
            end
            Srz = this.fun.srz.sym(this.z2event==ievent,:);
            for ip = 1:np
                Srz(:,ip) = subs(Srz(:,ip),this.fun.sx.sym(:,ip),Sx(:,ip));
            end
            znew = [this.event(ievent).z,reshape(Sz,[sum(this.z2event==ievent),np])];
        else
            znew = this.event(ievent).z;
        end
        tmp=subs(this.fun.deltasx.sym(:,:,ievent),this.fun.xdot.strsym_old,this.fun.xdot.sym);
        tmp=subs(tmp,this.fun.xdot.strsym,subs(this.fun.xdot.sym,this.fun.x.sym,this.fun.x.sym+this.event(ievent).bolus));
        tmp=subs(tmp,this.fun.stau.strsym,this.fun.stau.sym(ievent,:));
        for ip = 1:np
            tmp(:,ip)=subs(tmp(:,ip),this.fun.sx.sym(:,ip), Sx(:,ip));
        end
        bolusnew = [this.event(ievent).bolus;tmp(:)];
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
    this.getFun([],'Jy');
    tmp = arrayfun(@(x) sym(['var_y_' num2str(x)]),0:(augmodel.nytrue*(1+np)-1),'UniformOutput',false);
    tmp = transpose([tmp{:}]);
    aug_y_strsym =  reshape(tmp((augmodel.nytrue+1):end),[augmodel.nytrue,np]);% the update Jy must not contain any
    % references to x (otherwise we have to deal with sx in sJy), thus
    % we replace sy with the corresponding augmented y that we are about to
    % create
    SJy = jacobian(this.fun.Jy.sym,this.sym.p) ...
        + jacobian(this.fun.Jy.sym,this.fun.sigma_y.strsym)*this.fun.dsigma_ydp.sym ...
        + jacobian(this.fun.Jy.sym,this.fun.y.strsym)*aug_y_strsym;
    
    this.getFun([],'dsigma_zdp');
    this.getFun([],'rz');
    this.getFun([],'Jz');
    tmp = arrayfun(@(x) sym(['var_z_' num2str(x)]),0:(augmodel.nztrue*(1+np)-1),'UniformOutput',false);
    tmp = transpose([tmp{:}]);
    aug_z_strsym =  reshape(tmp((augmodel.nztrue+1):end),[augmodel.nztrue,np]);% the update Jz must not contain any
    % references to x (otherwise we have to deal with sx in sJy), thus
    % we replace sy with the corresponding augmented y that we are about to
    % create
    SJz = jacobian(this.fun.Jz.sym,this.sym.p);
    if(~isempty(this.fun.sigma_z.strsym))
        SJz = SJz + jacobian(this.fun.Jz.sym,this.fun.sigma_z.strsym)*this.fun.dsigma_zdp.sym ...
              + jacobian(this.fun.Jz.sym,this.fun.z.strsym)*aug_z_strsym;   
    end
    this.getFun([],'Jrz');
    tmp = arrayfun(@(x) sym(['var_rz_' num2str(x)]),0:(augmodel.nztrue*(1+np)-1),'UniformOutput',false);
    tmp = transpose([tmp{:}]);
    aug_rz_strsym =  reshape(tmp((augmodel.nztrue+1):end),[augmodel.nztrue,np]);% the update Jz must not contain any
    % references to x (otherwise we have to deal with sx in sJy), thus
    % we replace sy with the corresponding augmented y that we are about to
    % create
    SJrz = jacobian(this.fun.Jrz.sym,this.sym.p);
    if(~isempty(this.fun.sigma_z.strsym))
        SJrz = SJrz + jacobian(this.fun.Jrz.sym,this.fun.sigma_z.strsym)*this.fun.dsigma_zdp.sym ...
              + jacobian(this.fun.Jrz.sym,this.fun.rz.strsym)*aug_rz_strsym;   
    end
    
    % augment sigmas
    this.getFun([],'sigma_y');
    this.getFun([],'sigma_z');
    S0 = jacobian(this.sym.x0,this.sym.p);
    
    augmodel.sym.x = [this.sym.x;Sx(:)];
    augmodel.sym.xdot = [this.sym.xdot;Sdot(:)];
    augmodel.sym.f = augmodel.sym.xdot;
    augmodel.sym.y = [this.sym.y;Sy(:)];
    augmodel.sym.x0 = [this.sym.x0;S0(:)];
    augmodel.sym.Jy = [this.sym.Jy,SJy];
    augmodel.sym.Jz = [this.sym.Jz,SJz];
    augmodel.sym.Jrz = [this.sym.Jrz,SJrz];
    if(numel([this.event.z]))
        augmodel.sym.rz = [this.fun.rz.sym,Srz];
    end
    augmodel.sym.p = this.sym.p;
    augmodel.sym.k = this.sym.k;
    augmodel.sym.sigma_y = [transpose(this.sym.sigma_y(:)), reshape(transpose(this.fun.dsigma_ydp.sym), [1,numel(this.fun.dsigma_ydp.sym)])];
    augmodel.sym.sigma_z = [transpose(this.sym.sigma_z(:)), reshape(transpose(this.fun.dsigma_zdp.sym), [1,numel(this.fun.dsigma_zdp.sym)])];
    
    modelo2 = amimodel(augmodel,[this.modelname '_o2']);
    modelo2.o2flag = 1;
    modelo2.debug = this.debug;
    modelo2.forward = this.forward;
    modelo2.adjoint = this.adjoint;
    modelo2.param = this.param;
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
