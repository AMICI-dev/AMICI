function printLocalVars(this,model,fid)
    % printlocalvars prints the C code for the initialisation of local variables into the file specified by fid.
    %
    % Parameters:
    %  model: this struct must contain all necessary symbolic definitions
    %  @type amimodel
    %  fid: file id in which the final expression is written @type fileid
    %
    % Return values:
    %  Nothing
    
    
    nx = model.nx;
    ny = model.ny;
    
    for nvec = this.nvecs
        if(strcmp(nvec{1},'*sx'))
            fprintf(fid,'realtype *sx_tmp;\n');
        elseif(strcmp(nvec{1},'*sdx'))
            fprintf(fid,'realtype *sdx_tmp;\n');
        elseif(strcmp(nvec{1},'*sxdot'))
            fprintf(fid,'realtype *sxdot_tmp;\n');
        elseif(strcmp(nvec{1},'*sx0'))
            fprintf(fid,'realtype *sx0_tmp;\n');
        elseif(strcmp(nvec{1},'*sdx0'))
            fprintf(fid,'realtype *sdx0_tmp;\n');
        else
            if(or(strcmp(model.wtype,'iw'),~strcmp(nvec{1},'dx')))
                % check for nx ==0, otherwise we will access a NULL pointer and crash
                if(not(and(nx==0,strcmp(nvec{1},'x'))))
                    fprintf(fid,['realtype *' nvec{1} '_tmp = N_VGetArrayPointer(' nvec{1} ');\n']);
                end
            end
        end
    end
    
    if(this.sensiflag && not(strcmp(this.funstr,'sxdot')))    
        fprintf(fid,'int ip;\n');
    end
    
    switch(this.funstr)
        case 'xdot'
            fprintf(fid,'int ix;\n');
            fprintf(fid,['memset(xdot_tmp,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'xBdot'
            fprintf(fid,'int ix;\n');
            fprintf(fid,['memset(xBdot_tmp,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'qBdot'
            fprintf(fid,'memset(qBdot_tmp,0,sizeof(realtype)*udata->nplist*model->nJ);\n');
        case 'x0'
            fprintf(fid,['memset(x0_tmp,0,sizeof(realtype)*' num2str(nx) ');\n']);
            fprintf(fid,['realtype t = udata->tstart;\n']);
        case 'dx0'
            fprintf(fid,['memset(dx0_tmp,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'Jv'
            fprintf(fid,['memset(Jv_tmp,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'JvB'
            fprintf(fid,['memset(JvB_tmp,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'JBand'
            % nothing
        case 'J'
            fprintf(fid,'int ix;\n');
            fprintf(fid,['memset(J->data,0,sizeof(realtype)*' num2str(nx^2) ');\n']);
        case 'JDiag'
            fprintf(fid,'int ix;\n');
            fprintf(fid,['memset(JDiag_tmp,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'JSparse'
            fprintf(fid,'int inz;\n');
            fprintf(fid,'SparseSetMatToZero(J);\n');
            for i = 1:length(model.rowvals)
                fprintf(fid,['J->indexvals[' num2str(i-1) '] = ' num2str(model.rowvals(i)) ';\n']);
            end
            for i = 1:length(model.colptrs)
                fprintf(fid,['J->indexptrs[' num2str(i-1) '] = ' num2str(model.colptrs(i)) ';\n']);
            end
        case 'JBandB'
            % nothing
        case 'JB'
            fprintf(fid,['  memset(JB->data,0,sizeof(realtype)*' num2str(nx^2) ');\n']);
        case 'JSparseB'
            fprintf(fid,'  SparseSetMatToZero(JB);\n');
            for i = 1:length(model.rowvalsB)
                fprintf(fid,['  JB->indexvals[' num2str(i-1) '] = ' num2str(model.rowvalsB(i)) ';\n']);
            end
            for i = 1:length(model.colptrsB)
                fprintf(fid,['  JB->indexptrs[' num2str(i-1) '] = ' num2str(model.colptrsB(i)) ';\n']);
            end
        case 'sxdot'
            if(~strcmp(model.wtype,'iw'))
                fprintf(fid,['memset(sxdot_tmp,0,sizeof(realtype)*' num2str(nx) ');\n']);
            end
        case 'sx0'
            fprintf(fid,['realtype t = udata->tstart;\n']);
        case 's2x0'
            fprintf(fid,['realtype t = udata->tstart;\n']);
            % nothing
        case 'sdx0'
            fprintf(fid,['realtype t = udata->tstart;\n']);
            % nothing
        case 'y'
            % nothing
        case 'dydp'
            % nothing
        case 'dydx'
            % nothing
        case 'z'
            % nothing
        case 'rz'
            % nothing
        case 'sz'
            % nothing
        case 'srz'
            % nothing
        case 'dzdp'
            % nothing
        case 'dzdx'
            % nothing
        case 'drzdp'
            % nothing
        case 'drzdx'
            % nothing
        case 'deltax'
            fprintf(fid,['memset(tdata->deltax,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'deltaxB'
            fprintf(fid,['memset(tdata->deltaxB,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'deltaqB'
            fprintf(fid,['memset(tdata->deltaqB,0,sizeof(realtype)*udata->nplist*model->nJ);\n']);
        case 'deltasx'
            fprintf(fid,['memset(tdata->deltasx,0,sizeof(realtype)*' num2str(nx) '*udata->nplist);\n']);
        case 'stau'
            fprintf(fid,['memset(tdata->stau,0,sizeof(realtype)*udata->nplist);\n']);
        case 'dxdotdp'
            fprintf(fid,'int ix;\n');
            fprintf(fid,['memset(tdata->dxdotdp,0,sizeof(realtype)*' num2str(nx) '*udata->nplist);\n']);
        case 'root'
            % nothing
        case 'sigma_y'
            fprintf(fid,['memset(tdata->sigmay,0,sizeof(realtype)*' num2str(model.ny) ');\n']);
        case 'dsigma_ydp'
            fprintf(fid,['memset(tdata->dsigmaydp,0,sizeof(realtype)*' num2str(model.ny) '*udata->nplist);\n']);
        case 'sigma_z'
            fprintf(fid,['memset(tdata->sigmaz,0,sizeof(realtype)*' num2str(model.nz) ');\n']);
        case 'dsigma_zdp'
            fprintf(fid,['memset(tdata->dsigmazdp,0,sizeof(realtype)*' num2str(model.nz) '*udata->nplist);\n']);
        case 'Jy'
            % nothing
        case 'dJydy'
            fprintf(fid,['memset(tdata->dJydy,0,sizeof(realtype)*model->ny*model->nytrue*model->nJ);\n']);
        case 'dJydsigma'
            fprintf(fid,['memset(tdata->dJydsigma,0,sizeof(realtype)*model->nytrue*model->ny*model->nJ);\n']);
        case 'Jz'
            % nothing
        case 'Jrz'
            % nothing
        case 'dJzdz'
            fprintf(fid,['memset(tdata->dJzdz,0,sizeof(realtype)*model->nz*model->nztrue*model->nJ);\n']);
        case 'dJrzdz'
            fprintf(fid,['memset(tdata->dJrzdz,0,sizeof(realtype)*model->nz*model->nztrue*model->nJ);\n']);
        case 'dJzdsigma'
            fprintf(fid,['memset(tdata->dJzdsigma,0,sizeof(realtype)*model->nztrue*model->nz*model->nJ);\n']);
        case 'dJrzdsigma'
            fprintf(fid,['memset(tdata->dJrzdsigma,0,sizeof(realtype)*model->nztrue*model->nz*model->nJ);\n']);
        case 'w'
            fprintf(fid,['memset(tdata->w,0,sizeof(realtype)*' num2str(model.nw) ');\n']);
        case 'dwdx'
            fprintf(fid,['memset(tdata->dwdx,0,sizeof(realtype)*' num2str(model.ndwdx) ');\n']);
        case 'dwdp'
            fprintf(fid,['memset(tdata->dwdp,0,sizeof(realtype)*' num2str(model.ndwdp) ');\n']);
        case 'M'
            fprintf(fid,['memset(tdata->M,0,sizeof(realtype)*' num2str(model.nx^2) ');\n']);
        case 'dfdx'
            fprintf(fid,['memset(tdata->dfdx,0,sizeof(realtype)*' num2str(model.nx^2) ');\n']);
        case 'dJdx'
            fprintf(fid,'int ix;\n');
            fprintf(fid,['memset(tdata->dJdx,0,sizeof(realtype)*' num2str(nx^3) ');\n']);
        case 'dJdp'
            fprintf(fid,'int ix;\n');
            fprintf(fid,['memset(tdata->dJdx,0,sizeof(realtype)*' num2str(nx^2) '*udata->nplist);\n']);
        case 'ddxdotdpdp'
            fprintf(fid,'int ix;\n');
            fprintf(fid,['memset(tdata->ddxdotdpdp,0,sizeof(realtype)*' num2str(nx) '*udata->nplist*udata->nplist);\n']);
        case 'ddJydydy'
            fprintf(fid,'memset(tdata->ddJydydy,0,sizeof(realtype)*model->nytrue*model->ny*model->ny);\n');
        case 'ddJydsigmady'
            fprintf(fid,'memset(tdata->ddJydsigmady,0,sizeof(realtype)*model->nytrue*model->ny*model->ny);\n');
        case 'ddJydsigmadsigma'
            fprintf(fid,'memset(tdata->ddJydsigmadsigma,0,sizeof(realtype)*model->nytrue*model->ny*model->ny);\n');
        case 'ddJy_s2sigma'
            fprintf(fid,'memset(tdata->ddJy_s2sigma,0,sizeof(realtype)*model->nytrue*udata->nplist*udata->nplist);\n');
        case 'ddydpdp'
            % nothing
        case 'ddydpdx'
            % nothing
        case 'ddydxdx'
            % nothing
        otherwise
            error(['unkown function: ' this.funstr])
    end
    
end