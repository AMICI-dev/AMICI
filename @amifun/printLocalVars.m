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
            fprintf(fid,'memset(qBdot_tmp,0,sizeof(realtype)*nplist*ng);\n');
        case 'x0'
            fprintf(fid,['memset(x0_tmp,0,sizeof(realtype)*' num2str(nx) ');\n']);
            fprintf(fid,['realtype t = tstart;\n']);
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
            fprintf(fid,['realtype t = tstart;\n']);
            % nothing
        case 'sdx0'
            fprintf(fid,['realtype t = tstart;\n']);
            % nothing
        case 'y'
            % nothing
        case 'sy'
            % nothing
        case 'dydp'
            % nothing
        case 'dydx'
            % nothing
        case 'z'
            % nothing
        case 'sz'
            % nothing
        case 'sz_tf'
            % nothing
        case 'dzdp'
            % nothing
        case 'dzdx'
            % nothing
        case 'deltax'
            fprintf(fid,['memset(deltax,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'deltaxB'
            fprintf(fid,['memset(deltaxB,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'deltaqB'
            fprintf(fid,['memset(deltaqB,0,sizeof(realtype)*nplist*ng);\n']);
        case 'deltasx'
            fprintf(fid,['memset(deltasx,0,sizeof(realtype)*' num2str(nx) '*nplist);\n']);
        case 'stau'
            fprintf(fid,['memset(stau,0,sizeof(realtype)*nplist);\n']);
        case 'dxdotdp'
            fprintf(fid,'int ix;\n');
            fprintf(fid,['memset(dxdotdp,0,sizeof(realtype)*' num2str(nx) '*nplist);\n']);
        case 'root'
            % nothing
        case 'sroot'
            % nothing
        case 's2root'
            % nothing
        case 'sigma_y'
            fprintf(fid,['memset(sigma_y,0,sizeof(realtype)*' num2str(model.ny) ');\n']);
        case 'dsigma_ydp'
            fprintf(fid,['memset(dsigma_ydp,0,sizeof(realtype)*' num2str(model.ny) '*nplist);\n']);
        case 'sigma_z'
            fprintf(fid,['memset(sigma_z,0,sizeof(realtype)*' num2str(model.nz) ');\n']);
        case 'dsigma_zdp'
            fprintf(fid,['memset(dsigma_zdp,0,sizeof(realtype)*' num2str(model.nz) '*nplist);\n']);
        case 'Jy'
            % nothing
        case 'dJydx'
            % nothing
        case 'dJydy'
            fprintf(fid,['memset(dJydy,0,sizeof(realtype)*nytrue*nytrue*ng);\n']);
        case 'dJydp'
            fprintf(fid,['memset(dJydp,0,sizeof(realtype)*nytrue*nplist*ng);\n']);
        case 'sJy'
            % nothing
        case 'Jz'
            % nothing
        case 'dJzdx'
            % nothing
        case 'dJzdp'
            % nothing
        case 'sJz'
            % nothing
        case 'w'
            fprintf(fid,['memset(w_tmp,0,sizeof(realtype)*' num2str(model.nw) ');\n']);
        case 'dwdx'
            fprintf(fid,['memset(dwdx_tmp,0,sizeof(realtype)*' num2str(model.ndwdx) ');\n']);
        case 'dwdp'
            fprintf(fid,['memset(dwdp_tmp,0,sizeof(realtype)*' num2str(model.ndwdp) ');\n']);
        case 'M'
            fprintf(fid,['memset(M_tmp,0,sizeof(realtype)*' num2str(model.nx^2) ');\n']);
        case 'dfdx'
            fprintf(fid,['memset(dfdx_tmp,0,sizeof(realtype)*' num2str(model.nx^2) ');\n']);
        otherwise
            error(['unkown function: ' this.funstr])
    end
    
end