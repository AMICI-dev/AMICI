function this = getArgs(this,model)
    % getFArgs populates the fargstr property with the argument string of
    % the respective model function (if applicable). model functions are not
    % wrapped versions of functions which have a model specific name and
    % for which the call is solver specific.
    %
    % Parameters:
    %  model: model definition object @type amimodel
    %
    % Return values:
    %  this: updated function definition object @type amifun
    %
    
    if(strcmp(model.wtype,'iw'))
        xvec = ' N_Vector x,';
        dxvec = ' N_Vector dx,';
        sdxvec = ' N_Vector *sdx,';
        dxBvec = ' N_Vector dxB,';
        rtcj = ' realtype cj,';
        s = '*';
        intip = '';
        tmp3vec = ', N_Vector tmp3';
    else
        xvec = '';
        dxvec = '';
        sdxvec = '';
        dxBvec = '';
        rtcj = '';
        s = '';
        intip = 'int ip, ';
        tmp3vec = '';
    end
    
    
    switch(this.funstr)
        case 'xdot'
            this.argstr = ['(realtype t, N_Vector x,' dxvec ' N_Vector xdot, void *user_data)'];
            this.fargstr = '(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, void *user_data)';
        case 'xBdot'
            this.argstr = ['(realtype t, N_Vector x,' dxvec ' N_Vector xB,' dxBvec ' N_Vector xBdot, void *user_data)'];
        case 'qBdot'
            this.argstr = ['(realtype t, N_Vector x,' dxvec ' N_Vector xB,' dxBvec ' N_Vector qBdot, void *user_data)'];
        case 'x0'
            this.argstr = '(N_Vector x0, void *user_data)';
        case 'dx0'
            this.argstr = '(N_Vector x0, N_Vector dx0, void *user_data)';
        case 'Jv'
            if(strcmp(model.wtype,'iw'))
                this.argstr = ['(realtype t, N_Vector x,' dxvec ' N_Vector xdot, N_Vector v, N_Vector Jv,' rtcj ' void *user_data, N_Vector tmp1, N_Vector tmp2)'];
            else
                this.argstr = '(N_Vector v, N_Vector Jv, realtype t, N_Vector x, N_Vector xdot, void *user_data, N_Vector tmp)';
            end
        case 'JvB'
            if(strcmp(model.wtype,'iw'))
                this.argstr = ['(realtype t, N_Vector x,' dxvec ' N_Vector xB,' dxBvec ' N_Vector xBdot, N_Vector vB, N_Vector JvB,' rtcj ' void *user_data, N_Vector tmpB1, N_Vector tmpB2)'];
            else
                this.argstr = '(N_Vector vB, N_Vector JvB, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data, N_Vector tmpB)';
            end
        case 'JBand'
            this.argstr = ['(long int N, long int mupper, long int mlower, realtype t,' rtcj ' N_Vector x,' dxvec ' N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)'];
        case 'J'
            this.argstr = ['(long int N, realtype t,' rtcj ' N_Vector x,' dxvec ' N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)'];
            this.fargstr = '(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)';
        case 'JSparse'
            this.argstr = ['(realtype t,' rtcj ' N_Vector x,' dxvec ' N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)'];
        case 'JBandB'
            this.argstr = ['(long int NeqBdot, long int mupper, long int mlower, realtype t,' rtcj ' N_Vector x,' dxvec ' N_Vector xB,' dxBvec ' N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)'];
        case 'JB'
            this.argstr = ['(long int NeqBdot, realtype t,' rtcj ' N_Vector x,' dxvec ' N_Vector xB,' dxBvec ' N_Vector xBdot, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)'];
        case 'JSparseB'
            this.argstr = ['(realtype t,' rtcj ' N_Vector x,' dxvec ' N_Vector xB,' dxBvec ' N_Vector xBdot, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)'];
        case 'sxdot'
            this.argstr = ['(int Ns, realtype t, N_Vector x,' dxvec ' N_Vector xdot,' intip ' N_Vector ' s 'sx,' sdxvec ' N_Vector ' s 'sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2' tmp3vec ')'];
        case 'sx0'
            this.argstr = '(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data)';
        case 'sdx0'
            this.argstr = '(N_Vector *sdx0, N_Vector x, N_Vector dx, void *user_data)';
        case 'y'
            this.argstr = '(realtype t, int it, realtype *y, N_Vector x, void *user_data)';
        case 'sy'
            this.argstr = '(realtype t, int it, realtype *sy, realtype *dydx, realtype *dydp, N_Vector *sx, void *user_data)';
        case 'z'
            this.argstr = '(realtype t, int ie, int *nroots, realtype *z, N_Vector x, void *user_data)';
        case 'sz'
            this.argstr = '(realtype t, int ie, int *nroots, realtype *sz, N_Vector x, N_Vector *sx, void *user_data)';
        case 'sz_tf'
            this.argstr = '(realtype t, int ie, int *nroots, realtype *sz, N_Vector x, N_Vector *sx, void *user_data)';
        case 'dydp'
            this.argstr = '(realtype t, int it, realtype *dydp, N_Vector x, void *user_data)';
        case 'dydx'
            this.argstr = '(realtype t, int it, realtype *dydx, N_Vector x, void *user_data)';
        case 'dzdp'
            this.argstr = '(realtype t, int ie, realtype *dzdp, N_Vector x, void *user_data)';
        case 'dzdx'
            this.argstr = '(realtype t, int ie, realtype *dzdx, N_Vector x, void *user_data)';
        case 'deltax'
            this.argstr = '(realtype t, int ie, realtype *deltax, N_Vector x, N_Vector xdot, N_Vector xdot_old, void *user_data)';
        case 'deltaxB'
            this.argstr = '(realtype t, int ie, realtype *deltaxB, N_Vector x, N_Vector xB, N_Vector xdot, N_Vector xdot_old, void *user_data)';
        case 'deltaqB'
            this.argstr = '(realtype t, int ie, realtype *deltaqB, N_Vector x, N_Vector xB, N_Vector qBdot, N_Vector xdot, N_Vector xdot_old, void *user_data)';
        case 'deltasx'
            this.argstr = '(realtype t, int ie, realtype *deltasx, N_Vector x, N_Vector xdot, N_Vector xdot_old, N_Vector *sx, void *user_data)';
        case 'dxdotdp'
            this.argstr = ['(realtype t, realtype *dxdotdp, N_Vector x, N_Vector dx, void *user_data)'];
        case 'sigma_y'
            this.argstr = '(realtype t, realtype *sigma_y, void *user_data)';
        case 'dsigma_ydp'
            this.argstr = '(realtype t, realtype *dsigma_ydp, void *user_data)';
        case 'sigma_z'
            this.argstr = '(realtype t, int ie, realtype *sigma_z, void *user_data)';
        case 'dsigma_zdp'
            this.argstr = '(realtype t, int ie, realtype *dsigma_zdp, void *user_data)';
        case 'stau'
            this.argstr = '(realtype t, int ie, realtype *stau, N_Vector x, N_Vector *sx, void *user_data)';
        case 'root'
            this.argstr = ['(realtype t, N_Vector x,' dxvec ' realtype *root, void *user_data)'];
        case 'sroot'
            this.argstr = '(realtype t, int ie, int *nroots, realtype *sroot, N_Vector x, N_Vector *sx, void *user_data)';
        case 's2root'
            this.argstr = '(realtype t, int ie, int *nroots, realtype *s2root, N_Vector x, N_Vector *sx, void *user_data)';
        case 'Jy'
            this.argstr = '(realtype t, int it, realtype *Jy, realtype *y, N_Vector x, realtype *my, realtype *sigma_y, void *user_data)';
        case 'dJydx'
            this.argstr = '(realtype t, int it, realtype *dJydx, realtype *y, N_Vector x, realtype *dydx, realtype *my, realtype *sigma_y, void *user_data)';
        case 'dJydp'
            this.argstr = '(realtype t, int it, realtype *dJydp, realtype *y, N_Vector x, realtype *dydp, realtype *my, realtype *sigma_y, realtype *dsigma_ydp, void *user_data)';
        case 'sJy'
            this.argstr = '(realtype t, int it, realtype *sJy, realtype *dJydy, realtype *dJydp, realtype *sy, void *user_data)';
        case 'Jz'
            this.argstr = '(realtype t, int ie, realtype *Jz, realtype *z, N_Vector x, realtype *mz, realtype *sigma_z, void *user_data, void *temp_data)';
        case 'dJzdx'
            this.argstr = '(realtype t, int ie, realtype *dJzdx, realtype *z, N_Vector x, realtype *dzdx, realtype *mz, realtype *sigma_z, void *user_data, void *temp_data)';
        case 'dJzdp'
            this.argstr = '(realtype t, int ie, realtype *dJzdp, realtype *z, N_Vector x, realtype *dzdp, realtype *mz, realtype *sigma_z, realtype *dsigma_zdp, void *user_data, void *temp_data)';
        case 'sJz'
            this.argstr = '(realtype t, int ie, realtype *sJz, realtype *z, N_Vector x, realtype *dzdp, realtype *sz, realtype *mz, realtype *sigma_z, realtype *dsigma_zdp, void *user_data, void *temp_data)';
        case 'w'
            this.argstr = ['(realtype t, N_Vector x, N_Vector dx, void *user_data)'];
        case 'dwdp'
            this.argstr = ['(realtype t, N_Vector x, N_Vector dx, void *user_data)'];
        case 'dwdx'
            this.argstr = ['(realtype t, N_Vector x, N_Vector dx, void *user_data)'];
        case 'M'
            this.argstr = ['(realtype t, N_Vector x, N_Vector dx, void *user_data)'];
        case 'dfdx'
            this.argstr = ['(realtype t, N_Vector x, N_Vector dx, void *user_data)'];
        otherwise
            %nothing
    end
    
end