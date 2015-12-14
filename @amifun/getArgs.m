function this = getArgs(this,model)
    % getfunargs populates the argstr property with the argument string of
    % the respective function (if applicable)
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
            this.argstr = ['(N_Vector sx0,' xvec dxvec ' void *user_data)'];
        case 'sdx0'
            this.argstr = ['(N_Vector sdx0,' xvec dxvec ' void *user_data)'];
        case 'y'
            this.argstr = '(realtype t, int it, realtype *y, N_Vector x, void *user_data)';
        case 'sy'
            this.argstr = '(realtype t, int it, realtype *sy, N_Vector x, N_Vector *sx, void *user_data)';
        case 'z'
            this.argstr = '(realtype t, int ie, realtype *z, N_Vector x, void *user_data)';
        case 'sz'
            this.argstr = '(realtype t, int ie, realtype *sz, N_Vector x, N_Vector *sx, void *user_data)';
        case 'dydp'
            this.argstr = '(realtype t, int it, realtype *dydp, N_Vector x, void *user_data)';
        case 'dydx'
            this.argstr = '(realtype t, int it, realtype *dydx, N_Vector x, void *user_data)';
        case 'dzdp'
            this.argstr = '(realtype t, int ie, realtype *dzdp, N_Vector x, void *user_data)';
        case 'dzdx'
            this.argstr = '(realtype t, int ie, realtype *dzdx, N_Vector x, void *user_data)';
        case 'deltax'
            this.argstr = '(realtype t, realtype *deltax, N_Vector x, void *user_data)';
        case 'deltaxB'
            this.argstr = '(realtype t, realtype *deltaxB, N_Vector x, N_Vector xB, void *user_data)';
        case 'deltaqB'
            this.argstr = '(realtype t, realtype *deltaqB, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data)';
        case 'deltasx'
            this.argstr = '(realtype t, realtype *deltasx, N_Vector x, N_Vector *sx, void *user_data)';
        case 'dxdotdp'
            this.argstr = '(realtype t, realtype *dxdotdp, N_Vector x, void *user_data)';
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
        otherwise
            % do nothing
    end
    
end