function this = getArgs(this,model)
    % getfunargs returns a string which contains all input arguments for the function fun.
    %
    % Parameters:
    % funstr: function name @type string
    %
    % Return values:
    % str: string containing function arguments @type string
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
            this.argstr = ['(long int NeqBdot, long int mupper, long int mlower, realtype t,' rtcj ' N_Vector x,' dxvec ' N_Vector xB,' dxBvec ' N_Vector xdotB, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)'];
        case 'JB'
            this.argstr = ['(long int NeqBdot, realtype t,' rtcj ' N_Vector x,' dxvec ' N_Vector xB,' dxBvec ' N_Vector xdotB, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)'];
        case 'JSparseB'
            this.argstr = ['(realtype t,' rtcj ' N_Vector x,' dxvec ' N_Vector xB,' dxBvec ' N_Vector xdotB, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)'];
        case 'sxdot'
            this.argstr = ['(int Ns, realtype t, N_Vector x,' dxvec ' N_Vector xdot,' intip ' N_Vector ' s 'sx,' sdxvec ' N_Vector ' s 'sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2' tmp3vec ')'];
        case 'sx0'
            this.argstr = ['(int ip, N_Vector sx0,' xvec dxvec ' void *user_data)'];
        case 'sdx0'
            this.argstr = ['(int ip, N_Vector sdx0,' xvec dxvec ' void *user_data)'];
        case 'y'
            this.argstr = '(realtype t, int it, realtype *y, realtype *x, void *user_data)';
        case 'sy'
            this.argstr = '(realtype t, int it, realtype *sy, realtype *x, realtype *sx, void *user_data)';
        case 'dydp'
            this.argstr = '(realtype t, int it, realtype *dydp, realtype *y, realtype *x, void *user_data)';
        case 'dydx'
            this.argstr = '(realtype t, int it, realtype *dydx, realtype *y, realtype *x, void *user_data)';
        case 'dtdx'
            this.argstr = '(realtype t, realtype *dtdx, realtype *x, void *user_data)';
        case 'drvaldx'
            this.argstr = '(realtype t, realtype *drvaldx, realtype *x, void *user_data)';
        case 'dtdp'
            this.argstr = '(realtype t, realtype *dtdp, realtype *x, void *user_data)';
        case 'drvaldp'
            this.argstr = '(realtype t, realtype *drvaldp, realtype *x, void *user_data)';
        case 'root'
            this.argstr = ['(realtype t, N_Vector x,' dxvec ' realtype *root, void *user_data)'];
        case 'sroot'
            this.argstr = '(realtype t, int nroots, realtype *sroot, N_Vector x, N_Vector *sx, void *user_data)';
        case 's2root'
            this.argstr = '(realtype t, int nroots, realtype *s2root, N_Vector x, N_Vector *sx, void *user_data)';
        case 'rootval'
            this.argstr = '(realtype t, realtype *rootval, N_Vector x, void *user_data)';
        case 'srootval'
            this.argstr = '(realtype t, int nroots, realtype *srootval, N_Vector x, N_Vector *sx, void *user_data)';
        case 's2rootval'
            this.argstr = '(realtype t, int nroots, realtype *s2rootval, N_Vector x, N_Vector *sx, void *user_data)';
        case 'deltadisc'
            this.argstr = '(realtype t, int idisc, N_Vector x, void *user_data)';
        case 'bdeltadisc'
            this.argstr = '(realtype t, int idisc, N_Vector x, N_Vector xB, void *user_data)';
        case 'ideltadisc'
            this.argstr = '(realtype t, int idisc, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data)';
        case 'sdeltadisc'
            this.argstr = '(realtype t, int idisc, N_Vector x, N_Vector *sx, void *user_data)';
        case 'dxdotdp'
            this.argstr = '(realtype t, N_Vector x, realtype *dxdotdp, void *user_data)';
        case 'sigma_y'
            this.argstr = '(realtype t, realtype *sigma_y, void *user_data)';
        case 'dsigma_ydp'
            this.argstr = '(realtype t, realtype *dsigma_ydp, void *user_data)';
        case 'sigma_t'
            this.argstr = '(realtype t, realtype *sigma_t, void *user_data)';
        case 'dsigma_tdp'
            this.argstr = '(realtype t, realtype *dsigma_tdp, void *user_data)';
        case 'p'
            this.argstr = '';
        case 'x'
            this.argstr = '';
        case 'k'
            this.argstr = '';
        otherwise
            error(['unkown function: ' this.funstr])
    end
    
end