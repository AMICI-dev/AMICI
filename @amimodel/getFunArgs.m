function argstr = getFunArgs(this,funstr)
    % getfunargs returns a string which contains all input arguments for the function fun.
    %
    % Parameters:
    % funstr: function name @type string
    %
    % Return values:
    % str: string containing function arguments @type string
    %
    
    if(strcmp(this.wtype,'iw'))
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
    
    
    switch(funstr)
        case 'xdot'
            argstr = ['(realtype t, N_Vector x,' dxvec ' N_Vector xdot, void *user_data)'];
        case 'xBdot'
            argstr = ['(realtype t, N_Vector x,' dxvec ' N_Vector xB,' dxBvec ' N_Vector xBdot, void *user_data)'];
        case 'qBdot'
            argstr = ['(realtype t, N_Vector x,' dxvec ' N_Vector xB,' dxBvec ' N_Vector qBdot, void *user_data)'];
        case 'x0'
            argstr = '(N_Vector x0, void *user_data)';
        case 'dx0'
            argstr = '(N_Vector x0, N_Vector dx0, void *user_data)';
        case 'Jv'
            if(strcmp(this.wtype,'iw'))
                argstr = ['(realtype t, N_Vector x,' dxvec ' N_Vector xdot, N_Vector v, N_Vector Jv,' rtcj ' void *user_data, N_Vector tmp1, N_Vector tmp2)'];
            else
                argstr = '(N_Vector v, N_Vector Jv, realtype t, N_Vector x, N_Vector xdot, void *user_data, N_Vector tmp)';
            end
        case 'JvB'
            if(strcmp(this.wtype,'iw'))
                argstr = ['(realtype t, N_Vector x,' dxvec ' N_Vector xB,' dxBvec ' N_Vector xBdot, N_Vector vB, N_Vector JvB,' rtcj ' void *user_data, N_Vector tmpB1, N_Vector tmpB2)'];
            else
                argstr = '(N_Vector vB, N_Vector JvB, realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data, N_Vector tmpB)';
            end
        case 'JBand'
            argstr = ['(long int N, long int mupper, long int mlower, realtype t,' rtcj ' N_Vector x,' dxvec ' N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)'];
        case 'J'
            argstr = ['(long int N, realtype t,' rtcj ' N_Vector x,' dxvec ' N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)'];
        case 'JSparse'
            argstr = ['(realtype t,' rtcj ' N_Vector x,' dxvec ' N_Vector xdot, SlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)'];
        case 'JBandB'
            argstr = ['(long int NeqBdot, long int mupper, long int mlower, realtype t,' rtcj ' N_Vector x,' dxvec ' N_Vector xB,' dxBvec ' N_Vector xdotB, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)'];
        case 'JB'
            argstr = ['(long int NeqBdot, realtype t,' rtcj ' N_Vector x,' dxvec ' N_Vector xB,' dxBvec ' N_Vector xdotB, DlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)'];
        case 'JSparseB'
            argstr = ['(realtype t,' rtcj ' N_Vector x,' dxvec ' N_Vector xB,' dxBvec ' N_Vector xdotB, SlsMat JB, void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)'];
        case 'sxdot'
            argstr = ['(int Ns, realtype t, N_Vector x,' dxvec ' N_Vector xdot,' intip ' N_Vector ' s 'sx,' sdxvec ' N_Vector ' s 'sxdot, void *user_data, N_Vector tmp1, N_Vector tmp2' tmp3vec ')'];
        case 'sx0'
            argstr = ['(int ip, N_Vector sx0,' xvec dxvec ' void *user_data)'];
        case 'sdx0'
            argstr = ['(int ip, N_Vector sdx0,' xvec dxvec ' void *user_data)'];
        case 'y'
            argstr = '(realtype t, int it, realtype *y, realtype *x, void *user_data)';
        case 'sy'
            argstr = '(realtype t, int it, realtype *sy, realtype *x, realtype *sx, void *user_data)';
        case 'dydp'
            argstr = '(realtype t, int it, realtype *dydp, realtype *y, realtype *x, void *user_data)';
        case 'dydx'
            argstr = '(realtype t, int it, realtype *dydx, realtype *y, realtype *x, void *user_data)';
        case 'dtdx'
            argstr = '(realtype t, realtype *dtdx, realtype *x, void *user_data)';
        case 'drvaldx'
            argstr = '(realtype t, realtype *drvaldx, realtype *x, void *user_data)';
        case 'dtdp'
            argstr = '(realtype t, realtype *dtdp, realtype *x, void *user_data)';
        case 'drvaldp'
            argstr = '(realtype t, realtype *drvaldp, realtype *x, void *user_data)';
        case 'root'
            argstr = ['(realtype t, N_Vector x,' dxvec ' realtype *root, void *user_data)'];
        case 'sroot'
            argstr = '(realtype t, int nroots, realtype *sroot, N_Vector x, N_Vector *sx, void *user_data)';
        case 's2root'
            argstr = '(realtype t, int nroots, realtype *s2root, N_Vector x, N_Vector *sx, void *user_data)';
        case 'rootval'
            argstr = '(realtype t, realtype *rootval, N_Vector x, void *user_data)';
        case 'srootval'
            argstr = '(realtype t, int nroots, realtype *srootval, N_Vector x, N_Vector *sx, void *user_data)';
        case 's2rootval'
            argstr = '(realtype t, int nroots, realtype *s2rootval, N_Vector x, N_Vector *sx, void *user_data)';
        case 'deltadisc'
            argstr = '(realtype t, int idisc, N_Vector x, void *user_data)';
        case 'bdeltadisc'
            argstr = '(realtype t, int idisc, N_Vector x, N_Vector xB, void *user_data)';
        case 'ideltadisc'
            argstr = '(realtype t, int idisc, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data)';
        case 'sdeltadisc'
            argstr = '(realtype t, int idisc, N_Vector x, N_Vector *sx, void *user_data)';
        case 'dxdotdp'
            argstr = '(realtype t, N_Vector x, realtype *dxdotdp, void *user_data)';
        case 'sigma_y'
            argstr = '(realtype t, realtype *sigma_y, void *user_data)';
        case 'dsigma_ydp'
            argstr = '(realtype t, realtype *dsigma_ydp, void *user_data)';
        case 'sigma_t'
            argstr = '(realtype t, realtype *sigma_t, void *user_data)';
        case 'dsigma_tdp'
            argstr = '(realtype t, realtype *dsigma_tdp, void *user_data)';
        otherwise
            error(['unkown function: ' funstr])
    end
    
end