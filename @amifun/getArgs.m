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
        dxvec = ' N_Vector dx,';
        sdxvec = ' N_Vector sdx,';
        dxBvec = ' N_Vector dxB,';
        rtcj = ' realtype cj,';
        tmp3vec = ', N_Vector tmp3';
    else
        dxvec = '';
        sdxvec = '';
        dxBvec = '';
        rtcj = '';
        tmp3vec = '';
    end
    
    s = '';
    intip = 'int ip, ';
    
    switch(this.funstr)
        case 'xdot'
            this.argstr = ['(realtype t, N_Vector x,' dxvec ' N_Vector xdot, void *user_data)'];
        case 'xBdot'
            this.argstr = ['(realtype t, N_Vector x,' dxvec ' N_Vector xB,' dxBvec ' N_Vector xBdot, void *user_data)'];
        case 'qBdot'
            this.argstr = ['(realtype t, N_Vector x,' dxvec ' N_Vector xB,' dxBvec ' N_Vector qBdot, void *user_data)'];
        case 'x0'
            this.argstr = '(realtype *x0, const realtype t, const realtype *p, const realtype *k)';
        case 'dx0'
            this.argstr = '(N_Vector x0, N_Vector dx0, void *user_data)';
        case 'Jv'
            this.argstr = ['(realtype t, N_Vector x,' dxvec ' N_Vector xdot, N_Vector v, N_Vector Jv,' rtcj ' void *user_data, N_Vector tmp1, N_Vector tmp2)'];
        case 'JvB'
            this.argstr = ['(realtype t, N_Vector x,' dxvec ' N_Vector xB,' dxBvec ' N_Vector xBdot, N_Vector vB, N_Vector JvB,' rtcj ' void *user_data, N_Vector tmpB1, N_Vector tmpB2)'];
        case 'JBand'
            this.argstr = ['(long int N, long int mupper, long int mlower, realtype t,' rtcj ' N_Vector x,' dxvec ' N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)'];
        case 'J'
            this.argstr = ['(long int N, realtype t,' rtcj ' N_Vector x,' dxvec ' N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)'];
        case 'JDiag'
            this.argstr = ['(realtype t, N_Vector JDiag,' rtcj ' N_Vector x,' dxvec ' void *user_data)'];
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
            this.argstr = '(realtype *sx0, const realtype t,const realtype *x0, const realtype *p, const realtype *k, const int ip)';
        case 'root'
            this.argstr = ['(realtype t, N_Vector x,' dxvec ' realtype *root, void *user_data)'];
        case 'y'
            this.argstr = '(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k)';
        case 'z'
            this.argstr = '(double *z, const realtype t, const realtype *x, const realtype *p, const realtype *k)';
        case 'rz'
            this.argstr = '(double *rz, const realtype t, const realtype *x, const realtype *p, const realtype *k)';
        case 'sz'
            this.argstr = '(double *sz, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *sx, const int ip)';
        case 'srz'
            this.argstr = '(double *srz, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *sx, const int ip)';
        case 'dydp'
            this.argstr = '(double *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip)';
        case 'dydx'
            this.argstr = '(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k)';
        case 'dzdp'
            this.argstr = '(double *dzdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip)';
        case 'dzdx'
            this.argstr = '(double *dzdx, const realtype t, const realtype *x, const realtype *p, const realtype *k)';
        case 'drzdp'
            this.argstr = '(double *drzdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip)';
        case 'drzdx'
            this.argstr = '(double *drzdx, const realtype t, const realtype *x, const realtype *p, const realtype *k)';
        case 'deltax'
            this.argstr = '(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ie, const realtype *xdot, const realtype *xdot_old)';
        case 'deltaxB'
            this.argstr = '(double *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB)';
        case 'deltaqB'
            this.argstr = '(double *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB, const realtype *qBdot)';
        case 'deltasx'
            this.argstr = '(double *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau)';
        case 'dxdotdp'
            this.argstr = ['(realtype t, N_Vector x,' dxvec ', void *user_data)'];
        case 'sigma_y'
            this.argstr = '(double *sigmay, const realtype t, const realtype *p, const realtype *k)';
        case 'dsigma_ydp'
            this.argstr = '(double *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip)';
        case 'sigma_z'
            this.argstr = '(double *sigmaz, const realtype t, const realtype *p, const realtype *k)';
        case 'dsigma_zdp'
            this.argstr = '(double *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip)';
        case 'stau'
            this.argstr = '(double *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *sx, const int ip, const int ie)';
        case 'Jy'
            this.argstr = '(double *nllh, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my)';
        case 'dJydy'
            this.argstr = '(double *dJydy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my)';
        case 'dJydsigma'
            this.argstr = '(double *dJydsigma, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my)';
        case 'Jz'
            this.argstr = '(double *nllh, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz)';
        case 'Jrz'
            this.argstr = '(double *nllh, const realtype *p, const realtype *k, const double *z, const double *sigmaz)';
        case 'dJzdz'
            this.argstr = '(double *dJzdz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz)';
        case 'dJzdsigma'
            this.argstr = '(double *dJzdsigma, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz)';
        case 'dJrzdz'
            this.argstr = '(double *dJrzdz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz)';
        case 'dJrzdsigma'
            this.argstr = '(double *dJrzdsigma, const realtype *p, const realtype *k, const double *rz, const double *sigmaz)';
        case 'w'
            this.argstr = '(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k)';
        case 'dwdp'
            this.argstr = '(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *w)';
        case 'dwdx'
            this.argstr = '(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *w)';
        case 'M'
            this.argstr = '(realtype *M, const realtype t, const realtype *x, const realtype *p, const realtype *k)';
        otherwise
            %nothing
    end
    
end