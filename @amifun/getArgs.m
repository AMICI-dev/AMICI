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
            this.argstr = ['(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w)'];
        case 'xBdot'
            this.argstr = ['(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx)'];
        case 'qBdot'
            this.argstr = ['(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdp)'];
        case 'x0'
            this.argstr = '(realtype *x0, const realtype t, const realtype *p, const realtype *k)';
        case 'dx0'
            this.argstr = '(N_Vector x0, N_Vector dx0, void *user_data)';
        case 'Jv'
            this.argstr = ['(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *v, const realtype *w, const realtype *dwdx)'];
        case 'JvB'
            this.argstr = ['(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *xB, const realtype *dx, const realtype *dxB, const realtype *vB, const realtype *w, const realtype *dwdx)'];
        case 'J'
            this.argstr = ['(realtype *J, const realtype t, const realtype *x, const double *p, const double *k, const realtype *h, const realtype *w, const realtype *dwdx)'];
        case 'JDiag'
            this.argstr = ['(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx)'];
        case 'JSparse'
            this.argstr = ['(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx)'];
        case 'JB'
            this.argstr = ['(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx)'];
        case 'JSparseB'
            this.argstr = ['(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx)'];
        case 'sxdot'
            this.argstr = ['(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *sx, const realtype *w, const realtype *dwdx, const realtype *J, const realtype *dxdotdp)'];
        case 'sx0'
            this.argstr = '(realtype *sx0, const realtype t,const realtype *x0, const realtype *p, const realtype *k, const int ip)';
        case 'root'
            this.argstr = ['(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h)'];
        case 'y'
            this.argstr = '(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h)';
        case 'z'
            this.argstr = '(double *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h)';
        case 'rz'
            this.argstr = '(double *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h)';
        case 'sz'
            this.argstr = '(double *sz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip)';
        case 'srz'
            this.argstr = '(double *srz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip)';
        case 'dydp'
            this.argstr = '(double *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip)';
        case 'dydx'
            this.argstr = '(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h)';
        case 'dzdp'
            this.argstr = '(double *dzdp, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip)';
        case 'dzdx'
            this.argstr = '(double *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h)';
        case 'drzdp'
            this.argstr = '(double *drzdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip)';
        case 'drzdx'
            this.argstr = '(double *drzdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h)';
        case 'deltax'
            this.argstr = '(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old)';
        case 'deltaxB'
            this.argstr = '(double *deltaxB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB)';
        case 'deltaqB'
            this.argstr = '(double *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB, const realtype *qBdot)';
        case 'deltasx'
            this.argstr = '(double *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau)';
        case 'dxdotdp'
            this.argstr = ['(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp)'];
        case 'sigma_y'
            this.argstr = '(double *sigmay, const realtype t, const realtype *p, const realtype *k)';
        case 'dsigma_ydp'
            this.argstr = '(double *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip)';
        case 'sigma_z'
            this.argstr = '(double *sigmaz, const realtype t, const realtype *p, const realtype *k)';
        case 'dsigma_zdp'
            this.argstr = '(double *dsigmazdp, const realtype t, const realtype *p, const realtype *k, const int ip)';
        case 'stau'
            this.argstr = '(double *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip, const int ie)';
        case 'Jy'
            this.argstr = '(double *nllh, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my)';
        case 'dJydy'
            this.argstr = '(double *dJydy, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my)';
        case 'dJydsigma'
            this.argstr = '(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my)';
        case 'Jz'
            this.argstr = '(double *nllh, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz)';
        case 'Jrz'
            this.argstr = '(double *nllh, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz)';
        case 'dJzdz'
            this.argstr = '(double *dJzdz, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz)';
        case 'dJzdsigma'
            this.argstr = '(double *dJzdsigma, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz)';
        case 'dJrzdz'
            this.argstr = '(double *dJrzdz, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz)';
        case 'dJrzdsigma'
            this.argstr = '(double *dJrzdsigma, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz)';
        case 'w'
            this.argstr = '(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h)';
        case 'dwdp'
            this.argstr = '(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w)';
        case 'dwdx'
            this.argstr = '(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w)';
        case 'M'
            this.argstr = '(realtype *M, const realtype t, const realtype *x, const realtype *p, const realtype *k)';
        otherwise
            %nothing
    end
    
end