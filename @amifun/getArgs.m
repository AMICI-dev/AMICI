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
    
%     if(strcmp(model.wtype,'iw'))
        dxvec = ' N_Vector dx,';
        sdxvec = ' N_Vector sdx,';
        dxBvec = ' N_Vector dxB,';
        rtcj = ' realtype cj,';
%         s = '*';
%         intip = '';
        tmp3vec = ', N_Vector tmp3';
%     else
%         dxvec = '';
%         sdxvec = '';
%         dxBvec = '';
%         rtcj = '';
        s = '';
        intip = 'int ip, ';
%         tmp3vec = '';
%     end
    
    
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
            this.argstr = ['(realtype t, N_Vector x,' dxvec ' N_Vector xdot, N_Vector v, N_Vector Jv,' rtcj ' void *user_data, N_Vector tmp1, N_Vector tmp2)'];
        case 'JvB'
            this.argstr = ['(realtype t, N_Vector x,' dxvec ' N_Vector xB,' dxBvec ' N_Vector xBdot, N_Vector vB, N_Vector JvB,' rtcj ' void *user_data, N_Vector tmpB1, N_Vector tmpB2)'];
        case 'JBand'
            this.argstr = ['(long int N, long int mupper, long int mlower, realtype t,' rtcj ' N_Vector x,' dxvec ' N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)'];
        case 'J'
            this.argstr = ['(long int N, realtype t,' rtcj ' N_Vector x,' dxvec ' N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)'];
        case 'dJdx'
            this.argstr = '(realtype t, N_Vector x, N_Vector dx, void *user_data)';
        case 'dJdp'
            this.argstr = '(realtype t, N_Vector x, N_Vector dx, void *user_data)';
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
            this.argstr = '(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data)';
        case 's2x0'
            this.argstr = '(realtype *s2x0, N_Vector x, N_Vector dx, void *user_data)';
        case 'sdx0'
            this.argstr = '(N_Vector *sdx0, N_Vector x, N_Vector dx, void *user_data)';
        case 's2dx0'
            this.argstr = '(realtype *s2dx0, N_Vector x, N_Vector dx, void *user_data)';
        case 'root'
            this.argstr = ['(realtype t, N_Vector x,' dxvec ' realtype *root, void *user_data)'];
        case 'y'
            this.argstr = '(realtype t, int it, N_Vector x, void *user_data, ReturnData *rdata)';
        case 'z'
            this.argstr = '(realtype t, int ie, N_Vector x, TempData *tdata, ReturnData *rdata)';
        case 'rz'
            this.argstr = '(realtype t, int ie, N_Vector x, TempData *tdata, ReturnData *rdata)';
        case 'sz'
            this.argstr = '(realtype t, int ie, N_Vector x, N_Vector *sx, TempData *tdata, ReturnData *rdata)';
        case 'srz'
            this.argstr = '(realtype t, int ie, N_Vector x, N_Vector *sx, TempData *tdata, ReturnData *rdata)';
        case 'dydp'
            this.argstr = '(realtype t, int it, N_Vector x, TempData *tdata)';
        case 'dydx'
            this.argstr = '(realtype t, int it, N_Vector x, TempData *tdata)';
        case 'ddydpdp'
            this.argstr = '(realtype t, int it, N_Vector x, TempData *tdata)';
        case 'ddydpdx'
            this.argstr = '(realtype t, int it, N_Vector x, TempData *tdata)';
        case 'ddydxdx'
            this.argstr = '(realtype t, int it, N_Vector x, TempData *tdata)';
        case 'dzdp'
            this.argstr = '(realtype t, int ie, N_Vector x, TempData *tdata)';
        case 'dzdx'
            this.argstr = '(realtype t, int ie, N_Vector x, TempData *tdata)';
        case 'drzdp'
            this.argstr = '(realtype t, int ie, N_Vector x, TempData *tdata)';
        case 'drzdx'
            this.argstr = '(realtype t, int ie, N_Vector x, TempData *tdata)';
        case 'deltax'
            this.argstr = '(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, TempData *tdata)';
        case 'deltaxB'
            this.argstr = '(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector xdot, N_Vector xdot_old, TempData *tdata)';
        case 'deltaqB'
            this.argstr = '(realtype t, int ie, N_Vector x, N_Vector xB, N_Vector qBdot, N_Vector xdot, N_Vector xdot_old, TempData *tdata)';
        case 'deltasx'
            this.argstr = '(realtype t, int ie, N_Vector x, N_Vector xdot, N_Vector xdot_old, N_Vector *sx, TempData *tdata)';
        case 'dxdotdp'
            this.argstr = '(realtype t, N_Vector x, N_Vector dx, void *user_data)';
        case 'ddxdotdpdp'
            this.argstr = '(realtype t, N_Vector x, N_Vector dx, void *user_data)';
        case 'sigma_y'
            this.argstr = '(realtype t, TempData *tdata)';
        case 'dsigma_ydp'
            this.argstr = '(realtype t, TempData *tdata)';
        case 'sigma_z'
            this.argstr = '(realtype t, int ie, TempData *tdata)';
        case 'dsigma_zdp'
            this.argstr = '(realtype t, int ie, TempData *tdata)';
        case 'stau'
            this.argstr = '(realtype t, int ie, N_Vector x, N_Vector *sx, TempData *tdata)';
        case 'sroot'
            this.argstr = '(realtype t, int ie, N_Vector x, N_Vector *sx, TempData *tdata, ReturnData *rdata)';
        case 'Jy'
            this.argstr = '(realtype t, int it, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata)';
        case 'dJydy'
            this.argstr = '(realtype t, int it, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata)';
        case 'dJydsigma'
            this.argstr = '(realtype t, int it, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata)';
        case 'ddJydydy'
            this.argstr = '(realtype t, int it, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata)';
        case 'ddJydsigmady'
            this.argstr = '(realtype t, int it, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata)';
        case 'ddJydsigmadsigma'
            this.argstr = '(realtype t, int it, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata)';
        case 'ddJy_s2sigma'
            this.argstr = '(realtype t, int it, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata)';
        case 'Jz'
            this.argstr = '(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata)';
        case 'Jrz'
            this.argstr = '(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata)';
        case 'dJzdz'
            this.argstr = '(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata)';
        case 'dJzdsigma'
            this.argstr = '(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata)';
        case 'dJrzdz'
            this.argstr = '(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata)';
        case 'dJrzdsigma'
            this.argstr = '(realtype t, int ie, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata)';

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