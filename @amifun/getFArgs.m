function this = getFArgs(this)
    % getFArgs populates the fargstr property with the argument string of
    % the respective f-function (if applicable). f-function are wrapped
    % implementations of functions which no longer have a model specific
    % name and have solver independent calls.
    %
    % Parameters:
    %  model: model definition object @type amimodel
    %
    % Return values:
    %  this: updated function definition object @type amifun
    %
    
    switch(this.funstr)
        case 'xdot'
            this.fargstr = '(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, void *user_data)';
        case 'J'
            this.fargstr = '(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xdot, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)';
        case 'root'
            this.fargstr = '(realtype t, N_Vector x, N_Vector dx, realtype *root, void *user_data)';
        otherwise
            this.fargstr = this.argstr;
    end
    
end