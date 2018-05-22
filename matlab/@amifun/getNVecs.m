function this = getNVecs(this)
    % getfunargs populates the nvecs property with the names of the
    % N_Vector elements which are required in the execution of the function 
    % (if applicable). the information is directly extracted from the
    % argument string
    %
    % Parameters:
    %
    % Return values:
    %  this: updated function definition object @type amifun
    %
    
    vecs = {'x,','dx,','sx,','*sx,','sdx,','xB,','dxB,',...
        '*sxdot,','sxdot,','xdot,','xBdot,','qBdot,',...
        'x0,','dx0,','*sx0,','*sdx0,',...
        'v,','vB,','JDiag,','Jv,','JvB,',...
        'xdot_old,'};
    
    this.nvecs = {};
    for iv = 1:length(vecs)
        if strfind(this.argstr,['N_Vector ' vecs{iv}])
            this.nvecs = [this.nvecs,vecs{iv}(1:(end-1))];
        end
    end
    
end