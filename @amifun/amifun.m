%
% @file amifun
% @brief definition of amifun class
%
classdef amifun
    % the amifun class defines the prototype for all functions which later
    % on will be transformed into C code
    
    properties ( GetAccess = 'public', SetAccess = 'public' )
        % symbolic definition struct @type symbolic
        sym@sym;
        % short symbolic string which can be used for the reuse of precomputed values @type symbolic
        strsym@sym;
        % name of the model @type char
        funstr@char;
        % name of the c variable @type char
        cvar@char;
        % argument string @type char
        argstr@char;
        % dependencies on other functions @type cell
        deps@cell;
        % nvec dependencies
        nvecs@cell;
        % indicates whether the function is a sensitivity or derivative
        % with respect to parameters
        sensiflag@logical;
    end
    
    methods
        function AF = amifun(funstr,model)
            AF.funstr = funstr;
            AF = AF.getDeps(model);
            AF = AF.getArgs(model);
            AF = AF.getNVecs(model);
            AF = AF.getCVar();
            AF = AF.getSensiFlag();
        end

        printLocalVars(this,model,fid)
        
        writeCcode_sensi(this,model,fid)
        
        writeCcode(this,model,fid)
        
        gccode(this,fid)
        
        [ this ] = getDeps(this,model)
        
        [ this ] = getArgs(this,model)
        
        [ this ] = getNVecs(this,model)
        
        [ this ] = getCVar(this)
        
        [ this ] = getSensiFlag(this)
        
        [ this ] = checkDeps(this)
    end
end

