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
        % short symbolic string which can be used for the reuse of old values @type symbolic
        strsym_old@sym;
        % name of the model @type char
        funstr@char;
        % name of the c variable @type char
        cvar@char;
        % argument string (solver specific) @type char
        argstr@char;
        % argument string (solver unspecific) @type char
        fargstr@char;
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
            AF = AF.getFArgs();
            AF = AF.getNVecs();
            AF = AF.getCVar();
            AF = AF.getSensiFlag();
        end

        printLocalVars(this,model,fid)
        
        writeCcode_sensi(this,model,fid)
        
        writeCcode(this,model,fid)
        
        gccode(this,model,fid)
        
        [ this ] = getDeps(this,model)
        
        [ this ] = getArgs(this,model)
        
        [ this ] = getNVecs(this)
        
        [ this ] = getCVar(this)
        
        [ this ] = getSensiFlag(this)
        
        [ this ] = checkDeps(this)
    end
end

