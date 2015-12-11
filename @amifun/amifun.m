%
% @file amifun
% @brief definition of amifun class
%
classdef amifun
    % the amifun class defines the prototype for all functions which later
    % on will be transformed into C code
    
    properties ( GetAccess = 'public', SetAccess = 'public' )
        % symbolic definition struct @type symbolic
        sym;
        % short symbolic string which can be used for the reuse of precomputed values @type symbolic
        strsym;
        % flag which stores information for whether the symbolic defininition needs to be regenerated @type string
        generate;
        % name of the model @type char
        funstr;
        % name of the symbolic variable @type char
        svar;
        % name of the c variable @type char
        cvar;
        % argument string @type char
        argstr;
        % dependencies on other functions @type cell
        deps;
    end
    
    methods
        function AF = amifun(funstr,model)
            AF.funstr = funstr;
            AF = AF.getDeps(model);
            AF = AF.getArgs(model);
            AF = AF.getCVar();
            AF = AF.getSVar();
        end

        writeCcode_sensi(this,model,fid)
        
        writeCcode(this,model,fid,ip,jp)
        
        gccode(this,fid)
        
        [ this ] = getDeps(this,model)
        
        [ this ] = getArgs(this,model)
        
        [ this ] = getCVar(this)
        
        [ this ] = getSVar(this)
        
        [ this ] = checkDeps(this)

    end
end

