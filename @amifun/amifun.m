%
% @file amifun
% @brief definition of amifun class
%
classdef amifun
    % the amifun class defines the prototype for all functions which later
    % on will be transformed into 
    
    properties ( GetAccess = 'public', SetAccess = 'public' )
        % symbolic definition struct @type symbolic
        sym;
        % flag which stores information for whether the symbolic defininition needs to be regenerated @type string
        generate;
        % name of the model @type string
        funstr;
        % name of the symbolic variable @type string
        svar;
        % name of the c variable @type string
        cvar;
        % argument string
        args;
    end
    
    methods
        function AF = amifun(funstr)
            AF.funstr = funstr;
            AF = AF.getDeps();
            AF = AF.getArgs();
            AF = AF.getCvar();
            AF = AF.getSvar();
        end

        writeCcode_sensi(this,model,fid)
        
        writeCcode(this,model,fid,ip,jp)
        
        gccode(this,fid)
        
        [ this ] = getDeps(this)
        
        [ this ] = getArgs(this)
        
        [ this ] = getCVar(this)
        
        [ this ] = getSVar(this)
        
        [ this ] = checkDeps(this)

    end
end

