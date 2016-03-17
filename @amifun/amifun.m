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
        % short symbolic string which can be used for the reuse of old values @type symbolic
        strsym_old;
        % name of the model @type char
        funstr;
        % name of the c variable @type char
        cvar;
        % argument string (solver specific) @type char
        argstr;
        % argument string (solver unspecific) @type char
        fargstr;
        % dependencies on other functions @type cell
        deps;
        % nvec dependencies @type cell
        nvecs;
        % indicates whether the function is a sensitivity or derivative
        % with respect to parameters @type logical
        sensiflag;
    end
    
    methods
        function AF = amifun(funstr,model)
            % constructor of the amifun class. this function initializes the function object based on the provided
            % function name funstr and model definition object model
            %
            % Parameters:
            %  funstr: name of the function @type string
            %  model: model definition object @type amimodel
            % 
            % Return values:
            %  AM: model definition object
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
        
        [ this ] = getFArgs(this,model)
        
        [ this ] = getNVecs(this)
        
        [ this ] = getCVar(this)
        
        [ this ] = getSyms(this,model)
        
        [ this ] = getSensiFlag(this)
    end
end

