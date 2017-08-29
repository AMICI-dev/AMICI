function getFun(this,HTable,funstr)
    % getFun generates symbolic expressions for the requested function.
    %
    % Parameters:
    %  HTable: struct with hashes of symbolic definition from the previous
    %  compilation @type struct
    %  funstr: function for which symbolic expressions should be computed @type string
    %
    % Return values:
    %  void
    
    [wrap_path,~,~]=fileparts(which('amiwrap.m'));
    
    fun = amifun(funstr,this);
    

    if(~isfield(this.fun,funstr)) % check whether we already computed the respective fun
        if(~all(strcmp(fun.deps,funstr))) % prevent infinite loops
            % check whether the corresponding c file exists, if not we have to
            % force recompilation by passing an empty HTable
            if(~exist(fullfile(wrap_path,'models',this.modelname,[this.modelname '_' funstr '.cpp']),'file') ...
                    || and(~exist(fullfile(wrap_path,'models',this.modelname,[funstr  '_' this.modelname '.m'])),ismember(funstr,this.mfuns)))
                cflag = this.checkDeps([],fun.deps);
            else
                if(this.recompile)
                    cflag = this.checkDeps([],fun.deps);
                else
                    cflag = this.checkDeps(HTable,fun.deps);
                end
            end
        else
            if(~isempty(fun.deps))
                if(~isfield(this.strsym,[fun.deps{1} 's']))
                    cflag = 1;
                else
                    cflag = 0;
                end
            else
                cflag = 1;
            end
        end
    else
        cflag = 0;
    end
        
    
    if(cflag)
        fun = amifun(funstr,this);
        [fun,this] = fun.getSyms(this);
        this.fun(1).(funstr) = fun;
    end
    
end