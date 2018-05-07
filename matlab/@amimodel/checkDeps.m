function cflag = checkDeps(this,HTable,deps)
    % checkDeps checks the dependencies of functions and populates sym
    % fields if necessary
    %
    % Parameters:
    %  HTable: struct with reference hashes of functions in its fields
    %  @type struct
    %  deps: cell array with containing a list of dependencies @type cell
    %
    % Return values:
    %  cflag: boolean indicating whether any of the dependencies have 
    %  changed with respect to the hashes stored in HTable @type
    %  bool
    
    if(~isempty(HTable))
        cflags = zeros(length(deps),1);
        for id = 1:length(deps)
            if(ismember(deps{id},{'w','dwdx','dwdp'}))
                cflags(id) = 0;
            else
                if(~isfield(this.HTable,deps{id}))
                    % check subdependencies
                    fun = amifun(deps{id},this);
                    fun = fun.getDeps(this);
                    cflagdep = this.checkDeps(HTable,fun.deps);
                    cflags(id) = cflagdep;
                else
                    if(~isfield(HTable,deps{id})) %newly added field in HTable, need to recompile
                        cflags(id) = 1;
                    else
                        cflags(id) = ~strcmp(this.HTable.(deps{id}),HTable.(deps{id}));
                    end
                end
            end
        end
        cflag = any(cflags);
    else
        cflag = 1;
    end
    if(cflag)
        % if the function needs to be computed, we need to make sure that
        % all symbolic definitions are generated. this is forced by passing
        % an empty HTable to getFun
        for id = 1:length(deps)
            if(~isfield(this.fun,deps{id}))
                this.getFun([],deps{id});
            end
        end
    end
end