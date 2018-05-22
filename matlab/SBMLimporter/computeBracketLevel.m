function [ brl ] = computeBracketLevel( cstr )
    % Compute the bracket level for the input string cstr. The bracket
    % level is computed for every char in cstr and indicates how many
    % brackets have been opened up to this point. The bracket level is
    % useful to parse the arguments of functions in cstr. For this purpose
    % functions will have the same bracket level as the opening bracket of
    % the corresponding function call.
    %
    % Parameters:
    %  cstr: input string @type *char
    %
    % Return values:
    %  brl: bracket levels @type *int
    
    % compute bracket levels add one for each (, (before) remove 1 for each
    % ) (after)
    open = (cstr == '(');
    close = (cstr == ')');
    close = [0,close(1:end-1)];
    brl = cumsum(open) - cumsum(close);
    % take care of functions
    fun_startidx = regexp(cstr,'([\w_]+\()','start');
    fun_endidx = regexp(cstr,'([\w_]+\()','end');
    for ifun = 1:length(fun_startidx)
        brl(fun_startidx(ifun):(fun_endidx(ifun)-1)) = brl(fun_endidx(ifun));
    end
    
end

