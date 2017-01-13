function g = mfun(f,varargin)
% adapted from matlabFunction, follows the same syntax, but supports
% sparsity
%matlabFunction Generate a MATLAB file or anonymous function from a sym
%   G = matlabFunction(F) generates a MATLAB anonymous function from sym object
%   F. The free variables of F become the inputs for the resulting function
%   handle G. For example, if x is the free variable of F then G(z) computes in
%   MATLAB what subs(F,x,z) computes in the symbolic engine. The function handle
%   G can be used by functions in MATLAB like FZERO or by functions in other
%   toolboxes. The order of the inputs of G matches the order returned by
%   symvar(F).
%
%   G = matlabFunction(...,PARAM1,VALUE1,...) uses the specified parameter/value
%   pairs to customize the generated function. The parameters modify the
%   declaration that appears in the generated code. The template for the
%   generated code is
%      function [OUT] = NAME(IN1,IN2)
%   The parameter names can be any of the following
%
%     'file': The value, FILE, must be a string with a valid MATLAB function
%             name. If FILE does not end in '.m' it will be appended. If the
%             file already exists it will be overwritten.  A function handle to
%             the function will be returned in G. If the file name parameter is
%             empty an anonymous function is generated. The NAME in the function
%             template is the file name in FILE.
%
%     'vars': The value, IN, must be either
%                1) a cell array of M strings or sym arrays, or
%                2) a vector of M symbolic variables.
%             IN specifies the input variable names and their order IN1,
%             IN2,... INM in the function template. If INj is a sym array then
%             the name used in the function template is 'inj'.  The variables
%             listed in IN must be a superset of the free variables in all the
%             Fk. The default value for IN is the union of the free variables in
%             all the Fk.
%
%   Note: not all MuPAD expressions can be converted to a MATLAB function.
%   For example piecewise expressions and sets will not be converted.


narginchk(1,inf);

% process inputs
funs = sym(f);
args = varargin(1:end);
opts = getOptions(args);

% compute non-trivial defaults and verify inputs
funvars = getFunVars(funs);
vars = checkVars(funvars,opts);
varnames = format(opts.varnames);

% generate anonymous function or file
if ~isempty(opts.file)
    file = normalize(opts.file,'.m');
%     body = renameFileInputs(vars,inputs,funvars);
    clear(char(file)); % necessary to avoid problems in for-loops
    g = writeMATLAB(funs,file,varnames);
    tmp = exist(char(file),'file'); %#ok
                                    % necessary to make the function available
                                    % on the PATH right away
end

% compute the default value for 'vars'.
% returns the sorted union of the symvars of the funs
function funvars = getFunVars(f)
vars = symvar(f);
vars = unique(vars);
funparams = argnames(f);
funvars = unique([funparams, vars], 'stable');

% get the 'vars' value and check it for errors
function v = checkVars(funvars,opts)
if isempty(opts.vars)
    v = funvars;
else
    v = opts.vars;
end
[v,vexpanded] = var2cell(v);
if ~isempty(vexpanded)
    if ~iscellstr(vexpanded)
        error(message('symbolic:sym:matlabFunction:InvalidVars'));
    elseif ~all(cellfun(@(x)isvarname(x),vexpanded))
        error(message('symbolic:sym:matlabFunction:InvalidVarName'));
    end
end
checkVarsSubset(vexpanded,funvars);

% check that the funvars are a subset of the expanded vars.
function checkVarsSubset(vexpanded,funvars)
vars = var2cell(funvars);
missing = cellfun(@(x)~any(strcmp(char(x),vexpanded)),vars);
if any(missing)
    misvars = vars(missing);
    varnames = format(misvars);
    if length(misvars) > 1
        error(message('symbolic:sym:matlabFunction:FreeVariables', varnames));
    end
    error(message('symbolic:sym:matlabFunction:FreeVariable', varnames));
end

% convert a string or sym array into a 1-by-N cell array of strings
% also optionally return the cellstr of expanded sym array vars joined together
function [v,vexpand] = var2cell(v)
if isa(v,'sym')
    v = privsubsref(v,':');
    v = num2cell(v.');
    v = cellfun(@(x)char(x),v,'UniformOutput',false);
elseif ischar(v)
    v = {v};
end
if nargout > 1
    vexpand = cellfun(@var2cell,v,'UniformOutput',false);
    vexpand = [vexpand{:}];
end

function inputs = getInputs(v)
inputs = cell(size(v));
for k = 1:length(v)
        inputs{k} = sprintf('in%d',k);
end

% file = normalizeFile(file,ext) append extension ext if file doesn't
% have one already.
function file = normalize(file,ext)
[~,~,x] = fileparts(file);
if isempty(x)
    file = [file ext];
end

% Horzcat sym array or cell array v with commas
function varnames = format(v)
if isempty(v)
    varnames = '';
else
    if ~iscell(v)
        v = num2cell(v);
    end
    varnames = cellfun(@(x)[char(x) ','],v,'UniformOutput',false);
    varnames = [varnames{:}];
    varnames(end) = [];
end

% Generate MATLAB. f is the expr to generate. file is file name.
% varnames is the formatted input variables
% outputs is the cell array of output names
% mapping is string with input to variable mapping
function g = writeMATLAB(f,file,varnames)
[fid,msg] = fopen(file,'wt');
if fid == -1
    error(message('symbolic:sym:matlabFunction:FileError', file, msg));
end
tmp = onCleanup(@()fclose(fid));
% [f,tvalues,tnames] = optimize(f);
[~,fname] = fileparts(file);
writeHeader(fid,fname,varnames);
writeOutput(fid,f);
g = str2func(fname);

% write out the function declaration and help
function writeHeader(fid,fname,varnames)
symver = ver('symbolic');
if ~isempty(varnames)
    varnames = ['(' varnames ')'];
end
symver = symver(1);
fprintf(fid,'function out = %s%s\n',fname,varnames);
fprintf(fid,'%%%s\n',upper(fname));
fprintf(fid,'%%    out = %s%s\n\n',upper(fname),upper(varnames));
fprintf(fid,'%%    This function was generated by the Symbolic Math Toolbox version %s.\n',symver.Version);
fprintf(fid,'%%    %s\n\n',datestr(now));


% write the assignments to the output variables
function writeOutput(fid,expr)
idx = find(expr); % find nonzero entries
fprintf(fid,['out = sparse([],[],[],' num2str(size(expr,1)) ',' num2str(size(expr,2)) ',' num2str(length(idx)) ');\n']); % initialise with zero values

for k = 1:length(idx)
    str = char(expr(idx(k)));
    for var = {'p','x','k'}
        t = regexp(str,[var{1} '_([0-9]*)'],'tokens');
        ts = cellfun(@(x) [var{1} '_' x{1} ],t,'UniformOutput',false);
        td = cellfun(@(x) [var{1} '(' num2str(str2double(x)+1) ')'],t,'UniformOutput',false);
        ts = unique(ts);
        td = unique(td);
        [~,idx_ts] = sort(cellfun(@(x) length(x),ts),'descend');
        ts = ts(idx_ts);
        td = ts(idx_ts);
        for it = 1:length(ts)
            str = strrep(str,ts{it},td{it});
        end
    end
    fprintf(fid,['out(' num2str(idx(k)) ') = ' str ';\n']);
end

% validator for variable parameter
function t = isVars(x)
t = iscell(x) || (ischar(x)&&size(x,1)==1) || isa(x,'sym');

% validator for file parameter
function t = isFunc(x)
[~,file] = fileparts(x);
if length(file)>namelengthmax
    error('filename is longer than what matlab can handle.')
end
t = isempty(file) || isvarname(file);

% parse inputs and return option structure output
function opts = getOptions(args)
ip = inputParser;
verMAT = ver('MATLAB');
if(str2double(verMAT.Version)>=8.2)
    ip.addParameter('vars',{},@isVars);
    ip.addParameter('file','',@isFunc);
    ip.addParameter('outputs',{},@iscellstr);
    ip.addParameter('varnames',{},@iscellstr);
else
    ip.addParamValue('vars',{},@isVars);
    ip.addParamValue('file','',@isFunc);
    ip.addParamValue('outputs',{},@iscellstr);
    ip.addParamValue('varnames',{},@iscellstr);
end
ip.parse(args{:});
opts = ip.Results;
