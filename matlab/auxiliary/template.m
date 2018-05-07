classdef template < handle
    %TEMPLATE A class to replace strings in template files
    
    properties
        % strings in the template to be replaced
        templateStrings = {};
        % strings for replacement
        templateReplacements = {};
    end
    
    methods
        function replace(this, infile, outfile)
            % apply all provided template substitutions to infile and write
            % results to outfile
            fin = fopen(infile, 'r');
            fout = fopen(outfile, 'w');
            while(~feof(fin))
                s = fgetl(fin);
                s = this.replaceStr(s);
                fprintf(fout, '%s\n', s);
            end
            fclose(fin);
            fclose(fout);
        end
        
        function add(this, templateStr, replacementStr)
            % add a new template string and replacement
            nextIdx = numel(this.templateStrings);
            this.templateStrings{nextIdx + 1} = templateStr;
            this.templateReplacements{nextIdx + 1} = replacementStr;
        end
        
        function s = replaceStr(this, s)
            % apply all provided template substitutions to s 
            
            % do not use cellfun to guarantee order of replacements
            for n = 1:numel(this.templateStrings)
                s = strrep(s, this.templateStrings(n), this.templateReplacements(n));
                s = s{1};
            end
        end
    end
    
end

