function [ commit_hash,branch,url ] = getCommitHash( wrap_path )
    % extracts the commit hash from the git HEAD, when available. returns
    % empty array when no git directory is available.
    %
    % Parameters:
    %  wrap_path: path to AMICI folder @type char
    %
    % Return values:
    %  commit_hash: extracted hash value @type char
    %  branch: branch of the repository @type char
    %  url: employed remote origin @type char
    
    try
        fid = fopen(fullfile(wrap_path,'..','.git','FETCH_HEAD'));
        str = fgetl(fid);
        fclose(fid);
        if(str~=-1)
            t_hash = regexp(str,'^([\w]*)','tokens');
            commit_hash = t_hash{1}{1};
            t_branch = regexp(str,'branch ''([\w]*)''','tokens');
            branch = t_branch{1}{1};
            idx_url = strfind(str,'https://github.com');
            url = str(idx_url:end);
        else
            fid = fopen(fullfile(wrap_path,'.git','ORIG_HEAD'));
            commit_hash = ['dev_' fgetl(fid)];
            fclose(fid);
            
            fid = fopen(fullfile(wrap_path,'.git','HEAD'));
            branch = strrep(fgetl(fid),'ref: refs/heads/','');
            fclose(fid);
            
            url = 'local';
        end
    catch
        commit_hash = [];
        branch = 'unknown branch';
        url = 'unknown repository';
    end
end

