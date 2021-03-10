function generateRebuildM(this)
% generateRebuildM generates a Matlab script for recompilation of this
% model
%
% Return values:
%  void


filename = fullfile(this.wrap_path,'models',this.modelname,['rebuild_',this.modelname,'.m']);
% would require struct to string conversion, skipping for now
amimodelo2 = '[]';

fid = fopen(filename, 'w');

fprintf(fid, ['function ' ['rebuild_', this.modelname] '()\n']);
fprintf(fid, ['modelName = ''' this.modelname ''';\n']);
fprintf(fid, 'amimodel.compileAndLinkModel(modelName, '''', [], [], [], []);\n');
fprintf(fid, ['amimodel.generateMatlabWrapper(' num2str(this.nx) ', ' num2str(this.ny) ...
    ', ' num2str(this.np) ', ' num2str(this.nk) ', ' num2str(this.nz) ', ' num2str(this.o2flag) ', ' ...
    amimodelo2 ', [''simulate_'' modelName ''.m''], ''' ...
    this.modelname ''', ''' this.param ''', ' num2str(this.forward) ', ' num2str(this.adjoint) ');\n']);
fprintf(fid, 'end\n');

fclose(fid);
end

