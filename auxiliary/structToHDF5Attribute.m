function [ ] = structToHDF5Attribute( hdffile, path, mystruct )
%STRUCTTOHDF5ATTRIBUTE Write structs to HDF5 file
%   ...

try
    h5create(hdffile, path, 1);
catch me
    if ~strcmp(me.identifier,'MATLAB:imagesci:h5create:datasetAlreadyExists')
        rethrow(me)
    end
end

fields = fieldnames(mystruct);
for i = 1:numel(fields)
    f = fields{i};
    d = mystruct.(f);
    if isa(d, 'logical')
        d = int32(d);
    end
    h5writeatt(hdffile, path, f, d);
end
end

