function [ ] = structToHDF5Attribute( hdffile, path, mystruct )
%STRUCTTOHDF5ATTRIBUTE Write structs to HDF5 file
%   ...

try
    try
        fid = H5F.open(hdffile, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
    catch
        fid = H5F.create(hdffile, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
    end
    lcpl = H5P.create('H5P_LINK_CREATE');
    H5P.set_create_intermediate_group(lcpl,1);
    gid = H5G.create(fid, path, lcpl, 'H5P_DEFAULT', 'H5P_DEFAULT');
    H5G.close(gid);
    H5F.close(fid);
catch me
    if ~strcmp(me.identifier,'MATLAB:imagesci:h5create:datasetAlreadyExists')
        rethrow(me)
    end
end

fields = fieldnames(mystruct);
for i = 1:numel(fields)
    f = fields{i};
    d = mystruct.(f);
    if(isstruct(d))
        structToHDF5Attribute( hdffile, [path '/' f], d )
    else
        if isa(d, 'logical')
            d = int32(d);
        end
        if(numel(d)>0)
            if(numel(d)>1 && sum(size(d)>1)>1)
                h5create(hdffile, [path '/' f], size(d));
                h5write(hdffile, [path '/' f], d);
            else
                h5writeatt(hdffile, path, f, d);
            end
        end
    end
end
end