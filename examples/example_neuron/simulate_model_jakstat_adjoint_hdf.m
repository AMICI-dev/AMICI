    
    %% Write data that is passed to AMICI to HDF5
    global amiHDFfile;
    global amiHDFprefix;
    if ~isempty(amiHDFfile)
        amiHDFOptionPath = [amiHDFprefix '/options'];
        structToHDF5Attribute(amiHDFfile, amiHDFOptionPath, options_ami);
        h5writeatt(amiHDFfile, amiHDFOptionPath, 'ts', tout);
        h5writeatt(amiHDFfile, amiHDFOptionPath, 'nt', numel(tout));
        h5writeatt(amiHDFfile, amiHDFOptionPath, 'theta', theta);
        h5writeatt(amiHDFfile, amiHDFOptionPath, 'kappa', kappa);
        if(~isempty(data))
            structToHDF5Attribute(amiHDFfile, [amiHDFprefix '/data'], data);
        end
        structToHDF5Attribute(amiHDFfile, [amiHDFprefix '/results'], sol);
    end
end    
    %% Write data that is passed to AMICI to HDF5
    global amiHDFfile;
    global amiHDFprefix;
    if ~isempty(amiHDFfile)
        amiHDFOptionPath = [amiHDFprefix '/options'];
        structToHDF5Attribute(amiHDFfile, amiHDFOptionPath, options_ami);
        h5writeatt(amiHDFfile, amiHDFOptionPath, 'ts', tout);
        h5writeatt(amiHDFfile, amiHDFOptionPath, 'nt', numel(tout));
        h5writeatt(amiHDFfile, amiHDFOptionPath, 'theta', theta);
        h5writeatt(amiHDFfile, amiHDFOptionPath, 'kappa', kappa);
        if(~isempty(data))
            structToHDF5Attribute(amiHDFfile, [amiHDFprefix '/data'], data);
        end
        structToHDF5Attribute(amiHDFfile, [amiHDFprefix '/results'], sol);
    end
end    
    %% Write data that is passed to AMICI to HDF5
    global amiHDFfile;
    global amiHDFprefix;
    if ~isempty(amiHDFfile)
        amiHDFOptionPath = [amiHDFprefix '/options'];
        structToHDF5Attribute(amiHDFfile, amiHDFOptionPath, options_ami);
        h5writeatt(amiHDFfile, amiHDFOptionPath, 'ts', tout);
        h5writeatt(amiHDFfile, amiHDFOptionPath, 'nt', numel(tout));
        h5writeatt(amiHDFfile, amiHDFOptionPath, 'theta', theta);
        h5writeatt(amiHDFfile, amiHDFOptionPath, 'kappa', kappa);
        if(~isempty(data))
            structToHDF5Attribute(amiHDFfile, [amiHDFprefix '/data'], data);
        end
        structToHDF5Attribute(amiHDFfile, [amiHDFprefix '/results'], sol);
    end
end