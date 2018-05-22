function HTable = loadOldHashes(this)
    % loadOldHashes loads information from a previous compilation of the model.
    %
    % Return values:
    %  HTable: struct with hashes of symbolic definition from the previous
    %  compilation @type struct
    
    [wrap_path,~,~]=fileparts(fileparts(which('amiwrap.m')));
    try
        load(fullfile(wrap_path,'models',this.modelname,['hashes.mat']))
        try
            this.nw = nw;
            this.ndwdp = ndwdp;
            this.ndwdx = ndwdx;
            this.z2event = z2event;
            this.nnz = nnonzeros;
            this.id = id;
            this.ubw = ubw;
            this.lbw = lbw;
            this.colptrs = colptrs;
            this.rowvals = rowvals;
            this.sparseidx = sparseidx;
            this.colptrsB = colptrsB;
            this.rowvalsB = rowvalsB;
            this.sparseidxB = sparseidxB;
        catch err
        end
        
    catch
        HTable = struct();
    end
    DHTable.x = '';
    DHTable.k = '';

    DHTable.p = '';
    DHTable.x0 = '';
    DHTable.y = '';
    DHTable.sigma_y = '';
    DHTable.sigma_z = '';
    DHTable.z = '';
    DHTable.trigger = '';
    DHTable.deltax = '';
    DHTable.xdot = '';
    if(strcmp(this.wtype,'iw'))
        DHTable.dx0 = '';
        DHTable.M = '';
    end
    DHTable.Jy = '';
    DHTable.Jz = '';
    
    DHTable.generateC = '';
    DHTable.makeSyms = '';
    DHTable.makeEvents = '';
    DHTable.gccode = '';
    DHTable.getArgs = '';
    DHTable.getCVar = '';
    DHTable.getFArgs = '';
    DHTable.getSensiFlag = '';
    DHTable.getSyms = '';
    DHTable.printLocalVars = '';
    DHTable.writeCcode = '';
    DHTable.writeCcode_sensi = '';
    DHTable.tdata = '';

    
    HTable = am_setdefault(HTable,DHTable);
end
