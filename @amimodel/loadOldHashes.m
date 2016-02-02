function HTable = loadOldHashes(this)
    % loadOldHashes loads information from a previous compilation of the model.
    %
    % Return values:
    %  HTable: struct with hashes of symbolic definition from the previous
    %  compilation @type struct
    
    [wrap_path,~,~]=fileparts(which('amiwrap.m'));
    try
        load(fullfile(wrap_path,'models',this.modelname,['hashes.mat']))
        try
            if(this.compver>compver)
                this.recompile = true;
                error('catchmeifyoucan');
            end
        catch
            disp('Code changes made recompilation of model neventcessary!')
            this.recompile = true;
            error('catchmeifyoucan');
        end
        try
            this.nxtrue = nxtrue;
            this.nytrue = nytrue;
            this.nx = nx;
            this.ny = ny;
            this.np = np;
            this.nk = nk;
            this.nevent = nevent;
            this.nz = nz;
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
    DHTable.root = '';

    DHTable.p = '';
    DHTable.x0 = '';
    DHTable.y = '';
    DHTable.u = '';
    DHTable.sigma_y = '';
    DHTable.sigma_z = '';
    DHTable.z = '';
    DHTable.trigger = '';
    DHTable.bolus = '';
    if(strcmp(this.wtype,'iw'))
        DHTable.xdot = '';
        DHTable.dx0 = '';
        DHTable.M = '';
    else
        DHTable.xdot = '';
    end
    HTable = am_setdefault(HTable,DHTable);
end