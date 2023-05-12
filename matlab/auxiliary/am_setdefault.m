function robj = am_setdefault(obj,dobj)
    % sets the undefined fields of the struct obj to the default values specified in the default struct dobj.
    %
    % Parameters:
    %  obj: struct which is supposed to be updated based on default struct @type struct
    %  dobj: obj which carries d fields @type struct
    %
    % Return values:
    %  robj: updated obj @type struct

    fieldlist = fieldnames(obj);
    for i = 1:length(fieldlist)
        dobj.(fieldlist{i}) = obj.(fieldlist{i});
    end
    robj = dobj;
