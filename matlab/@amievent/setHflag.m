function this = setHflag(this,hflag)
    % setHflag sets the hflag property.
    %
    % Parameters:
    %  hflag: value for the hflag property, type double
    %
    % Return values:
    %  this: updated event definition object @type amievent

    try
        if(all(size(this.bolus) == size(hflag)))
            if(isa(hflag,'double'))
                this.hflag = hflag~=0;
            elseif(islogical(hflag))
                this.hflag = hflag;
            else
                error('provided hflag is not a double/logical value!');
            end
        else
            error('provided hflag does not match the bolus dimension!');
        end
    catch
        error('provided hflag does not match the bolus dimension!');
    end
