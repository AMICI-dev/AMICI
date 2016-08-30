function splinefun = am_spline_pos(varargin)
    n = nargin;
    str= '';
    if (round(n/2) - n/2 > 0.1)
        error('Input arguments of am_spline must have the following form: t, #ofNodes, t1, p1, ..., tn, pn, intss, dudt');
    end
    for i = 1 : n
        if (i < n)
            switch class(varargin{i})
                case 'sym'
                    str = strcat(strcat(str, char(varargin{i})), ',');
                case 'double'
                    str = strcat(strcat(str, num2str(varargin{i})), ',');
                otherwise
                    error(['Input argument ' num2str(i) ' of the splinefunction seems to be neither a symbolic nor a double. Please check, if it is correctly defined.']);
            end
        else
            switch class(varargin{i})
                case 'sym'
                    str = strcat(str, char(varargin{i}));
                case 'double'
                    str = strcat(str, num2str(varargin{i}));
                otherwise
                    error(['Input argument ' num2str(i) ' of the splinefunction seems to be neither a symbolic nor a double. Please check, if it is correctly defined.']);
            end
        end

    end
    str = strcat('(',strcat(strcat(str, char(varargin{n})), ')'));
    str = strrep(str, ' ', '');
    str = regexprep(str,'\,([0-9]*)\,','\,$1\.0\,'); % Black magic of sym
    str = regexprep(str,'\,([0-9]*)\,','\,$1\.0\,'); % Black magic of sym
    str = regexprep(str,'\,([0-9]*)\))','\,$1\.0\)'); % Black magic of sym
    
    splinefun = sym(strcat('am_spline_pos', str));
end
