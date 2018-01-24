function splinefun = am_spline(varargin)
    n = nargin;
    str= '';
    if (round(n/2) - n/2 < 0.1)
        error('Input arguments of am_spline must have the following form: t, t1, p2, ..., tn, pn, intss, dudt');
    end
    for i = 1 : n-1
        str = strcat(strcat(str, char(varargin{i})), ',');
    end
    str = strcat('(',strcat(strcat(str, char(varargin{n})), ')'));
    str = strrep(str, ' ', '');
    str = regexprep(str,'\,([0-9]*)\,','\,$1\.0\,'); % Black magic of sym
    splinefun = sym(strcat('spline', str));
end
