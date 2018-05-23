function csym = betterSym(str)
    matVer = ver('MATLAB');
    if(str2double(matVer.Version)>=9.4)
        csym = str2sym(str);
    else
        csym = sym(str);
    end
end