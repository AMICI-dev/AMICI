function fun = spline_pos5(t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, intss, dudt)
args = [ t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, intss, dudt];
str = [];
for j = 1:length(args)-1
    str = [ str char(args(j)) ', ' ];
end
str = [str char(args(end))];
str = ['(' strrep(str,' ','') ')'];
% we have to replace all numbers by number.0 in oder for sym to be able to
% generate a correctly differentiable function (this is matlab magic)
str = regexprep(str,'\,([0-9]*)\,','\,$1\.0\,');
str = regexprep(str,'\,([0-9]*)\)','\,$1\.0\)');
fun = sym(['spline_pos5' str '']);
end

function fun = spline5(t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, intss, dudt)
args = [ t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, intss, dudt];
str = [];
for j = 1:length(args)-1
str = [ str char(args(j)) ', ' ];
end
str = [str char(args(end))];
str = ['(' strrep(str,' ','') ')'];
% we have to replace all numbers by number.0 in oder for sym to be able to
% generate a correctly differentiable function (this is matlab magic)
str = regexprep(str,'\,([0-9]*)\,','\,$1\.0\,');
str = regexprep(str,'\,([0-9]*)\)','\,$1\.0\)');
fun = sym(['spline5' str '']);
end

function fun = spline_pos4(t, t1, p1, t2, p2, t3, p3, t4, p4, intss, dudt)
args = [ t, t1, p1, t2, p2, t3, p3, t4, p4, intss, dudt];
str = [];
for j = 1:length(args)-1
str = [ str char(args(j)) ', ' ];
end
str = [str char(args(end))];
str = ['(' strrep(str,' ','') ')'];
% we have to replace all numbers by number.0 in oder for sym to be able to
% generate a correctly differentiable function (this is matlab magic)
str = regexprep(str,'\,([0-9]*)\,','\,$1\.0\,');
str = regexprep(str,'\,([0-9]*)\)','\,$1\.0\)');
fun = sym(['spline_pos4' str '']);
end

function fun = spline4(t, t1, p1, t2, p2, t3, p3, t4, p4, intss, dudt)
args = [ t, t1, p1, t2, p2, t3, p3, t4, p4, intss, dudt];
str = [];
for j = 1:length(args)-1
str = [ str char(args(j)) ', ' ];
end
str = [str char(args(end))];
str = ['(' strrep(str,' ','') ')'];
% we have to replace all numbers by number.0 in oder for sym to be able to
% generate a correctly differentiable function (this is matlab magic)
str = regexprep(str,'\,([0-9]*)\,','\,$1\.0\,');
str = regexprep(str,'\,([0-9]*)\)','\,$1\.0\)');
fun = sym(['spline4' str '']);
end


function fun = spline_pos3(t, t1, p1, t2, p2, t3, p3, intss, dudt)
args = [ t, t1, p1, t2, p2, t3, p3, intss, dudt];
str = [];
for j = 1:length(args)-1
str = [ str char(args(j)) ', ' ];
end
str = [str char(args(end))];
str = ['(' strrep(str,' ','') ')'];
% we have to replace all numbers by number.0 in oder for sym to be able to
% generate a correctly differentiable function (this is matlab magic)
str = regexprep(str,'\,([0-9]*)\,','\,$1\.0\,');
str = regexprep(str,'\,([0-9]*)\)','\,$1\.0\)');
fun = sym(['spline_pos3' str '']);
end

function fun = spline3(t, t1, p1, t2, p2, t3, p3, intss, dudt)
args = [ t, t1, p1, t2, p2, t3, p3, intss, dudt];
str = [];
for j = 1:length(args)-1
str = [ str char(args(j)) ', ' ];
end
str = [str char(args(end))];
str = ['(' strrep(str,' ','') ')'];
% we have to replace all numbers by number.0 in oder for sym to be able to
% generate a correctly differentiable function (this is matlab magic)
str = regexprep(str,'\,([0-9]*)\,','\,$1\.0\,');
str = regexprep(str,'\,([0-9]*)\)','\,$1\.0\)');
fun = sym(['spline3' str '']);
end

function fun = spline_pos10(t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7, t8, p8, t9, p9, t10, p10, intss, dudt)
args = [ t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7, t8, p8, t9, p9, t10, p10, intss, dudt];
str = [];
for j = 1:length(args)-1
str = [ str char(args(j)) ', ' ];
end
str = [str char(args(end))];
str = ['(' strrep(str,' ','') ')'];
% we have to replace all numbers by number.0 in oder for sym to be able to
% generate a correctly differentiable function (this is matlab magic)
str = regexprep(str,'\,([0-9]*)\,','\,$1\.0\,');
str = regexprep(str,'\,([0-9]*)\)','\,$1\.0\)');
fun = sym(['spline_pos10' str '']);
end

function fun = spline10(t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7, t8, p8, t9, p9, t10, p10, intss, dudt)
args = [ t, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7, t8, p8, t9, p9, t10, p10, intss, dudt];
str = [];
for j = 1:length(args)-1
str = [ str char(args(j)) ', ' ];
end
str = [str char(args(end))];
str = ['(' strrep(str,' ','') ')'];
% we have to replace all numbers by number.0 in oder for sym to be able to
% generate a correctly differentiable function (this is matlab magic)
str = regexprep(str,'\,([0-9]*)\,','\,$1\.0\,');
str = regexprep(str,'\,([0-9]*)\)','\,$1\.0\)');
fun = sym(['spline10' str '']);
end