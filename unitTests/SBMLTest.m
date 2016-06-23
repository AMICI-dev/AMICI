function runSBMLTests
for iTest = 1:1196
    runSBMLTest(iTest)
end
end


function runSBMLTest(iTest)
    cd(fileparts(mfilename('fullpath')))
    curdir = pwd;
    testid = [repmat('0',1,4-floor(log10(iTest))),num2str(iTest)];
    cd(fullfile(pwd,'CustomSBMLTestsuite',testid))
    SBML2AMICI([testid '-sbml-l3v1'],['SBMLTEST_' testid])
    amiwrap(['SBMLTEST_' testid],['SBMLTEST_' testid '_syms'],pwd)
    load(['SBMLTEST_' testid '_knom.mat'])
    load(['SBMLTEST_' testid '_pnom.mat'])
    [t,options] = parseSettings(testid);
    eval(['sol = simulate_SBMLTEST_' testid '(t,pnom,knom,[],options);'])
    results = readtable([testid '-results.csv']);
    eval(['model = SBMLTEST_' testid '_syms;'])
    maxadev = 0;
    maxrdev = 0;
    for ispecies = 2:length(results.Properties.VariableNames)
        maxadevp = max(abs(sol.x(:,find(logical(sym(results.Properties.VariableNames{ispecies})==model.sym.x)))-results{:,ispecies}));
        rdev = abs((sol.x(:,find(logical(sym(results.Properties.VariableNames{ispecies})==model.sym.x)))-results{:,ispecies})./results{:,ispecies});
        rdev(isinf(rdev)) = 0;
        maxrdevp = max(rdev);
        if(maxadevp > maxadev)
            maxadev = maxadevp;
        end
        if(maxrdevp > maxrdev)
            maxrdev = maxrdevp;
        end
    end
    assert(not(and(maxadev>options.atol*1000,maxrdev>options.rtol*1000)))
    cd(curdir)
end

function [t,options] = parseSettings(testid)
fid = fopen([testid '-settings.txt']);
T = textscan(fid,'%s %f');
t = linspace(T{2}(1),T{2}(2),T{2}(3)+1);
tline = fgetl(fid);
T = textscan(fid,'%s %f');
options.atol = T{2}(1);
options.rtol = T{2}(2);
fclose(fid);

end