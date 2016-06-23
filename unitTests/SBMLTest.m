function runSBMLTests
for iTest = 1:1196
    runSBMLTest(iTest)
end
end


function runSBMLTest(iTest)
    cd(fileparts(mfilename('fullpath')))
    curdir = pwd;
    testid = [repmat('0',1,4-floor(log10(iTest))),num2str(iTest)];
    disp([' =================== ' testid ' =================== ']); 
    cd(fullfile(pwd,'CustomSBMLTestsuite',testid))
    SBML2AMICI([testid '-sbml-l3v1'],['SBMLTEST_' testid])
    amiwrap(['SBMLTEST_' testid],['SBMLTEST_' testid '_syms'],pwd)
    load(['SBMLTEST_' testid '_knom.mat'])
    load(['SBMLTEST_' testid '_pnom.mat'])
    [t,options] = parseSettings(testid);
    eval(['sol = simulate_SBMLTEST_' testid '(t,pnom,knom,[],options);'])
    results = readtable([testid '-results.csv']);
    eval(['model = SBMLTEST_' testid '_syms;'])
    adev = zeros(size(results{:,2:end}));
    rdev = zeros(size(results{:,2:end}));
    for ispecies = 2:length(results.Properties.VariableNames)
        if(any(logical(sym(results.Properties.VariableNames{ispecies})==model.sym.x)))
            adev(:,ispecies-1) = abs(sol.x(:,find(logical(sym(results.Properties.VariableNames{ispecies})==model.sym.x)))-results{:,ispecies});
            rdev(:,ispecies-1) = abs((sol.x(:,find(logical(sym(results.Properties.VariableNames{ispecies})==model.sym.x)))-results{:,ispecies})./results{:,ispecies});
        elseif(any(logical(sym(results.Properties.VariableNames{ispecies})==model.sym.k)))
            adev(:,ispecies-1) = abs(knom(find(logical(sym(results.Properties.VariableNames{ispecies})==model.sym.k)))-results{:,ispecies});
            rdev(:,ispecies-1) = abs((knom(find(logical(sym(results.Properties.VariableNames{ispecies})==model.sym.k)))-results{:,ispecies})./results{:,ispecies});
        end
    end
    rdev(isinf(rdev)) = 0;
    assert(not(any(any(and(adev>options.atol*1000,rdev>options.rtol*1000)))))
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