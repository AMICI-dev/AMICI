function runSBMLTests
for iTest = 175:1183
    try
    runSBMLTest(iTest)
    catch error
        disp(error.message);
    end
end
end


function runSBMLTest(iTest)
cd(fileparts(mfilename('fullpath')))
curdir = pwd;
testid = [repmat('0',1,4-floor(log10(iTest))),num2str(iTest)];
disp([' =================== ' testid ' =================== ']);
if(exist(fullfile(pwd,'CustomSBMLTestsuite',testid),'dir'))
    cd(fullfile(pwd,'CustomSBMLTestsuite',testid))
    try
        if(exist(fullfile(pwd,[testid '-sbml-l3v1.xml']),'file'))
            SBML2AMICI([testid '-sbml-l3v1'],['SBMLTEST_' testid])
        elseif(exist(fullfile(pwd,[testid '-sbml-l2v4.xml'])))
            SBML2AMICI([testid '-sbml-l2v4'],['SBMLTEST_' testid])
        end
        amiwrap(['SBMLTEST_' testid],['SBMLTEST_' testid '_syms'],pwd);
    catch error
        warning(['Test ' testid ' failed: ' error.message]);
        return
    end
    load(['SBMLTEST_' testid '_knom.mat'])
    load(['SBMLTEST_' testid '_pnom.mat'])
    load(['SBMLTEST_' testid '_vnom.mat'])
    [t,settings,concflag] = parseSettings(testid);
    options.sensi = 0;
    eval(['sol = simulate_SBMLTEST_' testid '(t,pnom,knom,[],options);'])
    results = readtable([testid '-results.csv']);
    eval(['model = SBMLTEST_' testid '_syms;'])
    adev = zeros(size(results{:,2:end}));
    rdev = zeros(size(results{:,2:end}));
    amiresults = results;
    for ispecies = 2:length(results.Properties.VariableNames)
        ix = find(logical(sym(results.Properties.VariableNames{ispecies})==model.sym.x));
        ik = find(logical(sym(results.Properties.VariableNames{ispecies})==model.sym.k));
        if(~isempty(ix))
            if(concflag)
                vol = vnom(ix);
                vol = subs(vol,model.sym.k(:),sym(knom(:)));
                vol = subs(vol,model.sym.p(:),sym(pnom(:)));
                vol = double(vol);
            else
                vol = 1;
            end
            amiresults{:,ispecies} = sol.x(:,ix)*vol;
            adev(:,ispecies-1) = abs(sol.x(:,ix)*vol-results{:,ispecies});
            rdev(:,ispecies-1) = abs((sol.x(:,ix)*vol-results{:,ispecies})./results{:,ispecies});
        elseif(~isempty(ik))
            amiresults{:,ispecies} = results{:,ispecies}*0 + knom(ik);
            
            adev(:,ispecies-1) = abs(knom(ik)-results{:,ispecies});
            rdev(:,ispecies-1) = abs((knom(ik)-results{:,ispecies})./results{:,ispecies});
        end
    end
    rdev(isinf(rdev)) = 0;
    assert(not(any(any(and(adev>settings.atol,rdev>settings.rtol)))))
    cd(curdir)
    mkdir(fullfile(pwd,'SBMLresults'))
    writetable(amiresults,fullfile(pwd,'SBMLresults',[testid '-results.csv']))
    clear all
end
end

function [t,options,concflag] = parseSettings(testid)
fid = fopen([testid '-settings.txt']);
T = textscan(fid,'%s %f');
t = linspace(T{2}(1),T{2}(2),T{2}(3)+1);
tline = fgetl(fid);
T = textscan(fid,'%s %f');
options.atol = T{2}(1);
options.rtol = T{2}(2);
tline = fgetl(fid);
if(strcmp(tline,''))
    concflag = false;
else
    concflag = true;
end
fclose(fid);

end