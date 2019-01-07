% Run SBML tests from sbml-semantic-test-cases

function runSBMLTests
exedir = fileparts(mfilename('fullpath'));
cd(exedir);
fid = fopen(['./SBMLTest_log_' date '.txt'],'w+');

for iTest = 21:1529
    cd(exedir);
    testid = [repmat('0',1,4-floor(log10(iTest))),num2str(iTest)];
    if(~exist(fullfile(pwd,'SBMLresults',[testid '-results.csv'])))
        try
            runSBMLTest(iTest,fid);
        catch error
            fprintf(fid,['Test ' testid ' failed: ' error.message '\n']);;
        end
    end
end
fclose(fid);
end


function runSBMLTest(iTest,fid)
cd(fileparts(mfilename('fullpath')))
curdir = pwd;
testid = [repmat('0',1,4-floor(log10(iTest))),num2str(iTest)];
disp([' =================== ' testid ' =================== ']);
if(exist(fullfile(pwd,'sbml-semantic-test-cases','cases','semantic',testid),'dir'))
    cd(fullfile(pwd,'sbml-semantic-test-cases','cases','semantic',testid))
    try
        if(exist(fullfile(pwd,[testid '-sbml-l3v1.xml']),'file'))
            SBML2AMICI([testid '-sbml-l3v1'],['SBMLTEST_' testid])
        elseif(exist(fullfile(pwd,[testid '-sbml-l2v4.xml'])))
            SBML2AMICI([testid '-sbml-l2v4'],['SBMLTEST_' testid])
        else
            return;
        end
        amiwrap(['SBMLTEST_' testid],['SBMLTEST_' testid '_syms'],pwd);
    catch error
        fprintf(fid,['Test ' testid ' failed: ' error.message '\n']);;
        return
    end
    load(['SBMLTEST_' testid '_knom.mat'])
    load(['SBMLTEST_' testid '_pnom.mat'])
    load(['SBMLTEST_' testid '_vnom.mat'])
    load(['SBMLTEST_' testid '_kvnom.mat'])
    [t,settings,concflag] = parseSettings(testid);
    options.sensi = 0;
    options.maxsteps = 1e6;
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
            if(~concflag)
                vol = vnom(ix);
                vol = subs(vol,model.sym.k(:),sym(knom(:)));
                vol = subs(vol,model.sym.p(:),sym(pnom(:)));
                vol = double(vol);
            else
                vol = 1;
            end
            amiresults{:,ispecies} = sol.x(:,ix)*vol;
        elseif(~isempty(ik))
            if(~concflag)
                vol = kvnom(ik);
                vol = subs(vol,model.sym.k(:),sym(knom(:)));
                vol = subs(vol,model.sym.p(:),sym(pnom(:)));
                vol = double(vol);
            else
                vol = 1;
            end
            amiresults{:,ispecies} = results{:,ispecies}*0 + knom(ik)*vol;
        end
        adev(:,ispecies-1) = abs(amiresults{:,ispecies}-results{:,ispecies});
        rdev(:,ispecies-1) = abs((amiresults{:,ispecies}-results{:,ispecies})./results{:,ispecies});
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
tline = fgetl(fid);
if(~strcmp(tline,'concentration: ') && ~strcmp(tline,'concentration:'))
    concflag = true;
else
    concflag = false;
end
fclose(fid);

end
