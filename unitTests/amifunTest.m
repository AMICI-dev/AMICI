function tests = amifunTest
    tests = functiontests(localfunctions);
end

%% Test Functions
function testConstructor(testCase)
    for ifun = 1:length(testCase.TestData.funs)
        try
            a = amifun(testCase.TestData.funs{ifun},testCase.TestData.model);
        catch
            a = [];
        end
        verifyEqual(testCase,isempty(a),false,'Generation of amifun failed');
    end
end

%% Optional fresh fixtures  
function setup(testCase)  % do not change function name
    
    testCase.TestData(1).model = amimodel('example_model_1_syms','test1');
    for itest = 1:length(testCase.TestData)
        testCase.TestData(itest).funs = {'xdot','w','dwdx','dwdp','J','x0','Jv','JBand','JSparse','y','z','deltax','dydp','dxdotdp','root','Jy','dJydx','dJydp','sJy','Jz','dJzdx','dJzdp','sJz'};
    end
    
end
