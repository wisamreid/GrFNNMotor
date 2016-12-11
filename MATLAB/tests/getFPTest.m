classdef getFPTest < matlab.unittest.Testcase
    % ExampleTest 
    %   Template for Matlab Unittests
    
    properties
        OriginalPath
    end
       
    methods (TestMethodSetup)
        function addLibToPath(testCase)
            testCase.OriginalPath = path;
            addpath(fullfile(pwd,'../lib'));
        end
    end
    
    methods (Test)
        function testFigure3A1(testCase)
            marginOfError = 0.01;
            expRegime = 1;
            expRstar = 0.125;
            expPsiStar = 0;
            [actRegime,actRstar, actPsiStar] = getFP(1, 1, 0, -100, 0, 0, 0.2);
            testCase.verifyEqual(actRegime,expRegime,'AbsTol',marginOfError);
            testCase.verifyEqual(actRstar,expRstar,'AbsTol',marginOfError);
            testCase.verifyEqual(actPsiStar,expPsiStar,'AbsTol',marginOfError);

        end
    end
    
end

