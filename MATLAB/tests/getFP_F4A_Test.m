classdef getFP_F4A_Test < matlab.unittest.TestCase
    % getFPTest 
    %   Fixed Point Location and Stability Tests from Kim & Large 2015 figures
    %   Note: The expected values are determined by inspection from the figures
    %   and are not ground truth
    
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
        
        %%%%%%%%%%%%%%%
        % Figure 4A
        %%%%%%%%%%%%%%%
        
        function testFigure4A1(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 2;
            actRegime = getFP(1, 0.6, 0, -100, 0, 0, 0.2);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        function testFigure4A2(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 1;
            actRegime = getFP(1, 1, 0, -100, 0, 0, 0.2);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        function testFigure4A3(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 2;
            actRegime = getFP(1, 1.4, 0, -100, 0, 0, 0.2);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
    end
    
end
