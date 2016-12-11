classdef getFP_F4C_Test < matlab.unittest.TestCase
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 4B :: F = 1.5
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testFigure4C1(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 2;
            actRegime = getFP(1, 0.25, -1, 4, -1, 1, 1.5);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        function testFigure4C2(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 4;
            actRegime = getFP(1, 0.5, -1, 4, -1, 1, 1.5);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        function testFigure4C3(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 1;
            actRegime = getFP(1, 1, -1, 4, -1, 1, 1.5);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        function testFigure4C4(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 4;
            actRegime = getFP(1, 1.5, -1, 4, -1, 1, 1.5);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        function testFigure4C5(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 2;
            actRegime = getFP(1, 1.75, -1, 4, -1, 1, 1.5);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 4B :: F = 0.3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
    end
    
end
