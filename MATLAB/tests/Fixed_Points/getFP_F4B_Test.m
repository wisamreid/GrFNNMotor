classdef getFP_F4B_Test < matlab.unittest.TestCase
    % getFP_F4B_Test 
    %   Stability Regime Tests from 
    %   Kim & Large 2015 figure 4B
    %
    %   Author: Wisam Reid
    %   Email: wisam@ccrma.stanford.edu
    
    properties
        OriginalPath
    end
       
    methods (TestMethodSetup)
        function addLibToPath(testCase)
            testCase.OriginalPath = path;
            addpath(fullfile(pwd,'../../lib'));
        end
    end
    
    methods (Test)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 4B :: F = 0.2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testFigure4B1(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 4;
            actRegime = getFP(1, 0.4, 1, -100, 0, 0, 0.2);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        function testFigure4B2(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 2;
            actRegime = getFP(1, 0.6, 1, -100, 0, 0, 0.2);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        function testFigure4B3(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 1;
            actRegime = getFP(1, 1, 1, -100, 0, 0, 0.2);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        function testFigure4B4(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 2;
            actRegime = getFP(1, 1.4, 1, -100, 0, 0, 0.2);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        function testFigure4B5(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 4;
            actRegime = getFP(1, 1.6, 1, -100, 0, 0, 0.2);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 4B :: F = 0.02
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function testFigure4B6(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 4;
            actRegime = getFP(1, 0.6, 1, -100, 0, 0, 0.02);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        function testFigure4B7(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 1;
            actRegime = getFP(1, 1, 1, -100, 0, 0, 0.02);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        function testFigure4B8(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 4;
            actRegime = getFP(1, 1.4, 1, -100, 0, 0, 0.02);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
    end
    
end
