classdef getFP_F4D_Test < matlab.unittest.TestCase
    % getFP_F4D_Test 
    %   Stability Regime Tests from 
    %   Kim & Large 2015 figure 4D
    %
    %   Author: Wisam Reid
    %   Email: wisam@ccrma.stanford.edu
    
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
        % Figure 4D :: F = 0.5
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testFigure4D1(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 2;
            actRegime = getFP(1, 0.8, -1, 2.5, -1, 1, 0.5);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        function testFigure4D2(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 4;
            actRegime = getFP(1, 0.87, -1, 2.5, -1, 1, 0.5);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        function testFigure4D3(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 2;
            actRegime = getFP(1, 0.89, -1, 2.5, -1, 1, 0.5);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        function testFigure4D4(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 1;
            actRegime = getFP(1, 1, -1, 2.5, -1, 1, 0.5);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        function testFigure4D5(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 2;
            actRegime = getFP(1, 1.11, -1, 2.5, -1, 1, 0.5);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        function testFigure4D6(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 4;
            actRegime = getFP(1, 1.13, -1, 2.5, -1, 1, 0.5);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        function testFigure4D7(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 2;
            actRegime = getFP(1, 1.2, -1, 2.5, -1, 1, 0.5);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 4D :: F = 0.2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testFigure4D8(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 2;
            actRegime = getFP(1, 0.8, -1, 4, -1, 1, 0.2);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        function testFigure4D9(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 5;
            actRegime = getFP(1, 1, -1, 4, -1, 1, 0.2);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))
            disp('WARNING! This Regime has 2 Fixed Points')
            
        end
        
        function testFigure4D10(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 2;
            actRegime = getFP(1, 1.2, -1, 4, -1, 1, 0.2);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 4D :: F = 0.1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testFigure4D11(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 2;
            actRegime = getFP(1, 0.8, -1, 4, -1, 1, 0.1);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end

        function testFigure4D12(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 1;
            actRegime = getFP(1, 1, -1, 4, -1, 1, 0.1);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end       
        
        function testFigure4D13(testCase)
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 2;
            actRegime = getFP(1, 1.2, -1, 4, -1, 1, 0.1);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))

        end
        
    end
    
end
