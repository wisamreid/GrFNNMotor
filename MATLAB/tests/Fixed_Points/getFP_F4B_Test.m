classdef getFP_F4B_Test < matlab.unittest.TestCase
    % getFP_F4B_Test 
    %   Stability Type Tests from 
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
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 4;
            [r_star, psi_star, actStability] = getFP(1, 0.4, 1, -100, 0, 0, 0.2);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4B2(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 0.6, 1, -100, 0, 0, 0.2);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4B3(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 1;
            [r_star, psi_star, actStability] = getFP(1, 1, 1, -100, 0, 0, 0.2);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4B4(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 1.4, 1, -100, 0, 0, 0.2);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4B5(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 4;
            [r_star, psi_star, actStability] = getFP(1, 1.6, 1, -100, 0, 0, 0.2);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 4B :: F = 0.02
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function testFigure4B6(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 4;
            [r_star, psi_star, actStability] = getFP(1, 0.6, 1, -100, 0, 0, 0.02);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4B7(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = [1; 5; 3];
            [r_star, psi_star, actStability] = getFP(1, 1, 1, -100, 0, 0, 0.02);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4B8(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 4;
            [r_star, psi_star, actStability] = getFP(1, 1.4, 1, -100, 0, 0, 0.02);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
    end
    
end
