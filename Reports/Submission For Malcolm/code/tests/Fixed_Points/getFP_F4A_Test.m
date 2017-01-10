classdef getFP_F4A_Test < matlab.unittest.TestCase
    % getFP_F4A_Test 
    %   Stability Type Tests from 
    %   Kim & Large 2015 figure 4A
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
        
        %%%%%%%%%%%%%%%
        % Figure 4A
        %%%%%%%%%%%%%%%
        
        function testFigure4A1(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 0.6, 0, -100, 0, 0, 0.2);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4A2(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 1;
            [r_star, psi_star, actStability] = getFP(1, 1, 0, -100, 0, 0, 0.2);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4A3(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 1.4, 0, -100, 0, 0, 0.2);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
    end
    
end
