classdef getFP_F4D_Test < matlab.unittest.TestCase
    % getFP_F4D_Test 
    %   Stability stability Tests from 
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
            addpath(fullfile(pwd,'../../lib'));
        end
    end
    
    methods (Test)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 4D :: F = 0.5
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testFigure4D1(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 0.8, -1, 2.5, -1, 1, 0.5);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4D2(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 4;
            [r_star, psi_star, actStability] = getFP(1, 0.87, -1, 2.5, -1, 1, 0.5);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4D3(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 0.89, -1, 2.5, -1, 1, 0.5);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4D4(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 1;
            [r_star, psi_star, actStability] = getFP(1, 1, -1, 2.5, -1, 1, 0.5);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4D5(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 1.11, -1, 2.5, -1, 1, 0.5);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4D6(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 4;
            [r_star, psi_star, actStability] = getFP(1, 1.13, -1, 2.5, -1, 1, 0.5);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4D7(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 1.2, -1, 2.5, -1, 1, 0.5);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 4D :: F = 0.2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testFigure4D8(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 0.8, -1, 4, -1, 1, 0.2);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4D9(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 5;
            [r_star, psi_star, actStability] = getFP(1, 1, -1, 4, -1, 1, 0.2);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))
            disp('WARNING! This stability has 2 Fixed Points')
            
        end
        
        function testFigure4D10(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 1.2, -1, 4, -1, 1, 0.2);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 4D :: F = 0.1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testFigure4D11(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 0.8, -1, 4, -1, 1, 0.1);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end

        function testFigure4D12(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 1;
            [r_star, psi_star, actStability] = getFP(1, 1, -1, 4, -1, 1, 0.1);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end       
        
        function testFigure4D13(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 1.2, -1, 4, -1, 1, 0.1);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
    end
    
end
