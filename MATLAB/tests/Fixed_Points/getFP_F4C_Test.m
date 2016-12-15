classdef getFP_F4C_Test < matlab.unittest.TestCase
    % getFP_F4C_Test 
    %   Stability Type Tests from 
    %   Kim & Large 2015 figure 4C
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
        % Figure 4C :: F = 1.5
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testFigure4C1(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 0.25, -1, 4, -1, 1, 1.5);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4C2(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 4;
            [r_star, psi_star, actStability] = getFP(1, 0.5, -1, 4, -1, 1, 1.5);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4C3(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 0.7, -1, 4, -1, 1, 1.5);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4C4(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 1;
            [r_star, psi_star, actStability] = getFP(1, 1, -1, 4, -1, 1, 1.5);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4C5(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 1.3, -1, 4, -1, 1, 1.5);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4C6(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 4;
            [r_star, psi_star, actStability] = getFP(1, 1.5, -1, 4, -1, 1, 1.5);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4C7(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 1.75, -1, 4, -1, 1, 1.5);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 4C :: F = 0.3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testFigure4C8(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 0.5, -1, 4, -1, 1, 0.3);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4C9(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 3;
            [r_star, psi_star, actStability] = getFP(1, 0.91, -1, 4, -1, 1, 0.3);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4C10(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = [1; 5; 3];
            [r_star, psi_star, actStability] = getFP(1, 1, -1, 4, -1, 1, 0.3);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4C11(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 3;
            [r_star, psi_star, actStability] = getFP(1, 1.09, -1, 4, -1, 1, 0.3);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        function testFigure4C12(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 1.5, -1, 4, -1, 1, 0.3);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Figure 4C :: F = 0.1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testFigure4C13(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 0.5, -1, 4, -1, 1, 0.1);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end

        function testFigure4C14(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = [1; 5; 3; 5; 1];
            [r_star, psi_star, actStability] = getFP(1, 1, -1, 4, -1, 1, 0.1);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end       
        
        function testFigure4C15(testCase)
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            [r_star, psi_star, actStability] = getFP(1, 1.5, -1, 4, -1, 1, 0.1);
            testCase.verifyEqual(actStability,expStability);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))

        end
        
    end
    
end
