classdef getFP_F8_Test < matlab.unittest.TestCase
    % getFP_F8_Test 
    %   Fixed Point Location and Stability Tests from 
    %   Kim & Large 2015 figure 8
    %   
    %   Note: The expected r_star and psi_star values are determined by 
    %   inspection from the figures and are not ground truth
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
        % Figure 8
        %%%%%%%%%%%%%%%
        
        function testFigure8A1(testCase)
            rError = 0.02; % error margin
            psiError = pi/12; % error margin
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 1;
            expRstar = 0.85;
            expPsiStar = pi/8;
            [actRstar, actPsiStar, actStability] = getFP(1, 0.98, -1, 4, -1, 1, 0.3);
            testCase.verifyEqual(actStability,expStability);
            testCase.verifyEqual(actRstar,expRstar,'AbsTol',rError);
            testCase.verifyEqual(actPsiStar,expPsiStar,'AbsTol',psiError);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))
            disp(['The expected R_Star is: ', num2str(expRstar)])
            disp(['The expected Psi_Star is: ', num2str(expPsiStar)])

        end
        
        function testFigure8A2(testCase)
            rError = 0.02; % error margin
            psiError = pi/12; % error margin
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 3;
            expRstar = 0.57;
            expPsiStar = 5*pi/8;
            [actRstar, actPsiStar, actStability] = getFP(1, 0.92, -1, 4, -1, 1, 0.3);
            testCase.verifyEqual(actStability,expStability);
            testCase.verifyEqual(actRstar,expRstar,'AbsTol',rError);
            testCase.verifyEqual(actPsiStar,expPsiStar,'AbsTol',psiError);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))
            disp(['The expected R_Star is: ', num2str(expRstar)])
            disp(['The expected Psi_Star is: ', num2str(expPsiStar)])

        end
        
        function testFigure8A3(testCase)
            rError = 0.02; % error margin
            psiError = pi/12; % error margin
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            expRstar = 0.2;
            expPsiStar = 3*pi/8;
            [actRstar, actPsiStar, actStability] = getFP(1, 0.8, -1, 4, -1, 1, 0.3);
            testCase.verifyEqual(actStability,expStability);
            testCase.verifyEqual(actRstar,expRstar,'AbsTol',rError);
            testCase.verifyEqual(actPsiStar,expPsiStar,'AbsTol',psiError);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))
            disp(['The expected R_Star is: ', num2str(expRstar)])
            disp(['The expected Psi_Star is: ', num2str(expPsiStar)])

        end
        
    end
end
