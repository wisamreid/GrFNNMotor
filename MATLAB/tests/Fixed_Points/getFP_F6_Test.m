classdef getFP_F6_Test < matlab.unittest.TestCase
    % getFP_F6_Test 
    %   Fixed Point Location and Stability Tests from 
    %   Kim & Large 2015 figure 6
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
        % Figure 6
        %%%%%%%%%%%%%%%
        
        function testFigure6A1(testCase)
            rError = 0.02; % error margin
            psiError = pi/12; % error margin
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 4;
            expRstar = 0.06;
            expPsiStar = pi/2;
            [actRstar, actPsiStar, actStability] = getFP(1, 0.5, 1, -100, 0, 0, 0.2);
            testCase.verifyEqual(actStability,expStability);
            testCase.verifyEqual(actRstar,expRstar,'AbsTol',rError);
            testCase.verifyEqual(actPsiStar,expPsiStar,'AbsTol',psiError);
            
            % Display
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability)))
            disp(['The expected R_Star is: ', num2str(expRstar)])
            disp(['The expected Psi_Stars is: ', num2str(expPsiStar)])

        end
        
        function testFigure6A2(testCase)
            rError = 0.02; % error margin
            psiError = pi/12; % error margin
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 4;
            expRstar = 0.04;
            expPsiStar = pi/2;
            [actRstar, actPsiStar, actStability] = getFP(1, 0.3, 1, -100, 0, 0, 0.2);
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
