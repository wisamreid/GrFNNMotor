classdef getFP_F7_Test < matlab.unittest.TestCase
    % getFP_F7_Test 
    %   Fixed Point Location and Stability Tests from 
    %   Kim & Large 2015 figure 7
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
        % Figure 7
        %%%%%%%%%%%%%%%
        
        function testFigure7A1(testCase)
            rError = 0.02; % error margin
            psiError = pi/12; % error margin
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = [1; 5; 3; 5; 2];
            expRstar = [0.85; 0.85; 0.57; 0.45; 0.11];
            expPsiStar = [pi/6; 5*pi/6; 7*pi/8; pi/8; 0.0];
            [actRstar, actPsiStar, actStability] = getFP(1, 0.99, -1, 4, -1, 1, 0.1);
            testCase.verifyEqual(actStability,expStability);
            testCase.verifyEqual(actRstar,expRstar,'AbsTol',rError);
            testCase.verifyEqual(actPsiStar,expPsiStar,'AbsTol',psiError);
            
            % Display
            disp('For Fixed Point 1:')
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability(1))))
            disp(['The expected R_Star is: ', num2str(expRstar(1))])
            disp(['The expected Psi_Star is: ', num2str(expPsiStar(1))])
            fprintf('\n')
            
            disp('For Fixed Point 1:')
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability(2))))
            disp(['The expected R_Star is: ', num2str(expRstar(2))])
            disp(['The expected Psi_Star is: ', num2str(expPsiStar(2))])
            fprintf('\n')
            
            disp('For Fixed Point 1:')
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability(3))))
            disp(['The expected R_Star is: ', num2str(expRstar(3))])
            disp(['The expected Psi_Star is: ', num2str(expPsiStar(3))])
            fprintf('\n')
            
            disp('For Fixed Point 1:')
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability(4))))
            disp(['The expected R_Star is: ', num2str(expRstar(4))])
            disp(['The expected Psi_Star is: ', num2str(expPsiStar(4))])
            fprintf('\n')
            
            disp('For Fixed Point 1:')
            disp(strcat('The expected stability is: a ', stabilityOptions(expStability(5))))
            disp(['The expected R_Star is: ', num2str(expRstar(5))])
            disp(['The expected Psi_Star is: ', num2str(expPsiStar(5))])
            fprintf('\n')

        end
        
        function testFigure7A2(testCase)
            rError = 0.02; % error margin
            psiError = pi/12; % error margin
            stabilityOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expStability = 2;
            expRstar = 0.1;
            expPsiStar = pi/8;
            [actRstar, actPsiStar, actStability] = getFP(1, 0.95, -1, 4, -1, 1, 0.1);
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
