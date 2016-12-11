classdef getFP_F5_Test < matlab.unittest.TestCase
    % getFP_F5_Test 
    %   Fixed Point Location and Stability Tests from 
    %   Kim & Large 2015 figure 5
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
        % Figure 5
        %%%%%%%%%%%%%%%
        
        function testFigure5A1(testCase)
            rError = 0.02; % error margin
            psiError = pi/12; % error margin
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 1;
            expRstar = 0.11;
            expPsiStar = pi/4;
            [actRegime, actRstar, actPsiStar] = getFP(1, 0.98, 1, -100, 0, 0, 0.02);
            testCase.verifyEqual(actRegime,expRegime);
            testCase.verifyEqual(actRstar,expRstar,'AbsTol',rError);
            testCase.verifyEqual(actPsiStar,expPsiStar,'AbsTol',psiError);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))
            disp(['The expected R_Star is: ', num2str(expRstar)])
            disp(['The expected Psi_Star is: ', num2str(expPsiStar)])

        end
        
        function testFigure5A2(testCase)
            rError = 0.02; % error margin
            psiError = pi/12; % error margin
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 4;
            expRstar = 0.02;
            expPsiStar = pi/2;
            [actRegime, actRstar, actPsiStar] = getFP(1, 0.96, 1, -100, 0, 0, 0.02);
            testCase.verifyEqual(actRegime,expRegime);
            testCase.verifyEqual(actRstar,expRstar,'AbsTol',rError);
            testCase.verifyEqual(actPsiStar,expPsiStar,'AbsTol',psiError);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))
            disp(['The expected R_Star is: ', num2str(expRstar)])
            disp(['The expected Psi_Star is: ', num2str(expPsiStar)])

        end
        
    end
end
