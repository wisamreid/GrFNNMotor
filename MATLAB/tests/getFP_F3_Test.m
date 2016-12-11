classdef getFP_F3_Test < matlab.unittest.TestCase
    % getFPTest 
    %   Fixed Point Location and Stability Tests from Kim & Large 2015 figures
    %   Note: The expected values are determined by inspection from the figures
    %   and are not ground truth
    
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
        
        %%%%%%%%%%%%%%%
        % Figure 3
        %%%%%%%%%%%%%%%
        
        function testFigure3A1(testCase)
            rError = 0.02;
            psiError = pi/12;
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 1;
            expRstar = 0.125;
            expPsiStar = pi/8;
            [actRegime, actRstar, actPsiStar] = getFP(1, 0.9, 0, -100, 0, 0, 0.2);
            testCase.verifyEqual(actRegime,expRegime);
            testCase.verifyEqual(actRstar,expRstar,'AbsTol',rError);
            testCase.verifyEqual(actPsiStar,expPsiStar,'AbsTol',psiError);
            
            % Display
            disp(strcat('The expected regime is: a ', regimeOptions(expRegime)))
            disp(['The expected R_Star is: ', num2str(expRstar)])
            disp(['The expected Psi_Star is: ', num2str(expPsiStar)])

        end
        
        function testFigure3A2(testCase)
            rError = 0.02;
            psiError = pi/12;
            regimeOptions = {' stable node',' stable spiral',' unstable node', ...
                ' unstable spiral',' saddle point'};
            expRegime = 2;
            expRstar = 0.06;
            expPsiStar = pi/2;
            [actRegime, actRstar, actPsiStar] = getFP(1, 0.5, 0, -100, 0, 0, 0.2);
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
