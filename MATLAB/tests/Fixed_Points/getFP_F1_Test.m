classdef getFP_F1_Test < matlab.unittest.TestCase
    % getFP_F1_Test 
    %   Fixed Point Location and Stability Tests from 
    %   Kim & Large 2015 figure 1
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
        % Figure 1
        %%%%%%%%%%%%%%%
        
        function testFigure1A(testCase)
            rError = 0.02; % error margin
            
            expRstar = 0.0;
            expNumFP = 1;
            
            f0 = 0; % We don't want the input involved
            F = 0; % We don't want the input involved
            f = 1; % osc freq
            
            % Bifurcation parameters
            alpha = 0;
            beta1 = -100;
            beta2 = 0;
            epsilon = 0;
            
            [actRstar, actPsiStar, actStability, actRegime] = ...
                                    getFP(f, f0,alpha, beta1, beta2, epsilon, F);
            
            actNumFP = length(actRstar);               
                               
            testCase.verifyEqual(actRstar,expRstar,'AbsTol',rError);
            testCase.verifyEqual(actNumFP,expNumFP);
            
            % Display
            disp(['The actual R_Star is: ', num2str(actRstar)])
            disp(['The expected R_Star is: ', num2str(actRstar)])
            disp(['The expected number of FP(s) is: ', num2str(expNumFP)])
            disp(['The actual number of FP(s) is: ', num2str(actRegime)])
            fprintf('\n')

        end
        
        function testFigure1B(testCase)
            
            % Bifurcation parameters
            alpha = 0.5;
            beta1 = -1;
            beta2 = 0;
            epsilon = 0;
            
            rError = 0.02; % error margin
            expRstar1 = 0.0;
            expRstar2 = sqrt(-alpha/beta1);
            expNumFP = 2; % expected number of fixed points
            
            f0 = 0; % We don't want the input involved
            F = 0; % We don't want the input involved
            f = 1; % osc freq
            
            [actRstar, actPsiStar, actStability, actRegime] = ...
                                    getFP(f, f0,alpha, beta1, beta2, epsilon, F);
            
            actNumFP = length(actRstar);               
                               
            testCase.verifyEqual(actRstar(1),expRstar1,'AbsTol',rError);
            testCase.verifyEqual(actRstar(2),expRstar2,'AbsTol',rError);
            testCase.verifyEqual(actNumFP,expNumFP);
            
            % Display
            disp(['The actual R_Star is: ', num2str(actRstar(1))])
            disp(['The actual R_Star is: ', num2str(actRstar(2))])
            disp(['The expected R_Star is: ', num2str(expRstar1)])
            disp(['The expected R_Star is: ', num2str(expRstar2)])
            disp(['The expected number of FP(s) is: ', num2str(expNumFP)])
            disp(['The actual number of FP(s) is: ', num2str(actRegime)])
            fprintf('\n')
            
        end
        
        function testFigure1C(testCase)
            
            % Bifurcation parameters
            alpha = -1;
            beta1 = 4;
            beta2 = -1;
            epsilon = 1;
            
            rError = 0.02; % error margin
            expRstar1 = 0.0;
            expRstar2 = sqrt(-alpha/beta1);
            expNumFP = 3; % expected number of fixed points
            
            f0 = 0; % We don't want the input involved
            F = 0; % We don't want the input involved
            f = 1; % osc freq
            
            [actRstar, actPsiStar, actStability, actRegime] = ...
                                    getFP(f, f0,alpha, beta1, beta2, epsilon, F);
            
            actNumFP = length(actRstar);               
                               
%             testCase.verifyEqual(actRstar(1),expRstar1,'AbsTol',rError);
%             testCase.verifyEqual(actRstar(2),expRstar2,'AbsTol',rError);
            testCase.verifyEqual(actNumFP,expNumFP);
            
            % Display
            disp(['The expected number of FP(s) is: ', num2str(expNumFP)])
            disp(['The actual number of FP(s) is: ', num2str(actRegime)])
            fprintf('\n')
            
        end
        
        function testFigure1D(testCase)
            
            % Bifurcation parameters
            alpha = -1;
            beta1 = 2.5;
            beta2 = -1;
            epsilon = 1;
            
            rError = 0.02; % error margin
            expRstar1 = 0.0;
            expRstar2 = sqrt(-alpha/beta1);
            expNumFP = 1; % expected number of fixed points
            
            f0 = 0; % We don't want the input involved
            F = 0; % We don't want the input involved
            f = 1; % osc freq
            
            [actRstar, actPsiStar, actStability, actRegime] = ...
                                    getFP(f, f0,alpha, beta1, beta2, epsilon, F);
            
            actNumFP = length(actRstar);               
                               
%             testCase.verifyEqual(actRstar(1),expRstar1,'AbsTol',rError);
%             testCase.verifyEqual(actRstar(2),expRstar2,'AbsTol',rError);
            testCase.verifyEqual(actNumFP,expNumFP);
            
            % Display
            disp(['The expected number of FP(s) is: ', num2str(expNumFP)])
            disp(['The actual number of FP(s) is: ', num2str(actNumFP)])
            fprintf('\n')
            
        end
        
    end
end
