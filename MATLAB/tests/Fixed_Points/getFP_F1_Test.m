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
            
            % Bifurcation parameters
            w = 1;
            w0 = 0; % We don't want the input involved
            alpha = 0;
            beta1 = -100;
            beta2 = 0;
            epsilon = 0;
            F = 0; % We don't want the input involved
            
            [actRstar, actPsiStar, actStability] = ...
                                    getFP(w, w0,alpha, beta1, beta2, epsilon, F);
            
            actNumFP = length(actRstar);               
                               
            testCase.verifyEqual(actRstar,expRstar,'AbsTol',rError);
            testCase.verifyEqual(actNumFP,expNumFP);
            
            % Display
            disp(['The expected R_Star is: ', num2str(expRstar)])
            disp(['The expected number of FP(s) is: ', num2str(expNumFP)])

        end
        
        function testFigure1B(testCase)
            
            rError = 0.02; % error margin
            expRstar1 = 0.0;
            expRstar2 = 0.1;
            expNumFP = 2; % expected number of fixed points
            
            w0 = 0; % We don't want the input involved
            F = 0; % We don't want the input involved
            w = 1;

            % Bifurcation parameters
            alpha = 1;
            beta1 = -100;
            beta2 = 0;
            epsilon = 0;
            
            [actRstar, actPsiStar, actStability] = ...
                                    getFP(w, w0,alpha, beta1, beta2, epsilon, F);
            
            actNumFP = size(actRstar);               
                               
            testCase.verifyEqual(actRstar(1),expRstar1,'AbsTol',rError);
            testCase.verifyEqual(actRstar(2),expRstar2,'AbsTol',rError);
            testCase.verifyEqual(actNumFP,expNumFP);
            
            % Display
            disp(['The expected R_Star is: ', num2str(expRstar)])
            disp(['The expected number of FP(s) is: ', num2str(expNumFP)])

        end
        
    end
end
