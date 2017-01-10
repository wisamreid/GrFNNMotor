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
            expStability = 1;
            
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
            testCase.verifyEqual(actStability,expStability);

            
            % Display
            disp(['The actual R_Star is: ', num2str(actRstar)])
            disp(['The expected R_Star is: ', num2str(actRstar)])
            disp(['The expected number of FP(s) is: ', num2str(expNumFP)])
            disp(['The actual number of FP(s) is: ', num2str(actRegime)])
            disp(['The expected stability of FP(s) is: ', num2str(expStability)])
            disp(['The actual stability of FP(s) is: ', num2str(actStability)])
            fprintf('\n')

        end
        
        function testFigure1B(testCase)
            
            % Bifurcation parameters
            alpha = 0.5;
            beta1 = -1;
            beta2 = 0;
            epsilon = 0;
            
            rError = 0.02; % error margin
            expRstar = [0.0; sqrt(-alpha/beta1)];
            expNumFP = 2; % expected number of fixed points
            expStability = [3; 1]; % expected FP stability
            
            f0 = 0; % We don't want the input involved
            F = 0; % We don't want the input involved
            f = 1; % osc freq
            
            [actRstar, actPsiStar, actStability, actRegime] = ...
                                    getFP(f, f0,alpha, beta1, beta2, epsilon, F);
            
            actNumFP = length(actRstar);
                               
            testCase.verifyEqual(actRstar(1),expRstar(1),'AbsTol',rError);
            testCase.verifyEqual(actRstar(2),expRstar(2),'AbsTol',rError);
            testCase.verifyEqual(actNumFP,expNumFP);
            testCase.verifyEqual(actStability(1),expStability(1),'AbsTol',rError);
            testCase.verifyEqual(actStability(2),expStability(2),'AbsTol',rError);
            
            % Display
            disp(['The actual R_Star is: ', num2str(actRstar(1))])
            disp(['The expected R_Star is: ', num2str(expRstar(1))])
            disp(['The actual R_Star is: ', num2str(actRstar(2))])
            disp(['The expected R_Star is: ', num2str(expRstar(2))])
            disp(['The expected number of FP(s) is: ', num2str(expNumFP)])
            disp(['The actual number of FP(s) is: ', num2str(actNumFP)])
            disp(['The expected stability of FP 1 is: ', num2str(expStability(1))])
            disp(['The actual stability of FP 1 is: ', num2str(actStability(1))])
            disp(['The expected stability of FP 2 is: ', num2str(expStability(2))])
            disp(['The actual stability of FP 2 is: ', num2str(actStability(2))])
            fprintf('\n')
            
        end
        
        function testFigure1C(testCase)
            
            % Bifurcation parameters
            alpha = -1;
            beta1 = 4;
            beta2 = -1;
            epsilon = 1;
            
            rError = 0.02; % error margin
            expRstar = [0.0; 0.0; 0.0];
            expNumFP = 3; % expected number of fixed points
            expStability = [1; 3; 1]; % expected FP stability
            
            f0 = 0; % We don't want the input involved
            F = 0; % We don't want the input involved
            f = 1; % osc freq
            
            [actRstar, actPsiStar, actStability, actRegime] = ...
                                    getFP(f, f0,alpha, beta1, beta2, epsilon, F);
            
            actNumFP = length(actRstar);               
                               
            testCase.verifyEqual(actRstar(1),expRstar(1),'AbsTol',rError);
%             testCase.verifyEqual(actRstar(2),expRstar(2),'AbsTol',rError);
%             testCase.verifyEqual(actRstar(3),expRstar(3),'AbsTol',rError);
            testCase.verifyEqual(actNumFP,expNumFP);
            testCase.verifyEqual(actStability(1),expStability(1),'AbsTol',rError);
            testCase.verifyEqual(actStability(2),expStability(2),'AbsTol',rError);
            testCase.verifyEqual(actStability(3),expStability(3),'AbsTol',rError);            
            
            % Display
            
            % non- zero r_star
            disp('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
            disp('Warning: Not Testing none zero r_star yet')
            disp('Values are needed to test against')
            disp('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')            
            fprintf('\n')
            disp(['The actual R_Star is: ', num2str(actRstar(2))])
            disp(['The expected R_Star is: ', num2str(expRstar(2))])
            disp(['The actual R_Star is: ', num2str(actRstar(3))])
            disp(['The expected R_Star is: ', num2str(expRstar(3))])
            fprintf('\n')
            
            % r_star @ zero
            disp(['The actual R_Star is: ', num2str(actRstar(1))])
            disp(['The expected R_Star is: ', num2str(expRstar(1))])
            
            % number of FPs
            disp(['The expected number of FP(s) is: ', num2str(expNumFP)])
            disp(['The actual number of FP(s) is: ', num2str(actRegime)])
            % stability
            disp(['The expected stability of FP 1 is: ', num2str(expStability(1))])
            disp(['The actual stability of FP 1 is: ', num2str(actStability(1))])
            disp(['The expected stability of FP 2 is: ', num2str(expStability(2))])
            disp(['The actual stability of FP 2 is: ', num2str(actStability(2))])
            fprintf('\n')
            
        end
        
        function testFigure1D(testCase)
            
            % Bifurcation parameters
            alpha = -1;
            beta1 = 2.5;
            beta2 = -1;
            epsilon = 1;
            
            rError = 0.02; % error margin
            expRstar = 0.0;
            expNumFP = 1; % expected number of fixed points
            expStability = 1;
            
            f0 = 0; % We don't want the input involved
            F = 0; % We don't want the input involved
            f = 1; % osc freq
            
            [actRstar, actPsiStar, actStability, actRegime] = ...
                                    getFP(f, f0,alpha, beta1, beta2, epsilon, F);
            
            actNumFP = length(actRstar);               
                               
            testCase.verifyEqual(actRstar,expRstar,'AbsTol',rError);
            testCase.verifyEqual(actNumFP,expNumFP);
            testCase.verifyEqual(actStability,expStability,'AbsTol',rError);
            
            % r_star
            disp(['The actual R_Star is: ', num2str(actRstar)])
            disp(['The expected R_Star is: ', num2str(expRstar)])
            
            % Display
            disp(['The expected number of FP(s) is: ', num2str(expNumFP)])
            disp(['The actual number of FP(s) is: ', num2str(actNumFP)])
            % stability
            disp(['The expected stability of FP 1 is: ', num2str(expStability)])
            disp(['The actual stability of FP 1 is: ', num2str(actStability)])

            fprintf('\n')
            
        end
        
    end
end
