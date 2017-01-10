classdef getOscParamRegime_Test < matlab.unittest.TestCase
    % getOscParamRegime_Test 
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
        % Critical Hopf : Regime 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testBasicCritical(testCase)

            alpha = 0;
            beta1 = -1;
            beta2 = 0;
            epsilon = 0;
            
            expRegime = 1;            
            [actRegime, actR, actDrdt] = getOscParamRegime(alpha,beta1,beta2,epsilon);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            fprintf('\n')
            disp(['alpha is: ', num2str(alpha)])
            disp(['beta1 is: ', num2str(beta1)])
            disp(['beta2 is: ', num2str(beta2)])
            disp(['epsilon is: ', num2str(epsilon)])
            fprintf('\n')
            disp(['r is: ', num2str(actR)])
            disp(['drdt is: ', num2str(actDrdt)])            
            disp(['The expected regime is: ', num2str(expRegime)])
            disp(['The actual regime is: ', num2str(actRegime)])

        end
        
        function testLargeNegBetaCritical(testCase)

            alpha = 0;
            beta1 = -1000;
            beta2 = 0;
            epsilon = 0;

            expRegime = 1;
            [actRegime, actR, actDrdt] = getOscParamRegime(alpha,beta1,beta2,epsilon);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            fprintf('\n')
            disp(['alpha is: ', num2str(alpha)])
            disp(['beta1 is: ', num2str(beta1)])
            disp(['beta2 is: ', num2str(beta2)])
            disp(['epsilon is: ', num2str(epsilon)])
            fprintf('\n')
            disp(['r is: ', num2str(actR)])
            disp(['drdt is: ', num2str(actDrdt)])            
            disp(['The expected regime is: ', num2str(expRegime)])
            disp(['The actual regime is: ', num2str(actRegime)])

        end
        
        function testWithNonZeroEpsilonSuperCritical(testCase)
            
            alpha = 0;
            beta1 = -1;
            beta2 = 0;
            epsilon = 1;
            
            expRegime = 1;
            [actRegime, actR, actDrdt] = getOscParamRegime(alpha,beta1,beta2,epsilon);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            fprintf('\n')
            disp(['alpha is: ', num2str(alpha)])
            disp(['beta1 is: ', num2str(beta1)])
            disp(['beta2 is: ', num2str(beta2)])
            disp(['epsilon is: ', num2str(epsilon)])
            fprintf('\n')
            disp(['r is: ', num2str(actR)])
            disp(['drdt is: ', num2str(actDrdt)])            
            disp(['The expected regime is: ', num2str(expRegime)])
            disp(['The actual regime is: ', num2str(actRegime)])

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Supercritical Hopf : Regime 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testBasicSuperCritical(testCase)
            
            alpha = 1;
            beta1 = -1;
            beta2 = 0;
            epsilon = 0;
            
            expRegime = 2; 
            [actRegime, actR, actDrdt] = getOscParamRegime(alpha,beta1,beta2,epsilon);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            fprintf('\n')
            disp(['alpha is: ', num2str(alpha)])
            disp(['beta1 is: ', num2str(beta1)])
            disp(['beta2 is: ', num2str(beta2)])
            disp(['epsilon is: ', num2str(epsilon)])
            fprintf('\n')
            disp(['r is: ', num2str(actR)])
            disp(['drdt is: ', num2str(actDrdt)])            
            disp(['The expected regime is: ', num2str(expRegime)])
            disp(['The actual regime is: ', num2str(actRegime)])

        end
        
        function testLargeAlphaSuperCritical(testCase)
            
            alpha = 1000;
            beta1 = -1;
            beta2 = 0;
            epsilon = 0;
            
            expRegime = 2;
            [actRegime, actR, actDrdt] = getOscParamRegime(alpha,beta1,beta2,epsilon);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            fprintf('\n')
            disp(['alpha is: ', num2str(alpha)])
            disp(['beta1 is: ', num2str(beta1)])
            disp(['beta2 is: ', num2str(beta2)])
            disp(['epsilon is: ', num2str(epsilon)])
            fprintf('\n')
            disp(['r is: ', num2str(actR)])
            disp(['drdt is: ', num2str(actDrdt)])            
            disp(['The expected regime is: ', num2str(expRegime)])
            disp(['The actual regime is: ', num2str(actRegime)])

        end
 
         function testLargeNegBeta1SuperCritical(testCase)
            
            alpha = 1;
            beta1 = -1000;
            beta2 = 0;
            epsilon = 0;
            
            expRegime = 2;
            [actRegime, actR, actDrdt] = getOscParamRegime(alpha,beta1,beta2,epsilon);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            fprintf('\n')
            disp(['alpha is: ', num2str(alpha)])
            disp(['beta1 is: ', num2str(beta1)])
            disp(['beta2 is: ', num2str(beta2)])
            disp(['epsilon is: ', num2str(epsilon)])
            fprintf('\n')
            disp(['r is: ', num2str(actR)])
            disp(['drdt is: ', num2str(actDrdt)])            
            disp(['The expected regime is: ', num2str(expRegime)])
            disp(['The actual regime is: ', num2str(actRegime)])

         end
         
         function testLargeAlphaAndBeta1SuperCritical(testCase)
            
            alpha = 500;
            beta1 = -1000;
            beta2 = 0;
            epsilon = 0;
            
            expRegime = 2;
            [actRegime, actR, actDrdt] = getOscParamRegime(alpha,beta1,beta2,epsilon);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            fprintf('\n')
            disp(['alpha is: ', num2str(alpha)])
            disp(['beta1 is: ', num2str(beta1)])
            disp(['beta2 is: ', num2str(beta2)])
            disp(['epsilon is: ', num2str(epsilon)])
            fprintf('\n')
            disp(['r is: ', num2str(actR)])
            disp(['drdt is: ', num2str(actDrdt)])            
            disp(['The expected regime is: ', num2str(expRegime)])
            disp(['The actual regime is: ', num2str(actRegime)])

         end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Supercritical Double Limit Cycle : Regime 3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testBasicSuperCriticalDLC(testCase)
            
            alpha = -1;
            beta1 = 4;
            beta2 = -1;
            epsilon = 1;
            
            expRegime = 3; 
            [actRegime, actR, actDrdt]  = getOscParamRegime(alpha,beta1,beta2,epsilon);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            fprintf('\n')
            disp(['alpha is: ', num2str(alpha)])
            disp(['beta1 is: ', num2str(beta1)])
            disp(['beta2 is: ', num2str(beta2)])
            disp(['epsilon is: ', num2str(epsilon)])
            fprintf('\n')
            disp(['r is: ', num2str(actR)])
            disp(['drdt is: ', num2str(actDrdt)])
            disp(['The expected regime is: ', num2str(expRegime)])
            disp(['The actual regime is: ', num2str(actRegime)])

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Subcritical Double Limit Cycle : Regime 4
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testBasicSubCriticalDLC(testCase)
            
            alpha = -1;
            beta1 = 2.5;
            beta2 = -1;
            epsilon = 1;
            
            expRegime = 4; % unreasonable epsilon value 
            [actRegime, actR, actDrdt] = getOscParamRegime(alpha,beta1,beta2,epsilon);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            fprintf('\n')
            disp(['alpha is: ', num2str(alpha)])
            disp(['beta1 is: ', num2str(beta1)])
            disp(['beta2 is: ', num2str(beta2)])
            disp(['epsilon is: ', num2str(epsilon)])
            fprintf('\n')
            disp(['r is: ', num2str(actR)])
            disp(['drdt is: ', num2str(actDrdt)])
            disp(['The expected regime is: ', num2str(expRegime)])
            disp(['The actual regime is: ', num2str(actRegime)])

       end
        
        function testLargeEpsilonSubCriticalDLC(testCase)
            
            alpha = -1;
            beta1 = 1;
            beta2 = -1;
            epsilon = 1000;
            
            expRegime = 0; % unreasonable epsilon value 
            [actRegime, actR, actDrdt] = getOscParamRegime(alpha,beta1,beta2,epsilon);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            fprintf('\n')
            disp(['alpha is: ', num2str(alpha)])
            disp(['beta1 is: ', num2str(beta1)])
            disp(['beta2 is: ', num2str(beta2)])
            disp(['epsilon is: ', num2str(epsilon)])
            fprintf('\n')
            disp(['r is: ', num2str(actR)])
            disp(['drdt is: ', num2str(actDrdt)])
            disp(['The expected regime is: ', num2str(expRegime)])
            disp(['The actual regime is: ', num2str(actRegime)])

       end
        
    end
end
