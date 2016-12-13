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
            actRegime = getOscParamRegime(alpha,beta1,beta2,epsilon);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            fprintf('\n')
            disp(['alpha is: ', num2str(alpha)])
            disp(['beta1 is: ', num2str(beta1)])
            disp(['beta2 is: ', num2str(beta2)])
            disp(['epsilon is: ', num2str(epsilon)])
            fprintf('\n')
            disp(['The expected regime is: ', num2str(actRegime)])
            disp(['The actual regime is: ', num2str(expRegime)])

        end
        
        function testLargeNegBetaCritical(testCase)

            alpha = 0;
            beta1 = -1000;
            beta2 = 0;
            epsilon = 0;

            expRegime = 1;
            actRegime = getOscParamRegime(alpha,beta1,beta2,epsilon);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            fprintf('\n')
            disp(['alpha is: ', num2str(alpha)])
            disp(['beta1 is: ', num2str(beta1)])
            disp(['beta2 is: ', num2str(beta2)])
            disp(['epsilon is: ', num2str(epsilon)])
            fprintf('\n')
            disp(['The expected regime is: ', num2str(actRegime)])
            disp(['The actual regime is: ', num2str(expRegime)])

        end
        
        function testWithNonZeroEpsilonSuperCritical(testCase)
            
            alpha = 0;
            beta1 = -1;
            beta2 = 0;
            epsilon = 1;
            
            expRegime = 1;
            actRegime = getOscParamRegime(alpha,beta1,beta2,epsilon);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            fprintf('\n')
            disp(['alpha is: ', num2str(alpha)])
            disp(['beta1 is: ', num2str(beta1)])
            disp(['beta2 is: ', num2str(beta2)])
            disp(['epsilon is: ', num2str(epsilon)])
            fprintf('\n')
            disp(['The expected regime is: ', num2str(actRegime)])
            disp(['The actual regime is: ', num2str(expRegime)])

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
            actRegime = getOscParamRegime(alpha,beta1,beta2,epsilon);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            fprintf('\n')
            disp(['alpha is: ', num2str(alpha)])
            disp(['beta1 is: ', num2str(beta1)])
            disp(['beta2 is: ', num2str(beta2)])
            disp(['epsilon is: ', num2str(epsilon)])
            fprintf('\n')
            disp(['The expected regime is: ', num2str(actRegime)])
            disp(['The actual regime is: ', num2str(expRegime)])

        end
        
        function testLargeAlphaSuperCritical(testCase)
            
            alpha = 1000;
            beta1 = -1;
            beta2 = 0;
            epsilon = 0;
            
            expRegime = 2;
            actRegime = getOscParamRegime(alpha,beta1,beta2,epsilon);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            fprintf('\n')
            disp(['alpha is: ', num2str(alpha)])
            disp(['beta1 is: ', num2str(beta1)])
            disp(['beta2 is: ', num2str(beta2)])
            disp(['epsilon is: ', num2str(epsilon)])
            fprintf('\n')
            disp(['The expected regime is: ', num2str(actRegime)])
            disp(['The actual regime is: ', num2str(expRegime)])

        end
 
         function testLargeNegBeta1SuperCritical(testCase)
            
            alpha = 1;
            beta1 = -1000;
            beta2 = 0;
            epsilon = 0;
            
            expRegime = 2;
            actRegime = getOscParamRegime(alpha,beta1,beta2,epsilon);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            fprintf('\n')
            disp(['alpha is: ', num2str(alpha)])
            disp(['beta1 is: ', num2str(beta1)])
            disp(['beta2 is: ', num2str(beta2)])
            disp(['epsilon is: ', num2str(epsilon)])
            fprintf('\n')
            disp(['The expected regime is: ', num2str(actRegime)])
            disp(['The actual regime is: ', num2str(expRegime)])

         end
         
         function testLargeAlphaAndBeta1SuperCritical(testCase)
            
            alpha = 500;
            beta1 = -1000;
            beta2 = 0;
            epsilon = 0;
            
            expRegime = 2;
            actRegime = getOscParamRegime(alpha,beta1,beta2,epsilon);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            fprintf('\n')
            disp(['alpha is: ', num2str(alpha)])
            disp(['beta1 is: ', num2str(beta1)])
            disp(['beta2 is: ', num2str(beta2)])
            disp(['epsilon is: ', num2str(epsilon)])
            fprintf('\n')
            disp(['The expected regime is: ', num2str(actRegime)])
            disp(['The actual regime is: ', num2str(expRegime)])

         end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Supercritical Double Limit Cycle : Regime 3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testBasicSuperCriticalDLC(testCase)
            
            alpha = -1;
            beta1 = 1;
            beta2 = -1;
            epsilon = 1;
            
            expRegime = 3; 
            actRegime = getOscParamRegime(alpha,beta1,beta2,epsilon);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            fprintf('\n')
            disp(['alpha is: ', num2str(alpha)])
            disp(['beta1 is: ', num2str(beta1)])
            disp(['beta2 is: ', num2str(beta2)])
            disp(['epsilon is: ', num2str(epsilon)])
            fprintf('\n')
            disp(['The expected regime is: ', num2str(actRegime)])
            disp(['The actual regime is: ', num2str(expRegime)])

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Subcritical Double Limit Cycle : Regime 4
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testLargeEpsilonSubCriticalDLC(testCase)
            
            alpha = -1;
            beta1 = 1;
            beta2 = -1;
            epsilon = 1000;
            
            expRegime = 4; 
            actRegime = getOscParamRegime(alpha,beta1,beta2,epsilon);
            testCase.verifyEqual(actRegime,expRegime);
            
            % Display
            fprintf('\n')
            disp(['alpha is: ', num2str(alpha)])
            disp(['beta1 is: ', num2str(beta1)])
            disp(['beta2 is: ', num2str(beta2)])
            disp(['epsilon is: ', num2str(epsilon)])
            fprintf('\n')
            disp(['The expected regime is: ', num2str(actRegime)])
            disp(['The actual regime is: ', num2str(expRegime)])

       end
        
    end
end
