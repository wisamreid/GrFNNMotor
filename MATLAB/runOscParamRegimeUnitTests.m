% A script for running unittests for the getOscParamRegime function
%      
% Author: Wisam Reid
% Email: wisam@ccrma.stanford.edu
%
%% CLEAN AND CLEAR

close all
clear
clc

%% Utilities

% test folder location
test_folder_path = 'tests/Param_Regime/';

% In order to test several classes together 
% we use the Test Suite Package
import matlab.unittest.TestSuite

% % Check for available methods
% methods('TestSuite')

% Only work on failed tests?
failed_only = 0;

%% Testing getFP function

% grab all the tests in the folder
getFPTestSuite = TestSuite.fromFolder(test_folder_path);

getFP_All_Results = run(getFPTestSuite);

% display test results
disp(getFP_All_Results)

%% reRun Failed Tests 

if failed_only
    clc 
    failedTests = getFPTestSuite([getFP_All_Results.Failed]);
    run(failedTests)
end