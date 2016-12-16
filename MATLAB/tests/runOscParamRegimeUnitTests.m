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
test_folder_path = 'Param_Regime/';

% In order to test several classes together 
% we use the Test Suite Package
import matlab.unittest.TestSuite

% % Check for available methods
% methods('TestSuite')

% Only work on failed tests?
failed_only = 1;

%% Testing getFP function

% grab all the tests in the folder
getFPTestSuite = TestSuite.fromFolder(test_folder_path);

getFP_All_Results = run(getFPTestSuite);

% display test results
disp(getFP_All_Results)

%% reRun Failed Tests 

if failed_only
    failedTests = getFPTestSuite([getFP_All_Results.Failed]);
    numFailed = size(failedTests);
    if numFailed(2) ~= 0
        clc 
        remaining_fails = run(failedTests);
        % display remaining fails
        disp(remaining_fails)
    end
end