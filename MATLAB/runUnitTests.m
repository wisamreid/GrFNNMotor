% A script for running unittests
% 
% Fixed Point Location and Stability Tests from Kim & Large 2015 figures
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
test_folder_path = 'tests/';

% In order to test several classes together 
% we use the Test Suite Package
import matlab.unittest.TestSuite

% % Check for available methods
% methods('TestSuite')

% Load tests
F3 = TestSuite.fromFile(strcat(test_folder_path,'getFP_F3_Test.m'));
F4A = TestSuite.fromFile(strcat(test_folder_path,'getFP_F4A_Test.m'));
F4B = TestSuite.fromFile(strcat(test_folder_path,'getFP_F4B_Test.m'));
F4C = TestSuite.fromFile(strcat(test_folder_path,'getFP_F4C_Test.m'));

%% Testing getFP function

% getFPTestSuite = [F3, F4A, F4B, F4C]; % select figures from Kim & Large 2015
% getFP_All_Results = run(getFPTestSuite);

% % display test results
% disp(getFP_All_Results)

% % reRun Failed Tests 
% failedTests = totalSuite([totalResult.Failed]);

%%

getFPTestSuite = [F4C];
getFP_All_Results = run(getFPTestSuite);


%% Testing getSS function

