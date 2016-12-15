% A script for running unittests for the getFP function
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
test_folder_path = 'tests/Fixed_Points/';

% In order to test several classes together 
% we use the Test Suite Package
import matlab.unittest.TestSuite

% % Check for available methods
% methods('TestSuite')

% Only work on failed tests?
failed_only = 1;

%% Load tests

% for loading tests from individual figures (See Kim & Large 2015)
F1 = TestSuite.fromFile(strcat(test_folder_path,'getFP_F1_Test.m'));
F3 = TestSuite.fromFile(strcat(test_folder_path,'getFP_F3_Test.m'));
F4A = TestSuite.fromFile(strcat(test_folder_path,'getFP_F4A_Test.m'));
F4B = TestSuite.fromFile(strcat(test_folder_path,'getFP_F4B_Test.m'));
F4C = TestSuite.fromFile(strcat(test_folder_path,'getFP_F4C_Test.m'));
F4D = TestSuite.fromFile(strcat(test_folder_path,'getFP_F4D_Test.m'));
F5 = TestSuite.fromFile(strcat(test_folder_path,'getFP_F5_Test.m'));
F6 = TestSuite.fromFile(strcat(test_folder_path,'getFP_F6_Test.m'));
F7 = TestSuite.fromFile(strcat(test_folder_path,'getFP_F7_Test.m'));
F8 = TestSuite.fromFile(strcat(test_folder_path,'getFP_F8_Test.m'));
F9 = TestSuite.fromFile(strcat(test_folder_path,'getFP_F9_Test.m'));
F10 = TestSuite.fromFile(strcat(test_folder_path,'getFP_F10_Test.m'));
F11 = TestSuite.fromFile(strcat(test_folder_path,'getFP_F11_Test.m'));

%% Testing getFP function

% grab all the tests in the folder
getFPTestSuite = TestSuite.fromFolder(test_folder_path);

% % run tests for a single figure
% getFPTestSuite = [F3];

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