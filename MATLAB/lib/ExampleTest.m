classdef ExampleTest < matlab.unittest.Testcase
    % ExampleTest 
    %   Template for Matlab Unittests
    
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
        function testName(testCase)
            actSolution = functionToTest(arg1, arg2);
            expSolution = 1;
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',sqrt(eps));
        end
    end
    
end

% % Run test file
% testCase1 = ExampleTest;
% result1 = run(testCase1);

% % Run a single test
% result2 = run(testCase1,'testName');

% % In order to test several classes together 
% % we use the Test Suite Package
% import matlab.unittest.TestSuite
% methods('TestSuite')

% % Methods for class TestSuite:
% % 
% % details      display      selectIf     
% % disp         run          
% % 
% % Static methods:
% % 
% % fromClass    fromMethod   
% % fromFile     fromName     
% % fromFolder   fromPackage 

% % From a folder 
% run(TestSuite.fromFolder('tests'));

% % 
% % From a file 
% run(TestSuite.fromFile('ExampleTest.m'));

% % combine them
% suite1 = TestSuite.fromFolder('tests');
% suite2 = TestSuite.fromFile('ExampleTest.m');
% totalSuite = [suite1, suite2];
% totalResult = run(totalSuite);
% disp(TotalResult)

% % reRun Failed Tests 
% failedTests = totalSuite([totalResult.Failed]);

