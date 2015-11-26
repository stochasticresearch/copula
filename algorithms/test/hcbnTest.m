classdef hcbnTest < matlab.unittest.TestCase
    %hcbnTest Test class for the HCBN definition (hcbn.m)
    
    properties
        bntPath;
    end
 
    methods(TestMethodSetup)
        function defineBntPath(testCase)
            % comment
            testCase.bntPath = '../../../bnt';
        end
    end
    
    methods(TestMethodTeardown)
        function closeFigure(testCase)
            close all;
        end
    end
    
    methods (Test)
        
        function testConstructor(testCase)
            % This test verifies that the correct node indices for the
            % discrete nodes are extracted from the provided constructor
            % information 
            
            % make bogus data
            D = 9;
            M = 100;
            X = randn(M,D);
            
            nodes = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'};
            discreteNodes = {'C', 'H', 'I'};
            X(:,3) = unidrnd(5,M,1);
            X(:,8) = unidrnd(5,M,1);
            X(:,9) = unidrnd(5,M,1);
            
            hcbnObj = hcbn(testCase.bntPath, X, nodes, discreteNodes);
            
            expectedDiscNodeIdxs = [3,8,9];
            testCase.verifyEqual(hcbnObj.discNodeIdxs, expectedDiscNodeIdxs);
        end
        
        function testCalcEmpInfo(testCase)
            % This test verifies that the empirical information calculated
            % is correct
            % make bogus data
            D = 9;
            M = 100;
            X = randn(M,D);
            
            nodes = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'};
            discreteNodes = {'C', 'H', 'I'};
            X(:,3) = unidrnd(5,M,1);
            X(:,8) = unidrnd(5,M,1);
            X(:,9) = unidrnd(5,M,1);
            
            hcbnObj = hcbn(testCase.bntPath, X, nodes, discreteNodes);
            
%             for ii=1:D
%                 subplot(1,3,1);
%                 histogram(X(:,ii))
%                 title('Histogram')
%                 subplot(1,3,2);
%                 plot(hcbnObj.empInfo{ii}.domain, hcbnObj.empInfo{ii}.density)
%                 grid on;
%                 title('Empirical Density Function')
%                 subplot(1,3,3)
%                 plot(hcbnObj.empInfo{ii}.domain, hcbnObj.empInfo{ii}.distribution)
%                 grid on;
%                 title('Empirical Distribution Function')
%                 pause;
%             end
            
        end
        
        function testGetParents(testCase)
            % make bogus data
            D = 4;
            M = 100;
            X = randn(M,D);
            
            nodes = {'C', 'S', 'R', 'W'};
            discreteNodes = {};
            
            % make a dag that is correct (not acyclic)
            dag = zeros(D,D);
            C = 1; S = 2; R = 3; W = 4;
            dag(C,[R S]) = 1;
            dag(R,W) = 1;
            dag(S,W)=1;
            hcbnObj = hcbn(testCase.bntPath, X, nodes, discreteNodes, dag);
            
            % query the parents
            [C_parent, C_parent_name] = hcbnObj.getParents(C);
            [S_parent, S_parent_name] = hcbnObj.getParents(S);
            [R_parent, R_parent_name] = hcbnObj.getParents(R);
            [W_parent, W_parent_name] = hcbnObj.getParents(W);
            C_parent_expect = zeros(1,0);    C_parent_name_expect = cell(1,0);
            S_parent_expect = C;     S_parent_name_expect = {'C'};
            R_parent_expect = C;     R_parent_name_expect = {'C'};
            W_parent_expect = [S R]; W_parent_name_expect = {'S', 'R'};
            
            testCase.verifyEqual(C_parent, C_parent_expect);
            testCase.verifyEqual(S_parent, S_parent_expect);
            testCase.verifyEqual(R_parent, R_parent_expect);
            testCase.verifyEqual(W_parent, W_parent_expect);
            testCase.verifyEqual(C_parent_name, C_parent_name_expect);
            testCase.verifyEqual(S_parent_name, S_parent_name_expect);
            testCase.verifyEqual(R_parent_name, R_parent_name_expect);
            testCase.verifyEqual(W_parent_name, W_parent_name_expect);
        end
        
        function testEstFamilyCopula(testCase)
            
        end
        
        function testHcbnLogLikelihood(testCase)
            
        end
        
        function testCopulall(testCase)
            
        end
        
        function testCopulaRatioVal(testCase)
            
        end
        
        function testLearnStruct_hc(testCase)
            
        end
    end
    
end

