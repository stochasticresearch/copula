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
            
            for ii=1:D
                subplot(1,3,1);
                histogram(X(:,ii))
                title('Histogram')
                subplot(1,3,2);
                plot(hcbnObj.empInfo{ii}.domain, hcbnObj.empInfo{ii}.density)
                grid on;
                title('Empirical Density Function')
                subplot(1,3,3)
                plot(hcbnObj.empInfo{ii}.domain, hcbnObj.empInfo{ii}.distribution)
                grid on;
                title('Empirical Distribution Function')
                pause;
            end
            
        end
        
        function testAcyclicCheck(testCase)
            
        end
        
        function testGetParents(testCase)
            
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

