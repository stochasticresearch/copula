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
            % make a hybrid network as follows w/ simulated data.
            %  A    B
            %   \ /  \ 
            %    C    D
            %    |
            %    E
            % All arrows point downwards.  A, B, and E will be continuous
            % nodes, C and D will be discrete nodes.  We will induce a 
            % Gaussian dependency between the nodes. of different
            % correlation amounts to simulate different dependencies
            % between the nodes.
            M = 1000;
            D = 5;
            
            % Generate samples from C1 (A,B,C) [Gaussian Copula]
            Rho_C1 = [1 .4 .2; .4 1 -.8; .2 -.8 1];
            Z = mvnrnd([0 0 0], Rho_C1, M);
            U_C1 = normcdf(Z,0,1);
            
            % Generate samples from C2 (B,D) [Clayton Copula]
            U_C2_1 = U_C1(:,2); c2_alpha = 2; p = rand(M,1);
            U_C2_2 = U_C2_1.*(p.^(-c2_alpha./(1+c2_alpha)) - 1 + U_C2_1.^c2_alpha).^(-1./c2_alpha);
            U_C2 = [U_C2_1 U_C2_2];
            
            % Generate samples from C3 (C,E) [Clayton Copula]
            U_C3_1 = U_C1(:,3); c3_alpha = 4; p = rand(M,1);
            U_C3_2 = U_C3_1.*(p.^(-c3_alpha./(1+c3_alpha)) - 1 + U_C3_1.^c3_alpha).^(-1./c3_alpha);
            U_C3 = [U_C3_1 U_C3_2];
            
            U = [U_C1 U_C2(:,2) U_C3(:,2)];
            
            X = [gaminv(U(:,1),2,1) ...
                   betainv(U(:,2),2,2) ...
                   unidinv(U(:,3),5) ...
                   unidinv(U(:,4),3) ...
                   norminv(U(:,5),0,1)];
                        
            nodes = {'A', 'B', 'C', 'D', 'E'};
            discreteNodes = {'C', 'D'};
            
            % make a dag that is correct (acyclic)
            dag = zeros(D,D);
            aa = 1; bb = 2; cc = 3; dd = 4; ee = 5;
            dag(aa,cc) = 1;
            dag(bb,[cc dd]) = 1;
            dag(cc, ee) = 1;
            hcbnObj = hcbn(testCase.bntPath, X, nodes, discreteNodes, dag);
            
            % ensure that we have no families for the root nodes
            A_family_expect = []; A_family_actual = hcbnObj.copulaFamilies{1};
            B_family_expect = []; B_family_actual = hcbnObj.copulaFamilies{2};
            testCase.verifyEqual(A_family_actual, A_family_expect);
            testCase.verifyEqual(B_family_actual, B_family_expect);
            
            % we compare the empirically calculated copulas to the known
            % form values
            D_family_actual = hcbnObj.copulaFamilies{4}.C;
            D_family_U = hcbnObj.copulaFamilies{4}.U;
            D_family_expect = zeros(size(D_family_actual));
            for ii=1:length(D_family_expect)
                for jj=1:length(D_family_expect)
                    u = [D_family_U{1}(ii,jj) D_family_U{2}(ii,jj)];
                    D_family_expect(ii,jj) = copulacdf('Clayton', u, c2_alpha);
                end
            end
            subplot(1,2,1); contour(D_family_U{1},D_family_U{2}, D_family_actual);
            grid on; title('Empirical C2')
            subplot(1,2,2); contour(D_family_U{1},D_family_U{2}, D_family_expect);
            grid on; title('Theoretical C2');
            pause;
            
            E_family_actual = hcbnObj.copulaFamilies{5}.C;
            E_family_U = hcbnObj.copulaFamilies{5}.U;
            E_family_expect = zeros(size(E_family_actual));
            for ii=1:length(E_family_expect)
                for jj=1:length(E_family_expect)
                    u = [E_family_U{1}(ii,jj) E_family_U{2}(ii,jj)];
                    E_family_expect(ii,jj) = copulacdf('Clayton', u, c3_alpha);
                end
            end
            subplot(1,2,1); contour(E_family_U{1},E_family_U{2}, E_family_actual);
            grid on; title('Empirical C3')
            subplot(1,2,2); contour(E_family_U{1},E_family_U{2}, E_family_expect);
            grid on; title('Theoretical C3');
            pause;

            % TODO: would it be useful to add an all continuous network
            % here?
        end
        
        function testCopulaRatioVal_allModel(testCase)
            % make a hybrid network as follows w/ simulated data.
            %  A    B
            %   \ /  \ 
            %    C    D
            %    |
            %    E
            % All arrows point downwards.  All nodes will be continuous to
            % ensure that we calculate the copula ratio properly (easy to
            % check w/ closed form values provided by matlab function
            % copulapdf this way :D )
            M = 1000;
            D = 5;
            
            % Generate samples from C1 (A,B,C) [Gaussian Copula]
            Rho_C1 = [1 .4 .2; .4 1 -.8; .2 -.8 1];
            Z = mvnrnd([0 0 0], Rho_C1, M);
            U_C1 = normcdf(Z,0,1);
            
            % Generate samples from C2 (B,D) [Clayton Copula]
            U_C2_1 = U_C1(:,2); c2_alpha = 2; p = rand(M,1);
            U_C2_2 = U_C2_1.*(p.^(-c2_alpha./(1+c2_alpha)) - 1 + U_C2_1.^c2_alpha).^(-1./c2_alpha);
            U_C2 = [U_C2_1 U_C2_2];
            
            % Generate samples from C3 (C,E) [Clayton Copula]
            U_C3_1 = U_C1(:,3); c3_alpha = 4; p = rand(M,1);
            U_C3_2 = U_C3_1.*(p.^(-c3_alpha./(1+c3_alpha)) - 1 + U_C3_1.^c3_alpha).^(-1./c3_alpha);
            U_C3 = [U_C3_1 U_C3_2];
            
            U = [U_C1 U_C2(:,2) U_C3(:,2)];
            
            X = [gaminv(U(:,1),2,1) ...
                   betainv(U(:,2),2,2) ...
                   expinv(U(:,3),5) ...
                   expinv(U(:,4),3) ...
                   norminv(U(:,5),0,1)];
                        
            nodes = {'A', 'B', 'C', 'D', 'E'};
            discreteNodes = {};
            
            % make a dag that is correct (not acyclic)
            dag = zeros(D,D);
            aa = 1; bb = 2; cc = 3; dd = 4; ee = 5;
            dag(aa,cc) = 1;
            dag(bb,[cc dd]) = 1;
            dag(cc, ee) = 1;
            hcbnObj = hcbn(testCase.bntPath, X, nodes, discreteNodes, dag);
            
            % the ratio values for nodes A,B should be 1 b/c they have no
            % parents (defined by Elidan 2010). (the value of u shouldn't
            % matter here)
            nodeA_Rc_u1_expect = 1; 
            nodeA_Rc_u1_actual = hcbnObj.copulaRatioVal(aa, [.1 .1], []);
            testCase.verifyEqual(nodeA_Rc_u1_actual, nodeA_Rc_u1_expect);
            
            nodeB_Rc_u1_expect = 1; 
            nodeB_Rc_u1_actual = hcbnObj.copulaRatioVal(bb, [.3 .3], []);
            testCase.verifyEqual(nodeB_Rc_u1_actual, nodeB_Rc_u1_expect);
            
            % Test node's D and E.  For this, we need to fit the families
            % to a Gaussian copula, b/c the hcbn code currently fits all
            % continuous nodes to the Gaussian copula even though they may
            % have originated from a different copula (Clayton in the case
            % above)
            % TODO: incorporate this into test also
            U_C1_test = [U_C1(:,3) U_C1(:,1) U_C1(:,2)];   % permute U_C1 to match child first, then parents
            Rho_C1_expect = copulafit('Gaussian', U_C1_test);
            Rho_C1_actual = hcbnObj.copulaFamilies{cc}.Rho;
            testCase.verifyEqual(Rho_C1_actual, Rho_C1_expect, 'AbsTol', 0.1);
            
            Rho_C1_1_expect = copulafit('Gaussian', U_C1_test(:,2:end));
            Rho_C1_1_actual = hcbnObj.copulaFamilies{cc}.Rho_parents;
            testCase.verifyEqual(Rho_C1_1_actual, Rho_C1_1_expect, 'AbsTol', 0.1);
            
            Rho_C2_expect = copulafit('Gaussian', U_C2);
            Rho_C2_actual = hcbnObj.copulaFamilies{dd}.Rho;
            testCase.verifyEqual(Rho_C2_actual, Rho_C2_expect, 'AbsTol', 0.1);
            
            Rho_C3_expect = copulafit('Gaussian', U_C3);
            Rho_C3_actual = hcbnObj.copulaFamilies{ee}.Rho;
            testCase.verifyEqual(Rho_C3_actual, Rho_C3_expect, 'AbsTol', 0.1);
            
            % Test Copula familes C1, C2, and C3
            for ii=randperm(M, 100)       % randomely test 100 of them
                x_c = [X(ii,cc) X(ii,aa) X(ii,bb)];
                u_c = [hcbnObj.empInfo{dd}.queryDistribution(X(ii,cc)) ...
                        hcbnObj.empInfo{aa}.queryDistribution(X(ii,aa)) ...
                        hcbnObj.empInfo{bb}.queryDistribution(X(ii,bb))];
                x_d = [X(ii,dd) X(ii,bb)];
                u_d = [hcbnObj.empInfo{dd}.queryDistribution(X(ii,dd)) ...
                        hcbnObj.empInfo{bb}.queryDistribution(X(ii,bb))];
                x_e = [X(ii,ee) X(ii,cc)];
                u_e = [hcbnObj.empInfo{ee}.queryDistribution(X(ii,ee)) ...
                        hcbnObj.empInfo{cc}.queryDistribution(X(ii,cc))];
                
                nodeC_Rc_actual = hcbnObj.copulaRatioVal(cc, u_c, x_c);
                nodeD_Rc_actual = hcbnObj.copulaRatioVal(dd, u_d, x_d);
                nodeE_Rc_actual = hcbnObj.copulaRatioVal(ee, u_e, x_e);
                
                % calculate expected copula ratio
                nodeC_Rc_numerator = copulapdf('Gaussian', u_c, hcbnObj.copulaFamilies{cc}.Rho);
                nodeD_Rc_numerator = copulapdf('Gaussian', u_d, hcbnObj.copulaFamilies{dd}.Rho);
                nodeE_Rc_numerator = copulapdf('Gaussian', u_e, hcbnObj.copulaFamilies{ee}.Rho);

                % calculate denominator manually
                nodeC_Rc_denominator = copulapdf('Gaussian', u_c(2:end), hcbnObj.copulaFamilies{cc}.Rho_parents);
                nodeD_Rc_denominator = hcbnObj.empInfo{bb}.queryDensity(x_d(2));
                nodeE_Rc_denominator = hcbnObj.empInfo{cc}.queryDensity(X(ii,cc));

                % compare 
                nodeC_Rc_expected = nodeC_Rc_numerator/nodeC_Rc_denominator;
                nodeD_Rc_expected = nodeD_Rc_numerator/nodeD_Rc_denominator;
                nodeE_Rc_expected = nodeE_Rc_numerator/nodeE_Rc_denominator;
                
                testCase.verifyEqual(nodeC_Rc_actual, nodeC_Rc_expected, 'AbsTol', 0.01);
                testCase.verifyEqual(nodeD_Rc_actual, nodeD_Rc_expected, 'AbsTol', 0.01);
                testCase.verifyEqual(nodeE_Rc_actual, nodeE_Rc_expected, 'AbsTol', 0.01);
            end
        end
        
        function testCopulaRatioVal_allEmpirical(testCase)
            % make a hybrid network as follows w/ simulated data.
            %  A    B
            %   \ /  \ 
            %    C    D
            %    |
            %    E
            % All arrows point downwards.  A, B, and E will be continuous
            % nodes, C and D will be discrete nodes.  We will induce a 
            % Gaussian dependency between the nodes. of different
            % correlation amounts to simulate different dependencies
            % between the nodes.
            M = 1000;
            D = 5;
            
            % Generate samples from C1 (A,B,C) [Gaussian Copula]
            Rho_C1 = [1 .4 .2; .4 1 -.8; .2 -.8 1];
            Z = mvnrnd([0 0 0], Rho_C1, M);
            U_C1 = normcdf(Z,0,1);
            
            % Generate samples from C2 (B,D) [Clayton Copula]
            U_C2_1 = U_C1(:,2); c2_alpha = 2; p = rand(M,1);
            U_C2_2 = U_C2_1.*(p.^(-c2_alpha./(1+c2_alpha)) - 1 + U_C2_1.^c2_alpha).^(-1./c2_alpha);
            U_C2 = [U_C2_1 U_C2_2];
            
            % Generate samples from C3 (C,E) [Clayton Copula]
            U_C3_1 = U_C1(:,3); c3_alpha = 4; p = rand(M,1);
            U_C3_2 = U_C3_1.*(p.^(-c3_alpha./(1+c3_alpha)) - 1 + U_C3_1.^c3_alpha).^(-1./c3_alpha);
            U_C3 = [U_C3_1 U_C3_2];
            
            U = [U_C1 U_C2(:,2) U_C3(:,2)];
            
            X = [gaminv(U(:,1),2,1) ...
                   betainv(U(:,2),2,2) ...
                   unidinv(U(:,3),5) ...
                   unidinv(U(:,4),3) ...
                   norminv(U(:,5),0,1)];
                        
            nodes = {'A', 'B', 'C', 'D', 'E'};
            discreteNodes = {'C', 'D'};
            
            % make a dag that is correct (acyclic)
            dag = zeros(D,D);
            aa = 1; bb = 2; cc = 3; dd = 4; ee = 5;
            dag(aa,cc) = 1;
            dag(bb,[cc dd]) = 1;
            dag(cc, ee) = 1;
            hcbnObj = hcbn(testCase.bntPath, X, nodes, discreteNodes, dag);
            
            % the ratio values for nodes A,B should be 1 b/c they have no
            % parents (defined by Elidan 2010). (the value of u shouldn't
            % matter here)
            nodeA_Rc_u1_expect = 1; 
            nodeA_Rc_u1_actual = hcbnObj.copulaRatioVal(aa, [.1 .1], []);
            testCase.verifyEqual(nodeA_Rc_u1_actual, nodeA_Rc_u1_expect);
            
            nodeB_Rc_u1_expect = 1; 
            nodeB_Rc_u1_actual = hcbnObj.copulaRatioVal(bb, [.3 .3], []);
            testCase.verifyEqual(nodeB_Rc_u1_actual, nodeB_Rc_u1_expect);
            
            % Test Copula familes C1, C2, and C3
            for ii=randperm(M, 100)       % randomely test 100 of them
                x_c = [X(ii,cc) X(ii,aa) X(ii,bb)];
                u_c = [hcbnObj.empInfo{dd}.queryDistribution(X(ii,cc)) ...
                        hcbnObj.empInfo{aa}.queryDistribution(X(ii,aa)) ...
                        hcbnObj.empInfo{bb}.queryDistribution(X(ii,bb))];
                x_d = [X(ii,dd) X(ii,bb)];
                u_d = [hcbnObj.empInfo{dd}.queryDistribution(X(ii,dd)) ...
                        hcbnObj.empInfo{bb}.queryDistribution(X(ii,bb))];
                x_e = [X(ii,ee) X(ii,cc)];
                u_e = [hcbnObj.empInfo{ee}.queryDistribution(X(ii,ee)) ...
                        hcbnObj.empInfo{cc}.queryDistribution(X(ii,cc))];
                
                nodeC_Rc_actual = hcbnObj.copulaRatioVal(cc, u_c, x_c);
                nodeD_Rc_actual = hcbnObj.copulaRatioVal(dd, u_d, x_d);
                nodeE_Rc_actual = hcbnObj.copulaRatioVal(ee, u_e, x_e);
                
                % calculate expected copula ratio
                [~,nodeC_Rc_numerator] = empcopula_val(hcbnObj.copulaFamilies{cc}.C, hcbnObj.copulaFamilies{cc}.c{end}, u_c);
                [~,nodeD_Rc_numerator] = empcopula_val(hcbnObj.copulaFamilies{dd}.C, hcbnObj.copulaFamilies{dd}.c{end}, u_d);
                [~,nodeE_Rc_numerator] = empcopula_val(hcbnObj.copulaFamilies{ee}.C, hcbnObj.copulaFamilies{ee}.c{end}, u_e);

                % calculate denominator manually
                [~,nodeC_Rc_denominator] = empcopula_val(hcbnObj.copulaFamilies{cc}.C_parents, hcbnObj.copulaFamilies{cc}.c_parents, u_c(2:end));
                nodeD_Rc_denominator = hcbnObj.empInfo{dd}.queryDensity(x_d(2));
                nodeE_Rc_denominator = hcbnObj.empInfo{dd}.queryDensity(x_e(2));

                % compare 
                nodeC_Rc_expected = nodeC_Rc_numerator/nodeC_Rc_denominator;
                nodeD_Rc_expected = nodeD_Rc_numerator/nodeD_Rc_denominator;
                nodeE_Rc_expected = nodeE_Rc_numerator/nodeE_Rc_denominator;
                
                testCase.verifyEqual(nodeC_Rc_actual, nodeC_Rc_expected, 'AbsTol', 0.01);
                testCase.verifyEqual(nodeD_Rc_actual, nodeD_Rc_expected, 'AbsTol', 0.01);
                testCase.verifyEqual(nodeE_Rc_actual, nodeE_Rc_expected, 'AbsTol', 0.01);
            end
        end
            
        function testCopulall(testCase)
            
        end
        
        function testHcbnLogLikelihood(testCase)
            
        end
        
        function testLearnStruct_hc(testCase)
            
        end
    end
end

