%**************************************************************************
%* 
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>
%*
%* This program is free software: you can redistribute it and/or modify
%* it under the terms of the GNU General Public License as published by
%* the Free Software Foundation, either version 3 of the License, or
%* (at your option) any later version.
%*
%* This program is distributed in the hope that it will be useful,
%* but WITHOUT ANY WARRANTY; without even the implied warranty of
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%* GNU General Public License for more details.
%*
%* You should have received a copy of the GNU General Public License
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
            
            % TODO: put in test case for node C
            
            % we compare the empirically calculated copulas to the known
            % form values
            D_family_actual = hcbnObj.copulaFamilies{4}.c;
            [D_family_U1,D_family_U2] = ndgrid(linspace(0,1,hcbnObj.K));
            D_family_expect = copulapdf('Clayton', [D_family_U1(:) D_family_U2(:)],c2_alpha);
            D_family_expect = reshape(D_family_expect,hcbnObj.K,hcbnObj.K);
            h1 = subplot(1,2,1); surf(D_family_U1,D_family_U2, D_family_actual);
            grid on; title('Empirical C2')
            h2 = subplot(1,2,2); surf(D_family_U1,D_family_U2, D_family_expect);
            grid on; title('Theoretical C2');
            linkprop([h1,h2],{'CameraPosition','CameraUpVector'}); rotate3d on; pause;
            
            E_family_actual = hcbnObj.copulaFamilies{5}.c;
            [E_family_U1,E_family_U2] = ndgrid(linspace(0,1,hcbnObj.K));
            E_family_expect = zeros(size(E_family_actual));
            for ii=1:length(E_family_expect)
                for jj=1:length(E_family_expect)
                    u = [E_family_U1(ii,jj) E_family_U2(ii,jj)];
                    E_family_expect(ii,jj) = copulapdf('Clayton', u, c3_alpha);
                end
            end
            h1 = subplot(1,2,1); surf(E_family_U1,E_family_U2, E_family_actual);
            grid on; title('Empirical C3')
            h2 = subplot(1,2,2); surf(E_family_U1,E_family_U2, E_family_expect);
            grid on; title('Theoretical C3');
            linkprop([h1,h2],{'CameraPosition','CameraUpVector'}); rotate3d on; pause;
            
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
            nodeA_Rc_u1_actual = hcbnObj.copulaRatioVal(aa, [.1 .1]);
            testCase.verifyEqual(nodeA_Rc_u1_actual, nodeA_Rc_u1_expect);
            
            nodeB_Rc_u1_expect = 1; 
            nodeB_Rc_u1_actual = hcbnObj.copulaRatioVal(bb, [.3 .3]);
            testCase.verifyEqual(nodeB_Rc_u1_actual, nodeB_Rc_u1_expect);
            
            % Test node's C, D and E.  For this, we need to fit the families
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
            
            % Test Copula familes for nodes C, D, and E
            for ii=randperm(M, 100)       % randomely test 100 of them
                u_c = [hcbnObj.empInfo{dd}.cdf(X(ii,cc)) ...
                        hcbnObj.empInfo{aa}.cdf(X(ii,aa)) ...
                        hcbnObj.empInfo{bb}.cdf(X(ii,bb))];
                u_d = [hcbnObj.empInfo{dd}.cdf(X(ii,dd)) ...
                        hcbnObj.empInfo{bb}.cdf(X(ii,bb))];
                u_e = [hcbnObj.empInfo{ee}.cdf(X(ii,ee)) ...
                        hcbnObj.empInfo{cc}.cdf(X(ii,cc))];
                
                nodeC_Rc_actual = hcbnObj.copulaRatioVal(cc, u_c);
                nodeD_Rc_actual = hcbnObj.copulaRatioVal(dd, u_d);
                nodeE_Rc_actual = hcbnObj.copulaRatioVal(ee, u_e);
                
                % calculate expected copula ratio
                nodeC_Rc_numerator = copulapdf('Gaussian', u_c, hcbnObj.copulaFamilies{cc}.Rho);
                nodeD_Rc_numerator = copulapdf('Gaussian', u_d, hcbnObj.copulaFamilies{dd}.Rho);
                nodeE_Rc_numerator = copulapdf('Gaussian', u_e, hcbnObj.copulaFamilies{ee}.Rho);

                % calculate denominator manually
                nodeC_Rc_denominator = copulapdf('Gaussian', u_c(2:end), hcbnObj.copulaFamilies{cc}.Rho_parents);
                nodeD_Rc_denominator = 1;
                nodeE_Rc_denominator = 1;

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
            nodeA_Rc_u1_actual = hcbnObj.copulaRatioVal(aa, [.1 .1]);
            testCase.verifyEqual(nodeA_Rc_u1_actual, nodeA_Rc_u1_expect);
            
            nodeB_Rc_u1_expect = 1; 
            nodeB_Rc_u1_actual = hcbnObj.copulaRatioVal(bb, [.3 .3]);
            testCase.verifyEqual(nodeB_Rc_u1_actual, nodeB_Rc_u1_expect);
            
            % Test Copula familes for nodes C, D, and E
            for ii=randperm(M, 100)       % randomely test 100 of them
                u_c = [hcbnObj.empInfo{dd}.cdf(X(ii,cc)) ...
                        hcbnObj.empInfo{aa}.cdf(X(ii,aa)) ...
                        hcbnObj.empInfo{bb}.cdf(X(ii,bb))];
                u_d = [hcbnObj.empInfo{dd}.cdf(X(ii,dd)) ...
                        hcbnObj.empInfo{bb}.cdf(X(ii,bb))];
                u_e = [hcbnObj.empInfo{ee}.cdf(X(ii,ee)) ...
                        hcbnObj.empInfo{cc}.cdf(X(ii,cc))];
                
                nodeC_Rc_actual = hcbnObj.copulaRatioVal(cc, u_c);
                nodeD_Rc_actual = hcbnObj.copulaRatioVal(dd, u_d);
                nodeE_Rc_actual = hcbnObj.copulaRatioVal(ee, u_e);
                
                % calculate expected copula ratio
                nodeC_Rc_numerator = empcopulaval(hcbnObj.copulaFamilies{cc}.c, u_c);
                nodeD_Rc_numerator = empcopulaval(hcbnObj.copulaFamilies{dd}.c, u_d);
                nodeE_Rc_numerator = empcopulaval(hcbnObj.copulaFamilies{ee}.c, u_e);

                % calculate denominator manually
                nodeC_Rc_denominator = empcopulaval(hcbnObj.copulaFamilies{cc}.c_parents, u_c(2:end));
                nodeD_Rc_denominator = 1;
                nodeE_Rc_denominator = 1;

                % compare 
                nodeC_Rc_expected = nodeC_Rc_numerator/nodeC_Rc_denominator;
                nodeD_Rc_expected = nodeD_Rc_numerator/nodeD_Rc_denominator;
                nodeE_Rc_expected = nodeE_Rc_numerator/nodeE_Rc_denominator;
                
                testCase.verifyEqual(nodeC_Rc_actual, nodeC_Rc_expected, 'AbsTol', 0.01);
                testCase.verifyEqual(nodeD_Rc_actual, nodeD_Rc_expected, 'AbsTol', 0.01);
                testCase.verifyEqual(nodeE_Rc_actual, nodeE_Rc_expected, 'AbsTol', 0.01);
            end
        end
        
        function testQueryDensity(testCase)
            M = 1000;
            X = normrnd(0,1,M,1);
            
            [F,x] = ecdf(X);
            F = F(2:end);
            x = x(2:end);
            f = ksdensity(X,x);
            isdiscrete = 0;
            empInfoObj = rvEmpiricalInfo(x, f, F, isdiscrete);
            
            for m=1:M
                f_expect = normpdf(X(m));
                f_actual = empInfoObj.pdf(X(m));
                
                F_expect = normcdf(X(m));
                F_actual = empInfoObj.cdf(X(m));
                
                testCase.verifyEqual(f_actual, f_expect, 'AbsTol', 0.1);
                testCase.verifyEqual(F_actual, F_expect, 'AbsTol', 0.1);
            end
        end

        function testCopulall_allContinuous(testCase)
            % make a hybrid network as follows w/ simulated data.
            %  A    B
            %   \ /  
            %    C 
            % All arrows point downwards.  All nodes will be continuous to
            % ensure that we calculate the copula ratio properly (easy to
            % check w/ closed form values provided by matlab function
            % copulapdf this way :D )
            M = 1000;
            D = 3;
            
            % Generate samples from C1 (A,B,C) [Gaussian Copula]
            Rho_C1 = [1 .4 .2; .4 1 -.8; .2 -.8 1];
            Z = mvnrnd([0 0 0], Rho_C1, M);
            U = normcdf(Z,0,1);
            
            X = [norminv(U(:,1),0,1) ...
                 norminv(U(:,2),0,1) ...
                 norminv(U(:,3),0,1)];
                        
            nodes = {'A', 'B', 'C'};
            discreteNodes = {};
            
            % make a dag that is correct (not acyclic)
            dag = zeros(D,D);
            aa = 1; bb = 2; cc = 3;
            dag(aa,cc) = 1;
            dag(bb,cc) = 1;
            hcbnObj = hcbn(testCase.bntPath, X, nodes, discreteNodes, dag);
            
            % for this test, we choose to use the exact same approximated
            % rho matrices to minimize propagation of error.  The above
            % tests compare how the fit of the hcbn code versus our fit
            % here match up.
            Rho_C1_est = hcbnObj.copulaFamilies{cc}.Rho;
            Rho_C1_parents_est = hcbnObj.copulaFamilies{cc}.Rho_parents;
            
            % first test to see if the copula ratio's are within the range
            % of expected values
            Rc_expect_nodeA = 1;
            Rc_expect_nodeB = 1;
            [~,Rc_val_actual_nodeA_vec] = hcbnObj.copulall(aa,X);
            [~,Rc_val_actual_nodeB_vec] = hcbnObj.copulall(bb,X);
            [~,Rc_val_actual_nodeC_vec] = hcbnObj.copulall(cc,X);
            Rc_expect_nodeC_vec = zeros(1,M);
            for m=1:M
                numVal = mvnpdf([X(m,3) X(m,1) X(m,2)], [0 0 0], Rho_C1_est);
                denVal = mvnpdf([X(m,1) X(m,2)], [0 0], Rho_C1_parents_est)*normpdf(X(m,3));
                Rc_expect_nodeC_vec(m) = numVal/denVal;
            end
            Rc_actual_nodeA = mean(Rc_val_actual_nodeA_vec);
            Rc_actual_nodeB = mean(Rc_val_actual_nodeB_vec);
            testCase.verifyEqual(Rc_actual_nodeA, Rc_expect_nodeA);
            testCase.verifyEqual(Rc_actual_nodeB, Rc_expect_nodeB);
            % NOTE: we expect some errors here, the conversion X->U->X->U
            % will have some effects on the error
%             idxVec = 1:M;
%             [idxVec' Rc_expect_nodeC_vec' Rc_val_actual_nodeC_vec' X]
            testCase.verifyEqual(Rc_val_actual_nodeC_vec, Rc_expect_nodeC_vec, 'RelTol', 0.25);
            
            % manually calculate copula likelihood value for nodes A,B,C
            ll_val_expect_nodeA = 0;
            ll_val_expect_nodeB = 0;
            ll_val_expect_nodeC = 0;
            for m=1:M
                ll_val_expect_nodeA = ll_val_expect_nodeA + log(normpdf(X(m,1)));
                ll_val_expect_nodeB = ll_val_expect_nodeB + log(normpdf(X(m,2)));
                % the order here matters b/c we defined Rho_C1 for this
                % order (X3,X1,X2)
                Rc = mvnpdf([X(m,3) X(m,1) X(m,2)], [0 0 0], Rho_C1_est)/( mvnpdf([X(m,1) X(m,2)], [0 0], Rho_C1_parents_est)*normpdf(X(m,3)) );
                ll_val_expect_nodeC = ll_val_expect_nodeC + log(normpdf(X(m,3))) + log(Rc);
            end
            ll_val_actual_nodeA = hcbnObj.copulall(aa,X);
            ll_val_actual_nodeB = hcbnObj.copulall(bb,X);
            ll_val_actual_nodeC = hcbnObj.copulall(cc,X);
            
            testCase.verifyEqual(ll_val_actual_nodeA, ll_val_expect_nodeA, 'RelTol', 0.25);
            testCase.verifyEqual(ll_val_actual_nodeB, ll_val_expect_nodeB, 'RelTol', 0.25);
            testCase.verifyEqual(ll_val_actual_nodeC, ll_val_expect_nodeC, 'RelTol', 0.25);
            % some experiments
%             for m=randperm(M,10)
%                 e1e = copulapdf('Gaussian', [U(m,3) U(m,1) U(m,2)], Rho_C1_est)/copulapdf('Gaussian', [U(m,1) U(m,2)], Rho_C1_parents_est);
%                 e2e = mvnpdf([X(m,3) X(m,1) X(m,2)], [0 0 0], Rho_C1_est)/( mvnpdf([X(m,1) X(m,2)], [0 0], Rho_C1_parents_est)*normpdf(X(m,3)) );
%                 fprintf('e1e=%d \t e2e=%d\n', e1e, e2e);
%             end
        end
            
        function testCopulall_allEmpirical(testCase)
            % make a hybrid network as follows w/ simulated data.
            %  A    B
            %   \ /  
            %    C 
            % All arrows point downwards.  All nodes are discrete;
            M = 1000;
            D = 3;
            
            % Generate samples from C1 (A,B,C) [Gaussian Copula]
            Rho_C1 = [1 .4 .2; .4 1 -.8; .2 -.8 1];
            Z = mvnrnd([0 0 0], Rho_C1, M);
            U = normcdf(Z,0,1);
            
            X = [unidinv(U(:,1),3) ...
                 unidinv(U(:,2),4) ...
                 unidinv(U(:,3),5)];
                        
            nodes = {'A', 'B', 'C'};
            discreteNodes = {'A', 'B', 'C'};
            
            % make a dag that is correct (not acyclic)
            dag = zeros(D,D);
            aa = 1; bb = 2; cc = 3;
            dag(aa,cc) = 1;
            dag(bb,cc) = 1;
            hcbnObj = hcbn(testCase.bntPath, X, nodes, discreteNodes, dag);
            
            % first test to see if the copula ratio's are within the range
            % of expected values
            Rc_expect_nodeA = 1;
            Rc_expect_nodeB = 1;
            [~,Rc_actual_nodeA_vec] = hcbnObj.copulall(aa,X);
            [~,Rc_actual_nodeB_vec] = hcbnObj.copulall(bb,X);
            Rc_actual_nodeA = mean(Rc_actual_nodeA_vec);
            Rc_actual_nodeB = mean(Rc_actual_nodeB_vec);
            testCase.verifyEqual(Rc_actual_nodeA, Rc_expect_nodeA);
            testCase.verifyEqual(Rc_actual_nodeB, Rc_expect_nodeB);
            
            % calculate Node C Rc and compare against expected Rc
            Rc_expect_nodeC_vec = zeros(1,M);
            Rc_expect_nodeC_num_vec = zeros(1,M);
            Rc_expect_nodeC_den_vec = zeros(1,M);
            [f_CAB, domain_CAB] = hist_discrete([X(:,3) X(:,1) X(:,2)]);
            [f_AB, domain_AB] = hist_discrete(X(:,1:2));
            [f_A, domain_A] = hist_discrete(X(:,1));
            [f_B, domain_B] = hist_discrete(X(:,2));
            [f_C, domain_C] = hist_discrete(X(:,3));
            for m=1:M
                % find which idx we need to query
                yy = [X(m,3) X(m,1) X(m,2)];
                CABidx = find(yy(1)==domain_CAB(1,:) & yy(2)==domain_CAB(2,:) & yy(3)==domain_CAB(3,:));
                
                yy = [X(m,1) X(m,2)];
                ABidx = find(yy(1)==domain_AB(1,:) & yy(2)==domain_AB(2,:));
                yy = [X(m,1)];
                Aidx = find(yy(1)==domain_A(1,:));
                yy = [X(m,2)];
                Bidx = find(yy(1)==domain_B(1,:));
                yy = [X(m,3)];
                Cidx = find(yy(1)==domain_C(1,:));
                
                numVal = f_CAB(CABidx);
                denVal1 = f_AB(ABidx);
                denVal2 = f_C(Cidx);
                
                Rc_expect_nodeC_vec(m) = numVal/(denVal1*denVal2);
                Rc_expect_nodeC_num_vec(m) = f_CAB(CABidx)/( f_A(Aidx)*f_B(Bidx)*f_C(Cidx) );
                Rc_expect_nodeC_den_vec(m) = f_AB(ABidx)/( f_A(Aidx)*f_B(Bidx) );
            end
            
            [~,Rc_actual_nodeC_vec, Rc_actual_nodeC_num_vec, Rc_actual_nodeC_den_vec] = hcbnObj.copulall(cc,X);
            
            numPtsToPlot = 50;
            subplot(3,1,1);
            plot(1:numPtsToPlot, Rc_expect_nodeC_vec(1:numPtsToPlot), 'b*', 1:numPtsToPlot, Rc_actual_nodeC_vec(1:numPtsToPlot), 'r*'); 
            title('Rc'); grid on; legend('Expected', 'Actual');
            
            subplot(3,1,2);
            plot(1:numPtsToPlot, Rc_expect_nodeC_num_vec(1:numPtsToPlot), 'b*', 1:numPtsToPlot, Rc_actual_nodeC_num_vec(1:numPtsToPlot), 'r*'); 
            title('c(A,B,C)'); grid on; legend('Expected', 'Actual');
            
            subplot(3,1,3);
            plot(1:numPtsToPlot, Rc_expect_nodeC_den_vec(1:numPtsToPlot), 'b*', 1:numPtsToPlot, Rc_actual_nodeC_den_vec(1:numPtsToPlot), 'r*'); 
            title('c(A,B)'); grid on; legend('Expected', 'Actual');
            pause;
            
        end
        
        function testLearnStruct_hc(testCase)
            
        end
    end
end

