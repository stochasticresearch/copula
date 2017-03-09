function out = pwcCopula_D_dim_sparse;
disp('**********************************************************************');
disp('**********************************************************************');

disp('This program is copyright the original authors and the University of Bonn,');
disp('and is published here under the GNU General Public   ');
disp('License. (See http://www.fsf.org/licenses/licenses.html)');


clear all
format compact;


D=input('Dimension ');

K=input('K, number of subsets of [0,1]}  ');
n_whole_sample = input('n_whole_sample');
%n_by_K = input('floor(n_whole_sample/K) ');
n_by_K = floor(n_whole_sample/K);
n = n_by_K * K              %sample size = floor(n_whole_sample/K)*K

%************************************************************************
%***********************   Sample:    ***********************************


% [Z,n] = IP_Daten(n);        % Returns a sample Z(1:D, 1:n)
Z = Stichprobe_Weibull_dim_2(3,1,n);

n

%********************   figure properties:   |   ************************%
%*****************************************   |   *************************%
%*****************************************   V   *************************%
  
bilder = input('Figures: No(0), some(1), many (2), all (D)  ');

PosFigure_1 = zeros(4,4,4);
PosFigure_2 = zeros(4,4,4);
PosFigure_3 = zeros(4,4,4);
PosFigure_4 = zeros(4,4,4);


delta_Pos = -20;
for d=1:D
   
% Bereich fuer die Figuren:
   maxZ = max(Z(d,:))*1.1;
   log_u = floor(log10(maxZ));
   normiert = maxZ/ 10^log_u;
   if normiert < 3 
       Achse(d) = ceil(10*normiert) * 10^(log_u-1);
   else
       Achse(d) = ceil(normiert) * 10^log_u;
   end
    for ds=d+1:D
        
        delta_Pos = delta_Pos+20;

        PosFigure_1(:,d,ds) = [10+delta_Pos 580-delta_Pos 480 360];
        PosFigure_2(:,d,ds) = [10+delta_Pos 80-delta_Pos 380 360];
        PosFigure_3(:,d,ds) = [550+delta_Pos 80-delta_Pos 380 360];
        PosFigure_4(:,d,ds) = [650+delta_Pos 580-delta_Pos 480 360];
    end
end
%*********************************************************************%
%*********************************************************************%
%*********************************************************************%

    nbr=n;

for d=1:bilder
    for ds=d+1:D
        figure ('name', 'Sample Z', 'Position', PosFigure_1(:,d,ds))
        scatter(Z(d,1:nbr),Z(ds,1:nbr),2.5,'*')
        axis([0 Achse(d) 0 Achse(ds)])
        xlabel(d); 
        ylabel(ds)
    end
end

%*********************************************************************%
%Empirical copula data (sample frequencies), and marginal distributions:

[Z_sorted, nU, ~] = empVF_v3(n,D,Z); 

% Z_sorted1=Z_sorted(1,1:5)
% nU1=nU(1,1:5) 
% F_emp1=F_emp(1,1:5)

%********************************************************************
%********************************************************************
%copula (marginal) densities:

delta = 1/K
j=zeros(1,D);
f2=zeros(K,K);
summand = 1*realpow(K,D)/n;

switch D
    case 4
        f4 = sparse(K*K+K, K*K+K); % Access f4(i,j,k,l) (not sparse) as follows:
                                   % f4(i*K + j, k*K + l) in the sparse
                                   % matrix!
        f3 = sparse(K*K+K,K);      % Access f3(i,j,k) (not sparse) as follows:
                                   % f3(i*K + j, l) in the sparse
                                   % matrix!
        for i=1:n
            for d=1:D
                j(d) = ceil( (nU(d,i)-0.00001)/n_by_K );   % rounding errors omitted
            end;
            f4(j(1)*K+j(2),j(3)*K+j(4)) = f4(j(1)*K+j(2),j(3)*K+j(4)) + summand;
            f3(j(1)*K+j(2),j(3)) = f3(j(1)*K+j(2),j(3)) + summand;
            f2(j(1),j(2)) = f2(j(1),j(2)) + summand;
        end;
        % max_f4 = max(max(max(max(f4))))
        % sum_f2 = sum(sum(f2)) / realpow(K,D)
        % sum_f3 = sum(sum(sum(f3)))/ realpow(K,D)
        % sum_f4 = sum(sum(sum(sum(f4))))/ realpow(K,D)

    case 3
        f3 = sparse(K*K+K,K);      % Access f3(i,j,k) (not sparse) as follows:
                                   % f3(i*K + j, l) in the sparse
                                   % matrix!
        for i=1:n
            for d=1:D
                j(d) = ceil( (nU(d,i)-0.00001)/n_by_K );   % rounding errors omitted
            end;
            f3(j(1)*K+j(2),j(3)) = f3(j(1)*K+j(2),j(3)) + summand;
            f2(j(1),j(2)) = f2(j(1),j(2)) + summand;
        end;
        % max_f3 = max(max(max(f3)))
        % sum_f2 = sum(sum(f2))/ realpow(K,D)
        % sum_f3 = sum(sum(sum(f3)))/ realpow(K,D)

    case 2
        for i=1:n
            for d=1:D
                j(d) = ceil( (nU(d,i)-0.00001)/n_by_K );   % rounding errors omitted
            end;
            f2(j(1),j(2)) = f2(j(1),j(2)) + summand;
        end;
        % max_f2 = max(max(f2))
        % sum_f2 = sum(sum(f2))/ realpow(K,D)

end; %switch


%*********************************************************************
if bilder>0
        U=nU/n;
end;
for d=1:bilder
    for ds=d+1:D
        %plot   F_emp_ds ( Z(ds,:) )  against F_emp_d ( Z(d,:) )
        figure ('name','U(d,:), U(ds,:)',...
        'Position',PosFigure_2(:,d,ds))
        scatter(U(d,1:nbr),U(ds,1:nbr),2.5,'*')
        axis([-.1 1.1 -.1 1.1])
        xlabel(d)
        ylabel(ds)
    end
end

%**********************************************************
%**********************************************************
% Generating, dependent uniform-distributed random numbers:

do_gen=1;
while(do_gen==1)
    ZZStrom = input('ZZStrom   ');
    rand('state',ZZStrom);   % 1. ZZStrom

    n_gen = input('generated sample size ');
    u_gen = zeros(D,n_gen);         %copula
    u_gen(1,:) = rand(1,n_gen);

    power_K = realpow(K,D-1);

    for i=1:n_gen
        u_gen_ = rand(1,1);
        u_power_K = u_gen_ * power_K;

        summ = 0;
        j2 = 0;
        u1_gen_up = ceil(u_gen(1,i)*K);   % delta = 1/K
        while (summ < u_power_K) 
            j2 = j2+1;
            summ_old = summ;
            summ = summ + f2(u1_gen_up, j2);
        end;
        u2_down = j2-1;
        u2_gen_up = j2;
        u_gen_rest = u_gen_ - summ_old/power_K;
        u2_gen_rest = u_gen_rest*realpow(K,D-2) / f2(u1_gen_up, u2_gen_up);
        u_gen(2,i) = u2_down*delta + u2_gen_rest;

        if D>2
        u_gen_ = rand(1,1);
        summ = 0;
        j3 = 0;
        u_gen_f2u1u2 = u_gen_*f2(u1_gen_up,u2_gen_up);

        while (summ < u_gen_f2u1u2) 
            j3 = j3+1;
            summ_old = summ;
            summ = summ + f3(u1_gen_up*K+u2_gen_up, j3);
        end;
        u3_down = j3-1;
        u3_gen_up = j3;
        u_gen_rest = delta*(u_gen_f2u1u2 - summ_old);
        u3_gen_rest = u_gen_rest / f3(u1_gen_up*K+u2_gen_up, u3_gen_up);
        u_gen(3,i) = u3_down*delta + u3_gen_rest;
        end; % D>2

        if D>3
        u_gen_ = rand(1,1);
        summ = 0;
        j4 = 0;
        u_gen_f3u1u2u3 = u_gen_*f3(u1_gen_up*K+u2_gen_up,u3_gen_up);

        while (summ < u_gen_f3u1u2u3) 
            j4 = j4+1;
            summ_old = summ;
            summ = summ + f4(u1_gen_up*K+u2_gen_up,u3_gen_up*K+ j4);
        end;
        u4_down = j4-1;
        u_gen_rest = delta*(u_gen_f3u1u2u3 - summ_old); 
        u4_gen_rest = u_gen_rest / f4(u1_gen_up*K+u2_gen_up,u3_gen_up*K+u4_down+1);
        u_gen(4,i) = u4_down*delta + u4_gen_rest;
        end; % D>3

    end;

    %*********************************************************************%
    %*********************************************************************%
    % Plot the simulated dependent uniform-distributed random numbers 
        nbr=n_gen;


    for d=1:bilder
        for ds=d+1:D
            %plot u_gen(ds,:) )  against   u_gen(d,:) 
             figure ('name','Gen. dependent uniform rv.',...
            'Position',PosFigure_3(:,d,ds))
            scatter(u_gen(d,1:nbr),u_gen(ds,1:nbr),2.5,'*')
            axis([-.1 1.1 -.1 1.1])
            xlabel(d)
            ylabel(ds)
        end
    end


    %***********************************************************
    %**********************************************************
    %**********************************************************
    % Generating, get dependent random numbers using inverse transfom:

    z_gen = zeros(D,n_gen);

    for d=1:D
        Dimension=d
        disp('Inverse transform u(d,:)');
        disp('   with linearly interpolated empirical CDF (1)');
        disp('   for discrete values with empirical CDF (4)');

        menu = input('   to Weibull (2), to log-logistic (3) ');

    switch menu
        case 1
            for j=1:n_gen
                un = u_gen(d,j)*n;
                i = ceil(un)-1;
                if i==0
                    z_gen(d,j) = (un-i)*Z_sorted(d,i+1);
                else 
                    z_gen(d,j) = Z_sorted(d,i) + (un-i)*(Z_sorted(d,i+1)-Z_sorted(d,i));
                end;
            end;
        case 2
            alpha  = input('shape parameter alpha ');
            beta  = input('scale parameter beta ');
            for i=1:n_gen
                z_gen(d,i) = beta * (-log(1-u_gen(d,i) ))^(1/alpha);
            end;
        case 3
            alpha  = input('shape parameter alpha ');
            beta  = input('scale parameter beta ');
            z_gen(d,:) = beta / (1./ u_gen(d,:) - 1).^(1/alpha);
        case 4
            for j=1:n_gen
                un = u_gen(d,j)*n;
                i = ceil(un);
                z_gen(d,j) = Z_sorted(d,i);
            end;

    end;
    end; %d

    size_z_gen = size(z_gen)

    %***********************************************************
    %***********************************************************

    r=1;
    for d=1:D
        %calculate means and the correlations.
       Mean_Z(d)=mean(Z(d,:));
       v = var(Z(d,:));
       VarKoeffZ(d) = sqrt(v)/Mean_Z(d);

       Mean_Z_generated(d)=mean(z_gen(d,:));
       Abweichung(r) = (Mean_Z(d) - Mean_Z_generated(d))/Mean_Z_generated(d);
       r=r+1;
       v = var(z_gen(d,:));
       VarKoeffZ_gen(d) = sqrt(v)/Mean_Z_generated(d);
       Abweichung(r) = (VarKoeffZ(d) - VarKoeffZ_gen(d));
                     %  Better: Relative deviation as well
       r=r+1;
    end

    for d=1:D   
        %estimate the correlations of the original time series:
       for ds=d+1:D
          Corr1=corrcoef(Z(d,:), Z(ds,:));
          Corr_Coeffizient(d,ds)=Corr1(1,2);
       end
       for ds=d+1:D
          Corr2=corrcoef(z_gen(d,:), z_gen(ds,:));
          Corr_Coeffizient_gen(d,ds)=Corr2(1,2);
          Abweichung(r) = (Corr_Coeffizient(d,ds)...
              - Corr_Coeffizient_gen(d,ds));
          r=r+1;
       end

    end

    Mean_Z
    Mean_Z_generated
    VarKoeffZ
    VarKoeffZ_gen
    Corr_Coeffizient
    Corr_Coeffizient_gen

    Abweichung = abs(Abweichung)
    display('(Relative deviations for the means, absolute otherwise');
    max_Abweichung  = max (Abweichung)

    %E_Z1_Z2_square = (Z(1,:).*(Z(2,:).*Z(2,:)))*ones(n,1)/n
    %E_Z1_Z2_square_generated = ( z_gen(1,:).*(z_gen(2,:).*z_gen(2,:)) )*ones(n_gen,1)/n_gen
    %display(' Abweichung (E_Z1_Z2_square - E_Z1_Z2_square_generated), relativ');
    %(E_Z1_Z2_square - E_Z1_Z2_square_generated)/E_Z1_Z2_square_generated

    %*********************************************************************%
    %*********************************************************************%

    %Plot the simulated dependent random variables


    for d=1:bilder
        for ds=d+1:D
            figure ('name', 'Generated Time Series', 'Position', PosFigure_4(:,d,ds))
            scatter(z_gen(d,1:nbr),z_gen(ds,1:nbr),2.5,'*')
            axis([0 Achse(d) 0 Achse(ds)])
            xlabel(d); 
            ylabel(ds)
        end
    end


    do_gen = input('another generating? (0 or 1) ');
end; %Loop for generating

end

%********************************************************************

function X = Stichprobe_Weibull_dim_2(alpha, beta,n)

display('Sample:');
display( 'X(1) ~ Weibull(alpha, beta)');
display('X(2,i) = X(1,i) + rand(1,1)*X(1,i)' );


alpha
beta
n   % Stichprobengroesse

X = zeros(2,n);

for i=1:n
    X(1,i) =  beta * (-log(rand(1,1)))^(1/alpha) ;
    X(2,i) = X(1,i) + rand(1,1)*X(1,i);
end

sample_mean_1 = mean(X(1,:));
sample_mean_2 = mean(X(2,:));

end

%**********************************************************************
function [Z_sorted, nU, F_emp] = empVF_v3(n,D,Z)

% D Dimensionen, in any empirical VF

%Empirical distribution:

step=[1/n:1/n:1]';  % [1/n, 2/n, ... ,1]

for d=1:D
   [Z_sorted(d,:), index_alt(d,:)] = sort(Z(d,:));
end 


U=zeros(D,n);     % to speed up processing

%Find  empirical CDF  F_emp
F_emp = zeros(4,n);
for d=1:D
    F_emp(d,:) = step';
end

  % Sampling values repeatedly have different values F i/n, (i+1)/n, ...
for i=1:n
    for d=1:D        
        nU(d,index_alt(d,i)) =  i; 
        % values n*U(d,index_alt) in order to omit rounding errors
    end                
end

end
