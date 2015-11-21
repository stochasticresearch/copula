% The main script, which performs simulations for both the BN and HCBN
% variants with the specified data

% TODO: [ ] - Figure out why copula generation balks when training size >
%             generation size?  Must be a bug in the regeneration part

clear;
clc;

dataset = 'adult';      % must be one of the following choices:
                        % {adult, bands, crx, imports}

adultDataFile = 'adult.data';
adultDataFmt = '%f%s%f%s%f%s%s%s%s%s%f%f%f%s%s';
adultContinuousIdx = [1,3,5,11,12,13];
adultDiscreteIdx = [2,4,6,7,8,9,10,14,15];

bandsDataFile = 'bands.data';
bandsDataFmt = '%f%s%s%f%s%s%s%s%s%s%s%s%s%s%s%f%f%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s';
bandsContinuousIdx = [1,21:39];
bandsDiscreteIdx = [2:20,40];

crxDataFile = 'crx.data';
crxDataFmt = '%s%f%f%s%s%s%s%f%s%s%f%s%s%f%f%s';
crxContinuousIdx = [2,3,8,11,14,15];
crxDiscreteIdx = [1,4,5,6,7,9,10,12,13,16];

importsDataFile = 'imports-85.data';
importsDataFmt = '%f%f%s%s%s%s%s%s%s%f%f%f%f%f%s%s%f%s%f%f%f%f%f%f%f%f';
importsContinuousIdx = [2,10,11,12,13,14,17,19:26];
importsDiscreteIdx = [1,3,4,5,6,7,8,9,15,16,18];

if(strcmpi(dataset,'adult'))
    datafile = adultDataFile;
    datafmt = adultDataFmt;
    continuousIdx = adultContinuousIdx;
    discreteIdx = adultDiscreteIdx;
elseif(strcmpi(dataset,'bands'))
    datafile = bandsDataFile;
    datafmt = bandsDataFmt;
    continuousIdx = bandsContinuousIdx;
    discreteIdx = bandsDiscreteIdx;
elseif(strcmpi(dataset,'crx'))
    datafile = crxDataFile;
    datafmt = crxDataFmt;
    continuousIdx = crxContinuousIdx;
    discreteIdx = crxDiscreteIdx;
elseif(strcmpi(dataset,'imports'))
    datafile = importsDataFile;
    datafmt = importsDataFmt;
    continuousIdx = importsContinuousIdx;
    discreteIdx = importsDiscreteIdx;
end

fprintf('Reading Dataset\n');
fid = fopen(datafile);
rawdata = textscan(fid, datafmt, 'delimiter', ',');
fclose(fid);

sz = size(rawdata);
numDataPts = length(rawdata{1});
numFeatures = sz(2);

% convert the discrete random variables, which are represented as strings,
% into ordinal numbers
fprintf('Processing Dataset for Model Generation\n');
dataMat = zeros(sz);
dataMatCtr = 1;
for ii=1:numDataPts
    datarow = zeros(1,numFeatures);
    for jj=1:numFeatures
        pt = rawdata{jj};
        if(isreal(pt))
            datarow(jj) = pt(ii);
        else
            x = mapStrToOrd(jj,pt{ii}, dataset);
            datarow(jj) = x;
        end
    end
    if(any(datarow==-999))
        % don't store this
    else
        dataMat(dataMatCtr,:) = datarow;
        dataMatCtr = dataMatCtr + 1;
    end
end

%% 
% Learn the Structure
% discretize the continuous nodes
dataTmp = dataMat;
for ii=continuousIdx
    x = dataMat(:,ii);
    % discretize
    edges = linspace(min(x),max(x),10);
    nbins = length(edges)-1;
    [~,bins] = histc(x,edges);
    bins(bins==nbins+1) = nbins;
    y = reshape(bins,size(x));
    dataTmp(:,ii) = y;
end
dataTmp = dataTmp';     % the learning algo likes the data in this format

% learn the structure using k2 algorithm
ns = zeros(1,size(dataTmp,1));
for ii=1:length(ns)
    ns(ii) = length(unique(dataTmp(ii,:)));
end

order = 1:size(dataTmp,1);
order_str = {'age','workclass','fnlwgt','education','education-num','martial-status','occupation',...
             'relationship', 'race', 'sex', 'capital-gain', 'capital-loss', 'hours-per-week', 'native-country', 'salary'};
max_fan_in = 2;

% print out structure
dag = learn_struct_K2(dataTmp, ns, order, 'max_fan_in', max_fan_in, 'verbose', 'no');
% print out the connections
for ii=1:length(order_str)
    fprintf('Parent = %s\n', order_str{ii});
    idxs = find(dag(ii,:));
    for zz=idxs
        fprintf(' ---> %s\n', order_str{zz});
    end
    fprintf('--------------------------------\n');
end

%% Extract a convenient child & parents for exposition of our technique
X_empirical = [dataMat(:,7) dataMat(:,11)];
X_empirical = X_empirical(1:10000,:);        % training data size

n = size(X_empirical,1);

n_gen = 1000;

fprintf('Applying CLG Model\n')
[x1_domain, x1_multinomial_est, x2_mle_params] = fit_clg_2D( X_empirical );
X_clg_gen = gen_samples_clg_2D(x1_domain, x1_multinomial_est, x2_mle_params, n_gen);

fprintf('Applying Copula Model\n')
X1_continued = X_empirical(:,1) + (rand(n,1)-1);
X_transformed = [X1_continued X_empirical(:,2)];
% K = n_gen;
K = 200;
D = size(X_transformed,2);

[ U_gen, ~, U_emp ] = emp_copularnd_old2( X_transformed, n_gen, K );
% [C,U,c] = empcopula(X_transformed,K);    
% U_gen = empcopularnd(c, n);

X_gen = empdistrnd(U_gen, X_empirical);

% Visualize Results
subplot(1,3,1);
scatter(X_empirical(1:1000,1),X_empirical(1:1000,2));
grid on
xlabel('Occupation')
ylabel('Capital Gains')
title('Empirical Data')

subplot(1,3,2);
scatter(X_gen(:,1),X_gen(:,2));
grid on
xlabel('Occupation')
ylabel('Capital Gains')
title('Copula Generative Model')

subplot(1,3,3);
scatter(X_clg_gen(:,1),X_clg_gen(:,2));
grid on
xlabel('Occupation')
ylabel('Capital Gains')
title('CLG Generative Model')

%% Some additional extractions for fun
% [race age fnlwgt]
X_empirical = [dataMat(:,9) dataMat(:,1) dataMat(:,3)];
X_empirical = X_empirical(1:10000,:);        % training data size

n = size(X_empirical,1);

n_gen = 1000;

fprintf('Applying CLG Model\n')
[ x1_domain, x1_multinomial_est, mle_params ] = fit_clg_3D( X_empirical );
[X1X2_clg_gen, X1X3_clg_gen, X2X3_clg_gen] = ...
    gen_samples_clg_3D(x1_domain, x1_multinomial_est, mle_params, n_gen);

fprintf('Applying Copula Model\n')
X1_continued = X_empirical(:,1) + (rand(n,1)-1);
X_transformed = [X1_continued X_empirical(:,2) X_empirical(:,3)];
K = n_gen;
D = size(X_transformed,2);
[ U_gen, ~, U_emp ] = emp_copularnd( X_transformed, n_gen, K );
X_gen = empdistrnd(U_gen, X_empirical);

% Visualize Results

figure;
subplot(2,1,1)
scatter(X_empirical(1:n_gen,1),X_empirical(1:n_gen,2))
xlabel('Race')
ylabel('Age')
grid on
title('Empirical Data')

subplot(2,1,2)
scatter(X_empirical(1:n_gen,2),X_empirical(1:n_gen,3))
xlabel('Age')
ylabel('fnlwgt')
grid on

figure;
subplot(2,1,1)
scatter(X1X2_clg_gen(:,1),X1X2_clg_gen(:,2))
xlabel('Race')
ylabel('Age')
grid on
title('CLG Generative Model')

subplot(2,1,2)
scatter(X2X3_clg_gen(:,1),X2X3_clg_gen(:,2))
xlabel('Age')
ylabel('fnlwgt')
grid on

figure;
subplot(2,1,1)
scatter(X_gen(:,1),X_gen(:,2))
xlabel('Race')
ylabel('Age')
grid on
title('HCBN Model')

subplot(2,1,2)
scatter(X_gen(:,2),X_gen(:,3))
xlabel('Age')
ylabel('fnlwgt')
grid on


%% R - Style Plots
figure;
% do an R-style "pairs" plot of the empirical data
subplot(3,3,1)
axis([0 1 0 1])
text(0.5,0.5,'X_1 (Race)')

subplot(3,3,2)
scatter(X_empirical(1:n_gen,1),X_empirical(1:n_gen,2))
grid on
title('Empirical Data (Census Data)')

subplot(3,3,3)
scatter(X_empirical(1:n_gen,1),X_empirical(1:n_gen,3))
grid on

subplot(3,3,4)
scatter(X_empirical(1:n_gen,2),X_empirical(1:n_gen,1))
grid on

subplot(3,3,5)
axis([0 1 0 1])
text(0.5,0.5,'X_2 (Age)')

subplot(3,3,6)
scatter(X_empirical(1:n_gen,2),X_empirical(1:n_gen,3))
grid on

subplot(3,3,7)
scatter(X_empirical(1:n_gen,3),X_empirical(1:n_gen,1))
grid on

subplot(3,3,8)
scatter(X_empirical(1:n_gen,3),X_empirical(1:n_gen,2))
grid on

subplot(3,3,9)
axis([0 1 0 1])
text(0.5,0.5,'X_3 (fnlwgt)')

% do an R-style "pairs" plot of the CLG Generated data
figure;
subplot(3,3,1)
axis([0 1 0 1])
text(0.5,0.5,'X_1 (Race)')

subplot(3,3,2)
scatter(X1X2_clg_gen(:,1),X1X2_clg_gen(:,2))
grid on
title('CLG Generative Model (Census Data)')

subplot(3,3,3)
scatter(X1X3_clg_gen(:,1),X1X3_clg_gen(:,2))
grid on

subplot(3,3,4)
scatter(X1X2_clg_gen(:,2),X1X2_clg_gen(:,1))
grid on

subplot(3,3,5)
axis([0 1 0 1])
text(0.5,0.5,'X_2 (Age)')

subplot(3,3,6)
scatter(X2X3_clg_gen(:,1),X2X3_clg_gen(:,2))
grid on

subplot(3,3,7)
scatter(X1X3_clg_gen(:,2),X1X3_clg_gen(:,1))
grid on

subplot(3,3,8)
scatter(X2X3_clg_gen(:,2),X2X3_clg_gen(:,1))
grid on

subplot(3,3,9)
axis([0 1 0 1])
text(0.5,0.5,'X_3 (fnlwgt)')

% do an R-style "pairs" plot of the Copula Generated data
figure;
subplot(3,3,1)
axis([0 1 0 1])
text(0.5,0.5,'X_1 (Race)')

subplot(3,3,2)
scatter(X_gen(:,1),X_gen(:,2))
grid on
title('Copula Generative Data (Census Data)')

subplot(3,3,3)
scatter(X_gen(:,1),X_gen(:,3))
grid on

subplot(3,3,4)
scatter(X_gen(:,2),X_gen(:,1))
grid on

subplot(3,3,5)
axis([0 1 0 1])
text(0.5,0.5,'X_2 (Age)')

subplot(3,3,6)
scatter(X_gen(:,2),X_gen(:,3))
grid on

subplot(3,3,7)
scatter(X_gen(:,3),X_gen(:,1))
grid on

subplot(3,3,8)
scatter(X_gen(:,3),X_gen(:,2))
grid on

subplot(3,3,9)
axis([0 1 0 1])
text(0.5,0.5,'X_3 (fnlwgt)')
