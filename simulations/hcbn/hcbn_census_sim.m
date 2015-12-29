% simulation w/ HCBN and census data
clear;
clc;
close all;

dataset = 'adult';

adultDataFile = 'adult.data';
adultDataFmt = '%f%s%f%s%f%s%s%s%s%s%f%f%f%s%s';
adultContinuousIdx = [1,3,5,11,12,13];
adultDiscreteIdx = [2,4,6,7,8,9,10,14,15];

datafile = adultDataFile;
datafmt = adultDataFmt;
continuousIdx = adultContinuousIdx;
discreteIdx = adultDiscreteIdx;

fprintf('Generating Dataset\n');
fid = fopen(datafile);
rawdata = textscan(fid, datafmt, 'delimiter', ',');
fclose(fid);

sz = size(rawdata);
numDataPts = length(rawdata{1});
numFeatures = sz(2);

% convert the discrete random variables, which are represented as strings,
% into ordinal numbers
fprintf('Processing Dataset for Model Generation\n');
X = zeros(sz);
idx = 1;
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
        X(idx,:) = datarow;
        idx = idx + 1;
    end
end

X_train = X(1:2000,:);
X_test = X(1:1000,:);
fprintf('Generating HCBN model\n');
% setup HCBN object
nodeNames = {'age', 'workclass', 'fnlwgt', 'education', 'education-num', ...
    'martial-status', 'occupation', 'relationship', 'race', 'sex', ...
    'capital-gain', 'capital-loss', 'hours-per-week', 'native-country', ...
    'salary-range'};
bntPath = '../bnt';
discreteNodeNames = {nodeNames{2}, nodeNames{4}, nodeNames{6}, nodeNames{7}, ...
    nodeNames{8}, nodeNames{9}, nodeNames{10}, nodeNames{14}, nodeNames{15}};
D = length(nodeNames);

dag = zeros(D,D);
% construct the DAG so we can show interesting graphs, we'd like to show
% the probabilistic modeling of age & occupation, and age & marital-status
% versus the collected data (this is a continuous parent w/ discrete
% children scenario)
dag(3,6) = 1; dag(3,8) = 1;

% instantiate hcbn object
hcbnObj = hcbn(bntPath, X_train, nodeNames, discreteNodeNames, dag);

fprintf('Generating samples from model\n');
tvec1 = hcbnObj.genFamilySamples('martial-status', 1000);
tvec2 = hcbnObj.genFamilySamples('relationship', 1000);

h1 = subplot(2,2,1); 
scatter(X_test(:,3),X_test(:,6)); grid on; xlabel('FNLWGT'); ylabel('Martial Status'); title('True Data')
h2 = subplot(2,2,2); 
scatter(tvec1(:,2),tvec1(:,1)); grid on; xlabel('FNLWGT'); ylabel('Martial Status'); title('HCBN Generative Model')
h3 = subplot(2,2,3); 
scatter(X_test(:,3),X_test(:,8)); grid on; xlabel('FNLWGT'); ylabel('Relationship'); title('True Data')
h4 = subplot(2,2,4); 
scatter(tvec2(:,2),tvec2(:,1)); grid on; xlabel('FNLWGT'); ylabel('Relationship'); title('HCBN Generative Model')

linkaxes([h1,h2],'xy')
linkaxes([h3,h4],'xy')