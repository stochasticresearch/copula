clear;
clc;

clear;
clc;

if(ispc)
    datafileDir = 'C:\Users\Kiran\Desktop\datasets\cu-antenna-data-200905';
elseif(ismac)
    datafileDir = '/Users/kiran/Desktop/datasets/cu-antenna-data-200905';
    saveDir = 'Users/kiran/stochasticresearch/copula/data';
elseif(isunix)
    datafileDir = '/home/kiran/Desktop/datasets/cu-antenna-data-200905';
    saveDir = '/home/kiran/stochasticresearch/copula/data';
else
end

%% Pre-Process CRAWDAD Antenna Dataset

% read in the crawdad cu antenna observations dataset into Matlab
filename = fullfile(datafileDir, 'packets.txt');
fid = fopen(filename,'r');
format = '%s %d %s %d %s %f %f';
pktsData = textscan(fid, format);

rss = pktsData{2};
antenna = pktsData{3};
azimuthAngle = pktsData{4};
tag = pktsData{5};
rssDiff = pktsData{7};

% convert antenna to ordinal
antennaTypes = unique(antenna);
valueSet = 1:length(antennaTypes);
mapObjAntenna = containers.Map(antennaTypes,valueSet);

antennaOrdinal = zeros(size(antenna));
for ii=1:length(antenna)
    antennaOrdinal(ii) = mapObjAntenna(antenna{ii});
end

% convert the tag to ordinal data
tagTypes = unique(tag);
valueSet = 1:length(tagTypes);
mapObjTag = containers.Map(tagTypes,valueSet);

tagOrdinal = zeros(size(tag));
for ii=1:length(tag)
    tagOrdinal(ii) = mapObjTag(tag{ii});
end

% save the data so we don't have to process the raw data everytime
save(fullfile(saveDir, 'cu_antenna_200905'));

%% Do fun stuff w/ the pre-processed data :D