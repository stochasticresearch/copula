function [] = testDetermineOrientation()

xMin = 0;
xMax = 1;

num_noise = 30;                    % The number of different noise levels used
noise = 3;                         % A constant to determine the amount of noise
M = 500;

numTestsForConfidence = 300;

%%%%%%%%%%%%%%%%%

noise_level = 1;

orientation0 = 0; orientation1 = 0;
for ii=1:numTestsForConfidence
    x = rand(M,1)*(xMax-xMin)+xMin;
    y = 4*(x-.5).^2 + noise*(noise_level/num_noise)*randn(M,1);
    x = pobs(x); y = pobs(y);
    orientation0 = orientation0 + determineOrientation(x,y);
    orientation1 = orientation1 + determineOrientation(y,x);
end
fprintf('quadratic(1) -- x/y correct = %0.02f %%, y/x correct = %0.02f %% \n', ...
   (numTestsForConfidence-orientation0)/numTestsForConfidence*100, ...
   orientation1/numTestsForConfidence*100);

orientation0 = 0; orientation1 = 0;
for ii=1:numTestsForConfidence
    x = rand(M,1)*(xMax-xMin)+xMin;
    y = sin(4*pi*x) + 2*noise*(noise_level/num_noise)*randn(M,1);
    x = pobs(x); y = pobs(y);
    orientation0 = orientation0 + determineOrientation(x,y);
    orientation1 = orientation1 + determineOrientation(y,x);
end
fprintf('sin(1) -- x/y correct = %0.02f %%, y/x correct = %0.02f %% \n', ...
   (numTestsForConfidence-orientation0)/numTestsForConfidence*100, ...
   orientation1/numTestsForConfidence*100);

%%%%%%%%%%%%%%%%%

noise_level = 10;

orientation0 = 0; orientation1 = 0;
for ii=1:numTestsForConfidence
    x = rand(M,1)*(xMax-xMin)+xMin;
    y = 4*(x-.5).^2 + noise*(noise_level/num_noise)*randn(M,1);
    x = pobs(x); y = pobs(y);
    orientation0 = orientation0 + determineOrientation(x,y);
    orientation1 = orientation1 + determineOrientation(y,x);
end
fprintf('quadratic(10) -- x/y correct = %0.02f %% , y/x correct = %0.02f %% \n', ...
   (numTestsForConfidence-orientation0)/numTestsForConfidence*100, ...
   orientation1/numTestsForConfidence*100);

orientation0 = 0; orientation1 = 0;
for ii=1:numTestsForConfidence
    x = rand(M,1)*(xMax-xMin)+xMin;
    y = sin(4*pi*x) + 2*noise*(noise_level/num_noise)*randn(M,1);
    x = pobs(x); y = pobs(y);
    orientation0 = orientation0 + determineOrientation(x,y);
    orientation1 = orientation1 + determineOrientation(y,x);
end
fprintf('sin(10) -- x/y correct = %0.02f %% , y/x correct = %0.02f %% \n', ...
   (numTestsForConfidence-orientation0)/numTestsForConfidence*100, ...
   orientation1/numTestsForConfidence*100);

%%%%%%%%%%%%%%%%%

noise_level = 15;

orientation0 = 0; orientation1 = 0;
for ii=1:numTestsForConfidence
    x = rand(M,1)*(xMax-xMin)+xMin;
    y = 4*(x-.5).^2 + noise*(noise_level/num_noise)*randn(M,1);
    x = pobs(x); y = pobs(y);
    orientation0 = orientation0 + determineOrientation(x,y);
    orientation1 = orientation1 + determineOrientation(y,x);
end
fprintf('quadratic(15) -- x/y correct = %0.02f %% , y/x correct = %0.02f %% \n', ...
   (numTestsForConfidence-orientation0)/numTestsForConfidence*100, ...
   orientation1/numTestsForConfidence*100);

orientation0 = 0; orientation1 = 0;
for ii=1:numTestsForConfidence
    x = rand(M,1)*(xMax-xMin)+xMin;
    y = sin(4*pi*x) + 2*noise*(noise_level/num_noise)*randn(M,1);
    x = pobs(x); y = pobs(y);
    orientation0 = orientation0 + determineOrientation(x,y);
    orientation1 = orientation1 + determineOrientation(y,x);
end
fprintf('sin(15) -- x/y correct = %0.02f %% , y/x correct = %0.02f %% \n', ...
   (numTestsForConfidence-orientation0)/numTestsForConfidence*100, ...
   orientation1/numTestsForConfidence*100);

%%%%%%%%%%%%%%%%%

noise_level = 20;

orientation0 = 0; orientation1 = 0;
for ii=1:numTestsForConfidence
    x = rand(M,1)*(xMax-xMin)+xMin;
    y = 4*(x-.5).^2 + noise*(noise_level/num_noise)*randn(M,1);
    x = pobs(x); y = pobs(y);
    orientation0 = orientation0 + determineOrientation(x,y);
    orientation1 = orientation1 + determineOrientation(y,x);
end
fprintf('quadratic(20) -- x/y correct = %0.02f %% , y/x correct = %0.02f %% \n', ...
   (numTestsForConfidence-orientation0)/numTestsForConfidence*100, ...
   orientation1/numTestsForConfidence*100);

orientation0 = 0; orientation1 = 0;
for ii=1:numTestsForConfidence
    x = rand(M,1)*(xMax-xMin)+xMin;
    y = sin(4*pi*x) + 2*noise*(noise_level/num_noise)*randn(M,1);
    x = pobs(x); y = pobs(y);
    orientation0 = orientation0 + determineOrientation(x,y);
    orientation1 = orientation1 + determineOrientation(y,x);
end
fprintf('sin(20) -- x/y correct = %0.02f %% , y/x correct = %0.02f %% \n', ...
   (numTestsForConfidence-orientation0)/numTestsForConfidence*100, ...
   orientation1/numTestsForConfidence*100);


end

function [orientation] = determineOrientation(ax1pts, ax2pts, checkBoxSize)
% if this function returns 0, then ax1pts is the independent variable. 
% else, ax2pts is the independent variable.

if(nargin<3)
    checkBoxSize = 0.1;
end

seePlots = 0;


range_ax1_ax2 = zeros(1,round(1/checkBoxSize));
range_ax2_ax1 = zeros(1,round(1/checkBoxSize));

% first check u/v
ax1min = 0; ax1max = checkBoxSize; ax2min = 0; ax2max = 1;
ii = 1;
if(seePlots)
    f = figure;
end
while(ax1max<1)
    ax1_match = find(ax1pts>ax1min & ax1pts<ax1max);
    ax2_match = find(ax2pts>ax2min & ax2pts<ax2max);
    matchIdxs = intersect(ax1_match,ax2_match);
    ax2matchPts = ax2pts(matchIdxs);
    
    % compute the number of unique points in ax1
    range_ax1_ax2(ii) = quantile(ax2matchPts, 0.75)-quantile(ax2matchPts,0.25);
    
    if(seePlots)
        clf(f);
        figure(f);
        scatter(ax1pts, ax2pts, 'b');
        hold on;
        scatter(ax1pts(matchIdxs), ax2matchPts, 'r'); grid on;
        title(sprintf('U/V -- Range = %0.02f\n', range_ax1_ax2(ii)));
        pause;
    end
    
    ii = ii + 1;
    
    ax1min = ax1min + checkBoxSize;
    ax1max = ax1max + checkBoxSize;
end

% now check v/u
tmp = ax1pts; ax1pts = ax2pts; ax2pts = tmp;
ax1min = 0; ax1max = checkBoxSize; ax2min = 0; ax2max = 1;
ii = 1;
while(ax1max<1)
    ax1_match = find(ax1pts>ax1min & ax1pts<ax1max);
    ax2_match = find(ax2pts>ax2min & ax2pts<ax2max);
    matchIdxs = intersect(ax1_match,ax2_match);
    ax2matchPts = ax2pts(matchIdxs);
    
    % compute the number of unique points in ax1
    range_ax2_ax1(ii) = quantile(ax2matchPts, 0.75)-quantile(ax2matchPts,0.25);
    
    if(seePlots)
        clf(f);
        figure(f);
        scatter(ax1pts, ax2pts, 'b');
        hold on;
        scatter(ax1pts(matchIdxs), ax2matchPts, 'r'); grid on;
        title(sprintf('V/U -- Range = %0.02f\n', range_ax2_ax1(ii)));
        pause;
    end
    
    ii = ii + 1;
    
    ax1min = ax1min + checkBoxSize;
    ax1max = ax1max + checkBoxSize;
end


if(sum(range_ax1_ax2) < sum(range_ax2_ax1))
    orientation = 0;
else
    orientation = 1;
end

end