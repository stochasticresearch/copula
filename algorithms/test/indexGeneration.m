clear;
clc;

x = [10 10];
y = zeros(prod(x), length(x));
for ii=1:prod(x)
    remVal = ii-1;
    divVal = prod(x(1:end-1));
    for jj=1:length(x)
        y(ii,jj) = floor( remVal/divVal );
%         fprintf('ii=%d remVal=%d divVal=%d y(%d,%d)=%d \n', ...
%             ii, remVal, divVal, ii, jj, y(ii,jj));
        
        remVal = remVal-y(ii,jj)*divVal;
        divVal = divVal / x(1);
    end
%     fprintf('*****************************\n');
end
y = y+1;
y = y / x(end);
% fliplr(y)