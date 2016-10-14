function r = rdc(x,y,k,s)
%RDC - Randomized Dependence Coefficient, implementation retrieved from
% http://is.tuebingen.mpg.de/fileadmin/user_upload/files/publications/2013/Lopez-Paz_talk_rdc.pdf
% See also: https://papers.nips.cc/paper/5138-the-randomized-dependence-coefficient.pdf
% Default parameters: k=20,s=1/6 (from the paper reference, where they used
%                                 these parameters to conduct the tests)

% The original code for RDC was obtained from David Lopez-Paz's GITHUB for RDC.  The URL is:
% git@github.com:lopezpaz/randomized_dependence_coefficient.git

% We put some stuff in here to supress warnings.  The warning ID was
% collected as follows:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> rdc(x,y+randn(500,1))
%    Warning: X is not full rank. 
%      In canoncorr (line 88)
%      In rdc (line 23) 
%    Warning: Y is not full rank. 
%      In canoncorr (line 96)
%      In rdc (line 23) 
%
%w = warning('query','last')
%
%w = 
%
%    identifier: 'stats:canoncorr:NotFullRank'
%         state: 'on'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch(nargin)
    case 4
    case 3
        s = 1/6;
    case 2
        k = 20; s = 1/6;
    otherwise
        error('Must specify atleast x & y!');
end

n = size(x,1);
x = [tiedrank(x)/n ones(n,1)];
y = [tiedrank(y)/n ones(n,1)];
x = sin(s/size(x,2)*x*randn(size(x,2),k));
y = sin(s/size(y,2)*y*randn(size(y,2),k));

warningID = 'stats:canoncorr:NotFullRank';
warning('off',warningID);
[~,~,r] = canoncorr([x ones(n,1)],[y ones(n,1)]);
warning('on',warningID);

r = max(r);     % should we return the maximum cannonical correlation? this
                % seems correct to me based on the paper ...

end
