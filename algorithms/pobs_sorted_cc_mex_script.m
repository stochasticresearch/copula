% POBS_SORTED_CC_MEX_SCRIPT   Generate MEX-function pobs_sorted_cc_mex from
%  pobs_sorted_cc.
% 
% Script generated from project 'pobs_sorted_cc.prj' on 10-May-2019.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.EnableJIT = true;

%% Define argument types for entry-point 'pobs_sorted_cc'.
ARGS = cell(1,1);
ARGS{1} = cell(2,1);
ARGS{1}{1} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(0,[Inf  1],[1 0]);

%% Invoke MATLAB Coder.
codegen -config cfg pobs_sorted_cc -args ARGS{1} -nargout 2

