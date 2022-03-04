%   Load MRI data and predefined options structure
load('InVivo_4p7T.mat');


%%% FIT SINGLE VOXEL %%%
S = squeeze(img(146,63,1,:));
opt.debug = 1;
opt.FitType = 'nnls';
[T2,B1,amp] = StimFitNNLS(S,opt);


%%% FIT ENTIRE IMAGE %%%
% opt.debug = 0;
% opt.th = 0;
% opt.FitType = 'nnls';
% parpool;
% [T2,B1,amp] = StimFitImgNNLS(img,opt);
