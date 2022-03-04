%   Load MRI data and predefined options structure
load('InVivo_4p7T.mat');



%%% FIT SINGLE VOXEL %%%
S = squeeze(img(91,69,1,:));
opt.debug = 1;
opt.FitType = 'lsq';
[T2,B1,amp] = StimFit(S,opt);


%%% FIT ENTIRE IMAGE %%%
%opt.debug = 0;
%opt.th = 0;
%opt.FitType = 'lsq';
%parpool;
%[T2,R2,B1,amp] = StimFitImg(img,opt);
