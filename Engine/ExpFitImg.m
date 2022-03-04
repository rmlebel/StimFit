function [T2, R2, AMP, opt] = ExpFitImg(img,opt)
%   T2 fitting routine with stimulated echo compensation
%   
%   Author: R. Marc Lebel
%   Date:   07/2011
%   
%   Usage:
%   [T2, R2, AMP, opt] = StimFitImg(img,opt)
%   
%   Input:
%   img: Img to be fit (size: NP x NV x NS x NE)
%   opt: Options structure defined by StimFit_optset (optional)
%   
%   Output:
%   T2:  Decay time (s) (size: NP x NV x NS)
%   R2:  Inverse of T2
%   AMP: Relative amplitude (au) (size: NP x NV x NS)
%   opt: Modified options structure


%   Check input arguments
if nargin < 1
    error('Function requires at least 1 input');
end
if nargin < 2 || isempty(opt) || ~isstruct(opt)
    opt = StimFit_optset;
end


%   Get and store image size, check for consistency
[opt.size(1),opt.size(2),opt.size(3),opt.size(4)] = size(img);
if opt.etl ~= opt.size(4)
    warning('StimFitImg:InconsistentParameters',...
        'opt.etl (=%g) has been modified to match image size (=%g)',opt.etl,opt.size(4));
    opt.etl = opt.size(4);
end


%   Check that image is real and double
if ~isreal(img)
    img = abs(img);
    warning('StimFitImg:ComplexImage','Image is complex. Taking absolute value');
end
if isa(img,'single')
    img = double(img);
end

%   Scale image. The fit routine is bounded, so arbitrarily scaled images
%   are not guaranteed to work - we want the signal to be within an order
%   of magnitude (or so) of unity
disp('Scaling Image');
i1 = img(:,:,:,1);
img = img ./ mean(i1(:));
clear i1

%   Remove noise baseline and define a fitting threshold
%   Estimate noise based on bottom portion of image
%   Threshold is based on noise STD
noise = max([mode(img(img~=0)) 0]);
if isempty(opt.th) || opt.th == 0
    opt.th = get_threshold(img(:,:,:,opt.th_te),3*noise);
end
clear noise

%   Compute mask
mask = img(:,:,:,opt.th_te) > opt.th;
mask = all(mask,4);
if opt.debug && ~(isempty(opt.th) || opt.th == 0)
    montage(reshape(mask,[opt.size(1) opt.size(2) 1 opt.size(3)]));
    title('Mask')
end

%   Initialize output variables
T2    = zeros(opt.size(1:3));
AMP = T2;

%   Fit a single voxel to get correct options structure, and plot if needed
S = img(1,1,1,:);
[~,~,opt] = ExpFit(S(:),opt);
opt.debug = 0;

%   Loop through voxels and fit
disp('Performing Stimulated Echo Compensated Fit');
t0 = tic;
parfor npi = 1:opt.size(1)  %   Parallel loop

    %   Some temporary variables needed for the parfor loop
    imgtmp = img(npi,:,:,:);
    opttmp = opt;
    t0tmp  = tic;
    masktmp = mask(npi,:,:);
    t2 = T2(npi,:,:,:);
    amp = AMP(npi,:,:,:);
    
    %   Loop through phase and slice points
    for nvi = 1:opttmp.size(2)
    for nsi = 1:opttmp.size(3)

        %   Fit if voxel is in mask
        if masktmp(1,nvi,nsi)
            
            %   Extract signal and fit
            S = imgtmp(1,nvi,nsi,:);
            [t2(1,nvi,nsi,:),amp(1,nvi,nsi,:)] = ExpFit(S(:),opttmp);
            
        end
    
    end
    end
    
    %   Save data in full variable
    T2(npi,:,:,:) = t2;
    AMP(npi,:,:,:) = amp;
    fprintf('Done point %3i (%g s)\n',npi,toc(t0tmp));
    
end
fprintf('\nDone all points (%g minutes)\n',toc(t0)/60);

%   Get R2
R2 = T2;
for i = 1:size(T2,4)
    R2tmp = R2(:,:,:,i);
    R2tmp(mask) = 1./R2tmp(mask);
    R2(:,:,:,i) = R2tmp;
end
clear R2tmp


%   Save Data
save('./ExpFit.mat','img','opt','T2','R2','AMP');
