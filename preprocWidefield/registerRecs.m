function [tform,regim,bregma,midline] = registerRecs(im,refim,refBregma,refMidline,type,saveFlag,regoptimizer)

% [tform,regim,bregma,midline] = registerRecs(im,refim,refBregma,refMidline,type,saveFlag,regoptimizer)
% registers im to refim using rigid transformation (including rotations)
% if refBregma and refMidline are provided, they will be transformed too;
% otherwise user should pass empty arrays (if not specifically passed will
% try to load from disk)
%
% INPUT
%   im: image to be registered, loaded by default
%   refim: reference image, loaded by default
%   refBregma: bregma data structure from setBregmaAndMidline(), loaded by default
%   refMidline: midline data structure from setBregmaAndMidline(), loaded by default
%   type: imregister input, default 'rigid'
%   saveFlag: true to save registration (default true)
%   regoptimizer: imregister input, default 'multimodal';
%
% OUTPUT
%    tform: imregister output
%    regim: registered image
%    bregma: registered bregma coordinates
%    midline: registered midline coordinates
%
% LP sep 2016

%% defaults
if nargin < 1 || isempty(im)
  load dff meanproj
  im = meanproj;
  clear meanproj
end
if nargin < 2 || isempty(refim)
  mouseID = mouseAndDateFromFileName(pwd);
  wf      = widefieldParams;
  root    = formatFilePath(wf.getRootDir(isThisSpock));
  load([root mouseID '/refIm.mat'],'meanproj')
  refim   = meanproj;
  clear meanproj
end
if nargin < 3 || isempty(refBregma) || isempty(refMidline) 
  mouseID = mouseAndDateFromFileName(pwd);
  wf      = widefieldParams;
  root    = formatFilePath(wf.getRootDir(isThisSpock));
  load([root mouseID '/refIm.mat'],'bregma','midline')
  refBregma  = bregma;
  refMidline = midline;
  clear bregma midline
end
if nargin < 5 || isempty(type)
  type    = 'rigid';
end
if nargin < 6 || isempty(saveFlag)
  saveFlag = true;
end
if nargin < 7 || isempty(regoptimizer)
  regoptimizer = 'multimodal';
end

%% make sure im and refim are the same size, otherwise downsample and rotate
if size(im,1) ~= size(refim,1) && size(im,1) == size(im,2)
  sf = size(refim,1)/size(im,1);
  im = fliplr(rot90(imresize(im,sf),3));
end
refim(isnan(refim)) = 0;
im(isnan(im))       = 0;

% 
%% normalize images
im    = im./max(im(:));
refim = refim./max(refim(:));


%% register images
[optimizer,metric]          = imregconfig(regoptimizer);
regim                       = imregister(im,refim,type,optimizer,metric);
tform                       = imregtform(im,refim,type,optimizer,metric);

%% transform reference points (bregma and midline)
if ~isempty(refBregma)
  [x,y]  = transformPointsInverse(tform,refBregma(1),refBregma(2));
  bregma = round([x y]);
end

if ~isempty(refMidline)
  [x,y] = transformPointsInverse(tform,refMidline.p(1,1),refMidline.p(1,2));
  midline.p(1,:) = [x y];
  [x,y] = transformPointsInverse(tform,refMidline.p(2,1),refMidline.p(2,2));
  midline.p(2,:) = [x y];
  
  % as in y = ax + b
  midline.a = (midline.p(1,2)-midline.p(2,2))/(midline.p(1,1)-midline.p(2,1));
  midline.b = midline.p(1,2) - midline.a*midline.p(1,1);
end

%% plot and save
f1 = figure;
subplot(1,3,1); imagesc(refim); axis image; colormap gray; axis off; title('reference')
if ~isempty(refBregma); hold on; plot(refBregma(2),refBregma(1),'+'); end
subplot(1,3,2); imagesc(im); axis image; colormap gray; axis off; title('image')
if ~isempty(refBregma); hold on; plot(bregma(2),bregma(1),'+'); end
subplot(1,3,3); imshowpair(refim,regim); axis image; axis off; title('reference + reg. image')
set(f1,'position',[100 200 700 250])

%%
if saveFlag; saveas(f1,'imreg'); end

clear f1 refIm refMidline refBregma sf x y
if saveFlag; save reg2ref; end