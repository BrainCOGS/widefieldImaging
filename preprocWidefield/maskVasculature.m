function mask = maskVasculature(meanproj, rec, cfg)

% mask = maskVasculature(meanproj, rec, cfg)
%
% creates a mask for pixels corresponding to detected vasculature using an
% adaptive median filtering algorithm
%   output mask is a logical matrix indicating vasculature pixels
%   meanproj is single image frame (typically mean Z projection)
%   rec is recording name (default current dir)
%   cfg is a structure containing analysis parameters, may contain any
%   subset of the following fields (missing fields will be added)
%
%         ALGORITHM PARAMETERS
%         cfg.imgScale         = 1;
%         cfg.maxIterations    = 1; *
%         cfg.globalBaseline   = false;
%         cfg.baseScale        = linspace(0.01, 1, 50); * % controls scales for filtering
%         cfg.pixBins          = linspace(-0.4, 0.3, 500);
%         cfg.pixMedian        = true(1, 5);
%         cfg.minFracPixels    = 0.5;
%         
%         PIXEL DISTRIBUTION PARAMETERS (CUTOFFS ETC)
%         cfg.maxQuantile      = 0.95; 
%         cfg.peakFrac         = 0.15;
%         cfg.minNSigmas       = 1.5; * 
%         cfg.maxNSigmas       = inf; *
%         cfg.blackTol         = 1-1e-4;
%
%         SPECIAL CASE FOR SAGITTAL SINUS / IMAGE CORNERS: ATTENTION:
%         DEFAULTS ASSUME 512 x 512 IMAGE
%         cfg.separateImCenter = true; *
%         cfg.subSearchPxls    = {180:320,1:230};{210:260,231:512}; * % do subsearches in
%                                                                       these {X,Y} image regions
%         cfg.subSearchBaseScale
%         cfg.subSearchMinSigmas
%         cfg.subSearchMaxSigmas
%         cfg.subSearchMaxIter
%         cfg.maxQuantile
%         cfg.subSearchFillInPxls = true; % fill mask (midline pixels)
%
%         LOAD OVAL MASK FROM DISK
%         cfg.removeOffBrain   = true; *
%
%         * denotes most typical fields to tweak for different datasets
%
% Lucas Pinto (lpinto@princeton.edu), modified from Sue Ann Koay

%% handle defaults and generate file name
if nargin < 2 || isempty(rec)
  rec      = pwd;
end
if nargin < 3 
  cfg      = struct([]);
end

% nested function to fill structure with all parameters
cfg        = populateCfg(cfg); 

% convert to rec to full path if necessary
if isempty(strfind(rec,'/')) && isempty(strfind(rec,'\'))
  rec      = [widefieldParams.savepath rec];
end
rec        = formatFilePath(rec,true);
outputFile = [rec 'vasculature'];


%% first run the regular code
fig             = gobjects(0);
cropFile        = '';
blackRemoval    = {cfg.blackTol, true};

% Set cropped pixels and extreme values to nan
data.vasculature  = meanproj;
[data.isCropped, data.moreCropped, data.lessCropped]    ...
                  = getCroppedPixels(cropFile, size(data.vasculature), cfg.imgScale, 0.9, 1.3);
data.maxF         = quantile(colvec(data.vasculature(~data.moreCropped)), cfg.maxQuantile);
data.outsideF     = quantile(colvec(data.vasculature( data.lessCropped)), cfg.maxQuantile);
if data.outsideF > data.maxF
  data.vasculature(data.isCropped | (data.moreCropped & data.vasculature > data.maxF))  = nan;
else
  data.vasculature(data.isCropped)  = nan;
end


% Assume that most of the image, i.e. the mode of the pixel value PDF, is brain and not
% vasculature. Run this procedure iteratively.
[masked, fig]   = maskOutliers(data.vasculature, cfg, 'vasculature', [1e3 6e3], fig);
data.noData     = masked(:,:,end);
mask1           = data.noData;
close(fig)

%% now just run center pxls corresponding to sinuses and/or image corners with special
% parameters
if cfg.separateImCenter

  fig             = gobjects(0);
  cropFile        = '';
  blackRemoval    = {cfg.blackTol, true};
  
  mask2           = zeros(size(meanproj));
  
  for iSub = 1:numel(cfg.subSearchPxls) 
    
    % this subregion
    centerPxlsX       = cfg.subSearchPxls{iSub}{1};
    centerPxlsY       = cfg.subSearchPxls{iSub}{2};
    cfg.minNSigmas    = cfg.subSearchMinSigmas(iSub);
    cfg.maxIterations = cfg.subSearchMaxIter(iSub);
    cfg.baseScale     = cfg.subSearchBaseScale{iSub};
    cfg.maxQuantile   = cfg.subSearchMaxQuantile(iSub);
    cfg.maxNSigmas    = cfg.subSearchMaxNSigmas(iSub);
    
    % Set cropped pixels and extreme values to nan
    data.vasculature  = meanproj(centerPxlsX,centerPxlsY);
    [data.isCropped, data.moreCropped, data.lessCropped]    ...
                      = getCroppedPixels(cropFile, size(data.vasculature), cfg.imgScale, 0.9, 1.3);
    data.maxF         = quantile(colvec(data.vasculature(~data.moreCropped)), cfg.maxQuantile);
    data.outsideF     = quantile(colvec(data.vasculature( data.lessCropped)), cfg.maxQuantile);
    if data.outsideF > data.maxF
      data.vasculature(data.isCropped | (data.moreCropped & data.vasculature > data.maxF))  = nan;
    else
      data.vasculature(data.isCropped)  = nan;
    end
    data.vasculature(data.isCropped)  = nan;
    % Assume that most of the image, i.e. the mode of the pixel value PDF, is brain and not
    % vasculature. Run this procedure iteratively.
    [masked, fig]                  = maskOutliers_LP(data.vasculature, cfg, 'vasculature', [1e3 6e6], fig);
    data.noData                    = masked(:,:,end);
    mask2(centerPxlsX,centerPxlsY) = data.noData;

  end
  
  % mask within-sinus pixels
  if cfg.subSearchFillInPxls
    se    = strel('line',10,6);
    mask2 = imdilate(mask2,se);
    se1   = strel('ball',4,4);
    for ii = 1:4; mask2 = imdilate(mask2,se1); end
    mask2 = mask2 > min(min(mask2));
  end
  
  % combine 2 masks
  mask = mask1+mask2;
  mask = mask>0;
  close(fig)
  save(outputFile, 'mask1', 'mask2')
  
else
  mask = mask1;
  save(outputFile, 'mask1')
end

if cfg.removeOffBrain
  if ~isempty(strfind(rec,'Volumes'))
    load(widefieldParams.brainMaskPath,'brainMask')
  else
    if ispc
      load(widefieldParams.brainMaskPath_pc,'brainMask')
    else
      load(widefieldParams.brainMaskPath_spock,'brainMask')
    end
  end
  mask = (mask + brainMask) > 0;
end

% Write output to disk
save(outputFile, '-struct', 'data', '-append');
save(outputFile, 'mask','-append');
fprintf(' ***  Saved to %s\n', outputFile);

% figure for the records
figure; 
subplot(1,2,1), imshow(uint8(meanproj)), 
subplot(1,2,2), imagesc(mask), colormap gray, axis off, axis image
saveas(gcf,'vasculature')
close

% Write additional copy and images
%   if ~exist(outputPath, 'dir')
%     mkdir(outputPath);
%   end
%   [~,name,ext]      = parsePath(outputFile);
%   copyfile(outputFile, fullfile(outputPath, [name ext]));

%   for iFig = 1:numel(fig)
%     figFile         = fullfile(outputPath, get(fig(iFig),'Name'));
%     saveas(fig(iFig), figFile)
%     close(fig(iFig))
% %     export_fig(fig(iFig), [figFile '.png']);
%   end

end

function cfg = populateCfg(cfg)

% ALGORITHM PARAMETERS
if ~isfield(cfg,'imgScale')
  cfg(1).imgScale         = 1;
end
if ~isfield(cfg,'maxIterations')
  cfg(1).maxIterations    = 1;
end
if ~isfield(cfg,'globalBaseline')
  cfg(1).globalBaseline   = false;
end
if ~isfield(cfg,'baseScale')
  cfg(1).baseScale        = linspace(0.01, .4, 30);%[linspace(0.01, .1, 10) linspace(0.15, 1, 10)]; % controls scales for filtering
end
if ~isfield(cfg,'pixBins')
  cfg(1).pixBins          = linspace(-0.4, 0.3, 500);
end
if ~isfield(cfg,'pixMedian')
  cfg(1).pixMedian        = true(1, 5);
end
if ~isfield(cfg,'minFracPixels')
  cfg(1).minFracPixels    = 0.5;
end

% PIXEL DISTRIBUTION PARAMETERS (CUTOFFS ETC)
if ~isfield(cfg,'maxQuantile')
  cfg(1).maxQuantile      = 0.95;
end
if ~isfield(cfg,'peakFrac')
  cfg(1).peakFrac         = 0.15;
end
if ~isfield(cfg,'minNSigmas')
  cfg(1).minNSigmas       = 1.5;
end
if ~isfield(cfg,'maxNSigmas')
  cfg(1).maxNSigmas       = 10;
end
if ~isfield(cfg,'blackTol')
  cfg(1).blackTol         = 1-1e-5;
end

% SPECIAL CASE FOR SAGITTAL SINUS / IMAGE CORNERS
if ~isfield(cfg,'separateImCenter')
  cfg(1).separateImCenter = true;
end
if ~isfield(cfg,'subSearchPxls') % separate image for special searches
  cfg(1).subSearchPxls{1} = {180:320,1:230}; % frontal part of the sinus
  cfg(1).subSearchPxls{2} = {210:270,231:512}; % posterior part of the sinus
%   % corners
%   cfg(1).subSearchPxls{3} = {1:512,1:512};
%   cfg(1).subSearchPxls{4} = {413:512,1:512};
%   cfg(1).subSearchPxls{5} = {413:512,1:100};
%   cfg(1).subSearchPxls{6} = {413:512,413:512};
else
  if ~iscell(cfg.subSearchPxls)
    cfg(1).subSearchPxls  = {cfg(1).subSearchPxls};
  end
end
if ~isfield(cfg,'subSearchMinSigmas') % separate image for special searches
  cfg(1).subSearchMinSigmas = [0.5; 5]; % posterior part of the sinus
end
if ~isfield(cfg,'subSearchMaxIter') % separate image for special searches
  cfg(1).subSearchMaxIter = [1; 1];  % posterior part of the sinus
end
if ~isfield(cfg,'subSearchMaxQuantile') % separate image for special searches
  cfg(1).subSearchMaxQuantile = [1; 1];  % posterior part of the sinus
end
if ~isfield(cfg,'subSearchMaxNSigmas') % separate image for special searches
  cfg(1).subSearchMaxNSigmas = [inf; inf];  % posterior part of the sinus
end
if ~isfield(cfg,'subSearchBaseScale') % separate image for special searches
  cfg(1).subSearchBaseScale{1} = linspace(0.15, 1, 10); % frontal part of the sinus
  cfg(1).subSearchBaseScale{2} = linspace(0.01, .5, 15); % posterior part of the sinus
%   cfg(1).subSearchBaseScale{3} = linspace(0.8, 1, 4); % frontal part of the sinus
%   cfg(1).subSearchBaseScale{4} = linspace(0.8, 1, 4); % frontal part of the sinus
%   cfg(1).subSearchBaseScale{5} = linspace(0.5, 1, 10); % frontal part of the sinus
%   cfg(1).subSearchBaseScale{6} = linspace(0.5, 1, 10); % frontal part of the sinus

else
  if ~iscell(cfg.subSearchBaseScale)
    cfg(1).subSearchBaseScale  = {cfg(1).subSearchBaseScale};
  end
end
if ~isfield(cfg,'subSearchFillInPxls')
  cfg(1).subSearchFillInPxls   = true;
end
if ~isfield(cfg,'removeOffBrain')
  cfg(1).removeOffBrain   = true;
end

end
