function dff = wideFieldPreProc(rec,cfg)

% dff = wideFieldPreProc(rec,cfg)
% applies full preprocessing pipeline to all .tifs in folder
% rec is full path to tif stacks, cfg is a strcuture with analysis
% parameters. refer to bottom of the m file for all fields
%
% Lucas Pinto (lpinto@princeton.edu)

%% defaults etc 
if nargin < 1 || isempty(rec)
  rec = pwd;
end
if nargin < 2
  cfg = struct([]);
end

try

% remove tankmousevr repo to avoid function name conflict  
if isThisSpock
  warning('off','all'); 
  rmpath(genpath('/usr/people/lpinto/code/tankmousevr/')); 
  warning('on','all'); 
end  

rec = formatFilePath(rec); % format rec path for this OS

cd(rec)
cfg = populateCfg(cfg); % get analysis config (nested fcn at bottom)
fl  = generateFilelist(rec); % stacks are typically in multiple files

%% motion correction
if cfg.motionCorrect
  mcorr   = getMotionCorrection(fl, false, false, 15, 5, false, 0.3, nan, 10);
end

%% extract time stamps from Tif headers
if isempty(dir('info.mat'))
  ts        = extractImageTimeStamps(fl);
  frameRate = 1/mode(diff(ts));
  recDur    = ts(end) - ts(1);
else
  load info ts frameRate recDur cfg
end
if ~isfield(cfg,'binTimeTo')
  cfg(1).binTimeTo = 20;
end

%% just load first movie segment for vascular detection and mean Z projection
if cfg.motionCorrect
  imstack   = cv.imreadx(fl{1}, mcorr(1).xShifts, mcorr(1).yShifts); % load image shifts and apply them
else
  imstack   = loadTiffStack(fl{1}); % load stack
end

% calculate mean projection just on blue LED frames
nZ          = size(imstack,3); 
if cfg.strobedLED
  % estimate LED illumination sequence from image contrast
  if isempty(cfg.strobeSeq)
    cfg.strobeSeq = estimateStrobeSeq(imstack(:,:,1:100));
  end
    
  strobeSeq = repmat(cfg.strobeSeq,[1 floor(nZ/numel(cfg.strobeSeq))]);
  excFrames = nZ - numel(strobeSeq);
  strobeSeq = [strobeSeq cfg.strobeSeq(1:excFrames)];
  idx       = strobeSeq == 1;
else
  idx       = 1:size(imstack,3);
end

meanproj    = mean(double(imstack(:,:,idx)),3); % Z projection
vascMask    = maskVasculature(meanproj, rec, cfg); % remove vasculature

save info frameRate rec fl cfg recDur ts

%% plot strobe seq as sanity check
if ~ispc
figure;
strobeSeq  = repmat(cfg.strobeSeq,[1 10]);
for iFrame = 1:20
  subplot(4,5,iFrame)
  imshow(imstack(:,:,iFrame)); axis image; axis off
  text(10,20,num2str(strobeSeq(iFrame)),'color','y','fontsize',14,'fontweight','bold')
end
set(gcf,'position',[10 10 1000 790])
export_fig strobeCheck.pdf; 
end
close
  
%% extract rawf (spatially bin and NaN masked pixels) 
rawfTemp    = cell(numel(fl),1);
rawfTemp{1} = extractF(double(imstack),vascMask);

%% now loop through other files and do the same
if cfg.parallelize
  parpool;
  parfor iFile = 2:numel(fl)
    if cfg.motionCorrect
      imstack = cv.imreadx(fl{iFile}, mcorr(iFile).xShifts, mcorr(iFile).yShifts); % load image shifts and apply them
    else
      imstack = loadTiffStack(fl{iFile}); % load stack
    end
    rawfTemp{iFile} = extractF(double(imstack),vascMask);
  end
  poolobj = gcp('nocreate');
  delete(poolobj);
else
  for iFile = 2:numel(fl)
    if cfg.motionCorrect
      imstack = cv.imreadx(fl{iFile}, mcorr(iFile).xShifts, mcorr(iFile).yShifts); % load image shifts and apply them
    else
      imstack = loadTiffStack(fl{iFile}); % load stack
    end
    rawfTemp{iFile} = extractF(double(imstack),vascMask);
    clear imstack
  end
end

%% make a single image stack of raw fluorescence
nFramesFl   = cellfun(@(x)(size(x,3)),rawfTemp);
nFrames     = sum(nFramesFl);
rawf        = zeros(size(rawfTemp{1},1),size(rawfTemp{1},2),nFrames);
fcount      = 0;
for iFile = 1:numel(rawfTemp)
  rawf(:,:,fcount+1:fcount+nFramesFl(iFile)) = rawfTemp{iFile};
  fcount                                     = fcount+nFramesFl(iFile);
end
clear rawfTemp fcount

% if there are nans in any frame, nan all frames (because motion correction)
rawf        = evenNan(rawf);

%% for VR widefield typically we want to rotate image / flip axis
if cfg.rotateImage
  rawf      = fliplr(rot90(rawf,3));
end

%% bin time if necessary
if cfg.binTimeTo < round(frameRate)
  timeBinFactor = round(frameRate)/cfg.binTimeTo;
  [nX,nY,~]     = size(rawf);
  if cfg.strobedLED
    [rf,bc]     = conditionDffMat(rawf);
    binB        = bintime(rf(strobeSeq==1,:)',[],timeBinFactor,'sum')';
    binV        = bintime(rf(strobeSeq==2,:)',[],timeBinFactor,'sum')';
    binR        = zeros(size(binV,1)+size(binB,1),size(binB,2));
    if strobeSeq(1) == 1
      binR(1:2:end,:) = binB;
      binR(2:2:end,:) = binV;
    else
      binR(1:2:end,:) = binV;
      binR(2:2:end,:) = binB;
    end
  else
    [rf,bc]     = conditionDffMat(rawf);
    binR        = bintime(rf',[],timeBinFactor,'sum')';
  end
  rawf          = conditionDffMat(binR,bc,[],[nX nY size(binR,1)]);
else
  timeBinFactor = 1;
end

%% save info
save info timeBinFactor nFrames -append

%% calculate dff and save (with or without hemodynamic correction)
if cfg.strobedLED
  [dff,rawfBlue,rawfViolet,dffBlue,dffViolet,cf] = dffFrom2Channels(rawf,frameRate,cfg); %#ok<ASGLU>
  
  tic; fprintf('saving... ')
  save('rawf','rawf','rawfBlue','rawfViolet','dffBlue','dffViolet','cf','-v7.3')
  save('dff','dff','vascMask','meanproj','-v7.3')
  fprintf('done after %1.1f min\n',toc/60)
else
  F0        = computeF0mode(rawf,cfg.winSizeModeSec,frameRate,true); % F0 is a runnig mode
  dff       = (rawf-F0)./F0;
  
  tic; fprintf('saving... ')
  save('rawf','rawf','F0','-v7.3')
  save('dff','dff','vascMask','meanproj','-v7.3')
  fprintf('done after %1.1f min\n',toc/60)
end

if isThisSpock; addpath(genpath('/usr/people/lpinto/code/tankmousevr/')); end  

catch err
  displayException(err);
end

end

%% --------------------------------------------------------------------
%%
function cfg = populateCfg(cfg)

%% PREPROCESSING PARAMS
if ~isfield(cfg,'motionCorrect')
  cfg(1).motionCorrect   = true;
end
if ~isfield(cfg,'parallelize')
  cfg(1).parallelize     = true;
end
if ~isfield(cfg,'strobedLED')
  cfg(1).strobedLED      = true;
end
if ~isfield(cfg,'rotateImage')
  cfg(1).rotateImage     = true;
end

%% 2-LED dff correction params
if ~isfield(cfg,'winSizeSec')
  cfg(1).winSizeSec            = 30;
end
if ~isfield(cfg,'dffMethod')
  cfg(1).dffMethod             = 'rawfFitDivide'; 
end
if ~isfield(cfg,'runWinMethod')
  cfg(1).runWinMethod          = 'overall'; 
end
if ~isfield(cfg,'smoothBaseline')
  cfg(1).smoothBaseline        = true;
end
if ~isfield(cfg,'smoothDff')
  cfg(1).smoothDff             = false;
end
if ~isfield(cfg,'smoothSigma')
  cfg(1).smoothSigma           = 3;
end
if ~isfield(cfg,'strobeSeq')
  cfg(1).strobeSeq             = [ ];
end
if ~isfield(cfg,'binTimeTo')
  cfg(1).binTimeTo             = 20;
end

%% VASCULAR DETECTION PARAMS
% ALGORITHM PARAMETERS
mid = mouseAndDateFromFileName(pwd); % vascular detection requires mouse-specific parameters

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
  cfg(1).baseScale        = linspace(0.03, .4, 30);%[linspace(0.01, .1, 10) linspace(0.15, 1, 10)]; % controls scales for filtering
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
  switch mid
    case {'ai2','ai3'}
      cfg(1).minNSigmas   = 1.8;
    case 'ai5'
      cfg(1).minNSigmas   = 1.2;
    otherwise
      cfg(1).minNSigmas   = 1.6;
  end
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
  switch mid
    case {'ai7','ai2'}
      cfg(1).subSearchPxls{1} = {200:300,1:210}; % frontal part of the sinus
      cfg(1).subSearchPxls{2} = {220:280,211:512}; % posterior part of the sinus
    otherwise
      cfg(1).subSearchPxls{1} = {180:320,1:230}; % frontal part of the sinus
      cfg(1).subSearchPxls{2} = {210:270,231:512}; % posterior part of the sinus
  end
else
  if ~iscell(cfg.subSearchPxls)
    cfg(1).subSearchPxls  = {cfg(1).subSearchPxls};
  end
end
if ~isfield(cfg,'subSearchMinSigmas') % separate image for special searches
  switch mid
    case 'ai2'
      cfg(1).subSearchMinSigmas = [1; 1]; % posterior part of the sinus
    case 'ai3'
      cfg(1).subSearchMinSigmas = [.75; .75]; % posterior part of the sinus
    case 'ai5'
      cfg(1).subSearchMinSigmas = [1.5; 1.5]; % posterior part of the sinus
    otherwise
      cfg(1).subSearchMinSigmas = [1.5; 1]; % posterior part of the sinus
  end
  
end
if ~isfield(cfg,'subSearchMaxIter') % separate image for special searches
  switch mid
    case 'ai2'
      cfg(1).subSearchMaxIter = [1; 2];  % posterior part of the sinus
    otherwise
      cfg(1).subSearchMaxIter = [1; 1];  % posterior part of the sinus
  end
  
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