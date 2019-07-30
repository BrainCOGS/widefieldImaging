function [dffROI,ROIlbl] = extractDffFromROI(rec,cfg,frameRate,ROIpath)

% dffROI = extractDffFromROI(rec,cfg,frameRate,ROIpath)
% extracts ROI-averaged dff by averaging F first and recalculating
% INPUT
%   rec is recording path
%   cfg is optional analysis params structure (loaded from info.mat by default)
%   frameRate is recording frame rate (loaded from info.mat by default)
%   ROI path is path for refROI.mat (assumed to be mouse root directory by default)
%
% OUTPUT
%   dffROI output is a frames x ROI matrix with dff traces
%   ROIlbl is cell array with corresponding ROI names (strings)

%%
if nargin < 1; rec = pwd;      end
if nargin < 2; cfg = [];       end
if nargin < 2; frameRate = []; end
if nargin < 4; ROIpath = [];   end

rec = formatFilePath(rec);
cd(rec)
if isempty(cfg);       load info cfg;       end
if isempty(frameRate); load info frameRate; end

%% get registered ROIs
[tform,~,~,~,meanproj] = registerRecs([],[],[],[],[],[],'monomodal'); close
[ROI,ROIlbl]           = getRefROI(ROIpath,tform,meanproj,true);

%% if same ROIs as saved, load and return
if ~isempty(dir('dffROI.mat'))
  thislbl = ROIlbl;
  load dffROI ROIlbl
  if numel(ROIlbl) == numel(thislbl)
    if sum(strcmpi(ROIlbl,thislbl)) == numel(ROIlbl)
      load dffROI dffROI
      return
    else
      ROIlbl = thislbl;
    end
  else
    ROIlbl = thislbl;
  end
end

%% avg rawf first, then calculate dff
tic
fprintf('extracting DFF for ROIs...')

%% calculate dff and save (with or without hemodynamic correction)
if cfg.strobedLED
  load rawf rawfViolet rawfBlue
  
  nZb      = size(rawfBlue,3);
  nZv      = size(rawfViolet,3);
  nZ       = min([nZb nZv]);
  nROI     = numel(ROI);
  rawfROIb = zeros(nZ,nROI); 
  rawfROIv = zeros(nZ,nROI); 

  for iROI = 1:nROI
    temp             = arrayfun(@(x,y)(squeeze(rawfBlue(x,y,1:nZ))),ROI{iROI}(:,1),ROI{iROI}(:,2),'UniformOutput',false);
    rawfROIb(:,iROI) = nanmean(cell2mat(temp'),2);
    temp             = arrayfun(@(x,y)(squeeze(rawfViolet(x,y,1:nZ))),ROI{iROI}(:,1),ROI{iROI}(:,2),'UniformOutput',false);
    rawfROIv(:,iROI) = nanmean(cell2mat(temp'),2);
  end
  
  % F0 is a runnig mode
  F0      = computeF0mode(rawfROIb,widefieldParams.winSizeModeSec,frameRate/2,false); 
  Rb      = rawfROIb./F0;
  F0      = computeF0mode(rawfROIv,widefieldParams.winSizeModeSec,frameRate/2,false); 
  Rv      = rawfROIv./F0;
  dffROI  = (Rb ./ Rv) - 1;
  
else 
  load rawf rawf
  nZ      = size(rawf,3);
  nROI    = numel(ROI);
  rawfROI = zeros(nZ,nROI); 

  for iROI = 1:nROI
    temp            = arrayfun(@(x,y)(squeeze(rawf(x,y,:))),ROI{iROI}(:,1),ROI{iROI}(:,2),'UniformOutput',false);
    rawfROI(:,iROI) = nanmean(cell2mat(temp'),2);
  end

  % F0 is a runnig mode
  F0      = computeF0mode(rawfROI,widefieldParams.winSizeModeSec,frameRate/2,false); 
  dffROI  = (rawfROI-F0)./F0;
end

save dffROI dffROI ROIlbl

fprintf(' done after %1.1f minutes\n',toc/60)