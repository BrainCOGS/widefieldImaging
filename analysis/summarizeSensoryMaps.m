function sensoryMaps = summarizeSensoryMaps(visPath,airpuffPath,manualVisFlag)

% summarizeSensoryMaps(visPath,airpuffPath,manualVisFlag)
% summarizes and registers visual and airpuff maps to each other
% manualVisFlag will prompt user to manually define visual ROIs (default true)
if nargin < 3
  manualVisFlag = true; % select visual areas manually
end

%% get reference image, bregma, midline
mouseID   = mouseAndDateFromFileName(visPath);
spockFlag = isThisSpock;
if spockFlag
  refPath     = [widefieldParams.serverpath_spock mouseID '/'];
else
  if ispc
    refPath   = [widefieldParams.serverpath_pc mouseID '/'];
  else
    if isempty(strfind(visPath,'Users'))
      refPath = [widefieldParams.serverpath_mac mouseID '/'];
    else
      refPath = [widefieldParams.savepath_mac mouseID '/'];
    end
  end
end

load(formatFilePath([refPath 'refIm.mat'],false))
refim = meanproj; clear meanproj
cd(refPath)

%% for some mice automatic registration fails, here are manual correct factors
% [colshift rowshift angle scale], typically from xRecAlign_manual
switch mouseID
  case 'ai2'
    cf_whiskers = [170 2 -1 1];
    cf_visual   = [];
    
  case 'ai3'
    cf_whiskers = [173 0  2  1];
    cf_visual   = [243 7 -5  1];

  case 'ai8'
    cf_whiskers = [163 8 0 1];
    cf_visual   = [163 8 0 1];
    
  otherwise
    cf_whiskers = [];
    cf_visual   = [];
end

%% get airpuffs map and register
load(formatFilePath(airpuffPath,false),'coord','meanproj','mapSmthTh')
if ~isempty(strfind(airpuffPath,'whiskersR'))
  ctxSide   = 'L';
  otherSide = 'R';
else
  ctxSide   = 'R';
  otherSide = 'L';
end

% spatially bin if necessary
if max([size(meanproj,1) size(meanproj,2)]) > size(refim,1)
  refim_temp = imresize(refim,widefieldParams.dsFactor);
  binFlag    = true;
else
  refim_temp = refim;
  binFlag    = false;
end

% register to reference image
if isempty(cf_whiskers)
  [tform_whiskers,regim_whiskers] = registerRecs(meanproj,refim_temp,[],[],'rigid',false); 
else
  tform_whiskers = geoTransfFromPxl(cf_whiskers);
  im             = meanproj;
  im(isnan(im))  = 0;
  Rfixed         = imref2d(size(refim_temp));
  regim_whiskers = imwarp(im,tform_whiskers,'OutputView',Rfixed); 
end

[c,r]            = transformPointsForward(tform_whiskers,coord(:,2),coord(:,1));
r                = round(r);
c                = round(c);

% spatially bin if necessary
if binFlag
  temp = zeros(size(refim_temp));
  for iPxl = 1:numel(r)
    temp(r(iPxl),c(iPxl)) = 1;
  end
  temp     = imBinSpace(temp);
  [r,c]    = find(temp > 0);
end

% pick contralateral hemisphere only
[Ridx,Lidx]                     = assignPxlToHem([r c],midline);
if strcmpi(ctxSide,'L')
  regcoord = [r(Lidx) c(Lidx)];
else
  regcoord = [r(Ridx) c(Ridx)];
end

% make sure pick only contiguous region with higher pxl denisty
temp = zeros(size(refim));
for iPxl = 1:size(regcoord,1)
  temp(regcoord(iPxl,1),regcoord(iPxl,2)) = 1;
end
centroid = nanmedian(regcoord);
temp     = grayconnected(temp,centroid(1),centroid(2));
[r,c]    = find(temp > 0);
regcoord = [r c];

% fill in ROI coordinates and copy over to contralateral hemisphere
[nXr,nYr]                       = size(refim);
areaCoord{1}                    = regcoord;
areaLbl{1}                      = ['SS' lower(ctxSide)];
areaCoord{end+1}                = reflectCoordOverMidline(regcoord,midline,otherSide,[nXr nYr]);
areaLbl{end+1}                  = ['SS' lower(otherSide)];

% put eveything in data structure
sensoryMaps.refIm               = refim;
sensoryMaps.refImIsBinned       = binFlag;
sensoryMaps.midline             = midline;
sensoryMaps.bregma              = bregma;
sensoryMaps.whiskers.ctxSide    = ctxSide;
sensoryMaps.whiskers.meanproj   = meanproj;
sensoryMaps.whiskers.tform      = tform_whiskers;
sensoryMaps.whiskers.regIm      = regim_whiskers;
sensoryMaps.whiskers.respMap    = mapSmthTh;
sensoryMaps.whiskers.path       = airpuffPath;

clear meanproj coord temp r c regcoord

%% get labels and visual area coordinates
load(formatFilePath(visPath,false),'visualAreas','noData')

% first adjust scale (gets compressed during segmentation)
imgScale  = size(noData,1) / size(visualAreas.vasculature,1);
meanproj  = imresize(visualAreas.vasculature,imgScale);
fieldSign = imresize(visualAreas.fieldSign,imgScale);

if ~isempty(strfind(visPath,'eyeR'))
  ctxSide   = 'L';
  otherSide = 'R';
elseif ~isempty(strfind(visPath,'eyeL'))
  ctxSide   = 'R';
  otherSide = 'L';
else
  ctxSide   = 'bilateral';
  otherSide = 'none';
end

% spatially bin if necessary
if max([size(meanproj,1) size(meanproj,2)]) > size(refim,1)
  refim_temp = imresize(refim,widefieldParams.dsFactor);
  binFlag    = true;
else
  refim_temp = refim;
  binFlag    = false;
end

% register to reference image
if isempty(cf_visual)
  [tform_visual,regim_visual] = registerRecs(meanproj,refim_temp,[],[],'rigid',false); 
else
  tform_visual  = geoTransfFromPxl(cf_visual);
  im            = meanproj;
  im(isnan(im)) = 0;
  Rfixed        = imref2d(size(refim_temp));
  regim_visual  = imwarp(im,tform_visual,'OutputView',Rfixed); 
end

% get ROIs and transfer coordinates
if manualVisFlag
  if isempty(dir([refPath 'manualROIfieldSign.mat']))
    [ROI,ROIlbl] = manualROIfromFieldSign(fieldSign,true);
  else
    load([refPath 'manualROIfieldSign.mat'],'ROI','ROIlbl') 
  end
else
  ROIlbl       = visualAreas.name;
  nROI         = numel(ROIlbl);
  ROI          = cell(1,nROI);
  for iROI = 1:nROI
    [r,c]      = find(visualAreas.areas == iROI);
    ROI{iROI}  = [r c];
  end
end

% ROIs in ref image space
for iROI = 1:numel(ROI)
  [c,r]     = transformPointsForward(tform_visual,ROI{iROI}(:,2),ROI{iROI}(:,1));
  r         = round(r);
  c         = round(c);
  delidx    = r < 1 | r > size(refim_temp,1) | c < 1 | c > size(refim_temp,2);
  r(delidx) = [];
  c(delidx) = [];
    
  % spatially bin if necessary
  if binFlag
    temp = zeros(size(refim_temp));
    for iPxl = 1:numel(r)
      temp(r(iPxl),c(iPxl)) = 1;
    end
    temp     = imBinSpace(temp);
    [r,c]    = find(temp > 0); 
  end
  
  % pick contralateral hemisphere only?
  if strcmpi(ctxSide,'bilateral')
    areaCoord{end+1}  = [r c];
    areaLbl{end+1}    = ROIlbl{iROI};
  else
    [Ridx,Lidx]       = assignPxlToHem([r c],midline);
    if strcmpi(ctxSide,'L')
      regcoord        = [r(Lidx) c(Lidx)];
    else
      regcoord        = [r(Ridx) c(Ridx)];
    end
    
    % fill in ROI coordinates and copy over to contralateral hemisphere
    [nXr,nYr]         = size(refim);
    areaCoord{end+1}  = regcoord;
    areaLbl{end+1}    = [ROIlbl{iROI} lower(ctxSide)];
    areaCoord{end+1}  = reflectCoordOverMidline(regcoord,midline,otherSide,[nXr nYr]);
    areaLbl{end+1}    = [ROIlbl{iROI} lower(otherSide)];
  end
end

% put eveything in data structure
sensoryMaps.visual.ctxSide     = ctxSide;
sensoryMaps.visual.meanproj    = meanproj;
sensoryMaps.visual.tform       = tform_visual;
sensoryMaps.visual.regIm       = regim_visual;
sensoryMaps.visual.fieldSign   = fieldSign;
sensoryMaps.visual.isManualROI = manualVisFlag;
sensoryMaps.visual.path        = visPath;

% flipped field sign for ai2 and 3 due to early bug
switch mouseID
  case {'ai2','ai3'}
    sensoryMaps.visual.fieldSign = -sensoryMaps.visual.fieldSign;
end

%% plot and save
I       = zeros(size(refim));
outline = I;
for iArea = 1:numel(areaCoord)
  subI = zeros(size(I));
  for iPxl = 1:size(areaCoord{iArea},1)
    I(areaCoord{iArea}(iPxl,1),areaCoord{iArea}(iPxl,2))    = iArea;
    subI(areaCoord{iArea}(iPxl,1),areaCoord{iArea}(iPxl,2)) = 1;
  end
  pm      = bwperim(subI);
 try outline = outline + pm; catch; keyboard; end
end
outline(outline >= 1) = 1;
ROIimage = I;
clear subI pm I

overlay = repmat(refim./max(max(refim)),[1 1 3]);
[ii,jj] = find(outline == 1);
for iP = 1:numel(ii)
  overlay(ii(iP),jj(iP),:) = [1 1 1];
end

% data structure
sensoryMaps.areaCoord    = areaCoord;
sensoryMaps.areaCentroid = cellfun(@nanmean,areaCoord,'uniformOutput',false);
sensoryMaps.areaLbl      = areaLbl;
sensoryMaps.ROIimage     = ROIimage;
sensoryMaps.ROIoutline   = outline;
sensoryMaps.vascOverlay  = overlay;

fh = figure; 
image(overlay); axis image; axis off
for iArea = 1:numel(areaCoord)
  centroid = mean(areaCoord{iArea});
  text(centroid(2),centroid(1),areaLbl{iArea},...
       'color','r','horizontalAlignment','center','fontsize',10)
end

%%
saveas(fh,formatFilePath([refPath 'sensoryMaps.fig'],false))
clear centroid fh

save(formatFilePath([refPath 'sensoryMaps.mat'],false),'sensoryMaps')
