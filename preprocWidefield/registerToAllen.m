function refROI = registerToAllen(mouseID,cfg)

% refROI = registerToAllen(mouseID,cfg)
% returns and saves reference ROI information registered to mouse mouseID
% cfg is optional analysis parameters structure
% requires https://github.com/BaselLaserMouse/AllenBrainAPI.git (Allen API interface)

%% defaults
if nargin < 1; mouseID = [];         end
if nargin < 2; cfg     = struct([]); end

%% initialize
if isempty(mouseID); mouseID = mouseAndDateFromFileName(pwd); end
wf   = widefieldParams;
root = wf.getRootDir(isThisSpock);

cd([root mouseID])

%%
cfg  = populateCfg(cfg);

%% make an image with right V1, bregma, approximate lambda as fiducial points
load sensoryMaps sensoryMaps
mouseIm  = zeros(size(sensoryMaps.refIm));
pxlPerMM = widefieldParams.pxlPerMM/widefieldParams.dsFactor;
bregma   = round(sensoryMaps.bregma);
lambdar  = round(bregma(1)+pxlPerMM*cfg.lambdaDist); % approximate lambda
if isfinite(sensoryMaps.midline.a)
  lambdac= (lambdar-sensoryMaps.midline.b)/sensoryMaps.midline.a;
else
  lambdac=bregma(2);
end

if cfg.includeBregmaAndLambda
  rowrange                   = bregma(1)-cfg.nPxlPoints:bregma(1)+cfg.nPxlPoints;
  colrange                   = bregma(2)-cfg.nPxlPoints:bregma(2)+cfg.nPxlPoints;
  mouseIm(rowrange,colrange) = 1;
  rowrange                   = lambdar-cfg.nPxlPoints:lambdar+cfg.nPxlPoints;
  colrange                   = lambdac-cfg.nPxlPoints:lambdac+cfg.nPxlPoints;
end

v1r      = strcmpi(sensoryMaps.areaLbl,'V1r');
v1r      = sensoryMaps.areaCoord{v1r};
v1l      = strcmpi(sensoryMaps.areaLbl,'V1l');
v1l      = sensoryMaps.areaCoord{v1l};
switch cfg.V1type
  
  case 'whole'
    for iPxl = 1:size(v1r,1)
      mouseIm(v1r(iPxl,1),v1r(iPxl,2)) = 1;
    end
    for iPxl = 1:size(v1l,1)
      mouseIm(v1l(iPxl,1),v1l(iPxl,2)) = 1;
    end
  case 'tip'
    subROI    = v1r((v1r(:,1) <= prctile(v1r(:,1),cfg.tipPrctile)),:);
    centr     = round(mean(subROI));
    rowrange  = centr(1)-cfg.nPxlPoints:centr(1)+cfg.nPxlPoints;
    colrange  = centr(2)-cfg.nPxlPoints:centr(2)+cfg.nPxlPoints;
    mouseIm(rowrange,colrange) = 1;
    
    subROI    = v1l((v1l(:,1) <= prctile(v1r(:,1),cfg.tipPrctile)),:);
    centr     = round(mean(subROI));
    rowrange  = centr(1)-cfg.nPxlPoints:centr(1)+cfg.nPxlPoints;
    colrange  = centr(2)-cfg.nPxlPoints:centr(2)+cfg.nPxlPoints;
    mouseIm(rowrange,colrange) = 1;
  case 'center'
    centr     = round(mean(v1r));
    rowrange  = centr(1)-cfg.nPxlPoints:centr(1)+cfg.nPxlPoints;
    colrange  = centr(2)-cfg.nPxlPoints:centr(2)+cfg.nPxlPoints;
    mouseIm(rowrange,colrange) = 1;
    
    centr     = round(mean(v1l));
    rowrange  = centr(1)-cfg.nPxlPoints:centr(1)+cfg.nPxlPoints;
    colrange  = centr(2)-cfg.nPxlPoints:centr(2)+cfg.nPxlPoints;
    mouseIm(rowrange,colrange) = 1;
end

%% do same for allen atlas
temp        = refROIfromAllen;
allenROI    = resizeAllenROI(temp,cfg);

%% allen refernce image
allenIm     = zeros(size(allenROI.refFig));

if cfg.includeBregmaAndLambda
  bregma      = round(allenROI.bregma);
  lambda      = bregma; 
  lambda(1)   = round(bregma(1)+(allenROI.pxlPerMM)*cfg.lambdaDist); % approximate lambda

  rowrange                   = bregma(1)-cfg.nPxlPoints:bregma(1)+cfg.nPxlPoints;
  colrange                   = bregma(2)-cfg.nPxlPoints:bregma(2)+cfg.nPxlPoints;
  allenIm(rowrange,colrange) = 1;
  rowrange                   = lambda(1)-cfg.nPxlPoints:lambda(1)+cfg.nPxlPoints;
  colrange                   = lambda(2)-cfg.nPxlPoints:lambda(2)+cfg.nPxlPoints;
end

v1          = strcmpi(allenROI.areaLbl,'VISp');
v1          = allenROI.ROI{v1};
v1r         = v1;
v1l         = v1;
v1r(v1r(:,2) < allenROI.bregma(2),:) = [];
v1l(v1l(:,2) > allenROI.bregma(2),:) = [];

switch cfg.V1type
  case 'whole'
    for iPxl = 1:size(v1r,1)
      allenIm(v1r(iPxl,1),v1r(iPxl,2)) = 1;
    end
    for iPxl = 1:size(v1l,1)
      allenIm(v1l(iPxl,1),v1l(iPxl,2)) = 1;
    end
  case 'tip'
    subROI    = v1r((v1r(:,1) <= prctile(v1r(:,1),cfg.tipPrctile)),:);
    centr     = round(mean(subROI));
    rowrange  = centr(1)-cfg.nPxlPoints:centr(1)+cfg.nPxlPoints;
    colrange  = centr(2)-cfg.nPxlPoints:centr(2)+cfg.nPxlPoints;
    allenIm(rowrange,colrange) = 1;
    
    subROI    = v1l((v1l(:,1) <= prctile(v1r(:,1),cfg.tipPrctile)),:);
    centr     = round(mean(subROI));
    rowrange  = centr(1)-cfg.nPxlPoints:centr(1)+cfg.nPxlPoints;
    colrange  = centr(2)-cfg.nPxlPoints:centr(2)+cfg.nPxlPoints;
    allenIm(rowrange,colrange) = 1;
  case 'center'
    centr     = round(mean(v1r));
    rowrange  = centr(1)-cfg.nPxlPoints:centr(1)+cfg.nPxlPoints;
    colrange  = centr(2)-cfg.nPxlPoints:centr(2)+cfg.nPxlPoints;
    allenIm(rowrange,colrange) = 1;
    
    centr     = round(mean(v1l));
    rowrange  = centr(1)-cfg.nPxlPoints:centr(1)+cfg.nPxlPoints;
    colrange  = centr(2)-cfg.nPxlPoints:centr(2)+cfg.nPxlPoints;
    allenIm(rowrange,colrange) = 1;
end

%% register, compile info
[tform,regim,br,md]      = registerRecs(allenIm,mouseIm,allenROI.bregma,allenROI.midline,cfg.transformType,false); close
allenROI.bregma          = br;
allenROI.midline         = md;
refROI.tform_allen2mouse = tform;
refROI.regIm_allen2mouse = regim;
refROI.sensoryMaps       = sensoryMaps;
refROI.allenROI          = allenROI;
refROI.refIm             = sensoryMaps.refIm;
refROI.bregma            = sensoryMaps.bregma;
refROI.midline           = sensoryMaps.midline;

%% apply registration to all ROIs
% ROIs in ref image space
areaCoord = {};
areaLbl   = {};
ROI       = allenROI.ROI;
preLbl    = allenROI.areaLbl;
[nX,nY]   = size(refROI.refIm);

if cfg.mergeMV2
  idx        = find(sum([strcmpi(preLbl,'VISam'); strcmpi(preLbl,'VISpm')]) > 0);
  mergeROI   = [];
  for iArea  = 1:numel(idx)
    mergeROI = [mergeROI; ROI{idx(iArea)}];
  end
  ROI(idx)      = [];
  preLbl(idx)   = [];
  ROI{end+1}    = mergeROI;
  preLbl{end+1} = 'mV2';
end

for iROI = 1:numel(cfg.ROIls)
  
%   if sum(strcmpi(preLbl{iROI},cfg.excludeROI)); continue; end
  
  ROIidx    = strcmpi(preLbl,cfg.ROIls{iROI});
  temp      = zeros(nX);
  for iPxl = 1:size(ROI{ROIidx},1)
    temp(ROI{ROIidx}(iPxl,1),ROI{ROIidx}(iPxl,2)) = 1;
  end
  Rfixed    = imref2d(size(temp));
  regim     = imwarp(temp,tform,'OutputView',Rfixed); 
  regim     = regim(1:nX,1:nY);
  [r,c]     = find(regim >= 1);
  ROIc      = cleanupROI([r c],nX);
  r         = ROIc(:,1);
  c         = ROIc(:,2);
  % apply transformation
%   [c,r]     = transformPointsForward(tform,ROI{ROIidx}(:,2),ROI{ROIidx}(:,1));
%   r         = round(r);
%   c         = round(c);
%   delidx    = r < 1 | r > nX | c < 1 | c > nY;
%   r(delidx) = [];
%   c(delidx) = [];
  
  % divide hemispheres, separate M2 if necessary
  [Ridx,Lidx]         = assignPxlToHem([r c],refROI.midline);
  if strcmpi(preLbl{ROIidx},'MOs') && cfg.divideM2
    if isnan(cfg.divideM2coord(1))
      coordL          = [0 refROI.bregma(1)-cfg.divideM2coord(2)*pxlPerMM];
      coordR          = [0 refROI.bregma(1)+cfg.divideM2coord(2)*pxlPerMM];
    elseif isnan(cfg.divideM2coord(2))
      coordL          = [refROI.bregma(2)+cfg.divideM2coord(1)*pxlPerMM bregma(1)];
      coordR          = [refROI.bregma(2)+cfg.divideM2coord(1)*pxlPerMM bregma(1)];
    else
      coordL          = [refROI.bregma(2)+cfg.divideM2coord(1)*pxlPerMM refROI.bregma(1)-cfg.divideM2coord(2)*pxlPerMM];
      coordR          = [refROI.bregma(2)+cfg.divideM2coord(1)*pxlPerMM  refROI.bregma(1)+cfg.divideM2coord(2)*pxlPerMM];
    end

    m2Idx             = find(r < max([nX coordL(1)]) & c < coordL(2));
    areaCoord{end+1}  = [r(intersect(Lidx,m2Idx)) c(intersect(Lidx,m2Idx))];
    areaLbl{end+1}    = ['a' preLbl{ROIidx} '-L'];
    
    m2Idx             = find(r < max([nX coordR(1)]) & c > coordR(2));
    areaCoord{end+1}  = [r(intersect(Ridx,m2Idx)) c(intersect(Ridx,m2Idx))];
    areaLbl{end+1}    = ['a' preLbl{ROIidx} '-R'];
    
    m2Idx             = find(r > min([0 coordL(1)]) & c > coordL(2));
    areaCoord{end+1}  = [r(intersect(Lidx,m2Idx)) c(intersect(Lidx,m2Idx))];
    areaLbl{end+1}    = ['m' preLbl{ROIidx} '-L'];
    
    m2Idx             = find(r > min([0 coordR(1)]) & c < coordR(2));
    areaCoord{end+1}  = [r(intersect(Ridx,m2Idx)) c(intersect(Ridx,m2Idx))];
    areaLbl{end+1}    = ['m' preLbl{ROIidx} '-R'];
    
  else
    areaCoord{end+1}  = [r(Lidx) c(Lidx)];
    areaLbl{end+1}    = [preLbl{ROIidx} '-L'];
    areaCoord{end+1}  = [r(Ridx) c(Ridx)];
    areaLbl{end+1}    = [preLbl{ROIidx} '-R'];
  end
end

refROI.areaCoord    = areaCoord;
refROI.areaLbl      = areaLbl;
refROI.areaCentroid = cellfun(@nanmean,areaCoord,'uniformOutput',false);

%% plot 
I       = zeros(size(refROI.refIm));
outline = I;
for iArea = 1:numel(areaCoord)
  subI = zeros(size(I));
  for iPxl = 1:size(areaCoord{iArea},1)
    I(areaCoord{iArea}(iPxl,1),areaCoord{iArea}(iPxl,2))    = iArea;
    subI(areaCoord{iArea}(iPxl,1),areaCoord{iArea}(iPxl,2)) = 1;
  end
  pm      = bwperim(subI);
  outline = outline + pm; 
end
outline(outline >= 1) = 1;

overlay = repmat(refROI.refIm./max(max(refROI.refIm)),[1 1 3]);
[ii,jj] = find(outline == 1);
for iP = 1:numel(ii)
  overlay(ii(iP),jj(iP),:) = [1 1 1];
end

refROI.ROIouline      = outline;
refROI.ROIoverlay     = overlay;

fh = figure; 
wf.applyFigDefaults(fh,[2 1],'k');
subplot(1,2,1); image(overlay); axis image; axis off
for iArea = 1:numel(areaCoord)
  centroid = nanmean(areaCoord{iArea});
  text(centroid(2),centroid(1),areaLbl{iArea},...
       'color','r','horizontalAlignment','center','fontsize',10)
end
title('Allen ROIs','color','w')

subplot(1,2,2); image(overlay); axis image; axis off
try; drawROI(sensoryMaps.ROIoutline,'y',gca); end
title('Vs. sensory maps','color','y')

%% save
save(cfg.fn,'refROI','cfg')
saveas(fh,cfg.figFn)
close(fh)

end

%% ---------------------------------------------------
function cfg = populateCfg(cfg)

if ~isfield(cfg,'transformType')
  cfg(1).transformType = 'similarity';
end
if ~isfield(cfg,'includeBregmaAndLambda')
  cfg(1).includeBregmaAndLambda = false;
end
if ~isfield(cfg,'ROIls')
  cfg(1).ROIls    = {'VISp','mV2','VISa','RSP','MOs','MOp','SS'}; 
end
if ~isfield(cfg,'V1type')
  cfg(1).V1type        = 'whole'; % tip, whole, center
end
if ~isfield(cfg,'tipPrctile')
  cfg(1).tipPrctile    = 2; % tip, whole, center
end
if ~isfield(cfg,'nPxlPoints')
  cfg(1).nPxlPoints    = 4;
end
if ~isfield(cfg,'lambdaDist')
  cfg(1).lambdaDist    = 3.1;
end
if ~isfield(cfg,'truncateRow')
  cfg(1).truncateRow   = [33 214];
end
if ~isfield(cfg,'truncateCol')
  cfg(1).truncateCol   = [11 218];
end
if ~isfield(cfg,'divideM2')
  cfg(1).divideM2      = true; 
end
if ~isfield(cfg,'excludeROI')
  cfg(1).excludeROI    = {'VISl','VISal','VISam','VISpm'}; 
end
if ~isfield(cfg,'mergeMV2')
  cfg(1).mergeMV2      = true; 
end
if ~isfield(cfg,'divideM2ap')
  cfg(1).divideM2coord = [nan 1.25]; % mm, [AP ML]
end
if ~isfield(cfg,'fn')
  cfg(1).fn = 'refROI.mat'; 
end
if ~isfield(cfg,'figFn')
  cfg(1).figFn = 'refROIAllen.fig'; 
end
end

%% ---------------------------------------------------
%% remove non-contiguous pixels and fill empty ones
function ROIc = cleanupROI(ROI,imsize,tol)

if nargin < 3; tol = 2; end

im = zeros(imsize);
for iPxl = 1:size(ROI,1)
  im(ROI(iPxl,1),ROI(iPxl,2)) = 1;
end

% remove orphan pixels (must have at least 2 in the neighborhood)
for iPxl = 1:size(ROI,1)
  subimrows = max([ROI(iPxl,1)-1 1]):min([ROI(iPxl,1)+1 imsize]);
  subimcols = max([ROI(iPxl,2)-1 1]):min([ROI(iPxl,2)+1 imsize]);
  subim     = im(subimrows,subimcols);
  if sum(subim(:)) < tol+1
    im(ROI(iPxl,1),ROI(iPxl,2)) = 0;
  end
end

% % fill in gaps, assuming single hemisphere
% im     = imfill(imgaussfilt(im));
% [r, c] = find(im ~= 0);
% im     = grayconnected(single(im),round(mean(r)),round(mean(c)),.9);
% 
% % reduce the border by 1 pxl (introduced by previous step)
% [r, ~] = find(im ~= 0);
% row    = round(mean(r));
% col    = find(im(row,:) > 0, 1, 'first');
% bound  = bwtraceboundary(im,[row col],'N');
% for iPxl = 1:size(bound,1)
%   im(bound(iPxl,1),bound(iPxl,2)) = 0;
% end

% return clean ROI
[r, c] = find(im ~= 0);
ROIc   = [r c];

end

%% ---------------------------------------------------
%% resize ROIs
function allenROI = resizeAllenROI(temp,cfg)

allenROI.pxlPerMM    = widefieldParams.pxlPerMM/widefieldParams.dsFactor;
allenROI.sizeRatio   = allenROI.pxlPerMM/temp.pxlPerMM;
allenROI.refFig      = temp.refFig(cfg.truncateRow(1):cfg.truncateRow(2),cfg.truncateCol(1):cfg.truncateCol(2));
allenROI.refFig      = round(imresize(allenROI.refFig,allenROI.sizeRatio));
allenROI.bregma      = round(allenROI.sizeRatio.*(temp.bregma - [cfg.truncateRow(1) cfg.truncateCol(1)]));
allenROI.areaLbl     = temp.areaLbl;

for iArea = 1:numel(temp.ROI)
  [r,c]               = find(allenROI.refFig == iArea);
  allenROI.ROI{iArea} = [r c];
end

allenROI.midline.a = -Inf;
allenROI.midline.b = Inf;
allenROI.midline.p = [allenROI.bregma; allenROI.bregma+[0 50]];

end