function airpuffsMap(rec,ISIflag,violetCorrectFlag)

if nargin < 1; rec = pwd;                 end
if nargin < 2; ISIflag = [];              end
if nargin < 3; violetCorrectFlag = false; end

tic
%% ------------------------------------------------------------------------
%% PREPROCESSING
%% look for stimulation file (deduct from rec name)
% mouse name (try 4 characters first)
if ~isempty(strfind(rec,'jukebox')) || ~isempty(strfind(rec,'mnt'))
  spockFlag = true;
else
  spockFlag = false;
end
[mouseID,rdate] = mouseAndDateFromFileName(rec);

% visual stim file path
if spockFlag
  serverpath   = widefieldParams.serverpathMap_spock;
else
  if ispc
    serverpath = widefieldParams.serverpathMap_pc;
  else
    serverpath = widefieldParams.serverpathMap_mac;
  end
end
fp  = formatFilePath(sprintf('%s%s/%s/',serverpath,mouseID,rdate));

% copy mapping file to local machine
cd(fp)
temp = dir('*whiskers*mat');
fn   = temp(1).name;
copyfile([fp fn],[formatFilePath(rec) fn]);

%% imaging frame rate
cd(rec)
if isempty(dir('info.mat'))
  try
    if ~isempty(dir('unprocessed'))
      [~,fl]   = generateFilelist([rec 'unprocessed/']); % imagej tends to mess up headers
    else
      [~,fl]   = generateFilelist(rec);
    end
    ts         = extractImageTimeStamps(fl{1});
    frameRate  = 1/mode(diff(ts));
    save info frameRate
  catch
    frameRate  = 20;
    save info frameRate
  end
else
  load info frameRate
end

%% motion correction / vascular detection
if isempty(dir('dff.mat'))
  fl         = generateFilelist(rec);
  try
    maskdata = removeVasculature_LP_visualMapping(fl);
    vascMask = maskdata.noData;
  catch
    vascMask = [];
  end
  close all
  mcorr    = getMotionCorrection(fl, false, false, 15, 5, false, 0.3, nan, 10);
  rawf = [];
  for iFile = 1:numel(fl)
    imstack = cv.imreadx(fl{iFile}, mcorr(iFile).xShifts, mcorr(iFile).yShifts); % load image shifts and apply them
    if isempty(vascMask); vascMask = false(size(imstack,1),size(imstack,2)); end
    rawf    = cat(3,rawf,extractF(double(imstack),vascMask,1));
%     clear imstack
  end
  
  if violetCorrectFlag
    frameRate  = frameRate /2;
    nZ         = size(rawf,3);
    strobeSeq  = estimateStrobeSeq(double(imstack(:,:,1:100)));
    strobeSeq  = repmat(strobeSeq,[1 floor(nZ/numel(strobeSeq))]);
    excFrames  = nZ - numel(strobeSeq);
    strobeSeq  = [strobeSeq strobeSeq(1:excFrames)];
    rawfViolet = rawf(:,:,strobeSeq==2);
    rawf       = rawf(:,:,strobeSeq==1);
    meanproj   = nanmean(rawf,3);
    
    save('dff','rawf','rawfViolet','vascMask','meanproj','-v7.3')
    save info frameRate strobeSeq
  else
    meanproj  = nanmean(rawf,3);
    save('dff','rawf','vascMask','meanproj','-v7.3')
  end
  
  
else
  fprintf('loading data...')
  load dff rawf meanproj
  fprintf('\n')
end

%% find and load file containing stim info
cd(rec)
fl  = dir('*mat');
fl  = {fl(:).name};

for iF = 1:numel(fl)
  fmatch = regexp(fl{iF},'[0-9]{8,}[_][0-9]{4,}','match');
  if ~isempty(fmatch); stimFile = fl{iF}; end
end

load(stimFile)
nPuffs = n;

%% ------------------------------------------------------------------------
%% ------------------------------------------------------------------------
%% ANALYSIS
%% parameters depend on whether it's Ca2+ or ISI
fprintf('computing map...')
if isempty(ISIflag)
  if strcmpi(mouseID(1:2),'vg')
    ISIflag = true;
  else
    ISIflag = false;
  end
end

if ISIflag
  win   = [1 5]; % sec, [pre post] airpuff
  th    = -1;   % z score threshold to consider response significant
  thsec = .5;    % amount
else
  win   = [.4 .4]; % sec, [pre post] airpuff
  th    = 1;   % z score threshold to consider response significant
  thsec = .1;    % amount
end

%% generate an average dff over time, where pre-airpuff period is F0
[nX,nY,~]  = size(rawf);

% in recs with violet correction, LED starts from frame 1, in others LED
% is triggered after a variable delay by the airpuff delivery code
if violetCorrectFlag
  f1       = 1;
else
  testF    = squeeze(nanmean(nanmean(rawf(round(nX/2):round(nX/2)+10,round(nY/2):round(nY/2)+10,:))));
  f1       = find(testF>prctile(testF,10),1,'first');
end

winFr      = round(win.*frameRate);
dffTrig    = zeros(nX,nY,winFr(2));
for iA = 1:nPuffs
  airpuffFrame = f1 + interval*frameRate;
  f1           = airpuffFrame;
  rf           = rawf(:,:,f1+1:f1+winFr(2)); % post-puff
  basel        = rawf(:,:,f1-winFr(1):f1); % pre-puff
  f0m          = repmat(nanmean(basel,3),[1 1 size(rf,3)]); % pre-puff mean
%   f0s          = repmat(nanstd(basel,0,3),[1 1 size(rf,3)]); % pre-
  dff          = (rf-f0m)./f0m;
  baseldff     = (basel - repmat(nanmean(basel,3),[1 1 size(basel,3)])) ./ repmat(nanmean(basel,3),[1 1 size(basel,3)]); % pre-puff
  dff0m        = repmat(nanmean(baseldff,3),[1 1 size(rf,3)]); % pre-puff mean
  dff0s        = repmat(nanstd(baseldff,0,3),[1 1 size(rf,3)]); % pre-puff std
  dffz         = (dff - dff0m) ./ dff0s;
  dffTrig      = dffTrig + dffz; % dff
end

dffTrig = dffTrig./nPuffs;
dffAvg  = nanmean(dffTrig,3);
th      = prctile(dffAvg(:),99)*sign(th);

%% apply pixel activity thresholds
% good pixels go beyond threshold for at thsec
[mat,bc,br] = conditionDffMat(dffTrig);
if th < 0
  thmat    = mat < th;
else
  thmat    = mat > th;
end

thmatTime   = cumsum(thmat);
valPxls     = sum((thmatTime > thsec*frameRate) & thmat) > 0;
map         = conditionDffMat(valPxls,bc,br,[nX nY 1]);

%% smooth and threshold again to obtain more continuous areas
smthSigma           = 10;
maptemp             = map;
maptemp(isnan(map)) = 0;
mapSmth             = imgaussfilt(maptemp,smthSigma);
mapSmthTh           = mapSmth;
if th > 0
  mapSmthTh(mapSmth < th/smthSigma) = 0;
else
  mapSmthTh(mapSmth > th/smthSigma) = 0;
end

%% segmentation
% (full disclosure this was more or less directly lifted from mathworks example)
I              = double(mapSmthTh);
I(isnan(I))    = 0;
I(I>0)         = 1;
[r,c]          = find(I>0);
coord          = [r c];

% edge detection
[~, threshold] = edge(I, 'sobel');
fudgeFactor    =  1;
outline        = imclearborder(edge(I,'sobel', threshold * fudgeFactor),4);

% % dilate borders to fill in gaps
% se90           = strel('line', 10, 90);
% se0            = strel('line', 10, 0);
% edgeDil        = imdilate(edges, [se90 se0]);
% % fill in gaps and smooth borders
% edgeFill = imfill(edgeDil, 'holes');
% seD = strel('diamond',2);
% edgeFill = imerode(edgeFill,seD);
% % detect outlines of putative ROIS
% coord   = bwboundaries(edgeFill);
% outline = imclearborder(edge(edgeFill,'sobel', threshold * fudgeFactor),4);

%% plot and save
overlay = repmat(meanproj./max(max(meanproj)),[1 1 3]);
[ii,jj] = find(outline == 1);
for iP = 1:numel(ii)
  overlay(ii(iP),jj(iP),:) = [1 0 0];
end

% [ii,jj] = find(I > 0);
% for iP = 1:numel(ii)
%   overlay(ii(iP),jj(iP),:) = [0 1 1];
% end

fprintf('\ndone after %1.1f min\n',toc/60)

fh = figure; 
subplot(1,2,1); imagesc(mapSmthTh); colorbar; axis image; axis off; title('Airpuff resp. z-score (smoothed)')
subplot(1,2,2); image(overlay); axis image; axis off; title('Area borders')
set(fh, 'position', [100 200 600 250])

%%
saveas(fh,'mapOnBrain.fig')
saveas(fh,'mapOnBrain.png','png')
close(fh)
save(['map_' fn],'coord','outline','map','dffTrig','mapSmth','mapSmthTh', ...
                 'meanproj','th','thsec','win','smthSigma','stimFile')

%% also for Violet, generate an average dff over time, where pre-airpuff period is F0
% here it will be just for comparison purposes, so i will just extract an
% average rawf for the region extracted above and compare triggered
% response with and without correction
if violetCorrectFlag
  if ~exist('rawfViolet','var'); load dff rawfViolet; end
  if ~exist('rawf','var');       load dff rawf;       end
  load info frameRate strobeSeq
  
  nZb        = size(rawf,3);
  nZv        = size(rawfViolet,3);
  mask       = mapSmthTh > 0;
  rawfb      = zeros(nZb,1);
  rawfv      = zeros(nZv,1);
  for iFrame = 1:nZb
    thisf          = rawf(:,:,iFrame);
    rawfb(iFrame)  = nanmean(thisf(mask));
  end
  for iFrame = 1:nZv
    thisf          = rawfViolet(:,:,iFrame);
    rawfv(iFrame)  = nanmean(thisf(mask));
  end
  rawfv      = interp1(find(strobeSeq==2)',rawfv,find(strobeSeq==1)','linear','extrap');
  
  win              = [1 5];
  winFr            = round(win.*frameRate);
  dffTrigBlue      = nan(sum(winFr),nPuffs);
  dffTrigViolet    = nan(sum(winFr),nPuffs);
  dffTrigVioletC   = nan(sum(winFr),nPuffs);
  dffTrigCorrected = nan(sum(winFr),nPuffs);
  f1               = 1;
  
  cf  = fit(rawfv,rawfb,'poly1');
  rawfvc = fiteval(cf,rawfv);
  
  for iA = 1:nPuffs
    airpuffFrame = f1 + round(interval*frameRate);
    f1           = airpuffFrame;
    
    % blue
    rf                  = rawfb(f1-winFr(1)+1:f1+winFr(2)); % post-puff
    basel               = rawfb(f1-winFr(1)+1:f1); % pre-puff
    f0m                 = repmat(nanmean(basel),[numel(rf) 1]); % pre-puff mean
    dff                 = (rf-f0m)./f0m;
    dffTrigBlue(:,iA)   = dff; % dff
    
    % violet
    rf                  = rawfv(f1-winFr(1)+1:f1+winFr(2)); % post-puff
    basel               = rawfv(f1-winFr(1)+1:f1); % pre-puff
    f0m                 = repmat(nanmean(basel),[numel(rf) 1]); % pre-puff mean
    dff                 = (rf-f0m)./f0m;
%     dff                 = approxGaussianBlur1D(dff,3,1);
    dffTrigViolet(:,iA) = dff; % dff
    
    % violet
    rf                  = rawfvc(f1-winFr(1)+1:f1+winFr(2)); % post-puff
    basel               = rawfvc(f1-winFr(1)+1:f1); % pre-puff
    f0m                 = repmat(nanmean(basel),[numel(rf) 1]); % pre-puff mean
    dff                 = (rf-f0m)./f0m;
%     dff                 = approxGaussianBlur1D(dff,3,1);
    dffTrigVioletC(:,iA) = dff; % dff
  end
  
  % corrected

  for iA = 1:nPuffs
    dffTrigCorrected(:,iA) = (dffTrigBlue(:,iA)+1)/(dffTrigVioletC(:,iA)+1) -1 ;
    dffTrigCorrected(:,iA) = dffTrigCorrected(:,iA) - mean(dffTrigCorrected(1:round(win*frameRate),iA));
  end
  
  dffMean  = mean(dffTrigCorrected,2);
  dffSem   = std(dffTrigCorrected,0,2)./sqrt(nPuffs-1);
  dffMeanB = mean(dffTrigBlue,2);
  dffSemB  = std(dffTrigBlue,0,2)./sqrt(nPuffs-1);
  dffMeanV = mean(dffTrigViolet,2);
  dffSemV  = std(dffTrigViolet,0,2)./sqrt(nPuffs-1);
  taxis    = linspace(-win(1)+1/frameRate,win(2),sum(winFr));
  
  save violetCorrection dffTrigBlue dffTrigViolet dffTrigCorrected dffMean dffSem ...
       dffMeanB dffSemB dffMeanV dffSemV taxis
  
  wf = widefieldParams;
  figure; wf.applyFigDefaults(gcf,[1 1],'w'); hold on
  plot(taxis, dffMeanV, '-', 'color', widefieldParams.mypurple, 'linewidth', 1)
  plot(taxis, dffMeanB, '-', 'color', widefieldParams.myblue,   'linewidth', 1)
  plot(taxis, dffMean,  '-', 'color', widefieldParams.darkgray, 'linewidth', 1)
  legend({'405','473','corrected'},'location','best'); legend('boxoff')
  plot(taxis, dffMeanV - dffSemV, '--', 'color', widefieldParams.mypurple, 'linewidth', .5)
  plot(taxis, dffMeanB - dffSemV, '--', 'color', widefieldParams.myblue,   'linewidth', .5)
  plot(taxis, dffMean - dffSem,   '--', 'color', widefieldParams.darkgray, 'linewidth', .5)
  plot(taxis, dffMeanV + dffSemV, '--', 'color', widefieldParams.mypurple, 'linewidth', .5)
  plot(taxis, dffMeanB + dffSemV, '--', 'color', widefieldParams.myblue,   'linewidth', .5)
  plot(taxis, dffMean + dffSem,   '--', 'color', widefieldParams.darkgray, 'linewidth', .5)
  axis tight; xl = get(gca,'xlim'); yl = get(gca,'ylim');
  plot(xl,[0 0],'--','color',widefieldParams.lightgray);
  plot([0 0],yl,'--','color',widefieldParams.lightgray);
  wf.applyAxisDefaults(gca,'k'); 
  wf.applyAxisLbls(gca,'Time from airpuff (s)','\DeltaF/F','Airpuff-triggered response')
  saveas(gcf,'violetCorrection')
  close
  
end