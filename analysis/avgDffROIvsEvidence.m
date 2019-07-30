function avgROI = avgDffROIvsEvidence(rec,doZScore,evBins,loadFlag,plotFlag)

% avgDff = avgDffROIvsEvidence(rec,doZScore,evBins)
% rec is path, doZscore logical flag, evBins are bins of \Delta towers 
% calculates average +/- sem dff for different amounts of final evidence
% LP dec 2018

%%
if nargin < 1; rec      = formatFilePath(pwd); end
if nargin < 2; doZScore = true;                end
if nargin < 3; evBins   = -15:5:15;            end
if nargin < 4; loadFlag = true;                end
if nargin < 5; plotFlag = true;                end
%%
tic
fprintf('calculating average dffs')

%% load if available and return
cd(rec)
if ~isempty(dir('avgROIvsEvidence.mat')) && loadFlag
  warning('off','all')
  load avgROIvsEvidence avgROI
  warning('on','all')
  if exist('avgROI','var')
    fprintf(' done after %1.1f min (loaded from disk)\n',toc/60)
    return
  end
else
  load dffROI dffROI ROIlbl
  dff = dffROI; clear dffROI
  load behavLog logSumm
end

%% some recs in blocks condensed have maze 12, this is just like maze 4 but with towers on both sides
% here they will be analyzed together for convenience
if contains(char(logSumm.info.protocol),'Condensed')
  logSumm.currMaze(logSumm.currMaze == 12) = 4;
end

%% analysis parameters
avgROI.ROIlbl             = ROIlbl;
avgROI.trialTypeList      = {'correct','error'};
avgROI.mazeList           = unique(logSumm.currMaze);
avgROI.RminusLBins        = evBins;
avgROI.RminusLvals        = toBinCenters(evBins);
avgROI.spaceBins          = 0:5:300;
avgROI.isZscored          = doZScore;
nSpace                    = numel(avgROI.spaceBins);
nEv                       = numel(avgROI.RminusLvals);
nROI                      = numel(ROIlbl);

%% z-score
if doZScore; dff = zscore(dff); end

%% loop over trial types, mazes and evidence bins to compile average dff 
for iType = 1:numel(avgROI.trialTypeList)
  fprintf('.')
  for iMaze = 1:numel(avgROI.mazeList)
    
    %% get list of relevant trials and trial-aligned dff
    [dffTrials,~,~,trialidx] = alignTrials(dff,logSumm,avgROI.trialTypeList{iType}, ...
                                           avgROI.mazeList(iMaze),'cueStart',avgROI.spaceBins,'space');    
    avgROI.maze(iMaze).(avgROI.trialTypeList{iType}).avg = nan(nSpace,nROI,nEv);
    avgROI.maze(iMaze).(avgROI.trialTypeList{iType}).sem = nan(nSpace,nROI,nEv);
    if isempty(trialidx); continue; end
    
    %% loop over evidence values
    RmL = logSumm.nCues_RminusL(trialidx);
    for iEv = 1:nEv
      thisidx = find(RmL >= evBins(iEv) & RmL < evBins(iEv+1));
      if isempty(thisidx); continue; end
      thismat = squeeze(nanmean(dffTrials(thisidx,:,:),1));
      avgROI.maze(iMaze).(avgROI.trialTypeList{iType}).avg(:,:,iEv) = thismat;
      thismat = squeeze(nanstd(dffTrials(thisidx,:,:),0,1))./sqrt(numel(thisidx)-1);
      avgROI.maze(iMaze).(avgROI.trialTypeList{iType}).sem(:,:,iEv) = thismat;
    end
  end
end

%% plot
if plotFlag
avgAT   = avgROI.maze(end).correct.avg;
avgVG   = avgROI.maze(1).correct.avg;
wf      = widefieldParams;
cl      = red2blue(nEv);
[nr,nc] = subplotOrg(nROI,4);

fh(1)   = figure;
fh(2)   = figure;
wf.applyFigDefaults(fh(1),[nc nr],'w')
wf.applyFigDefaults(fh(2),[nc nr],'w')
for iROI = 1:nROI
  figure(fh(1))
  subplot(nr,nc,iROI); hold on
  for iEv = 1:nEv
    plot(avgROI.spaceBins,avgAT(:,iROI,iEv),'-','color',cl(iEv,:),'linewidth',1.5)
  end
  wf.applyAxisDefaults(gca,'k')
  if iROI == 1
    wf.applyAxisLbls(gca,'y pos (cm)','\DeltaF/F',sprintf('%s - Towers',ROIlbl{iROI}));
  else
    wf.applyAxisLbls(gca,[],[],ROIlbl{iROI});
  end
  if iROI == nROI
    legend(cellfun(@num2str,num2cell(avgROI.RminusLvals),'UniformOutput',false),'location','best')
  end
  
  figure(fh(2))
  subplot(nr,nc,iROI); hold on
  for iEv = 1:nEv
    plot(avgROI.spaceBins,avgVG(:,iROI,iEv),'-','color',cl(iEv,:),'linewidth',1.5)
  end
  wf.applyAxisDefaults(gca,'k')
  if iROI == 1
    wf.applyAxisLbls(gca,'y pos (cm)','\DeltaF/F',sprintf('%s - VisGuide',ROIlbl{iROI}));
  else
    wf.applyAxisLbls(gca,[],[],ROIlbl{iROI});
  end
  if iROI == nROI
    legend(cellfun(@num2str,num2cell(avgROI.RminusLvals),'UniformOutput',false),'location','best')
  end
end

if ~isempty(dir('avgROIvsEvidence.pdf')); delete avgROIvsEvidence.pdf; end
figure(fh(1))
export_fig('avgROIvsEvidence.pdf','-append')
close
figure(fh(2))
export_fig('avgROIvsEvidence.pdf','-append')
close
end
%% save
save avgROIvsEvidence avgROI -v7.3
fprintf('\tdone after %1.1f min\n',toc/60)