function avgDffSumm = summarizeAvgDffvsEvidence(mice)

% avgDffSumm = summarizeAvgDff(mice)
% for accumulation maze and trial types compiles sensory evidence tuning 
% (space avg, COM, pxl sequences, ANOVA)
% mice is cell array with mouse names

%% ------------------------------------------------------------------------

if nargin < 1; mice  = {'ai2';'ai3';'ai5';'ai7';'ai9';'ai10'}; end

try
%% flag true to load from local disk instead of bucket
% (will have no effect if running on spock)
localFlag    = false; 

%% configuration for ROI corr analysis
cfg.trialTypes         = {'correct','error'};
cfg.bilateralROIflag   = true;
cfg.mice               = mice;
cfg.evBins             = -15:5:15;
wf                     = widefieldParams;

rootdir = wf.getRootDir(isThisSpock,localFlag);
cd(rootdir)

%% ------------------------------------------------------------------------
%% compile avg t4, t11 for each mouse (session)
%% ------------------------------------------------------------------------

fprintf('retrieving data...\n')
avgDffSumm.ROIlbl  = {};

% loop through recs
for iMouse  = 1:numel(mice)
  
  %% task (visual guide and accumulation)
  % average across recs for each animal first
  taskRecs                       = wf.getMouseRecs(mice{iMouse},'task',localFlag);
  avgDffSumm.mouse(iMouse).recls = taskRecs;
  for iRec = 1:numel(taskRecs)
    fprintf('\tmouse %d/%d, rec %d/%d\n',iMouse,numel(mice),iRec,numel(taskRecs))
    cd(taskRecs{iRec})
    
    %% collect ROI averages
%     load avgROIvsEvidence avgROI
%     avgDff = avgROI;
    avgDff = avgDffROIvsEvidence(pwd,true,cfg.evBins,false,false);
    
    avgDffSumm.ROIlbl      = avgDff.ROIlbl;
    avgDffSumm.spaceBins   = avgDff.spaceBins;
    avgDffSumm.RminusLvals = avgDff.RminusLvals;
    avgDffSumm.RminusLBins = avgDff.RminusLBins;
    
    for iType = 1:numel(cfg.trialTypes)
      mat    = avgDff.maze(end).(cfg.trialTypes{iType}).avg;
      avgDffSumm.mouse(iMouse).ROI.accumul.(cfg.trialTypes{iType}).recs(:,:,:,iRec)              ...
             = mat;
      mat    = avgDff.maze(1).(cfg.trialTypes{iType}).avg;
      avgDffSumm.mouse(iMouse).ROI.visGuide.(cfg.trialTypes{iType}).recs(:,:,:,iRec)             ...
             = mat;
    end

    clear avgDff avgROI
    
  end
  
  %% calculate averages, sem and corr across recs for this mouse
  for iType = 1:numel(cfg.trialTypes)
    mat    = avgDffSumm.mouse(iMouse).ROI.accumul.(cfg.trialTypes{iType}).recs;
    avgDffSumm.(['ROIavg_accumul_' (cfg.trialTypes{iType})])(:,:,:,iMouse)        = nanmean(mat,4);
    avgDffSumm.(['ROIavg_accumul_' (cfg.trialTypes{iType}) '_sem'])(:,:,:,iMouse) = nanstd(mat,0,4)./sqrt(size(mat,4)-1);
    
    mat    = avgDffSumm.mouse(iMouse).ROI.visGuide.(cfg.trialTypes{iType}).recs;
    avgDffSumm.(['ROIavg_visGuide_' (cfg.trialTypes{iType})])(:,:,:,iMouse)        = nanmean(mat,4);
    avgDffSumm.(['ROIavg_visGuide_' (cfg.trialTypes{iType}) '_sem'])(:,:,:,iMouse) = nanstd(mat,0,4)./sqrt(size(mat,4)-1);
  end
  
end

cd(rootdir)

%% ------------------------------------------------------------------------
%% stats
%% ------------------------------------------------------------------------


%% ------------------------------------------------------------------------
%% save and plot
%% ------------------------------------------------------------------------
cd(rootdir)
avgDffSumm.cfg = cfg;

save avgDffSummaryVsEvidence avgDffSumm cfg -v7.3
if ~isempty(dir('avgDffSummaryVsEvidence.pdf')); delete avgDffSummaryVsEvidence.pdf; end

%% plot ROI avgs: correct, error for both mazes, as a function of evidence
nROI    = numel(avgDffSumm.ROIlbl);
nEv     = numel(avgDffSumm.RminusLvals);
wf      = widefieldParams;
cl      = red2blue(nEv);
[nr,nc] = subplotOrg(nROI,4);
  
for iType = 1:numel(cfg.trialTypes)
  avgAT   = nanmean(avgDffSumm.(['ROIavg_accumul_' (cfg.trialTypes{iType})]),4);
  avgVG   = nanmean(avgDffSumm.(['ROIavg_visGuide_' (cfg.trialTypes{iType})]),4);

  fh(numel(cfg.trialTypes)*(iType-1)+1)   = figure;
  fh(numel(cfg.trialTypes)*(iType-1)+2)   = figure;
  wf.applyFigDefaults(fh(1),[nc nr],'w')
  wf.applyFigDefaults(fh(2),[nc nr],'w')
  for iROI = 1:nROI
    figure(fh(numel(cfg.trialTypes)*(iType-1)+1))
    subplot(nr,nc,iROI); hold on
    for iEv = 1:nEv
      plot(avgDffSumm.spaceBins,avgAT(:,iROI,iEv),'-','color',cl(iEv,:),'linewidth',1.5)
    end
    wf.applyAxisDefaults(gca,'k')
    if iROI == 1
      wf.applyAxisLbls(gca,'y pos (cm)','\DeltaF/F',sprintf('%s - Towers, %s',avgDffSumm.ROIlbl{iROI},cfg.trialTypes{iType}));
    else
      wf.applyAxisLbls(gca,[],[],avgDffSumm.ROIlbl{iROI});
    end
    if iROI == nROI
      legend(cellfun(@num2str,num2cell(avgDffSumm.RminusLvals),'UniformOutput',false),'location','best')
    end

    figure(fh(numel(cfg.trialTypes)*(iType-1)+2))
    subplot(nr,nc,iROI); hold on
    for iEv = 1:nEv
      plot(avgDffSumm.spaceBins,avgVG(:,iROI,iEv),'-','color',cl(iEv,:),'linewidth',1.5)
    end
    wf.applyAxisDefaults(gca,'k')
    if iROI == 1
      wf.applyAxisLbls(gca,'y pos (cm)','\DeltaF/F',sprintf('%s - VisGuide, %s',avgDffSumm.ROIlbl{iROI},cfg.trialTypes{iType}));
    else
      wf.applyAxisLbls(gca,[],[],avgDffSumm.ROIlbl{iROI});
    end
    if iROI == nROI
      legend(cellfun(@num2str,num2cell(avgDffSumm.RminusLvals),'UniformOutput',false),'location','best')
    end
  end
end

%% export
figls = get(0,'children'); % print to pdf
for ii = length(figls):-1:1
  figure(figls(ii))
  export_fig avgDffSummaryVsEvidence.pdf -q101 -append
end
close all

catch ME
  displayException(ME)
end

end