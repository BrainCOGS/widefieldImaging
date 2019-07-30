function avgDffSumm = summarizeAvgDffTrig(mice)

% avgDffSumm = summarizeAvgDffTrig(mice)
% for mazes and trial types compiles avg ROI tower- and tur-triggered resps
% mice is cell array with mouse names

%% ------------------------------------------------------------------------

if nargin < 1; mice  = {'ai2';'ai3';'ai5';'ai7';'ai9';'ai10'}; end

try
%% flag true to load from local disk instead of bucket
% (will have no effect if running on spock)
localFlag    = false; 

%% configuration for analysis
cfg.trialTypes         = {'correct','error'};
cfg.bilateralROIflag   = true;
cfg.mice               = mice;
cfg.summWin            = [0 2];
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
    load avgDffEventTrig_ROI avgDff
%     avgDff = avgDff_eventTriggered;
    avgDff.timeAxis            = avgDff.timeAxis + mode(diff(avgDff.timeAxis));
    avgDff.timeAxis_preturn    = avgDff.timeAxis_preturb + mode(diff(avgDff.timeAxis));
    avgDffSumm.ROIlbl          = avgDff.ROIlbl;
    avgDffSumm.isZscoredAll    = avgDff.isZscoredAll;
    avgDffSumm.isBaseSubt      = avgDff.isBaseSubt;
    avgDffSumm.isZscoredBasel  = avgDff.isZscoredBasel;
    if round(avgDff.fps) == 10
      avgDffSumm.timeAxis          = avgDff.timeAxis;  
      avgDffSumm.timeAxis_preturn  = avgDff.timeAxis_preturn;  
      downsample                   = false;
    else
      downsample                   = true;
    end
    
    for iType = 1:numel(cfg.trialTypes)
      % towers task
      turnR     = avgDff.maze(end).(cfg.trialTypes{iType}).turnR.avg;
      turnL     = avgDff.maze(end).(cfg.trialTypes{iType}).turnL.avg;
      preturnR  = avgDff.maze(end).(cfg.trialTypes{iType}).preturnR.avg;
      preturnL  = avgDff.maze(end).(cfg.trialTypes{iType}).preturnL.avg;
      towerR    = avgDff.maze(end).(cfg.trialTypes{iType}).towerR.avg;
      towerL    = avgDff.maze(end).(cfg.trialTypes{iType}).towerL.avg;
      
      if downsample
        turnR     = interp1(avgDff.timeAxis',turnR,avgDffSumm.timeAxis');
        turnL     = interp1(avgDff.timeAxis',turnL,avgDffSumm.timeAxis');
        preturnR  = interp1(avgDff.timeAxis',preturnR,avgDffSumm.timeAxis');
        preturnL  = interp1(avgDff.timeAxis',preturnL,avgDffSumm.timeAxis');
        towerR    = interp1(avgDff.timeAxis',towerR,avgDffSumm.timeAxis');
        towerL    = interp1(avgDff.timeAxis',towerL,avgDffSumm.timeAxis');
      end
      
      avgDffSumm.mouse(iMouse).accumul.(cfg.trialTypes{iType}).turnR.recs(:,:,iRec)     = turnR;
      avgDffSumm.mouse(iMouse).accumul.(cfg.trialTypes{iType}).turnL.recs(:,:,iRec)     = turnL;
      avgDffSumm.mouse(iMouse).accumul.(cfg.trialTypes{iType}).preturnR.recs(:,:,iRec)  = preturnR;
      avgDffSumm.mouse(iMouse).accumul.(cfg.trialTypes{iType}).preturnL.recs(:,:,iRec)  = preturnL;
      avgDffSumm.mouse(iMouse).accumul.(cfg.trialTypes{iType}).towerR.recs(:,:,iRec)    = towerR;
      avgDffSumm.mouse(iMouse).accumul.(cfg.trialTypes{iType}).towerL.recs(:,:,iRec)    = towerL;
           
      % visually-guided task
      turnR     = avgDff.maze(1).(cfg.trialTypes{iType}).turnR.avg;
      turnL     = avgDff.maze(1).(cfg.trialTypes{iType}).turnL.avg;
      preturnR  = avgDff.maze(1).(cfg.trialTypes{iType}).preturnR.avg;
      preturnL  = avgDff.maze(1).(cfg.trialTypes{iType}).preturnL.avg;
      towerR    = avgDff.maze(1).(cfg.trialTypes{iType}).towerR.avg;
      towerL    = avgDff.maze(1).(cfg.trialTypes{iType}).towerL.avg;
      
      if downsample
        turnR     = interp1(avgDff.timeAxis',turnR,avgDffSumm.timeAxis');
        turnL     = interp1(avgDff.timeAxis',turnL,avgDffSumm.timeAxis');
        preturnR  = interp1(avgDff.timeAxis',preturnR,avgDffSumm.timeAxis');
        preturnL  = interp1(avgDff.timeAxis',preturnL,avgDffSumm.timeAxis');
        towerR    = interp1(avgDff.timeAxis',towerR,avgDffSumm.timeAxis');
        towerL    = interp1(avgDff.timeAxis',towerL,avgDffSumm.timeAxis');
      end
      
      avgDffSumm.mouse(iMouse).visGuide.(cfg.trialTypes{iType}).turnR.recs(:,:,iRec)     = turnR;
      avgDffSumm.mouse(iMouse).visGuide.(cfg.trialTypes{iType}).turnL.recs(:,:,iRec)     = turnL;
      avgDffSumm.mouse(iMouse).visGuide.(cfg.trialTypes{iType}).preturnR.recs(:,:,iRec)  = preturnR;
      avgDffSumm.mouse(iMouse).visGuide.(cfg.trialTypes{iType}).preturnL.recs(:,:,iRec)  = preturnL;
      avgDffSumm.mouse(iMouse).visGuide.(cfg.trialTypes{iType}).towerR.recs(:,:,iRec)    = towerR;
      avgDffSumm.mouse(iMouse).visGuide.(cfg.trialTypes{iType}).towerL.recs(:,:,iRec)    = towerL;
    end

    clear avgDff 
  end
  
  %% calculate averages, sem across recs for this mouse
  vecs = {'turnR','turnL','preturnR','preturnL','towerR','towerL'};
  for iType = 1:numel(cfg.trialTypes)
    for iVec = 1:numel(vecs)
      mat    = avgDffSumm.mouse(iMouse).accumul.(cfg.trialTypes{iType}).(vecs{iVec}).recs;
      avgDffSumm.(['ROIavg_accumul_' (cfg.trialTypes{iType})]).(vecs{iVec})(:,:,iMouse)        = nanmean(mat,3);
      avgDffSumm.(['ROIavg_accumul_' (cfg.trialTypes{iType}) '_sem']).(vecs{iVec})(:,:,iMouse) = nanstd(mat,0,3)./sqrt(size(mat,3)-1);

      mat    = avgDffSumm.mouse(iMouse).visGuide.(cfg.trialTypes{iType}).(vecs{iVec}).recs;
      avgDffSumm.(['ROIavg_visGuide_' (cfg.trialTypes{iType})]).(vecs{iVec})(:,:,iMouse)        = nanmean(mat,3);
      avgDffSumm.(['ROIavg_visGuide_' (cfg.trialTypes{iType}) '_sem']).(vecs{iVec})(:,:,iMouse) = nanstd(mat,0,3)./sqrt(size(mat,3)-1);
    end
  end
  
end

cd(rootdir)

%% ------------------------------------------------------------------------
%% stats
%% ------------------------------------------------------------------------
events    = {'turn','preturn','tower'};
nROI      = numel(avgDffSumm.ROIlbl);
lblsAll   = unique(cellfun(@(x)(x(1:end-2)),avgDffSumm.ROIlbl,'UniformOutput',false),'stable');
[~ ,lbls] = getDefaultROIcl(lblsAll);
avgDffSumm.stats.ROIlbls_unil = lbls;
avgDffSumm.stats.testname     = 'ttest';
for iEv = 1:numel(events)
  % time averages
  if strcmp(events{iEv},'preturn')
    sumIdx  = find(avgDffSumm.timeAxis_preturn >= -cfg.summWin(2),1,'first'):find(avgDffSumm.timeAxis_preturn <= cfg.summWin(1),1,'last'); 
  else
    sumIdx  = find(avgDffSumm.timeAxis >= cfg.summWin(1),1,'first'):find(avgDffSumm.timeAxis <= cfg.summWin(2),1,'last'); 
  end
  avgATc_R      = avgDffSumm.ROIavg_accumul_correct.([events{iEv} 'R']);
  avgATc_R      = squeeze(nanmean(avgATc_R(sumIdx,:,:)))';
  avgATc_L      = avgDffSumm.ROIavg_accumul_correct.([events{iEv} 'L']);
  avgATc_L      = squeeze(nanmean(avgATc_L(sumIdx,:,:)))';
  avgATc_contra = (avgATc_R(:,1:2:end) + avgATc_L(:,2:2:end))./2; 
  avgATc_ipsi   = (avgATc_L(:,1:2:end) + avgATc_R(:,2:2:end,1))./2;
  
  avgATe_R      = avgDffSumm.ROIavg_accumul_error.([events{iEv} 'R']);
  avgATe_R      = squeeze(nanmean(avgATe_R(sumIdx,:,:)))';
  avgATe_L      = avgDffSumm.ROIavg_accumul_error.([events{iEv} 'L']);
  avgATe_L      = squeeze(nanmean(avgATe_L(sumIdx,:,:)))';
  avgATe_contra = (avgATe_R(:,1:2:nROI) + avgATe_L(:,2:2:end))./2;
  avgATe_ipsi   = (avgATe_L(:,1:2:nROI) + avgATe_R(:,2:2:end,1))./2;
  
  avgVGc_R      = avgDffSumm.ROIavg_visGuide_correct.([events{iEv} 'R']);
  avgVGc_R      = squeeze(nanmean(avgVGc_R(sumIdx,:,:)))';
  avgVGc_L      = avgDffSumm.ROIavg_visGuide_correct.([events{iEv} 'L']);
  avgVGc_L      = squeeze(nanmean(avgVGc_L(sumIdx,:,:)))';
  avgVGc_contra = (avgVGc_R(:,1:2:nROI) + avgVGc_L(:,2:2:end))./2;
  avgVGc_ipsi   = (avgVGc_L(:,1:2:nROI) + avgVGc_R(:,2:2:end,1))./2;
  
  avgVGe_R      = avgDffSumm.ROIavg_visGuide_error.([events{iEv} 'R']);
  avgVGe_R      = squeeze(nanmean(avgVGe_R(sumIdx,:,:)))';
  avgVGe_L      = avgDffSumm.ROIavg_visGuide_error.([events{iEv} 'L']);
  avgVGe_L      = squeeze(nanmean(avgVGe_L(sumIdx,:,:)))';
  avgVGe_contra = (avgVGe_R(:,1:2:nROI) + avgVGe_L(:,2:2:end))./2;
  avgVGe_ipsi   = (avgVGe_L(:,1:2:nROI) + avgVGe_R(:,2:2:end,1))./2;
  
  % compile time averages
  avgDffSumm.timeAvg_ROIavg_accumul_correct_contra.(events{iEv})  = avgATc_contra;
  avgDffSumm.timeAvg_ROIavg_accumul_correct_ipsi.(events{iEv})    = avgATc_ipsi;
  avgDffSumm.timeAvg_ROIavg_accumul_error_contra.(events{iEv})    = avgATe_contra;
  avgDffSumm.timeAvg_ROIavg_accumul_error_ipsi.(events{iEv})      = avgATe_ipsi;
  avgDffSumm.timeAvg_ROIavg_visGuide_correct_contra.(events{iEv}) = avgVGc_contra;
  avgDffSumm.timeAvg_ROIavg_visGuide_correct_ipsi.(events{iEv})   = avgVGc_ipsi;
  avgDffSumm.timeAvg_ROIavg_visGuide_error_contra.(events{iEv})   = avgVGe_contra;
  avgDffSumm.timeAvg_ROIavg_visGuide_error_ipsi.(events{iEv})     = avgVGe_ipsi;
  
  % 2-way rm ANOVA (ROI, mice and comparison -- either correct vs error, tasks, ipsi vs contra)
  [nmice,nROIu] = size(avgATc_contra);
  ROIvec        = repmat(1:nROIu,[nmice 1]);
  ROIvec        = repmat(ROIvec(:),[2 1]);
  mousevec      = repmat((1:nmice)',[1 nROIu]);
  mousevec      = repmat(mousevec(:),[2 1]);
  conditvec     = [ones(nmice*nROIu,1); 2*ones(nmice*nROIu,1)];
  
  [avgDffSumm.stats.(events{iEv}).p_ANOVA_accumul_correct_contraVSipsi,avgDffSumm.stats.(events{iEv}).table_ANOVA_accumul_correct_contraVSipsi,avgDffSumm.stats.(events{iEv}).stats_ANOVA_accumul_correct_contraVSipsi] ...
               = anovan([avgATc_contra(:); avgATc_ipsi(:)],{conditvec,ROIvec,mousevec},'varnames',{'side','ROI','mice'},'model','interaction','display','off');
  avgDffSumm.stats.(events{iEv}).multComp_ANOVA_accumul_correct_contraVSipsi      ...
               = multcompare(avgDffSumm.stats.(events{iEv}).stats_ANOVA_accumul_correct_contraVSipsi,'display','off');
  [avgDffSumm.stats.(events{iEv}).p_ANOVA_visGuide_correct_contraVSipsi,avgDffSumm.stats.(events{iEv}).table_ANOVA_visGuide_correct_contraVSipsi,avgDffSumm.stats.(events{iEv}).stats_ANOVA_visGuide_correct_contraVSipsi] ...
               = anovan([avgVGc_contra(:); avgVGc_ipsi(:)],{conditvec,ROIvec,mousevec},'varnames',{'side','ROI','mice'},'model','interaction','display','off');
  avgDffSumm.stats.(events{iEv}).multComp_ANOVA_visGuide_correct_contraVSipsi     ...
               = multcompare(avgDffSumm.stats.(events{iEv}).stats_ANOVA_visGuide_correct_contraVSipsi,'display','off');
  [avgDffSumm.stats.(events{iEv}).p_ANOVA_accumul_correctVSerror_contra,avgDffSumm.stats.(events{iEv}).table_ANOVA_accumul_correctVSerror_contra,avgDffSumm.stats.(events{iEv}).stats_ANOVA_accumul_correctVSerror_contra] ...
               = anovan([avgATc_contra(:); avgATe_contra(:)],{conditvec,ROIvec,mousevec},'varnames',{'trialType','ROI','mice'},'model','interaction','display','off');
  avgDffSumm.stats.(events{iEv}).multComp_ANOVA_accumul_correctVSerror_contra      ...
               = multcompare(avgDffSumm.stats.(events{iEv}).stats_ANOVA_accumul_correctVSerror_contra,'display','off');
  [avgDffSumm.stats.(events{iEv}).p_ANOVA_visGuide_correctVSerror_contra,avgDffSumm.stats.(events{iEv}).table_ANOVA_visGuide_correctVSerror_contra,avgDffSumm.stats.(events{iEv}).stats_ANOVA_visGuide_correctVSerror_contra] ...
               = anovan([avgVGc_contra(:); avgVGe_contra(:)],{conditvec,ROIvec,mousevec},'varnames',{'trialType','ROI','mice'},'model','interaction','display','off');
  avgDffSumm.stats.(events{iEv}).multComp_ANOVA_visGuide_correctVSerror_contra      ...
               = multcompare(avgDffSumm.stats.(events{iEv}).stats_ANOVA_visGuide_correctVSerror_contra,'display','off');
  [avgDffSumm.stats.(events{iEv}).p_ANOVA_accumulVSvisGuide_correct_contra,avgDffSumm.stats.(events{iEv}).table_ANOVA_accumulVSvisGuide_correct_contra,avgDffSumm.stats.(events{iEv}).stats_ANOVA_accumulVSvisGuide_correct_contra] ...
               = anovan([avgATc_contra(:); avgVGc_contra(:)],{conditvec,ROIvec,mousevec},'varnames',{'task','ROI','mice'},'model','interaction','display','off');
  avgDffSumm.stats.(events{iEv}).multComp_ANOVA_accumulVSvisGuide_correct_contra      ...
               = multcompare(avgDffSumm.stats.(events{iEv}).stats_ANOVA_accumulVSvisGuide_correct_contra,'display','off');
             
  % pairwise tests per ROI
  for iROI = 1:numel(lbls)
    switch avgDffSumm.stats.testname
      case 'signrank' 
        avgDffSumm.stats.(events{iEv}).p_accumul_correct_contraVSipsi(iROI,:)     = signrank(avgATc_contra(:,iROI),avgATc_ipsi(:,iROI));
        avgDffSumm.stats.(events{iEv}).p_accumul_error_contraVSipsi(iROI,:)       = signrank(avgATe_contra(:,iROI),avgATe_ipsi(:,iROI));
        avgDffSumm.stats.(events{iEv}).p_visGuide_correct_contraVSipsi(iROI,:)    = signrank(avgVGc_contra(:,iROI),avgVGc_ipsi(:,iROI));
        avgDffSumm.stats.(events{iEv}).p_visGuide_error_contraVSipsi(iROI,:)      = signrank(avgVGe_contra(:,iROI),avgVGe_ipsi(:,iROI));
        avgDffSumm.stats.(events{iEv}).p_accumulVSvisGuide_correct_contra(iROI,:) = signrank(avgATc_contra(:,iROI),avgVGc_contra(:,iROI));
        avgDffSumm.stats.(events{iEv}).p_accumulVSvisGuide_correct_ipsi(iROI,:)   = signrank(avgATc_ipsi(:,iROI),avgVGc_ipsi(:,iROI));
        avgDffSumm.stats.(events{iEv}).p_accumul_correctVSerror_contra(iROI,:)    = signrank(avgATc_contra(:,iROI),avgATe_contra(:,iROI));
        avgDffSumm.stats.(events{iEv}).p_accumul_correctVSerror_ipsi(iROI,:)      = signrank(avgATc_ipsi(:,iROI),avgATe_ipsi(:,iROI));
        avgDffSumm.stats.(events{iEv}).p_visGuide_correctVSerror_contra(iROI,:)   = signrank(avgVGc_contra(:,iROI),avgVGe_contra(:,iROI));
        avgDffSumm.stats.(events{iEv}).p_visGuide_correctVSerror_ipsi(iROI,:)     = signrank(avgVGc_ipsi(:,iROI),avgVGe_ipsi(:,iROI));
      case 'ttest'
        [~,avgDffSumm.stats.(events{iEv}).p_accumul_correct_contraVSipsi(iROI,:)]     = ttest(avgATc_contra(:,iROI),avgATc_ipsi(:,iROI));
        [~,avgDffSumm.stats.(events{iEv}).p_accumul_error_contraVSipsi(iROI,:)]       = ttest(avgATe_contra(:,iROI),avgATe_ipsi(:,iROI));
        [~,avgDffSumm.stats.(events{iEv}).p_visGuide_correct_contraVSipsi(iROI,:)]    = ttest(avgVGc_contra(:,iROI),avgVGc_ipsi(:,iROI));
        [~,avgDffSumm.stats.(events{iEv}).p_visGuide_error_contraVSipsi(iROI,:)]      = ttest(avgVGe_contra(:,iROI),avgVGe_ipsi(:,iROI));
        [~,avgDffSumm.stats.(events{iEv}).p_accumulVSvisGuide_correct_contra(iROI,:)] = ttest(avgATc_contra(:,iROI),avgVGc_contra(:,iROI));
        [~,avgDffSumm.stats.(events{iEv}).p_accumulVSvisGuide_correct_ipsi(iROI,:)]   = ttest(avgATc_ipsi(:,iROI),avgVGc_ipsi(:,iROI));
        [~,avgDffSumm.stats.(events{iEv}).p_accumul_correctVSerror_contra(iROI,:)]    = ttest(avgATc_contra(:,iROI),avgATe_contra(:,iROI));
        [~,avgDffSumm.stats.(events{iEv}).p_accumul_correctVSerror_ipsi(iROI,:)]      = ttest(avgATc_ipsi(:,iROI),avgATe_ipsi(:,iROI));
        [~,avgDffSumm.stats.(events{iEv}).p_visGuide_correctVSerror_contra(iROI,:)]   = ttest(avgVGc_contra(:,iROI),avgVGe_contra(:,iROI));
        [~,avgDffSumm.stats.(events{iEv}).p_visGuide_correctVSerror_ipsi(iROI,:)]     = ttest(avgVGc_ipsi(:,iROI),avgVGe_ipsi(:,iROI));

    end
  end
  
  % FDR correction
  pval_list = fieldnames(avgDffSumm.stats.(events{iEv}));
  for iTest = 1:numel(pval_list)
    if contains(pval_list{iTest},'ANOVA'); continue; end
    thisp = avgDffSumm.stats.(events{iEv}).(pval_list{iTest});
    avgDffSumm.stats.(events{iEv}).(['isSig_' pval_list{iTest}(3:end)]) = FDR(thisp,.05);
  end
end

%% ------------------------------------------------------------------------
%% save and plot
%% ------------------------------------------------------------------------
cd(rootdir)
avgDffSumm.cfg = cfg;
save avgDffSummary_triggered avgDffSumm cfg -v7.3

%% plot ROI avgs: correct, error for both mazes, as a function of evidence
nROI      = numel(avgDffSumm.ROIlbl);
wf        = widefieldParams;
[nr,nc]   = subplotOrg(nROI,4);
events    = {'turn','preturn','tower'};
lblsAll   = avgDffSumm.ROIlbl;
lbls      = unique(cellfun(@(x)(x(1:end-2)),lblsAll,'UniformOutput',false),'stable');
[cl,lbls] = getDefaultROIcl(lbls);

for iType = 1:numel(cfg.trialTypes)
  avgAT = []; avgVG = [];
  for iEv = 1:numel(events)
    avgAT(:,:,1)   = nanmean(avgDffSumm.(['ROIavg_accumul_' (cfg.trialTypes{iType})]).([events{iEv} 'R']),3);
    avgAT(:,:,2)   = nanmean(avgDffSumm.(['ROIavg_accumul_' (cfg.trialTypes{iType})]).([events{iEv} 'L']),3);
    avgVG(:,:,1)   = nanmean(avgDffSumm.(['ROIavg_visGuide_' (cfg.trialTypes{iType})]).([events{iEv} 'R']),3);
    avgVG(:,:,2)   = nanmean(avgDffSumm.(['ROIavg_visGuide_' (cfg.trialTypes{iType})]).([events{iEv} 'L']),3);
    
    if strcmp(events{iEv},'preturn')
      taxis  = avgDffSumm.timeAxis_preturn;
    else
      taxis  = avgDffSumm.timeAxis;
    end
  
    % ROI-by-ROI event-triggered average, towers
    figure;
    wf.applyFigDefaults(gcf,[nc nr],'w')
    for iROI = 1:nROI
      subplot(nr,nc,iROI); hold on
      plot([0 0],[-3 7],'--','color',[.5 .5 .5]);
      plot([taxis(1) taxis(end)],[0 0],'--','color',[.5 .5 .5])
      h(1) = plot(taxis,avgAT(:,iROI,1),'-','color',analysisParams.myBlue,'linewidth',1.5);
      h(2) = plot(taxis,avgAT(:,iROI,2),'-','color',analysisParams.myRed,'linewidth',1.5);
      wf.applyAxisDefaults(gca,'k')
      xlim([taxis(1) taxis(end)])
      ylim([-3 7])
      if iROI == 1
        wf.applyAxisLbls(gca,'time(s)','\DeltaF/F (z)',sprintf('%s - AT, %s, %s',lblsAll{iROI},events{iEv},cfg.trialTypes{iType}));
      else
        wf.applyAxisLbls(gca,[],[],avgDffSumm.ROIlbl{iROI});
      end
      if iROI == nROI
        legend(h,{'R','L'},'location','best')
      end
    end
    
    % single-panel summary, contra and ipsi, towers
    avgAT_contra = (avgAT(:,1:2:nROI,1) + avgAT(:,2:2:end,2))./2;
    avgAT_ipsi   = (avgAT(:,1:2:nROI,2) + avgAT(:,2:2:end,1))./2;
    figure;
    wf.applyFigDefaults(gcf,[4 2],'w')
    subplot(1,2,1); hold on 
    for iROI = 1:numel(lbls)
      plot(taxis, avgAT_contra(:,iROI),'-', 'LineWidth', 1.5, 'color', cl(iROI,:));
    end
    wf.applyAxisDefaults(gca,'k')
    xlim([taxis(1) 2])
    wf.applyAxisLbls(gca,'time(s)','\DeltaF/F (z)',sprintf('AT, contra %s, %s',events{iEv},cfg.trialTypes{iType}));
    subplot(1,2,2); hold on 
    for iROI = 1:numel(lbls)
      plot(taxis, avgAT_ipsi(:,iROI),'-', 'LineWidth', 1.5, 'color', cl(iROI,:));
    end
    wf.applyAxisDefaults(gca,'k')
    xlim([taxis(1) 2])
    legend(lbls,'location','best')
    wf.applyAxisLbls(gca,'time(s)','\DeltaF/F (z)',sprintf('AT, ipsi %s, %s',events{iEv},cfg.trialTypes{iType}));
    
    % ROI-by-ROI event-triggered average, visually-guided
    figure;
    wf.applyFigDefaults(gcf,[nc nr],'w')
    for iROI = 1:nROI
      subplot(nr,nc,iROI); hold on
      plot([0 0],[-3 7],'--','color',[.5 .5 .5]);
      plot([taxis(1) taxis(end)],[0 0],'--','color',[.5 .5 .5])
      h(1) = plot(taxis,avgVG(:,iROI,1),'-','color',analysisParams.myBlue,'linewidth',1.5);
      h(2) = plot(taxis,avgVG(:,iROI,2),'-','color',analysisParams.myRed,'linewidth',1.5);
      wf.applyAxisDefaults(gca,'k')
      xlim([taxis(1) taxis(end)])
      ylim([-3 7])
      if iROI == 1
        wf.applyAxisLbls(gca,'time(s)','\DeltaF/F (z)',sprintf('%s - VG, %s, %s',lblsAll{iROI},events{iEv},cfg.trialTypes{iType}));
      else
        wf.applyAxisLbls(gca,[],[],avgDffSumm.ROIlbl{iROI});
      end
      if iROI == nROI
        legend(h,{'R','L'},'location','best')
      end
    end
    
    % single-panel summary, contra and ipsi, vg
    avgVG_contra = (avgVG(:,1:2:nROI,1) + avgVG(:,2:2:end,2))./2;
    avgVG_ipsi   = (avgVG(:,1:2:nROI,2) + avgVG(:,2:2:end,1))./2;
    figure;
    wf.applyFigDefaults(gcf,[4 2],'w')
    subplot(1,2,1); hold on 
    for iROI = 1:numel(lbls)
      plot(taxis, avgVG_contra(:,iROI),'-', 'LineWidth', 1.5, 'color', cl(iROI,:));
    end
    wf.applyAxisDefaults(gca,'k')
    xlim([taxis(1) 2])
    wf.applyAxisLbls(gca,'time(s)','\DeltaF/F (z)',sprintf('VG, contra %s, %s',events{iEv},cfg.trialTypes{iType}));
    subplot(1,2,2); hold on 
    for iROI = 1:numel(lbls)
      plot(taxis, avgVG_ipsi(:,iROI),'-', 'LineWidth', 1.5, 'color', cl(iROI,:));
    end
    wf.applyAxisDefaults(gca,'k')
    xlim([taxis(1) 2])
    legend(lbls,'location','best')
    wf.applyAxisLbls(gca,'time(s)','\DeltaF/F (z)',sprintf('VG, ipsi %s, %s',events{iEv},cfg.trialTypes{iType}));
    
    % contra vs ipsi, within and across tasks
    if strcmp(events{iEv},'preturn')
      sumIdx  = find(avgDffSumm.timeAxis_preturn >= -cfg.summWin(2),1,'first'):find(avgDffSumm.timeAxis_preturn <= cfg.summWin(1),1,'last');
    else
      sumIdx  = find(avgDffSumm.timeAxis >= cfg.summWin(1),1,'first'):find(avgDffSumm.timeAxis <= cfg.summWin(2),1,'last');
    end
    ATcontra_timeAvg = nanmean(avgAT_contra(sumIdx,:),1);
    ATipsi_timeAvg   = nanmean(avgAT_ipsi(sumIdx,:),1);
    VGcontra_timeAvg = nanmean(avgVG_contra(sumIdx,:),1);
    VGipsi_timeAvg   = nanmean(avgVG_ipsi(sumIdx,:),1);
    
    figure;
    wf.applyFigDefaults(gcf,[2 2],'w')
    subplot(2,2,1); hold on
    for iROI = 1:numel(lbls)
      thish = plot([1 2], [ATcontra_timeAvg(iROI) ATipsi_timeAvg(iROI)],'o-', 'LineWidth', 1.5, 'color', cl(iROI,:));
      yl = get(gca,'ylim');
      text(2.05,yl(2)*.9,sprintf('P_cond = %1.1f\nP_ROI = %1.1f',avgDffSumm.stats.(events{iEv}).p_ANOVA_accumul_correct_contraVSipsi(1),avgDffSumm.stats.(events{iEv}).p_ANOVA_accumul_correct_contraVSipsi(2)))
      if avgDffSumm.stats.(events{iEv}).isSig_accumul_correct_contraVSipsi(iROI)
        set(thish, 'markerfacecolor', cl(iROI,:));
      end
    end
    wf.applyAxisDefaults(gca,'k')
    xlim([.75 2.25])
    set(gca,'xtick',1:2,'xticklabel',{'contra','ipsi'})
    rotateXLabels(gca,60)
    wf.applyAxisLbls(gca,[],'<\DeltaF/F (z)>',sprintf('AT, %s, %s',events{iEv},cfg.trialTypes{iType}));
    
    subplot(2,2,2); hold on
    for iROI = 1:numel(lbls)
      thish = plot([1 2], [VGcontra_timeAvg(iROI) VGipsi_timeAvg(iROI)],'o-', 'LineWidth', 1.5, 'color', cl(iROI,:));
      yl = get(gca,'ylim');
      text(2.05,yl(2)*.9,sprintf('P_cond = %1.1f\nP_ROI = %1.1f',avgDffSumm.stats.(events{iEv}).p_ANOVA_visGuide_correct_contraVSipsi(1),avgDffSumm.stats.(events{iEv}).p_ANOVA_visGuide_correct_contraVSipsi(2)))
      if avgDffSumm.stats.(events{iEv}).isSig_visGuide_correct_contraVSipsi(iROI)
        set(thish, 'markerfacecolor', cl(iROI,:));
      end
    end
    wf.applyAxisDefaults(gca,'k')
    xlim([.75 2.25])
    set(gca,'xtick',1:2,'xticklabel',{'contra','ipsi'})
    rotateXLabels(gca,60)
    wf.applyAxisLbls(gca,[],'<\DeltaF/F (z)>',sprintf('VG, %s, %s',events{iEv},cfg.trialTypes{iType}));
    
    subplot(2,2,3); hold on
    for iROI = 1:numel(lbls)
      thish = plot([1 2], [ATcontra_timeAvg(iROI) VGcontra_timeAvg(iROI)],'o-', 'LineWidth', 1.5, 'color', cl(iROI,:));
      yl = get(gca,'ylim');
      text(2.05,yl(2)*.9,sprintf('P_cond = %1.1f\nP_ROI = %1.1f',avgDffSumm.stats.(events{iEv}).p_ANOVA_accumulVSvisGuide_correct_contra(1),avgDffSumm.stats.(events{iEv}).p_ANOVA_accumulVSvisGuide_correct_contra(2)))
      if avgDffSumm.stats.(events{iEv}).isSig_accumulVSvisGuide_correct_contra(iROI)
        set(thish, 'markerfacecolor', cl(iROI,:));
      end
    end
    wf.applyAxisDefaults(gca,'k')
    xlim([.75 2.25])
    set(gca,'xtick',1:2,'xticklabel',{'AT','VG'})
    rotateXLabels(gca,60)
    wf.applyAxisLbls(gca,[],'<\DeltaF/F (z)>',sprintf('Contra %s, %s',events{iEv},cfg.trialTypes{iType}));
    
    subplot(2,2,4); hold on
    for iROI = 1:numel(lbls)
      thish = plot([1 2], [ATipsi_timeAvg(iROI) VGipsi_timeAvg(iROI)],'o-', 'LineWidth', 1.5, 'color', cl(iROI,:));
      if avgDffSumm.stats.(events{iEv}).isSig_accumulVSvisGuide_correct_ipsi(iROI)
        set(thish, 'markerfacecolor', cl(iROI,:));
      end
    end
    wf.applyAxisDefaults(gca,'k')
    xlim([.75 2.25])
    set(gca,'xtick',1:2,'xticklabel',{'AT','VG'})
    rotateXLabels(gca,60)
    wf.applyAxisLbls(gca,[],'<\DeltaF/F (z)>',sprintf('Ipsi %s, %s',events{iEv},cfg.trialTypes{iType}));
    
  end
end

%% correct vs error comparison
avgATc = []; avgVGc = [];
avgATe = []; avgVGe = [];
for iEv = 1:numel(events)
  avgATc(:,:,1)   = nanmean(avgDffSumm.ROIavg_accumul_correct.([events{iEv} 'R']),3);
  avgATc(:,:,2)   = nanmean(avgDffSumm.ROIavg_accumul_correct.([events{iEv} 'L']),3);
  avgVGc(:,:,1)   = nanmean(avgDffSumm.ROIavg_visGuide_correct.([events{iEv} 'R']),3);
  avgVGc(:,:,2)   = nanmean(avgDffSumm.ROIavg_visGuide_correct.([events{iEv} 'L']),3);
  avgATe(:,:,1)   = nanmean(avgDffSumm.ROIavg_accumul_correct.([events{iEv} 'R']),3);
  avgATe(:,:,2)   = nanmean(avgDffSumm.ROIavg_accumul_error.([events{iEv} 'L']),3);
  avgVGe(:,:,1)   = nanmean(avgDffSumm.ROIavg_visGuide_error.([events{iEv} 'R']),3);
  avgVGe(:,:,2)   = nanmean(avgDffSumm.ROIavg_visGuide_error.([events{iEv} 'L']),3);
  
  % single-panel summary, contra and ipsi, error vs correct, towers
  if strcmp(events{iEv},'preturn')
    sumIdx  = find(avgDffSumm.timeAxis_preturn >= -cfg.summWin(2),1,'first'):find(avgDffSumm.timeAxis_preturn <= cfg.summWin(1),1,'last'); 
  else
    sumIdx  = find(avgDffSumm.timeAxis >= cfg.summWin(1),1,'first'):find(avgDffSumm.timeAxis <= cfg.summWin(2),1,'last'); 
  end
  avgATc_contra = nanmean((avgATc(sumIdx,1:2:nROI,1) + avgATc(sumIdx,2:2:end,2))./2,1);
  avgATc_ipsi   = nanmean((avgATc(sumIdx,1:2:nROI,2) + avgATc(sumIdx,2:2:end,1))./2,1);
  avgATe_contra = nanmean((avgATe(sumIdx,1:2:nROI,1) + avgATe(sumIdx,2:2:end,2))./2,1);
  avgATe_ipsi   = nanmean((avgATe(sumIdx,1:2:nROI,2) + avgATe(sumIdx,2:2:end,1))./2,1);
  avgVGc_contra = nanmean((avgVGc(sumIdx,1:2:nROI,1) + avgVGc(sumIdx,2:2:end,2))./2,1);
  avgVGc_ipsi   = nanmean((avgVGc(sumIdx,1:2:nROI,2) + avgVGc(sumIdx,2:2:end,1))./2,1);
  avgVGe_contra = nanmean((avgVGe(sumIdx,1:2:nROI,1) + avgVGe(sumIdx,2:2:end,2))./2,1);
  avgVGe_ipsi   = nanmean((avgVGe(sumIdx,1:2:nROI,2) + avgVGe(sumIdx,2:2:end,1))./2,1);
 
  figure;
  wf.applyFigDefaults(gcf,[2 2],'w')
  subplot(2,2,1); hold on
  for iROI = 1:numel(lbls)
    thish = plot([1 2], [avgATc_contra(iROI) avgATe_contra(iROI)],'o-', 'LineWidth', 1.5, 'color', cl(iROI,:));
    yl = get(gca,'ylim');
    text(2.05,yl(2)*.9,sprintf('P_cond = %1.1f\nP_ROI = %1.1f',avgDffSumm.stats.(events{iEv}).p_ANOVA_accumul_correctVSerror_contra(1),avgDffSumm.stats.(events{iEv}).p_ANOVA_accumul_correctVSerror_contra(2)))
    if avgDffSumm.stats.(events{iEv}).isSig_accumul_correctVSerror_contra(iROI)
      set(thish, 'markerfacecolor', cl(iROI,:));
    end
  end
  wf.applyAxisDefaults(gca,'k')
  xlim([.75 2.25])
  set(gca,'xtick',1:2,'xticklabel',{'correct','error'})
  rotateXLabels(gca,60)
  wf.applyAxisLbls(gca ,[],'<\DeltaF/F (z)>',sprintf('AT, Contra %s',events{iEv}));

  subplot(2,2,2); hold on
  for iROI = 1:numel(lbls)
    thish = plot([1 2], [avgATc_ipsi(iROI) avgATe_ipsi(iROI)],'o-', 'LineWidth', 1.5, 'color', cl(iROI,:));
    if avgDffSumm.stats.(events{iEv}).isSig_accumul_correctVSerror_ipsi(iROI)
      set(thish, 'markerfacecolor', cl(iROI,:));
    end
  end
  wf.applyAxisDefaults(gca,'k')
  xlim([.75 2.25])
  set(gca,'xtick',1:2,'xticklabel',{'correct','error'})
  rotateXLabels(gca,60)
  wf.applyAxisLbls(gca,[],'<\DeltaF/F (z)>',sprintf('AT, Ipsi %s',events{iEv}));

  subplot(2,2,3); hold on
  for iROI = 1:numel(lbls)
    plot([1 2], [avgVGc_contra(iROI) avgVGe_contra(iROI)],'o-', 'LineWidth', 1.5, 'color', cl(iROI,:));
    yl = get(gca,'ylim');
    text(2.05,yl(2)*.9,sprintf('P_cond = %1.1f\nP_ROI = %1.1f',avgDffSumm.stats.(events{iEv}).p_ANOVA_visGuide_correctVSerror_contra(1),avgDffSumm.stats.(events{iEv}).p_ANOVA_visGuide_correctVSerror_contra(2)))
    if avgDffSumm.stats.(events{iEv}).isSig_visGuide_correctVSerror_contra(iROI)
      set(thish, 'markerfacecolor', cl(iROI,:));
    end
  end
  wf.applyAxisDefaults(gca,'k')
  xlim([.75 2.25])
  set(gca,'xtick',1:2,'xticklabel',{'correct','error'})
  rotateXLabels(gca,60)
  wf.applyAxisLbls(gca,[],'<\DeltaF/F (z)>',sprintf('VG, Contra %s',events{iEv}));

  subplot(2,2,4); hold on
  for iROI = 1:numel(lbls)
    thish = plot([1 2], [avgVGc_ipsi(iROI) avgVGe_ipsi(iROI)],'o-', 'LineWidth', 1.5, 'color', cl(iROI,:));
    if avgDffSumm.stats.(events{iEv}).isSig_visGuide_correctVSerror_ipsi(iROI)
      set(thish, 'markerfacecolor', cl(iROI,:));
    end
  end
  wf.applyAxisDefaults(gca,'k')
  xlim([.75 2.25])
  set(gca,'xtick',1:2,'xticklabel',{'correct','error'})
  rotateXLabels(gca,60)
  wf.applyAxisLbls(gca,[],'<\DeltaF/F (z)>',sprintf('VG, Ipsi %s',events{iEv}));
end

%% export
if ~isempty(dir('avgDffSummaryTriggered.pdf')); delete avgDffSummaryTriggered.pdf; end
figls = get(0,'children'); % print to pdf
for ii = length(figls):-1:1
  figure(figls(ii))
  export_fig avgDffSummaryTriggered.pdf -q101 -append
end
close all

catch ME
  displayException(ME)
end

end