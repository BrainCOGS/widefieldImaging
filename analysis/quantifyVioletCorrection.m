function quantifyVioletCorrection

% quantifyVioletCorrection;
% summarizes isosbestic correction data

%% analysis parameters
cfg.gcampMice      = {'ai2';'ai3';'ai5';'ai7';'ai9';'ai10'};
cfg.yfpMice        = {'ty1';'ty2'};
cfg.imageRegion{1} = {101:110,91:100}; % visual cortex
cfg.imageRegion{2} = {91:100,76:85}; % restrosplenial
cfg.imageRegion{3} = {46:55,81:90}; % mM2
cfg.imageRegion{4} = {26:35,91:100}; % aM2
cfg.recsPerAnimal  = 10; % some arbitrary high number to include all
cfg.egRecGCAMP     = 'ai10/20170814';%'ai7/20170814';
cfg.egRecYFP       = 'ty2/20170125';
cfg.egTimePts      = 13001:14200; % 2 min
cfg.histBins       = -.2:.002:.2; 
cfg.xcorr_taxis    = -5:.1:5; 
cfg.airpuffEgPath  = 'ai8/20170802_whiskersL_BV/';
BVrecs             = widefield_recLs.BVrecs;

wf      = widefieldParams;
rootdir = wf.getRootDir(isThisSpock);
BVrecs  = cellfun(@(x)(formatFilePath([rootdir x])),BVrecs,'Uniformoutput',false);
appath  = formatFilePath([rootdir cfg.airpuffEgPath]);
appath  = [appath 'violetCorrection.mat'];

%% get values for gcamp animals
gcamp_pxlVals      = cell(numel(cfg.gcampMice),cfg.recsPerAnimal,numel(cfg.imageRegion));
gcamp_pxlValsB     = cell(numel(cfg.gcampMice),cfg.recsPerAnimal,numel(cfg.imageRegion));
gcamp_pxlValsV     = cell(numel(cfg.gcampMice),cfg.recsPerAnimal,numel(cfg.imageRegion));
gcamp_egs_dff      = cell(numel(cfg.gcampMice),cfg.recsPerAnimal,numel(cfg.imageRegion));
gcamp_egs_dffB     = cell(numel(cfg.gcampMice),cfg.recsPerAnimal,numel(cfg.imageRegion));
gcamp_egs_dffV     = cell(numel(cfg.gcampMice),cfg.recsPerAnimal,numel(cfg.imageRegion));
gcamp_BVxcorr      = cell(numel(cfg.gcampMice),cfg.recsPerAnimal,numel(cfg.imageRegion));
gcamp_isPlotEg     = false(numel(cfg.gcampMice),cfg.recsPerAnimal,numel(cfg.imageRegion));

for iMouse = 1:numel(cfg.gcampMice)
  recIdx = cellfun(@(x)(~isempty(strfind(x,cfg.gcampMice{iMouse}))),BVrecs,'uniformOutput',false);
  recIdx = find([recIdx{:}],cfg.recsPerAnimal,'last'); 
  
  for iRec = 1:numel(recIdx)
    tic
    fprintf('analyzing %s...',BVrecs{recIdx(iRec)})
    cd(BVrecs{recIdx(iRec)})
    load dff dff 
    load rawf dffBlue dffViolet
    
    for iRegion = 1:numel(cfg.imageRegion)
      dffregion                   = squeeze(nanmean(nanmean(dff(cfg.imageRegion{iRegion}{1},cfg.imageRegion{iRegion}{2},:),1),2));
      dffregionB                  = squeeze(nanmean(nanmean(dffBlue(cfg.imageRegion{iRegion}{1},cfg.imageRegion{iRegion}{2},:),1),2));
      dffregionV                  = squeeze(nanmean(nanmean(dffViolet(cfg.imageRegion{iRegion}{1},cfg.imageRegion{iRegion}{2},:),1),2));
      gcamp_pxlVals{iMouse,iRec,iRegion}  = dffregion;
      gcamp_pxlValsB{iMouse,iRec,iRegion} = dffregionB;
      gcamp_pxlValsV{iMouse,iRec,iRegion} = dffregionV;
      gcamp_egs_dff{iMouse,iRec,iRegion}  = dffregion(cfg.egTimePts);
      gcamp_egs_dffB{iMouse,iRec,iRegion} = dffregionB(cfg.egTimePts);
      gcamp_egs_dffV{iMouse,iRec,iRegion} = dffregionV(cfg.egTimePts);
    end
    
    if ~isempty(strfind(BVrecs{recIdx(iRec)},cfg.egRecGCAMP))
      gcamp_isPlotEg(iMouse,iRec,1) = true;
    end
    try
      load BVxcorr xcorr_BvsV taxis
      gcamp_BVxcorr{iMouse,iRec,iRegion} = interp1(taxis',xcorr_BvsV,cfg.xcorr_taxis','linear','extrap');
    end
    clear dff dffBlue dffViolet
    fprintf(' done after %1.1f min\n',toc/60)
  end
  
end

%% get values for yfp animals
yfp_pxlVals      = cell(numel(cfg.yfpMice),cfg.recsPerAnimal,numel(cfg.imageRegion));
yfp_pxlValsB     = cell(numel(cfg.yfpMice),cfg.recsPerAnimal,numel(cfg.imageRegion));
yfp_pxlValsV     = cell(numel(cfg.yfpMice),cfg.recsPerAnimal,numel(cfg.imageRegion));
yfp_egs_dff      = cell(numel(cfg.yfpMice),cfg.recsPerAnimal,numel(cfg.imageRegion));
yfp_egs_dffB     = cell(numel(cfg.yfpMice),cfg.recsPerAnimal,numel(cfg.imageRegion));
yfp_egs_dffV     = cell(numel(cfg.yfpMice),cfg.recsPerAnimal,numel(cfg.imageRegion));
yfp_BVxcorr      = cell(numel(cfg.yfpMice),cfg.recsPerAnimal,numel(cfg.imageRegion));
yfp_isPlotEg     = false(numel(cfg.yfpMice),cfg.recsPerAnimal,numel(cfg.imageRegion));

for iMouse = 1:numel(cfg.yfpMice)
  recIdx = cellfun(@(x)(~isempty(strfind(x,cfg.yfpMice{iMouse}))),BVrecs,'uniformOutput',false);
  recIdx = find([recIdx{:}],cfg.recsPerAnimal,'last');
  
  for iRec = 1:numel(recIdx)
    tic
    fprintf('analyzing %s...',BVrecs{recIdx(iRec)})
    cd(BVrecs{recIdx(iRec)})
    load dff dff 
    load rawf dffBlue dffViolet
    
    for iRegion = 1:numel(cfg.imageRegion)
      dffregion                   = squeeze(nanmean(nanmean(dff(cfg.imageRegion{iRegion}{1},cfg.imageRegion{iRegion}{2},:),1),2));
      dffregionB                  = squeeze(nanmean(nanmean(dffBlue(cfg.imageRegion{iRegion}{1},cfg.imageRegion{iRegion}{2},:),1),2));
      dffregionV                  = squeeze(nanmean(nanmean(dffViolet(cfg.imageRegion{iRegion}{1},cfg.imageRegion{iRegion}{2},:),1),2));
      yfp_pxlVals{iMouse,iRec,iRegion}  = dffregion;
      yfp_pxlValsB{iMouse,iRec,iRegion} = dffregionB;
      yfp_pxlValsV{iMouse,iRec,iRegion} = dffregionV;
      yfp_egs_dff{iMouse,iRec,iRegion}  = dffregion(cfg.egTimePts);
      yfp_egs_dffB{iMouse,iRec,iRegion} = dffregionB(cfg.egTimePts);
      yfp_egs_dffV{iMouse,iRec,iRegion} = dffregionV(cfg.egTimePts);
    end
    
    if ~isempty(strfind(BVrecs{recIdx(iRec)},cfg.egRecYFP))
      yfp_isPlotEg(iMouse,iRec,1) = true;
    end
    
    load BVxcorr xcorr_BvsV taxis
    yfp_BVxcorr{iMouse,iRec}    = interp1(taxis',xcorr_BvsV,cfg.xcorr_taxis','linear','extrap');
    
    clear dff dffBlue dffViolet

    fprintf(' done after %1.1f min\n',toc/60)
  end
end

%% histograms and avg xcorr for overall distribution
gcamp = []; gcamp_xcorr = []; 
yfp   = []; yfp_xcorr   = [];
for iMouse = 1:numel(cfg.gcampMice)
  for iRec = 1:cfg.recsPerAnimal
    for iRegion = 1:numel(cfg.imageRegion)
      gcamp        = [gcamp; gcamp_pxlVals{iMouse,iRec,iRegion}];
      gcamp_xcorr  = [gcamp_xcorr; gcamp_BVxcorr{iMouse,iRec,iRegion}'./max(gcamp_BVxcorr{iMouse,iRec,iRegion})];
    end
  end
end
for iMouse = 1:numel(cfg.yfpMice)
  for iRec = 1:cfg.recsPerAnimal
    for iRegion = 1:numel(cfg.imageRegion)
      yfp          = [yfp; yfp_pxlVals{iMouse,iRec,iRegion}];
      yfp_xcorr    = [yfp_xcorr; yfp_BVxcorr{iMouse,iRec,iRegion}'./max(yfp_BVxcorr{iMouse,iRec,iRegion})];
    end
  end
end
hist_gcamp = histcounts(gcamp,cfg.histBins,'Normalization','probability');
hist_yfp   = histcounts(yfp,cfg.histBins,'Normalization','probability');

xcorr_mean_gcamp = nanmean(gcamp_xcorr);
xcorr_sem_gcamp  = nanstd(gcamp_xcorr)/sqrt(size(gcamp_xcorr,1));
xcorr_mean_yfp   = nanmean(yfp_xcorr);
xcorr_sem_yfp    = nanstd(yfp_xcorr)/sqrt(size(yfp_xcorr,1));

%% save
cd(rootdir)
clear iMouse iRec dff dffBlue dffViolet dffregion dffregionB dffregionV iRec iMouse gcamp yfp taxis xcorr_BvsV
save BVcorrectionQuantification

%% plot
wf = widefieldParams;

figure; wf.applyFigDefaults(gcf,[3 2],'w')
subplot(4,3,[1 2]); hold on
taxis = linspace(0,120,numel(cfg.egTimePts));
plot(taxis,gcamp_egs_dffB{gcamp_isPlotEg},'-','color',widefieldParams.myblue,'linewidth',1);
plot(taxis,gcamp_egs_dffV{gcamp_isPlotEg},'-','color',widefieldParams.mypurple,'linewidth',1);
legend({'470nm';'405nm'},'location','best'); legend('boxoff')
wf.applyAxisDefaults(gca,'k'); ylim([-.06 .06]); set(gca,'xtick',[],'xcolor','w')
wf.applyAxisLbls(gca,[],'\DeltaF/F','GCaMP6f')

subplot(4,3,[4 5]); hold on
plot(taxis,gcamp_egs_dff{gcamp_isPlotEg},'-','color',widefieldParams.darkgray,'linewidth',1);
legend({'corrected'},'location','best'); legend('boxoff')
wf.applyAxisDefaults(gca,'k'); ylim([-.06 .06]); set(gca,'xtick',[],'xcolor','w')
wf.applyAxisLbls(gca,[],'\DeltaF/F',[])

subplot(4,3,[7 8]); hold on
plot(taxis,yfp_egs_dffB{yfp_isPlotEg},'-','color',widefieldParams.myblue,'linewidth',1);
plot(taxis,yfp_egs_dffV{yfp_isPlotEg},'-','color',widefieldParams.mypurple,'linewidth',1);
wf.applyAxisDefaults(gca,'k'); ylim([-.06 .06]); set(gca,'xtick',[],'xcolor','w')
wf.applyAxisLbls(gca,[],'\DeltaF/F','YFP')

subplot(4,3,[10 11]); hold on
plot(taxis,yfp_egs_dff{yfp_isPlotEg},'-','color',widefieldParams.darkgray,'linewidth',1);
wf.applyAxisDefaults(gca,'k'); ylim([-.06 .06]); 
wf.applyAxisLbls(gca,'Time (s)','\DeltaF/F',[])

subplot(4,3,[3 6]); hold on
plot(cfg.xcorr_taxis, xcorr_mean_gcamp,'-','color',widefieldParams.mediumgreen,'linewidth',1)
plot(cfg.xcorr_taxis, xcorr_mean_yfp,'-','color',widefieldParams.darkgray,'linewidth',1)

legend({'GCaMP';'YFP'},'location','southwest'); legend('boxoff');

plot(cfg.xcorr_taxis, xcorr_mean_gcamp-xcorr_sem_gcamp,'--','color',widefieldParams.mediumgreen,'linewidth',.5)
plot(cfg.xcorr_taxis, xcorr_mean_gcamp+xcorr_sem_gcamp,'--','color',widefieldParams.mediumgreen,'linewidth',.5)

plot(cfg.xcorr_taxis, xcorr_mean_yfp-xcorr_sem_yfp,'--','color',widefieldParams.darkgray,'linewidth',.5)
plot(cfg.xcorr_taxis, xcorr_mean_yfp+xcorr_sem_yfp,'--','color',widefieldParams.darkgray,'linewidth',.5)
xlim([cfg.xcorr_taxis(1) cfg.xcorr_taxis(end)]); ylim([-.25 1])
xl = get(gca,'xlim'); yl = get(gca,'ylim');
plot(xl,[0 0],'--','color',widefieldParams.lightgray);
plot([0 0],yl,'--','color',widefieldParams.lightgray);
wf.applyAxisDefaults(gca,'k'); 
wf.applyAxisLbls(gca,'Blue -> Violet Lag (s)','Norm. x-corr',[])

subplot(4,3,[9 12]); hold on
load(appath)
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
set(gca,'ytick',-.01:.01:.02)
wf.applyAxisDefaults(gca,'k'); 
wf.applyAxisLbls(gca,'Time from airpuff (s)','\DeltaF/F',[])

%%
saveas(gcf,'violetCorrection')
export_fig violetCorrection.pdf
close