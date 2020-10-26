function glmSumm = summarizeGLM(mice,spockFlag,cfg,concatSessFlag,visGuideFlag,autoRegrFlag)

% glmSumm = summarizeGLM(mice,spockFlag,cfg,concatSessFlag)
%           mice: cell array with mouse name list
%      spockFlag: true to run on spock
%            cfg: structure with analysis config, pass empty to load defaults    
% concatSessFlag: true to load fits to concatenated sessions for each
%                 animal (not recommended) 

%% initialize
if nargin < 1 || isempty(mice);      mice      = {'ai2';'ai3';'ai5';'ai7';'ai9';'ai10'}; end
if nargin < 2 || isempty(spockFlag); spockFlag = true;                           end
if nargin < 3;                       cfg       = struct([]);                     end
if nargin < 4;                       concatSessFlag = false;                     end
if nargin < 5;                       visGuideFlag   = false;                     end
if nargin < 6;                       autoRegrFlag   = false;                     end

cfg                = populateCfg(cfg);
cfg.concatSessFlag = concatSessFlag;
fn                 = sprintf('dffGLM_%s_%s',cfg.timeOrSpace,cfg.whichMethod);
if cfg.ROIflag
  fn = [fn '_ROI'];
end
if visGuideFlag
  fn = [fn '_visGuide'];
end
if autoRegrFlag
  fn = [fn '_corr_autoRegr'];
end
if cfg.concatSessFlag
  fn = [fn '_concatSess.mat'];
else
  fn = [fn '.mat'];
end


glmSumm.mice = mice;
glmSumm.cfg  = cfg;

try
%% compile from saved data files
tic
fprintf('compiling data...\n')
recls = cell(numel(mice),1);
wf    = widefieldParams;
for iMouse = 1:numel(mice)
  fprintf('\tmouse %d / %d',iMouse,numel(mice))
  if concatSessFlag
    rootdir       = wf.getRootDir(spockFlag);
    recls{iMouse} = {[rootdir mice{iMouse}]};
  else
    [recls{iMouse},rootdir] = getFullRecPath(mice{iMouse},spockFlag);
  end
  glmSumm.mouse(iMouse).recls = recls{iMouse};
  
  for iRec = 1:numel(recls{iMouse})
    fprintf('.')
    cd(recls{iMouse}{iRec})
    load(fn,'dffFit')
    if ~isfield(dffFit,'ROIlbl'); warning(sprintf('outdated rec: rerun %s',pwd)); end
    glmSumm.mouse(iMouse).accuracy(iRec,:)         = dffFit.accuracy;
    glmSumm.mouse(iMouse).accuracy_isSig(iRec,:)   = dffFit.isSig;
    glmSumm.mouse(iMouse).accuracy_shuffle(iRec,:) = nanmean(dffFit.shuffle.accuracy);
    
    weights = dffFit.weights;
    if cfg.ROIflag
      if cfg.zscoreWeights; weights = zscore(weights')'; end
      glmSumm.mouse(iMouse).weights(:,:,iRec) = weights;
    else
      glmSumm.mouse(iMouse).weights{iRec}     = weights;
    end
    
    glmSumm.mouse(iMouse).weightCorr(:,:,iRec) = corr(weights');
    
    if ~isfield(glmSumm,'glmCfg')
      glmSumm.glmCfg   = dffFit.cfg;
    end
    if ~isfield(glmSumm,'predLbls')
      glmSumm.predLbls = dffFit.predLbls;
    end
    if ~isfield(glmSumm,'ROIlbl')
      glmSumm.ROIlbl   = dffFit.ROIlbl;
    end

    if ~exist('xaxis','var')
      switch cfg.timeOrSpace
        case 'space'
          xaxis = dffFit.cfg.posBins;
        case 'time'
          xaxis = dffFit.cfg.timeBins;
      end
      glmSumm.xaxis = xaxis;
    end
    
    if cfg.ROIflag && ~exist('ROIlbl','var')      
      ROIlbl = dffFit.ROIlbl;
      cd(rootdir)
    end

  end
  fprintf('\n')
  
  %% average for this mouse, x-rec corr etc
  % include only significant fits
  thisacc                              = glmSumm.mouse(iMouse).accuracy;
  thissig                              = glmSumm.mouse(iMouse).accuracy_isSig;
  thisacc(~thissig)                    = nan;
  glmSumm.accuracy(iMouse,:)           = nanmean(thisacc);
  glmSumm.accuracy_sem(iMouse,:)       = nanstd(thisacc)./sqrt(iRec-1);
  glmSumm.accuracyShuffle(iMouse,:)    = nanmean(glmSumm.mouse(iMouse).accuracy_shuffle);
  glmSumm.accuracySuffle_sem(iMouse,:) = nanstd(glmSumm.mouse(iMouse).accuracy_shuffle)./sqrt(iRec-1);
  
  if cfg.ROIflag
    weightmat                          = glmSumm.mouse(iMouse).weights;
    sig                                = zeros(size(weightmat));
    for iRec = 1:size(weightmat,3)
      sig(:,:,iRec)                    = repmat(thissig(iRec,:)',[1 size(weightmat,2)]);
    end
    weightmat(~sig)                    = nan;
    glmSumm.weights(:,:,iMouse)        = nanmean(weightmat,3);
    glmSumm.weights_sem(:,:,iMouse)    = nanstd(weightmat,0,3)./sqrt(iRec-1);
    glmSumm.weightCorr(:,:,iMouse)     = nanmean(glmSumm.mouse(iMouse).weightCorr,3);
    glmSumm.weightCorr_sem(:,:,iMouse) = nanstd(glmSumm.mouse(iMouse).weightCorr,0,3)./sqrt(iRec-1);
    
    mat         = glmSumm.mouse(iMouse).weights;
    [nx,ny,nz]  = size(mat);
    mat         = reshape(mat,[nx*ny nz]);
    thiscorr                                         = corr(mat);
    thiscorr(triu(true(size(thiscorr)),0))           = nan; % since corr is symmetrical dont double dip ROI pair
    thiscorr                                         = thiscorr(:);
    glmSumm.stats.xRecWeightCorr.mouse(iMouse).pairs = thiscorr;
    glmSumm.stats.xRecWeightCorr_avg(iMouse,:)       = nanmean(thiscorr);
  end
end
fprintf('\tdone after %1.1 min\n',toc/60)

%% stats 
glmSumm.stats.weight_avg            = nanmean(glmSumm.weights,3);
glmSumm.stats.weight_sem            = nanstd(glmSumm.weights,0,3)./sqrt(numel(mice)-1);
glmSumm.stats.accuracy_avg          = nanmean(glmSumm.accuracy);
glmSumm.stats.accuracy_sem          = nanstd(glmSumm.accuracy)./sqrt(numel(mice)-1);
glmSumm.stats.accuracy_shuffle_avg  = nanmean(glmSumm.accuracyShuffle);
glmSumm.stats.accuracy_shuffle_sem  = nanstd(glmSumm.accuracyShuffle)./sqrt(numel(mice)-1);

allmat = [];
for iMouse = 1:numel(glmSumm.mouse)
  mat         = glmSumm.mouse(iMouse).weights;
  [nx,ny,nz]  = size(mat);
  mat         = reshape(mat,[nx*ny nz]);
  allmat(:,end+1:end+nz) = mat; 
end
allcc                                     = corr(allmat); % rec vs rec corr
allcc(triu(true(size(thiscorr)),0))       = nan;
allcc                                     = allcc(:);
glmSumm.stats.xRecWeightCorr_allpairs     = allcc(~isnan(allcc));

%% save
cd(rootdir)
fn = 'glmSummary';
if cfg.ROIflag
  fn = [fn '_ROI'];
end
if visGuideFlag
  fn = [fn '_visGuide'];
end
if autoRegrFlag
  fn = [fn '_corr_autoRegr'];
end
save(fn,'glmSumm')

%% plot
if cfg.ROIflag 
  %%
  wf      = widefieldParams;
  nROI    = size(glmSumm.weights,1);
  nRec    = size(glmSumm.weights,3);
  ROIlbl  = glmSumm.ROIlbl;
  if sum(strcmpi(glmSumm.glmCfg.predList,'ROI')) > 0
    [nr,nc] = subplotOrg(numel(glmSumm.glmCfg.predList)+1,4);
  else
    [nr,nc] = subplotOrg(numel(glmSumm.glmCfg.predList),4);
  end
  colors  = jet(nROI);
  figure;
  wf.applyFigDefaults(gcf,[nc nr+2],'w')
  
  for iPred = 1:numel(glmSumm.glmCfg.predList)
    if strcmpi(glmSumm.glmCfg.predList{iPred},'ROI'); continue; end
    isPred = arrayfun(@(x)(~isempty(strmatch(glmSumm.glmCfg.predList{iPred},x))),glmSumm.predLbls);
    nlags  = sum(isPred);
    
    subplot(nr,nc,iPred); hold on
    if nlags > 1
      switch glmSumm.cfg.timeOrSpace
        case 'time'
          lags = linspace(-glmSumm.glmCfg.predLagSec{iPred}(1),glmSumm.glmCfg.predLagSec{iPred}(2),nlags);
          xlbl = 'Lag (s)';
        case 'space'
          lags = linspace(-glmSumm.glmCfg.predLagCm{iPred}(1),glmSumm.glmCfg.predLagCm{iPred}(2),nlags);
          xlbl = 'Lag (cm)';
      end
      for iROI = 1:nROI
        plot(lags,squeeze(mean(glmSumm.weights(iROI,isPred,:),3)),'-','linewidth',1,'color',colors(iROI,:));
      end
      wf.applyAxisLbls(gca,xlbl,'Weight (a.u.)',glmSumm.glmCfg.predList{iPred})
    else
      for iROI = 1:nROI
        bar(iROI,squeeze(mean(glmSumm.weights(iROI,isPred,:),3)),'edgecolor',colors(iROI,:),'facecolor',colors(iROI,:));
        errorbar(iROI,squeeze(mean(glmSumm.weights(iROI,isPred,:),3)),               ...
                 squeeze(std(glmSumm.weights(iROI,isPred,:),0,3))/sqrt(nRec-1),      ...
                 '-','color',colors(iROI,:));
      end
      set(gca,'xtick',1:nROI,'xticklabel',ROIlbl)
      rotateXLabels(gca,60);
      wf.applyAxisLbls(gca,[],'Weight (a.u.)',glmSumm.glmCfg.predList{iPred})
    end
    
    wf.applyAxisDefaults(gca,'k'); axis tight
%     if iPred == 1; legend(ROIlbl,'location','southoutside','position',[.01 .5 .12 .3]); end
  end
  
  if sum(strcmpi(glmSumm.glmCfg.predList,'ROI')) > 0
    subplot(nr,nc,numel(cfg.predList)); hold on
  else
    subplot(nr,nc,iPred+1); hold on
  end
  for iROI = 1:nROI
    bar(iROI,mean(glmSumm.accuracy(:,iROI)),...
        'edgecolor',colors(iROI,:),'facecolor',colors(iROI,:));
    errorbar(iROI,mean(glmSumm.accuracy(:,iROI)),...
             std(glmSumm.accuracy(:,iROI))/sqrt(nRec-1),...
             '-','color',colors(iROI,:));
  end
  set(gca,'xtick',1:nROI,'xticklabel',ROIlbl)
  rotateXLabels(gca,90);
  wf.applyAxisLbls(gca,[],'r (data vs. pred)','Model g.o.f.')
  axis tight
  yl = get(gca,'ylim'); 
  text(nROI,yl(1)-diff(yl)*.4,upper(sprintf('GLM, %s, %s',cfg.whichMethod,cfg.timeOrSpace)),...
       'fontsize',13,'horizontalAlignment','right')
%   saveas(gcf,fn);
end
%%
saveas(gcf,fn); close

catch ME
  displayException(ME)
end
end

%% get config
function cfg = populateCfg(cfg)

if ~isfield(cfg,'ROIflag')
  cfg(1).ROIflag = true;
end
if ~isfield(cfg,'zscoreWeights')
  cfg(1).zscoreWeights = false;
end
if ~isfield(cfg,'whichMethod')
  cfg(1).whichMethod = 'ridge';%{'ridge';'lasso'};
end
if ~isfield(cfg,'trialType')
  cfg(1).trialType = 'correct';
end
if ~isfield(cfg,'timeOrSpace')
  cfg(1).timeOrSpace = 'space';
end


end

