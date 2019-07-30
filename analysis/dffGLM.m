function dffFit = dffGLM(dff,logSumm,ROIflag,derivFlag,cfg,ROIlbl,saveFlag,visGuideFlag)

% dffFit = dffGLM(dff,logSumm,ROIflag,derivFlag,cfg,ROIlbl,saveFlag,visGuideFlag)
% fits GLMs to the activity of each ROI
%
% input:
%   dff: can be either string with recording path or a frames x ROI matrix
%        (default: pwd)
%   ROIflag: true for ROIs, false for pxls. 
%   derivFlag: true to take derivative of dff before decoding, as a dirty form
%              of "deconvolution" (not recommended)(default: false)
%   cfg: structure with analysis config, refer to bottom of function. If
%        empty will be filled with defaults (recommended)
%   ROIlbl: for ROIs, cell array with names (optional, if empty or not an
%           input labels will be loaded from disk)
%   saveFlag: true to save analysis results (default)
%   visGuideFlag: true will analyze visuall-guided instead of towers task
%   (default false)
%
% output is a data structure with analysis results

%%
if nargin < 2; logSumm      = [];         end
if nargin < 3; ROIflag      = true;       end
if nargin < 4; derivFlag    = false;      end
if nargin < 5; cfg          = struct([]); end
if nargin < 6; ROIlbl       = {};         end
if nargin < 7; saveFlag     = true;       end
if nargin < 8; visGuideFlag = false;      end

%%
if ischar(dff)
  cd(dff)
  if ROIflag
    load dffROI dffROI
    dff = dffROI;
    clear dffROI
  else
    load dff dff
  end
end
if isempty(logSumm)
  load behavLog logSumm
end

%% fn, config
cfg(1).diffFlag  = derivFlag;
cfg              = populateCfg(cfg);
fn               = sprintf('dffGLM_%s_%s',cfg.timeOrSpace,cfg.method);
if ROIflag;        fn  = [fn '_ROI']; end

if visGuideFlag
  mz = 4;
  fn = [fn '_visGuide'];
else
  mz = [];
end
cfg.visGuideFlag = visGuideFlag;

tic; fprintf('FITTING %s\n\tBuilding matrices...',upper(fn));

try
%% ROI labels
if ROIflag && isempty(ROIlbl)
  load dffROI ROIlbl
end
dffFit.ROIlbl = ROIlbl;

%% dff
if ~ROIflag && cfg.spatialBinFactor > 1
  dff       = imBinSpace(dff,[],cfg.spatialBinFactor,true);
  fn        = sprintf('%s_binned%dx',fn,cfg.spatialBinFactor);
end
if size(dff,3) > 1
  [dff,bc] = conditionDffMat(dff); %#ok<ASGLU>
end

if cfg.diffFlag % fast dirty deconvolution by derivative
  dff = [zeros(1,size(dff,2)); diff(dff,1,1)];
  fn  = [fn '_deriv']; 
end
if cfg.zscoreFlag; dff = zscore(dff); end

  
%% align trials (with parallelization)
poolobj = gcp('nocreate');
if isempty(poolobj); poolobj = parpool; end
addAttachedFiles(poolobj, 'behavLog.mat')

%% align trials (with parallelization)
dffTrialsTemp = cell(1,size(dff,2));

tt = cfg.trialType;
ap = cfg.alignPoint;
tb = cfg.timeBins;
pb = cfg.posBins;
ts = cfg.timeOrSpace;
switch cfg.timeOrSpace
  case 'time'
    [dffTrials,~,taxis,trialidx] = alignTrials(dff(:,1),logSumm,cfg.trialType,mz,cfg.alignPoint,cfg.timeBins,cfg.timeOrSpace);
    parfor iPxl = 2:size(dff,2)
      dffTrialsTemp{iPxl}    = alignTrials(dff(:,iPxl),[],tt,mz,ap,tb,ts);
    end
  case 'space'
    [dffTrials,~,taxis,trialidx] = alignTrials(dff(:,1),logSumm,cfg.trialType,mz,[],cfg.posBins,cfg.timeOrSpace);
    parfor iPxl = 2:size(dff,2)
      dffTrialsTemp{iPxl}    = alignTrials(dff(:,iPxl),[],tt,mz,[],pb,ts);
    end
end

% put back in matrix format
for iPxl = 2:size(dff,2)
 dffTrials = cat(3,dffTrials,dffTrialsTemp{iPxl});
end
clear dffTrialsTemp

%% generate X and y matrices for regression
if ROIflag
  [X, y, trialID, predLbls] = generateMats(dffTrials,logSumm,trialidx,taxis,cfg);
else
  [X, y, trialID, predLbls] = generateMats(dffTrials(:,:,1:2),logSumm,trialidx,taxis,cfg);
end
dffFit.predLbls             = predLbls;
dffFit.cfg                  = cfg;
nROI                        = size(dff,2);
dffFit.nROI                 = nROI;

fprintf(' done after %1.1f min\n',toc/60)

%% manually generate cross validation runs, sample by trial instead of time point
tic; fprintf('\tfitting and cross-validating models ...');
pseudoExp                 = generateXval(trialID,cfg);

%% fit / cross-validate model
jitter                = cfg.stateJitter;
cl                    = cfg.confidLevel;
l2                    = cfg.l2Penalty;
lambda                = cfg.lambdaRatio;
method                = cfg.method;

bestCoeff             = cell(1,dffFit.nROI);
accuracy              = cell(1,dffFit.nROI);
accInterval           = cell(1,dffFit.nROI);
coeff                 = cell(1,dffFit.nROI);
fitInfo               = cell(1,dffFit.nROI);
bestpred_y            = cell(1,dffFit.nROI);
bestpred_yhat         = cell(1,dffFit.nROI);

parfor iROI = 1:nROI %#ok<PFUIX>
%   try
  if ~ROIflag
    [pred, target]    = generateMats(dffTrials(:,:,[iROI iROI]),logSumm,trialidx,taxis,cfg);
    if size(target,2) > 1
      target          = target(:,1);
    end
    pred              = pred(:,:,1);
  else
    if size(y,2) > 1
      target          = y(:,iROI);
    else
      target          = y;
    end
    pred              = X(:,:,iROI);
  end
  
  % add jitter to all-zero pixels
  allzeros = sum(pred==0) == size(pred,1);
  pred(:,allzeros) = pred(:,allzeros) + jitter(1) * randn(size(pred(:,allzeros)));

  switch method
    case 'lasso'
      [coeff{iROI}, fitInfo{iROI}] ...
                     = elasticNetRegression_LP( pred, target, [], ...
                                            {'LeastR',jitter}, l2, lambda, pseudoExp, [], [] ); 
      fitInfo{iROI} = rmfield(fitInfo{iROI},{'X','y','w'});

    case 'ridge'
      [coeff{iROI}, fitInfo{iROI}] ...
                     = ridgeRegression( pred, target, lambda, pseudoExp, [], jitter); 
  end

  bestCoeff{iROI}  = coeff{iROI}(:,fitInfo{iROI}.Index1SE);

  % Compute accuracy as the fraction of correct predictions, with spreads across CV samples
  predAccuracy      = nan(size(fitInfo{iROI}.CVTestPrediction));

  % Cross-validated prediction accuracy
  prediction        = cell(size(fitInfo{iROI}.CVTestPrediction,1),size(fitInfo{iROI}.CVTestPrediction,2));
  predTruth         = prediction;
  for iMC = 1:size(fitInfo{iROI}.CVTestPrediction,1)
    for iFold = 1:size(fitInfo{iROI}.CVTestPrediction,2)
      % Compute optimal threshold for separating training set populations according to truth
      if isfield(fitInfo{iROI}.CVExperiments(iMC).test(iFold),'idx')
        testSel     = fitInfo{iROI}.CVExperiments(iMC).test(iFold).idx;
      else
        testSel     = fitInfo{iROI}.CVExperiments(iMC).test(iFold);
      end
      thispred                 = fitInfo{iROI}.CVTestPrediction{iMC,iFold}(:,fitInfo{iROI}.Index1SE);
      thistruth                = target(testSel);
      prediction{iMC,iFold}    = thispred;
      predTruth{iMC,iFold}     = thistruth;
      nanidx                   = isnan(thispred) | isnan(thistruth);
      if sum(nanidx) == numel(thispred)
        predAccuracy(iMC,iFold)  = nan;
      else
        predAccuracy(iMC,iFold)  = corr(thispred(~nanidx),thistruth(~nanidx));
      end
    end
  end
  
  
  accuracy{iROI}      = nanmean(predAccuracy(:));
  accInterval{iROI}   = quantile(predAccuracy(:), cl); 
  [row,col]           = find(predAccuracy==max(predAccuracy(:)),1,'first');
  bestpred_yhat{iROI} = prediction{row,col};
  bestpred_y{iROI}    = predTruth{row,col};
%   catch
%     fprintf('problem with ROI # %d',iROI)
%   end
end
fprintf(' done after %1.1f min\n',toc/60)

%% Consolidate values
dffFit.bestpred_y    = bestpred_y;
dffFit.bestpred_yhat = bestpred_yhat;
dffFit.allCoeffs     = coeff;
dffFit.fitInfo       = fitInfo;
dffFit.accuracy      = cell2mat(accuracy);
dffFit.weights       = cell2mat(bestCoeff)';

for iROI = 1:nROI
  dffFit.accuracyCI(iROI,:) = accInterval{iROI};
end

clear bestCoeff accuracy coeff 

%% do shuffling if necessary
if cfg.numShuffles > 1
  
  tic;
  fprintf('\tfitting and cross-validating models (shuffles)...');
  bestCoeff             = cell(cfg.numShuffles,nROI);
  accuracy              = cell(cfg.numShuffles,nROI);
  all_accuracy          = cell(nROI,1);
  numShuff              = cfg.numShuffles;
  numCVSamples          = cfg.numCVSamplesShuffle; % in each shuffle just do 5 runs of xval
  iPseudo               = pseudoExp(1:numCVSamples);
%   numCVFolds            = cfg.numCVFolds;
  
  
  % use best lambda from actual data
  bestLambdas           = zeros(1,iROI);
  for iROI = 1:nROI
    bestLambdas(iROI)  = lambda([fitInfo{iROI}.IndexMinMSE]);
  end

  parfor iROI = 1:nROI
    if ~ROIflag
      [pred, target]    = generateMats(dffTrials(:,:,[iROI iROI]),logSumm,trialidx,taxis,cfg);
      if size(target,2) > 1
        target          = target(:,1);
      end
      pred              = pred(:,:,1);
    else
      if size(y,2) > 1
        target          = y(:,iROI);
      else
        target          = y;
      end
      pred              = X(:,:,iROI);
    end
    
    % add jitter to all-zero pixels    
    allzeros = sum(pred==0) == size(pred,1);
    pred(:,allzeros) = pred(:,allzeros) + jitter(1) * randn(size(pred(:,allzeros))); %#ok<PFBNS>
    
    for iShuff = 1:numShuff 
      
      iTarget         = target(randperm(numel(target)));
      
      switch method 
        case 'lasso'
          [bestCoeff{iShuff,iROI}, sfitInfo] = elasticNetRegression_LP( pred, iTarget, [],                         ...
                                                                      {'LeastR',jitter}, l2, bestLambdas(iROI), ...
                                                                      iPseudo, [], [] );
          
        case 'ridge'
          [bestCoeff{iShuff,iROI}, sfitInfo] = ridgeRegression( pred, iTarget, bestLambdas(iROI), ...
                                                                 iPseudo, [], jitter);
          
      end

      % Compute accuracy as the fraction of correct predictions, with spreads across CV samples
      predAccuracy      = nan(size(sfitInfo.CVTestPrediction));
      
      % Cross-validated prediction accuracy
      for iMC = 1:size(sfitInfo.CVTestPrediction,1)
        for iFold = 1:size(sfitInfo.CVTestPrediction,2)
          % Compute accuracy by computing correlation between prediction
          % and data
          if isfield(sfitInfo.CVExperiments(iMC).test(iFold),'idx')
            testSel   = sfitInfo.CVExperiments(iMC).test(iFold).idx;
          else
            testSel   = sfitInfo.CVExperiments(iMC).test(iFold);
          end
          thispred                 = sfitInfo.CVTestPrediction{iMC,iFold}(:,sfitInfo.Index1SE);
          thistruth                = target(testSel);
          nanidx                   = isnan(thispred) | isnan(thistruth);
          if sum(nanidx) == numel(thispred)
            predAccuracy(iMC,iFold)  = nan;
          else
            predAccuracy(iMC,iFold)  = corr(thispred(~nanidx),thistruth(~nanidx));
          end
        end
      end
      
      accuracy{iShuff,iROI} = nanmean(predAccuracy(:));
      all_accuracy{iROI}    = [all_accuracy{iROI}; predAccuracy(:)];
    end
    all_accuracy{iROI}      = all_accuracy{iROI}(~isnan(all_accuracy{iROI}));
  end
  
  dffFit.shuffle.coeffs        = cell2mat(bestCoeff);
  dffFit.shuffle.accuracy      = cell2mat(accuracy);
  
  for iROI = 1:nROI
    dffFit.isSig(iROI) = dffFit.accuracy(iROI) > prctile(all_accuracy{iROI},95);
  end
  
  fprintf(' done after %1.1f min\n',toc/60)
end

%% shut down parallel pool
delete(poolobj);

%% save 
if ~saveFlag; return; end
save(fn,'dffFit','-v7.3')
if ~ROIflag; save(fn,'bc','-append'); end

%% plot
if ROIflag && cfg.plotFlag
  
  wf      = widefieldParams;
  if sum(strcmpi(dffFit.cfg.predList,'ROI')) > 0
    [nr,nc] = subplotOrg(numel(dffFit.cfg.predList)+1,4);
  else
    [nr,nc] = subplotOrg(numel(dffFit.cfg.predList),4);
  end
  colors  = jet(nROI);
  figure;
  wf.applyFigDefaults(gcf,[nc nr+2],'w')
  
  for iPred = 1:numel(dffFit.cfg.predList)
    if strcmpi(dffFit.cfg.predList{iPred},'ROI'); continue; end
    isPred = arrayfun(@(x)(~isempty(strmatch(dffFit.cfg.predList{iPred},x))),dffFit.predLbls);
    nlags  = sum(isPred);
    
    subplot(nr,nc,iPred); hold on
    if nlags > 1
      switch dffFit.cfg.timeOrSpace
        case 'time'
          lags = linspace(-dffFit.cfg.predLagSec{iPred}(1),dffFit.cfg.predLagSec{iPred}(2),nlags);
          xlbl = 'Lag (s)';
        case 'space'
          lags = linspace(-dffFit.cfg.predLagCm{iPred}(1),dffFit.cfg.predLagCm{iPred}(2),nlags);
          xlbl = 'Lag (cm)';
      end
      for iROI = 1:nROI
        plot(lags,dffFit.weights(iROI,isPred),'-','linewidth',1,'color',colors(iROI,:));
      end
      wf.applyAxisLbls(gca,xlbl,'Weight (a.u.)',dffFit.cfg.predList{iPred})
    else
      for iROI = 1:nROI
        bar(iROI,dffFit.weights(iROI,isPred),'edgecolor',colors(iROI,:),'facecolor',colors(iROI,:));
      end
      set(gca,'xtick',1:nROI,'xticklabel',ROIlbl)
      rotateXLabels(gca,60);
      wf.applyAxisLbls(gca,[],'Weight (a.u.)',dffFit.cfg.predList{iPred})
    end
    
    wf.applyAxisDefaults(gca,'k'); axis tight
%     if iPred == 1; legend(ROIlbl,'location','southoutside','position',[.01 .5 .12 .3]); end
  end
  
  if sum(strcmpi(dffFit.cfg.predList,'ROI')) > 0
    subplot(nr,nc,numel(cfg.predList)); hold on
  else
    subplot(nr,nc,iPred+1); hold on
  end
  for iROI = 1:nROI
    bar(iROI,dffFit.accuracy(iROI),'edgecolor',colors(iROI,:),'facecolor',colors(iROI,:));
  end
  set(gca,'xtick',1:nROI,'xticklabel',ROIlbl)
  rotateXLabels(gca,90);
  wf.applyAxisLbls(gca,[],'r (data vs. pred)','Model g.o.f.')
  axis tight
  yl = get(gca,'ylim'); 
  text(nROI,yl(1)-diff(yl)*.4,upper(sprintf('GLM, %s, %s',cfg.method,cfg.timeOrSpace)),...
       'fontsize',13,'horizontalAlignment','right')
  saveas(gcf,fn);
end

catch ME
%   keyboard
  displayException(ME)
end

end

%% ------------------------------------------------------------------------
%% GENERATE X and y
function [X,y,trialID,predLbls] = generateMats(dffTrials,logSumm,trialidx,taxis,cfg)

X         = [];
y         = [];
trialID   = [];
hasReward = strcmpi(cfg.predList,'rw');
predLbls  = {};

switch cfg.timeOrSpace
  case 'time'
    for iTrial = 1:numel(trialidx)
      %% chop trial to last valid time point, compile y for each ROI
      if sum(hasReward) > 1
        lastt = logSumm.keyFrames{trialidx(iTrial)}(end) + round(cfg.predLagSec{hasReward}(2)/cfg.logSumm.frameDtCam);
      else
        lastt = logSumm.keyFrames{trialidx(iTrial)}(end); % reward time is last point if no rw predictor
      end
      
      idx     = 1:find(taxis < logSumm.time{trialidx(iTrial)}(lastt), 1, 'last');
      thisy   = squeeze(dffTrials(iTrial,idx,:)); 
      trialID(end+1:end+numel(idx),:) = ones(numel(idx),1).*iTrial;
      y(end+1:end+numel(idx),:)       = thisy;
      
      %% add predictors to X, with lags if necessary
      thisX = [];
      for iPred = 1:numel(cfg.predList)
        switch cfg.predList{iPred}
          case {'tow_R','tow_L'}
            if ~isempty(strfind(cfg.predList{iPred},'L'))
              towers = logSumm.cueOnset_L{trialidx(iTrial)};
            else
              towers = logSumm.cueOnset_R{trialidx(iTrial)};
            end
            
            toweridx  = arrayfun(@(x)(find(taxis >= x,1,'first')),towers);
            thisvec   = zeros(numel(idx),1);
            thisvec(toweridx) = 1;
          
          case 'y'
            thisvec   = zeros(numel(idx),1);
            for iT = 1:numel(thisvec)
              iter        = find(logSumm.time{trialidx(iTrial)} >= taxis(iT),1,'first'): ...
                            min([size(logSumm.pos{trialidx(iTrial)},1) ...
                            find(logSumm.time{trialidx(iTrial)} < taxis(iT+1),1,'last')]);
              thisvec(iT) = mean(logSumm.pos{trialidx(iTrial)}(iter,2));
            end
            
          case 'speed'
            thisvec   = zeros(numel(idx),1);
            for iT = 1:numel(thisvec)
              iter        = find(logSumm.time{trialidx(iTrial)} >= taxis(iT),1,'first'): ...
                            min([size(logSumm.displ{trialidx(iTrial)},1) ...
                            find(logSumm.time{trialidx(iTrial)} < taxis(iT+1),1,'last')]);
              displ       = mean(sqrt(sum(logSumm.displ{trialidx(iTrial)}(iter,1:2).^2,2))); 
              thisvec(iT) = displ / logSumm.frameRateVirmen;
            end
          
          case 'd\theta/dt'
            thisvec   = zeros(numel(idx),1);
            for iT = 1:numel(thisvec)
              iter        = find(logSumm.time{trialidx(iTrial)} >= taxis(iT),1,'first'): ...
                            min([size(logSumm.displ{trialidx(iTrial)},1) ...
                            find(logSumm.time{trialidx(iTrial)} < taxis(iT+1),1,'last')]);
              displ       = mean(logSumm.displ{trialidx(iTrial)}(iter,3)); 
              thisvec(iT) = displ / logSumm.frameRateVirmen;
            end
            
          case '\theta'
            thisvec   = zeros(numel(idx),1);
            for iT = 1:numel(thisvec)
              iter        = find(logSumm.time{trialidx(iTrial)} >= taxis(iT),1,'first'): ...
                            min([size(logSumm.pos{trialidx(iTrial)},1) ...
                            find(logSumm.time{trialidx(iTrial)} < taxis(iT+1),1,'last')]);
              thisvec(iT) = mean(logSumm.pos{trialidx(iTrial)}(iter,3));
            end
          
          case '\Delta_bins'
            thisvec   = zeros(numel(idx),numel(cfg.evidenceBins)-1);
            for iT = 1:size(thisvec,1)
              thisev             = sum(logSumm.cueOnset_R{trialidx(iTrial)} <= taxis(iT)) - ...
                                   sum(logSumm.cueOnset_L{trialidx(iTrial)} <= taxis(iT));
              binidx             = find(cfg.evidenceBins <= thisev,1,'last');       
              thisvec(iT,binidx) = 1;            
            end
            
          case '\Delta'
            thisvec   = zeros(numel(idx),1);
            for iT = 1:numel(thisvec)
              thisvec(iT) = sum(logSumm.cueOnset_R{trialidx(iTrial)} <= taxis(iT)) - ...
                            sum(logSumm.cueOnset_L{trialidx(iTrial)} <= taxis(iT));
            end
          
          case 'abs(\Delta)'
            thisvec   = zeros(numel(idx),1);
            for iT = 1:numel(thisvec)
              thisvec(iT) = abs(sum(logSumm.cueOnset_R{trialidx(iTrial)} <= taxis(iT)) - ...
                                sum(logSumm.cueOnset_L{trialidx(iTrial)} <= taxis(iT)));
            end
            
          case 'ch'
            isLeft  = logSumm.choice(trialidx(iTrial)) == analysisParams.leftCode;
            thisvec = ones(numel(idx),1) - 2.*ones(numel(idx),1).*isLeft;
            
          case 'prevch'
            if trialidx(iTrial) == 1
              thisvec = zeros(numel(idx),1);
            else
              isLeft  = logSumm.choice(trialidx(iTrial)-1) == analysisParams.leftCode;
              thisvec = ones(numel(idx),1) - 2.*ones(numel(idx),1).*isLeft;
            end
            
          case 'prevrw'
            if trialidx(iTrial) == 1
              thisvec = zeros(numel(idx),1);
            else
              wasRW   = logSumm.choice(trialidx(iTrial)-1) == logSumm.trialType(trialidx(iTrial)-1);
              thisvec = ones(numel(idx),1) - 2.*ones(numel(idx),1).*wasRW;
            end
            
          case 'rw'
            thisvec = zeros(numel(idx),1);
            rwt     = logSumm.time{trialidx(iTrial)}(logSumm.keyFrames{trialidx(iTrial)}(end));
            idx     = find(taxis >= rwt,1,'first');
            thisvec(idx) = 1;
            
          case 'ROI'
            continue
        end
        
        % add lags
        nlags     = round(cfg.predLagSec{iPred} / mode(diff(taxis))); 
        lagct     = 1;
        vecs      = zeros(size(thisvec,1),sum(nlags)+1);
        for iLag = -nlags(1):nlags(2)
          if iLag < 0
            shiftvec = [thisvec(-iLag+1:end); zeros(abs(iLag),1)];
          elseif iLag > 0
            shiftvec = [zeros(abs(iLag),1); thisvec(1:end-iLag)];
          else
            shiftvec = thisvec;
          end
          
          vecs(:,lagct) = shiftvec;
          if iTrial == 1
            predLbls{end+1} = [cfg.predList{iPred} num2str(iLag)];
          end
          lagct = lagct+1;
        end
        
        thisX(:,end+1:end+size(vecs,2)) = vecs; 
      end
      
      X(end+1:end+size(thisX,1),:) = thisX;
      
    end
    
  case 'space'
    for iTrial = 1:numel(trialidx)
      
      idx     = 1:numel(cfg.posBins(1:end-1));
      taxis   = cfg.posBins;
      thisy   = squeeze(dffTrials(iTrial,idx,:));
      trialID(end+1:end+numel(idx),:) = ones(numel(idx),1).*iTrial;
      y(end+1:end+numel(idx),:)       = thisy;
      
      %% add predictors to X, with lags if necessary
      thisX = [];
      for iPred = 1:numel(cfg.predList)
        switch cfg.predList{iPred}
          case {'tow_R','tow_L'}
            if ~isempty(strfind(cfg.predList{iPred},'L'))
              towers = logSumm.cuePos_L{trialidx(iTrial)};
            else
              towers = logSumm.cuePos_R{trialidx(iTrial)};
            end
            
            toweridx  = arrayfun(@(x)(find(taxis >= x,1,'first')),towers);
            thisvec   = zeros(numel(idx),1);
            thisvec(toweridx) = 1;
            
          case 'speed'
            thisvec   = zeros(numel(idx),1);
            for iT = 1:numel(thisvec)
              iter        = find(logSumm.pos{trialidx(iTrial)}(:,2) >= taxis(iT),1,'first'): ...
                            find(logSumm.pos{trialidx(iTrial)}(:,2) < taxis(iT+1),1,'last');
              displ       = mean(sqrt(sum(logSumm.displ{trialidx(iTrial)}(iter,1:2).^2,2))); 
              thisvec(iT) = displ / logSumm.frameRateVirmen;
            end
          
          case 'd\theta/dt'
            thisvec   = zeros(numel(idx),1);
            for iT = 1:numel(thisvec)
              iter        = find(logSumm.pos{trialidx(iTrial)}(:,2) >= taxis(iT),1,'first'): ...
                            find(logSumm.pos{trialidx(iTrial)}(:,2) < taxis(iT+1),1,'last');
              displ       = mean(logSumm.displ{trialidx(iTrial)}(iter,3)); 
              thisvec(iT) = displ / logSumm.frameRateVirmen;
            end
            
          case '\theta'
            thisvec   = zeros(numel(idx),1);
            for iT = 1:numel(thisvec)
              iter        = find(logSumm.pos{trialidx(iTrial)}(:,2) >= taxis(iT),1,'first'): ...
                            find(logSumm.pos{trialidx(iTrial)}(:,2) < taxis(iT+1),1,'last');
              thisvec(iT) = mean(logSumm.pos{trialidx(iTrial)}(iter,3));
            end
          
          case '\Delta_bins'
            thisvec   = zeros(numel(idx),numel(cfg.evidenceBins)-1);
            for iT = 1:size(thisvec,1)
              thisev             = sum(logSumm.cuePos_R{trialidx(iTrial)} <= taxis(iT)) - ...
                                   sum(logSumm.cuePos_L{trialidx(iTrial)} <= taxis(iT));
              binidx             = find(cfg.evidenceBins <= thisev,1,'last');       
              thisvec(iT,binidx) = 1;            
            end
            
          case '\Delta'
            thisvec   = zeros(numel(idx),1);
            for iT = 1:numel(thisvec)
              thisvec(iT) = sum(logSumm.cuePos_R{trialidx(iTrial)} <= taxis(iT)) - ...
                            sum(logSumm.cuePos_L{trialidx(iTrial)} <= taxis(iT));
            end
          
          case 'abs(\Delta)'
            thisvec   = zeros(numel(idx),1);
            for iT = 1:numel(thisvec)
              thisvec(iT) = abs(sum(logSumm.cuePos_R{trialidx(iTrial)} <= taxis(iT)) - ...
                                sum(logSumm.cuePos_L{trialidx(iTrial)} <= taxis(iT)));
            end
            
          case 'y'
            thisvec   = zeros(numel(idx),1);
            for iT = 1:numel(thisvec)
              iter        = find(logSumm.pos{trialidx(iTrial)}(:,2) >= taxis(iT),1,'first'): ...
                            find(logSumm.pos{trialidx(iTrial)}(:,2) < taxis(iT+1),1,'last');
              thisvec(iT) = mean(logSumm.pos{trialidx(iTrial)}(iter,2));
            end
            
          case 'ch'
            isLeft  = logSumm.choice(trialidx(iTrial)) == analysisParams.leftCode;
            thisvec = ones(numel(idx),1) - 2.*ones(numel(idx),1).*isLeft;
            
          case 'prevch'
            if trialidx(iTrial) == 1
              thisvec = zeros(numel(idx),1);
            else
              isLeft  = logSumm.choice(trialidx(iTrial)-1) == analysisParams.leftCode;
              thisvec = ones(numel(idx),1) - 2.*ones(numel(idx),1).*isLeft;
            end
            
          case 'prevrw'
            if trialidx(iTrial) == 1
              thisvec = zeros(numel(idx),1);
            else
              wasRW   = logSumm.choice(trialidx(iTrial)-1) == logSumm.trialType(trialidx(iTrial)-1);
              thisvec = ones(numel(idx),1) - 2.*ones(numel(idx),1).*wasRW;
            end
            
          case 'rw'
            error('rw not supported for regression in space')
            
          case 'ROI'
            continue
        end
        
        % add lags
        nlags     = cfg.predLagCm{iPred} / mode(diff(taxis)); 
        lagct     = 1;
        ncols     = size(thisvec,2);
        vecs      = zeros(size(thisvec,1),(sum(nlags)+1)*ncols);
        for iLag = -nlags(1):nlags(2)
          if iLag < 0
            shiftvec = [thisvec(-iLag+1:end,:); zeros(abs(iLag),ncols)];
          elseif iLag > 0
            shiftvec = [zeros(abs(iLag),ncols); thisvec(1:end-iLag,:)];
          else
            shiftvec = thisvec;
          end
          
          vecs(:,(lagct-1)*ncols+1:lagct*ncols) = shiftvec;
          if iTrial == 1
            predLbls{end+1} = [cfg.predList{iPred} num2str(iLag)];
          end
          lagct = lagct+1;
        end
        
        thisX(:,end+1:end+size(vecs,2)) = vecs;
      end
      
      X(end+1:end+size(thisX,1),:) = thisX;
      
    end
end

%% add ROIs separately if necessary, in which case X is 3d (one per ROI)
isROIpred = strcmpi(cfg.predList,'ROI');
if sum(isROIpred) > 1
  baseX = X;
  X     = zeros(size(X,1),size(X,2),size(y,2));
  for iROI = 1:size(y,2)
    thisX = baseX;
    ct    = 1;
    for iOther = setdiff(1:size(y,2),iROI)
      thisvec   = y(:,iOther);
      % add lags
      nlags     = cfg.predLagSec{isROIpred} / mode(diff(taxis)); 
      lagct     = 1;
      vecs      = zeros(size(thisvec,1),sum(nlags)+1);
      for iLag = -nlags(1):nlags(2)
        if iLag < 0
          shiftvec = [thisvec(-iLag+1:end); zeros(abs(iLag),1)];
        elseif iLag > 0
          shiftvec = [zeros(abs(iLag),1); thisvec(1:end-iLag)];
        else
          shiftvec = thisvec;
        end

        vecs(:,lagct) = shiftvec;
        if iROI == 1
          predLbls{end+1} = sprintf('ROI%d-%d',ct,iLag);
        end
        lagct = lagct+1;
      end
      ct = ct +1;
    end
    thisX(:,end+1:end+size(vecs,2)) = vecs;
    X(:,:,iROI) = thisX;
  end
else
  X = repmat(X,[1 1 size(y,2)]);
end

end

%% ------------------------------------------------------------------------
%% SET UP XVAL
function pseudoExp = generateXval(trialID,cfg)

rng(723176)
trials    = unique(trialID);
nt        = numel(trials);
pseudoExp = struct([]);
for iRun = 1:cfg.numCVSamples
  pseudoExp(iRun).NumTestSets = cfg.numCVFolds;
  iTrials                     = trials(randperm(nt));
  
  for iFold = 1:cfg.numCVFolds
    testidx     = floor((iFold-1)*(1/cfg.numCVFolds))*nt+1:min([nt round(nt*iFold*(1/cfg.numCVFolds))]);
    trainTrials = iTrials(setdiff(1:nt,testidx));
    testTrials  = iTrials(testidx);
    pseudoExp(iRun).training(iFold).idx = ismember(trialID,trainTrials);
    pseudoExp(iRun).test(iFold).idx     = ismember(trialID,testTrials);
  end
end

end

%% ------------------------------------------------------------------------
%% CONFIGURATION
function cfg = populateCfg(cfg)
if ~isfield(cfg,'trialType')
  cfg(1).trialType     = 'all'; 
end
if ~isfield(cfg,'predList')
  cfg(1).predList     = {'tow_L','tow_R','\Delta','y','\theta', ...
                         'd\theta/dt','speed','ch','prevch','prevrw'}; 
end
if ~isfield(cfg,'predLagSec')
  cfg(1).predLagSec    = {[0 2],[0 2],[0 1],[0 0],[0 0],[1 1],[1 1],[0 0],[0 0],[0 0]}; 
end
if ~isfield(cfg,'predLagCm')
  cfg(1).predLagCm     = {[0 100],[0 100],[0 50],[0 0],[0 0],[30 30],[30 30],[0 0],[0 0],[0 0]}; 
end
if ~isfield(cfg,'alignPoint')
  cfg(1).alignPoint    = 'cueStart';
end
if ~isfield(cfg,'timeOrSpace')
  cfg(1).timeOrSpace   = 'space';
end
if ~isfield(cfg,'posBins')
  cfg(1).posBins       = 0:1:300; % in cm
end
if ~isfield(cfg,'timeBins')
  cfg(1).timeBins      = 0:0.1:60; % in sec 
end
if ~isfield(cfg,'evidenceBins')
  cfg(1).evidenceBins      = -14:4:14; % in delta towers
  cfg(1).evidenceBins(1)   = -20;
  cfg(1).evidenceBins(end) = 20;
end
if ~isfield(cfg,'zscoreFlag')
  cfg(1).zscoreFlag    = true; 
end
if ~isfield(cfg,'diffFlag')
  cfg(1).diffFlag      = false; 
end
if ~isfield(cfg,'plotFlag')
  cfg(1).plotFlag      = true; 
end
if ~isfield(cfg,'method')
  cfg(1).method        = 'ridge'; % ridge or lasso (if the latter setting l2penalty to true will do elastic net)
end
if ~isfield(cfg,'lambdaRatio')
  switch cfg.method
    case 'lasso'
      cfg(1).lambdaRatio   = flip(10.^(linspace(-3, log10(0.5), 15))); % l1/l2 penalty
    case 'ridge'
      cfg(1).lambdaRatio   = flip(10.^(linspace(-3, log10(1), 20))); % l1/l2 penalty
  end
end
if ~isfield(cfg,'stateJitter')
  cfg.stateJitter      = [1e-7,1e-6]; % to prevent degeneracies in population states e.g. induced by zero cue-locked amplitudes
end
if ~isfield(cfg,'l2penalty')
  cfg(1).l2Penalty     = 0; % Elastic net vs. lasso
end
if ~isfield(cfg,'numCVFolds')
  cfg(1).numCVFolds    = 3;
end
if ~isfield(cfg,'numCVSamples')
  cfg(1).numCVSamples  = 70;
end
if ~isfield(cfg,'numCVSamplesShuffle')
  cfg(1).numCVSamplesShuffle = 5;
end
if ~isfield(cfg,'numShuffles')
  cfg(1).numShuffles   = 10;
end
if ~isfield(cfg,'confidLevel')
  cfg(1).confidLevel   = normcdf([-1 1], 0, 1);
end
if ~isfield(cfg,'spatialBinFactor')
  cfg(1).spatialBinFactor = 2; % spatial binning
end
end