function decoder = dffPxlDecoder(rec,spockFlag,decodeWhat,method,correctViewAng,derivFlag,cfg,...
                                 ROIflag,ROIlbl,visGuideFlag,parsedROIFlag,spatialBinFactor)

% decoder = dffPxlDecoder(rec,spockFlag,decodeWhat,method,correctViewAng,derivFlag,cfg,ROIflag,ROIlbl)
% linear decoder of task variables based on simulatenous activity of pixels
% or ROIs
%
% input:
%   rec: string with path for recording session (default: pwd)
%   spockFlag: true for running on spock (default: true)
%   decodeWhat: 'choice' (default), 'evidence', 'viewangle', 'absevidence',
%               'prevchoice', or 'decisionPoint' (not recommended)
%   method: 'ridge' (default) or 'lasso'
%   correctViewAng: true to decode on residuals of regression against view
%                   angle (default: false)
%   derivFlag: true to take derivative of dff before decoding, as a dirty form
%              of "deconvolution" (not recommended)(default: false)
%   cfg: structure with analysis config, refer to bottom of function. If
%        empty will be filled with defaults (recommended)
%   ROIflag: true to use ROIs, false to use pixels (downsampled to a 32 x
%            32 image by default)
%   ROIlbl: for ROIs, cell array with names (optional, if empty or not an
%           input labels will be loaded from disk)
%   visGuideFlag: true to decode in control task instead of towers
%   parsedROIFlag: true to use subdivided ROIs
%   spatialBinFactor: how much to spatially bin the dff data by (default 4)
%
% output is a data structure with analysis results


try
%% defaults, paths etc
if nargin < 1  || isempty(rec);            rec        = pwd;         end
if nargin < 2  || isempty(spockFlag);      spockFlag  = true;        end
if nargin < 3  || isempty(decodeWhat);     decodeWhat = 'choice';    end
if nargin < 4  || isempty(method);         method     = 'ridge';     end
if nargin < 5  || isempty(correctViewAng); correctViewAng = false;   end
if nargin < 6  || isempty(derivFlag);      derivFlag  = false;       end
if nargin < 7  || isempty(cfg);            cfg        = struct([]);  end
if nargin < 8  || isempty(ROIflag);        ROIflag    = false;       end
if nargin < 9  || isempty(ROIlbl);         ROIlbl     = {};          end
if nargin < 10 || isempty(visGuideFlag);   visGuideFlag = false;     end
if nargin < 11 || isempty(parsedROIFlag);  parsedROIFlag = false;    end
if nargin < 12 || isempty(spatialBinFactor); spatialBinFactor = 4;   end

% remove tankmousevr repo to avoid function name conflict  
if spockFlag
  warning('off','all'); 
  rmpath(genpath('/usr/people/lpinto/code/tankmousevr/')); 
  warning('on','all'); 
end  

cfg(1).decodeWhat       = decodeWhat;
cfg(1).diffFlag         = derivFlag;
cfg(1).viewAngResiduals = correctViewAng;
cfg(1).spatialBinFactor = spatialBinFactor;

switch method
  case 'ridge'
    cfg(1).method    = method;
  case 'lasso'
    cfg(1).method    = method;
    cfg(1).l2Penalty = 0;
  case 'elasticNet'
    cfg(1).method    = 'lasso';
    cfg(1).l2Penalty = 1;
end

cfg               = populateCfg(cfg);
wf                = widefieldParams;

if ~isFullPath(rec)
  rootdir         = wf.getRootDir(spockFlag);
  rec             = formatFilePath([rootdir rec]);
end

fprintf('TRAINING %s DECODER IN %s\n%s, %s trials\n', ...
        upper(cfg.decodeWhat),upper(cfg.timeOrSpace),rec,cfg.trialType);
fn = sprintf('decoder_%s_%s_%sTrials_%s',cfg.decodeWhat,method,cfg.trialType,cfg.timeOrSpace);

if ROIflag && ~parsedROIFlag
  fn = [fn '_ROI']; 
elseif ROIflag && parsedROIFlag
  fn = [fn '_ROI_parsed']; 
end

if visGuideFlag
  mz = 4;
  fn = [fn '_visGuide'];
else
  mz = [];
end
cfg.visGuideFlag = visGuideFlag;

tic
fprintf('\tloading data...')

cd(rec)
if ROIflag && ~parsedROIFlag
  if ~isempty(dir('dffROI.mat'))
    load dffROI dffROI
    dff = dffROI; clear dffROI
    load behavLog logSumm
  elseif ~isempty(dir('concatSessionsROI.mat'))
    load concatSessionsROI dffROI logSumm
    save behavLog logSumm
    dff = dffROI; clear dffROI
  else
    error('no appropriate file to read dff from')
  end
elseif ROIflag && parsedROIFlag
  if ~isempty(dir('dffROI_parsed.mat'))
    load dffROI_parsed dffROI
    dff = dffROI; clear dffROI
    isbadROI = sum(isnan(dff)) > 0;
    dff(:,isbadROI) = [];
    load behavLog logSumm
  elseif ~isempty(dir('concatSessionsROI_parsed.mat'))
    load concatSessionsROI_parsed dffROI logSumm
    save behavLog logSumm
    dff = dffROI; clear dffROI
  else
    error('no appropriate file to read dff from')
  end
else
  if ~isempty(dir('dff.mat'))
    load dff dff
    load behavLog logSumm
  elseif ~isempty(dir('concatSessions.mat'))
    load concatSessions dff logSumm
    save behavLog logSumm
  else
    error('no appropriate file to read dff from')
  end
end

fprintf(' done after %1.1f min\n',toc/60)

%% condition matrix
tic
fprintf('\tgenerating matrices...')

if ~ROIflag
  if cfg.spatialBinFactor > 1
    dff       = imBinSpace(dff,[],cfg.spatialBinFactor,true);
    fn        = sprintf('%s_binned%dx',fn,cfg.spatialBinFactor);
  end
  
  [nX,nY,~]   = size(dff);
  [dff,bc,br] = conditionDffMat(dff);
end
if cfg.diffFlag % fast dirty deconvolution by derivative
  dff = [zeros(1,size(dff,2)); diff(dff,1,1)];
  fn  = [fn '_deriv']; 
end

dff   = zscore(dff);

%% align trials (with parallelization, takes forever)
% poolobj = gcp('nocreate');
% if isempty(poolobj) 
  poolobj = parpool; 
  addAttachedFiles(poolobj, 'behavLog.mat')
% end

dffTrialsTemp = cell(1,size(dff,2));

tt = cfg.trialType;
ap = cfg.alignPoint;
tb = cfg.timeBins;
pb = cfg.posBins;
ts = cfg.timeOrSpace;
switch cfg.timeOrSpace
  case 'time'
    [dffTrials,~,~,trialidx] = alignTrials(dff(:,1),logSumm,cfg.trialType,mz,cfg.alignPoint,cfg.timeBins,cfg.timeOrSpace);
    parfor iPxl = 2:size(dff,2)
      dffTrialsTemp{iPxl}    = alignTrials(dff(:,iPxl),[],tt,mz,ap,tb,ts);
    end
  case 'space'
    [dffTrials,~,~,trialidx] = alignTrials(dff(:,1),logSumm,cfg.trialType,mz,[],cfg.posBins,cfg.timeOrSpace);
    parfor iPxl = 2:size(dff,2)
      dffTrialsTemp{iPxl}    = alignTrials(dff(:,iPxl),[],tt,mz,[],pb,ts);
    end
end

% put back in matrix format
for iPxl = 2:size(dff,2)
 dffTrials = cat(3,dffTrials,dffTrialsTemp{iPxl});
end
clear dffTrialsTemp

%% set up response vector
switch cfg.decodeWhat
  case 'choice'
    y = logSumm.choice(trialidx)';
    y(y==0) = -1;
    
  case 'prevchoice'
    tix = trialidx(trialidx > 1);
    y   = logSumm.choice(tix-1)';
    y(y==0) = -1;
    
  case 'evidence'
    y = buildEvidenceY(logSumm,trialidx,cfg);
  
  case 'absevidence'
    y = abs(buildEvidenceY(logSumm,trialidx,cfg));
    
  case 'viewangle'
    y = buildVAngY(logSumm,trialidx,cfg);  
  
  case 'decisionPoint'
    cfg.decPt.posBins   = 0:300;
    cfg.decPt.basel     = 75;
    cfg.decPt.nSD       = 1.96;
    cfg.decPt.nConsec   = 4;
    cfg.decPt.trialType = cfg.trialType;
    cfg.decPt.secDeriv  = false;
    y = estimateDecisionPoint(logSumm,cfg.decPt);  
end

%% if desired, regress view angle away first
% dff is residuals from model that tries to explain dff from view angle
if cfg.viewAngResiduals
  fn        = [fn '_viewAngResid'];
  viewAng   = buildVAngY(logSumm,trialidx,cfg);
  dffTrials = dffViewAngResiduals(dffTrials,viewAng);
end

%% set up predictor matrix 
[nTrials,nTime,nROI]   = size(dffTrials);
if cfg.singleDecoder
  % if single decoder predictor matrix for desired time point (trials x pixels)
  switch cfg.timeOrSpace
    case 'time'
      ptidx            = find(cfg.posBins  <= cfg.singleDecoderTimePoint, 1, 'last');
    case 'space'
      ptidx            = find(cfg.timeBins <= cfg.singleDecoderTimePoint, 1, 'last');
  end
  predMat              = squeeze(dffTrials(:,ptidx,:));
  
else
  % if one decoder per time point, each slice of a 3d matrix is a predictor
  % matrix for that time point (trials x pixels)
  predMat              = nan(nTrials,nROI,nTime);
  for iTime = 1:nTime
    predMat(:,:,iTime) = squeeze(dffTrials(:,iTime,:));
  end
end

% clean up nans
nanidx              = isnan(y(:,1));
y(nanidx,:)         = [];
predMat(nanidx,:,:) = [];

% compile info
decoder.nPoints     = nTime;
decoder.nTrials     = nTrials;
decoder.nPxls       = nROI;
decoder.cfg         = cfg;

fprintf(' done after %1.1f min\n',toc/60)

%% set up xval
tic; fprintf('\tfitting and cross-validating models ...');

warning('off', 'stats:cvpartition:KFoldMissingGrp');
randGenerator         = RandStream('mt19937ar', 'Seed', 723176);
pseudoExp             = cvpartition(y(:,end), 'KFold', cfg.numCVFolds, randGenerator);
for iMC = 2:cfg.numCVSamples
  pseudoExp(iMC)      = repartition(pseudoExp(1));
end
warning('on', 'stats:cvpartition:KFoldMissingGrp');

%% train / test decoder

jitter                = cfg.stateJitter;
cl                    = cfg.confidLevel;
l2                    = cfg.l2Penalty;
lambda                = cfg.lambdaRatio;
method                = cfg.method;

bestCoeff             = cell(1,iTime);
accuracy              = cell(1,iTime);
accInterval           = cell(1,iTime);
coeff                 = cell(1,iTime);
fitInfo               = cell(1,iTime);

tic;

if ~cfg.singleDecoder
  parfor iTime = 1:nTime  %#ok<PFUIX>
    if size(y,2) > 1
      target          = y(:,iTime);
    else
      target          = y;
    end
    
    % add jitter to all-zero pixels 
    pred     = predMat(:,:,iTime);
    allzeros = sum(pred==0) == size(pred,1);
    pred(:,allzeros) = pred(:,allzeros) + jitter(1) * randn(size(pred(:,allzeros)));
    
    switch method
      case 'lasso'
        [coeff{iTime}, fitInfo{iTime}] ...
                       = elasticNetRegression( pred, target, [], ...
                                              {'LeastR',jitter}, l2, lambda, pseudoExp, [], [] ); 
        fitInfo{iTime} = rmfield(fitInfo{iTime},{'X','y','w'});
        
      case 'ridge'
        [coeff{iTime}, fitInfo{iTime}] ...
                       = ridgeRegression( pred, target, lambda, pseudoExp, [], jitter); 
    end
    
    bestCoeff{iTime}  = coeff{iTime}(:,fitInfo{iTime}.IndexMinMSE);
    
    % Compute accuracy as the fraction of correct predictions, with spreads across CV samples
    predAccuracy      = nan(size(fitInfo{iTime}.CVTestPrediction));

    % Cross-validated prediction accuracy
    for iMC = 1:size(fitInfo{iTime}.CVTestPrediction,1)
      for iFold = 1:size(fitInfo{iTime}.CVTestPrediction,2)
        % Compute optimal threshold for separating training set populations according to truth
        trainSel      = fitInfo{iTime}.CVExperiments(iMC).training(iFold);

        if strcmpi(decodeWhat,'evidence')  || strcmpi(decodeWhat,'absevidence') || ...
           strcmpi(decodeWhat,'viewangle') || strcmpi(decodeWhat,'decisionPoint')
          % Compute accuracy by computing correlation between prediction
          % and data
          testSel     = fitInfo{iTime}.CVExperiments(iMC).test(iFold);
          prediction  = fitInfo{iTime}.CVTestPrediction{iMC,iFold}(:,fitInfo{iTime}.Index1SE);
          predTruth   = target(testSel);
          predAccuracy(iMC,iFold)   = corr(prediction,predTruth);
          
        else
          prediction  = fitInfo{iTime}.CVTrainPrediction{iMC,iFold}(:,fitInfo{iTime}.Index1SE);
          try
            threshold   = optimalClassThreshold(prediction, target(trainSel)>0);
          catch
            threshold   = getThreshold(prediction,0);
          end
          % Compute accuracy by applying threshold to prediction in test set and comparing to truth
          testSel     = fitInfo{iTime}.CVExperiments(iMC).test(iFold);
          prediction  = fitInfo{iTime}.CVTestPrediction{iMC,iFold}(:,fitInfo{iTime}.Index1SE);
          predTruth   = target(testSel);
          prediction(prediction < threshold) = -1;
          prediction(prediction > threshold) =  1;
          predAccuracy(iMC,iFold) = sum(prediction == predTruth)./numel(predTruth); 

        end
      end
    end

    accuracy{iTime}    = nanmean(predAccuracy(:));
    accInterval{iTime} = quantile(predAccuracy(:), cl); 

  end
else
  error('single decoder not yet implemented')
end

fprintf(' done after %1.1f min\n',toc/60)

%% Consolidate values
decoder.allCoeffs  = coeff;
decoder.fitInfo    = fitInfo;
decoder.accuracy   = cell2mat(accuracy);
decoder.weights    = cell2mat(bestCoeff)';

for iTime = 1:nTime
  decoder.accuracyCI(iTime,:) = accInterval{iTime};
end

% weights in brain image x time points (3d matrix)
if ~ROIflag
  decoder.weights  = conditionDffMat(decoder.weights,bc,br,[nX nY nTime]);
end
clear bestCoeff accuracy coeff 

%% do shuffling if necessary
if cfg.numShuffles > 1
  
  tic;
  fprintf('\tfitting and cross-validating models (shuffles)...');
  bestCoeff             = cell(cfg.numShuffles,iTime);
  accuracy              = cell(cfg.numShuffles,iTime);
  numShuff              = cfg.numShuffles;
  numCVFolds            = cfg.numCVFolds;
  numCVSamples          = 5; % in each shuffle just do 5 runs of xval
  
  % use best lambda from actual data
  bestLambdas           = zeros(1,iTime);
  for iTime = 1:nTime
    bestLambdas(iTime)  = lambda([fitInfo{iTime}.IndexMinMSE]);
  end

  parfor iTime = 1:nTime
    if size(y,2) > 1
      target          = y(:,iTime);
    else
      target          = y;
    end
    
    % add jitter to all-zero pixels 
    pred     = predMat(:,:,iTime);
    allzeros = sum(pred==0) == size(pred,1);
    pred(:,allzeros) = pred(:,allzeros) + jitter(1) * randn(size(pred(:,allzeros))); %#ok<PFBNS>
    
    for iShuff = 1:numShuff 
      
      iTarget         = target(randperm(numel(target)));
      
      switch method 
        case 'lasso'
          [bestCoeff{iShuff,iTime}, sfitInfo] = elasticNetRegression( pred, iTarget, [],                         ...
                                                                      {'LeastR',jitter}, l2, bestLambdas(iTime), ...
                                                                      numCVFolds, numCVSamples, [] );
          
        case 'ridge'
          [bestCoeff{iShuff,iTime}, sfitInfo] = ridgeRegression( pred, iTarget, bestLambdas(iTime), ...
                                                                 numCVFolds, numCVSamples, jitter);
          
      end

      % Compute accuracy as the fraction of correct predictions, with spreads across CV samples
      predAccuracy      = nan(size(sfitInfo.CVTestPrediction));
      
      % Cross-validated prediction accuracy
      for iMC = 1:size(sfitInfo.CVTestPrediction,1)
        for iFold = 1:size(sfitInfo.CVTestPrediction,2)
          % Compute optimal threshold for separating training set populations according to truth
          trainSel      = sfitInfo.CVExperiments(iMC).training(iFold);
          
          if strcmpi(decodeWhat,'evidence')  || strcmpi(decodeWhat,'absevidence') || ...
             strcmpi(decodeWhat,'viewangle') || strcmpi(decodeWhat,'decisionPoint')
            % Compute accuracy by computing correlation between prediction
            % and data
            testSel     = sfitInfo.CVExperiments(iMC).test(iFold);
            prediction  = sfitInfo.CVTestPrediction{iMC,iFold};
            predTruth   = target(testSel);
            predAccuracy(iMC,iFold)   = corr(prediction,predTruth);
            
          else
            prediction  = sfitInfo.CVTrainPrediction{iMC,iFold};
            try
              threshold   = optimalClassThreshold(prediction, target(trainSel)>0);
            catch
              threshold   = getThreshold(prediction,0);
            end
            % Compute accuracy by applying threshold to prediction in test set and comparing to truth
            testSel     = sfitInfo.CVExperiments(iMC).test(iFold);
            prediction  = sfitInfo.CVTestPrediction{iMC,iFold};
            predTruth   = target(testSel);
            prediction(prediction < threshold) = -1;
            prediction(prediction > threshold) =  1;
            predAccuracy(iMC,iFold) = sum(prediction == predTruth)./numel(predTruth); 

          end
        end
      end
      
      accuracy{iShuff,iTime}    = nanmean(predAccuracy(:));
    end
  end
  
  decoder.shuffle.coeffs        = cell2mat(bestCoeff);
  decoder.shuffle.accuracy      = cell2mat(accuracy);
  
  fprintf(' done after %1.1f min\n',toc/60)
end

%% shut down parallel pool
delete(poolobj);

%% save 
save(fn,'decoder','cfg','-v7.3')

%% plot decoding accuracy
if ~isempty(dir([fn '.pdf'])); delete([fn '.pdf']); end

switch cfg.timeOrSpace
  case 'time'
    taxis = cfg.timeBins;
    xlbl  = 'Time (s)';
  case 'space'
    taxis = cfg.posBins;
    xlbl  = 'y position (cm)';
end

if ~ROIflag
  figure; wf.applyFigDefaults(gcf,[2 1],'w'); hold on
  plot(taxis, decoder.accuracy, '-', 'color', 'k', 'linewidth', 1.5)
  plot(taxis, decoder.accuracyCI(:,1), '--', 'color', 'k',   'linewidth', .75)
  plot(taxis, decoder.accuracyCI(:,2), '--', 'color', 'k',   'linewidth', .75)

  if isfield(decoder,'shuffle')
    m   = nanmean(decoder.shuffle.accuracy);
    sem = nanstd(decoder.shuffle.accuracy)./sqrt(cfg.numShuffles);

    plot(taxis, m, '-', 'color', wf.lightgray, 'linewidth', 1)
    plot(taxis, m+sem, '--', 'color', wf.lightgray,   'linewidth', .5)
    plot(taxis, m-sem, '--', 'color', wf.lightgray,   'linewidth', .5)
  end

  wf.applyAxisDefaults(gca,'k');
  wf.applyAxisLbls(gca,xlbl,'Decoding Acc.',sprintf('%s - %s trials',cfg.decodeWhat,cfg.trialType))
  saveas(gcf,[fn '_accuracy'])
  export_fig([fn '.pdf'],'-append')
  close

  %% plot decoding weights
  [nr,nc]   = subplotOrg(nTime,8);
  cmap      = colormap(red2blue);
  cmap(1,:) = [0 0 0];
  figure; wf.applyFigDefaults(gcf,[nc nr],'k'); hold on
  for iTime = 1:nTime
    subplot(nr,nc,iTime)
    thisweight = decoder.weights(:,:,iTime);
    maxc       = max(abs(thisweight(:)));
    if isnan(maxc); continue; end
    if maxc == 0; maxc = 0.01; end
    imagesc(thisweight,[-maxc maxc]); colormap(cmap)
    axis image; axis off
    title(['y = ' num2str(taxis(iTime))],'color','w','fontsize',14,'fontweight','bold')
  end
  saveas(gcf,[fn '_weights'])
  export_fig([fn '.pdf'],'-append')
  close
else
  %%
  figure; wf.applyFigDefaults(gcf,[4 2],'w'); 
  
  subplot(2,4,[1 2]); hold on
  plot(taxis, decoder.accuracy, '-', 'color', 'k', 'linewidth', 1.5)
  plot(taxis, decoder.accuracyCI(:,1), '--', 'color', 'k',   'linewidth', .75)
  plot(taxis, decoder.accuracyCI(:,2), '--', 'color', 'k',   'linewidth', .75)

  if isfield(decoder,'shuffle')
    m   = nanmean(decoder.shuffle.accuracy);
    sem = nanstd(decoder.shuffle.accuracy)./sqrt(cfg.numShuffles);

    plot(taxis, m, '-', 'color', wf.lightgray, 'linewidth', 1)
    plot(taxis, m+sem, '--', 'color', wf.lightgray,   'linewidth', .5)
    plot(taxis, m-sem, '--', 'color', wf.lightgray,   'linewidth', .5)
  end

  wf.applyAxisDefaults(gca,'k');
  wf.applyAxisLbls(gca,xlbl,'Decoding Acc.',sprintf('%s - %s trials',cfg.decodeWhat,cfg.trialType))
  
  subplot(2,4,[5 6]); hold on
  imagesc(taxis,1:nROI,decoder.weights',[-abs(max(decoder.weights(:))) abs(max(decoder.weights(:)))]); 
  axis tight; colormap red2blue
  colorbar('location','southoutside','position',[.1 .04 .1 .02])
  set(gca,'ytick',1:nROI,'yticklabel',ROIlbl)
  wf.applyAxisDefaults(gca,'k'); 
  wf.applyAxisLbls(gca,xlbl,[],'Decoding weights (a.u.)')
  
  subplot(2,4,[3 4]); hold on
  colors = jet(nROI);%feval(widefieldParams.colors,nROI);
  weightmean = mean(decoder.weights);
  for iROI = 1:nROI
    bar(iROI,weightmean(iROI),'edgecolor',colors(iROI,:),'facecolor',colors(iROI,:)); 
  end
  set(gca,'xtick',1:nROI,'xticklabel',ROIlbl)
  rotateXLabels(gca,30)
  wf.applyAxisDefaults(gca,'k'); axis tight
  wf.applyAxisLbls(gca,[],'Decoding weight (a.u.)','Average weight')
  
  subplot(2,4,[7 8]); hold on
  weightmean = decoder.weights(decoder.accuracy == max(decoder.accuracy),:);
  for iROI = 1:nROI
    bar(iROI,weightmean(iROI),'edgecolor',colors(iROI,:),'facecolor',colors(iROI,:)); 
  end
  set(gca,'xtick',1:nROI,'xticklabel',ROIlbl)
  rotateXLabels(gca,30)
  wf.applyAxisDefaults(gca,'k'); axis tight
  wf.applyAxisLbls(gca,[],'Decoding weight (a.u.)','Weight at max accuracy')
end

catch ME
  displayException(ME)
end
end

%% --------  build evidence response vector (per trial)
function y = buildEvidenceY(logSumm,trialIdx,cfg)

switch cfg.timeOrSpace
  case 'time'
    taxis = cfg.timeBins;
    
  case 'space'
    taxis = cfg.posBins;
end

%% compile cumulative evidence tally for each trial and position (time)
y = zeros(numel(trialIdx),numel(taxis));
for iTrial = 1:numel(trialIdx)
  for iTime = 1:numel(taxis)
    switch cfg.timeOrSpace
      case 'time'
        nR = sum(logSumm.cueOnset_R{trialIdx(iTrial)} <= taxis(iTime));
        nL = sum(logSumm.cueOnset_L{trialIdx(iTrial)} <= taxis(iTime));
        
      case 'space'
        nR = sum(logSumm.cuePos_R{trialIdx(iTrial)} <= taxis(iTime));
        nL = sum(logSumm.cuePos_L{trialIdx(iTrial)} <= taxis(iTime));
        
    end
    
    y(iTrial,iTime) = nR - nL;
  end
end

end

%% --------  build view angle response vector (per trial)
function y = buildVAngY(logSumm,trialIdx,cfg)

if strcmpi(cfg.timeOrSpace,'time')
    error('decoding view angle in time does not really make sense')
end

%% compile view angle for each trial and position 
y = sampleViewAngleVsY(logSumm.pos(trialIdx), cfg.posBins);
y = y';

end

%% --------  compile cfg
function cfg = populateCfg(cfg)

%% decoder config
if ~isfield(cfg,'trialType')
  cfg(1).trialType     = 'all'; % choice, prev. choice, evidence, position, view angle
end
if ~isfield(cfg,'decodeWhat')
  cfg(1).decodeWhat    = 'choice'; % choice, prev. choice, evidence, position, view angle
end
if ~isfield(cfg,'singleDecoder')
  cfg(1).singleDecoder = false; % true to train one decoder for all time points
end
if ~isfield(cfg,'singleDecoderTimePoint')
  cfg(1).singleDecoderTimePoint = 300; % if single decoder, which point to train it on?
end
if ~isfield(cfg,'singleCut')
  cfg(1).singleCut     = false; % if false will use different categorization boundaries per time point
end
if ~isfield(cfg,'alignPoint')
  cfg(1).alignPoint    = 'cueStart';
end
if ~isfield(cfg,'timeOrSpace')
  cfg(1).timeOrSpace   = 'space';
end
if ~isfield(cfg,'posBins')
  cfg(1).posBins       = 0:10:300; % in cm
end
if ~isfield(cfg,'timeBins')
  cfg(1).timeBins      = 0:0.2:10; % in sec 
end
if ~isfield(cfg,'evidenceBins')
  cfg(1).evidenceBins  = -12:12; % #R - #L towers
end
if ~isfield(cfg,'spatialBinFactor')
  cfg(1).spatialBinFactor = 4; % spatial binning
end
if ~isfield(cfg,'diffFlag')
  cfg(1).diffFlag      = false; 
end
if ~isfield(cfg,'viewAngResiduals')
  cfg(1).viewAngResiduals  = true; 
end

%% optimization config
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
  cfg.stateJitter      = [1e-6,1e-5]; % to prevent degeneracies in population states e.g. induced by zero cue-locked amplitudes
end
if ~isfield(cfg,'l2penalty')
  cfg(1).l2Penalty     = 0; % Elastic net vs. lasso
end
if ~isfield(cfg,'xvalNFold')
  cfg(1).numCVFolds    = 3;
end
if ~isfield(cfg,'xvalNSamples')
  cfg(1).numCVSamples  = 50;
end
if ~isfield(cfg,'numShuffles')
  cfg(1).numShuffles   = 10;
end
if ~isfield(cfg,'confidLevel')
  cfg(1).confidLevel   = normcdf([-1 1], 0, 1);
end

end

%%
function threshold = getThreshold(prediction,cut)

threshold = min(prediction) + (max(prediction)-min(prediction))/2;

if cut == 0
  if abs(threshold - cut) > .1;     threshold = cut; end
else
  if abs(threshold - cut) > .1*cut; threshold = cut; end
end

end

%% 
function dffTrials = dffViewAngResiduals(dffTrials,viewAng)

viewAng                = viewAng(:);
[nTrials,nTime,nROI]   = size(dffTrials);

parfor iPxl = 1:nROI
  dff           = dffTrials(:,:,iPxl);
  dff           = dff(:);
  b             = robustfit(viewAng, dff);
  yhat          = b(1) + viewAng*b(2);
  dffTemp{iPxl} = reshape(dff - yhat,[nTrials nTime]);
end

% put back in matrix format
dffTrials = zeros(size(dffTrials));
for iPxl = 1:nROI
 dffTrials(:,:,iPxl) = dffTemp{iPxl};
end
end