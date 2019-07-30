function ROIcorr = dffCorr(dff,logSumm,ROIlbl,cfg,subtractFlag)

% ROIcorr = dffCorr(dff,logSumm,ROIlbl,cfg,subtractFlag)
% calculates pairwise correlation between ROIs
%   dff is either frames x ROI matrix or dffROI.mat folder path
%   logSumm is flattened behavior log (default: [] to load)
%   ROIlbl is cell array with ROI names (default: {} to load)
%   cfg is optional analysis parameter structure
%   subtFlag to subtract ROI grand average dff (default false)

try
tic; fprintf('calculating correlations...\n')

%% defaults
if nargin < 1;      dff          = pwd;             end
if nargin < 2;      logSumm      = [];              end
if nargin < 3;      ROIlbl       = {};              end
if nargin < 4;      cfg          = struct([]);      end
if nargin < 5;      subtractFlag = false;           end

cfg              = populateCfg(cfg);
cfg.subtractFlag = subtractFlag; 

%% load if necessary
fn = 'ROIcorr';
if subtractFlag; fn = [fn '_subt.mat']; else; fn = [fn '.mat']; end

if ischar(dff)
  cd(dff)

  if cfg.ROIflag
    if isempty(dir('dffROI.mat'))
      [dff,ROIlbl] = extractDffFromROI(dff);
    else
      load dffROI dffROI ROIlbl
      dff = dffROI; clear dffROI
    end
  else
    load dff dff
  end
  
  if isempty(logSumm) && ~isempty(dir('behavLog.mat'))
  	load behavLog logSumm
  end
end

if isempty(logSumm)
  load info frameRate;
  fs = 1/(frameRate/2);
else
  fs = logSumm.frameDtCam;
end

if size(dff,3) > 1; dff = conditionDffMat(dff); end
if subtractFlag; dff = dff - repmat(mean(dff,2),[1 size(dff,2)]); end
ROIcorr.ROIlbl = ROIlbl;

%% some recs in blocks condensed have maze 12, this is just like maze 4 but with towers on both sides
% here they will be analyzed together for convenience
if ~isempty(logSumm)
  if contains(char(logSumm.info.protocol),'Condensed')
    logSumm.currMaze(logSumm.currMaze == 12) = 4;
  end
end
%%
nROI     = size(dff,2);
cfg      = populateCfg(cfg);
maxLag   = round(cfg.maxLagSec / fs);

switch cfg.timeOrSpace
  case 'time'
    maxLagEp = round(cfg.maxLagEpoch / fs);
  case 'space'
    maxLagEp = cfg.maxLagEpoch;
end
if cfg.zscoreFlag; dff = zscore(dff); end

%% start parallel pool
poolobj = gcp('nocreate');
if isempty(poolobj); poolobj = parpool; end
if ~isempty(logSumm); addAttachedFiles(poolobj, 'behavLog.mat'); end

%% overall correlation, x corr
% also calculates an index s.t. positive values mean ROI1 lags ROI2
ROIcorr.overall.cc                = corr(dff);
ROIcorr.overall.pair              = [];
ROIcorr.overall.xcorr             = [];
ROIcorr.overall.xcorrShuffle_avg  = [];
ROIcorr.overall.xcorrShuffle_std  = [];

posIdx   = maxLag+1:maxLag+round(cfg.maxLagSecIndex / fs);
negIdx   = maxLag-round(cfg.maxLagSecIndex / fs):maxLag-1;

for iROI1 = 1:nROI
  for iROI2 = 1:nROI
    if iROI2 > iROI1
      ROIcorr.overall.pair(end+1,:)  = [iROI1 iROI2];
      ROIcorr.overall.xcorr(end+1,:) = xcorr(dff(:,iROI1),dff(:,iROI2),maxLag,cfg.xCorrNorm);

      dff1       = dff(:,iROI1);
      dff2       = dff(:,iROI2);
      sigpos     = ROIcorr.overall.xcorr(end,posIdx);
      signeg     = ROIcorr.overall.xcorr(end,negIdx);
      parfor iShuff = 1:cfg.numShuff
        thisShuff{iShuff}          = ...
          xcorr(dff1(randperm(size(dff1,1))),dff2(randperm(size(dff2,1))),maxLag,cfg.xCorrNorm);
        
        % calculate index per iter for stats
        poslags = mean(sigpos - thisShuff{iShuff}(posIdx)');
        neglags = mean(signeg - thisShuff{iShuff}(negIdx)');
        indexShuff{iShuff} = (poslags - neglags)/(poslags+neglags);
      end
      
      ROIcorr.overall.xcorrShuffle_avg(end+1,:)  = mean(cell2mat(thisShuff)');
      ROIcorr.overall.xcorrShuffle_std(end+1,:)  = std(cell2mat(thisShuff)');
      
      poslags = mean(ROIcorr.overall.xcorr(end,posIdx) - ROIcorr.overall.xcorrShuffle_avg(end,posIdx));
      neglags = mean(ROIcorr.overall.xcorr(end,negIdx) - ROIcorr.overall.xcorrShuffle_avg(end,negIdx));
      ROIcorr.overall.directionIndex(iROI1,iROI2)      = (poslags - neglags)/(poslags+neglags);
      ROIcorr.overall.directionIndex_pval(iROI1,iROI2) = sum(sign(cell2mat(indexShuff)) ~= sign(ROIcorr.overall.directionIndex(iROI1,iROI2)))/cfg.numShuff;
      ROIcorr.overall.directionIndex(iROI2,iROI1)      = - (poslags - neglags)/(poslags+neglags);
      ROIcorr.overall.directionIndex_pval(iROI2,iROI1) = sum(sign(cell2mat(indexShuff)) ~= sign(ROIcorr.overall.directionIndex(iROI1,iROI2)))/cfg.numShuff;
    end
  end
end

% average cc per ROI excluding autocorr
cc                            = ROIcorr.overall.cc;
cc(logical(eye(size(cc))))    = nan;
ROIcorr.overall.conn_avgROIcc = nanmean(cc);

%% GLM coefficient correlations ("signal")
if ~isempty(dir(cfg.whichGLM))
  load(cfg.whichGLM,'dffFit')
  ROIcorr.glm.cc = corr(dffFit.weights');
elseif isempty(dir(cfg.whichGLM)) && ~isempty(logSumm)
  warning('GLM file not found')
end

%% epoch, maze specific, only if not spontaneous rec
if ~isempty(logSumm)  
  for iMz = 1:numel(cfg.mazes)
    
    %% all frames for that maze
    switch cfg.mazes{iMz}
      case 'visGuide'
        mid = min(logSumm.currMaze);
      case 'accumul'
       mid  = max(logSumm.currMaze);
    end

    idx    = logSumm.currMaze==mid;
    frames = logSumm.camFrameNum(idx);
    frames = cellfun(@(x)(x'),frames,'UniformOutput',false);
    fid    = unique([frames{:}]');
    
    %% overall correlation, x corr
    % also calculates an index s.t. positive values mean ROI1 lags ROI2
    ROIcorr.(cfg.mazes{iMz}).overall.cc                = corr(dff(fid,:));
    ROIcorr.(cfg.mazes{iMz}).overall.pair              = [];
    ROIcorr.(cfg.mazes{iMz}).overall.xcorr             = [];
    ROIcorr.(cfg.mazes{iMz}).overall.xcorrShuffle_avg  = [];
    ROIcorr.(cfg.mazes{iMz}).overall.xcorrShuffle_std  = [];
    
    % average cc per ROI excluding autocorr
    cc                            = ROIcorr.(cfg.mazes{iMz}).overall.cc;
    cc(logical(eye(size(cc))))    = nan;
    ROIcorr.(cfg.mazes{iMz}).overall.conn_avgROIcc = nanmean(cc);
    
    % xcorr
    posIdx   = maxLag+1:maxLag+round(cfg.maxLagSecIndex / fs);
    negIdx   = maxLag-round(cfg.maxLagSecIndex / fs):maxLag-1;

    for iROI1 = 1:nROI
      for iROI2 = 1:nROI
        if iROI2 > iROI1
          ROIcorr.(cfg.mazes{iMz}).overall.pair(end+1,:)  = [iROI1 iROI2];
          ROIcorr.(cfg.mazes{iMz}).overall.xcorr(end+1,:) = xcorr(dff(fid,iROI1),dff(fid,iROI2),maxLag,cfg.xCorrNorm);

          dff1       = dff(fid,iROI1);
          dff2       = dff(fid,iROI2);
          sigpos     = ROIcorr.(cfg.mazes{iMz}).overall.xcorr(end,posIdx);
          signeg     = ROIcorr.(cfg.mazes{iMz}).overall.xcorr(end,negIdx);
          parfor iShuff = 1:cfg.numShuff
            thisShuff{iShuff}          = ...
              xcorr(dff1(randperm(size(dff1,1))),dff2(randperm(size(dff2,1))),maxLag,cfg.xCorrNorm);

            % calculate index per iter for stats
            poslags = mean(sigpos - thisShuff{iShuff}(posIdx)');
            neglags = mean(signeg - thisShuff{iShuff}(negIdx)');
            indexShuff{iShuff} = (poslags - neglags)/(poslags+neglags);
          end

          ROIcorr.(cfg.mazes{iMz}).overall.xcorrShuffle_avg(end+1,:)  = mean(cell2mat(thisShuff)');
          ROIcorr.(cfg.mazes{iMz}).overall.xcorrShuffle_std(end+1,:)  = std(cell2mat(thisShuff)');

          poslags = mean(ROIcorr.(cfg.mazes{iMz}).overall.xcorr(end,posIdx) - ROIcorr.(cfg.mazes{iMz}).overall.xcorrShuffle_avg(end,posIdx));
          neglags = mean(ROIcorr.(cfg.mazes{iMz}).overall.xcorr(end,negIdx) - ROIcorr.(cfg.mazes{iMz}).overall.xcorrShuffle_avg(end,negIdx));
          ROIcorr.(cfg.mazes{iMz}).overall.directionIndex(iROI1,iROI2)      = (poslags - neglags)/(poslags+neglags);
          ROIcorr.(cfg.mazes{iMz}).overall.directionIndex_pval(iROI1,iROI2) = sum(sign(cell2mat(indexShuff)) ~= sign(ROIcorr.(cfg.mazes{iMz}).overall.directionIndex(iROI1,iROI2)))/cfg.numShuff;
          ROIcorr.(cfg.mazes{iMz}).overall.directionIndex(iROI2,iROI1)      = - (poslags - neglags)/(poslags+neglags);
          ROIcorr.(cfg.mazes{iMz}).overall.directionIndex_pval(iROI2,iROI1) = sum(sign(cell2mat(indexShuff)) ~= sign(ROIcorr.(cfg.mazes{iMz}).overall.directionIndex(iROI1,iROI2)))/cfg.numShuff;
        end
      end
    end

      
    for iType = 1:numel(cfg.trialType)

      %% align trials
      dffTrialsTemp = cell(1,size(dff,2));

      tt = cfg.trialType{iType};
      ap = cfg.alignPoint;
      tb = cfg.timeBins;
      pb = cfg.posBins;
      ts = cfg.timeOrSpace;
      switch cfg.timeOrSpace
        case 'time'
          dffTrials = alignTrials(dff(:,1),logSumm,cfg.trialType{iType},mid,cfg.alignPoint,cfg.timeBins,cfg.timeOrSpace);
          parfor iPxl = 2:size(dff,2)
            dffTrialsTemp{iPxl}    = alignTrials(dff(:,iPxl),[],tt,mid,ap,tb,ts);
          end
        case 'space'
          dffTrials = alignTrials(dff(:,1),logSumm,cfg.trialType{iType},mid,cfg.alignPoint,cfg.posBins,cfg.timeOrSpace);
          parfor iPxl = 2:size(dff,2)
            dffTrialsTemp{iPxl}    = alignTrials(dff(:,iPxl),[],tt,mid,ap,pb,ts);
          end
      end

      % put back in matrix format
      for iPxl = 2:size(dff,2)
       dffTrials = cat(3,dffTrials,dffTrialsTemp{iPxl});
      end
      clear dffTrialsTemp

      %% correlation by trial epoch
      
      posIdx     = maxLagEp+1:maxLagEp+round(cfg.maxLagSecIndex / fs);
      negIdx     = maxLagEp-round(cfg.maxLagSecIndex / fs):maxLagEp-1;

      for iEpoch = 1:numel(cfg.epochLbl)
        
        if isempty(dffTrials)
          ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).cc                  = nan(nROI);
          ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).pair                = [nan nan];
          ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).xcorr               = nan;
          ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).xcorrShuffle_avg    = nan;
          ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).xcorrShuffle_std    = nan;
          ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).directionIndex      = nan(nROI);
          ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).directionIndex_pval = nan(nROI);
          continue;
        end

        tidx    = cfg.posBins >= cfg.epochBins{iEpoch}(1) & cfg.posBins < cfg.epochBins{iEpoch}(2);
        thisdff = [];
        for iROI = 1:nROI; temp = dffTrials(:,tidx,iROI); thisdff(:,end+1) = temp(:); end

        ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).cc                = corr(thisdff);
        ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).pair              = [];
        ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).xcorr             = [];
        ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).xcorrShuffle_avg  = [];
        ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).xcorrShuffle_std  = [];
        
        % average cc per ROI excluding autocorr
        cc                            = ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).cc;
        cc(logical(eye(size(cc))))    = nan;
        ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).conn_avgROIcc = nanmean(cc);
        
        % xcorr
        for iROI1 = 1:nROI
          for iROI2 = 1:nROI
            if iROI2 > iROI1
              ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).pair(end+1,:)  = [iROI1 iROI2];
              ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).xcorr(end+1,:) = xcorr(thisdff(:,iROI1),thisdff(:,iROI2),maxLagEp,cfg.xCorrNorm);

              dff1       = thisdff(:,iROI1);
              dff2       = thisdff(:,iROI2);
              sigpos     = ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).xcorr(end,posIdx);
              signeg     = ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).xcorr(end,negIdx);
              parfor iShuff = 1:cfg.numShuff
                thisShuff{iShuff}          = ...
                  xcorr(dff1(randperm(size(dff1,1))),dff2(randperm(size(dff2,1))),maxLag,cfg.xCorrNorm);

                % calculate index per iter for stats
                poslags = mean(sigpos - thisShuff{iShuff}(posIdx)');
                neglags = mean(signeg - thisShuff{iShuff}(negIdx)');
                indexShuff{iShuff} = (poslags - neglags)/(poslags+neglags);
              end

              ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).xcorrShuffle_avg(end+1,:)  = mean(cell2mat(thisShuff)');
              ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).xcorrShuffle_std(end+1,:)  = std(cell2mat(thisShuff)');

              poslags = mean(ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).xcorr(end,posIdx) - ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).xcorrShuffle_avg(end,posIdx));
              neglags = mean(ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).xcorr(end,negIdx) - ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).xcorrShuffle_avg(end,negIdx));
              ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).directionIndex(iROI1,iROI2)      = (poslags - neglags)/(poslags+neglags);
              ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).directionIndex_pval(iROI1,iROI2) = sum(sign(cell2mat(indexShuff)) ~= sign(ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).directionIndex(iROI1,iROI2)))/cfg.numShuff;
              ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).directionIndex(iROI2,iROI1)      = - (poslags - neglags)/(poslags+neglags);
              ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).directionIndex_pval(iROI2,iROI1) = sum(sign(cell2mat(indexShuff)) ~= sign(ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(cfg.epochLbl{iEpoch}).directionIndex(iROI1,iROI2)))/cfg.numShuff;
            end
          end
        end
      end
    end
  end
end

%% shut down parallel pool
delete(poolobj);

%% save
ROIcorr.cfg = cfg;
fprintf('saving and plotting...\n')
save(fn,'ROIcorr','-v7.3')

%% plot correlations, cross-correlations
wf    = widefieldParams;
if isempty(logSumm)
  fh    = figure;
  plots = {'overall'};
  wf.applyFigDefaults(fh,[numel(plots)+1 2],'k')

  if ~isempty(cfg.sortBy)
    sortBy = strcmpi(ROIlbl,cfg.sortBy);
    if sum(sortBy) > 0
      sortIdx = find(sortBy);
    else
      sortIdx = 1;
    end
  end
  for iPlot = 1:numel(plots)
    thiscorr  = ROIcorr.(plots{iPlot}).cc;
    [~,order] = sort(thiscorr(:,sortIdx),'descend');

    subplot(2,numel(plots),iPlot)
    imagesc(thiscorr(order,order)); colormap red2blue
    set(gca,'xtick',1:nROI,'ytick',1:nROI,'xticklabel',ROIlbl(order),'yticklabel',ROIlbl(order))
    rotateXLabels(gca,90)
    wf.applyAxisDefaults(gca,'w'); axis tight
    wf.applyAxisLbls(gca,[],[],['corr - ' plots{iPlot}])

    if iPlot == numel(plots)
      ch = colorbar('location','eastoutside','position',[.92 .58 .01 .1],'Color','w');
      ch.Label.String = 'r';
    end

    thisindex  = ROIcorr.(plots{iPlot}).directionIndex;
    thisp      = ROIcorr.(plots{iPlot}).directionIndex_pval;
    thisindex(thisp > cfg.alpha) = 0; 
    [~,order] = sort(abs(thisindex(sortIdx,:)),'ascend');
    subplot(2,numel(plots),iPlot+numel(plots))
    imagesc(thisindex(order,order)); colormap red2blue
    set(gca,'xtick',1:nROI,'ytick',1:nROI,'xticklabel',ROIlbl(order),'yticklabel',ROIlbl(order))
    rotateXLabels(gca,90)
    wf.applyAxisDefaults(gca,'w'); axis tight
    wf.applyAxisLbls(gca,[],[],['x corr index - ' plots{iPlot}])

    if iPlot == numel(plots)
      ch = colorbar('location','eastoutside','position',[.92 .11 .01 .1],'Color','w');
      ch.Label.String = 'DI (<0 = leads)';
    end
  end
else
  %% compare visual guide and accumul
  
  % 1) overall mazes
  fh    = figure;
  plots = 'overall';
  wf.applyFigDefaults(fh,[numel(cfg.mazes)+1 2],'k')
  if ~isempty(cfg.sortBy)
    sortBy = strcmpi(ROIcorr.ROIlbl,cfg.sortBy);
    if sum(sortBy) > 0
      sortIdx = find(sortBy);
    else
      sortIdx = 1;
    end
  end
  
  for iMz = 1:numel(cfg.mazes)
    thiscorr  = ROIcorr.(cfg.mazes{iMz}).(plots).cc;
    [~,order] = sort(thiscorr(:,sortIdx),'descend');

    subplot(2,numel(cfg.mazes),iMz)
    imagesc(thiscorr(order,order),[0 1]); colormap red2blue
    set(gca,'xtick',1:nROI,'ytick',1:nROI,'xticklabel',ROIlbl(order),'yticklabel',ROIlbl(order))
    rotateXLabels(gca,90)
    wf.applyAxisDefaults(gca,'w'); axis tight
    wf.applyAxisLbls(gca,[],[],sprintf('corr %s %s',cfg.mazes{iMz},plots))

    if iMz == numel(cfg.mazes)
      ch = colorbar('location','eastoutside','position',[.92 .58 .01 .1],'Color','w');
      ch.Label.String = 'r';
    end

    thisindex  = ROIcorr.(cfg.mazes{iMz}).(plots).directionIndex;
    thisp      = ROIcorr.(cfg.mazes{iMz}).(plots).directionIndex_pval;
    thisindex(thisp > cfg.alpha) = 0; 
    [~,order] = sort(abs(thisindex(sortIdx,:)),'ascend');

    subplot(2,numel(cfg.mazes),iMz+numel(cfg.mazes))
    imagesc(thisindex(order,order),[-1.5 1.5]); colormap red2blue
    set(gca,'xtick',1:nROI,'ytick',1:nROI,'xticklabel',ROIlbl(order),'yticklabel',ROIlbl(order))
    rotateXLabels(gca,90)
    wf.applyAxisDefaults(gca,'w'); axis tight
    wf.applyAxisLbls(gca,[],[],sprintf('dir. index %s %s',cfg.mazes{iMz},plots))

    if iMz == numel(cfg.mazes)
      ch = colorbar('location','eastoutside','position',[.92 .11 .01 .1],'Color','w');
      ch.Label.String = 'DI (<0 = leads)';
    end
  end
  
  % 2) subtraction
  fh    = figure;
  wf.applyFigDefaults(fh,[2 2],'k')
  visGuide  = ROIcorr.visGuide.overall.cc;
  accumul   = ROIcorr.accumul.overall.cc;
  subt      = accumul-visGuide;
  [~,order] = sort(accumul(:,sortIdx),'descend');

  imagesc(subt(order,order),[-max(abs(subt(:))) max(abs(subt(:)))]); colormap red2blue
  set(gca,'xtick',1:nROI,'ytick',1:nROI,'xticklabel',ROIlbl(order),'yticklabel',ROIlbl(order))
  rotateXLabels(gca,90)
  wf.applyAxisDefaults(gca,'w'); axis tight
  wf.applyAxisLbls(gca,[],[],'accumul - visGuide, overall')

  ch = colorbar('location','eastoutside','position',[.92 .1 .03 .3],'Color','w');
  ch.Label.String = 'r';
  
  % 3) trial types  
  plots = cfg.epochLbl;
  for iMz = 1:numel(cfg.mazes)
    for iType = 1:numel(cfg.trialType)
      fh    = figure;
      wf.applyFigDefaults(fh,[numel(plots)+1 2],'k')

      if ~isempty(cfg.sortBy)
        sortBy = strcmpi(ROIlbl,cfg.sortBy);
        if sum(sortBy) > 0
          sortIdx = find(sortBy);
        else
          sortIdx = 1;
        end
      end
      for iPlot = 1:numel(plots)
        thiscorr  = ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(plots{iPlot}).cc;
        [~,order] = sort(thiscorr(:,sortIdx),'descend');

        subplot(2,numel(plots),iPlot)
        imagesc(thiscorr(order,order),[0 1]); colormap red2blue
        set(gca,'xtick',1:nROI,'ytick',1:nROI,'xticklabel',ROIlbl(order),'yticklabel',ROIlbl(order))
        rotateXLabels(gca,90)
        wf.applyAxisDefaults(gca,'w'); axis tight
        wf.applyAxisLbls(gca,[],[],sprintf('corr %s %s %s',cfg.mazes{iMz},cfg.trialType{iType},plots{iPlot}))

        if iPlot == numel(plots)
          ch = colorbar('location','eastoutside','position',[.92 .58 .01 .1],'Color','w');
          ch.Label.String = 'r';
        end

        thisindex  = ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(plots{iPlot}).directionIndex;
        thisp      = ROIcorr.(cfg.mazes{iMz}).(cfg.trialType{iType}).(plots{iPlot}).directionIndex_pval;
        thisindex(thisp > cfg.alpha) = 0; 
        [~,order] = sort(abs(thisindex(sortIdx,:)),'ascend');
        subplot(2,numel(plots),iPlot+numel(plots))
        imagesc(thisindex(order,order),[-1.5 1.5]); colormap red2blue
        set(gca,'xtick',1:nROI,'ytick',1:nROI,'xticklabel',ROIlbl(order),'yticklabel',ROIlbl(order))
        rotateXLabels(gca,90)
        wf.applyAxisDefaults(gca,'w'); axis tight
        wf.applyAxisLbls(gca,[],[],sprintf('dir. index %s %s %s',cfg.mazes{iMz},cfg.trialType{iType},plots{iPlot}))

        if iPlot == numel(plots)
          ch = colorbar('location','eastoutside','position',[.92 .11 .01 .1],'Color','w');
          ch.Label.String = 'DI (<0 = leads)';
        end
      end

    end
  end 
end

%% plot network metrics
if ~isempty(logSumm)
  fh    = figure;
  wf.applyFigDefaults(fh,[4 2],'w')
  
  % plot scatter plot comparisons between different mazes / epochs
  for iPlot = 1:3
    subplot(3,3,iPlot); hold on;
    switch iPlot
      case 1
        xlbl = 'visGuide overall'; ylbl = 'accumul overall';
      case 2 
        xlbl = 'accumul correct cueHalf1'; ylbl = 'accumul correct mem';
      case 3
        xlbl = 'accumul correct mem'; ylbl = 'accumul error mem';
    end
    
    vg   = degrees(strcmpi(conn_lbl,xlbl),:);
    acc  = degrees(strcmpi(conn_lbl,ylbl),:);
    
    if cfg.ROIflag 

      plotROIlbl  = {'V1','mV2','PPC','RSC','aM2','mM2','M1','SS'};
      plotROI     = {'VISp','mV2','VISa','RSP','aMOs','mMOs','MOp','SS'};
      
      for iROI = 1:numel(plotROI)
        idx = strcmpi(ROIcorr.ROIlbl,[plotROI{iROI} '-R']) | strcmpi(ROIcorr.ROIlbl,[plotROI{iROI} '-L']);
        plot(vg(idx),acc(idx),'o','color',wf.areaCl(iROI,:),'markerfacecolor',wf.areaCl(iROI,:),'markersize',8)
      end
      yl     = get(gca,'ylim');
      xl     = get(gca,'xlim');
      xlim([min([xl yl]) max([xl yl])]);
      ylim([min([xl yl]) max([xl yl])]);
      yl     = get(gca,'ylim');
      xl     = get(gca,'xlim');
      plot(xl,xl,'--','color',[.7 .7 .7])
      for iROI = 1:numel(plotROIlbl)
        text(xl(1)+1,yl(2)-.55*(iROI-1),plotROIlbl{iROI},'color',wf.areaCl(iROI,:))
      end
      
    else
      plot(vg,acc,'ko')
      yl     = get(gca,'ylim');
      xl     = get(gca,'xlim');
      xlim([min([xl yl]) max([xl yl])]);
      ylim([min([xl yl]) max([xl yl])]);
      yl     = get(gca,'ylim');
      xl     = get(gca,'xlim');
      
    end
    wf.applyAxisDefaults(gca,'k'); 
    wf.applyAxisLbls(gca,xlbl,ylbl,'connectivity (nodes)')
  end
  
  % plot averages for all comparisons
  subplot(3,3,4); hold on;
  wf.applyAxisDefaults(gca,'k'); 
  bar(mean(degrees,2),'facecolor',[1 1 1]/2)
  set(gca,'xtick',1:numel(conn_lbl),'xticklabel',conn_lbl)
  rotateXLabels(gca,90)
  wf.applyAxisLbls(gca,[],'avg. connectivity (nodes)')
  
  subplot(3,3,5); hold on;
  wf.applyAxisDefaults(gca,'k'); 
  bar(transit,'facecolor',[1 1 1]/2)
  set(gca,'xtick',1:numel(conn_lbl),'xticklabel',conn_lbl)
  rotateXLabels(gca,90)
  wf.applyAxisLbls(gca,[],'transitivity (clustering)')
  
  subplot(3,3,6); hold on;
  wf.applyAxisDefaults(gca,'k'); 
  bar(effic,'facecolor',[1 1 1]/2)
  set(gca,'xtick',1:numel(conn_lbl),'xticklabel',conn_lbl)
  rotateXLabels(gca,90)
  wf.applyAxisLbls(gca,[],'efficiency (integration)')
end

%% export to pdf
fn = 'ROIcorr';
if subtractFlag; fn = [fn '_subt.pdf']; else; fn = [fn '.pdf']; end
if ~isempty(dir(fn)); delete(fn); end
figls = get(0,'children'); % print to pdf
for ii = length(figls):-1:1
  figure(figls(ii))
  export_fig(fn, '-q101', '-append')
end
close all
  
fprintf(' done after %1.1f min\n',toc/60)
catch ME
  displayException(ME)
end
end

%% ------------------------------------------------------------------------
function cfg = populateCfg(cfg)
if ~isfield(cfg,'trialType')
  cfg(1).trialType     = {'all','correct','error','easy','hard','correct_nodistractor','correct_distractor'}; 
end
if ~isfield(cfg,'mazes')
  cfg(1).mazes = {'visGuide','accumul'};
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
if ~isfield(cfg,'zscoreFlag')
  cfg(1).zscoreFlag    = true; 
end
if ~isfield(cfg,'maxLagSec')
  cfg(1).maxLagSec     = 10; 
end
if ~isfield(cfg,'maxLagSecEpoch')
  cfg(1).maxLagEpoch   = 10; 
end
if ~isfield(cfg,'maxLagSecEpoch')
  cfg(1).maxLagSecIndex = .5; 
end
if ~isfield(cfg,'xCorrNorm')
  cfg(1).xCorrNorm = 'coeff'; 
end
if ~isfield(cfg,'numShuff')
  cfg(1).numShuff  = 2; 
end
if ~isfield(cfg,'sortBy')
  cfg(1).sortBy  = 'VISp-R'; 
end
if ~isfield(cfg,'alpha')
  cfg(1).alpha  = .05; 
end
if ~isfield(cfg,'ROIflag')
  cfg(1).ROIflag  = true; 
end
if ~isfield(cfg,'epochLbl')
  cfg(1).epochLbl  = {'cueHalf1','cueHalf2','mem','wholeTrial'};
end
if ~isfield(cfg,'epochBins')
  cfg(1).epochBins  = {[0 100],[100 200],[200 300],[0 300]};
end

%% which glm to use
if ~isfield(cfg,'whichGLM')
  cfg(1).whichGLM = 'dffGLM_space_ridge_ROI.mat';
end
%% for connectivity toolbox
if ~isfield(cfg,'fisherZ')
  cfg(1).fisherZ  = true; 
end
if ~isfield(cfg,'thCC')
  cfg(1).thCC    = 1; 
end
if ~isfield(cfg,'binarize')
  cfg(1).binarize = false;
end
end