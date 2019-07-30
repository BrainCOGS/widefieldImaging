function decodeSumm = summarizeDecoding_reduced(mice,cfg,whichDecoder,concatSessFlag)

% decodeSumm = summarizeDecoding(mice,cfg,whichDecoder)
%         mice: cell array with mouse name list
%          cfg: structure with analysis config, pass empty to load defaults    
% whichDecoder: 'evidence', 'choice', 'prevchoice' or 'viewangle'

%% initialize
if nargin < 1 || isempty(mice); mice           = {'ai2';'ai3';'ai5';'ai7';'ai9';'ai10'}; end
if nargin < 2 || isempty(cfg);  cfg            = struct([]);                end
if nargin < 3;                  whichDecoder   = [];                        end
if nargin < 4;                  concatSessFlag = false;                     end

cfg = populateCfg(cfg);
if ~isempty(whichDecoder); cfg.whichDecoder = whichDecoder; end
cfg.concatSessFlag = concatSessFlag;
warning off
try

%% which files to load?
decodeSumm.allROIs.fn               = sprintf('decoder_%s_%s_%sTrials_%s_ROI.mat',cfg.whichDecoder,                        ...
                                            cfg.whichMethod,cfg.trialType,cfg.timeOrSpace);
decodeSumm.allROIs_vAngResid.fn     = sprintf('decoder_%s_%s_%sTrials_%s_ROI_viewAngResid.mat',cfg.whichDecoder,           ...
                                          cfg.whichMethod,cfg.trialType,cfg.timeOrSpace);
decodeSumm.allPxls.fn               = sprintf('decoder_%s_%s_%sTrials_%s_binned4x.mat',cfg.whichDecoder,                      ...
                                          cfg.whichMethod,cfg.trialType,cfg.timeOrSpace);
decodeSumm.allPxls_vAngResid.fn     = sprintf('decoder_%s_%s_%sTrials_%s_binned4x_viewAngResid.mat',cfg.whichDecoder,         ...
                                          cfg.whichMethod,cfg.trialType,cfg.timeOrSpace);
decodeSumm.allPxls16x.fn            = sprintf('decoder_%s_%s_%sTrials_%s_binned16x.mat',cfg.whichDecoder,                      ...
                                          cfg.whichMethod,cfg.trialType,cfg.timeOrSpace);
decodeSumm.ROIpxls.fn               = sprintf('decoder_pxlsFromROIs_%s_%s_%sTrials_%s.mat',cfg.whichDecoder,               ...
                                          cfg.whichMethod,cfg.trialType,cfg.timeOrSpace);
decodeSumm.ROIpxls_vAngResid.fn     = sprintf('decoder_pxlsFromROIs_%s_%s_%sTrials_%s_viewAngResid.mat',cfg.whichDecoder,  ...
                                          cfg.whichMethod,cfg.trialType,cfg.timeOrSpace);

decoderList     = fieldnames(decodeSumm);

%% compile from saved data files
tic
fprintf('compiling data...\n')
recls = cell(numel(mice),1);
wf    = widefieldParams;
for iMouse = 1:numel(mice)
  
  fprintf('\tmouse %d / %d',iMouse,numel(mice))
  
  if concatSessFlag
    rootdir       = wf.getRootDir(isThisSpock);
    recls{iMouse} = {[rootdir mice{iMouse}]};
  else
    [recls{iMouse},rootdir] = getFullRecPath(mice{iMouse},isThisSpock);
  end
  for iRec = 1:numel(recls{iMouse})
    fprintf('.')
    cd(recls{iMouse}{iRec})
    if cfg.ROIflag && ~exist('ROIlbl','var')
    	load dffROI ROIlbl
      ROIlbl_unil = unique(cellfun(@(x)(x(1:end-2)),ROIlbl,'uniformOutput',false),'stable');
    end
    
    % compile accuracy of different types of decoder, weights for decoder
    % that takes all ROIs
    for iDecoder = 1:numel(decoderList) 
      load(decodeSumm.(decoderList{iDecoder}).fn,'decoder'); 
      
      switch decoderList{iDecoder}
        case {'allROIs','allROIs_vAngResid','allPxls','allPxls16x','allPxls_vAngResid'}
          decodeSumm.(decoderList{iDecoder}).mouse(iMouse).accuracy(:,:,iRec)                ...
                                        = decoder.accuracy;
          decodeSumm.(decoderList{iDecoder}).mouse(iMouse).accuracy_shuffle(:,:,iRec)        ...
                                        = nanmean(decoder.shuffle.accuracy);
          if ~isempty(strfind(decoderList{iDecoder},'ROI'))
            weights = decoder.weights;
            if cfg.zscoreWeights; weights = zscore(weights')'; end
            decodeSumm.(decoderList{iDecoder}).mouse(iMouse).weights(:,:,iRec) = weights;
          else
            weights = decoder.weights;
            decodeSumm.(decoderList{iDecoder}).mouse(iMouse).weights_im{iRec} = weights;
          end
          
          % collect weight spatial correlation
          if strcmpi(decoderList{iDecoder},'allPxls')
            decodeSumm.(decoderList{iDecoder}).pxlDistAxis = decoder.weightsSpatialCorr.distAxis;
            decodeSumm.(decoderList{iDecoder}).mouse(iMouse).weightSpatialCorr(:,iRec)         = decoder.weightsSpatialCorr.corrByDist_mean;
            decodeSumm.(decoderList{iDecoder}).mouse(iMouse).weightSpatialCorr_shuffle(:,iRec) = decoder.weightsSpatialCorr.corrByDist_shuffle_mean;
            decodeSumm.(decoderList{iDecoder}).mouse(iMouse).weightSpatialCorr_isSig(:,iRec)   = decoder.weightsSpatialCorr.corrByDist_isSig_VSshuffle;
          end
          
        case {'singROI','singROI_vAngResid','singROIbil','singROIbil_vAngResid'}
          decodeSumm.(decoderList{iDecoder}).mouse(iMouse).accuracy(:,:,iRec)                ...
                                        = decoder.accuracy';
          decodeSumm.(decoderList{iDecoder}).mouse(iMouse).accuracy_shuffle(:,:,iRec)        ...
                                        = squeeze(nanmean(decoder.shuffle.accuracy))';
        
        case {'ROIleaveOut','ROIleaveOut_vAngResid'}
          nROI = numel(decoder.ROIexclude);
          for iROI = 1:nROI
            decodeSumm.(decoderList{iDecoder}).mouse(iMouse).accuracy(iROI,:,iRec)                ...
                                          = decoder.ROIexclude{iROI}.accuracy';
            decodeSumm.(decoderList{iDecoder}).mouse(iMouse).accuracy_shuffle(iROI,:,iRec)        ...
                                          = squeeze(nanmean(decoder.ROIexclude{iROI}.shuffle.accuracy))';
          end          
        case {'ROIpxls','ROIpxls_vAngResid','allPxls_leaveOut','allPxls_leaveOut_vAngResid'}
          for iROI = 1:numel(ROIlbl_unil)
            idx = strcmpi(ROIlbl_unil{iROI},decoder.ROIlbl);
            decodeSumm.(decoderList{iDecoder}).mouse(iMouse).accuracy(iROI,:,iRec)           ...
                                          = decoder.ROI(idx).accuracy;
            decodeSumm.(decoderList{iDecoder}).mouse(iMouse).accuracy_shuffle(iROI,:,iRec)   ...
                                          = nanmean(decoder.ROI(idx).shuffle.accuracy);
          end
          
      end
      if ~exist('xaxis','var')
        switch cfg.timeOrSpace
          case 'space'
            xaxis = decoder.cfg.posBins;
          case 'time'
            xaxis = decoder.cfg.timeBins;
        end
      end   
    end
  end
  
  %% mouse average
  for iDecoder = 1:numel(decoderList) 
    decodeSumm.(decoderList{iDecoder}).accuracy(:,:,iMouse)              ...
          = nanmean(decodeSumm.(decoderList{iDecoder}).mouse(iMouse).accuracy,3);
    decodeSumm.(decoderList{iDecoder}).accuracy_sem(:,:,iMouse)          ...
          = nanstd(decodeSumm.(decoderList{iDecoder}).mouse(iMouse).accuracy,0,3)./ sqrt(iRec-1);
        
    decodeSumm.(decoderList{iDecoder}).accuracy_shuffle(:,:,iMouse)      ...
          = nanmean(decodeSumm.(decoderList{iDecoder}).mouse(iMouse).accuracy_shuffle,3);
    decodeSumm.(decoderList{iDecoder}).accuracy_shuffle_sem(:,:,iMouse)  ...
          = nanstd(decodeSumm.(decoderList{iDecoder}).mouse(iMouse).accuracy_shuffle,0,3)./ sqrt(iRec-1);
    
    if isfield(decodeSumm.(decoderList{iDecoder}).mouse(iMouse),'weights')
      decodeSumm.(decoderList{iDecoder}).weights(:,:,iMouse)             ...
            = nanmean(decodeSumm.(decoderList{iDecoder}).mouse(iMouse).weights,3);
      decodeSumm.(decoderList{iDecoder}).weights_sem(:,:,iMouse)         ...
            = nanstd(decodeSumm.(decoderList{iDecoder}).mouse(iMouse).weights,0,3)./ sqrt(iRec-1);
    end
    
    % collect weight spatial correlation
    if strcmpi(decoderList{iDecoder},'allPxls')
      decodeSumm.(decoderList{iDecoder}).weightSpatialCorr(:,iMouse)         = nanmean(decodeSumm.(decoderList{iDecoder}).mouse(iMouse).weightSpatialCorr(:,iRec),2);
      decodeSumm.(decoderList{iDecoder}).weightSpatialCorr_shuffle(:,iMouse) = nanmean(decodeSumm.(decoderList{iDecoder}).mouse(iMouse).weightSpatialCorr_shuffle,2);
    end
  end
  fprintf('\n')
  
end

%% compile more info
fprintf('\tdone after %1.1f min\n',toc/60)

if ~cfg.ROIflag; ROIlbl = {}; ROIlbl_unil = {}; end
decodeSumm.cfg         = cfg;
decodeSumm.xaxis       = xaxis;
decodeSumm.recls       = recls;
decodeSumm.mice        = mice;
decodeSumm.ROIlbl      = ROIlbl;
decodeSumm.ROIlbl_unil = ROIlbl_unil;

cd(rootdir)
save(['temp_decodeSumm_' whichDecoder])

%% significance of decoder predictions

% first build a vector including all time points from all decoder shuffles
% in order to obtain a single statistical threshold (95th percentile), then
% apply to decide significant prediction accuracies
shuff = [];
for iDecoder = 1:numel(decoderList)
  shuff = [shuff; decodeSumm.(decoderList{iDecoder}).accuracy_shuffle(:)];
end
decodeSumm.shufflePrctile95 = prctile(shuff,95);
for iDecoder = 1:numel(decoderList)
  decodeSumm.(decoderList{iDecoder}).accuracy_isSig = ...
    decodeSumm.(decoderList{iDecoder}).accuracy > decodeSumm.shufflePrctile95;
end

%% maximal / average accuracies for different decoding models
for iDecoder = 1:numel(decoderList) 
  nROI    = size(decodeSumm.(decoderList{iDecoder}).accuracy,1);
  accMax  = nan(numel(mice),nROI);
  accMaxT = nan(numel(mice),nROI);
  accMean = nan(numel(mice),nROI);
  accLat  = nan(numel(mice),nROI);
  for iMouse = 1:numel(mice)
    for iROI = 1:nROI
      thisacc               = squeeze(decodeSumm.(decoderList{iDecoder}).accuracy(iROI,:,iMouse));
      accMax(iMouse,iROI)   = nanmax(thisacc);
      if ~isnan(nanmax(thisacc))
        accMaxT(iMouse,iROI)= xaxis(find(thisacc == nanmax(thisacc), 1, 'first'));
      end
      accMean(iMouse,iROI)  = nanmean(thisacc);
      if size(thisacc,1) > 1; thisacc = thisacc'; end               
      sigIdx                = thisacc                  > decodeSumm.shufflePrctile95 & ...
                              [thisacc(2:end) nan]     > decodeSumm.shufflePrctile95 & ...
                              [thisacc(3:end) nan nan] > decodeSumm.shufflePrctile95;
      firstT                = find(sigIdx, 1, 'first');
      if ~isempty(firstT)
        accLat(iMouse,iROI) = xaxis(firstT);
      end
    end
  end
  
  decodeSumm.(decoderList{iDecoder}).accuracy_avg    = accMean;
  decodeSumm.(decoderList{iDecoder}).accuracy_max    = accMax;
  decodeSumm.(decoderList{iDecoder}).accuracy_maxy   = accMaxT;
  decodeSumm.(decoderList{iDecoder}).accuracy_firsty = accLat;
  
  % one-way repeated mesures ANOVA for prediction accuracies, times
  switch decoderList{iDecoder}
    case {'singROI','singROI_vAngResid','singROIbil','singROIbil_vAngResid',   ...
          'ROIpxls','ROIpxls_vAngResid','ROIleaveOut','ROIleaveOut_vAngResid', ...
          'allPxls_leaveOut','allPxls_leaveOut_vAngResid'}
      anova_types = {'accuracy_avg','accuracy_max','accuracy_maxy','accuracy_firsty'};
      for iType = 1:numel(anova_types)
        datavec      = decodeSumm.(decoderList{iDecoder}).(anova_types{iType});
        [nmice,nROI] = size(datavec);
        datavec      = datavec(:);
        micevec      = repmat((1:nmice)',[nROI 1]);
        roivec       = [];
        for iROI = 1:nROI
          roivec     = [roivec; ones(nmice,1)*iROI];
        end
        
        [decodeSumm.(decoderList{iDecoder}).stats.(anova_types{iType}).anovaP,     ...
         decodeSumm.(decoderList{iDecoder}).stats.(anova_types{iType}).anovaTable, ...
         decodeSumm.(decoderList{iDecoder}).stats.(anova_types{iType}).anovaStats] ...
                    = anovan(datavec(~isnan(datavec)),{roivec(~isnan(datavec)),micevec(~isnan(datavec))},'display','off');
        decodeSumm.(decoderList{iDecoder}).stats.(anova_types{iType}).multComp     ...
                     = multcompare(decodeSumm.(decoderList{iDecoder}).stats.(anova_types{iType}).anovaStats,'display','off');
      end
      
  end
end

%% pairwise tests for prediction accuracies, all ROIs vs individual ROIs
test_types   = {'accuracy_avg','accuracy_max','accuracy_maxy','accuracy_firsty'};

comps{1}     = {'ROIpxls','allPxls'};
% comps{end+1} = {'ROIleaveOut_vAngResid','allROIs_vAngResid'};
% comps{end+1} = {'allPxls_leaveOut','allPxls'};
% comps{end+1} = {'allPxls_vAngResid','allPxls_leaveOut_vAngResid'};

for iComp = 1:numel(comps)
  for iType = 1:numel(test_types)
    allROIvec    = decodeSumm.(comps{iComp}{2}).(test_types{iType});
    datavec      = decodeSumm.(comps{iComp}{1}).(test_types{iType});
    nROI         = size(datavec,2);
    pval         = zeros(nROI,1);
    for iROI = 1:nROI
      if sum(~isnan(datavec(:,iROI))) < 2
        pval(iROI) = nan;
      else
        pval(iROI) = signrank(allROIvec,datavec(:,iROI));
      end
    end
    [~,alphac] = FDR(pval(~isnan(pval)),.05); % multiple comparisons correction
    decodeSumm.(comps{iComp}{1}).stats.(test_types{iType}).vsAllROIs.test   = 'signrank';
    decodeSumm.(comps{iComp}{1}).stats.(test_types{iType}).vsAllROIs.pval   = pval';
    decodeSumm.(comps{iComp}{1}).stats.(test_types{iType}).vsAllROIs.alphac = alphac;
  end
end

%% cluster weights and absolute weights
if cfg.ROIflag
   for iDecoder = 1:numel(decoderList) 
     if ~isfield(decodeSumm.(decoderList{iDecoder}),'weights'); continue; end
     
     % include only weights for bins with significant predictions
     weights          = decodeSumm.(decoderList{iDecoder}).weights;
     isSig            = squeeze(decodeSumm.(decoderList{iDecoder}).accuracy_isSig);
     isSig            = repmat(reshape(isSig,[size(isSig,1) 1 size(isSig,2)]),[1 size(weights,2) 1]);
     weights(~isSig)  = nan;
     weights          = nanmean(weights,3);
%      weights          = nanmean(weights);
     [nanr,~]         = find(isnan(weights));
     weights(nanr,:)  = [];
     [~,nanc]         = find(isnan(weights));
     weights(:,nanc)  = [];
     cc               = corr(weights);
     [pcload,~,ev]    = pca(cc);
     decodeSumm.(decoderList{iDecoder}).clustStats.fracVarPCs  ...
                      = sum(ev(1:cfg.clustMaxPC))/sum(ev);
     [decodeSumm.(decoderList{iDecoder}).clustStats.clustID,~, ...
      decodeSumm.(decoderList{iDecoder}).clustStats.lnk,decodeSumm.(decoderList{iDecoder}).clustStats.CHindex] ...
                      = hierClustering(pcload(:,1:cfg.clustMaxPC),'eucl',cfg.clustMaxK,0);
     decodeSumm.(decoderList{iDecoder}).clustStats.nK          ...
                      = numel(unique(decodeSumm.(decoderList{iDecoder}).clustStats.clustID));
                    
     cc               = corr(abs(weights));
     [pcload,~,ev]    = pca(cc);
     decodeSumm.(decoderList{iDecoder}).clustStats_absW.fracVarPCs  ...
                      = sum(ev(1:cfg.clustMaxPC))/sum(ev);
     [decodeSumm.(decoderList{iDecoder}).clustStats_absW.clustID,~, ...
      decodeSumm.(decoderList{iDecoder}).clustStats_absW.lnk,decodeSumm.(decoderList{iDecoder}).clustStats_absW.CHindex] ...
                      = hierClustering(pcload(:,1:cfg.clustMaxPC),'eucl',cfg.clustMaxK,0);
     decodeSumm.(decoderList{iDecoder}).clustStats_absW.nK          ...
                      = numel(unique(decodeSumm.(decoderList{iDecoder}).clustStats_absW.clustID));
   end

end

%% save
cd(rootdir)
delete(['temp_decodeSumm_' whichDecoder '.mat'])
fn = ['decodingSummary_' whichDecoder]; 
if cfg.concatSessFlag; fn = [fn '_concatSess']; end
save(fn,'decodeSumm')

%% plot
  
wf          = widefieldParams;
if cfg.ROIflag
  
  %% all ROIs, accuracy and weights
  thisdecoderList = {'allROIs','allROIs_vAngResid','allPxls','allPxls_vAngResid'};
  [nr,nc] = subplotOrg(numel(thisdecoderList)*3,numel(thisdecoderList));
  figure; wf.applyFigDefaults(gcf,[nc nr],'w');

  for iDecoder = 1:numel(thisdecoderList)
    subplot(nr,nc,iDecoder); hold on
    % data
    acc     = decodeSumm.(thisdecoderList{iDecoder}).accuracy;
    accm    = squeeze(nanmean(acc,3));
    accs    = squeeze(nanstd(acc,0,3))./sqrt(size(acc,3)-1);
    thiscl  = eval(sprintf('cfg.%sCl',cfg.whichMethod));

    plot(xaxis,accm,'-','color',thiscl,'linewidth',1.5)
    plot(xaxis,accm-accs,'--','color',thiscl,'linewidth',.75)
    plot(xaxis,accm+accs,'--','color',thiscl,'linewidth',.75)

    % shuffle
    acc     = decodeSumm.(thisdecoderList{iDecoder}).accuracy_shuffle;
    accm    = squeeze(nanmean(acc,3));
    accs    = squeeze(nanstd(acc,0,3))./sqrt(size(acc,3)-1);
    thiscl  = eval(sprintf('cfg.%sCl',cfg.whichMethod));

    plot(xaxis,accm,'-','color',thiscl,'linewidth',1.5)
    plot(xaxis,accm-accs,'--','color',thiscl,'linewidth',.75)
    plot(xaxis,accm+accs,'--','color',thiscl,'linewidth',.75)


    % labels
    switch cfg.whichDecoder
      case 'evidence'
        titlestr = 'Evidence (Cumul. \Delta)';
        ylbl     = 'Cross-val. acc. (r)';
        ylim([-.05 .35])
      
      case 'absevidence'
        titlestr = 'Evidence (|Cumul. \Delta|)';
        ylbl     = 'Cross-val. acc. (r)';
        ylim([-.05 .35])
      
      case 'decisionPoint'
        titlestr = 'Decis. point (cm)';
        ylbl     = 'Cross-val. acc. (r)';
        ylim([-.05 .35])
        
      case 'viewangle'
        titlestr = 'View angle';
        ylbl     = 'Cross-val. acc. (r)';
        ylim([-.05 .35])

      case 'choice'
        titlestr = 'Upcoming choice';
        ylbl     = 'Cross-val. acc. (prop. correct)';
        ylim([.45 .85])

      case 'prevchoice'
        titlestr = 'Previous choice';
        ylbl     = 'Cross-val. acc. (prop. correct)';
        ylim([.45 .75])
    end

    switch cfg.timeOrSpace
      case 'time'
        xlbl     = 'Time (s)';

      case 'space'
        xlbl     = 'y position (cm)';

    end

    wf.applyAxisDefaults(gca,'k');
    wf.applyAxisLbls(gca,xlbl,ylbl,titlestr)
    set(gca,'ytick',0:.1:1)
    
    if ~isempty(strfind(thisdecoderList{iDecoder},'Pxls')); continue; end
    
    % plot avg weights
    subplot(nr,nc,iDecoder+numel(thisdecoderList)); hold on
    nROI    = numel(ROIlbl);
    weights = nanmean(decodeSumm.(thisdecoderList{iDecoder}).weights,3);
    imagesc(xaxis,1:nROI,nanmean(weights,3)',[-abs(max(weights(:))) abs(max(weights(:)))]); 
    axis tight; colormap red2blue
    set(gca,'ytick',1:nROI,'yticklabel',ROIlbl)
    wf.applyAxisDefaults(gca,'k'); 
    wf.applyAxisLbls(gca,xlbl,[],'Decode weights')
    smallcolorbar(gca);
    
    subplot(nr,nc,numel(thisdecoderList)*2+iDecoder); hold on
    weightmean = mean(weights);
    weightsem  = std(weights)./sqrt(numel(mice)-1);
    for iROI = 1:nROI
      [cl,lbls{iROI}] = getDefaultROIcl(ROIlbl{iROI});
      bar(iROI,weightmean(iROI),'edgecolor',cl,'facecolor',cl);
      errorbar(iROI,weightmean(iROI),weightsem(iROI),'-','color',cl);
    end
    set(gca,'xtick',1:nROI,'xticklabel',lbls)
    rotateXLabels(gca,90)
    wf.applyAxisDefaults(gca,'k'); axis tight
    wf.applyAxisLbls(gca,[],'Decoding weight','Average weight')
  end
  
  %% all ROIs, dendrograms
  figure; wf.applyFigDefaults(gcf,[4 1],'w');
  subplot(1,4,1)
  dendrogram(decodeSumm.allROIs.clustStats.lnk,'labels',ROIlbl)
  wf.applyAxisDefaults(gca,'k');
  wf.applyAxisLbls(gca,[],'eucl. dist','weights')
  rotateXLabels(gca,90)
  
  subplot(1,4,2)
  dendrogram(decodeSumm.allROIs.clustStats_absW.lnk,'labels',ROIlbl)
  wf.applyAxisDefaults(gca,'k');
  wf.applyAxisLbls(gca,[],'eucl. dist','abs(weights)')
  rotateXLabels(gca,90)
  
  subplot(1,4,3)
  dendrogram(decodeSumm.allROIs_vAngResid.clustStats.lnk,'labels',ROIlbl)
  wf.applyAxisDefaults(gca,'k');
  wf.applyAxisLbls(gca,[],'eucl. dist','weights')
  rotateXLabels(gca,90)
  
  subplot(1,4,4)
  dendrogram(decodeSumm.allROIs_vAngResid.clustStats_absW.lnk,'labels',ROIlbl)
  wf.applyAxisDefaults(gca,'k');
  wf.applyAxisLbls(gca,[],'eucl. dist','abs(weights)')
  rotateXLabels(gca,90)
  
  %% each ROI, accuracy
  thisdecoderList = {'ROIpxls','ROIpxls_vAngResid'};

  for iDecoder = 1:numel(thisdecoderList)
    nROI = size(decodeSumm.(thisdecoderList{iDecoder}).accuracy,1);
    if nROI == 8
      roilbls = decodeSumm.ROIlbl_unil;
    else
      roilbls = decodeSumm.ROIlbl;
    end
    [nr,nc] = subplotOrg(nROI,4);
    figure; wf.applyFigDefaults(gcf,[nc nr],'w');
    
    for iROI = 1:nROI
      subplot(nr,nc,iROI); hold on
      % data
      acc     = decodeSumm.(thisdecoderList{iDecoder}).accuracy(iROI,:,:);
      accm    = squeeze(nanmean(acc,3));
      accs    = squeeze(nanstd(acc,0,3))./sqrt(size(acc,3)-1);
      thiscl  = eval(sprintf('cfg.%sCl',cfg.whichMethod));

      plot(xaxis,accm,'-','color',thiscl,'linewidth',1.5)
      plot(xaxis,accm-accs,'--','color',thiscl,'linewidth',.75)
      plot(xaxis,accm+accs,'--','color',thiscl,'linewidth',.75)

      % shuffle
      acc     = decodeSumm.(thisdecoderList{iDecoder}).accuracy_shuffle(iROI,:,:);
      accm    = squeeze(nanmean(acc,3));
      accs    = squeeze(nanstd(acc,0,3))./sqrt(size(acc,3)-1);
      thiscl  = eval(sprintf('cfg.%sCl',cfg.whichMethod));

      plot(xaxis,accm,'-','color',thiscl,'linewidth',1.5)
      plot(xaxis,accm-accs,'--','color',thiscl,'linewidth',.75)
      plot(xaxis,accm+accs,'--','color',thiscl,'linewidth',.75)

      % labels
      switch cfg.whichDecoder
        case 'evidence'
          titlestr = 'Evidence (Cumul. \Delta)';
          ylbl     = 'Cross-val. acc. (r)';
          ylim([-.05 .35])
          
        case 'absevidence'
          titlestr = 'Evidence (|Cumul. \Delta|)';
          ylbl     = 'Cross-val. acc. (r)';
          ylim([-.05 .35])
        
      case 'decisionPoint'
          titlestr = 'Decis. point (cm)';
          ylbl     = 'Cross-val. acc. (r)';
          ylim([-.05 .35])
        
        case 'viewangle'
          titlestr = 'View angle';
          ylbl     = 'Cross-val. acc. (r)';
          ylim([-.05 .35])

        case 'choice'
          titlestr = 'Upcoming choice';
          ylbl     = 'Cross-val. acc. (prop. correct)';
          ylim([.45 .85])

        case 'prevchoice'
          titlestr = 'Previous choice';
          ylbl     = 'Cross-val. acc. (prop. correct)';
          ylim([.45 .75])
      end

      switch cfg.timeOrSpace
        case 'time'
          xlbl     = 'Time (s)';

        case 'space'
          xlbl     = 'y position (cm)';

      end

      wf.applyAxisDefaults(gca,'k');
      wf.applyAxisLbls(gca,xlbl,ylbl, ...
        sprintf('%s\n%s, %s',titlestr,removeUnderscores(thisdecoderList{iDecoder}),roilbls{iROI}))
      set(gca,'ytick',0:.1:1)
    end
  end
  
  %% decoder type comparisons
  test_types   = {'accuracy_avg','accuracy_max','accuracy_maxy','accuracy_firsty'};
  clear comps
  comps{1}     = {'ROIpxls','allPxls'};

  [nr,nc] = subplotOrg(numel(test_types)*numel(comps),numel(test_types));
  figure; wf.applyFigDefaults(gcf,[nc nr],'w');
  ct = 1;  
  for iComp = 1:numel(comps)
    for iType = 1:numel(test_types)
      subplot(nr,nc,ct); hold on
      allROIvec = decodeSumm.(comps{iComp}{2}).(test_types{iType});
      bar(1,nanmean(allROIvec),'facecolor','k','edgecolor','k');
      errorbar(1,nanmean(allROIvec),nanstd(allROIvec)./sqrt(numel(mice)-1),'k-');
      plotlbls{1}  = 'all';
      
      datavec      = decodeSumm.(comps{iComp}{1}).(test_types{iType});
      nROI         = size(datavec,2);
      if nROI == 8
        roilbls = decodeSumm.ROIlbl_unil;
      else
        roilbls = decodeSumm.ROIlbl;
      end
      for iROI = 1:nROI
        thismean = mean(datavec(:,iROI));
        thissem  = std(datavec(:,iROI))./sqrt(numel(mice)-1);
        [cl,plotlbls{end+1}] = getDefaultROIcl(roilbls{iROI});
        bar(iROI+1,thismean,'edgecolor',cl,'facecolor',cl);
        errorbar(iROI+1,thismean,thissem,'-','color',cl);
      end
      
      set(gca,'xtick',1:nROI+1,'xticklabel',plotlbls)
      wf.applyAxisDefaults(gca,'k');
      wf.applyAxisLbls(gca,[],removeUnderscores(test_types{iType}),removeUnderscores(comps{iComp}{1}));
      
      ct = ct+1;
    end
  end
  
else
  %% plot decoding accuracy
  [nr,nc] = subplotOrg(numel(decoderList),4);
  figure; wf.applyFigDefaults(gcf,[nc nr],'w');

  for iDecoder = 1:numel(decoderList)
    subplot(nr,nc,iDecoder); hold on
    % data
    acc     = decodeSumm.(decoderList{iDecoder}).accuracy;
    accm    = squeeze(nanmean(acc,3));
    accs    = squeeze(nanstd(acc,0,3))./sqrt(size(acc,3)-1);
    thiscl  = eval(sprintf('cfg.%sCl',cfg.whichMethod));

    plot(xaxis,accm,'-','color',thiscl,'linewidth',1.5)
    plot(xaxis,accm-accs,'--','color',thiscl,'linewidth',.75)
    plot(xaxis,accm+accs,'--','color',thiscl,'linewidth',.75)

    % shuffle
    acc     = decodeSumm.(decoderList{iDecoder}).accuracy_shuffle;
    accm    = squeeze(nanmean(acc,3));
    accs    = squeeze(nanstd(acc,0,3))./sqrt(size(acc,3)-1);
    thiscl  = eval(sprintf('cfg.%sCl',cfg.whichMethod));

    plot(xaxis,accm,'-','color',thiscl,'linewidth',1.5)
    plot(xaxis,accm-accs,'--','color',thiscl,'linewidth',.75)
    plot(xaxis,accm+accs,'--','color',thiscl,'linewidth',.75)


    % labels
    switch cfg.whichDecoder
      case 'evidence'
        titlestr = 'Evidence (Cumul. \Delta)';
        ylbl     = 'Cross-val. acc. (r)';
        ylim([-.05 .55])
      
      case 'absevidence'
        titlestr = 'Evidence (|Cumul. \Delta|)';
        ylbl     = 'Cross-val. acc. (r)';
        ylim([-.05 .55])
        
      
      case 'decisionPoint'
        titlestr = 'Decis. point (cm)';
        ylbl     = 'Cross-val. acc. (r)';
        ylim([-.05 .35])
      
      case 'viewangle'
        titlestr = 'View angle';
        ylbl     = 'Cross-val. acc. (r)';
        ylim([-.05 .55])

      case 'choice'
        titlestr = 'Upcoming choice';
        ylbl     = 'Cross-val. acc. (prop. correct)';
        ylim([.48 .75])

      case 'prevchoice'
        titlestr = 'Previous choice';
        ylbl     = 'Cross-val. acc. (prop. correct)';
        ylim([.48 .75])
    end

    switch cfg.timeOrSpace
      case 'time'
        xlbl     = 'Time (s)';

      case 'space'
        xlbl     = 'y position (cm)';

    end

    wf.applyAxisDefaults(gca,'k');
    wf.applyAxisLbls(gca,xlbl,ylbl,titlestr)
    set(gca,'ytick',0:.1:1)

    % legend
    if numel(cfg.whichMethod) > 1
      for iMethod = 1:numel(cfg.whichMethod)
        xl      = get(gca,'xlim');
        yl      = get(gca,'ylim');
        thiscl  = eval(sprintf('cfg.%sCl',cfg.whichMethod{iMethod}));
        text(xl(1)+diff(xl)*.05 , yl(2) - iMethod*diff(yl)*.1 , ...
          cfg.whichMethod, 'color', thiscl, 'fontsize', 14)
      end
    end
  end
end

%% plot spatial correlation of pxl weights
figure; wf.applyFigDefaults(gcf,[1 1],'w');
hold on; 
thisx = decodeSumm.allPxls.pxlDistAxis;
thism = nanmean(decodeSumm.allPxls.weightSpatialCorr,2);
thiss = nanstd(decodeSumm.allPxls.weightSpatialCorr,0,2)./sqrt(numel(mice)-1);
plot(thisx,thism,'k-','linewidth',1); 
plot(thisx,thism-thiss,'k--','linewidth',.5); 
plot(thisx,thism+thiss,'k--','linewidth',.5); 

thism = nanmean(decodeSumm.allPxls.weightSpatialCorr_shuffle,2);
thiss = nanstd(decodeSumm.allPxls.weightSpatialCorr_shuffle,0,2)./sqrt(numel(mice)-1);
plot(thisx,thism,'-','linewidth',1,'color',[.5 .5 .5]); 
plot(thisx,thism-thiss,'--','linewidth',.5,'color',[.5 .5 .5]); 
plot(thisx,thism+thiss,'--','linewidth',.5,'color',[.5 .5 .5]); 
wf.applyAxisDefaults(gca,'k');
wf.applyAxisLbls(gca,'Pxl dist (um)','CC','weight correlation, pxl decoder')
    
%%
%% export
if ~isempty(dir([fn '.pdf'])); delete([fn '.pdf']); end
figls = get(0,'children'); % print to pdf
for ii = length(figls):-1:1
  figure(figls(ii))
  export_fig([fn '.pdf'],'-q101','-append')
end
close all
warning on
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
if ~isfield(cfg,'whichDecoder')
  cfg(1).whichDecoder = 'evidence'; %{'evidence';'choice';'prevchoice';'viewangle'};
end
if ~isfield(cfg,'whichMethod')
  cfg(1).whichMethod = 'ridge';%{'ridge';'lasso'};
end
if ~isfield(cfg,'trialType')
  cfg(1).trialType = 'all';
end
if ~isfield(cfg,'timeOrSpace')
  cfg(1).timeOrSpace = 'space';
end
if ~isfield(cfg,'binx')
  cfg(1).binx         = 4;
end
if ~isfield(cfg,'ridgeCl')
  cfg(1).ridgeCl         = [0 0 0];
end
if ~isfield(cfg,'ridgeShuffleCl')
  cfg(1).ridgeShuffleCl  = widefieldParams.lightgray;
end
if ~isfield(cfg,'lassoCl')
  cfg(1).lassoCl         = widefieldParams.mypurple;
end
if ~isfield(cfg,'lassoShuffleCl')
  cfg(1).lassoShuffleCl  = widefieldParams.lightpurple;
end
if ~isfield(cfg,'clustMaxPC')
  cfg(1).clustMaxPC = 3;
end
if ~isfield(cfg,'clustMaxK')
  cfg(1).clustMaxK  = 5;
end

end