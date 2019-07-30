function avgDffSumm = summarizeAvgDff(mice)

% avgDffSumm = summarizeAvgDff(mice)
% for mazes and trial types compiles avg pxl/ROI analyses 
% (space avg, COM, pxl sequences, ANOVA)
% mice is cell array with mouse names

%% ------------------------------------------------------------------------

if nargin < 1; mice  = {'ai2';'ai3';'ai5';'ai7';'ai9';'ai10'}; end

try
%% flag true to load from local disk instead of bucket
% (will have no effect if running on spock)
localFlag    = false; 

%% configuration for ROI corr analysis
cfg.trialTypes         = { 'correctR','correctL','errorR','errorL',...
                           'correctR_easy','correctL_easy','errorR_easy','errorL_easy',...
                           'correctR_hard','correctL_hard','errorR_hard','errorL_hard',...
                         };
cfg.trialTypes_subt    = {{'correctR_hard','correctR_easy'}, ...
                          {'correctL_hard','correctL_easy'}, ...
                          {'correctR'     ,'errorR'       }, ...
                          {'correctL'     ,'errorL'       }, ...
                         };
cfg.maze_subt          = { 'correctR'     ,'correctL'     };
cfg.bilateralROIflag   = true;
cfg.maxBinSequence     = 300; 
cfg.plotSeqRec         = true;
cfg.mice               = mice;
wf                     = widefieldParams;

latencyCfg.trialType   = 'correct'; % pre-baseline before these trials
latencyCfg.time        = 0.5; % sec before trial start
latencyCfg.nSEM        = 3; % threshold: x s.e.m. beyond mean 
latencyCfg.nConsecBis  = 3; % num consecutive bins must cross threshold
cfg.latencyCfg         = latencyCfg;

rootdir = wf.getRootDir(isThisSpock,localFlag);
cd(rootdir)
if cfg.plotSeqRec && ~thy1Flag
  if ~isempty(dir('seqPlotsRecs.pdf')); end
end

%% ------------------------------------------------------------------------
%% compile avg t4, t11 for each mouse (session)
%% ------------------------------------------------------------------------

fprintf('retrieving data...\n')
avgDffSumm.ROIlbl  = {};

% loop through recs
for iMouse  = 1:numel(mice)
  
  %% task (visual guide and accumulation)
  % average across recs for each animal first, but keep track of x-rec
  % corrs
  taskRecs                   = wf.getMouseRecs(mice{iMouse},'task',localFlag);
  avgDffSumm.mouse(iMouse).recls = taskRecs;
  for iRec = 1:numel(taskRecs)
    fprintf('\tmouse %d/%d, rec %d/%d\n',iMouse,numel(mice),iRec,numel(taskRecs))
    cd(taskRecs{iRec})
    
    %% collect some pxl averages for the records (and examples)
    if ~thy1Flag
      load avgDff avgDff
      avgDffSumm.spaceBins = avgDff.spaceBins;
      for iType = 1:numel(cfg.trialTypes)
        avgDffSumm.mouse(iMouse).pxl.accumul.(cfg.trialTypes{iType}).recs{iRec}              ...
                           = avgDff.maze(end).(cfg.trialTypes{iType}).space.frames;
      end
      for iType = 1:numel(cfg.maze_subt)
        avgDffSumm.mouse(iMouse).pxl.visGuide.(cfg.maze_subt{iType}).recs{iRec}              ...
                           = avgDff.maze(1).(cfg.maze_subt{iType}).space.frames;
        avgDffSumm.mouse(iMouse).pxl.accumulMINUSvisGuide.(cfg.maze_subt{iType}).recs{iRec}  ...
                           = avgDff.maze(end).(cfg.maze_subt{iType}).space.frames -          ...
                             avgDff.maze(1).(cfg.maze_subt{iType}).space.frames;
      end
      for iType = 1:numel(cfg.trialTypes_subt)
        avgDffSumm.mouse(iMouse).pxl.accumul.([cfg.trialTypes_subt{iType}{1} 'MINUS' cfg.trialTypes_subt{iType}{2}]).recs{iRec}  ...
                           = avgDff.maze(end).(cfg.trialTypes_subt{iType}{1}).space.frames - ...
                             avgDff.maze(end).(cfg.trialTypes_subt{iType}{2}).space.frames;
      end


      %% compare avg pixel sequences for visGuide and accumul
      avgDffSumm.mouse(iMouse).pxl.seq.rec(iRec) = activitySequences(avgDff,cfg,cfg.plotSeqRec);
      if cfg.plotSeqRec; export_fig([rootdir 'seqPlotsRecs.pdf'],'-q101','-append'); close all; end
      clear avgDff
    end
    
    %% collect ROI averages, subtraction of averages, center of mass/peak
    load avgROI avgROI
    avgDff = avgROI;
    xaxis  = avgDff.spaceBins;
    isGoodBin = xaxis <= cfg.maxBinSequence;
    for iType = 1:numel(cfg.trialTypes)
      mat    = avgDff.maze(end).(cfg.trialTypes{iType}).space.frames;
      if isempty(mat)
        mat  = nan(size(avgDff.maze(end).(cfg.trialTypes{iType-1}).space.frames));
      end
      avgDffSumm.mouse(iMouse).ROI.accumul.(cfg.trialTypes{iType}).recs(:,:,iRec)              ...
             = mat;
      vals   = COM(mat(isGoodBin,:)',xaxis(isGoodBin));
      if sum(sum(isnan(mat))) > 1
        valsp  = nan(size(mat,2),1);
      else
        valsp  = arrayfun(@(x)(xaxis(find(mat(:,x)==max(mat(:,x)),1,'first'))),1:size(mat,2));
      end
      avgDffSumm.mouse(iMouse).ROI.accumul.(cfg.trialTypes{iType}).COM.recs(:,iRec)            ...
             = vals;
      avgDffSumm.mouse(iMouse).ROI.accumul.(cfg.trialTypes{iType}).peak.recs(:,iRec)           ...
             = valsp;
      
      if isempty(strfind(cfg.trialTypes{iType},'correct'))
        latencyCfg.trialType = 'error';
      else
        latencyCfg.trialType = 'correct';
      end
      avgDffSumm.mouse(iMouse).ROI.accumul.(cfg.trialTypes{iType}).latency.recs(:,iRec)           ...
             = calculateLatency(mat,xaxis,latencyCfg,pwd);
    end
    for iType = 1:numel(cfg.maze_subt)
      mat    = avgDff.maze(1).(cfg.maze_subt{iType}).space.frames;
      vals   = COM(mat(isGoodBin,:)',xaxis(isGoodBin));
      valsp  = arrayfun(@(x)(xaxis(find(mat(:,x)==max(mat(:,x)),1,'first'))),1:size(mat,2));
      avgDffSumm.mouse(iMouse).ROI.visGuide.(cfg.maze_subt{iType}).recs(:,:,iRec)              ...
             = mat;
      avgDffSumm.mouse(iMouse).ROI.visGuide.(cfg.maze_subt{iType}).COM.recs(:,iRec)            ...
             = vals;
      avgDffSumm.mouse(iMouse).ROI.visGuide.(cfg.maze_subt{iType}).peak.recs(:,iRec)           ...
             = valsp;
      latencyCfg.trialType = 'correct';     
      avgDffSumm.mouse(iMouse).ROI.visGuide.(cfg.maze_subt{iType}).latency.recs(:,iRec)        ...
               = calculateLatency(mat,xaxis,latencyCfg,pwd,4);
      avgDffSumm.mouse(iMouse).ROI.accumulMINUSvisGuide.(cfg.maze_subt{iType}).recs(:,:,iRec)  ...
             = avgDff.maze(end).(cfg.maze_subt{iType}).space.frames - mat;
    end
    for iType = 1:numel(cfg.trialTypes_subt)
      avgDffSumm.mouse(iMouse).ROI.accumul.([cfg.trialTypes_subt{iType}{1} 'MINUS' cfg.trialTypes_subt{iType}{2}]).recs(:,:,iRec)  ...
                         = avgDffSumm.mouse(iMouse).ROI.accumul.(cfg.trialTypes_subt{iType}{1}).recs(:,:,iRec) -     ...
                           avgDffSumm.mouse(iMouse).ROI.accumul.(cfg.trialTypes_subt{iType}{2}).recs(:,:,iRec);
    end
    
    load ROIfromRef ROI ROIlbl
    avgDffSumm.ROIlbl = ROIlbl;
    %% compare avg pixel sequences for visGuide and accumul
    if ~thy1Flag
      avgDffSumm.mouse(iMouse).ROI.seq.rec(iRec) = activitySequences(avgDff,cfg,cfg.plotSeqRec);
      if cfg.plotSeqRec; export_fig([rootdir 'seqPlotsRecs.pdf'],'-q101','-append'); close all; end
      clear avgDff avgROI

      %% compare % task-modulated pixels for each ROI (ANOVA), visGuide and accumul
      load dffANOVA3 anovaData
      ROIpercent        = getPercModPixels(anovaData.modMat,ROI,cfg.bilateralROIflag);
      avgDffSumm.mouse(iMouse).ROI.accumul.percModulatedPixels.rec(:,iRec)  = ROIpercent;

      load dffANOVA3_visGuide anovaData
      ROIpercent        = getPercModPixels(anovaData.modMat,ROI,cfg.bilateralROIflag);
      avgDffSumm.mouse(iMouse).ROI.visGuide.percModulatedPixels.rec(:,iRec) = ROIpercent;

      clear ROI ROIlbl anovaData
    end
  end
  
  %% calculate averages, sem and corr across recs for this mouse
  for iType = 1:numel(cfg.trialTypes)
    mat    = avgDffSumm.mouse(iMouse).ROI.accumul.(cfg.trialTypes{iType}).recs;
    avgDffSumm.(['ROIavg_accumul_' (cfg.trialTypes{iType})])(:,:,iMouse)        = nanmean(mat,3);
    avgDffSumm.(['ROIavg_accumul_' (cfg.trialTypes{iType}) '_sem'])(:,:,iMouse) = nanstd(mat,0,3)./sqrt(size(mat,3)-1);
    
    [nx,ny,nz]  = size(mat);
    mat     = reshape(mat,[nx*ny nz]);
    thiscorr                                    = corr(mat);
    thiscorr(triu(true(size(thiscorr)),0))      = nan; % since corr is symmetrical dont double dip ROI pair
    thiscorr                                    = thiscorr(:);
    avgDffSumm.stats.(['xRecCorr_accumul_' (cfg.trialTypes{iType})]).mouse(iMouse).pairs = thiscorr;
    avgDffSumm.stats.(['xRecCorr_accumul_' (cfg.trialTypes{iType})]).cc_avg(iMouse,:)    = nanmean(thiscorr);
  
    mat    = avgDffSumm.mouse(iMouse).ROI.accumul.(cfg.trialTypes{iType}).COM.recs;
    avgDffSumm.(['ROIavg_accumul_COM_' (cfg.trialTypes{iType})])(:,iMouse)        = nanmean(mat,2);
    avgDffSumm.(['ROIavg_accumul_COM_' (cfg.trialTypes{iType}) '_sem'])(:,iMouse) = nanstd(mat,0,2)./sqrt(size(mat,2)-1);
    thiscorr                                    = corr(mat);
    thiscorr(triu(true(size(thiscorr)),0))      = nan; % since corr is symmetrical dont double dip ROI pair
    thiscorr                                    = thiscorr(:);
    avgDffSumm.stats.(['xRecCorr_accumul_COM_' (cfg.trialTypes{iType})]).mouse(iMouse).pairs = thiscorr;
    avgDffSumm.stats.(['xRecCorr_accumul_COM_' (cfg.trialTypes{iType})]).cc_avg(iMouse,:)    = nanmean(thiscorr);
    
    mat    = avgDffSumm.mouse(iMouse).ROI.accumul.(cfg.trialTypes{iType}).latency.recs;
    avgDffSumm.(['ROIavg_accumul_latency_' (cfg.trialTypes{iType})])(:,iMouse)        = nanmean(mat,2);
    avgDffSumm.(['ROIavg_accumul_latency_' (cfg.trialTypes{iType}) '_sem'])(:,iMouse) = nanstd(mat,0,2)./sqrt(size(mat,2)-1);
    
    mat    = avgDffSumm.mouse(iMouse).ROI.accumul.(cfg.trialTypes{iType}).peak.recs;
    avgDffSumm.(['ROIavg_accumul_peak_' (cfg.trialTypes{iType})])(:,iMouse)        = nanmean(mat,2);
    avgDffSumm.(['ROIavg_accumul_peak_' (cfg.trialTypes{iType}) '_sem'])(:,iMouse) = nanstd(mat,0,2)./sqrt(size(mat,2)-1);
  end
  for iType = 1:numel(cfg.maze_subt)
    mat    = avgDffSumm.mouse(iMouse).ROI.visGuide.(cfg.maze_subt{iType}).recs;
    avgDffSumm.(['ROIavg_visGuide_' (cfg.maze_subt{iType})])(:,:,iMouse)        = nanmean(mat,3);
    avgDffSumm.(['ROIavg_visGuide_' (cfg.maze_subt{iType}) '_sem'])(:,:,iMouse) = nanstd(mat,0,3)./sqrt(size(mat,3)-1);
    
    [nx,ny,nz]  = size(mat);
    mat     = reshape(mat,[nx*ny nz]);
    thiscorr                                    = corr(mat);
    thiscorr(triu(true(size(thiscorr)),0))      = nan; % since corr is symmetrical dont double dip ROI pair
    thiscorr                                    = thiscorr(:);
    avgDffSumm.stats.(['xRecCorr_visGuide_' (cfg.maze_subt{iType})]).mouse(iMouse).pairs = thiscorr;
    avgDffSumm.stats.(['xRecCorr_visGuide_' (cfg.maze_subt{iType})]).cc_avg(iMouse,:)    = nanmean(thiscorr);
  
    mat    = avgDffSumm.mouse(iMouse).ROI.visGuide.(cfg.maze_subt{iType}).COM.recs;
    avgDffSumm.(['ROIavg_visGuide_COM_' (cfg.trialTypes{iType})])(:,iMouse)        = nanmean(mat,2);
    avgDffSumm.(['ROIavg_visGuide_COM_' (cfg.trialTypes{iType}) '_sem'])(:,iMouse) = nanstd(mat,0,2)./sqrt(size(mat,2)-1);
    thiscorr                                    = corr(mat);
    thiscorr(triu(true(size(thiscorr)),0))      = nan; % since corr is symmetrical dont double dip ROI pair
    thiscorr                                    = thiscorr(:);
    avgDffSumm.stats.(['xRecCorr_visGuide_COM_' (cfg.maze_subt{iType})]).mouse(iMouse).pairs = thiscorr;
    avgDffSumm.stats.(['xRecCorr_visGuide_COM_' (cfg.maze_subt{iType})]).cc_avg(iMouse,:)    = nanmean(thiscorr);
    
    mat    = avgDffSumm.mouse(iMouse).ROI.visGuide.(cfg.trialTypes{iType}).latency.recs;
    avgDffSumm.(['ROIavg_visGuide_latency_' (cfg.maze_subt{iType})])(:,iMouse)        = nanmean(mat,2);
    avgDffSumm.(['ROIavg_visGuide_latency_' (cfg.maze_subt{iType}) '_sem'])(:,iMouse) = nanstd(mat,0,2)./sqrt(size(mat,2)-1);

    mat    = avgDffSumm.mouse(iMouse).ROI.visGuide.(cfg.trialTypes{iType}).peak.recs;
    avgDffSumm.(['ROIavg_visGuide_peak_' (cfg.maze_subt{iType})])(:,iMouse)        = nanmean(mat,2);
    avgDffSumm.(['ROIavg_visGuide_peak_' (cfg.maze_subt{iType}) '_sem'])(:,iMouse) = nanstd(mat,0,2)./sqrt(size(mat,2)-1);
      
    mat    = avgDffSumm.mouse(iMouse).ROI.accumulMINUSvisGuide.(cfg.maze_subt{iType}).recs;
    avgDffSumm.(['ROIavg_accumulMINUSvisGuide_' (cfg.maze_subt{iType})])(:,:,iMouse)        = nanmean(mat,3);
    avgDffSumm.(['ROIavg_accumulMINUSvisGuide_' (cfg.maze_subt{iType}) '_sem'])(:,:,iMouse) = nanstd(mat,0,3)./sqrt(size(mat,3)-1);
  end
  for iType = 1:numel(cfg.trialTypes_subt)
    mat    = avgDffSumm.mouse(iMouse).ROI.accumul.([cfg.trialTypes_subt{iType}{1} 'MINUS' cfg.trialTypes_subt{iType}{2}]).recs;
    avgDffSumm.(['ROIavg_accumul_' ([cfg.trialTypes_subt{iType}{1} 'MINUS' cfg.trialTypes_subt{iType}{2}])])(:,:,iMouse)        = nanmean(mat,3);
    avgDffSumm.(['ROIavg_accumul_' ([cfg.trialTypes_subt{iType}{1} 'MINUS' cfg.trialTypes_subt{iType}{2}]) '_sem'])(:,:,iMouse) = nanstd(mat,0,3)./sqrt(size(mat,3)-1);
  end
  
  if ~thy1Flag
    avgDffSumm.percModulatedPixels_accumul(:,iMouse)      = nanmean(avgDffSumm.mouse(iMouse).ROI.accumul.percModulatedPixels.rec,2);
    avgDffSumm.percModulatedPixels_accumul_sem(:,iMouse)  = nanstd(avgDffSumm.mouse(iMouse).ROI.accumul.percModulatedPixels.rec,0,2)./sqrt(size(mat,3)-1);
    avgDffSumm.percModulatedPixels_visGuide(:,iMouse)     = mean(avgDffSumm.mouse(iMouse).ROI.visGuide.percModulatedPixels.rec,2);
    avgDffSumm.percModulatedPixels_visGuide_sem(:,iMouse) = nanstd(avgDffSumm.mouse(iMouse).ROI.visGuide.percModulatedPixels.rec,0,2)./sqrt(size(mat,3)-1);

    types = {'pxl','ROI'};
    nrecs = numel(avgDffSumm.mouse(iMouse).pxl.seq.rec);
    for iType = 1:numel(types)
      mat = zeros(1,nrecs);
      for iRec = 1:nrecs; mat(iRec) = avgDffSumm.mouse(iMouse).(types{iType}).seq.rec(iRec).COMcorr_accumulVSvisGuide; end
      avgDffSumm.([types{iType} 'Seq_COMcorr_accumulVSvisGuide'])(iMouse,:)        = nanmean(mat);
      avgDffSumm.([types{iType} 'Seq_COMcorr_accumulVSvisGuide_sem'])(iMouse,:)    = nanstd(mat)./sqrt(numel(mat)-1);

      mat = zeros(1,nrecs);
      for iRec = 1:nrecs; mat(iRec) = avgDffSumm.mouse(iMouse).(types{iType}).seq.rec(iRec).accumul.COMcorr_RvsL; end
      avgDffSumm.([types{iType} 'Seq_COMcorr_accumul_RvsL'])(iMouse,:)      = nanmean(mat);
      avgDffSumm.([types{iType} 'Seq_COMcorr_accumul_RvsL_sem'])(iMouse,:)  = nanstd(mat)./sqrt(numel(mat)-1);

      mat = zeros(1,nrecs);
      for iRec = 1:nrecs; mat(iRec) = avgDffSumm.mouse(iMouse).(types{iType}).seq.rec(iRec).visGuide.COMcorr_RvsL; end
      avgDffSumm.([types{iType} 'Seq_COMcorr_visGuide_RvsL'])(iMouse,:)     = nanmean(mat);
      avgDffSumm.([types{iType} 'Seq_COMcorr_visGuide_RvsL_sem'])(iMouse,:) = nanstd(mat)./sqrt(numel(mat)-1);

      mat = zeros(1,nrecs);
      for iRec = 1:nrecs; mat(iRec) = avgDffSumm.mouse(iMouse).(types{iType}).seq.rec(iRec).OVERALLcorr_accumulVSvisGuide; end
      avgDffSumm.([types{iType} 'Seq_OVERALLcorr_accumulVSvisGuide'])(iMouse,:)        = nanmean(mat);
      avgDffSumm.([types{iType} 'Seq_OVERALLcorr_accumulVSvisGuide_sem'])(iMouse,:)    = nanstd(mat)./sqrt(numel(mat)-1);

      mat = zeros(1,nrecs);
      for iRec = 1:nrecs; mat(iRec) = avgDffSumm.mouse(iMouse).(types{iType}).seq.rec(iRec).accumul.OVERALLcorr_RvsL; end
      avgDffSumm.([types{iType} 'Seq_OVERALLcorr_accumul_RvsL'])(iMouse,:)      = nanmean(mat);
      avgDffSumm.([types{iType} 'Seq_OVERALLcorr_accumul_RvsL_sem'])(iMouse,:)  = nanstd(mat)./sqrt(numel(mat)-1);

      mat = zeros(1,nrecs);
      for iRec = 1:nrecs; mat(iRec) = avgDffSumm.mouse(iMouse).(types{iType}).seq.rec(iRec).visGuide.OVERALLcorr_RvsL; end
      avgDffSumm.([types{iType} 'Seq_OVERALLcorr_visGuide_RvsL'])(iMouse,:)     = nanmean(mat);
      avgDffSumm.([types{iType} 'Seq_OVERALLcorr_visGuide_RvsL_sem'])(iMouse,:) = nanstd(mat)./sqrt(numel(mat)-1);
    end
  end
end

cd(rootdir)
save temp_avgSumm.mat -v7.3

%% ------------------------------------------------------------------------
%% stats
%% ------------------------------------------------------------------------

%% x-rec correlations across ROI avg activity and COM, accumulation
if ~thy1Flag
  allmat = [];
  allcom = [];
  for iMouse = 1:numel(avgDffSumm.mouse)
    matR   = avgDffSumm.mouse(iMouse).ROI.accumul.correctR.COM.recs;
    matL   = avgDffSumm.mouse(iMouse).ROI.accumul.correctR.COM.recs;
    allcom(:,end+1:end+size(matR,2)) = [matR; matL]; 

    matR            = avgDffSumm.mouse(iMouse).ROI.accumul.correctR.recs;
    matL            = avgDffSumm.mouse(iMouse).ROI.accumul.correctL.recs;
    mat         = [matR; matL];
    [nx,ny,nz]      = size(mat);
    allmat(:,end+1:end+nz) = reshape(mat,[nx*ny nz]); 
  end
  allcc                                          = corr(allmat); % rec vs rec corr
  allcc(triu(true(size(thiscorr)),0))            = nan;
  allcc                                          = allcc(:);
  avgDffSumm.stats.cc_ROIavg_accumul_allrecs     = allcc(~isnan(allcc));

  allcc                                          = corr(allcom); % rec vs rec corr
  allcc(triu(true(size(thiscorr)),0))            = nan;
  allcc                                          = allcc(:);
  avgDffSumm.stats.cc_ROIavg_COM_accumul_allrecs = allcc(~isnan(allcc));


  %% x-animal correlations (ie consistency of finding)
  % calculate average and corr across recs for this mouse
  fls = fields(avgDffSumm);
  for iF = 1:numel(fls)
    if ~isempty(strfind(fls{iF},'COM_')) && isempty(strfind(fls{iF},'sem')) && ~isempty(strfind(fls{iF},'ROIavg'))
      mat  = avgDffSumm.(fls{iF});

    elseif isempty(strfind(fls{iF},'COM_')) && isempty(strfind(fls{iF},'sem')) && ~isempty(strfind(fls{iF},'ROIavg'))
      mat        = avgDffSumm.(fls{iF});
      [nx,ny,nz] = size(mat);
      mat        = reshape(mat,[nx*ny nz]);

    else
      continue
    end

    mat(isnan(sum(mat,2)),:)                    = [];
    if isempty(mat)
      avgDffSumm.stats.xMouseCorr.pairs.(fls{iF})       = nan;
      avgDffSumm.stats.xMouseCorr.(fls{iF})             = nan;
      continue
    end
    thiscorr                                            = corr(mat);
    thiscorr(triu(true(size(thiscorr)),0))              = nan; % since corr is symmetrical dont double dip ROI pair
    thiscorr                                            = thiscorr(:);
    avgDffSumm.stats.xMouseCorr.pairs.(fls{iF})         = thiscorr;
    avgDffSumm.stats.xMouseCorr.(fls{iF})               = nanmean(thiscorr(:));

  end

  %% tests, t4 vs t11, correct vs error etc
  compwhat = {{'percModulatedPixels_accumul','percModulatedPixels_visGuide'},  ...
              {'pxlSeq_COMcorr_accumul_RvsL','pxlSeq_COMcorr_visGuide_RvsL'},  ...
              {'pxlSeq_OVERALLcorr_accumul_RvsL','pxlSeq_OVERALLcorr_visGuide_RvsL'},  ...
              {'ROIavg_accumul_COM_correctR','ROIavg_visGuide_COM_correctR'},  ...
              {'ROIavg_accumul_COM_correctL','ROIavg_visGuide_COM_correctL'},  ...
              {'ROIavg_accumul_peak_correctR','ROIavg_visGuide_peak_correctR'},  ...
              {'ROIavg_accumul_peak_correctL','ROIavg_visGuide_peak_correctL'},  ...
              {'ROIavg_accumul_latency_correctR','ROIavg_visGuide_latency_correctR'},  ...
              {'ROIavg_accumul_latency_correctL','ROIavg_visGuide_latency_correctL'},  ...
              {'ROIavg_accumul_COM_correctR','ROIavg_accumul_COM_errorR'},     ...
              {'ROIavg_accumul_COM_correctL','ROIavg_accumul_COM_errorL'},     ...
              {'ROIavg_accumul_COM_correctR_easy','ROIavg_accumul_COM_correctR_hard'},  ...
              {'ROIavg_accumul_COM_correctL_easy','ROIavg_accumul_COM_correctL_hard'},  ...
              {'ROIavg_accumul_correctR','ROIavg_visGuide_correctR'},  ...
              {'ROIavg_accumul_correctL','ROIavg_visGuide_correctL'},  ...
              {'ROIavg_accumul_correctR','ROIavg_accumul_errorR'},     ...
              {'ROIavg_accumul_correctL','ROIavg_accumul_errorL'},     ...
              {'ROIavg_accumul_correctR_easy','ROIavg_accumul_correctR_hard'},  ...
              {'ROIavg_accumul_correctL_easy','ROIavg_accumul_correctL_hard'},  ...
              };
  complbls  = {'percModulatedPixels_accumulVSvisGuide' , ...
               'pxlSeqCOM_RvsLcorr_accumulVSvisGuide'  , ...
               'pxlSeqOVERALL_RvsLcorr_accumulVSvisGuide'  , ...
               'ROIavg_correctR_COM_accumulVSvisGuide' , ...
               'ROIavg_correctL_COM_accumulVSvisGuide' , ...
               'ROIavg_correctR_peak_accumulVSvisGuide' , ...
               'ROIavg_correctL_peak_accumulVSvisGuide' , ...
               'ROIavg_correctR_latency_accumulVSvisGuide' , ...
               'ROIavg_correctL_latency_accumulVSvisGuide' , ...
               'ROIavg_COM_R_correctVSerror'           , ...
               'ROIavg_COM_L_correctVSerror'           , ...
               'ROIavg_COM_R_easyVShard'               , ...
               'ROIavg_COM_L_easyVShard'               , ...
               'ROIavg_correctR_avg_accumulVSvisGuide' , ...
               'ROIavg_correctL_avg_accumulVSvisGuide' , ...
               'ROIavg_avg_R_correctVSerror'           , ...
               'ROIavg_avg_L_correctVSerror'           , ...
               'ROIavg_avg_R_easyVShard'               , ...
               'ROIavg_avg_L_easyVShard'               , ...
              };

  for iComp = 1:numel(compwhat)
    c1  = avgDffSumm.(compwhat{iComp}{1});
    c2  = avgDffSumm.(compwhat{iComp}{2});
    try
      c1s = avgDffSumm.([compwhat{iComp}{1} '_sem']);
      c2s = avgDffSumm.([compwhat{iComp}{2} '_sem']);
    catch
      c1s = nan(size(c1)); c2s = c1s;
    end

    if size(c1,3) > 1 % for time series, collapse across time
      c1 = squeeze(nanmean(c1));
      c2 = squeeze(nanmean(c2));
    end

    if size(c1,2) > 1 % for ROI avgs, test both across mice and across ROIs
      c1m = mean(c1,1)';
      c2m = mean(c2,1)';
      c1r = mean(c1,2);
      c2r = mean(c2,2);
      c1rs = std(c1,0,2)./sqrt(size(c1,2)-1);
      c2rs = std(c2,0,2)./sqrt(size(c1,2)-1);
      c1ms = std(c1)'./sqrt(size(c1,1)-1);
      c2ms = std(c2)'./sqrt(size(c1,1)-1);

      avgDffSumm.stats.xmouse.(complbls{iComp}).avgs   = [c1m c2m];
      avgDffSumm.stats.xmouse.(complbls{iComp}).sems   = [c1ms  c2ms];
      avgDffSumm.stats.xmouse.(complbls{iComp}).lbls   = compwhat{iComp};

      avgDffSumm.stats.xroi.(complbls{iComp}).avgs     = [c1r c2r];
      avgDffSumm.stats.xroi.(complbls{iComp}).sems     = [c1rs  c2rs];
      avgDffSumm.stats.xroi.(complbls{iComp}).lbls     = compwhat{iComp};

      try
        if lillietest(c1m) || lillietest(c2m)
          avgDffSumm.stats.xmouse.(complbls{iComp}).pval = signrank(c1m,c2m);
          avgDffSumm.stats.xmouse.(complbls{iComp}).test = 'signrank';   
        else
          [~,avgDffSumm.stats.xmouse.(complbls{iComp}).pval] = ttest(c1m,c2m);
          avgDffSumm.stats.xmouse.(complbls{iComp}).test     = 'ttest';
        end
      catch
        avgDffSumm.stats.xmouse.(complbls{iComp}).pval = nan;
        avgDffSumm.stats.xmouse.(complbls{iComp}).test = 'insufficient data';
      end
      try
        if lillietest(c1r) || lillietest(c2r)
          avgDffSumm.stats.xroi.(complbls{iComp}).pval  = signrank(c1r,c2r);
          avgDffSumm.stats.xroi.(complbls{iComp}).test  = 'signrank';
        else
          [~,avgDffSumm.stats.xroi.(complbls{iComp}).pval]   = ttest(c1r,c2r);
          avgDffSumm.stats.xroi.(complbls{iComp}).test       = 'ttest';
        end
      catch
        avgDffSumm.stats.xroi.(complbls{iComp}).pval = nan;
        avgDffSumm.stats.xroi.(complbls{iComp}).test = 'insufficient data';
      end

    else
      avgDffSumm.stats.xmouse.(complbls{iComp}).avgs   = [c1  c2];
      avgDffSumm.stats.xmouse.(complbls{iComp}).sems   = [c1s  c2s];
      avgDffSumm.stats.xmouse.(complbls{iComp}).lbls   = compwhat{iComp};

      if lillietest(c1) || lillietest(c2)
        avgDffSumm.stats.xmouse.(complbls{iComp}).pval = signrank(c1,c2);
        avgDffSumm.stats.xmouse.(complbls{iComp}).test = 'signrank';
      else
        [~,avgDffSumm.stats.xmouse.(complbls{iComp}).pval] = ttest(c1,c2);
        avgDffSumm.stats.xmouse.(complbls{iComp}).test     = 'ttest';
      end
    end
  end
end

%% ------------------------------------------------------------------------
%% save and plot
%% ------------------------------------------------------------------------
cd(rootdir)
avgDffSumm.cfg = cfg;

save avgDffSummary avgDffSumm cfg -v7.3
delete temp_avgSumm.mat
if ~isempty(dir('avgDffSummary.pdf')); delete avgDffSummary.pdf; end

%% plot ROI avgs: correct vs error, easy vs hard, visual guide vs accumul (each + subtraction)
clear plots
plots{1} = {'accumul_correctL','visGuide_correctL'};
plots{2} = {'accumul_correctL_hard','accumul_correctL_easy'};
plots{3} = {'accumul_correctL','accumul_errorL'};

fh = figure;
wf.applyFigDefaults(fh,[5 numel(plots)+1],'w')
for iPlot = 1:numel(plots)
  avg1  = mean(avgDffSumm.(['ROIavg_' plots{iPlot}{1}]),3)';
  avg2  = mean(avgDffSumm.(['ROIavg_' plots{iPlot}{2}]),3)';
  avgs  = avg1 - avg2;
  nROI  = size(avg1,2);
  clim  = [min([avg1(:); avg2(:)]) max([avg1(:); avg2(:)])];
  clims = [-max(abs(avgs(:))) max(abs(avgs(:)))];
  
  subplot(numel(plots),3,(iPlot-1)*3+1);
  imagesc(avgDffSumm.spaceBins,1:nROI,avg1,clim); colormap(gca,parula)
  set(gca,'ytick',1:nROI,'yticklabel',avgDffSumm.ROIlbl)
  wf.applyAxisDefaults(gca,'k'); axis tight
  wf.applyAxisLbls(gca,'y pos (cm)',[],removeUnderscores(plots{iPlot}{1}))
  
  subplot(numel(plots),3,(iPlot-1)*3+2)
  imagesc(avgDffSumm.spaceBins,1:nROI,avg2,clim); colormap(gca,parula)
  set(gca,'ytick',1:nROI,'yticklabel',avgDffSumm.ROIlbl)
  wf.applyAxisDefaults(gca,'k'); axis tight
  wf.applyAxisLbls(gca,'y pos (cm)',[],removeUnderscores(plots{iPlot}{2}))
  
  ch1 = colorbar('location','eastoutside','position',[.63 .7-.3*(iPlot-1) .01 .08],'Color','k');
  ch1.Label.String = '\DeltaF/F';
  
  subplot(numel(plots),3,(iPlot-1)*3+3)
  imagesc(avgDffSumm.spaceBins,1:nROI,avgs,clims); colormap(gca,red2blue)
  set(gca,'ytick',1:nROI,'yticklabel',avgDffSumm.ROIlbl)
  wf.applyAxisDefaults(gca,'k'); axis tight
  wf.applyAxisLbls(gca,'y pos (cm)',[],'left - middle')
  
  ch2 = colorbar('location','eastoutside','position',[.93 .7-.3*(iPlot-1) .01 .08],'Color','k');
  ch2.Label.String = '\Delta(\DeltaF/F)';
end

%% plot ROI COM
clear plots
plots = {'accumul_COM_correctL','accumul_COM_correctR','visGuide_COM_correctL','visGuide_COM_correctL'};

fh = figure;
wf.applyFigDefaults(fh,[3 2],'w')
for iPlot = 1:numel(plots)
  subplot(2,2,iPlot); hold on
  
  commat = avgDffSumm.(['ROIavg_' plots{iPlot}]);
  mcom   = mean(commat,2);
  scom   = std(commat,0,2)/sqrt(size(commat,1)-1);

  for iR = 1:numel(avgDffSumm.ROIlbl)
    thism = mcom(iR);
    thiss = scom(iR);
    plot([thism-thiss thism+thiss],[iR iR],'-','color',[.8 .8 .8],'linewidth',3)
  end
  plot(mcom,1:iR,'k+')
  set(gca,'ytick',1:iR,'yticklabel',avgDffSumm.ROIlbl)
  wf.applyAxisDefaults(gca,'k'); axis ij; xlim([0 300]); ylim([.5 iR+.5])
  wf.applyAxisLbls(gca,'COM (cm)',[],removeUnderscores(plots{iPlot}))

end

%% percent modulated pixels, sequence correlations
fh      = figure;
wf.applyFigDefaults(fh,[2 2],'w')

% % modulated by ROI
subplot(2,2,[1 2]); hold on
ROIlbl     = unique(cellfun(@(x)(x(1:end-2)),avgDffSumm.ROIlbl,'uniformOutput',false));
avgAcc     = avgDffSumm.percModulatedPixels_accumul;
semAcc     = avgDffSumm.percModulatedPixels_accumul_sem;
avgVis     = avgDffSumm.percModulatedPixels_visGuide;
semVis     = avgDffSumm.percModulatedPixels_visGuide_sem;
plotROIlbl = {'V1','mV2','PPC','RSC','aM2','mM2','M1','SS'};
plotROI    = {'VISp','mV2','VISa','RSP','aMOs','mMOs','MOp','SS'};

xticks = [];
xlbls  = {};
for iROI = 1:numel(plotROI)
  idx = strcmpi(ROIlbl,plotROI{iROI});
  xt  = [(iROI-1)*2+1 (iROI-1)*2+2] + (iROI-1);
  xticks(end+1:end+2) = xt;
  xlbls(end+1:end+2)  = {[plotROIlbl{iROI} ' accumul'],[plotROIlbl{iROI} ' visGuide']};
  for iMouse = 1:size(avgAcc,2)
    errorbar(xt,[avgAcc(idx,iMouse) avgVis(idx,iMouse)],[semAcc(idx,iMouse) semVis(idx,iMouse)],...
            '-','linewidth',.75,'color',wf.areaCl(iROI,:))
  end
  pval = signrank(avgAcc(idx,:)',avgVis(idx,:)');
  yl   = get(gca,'ylim');
  text(mean(xt),yl(2)*.95,sprintf('P = %1.2g',pval),'color','k','fontsize',9,'horizontalAlignment','center')
end
xlim([xticks(1)-.5 xticks(end)+.5]); ylim([0 100])
set(gca,'xtick',xticks,'xticklabel',xlbls);
rotateXLabels(gca,60);
wf.applyAxisDefaults(gca,'k'); axis tight
wf.applyAxisLbls(gca,[],'Modulated pixels (%)')

% sequence similarities
types = {'pxl','ROI'};
for iType = 1:numel(types)
  subplot(2,2,2+iType); hold on
  
  thismean = mean(avgDffSumm.([types{iType} 'Seq_COMcorr_accumul_RvsL']));
  thissem  = std(avgDffSumm.([types{iType} 'Seq_COMcorr_accumul_RvsL'])) ./ sqrt(numel(cfg.mice)-1);
  bar(1,thismean,'edgecolor',wf.mediumgray,'facecolor',wf.mediumgray)
  errorbar(1,thismean,thissem,'-','color',wf.mediumgray)
  
  thismean = mean(avgDffSumm.([types{iType} 'Seq_COMcorr_visGuide_RvsL']));
  thissem  = std(avgDffSumm.([types{iType} 'Seq_COMcorr_visGuide_RvsL'])) ./ sqrt(numel(cfg.mice)-1);
  bar(2,thismean,'edgecolor',wf.mediumgray,'facecolor',wf.mediumgray)
  errorbar(2,thismean,thissem,'-','color',wf.mediumgray)
  
  thismean = mean(avgDffSumm.([types{iType} 'Seq_COMcorr_accumulVSvisGuide']));
  thissem  = std(avgDffSumm.([types{iType} 'Seq_COMcorr_accumulVSvisGuide'])) ./ sqrt(numel(cfg.mice)-1);
  bar(3,thismean,'edgecolor',wf.mediumgray,'facecolor',wf.mediumgray)
  errorbar(3,thismean,thissem,'-','color',wf.mediumgray)
  
  set(gca,'xtick',1:3,'xticklabel',{'accumul R vs L','visGuide R vs L','visGuide vs accumul'})
  rotateXLabels(gca,60);
  xlim([.5 3.5]); 
  wf.applyAxisDefaults(gca,'k'); axis tight
  wf.applyAxisLbls(gca,[],'COM corr. (r)',[types{iType} ' sequences'])
end

%% scatter plots (stats)
% plot crosses for mouse data (avg+/- sem)
[nr,nc] = subplotOrg(numel(compwhat),4);
fh      = figure;
wf.applyFigDefaults(fh,[nc nr],'w')

for iComp = 1:numel(compwhat)
  subplot(nr,nc,iComp); hold on
  avg  = avgDffSumm.stats.xmouse.(complbls{iComp}).avgs;
  sem  = avgDffSumm.stats.xmouse.(complbls{iComp}).sems;
  pval = avgDffSumm.stats.xmouse.(complbls{iComp}).pval;
  lbls = avgDffSumm.stats.xmouse.(complbls{iComp}).lbls;

  for iP = 1:size(avg,1)
    plot([avg(iP,1)-sem(iP,1) avg(iP,1)+sem(iP,1)],[avg(iP,2) avg(iP,2)],'k-','linewidth',1.5)
    plot([avg(iP,1) avg(iP,1)],[avg(iP,2)-sem(iP,2) avg(iP,2)+sem(iP,2)],'k-','linewidth',1.5)
  end
  yl     = get(gca,'ylim');
  xl     = get(gca,'xlim');
  xlim([min([xl yl]) max([xl yl])]);
  ylim([min([xl yl]) max([xl yl])]);
  yl     = get(gca,'ylim');
  plot(yl,yl,'--','color',[.7 .7 .7])

  text(yl(1)*1.1,yl(2)*.9,sprintf('p = %1.3g',pval))

  wf.applyAxisDefaults(gca,'k'); 
  wf.applyAxisLbls(gca,removeUnderscores(lbls{1}),removeUnderscores(lbls{2}))
end

%% plot crosses for roi data (avg+/- sem)
compwhat(2) = [];
complbls(2) = [];

[nr,nc] = subplotOrg(numel(compwhat),4);
fh      = figure;
wf.applyFigDefaults(fh,[nc nr],'w')

for iComp = 1:numel(compwhat)
  subplot(nr,nc,iComp); hold on
  avg  = avgDffSumm.stats.xroi.(complbls{iComp}).avgs;
  sem  = avgDffSumm.stats.xroi.(complbls{iComp}).sems;
  pval = avgDffSumm.stats.xroi.(complbls{iComp}).pval;
  lbls = avgDffSumm.stats.xroi.(complbls{iComp}).lbls;

  plotROIlbl  = {'V1','mV2','PPC','RSC','aM2','mM2','M1','SS'};
  plotROI     = {'VISp','mV2','VISa','RSP','aMOs','mMOs','MOp','SS'};

  for iROI = 1:numel(plotROI)
    if numel(plotROI) == numel(avgDffSumm.ROIlbl)
      idx = find(strcmpi(avgDffSumm.ROIlbl,[plotROI{iROI} '-R']) | strcmpi(avgDffSumm.ROIlbl,[plotROI{iROI} '-L']));
    else
      idx = find(strcmpi(unique(cellfun(@(x)(x(1:end-2)),avgDffSumm.ROIlbl,'UniformOutput',false)),plotROI{iROI}));
    end
    for iP = 1:numel(idx)
      plot([avg(idx(iP),1)-sem(idx(iP),1) avg(idx(iP),1)+sem(idx(iP),1)],[avg(idx(iP),2) avg(idx(iP),2)], ...
           '-','linewidth',1.5,'color',wf.areaCl(iROI,:))
      plot([avg(idx(iP),1) avg(idx(iP),1)],[avg(idx(iP),2)-sem(idx(iP),2) avg(idx(iP),2)+sem(idx(iP),2)], ...
           '-','linewidth',1.5,'color',wf.areaCl(iROI,:))
    end
  end

  yl     = get(gca,'ylim');
  xl     = get(gca,'xlim');
  xlim([min([xl yl]) max([xl yl])]);
  ylim([min([xl yl]) max([xl yl])]);
  yl     = get(gca,'ylim');
  plot(yl,yl,'--','color',[.7 .7 .7])

  for iROI = 1:numel(plotROIlbl)
    text(yl(2)*.9,yl(2)-.08.*diff(yl).*(iROI-1),plotROIlbl{iROI},'color',wf.areaCl(iROI,:))
  end

  text(yl(1)*1.1,yl(2)*.9,sprintf('p = %1.3g',pval))

  wf.applyAxisDefaults(gca,'k'); 
  wf.applyAxisLbls(gca,removeUnderscores(lbls{1}),removeUnderscores(lbls{2}))

end
  
%% overall x-rec corr 
fh = figure;
wf.applyFigDefaults(fh,[2 1],'w')

subplot(1,2,1)
h1 = histogram(avgDffSumm.stats.cc_ROIavg_accumul_allrecs,-1:.1:1);
h1.EdgeColor = wf.mediumgray;
wf.applyAxisDefaults(gca,'k'); xlim([-1 1])
wf.applyAxisLbls(gca,'r (accumul overall)','rec pairs','avg activity')

subplot(1,2,2)
h1 = histogram(avgDffSumm.stats.cc_ROIavg_COM_accumul_allrecs,-1:.1:1);
h1.EdgeColor = wf.mediumgray;
wf.applyAxisDefaults(gca,'k'); xlim([-1 1])
wf.applyAxisLbls(gca,'r (accumul overall)','rec pairs','COM')

%% export
figls = get(0,'children'); % print to pdf
for ii = length(figls):-1:1
  figure(figls(ii))
  export_fig avgDffSummary.pdf -q101 -append
end
close all

catch ME
  displayException(ME)
end

end

%% ------------------------------------------------------------------------
%% percentage of modulated pixels in ROI
%% ------------------------------------------------------------------------
function ROIpercent = getPercModPixels(modmat,ROI,bilateralFlag)

% dffSingleROI = getSingleROIpxlDff(dff,targetROIlbl,bilateralFlag)
% returns time x pxl matrix for subset of pixels belonging to ROI (from
% both hemispheres if bilateral flag is set to true (default)
% if single ROI returns matrix, if multiple, returns cell array of matrices
% LP dec 2017 lpinto@princeton.edu

%% input defaults
if nargin < 3; bilateralFlag = true; end

%% get relevant (bilateral) ROIs
if bilateralFlag
  unilROI = ROI;
  clear ROI; ct = 1;
  for iROI = 1:2:numel(unilROI)
    ROI{ct} = [unilROI{iROI}; unilROI{iROI+1}];
    ct      = ct+1;
  end
end

%% extract relevant pixels
ROIpercent  = zeros(numel(ROI),1);
for iROI = 1:numel(ROI)
  temp      = arrayfun(@(x,y)(squeeze(modmat(x,y))),ROI{iROI}(:,1),ROI{iROI}(:,2),'UniformOutput',false);
  temp      = cell2mat(temp');
  temp      = temp(:);
  ROIpercent(iROI) = sum(temp==1)./sum(~isnan(temp)) * 100;
end
end
