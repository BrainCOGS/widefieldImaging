function corrSumm = summarizeROIcorr(mice,subtractFlag)

% corrSumm = summarizeROIcorr(mice,subtractFlag)
% compiles correlation analyses from spontaneous and task recs
% mice is cell array with mouse names
% subtFlag to load analyses where grand ROI activity averages were
% subtracted (default false)

%% ------------------------------------------------------------------------

if nargin < 1 || isempty(mice); mice  = {'ai2';'ai3';'ai5';'ai7';'ai9';'ai10'}; end
if nargin < 2; subtractFlag = false; end

try
%% flag true to load from local disk instead of bucket
% (will have no effect if running on spock)
localFlag    = false; 

%% configuration for ROI corr analysis
cfg.ROIflag  = true;  % use ROIs
cfg.fisherZ  = true;  % for connectivity measures, transform cc mat
cfg.thCC     = 1;     % threshold z scored cc at this
cfg.binarize = false; % don't binarize
cfg.whichGLM = 'dffGLM_space_ridge_ROI.mat'; % for "signal corr"
cfg.mice     = mice;
wf           = widefieldParams;

%% ------------------------------------------------------------------------
%% compile spontaneous, t4, t11 for each mouse (session)
%% ------------------------------------------------------------------------

corrSumm.ROIlbl               = {};
corrSumm.allrecs_perf         = [];
corrSumm.allrecs_avgCCaccumul = [];

% loop through recs
for iMouse  = 1:numel(mice)
  %% spontaneous activity
  spontRec                           = wf.getMouseRecs(mice{iMouse},'spont',localFlag);
  corrSumm.mouse(iMouse).recls_spont = spontRec;
  spontCorr                          = dffCorr(spontRec,[],{},cfg,subtractFlag);
  ROIlbl                             = spontCorr.ROIlbl;
  if isempty(corrSumm.ROIlbl)
    corrSumm.ROIlbl = ROIlbl;
  else
    if sum(strcmpi(corrSumm.ROIlbl,ROIlbl)) ~= numel(ROIlbl)
      error('found incompatible ROI list')
    end
  end
  
  corrSumm.cc_spont_overall(:,:,iMouse)           = spontCorr.overall.cc;
  corrSumm.directionIndex_spont_overall(:,:,iMouse)     = spontCorr.overall.directionIndex;
  corrSumm.conn_avgROIcc_spont_overall(iMouse,:)  = spontCorr.overall.conn_avgROIcc;

  
  %% task (visual guide and accumulation)
  % average across recs for each animal first, but keep track of x-rec
  % corrs
  taskRecs                          = wf.getMouseRecs(mice{iMouse},'task',localFlag);
  corrSumm.mouse(iMouse).recls_task = taskRecs;
  
  for iRec = 1:numel(taskRecs)
    
    % calculate overall performance
    cd(taskRecs{iRec})
    load behavLog logSumm
    % include all trials (even below threshold)
    logSumm.currMaze(logSumm.currMaze == 12) = 4;
    trialidx = ~isBadTrial(logSumm,[],false,0) & logSumm.currMaze == max(logSumm.currMaze); 
    corrSumm.mouse(iMouse).perfOverall(iRec) = sum(logSumm.choice(trialidx) == logSumm.trialType(trialidx))/sum(trialidx);
    corrSumm.allrecs_perf(end+1,:)           = corrSumm.mouse(iMouse).perfOverall(iRec);
    
    % compile correlation measures
    taskCorr = dffCorr(taskRecs{iRec},[],{},cfg,subtractFlag);
       
    corrSumm.mouse(iMouse).cc_visGuide_overall(:,:,iRec)              = taskCorr.visGuide.overall.cc;
    corrSumm.mouse(iMouse).directionIndex_visGuide_overall(:,:,iRec)  = taskCorr.visGuide.overall.directionIndex;
    corrSumm.mouse(iMouse).conn_degrees_visGuide_overall(iRec,:)      = taskCorr.visGuide.overall.conn_degrees;
    corrSumm.mouse(iMouse).conn_avgROIcc_visGuide_overall(iRec,:)     = taskCorr.visGuide.overall.conn_avgROIcc;
    
    corrSumm.mouse(iMouse).cc_visGuide_wholeTrial_correct(:,:,iRec)             = taskCorr.visGuide.correct.wholeTrial.cc;
    corrSumm.mouse(iMouse).directionIndex_visGuide_wholeTrial_correct(:,:,iRec) = taskCorr.visGuide.correct.wholeTrial.directionIndex;
    corrSumm.mouse(iMouse).conn_degrees_visGuide_wholeTrial_correct(iRec,:)     = taskCorr.visGuide.correct.wholeTrial.conn_degrees;
    corrSumm.mouse(iMouse).conn_avgROIcc_visGuide_wholeTrial_correct(iRec,:)    = taskCorr.visGuide.correct.wholeTrial.conn_avgROIcc;
    
    
    if isfield(taskCorr,'glm')
      corrSumm.mouse(iMouse).cc_glm(:,:,iRec)                      = taskCorr.glm.cc;
    else
      corrSumm.mouse(iMouse).cc_glm(:,:,iRec)                      = nan(size(taskCorr.visGuide.overall.cc));
    end
    
    cc                                     = taskCorr.accumul.overall.cc;
    cc(logical(eye(size(cc))))             = nan;
    corrSumm.allrecs_avgCCaccumul(end+1,:) = nanmean(cc(:));
    
    corrSumm.mouse(iMouse).cc_accumul_overall(:,:,iRec)             = taskCorr.accumul.overall.cc;
    corrSumm.mouse(iMouse).directionIndex_accumul_overall(:,:,iRec) = taskCorr.accumul.overall.directionIndex;
    corrSumm.mouse(iMouse).conn_degrees_accumul_overall(iRec,:)     = taskCorr.accumul.overall.conn_degrees;
    corrSumm.mouse(iMouse).conn_avgROIcc_accumul_overall(iRec,:)    = taskCorr.accumul.overall.conn_avgROIcc;
    
    corrSumm.mouse(iMouse).cc_accumul_wholeTrial_correct(:,:,iRec)             = taskCorr.accumul.correct.wholeTrial.cc;
    corrSumm.mouse(iMouse).directionIndex_accumul_wholeTrial_correct(:,:,iRec) = taskCorr.accumul.correct.wholeTrial.directionIndex;
    corrSumm.mouse(iMouse).conn_degrees_accumul_wholeTrial_correct(iRec,:)     = taskCorr.accumul.correct.wholeTrial.conn_degrees;
    corrSumm.mouse(iMouse).conn_avgROIcc_accumul_wholeTrial_correct(iRec,:)    = taskCorr.accumul.correct.wholeTrial.conn_avgROIcc;
    
    corrSumm.mouse(iMouse).cc_accumul_wholeTrial_error(:,:,iRec)             = taskCorr.accumul.error.wholeTrial.cc;
    corrSumm.mouse(iMouse).directionIndex_accumul_wholeTrial_error(:,:,iRec) = taskCorr.accumul.error.wholeTrial.directionIndex;
    corrSumm.mouse(iMouse).conn_degrees_accumul_wholeTrial_error(iRec,:)     = taskCorr.accumul.error.wholeTrial.conn_degrees;
    corrSumm.mouse(iMouse).conn_avgROIcc_accumul_wholeTrial_error(iRec,:)    = taskCorr.accumul.error.wholeTrial.conn_avgROIcc;
    
    corrSumm.mouse(iMouse).cc_accumul_wholeTrial_easy(:,:,iRec)              = taskCorr.accumul.easy.wholeTrial.cc;
    corrSumm.mouse(iMouse).directionIndex_accumul_wholeTrial_easy(:,:,iRec)  = taskCorr.accumul.easy.wholeTrial.directionIndex;
    corrSumm.mouse(iMouse).conn_degrees_accumul_wholeTrial_easy(iRec,:)      = taskCorr.accumul.easy.wholeTrial.conn_degrees;
    corrSumm.mouse(iMouse).conn_avgROIcc_accumul_wholeTrial_easy(iRec,:)     = taskCorr.accumul.easy.wholeTrial.conn_avgROIcc;
    
    corrSumm.mouse(iMouse).cc_accumul_wholeTrial_hard(:,:,iRec)              = taskCorr.accumul.hard.wholeTrial.cc;
    corrSumm.mouse(iMouse).directionIndex_accumul_wholeTrial_hard(:,:,iRec)  = taskCorr.accumul.hard.wholeTrial.directionIndex;
    corrSumm.mouse(iMouse).conn_degrees_accumul_wholeTrial_hard(iRec,:)      = taskCorr.accumul.hard.wholeTrial.conn_degrees;
    corrSumm.mouse(iMouse).conn_avgROIcc_accumul_wholeTrial_hard(iRec,:)     = taskCorr.accumul.hard.wholeTrial.conn_avgROIcc;

    
    % some recs dont have distractor trials
    try
      corrSumm.mouse(iMouse).cc_accumul_wholeTrial_correct_nodistractor(:,:,iRec)             = taskCorr.accumul.correct_nodistractor.wholeTrial.cc;
      corrSumm.mouse(iMouse).directionIndex_accumul_wholeTrial_correct_nodistractor(:,:,iRec) = taskCorr.accumul.correct_nodistractor.wholeTrial.directionIndex;
      corrSumm.mouse(iMouse).conn_degrees_accumul_wholeTrial_correct_nodistractor(iRec,:)     = taskCorr.accumul.correct_nodistractor.wholeTrial.conn_degrees;
      corrSumm.mouse(iMouse).conn_avgROIcc_accumul_wholeTrial_correct_nodistractor(iRec,:)    = taskCorr.accumul.correct_nodistractor.wholeTrial.conn_avgROIcc;
    catch
      nROI = numel(corrSumm.ROIlbl);
      corrSumm.mouse(iMouse).cc_accumul_wholeTrial_correct_nodistractor(:,:,iRec)             = nan(nROI,nROI);
      corrSumm.mouse(iMouse).directionIndex_accumul_wholeTrial_correct_nodistractor(:,:,iRec) = nan(nROI,nROI);
      corrSumm.mouse(iMouse).conn_degrees_accumul_wholeTrial_correct_nodistractor(iRec,:)     = nan(1,nROI);
      corrSumm.mouse(iMouse).conn_avgROIcc_accumul_wholeTrial_correct_nodistractor(iRec,:)    = nan(1,nROI);
    end
    
    corrSumm.mouse(iMouse).cc_accumul_wholeTrial_correct_distractor(:,:,iRec)             = taskCorr.accumul.correct_distractor.wholeTrial.cc;
    corrSumm.mouse(iMouse).directionIndex_accumul_wholeTrial_correct_distractor(:,:,iRec) = taskCorr.accumul.correct_distractor.wholeTrial.directionIndex;
    corrSumm.mouse(iMouse).conn_degrees_accumul_wholeTrial_correct_distractor(iRec,:)     = taskCorr.accumul.correct_distractor.wholeTrial.conn_degrees;
    corrSumm.mouse(iMouse).conn_avgROIcc_accumul_wholeTrial_correct_distractor(iRec,:)    = taskCorr.accumul.correct_distractor.wholeTrial.conn_avgROIcc;
    
    corrSumm.mouse(iMouse).cc_accumul_cueHalf1_correct(:,:,iRec)             = taskCorr.accumul.correct.cueHalf1.cc;
    corrSumm.mouse(iMouse).directionIndex_accumul_cueHalf1_correct(:,:,iRec) = taskCorr.accumul.correct.cueHalf1.directionIndex;
    corrSumm.mouse(iMouse).conn_degrees_accumul_cueHalf1_correct(iRec,:)     = taskCorr.accumul.correct.cueHalf1.conn_degrees;
    corrSumm.mouse(iMouse).conn_avgROIcc_accumul_cueHalf1_correct(iRec,:)    = taskCorr.accumul.correct.cueHalf1.conn_avgROIcc;
    
    corrSumm.mouse(iMouse).cc_accumul_cueHalf2_correct(:,:,iRec)             = taskCorr.accumul.correct.cueHalf2.cc;
    corrSumm.mouse(iMouse).directionIndex_accumul_cueHalf2_correct(:,:,iRec) = taskCorr.accumul.correct.cueHalf2.directionIndex;
    corrSumm.mouse(iMouse).conn_degrees_accumul_cueHalf2_correct(iRec,:)     = taskCorr.accumul.correct.cueHalf2.conn_degrees;
    corrSumm.mouse(iMouse).conn_avgROIcc_accumul_cueHalf2_correct(iRec,:)    = taskCorr.accumul.correct.cueHalf2.conn_avgROIcc;
    
    corrSumm.mouse(iMouse).cc_accumul_delay_correct(:,:,iRec)                = taskCorr.accumul.correct.mem.cc;
    corrSumm.mouse(iMouse).directionIndex_accumul_delay_correct(:,:,iRec)    = taskCorr.accumul.correct.mem.directionIndex;
    corrSumm.mouse(iMouse).conn_degrees_accumul_delay_correct(iRec,:)        = taskCorr.accumul.correct.mem.conn_degrees;
    corrSumm.mouse(iMouse).conn_avgROIcc_accumul_delay_correct(iRec,:)       = taskCorr.accumul.correct.mem.conn_avgROIcc;
    
    corrSumm.mouse(iMouse).cc_visGuide_cueHalf1_correct(:,:,iRec)             = taskCorr.visGuide.correct.cueHalf1.cc;
    corrSumm.mouse(iMouse).directionIndex_visGuide_cueHalf1_correct(:,:,iRec) = taskCorr.visGuide.correct.cueHalf1.directionIndex;
    corrSumm.mouse(iMouse).conn_degrees_visGuide_cueHalf1_correct(iRec,:)     = taskCorr.visGuide.correct.cueHalf1.conn_degrees;
    corrSumm.mouse(iMouse).conn_avgROIcc_visGuide_cueHalf1_correct(iRec,:)    = taskCorr.visGuide.correct.cueHalf1.conn_avgROIcc;
    
    corrSumm.mouse(iMouse).cc_visGuide_cueHalf2_correct(:,:,iRec)             = taskCorr.visGuide.correct.cueHalf2.cc;
    corrSumm.mouse(iMouse).directionIndex_visGuide_cueHalf2_correct(:,:,iRec) = taskCorr.visGuide.correct.cueHalf2.directionIndex;
    corrSumm.mouse(iMouse).conn_degrees_visGuide_cueHalf2_correct(iRec,:)     = taskCorr.visGuide.correct.cueHalf2.conn_degrees;
    corrSumm.mouse(iMouse).conn_avgROIcc_visGuide_cueHalf2_correct(iRec,:)    = taskCorr.visGuide.correct.cueHalf2.conn_avgROIcc;
    
    corrSumm.mouse(iMouse).cc_visGuide_delay_correct(:,:,iRec)                = taskCorr.visGuide.correct.mem.cc;
    corrSumm.mouse(iMouse).directionIndex_visGuide_delay_correct(:,:,iRec)    = taskCorr.visGuide.correct.mem.directionIndex;
    corrSumm.mouse(iMouse).conn_degrees_visGuide_delay_correct(iRec,:)        = taskCorr.visGuide.correct.mem.conn_degrees;
    corrSumm.mouse(iMouse).conn_avgROIcc_visGuide_delay_correct(iRec,:)       = taskCorr.visGuide.correct.mem.conn_avgROIcc;
  end
  
  %% calculate average and corr across recs for this mouse
  fls = fields(corrSumm.mouse(iMouse));
  for iF = 1:numel(fls)  
    if ~isempty(strfind(fls{iF},'mouse'))       || ~isempty(strfind(fls{iF},'ROIlbl')) || ...
       ~isempty(strfind(fls{iF},'spont'))       || ~isempty(strfind(fls{iF},'stats'))  || ...
       ~isempty(strfind(fls{iF},'perfOverall')) || ~isempty(strfind(fls{iF},'recls'))
      continue
      
    elseif (~isempty(strfind(fls{iF},'cc_')) && isempty(strfind(fls{iF},'avgROI'))) ...
        || ~isempty(strfind(fls{iF},'directionIndex_'))
      thismat                        = corrSumm.mouse(iMouse).(fls{iF});
      [nx,ny,nz]                     = size(thismat);
      corrSumm.(fls{iF})(:,:,iMouse) = nanmean(thismat,3);
      
      thismat(repmat(logical(eye(nx,ny)),[1 1 nz]))       = nan; % i = j should be nan for this
      thismat                                             = reshape(thismat,[nx*ny nz]);
      thismat(isnan(sum(thismat,2)),:)                    = [];
      if isempty(thismat)
        corrSumm.stats.xRecCorr.pairs.(fls{iF}){:,iMouse} = nan;
        corrSumm.stats.xRecCorr.(fls{iF})(iMouse,:)       = nan;
        continue
      end
      thiscorr                                            = corr(thismat);
      thiscorr(triu(true(size(thiscorr)),0))              = nan; % since corr is symmetrical dont double dip ROI pair
      thiscorr                                            = thiscorr(:);
      corrSumm.stats.xRecCorr.pairs.(fls{iF}){:,iMouse}   = thiscorr;
      corrSumm.stats.xRecCorr.(fls{iF})(iMouse,:)         = nanmean(thiscorr(:));
    else
      corrSumm.(fls{iF})(iMouse,:)            = nanmean(corrSumm.mouse(iMouse).(fls{iF}));
      corrSumm.([fls{iF} '_sem'])(iMouse,:)   = nanstd(corrSumm.mouse(iMouse).(fls{iF}))./sqrt(iRec-1);

    end
  end
  
end

%% ------------------------------------------------------------------------
%% stats
%% ------------------------------------------------------------------------

%% x-rec correlations and relatinship to performance, just overall accumul cc
allcc = [];
for iMouse = 1:numel(corrSumm.mouse)
  for iRec = 1:size(corrSumm.mouse(iMouse).cc_accumul_overall,3)
    thiscc                             = corrSumm.mouse(iMouse).cc_accumul_overall(:,:,iRec);
    thiscc(triu(true(size(thiscc)),0)) = nan;
    thiscc                             = thiscc(:);
    allcc(:,end+1)                     = thiscc(~isnan(thiscc));
  end
end
allcc                                     = corr(allcc); % rec vs rec corr
allcc(triu(true(size(thiscorr)),0))       = nan;
allcc                                     = allcc(:);
corrSumm.stats.cc_accumul_overall_allrecs = allcc(~isnan(allcc));

[corrSumm.stats.cc_perfVSccAccumul,corrSumm.stats.p_perfVSccAccumul] ...
       = corr(corrSumm.allrecs_avgCCaccumul,corrSumm.allrecs_perf);

%% x-animal correlations (ie consistency of finding)
% calculate average and corr across recs for this mouse
fls = fields(corrSumm);
for iF = 1:numel(fls)
  if ~isempty(strfind(fls{iF},'mouse')) || ~isempty(strfind(fls{iF},'ROIlbl')) ||...
     ~isempty(strfind(fls{iF},'stats')) || ~isempty(strfind(fls{iF},'allrecs'))
    continue

  elseif (~isempty(strfind(fls{iF},'cc_')) && isempty(strfind(fls{iF},'avgROI'))) ...
        || ~isempty(strfind(fls{iF},'directionIndex_'))
    thismat                        = corrSumm.(fls{iF});
    [nx,ny,nz]                     = size(thismat);

    thismat(repmat(logical(eye(nx,ny)),[1 1 nz]))       = nan; % i = j should be nan for this
    thismat                                             = reshape(thismat,[nx*ny nz]);

  else
    thismat                        = corrSumm.(fls{iF})';

  end

  thismat(isnan(sum(thismat,2)),:)                    = [];
  if isempty(thismat)
    corrSumm.stats.xMouseCorr.pairs.(fls{iF})         = nan;
    corrSumm.stats.xMouserr.(fls{iF})                 = nan;
    continue
  end
  thiscorr                                            = corr(thismat);
  thiscorr(triu(true(size(thiscorr)),0))              = nan; % since corr is symmetrical dont double dip ROI pair
  thiscorr                                            = thiscorr(:);
  corrSumm.stats.xMouseCorr.pairs.(fls{iF})           = thiscorr;
  corrSumm.stats.xMouseCorr.(fls{iF})                 = nanmean(thiscorr(:));

end

%% tests, t4 vs t11, correct vs error etc
testwhat = {'conn_avgROIcc','conn_degrees'};
compwhat = {{'visGuide_overall','accumul_overall'},                {'visGuide_overall','spont_overall'}                       ...
            {'spont_overall','accumul_overall'},                   {'accumul_wholeTrial_correct','accumul_wholeTrial_error'}  ...
            {'accumul_wholeTrial_easy','accumul_wholeTrial_hard'}, {'accumul_cueHalf1_correct','accumul_cueHalf2_correct'}    ...
            {'accumul_cueHalf1_correct','accumul_delay_correct'},  {'accumul_cueHalf2_correct','accumul_delay_correct'}       ...
            {'accumul_wholeTrial_correct','visGuide_wholeTrial_correct'} , {'accumul_wholeTrial_correct_nodistractor','visGuide_wholeTrial_correct_nodistractor'}  ...
            };

for iComp = 1:numel(compwhat)
  for iTest = 1:numel(testwhat)
    c1  = corrSumm.([testwhat{iTest} '_' compwhat{iComp}{1}]);
    c2  = corrSumm.([testwhat{iTest} '_' compwhat{iComp}{2}]);
    try
    c1s = corrSumm.([testwhat{iTest} '_' compwhat{iComp}{1} '_sem']);
    c2s = corrSumm.([testwhat{iTest} '_' compwhat{iComp}{2} '_sem']);
    catch
      c1s = nan(size(c1)); c2s = c1s;
    end

    if size(c1,2) > 1 % for ROI avgs, test both across mice and across ROIs
      c1m = mean(c1,2);
      c2m = mean(c2,2);
      c1r = mean(c1,1)';
      c2r = mean(c2,1)';
      c1ms = std(c1,0,2)./sqrt(size(c1,2)-1);
      c2ms = std(c2,0,2)./sqrt(size(c1,2)-1);
      c1rs = std(c1)'./sqrt(size(c1,1)-1);
      c2rs = std(c2)'./sqrt(size(c1,1)-1);

      corrSumm.stats.xmouse.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).avgs     ...
                = [c1m c2m];
      corrSumm.stats.xmouse.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).sems     ...
                = [c1ms  c2ms];
      corrSumm.stats.xmouse.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).lbls     ...
                = compwhat{iComp};

      corrSumm.stats.xroi.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).avgs     ...
                = [c1r c2r];
      corrSumm.stats.xroi.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).sems     ...
                = [c1rs  c2rs];
      corrSumm.stats.xroi.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).lbls     ...
                = compwhat{iComp};

      if lillietest(c1m) || lillietest(c2m)
        corrSumm.stats.xmouse.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).pval     ...
                = signrank(c1m,c2m);
        corrSumm.stats.xmouse.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).test     ...
                = 'signrank';   
      else
        [~,corrSumm.stats.xmouse.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).pval] ...
                = ttest(c1m,c2m);
        corrSumm.stats.xmouse.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).test     ...
                = 'ttest';
      end
      if lillietest(c1r) || lillietest(c2r)
        corrSumm.stats.xroi.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).pval     ...
                = signrank(c1r,c2r);
        corrSumm.stats.xroi.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).test     ...
                = 'signrank';
      else
        [~,corrSumm.stats.xroi.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).pval] ...
                = ttest(c1r,c2r);
        corrSumm.stats.xroi.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).test     ...
                = 'ttest';
      end

    else
      corrSumm.stats.xmouse.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).avgs     ...
                = [c1  c2];
      corrSumm.stats.xmouse.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).sems     ...
                = [c1s  c2s];
      corrSumm.stats.xmouse.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).lbls     ...
                = compwhat{iComp};

      if lillietest(c1) || lillietest(c2)
        corrSumm.stats.xmouse.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).pval     ...
                = signrank(c1,c2);
        corrSumm.stats.xmouse.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).test     ...
                = 'signrank';
      else
        [~,corrSumm.stats.xmouse.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).pval] ...
                = ttest(c1,c2);
        corrSumm.stats.xmouse.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).test     ...
                = 'ttest';
      end
    end
  end
end

%% ------------------------------------------------------------------------
%% save and plot
%% ------------------------------------------------------------------------
rootdir      = wf.getRootDir(isThisSpock,localFlag);
corrSumm.cfg = cfg;
cd(rootdir)
if subtractFlag
  save ROIcorrSummary_subt corrSumm cfg
  if ~isempty(dir('ROIcorrSummary_subt.pdf')); delete ROIcorrSummary_subt.pdf; end
else
  save ROIcorrSummary corrSumm cfg
  if ~isempty(dir('ROIcorrSummary.pdf')); delete ROIcorrSummary.pdf; end
end
%% plot overall cc: spont, visual guide, accumulation
plots = {'spont_overall','visGuide_overall','accumul_overall'};
fh    = figure;
wf.applyFigDefaults(fh,[numel(plots)+1 1],'w')

for iPlot = 1:numel(plots)
  cc   = mean(corrSumm.(['cc_' plots{iPlot}]),3);
  nROI = size(cc,1);
  subplot(1,numel(plots),iPlot)
  imagesc(cc,[0 1]); colormap red2blue
  set(gca,'xtick',1:nROI,'ytick',1:nROI,'xticklabel',corrSumm.ROIlbl,'yticklabel',corrSumm.ROIlbl)
  rotateXLabels(gca,90)
  wf.applyAxisDefaults(gca,'k'); axis tight
  wf.applyAxisLbls(gca,[],[],removeUnderscores(plots{iPlot}))

  if iPlot == numel(plots)
    ch = colorbar('location','eastoutside','position',[.92 .18 .015 .15],'Color','k');
    ch.Label.String = 'r';
  end
end

%% plot cc: correct vs error, easy vs hard, visual guide vs accumul (each + subtraction)
clear plots
plots{1} = {'accumul_wholeTrial_correct','visGuide_wholeTrial_correct'};
plots{2} = {'accumul_wholeTrial_hard','accumul_wholeTrial_easy'};
plots{3} = {'accumul_wholeTrial_correct','accumul_wholeTrial_error'};
plots{3} = {'accumul_wholeTrial_correct_nodistractor','visGuide_wholeTrial_correct_nodistractor'};

fh = figure;
wf.applyFigDefaults(fh,[3 numel(plots)+1],'w')
for iPlot = 1:numel(plots)
  cc1  = mean(corrSumm.(['cc_' plots{iPlot}{1}]),3);
  cc2  = mean(corrSumm.(['cc_' plots{iPlot}{2}]),3);
  ccs  = cc1 - cc2;
  nROI = size(cc1,1);
  clim = [min([cc1(:); cc2(:)]) max([cc1(:); cc2(:)])];
  
  subplot(numel(plots),3,(iPlot-1)*3+1)
  imagesc(cc1,clim); colormap red2blue
  set(gca,'xtick',1:nROI,'ytick',1:nROI,'xticklabel',corrSumm.ROIlbl,'yticklabel',corrSumm.ROIlbl)
  rotateXLabels(gca,90)
  wf.applyAxisDefaults(gca,'k'); axis tight
  wf.applyAxisLbls(gca,[],[],removeUnderscores(plots{iPlot}{1}))
  
  subplot(numel(plots),3,(iPlot-1)*3+2)
  imagesc(cc2,clim); colormap red2blue
  set(gca,'xtick',1:nROI,'ytick',1:nROI,'xticklabel',corrSumm.ROIlbl,'yticklabel',corrSumm.ROIlbl)
  rotateXLabels(gca,90)
  wf.applyAxisDefaults(gca,'k'); axis tight
  wf.applyAxisLbls(gca,[],[],removeUnderscores(plots{iPlot}{2}))
  
  ch1 = colorbar('location','eastoutside','position',[.63 .7-.3*(iPlot-1) .01 .08],'Color','k');
  ch1.Label.String = 'r';
  
  subplot(numel(plots),3,(iPlot-1)*3+3)
  imagesc(ccs,[-max(abs(ccs(:))) max(abs(ccs(:)))]); colormap red2blue
  set(gca,'xtick',[],'ytick',[],'xticklabel',corrSumm.ROIlbl,'yticklabel',corrSumm.ROIlbl)
  rotateXLabels(gca,90)
  wf.applyAxisDefaults(gca,'k'); axis tight
  wf.applyAxisLbls(gca,[],[],'left - middle')
  
  ch2 = colorbar('location','eastoutside','position',[.93 .7-.3*(iPlot-1) .01 .08],'Color','k');
  ch2.Label.String = 'r';
end

%% plot cc, epochs (each + subtractions)
cc1   = mean(corrSumm.cc_accumul_cueHalf1_correct,3);
cc2   = mean(corrSumm.cc_accumul_cueHalf2_correct,3);
cc3   = mean(corrSumm.cc_accumul_delay_correct,3);
cc12  = cc1-cc2;
cc13  = cc1-cc3;
cc23  = cc2-cc3;
clim  = [min([cc1(:);cc2(:);cc3(:)]) max([cc1(:);cc2(:);cc3(:)])];
clims = [-max(abs([cc12(:);cc13(:);cc23(:)])) max(abs([cc12(:);cc13(:);cc23(:)]))];

fh = figure;
wf.applyFigDefaults(fh,[3 2],'w')

subplot(2,3,1)
imagesc(cc1,clim); colormap red2blue
set(gca,'xtick',1:nROI,'ytick',1:nROI,'xticklabel',corrSumm.ROIlbl,'yticklabel',corrSumm.ROIlbl)
rotateXLabels(gca,90)
wf.applyAxisDefaults(gca,'k'); axis tight
wf.applyAxisLbls(gca,[],[],'cueHalf1')

subplot(2,3,2)
imagesc(cc2,clim); colormap red2blue
set(gca,'xtick',1:nROI,'ytick',1:nROI,'xticklabel',corrSumm.ROIlbl,'yticklabel',corrSumm.ROIlbl)
rotateXLabels(gca,90)
wf.applyAxisDefaults(gca,'k'); axis tight
wf.applyAxisLbls(gca,[],[],'cueHalf2')

subplot(2,3,3)
imagesc(cc3,clim); colormap red2blue
set(gca,'xtick',1:nROI,'ytick',1:nROI,'xticklabel',corrSumm.ROIlbl,'yticklabel',corrSumm.ROIlbl)
rotateXLabels(gca,90)
wf.applyAxisDefaults(gca,'k'); axis tight
wf.applyAxisLbls(gca,[],[],'delay')

ch1 = colorbar('location','eastoutside','position',[.92 .6 .01 .08],'Color','k');
ch1.Label.String = 'r';

subplot(2,3,4)
imagesc(cc12,clims); colormap red2blue
set(gca,'xtick',1:nROI,'ytick',1:nROI,'xticklabel',corrSumm.ROIlbl,'yticklabel',corrSumm.ROIlbl)
rotateXLabels(gca,90)
wf.applyAxisDefaults(gca,'k'); axis tight
wf.applyAxisLbls(gca,[],[],'cueHalf1 - cueHalf2')

subplot(2,3,5)
imagesc(cc13,clims); colormap red2blue
set(gca,'xtick',1:nROI,'ytick',1:nROI,'xticklabel',corrSumm.ROIlbl,'yticklabel',corrSumm.ROIlbl)
rotateXLabels(gca,90)
wf.applyAxisDefaults(gca,'k'); axis tight
wf.applyAxisLbls(gca,[],[],'cueHalf1 - delay')

subplot(2,3,6)
imagesc(cc23,clims); colormap red2blue
set(gca,'xtick',1:nROI,'ytick',1:nROI,'xticklabel',corrSumm.ROIlbl,'yticklabel',corrSumm.ROIlbl)
rotateXLabels(gca,90)
wf.applyAxisDefaults(gca,'k'); axis tight
wf.applyAxisLbls(gca,[],[],'cueHalf2 - delay')

ch1 = colorbar('location','eastoutside','position',[.92 .1 .01 .08],'Color','k');
ch1.Label.String = 'r';

%% connectivity measures scatter plots
testwhat = {'conn_avgROIcc','conn_degrees'};
compwhat = {{'visGuide_overall','accumul_overall'},                {'visGuide_overall','spont_overall'}                       ...
            {'spont_overall','accumul_overall'},                   {'accumul_wholeTrial_correct','accumul_wholeTrial_error'}  ...
            {'accumul_wholeTrial_easy','accumul_wholeTrial_hard'}, {'accumul_cueHalf1_correct','accumul_cueHalf2_correct'}    ...
            {'accumul_cueHalf1_correct','accumul_delay_correct'},  {'accumul_cueHalf2_correct','accumul_delay_correct'}       ...
            {'accumul_wholeTrial_correct','visGuide_wholeTrial_correct'}  ...
            };

for iTest = 1:numel(testwhat)
  
  % plot crosses for mouse data (avg+/- sem)
  [nr,nc] = subplotOrg(numel(compwhat),4);
  fh      = figure;
  wf.applyFigDefaults(fh,[nc nr],'w')

  for iComp = 1:numel(compwhat)
    subplot(nr,nc,iComp); hold on
    avg  = corrSumm.stats.xmouse.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).avgs;
    sem  = corrSumm.stats.xmouse.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).sems;
    pval = corrSumm.stats.xmouse.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).pval;
    lbls = corrSumm.stats.xmouse.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).lbls;

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
    wf.applyAxisLbls(gca,removeUnderscores(lbls{1}),removeUnderscores(lbls{2}),removeUnderscores(testwhat{iTest}))
  end
  
  switch testwhat{iTest}
    % plot crosses for roi data (avg+/- sem)
    case {'conn_avgROIcc','conn_degrees'}
      [nr,nc] = subplotOrg(numel(compwhat),4);
      fh      = figure;
      wf.applyFigDefaults(fh,[nc nr],'w')

      for iComp = 1:numel(compwhat)
        subplot(nr,nc,iComp); hold on
        avg  = corrSumm.stats.xroi.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).avgs;
        sem  = corrSumm.stats.xroi.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).sems;
        pval = corrSumm.stats.xroi.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).pval;
        lbls = corrSumm.stats.xroi.(testwhat{iTest}).([compwhat{iComp}{1} '_' compwhat{iComp}{2}]).lbls;
        
        plotROIlbl  = {'V1','mV2','PPC','RSC','aM2','mM2','M1','SS'};
        plotROI     = {'VISp','mV2','VISa','RSP','aMOs','mMOs','MOp','SS'};
      
        for iROI = 1:numel(plotROI)
          idx = find(strcmpi(corrSumm.ROIlbl,[plotROI{iROI} '-R']) | strcmpi(corrSumm.ROIlbl,[plotROI{iROI} '-L']));
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
        wf.applyAxisLbls(gca,removeUnderscores(lbls{1}),removeUnderscores(lbls{2}),removeUnderscores(testwhat{iTest}))
    
      end
  end
end

%% overall x-rec corr and corr between corr and performance
fh = figure;
wf.applyFigDefaults(fh,[2 1],'w')

subplot(1,2,1)
h1 = histogram(corrSumm.stats.cc_accumul_overall_allrecs,0:.1:1);
h1.EdgeColor = wf.mediumgray;
wf.applyAxisDefaults(gca,'k'); xlim([0 1])
wf.applyAxisLbls(gca,'r (accumul overall)','rec pairs')

subplot(1,2,2)
plot(corrSumm.allrecs_perf,corrSumm.allrecs_avgCCaccumul,'o','color',wf.mediumgray)
xlim([.2 1]); ylim([.5 1])
text(.3,.9,sprintf('r = %1.2f\np = %1.2g',corrSumm.stats.cc_perfVSccAccumul,corrSumm.stats.p_perfVSccAccumul))
wf.applyAxisDefaults(gca,'k'); 
wf.applyAxisLbls(gca,'performance','<r (accumul overall)>')

%% export
figls = get(0,'children'); % print to pdf
for ii = length(figls):-1:1
  figure(figls(ii))
  if subtractFlag
    export_fig ROIcorrSummary_subt.pdf -q101 -append
  else
    export_fig ROIcorrSummary.pdf -q101 -append
  end
end
close all

catch ME
  displayException(ME)
end