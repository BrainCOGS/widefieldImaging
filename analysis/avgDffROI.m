function avgROI = avgDffROI(dff,logSumm,doZScore)

% avgDff = avgDffByTrialEpoch(logSumm,dff)
% logSumm is output by summarizeVirmenLog_widefield, dff is full frames x ROIs matrix
% calculates average +/- sem dff for combinations of mazes, trial types and
% trial epochs. Look at function body for details
% LP july 2016

%%
if nargin < 3; doZScore = false; end

%%
tic
fprintf('calculating average dffs')

%% load if available and return
if ~isempty(dir('avgROI.mat'))
  warning('off','all')
  load avgROI avgROI
  warning('on','all')
  if exist('avgDff','var')
    fprintf(' done after %1.1f min (loaded from disk)\n',toc/60)
    return
  end
end

%% some recs in blocks condensed have maze 12, this is just like maze 4 but with towers on both sides
% here they will be analyzed together for convenience
if contains(char(logSumm.info.protocol),'Condensed')
  logSumm.currMaze(logSumm.currMaze == 12) = 4;
end

%% analysis parameters
% list type, presence of distractors, maze
avgROI.trialTypeList      = {'correct','error','correctR','correctL','errorR','errorL',               ...
  'correct_distractor','correctR_distractor','correctL_distractor',        ...
  'errorR_distractor','errorL_distractor','correct_nodistractor',          ...
  'correctR_nodistractor','correctL_nodistractor','error_nodistractor',    ...
  'errorR_nodistractor','errorL_nodistractor','correct_hard','error_hard', ...
  'correctR_hard','correctL_hard','errorR_hard','errorL_hard',             ...
  'correct_easy','error_easy','correctR_easy','correctL_easy',             ...
  'errorR_easy','errorL_easy','error_distractor'};

avgROI.mazeList           = unique(logSumm.currMaze);
avgROI.epochList          = {'start','cueHalf1','cueHalf2','mem','choice','rw','iti'};
avgROI.towers.RminusLBins = -15:5:15;
avgROI.towers.RminusLvals = toBinCenters(avgROI.towers.RminusLBins);
avgROI.towers.ntowersBins = 0:3:15;
avgROI.towers.ntowers     = toBinCenters(avgROI.towers.ntowersBins);
avgROI.fps                = 1/logSumm.frameDtCam;
avgROI.nsec               = 8;
avgROI.nsecTower          = 4;
avgROI.nsecTowerResp      = .4;
avgROI.towerRespNFrames   = round(avgROI.nsecTowerResp * avgROI.fps);
avgROI.spaceBins          = 0:5:300;
avgROI.isZscored          = doZScore;

%% z-score
if doZScore; dff = zscore(dff); end

%% initialize matrices
[nZ,nROI]                 = size(dff);
tempdff_tR                = [];
tempdff_tL                = [];
tempdff_tR_fr             = zeros(floor(avgROI.nsecTower*avgROI.fps)+1,nROI);
tempdff_tL_fr             = zeros(floor(avgROI.nsecTower*avgROI.fps)+1,nROI);

for ii = 1:numel(avgROI.towers.ntowers)
  tempdff_nt(ii).R        = [];
  tempdff_nt(ii).L        = [];
end
for ii = 1:numel(avgROI.towers.RminusLvals)
  tempdff_dt(ii).RmL      = [];
end

%% loop over trial types, mazes and epochs to compile average dff for
% combinations thereof
for iType = 1:numel(avgROI.trialTypeList)
  fprintf('.')
  for iMaze = 1:numel(avgROI.mazeList)
    
    %% get list of relevant trials
    trialidx   = getTrialIdx(logSumm,avgROI.trialTypeList{iType},avgROI.mazeList(iMaze),widefieldParams.perfTh);
    
    if isempty(trialidx)
      avgROI.maze(iMaze).(avgROI.trialTypeList{iType}).space.frames = [];
      avgROI.maze(iMaze).(avgROI.trialTypeList{iType}).space.avg    = [];
      avgROI.maze(iMaze).(avgROI.trialTypeList{iType}).space.sem    = [];
      for iEpoch = 1:numel(avgROI.epochList)
        avgROI.maze(iMaze).(avgROI.trialTypeList{iType}).(avgROI.epochList{iEpoch}).frames = [];
        avgROI.maze(iMaze).(avgROI.trialTypeList{iType}).(avgROI.epochList{iEpoch}).avg    = [];
        avgROI.maze(iMaze).(avgROI.trialTypeList{iType}).(avgROI.epochList{iEpoch}).sem    = [];
      end
      continue
    end
    
    %% do spatial averages
    tempfr  = zeros(numel(avgROI.spaceBins),nROI);
    for iTrial = 1:numel(trialidx)
      idx     = arrayfun(@(x)(find(logSumm.binned.pos{trialidx(iTrial)}(:,2) >= x,1,'first')),avgROI.spaceBins);
      fidx    = logSumm.binned.camFrameID{trialidx(iTrial)}(idx);
      tempfr  = tempfr + dff(fidx,:);% preserves frames from iTrial = 0 to 1 s
    end
    
    avgROI.maze(iMaze).(avgROI.trialTypeList{iType}).space.frames = tempfr ./ numel(trialidx);
    %% do time averages triggered by different epochs
    for iEpoch = 1:numel(avgROI.epochList)
      tempdff = [];
      tempfr  = zeros(round(avgROI.fps*avgROI.nsec),nROI);
      for iTrial = 1:numel(trialidx)
        % definitions of trial epochs
        switch avgROI.epochList{iEpoch}
          case 'start'
            fidx = logSumm.binned.camFrameID{trialidx(iTrial)}(logSumm.binned.pos{trialidx(iTrial)}(:,2) <= 0);
          case 'cueHalf1'
            fidx = logSumm.binned.camFrameID{trialidx(iTrial)}(logSumm.binned.pos{trialidx(iTrial)}(:,2) > 0 & logSumm.binned.pos{trialidx(iTrial)}(:,2) <= 100);
          case 'cueHalf2'
            fidx = logSumm.binned.camFrameID{trialidx(iTrial)}(logSumm.binned.pos{trialidx(iTrial)}(:,2) > 100 & logSumm.binned.pos{trialidx(iTrial)}(:,2) <= 200);
          case 'mem' % first 80 cm of mem period
            fidx = logSumm.binned.camFrameID{trialidx(iTrial)}(logSumm.binned.pos{trialidx(iTrial)}(:,2) > 200 & logSumm.binned.pos{trialidx(iTrial)}(:,2) <= 280);
          case 'choice' % last 20 cm of memory until reward, defined by rotational velocity later
            fidx = logSumm.binned.camFrameID{trialidx(iTrial)}(logSumm.binned.pos{trialidx(iTrial)}(:,2) > 280);
            fidx(fidx>=logSumm.binned.keyFrames{trialidx(iTrial)}(end)) = [];
          case 'rw' % 2s after reward
            fidx = logSumm.binned.keyFrames{trialidx(iTrial)}(end):logSumm.binned.keyFrames{trialidx(iTrial)}(end)+floor(2*avgROI.fps);
            fidx(fidx>logSumm.binned.camFrameID{trialidx(iTrial)}(end)) = [];
          case 'iti' % remaining period after 2s
            fidx = logSumm.binned.keyFrames{trialidx(iTrial)}(end)+floor(2*(1/logSumm.frameDtVirmen))+1:logSumm.binned.camFrameID{trialidx(iTrial)}(end);
        end
        tempdff     = [tempdff; nanmean(dff(fidx,:))]; % will be averaged across frames
        theseframes = fidx(1):fidx(1)+round(avgROI.fps*avgROI.nsec)-1;
        if theseframes(end) > nZ; continue; end
        tempfr      = tempfr + dff(theseframes,:);
      end
      
      avgROI.maze(iMaze).(avgROI.trialTypeList{iType}).(avgROI.epochList{iEpoch}).frames = tempfr ./ numel(trialidx);
      avgROI.maze(iMaze).(avgROI.trialTypeList{iType}).(avgROI.epochList{iEpoch}).avg    = nanmean(tempdff);
      avgROI.maze(iMaze).(avgROI.trialTypeList{iType}).(avgROI.epochList{iEpoch}).sem    = nanstd(tempdff)./sqrt(numel(trialidx)-1);
    end
    
    %% add tower-triggered averages
    if iMaze > 1 && (strcmpi(avgROI.trialTypeList{iType},'correctL') || strcmpi(avgROI.trialTypeList{iType},'correctR'))
      for iTrial = 1:numel(trialidx)
        
        % basline-subtracted response for each R and L tower
        nr  = numel(logSumm.cuePos_R{trialidx(iTrial)});
        nl  = numel(logSumm.cuePos_L{trialidx(iTrial)});
        rml = nr - nl;
        for rr = 1:nr
          fidx          = logSumm.binned.camFrameID{trialidx(iTrial)} ...
            (find(logSumm.binned.time{trialidx(iTrial)} >= logSumm.cueOnset_R{trialidx(iTrial)}(rr),1,'first'));
          fidx          = fidx:fidx+floor(avgROI.nsecTower*avgROI.fps);
          thisdff       = dff(fidx,:);
          basel         = repmat(mean(dff(fidx(1)-2:fidx(1)-1,:)),[size(thisdff,1) 1]);
          thisdff       = thisdff - basel;
          tempdff_tR    = [tempdff_tR; nanmean(dff(1:avgROI.towerRespNFrames,:))];
          tempdff_tR_fr = tempdff_tR_fr + thisdff;
        end
        for ll = 1:nl
          fidx          = logSumm.binned.camFrameID{trialidx(iTrial)} ...
            (find(logSumm.binned.time{trialidx(iTrial)} >= logSumm.cueOnset_L{trialidx(iTrial)}(ll),1,'first'));
          fidx          = fidx:fidx+floor(avgROI.nsecTower*avgROI.fps);
          thisdff       = dff(fidx,:);
          basel         = repmat(mean(dff(fidx(1)-2:fidx(1)-1,:)),[size(thisdff,1) 1]);
          thisdff       = thisdff - basel;
          tempdff_tL    = [tempdff_tL; nanmean(dff(1:avgROI.towerRespNFrames,:))];
          tempdff_tL_fr = tempdff_tL_fr + thisdff;
        end
        
        % binned by #towers or delta, averaged between 180 & 230
        if nr >= max(avgROI.towers.ntowersBins); nr = max(avgROI.towers.ntowersBins)-1; end
        if nl >= max(avgROI.towers.ntowersBins); nl = max(avgROI.towers.ntowersBins)-1; end
        if abs(rml) >= max(abs(avgROI.towers.RminusLBins)); rml = sign(rml)*(max(abs(avgROI.towers.RminusLBins))-1); end
        
        fidx = logSumm.binned.camFrameID{trialidx(iTrial)}(logSumm.binned.pos{trialidx(iTrial)}(:,2) > 180 & logSumm.binned.pos{trialidx(iTrial)}(:,2) <= 230);
        
        binidx = find(avgROI.towers.ntowersBins <= nr,1,'last');
        tempdff_nt(binidx).R      = [tempdff_nt(binidx).R; nanmean(dff(fidx,:))];
        
        binidx = find(avgROI.towers.ntowersBins <= nl,1,'last');
        tempdff_nt(binidx).L      = [tempdff_nt(binidx).L; nanmean(dff(fidx,:))];
        
        binidx = find(avgROI.towers.RminusLBins <= rml,1,'last');
        tempdff_dt(binidx).RmL    = [tempdff_dt(binidx).RmL; nanmean(dff(fidx,:))];
      end
    end
  end
  
end

%% averages and SEM for towers
avgROI.towersR.frames                  = tempdff_tR_fr ./ size(tempdff_tR,1);
avgROI.towersR.avg                     = nanmean(tempdff_tR);
avgROI.towersR.sem                     = nanstd(tempdff_tR)./sqrt(size(tempdff_tR,1)-1);
avgROI.towersL.frames                  = tempdff_tL_fr ./ size(tempdff_tL,1);
avgROI.towersL.avg                     = nanmean(tempdff_tL);
avgROI.towersL.sem                     = nanstd(tempdff_tL)./sqrt(size(tempdff_tL,1)-1);

for ii = 1:numel(avgROI.towers.ntowers)
  avgROI.towers.numtowers(ii).R.avg    = nanmean(tempdff_nt(ii).R);
  avgROI.towers.numtowers(ii).R.sem    = nanstd(tempdff_nt(ii).R)./sqrt(size(tempdff_nt(ii).R,1)-1);
  
  avgROI.towers.numtowers(ii).L.avg    = nanmean(tempdff_nt(ii).L);
  avgROI.towers.numtowers(ii).L.sem    = nanstd(tempdff_nt(ii).L)./sqrt(size(tempdff_nt(ii).L,1)-1);
end
for ii = 1:numel(avgROI.towers.RminusLvals)
  avgROI.towers.RminusL(ii).avg        = nanmean(tempdff_dt(ii).RmL);
  avgROI.towers.RminusL(ii).sem        = nanstd(tempdff_dt(ii).RmL)./sqrt(size(tempdff_dt(ii).RmL,1)-1);
end

%% save
save avgROI avgROI -v7.3
fprintf('\tdone after %1.1f min\n',toc/60)