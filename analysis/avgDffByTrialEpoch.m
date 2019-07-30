function avgDff = avgDffByTrialEpoch(logSumm,dff)

% avgDff = avgDffByTrialEpoch(logSumm,dff)
% logSumm is output by summarizeVirmenLog_widefield, dff is full frames x pixels matrix
% calculates average +/- sem dff for combinations of mazes, trial types and
% trial epochs. Look at function body for details
% LP july 2016

tic
fprintf('calculating average dffs')

%% load if available and return
if ~isempty(dir('avgDff.mat'))
  warning('off','all')
  load avgDff avgDff
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
avgDff.trialTypeList      = {'correct','error','correctR','correctL','errorR','errorL',               ...
  'correct_distractor','correctR_distractor','correctL_distractor',        ...
  'errorR_distractor','errorL_distractor','correct_nodistractor',          ...
  'correctR_nodistractor','correctL_nodistractor','error_nodistractor',    ...
  'errorR_nodistractor','errorL_nodistractor','correct_hard','error_hard', ...
  'correctR_hard','correctL_hard','errorR_hard','errorL_hard',             ...
  'correct_easy','error_easy','correctR_easy','correctL_easy',             ...
  'errorR_easy','errorL_easy','error_distractor'};

avgDff.mazeList           = unique(logSumm.currMaze);
avgDff.epochList          = {'start','cueHalf1','cueHalf2','mem','choice','rw','iti'};
avgDff.towers.RminusLBins = -15:5:15;
avgDff.towers.RminusLvals = toBinCenters(avgDff.towers.RminusLBins);
avgDff.towers.ntowersBins = 0:3:15;
avgDff.towers.ntowers     = toBinCenters(avgDff.towers.ntowersBins);
avgDff.fps                = 1/logSumm.frameDtCam;
avgDff.nsec               = 8;
avgDff.nsecTower          = 4;
avgDff.nsecTowerResp      = .4;
avgDff.towerRespNFrames   = round(avgDff.nsecTowerResp * avgDff.fps);
avgDff.spaceBins          = 0:5:300;

%% initialize matrices
[nX,nY,nZ]                = size(dff);
tempdff_tR                = [];
tempdff_tL                = [];
tempdff_tR_fr             = zeros(nX,nY,floor(avgDff.nsecTower*avgDff.fps)+1);
tempdff_tL_fr             = zeros(nX,nY,floor(avgDff.nsecTower*avgDff.fps)+1);

for ii = 1:numel(avgDff.towers.ntowers)
  tempdff_nt(ii).R        = [];
  tempdff_nt(ii).L        = [];
end
for ii = 1:numel(avgDff.towers.RminusLvals)
  tempdff_dt(ii).RmL      = [];
end

%% loop over trial types, mazes and epochs to compile average dff for
% combinations thereof
for iType = 1:numel(avgDff.trialTypeList)
  fprintf('.')
  for iMaze = 1:numel(avgDff.mazeList)
    
    %% get list of relevant trials
    trialidx   = getTrialIdx(logSumm,avgDff.trialTypeList{iType},avgDff.mazeList(iMaze),widefieldParams.perfTh);
    
    if isempty(trialidx)
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).space.frames = [];
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).space.avg    = [];
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).space.sem    = [];
      for iEpoch = 1:numel(avgDff.epochList)
        avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).(avgDff.epochList{iEpoch}).frames = [];
        avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).(avgDff.epochList{iEpoch}).avg    = [];
        avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).(avgDff.epochList{iEpoch}).sem    = [];
      end
      continue
    end
    
    %% do spatial averages
    tempfr  = zeros(nX,nY,numel(avgDff.spaceBins));
    for iTrial = 1:numel(trialidx)
      idx     = arrayfun(@(x)(find(logSumm.binned.pos{trialidx(iTrial)}(:,2) >= x,1,'first')),avgDff.spaceBins);
      fidx    = logSumm.binned.camFrameID{trialidx(iTrial)}(idx);
      tempfr  = tempfr + dff(:,:,fidx);% preserves frames from iTrial = 0 to 1 s
    end
    
    avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).space.frames = tempfr ./ numel(trialidx);
    %% do time averages triggered by different epochs
    for iEpoch = 1:numel(avgDff.epochList)
      tempdff = [];
      tempfr  = zeros(nX,nY,round(avgDff.fps*avgDff.nsec));
      for iTrial = 1:numel(trialidx)
        % definitions of trial epochs
        switch avgDff.epochList{iEpoch}
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
            fidx = logSumm.binned.keyFrames{trialidx(iTrial)}(end):logSumm.binned.keyFrames{trialidx(iTrial)}(end)+floor(2*avgDff.fps);
            fidx(fidx>logSumm.binned.camFrameID{trialidx(iTrial)}(end)) = [];
          case 'iti' % remaining period after 2s
            fidx = logSumm.binned.keyFrames{trialidx(iTrial)}(end)+floor(2*(1/logSumm.frameDtVirmen))+1:logSumm.binned.camFrameID{trialidx(iTrial)}(end);
        end
        tempdff     = cat(3,tempdff,nanmean(dff(:,:,fidx),3)); % will be averaged across frames
        theseframes = fidx(1):fidx(1)+round(avgDff.fps*avgDff.nsec)-1;
        if theseframes(end) > nZ; continue; end
        tempfr      = tempfr + dff(:,:,theseframes);
      end
      
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).(avgDff.epochList{iEpoch}).frames = tempfr ./ numel(trialidx);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).(avgDff.epochList{iEpoch}).avg    = nanmean(tempdff,3);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).(avgDff.epochList{iEpoch}).sem    = nanstd(tempdff,0,3)./sqrt(numel(trialidx)-1);
    end
    
    %% add tower-triggered averages
    if iMaze > 1 && (strcmpi(avgDff.trialTypeList{iType},'correctL') || strcmpi(avgDff.trialTypeList{iType},'correctR'))
      for iTrial = 1:numel(trialidx)
        
        % basline-subtracted response for each R and L tower
        nr  = numel(logSumm.cuePos_R{trialidx(iTrial)});
        nl  = numel(logSumm.cuePos_L{trialidx(iTrial)});
        rml = nr - nl;
        for rr = 1:nr
          fidx          = logSumm.binned.camFrameID{trialidx(iTrial)} ...
            (find(logSumm.binned.time{trialidx(iTrial)} >= logSumm.cueOnset_R{trialidx(iTrial)}(rr),1,'first'));
          fidx          = fidx:fidx+floor(avgDff.nsecTower*avgDff.fps);
          thisdff       = dff(:,:,fidx);
          basel         = repmat(mean(dff(:,:,fidx(1)-2:fidx(1)-1),3),[1 1 size(thisdff,3)]);
          thisdff       = thisdff - basel;
          tempdff_tR    = cat(3,tempdff_tR,nanmean(dff(:,:,1:avgDff.towerRespNFrames),3));
          tempdff_tR_fr = tempdff_tR_fr + thisdff;
        end
        for ll = 1:nl
          fidx          = logSumm.binned.camFrameID{trialidx(iTrial)} ...
            (find(logSumm.binned.time{trialidx(iTrial)} >= logSumm.cueOnset_L{trialidx(iTrial)}(ll),1,'first'));
          fidx          = fidx:fidx+floor(avgDff.nsecTower*avgDff.fps);
          thisdff       = dff(:,:,fidx);
          basel         = repmat(mean(dff(:,:,fidx(1)-2:fidx(1)-1),3),[1 1 size(thisdff,3)]);
          thisdff       = thisdff - basel;
          tempdff_tL    = cat(3,tempdff_tL,nanmean(dff(:,:,1:avgDff.towerRespNFrames),3));
          tempdff_tL_fr = tempdff_tL_fr + thisdff;
        end
        
        % binned by #towers or delta, averaged between 180 & 230
        if nr >= max(avgDff.towers.ntowersBins); nr = max(avgDff.towers.ntowersBins)-1; end
        if nl >= max(avgDff.towers.ntowersBins); nl = max(avgDff.towers.ntowersBins)-1; end
        if abs(rml) >= max(abs(avgDff.towers.RminusLBins)); rml = sign(rml)*(max(abs(avgDff.towers.RminusLBins))-1); end
        
        fidx = logSumm.binned.camFrameID{trialidx(iTrial)}(logSumm.binned.pos{trialidx(iTrial)}(:,2) > 180 & logSumm.binned.pos{trialidx(iTrial)}(:,2) <= 230);
        
        binidx = find(avgDff.towers.ntowersBins <= nr,1,'last');
        tempdff_nt(binidx).R      = cat(3,tempdff_nt(binidx).R,nanmean(dff(:,:,fidx),3));
        
        binidx = find(avgDff.towers.ntowersBins <= nl,1,'last');
        tempdff_nt(binidx).L      = cat(3,tempdff_nt(binidx).L,nanmean(dff(:,:,fidx),3));
        
        binidx = find(avgDff.towers.RminusLBins <= rml,1,'last');
        tempdff_dt(binidx).RmL    = cat(3,tempdff_dt(binidx).RmL,nanmean(dff(:,:,fidx),3));
      end
    end
  end
  
end

%% averages and SEM for towers
avgDff.towersR.frames                  = tempdff_tR_fr ./ size(tempdff_tR,3);
avgDff.towersR.avg                     = nanmean(tempdff_tR,3);
avgDff.towersR.sem                     = nanstd(tempdff_tR,0,3)./sqrt(size(tempdff_tR,3)-1);
avgDff.towersL.frames                  = tempdff_tL_fr ./ size(tempdff_tL,3);
avgDff.towersL.avg                     = nanmean(tempdff_tL,3);
avgDff.towersL.sem                     = nanstd(tempdff_tL,0,3)./sqrt(size(tempdff_tL,3)-1);

for ii = 1:numel(avgDff.towers.ntowers)
  avgDff.towers.numtowers(ii).R.avg    = nanmean(tempdff_nt(ii).R,3);
  avgDff.towers.numtowers(ii).R.sem    = nanstd(tempdff_nt(ii).R,0,3)./sqrt(size(tempdff_nt(ii).R,3)-1);
  
  avgDff.towers.numtowers(ii).L.avg    = nanmean(tempdff_nt(ii).L,3);
  avgDff.towers.numtowers(ii).L.sem    = nanstd(tempdff_nt(ii).L,0,3)./sqrt(size(tempdff_nt(ii).L,3)-1);
end
for ii = 1:numel(avgDff.towers.RminusLvals)
  avgDff.towers.RminusL(ii).avg        = nanmean(tempdff_dt(ii).RmL,3);
  avgDff.towers.RminusL(ii).sem        = nanstd(tempdff_dt(ii).RmL,0,3)./sqrt(size(tempdff_dt(ii).RmL,3)-1);
end

%% save
save avgDff avgDff -v7.3
fprintf('\tdone after %1.1f min\n',toc/60)