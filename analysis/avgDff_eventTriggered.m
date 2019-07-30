function avgDff = avgDff_eventTriggered(rec,ROIflag,doZScore,baselSubt,baselZ)

% avgDff = avgDff_eventTriggered(rec,ROIflag,doZScore,baselSubt)
% turn and tower triggered Dff
%   rec: path for recording to be analyzed (default pwd)
%   ROIflag: true (default) to analyze ROI, false for pixels
%   doZscore: flag to zscore data (default true)
%   baselSubt: flag to subtract pre-event baseline (default false)
%   baselZ: flag to zscore to pre-event baseline (default true), overrides baselSubt
% LP jan 2019

%%
if nargin < 1 || isempty(rec);       rec       = pwd;   end
if nargin < 2 || isempty(ROIflag);   ROIflag   = true;  end
if nargin < 3 || isempty(doZscore);  doZScore  = false; end
if nargin < 4 || isempty(baselSubt); baselSubt = false; end
if nargin < 5 || isempty(baselz);    baselZ    = true;  end

%%
tic
fprintf('calculating average dffs')

%% load if available and return
rec = formatFilePath(rec);
cd(rec)
fn = 'avgDffEventTrig';
if ROIflag; fn = [fn '_ROI']; end

%% or else load data
load behavLog logSumm
if ROIflag
  load dffROI dffROI ROIlbl
  dff    = dffROI; 
  clear dffROI
else
  load dff dff
  ROIlbl     = {};
  [nX,nY,~]  = size(dff);
  [dff,bc]   = conditionDffMat(dff);
end
[nZ,nROI]    = size(dff);

%% some recs in blocks condensed have maze 12, this is just like maze 4 but with towers on both sides
% here they will be analyzed together for convenience
if contains(char(logSumm.info.protocol),'Condensed')
  logSumm.currMaze(logSumm.currMaze == 12) = 4;
end

%% get turns
[turnPoint,~,~,~,~,turnCfg] = estimateDecisionPoint(logSumm,[],[],false,'derivative');
avgDff.turnCfg              = turnCfg;

%% analysis parameters
% list type, presence of distractors, maze
avgDff.trialTypeList      = {'correct','error','correctR','correctL','errorR','errorL'};
avgDff.ROIlbl             = ROIlbl;
avgDff.mazeList           = unique(logSumm.currMaze);
avgDff.fps                = 1/logSumm.frameDtCam;
avgDff.respWinSec         = [.5 2];
avgDff.respNFrames        = round(avgDff.respWinSec .* avgDff.fps);
avgDff.timeAxis           = -avgDff.respNFrames(1)*(1/avgDff.fps):1/avgDff.fps:avgDff.respNFrames(2)*(1/avgDff.fps);
avgDff.timeAxis_preturb   = -avgDff.respNFrames(2)*(1/avgDff.fps):1/avgDff.fps:avgDff.respNFrames(1)*(1/avgDff.fps);
avgDff.isZscoredAll       = doZScore;
avgDff.isBaseSubt         = baselSubt;
avgDff.isZscoredBasel     = baselZ;

%% z-score
if doZScore; dff = zscore(dff); end

%% loop over trial types, mazes and epochs to compile average dff for
% combinations thereof
for iType = 1:numel(avgDff.trialTypeList)
  fprintf('.')
  for iMaze = 1:numel(avgDff.mazeList)
    
    %% initialize matrices
    if ROIflag
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).towerR.avg = nan(sum(avgDff.respNFrames)+1,nROI);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).towerR.sem = nan(sum(avgDff.respNFrames)+1,nROI);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).towerL.avg = nan(sum(avgDff.respNFrames)+1,nROI);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).towerL.sem = nan(sum(avgDff.respNFrames)+1,nROI);
      
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).turnR.avg  = nan(sum(avgDff.respNFrames)+1,nROI);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).turnR.sem  = nan(sum(avgDff.respNFrames)+1,nROI);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).turnL.avg  = nan(sum(avgDff.respNFrames)+1,nROI);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).turnL.sem  = nan(sum(avgDff.respNFrames)+1,nROI);
      
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).preturnR.avg  = nan(sum(avgDff.respNFrames)+1,nROI);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).preturnR.sem  = nan(sum(avgDff.respNFrames)+1,nROI);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).preturnL.avg  = nan(sum(avgDff.respNFrames)+1,nROI);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).preturnL.sem  = nan(sum(avgDff.respNFrames)+1,nROI);
    else
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).towerR.avg = nan(nX,nY,sum(avgDff.respNFrames)+1);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).towerR.sem = nan(nX,nY,sum(avgDff.respNFrames)+1);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).towerL.avg = nan(nX,nY,sum(avgDff.respNFrames)+1);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).towerL.sem = nan(nX,nY,sum(avgDff.respNFrames)+1);
      
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).turnR.avg  = nan(nX,nY,sum(avgDff.respNFrames)+1);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).turnR.sem  = nan(nX,nY,sum(avgDff.respNFrames)+1);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).turnL.avg  = nan(nX,nY,sum(avgDff.respNFrames)+1);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).turnL.sem  = nan(nX,nY,sum(avgDff.respNFrames)+1);
      
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).preturnR.avg  = nan(nX,nY,sum(avgDff.respNFrames)+1);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).preturnR.sem  = nan(nX,nY,sum(avgDff.respNFrames)+1);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).preturnL.avg  = nan(nX,nY,sum(avgDff.respNFrames)+1);
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).preturnL.sem  = nan(nX,nY,sum(avgDff.respNFrames)+1);
    end
    
    %% get list of relevant trials
    trialidx   = getTrialIdx(logSumm,avgDff.trialTypeList{iType},avgDff.mazeList(iMaze),widefieldParams.perfTh);
    if isempty(trialidx); continue; end
    
    %% triggered averages
    
    % compile trials
    turnDff_trials    = nan(sum(avgDff.respNFrames)+1,nROI,numel(trialidx));
    preturnDff_trials = nan(sum(avgDff.respNFrames)+1,nROI,numel(trialidx));
    tempdff_tR        = [];
    tempdff_tL        = [];
    
    for iTrial = 1:numel(trialidx)
      % turns
      if ~isnan(turnPoint(trialidx(iTrial)))
        fidx = logSumm.binned.camFrameID{trialidx(iTrial)} ...
                        (find(logSumm.binned.pos{trialidx(iTrial)}(:,2) <= turnPoint(trialidx(iTrial)),1,'last'));
        if fidx-avgDff.respNFrames(1) >= 1 || fidx+avgDff.respNFrames(2) <= nZ
          turnDff_trials(:,:,iTrial) = dff(fidx-avgDff.respNFrames(1):fidx+avgDff.respNFrames(2),:);
        end
      end
      
      % pre-turns (here the window is reversed, i.e. the baseline is
      % post-turn)
      if ~isnan(turnPoint(trialidx(iTrial)))
        fidx = logSumm.binned.camFrameID{trialidx(iTrial)} ...
                        (find(logSumm.binned.pos{trialidx(iTrial)}(:,2) <= turnPoint(trialidx(iTrial)),1,'last'));
        if fidx-avgDff.respNFrames(2) >= 1 || fidx+avgDff.respNFrames(1) <= nZ
          preturnDff_trials(:,:,iTrial) = dff(fidx-avgDff.respNFrames(2):fidx+avgDff.respNFrames(1),:);
        end
      end
      
      % towers
      nr  = numel(logSumm.cuePos_R{trialidx(iTrial)});
      nl  = numel(logSumm.cuePos_L{trialidx(iTrial)});
      for rr = 1:nr
        fidx          = logSumm.binned.camFrameID{trialidx(iTrial)} ...
                        (find(logSumm.binned.time{trialidx(iTrial)} <= logSumm.cueOnset_R{trialidx(iTrial)}(rr),1,'last'));
        if fidx-avgDff.respNFrames(1) < 1 || fidx+avgDff.respNFrames(2) > nZ; continue; end
        thisdff       = dff(fidx-avgDff.respNFrames(1):fidx+avgDff.respNFrames(2),:);
        tempdff_tR    = cat(3,tempdff_tR,thisdff);
      end
      for ll = 1:nl
        fidx          = logSumm.binned.camFrameID{trialidx(iTrial)} ...
                        (find(logSumm.binned.time{trialidx(iTrial)} <= logSumm.cueOnset_L{trialidx(iTrial)}(ll),1,'last'));
        if fidx-avgDff.respNFrames(1) < 1 || fidx+avgDff.respNFrames(2) > nZ; continue; end
        thisdff       = dff(fidx-avgDff.respNFrames(1):fidx+avgDff.respNFrames(2),:);
        tempdff_tL    = cat(3,tempdff_tL,thisdff);
      end
    end
    
    % subtract baseline
    if baselSubt && ~baselZ
      basel          = nanmean(turnDff_trials(1:avgDff.respNFrames(1),:,:),1);
      turnDff_trials = turnDff_trials - repmat(basel,[sum(avgDff.respNFrames)+1 1 1]);
      
      basel             = nanmean(preturnDff_trials(avgDff.respNFrames(2)+1:end,:,:),1);
      preturnDff_trials = preturnDff_trials - repmat(basel,[sum(avgDff.respNFrames)+1 1 1]);
      
      if ~isempty(tempdff_tR)
        basel        = nanmean(tempdff_tR(1:avgDff.respNFrames(1),:,:),1);
        tempdff_tR   = tempdff_tR - repmat(basel,[sum(avgDff.respNFrames)+1 1 1]);
      end
      if ~isempty(tempdff_tL)
        basel        = nanmean(tempdff_tL(1:avgDff.respNFrames(1),:,:),1);
        tempdff_tL   = tempdff_tL - repmat(basel,[sum(avgDff.respNFrames)+1 1 1]);
      end
    elseif baselZ
      basel          = nanmean(turnDff_trials(1:avgDff.respNFrames(1),:,:),1);
      baselSD        = nanstd(turnDff_trials(1:avgDff.respNFrames(1),:,:),0,1);
      turnDff_trials = turnDff_trials - repmat(basel,[sum(avgDff.respNFrames)+1 1 1]);
      turnDff_trials = turnDff_trials ./ repmat(baselSD,[sum(avgDff.respNFrames)+1 1 1]);
      
      basel          = nanmean(preturnDff_trials(avgDff.respNFrames(2)+1:end,:,:),1);
      baselSD        = nanstd(preturnDff_trials(avgDff.respNFrames(2)+1:end,:,:),0,1);
      preturnDff_trials = preturnDff_trials - repmat(basel,[sum(avgDff.respNFrames)+1 1 1]);
      preturnDff_trials = preturnDff_trials ./ repmat(baselSD,[sum(avgDff.respNFrames)+1 1 1]);
      
      if ~isempty(tempdff_tR)
        basel        = nanmean(tempdff_tR(1:avgDff.respNFrames(1),:,:),1);
        baselSD      = nanstd(tempdff_tR(1:avgDff.respNFrames(1),:,:),0,1);
        tempdff_tR   = tempdff_tR - repmat(basel,[sum(avgDff.respNFrames)+1 1 1]);
        tempdff_tR   = tempdff_tR ./ repmat(baselSD,[sum(avgDff.respNFrames)+1 1 1]);
      end
      if ~isempty(tempdff_tL)
        basel        = nanmean(tempdff_tL(1:avgDff.respNFrames(1),:,:),1);
        baselSD      = nanstd(tempdff_tL(1:avgDff.respNFrames(1),:,:),0,1);
        tempdff_tL   = tempdff_tL - repmat(basel,[sum(avgDff.respNFrames)+1 1 1]);
        tempdff_tL   = tempdff_tL ./ repmat(baselSD,[sum(avgDff.respNFrames)+1 1 1]);
      end
    end
    
    % avg and sem turns
    rightIdx = logSumm.choice(trialidx) == analysisParams.rightCode;
    turnR    = nanmean(turnDff_trials(:,:,rightIdx),3);
    turnRsem = nanstd(turnDff_trials(:,:,rightIdx),0,3)/sqrt(sum(rightIdx)-1);
    turnL    = nanmean(turnDff_trials(:,:,~rightIdx),3);
    turnLsem = nanstd(turnDff_trials(:,:,~rightIdx),0,3)/sqrt(sum(~rightIdx)-1);
    
    % avg and sem pre-turns
    preturnR    = nanmean(preturnDff_trials(:,:,rightIdx),3);
    preturnRsem = nanstd(preturnDff_trials(:,:,rightIdx),0,3)/sqrt(sum(rightIdx)-1);
    preturnL    = nanmean(preturnDff_trials(:,:,~rightIdx),3);
    preturnLsem = nanstd(preturnDff_trials(:,:,~rightIdx),0,3)/sqrt(sum(~rightIdx)-1);
    
    % avg and sem towers
    towR     = nanmean(tempdff_tR,3);
    towRsem  = nanstd(tempdff_tR,0,3)/sqrt(size(tempdff_tR,3)-1);
    towL     = nanmean(tempdff_tL,3);
    towLsem  = nanstd(tempdff_tL,0,3)/sqrt(size(tempdff_tL,3)-1);
    
    % put back in image format if necessary
    if ~ROIflag
      turnR    = conditionDffMat(turnR,bc,[],[nX nY sum(avgDff.respNFrames)+1]);
      turnRsem = conditionDffMat(turnRsem,bc,[],[nX nY sum(avgDff.respNFrames)+1]);
      turnL    = conditionDffMat(turnL,bc,[],[nX nY sum(avgDff.respNFrames)+1]);
      turnLsem = conditionDffMat(turnLsem,bc,[],[nX nY sum(avgDff.respNFrames)+1]);
      
      preturnR    = conditionDffMat(preturnR,bc,[],[nX nY sum(avgDff.respNFrames)+1]);
      preturnRsem = conditionDffMat(preturnRsem,bc,[],[nX nY sum(avgDff.respNFrames)+1]);
      preturnL    = conditionDffMat(preturnL,bc,[],[nX nY sum(avgDff.respNFrames)+1]);
      preturnLsem = conditionDffMat(preturnLsem,bc,[],[nX nY sum(avgDff.respNFrames)+1]);
      
      if ~isempty(towR)
        towR     = conditionDffMat(towR,bc,[],[nX nY sum(avgDff.respNFrames)+1]);
        towRsem  = conditionDffMat(towRsem,bc,[],[nX nY sum(avgDff.respNFrames)+1]);
      end
      if ~isempty(towL)
        towL     = conditionDffMat(towL,bc,[],[nX nY sum(avgDff.respNFrames)+1]);
        towLsem  = conditionDffMat(towLsem,bc,[],[nX nY sum(avgDff.respNFrames)+1]);
      end
    end
    
    % collect results
    avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).turnR.avg  = turnR;
    avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).turnR.sem  = turnRsem;
    avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).turnL.avg  = turnL;
    avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).turnL.sem  = turnLsem;
    
    avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).preturnR.avg  = preturnR;
    avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).preturnR.sem  = preturnRsem;
    avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).preturnL.avg  = preturnL;
    avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).preturnL.sem  = preturnLsem;
    
    if ~isempty(towR)
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).towerR.avg = towR;
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).towerR.sem = towRsem;
    end
    if ~isempty(towL)
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).towerL.avg = towL;
      avgDff.maze(iMaze).(avgDff.trialTypeList{iType}).towerL.sem = towLsem;
    end
    
  end
end

%% save
save(fn,'avgDff','-v7.3')
fprintf('\tdone after %1.1f min\n',toc/60)