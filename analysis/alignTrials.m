function [dffTrials,LorR,taxis,trialidx] = alignTrials(dff,logSumm,whichTrialType,maze,alignPoint,trange,timeOrSpace,combineVisGuide,avgPos)

% [dffTrials,LorR,taxis] = alignTrials(dff,logSumm,whichTrialType,maze,alignPoint,trange,timeOrSpace,combineVisGuide,avgPos)
% INPUT
%   dff is a frames x pxl(ROI) matrix
%   logSumm is flattened virmen log
%   whichTrialType: string for trial type, eg 'correct' (default), 'all'; see getTrialIdx() for full list
%   maze: maze number
%   alignPoint: affects time aignment only; 'cueStart' (default), 'cueHalf',
%               'memStart', 'memHalf', 'startStart', 'startHalf', 'turnStart',
%               'rewardStart', 'itiStart'
%   trange: how much to include in trials, in sec, eg [0 4] (default)
%   timeOrSpace: align in 'space' (default) or 'time'
%   combineVisGuide: flag to combine visually-guided mazes with and
%                    without distractors (default true)
%   avgPos: true to average over entire position bin, false (default) to
%           just sample first time frame >= that position
% OUTPUT dffTrials is a trial x time x pxl (ROI) matrix

%% defaults
if nargin < 2 || isempty(logSumm)
  load behavLog logSumm
end
if nargin < 3 || isempty(whichTrialType)
  whichTrialType = 'correct';
end
if nargin < 4 || isempty(maze)
  maze       = [];
end
if nargin < 5 || isempty(alignPoint)
  alignPoint = 'cueStart';
end
if nargin < 6 || isempty(trange)
  trange = [0 4];
end
if nargin < 7 || isempty(timeOrSpace)
  timeOrSpace = 'space';
end
if nargin < 8
  combineVisGuide = true;
end
if nargin < 9
  avgPos      = false;
end

if combineVisGuide
  logSumm.currMaze(logSumm.currMaze == 12) = 4;
end
if isempty(maze)
  maze = max(logSumm.currMaze);
end
%% get trial IDs, initialize matrices
trialidx   = getTrialIdx(logSumm,whichTrialType,maze,widefieldParams.perfTh);
LorR       = logSumm.trialType(trialidx);

if numel(trange) > 2 % if trange has multiople entries assume it's actually taxis
  taxis  = trange;
  frames = round(trange./logSumm.frameDtCam(1));
else
  switch timeOrSpace
    case 'time'
      taxis  = trange(1):logSumm.frameDtCam:trange(2);
      frames = round(trange./logSumm.frameDtCam(1));
    case 'space'
      taxis  = trange(1):trange(2);
  end
end

if sum(size(dff) == 1) < 1
  nROI = size(dff,2);
else
  nROI = 1;
  if size(dff,1) < size(dff,2); dff = dff'; end
end
dffTrials  = zeros(numel(trialidx),numel(taxis),nROI);
if avgPos && strcmpi(timeOrSpace,'space')
  dffTrials(:,end,:) = [];
end

%% align Trials
for nR = 1:nROI
  for tt = 1:numel(trialidx)
    switch timeOrSpace
      case 'time'
        switch alignPoint
          case 'startStart'
            f0idx = logSumm.binned.camFrameID{trialidx(tt)}(1);
          case 'startHalf'
            f0idx = logSumm.binned.camFrameID{trialidx(tt)}(find(logSumm.binned.pos{trialidx(tt)}(:,2)>=-15,1,'first'));
          case 'cueStart'
            f0idx = logSumm.binned.keyFrames{trialidx(tt)}(strcmpi(logSumm.keyFrameLabels,'cue'));
          case 'cueHalf'
            f0idx = logSumm.binned.camFrameID{trialidx(tt)}(find(logSumm.binned.pos{trialidx(tt)}(:,2)>=100,1,'first'));
          case 'memStart'
            f0idx = logSumm.binned.keyFrames{trialidx(tt)}(strcmpi(logSumm.keyFrameLabels,'mem'));
          case 'memHalf'
            f0idx = logSumm.binned.camFrameID{trialidx(tt)}(find(logSumm.binned.pos{trialidx(tt)}(:,2)>=250,1,'first'));
          case 'turnStart'
            f0idx = logSumm.binned.camFrameID{trialidx(tt)}(find(logSumm.binned.pos{trialidx(tt)}(:,2)>=280,1,'first'));
          case 'rewardStart'
            f0idx = logSumm.binned.keyFrames{trialidx(tt)}(strcmpi(logSumm.keyFrameLabels,'rew'));
          case 'itiStart'
            f0idx = logSumm.binned.camFrameID{trialidx(tt)}(size(logSumm.binned.pos{trialidx(tt)},1)+1);
        end
        idx       = f0idx+frames(1):min([f0idx+frames(2) numel(dff)]);
        dffTrials(tt,1:numel(idx),nR) = dff(idx,nR);
        
      case 'space'
        
        if avgPos
          for iPos = 1:numel(taxis)-1
            idx1 = find(logSumm.binned.pos{trialidx(tt)}(:,2) >= taxis(iPos),1,'first');
            idx2 = find(logSumm.binned.pos{trialidx(tt)}(:,2) <  taxis(iPos+1),1,'last');
            idx  = logSumm.binned.camFrameID{trialidx(tt)}(idx1):logSumm.binned.camFrameID{trialidx(tt)}(idx2);
            dffTrials(tt,iPos,nR) = nanmean(dff(idx,nR));
          end
          
        else
          idx = arrayfun(@(x)(find(logSumm.binned.pos{trialidx(tt)}(:,2) >= x,1,'first')),taxis);
          idx = logSumm.binned.camFrameID{trialidx(tt)}(idx);
          dffTrials(tt,1:numel(idx),nR) = dff(idx,nR);
        end
        
    end
  end
end

if avgPos && strcmpi(timeOrSpace,'space')
  taxis = taxis(1:end-1)+mode(diff(taxis))/2;
end