function [avg,sem] = getPreTrialBaseline(dff,logSumm,maze,ROIflag,cfg)

% [avg,sem] = getPreTrialBaseline(dff,logSumm,maze,ROIflag,cfg)
% for each ROI, get avg and sem over pre-trial baseline (def. 0.5 sec)

if nargin < 1 || isempty(dff);     dff     = pwd;  end
if nargin < 2 || isempty(logSumm); logSumm = [];   end
if nargin < 3 || isempty(maze);    maze    = [];   end
if nargin < 4 || isempty(ROIflag); ROIflag = true; end
if nargin < 5 || isempty(cfg);     cfg     = [];   end

%% load stuff if necessary, condition matrix
if isempty(cfg)
  cfg.trialType = 'correct'; % pre-baseline before these trials
  cfg.time      = 0.5; % sec before trial start
end

if ischar(dff)
  dff = formatFilePath(dff);
  if ROIflag
    load([dff 'dffROI.mat'],'dffROI');
    dff = dffROI;
  else
    load([dff 'dff.mat'],'dff');
  end
end

if size(dff,3) > 1
  [nX,nY,~] = size(dff);
  [dff,bc]  = conditionDffMat(dff);
  doReshape = true;
else
  doReshape = false;
end

if isempty(logSumm)
  load behavLog logSumm
end

% some recs in blocks condensed have maze 12, this is just like maze 4 but with towers on both sides
% here they will be analyzed together for convenience
if ~isempty(strfind(char(logSumm.info.protocol),'Condensed'))
  logSumm.currMaze(logSumm.currMaze == 12) = 4;
end

if isempty(maze)
  maze = max(logSumm.currMaze);
end

load info frameRate
frameRate = round(frameRate/2);
nframes   = round(frameRate*cfg.time);

%% get trial indices
trialidx  = getTrialIdx(logSumm,cfg.trialType,maze);

%% for each ROI (pxl) claculate baseline
nROI      = size(dff,2);
avg       = zeros(nROI,1);
sem       = zeros(nROI,1);

for iROI = 1:nROI
  firstFrames = cellfun(@(x)(x(1)),logSumm.camFrameNum(trialidx));
  frames      = zeros(numel(trialidx),nframes);
  for iTrial = 1:numel(trialidx)
    if firstFrames(iTrial)-nframes < 1
      excessF          = abs(firstFrames(iTrial)-nframes)+1;
      frames(iTrial,:) = [nan(1,excessF) dff(1:firstFrames(iTrial)-1,iROI)'];
    else
      frames(iTrial,:) = dff(firstFrames(iTrial)-nframes:firstFrames(iTrial)-1,iROI)';
    end
  end
  mframe    = frames(:);
%   mframe    = nanmean(frames); % average over trials first, sem will be over frames
  avg(iROI) = mean(mframe);
  sem(iROI) = std(mframe)./sqrt(numel(mframe)-1);
end

%% if in pxl return matrices in image format
if doReshape
  avg = conditionDffMat(avg',bc,[],[nX nY 1]);
  sem = conditionDffMat(sem',bc,[],[nX nY 1]);
end

