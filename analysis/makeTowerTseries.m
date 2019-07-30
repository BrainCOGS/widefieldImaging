function [tL,tR] = makeTowerTseries(logSumm,tidx,alignPoint,lengthSec,frameRate)

% [tL,tR] = makeTowerTseries(logSumm,tidx,alignPoint,lengthSec,frameRate)
% generate matrices where cols are delta functions for each tower in each
% trial, for each side

% handle input
if nargin < 2 || isempty(tidx)
    tidx = getTrialIdx(logSumm,'correct');
end
if nargin < 3 || isempty(alignPoint)
    alignPoint = 'cueStart'; % only one supported for now
end
if nargin < 4 || isempty(lengthSec)
    lengthSec = 15;
end
if nargin < 5 || isempty(frameRate)
    try
        load info frameRate
    catch
        frameRate = widefieldParams.frameRate;
    end
end

% initialize matrices
ntrials = numel(tidx);
tL      = zeros(round(lengthSec*frameRate),ntrials);
tR      = tL;

% find frames when towers occur
for ii = 1:ntrials
    
    nR  = numel(logSumm.cuePos_R{tidx(ii)});
    nL  = numel(logSumm.cuePos_L{tidx(ii)});
    ttl = zeros(nL,1); 
    ttr = zeros(nR,1);
    
    % find frame where each tower occurs
    for jj = 1:nR
        idx     = find(logSumm.binned.pos{tidx(ii)}(:,2) >= logSumm.cuePos_R{tidx(ii)}(jj),1,'first');
        ttr(jj) = logSumm.binned.camFrameID{tidx(ii)}(idx);
    end
    for jj = 1:nL
        idx     = find(logSumm.binned.pos{tidx(ii)}(:,2) >= logSumm.cuePos_L{tidx(ii)}(jj),1,'first');
        ttl(jj) = logSumm.binned.camFrameID{tidx(ii)}(idx);
    end
    
    % subtract align point to fit in desired time axis
    switch alignPoint
        case 'startStart'
            offset = logSumm.binned.camFrameID{tidx(ii)}(1);
        case 'cueStart'
            offset = logSumm.binned.keyFrames{tidx(ii)}(strcmpi(logSumm.keyFrameLabels,'cue'));
    end
    
    % times where towers occur = 1
    ttr        = ttr - offset + 1;
    ttl        = ttl - offset + 1;
    tL(ttl,ii) = 1;
    tR(ttr,ii) = 1;
end