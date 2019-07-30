function [dff,B,V,dffBlue,dffViolet,cf,cfg,Vorig] = dffFrom2Channels(rawf,frameRate,cfg)

% [dff,F0,B,V,dffBlue,dffViolet,cf,cfg] = dffFrom2Channels(rawf,frameRate,cfg)
% implements isosbestic hemodynamic correction
% 
% INPUT
%  rawf is 3D image stack with raw florescence values
%  frameRate is imaging rate (default 20)
%  cfg is analysis parameter structure with the following fields
%     cfg.winSizeSec     = 30; % window for running F0
%     cfg.dffMethod      = 'rawfFitDivide'; % method for correction
%     cfg.runWinMethod   = 'overall'; % linear fit of violet (green) channel in chunks or running baseline
%     cfg.smoothBaseline = true; % smooth F0 before subtracting? (relevant for dffMethod 'channelF0')
%     cfg.smoothDff      = false; % smooth dff before subtracting?
%     cfg.smoothSigma    = 3; % gaussian filter sigma for smoothing
%     cfg.strobeSeq      = [1 2]; sequence of blue (1) and other LED (2)
% OUPUT
%       dff: corrected image stack
%         B: substack with rawf for blue illumination
%         V: substack with rawf for violet illumination
%   dffBlue: deltaF/F for B
% dffViolet: deltaF/F for V
%        cf: heuristic correction factor
%       cfg: copy of analysis parameter structure
%     Vorig: V before smoothing if applicable
%
% Lucas Pinto (lpinto@princeton.edu)

%% defaults
if nargin < 2
  try
    load info frameRate
  catch
    frameRate = 20;
  end
end

if nargin < 3
  cfg = struct([]);
end

cfg   = populateCfg(cfg);

tic
fprintf('calculating hemodynamic-corrected dff, method: %s...\n ',cfg.dffMethod)

%% image size, reshape matrix for convenience
[nX,nY,nZ]   = size(rawf);
[temp,bc,br] = conditionDffMat(double(rawf));

%% select frames of different LEDs and interpolate
strobeSeq    = repmat(cfg.strobeSeq,[1 floor(nZ/numel(cfg.strobeSeq))]);
excFrames    = nZ - numel(strobeSeq);
strobeSeq    = [strobeSeq cfg.strobeSeq(1:excFrames)];
B            = temp(strobeSeq==1,:);
Vpre         = temp(strobeSeq==2,:);
V            = interp1(find(strobeSeq==2)',Vpre,find(strobeSeq==1)','linear','extrap');
V(1,:)       = Vpre(1,:);
nZtrue       = size(B,1);

%% calculate dff
switch cfg.dffMethod
  % green/violet channel is taken to be F0
  case 'channelF0'
    % linear fit for magnitude matching
    F0   = runningBaselineFit(B,V,cfg.winSizeSec,frameRate,cfg.runWinMethod,cfg.smoothBaseline);
    dff  = (B-F0)./F0;
    dffb = nan(size(dff));
    dffv = nan(size(dff));
    Vorig= V;
    % calculate dff separetely and subtract (divide) other channel from blue
  case {'channelSubtract','channelDivide'}
    fprintf('\t')
    F0   = computeF0mode(B,cfg.winSizeSec,frameRate,false); % F0 is a runnig mode
    dffb = (B-F0)./F0;
    fprintf('\t')
    F0v  = computeF0mode(V,cfg.winSizeSec,frameRate,false); % F0 is a runnig mode
    dffv = (V-F0v)./F0v;
    
    if cfg.smoothDff % smooth violet (green) channel before subtracting?
      dffv  = approxGaussianBlur1D(dffv,cfg.smoothSigma,1);
    end
    
    % linear fit for magnitude matching
    fprintf('\t')
    [dffv,cf] = runningBaselineFit(dffb,dffv,cfg.winSizeSec*2,frameRate,cfg.runWinMethod,cfg.smoothBaseline);
    
    % subtract or divide
    if strcmpi(cfg.dffMethod,'channelSubtract')
      dff  = dffb - dffv;
    elseif strcmpi(cfg.dffMethod,'channelDivide')
      dff  = dffb ./ dffv;
    end
    Vorig  = V;
    
    % new method: fit rawfV to rawfB, calculate ratio, divide, subtract 1
  case 'rawfFitDivide'
    
    % linear fit for magnitude matching
    fprintf('\t')
    Vorig  = V;
    [V,cf] = runningBaselineFit(B,V,cfg.winSizeSec,frameRate,'overall',cfg.smoothBaseline);
    fprintf('\t')
    F0     = computeF0mode(B,cfg.winSizeSec,frameRate,false); % F0 is a runnig mode
    rb     = B./F0;
    fprintf('\t')
    F0v    = computeF0mode(V,cfg.winSizeSec,frameRate,false); % F0 is a runnig mode
    rv     = V./F0v;
    
    dff    = (rb ./ rv) - 1;
    dffb   = rb - 1;
    dffv   = rv - 1;
    
  case 'rawfOffsetDivide'
    
    % linear fit for magnitude matching
    fprintf('\t')
    Vorig  = V;
    [V,cf] = runningBaselineFit(B,V,cfg.winSizeSec,frameRate,'offset',cfg.smoothBaseline);
    fprintf('\t')
    F0     = computeF0mode(B,cfg.winSizeSec,frameRate,false); % F0 is a runnig mode
    rb     = B./F0;
    fprintf('\t')
    F0v    = computeF0mode(V,cfg.winSizeSec,frameRate,false); % F0 is a runnig mode
    rv     = V./F0v;
    
    dff    = (rb ./ rv) - 1;
    dffb   = rb - 1;
    dffv   = rv - 1;
end

%% matrices back into original shape / size
% F0           = conditionDffMat(F0,bc,br,[nX nY nZtrue]);
dff          = conditionDffMat(dff,bc,br,[nX nY nZtrue]);
dffBlue      = conditionDffMat(dffb,bc,br,[nX nY nZtrue]);
dffViolet    = conditionDffMat(dffv,bc,br,[nX nY nZtrue]);
B            = conditionDffMat(B,bc,br,[nX nY nZtrue]);
V            = conditionDffMat(V,bc,br,[nX nY nZtrue]);
Vorig        = conditionDffMat(Vorig,bc,br,[nX nY nZtrue]);

fprintf('\tdone after %1.1f min\n',toc/60)

end

%% cfg defaults
function cfg = populateCfg(cfg)

if ~isfield(cfg,'winSizeSec')
  cfg(1).winSizeSec            = 30;
end
if ~isfield(cfg,'dffMethod')
  cfg(1).dffMethod             = 'rawfFitDivide';
end
if ~isfield(cfg,'runWinMethod')
  cfg(1).runWinMethod          = 'overall';
end
if ~isfield(cfg,'smoothBaseline')
  cfg(1).smoothBaseline        = true;
end
if ~isfield(cfg,'smoothDff')
  cfg(1).smoothDff             = false;
end
if ~isfield(cfg,'smoothSigma')
  cfg(1).smoothSigma           = 3;
end
if ~isfield(cfg,'strobeSeq')
  cfg(1).strobeSeq             = [1 2];
end

end