function F0 = computeF0mode(rawf,ws,frameRate,reshapeMat)

% F0 = computeF0mode(rawf,ws,frameRate, reshapeMat)
% computes mode of dff values from raw flourescence data F
% INPUT 
%   rawf: F matrix (pixels x pixels x frames)
%   ws: running window size in sec
%   frameRate: recording frame rate
%   reshape: true to return output in pixels x pixels x frames (default)
%   uses Sue Ann's mex file halfSampleMode

%% try loading info by default
if nargin < 2 || isempty(ws)
  ws  = widefieldParams.winSizeModeSec; % window size in sec
end
if nargin < 3
  
  try
    load info frameRate
  catch
    frameRate = widefieldParams.frameRate;
  end
  
end
if nargin < 4
  reshapeMat = 1;
end

tic

fprintf('calculating F0...') 

wsSamp     = ws*frameRate;

% estimate F0 as running mode
if reshapeMat
    temp = reshape(rawf,[],size(rawf,3));
    F0   = halfSampleMode(temp',wsSamp); % mex file
    F0   = reshape(F0',size(rawf));
else
    F0   = halfSampleMode(rawf,wsSamp);
end

t = toc;
fprintf(' done after %1.2f minutes\n',t/60)