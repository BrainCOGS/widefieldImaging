function [ROIperim,bregma,midline] = getSensoryMapOutlines(rec)

% ROIperim = getSensoryMapOutlines(rec)

if nargin < 1; rec = pwd; end
rec = formatFilePath(rec);

%% paths
[mouse,recdate]  = mouseAndDateFromFileName(rec);

if isempty(recdate) 
  % if no rec date, assume it's concat sessions
  load([rec 'sensoryMaps.mat'])
  ROIperim     = sensoryMaps.ROIoutline;
  bregma       = sensoryMaps.bregma;
  midline      = sensoryMaps.midline;
  
else
%   stopAt       = strfind(rec,recdate)-1;
%   rootdir      = rec(1:stopAt);
  wf           = widefieldParams;
  rootdir      = formatFilePath([wf.getRootDir(isThisSpock) mouse]);
  refFn        = [rootdir 'refIm.mat'];
  mapFn        = [rootdir 'sensoryMaps.mat'];

  %% register image to reference if necessary
  if isempty(dir('reg2ref.mat')) && ~isempty(dir(refFn)) && ~isempty(dir(mapFn))
    tic; fprintf('registering rec... ')
    load dff meanproj
    im                       = imBinSpace(fliplr(rot90(meanproj,3)));
    load(refFn,'meanproj','bregma','midline')
    [tform,~,bregma,midline] = registerRecs(im,meanproj,bregma,midline,'rigid',false);

    load(mapFn,'sensoryMaps')
    [~,~,~,~,ROIperim] = getRefROI(sensoryMaps.areaCoord,[],tform,im,false);
    fprintf('done after %1.1fmin\n',toc/60)

  elseif ~isempty(dir('reg2ref.mat')) && ~isempty(dir(refFn)) && ~isempty(dir(mapFn))
    load reg2ref tform im bregma midline
    load(mapFn,'sensoryMaps')
    [~,~,~,~,ROIperim] = getRefROI(sensoryMaps.areaCoord,[],tform,im,false);

  else
    ROIperim             = [];
    bregma               = [];
    midline              = [];
%     load dff meanproj
%     im                   = imBinSpace(fliplr(rot90(meanproj,3)));
%     load(refFn,'meanproj','bregma','midline')
%     [~,~,bregma,midline] = registerRecs(im,meanproj,bregma,midline,false);

  end
end