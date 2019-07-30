function [ROI,ROIlbl,ROIplotOrder,ROIcentroids,ROIperim] = getRefROI(ROIpath,tform,im,saveFlag)

% [ROI,ROIlbl,ROIcentroids,ROIperim] = getRefROI(ROIref,tform,im,saveFlag)
% transform reference ROIs (normally from Allen Brain Atlas) into
% present coordinates given by imaging registration (see registerRecs)
% will load everything from disk if no input is provided
% INPUT
%   ROIpath is folder containing refROI.mat
%   tform is registerRecs() output (im transformation object)
%   im is refernce image (typically meanProj)
%   saveFlag: true to save ROIfromRef.mat (default true)
% 
% OUTPUT
%   ROI is cell array with pixel coordinates belonging to each ROI
%   ROIlbl is cell array with corresponding ROI names
%   ROIplotOrder is loaded from refROI  and passed along, list of indices
%   determining plotting order for downstream functions
%   ROIcentroids in pxls: ROIs x [row col] matrix
%   ROIperim: cell array with outlines of the ROIs
%
% LP sep 2016

%%
if nargin < 1 || isempty(ROIpath)
  mouseID = mouseAndDateFromFileName(pwd);
  wf      = widefieldParams;
  root    = wf.getRootDir(isThisSpock,~contains(pwd,'/Volumes'));
  ROIpath = [formatFilePath(root) mouseID '/refROI.mat'];
end
if nargin < 2 || isempty(tform)
  load reg2ref tform
end
if nargin < 3 || isempty(im)
  load reg2ref im
end
if nargin < 4
  saveFlag = true;
end

%% 
load(ROIpath,'refROI')
ROIref   = refROI.areaCoord;
ROIlbl   = refROI.areaLbl;

%% apply transformation
[nX,nY]      = size(im);
ROI          = cell(size(ROIref));
ROIcentroids = zeros(length(ROIref),2);
ROIperim     = zeros(nX,nY);
for ii = 1:length(ROI)
  [ROI{ii}(:,2),ROI{ii}(:,1)] = ...
    transformPointsInverse(tform,ROIref{ii}(:,2),ROIref{ii}(:,1));
  ROI{ii}            = round(ROI{ii});
  delidx = ROI{ii}(:,1) < 1 | ROI{ii}(:,1) > nY | ...
           ROI{ii}(:,2) < 1 | ROI{ii}(:,2) > nX;
  ROI{ii}(delidx,:)  = [];
  ROIcentroids(ii,:) = nanmean(ROI{ii});
  
  temp = zeros(nX,nY);
  for jj = 1:size(ROI{ii},1)
    temp(ROI{ii}(jj,1),ROI{ii}(jj,2)) = 1;
  end
  ROIperim = bwperim(temp) + ROIperim;
end

ROIperim(ROIperim > 1) = 1;

if ~exist('ROIplotOrder','var')
  ROIplotOrder = 1:numel(ROI);
end

if saveFlag
  save ROIfromRef ROI ROIlbl ROIperim ROIcentroids ROIplotOrder
end


