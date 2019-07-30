function dffSingleROI = getSingleROIpxlDff(dff,targetROIlbl,bilateralFlag,returnNan)

% dffSingleROI = getSingleROIpxlDff(dff,targetROIlbl,bilateralFlag)
% returns time x pxl matrix for subset of pixels belonging to ROI (from
% both hemispheres if bilateral flag is set to true (default)
% if single ROI returns matrix, if multiple, returns cell array of matrices
% dff can be either a 3D matrix or recording path
% LP dec 2017 lpinto@princeton.edu

%% input defaults
if nargin < 1 || isempty(dff); dff = pwd;   end
if nargin < 2; targetROIlbl        = {};    end
if nargin < 3; bilateralFlag       = true;  end
if nargin < 4; returnNan           = false; end

%% handle dff
if ischar(dff)
  rec = dff; clear dff 
  cd(rec)
  load dff dff
end

%% get relevant (bilateral) ROIs
load ROIfromRef ROIlbl ROI
if bilateralFlag
  ROIlbl  = unique(cellfun(@(x)(x(1:end-2)),ROIlbl,'UniformOutput',false),'stable');
  unilROI = ROI;
  clear ROI; ct = 1;
  for iROI = 1:2:numel(unilROI)
    ROI{ct} = [unilROI{iROI}; unilROI{iROI+1}];
    ct      = ct+1;
  end
end

if isempty(targetROIlbl); targetROIlbl = ROIlbl; end

%% downsample ROIs if necessary
if size(dff,1) < 128
  dsFactor = size(dff,1) ./ 128;
  
  for iROI = 1:numel(ROI)
    thisim = zeros(128,128);
    for iPxl = 1:size(ROI{iROI},1)
      thisim(ROI{iROI}(iPxl,1),ROI{iROI}(iPxl,2)) = 1;
    end
    thisim    = imresize(thisim,dsFactor);
    [r,c]     = find(thisim > .5);
    ROI{iROI} = [r c];
  end
  
end

%% extract relevant pixels
if ~iscell(targetROIlbl)
  targetIdx    = strcmpi(ROIlbl,targetROIlbl);
  temp         = arrayfun(@(x,y)(squeeze(dff(x,y,:))),ROI{targetIdx}(:,1),ROI{targetIdx}(:,2),'UniformOutput',false);
  dffSingleROI = cell2mat(temp');
  nancol       = sum(isnan(dffSingleROI)) > 0;
  dffSingleROI(:,nancol) =[];
else
  dffSingleROI = cell(1,numel(targetROIlbl));
  for iROI = 1:numel(targetROIlbl)
    targetIdx = strcmpi(ROIlbl,targetROIlbl{iROI});
    temp      = arrayfun(@(x,y)(squeeze(dff(x,y,:))),ROI{targetIdx}(:,1),ROI{targetIdx}(:,2),'UniformOutput',false);
    temp      = cell2mat(temp');
    if ~returnNan
      nancol    = sum(isnan(temp)) > 0;
      temp(:,nancol)     = [];
    end
    dffSingleROI{iROI} = temp;
  end
end