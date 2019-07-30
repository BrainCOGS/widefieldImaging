function dffout = nanROIfromDff(dff,targetROIlbl,bilateralFlag)

% dffout = nanROIfromDff(dff,targetROIlbl,bilateralFlag)
% returns pxl x pxl x time matrix with nan masks for pixels belonging to ROI 
% (fromboth hemispheres if bilateral flag is set to true (default)
% if single ROI returns matrix, if multiple, returns cell array of matrices
% LP april 2018 lpinto@princeton.edu

%% input defaults
if nargin < 2; targetROIlbl  = {};   end
if nargin < 3; bilateralFlag = true; end

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

%% extract relevant pixels
if ~iscell(targetROIlbl)
  targetIdx    = strcmpi(ROIlbl,targetROIlbl);
  dffout       = dff;
  for iPxl = 1:size(ROI{targetIdx},1)
  	dffout(ROI{targetIdx}(iPxl,1),ROI{targetIdx}(iPxl,2),:) = nan;
  end
else
  dffout = cell(1,numel(targetROIlbl));
  for iROI = 1:numel(targetROIlbl)
    dffout{iROI} = dff;
    targetIdx = strcmpi(ROIlbl,targetROIlbl{iROI});
    for iPxl = 1:size(ROI{targetIdx},1)
      dffout{iROI}(ROI{targetIdx}(iPxl,1),ROI{targetIdx}(iPxl,2),:) = nan;
    end
  end
end