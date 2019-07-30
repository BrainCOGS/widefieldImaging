function [cl,lbl] = getDefaultROIcl(ROIlbl)

% [cl,lbl] = getDefaultROIcl(ROIlbl)
% return default RGB color and labels according to my convention from cell
% array with strings of ABI names

% compare to defaults
ROIlbls_allen = {'VISp','mV2','VISa','RSP','mMOs','aMOs','MOp','SS'};
ROIlbls_mine  = {'V1','mV2','PPC','RSC','mM2','aM2','M1','SS'};
colors        = widefieldParams.areaCl;

% find color and new label
doRevert = false;
if ~iscell(ROIlbl); ROIlbl = {ROIlbl}; doRevert = true; end
for iROI = 1:numel(ROIlbl)
  % remove side lbl temporarily
  thislbl = ROIlbl{iROI};
  if ~isempty(strfind(thislbl,'-'))
    addback = thislbl(end-1:end);
    thislbl = thislbl(1:end-2);
  else
    addback = '';
  end
  idx         = strcmpi(ROIlbls_allen,thislbl);
  cl(iROI,:)  = colors(idx,:);
  lbl{iROI}   = [ROIlbls_mine{idx} addback];
end

if doRevert; lbl = lbl{1}; end