function [ROIlbl,ROImm] = labelROIs(ROIcentroids,bregma,lblmethod,ROIfn,saveFlag)

% [ROIlbl,ROImm] = labelROIs(ROIcentroids,bregma,lblmethod,ROIfn)
% automatically or manually labels ROIs
% ROIcentroids is nROI x 2 vector specifying centroids in cartesian space
% bregma is cartesian coordinate for bregma pixel
% lblmethod: 'auto' (from atlas coordinates, usually fails for very large ROIs)
%            or 'manual'
% ROIfn:     file name containing ROIs, to both load vars and save ROIlbl
%            in case saveFlag is set to 1 
% output 
% ROIlbl: cell with string labels for ROIs
% ROImm:  conversion into [ML AP] distance from bregma
%
% LP aug 2016

if nargin < 2 || isempty(bregma)
    load bregma bregma
end
if nargin < 3 || isempty(lblmethod)
    lblmethod = 'auto';
end
if nargin < 4 || isempty(ROIfn)
    ROIfn = 'manualROI';
end
if nargin < 5
    saveFlag = 0;
end

nROI   = size(ROIcentroids,1);
ROIlbl = cell(nROI,1);
ROImm  = pxlsToMMFromBregma(ROIcentroids,bregma);

switch lblmethod
    case 'auto'
        for rr = 1:nROI
            ROIlbl{rr} = getLocationLabel(ROImm(rr,:));
        end
        
    case 'manual'
        load(ROIfn,'meanproj','ROI')
        for rr = 1:nROI
            figure; hold on
            imagesc(meanproj); colormap gray
            axis image
            axis ij
            axis off
            for ii = 1:size(ROI{rr})
                plot(ROI{rr}(ii,2),ROI{rr}(ii,1),'rs','markersize',1.2,'markerfacecolor','r');
            end
            ROIlbl{rr} = input('ROI label:  ');
            close
        end
end

if saveFlag
    save(ROIfn,'ROIlbl','ROImm','lblmethod','-append')
end