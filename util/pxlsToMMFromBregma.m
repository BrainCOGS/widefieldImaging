function ROImm = pxlsToMMFromBregma(ROIcentroids,bregma,scale)

% ROImm = pxlsToMMFromBregma(ROIcentroids,bregma,scale)
% converts ROIcentroids from pixels to distances from bregma in mm

if nargin < 2 || isempty(bregma)
    load bregma bregma
end
if nargin < 3 || isempty(scale)
    scale = widefieldParams.pxlPerMM/widefieldParams.dsFactor;
end

nROI  = size(ROIcentroids,1);
ROImm = zeros(nROI,2);
for rr = 1:nROI
    ROImm(rr,:) = fliplr((bregma-ROIcentroids(rr,:))/scale);
    ROImm(rr,1) = - ROImm(rr,1);
end
      