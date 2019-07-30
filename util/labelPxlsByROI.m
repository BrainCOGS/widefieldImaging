function imlbl = labelPxlsByROI(ROI,binSize)

% imlbl = labelPxlsByROI(ROI,binSize)
% return image with number labels corresponding to ROIs

if nargin < 1; ROI     = []; end
if nargin < 2; binSize = 1;  end

%%
if isempty(ROI); load ROIfromRef ROI; end

%%
im = nan(128,128);
for iROI = 1:numel(ROI)
  for iPxl = 1:size(ROI{iROI},1)
    im(ROI{iROI}(iPxl,1),ROI{iROI}(iPxl,2)) = iROI;
  end
end

if binSize == 1
  imlbl = im;
else
  temp  = imresize(im,1/binSize,'nearest');
  imlbl = zeros(size(temp));

  for iROI = 1:numel(ROI)
    if iROI == 1
      [r,c] = find(temp >= .5 & temp < iROI+1);
    elseif iROI == numel(ROI)
      [r,c] = find(temp >= iROI);
    else
      [r,c] = find(temp >= iROI & temp < iROI+1);
    end
    for iPxl = 1:numel(r)
      imlbl(r(iPxl),c(iPxl)) = iROI;
    end
  end
end
