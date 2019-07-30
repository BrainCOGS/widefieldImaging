function axisHandle = drawROI(ROI,cl,axisHandle)

% axisHandle = drawROI(ROI,cl,axisHandle);
% ROI is cell array with [row col] coortdinates of each pixel (loads by default)
% cl is color of oulines (defaault 'w')
% axisHandle is axis handle

if nargin < 1 || isempty(ROI)
  if ~isempty(dir('ROIfromRef.mat')) 
    load ROIfromRef ROI
  end
end
if nargin < 2 || isempty(cl); cl = 'w'; end
if nargin < 3 || isempty(axisHandle); axisHandle = gca; end

if isempty(ROI)
  ROIoutline = [];
else
  ROIoutline = getPerimCoords(ROI,[128 128]);
end

axis(axisHandle);
hold on

for iROI = 1:numel(ROIoutline)
  plot(ROIoutline{iROI}(:,2),ROIoutline{iROI}(:,1),'-','linewidth',.25,'color',cl);
end

end

%% smooth area borders to be used with plot function
function bound = getPerimCoords(ROI,imsize)

bound = cell(numel(ROI),1);
for iROI = 1:numel(ROI)
  im       = zeros(imsize);
  for iPxl = 1:size(ROI{iROI},1); im(ROI{iROI}(iPxl,1),ROI{iROI}(iPxl,2)) = 1; end
  if iROI == 3 || iROI == 4 % v2 is disjoint, requires special treatment
    se = strel('line',5,1);
    im = imdilate(im,se);
  end
  centroid = round(mean(ROI{iROI}));
  row      = centroid(1);
  col      = find(im(row,:) > 0, 1, 'first');
  try
    bound{iROI} = bwtraceboundary(im,[row col],'N');
  catch
    col         = find(im(row,:) > 0, 1, 'last');
    bound{iROI} = bwtraceboundary(im,[row col],'S');
  end
  if numel(bound{iROI}(:,1)) < 20
    col         = find(im(row,:) > 0, 1, 'last');
    bound{iROI} = bwtraceboundary(im,[row col],'S');
  end
  if numel(bound{iROI}(:,1)) < 20 
    col         = find(im(row,:) > 0, 2, 'first');
    bound{iROI} = bwtraceboundary(im,[row col(2)],'N');
  end
end

end