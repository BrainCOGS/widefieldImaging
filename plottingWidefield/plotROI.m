function fh = plotROI(ROI,imsize2,fh)

% fh = plotROI(ROI,imsize2,fh)
% plots color-coded ROIs whose pixels are provided in ROI
% imsize2 is a 2-element vector with image size
% fh is optional figure handle
%
% LP sep 2016

if nargin < 2
  load dff meanproj
  [nX,nY] = size(meanproj);
  imsize2 = [nX nY];
end
if nargin < 3
  figure;
  fh = gca;
end

im     = zeros(imsize2(1),imsize2(2),3);
colors = feval(widefieldParams.colors,numel(ROI)+1);
colors(sum(colors,2)==0,:) = [];
for iROI = 1:numel(ROI)
  for iPxl = 1:size(ROI{iROI},1)
    im(ROI{iROI}(iPxl,1),ROI{iROI}(iPxl,2),:) = colors(iROI,:);
  end
end

axis(fh);
image(im);
axis image; axis off;
for iROI = 1:numel(ROI)
  text(median(ROI{iROI}(:,2)),median(ROI{iROI}(:,1)),...
    num2str(iROI),'fontsize',14,'color',[.7 .7 .7],'fontweight','bold');
end
