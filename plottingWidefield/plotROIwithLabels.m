function fh = plotROIwithLabels(ROI,ROIlbl,imsize2,fh)

if nargin < 3
    imsize2 = [128 128];
end
if nargin < 4
    fh = figure;
end

figure(fh)
set(fh,'color','k','position',[360 380 810 320])
ah = subplot(1,2,1); plotROI(ROI,imsize2,ah);
subplot(1,2,2)
set(gca,'color','k','xcolor','k','ycolor','k')
nROI = length(ROI);
for ii = 1:2:length(ROI)-1 
    text(.01,1-1/nROI*ii,sprintf('%01d. %s',ii,ROIlbl{ii}),'fontsize',14,'color','w'); 
    text(.51,1-1/nROI*ii,sprintf('%01d. %s',ii+1,ROIlbl{ii+1}),'fontsize',14,'color','w'); 
end
