function [ROI,ROIlbl] = manualROIfromFieldSign(I,saveFlag)

% [ROI,ROIlbl] = manualROIfromFieldSign(I,saveFlag)
% utility to manually outline ROIs on field sign maps


if nargin < 2
  saveFlag = true;
end


if ~isempty(dir('manualROIfieldSign.mat'))
  load manualROIfieldSign ROI ROIlbl
  return
end

% manually select and label ROIs
ROIoutline   = zeros(size(I));
ROIcentroids = [];
figure; hold on; 
imagesc(I); colormap(red2blue)
axis image
axis off
stopFlag = 0; count = 1;
while ~stopFlag
    bw = roipoly;
    [r,c] = find(bw==1);
    ROIoutline            = ROIoutline+bw;
    ROI{count}            = [r c];
    ROIcentroids(count,:) = mean([r c]);
    plot(c,r,'k.')
    ROIlbl{count} = char(input('ROI label:  '));
    user_in = questdlg('select more ROIs?');
    if strcmpi(user_in,'yes')
        stopFlag = 0;
    else 
        stopFlag = 1;
    end
    count = count+1;
end
close

ROIperim = bwperim(ROIoutline);

if ~saveFlag; return; end

save manualROIfieldSign

figure;
hold on
imagesc(I); colormap red2blue
axis image
axis ij
axis off
[r,c]=find(ROIperim == 1);
for ii = 1:length(r)
    plot(c(ii),r(ii),'ks','markersize',1.2,'markerfacecolor','k');
end
for ii = 1:length(ROIlbl)
    text(ROIcentroids(ii,2),ROIcentroids(ii,1),ROIlbl{ii},'fontsize',13,'color','k');
end
saveas(gcf,'manualROIfieldSign')