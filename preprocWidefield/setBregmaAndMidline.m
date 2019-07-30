function [bregma,midline] = setBregmaAndMidline(meanproj,binFlag)

% [bregma,midline] = setBregmaAndMidline(meanproj,binFlag)
% user-set bregma coordinates + two points to specify the line going
% through the midline (will calculate the equation for the line)
% meanproj is refernce image
% binFlag to bin to 128 x 128 pixels (default true)
%
% LP sep 2016

if nargin < 1
    load dff meanproj 
end
if nargin < 2
    binFlag = 1;
end

if binFlag
    meanproj = fliplr(rot90(imBinSpace(meanproj),3));
end

% user select bregma
figure('name','select bregma');
imagesc(meanproj); colormap gray; axis image
bregma = round(ginput(1));
close

% user select two points to fit midline
figure('name','select two points for midline');
imagesc(meanproj); colormap gray; axis image
midline.p(1,:) = round(ginput(1));
midline.p(2,:) = round(ginput(1));
close

% as in y = ax + b
midline.a = (midline.p(1,2)-midline.p(2,2))/(midline.p(1,1)-midline.p(2,1));
midline.b = midline.p(1,2) - midline.a*midline.p(1,1);

save bregma
