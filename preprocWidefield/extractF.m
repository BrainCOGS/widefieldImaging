function rawf = extractF(imstack,vmask,ds)

% rawf = extractRawF(imstack,vmask)
% output is a nROI x frame matrix with rawf traces
% inputs:   imstack: 3d matrix or filename for loading
%           logical matrix where 1's indicate vasculature pixels
%           ds: downsampling factor
%
% LP may 2016

if nargin < 3
    ds = widefieldParams.dsFactor;
end

tic
fprintf('extracting raw fluorescence...')
% if imstack is a file name load mat file
if ischar(imstack)
    load imstackAlign imstack % assumes aligned stack is saved to mat file
end

% downsample
rawf = imBinSpace(imstack,vmask,ds);
rawf = rawf - widefieldParams.f_offset*ds; % because i'm summing instead of averaging

% save dff rawf dff dffCorrect correctFlag normMethod
t = toc;
fprintf(' done after %1.2f minutes\n',t/60)