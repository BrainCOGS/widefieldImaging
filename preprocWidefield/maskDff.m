function dff = maskDff(dff,doMask,rotateFlag,saveFlag)

% dff = maskDff(dff,doMask,rotateFlag,saveFlag)
% apply NaNs to dff image stack in pxls corresponding to true entries in
% doMask logical matrix. rotateFlag to rotate image by 270 degrees,
% saveFlag to save masked dff

if nargin < 1
    fprintf('loading dff...\n')
    load dff dff
end
if nargin < 2 || isempty(doMask)
    try
        load mask doMask
    catch
        err('found no mask, please provide as input');
    end
end
if nargin < 3
    rotateFlag = 1;
end
if nargin < 4
    saveFlag = 1;
end

if size(doMask,1) > size(dff,1)
    ds       = size(dff,1)/size(doMask,1);
    origMask = doMask;
    doMask   = imresize(doMask,ds);
else
    origMask = [];
end

% using repmat blows up memory
fprintf('applying mask...\n')
for ii = 1:size(dff,3); 
    frame         = dff(:,:,ii);
    frame(doMask) = nan;
    dff(:,:,ii)   = frame;
end

if rotateFlag
    fprintf('rotating...\n')
    dff    = fliplr(rot90(dff,3));
    doMask = fliplr(rot90(doMask,3));
end

if saveFlag
    fprintf('saving dff...\n')
    save('dff','dff','-append')
    save mask doMask origMask
end