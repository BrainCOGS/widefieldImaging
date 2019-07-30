function imstack = loadTiffStack(filelist,frameRange)

% imstack = loadTiffStack(filebasename,frameRange)
% loads tif stack from 1 or multiple files
% LP may 2013

warning off

if nargin == 1
    frameRange = [];
end

if ~iscell(filelist); temp{1} = filelist; filelist = temp; end

% figure out how many frames per file and image size
nZ = zeros(1,length(filelist));
tic;
fprintf('loading tiff stack...')
for i = 1:length(filelist)
    nZ(i) = tiff_Nframes(filelist{i});
end
temp = imread(filelist{i},1);
[nX nY] = size(temp); clear temp

if isempty(frameRange); frameRange = [1 sum(nZ)]; end
    
if length(filelist) == 1
    imstack = TIFFStack(filelist{i});
    if length(frameRange)>2
        imstack = imstack(:,:,frameRange);
    else
        imstack = imstack(:,:,frameRange(1):frameRange(end));
    end
else
    if length(frameRange)>2
        imstack = zeros(nX,nY,length(frameRange)); count = 1;
        step = mode(diff(frameRange));
        for i = 1:length(filelist)
            temp=TIFFStack(filelist{i});
            if count+nZ(i)-1 <= frameRange(end)
                imstack(:,:,count:count+nZ(i)/2-1)=temp(:,:,1:step:end);
            else imstack(:,:,count:frameRange(end)/2)=temp(:,:,1:step:frameRange(end)-count+1);
            end
            clear temp
            count = count + nZ(i);
        end
    else
        imstack = zeros(nX,nY,frameRange(end)); count = 1;
        for i = 1:length(filelist)
            temp=TIFFStack(filelist{i});
            if count+nZ(i)-1 <= frameRange(end)
                imstack(:,:,count:count+nZ(i)-1)=temp(:,:,:);
            else imstack(:,:,count:frameRange(end))=temp(:,:,1:frameRange(end)-count+1);
            end
            clear temp
            count = count + nZ(i);
        end
    end
end
t = toc; tic
fprintf(' done after %1.2f minutes\n',t/60)
end
