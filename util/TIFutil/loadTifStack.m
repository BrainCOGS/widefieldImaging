function thisstack = loadTifStack(imFile)

imheader   = imfinfo(imFile);
    
readObj    = Tiff(imFile,'r');
thisstack  = zeros(imheader(1).Height,imheader(1).Width,numel(imheader),'uint16');
for iFrame = 1:numel(imheader)
  readObj.setDirectory(iFrame);
  thisstack(:,:,iFrame) = readObj.read();
end
readObj.close();