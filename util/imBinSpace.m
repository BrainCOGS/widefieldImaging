function imds = imBinSpace(imstack,vmask,dsFactor,doNanAvg)

% imds = imBinSpace(imstack,vmask,dsFactor,doNanAvg)
% spatially bins image by a factor of dsFactor (default in widefieldParams)
% vmask is optional input to take into account vasculature pixels (NaN)
% doNaNaveraging: if true  when binning excludes nans
%
% LP sep 2016

if nargin < 2
  vmask    = [];
end
if nargin < 3
  dsFactor = widefieldParams.dsFactor;
end
if nargin < 4
  doNanAvg = false;
end

% size
[nX,nY,nZ]  = size(imstack);
imds        = zeros(ceil(nX/dsFactor),ceil(nY/dsFactor),nZ);
xbins       = unique([1:dsFactor:nX nX+1]);
ybins       = unique([1:dsFactor:nY nY+1]);

if isempty(vmask)
  vmask   = false(nX,nY);
end

% downsample and ignore vasculature pixels
for zz = 1:nZ
  im        = imstack(:,:,zz);
  im(vmask) = nan;
  for xx = 1:length(xbins)-1
    for yy = 1:length(ybins)-1
      vals              = im(xbins(xx):xbins(xx+1)-1,ybins(yy):ybins(yy+1)-1);
      if doNanAvg
        imds(xx,yy,zz)  = nansum(vals(:));
      else
        imds(xx,yy,zz)  = sum(vals(:));
      end
    end
  end
end