function spatialCorr = weightSpatialCorr(weightMat,cfg)

% spatialCorr = weightSpatialCorr(weightMat,cfg)
% calculates pairwise correlation by distance for decoding weights
% weightMat is pxl x pxl x decoder matrix
% cfg is optional analysis parameters strcuture

if nargin < 2
  cfg.pxlSize   = 1000 / (widefieldParams.pxlPerMM / widefieldParams.dsFactor / 4);
  cfg.spaceBins = 0:cfg.pxlSize:cfg.pxlSize*15;
  cfg.nShuffles = 50;
end

try
rng(0)
tic; fprintf('calculating weight spatial correlation...')
%% here each pixel is compared against the others using all spatial positions to calculate
% a single correlation per pair
[nX,nY,~ ]     = size(weightMat);
nanpxl         = isnan(weightMat(:,:,1));
weightMat      = conditionDffMat(weightMat); 
colmat         = ones(nX,nY)*diag(1:nX);     
rowmat         = colmat';
colmat(nanpxl) = nan;
rowmat(nanpxl) = nan;
colvec         = conditionDffMat(repmat(colmat,[1 1 2])); 
colvec         = colvec(1,:);
rowvec         = conditionDffMat(repmat(rowmat,[1 1 2])); 
rowvec         = rowvec(1,:);


xcorrs   = cell(1,size(weightMat,2));
distBins = cfg.spaceBins;
pxlSize  = cfg.pxlSize;

for iPxl = 1:numel(xcorrs)
  thisrow  = rowvec - rowvec(iPxl);
  thiscol  = colvec - colvec(iPxl);
  thisdist = pxlSize .* (sqrt(thisrow.^2 + thiscol.^2));
  corrvec  = zeros(numel(distBins)-1,1);
  
  for iBin = 1:numel(distBins)-1
    idx           = thisdist >= distBins(iBin) & thisdist < distBins(iBin+1);
    corrvec(iBin) = nanmean(corr(weightMat(:,iPxl),weightMat(:,idx)));
  end
  
  xcorrs{iPxl}    = corrvec;
end

spatialCorr.distAxis        = distBins(1:end-1);
spatialCorr.corrByDist      = cell2mat(xcorrs);
spatialCorr.corrByDist_mean = nanmean(spatialCorr.corrByDist,2);

%% now shuffle pixels
xcorrs   = cell(cfg.nShuffles,size(weightMat,2));
for iSh = 1:cfg.nShuffles
  thisShuffle = weightMat(:,randperm(size(weightMat,2)));
  for iPxl = 1:size(xcorrs,2)
    thisrow  = rowvec - rowvec(iPxl);
    thiscol  = colvec - colvec(iPxl);
    thisdist = pxlSize .* (sqrt(thisrow.^2 + thiscol.^2));
    corrvec  = zeros(numel(distBins)-1,1);

    for iBin = 1:numel(distBins)-1
      idx            = thisdist >= distBins(iBin) & thisdist < distBins(iBin+1);
      corrvec(iBin)  = nanmean(corr(thisShuffle(:,iPxl),thisShuffle(:,idx)));
    end

    xcorrs{iSh,iPxl} = corrvec;
  end
end

shuffle                                = cell2mat(reshape(xcorrs,[1 size(weightMat,2) cfg.nShuffles]));    
spatialCorr.corrByDist_shuffle_mean    = nanmean(nanmean(shuffle,3),2);

%% significance vs shuffle
spatialCorr.corrByDist_isSig_VSshuffle = false(size(spatialCorr.corrByDist_mean));
for iBin = 1:numel(spatialCorr.corrByDist_isSig_VSshuffle)
  thiss = squeeze(shuffle(iBin,:,:));
  spatialCorr.corrByDist_isSig_VSshuffle(iBin) = spatialCorr.corrByDist_mean(iBin) > prctile(thiss(:),95);
end

fprintf(' done after %1.1f min\n',toc/60)

catch ME
  displayExcpetion(ME)
end