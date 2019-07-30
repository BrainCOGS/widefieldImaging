function weightSpatialCorr_allDecoders(rec,cfg)

if nargin < 1 || isempty(rec)
  rec           = formatFilePath(pwd); 
end
if nargin < 2 || isempty(cfg)
  cfg.pxlSize   = 1000 / (widefieldParams.pxlPerMM / widefieldParams.dsFactor / 4);
  cfg.spaceBins = 0:cfg.pxlSize:cfg.pxlSize*15;
  cfg.nShuffles = 50;
end

%%
cd(rec)
filelist = dir('decode*binned4x*mat');
filelist = {filelist(:).name}';

for iFile = 1:numel(filelist)
  load(filelist{iFile},'decoder')
  decoder.weightsSpatialCorr = weightSpatialCorr(decoder.weights,cfg);
  save(filelist{iFile},'decoder','-append')
end
