function allenROI = refROIfromAllen(areaLbl)

% [refFig,areaLbl,pxlSize,ROI,metaData] = refROIfromAllen(areaLbl)
% gets allen brain atlas structures specificied in areaLbl
% uses https://github.com/BaselLaserMouse/AllenBrainAPI.git repo, which I
% updated to reflect the CCF v3 
% also needs nrrdread.m, from file exchange
% names should be in allen conventions

if nargin < 1
    areaLbl = {'VISa','VISrl','VISpm','VISpl','VISal','VISam','VISl','VISp','RSP','SS','MOp','MOs'};
end

wf      = widefieldParams;
root    = formatFilePath(wf.getRootDir(isThisSpock));
fn      = [root '/allenROI.mat'];

%%
if ~isempty(dir(fn))
  load(fn,'allenROI')
  return
end

%% get reference atlas at 50 um resolution
[sections,metaData] = loadRefAtlas;
pxlSize             = 0.05; % mm
pxlPerMM            = 1/pxlSize;
[~,AP,ML]           = size(sections);
refFig              = zeros(AP,ML);
ROI                 = cell(1,numel(areaLbl));
bregma              = [106 114.5];

%% for each area, 
tic; fprintf('finding area coordinates')
for iArea = 1:numel(areaLbl)
  fprintf('.')
  %% get children from this parent structure, exclude L5/6
  parentID  = acronym2structureID(areaLbl{iArea});
  structLs  = getAllenStructureList('childrenOf',structureID2name(parentID));
  structLs  = structLs.acronym;
  if ~strcmpi(areaLbl{iArea},'PTLp')
    validLs = cellfun(@(x)(~isempty(strfind(x,areaLbl{iArea})) & ...
                           (isempty(strfind(x,'5')) | isempty(strfind(x,'6')))),structLs); 
  else
    validLs = cellfun(@(x)((~isempty(strfind(x,areaLbl{iArea}))  | ...
                            ~isempty(strfind(x,'VISa'))          | ...
                            ~isempty(strfind(x,'VISrl')))        & ...
                            (isempty(strfind(x,'5')) | isempty(strfind(x,'6')))),structLs); 
  end
  areaIDs   = acronym2structureID(structLs(validLs));
  
  for iCoronal = 1:AP
    thisc   = squeeze(sections(:,iCoronal,:));
    [~,mls] = find(ismember(thisc,areaIDs));
    mls     = unique(mls);
    prev    = refFig(iCoronal,:);
    refFig(iCoronal,mls)     = iArea;
    refFig(iCoronal,prev~=0) = prev(prev~=0);
  end
  
  %% get coordinates
  [r,c]      = find(refFig == iArea);
  ROI{iArea} = [r c];
  
end

%% save
fieldls = {'sections','metaData','pxlSize','pxlPerMM','AP','ML','refFig','ROI','bregma','areaLbl'};
for iF = 1:numel(fieldls)
  allenROI.(fieldls{iF}) = eval(fieldls{iF});
end
save(fn,'allenROI')

fprintf(' done after %1.1f min\n',toc/60)
end

%% load atlas data
function [sections,metaData] = loadRefAtlas(fpath)

if nargin < 1
  fpath = widefieldParams.allenRefPath;
end

%3-D matrix of grid-level annotation labels
[sections, metaData] = nrrdread(fpath);
sections             = double(sections);

end