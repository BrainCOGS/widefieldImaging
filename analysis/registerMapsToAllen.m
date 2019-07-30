function regMaps = registerMapsToAllen(mapPath)

% regMaps = registerMapsToAllen(mapPath)
% registers measured sensory maps to allen ROIs

if nargin < 1; mapPath   = [widefieldParams.savepath_mac 'maps_all.mat']; end

rowCutoffFieldSign = 70; % lot f junk anterior to this
%% load data
mapData   = load(mapPath);

%% retrieve registered field signs and airpuffs, gcamp
gcampRecs = mapData.visMap_gcamp.recls;
mice      = cellfun(@mouseAndDateFromFileName,gcampRecs,'UniformOutput',false);
wf        = widefieldParams;
rootdir   = wf.getRootDir(isThisSpock,true);
fieldSign = cell(numel(mice),1);
bregma    = cell(numel(mice),1);
midline   = cell(numel(mice),1);
meanproj  = cell(numel(mice),1);
allenROI  = cell(numel(mice),1);
overlay   = cell(numel(mice),1);
outline   = cell(numel(mice),1);
airpuff   = cell(numel(mice),1);

for iMouse = 1:numel(mice)
  load([rootdir mice{iMouse} '/refROI.mat'])
  
  % register field sign
  if iMouse == 1; regMaps.gcamp_avg.ROIlbl = refROI.areaLbl; end
  im                = refROI.sensoryMaps.visual.fieldSign;
  refim             = imresize(refROI.refIm,widefieldParams.dsFactor);
  im(isnan(im))     = 0;
  Rfixed            = imref2d(size(refim));
  fieldSign{iMouse} = imresize(imwarp(im,refROI.sensoryMaps.visual.tform,'OutputView',Rfixed),1/widefieldParams.dsFactor);
  
  % register airpuff
  im                = refROI.sensoryMaps.whiskers.respMap;
  im(isnan(im))     = 0;
  airpuff{iMouse}   = imresize(imwarp(im,refROI.sensoryMaps.whiskers.tform,'OutputView',Rfixed),1/widefieldParams.dsFactor);
  
  % compile variables
  meanproj{iMouse}  = repmat(255.*refROI.refIm./max(refROI.refIm(:)),[1 1 3]);
  bregma{iMouse}    = refROI.bregma;
  midline{iMouse}   = refROI.midline;
  allenROI{iMouse}  = refROI.areaCoord;
  outline{iMouse}   = getPerimCoords(allenROI{iMouse},size(meanproj{iMouse},1));
  
  % cleanup hemisphere artifacts
  fieldSign{iMouse} = flipLeftHemSign(fieldSign{iMouse},midline{iMouse},rowCutoffFieldSign);
  airpuff{iMouse}   = deleteRightHem(airpuff{iMouse},midline{iMouse});
  
  % generate combined figure
  overlay{iMouse}   = overlayFieldSign(meanproj{iMouse},fieldSign{iMouse},airpuff{iMouse});
end

%% compile data structure
regMaps.gcamp.mice = mice;
regMaps.gcamp.recs = gcampRecs;
fieldLs = {'fieldSign','meanproj','bregma','midline','allenROI','outline','overlay','airpuff'};
for iF = 1:numel(fieldLs)
  regMaps.gcamp.(fieldLs{iF}) = eval(fieldLs{iF});
end

%% average maps
refID         = 2;
avgfs         = fieldSign{refID};
avgpuff       = airpuff{refID};
avgbregma     = bregma{refID};
avgmid        = midline{refID}.p;
regROI        = cell(1,numel(mice));
regROI{refID} = allenROI{refID};
for iMouse = 1:numel(mice)
  if iMouse == refID; continue; end
  if strcmpi(mice{iMouse},'ai9')
    cf     = [2 3 -3 .99];
    tform  = geoTransfFromPxl(cf);
    Rfixed = imref2d(size(fieldSign{refID}));
    regfs  = imwarp(fieldSign{iMouse},tform,'OutputView',Rfixed); 
    [x,y]  = transformPointsInverse(tform,bregma{iMouse}(1),bregma{iMouse}(2));
    regbregma = round([x y]);
    [x,y] = transformPointsInverse(tform,midline{iMouse}.p(1,1),midline{iMouse}.p(1,2));
    regmidline.p(1,:) = [x y];
    [x,y] = transformPointsInverse(tform,midline{iMouse}.p(2,1),midline{iMouse}.p(2,2));
    regmidline.p(2,:) = [x y];
  
    % as in y = ax + b
    regmidline.a = (regmidline.p(1,2)- regmidline.p(2,2))/(regmidline.p(1,1)-regmidline.p(2,1));
    regmidline.b = regmidline.p(1,2) - regmidline.a*regmidline.p(1,1);

  else
    [tform,regfs,regbregma,regmidline] = registerRecs(fieldSign{iMouse},fieldSign{refID}, ...
                                                      bregma{iMouse},midline{iMouse},'rigid',false,'multimodal');
  end
  close;
  [nX,nY]   = size(airpuff{iMouse});
  Rfixed    = imref2d([nX nY]);
  regap     = imwarp(airpuff{iMouse},tform,'OutputView',Rfixed); 
  
  avgfs     = avgfs + regfs;
  avgpuff   = avgpuff + regap;
  avgbregma = avgbregma + regbregma;
  avgmid    = avgmid + regmidline.p;
  
  
  for iROI = 1:numel(allenROI{iMouse})
    [c,r]     = transformPointsForward(tform,allenROI{iMouse}{iROI}(:,2),allenROI{iMouse}{iROI}(:,1));
    r         = round(r);
    c         = round(c);
    delidx    = r < 1 | r > nX | c < 1 | c > nY;
    r(delidx) = [];
    c(delidx) = [];
    regROI{iMouse}{iROI}    = [r c];
  end
end

regMaps.gcamp_avg.fieldSign = avgfs./numel(mice);
regMaps.gcamp_avg.airpuff   = avgpuff./numel(mice);
regMaps.gcamp_avg.bregma    = avgbregma./numel(mice);
mid.p                       = avgmid./numel(mice);
mid.a                       = (mid.p(1,2)- mid.p(2,2)) / (mid.p(1,1) - mid.p(2,1));
mid.b                       = mid.p(1,2) - mid.a*mid.p(1,1);
regMaps.gcamp_avg.midline   = mid;

%% average allen
meanROI = cell(1,numel(regROI{iMouse}));
for iROI = 1:numel(regROI{iMouse})
  im = zeros(size(fieldSign{1},1),size(fieldSign{1},1),numel(mice));
  for iMouse = 1:numel(mice)
    for iPxl = 1:size(regROI{iMouse}{iROI},1) 
      im(regROI{iMouse}{iROI}(iPxl,1),regROI{iMouse}{iROI}(iPxl,2),iMouse) = 1; 
    end
  end
  im            = mean(im,3);
  [r,c]         = find(im > 1/numel(mice));
  meanROI{iROI} = [r c]; 
end

regMaps.gcamp_avg.ROI     = meanROI;
regMaps.gcamp_avg.outline = getPerimCoords(meanROI,nX);

end

%% left hemisphere pixels have inverted sign for analysis convenience, revert
function fieldSign = flipLeftHemSign(fieldSign,midline,rowCutoffFieldSign)

[r,c]      = find(~isnan(fieldSign));
[~,isLeft] = assignPxlToHem([r c],midline);
r          = r(isLeft);
c          = c(isLeft);
for iL = 1:numel(r)
  fieldSign(r(iL),c(iL)) = -fieldSign(r(iL),c(iL));
end
fieldSign(1:rowCutoffFieldSign,:) = 0;
end

%% keep just contralateral (L hem) maps
function im = deleteRightHem(im,midline)

[r,c]      = find(~isnan(im));
isRight    = assignPxlToHem([r c],midline);
r          = r(isRight);
c          = c(isRight);
for iL = 1:numel(r)
  im(r(iL),c(iL)) = 0;
end
end

%% smooth area borders to be used with plot function
function bound = getPerimCoords(ROI,imsize)

bound = cell(numel(ROI),1);
for iROI = 1:numel(ROI)
  im       = zeros(imsize);
  for iPxl = 1:size(ROI{iROI},1); im(ROI{iROI}(iPxl,1),ROI{iROI}(iPxl,2)) = 1; end
  if iROI == 3 || iROI == 4 % v2 is disjoint, requires special treatment
    se = strel('line',5,1);
    im = imdilate(im,se);
  end
  centroid = round(mean(ROI{iROI}));
  row      = centroid(1);
  col      = find(im(row,:) > 0, 1, 'first');
  bound{iROI} = bwtraceboundary(im,[row col],'N');
  if numel(bound{iROI}(:,1)) < 20
    col         = find(im(row,:) > 0, 1, 'last');
    bound{iROI} = bwtraceboundary(im,[row col],'S');
  end
  if numel(bound{iROI}(:,1)) < 20 
    col         = find(im(row,:) > 0, 2, 'first');
    bound{iROI} = bwtraceboundary(im,[row col(2)],'N');
  end
  if numel(bound{iROI}(:,1)) < 20
    keyboard; 
  end
end

end

%% create combined image
function meanproj = overlayFieldSign(meanproj,fieldSign,airpuff)

[r,c]      = find(fieldSign<0);
for iL = 1:numel(r)
  meanproj(r(iL),c(iL),1) = abs(fieldSign(r(iL),c(iL)))*255;
end

[r,c]      = find(fieldSign>0);
for iL = 1:numel(r)
  meanproj(r(iL),c(iL),3) = abs(fieldSign(r(iL),c(iL)))*255;
end

airpuff    = airpuff./max(airpuff(:));
[r,c]      = find(airpuff>0);
for iL = 1:numel(r)
  meanproj(r(iL),c(iL),2) = airpuff(r(iL),c(iL))*255;
end


meanproj   = uint8(meanproj);
end