function imstackReg = applyRegistrationToStack(imstack,tform,refim)

% imstackReg = applyRegistrationToStack(imstack,tform,refim)
% applies image transformation in tform to entire imstack

Rfixed     = imref2d(size(refim));
imstackReg = zeros(size(imstack));

for iZ = 1:size(imstack,3)
  % separate vasculature from brain image
  thisFrame                   = imstack(:,:,iZ) .* 100; % scale to decrease precision problems
  vasc                        = thisFrame;
  vasc(~isnan(vasc))          = 0;
  vasc(isnan(vasc))           = 10;
  thisFrame(isnan(thisFrame)) = nanmean(thisFrame(:));
	frameOut                    = imwarp(thisFrame,tform,'nearest','OutputView',Rfixed,'FillValues',nan); 
  vasc                        = imwarp(vasc,tform,'nearest','OutputView',Rfixed,'FillValues',nan); 
  vasc                        = vasc > 5; % account for some interpolation-related change in value
  frameOut(vasc)              = nan; 
  imstackReg(:,:,iZ)          = frameOut ./ 100;
end